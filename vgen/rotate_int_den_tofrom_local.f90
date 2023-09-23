!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rotate_int_den_tofrom_local
   USE m_juDFT
   USE m_fft2d
   USE m_fft3d
   USE m_types
   
   IMPLICIT NONE

CONTAINS
   
   SUBROUTINE rotate_int_den_to_local(sym, stars, atoms, sphhar, vacuum, cell, &
                                      input, noco,   den)

      !--------------------------------------------------------------------------
      ! This subroutine calculates the spin-up and -down density in the intersti-
      ! tial region i.e. it takes the non-collinear density and rotates it into
      ! a local spin-frame, making it spin-diagonal.
      ! 
      ! The rotated density is needed to calculate the potential-energy integrals
      ! in vgen_xcpot. For accuracy reasons, the magnetisation for the potential
      ! itself is regenerated from the unrotated densities.
      ! 
      ! In addition this routine stores the angles used in the rotation. These
      ! angles are later needed to rotate the up- and down-potentials back to the
      ! global frame. DW 2018
      ! 
      ! Based on rhodirgen by
      ! Philipp Kurz 99/11/01
      !--------------------------------------------------------------------------

      !-------Important variables:----------------------------------------------- 
      ! ifft3: size of the 3d real space mesh
      ! ifft2: size of the 2d real space mesh
      ! ris:   first components of the density matrix
      !        later interstitial spin-up and -down density and direction of mag-
      !        netic field (theta and phi) all stored on real space mesh
      !--------------------------------------------------------------------------

      USE m_constants
      USE m_polangle
    
      TYPE(t_noco),   INTENT(IN)    :: noco
       
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_stars),  INTENT(IN)    :: stars
      TYPE(t_cell),   INTENT(IN)    :: cell
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_potden), INTENT(INOUT) :: den

      INTEGER                       :: iden, jspin, ivac, ifft2, ifft3
      INTEGER                       :: imz, ityp, iri, ilh, imesh, iq2, iq3
      REAL                          :: rho_11, rho_22, rho_21r, rho_21i
      REAL                          :: mx, my, mz, magmom, vz_r, vz_i, rziw
      REAL                          :: rhotot, rho_up, rho_down, theta, phi
      REAL                          :: eps=1E-20

      REAL, ALLOCATABLE             :: rz(:,:,:), rvacxy(:,:,:,:), ris(:,:)
      REAL, ALLOCATABLE             :: fftwork(:)

      ! Initialize arrays for the density matrix:

      ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
      IF (input%film) THEN
         ifft2 = 9*stars%mx1*stars%mx2
      ELSE
         ifft2=0
      END IF

      IF (ALLOCATED(den%phi_pw)) THEN
         DEALLOCATE(den%phi_pw)!,den%phi_vacz,den%phi_vacxy)
         DEALLOCATE(den%theta_pw)!,den%theta_vacz,den%theta_vacxy)
         DEALLOCATE(den%theta_vac,den%phi_vac)
      END IF

      ALLOCATE(den%phi_pw(ifft3),den%theta_pw(ifft3))
      !ALLOCATE(den%phi_vacz(vacuum%nmzd,2),den%theta_vacz(vacuum%nmzd,2))
      ALLOCATE(den%phi_vac(ifft2,vacuum%nmzd,2),den%theta_vac(ifft2,vacuum%nmzd,2))
      !ALLOCATE(den%phi_vacxy(ifft2,vacuum%nmzxyd,2),den%theta_vacxy(ifft2,vacuum%nmzxyd,2))

      ALLOCATE (ris(ifft3,4),fftwork(ifft3))
 
      ! Interstitial part:

      ! Fourier transform the diagonal part of the density matrix (den%pw)
      ! to real space (ris):
      DO iden = 1,2
         CALL fft3d(ris(:,iden),fftwork,den%pw(:,iden),stars,+1)
      END DO

      ! Fourier transform the off-diagonal part of the density matrix:
      CALL fft3d(ris(:,3),ris(:,4),den%pw(:,3),stars,+1)

      ! Calculate the charge and magnetization densities on the real space mesh:
      DO imesh = 1,ifft3
         rho_11   = ris(imesh,1)
         rho_22   = ris(imesh,2)
         rho_21r  = ris(imesh,3)
         rho_21i  = ris(imesh,4)
         mx       =  2*rho_21r
         my       = -2*rho_21i ! TODO: This is a magic minus.
         mz       = rho_11 - rho_22
         magmom   = SQRT(mx**2 + my**2 + mz**2)
         rhotot   = rho_11 + rho_22
         rho_up   = (rhotot + magmom)/2
         rho_down = (rhotot - magmom)/2

         CALL pol_angle(mx,my,mz,theta,phi)

         ris(imesh,1) = rho_up
         ris(imesh,2) = rho_down
         den%theta_pw(imesh) = theta
         den%phi_pw(imesh) = phi
      END DO

      ! Fourier transform the matrix potential back to reciprocal space:
      DO jspin = 1, input%jspins
         fftwork=0.0
         CALL fft3d(ris(:,jspin),fftwork,den%pw(:,jspin),stars,-1)
      END DO

      IF (.NOT.input%film) RETURN

      ! Now the vacuum part starts:
   
      ALLOCATE(rvacxy(ifft2,vacuum%nmzxyd,2,4))
      ALLOCATE (rz(vacuum%nmzd,2,2))

      ! Fourier transform the diagonal part of the density matrix (den%vacz and
      ! den%vacxy) to real space (rvacxy):
      DO iden = 1,2
         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               rziw = 0.0
                
                  !CALL fft2d(stars,rvacxy(:,imz,ivac,iden),fftwork,&
                  !     den%vacz(imz,ivac,iden),rziw,den%vacxy(imz,:,ivac,iden),&
                  !     1)
                  CALL fft2d(stars,rvacxy(:,imz,ivac,iden),fftwork,&
                       REAL(den%vac(imz,1,ivac,iden)),rziw,den%vac(imz,2:,ivac,iden),&
                       1)
            END DO
         END DO
      END DO

      ! Fourier transform the off-diagonal part of the density matrix:
      DO ivac = 1,vacuum%nvac
         DO imz = 1,vacuum%nmzxyd
            rziw = 0.0
            !vz_r = den%vacz(imz,ivac,3)
            !vz_i = den%vacz(imz,ivac,4)
            vz_r = REAL(den%vac(imz,1,ivac,3))
            vz_i = AIMAG(den%vac(imz,1,ivac,3))
             
               !CALL fft2d(stars,rvacxy(:,imz,ivac,3),rvacxy(:,imz,ivac,4),&
               !     vz_r,vz_i,den%vacxy(imz,:,ivac,3),1)
               CALL fft2d(stars,rvacxy(:,imz,ivac,3),rvacxy(:,imz,ivac,4),&
                    vz_r,vz_i,den%vac(imz,2:,ivac,3),1)
            
         END DO
      END DO

      ! Calculate the charge and magnetization densities on the real space mesh:
      DO ivac = 1,vacuum%nvac
         DO imz = 1,vacuum%nmzxyd
            DO imesh = 1,ifft2
               rho_11   = rvacxy(imesh,imz,ivac,1)
               rho_22   = rvacxy(imesh,imz,ivac,2)
               rho_21r  = rvacxy(imesh,imz,ivac,3)
               rho_21i  = rvacxy(imesh,imz,ivac,4)
               mx       =  2*rho_21r
               my       = -2*rho_21i
               mz       = rho_11 - rho_22
               magmom   = SQRT(mx**2 + my**2 + mz**2)
               rhotot   = rho_11 + rho_22
               rho_up   = (rhotot + magmom)/2
               rho_down = (rhotot - magmom)/2

               CALL pol_angle(mx,my,mz,theta,phi)

               rvacxy(imesh,imz,ivac,1) = rho_up
               rvacxy(imesh,imz,ivac,2) = rho_down
               !den%theta_vacxy(imesh,imz,ivac) = theta
               !den%phi_vacxy(imesh,imz,ivac) = phi
               den%theta_vac(imesh,imz,ivac) = theta
               den%phi_vac(imesh,imz,ivac) = phi
            END DO
         END DO
       
         DO imz = vacuum%nmzxyd+1,vacuum%nmzd
            !rho_11   = den%vacz(imz,ivac,1)
            !rho_22   = den%vacz(imz,ivac,2)
            !rho_21r  = den%vacz(imz,ivac,3)
            !rho_21i  = den%vacz(imz,ivac,4)
            rho_11   = REAL(den%vac(imz,1,ivac,1))
            rho_22   = REAL(den%vac(imz,1,ivac,2))
            rho_21r  = REAL(den%vac(imz,1,ivac,3))
            rho_21i  = AIMAG(den%vac(imz,1,ivac,3))
            mx       =  2*rho_21r
            my       = -2*rho_21i
            mz       = rho_11 - rho_22
            magmom   = SQRT(mx**2 + my**2 + mz**2)
            rhotot   = rho_11 + rho_22
            rho_up   = (rhotot + magmom)/2
            rho_down = (rhotot - magmom)/2

            CALL pol_angle(mx,my,mz,theta,phi)

            !den%vacz(imz,ivac,1) = rho_up
            !den%vacz(imz,ivac,2) = rho_down
            den%vac(imz,1,ivac,1) = rho_up
            den%vac(imz,1,ivac,2) = rho_down
            !den%theta_vacz(imz,ivac) = theta
            !den%phi_vacz(imz,ivac) = phi
            den%theta_vac(1,imz,ivac) = theta
            den%phi_vac(1,imz,ivac) = phi
         END DO
      END DO
    
      ! Fourier transform the matrix potential back to reciprocal space:
      DO jspin = 1,input%jspins
         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               fftwork=0.0
                
                  !CALL fft2d(stars,rvacxy(:,imz,ivac,jspin),fftwork,&
                  !     den%vacz(imz,ivac,jspin),rziw,den%vacxy(imz,:,ivac,jspin),&
                  !     -1)
                  CALL fft2d(stars,rvacxy(:,imz,ivac,jspin),fftwork,&
                       REAL(den%vac(imz,1,ivac,jspin)),rziw,den%vac(imz,2:,ivac,jspin),&
                       -1) ! TODO: AN, TB: This likely won't work. Dummygröße schreiben/fft2d fixen!
               
            END DO
         END DO
      END DO
    
      RETURN

   END SUBROUTINE rotate_int_den_to_local

   SUBROUTINE rotate_int_den_from_local(stars,atoms,vacuum,sym,input,den,vTot)

      !--------------------------------------------------------------------------
      ! This subroutine prepares the spin dependent 2x2 matrix potential for the 
      ! Hamiltonian setup. This is done in 4 steps.
      ! 
      ! i)   The spin up and down potential and the angles of the magentic field, 
      !      theta and phi, are reloaded from den.
      ! ii)  The spin up and down potential is Fourier transformed to real space 
      !      (theta and phi are also stored on the real space grid).
      ! iii) The four components of the matrix potential are calculated on the
      !      real space mesh.
      ! iv)  The matrix potential is Fourier transformed, stored in terms of
      !      stars and written to vTot%pw(_w).
      ! 
      ! Philipp Kurz 99/11/01
      !--------------------------------------------------------------------------

      !-------Important variables:----------------------------------------------- 
      ! ifft3: size of the 3d real space mesh
      ! ifft2: size of the 2d real space mesh
      ! vis: first interstitial spin up and down potential and angles of magnetic
      !      field (theta and phi)
      !      later four components of matrix potential all stored in real space
      !--------------------------------------------------------------------------

      ! 
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_stars),  INTENT(IN)    :: stars
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_potden), INTENT(IN)    :: den
      TYPE(t_potden), INTENT(INOUT) :: vTot
 
      INTEGER                       :: imeshpt, ipot, jspin, ig2, ig3, ivac
      INTEGER                       :: ifft2, ifft3, imz, iter, i
      REAL                          :: vup, vdown, veff, beff, vziw, theta, phi

      REAL, ALLOCATABLE             :: vvacxy(:,:,:,:), vis(:,:), vis2(:,:)
      REAL, ALLOCATABLE             :: fftwork(:)

      ! Initialize arrays for the potential matrix:

      ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
      IF (ifft3.NE.SIZE(den%theta_pw)) CALL judft_error("Wrong size of angles")
      ifft2 = SIZE(den%phi_vac,1) 
    
      ALLOCATE ( vis(ifft3,4),fftwork(ifft3),vis2(ifft3,4))
    
      ! Interstitial part:

      ! Fourier transform the diagonal part of the potential matrix (vTot%pw)
      ! to real space (vis):
      DO jspin = 1,input%jspins
         CALL fft3d(vis(:,jspin),fftwork, vTot%pw(:,jspin), stars,+1)
      END DO

      ! Calculate the four components of the matrix potential on the real space
      ! mesh:
      DO imeshpt = 1, ifft3
         vup   = vis(imeshpt,1)
         vdown = vis(imeshpt,2)
         theta = den%theta_pw(imeshpt)
         phi   = den%phi_pw(imeshpt)

         veff  = (vup + vdown)/2.0
         beff  = (vup - vdown)/2.0

         vis(imeshpt,1) = veff + beff*COS(theta) ! V_(1,1) [V+B_z]
         vis(imeshpt,2) = veff - beff*COS(theta) ! V_(2,2) [V-B_z]
         vis(imeshpt,3) = beff*SIN(theta)*COS(phi) ! Re(V_(2,1)) [B_x]
         vis(imeshpt,4) = beff*SIN(theta)*SIN(phi) ! Im(V_(2,1)) [B_y]

         DO ipot = 1,4
            vis2(imeshpt,ipot) =  vis(imeshpt,ipot) * stars%ufft(imeshpt-1)
         END DO
      END DO

      ! Fourier transform the matrix potential back to reciprocal space:
      DO ipot = 1,2
         fftwork=0.0
         CALL fft3d(vis(:,ipot),fftwork, vTot%pw(1,ipot), stars,-1)
         fftwork=0.0
         CALL fft3d(vis2(:,ipot),fftwork, vTot%pw_w(1,ipot), stars,-1)
      END DO
    
      CALL fft3d(vis(:,3),vis(:,4), vTot%pw(1,3), stars,-1)
      CALL fft3d(vis2(:,3),vis2(:,4), vTot%pw_w(1,3), stars,-1)

      IF (.NOT. input%film) RETURN

      ! Now the vacuum part starts:

      ALLOCATE(vvacxy(ifft2,vacuum%nmzxyd,2,4))
    
      ! Fourier transform the up and down potentials (vTot%vacz and vTot%vacxy)
      ! to real space (vvacxy):
      DO jspin = 1,input%jspins
         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               vziw = 0.0
               ! 
                  CALL fft2d(stars, vvacxy(:,imz,ivac,jspin),fftwork,&
                       REAL(vTot%vac(imz,1,ivac,jspin)),vziw,vTot%vac(imz,2:,ivac,jspin), 1)
               
            END DO
         END DO
      END DO

      ! Calculate the four components of the matrix potential in real space:
      DO ivac = 1,vacuum%nvac
         DO imz = 1,vacuum%nmzxyd
            DO imeshpt = 1,ifft2
               vup   = vvacxy(imeshpt,imz,ivac,1)
               vdown = vvacxy(imeshpt,imz,ivac,2)
               !theta = den%theta_vacxy(imeshpt,imz,ivac)
               !phi   = den%phi_vacxy(imeshpt,imz,ivac)
               theta = den%theta_vac(imeshpt,imz,ivac)
               phi   = den%phi_vac(imeshpt,imz,ivac)

               veff  = (vup + vdown)/2.0
               beff  = (vup - vdown)/2.0
               vvacxy(imeshpt,imz,ivac,1) = veff + beff*COS(theta)
               vvacxy(imeshpt,imz,ivac,2) = veff - beff*COS(theta)
               vvacxy(imeshpt,imz,ivac,3) = beff*SIN(theta)*COS(phi)
               vvacxy(imeshpt,imz,ivac,4) = beff*SIN(theta)*SIN(phi)
            END DO
         END DO
          
         DO imz = vacuum%nmzxyd+1,vacuum%nmzd
            vup   = REAL(vTot%vac(imz,1,ivac,1))
            vdown = REAL(vTot%vac(imz,1,ivac,2))
            !theta = den%theta_vacz(imz,ivac)
            !phi   = den%phi_vacz(imz,ivac)
            theta = den%theta_vac(1,imz,ivac)
            phi   = den%phi_vac(1,imz,ivac)
            veff  = (vup + vdown)/2.0
            beff  = (vup - vdown)/2.0
            vTot%vac(imz,1,ivac,1) = veff + beff*COS(theta)
            vTot%vac(imz,1,ivac,2) = veff - beff*COS(theta)
            vTot%vac(imz,1,ivac,3) = beff*SIN(theta)*COS(phi)+ImagUnit*beff*SIN(theta)*SIN(phi)
         END DO
      END DO

      ! Fourier transform the matrix potential back to reciprocal space:
      DO ipot = 1,2
         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               fftwork=0.0
               ! 
                  CALL fft2d(stars, vvacxy(:,imz,ivac,ipot),fftwork,&
                       REAL(vTot%vac(imz,1,ivac,ipot)),vziw,vTot%vac(imz,2:,ivac,ipot),-1)
               
            END DO
         END DO
      END DO

      DO ivac = 1,vacuum%nvac
         DO imz = 1,vacuum%nmzxyd
            fftwork=0.0
            ! 
               CALL fft2d(stars, vvacxy(:,imz,ivac,3),vvacxy(:,imz,ivac,4),&
                    REAL(vTot%vac(imz,1,ivac,3)),AIMAG(vTot%vac(imz,1,ivac,3)),vTot%vac(imz,2:,ivac,3),-1)
            
         END DO
      END DO

      RETURN

   END SUBROUTINE rotate_int_den_from_local

END MODULE m_rotate_int_den_tofrom_local
