!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rotate_int_den_to_local
  USE m_juDFT
  !**********************************************************************
  !     This subroutine calculates the spin-up and -down density, in the INT-region,
  !     e.i. it take the non-colinear density and rotates it locally into the
  !     spin-frame that make it spin-diagonal.
  !     The rotated density is needed to calculate the potential-energy integrals 
  !     in vgen_xcpot. For accuracy reasons, the magnetisation for the potential
  !     itself is regeneated from the unrotated densities.
  !     In addition this routine stores the angle used in the rotation.
  !     These angles are needed in vgen->vgen_finalize->vmatgen to rotate the up- and down-
  !     potentials back to the global frame.   DW 2018
  !     Based on rhodirgen by
  !     Philipp Kurz 99/11/01
  !**********************************************************************
CONTAINS
  SUBROUTINE rotate_int_den_to_local(DIMENSION,sym,stars,atoms,sphhar,vacuum,&
       cell,input,noco,oneD,den)
    !******** ABBREVIATIONS ***********************************************
    !     ifft3    : size of the 3d real space mesh
    !     ifft2    : size of the 2d real space mesh
    !     rpw      : diagonal components of the density matrix (rho_11 ,
    !                rho_22)
    !                later interstitial spin-up and -down density
    !                all stored in terms of 3d-stars
    !     ris      : first components of the density matrix
    !                later interstitial spin-up and -down density and
    !                direction of magnetic field (theta and phi)
    !                all stored on real space mesh
    !**********************************************************************

    USE m_constants
    USE m_fft2d
    USE m_fft3d
    USE m_types
    IMPLICIT NONE

    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_potden),INTENT(INOUT)   :: den

    !     .. Local Scalars ..
    INTEGER iden,jspin,ivac,ifft2,ifft3
    INTEGER imz,ityp,iri,ilh,imesh,iq2,iq3
    REAL   rho_11,rho_22,rho_21r,rho_21i,rhotot,magmom,phi
    REAL rho_up,rho_down,mx,my,mz,eps,vz_r,vz_i,rziw,theta
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: rz(:,:,:)
    REAL,    ALLOCATABLE :: rvacxy(:,:,:,:),ris(:,:),fftwork(:)
    !     ..
    eps = 1.0e-20



    !
    !---> initialize arrays for the density matrix
    !

    ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
    IF (input%film) THEN
       ifft2 = 9*stars%mx1*stars%mx2
       IF (oneD%odi%d1) ifft2 = 9*stars%mx3*oneD%odi%M
    ELSE
       ifft2=0
    END IF

    IF (ALLOCATED(den%phi_pw)) THEN
       DEALLOCATE(den%phi_pw,den%phi_vacz,den%phi_vacxy)
       DEALLOCATE(den%theta_pw,den%theta_vacz,den%theta_vacxy)
    ENDIF
    ALLOCATE(den%phi_pw(ifft3),den%theta_pw(ifft3))
    ALLOCATE(den%phi_vacz(vacuum%nmzd,2),den%theta_vacz(vacuum%nmzd,2))
    ALLOCATE(den%phi_vacxy(ifft2,vacuum%nmzxyd,2),den%theta_vacxy(ifft2,vacuum%nmzxyd,2))

     
    ALLOCATE (ris(ifft3,4),fftwork(ifft3))
    !---> fouriertransform the diagonal part of the density matrix
    !---> in the interstitial, den%pw, to real space (ris)
    DO iden = 1,2
       CALL fft3d(ris(:,iden),fftwork,den%pw(:,iden),stars,+1)
    ENDDO
    !---> fouriertransform the off-diagonal part of the density matrix
    CALL fft3d(ris(:,3),ris(:,4),den%pw(:,3),stars,+1)

    !test
    !      DO iden=1,4
    !         write(*,*)'iden=',iden
    !         write(*,8500)(ris(imesh,iden),imesh=0,ifft3-1)
    !      enddo
    !test
    !---> calculate the charge and magnetization density on the
    !---> real space mesh
    DO imesh = 1,ifft3
       rho_11  = ris(imesh,1)
       rho_22  = ris(imesh,2)
       rho_21r = ris(imesh,3)
       rho_21i = ris(imesh,4)
       mx      =  2*rho_21r
       my      = -2*rho_21i
       mz      = (rho_11-rho_22)
       magmom  = SQRT(mx**2 + my**2 + mz**2)
       rhotot  = rho_11 + rho_22
       rho_up  = (rhotot + magmom)/2
       rho_down= (rhotot - magmom)/2

       IF (ABS(mz) .LE. eps) THEN
          theta = pi_const/2
       ELSEIF (mz .GE. 0.0) THEN
          theta = ATAN(SQRT(mx**2 + my**2)/mz)
       ELSE
          theta = ATAN(SQRT(mx**2 + my**2)/mz) + pi_const
       ENDIF

       IF (ABS(mx) .LE. eps) THEN
          IF (ABS(my) .LE. eps) THEN
             phi = 0.0
          ELSEIF (my .GE. 0.0) THEN
             phi = pi_const/2
          ELSE
             phi = -pi_const/2
          ENDIF
       ELSEIF (mx .GE. 0.0) THEN
          phi = ATAN(my/mx)
       ELSE
          IF (my .GE. 0.0) THEN
             phi = ATAN(my/mx) + pi_const
          ELSE
             phi = ATAN(my/mx) - pi_const
          ENDIF
       ENDIF

       !         write(36,'(i4,2f12.6)') mod(imesh,33),rho_11,rho_22
       ris(imesh,1) = rho_up
       ris(imesh,2) = rho_down
       den%theta_pw(imesh) = theta
       den%phi_pw(imesh) = phi
    ENDDO

    DO jspin = 1,input%jspins
       fftwork=0.0
       CALL fft3d(ris(:,jspin),fftwork,den%pw(:,jspin),stars,-1)
    ENDDO

    IF (.NOT.input%film) RETURN


     !Now the vacuum part starts

   
    ALLOCATE(rvacxy(ifft2,vacuum%nmzxyd,2,4))
    ALLOCATE (rz(vacuum%nmzd,2,2))
    !---> fouriertransform the diagonal part of the density matrix
    !---> in the vacuum, rz & rxy, to real space (rvacxy)
    DO iden = 1,2
       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             rziw = 0.0
             IF (oneD%odi%d1) THEN
                CALL judft_error("oneD not implemented",calledby="rhodirgen")
                !CALL fft2d(oneD%k3,odi%M,odi%n2d,rvacxy(0,imz,ivac,iden),fftwork,&
                !           den%vacz(imz,ivac,iden),rziw,den%vacxy(imz,1,ivac,iden),&
                !           vacuum,odi%nq2,odi%kimax2,1,&
                !     &                  %igf,odl%pgf,odi%nst2)
             ELSE
                CALL fft2d(stars,rvacxy(:,imz,ivac,iden),fftwork,&
                     den%vacz(imz,ivac,iden),rziw,den%vacxy(imz,1,ivac,iden),&
                     vacuum%nmzxyd,1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !--->    fouriertransform the off-diagonal part of the density matrix
    DO ivac = 1,vacuum%nvac
       DO imz = 1,vacuum%nmzxyd
          rziw = 0.0
          vz_r = den%vacz(imz,ivac,3)
          vz_i = den%vacz(imz,ivac,4)
          IF (oneD%odi%d1) THEN
             CALL judft_error("oneD not implemented",calledby="rhodirgen")
             !CALL fft2d(oneD%k3,odi%M,odi%n2d,&
             !           rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),&
             !           vz_r,vz_i,den%vacxy(imz,1,ivac,3),&
             !           vacuum,odi%nq2,odi%kimax2,1,&
             !     &               %igf,odl%pgf,odi%nst2)
          ELSE
             CALL fft2d(stars,rvacxy(:,imz,ivac,3),rvacxy(:,imz,ivac,4),&
                  vz_r,vz_i,den%vacxy(imz,1,ivac,3),vacuum%nmzxyd,1)
          ENDIF
       ENDDO
    ENDDO

    !--->    calculate the four components of the matrix potential on
    !--->    real space mesh
    DO ivac = 1,vacuum%nvac
       DO imz = 1,vacuum%nmzxyd
          DO imesh = 1,ifft2
             rho_11  = rvacxy(imesh,imz,ivac,1)
             rho_22  = rvacxy(imesh,imz,ivac,2)
             rho_21r = rvacxy(imesh,imz,ivac,3)
             rho_21i = rvacxy(imesh,imz,ivac,4)
             mx      =  2*rho_21r
             my      = -2*rho_21i
             mz      = (rho_11-rho_22)
             magmom  = SQRT(mx**2 + my**2 + mz**2)
             rhotot  = rho_11 + rho_22
             rho_up  = (rhotot + magmom)/2
             rho_down= (rhotot - magmom)/2

             IF (ABS(mz) .LE. eps) THEN
                theta = pi_const/2
             ELSEIF (mz .GE. 0.0) THEN
                theta = ATAN(SQRT(mx**2 + my**2)/mz)
             ELSE
                theta = ATAN(SQRT(mx**2 + my**2)/mz) + pi_const
             ENDIF

             IF (ABS(mx) .LE. eps) THEN
                IF (ABS(my) .LE. eps) THEN
                   phi = 0.0
                ELSEIF (my .GE. 0.0) THEN
                   phi = pi_const/2
                ELSE
                   phi = -pi_const/2
                ENDIF
             ELSEIF (mx .GE. 0.0) THEN
                phi = ATAN(my/mx)
             ELSE
                IF (my .GE. 0.0) THEN
                   phi = ATAN(my/mx) + pi_const
                ELSE
                   phi = ATAN(my/mx) - pi_const
                ENDIF
             ENDIF

             rvacxy(imesh,imz,ivac,1) = rho_up
             rvacxy(imesh,imz,ivac,2) = rho_down
             den%theta_vacxy(imesh,imz,ivac) = theta
             den%phi_vacxy(imesh,imz,ivac) = phi
          ENDDO
       ENDDO
       DO imz = vacuum%nmzxyd+1,vacuum%nmzd
          rho_11  = den%vacz(imz,ivac,1)
          rho_22  = den%vacz(imz,ivac,2)
          rho_21r = den%vacz(imz,ivac,3)
          rho_21i = den%vacz(imz,ivac,4)
          mx      =  2*rho_21r
          my      = -2*rho_21i
          mz      = (rho_11-rho_22)
          magmom  = SQRT(mx**2 + my**2 + mz**2)
          rhotot  = rho_11 + rho_22
          rho_up  = (rhotot + magmom)/2
          rho_down= (rhotot - magmom)/2

          IF (ABS(mz) .LE. eps) THEN
             theta = pi_const/2
          ELSEIF (mz .GE. 0.0) THEN
             theta = ATAN(SQRT(mx**2 + my**2)/mz)
          ELSE
             theta = ATAN(SQRT(mx**2 + my**2)/mz) + pi_const
          ENDIF

          IF (ABS(mx) .LE. eps) THEN
             IF (ABS(my) .LE. eps) THEN
                phi = 0.0
             ELSEIF (my .GE. 0.0) THEN
                phi = pi_const/2
             ELSE
                phi = -pi_const/2
             ENDIF
          ELSEIF (mx .GE. 0.0) THEN
             phi = ATAN(my/mx)
          ELSE
             IF (my .GE. 0.0) THEN
                phi = ATAN(my/mx) + pi_const
             ELSE
                phi = ATAN(my/mx) - pi_const
             ENDIF
          ENDIF

          den%vacz(imz,ivac,1) = rho_up
          den%vacz(imz,ivac,2) = rho_down
          den%theta_vacz(imz,ivac) = theta
          den%phi_vacz(imz,ivac) = phi
       ENDDO
    ENDDO
    !--->    Fouriertransform the matrix potential back to reciprocal space
    DO jspin = 1,input%jspins
       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             fftwork=0.0
             IF (oneD%odi%d1) THEN
                call judft_error("oneD not implemented",calledby="rhodirgen")
                !CALL fft2d(oneD%k3,odi%M,odi%n2d,&
                !           rvacxy(0,imz,ivac,jspin),fftwork,&
                !           den%vacz(imz,ivac,jspin),rziw,den%vacxy(imz,1,ivac,jspin),&
                !           vacuum,odi%nq2,odi%kimax2,-1,&
                !     &                  %igf,odl%pgf,odi%nst2)
             ELSE
                CALL fft2d(stars,rvacxy(:,imz,ivac,jspin),fftwork,&
                     den%vacz(imz,ivac,jspin),rziw,den%vacxy(imz,1,ivac,jspin),&
                     vacuum%nmzxyd,-1)
             END IF
          ENDDO
       ENDDO
    ENDDO
    
    RETURN
  END SUBROUTINE rotate_int_den_to_local
END MODULE m_rotate_int_den_to_local
