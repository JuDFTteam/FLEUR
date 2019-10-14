!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vmatgen
  USE m_juDFT
  !**********************************************************************
  !     This subroutine prepares the spin dependent 2x2 matrix potential
  !     for the Hamiltonian setup. This is done in 4 steps.
  !
  !    i) The spin up and down potential and the direction of the
  !     magentic field, theta and phi, are reloaded from files nrp,
  !     dirofmag.
  !    ii) The spin up and down potential is Fouriertransformed to real
  !     space (theta and phi are stored in real space).
  !    iii) The four components of the matrix potential are calculated on
  !     the real space mesh.
  !    iv) The matrix potential is Fouriertransformed, stored in terms of
  !     stars and written to file potmat.
  !
  !     Philipp Kurz 99/11/01
  !   
  !**********************************************************************
CONTAINS
  SUBROUTINE vmatgen(stars,atoms,vacuum,sym,input,den,vTot)

    !******** ABBREVIATIONS ***********************************************
    !     ifft3    : size of the 3d real space mesh
    !     ifft2    : size of the 2d real space mesh
    !     vpw      : first interstitial spin up and down potential
    !                later four components of matrix potential
    !                all stored in terms of 3d-stars
    !     vis      : first interstitial spin up and down potential and
    !                direction of magnetic field (theta and phi)
    !                later four components of matrix potential
    !                all stored on real space mesh
    !**********************************************************************

    USE m_fft2d
    USE m_fft3d
    USE m_types
    IMPLICIT NONE
!    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_vacuum),INTENT(IN) :: vacuum
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_stars),INTENT(IN)  :: stars
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_potden),INTENT(IN) :: den
    TYPE(t_potden),INTENT(INOUT):: vTot
 
    !     ..
    !     .. Local Scalars ..
    INTEGER imeshpt,ipot,jspin,ig2 ,ig3,ivac,ifft2,ifft3,imz,iter,i
    REAL    vup,vdown,veff,beff,vziw,theta,phi
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: vvacxy(:,:,:,:),vis(:,:),fftwork(:)

    ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
    IF (ifft3.NE.SIZE(den%theta_pw)) CALL judft_error("Wrong size of angles")
    ifft2 = SIZE(den%phi_vacxy,1) 
    
    ALLOCATE ( vis(ifft3,4),fftwork(ifft3))
      
    !---> fouriertransform the spin up and down potential
    !---> in the interstitial, vpw, to real space (vis)
    DO jspin = 1,input%jspins
       CALL fft3d(vis(:,jspin),fftwork, vTot%pw(:,jspin), stars,+1)
    ENDDO

    !---> calculate the four components of the matrix potential on
    !---> real space mesh
    DO imeshpt = 1,ifft3
       vup   = vis(imeshpt,1)
       vdown = vis(imeshpt,2)
       theta = den%theta_pw(imeshpt)
       phi   = den%phi_pw(imeshpt)
       !--->    at first determine the effective potential and magnetic fields
       veff  = (vup + vdown)/2.0
       beff  = (vup - vdown)/2.0
       !--->    now calculate the matrix potential, which is hermitian.
       !--->    thus calculate the diagonal elements:
       !--->    V_11
       vis(imeshpt,1) = veff + beff*COS(theta)
       !--->    V_22
       vis(imeshpt,2) = veff - beff*COS(theta)
       !--->    the real part of V_21
       vis(imeshpt,3) = beff*SIN(theta)*COS(phi)
       !--->    and the imaginary part of V_21
       vis(imeshpt,4) = beff*SIN(theta)*SIN(phi)

       DO ipot = 1,4
          vis(imeshpt,ipot) =  vis(imeshpt,ipot) * stars%ufft(imeshpt-1)
       ENDDO

    ENDDO

    !---> Fouriertransform the matrix potential back to reciprocal space
    DO ipot = 1,2
       fftwork=0.0
       CALL fft3d(vis(:,ipot),fftwork, vTot%pw_w(1,ipot), stars,-1)
    ENDDO
    
    CALL fft3d(vis(:,3),vis(:,4), vTot%pw_w(1,3), stars,-1)

    IF (.NOT. input%film) RETURN

    !Now the vacuum part starts

 
    ALLOCATE(vvacxy(ifft2,vacuum%nmzxyd,2,4))
    
       !--->    fouriertransform the spin up and down potential
       !--->    in the vacuum, vz & vxy, to real space (vvacxy)
       DO jspin = 1,input%jspins
          DO ivac = 1,vacuum%nvac
             DO imz = 1,vacuum%nmzxyd
                vziw = 0.0
                !IF (oneD%odi%d1) THEN
                IF (.FALSE.) THEN
                   CALL judft_error("oneD not implemented",calledby="vmatgen")
                   !                  CALL fft2d(&
                   !     &                 oneD%k3,odi%M,odi%n2d,&
                   !     &                 vvacxy(0,imz,ivac,jspin),fftwork,&
                   !     &                 vz(imz,ivac,jspin),vziw,vxy(imz,1,ivac,jspin),&
                   !     &                 vacuum,odi%nq2,odi%kimax2,1,&
                   !     &                  %igf,odl%pgf,odi%nst2)
                ELSE
                   CALL fft2d(stars, vvacxy(:,imz,ivac,jspin),fftwork,&
                        vTot%vacz(imz,ivac,jspin),vziw,vTot%vacxy(imz,1,ivac,jspin), vacuum%nmzxyd,1)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       !--->    calculate the four components of the matrix potential on
       !--->    real space mesh
       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             DO imeshpt = 1,ifft2
                vup   = vvacxy(imeshpt,imz,ivac,1)
                vdown = vvacxy(imeshpt,imz,ivac,2)
                theta = den%theta_vacxy(imeshpt,imz,ivac)
                phi   = den%phi_vacxy(imeshpt,imz,ivac)
                veff  = (vup + vdown)/2.0
                beff  = (vup - vdown)/2.0
                vvacxy(imeshpt,imz,ivac,1) = veff + beff*COS(theta)
                vvacxy(imeshpt,imz,ivac,2) = veff - beff*COS(theta)
                vvacxy(imeshpt,imz,ivac,3) = beff*SIN(theta)*COS(phi)
                vvacxy(imeshpt,imz,ivac,4) = beff*SIN(theta)*SIN(phi)
             ENDDO
          ENDDO
          DO imz = vacuum%nmzxyd+1,vacuum%nmzd
             vup   = vTot%vacz(imz,ivac,1)
             vdown = vTot%vacz(imz,ivac,2)
             theta = den%theta_vacz(imz,ivac)
             phi   = den%phi_vacz(imz,ivac)
             veff  = (vup + vdown)/2.0
             beff  = (vup - vdown)/2.0
             vTot%vacz(imz,ivac,1) = veff + beff*COS(theta)
             vTot%vacz(imz,ivac,2) = veff - beff*COS(theta)
             vTot%vacz(imz,ivac,3) = beff*SIN(theta)*COS(phi)
             vTot%vacz(imz,ivac,4) = beff*SIN(theta)*SIN(phi)
          ENDDO
       ENDDO

       !--->    Fouriertransform the matrix potential back to reciprocal space
       DO ipot = 1,2
          DO ivac = 1,vacuum%nvac
             DO imz = 1,vacuum%nmzxyd
                fftwork=0.0
                !IF (oneD%odi%d1) THEN
                IF (.FALSE.) THEN
                   CALL judft_error("oneD not implemented",calledby="vmatgen")
                   !                CALL fft2d(&
                   !     &                 oneD%k3,odi%M,odi%n2d,&
                   !     &                 vvacxy(0,imz,ivac,ipot),fftwork,&
                   !     &                 vz(imz,ivac,ipot),vziw,vxy(imz,1,ivac,ipot),&
                   !     &                 vacuum,odi%nq2,odi%kimax2,-1,&
                   !     &                  %igf,odl%pgf,odi%nst2)
                ELSE
                   CALL fft2d(stars, vvacxy(:,imz,ivac,ipot),fftwork,&
                        vTot%vacz(imz,ivac,ipot),vziw,vTot%vacxy(imz,1,ivac,ipot), vacuum%nmzxyd,-1)
                END IF
             ENDDO
          ENDDO
       ENDDO

       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             fftwork=0.0
             !IF (oneD%odi%d1) THEN
             IF (.FALSE.) THEN
             CALL judft_error("oneD not implemented",calledby="vmatgen")
                !              CALL fft2d(&
                !   &              oneD%k3,odi%M,odi%n2d,&
                !   &              vvacxy(0,imz,ivac,3),vvacxy(0,imz,ivac,4),&
                !   &              vz(imz,ivac,3),vz(imz,ivac,4),vxy(imz,1,ivac,3),&
                !   &              vacuum,odi%nq2,odi%kimax2,-1,&
                !   &               %igf,odl%pgf,odi%nst2)
             ELSE
                CALL fft2d(stars, vvacxy(:,imz,ivac,3),vvacxy(:,imz,ivac,4),&
                     vTot%vacz(imz,ivac,3),vTot%vacz(imz,ivac,4),vTot%vacxy(imz,1,ivac,3), vacuum%nmzxyd,-1)
             END IF
          ENDDO
       ENDDO

    RETURN
  END SUBROUTINE vmatgen
END MODULE m_vmatgen
