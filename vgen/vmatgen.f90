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
  !**********************************************************************
CONTAINS
  SUBROUTINE vmatgen(&
       &                   stars,&
       &                   atoms,sphhar,vacuum,&
       &                   sym,input,oneD,&
       &                   nu,ndomfile,npotmatfile)

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

    USE m_loddop
    USE m_fft2d
    USE m_fft3d
    USE m_types
    IMPLICIT NONE
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_vacuum),INTENT(IN) :: vacuum
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_stars),INTENT(IN)  :: stars
    TYPE(t_sphhar),INTENT(IN) :: sphhar
    TYPE(t_atoms),INTENT(IN)  :: atoms

    !     .. Scalar Arguments ..    
    INTEGER, INTENT (IN) :: nu,ndomfile,npotmatfile  

    !     ..
    !     .. Local Scalars ..
    INTEGER imeshpt,ipot,jspin,ig2 ,ig3,ivac,ifft2,ifft3,imz,iter
    REAL    vup,vdown,veff,beff  ,zero,vziw,theta,phi
    LOGICAL l_domfexst
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: vpw(:,:),vxy(:,:,:,:)
    REAL,    ALLOCATABLE :: vr(:,:,:,:),vz(:,:,:)
    REAL,    ALLOCATABLE :: vvacxy(:,:,:,:),vis(:,:),fftwork(:)

    zero = 0.0
    ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
    ifft2 = 9*stars%mx1*stars%mx2
    IF (oneD%odi%d1) ifft2 = 9*stars%mx3*oneD%odi%M
    IF (input%film) ALLOCATE(vvacxy(0:ifft2-1,vacuum%nmzxyd,2,4))

    IF (input%jspins .NE. 2) THEN
       WRITE (6,*) 'This is the non-collinear version of the flapw-'
       WRITE (6,*) 'program. It can only perform spin-polarized'
       WRITE (6,*) 'calculations.'
       CALL juDFT_error("jspins not equal 2",calledby="vmatgen")
    ENDIF

    ALLOCATE ( vpw(stars%ng3,3),vis(0:27*stars%mx1*stars%mx2*stars%mx3-1,4),&
         &           vxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,3),vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),&
         &           vz(vacuum%nmzd,2,4),fftwork(0:27*stars%mx1*stars%mx2*stars%mx3-1) )

    !---> reload the spin up and down potential
    !      OPEN (nu,file='pottot',form='unformatted',status='old')
    OPEN (nu,file='nrp',form='unformatted',status='old')
    CALL loddop(stars,vacuum,atoms,sphhar,&
         &            input,sym,&
         &            nu,&
         &            iter,vr,vpw(1,1),vz(1,1,1),vxy(1,1,1,1))
    CLOSE(nu)

    !---> check, whether the direction of magnetic field file exists
    INQUIRE (FILE='dirofmag',EXIST=l_domfexst)
    IF (l_domfexst) THEN
       !--->    if it does, read the theta and phi values
       OPEN (ndomfile,FILE='dirofmag',FORM='unformatted',&
            &        STATUS='unknown')
       READ (ndomfile) (vis(imeshpt,3),imeshpt=0,ifft3-1)
       READ (ndomfile) (vis(imeshpt,4),imeshpt=0,ifft3-1)
       IF (input%film) THEN
          READ (ndomfile) ((vz(imz,ivac,3),imz=vacuum%nmzxyd+1,vacuum%nmzd),&
               &                       ivac=1,vacuum%nvac)
          READ (ndomfile) ((vz(imz,ivac,4),imz=vacuum%nmzxyd+1,vacuum%nmzd),&
               &                       ivac=1,vacuum%nvac)
          READ (ndomfile) (((vvacxy(imeshpt,imz,ivac,3),&
               &                   imeshpt=0,ifft2-1),imz=1,vacuum%nmzxyd),ivac=1,vacuum%nvac)
          READ (ndomfile) (((vvacxy(imeshpt,imz,ivac,4),&
               &                   imeshpt=0,ifft2-1),imz=1,vacuum%nmzxyd),ivac=1,vacuum%nvac)
       ENDIF
       CLOSE (ndomfile)
    ELSE
       !--->    if it doesn't, set all angles to zero
       vis(:,3:4)=0.0
       IF (input%film) THEN
          DO ivac = 1,2
             DO imz = vacuum%nmzxyd+1,vacuum%nmzd
                vz(imz,ivac,3) = 0.0
                vz(imz,ivac,4) = 0.0
             ENDDO
             DO imz = 1,vacuum%nmzxyd
                DO imeshpt = 0,ifft2-1
                   vvacxy(imeshpt,imz,ivac,3) = 0.0
                   vvacxy(imeshpt,imz,ivac,4) = 0.0
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    !---> fouriertransform the spin up and down potential
    !---> in the interstitial, vpw, to real space (vis)
    DO jspin = 1,input%jspins
       CALL fft3d(&
            &               vis(0,jspin),fftwork,&
            &               vpw(1,jspin),&
            &               stars,+1)
    ENDDO

    !---> calculate the four components of the matrix potential on
    !---> real space mesh
    DO imeshpt = 0,ifft3-1
       vup   = vis(imeshpt,1)
       vdown = vis(imeshpt,2)
       theta = vis(imeshpt,3)
       phi   = vis(imeshpt,4)
       !         write (35,'(i4,4f12.6)') mod(imeshpt,33),vup,vdown,theta,phi
       !--->    at first determine the effective potential and magnetic field
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
          vis(imeshpt,ipot) =  vis(imeshpt,ipot) * stars%ufft(imeshpt)
       ENDDO
    ENDDO

    !---> Fouriertransform the matrix potential back to reciprocal space
    DO ipot = 1,2
       fftwork=0.0
       CALL fft3d(&
            &               vis(0,ipot),fftwork,&
            &               vpw(1,ipot),&
            &               stars,-1)
    ENDDO
    CALL fft3d(&
         &           vis(0,3),vis(0,4),&
         &           vpw(1,3),&
         &           stars,-1)

    IF (input%film) THEN
       !--->    fouriertransform the spin up and down potential
       !--->    in the vacuum, vz & vxy, to real space (vvacxy)
       DO jspin = 1,input%jspins
          DO ivac = 1,vacuum%nvac
             DO imz = 1,vacuum%nmzxyd
                vziw = 0.0
                IF (oneD%odi%d1) THEN

                   CALL judft_error("oneD not implemented",calledby="vmatgen")
                   !                  CALL fft2d(&
                   !     &                 oneD%k3,odi%M,odi%n2d,&
                   !     &                 vvacxy(0,imz,ivac,jspin),fftwork,&
                   !     &                 vz(imz,ivac,jspin),vziw,vxy(imz,1,ivac,jspin),&
                   !     &                 vacuum,odi%nq2,odi%kimax2,1,&
                   !     &                  %igf,odl%pgf,odi%nst2)
                ELSE
                   CALL fft2d(&
                        &                 stars,&
                        &                 vvacxy(0,imz,ivac,jspin),fftwork,&
                        &                 vz(imz,ivac,jspin),vziw,vxy(imz,1,ivac,jspin),&
                        &                 vacuum%nmzxyd,1)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       !--->    calculate the four components of the matrix potential on
       !--->    real space mesh
       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             DO imeshpt = 0,ifft2-1
                vup   = vvacxy(imeshpt,imz,ivac,1)
                vdown = vvacxy(imeshpt,imz,ivac,2)
                theta = vvacxy(imeshpt,imz,ivac,3)
                phi   = vvacxy(imeshpt,imz,ivac,4)
                veff  = (vup + vdown)/2.0
                beff  = (vup - vdown)/2.0
                vvacxy(imeshpt,imz,ivac,1) = veff + beff*COS(theta)
                vvacxy(imeshpt,imz,ivac,2) = veff - beff*COS(theta)
                vvacxy(imeshpt,imz,ivac,3) = beff*SIN(theta)*COS(phi)
                vvacxy(imeshpt,imz,ivac,4) = beff*SIN(theta)*SIN(phi)
             ENDDO
          ENDDO
          DO imz = vacuum%nmzxyd+1,vacuum%nmzd
             vup   = vz(imz,ivac,1)
             vdown = vz(imz,ivac,2)
             theta = vz(imz,ivac,3)
             phi   = vz(imz,ivac,4)
             veff  = (vup + vdown)/2.0
             beff  = (vup - vdown)/2.0
             vz(imz,ivac,1) = veff + beff*COS(theta)
             vz(imz,ivac,2) = veff - beff*COS(theta)
             vz(imz,ivac,3) = beff*SIN(theta)*COS(phi)
             vz(imz,ivac,4) = beff*SIN(theta)*SIN(phi)
          ENDDO
       ENDDO

       !--->    Fouriertransform the matrix potential back to reciprocal space
       DO ipot = 1,2
          DO ivac = 1,vacuum%nvac
             DO imz = 1,vacuum%nmzxyd
                fftwork=0.0
                IF (oneD%odi%d1) THEN

                   CALL judft_error("oneD not implemented",calledby="vmatgen")
                   !                CALL fft2d(&
                   !     &                 oneD%k3,odi%M,odi%n2d,&
                   !     &                 vvacxy(0,imz,ivac,ipot),fftwork,&
                   !     &                 vz(imz,ivac,ipot),vziw,vxy(imz,1,ivac,ipot),&
                   !     &                 vacuum,odi%nq2,odi%kimax2,-1,&
                   !     &                  %igf,odl%pgf,odi%nst2)
                ELSE
                   CALL fft2d(&
                        &                 stars,&
                        &                 vvacxy(0,imz,ivac,ipot),fftwork,&
                        &                 vz(imz,ivac,ipot),vziw,vxy(imz,1,ivac,ipot),&
                        &                 vacuum%nmzxyd,-1)
                END IF
             ENDDO
          ENDDO
       ENDDO

       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             fftwork=0.0
             IF (oneD%odi%d1) THEN
                CALL judft_error("oneD not implemented",calledby="vmatgen")
                !              CALL fft2d(&
                !   &              oneD%k3,odi%M,odi%n2d,&
                !   &              vvacxy(0,imz,ivac,3),vvacxy(0,imz,ivac,4),&
                !   &              vz(imz,ivac,3),vz(imz,ivac,4),vxy(imz,1,ivac,3),&
                !   &              vacuum,odi%nq2,odi%kimax2,-1,&
                !   &               %igf,odl%pgf,odi%nst2)
             ELSE
                CALL fft2d(&
                     &              stars,&
                     &              vvacxy(0,imz,ivac,3),vvacxy(0,imz,ivac,4),&
                     &              vz(imz,ivac,3),vz(imz,ivac,4),vxy(imz,1,ivac,3),&
                     &              vacuum%nmzxyd,-1)
             END IF
          ENDDO
       ENDDO

    ENDIF
    !
    !---> save matrix potential to file potmat
    !
    OPEN (npotmatfile,FILE='potmat',FORM='unformatted',&
         &     STATUS='unknown')
    DO ipot = 1,3
       WRITE (npotmatfile) (vpw(ig3,ipot),ig3=1,stars%ng3)
    ENDDO
    IF (input%film) THEN
       DO ivac = 1,vacuum%nvac
          WRITE (npotmatfile)((vz(imz,ivac,ipot),imz=1,vacuum%nmzd),ipot=1,4)
          DO ipot = 1,3
             WRITE (npotmatfile)((vxy(imz,ig2,ivac,ipot),&
                  &                      imz=1,vacuum%nmzxyd),ig2=1,oneD%odi%nq2-1)
          ENDDO
       ENDDO
    ENDIF
8000 FORMAT(6f16.10)
    CLOSE (npotmatfile)

    DEALLOCATE ( vpw,vis,vxy,vr,vz,fftwork)
    IF (input%film) DEALLOCATE (vvacxy)
    RETURN
  END SUBROUTINE vmatgen
END MODULE m_vmatgen
