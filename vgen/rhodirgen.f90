!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rhodirgen
  USE m_juDFT
  !**********************************************************************
  !     This subroutine calculates the spin-up and -down density, which
  !     are needed to calculate the potential and writes them to the file
  !     cdn. The local angle of the magnetization is kept in real space
  !     and written to the file dirofmag. This is done in four steps.
  !
  !    i) The components of the hermitian density matrix (rho_11, rho_22,
  !     rho_21) are reloaded from the file rhomat_inp.
  !    ii) The density matrix in fouriertransformed to real space.
  !    iii) The spin-up and -down densities and the local angle of the
  !     magnetization are calculated on the real space mesh.    
  !    iv) The spin-up and -down densities are Fouriertransformed, stored
  !     in terms of stars and written to the file cdn. The local angle of
  !     magnetization is kept on the real space mesh and written to the
  !     file dirofmag.
  !
  !     Philipp Kurz 99/11/01
  !**********************************************************************
CONTAINS
  SUBROUTINE rhodirgen(DIMENSION,sym,stars,atoms,sphhar,vacuum,&
                       ndomfile,cell,input,noco,oneD)

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
    USE m_cdn_io
    USE m_qfix
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

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ndomfile    
    !     ..
    !     ..
    !-odim
    !+odim

    ! Local type instances
    TYPE(t_potden) :: den

    !     .. Local Scalars ..
    INTEGER iden,jspin,ivac,ifft2,ifft3
    INTEGER imz,ityp,iri,ilh,imesh,iq2,iq3
    REAL   zero,rho_11,rho_22,rho_21r,rho_21i,rhotot,magmom,phi
    REAL rho_up,rho_down,mx,my,mz,eps,pi,fix,vz_r,vz_i,rziw,theta
    REAL fermiEnergyTemp
    COMPLEX czero
    CHARACTER*8 dop,iop,name(10)
    LOGICAL l_qfix
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: rz(:,:,:)
    REAL,    ALLOCATABLE :: rvacxy(:,:,:,:),ris(:,:),fftwork(:)
    !     ..
    zero = 0.0 ; czero = CMPLX(0.0,0.0) 
    eps = 1.0e-20
  
    CALL den%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,.FALSE.,POTDEN_TYPE_DEN)

    ALLOCATE (ris(0:27*stars%mx1*stars%mx2*stars%mx3-1,4),fftwork(0:27*stars%mx1*stars%mx2*stars%mx3-1),&
              rz(vacuum%nmzd,2,2))
    !
    !---> initialize arrays for the density matrix
    !
    den%cdom(:) = czero
    den%cdomvz(:,:) = czero
    den%cdomvxy(:,:,:) = czero
    IF (input%film) THEN
       rz(:,:,:) = zero
    ENDIF

    ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
    ifft2 = 9*stars%mx1*stars%mx2
    IF (oneD%odi%d1) ifft2 = 9*stars%mx3*oneD%odi%M
    IF (input%film) ALLOCATE(rvacxy(0:ifft2-1,vacuum%nmzxyd,2,4))

    IF (input%jspins .NE. 2) THEN
       WRITE (6,*) 'This is the non-collinear version of the flapw-'
       WRITE (6,*) 'program. It can only perform spin-polarized'
       WRITE (6,*) 'calculations.'
       CALL juDFT_error("jspins not equal 2",calledby="rhodirgen")
    ENDIF

    !---> reload the density matrix from file rhomat_inp
    CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_NOCO_const,CDN_INPUT_DEN_const,&
                     0,fermiEnergyTemp,l_qfix,den)

    CALL qfix(stars,atoms,sym,vacuum,sphhar,input,cell,oneD,&
              den%pw,den%vacxy,den%mt,den%vacz,.FALSE.,.FALSE.,fix)

    !---> fouriertransform the diagonal part of the density matrix
    !---> in the interstitial, den%pw, to real space (ris)
    DO iden = 1,2
       CALL fft3d(ris(0,iden),fftwork,den%pw(1,iden),stars,+1)
    ENDDO
    !---> fouriertransform the off-diagonal part of the density matrix
    CALL fft3d(ris(0,3),ris(0,4),den%cdom(1),stars,+1)

    !test
    !      DO iden=1,4
    !         write(*,*)'iden=',iden
    !         write(*,8500)(ris(imesh,iden),imesh=0,ifft3-1)
    !      enddo
    !test
    !---> calculate the charge and magnetization density on the
    !---> real space mesh
    DO imesh = 0,ifft3-1
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
       ris(imesh,3) = theta
       ris(imesh,4) = phi
    ENDDO
    !test
    !      DO iden=1,4
    !         write(*,*)'iden=',iden
    !         write(*,8500)(ris(imesh,iden),imesh=0,ifft3-1)
    ! 8500    format(10e13.5)
    !      enddo
    !test
    !---> Fouriertransform the density matrix back to reciprocal space
    DO jspin = 1,input%jspins
       fftwork=0.0
       CALL fft3d(ris(0,jspin),fftwork,den%pw(1,jspin),stars,-1)
    ENDDO

    !---> fouriertransform the diagonal part of the density matrix
    !---> in the vacuum, rz & rxy, to real space (rvacxy)
    IF (input%film) THEN
       DO iden = 1,2
          DO ivac = 1,vacuum%nvac
             DO imz = 1,vacuum%nmzxyd
                rziw = 0.0
                IF (oneD%odi%d1) THEN
                   call judft_error("oneD not implemented",calledby="rhodirgen")
                   !CALL fft2d(oneD%k3,odi%M,odi%n2d,rvacxy(0,imz,ivac,iden),fftwork,&
                   !           den%vacz(imz,ivac,iden),rziw,den%vacxy(imz,1,ivac,iden),&
                   !           vacuum,odi%nq2,odi%kimax2,1,&
                   !     &                  %igf,odl%pgf,odi%nst2)
                ELSE
                   CALL fft2d(stars,rvacxy(0,imz,ivac,iden),fftwork,&
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
             vz_r = REAL(den%cdomvz(imz,ivac))
             vz_i = AIMAG(den%cdomvz(imz,ivac))
             IF (oneD%odi%d1) THEN
                   call judft_error("oneD not implemented",calledby="rhodirgen")
                !CALL fft2d(oneD%k3,odi%M,odi%n2d,&
                !           rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),&
                !           vz_r,vz_i,den%cdomvxy(imz,1,ivac),&
                !           vacuum,odi%nq2,odi%kimax2,1,&
                !     &               %igf,odl%pgf,odi%nst2)
             ELSE
                CALL fft2d(stars,rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),&
                           vz_r,vz_i,den%cdomvxy(imz,1,ivac),vacuum%nmzxyd,1)
             ENDIF
          ENDDO
       ENDDO

       !--->    calculate the four components of the matrix potential on
       !--->    real space mesh
       DO ivac = 1,vacuum%nvac
          DO imz = 1,vacuum%nmzxyd
             DO imesh = 0,ifft2-1
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
                rvacxy(imesh,imz,ivac,3) = theta
                rvacxy(imesh,imz,ivac,4) = phi
             ENDDO
          ENDDO
          DO imz = vacuum%nmzxyd+1,vacuum%nmzd
             rho_11  = den%vacz(imz,ivac,1)
             rho_22  = den%vacz(imz,ivac,2)
             rho_21r = REAL(den%cdomvz(imz,ivac))
             rho_21i = AIMAG(den%cdomvz(imz,ivac))
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
             rz(imz,ivac,1) = theta
             rz(imz,ivac,2) = phi
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
                   CALL fft2d(stars,rvacxy(0,imz,ivac,jspin),fftwork,&
                              den%vacz(imz,ivac,jspin),rziw,den%vacxy(imz,1,ivac,jspin),&
                              vacuum%nmzxyd,-1)
                END IF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !---> ndomfile is the dirofmag-file
    OPEN (ndomfile,FILE='dirofmag',FORM='unformatted',STATUS='unknown')
    WRITE (ndomfile) (ris(imesh,3),imesh=0,ifft3-1)
    WRITE (ndomfile) (ris(imesh,4),imesh=0,ifft3-1)
    IF (input%film) THEN
       WRITE (ndomfile) ((rz(imz,ivac,1),imz=vacuum%nmzxyd+1,vacuum%nmzd),ivac=1,vacuum%nvac)
       WRITE (ndomfile) ((rz(imz,ivac,2),imz=vacuum%nmzxyd+1,vacuum%nmzd),ivac=1,vacuum%nvac)
       WRITE (ndomfile) (((rvacxy(imesh,imz,ivac,3),imesh=0,ifft2-1),imz=1,vacuum%nmzxyd),ivac=1,vacuum%nvac)
       WRITE (ndomfile) (((rvacxy(imesh,imz,ivac,4),imesh=0,ifft2-1),imz=1,vacuum%nmzxyd),ivac=1,vacuum%nvac)
    ENDIF
    CLOSE (ndomfile)

    !---> write spin-up and -down density on file cdn

    CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                      0,-1.0,0.0,.FALSE.,den)

    DEALLOCATE (ris,fftwork,rz)
    IF (input%film) DEALLOCATE(rvacxy)
    RETURN
  END SUBROUTINE rhodirgen
END MODULE m_rhodirgen
