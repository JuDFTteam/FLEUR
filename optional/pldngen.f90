!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_pldngen

!**********************************************************************
!     This subroutine generates the charge and magetization densities
!     (mx,my,mz) and writes them to the files cdn, mdnx, mdny, mdnz.
!     These files are needed to generate plots of the density.
!
!    i) The components of the hermitian density matrix (rho_11, rho_22,
!     rho_21) are reloaded from the file rhomat_inp.
!    ii) The density matrix in fouriertransformed to real space.
!    iii) The charge and magnetization density (n, mx, my, mz) are
!     calculated on the real space mesh.
!    iv) n, mx, my, and mz are Fouriertransformed and stored in terms
!     of stars.
!
!     Philipp Kurz 99/10/29
!**********************************************************************

CONTAINS

SUBROUTINE pldngen(sym,stars,atoms,sphhar,vacuum,&
                   cell,input,noco,oneD,sliceplot)

          !******** ABBREVIATIONS ***********************************************
          !     ifft3    : size of the 3d real space mesh
          !     ifft2    : size of the 2d real space mesh
          !     rpw      : first diagonal components of the interstitial density
          !                matrix
          !                later charge and mag. density (n, mx, my, mz)
          !                all stored in terms of 3d-stars
          !     ris      : first componets of the density matrix
          !                later charge and mag. density (n, mx, my, mz)
          !                all stored on real space mesh
          !**********************************************************************

   USE m_juDFT
   USE m_constants
   USE m_cdn_io
   USE m_loddop
   USE m_wrtdop
   USE m_qfix
   USE m_fft2d
   USE m_fft3d
   USE m_types
   USE m_rotdenmat 

   IMPLICIT NONE

   TYPE(t_sym),INTENT(IN)    :: sym
   TYPE(t_stars),INTENT(IN)  :: stars
   TYPE(t_vacuum),INTENT(IN) :: vacuum
   TYPE(t_atoms),INTENT(IN)  :: atoms
   TYPE(t_sphhar),INTENT(IN) :: sphhar
   TYPE(t_input),INTENT(IN)  :: input
   TYPE(t_cell),INTENT(IN)   :: cell
   TYPE(t_oneD),INTENT(IN)    :: oneD
   TYPE(t_noco),INTENT(IN)   :: noco
   TYPE(t_sliceplot),INTENT(IN):: sliceplot

   ! Local type instances
   TYPE(t_input)  :: inp
   TYPE(t_potden) :: den

   ! Local Scalars
   INTEGER :: nrhomfile=26   
   INTEGER iden,ivac,ifft2,ifft3
   INTEGER imz,ityp,iri,ilh,imesh,lh,iq2,iq3,iter
   REAL cdnup,cdndown,chden,mgden,theta,phi,zero,rho_11,rziw
   REAL rho_22,rho_21r,rho_21i,rhotot,mx,my,mz,fix,vz_r,vz_i
   COMPLEX czero
   CHARACTER*8 dop,iop,name(10)

   ! Local Arrays
   !---> off-diagonal part of the density matrix
   COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
   COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
   REAL,    ALLOCATABLE :: rht(:,:,:),rho(:,:,:,:)
   REAL,    ALLOCATABLE :: rvacxy(:,:,:,:),ris(:,:),fftwork(:)

   !---> for testing: output of offdiag. output density matrix. to plot the
   !---> offdiag. part of the output density matrix, that part has to be
   !---> written the file rhomt21 in cdnmt.
   LOGICAL :: l_fmpl2
   REAL    :: cdn11, cdn22  
   COMPLEX :: cdn21 
   COMPLEX, ALLOCATABLE :: rho21(:,:,:)
   !---> end of test part

   zero = 0.0 ; czero = CMPLX(0.0,0.0)
   ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
   ifft2 = 9*stars%mx1*stars%mx2

   ALLOCATE (qpw(stars%ng3,4),rhtxy(vacuum%nmzxyd,stars%ng2-1,2,4),&
             cdom(stars%ng3),cdomvz(vacuum%nmzd,2),cdomvxy(vacuum%nmzxyd,stars%ng2-1,2),&
             ris(0:27*stars%mx1*stars%mx2*stars%mx3-1,4),fftwork(0:27*stars%mx1*stars%mx2*stars%mx3-1),&
             rvacxy(0:9*stars%mx1*stars%mx2-1,vacuum%nmzxyd,2,4),&
             rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,4),rht(vacuum%nmzd,2,4) )

   !---> initialize arrays for the density matrix
   rho(:,:,:,:) = zero ; qpw(:,:) = czero ; cdom(:) = czero
   IF (input%film) THEN
      cdomvz(:,:) = czero ;    rhtxy(:,:,:,:) = czero
      cdomvxy(:,:,:) = czero ; rht(:,:,:) = zero
   END IF

   IF (input%jspins .NE. 2) THEN
      WRITE (6,*) 'This is the non-collinear version of the flapw-'
      WRITE (6,*) 'program. It can only perform spin-polarized'
      WRITE (6,*) 'calculations.'
      CALL juDFT_error("jspins not equal 2",calledby = "pldngen",hint=&
                       "This is the non-collinear version of the flapw-"//&
                       ' PROGRAM. It can ONLY perform spin-polarized '//&
                       'calculations.')
   END IF

   !---> reload the density matrix from file rhomat_inp
   IF (.NOT. sliceplot%slice) THEN  
      OPEN (nrhomfile,FILE='rhomat_inp',FORM='unformatted',STATUS='unknown')
   ELSE
      OPEN (nrhomfile,FILE='cdn_slice',FORM='unformatted',STATUS='unknown')
   END IF
   !---> first the diagonal elements of the density matrix
   CALL loddop(stars,vacuum,atoms,sphhar,input,sym,nrhomfile,&
               iter,rho,qpw,rht,rhtxy)
   !---> and then the off-diagonal part
   READ (nrhomfile,END=100,ERR=50) (cdom(iq3),iq3=1,stars%ng3)
   IF (input%film) THEN
      READ (nrhomfile,END=75,ERR=50) ((cdomvz(imz,ivac),imz=1,vacuum%nmz),ivac=1,vacuum%nvac)
      READ (nrhomfile,END=75,ERR=50) (((cdomvxy(imz,iq2-1,ivac),imz=1,vacuum%nmzxy),iq2=2,stars%ng2),ivac=1,vacuum%nvac)
   END IF
   GOTO 150
   50 WRITE(6,*)'rhodirgen: ERROR: Problems while reading density'
   WRITE(6,*)'matrix from file rhomat_inp.'
   CALL juDFT_error("rhomatdir: ERROR while reading file rhomat_inp",calledby ="pldngen")
   75 WRITE(6,*)'rhomatdir: ERROR: reached end of file rhomat_inp'
   WRITE(6,*)'while reading the vacuum part of the off-diagonal'
   WRITE(6,*)'element of the desity matrix.'
   CALL juDFT_error("rhomatdir: ERROR while reading file rhomat_inp",calledby ="pldngen")
   100 WRITE(6,*)'rhodirgen: WARNING: The file rhomat_inp does not'
   WRITE(6,*)'contain off-diagonal part of the density matrix.'
   WRITE(6,*)'Assuming collinear magnetization.'
   150 CLOSE (nrhomfile)
   IF (.NOT. sliceplot%slice) THEN 
      CALL qfix(stars,atoms,sym,vacuum,sphhar,input,cell,oneD,&
                qpw,rhtxy,rho,rht,.FALSE.,.true.,fix)
   END IF

   !---> for testing: read offdiag. output density matrix
   INQUIRE (file= 'rhomt21', exist= l_fmpl2)
   IF (l_fmpl2) THEN
      ALLOCATE( rho21(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) )
      OPEN (26,file='rhomt21',form='unformatted',status='unknown')
      READ (26) rho21
      CLOSE (26)
   END IF
   !---> end of test output

   !---> calculate the charge and magnetization density in the muffin tins
   DO ityp = 1,atoms%ntype
      DO ilh = 0,sphhar%nlh(atoms%ntypsy(ityp))
         DO iri = 1,atoms%jri(ityp)
            IF (.NOT. l_fmpl2) THEN 
               cdnup   = rho(iri,ilh,ityp,1)
               cdndown = rho(iri,ilh,ityp,2)
               theta = noco%beta(ityp)
               phi   = noco%alph(ityp)
               chden  = cdnup + cdndown
               mgden  = cdnup - cdndown
               rho(iri,ilh,ityp,1) = chden
               rho(iri,ilh,ityp,2) = mgden*COS(phi)*SIN(theta)
               rho(iri,ilh,ityp,3) = mgden*SIN(phi)*SIN(theta)
               rho(iri,ilh,ityp,4) = mgden*COS(theta)
            ELSE 
               !--->            for testing: output of offdiag. output density matrix
               cdn11 = rho(iri,ilh,ityp,1)
               cdn22 = rho(iri,ilh,ityp,2)
               cdn21 = rho21(iri,ilh,ityp)
               CALL rot_den_mat(noco%alph(ityp),noco%beta(ityp),cdn11,cdn22,cdn21)
               rho(iri,ilh,ityp,1) = cdn11 + cdn22
               rho(iri,ilh,ityp,2) = 2*REAL(cdn21)
               rho(iri,ilh,ityp,3) = 2*AIMAG(cdn21)
               rho(iri,ilh,ityp,4) = cdn11 - cdn22
               !--->            end of test part
            END IF
         END DO
      END DO
   END DO

   IF (l_fmpl2) THEN
      DEALLOCATE( rho21 )
   END IF

   !---> fouriertransform the diagonal part of the density matrix
   !---> in the interstitial, qpw, to real space (ris)
   DO iden = 1,2
      CALL fft3d(ris(0,iden),fftwork,qpw(1,iden),stars,1)
   END DO
   !---> fouriertransform the off-diagonal part of the density matrix
   CALL fft3d(ris(0,3),ris(0,4),cdom(1),stars,+1)

   !---> calculate the charge and magnetization density on the
   !---> real space mesh
   DO imesh = 0,ifft3-1
      rho_11  = ris(imesh,1)
      rho_22  = ris(imesh,2)
      rho_21r = ris(imesh,3)
      rho_21i = ris(imesh,4)
      rhotot  = rho_11 + rho_22
      mx      =  2*rho_21r
      my      = -2*rho_21i
      mz      = (rho_11-rho_22)

      ris(imesh,1) = rhotot
      ris(imesh,2) = mx
      ris(imesh,3) = my
      ris(imesh,4) = mz
   END DO

   !---> Fouriertransform the density matrix back to reciprocal space
   DO iden = 1,4
      fftwork=zero
      CALL fft3d(ris(0,iden),fftwork,qpw(1,iden),stars,-1)
   END DO

   !---> fouriertransform the diagonal part of the density matrix
   !---> in the vacuum, rz & rxy, to real space (rvacxy)
   IF (input%film) THEN
      DO iden = 1,2
         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               rziw = 0.0
               CALL fft2d(stars,rvacxy(0,imz,ivac,iden),fftwork,rht(imz,ivac,iden),&
                          rziw,rhtxy(imz,1,ivac,iden),vacuum%nmzxyd,1)
            END DO
         END DO
      END DO
      !--->    fouriertransform the off-diagonal part of the density matrix
      DO ivac = 1,vacuum%nvac
         DO imz = 1,vacuum%nmzxyd
            rziw = 0.0
            vz_r = REAL(cdomvz(imz,ivac))
            vz_i = AIMAG(cdomvz(imz,ivac))
            CALL fft2d(stars,rvacxy(0,imz,ivac,3),rvacxy(0,imz,ivac,4),&
                       vz_r,vz_i,cdomvxy(imz,1,ivac),vacuum%nmzxyd,1)
         END DO
      END DO

      !--->    calculate the four components of the matrix potential on
      !--->    real space mesh
      DO ivac = 1,vacuum%nvac
         DO imz = 1,vacuum%nmzxyd
            DO imesh = 0,ifft2-1
               rho_11  = rvacxy(imesh,imz,ivac,1)
               rho_22  = rvacxy(imesh,imz,ivac,2)
               rho_21r = rvacxy(imesh,imz,ivac,3)
               rho_21i = rvacxy(imesh,imz,ivac,4)
               rhotot  = rho_11 + rho_22
               mx      =  2*rho_21r
               my      = -2*rho_21i
               mz      = (rho_11-rho_22)

               rvacxy(imesh,imz,ivac,1) = rhotot
               rvacxy(imesh,imz,ivac,2) = mx
               rvacxy(imesh,imz,ivac,3) = my
               rvacxy(imesh,imz,ivac,4) = mz
            END DO
         END DO
         DO imz = vacuum%nmzxyd+1,vacuum%nmzd
            rho_11  = rht(imz,ivac,1)
            rho_22  = rht(imz,ivac,2)
            rho_21r = REAL(cdomvz(imz,ivac))
            rho_21i = AIMAG(cdomvz(imz,ivac))
            rhotot  = rho_11 + rho_22
            mx      =  2*rho_21r
            my      = -2*rho_21i
            mz      = (rho_11-rho_22)

            rht(imz,ivac,1) = rhotot
            rht(imz,ivac,2) = mx
            rht(imz,ivac,3) = my
            rht(imz,ivac,4) = mz
         END DO
      END DO
      !--->    Fouriertransform the matrix potential back to reciprocal space
      DO iden = 1,4
         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               fftwork=zero
               CALL fft2d(stars,rvacxy(0,imz,ivac,iden),fftwork,rht(imz,ivac,iden),&
                          rziw,rhtxy(imz,1,ivac,iden),vacuum%nmzxyd,-1)
            END DO
         END DO
      END DO
   END IF

   !---> save charge density to file cdn
   inp=input
   inp%jspins=1

   CALL den%init(stars,atoms,sphhar,vacuum,noco,oneD,inp%jspins,.FALSE.,POTDEN_TYPE_DEN)
   den%iter = iter
   den%mt(:,0:,1:,1:1) = rho(:,0:,1:,1:1)
   den%pw(1:,1:1) = qpw(1:,1:1)
   den%vacz(1:,1:,1:1) = rht(1:,1:,1:1)
   den%vacxy(1:,1:,1:,1:1) = rhtxy(1:,1:,1:,1:1)
   den%cdom = cdom
   den%cdomvz = cdomvz
   den%cdomvxy = cdomvxy

   CALL writeDensity(stars,vacuum,atoms,cell,sphhar,inp,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                     0,-1.0,0.0,.FALSE.,den)

   !---> save mx to file mdnx
   OPEN (72,FILE='mdnx',FORM='unformatted',STATUS='unknown')
   CALL wrtdop(stars,vacuum,atoms,sphhar,&
               inp,sym,72,iter,rho(:,0:,1:,2:2),qpw(1:,2:2),rht(1:,1:,2:2),&
               rhtxy(1:,1:,1:,2:2))
   CLOSE (72)

   !---> save my to file mdny
   OPEN (72,FILE='mdny',FORM='unformatted',STATUS='unknown')
   CALL wrtdop(stars,vacuum,atoms,sphhar,&
               inp,sym,72,iter,rho(:,0:,1:,3:3),qpw(1:,3:3),rht(1:,1:,3:3),&
               rhtxy(1:,1:,1:,3:3))
   CLOSE (72)

   !---> save mz to file mdnz
   OPEN (72,FILE='mdnz',FORM='unformatted',STATUS='unknown')
   CALL wrtdop(stars,vacuum,atoms,sphhar,&
               inp,sym,72,iter,rho(:,0:,1:,4:4),qpw(1:,4:4),rht(1:,1:,4:4),&
               rhtxy(1:,1:,1:,4:4))
   CLOSE (72)

   DEALLOCATE (qpw,rhtxy,cdom,cdomvz,cdomvxy,ris,fftwork,rvacxy,rho,rht)

END SUBROUTINE pldngen

END MODULE m_pldngen
