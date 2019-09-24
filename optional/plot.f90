!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_plot
   USE m_types
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

   PRIVATE
   !------------------------------------------------
   ! A general purpose plotting routine for FLEUR.
   ! 
   ! Based on older plotting routines in pldngen.f90
   ! and plotdop.f90 called by optional.F90 and now
   ! called within a scf-loop instead of as a post
   ! process functionality.
   ! 
   ! A. Neukirchen, September 2019 
   !------------------------------------------------
   INTEGER, PARAMETER :: PLOT_INPDEN_const=2        !ind_plot= 1
   INTEGER, PARAMETER :: PLOT_OUTDEN_Y_CORE_const=4 !ind_plot= 2
   INTEGER, PARAMETER :: PLOT_INPDEN_N_CORE_const=8 !ind_plot= 4
   INTEGER, PARAMETER :: PLOT_POT_TOT_const=128     !ind_plot= 7
   INTEGER, PARAMETER :: PLOT_POT_EXT_const=256     !ind_plot= 8
   INTEGER, PARAMETER :: PLOT_POT_COU_const=512     !ind_plot= 9
   INTEGER, PARAMETER     :: PLOT_POT_VXC_const=1024    !ind_plot=10
   !CHARACTER, PARAMETER, DIMENSION(7) :: filenames=(/''/)
   PUBLIC             :: checkplotinp, genplotinp, doplots, vectorsplit, matrixsplit, scalarplot, vectorplot, matrixplot

CONTAINS

   SUBROUTINE checkplotinp()
      ! Checks for existing plot input. If an ancient plotin file is used, an
      ! error is called. If no usable plot_inp exists, a new one is generated. 

      oldform = .false.
      INQUIRE(file = "plotin", exist = oldform) 
      IF (oldform) THEN 
         CALL juDFT_error("Use of plotin file no longer supported",calledby = "plotdop")
      END IF

      INQUIRE(file = "plot_inp", exist = newform)
      IF (.NOT.newform) THEN CALL genplotinp()

   END SUBROUTINE checkplotinp

   SUBROUTINE genplotinp()
   ! Generates the necessary plot_inp file that dictates the parameters of all
   ! generated .xsf files.

      OPEN(20,file ="plot_inp")
      WRITE(20,'(i2,a5,l1)') 2,",xsf=",.true.
      WRITE(20,*) "&PLOT twodim=t,cartesian=t"
      WRITE(20,*) "  vec1(1)=10.0 vec2(2)=10.0"
      WRITE(20,*) "  filename='plot1' /"
      WRITE(20,*) "&PLOT twodim=f,cartesian=f"
      WRITE(20,*) "  vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0 "
      WRITE(20,*) "  vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0 "
      WRITE(20,*) "  vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0 "
      WRITE(20,*) "  grid(1)=30  grid(2)=30  grid(3)=30  "
      WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.5 "
      WRITE(20,*) "  filename ='plot2' /"
      CLOSE(20)

   END SUBROUTINE genplotinp
   
   SUBROUTINE doplots(jspins,noco,iplot,ind_plot,den)   
      INTEGER, INTENT(IN) :: iplot
      INTEGER, INTENT(IN) :: ind_plot !Index of the plot according to the constants set above
      INTEGER :: jplot
      
      IF btest(iplot,ind_plot) THEN
         jplot=2**ind_plot   
         CALL checkplotinp()
         CALL doplot(jspins,noco,jplot,den)
      END IF
   END SUBROUTINE doplots

   SUBROUTINE doplot(jspins,noco,iplot,ind_plot,den)   
      INTEGER, INTENT(IN) :: jplot
      CHARACTER (len=15), ALLOCATABLE :: outFilenames(:)
      INTEGER :: i
      
      ! Plotting the density matrix as n or n,m or n,mx,my,mz 
      IF jplot==2 THEN
         IF jspins==2 THEN
            IF noco%l_noco THEN
               ALLOCATE(outFilenames(4))
               outFilenames(1)='cden'
               outFilenames(2)='mdnx'
               outFilenames(3)='mdny'
               outFilenames(4)='mdnz'
               !--> matrixplot mit factor=1
            ELSE
               ALLOCATE(outFilenames(2))
               outFilenames(1)='cden'
               outFilenames(2)='mden'
            END IF
         ELSE
            ALLOCATE(outFilenames(1))
            outFilenames(1)='cden'
         END IF
      END IF
   END SUBROUTINE doplot
   
   
   SUBROUTINE vectorsplit(den,den1,den2)
   ! Takes a 2D potential/density vector and rearanges it into two plottable
   ! seperate ones (e.g. [rho_up, rho_down] ---> n, m).
      TYPE(t_potden) :: den1
      CALL den1%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)

   END SUBROUTINE vectorsplit

   SUBROUTINE matrixsplit(mpi,sym,stars,atoms,sphhar,vacuum,cell,input,noco,oneD,sliceplot,factor,denmat,cden,mxden,myden,mzden)
   ! Takes a 2x2 potential/density matrix and rearanges it into four plottable
   ! seperate ones (e.g. rho_mat ---> n, mx, my, mz).

      IMPLICIT NONE

      TYPE(t_mpi),       INTENT(IN)    :: mpi
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_sliceplot), INTENT(IN)    :: sliceplot
      REAL,              INTENT(IN)    :: factor
      TYPE(t_potden),    INTENT(INOUT) :: denmat
      TYPE(t_potden),    INTENT(OUT)   :: cden, mxden, myden, mzden

      ! Local type instances
      TYPE(t_input)  :: inp
      TYPE(t_potden) :: den

      ! Local scalars
      INTEGER iden,ivac,ifft2,ifft3
      INTEGER imz,ityp,iri,ilh,imesh,iter
      REAL cdnup,cdndown,chden,mgden,theta,phi,zero,rho_11,rziw,fermiEnergyTemp
      REAL rho_22,rho_21r,rho_21i,rhotot,mx,my,mz,fix,vz_r,vz_i
      COMPLEX czero

      ! Local arrays
      !---> off-diagonal part of the density matrix
      COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
      COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
      REAL,    ALLOCATABLE :: rht(:,:,:),rho(:,:,:,:)
      REAL,    ALLOCATABLE :: rvacxy(:,:,:,:),ris(:,:),fftwork(:)

      !---> for testing: output of offdiag. output density matrix. to plot the
      !---> offdiag. part of the output density matrix, that part has to be
      !---> written the file rhomt21 in cdnmt.
      LOGICAL :: l_qfix
      REAL    :: cdn11, cdn22
      COMPLEX :: cdn21
      !---> end of test part

      iter = 0 ! This is not clean!

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

      ! Save the density matrix to a work density
      CALL den%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
      den=denmat

      rho(:,0:,1:,:input%jspins) = factor*den%mt(:,0:,1:,:input%jspins)
      qpw(1:,:input%jspins) = factor*den%pw(1:,:input%jspins)
      rht(1:,1:,:input%jspins) = factor*den%vacz(1:,1:,:input%jspins)
      rhtxy(1:,1:,1:,:input%jspins) = factor*den%vacxy(1:,1:,1:,:input%jspins)
      IF(noco%l_noco) THEN
         cdom = factor*den%pw(:,3)
         cdomvz(:,:) = CMPLX(factor*den%vacz(:,:,3),factor*den%vacz(:,:,4))
         cdomvxy = factor*den%vacxy(:,:,:,3)
      END IF

      IF (.NOT. sliceplot%slice) THEN
         CALL den%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
         den%iter = iter
         den%mt(:,0:,1:,:input%jspins) = rho(:,0:,1:,:input%jspins)
         den%pw(1:,:input%jspins) = qpw(1:,:input%jspins)
         den%vacz(1:,1:,:input%jspins) = rht(1:,1:,:input%jspins)
         den%vacxy(1:,1:,1:,:input%jspins) = rhtxy(1:,1:,1:,:input%jspins)
         IF(noco%l_noco) THEN
            den%pw(:,3) = cdom
            den%vacz(:,:,3) = REAL(cdomvz(:,:))
            den%vacz(:,:,4) = AIMAG(cdomvz(:,:))
            den%vacxy(:,:,:,3) = cdomvxy
         END IF
         CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,den,noco%l_noco,.FALSE.,.true.,fix)
         rho(:,0:,1:,:input%jspins) = den%mt(:,0:,1:,:input%jspins)
         qpw(1:,:input%jspins) = den%pw(1:,:input%jspins)
         rht(1:,1:,:input%jspins) = den%vacz(1:,1:,:input%jspins)
         rhtxy(1:,1:,1:,:input%jspins) = den%vacxy(1:,1:,1:,:input%jspins)
         IF(noco%l_noco) THEN
            cdom = den%pw(:,3)
            cdomvz(:,:) = CMPLX(den%vacz(:,:,3),den%vacz(:,:,4))
            cdomvxy = den%vacxy(:,:,:,3)
         END IF
      END IF
   
      !---> calculate the charge and magnetization density in the muffin tins
      DO ityp = 1,atoms%ntype
         DO ilh = 0,sphhar%nlh(atoms%ntypsy(ityp))
            DO iri = 1,atoms%jri(ityp)
               IF (SIZE(den%mt,4).LE.2) THEN 
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
                  cdn21 = CMPLX(den%mt(iri,ilh,ityp,3),den%mt(iri,ilh,ityp,4))
                  CALL rot_den_mat(noco%alph(ityp),noco%beta(ityp),cdn11,cdn22,cdn21)
                  rho(iri,ilh,ityp,1) = cdn11 + cdn22
                  rho(iri,ilh,ityp,2) = 2.0*REAL(cdn21)
                  ! Note: The minus sign in the following line is temporary to adjust for differences in the offdiagonal
                  !       part of the density between this fleur version and ancient (v0.26) fleur.
                  rho(iri,ilh,ityp,3) = -2.0*AIMAG(cdn21)
                  rho(iri,ilh,ityp,4) = cdn11 - cdn22
                  !--->            end of test part
               END IF
            END DO
         END DO
      END DO


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

      cden=den

      !---> save mx to file mdnx
      den%mt(:,0:,1:,1) = rho(:,0:,1:,2)
      den%pw(1:,1) = qpw(1:,2)
      den%vacz(1:,1:,1) = rht(1:,1:,2)
      den%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,2)

      mxden=den

      !---> save my to file mdny
      den%mt(:,0:,1:,1) = rho(:,0:,1:,3)
      den%pw(1:,1) = qpw(1:,3)
      den%vacz(1:,1:,1) = rht(1:,1:,3)
      den%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,3)

      myden=den
   
      !---> save mz to file mdnz
      den%mt(:,0:,1:,1) = rho(:,0:,1:,4)
      den%pw(1:,1) = qpw(1:,4)
      den%vacz(1:,1:,1) = rht(1:,1:,4)
      den%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,4)

      mzden=den

      DEALLOCATE (qpw,rhtxy,cdom,cdomvz,cdomvxy,ris,fftwork,rvacxy,rho,rht)

   END SUBROUTINE matrixsplit

   SUBROUTINE scalarplot(iplot,den,filename)
   !Takes a 1-component t_potden density, i.e. a scalar field in MT-sphere/star
   !representation and makes it into a plottable .xsf file according to a scheme
   !given in plot_inp.
   
   END SUBROUTINE scalarplot

   SUBROUTINE vectorplot()
   !Takes a spin-polarized t_potden density, i.e. a 2D vector in MT-sphere/star
   !representation and makes it into a plottable .xsf file according to a scheme
   !given in plot_inp.

      CALL vectorsplit(den,den1,den2)
      CALL scalarplot(den1,filenames(1))
      CALL scalarplot(den2,filenames(1))

   END SUBROUTINE vectorplot

   SUBROUTINE matrixplot(iplot,den,filenames,factor)
   !Takes a 2x2 t_potden density, i.e. a sum of Pauli matrices in MT-sphere/star
   !representation and makes it into 4 plottable .xsf files according to a scheme
   !given in plot_inp.
   
      !Local variables:

      CALL matrixsplit(den,den1,den2,den3,den4,factor)
      CALL scalarplot(den1,filenames(1))
      CALL scalarplot(den2,filenames(1))
      CALL scalarplot(den3,filenames(1))
      CALL scalarplot(den4,filenames(1))
      
   END SUBROUTINE matrixplot

END MODULE m_plot
