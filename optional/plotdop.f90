!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_plotdop

use m_juDFT
use m_types

!+++++++++++++++++++++++++++++++++++++++++++++++++
!  plot the charge density for fleur  code output
!     
!  if twodim = .false. a 3-D plot with nplot layers in z-direction
!  is constructed; the 3x3 matrix gives the 3 vectors of the cell ..
!  .gustav
!
!  Changed the input/output for easier use. 
!  This subroutine uses the file plot_inp for input. 
!  The old plotin-file is still supported by the old subroutine at the 
!  end of the module
!                     Juelich, 21.1.06 DW
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++

CONTAINS

SUBROUTINE plotdop(oneD,dimension,stars,vacuum,sphhar,atoms,&
                   input,sym,cell,sliceplot,noco,cdnfname)

   USE m_outcdn
   USE m_loddop
   USE m_xsf_io
   USE m_cdn_io
   USE m_constants

   IMPLICIT NONE

   TYPE(t_oneD),                INTENT(IN)    :: oneD
   TYPE(t_dimension),           INTENT(IN)    :: dimension
   TYPE(t_stars),               INTENT(IN)    :: stars
   TYPE(t_vacuum),              INTENT(IN)    :: vacuum
   TYPE(t_sphhar),              INTENT(IN)    :: sphhar
   TYPE(t_atoms),               INTENT(IN)    :: atoms
   TYPE(t_input),               INTENT(IN)    :: input
   TYPE(t_sym),                 INTENT(IN)    :: sym
   TYPE(t_cell),                INTENT(IN)    :: cell
   TYPE(t_sliceplot),           INTENT(IN)    :: sliceplot
   TYPE(t_noco),                INTENT(IN)    :: noco
   CHARACTER(len=10), OPTIONAL, INTENT(IN)    :: cdnfname

!  .. Local Scalars ..
   REAL          :: tec,qint,fermiEnergyTemp,phi0,angss
   INTEGER       :: i,j,ix,iy,iz,jsp,na,nplo,iv,iflag,nfile
   INTEGER       :: nplot,nt,jm,jspin,numInFiles,numOutFiles
   LOGICAL       :: twodim,oldform,newform,l_qfix
   LOGICAL       :: cartesian,xsf,unwind,polar

!  .. Local Arrays ..
   TYPE(t_potden), ALLOCATABLE :: den(:)
   REAL, ALLOCATABLE    :: xdnout(:)
   REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3),help(3),qssc(3)
   INTEGER :: grid(3)
   REAL    :: rhocc(atoms%jmtd)
   REAL    :: point(3)
   CHARACTER (len=10), ALLOCATABLE :: cdnFilenames(:)
   CHARACTER (len=15), ALLOCATABLE :: outFilenames(:)
   CHARACTER (len=30)              :: filename
   CHARACTER (len=7)               :: textline

   REAL, PARAMETER :: eps = 1.0e-15

   NAMELIST /plot/twodim,cartesian,unwind,vec1,vec2,vec3,grid,zero,phi0,filename

   oldform = .false.
   INQUIRE(file ="plotin",exist = oldform) 
   IF (oldform) THEN 
      CALL juDFT_error("Use of plotin file no longer supported",calledby = "plotdop")
   END IF

   INQUIRE(file ="plot_inp",exist= newform)
   IF (.NOT.newform) THEN !no input file exists, create a template and exit
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
      WRITE(*,*) "No plot_inp file found. Created a template"
      CALL juDFT_error("Missing input for plot; modify plot_inp",calledby ="plotdop")
   END IF

   nfile = 120
   numInFiles = 0
   numOutFiles = 0
   IF(PRESENT(cdnfname)) THEN
      numInFiles = 1
      numOutFiles = 1
   ELSE
      IF(noco%l_noco) THEN
         numInFiles = 4
         numOutFiles = 4
      ELSE
         numInFiles = 1
         numOutFiles = 1
      END IF
   END IF
   ALLOCATE(den(numInFiles))
   ALLOCATE(cdnFilenames(numInFiles))
   IF(PRESENT(cdnfname)) THEN
      cdnFilenames(1) = cdnfname
   ELSE
      IF(noco%l_noco) THEN
        cdnFilenames(1)='cdn'
        cdnFilenames(2)='mdnx'
        cdnFilenames(3)='mdny'
        cdnFilenames(4)='mdnz'
      ELSE
         IF (sliceplot%slice) THEN
            cdnFilenames(1)='cdn_slice'
         ELSE
            cdnFilenames(1)='cdn1'
         END IF
      END IF
   END IF

   ! Read in charge/potential
   DO i = 1, numInFiles
      CALL den(i)%init(stars,atoms,sphhar,vacuum,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)
      IF(TRIM(ADJUSTL(cdnFilenames(i))).EQ.'cdn1') THEN
         CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,&
                          CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den(i))
      ELSE IF(TRIM(ADJUSTL(cdnFilenames(i))).EQ.'cdn') THEN
         CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,&
                          CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den(i))
      ELSE
         CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,&
                          CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den(i),TRIM(ADJUSTL(cdnFilenames(i))))
      END IF

      ! Subtract core charge if input%score is set
      IF ((.NOT.noco%l_noco).AND.(input%score)) THEN
         OPEN (17,file='cdnc',form='unformatted',status='old')
         REWIND 17
         DO jspin = 1, input%jspins
            DO nt = 1, atoms%ntype
               jm = atoms%jri(nt)
               READ (17) (rhocc(j),j=1,jm)
               DO j = 1, atoms%jri(nt)
                  den(i)%mt(j,0,nt,jspin) = den(i)%mt(j,0,nt,jspin) - rhocc(j)/2.0/SQRT(pi_const)
               END DO
               READ (17) tec
            END DO
            READ (17) qint
            den(i)%pw(1,jspin) = den(i)%pw(1,jspin) - qint/cell%volint
         END DO
         CLOSE (17)
      ELSE IF (input%score) THEN
         CALL juDFT_error('Subtracting core charge in noco calculations not supported', calledby = 'plotdop')
      END IF
   END DO

   IF (noco%l_ss) THEN 
      qssc = MATMUL(TRANSPOSE(cell%bmat),noco%qss) 
   END IF 

   ! Open the plot_inp file for input
   OPEN (18,file='plot_inp')
   READ(18,'(i2,5x,l1,a)') nplot,xsf,textline
   polar = .FALSE.
   IF ((noco%l_noco).AND.(numInFiles.EQ.4)) THEN
      polar = (textline(1:7)=='polar=T').OR.(textline(1:7)=='polar=t')
      IF (polar) THEN
         numOutFiles = 7
      END IF
   END IF
   ALLOCATE(outFilenames(numOutFiles))
   ALLOCATE(xdnout(numOutFiles))
   IF(numOutFiles.EQ.1) THEN
      outFilenames(1) = 'plot'
   ELSE
      DO i = 1, numInFiles
         outFilenames(i) = cdnFilenames(i)//'_pl'
      END DO
      IF (polar) THEN
         outFilenames(5) = 'mabs_pl'
         outFilenames(6) = 'mtha_pl'
         outFilenames(7) = 'mphi_pl'
      END IF
   END IF

   ! If xsf is specified we create input files for xcrysden
   IF (xsf) THEN
      DO i = 1, numOutFiles
         OPEN(nfile+i,file=TRIM(ADJUSTL(outFilenames(i)))//'.xsf',form='formatted')
         CALL xsf_WRITE_atoms(nfile+i,atoms,input%film,oneD%odi%d1,cell%amat)
      END DO
   END IF

   ! Loop over all plots
   DO nplo = 1, nplot

      ! the defaults
      twodim = .TRUE.
      cartesian = .TRUE.
      grid = (/100,100,100/)
      vec1 = (/0.,0.,0./)
      vec2 = (/0.,0.,0./)
      vec3 = (/0.,0.,0./)
      zero = (/0.,0.,0./)
      filename = "default"
      READ(18,plot)
      IF (twodim.AND.ANY(grid(1:2)<1)) &
         CALL juDFT_error("Illegal grid size in plot",calledby="plotdop")
      IF (.NOT.twodim.AND.ANY(grid<1)) &
         CALL juDFT_error("Illegal grid size in plot",calledby="plotdop")
      IF (twodim) grid(3) = 1

      !calculate cartesian coordinates if needed
      IF (.NOT.cartesian) THEN
         vec1=matmul(cell%amat,vec1)
         vec2=matmul(cell%amat,vec2)
         vec3=matmul(cell%amat,vec3)
         zero=matmul(cell%amat,zero)
      END IF

      !Open the file
      IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo
      DO i = 1, numOutFiles
         IF (xsf) THEN
            CALL xsf_WRITE_header(nfile+i,twodim,filename,vec1,vec2,vec3,zero,grid)
         ELSE
            IF (numOutFiles.NE.1) THEN
               OPEN (nfile+i,file = filename//outFilenames(i),form='formatted')
            ELSE
               OPEN (nfile+i,file = filename,form='formatted')
            END IF
         END IF
      END DO

      !loop over spins
      DO jsp = 1, input%jspins
         !loop over all points
         DO iz = 0, grid(3)-1
            DO iy = 0, grid(2)-1
               DO ix = 0, grid(1)-1

                  point = zero + vec1*REAL(ix)/(grid(1)-1) +&
                                 vec2*REAL(iy)/(grid(2)-1)
                  IF (.NOT.twodim) point = point + vec3*REAL(iz)/(grid(3)-1)

                  ! Set region specific parameters for point
                  
                  ! Get MT sphere for point if point is in MT sphere
                  CALL getMTSphere(input,cell,atoms,oneD,point,nt,na,pt)
                  IF (na.NE.0) THEN
                     ! In MT sphere
                     iv = 0
                     iflag = 1
                  ELSE IF (input%film.AND..NOT.oneD%odi%d1.AND.ABS(point(3))>=cell%z1) THEN
                     ! In vacuum in 2D system
                     iv = 1
                     iflag = 0
                     pt(:) = point(:)
                  ELSE IF ((oneD%odi%d1).AND.(SQRT((point(1))**2 + (point(2))**2)>=cell%z1)) THEN
                     ! In vacuum in 1D system
                     iv = 1
                     iflag = 0
                     pt(:) = point(:)
                  ELSE
                     ! In interstitial region
                     iv = 0
                     iflag = 2
                     pt(:) = point(:)
                  END IF

                  DO i = 1, numInFiles
                     CALL outcdn(pt,nt,na,iv,iflag,jsp,sliceplot,stars,&
                                 vacuum,sphhar,atoms,sym,cell,oneD,&
                                 den(i)%pw,den(i)%vacxy,den(i)%mt,&
                                 den(i)%vacz,xdnout(i))
                  END DO

                  IF (na.NE.0) THEN
                     IF (noco%l_ss) THEN 
                        ! rotate magnetization "backward"
                        angss = DOT_PRODUCT(qssc,pt-atoms%pos(:,na))
                        help(1) = xdnout(2)
                        help(2) = xdnout(3)
                        xdnout(2) = +help(1)*COS(angss)+help(2)*SIN(angss) 
                        xdnout(3) = -help(1)*SIN(angss)+help(2)*COS(angss) 
                        ! xdnout(2)=0. ; xdnout(3)=0. ; xdnout(4)=0. 
                     END IF
                  END IF

                  IF (noco%l_ss .AND. (.NOT. unwind)) THEN
                     ! rotate magnetization
                     angss = DOT_PRODUCT(qssc,point)
                     help(1) = xdnout(2)
                     help(2) = xdnout(3)
                     xdnout(2) = +help(1)*COS(angss) -help(2)*SIN(angss)
                     xdnout(3) = +help(1)*SIN(angss) +help(2)*COS(angss)
                  END IF

                  IF (polar) THEN
                     xdnout(5) = SQRT(ABS(xdnout(2)**2+xdnout(3)**2+xdnout(4)**2))
                     IF (xdnout(5)<eps) THEN
                        xdnout(5)= 0.0
                        xdnout(6)= -tpi_const
                        xdnout(7)= -tpi_const
                     ELSE
                        DO j = 1, 3
                           help(j) = xdnout(1+j)/xdnout(5) 
                        END DO
                        IF (help(3)<0.5) THEN
                           xdnout(6)= ACOS(help(3))
                        ELSE
                           xdnout(6)= pi_const/2.0-ASIN(help(3))
                        END IF
                        IF (SQRT(ABS(help(1)**2+help(2)**2)) < eps) THEN
                           xdnout(7)= -tpi_const
                        ELSE
                           IF ( ABS(help(1)) > ABS(help(2)) ) THEN
                              xdnout(7)= ABS(ATAN(help(2)/help(1)))
                           ELSE
                              xdnout(7)= pi_const/2.0-ABS(ATAN(help(1)/help(2)))
                           END IF
                           IF (help(2)<0.0) THEN
                              xdnout(7)= -xdnout(7)
                           END IF
                           IF (help(1)<0.0) THEN
                              xdnout(7)= pi_const-xdnout(7)
                           END IF
                           DO WHILE (xdnout(7)-pi_const*phi0 > +pi_const)
                              xdnout(7)= xdnout(7)-tpi_const
                           END DO
                           DO WHILE (xdnout(7)-pi_const*phi0 < -pi_const)
                              xdnout(7)= xdnout(7)+tpi_const
                           END DO
                        END IF
                     END IF
                     xdnout(6)= xdnout(6)/pi_const
                     xdnout(7)= xdnout(7)/pi_const
                  END IF ! (polar)

                  DO i = 1, numOutFiles
                     IF (xsf) THEN
                        WRITE(nfile+i,*) xdnout(i)
                     ELSE
                        WRITE(nfile+i,'(4e15.7)') point ,xdnout(i)
                     END IF
                  END DO

               END DO
            END DO
         END DO !z-loop
         DO i = 1, numOutFiles
            IF (xsf.AND.jsp /= input%jspins) &
               CALL xsf_WRITE_newblock(nfile+i,twodim,vec1,vec2,vec3,zero,grid)
         END DO
      END DO !Spin-loop

      DO i = 1, numOutFiles
         IF (xsf) THEN
            CALL xsf_WRITE_endblock(nfile+i,twodim)
         ELSE
            CLOSE(nfile+i)
         END IF
      END DO
   END DO !nplot  
    
   CLOSE(18)
   IF (xsf) THEN
      DO i = 1, numOutFiles
         CLOSE(nfile+i)
      END DO
   END IF

   DEALLOCATE(xdnout, cdnFilenames, outFilenames)

   END SUBROUTINE plotdop



   SUBROUTINE getMTSphere(input,cell,atoms,oneD,point,iType,iAtom,pt)

   IMPLICIT NONE

   TYPE(t_input), INTENT(IN)    :: input
   TYPE(t_cell),  INTENT(IN)    :: cell
   TYPE(t_atoms), INTENT(IN)    :: atoms
   TYPE(t_oneD),  INTENT(IN)    :: oneD

   INTEGER,       INTENT(OUT)   :: iType, iAtom
   REAL,          INTENT(OUT)   :: pt(3)
   REAL,          INTENT(IN)    :: point(3)

   INTEGER                      :: ii1, ii2, ii3, i1, i2, i3, nq
   REAL                         :: s

   ii1 = 3
   ii2 = 3
   ii3 = 3
   IF (input%film .AND. .NOT.oneD%odi%d1) ii3 = 0
   IF (oneD%odi%d1) THEN
      ii1 = 0
      ii2 = 0
   END IF

   DO i1 = -ii1, ii1
      DO i2 = -ii2, ii2
         DO i3 = -ii3, ii3
            pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
            iAtom = 0
            DO iType = 1, atoms%ntype
               DO nq = 1, atoms%neq(iType)
                  iAtom = iAtom + 1
                  s = SQRT(DOT_PRODUCT(atoms%pos(:,iAtom)-pt,atoms%pos(:,iAtom)-pt))
                  IF (s<atoms%rmsh(atoms%jri(iType),iType)) THEN
                     ! Return with the current iType, iAtom, pt
                     RETURN
                  END IF
               END DO
            END DO
         END DO
      END DO
   END DO !i1

   ! If no MT sphere encloses the point return 0 for iType, iAtom
   iType = 0
   iAtom = 0
   pt(:) = point(:)

   END SUBROUTINE getMTSphere

END MODULE m_plotdop
