!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_plotdop

use m_juDFT
use m_types
use m_constants,ONLY: pi_const

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

   IMPLICIT NONE

   TYPE(t_oneD),INTENT(IN)     :: oneD
   TYPE(t_dimension),INTENT(IN):: dimension
   TYPE(t_stars),INTENT(IN)    :: stars
   TYPE(t_vacuum),INTENT(IN)   :: vacuum
   TYPE(t_sphhar),INTENT(IN)   :: sphhar
   TYPE(t_atoms),INTENT(IN)    :: atoms
   TYPE(t_input),INTENT(IN)    :: input
   TYPE(t_sym),INTENT(IN)      :: sym
   TYPE(t_cell),INTENT(IN)     :: cell
   TYPE(t_sliceplot),INTENT(IN):: sliceplot
   TYPE(t_noco),INTENT(IN)     :: noco
   CHARACTER(len=10), INTENT (IN) :: cdnfname

!  .. Local Scalars ..
   REAL          :: s,tec,qint,xdnout,fermiEnergyTemp,phi0
   INTEGER       :: i,i1,i2,i3,ii3,ix,iy,iz,jsp,na,nplo
   INTEGER       :: nplot,nq,nt,jm,jspin,iter,ii1,ii2
   LOGICAL       :: twodim,oldform,newform,l_qfix
   LOGICAL       :: cartesian,xsf,unwind

!  .. Local Arrays ..
   COMPLEX :: qpw(stars%ng3,input%jspins),rhtxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
   REAL    :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),rht(vacuum%nmzd,2,input%jspins)
   REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
   COMPLEX :: cdom(1),cdomvz(1,1),cdomvxy(1,1,1)
   INTEGER :: grid(3)
   REAL    :: rhocc(atoms%jmtd)
   REAL    :: point(3)
   CHARACTER (len=30) :: filename
   CHARACTER (len=7)  :: append

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

   ! Read in charge/potential
   IF(TRIM(ADJUSTL(cdnfname)).EQ.'cdn1') THEN
      CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,&
                       CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
   ELSE IF(TRIM(ADJUSTL(cdnfname)).EQ.'cdn') THEN
      CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,&
                       CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,iter,rho,qpw,rht,rhtxy,cdom,cdomvz,cdomvxy)
   ELSE
      OPEN(20,file = cdnfname,form='unformatted',status='old')
      CALL loddop(stars,vacuum,atoms,sphhar,input,sym,20,&
                  iter,rho,qpw,rht,rhtxy)
      CLOSE(20)
   END IF

   ! Perhaps only the core charge should be plotted
   IF (input%score) THEN
      OPEN (17,file='cdnc',form='unformatted',status='old')
      REWIND 17
      DO jspin = 1, input%jspins
         DO nt = 1, atoms%ntype
            jm = atoms%jri(nt)
            READ (17) (rhocc(i),i=1,jm)
            DO i = 1, atoms%jri(nt)
               rho(i,0,nt,jspin) = rho(i,0,nt,jspin) - rhocc(i)/2.0/SQRT(pi_const)
            END DO
            READ (17) tec
         END DO
         READ (17) qint
         qpw(1,jspin) = qpw(1,jspin) - qint/cell%volint
      END DO
      CLOSE (17)
   END IF

   ! Open the plot_inp file for input
   OPEN (18,file='plot_inp')
   READ(18,'(i2,5x,l1)') nplot,xsf
   ! If xsf is specified we create an input file for xcrysden
   IF (xsf) THEN
      IF (noco%l_noco) THEN
         append = '_pl.xsf'
         OPEN (55,file = trim(cdnfname)//append,form='formatted')
      ELSE
         OPEN(55,file="plot.xsf")
      END IF
      CALL xsf_WRITE_atoms(55,atoms,input%film,oneD%odi%d1,cell%amat)
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
      IF (xsf) THEN
         CALL xsf_WRITE_header(55,twodim,filename,vec1,vec2,vec3,zero,grid)
      ELSE
         IF (noco%l_noco) THEN
            OPEN (55,file = filename//cdnfname,form='formatted')
         ELSE
            OPEN (55,file = filename,form='formatted')
         END IF
      END IF

      !loop over spins
      DO jsp = 1, input%jspins
         !loop over all points
         DO iz = 0, grid(3)-1
            DO iy = 0, grid(2)-1
               xloop:DO ix = 0, grid(1)-1

                  point = zero + vec1*REAL(ix)/(grid(1)-1) +&
                                 vec2*REAL(iy)/(grid(2)-1)
                  IF (.NOT.twodim) point = point + vec3*REAL(iz)/(grid(3)-1)
                  !Check if the point is in MT-sphere
                  ii1 = 3
                  ii2 = 3
                  ii3 = 3
                  IF (input%film .AND. .NOT.oneD%odi%d1) ii3 = 0
                  IF (oneD%odi%d1) THEN
                     ii1 = 0 ; ii2 = 0
                  END IF

                  DO i1 = -ii1, ii1
                     DO i2 = -ii2, ii2
                        DO i3 = -ii3, ii3
                           pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
                           na = 0
                           DO nt = 1, atoms%ntype
                              DO nq = 1, atoms%neq(nt)
                                 na   = na + 1
                                 s  = SQRT(dot_PRODUCT(atoms%pos(:,na)-pt,&
                                                       atoms%pos(:,na)-pt))
                                 IF (s<atoms%rmsh(atoms%jri(nt),nt)) THEN
                                    CALL outcdn(pt,nt,na,0,1,jsp,sliceplot,stars,&
                                                vacuum,sphhar,atoms,sym,cell,oneD,&
                                                qpw,rhtxy,rho,rht,xdnout)
                                    IF (xsf) THEN
                                       WRITE(55,*) xdnout
                                    ELSE
                                       WRITE(55,'(4e15.7)') point ,xdnout
                                    END IF
                                    CYCLE xloop
                                 END IF
                              END DO
                           END DO !nt
                        END DO
                     END DO
                  END DO !i1

                  !Check for point in vacuum
                  IF (input%film.AND..NOT.oneD%odi%d1.AND.ABS(point(3))>=cell%z1) THEN
                     CALL outcdn(point,0,0,1,0,jsp,sliceplot,stars,&
                                 vacuum,sphhar,atoms,sym,cell,oneD,&
                                 qpw,rhtxy,rho,rht,xdnout)
                     IF (xsf) THEN
                        WRITE(55,*) xdnout
                     ELSE
                        WRITE(55,'(4e15.7)') point(:), xdnout
                     END IF
                     CYCLE xloop
                  END IF

                  IF (oneD%odi%d1) THEN
                     IF (SQRT((pt(1))**2 + (pt(2))**2)>=cell%z1) THEN
                        CALL outcdn(pt,0,0,1,0,jsp,sliceplot,stars,&
                                    vacuum,sphhar,atoms,sym,cell,oneD,&
                                    qpw,rhtxy,rho,rht,xdnout)
                        IF (xsf) THEN
                           WRITE(55,*) xdnout
                        ELSE
                           WRITE (55,'(4e15.7)') point(:), xdnout
                        END IF
                        CYCLE xloop
                     END IF
                  END IF

                  CALL outcdn(point,0,0,0,2,jsp,sliceplot,stars,&
                              vacuum,sphhar,atoms,sym,cell,oneD,&
                              qpw,rhtxy,rho,rht,xdnout)
                  IF (xsf) THEN
                     WRITE(55,*) xdnout
                  ELSE
                     WRITE(55,'(4e15.7)') point(:), xdnout
                  END IF

               END DO xloop
            END DO
         END DO !z-loop
         IF (xsf.AND.jsp /= input%jspins) &
            CALL xsf_WRITE_newblock(55,twodim,vec1,vec2,vec3,zero,grid)
      END DO !Spin-loop

      IF (xsf) THEN
         CALL xsf_WRITE_endblock(55,twodim)
      ELSE
         CLOSE(55)
      END IF
   END DO !nplot  
    
   CLOSE(18)
   IF (xsf) CLOSE(55)
   END SUBROUTINE plotdop

END MODULE m_plotdop
