!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_plotdop
      use m_juDFT
      use m_types
      use m_constants,ONLY: pi_const
!     +++++++++++++++++++++++++++++++++++++++++++++++++
!     plot the charge density for fleur  code output
!     
!     if twodim = .false. a 3-D plot with nplot layers in z-direction
!     is constructed; the 3x3 matrix gives the 3 vectors of the cell ..
!     .gustav
!
!    Changed the input/output for easier use. 
!    This subroutine uses the file plot_inp for input. 
!    The old plotin-file is still supported by the old subroutine at the end of the module
!                      Juelich, 21.1.06 DW
!
!     +++++++++++++++++++++++++++++++++++++++++++++++++
      PRIVATE
      PUBLIC plotdop
      CONTAINS
      SUBROUTINE plotdop(oneD,dimension,&
     &     stars,vacuum,sphhar,atoms,&
     &     input,sym,cell,sliceplot,&
     &     l_noco,cdnfname)
!    *****************************************************
      USE m_outcdn
      USE m_loddop
      USE m_xsf_io
      USE m_cdn_io
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
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
      LOGICAL,INTENT(in)          :: l_noco
      CHARACTER(len=10), INTENT (IN) :: cdnfname

!     .. Local Scalars ..
      REAL          :: s,tec,qint,xdnout,fermiEnergyTemp
      INTEGER       :: i,i1,i2,i3,ii3,ix,iy,iz,jsp,na,nplo
      INTEGER       :: nplot,nq,nt,jm,jspin,iter,ii1,ii2
      LOGICAL       :: twodim,oldform,newform,l_qfix
!     ..
!     .. Local Arrays ..
      COMPLEX :: qpw(stars%ng3,input%jspins),rhtxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
      REAL    :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),rht(vacuum%nmzd,2,input%jspins)
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      COMPLEX :: cdom(1),cdomvz(1,1),cdomvxy(1,1,1)
      INTEGER :: grid(3)
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(atoms%jmtd)
      REAL    :: point(3)
      CHARACTER (len=30) :: filename
      CHARACTER (len=7)  :: append
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename


      oldform = .false.
      INQUIRE(file ="plotin",exist = oldform) 
      IF ( oldform ) THEN 
         CALL priv_old_plot(oneD,dimension,&
     &                stars,vacuum,sphhar,atoms,&
     &                input,sym,cell,sliceplot,&
     &                l_noco,cdnfname)
         RETURN
      ENDIF
      INQUIRE(file ="plot_inp",exist= newform)
      IF (.NOT.newform) THEN !no input file exists, create a template and
                            !exit
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
         CALL juDFT_error("Missing input for plot; modify plot_inp"&
     &        ,calledby ="plotdop")
      ENDIF

      ! new input
      !<-- Open the charge/potential file
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
      !<--Perhaps only the core charge should be plotted
      IF (input%score) THEN
         OPEN (17,file='cdnc',form='unformatted',status='old')
         REWIND 17
         DO jspin = 1,input%jspins
            DO nt = 1,atoms%ntype
               jm = atoms%jri(nt)
               READ (17) (rhocc(i),i=1,jm)
               DO i = 1,atoms%jri(nt)
                  rho(i,0,nt,jspin) = rho(i,0,nt,jspin) - rhocc(i)/2.0&
     &                 /SQRT( pi_const )
               ENDDO
               READ (17) tec
            ENDDO
            READ (17) qint
            qpw(1,jspin) = qpw(1,jspin) - qint/cell%volint
         ENDDO
         CLOSE (17)
      END IF
      !>
      !>
      !<-- Open the plot_inp file for input
      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)')    nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (xsf) THEN
         IF (l_noco) THEN
             append = '_pl.xsf'
             OPEN (55,file = trim(cdnfname)//append,form='formatted')
         ELSE
             OPEN(55,file="plot.xsf")
         ENDIF
         CALL xsf_WRITE_atoms(&
     &                        55,atoms,input%film,oneD%odi%d1,cell%amat)
      ENDIF
      !<-- Loop over all plots
      DO nplo=1,nplot
         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) &
     &        CALL juDFT_error("Illegal grid size in plot",calledby&
     &        ="plotdop")
         IF (.NOT.twodim.AND.ANY(grid<1)) &
     &        CALL juDFT_error("Illegal grid size in plot",calledby&
     &        ="plotdop")
         IF (twodim) grid(3) = 1
         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(cell%amat,vec1)
            vec2=matmul(cell%amat,vec2)
            vec3=matmul(cell%amat,vec3)
            zero=matmul(cell%amat,zero)
         ENDIF
         !Open the file
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo
         IF (xsf) THEN
            CALL xsf_WRITE_header(55,twodim,filename,vec1,vec2,vec3,zero&
     &           ,grid)
         ELSE
            IF (l_noco) THEN
               OPEN (55,file = filename//cdnfname,form='formatted')
            ELSE
               OPEN (55,file = filename,form='formatted')
            ENDIF
         ENDIF
         !loop over spins
         DO jsp = 1,input%jspins
            !loop over all points
            DO iz = 0,grid(3)-1
               DO iy = 0,grid(2)-1
                  xloop:DO ix = 0,grid(1)-1
                    point = zero + vec1*REAL(ix)/(grid(1)-1) +&
     &                             vec2*REAL(iy)/(grid(2)-1)
                    IF (.NOT.twodim) point = point +&
     &                             vec3*REAL(iz)/(grid(3)-1)
                    !Check if the point is in MT-sphere
                    ii1 = 3
                    ii2 = 3
                    ii3 = 3
                    IF (input%film .AND. .NOT.oneD%odi%d1) ii3 = 0
                    IF (oneD%odi%d1) THEN
                       ii1 = 0 ; ii2 = 0
                    END IF
                    DO  i1 = -ii1,ii1
                       DO  i2 = -ii2,ii2
                          DO  i3 = -ii3,ii3
                             pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
                             na = 0
                             DO nt = 1,atoms%ntype
                                DO nq = 1,atoms%neq(nt)
                                   na   = na + 1
                                   s  = SQRT(dot_PRODUCT(atoms%pos(:,na)&
     &                                  -pt,atoms%pos(:,na)-pt))
                                   IF (s<atoms%rmsh(atoms%jri(nt),nt)) THEN
                                      CALL outcdn(&
     &                                     pt,nt,na,0,1,jsp,sliceplot,stars,&
     &                                     vacuum,sphhar,atoms,sym,cell,oneD, &
     &                                     qpw,rhtxy,rho,rht,xdnout)
                                      IF (xsf) THEN
                                         write(55,*) xdnout
                                      ELSE
                                         IF (twodim) THEN
                                            WRITE(55,'(3e15.7)')&
     &                                           point(1:2),xdnout
                                         ELSE
                                            WRITE(55,'(4e15.7)') point&
     &                                           ,xdnout
                                         ENDIF
                                      ENDIF
                                      CYCLE xloop
                                   ENDIF
                                ENDDO
                             ENDDO !nt
                          ENDDO
                       ENDDO
                    ENDDO !i1
                    !Check for point in vacuum

                    IF (input%film.AND..NOT.oneD%odi%d1.AND.ABS(point(3))>=cell%z1) THEN
                       CALL outcdn(&
     &                      point,0,0,1,0,jsp,sliceplot,stars,&
     &                      vacuum,sphhar,atoms,sym,cell,oneD,&
     &                      qpw,rhtxy,rho,rht,xdnout)
                       IF (xsf) THEN
                          write(55,*) xdnout
                       ELSE
                          IF (twodim) THEN
                             WRITE(55,'(3e15.7)') point(1:2),xdnout
                          ELSE
                             WRITE(55,'(4e15.7)') point(:),xdnout
                          ENDIF
                       ENDIF
                       CYCLE xloop
                    END IF
                    IF (oneD%odi%d1) THEN
                       IF (SQRT((pt(1))**2 + (pt(2))**2)>=cell%z1) THEN
                          CALL outcdn(&
     &                         pt,0,0,1,0,jsp,sliceplot,stars,&
     &                         vacuum,sphhar,atoms,sym,cell,oneD, &
     &                         qpw,rhtxy,rho,rht,xdnout)
                          IF (xsf) THEN
                             WRITE(55,*) xdnout
                          ELSE
                             IF (twodim) THEN
                                WRITE (55,'(3e15.7)') point(1:2),xdnout
                             ELSE
                                WRITE (55,'(4e15.7)') point(:),xdnout
                             ENDIF
                          ENDIF
                          CYCLE xloop
                       END IF
                    END IF
                    CALL outcdn(&
     &                   point,0,0,0,2,jsp,sliceplot,stars,&
     &                   vacuum,sphhar,atoms,sym,cell,oneD,&
     &                   qpw,rhtxy,rho,rht,xdnout)
                    IF (xsf) THEN
                       WRITE(55,*) xdnout
                    ELSE
                       IF (twodim) THEN
                          WRITE(55,'(3e15.7)') point(1:2),xdnout
                       ELSE
                          WRITE(55,'(4e15.7)') point(:),xdnout
                       ENDIF
                    ENDIF
                 ENDDO xloop
              ENDDO
           ENDDO !z-loop
           IF (xsf.AND.jsp /= input%jspins) CALL xsf_WRITE_newblock(55,twodim&
     &          ,vec1,vec2,vec3,zero,grid)
        ENDDO !Spin-loop
        IF (xsf) THEN
           CALL xsf_WRITE_endblock(55,twodim)
        ELSE
           CLOSE(55)
        ENDIF
      ENDDO   !nplot      
      CLOSE(18)
      IF (xsf) CLOSE(55)
      RETURN
      END SUBROUTINE plotdop
!------------------------------------------
!     The old subroutine from Fleur is here
!------------------------------------------
      SUBROUTINE priv_old_plot(oneD,dimension,&
     &                stars,vacuum,sphhar,atoms,&
     &                input,sym,cell,sliceplot,&
     &                l_noco,cdnfname)
!
      USE m_outcdn
      USE m_loddop
      USE m_cdn_io
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
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
      LOGICAL,INTENT(IN)          :: l_noco
      CHARACTER(len=10), INTENT (IN) :: cdnfname
!     ..
!     .. Array Arguments ..
!-odim
!+odim
!     ..
!     .. Local Scalars ..
      REAL rx,ry,s,sl,sm,su,x,xm,y,ym,sfp,xdnout
      INTEGER i,i1,i2,i3,ii3,imshx,imshy,ix,iy,j,jsp,na,nfile,nplo,&
     &        nplot,nq,nt,nplott,jm,jspin,iter,ii1,ii2
      LOGICAL twodim
!     ..
!     .. Local Arrays ..
      COMPLEX qpw(stars%ng3,input%jspins),rhtxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
      REAL rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),rht(vacuum%nmzd,2,input%jspins)
      REAL ptp(3),rngl(3),rngm(3),rngu(3),tl(3),tm(3)
      REAL tu(3),vx1(3),vx2(3),tl_r(3),tm_r(3),tu_r(3),rhocc(atoms%jmtd,atoms%ntype,dimension%jspd)
      REAL tec(atoms%ntype,dimension%jspd),qint(atoms%ntype,dimension%jspd)
      REAL pt(3),a3(3)
      REAL, ALLOCATABLE :: cdn(:,:)
      CHARACTER(len=10), ALLOCATABLE :: plotname(:)
!     ..
      sfp = 2.0 * sqrt( pi_const )
      a3(3) = cell%amat(3,3)
!     ..
      OPEN (20,file=cdnfname,form='unformatted',status='old')
      CALL loddop(&
     &            stars,vacuum,atoms,sphhar,&
     &            input,sym,&
     &            20,&
     &            iter,rho,qpw,rht,rhtxy)
!
      IF (input%score) THEN
         CALL readCoreDensity(input,atoms,dimension,rhocc,tec,qint)
         DO jspin = 1,input%jspins
            DO nt = 1,atoms%ntype
               jm = atoms%jri(nt)
               DO i = 1,atoms%jri(nt)
                  rho(i,0,nt,jspin) = rho(i,0,nt,jspin) - rhocc(i,nt,jspin)/sfp
               ENDDO
             ENDDO
            qpw(1,jspin) = qpw(1,jspin) - qint(nt,jspin)/cell%volint
         ENDDO
      END IF
      OPEN (18,file='plotin')
      READ (18,FMT='(7x,l1)') twodim

      READ (18,FMT=8000) nplot
      ALLOCATE ( plotname(nplot) )
 8000 FORMAT (6x,i2)
      nplott=nplot
      if (nplot.eq.1) nplott=2
      DO 140 nplo = 1,nplot
         IF (twodim.OR.(nplo.eq.1)) THEN
           nfile = 55 + nplo
           READ (18,FMT='(a,3x,a)') plotname(nplo)
           IF (l_noco) THEN
              OPEN (nfile,file=plotname(nplo)//cdnfname,&
     &                    form='formatted')
           ELSE
              OPEN (nfile,file=plotname(nplo),form='formatted')
           ENDIF
           READ (18,FMT=8010) (tu_r(i),i=1,3)
           READ (18,FMT=8010) (tm_r(i),i=1,3)
           READ (18,FMT=8010) (tl_r(i),i=1,3)
 8010      FORMAT (4f10.6)
           READ (18,FMT=8020) imshx,imshy
 8020      FORMAT (6x,i5,7x,i5)
           ALLOCATE (cdn(imshx,imshy))
         ENDIF
         IF (twodim) THEN
           DO i=1,3
             tu(i) = tu_r(i)
             tm(i) = tm_r(i)
             tl(i) = tl_r(i)
           ENDDO 
         ELSE
           DO i=1,2
             tu(i) = tu_r(i)
             tm(i) = tm_r(i)
             tl(i) = tl_r(i)
           ENDDO 
           tu(3) = tu_r(3) + a3(3)*(nplo-1)/(2*nplott-2)
           tm(3) = tm_r(3) + a3(3)*(nplo-1)/(2*nplott-2)
           tl(3) = tl_r(3) + a3(3)*(nplo-1)/(2*nplott-2)
         ENDIF

         tu(3) = tu(3)/a3(3)
         tm(3) = tm(3)/a3(3)
         tl(3) = tl(3)/a3(3)
!--->    evaluate cartesian coordinates of positions
         DO 20 i = 1,3
            su = 0.
            sm = 0.
            sl = 0.
            DO 10 j = 1,3
               su = su + cell%amat(i,j)*tu(j)
               sm = sm + cell%amat(i,j)*tm(j)
               sl = sl + cell%amat(i,j)*tl(j)
   10       CONTINUE
            rngu(i) = su
            rngm(i) = sm
            rngl(i) = sl
   20    CONTINUE
         DO 30 i = 1,3
            vx1(i) = rngu(i) - rngm(i)
            vx2(i) = rngl(i) - rngm(i)
   30    CONTINUE
         rx = sqrt(vx1(1)*vx1(1)+vx1(2)*vx1(2)+vx1(3)*vx1(3))
         ry = sqrt(vx2(1)*vx2(1)+vx2(2)*vx2(2)+vx2(3)*vx2(3))
         DO 130 jsp = 1,input%jspins
            WRITE (16,FMT=8030) rx,ry
 8030       FORMAT (2f10.6)
            WRITE (nfile,FMT=8050) imshy,imshx,ry,rx
            WRITE (16,FMT=8050) imshx,imshy
            xm = imshx - 1
            ym = imshy - 1
            DO 120 ix = 1,imshx
               DO 110 iy = 1,imshy
                  x = ix - 1
                  y = iy - 1
                  pt(1) = rngm(1) + vx1(1)*x/xm + vx2(1)*y/ym
                  pt(2) = rngm(2) + vx1(2)*x/xm + vx2(2)*y/ym
                  pt(3) = rngm(3) + vx1(3)*x/xm + vx2(3)*y/ym
                  ii1 = 3
                  ii2 = 3
                  ii3 = 3
                  IF (input%film .AND. .NOT.oneD%odi%d1) ii3 = 0
                  IF (oneD%odi%d1) THEN
                     ii1 = 0 ; ii2 = 0
                  END IF
                  DO 100 i1 = -ii1,ii1
                     DO 90 i2 = -ii2,ii2
                        DO 80 i3 = -ii3,ii3
                           DO 40 i = 1,3
                              ptp(i) = pt(i) + i1*cell%amat(i,1) +&
     &                                 i2*cell%amat(i,2) + i3*cell%amat(i,3)
   40                      CONTINUE
                           na = 0
                           DO 70 nt = 1,atoms%ntype
                              DO 60 nq = 1,atoms%neq(nt)
                                 na = na + 1
                                 s = 0.
                                 DO 50 i = 1,3
                                    s = s + (ptp(i)-atoms%pos(i,na))**2
   50                            CONTINUE
                                 s = sqrt(s)
                                 IF (s.LT.atoms%rmsh(atoms%jri(nt),nt)) THEN
                                    CALL outcdn(&
     &                  ptp,nt,na,0,1,jsp,sliceplot,stars,&
     &                  vacuum,sphhar,atoms,sym,cell,oneD,&
     &                  qpw,rhtxy,rho,rht,xdnout)
                                    cdn(ix,iy) = xdnout
                                    GO TO 110
                                 END IF
   60                         CONTINUE
   70                      CONTINUE
   80                   CONTINUE
   90                CONTINUE
  100             CONTINUE
                  IF (input%film .AND. .NOT.oneD%odi%d1) THEN
                    IF (abs(pt(3)).GE.cell%z1) THEN
                      CALL outcdn(&
     &                  pt,0,0,1,0,jsp,sliceplot,stars,&
     &                  vacuum,sphhar,atoms,sym,cell,oneD,&
     &                  qpw,rhtxy,rho,rht,xdnout)
                        cdn(ix,iy) = xdnout
                      GO TO 110
                    END IF
                  END IF
!-odim
                  IF (oneD%odi%d1) THEN
                     IF (sqrt((pt(1))**2 + (pt(2))**2).GE.cell%z1) THEN
                      CALL outcdn(&
     &                  pt,0,0,1,0,jsp,sliceplot,stars,&
     &                  vacuum,sphhar,atoms,sym,cell,oneD,&
     &                  qpw,rhtxy,rho,rht,xdnout)
                        cdn(ix,iy) = xdnout
                      GO TO 110
                     END IF
                  END IF
!+odim
                  CALL outcdn(&
     &                  pt,0,0,0,2,jsp,sliceplot,stars,&
     &                  vacuum,sphhar,atoms,sym,cell,oneD,&
     &                  qpw,rhtxy,rho,rht,xdnout)
                  cdn(ix,iy) = xdnout
  110          CONTINUE
  120       CONTINUE
            WRITE (nfile,FMT=8040) ((cdn(ix,iy),ix=1,imshx),iy=1,imshy)
 8040       FORMAT (5e14.6)
 8050       FORMAT (2i5,2f12.6)
  130    CONTINUE
         IF (twodim) THEN
            DEALLOCATE (cdn)
            CLOSE (nfile)
         ENDIF
  140 CONTINUE
      DEALLOCATE ( plotname )
      IF (.not.twodim) CLOSE (nfile)
      CLOSE (18)
      CLOSE (20)
      END SUBROUTINE priv_old_plot

      END MODULE m_plotdop
