!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_plot_from_lapw
      use m_juDFT

c****************************************************************
c     read in WF1.lapw (WF2.lapw) and plot wannierfunctions
c     Frank Freimuth, October2006
c****************************************************************
      CONTAINS

      SUBROUTINE wann_plot_from_lapw(nv2d,jspins,odi,ods,n3d,nmzxyd,n2d,
     >     ntypsd,
     >     ntype,lmaxd,jmtd,natd,nmzd,neq,nq3,nvac,
     >     nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
     >     invs,invs2,z1,delz,ngopr,ntypsy,jri,pos,zatom,
     >     lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,nnne,kk,
     >     nlod,llod,lmd,omtil,nlo,llo)

      USE m_types, ONLY : od_inp, od_sym
      USE m_xsf_io
      USE m_wann_lapw_int_plot
      USE m_wann_lapw_sph_plot
 
      implicit none
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n3d,nmzxyd,n2d,ntypsd,nv2d
      INTEGER, INTENT (IN) :: lmaxd,jmtd,natd,nmzd
      INTEGER, INTENT (IN) :: nq3,nvac,nmz,nmzxy,nq2,nop,nop2,ntype
      INTEGER, INTENT (IN) :: lmd,llod,nlod
      INTEGER, INTENT (IN) :: nnne,kk,jspins
      LOGICAL, INTENT (IN) :: symor,invs,slice,invs2,film
      REAL,    INTENT (IN) :: z1,delz,volint,omtil
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ngopr(natd),ntypsy(natd),lmax(ntype)
      INTEGER, INTENT (IN) :: jri(ntype),neq(ntype),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invtab(nop),nlo(ntype),llo(nlod,ntype)
      REAL,    INTENT (IN) :: zatom(:),amat(3,3),bmat(3,3),pos(3,natd)
      REAL,    INTENT (IN) :: rmsh(jmtd,ntype),tau(3,nop)
      REAL,    INTENT (IN) :: bbmat(3,3)
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
      integer nplo,nplot,nbn,nbmin,nbmax,ix,iy,iz,ii1,ii2,ii3,i1,i2,i3
      integer nt,na,nq,ivac,jvac
      real amat_old(3,3)
      real pos_old(3,natd)
      real s
      logical l_file
      integer jspin
      integer num_wann,cell1,cell2,cell3
      character(len=3) :: spin12(2)
      data spin12/'WF1' , 'WF2'/
      complex,allocatable::wann_acof(:,:,:,:,:,:)
      complex,allocatable::wann_bcof(:,:,:,:,:,:)
      complex,allocatable::wann_ccof(:,:,:,:,:,:,:)
      real,allocatable::ff(:,:,:,:)
      real,allocatable::gg(:,:,:,:)
      real,allocatable::flo(:,:,:,:)
      complex xdnout
      LOGICAL       :: twodim
      LOGICAL :: cartesian,xsf
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3),v(3),point(3)
      INTEGER :: grid(3)
      CHARACTER(len=30):: filename
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename
      CHARACTER(len=20):: name1,name2,name3
      CHARACTER(len=10):: vandername
      CHARACTER*8      :: name(10)
      integer ntype2,jmtd2,lmaxd2,nlod2,natd2,lmd2,llod2
      integer unigrid(4)
      real tmp_zatom(1:ntype),tmp_dx(1:ntype),tmp_rmt(1:ntype)
      integer tmp_jri(1:ntype),tmp_neq(1:ntype),tmp_rmsh(jmtd,ntype)
      complex,allocatable::wannint(:,:,:,:)
      real delta1,delta2,plot_time_int,plot_read_data_time
      real plot_time_sph,time_sqrt
      plot_time_sph=0.0
      plot_time_int=0.0
      plot_read_data_time=0.0
      time_sqrt=0.0


      INQUIRE(file ="plot_inp",exist= twodim)
      IF (.NOT.twodim) THEN !no input file exists, create a template and
                            !exit
         OPEN(20,file ="plot_inp")
         WRITE(20,'(i2,a5,l1)') 1,",xsf=",.false.
c         WRITE(20,*) "&PLOT twodim=t,cartesian=t"
c         WRITE(20,*) "  vec1(1)=10.0 vec2(2)=10.0"
c         WRITE(20,*) "  filename='plot1' /"
         WRITE(20,*) "&PLOT twodim=f,cartesian=f"
         WRITE(20,*) "  vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0 "
         WRITE(20,*) "  vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0 "
         WRITE(20,*) "  vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0 "
         WRITE(20,*) "  grid(1)=30  grid(2)=30  grid(3)=30  "
c         WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.5 "
         WRITE(20,*) "  filename ='plot2' /"
         CLOSE(20)
         WRITE(*,*) "No plot_inp file found. Created a template"
         CALL juDFT_error("Missing input for plot; modify plot_inp"
     +        ,calledby ="wann_plot_from_lapw")
      ENDIF

      l_file=.false.
      do jspin=1,2  !spin loop
         call cpu_time(delta1)
c*******************************************************
c          read in data
c*******************************************************
       inquire(file=spin12(jspin)//'.lapw',exist=l_file)
       IF(.NOT.l_file) CALL juDFT_error("WF..lapw not found",calledby
     +      ="wann_plot_from_lapw")
       print*,"open file WF..lapw"
       open(344,file=spin12(jspin)//'.lapw',form='unformatted')
       read(344)num_wann,cell1,cell2,cell3
       read(344)amat_old(3,3)
       read(344)ntype2,jmtd2,lmaxd2,nlod2,natd2,lmd2,llod2
       IF(ntype2/=ntype) CALL juDFT_error("ntype2",calledby
     +      ="wann_plot_from_lapw")
       IF(jmtd2/=jmtd) CALL juDFT_error("jmtd2",calledby
     +      ="wann_plot_from_lapw")
       IF(lmaxd2/=lmaxd) CALL juDFT_error("lmaxd2",calledby
     +      ="wann_plot_from_lapw")
       IF(nlod2/=nlod) CALL juDFT_error("nlod2",calledby
     +      ="wann_plot_from_lapw")
       IF(natd2/=natd) CALL juDFT_error("natd2",calledby
     +      ="wann_plot_from_lapw")
       IF(lmd2/=lmd) CALL juDFT_error("lmd2",calledby
     +      ="wann_plot_from_lapw")
       IF(nlod2/=nlod) CALL juDFT_error("nlod2",calledby
     +      ="wann_plot_from_lapw")
       read(344)pos_old(1:3,1:natd)
       read(344)tmp_neq(1:ntype)
       read(344)tmp_zatom(1:ntype)
       read(344)tmp_dx(1:ntype)
       read(344)tmp_rmt(1:ntype)
       read(344)tmp_jri(1:ntype)
       read(344)tmp_rmsh(1:jmtd,1:ntype)
c*******************************************************
c         read in radial wavefunctions
c*******************************************************
       allocate(ff(1:ntype,1:jmtd,1:2,0:lmaxd))
       allocate(gg(1:ntype,1:jmtd,1:2,0:lmaxd))
       allocate(flo(1:ntype,1:jmtd,1:2,1:nlod))
       read(344)ff(1:ntype,1:jmtd,1:2,0:lmaxd)
       read(344)gg(1:ntype,1:jmtd,1:2,0:lmaxd)
       read(344)flo(1:ntype,1:jmtd,1:2,1:nlod)
c*******************************************************
c         read in a-,b-, and c-coefficients
c*******************************************************

       allocate(wann_acof(num_wann,0:lmd,natd,-cell1:cell1,-cell2:cell2,
     &   -cell3:cell3))
       allocate(wann_bcof(num_wann,0:lmd,natd,-cell1:cell1,-cell2:cell2,
     &   -cell3:cell3))
       allocate(wann_ccof(num_wann,-llod:llod,nlod,natd,-cell1:cell1,
     &   -cell2:cell2,-cell3:cell3))

       read(344)wann_acof(1:num_wann,0:lmd,1:natd,-cell1:cell1,
     &                    -cell2:cell2,-cell3:cell3)
       read(344)wann_bcof(1:num_wann,0:lmd,1:natd,-cell1:cell1,
     &                    -cell2:cell2,-cell3:cell3)
       read(344)wann_ccof(1:num_wann,-llod:llod,1:nlod,1:natd,
     &                -cell1:cell1,-cell2:cell2,-cell3:cell3)
c*********************************************************
c        read in planewave expansion
c*********************************************************
       read(344)unigrid(1:4)
       allocate(wannint(-unigrid(4):unigrid(4),-unigrid(4):unigrid(4),
     ,           -unigrid(4):unigrid(4),num_wann))
       read(344)wannint(:,:,:,:)

       print*,"num_wann=",num_wann
       print*,"wannier supercell: ",cell1,cell2,cell3
       close(344)
       call cpu_time(delta2)
       plot_read_data_time=delta2-delta1

c**********************************************************************
c       plot the wannierfunction
c**********************************************************************
      !<-- Open the plot_inp file for input
      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)') nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (nplot.ge.2) 
     &     CALL juDFT_error
     +     ("plots one by one, please, this is not charge density"
     +     ,calledby="wann_plot_from_lapw")
      !<-- Loop over all plots
      DO nplo=1,nplot
         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) 
     +        CALL juDFT_error("Illegal grid size in plot",calledby
     +        ="wann_plot_from_lapw")
         IF (.NOT.twodim.AND.ANY(grid<1)) 
     +        CALL juDFT_error("Illegal grid size in plot",calledby
     +        ="wann_plot_from_lapw")
         IF (twodim) grid(3) = 1
         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(amat,vec1)
            vec2=matmul(amat,vec2)
            vec3=matmul(amat,vec3)
            zero=matmul(amat,zero)
         ENDIF
         !Open the file
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo
c..loop by the wannierfunctions
         nbmin=1
         nbmax=num_wann
         bands:DO nbn = nbmin,nbmax

         IF (xsf) THEN
            write (name1,22) nbn,jspin
   22       format (i3.3,'.real.',i1,'.xsf')
            write (name2,23) nbn,jspin
   23       format (i3.3,'.imag.',i1,'.xsf')
            write (name3,24) nbn,jspin
   24       format (i3.3,'.absv.',i1,'.xsf')
            OPEN(55,file=name1)
!#if 1==1
            call judft_error("NOT INPLEMENTED")
c$$$#else
c$$$            CALL xsf_WRITE_atoms(
c$$$     >                        55,film,odi%d1,amat,neq(:ntype),
c$$$     >                        zatom(:ntype),pos)
c$$$            OPEN(56,file=name2)
c$$$            CALL xsf_WRITE_atoms(
c$$$     >                        56,film,odi%d1,amat,neq(:ntype),
c$$$     >                        zatom(:ntype),pos)
c$$$            OPEN(57,file=name3)
c$$$            CALL xsf_WRITE_atoms(
c$$$     >                        57,film,odi%d1,amat,neq(:ntype),
c$$$     >                        zatom(:ntype),pos)
c$$$            CALL xsf_WRITE_header(55,twodim,filename,(vec1),
c$$$     &       (vec2),(vec3),zero
c$$$     $           ,grid)
c$$$            CALL xsf_WRITE_header(56,twodim,filename,(vec1),
c$$$     &       (vec2),(vec3),zero
c$$$     $           ,grid)
c$$$            CALL xsf_WRITE_header(57,twodim,filename,(vec1),
c$$$     &       (vec2),(vec3),zero
c$$$     $           ,grid)
c$$$#endif
         ELSE
               WRITE (vandername,201) nbn,jspin
  201          FORMAT (i5.5,'.',i1)            
               OPEN(55,file=vandername)
               WRITE (55,7) grid(1),grid(2),grid(3)
    7          FORMAT (5i4)
         ENDIF

         DO iz = 0,grid(3)-1
          DO iy = 0,grid(2)-1
           xloop:DO ix = 0,grid(1)-1
            point = zero+vec1*REAL(ix)/(grid(1)-1)+vec2*REAL(iy)
     $                 /(grid(2)-1)
            IF (.NOT.twodim) point = point+vec3*REAL(iz)/(grid(3)-1)
!Check if the point is in MT-sphere

             ii1 = cell1
             ii2 = cell2
             ii3 = cell3
             IF (film .AND. .NOT.odi%d1) ii3 = 0
             IF (odi%d1) THEN
                ii1 = 0 ; ii2 = 0
             END IF
             DO  i1 = -ii1,ii1
              DO  i2 = -ii2,ii2
               DO  i3 = -ii3,ii3
                pt = point+MATMUL(amat,(/i1,i2,i3/))
                na = 0
                DO nt = 1,ntype
                 DO nq = 1,neq(nt)
                  na   = na + 1
                  s  = dot_PRODUCT(pos(:,na)-pt,pos(:,na)-pt)
                  IF (s<(rmsh(jri(nt),nt))**2) THEN
                   pt(:)=pt(:)-pos(:,na)
                   call wann_lapw_sph_plot(ff(nt,:,1,:),gg(nt,:,1,:),
     >                   flo(nt,:,1,:),wann_acof(nbn,:,na,-i1,-i2,-i3),
     >                   wann_bcof(nbn,:,na,-i1,-i2,-i3),
     >                   wann_ccof(nbn,:,:,na,-i1,-i2,-i3),pt,nlo(nt),
     >                   jmtd,lmaxd,nlod,llod,lmd,rmsh(:,nt),lmax(nt),
     >                   llo(:,nt),jri(nt),
     <                   xdnout)
                   IF (xsf) THEN
                      WRITE(55,*) real(xdnout)
                      WRITE(56,*) aimag(xdnout)
                      WRITE(57,*) real(xdnout*conjg(xdnout))
                   ELSE
                      WRITE(55,8) real(xdnout),aimag(xdnout)
                   ENDIF
                   CYCLE xloop
                  ENDIF
                 ENDDO
                ENDDO !nt
               ENDDO
              ENDDO
             ENDDO !i1
!Check for point in vacuum
             IF (film.AND..NOT.odi%d1.AND.ABS(point(3))>=z1) THEN
                ivac=1
                if (point(3).lt. 0.0)ivac=2
                jvac=ivac
                if(nvac==1)jvac=1
                xdnout=0.0
c     call wann_plot_vac
                if(real(xdnout).gt.9.0 .or.real(xdnout).lt.-9.0
     &        .or.aimag(xdnout).gt.9.0 .or. aimag(xdnout).lt.-9.0)then
                xdnout=cmplx(0.0,0.0)
                print*,"vac-problem at z=",point(3)
                endif

              IF (xsf) THEN
                 WRITE(55,*) real(xdnout)
                 WRITE(56,*) aimag(xdnout)
                 WRITE(57,*) real(xdnout*conjg(xdnout))
              ELSE
                 WRITE(55,8) real(xdnout),aimag(xdnout)
              ENDIF
              CYCLE xloop
             END IF
            
             IF (odi%d1) THEN
              IF (SQRT((pt(1))**2 + (pt(2))**2)>=z1) THEN
c     call wann_real
               IF (xsf) THEN
                  WRITE(55,*) real(xdnout)
                  WRITE(56,*) aimag(xdnout)
                  WRITE(57,*) real(xdnout*conjg(xdnout))
               ELSE
                  WRITE (55,8) real(xdnout),aimag(xdnout)
               ENDIF
               CYCLE xloop
              END IF
             END IF
               call wann_lapw_int_plot(point,bmat,unigrid,
     >            wannint(:,:,:,nbn),
     <            xdnout)
c               xdnout=cmplx(0.0,0.0)
             IF (xsf) THEN
                WRITE(55,*) real(xdnout)
                WRITE(56,*) aimag(xdnout)
                WRITE(57,*) real(xdnout*conjg(xdnout))
             ELSE
                WRITE(55,8) real(xdnout),aimag(xdnout)
             ENDIF
            ENDDO xloop
           ENDDO
          ENDDO !z-loop

          IF (xsf) THEN
              CALL xsf_WRITE_endblock(55,twodim)
              CALL xsf_WRITE_endblock(56,twodim)
              CALL xsf_WRITE_endblock(57,twodim)
              CLOSE (55) ; CLOSE (56) ; CLOSE (57)
          ENDIF
c..end of the loop by the bands
          ENDDO bands   
      ENDDO   !nplot      
      CLOSE(18)
      IF (.not.xsf) CLOSE(55)
    8 FORMAT (2f16.12)


       deallocate(wann_acof,wann_bcof,wann_ccof)

       if(jspins.eq.1)exit
      enddo !spinloop
      print*,"plot_time_int=",plot_time_int
      print*,"plot_time_sph=",plot_time_sph
      print*,"plot_read_data_time=",plot_read_data_time
      print*,"time_sqrt=",time_sqrt

      END SUBROUTINE wann_plot_from_lapw
      END MODULE m_wann_plot_from_lapw
