!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_plot
      use m_juDFT
!     +++++++++++++++++++++++++++++++++++++++++++++++++
!     plots the periodic part of the wavefunction at the given k-point
!       and specified bands. Needed for real-space
!       visualization of the Wannier Functions
!     
!     if twodim = .false. a 3-D plot with nplot layers in z-direction
!     is constructed; the 3x3 matrix gives the 3 vectors of the cell ..
!     .gustav
!
!    Changed the input/output for easier use. 
!    This subroutine uses the file plot_inp for input. 
!       based on the routine by:  Juelich, 21.1.06 D.Wortmann
!            Y.Mokrousov 16.8.6 Hamburg
!
!    First of all, if xsf then the files are written in the following
!       format: k.bnd.real.jsp.xsf and k.bnd.imag.jsp.xsf
!               for reading with XCrysDen,where 
!               k is the kpoint, bnd is the band index
!               imag/real stand for the real and imaginary parts of wavefn
!               jsp stands for spin, for all the k-points and bands
!    If the xsf is false, than the files are written in the
!       format suitable for the Vanderbilt's code:
!       UNK0000k.jsp, having all the bands in this file, real and imag
!
!    Then, if xsf=true, then:
! 
!    If slice=T, the kk specifies the k-point and the nnne
!       specifies the number of the band, the wavefunction is
!       written to the file kk.nnne.real.jsp.xsf and kk.nnne.imag.jsp.xsf
!     +++++++++++++++++++++++++++++++++++++++++++++++++
      CONTAINS
      SUBROUTINE wann_plot(
     >      vacuum,stars,cell,atoms,sym,
     >     nv2d,jspin,ntypsd,
     >     film,slice,ikpt,nnne,kk,
     >     nvd,nv,lmd,bkpt,k1,k2,k3,evac,vz,
     >     nslibd,nbasfcn,neigd,ff,gg,flo,acof,bcof,ccof,zMat,
     >     l_noco,l_ss,qss,addnoco,
     >     index_kq,l_sgwf)
!    *****************************************************
      USE m_constants
      USE m_types
      USE m_wann_real
      USE m_xsf_io
      USE m_wann_plot_vac
      USE m_wann_2dvacabcof
   
      IMPLICIT NONE

      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_stars),INTENT(IN)       :: stars
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_mat),INTENT(IN)         :: zMat

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypsd,ikpt,jspin,nv2d
      INTEGER, INTENT (IN) :: nvd,nv,lmd
      INTEGER, INTENT (IN) :: nslibd,nbasfcn,neigd
      INTEGER, INTENT (IN) :: nnne,kk,addnoco,index_kq
      LOGICAL, INTENT (IN) :: slice,film
      LOGICAL, INTENT (IN) :: l_noco,l_ss,l_sgwf
      REAL,    INTENT (IN) :: qss(3)
      real,    intent (in) :: vz(:,:)
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(:),k2(:),k3(:)
      REAL,    INTENT (IN) :: bkpt(3),evac(2)
      REAL,    INTENT (IN) :: ff(:,:,:,0:)
      REAL,    INTENT (IN) :: gg(:,:,:,0:)
      REAL,    INTENT (IN) :: flo(:,:,:,:)
      COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:,:)
      COMPLEX, INTENT (IN) :: acof(:,0:,:)
      COMPLEX, INTENT (IN) :: bcof(:,0:,:)
    
!     ..
!     .. Local Scalars ..
      integer n2,k,j,l
      REAL          :: s
      COMPLEX       :: xdnout
      INTEGER       :: i,i1,i2,i3,ii3,ix,iy,iz,na,nplo,nbn
      INTEGER       :: nbmin,nbmax,n
      INTEGER       :: nplot,nq,nt,jm,iter,ii1,ii2
      CHARACTER*8   :: dop,iop
      LOGICAL       :: twodim
!     ..
!     .. Local Arrays ..
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      INTEGER :: grid(3)
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(atoms%jmtd)
      REAL    :: point(3),post(3,atoms%nat)
      CHARACTER(len=30):: filename
      CHARACTER(len=20):: name1,name2,name3
      CHARACTER(len=10):: vandername
      CHARACTER*8      :: name(10)
c     ..basis wavefunctions in the vacuum
      complex, allocatable :: ac(:,:,:),bc(:,:,:)
      complex, allocatable :: ac_1(:,:,:),bc_1(:,:,:)
      real,    allocatable :: u(:,:,:),ue(:,:,:)
      real,    allocatable :: u_1(:,:,:),ue_1(:,:,:)

      integer ik,nv2,ivac,jvac
      complex factor
      real fas


      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename
      
      intrinsic real,aimag,conjg
c      external dotirp
c      real dotirp  !module now

      call timestart("wann_plot")
      IF (slice) THEN
        nbmin = nnne
        nbmax = nnne
      ELSE
        nbmin = 1
        nbmax = nslibd
      ENDIF

      !write(*,*)'in wann_plot'
      !write(*,*)'jspin:',jspin
      !write(*,*)'addnoco:',addnoco
      !write(*,*)'l_noco,l_ss: ',l_noco,l_ss
      !write(*,*)'qss:',qss

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
         WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.0 "
         WRITE(20,*) "  filename ='plot2' /"
         CLOSE(20)
         WRITE(*,*) "No plot_inp file found. Created a template"
         CALL juDFT_error("Missing input for plot; modify plot_inp"
     +        ,calledby ="wann_plot")
      ENDIF

c make preparations for plotting in vacuum

      if(film)then
         allocate ( ac(nv2d,nslibd,2),bc(nv2d,nslibd,2),
     +         u(vacuum%nmzd,nv2d,vacuum%nvac),
     +         ue(vacuum%nmzd,nv2d,vacuum%nvac))
         call wann_2dvacabcof(
     >         cell,vacuum,
     >         nv2d,nslibd,vz,nv,bkpt,
     >         nvd,k1,k2,k3,evac,
     >         nbasfcn,neigd,zMat,
     <         ac,bc,u,ue,addnoco,l_ss,qss,jspin)
      endif

     


      !<-- Open the plot_inp file for input
      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)') nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (nplot.ge.2) 
     &     CALL juDFT_error
     +     ("plots one by one, please, this is not charge density"
     +     ,calledby="wann_plot")
      !<-- Loop over all plots
      DO nplo=1,nplot
         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) 
     +        CALL juDFT_error("Illegal grid size in plot",calledby
     +        ="wann_plot")
         IF (.NOT.twodim.AND.ANY(grid<1)) 
     +        CALL juDFT_error("Illegal grid size in plot",calledby
     +        ="wann_plot")
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
c..loop by the bands
         bands:DO nbn = nbmin,nbmax

         IF (xsf) THEN
            write (name1,22) ikpt,nbn,jspin
   22       format (i5.5,'.',i3.3,'.real.',i1,'.xsf')
            write (name2,23) ikpt,nbn,jspin
   23       format (i5.5,'.',i3.3,'.imag.',i1,'.xsf')
            write (name3,24) ikpt,nbn,jspin
   24       format (i5.5,'.',i3.3,'.absv.',i1,'.xsf')
            OPEN(55,file=name1)
            CALL xsf_WRITE_atoms(55,atoms,film,cell%amat)
            OPEN(56,file=name2)
            CALL xsf_WRITE_atoms(56,atoms,film,cell%amat)
            OPEN(57,file=name3)
            CALL xsf_WRITE_atoms(57,atoms,film,cell%amat)
            CALL xsf_WRITE_header(55,twodim,filename,(vec1),
     &       (vec2),(vec3),zero
     $           ,grid)
            CALL xsf_WRITE_header(56,twodim,filename,(vec1),
     &       (vec2),(vec3),zero
     $           ,grid)
            CALL xsf_WRITE_header(57,twodim,filename,(vec1),
     &       (vec2),(vec3),zero
     $           ,grid)
         ELSE
            IF (nbn.EQ.nbmin) THEN
               WRITE (vandername,201) ikpt,jspin
               IF(l_noco) WRITE(vandername,202)ikpt,jspin
               IF(l_sgwf) WRITE(vandername,202)index_kq,jspin
  201          FORMAT ('UNK',i5.5,'.',i1)
  202          FORMAT ('RNK',i5.5,'.',i1)
               OPEN(55,file=vandername)
               IF(.NOT.l_sgwf) THEN
                  WRITE(55,7) grid(1),grid(2),grid(3),ikpt,nslibd
               ELSE
                  WRITE(55,7) grid(1),grid(2),grid(3),index_kq,nslibd
               ENDIF
    7          FORMAT (5i4)
            ENDIF
         ENDIF

         if(film)then
           fas=-bkpt(3)*cell%bmat(3,3)*cell%z1
           factor=cmplx(cos(fas),sin(fas))
         else
           factor=cmplx(1.0,0.0)
         endif   

         DO iz = 0,grid(3)-1
          DO iy = 0,grid(2)-1
           xloop:DO ix = 0,grid(1)-1
            point = zero+vec1*REAL(ix)/grid(1)+vec2*REAL(iy)
     $                 /grid(2)
            IF (.NOT.twodim) point = point+vec3*REAL(iz)/grid(3)
!Check if the point is in MT-sphere

             ii1 = 3
             ii2 = 3
             ii3 = 3
             IF (film ) ii3 = 0
           
             DO  i1 = -ii1,ii1
              DO  i2 = -ii2,ii2
               DO  i3 = -ii3,ii3
                pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
                na = 0
                DO nt = 1,atoms%ntype
                 DO nq = 1,atoms%neq(nt)
                  na   = na + 1
                  s  = SQRT(dot_PRODUCT(
     +                 atoms%pos(:,na)-pt,
     +                 atoms%pos(:,na)-pt))
                  IF (s<atoms%rmsh(atoms%jri(nt),nt)) THEN
                    CALL wann_real(
     >            pt,nt,na,0,1,bkpt,nbn,
     >            atoms,cell,sym,
     >            ff,gg,flo,
     >            acof(nbn,:,:),
     >            bcof(nbn,:,:),
     >            ccof(:,nbn,:,:),zMat,
     >            nv,k1,k2,k3,lmd,
     >            nbasfcn,l_ss,qss,
     >            jspin,addnoco,
     <            xdnout)
                    xdnout=xdnout*factor
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
             IF (film.AND.ABS(point(3))>=cell%z1) THEN
                ivac=1
                if (point(3).lt. 0.0)ivac=2
                jvac=ivac
                if(vacuum%nvac==1)jvac=1
                call wann_plot_vac(point,
     >            cell,vacuum,
     >            nv2d,evac,bkpt,vz,
     >            jspin,k1,k2,k3,nvd,
     >            nbasfcn,neigd,nv,
     >            nslibd,ac(:,nbn,ivac),
     &            bc(:,nbn,ivac),
     &            u(:,:,jvac),ue(:,:,jvac),xdnout)
                xdnout=xdnout*factor
                if(real(xdnout).gt.9.0 .or.real(xdnout).lt.-9.0
     &        .or.aimag(xdnout).gt.9.0 .or. aimag(xdnout).lt.-9.0)then
                xdnout=cmplx(0.0,0.0)
                print*,"vac-problem at z=",point(3)
                endif
c               CALL wann_real(
c     >              point,0,0,1,0,bkpt,
c     >              n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
c     >              natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
c     >              nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
c     >              lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
c     >              omtil,amat,bmat,nlod,llod,nlo,llo,
c     >              ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
c     >              ccof(:,nbn,:,:),z(:,nbn),
c     >              nv,k1,k2,k3,lmd,nbasfcn,
c     <              xdnout)
              IF (xsf) THEN
                 WRITE(55,*) real(xdnout)
                 WRITE(56,*) aimag(xdnout)
                 WRITE(57,*) real(xdnout*conjg(xdnout))
              ELSE
                 WRITE(55,8) real(xdnout),aimag(xdnout)
              ENDIF
              CYCLE xloop
             END IF
            
            
             CALL wann_real(
     >            point,0,0,0,2,bkpt,nbn,
     >            atoms,cell,sym,
     >            ff,gg,flo,
     >            acof(nbn,:,:),
     >            bcof(nbn,:,:),
     >            ccof(:,nbn,:,:),zMat,
     >            nv,k1,k2,k3,lmd,
     >            nbasfcn,l_ss,qss,
     >            jspin,addnoco,
     <            xdnout)
             xdnout=xdnout*factor
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
    8 FORMAT (f20.12,1x,f20.12)


      call timestop("wann_plot")
      RETURN
      END SUBROUTINE wann_plot
!------------------------------------------
      
      END MODULE m_wann_plot
