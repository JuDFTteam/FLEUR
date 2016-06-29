!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_plot_um_dat
      use m_juDFT
c******************************************************************
c       plot wannierfunctions directly within fleur
c       based on wann_plot
c       FF, September 2006
c******************************************************************
      CONTAINS
      SUBROUTINE wann_plot_um_dat(
     >          mpi_comm,sortrule,band_min,band_max,l_soc,
     >          l_dulo,l_noco,l_ss,lmaxd,
     >          ntypd,
     >          neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >          nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >          invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >          beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >          jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >          ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >          ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >          z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >          volint,symor,pos,ef,l_bzsym,
     >          l_proj_plot,irecl,wan90version)

      use m_wann_rw_eig
      use m_abcof
      use m_wann_2dvacabcof
      use m_dotir
      use m_radfun
      use m_radflo
      use m_cdnread, only : cdn_read0
      use m_od_types, only : od_inp, od_sym
      use m_loddop
      use m_constants, only : pimach
      use m_wann_real
      use m_xsf_io
      use m_cotra,only:cotra1
      use m_wann_read_umatrix
      use m_wann_kptsrotate
      use m_wann_plot_vac
      use m_wann_abinv


      implicit none
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer mpiierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)
#endif

      integer,intent(in)::band_min(2),band_max(2),mpi_comm
      logical,intent(in) :: l_soc
      logical, intent (in) :: invs,invs2,film,slice,symor
      integer, intent (in) :: lmaxd,ntypd,neigd,nkptd,kk,nnne
      integer, intent (in) :: natd,nop,nvd,jspd,nbasfcn,nq2,nop2
      integer, intent (in) :: llod,nlod,ntype,nwdd,n3d,n2d
      integer, intent (in) :: nmzxyd,nmzd,jmtd,nlhd,nq3,nvac 
      integer, intent (in) :: ntypsd,jspins,k1d,k2d,k3d
      real,    intent (in) :: omtil,e1s,e2s,delz,area,z1,volint
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      real,    intent (in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: nlh(ntypsd),jri(ntypd),ntypsy(natd)
      integer, intent (in) :: nlo(ntypd),llo(nlod,ntypd),lapw_l(ntypd)
      integer, intent (in) :: invtab(nop),mrot(3,3,nop),ngopr(natd)
      integer, intent (in) :: neq(ntypd),lmax(ntypd)
      integer, intent (in) :: wan90version
      integer, intent (in) :: invsat(natd),invsatnr(natd),nkpt(nwdd)
      integer, intent (in) :: irank,isize,nv2d,nmzxy,nmz
      integer, intent (in) :: ulo_der(nlod,ntypd),ig2(n3d)
      real,    intent (in) :: taual(3,natd),rmt(ntypd),dx(ntypd) 
      real,    intent (in) :: amat(3,3),bmat(3,3),bbmat(3,3)
      real,    intent (in) :: rmsh(jmtd,ntypd),tau(3,nop),zatom(ntype)
      real, intent (inout) :: alph(ntypd),beta(ntypd),qss(3)  
      real,    intent (in) :: pos(3,natd),ef
      complex, intent (in) :: ustep(n3d)
      logical, intent (in) :: l_dulo(nlod,ntypd),l_noco,l_ss,l_bzsym
      logical,intent(in)::l_proj_plot
      integer,intent(in)::sortrule
#ifndef CPP_MPI
      integer(8), intent (in) :: irecl
#else
      INTEGER(KIND=MPI_OFFSET_KIND),intent(in) :: irecl
#endif
c-odim
      real,    intent (in) :: sk2(n2d),phi2(n2d)
      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods
c+odim
      logical l_spreadcal
      complex, allocatable::spreads(:,:)
      real,allocatable:: centers(:,:)
cccccccccccccccccc   local variables   cccccccccccccccccccc
      integer lmd,nlotot,n,nmat,nw,ispin,iter,ikpt,ilo
      integer noccbd,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
      integer jsp_start,jsp_end,nrec,nrec1,nbands
      integer nodeu,noded,n_size,na,n_rank,nbnd,nkbnd,idum,jdum,kdum
      integer i1,i2,i3,in,ikpt_k,lda
      integer n_bands(0:neigd),nslibd
      character*8 dop,iop,name(10)
      real bkpt(3),sfp,tpi,wronk,wk,wk_b,phase
      real eig(neigd),cp_time(9)
      logical l_p0,l_bkpts,l_proj,l_amn,l_mmn,l_eig,wann,wann_plott
!!! energy window boundaries
      integer, allocatable :: kveclo(:),nv(:)
      integer, allocatable :: k1(:,:),k2(:,:),k3(:,:)
      real,    allocatable :: we(:)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex, allocatable :: z(:,:),zz(:,:)
#else
      real,    allocatable :: z(:,:),zz(:,:)
#endif
      real,    allocatable :: eigg(:)
      real kpoints(nkptd)
!!! a and b coeff. constructed for each k-point
      complex, allocatable :: acof(:,:,:)
      complex, allocatable :: bcof(:,:,:)
      complex, allocatable :: ccof(:,:,:,:)
!!! the parameters for the number of wfs
      integer :: nwfs
!!! the potential in the spheres and the vacuum
      real, allocatable :: vr(:,:,:),vz(:,:,:),vrf(:,:,:,:)
!!! auxiliary potentials
      complex, allocatable :: vpw(:,:),vzxy(:,:,:,:)
!!! bkpts data
      integer nntot,ikpt_help
      integer, allocatable :: gb(:,:,:),bpt(:,:)
!!! radial wavefunctions in the muffin-tins and more ...
      real,    allocatable :: flo(:,:,:,:)
      real,    allocatable :: ff(:,:,:,:),gg(:,:,:,:)
      real,    allocatable :: us(:,:),uds(:,:),ulos(:,:)
      real,    allocatable :: dus(:,:),duds(:,:),dulos(:,:)
      real,    allocatable :: ddn(:,:),uulon(:,:),dulon(:,:)
      real,    allocatable :: uloulopn(:,:,:)
      real     uuilon(nlod,ntypd),duilon(nlod,ntypd)
      real     ulouilopn(nlod,nlod,ntypd)
!!! energy parameters
      real    ello(nlod,ntypd,jspd),evac(2,jspd)
      real    epar(0:lmaxd,ntypd,jspd),evdu(2,jspd)
      character(len=3) :: spin12(2)
      data spin12/'WF1' , 'WF2'/
      character(len=2)::spinspin12(2)
      data spinspin12/'.1','.2'/
      complex,allocatable::wannierfunc(:,:)
      complex,allocatable::wannierfunc_temp(:,:)
      character(len=33):: header
      integer num_nnmax
      integer posi
      complex,allocatable::u_matrix_tmp(:,:,:)
      complex,allocatable::u_matrix(:,:,:)
      real tmp_omi
      integer kpt,oper
      real poinint(3)
      real phas,tmax
      real bkrot(3)
      integer j1,j2,j3
      logical um_format
      logical have_disentangled,l_chk,l_umdat
      integer,allocatable :: ndimwin(:)
      logical,allocatable :: lwindow(:,:)
      integer :: chk_unit,nkp,ntmp,ierr,fullnkpts
      integer,allocatable::irreduc(:),mapkoper(:)
      character(len=20)::checkpoint
      real :: tmp_latt(3,3)
      real,allocatable:: tmp_kpt_latt(:,:)
      real omega_invariant
      complex,allocatable::u_matrix_opt(:,:,:)
      logical l_file
      logical,allocatable::inc_band(:)
      integer num_inc,counter,kptibz
      logical l_byindex, l_byenergy, l_bynumber
      integer num_wann,num_bands,kpun
      complex d_wgn(-3:3,-3:3,3,nop)
      integer pos1,pos2,ato,loc,invop
c     ..basis wavefunctions in the vacuum
      complex, allocatable :: ac(:,:,:),bc(:,:,:)
      real,    allocatable :: dt(:),dte(:)
      real,    allocatable :: t(:),te(:),tei(:)
      real,    allocatable :: u(:,:,:),ue(:,:,:)
      real vz0(2)
      integer ik,nv2,ivac,jvac,symvac,symvacvac
      real evacp,sign,arg
      complex c_1
      integer kvac1(nv2d),kvac2(nv2d),map2(nvd)
      real fas,zks
      integer mesh
      integer n2
      real v(3),scale,ev
      complex av,bv
      real volume
c      external dotirp !module now
c      real dotirp
      REAL          :: s,const
      COMPLEX       :: xdnout,factor
      INTEGER       :: ii3,ix,iy,iz,nplo,nbn
      INTEGER       :: nbmin,nbmax
      INTEGER       :: nplot,nq,nt,jm,ii1,ii2
      LOGICAL       :: twodim
      real,allocatable::knorm(:,:)
      real,allocatable::wfnorm(:)
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      INTEGER :: grid(3),k      
      integer,allocatable :: shiftkpt(:,:)
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(jmtd),realpart,imagpart
      REAL    :: point(3),post(3,natd)
      CHARACTER(len=30):: filename
      CHARACTER(len=20):: name1,name2,name3
      CHARACTER(len=10):: vandername
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename 
      complex :: nsfactor
      integer :: ngopr1(natd)

      um_format=.false.
      l_byindex=.false.
      l_byenergy=.false.
      l_bynumber=.false.
      if(sortrule==1)l_byindex=.true.
      if(sortrule==2)l_bynumber=.true.
      if(sortrule==3)l_byenergy=.true.
      ngopr1(:)=1



c     read in plot_inp

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
     +        ,calledby ="wann_plot_um_dat")
      ENDIF

      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)') nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (nplot.ge.2) 
     &     CALL juDFT_error
     +     ("plots one by one, please, this is not charge density"
     +     ,calledby ="wann_plot_um_dat")
      twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
      vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
      zero = (/0.,0.,0./);filename="default"
      READ(18,plot)
      IF (twodim.AND.ANY(grid(1:2)<1)) 
     +     CALL juDFT_error("Illegal grid size in plot",calledby
     +     ="wann_plot_um_dat")
      IF (.NOT.twodim.AND.ANY(grid<1)) 
     +     CALL juDFT_error("Illegal grid size in plot",calledby
     +     ="wann_plot_um_dat")
      IF (twodim) grid(3) = 1
      !calculate cartesian coordinates if needed
      IF (.NOT.cartesian) THEN
         vec1=matmul(amat,vec1)
         vec2=matmul(amat,vec2)
         vec3=matmul(amat,vec3)
         zero=matmul(amat,zero)
      ENDIF
      Close(18)

      !calculate volume
      volume  = vec1(1)*vec2(2)*vec3(3) + vec2(1)*vec3(2)*vec1(3) +
     &          vec3(1)*vec1(2)*vec2(3) - vec1(3)*vec2(2)*vec3(1) -
     &          vec2(3)*vec3(2)*vec1(1) - vec3(3)*vec1(2)*vec2(1)
    


      sfp = 2* sqrt( pimach() )
      tpi = 2* pimach()
      lmd = lmaxd*(lmaxd+2)
      nkpts = maxval(nkpt(:))

      nrec = 0
      nlotot = 0
      do n = 1, ntype
        do l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        enddo
      enddo

cccccccccccccccc   initialize the potential   cccccccccccc

      allocate ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
      allocate ( vz(nmzd,2,jspd),vrf(jmtd,0:nlhd,ntypd,jspd) )
      allocate ( vr(jmtd,ntypd,jspd) )

      open (8,file='pottot',form='unformatted',status='old')
      rewind (8)

      call loddop(
     >      jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >      jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >      nlh,jri,ntypsd,ntypsy,8,natd,neq,
     <      iop,dop,iter,vrf,vpw,vz,vzxy,name)

      close (8)

      do jspin = 1,jspins
        do n = 1, ntype
          do j = 1,jri(n)
            vr(j,n,jspin) = vrf(j,0,n,jspin)
          enddo
        enddo
      enddo

      deallocate ( vpw,vzxy,vrf )
 
cccccccccccccccc   end of the potential part  ccccccccccc

      allocate ( kveclo(nlotot),nv(jspd) )
      allocate ( k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd) )
      allocate ( ff(ntypd,jmtd,2,0:lmaxd) )
      allocate ( gg(ntypd,jmtd,2,0:lmaxd) )
      allocate ( us(0:lmaxd,ntypd) )
      allocate ( uds(0:lmaxd,ntypd) )
      allocate ( dus(0:lmaxd,ntypd) )
      allocate ( duds(0:lmaxd,ntypd) )
      allocate ( ddn(0:lmaxd,ntypd) )
      allocate ( ulos(nlod,ntypd) )
      allocate ( dulos(nlod,ntypd) )
      allocate ( uulon(nlod,ntypd) )
      allocate ( dulon(nlod,ntypd) )
      allocate ( uloulopn(nlod,nlod,ntypd) )

c cycle by spins starts! NO NON-COLLINEAR WHATSOEVER!

      do 110 jspin=1,jspd   ! cycle by spins
         print*,"spin=",jspin
      jsp_start = jspin ; jsp_end = jspin

c*******************************************************
c      get num_bands and num_wann from WF1.amn (WF2.amn)
c*******************************************************
      l_file=.false.
      inquire(file=spin12(jspin)//'.amn',exist=l_file)
      open(355,file=spin12(jspin)//'.amn')
      read(355,*)
      read(355,*)num_bands,kpun,num_wann
      close(355)
      if(l_byindex.and..not.((1+band_max(jspin)-
     &  band_min(jspin)).eq.num_bands))
     &     CALL juDFT_error("1+band_max-band_min  /=  num_bands",calledby
     +     ="wann_plot_um_dat")

c**************************************************************
!   for bzsym = .true.: determine mapping between kpts and w90kpts
c**************************************************************
      if (l_bzsym) then
         l_file=.false.
         inquire(file='w90kpts',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("w90kpts not found, needed if bzsym",calledby
     +        ="wann_plot_um_dat")
         open(412,file='w90kpts',form='formatted')
         read(412,*)fullnkpts
         close(412)
         print*,"fullnkpts=",fullnkpts
         IF(fullnkpts<=nkpts) CALL juDFT_error("fullnkpts.le.nkpts"
     +        ,calledby ="wann_plot_um_dat")
         allocate(irreduc(fullnkpts),mapkoper(fullnkpts))
         allocate(shiftkpt(3,fullnkpts))
         l_file=.false.
         inquire(file='kptsmap',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("kptsmap not found, needed if bzsym",calledby
     +        ="wann_plot_um_dat")
         open(713,file='kptsmap')
         do i=1,fullnkpts
            read(713,*)kpt,irreduc(i),mapkoper(i),shiftkpt(:,i)
            IF(kpt/=i) CALL juDFT_error("kpt.ne.i",calledby
     +           ="wann_plot_um_dat")
            print*,i,irreduc(i),mapkoper(i)
         enddo   
         close(713)
         IF(MAXVAL(irreduc(:))/=nkpts) CALL juDFT_error
     +        ("max(irreduc(:))/=nkpts",calledby ="wann_plot_um_dat")
      else   
         fullnkpts=nkpts
      endif   

      IF(kpun/=fullnkpts) CALL juDFT_error
     +     ("mismatch in kpun and fullnkpts",calledby
     +     ="wann_plot_um_dat")

      if(.not.l_proj_plot)then
c**************************************************************
c        read in chk
c*************************************************************
         allocate( u_matrix(num_bands,num_wann,fullnkpts) )
         allocate( lwindow(num_bands,fullnkpts) )
         allocate( ndimwin(fullnkpts) )
         call wann_read_umatrix(
     >            fullnkpts,num_wann,num_bands,
     >            um_format,jspin,wan90version,
     <            have_disentangled,
     <            lwindow,ndimwin,
     <            u_matrix)
      else
c**************************************************************
c        read WF1.umn (WF2.umn) (if projmethod)
c**************************************************************
         have_disentangled=.false.
         l_file=.false.
         inquire(file=spin12(jspin)//'.umn',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error("no umn file foun
     +d",calledby ="wann_plot_um_dat")
         open(419,file=spin12(jspin)//'.umn')
         read(419,*)     !num_wann,num_bands
         allocate(u_matrix(num_bands,num_wann,fullnkpts))
         do ikpt=1,fullnkpts
            do j=1,num_wann
               do i=1,num_bands
                  read(419,*)idum,jdum,kdum,realpart,imagpart
                   u_matrix(i,j,ikpt)=cmplx(realpart,imagpart)
               enddo
            enddo   
         enddo   
         close(419)
      endif   

c      if(um_format)then
c         open(419,file='umatrix_formatted')
c         do ikpt=1,fullnkpts
c            do j=1,num_wann
c               do i=1,num_bands
c                  write(419,*)u_matrix(i,j,ikpt)
c               enddo
c            enddo   
c         enddo   
c         close(419)
c      endif


***********************************************************
***********************************************************

      print*,"num_wann=",num_wann
      print*,"num_bands=",num_bands
      allocate(wannierfunc(num_wann,
     &  (grid(1))*(grid(2))*(grid(3))))



      wannierfunc(:,:)=cmplx(0.0,0.0)


cccccccccccc   read in the eigenvalues and vectors   cccccc

      l_p0 = .false.
      if (irank.eq.0) l_p0 = .true.

      call cdn_read0(
     >               lmaxd,ntypd,nlod,neigd,jspd,
     >               irank,isize,jspin,jsp_start,jsp_end,
     >               l_noco,nrec,66,
     <               ello,evac,epar,bkpt,wk,n_bands,nrec1,n_size)


      allocate ( flo(ntypd,jmtd,2,nlod) )

      na = 1
      do 40 n = 1,ntype
       do 30 l = 0,lmax(n)
c...compute the l-dependent, k-independent radial MT- basis functions
         call radfun(
     >              l,epar(l,n,jspin),vr(1,n,jspin),jri(n),rmsh(1,n),
     >              dx(n),jmtd,
     <              ff(n,:,:,l),gg(n,:,:,l),us(l,n),
     <              dus(l,n),uds(l,n),duds(l,n),
     <              ddn(l,n),nodeu,noded,wronk)
   30  continue
c...and the local orbital radial functions
c       do ilo = 1, nlo(n)
         call radflo(
     >             ntypd,nlod,jspd,jmtd,lmaxd,n,jspin,
     >             ello(1,1,jspin),vr(1,n,jspin),
     >             jri(n),rmsh(1,n),dx(n),ff(n,1:,1:,0:),
     >             gg(n,1:,1:,0:),llo,nlo,l_dulo(1,n),irank,ulo_der,
     <             ulos(1,1),dulos(1,1),uulon(1,1),dulon(1,1),
     <             uloulopn(1,1,1),uuilon,duilon,ulouilopn,flo(n,:,:,:))
c       enddo
c       na = na + neq(n)
   40 continue
      i_rec = 0 ; n_rank = 0

c******************************************************************
c          beginning of k-point loop,each may be a separate task
c******************************************************************

      allocate(knorm(fullnkpts,num_bands))
      print*,"num_bands=",num_bands
      knorm(:,:)=0.0

      do ikpt = 1,fullnkpts  ! loop by k-points starts

        i_rec = i_rec + 1
        if (mod(i_rec-1,isize).eq.irank) then
        print*,"k-point=",ikpt
        kptibz=ikpt 
        if(l_bzsym) kptibz=irreduc(ikpt)
        if(l_bzsym) oper=mapkoper(ikpt)

       if(have_disentangled) then
          allocate(inc_band(size(lwindow,1)))
          inc_band(:)=lwindow(:,ikpt)
          num_inc=ndimwin(ikpt)
       end if

      allocate ( z(nbasfcn,neigd),we(neigd),
     &                   zz(nbasfcn,neigd),eigg(neigd) )



      call wann_read_eig(66,     
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eigg,zz,cp_time)


      nslibd = 0

c...we work only within the energy window

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      z(:,:) = cmplx(0.,0.)
#else
      z(:,:) = 0.
#endif
      eig(:) = 0.

      print*,"bands used"
      do i = 1,nbands
       if ((eigg(i).ge.e1s .and. nslibd.lt.num_bands.and.l_bynumber)
     &.or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.l_byenergy)
     &.or.(i.ge.band_min(jspin).and.i.le.band_max(jspin)
     &.and.l_byindex))then
          print*,i
        nslibd = nslibd + 1
        eig(nslibd) = eigg(i)
        we(nslibd) = we(i)
        do  j = 1,nv(jspin) + nlotot
         z(j,nslibd) = zz(j,i)
        enddo
       endif 
      enddo 
c***********************************************************
c              rotate the wavefunction
c***********************************************************
c      if (l_bzsym.and.oper.ne.1) then  !rotate bkpt
c         bkrot(:)=0.0
c         do k=1,3
c           bkrot(:)=bkrot(:)+mrot(k,:,oper)*bkpt(k)
c         enddo
c         bkpt(:)=bkrot(:)

c         jloop:do j=1,nv(jspin)
c               j1=mrot(1,1,oper)*k1(j,jspin)+
c     +             mrot(2,1,oper)*k2(j,jspin)+mrot(3,1,oper)*k3(j,jspin)
c               j2=mrot(1,2,oper)*k1(j,jspin)+
c     +             mrot(2,2,oper)*k2(j,jspin)+mrot(3,2,oper)*k3(j,jspin)  
c               j3=mrot(1,3,oper)*k1(j,jspin)+
c     +             mrot(2,3,oper)*k2(j,jspin)+mrot(3,3,oper)*k3(j,jspin)
c               k1(j,jspin)=j1
c               k2(j,jspin)=j2
c               k3(j,jspin)=j3
c         enddo jloop   

c      endif   

      if (l_bzsym.and.oper.ne.1) then  !rotate bkpt
         call wann_kptsrotate(
     >            natd,nlod,llod,
     >            ntypd,nlo,llo,invsat,
     >            l_noco,l_soc,
     >            ntype,neq,nlotot,
     >            kveclo,jspin,
     >            oper,nop,mrot,nvd,nv,
     >            shiftkpt(:,ikpt),
     >            tau,
     x            bkpt,k1(:,:),k2(:,:),k3(:,:),
     x            z,nsfactor)
      endif

      print*,"bkpt=",bkpt(:)
c******************************************************************


      noccbd = nslibd

      allocate ( acof(noccbd,0:lmd,natd),
     &           bcof(noccbd,0:lmd,natd),
     &           ccof(-llod:llod,noccbd,nlod,natd))

      acof(:,:,:) = cmplx(0.,0.) ; bcof(:,:,:) = cmplx(0.,0.)
      ccof(:,:,:,:) = cmplx(0.,0.)

c...generation of the A,B,C coefficients in the spheres 
c...for the lapws and local orbitals, summed by the basis functions 

      call abcof(
     >           lmaxd,ntypd,neigd,noccbd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngopr1,taual,neq,lmax,rmt,omtil,
     >           bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,noccbd,z,
     >           us(0,1),dus(0,1),uds(0,1),
     >           duds(0,1),ddn(0,1),invsat,invsatnr,
     >           ulos(1,1),uulon(1,1),dulon(1,1),
     >           dulos(1,1),llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,jspin,alph,beta,qss,kveclo,odi,ods,
     <           acof(1,0,1),bcof(1,0,1),
     <           ccof(-llod,1,1,1))

c      call abcrot(
c     >        ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
c     >        noccbd,lmax,nlo,llo,nop,ngopr,mrot,invsat,invsatnr,
c     >                 bmat,odi,ods,
c     X                 acof,bcof,ccof)

      call wann_abinv(
     >        ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
     >        noccbd,lmax,nlo,llo,invsat,invsatnr,bkpt,taual,
     X        acof,bcof,ccof)


c***********************************************************************
c make preparations for plotting in vacuum
c***********************************************************************
      if (film)then
         allocate ( ac(nv2d,nslibd,2),bc(nv2d,nslibd,2),
     +         u(nmzd,nv2d,nvac),ue(nmzd,nv2d,nvac))
         call wann_2dvacabcof(
     >         nv2d,nslibd,nvac,nmzd,nmz,omtil,vz(:,:,jspin),
     >         nv(jspin),bkpt,z1,
     >         nvd,k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
     >         bbmat,delz,bmat,nbasfcn,
     >         neigd,z,
     <         ac,bc,u,ue)
      endif !preparations for vacuum
c**************************************************************************
c**************************************************************************

      nbmin=1
      nbmax=nslibd
      counter=1

      band:DO nbn = nbmin,nbmax
       if(have_disentangled) then
             if(counter>num_inc) exit
             if(.not.inc_band(nbn))cycle band
       endif
          

          DO iz = 0,grid(3)-1
          DO iy = 0,grid(2)-1
           xloop:DO ix = 0,grid(1)-1
            posi=ix+1+iy*(grid(1))+iz*(grid(1))*(grid(2))
            point = zero+vec1*(ix+0.0)/(grid(1)-1)+vec2*(iy+0.0)
     $                 /(grid(2)-1)
            IF (.NOT.twodim) point = point+vec3*(iz+0.0)/(grid(3)-1)
            call cotra1(point,poinint,bmat)
            phas=tpi*(bkpt(1)*poinint(1)
     +      +bkpt(2)*poinint(2)+bkpt(3)*poinint(3))
            factor=cmplx(cos(phas),sin(phas))
!Check if the point is in MT-sphere
             ii1 = 3
             ii2 = 3
             ii3 = 3
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
                  s  = SQRT(dot_PRODUCT(pos(:,na)-pt,pos(:,na)-pt))
                  IF (s<rmsh(jri(nt),nt)) THEN
                    CALL wann_real(
     >                   pt,nt,na,0,1,bkpt,
     >                   n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >                   natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >                   nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >                   lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >                   omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >                   ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
     >                   ccof(:,nbn,:,:),z(:,nbn),
     >               nv(jspin),k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >                   lmd,nbasfcn,
     <                   xdnout)
                    wannierfunc(:,posi)=
     =   wannierfunc(:,posi)+xdnout*u_matrix(counter,:,ikpt)*factor
                    knorm(ikpt,nbn)=knorm(ikpt,nbn)+(abs(xdnout))**2
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
                call wann_plot_vac(point,z1,nmzd,nv2d,n3d,nvac,
     >            nmz,delz,bmat,bbmat,evac(:,jspin),bkpt,vz,jspin,
     >            k1(:,jspin),k2(:,jspin),k3(:,jspin),nvd, 
     >            nbasfcn,neigd,z(:,nbn),nv(jspin),omtil,nslibd,
     >            ac(:,nbn,ivac),
     &             bc(:,nbn,ivac),
     &            u(:,:,jvac),ue(:,:,jvac),xdnout)

               
c               CALL wann_real(
c     >              point,0,0,1,0,bkpt,
c     >              n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
c     >              natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
c     >              nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
c     >              lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
c     >              omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
c     >              ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
c     >              ccof(:,nbn,:,:),z(:,nbn),
c     >              nv,k1,k2,k3,lmd,nbasfcn,
c     <              xdnout)
                wannierfunc(:,posi)=
     =             wannierfunc(:,posi)+
     +                  xdnout*u_matrix(counter,:,ikpt)*factor
              CYCLE xloop
             END IF
            
             IF (odi%d1) THEN
              IF (SQRT((pt(1))**2 + (pt(2))**2)>=z1) THEN
                CALL wann_real(
     >               pt,0,0,1,0,bkpt,
     >               n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >               natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >               nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >               lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >               omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >               ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
     >               ccof(:,nbn,:,:),z(:,nbn),
     >               nv(jspin),k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >                lmd,nbasfcn,
     <               xdnout)
                wannierfunc(:,posi)=
     =wannierfunc(:,posi)+xdnout*u_matrix(counter,:,ikpt)*factor
                CYCLE xloop
              END IF
             END IF
             CALL wann_real(
     >             point,0,0,0,2,bkpt,
     >             n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >             natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >             nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >             lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >             omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >             ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
     >             ccof(:,nbn,:,:),z(:,nbn),
     >               nv(jspin),k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >             lmd,nbasfcn,
     <             xdnout)
             wannierfunc(:,posi)=
     =wannierfunc(:,posi)+xdnout*u_matrix(counter,:,ikpt)*factor
                    knorm(ikpt,nbn)=knorm(ikpt,nbn)+(abs(xdnout))**2
            ENDDO xloop
           ENDDO
          ENDDO !z-loop

c..end of the loop by the bands
      counter=counter+1


      ENDDO band   

      deallocate ( acof,bcof,ccof,z,we,zz,eigg )

      write (*,*) 'nslibd=',nslibd


      if(film)then
         deallocate(ac,bc,u,ue)
      endif   


      endif!processors

      enddo !loop over k-points
      mesh=grid(1)*grid(2)*grid(3)

#ifdef CPP_MPI
c      call MPI_BARRIER(mpi_comm,ierr)
      if(l_p0)then
         if(isize.ne.1)then
       allocate(wannierfunc_temp(num_wann,mesh))  
       do cpu_index=1,isize-1
        do ikpt=1,fullnkpts
         if(mod(ikpt-1,isize).eq.cpu_index)then
           call MPI_RECV(knorm(ikpt,1:num_bands),num_bands,
     &             CPP_MPI_REAL,cpu_index,
     &        ikpt,mpi_comm,stt,mpiierr)
         endif !processors
        enddo !ikpt
           call MPI_RECV(wannierfunc_temp(1:num_wann,1:mesh),
     &         num_wann*mesh,
     &             CPP_MPI_COMPLEX,cpu_index,
     &        cpu_index+fullnkpts,mpi_comm,stt,mpiierr)
           wannierfunc(:,:)=wannierfunc(:,:)+
     &                   wannierfunc_temp(:,:)



       enddo !cpu_index 
       deallocate(wannierfunc_temp)
       endif !isize
      else
       do ikpt=1,fullnkpts
        if(mod(ikpt-1,isize).eq.irank)then
             call MPI_SEND(knorm(ikpt,1:num_bands),num_bands,
     &           CPP_MPI_REAL,0,ikpt,mpi_comm,mpiierr)
        endif !processors  
       enddo !ikpt  
             call MPI_SEND(wannierfunc(1:num_wann,1:mesh),
     &        num_wann*mesh,
     & CPP_MPI_COMPLEX,0,fullnkpts+irank,mpi_comm,mpiierr)



      endif ! l_p0   



#endif


      deallocate(flo)
      if(l_p0)then
      wannierfunc(:,:)=wannierfunc(:,:)/real(fullnkpts)
      DO nplo=1,num_wann


c****************************************************************
c      make Wannier function real (as much as possible)
c****************************************************************
       phas=0.0
       do iz=0,grid(3)-1
          do iy=0,grid(2)-1
             do ix=0,grid(1)-1
                posi=ix+1+iy*(grid(1))+iz*(grid(1))*(grid(2))
               tmax=wannierfunc(nplo,posi)*conjg(wannierfunc(nplo,posi))

                if (tmax>phas) then
                   phas=tmax
                   factor=wannierfunc(nplo,posi)
                end if
             end do
          end do
       end do
       factor=factor/sqrt(real(factor)**2+aimag(factor)**2)
       wannierfunc(nplo,:)=wannierfunc(nplo,:)/factor


c***************************************************************
c       open files for plot and make headers
c***************************************************************
         IF (xsf) THEN
            write (name1,22) nplo,jspin
   22       format (i3.3,'.real.',i1,'.xsf')
            write (name2,23) nplo,jspin
   23       format (i3.3,'.imag.',i1,'.xsf')
            write (name3,24) nplo,jspin
   24       format (i3.3,'.absv.',i1,'.xsf')
            OPEN(55,file=name1)
            CALL xsf_WRITE_atoms(
     >                        55,film,odi%d1,amat,neq(:ntype),
     >                        zatom(:ntype),pos)
            OPEN(56,file=name2)
            CALL xsf_WRITE_atoms(
     >                        56,film,odi%d1,amat,neq(:ntype),
     >                        zatom(:ntype),pos)
            OPEN(57,file=name3)
            CALL xsf_WRITE_atoms(
     >                        57,film,odi%d1,amat,neq(:ntype),
     >                        zatom(:ntype),pos)
            CALL xsf_WRITE_header(55,twodim,filename,vec1,vec2,vec3,zero
     $           ,grid)
            CALL xsf_WRITE_header(56,twodim,filename,vec1,vec2,vec3,zero
     $           ,grid)
            CALL xsf_WRITE_header(57,twodim,filename,vec1,vec2,vec3,zero
     $           ,grid)
         ELSE
               WRITE (vandername,201) nplo,jspin
  201          FORMAT (i5.5,'.',i1)            
               OPEN(55,file=vandername)
               WRITE (55,7) grid(1),grid(2),grid(3),ikpt,nslibd
    7          FORMAT (5i4)
         ENDIF
c********************************************************************
c        write data to files
c********************************************************************
         DO iz = 0,grid(3)-1
          DO iy = 0,grid(2)-1
           DO ix = 0,grid(1)-1
              posi=ix+1+iy*grid(1)+iz*grid(1)*grid(2)
              IF (xsf) THEN
                 WRITE(55,*) real(wannierfunc(nplo,posi))
                 WRITE(56,*) aimag(wannierfunc(nplo,posi))
                 WRITE(57,*) abs(wannierfunc(nplo,posi))
              ELSE
                 WRITE(55,8) real(wannierfunc(nplo,posi))
              ENDIF
           enddo
          enddo
         enddo 
         IF (xsf) THEN
              CALL xsf_WRITE_endblock(55,twodim)
              CALL xsf_WRITE_endblock(56,twodim)
              CALL xsf_WRITE_endblock(57,twodim)
              CLOSE (55) ; CLOSE (56) ; CLOSE (57)
         ENDIF
           
      ENDDO   !nplo      
      IF (.not.xsf) CLOSE(55)
    8 FORMAT (2f7.3)

c*******************************************************************
c     determine spreads, centers, norms
c*******************************************************************
      
      l_spreadcal=.true.
      if(l_spreadcal)then !calculate spreads and centers from real space grid
         print*,"calculate spreads and centers"
         allocate (spreads(num_wann,num_wann))
         allocate (centers(3,num_wann))
         allocate(wfnorm(num_wann))
         centers(:,:)=0.0
         wfnorm(:)=0.0
         spreads(:,:)=cmplx(0.0,0.0)
         do nplo=1,num_wann
           do iz=0,grid(3)-1
            do iy=0,grid(2)-1
               do ix=0,grid(1)-1
                  posi=ix+1+iy*(grid(1))+iz*(grid(1))*(grid(2))
                  point = zero+vec1*(ix+0.0)/(grid(1)-1)+vec2*(iy+0.0)
     $                 /(grid(2)-1)
           IF (.NOT.twodim) point = point+vec3*(iz+0.0)/(grid(3)-1)
                  centers(:,nplo)=centers(:,nplo)
     +               +point(:)*(abs(wannierfunc(nplo,posi)))**2
                  wfnorm(nplo)=wfnorm(nplo)+
     +         (abs(wannierfunc(nplo,posi)))**2
                  do ii1=1,num_wann
                     spreads(nplo,ii1)=spreads(nplo,ii1)
     +  +wannierfunc(nplo,posi)*conjg(wannierfunc(ii1,posi))*
     *       dot_product(point,point)   
                  enddo   
               enddo
            enddo
           enddo   
        enddo

        do nplo=1,num_wann !normalize centers
           centers(:,nplo)=centers(:,nplo)/(mesh)*volume
        enddo   
 
        do nplo=1,num_wann  !normalize spreads
           do ii1=1,num_wann
             spreads(nplo,ii1)=spreads(nplo,ii1)/mesh*volume
           enddo        
           spreads(nplo,nplo)=spreads(nplo,nplo)
     &           -dot_product(centers(1:3,nplo),centers(1:3,nplo))      
        enddo   

        wfnorm(:)=wfnorm(:)/(mesh)*volume !normalize wfnorm

         knorm(:,:)=knorm(:,:)/mesh*volume !normalize knorm       

c***********************************************************
c        write spreads and so on to files
c***********************************************************

         open(518,file=spin12(jspin)//'.centers')
         do nplo=1,num_wann
            write(518,*)centers(:,nplo)
         enddo
         close(518)
         open(519,file=spin12(jspin)//'.spreads')
         do nplo=1,num_wann
            do ii1=1,num_wann
               write(519,*)nplo,ii1,spreads(nplo,ii1)
            enddo   
         enddo   
         close(519)
         open(521,file=spin12(jspin)//'.norm')
         do ii1=1,num_wann
            write(521,*)wfnorm(ii1)
         enddo   
         close(521)
         deallocate(centers)
         deallocate(spreads)
         deallocate(wfnorm)
      endif   

      open(611,file=spin12(jspin)//'.knorm')
      do nplo=1,num_bands
         do ikpt=1,fullnkpts
            write(611,*)ikpt,nplo,knorm(ikpt,nplo)
         enddo   
      enddo
      close(611)

c*****************************************************************
c*****************************************************************

      endif !l_p0

      deallocate(knorm)
      if(have_disentangled)deallocate(inc_band)
      if(.not.l_proj_plot)then
         deallocate(lwindow,ndimwin)
      endif
      deallocate(wannierfunc)
      nrec=nrec+nkpts
      deallocate(u_matrix)

#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,mpiierr)
#endif      

110   continue ! end of cycle by spins


      deallocate ( vr,vz,kveclo,nv,k1,k2,k3 )
      deallocate ( ff,gg,us,dus,duds,uds,ddn)
      deallocate ( ulos,dulos,uulon,dulon,uloulopn )
 
#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,mpiierr)
#endif      

      END SUBROUTINE wann_plot_um_dat
      END MODULE m_wann_plot_um_dat
