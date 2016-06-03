      MODULE m_wannier
      use m_juDFT
      CONTAINS
      SUBROUTINE wannier(
     >      mpi_comm,l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,l_soc,
     >      memd,lnonsph,clnu,lmplmd,mlh,nmem,llh,lo1l,
     >      theta,phi,soc_opt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Makes necessary for the construction of the wannier functions
c     (W90: Yates, Mostofi, Marzari, Souza, Vanderbilt '06 f90 code)
c     ab initio preliminaries: constructs the overlaps of the periodic
c     parts of the wavefunctions and the projections of the 
c     wavefunctions
c     onto a set of starting wfs, i.e. atomic-like orbitals.
c                                                            YM 06
c     Mmn(k,b) = <u_{nk}|u_{m(k+b)}>, u being a periodic part
c                        of the wavefunction psi_nk
c     A_mn^k = <psi_mk|g_n>, where g_n is a trial orbital
c     which are written into the files 'WF1.mmn' and 'WF1.amn'
c           Marzari Vanderbilt PRB 56,12847(1997)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parallelization, Optionals, Symmetry, Noco&Soc:
c     Frank Freimuth
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The routine reads the bkpts file, which contains the following
c     information:
c     1st line: nntot (INT) - number of the nearest neighbors for
c                             each k-point in the MP mesh
c     2-nkpts*nntot lines containing 5 integers i1,i2,i3,i4,i5:
c     i1 - the number of the k-point in the kpts file
c     i2 - number of the k-point, which is a periodic image of
c          k+b in the 1st BZ
c     i3-i5 - coordinates of the G-vector, connecting k-point
c             i2 with the actual k+b k-point
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     In general, the number of bands for each k-point is
c     different: N(k), and also differs from the number of bands
c     we are interested in: N (for instance 5 d-bands of Cu among
c     the 6 s- and d-bands). While matrices Mmn(k) are square
c     for each k-point, matrices Mmn(k,b) can be made so after
c     defining the maximum number of bands max(N(k)).
c     The matrix Amn is non-diagonal by default (N(k)*N).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Total number of wannier functions: nwfs
c
c     e1s,e2s: lower and upper boundaries of the energy window:
c     Needed for sorting by number and sorting by energy.
c     Not needed for sorting by index.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use m_types
      use m_wann_mmnk_symm
      use m_wann_rw_eig
      use m_abcof
      use m_radfun
      use m_radflo
      use m_cdnread, only : cdn_read0, cdn_read
      use m_od_types, only : od_inp, od_sym
      use m_loddop
      use m_constants, only : pimach
      use m_wann_mmk0_od_vac
      use m_wann_mmkb_od_vac
      use m_wann_mmk0_vac
      use m_wann_mmkb_vac
      use m_wann_updown
      use m_wann_mmk0_sph
      use m_wann_ujugaunt
      use m_wann_mmkb_sph
      use m_wann_projmethod
      use m_wann_amn
      use m_wann_abinv
      use m_wann_kptsrotate
      use m_wann_plot
      use m_wann_read_inp
      use m_wann_plot_symm
      use m_wann_mmkb_int
      use m_wann_postproc
      use m_matmul,only : matmul3,matmul3r
      use m_wann_write_mmnk
      use m_wann_write_amn
      use m_wann_write_nabla
      use m_vsoc
      use m_wann_write_matrix4
      use m_wann_write_matrix5
      use m_wann_orbcomp
      use m_wann_anglmom
#ifdef CPP_TOPO
      use m_wann_surfcurr
      use m_wann_surfcurr_int2
      use m_wann_nabla
      use m_wann_nabla_vac
      use m_wann_soc_to_mom
#endif

      IMPLICIT NONE
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer ierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)

#endif
      logical, intent (in) :: invs,invs2,film,slice,symor
      integer, intent (in) :: lmaxd,ntypd,neigd,nkptd,kk,nnne,mpi_comm
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
      integer, intent (in) :: invsat(natd),invsatnr(natd),nkpt(nwdd)
      integer, intent (in) :: irank,isize,nv2d,nmzxy,nmz
      integer, intent (in) :: ulo_der(nlod,ntypd),ig2(n3d)
      real,    intent (in) :: taual(3,natd),rmt(ntypd),dx(ntypd)
      real,    intent (in) :: amat(3,3),bmat(3,3),bbmat(3,3)
      real,    intent (in) :: rmsh(jmtd,ntypd),tau(3,nop),zatom(ntype)
      real, intent (inout) :: alph(ntypd),beta(ntypd),qss(3)
      real,    intent (in) :: pos(3,natd)
      real, intent(in)     :: ef
      complex, intent (in) :: ustep(n3d)
      logical, intent (in) :: l_dulo(nlod,ntypd),l_noco,l_ss
      logical, intent (in) :: l_soc
      integer, intent (in) :: memd
      integer, intent (in) :: lnonsph(ntypd)
      complex, intent (in) :: clnu(memd,0:nlhd,ntypsd)
      integer, intent (in) :: lmplmd
      integer, intent (in) :: mlh(memd,0:nlhd,ntypsd)
      integer, intent (in) :: nmem(0:nlhd,ntypsd)
      integer, intent (in) :: llh(0:nlhd,ntypsd)
      integer, intent (in) :: lo1l(0:llod,ntypd)
      real,    intent (in) :: sk2(n2d),phi2(n2d)
      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods
      real, intent(in) :: theta,phi
      logical,intent(in) :: soc_opt(ntype+2)

cccccccccccccccccc   local variables   cccccccccccccccccccc
      type(t_wann) :: wann
      integer :: lmd,nlotot,n,nmat,iter,ikpt,ikpt_b,nmat_b
      integer :: addnoco,funbas,loplod,addnoco2,igvm2
      integer :: noccbd,noccbd_b,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
      integer :: jsp_start,jsp_end,nrec,nrec1,nbands,nbands_b
      integer :: nodeu,noded,n_size,na,n_rank,nbnd,numbands
      integer :: i1,i2,i3,in,lda
      integer :: n_bands(0:neigd),nslibd,nslibd_b
      character(len=8) :: dop,iop,name(10)
      real    :: bkpt(3),bkpt_b(3),sfp,tpi,wronk,wk,wk_b,phase
      complex :: c_phase
      real    :: eig(neigd),eig_b(neigd),cp_time(9),efermi
      logical :: l_p0,l_bkpts,l_proj,l_amn,l_mmn
!!! energy window boundaries
      integer, allocatable :: kveclo(:),nv(:)
      integer, allocatable :: kveclo_b(:),nv_b(:)
      integer, allocatable :: k1(:,:),k2(:,:),k3(:,:)
      integer, allocatable :: k1_b(:,:),k2_b(:,:),k3_b(:,:)
      real,    allocatable :: we(:),we_b(:)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex, allocatable :: z(:,:),zz(:,:)
      complex, allocatable :: z_b(:,:)
      complex, allocatable :: ztrans(:,:)
#else
      real,    allocatable :: z(:,:),zz(:,:)
      real,    allocatable :: z_b(:,:)
      real,    allocatable :: ztrans(:,:)
#endif
      real,    allocatable :: eigg(:)
      real kpoints(nkptd)
!!! a and b coeff. constructed for each k-point
      complex, allocatable :: acof(:,:,:),acof_b(:,:,:)
      complex, allocatable :: bcof(:,:,:),bcof_b(:,:,:)
      complex, allocatable :: ccof(:,:,:,:),ccof_b(:,:,:,:)
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
      real,    allocatable :: flo(:,:,:,:),vso(:,:,:)
      real,    allocatable :: ff(:,:,:,:),gg(:,:,:,:)
      real,    allocatable :: us(:,:),uds(:,:),ulos(:,:)
      real,    allocatable :: dus(:,:),duds(:,:),dulos(:,:)
      real,    allocatable :: ddn(:,:),uulon(:,:),dulon(:,:)
      real,    allocatable :: uloulopn(:,:,:)
      real     :: uuilon(nlod,ntypd),duilon(nlod,ntypd)
      real     :: ulouilopn(nlod,nlod,ntypd)
!!! energy parameters
      real    :: ello(nlod,ntypd,max(2,jspd)),evac(2,max(2,jspd))
      real    :: epar(0:lmaxd,ntypd,max(2,jspd)),evdu(2,max(jspd,2))
!!! the Mmn matrices
      complex, allocatable :: mmnk(:,:,:,:),mmn(:,:,:)
      complex, allocatable :: amn(:,:,:),nablamat(:,:,:,:)
      complex, allocatable :: soctomom(:,:,:,:)
      complex, allocatable :: surfcurr(:,:,:,:)
      complex, allocatable :: socmmn(:,:,:)
      complex, allocatable :: a(:)
      complex, allocatable :: psiw(:,:,:)
      complex, allocatable :: anglmom(:,:,:,:)
      complex, allocatable :: orbcomp(:,:,:,:,:)
c..wf-hamiltonian in real space (hopping in the same unit cell)
      complex, allocatable :: hwfr(:,:),hwfr2(:,:)
c      real, allocatable :: ei(:)
      complex, allocatable :: work(:)
      real,allocatable::centers(:,:,:)
      logical :: l_file
      logical :: l_amn2, l_conjugate
      character(len=3) :: spin12(2)
      data   spin12/'WF1' , 'WF2'/
      character(len=30)  :: task
      integer,allocatable::irreduc(:)
      integer,allocatable::mapkoper(:)
      integer :: fullnkpts,kpt,kptibz,kptibz_b,j1,j2,j3,oper,oper_b,k
      real :: bkrot(3),dirfacs(3)
      integer :: ios,kplot,kplotoper,plotoper,gfcut
      complex :: phasust
      integer,allocatable::pair_to_do(:,:)
      integer :: ngopr1(natd)
      integer,allocatable::maptopair(:,:,:)
      integer :: wannierspin,jspin2
      real, allocatable :: rwork(:)
      real,allocatable::kdiff(:,:)
      integer,allocatable :: shiftkpt(:,:)
      integer :: unigrid(6),gfthick
      complex,allocatable::ujug(:,:,:,:),ujdg(:,:,:,:)
      complex,allocatable::djug(:,:,:,:),djdg(:,:,:,:)
      complex,allocatable::ujulog(:,:,:,:,:)
      complex,allocatable::djulog(:,:,:,:,:)
      complex,allocatable::ulojug(:,:,:,:,:)
      complex,allocatable::ulojdg(:,:,:,:,:)
      complex,allocatable::ulojulog(:,:,:,:,:,:)
      real(4) delta,delta1,time_interstitial,time_mmn
      real(4) time_total,delta2,delta3
      real(4) time_lapw_expand,time_rw,time_symm,time_film
      real(4) time_lapw_plot,time_ujugaunt,time_abcof
      integer :: n_start,n_end,mlotot,mlolotot,err
      integer :: mlot_d,mlolot_d,ilo,dir,length
      complex :: ci
      character(len=2) :: spin012(0:2)
      data spin012/'  ', '.1', '.2'/
      character(len=6) :: filename
      real :: arg,hescale
      complex :: nsfactor,nsfactor_b,value
      real :: b1(3),b2(3)
      real,parameter :: bohrtocm=0.529177e-8
      real,parameter :: condquant=7.7480917e-5
      integer :: npotmatfile,ig3,maxvac,irec,imz,ivac,ipot
      logical :: l_orbcompinp
      integer :: num_angl
      complex,allocatable :: vxy(:,:,:)


c-----initializations
      time_interstitial=0.0
      time_mmn=0.0
      time_total=0.0
      time_ujugaunt=0.0
      time_abcof=0.0
      time_rw=0.0
      time_symm=0.0
      time_film=0.0
      ngopr1(:)=1

      ci=cmplx(0.0,1.0)

      call cpu_time(time_total)

      l_p0 = .false.
      if (irank.eq.0) l_p0 = .true.

      sfp = 2* sqrt( pimach() )
      tpi = 2* pimach()
      lmd = lmaxd*(lmaxd+2)

!!!   should be changed in case the windows are really used
      nkpts = maxval(nkpt(:))

c-----read the input file to determine what to do
      wann%atomlist_num=natd
      wann%oc_num_orbs=natd
      call wann_read_inp(
     >         l_p0,
     <         wann)

c-----input file for orbital decomposition
c      if(wann%l_orbcomp.or.wann%l_orbcomprs)then
c         inquire(file='orbcomp_inp',exist=l_orbcompinp)
c         if(l_orbcompinp)then
c            open(159,file='orbcomp_inp')
c            read(159,*)wann%oc_num_orbs,wann%l_oc_f
c            allocate(wann%oc_orbs(wann%oc_num_orbs))
c            do n=1,wann%oc_num_orbs
c               read(159,*)wann%oc_orbs(n)
c            enddo
c            close(159)
c         else !default is all atoms including f
c            wann%oc_num_orbs=natd
c            wann%l_oc_f=.true.
c            allocate(wann%oc_orbs(wann%oc_num_orbs))
c            do n=1,wann%oc_num_orbs
c               wann%oc_orbs(n)=n
c            enddo          
c         endif
c      endif

      if(wann%l_updown)then
         call juDFT_error("TODO fix IO in Wannier")
         call wann_updown(
     >      mpi_comm,l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,0,l_soc, !irecl=0
     >      memd,lnonsph,clnu,lmplmd,mlh,nmem,llh,lo1l,
     >      theta,phi,soc_opt)
         CALL juDFT_end("updown done",irank)
      endif

      if(wann%l_byenergy.and.wann%l_byindex) CALL juDFT_error
     +     ("byenergy.and.byindex",calledby ="wannier")
      if(wann%l_byenergy.and.wann%l_bynumber) CALL juDFT_error
     +     ("byenergy.and.bynumber",calledby ="wannier")
      if(wann%l_bynumber.and.wann%l_byindex) CALL juDFT_error
     +     ("bynumber.and.byindex",calledby ="wannier")
      if(.not.(wann%l_bynumber.or.wann%l_byindex.or.wann%l_byenergy))
     &     CALL juDFT_error("no rule to sort bands",calledby ="wannier")

      efermi=ef
      if(.not.wann%l_fermi)efermi=0.0

#ifdef CPP_MPI
         call MPI_BARRIER(mpi_comm,ierr)
#endif

c**************************************************************
c   for bzsym=.true.: determine mapping between kpts and w90kpts
c**************************************************************
      if (wann%l_bzsym) then
         l_file=.false.
         inquire(file='w90kpts',exist=l_file)
         if(.not.l_file)  CALL juDFT_error
     +        ("w90kpts not found, needed if bzsym",calledby ="wannier")
         open(412,file='w90kpts',form='formatted')
         read(412,*)fullnkpts
         close(412)
         if(l_p0)print*,"fullnkpts=",fullnkpts
         if(fullnkpts<nkpts) CALL juDFT_error("fullnkpts.lt.nkpts"
     +        ,calledby ="wannier")
         allocate(irreduc(fullnkpts),mapkoper(fullnkpts))
         allocate(shiftkpt(3,fullnkpts))
         l_file=.false.
         inquire(file='kptsmap',exist=l_file)
         if(.not.l_file)  CALL juDFT_error
     +        ("kptsmap not found, needed if bzsym",calledby ="wannier")
         open(713,file='kptsmap')
         do i=1,fullnkpts
            read(713,*)kpt,irreduc(i),mapkoper(i),shiftkpt(:,i)
            if(kpt/=i) CALL juDFT_error("kpt.ne.i",calledby ="wannier")
            if(l_p0)print*,i,irreduc(i),mapkoper(i)
         enddo
         close(713)
         if(maxval(irreduc(:))/=nkpts) CALL juDFT_error
     +        ("max(irreduc(:))/=nkpts",calledby ="wannier")
      else
         fullnkpts=nkpts
      endif

      nrec = 0
      if(l_p0)then
      write (*,*) 'fermi energy:',efermi
      write (*,*) 'emin,emax=',e1s,e2s
      write (*,*) 'nbasfcn =',nbasfcn
      endif
      nlotot = 0
      do n = 1, ntype
        do l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        enddo
      enddo

      if((.not.wann%l_matrixmmn).and.(.not.wann%l_wann_plot).and.
     &   (.not.wann%l_matrixamn).and.(.not.wann%l_projmethod).and.
     &   (.not.wann%l_bestproj).and.(.not.wann%l_nabla).and.
     &   (.not.wann%l_mmn0).and.(.not.wann%l_surfcurr).and.
     &   (.not.wann%l_offdiposop).and.(.not.wann%l_anglmom).and.
     &   (.not.wann%l_orbcomp).and.(.not.wann%l_perturb) ) goto 1911

c**********************************************************
ccccccccccccccc   read in the bkpts file  ccccccccccccccccc
c**********************************************************
      if (wann%l_matrixmmn) then ! for Omega functional minimization
       l_bkpts = .false.
       inquire (file='bkpts',exist=l_bkpts)
       if (.not.l_bkpts)  CALL juDFT_error("need bkpts for matrixmmn"
     +      ,calledby ="wannier")
       open (202,file='bkpts',form='formatted',status='old')
       rewind (202)
       read (202,'(i4)') nntot
       if(l_p0)then
       write (*,*) 'nntot=',nntot
       write(*,*) 'fullnkpts=',fullnkpts
       write(*,*) 'nkpts=',nkpts
       endif
       allocate ( gb(1:3,1:nntot,1:fullnkpts),bpt(1:nntot,1:fullnkpts))
       do ikpt=1,fullnkpts
        do nn=1,nntot
         read (202,'(2i6,3x,3i4)')
     &     ikpt_help,bpt(nn,ikpt),(gb(i,nn,ikpt),i=1,3)
         if (ikpt/=ikpt_help)  CALL juDFT_error("ikpt.ne.ikpt_help"
     +        ,calledby ="wannier")
         if (bpt(nn,ikpt)>fullnkpts)  CALL juDFT_error("bpt.gt.fullnkpts"
     +        ,calledby ="wannier")
        enddo
       enddo
       close (202)
        allocate(kdiff(3,nntot))
      endif
c********************************************************
c      find symmetry-related elements in mmkb
c********************************************************
      if(wann%l_matrixmmn)then
         call cpu_time(delta)
         allocate(maptopair(3,fullnkpts,nntot))
         allocate(pair_to_do(fullnkpts,nntot))
         call wann_mmnk_symm(
     >     fullnkpts,nntot,bpt,gb,wann%l_bzsym,
     >     irreduc,mapkoper,l_p0,film,nop,invtab,mrot,odi%d1,
     >     tau,
     <     pair_to_do,maptopair,kdiff)
         call cpu_time(delta1)
         time_symm=delta1-delta
      endif
c*********************************************************
cccccccccccccccc   initialize the potential   cccccccccccc
c*********************************************************

      allocate ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
      allocate ( vz(nmzd,2,4),vrf(jmtd,0:nlhd,ntypd,jspd) )
      allocate ( vr(jmtd,ntypd,jspd) )
      allocate ( vso(jmtd,natd,2) )

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

      if(wann%l_soctomom)then
        call vsoc(
     >         jspins,ntype,jmtd,jri,
     >         dx,vr,rmsh,epar,.true.,
     <         vso)
      endif

      deallocate ( vpw,vzxy )
      deallocate(vrf)


      if(l_noco.and.film)then
         npotmatfile=25
         allocate(vpw(n3d,1))
         if(.not.odi%d1)
     &    allocate( vxy(nmzxyd,n2d-1,2) )

         OPEN (npotmatfile,FILE='potmat',FORM='unformatted',
     +                  STATUS='old')
         READ (npotmatfile) (vpw(ig3,1),ig3=1,n3d)
         READ (npotmatfile) (vpw(ig3,1),ig3=1,n3d)
         READ (npotmatfile) (vpw(ig3,1),ig3=1,n3d)
         maxvac=2
         if(odi%d1)maxvac=1
         DO ivac = 1,maxvac
c--->       if the two vacuua are equivalent, the potential file has to
c--->       be backspaced, because the potential is the same at both
c--->       surfaces of the film
            IF ((ivac.EQ.2) .AND. (nvac.EQ.1)) THEN
               DO irec = 1,4
                  BACKSPACE (npotmatfile)
               ENDDO
            ENDIF
c--->       load the non-warping part of the potential
            READ (npotmatfile)((vz(imz,ivac,ipot),imz=1,nmzd),ipot=1,4)

            if(.not.odi%d1)then
            DO ipot = 1,3
                  READ (npotmatfile)((vxy(imz,igvm2,ivac),
     +                           imz=1,nmzxy),igvm2=1,nq2-1)
            enddo    
            endif   
         enddo   
         CLOSE (npotmatfile)
         deallocate(vpw)
         if(.not.odi%d1)deallocate(vxy)
      endif   


cccccccccccccccc   end of the potential part  ccccccccccc
      wannierspin=jspd
      if(l_soc) wannierspin=2

      allocate ( kveclo(nlotot),nv(wannierspin) )
      allocate ( kveclo_b(nlotot),nv_b(wannierspin) )
      allocate ( k1(nvd,wannierspin),k2(nvd,wannierspin),
     &           k3(nvd,wannierspin) )
      allocate ( k1_b(nvd,wannierspin),k2_b(nvd,wannierspin),
     &           k3_b(nvd,wannierspin) )
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
c****************************************************
c cycle by spins starts! 
c****************************************************
      do 110 jspin=1,wannierspin   ! cycle by spins

c...read number of bands and wannier functions from file proj

c..reading the proj.1 / proj.2 / proj file
       l_proj=.false.  
       do j=jspin,0,-1
         inquire(file=trim('proj'//spin012(j)),exist=l_proj)
         if(l_proj)then
            filename='proj'//spin012(j)
            exit
         endif
       enddo

       if(l_proj)then
         open (203,file=trim(filename),status='old')
         rewind (203)
         read (203,*) nwfs,numbands
         rewind (203)
         close (203)
       elseif(wann%l_projmethod.or.wann%l_bestproj
     &                         .or.wann%l_matrixamn)then
         CALL juDFT_error("no proj/proj.1/proj.2"
     &             ,calledby ="wannier")
       endif  


       jspin2=jspin
       if(l_soc .and. jspins.eq.1)jspin2=1
       jsp_start = jspin ; jsp_end = jspin

cccccccccccc   read in the eigenvalues and vectors   cccccc

       call cdn_read0(
     >               lmaxd,ntypd,nlod,neigd,jspd,!wannierspin,
     >               irank,isize,jspin,jsp_start,jsp_end,
     >               l_noco,nrec,66,
     <               ello,evac,epar,bkpt,wk,n_bands,nrec1,n_size)

c.. now we want to define the maximum number of the bands by all kpts

      nbnd = 0

      i_rec = 0 ; n_rank = 0
c*************************************************************
c..writing down the eig.1 and/or eig.2 files
c*************************************************************
      if(l_p0)then
         CALL juDFT_error("TODO: implement correct IO for Wannier")
!         call wann_write_eig(
!     >            lmaxd,ntypd,nlod,neigd,nvd,wannierspin,
!     >            isize,jspin,nbasfcn,nlotot,
!     >            l_ss,l_noco,nrec,fullnkpts,irecl,
!     >            wann%l_bzsym,wann%l_byindex,wann%l_bynumber,
!     >            wann%l_byenergy,
!     >            irreduc,odi,wann%band_min(jspin),
!     >            wann%band_max(jspin),
!     >            numbands,
!     >            e1s,e2s,efermi,.false.,nkpts,
!     <            nbnd,kpoints)
         if(odi%d1)then
            kpoints(:)=kpoints(:)*bmat(3,3)
         endif
      endif!l_p0

! nbnd is calculated for process zero and is sent here to the others
#ifdef CPP_MPI
      if(l_p0)then
         do cpu_index=1,isize-1
      call MPI_SEND(nbnd,1,MPI_INTEGER,cpu_index,1,mpi_comm,ierr)
         enddo
      else
       call MPI_RECV(nbnd,1,MPI_INTEGER,0,1,mpi_comm,stt,ierr)
      endif
#endif

       print*,"process: ",irank," nbnd= ",nbnd
c##################################################################
      if(wann%l_mmn0)then
           allocate ( mmn(nbnd,nbnd,fullnkpts) )
           mmn(:,:,:) = cmplx(0.,0.)
           if((l_soc.or.l_noco) .and. (jspin.eq.1))
     &          allocate(socmmn(nbnd,nbnd,fullnkpts) )
      endif     
      if(wann%l_nabla)then
         allocate ( nablamat(3,nbnd,nbnd,fullnkpts) )
         nablamat = cmplx(0.,0.)
      endif

      if(wann%l_soctomom)then
         allocate ( soctomom(3,nbnd,nbnd,fullnkpts) )
         soctomom = cmplx(0.,0.)
      endif

      if(wann%l_surfcurr)then
         allocate ( surfcurr(3,nbnd,nbnd,fullnkpts) )
         surfcurr = cmplx(0.,0.)
      endif

      if(wann%l_anglmom)then
        if(.not.allocated(anglmom))then  
          allocate ( anglmom(3,nbnd,nbnd,fullnkpts) )
          anglmom=cmplx(0.,0.)
        endif  
      endif

      if(wann%l_orbcomp)then
         if(allocated(orbcomp))deallocate(orbcomp)
         if(wann%l_oc_f)then
            allocate(orbcomp(16,wann%oc_num_orbs,nbnd,nbnd,fullnkpts))
         else
            allocate(orbcomp(9,wann%oc_num_orbs,nbnd,nbnd,fullnkpts))
         endif
         orbcomp=cmplx(0.,0.)
      endif

      write (*,*) 'nwfs=',nwfs
      if(wann%l_projmethod.or.wann%l_bestproj.or.wann%l_matrixamn)then
        if(.not.allocated(amn))then
          allocate ( amn(nbnd,nwfs,fullnkpts) )
          amn(:,:,:) = cmplx(0.,0.)
        endif
      endif  

      if (wann%l_projmethod.or.wann%l_bestproj) then
        allocate ( psiw(nbnd,nwfs,fullnkpts) )
        psiw(:,:,:) = cmplx(0.,0.)
        if(.not.allocated(hwfr))then
           allocate ( hwfr(nwfs,nwfs) )
           hwfr(:,:) = cmplx(0.,0.)
        endif
      endif


      if (wann%l_matrixmmn) then
       if(.not.allocated(mmnk))then  
         allocate ( mmnk(nbnd,nbnd,nntot,fullnkpts) )
         mmnk = (0.,0.)
       endif  
      endif

      allocate ( flo(ntypd,jmtd,2,nlod) )

      na = 1
      do 40 n = 1,ntype
       do 30 l = 0,lmax(n)
c...compute the l-dependent, k-independent radial MT- basis functions
         call radfun(
     >              l,epar(l,n,jspin),vr(1,n,jspin2),jri(n),rmsh(1,n),
     >              dx(n),jmtd,
     <              ff(n,:,:,l),gg(n,:,:,l),us(l,n),
     <              dus(l,n),uds(l,n),duds(l,n),
     <              ddn(l,n),nodeu,noded,wronk)
   30  continue
c...and the local orbital radial functions
       do ilo = 1, nlo(n)
         call radflo(
     >             ntypd,nlod,jspd,jmtd,lmaxd,n,jspin,
     >             ello(1,1,jspin),vr(1,n,jspin2),
     >             jri(n),rmsh(1,n),dx(n),ff(n,1:,1:,0:),
     >             gg(n,1:,1:,0:),llo,nlo,l_dulo(1,n),irank,ulo_der,
     <             ulos(1,1),dulos(1,1),uulon(1,1),dulon(1,1),
     <             uloulopn(1,1,1),uuilon,duilon,ulouilopn,flo(n,:,:,:))
       enddo
c       na = na + neq(n)
   40 continue
c****************************************************************
c   calculate the k-independent uju*gaunt-matrix needed for
c   mmnmatrix
c****************************************************************
      if(wann%l_matrixmmn)then
         call cpu_time(delta)
         allocate(ujug(0:lmd,0:lmd,
     &              1:ntype,1:nntot))
         allocate(ujdg(0:lmd,0:lmd,
     &              1:ntype,1:nntot))
         allocate(djug(0:lmd,0:lmd,
     &              1:ntype,1:nntot))
         allocate(djdg(0:lmd,0:lmd,
     &              1:ntype,1:nntot))
         allocate(ujulog(0:lmd,1:nlod,-llod:llod,
     &              1:ntype,1:nntot))
         allocate(djulog(0:lmd,1:nlod,-llod:llod,
     &              1:ntype,1:nntot))
         allocate(ulojug(0:lmd,1:nlod,-llod:llod,
     &              1:ntype,1:nntot))
         allocate(ulojdg(0:lmd,1:nlod,-llod:llod,
     &              1:ntype,1:nntot))
         allocate(ulojulog(1:nlod,-llod:llod,1:nlod,-llod:llod,
     &              1:ntype,1:nntot))
         call wann_ujugaunt(
     >            llod,nntot,kdiff,lmax,ntype,ntypd,bbmat,
     >            bmat,nlod,nlo,llo,flo,ff,gg,jri,rmsh,dx,jmtd,
     >            lmaxd,lmd,
     <            ujug,ujdg,djug,djdg,
     <            ujulog,djulog,ulojug,ulojdg,ulojulog)
         call cpu_time(delta1)
         time_ujugaunt=delta1-delta
      endif


      i_rec = 0 ; n_rank = 0

c****************************************************************
c.. loop by kpoints starts!      each may be a separate task
c****************************************************************
      do 10 ikpt = wann%ikptstart,fullnkpts  ! loop by k-points starts
        kptibz=ikpt
        if(wann%l_bzsym) kptibz=irreduc(ikpt)
        if(wann%l_bzsym) oper=mapkoper(ikpt)

        i_rec = i_rec + 1
      if (mod(i_rec-1,isize).eq.irank) then

      allocate ( z(nbasfcn,neigd),we(neigd),
     &                   zz(nbasfcn,neigd),eigg(neigd) )

      call cpu_time(delta)
      n_start=1
      n_end=neigd
      call cdn_read(
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,!wannierspin,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,kptibz,66,neigd,
     >              n_start,n_end,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eigg,zz,cp_time)

      call cpu_time(delta1)
      time_rw=time_rw+delta1-delta
      nslibd = 0



c...we work only within the energy window

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      z(:,:) = cmplx(0.,0.)
#else
      z(:,:) = 0.
#endif
      eig(:) = 0.


      print*,"bands used:"
      do i = 1,nbands
        if ((eigg(i).ge.e1s.and.nslibd.lt.numbands.and.wann%l_bynumber)
     &    .or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.wann%l_byenergy)
     &    .or.(i.ge.wann%band_min(jspin).and.
     &    (i.le.wann%band_max(jspin))
     &       .and.wann%l_byindex))then

           print*,i
           nslibd = nslibd + 1
           eig(nslibd) = eigg(i)
           we(nslibd) = we(i)
           if(l_noco)then
             funbas=       nv(1)+nlotot
             funbas=funbas+nv(2)+nlotot
           else
             funbas=nv(jspin)+nlotot
           endif
           do  j = 1,funbas
             z(j,nslibd) = zz(j,i)
           enddo
        endif
      enddo

c***********************************************************
c              rotate the wavefunction
c***********************************************************
      if (wann%l_bzsym.and.oper.ne.1) then  !rotate bkpt
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
      else
         nsfactor=cmplx(1.0,0.0)
      endif
c      print*,"bkpt1=",bkpt

c******************************************************************

c...the overlap matrix Mmn which is computed for each k- and b-point

      noccbd = nslibd

      allocate ( acof(noccbd,0:lmd,natd),
     &           bcof(noccbd,0:lmd,natd),
     &           ccof(-llod:llod,noccbd,nlod,natd))

      acof(:,:,:) = cmplx(0.,0.) ; bcof(:,:,:) = cmplx(0.,0.)
      ccof(:,:,:,:) = cmplx(0.,0.)

c...generation the A,B,C coefficients in the spheres
c...for the lapws and local orbitals, summed by the basis functions

      call cpu_time(delta)
      call abcof(
     >           lmaxd,ntypd,neigd,noccbd,natd,nop,nvd,wannierspin,
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
      call cpu_time(delta1)
      time_abcof=time_abcof+delta1-delta

      call wann_abinv(
     >        ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
     >        noccbd,lmax,nlo,llo,invsat,invsatnr,bkpt,taual,
     X        acof,bcof,ccof)

      if(wann%l_orbcomp)then
         call wann_orbcomp(
     >            llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >            ntype,neq,nlo,llo,acof,bcof,ccof,
     >            ddn,uulon,dulon,uloulopn,
     >            wann%oc_num_orbs,
     >            wann%oc_orbs,
     >            wann%l_oc_f,
     =            orbcomp(:,:,:,:,ikpt))
      endif

      if(wann%l_anglmom)then
         call wann_anglmom(
     >                  llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >                  ntype,neq,nlo,llo,acof,bcof,ccof,
     >                  ddn,uulon,dulon,uloulopn,
     =                  anglmom(:,:,:,ikpt))
      endif

#ifdef CPP_TOPO
      if(wann%l_surfcurr)then
c         call wann_surfcurr_int(
c     >        nv2d,jspin,odi,ods,n3d,nmzxyd,n2d,ntypsd,
c     >        ntype,lmaxd,jmtd,ntypd,natd,nmzd,neq,nq3,nvac,
c     >        nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
c     >        invs,invs2,z1,delz,ngopr,ntypsy,jri,pos,zatom,
c     >        lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,ikpt,nnne,kk,
c     >        nvd,nlod,llod,nv(jspin),lmd,bkpt,omtil,nlo,llo,
c     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
c     >        vz(:,:,jspin2),
c     >        nslibd,nbasfcn,neigd,ff,gg,flo,acof,bcof,ccof,z,
c     >        surfcurr(:,:,:,ikpt))

         call wann_surfcurr_int2(
     >            nv2d,jspin,odi,ods,n3d,nmzxyd,n2d,ntypsd,
     >            ntype,lmaxd,jmtd,ntypd,natd,nmzd,neq,nq3,nvac,
     >            nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
     >            invs,invs2,z1,delz,ngopr,ntypsy,jri,
     >            pos,taual,zatom,rmt,
     >            lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,ikpt,
     >            nvd,nv(jspin),bkpt,omtil,
     >            k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >            nslibd,nbasfcn,neigd,z,
     <            dirfacs,
     >            surfcurr(:,:,:,ikpt))

         call wann_surfcurr(
     >            dirfacs,amat,
     >            jspin,ntype,lmaxd,lmax,natd,
     >            neq,noccbd,lmd,natd,llod,nlod,
     >            nlo,llo, 
     >            acof,bcof,ccof,
     >            us,dus,duds,uds,
     >            ulos(:,:),dulos(:,:),
     >            rmt,pos, 
     &            surfcurr(:,:,:,ikpt))
        write(6,*)"dirfacs=",dirfacs
      endif

      if(wann%l_soctomom)then
         call wann_soc_to_mom(
     >            jspin,ntype,lmaxd,lmax,natd,
     >            jmtd,jri,rmsh,dx,neq,
     >            noccbd,lmd,natd,llod,nlod,
     >            vso(:,:,1), 
     >            ff,gg,
     >            acof,bcof,ccof,
     &            soctomom(:,:,:,ikpt))
      endif

      if(wann%l_nabla)then
         call wann_nabla(
     >       nlo,llo,
     >       jspin,ntype,lmaxd,lmax,natd,
     >       jmtd,jri,rmsh,dx,neq,
     >       noccbd,lmd,natd,llod,nlod, 
     >       ff,gg,flo,
     >       acof,bcof,ccof,
     &       nablamat(:,:,:,ikpt))
         if(film.and..not.odi%d1)then
            call wann_nabla_vac(
     >          z1,nmzd,nv2d,k1d,k2d,k3d,n3d,nvac,
     >          ig,nmz,delz,ig2,area,bmat,
     >          bbmat,evac(:,jspin),bkpt,vz(:,:,jspin2),
     >          nslibd,jspin,k1,k2,k3,wannierspin,nvd,
     >          nbasfcn,neigd,z,nv,omtil,
     <          nablamat(:,:,:,ikpt))
         endif   
         addnoco=0
         do 41 i = n_rank+1,nv(jspin),n_size
           b1(1)=bkpt(1)+k1(i,jspin)
           b1(2)=bkpt(2)+k2(i,jspin)
           b1(3)=bkpt(3)+k3(i,jspin)
           b2(1)=b1(1)*bmat(1,1)+b1(2)*bmat(2,1)+b1(3)*bmat(3,1)
           b2(2)=b1(1)*bmat(1,2)+b1(2)*bmat(2,2)+b1(3)*bmat(3,2)
           b2(3)=b1(1)*bmat(1,3)+b1(2)*bmat(2,3)+b1(3)*bmat(3,3)
           do 42 j = n_rank+1,nv(jspin),n_size
c-->     determine index and phase factor
            i1 = k1(j,jspin) - k1(i,jspin)
            i2 = k2(j,jspin) - k2(i,jspin)
            i3 = k3(j,jspin) - k3(i,jspin)
            in = ig(i1,i2,i3)
            if (in.eq.0) goto 42
            phase   = rgphs(i1,i2,i3)
            phasust = cmplx(phase,0.0)*ustep(in)

            do m = 1,nslibd
             do n = 1,nslibd
              do dir=1,3  
                  
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
               value=phasust*z(i+addnoco,m)*conjg(z(j+addnoco,n))
               nablamat(dir,m,n,ikpt) =
     =         nablamat(dir,m,n,ikpt) -
     +          value*b2(dir)
#else
               value=phasust*cmplx(z(i+addnoco,m)*z(j+addnoco,n),0.0)
               nablamat(dir,m,n,ikpt) =
     =         nablamat(dir,m,n,ikpt) - 
     +          value*b2(dir)
#endif
              enddo 
             enddo
            enddo

 42        continue
 41      continue

      endif   
#endif
c      goto jump no longer needed?
c      if ((.not.wann%l_matrixmmn).and.(.not.wann%l_matrixamn).and.
c     &    (.not.wann%l_bestproj).and.(.not.wann%l_projmethod).and.
c     &    (.not.wann%l_mmn0)) goto 3

c------mmn0-matrix
       if(wann%l_mmn0)then
       addnoco=0
       if(l_noco.and.(jspin.eq.2))then
          addnoco=nv(1)+nlotot
       endif

c$$$       do 21 i = n_rank+1,nv(jspin),n_size
c$$$         do 22 j = n_rank+1,nv(jspin),n_size
c$$$
c$$$c-->     determine index and phase factor
c$$$            i1 = k1(j,jspin) - k1(i,jspin)
c$$$            i2 = k2(j,jspin) - k2(i,jspin)
c$$$            i3 = k3(j,jspin) - k3(i,jspin)
c$$$c            if(abs(i1).gt.k1d.or.abs(i2).gt.k2d
c$$$c     & .or. abs(i3).gt.k3d)print*,"interstitial warning"
c$$$            in = ig(i1,i2,i3)
c$$$            if (in.eq.0) goto 22
c$$$            phase   = rgphs(i1,i2,i3)
c$$$            phasust = cmplx(phase,0.0)*ustep(in)
c$$$            do m = 1,nslibd
c$$$             do n = 1,nslibd
c$$$#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
c$$$              mmn(m,n,ikpt) =
c$$$     =        mmn(m,n,ikpt) + 
c$$$     +          phasust*z(i+addnoco,m)*conjg(z(j+addnoco,n))
c$$$#else
c$$$              mmn(m,n,ikpt) =
c$$$     =        mmn(m,n,ikpt) + 
c$$$     +          phasust*cmplx(z(i+addnoco,m)*z(j+addnoco,n),0.0)
c$$$#endif
c$$$             enddo
c$$$            enddo
c$$$
c$$$  22     continue
c$$$  21   continue

c-----> interstitial contribution to mmn0-matrix
      call wann_mmkb_int(
     >         addnoco,addnoco,
     >         nvd,k1d,k2d,k3d,
     >         n3d,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),neigd,nbasfcn,z,nslibd,
     >         k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),z,nslibd,
     >         nbnd,
     >         rgphs,ustep,ig,(/ 0,0,0 /),
     <         mmn(:,:,ikpt))


c---> spherical contribution to mmn0-matrix

       call wann_mmk0_sph(
     >           llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >           ntype,neq,nlo,llo,acof(1:noccbd,:,:),
     >           bcof(1:noccbd,:,:),ccof(:,1:noccbd,:,:),
     >           ddn,uulon,dulon,uloulopn,
     =           mmn(:,:,ikpt))
c---> vacuum contribution to mmn0-matrix

       if (film .and. .not.odi%d1) then

         call wann_mmk0_vac(
     >           l_noco,nlotot,qss,
     >           z1,nmzd,nv2d,k1d,k2d,k3d,n3d,nvac,
     >           ig,nmz,delz,ig2,area,bmat,
     >           bbmat,evac(:,jspin),bkpt,vz(:,:,jspin2),
     >           nslibd,jspin,k1,k2,k3,wannierspin,nvd,
     >           nbasfcn,neigd,z,nv,omtil,
     <           mmn(:,:,ikpt))

       elseif (odi%d1) then

          call wann_mmk0_od_vac(
     >           l_noco,nlotot,
     >           z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >           ig,nmzxy,nmz,delz,ig2,odi%n2d,
     >           bbmat,evac(1,jspin),bkpt,odi%M,odi%mb,
     >           vz(:,1,jspin2),odi,
     >           nslibd,jspin,k1,k2,k3,wannierspin,nvd,area,
     >           nbasfcn,neigd,z,nv,sk2,phi2,omtil,qss,
     <           mmn(:,:,ikpt))

       endif
       endif !l_mmn0

c---> overlaps with the trial orbitals
      if (wann%l_projmethod.or.wann%l_bestproj.or.wann%l_matrixamn) then
         l_amn2=.false.
         call wann_amn (
     >            nslibd,nwfs,ntypd,nlod,llod,llo,nlo,
     >            lmaxd,jmtd,lmd,neq,natd,ikpt,nbnd,
     >            rmsh,rmt,jri,dx,lmax,
     >            us,dus,uds,duds,flo,
     >            ff,gg,acof,bcof,ccof,
     >            (l_soc.or.l_noco),jspin,
     &            l_amn2,
     <            amn(:,:,ikpt))
         if(l_amn2)then
          call wann_amn (
     >            nslibd,nwfs,ntypd,nlod,llod,llo,nlo,
     >            lmaxd,jmtd,lmd,neq,natd,ikpt,nbnd,
     >            rmsh,rmt,jri,dx,lmax,
     >            us,dus,uds,duds,flo,
     >            ff,gg,acof,bcof,ccof,
     >            (l_soc.or.l_noco),jspin,
     &            l_amn2,
     <            amn(:,:,ikpt),
     >            bkpt)
         endif
c         amn(ikpt,:,:)=amn(ikpt,:,:)*conjg(nsfactor)
      endif

c****************************************************************
c...         vanderbilt mmn matrix
c***************************************************************
      if (wann%l_matrixmmn) then   !  vanderbilt procedure Mmn matrix
      allocate ( z_b(nbasfcn,neigd), we_b(neigd) )


!!! the cycle by the nearest neighbors (nntot) for each kpoint

       do 15 ikpt_b = 1,nntot
          if(pair_to_do(ikpt,ikpt_b).eq.0)cycle !save time by symmetry
          kptibz_b=bpt(ikpt_b,ikpt)
          if(wann%l_bzsym) oper_b=mapkoper(kptibz_b)
          if (wann%l_bzsym) kptibz_b=irreduc(kptibz_b)



! now we need the wavefunctions for k_b kpoint

c        print*,"something to do"


          call cpu_time(delta)
          n_start=1
          n_end=neigd

        call cdn_read(
     >                lmaxd,ntypd,nlod,neigd,nvd,jspd,!wannierspin,
     >                irank,isize,kptibz_b,jspin,nbasfcn,nlotot,
     >                l_ss,l_noco,nrec,kptibz_b,66,
     >                neigd,n_start,n_end,
     <                nmat_b,nv_b,ello,evdu,epar,kveclo_b,
     <                k1_b,k2_b,k3_b,bkpt_b,wk_b,nbands_b,
     <                eigg,zz,cp_time)


        call cpu_time(delta1)
        time_rw=time_rw+delta1-delta
        nslibd_b = 0




#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      z_b(:,:) = cmplx(0.,0.)
#else
      z_b(:,:) = 0.
#endif
      eig_b(:) = 0.

        do i = 1,nbands_b
          if((eigg(i).ge.e1s.and.nslibd_b.lt.numbands
     & .and.wann%l_bynumber)
     &.or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.wann%l_byenergy)
     &.or.(i.ge.wann%band_min(jspin)
     &                 .and.
     &       (i.le.wann%band_max(jspin))
     &                 .and.
     &       wann%l_byindex))then
            nslibd_b = nslibd_b + 1
            eig_b(nslibd_b) = eigg(i)
            we_b(nslibd_b) = we_b(i)
            if(l_noco)then
               funbas =        nv_b(1) + nlotot
               funbas = funbas+nv_b(2) + nlotot
            else
               funbas = nv_b(jspin) + nlotot
            endif
            do  j = 1,funbas
              z_b(j,nslibd_b) = zz(j,i)
            enddo
          endif
        enddo

c***********************************************************
c              Rotate the wavefunction of next neighbor.
c***********************************************************
      if (wann%l_bzsym .and. (oper_b.ne.1)  ) then
         call wann_kptsrotate(
     >            natd,nlod,llod,
     >            ntypd,nlo,llo,invsat,
     >            l_noco,l_soc,
     >            ntype,neq,nlotot,
     >            kveclo_b,jspin,
     >            oper_b,nop,mrot,nvd,
     >            nv_b,
     >            shiftkpt(:,bpt(ikpt_b,ikpt)),
     >            tau,
     x            bkpt_b,k1_b(:,:),
     x            k2_b(:,:),k3_b(:,:),
     x            z_b,nsfactor_b)
      else
         nsfactor_b=cmplx(1.0,0.0)
      endif
c      print*,"kpt2=",bkpt_b

        noccbd_b = nslibd_b

ccccc   we start with the Mmn matrix   ccccccccccccc

!!! matrix elements of the interstitial overlap
!!! matrix with the weights given by products of
!!! the c-coeff. for different G-vectors, bands and k-points
!!! Mmn(k,b)(IR) = \sum(G,G')C_G^(k,n)*C_G'^(k+b,m)\theta_(G-G')

cccccccccccc  Spherical Contributions         ccccccccccccc

      allocate ( acof_b(noccbd_b,0:lmd,natd),
     &             bcof_b(noccbd_b,0:lmd,natd),
     &             ccof_b(-llod:llod,noccbd_b,nlod,natd))

!!! get the band-dependent k-dependent ab coeff.
      call cpu_time(delta)
      call abcof(
     >             lmaxd,ntypd,neigd,noccbd_b,natd,nop,nvd,wannierspin,
     >             lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >             ntype,mrot,ngopr1,taual,neq,lmax,rmt,omtil,
     >             bmat,bbmat,bkpt_b,
     >             k1_b,k2_b,k3_b,nv_b,nmat_b,noccbd_b,z_b,
     >             us(0,1),dus(0,1),uds(0,1),
     >             duds(0,1),ddn(0,1),invsat,invsatnr,
     >             ulos(1,1),uulon(1,1),dulon(1,1),
     >             dulos(1,1),llo,nlo,l_dulo,lapw_l,
     >             l_noco,l_ss,jspin,alph,beta,qss,kveclo_b,odi,ods,
     <             acof_b(1,0,1),bcof_b(1,0,1),
     <             ccof_b(-llod,1,1,1))
      call cpu_time(delta1)
      time_abcof=time_abcof+delta1-delta

      call wann_abinv(
     >         ntypd,natd,noccbd_b,lmaxd,lmd,llod,nlod,ntype,neq,
     >         noccbd_b,lmax,nlo,llo,invsat,invsatnr,bkpt_b,taual,
     X         acof_b,bcof_b,ccof_b)


ccccccccc  Interstitial  ccccccccccccccccccccccccc
!!! matrix elements of the interstitial overlap
!!! matrix with the weights given by products of
!!! the c-coeff. for different G-vectors, bands and k-points
!!! Mmn(k,b)(IR) = \sum(G,G')C_G^(k,n)*C_G'^(k+b,m)\theta_(G-G')
c... these overlaps are the same as in the mmk0 case, only differ
c... by the b-dependence of the C-coefficients
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call cpu_time(delta)

      addnoco=0
      addnoco2=0
      if(l_noco.and.(jspin.eq.2))then
         addnoco  = nv(1)   + nlotot
         addnoco2 = nv_b(1) + nlotot
      endif
      call wann_mmkb_int(
     >         addnoco,addnoco2,
     >         nvd,k1d,k2d,k3d,
     >         n3d,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),neigd,nbasfcn,z,nslibd,
     >         k1_b(:,jspin),k2_b(:,jspin),k3_b(:,jspin),
     >         nv_b(jspin),z_b,nslibd_b,
     >         nbnd,
     >         rgphs,ustep,ig,gb(:,ikpt_b,ikpt),
     <         mmnk(:,:,ikpt_b,ikpt))


        call cpu_time(delta1)
        time_interstitial=time_interstitial+delta1-delta

cccccccccccc  Spherical Contributions  ccccccccccccc
        call cpu_time(delta)
         call wann_mmkb_sph(
     >       nbnd,llod,nslibd,nslibd_b,nlod,natd,ntypd,lmd,jmtd,
     >       taual,nop,lmax,
     >       ntype,neq,nlo,llo,acof,bcof,ccof,bkpt_b,
     >       acof_b,bcof_b,ccof_b,gb(:,ikpt_b,ikpt),bkpt,ujug,ujdg,
     >       djug,djdg,ujulog,djulog,ulojug,ulojdg,ulojulog,kdiff,
     >       nntot,
     =       mmnk(:,:,ikpt_b,ikpt))

         call cpu_time(delta1)
         time_mmn=time_mmn+delta1-delta

c...vacuum contributions

         call cpu_time(delta)
        if (film .and. .not.odi%d1) then
           call wann_mmkb_vac(
     >              l_noco,nlotot,qss,
     >              nbnd,z1,nmzd,nv2d,k1d,k2d,k3d,n3d,nvac,
     >              ig,nmz,delz,ig2,area,bmat,
     >              bbmat,evac(:,jspin),bkpt,bkpt_b,vz(:,:,jspin2),
     >              nslibd,nslibd_b,jspin,k1,k2,k3,k1_b,k2_b,k3_b,
     >              wannierspin,nvd,
     >              nbasfcn,neigd,z,z_b,nv,nv_b,omtil,gb(:,ikpt_b,ikpt),
     <              mmnk(:,:,ikpt_b,ikpt))
         elseif (odi%d1) then

           call wann_mmkb_od_vac(
     >           l_noco,nlotot,
     >           nbnd,z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >           ig,nmzxy,nmz,delz,ig2,odi%n2d,
     >           bbmat,evac(1,jspin),bkpt,bkpt_b,odi%M,odi%mb,
     >           vz(:,1,jspin2),odi,
     >           nslibd,nslibd_b,jspin,k1,k2,k3,k1_b,k2_b,k3_b,
     >           wannierspin,nvd,area,nbasfcn,neigd,
     >           z,z_b,nv,nv_b,sk2,phi2,omtil,
     >           gb(:,ikpt_b,ikpt),qss,
     <           mmnk(:,:,ikpt_b,ikpt))
        endif
c        mmnk(:,:,ikpt_b,ikpt)=
c     &          mmnk(:,:,ikpt_b,ikpt)*nsfactor*conjg(nsfactor_b)
        call cpu_time(delta1)
        time_film=time_film+delta1-delta

        deallocate ( acof_b,bcof_b,ccof_b )

15      continue ! end of loop by the nearest k-neighbors

        deallocate ( z_b,we_b )

        endif    !   if wann%l_matrixmmn = true

c$$$  3     continue !goto jump no longer needed?

        if(.not.wann%l_bzsym)oper=0
        if(.not.wann%l_plot_symm.or.oper.eq.1)then
        if (wann%l_wann_plot) then
         if (slice) then
          if (ikpt.eq.kk) then

           write (6,*) 'nnne=',nnne
           write (6,*) 'eig(nnne)=',eig(nnne)
           write (6,*) 'we(nnne)=',we(nnne)

           call wann_plot(nv2d,jspin,odi,ods,n3d,nmzxyd,n2d,ntypsd,
     >        ntype,lmaxd,jmtd,ntypd,natd,nmzd,neq,nq3,nvac,
     >        nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
     >        invs,invs2,z1,delz,ngopr,ntypsy,jri,pos,zatom,
     >        lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,ikpt,nnne,kk,
     >        nvd,nlod,llod,nv(jspin),lmd,bkpt,omtil,nlo,llo,
     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
     >        vz(:,:,jspin2),
     >        nslibd,nbasfcn,neigd,ff,gg,flo,acof,bcof,ccof,z)

          endif
         else ! not slice

          call wann_plot(nv2d,jspin,odi,ods,n3d,nmzxyd,n2d,ntypsd,
     >        ntype,lmaxd,jmtd,ntypd,natd,nmzd,neq,nq3,nvac,
     >        nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
     >        invs,invs2,z1,delz,ngopr,ntypsy,jri,pos,zatom,
     >        lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,ikpt,nnne,kk,
     >        nvd,nlod,llod,nv(jspin),lmd,bkpt,omtil,nlo,llo,
     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
     >         vz(:,:,jspin2),
     >        nslibd,nbasfcn,neigd,ff,gg,flo,acof,bcof,ccof,z)

          if(wann%l_plot_symm.and.wann%l_bzsym)then
           do kplot=1,fullnkpts
            if(irreduc(kplot).eq.kptibz)then
              plotoper=mapkoper(kplot) 
              if(plotoper.lt.0)then
                 plotoper=-plotoper
                 l_conjugate=.true.
              else
                 l_conjugate=.false.
              endif
              kplotoper=invtab(plotoper)
              call wann_plot_symm(jspin,mrot(:,:,kplotoper),ikpt,kplot,
     &                                 l_conjugate)
            endif
           enddo
          endif
                
         endif
        endif
        endif !wann%l_plot_symm
      deallocate ( acof,bcof,ccof,z,we,zz,eigg )

      write (*,*) 'nslibd=',nslibd

      if(wann%l_projmethod.or.wann%l_bestproj)then
         call wann_projmethod(
     >            fullnkpts,
     >            wann%l_projmethod,wann%l_bestproj,
     >            ikpt,nwfs,nslibd,amn,eig,
     <            psiw,hwfr)
      endif ! projmethod

      endif   ! loop by processors

10    continue ! end of cycle by the k-points

      if(wann%l_matrixmmn)deallocate(ujug,ujdg,djug,djdg,
     &                 ujulog,djulog,ulojulog,ulojug,ulojdg)

#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,ierr)
#endif
c******************************************************
c     Write down the projections.
c******************************************************
  5   continue

      if(wann%l_nabla)then
        hescale=tpi*condquant/bohrtocm/omtil
        hescale=sqrt(hescale)
        nablamat=nablamat*hescale
        call wann_write_nabla(
     >           mpi_comm,l_p0,spin12(jspin)//'.nabl',
     >           'Matrix elements of nabla operator',
     >           nbnd,fullnkpts,nbnd,
     >           irank,isize,
     <           nablamat)
      endif

      if(wann%l_soctomom)then
        hescale=tpi*condquant/bohrtocm/omtil
        hescale=sqrt(hescale)
        soctomom=soctomom*hescale
        call wann_write_nabla(
     >           mpi_comm,l_p0,spin12(jspin)//'.stm',
     >           'Matrix elements of stm operator',
     >           nbnd,fullnkpts,nbnd,
     >           irank,isize,
     <           soctomom)
      endif

      if(wann%l_surfcurr)then
        surfcurr = conjg( surfcurr/cmplx(0.0,2.0) )
        hescale=tpi*condquant/bohrtocm/omtil
        hescale=sqrt(hescale)
        surfcurr=surfcurr*hescale
        call wann_write_nabla(
     >           mpi_comm,l_p0,spin12(jspin)//'.surfcurr',
     >           'Surface currents',
     >           nbnd,fullnkpts,nbnd,
     >           irank,isize,
     <           surfcurr)
      endif

      if((l_soc.or.l_noco).and.wann%l_mmn0)then
        call wann_write_amn(
     >           mpi_comm,l_p0,spin12(jspin)//'.socmmn0',
     >           'Overlaps of the wavefunct. at the same kpoint',
     >           nbnd,fullnkpts,nbnd,
     >           irank,isize,.false.,wann%l_unformatted,
     <           mmn)
      endif !l_soc and l_mmn0  

      if(wann%l_orbcomp)then
         num_angl=9
         if(wann%l_oc_f)num_angl=16
         call wann_write_matrix5(
     >            mpi_comm,l_p0,spin12(jspin)//'.orbcomp',
     >            'angular components',
     >            nbnd,nbnd,
     >            num_angl,wann%oc_num_orbs,fullnkpts,
     >            irank,isize,
     <            orbcomp)
      endif

      if((l_soc.or.l_noco) .and. (jspin.eq.1)) then
         if(wann%l_mmn0) socmmn(:,:,:)=mmn(:,:,:)
         goto 912
      endif

      if(l_soc.or.l_noco ) then
         jspin2=1
         if(wann%l_mmn0)      mmn(:,:,:)=socmmn(:,:,:)+mmn(:,:,:)
         if(wann%l_mmn0)     deallocate(socmmn)
      endif

      if (wann%l_matrixamn)then
         call wann_write_amn(
     >           mpi_comm,l_p0,spin12(jspin2)//'.amn',
     >           'Overlaps of the wavefunct. with the trial orbitals',
     >           nbnd,fullnkpts,nwfs,
     >           irank,isize,.false.,wann%l_unformatted,
     <           amn)
      endif !wann%l_matrixamn

      if(wann%l_anglmom)then
        call wann_write_matrix4(
     >              mpi_comm,l_p0,spin12(jspin2)//'.anglmom',
     >              'Matrix elements of angular momentum',
     >              nbnd,nbnd,3,fullnkpts,
     >              irank,isize,
     <              anglmom)
      endif

      if (l_proj) then
c**************************************************************
c            for projmethod: write down WF1.umn
c*************************************************************
       if((wann%l_projmethod.or.wann%l_bestproj))then
         call wann_write_amn(
     >           mpi_comm,l_p0,spin12(jspin2)//'.umn',
     >           'transformation to first guess Wannier functions',
     >           nbnd,fullnkpts,nwfs,
     >           irank,isize,.false.,wann%l_unformatted,
     <           psiw)
#ifdef CPP_MPI
         allocate( hwfr2(size(hwfr,1),size(hwfr,2)) )
         length=nwfs*nwfs
         CALL MPI_REDUCE(
     +           hwfr,hwfr2,length,
     +           CPP_MPI_COMPLEX,MPI_SUM,0,
     +           mpi_comm,ierr)
         hwfr=hwfr2
         deallocate(hwfr2)
#endif
c********************************************************
c        projmethod: hamiltonian matrix in real space
c********************************************************
        if(l_p0)then
        write (6,*) 'the hamiltonian matrix in real space:'
        do i = 1,nwfs
         do j = 1,nwfs
           write (6,*) '   WFs:',i,'and',j
           write (6,*) '     matrix element:',hwfr(i,j)
         enddo
        enddo
        endif !l_p0
        deallocate(hwfr)
       endif !wann%l_projmethod or wann%l_bestproj
      endif !l_proj

c*********************************************************
c.....write down the mmn0 matrix 
c*********************************************************
      if(wann%l_mmn0)then
         call wann_write_amn(
     >            mpi_comm,l_p0,spin12(jspin2)//'.mmn0',
     >            'Overlaps of the wavefunct. at the same kpoint',
     >            nbnd,fullnkpts,nbnd,
     >            irank,isize,.false.,wann%l_unformatted,
     <            mmn)
      endif !wann%l_mmn0   
      
c*****************************************************
c.....write down the matrix M^{k,b}_{mn} 
c*****************************************************
      if(wann%l_matrixmmn)then
         call wann_write_mmnk(
     >            mpi_comm,jspin2,l_p0,fullnkpts,nntot,wann,
     >            maptopair,pair_to_do,nbnd,bpt,gb,
     >            isize,irank,
     <            mmnk)
      endif !wann%l_matrixmmn

 912  continue

      if( allocated (nablamat) ) deallocate( nablamat )
      if( allocated (soctomom) ) deallocate( soctomom )
      if( allocated (surfcurr) ) deallocate( surfcurr )
      if( allocated( mmn ) ) deallocate(mmn)
      if( allocated( amn ) ) then
         if(.not.(l_soc.or.l_noco))deallocate(amn)
      endif
      if ( allocated (psiw) ) deallocate ( psiw )
      if (wann%l_matrixmmn) then
         if(.not.(l_soc.or.l_noco))deallocate (mmnk)
      endif   
      if (wann%l_anglmom)then
         if(.not.(l_soc.or.l_noco))deallocate (anglmom)
      endif
      deallocate(flo)

      if(.not.l_noco)nrec=nrec+nkpts

c#ifdef CPP_MPI
c      call MPI_BARRIER(mpi_comm,ierr)
c#endif

110   continue ! end of cycle by spins


      if (wann%l_matrixmmn) deallocate (gb,bpt)
      if(wann%l_matrixmmn)deallocate(pair_to_do,maptopair)

      deallocate ( vr,vz,kveclo,nv,k1,k2,k3 )
      deallocate ( ff,gg,us,dus,duds,uds,ddn)
      deallocate ( ulos,dulos,uulon,dulon,uloulopn )
      if (wann%l_matrixmmn) deallocate ( kveclo_b,nv_b,k1_b,k2_b,k3_b )
      if (wann%l_bzsym)  deallocate(irreduc,mapkoper)

 1911 continue

      call cpu_time(delta)
      time_total=delta-time_total
      write(6,*)"time_total=",time_total
      write(6,*)"time_mmn=",time_mmn
      write(6,*)"time_interstitial=",time_interstitial
      write(6,*)"time_ujugaunt=",time_ujugaunt
      write(6,*)"time_abcof=",time_abcof
      write(6,*)"time_symm=",time_symm
      write(6,*)"time_rw=",time_rw
      write(6,*)"time_film=",time_film

#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,ierr)
#endif
      CALL juDFT_error("TODO: fix IO in Wannier")
      call wann_postproc(
     >         mpi_comm,wann,l_p0,film,jspins,n2d,
     >         natd,pos,amat,bmat,ntype,neq,zatom,
     >         omtil,l_soc,l_noco,neigd,fullnkpts,
     >         odi,l_proj,
     >         l_dulo,l_ss,lmaxd,ntypd,nop,nvd,jspd,nbasfcn,
     >         llod,nlod,nwdd,nlo,llo,lapw_l,invtab,mrot,ngopr,
     >         lmax,invsat,invsatnr,nkpt,taual,rmt,bbmat,alph,
     >         beta,qss,sk2,phi2,ods,
     >         irank,isize,n3d,nmzxyd,nmzd,jmtd,nlhd,nq3,nvac,
     >         invs,invs2,nlh,jri,ntypsy,ntypsd,nkptd,
     >         dx,rmsh,
     >         e1s,e2s,ulo_der,ustep,k1d,k2d,k3d,nop2,
     >         ig,rgphs,slice,kk,nnne,0,nv2d,nmzxy,nmz,nq2,             !irecl=0 here
     >         delz,area,z1,volint,ig2,tau,symor,ef)

#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,ierr)
#endif
      if(.not.wann%l_ldauwan)  CALL juDFT_end("wannier good",irank)

      END SUBROUTINE wannier
      END MODULE m_wannier
