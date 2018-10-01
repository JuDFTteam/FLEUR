!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_updown
      use m_juDFT
      CONTAINS
      SUBROUTINE wann_updown(
     >      mpi,input,sym,atoms,stars,vacuum,sphhar,oneD,noco,cell,vTot,
     >      enpara,
     >      eig_id,l_real,mpi_comm,l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,l_soc,
     >      memd,lnonsph,clnu,lmplmd,mlh,nmem,llh,lo1l,
     >      theta,phi)
c****************************************************************************
c     Calculate spin-up-spin-down matrix elements.
c     Useful in the case of spin-orbit coupling and noncollinear magnetism.
c     Implemented matrix elements:
c         (i)   Overlaps => Pauli matrix elements.
c         (ii)  Momentum operator => Spin current.
c         (iii) H_{SOC} => Include SOC in the basis set of Wannier functions.
c         (iv)  Spin-current into the MT-spheres => Torque.
c         (v)   Matrix elements of exchange field => Torque.
c     Frank Freimuth
c****************************************************************************
      use m_types
      use m_juDFT
      use m_wann_rw_eig
      use m_abcof
      use m_radfun
      use m_radflo
      use m_cdnread, only : cdn_read0, cdn_read
      use m_constants, only : pimach
!      use m_wann_mmk0_od_vac
      use m_wann_mmk0_vac
      use m_wann_mmk0_updown_sph
      use m_wann_abinv
      use m_wann_kptsrotate
      use m_wann_read_inp
      use m_wann_radovlp_integrals
      use m_matmul,only : matmul3,matmul3r
      use m_wann_socmat
      use m_wann_write_hsomtx
      use m_wann_write_nabla
      use m_wann_write_amn
      use m_wann_write_nabla
      use m_wann_mmkb_int
#ifdef CPP_TOPO
      use m_wann_perpmag
      use m_wann_nabla_updown
      use m_wann_surfcurr_updown
      use m_wann_surfcurr_int_updown
      use m_wann_gabeffgagaunt
#endif
      IMPLICIT NONE
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer ierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)

#endif

      TYPE(t_mpi),INTENT(IN)    :: mpi
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_oneD),INTENT(IN)   :: oneD
      TYPE(t_noco),INTENT(IN)   :: noco
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_potden),INTENT(IN) :: vTot
      TYPE(t_enpara),INTENT(IN) :: enpara

      logical, intent (in) :: invs,invs2,film,slice,symor,l_real
      integer, intent (in) :: lmaxd,ntypd,neigd,nkptd,kk,nnne,mpi_comm
      integer, intent (in) :: natd,nop,nvd,jspd,nbasfcn,nq2,nop2,eig_id
      integer, intent (in) :: llod,nlod,ntype,n3d,n2d
      integer, intent (in) :: nmzxyd,nmzd,jmtd,nlhd,nq3,nvac
      integer, intent (in) :: ntypsd,jspins,k1d,k2d,k3d
      real,    intent (in) :: omtil,e1s,e2s,delz,area,z1,volint
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      complex, intent (in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: nlh(ntypsd),jri(ntypd),ntypsy(natd)
      integer, intent (in) :: nlo(ntypd),llo(nlod,ntypd),lapw_l(ntypd)
      integer, intent (in) :: invtab(nop),mrot(3,3,nop),ngopr(natd)
      integer, intent (in) :: neq(ntypd),lmax(ntypd)
      integer, intent (in) :: invsat(natd),invsatnr(natd),nkpt
      integer, intent (in) :: irank,isize,nv2d,nmzxy,nmz
      integer, intent (in) :: ulo_der(nlod,ntypd),ig2(n3d)
      real,    intent (in) :: taual(3,natd),rmt(ntypd),dx(ntypd)
      real,    intent (in) :: amat(3,3),bmat(3,3),bbmat(3,3)
      real,    intent (in) :: rmsh(jmtd,ntypd),tau(3,nop),zatom(ntype)
      real,    intent (in) :: alph(ntypd),beta(ntypd),qss(3)
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
     
cccccccccccccccccc   local variables   cccccccccccccccccccc
      type(t_wann) :: wann
      integer :: lmd,nlotot,n,nmat,iter,ikpt,ikpt_b,nmat_b
      integer :: addnoco,addnoco2,funbas,loplod
      integer :: noccbd,noccbd_b,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
      integer :: jsp_start,jsp_end,nrec,nrec1,nbands,nbands_b
      integer :: nodeu,noded,n_size,na,n_rank,nbnd,numbands
      integer :: i1,i2,i3,in,lda
      integer :: n_bands(0:neigd),nslibd,nslibd_b
      character(len=8) :: dop,iop,name(10)
      real    :: bkpt(3),bkpt_b(3),sfp,tpi,wronk,wk,wk_b,phase
      complex :: c_phase, zzConjgTemp
      real    :: eig(neigd),cp_time(9),efermi
      logical :: l_p0,l_bkpts,l_proj,l_amn,l_mmn
!!! energy window boundaries
      integer, allocatable :: nv(:)
      integer, allocatable :: nv_b(:)
      integer, allocatable :: k1(:,:),k2(:,:),k3(:,:)
      integer, allocatable :: k1_b(:,:),k2_b(:,:),k3_b(:,:)

      real,    allocatable :: eigg(:,:)
      real kpoints(nkptd)
!!! a and b coeff. constructed for each k-point
      complex, allocatable :: acof(:,:,:,:)
      complex, allocatable :: bcof(:,:,:,:)
      complex, allocatable :: ccof(:,:,:,:,:)
!!! the parameters for the number of wfs
      integer :: nwfs
!!! the potential in the spheres and the vacuum
      real, allocatable :: vr(:,:,:),vz(:,:,:)
!!! bkpts data
      integer nntot,ikpt_help
      integer, allocatable :: gb(:,:,:),bpt(:,:)
!!! radial wavefunctions in the muffin-tins and more ...
      real,    allocatable :: flo(:,:,:,:,:)
      real,    allocatable :: ff(:,:,:,:,:),gg(:,:,:,:,:)

      real     :: uuilon(nlod,ntypd),duilon(nlod,ntypd)
      real     :: ulouilopn(nlod,nlod,ntypd)
!!! energy parameters
      real    :: ello(nlod,ntypd,max(2,jspd)),evac(2,max(2,jspd))
      real    :: epar(0:lmaxd,ntypd,max(2,jspd)),evdu(2,max(jspd,2))
!!! the Mmn matrices
      complex, allocatable :: surfcurr(:,:,:,:)
      complex, allocatable :: hsomtx(:,:,:,:,:)
      complex, allocatable :: mmnk(:,:,:,:),mmn(:,:,:)
      complex, allocatable :: perpmag(:,:,:)
      complex, allocatable :: amn(:,:,:),nablamat(:,:,:,:)
      complex, allocatable :: socamn(:,:,:)
      complex, allocatable :: socmmn(:,:,:)
      complex, allocatable :: socmmnk(:,:,:,:)
      complex, allocatable :: omn(:,:,:),vec(:,:),a(:)
      complex, allocatable :: smn(:,:,:)
      complex, allocatable :: psiw(:,:,:)
c..wf-hamiltonian in reciprocal space Hwf(k)
      complex, allocatable :: hwfk(:,:,:)
c..egenvalues of the wf-hamiltonian Hwf(k)
      real, allocatable :: eigw(:,:)
c..wf-hamiltonian in real space (hopping in the same unit cell)
      complex, allocatable :: hwfr(:,:)
      real, allocatable :: ei(:)
      complex, allocatable :: work(:)
      real,allocatable::centers(:,:,:)
      integer, allocatable :: iwork(:)
      integer :: info,lwork,lrwork,liwork,lee
      character :: jobz, uplo
      logical :: l_file
      logical :: l_amn2, l_conjugate
      character(len=3) :: spin12(2)
      data   spin12/'WF1' , 'WF2'/
      character(len=30)  :: task
      integer,allocatable::irreduc(:)
      integer,allocatable::mapkoper(:)
      integer :: fullnkpts,kpt,kptibz,kptibz_b,j1,j2,j3,oper,oper_b,k
      real :: bkrot(3)
      integer :: ios,kplot,kplotoper,plotoper,gfcut
      integer :: funbas_start,funbas_stop
      complex :: phasust
      integer,allocatable::pair_to_do(:,:)
      integer :: ngopr1(natd)
      integer,allocatable::maptopair(:,:,:)
      integer :: ldu,ldv,wannierspin,jspin2
      complex, allocatable :: amn_copy(:,:)
      complex, allocatable :: sunit(:,:)
      real,    allocatable :: sdiag(:)
      complex, allocatable :: umn(:,:),vmn(:,:)
      real, allocatable :: rwork(:)
      character :: jobu*1,jobv*1
      real,allocatable::kdiff(:,:)
      integer,allocatable :: shiftkpt(:,:)
      integer :: unigrid(6),gfthick
      complex,allocatable::ujug(:,:,:,:),ujdg(:,:,:,:)
      complex,allocatable::djug(:,:,:,:),djdg(:,:,:,:)
      complex,allocatable::ujulog(:,:,:,:)
      complex,allocatable::djulog(:,:,:,:)
      complex,allocatable::ulojug(:,:,:,:)
      complex,allocatable::ulojdg(:,:,:,:)
      complex,allocatable::ulojulog(:,:,:,:)
      real(4) delta,delta1,time_interstitial,time_mmn
      real(4) time_total,delta2,delta3
      real(4) time_lapw_expand,time_rw,time_symm,time_film
      real(4) time_lapw_plot,time_ujugaunt,time_abcof
      integer :: n_start,n_end,mlotot,mlolotot,err
      integer :: mlot_d,mlolot_d,ilo,dir,l1,l2,ii,jj
      character(len=2)    :: spin012(0:2)
      data spin012/'  ', '.1', '.2'/
      character(len=6)    :: filename
      real                :: arg,hescale
      complex             :: nsfactor,nsfactor_b,value
      real                :: b1(3),b2(3),b3(3)
      real,allocatable    :: radial1_ff(:,:,:,:,:)
      real,allocatable    :: radial1_gg(:,:,:,:,:)
      real,allocatable    :: radial1_fg(:,:,:,:,:)
      real,allocatable    :: radial1_gf(:,:,:,:,:)
      real,allocatable    :: radial1_flo(:,:,:,:,:)
      real,allocatable    :: radial1_glo(:,:,:,:,:)
      real,allocatable    :: radial1_lof(:,:,:,:,:)
      real,allocatable    :: radial1_log(:,:,:,:,:)
      real,allocatable    :: radial1_lolo(:,:,:,:,:)

      real,parameter      :: bohrtocm=0.529177e-8
      real,parameter      :: condquant=7.7480917e-5
      real                :: dirfacs(3)
      complex,allocatable :: ubeffug(:,:,:) 
      complex,allocatable :: ubeffdg(:,:,:) 
      complex,allocatable :: dbeffug(:,:,:) 
      complex,allocatable :: dbeffdg(:,:,:) 
      complex,allocatable :: ubeffulog(:,:,:) 
      complex,allocatable :: dbeffulog(:,:,:) 
      complex,allocatable :: ulobeffug(:,:,:) 
      complex,allocatable :: ulobeffdg(:,:,:) 
      complex,allocatable :: ulobeffulog(:,:,:) 
      integer :: spinloop1,spinloop2,ilop

      TYPE(t_usdus) :: usdus
      TYPE(t_mat)   :: zMat(2),zzMat(2)
      TYPE(t_lapw)  :: lapw

      IF(l_ss) CALL juDFT_error("spin-spiral not yet in this version"
     +     ,calledby ="wann_updown")

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

      call cpu_time(time_total)

      l_p0 = .false.
      if (irank.eq.0) l_p0 = .true.

      sfp = 2* sqrt( pimach() )
      tpi = 2* pimach()
      lmd = lmaxd*(lmaxd+2)

!!!   should be changed in case the windows are really used
      nkpts = nkpt

c-----read the input file to determine what to do
      call wann_read_inp(input,l_p0,wann)

      IF(wann%l_byenergy.AND.wann%l_byindex) CALL juDFT_error
     +     ("byenergy.and.byindex",calledby ="wann_updown")
      IF(wann%l_byenergy.AND.wann%l_bynumber) CALL juDFT_error
     +     ("byenergy.and.bynumber",calledby ="wann_updown")
      IF(wann%l_bynumber.AND.wann%l_byindex) CALL juDFT_error
     +     ("bynumber.and.byindex",calledby ="wann_updown")
      IF(.NOT.(wann%l_bynumber.OR.wann%l_byindex.OR.wann%l_byenergy))
     &     CALL juDFT_error("no rule to sort bands",calledby
     +     ="wann_updown")

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
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("w90kpts not found, needed if bzsym",calledby
     +        ="wann_updown")
         open(412,file='w90kpts',form='formatted')
         read(412,*)fullnkpts
         close(412)
         if(l_p0)print*,"fullnkpts=",fullnkpts
         IF(fullnkpts<nkpts) CALL juDFT_error("fullnkpts.lt.nkpts"
     +        ,calledby ="wann_updown")
         allocate(irreduc(fullnkpts),mapkoper(fullnkpts))
         allocate(shiftkpt(3,fullnkpts))
         l_file=.false.
         inquire(file='kptsmap',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("kptsmap not found, needed if bzsym",calledby
     +        ="wann_updown")
         open(713,file='kptsmap')
         do i=1,fullnkpts
            read(713,*)kpt,irreduc(i),mapkoper(i),shiftkpt(:,i)
            IF(kpt/=i) CALL juDFT_error("kpt.ne.i",
     +                                  calledby ="wann_updown")
            if(l_p0)print*,i,irreduc(i),mapkoper(i)
         enddo
         close(713)
         IF(MAXVAL(irreduc(:))/=nkpts) CALL juDFT_error
     +        ("max(irreduc(:))/=nkpts",calledby ="wann_updown")
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

c*********************************************************
cccccccccccccccc   initialize the potential   cccccccccccc
c*********************************************************

      allocate (vz(nmzd,2,jspd))
      allocate (vr(jmtd,ntypd,jspd))

      vz = vTot%vacz

      do jspin = 1,jspins
        do n = 1, ntype
          do j = 1,jri(n)
            vr(j,n,jspin) = vTot%mt(j,0,n,jspin)
          enddo
        enddo
      enddo
    
cccccccccccccccc   end of the potential part  ccccccccccc
      wannierspin=jspd
      if(l_soc) wannierspin=2

      allocate ( nv(wannierspin) )
      allocate ( k1(nvd,wannierspin),k2(nvd,wannierspin),
     &           k3(nvd,wannierspin) )
      allocate ( ff(ntypd,jmtd,2,0:lmaxd,2) )
      allocate ( gg(ntypd,jmtd,2,0:lmaxd,2) )
      allocate ( usdus%us(0:lmaxd,ntypd,2) )
      allocate ( usdus%uds(0:lmaxd,ntypd,2) )
      allocate ( usdus%dus(0:lmaxd,ntypd,2) )
      allocate ( usdus%duds(0:lmaxd,ntypd,2) )
      allocate ( usdus%ddn(0:lmaxd,ntypd,2) )
      allocate ( usdus%ulos(nlod,ntypd,2) )
      allocate ( usdus%dulos(nlod,ntypd,2) )
      allocate ( usdus%uulon(nlod,ntypd,2) )
      allocate ( usdus%dulon(nlod,ntypd,2) )
      allocate ( usdus%uloulopn(nlod,nlod,ntypd,2) )
c****************************************************
c     jspin=1 
c****************************************************
      jspin=1

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
       else
          CALL juDFT_error("no proj/proj.1/proj.2",
     +                     calledby ="wann_updown")
       endif  


       jspin2=jspin
       if(l_soc .and. jspins.eq.1)jspin2=1
       jsp_start = jspin ; jsp_end = jspin

cccccccccccc   read in the eigenvalues and vectors   cccccc


       do jspin=1,jspins
          CALL judft_error("TODO:wann_updown")
!         call cdn_read0(eig_id,irank,isize,jspin,jspd,l_noco,
 !    <                  ello,evac,epar,wk,n_bands,n_size)
       enddo

       jspin=1

c.. now we want to define the maximum number of the bands by all kpts

      nbnd = 0

      i_rec = 0 ; n_rank = 0
c*************************************************************
c..writing down the eig.1 and/or eig.2 files
c*************************************************************
      if(l_p0)then
         call wann_write_eig(
     >            eig_id,l_real,
     >            lmaxd,ntypd,nlod,neigd,nvd,wannierspin,
     >            isize,jspin,nbasfcn,nlotot,
     >            l_ss,l_noco,nrec,fullnkpts,
     >            wann%l_bzsym,wann%l_byindex,wann%l_bynumber,
     >            wann%l_byenergy,
     >            irreduc,odi,wann%band_min(jspin),
     >            wann%band_max(jspin),
     >            numbands,
     >            e1s,e2s,efermi,wann%l_pauli,nkpts,
     <            nbnd,kpoints,.false.,-1)
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

      allocate ( mmn(nbnd,nbnd,fullnkpts) )
      mmn(:,:,:) = cmplx(0.,0.)
      if((l_soc.or.l_noco) .and. (jspin.eq.1))
     &   allocate( socmmn(nbnd,nbnd,fullnkpts) )

      if(wann%l_nabla)then
         allocate ( nablamat(3,nbnd,nbnd,fullnkpts) )
         nablamat = cmplx(0.,0.)
      endif

      if(wann%l_surfcurr)then
         allocate ( surfcurr(3,nbnd,nbnd,fullnkpts) )
         surfcurr = cmplx(0.,0.)
      endif

      if(wann%l_socmat)then
         allocate ( hsomtx(2,2,nbnd,nbnd,fullnkpts) )
      endif

      write (*,*) 'nwfs=',nwfs

      allocate ( flo(ntypd,jmtd,2,nlod,2) )

c***********************************************************
c     Calculate the radial wave functions for the two spins.
c***********************************************************
      do jspin=1,wannierspin
        jspin2=jspin
        if(l_soc .and. jspins.eq.1)jspin2=1
        do n = 1,ntype
          do l = 0,lmax(n)
             call radfun(
     >              l,n,jspin,epar(l,n,jspin2),vr(1,n,jspin2),atoms,
     <              ff(n,:,:,l,jspin),gg(n,:,:,l,jspin),usdus,
     <              nodeu,noded,wronk)
          enddo !l
          do ilo = 1, nlo(n)

             call radflo(
     >             atoms,n,jspin,ello(:,:,jspin2),vr(1,n,jspin2),
     >             ff(n,1:,1:,0:,jspin),gg(n,1:,1:,0:,jspin),mpi,
     <             usdus,uuilon,duilon,ulouilopn,flo(n,:,:,:,jspin))

          enddo !ilo
        enddo !n
      enddo !jspin

c****************************************************************
c     Calculate the radial overlap integrals.
c****************************************************************
      if(wann%l_mmn0)then
       allocate( radial1_ff(2,2,0:lmaxd,0:lmaxd,ntype) )
       allocate( radial1_fg(2,2,0:lmaxd,0:lmaxd,ntype) )
       allocate( radial1_gf(2,2,0:lmaxd,0:lmaxd,ntype) )
       allocate( radial1_gg(2,2,0:lmaxd,0:lmaxd,ntype) )

       allocate( radial1_flo(2,2,0:lmaxd,1:nlod,ntype) )
       allocate( radial1_glo(2,2,0:lmaxd,1:nlod,ntype) )
       allocate( radial1_lof(2,2,1:nlod,0:lmaxd,ntype) )
       allocate( radial1_log(2,2,1:nlod,0:lmaxd,ntype) )
       allocate( radial1_lolo(2,2,1:nlod,1:nlod,ntype) )

       do n=1,ntype
        do l1=0,lmaxd
         do l2=0,lmaxd
          do spinloop1=1,2
           do spinloop2=1,2
            call wann_radovlp_integrals(
     >         rmsh(:,n),jri(n),dx(n),
     >         ff(n,:,:,l1,spinloop2),
     >         ff(n,:,:,l2,spinloop1),
     <         radial1_ff(spinloop2,spinloop1,l1,l2,n))
            call wann_radovlp_integrals(
     >         rmsh(:,n),jri(n),dx(n),
     >         ff(n,:,:,l1,spinloop2),
     >         gg(n,:,:,l2,spinloop1),
     <         radial1_fg(spinloop2,spinloop1,l1,l2,n))
            call wann_radovlp_integrals(
     >         rmsh(:,n),jri(n),dx(n),
     >         gg(n,:,:,l1,spinloop2),
     >         ff(n,:,:,l2,spinloop1),
     <         radial1_gf(spinloop2,spinloop1,l1,l2,n))
            call wann_radovlp_integrals(
     >         rmsh(:,n),jri(n),dx(n),
     >         gg(n,:,:,l1,spinloop2),
     >         gg(n,:,:,l2,spinloop1),
     <         radial1_gg(spinloop2,spinloop1,l1,l2,n))
           enddo !spinloop2
          enddo !spinloop1
         enddo !l2
        enddo !l1
       enddo !n

       do n=1,ntype
        do ilo=1,nlo(n)
         do ilop=1,nlo(n)
          do spinloop1=1,2
           do spinloop2=1,2
             call wann_radovlp_integrals(
     > rmsh(:,n),jri(n),dx(n),
     > flo(n,:,:,ilo,spinloop2),
     > flo(n,:,:,ilop,spinloop1),
     < radial1_lolo(spinloop2,spinloop1,ilo,ilop,n))
           enddo !spinloop2
          enddo !spinloop1
         enddo!ilop
        enddo !ilo
       enddo !n

       do n=1,ntype
        do ilo=1,nlo(n)
         do l2=0,lmaxd
          do spinloop1=1,2
           do spinloop2=1,2
             call wann_radovlp_integrals(
     > rmsh(:,n),jri(n),dx(n),
     > flo(n,:,:,ilo,spinloop2),
     > ff(n,:,:,l2,spinloop1),
     < radial1_lof(spinloop2,spinloop1,ilo,l2,n))
           enddo !spinloop2
          enddo !spinloop1
         enddo!l2
        enddo !ilo
       enddo !n

       do n=1,ntype
        do ilo=1,nlo(n)
         do l2=0,lmaxd
          do spinloop1=1,2
           do spinloop2=1,2
             call wann_radovlp_integrals(
     > rmsh(:,n),jri(n),dx(n),
     > flo(n,:,:,ilo,spinloop2),
     > gg(n,:,:,l2,spinloop1),
     < radial1_log(spinloop2,spinloop1,ilo,l2,n))
           enddo !spinloop2
          enddo !spinloop1
         enddo!l2
        enddo !ilo
       enddo !n

       do n=1,ntype
        do ilo=1,nlo(n)
         do l2=0,lmaxd
          do spinloop1=1,2
           do spinloop2=1,2
             call wann_radovlp_integrals(
     > rmsh(:,n),jri(n),dx(n),
     > ff(n,:,:,l2,spinloop2),
     > flo(n,:,:,ilo,spinloop1),
     < radial1_flo(spinloop2,spinloop1,l2,ilo,n))
           enddo !spinloop2
          enddo !spinloop1
         enddo!l2
        enddo !ilo
       enddo !n

       do n=1,ntype
        do ilo=1,nlo(n)
         do l2=0,lmaxd
          do spinloop1=1,2
           do spinloop2=1,2
             call wann_radovlp_integrals(
     > rmsh(:,n),jri(n),dx(n),
     > gg(n,:,:,l2,spinloop2),
     > flo(n,:,:,ilo,spinloop1),
     < radial1_glo(spinloop2,spinloop1,l2,ilo,n))
           enddo !spinloop2
          enddo !spinloop1
         enddo!l2
        enddo !ilo
       enddo !n

      endif


      jspin=1

c****************************************************************
c     compute matrix elements of exchange field
c****************************************************************
#ifdef CPP_TOPO
      if(wann%l_perpmag)then
         allocate( perpmag(nbnd,nbnd,fullnkpts) )
         perpmag=0.0
         allocate( ubeffug(0:lmd,0:lmd,1:ntype) )
         allocate( ubeffdg(0:lmd,0:lmd,1:ntype) )
         allocate( dbeffug(0:lmd,0:lmd,1:ntype) )
         allocate( dbeffdg(0:lmd,0:lmd,1:ntype) )
         allocate( ubeffulog(0:lmd,0:15,1:ntype) )
         allocate( dbeffulog(0:lmd,0:15,1:ntype) )
         allocate( ulobeffug(0:lmd,0:15,1:ntype) )
         allocate( ulobeffdg(0:lmd,0:15,1:ntype) )
         allocate( ulobeffulog(0:15,0:15,1:ntype) )
         call wann_gabeffgagaunt(
     >            neq,mlh,nlhd,nlh,ntypsy,llh,lmax,
     >            nmem,ntype,ntypd,bbmat,bmat,
     >            nlod,nlo,llo,flo,ff,gg,jri,rmsh,dx,jmtd,
     >            lmaxd,lmd,clnu,
     <            ubeffug,ubeffdg,dbeffug,dbeffdg,ubeffulog,
     <            dbeffulog,
     <            ulobeffug,ulobeffdg,ulobeffulog)
      endif   
#endif      


      i_rec = 0 ; n_rank = 0

c****************************************************************
c..   loop by kpoints starts!      each may be a separate task
c****************************************************************
      do 10 ikpt = wann%ikptstart,fullnkpts  ! loop by k-points starts
        kptibz=ikpt
        if(wann%l_bzsym) kptibz=irreduc(ikpt)
        if(wann%l_bzsym) oper=mapkoper(ikpt)

        i_rec = i_rec + 1
      if (mod(i_rec-1,isize).eq.irank) then

c****************************************************************
c     Obtain eigenvectors for the two spin channels.
c****************************************************************
      
      allocate ( eigg(neigd,2) )
      call cpu_time(delta)
      do jspin=1,wannierspin
         if(jspin.eq.1)nrec=0
         if(jspin.eq.2)nrec=nkpts
         if(l_noco)nrec=0
         n_start=1
         n_end=neigd

      zzMat(jspin)%l_real = l_real
      zzMat(jspin)%matsize1 = nbasfcn
      zzMat(jspin)%matsize2 = neigd
      IF(l_real) THEN
         ALLOCATE (zzMat(jspin)%data_r(zzMat(jspin)%matsize1,
     +                                 zzMat(jspin)%matsize2))
      ELSE
         ALLOCATE (zzMat(jspin)%data_c(zzMat(jspin)%matsize1,
     +                                 zzMat(jspin)%matsize2))
      END IF

!         CALL cdn_read(
!     >              eig_id,
!     >              nvd,jspd,irank,isize,kptibz,jspin,nbasfcn,
!     >              l_ss,l_noco,neigd,n_start,n_end,
!     <              nmat,nv,ello,evdu,epar,
!     <              k1,k2,k3,bkpt,wk,nbands,eigg(:,jspin),zzMat(jspin))
      call judft_error("cdn_read in wann_updown")
         nrec=0
      enddo !jspin
      call cpu_time(delta1)
      time_rw=time_rw+delta1-delta
      nslibd = 0
      
      jspin=1

c...we work only within the energy window

      DO jspin = 1, wannierspin
         zMat(jspin)%l_real = zzMat(jspin)%l_real
         zMat(jspin)%matsize1 = zzMat(jspin)%matsize1
         zMat(jspin)%matsize2 = zzMat(jspin)%matsize2
         IF (zzMat(jspin)%l_real) THEN
	       IF(.not.allocated(zmat(jspin)%data_r))
     +      ALLOCATE (zMat(jspin)%data_r(zMat(jspin)%matsize1,
     +                                zMat(jspin)%matsize2))
            zMat(jspin)%data_r = 0.0
         ELSE
	       IF(.not.allocated(zMat(jspin)%data_c))
     +      ALLOCATE (zMat(jspin)%data_c(zMat(jspin)%matsize1,
     +                                zMat(jspin)%matsize2))
            zMat(jspin)%data_c = CMPLX(0.0,0.0)
         END IF
      END DO

      eig(:) = 0.

      do jspin=1,wannierspin
       nslibd=0
       !print*,"bands used:"
       do i = 1,nbands
        if ((eigg(i,1).ge.e1s.and.nslibd.lt.numbands.
     &       and.wann%l_bynumber)
     &    .or.(eigg(i,1).ge.e1s.and.eigg(i,1).le.e2s.
     &       and.wann%l_byenergy)
     &    .or.(i.ge.wann%band_min(jspin).and.
     &        (i.le.wann%band_max(jspin))
     &       .and.wann%l_byindex))then

           !print*,i
           nslibd = nslibd + 1
           eig(nslibd) = eigg(i,1)
           if(l_noco)then
             funbas=nv(1)+nlotot
             funbas=funbas+nv(2)+nlotot
           else
             funbas=nv(jspin)+nlotot
           endif
           IF (zzMat(jspin)%l_real) THEN
              do j = 1,funbas
                 zMat(jspin)%data_r(j,nslibd) = zzMat(jspin)%data_r(j,i)
              enddo
           ELSE
              do j = 1,funbas
                 zMat(jspin)%data_c(j,nslibd) = zzMat(jspin)%data_c(j,i)
              enddo
           END IF
        endif
       enddo
      enddo
      jspin=1

c***********************************************************
c              rotate the wavefunction
c***********************************************************
      do jspin=1,wannierspin
        if (wann%l_bzsym.and.oper.ne.1) then  !rotate bkpt
          call wann_kptsrotate(
     >            natd,nlod,llod,
     >            ntypd,nlo,llo,invsat,
     >            l_noco,l_soc,
     >            ntype,neq,nlotot,
     >            jspin,
     >            oper,nop,mrot,nvd,nv,
     >            shiftkpt(:,ikpt),
     >            tau,
     x            bkpt,k1(:,:),k2(:,:),k3(:,:),
     x            zMat(jspin),nsfactor)
        else
          nsfactor=cmplx(1.0,0.0)
        endif
      enddo !jspin 
c      print*,"bkpt1=",bkpt

c******************************************************************

c...the overlap matrix Mmn which is computed for each k- and b-point

      noccbd = nslibd

      allocate ( acof(noccbd,0:lmd,natd,wannierspin),
     &           bcof(noccbd,0:lmd,natd,wannierspin),
     &           ccof(-llod:llod,noccbd,nlod,natd,wannierspin))

      acof(:,:,:,:) = cmplx(0.,0.) ; bcof(:,:,:,:) = cmplx(0.,0.)
      ccof(:,:,:,:,:) = cmplx(0.,0.)

c...generation the A,B,C coefficients in the spheres
c...for the lapws and local orbitals, summed by the basis functions

      ALLOCATE(lapw%k1(SIZE(k1,1),SIZE(k1,2)))
      ALLOCATE(lapw%k2(SIZE(k1,1),SIZE(k1,2)))
      ALLOCATE(lapw%k3(SIZE(k1,1),SIZE(k1,2)))
      lapw%nv = nv
      lapw%nmat = nmat
      lapw%k1 = k1
      lapw%k2 = k2
      lapw%k3 = k3
      ! I think the other variables of lapw are not needed here.

      call cpu_time(delta)
      do jspin=1,wannierspin

         CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,
     +              noco,jspin,oneD,acof(1:,0:,1:,jspin),
     +              bcof(1:,0:,1:,jspin),ccof(-llod:,1:,1:,1:,jspin),
     +              zMat(jspin))

         call wann_abinv(
     >       ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
     >       noccbd,lmax,nlo,llo,invsat,invsatnr,bkpt,taual,
     X       acof(1,0,1,jspin),bcof(1,0,1,jspin),
     <       ccof(-llod,1,1,1,jspin))
      enddo !jspin

      DEALLOCATE(lapw%k1,lapw%k2,lapw%k3)

      jspin=1
      call cpu_time(delta1)
      time_abcof=time_abcof+delta1-delta

      if(wann%l_socmat)then
         call wann_socmat(
     >       mpi,enpara,input,noco,atoms,
     >       lmaxd,ntypd,nlod,natd,noccbd,
     >       llod,jmtd,jspd,nlhd,neq,
     >       ntype,theta,phi,jspins,irank,
     >       jri,lmax,dx,rmsh,epar,ello,nlo,llo,
     >       l_dulo,ulo_der,vTot%mt,
     >       acof,bcof,ccof,
     <       hsomtx(:,:,:,:,ikpt))
      endif

#ifdef CPP_TOPO

      if(wann%l_perpmag)then
         call wann_perpmag(
     >            nslibd,neq,mlh,nlhd,nlh,ntypsy,llh,lmax,
     >            nmem,ntype,ntypd,bbmat,bmat,
     >            nlod,nlo,llo,flo,ff,gg,jri,rmsh,dx,jmtd,
     >            lmaxd,lmd,clnu,
     >            ubeffug,ubeffdg,dbeffug,dbeffdg,ubeffulog,
     >            dbeffulog,
     >            ulobeffug,ulobeffdg,ulobeffulog,
     >            acof,bcof,
     <            perpmag(:,:,ikpt))
      endif   

      if(wann%l_surfcurr)then
         call wann_surfcurr_int_updown(
     >       nv2d,jspin,odi,ods,n3d,nmzxyd,n2d,ntypsd,
     >       ntype,lmaxd,jmtd,ntypd,natd,nmzd,neq,nq3,nvac,
     >       nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
     >       invs,invs2,z1,delz,ngopr,ntypsy,jri,
     >       pos,taual,zatom,rmt,
     >       lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,ikpt,
     >       nvd,nv(jspin),bkpt,omtil,
     >       k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >       nslibd,nbasfcn,neigd,zMat,
     >       dirfacs,
     >       surfcurr(:,:,:,ikpt))
         call wann_surfcurr_updown(
     >       dirfacs,amat,
     >       ntype,lmaxd,lmax,natd,
     >       neq,noccbd,lmd,natd,llod,nlod,
     >       nlo,llo, 
     >       acof,bcof,ccof,
     >       us,dus,duds,uds,
     >       ulos,dulos,
     >       rmt,pos, 
     &       surfcurr(:,:,:,ikpt))
        write(6,*)"dirfacs=",dirfacs
      endif

      if(wann%l_nabla)then
         call wann_nabla_updown(
     >       jspin,ntype,lmaxd,lmax,natd,
     >       jmtd,jri,rmsh,dx,neq,
     >       noccbd,lmd,natd,llod,nlod, 
     >       ff(:,:,:,:,1),gg(:,:,:,:,1),
     >       ff(:,:,:,:,2),gg(:,:,:,:,2),
     >       acof(:,:,:,1),bcof(:,:,:,1),
     >       ccof(:,:,:,:,1),
     >       acof(:,:,:,2),bcof(:,:,:,2),
     >       ccof(:,:,:,:,2),
     &       nablamat(:,:,:,ikpt))

         addnoco=0
         do 41 i = n_rank+1,nv(jspin),n_size
           b1(1)=bkpt(1)+k1(i,jspin)
           b1(2)=bkpt(2)+k2(i,jspin)
           b1(3)=bkpt(3)+k3(i,jspin)
           b2(1)=b1(1)*bmat(1,1)+b1(2)*bmat(2,1)+b1(3)*bmat(3,1)
           b2(2)=b1(1)*bmat(1,2)+b1(2)*bmat(2,2)+b1(3)*bmat(3,2)
           b2(3)=b1(1)*bmat(1,3)+b1(2)*bmat(2,3)+b1(3)*bmat(3,3)
           do 42 j = n_rank+1,nv(jspin),n_size
c            b1(1)=bkpt(1)+k1(j,jspin)
c            b1(2)=bkpt(2)+k2(j,jspin)
c            b1(3)=bkpt(3)+k3(j,jspin)
c            b3(1)=b1(1)*bmat(1,1)+b1(2)*bmat(2,1)+b1(3)*bmat(3,1)
c            b3(2)=b1(1)*bmat(1,2)+b1(2)*bmat(2,2)+b1(3)*bmat(3,2)
c            b3(3)=b1(1)*bmat(1,3)+b1(2)*bmat(2,3)+b1(3)*bmat(3,3)
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
               IF (zMat%l_real) THEN
                  zzConjgTemp = cmplx(zMat(1)%data_r(i+addnoco,m) *
     +                                zMat(2)%data_r(j+addnoco,n), 0.0)
               ELSE
                  zzConjgTemp =       zMat(1)%data_c(i+addnoco,m) *
     +                          CONJG(zMat(2)%data_c(j+addnoco,n))
               END IF
               value=phasust*zzConjgTemp
               nablamat(dir,m,n,ikpt) = nablamat(dir,m,n,ikpt) -
     +                                  value*b2(dir)
              enddo 
             enddo
            enddo

 42        continue
 41      continue
      endif   
#endif

c***************************************************
c     Calculate matrix elements of the Pauli matrix.
c***************************************************
      if(wann%l_mmn0)then

c------> interstitial contribution to updown.mmn0-matrix
       addnoco=0
       if(l_noco)then
          addnoco=nv(1)+nlotot
       endif

c$$$       do i = n_rank+1,nv(jspin),n_size
c$$$         do 22 j = n_rank+1,nv(jspin),n_size
c$$$
c$$$c-->     determine index and phase factor
c$$$            i1 = k1(j,jspin) - k1(i,jspin)
c$$$            i2 = k2(j,jspin) - k2(i,jspin)
c$$$            i3 = k3(j,jspin) - k3(i,jspin)
c$$$            in = ig(i1,i2,i3)
c$$$            if (in.eq.0) goto 22
c$$$            phase   = rgphs(i1,i2,i3)
c$$$            phasust = cmplx(phase,0.0)*ustep(in)
c$$$            do m = 1,nslibd
c$$$             do n = 1,nslibd
c$$$#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
c$$$              mmn(m,n,ikpt) =
c$$$     =        mmn(m,n,ikpt) + 
c$$$     +          phasust*z(i,m,1)*conjg(z(j+addnoco,n,2))
c$$$#else
c$$$              mmn(m,n,ikpt) =
c$$$     =        mmn(m,n,ikpt) + 
c$$$     +          phasust*cmplx(z(i,m,1)*z(j+addnoco,n,2),0.0)
c$$$#endif
c$$$             enddo !n
c$$$            enddo !m
c$$$
c$$$  22     continue
c$$$       enddo !i  

      addnoco2 = addnoco      ! TODO: correct for addnoco2
      call wann_mmkb_int(
     >         cmplx(1.,0.),addnoco,addnoco2,
     >         nvd,k1d,k2d,k3d,
     >         n3d,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),neigd,nbasfcn,zMat(1),nslibd,
     >         k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),zMat(2),nslibd,
     >         nbnd,
     >         rgphs,ustep,ig,(/ 0,0,0 /),
     <         mmn(:,:,ikpt))


c---> spherical contribution to updown.mmn0-matrix
       call wann_mmk0_updown_sph(
     >          l_noco,alph,beta,
     >          llod,noccbd,nlod,natd,ntypd,lmaxd,lmax,lmd,
     >          ntype,neq,nlo,llo,
     >          radial1_ff,radial1_gg,radial1_fg,radial1_gf,
     >          radial1_flo,radial1_glo,
     >          radial1_lof,radial1_log,
     >          radial1_lolo,
     >          acof(1:noccbd,:,:,:),
     >          bcof(1:noccbd,:,:,:),ccof(:,1:noccbd,:,:,:),
     >          usdus%ddn,usdus%uulon,usdus%dulon,usdus%uloulopn,
     =          mmn(1:noccbd,1:noccbd,ikpt))

c---> vacuum contribution to updown.mmn0-matrix

       if (film .and. .not.odi%d1) then

c         call wann_mmk0_vac(
c     >           z1,nmzd,nv2d,k1d,k2d,k3d,n3d,nvac,
c     >           ig,nmz,delz,ig2,area,bmat,
c     >           bbmat,evac(:,jspin),bkpt,vz(:,:,jspin2),
c     >           nslibd,jspin,k1,k2,k3,wannierspin,nvd,
c     >           nbasfcn,neigd,z,nv,omtil,
c     <           mmn(1:noccbd,1:noccbd,ikpt))

       elseif (odi%d1) then

c          call wann_mmk0_od_vac(
c     >           z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
c     >           ig,nmzxy,nmz,delz,ig2,odi%n2d,
c     >           bbmat,evac(1,jspin),bkpt,odi%M,odi%mb,
c     >           vz(:,1,jspin2),odi,
c     >           nslibd,jspin,k1,k2,k3,wannierspin,nvd,area,
c     >           nbasfcn,neigd,z,nv,sk2,phi2,omtil,qss,
c     <           mmn(1:noccbd,1:noccbd,ikpt))

       endif
      endif !l_mmn0

      deallocate ( acof,bcof,ccof,eigg )

      write (*,*) 'nslibd=',nslibd


      endif   ! loop by processors

10    continue ! end of cycle by the k-points

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
     >      mpi_comm,l_p0,'nabl.updown',
     >      'Matrix elements of nabla operator',
     >      nbnd,fullnkpts,nbnd,
     >      irank,isize,
     <      nablamat)
      endif

      if(wann%l_socmat)then
        call wann_write_hsomtx(
     >      mpi_comm,l_p0,'WF1.hsomtx',
     >      'Matrix elements of spin-orbit coupling',
     >      2,2,nbnd,nbnd,fullnkpts,
     >      irank,isize,
     <      hsomtx)
      endif

      if(wann%l_mmn0)then
        call wann_write_amn(
     >      mpi_comm,l_p0,'updown.mmn0',
     >      'Overlaps of the wavefunct. at the same kpoint',
     >      nbnd,fullnkpts,nbnd,
     >      irank,isize,.false.,
     <      mmn,wann%l_unformatted)
      endif !l_soc and l_mmn0  

      if(wann%l_perpmag)then
        call wann_write_amn(
     >      mpi_comm,l_p0,'updown.perpmag',
     >      'Overlaps of the wavefunct. with Beff',
     >      nbnd,fullnkpts,nbnd,
     >      irank,isize,.true.,
     <      perpmag,wann%l_unformatted)
      endif !l_soc and l_mmn0  

      if(wann%l_surfcurr)then
        surfcurr=conjg(surfcurr/cmplx(0.0,2.0))
        hescale=tpi*condquant/bohrtocm/omtil
        hescale=sqrt(hescale)
        surfcurr=surfcurr*hescale
        call wann_write_nabla(
     >      mpi_comm,l_p0,'surfcurr.updown',
     >      'Surface currents',
     >      nbnd,fullnkpts,nbnd,
     >      irank,isize,
     <      surfcurr)
      endif 

      if( allocated (nablamat) ) deallocate( nablamat )
      deallocate ( mmn )

      deallocate(flo)
 
      deallocate ( vr,vz,nv,k1,k2,k3 )
      deallocate ( ff,gg )
      if (wann%l_bzsym)  deallocate(irreduc,mapkoper)

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

#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,ierr)
#endif


      END SUBROUTINE wann_updown
      END MODULE m_wann_updown
