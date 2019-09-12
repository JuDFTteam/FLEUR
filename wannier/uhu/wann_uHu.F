c*******************************************c
c   Set up uHu matrix necessary for         c
c   Wannier-based calc. of orbital moment   c
c*******************************************c
c   keyword is 'matrixuhu' in wann_inp      c
c*******************************************c
c   uHu = < u_{k+b1} | H_{k} | u_{k+b2} >   c
c                                           c
c   Contributions to Hamiltonian:           c
c       (i)   interstitial                  c
c       (ii)  muffin tin  (a) spherical     c
c                         (b) non-sph.      c
c                         (c) SOC           c
c       (iii) vacuum                        c
c*******************************************c
c                  J.-P. Hanke, Dec. 2015   c
c*******************************************c
      MODULE m_wann_uHu
      USE m_juDFT
      CONTAINS
      SUBROUTINE wann_uHu(
     >      DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,
     >      banddos,oneD,noco,cell,vTot,wann,
     >      eig_idList,l_real,l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,zrfs,ig2,area,tau,zatom,nq2,kv2,nop2,
     >      volint,symor,pos,ef,l_soc,
     >      memd,lnonsph,clnu,lmplmd,mlh,nmem,llh,lo1l,
     >      theta,phi,
     >      l_ms,l_sgwf,l_socgwf,aux_latt_const,
     >      param_file,param_vec,nparampts,param_alpha,l_dim)

      use m_types
      use m_wann_mmnk_symm
      use m_wann_rw_eig
      use m_abcof
      use m_radfun
      use m_radflo
      use m_cdnread
      use m_constants, only : pimach
      use m_wann_projmethod
      use m_wann_abinv
      use m_wann_kptsrotate
      use m_wann_read_inp
      use m_matmul,only : matmul3,matmul3r
      use m_wann_maxbnd
      use m_wann_uHu_tlmplm
      use m_wann_uHu_sph
      use m_wann_uHu_int
      use m_wann_uHu_soc
      use m_wann_uHu_vac
      use m_wann_uHu_od_vac
      use m_wann_uHu_util
      use m_wann_uHu_commat
      use m_wann_write_uHu
      USE m_eig66_io

      IMPLICIT NONE
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer ierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)
#endif

      TYPE(t_dimension),INTENT(IN) :: DIMENSION
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_kpts),INTENT(IN)      :: kpts
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_banddos),INTENT(IN)   :: banddos
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_potden),INTENT(IN)    :: vTot
      TYPE(t_wann),INTENT(INOUT)   :: wann

c     ..scalar arguments..
      character(len=20),intent(in) :: param_file
      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods
      INTEGER, INTENT (IN) :: eig_idList(wann%nparampts)
      logical, intent (in) :: invs,invs2,film,slice,symor,zrfs
      logical, intent (in) :: l_real,l_noco,l_ss,l_soc
      logical, intent (in) :: l_ms,l_sgwf,l_socgwf
      integer, intent (in) :: lmaxd,ntypd,neigd,nkptd,kk,nnne
      integer, intent (in) :: natd,nop,nvd,jspd,nbasfcn,nq2,nop2
      integer, intent (in) :: llod,nlod,ntype,n3d,n2d
      integer, intent (in) :: nmzxyd,nmzd,jmtd,nlhd,nq3,nvac
      integer, intent (in) :: ntypsd,jspins,k1d,k2d,k3d
      integer, intent (in) :: irank,isize,nv2d,nmzxy,nmz
      integer, intent (in) :: memd,lmplmd,nparampts
      real,    intent (in) :: omtil,e1s,e2s,delz,area,z1,volint
      real,    intent (in) :: ef,theta,phi,aux_latt_const

c     ..array arguments..
      logical, intent (in) :: l_dulo(nlod,ntypd)
      logical, intent (in) :: l_dim(3)
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: nlh(ntypsd),jri(ntypd),ntypsy(natd)
      integer, intent (in) :: nlo(ntypd),llo(nlod,ntypd),lapw_l(ntypd)
      integer, intent (in) :: invtab(nop),mrot(3,3,nop),ngopr(natd)
      integer, intent (in) :: neq(ntypd),lmax(ntypd)
      integer, intent (in) :: invsat(natd),invsatnr(natd),nkpt
      integer, intent (in) :: ulo_der(nlod,ntypd),ig2(n3d),kv2(2,n2d)
      integer, intent (in) :: mlh(memd,0:nlhd,ntypsd)
      integer, intent (in) :: nmem(0:nlhd,ntypsd)
      integer, intent (in) :: llh(0:nlhd,ntypsd),lnonsph(ntypd)
      integer, intent (in) :: lo1l(0:llod,ntypd)
      complex, intent (in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      real,    intent (in) :: taual(3,natd),rmt(ntypd),dx(ntypd)
      real,    intent (in) :: amat(3,3),bmat(3,3),bbmat(3,3)
      real,    intent (in) :: rmsh(jmtd,ntypd),tau(3,nop),zatom(ntype)
      real,    intent (in) :: alph(ntypd),beta(ntypd),qss(3)
      real,    intent (in) :: pos(3,natd),sk2(n2d),phi2(n2d)
      real,    intent (in) :: param_vec(3,nparampts)
      real,    intent (in) :: param_alpha(ntypd,nparampts)
      complex, intent (in) :: ustep(n3d)
      complex, intent (in) :: clnu(memd,0:nlhd,ntypsd)

c     ..allocatable arrays..
      integer, allocatable :: kveclo(:)   , nv(:)
      integer, allocatable :: kveclo_b(:) , nv_b(:)
      integer, allocatable :: kveclo_b2(:), nv_b2(:)
      integer, allocatable :: k1(:,:)   , k2(:,:)   , k3(:,:)
      integer, allocatable :: k1_b(:,:) , k2_b(:,:) , k3_b(:,:)
      integer, allocatable :: k1_b2(:,:), k2_b2(:,:), k3_b2(:,:)
      integer, allocatable :: irreduc(:),mapkoper(:)
      integer, allocatable :: irreduc_q(:),mapqoper(:)        
      integer, allocatable :: shiftkpt(:,:),pair_to_do(:,:)
      integer, allocatable :: shiftqpt(:,:),pair_to_do_q(:,:)  
      integer, allocatable :: maptopair(:,:,:)
      integer, allocatable :: maptopair_q(:,:,:)              
      integer, allocatable :: counts(:),displs(:)
      integer, allocatable :: gb(:,:,:),bpt(:,:)
      integer, allocatable :: gb_q(:,:,:),bpt_q(:,:)
      INTEGER, ALLOCATABLE :: innerEig_idList(:)
      real,    allocatable :: we(:),we_b(:),we_b2(:)
      real,    allocatable :: eigg(:)
      real,    allocatable :: vr(:,:,:),vz(:,:,:)
      real,    allocatable :: flo(:,:,:,:,:)
      real,    allocatable :: ff(:,:,:,:,:),gg(:,:,:,:,:)
      real,    allocatable :: kdiff(:,:),qdiff(:,:),zero_qdiff(:,:)
      complex, allocatable :: vpw(:,:),vzxy(:,:,:,:)
      complex, allocatable :: uHu(:,:,:,:,:)
c      complex, allocatable :: uHuold(:,:)
      complex, allocatable :: acof_b(:,:,:),acof_b2(:,:,:)
      complex, allocatable :: bcof_b(:,:,:),bcof_b2(:,:,:)
      complex, allocatable :: ccof_b(:,:,:,:),ccof_b2(:,:,:,:)
      complex, allocatable :: tdd(:,:,:,:,:),tdu(:,:,:,:,:)
      complex, allocatable :: tud(:,:,:,:,:),tuu(:,:,:,:,:)
      complex, allocatable :: tdulo(:,:,:,:,:,:),tuulo(:,:,:,:,:,:)
      complex, allocatable :: tulod(:,:,:,:,:,:),tulou(:,:,:,:,:,:)
      complex, allocatable :: tuloulo(:,:,:,:,:,:,:)
      complex, allocatable :: tdd_soc(:,:,:,:),tdu_soc(:,:,:,:)
      complex, allocatable :: tud_soc(:,:,:,:),tuu_soc(:,:,:,:)
      complex, allocatable :: tdulo_soc(:,:,:,:,:)
      complex, allocatable :: tuulo_soc(:,:,:,:,:)
      complex, allocatable :: tulod_soc(:,:,:,:,:)
      complex, allocatable :: tulou_soc(:,:,:,:,:)
      complex, allocatable :: tuloulo_soc(:,:,:,:,:,:)

c     ..local arrays..
      character(len=2) :: spin012(0:2)
      data spin012/'  ', '.1', '.2'/
      character(len=3) :: spin12(2)
      data   spin12/'WF1' , 'WF2'/
      character(len=8) :: name(10)
      integer :: n_bands(0:neigd),ngopr1(natd)
      real    :: bkpt(3),bkpt_b(3),bkpt_b2(3),bkrot(3)
      real    :: eig(neigd),eig_b(neigd),eig_b2(neigd)
      real    :: uuilon(nlod,ntypd),duilon(nlod,ntypd)
      real    :: ulouilopn(nlod,nlod,ntypd)
      real    :: ello(nlod,ntypd,max(2,jspd)),evac(2,max(2,jspd))
      real    :: epar(0:lmaxd,ntypd,max(2,jspd)),evdu(2,max(jspd,2))
      real    :: qpt_i(3),qptb_i(3)
      real    :: alph_i(ntypd),alphb_i(ntypd)
      real    :: beta_i(ntypd),betab_i(ntypd)
      real    :: cp_time(9)

c     ..local scalars..
      character(len=6) :: filename
      character(len=8) :: dop,iop
      character(len=12) fending
      character(len=30) fstart
      logical :: l_p0,l_bkpts,l_proj,l_file
      logical :: l_bqpts,l_gwf,l_exist,l_nocosoc,l_symcheck
      logical :: l_skip_sph,l_skip_non,l_skip_soc
      logical :: l_skip_int,l_skip_vac,l_skip_loc
      integer :: lmd,nlotot,n,iter,ikpt,ikpt_b,ikpt_b2,iqpt,iqpt_b
      integer :: addnoco,addnoco2,funbas,loplod,igvm2
      integer :: nn,nkpts,i,j,l,i_rec,m,nwf,nwfp
      integer :: jsp_start,jsp_end,nrec,nrec_b,nrec1
      integer :: nodeu,noded,n_size,na,n_rank,nbnd,numbands
      integer :: i1,i2,i3,in,lda
      integer :: nmat,nmat_b,nmat_b2,nmat_qb
      integer :: nbands,nbands_b,nbands_b2,nbands_qb
      integer :: nslibd,nslibd_b,nslibd_b2,nslibd_qb
      integer :: noccbd,noccbd_b,noccbd_b2,noccbd_qb
      integer :: kptibz,kptibz_b,kptibz_b2
      integer :: qptibz, qptibz_b
      integer :: oper,oper_b,oper_b2,oper_q, oper_qb
      integer :: nwfs,nntot,nntot_q,fullnkpts,fullnqpts
      integer :: kpt,qpt,j1,j2,j3,k,ikpt_help,iqpt_help
      integer :: wannierspin,jspin,jspin_b,jspin2
      integer :: jspin3,jspin4_b,jspin4,jspin5,tspin,tspin2
      integer :: n_start,n_end,mlotot,mlolotot,err
      integer :: mlot_d,mlolot_d,ilo,dir,length
      integer :: npotmatfile,ig3,maxvac,irec,imz,ivac,ipot
      integer :: funit_start,band_help,sign2
      integer :: doublespin,doublespin_max,nrec5
      integer :: aoff,d1,d10,d100
      integer :: eig_id
      real    :: tpi,wronk,wk,wk_b,wk_b2
      real    :: t0,t00,t1,t_myTlmplm,t_init
      real    :: t_int,t_sph,t_vac,t_abcof,t_eig,t_total
      real    :: efermi,htr2ev
      real    :: theta_i, thetab_i, phi_i, phib_i
      complex :: nsfactor,nsfactor_b,nsfactor_b2

      TYPE(t_usdus) :: usdus
      TYPE(t_mat)  :: zMat, zzMat, zMat_b, zMat_b2
      TYPE(t_lapw)  :: lapw, lapw_b, lapw_b2


c     ..initializations..
      call cpu_time(t00)

      t_init = 0.
      t_myTlmplm = 0.
      t_eig = 0.
      t_abcof = 0.
      t_int = 0.
      t_sph = 0.
      t_vac = 0.
      t_total = 0.
      htr2ev = 27.2
      nntot_q = 1
      fullnqpts = 1
      funit_start = 5000
      
    ! In some modes the code skips the Hamiltonian setup and
    ! jumps directly into the Wannier branch.
      if(.not.allocated (usdus%us))then
         call usdus%init(atoms,input%jspins)
      endif

      aoff = iachar('1')-1
      d1  = mod(irank,10)
      IF (irank < 100) THEN
        d10 = int( (irank + 0.5)/10 )
        fstart = 'eig'//achar(d10+aoff)//achar(d1+aoff)
      ELSE
        d10 = mod((irank-d1)/10,10)
        d100 = (irank-10*d10-d1)/100
        IF ( d100.GE.10 ) d100 = d100 + iachar('7')
        fstart =
     +  'eig'//achar(d100+aoff)//achar(d10+aoff)//achar(d1+aoff)
      ENDIF


      ngopr1(:)=1

      l_p0 = .false.
      if (irank.eq.0) l_p0 = .true.
      l_nocosoc = l_noco.or.l_soc

      tpi = 2* pimach()
      lmd = lmaxd*(lmaxd+2)

!!!   should be changed in case the windows are really used
      nkpts = nkpt

      ! do we have to construct GWF ?
      l_gwf = .false.
      l_gwf = l_sgwf.or.l_socgwf 
      

c-----read the input file to determine what to do
      call wann_read_inp(input,l_p0,wann)

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
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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


      if(l_gwf) fullnqpts = nparampts


      nrec = 0
      if(l_p0)then
      write (*,*) 'fermi energy:',efermi
      write (*,*) 'emin,emax=',e1s,e2s
      write (*,*) 'nbasfcn =',nbasfcn
      endif
      nlotot = 0
      mlotot = 0
      mlolotot = 0
      do n = 1, ntype
        mlotot = mlotot + nlo(n)
        mlolotot = mlolotot + nlo(n)*(nlo(n)+1)/2
        do l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        enddo
      enddo


      allocate(counts(0:isize-1),displs(0:isize-1))
      call array_split(fullnkpts,isize,counts,displs)

c**********************************************************
ccccccccccccccc   read in the bkpts file  ccccccccccccccccc
c**********************************************************
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
         if (bpt(nn,ikpt)>fullnkpts) CALL juDFT_error("bpt.gt.fullnkpts"
     +        ,calledby ="wannier")
        enddo
       enddo
       close (202)
        allocate(kdiff(3,nntot))

c**********************************************************
ccccccccccccccc   read in the bqpts file  ccccccccccccccccc         
c**********************************************************
      if (l_gwf.or.l_ms) then ! for Omega functional minimization
       l_bqpts = .false.
       inquire (file='bqpts',exist=l_bqpts)
       if (.not.l_bqpts)  CALL juDFT_error("need bqpts for matrixmmn"
     +     ,calledby ="wannier")
       open (202,file='bqpts',form='formatted',status='old')
       rewind (202)
       read (202,'(i4)') nntot_q
       if(l_p0)then
       write (*,*) 'nntot_q=',nntot_q
       write(*,*) 'fullnqpts=',fullnqpts
       endif
       allocate ( gb_q(1:3,1:nntot_q,1:fullnqpts),
     &            bpt_q(1:nntot_q,1:fullnqpts))
       do iqpt=1,fullnqpts
        do nn=1,nntot_q
         read (202,'(2i6,3x,3i4)')
     &     iqpt_help,bpt_q(nn,iqpt),(gb_q(i,nn,iqpt),i=1,3)
         if (iqpt/=iqpt_help)  CALL juDFT_error("iqpt.ne.iqpt_help"
     +        ,calledby ="wannier")
         if (bpt_q(nn,iqpt)>fullnqpts)
     &        CALL juDFT_error("bpt_q.gt.fullnqpts",calledby ="wannier")
        enddo
       enddo
       close (202)
        allocate(qdiff(3,nntot_q))
        allocate(zero_qdiff(3,nntot_q))
        zero_qdiff=0.0
      endif                                                        


! when treating gen. WF for spin spirals, the Brillouin zone
! of q-points is twice as large compared to k-BZ. Thus,
! the G-vectors connecting neighbors across the boundary
! need to be doubled
      if(l_sgwf) gb_q = 2*gb_q    
      if(l_socgwf) gb_q = 2*gb_q 

      if(wann%l_finishgwf) goto 9110
c********************************************************
c      find symmetry-related elements in mmkb
c********************************************************
         allocate(maptopair(3,fullnkpts,nntot))
         allocate(pair_to_do(fullnkpts,nntot))
         call wann_mmnk_symm(input,kpts,
     >     fullnkpts,nntot,bpt,gb,wann%l_bzsym,
     >     irreduc,mapkoper,l_p0,film,nop,invtab,mrot,odi%d1,
     >     tau,
     <     pair_to_do,maptopair,kdiff,.false.,param_file)

      ! do the same for q-points to construct GWFs
      if(l_gwf)then 
         allocate(maptopair_q(3,fullnqpts,nntot_q))
         allocate(pair_to_do_q(fullnqpts,nntot_q))
         call wann_mmnk_symm(input,kpts,
     >     fullnqpts,nntot_q,bpt_q,gb_q,wann%l_bzsym,
     >     irreduc_q,mapqoper,l_p0,.false.,1,invtab(1),mrot(:,:,1),
     >     .false.,tau,
     <     pair_to_do_q,maptopair_q,qdiff,.true.,param_file)
      endif


c*********************************************************
cccccccccccccccc   initialize the potential   cccccccccccc
c*********************************************************

      if(.not. l_noco) then
       allocate ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
      else
       allocate ( vpw(n3d,4),vzxy(nmzxyd,odi%n2d-1,2,4) )
      endif

      allocate (vz(nmzd,2,4))
      allocate (vr(jmtd,ntypd,jspd))

      vpw(:,1:SIZE(vTot%pw,2)) = vTot%pw(:,1:SIZE(vTot%pw,2))
      IF(film) THEN
         vz(:,:,1:SIZE(vTot%vacz,3)) =
     &      vTot%vacz(:,:,1:SIZE(vTot%vacz,3))
         vzxy(:,:,:,1:SIZE(vTot%vacxy,4)) =
     &      vTot%vacxy(:,:,:,1:SIZE(vTot%vacxy,4))
      END IF

      do jspin = 1,jspins
        do n = 1, ntype
          do j = 1,jri(n)
            vr(j,n,jspin) = vTot%mt(j,0,n,jspin)
          enddo
        enddo
      enddo

      if(.not.film) deallocate(vz,vzxy)

      if(l_noco)then
         npotmatfile=25

         OPEN (npotmatfile,FILE='potmat',FORM='unformatted',
     +                  STATUS='old')
c--->    load the interstitial potential
         if(l_p0 .and. nq3.ne.n3d) then
          write(*,*)'WARNING for reading potmat: nq3.ne.n3d'
         endif
         vpw = cmplx(0.0, 0.0)
         READ (npotmatfile) (vpw(ig3,1),ig3=1,n3d)
         READ (npotmatfile) (vpw(ig3,2),ig3=1,n3d)
         READ (npotmatfile) (vpw(ig3,3),ig3=1,n3d)
         vpw(:,4) = conjg(vpw(:,3))
         if(film) then
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

c--->       load the warping part of the potential
            if(.not.odi%d1)then
               DO ipot = 1,3
                  READ (npotmatfile)((vzxy(imz,igvm2,ivac,ipot),
     +                           imz=1,nmzxy),igvm2=1,nq2-1)
               ENDDO   
            else
               DO ipot = 1,3
                  READ (npotmatfile)((vzxy(imz,igvm2,ivac,ipot),
     +                           imz=1,nmzxy),igvm2=1,odi%n2d-1)
               ENDDO   
            endif   
            vzxy(:,:,:,4) = conjg(vzxy(:,:,:,3))
          enddo   
         endif
         CLOSE (npotmatfile)
      endif   
   

      if(film .and. l_p0) write(*,*)'nvac',nvac

cccccccccccccccc   end of the potential part  ccccccccccc
      wannierspin=jspd
      if(l_soc) wannierspin=2
     
      allocate(flo(ntypd,jmtd,2,nlod,2))
      allocate ( ff(ntypd,jmtd,2,0:lmaxd,2) )
      allocate ( gg(ntypd,jmtd,2,0:lmaxd,2) )
      allocate ( kveclo(nlotot),nv(wannierspin) )
      allocate ( kveclo_b(nlotot),nv_b(wannierspin) )
      allocate ( kveclo_b2(nlotot),nv_b2(wannierspin) )
      allocate ( k1(nvd,wannierspin),k2(nvd,wannierspin),
     &           k3(nvd,wannierspin) )
      allocate ( k1_b(nvd,wannierspin),k2_b(nvd,wannierspin),
     &           k3_b(nvd,wannierspin) )
      allocate ( k1_b2(nvd,wannierspin),k2_b2(nvd,wannierspin),
     &           k3_b2(nvd,wannierspin) )

      if(l_nocosoc) then
         doublespin_max=4
      else
         doublespin_max=wannierspin
      endif
      
c      if(l_soc.and.(jspins.eq.1)) doublespin_max=wannierspin

      l_skip_int = .false.; l_skip_soc = .false.; l_skip_vac = .false.
      l_skip_sph = .false.; l_skip_non = .false.; l_skip_loc = .false.
      inquire(file='debug_uHu',exist=l_exist)
      if(l_exist) then
       open(888,file='debug_uHu')
       read(888,*)l_skip_int
       read(888,*)l_skip_sph
       read(888,*)l_skip_non
       read(888,*)l_skip_soc
       read(888,*)l_skip_loc
       read(888,*)l_skip_vac
       read(888,*)doublespin_max 
       close(888)
       if(l_p0) then
        write(*,*)'skip INT :',l_skip_int
        write(*,*)'skip SPH :',l_skip_sph
        write(*,*)'skip NON :',l_skip_non
        write(*,*)'skip SOC :',l_skip_soc
        write(*,*)'skip LOC :',l_skip_loc
        write(*,*)'skip VAC :',l_skip_vac
        write(*,*)'doublespin_max:',doublespin_max
       endif
      endif

      tspin = doublespin_max
      if(fullnqpts.eq.1) tspin=1
      if(l_p0) write(*,*)'tspin',tspin
      allocate( tdd(0:lmd,0:lmd,ntypd,nntot*nntot,tspin) )
      allocate( tdu(0:lmd,0:lmd,ntypd,nntot*nntot,tspin) )
      allocate( tud(0:lmd,0:lmd,ntypd,nntot*nntot,tspin) )
      allocate( tuu(0:lmd,0:lmd,ntypd,nntot*nntot,tspin) )
      allocate( tdulo(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot,tspin) )
      allocate( tuulo(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot,tspin) )
      allocate( tulou(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot,tspin) )
      allocate( tulod(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot,tspin) )
      allocate( tuloulo(nlod,-llod:llod,nlod,-llod:llod,
     >                  ntypd,nntot*nntot,tspin) )
      allocate( tdd_soc(0:lmd,0:lmd,ntypd,nntot*nntot) )
      allocate( tdu_soc(0:lmd,0:lmd,ntypd,nntot*nntot) )
      allocate( tud_soc(0:lmd,0:lmd,ntypd,nntot*nntot) )
      allocate( tuu_soc(0:lmd,0:lmd,ntypd,nntot*nntot) )
      allocate( tdulo_soc(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot) )
      allocate( tuulo_soc(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot) )
      allocate( tulou_soc(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot) )
      allocate( tulod_soc(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot) )
      allocate( tuloulo_soc(nlod,-llod:llod,nlod,-llod:llod,
     >                  ntypd,nntot*nntot) )

      tuu = cmplx(0.,0.); tdu = cmplx(0.,0.)
      tud = cmplx(0.,0.); tdd = cmplx(0.,0.)
      tuulo = cmplx(0.,0.); tdulo = cmplx(0.,0.)
      tulou = cmplx(0.,0.); tulod = cmplx(0.,0.)
      tuloulo = cmplx(0.,0.)
      tuu_soc = cmplx(0.,0.); tdu_soc = cmplx(0.,0.)
      tud_soc = cmplx(0.,0.); tdd_soc = cmplx(0.,0.)
      tuulo_soc = cmplx(0.,0.); tdulo_soc = cmplx(0.,0.)
      tulou_soc = cmplx(0.,0.); tulod_soc = cmplx(0.,0.)
      tuloulo_soc = cmplx(0.,0.)

      call cpu_time(t1)
      t_init = t1-t00

c*****************************************************************c
c                         START Q LOOP                            c
c*****************************************************************c
      do 314 iqpt = 1,fullnqpts  ! loop by q-points starts

       ALLOCATE(innerEig_idList(nntot_q))

        qptibz=iqpt                          
        if(wann%l_bzsym .AND. l_gwf) qptibz=irreduc_q(iqpt)
        if(wann%l_bzsym .AND. l_gwf) oper_q=mapqoper(iqpt)

       qpt_i = qss
       alph_i = alph
       beta_i = beta
       theta_i = theta
       phi_i = phi
       if(l_sgwf.or.l_ms) then
          qpt_i(:) = param_vec(:,qptibz)
          alph_i(:) = param_alpha(:,qptibz)
       elseif(l_socgwf) then 
          if(l_dim(2)) phi_i = tpi*param_vec(2,qptibz)
          if(l_dim(3)) theta_i = tpi*param_vec(3,qptibz)
       endif

       if(l_p0) then
        write(*,*)'qpt_i',qpt_i
        do n=1,ntype
         write(*,*)'n=',n,alph_i(n),beta_i(n)
        enddo
        write(*,*)'theta_i=',theta_i
        write(*,*)'phi_i=',phi_i
       endif

       IF (l_gwf) THEN
          do iqpt_b=1,nntot_q

            innerEig_idList(iqpt_b) = eig_idList(bpt_q(iqpt_b,iqpt))

!            WRITE(fending,'("_",i4.4)')bpt_q(iqpt_b,iqpt)
!            innerEig_idList(iqpt_b)=open_eig(mpi%mpi_comm,
!     +                  DIMENSION%nbasfcn,DIMENSION%neigd,
!     +                  nkpts,wannierspin,atoms%lmaxd,
!     +                  atoms%nlod,atoms%ntype,atoms%nlotot,
!     +                  l_noco,.FALSE.,l_real,l_soc,.FALSE.,
!     +                  mpi%n_size,filename=trim(fstart)//fending,
!     +                  layers=vacuum%layers,nstars=vacuum%nstars,
!     +                  ncored=DIMENSION%nstd,nsld=atoms%nat,
!     +                  nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,
!     +                  l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

          enddo

          eig_id = eig_idList(qptibz)

!          WRITE(fending,'("_",i4.4)')qptibz
!         eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,
!     +                   nkpts,wannierspin,atoms%lmaxd,
!     +                   atoms%nlod,atoms%ntype,atoms%nlotot,
!     +                   l_noco,.FALSE.,l_real,l_soc,.FALSE.,
!     +                   mpi%n_size,filename=trim(fstart)//fending,
!     +                   layers=vacuum%layers,nstars=vacuum%nstars,
!     +                   ncored=DIMENSION%nstd,nsld=atoms%nat,
!     +                   nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,
!     +                   l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

!        OPEN(66,file=trim(fstart)//fending,access='direct',
!     >          form='unformatted',recl=irecl,status='old')
       ELSEIF(l_ms) THEN

          eig_id = eig_idList(qptibz)

!        WRITE(fending,'("_",i4.4)')qptibz
!        eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,
!     +                  nkpts,wannierspin,atoms%lmaxd,
!     +                  atoms%nlod,atoms%ntype,atoms%nlotot,
!     +                  l_noco,.FALSE.,l_real,l_soc,.FALSE.,
!     +                  mpi%n_size,filename=trim(fstart)//fending,
!     +                  layers=vacuum%layers,nstars=vacuum%nstars,
!     +                  ncored=DIMENSION%nstd,nsld=atoms%nat,
!     +                  nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,
!     +                  l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

!        OPEN(66,file=trim(fstart)//fending,access='direct',
!     >          form='unformatted',recl=irecl,status='old')
       ELSE
          fending=''
       ENDIF ! l_gwf.or.l_ms
       nrec=0
       nrec_b=0

c****************************************************
c cycle by spins starts! 
c****************************************************
      do 110 doublespin=1,doublespin_max   ! cycle by spins
         if(l_p0) write(*,*)'spin loop:',doublespin

c         jspin=mod(doublespin+1,2)+1
c         jspin_b=jspin
c         if(doublespin.eq.3) jspin_b=2
c         if(doublespin.eq.4) jspin_b=1
         jspin_b=mod(doublespin+1,2)+1
         jspin=jspin_b
         if(doublespin.eq.3) jspin=2
         if(doublespin.eq.4) jspin=1

         tspin2 = doublespin
         if(fullnqpts.eq.1) tspin2 = 1

         nrec_b = nrec

         if(.not.l_noco) then
            nrec = (jspin-1)*nkpts
            nrec_b = (jspin_b-1)*nkpts
         endif

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
          CALL juDFT_error("no proj/proj.1/proj.2",calledby ="wannier")
       endif  


       jspin2=jspin
       if(l_soc .and. jspins.eq.1)jspin2=1
       jsp_start = jspin ; jsp_end = jspin

cccccccccccc   read in the eigenvalues and vectors   cccccc
       do jspin5=1,wannierspin!1,2
       jsp_start=jspin5; jsp_end=jspin5
       nrec5=0
       if(.not.l_noco) nrec5 = (jspin5-1)*nkpts
       CALL judft_error("TODO:wann_uHu")
!       call cdn_read0(eig_id,irank,isize,jspin5,wannierspin,l_noco,
!     <                ello,evac,epar,wk,n_bands,n_size)

       enddo

c..   now we want to define the maximum number of the bands by all kpts
      nbnd = 0
      i_rec = 0 ; n_rank = 0

      if(l_p0)then         
         call wann_maxbnd(
     >            eig_id,l_real,
     >            lmaxd,ntypd,nlod,neigd,nvd,wannierspin,
     >            isize,jspin,nbasfcn,nlotot,
     >            l_ss,l_noco,nrec,fullnkpts,
     >            wann%l_bzsym,wann%l_byindex,wann%l_bynumber,
     >            wann%l_byenergy,
     >            irreduc,odi,wann%band_min(jspin),
     >            wann%band_max(jspin),
     >            numbands,e1s,e2s,efermi,nkpts,
     <            nbnd,l_gwf,iqpt)       
      endif!l_p0

! nbnd is calculated for process zero and is sent here to the others
#ifdef CPP_MPI
      if(l_p0)then
         do cpu_index=1,isize-1
      call MPI_SEND(nbnd,1,MPI_INTEGER,cpu_index,1,MPI_COMM_WORLD,ierr)
         enddo
      else
       call MPI_RECV(nbnd,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,stt,ierr)
      endif
#endif
     
c##################################################################
         if(.not.allocated(uHu)) then
            allocate(uHu(nbnd,nbnd,nntot,nntot,counts(irank)))
c            allocate(uHuold(nbnd,nbnd))
            uHu = cmplx(0.,0.)
c            uHuold = cmplx(0.,0.)
         endif



! set up T-matrix for spherical and non-spherical part
! only once since it is k- and q-independent !
      if(iqpt.eq.1) then
         if(l_p0) write(*,*)'tspin2=',tspin2
         tuu(:,:,:,:,tspin2) = cmplx(0.,0.)
         tdu(:,:,:,:,tspin2) = cmplx(0.,0.)
         tud(:,:,:,:,tspin2) = cmplx(0.,0.)
         tdd(:,:,:,:,tspin2) = cmplx(0.,0.)
         tuulo(:,:,:,:,:,tspin2) = cmplx(0.,0.)
         tdulo(:,:,:,:,:,tspin2) = cmplx(0.,0.)
         tulou(:,:,:,:,:,tspin2) = cmplx(0.,0.)
         tulod(:,:,:,:,:,tspin2) = cmplx(0.,0.)
         tuloulo(:,:,:,:,:,:,tspin2) = cmplx(0.,0.)

      do jspin4=1,2
         jspin3=jspin4
         if(jspins.eq.1) jspin3=1
      na = 1
      do 40 n = 1,ntype
       do 30 l = 0,lmax(n)
c...compute the l-dependent, k-independent radial MT- basis functions

         call radfun(
     >              l,n,jspin4,epar(l,n,jspin3),vr(1,n,jspin3),atoms,
     <              ff(n,:,:,l,jspin4),gg(n,:,:,l,jspin4),usdus,
     <              nodeu,noded,wronk)

   30  continue
c...and the local orbital radial functions
       do ilo = 1, nlo(n)

         call radflo(
     >             atoms,n,jspin4,ello(:,:,jspin3),vr(1,n,jspin3),
     >             ff(n,1:,1:,0:,jspin4),gg(n,1:,1:,0:,jspin4),mpi,
     <             usdus,uuilon,duilon,ulouilopn,flo(n,:,:,:,jspin4))

       enddo
c       na = na + neq(n)
   40 continue
      enddo!jspin3


      if(jspin.eq.jspin_b) then
         if(l_p0) write(*,*)'wann_uHu_tlmplm'
         jspin3=jspin
         if(jspins.eq.1) jspin3=1
         call cpu_time(t0)
         call wann_uHu_tlmplm(
     >        memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd,ntype,dx,
     >        rmsh,jri,lmax,ntypsy,natd,lnonsph,lmd,lmplmd,clnu,
     >        mlh,nmem,llh,nlh,neq,irank,mlotot,mlolotot,
     >        vTot%mt(:,:,:,jspin3),nlod,llod,loplod,ello(1,1,jspin3),
     >        llo,nlo,lo1l,l_dulo,ulo_der,ff(:,:,:,:,jspin),
     >        gg(:,:,:,:,jspin),flo(:,:,:,:,jspin),
     >        ff(:,:,:,:,jspin_b),gg(:,:,:,:,jspin_b),
     >        flo(:,:,:,:,jspin_b),kdiff,kdiff,nntot,nntot,bmat,bbmat,
     >        vr(1,1,jspin3),epar(0,1,jspin3),invsat,
     >        l_skip_sph,l_skip_non,l_skip_loc,
     <        tuu(:,:,:,:,tspin2),tud(:,:,:,:,tspin2),
     >        tdu(:,:,:,:,tspin2),tdd(:,:,:,:,tspin2),
     >        tuulo(:,:,:,:,:,tspin2),tulou(:,:,:,:,:,tspin2),
     >        tdulo(:,:,:,:,:,tspin2),tulod(:,:,:,:,:,tspin2),
     >        tuloulo(:,:,:,:,:,:,tspin2))
         call cpu_time(t1)
         t_myTlmplm = t_myTlmplm + t1-t0
      endif

      endif!iqpt.eq.1


! compute SOC-contribution for each theta_i !
      if(l_soc.and. (.not.l_skip_soc)) then
         if(l_p0) write(*,*)'wann_uHu_soc'
         tuu_soc = cmplx(0.,0.); tdu_soc = cmplx(0.,0.)
         tud_soc = cmplx(0.,0.); tdd_soc = cmplx(0.,0.)
         tuulo_soc = cmplx(0.,0.); tdulo_soc = cmplx(0.,0.)
         tulou_soc = cmplx(0.,0.); tulod_soc = cmplx(0.,0.)
         tuloulo_soc = cmplx(0.,0.)
         call cpu_time(t0)
         call wann_uHu_soc(
     >     input,atoms,
     >     ntypd,jmtd,lmaxd,jspd,
     >     ntype,dx,rmsh,jri,lmax,natd,
     >     lmd,lmplmd,neq,irank,
     >     nlod,llod,loplod,ello,llo,nlo,lo1l,l_dulo,ulo_der,
     >     ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),flo(:,:,:,:,jspin),
     >     ff(:,:,:,:,jspin_b),gg(:,:,:,:,jspin_b),flo(:,:,:,:,jspin_b),
     >     kdiff,kdiff,nntot,nntot,bmat,bbmat,
     >     vr,epar,jspin,jspin_b,jspins,
     >     .true.,theta_i,phi_i,alph_i,beta_i,l_noco,l_skip_loc,
     <     tuu_soc,tud_soc,
     >     tdu_soc,tdd_soc,
     >     tuulo_soc,tulou_soc,
     >     tdulo_soc,tulod_soc,
     >     tuloulo_soc)
         call cpu_time(t1)
         t_myTlmplm = t_myTlmplm + t1-t0
      endif!l_soc


      i_rec = 0 ; n_rank = 0
      
c****************************************************************
c.. loop by kpoints starts!      each may be a separate task
c****************************************************************
      if(l_p0) write(*,*)'start k-loop'
      do 10 ikpt = wann%ikptstart,fullnkpts  ! loop by k-points starts
        kptibz=ikpt
        if(wann%l_bzsym) kptibz=irreduc(ikpt)
        if(wann%l_bzsym) oper=mapkoper(ikpt)

        if(kpt_on_node(ikpt,isize,counts,displs).eq.irank) then
c        write(*,*)irank,ikpt
        i_rec = i_rec + 1
c      if (mod(i_rec-1,isize).eq.irank) then

      allocate ( eigg(neigd) )

          n_start=1
          n_end=DIMENSION%neigd
 
      call cpu_time(t0)
      ! get current bkpt vector

          zzMat%l_real = l_real
          zzMat%matsize1 = DIMENSION%nbasfcn
          zzMat%matsize2 = DIMENSION%neigd
          IF(l_real) THEN
             IF(.NOT.ALLOCATED(zzMat%data_r))
     &               ALLOCATE (zzMat%data_r(zzMat%matsize1,zzMat%matsize2))
          ELSE
               IF(.NOT.ALLOCATED(zzMat%data_c))
     &          ALLOCATE (zzMat%data_c(zzMat%matsize1,zzMat%matsize2))
          END IF



    
          zMat%l_real = zzMat%l_real
          zMat%matsize1 = zzMat%matsize1
          zMat%matsize2 = zzMat%matsize2
          IF (zzMat%l_real) THEN
             IF(.NOT.ALLOCATED(zMat%data_r))
     &              ALLOCATE (zMat%data_r(zMat%matsize1,zMat%matsize2))
             zMat%data_r = 0.0
          ELSE
             IF(.NOT.ALLOCATED(zMat%data_c))
     &              ALLOCATE (zMat%data_c(zMat%matsize1,zMat%matsize2))
             zMat%data_c = CMPLX(0.0,0.0)
          END IF
    
          zMat_b%l_real = zzMat%l_real
          zMat_b%matsize1 = zzMat%matsize1
          zMat_b%matsize2 = zzMat%matsize2
          IF (zzMat%l_real) THEN
               IF(.NOT.ALLOCATED(zMat_b%data_r))
     &           ALLOCATE (zMat_b%data_r(zMat_b%matsize1,zMat_b%matsize2))
             zMat_b%data_r = 0.0
          ELSE
               IF(.NOT.ALLOCATED(zMat_b%data_c))
     &           ALLOCATE (zMat_b%data_c(zMat_b%matsize1,zMat_b%matsize2))
             zMat_b%data_c = CMPLX(0.0,0.0)
          END IF
    
    
          zMat_b2%l_real = zzMat%l_real
          zMat_b2%matsize1 = zzMat%matsize1
          zMat_b2%matsize2 = zzMat%matsize2
          IF (zzMat%l_real) THEN
               IF(.NOT.ALLOCATED(zMat_b2%data_r))
     &           ALLOCATE (zMat_b2%data_r(zMat_b%matsize1,zMat_b%matsize2))
             zMat_b2%data_r = 0.0
          ELSE
               IF(.NOT.ALLOCATED(zMat_b2%data_c))
     &           ALLOCATE (zMat_b2%data_c(zMat_b%matsize1,zMat_b%matsize2))
             zMat_b2%data_c = CMPLX(0.0,0.0)
          END IF
    
    
    

!      CALL cdn_read(
!     >              eig_id,
!     >              nvd,jspd,irank,isize,kptibz,jspin,nbasfcn, !wannierspin instead of jspd?
!     >              l_ss,l_noco,neigd,n_start,n_end,
!     <              nmat,nv,ello,evdu,epar,
!     <              wk,nbands,eigg,zzMat)



        CALL lapw%init(input,noco,kpts,atoms,sym,kptibz,cell,
     &    (sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9)
     &     .AND..NOT.noco%l_noco.and.mpi%n_size==1),mpi)

        CALL cdn_read(
     &      eig_id,
     &      DIMENSION%nvd,input%jspins,mpi%irank,mpi%isize, !wannierspin instead of DIMENSION%jspd?&
     &      kptibz,jspin,DIMENSION%nbasfcn,
     &      noco%l_ss,noco%l_noco,DIMENSION%neigd,
     &      n_start,n_end,
     &      nbands,eigg,zzMat)



        do i = 1,nbands
          if((eigg(i).ge.e1s.and.nslibd.lt.numbands
     & .and.wann%l_bynumber)
     &.or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.wann%l_byenergy)
     &.or.(i.ge.wann%band_min(jspin)
     &                 .and.
     &       (i.le.wann%band_max(jspin))
     &                 .and.
     &       wann%l_byindex))then
            nslibd = nslibd + 1
            eig(nslibd) = eigg(i)
            we(nslibd) = we(i)
            if(l_noco)then
               funbas =        nv(1) + nlotot
               funbas = funbas+nv(2) + nlotot
            else
               funbas = nv(jspin) + nlotot
            endif
            IF (zzMat%l_real) THEN
               do j = 1,funbas
                  zMat%data_r(j,nslibd) = zzMat%data_r(j,i)
               enddo
            ELSE
               do j = 1,funbas
                  zMat%data_c(j,nslibd) = zzMat%data_c(j,i)
               enddo
            END IF
          endif
        enddo




      call cpu_time(t1)
      t_eig = t_eig + t1 - t0

 
      allocate (we_b(neigd), we_b2(neigd))

  !!! the cycle by the nearest neighbors (nntot) for each kpoint

       do 15  ikpt_b = 1,nntot
          kptibz_b=bpt(ikpt_b,ikpt)

          if(wann%l_bzsym) oper_b=mapkoper(kptibz_b)
          if (wann%l_bzsym) kptibz_b=irreduc(kptibz_b)

          n_start=1
          n_end=neigd

        eigg = 0.
        call cpu_time(t0)
        call judft_error("BUG not implement")
!        CALL cdn_read(
!     >                eig_id,
!     >                nvd,jspd,irank,isize,kptibz_b,jspin,nbasfcn, !wannierspin instead of jspd?
!     >                l_ss,l_noco,neigd,n_start,n_end,
!     <                nmat_b,nv_b,ello,evdu,epar,
!     <                wk_b,nbands_b,eigg,zzMat)



        call lapw_b%init(input,noco,kpts,atoms,sym,kptibz_b,cell,
     &          (sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9)
     &           .AND..NOT.noco%l_noco.and.mpi%n_size==1),mpi)
        CALL cdn_read(
     &        eig_id,
     &      DIMENSION%nvd,input%jspins,mpi%irank,mpi%isize, !wannierspin instead of DIMENSION%jspd?
     &    kptibz_b,jspin,DIMENSION%nbasfcn,
     &    noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,n_end,
     &   nbands_b,eigg,zzMat)




        nslibd_b = 0

      IF(zzMat%l_real) THEN
         zMat_b%data_r = 0.0
      ELSE
         zMat_b%data_c = CMPLX(0.0,0.0)
      END IF

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
            IF (zzMat%l_real) THEN
               do j = 1,funbas
                  zMat_b%data_r(j,nslibd_b) = zzMat%data_r(j,i)
               enddo
            ELSE
               do j = 1,funbas
                  zMat_b%data_c(j,nslibd_b) = zzMat%data_c(j,i)
               enddo
            END IF
          endif
        enddo

c***********************************************************
c              Rotate the wavefunction of next neighbor.
c***********************************************************
      if (wann%l_bzsym .and. (oper_b.ne.1)  ) then
!         TODO
!         call wann_kptsrotate(
!     >            natd,nlod,llod,
!     >            ntypd,nlo,llo,invsat,
!     >            l_noco,l_soc,
!     >            ntype,neq,nlotot,
!     >            kveclo_b,jspin,
!     >            oper_b,nop,mrot,nvd,
!     >            nv_b,
!     >            shiftkpt(:,bpt(ikpt_b,ikpt)),
!     >            tau,
!     x            bkpt_b,k1_b(:,:),
!     x            k2_b(:,:),k3_b(:,:),
!     x            zMat_b,nsfactor_b)
      else
         nsfactor_b=cmplx(1.0,0.0)
      endif
        noccbd_b = nslibd_b
      call cpu_time(t1)
      t_eig = t_eig + t1 - t0

      addnoco = 0
      if(l_noco.and.jspin.eq.2) addnoco=nv_b(1)+nlotot

      ! set up a(k+b1),b(k+b1),c(k+b1)
      allocate( acof_b(noccbd_b,0:lmd,natd),
     >          bcof_b(noccbd_b,0:lmd,natd),
     >          ccof_b(-llod:llod,noccbd_b,nlod,natd) )

      call cpu_time(t0)

!      ALLOCATE(lapw_b%k1(SIZE(k1_b,1),SIZE(k1_b,2)))
!      ALLOCATE(lapw_b%k2(SIZE(k1_b,1),SIZE(k1_b,2)))
!      ALLOCATE(lapw_b%k3(SIZE(k1_b,1),SIZE(k1_b,2)))
!      lapw_b%k1 = k1_b
!      lapw_b%k2 = k2_b
!      lapw_b%k3 = k3_b
!      lapw_b%nmat = nmat_b
!      lapw_b%nv = nv_b
      ! I think the other variables of lapw are not needed here.
!TODO
!      CALL abcof(input,atoms,noccbd_b,sym,cell,bkpt_b,lapw_b,
!     +           noccbd_b,usdus,noco,jspin,kveclo_b,oneD,
!     +           acof_b,bcof_b,ccof_b,zMat_b)

!      DEALLOCATE(lapw_b%k1,lapw_b%k2,lapw_b%k3)

      call wann_abinv(atoms,
     X        acof_b,bcof_b,ccof_b)
      call cpu_time(t1)
      t_abcof = t_abcof + t1 - t0

       do 25  ikpt_b2 = 1,nntot!ikpt_b!nntot

          kptibz_b2=bpt(ikpt_b2,ikpt)
          if(wann%l_bzsym) oper_b2=mapkoper(kptibz_b2)
          if (wann%l_bzsym) kptibz_b2=irreduc(kptibz_b2)

          n_start=1
          n_end=DIMENSION%neigd

          eigg = 0.
          call cpu_time(t0)

          WRITE(*,*) 'Here probably the wrong record is read in'
          WRITE(*,*) 'Should eig_id not be dependent on iqpt_b?'
          WRITE(*,*) '(in wann_uHu)'
!          CALL cdn_read(
!     >                  eig_id,
!     >                  nvd,jspd,irank,isize,kptibz_b2,jspin,nbasfcn, !wannierspin instead of jspd?
!     >                  l_ss,l_noco,neigd,n_start,n_end,
!     <                  nmat_b2,nv_b2,ello,evdu,epar,kveclo_b2,
!     <                  k1_b2,k2_b2,k3_b2,bkpt_b2,wk_b2,nbands_b2,
!     <                  eigg,zzMat)

          call lapw_b2%init(input,noco,kpts,atoms,sym,kptibz_b2,
     &     cell,(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9)
     &     .AND..NOT.noco%l_noco.and.mpi%n_size==1),mpi)

          CALL cdn_read(
     &   eig_id,
     &   DIMENSION%nvd,input%jspins,mpi%irank,mpi%isize, !wannierspin instead of DIMENSION%jspd?&
     &   kptibz_b,jspin,DIMENSION%nbasfcn,
     &   noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,n_end,
     &   nbands_b2,eigg,zzMat)


        nslibd_b2 = 0

      IF(zzMat%l_real) THEN
         zMat_b2%z_r = 0.0
      ELSE
         zMat_b2%z_c = CMPLX(0.0,0.0)
      END IF

      eig_b2(:) = 0.

        do i = 1,nbands_b2
          if((eigg(i).ge.e1s.and.nslibd_b2.lt.numbands
     & .and.wann%l_bynumber)
     &.or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.wann%l_byenergy)
     &.or.(i.ge.wann%band_min(jspin_b)
     &                 .and.
     &       (i.le.wann%band_max(jspin_b))
     &                 .and.
     &       wann%l_byindex))then
            nslibd_b2 = nslibd_b2 + 1
            eig_b2(nslibd_b2) = eigg(i)
            we_b2(nslibd_b2) = we_b2(i)
            if(l_noco)then
               funbas =        nv_b2(1) + nlotot
               funbas = funbas+nv_b2(2) + nlotot
            else
               funbas = nv_b2(jspin_b) + nlotot
            endif
            IF (zzMat%l_real) THEN
               do j = 1,funbas
                  zMat_b2%data_r(j,nslibd_b2) = zzMat%data_r(j,i)
               enddo
            ELSE
               do j = 1,funbas
                  zMat_b2%data_c(j,nslibd_b2) = zzMat%data_c(j,i)
               enddo
            END IF
          endif
        enddo

c***********************************************************
c              Rotate the wavefunction of next neighbor.
c***********************************************************
      if (wann%l_bzsym .and. (oper_b2.ne.1)  ) then
!  TODO
!         call wann_kptsrotate(
!     >            natd,nlod,llod,
!     >            ntypd,nlo,llo,invsat,
!     >            l_noco,l_soc,
!     >            ntype,neq,nlotot,
!     >            kveclo_b2,jspin_b,
!     >            oper_b2,nop,mrot,nvd,
!     >            nv_b2,
!     >            shiftkpt(:,bpt(ikpt_b2,ikpt)),
!     >            tau,
!     x            bkpt_b2,k1_b2(:,:),
!     x            k2_b2(:,:),k3_b2(:,:),
!     x            zMat_b2,nsfactor_b2)
      else
         nsfactor_b2=cmplx(1.0,0.0)
      endif
        noccbd_b2 = nslibd_b2
      call cpu_time(t1)
      t_eig = t_eig + t1 - t0

      addnoco2 = 0
      if(l_noco.and.jspin_b.eq.2) addnoco2=nv_b2(1)+nlotot

      ! set up a(k+b2),b(k+b2),c(k+b2)
      allocate( acof_b2(noccbd_b2,0:lmd,natd),
     >          bcof_b2(noccbd_b2,0:lmd,natd),
     >          ccof_b2(-llod:llod,noccbd_b2,nlod,natd) )

      call cpu_time(t0)


!      ALLOCATE(lapw_b2%k1(SIZE(k1_b2,1),SIZE(k1_b2,2)))
!      ALLOCATE(lapw_b2%k2(SIZE(k1_b2,1),SIZE(k1_b2,2)))
!      ALLOCATE(lapw_b2%k3(SIZE(k1_b2,1),SIZE(k1_b2,2)))
!      lapw_b2%k1 = k1_b2
!      lapw_b2%k2 = k2_b2
!      lapw_b2%k3 = k3_b2
!      lapw_b2%nmat = nmat_b2
!      lapw_b2%nv = nv_b2
      ! I think the other variables of lapw are not needed here.
!TODO
!      CALL abcof(input,atoms,noccbd_b2,sym,cell,bkpt_b2,lapw_b2,
!     +           noccbd_b2,usdus,noco,jspin_b,kveclo_b2,oneD,
!     +           acof_b2,bcof_b2,ccof_b2,zMat_b2)

!      DEALLOCATE(lapw_b2%k1,lapw_b2%k2,lapw_b2%k3)

      call wann_abinv(atoms,
     X        acof_b2,bcof_b2,ccof_b2)
      call cpu_time(t1)
      t_abcof = t_abcof + t1 - t0


c**************************************c
c    calculate uHu matrix due to:      c
c             (i)    interstitial      c
c             (ii)   muffin tins       c
c             (iii)  vacuum            c 
c**************************************c
c      if(ANY(gb(3:3,ikpt_b,ikpt).ne.0) .or.
c     >   ANY(gb(3:3,ikpt_b2,ikpt).ne.0)) then
c       write(*,*)ikpt,ikpt_b,ikpt_b2
c       write(*,*)gb(:,ikpt_b,ikpt)
c       write(*,*)gb(:,ikpt_b2,ikpt)
c       uHu(:,:,ikpt_b2,ikpt_b,i_rec) = cmplx(0.,0.)
c       goto 848
c      endif

      ! MT (SPH, NON)
      if(.not.(l_skip_sph.and.l_skip_non)) then
         call cpu_time(t0)
         call wann_uHu_sph(
     >        cmplx(1.,0.),nbnd,llod,nslibd_b,nslibd_b2,nlod,natd,ntypd,
     >        lmd,jmtd,taual,nop,lmax,ntype,neq,nlo,llo,
     >        acof_b,bcof_b,ccof_b,bkpt_b2,
     >        acof_b2,bcof_b2,ccof_b2,bkpt_b,bkpt,
     >        gb(:,ikpt_b,ikpt),gb(:,ikpt_b2,ikpt),
     <        tuu(:,:,:,:,tspin2),tud(:,:,:,:,tspin2),
     >        tdu(:,:,:,:,tspin2),tdd(:,:,:,:,tspin2),
     >        tuulo(:,:,:,:,:,tspin2),tulou(:,:,:,:,:,tspin2),
     >        tdulo(:,:,:,:,:,tspin2),tulod(:,:,:,:,:,tspin2),
     >        tuloulo(:,:,:,:,:,:,tspin2),
     >        kdiff,kdiff,nntot,nntot,
     >        uHu(:,:,ikpt_b2,ikpt_b,i_rec))
         call cpu_time(t1)
         t_sph = t_sph + t1 - t0
      endif
c      if(ikpt.eq.1) then
c       write(*,*)'SPH,NON',ikpt_b2,ikpt_b
c       write(*,*)uHu(1,3,ikpt_b2,ikpt_b,i_rec)
c       write(*,*)uHu(3,1,ikpt_b2,ikpt_b,i_rec)
c       uHuold(:,:) = uHu(:,:,ikpt_b2,ikpt_b,i_rec)
c      endif

      ! MT (SOC)
      if(l_soc.and.(.not.l_skip_soc)) then
         call cpu_time(t0)
         call wann_uHu_sph(
     >        cmplx(1.,0.),nbnd,llod,nslibd_b,nslibd_b2,nlod,natd,ntypd,
     >        lmd,jmtd,taual,nop,lmax,ntype,neq,nlo,llo,
     >        acof_b,bcof_b,ccof_b,bkpt_b2,
     >        acof_b2,bcof_b2,ccof_b2,bkpt_b,bkpt,
     >        gb(:,ikpt_b,ikpt),gb(:,ikpt_b2,ikpt),
     <        tuu_soc,tud_soc,
     >        tdu_soc,tdd_soc,
     >        tuulo_soc,tulou_soc,
     >        tdulo_soc,tulod_soc,
     >        tuloulo_soc,
     >        kdiff,kdiff,nntot,nntot,
     >        uHu(:,:,ikpt_b2,ikpt_b,i_rec))
         call cpu_time(t1)
         t_sph = t_sph + t1 - t0
c      if(ikpt.eq.1) then
c       write(*,*)'SOC',ikpt_b2,ikpt_b
c       write(*,*)uHu(1,3,ikpt_b2,ikpt_b,i_rec)-uHuold(1,3)
c       write(*,*)uHu(3,1,ikpt_b2,ikpt_b,i_rec)-uHuold(3,1)
c       uHuold(:,:) = uHu(:,:,ikpt_b2,ikpt_b,i_rec)
c      endif
      endif


      ! In collinear calculation, the potential in interstital and
      ! vacuum region is always diagonal with respect to the spin
      IF( l_noco .OR. ((.NOT. l_noco) .AND. (doublespin.LT.3)) ) THEN

         jspin3 = doublespin
         if(jspins.eq.1) jspin3=1

         jspin4 = jspin
         jspin4_b=jspin_b
         if(jspins.eq.1) then
          jspin4=1
          jspin4_b=1
         endif

         sign2 = 1
         if(doublespin.EQ.4) sign2=-1

         ! INT
         if(.not.l_skip_int) then
            call cpu_time(t0)
            call wann_uHu_int(cmplx(1.,0.),nvd,k1d,k2d,k3d,n3d,
     >            nv_b(jspin),nv_b2(jspin_b),nbnd,neigd,
     >            nslibd_b,nslibd_b2,nbasfcn,addnoco,addnoco2,
     >            k1_b(:,jspin), k2_b(:,jspin), k3_b(:,jspin),
     >            gb(:,ikpt_b,ikpt),
     >            k1_b2(:,jspin_b),k2_b2(:,jspin_b),k3_b2(:,jspin_b),
     >            gb(:,ikpt_b2,ikpt),
     >            bkpt,bbmat,vpw(:,jspin3),zMat_b,zMat_b2,rgphs,
     >            ustep,ig,jspin.eq.jspin_b,sign2,
     >            uHu(:,:,ikpt_b2,ikpt_b,i_rec))
            call cpu_time(t1)
            t_int = t_int + t1 - t0
         endif
c      if(ikpt.eq.1) then
c       write(*,*)'INT',ikpt_b2,ikpt_b
c       write(*,*)uHu(1,3,ikpt_b2,ikpt_b,i_rec)-uHuold(1,3)
c       write(*,*)uHu(3,1,ikpt_b2,ikpt_b,i_rec)-uHuold(3,1)
c       uHuold(:,:) = uHu(:,:,ikpt_b2,ikpt_b,i_rec)
c      endif


         ! VAC
         if ((.not.l_skip_vac) .and. film .and. (.not.odi%d1)) then

            call cpu_time(t0)
            call wann_uHu_vac(
     >            cmplx(1.,0.),l_noco,l_soc,zrfs,jspins,nlotot,qpt_i,
     >            nbnd,z1,nmzxyd,nmzd,n2d,nv2d,k1d,k2d,k3d,n3d,nvac,ig,
     >            rgphs,nmzxy,nmz,delz,ig2,nq2,kv2,area,bmat,bbmat,
     >            evac(:,jspin4),evac(:,jspin4_b),bkpt_b,bkpt_b2,
     >            vzxy(:,:,:,jspin3),vz,nslibd_b,nslibd_b2,
     >            jspin,jspin_b,doublespin,k1_b,k2_b,k3_b,
     >            k1_b2,k2_b2,k3_b2,wannierspin,nvd,nbasfcn,neigd,
     >            zMat_b,zMat_b2,nv_b,nv_b2,omtil,gb(:,ikpt_b,ikpt),
     >            gb(:,ikpt_b2,ikpt),sign2,
     >            uHu(:,:,ikpt_b2,ikpt_b,i_rec))
            call cpu_time(t1)
            t_vac = t_vac + t1 - t0

         elseif ((.not.l_skip_vac) .and. odi%d1) then

            call cpu_time(t0)
            call wann_uHu_od_vac(
     >            DIMENSION,oneD,vacuum,stars,cell,
     >            cmplx(1.,0.),l_noco,l_soc,jspins,nlotot,nbnd,z1,
     >            nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,ig,nmzxy,nmz,
     >            delz,ig2,bbmat,evac(1,jspin4),evac(1,jspin4_b),bkpt_b,
     >            bkpt_b2,odi,vzxy(:,:,:,jspin3),vz,nslibd_b,nslibd_b2,
     >            jspin,jspin_b,doublespin,k1_b,k2_b,k3_b,k1_b2,k2_b2,
     >          k3_b2,wannierspin,nvd,area,nbasfcn,neigd,zMat_b,zMat_b2,
     >            nv_b,nv_b2,sk2,phi2,omtil,gb(:,ikpt_b,ikpt),
     >            gb(:,ikpt_b2,ikpt),qpt_i,sign2,
     >            uHu(:,:,ikpt_b2,ikpt_b,i_rec))
            call cpu_time(t1)
            t_vac = t_vac + t1 - t0

         endif
c      if(ikpt.eq.1) then
c       write(*,*)'VAC',ikpt_b2,ikpt_b
c       write(*,*)uHu(1,3,ikpt_b2,ikpt_b,i_rec)-uHuold(1,3)
c       write(*,*)uHu(3,1,ikpt_b2,ikpt_b,i_rec)-uHuold(3,1)
c       uHuold(:,:) = uHu(:,:,ikpt_b2,ikpt_b,i_rec)
c      endif
      ELSE
       IF(ikpt.EQ.1 .AND. ikpt_b.EQ.1 .AND. ikpt_b2.EQ.1) THEN
          if(l_p0) write(*,*)'skip INT and VAC'
       ENDIF
      ENDIF


848     continue

        deallocate (acof_b2,bcof_b2,ccof_b2)
25      continue ! end of loop by the nearest k-neighbors

        deallocate (acof_b,bcof_b,ccof_b)
15      continue ! end of loop by the nearest k-neighbors
        deallocate ( we_b,we_b2 )
      deallocate ( eigg )



      IF(l_real) THEN
         DEALLOCATE ( zzMat%z_r )
      ELSE
         DEALLOCATE ( zzMat%z_c )
      END IF



      IF(zzMat%l_real) THEN
         DEALLOCATE (zMat_b%z_r)
         DEALLOCATE (zMat_b2%z_r)
      ELSE
         DEALLOCATE (zMat_b%z_c)
         DEALLOCATE (zMat_b2%z_c)
      END IF



      endif   ! loop by processors






10    continue ! end of cycle by the k-points

#ifdef CPP_MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
  5   continue

      l_symcheck = (.not.l_nocosoc).or.(doublespin.eq.doublespin_max)
      jspin4=jspin
      if(l_nocosoc) jspin4=1
      if(l_p0) write(*,*)'write uHu file'
      uHu = uHu * htr2ev
      call wann_write_uHu(jspin4,l_p0,fullnkpts,nntot,nntot,wann,
     >          nbnd,bpt,gb,isize,irank,fending,'',uHu,
     >          counts(irank),counts,displs,isize,
     >          wann%l_unformatted,.true.,l_symcheck)
      uHu = uHu / htr2ev
      if(.not.l_nocosoc) deallocate(uHu)

      if(.not.l_noco)nrec=nrec+nkpts

110   continue ! end of cycle by spins

#ifdef CPP_MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      ! close eig files
      IF (l_gwf) THEN

!         CALL close_eig(eig_id)

!         CLOSE(66)
         DO iqpt_b=1,nntot_q
!            CALL close_eig(innerEig_idList(iqpt_b))
         ENDDO
      ENDIF


c      if(l_p0) write(*,*)'write uHu file'
c      uHu = htr2ev * uHu
c      call wann_write_uHu(1,l_p0,fullnkpts,nntot,nntot,wann,
c     >          nbnd,bpt,gb,isize,irank,fending,'',uHu,
c     >          counts(irank),counts,displs,isize,
c     >          wann%l_unformatted,.true.)
      if(allocated(uHu)) deallocate(uHu)
c      if(allocated(uHuold)) deallocate(uHuold)

      DEALLOCATE(innerEig_idList)

314   continue ! iqpt, q-points
c************************************************c
c               END Q LOOP                       c
c************************************************c


      deallocate ( kveclo,nv,k1,k2,k3 )
      deallocate(flo)
      deallocate ( ff,gg)
      if (allocated(nv_b))deallocate(kveclo_b,nv_b,k1_b,k2_b,k3_b)
      if (allocated(nv_b2))deallocate(kveclo_b2,nv_b2,k1_b2,k2_b2,k3_b2)
      if (wann%l_bzsym)deallocate(irreduc,mapkoper,shiftkpt)
      if (wann%l_bzsym.AND.l_gwf)deallocate(irreduc_q,mapqoper,shiftqpt)
      if (allocated(pair_to_do)) deallocate(pair_to_do,maptopair)
      if (allocated(pair_to_do_q)) deallocate(pair_to_do_q,maptopair_q)
      if (allocated(kdiff)) deallocate ( kdiff )
      if (allocated(qdiff)) deallocate(qdiff,zero_qdiff)
      if(allocated(vpw)) deallocate(vpw)
      if(allocated(vzxy)) deallocate(vzxy)
      if(allocated(vz)) deallocate(vz)
      deallocate(vr)
      deallocate(tdd,tdu,tud,tuu)
      deallocate(tdulo,tuulo)
      deallocate(tulou,tulod)
      deallocate(tuloulo)
      deallocate(tdd_soc,tdu_soc,tud_soc,tuu_soc)
      deallocate(tdulo_soc,tuulo_soc)
      deallocate(tulou_soc,tulod_soc)
      deallocate(tuloulo_soc)
      deallocate(counts,displs)   

9110  continue

      if(l_sgwf .or. l_socgwf) gb_q = gb_q/2

      if(l_p0.and.l_gwf) then
         call wann_uHu_commat(
     >              fullnkpts,nntot,bpt,fullnqpts,nntot_q,bpt_q,
     >              gb,gb_q,aux_latt_const,wann%l_unformatted,l_dim,
     >              nparampts,param_vec/2.0)
      endif

      if(allocated(gb)) deallocate(gb,bpt)
      if(allocated(gb_q)) deallocate(gb_q,bpt_q)

#ifdef CPP_MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      call cpu_time(t1)
      t_total = t1-t00
      if(l_p0) then
       write(*,900)'t_init  =',t_init
       write(*,900)'t_myTlmplm=',t_myTlmplm
       write(*,900)'t_eig   =',t_eig   
       write(*,900)'t_abcof =',t_abcof 
       write(*,900)'t_int   =',t_int   
       write(*,900)'t_sph   =',t_sph   
       write(*,900)'t_vac   =',t_vac   
       write(*,900)'t_total =',t_total
      endif

900   FORMAT(a,f14.4)

      END SUBROUTINE wann_uHu

 
      END MODULE m_wann_uHu


