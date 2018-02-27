!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wannier
      use m_juDFT
      CONTAINS
      SUBROUTINE wannier(
     >      DIMENSION,mpi,input,kpts,sym,atoms,stars,vacuum,sphhar,oneD,
     >      wann,noco,cell,enpara,banddos,sliceplot,vTot,results,
     >      eig_idList,l_real,nkpt)
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
c     In general, the number of bands for each k-poin t is
c     different: N(k), and also differs from the number of bands
c     we are interested in: N (for instance 5 d-bands of Cu among
c     the 6 s- and d-bands). While matrices Mmn(k) are square
c     for each k-point, matrices Mmn(k,b) can be made so after
c     defining the maximum number of bands max(N(k)).
c     The matrix Amn is non-diagonal by default (N(k)*N).
cccccc ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Total number of wannier functions: nwfs
c
c     sliceplot%e1s,sliceplot%e2s: lower and upper boundaries of the energy window:
c     Needed for sorting by number and sorting by energy.
c     Not needed for sorting by index.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Extension to case of higher-dimensional Wannier functions        
c     according to the formalism in PRB 91, 184413 (2015)
c     Jan-Philipp Hanke                                         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use m_types
      use m_wann_mmnk_symm
      use m_wann_rw_eig
      use m_abcof
      use m_radfun
      use m_radflo
      use m_cdnread
      use m_constants
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
      use m_wann_gwf_tools, only : get_index_kq, gwf_plottemplate
      use m_wann_gwf_commat
      use m_wann_gwf_anglmom
      use m_wann_write_mmnk2
      use m_wann_uHu
      use m_wann_uHu_dmi
      USE m_eig66_io

      IMPLICIT NONE
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer ierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)
#endif

      TYPE(t_dimension), INTENT(IN) :: DIMENSION
      TYPE(t_mpi),       INTENT(IN) :: mpi
      TYPE(t_input),     INTENT(IN) :: input
      TYPE(t_kpts),      INTENT(IN) :: kpts
      TYPE(t_sym),       INTENT(IN) :: sym
      TYPE(t_atoms),     INTENT(IN) :: atoms
      TYPE(t_stars),     INTENT(IN) :: stars
      TYPE(t_vacuum),    INTENT(IN) :: vacuum
      TYPE(t_sphhar),    INTENT(IN) :: sphhar
      TYPE(t_oneD),      INTENT(IN) :: oneD
      TYPE(t_noco),      INTENT(IN) :: noco
      TYPE(t_cell),      INTENT(IN) :: cell
      TYPE(t_enpara),    INTENT(IN) :: enpara
      TYPE(t_banddos),   INTENT(IN) :: banddos
      TYPE(t_sliceplot), INTENT(IN) :: sliceplot
      TYPE(t_potden),    INTENT(IN) :: vTot
      TYPE(t_results),   INTENT(IN) :: results
      TYPE(t_wann),      INTENT(INOUT) :: wann

      logical, intent (in) :: l_real
      integer, intent (in) :: nkpt
      INTEGER, INTENT (IN) :: eig_idList(wann%nparampts)

cccccccccccccccccc   local variables   cccccccccccccccccccc
      integer :: lmd,nlotot,n,nmat,iter,ikpt,ikpt_b,nmat_b, pc
      integer :: addnoco,funbas,loplod,addnoco2,igvm2,eig_id
      integer :: noccbd,noccbd_b,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
      integer :: jsp_start,jsp_end,nrec,nrec1,nrec_b,nbands,nbands_b
      integer :: nodeu,noded,n_size,na,n_rank,nbnd,numbands
      integer :: i1,i2,i3,in,lda
      integer :: n_bands(0:DIMENSION%neigd),nslibd,nslibd_b
      character(len=8) :: dop,iop,name(10)
      real    :: bkpt(3),bkpt_b(3),sfp,tpi,wronk,wk,wk_b,phase
      complex :: c_phase
      real    :: eig(DIMENSION%neigd),eig_b(DIMENSION%neigd),cp_time(9)
      real    :: efermi
      logical :: l_p0,l_bkpts,l_proj,l_amn,l_mmn
!!! energy window boundaries
      integer, allocatable :: kveclo(:),nv(:)
      integer, allocatable :: kveclo_b(:),nv_b(:)
      integer, allocatable :: k1(:,:),k2(:,:),k3(:,:)
      integer, allocatable :: k1_b(:,:),k2_b(:,:),k3_b(:,:)
      INTEGER, ALLOCATABLE :: innerEig_idList(:)
      real,    allocatable :: we(:),we_b(:)

      real,    allocatable :: eigg(:)
      real kpoints(nkpt)
!!! a and b coeff. constructed for each k-point
      complex, allocatable :: acof(:,:,:),acof_b(:,:,:)
      complex, allocatable :: bcof(:,:,:),bcof_b(:,:,:)
      complex, allocatable :: ccof(:,:,:,:),ccof_b(:,:,:,:)
!!! the parameters for the number of wfs
      integer :: nwfs
!!! the potential in the spheres and the vacuum
      real, allocatable :: vr(:,:,:),vz(:,:,:)
!!! auxiliary potentials
      complex, allocatable :: vpw(:,:)
!!! bkpts data
      integer nntot,ikpt_help
      integer, allocatable :: gb(:,:,:),bpt(:,:)
!!! radial wavefunctions in the muffin-tins and more ...
      real,    allocatable :: flo(:,:,:,:,:),vso(:,:,:)
      real,    allocatable :: ff(:,:,:,:,:),gg(:,:,:,:,:)

      real     :: uuilon(atoms%nlod,atoms%ntype)
      real     :: duilon(atoms%nlod,atoms%ntype)
      real     :: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
!!! energy parameters
      real    :: ello(atoms%nlod,atoms%ntype,max(2,DIMENSION%jspd))
      real    :: evac(2,max(2,DIMENSION%jspd))
      real    :: epar(0:atoms%lmaxd,atoms%ntype,max(2,DIMENSION%jspd))
      real    :: evdu(2,max(DIMENSION%jspd,2))
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
      integer :: ngopr1(atoms%nat)
      integer,allocatable::maptopair(:,:,:)
      integer :: wannierspin,jspin2,jspin7,jspin2_b
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


c---->gwf

      ! FURTHER VARIABLES
      real :: qpt_i(3),qptb_i(3)
      real :: alph_i(atoms%ntype),alphb_i(atoms%ntype)
      real :: beta_i(atoms%ntype),betab_i(atoms%ntype)
      real :: theta_i, thetab_i, phi_i, phib_i
      real :: dalph,db1,db2,coph,siph
      real :: zero_taual(3,atoms%nat),bqpt(3),wqb
      real :: tt1,tt2,tt3,tt4,tt5,tt6,t_it,t_kov,t_qov
      real :: eig_qb(DIMENSION%neigd)

      real,allocatable :: qdiff(:,:), we_qb(:)                
      real,allocatable :: energies(:,:,:)  
      real,allocatable :: zero_qdiff(:,:)


      integer,allocatable :: irreduc_q(:),mapqoper(:)        
      integer,allocatable :: shiftqpt(:,:),pair_to_do_q(:,:)  
      integer,allocatable :: maptopair_q(:,:,:)              
      integer,allocatable :: gb_q(:,:,:),bpt_q(:,:)         
      integer,allocatable :: k1_qb(:,:),k2_qb(:,:),k3_qb(:,:)
      integer,allocatable :: kveclo_qb(:),nv_qb(:)

      integer :: nntot_q = 1                               
      integer :: fullnqpts = 1                                
      integer :: funit_start = 5000
      integer :: qptibz, qptibz_b, oper_q, oper_qb
      integer :: qpt,iqpt_help, iqpt, iqpt_b
      integer :: nbands_qb, nmat_qb, nslibd_qb, noccbd_qb
      integer :: sign_q = 1,band_help
      integer :: doublespin,jspin_b,jspin3,jspin4,jspin5
      integer :: doublespin_max,nrec5
      integer :: count_i,count_j
      integer :: aoff,d1,d10,d100
      integer :: n1,n2,ii,jj

      complex :: interchi,vacchi,amnchi
      complex :: phasfac,phasfac2,cmplx_1                                 

      complex,allocatable :: chi(:)
      complex,allocatable :: acof_qb(:,:,:)                 
      complex,allocatable :: bcof_qb(:,:,:)                   
      complex,allocatable :: ccof_qb(:,:,:,:)                 
      complex,allocatable :: mmnk_q(:,:,:,:)                 
      complex,allocatable :: m_int(:,:,:,:)
      complex,allocatable :: m_sph(:,:,:,:)
      complex,allocatable :: m_vac(:,:,:,:)
      complex,allocatable :: ujug_q(:,:,:,:),ujdg_q(:,:,:,:) 
      complex,allocatable :: djug_q(:,:,:,:),djdg_q(:,:,:,:) 
      complex,allocatable :: ujulog_q(:,:,:,:,:)              
      complex,allocatable :: djulog_q(:,:,:,:,:)             
      complex,allocatable :: ulojug_q(:,:,:,:,:)             
      complex,allocatable :: ulojdg_q(:,:,:,:,:)             
      complex,allocatable :: ulojulog_q(:,:,:,:,:,:)         

      character(len=12) tmp_filename,fending
      character(len=30) fname,fstart,eigfile

      logical :: l_bqpts,l_gwf,l_nochi

      TYPE(t_usdus) :: usdus
      TYPE(t_zmat)  :: zMat, zzMat, zMat_b, zMat_qb
      TYPE(t_lapw)  :: lapw, lapw_b, lapw_qb
      TYPE(t_wann)  :: wannTemp

      eig_id = eig_idList(1)

c----<gwf


      aoff = iachar('1')-1
      d1  = mod(mpi%irank,10)
      IF (mpi%irank < 100) THEN
        d10 = int( (mpi%irank + 0.5)/10 )
        fstart = 'eig'//achar(d10+aoff)//achar(d1+aoff)
      ELSE
        d10 = mod((mpi%irank-d1)/10,10)
        d100 = (mpi%irank-10*d10-d1)/100
        IF ( d100.GE.10 ) d100 = d100 + iachar('7')
        fstart =
     +  'eig'//achar(d100+aoff)//achar(d10+aoff)//achar(d1+aoff)
      ENDIF


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
      zero_taual = 0.0

      ci=cmplx(0.0,1.0)
      cmplx_1 = cmplx(1.0,0.0)

      call cpu_time(time_total)

      l_p0 = .false.
      if (mpi%irank.eq.0) l_p0 = .true.

      sfp = 2* sqrt( pimach() )
      tpi = 2* pimach()
      lmd = atoms%lmaxd*(atoms%lmaxd+2)

!!!   should be changed in case the windows are really used
      nkpts = nkpt

      ! do we have to construct GWF ?
      l_gwf = .false.
      l_gwf = wann%l_sgwf.or.wann%l_socgwf 
      
      l_nochi = .false.
      inquire(file='l_nochi',exist=l_nochi)
      if(l_gwf.and.l_p0) write(*,*)'disable chi trafo: ',l_nochi

      if(l_gwf.and.l_p0) call gwf_plottemplate()
      allocate( chi(atoms%ntype) )

c-----read the input file to determine what to do
      wann%atomlist_num=atoms%nat
      wann%oc_num_orbs=atoms%nat

      CALL wann_read_inp(input,l_p0,wann)


c-----input file for orbital decomposition
      if(wann%l_orbcomp.or.wann%l_orbcomprs)then
         inquire(file='orbcomp_inp',exist=l_orbcompinp)
         if(l_orbcompinp)then
            open(159,file='orbcomp_inp')
            read(159,*)wann%oc_num_orbs,wann%l_oc_f
            allocate(wann%oc_orbs(wann%oc_num_orbs))
            do n=1,wann%oc_num_orbs
               read(159,*)wann%oc_orbs(n)
            enddo
            close(159)
         else !default is all atoms including f
            wann%oc_num_orbs=atoms%nat
            wann%l_oc_f=.true.
            allocate(wann%oc_orbs(wann%oc_num_orbs))
            do n=1,wann%oc_num_orbs
               wann%oc_orbs(n)=n
            enddo          
         endif
      endif

      if(wann%l_updown)then            
         call wann_updown(
     >      mpi,input,sym,atoms,stars,vacuum,sphhar,oneD,noco,cell,vTot,
     >      enpara,eig_idList(1),l_real,
     >      mpi%mpi_comm,atoms%l_dulo,noco%l_noco,noco%l_ss,
     >      atoms%lmaxd,atoms%ntype,DIMENSION%neigd,atoms%nat,sym%nop,
     >      DIMENSION%nvd,DIMENSION%jspd,DIMENSION%nbasfcn,atoms%llod,
     >      atoms%nlod,atoms%ntype,cell%omtil,atoms%nlo,atoms%llo,
     >      atoms%lapw_l,sym%invtab,sym%mrot,atoms%ngopr,atoms%neq,
     >      atoms%lmax,atoms%invsat,sym%invsatnr,nkpt,atoms%taual,
     >      atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,
     >      noco%beta,noco%qss,                    ! TODO: adapt if needed
     >      stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,
     >      stars%ng3,
     >      vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,sphhar%nlhd,stars%ng3,
     >      vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,
     >      atoms%jri,sphhar%ntypsd,atoms%ntypsy,input%jspins,nkpt,
     >      atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,
     >      atoms%ulo_der,stars%ustep,stars%ig,stars%mx1,
     >      stars%mx2,stars%mx3,stars%rgphs,
     >      sliceplot%slice,sliceplot%kk,sliceplot%nnne,cell%z1,
     >      DIMENSION%nv2d,vacuum%nmzxy,vacuum%nmz,vacuum%delz,
     >      stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,
     >      cell%volint,sym%symor,atoms%pos,results%ef,noco%l_soc,
     >      sphhar%memd,atoms%lnonsph,sphhar%clnu,DIMENSION%lmplmd,
     >      sphhar%mlh,sphhar%nmem,sphhar%llh,atoms%lo1l,
     >      noco%theta,noco%phi,noco%soc_opt)

         DO pc = 1, wann%nparampts
            CALL close_eig(eig_idList(pc))
         END DO

         CALL juDFT_end("updown done",mpi%irank)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  modern theory of orbital magnetization from Wannier functions
c  Jan-Philipp Hanke
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(wann%l_matrixuHu)then
         wannTemp = wann
         call wann_uHu(
     >      DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,
     >      banddos,oneD,noco,cell,vTot,wannTemp,eig_idList,
     >      l_real,atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,
     >      atoms%ntype,DIMENSION%neigd,atoms%nat,sym%nop,DIMENSION%nvd,
     >      DIMENSION%jspd,DIMENSION%nbasfcn,atoms%llod,atoms%nlod,
     >      atoms%ntype,cell%omtil,atoms%nlo,atoms%llo,
     >      atoms%lapw_l,sym%invtab,sym%mrot,atoms%ngopr,atoms%neq,
     >      atoms%lmax,atoms%invsat,sym%invsatnr,nkpt,atoms%taual,
     >      atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,
     >      noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,
     >      mpi%irank,
     >      mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,
     >      sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,
     >      input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,atoms%ntypsy,
     >      input%jspins,nkpt,atoms%dx,stars%ng2,atoms%rmsh,
     >      sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,stars%ustep,
     >      stars%ig,stars%mx1,stars%mx2,stars%mx3,
     >      stars%rgphs,sliceplot%slice,
     >      sliceplot%kk,sliceplot%nnne,cell%z1,DIMENSION%nv2d,
     >      vacuum%nmzxy,vacuum%nmz,vacuum%delz,sym%zrfs,stars%ig2,
     >      cell%area,sym%tau,atoms%zatom,stars%ng2,stars%kv2,sym%nop2,
     >      cell%volint,sym%symor,atoms%pos,results%ef,noco%l_soc,
     >      sphhar%memd,atoms%lnonsph,sphhar%clnu,DIMENSION%lmplmd,
     >      sphhar%mlh,sphhar%nmem,sphhar%llh,atoms%lo1l,
     >      noco%theta,noco%phi,noco%soc_opt,
     >      wann%l_ms,wann%l_sgwf,wann%l_socgwf,wann%aux_latt_const,
     >      wann%param_file,wann%param_vec,wann%nparampts,
     >      wann%param_alpha,wann%l_dim)

         DO pc = 1, wann%nparampts
            CALL close_eig(eig_idList(pc))
         END DO

        CALL juDFT_end("wann_uHu done",mpi%irank)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  modern theory of DMI from higher-dimensional Wannier functions
c  Jan-Philipp Hanke
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(wann%l_matrixuHu_dmi)then
         wannTemp = wann
         call wann_uHu_dmi(
     >      DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,
     >      banddos,oneD,noco,cell,vTot,wannTemp,eig_idList,
     >      l_real,atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,
     >      atoms%ntype,DIMENSION%neigd,atoms%nat,sym%nop,DIMENSION%nvd,
     >      DIMENSION%jspd,DIMENSION%nbasfcn,atoms%llod,atoms%nlod,
     >      atoms%ntype,cell%omtil,atoms%nlo,atoms%llo,
     >      atoms%lapw_l,sym%invtab,sym%mrot,atoms%ngopr,atoms%neq,
     >      atoms%lmax,atoms%invsat,sym%invsatnr,nkpt,atoms%taual,
     >      atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,
     >      noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,
     >      mpi%irank,
     >      mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,
     >      sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,
     >      input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,atoms%ntypsy,
     >      input%jspins,nkpt,atoms%dx,stars%ng2,atoms%rmsh,
     >      sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,stars%ustep,
     >      stars%ig,stars%mx1,stars%mx2,stars%mx3,
     >      stars%rgphs,sliceplot%slice,
     >      sliceplot%kk,sliceplot%nnne,cell%z1,DIMENSION%nv2d,
     >      vacuum%nmzxy,vacuum%nmz,vacuum%delz,sym%zrfs,stars%ig2,
     >      cell%area,sym%tau,atoms%zatom,stars%ng2,stars%kv2,sym%nop2,
     >      cell%volint,sym%symor,atoms%pos,results%ef,noco%l_soc,
     >      sphhar%memd,atoms%lnonsph,sphhar%clnu,DIMENSION%lmplmd,
     >      sphhar%mlh,sphhar%nmem,sphhar%llh,atoms%lo1l,
     >      noco%theta,noco%phi,noco%soc_opt,
     >      wann%l_ms,wann%l_sgwf,wann%l_socgwf,wann%aux_latt_const,
     >      wann%param_file,wann%param_vec,wann%nparampts,
     >      wann%param_alpha,wann%l_dim,l_nochi)

         DO pc = 1, wann%nparampts
            CALL close_eig(eig_idList(pc))
         END DO

         CALL juDFT_end("wann_uHu dmi done",mpi%irank)
      endif

      if(wann%l_byenergy.and.wann%l_byindex) CALL juDFT_error
     +     ("byenergy.and.byindex",calledby ="wannier")
      if(wann%l_byenergy.and.wann%l_bynumber) CALL juDFT_error
     +     ("byenergy.and.bynumber",calledby ="wannier")
      if(wann%l_bynumber.and.wann%l_byindex) CALL juDFT_error
     +     ("bynumber.and.byindex",calledby ="wannier")
      if(.not.(wann%l_bynumber.or.wann%l_byindex.or.wann%l_byenergy))
     &     CALL juDFT_error("no rule to sort bands",calledby ="wannier")


      efermi=results%ef
      if(.not.wann%l_fermi)efermi=0.0

#ifdef CPP_MPI
         call MPI_BARRIER(mpi%mpi_comm,ierr)
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
         allocate(irreduc(fullnkpts),mapkoper(fullnkpts))
         allocate(shiftkpt(3,fullnkpts))
      endif


      if(l_gwf) fullnqpts = wann%nparampts


      nrec = 0
      if(l_p0)then
      write (*,*) 'fermi energy:',efermi
      write (*,*) 'emin,emax=',sliceplot%e1s,sliceplot%e2s
      write (*,*) 'nbasfcn =',DIMENSION%nbasfcn
      endif
      nlotot = 0
      do n = 1, atoms%ntype
        do l = 1,atoms%nlo(n)
          nlotot = nlotot + atoms%neq(n) * ( 2*atoms%llo(l,n) + 1 )
        enddo
      enddo

      if((.not.wann%l_matrixmmn).and.(.not.wann%l_wann_plot).and.
     &   (.not.wann%l_matrixamn).and.(.not.wann%l_projmethod).and.
     &   (.not.wann%l_bestproj).and.(.not.wann%l_nabla).and.
     &   (.not.wann%l_mmn0).and.(.not.wann%l_surfcurr).and.
     &   (.not.wann%l_offdiposop).and.(.not.wann%l_anglmom).and.
     &   (.not.wann%l_orbcomp).and.(.not.wann%l_perturb) .and.
     &   (.not.wann%l_finishgwf) ) goto 1911

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
         if (bpt(nn,ikpt)>fullnkpts) CALL juDFT_error("bpt.gt.fullnkpts"
     +        ,calledby ="wannier")
        enddo
       enddo
       close (202)
        allocate(kdiff(3,nntot))
      endif

c**********************************************************
ccccccccccccccc   read in the bqpts file  ccccccccccccccccc         
c**********************************************************
      if ((wann%l_matrixmmn).AND.(l_gwf.or.wann%l_ms)) then
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
      if(wann%l_sgwf) gb_q = 2*gb_q    
      if(wann%l_socgwf) gb_q = 2*gb_q 

      if(wann%l_finishgwf) goto 9110
c********************************************************
c      find symmetry-related elements in mmkb
c********************************************************
      if(wann%l_matrixmmn)then
         call cpu_time(delta)
         allocate(maptopair(3,fullnkpts,nntot))
         allocate(pair_to_do(fullnkpts,nntot))
         call wann_mmnk_symm(input,kpts,
     >     fullnkpts,nntot,bpt,gb,wann%l_bzsym,
     >     irreduc,mapkoper,l_p0,input%film,sym%nop,sym%invtab,sym%mrot,
     >     oneD%odi%d1,sym%tau,
     <     pair_to_do,maptopair,kdiff,.false.,wann%param_file)
         call cpu_time(delta1)
         time_symm=delta1-delta
      endif

      ! do the same for q-points to construct GWFs
      if(wann%l_matrixmmn.AND.l_gwf)then 
         allocate(maptopair_q(3,fullnqpts,nntot_q))
         allocate(pair_to_do_q(fullnqpts,nntot_q))
         call wann_mmnk_symm(input,kpts,
     >     fullnqpts,nntot_q,bpt_q,gb_q,wann%l_bzsym,
     >     irreduc_q,mapqoper,l_p0,.false.,1,sym%invtab(1),
     >     sym%mrot(:,:,1),.false.,sym%tau,
     <     pair_to_do_q,maptopair_q,qdiff,.true.,wann%param_file)
      endif


c*********************************************************
cccccccccccccccc   initialize the potential   cccccccccccc
c*********************************************************

      allocate ( vz(vacuum%nmzd,2,4) )
      allocate ( vr(atoms%jmtd,atoms%ntype,DIMENSION%jspd) )
      allocate ( vso(atoms%jmtd,atoms%nat,2) )

      vz = 0.0
      vz(:,:,:SIZE(vTot%vacz,3)) = vTot%vacz(:,:,:)

      do jspin = 1,input%jspins
        do n = 1, atoms%ntype
          do j = 1,atoms%jri(n)
            vr(j,n,jspin) = vTot%mt(j,0,n,jspin)
          enddo
        enddo
      enddo

      if(wann%l_soctomom)then
        CALL vsoc(input,atoms,vr,epar,.TRUE., vso)
      endif

      if(noco%l_noco.and.input%film)then
         npotmatfile=25
         allocate(vpw(stars%ng3,1))
         if(.not.oneD%odi%d1)
     &    allocate( vxy(vacuum%nmzxyd,stars%ng2-1,2) )

         OPEN (npotmatfile,FILE='potmat',FORM='unformatted',
     +                  STATUS='old')
         READ (npotmatfile) (vpw(ig3,1),ig3=1,stars%ng3)
         READ (npotmatfile) (vpw(ig3,1),ig3=1,stars%ng3)
         READ (npotmatfile) (vpw(ig3,1),ig3=1,stars%ng3)
         maxvac=2
         if(oneD%odi%d1)maxvac=1
         DO ivac = 1,maxvac
c--->       if the two vacuua are equivalent, the potential file has to
c--->       be backspaced, because the potential is the same at both
c--->       surfaces of the film
            IF ((ivac.EQ.2) .AND. (vacuum%nvac.EQ.1)) THEN
               DO irec = 1,4
                  BACKSPACE (npotmatfile)
               ENDDO
            ENDIF
c--->       load the non-warping part of the potential
            READ (npotmatfile)
     +         ((vz(imz,ivac,ipot),imz=1,vacuum%nmzd),ipot=1,4)

            if(.not.oneD%odi%d1)then
            DO ipot = 1,3
               READ (npotmatfile)((vxy(imz,igvm2,ivac),
     +            imz=1,vacuum%nmzxy),igvm2=1,stars%ng2-1)
            enddo    
            endif   
         enddo   
         CLOSE (npotmatfile)
         deallocate(vpw)
         if(.not.oneD%odi%d1)then
            deallocate(vxy)
         endif
      endif   
    
cccccccccccccccc   end of the potential part  ccccccccccc
      wannierspin=DIMENSION%jspd
      if(noco%l_soc) wannierspin=2
     
      allocate ( kveclo(nlotot),nv(wannierspin) )
      allocate ( kveclo_b(nlotot),nv_b(wannierspin) )
      IF(l_gwf) allocate ( kveclo_qb(nlotot),nv_qb(wannierspin) )
      allocate ( k1(DIMENSION%nvd,wannierspin) )
      allocate ( k2(DIMENSION%nvd,wannierspin) )
      allocate ( k3(DIMENSION%nvd,wannierspin) )
      allocate ( k1_b(DIMENSION%nvd,wannierspin) )
      allocate ( k2_b(DIMENSION%nvd,wannierspin) )
      allocate ( k3_b(DIMENSION%nvd,wannierspin) )
      IF(l_gwf) THEN
         allocate ( k1_qb(DIMENSION%nvd,wannierspin) )
         allocate ( k2_qb(DIMENSION%nvd,wannierspin) )
         allocate ( k3_qb(DIMENSION%nvd,wannierspin) )
      END IF

      allocate ( ff(atoms%ntype,atoms%jmtd,2,0:atoms%lmaxd,2) )
      allocate ( gg(atoms%ntype,atoms%jmtd,2,0:atoms%lmaxd,2) )
      allocate ( usdus%us(0:atoms%lmaxd,atoms%ntype,2) )
      allocate ( usdus%uds(0:atoms%lmaxd,atoms%ntype,2) )
      allocate ( usdus%dus(0:atoms%lmaxd,atoms%ntype,2) )
      allocate ( usdus%duds(0:atoms%lmaxd,atoms%ntype,2) )
      allocate ( usdus%ddn(0:atoms%lmaxd,atoms%ntype,2) )
      allocate ( usdus%ulos(atoms%nlod,atoms%ntype,2) )
      allocate ( usdus%dulos(atoms%nlod,atoms%ntype,2) )
      allocate ( usdus%uulon(atoms%nlod,atoms%ntype,2) )
      allocate ( usdus%dulon(atoms%nlod,atoms%ntype,2) )
      allocate ( usdus%uloulopn(atoms%nlod,atoms%nlod,atoms%ntype,2) )

      if(l_gwf.and..not.(wann%l_wann_plot)) then
         doublespin_max=4!2
      else
         doublespin_max=wannierspin
      endif

      t_it=0.0
!      t_kov=0.0
      t_qov=0.0
c*****************************************************************c
c                         START Q LOOP                            c
c   standard functionality of code for fullnqpts = nntot_q = 1    c
c        and wann%l_ms = wann%l_sgwf = wann%l_socgwf = F          c
c*****************************************************************c
      do 314 iqpt = 1,fullnqpts  ! loop by q-points starts

       ALLOCATE(innerEig_idList(nntot_q))

       call cpu_time(tt1)
!       t_kov=0.0
       t_qov=0.0

        qptibz=iqpt                          
        if(wann%l_bzsym .AND. l_gwf) qptibz=irreduc_q(iqpt)
        if(wann%l_bzsym .AND. l_gwf) oper_q=mapqoper(iqpt)

       qpt_i = noco%qss
       alph_i = noco%alph
       beta_i = noco%beta
       theta_i = noco%theta
       phi_i = noco%phi
       if(wann%l_sgwf.or.wann%l_ms) then
          qpt_i(:) = wann%param_vec(:,qptibz)
          alph_i(:) = wann%param_alpha(:,qptibz)
       elseif(wann%l_socgwf) then 
          if(wann%l_dim(2)) phi_i = tpi*wann%param_vec(2,qptibz)
          if(wann%l_dim(3)) theta_i = tpi*wann%param_vec(3,qptibz)
       endif

       IF (l_gwf) THEN
        IF(wann%l_matrixmmn)THEN
          do iqpt_b=1,nntot_q

            innerEig_idList(iqpt_b) = eig_idList(bpt_q(iqpt_b,iqpt))

!            WRITE(fending,'("_",i4.4)')bpt_q(iqpt_b,iqpt)
!            innerEig_idList(iqpt_b)=open_eig(mpi%mpi_comm,
!     +                  DIMENSION%nbasfcn,DIMENSION%neigd,
!     +                  nkpts,wannierspin,atoms%lmaxd,
!     +                  atoms%nlod,atoms%ntype,atoms%nlotot,
!     +                  noco%l_noco,.FALSE.,l_real,noco%l_soc,.FALSE.,
!     +                  mpi%n_size,filename=trim(fstart)//fending,
!     +                  layers=vacuum%layers,nstars=vacuum%nstars,
!     +                  ncored=DIMENSION%nstd,nsld=atoms%nat,
!     +                  nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,
!     +                  l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

          enddo
        ENDIF

        eig_id = eig_idList(qptibz)

!        WRITE(fending,'("_",i4.4)')qptibz
!        eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,
!     +                  nkpts,wannierspin,atoms%lmaxd,
!     +                  atoms%nlod,atoms%ntype,atoms%nlotot,
!     +                  noco%l_noco,.FALSE.,l_real,noco%l_soc,.FALSE.,
!     +                  mpi%n_size,filename=trim(fstart)//fending,
!     +                  layers=vacuum%layers,nstars=vacuum%nstars,
!     +                  ncored=DIMENSION%nstd,nsld=atoms%nat,
!     +                  nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,
!     +                  l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

       ELSEIF(wann%l_ms) THEN

        eig_id = eig_idList(qptibz)

!        WRITE(fending,'("_",i4.4)')qptibz
!        eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,
!     +                  nkpts,wannierspin,atoms%lmaxd,
!     +                  atoms%nlod,atoms%ntype,atoms%nlotot,
!     +                  noco%l_noco,.FALSE.,l_real,noco%l_soc,.FALSE.,
!     +                  mpi%n_size,filename=trim(fstart)//fending,
!     +                  layers=vacuum%layers,nstars=vacuum%nstars,
!     +                  ncored=DIMENSION%nstd,nsld=atoms%nat,
!     +                  nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,
!     +                  l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

       ELSE
          fending=''
       ENDIF ! l_gwf.or.wann%l_ms
       nrec=0
       nrec_b=0


c****************************************************
c cycle by spins starts! 
c****************************************************
      do 110 doublespin=1,doublespin_max   ! cycle by spins

         jspin=mod(doublespin+1,2)+1
         jspin_b=jspin
         if(doublespin.eq.3) jspin_b=2
         if(doublespin.eq.4) jspin_b=1

         nrec_b = nrec

         if(.not.noco%l_noco) then
            nrec = (jspin-1)*nkpts
            nrec_b = (jspin_b-1)*nkpts
         endif

       ! spin-dependent sign of the q-dependent phase
       ! in the generalized Bloch theorem
       ! -1: spin up, +1: spin down
       sign_q = -sign_q

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
       if(noco%l_soc .and. input%jspins.eq.1)jspin2=1
       jspin2_b=jspin_b
       if(noco%l_soc .and. input%jspins.eq.1)jspin2_b=1

       jsp_start = jspin ; jsp_end = jspin

       call cpu_time(tt2)
cccccccccccc   read in the eigenvalues and vectors   cccccc
       write(*,*)'wannierspin',wannierspin
       do jspin5=1,wannierspin!1!2
!       jspin5=jspin
       jsp_start=jspin5; jsp_end=jspin5
       nrec5=0
       if(.not.noco%l_noco) nrec5 = (jspin5-1)*nkpts

       call cdn_read0(eig_id,mpi%irank,mpi%isize,jspin5,DIMENSION%jspd, !wannierspin instead of DIMENSION%jspd?
     >                noco%l_noco,
     <                ello,evac,epar,wk,n_bands,n_size)

       enddo
!       call cpu_time(tt3)
       !if(l_p0) write(*,*)'cdn_read0=',tt3-tt2

c..   now we want to define the maximum number of the bands by all kpts
      nbnd = 0

      i_rec = 0 ; n_rank = 0
c*************************************************************
c..writing down the eig.1 and/or eig.2 files

c..write individual files if multi-spiral mode wann%l_ms=T
c*************************************************************
      if(l_p0)then         
         call wann_write_eig(
     >            eig_id,l_real,
     >            atoms%lmaxd,atoms%ntype,atoms%nlod,DIMENSION%neigd,
     >            DIMENSION%nvd,wannierspin,
     >            mpi%isize,jspin,DIMENSION%nbasfcn,nlotot,
     >            noco%l_ss,noco%l_noco,nrec,fullnkpts,
     >            wann%l_bzsym,wann%l_byindex,wann%l_bynumber,
     >            wann%l_byenergy,
     >            irreduc,oneD%odi,wann%band_min(jspin),
     >            wann%band_max(jspin),
     >            numbands,
     >            sliceplot%e1s,sliceplot%e2s,efermi,.false.,nkpts,
     <            nbnd,kpoints,l_gwf,iqpt)       

         if(oneD%odi%d1)then
            kpoints(:)=kpoints(:)*cell%bmat(3,3)         
         endif
      endif!l_p0
!      call cpu_time(tt2)
      !if(l_p0)write(*,*)'write_eig=',tt2-tt3

     
! nbnd is calculated for process zero and is sent here to the others
#ifdef CPP_MPI
      if(l_p0)then
         do cpu_index=1,mpi%isize-1
      call MPI_SEND(nbnd,1,MPI_INTEGER,cpu_index,1,mpi%mpi_comm,ierr)
         enddo
      else
       call MPI_RECV(nbnd,1,MPI_INTEGER,0,1,mpi%mpi_comm,stt,ierr)
      endif
#endif
     
       print*,"process: ",mpi%irank," nbnd= ",nbnd
!       call cpu_time(tt3)
c##################################################################
      if(wann%l_mmn0)then
           allocate ( mmn(nbnd,nbnd,fullnkpts) )
           mmn(:,:,:) = cmplx(0.,0.)
           if((noco%l_soc.or.noco%l_noco) .and. (doublespin.eq.1))
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

      !write (*,*) 'nwfs=',nwfs
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

      if(wann%l_matrixmmn)then
         if(.not.allocated(mmnk_q).AND.l_gwf)then
             allocate ( mmnk_q (nbnd,nbnd,nntot_q,fullnkpts) )
             mmnk_q = (0.,0.)

!             allocate ( m_int(nbnd,nbnd,nntot_q,fullnkpts) )
!             allocate ( m_sph(nbnd,nbnd,nntot_q,fullnkpts) )
!             allocate ( m_vac(nbnd,nbnd,nntot_q,fullnkpts) )
!             m_int = cmplx(0.,0.)
!             m_sph = cmplx(0.,0.)
!             m_vac = cmplx(0.,0.)
         endif
      endif

 
      allocate ( flo(atoms%ntype,atoms%jmtd,2,atoms%nlod,2) )
!      call cpu_time(tt2)
      !if(l_p0) write(*,*)'allocate',tt2-tt3

      do jspin4=1,wannierspin!2
         jspin3=jspin4
         if(input%jspins.eq.1) jspin3=1
      na = 1
      do 40 n = 1,atoms%ntype
       do 30 l = 0,atoms%lmax(n)
c...compute the l-dependent, k-independent radial MT- basis functions

         call radfun(
     >              l,n,jspin4,epar(l,n,jspin3),vr(1,n,jspin3),atoms,
     <              ff(n,:,:,l,jspin4),gg(n,:,:,l,jspin4),usdus,
     <              nodeu,noded,wronk)

   30  continue
c...and the local orbital radial functions
       do ilo = 1, atoms%nlo(n)

         call radflo(
     >             atoms,n,jspin4,ello(:,:,jspin3),vr(1,n,jspin3),
     >             ff(n,1:,1:,0:,jspin4),gg(n,1:,1:,0:,jspin4),mpi,
     <             usdus,uuilon,duilon,ulouilopn,flo(n,:,:,:,jspin4))

       enddo
c       na = na + atoms%neq(n)
   40 continue
      enddo!jspin3
!      call cpu_time(tt3)
      !if(l_p0) write(*,*)'radfun=',tt3-tt2

c****************************************************************
c   calculate the k-independent uju*gaunt-matrix needed for
c   mmnmatrix
c****************************************************************
! TODO: make this more efficient (i.e., compute ujugaunt only once
! and not for all q-points).
      if(wann%l_matrixmmn)then
         call cpu_time(delta)
         allocate(ujug(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot))
         allocate(ujdg(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot))
         allocate(djug(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot))
         allocate(djdg(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot))
         allocate(ujulog(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot))
         allocate(djulog(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot))
         allocate(ulojug(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot))
         allocate(ulojdg(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot))
         allocate(ulojulog(1:atoms%nlod,-atoms%llod:atoms%llod,
     &                     1:atoms%nlod,-atoms%llod:atoms%llod,
     &                     1:atoms%ntype,1:nntot))

         call wann_ujugaunt(
     >            atoms%llod,nntot,kdiff,atoms%lmax,atoms%ntype,
     >            atoms%ntype,cell%bbmat,cell%bmat,atoms%nlod,atoms%nlo,
     >            atoms%llo,flo(:,:,:,:,jspin),
     >            flo(:,:,:,:,jspin),
     >            ff(:,:,:,:,jspin),
     >            ff(:,:,:,:,jspin),
     >            gg(:,:,:,:,jspin),
     >            gg(:,:,:,:,jspin),atoms%jri,atoms%rmsh,atoms%dx,
     >            atoms%jmtd,atoms%lmaxd,lmd,
     <            ujug,ujdg,djug,djdg,
     <            ujulog,djulog,ulojug,ulojdg,ulojulog,.false.,1)

! compute integrals of radial solution, according energy derivatives,
! the spherical Bessel function and the Gaunt coefficients in order
! to account for the overlap of the lattice periodic parts at
! neighboring q-points
         IF(l_gwf)THEN
            allocate(ujug_q(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot_q))
            allocate(ujdg_q(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot_q))
            allocate(djug_q(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot_q))
            allocate(djdg_q(0:lmd,0:lmd,
     &              1:atoms%ntype,1:nntot_q))
            allocate(ujulog_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot_q))
            allocate(djulog_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot_q))
            allocate(ulojug_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot_q))
            allocate(ulojdg_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,
     &              1:atoms%ntype,1:nntot_q))
            allocate(ulojulog_q(1:atoms%nlod,-atoms%llod:atoms%llod,
     &                          1:atoms%nlod,-atoms%llod:atoms%llod,
     &                          1:atoms%ntype,1:nntot_q))
            
            ! we need G(q+b)/2 as argument for the sph. Bessel func.
            ! and additionally a spin-dependent sign (-/+ 1)^{lpp}
            if(wann%l_sgwf) call wann_ujugaunt(
     >            atoms%llod,nntot_q,qdiff/2.0,atoms%lmax,atoms%ntype,
     >            atoms%ntype,cell%bbmat,cell%bmat,atoms%nlod,atoms%nlo,
     >            atoms%llo,flo(:,:,:,:,jspin),
     >            flo(:,:,:,:,jspin_b),
     >            ff(:,:,:,:,jspin),
     >            ff(:,:,:,:,jspin_b),
     >            gg(:,:,:,:,jspin),
     >            gg(:,:,:,:,jspin_b),atoms%jri,atoms%rmsh,atoms%dx,
     >            atoms%jmtd,atoms%lmaxd,lmd,
     <            ujug_q,ujdg_q,djug_q,djdg_q,
     <            ujulog_q,djulog_q,ulojug_q,ulojdg_q,ulojulog_q,.true.,
     >            sign_q)

            if(wann%l_socgwf) call wann_ujugaunt(
     >            atoms%llod,nntot_q,zero_qdiff,atoms%lmax,atoms%ntype,
     >            atoms%ntype,cell%bbmat,cell%bmat,atoms%nlod,atoms%nlo,
     >            atoms%llo,flo(:,:,:,:,jspin),
     >            flo(:,:,:,:,jspin_b),
     >            ff(:,:,:,:,jspin),
     >            ff(:,:,:,:,jspin_b),
     >            gg(:,:,:,:,jspin),
     >            gg(:,:,:,:,jspin_b),atoms%jri,atoms%rmsh,atoms%dx,
     >            atoms%jmtd,
     >            atoms%lmaxd,lmd,ujug_q,ujdg_q,djug_q,djdg_q,
     <            ujulog_q,djulog_q,ulojug_q,ulojdg_q,ulojulog_q,
     >            .false.,1)
      
         ENDIF ! l_gwf

         call cpu_time(delta1)
         time_ujugaunt=delta1-delta
      endif !l_matrixmmn
!      call cpu_time(tt2)
      !if(l_p0) write(*,*)'ujugaunt=',tt2-tt3

      zzMat%l_real = l_real
      zzMat%nbasfcn = DIMENSION%nbasfcn
      zzMat%nbands = DIMENSION%neigd
      IF(l_real) THEN
	    IF(.not.ALLOCATED(zzMat%z_r))
     >    ALLOCATE (zzMat%z_r(zzMat%nbasfcn,zzMat%nbands))
      ELSE
	    IF(.not.ALLOCATED(zzMat%z_c))
     >    ALLOCATE (zzMat%z_c(zzMat%nbasfcn,zzMat%nbands))
      END IF

      zMat%l_real = zzMat%l_real
      zMat%nbasfcn = zzMat%nbasfcn
      zMat%nbands = zzMat%nbands
      IF (zzMat%l_real) THEN
	    IF(.not.ALLOCATED(zMat%z_r))
     >    ALLOCATE (zMat%z_r(zMat%nbasfcn,zMat%nbands))
         zMat%z_r = 0.0
      ELSE
	    IF(.not.ALLOCATED(zMat%z_c))
     >    ALLOCATE (zMat%z_c(zMat%nbasfcn,zMat%nbands))
         zMat%z_c = CMPLX(0.0,0.0)
      END IF

      zMat_b%l_real = zzMat%l_real
      zMat_b%nbasfcn = zzMat%nbasfcn
      zMat_b%nbands = zzMat%nbands
      IF (zzMat%l_real) THEN
	    IF(.not.ALLOCATED(zMat_b%z_r))
     >    ALLOCATE (zMat_b%z_r(zMat_b%nbasfcn,zMat_b%nbands))
         zMat_b%z_r = 0.0
      ELSE
	    IF(.not.ALLOCATED(zMat_b%z_c))
     >    ALLOCATE (zMat_b%z_c(zMat_b%nbasfcn,zMat_b%nbands))
         zMat_b%z_c = CMPLX(0.0,0.0)
      END IF

      i_rec = 0 ; n_rank = 0
      
c****************************************************************
c.. loop by kpoints starts!      each may be a separate task
c****************************************************************
      do 10 ikpt = wann%ikptstart,fullnkpts  ! loop by k-points starts
!         call cpu_time(tt2)
!         call cpu_time(tt3)
        kptibz=ikpt
        if(wann%l_bzsym) kptibz=irreduc(ikpt)
        if(wann%l_bzsym) oper=mapkoper(ikpt)

        i_rec = i_rec + 1
      if (mod(i_rec-1,mpi%isize).eq.mpi%irank) then

      allocate ( we(DIMENSION%neigd),eigg(DIMENSION%neigd) )

      call cpu_time(delta)
      n_start=1
      n_end=DIMENSION%neigd


! read information of diagonalization for fixed q-point iqpt
! stored in the eig file on unit 66. the lattice respectively
! plane-wave vectors G(k,q) are saved in (k1,k2,k3).

      CALL cdn_read(
     >              eig_id,
     >              DIMENSION%nvd,DIMENSION%jspd,mpi%irank,mpi%isize, !wannierspin instead of DIMENSION%jspd?
     >              kptibz,jspin,DIMENSION%nbasfcn,
     >              noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,n_end,
     <              ello,evdu,epar,
     <              wk,nbands,eigg,zzMat)

!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'cdn_read=',tt4-tt3

      call cpu_time(delta1)
      time_rw=time_rw+delta1-delta
      nslibd = 0

c...we work only within the energy window

      eig(:) = 0.

!      print*,"bands used:"

      do i = 1,nbands
        if ((eigg(i).ge.sliceplot%e1s.and.nslibd.lt.numbands.and.
     &       wann%l_bynumber).or.
     &      (eigg(i).ge.sliceplot%e1s.and.eigg(i).le.sliceplot%e2s.and.
     &       wann%l_byenergy).or.(i.ge.wann%band_min(jspin).and.
     &      (i.le.wann%band_max(jspin)).and.wann%l_byindex))then

!           print*,i
           nslibd = nslibd + 1
           eig(nslibd) = eigg(i)
           we(nslibd) = we(i)
           if(noco%l_noco)then
             funbas=       nv(1)+nlotot
             funbas=funbas+nv(2)+nlotot
           else
             funbas=nv(jspin)+nlotot
           endif
           IF(zzMat%l_real) THEN
              do j = 1, funbas
                 zMat%z_r(j,nslibd) = zzMat%z_r(j,i)
              end do
           ELSE
              do j = 1, funbas
                 zMat%z_c(j,nslibd) = zzMat%z_c(j,i)
              end do
           END IF
        endif
      enddo

c***********************************************************
c              rotate the wavefunction
c***********************************************************
      if (wann%l_bzsym.and.oper.ne.1) then  !rotate bkpt
!         call wann_kptsrotate(
!     >            atoms%nat,atoms%nlod,atoms%llod,
!     >            atoms%ntype,atoms%nlo,atoms%llo,atoms%invsat,
!     >            noco%l_noco,noco%l_soc,
!     >            atoms%ntype,atoms%neq,nlotot,
!     >            kveclo,jspin,
!     >            oper,sym%nop,sym%mrot,DIMENSION%nvd,nv,
!     >            shiftkpt(:,ikpt),
!     >            sym%tau,
!     x            bkpt,k1(:,:),k2(:,:),k3(:,:),
!     x            zMat,nsfactor)
      else
         nsfactor=cmplx(1.0,0.0)
      endif
c      print*,"bkpt1=",bkpt
!      call cpu_time(tt3)
      !if(l_p0) write(*,*)'nbnd=',tt3-tt4

c******************************************************************

c...the overlap matrix Mmn which is computed for each k- and b-point

      noccbd = nslibd

      allocate(acof(noccbd,0:lmd,atoms%nat),
     &         bcof(noccbd,0:lmd,atoms%nat),
     &         ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))

      acof(:,:,:) = cmplx(0.,0.) ; bcof(:,:,:) = cmplx(0.,0.)
      ccof(:,:,:,:) = cmplx(0.,0.)

c...generation the A,B,C coefficients in the spheres
c...for the lapws and local orbitals, summed by the basis functions

      call cpu_time(delta)

      ALLOCATE(lapw%k1(SIZE(k1,1),SIZE(k1,2)))
      ALLOCATE(lapw%k2(SIZE(k1,1),SIZE(k1,2)))
      ALLOCATE(lapw%k3(SIZE(k1,1),SIZE(k1,2)))
      lapw%k1 = k1
      lapw%k2 = k2
      lapw%k3 = k3
      lapw%nmat = nmat
      lapw%nv = nv
      ! I think the other variables of lapw are not needed here.

!      CALL abcof(input,atoms,noccbd,sym,cell,bkpt,lapw,noccbd,usdus,
!     >           noco,jspin,kveclo,oneD,acof,bcof,ccof,zMat)

      DEALLOCATE(lapw%k1,lapw%k2,lapw%k3)

      call cpu_time(delta1)
      time_abcof=time_abcof+delta1-delta

      call wann_abinv(
     >        atoms%ntype,atoms%nat,noccbd,atoms%lmaxd,lmd,atoms%llod,
     >        atoms%nlod,atoms%ntype,atoms%neq,noccbd,atoms%lmax,
     >        atoms%nlo,atoms%llo,atoms%invsat,sym%invsatnr,bkpt,
     >        atoms%taual,
     X        acof,bcof,ccof)
!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'abcof=',tt4-tt3


      if((doublespin.eq.3).or.(doublespin.eq.4)) goto 9900


      if(wann%l_orbcomp)then
         call wann_orbcomp(
     >            atoms%llod,noccbd,atoms%nlod,atoms%nat,atoms%ntype,
     >            atoms%lmaxd,lmd,atoms%ntype,atoms%neq,atoms%nlo,
     >            atoms%llo,acof,bcof,ccof,
     >            usdus%ddn(:,:,jspin),usdus%uulon(:,:,jspin),
     >            usdus%dulon(:,:,jspin),usdus%uloulopn(:,:,:,jspin),
     >            wann%oc_num_orbs,
     >            wann%oc_orbs,
     >            wann%l_oc_f,
     =            orbcomp(:,:,:,:,ikpt))
      endif

      if(wann%l_anglmom)then
         call wann_anglmom(
     >                  atoms%llod,noccbd,atoms%nlod,atoms%nat,
     >                  atoms%ntype,atoms%lmax,lmd,!atoms%lmaxd,lmd,
     >                  atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,
     >                  acof,bcof,ccof,
     >                  usdus%ddn(:,:,jspin),usdus%uulon(:,:,jspin),
     >                  usdus%dulon(:,:,jspin),
     >                  usdus%uloulopn(:,:,:,jspin),
     =                  anglmom(:,:,:,ikpt))
      endif

#ifdef CPP_TOPO
      if(wann%l_surfcurr)then
c         call wann_surfcurr_int(
c     >        DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,vacuum%nmzxyd,stars%ng2,sphhar%ntypsd,
c     >        atoms%ntype,atoms%lmaxd,atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,atoms%neq,stars%ng3,vacuum%nvac,
c     >        vacuum%nmz,vacuum%nmzxy,stars%ng2,sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,sym%symor,
c     >        sym%invs,sym%invs2,cell%z1,vacuum%delz,atoms%ngopr,atoms%ntypsy,atoms%jri,atoms%pos,atoms%zatom,
c     >        atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,cell%amat,cell%bmat,cell%bbmat,ikpt,sliceplot%nnne,sliceplot%kk,
c     >        DIMENSION%nvd,atoms%nlod,atoms%llod,nv(jspin),lmd,bkpt,cell%omtil,atoms%nlo,atoms%llo,
c     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
c     >        vz(:,:,jspin2),
c     >        nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,ff,gg,flo,acof,bcof,ccof,z,
c     >        surfcurr(:,:,:,ikpt))

         call wann_surfcurr_int2(
     >            DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,
     >            vacuum%nmzxyd,
     >            stars%ng2,sphhar%ntypsd,atoms%ntype,atoms%lmaxd,
     >            atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,
     >            atoms%neq,stars%ng3,vacuum%nvac,vacuum%nmz,
     >            vacuum%nmzxy,stars%ng2,sym%nop,sym%nop2,cell%volint,
     >            input%film,sliceplot%slice,sym%symor,
     >            sym%invs,sym%invs2,cell%z1,vacuum%delz,atoms%ngopr,
     >            atoms%ntypsy,atoms%jri,atoms%pos,atoms%taual,
     >            atoms%zatom,atoms%rmt,atoms%lmax,sym%mrot,sym%tau,
     >            atoms%rmsh,sym%invtab,cell%amat,cell%bmat,cell%bbmat,
     >            ikpt,DIMENSION%nvd,nv(jspin),bkpt,cell%omtil,
     >            k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >            nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,z,
     <            dirfacs,
     >            surfcurr(:,:,:,ikpt))

         call wann_surfcurr(
     >            dirfacs,cell%amat,
     >            jspin,atoms%ntype,atoms%lmaxd,atoms%lmax,atoms%nat,
     >            atoms%neq,noccbd,lmd,atoms%nat,atoms%llod,atoms%nlod,
     >            atoms%nlo,atoms%llo, 
     >            acof,bcof,ccof,
     >            us(:,:,jspin),dus(:,:,jspin),duds(:,:,jspin),
     >            uds(:,:,jspin),
     >            ulos(:,:,jspin),dulos(:,:,jspin),
     >            atoms%rmt,atoms%pos, 
     &            surfcurr(:,:,:,ikpt))
        write(6,*)"dirfacs=",dirfacs
      endif

      if(wann%l_soctomom)then
         call wann_soc_to_mom(
     >            jspin,atoms%ntype,atoms%lmaxd,atoms%lmax,atoms%nat,
     >            atoms%jmtd,atoms%jri,atoms%rmsh,atoms%dx,atoms%neq,
     >            noccbd,lmd,atoms%nat,atoms%llod,atoms%nlod,
     >            vso(:,:,1), 
     >            ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),
     >            acof,bcof,ccof,
     &            soctomom(:,:,:,ikpt))
      endif

      if(wann%l_nabla)then
         call wann_nabla(
     >       atoms%nlo,atoms%llo,
     >       jspin,atoms%ntype,atoms%lmaxd,atoms%lmax,atoms%nat,
     >       atoms%jmtd,atoms%jri,atoms%rmsh,atoms%dx,atoms%neq,
     >       noccbd,lmd,atoms%nat,atoms%llod,atoms%nlod, 
     >       ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),flo(:,:,:,:,jspin),
     >       acof,bcof,ccof,
     &       nablamat(:,:,:,ikpt))
         if(input%film.and..not.oneD%odi%d1)then
            call wann_nabla_vac(
     >          cell%z1,vacuum%nmzd,DIMENSION%nv2d,
     >          stars%mx1,stars%mx2,stars%mx3,
     >          stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,vacuum%delz,
     >          stars%ig2,cell%area,cell%bmat,cell%bbmat,evac(:,jspin),
     >          bkpt,vz(:,:,jspin2),nslibd,jspin,k1,k2,k3,wannierspin,
     >          DIMENSION%nvd,DIMENSION%nbasfcn,DIMENSION%neigd,z,nv,
     >          cell%omtil,
     <          nablamat(:,:,:,ikpt))
         endif   
         addnoco=0
         do 41 i = n_rank+1,nv(jspin),n_size
           b1(1)=bkpt(1)+k1(i,jspin)
           b1(2)=bkpt(2)+k2(i,jspin)
           b1(3)=bkpt(3)+k3(i,jspin)
           b2(1)=b1(1)*cell%bmat(1,1)+b1(2)*cell%bmat(2,1)+
     +           b1(3)*cell%bmat(3,1)
           b2(2)=b1(1)*cell%bmat(1,2)+b1(2)*cell%bmat(2,2)+
     +           b1(3)*cell%bmat(3,2)
           b2(3)=b1(1)*cell%bmat(1,3)+b1(2)*cell%bmat(2,3)+
     +           b1(3)*cell%bmat(3,3)
           do 42 j = n_rank+1,nv(jspin),n_size
c-->     determine index and phase factor
            i1 = k1(j,jspin) - k1(i,jspin)
            i2 = k2(j,jspin) - k2(i,jspin)
            i3 = k3(j,jspin) - k3(i,jspin)
            in = stars%ig(i1,i2,i3)
            if (in.eq.0) goto 42
            phase   = stars%rgphs(i1,i2,i3)
            phasust = cmplx(phase,0.0)*stars%ustep(in)

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
       if(noco%l_noco.and.(jspin.eq.2))then
          addnoco=nv(1)+nlotot
       endif
       
c$$$       do 21 i = n_rank+1,nv(jspin),n_size
c$$$         do 22 j = n_rank+1,nv(jspin),n_size
c$$$
c$$$c-->     determine index and phase factor
c$$$            i1 = k1(j,jspin) - k1(i,jspin)
c$$$            i2 = k2(j,jspin) - k2(i,jspin)
c$$$            i3 = k3(j,jspin) - k3(i,jspin)
c$$$c            if(abs(i1).gt.stars%mx1.or.abs(i2).gt.stars%mx2
c$$$c     & .or. abs(i3).gt.stars%mx3)print*,"interstitial warning"
c$$$            in = stars%ig(i1,i2,i3)
c$$$            if (in.eq.0) goto 22
c$$$            phase   = stars%rgphs(i1,i2,i3)
c$$$            phasust = cmplx(phase,0.0)*stars%ustep(in)
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
     >         cmplx_1,addnoco,addnoco,
     >         DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,
     >         stars%ng3,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,
     >         k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),zMat,nslibd,
     >         nbnd,
     >         stars%rgphs,stars%ustep,stars%ig,(/ 0,0,0 /),
     <         mmn(:,:,ikpt))

c---> spherical contribution to mmn0-matrix

       call wann_mmk0_sph(
     >           atoms%llod,noccbd,atoms%nlod,atoms%nat,atoms%ntype,
     >           atoms%lmaxd,atoms%lmax,lmd,atoms%ntype,atoms%neq,
     >           atoms%nlo,atoms%llo,acof(1:noccbd,:,:),
     >           bcof(1:noccbd,:,:),ccof(:,1:noccbd,:,:),
     >           usdus%ddn(:,:,jspin),usdus%uulon(:,:,jspin),
     >           usdus%dulon(:,:,jspin),usdus%uloulopn,
     =           mmn(:,:,ikpt))
c---> vacuum contribution to mmn0-matrix

       if (input%film .and. .not.oneD%odi%d1) then

            call wann_mmk0_vac(
     >           noco%l_noco,nlotot,qpt_i,
     >           cell%z1,vacuum%nmzd,DIMENSION%nv2d,
     >           stars%mx1,stars%mx2,stars%mx3,
     >           stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,vacuum%delz,
     >           stars%ig2,cell%area,cell%bmat,
     >           cell%bbmat,evac(:,jspin),bkpt,vz(:,:,jspin2),
     >           nslibd,jspin,k1,k2,k3,wannierspin,DIMENSION%nvd,
     >           DIMENSION%nbasfcn,DIMENSION%neigd,zMat,nv,cell%omtil,
     <           mmn(:,:,ikpt))
       elseif (oneD%odi%d1) then

             call wann_mmk0_od_vac(
     >           DIMENSION, oneD, vacuum, stars, cell,
     >           noco%l_noco,nlotot,
     >           cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,
     >           stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,
     >           stars%ig,vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,
     >           oneD%odi%n2d,cell%bbmat,evac(1,jspin),bkpt,oneD%odi%M,
     >           oneD%odi%mb,vz(:,1,jspin2),oneD%odi,
     >           nslibd,jspin,k1,k2,k3,wannierspin,DIMENSION%nvd,
     >           cell%area,DIMENSION%nbasfcn,DIMENSION%neigd,zMat,nv,
     >           stars%sk2,stars%phi2,cell%omtil,qpt_i,
     <           mmn(:,:,ikpt))

       endif
       endif !l_mmn0

!       call cpu_time(tt3)

c---> overlaps with the trial orbitals
      if (wann%l_projmethod.or.wann%l_bestproj.or.wann%l_matrixamn) then
         l_amn2=.false.
         amnchi = cmplx_1!cmplx(1.,0.)

         call wann_amn (
     >            amnchi,nslibd,nwfs,atoms%ntype,atoms%nlod,atoms%llod,
     >            atoms%llo,atoms%nlo,atoms%lmaxd,atoms%jmtd,lmd,
     >            atoms%neq,atoms%nat,ikpt,nbnd,
     >            atoms%rmsh,atoms%rmt,atoms%jri,atoms%dx,atoms%lmax,
     >            usdus%us(:,:,jspin),usdus%dus(:,:,jspin),
     >            usdus%uds(:,:,jspin),
     >            usdus%duds(:,:,jspin),flo(:,:,:,:,jspin),
     >            ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),acof,bcof,ccof,
     >            (noco%l_soc.or.noco%l_noco),jspin,
     &            l_amn2,amn(:,:,ikpt))
         if(l_amn2)then
          call wann_amn (
     >            amnchi,nslibd,nwfs,atoms%ntype,atoms%nlod,atoms%llod,
     >            atoms%llo,atoms%nlo,atoms%lmaxd,atoms%jmtd,lmd,
     >            atoms%neq,atoms%nat,ikpt,nbnd,
     >            atoms%rmsh,atoms%rmt,atoms%jri,atoms%dx,atoms%lmax,
     >            usdus%us(:,:,jspin),usdus%dus(:,:,jspin),
     >            usdus%uds(:,:,jspin),
     >            usdus%duds(:,:,jspin),flo(:,:,:,:,jspin),
     >            ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),acof,bcof,ccof,
     >            (noco%l_soc.or.noco%l_noco),jspin,
     &            l_amn2,amn(:,:,ikpt),bkpt)
         endif
c         amn(ikpt,:,:)=amn(ikpt,:,:)*conjg(nsfactor)
      endif

!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'amn',tt4-tt3


c****************************************************************
c...         vanderbilt mmn matrix
c***************************************************************
      if (wann%l_matrixmmn .and.
     >     (.not.wann%l_skipkov)) then   !  vanderbilt procedure Mmn matrix
      allocate ( we_b(DIMENSION%neigd) )

  !!! the cycle by the nearest neighbors (nntot) for each kpoint

       do 15  ikpt_b = 1,nntot
!      call cpu_time(tt3)

          if(pair_to_do(ikpt,ikpt_b).eq.0)cycle !save time by symmetry
          kptibz_b=bpt(ikpt_b,ikpt)
          if(wann%l_bzsym) oper_b=mapkoper(kptibz_b)
          if (wann%l_bzsym) kptibz_b=irreduc(kptibz_b)



! now we need the wavefunctions for k_b kpoint

c        print*,"something to do"


          call cpu_time(delta)
          n_start=1
          n_end=DIMENSION%neigd

      CALL cdn_read(
     >              eig_id,
     >              DIMENSION%nvd,DIMENSION%jspd,mpi%irank,mpi%isize, !wannierspin instead of DIMENSION%jspd?
     >              kptibz_b,jspin,DIMENSION%nbasfcn,
     >              noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,n_end,
     <              ello,evdu,epar,
     <              wk_b,nbands_b,eigg,zzMat)

!        call cpu_time(tt4)
        !if(l_p0) write(*,*)'read_eig=',tt4-tt3

        call cpu_time(delta1)
        time_rw=time_rw+delta1-delta
        nslibd_b = 0

      eig_b(:) = 0.

        do i = 1,nbands_b
          if((eigg(i).ge.sliceplot%e1s.and.nslibd_b.lt.numbands
     &        .and.wann%l_bynumber).or.
     &       (eigg(i).ge.sliceplot%e1s.and.eigg(i).le.sliceplot%e2s.and.
     &        wann%l_byenergy).or.(i.ge.wann%band_min(jspin).and.
     &        (i.le.wann%band_max(jspin)).and.
     &       wann%l_byindex))then
            nslibd_b = nslibd_b + 1
            eig_b(nslibd_b) = eigg(i)
            we_b(nslibd_b) = we_b(i)
            if(noco%l_noco)then
               funbas =        nv_b(1) + nlotot
               funbas = funbas+nv_b(2) + nlotot
            else
               funbas = nv_b(jspin) + nlotot
            endif
            IF (zzMat%l_real) THEN
               do j = 1,funbas
                  zMat_b%z_r(j,nslibd_b) = zzMat%z_r(j,i)
               enddo
            ELSE
               do j = 1,funbas
                  zMat_b%z_c(j,nslibd_b) = zzMat%z_c(j,i)
               enddo
            END IF
          endif
        enddo

c***********************************************************
c              Rotate the wavefunction of next neighbor.
c***********************************************************
      if (wann%l_bzsym .and. (oper_b.ne.1)  ) then
!         call wann_kptsrotate(
!     >            atoms%nat,atoms%nlod,atoms%llod,
!     >            atoms%ntype,atoms%nlo,atoms%llo,atoms%invsat,
!     >            noco%l_noco,noco%l_soc,
!     >            atoms%ntype,atoms%neq,nlotot,
!     >            kveclo_b,jspin,
!     >            oper_b,sym%nop,sym%mrot,DIMENSION%nvd,
!     >            nv_b,
!     >            shiftkpt(:,bpt(ikpt_b,ikpt)),
!     >            sym%tau,
!     x            bkpt_b,k1_b(:,:),
!     x            k2_b(:,:),k3_b(:,:),
!     x            zMat_b,nsfactor_b)
      else
         nsfactor_b=cmplx(1.0,0.0)
      endif
c      print*,"kpt2=",bkpt_b
!      call cpu_time(tt5)
      !if(l_p0) write(*,*)'nbnd_b=',tt5-tt4

        noccbd_b = nslibd_b

ccccc   we start with the Mmn matrix   ccccccccccccc

!!! matrix elements of the interstitial overlap
!!! matrix with the weights given by products of
!!! the c-coeff. for different G-vectors, bands and k-points
!!! Mmn(k,b)(IR) = \sum(G,G')C_G^(k,n)*C_G'^(k+b,m)\theta_(G-G')

cccccccccccc  Spherical Contributions         ccccccccccccc

      allocate (acof_b(noccbd_b,0:lmd,atoms%nat),
     &          bcof_b(noccbd_b,0:lmd,atoms%nat),
     &          ccof_b(-atoms%llod:atoms%llod,noccbd_b,atoms%nlod,
     &          atoms%nat))
      
!!! get the band-dependent k-dependent ab coeff.
      call cpu_time(delta)

      ALLOCATE(lapw_b%k1(SIZE(k1_b,1),SIZE(k1_b,2)))
      ALLOCATE(lapw_b%k2(SIZE(k1_b,1),SIZE(k1_b,2)))
      ALLOCATE(lapw_b%k3(SIZE(k1_b,1),SIZE(k1_b,2)))
      lapw_b%k1 = k1_b
      lapw_b%k2 = k2_b
      lapw_b%k3 = k3_b
      lapw_b%nmat = nmat_b
      lapw_b%nv = nv_b
      ! I think the other variables of lapw are not needed here.

!      CALL abcof(input,atoms,noccbd_b,sym,cell,bkpt_b,lapw_b,
!     +           noccbd_b,usdus,noco,jspin,kveclo_b,oneD,
!     +           acof_b,bcof_b,ccof_b,zMat_b)

      DEALLOCATE(lapw_b%k1,lapw_b%k2,lapw_b%k3)

      call cpu_time(delta1)
      time_abcof=time_abcof+delta1-delta

      call wann_abinv(
     >         atoms%ntype,atoms%nat,noccbd_b,atoms%lmaxd,lmd,
     >         atoms%llod,atoms%nlod,atoms%ntype,atoms%neq,
     >         noccbd_b,atoms%lmax,atoms%nlo,atoms%llo,atoms%invsat,
     >         sym%invsatnr,bkpt_b,atoms%taual,
     X         acof_b,bcof_b,ccof_b)
!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'abcof_b=',tt4-tt5

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
      if(noco%l_noco.and.(jspin.eq.2))then
         addnoco  = nv(1)   + nlotot
         addnoco2 = nv_b(1) + nlotot
      endif

      call wann_mmkb_int(
     >         cmplx_1,addnoco,addnoco2,
     >         DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,
     >         stars%ng3,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,
     >         k1_b(:,jspin),k2_b(:,jspin),k3_b(:,jspin),
     >         nv_b(jspin),zMat_b,nslibd_b,
     >         nbnd,
     >         stars%rgphs,stars%ustep,stars%ig,gb(:,ikpt_b,ikpt),
     <         mmnk(:,:,ikpt_b,ikpt))

!      call cpu_time(tt5)
      !if(l_p0) write(*,*)'int=',tt5-tt4

        call cpu_time(delta1)
        time_interstitial=time_interstitial+delta1-delta

cccccccccccc  Spherical Contributions  ccccccccccccc
        call cpu_time(delta)

        chi = cmplx_1
         call wann_mmkb_sph(
     >       nbnd,atoms%llod,nslibd,nslibd_b,atoms%nlod,atoms%nat,
     >       atoms%ntype,lmd,atoms%jmtd,atoms%taual,sym%nop,atoms%lmax,
     >       atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,acof,bcof,ccof,
     >       bkpt_b,acof_b,bcof_b,ccof_b,gb(:,ikpt_b,ikpt),bkpt,
     >       ujug,ujdg,
     >       djug,djdg,ujulog,djulog,ulojug,ulojdg,ulojulog,kdiff,
     >       nntot,chi,
     =       mmnk(:,:,ikpt_b,ikpt))


!         call cpu_time(tt4)
         !if(l_p0) write(*,*)'sph=',tt4-tt5

         call cpu_time(delta1)
         time_mmn=time_mmn+delta1-delta

c...vacuum contributions
         call cpu_time(delta)
        if (input%film .and. .not.oneD%odi%d1) then

           call wann_mmkb_vac(
     >              cmplx_1,noco%l_noco,nlotot,qpt_i,
     >              nbnd,cell%z1,vacuum%nmzd,DIMENSION%nv2d,
     >              stars%mx1,stars%mx2,stars%mx3,
     >              stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,
     >              vacuum%delz,stars%ig2,cell%area,cell%bmat,
     >              cell%bbmat,evac(:,jspin),evac(:,jspin_b),
     >              bkpt,bkpt_b,vz(:,:,jspin2),vz(:,:,jspin2_b),
     >              nslibd,nslibd_b,jspin,jspin_b,
     >              k1,k2,k3,k1_b,k2_b,k3_b,
     >              wannierspin,DIMENSION%nvd,
     >              DIMENSION%nbasfcn,DIMENSION%neigd,zMat,zMat_b,
     >              nv,nv_b,cell%omtil,
     >              gb(:,ikpt_b,ikpt),
     <              mmnk(:,:,ikpt_b,ikpt))
         elseif (oneD%odi%d1) then

           call wann_mmkb_od_vac(
     >           DIMENSION,oneD,vacuum,stars,cell,
     >           cmplx_1,noco%l_noco,nlotot,
     >           nbnd,cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,
     >           stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,
     >           stars%ig,vacuum%nmzxy,
     >           vacuum%nmz,vacuum%delz,stars%ig2,oneD%odi%n2d,
     >           cell%bbmat,evac(1,jspin),evac(1,jspin_b),
     >           bkpt,bkpt_b,oneD%odi%M,oneD%odi%mb,
     >           vz(:,1,jspin2),vz(:,1,jspin2_b),oneD%odi,
     >           nslibd,nslibd_b,jspin,jspin_b,k1,k2,k3,k1_b,k2_b,k3_b,
     >           wannierspin,DIMENSION%nvd,cell%area,DIMENSION%nbasfcn,
     >           DIMENSION%neigd,
     >           zMat,zMat_b,nv,nv_b,stars%sk2,stars%phi2,cell%omtil,
     >           gb(:,ikpt_b,ikpt),qpt_i,
     >           .false.,1,
     <           mmnk(:,:,ikpt_b,ikpt))
        endif

c        mmnk(:,:,ikpt_b,ikpt)=
c     &          mmnk(:,:,ikpt_b,ikpt)*nsfactor*conjg(nsfactor_b)
!        call cpu_time(tt5)
        !if(l_p0) write(*,*)'film=',tt5-tt4

        call cpu_time(delta1)
        time_film=time_film+delta1-delta

        deallocate ( acof_b,bcof_b,ccof_b )

!        call cpu_time(tt4)
        !if(l_p0) write(*,*)'k-neig',ikpt_b,' total=',tt4-tt3
!        t_kov=t_kov+tt4-tt3

15      continue ! end of loop by the nearest k-neighbors

        deallocate ( we_b )
      
      endif!l_matrixmmn=.true.

 9900 continue ! jump for doublespin loop

      if (wann%l_matrixmmn) then   !  vanderbilt procedure Mmn matrix

c*******************************************c 
c          START Q-NEIGHBOR LOOP            c
c*******************************************c      
      allocate ( we_qb(DIMENSION%neigd) )

      DO iqpt_b=1,nntot_q 
!         call cpu_time(tt3)
        IF(.NOT.l_gwf) EXIT              ! old functionality

        qptibz_b = bpt_q(iqpt_b,iqpt)
        IF(qptibz_b.EQ.qptibz) CYCLE     ! no need to compute overlaps
                                         ! with periodic images for now

        qptb_i = noco%qss
        alphb_i = noco%alph
        betab_i = noco%beta
        thetab_i = noco%theta
        phib_i = noco%phi
        if(wann%l_sgwf) then
           qptb_i(:) = wann%param_vec(:,qptibz_b)
           alphb_i(:) = wann%param_alpha(:,qptibz_b)
        elseif(wann%l_socgwf) then
           if(wann%l_dim(2)) phib_i = tpi*wann%param_vec(2,qptibz_b)
           if(wann%l_dim(3)) thetab_i = tpi*wann%param_vec(3,qptibz_b)
        endif

        !if(pair_to_do_q(iqpt,iqpt_b).eq.0)cycle    ! TODO: use symmetry
        if(wann%l_bzsym) oper_qb=mapqoper(qptibz_b)
        if (wann%l_bzsym) qptibz_b=irreduc_q(qptibz_b)

        n_start=1
        n_end=DIMENSION%neigd
        
! read in diagonalization information from corresponding
! eig file to q-point iqpt_b at a given k-point ikpt.
! as a check verify that bkpt.eq.bqpt (same k).
! moreover, the plane-wave vectors G(k,q+b) are stored
! in (k1_qb,k2_qb,k3_qb) for later use.

       CALL cdn_read(
     >               innerEig_idList(iqpt_b),
     >               DIMENSION%nvd,DIMENSION%jspd,mpi%irank,mpi%isize, !wannierspin instead of DIMENSION%jspd? !kptibz_b2?
     >               kptibz,jspin_b,DIMENSION%nbasfcn,
     >               noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,
     >               n_end,
     <               ello,evdu,epar,
     <               wqb,nbands_qb,eigg,   !we_qb instead of wqb?
     <               zzMat)

!        call cpu_time(tt4)
        !if(l_p0) write(*,*)'read_eig=',tt4-tt3

       ! are we dealing with the same k-point at which
       ! we want to construct A,B,C coefficients etc. ?
       IF(ANY(bqpt.NE.bkpt)) CALL juDFT_error("bqpt.ne.bkpt",
     >                             calledby="wannier")

       zMat_qb%l_real = zzMat%l_real
       zMat_qb%nbasfcn = zzMat%nbasfcn
       zMat_qb%nbands = zzMat%nbands
       IF (zzMat%l_real) THEN
          ALLOCATE (zMat_qb%z_r(zMat%nbasfcn,zMat%nbands))
          zMat_qb%z_r = 0.0
       ELSE
          ALLOCATE (zMat_qb%z_c(zMat%nbasfcn,zMat%nbands))
          zMat_qb%z_c = CMPLX(0.0,0.0)
       END IF

       eig_qb(:) = 0.

       nslibd_qb = 0
       do i = 1,nbands_qb
          if((eigg(i).ge.sliceplot%e1s.and.nslibd_qb.lt.numbands
     &       .and.wann%l_bynumber).or.
     &       (eigg(i).ge.sliceplot%e1s.and.eigg(i).le.sliceplot%e2s.and.
     &        wann%l_byenergy).or.(i.ge.wann%band_min(jspin).and.
     &       (i.le.wann%band_max(jspin)).and.wann%l_byindex)) then
            nslibd_qb = nslibd_qb + 1
            eig_qb(nslibd_qb) = eigg(i)
            we_qb(nslibd_qb) = we_qb(i)
            if(noco%l_noco)then
               funbas =        nv_qb(1) + nlotot
               funbas = funbas+nv_qb(2) + nlotot
            else
               funbas = nv_qb(jspin_b) + nlotot      
            endif
            IF (zzMat%l_real) THEN
               do j = 1,funbas
                  zMat_qb%z_r(j,nslibd_qb) = zzMat%z_r(j,i)
               enddo
            ELSE
               do j = 1,funbas
                  zMat_qb%z_c(j,nslibd_qb) = zzMat%z_c(j,i)
               enddo
            END IF
          endif
        enddo

! check that eigenvectors and -values are identical if q=q+b
       if(iqpt.eq.qptibz_b .and. jspin.eq.jspin_b) then
          IF(zMat%l_real) THEN
             if(any(zMat%z_r.ne.zMat_qb%z_r)) 
     +          write(*,*)'z.ne.z_qb',iqpt,ikpt
          ELSE
             if(any(zMat%z_c.ne.zMat_qb%z_c)) 
     +          write(*,*)'z.ne.z_qb',iqpt,ikpt
          END IF
          if(any(eig.ne.eig_qb)) write(*,*)'eig.ne.eiq_qb',iqpt,ikpt
          if(nv(jspin).ne.nv_qb(jspin)) write(*,*)'nv!=nv_qb',iqpt,ikpt
       endif

! check that number of bands are the same at (k,q) and (k,q+b)
        if(nslibd.ne.nslibd_qb)
     >      write(*,*)'nslibd.ne.nslibd_qb',ikpt,iqpt,iqpt_b


        noccbd_qb = nslibd_qb
        nsfactor_b=cmplx(1.0,0.0)

!         call cpu_time(tt5)
         !if(l_p0) write(*,*)'nbnd_qb=',tt5-tt4

        allocate (acof_qb(noccbd_qb,0:lmd,atoms%nat),
     &            bcof_qb(noccbd_qb,0:lmd,atoms%nat),
     &            ccof_qb(-atoms%llod:atoms%llod,noccbd_qb,atoms%nlod,
     &            atoms%nat))

        acof_qb(:,:,:) = cmplx(0.,0.) ; bcof_qb(:,:,:) = cmplx(0.,0.)
        ccof_qb(:,:,:,:) = cmplx(0.,0.)

! construct the A,B,C coefficients of the wave function
! at the point (k,q+b) using previously read information

      ALLOCATE(lapw_qb%k1(SIZE(k1_qb,1),SIZE(k1_qb,2)))
      ALLOCATE(lapw_qb%k2(SIZE(k1_qb,1),SIZE(k1_qb,2)))
      ALLOCATE(lapw_qb%k3(SIZE(k1_qb,1),SIZE(k1_qb,2)))
      lapw_qb%k1 = k1_qb
      lapw_qb%k2 = k2_qb
      lapw_qb%k3 = k3_qb
      lapw_qb%nmat = nmat_qb
      lapw_qb%nv = nv_qb
      ! I think the other variables of lapw are not needed here.

!        CALL abcof(input,atoms,noccbd_qb,sym,cell,bkpt,lapw_qb,
!     +             noccbd_qb,usdus,noco,jspin_b,kveclo_qb,oneD,
!     +             acof_qb,bcof_qb,ccof_qb,zMat_qb)

      DEALLOCATE(lapw_qb%k1,lapw_qb%k2,lapw_qb%k3)

        call wann_abinv(
     >        atoms%ntype,atoms%nat,noccbd_qb,atoms%lmaxd,lmd,
     >        atoms%llod,atoms%nlod,atoms%ntype,atoms%neq,      
     >        noccbd_qb,atoms%lmax,atoms%nlo,atoms%llo,atoms%invsat,
     >        sym%invsatnr,
     >        bqpt,  ! does not enter            
     >        atoms%taual, ! does not enter
     X        acof_qb,bcof_qb,ccof_qb)

! check that A,B,C coefficients are the same if q+b = q
        if(l_gwf.and.(iqpt.eq.qptibz_b).and.(jspin.eq.jspin_b)) then
           if(any(acof_qb.ne.acof)) write(*,*)'acof',iqpt,ikpt
           if(any(bcof_qb.ne.bcof)) write(*,*)'bcof',iqpt,ikpt
           if(any(ccof_qb.ne.ccof)) write(*,*)'ccof',iqpt,ikpt
        endif  

!        call cpu_time(tt4)
        !if(l_p0) write(*,*)'abcof_qb=',tt4-tt5

        addnoco=0
        addnoco2=0
        if(noco%l_noco.and.(jspin.eq.2))then
           addnoco  = nv(1)   + nlotot
        endif
        if(noco%l_noco.and.(jspin_b.eq.2))then
           addnoco2 = nv_qb(1) + nlotot
        endif


! set up local->global transformation for overlaps
        do n=1,atoms%ntype
           if(wann%l_sgwf) then
              dalph = alph_i(n)-alphb_i(n)
              db1 = beta_i(n)/2.
              db2 = betab_i(n)/2.
           elseif(wann%l_socgwf) then
              dalph = phi_i-phib_i
              db1 = theta_i/2.
              db2 = thetab_i/2.
           endif
           coph = cos(dalph)    
           siph = sin(dalph)     
           phasfac = cmplx(coph,siph)
           phasfac2= cmplx(coph,-siph)

           if(l_p0 .and. dalph.ne.0.0) then
            write(*,*)'WARNING: include dalph in chi trafo!'
           endif

           if( (jspin.eq.1) .and. (jspin_b.eq.1) ) then ! uu
!              chi(n) = cos(db1)*cos(db2)*phasfac
!     >               + sin(db1)*sin(db2)*phasfac2
            chi(n) = cos(db2-db1)
           elseif( (jspin.eq.2) .and. (jspin_b.eq.2) ) then !dd
!              chi(n) = cos(db1)*cos(db2)*phasfac2
!     >               + sin(db1)*sin(db2)*phasfac
            chi(n) = cos(db2-db1)
           elseif( (jspin.eq.1) .and. (jspin_b.eq.2) ) then ! ud
!              chi(n) = sin(db1)*cos(db2)*phasfac2
!     >               - cos(db1)*sin(db2)*phasfac
            chi(n) = -sin(db2-db1)
           elseif( (jspin.eq.2) .and. (jspin_b.eq.1) ) then ! du
!              chi(n) = cos(db1)*sin(db2)*phasfac2
!     >               - sin(db1)*cos(db2)*phasfac
            chi(n) = sin(db2-db1)
           else
            stop 'problem setting up chi: jspin,jspin_b'
           endif
        enddo
        chi = conjg(chi)
        !chi = cmplx_1

        ! optional: disable chi transformation
        ! instead of computing overlap w.r.t. global frame
        ! only consider wave function overlaps in local frames
        if(l_nochi) then
         if(doublespin.lt.3) then
          chi = cmplx(1.0,0.0)
         else
          chi = cmplx(0.0,0.0)
         endif
        endif

        if((iqpt.eq.1).and.(ikpt.eq.1).and.(iqpt_b.eq.1)) then
           write(*,*)'dbs',doublespin,'chi',chi(1)
        endif

! muffin tin contribution to overlap is computed taking into account
! the spin-dependent phase exp(+/- tau*b/2) and the ujugaunt integrals
! calculated especially for the q-points before. then, q and q+b
! take the role of k and k+b and the same for G(q+b) and G(k+b)
        if(wann%l_sgwf) call wann_mmkb_sph(        
     >       nbnd,atoms%llod,nslibd,nslibd_qb,atoms%nlod,atoms%nat,
     >       atoms%ntype,lmd,atoms%jmtd,sign_q*atoms%taual/2.0,sym%nop,
     >       atoms%lmax,atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,
     >       acof,bcof,ccof,qptb_i,
     >       acof_qb,bcof_qb,ccof_qb,gb_q(:,iqpt_b,iqpt),qpt_i,
     >       ujug_q,ujdg_q,
     >       djug_q,djdg_q,ujulog_q,djulog_q,ulojug_q,ulojdg_q,
     >       ulojulog_q,qdiff,      
     >       nntot_q,chi,
     =       mmnk_q(:,:,iqpt_b,ikpt))
        if(wann%l_socgwf) call wann_mmkb_sph(        
     >       nbnd,atoms%llod,nslibd,nslibd_qb,atoms%nlod,atoms%nat,
     >       atoms%ntype,lmd,atoms%jmtd,
     >       zero_taual,sym%nop,atoms%lmax,                        
     >       atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,acof,bcof,ccof,
     >       (/ 0.0, phib_i/tpi, thetab_i/tpi /),
     >       acof_qb,bcof_qb,ccof_qb,gb_q(:,iqpt_b,iqpt),
     >       (/ 0.0, phi_i/tpi, theta_i/tpi /),
     >       ujug_q,ujdg_q,
     >       djug_q,djdg_q,ujulog_q,djulog_q,ulojug_q,ulojdg_q,
     >       ulojulog_q,qdiff,      
     >       nntot_q,chi,
     =       mmnk_q(:,:,iqpt_b,ikpt))

!        call cpu_time(tt5)
        !if(l_p0) write(*,*)'sph=',tt5-tt4


       if(((doublespin.ne.3).and.(doublespin.ne.4))
     >     .or.(.not.noco%l_noco)) then

          if(.not.noco%l_noco)then
             interchi=chi(1)
             vacchi=chi(1)
          else
             interchi=cmplx_1
             vacchi=cmplx_1
          endif
          
!          call cpu_time(tt4)
! interstitial contribution to overlap is computed using
! (-/+ 1)*G(q+b)/2 as G-vector connecting neighbors
! and lattice vectors G(k,q) (k1...) and G(k,q+b) (k1_qb...)
        if(wann%l_sgwf) call wann_mmkb_int(
     >         interchi,addnoco,addnoco2,
     >         DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,
     >         stars%ng3,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,
     >         k1_qb(:,jspin_b),k2_qb(:,jspin_b),k3_qb(:,jspin_b),
     >         nv_qb(jspin_b),zMat_qb,nslibd_qb,
     >         nbnd,
     >         stars%rgphs,stars%ustep,stars%ig,
     >         sign_q*gb_q(:,iqpt_b,iqpt)/2,     
     =         mmnk_q(:,:,iqpt_b,ikpt))                        
        if(wann%l_socgwf) call wann_mmkb_int(
     >         interchi,addnoco,addnoco2,
     >         DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,
     >         stars%ng3,k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >         nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,
     >         k1_qb(:,jspin_b),k2_qb(:,jspin_b),k3_qb(:,jspin_b),
     >         nv_qb(jspin_b),zMat_qb,nslibd_qb,
     >         nbnd,
     >         stars%rgphs,stars%ustep,stars%ig,(/ 0, 0, 0 /),     
     =         mmnk_q(:,:,iqpt_b,ikpt))!m_int(:,:,iqpt_b,ikpt))      

!        call cpu_time(tt5)
        !if(l_p0) write(*,*)'int=',tt5-tt4

! vacuum contribution in film calculation
        if (input%film .and. .not.oneD%odi%d1) then
           if(wann%l_sgwf) call wann_mmkb_vac(                  
     >              vacchi,noco%l_noco,nlotot,sign_q*2.*bkpt,
     >              nbnd,cell%z1,vacuum%nmzd,DIMENSION%nv2d,
     >              stars%mx1,stars%mx2,stars%mx3,
     >              stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,
     >              vacuum%delz,stars%ig2,cell%area,cell%bmat,
     >              cell%bbmat,evac(:,jspin),evac(:,jspin_b),
     >              sign_q*qpt_i/2.,
     >              sign_q*qptb_i/2.,
     >              vz(:,:,jspin2),vz(:,:,jspin2_b),
     >              nslibd,nslibd_qb,jspin,jspin_b,
     >              k1,k2,k3,k1_qb,k2_qb,k3_qb,
     >              wannierspin,DIMENSION%nvd,
     >              DIMENSION%nbasfcn,DIMENSION%neigd,zMat,zMat_qb,nv,
     >              nv_qb,cell%omtil,
     >              sign_q*gb_q(:,iqpt_b,iqpt)/2,
     =              mmnk_q(:,:,iqpt_b,ikpt))
           if(wann%l_socgwf) call wann_mmkb_vac(                    
     >              vacchi,noco%l_noco,nlotot,qpt_i,
     >              nbnd,cell%z1,vacuum%nmzd,DIMENSION%nv2d,
     >              stars%mx1,stars%mx2,stars%mx3,
     >              stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,
     >              vacuum%delz,stars%ig2,cell%area,cell%bmat,
     >              cell%bbmat,evac(:,jspin),evac(:,jspin_b),
     >              bqpt,bqpt,
     >              vz(:,:,jspin2),vz(:,:,jspin2_b),
     >              nslibd,nslibd_qb,jspin,jspin_b,
     >              k1,k2,k3,k1_qb,k2_qb,k3_qb,
     >              wannierspin,DIMENSION%nvd,
     >              DIMENSION%nbasfcn,DIMENSION%neigd,zMat,zMat_qb,nv,
     >              nv_qb,cell%omtil,
     >              (/ 0, 0, 0 /),
     =              mmnk_q(:,:,iqpt_b,ikpt))

! vacuum contribution in one-dimensional chain calculation where
! q-point plays role of k-point with proper prefactor of (-/+ 1/2).
! moreover, the k-point ikpt is treated like qss in the subroutine
! such that a correction for sign and factor has to appear as well.
! lattice vectors G(k,q) (k1...) and G(k,q+b) (k1_qb...) need to 
! be provided.
         elseif (oneD%odi%d1) then
           if(wann%l_sgwf) call wann_mmkb_od_vac(
     >           DIMENSION,oneD,vacuum,stars,cell,
     >           vacchi,noco%l_noco,nlotot,          
     >           nbnd,cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,
     >           stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,
     >           stars%ig,vacuum%nmzxy,
     >           vacuum%nmz,vacuum%delz,stars%ig2,oneD%odi%n2d,
     >           cell%bbmat,evac(1,jspin),evac(1,jspin_b),
     >           sign_q*qpt_i/2.,sign_q*qptb_i/2.,
     >           oneD%odi%M,oneD%odi%mb,
     >           vz(:,1,jspin2),vz(:,1,jspin2_b),oneD%odi,
     >           nslibd,nslibd_qb,jspin,jspin_b,
     >           k1,k2,k3,k1_qb,k2_qb,k3_qb,
     >           wannierspin,DIMENSION%nvd,cell%area,DIMENSION%nbasfcn,
     >           DIMENSION%neigd,zMat,zMat_qb,nv,nv_qb,stars%sk2,
     >           stars%phi2,cell%omtil,
     >           sign_q*gb_q(:,iqpt_b,iqpt)/2,sign_q*2.*bkpt, 
     >           .true.,sign_q,
     =           mmnk_q(:,:,iqpt_b,ikpt))
           if(wann%l_socgwf) call wann_mmkb_od_vac(   
     >           DIMENSION,oneD,vacuum,stars,cell,               
     >           vacchi,noco%l_noco,nlotot,          
     >           nbnd,cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,
     >           stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,
     >           stars%ig,vacuum%nmzxy,
     >           vacuum%nmz,vacuum%delz,stars%ig2,oneD%odi%n2d,
     >           cell%bbmat,evac(1,jspin),evac(1,jspin_b),
     >           bqpt,bqpt,
     >           oneD%odi%M,oneD%odi%mb,
     >           vz(:,1,jspin2),vz(:,1,jspin2_b),oneD%odi,
     >           nslibd,nslibd_qb,jspin,jspin_b,
     >           k1,k2,k3,k1_qb,k2_qb,k3_qb,
     >           wannierspin,DIMENSION%nvd,cell%area,DIMENSION%nbasfcn,
     >           DIMENSION%neigd,zMat,zMat_qb,nv,nv_qb,stars%sk2,
     >           stars%phi2,cell%omtil,
     >           (/ 0, 0, 0 /),qpt_i, 
     >           .false.,1,
     =           mmnk_q(:,:,iqpt_b,ikpt))
        endif!film resp. odi

!        call cpu_time(tt4)
        !if(l_p0) write(*,*)'film=',tt4-tt5

      endif!doublespin

        deallocate ( acof_qb,bcof_qb,ccof_qb )

!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'q-neig',iqpt_b,' total=',tt4-tt3
      t_qov=t_qov+tt4-tt3

      ENDDO !iqpt_b, q-neighbors
c**************************************************c
c              END Q-NEIGHBOR LOOP                 c
c**************************************************c     

      deallocate ( we_qb )

      endif    !   if wann%l_matrixmmn = true

c$$$  3     continue !goto jump no longer needed?

        if(.not.wann%l_bzsym)oper=0
        if(.not.wann%l_plot_symm.or.oper.eq.1)then
        if (wann%l_wann_plot .and. 
     >        (doublespin.eq.1 .or. doublespin.eq.2)) then

         addnoco=0
         if(noco%l_noco.and.(jspin.eq.2))then
            addnoco=nv(1)+nlotot
         endif

         if (sliceplot%slice) then
          if (ikpt.eq.sliceplot%kk) then

           write (6,*) 'nnne=',sliceplot%nnne
           write (6,*) 'eig(nnne)=',eig(sliceplot%nnne)
           write (6,*) 'we(nnne)=',we(sliceplot%nnne)

           call wann_plot(
     >        DIMENSION,oneD,vacuum,stars,cell,atoms,
     >        DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,
     >        vacuum%nmzxyd,
     >        stars%ng2,sphhar%ntypsd,atoms%ntype,atoms%lmaxd,
     >        atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,atoms%neq,
     >        stars%ng3,vacuum%nvac,vacuum%nmz,vacuum%nmzxy,stars%ng2,
     >        sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,
     >        sym%symor,sym%invs,sym%invs2,cell%z1,vacuum%delz,
     >        atoms%ngopr,atoms%ntypsy,atoms%jri,atoms%pos,atoms%zatom,
     >        atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,
     >        cell%amat,cell%bmat,cell%bbmat,ikpt,sliceplot%nnne,
     >        sliceplot%kk,DIMENSION%nvd,atoms%nlod,atoms%llod,
     >        nv(jspin),lmd,bkpt,cell%omtil,atoms%nlo,atoms%llo,
     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
     >        vz(:,:,jspin2),
     >        nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,
     >        ff(:,:,:,:,jspin),
     >        gg(:,:,:,:,jspin),flo,acof,bcof,ccof,zMat,
     >        stars%mx1,stars%mx2,stars%mx3,stars%ig,stars%ig2,
     >        stars%sk2,stars%phi2,
     >        noco%l_noco,noco%l_ss,qpt_i,
     >        addnoco,get_index_kq(ikpt,iqpt,fullnkpts),wann%l_sgwf)

          endif
         else ! not sliceplot%slice

          call wann_plot(
     >        DIMENSION,oneD,vacuum,stars,cell,atoms,
     >        DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,
     >        vacuum%nmzxyd,
     >        stars%ng2,sphhar%ntypsd,atoms%ntype,atoms%lmaxd,
     >        atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,atoms%neq,
     >        stars%ng3,vacuum%nvac,vacuum%nmz,vacuum%nmzxy,stars%ng2,
     >        sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,
     >        sym%symor,sym%invs,sym%invs2,cell%z1,vacuum%delz,
     >        atoms%ngopr,atoms%ntypsy,atoms%jri,atoms%pos,atoms%zatom,
     >        atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,
     >        cell%amat,cell%bmat,cell%bbmat,ikpt,sliceplot%nnne,
     >        sliceplot%kk,DIMENSION%nvd,atoms%nlod,atoms%llod,
     >        nv(jspin),lmd,bkpt,cell%omtil,atoms%nlo,atoms%llo,
     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
     >         vz(:,:,jspin2),
     >        nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,
     >        ff(:,:,:,:,jspin),
     >        gg(:,:,:,:,jspin),flo,acof,bcof,ccof,zMat,
     >        stars%mx1,stars%mx2,stars%mx3,stars%ig,stars%ig2,
     >        stars%sk2,stars%phi2,
     >        noco%l_noco,noco%l_ss,qpt_i,
     >        addnoco,get_index_kq(ikpt,iqpt,fullnkpts),wann%l_sgwf)

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
              kplotoper=sym%invtab(plotoper)
              call wann_plot_symm(jspin,sym%mrot(:,:,kplotoper),ikpt,
     &                            kplot,l_conjugate)
            endif
           enddo
          endif
       
         endif
        endif
        endif !wann%l_plot_symm
      deallocate ( acof,bcof,ccof,we,eigg )

      if(wann%l_projmethod.or.wann%l_bestproj)then
         call wann_projmethod(
     >            fullnkpts,
     >            wann%l_projmethod,wann%l_bestproj,
     >            ikpt,nwfs,nslibd,amn,eig,
     <            psiw,hwfr)
      endif ! projmethod

      endif   ! loop by processors

!      call cpu_time(tt3)
      !if(l_p0) write(*,*)'k-loop',ikpt,' total=',tt3-tt2
10    continue ! end of cycle by the k-points


      if(wann%l_matrixmmn)
     &      deallocate(ujug,ujdg,djug,djdg,
     &                 ujulog,djulog,ulojulog,ulojug,ulojdg)
      if(wann%l_matrixmmn.AND.l_gwf)
     &      deallocate(ujug_q,ujdg_q,djug_q,
     &                 djdg_q,ujulog_q,djulog_q,ulojulog_q,ulojug_q,
     &                 ulojdg_q)

#ifdef CPP_MPI
      call MPI_BARRIER(mpi%mpi_comm,ierr)
#endif


c******************************************************
c     Write down the projections.
c******************************************************
  5   continue

      if(doublespin.eq.3 .or. doublespin.eq.4) goto 912

      if(wann%l_nabla)then
        hescale=tpi*condquant/bohrtocm/cell%omtil
        hescale=sqrt(hescale)
        nablamat=nablamat*hescale
        call wann_write_nabla(
     >           mpi%mpi_comm,l_p0,spin12(jspin)//'.nabl',
     >           'Matrix elements of nabla operator',
     >           nbnd,fullnkpts,nbnd,
     >           mpi%irank,mpi%isize,
     <           nablamat)
      endif

      if(wann%l_soctomom)then
        hescale=tpi*condquant/bohrtocm/cell%omtil
        hescale=sqrt(hescale)
        soctomom=soctomom*hescale
        call wann_write_nabla(
     >           mpi%mpi_comm,l_p0,spin12(jspin)//'.stm',
     >           'Matrix elements of stm operator',
     >           nbnd,fullnkpts,nbnd,
     >           mpi%irank,mpi%isize,
     <           soctomom)
      endif

      if(wann%l_surfcurr)then
        surfcurr = conjg( surfcurr/cmplx(0.0,2.0) )
        hescale=tpi*condquant/bohrtocm/cell%omtil
        hescale=sqrt(hescale)
        surfcurr=surfcurr*hescale
        call wann_write_nabla(
     >           mpi%mpi_comm,l_p0,spin12(jspin)//'.surfcurr',
     >           'Surface currents',
     >           nbnd,fullnkpts,nbnd,
     >           mpi%irank,mpi%isize,
     <           surfcurr)
      endif

      if((noco%l_soc.or.noco%l_noco).and.wann%l_mmn0)then
        call wann_write_amn(
     >           mpi%mpi_comm,
     >           l_p0,spin12(jspin)//trim(fending)//'.socmmn0',
     >           'Overlaps of the wavefunct. at the same kpoint',
     >           nbnd,fullnkpts,nbnd,
     >           mpi%irank,mpi%isize,.false.,
     <           mmn,.false.)
      endif !noco%l_soc and l_mmn0  

      if(wann%l_orbcomp)then
         num_angl=9
         if(wann%l_oc_f)num_angl=16
         call wann_write_matrix5(
     >            mpi%mpi_comm,l_p0,spin12(jspin)//'.orbcomp',
     >            'angular components',
     >            nbnd,nbnd,
     >            num_angl,wann%oc_num_orbs,fullnkpts,
     >            mpi%irank,mpi%isize,
     <            orbcomp)
      endif

      if((noco%l_soc.or.noco%l_noco) .and. (doublespin.eq.1)) then
         if(wann%l_mmn0) socmmn(:,:,:)=mmn(:,:,:)
         goto 912
      endif

      if(noco%l_soc.or.noco%l_noco) then
         jspin2=1
         if(wann%l_mmn0)      mmn(:,:,:)=socmmn(:,:,:)+mmn(:,:,:)
         if(wann%l_mmn0)      deallocate(socmmn)
      endif

!      call cpu_time(tt3)
      if (wann%l_matrixamn)then
         call wann_write_amn(
     >           mpi%mpi_comm,
     >           l_p0,spin12(jspin2)//trim(fending)//'.amn',
     >           'Overlaps of the wavefunct. with the trial orbitals',
     >           nbnd,fullnkpts,nwfs,
     >           mpi%irank,mpi%isize,.false.,
     <           amn(:,:,:),wann%l_unformatted)
      endif !wann%l_matrixamn
!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'write amn=',tt4-tt3

      if(wann%l_anglmom)then
        call wann_write_matrix4(
     >          mpi%mpi_comm,
     >          l_p0,spin12(jspin2)//trim(fending)//'.anglmom',
     >          'Matrix elements of angular momentum',
     >          nbnd,nbnd,3,fullnkpts,
     >          mpi%irank,mpi%isize,
     <          anglmom)
      endif

      if (l_proj) then
c**************************************************************
c            for projmethod: write down WF1.umn
c*************************************************************
       if((wann%l_projmethod.or.wann%l_bestproj))then
         call wann_write_amn(
     >           mpi%mpi_comm,l_p0,spin12(jspin2)//'.umn',
     >           'transformation to first guess Wannier functions',
     >           nbnd,fullnkpts,nwfs,
     >           mpi%irank,mpi%isize,.false.,
     <           psiw,.false.)
#ifdef CPP_MPI
         allocate( hwfr2(size(hwfr,1),size(hwfr,2)) )
         length=nwfs*nwfs
         CALL MPI_REDUCE(
     +           hwfr,hwfr2,length,
     +           CPP_MPI_COMPLEX,MPI_SUM,0,
     +           mpi%mpi_comm,ierr)
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
     >            mpi%mpi_comm,
     >            l_p0,spin12(jspin2)//trim(fending)//'.mmn0',
     >            'Overlaps of the wavefunct. at the same kpoint',
     >            nbnd,fullnkpts,nbnd,
     >            mpi%irank,mpi%isize,.false.,
     <            mmn,.false.)
      endif !wann%l_mmn0  

      
c*****************************************************
c.....write down the matrix M^{k,b}_{mn} 
c*****************************************************

! 912   continue
!      call cpu_time(tt3)
      if(wann%l_matrixmmn.and.(.not.wann%l_skipkov))then
         call wann_write_mmnk(
     >            mpi%mpi_comm,jspin2,l_p0,fullnkpts,nntot,wann,
     >            maptopair,pair_to_do,nbnd,bpt,gb,
     >            mpi%isize,mpi%irank,fending,     
     <            mmnk,wann%l_unformatted)
      endif !wann%l_matrixmmn
!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'write mmnk=',tt4-tt3

 912  continue


      call cpu_time(tt3)
      if(l_gwf .and. wann%l_matrixmmn) then
!      mmnk_q = mmnk_q + m_sph+m_int+m_vac
      if(doublespin.eq.doublespin_max) then
         write(fname,'("param_",i4.4,".mmn")')iqpt
         call wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,
     >                         nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,
     >                         mpi%isize,mpi%irank,fname,mmnk_q,
     >                         wann%l_unformatted)
      endif

      if(.false.) then
         write(fname,'("param_",i4.4,"_",i1,".mmn")')iqpt,doublespin
         call wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,
     >                         nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,
     >                         mpi%isize,mpi%irank,fname,
     >                         m_sph+m_int+m_vac,wann%l_unformatted)

         write(fname,'("param_",i4.4,"_",i1,"_int.mmn")')iqpt,doublespin
         call wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,
     >                         nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,
     >                         mpi%isize,mpi%irank,fname,m_int,
     >                         wann%l_unformatted)
         write(fname,'("param_",i4.4,"_",i1,"_sph.mmn")')iqpt,doublespin
         call wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,
     >                         nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,
     >                         mpi%isize,mpi%irank,fname,m_sph,
     >                         wann%l_unformatted)
         write(fname,'("param_",i4.4,"_",i1,"_vac.mmn")')iqpt,doublespin
         call wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,
     >                         nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,
     >                         mpi%isize,mpi%irank,fname,m_vac,
     >                         wann%l_unformatted)
      endif
!      m_int = cmplx(0.,0.)
!      m_sph = cmplx(0.,0.)
!      m_vac = cmplx(0.,0.)

      endif

!      call cpu_time(tt4)
      !if(l_p0) write(*,*)'write mmnq=',tt4-tt3

      if( allocated (nablamat) ) deallocate( nablamat )
      if( allocated (soctomom) ) deallocate( soctomom )

      if( allocated (surfcurr) ) deallocate( surfcurr )
      if( allocated( mmn ) ) deallocate(mmn)
      if( allocated( amn ) ) then
         if(.not.(noco%l_soc.or.noco%l_noco))deallocate(amn)
      endif
      if ( allocated (psiw) ) deallocate ( psiw )
      if (wann%l_matrixmmn) then
         if(.not.(noco%l_soc.or.noco%l_noco))deallocate (mmnk)
      endif  
      if (wann%l_anglmom .and. allocated(anglmom))then
         if(.not.(noco%l_soc.or.noco%l_noco))deallocate (anglmom)
      endif
      deallocate(flo)

      if(.not.noco%l_noco)nrec=nrec+nkpts

c#ifdef CPP_MPI
c      call MPI_BARRIER(mpi%mpi_comm,ierr)
c#endif

110   continue ! end of cycle by spins


#ifdef CPP_MPI
      call MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

      ! close eig files
      IF (l_gwf) THEN
!         CALL close_eig(eig_id)
         IF(wann%l_matrixmmn)THEN
            DO iqpt_b=1,nntot_q
!               CALL close_eig(innerEig_idList(iqpt_b))
            ENDDO
         ENDIF
      ENDIF

      if (wann%l_matrixmmn.AND.allocated(mmnk))then 
         deallocate ( mmnk )
      endif

      if(allocated(mmnk_q)) deallocate(mmnk_q)
      if(allocated(m_int)) deallocate(m_int)
      if(allocated(m_sph)) deallocate(m_sph)
      if(allocated(m_vac)) deallocate(m_vac)

      if ((wann%l_projmethod.or.wann%l_bestproj.or.wann%l_matrixamn)
     >       .AND.allocated(amn))then
         deallocate ( amn )
      endif

      if(wann%l_anglmom.AND.allocated(anglmom))then
         deallocate ( anglmom )  
      endif

      if ((wann%l_projmethod.or.wann%l_bestproj)
     >       .AND.allocated(hwfr)) then
         deallocate ( hwfr )
      endif

!      call cpu_time(tt2)
      !if(l_p0) write(*,*)'q-loop',iqpt,' total=',tt2-tt1
      t_it = tt2-tt1
      if(l_p0) then
!         write(*,*)t_kov,t_qov,t_it
      endif

      DEALLOCATE(innerEig_idList)

314   continue ! iqpt, q-points
c************************************************c
c               END Q LOOP                       c
c************************************************c

      if(allocated(pair_to_do))deallocate(pair_to_do,maptopair)

      deallocate ( vr,vz,kveclo,nv,k1,k2,k3 )
      deallocate ( ff,gg )
      if (wann%l_matrixmmn) deallocate ( kveclo_b,nv_b,k1_b,k2_b,k3_b )
      if (wann%l_matrixmmn.AND.l_gwf) 
     >     deallocate ( kveclo_qb,nv_qb,k1_qb,k2_qb,k3_qb ) 
      if (wann%l_bzsym)  deallocate(irreduc,mapkoper)
      if (wann%l_bzsym.AND.l_gwf)  deallocate(irreduc_q,mapqoper)
      if(allocated(pair_to_do_q))
     >     deallocate(pair_to_do_q,maptopair_q)

      if (allocated(kdiff)) deallocate ( kdiff )
      if (allocated(qdiff)) deallocate(qdiff,zero_qdiff)


9110  continue  ! jump for l_finishgwf=T

! correct for previously introduced factor of 2 in the
! G-vectors connecting neighbors across the BZ boundary
      if(wann%l_sgwf) gb_q = gb_q/2 
      if(wann%l_socgwf) gb_q = gb_q/2

! set up input files for wannier90 --> HDWFs
      if(l_p0.and.l_gwf.and.(wann%l_matrixmmn.or.wann%l_matrixamn)
     >   .and.(.not.wann%l_skipkov)) then
         call wann_gwf_commat(fullnkpts,nntot,bpt,fullnqpts,
     >                        nntot_q,bpt_q,gb,gb_q,wann%aux_latt_const,
     >                        wann%l_unformatted,wann%l_matrixamn,
     >                        wann%l_matrixmmn,wann%l_dim,
     >                        wann%nparampts,wann%param_vec/2.0)
      endif

      if(l_p0.and.l_gwf.and.wann%l_anglmom) then
       call wann_gwf_anglmom(fullnkpts,fullnqpts,wann%l_unformatted)
      endif

      if (wann%l_matrixmmn) THEN
            deallocate (gb,bpt)
            IF(l_gwf) deallocate (gb_q,bpt_q)
      endif


 1911 continue

      if(allocated(chi)) deallocate(chi)

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

      write(*,*)'time_total wannier=',time_total

      call wann_postproc(
     >         DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,
     >         lapw,oneD,noco,cell,vTot,eig_id,l_real,                !eig_id is used here after closing the files?!
     >         mpi%mpi_comm,wann,l_p0,input%film,input%jspins,stars%ng2,
     >         atoms%nat,atoms%pos,cell%amat,cell%bmat,atoms%ntype,
     >         atoms%neq,atoms%zatom,cell%omtil,noco%l_soc,noco%l_noco,
     >         DIMENSION%neigd,fullnkpts,oneD%odi,l_proj,atoms%l_dulo,
     >         noco%l_ss,atoms%lmaxd,atoms%ntype,sym%nop,DIMENSION%nvd,
     >         DIMENSION%jspd,DIMENSION%nbasfcn,atoms%llod,atoms%nlod,
     >         atoms%nlo,atoms%llo,atoms%lapw_l,sym%invtab,
     >         sym%mrot,atoms%ngopr,atoms%lmax,atoms%invsat,
     >         sym%invsatnr,nkpt,atoms%taual,atoms%rmt,cell%bbmat,
     >         noco%alph,noco%beta,noco%qss,stars%sk2,stars%phi2,
     >         oneD%ods,
     >         mpi%irank,mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,
     >         atoms%jmtd,sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,
     >         sym%invs2,sphhar%nlh,atoms%jri,atoms%ntypsy,
     >         sphhar%ntypsd,nkpt,atoms%dx,atoms%rmsh,
     >         sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,stars%ustep,
     >         stars%mx1,stars%mx2,stars%mx3,sym%nop2,
     >         stars%ig,stars%rgphs,
     >         sliceplot%slice,sliceplot%kk,sliceplot%nnne,
     >         DIMENSION%nv2d,vacuum%nmzxy,vacuum%nmz,stars%ng2,
     >         vacuum%delz,cell%area,cell%z1,cell%volint,stars%ig2,
     >         sym%tau,sym%symor,results%ef,wann%l_sgwf,fullnqpts)

#ifdef CPP_MPI
      call MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

      IF(.not.wann%l_ldauwan) THEN
         DO pc = 1, wann%nparampts
            CALL close_eig(eig_idList(pc))
         END DO

         CALL juDFT_end("wannier good",mpi%irank)
      END IF

      END SUBROUTINE wannier


      END MODULE m_wannier
