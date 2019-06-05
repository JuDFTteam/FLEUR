!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_wannier
  USE m_juDFT
CONTAINS
  SUBROUTINE wannier(&
       DIMENSION,mpi,input,kpts,sym,atoms,stars,vacuum,sphhar,oneD,&
       wann,noco,cell,enpara,banddos,sliceplot,vTot,results,&
       eig_idList,l_real,nkpt)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     Makes necessary for the construction of the wannier functions
    !     (W90: Yates, Mostofi, Marzari, Souza, Vanderbilt '06 f90 code)
    !     ab initio preliminaries: constructs the overlaps of the periodic
    !     parts of the wavefunctions and the projections of the 
    !     wavefunctions
    !     onto a set of starting wfs, i.e. atomic-like orbitals.
    !                                                            YM 06
    !     Mmn(k,b) = <u_{nk}|u_{m(k+b)}>, u being a periodic part
    !                        of the wavefunction psi_nk
    !     A_mn^k = <psi_mk|g_n>, where g_n is a trial orbital
    !     which are written into the files 'WF1.mmn' and 'WF1.amn'
    !           Marzari Vanderbilt PRB 56,12847(1997)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     Parallelization, Optionals, Symmetry, Noco&Soc:
    !     Frank Freimuth
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     The routine reads the bkpts file, which contains the following
    !     information:
    !     1st line: nntot (INT) - number of the nearest neighbors for
    !                             each k-point in the MP mesh
    !     2-nkpts*nntot lines containing 5 integers i1,i2,i3,i4,i5:
    !     i1 - the number of the k-point in the kpts file
    !     i2 - number of the k-point, which is a periodic image of
    !          k+b in the 1st BZ
    !     i3-i5 - coordinates of the G-vector, connecting k-point
    !             i2 with the actual k+b k-point
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     In general, the number of bands for each k-poin t is
    !     different: N(k), and also differs from the number of bands
    !     we are interested in: N (for instance 5 d-bands of Cu among
    !     the 6 s- and d-bands). While matrices Mmn(k) are square
    !     for each k-point, matrices Mmn(k,b) can be made so after
    !     defining the maximum number of bands max(N(k)).
    !     The matrix Amn is non-diagonal by default (N(k)*N).
    !ccccc ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     Total number of wannier functions: nwfs
    !     sliceplot%e1s,sliceplot%e2s: lower and upper boundaries of the energy window:
    !     Needed for sorting by number and sorting by energy.
    !     Not needed for sorting by index.
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     Extension to case of higher-dimensional Wannier functions        
    !     according to the formalism in PRB 91, 184413 (2015)
    !     Jan-Philipp Hanke                                         
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    USE m_types
    USE m_wann_mmnk_symm
    USE m_wann_rw_eig
    USE m_abcof
    USE m_radfun
    USE m_radflo
    USE m_cdnread
    USE m_constants
    USE m_wann_mmk0_od_vac
    USE m_wann_mmkb_od_vac
    USE m_wann_mmk0_vac
    USE m_wann_mmkb_vac
    USE m_wann_updown
    USE m_wann_mmk0_sph
    USE m_wann_ujugaunt
    USE m_wann_mmkb_sph
    USE m_wann_projmethod
    USE m_wann_amn
    USE m_wann_abinv
    USE m_wann_kptsrotate
    USE m_wann_plot
    USE m_wann_read_inp
    USE m_wann_plot_symm
    USE m_wann_mmkb_int
    USE m_wann_postproc
    USE m_matmul,ONLY : matmul3,matmul3r
    USE m_wann_write_mmnk
    USE m_wann_write_amn
    USE m_wann_write_nabla
    USE m_vsoc
    USE m_wann_write_matrix4
    USE m_wann_write_matrix5
    USE m_wann_orbcomp
    USE m_wann_anglmom
#ifdef CPP_TOPO
    USE m_wann_surfcurr
    USE m_wann_surfcurr_int2
    USE m_wann_nabla
    USE m_wann_nabla_vac
    USE m_wann_soc_to_mom
#endif
    USE m_wann_gwf_tools, ONLY : get_index_kq, gwf_plottemplate
    USE m_wann_gwf_commat
    USE m_wann_gwf_anglmom
    USE m_wann_write_mmnk2
    USE m_wann_uHu
    USE m_wann_uHu_dmi
    USE m_eig66_io

    IMPLICIT NONE
#include "cpp_double.h"
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER ierr(3)
    INTEGER cpu_index
    INTEGER stt(MPI_STATUS_SIZE)
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

    LOGICAL, INTENT (in) :: l_real
    INTEGER, INTENT (in) :: nkpt
    INTEGER, INTENT (IN) :: eig_idList(wann%nparampts)

    !ccccccccccccccccc   local variables   cccccccccccccccccccc
    INTEGER :: lmd,n,nmat,iter,ikpt,ikpt_b, pc
    INTEGER :: addnoco,loplod,addnoco2,igvm2,eig_id
    INTEGER :: noccbd,noccbd_b,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
    INTEGER :: jsp_start,jsp_end,nrec,nrec1,nrec_b,nbands,nbands_b
    INTEGER :: nodeu,noded,n_size,na,n_rank,nbnd,numbands
    INTEGER :: i1,i2,i3,in,lda
    INTEGER :: n_bands(0:DIMENSION%neigd),nslibd,nslibd_b
    CHARACTER(len=8) :: dop,iop,name(10)
    REAL    :: wronk,phase
    COMPLEX :: c_phase
    REAL    :: eig(DIMENSION%neigd),eig_b(DIMENSION%neigd)
    REAL    :: efermi
    LOGICAL :: l_p0,l_bkpts,l_proj,l_amn,l_mmn
!!! energy window boundaries
    INTEGER, ALLOCATABLE :: innerEig_idList(:)
    REAL,    ALLOCATABLE :: we(:),we_b(:)

    REAL,    ALLOCATABLE :: eigg(:)
    REAL kpoints(nkpt)
!!! a and b coeff. constructed for each k-point
    COMPLEX, ALLOCATABLE :: acof(:,:,:),acof_b(:,:,:)
    COMPLEX, ALLOCATABLE :: bcof(:,:,:),bcof_b(:,:,:)
    COMPLEX, ALLOCATABLE :: ccof(:,:,:,:),ccof_b(:,:,:,:)
!!! the parameters for the number of wfs
    INTEGER :: nwfs
!!! the potential in the spheres and the vacuum
    REAL, ALLOCATABLE :: vr(:,:,:),vz(:,:,:)
!!! auxiliary potentials
    COMPLEX, ALLOCATABLE :: vpw(:,:)
!!! bkpts data
    INTEGER nntot,ikpt_help
    INTEGER, ALLOCATABLE :: gb(:,:,:),bpt(:,:)
!!! radial wavefunctions in the muffin-tins and more ...
    REAL,    ALLOCATABLE :: flo(:,:,:,:,:),vso(:,:,:)
    REAL,    ALLOCATABLE :: ff(:,:,:,:,:),gg(:,:,:,:,:)

    REAL     :: uuilon(atoms%nlod,atoms%ntype)
    REAL     :: duilon(atoms%nlod,atoms%ntype)
    REAL     :: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
!!! the Mmn matrices
    COMPLEX, ALLOCATABLE :: mmnk(:,:,:,:),mmn(:,:,:)           
    COMPLEX, ALLOCATABLE :: amn(:,:,:),nablamat(:,:,:,:)        
    COMPLEX, ALLOCATABLE :: soctomom(:,:,:,:)
    COMPLEX, ALLOCATABLE :: surfcurr(:,:,:,:)
    COMPLEX, ALLOCATABLE :: socmmn(:,:,:)
    COMPLEX, ALLOCATABLE :: a(:)
    COMPLEX, ALLOCATABLE :: psiw(:,:,:)
    COMPLEX, ALLOCATABLE :: anglmom(:,:,:,:)
    COMPLEX, ALLOCATABLE :: orbcomp(:,:,:,:,:)
    !..wf-hamiltonian in real space (hopping in the same unit cell)
    COMPLEX, ALLOCATABLE :: hwfr(:,:),hwfr2(:,:)
    !      real, allocatable :: ei(:)
    COMPLEX, ALLOCATABLE :: work(:)
    REAL,ALLOCATABLE::centers(:,:,:)
    LOGICAL :: l_file
    LOGICAL :: l_amn2, l_conjugate
    CHARACTER(len=3) :: spin12(2)
    DATA   spin12/'WF1' , 'WF2'/
    CHARACTER(len=30)  :: task
    INTEGER,ALLOCATABLE::irreduc(:)
    INTEGER,ALLOCATABLE::mapkoper(:)
    INTEGER :: fullnkpts,kpt,kptibz,kptibz_b,j1,j2,j3,oper,oper_b,k
    REAL :: bkrot(3),dirfacs(3)
    INTEGER :: ios,kplot,kplotoper,plotoper,gfcut
    COMPLEX :: phasust
    INTEGER,ALLOCATABLE::pair_to_do(:,:)
    INTEGER :: ngopr1(atoms%nat)
    INTEGER,ALLOCATABLE::maptopair(:,:,:)
    INTEGER :: wannierspin,jspin2,jspin7,jspin2_b
    REAL, ALLOCATABLE :: rwork(:)
    REAL,ALLOCATABLE::kdiff(:,:)
    INTEGER,ALLOCATABLE :: shiftkpt(:,:)
    INTEGER :: unigrid(6),gfthick
    COMPLEX,ALLOCATABLE::ujug(:,:,:,:),ujdg(:,:,:,:)
    COMPLEX,ALLOCATABLE::djug(:,:,:,:),djdg(:,:,:,:)
    COMPLEX,ALLOCATABLE::ujulog(:,:,:,:,:)
    COMPLEX,ALLOCATABLE::djulog(:,:,:,:,:)
    COMPLEX,ALLOCATABLE::ulojug(:,:,:,:,:)
    COMPLEX,ALLOCATABLE::ulojdg(:,:,:,:,:)
    COMPLEX,ALLOCATABLE::ulojulog(:,:,:,:,:,:)
    INTEGER :: n_start,n_end,mlotot,mlolotot,err
    INTEGER :: mlot_d,mlolot_d,ilo,dir,length
    CHARACTER(len=2) :: spin012(0:2)
    DATA spin012/'  ', '.1', '.2'/
    CHARACTER(len=6) :: filename
    REAL :: arg,hescale
    COMPLEX :: nsfactor,nsfactor_b,VALUE
    REAL :: b1(3),b2(3)
    REAL,PARAMETER :: bohrtocm=0.529177e-8
    REAL,PARAMETER :: condquant=7.7480917e-5
    INTEGER :: npotmatfile,ig3,maxvac,irec,imz,ivac,ipot
    LOGICAL :: l_orbcompinp
    INTEGER :: num_angl
    COMPLEX,ALLOCATABLE :: vxy(:,:,:)


    !---->gwf

    ! FURTHER VARIABLES
    REAL :: qpt_i(3),qptb_i(3)
    REAL :: alph_i(atoms%ntype),alphb_i(atoms%ntype)
    REAL :: beta_i(atoms%ntype),betab_i(atoms%ntype)
    REAL :: theta_i, thetab_i, phi_i, phib_i
    REAL :: dalph,db1,db2,coph,siph
    REAL :: zero_taual(3,atoms%nat),bqpt(3)
    REAL :: eig_qb(DIMENSION%neigd)

    REAL,ALLOCATABLE :: qdiff(:,:), we_qb(:)                
    REAL,ALLOCATABLE :: energies(:,:,:)  
    REAL,ALLOCATABLE :: zero_qdiff(:,:)


    INTEGER,ALLOCATABLE :: irreduc_q(:),mapqoper(:)        
    INTEGER,ALLOCATABLE :: shiftqpt(:,:),pair_to_do_q(:,:)  
    INTEGER,ALLOCATABLE :: maptopair_q(:,:,:)              
    INTEGER,ALLOCATABLE :: gb_q(:,:,:),bpt_q(:,:)         

    INTEGER :: nntot_q = 1                               
    INTEGER :: fullnqpts = 1                                
    INTEGER :: funit_start = 5000
    INTEGER :: qptibz, qptibz_b, oper_q, oper_qb
    INTEGER :: qpt,iqpt_help, iqpt, iqpt_b
    INTEGER :: nbands_qb, nmat_qb, nslibd_qb, noccbd_qb
    INTEGER :: sign_q = 1,band_help
    INTEGER :: doublespin,jspin_b,jspin3,jspin4,jspin5
    INTEGER :: doublespin_max,nrec5
    INTEGER :: count_i,count_j
    INTEGER :: n1,n2,ii,jj

    COMPLEX :: interchi,vacchi,amnchi
    COMPLEX :: phasfac,phasfac2,cmplx_1                                 

    COMPLEX,ALLOCATABLE :: chi(:)
    COMPLEX,ALLOCATABLE :: acof_qb(:,:,:)                 
    COMPLEX,ALLOCATABLE :: bcof_qb(:,:,:)                   
    COMPLEX,ALLOCATABLE :: ccof_qb(:,:,:,:)                 
    COMPLEX,ALLOCATABLE :: mmnk_q(:,:,:,:)                 
    COMPLEX,ALLOCATABLE :: m_int(:,:,:,:)
    COMPLEX,ALLOCATABLE :: m_sph(:,:,:,:)
    COMPLEX,ALLOCATABLE :: m_vac(:,:,:,:)
    COMPLEX,ALLOCATABLE :: ujug_q(:,:,:,:),ujdg_q(:,:,:,:) 
    COMPLEX,ALLOCATABLE :: djug_q(:,:,:,:),djdg_q(:,:,:,:) 
    COMPLEX,ALLOCATABLE :: ujulog_q(:,:,:,:,:)              
    COMPLEX,ALLOCATABLE :: djulog_q(:,:,:,:,:)             
    COMPLEX,ALLOCATABLE :: ulojug_q(:,:,:,:,:)             
    COMPLEX,ALLOCATABLE :: ulojdg_q(:,:,:,:,:)             
    COMPLEX,ALLOCATABLE :: ulojulog_q(:,:,:,:,:,:)         

    CHARACTER(len=30) fname

    LOGICAL :: l_bqpts,l_gwf,l_nochi

    TYPE(t_usdus) :: usdus
    TYPE(t_mat)   :: zMat, zzMat, zMat_b, zMat_qb
    TYPE(t_lapw)  :: lapw, lapw_b, lapw_qb
    TYPE(t_wann)  :: wannTemp

    eig_id = eig_idList(1)

    !----<gwf



    !-----initializations
    ngopr1(:)=1
    zero_taual = 0.0

    hescale=sqrt(tpi_const*condquant/bohrtocm/cell%omtil)

    cmplx_1 = CMPLX(1.0,0.0)

    CALL timestart("Wannier total")

    l_p0 = .FALSE.
    IF (mpi%irank.EQ.0) l_p0 = .TRUE.

    lmd = atoms%lmaxd*(atoms%lmaxd+2)

!!!   should be changed in case the windows are really used
    nkpts = nkpt

    ! do we have to construct GWF ?
    l_gwf = .FALSE.
    l_gwf = wann%l_sgwf.OR.wann%l_socgwf 

    l_nochi = .FALSE.
    INQUIRE(file='l_nochi',exist=l_nochi)
    IF(l_gwf.AND.l_p0) WRITE(*,*)'disable chi trafo: ',l_nochi

    IF(l_gwf.AND.l_p0) CALL gwf_plottemplate()
    ALLOCATE( chi(atoms%ntype) )

    !-----read the input file to determine what to do
    wann%atomlist_num=atoms%nat
    wann%oc_num_orbs=atoms%nat

    CALL wann_read_inp(input,l_p0,wann)


    !-----input file for orbital decomposition
    IF(wann%l_orbcomp.OR.wann%l_orbcomprs)THEN
       INQUIRE(file='orbcomp_inp',exist=l_orbcompinp)
       IF(l_orbcompinp)THEN
          OPEN(159,file='orbcomp_inp')
          READ(159,*)wann%oc_num_orbs,wann%l_oc_f
          ALLOCATE(wann%oc_orbs(wann%oc_num_orbs))
          DO n=1,wann%oc_num_orbs
             READ(159,*)wann%oc_orbs(n)
          ENDDO
          CLOSE(159)
       ELSE !default is all atoms including f
          wann%oc_num_orbs=atoms%nat
          wann%l_oc_f=.TRUE.
          ALLOCATE(wann%oc_orbs(wann%oc_num_orbs))
          DO n=1,wann%oc_num_orbs
             wann%oc_orbs(n)=n
          ENDDO
       ENDIF
    ENDIF

    IF(wann%l_updown)THEN            
       CALL wann_updown(&
            mpi,input,sym,atoms,stars,vacuum,sphhar,oneD,noco,cell,vTot,&
            enpara,eig_idList(1),l_real,&
            mpi%mpi_comm,atoms%l_dulo,noco%l_noco,noco%l_ss,&
            atoms%lmaxd,atoms%ntype,DIMENSION%neigd,atoms%nat,sym%nop,&
            DIMENSION%nvd,input%jspins,DIMENSION%nbasfcn,atoms%llod,&
            atoms%nlod,atoms%ntype,cell%omtil,atoms%nlo,atoms%llo,&
            atoms%lapw_l,sym%invtab,sym%mrot,atoms%ngopr,atoms%neq,&
            atoms%lmax,atoms%invsat,sym%invsatnr,nkpt,atoms%taual,&
            atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,&                    ! TODO: adapt if needed&
            stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,&
            stars%ng3,&
            vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,sphhar%nlhd,stars%ng3,&
            vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,&
            atoms%jri,sphhar%ntypsd,atoms%ntypsy,input%jspins,nkpt,&
            atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,&
            atoms%ulo_der,stars%ustep,stars%ig,stars%mx1,&
            stars%mx2,stars%mx3,stars%rgphs,&
            sliceplot%slice,sliceplot%kk,sliceplot%nnne,cell%z1,&
            DIMENSION%nv2d,vacuum%nmzxy,vacuum%nmz,vacuum%delz,&
            stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,results%ef,noco%l_soc,&
            sphhar%memd,atoms%lnonsph,sphhar%clnu,DIMENSION%lmplmd,&
            sphhar%mlh,sphhar%nmem,sphhar%llh,atoms%lo1l,&
            noco%theta,noco%phi)

       DO pc = 1, wann%nparampts
          CALL close_eig(eig_idList(pc))
       END DO

       CALL juDFT_end("updown done",mpi%irank)
    ENDIF

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  modern theory of orbital magnetization from Wannier functions
    !  Jan-Philipp Hanke
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IF(wann%l_matrixuHu)THEN
       wannTemp = wann
       CALL wann_uHu(&
            DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,&
            banddos,oneD,noco,cell,vTot,wannTemp,eig_idList,&
            l_real,atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,&
            atoms%ntype,DIMENSION%neigd,atoms%nat,sym%nop,DIMENSION%nvd,&
            input%jspins,DIMENSION%nbasfcn,atoms%llod,atoms%nlod,&
            atoms%ntype,cell%omtil,atoms%nlo,atoms%llo,&
            atoms%lapw_l,sym%invtab,sym%mrot,atoms%ngopr,atoms%neq,&
            atoms%lmax,atoms%invsat,sym%invsatnr,nkpt,atoms%taual,&
            atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,&
            mpi%irank,&
            mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,&
            sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,&
            input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,atoms%ntypsy,&
            input%jspins,nkpt,atoms%dx,stars%ng2,atoms%rmsh,&
            sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,stars%ustep,&
            stars%ig,stars%mx1,stars%mx2,stars%mx3,&
            stars%rgphs,sliceplot%slice,&
            sliceplot%kk,sliceplot%nnne,cell%z1,DIMENSION%nv2d,&
            vacuum%nmzxy,vacuum%nmz,vacuum%delz,sym%zrfs,stars%ig2,&
            cell%area,sym%tau,atoms%zatom,stars%ng2,stars%kv2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,results%ef,noco%l_soc,&
            sphhar%memd,atoms%lnonsph,sphhar%clnu,DIMENSION%lmplmd,&
            sphhar%mlh,sphhar%nmem,sphhar%llh,atoms%lo1l,&
            noco%theta,noco%phi,&
            wann%l_ms,wann%l_sgwf,wann%l_socgwf,wann%aux_latt_const,&
            wann%param_file,wann%param_vec,wann%nparampts,&
            wann%param_alpha,wann%l_dim)

       DO pc = 1, wann%nparampts
          CALL close_eig(eig_idList(pc))
       END DO

       CALL juDFT_end("wann_uHu done",mpi%irank)
    ENDIF

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !  modern theory of DMI from higher-dimensional Wannier functions
    !  Jan-Philipp Hanke
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IF(wann%l_matrixuHu_dmi)THEN
       wannTemp = wann
       CALL wann_uHu_dmi(&
            DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,&
            banddos,oneD,noco,cell,vTot,wannTemp,eig_idList,&
            l_real,atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,&
            atoms%ntype,DIMENSION%neigd,atoms%nat,sym%nop,DIMENSION%nvd,&
            input%jspins,DIMENSION%nbasfcn,atoms%llod,atoms%nlod,&
            atoms%ntype,cell%omtil,atoms%nlo,atoms%llo,&
            atoms%lapw_l,sym%invtab,sym%mrot,atoms%ngopr,atoms%neq,&
            atoms%lmax,atoms%invsat,sym%invsatnr,nkpt,atoms%taual,&
            atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,&
            mpi%irank,&
            mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,atoms%jmtd,&
            sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,&
            input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,atoms%ntypsy,&
            input%jspins,nkpt,atoms%dx,stars%ng2,atoms%rmsh,&
            sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,stars%ustep,&
            stars%ig,stars%mx1,stars%mx2,stars%mx3,&
            stars%rgphs,sliceplot%slice,&
            sliceplot%kk,sliceplot%nnne,cell%z1,DIMENSION%nv2d,&
            vacuum%nmzxy,vacuum%nmz,vacuum%delz,sym%zrfs,stars%ig2,&
            cell%area,sym%tau,atoms%zatom,stars%ng2,stars%kv2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,results%ef,noco%l_soc,&
            sphhar%memd,atoms%lnonsph,sphhar%clnu,DIMENSION%lmplmd,&
            sphhar%mlh,sphhar%nmem,sphhar%llh,atoms%lo1l,&
            noco%theta,noco%phi,&
            wann%l_ms,wann%l_sgwf,wann%l_socgwf,wann%aux_latt_const,&
            wann%param_file,wann%param_vec,wann%nparampts,&
            wann%param_alpha,wann%l_dim,l_nochi)

       DO pc = 1, wann%nparampts
          CALL close_eig(eig_idList(pc))
       END DO

       CALL juDFT_end("wann_uHu dmi done",mpi%irank)
    ENDIF

    IF(wann%l_byenergy.AND.wann%l_byindex) CALL juDFT_error&
         ("byenergy.and.byindex",calledby ="wannier")
    IF(wann%l_byenergy.AND.wann%l_bynumber) CALL juDFT_error&
         ("byenergy.and.bynumber",calledby ="wannier")
    IF(wann%l_bynumber.AND.wann%l_byindex) CALL juDFT_error&
         ("bynumber.and.byindex",calledby ="wannier")
    IF(.NOT.(wann%l_bynumber.OR.wann%l_byindex.OR.wann%l_byenergy))&
         CALL juDFT_error("no rule to sort bands",calledby ="wannier")


    efermi=results%ef
    IF(.NOT.wann%l_fermi)efermi=0.0

#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

    !**************************************************************
    !   for bzsym=.true.: determine mapping between kpts and w90kpts
    !**************************************************************
    IF (wann%l_bzsym) THEN
       l_file=.FALSE.
       INQUIRE(file='w90kpts',exist=l_file)
       IF(.NOT.l_file)  CALL juDFT_error&
            ("w90kpts not found, needed if bzsym",calledby ="wannier")
       OPEN(412,file='w90kpts',form='formatted')
       READ(412,*)fullnkpts
       CLOSE(412)
       IF(l_p0)PRINT*,"fullnkpts=",fullnkpts
       IF(fullnkpts<nkpts) CALL juDFT_error("fullnkpts.lt.nkpts"&
            ,calledby ="wannier")
       ALLOCATE(irreduc(fullnkpts),mapkoper(fullnkpts))
       ALLOCATE(shiftkpt(3,fullnkpts))
       l_file=.FALSE.
       INQUIRE(file='kptsmap',exist=l_file)
       IF(.NOT.l_file)  CALL juDFT_error&
            ("kptsmap not found, needed if bzsym",calledby ="wannier")
       OPEN(713,file='kptsmap')
       DO i=1,fullnkpts
          READ(713,*)kpt,irreduc(i),mapkoper(i),shiftkpt(:,i)
          IF(kpt/=i) CALL juDFT_error("kpt.ne.i",calledby ="wannier")
          IF(l_p0)PRINT*,i,irreduc(i),mapkoper(i)
       ENDDO
       CLOSE(713)
       IF(MAXVAL(irreduc(:))/=nkpts) CALL juDFT_error&
            ("max(irreduc(:))/=nkpts",calledby ="wannier")
    ELSE
       fullnkpts=nkpts
       ALLOCATE(irreduc(fullnkpts),mapkoper(fullnkpts))
       ALLOCATE(shiftkpt(3,fullnkpts))
    ENDIF


    IF(l_gwf) fullnqpts = wann%nparampts


    nrec = 0
    IF(l_p0)THEN
       WRITE (*,*) 'fermi energy:',efermi
       WRITE (*,*) 'emin,emax=',sliceplot%e1s,sliceplot%e2s
       WRITE (*,*) 'nbasfcn =',DIMENSION%nbasfcn
    ENDIF

    IF((.NOT.wann%l_matrixmmn).AND.(.NOT.wann%l_wann_plot).AND.&
         (.NOT.wann%l_matrixamn).AND.(.NOT.wann%l_projmethod).AND.&
         (.NOT.wann%l_bestproj).AND.(.NOT.wann%l_nabla).AND.&
         (.NOT.wann%l_mmn0).AND.(.NOT.wann%l_surfcurr).AND.&
         (.NOT.wann%l_offdiposop).AND.(.NOT.wann%l_anglmom).AND.&
         (.NOT.wann%l_orbcomp).AND.(.NOT.wann%l_perturb) .AND.&
         (.NOT.wann%l_finishgwf) ) GOTO 1911

    !**********************************************************
    !cccccccccccccc   read in the bkpts file  ccccccccccccccccc
    !**********************************************************
    IF (wann%l_matrixmmn) THEN ! for Omega functional minimization
       l_bkpts = .FALSE.
       INQUIRE (file='bkpts',exist=l_bkpts)
       IF (.NOT.l_bkpts)  CALL juDFT_error("need bkpts for matrixmmn"&
            ,calledby ="wannier")
       OPEN (202,file='bkpts',form='formatted',status='old')
       REWIND (202)
       READ (202,'(i4)') nntot
       IF(l_p0)THEN
          WRITE (*,*) 'nntot=',nntot
          WRITE(*,*) 'fullnkpts=',fullnkpts
          WRITE(*,*) 'nkpts=',nkpts
       ENDIF
       ALLOCATE ( gb(1:3,1:nntot,1:fullnkpts),bpt(1:nntot,1:fullnkpts))
       DO ikpt=1,fullnkpts
          DO nn=1,nntot
             READ (202,'(2i6,3x,3i4)')&
                  ikpt_help,bpt(nn,ikpt),(gb(i,nn,ikpt),i=1,3)
             IF (ikpt/=ikpt_help)  CALL juDFT_error("ikpt.ne.ikpt_help"&
                  ,calledby ="wannier")
             IF (bpt(nn,ikpt)>fullnkpts) CALL juDFT_error("bpt.gt.fullnkpts"&
                  ,calledby ="wannier")
          ENDDO
       ENDDO
       CLOSE (202)
       ALLOCATE(kdiff(3,nntot))
    ENDIF

    !**********************************************************
    !cccccccccccccc   read in the bqpts file  ccccccccccccccccc         
    !**********************************************************
    IF ((wann%l_matrixmmn).AND.(l_gwf.OR.wann%l_ms)) THEN
       l_bqpts = .FALSE.
       INQUIRE (file='bqpts',exist=l_bqpts)
       IF (.NOT.l_bqpts)  CALL juDFT_error("need bqpts for matrixmmn"&
            ,calledby ="wannier")
       OPEN (202,file='bqpts',form='formatted',status='old')
       REWIND (202)
       READ (202,'(i4)') nntot_q
       IF(l_p0)THEN
          WRITE (*,*) 'nntot_q=',nntot_q
          WRITE(*,*) 'fullnqpts=',fullnqpts
       ENDIF
       ALLOCATE ( gb_q(1:3,1:nntot_q,1:fullnqpts),&
            bpt_q(1:nntot_q,1:fullnqpts))
       DO iqpt=1,fullnqpts
          DO nn=1,nntot_q
             READ (202,'(2i6,3x,3i4)')&
                  iqpt_help,bpt_q(nn,iqpt),(gb_q(i,nn,iqpt),i=1,3)
             IF (iqpt/=iqpt_help)  CALL juDFT_error("iqpt.ne.iqpt_help"&
                  ,calledby ="wannier")
             IF (bpt_q(nn,iqpt)>fullnqpts)&
                  CALL juDFT_error("bpt_q.gt.fullnqpts",calledby ="wannier")
          ENDDO
       ENDDO
       CLOSE (202)
       ALLOCATE(qdiff(3,nntot_q))
       ALLOCATE(zero_qdiff(3,nntot_q))
       zero_qdiff=0.0
    ENDIF


    ! when treating gen. WF for spin spirals, the Brillouin zone
    ! of q-points is twice as large compared to k-BZ. Thus,
    ! the G-vectors connecting neighbors across the boundary
    ! need to be doubled
    IF(wann%l_sgwf) gb_q = 2*gb_q    
    IF(wann%l_socgwf) gb_q = 2*gb_q 

    IF(wann%l_finishgwf) GOTO 9110
    !********************************************************
    !      find symmetry-related elements in mmkb
    !********************************************************
    IF(wann%l_matrixmmn)THEN
       ALLOCATE(maptopair(3,fullnkpts,nntot))
       ALLOCATE(pair_to_do(fullnkpts,nntot))
       CALL wann_mmnk_symm(input,kpts,&
            fullnkpts,nntot,bpt,gb,wann%l_bzsym,&
            irreduc,mapkoper,l_p0,input%film,sym%nop,sym%invtab,sym%mrot,&
            oneD%odi%d1,sym%tau,&
            pair_to_do,maptopair,kdiff,.FALSE.,wann%param_file)
    ENDIF

    ! do the same for q-points to construct GWFs
    IF(wann%l_matrixmmn.AND.l_gwf)THEN 
       ALLOCATE(maptopair_q(3,fullnqpts,nntot_q))
       ALLOCATE(pair_to_do_q(fullnqpts,nntot_q))
       CALL wann_mmnk_symm(input,kpts,&
            fullnqpts,nntot_q,bpt_q,gb_q,wann%l_bzsym,&
            irreduc_q,mapqoper,l_p0,.FALSE.,1,sym%invtab(1),&
            sym%mrot(:,:,1),.FALSE.,sym%tau,&
            pair_to_do_q,maptopair_q,qdiff,.TRUE.,wann%param_file)
    ENDIF


    !*********************************************************
    !ccccccccccccccc   initialize the potential   cccccccccccc
    !*********************************************************

    ALLOCATE ( vz(vacuum%nmzd,2,4) )
    ALLOCATE ( vr(atoms%jmtd,atoms%ntype,input%jspins) )
    ALLOCATE ( vso(atoms%jmtd,atoms%nat,2) )

    vz = 0.0
    vz(:,:,:SIZE(vTot%vacz,3)) = vTot%vacz(:,:,:)

    DO jspin = 1,input%jspins
       DO n = 1, atoms%ntype
          DO j = 1,atoms%jri(n)
             vr(j,n,jspin) = vTot%mt(j,0,n,jspin)
          ENDDO
       ENDDO
    ENDDO

    IF(wann%l_soctomom)THEN
       CALL vsoc(input,atoms,vr,enpara%el0,.TRUE., vso)
    ENDIF

    IF(noco%l_noco.AND.input%film)THEN
       npotmatfile=25
       ALLOCATE(vpw(stars%ng3,1))
       IF(.NOT.oneD%odi%d1)&
            ALLOCATE( vxy(vacuum%nmzxyd,stars%ng2-1,2) )

       OPEN (npotmatfile,FILE='potmat',FORM='unformatted',&
            STATUS='old')
       READ (npotmatfile) (vpw(ig3,1),ig3=1,stars%ng3)
       READ (npotmatfile) (vpw(ig3,1),ig3=1,stars%ng3)
       READ (npotmatfile) (vpw(ig3,1),ig3=1,stars%ng3)
       maxvac=2
       IF(oneD%odi%d1)maxvac=1
       DO ivac = 1,maxvac
          !--->       if the two vacuua are equivalent, the potential file has to
          !--->       be backspaced, because the potential is the same at both
          !--->       surfaces of the film
          IF ((ivac.EQ.2) .AND. (vacuum%nvac.EQ.1)) THEN
             DO irec = 1,4
                BACKSPACE (npotmatfile)
             ENDDO
          ENDIF
          !--->       load the non-warping part of the potential
          READ (npotmatfile)&
               ((vz(imz,ivac,ipot),imz=1,vacuum%nmzd),ipot=1,4)

          IF(.NOT.oneD%odi%d1)THEN
             DO ipot = 1,3
                READ (npotmatfile)((vxy(imz,igvm2,ivac),&
                     imz=1,vacuum%nmzxy),igvm2=1,stars%ng2-1)
             ENDDO
          ENDIF
       ENDDO
       CLOSE (npotmatfile)
    ENDIF

    !ccccccccccccccc   end of the potential part  ccccccccccc
    wannierspin=input%jspins
    IF(noco%l_soc) wannierspin=2


    ALLOCATE ( ff(atoms%ntype,atoms%jmtd,2,0:atoms%lmaxd,2) )
    ALLOCATE ( gg(atoms%ntype,atoms%jmtd,2,0:atoms%lmaxd,2) )
    ALLOCATE ( usdus%us(0:atoms%lmaxd,atoms%ntype,2) )
    ALLOCATE ( usdus%uds(0:atoms%lmaxd,atoms%ntype,2) )
    ALLOCATE ( usdus%dus(0:atoms%lmaxd,atoms%ntype,2) )
    ALLOCATE ( usdus%duds(0:atoms%lmaxd,atoms%ntype,2) )
    ALLOCATE ( usdus%ddn(0:atoms%lmaxd,atoms%ntype,2) )
    ALLOCATE ( usdus%ulos(atoms%nlod,atoms%ntype,2) )
    ALLOCATE ( usdus%dulos(atoms%nlod,atoms%ntype,2) )
    ALLOCATE ( usdus%uulon(atoms%nlod,atoms%ntype,2) )
    ALLOCATE ( usdus%dulon(atoms%nlod,atoms%ntype,2) )
    ALLOCATE ( usdus%uloulopn(atoms%nlod,atoms%nlod,atoms%ntype,2) )

    IF(l_gwf.AND..NOT.(wann%l_wann_plot)) THEN
       doublespin_max=4!2
    ELSE
       doublespin_max=wannierspin
    ENDIF

    !*****************************************************************c
    !                         START Q LOOP                            c
    !   standard functionality of code for fullnqpts = nntot_q = 1    c
    !        and wann%l_ms = wann%l_sgwf = wann%l_socgwf = F          c
    !*****************************************************************c
    DO iqpt = 1,fullnqpts  ! loop by q-points starts

       ALLOCATE(innerEig_idList(nntot_q))

       qptibz=iqpt                          
       IF(wann%l_bzsym .AND. l_gwf) qptibz=irreduc_q(iqpt)
       IF(wann%l_bzsym .AND. l_gwf) oper_q=mapqoper(iqpt)

       qpt_i = noco%qss
       alph_i = noco%alph
       beta_i = noco%beta
       theta_i = noco%theta
       phi_i = noco%phi
       IF(wann%l_sgwf.OR.wann%l_ms) THEN
          qpt_i(:) = wann%param_vec(:,qptibz)
          alph_i(:) = wann%param_alpha(:,qptibz)
       ELSEIF(wann%l_socgwf) THEN 
          IF(wann%l_dim(2)) phi_i = tpi_const*wann%param_vec(2,qptibz)
          IF(wann%l_dim(3)) theta_i = tpi_const*wann%param_vec(3,qptibz)
       ENDIF

       IF (l_gwf) THEN
          IF(wann%l_matrixmmn)THEN
             DO iqpt_b=1,nntot_q

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

             ENDDO
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

       ENDIF ! l_gwf.or.wann%l_ms
       nrec=0
       nrec_b=0


       !****************************************************
       ! cycle by spins starts! 
       !****************************************************
       DO doublespin=1,doublespin_max   ! cycle by spins

          jspin=MOD(doublespin+1,2)+1
          jspin_b=jspin
          IF(doublespin.EQ.3) jspin_b=2
          IF(doublespin.EQ.4) jspin_b=1

          nrec_b = nrec

          IF(.NOT.noco%l_noco) THEN
             nrec = (jspin-1)*nkpts
             nrec_b = (jspin_b-1)*nkpts
          ENDIF

          ! spin-dependent sign of the q-dependent phase
          ! in the generalized Bloch theorem
          ! -1: spin up, +1: spin down
          sign_q = -sign_q

          !...read number of bands and wannier functions from file proj

          !..reading the proj.1 / proj.2 / proj file
          l_proj=.FALSE.  
          DO j=jspin,0,-1
             INQUIRE(file=TRIM('proj'//spin012(j)),exist=l_proj)
             IF(l_proj)THEN
                filename='proj'//spin012(j)
                EXIT
             ENDIF
          ENDDO

          IF(l_proj)THEN
             OPEN (203,file=TRIM(filename),status='old')
             REWIND (203)
             READ (203,*) nwfs,numbands
             REWIND (203)
             CLOSE (203)
          ELSEIF(wann%l_projmethod.OR.wann%l_bestproj&
               .OR.wann%l_matrixamn)THEN
             CALL juDFT_error("no proj/proj.1/proj.2"&
                  ,calledby ="wannier")
          ENDIF


          jspin2=jspin
          IF(noco%l_soc .AND. input%jspins.EQ.1)jspin2=1
          jspin2_b=jspin_b
          IF(noco%l_soc .AND. input%jspins.EQ.1)jspin2_b=1

          jsp_start = jspin ; jsp_end = jspin

          !ccccccccccc   read in the eigenvalues and vectors   cccccc
          WRITE(*,*)'wannierspin',wannierspin
          DO jspin5=1,wannierspin!1!2
             !       jspin5=jspin
             jsp_start=jspin5; jsp_end=jspin5
             nrec5=0
             IF(.NOT.noco%l_noco) nrec5 = (jspin5-1)*nkpts

             CALL cdn_read0(eig_id,mpi%irank,mpi%isize,jspin5,input%jspins, &!wannierspin instead of DIMENSION%jspd?&
                  noco%l_noco, n_bands,n_size)

          ENDDO
          !..   now we want to define the maximum number of the bands by all kpts
          nbnd = 0

          i_rec = 0 ; n_rank = 0
          !*************************************************************
          !..writing down the eig.1 and/or eig.2 files

          !..write individual files if multi-spiral mode wann%l_ms=T
          !*************************************************************
          IF(l_p0)THEN         
             CALL wann_write_eig(&
                  eig_id,l_real,&
                  atoms%lmaxd,atoms%ntype,atoms%nlod,DIMENSION%neigd,&
                  DIMENSION%nvd,wannierspin,&
                  mpi%isize,jspin,DIMENSION%nbasfcn,atoms%nlotot,&
                  noco%l_ss,noco%l_noco,nrec,fullnkpts,&
                  wann%l_bzsym,wann%l_byindex,wann%l_bynumber,&
                  wann%l_byenergy,&
                  irreduc,oneD%odi,wann%band_min(jspin),&
                  wann%band_max(jspin),&
                  numbands,&
                  sliceplot%e1s,sliceplot%e2s,efermi,.FALSE.,nkpts,&
                  nbnd,kpoints,l_gwf,iqpt)       

             IF(oneD%odi%d1)THEN
                kpoints(:)=kpoints(:)*cell%bmat(3,3)         
             ENDIF
          ENDIF!l_p0

          ! nbnd is calculated for process zero and is sent here to the others
#ifdef CPP_MPI
          IF(l_p0)THEN
             DO cpu_index=1,mpi%isize-1
                CALL MPI_SEND(nbnd,1,MPI_INTEGER,cpu_index,1,mpi%mpi_comm,ierr)
             ENDDO
          ELSE
             CALL MPI_RECV(nbnd,1,MPI_INTEGER,0,1,mpi%mpi_comm,stt,ierr)
          ENDIF
#endif

          PRINT*,"process: ",mpi%irank," nbnd= ",nbnd
          !##################################################################
          IF(wann%l_mmn0)THEN
             ALLOCATE ( mmn(nbnd,nbnd,fullnkpts) )
             mmn(:,:,:) = CMPLX(0.,0.)
             IF((noco%l_soc.OR.noco%l_noco) .AND. (doublespin.EQ.1))&
                  ALLOCATE(socmmn(nbnd,nbnd,fullnkpts) )
          ENDIF
          IF(wann%l_nabla)THEN
             ALLOCATE ( nablamat(3,nbnd,nbnd,fullnkpts) )
             nablamat = CMPLX(0.,0.)
          ENDIF

          IF(wann%l_soctomom)THEN
             ALLOCATE ( soctomom(3,nbnd,nbnd,fullnkpts) )
             soctomom = CMPLX(0.,0.)
          ENDIF

          IF(wann%l_surfcurr)THEN
             ALLOCATE ( surfcurr(3,nbnd,nbnd,fullnkpts) )
             surfcurr = CMPLX(0.,0.)
          ENDIF

          IF(wann%l_anglmom)THEN
             IF(.NOT.ALLOCATED(anglmom))THEN  
                ALLOCATE ( anglmom(3,nbnd,nbnd,fullnkpts) )
                anglmom=CMPLX(0.,0.)
             ENDIF
          ENDIF

          IF(wann%l_orbcomp)THEN
             IF(ALLOCATED(orbcomp))DEALLOCATE(orbcomp)
             IF(wann%l_oc_f)THEN
                ALLOCATE(orbcomp(16,wann%oc_num_orbs,nbnd,nbnd,fullnkpts))
             ELSE
                ALLOCATE(orbcomp(9,wann%oc_num_orbs,nbnd,nbnd,fullnkpts))
             ENDIF
             orbcomp=CMPLX(0.,0.)
          ENDIF

          !write (*,*) 'nwfs=',nwfs
          IF(wann%l_projmethod.OR.wann%l_bestproj.OR.wann%l_matrixamn)THEN
             IF(.NOT.ALLOCATED(amn))THEN
                ALLOCATE ( amn(nbnd,nwfs,fullnkpts) )
                amn(:,:,:) = CMPLX(0.,0.)
             ENDIF
          ENDIF

          IF (wann%l_projmethod.OR.wann%l_bestproj) THEN
             ALLOCATE ( psiw(nbnd,nwfs,fullnkpts) )
             psiw(:,:,:) = CMPLX(0.,0.)
             IF(.NOT.ALLOCATED(hwfr))THEN
                ALLOCATE ( hwfr(nwfs,nwfs) )
                hwfr(:,:) = CMPLX(0.,0.)
             ENDIF
          ENDIF


          IF (wann%l_matrixmmn) THEN
             IF(.NOT.ALLOCATED(mmnk))THEN
                ALLOCATE ( mmnk(nbnd,nbnd,nntot,fullnkpts) )
                mmnk = (0.,0.)
             ENDIF
          ENDIF

          IF(wann%l_matrixmmn)THEN
             IF(.NOT.ALLOCATED(mmnk_q).AND.l_gwf)THEN
                ALLOCATE ( mmnk_q (nbnd,nbnd,nntot_q,fullnkpts) )
                mmnk_q = (0.,0.)

                !             allocate ( m_int(nbnd,nbnd,nntot_q,fullnkpts) )
                !             allocate ( m_sph(nbnd,nbnd,nntot_q,fullnkpts) )
                !             allocate ( m_vac(nbnd,nbnd,nntot_q,fullnkpts) )
                !             m_int = cmplx(0.,0.)
                !             m_sph = cmplx(0.,0.)
                !             m_vac = cmplx(0.,0.)
             ENDIF
          ENDIF


          ALLOCATE ( flo(atoms%ntype,atoms%jmtd,2,atoms%nlod,2) )

          DO jspin4=1,wannierspin!2
             jspin3=jspin4
             IF(input%jspins.EQ.1) jspin3=1
             na = 1
             DO  n = 1,atoms%ntype
                DO  l = 0,atoms%lmax(n)
                   !...compute the l-dependent, k-independent radial MT- basis functions

                   CALL radfun(&
                        l,n,jspin4,enpara%el0(l,n,jspin3),vr(1,n,jspin3),atoms,&
                        ff(n,:,:,l,jspin4),gg(n,:,:,l,jspin4),usdus,&
                        nodeu,noded,wronk)

                ENDDO
                !...and the local orbital radial functions
                DO ilo = 1, atoms%nlo(n)

                   CALL radflo(&
                        atoms,n,jspin4,enpara%ello0(:,:,jspin3),vr(1,n,jspin3),&
                        ff(n,1:,1:,0:,jspin4),gg(n,1:,1:,0:,jspin4),mpi,&
                        usdus,uuilon,duilon,ulouilopn,flo(n,:,:,:,jspin4))

                ENDDO
                !       na = na + atoms%neq(n)
             ENDDO
          ENDDO!jspin3
          !****************************************************************
          !   calculate the k-independent uju*gaunt-matrix needed for
          !   mmnmatrix
          !****************************************************************
          ! TODO: make this more efficient (i.e., compute ujugaunt only once
          ! and not for all q-points).
          IF(wann%l_matrixmmn)THEN
             ALLOCATE(ujug(0:lmd,0:lmd,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(ujdg(0:lmd,0:lmd,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(djug(0:lmd,0:lmd,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(djdg(0:lmd,0:lmd,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(ujulog(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(djulog(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(ulojug(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(ulojdg(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                  1:atoms%ntype,1:nntot))
             ALLOCATE(ulojulog(1:atoms%nlod,-atoms%llod:atoms%llod,&
                  1:atoms%nlod,-atoms%llod:atoms%llod,&
                  1:atoms%ntype,1:nntot))

             CALL wann_ujugaunt(&
                  atoms%llod,nntot,kdiff,atoms%lmax,atoms%ntype,&
                  atoms%ntype,cell%bbmat,cell%bmat,atoms%nlod,atoms%nlo,&
                  atoms%llo,flo(:,:,:,:,jspin),&
                  flo(:,:,:,:,jspin),&
                  ff(:,:,:,:,jspin),&
                  ff(:,:,:,:,jspin),&
                  gg(:,:,:,:,jspin),&
                  gg(:,:,:,:,jspin),atoms%jri,atoms%rmsh,atoms%dx,&
                  atoms%jmtd,atoms%lmaxd,lmd,&
                  ujug,ujdg,djug,djdg,&
                  ujulog,djulog,ulojug,ulojdg,ulojulog,.FALSE.,1)

             ! compute integrals of radial solution, according energy derivatives,
             ! the spherical Bessel function and the Gaunt coefficients in order
             ! to account for the overlap of the lattice periodic parts at
             ! neighboring q-points
             IF(l_gwf)THEN
                ALLOCATE(ujug_q(0:lmd,0:lmd,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(ujdg_q(0:lmd,0:lmd,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(djug_q(0:lmd,0:lmd,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(djdg_q(0:lmd,0:lmd,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(ujulog_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(djulog_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(ulojug_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(ulojdg_q(0:lmd,1:atoms%nlod,-atoms%llod:atoms%llod,&
                     1:atoms%ntype,1:nntot_q))
                ALLOCATE(ulojulog_q(1:atoms%nlod,-atoms%llod:atoms%llod,&
                     1:atoms%nlod,-atoms%llod:atoms%llod,&
                     1:atoms%ntype,1:nntot_q))

                ! we need G(q+b)/2 as argument for the sph. Bessel func.
                ! and additionally a spin-dependent sign (-/+ 1)^{lpp}
                IF(wann%l_sgwf) CALL wann_ujugaunt(&
                     atoms%llod,nntot_q,qdiff/2.0,atoms%lmax,atoms%ntype,&
                     atoms%ntype,cell%bbmat,cell%bmat,atoms%nlod,atoms%nlo,&
                     atoms%llo,flo(:,:,:,:,jspin),&
                     flo(:,:,:,:,jspin_b),&
                     ff(:,:,:,:,jspin),&
                     ff(:,:,:,:,jspin_b),&
                     gg(:,:,:,:,jspin),&
                     gg(:,:,:,:,jspin_b),atoms%jri,atoms%rmsh,atoms%dx,&
                     atoms%jmtd,atoms%lmaxd,lmd,&
                     ujug_q,ujdg_q,djug_q,djdg_q,&
                     ujulog_q,djulog_q,ulojug_q,ulojdg_q,ulojulog_q,.TRUE.,&
                     sign_q)

                IF(wann%l_socgwf) CALL wann_ujugaunt(&
                     atoms%llod,nntot_q,zero_qdiff,atoms%lmax,atoms%ntype,&
                     atoms%ntype,cell%bbmat,cell%bmat,atoms%nlod,atoms%nlo,&
                     atoms%llo,flo(:,:,:,:,jspin),&
                     flo(:,:,:,:,jspin_b),&
                     ff(:,:,:,:,jspin),&
                     ff(:,:,:,:,jspin_b),&
                     gg(:,:,:,:,jspin),&
                     gg(:,:,:,:,jspin_b),atoms%jri,atoms%rmsh,atoms%dx,&
                     atoms%jmtd,&
                     atoms%lmaxd,lmd,ujug_q,ujdg_q,djug_q,djdg_q,&
                     ujulog_q,djulog_q,ulojug_q,ulojdg_q,ulojulog_q,&
                     .FALSE.,1)

             ENDIF ! l_gwf

          ENDIF !l_matrixmmn
          zzMat%l_real = l_real
          zzMat%matsize1 = DIMENSION%nbasfcn
          zzMat%matsize2 = DIMENSION%neigd
          IF(l_real) THEN
             IF(.NOT.ALLOCATED(zzMat%data_r))&
               ALLOCATE (zzMat%data_r(zzMat%matsize1,zzMat%matsize2))
          ELSE
               IF(.NOT.ALLOCATED(zzMat%data_c))&
               ALLOCATE (zzMat%data_c(zzMat%matsize1,zzMat%matsize2))
          END IF

          zMat%l_real = zzMat%l_real
          zMat%matsize1 = zzMat%matsize1
          zMat%matsize2 = zzMat%matsize2
          IF (zzMat%l_real) THEN
             IF(.NOT.ALLOCATED(zMat%data_r))&
                  ALLOCATE (zMat%data_r(zMat%matsize1,zMat%matsize2))
             zMat%data_r = 0.0
          ELSE
             IF(.NOT.ALLOCATED(zMat%data_c))&
                  ALLOCATE (zMat%data_c(zMat%matsize1,zMat%matsize2))
             zMat%data_c = CMPLX(0.0,0.0)
          END IF

          zMat_b%l_real = zzMat%l_real
          zMat_b%matsize1 = zzMat%matsize1
          zMat_b%matsize2 = zzMat%matsize2
          IF (zzMat%l_real) THEN
               IF(.NOT.ALLOCATED(zMat_b%data_r))&
               ALLOCATE (zMat_b%data_r(zMat_b%matsize1,zMat_b%matsize2))
             zMat_b%data_r = 0.0
          ELSE
               IF(.NOT.ALLOCATED(zMat_b%data_c))&
               ALLOCATE (zMat_b%data_c(zMat_b%matsize1,zMat_b%matsize2))
             zMat_b%data_c = CMPLX(0.0,0.0)
          END IF

          i_rec = 0 ; n_rank = 0

          !****************************************************************
          !.. loop by kpoints starts!      each may be a separate task
          !****************************************************************
          DO ikpt = wann%ikptstart,fullnkpts  ! loop by k-points starts
             kptibz=ikpt
             IF(wann%l_bzsym) kptibz=irreduc(ikpt)
             IF(wann%l_bzsym) oper=mapkoper(ikpt)

             i_rec = i_rec + 1
             IF (MOD(i_rec-1,mpi%isize).EQ.mpi%irank) THEN

                ALLOCATE ( we(DIMENSION%neigd),eigg(DIMENSION%neigd) )

                n_start=1
                n_end=DIMENSION%neigd


                ! read information of diagonalization for fixed q-point iqpt
                ! stored in the eig file on unit 66. the lattice respectively
                ! plane-wave vectors G(k,q) are saved in (k1,k2,k3).

                CALL lapw%init(input,noco,kpts,atoms,sym,kptibz,cell,(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco.and.mpi%n_size==1),mpi)

                CALL cdn_read(&
                     eig_id,&
                     DIMENSION%nvd,input%jspins,mpi%irank,mpi%isize, &!wannierspin instead of DIMENSION%jspd?&
                     kptibz,jspin,DIMENSION%nbasfcn,&
                     noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,n_end,&
                     nbands,eigg,zzMat)


                nslibd = 0

                !...we work only within the energy window

                eig(:) = 0.

                !      print*,"bands used:"

                DO i = 1,nbands
                   IF ((eigg(i).GE.sliceplot%e1s.AND.nslibd.LT.numbands.AND.&
                        wann%l_bynumber).OR.&
                        (eigg(i).GE.sliceplot%e1s.AND.eigg(i).LE.sliceplot%e2s.AND.&
                        wann%l_byenergy).OR.(i.GE.wann%band_min(jspin).AND.&
                        (i.LE.wann%band_max(jspin)).AND.wann%l_byindex))THEN

                      !           print*,i
                      nslibd = nslibd + 1
                      eig(nslibd) = eigg(i)
                      we(nslibd) = we(i)
                      IF(zzMat%l_real) THEN
                         zMat%data_r(:,nslibd) = zzMat%data_r(:,i)
                      ELSE
                         zMat%data_c(:,nslibd) = zzMat%data_c(:,i)
                      END IF
                   ENDIF
                ENDDO

                !***********************************************************
                !              rotate the wavefunction
                !***********************************************************
                IF (wann%l_bzsym.AND.oper.NE.1) THEN  !rotate bkpt
                   !         call wann_kptsrotate(
                   !     >            atoms%nat,atoms%nlod,atoms%llod,
                   !     >            atoms%ntype,atoms%nlo,atoms%llo,atoms%invsat,
                   !     >            noco%l_noco,noco%l_soc,
                   !     >            atoms%ntype,atoms%neq,atoms%nlotot,
                   !     >            kveclo,jspin,
                   !     >            oper,sym%nop,sym%mrot,DIMENSION%nvd,nv,
                   !     >            shiftkpt(:,ikpt),
                   !     >            sym%tau,
                   !     x            lapw%bkpt,k1(:,:),k2(:,:),k3(:,:),
                   !     x            zMat,nsfactor)
                ELSE
                   nsfactor=CMPLX(1.0,0.0)
                ENDIF

                !******************************************************************

                !...the overlap matrix Mmn which is computed for each k- and b-point

                noccbd = nslibd

                ALLOCATE(acof(noccbd,0:lmd,atoms%nat),&
                     bcof(noccbd,0:lmd,atoms%nat),&
                     ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))

                acof(:,:,:) = CMPLX(0.,0.) ; bcof(:,:,:) = CMPLX(0.,0.)
                ccof(:,:,:,:) = CMPLX(0.,0.)

                !...generation the A,B,C coefficients in the spheres
                !...for the lapws and local orbitals, summed by the basis functions


                CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,&
                     noco,jspin,oneD,acof,bcof,ccof,zMat)


                CALL wann_abinv(atoms, acof,bcof,ccof)


                IF((doublespin.EQ.3).OR.(doublespin.EQ.4)) GOTO 9900


                IF(wann%l_orbcomp)THEN
                   CALL wann_orbcomp(atoms,usdus,jspin,acof,bcof,ccof,&
                        wann%oc_num_orbs,wann%oc_orbs,wann%l_oc_f,orbcomp(:,:,:,:,ikpt))
                ENDIF

                IF(wann%l_anglmom)THEN
                   CALL wann_anglmom(atoms,usdus,jspin, acof,bcof,ccof,anglmom(:,:,:,ikpt))
                ENDIF

#ifdef CPP_TOPO
                IF(wann%l_surfcurr)THEN
                   !         call wann_surfcurr_int(
                   !     >        DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,vacuum%nmzxyd,stars%ng2,sphhar%ntypsd,
                   !     >        atoms%ntype,atoms%lmaxd,atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,atoms%neq,stars%ng3,vacuum%nvac,
                   !     >        vacuum%nmz,vacuum%nmzxy,stars%ng2,sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,sym%symor,
                   !     >        sym%invs,sym%invs2,cell%z1,vacuum%delz,atoms%ngopr,atoms%ntypsy,atoms%jri,atoms%pos,atoms%zatom,
                   !     >        atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,cell%amat,cell%bmat,cell%bbmat,ikpt,sliceplot%nnne,sliceplot%kk,
                   !     >        DIMENSION%nvd,atoms%nlod,atoms%llod,nv(jspin),lmd,lapw%bkpt,cell%omtil,atoms%nlo,atoms%llo,
                   !     >        k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
                   !     >        vz(:,:,jspin2),
                   !     >        nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,ff,gg,flo,acof,bcof,ccof,z,
                   !     >        surfcurr(:,:,:,ikpt))

                   CALL wann_surfcurr_int2(&
                        DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,&
                        vacuum%nmzxyd,&
                        stars%ng2,sphhar%ntypsd,atoms%ntype,atoms%lmaxd,&
                        atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,&
                        atoms%neq,stars%ng3,vacuum%nvac,vacuum%nmz,&
                        vacuum%nmzxy,stars%ng2,sym%nop,sym%nop2,cell%volint,&
                        input%film,sliceplot%slice,sym%symor,&
                        sym%invs,sym%invs2,cell%z1,vacuum%delz,atoms%ngopr,&
                        atoms%ntypsy,atoms%jri,atoms%pos,atoms%taual,&
                        atoms%zatom,atoms%rmt,atoms%lmax,sym%mrot,sym%tau,&
                        atoms%rmsh,sym%invtab,cell%amat,cell%bmat,cell%bbmat,&
                        ikpt,DIMENSION%nvd,lapw%nv(jspin),lapw%bkpt,cell%omtil,&
                        lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),&
                        nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,z,&
                        dirfacs,&
                        surfcurr(:,:,:,ikpt))

                   CALL wann_surfcurr(&
                        dirfacs,cell%amat,&
                        jspin,atoms%ntype,atoms%lmaxd,atoms%lmax,atoms%nat,&
                        atoms%neq,noccbd,lmd,atoms%nat,atoms%llod,atoms%nlod,&
                        atoms%nlo,atoms%llo, &
                        acof,bcof,ccof,&
                        us(:,:,jspin),dus(:,:,jspin),duds(:,:,jspin),&
                        uds(:,:,jspin),&
                        ulos(:,:,jspin),dulos(:,:,jspin),&
                        atoms%rmt,atoms%pos, &
                        surfcurr(:,:,:,ikpt))
                   WRITE(6,*)"dirfacs=",dirfacs
                ENDIF

                IF(wann%l_soctomom)THEN
                   CALL wann_soc_to_mom(&
                        jspin,atoms%ntype,atoms%lmaxd,atoms%lmax,atoms%nat,&
                        atoms%jmtd,atoms%jri,atoms%rmsh,atoms%dx,atoms%neq,&
                        noccbd,lmd,atoms%nat,atoms%llod,atoms%nlod,&
                        vso(:,:,1), &
                        ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),&
                        acof,bcof,ccof,&
                        soctomom(:,:,:,ikpt))
                ENDIF

                IF(wann%l_nabla)THEN
                   CALL wann_nabla(&
                        atoms%nlo,atoms%llo,&
                        jspin,atoms%ntype,atoms%lmaxd,atoms%lmax,atoms%nat,&
                        atoms%jmtd,atoms%jri,atoms%rmsh,atoms%dx,atoms%neq,&
                        noccbd,lmd,atoms%nat,atoms%llod,atoms%nlod, &
                        ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),flo(:,:,:,:,jspin),&
                        acof,bcof,ccof,&
                        nablamat(:,:,:,ikpt))
                   IF(input%film.AND..NOT.oneD%odi%d1)THEN
                      CALL wann_nabla_vac(&
                           cell%z1,vacuum%nmzd,DIMENSION%nv2d,&
                           stars%mx1,stars%mx2,stars%mx3,&
                           stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,vacuum%delz,&
                           stars%ig2,cell%area,cell%bmat,cell%bbmat,enpara%evac0(:,jspin),&
                           lapw%bkpt,vz(:,:,jspin2),nslibd,jspin,lapw%k1,lapw%k2,lapw%k3,&
                           wannierspin,DIMENSION%nvd,DIMENSION%nbasfcn,DIMENSION%neigd,z,nv,&
                           cell%omtil,&
                           nablamat(:,:,:,ikpt))
                   ENDIF
                   addnoco=0
                   DO  i = n_rank+1,lapw%nv(jspin),n_size
                      b1 = lapw%bkpt+lapw%gvec(:,i,jspin)
                      b2(1)=b1(1)*cell%bmat(1,1)+b1(2)*cell%bmat(2,1)+&
                           b1(3)*cell%bmat(3,1)
                      b2(2)=b1(1)*cell%bmat(1,2)+b1(2)*cell%bmat(2,2)+&
                           b1(3)*cell%bmat(3,2)
                      b2(3)=b1(1)*cell%bmat(1,3)+b1(2)*cell%bmat(2,3)+&
                           b1(3)*cell%bmat(3,3)
                      DO  j = n_rank+1,lapw%nv(jspin),n_size
                         !-->     determine index and phase factor
                         i1 = lapw%k1(j,jspin) - lapw%k1(i,jspin)
                         i2 = lapw%k2(j,jspin) - lapw%k2(i,jspin)
                         i3 = lapw%k3(j,jspin) - lapw%k3(i,jspin)
                         in = stars%ig(i1,i2,i3)
                         IF (in.EQ.0) CYCLE
                         phase   = stars%rgphs(i1,i2,i3)
                         phasust = CMPLX(phase,0.0)*stars%ustep(in)

                         DO m = 1,nslibd
                            DO n = 1,nslibd
                               DO dir=1,3  

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
                                  VALUE=phasust*z(i+addnoco,m)*CONJG(z(j+addnoco,n))
                                  nablamat(dir,m,n,ikpt) =&
                                       nablamat(dir,m,n,ikpt) -&
                                       VALUE*b2(dir)
#else
                                  VALUE=phasust*CMPLX(z(i+addnoco,m)*z(j+addnoco,n),0.0)
                                  nablamat(dir,m,n,ikpt) =&
                                       nablamat(dir,m,n,ikpt) - &
                                       VALUE*b2(dir)
#endif
                               ENDDO
                            ENDDO
                         ENDDO

                      ENDDO
                   ENDDO

                ENDIF
#endif
                !      goto jump no longer needed?
                !      if ((.not.wann%l_matrixmmn).and.(.not.wann%l_matrixamn).and.
                !          (.not.wann%l_bestproj).and.(.not.wann%l_projmethod).and.
                !          (.not.wann%l_mmn0)) goto 3


                !------mmn0-matrix
                IF(wann%l_mmn0)THEN
                   addnoco=0
                   IF(noco%l_noco.AND.(jspin.EQ.2))THEN
                      addnoco=lapw%nv(1)+atoms%nlotot
                   ENDIF

                   !-----> interstitial contribution to mmn0-matrix

                   CALL wann_mmkb_int(&
                        cmplx_1,addnoco,addnoco,&
                        DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,&
                        stars%ng3,lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),&
                        lapw%nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,&
                        lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),&
                        lapw%nv(jspin),zMat,nslibd,&
                        nbnd,&
                        stars%rgphs,stars%ustep,stars%ig,(/ 0,0,0 /),&
                        mmn(:,:,ikpt))

                   !---> spherical contribution to mmn0-matrix

                   CALL wann_mmk0_sph(&
                        atoms%llod,noccbd,atoms%nlod,atoms%nat,atoms%ntype,&
                        atoms%lmaxd,atoms%lmax,lmd,atoms%ntype,atoms%neq,&
                        atoms%nlo,atoms%llo,acof(1:noccbd,:,:),&
                        bcof(1:noccbd,:,:),ccof(:,1:noccbd,:,:),&
                        usdus%ddn(:,:,jspin),usdus%uulon(:,:,jspin),&
                        usdus%dulon(:,:,jspin),usdus%uloulopn,&
                        mmn(:,:,ikpt))
                   !---> vacuum contribution to mmn0-matrix

                   IF (input%film .AND. .NOT.oneD%odi%d1) THEN

                      CALL wann_mmk0_vac(&
                           noco%l_noco,atoms%nlotot,qpt_i,&
                           cell%z1,vacuum%nmzd,DIMENSION%nv2d,&
                           stars%mx1,stars%mx2,stars%mx3,&
                           stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,vacuum%delz,&
                           stars%ig2,cell%area,cell%bmat,&
                           cell%bbmat,enpara%evac0(:,jspin),lapw%bkpt,vz(:,:,jspin2),&
                           nslibd,jspin,lapw%k1,lapw%k2,lapw%k3,wannierspin,DIMENSION%nvd,&
                           DIMENSION%nbasfcn,DIMENSION%neigd,zMat,lapw%nv,cell%omtil,&
                           mmn(:,:,ikpt))
                   ELSEIF (oneD%odi%d1) THEN

                      CALL wann_mmk0_od_vac(&
                           DIMENSION, oneD, vacuum, stars, cell,&
                           noco%l_noco,atoms%nlotot,&
                           cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,&
                           stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,&
                           stars%ig,vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,&
                           oneD%odi%n2d,cell%bbmat,enpara%evac0(1,jspin),lapw%bkpt,oneD%odi%M,&
                           oneD%odi%mb,vz(:,1,jspin2),oneD%odi,&
                           nslibd,jspin,lapw%k1,lapw%k2,lapw%k3,wannierspin,DIMENSION%nvd,&
                           cell%area,DIMENSION%nbasfcn,DIMENSION%neigd,zMat,lapw%nv,&
                           stars%sk2,stars%phi2,cell%omtil,qpt_i,&
                           mmn(:,:,ikpt))

                   ENDIF
                ENDIF !l_mmn0


                !---> overlaps with the trial orbitals
                IF (wann%l_projmethod.OR.wann%l_bestproj.OR.wann%l_matrixamn) THEN
                   l_amn2=.FALSE.
                   amnchi = cmplx_1!cmplx(1.,0.)

                   CALL wann_amn (&
                        amnchi,nslibd,nwfs,atoms%ntype,atoms%nlod,atoms%llod,&
                        atoms%llo,atoms%nlo,atoms%lmaxd,atoms%jmtd,lmd,&
                        atoms%neq,atoms%nat,ikpt,nbnd,&
                        atoms%rmsh,atoms%rmt,atoms%jri,atoms%dx,atoms%lmax,&
                        usdus%us(:,:,jspin),usdus%dus(:,:,jspin),&
                        usdus%uds(:,:,jspin),&
                        usdus%duds(:,:,jspin),flo(:,:,:,:,jspin),&
                        ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),acof,bcof,ccof,&
                        (noco%l_soc.OR.noco%l_noco),jspin,&
                        l_amn2,amn(:,:,ikpt))
                   IF(l_amn2)THEN
                      CALL wann_amn (&
                           amnchi,nslibd,nwfs,atoms%ntype,atoms%nlod,atoms%llod,&
                           atoms%llo,atoms%nlo,atoms%lmaxd,atoms%jmtd,lmd,&
                           atoms%neq,atoms%nat,ikpt,nbnd,&
                           atoms%rmsh,atoms%rmt,atoms%jri,atoms%dx,atoms%lmax,&
                           usdus%us(:,:,jspin),usdus%dus(:,:,jspin),&
                           usdus%uds(:,:,jspin),&
                           usdus%duds(:,:,jspin),flo(:,:,:,:,jspin),&
                           ff(:,:,:,:,jspin),gg(:,:,:,:,jspin),acof,bcof,ccof,&
                           (noco%l_soc.OR.noco%l_noco),jspin,&
                           l_amn2,amn(:,:,ikpt),lapw%bkpt)
                   ENDIF
                   !         amn(ikpt,:,:)=amn(ikpt,:,:)*conjg(nsfactor)
                ENDIF



                !****************************************************************
                !...         vanderbilt mmn matrix
                !***************************************************************
                IF (wann%l_matrixmmn .AND.&
                     (.NOT.wann%l_skipkov)) THEN   !  vanderbilt procedure Mmn matrix
                   ALLOCATE ( we_b(DIMENSION%neigd) )

!!! the cycle by the nearest neighbors (nntot) for each kpoint

                   DO   ikpt_b = 1,nntot

                      IF(pair_to_do(ikpt,ikpt_b).EQ.0)CYCLE !save time by symmetry
                      kptibz_b=bpt(ikpt_b,ikpt)
                      IF(wann%l_bzsym) oper_b=mapkoper(kptibz_b)
                      IF (wann%l_bzsym) kptibz_b=irreduc(kptibz_b)



                      ! now we need the wavefunctions for k_b kpoint

                      !        print*,"something to do"


                      n_start=1
                      n_end=DIMENSION%neigd
                      call lapw_b%init(input,noco,kpts,atoms,sym,kptibz_b,cell,(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco.and.mpi%n_size==1),mpi)
                      CALL cdn_read(&
                           eig_id,&
                           DIMENSION%nvd,input%jspins,mpi%irank,mpi%isize, &!wannierspin instead of DIMENSION%jspd?&
                           kptibz_b,jspin,DIMENSION%nbasfcn,&
                           noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,n_end,&
                           nbands_b,eigg,zzMat)


                      nslibd_b = 0

                      eig_b(:) = 0.

                      DO i = 1,nbands_b
                         IF((eigg(i).GE.sliceplot%e1s.AND.nslibd_b.LT.numbands&
                              .AND.wann%l_bynumber).OR.&
                              (eigg(i).GE.sliceplot%e1s.AND.eigg(i).LE.sliceplot%e2s.AND.&
                              wann%l_byenergy).OR.(i.GE.wann%band_min(jspin).AND.&
                              (i.LE.wann%band_max(jspin)).AND.&
                              wann%l_byindex))THEN
                            nslibd_b = nslibd_b + 1
                            eig_b(nslibd_b) = eigg(i)
                            we_b(nslibd_b) = we_b(i)
                            IF (zzMat%l_real) THEN
                               zMat_b%data_r(:,nslibd_b) = zzMat%data_r(:,i)
                            ELSE
                               zMat_b%data_c(:,nslibd_b) = zzMat%data_c(:,i)
                            END IF
                         ENDIF
                      ENDDO

                      !***********************************************************
                      !              Rotate the wavefunction of next neighbor.
                      !***********************************************************
                      IF (wann%l_bzsym .AND. (oper_b.NE.1)  ) THEN
                         !         call wann_kptsrotate(
                         !     >            atoms%nat,atoms%nlod,atoms%llod,
                         !     >            atoms%ntype,atoms%nlo,atoms%llo,atoms%invsat,
                         !     >            noco%l_noco,noco%l_soc,
                         !     >            atoms%ntype,atoms%neq,atoms%nlotot,
                         !     >            kveclo_b,jspin,
                         !     >            oper_b,sym%nop,sym%mrot,DIMENSION%nvd,
                         !     >            nv_b,
                         !     >            shiftkpt(:,bpt(ikpt_b,ikpt)),
                         !     >            sym%tau,
                         !     x            bkpt_b,k1_b(:,:),
                         !     x            k2_b(:,:),k3_b(:,:),
                         !     x            zMat_b,nsfactor_b)
                      ELSE
                         nsfactor_b=CMPLX(1.0,0.0)
                      ENDIF
                      !      print*,"kpt2=",bkpt_b

                      noccbd_b = nslibd_b

                      !cccc   we start with the Mmn matrix   ccccccccccccc

!!! matrix elements of the interstitial overlap
!!! matrix with the weights given by products of
!!! the c-coeff. for different G-vectors, bands and k-points
!!! Mmn(k,b)(IR) = \sum(G,G')C_G^(k,n)*C_G'^(k+b,m)\theta_(G-G')

                      !ccccccccccc  Spherical Contributions         ccccccccccccc

                      ALLOCATE (acof_b(noccbd_b,0:lmd,atoms%nat),&
                           bcof_b(noccbd_b,0:lmd,atoms%nat),&
                           ccof_b(-atoms%llod:atoms%llod,noccbd_b,atoms%nlod,&
                           atoms%nat))

!!! get the band-dependent k-dependent ab coeff.

                      CALL abcof(input,atoms,sym,cell,lapw_b,&
                           noccbd_b,usdus,noco,jspin,oneD,&
                           acof_b,bcof_b,ccof_b,zMat_b)


                      CALL wann_abinv(atoms,&
                           acof_b,bcof_b,ccof_b)

                      !cccccccc  Interstitial  ccccccccccccccccccccccccc
!!! matrix elements of the interstitial overlap
!!! matrix with the weights given by products of
!!! the c-coeff. for different G-vectors, bands and k-points
!!! Mmn(k,b)(IR) = \sum(G,G')C_G^(k,n)*C_G'^(k+b,m)\theta_(G-G')
                      !... these overlaps are the same as in the mmk0 case, only differ
                      !... by the b-dependence of the C-coefficients
                      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

                      addnoco=0
                      addnoco2=0
                      IF(noco%l_noco.AND.(jspin.EQ.2))THEN
                         addnoco  = lapw%nv(1)   + atoms%nlotot
                         addnoco2 = lapw_b%nv(1) + atoms%nlotot
                      ENDIF

                      CALL wann_mmkb_int(&
                           cmplx_1,addnoco,addnoco2,&
                           DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,&
                           stars%ng3,lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),&
                           lapw%nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,&
                           lapw_b%k1(:,jspin),lapw_b%k2(:,jspin),lapw_b%k3(:,jspin),&
                           lapw_b%nv(jspin),zMat_b,nslibd_b,&
                           nbnd,&
                           stars%rgphs,stars%ustep,stars%ig,gb(:,ikpt_b,ikpt),&
                           mmnk(:,:,ikpt_b,ikpt))



                      !ccccccccccc  Spherical Contributions  ccccccccccccc

                      chi = cmplx_1
                      CALL wann_mmkb_sph(&
                           nbnd,atoms%llod,nslibd,nslibd_b,atoms%nlod,atoms%nat,&
                           atoms%ntype,lmd,atoms%jmtd,atoms%taual,sym%nop,atoms%lmax,&
                           atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,acof,bcof,ccof,&
                           lapw_b%bkpt,acof_b,bcof_b,ccof_b,gb(:,ikpt_b,ikpt),lapw%bkpt,&
                           ujug,ujdg,&
                           djug,djdg,ujulog,djulog,ulojug,ulojdg,ulojulog,kdiff,&
                           nntot,chi,&
                           mmnk(:,:,ikpt_b,ikpt))




                      !...vacuum contributions
                      IF (input%film .AND. .NOT.oneD%odi%d1) THEN

                         CALL wann_mmkb_vac(&
                              cmplx_1,noco%l_noco,atoms%nlotot,qpt_i,&
                              nbnd,cell%z1,vacuum%nmzd,DIMENSION%nv2d,&
                              stars%mx1,stars%mx2,stars%mx3,&
                              stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,&
                              vacuum%delz,stars%ig2,cell%area,cell%bmat,&
                              cell%bbmat,enpara%evac0(:,jspin),&
                              enpara%evac0(:,jspin_b),&
                              lapw%bkpt,lapw_b%bkpt,vz(:,:,jspin2),vz(:,:,jspin2_b),&
                              nslibd,nslibd_b,jspin,jspin_b,&
                              lapw%k1,lapw%k2,lapw%k3,lapw_b%k1,lapw_b%k2,lapw_b%k3,&
                              wannierspin,DIMENSION%nvd,&
                              DIMENSION%nbasfcn,DIMENSION%neigd,zMat,zMat_b,&
                              lapw%nv,lapw_b%nv,cell%omtil,&
                              gb(:,ikpt_b,ikpt),&
                              mmnk(:,:,ikpt_b,ikpt))
                      ELSEIF (oneD%odi%d1) THEN

                         CALL wann_mmkb_od_vac(&
                              DIMENSION,oneD,vacuum,stars,cell,&
                              cmplx_1,noco%l_noco,atoms%nlotot,&
                              nbnd,cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,&
                              stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,&
                              stars%ig,vacuum%nmzxy,&
                              vacuum%nmz,vacuum%delz,stars%ig2,oneD%odi%n2d,&
                              cell%bbmat,enpara%evac0(1,jspin),enpara%evac0(1,jspin_b),&
                              lapw%bkpt,lapw_b%bkpt,oneD%odi%M,oneD%odi%mb,&
                              vz(:,1,jspin2),vz(:,1,jspin2_b),oneD%odi,&
                              nslibd,nslibd_b,jspin,jspin_b,lapw%k1,lapw%k2,lapw%k3,lapw_b%k1,lapw_b%k2,lapw_b%k3,&
                              wannierspin,DIMENSION%nvd,cell%area,DIMENSION%nbasfcn,&
                              DIMENSION%neigd,&
                              zMat,zMat_b,lapw%nv,lapw_b%nv,stars%sk2,stars%phi2,cell%omtil,&
                              gb(:,ikpt_b,ikpt),qpt_i,&
                              .FALSE.,1,&
                              mmnk(:,:,ikpt_b,ikpt))
                      ENDIF

                      !        mmnk(:,:,ikpt_b,ikpt)=
                      !                mmnk(:,:,ikpt_b,ikpt)*nsfactor*conjg(nsfactor_b)


                      DEALLOCATE ( acof_b,bcof_b,ccof_b )


                   enddo ! end of loop by the nearest k-neighbors

                   DEALLOCATE ( we_b )

                ENDIF!l_matrixmmn=.true.

9900            CONTINUE ! jump for doublespin loop

                IF (wann%l_matrixmmn) THEN   !  vanderbilt procedure Mmn matrix

                   !*******************************************c 
                   !          START Q-NEIGHBOR LOOP            c
                   !*******************************************c      
                   ALLOCATE ( we_qb(DIMENSION%neigd) )

                   DO iqpt_b=1,nntot_q 
                      IF(.NOT.l_gwf) EXIT              ! old functionality

                      qptibz_b = bpt_q(iqpt_b,iqpt)
                      IF(qptibz_b.EQ.qptibz) CYCLE     ! no need to compute overlaps
                      ! with periodic images for now

                      qptb_i = noco%qss
                      alphb_i = noco%alph
                      betab_i = noco%beta
                      thetab_i = noco%theta
                      phib_i = noco%phi
                      IF(wann%l_sgwf) THEN
                         qptb_i(:) = wann%param_vec(:,qptibz_b)
                         alphb_i(:) = wann%param_alpha(:,qptibz_b)
                      ELSEIF(wann%l_socgwf) THEN
                         IF(wann%l_dim(2)) phib_i = tpi_const*wann%param_vec(2,qptibz_b)
                         IF(wann%l_dim(3)) thetab_i = tpi_const*wann%param_vec(3,qptibz_b)
                      ENDIF

                      !if(pair_to_do_q(iqpt,iqpt_b).eq.0)cycle    ! TODO: use symmetry
                      IF(wann%l_bzsym) oper_qb=mapqoper(qptibz_b)
                      IF (wann%l_bzsym) qptibz_b=irreduc_q(qptibz_b)

                      n_start=1
                      n_end=DIMENSION%neigd

                      ! read in diagonalization information from corresponding
                      ! eig file to q-point iqpt_b at a given k-point ikpt.
                      ! as a check verify that bkpt.eq.bqpt (same k).
                      ! moreover, the plane-wave vectors G(k,q+b) are stored
                      ! in (k1_qb,k2_qb,k3_qb) for later use.

                      CALL lapw_qb%init(input,noco,kpts,atoms,sym,kptibz,cell,(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco.and.mpi%n_size==1),mpi)
                      CALL cdn_read(&
                           innerEig_idList(iqpt_b),&
                           DIMENSION%nvd,input%jspins,mpi%irank,mpi%isize, &!wannierspin instead of DIMENSION%jspd? !kptibz_b2?&
                           kptibz,jspin_b,DIMENSION%nbasfcn,&
                           noco%l_ss,noco%l_noco,DIMENSION%neigd,n_start,&
                           n_end,&
                           nbands_qb,eigg,   &
                           zzMat)


                      ! are we dealing with the same k-point at which
                      ! we want to construct A,B,C coefficients etc. ?
                      IF(ANY(bqpt.NE.lapw%bkpt)) CALL juDFT_error("bqpt.ne.bkpt",&
                           calledby="wannier")

                      zMat_qb%l_real = zzMat%l_real
                      zMat_qb%matsize1 = zzMat%matsize1
                      zMat_qb%matsize2 = zzMat%matsize2
                      IF (zzMat%l_real) THEN
                         ALLOCATE (zMat_qb%data_r(zMat%matsize1,zMat%matsize2))
                         zMat_qb%data_r = 0.0
                      ELSE
                         ALLOCATE (zMat_qb%data_c(zMat%matsize1,zMat%matsize2))
                         zMat_qb%data_c = CMPLX(0.0,0.0)
                      END IF

                      eig_qb(:) = 0.

                      nslibd_qb = 0
                      DO i = 1,nbands_qb
                         IF((eigg(i).GE.sliceplot%e1s.AND.nslibd_qb.LT.numbands&
                              .AND.wann%l_bynumber).OR.&
                              (eigg(i).GE.sliceplot%e1s.AND.eigg(i).LE.sliceplot%e2s.AND.&
                              wann%l_byenergy).OR.(i.GE.wann%band_min(jspin).AND.&
                              (i.LE.wann%band_max(jspin)).AND.wann%l_byindex)) THEN
                            nslibd_qb = nslibd_qb + 1
                            eig_qb(nslibd_qb) = eigg(i)
                            we_qb(nslibd_qb) = we_qb(i)
                            IF (zzMat%l_real) THEN
                               zMat_qb%data_r(:,nslibd_qb) = zzMat%data_r(:,i)
                            ELSE
                               zMat_qb%data_c(:,nslibd_qb) = zzMat%data_c(:,i)
                            END IF
                         ENDIF
                      ENDDO

                      ! check that eigenvectors and -values are identical if q=q+b
                      IF(iqpt.EQ.qptibz_b .AND. jspin.EQ.jspin_b) THEN
                         IF(zMat%l_real) THEN
                            IF(ANY(zMat%data_r.NE.zMat_qb%data_r)) &
                                 WRITE(*,*)'z.ne.z_qb',iqpt,ikpt
                         ELSE
                            IF(ANY(zMat%data_c.NE.zMat_qb%data_c)) &
                                 WRITE(*,*)'z.ne.z_qb',iqpt,ikpt
                         END IF
                         IF(ANY(eig.NE.eig_qb)) WRITE(*,*)'eig.ne.eiq_qb',iqpt,ikpt
                         IF(lapw%nv(jspin).NE.lapw_qb%nv(jspin)) WRITE(*,*)'nv!=nv_qb',iqpt,ikpt
                      ENDIF

                      ! check that number of bands are the same at (k,q) and (k,q+b)
                      IF(nslibd.NE.nslibd_qb)&
                           WRITE(*,*)'nslibd.ne.nslibd_qb',ikpt,iqpt,iqpt_b


                      noccbd_qb = nslibd_qb
                      nsfactor_b=CMPLX(1.0,0.0)


                      ALLOCATE (acof_qb(noccbd_qb,0:lmd,atoms%nat),&
                           bcof_qb(noccbd_qb,0:lmd,atoms%nat),&
                           ccof_qb(-atoms%llod:atoms%llod,noccbd_qb,atoms%nlod,&
                           atoms%nat))

                      acof_qb(:,:,:) = CMPLX(0.,0.) ; bcof_qb(:,:,:) = CMPLX(0.,0.)
                      ccof_qb(:,:,:,:) = CMPLX(0.,0.)

                      ! construct the A,B,C coefficients of the wave function
                      ! at the point (k,q+b) using previously read information

                      CALL abcof(input,atoms,sym,cell,lapw_qb,&
                           noccbd_qb,usdus,noco,jspin_b,oneD,&
                           acof_qb,bcof_qb,ccof_qb,zMat_qb)

                      CALL wann_abinv(atoms,&
                           acof_qb,bcof_qb,ccof_qb)

                      ! check that A,B,C coefficients are the same if q+b = q
                      IF(l_gwf.AND.(iqpt.EQ.qptibz_b).AND.(jspin.EQ.jspin_b)) THEN
                         IF(ANY(acof_qb.NE.acof)) WRITE(*,*)'acof',iqpt,ikpt
                         IF(ANY(bcof_qb.NE.bcof)) WRITE(*,*)'bcof',iqpt,ikpt
                         IF(ANY(ccof_qb.NE.ccof)) WRITE(*,*)'ccof',iqpt,ikpt
                      ENDIF

                      addnoco=0
                      addnoco2=0
                      IF(noco%l_noco.AND.(jspin.EQ.2))THEN
                         addnoco  = lapw%nv(1)   + atoms%nlotot
                      ENDIF
                      IF(noco%l_noco.AND.(jspin_b.EQ.2))THEN
                         addnoco2 = lapw_qb%nv(1) + atoms%nlotot
                      ENDIF


                      ! set up local->global transformation for overlaps
                      DO n=1,atoms%ntype
                         IF(wann%l_sgwf) THEN
                            dalph = alph_i(n)-alphb_i(n)
                            db1 = beta_i(n)/2.
                            db2 = betab_i(n)/2.
                         ELSEIF(wann%l_socgwf) THEN
                            dalph = phi_i-phib_i
                            db1 = theta_i/2.
                            db2 = thetab_i/2.
                         ENDIF
                         coph = COS(dalph)    
                         siph = SIN(dalph)     
                         phasfac = CMPLX(coph,siph)
                         phasfac2= CMPLX(coph,-siph)

                         IF(l_p0 .AND. dalph.NE.0.0) THEN
                            WRITE(*,*)'WARNING: include dalph in chi trafo!'
                         ENDIF

                         IF( (jspin.EQ.1) .AND. (jspin_b.EQ.1) ) THEN ! uu
                            !              chi(n) = cos(db1)*cos(db2)*phasfac
                            !     >               + sin(db1)*sin(db2)*phasfac2
                            chi(n) = COS(db2-db1)
                         ELSEIF( (jspin.EQ.2) .AND. (jspin_b.EQ.2) ) THEN !dd
                            !              chi(n) = cos(db1)*cos(db2)*phasfac2
                            !     >               + sin(db1)*sin(db2)*phasfac
                            chi(n) = COS(db2-db1)
                         ELSEIF( (jspin.EQ.1) .AND. (jspin_b.EQ.2) ) THEN ! ud
                            !              chi(n) = sin(db1)*cos(db2)*phasfac2
                            !     >               - cos(db1)*sin(db2)*phasfac
                            chi(n) = -SIN(db2-db1)
                         ELSEIF( (jspin.EQ.2) .AND. (jspin_b.EQ.1) ) THEN ! du
                            !              chi(n) = cos(db1)*sin(db2)*phasfac2
                            !     >               - sin(db1)*cos(db2)*phasfac
                            chi(n) = SIN(db2-db1)
                         ELSE
                            STOP 'problem setting up chi: jspin,jspin_b'
                         ENDIF
                      ENDDO
                      chi = CONJG(chi)
                      !chi = cmplx_1

                      ! optional: disable chi transformation
                      ! instead of computing overlap w.r.t. global frame
                      ! only consider wave function overlaps in local frames
                      IF(l_nochi) THEN
                         IF(doublespin.LT.3) THEN
                            chi = CMPLX(1.0,0.0)
                         ELSE
                            chi = CMPLX(0.0,0.0)
                         ENDIF
                      ENDIF

                      IF((iqpt.EQ.1).AND.(ikpt.EQ.1).AND.(iqpt_b.EQ.1)) THEN
                         WRITE(*,*)'dbs',doublespin,'chi',chi(1)
                      ENDIF

                      ! muffin tin contribution to overlap is computed taking into account
                      ! the spin-dependent phase exp(+/- tau*b/2) and the ujugaunt integrals
                      ! calculated especially for the q-points before. then, q and q+b
                      ! take the role of k and k+b and the same for G(q+b) and G(k+b)
                      IF(wann%l_sgwf) CALL wann_mmkb_sph( &        
                           nbnd,atoms%llod,nslibd,nslibd_qb,atoms%nlod,atoms%nat,&
                           atoms%ntype,lmd,atoms%jmtd,sign_q*atoms%taual/2.0,sym%nop,&
                           atoms%lmax,atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,&
                           acof,bcof,ccof,qptb_i,&
                           acof_qb,bcof_qb,ccof_qb,gb_q(:,iqpt_b,iqpt),qpt_i,&
                           ujug_q,ujdg_q,&
                           djug_q,djdg_q,ujulog_q,djulog_q,ulojug_q,ulojdg_q,&
                           ulojulog_q,qdiff,       &
                           nntot_q,chi,&
                           mmnk_q(:,:,iqpt_b,ikpt))
                      IF(wann%l_socgwf) CALL wann_mmkb_sph( &        
                           nbnd,atoms%llod,nslibd,nslibd_qb,atoms%nlod,atoms%nat,&
                           atoms%ntype,lmd,atoms%jmtd,&
                           zero_taual,sym%nop,atoms%lmax,         &                
                           atoms%ntype,atoms%neq,atoms%nlo,atoms%llo,acof,bcof,ccof,&
                           (/ 0.0, phib_i/tpi_const, thetab_i/tpi_const /),&
                           acof_qb,bcof_qb,ccof_qb,gb_q(:,iqpt_b,iqpt),&
                           (/ 0.0, phi_i/tpi_const, theta_i/tpi_const /),&
                           ujug_q,ujdg_q,&
                           djug_q,djdg_q,ujulog_q,djulog_q,ulojug_q,ulojdg_q,&
                           ulojulog_q,qdiff,       &
                           nntot_q,chi,&
                           mmnk_q(:,:,iqpt_b,ikpt))


                      IF(((doublespin.NE.3).AND.(doublespin.NE.4))&
                           .OR.(.NOT.noco%l_noco)) THEN

                         IF(.NOT.noco%l_noco)THEN
                            interchi=chi(1)
                            vacchi=chi(1)
                         ELSE
                            interchi=cmplx_1
                            vacchi=cmplx_1
                         ENDIF

                         ! interstitial contribution to overlap is computed using
                         ! (-/+ 1)*G(q+b)/2 as G-vector connecting neighbors
                         ! and lattice vectors G(k,q) (k1...) and G(k,q+b) (k1_qb...)
                         IF(wann%l_sgwf) CALL wann_mmkb_int(&
                              interchi,addnoco,addnoco2,&
                              DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,&
                              stars%ng3,lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),&
                              lapw%nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,&
                              lapw_qb%k1(:,jspin_b),lapw_qb%k2(:,jspin_b),lapw_qb%k3(:,jspin_b),&
                              lapw_qb%nv(jspin_b),zMat_qb,nslibd_qb,&
                              nbnd,&
                              stars%rgphs,stars%ustep,stars%ig,&
                              sign_q*gb_q(:,iqpt_b,iqpt)/2,   &   
                              mmnk_q(:,:,iqpt_b,ikpt))                        
                         IF(wann%l_socgwf) CALL wann_mmkb_int(&
                              interchi,addnoco,addnoco2,&
                              DIMENSION%nvd,stars%mx1,stars%mx2,stars%mx3,&
                              stars%ng3,lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),&
                              lapw%nv(jspin),DIMENSION%neigd,DIMENSION%nbasfcn,zMat,nslibd,&
                              lapw_qb%k1(:,jspin_b),lapw_qb%k2(:,jspin_b),lapw_qb%k3(:,jspin_b),&
                              lapw_qb%nv(jspin_b),zMat_qb,nslibd_qb,&
                              nbnd,&
                              stars%rgphs,stars%ustep,stars%ig,(/ 0, 0, 0 /),  &    
                              mmnk_q(:,:,iqpt_b,ikpt))!m_int(:,:,iqpt_b,ikpt))      


                         ! vacuum contribution in film calculation
                         IF (input%film .AND. .NOT.oneD%odi%d1) THEN
                            IF(wann%l_sgwf) CALL wann_mmkb_vac(&                   
                                 vacchi,noco%l_noco,atoms%nlotot,sign_q*2.*lapw%bkpt,&
                                 nbnd,cell%z1,vacuum%nmzd,DIMENSION%nv2d,&
                                 stars%mx1,stars%mx2,stars%mx3,&
                                 stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,&
                                 vacuum%delz,stars%ig2,cell%area,cell%bmat,&
                                 cell%bbmat,enpara%evac0(:,jspin),enpara%evac0(:,jspin_b),&
                                 sign_q*qpt_i/2.,&
                                 sign_q*qptb_i/2.,&
                                 vz(:,:,jspin2),vz(:,:,jspin2_b),&
                                 nslibd,nslibd_qb,jspin,jspin_b,&
                                 lapw%k1,lapw%k2,lapw%k3,lapw_qb%k1,lapw_qb%k2,lapw_qb%k3,&
                                 wannierspin,DIMENSION%nvd,&
                                 DIMENSION%nbasfcn,DIMENSION%neigd,zMat,zMat_qb,lapw%nv,&
                                 lapw_qb%nv,cell%omtil,&
                                 sign_q*gb_q(:,iqpt_b,iqpt)/2,&
                                 mmnk_q(:,:,iqpt_b,ikpt))
                            IF(wann%l_socgwf) CALL wann_mmkb_vac(&  
                                 vacchi,noco%l_noco,atoms%nlotot,qpt_i,&
                                 nbnd,cell%z1,vacuum%nmzd,DIMENSION%nv2d,&
                                 stars%mx1,stars%mx2,stars%mx3,&
                                 stars%ng3,vacuum%nvac,stars%ig,vacuum%nmz,&
                                 vacuum%delz,stars%ig2,cell%area,cell%bmat,&
                                 cell%bbmat,enpara%evac0(:,jspin),enpara%evac0(:,jspin_b),&
                                 bqpt,bqpt,&
                                 vz(:,:,jspin2),vz(:,:,jspin2_b),&
                                 nslibd,nslibd_qb,jspin,jspin_b,&
                                 lapw%k1,lapw%k2,lapw%k3,lapw_qb%k1,lapw_qb%k2,lapw_qb%k3,&
                                 wannierspin,DIMENSION%nvd,&
                                 DIMENSION%nbasfcn,DIMENSION%neigd,zMat,zMat_qb,lapw%nv,&
                                 lapw_qb%nv,cell%omtil,&
                                 (/ 0, 0, 0 /),&
                                 mmnk_q(:,:,iqpt_b,ikpt))

                            ! vacuum contribution in one-dimensional chain calculation where
                            ! q-point plays role of k-point with proper prefactor of (-/+ 1/2).
                            ! moreover, the k-point ikpt is treated like qss in the subroutine
                            ! such that a correction for sign and factor has to appear as well.
                            ! lattice vectors G(k,q) (k1...) and G(k,q+b) (k1_qb...) need to 
                            ! be provided.
                         ELSEIF (oneD%odi%d1) THEN
                            IF(wann%l_sgwf) CALL wann_mmkb_od_vac(&
                                 DIMENSION,oneD,vacuum,stars,cell,&
                                 vacchi,noco%l_noco,atoms%nlotot, &          
                                 nbnd,cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,&
                                 stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,&
                                 stars%ig,vacuum%nmzxy,&
                                 vacuum%nmz,vacuum%delz,stars%ig2,oneD%odi%n2d,&
                                 cell%bbmat,enpara%evac0(1,jspin),enpara%evac0(1,jspin_b),&
                                 sign_q*qpt_i/2.,sign_q*qptb_i/2.,&
                                 oneD%odi%M,oneD%odi%mb,&
                                 vz(:,1,jspin2),vz(:,1,jspin2_b),oneD%odi,&
                                 nslibd,nslibd_qb,jspin,jspin_b,&
                                 lapw%k1,lapw%k2,lapw%k3,lapw_qb%k1,lapw_qb%k2,lapw_qb%k3,&
                                 wannierspin,DIMENSION%nvd,cell%area,DIMENSION%nbasfcn,&
                                 DIMENSION%neigd,zMat,zMat_qb,lapw%nv,lapw_qb%nv,stars%sk2,&
                                 stars%phi2,cell%omtil,&
                                 sign_q*gb_q(:,iqpt_b,iqpt)/2,sign_q*2.*lapw%bkpt, &
                                 .TRUE.,sign_q,&
                                 mmnk_q(:,:,iqpt_b,ikpt))
                            IF(wann%l_socgwf) CALL wann_mmkb_od_vac(   &
                                 DIMENSION,oneD,vacuum,stars,cell,&                
                                 vacchi,noco%l_noco,atoms%nlotot,&           
                                 nbnd,cell%z1,vacuum%nmzxyd,vacuum%nmzd,DIMENSION%nv2d,&
                                 stars%mx1,stars%mx2,stars%mx3,stars%ng2,stars%ng3,&
                                 stars%ig,vacuum%nmzxy,&
                                 vacuum%nmz,vacuum%delz,stars%ig2,oneD%odi%n2d,&
                                 cell%bbmat,enpara%evac0(1,jspin),enpara%evac0(1,jspin_b),&
                                 bqpt,bqpt,&
                                 oneD%odi%M,oneD%odi%mb,&
                                 vz(:,1,jspin2),vz(:,1,jspin2_b),oneD%odi,&
                                 nslibd,nslibd_qb,jspin,jspin_b,&
                                 lapw%k1,lapw%k2,lapw%k3,lapw_qb%k1,lapw_qb%k2,lapw_qb%k3,&
                                 wannierspin,DIMENSION%nvd,cell%area,DIMENSION%nbasfcn,&
                                 DIMENSION%neigd,zMat,zMat_qb,lapw%nv,lapw_qb%nv,stars%sk2,&
                                 stars%phi2,cell%omtil,&
                                 (/ 0, 0, 0 /),qpt_i, &
                                 .FALSE.,1,&
                                 mmnk_q(:,:,iqpt_b,ikpt))
                         ENDIF!film resp. odi


                      ENDIF!doublespin

                      DEALLOCATE ( acof_qb,bcof_qb,ccof_qb )


                   ENDDO !iqpt_b, q-neighbors
                   !**************************************************c
                   !              END Q-NEIGHBOR LOOP                 c
                   !**************************************************c     

                   DEALLOCATE ( we_qb )

                ENDIF    !   if wann%l_matrixmmn = true

                IF(.NOT.wann%l_bzsym)oper=0
                IF(.NOT.wann%l_plot_symm.OR.oper.EQ.1)THEN
                   IF (wann%l_wann_plot .AND. &
                        (doublespin.EQ.1 .OR. doublespin.EQ.2)) THEN

                      addnoco=0
                      IF(noco%l_noco.AND.(jspin.EQ.2))THEN
                         addnoco=lapw%nv(1)+atoms%nlotot
                      ENDIF

                      IF (sliceplot%slice) THEN
                         IF (ikpt.EQ.sliceplot%kk) THEN

                            WRITE (6,*) 'nnne=',sliceplot%nnne
                            WRITE (6,*) 'eig(nnne)=',eig(sliceplot%nnne)
                            WRITE (6,*) 'we(nnne)=',we(sliceplot%nnne)

                            CALL wann_plot(&
                                 DIMENSION,oneD,vacuum,stars,cell,atoms,&
                                 DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,&
                                 vacuum%nmzxyd,&
                                 stars%ng2,sphhar%ntypsd,atoms%ntype,atoms%lmaxd,&
                                 atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,atoms%neq,&
                                 stars%ng3,vacuum%nvac,vacuum%nmz,vacuum%nmzxy,stars%ng2,&
                                 sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,&
                                 sym%symor,sym%invs,sym%invs2,cell%z1,vacuum%delz,&
                                 atoms%ngopr,atoms%ntypsy,atoms%jri,atoms%pos,atoms%zatom,&
                                 atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,&
                                 cell%amat,cell%bmat,cell%bbmat,ikpt,sliceplot%nnne,&
                                 sliceplot%kk,DIMENSION%nvd,atoms%nlod,atoms%llod,&
                                 lapw%nv(jspin),lmd,lapw%bkpt,cell%omtil,atoms%nlo,atoms%llo,&
                                 lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),enpara%evac0(:,jspin),&
                                 vz(:,:,jspin2),&
                                 nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,&
                                 ff(:,:,:,:,jspin),&
                                 gg(:,:,:,:,jspin),flo,acof,bcof,ccof,zMat,&
                                 stars%mx1,stars%mx2,stars%mx3,stars%ig,stars%ig2,&
                                 stars%sk2,stars%phi2,&
                                 noco%l_noco,noco%l_ss,qpt_i,&
                                 addnoco,get_index_kq(ikpt,iqpt,fullnkpts),wann%l_sgwf)

                         ENDIF
                      ELSE ! not sliceplot%slice

                         CALL wann_plot(&
                              DIMENSION,oneD,vacuum,stars,cell,atoms,&
                              DIMENSION%nv2d,jspin,oneD%odi,oneD%ods,stars%ng3,&
                              vacuum%nmzxyd,&
                              stars%ng2,sphhar%ntypsd,atoms%ntype,atoms%lmaxd,&
                              atoms%jmtd,atoms%ntype,atoms%nat,vacuum%nmzd,atoms%neq,&
                              stars%ng3,vacuum%nvac,vacuum%nmz,vacuum%nmzxy,stars%ng2,&
                              sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,&
                              sym%symor,sym%invs,sym%invs2,cell%z1,vacuum%delz,&
                              atoms%ngopr,atoms%ntypsy,atoms%jri,atoms%pos,atoms%zatom,&
                              atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,&
                              cell%amat,cell%bmat,cell%bbmat,ikpt,sliceplot%nnne,&
                              sliceplot%kk,DIMENSION%nvd,atoms%nlod,atoms%llod,&
                              lapw%nv(jspin),lmd,lapw%bkpt,cell%omtil,atoms%nlo,atoms%llo,&
                              lapw%k1(:,jspin),lapw%k2(:,jspin),lapw%k3(:,jspin),enpara%evac0(:,jspin),&
                              vz(:,:,jspin2),&
                              nslibd,DIMENSION%nbasfcn,DIMENSION%neigd,&
                              ff(:,:,:,:,jspin),&
                              gg(:,:,:,:,jspin),flo,acof,bcof,ccof,zMat,&
                              stars%mx1,stars%mx2,stars%mx3,stars%ig,stars%ig2,&
                              stars%sk2,stars%phi2,&
                              noco%l_noco,noco%l_ss,qpt_i,&
                              addnoco,get_index_kq(ikpt,iqpt,fullnkpts),wann%l_sgwf)

                         IF(wann%l_plot_symm.AND.wann%l_bzsym)THEN
                            DO kplot=1,fullnkpts
                               IF(irreduc(kplot).EQ.kptibz)THEN
                                  plotoper=mapkoper(kplot) 
                                  IF(plotoper.LT.0)THEN
                                     plotoper=-plotoper
                                     l_conjugate=.TRUE.
                                  ELSE
                                     l_conjugate=.FALSE.
                                  ENDIF
                                  kplotoper=sym%invtab(plotoper)
                                  CALL wann_plot_symm(jspin,sym%mrot(:,:,kplotoper),ikpt,&
                                       kplot,l_conjugate)
                               ENDIF
                            ENDDO
                         ENDIF

                      ENDIF
                   ENDIF
                ENDIF !wann%l_plot_symm
                DEALLOCATE ( acof,bcof,ccof,we,eigg )

                IF(wann%l_projmethod.OR.wann%l_bestproj)THEN
                   CALL wann_projmethod(&
                        fullnkpts,&
                        wann%l_projmethod,wann%l_bestproj,&
                        ikpt,nwfs,nslibd,amn,eig,&
                        psiw,hwfr)
                ENDIF ! projmethod

             ENDIF   ! loop by processors

          ENDDO ! end of cycle by the k-points


          IF(wann%l_matrixmmn)&
               DEALLOCATE(ujug,ujdg,djug,djdg,&
               ujulog,djulog,ulojulog,ulojug,ulojdg)
          IF(wann%l_matrixmmn.AND.l_gwf)&
               DEALLOCATE(ujug_q,ujdg_q,djug_q,&
               djdg_q,ujulog_q,djulog_q,ulojulog_q,ulojug_q,&
               ulojdg_q)

#ifdef CPP_MPI
          CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif


          !******************************************************
          !     Write down the projections.
          !******************************************************
5         CONTINUE

          IF(doublespin.EQ.3 .OR. doublespin.EQ.4) GOTO 912

          IF(wann%l_nabla)THEN
            
             nablamat=nablamat*hescale
             CALL wann_write_nabla(&
                  mpi%mpi_comm,l_p0,spin12(jspin)//'.nabl',&
                  'Matrix elements of nabla operator',&
                  nbnd,fullnkpts,nbnd,&
                  mpi%irank,mpi%isize,&
                  nablamat)
          ENDIF

          IF(wann%l_soctomom)THEN
             soctomom=soctomom*hescale
             CALL wann_write_nabla(&
                  mpi%mpi_comm,l_p0,spin12(jspin)//'.stm',&
                  'Matrix elements of stm operator',&
                  nbnd,fullnkpts,nbnd,&
                  mpi%irank,mpi%isize,&
                  soctomom)
          ENDIF

          IF(wann%l_surfcurr)THEN
             surfcurr=surfcurr*hescale
             CALL wann_write_nabla(&
                  mpi%mpi_comm,l_p0,spin12(jspin)//'.surfcurr',&
                  'Surface currents',&
                  nbnd,fullnkpts,nbnd,&
                  mpi%irank,mpi%isize,&
                  surfcurr)
          ENDIF

          IF((noco%l_soc.OR.noco%l_noco).AND.wann%l_mmn0)THEN
             CALL wann_write_amn(&
                  mpi%mpi_comm,&
                  l_p0,spin12(jspin)//'.socmmn0',&
                  'Overlaps of the wavefunct. at the same kpoint',&
                  nbnd,fullnkpts,nbnd,&
                  mpi%irank,mpi%isize,.FALSE.,&
                  mmn,.FALSE.)
          ENDIF !noco%l_soc and l_mmn0  

          IF(wann%l_orbcomp)THEN
             num_angl=9
             IF(wann%l_oc_f)num_angl=16
             CALL wann_write_matrix5(&
                  mpi%mpi_comm,l_p0,spin12(jspin)//'.orbcomp',&
                  'angular components',&
                  nbnd,nbnd,&
                  num_angl,wann%oc_num_orbs,fullnkpts,&
                  mpi%irank,mpi%isize,&
                  orbcomp)
          ENDIF

          IF((noco%l_soc.OR.noco%l_noco) .AND. (doublespin.EQ.1)) THEN
             IF(wann%l_mmn0) socmmn(:,:,:)=mmn(:,:,:)
             GOTO 912
          ENDIF

          IF(noco%l_soc.OR.noco%l_noco) THEN
             jspin2=1
             IF(wann%l_mmn0)      mmn(:,:,:)=socmmn(:,:,:)+mmn(:,:,:)
             IF(wann%l_mmn0)      DEALLOCATE(socmmn)
          ENDIF

          IF (wann%l_matrixamn)THEN
             CALL wann_write_amn(&
                  mpi%mpi_comm,&
                  l_p0,spin12(jspin2)//'.amn',&
                  'Overlaps of the wavefunct. with the trial orbitals',&
                  nbnd,fullnkpts,nwfs,&
                  mpi%irank,mpi%isize,.FALSE.,&
                  amn(:,:,:),wann%l_unformatted)
          ENDIF !wann%l_matrixamn

          IF(wann%l_anglmom)THEN
             CALL wann_write_matrix4(&
                  mpi%mpi_comm,&
                  l_p0,spin12(jspin2)//'.anglmom',&
                  'Matrix elements of angular momentum',&
                  nbnd,nbnd,3,fullnkpts,&
                  mpi%irank,mpi%isize,&
                  anglmom)
          ENDIF

          IF (l_proj) THEN
             !**************************************************************
             !            for projmethod: write down WF1.umn
             !*************************************************************
             IF((wann%l_projmethod.OR.wann%l_bestproj))THEN
                CALL wann_write_amn(&
                     mpi%mpi_comm,l_p0,spin12(jspin2)//'.umn',&
                     'transformation to first guess Wannier functions',&
                     nbnd,fullnkpts,nwfs,&
                     mpi%irank,mpi%isize,.FALSE.,&
                     psiw,.FALSE.)
#ifdef CPP_MPI
                ALLOCATE( hwfr2(SIZE(hwfr,1),SIZE(hwfr,2)) )
                length=nwfs*nwfs
                CALL MPI_REDUCE(&
                     hwfr,hwfr2,length,&
                     CPP_MPI_COMPLEX,MPI_SUM,0,&
                     mpi%mpi_comm,ierr)
                hwfr=hwfr2
                DEALLOCATE(hwfr2)
#endif
                !********************************************************
                !        projmethod: hamiltonian matrix in real space
                !********************************************************
                IF(l_p0)THEN
                   WRITE (6,*) 'the hamiltonian matrix in real space:'
                   DO i = 1,nwfs
                      DO j = 1,nwfs
                         WRITE (6,*) '   WFs:',i,'and',j
                         WRITE (6,*) '     matrix element:',hwfr(i,j)
                      ENDDO
                   ENDDO
                ENDIF !l_p0
                DEALLOCATE(hwfr)
             ENDIF !wann%l_projmethod or wann%l_bestproj
          ENDIF !l_proj

          !*********************************************************
          !.....write down the mmn0 matrix 
          !*********************************************************

          IF(wann%l_mmn0)THEN
             CALL wann_write_amn(&
                  mpi%mpi_comm,&
                  l_p0,spin12(jspin2)//'.mmn0',&
                  'Overlaps of the wavefunct. at the same kpoint',&
                  nbnd,fullnkpts,nbnd,&
                  mpi%irank,mpi%isize,.FALSE.,&
                  mmn,.FALSE.)
          ENDIF !wann%l_mmn0  


          !*****************************************************
          !.....write down the matrix M^{k,b}_{mn} 
          !*****************************************************

          ! 912   continue
          IF(wann%l_matrixmmn.AND.(.NOT.wann%l_skipkov))THEN
             CALL wann_write_mmnk(&
                  mpi%mpi_comm,jspin2,l_p0,fullnkpts,nntot,wann,&
                  maptopair,pair_to_do,nbnd,bpt,gb,&
                  mpi%isize,mpi%irank,"            ",&      
                  mmnk,wann%l_unformatted)
          ENDIF !wann%l_matrixmmn

912       CONTINUE

          IF(l_gwf .AND. wann%l_matrixmmn) THEN
             !      mmnk_q = mmnk_q + m_sph+m_int+m_vac
             IF(doublespin.EQ.doublespin_max) THEN
                WRITE(fname,'("param_",i4.4,".mmn")')iqpt
                CALL wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,&
                     nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,&
                     mpi%isize,mpi%irank,fname,mmnk_q,&
                     wann%l_unformatted)
             ENDIF

             IF(.FALSE.) THEN
                WRITE(fname,'("param_",i4.4,"_",i1,".mmn")')iqpt,doublespin
                CALL wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,&
                     nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,&
                     mpi%isize,mpi%irank,fname,&
                     m_sph+m_int+m_vac,wann%l_unformatted)

                WRITE(fname,'("param_",i4.4,"_",i1,"_int.mmn")')iqpt,doublespin
                CALL wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,&
                     nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,&
                     mpi%isize,mpi%irank,fname,m_int,&
                     wann%l_unformatted)
                WRITE(fname,'("param_",i4.4,"_",i1,"_sph.mmn")')iqpt,doublespin
                CALL wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,&
                     nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,&
                     mpi%isize,mpi%irank,fname,m_sph,&
                     wann%l_unformatted)
                WRITE(fname,'("param_",i4.4,"_",i1,"_vac.mmn")')iqpt,doublespin
                CALL wann_write_mmnk2(l_p0,fullnkpts,nntot_q,wann,&
                     nbnd,bpt_q(:,iqpt),gb_q(:,:,iqpt)/2,&
                     mpi%isize,mpi%irank,fname,m_vac,&
                     wann%l_unformatted)
             ENDIF
             !      m_int = cmplx(0.,0.)
             !      m_sph = cmplx(0.,0.)
             !      m_vac = cmplx(0.,0.)

          ENDIF


          IF( ALLOCATED (nablamat) ) DEALLOCATE( nablamat )
          IF( ALLOCATED (soctomom) ) DEALLOCATE( soctomom )

          IF( ALLOCATED (surfcurr) ) DEALLOCATE( surfcurr )
          IF( ALLOCATED( mmn ) ) DEALLOCATE(mmn)
          IF( ALLOCATED( amn ) ) THEN
             IF(.NOT.(noco%l_soc.OR.noco%l_noco))DEALLOCATE(amn)
          ENDIF
          IF ( ALLOCATED (psiw) ) DEALLOCATE ( psiw )
          IF (wann%l_matrixmmn) THEN
             IF(.NOT.(noco%l_soc.OR.noco%l_noco))DEALLOCATE (mmnk)
          ENDIF
          IF (wann%l_anglmom .AND. ALLOCATED(anglmom))THEN
             IF(.NOT.(noco%l_soc.OR.noco%l_noco))DEALLOCATE (anglmom)
          ENDIF
          DEALLOCATE(flo)

          IF(.NOT.noco%l_noco)nrec=nrec+nkpts

          !#ifdef CPP_MPI
          !      call MPI_BARRIER(mpi%mpi_comm,ierr)
          !#endif

       ENDDO ! end of cycle by spins


#ifdef CPP_MPI
       CALL MPI_BARRIER(mpi%mpi_comm,ierr)
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

       IF (wann%l_matrixmmn.AND.ALLOCATED(mmnk))THEN 
          DEALLOCATE ( mmnk )
       ENDIF

       IF(ALLOCATED(mmnk_q)) DEALLOCATE(mmnk_q)
       IF(ALLOCATED(m_int)) DEALLOCATE(m_int)
       IF(ALLOCATED(m_sph)) DEALLOCATE(m_sph)
       IF(ALLOCATED(m_vac)) DEALLOCATE(m_vac)

       IF ((wann%l_projmethod.OR.wann%l_bestproj.OR.wann%l_matrixamn)&
            .AND.ALLOCATED(amn))THEN
          DEALLOCATE ( amn )
       ENDIF

       IF(wann%l_anglmom.AND.ALLOCATED(anglmom))THEN
          DEALLOCATE ( anglmom )  
       ENDIF

       IF ((wann%l_projmethod.OR.wann%l_bestproj)&
            .AND.ALLOCATED(hwfr)) THEN
          DEALLOCATE ( hwfr )
       ENDIF

       DEALLOCATE(innerEig_idList)
       
    ENDDO ! iqpt, q-points
       !************************************************c
       !               END Q LOOP                       c
       !************************************************c

       IF(ALLOCATED(pair_to_do))DEALLOCATE(pair_to_do,maptopair)

       DEALLOCATE ( vr,vz)
       DEALLOCATE ( ff,gg ) 
       IF (wann%l_bzsym)  DEALLOCATE(irreduc,mapkoper)
       IF (wann%l_bzsym.AND.l_gwf)  DEALLOCATE(irreduc_q,mapqoper)
       IF(ALLOCATED(pair_to_do_q))&
            DEALLOCATE(pair_to_do_q,maptopair_q)

       IF (ALLOCATED(kdiff)) DEALLOCATE ( kdiff )
       IF (ALLOCATED(qdiff)) DEALLOCATE(qdiff,zero_qdiff)


9110   CONTINUE  ! jump for l_finishgwf=T

       ! correct for previously introduced factor of 2 in the
       ! G-vectors connecting neighbors across the BZ boundary
       IF(wann%l_sgwf) gb_q = gb_q/2 
       IF(wann%l_socgwf) gb_q = gb_q/2

       ! set up input files for wannier90 --> HDWFs
       IF(l_p0.AND.l_gwf.AND.(wann%l_matrixmmn.OR.wann%l_matrixamn)&
            .AND.(.NOT.wann%l_skipkov)) THEN
          CALL wann_gwf_commat(fullnkpts,nntot,bpt,fullnqpts,&
               nntot_q,bpt_q,gb,gb_q,wann%aux_latt_const,&
               wann%l_unformatted,wann%l_matrixamn,&
               wann%l_matrixmmn,wann%l_dim,&
               wann%nparampts,wann%param_vec/2.0)
       ENDIF

       IF(l_p0.AND.l_gwf.AND.wann%l_anglmom) THEN
          CALL wann_gwf_anglmom(fullnkpts,fullnqpts,wann%l_unformatted)
       ENDIF

       IF (wann%l_matrixmmn) THEN
          DEALLOCATE (gb,bpt)
          IF(l_gwf) DEALLOCATE (gb_q,bpt_q)
       ENDIF


1911   CONTINUE

       IF(ALLOCATED(chi)) DEALLOCATE(chi)

       CALL timeStop("Wannier total")
       CALL wann_postproc(&
            DIMENSION,stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,&
            lapw,oneD,noco,cell,vTot,enpara,sliceplot,eig_id,l_real,&                !eig_id is used here after closing the files?!&
            wann,fullnkpts,l_proj,results%ef,wann%l_sgwf,fullnqpts)

#ifdef CPP_MPI
       CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

       IF(.NOT.wann%l_ldauwan) THEN
          DO pc = 1, wann%nparampts
             CALL close_eig(eig_idList(pc))
          END DO

          CALL juDFT_end("wannier good",mpi%irank)
       END IF

     END SUBROUTINE wannier


   END MODULE m_wannier
