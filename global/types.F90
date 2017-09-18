!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types
  USE m_types_rcmat
  USE m_types_xcpot
  !*************************************************************
  !     This module contains definitions for all kind of types
  !*************************************************************
  !
  ! Types for orbital moment calculation:
  !
  !
  ! Types for spin-off-diagonal charge density:
  !
  TYPE t_mt21                          ! 'normal' contributions
     SEQUENCE
     REAL ::  uun,udn,dun,ddn           ! normes of radial overlaps
     COMPLEX :: uu,ud,du,dd             ! values
  END TYPE t_mt21

  TYPE t_lo21                          ! ocal orbitals & (u,d)
     SEQUENCE
     REAL ::  uulon,dulon,uloun,ulodn   ! normes of radial overlaps
     COMPLEX :: uulo,dulo,ulou,ulod     ! values
  END TYPE t_lo21

  TYPE t_orb                           ! 'normal' contributions
     SEQUENCE
     REAL :: uu,dd                      ! z   component
     COMPLEX :: uup,uum,ddp,ddm         ! +/- component
  END TYPE t_orb

  TYPE t_orbl                          ! local orbitals & (u,d)
     SEQUENCE
     REAL :: uulo,dulo
     COMPLEX :: uulop,uulom,dulop,dulom
  END TYPE t_orbl

  TYPE t_orblo                         ! lo,lo' contributions
     SEQUENCE
     REAL :: z
     COMPLEX :: p,m
  END TYPE t_orblo
  TYPE t_lapw
     INTEGER :: nv(2)
     INTEGER :: nv_tot
     INTEGER :: nmat
     INTEGER,ALLOCATABLE:: k1(:,:)
     INTEGER,ALLOCATABLE:: k2(:,:)
     INTEGER,ALLOCATABLE:: k3(:,:)
     INTEGER,ALLOCATABLE:: kp(:,:)
     REAL,ALLOCATABLE::rk(:,:)
  END TYPE t_lapw

  TYPE t_tlmplm
     COMPLEX,ALLOCATABLE :: tdd(:,:,:)
     COMPLEX,ALLOCATABLE :: tdu(:,:,:)
     !(0:lmplmd,ntypd,tspin)
     COMPLEX,ALLOCATABLE :: tud(:,:,:)
     COMPLEX,ALLOCATABLE :: tuu(:,:,:)
     !(0:lmplmd,ntypd,tspin)
     INTEGER,ALLOCATABLE :: ind(:,:,:,:)
     !(0:lmd,0:lmd,ntypd,tspin)
     COMPLEX,ALLOCATABLE :: tdulo(:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuulo(:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuloulo(:,:,:,:)
     !(-llod:llod,-llod:llod,mlolotot,tspin)
  END TYPE t_tlmplm

  TYPE t_usdus
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: us
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: dus
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: uds
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: duds !(0:lmaxd,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ddn  !(0:lmaxd,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ulos
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: dulos
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: uulon
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: dulon !(nlod,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: uloulopn!  (nlod,nlod,ntypd,jspd)
   CONTAINS
     PROCEDURE :: init => usdus_init
  END TYPE t_usdus



  ! types for 1D calculations
  TYPE od_dim
     LOGICAL :: d1
     INTEGER :: mb,M,k3,m_cyl
     INTEGER :: chi,rot
     LOGICAL :: invs,zrfs
     INTEGER :: n2d,nq2,nn2d
     INTEGER :: kimax2
     INTEGER :: nop,nat
  END TYPE od_dim

  TYPE od_inp
     LOGICAL :: d1
     INTEGER :: mb,M,k3,m_cyl
     INTEGER :: chi,rot
     LOGICAL :: invs,zrfs
     INTEGER :: n2d,nq2,nn2d
     INTEGER :: kimax2
     INTEGER, POINTER :: ig(:,:)  !(-k3:k3,-M:M)
     INTEGER, POINTER :: kv(:,:)        !(2,n2d)
     INTEGER, POINTER :: nst2(:)        !(n2d)
  END TYPE od_inp

  TYPE od_sym
     INTEGER :: nop,nat
     INTEGER, POINTER :: ngopr(:)     !(nat)
     REAL   , POINTER :: mrot(:,:,:)  !(3,3,nop)
     REAL   , POINTER :: tau(:,:)     !(3,nop)
     INTEGER, POINTER :: invtab(:)    !(nop)
     INTEGER, POINTER :: multab(:,:)  !(nop,nop)
  END TYPE od_sym

  TYPE od_lda
     INTEGER :: nn2d
     INTEGER, POINTER :: igf(:,:)  !(0:nn2d-1,2)
     REAL   , POINTER :: pgf(:)    !(0:nn2d-1)
  END TYPE od_lda

  TYPE od_gga
     INTEGER          :: nn2d
     REAL, POINTER    :: pgfx(:)  ! (0:nn2d-1)
     REAL, POINTER    :: pgfy(:)
     REAL, POINTER    :: pgfxx(:)
     REAL, POINTER    :: pgfyy(:)
     REAL, POINTER    :: pgfxy(:)
  END TYPE od_gga


  !
  ! Type for LDA+U:
  !
  TYPE t_utype
     SEQUENCE
     REAL u,j         ! the actual U and J parameters
     INTEGER l        ! the l quantum number to which this U parameter belongs
     INTEGER atomType ! The atom type to which this U parameter belongs
     LOGICAL :: l_amf ! logical switch to choose the "around mean field" LDA+U limit
  END TYPE t_utype
  !
  ! Type for the HF total energy
  !
  TYPE t_energy_hf
     REAL :: valence
     REAL :: core
  END TYPE t_energy_hf

  !
  ! Type for the electric field
  !
  TYPE t_efield
     REAL    :: zsigma  = 10.0  ! Distance to the charged plates
     REAL    :: sigma   =  0.0  ! charge at the plates
     REAL    :: sig_b(2)=  0.0  ! Extra charge for the top/bottom plate
     COMPLEX :: vslope  =  0.0  ! Dirichlet bnd. cond.: Slope
     REAL,    ALLOCATABLE :: sigEF(:,:,:) ! (nx, ny, nvac)
     COMPLEX, ALLOCATABLE :: rhoEF(:,:)   ! (g_||, nvac)
     COMPLEX, ALLOCATABLE :: C1(:), C2(:) ! Coeff. for Dirichlet bnd.cond.
     LOGICAL :: l_segmented = .FALSE.
     LOGICAL :: plot_charge = .FALSE. ! Plot charge as inputted
     LOGICAL :: plot_rho    = .FALSE. ! Plot Fourier-transformed charge
     LOGICAL :: autocomp    = .TRUE.  ! Auto-compensate film charge
     LOGICAL :: dirichlet = .FALSE. ! Dirichlet vs. Neumann boundary cond.
     LOGICAL :: l_dirichlet_coeff = .FALSE. ! For MPI, true if C1/C2 set
  END TYPE t_efield


  TYPE t_atoms
     !<no of types
     INTEGER :: ntype
     !<total-no of atoms
     INTEGER :: nat
     !<dimensions of LO's
     INTEGER ::nlod
     INTEGER ::llod
     INTEGER ::nlotot
     !lmaxd=maxval(lmax)
     INTEGER:: lmaxd
     ! no of lda+us
     INTEGER ::n_u
     ! dimensions
     INTEGER :: jmtd
     !No of element
     INTEGER,ALLOCATABLE ::nz(:)
     !atoms per type
     INTEGER,ALLOCATABLE::neq(:)
     !radial grid points
     INTEGER,ALLOCATABLE::jri(:)
     !core states
     INTEGER,ALLOCATABLE::ncst(:)
     !How many states are explicitely provided?
     INTEGER,ALLOCATABLE::numStatesProvided(:)
     !core state occupations
     REAL,ALLOCATABLE::coreStateOccs(:,:,:)
     !core state nprnc
     INTEGER,ALLOCATABLE::coreStateNprnc(:,:)
     !core state kappa
     INTEGER,ALLOCATABLE::coreStateKappa(:,:)
     !lmax
     INTEGER,ALLOCATABLE::lmax(:)
     !lmax non-spherical
     INTEGER,ALLOCATABLE::lnonsph(:)
     !expansion of pseudo-charge
     INTEGER,ALLOCATABLE::ncv(:)
     !no of LO
     INTEGER,ALLOCATABLE::nlo(:)
     !l of LO (nlo,ntype)
     INTEGER,ALLOCATABLE::llo(:,:)
     !lmax for lapw (ntype)
     INTEGER,ALLOCATABLE::lapw_l(:)
     !first LO with a given l (max(nlo
     INTEGER,ALLOCATABLE::lo1l(:,:)
     !??
     INTEGER,ALLOCATABLE::ulo_der(:,:)
     !no of LOs per l (max(nlo1),ntype
     INTEGER,ALLOCATABLE::nlol(:,:)
     !true if LO is formed by \dot u (
     LOGICAL,ALLOCATABLE::l_dulo(:,:)
     !no of op that maps atom into
     INTEGER,ALLOCATABLE::ngopr(:)
     !symetry of atom (nat)
     INTEGER,ALLOCATABLE::ntypsy(:)
     !no of sphhar for atom type(ntype
     INTEGER,ALLOCATABLE ::nlhtyp(:)
     !atom mapped to by inversion (nat
     INTEGER,ALLOCATABLE ::invsat(:)
     !Calaculate forces for this atom?
     LOGICAL,ALLOCATABLE :: l_geo(:)
     !MT-Radius (ntype)
     REAL,ALLOCATABLE::rmt(:)
     !log increment(ntype)
     REAL,ALLOCATABLE::dx(:)
     !vol of MT(ntype)
     REAL,ALLOCATABLE::volmts(:)
     !radial grid points(max(jri),ntyp
     REAL,ALLOCATABLE::rmsh(:,:)
     !charge of nucleus(ntype)
     REAL,ALLOCATABLE::zatom(:)
     !initial mag moment(ntype)
     REAL,ALLOCATABLE::bmu(:)
     !pos of atom (absol) (3,nat)
     REAL,ALLOCATABLE::pos(:,:)
     !pos of atom (relat)(3,nat)
     REAL,ALLOCATABLE::taual(:,:)  
     !labels
     CHARACTER(LEN=20), ALLOCATABLE :: label(:)
     CHARACTER(len=20), ALLOCATABLE :: speciesName(:)
     !name and other data of explicitely provided xc functional
     CHARACTER(len=4), ALLOCATABLE :: namex(:)
     INTEGER,          ALLOCATABLE :: icorr(:)
     INTEGER,          ALLOCATABLE :: igrd(:)
     INTEGER,          ALLOCATABLE :: krla(:)
     LOGICAL,          ALLOCATABLE :: relcor(:)
     !lda_u information(ntype)
     TYPE(t_utype),ALLOCATABLE::lda_u(:)
     INTEGER,ALLOCATABLE :: relax(:,:) !<(3,ntype)
     INTEGER, ALLOCATABLE :: nflip(:) !<flip magnetisation of this atom
     REAL,ALLOCATABLE:: vr0(:) !< Average Coulomb potential for atoms
  END TYPE t_atoms

  TYPE t_kpts
     INTEGER :: specificationType
     !no
     INTEGER :: nkpt
     INTEGER :: ntet
     REAL    :: posScale
     LOGICAL :: l_gamma
     !(3,nkpt) k-vectors internal units
     REAL,ALLOCATABLE ::bk(:,:)
     !(nkpts) weights
     REAL,ALLOCATABLE ::wtkpt(:)
     INTEGER               ::  nkptf !<k-vectors in full BZ
     INTEGER               ::  nkpt3(3)
     REAL   ,ALLOCATABLE   ::  bkf(:,:)
     INTEGER,ALLOCATABLE   ::  bkp(:)
     INTEGER,ALLOCATABLE   ::  bksym(:)
     INTEGER                       :: numSpecialPoints
     CHARACTER(LEN=50),ALLOCATABLE :: specialPointNames(:)
     REAL   ,ALLOCATABLE           :: specialPoints(:,:)
     INTEGER,ALLOCATABLE           :: ntetra(:,:)
     REAL   ,ALLOCATABLE           :: voltet(:)
  ENDTYPE t_kpts


  TYPE t_cell
     !name of 2D-lattice type
     CHARACTER*3::latnam
     !vol of dtilde box
     REAL::omtil
     !2D area
     REAL::area
     !bravais matrix
     REAL::amat(3,3)
     !rez. bravais matrx
     REAL::bmat(3,3)
     !square of bbmat
     REAL::bbmat(3,3)
     !d-value
     REAL::z1
     !volume of cell
     REAL::vol
     !volume of interstitial
     REAL::volint
     REAL:: c
  END TYPE t_cell

  !The stars
  TYPE t_stars
     !max-length of star
     REAL :: gmax
     REAL :: gmaxInit
     !no of 3d-stars
     INTEGER :: ng3
     !no of 2d-stars
     INTEGER :: ng2
     !dim of box
     INTEGER ::mx1
     INTEGER ::mx2
     INTEGER ::mx3
     !No of elements in FFT
     INTEGER ::kimax
     !No of elements in 2D-FFT
     INTEGER ::kimax2

     !Box for FFT in pwden
     INTEGER :: kq1_fft
     INTEGER :: kq2_fft
     INTEGER :: kq3_fft
     INTEGER :: kmxq_fft !no of g-vectors in sphere

     !fft box for xc-pot
     INTEGER :: kxc1_fft
     INTEGER :: kxc2_fft
     INTEGER :: kxc3_fft

     INTEGER :: ng3_fft
     INTEGER :: kmxxc_fft !<number of g-vectors forming the nxc3_fft stars in the charge density or xc-density sphere

     INTEGER :: nxc3_fft !< number of stars in the  charge density  fft-box

     !rep. g-vector of star
     INTEGER,ALLOCATABLE ::kv3(:,:)
     !length of star
     REAL,ALLOCATABLE    ::sk3(:)
     !mapping of g-vectors to stars
     INTEGER,ALLOCATABLE ::ig(:,:,:)
     !No of g-vectors in star
     INTEGER,ALLOCATABLE ::nstr(:)
     !rep. g-vector of 2D-star
     INTEGER,ALLOCATABLE ::kv2(:,:)
     !length of 2D-star
     REAL,ALLOCATABLE    ::sk2(:)
     !No of g-vecs in 2D-star
     INTEGER,ALLOCATABLE ::nstr2(:)
     !mapping of
     INTEGER,ALLOCATABLE ::ig2(:)
     !
     REAL,ALLOCATABLE:: phi2(:) !<(n2d)
     !phase phactor of g-vector
     COMPLEX,ALLOCATABLE    ::rgphs(:,:,:)
     !mapping of stars to FFT-box
     INTEGER, ALLOCATABLE :: igfft(:,:)
     !same for 2D
     INTEGER, ALLOCATABLE :: igfft2(:,:)
     !phasefactors for mapping
     COMPLEX,ALLOCATABLE  :: pgfft(:)
     !same of 2D
     COMPLEX,ALLOCATABLE  :: pgfft2(:)
     !
     REAL,ALLOCATABLE     :: ft2_gfx(:),ft2_gfy(:)
     COMPLEX, ALLOCATABLE :: ustep(:)
     REAL, ALLOCATABLE    :: ufft(:)
  END TYPE t_stars

  TYPE t_oneD
     TYPE (od_dim) :: odd
     TYPE (od_inp) :: odi
     TYPE (od_sym) :: ods
     TYPE (od_lda) :: odl
     TYPE (od_gga) :: odg
     INTEGER,  POINTER :: ig1(:,:)
     INTEGER,  POINTER :: kv1(:,:)
     INTEGER,  POINTER :: nstr1(:)
     INTEGER,  POINTER :: ngopr1(:)
     REAL,     POINTER :: mrot1(:,:,:)
     REAL,     POINTER :: tau1(:,:)
     INTEGER,  POINTER :: invtab1(:)
     INTEGER,  POINTER :: multab1(:,:)
     INTEGER,  POINTER :: igfft1(:,:)
     REAL,     POINTER :: pgfft1(:)
     REAL,     POINTER :: pgft1x(:)
     REAL,     POINTER :: pgft1y(:)
     REAL,     POINTER :: pgft1xx(:)
     REAL,     POINTER :: pgft1yy(:)
     REAL,     POINTER :: pgft1xy(:)
  END TYPE t_oneD

  TYPE t_hybrid
     LOGICAL               ::  l_hybrid
     LOGICAL               ::  l_subvxc
     LOGICAL               ::  l_calhf
     LOGICAL               ::  l_addhf
     INTEGER               ::  ewaldlambda
     INTEGER               ::  lexp
     INTEGER               ::  bands1 !Only read in
     INTEGER               ::  nbasp
     INTEGER               ::  maxlcutm1
     INTEGER               ::  maxindxm1
     INTEGER               ::  maxbasm1
     INTEGER               ::  maxindxp1
     INTEGER               ::  maxgptm
     INTEGER               ::  maxgptm1
     INTEGER               ::  maxindx
     INTEGER               ::  maxlmindx
     INTEGER               ::  gptmd
     INTEGER,ALLOCATABLE   ::  nindx(:,:)
     INTEGER,ALLOCATABLE   ::  select1(:,:)
     INTEGER,ALLOCATABLE   ::  lcutm1(:)
     INTEGER,ALLOCATABLE   ::  nindxm1(:,:)
     INTEGER,ALLOCATABLE   ::  gptm(:,:)
     INTEGER,ALLOCATABLE   ::  ngptm1(:)
     INTEGER,ALLOCATABLE   ::  pgptm1(:,:)
     INTEGER,ALLOCATABLE   ::  ngptm (:)
     INTEGER,ALLOCATABLE   ::  pgptm (:,:)
     INTEGER,ALLOCATABLE   ::  lcutwf(:)
     INTEGER,ALLOCATABLE   ::  map(:,:)
     INTEGER,ALLOCATABLE   ::  tvec(:,:,:)
     INTEGER , ALLOCATABLE ::  nbasm(:)                                     
     REAL                  ::  gcutm1
     REAL                  ::  tolerance1  !only read in
     REAL   ,ALLOCATABLE   ::  basm1(:,:,:,:)
     COMPLEX,ALLOCATABLE   ::  d_wgn2(:,:,:,:)
     INTEGER,ALLOCATABLE   ::  ne_eig(:),nbands(:),nobd(:)                   !alloc in eigen_HF_init
     REAL   ,ALLOCATABLE   ::  div_vv(:,:,:)
  END TYPE t_hybrid

  TYPE prodtype
     INTEGER :: l1,l2,n1,n2
  END TYPE prodtype
  TYPE t_hybdat
     INTEGER              :: lmaxcd,maxindxc
     REAL,  ALLOCATABLE   ::  gridf(:,:)                                    !alloc in util.F
     INTEGER , ALLOCATABLE::  nindxc(:,:)                                   !alloc in eigen_HF_init
     INTEGER,ALLOCATABLE  :: lmaxc(:)                                       !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  core1(:,:,:,:),core2(:,:,:,:)                 !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  eig_c(:,:,:)                                  !alloc in eigen_HF_init
     INTEGER , ALLOCATABLE::  kveclo_eig(:,:)                               !alloc in eigen_HF_setup
     INTEGER              ::  maxfac
     REAL,    ALLOCATABLE ::  sfac(:),fac(:)                                !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  gauntarr(:,:,:,:,:,:)                         !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  bas1(:,:,:,:),bas2(:,:,:,:)                   !alloc in eigen_HF_init
     REAL ,   ALLOCATABLE ::  bas1_MT(:,:,:),drbas1_MT(:,:,:)               !alloc in eigen_HF_init
     REAL, ALLOCATABLE    ::  prodm(:,:,:,:)                                !alloc in eigen_HF_setup
     TYPE(PRODTYPE),ALLOCATABLE :: prod(:,:,:)                              !alloc in eigen_HF_setup
     INTEGER, ALLOCATABLE :: pntgptd(:)                                     !alloc in eigen_HF_setup
     INTEGER, ALLOCATABLE :: pntgpt(:,:,:,:)                                !alloc in eigen_HF_setup
     INTEGER,ALLOCATABLE   ::  nindxp1(:,:)
  END TYPE t_hybdat

  TYPE t_dimension
     INTEGER :: jspd
     INTEGER :: nspd
     INTEGER :: nvd
     INTEGER :: nv2d
     INTEGER :: neigd
     INTEGER :: neigd2
     INTEGER :: ncvd
     INTEGER :: nn2d
     INTEGER :: nn3d
     INTEGER :: nstd
     INTEGER :: msh
     INTEGER :: lmd
     INTEGER :: lmplmd
     INTEGER :: nbasfcn
  END TYPE t_dimension

  TYPE t_Jij
     LOGICAL :: l_J
     INTEGER :: nqpt
     INTEGER :: nqptd
     INTEGER ::phnd
     INTEGER ::nsh
     INTEGER ::mtypes
     INTEGER :: nmopq(3)
     REAL    :: thetaJ
     REAL    :: qn
     INTEGER :: nmagn
     INTEGER :: nkpt_l
     LOGICAL :: l_disp
     LOGICAL :: l_wr
     LOGICAL :: l_jenerg
     REAL, ALLOCATABLE :: qj(:,:)
     LOGICAL, ALLOCATABLE:: l_magn(:)
     REAL, ALLOCATABLE   :: M(:)
     INTEGER, ALLOCATABLE:: magtype(:)
     INTEGER, ALLOCATABLE::nmagtype(:)
     REAL, ALLOCATABLE :: eig_l(:,:)
     REAL, ALLOCATABLE :: alph1(:)
  END TYPE t_Jij

  TYPE t_noco
     LOGICAL:: l_noco
     LOGICAL:: l_ss
     LOGICAL:: l_mperp
     LOGICAL:: l_constr
     REAL:: qss(3)
     REAL:: mix_b
     LOGICAL, ALLOCATABLE :: l_relax(:)
     REAL, ALLOCATABLE :: alphInit(:)
     REAL, ALLOCATABLE :: alph(:)
     REAL, ALLOCATABLE :: beta(:)
     REAL, ALLOCATABLE :: b_con(:,:)
     LOGICAL              :: l_soc
     LOGICAL, ALLOCATABLE :: soc_opt(:)
     REAL                 :: theta
     REAL                 :: phi
  END TYPE t_noco

  TYPE t_input
     LOGICAL :: strho
     LOGICAL :: cdinf
     LOGICAL :: vchk
     LOGICAL :: l_f
     LOGICAL :: eonly
     LOGICAL :: film
     LOGICAL :: ctail
     INTEGER :: coretail_lmax
     INTEGER :: itmax
     REAL    :: minDistance
     INTEGER :: maxiter
     INTEGER :: imix
     INTEGER :: gw
     INTEGER :: gw_neigd
     INTEGER :: qfix
     REAL    :: xa !< mixing parameter for geometry optimzer
     REAL    :: thetad !< Debey temperature for first step of geometry optimzer
     REAL    :: epsdisp !< minimal displacement. If all displacements are < epsdisp stop
     REAL    :: epsforce !< minimal force. If all forces <epsforce stop
     INTEGER :: isec1
     REAL    :: delgau
     REAL    :: alpha
     REAL    :: spinf
     REAL    :: tkb
     LOGICAL :: gauss
     LOGICAL :: l_bmt
     !INTEGER:: scale
     INTEGER:: jspins
     INTEGER:: kcrel
     LOGICAL:: frcor
     LOGICAL:: lflip
     LOGICAL:: score
     LOGICAL:: swsp
     LOGICAL:: tria
     LOGICAL:: integ
     LOGICAL:: pallst
     LOGICAL:: l_wann
     LOGICAL:: secvar
     LOGICAL:: evonly(2)
     LOGICAL:: eigvar(3)
     LOGICAL:: sso_opt(2)
     LOGICAL:: total
     LOGICAL:: l_inpXML
     REAL :: ellow
     REAL :: elup
     REAL :: rkmax
     REAL :: zelec
     CHARACTER(LEN=8) :: comment(10)
     TYPE(t_efield)::efield
     LOGICAL :: l_core_confpot
     LOGICAL :: l_useapw
  END TYPE t_input

  TYPE t_sliceplot
     LOGICAL :: iplot
     LOGICAL :: slice
     LOGICAL :: plpot
     INTEGER :: kk
     INTEGER :: nnne
     REAL    :: e1s
     REAL    :: e2s
  END TYPE t_sliceplot
  TYPE t_banddos
     LOGICAL :: dos
     LOGICAL :: band
     LOGICAL :: l_mcd
     LOGICAL :: l_orb
     LOGICAL :: vacdos
     INTEGER :: ndir
     INTEGER :: orbCompAtom
     REAL    :: e1_dos
     REAL    :: e2_dos
     REAL    :: sig_dos
  END TYPE t_banddos

  TYPE t_obsolete
     INTEGER:: lepr !floating energy parameters...
     LOGICAL:: disp
     INTEGER:: ndvgrd
     REAL   :: chng
     LOGICAL :: lwb
     LOGICAL:: l_u2f
     LOGICAL:: l_f2u
     LOGICAL :: pot8
  END TYPE t_obsolete


  TYPE t_enpara
     REAL, ALLOCATABLE :: el0(:,:,:)
     REAL, ALLOCATABLE :: evac0(:,:)
     REAL, ALLOCATABLE :: ello0(:,:,:)
     REAL, ALLOCATABLE :: enmix(:)
     INTEGER, ALLOCATABLE :: skiplo(:,:)
     LOGICAL, ALLOCATABLE :: lchange(:,:,:)
     LOGICAL, ALLOCATABLE :: lchg_v(:,:)
     LOGICAL, ALLOCATABLE :: llochg(:,:,:)
     REAL                 :: epara_min
  END TYPE t_enpara

  TYPE t_vacuum
     !Stuff for the vacuum
     INTEGER ::nmz
     INTEGER ::nmzd
     INTEGER ::nmzxy
     INTEGER ::nmzxyd
     INTEGER :: layerd
     INTEGER :: layers
     INTEGER :: nvac
     INTEGER :: nvacd
     REAL :: delz
     REAL :: dvac
     INTEGER::nstars
     INTEGER:: nstm
     REAL :: tworkf
     REAL :: locx(2)
     REAL :: locy(2)
     LOGICAL ::starcoeff
     INTEGER, ALLOCATABLE :: izlay(:,:)
  END TYPE t_vacuum


  !Data for the spherical harmonics
  TYPE t_sphhar
     !No of symmetry types (must
     !equal maxval(atoms%ntypsy)
     INTEGER ::ntypsd
     !Max no of members of sphhar
     INTEGER ::memd
     !max of nlh
     INTEGER ::nlhd
     !No of sphhar (ntypsd)
     INTEGER,ALLOCATABLE ::nlh(:)
     !l's of sphhar (0:nlhd,ntypsd)
     INTEGER,ALLOCATABLE ::llh(:,:)
     !No of members in sphhar (0:nlh
     INTEGER,ALLOCATABLE ::nmem(:,:)
     !lm's of of members (max(nmem),
     INTEGER,ALLOCATABLE ::mlh(:,:,:)
     !phasefactors (max(nmem),0:nlhd
     COMPLEX,ALLOCATABLE ::clnu(:,:,:)
  END TYPE t_sphhar

  !symmetry information
  TYPE t_sym
     INTEGER :: symSpecType
     !Symophic group
     LOGICAL ::symor
     INTEGER ::nsymt
     INTEGER :: nsym
     COMPLEX,ALLOCATABLE:: d_wgn(:,:,:,:)
     !2D-inv-sym
     LOGICAL ::invs2
     !Inversion-sym
     LOGICAL ::invs
     !Z-refls. sym
     LOGICAL ::zrfs
     !No of sym ops
     INTEGER ::nop
     !No of 2D-sym ops
     INTEGER ::nop2
     !Rot-matrices (3,3,nop)
     INTEGER,ALLOCATABLE::mrot(:,:,:)
     !inverse operation (nop)
     INTEGER,ALLOCATABLE::invtab(:)
     !translation vectors (3,nop)
     REAL,ALLOCATABLE::tau(:,:)
     !Name of lattice type
     CHARACTER*3   :: latnam
     !Name of sym
     CHARACTER*4   :: namgrp
     INTEGER, ALLOCATABLE :: multab(:,:)
     INTEGER, ALLOCATABLE :: invsatnr(:)
     INTEGER, ALLOCATABLE :: invarop(:,:)
     INTEGER, ALLOCATABLE :: invarind(:)

  END TYPE t_sym

  TYPE t_results
     REAL, ALLOCATABLE :: force(:,:,:)   !< Forces calculated on all atoms (for each spin)
     REAL, ALLOCATABLE :: force_old(:,:) !< Forces on all atoms from last iteration
     REAL              :: ef        !<Fermie energy
     REAL              :: seigc     !<sum of the core eigenvalues
     REAL              :: seigsc    !<weighted sum of the semi-core eigenvalues
     REAL              :: seigv     !<weighted sum of the occupied valence eigenvalues
     REAL              :: seigscv   !<sum of seigv and seigsc
     REAL              :: ts        !<entropy contribution to the free energy
     REAL              :: te_vcoul  !<charge density-coulomb potential integral
     REAL              :: te_veff   !<charge density-effective potential integral
     REAL              :: te_exc    !<charge density-ex-corr.energy density integral
     REAL              :: e_ldau    !<total energy contribution of LDA+U
     REAL              :: tote
     REAL              :: last_distance
     REAL              :: bandgap
     TYPE(t_energy_hf) ::  te_hfex
     REAL              ::  te_hfex_loc(2)
     REAL, ALLOCATABLE :: w_iks(:,:,:)
  END TYPE t_results


  TYPE t_mpi
     INTEGER :: mpi_comm !< replaces MPI_COMM_WORLD
     INTEGER :: irank    !< rank of task in mpi_comm
     INTEGER :: isize    !< no of tasks in mpi_comm
     INTEGER :: n_start  !< no of first k-point to calculate on this PE
     INTEGER :: n_stride !< stride for k-loops
     INTEGER :: n_size   !< PE per kpoint, i.e. "isize" for eigenvalue parallelization
     INTEGER :: n_groups !< No of k-loops per PE
     INTEGER :: sub_comm !< Sub-Communicator for eigenvalue parallelization (all PE working on same k-point)
     INTEGER :: n_rank   !< rank in sub_comm
  END TYPE t_mpi

  TYPE t_zMat
     LOGICAL              :: l_real
     INTEGER              :: nbasfcn
     INTEGER              :: nbands
     REAL,    ALLOCATABLE :: z_r(:,:) ! z_r(nbasfcn,nbands)
     COMPLEX, ALLOCATABLE :: z_c(:,:) ! z_c(nbasfcn,nbands)
  END TYPE t_zMat

  TYPE t_hamOvlp
     LOGICAL              :: l_real
     INTEGER              :: matsize
     REAL,    ALLOCATABLE :: a_r(:), b_r(:)
     COMPLEX, ALLOCATABLE :: a_c(:), b_c(:)
  END TYPE t_hamOvlp


  !
  ! type for wannier-functions
  !
  TYPE t_wann
     INTEGER :: wan90version
     INTEGER :: oc_num_orbs
     INTEGER,ALLOCATABLE :: oc_orbs(:)
     LOGICAL :: l_unformatted
     LOGICAL :: l_oc_f
     LOGICAL :: l_ndegen
     LOGICAL :: l_orbitalmom
     LOGICAL :: l_orbcomp
     LOGICAL :: l_orbcomprs
     LOGICAL :: l_denmat
     LOGICAL :: l_perturbrs
     LOGICAL :: l_perturb
     LOGICAL :: l_nedrho
     LOGICAL :: l_anglmomrs
     LOGICAL :: l_anglmom
     LOGICAL :: l_spindisp
     LOGICAL :: l_spindisprs
     LOGICAL :: l_socspicom
     LOGICAL :: l_socspicomrs
     LOGICAL :: l_offdiposoprs
     LOGICAL :: l_offdiposop
     LOGICAL :: l_torque
     LOGICAL :: l_torquers
     LOGICAL :: l_atomlist
     INTEGER :: atomlist_num
     INTEGER,ALLOCATABLE :: atomlist(:)
     LOGICAL :: l_berry
     LOGICAL :: l_perpmagrs
     LOGICAL :: l_perpmag
     LOGICAL :: l_perpmagat
     LOGICAL :: l_perpmagatrs
     LOGICAL :: l_socmatrs
     LOGICAL :: l_socmat
     LOGICAL :: l_soctomom
     LOGICAL :: l_kptsreduc2
     LOGICAL :: l_nablapaulirs
     LOGICAL :: l_nablars
     LOGICAL :: l_surfcurr
     LOGICAL :: l_updown
     LOGICAL :: l_ahe
     LOGICAL :: l_she
     LOGICAL :: l_rmat
     LOGICAL :: l_nabla
     LOGICAL :: l_socodi
     LOGICAL :: l_pauli
     LOGICAL :: l_pauliat
     LOGICAL :: l_potmat
     LOGICAL :: l_projgen
     LOGICAL :: l_plot_symm
     LOGICAL :: l_socmmn0
     LOGICAL :: l_bzsym
     LOGICAL :: l_hopping
     LOGICAL :: l_kptsreduc
     LOGICAL :: l_prepwan90
     LOGICAL :: l_plot_umdat
     LOGICAL :: l_wann_plot
     LOGICAL :: l_bynumber
     LOGICAL :: l_stopopt
     LOGICAL :: l_matrixmmn
     LOGICAL :: l_matrixamn
     LOGICAL :: l_projmethod
     LOGICAL :: l_wannierize
     LOGICAL :: l_plotw90
     LOGICAL :: l_byindex
     LOGICAL :: l_byenergy
     LOGICAL :: l_proj_plot
     LOGICAL :: l_bestproj
     LOGICAL :: l_ikptstart
     LOGICAL :: l_lapw
     LOGICAL :: l_plot_lapw
     LOGICAL :: l_fermi
     LOGICAL :: l_dipole
     LOGICAL :: l_dipole2
     LOGICAL :: l_dipole3
     LOGICAL :: l_mmn0
     LOGICAL :: l_mmn0at
     LOGICAL :: l_manyfiles
     LOGICAL :: l_collectmanyfiles
     LOGICAL :: l_ldauwan
     LOGICAL :: l_lapw_kpts
     LOGICAL :: l_lapw_gfleur
     LOGICAL :: l_kpointgen
     LOGICAL :: l_w90kpointgen
     LOGICAL :: l_finishnocoplot
     LOGICAL :: l_finishgwf
     LOGICAL :: l_skipkov
     LOGICAL :: l_matrixuHu
     LOGICAL :: l_matrixuHu_dmi
     INTEGER :: ikptstart
     INTEGER :: band_min(1:2)
     INTEGER :: band_max(1:2)
     INTEGER :: gfthick
     INTEGER :: gfcut
     INTEGER :: unigrid(6)
     INTEGER :: mhp(3)
     !---> gwf
     LOGICAL :: l_ms
     LOGICAL :: l_sgwf
     LOGICAL :: l_socgwf
     LOGICAL :: l_gwf
     LOGICAL :: l_bs_comf
     LOGICAL :: l_exist
     LOGICAL :: l_opened
     LOGICAL :: l_cleverskip
     LOGICAL :: l_dim(3)
     REAL    :: scale_param
     REAL    :: aux_latt_const
     REAL    :: hdwf_t1
     REAL    :: hdwf_t2
     INTEGER :: nparampts
     CHARACTER(len=20) :: fn_eig
     CHARACTER(len=20) :: param_file
     REAL,ALLOCATABLE :: param_vec(:,:)
     REAL,ALLOCATABLE :: param_alpha(:,:)
     CHARACTER(LEN=20), ALLOCATABLE :: jobList(:)
     !---> gwf

  END TYPE t_wann


  TYPE t_potden
     INTEGER             :: iter
     COMPLEX,ALLOCATABLE :: pw(:,:)
     REAL,ALLOCATABLE    :: mt(:,:,:,:)
     REAL,ALLOCATABLE    :: vacz(:,:,:)
     COMPLEX,ALLOCATABLE :: vacxy(:,:,:,:)
     !this type contains two init routines that should be used to allocate
     !memory. You can either specify the datatypes or give the dimensions as integers
     !See implementation below!
   CONTAINS
     PROCEDURE :: init_potden_types
     PROCEDURE :: init_potden_simple
     GENERIC   :: init=>init_potden_types,init_potden_simple
  END TYPE t_potden
CONTAINS
  SUBROUTINE usdus_init(ud,atoms,jsp)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_usdus)           :: ud
    TYPE(t_atoms),INTENT(IN) :: atoms
    INTEGER,INTENT(IN)       :: jsp

    INTEGER :: err(10)
    ALLOCATE ( ud%uloulopn(atoms%nlod,atoms%nlod,atoms%ntype,jsp),stat=err(1) )
    ALLOCATE ( ud%ddn(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(2) )
    ALLOCATE ( ud%us(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(3))
    ALLOCATE ( ud%uds(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(4) )
    ALLOCATE ( ud%dus(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(5))
    ALLOCATE ( ud%duds(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(6))
    ALLOCATE ( ud%ulos(atoms%nlod,atoms%ntype,jsp ),stat=err(7))
    ALLOCATE (ud%dulos(atoms%nlod,atoms%ntype,jsp ),stat=err(8) )
    ALLOCATE ( ud%uulon(atoms%nlod,atoms%ntype,jsp ),stat=err(9))
    ALLOCATE (ud%dulon(atoms%nlod,atoms%ntype,jsp) ,stat=err(10))

    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating usdus datatype")

  END SUBROUTINE usdus_init

  SUBROUTINE init_potden_types(pd,stars,atoms,sphhar,vacuum,oneD,jsp,l_noco)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT):: pd
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_sphhar),INTENT(IN):: sphhar
    TYPE(t_vacuum),INTENT(IN):: vacuum
    TYPE(t_oneD),INTENT(IN)  :: oneD
    INTEGER,INTENT(IN)       :: jsp
    LOGICAL,INTENT(IN)       :: l_noco
    CALL  init_potden_simple(pd,stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,jsp,l_noco,vacuum%nmzd,vacuum%nmzxyd,oneD%odi%n2d)
  END SUBROUTINE init_potden_types

  SUBROUTINE init_potden_simple(pd,ng3,jmtd,nlhd,ntype,jsp,l_noco,nmzd,nmzxyd,n2d)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT) :: pd
    INTEGER,INTENT(IN)          :: ng3,jmtd,nlhd,ntype,jsp
    LOGICAL,INTENT(IN)          :: l_noco
    INTEGER,INTENT(IN),OPTIONAL :: nmzd,nmzxyd,n2d

    INTEGER:: err(4)

    err=0
    pd%iter=0
    ALLOCATE(pd%pw(ng3,jsp),stat=err(1))
    ALLOCATE(pd%mt(jmtd,0:nlhd,ntype,jsp),stat=err(2))
    IF (PRESENT(nmzd)) THEN
       ALLOCATE(pd%vacz(nmzd,2,MERGE(jsp,4,l_noco)),stat=err(3))
       ALLOCATE(pd%vacxy(nmzxyd,n2d-1,2,jsp),stat=err(4))
    ENDIF
    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating potential or density")
    pd%pw=0.0
    pd%mt=0.0
    IF (PRESENT(nmzd)) THEN
       pd%vacz=0.0
       pd%vacxy=0.0
    ENDIF
  END SUBROUTINE init_potden_simple

 
END MODULE m_types
