!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

       MODULE m_types
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
   end type

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
   END TYPE

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
   END TYPE



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
        REAL u,j
        INTEGER l
        LOGICAL :: l_amf
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
        LOGICAL :: l_segmented = .false.
        LOGICAL :: plot_charge = .false. ! Plot charge as inputted
        LOGICAL :: plot_rho    = .false. ! Plot Fourier-transformed charge
        LOGICAL :: autocomp    = .true.  ! Auto-compensate film charge
        LOGICAL :: dirichlet = .false. ! Dirichlet vs. Neumann boundary cond.
        LOGICAL :: l_dirichlet_coeff = .false. ! For MPI, true if C1/C2 set
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
       !Claculate forces for this atom?
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
       !lda_u information(ntype)
       TYPE(t_utype),ALLOCATABLE::lda_u(:)
       INTEGER,ALLOCATABLE :: relax(:,:) !<(3,ntype)
       INTEGER, ALLOCATABLE :: nflip(:) !<flip magnetisation of this atom
       REAL,ALLOCATABLE:: vr0(:) !< Average Coulomb potential for atoms
      END TYPE

      TYPE t_kpts
       !no
       INTEGER :: nkpt
       INTEGER :: ntet
       REAL    :: posScale
       LOGICAL :: l_gamma
       INTEGER :: nmop(3) !<number of k-points in 3 directions
       !(3,nkpts) k-vectors internal units
       REAL,ALLOCATABLE ::bk(:,:)
       !(nkpts) weights
       REAL,ALLOCATABLE ::wtkpt(:)
       INTEGER, ALLOCATABLE :: pntgptd(:)
       INTEGER, ALLOCATABLE :: pntgpt(:,:,:,:,:)
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
      ENDTYPE


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
      END TYPE

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
        integer :: kq1_fft
        integer :: kq2_fft
        integer :: kq3_fft
        integer :: kmxq_fft !no of g-vectors in sphere

        !fft box for xc-pot
        integer :: kxc1_fft
        integer :: kxc2_fft
        integer :: kxc3_fft

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
     END TYPE

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
     END TYPE

     TYPE t_hybrid
      INTEGER               ::  ewaldlambda
      INTEGER               ::  lexp
      INTEGER               ::  bands1
      INTEGER               ::  bands2
      INTEGER(4)            ::  maxlcutm1
      INTEGER(4)            ::  maxindxm1
      INTEGER(4)            ::  maxbasm1
      INTEGER               ::  maxlcutm2
      INTEGER               ::  maxindxm2
      INTEGER               ::  maxbasm2
      INTEGER               ::  maxindxp1
      INTEGER               ::  maxindxp2
      INTEGER               ::  maxgptm
      INTEGER               ::  maxgptm1
      INTEGER               ::  maxgptm2
      INTEGER               ::  maxindx
      INTEGER               ::  maxlmindx
      INTEGER               ::  gptmd
      INTEGER,ALLOCATABLE   ::  nindx(:,:)
      INTEGER(4),ALLOCATABLE::  select1(:,:)
      INTEGER(4),ALLOCATABLE::  lcutm1(:)
      INTEGER(4),ALLOCATABLE::  select2(:,:)
      INTEGER,ALLOCATABLE   ::  lcutm2(:)
      INTEGER,ALLOCATABLE   ::  nindxm1(:,:)
      INTEGER,ALLOCATABLE   ::  nindxm2(:,:)
      INTEGER,ALLOCATABLE   ::  gptm(:,:)
      INTEGER,ALLOCATABLE   ::  ngptm1(:)
      INTEGER,ALLOCATABLE   ::  pgptm1(:,:)
      INTEGER,ALLOCATABLE   ::  ngptm2(:)
      INTEGER,ALLOCATABLE   ::  pgptm2(:,:)
      INTEGER,ALLOCATABLE   ::  ngptm (:)
      INTEGER,ALLOCATABLE   ::  pgptm (:,:)
      INTEGER,ALLOCATABLE   ::  nindxp1(:,:)
      INTEGER,ALLOCATABLE   ::  nindxp2(:,:)
      INTEGER,ALLOCATABLE   ::  lcutwf(:)
      INTEGER,ALLOCATABLE   ::  map(:,:)
      INTEGER,ALLOCATABLE   ::  tvec(:,:,:)
      REAL(8)               ::  radshmin
      REAL                  ::  gcutm1
      REAL                  ::  gcutm2
      REAL                  ::  tolerance1
      REAL                  ::  tolerance2
      REAL(8),ALLOCATABLE   ::  ddist(:)
      REAL   ,ALLOCATABLE   ::  basm1(:,:,:,:)
      REAL   ,ALLOCATABLE   ::  basm2(:,:,:,:)
      COMPLEX,ALLOCATABLE   ::  d_wgn2(:,:,:,:)
      LOGICAL               ::  l_subvxc
      LOGICAL               ::  l_calhf
      LOGICAL,ALLOCATABLE   ::  l_exxc(:,:)
     END TYPE

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
     END TYPE

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
     END TYPE

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
     END TYPE

     TYPE t_xcpot
        INTEGER :: icorr
        INTEGER :: igrd
        REAL    :: gmaxxc
     END TYPE


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
        INTEGER :: krla
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
     END TYPE

     type t_sliceplot
        logical :: iplot
        logical :: slice
        logical :: plpot
        integer :: kk
        integer :: nnne
        real    :: e1s
        real    :: e2s
     end type
     TYPE t_banddos
       LOGICAL :: dos
       LOGICAL :: band
       LOGICAL :: l_mcd
       LOGICAL :: l_orb
       LOGICAL :: vacdos
       INTEGER :: ndir
       REAL    :: e1_dos
       real    :: e2_dos
       real    :: sig_dos
     END TYPE

     TYPE t_obsolete
        INTEGER:: lepr !floating energy parameters...
        LOGICAL:: disp
        INTEGER:: ndvgrd
        REAL   :: chng
        LOGICAL :: lwb
        LOGICAL:: l_u2f
        LOGICAL:: l_f2u
        LOGICAL :: pot8
     END TYPE


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
     end type

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
      END TYPE


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
      END TYPE

      !symmetry information
      TYPE t_sym
       !Symophic group
       LOGICAL ::symor
       INTEGER ::nsymt
       INTEGER               ::  nsym

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

      END TYPE

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
      END TYPE


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
      END TYPE

      TYPE t_zMat
        LOGICAL              :: l_real
        INTEGER              :: nbasfcn
        INTEGER              :: nbands
        REAL,    ALLOCATABLE :: z_r(:,:) ! z_r(nbasfcn,nbands)
        COMPLEX, ALLOCATABLE :: z_c(:,:) ! z_c(nbasfcn,nbands)
      END TYPE

      TYPE t_hamOvlp
        LOGICAL              :: l_real
        INTEGER              :: matsize
        REAL,    ALLOCATABLE :: a_r(:), b_r(:)
        COMPLEX, ALLOCATABLE :: a_c(:), b_c(:)
      END TYPE

      END
