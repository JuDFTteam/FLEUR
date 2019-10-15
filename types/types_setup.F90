!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_setup
   !*************************************************************
   !     This module contains definitions for all kind of types
   !*************************************************************

   ! types for 1D calculations
   TYPE od_dim
      LOGICAL :: d1
      INTEGER :: mb, M, k3, m_cyl
      INTEGER :: chi, rot
      LOGICAL :: invs, zrfs
      INTEGER :: n2d, nq2, nn2d
      INTEGER :: kimax2
      INTEGER :: nop, nat
   END TYPE od_dim

   TYPE od_inp
      LOGICAL :: d1
      INTEGER :: mb, M, k3, m_cyl
      INTEGER :: chi, rot
      LOGICAL :: invs, zrfs
      INTEGER :: n2d, nq2, nn2d
      INTEGER :: kimax2
      INTEGER, POINTER :: ig(:, :)  !(-k3:k3,-M:M)
      INTEGER, POINTER :: kv(:, :)        !(2,n2d)
      INTEGER, POINTER :: nst2(:)        !(n2d)
   END TYPE od_inp

   TYPE od_sym
      INTEGER :: nop, nat
      INTEGER, POINTER :: ngopr(:)     !(nat)
      REAL, POINTER :: mrot(:, :, :)  !(3,3,nop)
      REAL, POINTER :: tau(:, :)     !(3,nop)
      INTEGER, POINTER :: invtab(:)    !(nop)
      INTEGER, POINTER :: multab(:, :)  !(nop,nop)
   END TYPE od_sym

   TYPE od_lda
      INTEGER :: nn2d
      INTEGER, POINTER :: igf(:, :)  !(0:nn2d-1,2)
      REAL, POINTER :: pgf(:)    !(0:nn2d-1)
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
      REAL :: u, j         ! the actual U and J parameters
      REAL :: theta,phi   !the rotation angles by which the density metrics is rotated
      INTEGER :: l        ! the l quantum number to which this U parameter belongs
      INTEGER :: atomType ! The atom type to which this U parameter belongs
      LOGICAL :: l_amf ! logical switch to choose the "around mean field" LDA+U limit
   END TYPE t_utype

   !
   ! Type for the electric field
   !

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
      INTEGER, ALLOCATABLE ::nz(:)
      !atoms per type
      INTEGER, ALLOCATABLE::neq(:)
      !radial grid points
      INTEGER, ALLOCATABLE::jri(:)
      !core states
      INTEGER, ALLOCATABLE::ncst(:)
      !How many states are explicitely provided?
      INTEGER, ALLOCATABLE::numStatesProvided(:)
      !core state occupations
      REAL, ALLOCATABLE::coreStateOccs(:, :, :)
      !core state nprnc
      INTEGER, ALLOCATABLE::coreStateNprnc(:, :)
      !core state kappa
      INTEGER, ALLOCATABLE::coreStateKappa(:, :)
      !lmax
      INTEGER, ALLOCATABLE::lmax(:)
      !lmax non-spherical
      INTEGER, ALLOCATABLE::lnonsph(:)
      !expansion of pseudo-charge
      INTEGER, ALLOCATABLE::ncv(:)
      !no of LO
      INTEGER, ALLOCATABLE::nlo(:)
      !l of LO (nlo,ntype)
      INTEGER, ALLOCATABLE::llo(:, :)
      !lmax for lapw (ntype)
      INTEGER, ALLOCATABLE::lapw_l(:)
      !first LO with a given l (max(nlo
      INTEGER, ALLOCATABLE::lo1l(:, :)
      !??
      INTEGER, ALLOCATABLE::ulo_der(:, :)
      !no of LOs per l (max(nlo1),ntype
      INTEGER, ALLOCATABLE::nlol(:, :)
      !true if LO is formed by \dot u (
      LOGICAL, ALLOCATABLE::l_dulo(:, :)
      !no of op that maps atom into
      INTEGER, ALLOCATABLE::ngopr(:)
      !symetry of atom (nat)
      INTEGER, ALLOCATABLE::ntypsy(:)
      !no of sphhar for atom type(ntype
      INTEGER, ALLOCATABLE ::nlhtyp(:)
      !atom mapped to by inversion (nat
      INTEGER, ALLOCATABLE ::invsat(:)
      !Calaculate forces for this atom?
      LOGICAL, ALLOCATABLE :: l_geo(:)
      !MT-Radius (ntype)
      REAL, ALLOCATABLE CPP_MANAGED::rmt(:)
      !log increment(ntype)
      REAL, ALLOCATABLE::dx(:)
      !vol of MT(ntype)
      REAL, ALLOCATABLE::volmts(:)
      !radial grid points(max(jri),ntyp
      REAL, ALLOCATABLE::rmsh(:, :)
      !charge of nucleus(ntype)
      REAL, ALLOCATABLE::zatom(:)
      !initial mag moment(ntype)
      REAL, ALLOCATABLE::bmu(:)
      !pos of atom (absol) (3,nat)
      REAL, ALLOCATABLE::pos(:, :)
      !pos of atom (relat)(3,nat)
      REAL, ALLOCATABLE CPP_MANAGED::taual(:, :)
      !labels
      CHARACTER(LEN=20), ALLOCATABLE :: label(:)
      CHARACTER(len=20), ALLOCATABLE :: speciesName(:)
      !name and other data of explicitely provided xc functional
      CHARACTER(len=4), ALLOCATABLE :: namex(:)
      INTEGER, ALLOCATABLE :: icorr(:)
      INTEGER, ALLOCATABLE :: igrd(:)
      INTEGER, ALLOCATABLE :: krla(:)
      LOGICAL, ALLOCATABLE :: relcor(:)
      !lda_u information(ntype)
      TYPE(t_utype), ALLOCATABLE::lda_u(:)
      INTEGER, ALLOCATABLE :: relax(:, :) !<(3,ntype)
      REAL, ALLOCATABLE :: flipSpinPhi(:) !<flip magnetisation of this atom by angle phi
      REAL, ALLOCATABLE :: flipSpinTheta(:)
      LOGICAL :: l_flipSpinScale
   CONTAINS
      procedure :: nsp => calc_nsp_atom
   END TYPE t_atoms

   TYPE t_cell
      !name of 2D-lattice type
      CHARACTER*3::latnam
      !vol of dtilde box
      REAL::omtil
      !2D area
      REAL::area
      !bravais matrix
      REAL::amat(3, 3)
      !rez. bravais matrx
      REAL::bmat(3, 3)
      !square of bbmat
      REAL::bbmat(3, 3)
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
      INTEGER, ALLOCATABLE :: igq_fft(:)
      INTEGER, ALLOCATABLE :: igq2_fft(:)

      !fft box for xc-pot
      INTEGER :: kxc1_fft
      INTEGER :: kxc2_fft
      INTEGER :: kxc3_fft

      INTEGER :: ng3_fft
      INTEGER :: kmxxc_fft !<number of g-vectors forming the nxc3_fft stars in the charge density or xc-density sphere

      INTEGER :: nxc3_fft !< number of stars in the  charge density  fft-box

      !rep. g-vector of star
      INTEGER, ALLOCATABLE ::kv3(:, :)
      !length of star
      REAL, ALLOCATABLE    ::sk3(:)
      !mapping of g-vectors to stars
      INTEGER, ALLOCATABLE ::ig(:, :, :)
      !No of g-vectors in star
      INTEGER, ALLOCATABLE ::nstr(:)
      !rep. g-vector of 2D-star
      INTEGER, ALLOCATABLE ::kv2(:, :)
      !length of 2D-star
      REAL, ALLOCATABLE    ::sk2(:)
      !No of g-vecs in 2D-star
      INTEGER, ALLOCATABLE ::nstr2(:)
      !mapping of
      INTEGER, ALLOCATABLE ::ig2(:)
      !
      REAL, ALLOCATABLE:: phi2(:) !<(n2d)
      !phase phactor of g-vector
      COMPLEX, ALLOCATABLE    ::rgphs(:, :, :)
      !mapping of stars to FFT-box
      INTEGER, ALLOCATABLE :: igfft(:, :)
      !same for 2D
      INTEGER, ALLOCATABLE :: igfft2(:, :)
      !phasefactors for mapping
      COMPLEX, ALLOCATABLE  :: pgfft(:)
      !same of 2D
      COMPLEX, ALLOCATABLE  :: pgfft2(:)
      !
      REAL, ALLOCATABLE     :: ft2_gfx(:), ft2_gfy(:)
      COMPLEX, ALLOCATABLE :: ustep(:)
      REAL, ALLOCATABLE    :: ufft(:)
   END TYPE t_stars

   TYPE t_oneD
      TYPE(od_dim) :: odd
      TYPE(od_inp) :: odi
      TYPE(od_sym) :: ods
      TYPE(od_lda) :: odl
      TYPE(od_gga) :: odg
      INTEGER, POINTER :: ig1(:, :)
      INTEGER, POINTER :: kv1(:, :)
      INTEGER, POINTER :: nstr1(:)
      INTEGER, POINTER :: ngopr1(:)
      REAL, POINTER :: mrot1(:, :, :)
      REAL, POINTER :: tau1(:, :)
      INTEGER, POINTER :: invtab1(:)
      INTEGER, POINTER :: multab1(:, :)
      INTEGER, POINTER :: igfft1(:, :)
      REAL, POINTER :: pgfft1(:)
      REAL, POINTER :: pgft1x(:)
      REAL, POINTER :: pgft1y(:)
      REAL, POINTER :: pgft1xx(:)
      REAL, POINTER :: pgft1yy(:)
      REAL, POINTER :: pgft1xy(:)
   END TYPE t_oneD

   TYPE t_hybrid
      LOGICAL               ::  l_hybrid = .false.
      LOGICAL               ::  l_subvxc = .false.
      LOGICAL               ::  l_calhf = .false.
      LOGICAL               ::  l_addhf = .false.
      INTEGER               ::  ewaldlambda
      INTEGER               ::  lexp = 0
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
      INTEGER, ALLOCATABLE   ::  nindx(:, :)
      INTEGER, ALLOCATABLE   ::  select1(:, :)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  nindxm1(:, :)
      INTEGER, ALLOCATABLE   ::  gptm(:, :)
      INTEGER, ALLOCATABLE   ::  ngptm1(:)
      INTEGER, ALLOCATABLE   ::  pgptm1(:, :)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER, ALLOCATABLE   ::  pgptm(:, :)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:, :)
      INTEGER, ALLOCATABLE   ::  tvec(:, :, :)
      INTEGER, ALLOCATABLE ::  nbasm(:)
      REAL                  ::  gcutm1
      REAL                  ::  tolerance1  !only read in
      REAL, ALLOCATABLE   ::  basm1(:, :, :, :)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:, :, :, :)
      INTEGER, ALLOCATABLE   ::  ne_eig(:), nbands(:), nobd(:)                   !alloc in eigen_HF_init
      REAL, ALLOCATABLE   ::  div_vv(:, :, :)
   END TYPE t_hybrid

   TYPE t_dimension
      INTEGER :: nspd
      INTEGER :: nvd
      INTEGER :: nv2d
      INTEGER :: neigd
      INTEGER :: neigd2
      INTEGER :: ncvd
      INTEGER :: nstd
      INTEGER :: msh
      INTEGER :: lmd
      INTEGER :: lmplmd
      INTEGER :: nbasfcn
   END TYPE t_dimension

   TYPE t_noco
      LOGICAL:: l_noco
      LOGICAL:: l_ss
      LOGICAL:: l_mperp
      LOGICAL:: l_constr
      LOGICAL:: l_mtNocoPot
      REAL:: qss(3)
      REAL:: mix_b
      LOGICAL, ALLOCATABLE :: l_relax(:)
      REAL, ALLOCATABLE :: alphInit(:)
      REAL, ALLOCATABLE :: alph(:)
      REAL, ALLOCATABLE :: beta(:)
      REAL, ALLOCATABLE :: b_con(:, :)
      LOGICAL           :: l_soc
      LOGICAL           :: l_spav
      REAL              :: theta
      REAL              :: phi
      REAL, ALLOCATABLE  :: socscale(:)
   END TYPE t_noco

   TYPE t_input
      LOGICAL :: eig66(2)
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
      INTEGER :: numBandsKPoints
      REAL    :: forcealpha !< mixing parameter for geometry optimzer
      REAL    :: epsdisp !< minimal displacement. If all displacements are < epsdisp stop
      REAL    :: epsforce !< minimal force. If all forces <epsforce stop
      REAL    :: force_converged=0.00001
      INTEGER :: forcemix=3
      REAL    :: delgau
      REAL    :: alpha
      REAL    :: preconditioning_param
      REAL    :: spinf
      REAL    :: tkb
      LOGICAL :: gauss
      LOGICAL :: l_bmt
      !INTEGER:: scale
      INTEGER:: jspins
      INTEGER:: kcrel
      LOGICAL:: frcor
      LOGICAL:: lflip
      LOGICAL:: swsp
      LOGICAL:: tria
      LOGICAL:: integ
      LOGICAL:: pallst
      LOGICAL:: l_coreSpec
      LOGICAL:: l_wann
      LOGICAL:: secvar
      LOGICAL:: evonly
      LOGICAL:: total
      LOGICAL:: l_inpXML
      REAL :: scaleCell
      REAL :: scaleA1
      REAL :: scaleA2
      REAL :: scaleC
      REAL :: ellow
      REAL :: elup
      REAL :: rkmax
      REAL :: zelec
      REAL :: fixed_moment = 0.0
      CHARACTER(LEN=8) :: comment(10)
      REAL, POINTER :: sigma !this is the difference in charge due to the electric field it points to the value stored in t_efield
      LOGICAL :: l_core_confpot
      LOGICAL :: l_useapw
      LOGICAL :: ldauLinMix
      REAL    :: ldauMixParam
      REAL    :: ldauSpinf
      LOGICAL :: l_rdmft
      REAL    :: rdmftOccEps
      INTEGER :: rdmftStatesBelow
      INTEGER :: rdmftStatesAbove
      INTEGER :: rdmftFunctional
   END TYPE t_input

   TYPE t_sliceplot
      INTEGER :: iplot
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
      REAL    :: e_mcd_lo
      REAL    :: e_mcd_up
      LOGICAL :: unfoldband
      INTEGER :: s_cell_x
      INTEGER :: s_cell_y
      INTEGER :: s_cell_z
      REAL    :: alpha,beta,gamma !For orbital decomp. (was orbcomprot)
   END TYPE t_banddos

   TYPE t_obsolete
      INTEGER:: lepr !floating energy parameters...
      INTEGER:: ndvgrd !remove
      REAL   :: chng   !remove
      LOGICAL :: lwb   !remove
   END TYPE t_obsolete

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
      INTEGER, ALLOCATABLE :: izlay(:, :)
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
      INTEGER, ALLOCATABLE ::nlh(:)
      !l's of sphhar (0:nlhd,ntypsd)
      INTEGER, ALLOCATABLE ::llh(:, :)
      !No of members in sphhar (0:nlh
      INTEGER, ALLOCATABLE ::nmem(:, :)
      !lm's of of members (max(nmem),
      INTEGER, ALLOCATABLE ::mlh(:, :, :)
      !phasefactors (max(nmem),0:nlhd
      COMPLEX, ALLOCATABLE ::clnu(:, :, :)
   END TYPE t_sphhar

   !symmetry information
   TYPE t_sym
      INTEGER :: symSpecType
      !Symophic group
      LOGICAL ::symor
      INTEGER ::nsymt
      INTEGER :: nsym
      COMPLEX, ALLOCATABLE:: d_wgn(:, :, :, :)
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
      INTEGER, ALLOCATABLE::mrot(:, :, :)
      !inverse operation (nop)
      INTEGER, ALLOCATABLE::invtab(:)
      !translation vectors (3,nop)
      REAL, ALLOCATABLE::tau(:, :)
      !Name of lattice type
      CHARACTER*3   :: latnam
      !Name of sym
      CHARACTER*4   :: namgrp
      INTEGER, ALLOCATABLE :: multab(:, :)
      INTEGER, ALLOCATABLE :: invsatnr(:)
      INTEGER, ALLOCATABLE :: invarop(:, :)
      INTEGER, ALLOCATABLE :: invarind(:)

   END TYPE t_sym

   ! type for the input to the calculation of the core spectrum (EELS)
   TYPE t_coreSpecInput
      integer :: verb  ! output verbosity
      integer :: atomType  ! atomic type used for calculation of core spectra
      character(LEN=1) :: edge  ! edge character (K,L,M,N,O,P)
      integer :: edgeidx(11)  ! l-j edges
      integer :: lx  ! maximum lmax considered in spectra calculation
      real :: ek0  ! kinetic energy of incoming electrons
      real :: emn  ! energy spectrum lower bound
      real :: emx  ! energy spectrum upper bound
      real :: ein  ! energy spectrum increment
      integer :: nqphi ! no. of angle-sectors for integral over q vectors
      integer :: nqr   ! no. of radial-sectors for integral over q vectors
      real :: alpha_ex  ! maximal angle of incoming electrons
      real :: beta_ex   ! maximal (measured) angle of outcoming electrons
      real :: I0        ! incoming intensity
   END TYPE t_coreSpecInput

   !
   ! type for wannier-functions
   !
   TYPE t_wann
      INTEGER :: wan90version
      INTEGER :: oc_num_orbs
      INTEGER, ALLOCATABLE :: oc_orbs(:)
      LOGICAL :: l_mmn0_unf_to_spn_unf
      LOGICAL :: l_mmn0_to_spn_unf
      LOGICAL :: l_mmn0_to_spn
      LOGICAL :: l_mmn0_to_spn2
      LOGICAL :: l_mmn0_unf_to_spn
      LOGICAL :: l_perpmag_unf_to_tor_unf 
      LOGICAL :: l_perpmag_to_tor_unf
      LOGICAL :: l_perpmag_to_tor
      LOGICAL :: l_perpmag_unf_to_tor
      LOGICAL :: l_hsomtxvec_unf_to_lmpzsoc_unf
      LOGICAL :: l_hsomtxvec_to_lmpzsoc_unf
      LOGICAL :: l_hsomtxvec_to_lmpzsoc
      LOGICAL :: l_hsomtxvec_unf_to_lmpzsoc
      LOGICAL :: l_hsomtx_unf_to_hsoc_unf
      LOGICAL :: l_hsomtx_to_hsoc_unf
      LOGICAL :: l_hsomtx_to_hsoc
      LOGICAL :: l_hsomtx_unf_to_hsoc
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
      INTEGER, ALLOCATABLE :: atomlist(:)
      LOGICAL :: l_berry
      LOGICAL :: l_perpmagrs
      LOGICAL :: l_perpmag
      LOGICAL :: l_perpmagat
      INTEGER :: perpmagl
      LOGICAL :: l_perpmagatlres
      LOGICAL :: l_perpmagatrs
      LOGICAL :: l_socmatrs
      LOGICAL :: l_socmat
      LOGICAL :: l_socmatvec
      LOGICAL :: l_socmatvecrs
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
      REAL, ALLOCATABLE :: param_vec(:, :)
      REAL, ALLOCATABLE :: param_alpha(:, :)
      CHARACTER(LEN=20), ALLOCATABLE :: jobList(:)
      !---> gwf

   END TYPE t_wann
CONTAINS
   pure function calc_nsp_atom(self) result(nsp) 
      implicit none
      CLASS(t_atoms),INTENT(IN)      :: self
      INTEGER                        :: nsp

      nsp = (self%lmaxd+1+MOD(self%lmaxd+1,2))*(2*self%lmaxd+1)
   end function
END MODULE m_types_setup
