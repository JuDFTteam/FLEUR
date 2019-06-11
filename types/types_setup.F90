!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_setup
   !*************************************************************
   !     This module contains definitions for all kind of types
   !*************************************************************
  USE m_types_cell
  USE m_types_sym
  USE m_types_banddos
  USE m_types_econfig
  USE m_types_input
  USE m_types_sliceplot
  USE m_types_oneD
  USE m_types_hybrid
  USE m_types_noco
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
      TYPE(t_econfig),ALLOCATABLE::econf(:)
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
      INTEGER, ALLOCATABLE :: nflip(:) !<flip magnetisation of this atom
   CONTAINS
      procedure :: nsp => calc_nsp_atom
   END TYPE t_atoms

  
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

  

   TYPE t_obsolete
      INTEGER:: lepr !floating energy parameters...
      INTEGER:: ndvgrd !remove
      REAL   :: chng   !remove
      LOGICAL :: lwb   !remove
   END TYPE t_obsolete

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
