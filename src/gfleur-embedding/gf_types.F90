!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_types
   !*************************************************************
   !  Types of the Green's-function embedding code.
   !
   !  The standard FLEUR types used by the gfleur sources (t_atoms,
   !  t_cell, t_sym, t_stars, t_sphhar, t_enpara, t_lapw, t_mpi,
   !  t_potden, t_mat, ...) are re-exported from here; only genuinely
   !  GF-specific types are defined here. The historic local copies of FLEUR types were
   !  removed in the port to the modern FLEUR datatypes.
   !*************************************************************
   USE m_types_atoms
   USE m_types_cell
   USE m_types_enpara
   USE m_types_enparaXML
   USE m_types_fftGrid
   USE m_types_fleurinput
   USE m_types_forcetheo
   USE m_types_forcetheo_data
   USE m_types_hub1data
   USE m_types_hybdat
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_mat
   USE m_types_misc
   USE m_types_mpi
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_potden
   USE m_types_sphhar
   USE m_types_stars
   USE m_types_sym
   USE m_types_tlmplm
   USE m_types_usdus
   USE m_types_vacuum
   USE m_types_xcpot
   USE m_types_xcpot_inbuild
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_embinp, t_layers, t_gfmix, t_gfmpi, t_lapw_gf, t_stars_gf, t_gflayer, gf_l, gf_r
   !FLEUR types re-exported for the gfleur sources
   PUBLIC :: t_atoms, t_cell, t_enpara, t_enparaXML, t_fftGrid, t_fleurinput, t_forcetheo, t_forcetheo_data, t_hub1data, &
             t_hybdat, t_input, t_kpts, t_lapw, t_mat, t_results, t_mpi, t_noco, t_nococonv, &
             t_potden, t_sphhar, t_stars, t_sym, t_tlmplm, t_usdus, t_vacuum, t_xcpot, t_xcpot_inbuild

   !Two constants for the left and right side of a layer
   INTEGER, PARAMETER :: GF_L = 1
   INTEGER, PARAMETER :: GF_R = 2

   !Control parameters of the embedding calculation, read from the gf_inp
   !file. (Was called t_gfinp; renamed because modern FLEUR has an
   !unrelated t_gfinp for its onsite Green's functions.)
   TYPE t_embinp
      LOGICAL :: l_gf, l_eigen, l_surface, l_inp2plot, l_savemem, l_writeHS
      LOGICAL :: l_tmat, l_gmat, l_CBS, l_solwil, l_gproj, l_band, l_embmt
      LOGICAL :: l_charge, l_writeT, l_addemb, l_dos, l_potmix, l_adv, l_IEC
      LOGICAL :: l_totalmix, l_hdfio
      !use spectral representation to get Green function
      LOGICAL :: l_spectral, l_embspectral, l_simplevacuum
      !vacuum embedding below layer 1 (l_surface: above the last layer)
      LOGICAL :: l_vacleft = .FALSE.
      !calculate full Green function instead of only projections
      LOGICAL :: l_fullgreen
      LOGICAL :: l_nogno
      LOGICAL :: l_nohelpregion
      LOGICAL :: l_addselfen
      LOGICAL :: l_intdos
      CHARACTER(LEN=4) :: kpts
      !c1,c2, info about embedding plane
      REAL    :: dp1, dp2, CBS_bz
      REAL    :: eps_current, eps_non_bloch, charge_limit, kappa_max
      REAL    :: imag_broad   !broadening by imaginary energy
      !should an E-field be applied
      REAL    :: Efield, vacuum_energy
      !min and max of z-positions of dos
      REAL    :: z_min, z_max
      !relative threshold below which eigenvalues of S are dropped
      !in the spectral representation (0: plain generalized EVP)
      REAL    :: s_eps = 0.0
      !no of dos-planes
      INTEGER :: nz_dos
      INTEGER :: npw, curr
      INTEGER, ALLOCATABLE :: napw(:)
      LOGICAL, ALLOCATABLE :: l_doslayer(:)
      REAL    :: trans
      REAL    :: c_sc
      !cutoffs for the corrections of the Coulomb potential
      !(distributed into the per-layer t_stars_gf)
      REAL    :: gmax_pot, gmax_decouple
      !occupation/smearing parameters of the energy contour
      !(was the separate t_fermi type)
      LOGICAL :: l_gauss, l_tria
      REAL    :: tkb, delgau
   END TYPE t_embinp

   !Geometry of the layer decomposition along z
   TYPE t_layers
      INTEGER :: num_layers
      !prefix of the per-layer input files: <prefix><ilayer>_inp.xml
      CHARACTER(LEN=20) :: prefix = "layer"
      REAL, ALLOCATABLE :: c1(:), c2(:), d(:), dt(:), c(:)
   END TYPE t_layers

   !Mixing and self-consistency parameters of the embedding SCF
   !(was called t_mix)
   TYPE t_gfmix
      !no of iterations
      INTEGER :: iter
      !Mixing method
      INTEGER :: imix
      !Max no of iterations
      INTEGER :: maxiter
      !Mixing parameter
      REAL    :: alpha
      !Spin-mixing parameter
      REAL    :: spinf
      REAL    :: k_kerker, g0max, g0scale
      INTEGER :: precond
   END TYPE t_gfmix

   !MPI setup of gfleur: the standard FLEUR t_mpi plus the GF-specific
   !sub-communicators for the k-point x layer x energy decomposition
   !(set up in gf_mpi_groups).
   TYPE t_gfmpi
      TYPE(t_mpi) :: fmpi
      LOGICAL     :: pe0 = .FALSE.  !this PE does I/O (fmpi%irank==0)
      !k-parallelization
      INTEGER              :: k_kperPE, k_PEperK, k_rank
      INTEGER, ALLOCATABLE :: k_kpts(:)
      INTEGER              :: self_subcom
      !MPI-communicator for I/O
      INTEGER              :: iodop_subcom
      !communicators of all PEs working on the same layer/k
      INTEGER              :: samelayer_subcom, samek_subcom
      !k+layer parallelization
      INTEGER              :: kl_LayerPerPE
      INTEGER, ALLOCATABLE :: kl_layers(:)
      !k+energy parallelization
      INTEGER              :: ke_ENPerPE
      INTEGER, ALLOCATABLE :: ke_energies(:)
   END TYPE t_gfmpi

   !GF-specific extension of the LAPW basis (per layer and k-point).
   !The standard basis data lives in a modern t_lapw; this type carries
   !what the embedding needs on top: the cylinder-basis bookkeeping, the
   !2D (in-plane) basis with its mapping to a global list used for the
   !embedding-potential planes, and the sphere-restricted counts.
   TYPE t_lapw_gf
      !max length of lapw for this k / input value
      REAL    :: rkmax
      !matrix dimension and its restriction to the |k+G|<rkmax sphere
      INTEGER :: nmat, nmat_sphere
      !basis sizes (2*nv for noco)
      INTEGER :: nv_tot, nv_tot_sphere
      INTEGER :: nv2_tot
      !actual no of G's inside the rkmax sphere for this k (jspin)
      INTEGER :: nv_sphere(2)
      !actual no of 2D-G's for this k / the dimensions
      INTEGER :: nv2(2)
      INTEGER :: nvd, nv2d
      !cylinder basis (all G_z for the in-plane stars) in use?
      LOGICAL :: l_cylinder
      !3D basis index -> 2D basis index (per G, per spin)
      INTEGER, ALLOCATABLE :: kp(:, :)
      !the 2D g-vectors (sorted) and their lengths
      INTEGER, ALLOCATABLE :: k1p(:, :), k2p(:, :)
      REAL, ALLOCATABLE    :: rkp(:, :)
      !The 2D g-vectors of this k-point mapped to a global list
      INTEGER, ALLOCATABLE :: global2Dlist(:, :)
      INTEGER, ALLOCATABLE :: global2Dmap(:, :)
      INTEGER, ALLOCATABLE :: g2map(:)
      INTEGER              :: g_MAX(2)
   END TYPE t_lapw_gf

   !GF-specific extension of the star lists: additional cutoffs and the
   !"boxed" star sets filling the full FFT box that are used for the
   !embedding-potential planes and decoupled potentials, plus the
   !star<->FFT-grid maps the GF FFT code uses. Filled in gf_stars /
   !gf_stepsanaly.
   TYPE t_stars_gf
      !input cutoff
      REAL :: gmax_inp
      !cutoff for the decoupling correction of the Coulomb potential
      REAL :: gmax_decouple
      !cutoff for the correction of the Coulomb potential
      REAL :: gmax_pot
      !no of elements in the (2D-)FFT
      INTEGER :: kimax, kimax2
      !no of elements in z-direction
      INTEGER :: ngz, izmin, izmax
      INTEGER, ALLOCATABLE :: igz(:)
      !mapping of stars to the FFT-box (and phases)
      INTEGER, ALLOCATABLE :: igfft(:, :), igfft2(:, :)
      REAL, ALLOCATABLE    :: pgfft(:), pgfft2(:)
      REAL, ALLOCATABLE    :: pgft2x(:), pgft2y(:), pgft2xy(:)
      REAL, ALLOCATABLE    :: pgft2xx(:), pgft2yy(:)
   END TYPE t_stars_gf

   !Everything gfleur holds per layer: the complete FLEUR input of the
   !layer (from its own inp.xml) plus the derived per-layer objects.
   TYPE t_gflayer
      TYPE(t_fleurinput) :: fi
      TYPE(t_sphhar)     :: sphhar
      TYPE(t_stars)      :: stars
      TYPE(t_stars_gf)   :: stars_gf
      TYPE(t_enpara)     :: enpara
      TYPE(t_nococonv)   :: nococonv
      TYPE(t_hub1data)   :: hub1data
      CLASS(t_xcpot), ALLOCATABLE :: xcpot
      !per-layer snapshots of the global LAPW dimensions
      !(lapw_dim_nvd/nv2d/nbasfcn are module-global in m_types_lapw and
      !must be restored via lapw%init_dim before per-layer kernels run)
      INTEGER :: nvd = 0, nv2d = 0, nbasfcn = 0
      !the GF-specific LAPW bookkeeping of this layer (dimensions,
      !cylinder mode, global 2D lists; the per-k parts are refilled by
      !gf_apws inside the k-loops)
      TYPE(t_lapw_gf) :: lapw_gf
      !potential and density of the layer
      TYPE(t_potden) :: vTot
      TYPE(t_potden) :: cdn_new
      REAL, ALLOCATABLE :: qmtl_new(:, :)
      !radial functions and local Hamiltonian of the layer
      TYPE(t_tlmplm) :: tlmplm
   END TYPE t_gflayer

END MODULE m_gf_types
