!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_driver
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_UNDERFLOW, IEEE_GET_FLAG, IEEE_SET_FLAG
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_eig66_io, ONLY: read_eig
   USE m_genMTBasis, ONLY: genMTBasis
   USE m_juDFT, ONLY: juDFT_error
   USE m_mpi_reduce_tool, ONLY: mpi_sum_reduce
   USE m_types_abc, ONLY: t_abc
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_cell, ONLY: t_cell
   USE m_types_enpara, ONLY: t_enpara
   USE m_types_input, ONLY: t_input
   USE m_types_kpts, ONLY: t_kpts
   USE m_types_lapw, ONLY: t_lapw
   USE m_types_mat, ONLY: t_mat
   USE m_types_misc, ONLY: t_results
   USE m_types_mpi, ONLY: t_mpi
   USE m_types_noco, ONLY: t_noco
   USE m_types_nococonv, ONLY: t_nococonv
   USE m_types_potden, ONLY: t_potden
   USE m_types_radfun, ONLY: t_radfun
   USE m_types_sym, ONLY: t_sym
   USE m_types_usdus, ONLY: t_usdus
   USE m_xas_angular, ONLY: xas_cartesian_to_spherical, xas_print_angular_sumrule
   USE m_xas_core, ONLY: t_xas_core_state, xas_extract_core_states
   USE m_xas_io, ONLY: xas_write_spectrum_text
   USE m_xas_matrixelements, ONLY: xas_core_band_matrixelements
   USE m_xas_radial, ONLY: xas_radial_dipole_integrals
   USE m_xas_spectrum, ONLY: xas_accumulate_matrix_spectrum
	   USE m_xas_symmetry, ONLY: xas_count_star_members, xas_star_member_weight, xas_rotate_abc_star_member, &
	                             xas_rotate_abc_star_member_spinor, xas_su2_from_sym, xas_local_spin_transform, &
	                             xas_cart_rotation_from_sym, xas_star_operation, xas_print_symmetry_rotation_diagnostics
   IMPLICIT NONE
   PRIVATE

   ! Temporary hardwired XAS development parameters. Fe L3 remains the
   ! default validation case, but these choices are no longer baked into the
   ! driver logic and should become XML/input controlled later.
   INTEGER, PARAMETER :: xas_absorber_z = 26
   CHARACTER(LEN=2), PARAMETER :: xas_edge_name = "L3"
   INTEGER, PARAMETER :: xas_test_n_grid = 401
   ! Hardwired Gaussian broadening for the smoke test, in Hartree.
   ! Useful trial values: 0.01 Ha (sharper), 0.03 Ha (default), 0.05 Ha (smoother).
   REAL,    PARAMETER :: xas_test_eta = 0.03
   INTEGER, PARAMETER :: xas_debug_n_pol = 3
   ! Debug verbosity:
   !   0: only final spectrum filename
   !   1: compact atom/spin/total/l-channel strength summary
   !   2: level 1 plus compact k-point contribution table
   !   3: level 2 plus floating-point underflow diagnostics
   INTEGER, PARAMETER :: xas_debug_verbosity = 1
   ! Default star treatment: reuse parent irreducible-k eigenvectors and rotate
   ! MT abc coefficients with Wigner-D matrices for each full-zone star member.
   ! The lab-frame polarization is kept fixed. Older diagnostic paths that only
   ! rotated polarization, or used full-zone read_eig, were not production-valid.
   LOGICAL, PARAMETER :: xas_use_spatial_star = .TRUE.
   ! Development diagnostic for the SOC/noncollinear spinor part of the star
   ! transform. This prints SU(2) and local-spin-frame matrices only; it is not
   ! connected to spectrum accumulation yet.
   LOGICAL, PARAMETER :: xas_debug_spinor_star = .FALSE.
   ! Temporary pure-angular L3 sum-rule diagnostic. Enable manually when
   ! checking angular isotropy; it does not affect spectra.
   LOGICAL, PARAMETER :: xas_debug_angular_sumrule = .FALSE.
   ! Temporary symmetry-basis diagnostic. Enable manually to print the raw
	   ! lattice/fractional operations and the Cartesian/proper rotations used by
	   ! XAS Wigner-D and SU(2) star reconstruction.
	   LOGICAL, PARAMETER :: xas_debug_symmetry_rotations = .FALSE.
   ! Temporary abc fingerprint diagnostic for spinor spatial-star validation.
   ! Explicit full-k runs write DIRECT records; symmetry-reduced runs write
   ! STAR records keyed by the same full-zone k coordinates for comparison.
   LOGICAL, PARAMETER :: xas_debug_abc_star_compare = .FALSE.
   INTEGER, PARAMETER :: xas_debug_abc_star_max_records = 20
	   CHARACTER(LEN=1), PARAMETER :: xas_debug_pol_label(xas_debug_n_pol) = [CHARACTER(LEN=1) :: "x", "y", "z"]

   PUBLIC :: xas_hardwired_test_driver

CONTAINS

   SUBROUTINE xas_hardwired_test_driver(eig_id, fmpi, input, kpts, atoms, sym, cell, noco, nococonv, enpara, vTot, results)
      ! Temporary internal XAS smoke-test driver.
      !
      ! Hardwired choices:
      !   all atom types with Z=xas_absorber_z, edge=xas_edge_name, z
      !   polarization, Gaussian eta=xas_test_eta Ha, output file labelled by
      !   edge, polarization, and broadening. The absorbing-atom selection is
      !   temporary hardwired behavior and must become XML/input controlled in
      !   the real driver.
      !
      ! The current test spectrum is a per-selected-absorbing-atom local
      ! muffin-tin XAS signal. It is k-weighted and written in arbitrary units,
      ! but it is not normalized per cell volume, film area, film thickness, or
      ! number of equivalent atoms. Additive spectral and diagnostic quantities
      ! are accumulated over the rank-local fmpi%k_list and reduced to rank 0
      ! before output.
      !
      ! The local dipole approximation intentionally neglects interstitial and
      ! vacuum contributions. This is appropriate for core-level local XAS as a
      ! first implementation.
      !
      ! Spatial-star handling reuses parent irreducible-k eigenvectors, rotates
      ! the MT angular coefficients with Wigner-D matrices, and keeps the
      ! laboratory polarization fixed for each full k-star member. The current
      ! atom-mapping helper requires mapped absorbers to stay in the same atom
      ! type. Scalar calculations use orbital Wigner-D only. SOC/noncollinear
      ! calculations use orbital Wigner-D plus local-frame SU(2) for unitary
      ! spatial operations; time-reversal star handling is still guarded.
      !
      ! For k-mesh and Kmax convergence tests, reuse the same converged cdn.hdf
      ! (for example from a dense 12x12x12 SCF run), then vary only the XAS
      ! evaluation k mesh and basis settings in one-shot runs. This separates
      ! spectral k/basis convergence from self-consistent density changes.
      INTEGER,             INTENT(IN) :: eig_id
      TYPE(t_mpi),         INTENT(IN) :: fmpi
      TYPE(t_input),       INTENT(IN) :: input
      TYPE(t_kpts),        INTENT(IN) :: kpts
      TYPE(t_atoms),       INTENT(IN) :: atoms
      TYPE(t_sym),         INTENT(IN) :: sym
      TYPE(t_cell),        INTENT(IN) :: cell
      TYPE(t_noco),        INTENT(IN) :: noco
      TYPE(t_nococonv),    INTENT(IN) :: nococonv
      TYPE(t_enpara),      INTENT(IN) :: enpara
      TYPE(t_potden),      INTENT(IN) :: vTot
      TYPE(t_results),     INTENT(IN) :: results

      TYPE(t_usdus) :: usdus
      TYPE(t_radfun) :: radfun
      TYPE(t_xas_core_state), ALLOCATABLE :: core_states(:)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      TYPE(t_abc), ALLOCATABLE :: abc_spin(:)
      TYPE(t_abc), ALLOCATABLE :: abc_star_spin(:)

      COMPLEX :: eps_cart(3), eps_sph(-1:1), eps_cart_debug(3), eps_sph_debug(-1:1)
      COMPLEX :: spin_frame_transform(2, 2)
      COMPLEX, ALLOCATABLE :: matrix(:, :), matrix_debug(:, :), matrix_l0(:, :), matrix_l2(:, :)
      REAL, ALLOCATABLE :: energy_grid(:), intensity(:), radial_xas(:, :, :)
      REAL, ALLOCATABLE :: intensity_reduced(:)
      REAL, ALLOCATABLE :: f(:, :, :, :), g(:, :, :, :), flo(:, :, :, :)
      REAL, ALLOCATABLE :: eig_band(:), occ_band(:)
      REAL, ALLOCATABLE :: xas_debug_strength_kpt(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_kpt_reduced(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_spin(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_spin_reduced(:, :)
      INTEGER, ALLOCATABLE :: ev_list(:)

      CHARACTER(LEN=1)   :: pol_label
      CHARACTER(LEN=5)   :: eta_label
      CHARACTER(LEN=128) :: output_filename
      CHARACTER(LEN=128) :: xas_debug_filename
      CHARACTER(LEN=200) :: error_message
      INTEGER :: ikpt_i, ikpt, ikptf, bksym, jsp_loop, jsp, ispin, nbands, nbands_read, itype
      INTEGER :: n_spin_channels, n_local_spins, max_order, nbasfcn, lmax_xas, i_band
	      INTEGER :: n_underflow_spectrum, i_char, i_pol, iatom_l, n_absorber_types, n_absorber_atoms, nstar
	      INTEGER :: n_abc_debug_direct_records, n_abc_debug_star_records
	      INTEGER :: xas_debug_unit, xas_abc_debug_unit
      REAL    :: xas_debug_strength_l0(xas_debug_n_pol), xas_debug_strength_l2(xas_debug_n_pol)
      REAL    :: xas_debug_strength_total(xas_debug_n_pol)
      REAL    :: xas_debug_strength_l0_reduced(xas_debug_n_pol), xas_debug_strength_l2_reduced(xas_debug_n_pol)
      REAL    :: xas_debug_strength_total_reduced(xas_debug_n_pol), xas_debug_strength_cross_reduced(xas_debug_n_pol)
      REAL    :: weight_sums_local(2), weight_sums_reduced(2)
      REAL    :: transition_min, transition_max, transition_padding, occ, debug_strength, debug_l0, debug_l2
      REAL    :: debug_avg, debug_rel_xz
      REAL    :: wk_current, wk_star, weight_sum_parent, weight_sum_star
      INTEGER :: underflow_local(1), underflow_reduced(1)
      LOGICAL :: l_real, l_xas_debug_fp, l_xas_debug_strength, l_xas_debug_kpt_strength, l_root
      LOGICAL :: l_spinor_abc, l_xas_angular_sumrule_printed

      l_root = fmpi%irank == 0
      l_xas_debug_fp = xas_debug_verbosity >= 3
      l_xas_debug_strength = xas_debug_verbosity >= 1
      l_xas_debug_kpt_strength = xas_debug_verbosity >= 2
      n_underflow_spectrum = 0
      xas_debug_strength_l0 = 0.0
      xas_debug_strength_l2 = 0.0
      xas_debug_strength_total = 0.0
      n_absorber_types = 0
      n_absorber_atoms = 0
	      weight_sum_parent = 0.0
	      weight_sum_star = 0.0
	      xas_debug_unit = -1
	      xas_abc_debug_unit = -1
	      n_abc_debug_direct_records = 0
	      n_abc_debug_star_records = 0
	      l_xas_angular_sumrule_printed = .FALSE.

      IF (.NOT. ALLOCATED(results%w_iks)) THEN
         CALL juDFT_error("results%w_iks is not allocated in xas_hardwired_test_driver", calledby="m_xas_driver")
      END IF

	      xas_debug_filename = ""
	      IF (l_root) CALL xas_debug_open_log(kpts, xas_use_spatial_star, xas_debug_unit, xas_debug_filename)
	      IF (l_root .AND. xas_debug_abc_star_compare) THEN
	         OPEN(NEWUNIT=xas_abc_debug_unit, FILE="xas_abc_star_compare.dat", STATUS="REPLACE", ACTION="WRITE")
	         CALL xas_debug_write_abc_header(xas_abc_debug_unit)
	      END IF

      IF (l_root .AND. xas_debug_verbosity >= 3) THEN
         WRITE(*, '(a,i0)') "XAS DEBUG atom types: ntype=", atoms%ntype
         WRITE(xas_debug_unit, '(a,i0)') "XAS DEBUG atom types: ntype=", atoms%ntype
      END IF
      DO itype = 1, atoms%ntype
         IF (l_root .AND. xas_debug_verbosity >= 3) THEN
            WRITE(*, '(a,i0,a,i0,a,i0,a,a)') "XAS DEBUG atom type ", itype, " Z=", atoms%nz(itype), &
                                              " neq=", atoms%neq(itype), " species=", TRIM(atoms%speciesName(itype))
            WRITE(xas_debug_unit, '(a,i0,a,i0,a,i0,a,a)') "XAS DEBUG atom type ", itype, " Z=", atoms%nz(itype), &
                                                          " neq=", atoms%neq(itype), " species=", TRIM(atoms%speciesName(itype))
         END IF
         IF (atoms%nz(itype) == xas_absorber_z) THEN
            n_absorber_types = n_absorber_types + 1
            n_absorber_atoms = n_absorber_atoms + atoms%neq(itype)
         END IF
      END DO
      IF (n_absorber_types == 0) THEN
         WRITE(error_message, '(a,i0)') "No atom types found for requested XAS absorber Z=", xas_absorber_z
         CALL juDFT_error(TRIM(error_message), calledby="m_xas_driver")
      END IF
      IF (l_root .AND. xas_debug_spinor_star) THEN
         CALL xas_debug_print_spinor_star(sym, cell, atoms, nococonv, xas_absorber_z, xas_debug_unit)
      END IF
      IF (l_root .AND. xas_debug_symmetry_rotations) THEN
         CALL xas_print_symmetry_rotation_diagnostics(sym, cell, xas_debug_unit)
      END IF

      CALL usdus%init(atoms, input%jspins)
      ALLOCATE(f(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
      ALLOCATE(g(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
      ALLOCATE(flo(atoms%jmtd, 2, atoms%nlod, input%jspins))
      IF (l_xas_debug_kpt_strength) THEN
         ALLOCATE(xas_debug_strength_kpt(xas_debug_n_pol, kpts%nkpt), SOURCE=0.0)
      END IF

      eps_cart = CMPLX(0.0, 0.0)
      ! Hardwired linear polarization vector in Cartesian components:
      !   x: eps_cart(1) = CMPLX(1.0, 0.0); pol_label = "x"
      !   y: eps_cart(2) = CMPLX(1.0, 0.0); pol_label = "y"
      !   z: eps_cart(3) = CMPLX(1.0, 0.0); pol_label = "z"  ! default
      pol_label = "z"
      eps_cart(3) = CMPLX(1.0, 0.0)
      CALL xas_cartesian_to_spherical(eps_cart, eps_sph)

      WRITE(eta_label, '(f5.3)') xas_test_eta
      DO i_char = 1, LEN(eta_label)
         IF (eta_label(i_char:i_char) == ".") eta_label(i_char:i_char) = "p"
      END DO
      output_filename = "xas_test_"//TRIM(xas_edge_name)//"_"//TRIM(pol_label)//"_eta"//TRIM(eta_label)//".dat"

      ! Keep FLEUR eig/occupation spin channels separate from the local MT
      ! spinor components used by the XAS dipole matrix element. In noco
      ! calculations the eig file is addressed with jsp=1, while abc%calc_abc
      ! must still be called for the two local components ispin=1,2.
      l_spinor_abc = noco%l_noco
      n_spin_channels = MERGE(1, input%jspins, l_spinor_abc)
      IF (l_xas_debug_strength) THEN
         ALLOCATE(xas_debug_strength_spin(xas_debug_n_pol, n_spin_channels), SOURCE=0.0)
      END IF
      DO ikpt_i = 1, SIZE(fmpi%k_list)
         ikpt = fmpi%k_list(ikpt_i)
         IF (kpts%wtkpt(ikpt) <= 0.0) CYCLE
         weight_sum_parent = weight_sum_parent + kpts%wtkpt(ikpt)
         IF (xas_use_spatial_star) THEN
            CALL xas_count_star_members(kpts, ikpt, nstar)
            CALL xas_star_member_weight(kpts, ikpt, wk_star)
            weight_sum_star = weight_sum_star + REAL(nstar)*wk_star
         ELSE
            weight_sum_star = weight_sum_star + kpts%wtkpt(ikpt)
         END IF
      END DO
      transition_min = HUGE(1.0)
      transition_max = -HUGE(1.0)
      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) /= xas_absorber_z) CYCLE
         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_extract_core_states(atoms, itype, xas_edge_name, vTot%mt(1:atoms%jri(itype), 0, itype, 1), core_states)
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_extract_core_states", unit=xas_debug_unit)
         IF (SIZE(core_states) < 1) THEN
            WRITE(error_message, '(a,a,a,i0,a,i0)') "No core state found for requested XAS edge ", TRIM(xas_edge_name), &
                                                    " in absorber Z=", xas_absorber_z, " atom type ", itype
            CALL juDFT_error(TRIM(error_message), calledby="m_xas_driver")
         END IF
         DO jsp_loop = 1, n_spin_channels
            jsp = MERGE(1, jsp_loop, noco%l_noco)
            DO ikpt_i = 1, SIZE(fmpi%k_list)
               ikpt = fmpi%k_list(ikpt_i)
               IF (kpts%wtkpt(ikpt) <= 0.0) CYCLE
               nbands = results%neig(ikpt, jsp)
               IF (nbands <= 0) CYCLE
               DO i_band = 1, nbands
                  occ = results%w_iks(i_band, ikpt, jsp)/kpts%wtkpt(ikpt)
                  IF (1.0 - occ <= 1.0e-10) CYCLE
                  transition_min = MIN(transition_min, results%eig(i_band, ikpt, jsp) - core_states(1)%energy)
                  transition_max = MAX(transition_max, results%eig(i_band, ikpt, jsp) - core_states(1)%energy)
               END DO
            END DO
         END DO
         DEALLOCATE(core_states)
      END DO
      CALL xas_allreduce_transition_window(fmpi, transition_min, transition_max)
      IF (transition_min > transition_max) THEN
         CALL juDFT_error("No empty final-state bands found for hardwired XAS test", calledby="m_xas_driver")
      END IF
      transition_padding = MAX(5.0*xas_test_eta, 0.05)
      transition_min = transition_min - transition_padding
      transition_max = transition_max + transition_padding

      ALLOCATE(energy_grid(xas_test_n_grid), intensity(xas_test_n_grid), SOURCE=0.0)
      DO i_band = 1, xas_test_n_grid
         energy_grid(i_band) = transition_min + (transition_max - transition_min)*REAL(i_band - 1)/REAL(xas_test_n_grid - 1)
      END DO

      l_real = sym%invs .AND. (.NOT. noco%l_soc) .AND. (.NOT. noco%l_noco) .AND. atoms%n_hia == 0

      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) /= xas_absorber_z) CYCLE
         DO ispin = 1, input%jspins
            CALL genMTBasis(atoms, enpara, vTot, fmpi, itype, ispin, usdus, &
                            f(:, :, 0:, ispin), g(:, :, 0:, ispin), flo(:, :, :, ispin), l_writeArg=.FALSE.)
         END DO

         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_extract_core_states(atoms, itype, xas_edge_name, vTot%mt(1:atoms%jri(itype), 0, itype, 1), core_states)
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_extract_core_states", unit=xas_debug_unit)
         IF (SIZE(core_states) < 1) THEN
            WRITE(error_message, '(a,a,a,i0,a,i0)') "No core state found for requested XAS edge ", TRIM(xas_edge_name), &
                                                    " in absorber Z=", xas_absorber_z, " atom type ", itype
            CALL juDFT_error(TRIM(error_message), calledby="m_xas_driver")
         END IF
         IF (l_root .AND. xas_debug_angular_sumrule .AND. (.NOT. l_xas_angular_sumrule_printed)) THEN
            CALL xas_print_angular_sumrule(core_states(1)%lc, core_states(1)%twice_j, xas_debug_unit)
            l_xas_angular_sumrule_printed = .TRUE.
         END IF

         CALL radfun%generate_radial_functions(atoms, input, enpara, fmpi, vTot, itype)
         max_order = MAXVAL(radfun%n_r(0:atoms%lmax(itype)))
         ALLOCATE(radial_xas(max_order, 0:atoms%lmaxd, input%jspins), SOURCE=0.0)
         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_radial_dipole_integrals(atoms, itype, radfun, core_states(1)%p_core, radial_xas)
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_radial_dipole_integrals", unit=xas_debug_unit)

         lmax_xas = atoms%lmax(itype)
         IF (l_spinor_abc) THEN
            ! abc%calc_abc returns local MT spin components in noco:
            ! abc_local(tau) = sum_s conjg(U(s,tau)) abc_global(s).
            ! Rotate the global core spin-angular coefficients with the same
            ! U^\dagger convention before contracting with abc_local.
            spin_frame_transform = CONJG(TRANSPOSE(nococonv%umat(itype)))
         END IF
         DO jsp_loop = 1, n_spin_channels
            jsp = MERGE(1, jsp_loop, l_spinor_abc)
            IF (l_spinor_abc) THEN
               n_local_spins = 2
            ELSE IF (input%jspins == 2) THEN
               n_local_spins = 2
            ELSE
               n_local_spins = 1
            END IF
            DO ikpt_i = 1, SIZE(fmpi%k_list)
               ikpt = fmpi%k_list(ikpt_i)
               IF (kpts%wtkpt(ikpt) <= 0.0) CYCLE
               nbands = results%neig(ikpt, jsp)
               IF (nbands <= 0) CYCLE

               ALLOCATE(ev_list(nbands))
               ev_list = [(i_band, i_band=1, nbands)]
               ALLOCATE(eig_band(nbands), occ_band(nbands))
               eig_band = results%eig(1:nbands, ikpt, jsp)
               occ_band = results%w_iks(1:nbands, ikpt, jsp)/kpts%wtkpt(ikpt)

               CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell, fmpi)
               nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
               CALL zMat%init(l_real, nbasfcn, nbands)
               CALL read_eig(eig_id, ikpt, jsp, list=ev_list, neig=nbands_read, zmat=zMat)
               IF (nbands_read < nbands) THEN
                  CALL juDFT_error("read_eig returned fewer bands than requested in hardwired XAS test", calledby="m_xas_driver")
               END IF

               ALLOCATE(abc_spin(n_local_spins))
               DO ispin = 1, n_local_spins
                  CALL abc_spin(ispin)%init(input, atoms, radfun%n_r, nbands, itype)
               END DO
               IF (l_spinor_abc) THEN
                  DO ispin = 1, n_local_spins
                     CALL abc_spin(ispin)%calc_abc(input, atoms, sym, cell, lapw, nbands, usdus, noco, nococonv, &
                                                   ispin, itype, zMat)
                  END DO
               ELSE IF (input%jspins == 2) THEN
                  CALL abc_spin(jsp_loop)%calc_abc(input, atoms, sym, cell, lapw, nbands, usdus, noco, nococonv, &
                                                   jsp_loop, itype, zMat)
	               ELSE
	                  CALL abc_spin(1)%calc_abc(input, atoms, sym, cell, lapw, nbands, usdus, noco, nococonv, &
	                                            1, itype, zMat)
	               END IF
               IF (l_root .AND. xas_debug_abc_star_compare .AND. l_spinor_abc .AND. &
                   xas_debug_is_direct_target_k(kpts%bk(:, ikpt)) .AND. &
                   n_abc_debug_direct_records < xas_debug_abc_star_max_records) THEN
	                  DO iatom_l = 1, atoms%neq(itype)
	                     IF (n_abc_debug_direct_records >= xas_debug_abc_star_max_records) EXIT
	                     CALL xas_debug_dump_abc_fingerprint(xas_abc_debug_unit, "DIRECT", abc_spin, atoms, sym, cell, &
	                                                         nococonv, itype, iatom_l, ikpt, ikpt, 1, kpts%bk(:, ikpt), &
	                                                         lmax_xas, eig_band)
	                     n_abc_debug_direct_records = n_abc_debug_direct_records + 1
	                  END DO
	               END IF

	               ALLOCATE(matrix(nbands, SIZE(core_states(1)%twice_mj)))
               IF (l_xas_debug_strength) THEN
                  ALLOCATE(matrix_debug(nbands, SIZE(core_states(1)%twice_mj)))
                  ALLOCATE(matrix_l0(nbands, SIZE(core_states(1)%twice_mj)))
                  ALLOCATE(matrix_l2(nbands, SIZE(core_states(1)%twice_mj)))
               END IF

               IF (xas_use_spatial_star) THEN
                  CALL xas_count_star_members(kpts, ikpt, nstar)
                  CALL xas_star_member_weight(kpts, ikpt, wk_star)
                  DO ikptf = 1, kpts%nkptf
                     IF (kpts%bkp(ikptf) /= ikpt) CYCLE
                     bksym = kpts%bksym(ikptf)
                     wk_current = wk_star

                     ALLOCATE(abc_star_spin(n_local_spins))
	                     IF (l_spinor_abc) THEN
	                        CALL xas_rotate_abc_star_member_spinor(abc_spin, atoms, sym, cell, nococonv, itype, bksym, &
	                                                               lmax_xas, abc_star_spin)
	                        IF (l_root .AND. xas_debug_abc_star_compare .AND. kpts%nkpt /= kpts%nkptf .AND. &
	                            bksym /= 1 .AND. n_abc_debug_star_records < xas_debug_abc_star_max_records) THEN
	                           DO iatom_l = 1, atoms%neq(itype)
	                              IF (n_abc_debug_star_records >= xas_debug_abc_star_max_records) EXIT
	                              CALL xas_debug_dump_abc_fingerprint(xas_abc_debug_unit, "STAR", abc_star_spin, atoms, sym, cell, &
	                                                                  nococonv, itype, iatom_l, ikpt, ikptf, bksym, &
	                                                                  kpts%bkf(:, ikptf), lmax_xas, eig_band)
	                              n_abc_debug_star_records = n_abc_debug_star_records + 1
	                           END DO
	                        END IF
	                     ELSE IF (noco%l_soc) THEN
	                        CALL xas_abort_missing_spinor_abc(input, noco, n_local_spins)
                     ELSE
                        DO ispin = 1, n_local_spins
                           CALL xas_rotate_abc_star_member(abc_spin(ispin), atoms, sym, cell, itype, bksym, lmax_xas, &
                                                           abc_star_spin(ispin))
                        END DO
                     END IF

                     IF (l_xas_debug_strength) THEN
                        DO iatom_l = 1, atoms%neq(itype)
                           DO i_pol = 1, xas_debug_n_pol
                              eps_cart_debug = CMPLX(0.0, 0.0)
                              eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                              CALL xas_cartesian_to_spherical(eps_cart_debug, eps_sph_debug)

                              IF (l_spinor_abc) THEN
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_debug, &
                                                                   spin_frame_transform=spin_frame_transform)
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_l0, final_l=0, &
                                                                   spin_frame_transform=spin_frame_transform)
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_l2, final_l=2, &
                                                                   spin_frame_transform=spin_frame_transform)
                              ELSE
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_debug)
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_l0, final_l=0)
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_l2, final_l=2)
                              END IF

                              debug_strength = xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                              debug_l0 = xas_debug_transition_strength(matrix_l0, occ_band, wk_current)
                              debug_l2 = xas_debug_transition_strength(matrix_l2, occ_band, wk_current)
                              xas_debug_strength_total(i_pol) = xas_debug_strength_total(i_pol) + debug_strength
                              xas_debug_strength_spin(i_pol, jsp_loop) = xas_debug_strength_spin(i_pol, jsp_loop) + debug_strength
                              IF (l_xas_debug_kpt_strength) THEN
                                 xas_debug_strength_kpt(i_pol, ikpt) = xas_debug_strength_kpt(i_pol, ikpt) + debug_strength
                              END IF
                              xas_debug_strength_l0(i_pol) = xas_debug_strength_l0(i_pol) + debug_l0
                              xas_debug_strength_l2(i_pol) = xas_debug_strength_l2(i_pol) + debug_l2
                           END DO
                        END DO
                     END IF

                     DO iatom_l = 1, atoms%neq(itype)
                        CALL xas_debug_clear_underflow(l_xas_debug_fp)
                        IF (l_spinor_abc) THEN
                           CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                             iatom_l, lmax_xas, matrix, &
                                                             spin_frame_transform=spin_frame_transform)
                        ELSE
                           CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                             iatom_l, lmax_xas, matrix)
                        END IF
                        CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_core_band_matrixelements", unit=xas_debug_unit)
                        CALL xas_debug_clear_underflow(l_xas_debug_fp)
                        CALL xas_accumulate_matrix_spectrum(energy_grid, eig_band, occ_band, wk_current, &
                                                            core_states(1)%energy, matrix, xas_test_eta, intensity, "gaussian")
                        ! The hardwired test can raise harmless underflow from negligible
                        ! Gaussian tails or tiny spectral products. Count it, then clear
                        ! only IEEE_UNDERFLOW so the temporary debug path does not leave
                        ! this benign flag set at program exit.
                        CALL xas_debug_count_and_clear_underflow(l_xas_debug_fp, n_underflow_spectrum)
                     END DO
                     DEALLOCATE(abc_star_spin)
                  END DO
               ELSE
                  wk_current = kpts%wtkpt(ikpt)
                  IF (l_xas_debug_strength) THEN
                     DO iatom_l = 1, atoms%neq(itype)
                        DO i_pol = 1, xas_debug_n_pol
                           eps_cart_debug = CMPLX(0.0, 0.0)
                           eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                           CALL xas_cartesian_to_spherical(eps_cart_debug, eps_sph_debug)

                           IF (l_spinor_abc) THEN
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_debug, &
                                                                spin_frame_transform=spin_frame_transform)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l0, final_l=0, &
                                                                spin_frame_transform=spin_frame_transform)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l2, final_l=2, &
                                                                spin_frame_transform=spin_frame_transform)
                           ELSE
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_debug)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l0, final_l=0)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l2, final_l=2)
                           END IF

                           debug_strength = xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                           debug_l0 = xas_debug_transition_strength(matrix_l0, occ_band, wk_current)
                           debug_l2 = xas_debug_transition_strength(matrix_l2, occ_band, wk_current)
                           xas_debug_strength_total(i_pol) = xas_debug_strength_total(i_pol) + debug_strength
                           xas_debug_strength_spin(i_pol, jsp_loop) = xas_debug_strength_spin(i_pol, jsp_loop) + debug_strength
                           IF (l_xas_debug_kpt_strength) THEN
                              xas_debug_strength_kpt(i_pol, ikpt) = xas_debug_strength_kpt(i_pol, ikpt) + debug_strength
                           END IF
                           xas_debug_strength_l0(i_pol) = xas_debug_strength_l0(i_pol) + debug_l0
                           xas_debug_strength_l2(i_pol) = xas_debug_strength_l2(i_pol) + debug_l2
                        END DO
                     END DO
                  END IF

                  DO iatom_l = 1, atoms%neq(itype)
                     CALL xas_debug_clear_underflow(l_xas_debug_fp)
                     IF (l_spinor_abc) THEN
                        CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                          iatom_l, lmax_xas, matrix, &
                                                          spin_frame_transform=spin_frame_transform)
                     ELSE
                        CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                          iatom_l, lmax_xas, matrix)
                     END IF
                     CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_core_band_matrixelements", unit=xas_debug_unit)
                     CALL xas_debug_clear_underflow(l_xas_debug_fp)
                     CALL xas_accumulate_matrix_spectrum(energy_grid, eig_band, occ_band, wk_current, &
                                                         core_states(1)%energy, matrix, xas_test_eta, intensity, "gaussian")
                     ! The hardwired test can raise harmless underflow from negligible
                     ! Gaussian tails or tiny spectral products. Count it, then clear
                     ! only IEEE_UNDERFLOW so the temporary debug path does not leave
                     ! this benign flag set at program exit.
                     CALL xas_debug_count_and_clear_underflow(l_xas_debug_fp, n_underflow_spectrum)
                  END DO
               END IF

               IF (l_xas_debug_strength) DEALLOCATE(matrix_debug, matrix_l0, matrix_l2)

               DEALLOCATE(matrix, abc_spin, eig_band, occ_band, ev_list)
            END DO
         END DO
         DEALLOCATE(radial_xas, core_states)
      END DO

      ! Rank-local k-point/star contributions are additive. Reduce them once
      ! after the k loops, then let rank 0 do all text/debug output.
      ALLOCATE(intensity_reduced(SIZE(intensity)), SOURCE=0.0)
      CALL mpi_sum_reduce(intensity, intensity_reduced, fmpi%mpi_comm)

      weight_sums_local = [weight_sum_parent, weight_sum_star]
      weight_sums_reduced = 0.0
      CALL mpi_sum_reduce(weight_sums_local, weight_sums_reduced, fmpi%mpi_comm)

      underflow_local = [n_underflow_spectrum]
      underflow_reduced = 0
      CALL mpi_sum_reduce(underflow_local, underflow_reduced, fmpi%mpi_comm)

      IF (l_xas_debug_strength) THEN
         ALLOCATE(xas_debug_strength_spin_reduced(SIZE(xas_debug_strength_spin, 1), &
                                                  SIZE(xas_debug_strength_spin, 2)), SOURCE=0.0)
         CALL mpi_sum_reduce(xas_debug_strength_spin, xas_debug_strength_spin_reduced, fmpi%mpi_comm)
         xas_debug_strength_total_reduced = 0.0
         xas_debug_strength_l0_reduced = 0.0
         xas_debug_strength_l2_reduced = 0.0
         CALL mpi_sum_reduce(xas_debug_strength_total, xas_debug_strength_total_reduced, fmpi%mpi_comm)
         CALL mpi_sum_reduce(xas_debug_strength_l0, xas_debug_strength_l0_reduced, fmpi%mpi_comm)
         CALL mpi_sum_reduce(xas_debug_strength_l2, xas_debug_strength_l2_reduced, fmpi%mpi_comm)
      END IF
      IF (l_xas_debug_kpt_strength) THEN
         ALLOCATE(xas_debug_strength_kpt_reduced(SIZE(xas_debug_strength_kpt, 1), &
                                                 SIZE(xas_debug_strength_kpt, 2)), SOURCE=0.0)
         CALL mpi_sum_reduce(xas_debug_strength_kpt, xas_debug_strength_kpt_reduced, fmpi%mpi_comm)
      END IF

      IF (l_root .AND. xas_debug_verbosity >= 1) THEN
         WRITE(*, '(a,i0,a,i0,a,l1,a,es12.4)') "XAS DEBUG star setup: nkpt=", kpts%nkpt, &
                                      " nkptf=", kpts%nkptf, " use_spatial_star=", xas_use_spatial_star, &
                                      " reduction_factor=", REAL(kpts%nkptf)/REAL(kpts%nkpt)
         WRITE(*, '(a,es18.10,a,es18.10,a,es12.4)') "XAS DEBUG star weight sums: selected=", weight_sums_reduced(1), &
                                            " expanded=", weight_sums_reduced(2), &
                                            " diff=", weight_sums_reduced(2) - weight_sums_reduced(1)
         WRITE(xas_debug_unit, '(a,i0,a,i0,a,l1,a,es12.4)') "XAS DEBUG star setup: nkpt=", kpts%nkpt, &
                                      " nkptf=", kpts%nkptf, " use_spatial_star=", xas_use_spatial_star, &
                                      " reduction_factor=", REAL(kpts%nkptf)/REAL(kpts%nkpt)
         WRITE(xas_debug_unit, '(a,es18.10,a,es18.10,a,es12.4)') "XAS DEBUG star weight sums: selected=", weight_sums_reduced(1), &
                                            " expanded=", weight_sums_reduced(2), &
                                            " diff=", weight_sums_reduced(2) - weight_sums_reduced(1)
         WRITE(*, '(a,i0,a,i0,a,i0,a)') "XAS DEBUG: summed ", n_absorber_atoms, &
                                        " total absorber atoms with Z=", xas_absorber_z, &
                                        " over ", n_absorber_types, " atom types"
         WRITE(xas_debug_unit, '(a,i0,a,i0,a,i0,a)') "XAS DEBUG: summed ", n_absorber_atoms, &
                                                     " total absorber atoms with Z=", xas_absorber_z, &
                                                     " over ", n_absorber_types, " atom types"
      END IF
      IF (l_root .AND. l_xas_debug_kpt_strength) THEN
         WRITE(*, '(a)') "XAS DEBUG kpt strength columns: ikpt weight Sx Sy Sz rel_xz=(Sx-Sz)/avg"
         WRITE(xas_debug_unit, '(a)') "XAS DEBUG kpt strength columns: ikpt weight Sx Sy Sz rel_xz=(Sx-Sz)/avg"
         DO ikpt = 1, kpts%nkpt
            debug_avg = SUM(xas_debug_strength_kpt_reduced(:, ikpt))/REAL(xas_debug_n_pol)
            debug_rel_xz = 0.0
            IF (ABS(debug_avg) > TINY(debug_avg)) THEN
               debug_rel_xz = (xas_debug_strength_kpt_reduced(1, ikpt) - xas_debug_strength_kpt_reduced(3, ikpt))/debug_avg
            END IF
            WRITE(*, '(a,i0,1x,5es16.8)') "XAS DEBUG kpt strength ", ikpt, kpts%wtkpt(ikpt), &
               xas_debug_strength_kpt_reduced(1, ikpt), xas_debug_strength_kpt_reduced(2, ikpt), &
               xas_debug_strength_kpt_reduced(3, ikpt), debug_rel_xz
            WRITE(xas_debug_unit, '(a,i0,1x,5es16.8)') "XAS DEBUG kpt strength ", ikpt, kpts%wtkpt(ikpt), &
               xas_debug_strength_kpt_reduced(1, ikpt), xas_debug_strength_kpt_reduced(2, ikpt), &
               xas_debug_strength_kpt_reduced(3, ikpt), debug_rel_xz
         END DO
      END IF
      IF (l_root .AND. l_xas_debug_strength) THEN
         DO jsp_loop = 1, n_spin_channels
            CALL xas_debug_print_strength_summary("spin-resolved", jsp_loop, &
                                                  xas_debug_strength_spin_reduced(:, jsp_loop), xas_debug_unit)
         END DO
         DO i_pol = 1, xas_debug_n_pol
            xas_debug_strength_cross_reduced(i_pol) = xas_debug_strength_total_reduced(i_pol) - &
                                                      xas_debug_strength_l0_reduced(i_pol) - &
                                                      xas_debug_strength_l2_reduced(i_pol)
            WRITE(*, '(a,a,a,es18.10,a,es18.10,a,es18.10,a,es18.10)') &
               "XAS DEBUG strength pol=", xas_debug_pol_label(i_pol), &
               " total=", xas_debug_strength_total_reduced(i_pol), &
               " l0=", xas_debug_strength_l0_reduced(i_pol), &
               " l2=", xas_debug_strength_l2_reduced(i_pol), &
               " cross_l0_l2=", xas_debug_strength_cross_reduced(i_pol)
            WRITE(xas_debug_unit, '(a,a,a,es18.10,a,es18.10,a,es18.10,a,es18.10)') &
               "XAS DEBUG strength pol=", xas_debug_pol_label(i_pol), &
               " total=", xas_debug_strength_total_reduced(i_pol), &
               " l0=", xas_debug_strength_l0_reduced(i_pol), &
               " l2=", xas_debug_strength_l2_reduced(i_pol), &
               " cross_l0_l2=", xas_debug_strength_cross_reduced(i_pol)
         END DO
         CALL xas_debug_print_strength_summary("total", 0, xas_debug_strength_total_reduced, xas_debug_unit)
      END IF

      IF (l_root) THEN
         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_write_spectrum_text(TRIM(output_filename), energy_grid, intensity_reduced, "Hartree")
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_write_spectrum_text", unit=xas_debug_unit)
         WRITE(*, '(a,a)') "XAS hardwired test wrote spectrum to ", TRIM(output_filename)
         WRITE(xas_debug_unit, '(a,a)') "XAS hardwired test wrote spectrum to ", TRIM(output_filename)
         WRITE(*, '(a,a)') "XAS hardwired test wrote debug log to ", TRIM(xas_debug_filename)
         WRITE(xas_debug_unit, '(a,a)') "XAS hardwired test wrote debug log to ", TRIM(xas_debug_filename)
         IF (xas_debug_verbosity >= 3 .AND. underflow_reduced(1) > 0) THEN
            WRITE(*, '(a,i0,a)') "XAS DEBUG SUMMARY: underflow occurred in xas_accumulate_matrix_spectrum ", &
                                 underflow_reduced(1), &
                                 " time(s); cleared as harmless Gaussian-tail/tiny-product underflow."
            WRITE(xas_debug_unit, '(a,i0,a)') "XAS DEBUG SUMMARY: underflow occurred in xas_accumulate_matrix_spectrum ", &
                                             underflow_reduced(1), &
                                             " time(s); cleared as harmless Gaussian-tail/tiny-product underflow."
	         END IF
	         CLOSE(xas_debug_unit)
	         IF (xas_abc_debug_unit /= -1) CLOSE(xas_abc_debug_unit)
	      END IF
	   END SUBROUTINE xas_hardwired_test_driver

   SUBROUTINE xas_abort_missing_spinor_abc(input, noco, n_local_spins)
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_noco),  INTENT(IN) :: noco
      INTEGER,       INTENT(IN) :: n_local_spins

      CHARACTER(LEN=320) :: error_message

      WRITE(error_message, '(a,i0,a,l1,a,l1,a,i0,a)') &
         "XAS spinor star rotation needs two local MT spinor abc components from abc%calc_abc; found input%jspins=", &
         input%jspins, ", noco%l_soc=", noco%l_soc, ", noco%l_noco=", noco%l_noco, &
         ", constructed local components=", n_local_spins, &
         ". The connected XAS spinor-star path is currently implemented for noco%l_noco; SOC without noco needs a separate abc construction."
	      CALL juDFT_error(TRIM(error_message), calledby="m_xas_driver")
	   END SUBROUTINE xas_abort_missing_spinor_abc

   SUBROUTINE xas_debug_write_abc_header(unit)
      INTEGER, INTENT(IN) :: unit

      WRITE(unit, '(a)') "# XAS abc star comparison fingerprints"
      WRITE(unit, '(a)') "# Run explicit full-k/no-symmetry with this flag enabled to get DIRECT records."
      WRITE(unit, '(a)') "# Run symmetry-reduced with the same setup to get STAR records."
      WRITE(unit, '(a)') "# Match records by full-zone k coordinates; compare l/tau/lm norms and complex checksums."
      WRITE(unit, '(a)') "# ABC_RECORD mode parent_ikpt full_ikpt bksym spatial_iop det_int itype atom_l time_reversal mapped_atom mapped_l Z kx ky kz"
      WRITE(unit, '(a)') "# ABC_BAND mode band eigenvalue_ha"
      WRITE(unit, '(a)') "# ABC_NORM mode tau l norm sum_re sum_im max_abs"
      WRITE(unit, '(a)') "# ABC_LM mode tau l m lm norm sum_re sum_im max_abs"
      WRITE(unit, '(a)') "# ABC_COEF mode tau l m lm iord band re im"
      WRITE(unit, '(a)') "# ABC_SU2 row re1 im1 re2 im2"
      WRITE(unit, '(a)') "# ABC_SLOCAL row re1 im1 re2 im2"
   END SUBROUTINE xas_debug_write_abc_header

   SUBROUTINE xas_debug_dump_abc_fingerprint(unit, mode, abc_spin, atoms, sym, cell, nococonv, itype, iatom_l, &
                                             ikpt_parent, ikpt_full, bksym, k_full, lmax_xas, eig_band)
      INTEGER,             INTENT(IN) :: unit
      CHARACTER(LEN=*),    INTENT(IN) :: mode
      TYPE(t_abc),         INTENT(IN) :: abc_spin(:)
      TYPE(t_atoms),       INTENT(IN) :: atoms
      TYPE(t_sym),         INTENT(IN) :: sym
      TYPE(t_cell),        INTENT(IN) :: cell
      TYPE(t_nococonv),    INTENT(IN) :: nococonv
      INTEGER,             INTENT(IN) :: itype, iatom_l, ikpt_parent, ikpt_full, bksym, lmax_xas
      REAL,                INTENT(IN) :: k_full(3)
      REAL,                INTENT(IN) :: eig_band(:)

      INTEGER :: spatial_iop, iatom, mapped_atom, mapped_l, tau, l, lm, m, det_int, irow, band
      LOGICAL :: l_time_reversal
      COMPLEX :: su_global(2, 2), su_local(2, 2)

      IF (unit == -1) RETURN
      IF (SIZE(abc_spin) < 2) THEN
         WRITE(unit, '(a,1x,a,1x,a,i0)') "# ABC_SKIP", TRIM(mode), "abc_spin_size=", SIZE(abc_spin)
         RETURN
      END IF
      IF (.NOT. ALLOCATED(abc_spin(1)%cof) .OR. .NOT. ALLOCATED(abc_spin(2)%cof)) THEN
         WRITE(unit, '(a,1x,a)') "# ABC_SKIP", TRIM(mode)//" missing allocated abc cof"
         RETURN
      END IF

      iatom = atoms%firstAtom(itype) + iatom_l - 1
      mapped_atom = iatom
      CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
      IF (ALLOCATED(sym%mapped_atom)) THEN
         IF (iatom >= 1 .AND. iatom <= atoms%nat) mapped_atom = sym%mapped_atom(spatial_iop, iatom)
      END IF
      IF (mapped_atom < 1 .OR. mapped_atom > atoms%nat) mapped_atom = iatom
      IF (mapped_atom >= 1 .AND. mapped_atom <= atoms%nat) THEN
         mapped_l = mapped_atom - atoms%firstAtom(atoms%itype(mapped_atom)) + 1
      ELSE
         mapped_l = -1
      END IF
      det_int = xas_debug_det3_int(sym%mrot(:, :, spatial_iop))

      WRITE(unit, '(a,1x,a,1x,7(i0,1x),l1,1x,i0,1x,i0,1x,i0,1x,3es22.14)') &
         "ABC_RECORD", TRIM(mode), ikpt_parent, ikpt_full, bksym, spatial_iop, det_int, itype, iatom_l, &
         l_time_reversal, mapped_atom, mapped_l, atoms%nz(itype), k_full

      CALL xas_su2_from_sym(sym, cell, spatial_iop, su_global)
      CALL xas_local_spin_transform(atoms, nococonv, iatom, mapped_atom, su_global, su_local)
      DO irow = 1, 2
         WRITE(unit, '(a,1x,i0,1x,4es22.14)') "ABC_SU2", irow, REAL(su_global(irow, 1)), AIMAG(su_global(irow, 1)), &
                                              REAL(su_global(irow, 2)), AIMAG(su_global(irow, 2))
         WRITE(unit, '(a,1x,i0,1x,4es22.14)') "ABC_SLOCAL", irow, REAL(su_local(irow, 1)), AIMAG(su_local(irow, 1)), &
                                                 REAL(su_local(irow, 2)), AIMAG(su_local(irow, 2))
      END DO
      DO band = 1, MIN(SIZE(eig_band), SIZE(abc_spin(1)%cof, 1))
         WRITE(unit, '(a,1x,a,1x,i0,1x,es22.14)') "ABC_BAND", TRIM(mode), band, eig_band(band)
      END DO

      DO tau = 1, 2
         DO l = 0, MIN(lmax_xas, UBOUND(abc_spin(tau)%n_r, 1))
            IF (l /= 0 .AND. l /= 2) CYCLE
            CALL xas_debug_write_abc_block(unit, mode, abc_spin(tau), tau, l, iatom_l)
            DO lm = l*l, l*(l + 2)
               m = lm - l*(l + 1)
               CALL xas_debug_write_abc_lm(unit, mode, abc_spin(tau), tau, l, m, lm, iatom_l)
            END DO
            CALL xas_debug_write_abc_coefficients(unit, mode, abc_spin(tau), tau, l, iatom_l)
         END DO
      END DO
      WRITE(unit, '(a)') "ABC_END"
   END SUBROUTINE xas_debug_dump_abc_fingerprint

   SUBROUTINE xas_debug_write_abc_block(unit, mode, abc, tau, l, iatom_l)
      INTEGER,          INTENT(IN) :: unit, tau, l, iatom_l
      CHARACTER(LEN=*), INTENT(IN) :: mode
      TYPE(t_abc),      INTENT(IN) :: abc

      INTEGER :: lm1, lm2
      REAL    :: norm_val, max_abs
      COMPLEX :: checksum

      lm1 = l*l
      lm2 = l*(l + 2)
      IF (lm2 > UBOUND(abc%cof, 2)) RETURN
      norm_val = SQRT(SUM(ABS(abc%cof(:, lm1:lm2, :, iatom_l))**2))
      max_abs = MAXVAL(ABS(abc%cof(:, lm1:lm2, :, iatom_l)))
      checksum = SUM(abc%cof(:, lm1:lm2, :, iatom_l))
      WRITE(unit, '(a,1x,a,1x,i0,1x,i0,1x,4es22.14)') "ABC_NORM", TRIM(mode), tau, l, norm_val, &
                                                        REAL(checksum), AIMAG(checksum), max_abs
   END SUBROUTINE xas_debug_write_abc_block

   SUBROUTINE xas_debug_write_abc_lm(unit, mode, abc, tau, l, m, lm, iatom_l)
      INTEGER,          INTENT(IN) :: unit, tau, l, m, lm, iatom_l
      CHARACTER(LEN=*), INTENT(IN) :: mode
      TYPE(t_abc),      INTENT(IN) :: abc

      REAL    :: norm_val, max_abs
      COMPLEX :: checksum

      IF (lm > UBOUND(abc%cof, 2)) RETURN
      norm_val = SQRT(SUM(ABS(abc%cof(:, lm, :, iatom_l))**2))
      max_abs = MAXVAL(ABS(abc%cof(:, lm, :, iatom_l)))
      checksum = SUM(abc%cof(:, lm, :, iatom_l))
      WRITE(unit, '(a,1x,a,1x,5(i0,1x),4es22.14)') "ABC_LM", TRIM(mode), tau, l, m, lm, iatom_l, &
                                                     norm_val, REAL(checksum), AIMAG(checksum), max_abs
   END SUBROUTINE xas_debug_write_abc_lm

   SUBROUTINE xas_debug_write_abc_coefficients(unit, mode, abc, tau, l, iatom_l)
      INTEGER,          INTENT(IN) :: unit, tau, l, iatom_l
      CHARACTER(LEN=*), INTENT(IN) :: mode
      TYPE(t_abc),      INTENT(IN) :: abc

      INTEGER :: lm, m, iord, band
      COMPLEX :: coeff

      IF (l > UBOUND(abc%n_r, 1)) RETURN
      DO lm = l*l, l*(l + 2)
         IF (lm > UBOUND(abc%cof, 2)) CYCLE
         m = lm - l*(l + 1)
         DO iord = 1, MIN(abc%n_r(l), SIZE(abc%cof, 3))
            DO band = 1, SIZE(abc%cof, 1)
               coeff = abc%cof(band, lm, iord, iatom_l)
               WRITE(unit, '(a,1x,a,1x,6(i0,1x),2es22.14)') "ABC_COEF", TRIM(mode), tau, l, m, lm, iord, band, &
                                                             REAL(coeff), AIMAG(coeff)
            END DO
         END DO
      END DO
   END SUBROUTINE xas_debug_write_abc_coefficients

   LOGICAL FUNCTION xas_debug_is_direct_target_k(kvec) RESULT(l_match)
      REAL, INTENT(IN) :: kvec(3)

      REAL, PARAMETER :: tol = 1.0e-7
      REAL, PARAMETER :: targets(3, 20) = RESHAPE([REAL :: &
         3.75000000000000E-01, 5.41666666666667E-01, 5.41666666666667E-01, &
         4.58333333333333E-01, 4.58333333333333E-01, 6.25000000000000E-01, &
         5.41666666666667E-01, 5.41666666666667E-01, 5.41666666666667E-01, &
         5.41666666666667E-01, 3.75000000000000E-01, 5.41666666666667E-01, &
         4.58333333333333E-01, 6.25000000000000E-01, 4.58333333333333E-01, &
         4.58333333333333E-01, 4.58333333333333E-01, 4.58333333333333E-01, &
         5.41666666666667E-01, 5.41666666666667E-01, 3.75000000000000E-01, &
         4.58333333333333E-01, 5.41666666666667E-01, 5.41666666666667E-01, &
         3.75000000000000E-01, 4.58333333333333E-01, 7.08333333333333E-01, &
         5.41666666666667E-01, 5.41666666666667E-01, 6.25000000000000E-01, &
         4.58333333333333E-01, 7.08333333333333E-01, 4.58333333333333E-01, &
         7.08333333333333E-01, 3.75000000000000E-01, 4.58333333333333E-01, &
         2.91666666666667E-01, 6.25000000000000E-01, 5.41666666666667E-01, &
         5.41666666666667E-01, 2.91666666666667E-01, 5.41666666666667E-01, &
         4.58333333333333E-01, 4.58333333333333E-01, 3.75000000000000E-01, &
         3.75000000000000E-01, 5.41666666666667E-01, 6.25000000000000E-01, &
         4.58333333333333E-01, 5.41666666666667E-01, 6.25000000000000E-01, &
         4.58333333333333E-01, 6.25000000000000E-01, 5.41666666666667E-01, &
         6.25000000000000E-01, 3.75000000000000E-01, 4.58333333333333E-01, &
         3.75000000000000E-01, 6.25000000000000E-01, 5.41666666666667E-01], [3, 20])

      INTEGER :: i
      REAL    :: diff(3)

      l_match = .FALSE.
      DO i = 1, SIZE(targets, 2)
         diff = kvec - targets(:, i)
         diff = diff - REAL(NINT(diff))
         IF (MAXVAL(ABS(diff)) < tol) THEN
            l_match = .TRUE.
            EXIT
         END IF
      END DO
   END FUNCTION xas_debug_is_direct_target_k

   SUBROUTINE xas_allreduce_transition_window(fmpi, transition_min, transition_max)
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      REAL,        INTENT(INOUT) :: transition_min
      REAL,        INTENT(INOUT) :: transition_max

#ifdef CPP_MPI
      INTEGER :: ierr
      REAL    :: local_min, local_max

      local_min = transition_min
      local_max = transition_max
      CALL MPI_ALLREDUCE(local_min, transition_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, fmpi%mpi_comm, ierr)
      IF (ierr /= 0) CALL juDFT_error("MPI_ALLREDUCE failed for XAS transition minimum", calledby="m_xas_driver")
      CALL MPI_ALLREDUCE(local_max, transition_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, fmpi%mpi_comm, ierr)
      IF (ierr /= 0) CALL juDFT_error("MPI_ALLREDUCE failed for XAS transition maximum", calledby="m_xas_driver")
#else
      ! Serial/non-MPI builds already have the complete transition window.
      transition_min = transition_min
      transition_max = transition_max
#endif
   END SUBROUTINE xas_allreduce_transition_window

   SUBROUTINE xas_debug_open_log(kpts, l_use_spatial_star, unit, filename)
      TYPE(t_kpts),      INTENT(IN)  :: kpts
      LOGICAL,           INTENT(IN)  :: l_use_spatial_star
      INTEGER,           INTENT(OUT) :: unit
      CHARACTER(LEN=*),  INTENT(OUT) :: filename

      CHARACTER(LEN=16) :: reduction_label
      INTEGER           :: i_char
      REAL              :: reduction_factor

      IF (kpts%nkpt <= 0) THEN
         CALL juDFT_error("Invalid k-point count in xas_debug_open_log", calledby="m_xas_driver")
      END IF

      reduction_factor = REAL(kpts%nkptf)/REAL(kpts%nkpt)
      WRITE(reduction_label, '(f8.4)') reduction_factor
      reduction_label = ADJUSTL(reduction_label)
      DO i_char = 1, LEN_TRIM(reduction_label)
         IF (reduction_label(i_char:i_char) == ".") reduction_label(i_char:i_char) = "p"
      END DO

      WRITE(filename, '("xas_debug_nkpt",i0,"_nkptf",i0,"_star",l1,"_red",a,".txt")') &
         kpts%nkpt, kpts%nkptf, l_use_spatial_star, TRIM(reduction_label)
      OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE")
   END SUBROUTINE xas_debug_open_log

   SUBROUTINE xas_debug_print_spinor_star(sym, cell, atoms, nococonv, absorber_z, unit)
      TYPE(t_sym),      INTENT(IN) :: sym
      TYPE(t_cell),     INTENT(IN) :: cell
      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      INTEGER,          INTENT(IN) :: absorber_z
      INTEGER,          INTENT(IN) :: unit

      COMPLEX :: su_global(2, 2), su_local(2, 2), su_check(2, 2), det_su
      REAL    :: r_cart(3, 3), unitarity_norm
      INTEGER :: identity_rot(3, 3), iop_list(2), n_iop, iop, itype, det_rot
      INTEGER :: parent_atom, mapped_atom, i, j
      LOGICAL :: found_absorber

      identity_rot = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      found_absorber = .FALSE.
      parent_atom = 1
      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) == absorber_z) THEN
            parent_atom = atoms%firstAtom(itype)
            found_absorber = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT. found_absorber) RETURN

      iop_list = 1
      n_iop = 1
      DO iop = 1, sym%nop
         det_rot = xas_debug_det3_int(sym%mrot(:, :, iop))
         IF (det_rot == 1 .AND. ANY(sym%mrot(:, :, iop) /= identity_rot)) THEN
            n_iop = 2
            iop_list(2) = iop
            EXIT
         END IF
      END DO
      IF (n_iop == 1) THEN
         DO iop = 1, sym%nop
            IF (ANY(sym%mrot(:, :, iop) /= identity_rot)) THEN
               n_iop = 2
               iop_list(2) = iop
               EXIT
            END IF
         END DO
      END IF

      CALL xas_debug_write_line(unit, "XAS DEBUG SPINOR STAR: diagnostic-only SU(2) matrices")
      DO i = 1, n_iop
         iop = iop_list(i)
         mapped_atom = parent_atom
         IF (ALLOCATED(sym%mapped_atom)) mapped_atom = sym%mapped_atom(iop, parent_atom)

         CALL xas_cart_rotation_from_sym(sym, cell, iop, r_cart)
         CALL xas_su2_from_sym(sym, cell, iop, su_global)
         CALL xas_local_spin_transform(atoms, nococonv, parent_atom, mapped_atom, su_global, su_local)

         su_check = MATMUL(CONJG(TRANSPOSE(su_global)), su_global)
         su_check(1, 1) = su_check(1, 1) - CMPLX(1.0, 0.0)
         su_check(2, 2) = su_check(2, 2) - CMPLX(1.0, 0.0)
         unitarity_norm = MAXVAL(ABS(su_check))
         det_su = su_global(1, 1)*su_global(2, 2) - su_global(1, 2)*su_global(2, 1)

         WRITE(*, '(a,i0,a,i0,a,i0)') "XAS DEBUG SPINOR STAR iop=", iop, " parent_atom=", parent_atom, &
                                      " mapped_atom=", mapped_atom
         IF (unit > 0) WRITE(unit, '(a,i0,a,i0,a,i0)') "XAS DEBUG SPINOR STAR iop=", iop, &
                                                        " parent_atom=", parent_atom, " mapped_atom=", mapped_atom
         DO j = 1, 3
            WRITE(*, '(a,i0,a,3es16.8)') "XAS DEBUG SPINOR STAR R_cart row ", j, ":", r_cart(j, :)
            IF (unit > 0) WRITE(unit, '(a,i0,a,3es16.8)') "XAS DEBUG SPINOR STAR R_cart row ", j, ":", r_cart(j, :)
         END DO
         CALL xas_debug_write_complex2(unit, "XAS DEBUG SPINOR STAR SU2 global row 1:", su_global(1, :))
         CALL xas_debug_write_complex2(unit, "XAS DEBUG SPINOR STAR SU2 global row 2:", su_global(2, :))
         WRITE(*, '(a,es12.4,a,2es16.8)') "XAS DEBUG SPINOR STAR SU2 unitarity_norm=", unitarity_norm, &
                                          " det=", REAL(det_su), AIMAG(det_su)
         IF (unit > 0) WRITE(unit, '(a,es12.4,a,2es16.8)') &
            "XAS DEBUG SPINOR STAR SU2 unitarity_norm=", unitarity_norm, " det=", REAL(det_su), AIMAG(det_su)
         CALL xas_debug_write_complex2(unit, "XAS DEBUG SPINOR STAR SU2 local row 1:", su_local(1, :))
         CALL xas_debug_write_complex2(unit, "XAS DEBUG SPINOR STAR SU2 local row 2:", su_local(2, :))
      END DO
   END SUBROUTINE xas_debug_print_spinor_star

   INTEGER FUNCTION xas_debug_det3_int(mat)
      INTEGER, INTENT(IN) :: mat(3, 3)

      xas_debug_det3_int = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(3, 2)*mat(2, 3)) + &
                           mat(1, 2)*(mat(2, 3)*mat(3, 1) - mat(2, 1)*mat(3, 3)) + &
                           mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))
   END FUNCTION xas_debug_det3_int

   SUBROUTINE xas_debug_write_line(unit, text)
      INTEGER,          INTENT(IN) :: unit
      CHARACTER(LEN=*), INTENT(IN) :: text

      WRITE(*, '(a)') TRIM(text)
      IF (unit > 0) WRITE(unit, '(a)') TRIM(text)
   END SUBROUTINE xas_debug_write_line

   SUBROUTINE xas_debug_write_complex2(unit, label, values)
      INTEGER,          INTENT(IN) :: unit
      CHARACTER(LEN=*), INTENT(IN) :: label
      COMPLEX,          INTENT(IN) :: values(2)

      WRITE(*, '(a,2(1x,"(",es13.5,",",es13.5,")"))') TRIM(label), &
         REAL(values(1)), AIMAG(values(1)), REAL(values(2)), AIMAG(values(2))
      IF (unit > 0) WRITE(unit, '(a,2(1x,"(",es13.5,",",es13.5,")"))') TRIM(label), &
         REAL(values(1)), AIMAG(values(1)), REAL(values(2)), AIMAG(values(2))
   END SUBROUTINE xas_debug_write_complex2

   SUBROUTINE xas_debug_print_strength_summary(label, jsp_loop, strength, unit)
      CHARACTER(LEN=*), INTENT(IN) :: label
      INTEGER,          INTENT(IN) :: jsp_loop
      REAL,             INTENT(IN) :: strength(:)
      INTEGER,          INTENT(IN) :: unit

      REAL :: avg, rel_xy, rel_xz, rel_yz

      avg = SUM(strength(1:xas_debug_n_pol))/REAL(xas_debug_n_pol)
      rel_xy = 0.0
      rel_xz = 0.0
      rel_yz = 0.0
      IF (ABS(avg) > TINY(avg)) THEN
         rel_xy = (strength(1) - strength(2))/avg
         rel_xz = (strength(1) - strength(3))/avg
         rel_yz = (strength(2) - strength(3))/avg
      END IF

      IF (jsp_loop > 0) THEN
         WRITE(*, '(a,a,a,i0,a,7es18.10)') "XAS DEBUG strength summary ", TRIM(label), &
            " jsp=", jsp_loop, " Sx Sy Sz avg x_minus_y_over_avg x_minus_z_over_avg y_minus_z_over_avg ", &
            strength(1), strength(2), strength(3), avg, rel_xy, rel_xz, rel_yz
         WRITE(unit, '(a,a,a,i0,a,7es18.10)') "XAS DEBUG strength summary ", TRIM(label), &
            " jsp=", jsp_loop, " Sx Sy Sz avg x_minus_y_over_avg x_minus_z_over_avg y_minus_z_over_avg ", &
            strength(1), strength(2), strength(3), avg, rel_xy, rel_xz, rel_yz
      ELSE
         WRITE(*, '(a,a,a,7es18.10)') "XAS DEBUG strength summary ", TRIM(label), &
            " Sx Sy Sz avg x_minus_y_over_avg x_minus_z_over_avg y_minus_z_over_avg ", &
            strength(1), strength(2), strength(3), avg, rel_xy, rel_xz, rel_yz
         WRITE(unit, '(a,a,a,7es18.10)') "XAS DEBUG strength summary ", TRIM(label), &
            " Sx Sy Sz avg x_minus_y_over_avg x_minus_z_over_avg y_minus_z_over_avg ", &
            strength(1), strength(2), strength(3), avg, rel_xy, rel_xz, rel_yz
      END IF
   END SUBROUTINE xas_debug_print_strength_summary

   REAL FUNCTION xas_debug_transition_strength(matrix, occ_band, wk) RESULT(strength)
      COMPLEX, INTENT(IN) :: matrix(:, :)
      REAL,    INTENT(IN) :: occ_band(:)
      REAL,    INTENT(IN) :: wk

      INTEGER :: i_band, i_mj
      REAL    :: empty_factor

      strength = 0.0
      DO i_band = 1, MIN(SIZE(matrix, 1), SIZE(occ_band))
         empty_factor = 1.0 - occ_band(i_band)
         IF (empty_factor <= 1.0e-10) CYCLE
         DO i_mj = 1, SIZE(matrix, 2)
            strength = strength + wk*empty_factor*ABS(matrix(i_band, i_mj))**2
         END DO
      END DO
   END FUNCTION xas_debug_transition_strength

   SUBROUTINE xas_debug_clear_underflow(l_debug)
      LOGICAL, INTENT(IN) :: l_debug

      IF (l_debug) CALL IEEE_SET_FLAG(IEEE_UNDERFLOW, .FALSE.)
   END SUBROUTINE xas_debug_clear_underflow

   SUBROUTINE xas_debug_report_underflow(l_debug, stage_name, l_clear_after, unit)
      LOGICAL,          INTENT(IN) :: l_debug
      CHARACTER(LEN=*), INTENT(IN) :: stage_name
      LOGICAL, OPTIONAL, INTENT(IN) :: l_clear_after
      INTEGER, OPTIONAL, INTENT(IN) :: unit

      LOGICAL :: l_underflow
      LOGICAL :: l_write

      IF (.NOT. l_debug) RETURN

      CALL IEEE_GET_FLAG(IEEE_UNDERFLOW, l_underflow)
      IF (l_underflow) THEN
         l_write = .TRUE.
         IF (PRESENT(unit)) l_write = unit > 0
         IF (l_write) WRITE(*, '(a,a)') "XAS DEBUG: underflow occurred in ", TRIM(stage_name)
         IF (PRESENT(unit)) THEN
            IF (unit > 0) WRITE(unit, '(a,a)') "XAS DEBUG: underflow occurred in ", TRIM(stage_name)
         END IF
         IF (PRESENT(l_clear_after)) THEN
            IF (l_clear_after) CALL IEEE_SET_FLAG(IEEE_UNDERFLOW, .FALSE.)
         END IF
      END IF
   END SUBROUTINE xas_debug_report_underflow

   SUBROUTINE xas_debug_count_and_clear_underflow(l_debug, n_underflow)
      LOGICAL, INTENT(IN)    :: l_debug
      INTEGER, INTENT(INOUT) :: n_underflow

      LOGICAL :: l_underflow

      IF (.NOT. l_debug) RETURN

      CALL IEEE_GET_FLAG(IEEE_UNDERFLOW, l_underflow)
      IF (l_underflow) THEN
         n_underflow = n_underflow + 1
         CALL IEEE_SET_FLAG(IEEE_UNDERFLOW, .FALSE.)
      END IF
   END SUBROUTINE xas_debug_count_and_clear_underflow

END MODULE m_xas_driver
