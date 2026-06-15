!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_driver
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_UNDERFLOW, IEEE_GET_FLAG, IEEE_SET_FLAG
   USE m_eig66_io, ONLY: read_eig
   USE m_genMTBasis, ONLY: genMTBasis
   USE m_juDFT, ONLY: juDFT_error
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
   USE m_xas_angular, ONLY: xas_cartesian_to_spherical
   USE m_xas_core, ONLY: t_xas_core_state, xas_extract_core_states
   USE m_xas_io, ONLY: xas_write_spectrum_text
   USE m_xas_matrixelements, ONLY: xas_core_band_matrixelements
   USE m_xas_radial, ONLY: xas_radial_dipole_integrals
   USE m_xas_spectrum, ONLY: xas_accumulate_matrix_spectrum
   USE m_xas_symmetry, ONLY: xas_count_star_members, xas_star_member_weight, xas_star_operation, &
                             xas_cart_rotation_from_sym, xas_rotate_lab_polarization_for_parent, &
                             xas_rotate_abc_star_member
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: xas_test_atomic_number = 26
   INTEGER, PARAMETER :: xas_test_n_grid = 401
   INTEGER, PARAMETER :: xas_debug_n_pol = 3
   ! Debug verbosity:
   !   0: only final spectrum filename
   !   1: compact atom/spin/total/l-channel strength summary
   !   2: level 1 plus compact k-point contribution table
   !   3: level 2 plus floating-point underflow diagnostics
   INTEGER, PARAMETER :: xas_debug_verbosity = 1
   ! Temporary first spatial-star implementations for the hardwired driver.
   !   0: no star reconstruction, selected k only
   !   1: diagnostic/incomplete, parent eigenvectors with rotated polarization only
   !   2: full-zone read_eig transformation with lab polarization
   !      (diagnostic only: read_eig's current wavefunction transform appears
   !       incomplete for LO/abc usage and is not production-safe)
   !   3: preferred current implementation for collinear no-SOC spatial-star
   !      XAS: parent eigenvectors with Wigner-D transformed MT coefficients
   !      and lab polarization
   ! SOC and noncollinear spinors require additional spinor-rotation considerations.
   INTEGER, PARAMETER :: xas_star_mode = 3
   LOGICAL, PARAMETER :: l_xas_use_spatial_star = xas_star_mode > 0
   LOGICAL, PARAMETER :: l_xas_star_sanity_print = .TRUE.
   LOGICAL, PARAMETER :: l_xas_debug_fullk_contrib_table = .TRUE.
   LOGICAL, PARAMETER :: l_xas_debug_compare_rotation_modes = .FALSE.
   LOGICAL, PARAMETER :: l_xas_debug_mode2_norms = .TRUE.
   CHARACTER(LEN=1), PARAMETER :: xas_debug_pol_label(xas_debug_n_pol) = [CHARACTER(LEN=1) :: "x", "y", "z"]
   ! Hardwired Gaussian broadening for the smoke test, in Hartree.
   ! Useful trial values: 0.01 Ha (sharper), 0.03 Ha (default), 0.05 Ha (smoother).
   REAL,    PARAMETER :: xas_test_eta = 0.03

   PUBLIC :: xas_hardwired_test_driver

CONTAINS

   SUBROUTINE xas_hardwired_test_driver(eig_id, fmpi, input, kpts, atoms, sym, cell, noco, nococonv, enpara, vTot, results)
      ! Temporary internal XAS smoke-test driver.
      !
      ! Hardwired choices:
      !   all Fe atom types (Z=26), L3 edge, z polarization, Gaussian
      !   eta = xas_test_eta Ha, output file labelled by edge, polarization,
      !   and broadening. The absorbing-atom selection is temporary Fe-test
      !   behavior and must become XML/input controlled in the real driver.
      !
      ! The current test spectrum is a per-selected-absorbing-atom local
      ! muffin-tin XAS signal. It is k-weighted and written in arbitrary units,
      ! but it is not normalized per cell volume, film area, film thickness, or
      ! number of equivalent atoms. The routine is serial-only and is not yet
      ! suitable for quantitative comparison between bulk and film calculations.
      !
      ! The local dipole approximation intentionally neglects interstitial and
      ! vacuum contributions. This is appropriate for core-level local XAS as a
      ! first implementation.
      !
      ! Spatial-star handling is currently guarded by xas_star_mode. Mode 3
      ! reuses parent irreducible-k eigenvectors, rotates the MT angular
      ! coefficients with Wigner-D matrices, and keeps the lab polarization
      ! fixed for each full k-star member. SOC/noncollinear calculations need
      ! additional spinor rotations before this can become production machinery.
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
      COMPLEX :: eps_parent_cart(3), eps_parent_sph(-1:1), eps_star_cart_debug(3), eps_star_sph_debug(-1:1)
      COMPLEX, ALLOCATABLE :: matrix(:, :), matrix_debug(:, :), matrix_l0(:, :), matrix_l2(:, :)
      REAL, ALLOCATABLE :: energy_grid(:), intensity(:), radial_xas(:, :, :)
      REAL, ALLOCATABLE :: f(:, :, :, :), g(:, :, :, :), flo(:, :, :, :)
      REAL, ALLOCATABLE :: eig_band(:), occ_band(:)
      REAL, ALLOCATABLE :: xas_debug_strength_kpt(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_spin(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_fullk(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_fullk_l0(:, :)
      REAL, ALLOCATABLE :: xas_debug_strength_fullk_l2(:, :)
      INTEGER, ALLOCATABLE :: ev_list(:)

      CHARACTER(LEN=2)   :: edge_name
      CHARACTER(LEN=1)   :: pol_label
      CHARACTER(LEN=5)   :: eta_label
      CHARACTER(LEN=128) :: output_filename
      CHARACTER(LEN=128) :: xas_debug_filename
      CHARACTER(LEN=160) :: xas_contrib_filename
      INTEGER :: ikpt_i, ikpt, ikptf, bksym, spatial_iop, jsp_loop, jsp, ispin, nbands, nbands_read, itype
      INTEGER :: n_spin_channels, n_local_spins, max_order, nbasfcn, n_apw_rows, lmax_xas, i_band
      INTEGER :: n_underflow_spectrum, i_char, i_pol, iatom_l, n_fe_types, n_fe_atoms, nstar, n_sanity_print
      INTEGER :: xas_debug_unit
      REAL    :: xas_debug_strength_l0(xas_debug_n_pol), xas_debug_strength_l2(xas_debug_n_pol)
      REAL    :: xas_debug_strength_total(xas_debug_n_pol), xas_debug_strength_cross(xas_debug_n_pol)
      REAL    :: transition_min, transition_max, transition_padding, occ, debug_strength, debug_l0, debug_l2
      REAL    :: debug_avg, debug_rel_xz, z_norm, z_apw_norm, z_lo_norm, z_max_abs, abc_max_abs, abc_sum_abs2
      REAL    :: wk_current, wk_star, weight_sum_parent, weight_sum_star, reconstructed_k(3), rrot(3, 3)
      REAL    :: xas_debug_strength_mode_b(xas_debug_n_pol)
      LOGICAL :: l_real, l_xas_debug_fp, l_xas_debug_strength, l_xas_debug_kpt_strength, l_time_reversal

      l_xas_debug_fp = xas_debug_verbosity >= 3
      l_xas_debug_strength = xas_debug_verbosity >= 1
      l_xas_debug_kpt_strength = xas_debug_verbosity >= 2
      n_underflow_spectrum = 0
      edge_name = "L3"
      xas_debug_strength_l0 = 0.0
      xas_debug_strength_l2 = 0.0
      xas_debug_strength_total = 0.0
      xas_debug_strength_cross = 0.0
      xas_debug_strength_mode_b = 0.0
      n_fe_types = 0
      n_fe_atoms = 0
      n_sanity_print = 0
      weight_sum_parent = 0.0
      weight_sum_star = 0.0
      xas_debug_unit = -1

      IF (fmpi%isize /= 1) THEN
         CALL juDFT_error("xas_hardwired_test_driver is serial-only for now", calledby="m_xas_driver")
      END IF
      IF (.NOT. ALLOCATED(results%w_iks)) THEN
         CALL juDFT_error("results%w_iks is not allocated in xas_hardwired_test_driver", calledby="m_xas_driver")
      END IF

      CALL xas_debug_open_log(kpts, l_xas_use_spatial_star, xas_debug_unit, xas_debug_filename)
      xas_contrib_filename = "xas_kcontrib_"//TRIM(xas_debug_filename)

      IF (xas_debug_verbosity >= 3) THEN
         WRITE(*, '(a,i0)') "XAS DEBUG atom types: ntype=", atoms%ntype
         WRITE(xas_debug_unit, '(a,i0)') "XAS DEBUG atom types: ntype=", atoms%ntype
      END IF
      DO itype = 1, atoms%ntype
         IF (xas_debug_verbosity >= 3) THEN
            WRITE(*, '(a,i0,a,i0,a,i0,a,a)') "XAS DEBUG atom type ", itype, " Z=", atoms%nz(itype), &
                                              " neq=", atoms%neq(itype), " species=", TRIM(atoms%speciesName(itype))
            WRITE(xas_debug_unit, '(a,i0,a,i0,a,i0,a,a)') "XAS DEBUG atom type ", itype, " Z=", atoms%nz(itype), &
                                                          " neq=", atoms%neq(itype), " species=", TRIM(atoms%speciesName(itype))
         END IF
         IF (atoms%nz(itype) == xas_test_atomic_number) THEN
            n_fe_types = n_fe_types + 1
            n_fe_atoms = n_fe_atoms + atoms%neq(itype)
         END IF
      END DO
      IF (n_fe_types == 0) THEN
         CALL juDFT_error("No Fe atom types found for hardwired XAS test", calledby="m_xas_driver")
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
      output_filename = "xas_test_"//TRIM(edge_name)//"_"//TRIM(pol_label)//"_eta"//TRIM(eta_label)//".dat"

      n_spin_channels = MERGE(1, input%jspins, noco%l_noco)
      IF (l_xas_debug_strength) THEN
         ALLOCATE(xas_debug_strength_spin(xas_debug_n_pol, n_spin_channels), SOURCE=0.0)
         IF (l_xas_debug_fullk_contrib_table .AND. xas_debug_verbosity >= 2) THEN
            ALLOCATE(xas_debug_strength_fullk(xas_debug_n_pol, kpts%nkptf), SOURCE=0.0)
            ALLOCATE(xas_debug_strength_fullk_l0(xas_debug_n_pol, kpts%nkptf), SOURCE=0.0)
            ALLOCATE(xas_debug_strength_fullk_l2(xas_debug_n_pol, kpts%nkptf), SOURCE=0.0)
         END IF
      END IF
      DO ikpt_i = 1, SIZE(fmpi%k_list)
         ikpt = fmpi%k_list(ikpt_i)
         IF (kpts%wtkpt(ikpt) <= 0.0) CYCLE
         weight_sum_parent = weight_sum_parent + kpts%wtkpt(ikpt)
         IF (l_xas_use_spatial_star) THEN
            CALL xas_count_star_members(kpts, ikpt, nstar)
            CALL xas_star_member_weight(kpts, ikpt, wk_star)
            weight_sum_star = weight_sum_star + REAL(nstar)*wk_star
         ELSE
            weight_sum_star = weight_sum_star + kpts%wtkpt(ikpt)
         END IF
      END DO
      IF (xas_debug_verbosity >= 1) THEN
         WRITE(*, '(a,i0,a,i0,a,l1,a,es12.4)') "XAS DEBUG star setup: nkpt=", kpts%nkpt, &
                                      " nkptf=", kpts%nkptf, " use_spatial_star=", l_xas_use_spatial_star, &
                                      " reduction_factor=", REAL(kpts%nkptf)/REAL(kpts%nkpt)
         WRITE(*, '(a,i0)') "XAS DEBUG star mode = ", xas_star_mode
         WRITE(*, '(a,es18.10,a,es18.10,a,es12.4)') "XAS DEBUG star weight sums: selected=", weight_sum_parent, &
                                            " expanded=", weight_sum_star, &
                                            " diff=", weight_sum_star - weight_sum_parent
         WRITE(xas_debug_unit, '(a,i0,a,i0,a,l1,a,es12.4)') "XAS DEBUG star setup: nkpt=", kpts%nkpt, &
                                      " nkptf=", kpts%nkptf, " use_spatial_star=", l_xas_use_spatial_star, &
                                      " reduction_factor=", REAL(kpts%nkptf)/REAL(kpts%nkpt)
         WRITE(xas_debug_unit, '(a,i0)') "XAS DEBUG star mode = ", xas_star_mode
         WRITE(xas_debug_unit, '(a,es18.10,a,es18.10,a,es12.4)') "XAS DEBUG star weight sums: selected=", weight_sum_parent, &
                                            " expanded=", weight_sum_star, &
                                            " diff=", weight_sum_star - weight_sum_parent
      END IF
      transition_min = HUGE(1.0)
      transition_max = -HUGE(1.0)
      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) /= xas_test_atomic_number) CYCLE
         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_extract_core_states(atoms, itype, edge_name, vTot%mt(1:atoms%jri(itype), 0, itype, 1), core_states)
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_extract_core_states", unit=xas_debug_unit)
         IF (SIZE(core_states) < 1) THEN
            CALL juDFT_error("No L3 core state found for selected Fe atom type in hardwired XAS test", calledby="m_xas_driver")
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
         IF (atoms%nz(itype) /= xas_test_atomic_number) CYCLE
         DO ispin = 1, input%jspins
            CALL genMTBasis(atoms, enpara, vTot, fmpi, itype, ispin, usdus, &
                            f(:, :, 0:, ispin), g(:, :, 0:, ispin), flo(:, :, :, ispin), l_writeArg=.FALSE.)
         END DO

         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_extract_core_states(atoms, itype, edge_name, vTot%mt(1:atoms%jri(itype), 0, itype, 1), core_states)
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_extract_core_states", unit=xas_debug_unit)
         IF (SIZE(core_states) < 1) THEN
            CALL juDFT_error("No L3 core state found for selected Fe atom type in hardwired XAS test", calledby="m_xas_driver")
         END IF

         CALL radfun%generate_radial_functions(atoms, input, enpara, fmpi, vTot, itype)
         max_order = MAXVAL(radfun%n_r(0:atoms%lmax(itype)))
         ALLOCATE(radial_xas(max_order, 0:atoms%lmaxd, input%jspins), SOURCE=0.0)
         CALL xas_debug_clear_underflow(l_xas_debug_fp)
         CALL xas_radial_dipole_integrals(atoms, itype, radfun, core_states(1)%p_core, radial_xas)
         CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_radial_dipole_integrals", unit=xas_debug_unit)

         lmax_xas = atoms%lmax(itype)
         DO jsp_loop = 1, n_spin_channels
            jsp = MERGE(1, jsp_loop, noco%l_noco)
            n_local_spins = MERGE(2, 1, noco%l_noco .OR. input%jspins == 2)
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

               IF (xas_star_mode == 2) THEN
                  CALL xas_count_star_members(kpts, ikpt, nstar)
                  CALL xas_star_member_weight(kpts, ikpt, wk_star)
                  ALLOCATE(matrix(nbands, SIZE(core_states(1)%twice_mj)))
                  IF (l_xas_debug_strength) THEN
                     ALLOCATE(matrix_debug(nbands, SIZE(core_states(1)%twice_mj)))
                     ALLOCATE(matrix_l0(nbands, SIZE(core_states(1)%twice_mj)))
                     ALLOCATE(matrix_l2(nbands, SIZE(core_states(1)%twice_mj)))
                  END IF

                  DO ikptf = 1, kpts%nkptf
                     IF (kpts%bkp(ikptf) /= ikpt) CYCLE
                     wk_current = wk_star

                     CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikptf, cell, fmpi)
                     nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
                     n_apw_rows = MERGE(lapw%nv(1) + lapw%nv(2), lapw%nv(1), noco%l_noco)
                     CALL zMat%init(l_real, nbasfcn, nbands)
                     CALL read_eig(eig_id, ikptf, jsp, list=ev_list, neig=nbands_read, zmat=zMat, &
                                   kpts=kpts, input=input, noco=noco, nococonv=nococonv, sym=sym, atoms=atoms, cell=cell)
                     IF (nbands_read < nbands) THEN
                        CALL juDFT_error("read_eig returned fewer bands than requested in hardwired XAS full-zone test", &
                                         calledby="m_xas_driver")
                     END IF
                     IF (l_xas_debug_mode2_norms .AND. xas_debug_verbosity >= 3) THEN
                        CALL xas_debug_zmat_stats(zMat, n_apw_rows, z_norm, z_apw_norm, z_lo_norm, z_max_abs)
                        WRITE(*, '(a,3i8,4es16.8)') "XAS DEBUG mode2 zmat ikptf parent bksym norm apw lo maxabs ", &
                           ikptf, ikpt, kpts%bksym(ikptf), z_norm, z_apw_norm, z_lo_norm, z_max_abs
                        WRITE(xas_debug_unit, '(a,3i8,4es16.8)') "XAS DEBUG mode2 zmat ikptf parent bksym norm apw lo maxabs ", &
                           ikptf, ikpt, kpts%bksym(ikptf), z_norm, z_apw_norm, z_lo_norm, z_max_abs
                     END IF

                     ALLOCATE(abc_spin(n_local_spins))
                     DO ispin = 1, n_local_spins
                        CALL abc_spin(ispin)%init(input, atoms, radfun%n_r, nbands, itype)
                     END DO
                     IF (noco%l_noco) THEN
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
                     IF (l_xas_debug_mode2_norms .AND. xas_debug_verbosity >= 3) THEN
                        DO ispin = 1, n_local_spins
                           CALL xas_debug_abc_stats(abc_spin(ispin), abc_max_abs, abc_sum_abs2)
                           WRITE(*, '(a,4i8,2es16.8)') "XAS DEBUG mode2 abc ikptf parent bksym ispin maxabs sumabs2 ", &
                              ikptf, ikpt, kpts%bksym(ikptf), ispin, abc_max_abs, abc_sum_abs2
                           WRITE(xas_debug_unit, '(a,4i8,2es16.8)') "XAS DEBUG mode2 abc ikptf parent bksym ispin maxabs sumabs2 ", &
                              ikptf, ikpt, kpts%bksym(ikptf), ispin, abc_max_abs, abc_sum_abs2
                        END DO
                     END IF

                     IF (l_xas_debug_strength) THEN
                        DO iatom_l = 1, atoms%neq(itype)
                           DO i_pol = 1, xas_debug_n_pol
                              eps_cart_debug = CMPLX(0.0, 0.0)
                              eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                              CALL xas_cartesian_to_spherical(eps_cart_debug, eps_sph_debug)

                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_debug)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l0, final_l=0)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l2, final_l=2)

                              debug_strength = xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                              debug_l0 = xas_debug_transition_strength(matrix_l0, occ_band, wk_current)
                              debug_l2 = xas_debug_transition_strength(matrix_l2, occ_band, wk_current)
                              xas_debug_strength_total(i_pol) = xas_debug_strength_total(i_pol) + debug_strength
                              xas_debug_strength_spin(i_pol, jsp_loop) = xas_debug_strength_spin(i_pol, jsp_loop) + debug_strength
                              IF (l_xas_debug_kpt_strength) THEN
                                 xas_debug_strength_kpt(i_pol, ikpt) = xas_debug_strength_kpt(i_pol, ikpt) + debug_strength
                              END IF
                              IF (l_xas_debug_fullk_contrib_table .AND. xas_debug_verbosity >= 2) THEN
                                 xas_debug_strength_fullk(i_pol, ikptf) = xas_debug_strength_fullk(i_pol, ikptf) + debug_strength
                                 xas_debug_strength_fullk_l0(i_pol, ikptf) = xas_debug_strength_fullk_l0(i_pol, ikptf) + debug_l0
                                 xas_debug_strength_fullk_l2(i_pol, ikptf) = xas_debug_strength_fullk_l2(i_pol, ikptf) + debug_l2
                              END IF
                              xas_debug_strength_l0(i_pol) = xas_debug_strength_l0(i_pol) + debug_l0
                              xas_debug_strength_l2(i_pol) = xas_debug_strength_l2(i_pol) + debug_l2
                           END DO
                        END DO
                     END IF

                     DO iatom_l = 1, atoms%neq(itype)
                        CALL xas_debug_clear_underflow(l_xas_debug_fp)
                        CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                          iatom_l, lmax_xas, matrix)
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

                     DEALLOCATE(abc_spin)
                  END DO

                  IF (l_xas_debug_strength) DEALLOCATE(matrix_debug, matrix_l0, matrix_l2)
                  DEALLOCATE(matrix, eig_band, occ_band, ev_list)
                  CYCLE
               END IF

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
               IF (noco%l_noco) THEN
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

               ALLOCATE(matrix(nbands, SIZE(core_states(1)%twice_mj)))
               IF (l_xas_debug_strength) THEN
                  ALLOCATE(matrix_debug(nbands, SIZE(core_states(1)%twice_mj)))
                  ALLOCATE(matrix_l0(nbands, SIZE(core_states(1)%twice_mj)))
                  ALLOCATE(matrix_l2(nbands, SIZE(core_states(1)%twice_mj)))
               END IF

               IF (l_xas_use_spatial_star) THEN
                  CALL xas_count_star_members(kpts, ikpt, nstar)
                  CALL xas_star_member_weight(kpts, ikpt, wk_star)
                  DO ikptf = 1, kpts%nkptf
                     IF (kpts%bkp(ikptf) /= ikpt) CYCLE
                     bksym = kpts%bksym(ikptf)
                     wk_current = wk_star

                     IF (xas_star_mode == 3) THEN
                        CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
                        ALLOCATE(abc_star_spin(n_local_spins))
                        DO ispin = 1, n_local_spins
                           CALL xas_rotate_abc_star_member(abc_spin(ispin), atoms, sym, cell, itype, bksym, lmax_xas, &
                                                           abc_star_spin(ispin))
                        END DO

                        IF (xas_debug_verbosity >= 3 .AND. bksym /= 1 .AND. n_sanity_print < 10) THEN
                           DO ispin = 1, n_local_spins
                              CALL xas_debug_abc_stats(abc_spin(ispin), abc_max_abs, abc_sum_abs2)
                              CALL xas_debug_abc_stats(abc_star_spin(ispin), z_max_abs, z_norm)
                              WRITE(*, '(a,5i8,4es16.8)') &
                                 "XAS DEBUG mode3 abc ikptf parent bksym spatial_iop ispin parent_max star_max parent_sum star_sum ", &
                                 ikptf, ikpt, bksym, spatial_iop, ispin, abc_max_abs, z_max_abs, abc_sum_abs2, z_norm
                              WRITE(xas_debug_unit, '(a,5i8,4es16.8)') &
                                 "XAS DEBUG mode3 abc ikptf parent bksym spatial_iop ispin parent_max star_max parent_sum star_sum ", &
                                 ikptf, ikpt, bksym, spatial_iop, ispin, abc_max_abs, z_max_abs, abc_sum_abs2, z_norm
                           END DO
                           n_sanity_print = n_sanity_print + 1
                        END IF

                        IF (l_xas_debug_strength) THEN
                           DO iatom_l = 1, atoms%neq(itype)
                              DO i_pol = 1, xas_debug_n_pol
                                 eps_cart_debug = CMPLX(0.0, 0.0)
                                 eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                                 CALL xas_cartesian_to_spherical(eps_cart_debug, eps_sph_debug)

                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_debug)
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_l0, final_l=0)
                                 CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                                   iatom_l, lmax_xas, matrix_l2, final_l=2)

                                 debug_strength = xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                                 debug_l0 = xas_debug_transition_strength(matrix_l0, occ_band, wk_current)
                                 debug_l2 = xas_debug_transition_strength(matrix_l2, occ_band, wk_current)
                                 xas_debug_strength_total(i_pol) = xas_debug_strength_total(i_pol) + debug_strength
                                 xas_debug_strength_spin(i_pol, jsp_loop) = xas_debug_strength_spin(i_pol, jsp_loop) + debug_strength
                                 IF (l_xas_debug_kpt_strength) THEN
                                    xas_debug_strength_kpt(i_pol, ikpt) = xas_debug_strength_kpt(i_pol, ikpt) + debug_strength
                                 END IF
                                 IF (l_xas_debug_fullk_contrib_table .AND. xas_debug_verbosity >= 2) THEN
                                    xas_debug_strength_fullk(i_pol, ikptf) = xas_debug_strength_fullk(i_pol, ikptf) + debug_strength
                                    xas_debug_strength_fullk_l0(i_pol, ikptf) = xas_debug_strength_fullk_l0(i_pol, ikptf) + debug_l0
                                    xas_debug_strength_fullk_l2(i_pol, ikptf) = xas_debug_strength_fullk_l2(i_pol, ikptf) + debug_l2
                                 END IF
                                 xas_debug_strength_l0(i_pol) = xas_debug_strength_l0(i_pol) + debug_l0
                                 xas_debug_strength_l2(i_pol) = xas_debug_strength_l2(i_pol) + debug_l2
                              END DO
                           END DO
                        END IF

                        DO iatom_l = 1, atoms%neq(itype)
                           CALL xas_debug_clear_underflow(l_xas_debug_fp)
                           CALL xas_core_band_matrixelements(abc_star_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                             iatom_l, lmax_xas, matrix)
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
                        CYCLE
                     END IF

                     CALL xas_rotate_lab_polarization_for_parent(sym, cell, bksym, eps_cart, eps_parent_cart)
                     CALL xas_cartesian_to_spherical(eps_parent_cart, eps_parent_sph)

                     IF (l_xas_star_sanity_print .AND. xas_debug_verbosity >= 2 .AND. n_sanity_print < 10) THEN
                        CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
                        IF (bksym <= sym%nop) THEN
                           rrot = TRANSPOSE(REAL(sym%mrot(:, :, sym%invtab(bksym))))
                        ELSE
                           rrot = -TRANSPOSE(REAL(sym%mrot(:, :, sym%invtab(bksym - sym%nop))))
                        END IF
                        reconstructed_k = kpts%to_first_bz(MATMUL(rrot, kpts%bk(:, ikpt)))
                        WRITE(*, '(a,i0,a,i0,a,i0,a,i0,a,l1)') "XAS DEBUG star member parent=", ikpt, &
                           " ikptf=", ikptf, " bksym=", bksym, " spatial_iop=", spatial_iop, &
                           " time_reversal=", l_time_reversal
                        WRITE(*, '(a,3f12.7,a,3f12.7,a,3f12.7)') "XAS DEBUG star k parent=", kpts%bk(:, ikpt), &
                           " full=", kpts%bkf(:, ikptf), " reconstructed=", reconstructed_k
                        WRITE(xas_debug_unit, '(a,i0,a,i0,a,i0,a,i0,a,l1)') "XAS DEBUG star member parent=", ikpt, &
                           " ikptf=", ikptf, " bksym=", bksym, " spatial_iop=", spatial_iop, &
                           " time_reversal=", l_time_reversal
                        WRITE(xas_debug_unit, '(a,3f12.7,a,3f12.7,a,3f12.7)') "XAS DEBUG star k parent=", kpts%bk(:, ikpt), &
                           " full=", kpts%bkf(:, ikptf), " reconstructed=", reconstructed_k
                        DO i_pol = 1, xas_debug_n_pol
                           eps_cart_debug = CMPLX(0.0, 0.0)
                           eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                           CALL xas_rotate_lab_polarization_for_parent(sym, cell, bksym, eps_cart_debug, eps_star_cart_debug)
                           WRITE(*, '(a,a,a,6f10.5)') "XAS DEBUG star eps_parent pol=", &
                              xas_debug_pol_label(i_pol), " re/im(x,y,z)=", &
                              REAL(eps_star_cart_debug(1)), AIMAG(eps_star_cart_debug(1)), &
                              REAL(eps_star_cart_debug(2)), AIMAG(eps_star_cart_debug(2)), &
                              REAL(eps_star_cart_debug(3)), AIMAG(eps_star_cart_debug(3))
                           WRITE(xas_debug_unit, '(a,a,a,6f10.5)') "XAS DEBUG star eps_parent pol=", &
                              xas_debug_pol_label(i_pol), " re/im(x,y,z)=", &
                              REAL(eps_star_cart_debug(1)), AIMAG(eps_star_cart_debug(1)), &
                              REAL(eps_star_cart_debug(2)), AIMAG(eps_star_cart_debug(2)), &
                              REAL(eps_star_cart_debug(3)), AIMAG(eps_star_cart_debug(3))
                        END DO
                        n_sanity_print = n_sanity_print + 1
                     END IF

                     IF (l_xas_debug_strength) THEN
                        DO iatom_l = 1, atoms%neq(itype)
                           DO i_pol = 1, xas_debug_n_pol
                              eps_cart_debug = CMPLX(0.0, 0.0)
                              eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                              CALL xas_rotate_lab_polarization_for_parent(sym, cell, bksym, eps_cart_debug, eps_star_cart_debug)
                              CALL xas_cartesian_to_spherical(eps_star_cart_debug, eps_star_sph_debug)

                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_star_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_debug)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_star_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l0, final_l=0)
                              CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_star_sph_debug, &
                                                                iatom_l, lmax_xas, matrix_l2, final_l=2)

                              debug_strength = xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                              debug_l0 = xas_debug_transition_strength(matrix_l0, occ_band, wk_current)
                              debug_l2 = xas_debug_transition_strength(matrix_l2, occ_band, wk_current)
                              xas_debug_strength_total(i_pol) = xas_debug_strength_total(i_pol) + debug_strength
                              xas_debug_strength_spin(i_pol, jsp_loop) = xas_debug_strength_spin(i_pol, jsp_loop) + debug_strength
                              IF (l_xas_debug_kpt_strength) THEN
                                 xas_debug_strength_kpt(i_pol, ikpt) = xas_debug_strength_kpt(i_pol, ikpt) + debug_strength
                              END IF
                              IF (l_xas_debug_fullk_contrib_table .AND. xas_debug_verbosity >= 2) THEN
                                 xas_debug_strength_fullk(i_pol, ikptf) = xas_debug_strength_fullk(i_pol, ikptf) + debug_strength
                                 xas_debug_strength_fullk_l0(i_pol, ikptf) = xas_debug_strength_fullk_l0(i_pol, ikptf) + debug_l0
                                 xas_debug_strength_fullk_l2(i_pol, ikptf) = xas_debug_strength_fullk_l2(i_pol, ikptf) + debug_l2
                              END IF
                              xas_debug_strength_l0(i_pol) = xas_debug_strength_l0(i_pol) + debug_l0
                              xas_debug_strength_l2(i_pol) = xas_debug_strength_l2(i_pol) + debug_l2

                              IF (l_xas_debug_compare_rotation_modes .AND. xas_star_mode == 1) THEN
                                 CALL xas_debug_rotate_lab_polarization_mode(sym, cell, bksym, eps_cart_debug, &
                                                                             eps_star_cart_debug, l_transpose=.FALSE.)
                                 CALL xas_cartesian_to_spherical(eps_star_cart_debug, eps_star_sph_debug)
                                 CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), &
                                                                   eps_star_sph_debug, iatom_l, lmax_xas, matrix_debug)
                                 xas_debug_strength_mode_b(i_pol) = xas_debug_strength_mode_b(i_pol) + &
                                                                    xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                              END IF
                           END DO
                        END DO
                     END IF

                     DO iatom_l = 1, atoms%neq(itype)
                        CALL xas_debug_clear_underflow(l_xas_debug_fp)
                        CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_parent_sph, &
                                                          iatom_l, lmax_xas, matrix)
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
                  END DO
               ELSE
                  wk_current = kpts%wtkpt(ikpt)
                  IF (l_xas_debug_strength) THEN
                     DO iatom_l = 1, atoms%neq(itype)
                        DO i_pol = 1, xas_debug_n_pol
                           eps_cart_debug = CMPLX(0.0, 0.0)
                           eps_cart_debug(i_pol) = CMPLX(1.0, 0.0)
                           CALL xas_cartesian_to_spherical(eps_cart_debug, eps_sph_debug)

                           CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                             iatom_l, lmax_xas, matrix_debug)
                           CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                             iatom_l, lmax_xas, matrix_l0, final_l=0)
                           CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph_debug, &
                                                             iatom_l, lmax_xas, matrix_l2, final_l=2)

                           debug_strength = xas_debug_transition_strength(matrix_debug, occ_band, wk_current)
                           debug_l0 = xas_debug_transition_strength(matrix_l0, occ_band, wk_current)
                           debug_l2 = xas_debug_transition_strength(matrix_l2, occ_band, wk_current)
                           xas_debug_strength_total(i_pol) = xas_debug_strength_total(i_pol) + debug_strength
                           xas_debug_strength_spin(i_pol, jsp_loop) = xas_debug_strength_spin(i_pol, jsp_loop) + debug_strength
                           IF (l_xas_debug_kpt_strength) THEN
                              xas_debug_strength_kpt(i_pol, ikpt) = xas_debug_strength_kpt(i_pol, ikpt) + debug_strength
                           END IF
                           IF (l_xas_debug_fullk_contrib_table .AND. xas_debug_verbosity >= 2) THEN
                              xas_debug_strength_fullk(i_pol, ikpt) = xas_debug_strength_fullk(i_pol, ikpt) + debug_strength
                              xas_debug_strength_fullk_l0(i_pol, ikpt) = xas_debug_strength_fullk_l0(i_pol, ikpt) + debug_l0
                              xas_debug_strength_fullk_l2(i_pol, ikpt) = xas_debug_strength_fullk_l2(i_pol, ikpt) + debug_l2
                           END IF
                           xas_debug_strength_l0(i_pol) = xas_debug_strength_l0(i_pol) + debug_l0
                           xas_debug_strength_l2(i_pol) = xas_debug_strength_l2(i_pol) + debug_l2
                        END DO
                     END DO
                  END IF

                  DO iatom_l = 1, atoms%neq(itype)
                     CALL xas_debug_clear_underflow(l_xas_debug_fp)
                     CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_sph, &
                                                       iatom_l, lmax_xas, matrix)
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

      IF (xas_debug_verbosity >= 1) THEN
         WRITE(*, '(a,i0,a,i0,a)') "XAS DEBUG: summed ", n_fe_atoms, &
                                   " total Fe atoms over ", n_fe_types, " atom types"
         WRITE(xas_debug_unit, '(a,i0,a,i0,a)') "XAS DEBUG: summed ", n_fe_atoms, &
                                                " total Fe atoms over ", n_fe_types, " atom types"
      END IF
      IF (l_xas_debug_kpt_strength) THEN
         WRITE(*, '(a)') "XAS DEBUG kpt strength columns: ikpt weight Sx Sy Sz rel_xz=(Sx-Sz)/avg"
         WRITE(xas_debug_unit, '(a)') "XAS DEBUG kpt strength columns: ikpt weight Sx Sy Sz rel_xz=(Sx-Sz)/avg"
         DO ikpt_i = 1, SIZE(fmpi%k_list)
            ikpt = fmpi%k_list(ikpt_i)
            debug_avg = SUM(xas_debug_strength_kpt(:, ikpt))/REAL(xas_debug_n_pol)
            debug_rel_xz = 0.0
            IF (ABS(debug_avg) > TINY(debug_avg)) THEN
               debug_rel_xz = (xas_debug_strength_kpt(1, ikpt) - xas_debug_strength_kpt(3, ikpt))/debug_avg
            END IF
            WRITE(*, '(a,i0,1x,5es16.8)') "XAS DEBUG kpt strength ", ikpt, kpts%wtkpt(ikpt), &
               xas_debug_strength_kpt(1, ikpt), xas_debug_strength_kpt(2, ikpt), &
               xas_debug_strength_kpt(3, ikpt), debug_rel_xz
            WRITE(xas_debug_unit, '(a,i0,1x,5es16.8)') "XAS DEBUG kpt strength ", ikpt, kpts%wtkpt(ikpt), &
               xas_debug_strength_kpt(1, ikpt), xas_debug_strength_kpt(2, ikpt), &
               xas_debug_strength_kpt(3, ikpt), debug_rel_xz
         END DO
      END IF
      IF (l_xas_debug_strength) THEN
         DO jsp_loop = 1, n_spin_channels
            CALL xas_debug_print_strength_summary("spin-resolved", jsp_loop, xas_debug_strength_spin(:, jsp_loop), xas_debug_unit)
         END DO
         DO i_pol = 1, xas_debug_n_pol
            xas_debug_strength_cross(i_pol) = xas_debug_strength_total(i_pol) - &
                                              xas_debug_strength_l0(i_pol) - xas_debug_strength_l2(i_pol)
            WRITE(*, '(a,a,a,es18.10,a,es18.10,a,es18.10,a,es18.10)') &
               "XAS DEBUG strength pol=", xas_debug_pol_label(i_pol), &
               " total=", xas_debug_strength_total(i_pol), &
               " l0=", xas_debug_strength_l0(i_pol), &
               " l2=", xas_debug_strength_l2(i_pol), &
               " cross_l0_l2=", xas_debug_strength_cross(i_pol)
            WRITE(xas_debug_unit, '(a,a,a,es18.10,a,es18.10,a,es18.10,a,es18.10)') &
               "XAS DEBUG strength pol=", xas_debug_pol_label(i_pol), &
               " total=", xas_debug_strength_total(i_pol), &
               " l0=", xas_debug_strength_l0(i_pol), &
               " l2=", xas_debug_strength_l2(i_pol), &
               " cross_l0_l2=", xas_debug_strength_cross(i_pol)
         END DO
         CALL xas_debug_print_strength_summary("total", 0, xas_debug_strength_total, xas_debug_unit)
         IF (l_xas_debug_compare_rotation_modes .AND. xas_star_mode == 1) THEN
            CALL xas_debug_print_strength_summary("rotation-mode-B-Rcart-eps", 0, xas_debug_strength_mode_b, xas_debug_unit)
         END IF
         IF (l_xas_debug_fullk_contrib_table .AND. xas_debug_verbosity >= 2) THEN
            CALL xas_debug_write_fullk_contrib_table(TRIM(xas_contrib_filename), kpts, sym, &
                                                     xas_debug_strength_fullk, xas_debug_strength_fullk_l0, &
                                                     xas_debug_strength_fullk_l2)
            WRITE(*, '(a,a)') "XAS hardwired test wrote k contribution table to ", TRIM(xas_contrib_filename)
            WRITE(xas_debug_unit, '(a,a)') "XAS hardwired test wrote k contribution table to ", TRIM(xas_contrib_filename)
         END IF
      END IF

      CALL xas_debug_clear_underflow(l_xas_debug_fp)
      CALL xas_write_spectrum_text(TRIM(output_filename), energy_grid, intensity, "Hartree")
      CALL xas_debug_report_underflow(l_xas_debug_fp, "xas_write_spectrum_text", unit=xas_debug_unit)
      WRITE(*, '(a,a)') "XAS hardwired test wrote spectrum to ", TRIM(output_filename)
      WRITE(xas_debug_unit, '(a,a)') "XAS hardwired test wrote spectrum to ", TRIM(output_filename)
      WRITE(*, '(a,a)') "XAS hardwired test wrote debug log to ", TRIM(xas_debug_filename)
      WRITE(xas_debug_unit, '(a,a)') "XAS hardwired test wrote debug log to ", TRIM(xas_debug_filename)
      IF (xas_debug_verbosity >= 3 .AND. n_underflow_spectrum > 0) THEN
         WRITE(*, '(a,i0,a)') "XAS DEBUG SUMMARY: underflow occurred in xas_accumulate_matrix_spectrum ", &
                              n_underflow_spectrum, &
                              " time(s); cleared as harmless Gaussian-tail/tiny-product underflow."
         WRITE(xas_debug_unit, '(a,i0,a)') "XAS DEBUG SUMMARY: underflow occurred in xas_accumulate_matrix_spectrum ", &
                                          n_underflow_spectrum, &
                                          " time(s); cleared as harmless Gaussian-tail/tiny-product underflow."
      END IF
      CLOSE(xas_debug_unit)
   END SUBROUTINE xas_hardwired_test_driver

   SUBROUTINE xas_debug_write_fullk_contrib_table(filename, kpts, sym, strength, strength_l0, strength_l2)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      TYPE(t_kpts),     INTENT(IN) :: kpts
      TYPE(t_sym),      INTENT(IN) :: sym
      REAL,             INTENT(IN) :: strength(:, :)
      REAL,             INTENT(IN) :: strength_l0(:, :)
      REAL,             INTENT(IN) :: strength_l2(:, :)

      INTEGER :: unit, ikptf, parent, bksym, spatial_iop
      LOGICAL :: l_time_reversal
      REAL    :: cross(xas_debug_n_pol)

      OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE")
      WRITE(unit, '(a)') "# XAS per-full-k pre-broadening contribution table"
      WRITE(unit, '(a)') "# columns:"
      WRITE(unit, '(a)') "# ikptf parent bksym spatial_iop time_reversal kx ky kz"
      WRITE(unit, '(a)') "# Sx Sy Sz l0x l0y l0z l2x l2y l2z crossx crossy crossz"
      DO ikptf = 1, kpts%nkptf
         parent = kpts%bkp(ikptf)
         bksym = kpts%bksym(ikptf)
         CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
         cross = strength(:, ikptf) - strength_l0(:, ikptf) - strength_l2(:, ikptf)
         WRITE(unit, '(4i8,1x,l1,15es18.10)') ikptf, parent, bksym, spatial_iop, l_time_reversal, &
            kpts%bkf(:, ikptf), strength(:, ikptf), strength_l0(:, ikptf), &
            strength_l2(:, ikptf), cross
      END DO
      CLOSE(unit)
   END SUBROUTINE xas_debug_write_fullk_contrib_table

   SUBROUTINE xas_debug_rotate_lab_polarization_mode(sym, cell, bksym, eps_lab_cart, eps_parent_cart, l_transpose)
      TYPE(t_sym),  INTENT(IN)  :: sym
      TYPE(t_cell), INTENT(IN)  :: cell
      INTEGER,      INTENT(IN)  :: bksym
      COMPLEX,      INTENT(IN)  :: eps_lab_cart(3)
      COMPLEX,      INTENT(OUT) :: eps_parent_cart(3)
      LOGICAL,      INTENT(IN)  :: l_transpose

      REAL    :: r_cart(3, 3)
      COMPLEX :: eps_source(3)
      INTEGER :: spatial_iop
      LOGICAL :: l_time_reversal

      CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
      CALL xas_cart_rotation_from_sym(sym, cell, spatial_iop, r_cart)

      eps_source = eps_lab_cart
      IF (l_time_reversal) eps_source = CONJG(eps_source)

      IF (l_transpose) THEN
         eps_parent_cart = MATMUL(CMPLX(TRANSPOSE(r_cart), 0.0), eps_source)
      ELSE
         eps_parent_cart = MATMUL(CMPLX(r_cart, 0.0), eps_source)
      END IF
   END SUBROUTINE xas_debug_rotate_lab_polarization_mode

   SUBROUTINE xas_debug_zmat_stats(zMat, n_apw_rows, norm_total, norm_apw, norm_lo, max_abs)
      TYPE(t_mat), INTENT(IN)  :: zMat
      INTEGER,     INTENT(IN)  :: n_apw_rows
      REAL,        INTENT(OUT) :: norm_total
      REAL,        INTENT(OUT) :: norm_apw
      REAL,        INTENT(OUT) :: norm_lo
      REAL,        INTENT(OUT) :: max_abs

      INTEGER :: n_apw

      n_apw = MIN(MAX(n_apw_rows, 0), zMat%matsize1)
      IF (zMat%l_real) THEN
         norm_total = SUM(zMat%data_r(1:zMat%matsize1, 1:zMat%matsize2)**2)
         max_abs = MAXVAL(ABS(zMat%data_r(1:zMat%matsize1, 1:zMat%matsize2)))
         IF (n_apw > 0) THEN
            norm_apw = SUM(zMat%data_r(1:n_apw, 1:zMat%matsize2)**2)
         ELSE
            norm_apw = 0.0
         END IF
         IF (n_apw < zMat%matsize1) THEN
            norm_lo = SUM(zMat%data_r(n_apw + 1:zMat%matsize1, 1:zMat%matsize2)**2)
         ELSE
            norm_lo = 0.0
         END IF
      ELSE
         norm_total = SUM(ABS(zMat%data_c(1:zMat%matsize1, 1:zMat%matsize2))**2)
         max_abs = MAXVAL(ABS(zMat%data_c(1:zMat%matsize1, 1:zMat%matsize2)))
         IF (n_apw > 0) THEN
            norm_apw = SUM(ABS(zMat%data_c(1:n_apw, 1:zMat%matsize2))**2)
         ELSE
            norm_apw = 0.0
         END IF
         IF (n_apw < zMat%matsize1) THEN
            norm_lo = SUM(ABS(zMat%data_c(n_apw + 1:zMat%matsize1, 1:zMat%matsize2))**2)
         ELSE
            norm_lo = 0.0
         END IF
      END IF
   END SUBROUTINE xas_debug_zmat_stats

   SUBROUTINE xas_debug_abc_stats(abc, max_abs, sum_abs2)
      TYPE(t_abc), INTENT(IN)  :: abc
      REAL,        INTENT(OUT) :: max_abs
      REAL,        INTENT(OUT) :: sum_abs2

      IF (.NOT. ALLOCATED(abc%cof)) THEN
         max_abs = 0.0
         sum_abs2 = 0.0
         RETURN
      END IF

      max_abs = MAXVAL(ABS(abc%cof))
      sum_abs2 = SUM(ABS(abc%cof)**2)
   END SUBROUTINE xas_debug_abc_stats

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

      IF (.NOT. l_debug) RETURN

      CALL IEEE_GET_FLAG(IEEE_UNDERFLOW, l_underflow)
      IF (l_underflow) THEN
         WRITE(*, '(a,a)') "XAS DEBUG: underflow occurred in ", TRIM(stage_name)
         IF (PRESENT(unit)) WRITE(unit, '(a,a)') "XAS DEBUG: underflow occurred in ", TRIM(stage_name)
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
