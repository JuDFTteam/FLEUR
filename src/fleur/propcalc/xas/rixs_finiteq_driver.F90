!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_finiteq_driver
   USE m_eig66_io, ONLY: read_eig
   USE m_genMTBasis, ONLY: genMTBasis
   USE m_juDFT, ONLY: juDFT_error
   USE m_mpi_reduce_tool, ONLY: mpi_sum_reduce
   USE m_rixs_finiteq_io, ONLY: rixs_finiteq_label, rixs_open_finiteq_pair_table, rixs_open_finiteq_site_table, &
                                rixs_write_finiteq_pair_rows, rixs_write_finiteq_site_rows
   USE m_rixs_io, ONLY: rixs_close_contribution_table, rixs_energy_label, rixs_print_contribution_check, &
                        rixs_print_pair_summary, rixs_write_spectrum_text
   USE m_rixs_momentum, ONLY: rixs_build_kq_map, rixs_fold_rlu, rixs_site_phase
   USE m_rixs_spectrum, ONLY: rixs_accumulate_finiteq_spinor_spectrum, &
                              rixs_add_finiteq_spinor_site_amplitudes, rixs_occupation_tolerance
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
   USE m_types_xas, ONLY: t_xas
   USE m_xas_angular, ONLY: xas_cartesian_to_spherical
   USE m_xas_core, ONLY: t_xas_core_state, xas_extract_core_states
   USE m_xas_matrixelements, ONLY: xas_band_core_emission_matrixelements, xas_core_band_matrixelements
   USE m_xas_radial, ONLY: xas_radial_dipole_integrals
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: n_pol = 3
   CHARACTER(LEN=1), PARAMETER :: pol_label(n_pol) = [CHARACTER(LEN=1) :: "x", "y", "z"]

   TYPE t_absorber_context
      INTEGER :: itype = 0
      TYPE(t_radfun) :: radfun
      TYPE(t_xas_core_state) :: core_state
      REAL, ALLOCATABLE :: radial_xas(:, :, :)
      COMPLEX :: spin_frame_transform(2, 2) = CMPLX(0.0, 0.0)
   END TYPE t_absorber_context

   PUBLIC :: rixs_run_finiteq_spinor

CONTAINS

   SUBROUTINE rixs_run_finiteq_spinor(eig_id, fmpi, input, rixs, kpts, atoms, sym, cell, noco, nococonv, enpara, vTot, results)
      INTEGER, INTENT(IN) :: eig_id
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_xas), INTENT(IN) :: rixs
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_potden), INTENT(IN) :: vTot
      TYPE(t_results), INTENT(IN) :: results

      TYPE(t_absorber_context), ALLOCATABLE :: contexts(:)
      TYPE(t_usdus) :: usdus
      TYPE(t_lapw) :: lapw_v, lapw_n
      TYPE(t_mat) :: zmat_v, zmat_n
      TYPE(t_abc), ALLOCATABLE :: abc_v(:), abc_n(:)
      COMPLEX, ALLOCATABLE :: matrix_abs(:, :), matrix_emit(:, :), amplitude_vn(:, :, :, :), site_partial_vn(:, :)
      COMPLEX :: eps_cart(3), eps_in_sph(-1:1), eps_out_sph(-1:1), site_phase
      REAL, ALLOCATABLE :: loss_grid(:), intensity(:, :, :), intensity_reduced(:, :, :)
      REAL, ALLOCATABLE :: contribution_intensity(:, :, :), contribution_intensity_reduced(:, :, :)
      REAL :: q_reduced_rlu(3)
      REAL, ALLOCATABLE :: f(:, :, :, :), g(:, :, :, :), flo(:, :, :, :)
      REAL, ALLOCATABLE :: eig_v(:), eig_n(:), occ_v(:), occ_n(:)
      INTEGER, ALLOCATABLE :: ev_list_v(:), ev_list_n(:), kq_index(:), reciprocal_shift(:, :)
      INTEGER :: pair_units(n_pol, n_pol), site_units(n_pol, n_pol)
      INTEGER :: ikpt_i, ikpt_v, ikpt_n, nbands_v, nbands_n, nbands_read, nbas_v, nbas_n
      INTEGER :: val_min, val_max, int_min, int_max, i_grid, i_band, i_context, iatom_l, iatom, ispin
      INTEGER :: i_pin, i_pout, n_contexts
      CHARACTER(LEN=4) :: rank_label
      CHARACTER(LEN=12) :: omega_label
      CHARACTER(LEN=60) :: q_full_label
      CHARACTER(LEN=240) :: filename
      LOGICAL :: l_root, l_kpt_group_root

      IF (.NOT. noco%l_noco .OR. input%jspins /= 2) THEN
         CALL juDFT_error("Finite-Q RIXS currently requires first-variation noco spinors with jspins=2.", &
                          calledby="m_rixs_finiteq_driver")
      END IF
      IF (rixs%rixs_write_state_character) THEN
         CALL juDFT_error("writeStateCharacter is not yet defined for distinct finite-Q valence/intermediate k points.", &
                          calledby="m_rixs_finiteq_driver")
      END IF

      l_root = fmpi%irank == 0
      l_kpt_group_root = fmpi%n_rank == 0
      pair_units = -1
      site_units = -1
      q_reduced_rlu = rixs_fold_rlu(rixs%rixs_q_full_rlu)
      CALL rixs_build_kq_map(kpts, rixs%rixs_q_full_rlu, kq_index, reciprocal_shift)

      ALLOCATE(loss_grid(rixs%rixs_n_loss))
      DO i_grid = 1, rixs%rixs_n_loss
         loss_grid(i_grid) = rixs%rixs_loss_min + (rixs%rixs_loss_max - rixs%rixs_loss_min) &
                             *REAL(i_grid - 1)/REAL(rixs%rixs_n_loss - 1)
      END DO
      ALLOCATE(intensity(rixs%rixs_n_loss, n_pol, n_pol), SOURCE=0.0)
      IF (rixs%rixs_write_contributions) THEN
         ALLOCATE(contribution_intensity(rixs%rixs_n_loss, n_pol, n_pol), SOURCE=0.0)
      END IF

      omega_label = rixs_energy_label(rixs%rixs_omega_in)
      q_full_label = rixs_finiteq_label(rixs%rixs_q_full_rlu)
      IF (rixs%rixs_write_contributions) THEN
         WRITE(rank_label, '(i4.4)') fmpi%irank
         DO i_pin = 1, n_pol
            IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
            DO i_pout = 1, n_pol
               IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
               filename = TRIM(rixs%rixs_output_prefix)//"_"//TRIM(rixs%rixs_edge)//"_"//pol_label(i_pin)//"_"// &
                          pol_label(i_pout)//"_"//TRIM(q_full_label)//"_omega"//TRIM(omega_label)// &
                          "_finiteq_pair_rank"//rank_label//".dat"
               CALL rixs_open_finiteq_pair_table(TRIM(filename), rixs%rixs_edge, rixs%rixs_absorber_z, pol_label(i_pin), &
                  pol_label(i_pout), rixs%rixs_omega_in, rixs%rixs_q_full_rlu, q_reduced_rlu, &
                  fmpi%irank, pair_units(i_pin, i_pout))
               filename = TRIM(rixs%rixs_output_prefix)//"_"//TRIM(rixs%rixs_edge)//"_"//pol_label(i_pin)//"_"// &
                          pol_label(i_pout)//"_"//TRIM(q_full_label)//"_omega"//TRIM(omega_label)// &
                          "_finiteq_site_rank"//rank_label//".dat"
               CALL rixs_open_finiteq_site_table(TRIM(filename), rixs%rixs_edge, rixs%rixs_absorber_z, pol_label(i_pin), &
                  pol_label(i_pout), rixs%rixs_omega_in, rixs%rixs_q_full_rlu, q_reduced_rlu, &
                  fmpi%irank, site_units(i_pin, i_pout))
            END DO
         END DO
      END IF

      IF (l_kpt_group_root) THEN
         CALL prepare_absorbers(contexts, usdus, f, g, flo, fmpi, input, rixs, atoms, enpara, vTot, nococonv)
         n_contexts = SIZE(contexts)
         DO ikpt_i = 1, SIZE(fmpi%k_list)
            ikpt_v = fmpi%k_list(ikpt_i)
            ikpt_n = kq_index(ikpt_v)
            IF (kpts%wtkpt(ikpt_v) <= 0.0) CYCLE
            nbands_v = results%neig(ikpt_v, 1)
            nbands_n = results%neig(ikpt_n, 1)
            IF (nbands_v <= 0 .OR. nbands_n <= 0) CYCLE
            CALL finiteq_band_bounds(rixs, nbands_v, nbands_n, val_min, val_max, int_min, int_max)
            IF (val_min > val_max .OR. int_min > int_max) CYCLE

            ALLOCATE(eig_v(nbands_v), occ_v(nbands_v), eig_n(nbands_n), occ_n(nbands_n))
            eig_v = results%eig(1:nbands_v, ikpt_v, 1)
            eig_n = results%eig(1:nbands_n, ikpt_n, 1)
            occ_v = results%w_iks(1:nbands_v, ikpt_v, 1)/kpts%wtkpt(ikpt_v)
            occ_n = results%w_iks(1:nbands_n, ikpt_n, 1)/kpts%wtkpt(ikpt_n)
            IF (.NOT. ANY(occ_v(val_min:val_max) > rixs_occupation_tolerance) .OR. &
                .NOT. ANY(1.0 - occ_n(int_min:int_max) > rixs_occupation_tolerance)) THEN
               DEALLOCATE(eig_v, occ_v, eig_n, occ_n)
               CYCLE
            END IF

            ALLOCATE(ev_list_v(nbands_v), ev_list_n(nbands_n))
            ev_list_v = [(i_band, i_band=1, nbands_v)]
            ev_list_n = [(i_band, i_band=1, nbands_n)]
            CALL lapw_v%init(input, noco, nococonv, kpts, atoms, sym, ikpt_v, cell, fmpi)
            CALL lapw_n%init(input, noco, nococonv, kpts, atoms, sym, ikpt_n, cell, fmpi)
            nbas_v = lapw_v%nv(1) + lapw_v%nv(2) + 2*atoms%nlotot
            nbas_n = lapw_n%nv(1) + lapw_n%nv(2) + 2*atoms%nlotot
            CALL zmat_v%init(.FALSE., nbas_v, nbands_v)
            CALL zmat_n%init(.FALSE., nbas_n, nbands_n)
            CALL read_eig(eig_id, ikpt_v, 1, list=ev_list_v, neig=nbands_read, zmat=zmat_v)
            IF (nbands_read < nbands_v) CALL juDFT_error("Finite-Q RIXS valence eigenvector read is incomplete.", &
                                                         calledby="m_rixs_finiteq_driver")
            CALL read_eig(eig_id, ikpt_n, 1, list=ev_list_n, neig=nbands_read, zmat=zmat_n)
            IF (nbands_read < nbands_n) CALL juDFT_error("Finite-Q RIXS intermediate eigenvector read is incomplete.", &
                                                         calledby="m_rixs_finiteq_driver")
            ALLOCATE(amplitude_vn(nbands_v, nbands_n, n_pol, n_pol), SOURCE=CMPLX(0.0, 0.0))
            ALLOCATE(site_partial_vn(nbands_v, nbands_n))

            DO i_context = 1, n_contexts
               ALLOCATE(abc_v(2), abc_n(2))
               DO ispin = 1, 2
                  CALL abc_v(ispin)%init(input, atoms, nbands_v, contexts(i_context)%itype)
                  CALL abc_v(ispin)%calc_abc(input, atoms, sym, cell, lapw_v, nbands_v, usdus, noco, nococonv, ispin, &
                                             contexts(i_context)%itype, zmat_v)
                  CALL abc_n(ispin)%init(input, atoms, nbands_n, contexts(i_context)%itype)
                  CALL abc_n(ispin)%calc_abc(input, atoms, sym, cell, lapw_n, nbands_n, usdus, noco, nococonv, ispin, &
                                             contexts(i_context)%itype, zmat_n)
               END DO
               ALLOCATE(matrix_abs(nbands_n, SIZE(contexts(i_context)%core_state%twice_mj)))
               ALLOCATE(matrix_emit(nbands_v, SIZE(contexts(i_context)%core_state%twice_mj)))
               DO iatom_l = 1, atoms%neq(contexts(i_context)%itype)
                  iatom = atoms%firstAtom(contexts(i_context)%itype) + iatom_l - 1
                  site_phase = rixs_site_phase(rixs%rixs_q_full_rlu, atoms%taual(:, iatom))
                  DO i_pin = 1, n_pol
                     IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
                     eps_cart = CMPLX(0.0, 0.0)
                     eps_cart(i_pin) = CMPLX(1.0, 0.0)
                     CALL xas_cartesian_to_spherical(eps_cart, eps_in_sph)
                     CALL xas_core_band_matrixelements(abc_n, contexts(i_context)%radfun, contexts(i_context)%radial_xas, &
                        contexts(i_context)%core_state, eps_in_sph, iatom_l, atoms%lmax(contexts(i_context)%itype), &
                        matrix_abs, spin_frame_transform=contexts(i_context)%spin_frame_transform)
                     DO i_pout = 1, n_pol
                        IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
                        eps_cart = CMPLX(0.0, 0.0)
                        eps_cart(i_pout) = CMPLX(1.0, 0.0)
                        CALL xas_cartesian_to_spherical(eps_cart, eps_out_sph)
                        ! Pass epsilon_out itself.  The emission helper performs
                        ! the Hermitian conjugation exactly once.
                        CALL xas_band_core_emission_matrixelements(abc_v, contexts(i_context)%radfun, &
                           contexts(i_context)%radial_xas, contexts(i_context)%core_state, eps_out_sph, iatom_l, &
                           atoms%lmax(contexts(i_context)%itype), matrix_emit, &
                           spin_frame_transform=contexts(i_context)%spin_frame_transform)
                        CALL rixs_add_finiteq_spinor_site_amplitudes(eig_n, contexts(i_context)%core_state%energy, &
                           rixs%rixs_omega_in, rixs%rixs_gamma_core, matrix_abs, matrix_emit, site_phase, val_min, val_max, &
                           int_min, int_max, amplitude_vn(:, :, i_pin, i_pout), site_partial_vn)
                        IF (rixs%rixs_write_contributions) THEN
                           CALL rixs_write_finiteq_site_rows(site_units(i_pin, i_pout), ikpt_v, ikpt_n, iatom, &
                              contexts(i_context)%itype, atoms%taual(:, iatom), site_phase, eig_n, occ_v, occ_n, &
                              contexts(i_context)%core_state%energy, rixs%rixs_omega_in, rixs%rixs_gamma_core, &
                              site_partial_vn, val_min, val_max, int_min, int_max)
                        END IF
                     END DO
                  END DO
               END DO
               DEALLOCATE(matrix_abs, matrix_emit, abc_v, abc_n)
            END DO

            DO i_pin = 1, n_pol
               IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
               DO i_pout = 1, n_pol
                  IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
                  CALL rixs_accumulate_finiteq_spinor_spectrum(loss_grid, eig_v, occ_v, eig_n, occ_n, &
                     kpts%wtkpt(ikpt_v), rixs%rixs_eta_loss, amplitude_vn(:, :, i_pin, i_pout), val_min, val_max, &
                     int_min, int_max, intensity(:, i_pin, i_pout))
                  IF (rixs%rixs_write_contributions) THEN
                     CALL rixs_write_finiteq_pair_rows(pair_units(i_pin, i_pout), ikpt_v, ikpt_n, kpts%bk(:, ikpt_v), &
                        kpts%bk(:, ikpt_n), reciprocal_shift(:, ikpt_v), eig_v, occ_v, eig_n, occ_n, kpts%wtkpt(ikpt_v), &
                        amplitude_vn(:, :, i_pin, i_pout), loss_grid, rixs%rixs_eta_loss, val_min, val_max, int_min, int_max, &
                        contribution_intensity(:, i_pin, i_pout))
                  END IF
               END DO
            END DO
            DEALLOCATE(amplitude_vn, site_partial_vn, ev_list_v, ev_list_n, eig_v, occ_v, eig_n, occ_n)
         END DO
      END IF

      IF (rixs%rixs_write_contributions) THEN
         DO i_pin = 1, n_pol
            DO i_pout = 1, n_pol
               CALL rixs_close_contribution_table(pair_units(i_pin, i_pout))
               CALL rixs_close_contribution_table(site_units(i_pin, i_pout))
            END DO
         END DO
      END IF
      ALLOCATE(intensity_reduced(rixs%rixs_n_loss, n_pol, n_pol), SOURCE=0.0)
      CALL mpi_sum_reduce(intensity, intensity_reduced, fmpi%mpi_comm)
      IF (rixs%rixs_write_contributions) THEN
         ALLOCATE(contribution_intensity_reduced(rixs%rixs_n_loss, n_pol, n_pol), SOURCE=0.0)
         CALL mpi_sum_reduce(contribution_intensity, contribution_intensity_reduced, fmpi%mpi_comm)
      END IF

      IF (l_root) THEN
         WRITE(*, '(a,3f12.6)') " Finite-Q RIXS full Q (FLEUR RLU): ", rixs%rixs_q_full_rlu
         WRITE(*, '(a,3f12.6)') " Finite-Q RIXS reduced q          : ", q_reduced_rlu
         WRITE(*, '(a)') " Finite-Q absorber-site sum       : coherent complex amplitudes before squaring"
         DO i_pin = 1, n_pol
            IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
            DO i_pout = 1, n_pol
               IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
               CALL rixs_print_pair_summary(pol_label(i_pin), pol_label(i_pout), loss_grid, &
                                            intensity_reduced(:, i_pin, i_pout))
               IF (rixs%rixs_write_contributions) THEN
                  CALL rixs_print_contribution_check(pol_label(i_pin), pol_label(i_pout), &
                     intensity_reduced(:, i_pin, i_pout), contribution_intensity_reduced(:, i_pin, i_pout), l_spinor=.TRUE.)
               END IF
               filename = TRIM(rixs%rixs_output_prefix)//"_"//TRIM(rixs%rixs_edge)//"_"//pol_label(i_pin)//"_"// &
                          pol_label(i_pout)//"_"//TRIM(q_full_label)//"_omega"//TRIM(omega_label)//"_finiteq_rixs.dat"
               CALL rixs_write_spectrum_text(TRIM(filename), loss_grid, intensity_reduced(:, i_pin, i_pout))
               WRITE(*, '(a,a)') "RIXS wrote finite-Q spectrum to ", TRIM(filename)
            END DO
         END DO
         IF (rixs%rixs_in_polarizations(3) .AND. rixs%rixs_out_polarizations(1) .AND. &
             rixs%rixs_out_polarizations(2)) THEN
            WRITE(*, '(a)') " Unanalysed outgoing signal: form I_zx + I_zy from the separate spectra after squaring."
         END IF
      END IF
   END SUBROUTINE rixs_run_finiteq_spinor

   SUBROUTINE prepare_absorbers(contexts, usdus, f, g, flo, fmpi, input, rixs, atoms, enpara, vTot, nococonv)
      TYPE(t_absorber_context), ALLOCATABLE, INTENT(OUT) :: contexts(:)
      TYPE(t_usdus), INTENT(INOUT) :: usdus
      REAL, ALLOCATABLE, INTENT(OUT) :: f(:, :, :, :), g(:, :, :, :), flo(:, :, :, :)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_xas), INTENT(IN) :: rixs
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_potden), INTENT(IN) :: vTot
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_xas_core_state), ALLOCATABLE :: core_states(:)
      INTEGER :: itype, ispin, n_contexts, i_context, max_order

      n_contexts = COUNT(atoms%nz == rixs%rixs_absorber_z)
      IF (n_contexts < 1) CALL juDFT_error("Finite-Q RIXS found no matching absorber type.", &
                                           calledby="m_rixs_finiteq_driver")
      ALLOCATE(contexts(n_contexts))
      CALL usdus%init(atoms, input%jspins)
      ALLOCATE(f(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
      ALLOCATE(g(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
      ALLOCATE(flo(atoms%jmtd, 2, atoms%nlod, input%jspins))
      i_context = 0
      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) /= rixs%rixs_absorber_z) CYCLE
         i_context = i_context + 1
         contexts(i_context)%itype = itype
         DO ispin = 1, input%jspins
            CALL genMTBasis(atoms, enpara, vTot, fmpi, itype, ispin, usdus, f(:, :, 0:, ispin), g(:, :, 0:, ispin), &
                            flo(:, :, :, ispin), l_writeArg=.FALSE.)
         END DO
         CALL xas_extract_core_states(atoms, itype, rixs%rixs_edge, vTot%mt(1:atoms%jri(itype), 0, itype, 1), core_states)
         IF (SIZE(core_states) < 1) CALL juDFT_error("Finite-Q RIXS could not extract the requested core edge.", &
                                                    calledby="m_rixs_finiteq_driver")
         contexts(i_context)%core_state = core_states(1)
         CALL contexts(i_context)%radfun%generate_radial_functions(atoms, input, enpara, fmpi, vTot, itype)
         max_order = MAXVAL(contexts(i_context)%radfun%n_r(0:atoms%lmax(itype)))
         ALLOCATE(contexts(i_context)%radial_xas(max_order, 0:atoms%lmaxd, input%jspins), SOURCE=0.0)
         CALL xas_radial_dipole_integrals(atoms, itype, contexts(i_context)%radfun, &
                                          contexts(i_context)%core_state%p_core, contexts(i_context)%radial_xas)
         contexts(i_context)%spin_frame_transform = CONJG(TRANSPOSE(nococonv%umat(itype)))
         DEALLOCATE(core_states)
      END DO
   END SUBROUTINE prepare_absorbers

   SUBROUTINE finiteq_band_bounds(rixs, nbands_v, nbands_n, val_min, val_max, int_min, int_max)
      TYPE(t_xas), INTENT(IN) :: rixs
      INTEGER, INTENT(IN) :: nbands_v, nbands_n
      INTEGER, INTENT(OUT) :: val_min, val_max, int_min, int_max

      val_min = 1
      IF (rixs%l_rixs_valence_band_min) val_min = rixs%rixs_valence_band_min
      val_max = nbands_v
      IF (rixs%l_rixs_valence_band_max) val_max = MIN(rixs%rixs_valence_band_max, nbands_v)
      int_min = 1
      IF (rixs%l_rixs_intermediate_band_min) int_min = rixs%rixs_intermediate_band_min
      int_max = nbands_n
      IF (rixs%l_rixs_intermediate_band_max) int_max = MIN(rixs%rixs_intermediate_band_max, nbands_n)
   END SUBROUTINE finiteq_band_bounds

END MODULE m_rixs_finiteq_driver
