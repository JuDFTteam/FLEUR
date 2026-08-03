!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_driver
   USE m_constants, ONLY: hartree_to_ev_const
   USE m_eig66_io, ONLY: read_eig
   USE m_genMTBasis, ONLY: genMTBasis
   USE m_juDFT, ONLY: juDFT_error
   USE m_mpi_reduce_tool, ONLY: mpi_sum_reduce
   USE m_rixs_io, ONLY: rixs_close_contribution_table, rixs_energy_label, rixs_open_contribution_table, &
                        rixs_open_spinor_contribution_table, &
                        rixs_print_contribution_check, rixs_print_pair_summary, rixs_print_setup_summary, &
                        rixs_write_contribution_rows, rixs_write_spinor_contribution_rows, rixs_write_spectrum_text
   USE m_rixs_spectrum, ONLY: rixs_accumulate_scalar_spin_trace_spectrum, rixs_accumulate_spinor_spectrum
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
   USE m_xas_matrixelements, ONLY: xas_band_core_emission_matrixelements, xas_band_core_emission_matrixelements_one_spin, &
                                  xas_core_band_matrixelements, xas_core_band_matrixelements_one_spin
   USE m_xas_radial, ONLY: xas_radial_dipole_integrals
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: rixs_n_pol = 3
   REAL, PARAMETER :: rixs_occ_tol = 1.0e-10
   CHARACTER(LEN=1), PARAMETER :: rixs_pol_label(rixs_n_pol) = [CHARACTER(LEN=1) :: "x", "y", "z"]

   PUBLIC :: rixs_run_driver

CONTAINS

   SUBROUTINE rixs_run_driver(eig_id, fmpi, input, rixs, kpts, atoms, sym, cell, noco, nococonv, enpara, vTot, results)
      INTEGER,             INTENT(IN) :: eig_id
      TYPE(t_mpi),         INTENT(IN) :: fmpi
      TYPE(t_input),       INTENT(IN) :: input
      TYPE(t_xas),         INTENT(IN) :: rixs
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
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      TYPE(t_xas_core_state), ALLOCATABLE :: core_states(:)
      TYPE(t_abc), ALLOCATABLE :: abc_spin(:)

      COMPLEX :: eps_cart(3), eps_in_sph(-1:1), eps_out_sph(-1:1)
      COMPLEX :: spin_frame_transform(2, 2)
      COMPLEX, ALLOCATABLE :: matrix_abs(:, :), matrix_emit(:, :)
      COMPLEX, ALLOCATABLE :: matrix_abs_spin(:, :, :), matrix_emit_spin(:, :, :)
      REAL, ALLOCATABLE :: loss_grid(:), intensity(:, :, :), intensity_reduced(:, :, :)
      REAL, ALLOCATABLE :: contribution_intensity(:, :, :), contribution_intensity_reduced(:, :, :)
      REAL, ALLOCATABLE :: radial_xas(:, :, :)
      REAL, ALLOCATABLE :: f(:, :, :, :), g(:, :, :, :), flo(:, :, :, :)
      REAL, ALLOCATABLE :: eig_band(:), occ_band(:)
      INTEGER, ALLOCATABLE :: ev_list(:)

      CHARACTER(LEN=4) :: rank_label
      CHARACTER(LEN=12) :: omega_label
      CHARACTER(LEN=160) :: output_filename, contribution_filename
      CHARACTER(LEN=240) :: error_message
      INTEGER :: ikpt_i, ikpt, jsp, nbands, nbands_read, itype, ispin, max_order, nbasfcn, lmax_rixs
      INTEGER :: i_grid, i_band, iatom_l, i_pin, i_pout, n_absorber_types, n_absorber_atoms
      INTEGER :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max
      INTEGER :: contribution_units(rixs_n_pol, rixs_n_pol)
      LOGICAL :: l_root, l_kpt_group_root, l_real, l_spinor_rixs

      IF (.NOT. rixs%l_rixs) RETURN

      l_root = fmpi%irank == 0
      l_kpt_group_root = fmpi%n_rank == 0
      contribution_units = -1
      CALL rixs_check_supported_input(input, rixs, kpts, noco)
      l_spinor_rixs = noco%l_noco
      IF (.NOT. ALLOCATED(results%w_iks)) THEN
         CALL juDFT_error("results%w_iks is not allocated in rixs_run_driver", calledby="m_rixs_driver")
      END IF

      IF (l_root) CALL rixs_print_setup_summary(rixs, noco%l_noco, noco%l_soc)
      omega_label = rixs_energy_label(rixs%rixs_omega_in)
      IF (rixs%rixs_write_contributions) THEN
         WRITE(rank_label, '(i4.4)') fmpi%irank
         DO i_pin = 1, rixs_n_pol
            IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
            DO i_pout = 1, rixs_n_pol
               IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
               contribution_filename = TRIM(rixs%rixs_output_prefix)//"_"//TRIM(rixs%rixs_edge)//"_"// &
                                       rixs_pol_label(i_pin)//"_"//rixs_pol_label(i_pout)// &
                                       "_omega"//TRIM(omega_label)//"_contrib_rank"//rank_label//".dat"
               IF (l_spinor_rixs) THEN
                  CALL rixs_open_spinor_contribution_table(TRIM(contribution_filename), rixs%rixs_edge, &
                     rixs%rixs_absorber_z, rixs_pol_label(i_pin), rixs_pol_label(i_pout), rixs%rixs_omega_in, &
                     fmpi%irank, contribution_units(i_pin, i_pout))
               ELSE
                  CALL rixs_open_contribution_table(TRIM(contribution_filename), rixs%rixs_edge, rixs%rixs_absorber_z, &
                                                    rixs_pol_label(i_pin), rixs_pol_label(i_pout), rixs%rixs_omega_in, &
                                                    fmpi%irank, contribution_units(i_pin, i_pout))
               END IF
            END DO
         END DO
      END IF

      n_absorber_types = 0
      n_absorber_atoms = 0
      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) == rixs%rixs_absorber_z) THEN
            n_absorber_types = n_absorber_types + 1
            n_absorber_atoms = n_absorber_atoms + atoms%neq(itype)
         END IF
      END DO
      IF (n_absorber_types == 0) THEN
         WRITE(error_message, '(a,i0)') "No atom types found for requested RIXS absorber Z=", rixs%rixs_absorber_z
         CALL juDFT_error(TRIM(error_message), calledby="m_rixs_driver")
      END IF

      ALLOCATE(loss_grid(rixs%rixs_n_loss))
      DO i_grid = 1, rixs%rixs_n_loss
         loss_grid(i_grid) = rixs%rixs_loss_min + (rixs%rixs_loss_max - rixs%rixs_loss_min) &
                             *REAL(i_grid - 1)/REAL(rixs%rixs_n_loss - 1)
      END DO
      ALLOCATE(intensity(rixs%rixs_n_loss, rixs_n_pol, rixs_n_pol), SOURCE=0.0)
      IF (rixs%rixs_write_contributions) THEN
         ALLOCATE(contribution_intensity(rixs%rixs_n_loss, rixs_n_pol, rixs_n_pol), SOURCE=0.0)
      END IF

      ! fmpi%k_list is shared by all ranks in a k-point subgroup. RIXS currently
      ! performs serial work for each k-point, so only subgroup roots calculate it.
      IF (l_kpt_group_root) THEN
      CALL usdus%init(atoms, input%jspins)
      ALLOCATE(f(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
      ALLOCATE(g(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
      ALLOCATE(flo(atoms%jmtd, 2, atoms%nlod, input%jspins))

      jsp = 1
      l_real = sym%invs .AND. (.NOT. noco%l_soc) .AND. (.NOT. noco%l_noco) .AND. atoms%n_hia == 0
      DO itype = 1, atoms%ntype
         IF (atoms%nz(itype) /= rixs%rixs_absorber_z) CYCLE
         IF (l_spinor_rixs) THEN
            DO ispin = 1, input%jspins
               CALL genMTBasis(atoms, enpara, vTot, fmpi, itype, ispin, usdus, &
                               f(:, :, 0:, ispin), g(:, :, 0:, ispin), flo(:, :, :, ispin), l_writeArg=.FALSE.)
            END DO
         ELSE
            CALL genMTBasis(atoms, enpara, vTot, fmpi, itype, jsp, usdus, &
                            f(:, :, 0:, jsp), g(:, :, 0:, jsp), flo(:, :, :, jsp), l_writeArg=.FALSE.)
         END IF
         CALL xas_extract_core_states(atoms, itype, rixs%rixs_edge, vTot%mt(1:atoms%jri(itype), 0, itype, 1), core_states)
         IF (SIZE(core_states) < 1) THEN
            WRITE(error_message, '(a,a,a,i0,a,i0)') "No core state found for requested RIXS edge ", TRIM(rixs%rixs_edge), &
                                                    " in absorber Z=", rixs%rixs_absorber_z, " atom type ", itype
            CALL juDFT_error(TRIM(error_message), calledby="m_rixs_driver")
         END IF

         CALL radfun%generate_radial_functions(atoms, input, enpara, fmpi, vTot, itype)
         max_order = MAXVAL(radfun%n_r(0:atoms%lmax(itype)))
         ALLOCATE(radial_xas(max_order, 0:atoms%lmaxd, input%jspins), SOURCE=0.0)
         CALL xas_radial_dipole_integrals(atoms, itype, radfun, core_states(1)%p_core, radial_xas)
         lmax_rixs = atoms%lmax(itype)
         IF (l_spinor_rixs) THEN
            ! abc%calc_abc returns the two first-variation components in the
            ! local MT spin frame. Rotate the global core spin-angular
            ! coefficients with the identical validated XAS U^\dagger
            ! convention before contracting them with these local components.
            spin_frame_transform = CONJG(TRANSPOSE(nococonv%umat(itype)))
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
            CALL rixs_effective_band_bounds(rixs, nbands, valence_band_min, valence_band_max, &
                                            intermediate_band_min, intermediate_band_max)
            IF (valence_band_min > valence_band_max .OR. intermediate_band_min > intermediate_band_max) THEN
               DEALLOCATE(ev_list, eig_band, occ_band)
               CYCLE
            END IF
            IF (.NOT. rixs_kpoint_has_valence_and_empty(occ_band, valence_band_min, valence_band_max, &
                                                        intermediate_band_min, intermediate_band_max)) THEN
               DEALLOCATE(ev_list, eig_band, occ_band)
               CYCLE
            END IF

            CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell, fmpi)
            nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, &
                            lapw%nv(1) + atoms%nlotot, l_spinor_rixs)
            CALL zMat%init(l_real, nbasfcn, nbands)
            CALL read_eig(eig_id, ikpt, jsp, list=ev_list, neig=nbands_read, zmat=zMat)
            IF (nbands_read < nbands) THEN
               CALL juDFT_error("read_eig returned fewer bands than requested in RIXS", calledby="m_rixs_driver")
            END IF

            IF (l_spinor_rixs) THEN
               ALLOCATE(abc_spin(2))
               DO ispin = 1, 2
                  CALL abc_spin(ispin)%init(input, atoms, nbands, itype)
                  CALL abc_spin(ispin)%calc_abc(input, atoms, sym, cell, lapw, nbands, usdus, noco, nococonv, &
                                                ispin, itype, zMat)
               END DO
               ALLOCATE(matrix_abs(nbands, SIZE(core_states(1)%twice_mj)))
               ALLOCATE(matrix_emit(nbands, SIZE(core_states(1)%twice_mj)))
            ELSE
               ALLOCATE(abc_spin(1))
               CALL abc_spin(1)%init(input, atoms, nbands, itype)
               CALL abc_spin(1)%calc_abc(input, atoms, sym, cell, lapw, nbands, usdus, noco, nococonv, jsp, itype, zMat)
               ALLOCATE(matrix_abs_spin(nbands, SIZE(core_states(1)%twice_mj), 2))
               ALLOCATE(matrix_emit_spin(nbands, SIZE(core_states(1)%twice_mj), 2))
            END IF

            DO iatom_l = 1, atoms%neq(itype)
               DO i_pin = 1, rixs_n_pol
                  IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
                  eps_cart = CMPLX(0.0, 0.0)
                  eps_cart(i_pin) = CMPLX(1.0, 0.0)
                  CALL xas_cartesian_to_spherical(eps_cart, eps_in_sph)
                  IF (l_spinor_rixs) THEN
                     CALL xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_states(1), eps_in_sph, &
                                                       iatom_l, lmax_rixs, matrix_abs, &
                                                       spin_frame_transform=spin_frame_transform)
                  ELSE
                     DO ispin = 1, 2
                        CALL xas_core_band_matrixelements_one_spin(abc_spin(1), radfun, radial_xas(:, 0:, 1), &
                           core_states(1), eps_in_sph, iatom_l, lmax_rixs, ispin, matrix_abs_spin(:, :, ispin))
                     END DO
                  END IF

                  DO i_pout = 1, rixs_n_pol
                     IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
                     eps_cart = CMPLX(0.0, 0.0)
                     eps_cart(i_pout) = CMPLX(1.0, 0.0)
                     CALL xas_cartesian_to_spherical(eps_cart, eps_out_sph)
                     IF (l_spinor_rixs) THEN
                        CALL xas_band_core_emission_matrixelements(abc_spin, radfun, radial_xas, core_states(1), &
                           eps_out_sph, iatom_l, lmax_rixs, matrix_emit, spin_frame_transform=spin_frame_transform)
                        CALL rixs_accumulate_spinor_spectrum(loss_grid, eig_band, occ_band, kpts%wtkpt(ikpt), &
                           core_states(1)%energy, rixs%rixs_omega_in, rixs%rixs_gamma_core, rixs%rixs_eta_loss, &
                           matrix_abs, matrix_emit, valence_band_min, valence_band_max, intermediate_band_min, &
                           intermediate_band_max, intensity(:, i_pin, i_pout))
                     ELSE
                        DO ispin = 1, 2
                           CALL xas_band_core_emission_matrixelements_one_spin(abc_spin(1), radfun, radial_xas(:, 0:, 1), &
                              core_states(1), eps_out_sph, iatom_l, lmax_rixs, ispin, matrix_emit_spin(:, :, ispin))
                        END DO
                        CALL rixs_accumulate_scalar_spin_trace_spectrum(loss_grid, eig_band, occ_band, kpts%wtkpt(ikpt), &
                           core_states(1)%energy, rixs%rixs_omega_in, rixs%rixs_gamma_core, rixs%rixs_eta_loss, &
                           matrix_abs_spin, matrix_emit_spin, valence_band_min, valence_band_max, intermediate_band_min, &
                           intermediate_band_max, intensity(:, i_pin, i_pout))
                     END IF
                     IF (rixs%rixs_write_contributions) THEN
                        IF (l_spinor_rixs) THEN
                           CALL rixs_write_spinor_contribution_rows(contribution_units(i_pin, i_pout), ikpt, &
                              atoms%firstAtom(itype) + iatom_l - 1, itype, eig_band, occ_band, kpts%wtkpt(ikpt), &
                              core_states(1)%energy, rixs%rixs_omega_in, rixs%rixs_gamma_core, matrix_abs, matrix_emit, &
                              loss_grid, rixs%rixs_eta_loss, valence_band_min, valence_band_max, intermediate_band_min, &
                              intermediate_band_max, contribution_intensity(:, i_pin, i_pout))
                        ELSE
                           CALL rixs_write_contribution_rows(contribution_units(i_pin, i_pout), ikpt, &
                              atoms%firstAtom(itype) + iatom_l - 1, itype, eig_band, occ_band, kpts%wtkpt(ikpt), &
                              core_states(1)%energy, rixs%rixs_omega_in, rixs%rixs_gamma_core, matrix_abs_spin, &
                              matrix_emit_spin, loss_grid, rixs%rixs_eta_loss, valence_band_min, valence_band_max, &
                              intermediate_band_min, intermediate_band_max, contribution_intensity(:, i_pin, i_pout))
                        END IF
                     END IF
                  END DO
               END DO
            END DO
            IF (l_spinor_rixs) THEN
               DEALLOCATE(matrix_abs, matrix_emit)
            ELSE
               DEALLOCATE(matrix_abs_spin, matrix_emit_spin)
            END IF
            DEALLOCATE(abc_spin, ev_list, eig_band, occ_band)
         END DO
         DEALLOCATE(radial_xas, core_states)
      END DO
      END IF

      IF (rixs%rixs_write_contributions) THEN
         DO i_pin = 1, rixs_n_pol
            DO i_pout = 1, rixs_n_pol
               CALL rixs_close_contribution_table(contribution_units(i_pin, i_pout))
            END DO
         END DO
      END IF

      ALLOCATE(intensity_reduced(SIZE(intensity, 1), SIZE(intensity, 2), SIZE(intensity, 3)), SOURCE=0.0)
      CALL mpi_sum_reduce(intensity, intensity_reduced, fmpi%mpi_comm)
      IF (rixs%rixs_write_contributions) THEN
         ALLOCATE(contribution_intensity_reduced(SIZE(contribution_intensity, 1), SIZE(contribution_intensity, 2), &
                                                 SIZE(contribution_intensity, 3)), SOURCE=0.0)
         CALL mpi_sum_reduce(contribution_intensity, contribution_intensity_reduced, fmpi%mpi_comm)
      END IF

      IF (l_root) THEN
         WRITE(*, '(a,i0,a,i0,a)') "RIXS summed ", n_absorber_atoms, " absorber atoms over ", &
                                   n_absorber_types, " atom types."
         DO i_pin = 1, rixs_n_pol
            IF (.NOT. rixs%rixs_in_polarizations(i_pin)) CYCLE
            DO i_pout = 1, rixs_n_pol
               IF (.NOT. rixs%rixs_out_polarizations(i_pout)) CYCLE
               output_filename = TRIM(rixs%rixs_output_prefix)//"_"//TRIM(rixs%rixs_edge)//"_"// &
                                 rixs_pol_label(i_pin)//"_"//rixs_pol_label(i_pout)// &
                                 "_omega"//TRIM(omega_label)//"_rixs.dat"
               CALL rixs_print_pair_summary(rixs_pol_label(i_pin), rixs_pol_label(i_pout), loss_grid, &
                                            intensity_reduced(:, i_pin, i_pout))
               IF (rixs%rixs_write_contributions) THEN
                  CALL rixs_print_contribution_check(rixs_pol_label(i_pin), rixs_pol_label(i_pout), &
                                                     intensity_reduced(:, i_pin, i_pout), &
                                                     contribution_intensity_reduced(:, i_pin, i_pout), &
                                                     l_spinor=l_spinor_rixs)
               END IF
               CALL rixs_write_spectrum_text(TRIM(output_filename), loss_grid, intensity_reduced(:, i_pin, i_pout))
               WRITE(*, '(a,a)') "RIXS wrote spectrum to ", TRIM(output_filename)
            END DO
         END DO
      END IF
   END SUBROUTINE rixs_run_driver

   SUBROUTINE rixs_check_supported_input(input, rixs, kpts, noco)
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_xas),   INTENT(IN) :: rixs
      TYPE(t_kpts),  INTENT(IN) :: kpts
      TYPE(t_noco),  INTENT(IN) :: noco

      IF (noco%l_soc .AND. .NOT. noco%l_noco) THEN
         CALL juDFT_error("Second-variation SOC RIXS is unsupported because the validated first-variation two-component "// &
                          "local spinor abc representation is unavailable.", calledby="m_rixs_driver")
      END IF
      IF (.NOT. noco%l_noco .AND. input%jspins /= 1) THEN
         CALL juDFT_error("Scalar RIXS requires input%jspins=1; collinear spin-polarized RIXS is not supported.", &
                          calledby="m_rixs_driver")
      END IF
      IF (noco%l_noco .AND. input%jspins /= 2) THEN
         CALL juDFT_error("First-variation spinor RIXS requires input%jspins=2 to construct both local spinor abc components.", &
                          calledby="m_rixs_driver")
      END IF
      IF (kpts%nkpt /= kpts%nkptf) THEN
         ! Keep RIXS on explicit full-k points: spinor time-reversal star
         ! reconstruction is not implemented, and current XAS star helpers omit
         ! non-symmorphic translation phases that coherent RIXS amplitudes need.
         CALL juDFT_error("RIXS requires a full-k/no-star k-point set (nkpt == nkptf); symmetry-star RIXS is unsupported.", &
                          calledby="m_rixs_driver")
      END IF
      SELECT CASE (TRIM(ADJUSTL(rixs%rixs_edge)))
      CASE ("K", "k", "1s1/2", "1S1/2", "L2", "l2", "2p1/2", "2P1/2", "L3", "l3", "2p3/2", "2P3/2")
      CASE DEFAULT
         CALL juDFT_error("Unsupported RIXS edge. Supported first-prototype edges are K, L2, and L3.", calledby="m_rixs_driver")
      END SELECT
   END SUBROUTINE rixs_check_supported_input

   SUBROUTINE rixs_effective_band_bounds(rixs, nbands, valence_band_min, valence_band_max, &
                                         intermediate_band_min, intermediate_band_max)
      TYPE(t_xas), INTENT(IN) :: rixs
      INTEGER,     INTENT(IN) :: nbands
      INTEGER,     INTENT(OUT) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max

      valence_band_min = 1
      IF (rixs%l_rixs_valence_band_min) valence_band_min = rixs%rixs_valence_band_min
      valence_band_max = nbands
      IF (rixs%l_rixs_valence_band_max) valence_band_max = MIN(rixs%rixs_valence_band_max, nbands)

      intermediate_band_min = 1
      IF (rixs%l_rixs_intermediate_band_min) intermediate_band_min = rixs%rixs_intermediate_band_min
      intermediate_band_max = nbands
      IF (rixs%l_rixs_intermediate_band_max) intermediate_band_max = MIN(rixs%rixs_intermediate_band_max, nbands)
   END SUBROUTINE rixs_effective_band_bounds

   LOGICAL FUNCTION rixs_kpoint_has_valence_and_empty(occ_band, valence_band_min, valence_band_max, &
                                                     intermediate_band_min, intermediate_band_max) RESULT(ok)
      REAL, INTENT(IN) :: occ_band(:)
      INTEGER, INTENT(IN) :: valence_band_min, valence_band_max, intermediate_band_min, intermediate_band_max

      ok = ANY(occ_band(valence_band_min:valence_band_max) > rixs_occ_tol) .AND. &
           ANY(1.0 - occ_band(intermediate_band_min:intermediate_band_max) > rixs_occ_tol)
   END FUNCTION rixs_kpoint_has_valence_and_empty

END MODULE m_rixs_driver
