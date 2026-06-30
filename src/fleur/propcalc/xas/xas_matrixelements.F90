!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_matrixelements
   USE m_juDFT, ONLY: juDFT_error
   USE m_types_abc, ONLY: t_abc
   USE m_types_radfun, ONLY: t_radfun
   USE m_xas_angular, ONLY: xas_dipole_angular_coeff, xas_sigma_down, xas_sigma_up
   USE m_xas_core, ONLY: t_xas_core_state
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: xas_core_band_matrixelements
   PUBLIC :: xas_print_largest_matrixelement_partials

CONTAINS

   SUBROUTINE xas_core_band_matrixelements(abc_spin, radfun, radial_xas, core_state, eps_sph, iAtom_l, lmax, matrix, &
                                           final_l, spin_frame_transform)
      ! Computes M(band,mj_index) using FLEUR lm indexing lm=l*(l+1)+m.
      ! abc_spin(1) is the local MT spin-up component and maps to sigma=+1;
      ! abc_spin(2) is the local MT spin-down component and maps to sigma=-1.
      ! In noco calculations abc_spin(:) is in the local MT spin frame. The
      ! core spin-angular coefficients are constructed in the global spin
      ! basis and must be rotated into the same local frame with the optional
      ! spin_frame_transform(tau,s) before contraction.
      ! If final_l is present, only that final-state angular momentum channel
      ! is included. The driver uses this for compact channel summaries.
      TYPE(t_abc),            INTENT(IN)  :: abc_spin(:)
      TYPE(t_radfun),         INTENT(IN)  :: radfun
      REAL,                   INTENT(IN)  :: radial_xas(:, 0:, :)
      TYPE(t_xas_core_state), INTENT(IN)  :: core_state
      COMPLEX,                INTENT(IN)  :: eps_sph(-1:1)
      INTEGER,                INTENT(IN)  :: iAtom_l, lmax
      COMPLEX,                INTENT(OUT) :: matrix(:, :)
      INTEGER, OPTIONAL,      INTENT(IN)  :: final_l
      COMPLEX, OPTIONAL,      INTENT(IN)  :: spin_frame_transform(:, :)

      INTEGER :: band, i_mj, ispin, sigma, l, m, lm, iOrd, nbands, s_global
      COMPLEX :: angular, radial_sum
      COMPLEX :: angular_global(2), angular_local(2)
      LOGICAL :: l_rotate_spin_frame

      CALL xas_check_matrixelement_inputs(abc_spin, radfun, radial_xas, core_state, iAtom_l, lmax, matrix)
      IF (PRESENT(final_l)) THEN
         IF (final_l < 0 .OR. final_l > lmax) THEN
            CALL juDFT_error("Invalid final_l in XAS matrix elements", calledby="m_xas_matrixelements")
         END IF
      END IF
      l_rotate_spin_frame = PRESENT(spin_frame_transform)
      IF (l_rotate_spin_frame) THEN
         IF (SIZE(abc_spin) /= 2) THEN
            CALL juDFT_error("XAS spin-frame transform requires two abc spin components", calledby="m_xas_matrixelements")
         END IF
         IF (SIZE(spin_frame_transform, 1) < 2 .OR. SIZE(spin_frame_transform, 2) < 2) THEN
            CALL juDFT_error("XAS spin-frame transform must be at least 2x2", calledby="m_xas_matrixelements")
         END IF
      END IF

      nbands = SIZE(abc_spin(1)%cof, 1)
      matrix = CMPLX(0.0, 0.0)

      DO i_mj = 1, SIZE(core_state%twice_mj)
         DO l = 0, lmax
            IF (PRESENT(final_l)) THEN
               IF (l /= final_l) CYCLE
            END IF
            DO m = -l, l
               lm = l*(l + 1) + m
               IF (l_rotate_spin_frame) THEN
                  DO s_global = 1, 2
                     sigma = xas_sigma_from_spin_index(s_global)
                     angular_global(s_global) = xas_dipole_angular_coeff(l, m, core_state%lc, core_state%twice_j, &
                                                                          core_state%twice_mj(i_mj), sigma, eps_sph)
                  END DO
                  DO ispin = 1, 2
                     angular_local(ispin) = spin_frame_transform(ispin, 1)*angular_global(1) + &
                                            spin_frame_transform(ispin, 2)*angular_global(2)
                  END DO
               ELSE
                  DO ispin = 1, SIZE(abc_spin)
                     sigma = xas_sigma_from_spin_index(ispin)
                     angular_local(ispin) = xas_dipole_angular_coeff(l, m, core_state%lc, core_state%twice_j, &
                                                                     core_state%twice_mj(i_mj), sigma, eps_sph)
                  END DO
               END IF
               DO ispin = 1, SIZE(abc_spin)
                  angular = angular_local(ispin)
                  IF (ABS(angular) == 0.0) CYCLE
                  DO band = 1, nbands
                     radial_sum = CMPLX(0.0, 0.0)
                     DO iOrd = 1, radfun%n_r(l)
                        radial_sum = radial_sum + CONJG(abc_spin(ispin)%cof(band, lm, iOrd, iAtom_l)) &
                                                   * radial_xas(iOrd, l, ispin)
                     END DO
                     matrix(band, i_mj) = matrix(band, i_mj) + angular*radial_sum
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE xas_core_band_matrixelements

   SUBROUTINE xas_print_largest_matrixelement_partials(abc_spin, radfun, radial_xas, core_state, eps_sph, iAtom_l, lmax, &
                                                       n_print, out_unit)
      TYPE(t_abc),            INTENT(IN)           :: abc_spin(:)
      TYPE(t_radfun),         INTENT(IN)           :: radfun
      REAL,                   INTENT(IN)           :: radial_xas(:, 0:, :)
      TYPE(t_xas_core_state), INTENT(IN)           :: core_state
      COMPLEX,                INTENT(IN)           :: eps_sph(-1:1)
      INTEGER,                INTENT(IN)           :: iAtom_l, lmax
      INTEGER,                INTENT(IN), OPTIONAL :: n_print
      INTEGER,                INTENT(IN), OPTIONAL :: out_unit

      INTEGER :: nout, unit, nbands, total_partials, printed
      INTEGER :: band, i_mj, ispin, sigma, l, m, lm, iOrd
      INTEGER :: best_band, best_mj, best_spin, best_l, best_m, best_iOrd
      LOGICAL, ALLOCATABLE :: used(:, :, :, :, :, :)
      REAL    :: best_abs, abs_partial
      COMPLEX :: angular, partial

      CALL xas_check_matrixelement_inputs(abc_spin, radfun, radial_xas, core_state, iAtom_l, lmax)

      nout = 20
      IF (PRESENT(n_print)) nout = n_print
      unit = 6
      IF (PRESENT(out_unit)) unit = out_unit
      nbands = SIZE(abc_spin(1)%cof, 1)
      total_partials = nbands*SIZE(core_state%twice_mj)*SIZE(abc_spin)*(lmax + 1)**2*SIZE(radial_xas, 1)
      ALLOCATE(used(nbands, SIZE(core_state%twice_mj), SIZE(abc_spin), 0:lmax, -lmax:lmax, SIZE(radial_xas, 1)), SOURCE=.FALSE.)

      WRITE (unit, '(a)') "Largest XAS matrix-element partial contributions:"
      DO printed = 1, MIN(nout, total_partials)
         best_abs = -1.0
         best_band = 0
         best_mj = 0
         best_spin = 0
         best_l = 0
         best_m = 0
         best_iOrd = 0

         DO i_mj = 1, SIZE(core_state%twice_mj)
            DO ispin = 1, SIZE(abc_spin)
               sigma = xas_sigma_from_spin_index(ispin)
               DO l = 0, lmax
                  DO m = -l, l
                     lm = l*(l + 1) + m
                     angular = xas_dipole_angular_coeff(l, m, core_state%lc, core_state%twice_j, &
                                                         core_state%twice_mj(i_mj), sigma, eps_sph)
                     IF (ABS(angular) == 0.0) CYCLE
                     DO band = 1, nbands
                        DO iOrd = 1, radfun%n_r(l)
                           IF (used(band, i_mj, ispin, l, m, iOrd)) CYCLE
                           partial = angular*CONJG(abc_spin(ispin)%cof(band, lm, iOrd, iAtom_l))*radial_xas(iOrd, l, ispin)
                           abs_partial = ABS(partial)
                           IF (abs_partial > best_abs) THEN
                              best_abs = abs_partial
                              best_band = band
                              best_mj = i_mj
                              best_spin = ispin
                              best_l = l
                              best_m = m
                              best_iOrd = iOrd
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         IF (best_abs < 0.0) EXIT
         used(best_band, best_mj, best_spin, best_l, best_m, best_iOrd) = .TRUE.
         sigma = xas_sigma_from_spin_index(best_spin)
         angular = xas_dipole_angular_coeff(best_l, best_m, core_state%lc, core_state%twice_j, &
                                             core_state%twice_mj(best_mj), sigma, eps_sph)
         partial = angular*CONJG(abc_spin(best_spin)%cof(best_band, best_l*(best_l + 1) + best_m, best_iOrd, iAtom_l)) &
                   * radial_xas(best_iOrd, best_l, best_spin)
         WRITE (unit, '(a,i5,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,2es16.8,a,es16.8)') &
            "  band=", best_band, " mj_index=", best_mj, " twice_mj=", core_state%twice_mj(best_mj), &
            " spin=", best_spin, " l=", best_l, " m=", best_m, " iOrd=", best_iOrd, &
            " value=", REAL(partial), AIMAG(partial), " abs=", best_abs
      END DO
   END SUBROUTINE xas_print_largest_matrixelement_partials

   SUBROUTINE xas_check_matrixelement_inputs(abc_spin, radfun, radial_xas, core_state, iAtom_l, lmax, matrix)
      TYPE(t_abc),            INTENT(IN)           :: abc_spin(:)
      TYPE(t_radfun),         INTENT(IN)           :: radfun
      REAL,                   INTENT(IN)           :: radial_xas(:, 0:, :)
      TYPE(t_xas_core_state), INTENT(IN)           :: core_state
      INTEGER,                INTENT(IN)           :: iAtom_l, lmax
      COMPLEX,                INTENT(IN), OPTIONAL :: matrix(:, :)

      INTEGER :: ispin, l, lm_max_needed, nbands

      IF (SIZE(abc_spin) < 1 .OR. SIZE(abc_spin) > 2) THEN
         CALL juDFT_error("abc_spin must contain one or two local spin components in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF
      IF (.NOT. ALLOCATED(radfun%n_r)) THEN
         CALL juDFT_error("radfun%n_r is not allocated in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF
      IF (.NOT. ALLOCATED(core_state%twice_mj)) THEN
         CALL juDFT_error("core_state%twice_mj is not allocated in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF
      IF (lmax < 0 .OR. lmax > UBOUND(radfun%n_r, 1)) THEN
         CALL juDFT_error("Invalid lmax in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF
      IF (UBOUND(radial_xas, 2) < lmax) THEN
         CALL juDFT_error("radial_xas l dimension is too small in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF
      IF (SIZE(radial_xas, 3) < SIZE(abc_spin)) THEN
         CALL juDFT_error("radial_xas spin dimension is too small in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF
      IF (SIZE(radial_xas, 1) < MAXVAL(radfun%n_r(0:lmax))) THEN
         CALL juDFT_error("radial_xas radial-order dimension is too small in XAS matrix elements", calledby="m_xas_matrixelements")
      END IF

      DO ispin = 1, SIZE(abc_spin)
         IF (.NOT. ALLOCATED(abc_spin(ispin)%cof)) THEN
            CALL juDFT_error("abc_spin%cof is not allocated in XAS matrix elements", calledby="m_xas_matrixelements")
         END IF
      END DO
      nbands = SIZE(abc_spin(1)%cof, 1)
      DO ispin = 1, SIZE(abc_spin)
         IF (SIZE(abc_spin(ispin)%cof, 1) /= nbands) THEN
            CALL juDFT_error("abc_spin band dimensions differ in XAS matrix elements", calledby="m_xas_matrixelements")
         END IF
         IF (iAtom_l < LBOUND(abc_spin(ispin)%cof, 4) .OR. iAtom_l > UBOUND(abc_spin(ispin)%cof, 4)) THEN
            CALL juDFT_error("iAtom_l outside abc_spin atom dimension in XAS matrix elements", calledby="m_xas_matrixelements")
         END IF
         lm_max_needed = lmax*(lmax + 2)
         IF (UBOUND(abc_spin(ispin)%cof, 2) < lm_max_needed) THEN
            CALL juDFT_error("lmax exceeds abc_spin lm dimension in XAS matrix elements", calledby="m_xas_matrixelements")
         END IF
         DO l = 0, lmax
            IF (SIZE(abc_spin(ispin)%cof, 3) < radfun%n_r(l)) THEN
               CALL juDFT_error("abc_spin radial-order dimension is too small in XAS matrix elements", calledby="m_xas_matrixelements")
            END IF
         END DO
      END DO
      IF (PRESENT(matrix)) THEN
         IF (SIZE(matrix, 1) < nbands .OR. SIZE(matrix, 2) < SIZE(core_state%twice_mj)) THEN
            CALL juDFT_error("matrix output dimensions are too small in XAS matrix elements", calledby="m_xas_matrixelements")
         END IF
      END IF
   END SUBROUTINE xas_check_matrixelement_inputs

   INTEGER FUNCTION xas_sigma_from_spin_index(ispin) RESULT(sigma)
      INTEGER, INTENT(IN) :: ispin

      SELECT CASE (ispin)
      CASE (1)
         sigma = xas_sigma_up
      CASE (2)
         sigma = xas_sigma_down
      CASE DEFAULT
         CALL juDFT_error("Unsupported local spin index in XAS matrix elements", calledby="m_xas_matrixelements")
      END SELECT
   END FUNCTION xas_sigma_from_spin_index

END MODULE m_xas_matrixelements
