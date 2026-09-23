!--------------------------------------------------------------------------------
! Projector-only validation for m_local_d_spin_projector_core.
!--------------------------------------------------------------------------------

PROGRAM local_d_spin_projector_test
   USE m_intgr, ONLY: intgr0
   USE m_local_d_spin_projector_core, ONLY: local_d_l, local_d_n_spins, &
                                            local_d_spin_contract_density, &
                                            local_d_spin_ordered_radial_metric
   IMPLICIT NONE

   INTEGER, PARAMETER :: n_grid = 61, n_radial = 3, n_state = 10
   REAL, PARAMETER :: r0 = 1.0e-4, dx = LOG(2.8/r0)/REAL(n_grid - 1)
   REAL, PARAMETER :: tolerance = 5.0e-12
   REAL :: radial(n_grid, 2, n_radial, local_d_n_spins)
   REAL :: radial_identical(n_grid, 2, n_radial, local_d_n_spins)
   REAL :: metric(n_radial, n_radial, local_d_n_spins, local_d_n_spins)
   REAL :: metric_identical(n_radial, n_radial, local_d_n_spins, local_d_n_spins)
   COMPLEX :: coefficients(-local_d_l:local_d_l, local_d_n_spins, n_radial)
   COMPLEX :: coefficients_phase(-local_d_l:local_d_l, local_d_n_spins, n_radial)
   COMPLEX :: rho(-local_d_l:local_d_l, local_d_n_spins, -local_d_l:local_d_l, local_d_n_spins)
   COMPLEX :: rho_reference(-local_d_l:local_d_l, local_d_n_spins, -local_d_l:local_d_l, local_d_n_spins)
   COMPLEX :: rho_phase(-local_d_l:local_d_l, local_d_n_spins, -local_d_l:local_d_l, local_d_n_spins)
   REAL :: direct_metric_12, direct_metric_21, norm, trace_rho, trace_direct
   REAL :: hermiticity_error, phase_error, direct_error, cross_spin_error
   REAL :: eigenvalues(n_state), min_eigenvalue, sigma_y
   REAL :: max_abs_error, max_relative_error, rho_scale
   COMPLEX :: phase
   INTEGER :: failures, i, m, ispin

   failures = 0
   CALL initialize_spin_dependent_radial_functions(radial)
   CALL local_d_spin_ordered_radial_metric(radial, r0, dx, n_grid, metric)

   CALL direct_radial_overlap(radial, 1, 2, 1, 2, direct_metric_12)
   CALL direct_radial_overlap(radial, 2, 1, 1, 2, direct_metric_21)
   CALL check_real("ordered S12(1,2)", metric(1, 2, 1, 2), direct_metric_12, tolerance, failures)
   CALL check_real("ordered S12(2,1)", metric(2, 1, 1, 2), direct_metric_21, tolerance, failures)
   IF (ABS(metric(1, 2, 1, 2) - metric(2, 1, 1, 2)) < 1.0e-5) THEN
      WRITE (*, '(a)') "FAIL: critical cross-spin radial functions did not produce distinct ordered overlaps"
      failures = failures + 1
   END IF
   cross_spin_error = MAXVAL(ABS(metric(:, :, 1, 2) - TRANSPOSE(metric(:, :, 2, 1))))
   CALL check_bound("S12-transpose(S21)", cross_spin_error, tolerance, failures)
   IF (MIN(ABS(metric(1, 2, 1, 1)), ABS(metric(1, 3, 1, 1)), ABS(metric(2, 3, 1, 1))) < 1.0e-5) THEN
      WRITE (*, '(a)') "FAIL: synthetic u/udot/LO basis does not have the required nonzero cross terms"
      failures = failures + 1
   ELSE
      WRITE (*, '(a)') "PASS: synthetic u/udot/LO basis has nonzero cross terms"
   END IF

   DO m = -local_d_l, local_d_l
      DO ispin = 1, local_d_n_spins
         DO i = 1, n_radial
            coefficients(m, ispin, i) = CMPLX(0.07*REAL(m + 3) + 0.03*REAL(ispin) + 0.02*REAL(i), &
                                               -0.04*REAL(m) + 0.05*REAL(ispin) - 0.01*REAL(i))
         END DO
      END DO
   END DO
   CALL local_d_spin_contract_density(coefficients, metric, rho)
   CALL direct_density(coefficients, radial, rho_reference)
   rho_scale = MAXVAL(ABS(rho_reference))
   direct_error = MAXVAL(ABS(rho - rho_reference))
   CALL check_bound("complete nonorthogonal u/udot/LO contraction", direct_error, tolerance, failures)

   hermiticity_error = local_hermiticity_error(rho)
   CALL check_bound("rho Hermiticity", hermiticity_error, tolerance, failures)
   CALL hermitian_eigenvalues(rho, eigenvalues)
   min_eigenvalue = MINVAL(eigenvalues)
   IF (min_eigenvalue < -tolerance) THEN
      WRITE (*, '(a,es24.16)') "FAIL: minimum rho eigenvalue = ", min_eigenvalue
      failures = failures + 1
   ELSE
      WRITE (*, '(a,es24.16)') "PASS: minimum rho eigenvalue = ", min_eigenvalue
   END IF

   trace_rho = local_trace(rho)
   trace_direct = local_trace(rho_reference)
   CALL check_real("trace/direct MT d weight", trace_rho, trace_direct, tolerance, failures)

   phase = EXP(CMPLX(0.0, 0.731))
   coefficients_phase = phase*coefficients
   CALL local_d_spin_contract_density(coefficients_phase, metric, rho_phase)
   phase_error = MAXVAL(ABS(rho_phase - rho))
   CALL check_bound("global band phase invariance", phase_error, tolerance, failures)

   radial_identical = radial
   radial_identical(:, :, :, 2) = radial_identical(:, :, :, 1)
   CALL local_d_spin_ordered_radial_metric(radial_identical, r0, dx, n_grid, metric_identical)
   norm = SQRT(metric_identical(1, 1, 1, 1))
   radial_identical(:, :, 1, :) = radial_identical(:, :, 1, :)/norm
   CALL local_d_spin_ordered_radial_metric(radial_identical, r0, dx, n_grid, metric_identical)

   DO m = -local_d_l, local_d_l
      coefficients = CMPLX(0.0, 0.0)
      coefficients(m, 1, 1) = CMPLX(1.0, 0.0)
      CALL local_d_spin_contract_density(coefficients, metric_identical, rho)
      rho_reference = CMPLX(0.0, 0.0)
      rho_reference(m, 1, m, 1) = CMPLX(1.0, 0.0)
      CALL check_bound("pure orbital m channel", MAXVAL(ABS(rho - rho_reference)), tolerance, failures)
   END DO

   coefficients = CMPLX(0.0, 0.0)
   coefficients(0, 1, 1) = CMPLX(1.0, 0.0)
   CALL local_d_spin_contract_density(coefficients, metric_identical, rho)
   CALL check_real("pure up weight", REAL(rho(0, 1, 0, 1)), 1.0, tolerance, failures)
   CALL check_bound("pure up spin leakage", MAXVAL(ABS(rho(:, 2, :, :))), tolerance, failures)

   coefficients = CMPLX(0.0, 0.0)
   coefficients(0, 2, 1) = CMPLX(1.0, 0.0)
   CALL local_d_spin_contract_density(coefficients, metric_identical, rho)
   CALL check_real("pure down weight", REAL(rho(0, 2, 0, 2)), 1.0, tolerance, failures)
   CALL check_bound("pure down spin leakage", MAXVAL(ABS(rho(:, 1, :, :))), tolerance, failures)

   coefficients = CMPLX(0.0, 0.0)
   coefficients(0, 1, 1) = CMPLX(SQRT(0.5), 0.0)
   coefficients(0, 2, 1) = CMPLX(SQRT(0.5), 0.0)
   CALL local_d_spin_contract_density(coefficients, metric_identical, rho)
   CALL check_complex("(up+down)/sqrt(2) coherence", rho(0, 1, 0, 2), CMPLX(0.5, 0.0), tolerance, failures)

   coefficients(0, 2, 1) = CMPLX(0.0, SQRT(0.5))
   CALL local_d_spin_contract_density(coefficients, metric_identical, rho)
   CALL check_complex("(up+i down)/sqrt(2) coherence", rho(0, 1, 0, 2), CMPLX(0.0, -0.5), tolerance, failures)
   sigma_y = 2.0*AIMAG(rho(0, 2, 0, 1))
   CALL check_real("sigma_y expectation", sigma_y, 1.0, tolerance, failures)

   max_abs_error = MAX(direct_error, MAX(hermiticity_error, phase_error))
   max_relative_error = max_abs_error/MAX(1.0, rho_scale)
   WRITE (*, '(a,es24.16)') "Critical S12(1,2)               = ", metric(1, 2, 1, 2)
   WRITE (*, '(a,es24.16)') "Critical S12(2,1)               = ", metric(2, 1, 1, 2)
   WRITE (*, '(a,es24.16)') "Critical ordered-overlap split  = ", &
      ABS(metric(1, 2, 1, 2) - metric(2, 1, 1, 2))
   WRITE (*, '(a,es24.16)') "Maximum selected absolute error = ", max_abs_error
   WRITE (*, '(a,es24.16)') "Maximum selected relative error = ", max_relative_error

   IF (failures /= 0) THEN
      WRITE (*, '(a,i0)') "LOCAL D-SPIN PROJECTOR TEST: FAIL; failures = ", failures
      ERROR STOP 1
   END IF
   WRITE (*, '(a)') "LOCAL D-SPIN PROJECTOR TEST: PASS"

CONTAINS

   SUBROUTINE initialize_spin_dependent_radial_functions(values)
      REAL, INTENT(OUT) :: values(:, :, :, :)
      REAL :: r
      INTEGER :: j

      values = 0.0
      DO j = 1, n_grid
         r = r0*EXP(dx*REAL(j - 1))
         values(j, 1, 1, 1) = r**3*EXP(-0.70*r)
         values(j, 1, 2, 1) = r**3*(1.0 + 0.20*r)*EXP(-1.10*r)
         values(j, 1, 3, 1) = r**4*EXP(-0.90*r)
         values(j, 2, :, 1) = [0.05*values(j, 1, 1, 1), -0.03*values(j, 1, 2, 1), &
                                0.04*values(j, 1, 3, 1)]

         values(j, 1, 1, 2) = 0.80*r**3*EXP(-0.50*r)
         values(j, 1, 2, 2) = r**4*(1.0 + 0.10*r)*EXP(-1.30*r)
         values(j, 1, 3, 2) = r**3*(0.30 + r)*EXP(-0.80*r)
         values(j, 2, :, 2) = [-0.02*values(j, 1, 1, 2), 0.06*values(j, 1, 2, 2), &
                                0.03*values(j, 1, 3, 2)]
      END DO
   END SUBROUTINE initialize_spin_dependent_radial_functions

   SUBROUTINE direct_radial_overlap(values, i_radial, j_radial, ispin_in, jspin_in, overlap)
      REAL, INTENT(IN) :: values(:, :, :, :)
      INTEGER, INTENT(IN) :: i_radial, j_radial, ispin_in, jspin_in
      REAL, INTENT(OUT) :: overlap
      REAL :: integrand(n_grid)

      integrand = values(:, 1, i_radial, ispin_in)*values(:, 1, j_radial, jspin_in) &
                + values(:, 2, i_radial, ispin_in)*values(:, 2, j_radial, jspin_in)
      CALL intgr0(integrand, r0, dx, n_grid, overlap)
   END SUBROUTINE direct_radial_overlap

   SUBROUTINE direct_density(coeff, values, density)
      COMPLEX, INTENT(IN) :: coeff(-local_d_l:, :, :)
      REAL, INTENT(IN) :: values(:, :, :, :)
      COMPLEX, INTENT(OUT) :: density(-local_d_l:local_d_l, local_d_n_spins, &
                                      -local_d_l:local_d_l, local_d_n_spins)
      COMPLEX :: wave(n_grid, 2, -local_d_l:local_d_l, local_d_n_spins)
      COMPLEX :: product(n_grid)
      REAL :: product_part(n_grid), integral_part
      INTEGER :: i_radial, component, m1, m2, s1, s2

      wave = CMPLX(0.0, 0.0)
      DO m1 = -local_d_l, local_d_l
         DO s1 = 1, local_d_n_spins
            DO i_radial = 1, n_radial
               DO component = 1, 2
                  wave(:, component, m1, s1) = wave(:, component, m1, s1) &
                     + coeff(m1, s1, i_radial)*values(:, component, i_radial, s1)
               END DO
            END DO
         END DO
      END DO

      density = CMPLX(0.0, 0.0)
      DO m1 = -local_d_l, local_d_l
         DO s1 = 1, local_d_n_spins
            DO m2 = -local_d_l, local_d_l
               DO s2 = 1, local_d_n_spins
                  product = wave(:, 1, m1, s1)*CONJG(wave(:, 1, m2, s2)) &
                          + wave(:, 2, m1, s1)*CONJG(wave(:, 2, m2, s2))
                  product_part = REAL(product)
                  CALL intgr0(product_part, r0, dx, n_grid, integral_part)
                  density(m1, s1, m2, s2) = CMPLX(integral_part, 0.0)
                  product_part = AIMAG(product)
                  CALL intgr0(product_part, r0, dx, n_grid, integral_part)
                  density(m1, s1, m2, s2) = density(m1, s1, m2, s2) + CMPLX(0.0, integral_part)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE direct_density

   REAL FUNCTION local_hermiticity_error(density)
      COMPLEX, INTENT(IN) :: density(-local_d_l:local_d_l, local_d_n_spins, &
                                     -local_d_l:local_d_l, local_d_n_spins)
      INTEGER :: m1, m2, s1, s2

      local_hermiticity_error = 0.0
      DO m1 = -local_d_l, local_d_l
         DO s1 = 1, local_d_n_spins
            DO m2 = -local_d_l, local_d_l
               DO s2 = 1, local_d_n_spins
                  local_hermiticity_error = MAX(local_hermiticity_error, &
                     ABS(density(m1, s1, m2, s2) - CONJG(density(m2, s2, m1, s1))))
               END DO
            END DO
         END DO
      END DO
   END FUNCTION local_hermiticity_error

   REAL FUNCTION local_trace(density)
      COMPLEX, INTENT(IN) :: density(-local_d_l:local_d_l, local_d_n_spins, &
                                     -local_d_l:local_d_l, local_d_n_spins)
      INTEGER :: m1, s1

      local_trace = 0.0
      DO m1 = -local_d_l, local_d_l
         DO s1 = 1, local_d_n_spins
            local_trace = local_trace + REAL(density(m1, s1, m1, s1))
         END DO
      END DO
   END FUNCTION local_trace

   SUBROUTINE hermitian_eigenvalues(density, values)
      COMPLEX, INTENT(IN) :: density(-local_d_l:local_d_l, local_d_n_spins, &
                                     -local_d_l:local_d_l, local_d_n_spins)
      REAL, INTENT(OUT) :: values(n_state)
      COMPLEX :: matrix(n_state, n_state), work(4*n_state)
      REAL :: rwork(4*n_state)
      INTEGER :: m1, m2, s1, s2, row, column, info

      DO m1 = -local_d_l, local_d_l
         DO s1 = 1, local_d_n_spins
            row = (m1 + local_d_l)*local_d_n_spins + s1
            DO m2 = -local_d_l, local_d_l
               DO s2 = 1, local_d_n_spins
                  column = (m2 + local_d_l)*local_d_n_spins + s2
                  matrix(row, column) = density(m1, s1, m2, s2)
               END DO
            END DO
         END DO
      END DO
      CALL zheev('N', 'U', n_state, matrix, n_state, values, work, SIZE(work), rwork, info)
      IF (info /= 0) THEN
         WRITE (*, '(a,i0)') "FAIL: ZHEEV returned info = ", info
         ERROR STOP 2
      END IF
   END SUBROUTINE hermitian_eigenvalues

   SUBROUTINE check_bound(name, value, bound, n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value, bound
      INTEGER, INTENT(INOUT) :: n_failures

      IF (value > bound) THEN
         WRITE (*, '(a,a,a,es24.16,a,es24.16)') "FAIL: ", TRIM(name), " error ", value, " > ", bound
         n_failures = n_failures + 1
      ELSE
         WRITE (*, '(a,a,a,es24.16)') "PASS: ", TRIM(name), " error = ", value
      END IF
   END SUBROUTINE check_bound

   SUBROUTINE check_real(name, value, expected, bound, n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value, expected, bound
      INTEGER, INTENT(INOUT) :: n_failures

      CALL check_bound(name, ABS(value - expected), bound, n_failures)
   END SUBROUTINE check_real

   SUBROUTINE check_complex(name, value, expected, bound, n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      COMPLEX, INTENT(IN) :: value, expected
      REAL, INTENT(IN) :: bound
      INTEGER, INTENT(INOUT) :: n_failures

      CALL check_bound(name, ABS(value - expected), bound, n_failures)
   END SUBROUTINE check_complex

END PROGRAM local_d_spin_projector_test
