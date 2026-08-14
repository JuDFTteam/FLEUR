PROGRAM radfun_integrals_test
   USE m_intgr
   USE m_radfun_integrals
   IMPLICIT NONE

   INTEGER, PARAMETER :: jri = 13, nr = 2, lmax = 0, jspins = 2
   REAL, PARAMETER :: r0 = 0.02, dx = 0.15, tolerance = 2.0e-6
   INTEGER :: ir
   INTEGER :: n_r(0:lmax)
   REAL :: rmsh(jri), radial(jri, 2, nr, 0:lmax, jspins)
   REAL :: overlap(nr, nr, 0:lmax, jspins, jspins)
   REAL :: direct12, direct12_transposed, direct11, direct22
   REAL :: integrand(jri)

   n_r = nr
   overlap = 0.0
   DO ir = 1, jri
      rmsh(ir) = r0*EXP(dx*(ir - 1))
      radial(ir, 1, 1, 0, 1) = 0.7 + 0.8*rmsh(ir)
      radial(ir, 2, 1, 0, 1) = 0.2 + 0.3*rmsh(ir)**2
      radial(ir, 1, 2, 0, 1) = 1.1 + 1.4*rmsh(ir)**2
      radial(ir, 2, 2, 0, 1) = 0.4 + 0.9*rmsh(ir)
      radial(ir, 1, 1, 0, 2) = 1.8 + 0.5*rmsh(ir)**2
      radial(ir, 2, 1, 0, 2) = 0.6 + 1.2*rmsh(ir)
      radial(ir, 1, 2, 0, 2) = 0.3 + 2.1*rmsh(ir)
      radial(ir, 2, 2, 0, 2) = 1.5 + 0.7*rmsh(ir)**2
   END DO

   integrand = radial(:, 1, 1, 0, 1)*radial(:, 1, 2, 0, 2) &
             + radial(:, 2, 1, 0, 1)*radial(:, 2, 2, 0, 2)
   CALL intgr0(integrand, r0, dx, jri, direct12)
   integrand = radial(:, 1, 2, 0, 1)*radial(:, 1, 1, 0, 2) &
             + radial(:, 2, 2, 0, 1)*radial(:, 2, 1, 0, 2)
   CALL intgr0(integrand, r0, dx, jri, direct12_transposed)
   integrand = radial(:, 1, 1, 0, 1)*radial(:, 1, 2, 0, 1) &
             + radial(:, 2, 1, 0, 1)*radial(:, 2, 2, 0, 1)
   CALL intgr0(integrand, r0, dx, jri, direct11)
   integrand = radial(:, 1, 1, 0, 2)*radial(:, 1, 2, 0, 2) &
             + radial(:, 2, 1, 0, 2)*radial(:, 2, 2, 0, 2)
   CALL intgr0(integrand, r0, dx, jri, direct22)

   CALL calculate_radial_integrals(radial, overlap, n_r, rmsh, dx, jri, lmax, jspins)

   WRITE (*, '(A,ES24.16)') 'direct S12(1,2): ', direct12
   WRITE (*, '(A,ES24.16)') 'direct S12(2,1): ', direct12_transposed
   WRITE (*, '(A,ES24.16)') 'absolute difference: ', ABS(direct12 - direct12_transposed)
   WRITE (*, '(A,ES24.16)') 'cached S12(1,2): ', overlap(1, 2, 0, 1, 2)
   WRITE (*, '(A,ES24.16)') 'cached S12(2,1): ', overlap(2, 1, 0, 1, 2)

   CALL assert_close(overlap(1, 2, 0, 1, 2), direct12, 'ordered S12(1,2)')
   CALL assert_close(overlap(2, 1, 0, 1, 2), direct12_transposed, 'ordered S12(2,1)')
   CALL assert_close(overlap(2, 1, 0, 2, 1), direct12, 'S21(2,1) = S12(1,2)')
   CALL assert_close(overlap(1, 2, 0, 2, 1), direct12_transposed, 'S21(1,2) = S12(2,1)')
   CALL assert_close(overlap(1, 2, 0, 1, 1), direct11, 'direct S11(1,2)')
   CALL assert_close(overlap(2, 1, 0, 1, 1), direct11, 'symmetric S11')
   CALL assert_close(overlap(1, 2, 0, 2, 2), direct22, 'direct S22(1,2)')
   CALL assert_close(overlap(2, 1, 0, 2, 2), direct22, 'symmetric S22')
   IF (ABS(direct12 - direct12_transposed) < 100.0*tolerance) &
      ERROR STOP 'synthetic cross-spin overlaps are not sufficiently asymmetric'
   WRITE (*, '(A)') 'radfun_integrals_test: PASS'

CONTAINS
   SUBROUTINE assert_close(actual, expected, label)
      REAL, INTENT(IN) :: actual, expected
      CHARACTER(*), INTENT(IN) :: label
      IF (ABS(actual - expected) > tolerance*MAX(1.0, ABS(expected))) THEN
         WRITE (*, '(3A,2ES24.16)') 'FAIL: ', label, ': actual/expected = ', actual, expected
         ERROR STOP 1
      END IF
   END SUBROUTINE assert_close
END PROGRAM radfun_integrals_test
