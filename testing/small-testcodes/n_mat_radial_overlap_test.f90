PROGRAM n_mat_radial_overlap_test
   USE m_constants, ONLY: lmaxU_const
   USE m_intgr
   USE m_radfun_integrals
   USE m_types
   USE m_types_abc
   USE m_types_radfun
   USE m_nmat
   IMPLICIT NONE

   INTEGER, PARAMETER :: jri = 13, nr = 2, lmax = 2, jspins = 2, l_u = 2
   REAL, PARAMETER :: r0 = 0.02, dx = 0.15, tolerance = 2.0e-6
   TYPE(t_atoms) :: atoms
   TYPE(t_sym) :: sym
   TYPE(t_radfun) :: radfun
   TYPE(t_abc) :: abc_up, abc_down, abc_cross_left, abc_cross_right
   INTEGER :: ir, j, jj
   INTEGER :: n_r(0:lmax)
   REAL :: rmsh(jri), radial(jri, 2, nr, 0:lmax, jspins)
   REAL :: direct(2, 2, 2, 2), integrand(jri), weights(1)
   REAL :: expected_cross, transposed_cross, expected_up, expected_down
   COMPLEX :: mmpMat(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, 1, 3)

   n_r = nr
   radial = 0.0
   DO ir = 1, jri
      rmsh(ir) = r0*EXP(dx*(ir - 1))
      radial(ir, 1, 1, l_u, 1) = 0.7 + 0.8*rmsh(ir)
      radial(ir, 2, 1, l_u, 1) = 0.2 + 0.3*rmsh(ir)**2
      radial(ir, 1, 2, l_u, 1) = 1.1 + 1.4*rmsh(ir)**2
      radial(ir, 2, 2, l_u, 1) = 0.4 + 0.9*rmsh(ir)
      radial(ir, 1, 1, l_u, 2) = 1.8 + 0.5*rmsh(ir)**2
      radial(ir, 2, 1, l_u, 2) = 0.6 + 1.2*rmsh(ir)
      radial(ir, 1, 2, l_u, 2) = 0.3 + 2.1*rmsh(ir)
      radial(ir, 2, 2, l_u, 2) = 1.5 + 0.7*rmsh(ir)**2
   END DO

   DO j = 1, nr
      DO jj = 1, nr
         integrand = radial(:, 1, j, l_u, 1)*radial(:, 1, jj, l_u, 1) &
                   + radial(:, 2, j, l_u, 1)*radial(:, 2, jj, l_u, 1)
         CALL intgr0(integrand, r0, dx, jri, direct(j, jj, 1, 1))
         integrand = radial(:, 1, j, l_u, 2)*radial(:, 1, jj, l_u, 2) &
                   + radial(:, 2, j, l_u, 2)*radial(:, 2, jj, l_u, 2)
         CALL intgr0(integrand, r0, dx, jri, direct(j, jj, 2, 2))
         integrand = radial(:, 1, j, l_u, 1)*radial(:, 1, jj, l_u, 2) &
                   + radial(:, 2, j, l_u, 1)*radial(:, 2, jj, l_u, 2)
         CALL intgr0(integrand, r0, dx, jri, direct(j, jj, 1, 2))
         direct(jj, j, 2, 1) = direct(j, jj, 1, 2)
      END DO
   END DO

   ALLOCATE(radfun%integral(nr, nr, 0:lmax, jspins, jspins), source=0.0)
   CALL calculate_radial_integrals(radial, radfun%integral, n_r, rmsh, dx, jri, lmax, jspins)
   atoms%n_u = 1
   ALLOCATE(atoms%neq(1), atoms%firstAtom(1), atoms%lda_u(1))
   atoms%neq = 1
   atoms%firstAtom = 0
   atoms%lda_u(1)%atomType = 1
   atoms%lda_u(1)%l = l_u
   CALL allocate_abc(abc_up)
   CALL allocate_abc(abc_down)
   CALL allocate_abc(abc_cross_left)
   CALL allocate_abc(abc_cross_right)
   abc_up%cof(1, 6, :, 1) = [CMPLX(0.6, 0.0), CMPLX(-1.1, 0.0)]
   abc_down%cof(1, 6, :, 1) = [CMPLX(0.9, 0.0), CMPLX(0.4, 0.0)]
   abc_cross_left%cof(1, 7, 1, 1) = CMPLX(1.25, 0.0)
   abc_cross_right%cof(1, 6, 2, 1) = CMPLX(-0.75, 0.0)
   weights = 0.8
   mmpMat = CMPLX(0.0, 0.0)
   CALL n_mat(atoms, radfun, sym, 1, weights, abc_up, abc_up, mmpMat(:, :, :, 1), 1, 1, 1)
   CALL n_mat(atoms, radfun, sym, 1, weights, abc_down, abc_down, mmpMat(:, :, :, 2), 1, 2, 2)
   CALL n_mat(atoms, radfun, sym, 1, weights, abc_cross_left, abc_cross_right, mmpMat(:, :, :, 3), 1, 2, 1)

   expected_up = direct_contraction(abc_up, abc_up, direct(:, :, 1, 1), weights(1))
   expected_down = direct_contraction(abc_down, abc_down, direct(:, :, 2, 2), weights(1))
   expected_cross = weights(1)*1.25*(-0.75)*direct(1, 2, 2, 1)
   transposed_cross = weights(1)*1.25*(-0.75)*direct(2, 1, 2, 1)
   WRITE (*, '(A,ES24.16)') 'direct ordered S21(1,2): ', direct(1, 2, 2, 1)
   WRITE (*, '(A,ES24.16)') 'cached ordered S21(1,2): ', radfun%integral(1, 2, l_u, 2, 1)
   WRITE (*, '(A,ES24.16)') 'incorrect transposed S21(2,1): ', direct(2, 1, 2, 1)
   WRITE (*, '(A,ES24.16)') 'expected mmpMat(0,1,1,3): ', expected_cross
   WRITE (*, '(A,2ES24.16)') 'actual mmpMat(0,1,1,3): ', mmpMat(0, 1, 1, 3)
   WRITE (*, '(A,ES24.16)') 'incorrect transposed-overlap result: ', transposed_cross
   WRITE (*, '(A,2ES24.16)') 'same-spin up actual/expected: ', REAL(mmpMat(0, 0, 1, 1)), expected_up
   WRITE (*, '(A,2ES24.16)') 'same-spin down actual/expected: ', REAL(mmpMat(0, 0, 1, 2)), expected_down
   CALL assert_close(REAL(mmpMat(0, 1, 1, 3)), expected_cross, 'ordered cross-spin n_mat contraction')
   CALL assert_close(AIMAG(mmpMat(0, 1, 1, 3)), 0.0, 'cross-spin n_mat imaginary part')
   IF (ABS(REAL(mmpMat(0, 1, 1, 3)) - transposed_cross) < 100.0*tolerance) &
      ERROR STOP 'cross-spin n_mat result does not distinguish the transposed overlap'
   CALL assert_close(REAL(mmpMat(0, 0, 1, 1)), expected_up, 'same-spin-up n_mat contraction')
   CALL assert_close(REAL(mmpMat(0, 0, 1, 2)), expected_down, 'same-spin-down n_mat contraction')
   WRITE (*, '(A)') 'n_mat_radial_overlap_test: PASS'

CONTAINS
   SUBROUTINE allocate_abc(abc)
      TYPE(t_abc), INTENT(OUT) :: abc
      ALLOCATE(abc%cof(1, 0:l_u*(l_u + 2), nr, 1), source=CMPLX(0.0, 0.0))
   END SUBROUTINE allocate_abc
   REAL FUNCTION direct_contraction(left, right, ordered_overlap, weight)
      TYPE(t_abc), INTENT(IN) :: left, right
      REAL, INTENT(IN) :: ordered_overlap(nr, nr), weight
      INTEGER :: i, k
      direct_contraction = 0.0
      DO i = 1, nr
         DO k = 1, nr
            direct_contraction = direct_contraction &
                               + weight*REAL(CONJG(left%cof(1, 6, i, 1))*right%cof(1, 6, k, 1))*ordered_overlap(i, k)
         END DO
      END DO
   END FUNCTION direct_contraction
   SUBROUTINE assert_close(actual, expected, label)
      REAL, INTENT(IN) :: actual, expected
      CHARACTER(*), INTENT(IN) :: label
      IF (ABS(actual - expected) > tolerance*MAX(1.0, ABS(expected))) THEN
         WRITE (*, '(3A,2ES24.16)') 'FAIL: ', label, ': actual/expected = ', actual, expected
         ERROR STOP 1
      END IF
   END SUBROUTINE assert_close
END PROGRAM n_mat_radial_overlap_test
