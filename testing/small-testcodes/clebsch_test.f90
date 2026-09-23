!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

! Regression test for the Clebsch-Gordan selection rule, see issue #804.
!
! <j1 m1; j2 m2 | J M> vanishes unless m1+m2 = M.  A truncating test on the real
! difference used to let a mismatch of -1 through and returned a non-zero
! coefficient; one of the cases below returned sqrt(2), which is impossible for
! a Clebsch-Gordan coefficient.

MODULE m_clebsch_test
   USE m_clebsch
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_clebsch_checks

   INTEGER :: n_fail = 0

CONTAINS

   !> Returns the number of failed cases.
   INTEGER FUNCTION run_clebsch_checks()

      n_fail = 0

      ! Forbidden couplings: m1+m2 /= M, must be exactly zero.  The first two
      ! are the regression; before the fix they returned 5.7735026919E-01 and
      ! 1.4142135624E+00.
      CALL expect(1.0, 0.5, -1.0, 0.5, 1.5,  0.5, 0.0, "d=-1 (was 1/sqrt3)")
      CALL expect(1.0, 0.5,  0.0, 0.5, 1.5,  1.5, 0.0, "d=-1 (was sqrt2)")
      CALL expect(1.0, 0.5,  1.0, 0.5, 1.5, -0.5, 0.0, "d=+2")
      CALL expect(1.0, 0.5, -1.0, 0.5, 1.5,  1.5, 0.0, "d=-2")
      CALL expect(1.0, 0.5,  0.0, 0.5, 1.5,  1.0, 0.0, "d=-0.5")

      ! Allowed couplings must keep their values: the guard must not
      ! over-reject.  Signs follow the Condon-Shortley convention.  The J=3/2
      ! and J=1/2 entries for m1=-1, m2=+1/2 satisfy the sum rule
      ! 1/3 + 2/3 = 1, which fixes the magnitudes independently of the sign
      ! convention.
      CALL expect(1.0, 0.5,  1.0, 0.5, 1.5,  1.5,  1.0,            "stretched state")
      CALL expect(1.0, 0.5, -1.0, 0.5, 1.5, -0.5,  SQRT(1.0/3.0),  "d=0, J=3/2")
      CALL expect(1.0, 0.5,  0.0, 0.5, 1.5,  0.5,  SQRT(2.0/3.0),  "d=0, sqrt(2/3)")
      CALL expect(1.0, 0.5, -1.0, 0.5, 0.5, -0.5, -SQRT(2.0/3.0),  "d=0, J=1/2")

      run_clebsch_checks = n_fail
   END FUNCTION run_clebsch_checks

   SUBROUTINE expect(aj, bj, am, bm, cj, cm, want, label)
      REAL, INTENT(IN) :: aj, bj, am, bm, cj, cm, want
      CHARACTER(LEN=*), INTENT(IN) :: label

      REAL, PARAMETER :: tol = 1.0e-10
      REAL :: got

      got = clebsch(aj, bj, am, bm, cj, cm)
      IF (ABS(got-want) > tol) THEN
         n_fail = n_fail + 1
         WRITE (*,'(a,a,a,6f6.1,a,es18.10,a,es18.10)') "FAIL ", label, &
            " args=", aj, bj, am, bm, cj, cm, " got=", got, " want=", want
      ELSE
         WRITE (*,'(a,a)') "ok   ", label
      END IF
   END SUBROUTINE expect

END MODULE m_clebsch_test

PROGRAM clebsch_test
   USE m_clebsch_test, ONLY: run_clebsch_checks
   IMPLICIT NONE

   INTEGER :: n_fail

   n_fail = run_clebsch_checks()

   IF (n_fail == 0) THEN
      WRITE (*,'(a)') "CLEBSCH SELECTION RULE TEST: PASS"
   ELSE
      WRITE (*,'(a,i0,a)') "CLEBSCH SELECTION RULE TEST: FAIL (", n_fail, " case(s))"
      ERROR STOP 1
   END IF

END PROGRAM clebsch_test
