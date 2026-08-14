!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_radfun_integrals
   USE m_intgr
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: calculate_radial_integrals

CONTAINS

   SUBROUTINE calculate_radial_integrals(r, integral, n_r, rmsh, dx, jri, lmax, jspins)
      REAL,    INTENT(IN)    :: r(:, :, :, 0:, :)
      REAL,    INTENT(INOUT) :: integral(:, :, 0:, :, :)
      INTEGER, INTENT(IN)    :: n_r(0:)
      REAL,    INTENT(IN)    :: rmsh(:), dx
      INTEGER, INTENT(IN)    :: jri, lmax, jspins

      INTEGER :: ispin, jspin, l, i, j
      REAL, ALLOCATABLE :: rf(:)
      REAL :: ovlp

      DO ispin = 1, jspins
         DO jspin = 1, jspins
            DO l = 0, lmax
               IF (ispin == jspin) THEN
                  DO i = 1, n_r(l)
                     DO j = 1, i
                        rf = r(1:jri, 1, i, l, ispin)*r(1:jri, 1, j, l, jspin) &
                           + r(1:jri, 2, i, l, ispin)*r(1:jri, 2, j, l, jspin)
                        CALL intgr0(rf, rmsh(1), dx, jri, ovlp)
                        integral(i, j, l, ispin, jspin) = ovlp
                        integral(j, i, l, ispin, jspin) = ovlp
                     END DO
                  END DO
               ELSE IF (ispin < jspin) THEN
                  DO i = 1, n_r(l)
                     DO j = 1, n_r(l)
                        rf = r(1:jri, 1, i, l, ispin)*r(1:jri, 1, j, l, jspin) &
                           + r(1:jri, 2, i, l, ispin)*r(1:jri, 2, j, l, jspin)
                        CALL intgr0(rf, rmsh(1), dx, jri, ovlp)
                        integral(i, j, l, ispin, jspin) = ovlp
                        integral(j, i, l, jspin, ispin) = ovlp
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE calculate_radial_integrals

END MODULE m_radfun_integrals
