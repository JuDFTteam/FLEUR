module m_rorder
   implicit none

CONTAINS
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !     Orders iarr(1:n) according to size and returns a correspondingly defined pointer in pnt
   SUBROUTINE iorderp(pnt, iarr, n)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: pnt(1:n)
      INTEGER, INTENT(IN)  :: iarr(1:n)
      INTEGER                :: i, j, k

      DO i = 1, n
         pnt(i) = i
         DO j = 1, i - 1
            IF (iarr(pnt(j)) > iarr(i)) THEN
               DO k = i, j + 1, -1
                  pnt(k) = pnt(k - 1)
               END DO
               pnt(j) = i
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE iorderp

end module m_rorder
