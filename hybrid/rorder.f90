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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Orders rarr(1:n) according to size and returns a correspondingly defined pointer in pnt

   SUBROUTINE rorderp(pnt, rarr, n)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: pnt(1:n)
      REAL, INTENT(IN)  :: rarr(1:n)
      INTEGER                :: i, j, k

      DO i = 1, n
         pnt(i) = i
         DO j = 1, i - 1
            IF (rarr(pnt(j)) > rarr(i)) THEN
               DO k = i, j + 1, -1
                  pnt(k) = pnt(k - 1)
               END DO
               pnt(j) = i
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE rorderp

!     Same as rorderp but divides the problem in halves np times (leading to 2**np intervals) and is
!     much faster than rorderp (devide and conquer algorithm).
!     There is an optimal np, while for larger np the overhead (also memory-wise) outweights the speed-up.
!     np = max(0,int(log(n*0.001)/log(2.0))) should be a safe choice.
   RECURSIVE SUBROUTINE rorderpf(pnt, rarr, n, np)
      use m_judft
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n, np
      INTEGER, INTENT(OUT) :: pnt(n)
      REAL  , INTENT(IN)  :: rarr(n)
      REAL                   :: rarr1(n)
      REAL                   :: ravg
      INTEGER                :: pnt1(n), pnt2(n)
      INTEGER                :: n1, n2, i

      IF (np == 0) THEN
         CALL rorderp(pnt, rarr, n)
         RETURN
      ELSE IF (np < 0) THEN
         call judft_error("rorderpf: fourth argument must be non-negative (bug?).")
      END IF
      ravg = sum(rarr)/n
! first half
      n1 = 0
      DO i = 1, n
         IF (rarr(i) <= ravg) THEN
            n1 = n1 + 1
            rarr1(n1) = rarr(i)
            pnt1(n1) = i
         END IF
      END DO
      CALL rorderpf(pnt2, rarr1, n1, np - 1)
      pnt(:n1) = pnt1(pnt2(:n1))
! second half
      n2 = 0
      DO i = 1, n
         IF (rarr(i) > ravg) THEN
            n2 = n2 + 1
            rarr1(n2) = rarr(i)
            pnt1(n2) = i
         END IF
      END DO
      CALL rorderpf(pnt2, rarr1, n2, np - 1)
      pnt(n1 + 1:) = pnt1(pnt2(:n2))
   END SUBROUTINE rorderpf


end module m_rorder
