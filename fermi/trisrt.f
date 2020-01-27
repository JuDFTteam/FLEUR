      MODULE m_trisrt
!
!     orders (k1,x1), (k2,x2), (k3,x3) such that x1 < x2 < x3
!
      CONTAINS
      SUBROUTINE trisrt(x1,x2,x3,k1,k2,k3)

      IMPLICIT NONE
      INTEGER, INTENT (INOUT) :: k1,k2,k3
      REAL,    INTENT (INOUT) :: x1,x2,x3
      INTEGER k
      REAL    x

      IF (x2 < x1) THEN             ! interchange x1 and x2
        x = x1 ; x1 = x2 ; x2 = x
        k = k1 ; k1 = k2 ; k2 = k
      ENDIF
      IF (x3 < x1) THEN             ! interchange x1 and x3
        x = x1 ; x1 = x3 ; x3 = x
        k = k1 ; k1 = k3 ; k3 = k
      ENDIF
      IF (x3 < x2) THEN             ! interchange x2 and x3
        x = x2 ; x2 = x3 ; x3 = x
        k = k2 ; k2 = k3 ; k3 = k
      ENDIF

      END SUBROUTINE trisrt
      END MODULE m_trisrt
