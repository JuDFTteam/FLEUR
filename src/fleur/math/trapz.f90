MODULE m_trapz
   !General Purpose trapezian method integration
   !Used in green's function calculations because the
   !integrands are very spiky

   IMPLICIT NONE

   INTERFACE trapz
      PROCEDURE :: trapzr, trapzc
   END INTERFACE

   CONTAINS

   PURE REAL FUNCTION trapzr(y,h,n)

      REAL,          INTENT(IN)     :: y(:)

      INTEGER,       INTENT(IN)     :: n
      REAL,          INTENT(IN)     :: h

      INTEGER i

      trapzr = y(1)
      DO i = 2, n-1
         trapzr = trapzr + 2*y(i)
      ENDDO
      trapzr = trapzr + y(n)

      trapzr = trapzr*h/2.0

   END FUNCTION trapzr

   PURE COMPLEX FUNCTION trapzc(y,h,n)

      COMPLEX,       INTENT(IN)     :: y(:)

      INTEGER,       INTENT(IN)     :: n
      REAL,          INTENT(IN)     :: h

      INTEGER i

      trapzc = y(1)
      DO i = 2, n-1
         trapzc = trapzc + 2*y(i)
      ENDDO
      trapzc = trapzc + y(n)

      trapzc = trapzc*h/2.0

   END FUNCTION trapzc


END MODULE m_trapz