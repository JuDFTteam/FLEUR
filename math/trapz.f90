MODULE m_trapz
   !General Purpose trapezian method integration
   !Used in green's function calculations because the
   !integrands are very spiky

   IMPLICIT NONE

   CONTAINS

   PURE REAL FUNCTION trapz(y,h,n)

      REAL,          INTENT(IN)     :: y(:)

      INTEGER,       INTENT(IN)     :: n
      REAL,          INTENT(IN)     :: h


      INTEGER i

      trapz = y(1)
      DO i = 2, n-1
         trapz = trapz + 2*y(i)
      ENDDO
      trapz = trapz + y(n)

      trapz = trapz*h/2.0

   END FUNCTION trapz

END MODULE m_trapz