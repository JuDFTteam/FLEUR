!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_trapz
   !General Purpose trapezian method integration
   !Used in green's function calculations because the
   !integrands are very spiky

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: trapzr, trapzc, trapz

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