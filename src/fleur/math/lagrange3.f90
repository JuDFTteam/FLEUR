!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_lagrange3

   !! Quadratic (3-point) Lagrange interpolation on an equidistant mesh.
   !!
   !! Given the values fLeft, fCenter, fRight of a function at three consecutive
   !! mesh points z, z+delz, z+2*delz, lagrange3 returns its interpolation at
   !! z + relPos*delz. relPos is thus the distance from the left mesh point in
   !! units of the mesh spacing; the three points it interpolates between are at
   !! relPos = 0, 1 and 2.

   IMPLICIT NONE
   PRIVATE

   INTERFACE lagrange3
      PROCEDURE :: lagrange3Real, lagrange3Complex
   END INTERFACE

   PUBLIC :: lagrange3

CONTAINS

   PURE REAL FUNCTION lagrange3Real(relPos, fLeft, fCenter, fRight)
      REAL, INTENT(IN) :: relPos
      REAL, INTENT(IN) :: fLeft, fCenter, fRight

      lagrange3Real = 0.5*(relPos-1.)*(relPos-2.)*fLeft &
                    - relPos*(relPos-2.)*fCenter &
                    + 0.5*relPos*(relPos-1.)*fRight

   END FUNCTION lagrange3Real

   PURE COMPLEX FUNCTION lagrange3Complex(relPos, fLeft, fCenter, fRight)
      REAL,    INTENT(IN) :: relPos
      COMPLEX, INTENT(IN) :: fLeft, fCenter, fRight

      lagrange3Complex = 0.5*(relPos-1.)*(relPos-2.)*fLeft &
                       - relPos*(relPos-2.)*fCenter &
                       + 0.5*relPos*(relPos-1.)*fRight

   END FUNCTION lagrange3Complex

END MODULE m_lagrange3
