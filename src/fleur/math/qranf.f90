!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_qranf
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qranf
CONTAINS
   REAL FUNCTION qranf(x,j)
!     **********************************************************
!     quasi random generator in the interval (0.,1.)
!     **********************************************************

      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL x
      INTEGER j
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC aint
!     ..
      j = j + 1
      qranf = j*x
      qranf = qranf - aint(qranf)
      RETURN
   END FUNCTION
END MODULE m_qranf
