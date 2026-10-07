!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_coreerr

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: coreerr
   CONTAINS
   
   SUBROUTINE coreerr(err,var,s,nsol,pow,qow,piw,qiw)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE THE MISMATCH OF THE RADIAL WAVE FUNCTIONS AT THE     *
!   *   POINT  NMATCH  FOR OUT- AND INWARD INTEGRATION                 *
!   *                                                                  *
!   ********************************************************************

      IMPLICIT NONE

      INTEGER nsol,s

      REAL err(4),piw(2,2),pow(2,2),qiw(2,2),qow(2,2),var(4)

      INTEGER t

      err(1) = pow(s,s) - piw(s,s)*var(2)
      err(2) = qow(s,s) - qiw(s,s)*var(2)

      IF (nsol.EQ.1) RETURN

      t = 3 - s

      err(1) = err(1) + pow(s,t)*var(3) - piw(s,t)*var(2)*var(4)
      err(2) = err(2) + qow(s,t)*var(3) - qiw(s,t)*var(2)*var(4)
      err(3) = pow(t,s) - piw(t,s)*var(2) + pow(t,t)*var(3) - piw(t,t)*var(2)*var(4)
      err(4) = qow(t,s) - qiw(t,s)*var(2) + qow(t,t)*var(3) - qiw(t,t)*var(2)*var(4)

   END SUBROUTINE coreerr

END MODULE m_coreerr
