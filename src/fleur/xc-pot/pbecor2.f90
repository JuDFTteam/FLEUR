!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_pbecor2
!---------------------
! slimmed down version of gcor used in pw91 routines, to interpolate
! lsd correlation energy, as given by (10) of
! j. p. perdew and y. wang, phys. rev. b {\bf 45}, 13244 (1992).
! k. burke, may 11, 1996.
!---------------------
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: pbecor2
CONTAINS
   SUBROUTINE pbecor2( &
      a,a1,b1,b2,b3,b4,rtrs, &
      gg,ggrs)
      IMPLICIT NONE

      REAL, INTENT (IN)  :: a,a1,b1,b2,b3,b4,rtrs
      REAL, INTENT (OUT) :: gg,ggrs

      REAL :: q0,q1,q2,q3
!     ..
      q0 = -2.e0*a* (1.e0+a1*rtrs*rtrs)
      q1 = 2.e0*a*rtrs* (b1+rtrs* (b2+rtrs* (b3+b4*rtrs)))
      q2 = log(1.e0+1.e0/q1)
      gg = q0*q2
      q3 = a* (b1/rtrs+2.e0*b2+rtrs* (3.e0*b3+4.e0*b4*rtrs))
      ggrs = -2.e0*a*a1*q2 - q0*q3/ (q1* (1.e0+q1))

   END SUBROUTINE pbecor2
END MODULE m_pbecor2
