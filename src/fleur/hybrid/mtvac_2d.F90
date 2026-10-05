!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>MT-VAC element for in-plane |q+G||| > 0: the (l,m) multipole of an MT carrier seen by a
!>vacuum function.  mtvac_shape excludes moment, mu and the decay exp(-k(z1 - sigma z_a)).
MODULE m_mtvac_2d
   USE m_juDFT
   USE m_constants
   USE m_ylm
   USE m_multipole_2d, ONLY: multipole_2d_pot

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mtvac_shape

CONTAINS

   COMPLEX FUNCTION mtvac_shape(l, m, kv, sigma, area)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, m
      REAL, INTENT(IN)    :: kv(2), sigma, area

      INTEGER :: lm
      REAL    :: kk, dfac, ii
      COMPLEX, ALLOCATABLE :: w(:)

      kk = SQRT(kv(1)**2 + kv(2)**2)
      ALLOCATE (w((l + 1)**2))
      CALL multipole_2d_pot(l, kv, sigma, w)
      dfac = 1.0
      DO lm = 1, l
         dfac = dfac*(2*lm - 1)
      END DO
      lm = l**2 + l - m + 1
      ii = (-1)**m/((-1)**l*dfac)
      mtvac_shape = (fpi_const/(2*l + 1))*ii*w(lm)*EXP(kk*ABS(sigma))/SQRT(area)
   END FUNCTION mtvac_shape

END MODULE m_mtvac_2d
