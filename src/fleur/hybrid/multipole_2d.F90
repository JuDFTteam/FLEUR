!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>In-plane Fourier transform of a point-multipole potential at height dz above it:
!>W_lm = (2 pi/k) exp(-k|dz|) exp(-i m phi_k) sum_g d(l,m,g) (ik)^(l-g) (-k sgn dz)^g.
MODULE m_multipole_2d
   USE m_juDFT
   USE m_constants
   USE m_ylm
   USE m_solid_harmonics, ONLY: solidharm_cyl_coeff

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: multipole_2d_pot

CONTAINS

   SUBROUTINE multipole_2d_pot(lmax, kv, dz, w)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: lmax
      REAL, INTENT(IN)     :: kv(2), dz
      COMPLEX, INTENT(OUT) :: w((lmax + 1)**2)

      INTEGER :: l, m, g, lm
      REAL    :: kk, phik, sgz, base
      REAL    :: d(0:lmax, -lmax:lmax, 0:lmax)
      COMPLEX :: csum

      w = cmplx_0
      kk = SQRT(kv(1)**2 + kv(2)**2)
      IF (kk < 1e-12) RETURN

      phik = ATAN2(kv(2), kv(1))
      sgz = SIGN(1.0, dz)
      base = (tpi_const/kk)*EXP(-kk*ABS(dz))

      CALL solidharm_cyl_coeff(lmax, d)

      DO l = 0, lmax
         lm = l**2
         DO m = -l, l
            lm = lm + 1
            csum = cmplx_0
            DO g = 0, l
               IF (d(l, m, g) == 0.0) CYCLE
               csum = csum + d(l, m, g)*(ImagUnit*kk)**(l - g)*(-kk*sgz)**g
            END DO
            w(lm) = base*EXP(-ImagUnit*m*phik)*csum
         END DO
      END DO
   END SUBROUTINE multipole_2d_pot

END MODULE m_multipole_2d
