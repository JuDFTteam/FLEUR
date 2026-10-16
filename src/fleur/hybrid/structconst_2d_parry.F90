!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Low-l 2D structure constants by Ewald splitting (Parry): real-space sum with incomplete
!>gamma functions plus a reciprocal sum of z-derivatives of the 2D kernel.  The result must
!>not depend on the splitting parameter eta.
MODULE m_structconst_2d_parry
   USE m_juDFT
   USE m_constants
   USE m_ylm
   USE m_solid_harmonics, ONLY: solidharm_cyl_coeff
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: parry_structconst_2d, parry_default_eta

CONTAINS

   PURE REAL FUNCTION ew_term(a, b)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a, b

      REAL :: x

      x = a + b
      IF (x >= 0.0) THEN
         ew_term = EXP(-a*a - b*b)*ERFC_SCALED(x)
      ELSE
         ew_term = 2.0*EXP(2.0*a*b) - EXP(-a*a - b*b)*ERFC_SCALED(-x)
      END IF
   END FUNCTION ew_term

   PURE SUBROUTINE gamma_ratio(lmax, x, gr)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: lmax
      REAL, INTENT(IN)    :: x
      REAL, INTENT(OUT)   :: gr(0:lmax)

      INTEGER :: l
      REAL    :: gam, gfull, a

      gam = SQRT(pi_const)*ERFC(SQRT(x))
      gfull = SQRT(pi_const)
      gr(0) = gam/gfull

      DO l = 0, lmax - 1
         a = l + 0.5
         gam = a*gam + x**a*EXP(-x)
         gfull = a*gfull
         gr(l + 1) = gam/gfull
      END DO
   END SUBROUTINE gamma_ratio

   PURE SUBROUTINE hermite(n, u, h)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN)    :: u
      REAL, INTENT(OUT)   :: h(0:n)

      INTEGER :: i

      h(0) = 1.0
      IF (n >= 1) h(1) = 2.0*u
      DO i = 1, n - 1
         h(i + 1) = 2.0*u*h(i) - 2.0*i*h(i - 1)
      END DO
   END SUBROUTINE hermite

   PURE REAL FUNCTION parry_default_eta(area)
      IMPLICIT NONE
      REAL, INTENT(IN) :: area

      parry_default_eta = SQRT(pi_const)/SQRT(area)
   END FUNCTION parry_default_eta

   SUBROUTINE parry_structconst_2d(lmax, qrel, r0, amat, bmat, eta_in, s)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: lmax
      REAL, INTENT(IN)     :: qrel(3), r0(3), amat(3, 3), bmat(3, 3), eta_in
      COMPLEX, INTENT(OUT) :: s((lmax + 1)**2)

      INTEGER :: n1, n2, l, m, g, lm, nm1, nm2
      REAL    :: eta, area, rcut, kcut, r, z, a, b, kk, phik, cross(3)
      LOGICAL :: l_self
      REAL    :: a1(3), a2(3), b1(3), b2(3), rv(3), kv(3), h1, h2
      REAL    :: dfac(0:lmax), gr(0:lmax), herm(0:lmax), wn
      REAL    :: hn(0:lmax), pn(0:lmax), hd(0:lmax)
      REAL    :: d(0:lmax, -lmax:lmax, 0:lmax)
      COMPLEX :: y((lmax + 1)**2), cexp, csum, opref

      s = cmplx_0
      z = r0(3)

      a1 = amat(:, 1); a2 = amat(:, 2)
      b1 = bmat(1, :); b2 = bmat(2, :)
      cross(1) = a1(2)*a2(3) - a1(3)*a2(2)
      cross(2) = a1(3)*a2(1) - a1(1)*a2(3)
      cross(3) = a1(1)*a2(2) - a1(2)*a2(1)
      area = NORM2(cross)
      IF (area < 1e-12) CALL juDFT_error("degenerate in-plane cell", &
                                         calledby="parry_structconst_2d")

      eta = eta_in
      IF (eta <= 0.0) eta = parry_default_eta(area)
      l_self = .FALSE.

      dfac(0) = 1.0
      DO l = 1, lmax
         dfac(l) = dfac(l - 1)*(2*l - 1)
      END DO

      CALL solidharm_cyl_coeff(lmax, d)

      rcut = 6.0/eta + NORM2(r0(1:2))
      kcut = 12.0*eta

      h1 = area/MAX(NORM2(a2), 1e-30)
      h2 = area/MAX(NORM2(a1), 1e-30)
      nm1 = CEILING(rcut/h1) + 1
      nm2 = CEILING(rcut/h2) + 1

      DO n1 = -nm1, nm1
         DO n2 = -nm2, nm2
            rv = n1*a1 + n2*a2 + r0
            r = NORM2(rv)
            IF (r < 1e-12) THEN
               l_self = .TRUE.
               CYCLE
            END IF
            IF (r > rcut) CYCLE

            CALL gamma_ratio(lmax, (eta*r)**2, gr)
            cexp = EXP(ImagUnit*tpi_const*(qrel(1)*n1 + qrel(2)*n2))
            CALL ylm4(lmax, rv, y)
            y = CONJG(y)

            DO l = 0, lmax
               lm = l**2
               DO m = -l, l
                  lm = lm + 1
                  s(lm) = s(lm) + cexp*y(lm)*r**(-l - 1)*gr(l)
               END DO
            END DO
         END DO
      END DO

      ! remove the self term of the Gaussian at T+r0 = 0
      IF (l_self) s(1) = s(1) - (2.0*eta/SQRT(pi_const))/sfp_const

      h1 = tpi_const/MAX(NORM2(a1), 1e-30)
      h2 = tpi_const/MAX(NORM2(a2), 1e-30)
      nm1 = CEILING(kcut/h1) + 1
      nm2 = CEILING(kcut/h2) + 1

      DO n1 = -nm1, nm1
         DO n2 = -nm2, nm2
            kv = (n1 - qrel(1))*b1 + (n2 - qrel(2))*b2
            kk = NORM2(kv)

            ! G|| = 0 term without the divergent 2 pi/(A k): finite remainder and its z-derivatives
            IF (kk < 1e-10) THEN
               CALL hermite(MAX(lmax, 1), eta*z, herm)
               DO l = 0, lmax
                  IF (l == 0) THEN
                     wn = -(tpi_const/area)*(z*ERF(eta*z) + EXP(-(eta*z)**2)/(eta*SQRT(pi_const)))
                  ELSE IF (l == 1) THEN
                     wn = -(tpi_const/area)*ERF(eta*z)
                  ELSE
                     wn = -(tpi_const/area)*(2.0*eta/SQRT(pi_const)) &
                          *(-eta)**(l - 2)*herm(l - 2)*EXP(-(eta*z)**2)
                  END IF
                  lm = l**2 + l + 1
                  s(lm) = s(lm) + ((-1)**l/dfac(l))*d(l, 0, l)*wn
               END DO
               CYCLE
            END IF

            IF (kk > kcut) CYCLE

            a = kk/(2.0*eta)
            b = eta*z
            phik = ATAN2(kv(2), kv(1))
            cexp = EXP(ImagUnit*(kv(1)*r0(1) + kv(2)*r0(2)))

            hn(0) = ew_term(a, b)
            pn(0) = ew_term(a, -b)
            CALL hermite(lmax, eta*z, herm)
            DO g = 0, lmax - 1
               wn = EXP(-a*a)*(-eta)**g*herm(g)*EXP(-(eta*z)**2)
               hn(g + 1) = kk*hn(g) - (2.0*eta/SQRT(pi_const))*wn
               pn(g + 1) = -kk*pn(g) + (2.0*eta/SQRT(pi_const))*wn
            END DO
            hd = hn + pn

            DO l = 0, lmax
               opref = ((-1)**l/dfac(l))*(pi_const/area)*cexp
               lm = l**2
               DO m = -l, l
                  lm = lm + 1
                  csum = cmplx_0
                  DO g = 0, l
                     IF (d(l, m, g) == 0.0) CYCLE
                     csum = csum + d(l, m, g)*(ImagUnit*kk)**(l - g)*hd(g)/kk
                  END DO
                  s(lm) = s(lm) + opref*EXP(-ImagUnit*m*phik)*csum
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE parry_structconst_2d

END MODULE m_structconst_2d_parry
