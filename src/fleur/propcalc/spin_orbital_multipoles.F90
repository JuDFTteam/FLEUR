!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Nordstrom/Bultmark complex spin-orbital multipole transformation.
!
! The density is the explicit physical ket-bra density
!
!   rho(m,s,m',s') = C(m,s) CONJG(C(m',s')),
!
! in the complex spherical-harmonic basis m=-l,...,+l and spin order
! s=1 (+1/2), s=2 (-1/2).  In particular, this is not the historical FLEUR
! mmpMat convention, which contains an orbital partial transpose and must not be
! passed directly to this module.
!
! The dense coupled array uses w(k,p,r,t), with invalid combinations set to
! exactly zero.  Valid indices are k=0,...,2*l, p=0,1,
! r=ABS(k-p),...,k+p and t=-r,...,+r.  The public convention is the complex
! convention of Nordstrom (2020), Eqs. (45), (48)-(52), and (59)-(62).
!--------------------------------------------------------------------------------
MODULE m_spin_orbital_multipoles
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: spin_orbital_density_to_multipoles
   PUBLIC :: spin_orbital_multipoles_to_density

   ! The first public API version is deliberately limited to the shell ranks
   ! covered by repository validation: s, p, d, and f (l=0,...,3).  Higher l
   ! values are mathematically possible but are not advertised for the private
   ! direct-factorial Wigner evaluator used here.
   INTEGER, PARAMETER :: supported_l_max = 3

CONTAINS

   SUBROUTINE spin_orbital_density_to_multipoles(l, rho, w)
      INTEGER, INTENT(IN) :: l
      COMPLEX, INTENT(IN) :: rho(-l:, 1:, -l:, 1:)
      COMPLEX, INTENT(OUT) :: w(0:, 0:, 0:, -(2*l+1):)

      COMPLEX :: wxy(0:2*supported_l_max, 0:1, -2*supported_l_max:2*supported_l_max, -1:1)
      COMPLEX :: chi(2, 2, 0:1, -1:1)
      COMPLEX :: bracket, value
      REAL :: orbital_factor
      INTEGER :: k, p, r, t, x, y, ma, mb, sa, sb

      CALL validate_l(l)
      IF (ANY(SHAPE(rho) /= [2*l+1, 2, 2*l+1, 2])) THEN
         CALL juDFT_error("rho must have shape (2*l+1,2,2*l+1,2)", &
                          calledby="spin_orbital_density_to_multipoles")
      END IF
      IF (ANY(SHAPE(w) /= [2*l+1, 2, 2*l+2, 4*l+3])) THEN
         CALL juDFT_error("w must have shape (2*l+1,2,2*l+2,4*l+3)", &
                          calledby="spin_orbital_density_to_multipoles")
      END IF

      CALL spherical_pauli_matrices(chi)
      wxy = CMPLX(0.0, 0.0)
      w = CMPLX(0.0, 0.0)

      ! Nordstrom Eqs. (48)-(50): uncoupled orbital-spin double tensor.
      DO k = 0, 2*l
         DO p = 0, 1
            DO x = -k, k
               DO y = -p, p
                  value = CMPLX(0.0, 0.0)
                  DO ma = -l, l
                     DO mb = -l, l
                        orbital_factor = phase_factor(l-mb)*integer_wigner_3j(l, k, l, -mb, x, ma) &
                                         / orbital_normalization(l, k)
                        DO sa = 1, 2
                           DO sb = 1, 2
                              value = value + orbital_factor*chi(sb, sa, p, y)*rho(ma, sa, mb, sb)
                           END DO
                        END DO
                     END DO
                  END DO
                  wxy(k, p, x, y) = value
               END DO
            END DO
         END DO
      END DO

      ! Nordstrom Eq. (59): couple the orbital and spin tensor ranks.
      DO k = 0, 2*l
         DO p = 0, 1
            DO r = ABS(k-p), k+p
               bracket = generalized_three_bracket(k, p, r)
               DO t = -r, r
                  value = CMPLX(0.0, 0.0)
                  DO x = -k, k
                     DO y = -p, p
                        value = value + wxy(k, p, x, y)*integer_wigner_3j(k, r, p, -x, t, -y) &
                                        *phase_factor(-x-y)
                     END DO
                  END DO
                  w(k, p, r, t) = phase_factor(k+p)*value/bracket
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE spin_orbital_density_to_multipoles

   SUBROUTINE spin_orbital_multipoles_to_density(l, w, rho)
      INTEGER, INTENT(IN) :: l
      COMPLEX, INTENT(IN) :: w(0:, 0:, 0:, -(2*l+1):)
      COMPLEX, INTENT(OUT) :: rho(-l:, 1:, -l:, 1:)

      COMPLEX :: wxy(0:2*supported_l_max, 0:1, -2*supported_l_max:2*supported_l_max, -1:1)
      COMPLEX :: chi(2, 2, 0:1, -1:1)
      COMPLEX :: bracket, value, spin_dual
      REAL :: orbital_dual
      INTEGER :: k, p, r, t, x, y, ma, mc, sa, sc

      CALL validate_l(l)
      IF (ANY(SHAPE(w) /= [2*l+1, 2, 2*l+2, 4*l+3])) THEN
         CALL juDFT_error("w must have shape (2*l+1,2,2*l+2,4*l+3)", &
                          calledby="spin_orbital_multipoles_to_density")
      END IF
      IF (ANY(SHAPE(rho) /= [2*l+1, 2, 2*l+1, 2])) THEN
         CALL juDFT_error("rho must have shape (2*l+1,2,2*l+1,2)", &
                          calledby="spin_orbital_multipoles_to_density")
      END IF

      wxy = CMPLX(0.0, 0.0)
      ! Nordstrom Eq. (62): analytic inverse of the rank coupling.
      DO k = 0, 2*l
         DO p = 0, 1
            DO x = -k, k
               DO y = -p, p
                  value = CMPLX(0.0, 0.0)
                  DO r = ABS(k-p), k+p
                     bracket = generalized_three_bracket(k, p, r)
                     DO t = -r, r
                        value = value + REAL(2*r+1)*integer_wigner_3j(k, r, p, -x, t, -y) &
                                        *phase_factor(x-k+y-p)*bracket*w(k, p, r, t)
                     END DO
                  END DO
                  wxy(k, p, x, y) = value
               END DO
            END DO
         END DO
      END DO

      CALL spherical_pauli_matrices(chi)
      rho = CMPLX(0.0, 0.0)
      ! Nordstrom Eq. (52).  For explicit spin 1/2 the spin part is written as
      ! the exact spherical-Pauli dual, (-1)^y chi_-y/2, avoiding any
      ! half-integer Wigner dependency.
      DO ma = -l, l
         DO mc = -l, l
            DO sa = 1, 2
               DO sc = 1, 2
                  value = CMPLX(0.0, 0.0)
                  DO k = 0, 2*l
                     DO x = -k, k
                        orbital_dual = REAL(2*k+1)*orbital_normalization(l, k)*phase_factor(mc-l) &
                                       *integer_wigner_3j(l, k, l, -mc, x, ma)
                        DO p = 0, 1
                           DO y = -p, p
                              spin_dual = 0.5*phase_factor(y)*chi(sa, sc, p, -y)
                              value = value + orbital_dual*spin_dual*wxy(k, p, x, y)
                           END DO
                        END DO
                     END DO
                  END DO
                  rho(ma, sa, mc, sc) = value
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE spin_orbital_multipoles_to_density

   SUBROUTINE validate_l(l)
      INTEGER, INTENT(IN) :: l

      IF (l < 0 .OR. l > supported_l_max) THEN
         CALL juDFT_error("this API version supports only l=0,...,3", &
                          calledby="spin_orbital_multipoles")
      END IF
   END SUBROUTINE validate_l

   PURE SUBROUTINE spherical_pauli_matrices(chi)
      COMPLEX, INTENT(OUT) :: chi(2, 2, 0:1, -1:1)

      chi = CMPLX(0.0, 0.0)
      chi(1, 1, 0, 0) = 1.0
      chi(2, 2, 0, 0) = 1.0
      chi(1, 1, 1, 0) = 1.0
      chi(2, 2, 1, 0) = -1.0
      chi(1, 2, 1, 1) = -SQRT(2.0)
      chi(2, 1, 1, -1) = SQRT(2.0)
   END SUBROUTINE spherical_pauli_matrices

   PURE REAL FUNCTION orbital_normalization(l, k)
      INTEGER, INTENT(IN) :: l, k

      orbital_normalization = factorial_real(2*l) &
                              /SQRT(factorial_real(2*l-k)*factorial_real(2*l+k+1))
   END FUNCTION orbital_normalization

   PURE COMPLEX FUNCTION generalized_three_bracket(k, p, r)
      INTEGER, INTENT(IN) :: k, p, r

      INTEGER :: g
      REAL :: magnitude

      g = k + p + r
      IF (MOD(g, 2) == 0) THEN
         generalized_three_bracket = CMPLX(integer_wigner_3j(k, p, r, 0, 0, 0), 0.0)
      ELSE
         ! Nordstrom Eq. (61), including the i^g phase for odd g.
         magnitude = SQRT(factorial_real(g-2*k)*factorial_real(g-2*p)*factorial_real(g-2*r) &
                          /factorial_real(g+1))
         magnitude = magnitude*double_factorial_real(g) &
                     /(double_factorial_real(g-2*k)*double_factorial_real(g-2*p) &
                       *double_factorial_real(g-2*r))
         generalized_three_bracket = CMPLX(0.0, 1.0)**g*magnitude
      END IF
   END FUNCTION generalized_three_bracket

   PURE REAL FUNCTION integer_wigner_3j(j1, j2, j3, m1, m2, m3)
      INTEGER, INTENT(IN) :: j1, j2, j3, m1, m2, m3

      INTEGER :: z, zmin, zmax
      REAL :: delta_factor, magnetic_factor, racah_sum

      ! Integer-only Racah factorial expression.  Keeping this evaluator
      ! private avoids the caller-owned factorial tables and unrelated module
      ! dependencies of FLEUR's legacy m_util implementation.
      integer_wigner_3j = 0.0
      IF (m1+m2+m3 /= 0) RETURN
      IF (ABS(m1) > j1 .OR. ABS(m2) > j2 .OR. ABS(m3) > j3) RETURN
      IF (j3 < ABS(j1-j2) .OR. j3 > j1+j2) RETURN

      delta_factor = SQRT(factorial_real(j1+j2-j3)*factorial_real(j1-j2+j3) &
                          *factorial_real(-j1+j2+j3)/factorial_real(j1+j2+j3+1))
      magnetic_factor = SQRT(factorial_real(j1+m1)*factorial_real(j1-m1) &
                             *factorial_real(j2+m2)*factorial_real(j2-m2) &
                             *factorial_real(j3+m3)*factorial_real(j3-m3))
      zmin = MAX(0, j2-j3-m1, j1-j3+m2)
      zmax = MIN(j1+j2-j3, j1-m1, j2+m2)
      racah_sum = 0.0
      DO z = zmin, zmax
         racah_sum = racah_sum + phase_factor(z) &
            /(factorial_real(z)*factorial_real(j1+j2-j3-z) &
              *factorial_real(j1-m1-z)*factorial_real(j2+m2-z) &
              *factorial_real(j3-j2+m1+z)*factorial_real(j3-j1-m2+z))
      END DO
      integer_wigner_3j = phase_factor(j1-j2-m3)*delta_factor*magnetic_factor*racah_sum
   END FUNCTION integer_wigner_3j

   PURE REAL FUNCTION factorial_real(n)
      INTEGER, INTENT(IN) :: n

      INTEGER :: i

      factorial_real = 1.0
      DO i = 2, n
         factorial_real = factorial_real*REAL(i)
      END DO
   END FUNCTION factorial_real

   PURE REAL FUNCTION double_factorial_real(n)
      INTEGER, INTENT(IN) :: n

      INTEGER :: i

      double_factorial_real = 1.0
      DO i = n, 1, -2
         double_factorial_real = double_factorial_real*REAL(i)
      END DO
   END FUNCTION double_factorial_real

   PURE REAL FUNCTION phase_factor(n)
      INTEGER, INTENT(IN) :: n

      phase_factor = MERGE(1.0, -1.0, MOD(ABS(n), 2) == 0)
   END FUNCTION phase_factor

END MODULE m_spin_orbital_multipoles
