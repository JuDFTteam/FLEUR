!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Solid harmonics r^l Y_lm as polynomials in x, y, z (fitted once), evaluated at complex
!>null vectors (analytic continuation of conjg(Y_lm)); also their cylindrical coefficients.
MODULE m_solid_harmonics
   USE m_juDFT
   USE m_constants
   USE m_ylm
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: solidharm_init, solidharm_eval_conj, solidharm_cyl_coeff

   INTEGER, ALLOCATABLE, SAVE :: mono_a(:, :), mono_b(:, :), mono_c(:, :)
   COMPLEX, ALLOCATABLE, SAVE :: mono_coeff(:, :)
   INTEGER, SAVE              :: solidharm_lmax = -1

CONTAINS

   PURE INTEGER FUNCTION nmono(l)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l
      nmono = (l + 1)*(l + 2)/2
   END FUNCTION nmono

   SUBROUTINE solidharm_init(lmax)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: lmax

      INTEGER :: l, m, lm, i, j, a, b, npt, nmx, nm, info, rank, lwork
      REAL    :: gold, phi, ct, st, kv(3)
      REAL,    ALLOCATABLE :: amat(:, :), rhs(:, :), acopy(:, :), sv(:), work(:)
      REAL,    ALLOCATABLE :: ksamp(:, :)
      COMPLEX, ALLOCATABLE :: ysamp(:, :), yl(:)

      IF (solidharm_lmax >= lmax) RETURN

      nmx = nmono(lmax)
      npt = MAX(4*nmx, 40)

      IF (ALLOCATED(mono_a)) DEALLOCATE (mono_a, mono_b, mono_c, mono_coeff)

      ALLOCATE (mono_a(nmx, 0:lmax), mono_b(nmx, 0:lmax), mono_c(nmx, 0:lmax), source=0)
      ALLOCATE (mono_coeff(nmx, (lmax + 1)**2), source=cmplx_0)
      ALLOCATE (ksamp(3, npt), ysamp((lmax + 1)**2, npt), yl((lmax + 1)**2))

      gold = pi_const*(3.0 - SQRT(5.0))
      DO i = 1, npt
         ct = 1.0 - 2.0*(i - 0.5)/npt
         st = SQRT(MAX(0.0, 1.0 - ct*ct))
         phi = gold*i
         kv = [st*COS(phi), st*SIN(phi), ct]
         ksamp(:, i) = kv
         CALL ylm4(lmax, kv, yl)
         ysamp(:, i) = yl
      END DO

      DO l = 0, lmax
         j = 0
         DO a = l, 0, -1
            DO b = l - a, 0, -1
               j = j + 1
               mono_a(j, l) = a; mono_b(j, l) = b; mono_c(j, l) = l - a - b
            END DO
         END DO
      END DO

      DO l = 0, lmax
         nm = nmono(l)
         ALLOCATE (amat(npt, nm), acopy(npt, nm), rhs(npt, 2), sv(nm), work(1))

         DO i = 1, npt
            DO j = 1, nm
               acopy(i, j) = ksamp(1, i)**mono_a(j, l)*ksamp(2, i)**mono_b(j, l) &
                             *ksamp(3, i)**mono_c(j, l)
            END DO
         END DO

         amat = acopy; rhs = 0.0
         CALL DGELSS(npt, nm, 2, amat, npt, rhs, npt, sv, -1.0, rank, work, -1, info)
         lwork = MAX(1, NINT(work(1)))
         DEALLOCATE (work); ALLOCATE (work(lwork))

         DO m = -l, l
            lm = l**2 + l + m + 1
            amat = acopy
            rhs(:, 1) = REAL(ysamp(lm, :)); rhs(:, 2) = AIMAG(ysamp(lm, :))
            CALL DGELSS(npt, nm, 2, amat, npt, rhs, npt, sv, -1.0, rank, work, lwork, info)
            IF (info /= 0) CALL juDFT_error("DGELSS failed", calledby="solidharm_init")
            IF (rank < nm) CALL juDFT_error("rank-deficient monomial fit", &
                                            calledby="solidharm_init", &
                                            hint="sample directions degenerate for this degree")
            mono_coeff(1:nm, lm) = CMPLX(rhs(1:nm, 1), rhs(1:nm, 2))
         END DO
         DEALLOCATE (amat, acopy, rhs, sv, work)
      END DO

      solidharm_lmax = lmax
   END SUBROUTINE solidharm_init

   SUBROUTINE solidharm_cyl_coeff(lmax, d)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: lmax
      REAL, INTENT(OUT)   :: d(0:lmax, -lmax:lmax, 0:lmax)

      INTEGER :: l, m, g, j, n, info
      INTEGER, ALLOCATABLE :: ipiv(:)
      REAL, ALLOCATABLE    :: vdm(:, :), rhs(:, :), znod(:)
      REAL    :: v(3)
      COMPLEX :: y((lmax + 1)**2)

      d = 0.0

      DO l = 0, lmax
         n = l + 1
         ALLOCATE (vdm(n, n), rhs(n, 2*l + 1), znod(n), ipiv(n))

         DO j = 1, n
            znod(j) = COS(pi_const*(j - 0.5)/n)
            DO g = 0, l
               vdm(j, g + 1) = znod(j)**g
            END DO
            v = [1.0, 0.0, znod(j)]
            CALL ylm4(l, v, y)
            DO m = -l, l
               rhs(j, m + l + 1) = NORM2(v)**l*REAL(y(l**2 + l + m + 1))
            END DO
         END DO

         CALL dgesv(n, 2*l + 1, vdm, n, ipiv, rhs, n, info)
         IF (info /= 0) CALL juDFT_error("solid harmonic Vandermonde solve failed", &
                                         calledby="solidharm_cyl_coeff")
         DO m = -l, l
            DO g = 0, l
               d(l, m, g) = rhs(g + 1, m + l + 1)
            END DO
         END DO

         DEALLOCATE (vdm, rhs, znod, ipiv)
      END DO
   END SUBROUTINE solidharm_cyl_coeff

   COMPLEX FUNCTION solidharm_eval_conj(l, m, k) RESULT(s)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, m
      COMPLEX, INTENT(IN) :: k(3)

      INTEGER :: lm, j

      IF (l > solidharm_lmax) CALL juDFT_error("table not built to this l", &
                                               calledby="solidharm_eval_conj")
      lm = l**2 + l + m + 1
      s = cmplx_0
      DO j = 1, nmono(l)
         s = s + CONJG(mono_coeff(j, lm))*k(1)**mono_a(j, l)*k(2)**mono_b(j, l) &
             *k(3)**mono_c(j, l)
      END DO
   END FUNCTION solidharm_eval_conj

END MODULE m_solid_harmonics
