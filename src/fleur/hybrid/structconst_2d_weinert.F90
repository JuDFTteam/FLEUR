!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Low-l 2D structure constants with a compact Weinert pseudocharge: one reciprocal sum and
!>no splitting parameter.  Independent of m_structconst_2d_parry; used as a cross-check.
MODULE m_structconst_2d_weinert
   USE m_juDFT
   USE m_constants
   USE m_ylm
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: weinert_structconst_2d

   INTEGER, PARAMETER :: NPC_DEFAULT = 8

   INTEGER, PARAMETER :: NTHETA = 16, NPHI = 32

   REAL, PARAMETER :: F_RPC = 0.4, F_R0 = 0.4

   INTEGER, PARAMETER :: NZQ = 24

   INTEGER, PARAMETER :: MAXSHELL = 200
   REAL, PARAMETER    :: GTOL = 1e-13

   COMPLEX, ALLOCATABLE, SAVE :: ygrid(:, :)
   INTEGER, SAVE              :: ygrid_lmax = -1

CONTAINS

   PURE REAL FUNCTION bes_j0(x)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL :: ax, z, xx, y, p, q
      ax = ABS(x)
      IF (ax < 8.0) THEN
         y = x*x
         p = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7 + &
             y*(-11214424.18 + y*(77392.33017 + y*(-184.9052456)))))
         q = 57568490411.0 + y*(1029532985.0 + y*(9494680.718 + &
             y*(59272.64853 + y*(267.8532712 + y))))
         bes_j0 = p/q
      ELSE
         z = 8.0/ax; y = z*z; xx = ax - 0.785398164
         p = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4 + &
             y*(-0.2073370639e-5 + y*0.2093887211e-6)))
         q = -0.1562499995e-1 + y*(0.1430488765e-3 + y*(-0.6911147651e-5 + &
             y*(0.7621095161e-6 + y*(-0.934935152e-7))))
         bes_j0 = SQRT(0.636619772/ax)*(COS(xx)*p - z*SIN(xx)*q)
      END IF
   END FUNCTION bes_j0

   PURE REAL FUNCTION bes_j1(x)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL :: ax, z, xx, y, p, q, r
      ax = ABS(x)
      IF (ax < 8.0) THEN
         y = x*x
         p = x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1 + &
             y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))))
         q = 144725228442.0 + y*(2300535178.0 + y*(18583304.74 + &
             y*(99447.43394 + y*(376.9991397 + y))))
         bes_j1 = p/q
      ELSE
         z = 8.0/ax; y = z*z; xx = ax - 2.356194491
         p = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4 + &
             y*(0.2457520174e-5 + y*(-0.240337019e-6))))
         q = 0.04687499995 + y*(-0.2002690873e-3 + y*(0.8449199096e-5 + &
             y*(-0.88228987e-6 + y*0.105787412e-6)))
         r = SQRT(0.636619772/ax)*(COS(xx)*p - z*SIN(xx)*q)
         bes_j1 = MERGE(-r, r, x < 0.0)
      END IF
   END FUNCTION bes_j1

   PURE REAL FUNCTION bes_jn(n, x)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN)    :: x

      INTEGER, PARAMETER :: IACC = 40
      REAL, PARAMETER    :: BIGNO = 1.0e10, BIGNI = 1.0e-10
      INTEGER :: j, mm
      REAL    :: ax, bj, bjm, bjp, sm, tox, ans

      ax = ABS(x)
      IF (n == 0) THEN
         bes_jn = bes_j0(x); RETURN
      ELSE IF (n == 1) THEN
         bes_jn = bes_j1(x); RETURN
      END IF
      IF (ax < 1e-30) THEN
         bes_jn = 0.0; RETURN
      END IF

      tox = 2.0/ax
      IF (ax > REAL(n)) THEN
         bjm = bes_j0(ax); bj = bes_j1(ax)
         DO j = 1, n - 1
            bjp = j*tox*bj - bjm
            bjm = bj; bj = bjp
         END DO
         ans = bj
      ELSE
         mm = 2*((n + INT(SQRT(REAL(IACC*n))))/2)
         ans = 0.0; sm = 0.0
         bjp = 0.0; bj = 1.0
         DO j = mm, 1, -1
            bjm = j*tox*bj - bjp
            bjp = bj; bj = bjm
            IF (ABS(bj) > BIGNO) THEN
               bj = bj*BIGNI; bjp = bjp*BIGNI
               ans = ans*BIGNI; sm = sm*BIGNI
            END IF
            IF (MOD(j, 2) == 1) sm = sm + bj
            IF (j == n) ans = bjp
         END DO
         sm = 2.0*sm - bj
         ans = ans/sm
      END IF
      bes_jn = MERGE(-ans, ans, x < 0.0 .AND. MOD(n, 2) == 1)
   END FUNCTION bes_jn

   SUBROUTINE gauleg_unit(n, x, w)
      USE m_grule, ONLY: grule
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(OUT)   :: x(n), w(n)

      INTEGER :: i, m
      REAL    :: xh((n + 1)/2), wh((n + 1)/2)

      m = (n + 1)/2
      CALL grule(n, xh, wh)
      DO i = 1, m
         x(n + 1 - i) = xh(i);  x(i) = -xh(i)
         w(n + 1 - i) = wh(i);  w(i) = wh(i)
      END DO
   END SUBROUTINE gauleg_unit

   PURE REAL FUNCTION pc_norm(rpc, npc)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: rpc
      INTEGER, INTENT(IN) :: npc
      INTEGER :: j
      REAL    :: ii
      ii = 1.0/3.0
      DO j = 1, npc
         ii = ii*(2.0*j)/(2.0*j + 3.0)
      END DO
      pc_norm = 1.0/(fpi_const*rpc**3*ii)
   END FUNCTION pc_norm

   PURE REAL FUNCTION pc_sigma(k, zeta, rpc, npc, cnrm)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: k, zeta, rpc, cnrm
      INTEGER, INTENT(IN) :: npc

      INTEGER :: j
      REAL    :: b2, b, fac, x

      b2 = rpc*rpc - zeta*zeta
      IF (b2 <= 0.0) THEN
         pc_sigma = 0.0; RETURN
      END IF
      b = SQRT(b2)
      fac = 1.0
      DO j = 1, npc
         fac = fac*2.0*j
      END DO
      fac = tpi_const*cnrm*rpc**(-2*npc)*fac

      x = k*b
      IF (x < 1e-8) THEN
         pc_sigma = fac*b**(2*npc + 2)
         DO j = 1, npc + 1
            pc_sigma = pc_sigma/(2.0*j)
         END DO
      ELSE
         pc_sigma = fac*b**(npc + 1)*bes_jn(npc + 1, x)/k**(npc + 1)
      END IF
   END FUNCTION pc_sigma

   PURE REAL FUNCTION pc_pot_k(k, z, zr, rpc, npc, cnrm, xq, wq)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: k, z, zr, rpc, cnrm, xq(:), wq(:)
      INTEGER, INTENT(IN) :: npc

      INTEGER :: i, ns
      REAL    :: cut, acc, zp, ker, a1, b1, hw, mid, seg(3)

      cut = MIN(MAX(z, zr - rpc), zr + rpc)
      seg = [zr - rpc, cut, zr + rpc]

      acc = 0.0
      DO ns = 1, 2
         a1 = seg(ns); b1 = seg(ns + 1)
         IF (b1 - a1 < 1e-14) CYCLE
         hw = 0.5*(b1 - a1); mid = 0.5*(b1 + a1)
         DO i = 1, SIZE(xq)
            zp = mid + hw*xq(i)
            IF (k > 0.0) THEN
               ker = EXP(-k*ABS(z - zp))/k
            ELSE
               ker = -ABS(z - zp)
            END IF
            acc = acc + hw*wq(i)*ker*pc_sigma(k, zp - zr, rpc, npc, cnrm)
         END DO
      END DO
      pc_pot_k = tpi_const*acc
   END FUNCTION pc_pot_k

   PURE REAL FUNCTION pc_pot_r(r, rpc, npc, cnrm, xq, wq)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: r, rpc, cnrm, xq(:), wq(:)
      INTEGER, INTENT(IN) :: npc

      INTEGER :: i
      REAL    :: acc1, acc2, hw, mid, s

      acc1 = 0.0
      hw = 0.5*r; mid = 0.5*r
      DO i = 1, SIZE(xq)
         s = mid + hw*xq(i)
         acc1 = acc1 + hw*wq(i)*cnrm*(1.0 - (s/rpc)**2)**npc*s*s
      END DO
      acc2 = 0.0
      hw = 0.5*(rpc - r); mid = 0.5*(rpc + r)
      DO i = 1, SIZE(xq)
         s = mid + hw*xq(i)
         acc2 = acc2 + hw*wq(i)*cnrm*(1.0 - (s/rpc)**2)**npc*s
      END DO
      pc_pot_r = fpi_const*(acc1/MAX(r, 1e-30) + acc2)
   END FUNCTION pc_pot_r

   SUBROUTINE build_ygrid(lmax, ct, cos_phi, sin_phi)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: lmax
      REAL, INTENT(IN)    :: ct(:), cos_phi(:), sin_phi(:)

      INTEGER :: it, ip, k
      REAL    :: st, rvec(3)
      COMPLEX :: y((lmax + 1)**2)

      IF (ygrid_lmax >= lmax .AND. ALLOCATED(ygrid)) RETURN
      IF (ALLOCATED(ygrid)) DEALLOCATE (ygrid)
      ALLOCATE (ygrid(NTHETA*NPHI, (lmax + 1)**2))

      DO it = 1, NTHETA
         st = SQRT(MAX(0.0, 1.0 - ct(it)**2))
         DO ip = 1, NPHI
            rvec = [st*cos_phi(ip), st*sin_phi(ip), ct(it)]
            CALL ylm4(lmax, rvec, y)
            k = (ip - 1)*NTHETA + it
            ygrid(k, :) = y
         END DO
      END DO
      ygrid_lmax = lmax
   END SUBROUTINE build_ygrid

   SUBROUTINE weinert_structconst_2d(lmax, qrel, r0, amat, bmat, npc_in, s)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: lmax
      REAL, INTENT(IN)     :: qrel(3), r0(3), amat(3, 3), bmat(3, 3)
      INTEGER, INTENT(IN)  :: npc_in
      COMPLEX, INTENT(OUT) :: s((lmax + 1)**2)

      INTEGER :: npc, n1, n2, ish, it, ip, l, m, lm, nquiet
      REAL    :: qc(3), gv(3), kv(3), kk, area, dmin, rpc, rr0, cnrm
      REAL    :: ct(NTHETA), wt(NTHETA), xz(NZQ), wz(NZQ)
      REAL    :: st, phi, shell, running, rvec(3), vth(NTHETA)
      REAL    :: cos_phi(NPHI), sin_phi(NPHI)
      COMPLEX :: vgrid(NTHETA, NPHI), ph, y((lmax + 1)**2)
      COMPLEX :: vw(NTHETA*NPHI)
      LOGICAL :: l_onsite

      npc = MERGE(NPC_DEFAULT, npc_in, npc_in <= 0)
      s = cmplx_0

      area = ABS(amat(1, 1)*amat(2, 2) - amat(1, 2)*amat(2, 1))
      IF (area < 1e-12) CALL juDFT_error("degenerate in-plane cell", &
                                         calledby="weinert_structconst_2d")
      qc = MATMUL(qrel, bmat)
      qc(3) = 0.0
      l_onsite = NORM2(r0) < 1e-10

      dmin = HUGE(1.0)
      DO n1 = -3, 3
         DO n2 = -3, 3
            rvec = n1*amat(:, 1) + n2*amat(:, 2) + r0
            IF (NORM2(rvec) > 1e-10) dmin = MIN(dmin, NORM2(rvec))
         END DO
      END DO
      IF (dmin > 1e29) CALL juDFT_error("no lattice image found", &
                                        calledby="weinert_structconst_2d")
      rpc = F_RPC*dmin
      rr0 = F_R0*dmin

      CALL gauleg_unit(NTHETA, ct, wt)
      CALL gauleg_unit(NZQ, xz, wz)
      cnrm = pc_norm(rpc, npc)
      DO ip = 1, NPHI
         phi = tpi_const*(ip - 1)/NPHI
         cos_phi(ip) = COS(phi); sin_phi(ip) = SIN(phi)
      END DO

      vgrid = cmplx_0
      running = 0.0
      nquiet = 0
      DO ish = 0, MAXSHELL
         shell = 0.0
         DO n1 = -ish, ish
            DO n2 = -ish, ish
               IF (MAX(ABS(n1), ABS(n2)) /= ish) CYCLE
               gv = n1*bmat(1, :) + n2*bmat(2, :)
               kv = qc + gv
               kv(3) = 0.0
               kk = NORM2(kv)
               IF (kk < 1e-10) kk = 0.0

               DO it = 1, NTHETA
                  vth(it) = pc_pot_k(kk, rr0*ct(it), r0(3), rpc, npc, cnrm, xz, wz)/area
               END DO

               DO it = 1, NTHETA
                  st = SQRT(MAX(0.0, 1.0 - ct(it)**2))
                  DO ip = 1, NPHI
                     ph = EXP(ImagUnit*(kv(1)*(rr0*st*cos_phi(ip) - r0(1)) &
                                        + kv(2)*(rr0*st*sin_phi(ip) - r0(2))))
                     vgrid(it, ip) = vgrid(it, ip) + ph*vth(it)
                  END DO
               END DO
               shell = MAX(shell, MAXVAL(ABS(vth)))
            END DO
         END DO
         running = MAX(running, shell)
         IF (ish >= 2 .AND. shell < GTOL*MAX(running, 1e-30)) THEN
            nquiet = nquiet + 1
            IF (nquiet >= 2) EXIT
         ELSE
            nquiet = 0
         END IF
      END DO

      IF (l_onsite) THEN
         vgrid = vgrid - pc_pot_r(rr0, rpc, npc, cnrm, xz, wz)
      END IF

      CALL build_ygrid(lmax, ct, cos_phi, sin_phi)
      DO it = 1, NTHETA
         DO ip = 1, NPHI
            vw((ip - 1)*NTHETA + it) = wt(it)*(tpi_const/NPHI)*vgrid(it, ip)
         END DO
      END DO
      CALL zgemv('C', NTHETA*NPHI, (lmax + 1)**2, cmplx_1, ygrid, NTHETA*NPHI, &
                 vw, 1, cmplx_0, s, 1)

      DO l = 0, lmax
         lm = l**2
         DO m = -l, l
            lm = lm + 1
            s(lm) = s(lm)*(2*l + 1)/(fpi_const*rr0**l)
         END DO
      END DO
   END SUBROUTINE weinert_structconst_2d

END MODULE m_structconst_2d_weinert
