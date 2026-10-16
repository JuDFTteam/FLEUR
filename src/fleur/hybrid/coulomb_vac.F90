!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>z-integrals of the 2D Coulomb kernel in the vacuum (films).
!>v(g;z,z') = (2 pi/g) exp(-g|z-z'|) for g > 0; at g = 0 the divergent head is left to
!>m_gamma_2d and the remainder is -2 pi |z-z'|.  The cell area enters in m_coulomb_vac_blocks.
MODULE m_coulomb_vac
   USE m_juDFT
   USE m_constants
   USE m_grule, ONLY: grule
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: vac_int_samevac, vac_int_samevac_g0, vac_exp_mom, vac_mom_g0
   PUBLIC :: vac_pref
   PUBLIC :: vac_elem_gram, vac_elem_gram_cached
   PUBLIC :: zquad
   PUBLIC :: vac_zint

CONTAINS

   PURE REAL FUNCTION vac_pref(g)
      IMPLICIT NONE
      REAL, INTENT(IN) :: g

      IF (g > 1e-12) THEN
         vac_pref = tpi_const/g
      ELSE
         vac_pref = tpi_const
      END IF
   END FUNCTION vac_pref

   !>Vacuum z-integral; the same rule as in the Coulomb elements keeps the basis orthonormal
   !>in the metric of the operator.  No tail correction: the moment elimination needs a linear rule.
   REAL FUNCTION vac_zint(prod, nz, delz)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: prod(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz

      vac_zint = zquad(prod, nz, delz)
   END FUNCTION vac_zint

   PURE REAL FUNCTION zquad(y, n, h)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: y(:)
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN)    :: h

      INTEGER :: i, nodd
      REAL    :: acc

      IF (n <= 1) THEN
         zquad = 0.0
         RETURN
      ELSE IF (n == 2) THEN
         zquad = 0.5*h*(y(1) + y(2))
         RETURN
      END IF

      nodd = n - 1 + MOD(n, 2)
      acc = y(1) + y(nodd)
      DO i = 2, nodd - 1, 2
         acc = acc + 4.0*y(i)
      END DO
      DO i = 3, nodd - 2, 2
         acc = acc + 2.0*y(i)
      END DO
      zquad = h*acc/3.0
      IF (nodd < n) zquad = zquad + 0.5*h*(y(nodd) + y(n))
   END FUNCTION zquad

   !>Galerkin matrix of the kernel on piecewise-quadratic elements over the z mesh.
   !>Symmetric and positive definite for g > 0; the kink at z = z' lies on cell boundaries.
   SUBROUTINE vac_elem_gram(nz, delz, g, a)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz, g
      REAL, INTENT(OUT)   :: a(:, :)

      INTEGER, PARAMETER :: NG = 6
      REAL    :: TT(NG), WT(NG), xh((NG + 1)/2), wh((NG + 1)/2)
      INTEGER :: ncell, nelem, ic, jc, ip, iq, iside, m, n, nni, nnj
      INTEGER :: idxi(3), idxj(3)
      REAL    :: z, zp, z0, jac, sl, kk, ww, vali(3), valj(3)

      CALL grule(NG, xh, wh)
      DO ip = 1, (NG + 1)/2
         TT(NG + 1 - ip) = 0.5*(1.0 + xh(ip));  TT(ip) = 0.5*(1.0 - xh(ip))
         WT(NG + 1 - ip) = 0.5*wh(ip);          WT(ip) = 0.5*wh(ip)
      END DO

      a(:nz, :nz) = 0.0
      ncell = nz - 1
      IF (ncell < 1) RETURN
      nelem = ncell/2

      DO ic = 1, ncell
         DO jc = 1, ncell
            IF (ic /= jc) THEN
               DO ip = 1, NG
                  z = (ic - 1 + TT(ip))*delz
                  CALL vac_elem_shapes(ic, TT(ip), nelem, nni, idxi, vali)
                  DO iq = 1, NG
                     zp = (jc - 1 + TT(iq))*delz
                     CALL vac_elem_shapes(jc, TT(iq), nelem, nnj, idxj, valj)
                     kk = vac_elem_kernel(z - zp, g)
                     ww = WT(ip)*WT(iq)*delz*delz*kk
                     DO m = 1, nni
                        DO n = 1, nnj
                           a(idxi(m), idxj(n)) = a(idxi(m), idxj(n)) + ww*vali(m)*valj(n)
                        END DO
                     END DO
                  END DO
               END DO
            ELSE
               z0 = (ic - 1)*delz
               DO ip = 1, NG
                  z = z0 + TT(ip)*delz
                  CALL vac_elem_shapes(ic, TT(ip), nelem, nni, idxi, vali)
                  DO iside = 1, 2
                     IF (iside == 1) THEN
                        jac = z - z0
                     ELSE
                        jac = z0 + delz - z
                     END IF
                     IF (jac <= 0.0) CYCLE
                     DO iq = 1, NG
                        IF (iside == 1) THEN
                           zp = z0 + jac*TT(iq)
                        ELSE
                           zp = z + jac*TT(iq)
                        END IF
                        sl = (zp - z0)/delz
                        CALL vac_elem_shapes(ic, sl, nelem, nnj, idxj, valj)
                        kk = vac_elem_kernel(z - zp, g)
                        ww = WT(ip)*delz*WT(iq)*jac*kk
                        DO m = 1, nni
                           DO n = 1, nnj
                              a(idxi(m), idxj(n)) = a(idxi(m), idxj(n)) + ww*vali(m)*valj(n)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO

      DO jc = 1, nz
         DO ic = 1, jc - 1
            a(ic, jc) = 0.5*(a(ic, jc) + a(jc, ic))
            a(jc, ic) = a(ic, jc)
         END DO
      END DO
   END SUBROUTINE vac_elem_gram

   PURE SUBROUTINE vac_elem_shapes(ic, t, nelem, nn, idx, val)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ic, nelem
      REAL, INTENT(IN)     :: t
      INTEGER, INTENT(OUT) :: nn, idx(3)
      REAL, INTENT(OUT)    :: val(3)

      INTEGER :: ielem, ibase
      REAL    :: xi

      IF (ic <= 2*nelem) THEN
         ielem = (ic + 1)/2
         ibase = 2*ielem - 1
         IF (ic == 2*ielem - 1) THEN
            xi = t
         ELSE
            xi = 1.0 + t
         END IF
         nn = 3
         idx = [ibase, ibase + 1, ibase + 2]
         val = [(xi - 1.0)*(xi - 2.0)/2.0, xi*(2.0 - xi), xi*(xi - 1.0)/2.0]
      ELSE
         nn = 2
         idx = [ic, ic + 1, 1]
         val = [1.0 - t, t, 0.0]
      END IF
   END SUBROUTINE vac_elem_shapes

   PURE REAL FUNCTION vac_elem_kernel(dz, g)
      IMPLICIT NONE
      REAL, INTENT(IN) :: dz, g

      IF (g > 1e-12) THEN
         vac_elem_kernel = EXP(-g*ABS(dz))
      ELSE
         vac_elem_kernel = -ABS(dz)
      END IF
   END FUNCTION vac_elem_kernel

   SUBROUTINE vac_elem_gram_cached(nz, delz, g, a)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz, g
      REAL, INTENT(OUT)   :: a(:, :)

      INTEGER, SAVE :: c_nz = -1
      REAL, SAVE    :: c_delz = -1.0, c_g = -2.0
      REAL, ALLOCATABLE, SAVE :: c_a(:, :)
      !$OMP THREADPRIVATE(c_nz, c_delz, c_g, c_a)

      IF (c_nz /= nz .OR. c_delz /= delz .OR. c_g /= g .OR. .NOT. ALLOCATED(c_a)) THEN
         IF (ALLOCATED(c_a)) DEALLOCATE (c_a)
         ALLOCATE (c_a(nz, nz))
         CALL vac_elem_gram(nz, delz, g, c_a)
         c_nz = nz; c_delz = delz; c_g = g
      END IF
      a(:nz, :nz) = c_a(:nz, :nz)
   END SUBROUTINE vac_elem_gram_cached

   REAL FUNCTION vac_int_samevac(f, h, nz, delz, g)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:), h(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz, g

      REAL :: a(nz, nz)

      CALL vac_elem_gram_cached(nz, delz, g, a)
      vac_int_samevac = DOT_PRODUCT(f(:nz), MATMUL(a, h(:nz)))
   END FUNCTION vac_int_samevac

   REAL FUNCTION vac_int_samevac_g0(f, h, nz, delz)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:), h(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz

      REAL :: a(nz, nz)

      CALL vac_elem_gram_cached(nz, delz, 0.0, a)
      vac_int_samevac_g0 = DOT_PRODUCT(f(:nz), MATMUL(a, h(:nz)))
   END FUNCTION vac_int_samevac_g0

   !>int f(s) exp(-g s) ds: the only functional coupling a vacuum function across the slab.
   REAL FUNCTION vac_exp_mom(f, nz, delz, g)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz, g

      INTEGER :: i
      REAL    :: pr(nz)

      DO i = 1, nz
         pr(i) = f(i)*EXP(-g*(i - 1)*delz)
      END DO
      vac_exp_mom = zquad(pr, nz, delz)
   END FUNCTION vac_exp_mom

   !>Charge and first z-moment, which carry the g = 0 coupling across the slab.
   SUBROUTINE vac_mom_g0(f, nz, delz, q0, q1)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz
      REAL, INTENT(OUT)   :: q0, q1

      INTEGER :: i
      REAL    :: pr(nz)

      pr = f(:nz)
      q0 = zquad(pr, nz, delz)
      DO i = 1, nz
         pr(i) = f(i)*(i - 1)*delz
      END DO
      q1 = zquad(pr, nz, delz)
   END SUBROUTINE vac_mom_g0

END MODULE m_coulomb_vac
