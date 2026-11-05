!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfprep_potconv
   !Conversion of a bulk potential/density to the setup of a layer.
   !
   !Interstitial: the bulk function, taken in the window |z-zc|<dt/2 and
   !smoothly made periodic in the buffer d/2<|z-zc|<dt/2, is expanded with
   !the layer period dt and projected onto the layer stars.
   !For a film the function is the plane-wave expansion for |z|<z1 and the
   !vacuum expansion beyond; the window transform is done by quadrature.
   !Muffin-tins: the lattice-harmonic expansion of the corresponding bulk
   !atom is projected onto the lattice harmonics of the layer (lower site
   !symmetry, possibly different local frame) by angular quadrature.
   USE m_juDFT
   USE m_constants
   USE m_types_sphhar
   USE m_types_stars
   USE m_types_cell
   USE m_types_sym
   USE m_types_atoms
   USE m_types_vacuum
   USE m_ylm
   USE m_gfprep_bulk, ONLY: t_gfprep_struct
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gfprep_pwconv, gfprep_pwconv_film, gfprep_vacemb, gfprep_mtconv, gfprep_check

   !composite Gauss-Legendre quadrature of the film transform: points per
   !panel and largest panel length (a.u.)
   INTEGER, PARAMETER :: nglf = 10
   REAL, PARAMETER    :: hpanel = 0.2

CONTAINS

   SUBROUTINE gfprep_pwconv(bs, bcell, ls, lcell, zc, d, pw_in, pw_out)
      TYPE(t_stars), INTENT(IN) :: bs, ls       !bulk and layer stars
      TYPE(t_cell), INTENT(IN)  :: bcell, lcell
      REAL, INTENT(IN)          :: zc           !layer center in the bulk frame
      REAL, INTENT(IN)          :: d            !region kept unchanged, |z|<d/2
      COMPLEX, INTENT(IN)       :: pw_in(:)
      COMPLEX, INTENT(OUT)      :: pw_out(:)

      COMPLEX :: cb(-bs%mx3:bs%mx3), tz(-bs%mx3:bs%mx3, -ls%mx3:ls%mx3), c
      INTEGER :: k1, k2, k3, k3l, is

      CALL priv_window(bs%mx3, bcell%bmat(3, 3), ls%mx3, lcell%amat(3, 3), d, tz)
      pw_out = CMPLX(0.0, 0.0)
      DO k2 = -MIN(ls%mx2, bs%mx2), MIN(ls%mx2, bs%mx2)
         DO k1 = -MIN(ls%mx1, bs%mx1), MIN(ls%mx1, bs%mx1)
            IF (ALL(ls%ig(k1, k2, :) == 0)) CYCLE
            !bulk coefficients of the individual G, shifted to the layer center
            cb = CMPLX(0.0, 0.0)
            DO k3 = -bs%mx3, bs%mx3
               is = bs%ig(k1, k2, k3)
               IF (is == 0) CYCLE
               cb(k3) = pw_in(is)*bs%rgphs(k1, k2, k3)*EXP(CMPLX(0.0, k3*bcell%bmat(3, 3)*zc))
            END DO
            DO k3l = -ls%mx3, ls%mx3
               is = ls%ig(k1, k2, k3l)
               IF (is == 0) CYCLE
               c = SUM(cb*tz(:, k3l))
               pw_out(is) = pw_out(is) + c*CONJG(ls%rgphs(k1, k2, k3l))
            END DO
         END DO
      END DO
      pw_out(:ls%ng3) = pw_out(:ls%ng3)/ls%nstr(:ls%ng3)
   END SUBROUTINE gfprep_pwconv

   SUBROUTINE priv_window(mx3, b3, mx3l, dt, d, tz)
      !tz(k,k') = 1/dt * int_{-dt/2}^{dt/2} g_k(z) exp(-i k' b3' z) dz with the
      !layer reciprocal vector b3'=2pi/dt and g_k the bulk plane wave
      !exp(i k b3 z) for |z|<d/2. In the buffer d/2<|z|<dt/2 g_k is blended
      !smoothly into its periodic image, so that the layer function is
      !continuous and smooth at z=+-dt/2 and does not suffer from Gibbs
      !oscillations of a hard cut.
      INTEGER, INTENT(IN)  :: mx3, mx3l
      REAL, INTENT(IN)     :: b3, dt, d
      COMPLEX, INTENT(OUT) :: tz(-mx3:mx3, -mx3l:mx3l)

      INTEGER, PARAMETER :: ngl = 64
      REAL, ALLOCATABLE  :: xg(:), wg(:)
      INTEGER :: k, kl, ig, side
      REAL    :: q, b3l, z, u, w, h, z0

      b3l = tpi_const/dt
      CALL priv_gauleg(ngl, xg, wg)
      h = 0.5*(dt - d)
      DO kl = -mx3l, mx3l
         DO k = -mx3, mx3
            q = k*b3 - kl*b3l
            !inner region
            IF (ABS(q) > 1.0E-12) THEN
               tz(k, kl) = 2.0*SIN(0.5*q*d)/q
            ELSE
               tz(k, kl) = d
            END IF
            !buffers
            IF (h > 1.0E-8) THEN
               DO side = -1, 1, 2
                  z0 = side*0.5*(d + h)
                  DO ig = 1, ngl
                     z = z0 + 0.5*h*xg(ig)
                     u = (ABS(z) - 0.5*d)/h
                     w = 1.0 - 0.5*SIN(0.5*pi_const*u)**2
                     tz(k, kl) = tz(k, kl) + 0.5*h*wg(ig)*EXP(CMPLX(0.0, -kl*b3l*z))* &
                                 (w*EXP(CMPLX(0.0, k*b3*z)) + (1.0 - w)*EXP(CMPLX(0.0, k*b3*(z - side*dt))))
                  END DO
               END DO
            END IF
            tz(k, kl) = tz(k, kl)/dt
         END DO
      END DO
   END SUBROUTINE priv_window

   SUBROUTINE gfprep_pwconv_film(fs, fcell, fvac, ls, lcell, zc, d, pw_in, vac_in, pw_out)
      !Film version of gfprep_pwconv: the film function is sampled in the
      !layer window |z-zc|<dt/2, blended in the buffer as in priv_window, and
      !projected onto the layer plane waves by quadrature
      TYPE(t_stars), INTENT(IN)  :: fs, ls       !film and layer stars
      TYPE(t_cell), INTENT(IN)   :: fcell, lcell
      TYPE(t_vacuum), INTENT(IN) :: fvac
      REAL, INTENT(IN)           :: zc, d
      COMPLEX, INTENT(IN)        :: pw_in(:), vac_in(:, :, :)   !vac(nmzd,ng2,2)
      COMPLEX, INTENT(OUT)       :: pw_out(:)

      REAL, ALLOCATABLE    :: zq(:), wq(:), za(:), zb(:), wa(:)
      COMPLEX, ALLOCATABLE :: eq(:, :), g(:)
      COMPLEX :: c
      INTEGER :: k1, k2, k3l, is, nq, q
      REAL    :: dt, b3l

      dt = lcell%amat(3, 3)
      b3l = tpi_const/dt
      CALL priv_film_nodes(fcell%z1, zc, d, dt, zq, wq, za, zb, wa)
      nq = SIZE(zq)
      ALLOCATE (eq(nq, -ls%mx3:ls%mx3), g(nq))
      DO k3l = -ls%mx3, ls%mx3
         eq(:, k3l) = wq*EXP(CMPLX(0.0, -k3l*b3l*zq))/dt
      END DO
      pw_out = CMPLX(0.0, 0.0)
      DO k2 = -ls%mx2, ls%mx2
         DO k1 = -ls%mx1, ls%mx1
            IF (ALL(ls%ig(k1, k2, :) == 0)) CYCLE
            DO q = 1, nq
               g(q) = wa(q)*priv_filmcoef(fs, fcell, fvac, pw_in, vac_in, k1, k2, za(q))
               IF (wa(q) < 1.0) g(q) = g(q) + (1.0 - wa(q))*priv_filmcoef(fs, fcell, fvac, pw_in, vac_in, k1, k2, zb(q))
            END DO
            DO k3l = -ls%mx3, ls%mx3
               is = ls%ig(k1, k2, k3l)
               IF (is == 0) CYCLE
               c = SUM(g*eq(:, k3l))
               pw_out(is) = pw_out(is) + c*CONJG(ls%rgphs(k1, k2, k3l))
            END DO
         END DO
      END DO
      pw_out(:ls%ng3) = pw_out(:ls%ng3)/ls%nstr(:ls%ng3)
   END SUBROUTINE gfprep_pwconv_film

   SUBROUTINE priv_film_nodes(z1, zc, d, dt, zq, wq, za, zb, wa)
      !Quadrature nodes zq (relative to zc) and weights wq on [-dt/2,dt/2].
      !The layer function at zq is wa*f(za)+(1-wa)*f(zb) with absolute film
      !positions za, zb (zb: periodic image used in the buffer). The panels
      !end at the kinks of f (|z|=z1) and of the blending (|z|=d/2).
      REAL, INTENT(IN)               :: z1, zc, d, dt
      REAL, ALLOCATABLE, INTENT(OUT) :: zq(:), wq(:), za(:), zb(:), wa(:)

      REAL, ALLOCATABLE :: xg(:), wg(:)
      REAL    :: bp(14), a, b, h, z, u
      INTEGER :: nbp, i, j, np, ip, n, side, is

      nbp = 0
      CALL priv_addbp(-0.5*dt); CALL priv_addbp(-0.5*d); CALL priv_addbp(0.5*d); CALL priv_addbp(0.5*dt)
      DO is = -1, 1, 2
         CALL priv_addbp(is*z1 - zc)
         DO side = -1, 1, 2
            !the image zq-side*dt crosses |z|=z1 inside the buffer
            z = is*z1 - zc + side*dt
            IF (ABS(z) > 0.5*d) CALL priv_addbp(z)
         END DO
      END DO
      !sort
      DO i = 2, nbp
         DO j = i, 2, -1
            IF (bp(j) >= bp(j - 1)) EXIT
            bp(j - 1:j) = [bp(j), bp(j - 1)]
         END DO
      END DO
      CALL priv_gauleg(nglf, xg, wg)
      n = 0
      DO i = 1, nbp - 1
         IF (bp(i + 1) - bp(i) > 1.0E-10) n = n + nglf*CEILING((bp(i + 1) - bp(i))/hpanel)
      END DO
      ALLOCATE (zq(n), wq(n), za(n), zb(n), wa(n))
      n = 0
      DO i = 1, nbp - 1
         IF (bp(i + 1) - bp(i) <= 1.0E-10) CYCLE
         np = CEILING((bp(i + 1) - bp(i))/hpanel)
         h = (bp(i + 1) - bp(i))/np
         DO ip = 1, np
            a = bp(i) + (ip - 1)*h
            b = a + h
            DO j = 1, nglf
               n = n + 1
               zq(n) = 0.5*(a + b) + 0.5*h*xg(j)
               wq(n) = 0.5*h*wg(j)
            END DO
         END DO
      END DO
      DO n = 1, SIZE(zq)
         za(n) = zc + zq(n)
         zb(n) = za(n)
         wa(n) = 1.0
         IF (ABS(zq(n)) > 0.5*d) THEN
            u = (ABS(zq(n)) - 0.5*d)/(0.5*(dt - d))
            wa(n) = 1.0 - 0.5*SIN(0.5*pi_const*u)**2
            zb(n) = za(n) - SIGN(dt, zq(n))
         END IF
      END DO
   CONTAINS
      SUBROUTINE priv_addbp(z)
         REAL, INTENT(IN) :: z
         IF (z < -0.5*dt - 1.0E-10 .OR. z > 0.5*dt + 1.0E-10) RETURN
         nbp = nbp + 1
         bp(nbp) = MAX(-0.5*dt, MIN(0.5*dt, z))
      END SUBROUTINE priv_addbp
   END SUBROUTINE priv_film_nodes

   COMPLEX FUNCTION priv_filmcoef(fs, fcell, fvac, pw, vac, k1, k2, z)
      !Coefficient of the single in-plane G=(k1,k2) of the film function at
      !the cartesian height z: plane waves for |z|<=z1, vacuum beyond. The
      !warping is zero beyond nmzxy, the planar average is constant beyond nmz.
      TYPE(t_stars), INTENT(IN)  :: fs
      TYPE(t_cell), INTENT(IN)   :: fcell
      TYPE(t_vacuum), INTENT(IN) :: fvac
      COMPLEX, INTENT(IN)        :: pw(:), vac(:, :, :)
      INTEGER, INTENT(IN)        :: k1, k2
      REAL, INTENT(IN)           :: z

      INTEGER :: k3, is, s2, ivac, nz
      REAL    :: x

      priv_filmcoef = CMPLX(0.0, 0.0)
      IF (ABS(k1) > fs%mx1 .OR. ABS(k2) > fs%mx2) RETURN
      IF (ABS(z) <= fcell%z1) THEN
         DO k3 = -fs%mx3, fs%mx3
            is = fs%ig(k1, k2, k3)
            IF (is == 0) CYCLE
            priv_filmcoef = priv_filmcoef + pw(is)*fs%rgphs(k1, k2, k3)*EXP(CMPLX(0.0, k3*fcell%bmat(3, 3)*z))
         END DO
         RETURN
      END IF
      s2 = fs%i2g(k1, k2)
      IF (s2 == 0) RETURN
      ivac = MERGE(1, 2, z > 0.0)
      nz = MERGE(fvac%nmz, fvac%nmzxy, s2 == 1)
      x = (ABS(z) - fcell%z1)/fvac%delz + 1.0
      IF (x > nz) THEN
         IF (s2 == 1) priv_filmcoef = vac(fvac%nmz, 1, ivac)
         RETURN
      END IF
      priv_filmcoef = priv_interp(vac(:nz, s2, ivac), x)*fs%r2gphs(k1, k2)
   END FUNCTION priv_filmcoef

   COMPLEX FUNCTION priv_interp(f, x)
      !cubic Lagrange interpolation of f at the (real) grid index x
      COMPLEX, INTENT(IN) :: f(:)
      REAL, INTENT(IN)    :: x
      INTEGER :: i0, i, j
      REAL    :: l
      i0 = MIN(MAX(FLOOR(x) - 1, 1), SIZE(f) - 3)
      priv_interp = CMPLX(0.0, 0.0)
      DO i = i0, i0 + 3
         l = 1.0
         DO j = i0, i0 + 3
            IF (j /= i) l = l*(x - j)/(i - j)
         END DO
         priv_interp = priv_interp + l*f(i)
      END DO
   END FUNCTION priv_interp

   SUBROUTINE gfprep_vacemb(fs, fcell, fvac, side, zs, pw, vac, delz, vz)
      !Planar average of the film potential from the vacuum embedding plane
      !at |z|=zs (side 1: z=-zs, side 2: z=+zs) outwards on the vacuum grid,
      !up to the end of the FLEUR vacuum (vz(1) is the value at the plane)
      TYPE(t_stars), INTENT(IN)      :: fs
      TYPE(t_cell), INTENT(IN)       :: fcell
      TYPE(t_vacuum), INTENT(IN)     :: fvac
      INTEGER, INTENT(IN)            :: side
      REAL, INTENT(IN)               :: zs
      COMPLEX, INTENT(IN)            :: pw(:, :), vac(:, :, :, :)   !pw(ng3,js), vac(nmzd,ng2,2,js)
      REAL, INTENT(OUT)              :: delz
      REAL, ALLOCATABLE, INTENT(OUT) :: vz(:, :)

      INTEGER :: n, i, js
      REAL    :: zend, z

      delz = fvac%delz
      zend = fcell%z1 + (fvac%nmz - 1)*fvac%delz
      IF (zs >= zend - 3*delz) CALL juDFT_error("Vacuum embedding plane beyond the FLEUR vacuum grid", &
                                                hint="decrease zvac", calledby="gfprep_vacemb")
      n = FLOOR((zend - zs)/delz + 1.0E-6) + 1
      ALLOCATE (vz(n, SIZE(pw, 2)))
      DO js = 1, SIZE(pw, 2)
         DO i = 1, n
            z = (zs + (i - 1)*delz)*(2*side - 3)
            vz(i, js) = REAL(priv_filmcoef(fs, fcell, fvac, pw(:, js), vac(:, :, :, js), 0, 0, z))
         END DO
      END DO
   END SUBROUTINE gfprep_vacemb

   SUBROUTINE gfprep_mtconv(b, l, typmap, atmap, mt_in, mt_out)
      !mt_in/mt_out: (r,lh,type,spin) of bulk and layer
      TYPE(t_gfprep_struct), INTENT(IN) :: b, l
      INTEGER, INTENT(IN)               :: typmap(:), atmap(:)
      REAL, INTENT(IN)                  :: mt_in(:, 0:, :, :)
      REAL, INTENT(OUT)                 :: mt_out(:, 0:, :, :)

      REAL, ALLOCATABLE :: tm(:, :)
      INTEGER :: t, n, ab, al, nsb, nsl, js, jri
      REAL    :: rb(3, 3), rl(3, 3), loss

      mt_out = 0.0
      DO t = 1, l%fi%atoms%ntype
         n = typmap(t)
         al = l%fi%atoms%firstAtom(t)
         ab = atmap(al)
         nsb = b%fi%sym%ntypsy(ab)
         nsl = l%fi%sym%ntypsy(al)
         rb = priv_cartrot(b%fi%cell, b%fi%sym%mrot(:, :, b%fi%sym%ngopr(ab)))
         rl = priv_cartrot(l%fi%cell, l%fi%sym%mrot(:, :, l%fi%sym%ngopr(al)))
         CALL priv_transform(b%sphhar, nsb, l%sphhar, nsl, MATMUL(rb, TRANSPOSE(rl)), tm, loss)
         jri = l%fi%atoms%jri(t)
         IF (jri /= b%fi%atoms%jri(n)) CALL juDFT_error("BUG: radial meshes differ", calledby="gfprep_mtconv")
         DO js = 1, SIZE(mt_out, 4)
            mt_out(:jri, 0:l%sphhar%nlh(nsl), t, js) = MATMUL(mt_in(:jri, 0:b%sphhar%nlh(nsb), n, js), TRANSPOSE(tm))
         END DO
         IF (loss > 1.0E-8) THEN
            WRITE (oUnit, "(a,i0,a,es10.2)") "gfprep: WARNING: MT expansion of layer type ", t, &
               " not complete, lost norm fraction:", loss
            CALL juDFT_warn("Lattice harmonics of a layer atom do not span the bulk function", &
                            calledby="gfprep_mtconv")
         END IF
      END DO
   END SUBROUTINE gfprep_mtconv

   SUBROUTINE priv_transform(sb, nsb, sl, nsl, rot, tm, loss)
      !tm(lh',lh) = <K'_lh'(y)|K_lh(rot y)>/<K'_lh'|K'_lh'> by quadrature;
      !loss is the largest relative norm of a bulk harmonic not
      !representable by the layer harmonics
      TYPE(t_sphhar), INTENT(IN)      :: sb, sl
      INTEGER, INTENT(IN)             :: nsb, nsl
      REAL, INTENT(IN)                :: rot(3, 3)
      REAL, ALLOCATABLE, INTENT(OUT)  :: tm(:, :)
      REAL, INTENT(OUT)               :: loss

      REAL, ALLOCATABLE    :: xg(:), wg(:), kb(:), kl(:), nb(:), nl(:)
      COMPLEX, ALLOCATABLE :: ylm(:)
      INTEGER :: lq, nth, nph, it, ip, lh
      REAL    :: y(3), w, phi, st
      REAL, ALLOCATABLE :: proj(:)

      lq = MAX(MAXVAL(sb%llh(0:sb%nlh(nsb), nsb)), MAXVAL(sl%llh(0:sl%nlh(nsl), nsl)))
      nth = lq + 2
      nph = 2*lq + 3
      CALL priv_gauleg(nth, xg, wg)
      ALLOCATE (tm(0:sl%nlh(nsl), 0:sb%nlh(nsb)), kb(0:sb%nlh(nsb)), kl(0:sl%nlh(nsl)))
      ALLOCATE (nb(0:sb%nlh(nsb)), nl(0:sl%nlh(nsl)), ylm((lq + 1)**2))
      tm = 0.0; nb = 0.0; nl = 0.0
      DO it = 1, nth
         st = SQRT(MAX(0.0, 1.0 - xg(it)**2))
         DO ip = 1, nph
            phi = tpi_const*(ip - 1)/nph
            y = [st*COS(phi), st*SIN(phi), xg(it)]
            w = wg(it)*tpi_const/nph
            CALL ylm4(lq, y, ylm)
            CALL priv_latharm(sl, nsl, ylm, kl)
            CALL ylm4(lq, MATMUL(rot, y), ylm)
            CALL priv_latharm(sb, nsb, ylm, kb)
            DO lh = 0, sl%nlh(nsl)
               tm(lh, :) = tm(lh, :) + w*kl(lh)*kb
            END DO
            nb = nb + w*kb**2
            nl = nl + w*kl**2
         END DO
      END DO
      !harmonics of different l are orthogonal: the quadrature leaves only rounding
      !noise (~1e-16) in these entries, which would pollute the high-l components
      !of the density (the Poisson kernel r^(-l-1) amplifies them enormously)
      DO lh = 0, sl%nlh(nsl)
         DO it = 0, sb%nlh(nsb)
            IF (sl%llh(lh, nsl) /= sb%llh(it, nsb)) tm(lh, it) = 0.0
         END DO
      END DO
      !norm of the bulk harmonics captured by the layer harmonics
      ALLOCATE (proj(0:sb%nlh(nsb)))
      DO lh = 0, sb%nlh(nsb)
         proj(lh) = SUM(tm(:, lh)**2/nl)
      END DO
      loss = MAXVAL(ABS(1.0 - proj/nb))
      DO lh = 0, sl%nlh(nsl)
         tm(lh, :) = tm(lh, :)/nl(lh)
      END DO
   END SUBROUTINE priv_transform

   SUBROUTINE priv_latharm(s, ns, ylm, k)
      !real lattice harmonics of symmetry type ns from the Ylm
      TYPE(t_sphhar), INTENT(IN) :: s
      INTEGER, INTENT(IN)        :: ns
      COMPLEX, INTENT(IN)        :: ylm(:)
      REAL, INTENT(OUT)          :: k(0:)
      INTEGER :: lh, mem, l
      COMPLEX :: c
      DO lh = 0, s%nlh(ns)
         l = s%llh(lh, ns)
         c = CMPLX(0.0, 0.0)
         DO mem = 1, s%nmem(lh, ns)
            c = c + s%clnu(mem, lh, ns)*ylm(l*(l + 1) + s%mlh(mem, lh, ns) + 1)
         END DO
         k(lh) = REAL(c)
      END DO
   END SUBROUTINE priv_latharm

   FUNCTION priv_cartrot(cell, mrot) RESULT(r)
      !cartesian matrix of a rotation given in internal coordinates
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER, INTENT(IN)      :: mrot(3, 3)
      REAL :: r(3, 3)
      r = MATMUL(cell%amat, MATMUL(REAL(mrot), cell%bmat))/tpi_const
   END FUNCTION priv_cartrot

   SUBROUTINE priv_gauleg(n, x, w)
      !Gauss-Legendre nodes and weights on [-1,1]
      INTEGER, INTENT(IN)            :: n
      REAL, ALLOCATABLE, INTENT(OUT) :: x(:), w(:)
      INTEGER :: i, j, it
      REAL    :: z, p0, p1, p2, dp
      ALLOCATE (x(n), w(n))
      DO i = 1, (n + 1)/2
         z = COS(pi_const*(i - 0.25)/(n + 0.5))
         DO it = 1, 100
            p0 = 1.0; p1 = z
            DO j = 2, n
               p2 = ((2*j - 1)*z*p1 - (j - 1)*p0)/j
               p0 = p1; p1 = p2
            END DO
            IF (n == 1) THEN
               p1 = z; p0 = 1.0
            END IF
            dp = n*(z*p1 - p0)/(z*z - 1.0)
            z = z - p1/dp
            IF (ABS(p1/dp) < 1.0E-15) EXIT
         END DO
         x(i) = -z; x(n + 1 - i) = z
         w(i) = 2.0/((1.0 - z*z)*dp*dp); w(n + 1 - i) = w(i)
      END DO
   END SUBROUTINE priv_gauleg

   SUBROUTINE gfprep_check(name, b, l, ilay, zc, c, typmap, atmap, bpw, lpw, bmt, lmt, bvac)
      !Compares the bulk (film) and the converted layer function at test
      !points inside the physical region |z|<c/2 (interstitial and, for a
      !film, vacuum) and on spheres of half the MT radius (muffin-tins) for
      !the first spin
      CHARACTER(LEN=*), INTENT(IN)      :: name
      TYPE(t_gfprep_struct), INTENT(IN) :: b, l
      INTEGER, INTENT(IN)               :: ilay, typmap(:), atmap(:)
      REAL, INTENT(IN)                  :: zc, c
      COMPLEX, INTENT(IN)               :: bpw(:), lpw(:)
      REAL, INTENT(IN)                  :: bmt(:, 0:, :, :), lmt(:, 0:, :, :)
      COMPLEX, INTENT(IN), OPTIONAL     :: bvac(:, :, :)   !film vacuum (nmzd,ng2,2)

      INTEGER :: i1, i2, i3, t, al, ab, n, j, ip
      REAL    :: r(3), fb, fl, dmax, fmax, y(3), rb(3, 3), rl(3, 3)
      REAL, PARAMETER :: dirs(3, 3) = RESHAPE([0.3, 0.5, 0.81, -0.7, 0.2, 0.68, 0.1, -0.9, -0.42], [3, 3])

      !interstitial points only
      dmax = 0.0; fmax = 0.0
      DO i3 = -4, 4
         DO i2 = 0, 5
            DO i1 = 0, 5
               r = MATMUL(l%fi%cell%amat(:, 1:2), [i1/6.0, i2/6.0])
               r(3) = 0.999*0.5*c*i3/4.0 + zc
               IF (priv_in_mt(b, r)) CYCLE
               IF (PRESENT(bvac)) THEN
                  !the plane-wave and vacuum expansions of FLEUR do not match
                  !exactly at |z|=z1
                  IF (ABS(ABS(r(3)) - b%fi%cell%z1) < 0.1) CYCLE
                  fb = priv_filmeval(b, bpw, bvac, r)
               ELSE
                  fb = priv_pweval(b%stars, b%fi%cell, bpw, r)
               END IF
               r(3) = r(3) - zc
               fl = priv_pweval(l%stars, l%fi%cell, lpw, r)
               dmax = MAX(dmax, ABS(fb - fl)); fmax = MAX(fmax, ABS(fb))
            END DO
         END DO
      END DO
      WRITE (oUnit, "(a,a,a,i0,a,2es11.3)") "gfprep check ", name, " layer ", ilay, &
         " interstitial: max|f|, max deviation:", fmax, dmax
      WRITE (*, "(a,a,a,i0,a,es11.3)") "gfprep check ", name, " layer ", ilay, &
         " interstitial: relative deviation", dmax/MAX(fmax, 1.0E-12)

      dmax = 0.0; fmax = 0.0
      DO t = 1, l%fi%atoms%ntype
         n = typmap(t)
         al = l%fi%atoms%firstAtom(t)
         ab = atmap(al)
         rb = priv_cartrot(b%fi%cell, b%fi%sym%mrot(:, :, b%fi%sym%ngopr(ab)))
         rl = priv_cartrot(l%fi%cell, l%fi%sym%mrot(:, :, l%fi%sym%ngopr(al)))
         j = l%fi%atoms%jri(t)/2
         DO ip = 1, 3
            y = dirs(:, ip)/NORM2(dirs(:, ip))
            fb = priv_mteval(b%sphhar, b%fi%sym%ntypsy(ab), bmt(j, :, n, 1), MATMUL(rb, y))
            fl = priv_mteval(l%sphhar, l%fi%sym%ntypsy(al), lmt(j, :, t, 1), MATMUL(rl, y))
            dmax = MAX(dmax, ABS(fb - fl)); fmax = MAX(fmax, ABS(fb))
         END DO
      END DO
      WRITE (oUnit, "(a,a,a,i0,a,2es11.3)") "gfprep check ", name, " layer ", ilay, &
         " muffin-tins: max|f|, max deviation:", fmax, dmax
      WRITE (*, "(a,a,a,i0,a,es11.3)") "gfprep check ", name, " layer ", ilay, &
         " muffin-tins: relative deviation", dmax/MAX(fmax, 1.0E-12)
   END SUBROUTINE gfprep_check

   LOGICAL FUNCTION priv_in_mt(b, r)
      TYPE(t_gfprep_struct), INTENT(IN) :: b
      REAL, INTENT(IN)                  :: r(3)
      INTEGER :: na, i1, i2, i3, n3
      priv_in_mt = .FALSE.
      n3 = MERGE(0, 2, b%fi%input%film)
      DO na = 1, b%fi%atoms%nat
         DO i3 = -n3, n3
            DO i2 = -1, 1
               DO i1 = -1, 1
                  IF (NORM2(r - b%fi%atoms%pos(:, na) - MATMUL(b%fi%cell%amat, REAL([i1, i2, i3]))) &
                      < b%fi%atoms%rmt(b%fi%atoms%itype(na))) priv_in_mt = .TRUE.
               END DO
            END DO
         END DO
      END DO
   END FUNCTION priv_in_mt

   REAL FUNCTION priv_pweval(s, cell, pw, r)
      TYPE(t_stars), INTENT(IN) :: s
      TYPE(t_cell), INTENT(IN)  :: cell
      COMPLEX, INTENT(IN)       :: pw(:)
      REAL, INTENT(IN)          :: r(3)
      INTEGER :: k1, k2, k3
      COMPLEX :: f
      f = CMPLX(0.0, 0.0)
      DO k3 = -s%mx3, s%mx3
         DO k2 = -s%mx2, s%mx2
            DO k1 = -s%mx1, s%mx1
               IF (s%ig(k1, k2, k3) == 0) CYCLE
               f = f + pw(s%ig(k1, k2, k3))*s%rgphs(k1, k2, k3)* &
                   EXP(CMPLX(0.0, DOT_PRODUCT(MATMUL(REAL([k1, k2, k3]), cell%bmat), r)))
            END DO
         END DO
      END DO
      priv_pweval = REAL(f)
   END FUNCTION priv_pweval

   REAL FUNCTION priv_filmeval(b, pw, vac, r)
      !film function at r: plane waves for |z|<=z1, vacuum beyond
      TYPE(t_gfprep_struct), INTENT(IN) :: b
      COMPLEX, INTENT(IN)               :: pw(:), vac(:, :, :)
      REAL, INTENT(IN)                  :: r(3)
      INTEGER :: k1, k2
      COMPLEX :: f
      IF (ABS(r(3)) <= b%fi%cell%z1) THEN
         priv_filmeval = priv_pweval(b%stars, b%fi%cell, pw, r)
         RETURN
      END IF
      f = CMPLX(0.0, 0.0)
      DO k2 = -b%stars%mx2, b%stars%mx2
         DO k1 = -b%stars%mx1, b%stars%mx1
            IF (b%stars%i2g(k1, k2) == 0) CYCLE
            f = f + priv_filmcoef(b%stars, b%fi%cell, b%fi%vacuum, pw, vac, k1, k2, r(3))* &
                EXP(CMPLX(0.0, DOT_PRODUCT(MATMUL(REAL([k1, k2]), b%fi%cell%bmat(1:2, 1:2)), r(1:2))))
         END DO
      END DO
      priv_filmeval = REAL(f)
   END FUNCTION priv_filmeval

   REAL FUNCTION priv_mteval(s, ns, f, y)
      TYPE(t_sphhar), INTENT(IN) :: s
      INTEGER, INTENT(IN)        :: ns
      REAL, INTENT(IN)           :: f(0:), y(3)
      COMPLEX :: ylm((MAXVAL(s%llh(0:s%nlh(ns), ns)) + 1)**2)
      REAL    :: k(0:s%nlh(ns))
      CALL ylm4(MAXVAL(s%llh(0:s%nlh(ns), ns)), y, ylm)
      CALL priv_latharm(s, ns, ylm, k)
      priv_mteval = DOT_PRODUCT(f(0:s%nlh(ns)), k)
   END FUNCTION priv_mteval

END MODULE m_gfprep_potconv
