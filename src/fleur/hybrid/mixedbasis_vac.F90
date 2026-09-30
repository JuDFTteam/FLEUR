!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>VAC block of the mixed product basis: M = U_i(z) exp(i(q+G||).r||)/sqrt(A) for |z| > z1.
!>U_i are orthonormalised products of the vacuum functions u, udot, taken on a ladder of
!>individual in-plane momenta (the product of k+G and k+q+G' decays with both), per length |q+G||.
MODULE m_mixedbasis_vac
   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_vac_rows, ONLY: NVAC_MPB, vac_src
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: gen_vac_basis

   ! rungs of the momentum ladder, uniform in kappa = sqrt(h^2 + 2(V_inf - E))
   INTEGER, PARAMETER :: NH_LADDER = 6

CONTAINS

   SUBROUTINE vac_rung(vacuum, vz, vzero, evac, h, u, ud)
      USE m_vacuz
      USE m_vacudz
      IMPLICIT NONE
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      REAL, INTENT(IN)           :: vz(:), vzero, evac, h
      REAL, INTENT(OUT)          :: u(:), ud(:)

      REAL :: ev, uz, duz, udz, dudz, ddnv, wronk, scale

      ev = evac - 0.5*h**2
      wronk = 2.0
      CALL vacuz(ev, vz, vzero, vacuum%nmz, vacuum%delz, uz, duz, u)
      CALL vacudz(ev, vz, vzero, vacuum%nmz, vacuum%delz, udz, dudz, ddnv, ud, duz, u)
      scale = wronk/(udz*duz - dudz*uz)
      ud = scale*ud
   END SUBROUTINE vac_rung

   SUBROUTINE vac_mom_ladder(rkmax, cc, nh, h)
      IMPLICIT NONE
      REAL, INTENT(IN)               :: rkmax, cc
      INTEGER, INTENT(OUT)           :: nh
      REAL, ALLOCATABLE, INTENT(OUT) :: h(:)

      INTEGER :: i
      REAL    :: k0, k1, kk

      nh = NH_LADDER
      ALLOCATE (h(nh))
      k0 = SQRT(cc)
      k1 = SQRT(rkmax**2 + cc)
      DO i = 1, nh
         kk = k0 + (k1 - k0)*(i - 1)/REAL(nh - 1)
         h(i) = SQRT(MAX(0.0, kk**2 - cc))
      END DO
   END SUBROUTINE vac_mom_ladder

   SUBROUTINE gen_vac_basis(vacuum, cell, kpts, input, mpinp, mpdata, enpara, v, fmpi)
      USE m_vacuz
      USE m_vacudz
      IMPLICIT NONE
      TYPE(t_vacuum), INTENT(IN)        :: vacuum
      TYPE(t_cell), INTENT(IN)          :: cell
      TYPE(t_kpts), INTENT(IN)          :: kpts
      TYPE(t_input), INTENT(IN)         :: input
      TYPE(t_mpinp), INTENT(IN)         :: mpinp
      TYPE(t_mpdata), INTENT(INOUT)     :: mpdata
      TYPE(t_enpara), INTENT(IN)        :: enpara
      TYPE(t_potden), INTENT(IN)        :: v
      TYPE(t_mpi), INTENT(IN)           :: fmpi

      INTEGER :: ivac, ilen, nlen, nz, i, j, nraw, nrawmx, nkeep
      INTEGER :: nh, ih1, ih2, nkmin, nkmax, nhmx, nhg
      REAL    :: g, vzero
      REAL    :: tol, cc, hnew

      REAL, ALLOCATABLE :: u(:), ud(:)
      REAL, ALLOCATABLE :: raw(:, :), evsc(:, :)
      REAL, ALLOCATABLE :: olap(:, :), eig(:), eigv(:, :)
      REAL, ALLOCATABLE :: vz(:)
      REAL, ALLOCATABLE :: hlad(:)
      REAL, ALLOCATABLE :: hgrp(:)
      REAL, ALLOCATABLE :: ulad(:, :), udlad(:, :)

      CALL timestart("vacuum mixed basis")

      tol = mpinp%linear_dep_tol

      CALL gen_vac_gvec(cell, kpts, mpinp, mpdata)
      nlen = SIZE(mpdata%glen_vac)

      ALLOCATE (vz(vacuum%nmz))
      vz = REAL(v%vac(:vacuum%nmz, 1, 1, 1))

      CALL set_vac_zrange(vacuum, mpinp, enpara, vz, mpdata)
      nz = mpdata%nmz_vac

      cc = 2.0*(vz(vacuum%nmz) - enpara%evac(1, 1))
      CALL vac_mom_ladder(input%rkmax, MAX(cc, 1e-8), nh, hlad)

      nhmx = nh + 2

      nrawmx = 3*nhmx + 4*(nhmx*(nhmx - 1))/2

      ALLOCATE (mpdata%num_zbasfn_vac(nlen, NVAC_MPB), source=0)
      ALLOCATE (mpdata%zbasfn_vac(nz, nrawmx + 1, nlen, NVAC_MPB), source=0.0)
      ALLOCATE (u(vacuum%nmz), ud(vacuum%nmz))
      ALLOCATE (ulad(vacuum%nmz, nhmx), udlad(vacuum%nmz, nhmx))
      ALLOCATE (hgrp(nhmx))

      nkmin = HUGE(nkmin); nkmax = 0
      DO ivac = 1, NVAC_MPB
         vz = REAL(v%vac(:vacuum%nmz, 1, vac_src(vacuum%nvac, ivac), 1))
         vzero = vz(vacuum%nmz)

         ! identical vacua share the same functions (needed by the invs path)
         IF (ivac > 1) THEN
            IF (MAXVAL(ABS(vz - REAL(v%vac(:vacuum%nmz, 1, 1, 1)))) < 1e-10) THEN
               mpdata%num_zbasfn_vac(:, ivac) = mpdata%num_zbasfn_vac(:, 1)
               mpdata%zbasfn_vac(:, :, :, ivac) = mpdata%zbasfn_vac(:, :, :, 1)
               CYCLE
            END IF
         END IF

         DO i = 1, nh
            CALL vac_rung(vacuum, vz, vzero, enpara%evac(vac_src(vacuum%nvac, ivac), 1), hlad(i), u, ud)
            ulad(:, i) = u
            udlad(:, i) = ud
         END DO
         hgrp(:nh) = hlad

         DO ilen = 1, nlen
            g = mpdata%glen_vac(ilen)

            ! add the slowest admissible momenta g/2 and g of this group
            nhg = nh
            DO i = 1, 2
               hnew = MERGE(0.5*g, g, i == 1)
               IF (ANY(ABS(hgrp(:nhg) - hnew) < 1e-6)) CYCLE
               nhg = nhg + 1
               hgrp(nhg) = hnew
               CALL vac_rung(vacuum, vz, vzero, enpara%evac(vac_src(vacuum%nvac, ivac), 1), hnew, u, ud)
               ulad(:, nhg) = u
               udlad(:, nhg) = ud
            END DO

            ALLOCATE (raw(nz, nrawmx))
            ! pairs allowed by the triangle inequality |h_a - h_b| <= g <= h_a + h_b
            nraw = 0
            DO ih1 = 1, nhg
               DO ih2 = ih1, nhg
                  IF (ABS(hgrp(ih1) - hgrp(ih2)) > g + 1e-10) CYCLE
                  IF (g > hgrp(ih1) + hgrp(ih2) + 1e-10) CYCLE
                  nraw = nraw + 1
                  raw(:nz, nraw) = ulad(:nz, ih1)*ulad(:nz, ih2)
                  nraw = nraw + 1
                  raw(:nz, nraw) = ulad(:nz, ih1)*udlad(:nz, ih2)
                  IF (ih1 /= ih2) THEN
                     nraw = nraw + 1
                     raw(:nz, nraw) = udlad(:nz, ih1)*ulad(:nz, ih2)
                  END IF
                  nraw = nraw + 1
                  raw(:nz, nraw) = udlad(:nz, ih1)*udlad(:nz, ih2)
               END DO
            END DO
            IF (nraw == 0) CALL juDFT_error( &
               "vacuum mixed basis: no admissible momentum pair for a length group", &
               calledby="gen_vac_basis", &
               hint="the ladder does not reach this transferred |q+G|||")

            ALLOCATE (olap(nraw, nraw))
            DO i = 1, nraw
               DO j = 1, i
                  olap(i, j) = z_dot(raw(:, i), raw(:, j), nz, vacuum%delz)
                  olap(j, i) = olap(i, j)
               END DO
            END DO

            CALL diag_sym(olap, eig, eigv)

            nkeep = COUNT(eig > tol)
            IF (nkeep == 0) THEN
               CALL juDFT_error("vacuum mixed basis: all product functions filtered out", &
                                calledby="gen_vac_basis", &
                                hint="check the vacuum potential and prodBasis/@tolerance")
            END IF

            ALLOCATE (evsc(nraw, nkeep))
            DO i = 1, nkeep
               j = nraw - i + 1
               evsc(:, i) = eigv(:, j)/SQRT(eig(j))
            END DO
            CALL dgemm('N', 'N', nz, nkeep, nraw, 1.0, raw, SIZE(raw, 1), &
                       evsc, nraw, 0.0, mpdata%zbasfn_vac(:, :, ilen, ivac), &
                       SIZE(mpdata%zbasfn_vac, 1))
            DEALLOCATE (evsc)
            mpdata%num_zbasfn_vac(ilen, ivac) = nkeep
            nkmin = MIN(nkmin, nkeep); nkmax = MAX(nkmax, nkeep)

            DEALLOCATE (raw, olap, eig, eigv)
         END DO

         ! one function per group carries the exponential moment, at g = 0 a second the first moment;
         ! only these couple outside their own vacuum
         CALL add_vac_moment_fun(mpdata, vacuum, ivac, tol)

         DO ilen = 1, nlen
            CALL concentrate_exp_moment(mpdata, vacuum, ilen, ivac, &
                                        mpdata%glen_vac(ilen), fmpi)
            IF (mpdata%glen_vac(ilen) < 1e-12) &
               CALL concentrate_first_moment(mpdata, vacuum, ilen, ivac)
         END DO
      END DO

      IF (fmpi%irank == 0) WRITE (oUnit, '(A,I0,A,I0,A,I0)') &
         "  vacuum product basis: z-functions per group ", nkmin, " .. ", nkmax, ", total ", &
         SUM(mpdata%num_zbasfn_vac)
      CALL check_exp_moment(mpdata, vacuum, i)
      IF (i /= 0) CALL juDFT_error("exponential moment is not concentrated in one &
                                   &function in "//int2str(i)//" groups", &
                                   calledby="gen_vac_basis")

      CALL timestop("vacuum mixed basis")
   END SUBROUTINE gen_vac_basis

   REAL FUNCTION z_dot(f, h, nz, delz)
      USE m_coulomb_vac, ONLY: vac_zint
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:), h(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz

      z_dot = vac_zint(f(:nz)*h(:nz), nz, delz)
   END FUNCTION z_dot

   !>In-plane vectors of the VAC block (cutoff as for IR) and their length groups.
   SUBROUTINE gen_vac_gvec(cell, kpts, mpinp, mpdata)
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN)       :: cell
      TYPE(t_kpts), INTENT(IN)       :: kpts
      TYPE(t_mpinp), INTENT(IN)      :: mpinp
      TYPE(t_mpdata), INTENT(INOUT)  :: mpdata

      INTEGER :: ik, ig, i, j, n, nmax, cnt
      INTEGER :: gv(2)
      REAL    :: q(3), qn
      INTEGER, ALLOCATABLE :: glist(:, :)
      REAL, ALLOCATABLE    :: nrm(:)

      nmax = mpdata%num_gpts()
      ALLOCATE (glist(2, nmax), source=0)
      n = 0
      DO ig = 1, nmax
         gv = mpdata%g(1:2, ig)
         DO i = 1, n
            IF (ALL(glist(:, i) == gv)) GOTO 10
         END DO
         n = n + 1
         glist(:, n) = gv
10       CONTINUE
      END DO
      ALLOCATE (mpdata%g_vac(2, n))
      mpdata%g_vac = glist(:, :n)
      DEALLOCATE (glist)

      ALLOCATE (mpdata%n_g_vac(kpts%nkptf), source=0)
      ALLOCATE (mpdata%gptm_ptr_vac(n, kpts%nkptf), source=0)
      ALLOCATE (mpdata%glen_ptr_vac(n, kpts%nkptf), source=0)
      ALLOCATE (nrm(n*kpts%nkptf), source=0.0)
      cnt = 0

      DO ik = 1, kpts%nkptf
         DO i = 1, n
            q = MATMUL([kpts%bkf(1, ik) + mpdata%g_vac(1, i), &
                        kpts%bkf(2, ik) + mpdata%g_vac(2, i), 0.0], cell%bmat)
            qn = SQRT(q(1)**2 + q(2)**2)
            IF (qn > mpinp%g_cutoff) CYCLE
            mpdata%n_g_vac(ik) = mpdata%n_g_vac(ik) + 1
            mpdata%gptm_ptr_vac(mpdata%n_g_vac(ik), ik) = i
            DO j = 1, cnt
               IF (ABS(nrm(j) - qn) < 1e-5*MAX(1.0, qn)) THEN
                  mpdata%glen_ptr_vac(mpdata%n_g_vac(ik), ik) = j
                  GOTO 20
               END IF
            END DO
            cnt = cnt + 1
            nrm(cnt) = qn
            mpdata%glen_ptr_vac(mpdata%n_g_vac(ik), ik) = cnt
20          CONTINUE
         END DO
      END DO

      ALLOCATE (mpdata%glen_vac(cnt))
      mpdata%glen_vac = nrm(:cnt)
   END SUBROUTINE gen_vac_gvec

   !>z range: until the slowest product drops below the linear-dependence tolerance.
   SUBROUTINE set_vac_zrange(vacuum, mpinp, enpara, vz, mpdata)
      USE m_vacuz
      IMPLICIT NONE
      TYPE(t_vacuum), INTENT(IN)     :: vacuum
      TYPE(t_mpinp), INTENT(IN)      :: mpinp
      TYPE(t_enpara), INTENT(IN)     :: enpara
      REAL, INTENT(IN)               :: vz(:)
      TYPE(t_mpdata), INTENT(INOUT)  :: mpdata

      REAL    :: u(vacuum%nmz), uz, duz, ev, vzero, ref
      INTEGER :: iz, nz

      ev = enpara%evac(1, 1)
      vzero = vz(vacuum%nmz)
      CALL vacuz(ev, vz, vzero, vacuum%nmz, vacuum%delz, uz, duz, u)

      ref = u(1)**2
      nz = vacuum%nmz
      IF (ref > 0.0) THEN
         DO iz = 1, vacuum%nmz
            IF (u(iz)**2 < mpinp%linear_dep_tol*ref) THEN
               nz = iz
               EXIT
            END IF
         END DO
      END IF
      mpdata%nmz_vac = MAX(MIN(nz, vacuum%nmz), 9)

   END SUBROUTINE set_vac_zrange

   !>Add the constant to the g = 0 group, as add_l0_fun does for the MT basis.
   SUBROUTINE add_vac_moment_fun(mpdata, vacuum, ivac, tol)
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(INOUT)  :: mpdata
      TYPE(t_vacuum), INTENT(IN)     :: vacuum
      INTEGER, INTENT(IN)            :: ivac
      REAL, INTENT(IN)               :: tol

      INTEGER :: ilen, i, nn, nz, ilo, ihi
      REAL    :: ovl, nrm, g
      REAL, ALLOCATABLE :: f(:), seed(:)

      nz = mpdata%nmz_vac
      ALLOCATE (f(nz), seed(nz))

      ilo = MINLOC(ABS(mpdata%glen_vac), DIM=1); ihi = ilo
      IF (ABS(mpdata%glen_vac(ilo)) > 1e-10) THEN
         DEALLOCATE (f, seed); RETURN
      END IF

      DO ilen = ilo, ihi
         g = mpdata%glen_vac(ilen)
         nn = mpdata%num_zbasfn_vac(ilen, ivac)
         IF (nn == 0) CYCLE
         IF (nn + 1 > SIZE(mpdata%zbasfn_vac, 2)) CYCLE

         DO i = 1, nz
            seed(i) = EXP(-g*(i - 1)*vacuum%delz)
         END DO
         f = seed
         DO i = 1, nn
            ovl = z_dot(mpdata%zbasfn_vac(:, i, ilen, ivac), seed, nz, vacuum%delz)
            f = f - ovl*mpdata%zbasfn_vac(:nz, i, ilen, ivac)
         END DO
         nrm = z_dot(f, f, nz, vacuum%delz)
         IF (nrm > tol) THEN
            mpdata%zbasfn_vac(:nz, nn + 1, ilen, ivac) = f/SQRT(nrm)
            mpdata%num_zbasfn_vac(ilen, ivac) = nn + 1
         END IF
      END DO
      DEALLOCATE (f, seed)
   END SUBROUTINE add_vac_moment_fun

   !>Rotate the first z-moment of the charge-free functions into slot nn-1 (g = 0 only).
   SUBROUTINE concentrate_first_moment(mpdata, vacuum, ilen, ivac)
      USE m_coulomb_vac, ONLY: vac_mom_g0
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(INOUT) :: mpdata
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      INTEGER, INTENT(IN)           :: ilen, ivac

      INTEGER :: nz, nn, i, ilast, imax
      REAL    :: mu_i, mu_n, nrm, mmax
      REAL, ALLOCATABLE :: hlp(:), swp(:)

      nz = mpdata%nmz_vac
      nn = mpdata%num_zbasfn_vac(ilen, ivac)
      IF (nn <= 2) RETURN

      ALLOCATE (hlp(nz), swp(nz))
      ilast = nn - 1

      mmax = -1.0; imax = ilast
      DO i = 1, ilast
         mu_i = ABS(zmom(mpdata%zbasfn_vac(:, i, ilen, ivac)))
         IF (mu_i > mmax) THEN
            mmax = mu_i; imax = i
         END IF
      END DO
      IF (imax /= ilast) THEN
         swp = mpdata%zbasfn_vac(:nz, ilast, ilen, ivac)
         mpdata%zbasfn_vac(:nz, ilast, ilen, ivac) = mpdata%zbasfn_vac(:nz, imax, ilen, ivac)
         mpdata%zbasfn_vac(:nz, imax, ilen, ivac) = swp
      END IF

      DO i = 1, ilast - 1
         mu_i = zmom(mpdata%zbasfn_vac(:, i, ilen, ivac))
         mu_n = zmom(mpdata%zbasfn_vac(:, ilast, ilen, ivac))
         nrm = SQRT(mu_n**2 + mu_i**2)
         IF (nrm < 1e-12) CYCLE
         hlp = mpdata%zbasfn_vac(:nz, ilast, ilen, ivac)
         mpdata%zbasfn_vac(:nz, ilast, ilen, ivac) = &
            mu_n/nrm*hlp + mu_i/nrm*mpdata%zbasfn_vac(:nz, i, ilen, ivac)
         mpdata%zbasfn_vac(:nz, i, ilen, ivac) = &
            mu_i/nrm*hlp - mu_n/nrm*mpdata%zbasfn_vac(:nz, i, ilen, ivac)

         mu_i = zmom(mpdata%zbasfn_vac(:, i, ilen, ivac))
         IF (ABS(mu_i) > 1e-10*MAX(1.0, mmax)) &
            CALL juDFT_error("first z-moment of a vacuum function does not vanish", &
                             calledby="concentrate_first_moment")
      END DO
      DEALLOCATE (hlp, swp)

   CONTAINS
      REAL FUNCTION zmom(f)
         IMPLICIT NONE
         REAL, INTENT(IN) :: f(:)
         REAL :: q0, q1
         CALL vac_mom_g0(f, nz, vacuum%delz, q0, q1)
         zmom = q1
      END FUNCTION zmom
   END SUBROUTINE concentrate_first_moment

   !>Givens rotations moving the exponential moment into the last function of the group.
   SUBROUTINE concentrate_exp_moment(mpdata, vacuum, ilen, ivac, g, fmpi)
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(INOUT) :: mpdata
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      INTEGER, INTENT(IN)           :: ilen, ivac
      REAL, INTENT(IN)              :: g
      TYPE(t_mpi), INTENT(IN)       :: fmpi

      INTEGER :: nz, nn, i, ilast, imax, iz
      REAL    :: mu_i, mu_n, nrm, mmax
      REAL, ALLOCATABLE :: w(:), hlp(:), swp(:)

      nz = mpdata%nmz_vac
      nn = mpdata%num_zbasfn_vac(ilen, ivac)
      IF (nn <= 1) RETURN

      ALLOCATE (w(nz), hlp(nz), swp(nz))
      DO iz = 1, nz
         w(iz) = EXP(-g*(iz - 1)*vacuum%delz)
      END DO

      mmax = -1.0; imax = nn
      DO i = 1, nn
         mu_i = ABS(mom(mpdata%zbasfn_vac(:, i, ilen, ivac)))
         IF (mu_i > mmax) THEN
            mmax = mu_i; imax = i
         END IF
      END DO
      IF (imax /= nn) THEN
         swp = mpdata%zbasfn_vac(:nz, nn, ilen, ivac)
         mpdata%zbasfn_vac(:nz, nn, ilen, ivac) = mpdata%zbasfn_vac(:nz, imax, ilen, ivac)
         mpdata%zbasfn_vac(:nz, imax, ilen, ivac) = swp
      END IF

      ilast = nn
      DO i = 1, nn - 1
         mu_i = mom(mpdata%zbasfn_vac(:, i, ilen, ivac))
         mu_n = mom(mpdata%zbasfn_vac(:, ilast, ilen, ivac))
         nrm = SQRT(mu_n**2 + mu_i**2)
         IF (nrm < 1e-12) CYCLE
         hlp = mpdata%zbasfn_vac(:nz, ilast, ilen, ivac)
         mpdata%zbasfn_vac(:nz, ilast, ilen, ivac) = &
            mu_n/nrm*hlp + mu_i/nrm*mpdata%zbasfn_vac(:nz, i, ilen, ivac)
         mpdata%zbasfn_vac(:nz, i, ilen, ivac) = &
            mu_i/nrm*hlp - mu_n/nrm*mpdata%zbasfn_vac(:nz, i, ilen, ivac)

         mu_i = mom(mpdata%zbasfn_vac(:, i, ilen, ivac))
         IF (ABS(mu_i) > 1e-10*MAX(1.0, mmax)) THEN
            CALL juDFT_error("exponential moment of a vacuum function does not vanish", &
                             calledby="concentrate_exp_moment")
         END IF
      END DO

      DEALLOCATE (w, hlp, swp)

   CONTAINS
      REAL FUNCTION mom(f)
         USE m_coulomb_vac, ONLY: vac_exp_mom
         IMPLICIT NONE
         REAL, INTENT(IN) :: f(:)
         mom = vac_exp_mom(f, nz, vacuum%delz, g)
      END FUNCTION mom
   END SUBROUTINE concentrate_exp_moment

   SUBROUTINE check_exp_moment(mpdata, vacuum, nmore)
      USE m_coulomb_vac, ONLY: vac_exp_mom
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(IN) :: mpdata
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      INTEGER, INTENT(OUT)       :: nmore

      INTEGER :: ivac, ilen, i, iz, nz, nn, ncarry
      REAL    :: g, r, mmax
      REAL, ALLOCATABLE :: w(:), pr(:)

      nz = mpdata%nmz_vac
      ALLOCATE (w(nz), pr(nz))
      nmore = 0
      DO ivac = 1, NVAC_MPB
         DO ilen = 1, SIZE(mpdata%glen_vac)
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            g = mpdata%glen_vac(ilen)
            DO iz = 1, nz
               w(iz) = EXP(-g*(iz - 1)*vacuum%delz)
            END DO
            mmax = 0.0
            DO i = 1, nn
               r = vac_exp_mom(mpdata%zbasfn_vac(:, i, ilen, ivac), nz, vacuum%delz, g)
               mmax = MAX(mmax, ABS(r))
            END DO
            ncarry = 0
            DO i = 1, nn
               r = vac_exp_mom(mpdata%zbasfn_vac(:, i, ilen, ivac), nz, vacuum%delz, g)
               IF (ABS(r) > 1e-8*MAX(mmax, 1e-30)) ncarry = ncarry + 1
            END DO
            IF (ncarry > 1) nmore = nmore + 1
         END DO
      END DO
   END SUBROUTINE check_exp_moment
   SUBROUTINE diag_sym(a, eig, eigv)
      IMPLICIT NONE
      REAL, INTENT(IN)               :: a(:, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: eig(:), eigv(:, :)

      INTEGER :: n, info, lwork
      REAL :: query(1)
      REAL, ALLOCATABLE :: work(:)

      n = SIZE(a, 1)
      ALLOCATE (eig(n), eigv(n, n))
      eigv = a
      lwork = -1
      CALL dsyev('V', 'U', n, eigv, n, eig, query, lwork, info)
      lwork = MAX(3*n, NINT(query(1)))
      ALLOCATE (work(lwork))
      CALL dsyev('V', 'U', n, eigv, n, eig, work, lwork, info)
      IF (info /= 0) CALL juDFT_error("dsyev failed in the vacuum mixed basis", &
                                      calledby="diag_sym")
      DEALLOCATE (work)
   END SUBROUTINE diag_sym

END MODULE m_mixedbasis_vac
