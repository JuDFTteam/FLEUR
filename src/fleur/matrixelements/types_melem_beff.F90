!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Muffin-tin matrix elements of the exchange field B_xc = (V_up - V_dn)/2, as a
!>  table per atom type that does not depend on k:
!>
!>      g(lm_a, i, lm_b, j, a, b, n) = < r_i^a Y_lm_a |  B_xc  | r_j^b Y_lm_b >
!>
!>  with a, b the spin channel the radial function is taken from. It is the weight
!>  the exchange-torque operator contracts the abc coefficients with, in the place
!>  where the spin operator uses radfun%integral.
!>
!>  The two differ in one way that decides the shape of this type: radfun%integral
!>  weighs a plain overlap, which the orthogonality of the Y_lm makes diagonal in
!>  lm, whereas B_xc is NOT spherical. It carries the lattice-harmonic expansion
!>  of the potential, so it couples lm_a to lm_b through Gaunt coefficients and the
!>  table needs both indices.
!>
!>  (lm, radial slot) are adjacent on each side on purpose: that is what lets the
!>  consumer contract the table in two passes over a composite index instead of
!>  running the nested lm sums band by band, which costs the square of what the
!>  spin operator does.
!>
!>  One table covers every radial function: the unified slot basis of t_radfun has
!>  u, udot and the local orbitals of an l as slots 1, 2, 3... of the same array,
!>  so the nine separate tables of the original (uju, ujd, dju, djd, ujulo, djulo,
!>  loju, lojd, ulojulo) collapse into the (i,j) pair of indices here.
MODULE m_types_melem_beff
   USE m_judft
   USE m_constants, ONLY: sfp_const, ImagUnit
   USE m_types_atoms
   USE m_types_sym
   USE m_types_sphhar
   USE m_types_potden
   USE m_types_radfun
   USE m_intgr, ONLY: intgr3
   USE m_gaunt, ONLY: gaunt1
   IMPLICIT NONE
   PRIVATE

   TYPE :: t_melem_beff
      INTEGER :: lmd  = -1                     !> largest lm index, (lmaxd+1)^2 - 1
      INTEGER :: n_rd = 0                      !> largest number of radial slots of any l
      COMPLEX, ALLOCATABLE :: g(:, :, :, :, :, :, :)
      !> (0:lmd, n_rd, 0:lmd, n_rd, 2, 2, ntype) = (lm_a, slot_a, lm_b, slot_b, spin_a, spin_b, type)
   CONTAINS
      PROCEDURE :: init => melem_beff_init
   END TYPE t_melem_beff

   PUBLIC :: t_melem_beff

CONTAINS

   !> Build the table for every atom type. radfun must already hold the radial
   !> functions of both spin channels, which is what makes the spin pairing below
   !> meaningful.
   SUBROUTINE melem_beff_init(this, atoms, sym, sphhar, vtot, radfun)
      CLASS(t_melem_beff), INTENT(INOUT) :: this
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_sym),         INTENT(IN)    :: sym
      TYPE(t_sphhar),      INTENT(IN)    :: sphhar
      TYPE(t_potden),      INTENT(IN)    :: vtot
      TYPE(t_radfun),      INTENT(IN)    :: radfun(:)   !> (ntype)

      REAL,    ALLOCATABLE :: beff(:, :, :)      ! (jmtd, 0:nlhd, ntype)
      REAL,    ALLOCATABLE :: rad(:, :, :, :, :) ! (n_rd, n_rd, 0:lmaxd, 0:lmaxd, 0:nlhd)
      REAL,    ALLOCATABLE :: x(:)
      INTEGER :: n, ns, lh, lpp, la, lb, ma, mb, mpp, lma, lmb
      INTEGER :: jmem, jj, ia, ib, lmini, lmaxi, jri, sa, sb
      COMPLEX :: factor
      LOGICAL :: l_doit

      IF (SIZE(radfun) /= atoms%ntype) &
         CALL judft_bug("melem_beff_init: radfun must carry one entry per atom type")
      IF (SIZE(vtot%mt, 4) < 2) &
         CALL judft_error('the exchange torque needs a spin-polarised potential', &
                          hint='B_xc = (V_up - V_dn)/2 is identically zero for jspins=1, '// &
                               'so the operator would return nothing but noise', &
                          calledby='melem_beff_init')

      this%lmd  = (atoms%lmaxd + 1)**2 - 1
      this%n_rd = 0
      DO n = 1, atoms%ntype
         this%n_rd = MAX(this%n_rd, MAXVAL(radfun(n)%n_r(0:atoms%lmax(n))))
      END DO

      IF (ALLOCATED(this%g)) DEALLOCATE(this%g)
      ALLOCATE(this%g(0:this%lmd, this%n_rd, 0:this%lmd, this%n_rd, 2, 2, atoms%ntype))
      this%g = CMPLX(0.0, 0.0)

      ! ---- the exchange field in the lattice-harmonic expansion ----
      !> vgen_finalize stores the spherical channel as r*V/sqrt(4 pi); the factor is
      !> taken out here so that every lh channel is a plain V(r) again.
      ALLOCATE(beff(SIZE(vtot%mt, 1), 0:sphhar%nlhd, atoms%ntype))
      beff = (vtot%mt(:, :, :, 1) - vtot%mt(:, :, :, 2)) / 2.0
      DO n = 1, atoms%ntype
         jri = atoms%jri(n)
         beff(:jri, 0, n) = beff(:jri, 0, n) * sfp_const / atoms%rmsh(:jri, n)
      END DO

      ALLOCATE(rad(this%n_rd, this%n_rd, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:sphhar%nlhd))
      ALLOCATE(x(SIZE(vtot%mt, 1)))

      DO n = 1, atoms%ntype
         !> the site symmetry of this type, named the way the rest of the tree names it
         ns  = sym%ntypsy(atoms%firstAtom(n))
         jri = atoms%jri(n)

         DO sb = 1, 2
            DO sa = 1, 2
               rad = 0.0

               ! ---- radial integrals, one per lattice harmonic and pair of l ----
               DO lh = 0, sphhar%nlh(ns)
                  lpp = sphhar%llh(lh, ns)
                  DO lb = 0, atoms%lmax(n)
                     DO la = 0, atoms%lmax(n)
                        !> The Gaunt coefficient vanishes outside the triangle and for odd
                        !> la + lb + lpp, so the integral is not worth doing there.
                        lmini = ABS(la - lb); lmaxi = la + lb
                        IF (MOD(la + lb + lpp, 2) == 1 .OR. lpp < lmini .OR. lpp > lmaxi) CYCLE
                        DO ib = 1, radfun(n)%n_r(lb)
                           DO ia = 1, radfun(n)%n_r(la)
                              !> both components of the scalar-relativistic pair are summed
                              x(:jri) = ( radfun(n)%r(:jri, 1, ia, la, sa)    &
                                        * radfun(n)%r(:jri, 1, ib, lb, sb)    &
                                        + radfun(n)%r(:jri, 2, ia, la, sa)    &
                                        * radfun(n)%r(:jri, 2, ib, lb, sb) )  &
                                        * beff(:jri, lh, n)
                              CALL intgr3(x, atoms%rmsh(1:, n), atoms%dx(n), jri, &
                                          rad(ia, ib, la, lb, lh))
                           END DO
                        END DO
                     END DO
                  END DO
               END DO

               ! ---- dress them with the Gaunt coefficients ----
               DO lb = 0, atoms%lmax(n)
                  DO mb = -lb, lb
                     lmb = lb*(lb + 1) + mb
                     DO la = 0, atoms%lmax(n)
                        DO ma = -la, la
                           lma = la*(la + 1) + ma
                           DO lh = 0, sphhar%nlh(ns)
                              lpp = sphhar%llh(lh, ns)
                              mpp = ma - mb
                              !> which member of this lattice harmonic carries mpp, if any
                              jmem = 0
                              DO jj = 1, sphhar%nmem(lh, ns)
                                 IF (sphhar%mlh(jj, lh, ns) == mpp) THEN
                                    jmem = jj; EXIT
                                 END IF
                              END DO
                              lmini = ABS(lb - lpp); lmaxi = lb + lpp
                              l_doit =        lmini <= la
                              l_doit = l_doit .AND. la <= lmaxi
                              l_doit = l_doit .AND. MOD(la + lb + lpp, 2) == 0
                              l_doit = l_doit .AND. ABS(mpp) <= lpp
                              l_doit = l_doit .AND. jmem /= 0
                              IF (.NOT. l_doit) CYCLE
                              factor = ImagUnit**(lb - la) &
                                     * gaunt1(la, lpp, lb, ma, mpp, mb, atoms%lmaxd) &
                                     * sphhar%clnu(jmem, lh, ns)
                              DO ib = 1, radfun(n)%n_r(lb)
                                 DO ia = 1, radfun(n)%n_r(la)
                                    this%g(lma, ia, lmb, ib, sa, sb, n) =      &
                                       this%g(lma, ia, lmb, ib, sa, sb, n)     &
                                       + factor * rad(ia, ib, la, lb, lh)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO

            END DO
         END DO
      END DO

      DEALLOCATE(beff, rad, x)
   END SUBROUTINE melem_beff_init

END MODULE m_types_melem_beff
