!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Groups the k-points of the full-zone mesh into stars under the space group. This is what
!> an external symmetry consumer needs before any representation matrix can be written: both
!> Wannier90's .dmn and WannierBerri's symmetrizer are indexed by a star representative and
!> an operation.
!>
!> WHY THIS IS NOT TAKEN FROM kpts. t_kpts already carries bkp and bksym -- the irreducible
!> parent of each k, and the operation that reaches it. But a wannierisation runs on the
!> FULL-zone mesh, because Wannier90 needs every k, and gen_bz then builds the zone from a
!> list that is already complete: it reaches every point with the identity and leaves
!> bkp(ik) = ik and bksym(ik) = 1. The tables exist and they are trivial, so the stars are
!> rebuilt here from the mesh and sym%mrot.
!>
!> TIME REVERSAL IS OPTIONAL, AND WHY. A .dmn has no antiunitary slot: its kptsym assigns
!> every k to a rotation, so writing an antiunitary operation there would assign k-points to
!> operations that do not generate them, with nothing to complain. FLEUR's own file has room
!> for it, so the caller decides, and l_tr says which: under the unitary group alone every
!> image is an ordinary rotation and the irreducible set is larger; with time reversal the
!> operations double, exactly as gen_bz doubles them (rrot = -rrot for the second half,
!> since T sends k to -k), and the set shrinks.
!>
!> It is a symmetry of the HAMILTONIAN only without magnetism. T alone maps a magnet onto
!> the state with the moment reversed: measured on bcc Fe with spin-orbit, where T by itself
!> misses by 7.2e-03 while the eight operations that preserve M close at 2.9e-12. So the
!> caller gates it on the run being non-magnetic, and a magnetic run gets the unitary group,
!> which is right rather than merely smaller.
!>
!> NOT EVERY OPERATION OF sym CAN BE USED, AND THE TWO HALVES DO NOT HOLD THE SAME ONES.
!> mapatom leaves sym%mapped_atom at zero wherever it could not place an atom: soc_sym
!> rejects the operations that do not hold the spin quantisation axis, and the search runs
!> inside one type, so one that carries an atom into another type leaves a zero as well.
!> Those are operations the muffin-tin rotation has no destination for, so the unitary half
!> is what mapatom kept.
!>
!> The antiunitary half is a DIFFERENT set. Time reversal flips every moment, so on a
!> magnet with spin-orbit the operations that pair with it are the ones that REVERSE the
!> moments -- precisely the ones soc_sym struck out -- and they need an atom map mapatom
!> never stored. Measured on AFM Fe: those eight close at 2.7e-12, the same tier as the
!> eight unitary ones, because both come from rotating a k-point of the mesh.
!>
!> So a slot is read through op, which holds FLEUR's own index of its rotation, and never
!> by offsetting nop. The matrices written beside the stars are the ones op selects, so a
!> consumer lines its operations up with those and never with sym%mrot.
!>
!> The reciprocal-space rotation is built exactly as gen_bz builds it,
!> TRANSPOSE(mrot(invtab)), rather than from mrot directly: a k vector transforms under the
!> inverse transpose, and getting that wrong permutes the star without raising anything.
MODULE m_wannierlib_sym_kstar
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_wl_kstar, wl_build_kstar

   !> One entry per star, plus the map back from the mesh. Index conventions, since they are
   !> the thing a consumer gets wrong silently: kirr and kirr2k hold indices into the FULL
   !> mesh (1..nk, i.e. kpts%bkf), while k2kirr holds a star number (1..nkirr).
   TYPE t_wl_kstar
      INTEGER :: nkirr = 0                  !< stars found
      INTEGER :: nsym = 0                   !< operations used, antiunitary half included
      INTEGER :: nk = 0                     !< points of the full mesh, = kpts%nkptf
      INTEGER, ALLOCATABLE :: kirr(:)       !< (nkirr) the representative of each star
      INTEGER, ALLOCATABLE :: k2kirr(:)     !< (nk) the star this k belongs to
      INTEGER, ALLOCATABLE :: k2kirr_sym(:) !< (nk) an operation taking the representative here
      INTEGER, ALLOCATABLE :: kirr2k(:, :)  !< (nkirr,nsym) where each operation sends it
      INTEGER, ALLOCATABLE :: weight(:)     !< (nkirr) points in the star
      !> (nsym) .TRUE. where the operation carries time reversal, i.e. acts as
      !> k -> -R k and conjugates. The unitary slots come first, but the antiunitary ones
      !> are not the same rotations shifted by nop: read op, never iop - nop.
      LOGICAL, ALLOCATABLE :: antiunit(:)
      INTEGER :: nop = 0                    !< unitary operations used, <= sym%nop
      !> operations sym offers that the MOMENTS rule out: they neither hold nor
      !> reverse them, so they are no symmetry with or without time reversal.
      INTEGER :: ntilt = 0
      !> (nsym) FLEUR's own index of the rotation each slot carries, and (nsym) where the
      !> inverse of each slot sits in THIS list. The second is what the written file needs:
      !> an invtab taken from sym would point at an operation that was not written.
      INTEGER, ALLOCATABLE :: op(:)
      INTEGER, ALLOCATABLE :: invtab(:)
      !> (nat,nsym) where each operation sends each atom. Rebuilt rather than taken from
      !> sym%mapped_atom, which is zero for every operation soc_sym vetoed -- and the
      !> antiunitary half is made of exactly those.
      INTEGER, ALLOCATABLE :: amap(:, :)
   END TYPE t_wl_kstar

CONTAINS

   !> Every rank builds the tables: they are integers, the loop is cheap next to anything
   !> else in the run, and d_band will be wanted per k on whichever rank owns it. Only the
   !> report is rank 0, because oUnit is open there.
   SUBROUTINE wl_build_kstar(kpts, sym, atoms, cell, nococonv, fmpi, star, l_tr, l_mfilter)
      USE m_types_kpts; USE m_types_sym; USE m_types_mpi
      USE m_types_atoms; USE m_types_cell; USE m_types_nococonv
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_wl_kstar), INTENT(OUT) :: star
      !> Whether time reversal may be used at all, and whether the moments decide which
      !> operations it pairs with. Both are the caller's: l_mfilter asks for the magnetic
      !> rule, and without it every unitary operation gets an antiunitary twin, which is
      !> what a run with nothing to reverse wants.
      LOGICAL, INTENT(IN) :: l_tr, l_mfilter

      REAL, ALLOCATABLE :: rrot(:, :, :)
      INTEGER :: ik, jk, iop, ir, id(3, 3)

      !> The eigenvectors of a k that is not in the mesh cannot be rotated into place here,
      !> and the stars would be built over a wedge while the Wannier matrices live on the
      !> full zone: two different k indexings, which is the one error no later check sees.
      IF (kpts%nkpt /= kpts%nkptf) CALL juDFT_error( &
         "wannierlib: the symmetry export needs the full-zone k-mesh", &
         hint="generate it with inpgen -noKsym", calledby="wl_build_kstar")

      !> Both consumers require the first operation to be the identity, and so does the gate
      !> below that the representative maps to itself under it. gen_bz already aborts if it
      !> is not, but this module is reachable without having gone through gen_bz.
      id = 0; id(1, 1) = 1; id(2, 2) = 1; id(3, 3) = 1
      IF (ANY(sym%mrot(:, :, 1) /= id)) CALL juDFT_error( &
         "wannierlib: the first symmetry operation is not the identity", &
         calledby="wl_build_kstar")

      star%nk = kpts%nkptf
      CALL wl_build_ops(sym, atoms, cell, nococonv, star, l_tr, l_mfilter)
      ALLOCATE (rrot(3, 3, star%nsym))
      DO iop = 1, star%nsym
         rrot(:, :, iop) = TRANSPOSE(sym%mrot(:, :, sym%invtab(star%op(iop))))
         !> T sends k to -k, so the antiunitary slots carry the same rotation with the
         !> sign flipped, exactly as gen_bz doubles them.
         IF (star%antiunit(iop)) rrot(:, :, iop) = -rrot(:, :, iop)
      END DO

      ALLOCATE (star%k2kirr(star%nk), star%k2kirr_sym(star%nk))
      ALLOCATE (star%kirr(star%nk), star%kirr2k(star%nk, star%nsym), star%weight(star%nk))
      star%k2kirr = 0; star%k2kirr_sym = 0; star%weight = 0
      star%nkirr = 0

      !> Walk the mesh in order; the first point of a star that has not been claimed yet
      !> becomes its representative. Taking the lowest index keeps the set reproducible,
      !> which matters because a consumer stores it and compares runs.
      DO ik = 1, star%nk
         IF (star%k2kirr(ik) /= 0) CYCLE
         star%nkirr = star%nkirr + 1
         ir = star%nkirr
         star%kirr(ir) = ik
         DO iop = 1, star%nsym
            jk = kpts%get_nk(MATMUL(rrot(:, :, iop), kpts%bkf(:, ik)))
            IF (jk == 0) CALL juDFT_error( &
               "wannierlib: a rotated k-point is not in the mesh", &
               hint="the mesh is not invariant under the space group of inp.xml", &
               calledby="wl_build_kstar")
            star%kirr2k(ir, iop) = jk
            !> Several operations reach the same point when the star is shorter than the
            !> group; the first one is kept, and which one it is does not matter as long as
            !> the same table is handed out with it.
            IF (star%k2kirr(jk) == 0) THEN
               star%k2kirr(jk) = ir
               star%k2kirr_sym(jk) = iop
               star%weight(ir) = star%weight(ir) + 1
            END IF
         END DO
      END DO

      CALL wl_check_kstar(star)
      IF (fmpi%irank == 0) CALL wl_report_kstar(star, kpts%nkpt3, sym%nop)
   END SUBROUTINE wl_build_kstar

   !> The cartesian rotation of a lattice operation. x_cart = amat.x_frac is FLEUR's own
   !> idiom, so R = amat.mrot.amat^-1, and the inverse is read off bmat rather than built
   !> here, for the reason the body gives. A lattice operation is an isometry, so R comes
   !> out orthogonal, and that is checked rather than assumed: it is the only thing that
   !> says the convention is the right way round.
   SUBROUTINE wl_rcart(cell, mrot, r, det)
      USE m_types_cell
      USE m_constants, ONLY: tpi_const
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER, INTENT(IN) :: mrot(3, 3)
      REAL, INTENT(OUT) :: r(3, 3), det

      REAL :: p(3, 3)
      INTEGER :: i

      !> The inverse is not computed here: cell%bmat is 2.pi times it. MEASURED on wurtzite
      !> AlN, which is as far from a symmetric amat as a cell gets -- MATMUL(amat, bmat/2pi)
      !> is the identity to 2e-16 while the same product against its transpose is off by
      !> 1.4, so the two are not interchangeable and this is the one that holds.
      !>
      !> xas_cart_rotation_from_sym is the same expression and has carried it since June;
      !> over the twelve operations of that cell the two agree to 3e-15. They are not one
      !> routine because the gate below is what this caller needs and that one does not
      !> have it, and because a shared one belongs beside euler in m_grp_k rather than in
      !> either caller -- which is a change to code neither of them owns.
      r = MATMUL(cell%amat, MATMUL(REAL(mrot), cell%bmat/tpi_const))

      p = MATMUL(r, TRANSPOSE(r))
      DO i = 1, 3
         p(i, i) = p(i, i) - 1.0
      END DO
      IF (MAXVAL(ABS(p)) > 1e-6) CALL juDFT_error( &
         "wannierlib: the cartesian rotation is not orthogonal", &
         hint="amat does not hold the lattice vectors in its columns here", &
         calledby="wl_rcart")
      det = r(1, 1)*(r(2, 2)*r(3, 3) - r(2, 3)*r(3, 2)) &
            - r(1, 2)*(r(2, 1)*r(3, 3) - r(2, 3)*r(3, 1)) &
            + r(1, 3)*(r(2, 1)*r(3, 2) - r(2, 2)*r(3, 1))
   END SUBROUTINE wl_rcart

   !> The direction of each atom's moment in the global frame, from the non-collinear angles
   !> of its type: m = (sin b cos a, sin b sin a, cos b). Only the direction is wanted, since
   !> the question is whether an operation carries one atom's onto its image's.
   SUBROUTINE wl_moments(atoms, nococonv, m)
      USE m_types_atoms; USE m_types_nococonv
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      REAL, INTENT(OUT) :: m(:, :)

      INTEGER :: n, na, ia
      REAL :: al, be

      DO n = 1, atoms%ntype
         al = nococonv%alph(n)
         be = nococonv%beta(n)
         DO na = 1, atoms%neq(n)
            ia = atoms%firstAtom(n) + na - 1
            m(:, ia) = [SIN(be)*COS(al), SIN(be)*SIN(al), COS(be)]
         END DO
      END DO
   END SUBROUTINE wl_moments

   !> Where an operation sends each atom, by the test mapatom uses: R.taual(a) + tau =
   !> taual(a') + L with L a lattice vector and a' of the same type. Rebuilt here because
   !> mapatom stores a zero wherever soc_sym vetoed the operation, and the antiunitary half
   !> is made of exactly those; l_ok is false where no image was found.
   SUBROUTINE wl_atom_map(atoms, sym, iop, map, l_ok)
      USE m_types_atoms; USE m_types_sym
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      INTEGER, INTENT(IN) :: iop
      INTEGER, INTENT(OUT) :: map(:)
      LOGICAL, INTENT(OUT) :: l_ok

      INTEGER :: n, na, nb, ia, ja
      REAL :: x(3), d(3)

      map = 0
      DO n = 1, atoms%ntype
         DO na = 1, atoms%neq(n)
            ia = atoms%firstAtom(n) + na - 1
            x = MATMUL(REAL(sym%mrot(:, :, iop)), atoms%taual(:, ia)) + sym%tau(:, iop)
            DO nb = 1, atoms%neq(n)
               ja = atoms%firstAtom(n) + nb - 1
               d = x - atoms%taual(:, ja)
               IF (ANY(ABS(d - ANINT(d)) > 1e-6)) CYCLE
               map(ia) = ja
               EXIT
            END DO
         END DO
      END DO
      l_ok = ALL(map > 0)
   END SUBROUTINE wl_atom_map

   !> The two halves of the operation list. The unitary one is what mapatom decided and is
   !> not revisited: a zero in sym%mapped_atom is FLEUR saying the operation is not a
   !> symmetry of this state, and dropping those leaves a group.
   !>
   !> The antiunitary one is chosen here. Time reversal flips every moment, so g x T is a
   !> symmetry when g REVERSES them: det(R).R.m(a) = -m(g a) for every atom, with the axial
   !> sign because a moment is a pseudovector. Without the magnetic rule every unitary
   !> operation gets a twin, which is the non-magnetic case and what l_tr alone used to mean.
   !>
   !> It is not taken on trust. The inverse of every slot has to be in the list, the rebuilt
   !> map has to resolve every operation mapatom resolved, and the operations mapatom kept
   !> have to HOLD the moments -- a control the unitary half is asked for nowhere else.
   SUBROUTINE wl_build_ops(sym, atoms, cell, nococonv, star, l_tr, l_mfilter)
      USE m_types_sym; USE m_types_atoms; USE m_types_cell; USE m_types_nococonv
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_wl_kstar), INTENT(INOUT) :: star
      LOGICAL, INTENT(IN) :: l_tr, l_mfilter

      LOGICAL :: l_u(sym%nop), l_a(sym%nop), l_ok
      INTEGER :: amap(atoms%nat, sym%nop)
      REAL :: m(3, atoms%nat), r(3, 3), det
      INTEGER :: i, j, ia, ntilt

      ntilt = 0
      l_u = .TRUE.
      IF (ALLOCATED(sym%mapped_atom)) THEN
         DO i = 1, sym%nop
            l_u(i) = ALL(sym%mapped_atom(i, :) > 0)
         END DO
      END IF

      DO i = 1, sym%nop
         CALL wl_atom_map(atoms, sym, i, amap(:, i), l_ok)
         IF (l_u(i) .AND. .NOT. l_ok) CALL juDFT_error( &
            "wannierlib: the rebuilt atom map misses an operation mapatom resolved", &
            calledby="wl_build_ops")
      END DO

      l_a = .FALSE.
      IF (l_tr .AND. .NOT. l_mfilter) THEN
         l_a = l_u
      ELSE IF (l_tr) THEN
         CALL wl_moments(atoms, nococonv, m)
         DO i = 1, sym%nop
            IF (ANY(amap(:, i) == 0)) CYCLE
            CALL wl_rcart(cell, sym%mrot(:, :, i), r, det)
            l_a(i) = .TRUE.
            DO ia = 1, atoms%nat
               IF (ANY(ABS(det*MATMUL(r, m(:, ia)) + m(:, amap(ia, i))) > 1e-6)) &
                  l_a(i) = .FALSE.
            END DO
            IF (.NOT. l_u(i)) CYCLE
            !> An operation that REVERSES every moment is a symmetry with time reversal and
            !> not without it, so it belongs in the antiunitary half. mapped_atom filters on
            !> the spin quantisation axis alone, which does not see a moment turned in the
            !> plane, and keeps such an operation in the unitary half; the two halves then
            !> overlap and the operation would be counted twice. It is moved rather than
            !> treated as a contradiction. The operations that hold the moments are a
            !> subgroup and the ones that reverse them its coset, so the unitary half stays
            !> closed -- which the inversion gate below checks and does not take on trust.
            IF (l_a(i)) THEN
               l_u(i) = .FALSE.
               CYCLE
            END IF
            !> Neither held nor reversed: the operation TILTS the moments, so it is no
            !> symmetry of this state with or without time reversal and it goes into
            !> neither half. mapped_atom does not see this, because it filters on the spin
            !> quantisation axis alone and a moment turned inside the plane keeps its atom
            !> map, so sym offers operations the magnetic state does not have.
            !>
            !> MEASURED on hcp Co with a basal moment: of twelve operations two hold the
            !> moments, two reverse them and eight tilt them, the threefold ones. The count
            !> is reported rather than swallowed, because a group cut from twelve to four
            !> in silence is what later reads as a bug. What is left is still checked: the
            !> ones that hold the moments are a subgroup and the ones that reverse them its
            !> coset, and the inversion gate below fails if that is not what came out.
            DO ia = 1, atoms%nat
               IF (ANY(ABS(det*MATMUL(r, m(:, ia)) - m(:, amap(ia, i))) > 1e-6)) THEN
                  l_u(i) = .FALSE.
                  ntilt = ntilt + 1
                  EXIT
               END IF
            END DO
         END DO
      END IF
      star%ntilt = ntilt

      star%nop = COUNT(l_u)
      star%nsym = star%nop + COUNT(l_a)
      ALLOCATE (star%op(star%nsym), star%antiunit(star%nsym))
      ALLOCATE (star%amap(atoms%nat, star%nsym))
      j = 0
      DO i = 1, sym%nop
         IF (.NOT. l_u(i)) CYCLE
         j = j + 1
         star%op(j) = i
         star%antiunit(j) = .FALSE.
         star%amap(:, j) = amap(:, i)
      END DO
      DO i = 1, sym%nop
         IF (.NOT. l_a(i)) CYCLE
         j = j + 1
         star%op(j) = i
         star%antiunit(j) = .TRUE.
         star%amap(:, j) = amap(:, i)
      END DO
      !> Both consumers want the identity first, and the caller has already checked that
      !> sym itself starts with it.
      IF (star%op(1) /= 1 .OR. star%antiunit(1)) CALL juDFT_error( &
         "wannierlib: the identity is not the first usable operation", &
         calledby="wl_build_ops")

      ALLOCATE (star%invtab(star%nsym), source=0)
      DO i = 1, star%nsym
         DO j = 1, star%nsym
            IF (star%antiunit(j) .NEQV. star%antiunit(i)) CYCLE
            IF (star%op(j) == sym%invtab(star%op(i))) star%invtab(i) = j
         END DO
      END DO
      IF (ANY(star%invtab == 0)) CALL juDFT_error( &
         "wannierlib: the usable operations are not closed under inversion", &
         calledby="wl_build_ops")
   END SUBROUTINE wl_build_ops

   !> The gates. Each one fails on a different mistake and none of them needs any physics,
   !> which is the point of doing this before the representation matrices exist.
   SUBROUTINE wl_check_kstar(star)
      TYPE(t_wl_kstar), INTENT(IN) :: star
      INTEGER :: ir, iop

      IF (ANY(star%k2kirr(:star%nk) == 0)) CALL juDFT_error( &
         "wannierlib: a k-point of the mesh is in no star", calledby="wl_check_kstar")
      IF (SUM(star%weight(:star%nkirr)) /= star%nk) CALL juDFT_error( &
         "wannierlib: the stars do not add up to the mesh", calledby="wl_check_kstar")
      DO ir = 1, star%nkirr
         !> Under the identity the representative has to be itself, which catches an
         !> operation index out of step: a kirr2k filled off by one column, or a first
         !> operation that is not the identity. It does NOT catch a reciprocal rotation
         !> built from mrot instead of its inverse transpose, and nothing here does -- the
         !> identity is the identity under either convention, and measurement says the wrong
         !> one passes every gate below with only the star-size histogram differing. That
         !> convention is secured by sharing gen_bz_internal's own line, not by a gate.
         IF (star%kirr2k(ir, 1) /= star%kirr(ir)) CALL juDFT_error( &
            "wannierlib: the identity does not fix the star representative", &
            calledby="wl_check_kstar")
         !> |star| * |little group| = |group|, so the star size divides the group order. A
         !> count that does not divide it means points were merged or missed.
         IF (MOD(star%nsym, star%weight(ir)) /= 0) CALL juDFT_error( &
            "wannierlib: a star size does not divide the number of operations", &
            calledby="wl_check_kstar")
         !> Every image of a representative has to report that same star. This is the gate
         !> that sees a set of operations not closed as a group: one operation outside it
         !> reaches points an earlier star already claimed, while every k is still in
         !> exactly one star and the weights still add up to the mesh. Measured by injecting
         !> a single integer matrix from outside a 48-operation group: all 47 injections
         !> fire here, and on none of them does any other gate fire.
         DO iop = 1, star%nsym
            IF (star%k2kirr(star%kirr2k(ir, iop)) /= ir) CALL juDFT_error( &
               "wannierlib: a rotated k-point belongs to another star", &
               calledby="wl_check_kstar")
         END DO
      END DO
   END SUBROUTINE wl_check_kstar

   !> Report the stars to the out file: the mesh, how many operations were usable, how many
   !> of sym were dropped for having no atom map, whether time reversal joined them as a set
   !> of its own, and how far the reduction went. This is the only place a reader of the out
   !> file learns that operations were dropped at all, which is why the count that was NOT
   !> used is printed beside the one that was.
   SUBROUTINE wl_report_kstar(star, nkpt3, nop_all)
      TYPE(t_wl_kstar), INTENT(IN) :: star
      INTEGER, INTENT(IN) :: nkpt3(3), nop_all
      INTEGER :: ir, w, wmax

      WRITE (oUnit, '(/,a)') "Wannier90/WannierBerri symmetry export: k-point stars"
      WRITE (oUnit, '(a,2(i0,a),i0,a,i0)') "  mesh ", nkpt3(1), "x", nkpt3(2), "x", &
         nkpt3(3), ", points ", star%nk
      WRITE (oUnit, '(a,i0)') "  unitary operations used  ", star%nop
      IF (star%nop < nop_all) &
         WRITE (oUnit, '(a,i0,a)') "  of ", nop_all, &
         " in sym; the rest are not symmetries of this state"
      IF (star%ntilt > 0) &
         WRITE (oUnit, '(a,i0,a)') "  of those, by the moments ", star%ntilt, &
         "; the rotation turns them, so no operation of the group holds them"
      IF (star%nsym > star%nop) &
         WRITE (oUnit, '(a,i0,a)') "  with time reversal       ", star%nsym, " in all"

      WRITE (oUnit, '(a,i0)') "  stars (irreducible k)    ", star%nkirr
      WRITE (oUnit, '(a,f6.2)') "  reduction factor         ", REAL(star%nk)/star%nkirr
      wmax = MAXVAL(star%weight(:star%nkirr))
      WRITE (oUnit, '(a)') "  star sizes:"
      DO w = 1, wmax
         ir = COUNT(star%weight(:star%nkirr) == w)
         IF (ir > 0) WRITE (oUnit, '(a,i0,a,i0,a)') "    ", ir, " star(s) of ", w, " points"
      END DO
      ! The representatives themselves, so that the partition can be compared against an
      ! independently built one star by star, not just through the histogram above: the
      ! size histogram alone does not pin down which k is in which star.
      IF (star%nkirr <= 200) THEN
         WRITE (oUnit, '(a)') "  stars as <k in full mesh>:<size>"
         WRITE (oUnit, '(a)', ADVANCE='no') "   "
         DO ir = 1, star%nkirr
            WRITE (oUnit, '(2(a,i0))', ADVANCE='no') "  ", star%kirr(ir), ":", star%weight(ir)
            IF (MOD(ir, 8) == 0 .AND. ir < star%nkirr) WRITE (oUnit, '(/,a)', ADVANCE='no') "   "
         END DO
         WRITE (oUnit, '(a)') ""
      END IF
   END SUBROUTINE wl_report_kstar

END MODULE m_wannierlib_sym_kstar
