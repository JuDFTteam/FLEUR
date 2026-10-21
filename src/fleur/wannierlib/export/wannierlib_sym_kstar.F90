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
!> The operation index is the one FLEUR uses, so a consumer that is also handed sym%mrot and
!> sym%tau can line the two up. The reciprocal-space rotation is built exactly as gen_bz
!> builds it, TRANSPOSE(mrot(invtab)), rather than from mrot directly: a k vector transforms
!> under the inverse transpose, and getting that wrong permutes the star without raising
!> anything.
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
      INTEGER :: nsym = 0                   !< unitary operations used, = sym%nop
      INTEGER :: nk = 0                     !< points of the full mesh, = kpts%nkptf
      INTEGER, ALLOCATABLE :: kirr(:)       !< (nkirr) the representative of each star
      INTEGER, ALLOCATABLE :: k2kirr(:)     !< (nk) the star this k belongs to
      INTEGER, ALLOCATABLE :: k2kirr_sym(:) !< (nk) an operation taking the representative here
      INTEGER, ALLOCATABLE :: kirr2k(:, :)  !< (nkirr,nsym) where each operation sends it
      INTEGER, ALLOCATABLE :: weight(:)     !< (nkirr) points in the star
      !> (nsym) .TRUE. where the operation carries time reversal, i.e. acts as
      !> k -> -R k and conjugates. The first nop are always unitary and in FLEUR's own
      !> order, so an index below nop means the same operation with or without l_tr.
      LOGICAL, ALLOCATABLE :: antiunit(:)
      INTEGER :: nop = 0                    !< unitary operations, = sym%nop
   END TYPE t_wl_kstar

CONTAINS

   !> Every rank builds the tables: they are integers, the loop is cheap next to anything
   !> else in the run, and d_band will be wanted per k on whichever rank owns it. Only the
   !> report is rank 0, because oUnit is open there.
   SUBROUTINE wl_build_kstar(kpts, sym, fmpi, star, l_tr)
      USE m_types_kpts; USE m_types_sym; USE m_types_mpi
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_wl_kstar), INTENT(OUT) :: star
      LOGICAL, INTENT(IN) :: l_tr   !> add the antiunitary half

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
      star%nop = sym%nop
      star%nsym = MERGE(2*sym%nop, sym%nop, l_tr)
      ALLOCATE (rrot(3, 3, star%nsym), star%antiunit(star%nsym))
      star%antiunit = .FALSE.
      DO iop = 1, sym%nop
         rrot(:, :, iop) = TRANSPOSE(sym%mrot(:, :, sym%invtab(iop)))
      END DO
      IF (l_tr) THEN
         DO iop = 1, sym%nop
            rrot(:, :, sym%nop + iop) = -rrot(:, :, iop)
            star%antiunit(sym%nop + iop) = .TRUE.
         END DO
      END IF

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
      IF (fmpi%irank == 0) CALL wl_report_kstar(star, kpts%nkpt3)
   END SUBROUTINE wl_build_kstar

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

   SUBROUTINE wl_report_kstar(star, nkpt3)
      TYPE(t_wl_kstar), INTENT(IN) :: star
      INTEGER, INTENT(IN) :: nkpt3(3)
      INTEGER :: ir, w, wmax

      WRITE (oUnit, '(/,a)') "Wannier90/WannierBerri symmetry export: k-point stars"
      WRITE (oUnit, '(a,2(i0,a),i0,a,i0)') "  mesh ", nkpt3(1), "x", nkpt3(2), "x", &
         nkpt3(3), ", points ", star%nk
      WRITE (oUnit, '(a,i0)') "  unitary operations used  ", star%nop
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
