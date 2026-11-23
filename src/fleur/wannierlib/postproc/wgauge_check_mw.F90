!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Wannier90's overlap matrix against ours, measured rather than assumed.
!>
!> Both sides hold the same quantity by two different routes. Wannier90 starts from the
!> Bloch overlaps it was handed and rotates them in place: the disentanglement writes
!> u_opt^dagger M u_opt into the num_wann corner (disentangle.F90:383), the projection then
!> applies u_matrix (overlap.F90:999), and the minimisation rotates M in step with every
!> update of u_matrix (wannierise.F90:1761). What it leaves behind is therefore
!> M^(W,b)(k) = V^dagger(k) M^(0,b)(k) V(k_b) with V = u_opt.u_matrix -- which is the same
!> line wgauge_coeff_a builds from the untouched Bloch overlaps and the final gauge.
!>
!> WHY THIS IS WORTH A CHECK. It is the one place where an independent implementation of
!> the Wannier gauge is available to compare against, and it covers the whole chain at
!> once: the b-shell neighbour list, the order of the two rotations, the conjugation of
!> the overlaps, and which of u_opt and u_matrix goes on which side. A sign or a transpose
!> anywhere in that chain gives a result that still looks like an overlap matrix.
!>
!> AND WHY IT IS WORTH HAVING THE NUMBER. The Wannier90 side is what leaves this run: it
!> is the nw x nw block of the checkpoint handed to whoever carries on. That reader has
!> the block and no Bloch overlaps at all, so it cannot rebuild the quantity and has to
!> take the one in the file. This is the only run in which both sides exist, so it is the
!> only place the two can be shown to be interchangeable.
MODULE m_wgauge_check_mw
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_mpi, ONLY: t_mpi
   USE m_types_wgauge_bmesh
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wgauge_check_mw
CONTAINS

   !> tol is relative to the largest element of the overlap. The two sides are the same
   !> matrix products in a different order -- zgemm against MATMUL, and the rotations
   !> accumulated one minimisation step at a time against applied once at the end -- so
   !> they agree to rounding and not exactly, and demanding exactness here would be
   !> demanding that two implementations round the same way.
   !>
   !> MEASURED, on Pt with 108 bands disentangled to 54 Wannier functions over 100 k:
   !> 2.1E-15 after one minimisation step, 1.1E-14 after fifty and 2.9E-14 after three
   !> thousand. The drift is accumulation, so it tracks the step count rather than the
   !> system, and 1.0E-6 leaves it eight decades of room.
   SUBROUTINE wgauge_check_mw(mmn_loc, distk, u_opt, u_matrix, bmesh, mw_chk, fmpi, tol)
#ifdef CPP_MPI
      USE mpi
#endif
      COMPLEX,                 INTENT(IN) :: mmn_loc(:, :, :, :)   !< (nb,nb,nntot,nk_loc) Bloch
      INTEGER,                 INTENT(IN) :: distk(:)              !< (nk) owning rank of each k
      COMPLEX,                 INTENT(IN) :: u_opt(:, :, :)        !< (nb,nw,nk)
      COMPLEX,                 INTENT(IN) :: u_matrix(:, :, :)     !< (nw,nw,nk)
      TYPE(t_wgauge_bmesh),    INTENT(IN) :: bmesh
      COMPLEX,                 INTENT(IN) :: mw_chk(:, :, :, :)    !< (nw,nw,nntot,nk_loc) w90's
      TYPE(t_mpi),             INTENT(IN) :: fmpi
      REAL,                    INTENT(IN) :: tol

      COMPLEX, ALLOCATABLE :: Vk(:, :), Vkb(:, :), Mw(:, :)
      REAL    :: ds(2), dloc, dmine
      !> Where the worst element sits, carried alongside its size. A number on its own
      !> says the two sides disagree; the band pair, the neighbour and the k say which
      !> of the things that could be wrong actually is -- a disagreement confined to the
      !> off-diagonal is not the same defect as one that reaches the diagonal, and one at
      !> a single neighbour is not the same as one spread over all of them.
      INTEGER :: worst(4)
      INTEGER :: nb, nw, nnt, nk, nkl, kl, k, kb, nn, ierr, rk, ij(2)

      !> Nothing to compare is said out loud and not skipped quietly: this is the only
      !> cross-check of the gauge chain there is, and a run that lost it should say so.
      IF (SIZE(mw_chk, 1) == 0 .OR. SIZE(mmn_loc, 1) == 0) THEN
         IF (fmpi%irank == 0) WRITE (oUnit, '(a)') 'wannierlib: the Wannier90 overlaps and '// &
            'ours were not compared, one of the two is not there'
         RETURN
      END IF

      nb = SIZE(mmn_loc, 1); nw = SIZE(mw_chk, 1)
      nnt = SIZE(mmn_loc, 3); nk = SIZE(u_matrix, 3)
      nkl = COUNT(distk == fmpi%irank)

      IF (SIZE(u_opt, 1) /= nb .OR. SIZE(u_opt, 2) /= nw) &
         CALL juDFT_error("wannierlib: the gauge does not match the overlaps it is to "// &
                          "rotate", calledby="wgauge_check_mw")
      IF (SIZE(mw_chk, 3) /= nnt) &
         CALL juDFT_error("wannierlib: the two overlaps have different neighbour counts", &
                          calledby="wgauge_check_mw")

      ALLOCATE (Vk(nb, nw), Vkb(nb, nw), Mw(nw, nw))
      ds = 0.0
      worst = 0
      dloc = 0.0
      !> Same traversal as the operator: the slice is this rank's k in the order distk
      !> gives them, and the gauge is indexed by the GLOBAL k. Mixing the two up is the
      !> mistake this check exists to catch, so it is spelled out rather than reused.
      kl = 0
      DO k = 1, nk
         IF (distk(k) /= fmpi%irank) CYCLE
         kl = kl + 1
         IF (kl > nkl .OR. kl > SIZE(mmn_loc, 4)) EXIT
         Vk = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
         DO nn = 1, nnt
            kb = bmesh%nnlist(k, nn)
            IF (kb < 1 .OR. kb > nk) CALL juDFT_error("wannierlib: neighbour index out of "// &
                                                      "range", calledby="wgauge_check_mw")
            Vkb = MATMUL(u_opt(:, :, kb), u_matrix(:, :, kb))
            Mw = MATMUL(CONJG(TRANSPOSE(Vk)), MATMUL(mmn_loc(:, :, nn, kl), Vkb))
            dloc = MAXVAL(ABS(Mw - mw_chk(:, :, nn, kl)))
            IF (dloc > ds(1)) THEN
               ij = MAXLOC(ABS(Mw - mw_chk(:, :, nn, kl)))
               worst = [ij(1), ij(2), nn, k]
            END IF
            ds(1) = MAX(ds(1), dloc)
            ds(2) = MAX(ds(2), MAXVAL(ABS(Mw)))
         END DO
      END DO

      !> Reduced before it is judged and before it is printed: each rank sees only its own
      !> k, so a verdict taken on one of them is a verdict about a sixteenth of the mesh.
      dmine = ds(1)
#ifdef CPP_MPI
      IF (fmpi%isize > 1) CALL MPI_ALLREDUCE(MPI_IN_PLACE, ds, 2, MPI_DOUBLE_PRECISION, &
                                             MPI_MAX, fmpi%mpi_comm, ierr)
      !> Whose worst element it was, decided before it is reported. Taking the lowest rank
      !> among those holding the maximum keeps the line the same from one run to the next
      !> even when two ranks tie, which they do whenever the difference is exactly zero.
      rk = fmpi%isize
      IF (dmine == ds(1)) rk = fmpi%irank
      IF (fmpi%isize > 1) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, rk, 1, MPI_INTEGER, MPI_MIN, fmpi%mpi_comm, ierr)
         CALL MPI_BCAST(worst, 4, MPI_INTEGER, rk, fmpi%mpi_comm, ierr)
      END IF
#endif
      IF (fmpi%irank == 0) WRITE (oUnit, '(a,es10.2,a,es10.2,a,es8.1,a,3(i0,a),i0,a)') &
         'wannierlib: Wannier90 overlaps against ours, |w90 - wgauge| ', ds(1), &
         '   over a largest element of ', ds(2), '   allowed, relative ', tol, &
         '   worst at (', worst(1), ',', worst(2), ') neighbour ', worst(3), ' k ', &
         worst(4), MERGE('  ON the diagonal ', '  off the diagonal', worst(1) == worst(2))

      IF (ds(2) == 0.0) &
         CALL juDFT_error("wannierlib: the overlap matrix in the Wannier gauge is "// &
                          "identically zero", hint='nothing was compared, so the gauge '// &
                          'chain is unchecked', calledby="wgauge_check_mw")
      IF (ds(1) > tol*ds(2)) &
         CALL juDFT_error("wannierlib: Wannier90 and wgauge disagree on the overlap "// &
                          "matrix in the Wannier gauge", hint='a rotation, a conjugation '// &
                          'or a neighbour index is wrong on one of the two sides', &
                          calledby="wgauge_check_mw")
   END SUBROUTINE wgauge_check_mw

END MODULE m_wgauge_check_mw
