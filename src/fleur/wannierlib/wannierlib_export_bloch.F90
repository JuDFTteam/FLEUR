!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The files written from the BLOCH basis, before any gauge exists: the eigenvalues in
!> Wannier90's .eig format, and the coarse spin operator.
!>
!> They are here and not in postproc/wgauge_io because that file is the layout of the
!> real-space operators O(R), which only exist once the gauge does. Keeping the two apart
!> is what lets either be described in one sentence: before the gauge, and after it.
MODULE m_wannierlib_export_bloch
   USE m_juDFT
   USE m_constants, ONLY: hartree_to_ev_const
   USE m_types_mpi, ONLY: t_mpi
   USE m_types_kpts
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannierlib_write_eig, wannierlib_write_s0, &
             wannierlib_write_mmn

CONTAINS

   ! Write the eigenvalues in the Wannier90 file format: one line per (band, k), in eV.
   ! The library hands w90 the same numbers through the API and never puts them on disk, so
   ! the .amn and .mmn alone cannot be replayed -- a standalone run stops at "No .eig file
   ! found. Needed for disentanglement". FLEUR keeps eigenvalues in Hartree; the file is in
   ! eV, the unit the energy windows are also written in.
   SUBROUTINE wannierlib_write_eig(eig, filename)
      REAL, INTENT(IN) :: eig(:, :)            ! (num_bands, nkptf), Hartree
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER :: iunit, iband, ikpt
      OPEN (NEWUNIT=iunit, FILE=filename, FORM='formatted', STATUS='replace')
      DO ikpt = 1, SIZE(eig, 2)
         DO iband = 1, SIZE(eig, 1)
            WRITE (iunit, '(2i12,f18.13)') iband, ikpt, hartree_to_ev_const*eig(iband, ikpt)
         END DO
      END DO
      CLOSE (iunit)
   END SUBROUTINE wannierlib_write_eig

   ! Write the spin operator in the BLOCH basis, before any gauge is applied to it.
   ! It is the one ingredient of S(R) that no file carries: the .amn and .mmn fix the
   ! gauge and the .eig the energies, but the operator itself only ever exists in memory,
   ! so a disagreement between it and the basis the gauge was built for cannot be seen
   ! from outside. The classic route publishes the same quantity as WF1.mmn0 (diagonal
   ! spin blocks) and updown.mmn0 (cross block), which is what this can be compared with.
   ! ONE file, gathered, like every other export: the rank's own k-slice is completed by
   ! the same reduce the overlaps use and rank 0 writes. A file per rank would have been
   ! cheaper here, but then the number of output files would depend on how the run was
   ! parallelised, which is not something a reader should have to know.
   SUBROUTINE wannierlib_write_s0(fmpi, distk, nkptf, s0_loc, filename)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)             ! (nkptf) owning rank of each k
      INTEGER, INTENT(IN) :: nkptf
      COMPLEX, INTENT(IN) :: s0_loc(:, :, :, :)   ! (nb, nb, 3, nk_local), ascending global k
      CHARACTER(LEN=*), INTENT(IN) :: filename
      COMPLEX, ALLOCATABLE :: s0_full(:, :, :, :)
      INTEGER :: iunit, ikpt, ik_local, ib, jb, ic, ierr
#ifdef CPP_MPI
      INTEGER :: mpi_err
#endif

      ALLOCATE (s0_full(SIZE(s0_loc, 1), SIZE(s0_loc, 2), SIZE(s0_loc, 3), nkptf), &
                stat=ierr, source=CMPLX(0.0, 0.0))
      IF (ierr /= 0) CALL juDFT_error('wannierlib_write_s0: allocation failed', &
                                      calledby='wannierlib_write_s0')
      !> Indexed by GLOBAL k and zero everywhere else, so the sum below completes it.
      !> The local array is indexed by position in this rank's slice, which is why it
      !> cannot be reduced as it stands: the same slot means a different k on each rank.
      ik_local = 0
      DO ikpt = 1, nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE
         ik_local = ik_local + 1
         IF (ik_local > SIZE(s0_loc, 4)) EXIT
         s0_full(:, :, :, ikpt) = s0_loc(:, :, :, ik_local)
      END DO
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, s0_full, SIZE(s0_full), MPI_DOUBLE_COMPLEX, &
                         MPI_SUM, fmpi%mpi_comm, mpi_err)
#endif
      IF (fmpi%irank == 0) THEN
         OPEN (NEWUNIT=iunit, FILE=filename, FORM='formatted', STATUS='replace')
         WRITE (iunit, '(a)') '# spin operator in the Bloch basis, before the gauge'
         WRITE (iunit, '(a)') '# ikpt  iband  jband  comp(1=x,2=y,3=z)  Re  Im'
         DO ikpt = 1, nkptf
            DO ic = 1, SIZE(s0_full, 3)
               DO jb = 1, SIZE(s0_full, 2)
                  DO ib = 1, SIZE(s0_full, 1)
                     WRITE (iunit, '(4i7,2es24.15)') ikpt, ib, jb, ic, &
                        REAL(s0_full(ib, jb, ic, ikpt)), AIMAG(s0_full(ib, jb, ic, ikpt))
                  END DO
               END DO
            END DO
         END DO
         CLOSE (iunit)
      END IF
      DEALLOCATE (s0_full)
   END SUBROUTINE wannierlib_write_s0


   subroutine wannierlib_write_mmn(this, mmnk, kpts, nnkp, gkpb, jspin, fending)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      COMPLEX, INTENT(IN) :: mmnk(:, :, :, :)
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: nnkp(:, :)
      INTEGER, INTENT(IN) :: gkpb(:, :, :)
      INTEGER, INTENT(IN), OPTIONAL :: jspin
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fending

      INTEGER :: ikpt, ikpt_b, i, j, nntot, nbnd, jspin2
      CHARACTER(LEN=3) :: spin12(2)
      CHARACTER(LEN=12) :: file_ending

      spin12 = (/'WF1', 'WF2'/)
      jspin2 = 1
      IF (PRESENT(jspin)) jspin2 = jspin
      IF (PRESENT(fending)) THEN
         file_ending = TRIM(fending)
      ELSE
         file_ending = ''
      END IF

      nntot = SIZE(mmnk, 3)
      nbnd = SIZE(mmnk, 1)

      IF (SIZE(mmnk, 2) /= nbnd) THEN
         CALL juDFT_error('wannierlib_write_mmn: mmnk must be square in band indices', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (SIZE(mmnk, 4) /= kpts%nkptf) THEN
         CALL juDFT_error('wannierlib_write_mmn: k-point dimension mismatch', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (SIZE(nnkp, 1) /= kpts%nkptf .OR. SIZE(nnkp, 2) /= nntot) THEN
         CALL juDFT_error('wannierlib_write_mmn: nnkp dimension mismatch', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (SIZE(gkpb, 1) /= 3 .OR. SIZE(gkpb, 2) /= kpts%nkptf .OR. SIZE(gkpb, 3) /= nntot) THEN
         CALL juDFT_error('wannierlib_write_mmn: gkpb dimension mismatch', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (jspin2 < 1 .OR. jspin2 > 2) THEN
         CALL juDFT_error('wannierlib_write_mmn: jspin must be 1 or 2', &
                          calledby='wannierlib_write_mmn')
      END IF

      open (305, file=spin12(jspin2)//TRIM(file_ending)//'.mmn')
      write (305, *) 'Overlaps of the wavefunct. the k- and b-points'
      write (305, '(3i5)') nbnd, kpts%nkptf, nntot
      DO ikpt = 1, kpts%nkptf
         DO ikpt_b = 1, nntot
            write (305, '(2i5,3x,3i4)') ikpt, nnkp(ikpt, ikpt_b), gkpb(1:3, ikpt, ikpt_b)
            DO i = 1, nbnd
               DO j = 1, nbnd
                  write (305, '(2f24.18)') real(mmnk(j, i, ikpt_b, ikpt)), &
                     aimag(mmnk(j, i, ikpt_b, ikpt))
               END DO
            END DO
         END DO
      END DO
      close (305)
   END SUBROUTINE wannierlib_write_mmn
END MODULE m_wannierlib_export_bloch
