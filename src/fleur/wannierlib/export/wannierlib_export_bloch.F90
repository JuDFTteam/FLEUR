!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The files written from the BLOCH basis, before any gauge exists: the eigenvalues in
!> Wannier90's .eig format, and the coarse spin operator.
!>
!> Everything here is written before any gauge exists; the real-space operators O(R),
!> which need one, are laid out elsewhere.
MODULE m_wannierlib_export_bloch
   USE m_juDFT
   USE m_constants, ONLY: hartree_to_ev_const
   USE m_types_mpi, ONLY: t_mpi
   USE m_types_kpts
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
#ifdef CPP_MPI
   USE mpi
#endif
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannierlib_write_eig, wannierlib_write_mmn, wannierlib_write_amn

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
   ! Write the projections in the Wannier90 .amn format: one line per (band, Wannier
   ! function, k-point), real and imaginary part. The library hands w90 the same matrix
   ! through the API and never puts it on disk, so without this file a standalone
   ! wannier90.x cannot repeat the run -- the sibling of the .mmn and .eig above.
   !
   ! The matrix arrives complete. The k-distributed contributions are reduced before this
   ! is reached, so there is nothing to gather here and only rank 0 writes.
   SUBROUTINE wannierlib_write_amn(amn, filename, title)
      COMPLEX, INTENT(IN) :: amn(:, :, :)      ! (num_bands, num_wann, nkptf)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      CHARACTER(LEN=*), INTENT(IN) :: title
      INTEGER :: iunit, iband, iwf, ikpt
      OPEN (NEWUNIT=iunit, FILE=filename, FORM='formatted', STATUS='replace')
      WRITE (iunit, *) title
      WRITE (iunit, '(i5,i7,i5)') SIZE(amn, 1), SIZE(amn, 3), SIZE(amn, 2)
      DO ikpt = 1, SIZE(amn, 3)
         DO iwf = 1, SIZE(amn, 2)
            DO iband = 1, SIZE(amn, 1)
               WRITE (iunit, '(i5,i5,i7,3x,2f18.12)') iband, iwf, ikpt, &
                  REAL(amn(iband, iwf, ikpt)), AIMAG(amn(iband, iwf, ikpt))
            END DO
         END DO
      END DO
      CLOSE (iunit)
   END SUBROUTINE wannierlib_write_amn

END MODULE m_wannierlib_export_bloch
