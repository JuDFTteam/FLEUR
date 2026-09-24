!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Writes the gauge Wannier90 returned: the disentanglement factor and the MLWF rotation,
!> one file per wannierised channel.
!>
!> WHY IT IS NEEDED. Every operator this code exports -- H(R), S(R), L(R), the SOC matrix --
!> is the Bloch-basis operator sandwiched between these two matrices, so the gauge is what
!> decides whether the result can be interpolated at all. It is also the only ingredient of
!> that sandwich which leaves no trace: the overlaps are in the .mmn, the projections in the
!> .amn, the eigenvalues in out.xml, but U(k) lives and dies inside the run. Anyone studying
!> a gauge-dependent quantity from outside has had to rebuild it with their own minimiser,
!> which reaches a different minimum and therefore answers a different question.
!>
!> THE TWO FACTORS ARE WRITTEN SEPARATELY, not their product, because they answer different
!> questions: u_opt fixes WHICH subspace was selected, u_mlwf fixes HOW it was rotated inside
!> that subspace. Effects that look alike from the outside are told apart only by attributing
!> them to one factor or the other. The total gauge is the product, in this order:
!>
!>       U(k) = u_opt(k) . u_mlwf(k)          (num_bands x num_wann)
!>
!> Multiplying them is one line for the reader; separating them again afterwards is not
!> possible. Nothing else is interpreted here -- no reordering, no units, no convention
!> beyond the shapes the attributes state.
MODULE m_wannierlib_export_gauge
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_gauge
CONTAINS

   SUBROUTINE wannierlib_export_gauge(kpts, fmpi, jspin, u_opt, u_matrix)
      USE m_types_kpts
      USE m_types_mpi
#ifdef CPP_HDF
      USE hdf5
      USE m_hdf_tools
#endif
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: jspin
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)      !> (num_bands, num_wann, nkptf)
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)   !> (num_wann,  num_wann, nkptf)
#ifdef CPP_HDF
      INTEGER(HID_T) :: fid, gid
      INTEGER :: err, nb, nw, ik, i, j, n(4)
      REAL, ALLOCATABLE :: buf(:, :, :, :)
      CHARACTER(LEN=64) :: filename
      LOGICAL :: l_ex

      !> Rank 0 only: w90 runs in library mode on the master, where both factors are whole.
      !> Writing from every rank would have each of them truncate the same file.
      IF (fmpi%irank /= 0) RETURN

      nb = SIZE(u_opt, 1); nw = SIZE(u_opt, 2)
      IF (SIZE(u_opt, 3) < kpts%nkptf .OR. SIZE(u_matrix, 3) < kpts%nkptf) &
         CALL juDFT_error("wannierlib: the gauge does not cover the full k-set", &
                          calledby="wannierlib_export_gauge")

      WRITE (filename, '(a,i0,a)') 'WF', jspin, '_gauge.hdf'
      INQUIRE (file=TRIM(filename), exist=l_ex)
      IF (l_ex) CALL system('rm '//TRIM(filename))
      CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, fid, err, H5P_DEFAULT_F, H5P_DEFAULT_F)
      CALL h5gopen_f(fid, '/', gid, err)

      CALL io_write_attint0(gid, 'version', 1)
      CALL io_write_attint0(gid, 'nkpt', kpts%nkptf)
      CALL io_write_attint0(gid, 'num_bands', nb)
      CALL io_write_attint0(gid, 'num_wann', nw)
      CALL io_write_attint0(gid, 'spin', jspin)

      ALLOCATE (buf(2, nb, nw, kpts%nkptf))
      DO ik = 1, kpts%nkptf
         DO j = 1, nw
            DO i = 1, nb
               buf(1, i, j, ik) = REAL(u_opt(i, j, ik))
               buf(2, i, j, ik) = AIMAG(u_opt(i, j, ik))
            END DO
         END DO
      END DO
      n = (/2, nb, nw, kpts%nkptf/)
      CALL wr_r4(gid, 'u_opt', n, buf)
      DEALLOCATE (buf)

      ALLOCATE (buf(2, nw, nw, kpts%nkptf))
      DO ik = 1, kpts%nkptf
         DO j = 1, nw
            DO i = 1, nw
               buf(1, i, j, ik) = REAL(u_matrix(i, j, ik))
               buf(2, i, j, ik) = AIMAG(u_matrix(i, j, ik))
            END DO
         END DO
      END DO
      n = (/2, nw, nw, kpts%nkptf/)
      CALL wr_r4(gid, 'u_mlwf', n, buf)
      DEALLOCATE (buf)

      CALL h5gclose_f(gid, err)
      CALL h5fclose_f(fid, err)
      WRITE (oUnit, '(a,i0,a,i0,a,i0,a)') 'wannierlib: wrote '//TRIM(filename)//' (', &
         kpts%nkptf, ' k-points, ', nb, ' bands -> ', nw, &
         ' Wannier functions; U = u_opt . u_mlwf)'
#else
      CALL juDFT_error("wannierlib: the gauge export needs FLEUR built with HDF5", &
                       calledby="wannierlib_export_gauge")
#endif
   END SUBROUTINE wannierlib_export_gauge

#ifdef CPP_HDF
   !> Dataspace, dataset, write, close -- the same four steps every time.
   SUBROUTINE wr_r4(gid, name, n, dat)
      USE hdf5
      USE m_hdf_tools
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:)
      REAL, INTENT(IN) :: dat(:, :, :, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(4, INT(n(:4), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_real4(did, (/1, 1, 1, 1/), n(:4), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_r4
#endif

END MODULE m_wannierlib_export_gauge
