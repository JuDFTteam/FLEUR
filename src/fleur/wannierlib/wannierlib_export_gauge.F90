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
   USE m_types, ONLY: t_results
   USE m_types_kpts
   USE m_types_input
   USE m_types_mpi
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_gauge, &
             wannierlib_store_gauge
CONTAINS

   SUBROUTINE wannierlib_export_gauge(kpts, fmpi, jspin, u_opt, u_matrix)
      USE m_types_kpts
      USE m_types_mpi
#ifdef CPP_HDF
      USE hdf5
      USE m_hdf_tools
      USE m_wannierlib_hdf_util, ONLY: wr_r4
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


   !> Copy the Wannier gauge into results, completed over the whole mesh.
   !>
   !> NOT COVERED BY ANY TEST, and it cannot be from here: nothing in the wannier testset
   !> reads results%U_mat, and the only consumer is interpolate_bandstructure, reached from
   !> the electron-phonon path in dfpt/superconductivity. Both of its call sites ask it not
   !> to write, so there is no file to compare either. What can be checked without running
   !> that path: <export gauge="T"/> writes the same two matrices to WF<n>_gauge.hdf, so the
   !> values this stores are the ones in that file.
   !>
   !> The electron-phonon interpolation reads results%U_mat and results%U_dis; what run_w90
   !> returns is each rank's own k-points with the rest left at zero, which is what the
   !> operator pass expects, so the sum that completes the mesh is done on a copy and the
   !> arrays the caller holds are left exactly as they were.
   SUBROUTINE wannierlib_store_gauge(fmpi, distk, kpts, input, num_bands, num_wann, jspin, &
                                     u_matrix, u_opt, results)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_input), INTENT(IN) :: input
      INTEGER, INTENT(IN) :: num_bands, num_wann, jspin
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (nw, nw, nk), this rank's k-points
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (nb, nw, nk), this rank's k-points
      TYPE(t_results), INTENT(INOUT) :: results

      COMPLEX, ALLOCATABLE :: buf(:, :, :)
      INTEGER :: ikpt, ierr

      IF (SIZE(u_matrix, 3) == kpts%nkptf) THEN
         buf = u_matrix
         DO ikpt = 1, kpts%nkptf
            IF (distk(ikpt) /= fmpi%irank) buf(:, :, ikpt) = CMPLX(0.0, 0.0)
         END DO
#ifdef CPP_MPI
         IF (fmpi%isize > 1) &
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, buf, SIZE(buf), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                               fmpi%mpi_comm, ierr)
#endif
         IF (.NOT. ALLOCATED(results%U_mat)) THEN
            ALLOCATE (results%U_mat(num_wann, num_wann, kpts%nkptf, input%jspins), &
                      source=CMPLX(0.0, 0.0))
         END IF
         results%U_mat(:, :, :, jspin) = buf
         DEALLOCATE (buf)
      END IF

      !> Only meaningful when the subspace was chosen, i.e. more bands than functions;
      !> with num_bands == num_wann u_opt is the projection itself and carries nothing.
      !> That is the condition the routine this replaces used, expressed as a test because
      !> here u_opt is always allocated.
      IF (num_bands > num_wann .AND. SIZE(u_opt, 3) == kpts%nkptf) THEN
         buf = u_opt
         DO ikpt = 1, kpts%nkptf
            IF (distk(ikpt) /= fmpi%irank) buf(:, :, ikpt) = CMPLX(0.0, 0.0)
         END DO
#ifdef CPP_MPI
         IF (fmpi%isize > 1) &
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, buf, SIZE(buf), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                               fmpi%mpi_comm, ierr)
#endif
         IF (.NOT. ALLOCATED(results%U_dis)) THEN
            ALLOCATE (results%U_dis(num_bands, num_wann, kpts%nkptf, input%jspins), &
                      source=CMPLX(0.0, 0.0))
         END IF
         results%U_dis(:, :, :, jspin) = buf
         DEALLOCATE (buf)
      END IF
   END SUBROUTINE wannierlib_store_gauge
END MODULE m_wannierlib_export_gauge
