!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Read WF<spin>_operators_ks.hdf back into the coarse operator slices it was written from.
!>
!> The file holds the operator matrices in the Bloch basis, which is the one thing a later
!> run cannot rebuild without the eigenvectors. Everything else a gauge run consumes is
!> either in the Wannier90 checkpoint (the gauge and the overlaps) or costs nothing to
!> recompute, so this reader is what makes it possible to resume from the operator stage.
!>
!> WHAT IS CHECKED. The layout version, the number of k-points, and the band, atom and
!> channel counts of the dims record, each against this run. A file written for a different
!> window or a different structure is refused by number rather than discovered later as a
!> shape mismatch. What NO check here can see is the density: a file made on the same mesh
!> from a different SCF matches on every count above and its matrices belong to states this
!> run never computed. That is the same blind spot the Wannier90 checkpoint has.
!>
!> WHAT IS NOT IN THE FILE. soc0, the sum over the four spin-orbit blocks, is not stored
!> because it is SUM(soc4, DIM=3) and the blocks are. The cross-spin blocks x0 and tx0 are
!> not stored either: they only exist in a collinear two-channel run, which the exporter
!> refuses outright.
!>
!> DISTRIBUTION. Mirror of the exporter: rank 0 holds the file and reads one component over
!> the whole mesh, it is broadcast, and every rank keeps the k-points distk gives it. The
!> full-mesh buffer therefore exists once per rank, as it already does on the way out.
MODULE m_wannierlib_import_ks
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_mpi, ONLY: t_mpi
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_import_ks, wannierlib_ks_roundtrip
CONTAINS

   !> Fill the coarse slices from the file. l_found comes back .FALSE. when there is no
   !> file at all, which is a question for the caller and not an error here; anything else
   !> that is wrong with it stops the run.
   SUBROUTINE wannierlib_import_ks(l_spin, l_orbital, l_spinorbit, l_torque, &
                                   s0, l0, soc4, t0, n_channels, nat, &
                                   distk, nkptf, nb, jspin, fmpi, l_found)
#ifdef CPP_HDF
      USE hdf5
      USE m_wannierlib_hdf_util
#endif
#ifdef CPP_MPI
      USE mpi
#endif
      LOGICAL,              INTENT(IN)  :: l_spin, l_orbital, l_spinorbit, l_torque
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: s0(:, :, :, :)          !< (nb,nb,3,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: l0(:, :, :, :, :, :)    !< (nb,nb,3,nat,channel,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: soc4(:, :, :, :)        !< (nb,nb,4,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: t0(:, :, :, :)          !< (nb,nb,3,nk_loc)
      INTEGER,              INTENT(IN)  :: n_channels, nat
      INTEGER,              INTENT(IN)  :: distk(:)                !< (nkptf) owning rank of each k
      INTEGER,              INTENT(IN)  :: nkptf, nb, jspin
      TYPE(t_mpi),          INTENT(IN)  :: fmpi
      LOGICAL,              INTENT(OUT) :: l_found
#ifdef CPP_HDF
      INTEGER(HID_T)    :: fid, gid, did
      CHARACTER(LEN=80) :: filename
      INTEGER           :: hdr(7), dims(4), n5(5)
      INTEGER           :: ic, iat, ich, icomp, ncomp, nkl, e, ierr
      LOGICAL           :: l_root

      l_root = (fmpi%irank == 0)
      fid = -1; gid = -1; did = -1
      hdr = 0

      !> Rank 0 reads the header; every rank judges it. Having one rank decide and the
      !> others carry on is how half a communicator ends up waiting for the other half.
      IF (l_root) THEN
         CALL wl_hdf_open('operators_ks', jspin, fid, filename, l_found)
         IF (l_found) THEN
            CALL wl_hdf_root_read(fid, hdr(2), hdr(3), gid)
            IF (.NOT. wl_ds_exists(gid, 'dims')) &
               CALL juDFT_error("wannierlib: "//TRIM(filename)//" has no dims record", &
                                calledby="wannierlib_import_ks")
            CALL rd_i1(gid, 'dims', [4], dims)
            hdr(1) = 1
            hdr(4:7) = dims
         END IF
      END IF
#ifdef CPP_MPI
      IF (fmpi%isize > 1) CALL MPI_BCAST(hdr, SIZE(hdr), MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
#endif
      l_found = (hdr(1) == 1)
      IF (.NOT. l_found) RETURN

      IF (hdr(2) /= 1) CALL juDFT_error("wannierlib: the Bloch-basis operator file has "// &
                                        "layout version "//i2s(hdr(2))//", not 1", &
                                        calledby="wannierlib_import_ks")
      !> The number of k-points is written TWICE, in the root attribute and in the dims
      !> record. Both are checked, separately and under different names, so that the day
      !> they disagree with each other the message can say which of the two it was.
      CALL assert_same('k-points in its root attribute', hdr(3), nkptf)
      CALL assert_same('k-points in its dims record',    hdr(5), nkptf)
      CALL assert_same('bands',     hdr(4), nb)
      CALL assert_same('atoms',     hdr(6), nat)
      CALL assert_same('channels',  hdr(7), n_channels)

      nkl = MAX(1, COUNT(distk == fmpi%irank))

      IF (l_spin) THEN
         CALL open_ds(gid, 'spin', l_root, did, n5)
         ncomp = want_ncomp(n5, 3, 'spin', l_root)
         ALLOCATE (s0(nb, nb, ncomp, nkl), source=CMPLX(0.0, 0.0))
         DO ic = 1, ncomp
            CALL get_component(did, 'spin', ic, l_root, fmpi, distk, nkptf, nb, s0(:, :, ic, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_orbital) THEN
         CALL open_ds(gid, 'orbital', l_root, did, n5)
         ncomp = want_ncomp(n5, 3*nat*n_channels, 'orbital', l_root)
         ALLOCATE (l0(nb, nb, 3, nat, n_channels, nkl), source=CMPLX(0.0, 0.0))
         !> The exporter numbers the components with the cartesian index fastest, then the
         !> atom, then the channel. Undone here in the same order and nowhere else, so the
         !> two loops can only disagree by being changed together.
         icomp = 0
         DO ich = 1, n_channels
            DO iat = 1, nat
               DO ic = 1, 3
                  icomp = icomp + 1
                  CALL get_component(did, 'orbital', icomp, l_root, fmpi, distk, nkptf, nb, &
                                     l0(:, :, ic, iat, ich, :))
               END DO
            END DO
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_spinorbit) THEN
         CALL open_ds(gid, 'spinOrbit', l_root, did, n5)
         ncomp = want_ncomp(n5, 4, 'spinOrbit', l_root)
         ALLOCATE (soc4(nb, nb, ncomp, nkl), source=CMPLX(0.0, 0.0))
         DO ic = 1, ncomp
            CALL get_component(did, 'spinOrbit', ic, l_root, fmpi, distk, nkptf, nb, soc4(:, :, ic, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_torque) THEN
         CALL open_ds(gid, 'exchangeTorque', l_root, did, n5)
         ncomp = want_ncomp(n5, 3, 'exchangeTorque', l_root)
         ALLOCATE (t0(nb, nb, ncomp, nkl), source=CMPLX(0.0, 0.0))
         DO ic = 1, ncomp
            CALL get_component(did, 'exchangeTorque', ic, l_root, fmpi, distk, nkptf, nb, t0(:, :, ic, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_root) THEN
         CALL h5gclose_f(gid, e)
         CALL h5fclose_f(fid, e)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: read '//TRIM(filename)//' (', nb, &
            ' bands, ', nkptf, ' k-points, no gauge applied)'
      END IF
#else
      l_found = .FALSE.
      CALL juDFT_error("wannierlib: the Bloch-basis operator import needs FLEUR built "// &
                       "with HDF5", calledby="wannierlib_import_ks")
#endif
   END SUBROUTINE wannierlib_import_ks

   !> Read the file back and require it to reproduce, exactly, the matrices it was written
   !> from. A zero return code is not evidence of having written: a full quota, a short
   !> write or a component sent to the wrong slab all leave a file that opens cleanly and
   !> answers a different question, and the run that reads it next has no way to tell.
   !> Exact is the right demand -- both sides are the same IEEE doubles through the same
   !> slab calls -- so any non-zero residual is a defect and not a tolerance to widen.
   SUBROUTINE wannierlib_ks_roundtrip(l_spin, l_orbital, l_spinorbit, l_torque, &
                                      s0, l0, soc4, t0, n_channels, nat, &
                                      distk, nkptf, nb, jspin, fmpi)
#ifdef CPP_MPI
      USE mpi
#endif
      LOGICAL,              INTENT(IN) :: l_spin, l_orbital, l_spinorbit, l_torque
      COMPLEX, ALLOCATABLE, INTENT(IN) :: s0(:, :, :, :)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: l0(:, :, :, :, :, :)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: soc4(:, :, :, :)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: t0(:, :, :, :)
      INTEGER,              INTENT(IN) :: n_channels, nat, distk(:), nkptf, nb, jspin
      TYPE(t_mpi),          INTENT(IN) :: fmpi

      COMPLEX, ALLOCATABLE :: r_s0(:, :, :, :), r_soc4(:, :, :, :), r_t0(:, :, :, :)
      COMPLEX, ALLOCATABLE :: r_l0(:, :, :, :, :, :)
      REAL    :: dmax
      LOGICAL :: l_found

      CALL wannierlib_import_ks(l_spin, l_orbital, l_spinorbit, l_torque, &
                                r_s0, r_l0, r_soc4, r_t0, n_channels, nat, &
                                distk, nkptf, nb, jspin, fmpi, l_found)
      IF (.NOT. l_found) CALL juDFT_error("wannierlib: the file just written is not there", &
                                          calledby="wannierlib_ks_roundtrip")
      !> Asked for and not there on one of the two sides is a defect of its own, and one
      !> that would otherwise leave the operator silently unchecked.
      CALL assert_pair('spin',           l_spin,      ALLOCATED(s0),   ALLOCATED(r_s0))
      CALL assert_pair('orbital',        l_orbital,   ALLOCATED(l0),   ALLOCATED(r_l0))
      CALL assert_pair('spinOrbit',      l_spinorbit, ALLOCATED(soc4), ALLOCATED(r_soc4))
      CALL assert_pair('exchangeTorque', l_torque,    ALLOCATED(t0),   ALLOCATED(r_t0))

      dmax = 0.0
      IF (l_spin)      CALL one_op('spin',    MAXVAL(ABS(r_s0 - s0)),     MAXVAL(ABS(s0)),   fmpi, dmax)
      IF (l_orbital)   CALL one_op('orbital', MAXVAL(ABS(r_l0 - l0)),     MAXVAL(ABS(l0)),   fmpi, dmax)
      IF (l_spinorbit) CALL one_op('soc',     MAXVAL(ABS(r_soc4 - soc4)), MAXVAL(ABS(soc4)), fmpi, dmax)
      IF (l_torque)    CALL one_op('torque',  MAXVAL(ABS(r_t0 - t0)),     MAXVAL(ABS(t0)),   fmpi, dmax)

      !> dmax is already the worst over the ranks: one_op reduces each operator before it
      !> adds it in. Every rank therefore reaches the same verdict, which is what keeps
      !> the abort from being taken by one rank while the others carry on.
      IF (dmax /= 0.0) CALL juDFT_error("wannierlib: the Bloch-basis operator file does "// &
                                        "not read back as it was written", &
                                        calledby="wannierlib_ks_roundtrip")
   END SUBROUTINE wannierlib_ks_roundtrip

   !> One operator: the residual, and the size of what it was compared against. The
   !> second number is not decoration. An operator that was never asked for, or that came
   !> back empty on both sides, also has a residual of zero, and a line that printed only
   !> the zero would read exactly like a check that passed -- a test that cannot fail is
   !> not a test. With the scale beside it, a zero residual over a zero scale is visibly
   !> nothing having been compared.
   SUBROUTINE one_op(name, resid, scale, fmpi, dmax)
#ifdef CPP_MPI
      USE mpi
#endif
      CHARACTER(LEN=*), INTENT(IN)    :: name
      REAL,             INTENT(IN)    :: resid, scale
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      REAL,             INTENT(INOUT) :: dmax
      REAL    :: r, s
      INTEGER :: ierr
      !> Reduced BEFORE it is printed, not after. Each rank only ever sees its own
      !> k-points, so a line written by rank 0 with rank 0's number reads as a clean zero
      !> while the run aborts on a mismatch that happened on another rank -- which is the
      !> one case where the reader of the log most needs the number to be the true one.
      r = resid; s = scale
#ifdef CPP_MPI
      IF (fmpi%isize > 1) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, r, 1, MPI_DOUBLE_PRECISION, MPI_MAX, fmpi%mpi_comm, ierr)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, s, 1, MPI_DOUBLE_PRECISION, MPI_MAX, fmpi%mpi_comm, ierr)
      END IF
#endif
      dmax = MAX(dmax, r)
      IF (fmpi%irank /= 0) RETURN
      !> s, not scale: the verdict rests on the maximum over the ranks. With the local
      !> one, a rank 0 whose own slice happens to be zero announces that nothing was
      !> compared while the others were comparing -- the very thing the reduction fixes.
      IF (s == 0.0) THEN
         WRITE (oUnit, '(a,es10.2,a)') &
            'wannierlib: operators_ks round trip, '//name//': |read - written| ', r, &
            '   but the operator is identically zero: NOTHING WAS COMPARED'
      ELSE
         WRITE (oUnit, '(a,es10.2,a,es10.2)') &
            'wannierlib: operators_ks round trip, '//name//': |read - written| ', r, &
            '   over a largest element of ', s
      END IF
   END SUBROUTINE one_op

#ifdef CPP_HDF
   !> One component over the whole mesh: read on rank 0, broadcast, and kept where distk
   !> says. Zeroed first, because a rank that owns no k-point still has an array and an
   !> INTENT(OUT) one would otherwise be whatever was in that memory.
   SUBROUTINE get_component(did, name, icomp, l_root, fmpi, distk, nkptf, nb, loc)
      USE hdf5
      USE m_wannierlib_hdf_util, ONLY: wl_ds5_read
#ifdef CPP_MPI
      USE mpi
#endif
      INTEGER(HID_T),   INTENT(IN)  :: did
      CHARACTER(LEN=*), INTENT(IN)  :: name
      INTEGER,          INTENT(IN)  :: icomp, nkptf, nb
      LOGICAL,          INTENT(IN)  :: l_root
      TYPE(t_mpi),      INTENT(IN)  :: fmpi
      INTEGER,          INTENT(IN)  :: distk(:)
      COMPLEX,          INTENT(OUT) :: loc(:, :, :)     !< (nb,nb,nk_loc)
      COMPLEX, ALLOCATABLE :: full(:, :, :)
      REAL,    ALLOCATABLE :: buf(:, :, :, :, :)
      INTEGER :: ikpt, ik_local, ierr

      loc = CMPLX(0.0, 0.0)
      ALLOCATE (full(nb, nb, nkptf), source=CMPLX(0.0, 0.0))
      IF (l_root) THEN
         ALLOCATE (buf(2, nb, nb, 1, nkptf))
         CALL wl_ds5_read(did, name, [1, 1, 1, icomp, 1], [2, nb, nb, 1, nkptf], buf)
         full = CMPLX(buf(1, :, :, 1, :), buf(2, :, :, 1, :))
         DEALLOCATE (buf)
      END IF
#ifdef CPP_MPI
      IF (fmpi%isize > 1) CALL MPI_BCAST(full, SIZE(full), MPI_DOUBLE_COMPLEX, 0, &
                                         fmpi%mpi_comm, ierr)
#endif
      ik_local = 0
      DO ikpt = 1, nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE
         ik_local = ik_local + 1
         IF (ik_local > SIZE(loc, 3)) &
            CALL juDFT_error("wannierlib: more k-points assigned to this rank than the "// &
                             "operator holds", calledby="get_component")
         loc(:, :, ik_local) = full(:, :, ikpt)
      END DO
   END SUBROUTINE get_component

   !> Open a named operator dataset on rank 0 and take its shape. Absence stops the run:
   !> the caller only asks for what inp.xml named, so a missing dataset means the file was
   !> written by a run that was asked for less.
   SUBROUTINE open_ds(gid, name, l_root, did, n5)
      USE hdf5
      USE m_wannierlib_hdf_util, ONLY: wl_ds5_open_read, wl_ds5_dims, wl_ds_exists
      INTEGER(HID_T),   INTENT(IN)  :: gid
      CHARACTER(LEN=*), INTENT(IN)  :: name
      LOGICAL,          INTENT(IN)  :: l_root
      INTEGER(HID_T),   INTENT(OUT) :: did
      INTEGER,          INTENT(OUT) :: n5(5)
      did = -1
      n5 = 0
      IF (.NOT. l_root) RETURN
      IF (.NOT. wl_ds_exists(gid, name)) &
         CALL juDFT_error('wannierlib: the Bloch-basis operator file has no "'//name// &
                          '"; it was written by a run that was not asked for it', &
                          calledby="wannierlib_import_ks")
      CALL wl_ds5_open_read(gid, name, did)
      CALL wl_ds5_dims(did, n5)
   END SUBROUTINE open_ds

   !> The component count comes from the file, and has to be the one this run would have
   !> written. Only rank 0 has the shape, so the judgement is made there and the number is
   !> the expected one everywhere else -- the ranks that do not read cannot disagree with
   !> a file they have not opened.
   INTEGER FUNCTION want_ncomp(n5, expect, name, l_root)
      INTEGER,          INTENT(IN) :: n5(5), expect
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL,          INTENT(IN) :: l_root
      want_ncomp = expect
      IF (.NOT. l_root) RETURN
      IF (n5(4) == expect) RETURN
      CALL juDFT_error('wannierlib: "'//name//'" in the Bloch-basis operator file has '// &
                       i2s(n5(4))//' components, this run expects '//i2s(expect), &
                       calledby="wannierlib_import_ks")
   END FUNCTION want_ncomp
#endif

   !> An operator that was asked for has to be on both sides of the comparison.
   SUBROUTINE assert_pair(name, l_asked, l_mem, l_file)
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL,          INTENT(IN) :: l_asked, l_mem, l_file
      IF (.NOT. l_asked) RETURN
      IF (l_mem .AND. l_file) RETURN
      CALL juDFT_error('wannierlib: "'//name//'" was asked for but is missing '// &
                       MERGE('in memory', 'in a file', .NOT. l_mem)//', so the round '// &
                       'trip would have passed without comparing it', &
                       calledby="wannierlib_ks_roundtrip")
   END SUBROUTINE assert_pair

   !> One count from the file against the same count in this run.
   SUBROUTINE assert_same(what, in_file, in_run)
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER,          INTENT(IN) :: in_file, in_run
      IF (in_file == in_run) RETURN
      CALL juDFT_error('wannierlib: the Bloch-basis operator file was written with '// &
                       i2s(in_file)//' '//what//', this run has '//i2s(in_run), &
                       hint='it belongs to another window, mesh or structure', &
                       calledby="wannierlib_import_ks")
   END SUBROUTINE assert_same

   !> A small integer as text, for messages that have to carry the two numbers that
   !> disagree: a message that says only that something does not match sends the reader
   !> back to the file to find out what.
   FUNCTION i2s(i) RESULT(s)
      INTEGER, INTENT(IN) :: i
      CHARACTER(LEN=:), ALLOCATABLE :: s
      CHARACTER(LEN=20) :: b
      WRITE (b, '(i0)') i
      s = TRIM(b)
   END FUNCTION i2s

END MODULE m_wannierlib_import_ks
