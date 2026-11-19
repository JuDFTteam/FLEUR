!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Dataspace, dataset, write, close: the same four steps every array needs, one call per
!> rank and type, so a caller states what it writes and not how HDF5 is driven. And the two
!> steps that come before any of them: making the file and opening its root group.
!>
!> The file is WF<spin>_<what>.hdf and its root carries a layout version and the number of
!> k-points. That is the part of the contract the readers share, so it lives here rather
!> than being spelled out once per exporter.
!>
!> Present only in an HDF5 build. A caller reaches these from inside its own CPP_HDF
!> branch, which is where it has a group to write into at all.
MODULE m_wannierlib_hdf_util
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wr_r1, wr_r2, wr_r3, wr_r4, wr_i1, wr_i2, wr_i3, wl_hdf_create, wl_hdf_root
   PUBLIC :: wl_ds5_open, wl_ds5_slab, wl_ds_close
   PUBLIC :: wl_hdf_name, wl_hdf_open, wl_hdf_root_read, rd_i1, wl_ds_exists
   PUBLIC :: wl_ds5_open_read, wl_ds5_read, wl_ds5_dims, rd_r2
CONTAINS

   !> Make WF<jspin>_<what>.hdf and return it open. A stale file is removed rather than
   !> written over, so a shorter run cannot leave the tail of a longer one behind.
   SUBROUTINE wl_hdf_create(what, jspin, fid, filename)
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN) :: jspin
      INTEGER(HID_T), INTENT(OUT) :: fid
      CHARACTER(LEN=*), INTENT(OUT) :: filename
      LOGICAL :: l_ex
      INTEGER :: e
      WRITE (filename, '(a,i0,a)') 'WF', jspin, '_'//TRIM(what)//'.hdf'
      INQUIRE (file=TRIM(filename), exist=l_ex)
      IF (l_ex) CALL system('rm '//TRIM(filename))
      CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, fid, e, H5P_DEFAULT_F, H5P_DEFAULT_F)
   END SUBROUTINE wl_hdf_create

   !> The root group with the two attributes a reader needs before it can read anything
   !> else. Returned OPEN: what else belongs on the root is the exporter's own business,
   !> and it closes the group when it has written it.
   SUBROUTINE wl_hdf_root(fid, version, nkpt, gid)
      INTEGER(HID_T), INTENT(IN) :: fid
      INTEGER, INTENT(IN) :: version, nkpt
      INTEGER(HID_T), INTENT(OUT) :: gid
      INTEGER :: e
      CALL h5gopen_f(fid, '/', gid, e)
      CALL io_write_attint0(gid, 'version', version)
      CALL io_write_attint0(gid, 'nkpt', nkpt)
   END SUBROUTINE wl_hdf_root

   !> A rank-5 dataset filled one slab at a time. The arrays this serves are an operator
   !> over every band pair, component and k-point at once, and they do not fit beside the
   !> buffer they are gathered through, so the caller opens the dataset, writes a component
   !> per call and closes it, instead of handing over a whole array.
   SUBROUTINE wl_ds5_open(gid, name, n, did)
      INTEGER(HID_T),   INTENT(IN)  :: gid
      CHARACTER(LEN=*), INTENT(IN)  :: name
      INTEGER,          INTENT(IN)  :: n(5)
      INTEGER(HID_T),   INTENT(OUT) :: did
      INTEGER(HID_T) :: sid
      INTEGER :: e
      CALL h5screate_simple_f(5, INT(n, HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, e)
      CALL h5sclose_f(sid, e)
   END SUBROUTINE wl_ds5_open

   SUBROUTINE wl_ds5_slab(did, name, start, cnt, dat)
      INTEGER(HID_T),   INTENT(IN) :: did
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER,          INTENT(IN) :: start(5), cnt(5)
      REAL,             INTENT(IN) :: dat(:, :, :, :, :)
      CALL io_write_real5(did, start, cnt, name, dat)
   END SUBROUTINE wl_ds5_slab

   SUBROUTINE wl_ds_close(did)
      INTEGER(HID_T), INTENT(IN) :: did
      INTEGER :: e
      CALL h5dclose_f(did, e)
   END SUBROUTINE wl_ds_close

   SUBROUTINE wr_r2(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:)
      REAL, INTENT(IN) :: dat(:, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(2, INT(n(:2), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_real2(did, (/1, 1/), n(:2), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_r2

   SUBROUTINE wr_r1(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:)
      REAL, INTENT(IN) :: dat(:)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(1, INT(n(:1), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_real1(did, (/1/), n(:1), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_r1

   SUBROUTINE wr_r3(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:)
      REAL, INTENT(IN) :: dat(:, :, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(3, INT(n(:3), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_real3(did, (/1, 1, 1/), n(:3), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_r3

   SUBROUTINE wr_r4(gid, name, n, dat)
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

   SUBROUTINE wr_i1(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:), dat(:)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(1, INT(n(:1), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_integer1(did, (/1/), n(:1), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_i1

   SUBROUTINE wr_i2(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:), dat(:, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(2, INT(n(:2), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_integer2(did, (/1, 1/), n(:2), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_i2

   SUBROUTINE wr_i3(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:), dat(:, :, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(3, INT(n(:3), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_integer3(did, (/1, 1, 1/), n(:3), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_i3

   !> ---- reading back what the four above wrote --------------------------------------
   !>
   !> Symmetric on purpose: a file written through wl_ds5_open/wl_ds5_slab is read through
   !> wl_ds5_open_read/wl_ds5_read with the SAME start and count, and the round trip is the
   !> identity. In particular there is no transpose to undo here. The warning in the
   !> operator exporters is about a reader that sees HDF's own C ordering -- python, or
   !> anything not going through this module; a Fortran reader handing back the shape the
   !> writer declared gets the matrix the way it was stored. Transposing "to be safe" is
   !> the silent bug, not the cure: a Hermitian matrix transposed is still Hermitian and
   !> keeps its norm and its eigenvalues.

   !> Open WF<jspin>_<what>.hdf for reading. A missing file is NOT an error here: the
   !> caller may know how to produce what it holds, and only it can decide.
   !> The one place the file name is spelled, for the writer and the reader alike:
   !> two spellings of one convention drift apart without anything noticing.
   FUNCTION wl_hdf_name(what, jspin) RESULT(filename)
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER,          INTENT(IN) :: jspin
      CHARACTER(LEN=80) :: filename
      WRITE (filename, '(a,i0,a)') 'WF', jspin, '_'//TRIM(what)//'.hdf'
   END FUNCTION wl_hdf_name

   SUBROUTINE wl_hdf_open(what, jspin, fid, filename, l_found)
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN) :: jspin
      INTEGER(HID_T), INTENT(OUT) :: fid
      CHARACTER(LEN=*), INTENT(OUT) :: filename
      LOGICAL, INTENT(OUT) :: l_found
      INTEGER :: e
      fid = -1
      filename = wl_hdf_name(what, jspin)
      INQUIRE (file=TRIM(filename), exist=l_found)
      IF (.NOT. l_found) RETURN
      CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, fid, e)
      !> A file that is there and does not open is not the same thing as a file that is
      !> not there. Dropping the status here hands the caller an invalid handle together
      !> with l_found true, and the run dies several calls later inside HDF5 instead.
      IF (e < 0) CALL juDFT_error("wannierlib: "//TRIM(filename)//" exists but will not "// &
                                  "open as HDF5", hint='it may be truncated, or written '// &
                                  'by a run that was killed', calledby="wl_hdf_open")
   END SUBROUTINE wl_hdf_open

   !> The two attributes wl_hdf_root wrote. Read before anything else, because what they
   !> say decides whether the rest of the file means what this run thinks it means.
   SUBROUTINE wl_hdf_root_read(fid, version, nkpt, gid)
      INTEGER(HID_T), INTENT(IN) :: fid
      INTEGER, INTENT(OUT) :: version, nkpt
      INTEGER(HID_T), INTENT(OUT) :: gid
      INTEGER :: e
      CALL h5gopen_f(fid, '/', gid, e)
      IF (e < 0) CALL juDFT_error("wannierlib: the operator file has no root group", &
                                  calledby="wl_hdf_root_read")
      CALL io_read_attint0(gid, 'version', version)
      CALL io_read_attint0(gid, 'nkpt', nkpt)
   END SUBROUTINE wl_hdf_root_read

   LOGICAL FUNCTION wl_ds_exists(gid, name)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      wl_ds_exists = io_dataexists(gid, name)
   END FUNCTION wl_ds_exists

   SUBROUTINE rd_i1(gid, name, n, dat)
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:)
      INTEGER, INTENT(OUT) :: dat(:)
      INTEGER(HID_T) :: did
      INTEGER :: e
      CALL h5dopen_f(gid, name, did, e)
      CALL io_read_integer1(did, (/1/), n(:1), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE rd_i1

   !> The counterpart of wr_r2, with the same shape in the same order.
   SUBROUTINE rd_r2(gid, name, n, dat)
      INTEGER(HID_T),   INTENT(IN)  :: gid
      CHARACTER(LEN=*), INTENT(IN)  :: name
      INTEGER,          INTENT(IN)  :: n(:)
      REAL,             INTENT(OUT) :: dat(:, :)
      INTEGER(HID_T) :: did
      INTEGER :: e
      CALL h5dopen_f(gid, name, did, e)
      IF (e < 0) CALL juDFT_error('wannierlib: "'//name//'" will not open', &
                                  calledby="rd_r2")
      CALL io_read_real2(did, (/1, 1/), n(:2), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE rd_r2

   SUBROUTINE wl_ds5_open_read(gid, name, did)
      INTEGER(HID_T),   INTENT(IN)  :: gid
      CHARACTER(LEN=*), INTENT(IN)  :: name
      INTEGER(HID_T),   INTENT(OUT) :: did
      INTEGER :: e
      CALL h5dopen_f(gid, name, did, e)
      IF (e < 0) CALL juDFT_error('wannierlib: "'//name//'" will not open', &
                                  calledby="wl_ds5_open_read")
   END SUBROUTINE wl_ds5_open_read

   !> The shape the dataset was created with, in the SAME order wl_ds5_open was given it:
   !> the Fortran interface does not reverse it. So n(4) is the component count, and a
   !> reader does not have to be told it separately -- which is what keeps the component
   !> count of an operator a property of the file and not of the caller's expectation.
   SUBROUTINE wl_ds5_dims(did, n)
      INTEGER(HID_T), INTENT(IN)  :: did
      INTEGER,        INTENT(OUT) :: n(5)
      INTEGER(HID_T)   :: sid
      INTEGER(HSIZE_T) :: d(5), m(5)
      INTEGER :: e
      CALL h5dget_space_f(did, sid, e)
      !> This one returns the RANK in its status, not zero. So the rank comes for free,
      !> and taking it means n is never read past what the dataspace filled: a rank-4
      !> dataset would otherwise leave n(5) as whatever was on the stack.
      CALL h5sget_simple_extent_dims_f(sid, d, m, e)
      IF (e /= 5) CALL juDFT_error("wannierlib: an operator dataset has rank "// &
                                   int2str(e)//", not 5", calledby="wl_ds5_dims")
      CALL h5sclose_f(sid, e)
      n = INT(d)
   END SUBROUTINE wl_ds5_dims

   SUBROUTINE wl_ds5_read(did, name, start, cnt, dat)
      INTEGER(HID_T),   INTENT(IN)  :: did
      CHARACTER(LEN=*), INTENT(IN)  :: name
      INTEGER,          INTENT(IN)  :: start(5), cnt(5)
      REAL,             INTENT(OUT) :: dat(:, :, :, :, :)
      CALL io_read_real5(did, start, cnt, name, dat)
   END SUBROUTINE wl_ds5_read

#else
   IMPLICIT NONE
   PRIVATE
#endif
END MODULE m_wannierlib_hdf_util
