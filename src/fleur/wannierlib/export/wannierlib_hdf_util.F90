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
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wr_r1, wr_r2, wr_r3, wr_r4, wr_i1, wr_i2, wr_i3, wl_hdf_create, wl_hdf_root
   PUBLIC :: wl_ds5_open, wl_ds5_slab, wl_ds_close
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

#else
   IMPLICIT NONE
   PRIVATE
#endif
END MODULE m_wannierlib_hdf_util
