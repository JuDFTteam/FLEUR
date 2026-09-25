!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Dataspace, dataset, write, close: the same four steps every array needs, one call per
!> rank and type, so a caller states what it writes and not how HDF5 is driven.
!>
!> Present only in an HDF5 build. A caller reaches these from inside its own CPP_HDF
!> branch, which is where it has a group to write into at all.
MODULE m_wannierlib_hdf_util
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wr_r4, wr_i3
CONTAINS

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
