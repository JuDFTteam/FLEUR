!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_io_matrix
  USE m_types_mat
  USE m_types_mpimat
  USE m_judft
  USE m_iomatrix_hdf
#ifdef CPP_HDF
  USE hdf5
#endif
  IMPLICIT NONE
  PRIVATE
  TYPE t_iomatrix_handle
     INTEGER:: mode=0 !can be 1 for DA or 2 for HDF
     INTEGER:: id !file ID in direct-access mode
     INTEGER(hid_t):: fid,did !file-handle in hdf mode
  END TYPE t_iomatrix_handle

  TYPE(t_iomatrix_handle)::fh(10)

  PUBLIC:: t_iomatrix_handle,open_matrix,read_matrix,write_matrix,close_matrix
CONTAINS
  INTEGER FUNCTION OPEN_matrix(l_real,matsize,mode,no_rec,filename)
    LOGICAL,INTENT(IN)          :: l_real
    INTEGER,INTENT(in)          :: matsize,no_rec,mode
    CHARACTER(len=*),INTENT(in) :: filename
    !Find free handle
    DO open_matrix=1,SIZE(fh)
       IF (fh(open_matrix)%mode==0) EXIT
    ENDDO
    IF (open_matrix>SIZE(fh)) CALL judft_error("Too many filehandles for matrix IO")

    SELECT CASE (mode)
    CASE (1)
       fh(open_matrix)%mode=1
       fh(OPEN_matrix)%id=open_DA(l_real,matsize,no_rec,filename)
    CASE(2)
       fh(open_matrix)%mode=2
       CALL iomatrix_hdf_open(l_real,matsize,no_rec,filename,fh(open_matrix)%fid,fh(open_matrix)%did)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
  END FUNCTION OPEN_MATRIX

  SUBROUTINE read_matrix(mat,rec,id)
    CLASS(t_Mat),INTENT(INOUT)  :: mat
    INTEGER,INTENT(IN)          :: rec,id

    CALL mat%alloc()
    SELECT CASE (fh(id)%mode)
    CASE (1)
       SELECT TYPE(mat)
       TYPE is (t_mat)
          CALL read_matrix_DA(mat,rec,fh(id)%id)
       TYPE is (t_mpimat)
          CALL judft_error("Matrix IO for parallel matrix only with HDF5")
       END SELECT  
    CASE(2)
       CALL iomatrix_hdf_read(mat,rec,fh(id)%did)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
  END SUBROUTINE read_matrix

  SUBROUTINE write_matrix(mat,rec,id)
    CLASS(t_Mat),INTENT(IN)  :: mat
    INTEGER,INTENT(IN)       :: rec,id

    SELECT CASE (fh(id)%mode)
    CASE (1)
       SELECT TYPE(mat)
       TYPE is (t_mat)
          CALL write_matrix_DA(mat,rec,fh(id)%id)
       TYPE is (t_mpimat)
          CALL judft_error("Matrix IO for parallel matrix only with HDF5")
       END SELECT
    CASE(2)
       CALL iomatrix_hdf_write(mat,rec,fh(id)%did)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
  END SUBROUTINE write_matrix
  
  SUBROUTINE close_matrix(id)
    INTEGER,INTENT(IN):: id
    SELECT CASE (fh(id)%mode)
    CASE (1)
       CALL close_matrix_DA(fh(id)%id)
    CASE (2)
       CALL iomatrix_hdf_close(fh(id)%fid,fh(id)%did)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
       fh(id)%mode=0
  END SUBROUTINE CLOSE_MATRIX

  !Now the implementation in terms of fortran DA-files
  INTEGER FUNCTION open_DA(l_real,matsize,no_rec,filename)
    LOGICAL,INTENT(IN)           :: l_real
    INTEGER,INTENT(in)          :: matsize,no_rec
    CHARACTER(len=*),INTENT(in) :: filename

    LOGICAL :: used_unit
    REAL    :: r
    COMPLEX :: c
    INTEGER :: datasize

    !Determine size of data
    IF (l_real) THEN
       INQUIRE(IOLENGTH=datasize) r
    ELSE
       INQUIRE(IOLENGTH=datasize) c
    END IF

    !find free unit starting at 901
    open_DA=901
    DO
       INQUIRE(unit=open_DA,opened=used_unit)
       IF (.NOT.used_unit) EXIT
       open_DA=open_DA+1
    END DO
    !openfile
    OPEN(unit=open_DA,file=filename,access='direct',recl=datasize*(matsize*matsize+6))!Three to include matsize


  END FUNCTION open_DA

  SUBROUTINE read_matrix_DA(mat,rec,id)
    TYPE(t_Mat),INTENT(OUT):: mat
    INTEGER,INTENT(IN)           :: rec,id
    LOGICAL :: l_real
    INTEGER:: err,matsize1,matsize2
    l_real=mat%l_real
   
    READ(id,rec=rec,iostat=err) l_real,matsize1,matsize2
    IF (err.NE.0) CALL judft_error("Data not found in file")
    CALL mat%init(l_real,matsize1,matsize2)

    IF (mat%l_real) THEN
       READ(id,rec=rec,iostat=err) l_real,matsize1,matsize2,mat%data_r
    ELSE
       READ(id,rec=rec,iostat=err) l_real,matsize1,matsize2,mat%data_c
    END IF
    IF (err.NE.0) CALL judft_error("Failed in reading of matrix")
  END SUBROUTINE read_matrix_DA

  SUBROUTINE write_matrix_DA(mat,rec,id)
    TYPE(t_Mat),INTENT(IN):: mat
    INTEGER,INTENT(IN)        :: rec,id
    INTEGER:: err
    IF (mat%l_real) THEN
       WRITE(id,rec=rec,iostat=err) mat%l_real,mat%matsize1,mat%matsize2,mat%data_r
    ELSE
       WRITE(id,rec=rec,iostat=err) mat%l_real,mat%matsize1,mat%matsize2,mat%data_c
    END IF
    IF (err.NE.0) CALL judft_error("Failed in writing of matrix")
  END SUBROUTINE write_matrix_DA

  SUBROUTINE close_matrix_DA(id)
    INTEGER,INTENT(IN)        :: id
    INTEGER:: err

    close(id)
  END SUBROUTINE close_matrix_DA

END MODULE m_io_matrix
