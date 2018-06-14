!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_io_matrix
  USE m_types
  IMPLICIT NONE
  private
  INTEGER:: mode=1
  public:: open_matrix,read_matrix,write_matrix
CONTAINS
  INTEGER FUNCTION OPEN_matrix(l_real,matsize,no_rec,filename)
    LOGICAL,INTENT(IN)          :: l_real
    INTEGER,INTENT(in)          :: matsize,no_rec
    CHARACTER(len=*),INTENT(in) :: filename
    SELECT CASE (mode)
    CASE (1)
       OPEN_matrix=open_DA(l_real,matsize,no_rec,filename)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
  END FUNCTION OPEN_MATRIX

  SUBROUTINE read_matrix(mat,rec,id)
    TYPE(t_Mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec,id

    CALL mat%alloc()
    SELECT CASE (mode)
    CASE (1)
       CALL read_matrix_DA(mat,rec,id)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
  END SUBROUTINE read_matrix

  SUBROUTINE write_matrix(mat,rec,id)
    TYPE(t_Mat),INTENT(IN):: mat
    INTEGER,INTENT(IN)        :: rec,id

    SELECT CASE (mode)
    CASE (1)
       CALL write_matrix_DA(mat,rec,id)
    CASE default
       CALL judft_error("BUG in io_matrix")
    END SELECT
  END SUBROUTINE write_matrix

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
    OPEN(unit=open_DA,file=filename,access='direct',recl=datasize*(matsize*matsize+4))!Three to include matsize


  END FUNCTION open_DA

  SUBROUTINE read_matrix_DA(mat,rec,id)
    TYPE(t_Mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec,id
    LOGICAL :: l_real
    INTEGER:: err,matsize1,matsize2
    l_real=mat%l_real
    READ(id,rec=rec,iostat=err) matsize1,matsize2
    if (matsize1<1) call judft_error("Data not found in file")
    IF (mat%matsize1.NE.matsize1.OR.mat%matsize2.NE.matsize2) CALL mat%alloc(l_real,matsize1,matsize2)
    IF (mat%l_real) THEN
       READ(id,rec=rec,iostat=err) matsize1,matsize2,mat%data_r
    ELSE
       READ(id,rec=rec,iostat=err) matsize1,matsize2,mat%data_c
    END IF
    IF (err.NE.0) CALL judft_error("Failed in reading of matrix")
  END SUBROUTINE read_matrix_DA

  SUBROUTINE write_matrix_DA(mat,rec,id)
    TYPE(t_Mat),INTENT(IN):: mat
    INTEGER,INTENT(IN)        :: rec,id
    INTEGER:: err
    IF (mat%l_real) THEN
       WRITE(id,rec=rec,iostat=err) mat%matsize1,mat%matsize2,mat%data_r
    ELSE
       WRITE(id,rec=rec,iostat=err) mat%matsize1,mat%matsize2,mat%data_c
    END IF
    IF (err.NE.0) CALL judft_error("Failed in writing of matrix")
  END SUBROUTINE write_matrix_DA

END MODULE m_io_matrix
