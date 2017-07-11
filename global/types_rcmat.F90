module m_types_rcmat
  use m_judft
    IMPLICIT NONE
    TYPE t_mat
     LOGICAL :: l_real
     INTEGER :: matsize1=-1
     INTEGER :: matsize2=-1
     REAL,ALLOCATABLE    :: data_r(:,:)
     COMPLEX,ALLOCATABLE :: data_c(:,:)
   CONTAINS
     PROCEDURE        :: alloc => t_mat_alloc
     PROCEDURE        :: multiply=>t_mat_multiply
     PROCEDURE        :: transpose=>t_mat_transpose
     PROCEDURE        :: from_packed=>t_mat_from_packed
   END type t_mat

 CONTAINS
  SUBROUTINE t_mat_alloc(mat,l_real,matsize1,matsize2)
    CLASS(t_mat) :: mat
    LOGICAL,INTENT(IN),OPTIONAL:: l_real
    INTEGER,INTENT(IN),OPTIONAL:: matsize1,matsize2

    INTEGER:: err
    IF (present(l_real)) mat%l_real=l_real
    IF (present(matsize1)) mat%matsize1=matsize1
    IF (present(matsize2)) mat%matsize2=matsize2

    IF (mat%matsize1<0) CALL judft_error("Cannot allocate memory for mat datatype that is not initialized",hint="This is a BUG in FLEUR, please report")
    IF (mat%matsize2<0)  mat%matsize2=mat%matsize1 !Square by default
    
    IF (allocated(mat%data_r)) deallocate(mat%data_r)
    IF (allocated(mat%data_c)) deallocate(mat%data_c)
    
    IF (mat%l_real) THEN
       ALLOCATE(mat%data_r(mat%matsize1,mat%matsize2),STAT=err)
    ELSE
       ALLOCATE(mat%data_c(mat%matsize1,mat%matsize2),STAT=err)
    ENDIF

    IF (err>0) CALL judft_error("Allocation of memmory failed for mat datatype",hint="You probably run out of memory")
  END SUBROUTINE t_mat_alloc

  SUBROUTINE t_mat_multiply(mat1,mat2,res)
    CLASS(t_mat),INTENT(INOUT)       ::mat1
    TYPE(t_mat),INTENT(IN)           ::mat2
    TYPE(t_mat),INTENT(OUT),OPTIONAL ::res

    if (mat1%matsize2.ne.mat2%matsize1) CALL judft_error("Cannot multiply matrices because of non-matching dimensions",hint="This is a BUG in FLEUR, please report")
    
    IF (present(res)) THEN
       call res%alloc(mat1%l_real,mat1%matsize1,mat2%matsize2)
       IF (mat1%l_real) THEN
          res%data_r=matmul(mat1%data_r(:mat1%matsize1,:mat1%matsize2),mat2%data_r(:mat2%matsize1,:mat2%matsize2))
       ELSE
          res%data_c=matmul(mat1%data_c(:mat1%matsize1,:mat1%matsize2),mat2%data_c(:mat2%matsize1,:mat2%matsize2))
       ENDIF
    else
       if (mat1%matsize1.ne.mat1%matsize2) CALL judft_error("Cannot multiply matrices inplace because of non-matching dimensions",hint="This is a BUG in FLEUR, please report")
       if (mat1%l_real) THEN
          mat1%data_r(:mat1%matsize1,:mat1%matsize2)=matmul(mat1%data_r(:mat1%matsize1,:mat1%matsize2),mat2%data_r(:mat2%matsize1,:mat2%matsize2))
       ELSE
          mat1%data_c(:mat1%matsize1,:mat1%matsize2)=matmul(mat1%data_c(:mat1%matsize1,:mat1%matsize2),mat2%data_c(:mat2%matsize1,:mat2%matsize2))
       ENDIF
    end IF
  end SUBROUTINE t_mat_multiply


  SUBROUTINE t_mat_transpose(mat1,res)
    CLASS(t_mat),INTENT(INOUT)       ::mat1
    TYPE(t_mat),INTENT(OUT),OPTIONAL ::res
    
    IF (present(res)) THEN
       call res%alloc(mat1%l_real,mat1%matsize2,mat1%matsize1)
       IF (mat1%l_real) THEN
          res%data_r=transpose(mat1%data_r(:mat1%matsize1,:mat1%matsize2))
       ELSE
          res%data_c=transpose(mat1%data_c(:mat1%matsize1,:mat1%matsize2))
       ENDIF
    else
       if (mat1%matsize1.ne.mat1%matsize2) CALL judft_error("Cannot transpose matrices inplace because of non-matching dimensions",hint="This is a BUG in FLEUR, please report")
       IF (mat1%l_real) THEN
          mat1%data_r(:mat1%matsize1,:mat1%matsize2)=transpose(mat1%data_r(:mat1%matsize1,:mat1%matsize2))
       ELSE
          mat1%data_c(:mat1%matsize1,:mat1%matsize2)=transpose(mat1%data_c(:mat1%matsize1,:mat1%matsize2))
       ENDIF
    end IF
  end SUBROUTINE t_mat_transpose

  SUBROUTINE t_mat_from_packed(mat1,matsize,l_real,packed_r,packed_c)
    CLASS(t_mat),INTENT(INOUT)       :: mat1
    INTEGER,INTENT(IN)               :: matsize
    LOGICAL,INTENT(IN)               :: l_real
    REAL,INTENT(IN)                  :: packed_r(:)
    COMPLEX,INTENT(IN)               :: packed_c(:)

    INTEGER:: n,nn,i
    call mat1%alloc(l_real,matsize)
    i=1
    DO n=1,matsize
       DO nn=1,n
          if (l_real) THEN
             mat1%data_r(n,nn)=packed_r(i)
          else
             mat1%data_c(n,nn)=packed_c(i)
          end if
          i=i+1
       end DO
    end DO
  end SUBROUTINE t_mat_from_packed
end module m_types_rcmat
