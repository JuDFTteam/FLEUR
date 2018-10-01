MODULE m_types_mat
  USE m_judft
  IMPLICIT NONE

  !<This is the basic type to store and manipulate real/complex rank-2 matrices
  !!
  !! In its simple implementation here, it contains a fields for the matrix-size and
  !! a real and complex array for the data
  !! This data-type will be overwritten for distributed matrixes by t_mpimat as defined in types_mpimat.F90
  
   TYPE :: t_mat
     LOGICAL :: l_real                     !>Store either real or complex data
     INTEGER :: matsize1=-1                !> matsize1=size(data_?,1),i.e. no of rows
     INTEGER :: matsize2=-1                !> matsize2=size(data_?,2),i.e. no of columns
     REAL,ALLOCATABLE    :: data_r(:,:)
     COMPLEX,ALLOCATABLE :: data_c(:,:)
   CONTAINS
     PROCEDURE        :: alloc => t_mat_alloc                !> allocate the data-arrays
     PROCEDURE        :: multiply=>t_mat_multiply            !> do a matrix-matrix multiply
     PROCEDURE        :: transpose=>t_mat_transpose          !> transpose the matrix
     PROCEDURE        :: from_packed=>t_mat_from_packed      !> initialized from a packed-storage matrix
     PROCEDURE        :: inverse =>t_mat_inverse             !> invert the matrix
     PROCEDURE        :: to_packed=>t_mat_to_packed          !> convert to packed-storage matrix
     PROCEDURE        :: clear => t_mat_clear                !> set data arrays to zero
     PROCEDURE        :: copy => t_mat_copy                  !> copy into another t_mat (overloaded for t_mpimat)
     PROCEDURE        :: move => t_mat_move                  !> move data into another t_mat (overloaded for t_mpimat)
     PROCEDURE        :: init_details => t_mat_init
     PROCEDURE        :: init_template => t_mat_init_template              !> initalize the matrix(overloaded for t_mpimat)
     GENERIC          :: init => init_details,init_template
     PROCEDURE        :: free => t_mat_free                  !> dealloc the data (overloaded for t_mpimat)
     PROCEDURE        :: add_transpose => t_mat_add_transpose!> add the tranpose/Hermitian conjg. without the diagonal (overloaded for t_mpimat)
   END type t_mat
   
 CONTAINS
   
   SUBROUTINE t_mat_free(mat)
     CLASS(t_mat),INTENT(INOUT)::mat
     IF (ALLOCATED(mat%data_c)) DEALLOCATE(mat%data_c)
     IF (ALLOCATED(mat%data_r)) DEALLOCATE(mat%data_r)
   END SUBROUTINE t_mat_free
   
   SUBROUTINE t_mat_add_transpose(mat,mat1)
    CLASS(t_mat),INTENT(INOUT)::mat,mat1
    INTEGER::i,j
    IF ((mat%matsize1.NE.mat1%matsize2).OR. &
         (mat%matsize2.NE.mat1%matsize1)) &
         CALL judft_error("Matrix sizes missmatch in add_transpose")
    IF (mat%l_real.AND.mat1%l_real) THEN
       DO i=1,mat%matsize2
          DO j=i+1,mat%matsize1
             mat%data_r(j,i)=mat1%data_r(i,j)
          ENDDO
       ENDDO
    ELSEIF((.NOT.mat%l_real).AND.(.NOT.mat1%l_real)) THEN
       DO i=1,mat%matsize2
          DO j=i+1,mat%matsize1
             mat%data_c(j,i)=CONJG(mat1%data_c(i,j))
          ENDDO
       ENDDO
    ELSE
       call judft_error("Inconsistency between data types in m_mat")
    END IF
  END SUBROUTINE t_mat_add_transpose

  
   
   SUBROUTINE t_mat_init(mat,l_real,matsize1,matsize2,mpi_subcom,l_2d,nb_x,nb_y)
     CLASS(t_mat) :: mat
     LOGICAL,INTENT(IN),OPTIONAL:: l_real
     INTEGER,INTENT(IN),OPTIONAL:: matsize1,matsize2
     INTEGER,INTENT(IN),OPTIONAL:: mpi_subcom,nb_x,nb_y !not needed here, only for allowing overloading this in t_mpimat
     LOGICAL,INTENT(IN),OPTIONAL:: l_2d                 !not needed here either

     CALL mat%alloc(l_real,matsize1,matsize2)
   END SUBROUTINE t_mat_init
   SUBROUTINE t_mat_init_template(mat,templ)
     IMPLICIT NONE
     CLASS(t_mat),INTENT(INOUT) :: mat
     CLASS(t_mat),INTENT(IN)    :: templ
     SELECT TYPE(templ)
     TYPE is(t_mat)
        mat%l_real=templ%l_real
        mat%matsize1=templ%matsize1
        mat%matsize2=templ%matsize2
        IF (mat%l_real) THEN
           ALLOCATE(mat%data_r(mat%matsize1,mat%matsize2))
           ALLOCATE(mat%data_c(1,1))
           mat%data_r=0.0
        ELSE
           ALLOCATE(mat%data_c(mat%matsize1,mat%matsize2))
           ALLOCATE(mat%data_r(1,1))
           mat%data_c=0.0
        END IF
     CLASS default
        CALL judft_error("Mixed initialization in t_mat not possible(BUG)")
     END SELECT
   END SUBROUTINE t_mat_init_template
     
  SUBROUTINE t_mat_alloc(mat,l_real,matsize1,matsize2,init)
    CLASS(t_mat) :: mat
    LOGICAL,INTENT(IN),OPTIONAL:: l_real
    INTEGER,INTENT(IN),OPTIONAL:: matsize1,matsize2
    REAL,INTENT(IN),OPTIONAL   :: init

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
       ALLOCATE(mat%data_c(0,0))
       mat%data_r=0.0
       if (present(init)) mat%data_r=init
    ELSE
       ALLOCATE(mat%data_r(0,0))
       ALLOCATE(mat%data_c(mat%matsize1,mat%matsize2),STAT=err)
       mat%data_c=0.0
       IF (PRESENT(init)) mat%data_c=init
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

  SUBROUTINE t_mat_from_packed(mat1,l_real,matsize,packed_r,packed_c)
    CLASS(t_mat),INTENT(INOUT)       :: mat1
    INTEGER,INTENT(IN)               :: matsize
    LOGICAL,INTENT(IN)               :: l_real
    REAL,INTENT(IN)                  :: packed_r(:)
    COMPLEX,INTENT(IN)               :: packed_c(:)

    INTEGER:: n,nn,i
    call mat1%alloc(l_real,matsize,matsize)
    i=1
    DO n=1,matsize
       DO nn=1,n
          if (l_real) THEN
             mat1%data_r(n,nn)=packed_r(i)
             mat1%data_r(nn,n)=packed_r(i)
          else
             mat1%data_c(n,nn)=conjg(packed_c(i))
             mat1%data_c(nn,n)=packed_c(i)
          end if
          i=i+1
       end DO
    end DO
  end SUBROUTINE t_mat_from_packed

  function t_mat_to_packed(mat)result(packed)
    CLASS(t_mat),INTENT(IN)       :: mat
    COMPLEX                       :: packed(mat%matsize1*(mat%matsize1+1)/2)
    integer :: n,nn,i
    real,parameter :: tol=1e-5
    if (mat%matsize1.ne.mat%matsize2) call judft_error("Could not pack no-square matrix",hint='This is a BUG, please report')
    i=1
    DO n=1,mat%matsize1
       DO nn=1,n
          if (mat%l_real) THEN
             packed(i)=(mat%data_r(n,nn)+mat%data_r(nn,n))/2.
             if (abs(mat%data_r(n,nn)-mat%data_r(nn,n))>tol) call judft_warn("Large unsymmetry in matrix packing")
          else
             packed(i)=(conjg(mat%data_c(n,nn))+mat%data_c(nn,n))/2.
             if (abs(conjg(mat%data_c(n,nn))-mat%data_c(nn,n))>tol) call judft_warn("Large unsymmetry in matrix packing")
          endif
          i=i+1
       end DO
    end DO
  end function t_mat_to_packed

  subroutine t_mat_inverse(mat)
    implicit none
    CLASS(t_mat),INTENT(INOUT)       :: mat
    integer                :: info
    real, allocatable      :: work_r(:)
    integer, allocatable   :: ipiv(:)
    complex,allocatable    :: work_c(:)
    
    
    if (mat%matsize1.ne.mat%matsize2) call judft_error("Can only invert square matrices",hint="This is a BUG in FLEUR, please report")
    ALLOCATE(ipiv(mat%matsize1))

    if (mat%l_real) THEN
       ALLOCATE(work_r(mat%matsize1))
       call dgetrf(mat%matsize1,mat%matsize1,mat%data_r,size(mat%data_r,1),ipiv,info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
       call dgetri(mat%matsize1,mat%data_r,size(mat%data_r,1),ipiv,work_r,size(work_r),info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
    else
       ALLOCATE(work_c(mat%matsize1))
       call zgetrf(mat%matsize1,mat%matsize1,mat%data_c,size(mat%data_c,1),ipiv,info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
       call zgetri(mat%matsize1,mat%data_c,size(mat%data_c,1),ipiv,work_c,size(work_c),info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
    end if
  end subroutine t_mat_inverse

  SUBROUTINE t_mat_move(mat,mat1)
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT):: mat
    CLASS(t_mat),INTENT(INOUT):: mat1
    !Special case, the full matrix is copied. Then use move alloc
    IF (mat%l_real) THEN
       CALL move_ALLOC(mat1%data_r,mat%data_r)
    ELSE
       CALL move_ALLOC(mat1%data_c,mat%data_c)
    END IF
  END SUBROUTINE t_mat_move
  
  SUBROUTINE t_mat_copy(mat,mat1,n1,n2)
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT):: mat
    CLASS(t_mat),INTENT(IN)   :: mat1
    INTEGER,INTENT(IN)        :: n1,n2

    INTEGER:: i1,i2

    i1=mat%matsize1-n1+1  !space available for first dimension
    i2=mat%matsize2-n1+1
    i1=MIN(i1,mat1%matsize1)
    i2=MIN(i2,mat1%matsize2)
    IF (mat%l_real) THEN
       mat%data_r(n1:n1+i1-1,n2:n2+i2-1)=mat1%data_r(:i1,:i2)
    ELSE
       mat%data_c(n1:n1+i1-1,n2:n2+i2-1)=mat1%data_c(:i1,:i2)
    END IF
       
  END SUBROUTINE t_mat_copy
 
  SUBROUTINE t_mat_clear(mat)
    IMPLICIT NONE
    CLASS(t_mat),INTENT(INOUT):: mat

    IF (mat%l_real) THEN
       mat%data_r=0.0
    ELSE
       mat%data_c=0.0
    ENDIF
  END SUBROUTINE t_mat_clear
END MODULE m_types_mat

MODULE m_types_rcmat
  USE m_types_mat
END MODULE m_types_rcmat
