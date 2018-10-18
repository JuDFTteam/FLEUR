MODULE m_types_gpumat
  USE m_judft
  USE m_types_mat
#ifdef CPP_GPU    
  USE cudafor
#endif  
  IMPLICIT NONE

  !<Some routines are overwritten for GPU handling
  !!
  !! Several select type constructs are only needed to make the PGI compiler happy...
  TYPE,EXTENDS(t_mat) :: t_gpumat
#ifdef CPP_GPU     
   CONTAINS
     PROCEDURE        :: alloc =>t_gpumat_alloc
     PROCEDURE        :: multiply=>t_gpumat_multiply            !> do a matrix-matrix multiply
     PROCEDURE        :: transpose=>t_gpumat_transpose          !> transpose the matrix
     PROCEDURE        :: from_packed=>t_gpumat_from_packed      !> initialized from a packed-storage matrix
     PROCEDURE        :: inverse =>t_gpumat_inverse             !> invert the matrix
     PROCEDURE        :: to_packed=>t_gpumat_to_packed          !> convert to packed-storage matrix
     PROCEDURE        :: clear => t_gpumat_clear                !> set data arrays to zero
     PROCEDURE        :: copy => t_gpumat_copy                  !> copy into another t_mat 
     PROCEDURE        :: move => t_gpumat_move                  !> move data into another t_mat 
     PROCEDURE        :: add_transpose => t_gpumat_add_transpose!> add the tranpose/Hermitian conjg. without the diagonal 
     PROCEDURE        :: init_from                              !> Initialize from an existing matrix
#endif
  END TYPE t_gpumat

CONTAINS
#ifdef CPP_GPU
  SUBROUTINE t_gpumat_alloc(mat,l_real,matsize1,matsize2,init)
    CLASS(t_gpumat) :: mat
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
    ELSE
       ALLOCATE(mat%data_r(0,0))
       ALLOCATE(mat%data_c(mat%matsize1,mat%matsize2),STAT=err)
    ENDIF

    IF (err>0) CALL judft_error("Allocation of memmory failed for mat datatype",hint="You probably run out of memory")
  END SUBROUTINE t_gpumat_alloc

  SUBROUTINE init_from(this,mat)
    CLASS(t_gpumat),INTENT(INOUT):: this
    CLASS(t_mat),INTENT(IN)      :: mat

    SELECT TYPE(mat)
    TYPE IS (t_mat)
       call this%alloc(mat%l_real,mat%matsize1,mat%matsize2)
       IF (mat%l_real) THEN
          this%data_r=mat%data_r
       ELSE
          this%data_c=mat%data_c
       END IF
       CLASS default
       CALL judft_error("Initializiation not implemented in t_gpumat")
    END SELECT
  END SUBROUTINE init_from


  SUBROUTINE t_gpumat_add_transpose(mat,mat1)
    USE m_types
    IMPLICIT NONE
    CLASS(t_gpumat),INTENT(INOUT):: mat
    CLASS(t_mat),INTENT(IN)    :: mat1
    INTEGER::i,j
    SELECT TYPE(mat)
    TYPE is (t_gpumat)
       SELECT TYPE(mat1)
       TYPE is (t_gpumat)
          IF ((mat%matsize1.NE.mat1%matsize2).OR. &
               (mat%matsize2.NE.mat1%matsize1)) &
               CALL judft_error("Matrix sizes missmatch in add_transpose")
          IF (mat%l_real.AND.mat1%l_real) THEN
             !$cuf kernel do <<<*,256>>>
             DO i=1,mat%matsize2
                DO j=i+1,mat%matsize1
                   mat%data_r(j,i)=mat1%data_r(i,j)
                ENDDO
             ENDDO
          ELSEIF((.NOT.mat%l_real).AND.(.NOT.mat1%l_real)) THEN
             !$cuf kernel do <<<*,256>>>
             DO i=1,mat%matsize2
                DO j=i+1,mat%matsize1
                   mat%data_c(j,i)=CONJG(mat1%data_c(i,j))
                ENDDO
             ENDDO
          ELSE
             CALL judft_error("Inconsistency between data types in m_mat")
          END IF
       END SELECT
       i=cudaDeviceSynchronize()
    END SELECT
  END SUBROUTINE t_gpumat_add_transpose



  SUBROUTINE t_gpumat_multiply(mat1,mat2,res)
    CLASS(t_gpumat),INTENT(INOUT)       ::mat1
    CLASS(t_mat),INTENT(IN)           ::mat2
    CLASS(t_mat),INTENT(OUT),OPTIONAL ::res

    if (mat1%matsize2.ne.mat2%matsize1) CALL judft_error("Cannot multiply matrices because of non-matching dimensions",hint="This is a BUG in FLEUR, please report")
    CALL judft_error(" multiply method not implemented for t_gpumat")

  end SUBROUTINE t_gpumat_multiply


  SUBROUTINE t_gpumat_transpose(mat1,res)
    CLASS(t_gpumat),INTENT(INOUT)       ::mat1
    TYPE(t_gpumat),INTENT(OUT),OPTIONAL ::res

    INTEGER :: i,j

    IF (PRESENT(res)) THEN
       CALL res%alloc(mat1%l_real,mat1%matsize2,mat1%matsize1)
       IF (mat1%l_real) THEN
          DO i=1,mat1%matsize1
             DO j=1,mat1%matsize2
                res%data_r(i,j)=mat1%data_r(j,i)
             ENDDO
          ENDDO
       ELSE
          DO i=1,mat1%matsize1
             DO j=1,mat1%matsize2
                res%data_c(i,j)=mat1%data_c(j,i)
             ENDDO
          ENDDO
       ENDIF
    ELSE
       SELECT TYPE(mat1)
       TYPE is (t_gpumat)
          IF (mat1%matsize1.NE.mat1%matsize2) CALL judft_error("Cannot transpose matrices inplace because of non-matching dimensions",hint="This is a BUG in FLEUR, please report")
          IF (mat1%l_real) THEN
             !$cuf kernel do <<<*,256>>>
             DO i=1,mat1%matsize1
                DO j=1,mat1%matsize2
                   mat1%data_r(i,j)=mat1%data_r(j,i)
                ENDDO
             ENDDO
          ELSE
             !$cuf kernel do <<<*,256>>>
             DO i=1,mat1%matsize1
                DO j=1,mat1%matsize2
                   mat1%data_c(i,j)=mat1%data_c(j,i)
                ENDDO
             ENDDO
          ENDIF
       END SELECT

       i=cudaDeviceSynchronize()
    END IF
  end SUBROUTINE t_gpumat_transpose

  SUBROUTINE t_gpumat_from_packed(mat1,l_real,matsize,packed_r,packed_c)
    CLASS(t_gpumat),INTENT(INOUT)       :: mat1
    INTEGER,INTENT(IN)               :: matsize
    LOGICAL,INTENT(IN)               :: l_real
    REAL,INTENT(IN)                  :: packed_r(:)
    COMPLEX,INTENT(IN)               :: packed_c(:)
    CALL judft_error(" mat_from_packed method not implemented for t_gpumat")
  end SUBROUTINE t_gpumat_from_packed

  function t_gpumat_to_packed(mat)result(packed)
    CLASS(t_gpumat),INTENT(IN)       :: mat
    COMPLEX                       :: packed(mat%matsize1*(mat%matsize1+1)/2)
    integer :: n,nn,i
    real,parameter :: tol=1e-5
    if (mat%matsize1.ne.mat%matsize2) call judft_error("Could not pack no-square matrix",hint='This is a BUG, please report')
    CALL judft_error(" mat_to_packed method not implemented for t_gpumat")
  end function t_gpumat_to_packed

  SUBROUTINE t_gpumat_inverse(mat)
    USE cublas
    implicit none
    CLASS(t_gpumat),INTENT(INOUT)       :: mat
    integer                :: info
    REAL, ALLOCATABLE,device      :: work_r(:)
    INTEGER, ALLOCATABLE,device   :: ipiv(:)
    COMPLEX,ALLOCATABLE,device    :: work_c(:)


    if (mat%matsize1.ne.mat%matsize2) call judft_error("Can only invert square matrices",hint="This is a BUG in FLEUR, please report")
    ALLOCATE(ipiv(mat%matsize1))

    if (mat%l_real) THEN
       ALLOCATE(work_r(mat%matsize1))
       call dgetrf(mat%matsize1,mat%matsize1,mat%data_r,mat%matsize1,ipiv,info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
       call dgetri(mat%matsize1,mat%data_r,mat%matsize1,ipiv,work_r,mat%matsize1,info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
    else
       ALLOCATE(work_c(mat%matsize1))
       call zgetrf(mat%matsize1,mat%matsize1,mat%data_c,mat%matsize1,ipiv,info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
       call zgetri(mat%matsize1,mat%data_c,size(mat%data_c,1),ipiv,work_c,mat%matsize1,info)
       if(info.ne.0) call judft_error("Failed to invert matrix: dpotrf failed.")
    end if
  end subroutine t_gpumat_inverse

  SUBROUTINE t_gpumat_move(mat,mat1)
    IMPLICIT NONE
    CLASS(t_gpumat),INTENT(INOUT):: mat
    CLASS(t_gpumat),INTENT(INOUT):: mat1
    !Special case, the full matrix is copied. Then use move alloc
    CALL mat%copy(mat1,1,1)
  END SUBROUTINE t_gpumat_move

  SUBROUTINE t_gpumat_copy(mat,mat1,n1,n2)
    IMPLICIT NONE
    CLASS(t_gpumat),INTENT(INOUT):: mat
    CLASS(t_mat),INTENT(IN)   :: mat1
    INTEGER,INTENT(IN)        :: n1,n2

    INTEGER:: i1,i2,i,j
    LOGICAL no_gpu

    no_gpu=.FALSE.
    i1=mat%matsize1-n1+1  !space available for first dimension
    i2=mat%matsize2-n1+1
    i1=MIN(i1,mat1%matsize1)
    i2=MIN(i2,mat1%matsize2)

    SELECT TYPE(mat)
    TYPE IS (t_gpumat)
       SELECT TYPE(mat1)
       TYPE IS (t_gpumat)
          IF (mat%l_real) THEN
             !$cuf kernel do <<<*,256>>>
             DO i=1,i1
                DO j=1,i2
                   mat%data_r(n1+i-1,n2+j-1)=mat1%data_r(i,j)
                ENDDO
             ENDDO
          ELSE
             !$cuf kernel do <<<*,256>>>
             DO i=1,i1
                DO j=1,i2
                   mat%data_c(n1+i-1,n2+j-1)=mat1%data_c(i,j)
                ENDDO
             ENDDO
          END IF
          i=cudaDeviceSynchronize()
          CLASS default
          no_gpu=.TRUE.
       END SELECT
       CLASS default
       no_gpu=.TRUE.
    END SELECT
    IF (no_gpu)THEN
       IF (mat%l_real) THEN
          DO i=1,i1
             DO j=1,i2
                mat%data_r(n1+i-1,n2+j-1)=mat1%data_r(i,j)
             ENDDO
          ENDDO
       ELSE
          DO i=1,i1
             DO j=1,i2
                mat%data_c(n1+i-1,n2+j-1)=mat1%data_c(i,j)
             ENDDO
          ENDDO
       END IF
    END IF
  END SUBROUTINE t_gpumat_copy

  SUBROUTINE t_gpumat_clear(mat)
    IMPLICIT NONE
    CLASS(t_gpumat),INTENT(INOUT):: mat
    INTEGER:: n,nn
    SELECT TYPE(mat)
    TYPE IS (t_gpumaT)
       IF (mat%l_real) THEN
          !$cuf kernel do <<<*,256>>>
          DO n=1,mat%matsize1
             DO nn=1,mat%matsize2
                mat%data_r(n,nn)=0.0
             ENDDO
          ENDDO
       ELSE
          !$cuf kernel do <<<*,256>>>
          DO n=1,mat%matsize1
             DO nn=1,mat%matsize2
                mat%data_c(n,nn)=0.0
             ENDDO
          ENDDO
       ENDIF
       n=cudaDeviceSynchronize()
    END SELECT
  END SUBROUTINE t_gpumat_clear
#endif  
END MODULE m_types_gpumat
