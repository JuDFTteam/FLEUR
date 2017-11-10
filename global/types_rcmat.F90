MODULE m_types_rcmat
  use m_judft
    IMPLICIT NONE
    TYPE :: t_mat
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
     PROCEDURE        :: inverse =>t_mat_inverse
     PROCEDURE        :: to_packed=>t_mat_to_packed
   END type t_mat

   TYPE,EXTENDS(t_mat) :: t_lapwmat
      INTEGER,ALLOCATABLE :: local_rk_map(:,:) !<Maps a global k+g value to a local one in MPI case
      INTEGER             :: matsize_half      !<Half of matrix in noco case
    CONTAINS
      PROCEDURE,PASS :: init => t_lapwmat_init
   END TYPE t_lapwmat
   
 CONTAINS
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
       if (present(init)) mat%data_r=init
    ELSE
       ALLOCATE(mat%data_r(0,0))
       ALLOCATE(mat%data_c(mat%matsize1,mat%matsize2),STAT=err)
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


  SUBROUTINE   t_lapwmat_init(mat,lapw,mpi,nlotot,l_noco,jspin,l_real)
    USE m_types_lapw
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_lapwmat),INTENT(OUT):: mat
    TYPE(t_lapw),INTENT(in)     :: lapw
    TYPE(t_mpi),INTENT(IN)      :: mpi
    LOGICAL,INTENT(IN)          :: l_noco
    INTEGER,INTENT(IN)          :: nlotot
    INTEGER,INTENT(IN)          :: jspin
    LOGICAL,INTENT(IN)          :: l_real

    INTEGER :: matsize,i,n
    !Fist calculate size of matrix
    IF (l_noco) THEN
       matsize=SUM(lapw%nv)+2*nlotot
       mat%matsize_half=lapw%nv(1)+nlotot
    ELSE
       matsize=lapw%nv(jspin)+nlotot
       mat%matsize_half=0
    ENDIF
    !Calculate Matrix distribution
    IF (ALLOCATED(mat%local_rk_map)) DEALLOCATE(mat%local_rk_map)
    ALLOCATE(mat%local_rk_map(MAXVAL(lapw%nv)+nlotot,2))
    mat%local_rk_map=0
    i=0
    DO n=mpi%n_rank+1, lapw%nv(jspin)+nlotot, mpi%n_size
       i=i+1
       mat%local_rk_map(n,1)=i
    ENDDO
    IF (l_noco) THEN ! In noco case do it again
       DO n=mpi%n_rank+1, lapw%nv(2)+nlotot, mpi%n_size
          i=i+1
          mat%local_rk_map(n,2)=i
       ENDDO
    ELSE
       mat%local_rk_map(:,2)=mat%local_rk_map(:,1)
    ENDIF
    
    !Allocate space
    CALL mat%alloc(l_real,matsize,i)
  END SUBROUTINE t_lapwmat_init
END MODULE m_types_rcmat
