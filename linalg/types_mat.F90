MODULE m_types_mat
#ifdef _OPENACC
   use openacc
   use cusolverDn
   use cublas
#endif
   USE m_judft
   use m_constants
   IMPLICIT NONE
   PRIVATE
   !<This is the basic type to store and manipulate real/complex rank-2 matrices
   !!
   !! In its simple implementation here, it contains a fields for the matrix-size and
   !! a real and complex array for the data
   !! This data-type will be overwritten for distributed matrixes by t_mpimat as defined in types_mpimat.F90

   TYPE :: t_mat
      LOGICAL :: l_real                       !>Store either real or complex data
      INTEGER :: matsize1 = -1                !> matsize1=size(data_?,1),i.e. no of rows
      INTEGER :: matsize2 = -1                !> matsize2=size(data_?,2),i.e. no of columns
      LOGICAL :: hidden_transpose=.false.
      REAL    :: hidden_scaled=1.0
      LOGICAL :: data_on_device=.false.
      REAL, ALLOCATABLE    :: data_r(:, :)
      COMPLEX, ALLOCATABLE :: data_c(:, :)
   CONTAINS
      !Initialization
      PROCEDURE        :: t_mat_init
      PROCEDURE        :: t_mat_init_template              !> initalize the matrix(overloaded for t_mpimat)
      PROCEDURE        :: t_mat_init_rmatrix,t_mat_init_cmatrix
      GENERIC          :: init => init_details, init_template,t_mat_init_rmatrix,t_mat_init_cmatrix,t_mat_init_cmatrix_alloc,t_mat_init_rmatrix_alloc
      !Data movement and initialization
      PROCEDURE        :: clear => t_mat_clear                !> set data arrays to zero
      PROCEDURE        :: copy => t_mat_copy                  !> copy into another t_mat (overloaded for t_mpimat)
      PROCEDURE        :: move => t_mat_move                  !> move data into another t_mat (overloaded for t_mpimat)
      PROCEDURE        :: free => t_mat_free                  !> dealloc the data (overloaded for t_mpimat)
      !Calculations
      PROCEDURE        :: multiply => t_mat_multiply            !> do a matrix-matrix multiply
      procedure        :: add=>t_mat_add
      PROCEDURE        :: transpose => t_mat_transpose          !> transpose the matrix
      PROCEDURE        :: evaluate =>t_mat_evaluate
      PROCEDURE        :: scale =>t_mat_scale
      !GPU stuff
      PROCEDURE        :: tohost => t_mat_data2host
      PROCEDURE        :: todevice => t_mat_data2device
      
   END type t_mat

   interface assignment(=)
       module procedure t_mat_init_matrix
   end interface

   interface operator (+)
       module procedure t_mat_multiply
   end interface
   
   PUBLIC t_mat
CONTAINS
SUBROUTINE t_mat_init(mat, l_real, matsize1, matsize2, mpi_subcom, l_2d, nb_x, nb_y, mat_name)
   CLASS(t_mat) :: mat
   LOGICAL, INTENT(IN), OPTIONAL        :: l_real
   INTEGER, INTENT(IN), OPTIONAL        :: matsize1, matsize2
   INTEGER, INTENT(IN), OPTIONAL        :: mpi_subcom, nb_x, nb_y !not needed here, only for allowing overloading this in t_mpimat
   LOGICAL, INTENT(IN), OPTIONAL        :: l_2d                 !not needed here either
   character(len=*),intent(in),optional :: mat_name

   CALL mat%alloc(l_real, matsize1, matsize2, mat_name=mat_name)
END SUBROUTINE t_mat_init

SUBROUTINE t_mat_init_template(mat, templ, global_size1, global_size2, mat_name)
   IMPLICIT NONE
   CLASS(t_mat), INTENT(INOUT) :: mat
   CLASS(t_mat), INTENT(IN)    :: templ
   INTEGER, INTENT(IN), OPTIONAL:: global_size1, global_size2
   character(len=*),intent(in),optional :: mat_name

   integer :: ierr

   IF (PRESENT(global_size1) .AND. PRESENT(global_size2)) THEN
      IF ((global_size1 .NE. templ%matsize1) .OR. (global_size2 .NE. templ%matsize2)) CALL judft_error("BUG:Invalid change of size in init by template")
   END IF
   mat%l_real = templ%l_real
   mat%matsize1 = templ%matsize1
   mat%matsize2 = templ%matsize2
   IF (mat%l_real) THEN
      ALLOCATE (mat%data_r(mat%matsize1, mat%matsize2), source=0.0, stat=ierr)
      ALLOCATE (mat%data_c(1, 1))
   ELSE
      ALLOCATE (mat%data_c(mat%matsize1, mat%matsize2), source=(0.0,0.0), stat=ierr)
      ALLOCATE (mat%data_r(1, 1))
   END IF
   if(ierr /= 0) then
      if(present(mat_name)) then
         call judft_error("can't alloc matrix of size: [" // &
            int2str(mat%matsize1) // ", " // int2str(mat%matsize2) // "]. Name: " // trim(mat_name))
      else
         call judft_error("can't alloc matrix of size: [" // &
            int2str(mat%matsize1) // ", " // int2str(mat%matsize2) // "].")
      endif
   endif
END SUBROUTINE t_mat_init_template

SUBROUTINE t_mat_alloc(mat, l_real, matsize1, matsize2, init, mat_name)
   use m_judft
   CLASS(t_mat) :: mat
   LOGICAL, INTENT(IN), OPTIONAL:: l_real
   INTEGER, INTENT(IN), OPTIONAL:: matsize1, matsize2
   REAL, INTENT(IN), OPTIONAL   :: init
   character(len=*), intent(in), optional :: mat_name
   character(len=300)           :: errmsg

   INTEGER:: err

   call timestart("t_mat_alloc")
   IF (present(l_real)) mat%l_real = l_real
   IF (present(matsize1)) mat%matsize1 = matsize1
   IF (present(matsize2)) mat%matsize2 = matsize2

   IF (mat%matsize1 < 0) CALL judft_error("Cannot allocate memory for mat datatype that is not initialized", hint="This is a BUG in FLEUR, please report")
   IF (mat%matsize2 < 0) mat%matsize2 = mat%matsize1 !Square by default

   IF (allocated(mat%data_r)) deallocate (mat%data_r)
   IF (allocated(mat%data_c)) deallocate (mat%data_c)

   IF (mat%l_real) THEN
      ALLOCATE (mat%data_r(mat%matsize1, mat%matsize2), STAT=err, errmsg=errmsg)
      ALLOCATE (mat%data_c(0, 0))
      IF (err /= 0) then
         write (*,*) "Failed to allocate mem of shape: [" &
                    // int2str(mat%matsize1) // ", " //  int2str(mat%matsize2) // "]"
         if(present(mat_name)) then
            CALL judft_error("Allocation of memmory failed for mat datatype. Name:" // trim(mat_name), &
                                    hint="Errormessage: " // trim(errmsg))
         else
            CALL judft_error("Allocation of memmory failed for mat datatype", &
                                    hint="Errormessage: " // trim(errmsg))
         endif
      endif
      mat%data_r = 0.0
      if (present(init)) mat%data_r = init
   ELSE
      ALLOCATE (mat%data_r(0, 0))
      ALLOCATE (mat%data_c(mat%matsize1, mat%matsize2), STAT=err, errmsg=errmsg)
      IF (err /= 0) CALL judft_error("Allocation of memmory failed for mat datatype", &
                                     hint="Errormessage: " // trim(errmsg))
      mat%data_c = 0.0
      IF (PRESENT(init)) mat%data_c = init
   ENDIF
   call timestop("t_mat_alloc")
END SUBROUTINE t_mat_alloc

subroutine t_mat_init_rmatrix(mat,matrix)
   CLASS(t_mat), INTENT(OUT)      :: mat
   REAL,INTENT(IN)                :: matrix(:,:)

   mat%l_real=.true.
   mat%matsize1=size(matrix,1)
   mat%matsize2=size(matrix,2)
   mat%data_r=matrix
   allocate(mat%data_c(0,0))
end subroutine
subroutine t_mat_init_cmatrix(mat,matrix)
   CLASS(t_mat), INTENT(OUT)      :: mat
   COMPLEX,INTENT(IN)             :: matrix(:,:)

   mat%l_real=.false.
   mat%matsize1=size(matrix,1)
   mat%matsize2=size(matrix,2)
   mat%data_c=matrix
   allocate(mat%data_r(0,0))
end subroutine
subroutine t_mat_init_rmatrix_alloc(mat,matrix,alloc)
   CLASS(t_mat), INTENT(OUT)      :: mat
   REAL,INTENT(INOUT),ALLOCATABLE :: matrix(:,:)
   LOGICAL,INTENT(IN)             :: alloc

   mat%l_real=.true.
   mat%matsize1=size(matrix,1)
   mat%matsize2=size(matrix,2)
   if (alloc) THEN
      move_alloc(matrix,mat%data_r)
   else   
      mat%data_r=matrix
   endif   
   allocate(mat%data_c(0,0))
end subroutine
subroutine t_mat_init_cmatrix_alloc(mat,matrix,alloc)
   CLASS(t_mat), INTENT(OUT)      :: mat
   COMPLEX,INTENT(INOUT),ALLOCATABLE:: matrix(:,:)
   LOGICAL,INTENT(IN)             :: alloc

   mat%l_real=.false.
   mat%matsize1=size(matrix,1)
   mat%matsize2=size(matrix,2)
   if (alloc) THEN
      move_alloc(matrix,mat%data_c)
   else   
      mat%data_c=matrix
   endif   
   allocate(mat%data_r(0,0))
end subroutine

!Data movement and initialization
SUBROUTINE t_mat_clear(mat)
   IMPLICIT NONE
   CLASS(t_mat), INTENT(INOUT):: mat
   INTEGER :: i

   IF (mat%l_real) THEN
      call dlaset("A",mat%matsize1,mat%matsize2,0.0,0.0,mat%data_r,mat%matsize1)
   ELSE
      call zlaset("A",mat%matsize1,mat%matsize2,cmplx(0.0,0.0),cmplx(0.0,0.0),mat%data_c,mat%matsize1)
   ENDIF
#ifdef _OPENACC
   IF (mat%l_real) THEN
     if (acc_is_present(mat%data_r)) Then
       !$acc kernels present(mat%data_r)
       mat%data_r(:,:)=0.0
       !$acc end kernels
     endif
   ELSE
     if (acc_is_present(mat%data_c)) Then
       !$acc kernels present(mat%data_c)
       mat%data_c(:,:)=0.0
       !$acc end kernels
     endif
   ENDIF
#endif
END SUBROUTINE t_mat_clear

SUBROUTINE t_mat_copy(mat, mat1, n1, n2)
   IMPLICIT NONE
   CLASS(t_mat), INTENT(INOUT):: mat
   class(t_mat), INTENT(IN)   :: mat1
   INTEGER, INTENT(IN)        :: n1, n2

   INTEGER:: i1, i2, j1, j2
   logical:: both_on_gpu, tmp

   call timestart("t_mat_copy")

   if(.not. mat%allocated()) then
#ifdef _OPENACC
      tmp = acc_is_present(mat1%data_c)
#else 
      tmp = .False.
#endif
      if(tmp) then
         call judft_error("can't use copy alloc on GPU")
      else
         call mat%init(mat1)
      endif
   endif

#ifdef _OPENACC
   if(mat1%l_real) then
      both_on_gpu = acc_is_present(mat%data_r) .and. acc_is_present(mat1%data_r)
   else
      both_on_gpu = acc_is_present(mat%data_c) .and. acc_is_present(mat1%data_c)
   endif
#else 
   both_on_GPU = .False.
#endif 

   select type (mat1)
   type is(t_mat)

   class default
      call judft_error("you can only copy a t_mat to a t_mat")
   end select

   i1 = mat%matsize1 - n1 + 1  !space available for first dimension
   i2 = mat%matsize2 - n2 + 1
   i1 = MIN(i1, mat1%matsize1)
   i2 = MIN(i2, mat1%matsize2)

   if(both_on_GPU )then
      if(mat%l_real) then
         !$acc kernels present(mat, mat%data_r, mat1, mat1%data_r)
         do j1 = 1,i1 
            do j2 = 1,i2 
               mat%data_r(n1+j1-1, n2+j2-1) = mat1%data_r(j1,j2)
            enddo
         enddo
         !$acc end kernels
      else 
         !$acc kernels present(mat, mat%data_c, mat1, mat1%data_c)
         do j1 = 1,i1 
            do j2 = 1,i2 
               mat%data_c(n1+j1-1, n2+j2-1) = mat1%data_c(j1,j2)
            enddo
         enddo
         !$acc end kernels
      endif
   else
      IF (mat%l_real) THEN
         call dlacpy("N", i1, i2, mat1%data_r, size(mat1%data_r, 1),  mat%data_r(n1,n2), size(mat%data_r,1) )
      ELSE
         call zlacpy("N", i1, i2, mat1%data_c, size(mat1%data_c, 1),  mat%data_c(n1,n2), size(mat%data_c,1) )
      END IF
   endif

   call timestop("t_mat_copy")
END SUBROUTINE t_mat_copy


SUBROUTINE t_mat_move(mat, mat1)
   IMPLICIT NONE
   CLASS(t_mat), INTENT(INOUT):: mat
   CLASS(t_mat), INTENT(INOUT):: mat1
   !Special case, the full matrix is copied. Then use move alloc
   IF (mat%l_real) THEN
      CALL move_ALLOC(mat1%data_r, mat%data_r)
   ELSE
      CALL move_ALLOC(mat1%data_c, mat%data_c)
   END IF
END SUBROUTINE t_mat_move

SUBROUTINE t_mat_free(mat)
   use m_judft
   CLASS(t_mat), INTENT(INOUT)::mat
   call timestart("t_mat_free")
   IF (ALLOCATED(mat%data_c)) DEALLOCATE (mat%data_c)
   IF (ALLOCATED(mat%data_r)) DEALLOCATE (mat%data_r)
   mat%matsize1 = -1
   mat%matsize2 = -1
   call timestop("t_mat_free")
END SUBROUTINE t_mat_free


 !Calculations

   subroutine add(mat,mat2,alpha_c,alpha_r)
#ifdef _OPENACC
         use openacc
#endif            
         IMPLICIT NONE 
         CLASS(t_mat), INTENT(INOUT)      :: mat
         class(t_mat), INTENT(IN)         :: mat2
         complex,intent(in),optional      :: alpha_c
         complex,intent(in),optional      :: alpha_r
   
         real:: a_r
         complex:: a_c
         INTEGER:: i,j

         INTEGER,PARAMETER :: gpu_mode=1,omp_mode=2,caxpy_mode=2
         INTEGER :: mode=caxpy_mode
         
         a_c=1.0
         if (present(alpha_c)) a_c=alpha_c
      
         a_r=1.0
         if(present(alpha_r)) a_r=alpha_r
         
#ifdef _OPENACC
         if (mat%l_real) THEN 
            mode=merge(gpu_mode,mode,acc_is_present(mat%data_r).and.acc_is_present(mat2%data_r))
         else
            mode=merge(gpu_mode,mode,acc_is_present(mat%data_c).and.acc_is_present(mat2%data_c))
         endif
#endif         

         if (mat%l_real) THEN
            select case(mode)
               CASE(GPU_MODE)
                  !Data is on Device, hence we can operate on GPU
                  !$acc kernels present(mat%data_r,mat2%data_r)
                  mat%data_r=mat%data_r+a_r*mat2%data_r
                  !$acc end kernels
               CASE(OMP_MODE)   
                  !$OMP parallel do collapse(2) shared (mat,mat2,a_r) default(none)
                  DO j=1,mat%matsize2
                     DO i=1,mat%matsize1
                        mat%data_r(i,j)=mat%data_r(i,j)+a_r*mat2%data_r(i,j)
                     ENDDO
                  ENDDO
                  !$OMP end parallel do
               CASE default
               call daxpy(size(mat%data_r),a_r,mat2%data_r(1,1),1,mat%data_r(1,1),1)    
            end select   
         ELSE
            select case(mode)
               case(GPU_MODE)
                  !Data is on Device, hence we can operate on GPU
                  !$acc kernels present(mat%data_c,mat2%data_c)
                  mat%data_c=mat%data_c+a_c*mat2%data_c
                  !$acc end kernels
               case(OMP_MODE)
                  !$OMP parallel do collapse(2) shared (mat,mat2,a_c) default(none)
                  DO j=1,mat%matsize2
                     DO i=1,mat%matsize1
                      mat%data_c(i,j)=mat%data_c(i,j)+a_c*mat2%data_c(i,j)
                     ENDDO
                  ENDDO
                  !$OMP end parallel do
               CASE DEFAULT
                  call zaxpy(size(mat%data_c),a_c,mat2%data_c(1,1),1,mat%data_c(1,1),1)              
            END SELECT
         ENDIF
   END SUBROUTINE


   function t_mat_size_mb(mat) result(mb_size)
      implicit none
      class(t_mat), intent(inout) :: mat
      real :: mb_size

      if(mat%l_real) then
         mb_size =  8e-6 * size(mat%data_r)
      else
         mb_size = 16e-6 * size(mat%data_c)
      endif
   end function t_mat_size_mb

   function t_mat_norm2(mat) result(norm)
      implicit none
      class(t_mat), intent(in) :: mat
      real :: norm
      real, external :: dnrm2, dznrm2

      if (mat%l_real) then
         norm = dnrm2(size(mat%data_r), mat%data_r, 1)
      else
         norm = dznrm2(size(mat%data_c), mat%data_c, 1)
      endif
   end function t_mat_norm2

   function t_mat_allocated(mat) result(var_alloc)
      implicit none
      class(t_mat), intent(in) :: mat
      logical :: var_alloc

      var_alloc = allocated(mat%data_r) .or. allocated(mat%data_c)
   end function t_mat_allocated

   SUBROUTINE t_mat_multiply(mat1, mat2, res, transA, transB)
      use m_judft
      CLASS(t_mat), INTENT(INOUT)            :: mat1
      CLASS(t_mat), INTENT(IN)               :: mat2
      CLASS(t_mat), INTENT(INOUT), OPTIONAL    :: res
      character(len=1), intent(in), optional :: transA, transB

      integer           :: m,n,k, lda, ldb, ldc
      character(len=1)  :: transA_i, transB_i
      type(t_mat)       :: tmp
      logical           :: run_on_gpu 

      call timestart("t_mat_multiply")

      if(mat1%matsize1 == -1 .or. mat1%matsize2 == -1) call judft_error("mat1 not initialized")
      if(mat2%matsize1 == -1 .or. mat2%matsize2 == -1) call judft_error("mat2 not initialized")

      transA_i = "N"
      if(present(transA)) transA_i = transA
      transB_i = "N"
      if(present(transB)) transB_i = transB

      if(transA_i == "N") then
         m = mat1%matsize1
         k = mat1%matsize2
      else
         m = mat1%matsize2
         k = mat1%matsize1
      endif

      if(mat1%l_real .neqv. mat2%l_real) call judft_error("can only multiply matricies of the same type")

      if(mat1%l_real) then
#ifdef _OPENACC
         run_on_gpu = acc_is_present(mat1%data_r) .and. acc_is_present(mat2%data_r)
         if(present(res)) then
            run_on_gpu = run_on_gpu .and. acc_is_present(res%data_r)
         endif 
#else
         run_on_gpu = .False.
#endif
      else
#ifdef _OPENACC
         run_on_gpu = acc_is_present(mat1%data_c) .and. acc_is_present(mat2%data_c)
         if(present(res)) then
            run_on_gpu = run_on_gpu .and. acc_is_present(res%data_c)
         endif 
#else
         run_on_gpu = .False.
#endif
      endif

      if(transB_i == "N" ) then
         if(k /= mat2%matsize1) call judft_error("dimensions don't agree for matmul")
         n = mat2%matsize2
      else
         if(k /= mat2%matsize2) call judft_error("dimensions don't agree for matmul")
         n = mat2%matsize1
      endif

      lda = merge(size(mat1%data_r, dim=1), size(mat1%data_c, dim=1), mat1%l_real)
      if(transA_i == "N") then
         if(lda < max(1,m)) call judft_error("problem with lda")
      else
         if(lda < max(1,k)) call judft_error("problem with lda")
      endif

      ldb = merge(size(mat2%data_r, dim=1), size(mat2%data_c, dim=1), mat2%l_real)
      if(transB_i == "N") then
         if(ldb < max(1,k)) call judft_error("problem with ldb")
      else
         if(ldb < max(1,n)) call judft_error("problem with ldb")
      endif

      IF (present(res)) THEN
         ! prepare res matrix
         if(res%allocated()) then
            if(res%l_real .neqv. mat1%l_real) then
               call juDFT_error("res must be of the correct type")
            else
               if(res%l_real) then
                  if(any(shape(res%data_r) /= [m,n])) then
                     call juDFT_error("res must be of the correct size!")
                  endif
               else
                  if(any(shape(res%data_c) /= [m,n])) then
                     call juDFT_error("res must be of the correct size!")
                  endif
               endif
            endif
         else
            call juDFT_error("res must be allocated")
         endif

         ldc = merge(size(res%data_r, dim=1), size(res%data_c, dim=1), mat2%l_real)
         if(ldc < max(1,m)) call judft_error("problem with ldc")

         if(run_on_gpu) then
            call perform_cublas_gemm(mat1%l_real, transA_i,transB_i,m,n,k, lda, ldb, ldc,&
                                    mat1%data_r, mat1%data_c, mat2%data_r, mat2%data_c, res%data_r, res%data_c)
         else
            IF (mat1%l_real) THEN
               call dgemm(transA_i,transB_i,m,n,k, 1.0, mat1%data_r, lda, mat2%data_r, ldb, 0.0, res%data_r, ldc)
            ELSE
               call zgemm(transA_i,transB_i,m,n,k,cmplx_1, mat1%data_c, lda, mat2%data_c, ldb, cmplx_0,res%data_c, ldc)
            ENDIF
         endif
      else
         if (mat1%matsize1  /= mat1%matsize2 .or. mat2%matsize2 /= mat2%matsize1)&
            CALL judft_error("Cannot multiply matrices inplace because of non-matching dimensions", hint="This is a BUG in FLEUR, please report")

         call tmp%alloc(mat1%l_real, n,n)
         ldc = merge(size(tmp%data_r, dim=1), size(tmp%data_c, dim=1), tmp%l_real)
         if(ldc < max(1,m)) call judft_error("problem with ldc")

         if(run_on_gpu) then
            !$acc data create(tmp, tmp%data_r, tmp%data_c)
            call perform_cublas_gemm(mat1%l_real, transA_i, transB_i, m,n,k, lda, ldb, ldc,& 
                                     mat1%data_r, mat1%data_c, mat2%data_r, mat2%data_c, tmp%data_r, tmp%data_c)
            call mat1%copy(tmp,1,1)
            !$acc end data
         else
            if (mat1%l_real) THEN
               call dgemm(transA_i,transB_i,n,n,n, 1.0, mat1%data_r, lda, mat2%data_r, ldb, 0.0, tmp%data_r, ldc)
            ELSE
               call zgemm(transA_i,transB_i,n,n,n,cmplx_1, mat1%data_c, lda, mat2%data_c, ldb, cmplx_0, tmp%data_c, ldc)
            ENDIF
            call mat1%copy(tmp,1,1)
         endif
         call tmp%free()
      end IF
      call timestop("t_mat_multiply")
   end SUBROUTINE t_mat_multiply

   subroutine perform_cublas_gemm(l_real, transA, transB, m,n,k, lda, ldb, ldc,& 
                                  mat1_data_r, mat1_data_c, mat2_data_r, mat2_data_c, res_data_r, res_data_c)
      implicit none 
      logical, intent(in)           :: l_real
      character(len=*), intent(in)  :: transA, transB 
      integer, intent(in)           :: m, n, k, lda, ldb, ldc
      real, intent(in)              :: mat1_data_r(:,:), mat2_data_r(:,:)
      complex, intent(in)           :: mat1_data_c(:,:), mat2_data_c(:,:)
      real, intent(inout)           :: res_data_r(:,:)
      complex, intent(inout)        :: res_data_c(:,:)

#ifdef _OPENACC
      if(l_real) then
         !$acc host_data use_device(mat1_data_r, mat2_data_r, res_data_r)
         call cublasDgemm(transA, transB, m, n, k, 1.0, mat1_data_r, lda, mat2_data_r, ldb, 0.0, res_data_r, ldc)
         !$acc end host_data
      else
         !$acc host_data use_device(mat1_data_c, mat2_data_c, res_data_c)
         call cublasZgemm(transA, transB, m, n, k, cmplx_1, mat1_data_c, lda, mat2_data_c, ldb, cmplx_0, res_data_c, ldc)
         !$acc end host_data
      endif
#endif
   end subroutine perform_cublas_gemm


   SUBROUTINE t_mat_transpose(mat1, res,hidden)
      CLASS(t_mat), INTENT(INOUT)         :: mat1
      CLASS(t_mat), INTENT(OUT), OPTIONAL :: res
      LOGICAL,INTENT(IN),OPTIONAL         :: hidden

      if (present(hidden)) THEN
         if (hidden) THEN
          if (present(res)) THEN
            call res%init(mat1)
            call res%copy(mat1)
            res%hidden_transpose=.true.
          ELSE
            mat1%hidden_transpose=.not.mat1%hidden_transpose
          ENDIF
          return
         ENDIF 
      ENDIF         
         
      !$acc kernels present(mat%data_r,mat%data_c) if(mat%data_on_device)
      IF (present(res)) THEN
         call res%alloc(mat1%l_real, mat1%matsize2, mat1%matsize1)
         IF (mat1%l_real) THEN
            res%data_r = transpose(mat1%data_r(:mat1%matsize1, :mat1%matsize2))
         ELSE
            res%data_c = conjg(transpose(mat1%data_c(:mat1%matsize1, :mat1%matsize2)))
         ENDIF
      else
         IF (mat1%l_real) THEN
            mat1%data_r = transpose(mat1%data_r(:mat1%matsize1, :mat1%matsize2))
         ELSE
            mat1%data_c = conjg(transpose(mat1%data_c(:mat1%matsize1, :mat1%matsize2)))
         ENDIF
      end IF
      !$acc end kernels
   end SUBROUTINE t_mat_transpose

   subroutine t_mat_scale(mat,scale,hidden)
      CLASS(t_mat), INTENT(INOUT)         :: mat
      REAL,INTENT(IN)                     :: scale
      LOGICAL,INTENT(IN),OPTIONAL         :: hidden

      if (present(hidden)) THEN
         if (hidden) THEN
            mat%hidden_scaled=mat%hidden_scaled*scale
            return
         endif
      endif
      
      !$acc kernels present(mat%data_r,mat%data_c) if (mat%data_on_device)
      if (mat%l_real) THEN
         mat%data_r=mat%hidden_scaled*mat%data_r
      else
         mat%data_c=mat%hidden_scaled*mat%data_c
      endif   
      !$acc end kernels
   end subroutine t_mat_scale 

   subroutine t_mat_evaluate(mat)
      CLASS(t_mat), INTENT(INOUT)         :: mat

      IF (abs(mat%hidden_scaled-1.0)>1E-20) THEN
         mat%scale(mat%hidden_scaled)
         mat%hidden_scaled=1.0
      ENDIF
      
      IF (mat%hidden_transpose) THEN
         call mat%transpose()
         mat%hidden_transpose=.false.
      ENDIF
   END subroutine t_mat_evaluate

   subroutine t_mat_data2device(mat)
      CLASS(t_mat), INTENT(INOUT)         :: mat
    
      if (mat%l_real)THEN
         if (acc_is_present(mat%data_r))THEN
            !$acc update device(mat%data_r)
         else   
            !$acc data copyin(mat%data_r)
         endif   
      else
         if (acc_is_present(mat%data_c))THEN
            !$acc update device(mat%data_c)
         else   
            !$acc data copyin(mat%data_c)
         endif  
      endif
      mat%data_on_device=.true.
   end subroutine
   
   subroutine t_mat_data2host(mat,delete)
      CLASS(t_mat), INTENT(INOUT)         :: mat
      LOGICAL,INTENT(IN),OPTIONAL         :: delete
      if (mat%data_on_device) THEN
         if (mat%l_real)THEN
          !$acc update self(mat%data_r)
         else
          !$acc update self(mat%data_c)
         endif
      endif
      if (present(delete)) THEN
         if (delete) THEN
            if (mat%l_real) THEN
               if (acc_is_present(mat%data_r)) THEN
                  !$acc data delete(mat%data_r)
               endif
            else    
               if (acc_is_present(mat%data_c)) THEN
                  !$acc data delete(mat%data_c)
               endif
            endif
         endif
      endif           
      mat%data_on_device=.false.
   end subroutine
   

END MODULE m_types_mat
