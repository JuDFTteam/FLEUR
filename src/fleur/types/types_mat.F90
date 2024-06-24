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
      LOGICAL :: l_real                     !>Store either real or complex data
      INTEGER :: matsize1 = -1                !> matsize1=size(data_?,1),i.e. no of rows
      INTEGER :: matsize2 = -1                !> matsize2=size(data_?,2),i.e. no of columns
      REAL, ALLOCATABLE    :: data_r(:, :)
      COMPLEX, ALLOCATABLE :: data_c(:, :)
   CONTAINS
      PROCEDURE        :: alloc => t_mat_alloc                !> allocate the data-arrays
      PROCEDURE        :: multiply => t_mat_multiply            !> do a matrix-matrix multiply
      PROCEDURE        :: transpose => t_mat_transpose          !> transpose the matrix
      PROCEDURE        :: from_packed_real => t_mat_from_packed_real
      PROCEDURE        :: from_packed_cmplx => t_mat_from_packed_cmplx !> initialized from a packed-storage matrix
      generic          :: from_packed => from_packed_real, from_packed_cmplx
      PROCEDURE        :: inverse => t_mat_inverse             !> invert the matrix
      PROCEDURE        :: linear_problem => t_mat_lproblem    !> Solve linear equation
      PROCEDURE        :: to_packed => t_mat_to_packed          !> convert to packed-storage matrix
      PROCEDURE        :: clear => t_mat_clear                !> set data arrays to zero
      PROCEDURE        :: copy => t_mat_copy                  !> copy into another t_mat (overloaded for t_mpimat)
      PROCEDURE        :: move => t_mat_move                  !> move data into another t_mat (overloaded for t_mpimat)
      PROCEDURE        :: save_npy => t_mat_save_npy
      procedure        :: allocated => t_mat_allocated
      PROCEDURE        :: init_details => t_mat_init
      PROCEDURE        :: init_template => t_mat_init_template              !> initalize the matrix(overloaded for t_mpimat)
      GENERIC          :: init => init_details, init_template
      PROCEDURE        :: free => t_mat_free                  !> dealloc the data (overloaded for t_mpimat)
      PROCEDURE        :: add_transpose => t_mat_add_transpose!> add the tranpose/Hermitian conjg. without the diagonal (overloaded for t_mpimat)
      PROCEDURE        :: unsymmetry => t_mat_unsymmetry
      procedure        :: norm2 => t_mat_norm2
      procedure        :: subtract => t_mat_subtract
      procedure        :: u2l => t_mat_u2l
      procedure        :: l2u => t_mat_l2u
      procedure        :: size_mb => t_mat_size_mb
      procedure        :: print_type => t_mat_print_type
      procedure        :: conjugate => t_mat_conjg
      procedure        :: reset => t_mat_reset
      procedure        :: bcast => t_mat_bcast
      procedure        :: pos_eigvec_sum => t_mat_pos_eigvec_sum
      procedure        :: leastsq => t_mat_leastsq
      procedure        :: add
   END type t_mat
   PUBLIC t_mat
CONTAINS
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

         INTEGER,PARAMETER :: gpu_mode=1,omp_mode=2,caxpy_mode=3
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

   subroutine t_mat_leastsq(A, b)
      use m_constants
      implicit none
      class(t_mat), intent(inout) :: A
      type(t_mat), intent(inout)  :: b

      type(t_mat) :: tmp
      integer              :: m, n, nrhs, lda, ldb, info, lwork

      real    :: rwork_req(1)
      complex :: cwork_req(1)

      real, allocatable    :: rwork(:)
      complex, allocatable :: cwork(:)

      if(A%matsize2 /= b%matsize2) call judft_error("least-squares dimension problem")
      if(A%l_real .neqv. b%l_real) call judft_error("least-squares kind problem")

      m = A%matsize1
      n = A%matsize2
      nrhs = b%matsize2
      if(A%l_real) then
         lda = size(A%data_r,1)
         ldb = size(b%data_r,1)

         call dgels("N", m, n, nrhs, A%data_r, lda, b%data_r, ldb, rwork_req, -1, info)
         lwork = int(rwork_req(1))
         allocate(rwork(lwork), source=0.0)

         call dgels("N", m, n, nrhs, A%data_r, lda, b%data_r, ldb, rwork, lwork, info)
      else
         lda = size(A%data_c,1)
         ldb = size(b%data_c,1)

         call zgels("N", m, n, nrhs, A%data_c, lda, b%data_c, ldb, cwork_req, -1, info)
         lwork = int(cwork_req(1))
         allocate(cwork(lwork), source=cmplx_0)

         call zgels("N", m, n, nrhs, A%data_c, lda, b%data_c, ldb, cwork, lwork, info)
      endif

      if(info /= 0) call judft_error("least squares failed.")

      call tmp%init(A%l_real, n, nrhs)

      if(tmp%l_real) then
         tmp%data_r = b%data_r(:n,:)
      else
         tmp%data_c = b%data_c(:n,:)
      endif

      call b%free()
      call b%init(tmp)
      call b%copy(tmp, 1,1)
      call tmp%free()
   end subroutine t_mat_leastsq

   subroutine t_mat_pos_eigvec_sum(mat)
      implicit none
      CLASS(t_mat), INTENT(INOUT)   :: mat
      integer :: ne, i
      real    :: sum_sign_r

      do ne = 1,size(mat%data_r,2)
         i = -1
         sum_sign_r = 0.0
         do while (abs(sum_sign_r) < 1e-8)
            i = i + 1
            if(mat%matsize1-i < 1) exit
            sum_sign_r = sum(mat%data_r(:mat%matsize1-i,ne))
         enddo
         if(mat%matsize1-i >= 1) then
            mat%data_r(:,ne) = mat%data_r(:,ne) / sign(1.0, sum_sign_r)
         endif
      enddo
   end subroutine t_mat_pos_eigvec_sum

   subroutine t_mat_bcast(mat, root, comm)
      use m_divide_most_evenly 
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      CLASS(t_mat), INTENT(INOUT)   :: mat
      integer, intent(in)           :: root, comm

      integer              :: ierr, full_shape(2), me, n_parts, i
      integer, allocatable :: start_idx(:), psize(:)
      integer(8) :: sz_in_byte

#ifdef CPP_MPI
      call MPI_Comm_rank(comm, me, ierr)
      call MPI_Bcast(mat%l_real, 1, MPI_LOGICAL, root, comm, ierr)
      !alloc mat same as root
      if(me == root) then
         full_shape = merge(shape(mat%data_r), shape(mat%data_c), mat%l_real)
         call MPI_Bcast(full_shape, 2, MPI_INTEGER, root, comm, ierr)
      else
         call MPI_Bcast(full_shape, 2, MPI_INTEGER, root, comm, ierr)
         call mat%alloc(mat%l_real, full_shape(1), full_shape(2))
      endif

      ! overwrite matsize as needed
      call MPI_Bcast(mat%matsize1, 1, MPI_INTEGER, root, comm, ierr)
      call MPI_Bcast(mat%matsize2, 1, MPI_INTEGER, root, comm, ierr)

      sz_in_byte = full_shape(1)
      sz_in_byte = sz_in_byte * full_shape(2) 
      sz_in_byte = sz_in_byte * merge(8, 16, mat%l_real)
      !make sure everything is smaller than 4 GB
      n_parts = ceiling(sz_in_byte / 4e9) 
      call divide_most_evenly(mat%matsize2, n_parts, start_idx, psize)

      do i = 1,n_parts 
         if(mat%l_real) then
            call MPI_bcast(mat%data_r(:,start_idx(i)), full_shape(1)*psize(i), MPI_DOUBLE_PRECISION, root, comm, ierr)
         else
            call MPI_bcast(mat%data_c(:,start_idx(i)), full_shape(1)*psize(i), MPI_DOUBLE_COMPLEX, root, comm, ierr)
         endif
      enddo
#endif
   end subroutine t_mat_bcast

   subroutine t_mat_reset(mat, val)
      implicit none
      CLASS(t_mat), INTENT(INOUT)   :: mat
      complex, intent(in)           :: val

      if(mat%l_real) then
         mat%data_r = real(val)
      else
         mat%data_c = val
      endif
   end subroutine t_mat_reset

   subroutine t_mat_conjg(mat)
      implicit none
      CLASS(t_mat), INTENT(INOUT) :: mat
      integer :: i,j

      if(.not. mat%l_real) then
         if(mat%matsize1 == size(mat%data_c,1) .and. mat%matsize2 == size(mat%data_c,2)) then
            call zlacgv(mat%matsize1 * mat%matsize2, mat%data_c, 1)
         else
            !$OMP parallel do default(none) private(i) shared(mat)
            do i =1,mat%matsize2
               call zlacgv(mat%matsize1, mat%data_c(1,i), 1)
            enddo
            !$OMP end parallel do
         endif
      endif
   end subroutine t_mat_conjg


   subroutine t_mat_print_type(mat)
      implicit none
      CLASS(t_mat), INTENT(IN)     :: mat

      write (*,*) "type -> t_mat"
   end subroutine t_mat_print_type

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
   ! copy lower triangle to upper triangle
   subroutine t_mat_l2u(mat)
      implicit none
      class(t_mat), intent(inout) :: mat
      integer :: i,j

      call timestart("copy lower to upper matrix")
      if(mat%matsize1 /= mat%matsize2) call judft_error("l2u only works for square matricies")

      if(mat%l_real) then
         do i = 1,mat%matsize1
            do j = 1,i-1
               mat%data_r(j,i) = mat%data_r(i,j)
            enddo
         enddo
      else
         do i = 1,mat%matsize1
            do j = 1,i-1
               mat%data_c(j,i) = conjg(mat%data_c(i,j))
            enddo
         enddo
      endif
      call timestop("copy lower to upper matrix")
   end subroutine t_mat_l2u

   ! copy upper triangle to lower triangle
   subroutine t_mat_u2l(mat)
      use m_judft
      implicit none
      class(t_mat), intent(inout) :: mat
      integer :: i,j

      call timestart("copy upper to lower matrix")
      if(mat%matsize1 /= mat%matsize2) call judft_error("l2u only works for square matricies")
      if(mat%l_real) then
         do i = 1,mat%matsize1
            do j = 1,i-1
               mat%data_r(i,j) = mat%data_r(j,i)
            enddo
         enddo
      else
         do i = 1,mat%matsize1
            do j = 1,i-1
               mat%data_c(i,j) = conjg(mat%data_c(j,i))
            enddo
         enddo
      endif
      call timestop("copy upper to lower matrix")
   end subroutine t_mat_u2l

   subroutine t_mat_subtract(res_mat, mat1, mat2)
      implicit none
      class(t_mat), intent(inout) :: res_mat
      type(t_mat), intent(in)     :: mat1, mat2
      logical :: real_res
      integer :: s1, s2

      ! check dimensions
      if(mat1%matsize1 /= mat2%matsize1) call judft_error("matsize 1 doesn't agree")
      s1 = mat1%matsize1
      if(mat1%matsize2 /= mat2%matsize2) call judft_error("matsize 2 doesn't agree")
      s2 = mat1%matsize2

      ! check real/cmplx
      real_res = mat1%l_real .and. mat2%l_real
      if(res_mat%l_real .neqv. real_res) then
         call res_mat%free()
      endif
      if(.not. res_mat%allocated())   call res_mat%alloc(real_res, s1, s2)

      if(res_mat%l_real) then
         res_mat%data_r = mat1%data_r(:s1,:s2) - mat2%data_r(:s1,:s2)
      elseif(mat1%l_real .and. (.not. mat2%l_real)) then
         res_mat%data_c = mat1%data_r(:s1,:s2) - mat2%data_c(:s1,:s2)
      elseif((.not. mat1%l_real) .and. mat2%l_real) then
         res_mat%data_c = mat1%data_c(:s1,:s2) - mat2%data_r(:s1,:s2)
      else
         res_mat%data_c(:s1,:s2) = mat1%data_c(:s1,:s2) - mat2%data_c(:s1,:s2)
      endif
   end subroutine t_mat_subtract

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

   SUBROUTINE t_mat_lproblem(mat, vec)
      IMPLICIT NONE
      CLASS(t_mat), INTENT(IN)     :: mat
      class(t_mat), INTENT(INOUT)   :: vec

      INTEGER:: lwork, info
      REAL, ALLOCATABLE:: work(:)
      INTEGER, allocatable::ipiv(:)
      logical :: both_on_gpu
#ifdef _OPENACC 
      integer :: ierr, sz
      real(8), dimension(100,100) :: A
      type(cusolverDnHandle) :: handle
#endif

      select type (vec) 
      class is (t_mat)
      class default
         call judft_error("lproblem can only be solved if vec and mat are the same class")
      end select

      IF ((mat%l_real .NEQV. vec%l_real) .OR. (mat%matsize1 .NE. mat%matsize2) .OR. (mat%matsize1 .NE. vec%matsize1)) then
         CALL judft_error("Invalid matices in t_mat_lproblem")
      endif 

#ifdef _OPENACC
      if(mat%l_real) then
         both_on_gpu = acc_is_present(mat%data_r) .and. acc_is_present(vec%data_r)
      else
         both_on_gpu = acc_is_present(mat%data_c) .and. acc_is_present(vec%data_c)
      endif
#else 
      both_on_gpu = .False.
#endif

      if(both_on_gpu) then
#ifdef _OPENACC       
         allocate(ipiv(mat%matsize1))
         sz = mat%matsize1
         ierr = cusolverDnCreate(handle)
         if(ierr /= 0) call juDFT_error("can't create handle")
         
         !$acc data create(ipiv)
            call perform_LU_cusolver(handle, mat%l_real, mat%data_r, mat%data_c, ipiv)
            call calc_rhs(handle, mat%l_real, mat%data_r, mat%data_c, vec%data_r, vec%data_c, ipiv)
         !$acc end data
         
         ierr = cusolverDnDestroy(handle)
         deallocate(ipiv)
         if(ierr /= 0) call juDFT_error("can't destroy handle")
#endif
      else
         IF (mat%l_real) THEN
            call timestart("solve real linear problem")
            IF (mat%unsymmetry() < 1E-8) THEN
               !Matrix is symmetric
               CALL DPOSV('Upper', mat%matsize1, vec%matsize2, mat%data_r, mat%matsize1, &
                        vec%data_r, vec%matsize1, INFO)
               IF (INFO > 0) THEN
                  !Matrix was not positive definite
                  lwork = -1; ALLOCATE (work(1))
                  CALL DSYSV('Upper', mat%matsize1, vec%matsize2, mat%data_r, mat%matsize1, IPIV, &
                           vec%data_r, vec%matsize1, WORK, LWORK, INFO)
                  lwork = INT(work(1))
                  DEALLOCATE (work); ALLOCATE (ipiv(mat%matsize1), work(lwork))
                  CALL DSYSV('Upper', mat%matsize1, vec%matsize2, mat%data_r, mat%matsize1, IPIV, &
                           vec%data_r, vec%matsize1, WORK, LWORK, INFO)
                  IF (info .NE. 0) CALL judft_error("Could not solve linear equation, matrix singular")
               END IF
            ELSE
               allocate (ipiv(mat%matsize1))
               call dgesv(mat%matsize1, vec%matsize2, mat%data_r, mat%matsize1, ipiv, vec%data_r, vec%matsize1, info)
               if (info /= 0) call judft_error("Error in dgesv for lproblem")
            END IF
            call timestop("solve real linear problem")
         ELSE
            call timestart("solve cmplx linear problem")
            ! I don't to do the whole testing for symmetry:
            allocate (ipiv(mat%matsize1))
            call zgesv(mat%matsize1, vec%matsize2, mat%data_c, mat%matsize1, ipiv, vec%data_c, vec%matsize1, info)
            if (info /= 0) call judft_error("Error in zgesv for lproblem")
            call timestop("solve cmplx linear problem")
         ENDIF
      endif
   END SUBROUTINE t_mat_lproblem

#ifdef _OPENACC  
   subroutine perform_LU_cusolver(handle, l_real, data_r, data_c, ipiv)
      implicit none
      type(cusolverDnHandle), intent(in) :: handle
      logical, intent(in)                :: l_real 
      real, intent(inout)                :: data_r(:,:)
      complex, intent(inout)             :: data_c(:,:)
      integer, intent(inout)             :: ipiv(:)

      real, allocatable    :: r_work(:)
      complex, allocatable :: c_work(:)
      integer :: lwork, devinfo, sz, ierr

      lwork = get_lwork_cusolver(handle, l_real, data_r, data_c) 
      if(l_real) then 
         sz = size(data_r,1)
         allocate(r_work(lwork), stat=ierr)
         if(ierr /= 0) call juDFT_error("cant' alloc r_work")

         !$acc data create(r_work) copyout(devinfo)
            !$acc host_data use_device(data_r, r_work, ipiv, devinfo)
            ierr = cusolverDnDgetrf(handle, sz, sz, data_r, sz, r_work, ipiv, devinfo)
            !$acc end host_data 
         !$acc end data
         if(ierr /= 0) call juDFT_error("cusolver failed R")
         deallocate(r_work)
      else
         sz = size(data_c,1)
         allocate(c_work(lwork), stat=ierr)
         if(ierr /= 0) call juDFT_error("cant' alloc c_work")

         !$acc data create(c_work) copyout(devinfo)
            !$acc host_data use_device(data_c, c_work, ipiv, devinfo)
            ierr = cusolverDnZgetrf(handle, sz, sz, data_c, sz, c_work, ipiv, devinfo)
            !$acc end host_data 
         !$acc end data
         if(ierr /= 0) call juDFT_error("cusolver failed C")
         deallocate(c_work)
      endif 
   end subroutine perform_LU_cusolver
   
   subroutine calc_rhs(handle, l_real, mat_r, mat_c, vec_r, vec_c, ipiv)
      implicit none
      type(cusolverDnHandle), intent(in) :: handle
      logical, intent(in)                :: l_real 
      real, intent(inout)                :: mat_r(:,:), vec_r(:,:)
      complex, intent(inout)             :: mat_c(:,:), vec_c(:,:)
      integer, intent(in)                :: ipiv(:)

      integer :: ierr, n, nrhs, devinfo

      !$acc data copyout(devinfo)
         if(l_real) then 
            n = size(mat_r, 1) 
            nrhs = size(vec_r, 2)
            !$acc host_data use_device(mat_r, vec_r, ipiv, devinfo)
            ierr = cusolverDnDgetrs(handle, CUBLAS_OP_N, n, nrhs, mat_r, n, ipiv, vec_r, n, devinfo)
            !$acc end host_data
         else 
            n = size(mat_c, 1) 
            nrhs = size(vec_c, 2)
            !$acc host_data use_device(mat_c, vec_c, ipiv, devinfo)
            ierr = cusolverDnZgetrs(handle, CUBLAS_OP_N, n, nrhs, mat_c, n, ipiv, vec_c, n, devinfo)
            !$acc end host_data
         endif
      !$acc end data
      if(ierr /= 0) call judft_error("rhs failed")
   end subroutine calc_rhs 
      
   function get_lwork_cusolver(handle, l_real, data_r, data_c) result(lwork)
      implicit none 
      type(cusolverDnHandle), intent(in) :: handle
      logical, intent(in)                :: l_real 
      real, intent(in)                   :: data_r(:,:)
      complex, intent(in)                :: data_c(:,:)
      integer :: lwork, sz, ierr

      if(l_real) then
         sz = size(data_r,1)
         !$acc host_data use_device(data_r)
         ierr = cusolverDnDgetrf_buffersize(handle, sz, sz, data_r, sz, lwork)
         !$acc end host_data
      else
         sz = size(data_c,1)
         !$acc host_data use_device(data_c)
         ierr = cusolverDnZgetrf_buffersize(handle, sz, sz, data_c, sz, lwork)
         !$acc end host_data
      endif

      if(ierr /= 0) call juDFT_error("can't get lwork")
   end function get_lwork_cusolver
#endif

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

   SUBROUTINE t_mat_add_transpose(mat, mat1)
      CLASS(t_mat), INTENT(INOUT)::mat, mat1
      INTEGER::i, j
      IF ((mat%matsize1 .NE. mat1%matsize2) .OR. &
          (mat%matsize2 .NE. mat1%matsize1)) &
         CALL judft_error("Matrix sizes missmatch in add_transpose")
      IF (mat%l_real .AND. mat1%l_real) THEN
         DO i = 1, mat%matsize2
            DO j = i + 1, mat%matsize1
               mat%data_r(j, i) = mat1%data_r(i, j)
            ENDDO
         ENDDO
      ELSEIF ((.NOT. mat%l_real) .AND. (.NOT. mat1%l_real)) THEN
         DO i = 1, mat%matsize2
            DO j = i + 1, mat%matsize1
               mat%data_c(j, i) = CONJG(mat1%data_c(i, j))
            ENDDO
         ENDDO
      ELSE
         call judft_error("Inconsistency between data types in m_mat")
      END IF
   END SUBROUTINE t_mat_add_transpose

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


   SUBROUTINE t_mat_transpose(mat1, res)
      CLASS(t_mat), INTENT(INOUT)       ::mat1
      CLASS(t_mat), INTENT(OUT), OPTIONAL ::res

      IF (present(res)) THEN
         call res%alloc(mat1%l_real, mat1%matsize2, mat1%matsize1)
         IF (mat1%l_real) THEN
            res%data_r = transpose(mat1%data_r(:mat1%matsize1, :mat1%matsize2))
         ELSE
            res%data_c = conjg(transpose(mat1%data_c(:mat1%matsize1, :mat1%matsize2)))
         ENDIF
      else
         if (mat1%matsize1 .ne. mat1%matsize2) CALL judft_error("Cannot transpose matrices inplace because of non-matching dimensions", hint="This is a BUG in FLEUR, please report")
         IF (mat1%l_real) THEN
            mat1%data_r(:mat1%matsize1, :mat1%matsize2) = transpose(mat1%data_r(:mat1%matsize1, :mat1%matsize2))
         ELSE
            mat1%data_c(:mat1%matsize1, :mat1%matsize2) = conjg(transpose(mat1%data_c(:mat1%matsize1, :mat1%matsize2)))
         ENDIF
      end IF
   end SUBROUTINE t_mat_transpose

   SUBROUTINE t_mat_from_packed_real(mat1, matsize, packed_r)
      use m_judft
      CLASS(t_mat), INTENT(INOUT)       :: mat1
      INTEGER, INTENT(IN)               :: matsize
      REAL, INTENT(IN)                  :: packed_r(:)

      INTEGER:: n, nn, i
      call timestart("t_mat_from_packed_real")
      call mat1%alloc(.true., matsize, matsize)

      !$OMP PARALLEL DO default(none) &
      !$OMP shared(matsize, mat1, packed_r) private(n, nn, i) &
      !$OMP schedule(dynamic, 10)
      DO n = 1, matsize
         DO nn = 1, n
            i = ((n - 1)*n)/2 + nn
            mat1%data_r(n, nn) = packed_r(i)
            mat1%data_r(nn, n) = packed_r(i)
         end DO
      end DO
      !$OMP END PARALLEL DO
      call timestop("t_mat_from_packed_real")
   end SUBROUTINE t_mat_from_packed_real

   SUBROUTINE t_mat_from_packed_cmplx(mat1, matsize, packed_c)
      use m_judft
      CLASS(t_mat), INTENT(INOUT)       :: mat1
      INTEGER, INTENT(IN)               :: matsize
      COMPLEX, INTENT(IN)               :: packed_c(:)

      INTEGER:: n, nn, i
      call timestart("t_mat_from_packed_cmplx")
      call mat1%alloc(.false., matsize, matsize)

      !$OMP PARALLEL DO default(none) &
      !$OMP shared(matsize, mat1, packed_c) private(n, nn, i) &
      !$OMP schedule(dynamic, 10)
      DO n = 1, matsize
         DO nn = 1, n
            i = ((n - 1)*n)/2 + nn
            mat1%data_c(n, nn) = conjg(packed_c(i))
            mat1%data_c(nn, n) = packed_c(i)
            i = i + 1
         end DO
      end DO
      !$OMP END PARALLEL DO
      call timestop("t_mat_from_packed_cmplx")
   end SUBROUTINE t_mat_from_packed_cmplx

   function t_mat_to_packed(mat) result(packed)
      CLASS(t_mat), INTENT(IN)      :: mat
      COMPLEX                       :: packed(mat%matsize1*(mat%matsize1 + 1)/2)
      integer :: n, nn, i
      real, parameter :: tol = 1e-5
      if (mat%matsize1 .ne. mat%matsize2) call judft_error("Could not pack no-square matrix", hint='This is a BUG, please report')

      if (mat%l_real) THEN
         !$OMP PARALLEL DO default(none) &
         !$OMP shared(mat, packed) private(n, nn, i) &
         !$OMP schedule(dynamic, 10)
         DO n = 1, mat%matsize1
            DO nn = 1, n
               i = ((n - 1)*n)/2 + nn
               packed(i) = (mat%data_r(n, nn) + mat%data_r(nn, n))/2.
            end DO
         end DO
         !$OMP END PARALLEL DO
      else
         !$OMP PARALLEL DO default(none) &
         !$OMP shared(mat, packed) private(n, nn, i) &
         !$OMP schedule(dynamic, 10)
         DO n = 1, mat%matsize1
            DO nn = 1, n
               i = ((n - 1)*n)/2 + nn
               packed(i) = (conjg(mat%data_c(n, nn)) + mat%data_c(nn, n))/2.
            end DO
         end DO
         !$OMP END PARALLEL DO
      endif
   end function t_mat_to_packed

   subroutine t_mat_inverse(mat)
      implicit none
      CLASS(t_mat), INTENT(INOUT)       :: mat
      integer                :: info
      real, allocatable      :: work_r(:)
      integer, allocatable   :: ipiv(:)
      complex, allocatable    :: work_c(:)

      call timestart("invert matrix")
      if (mat%matsize1 .ne. mat%matsize2) then
         call judft_error("Can only invert square matrices", &
                          hint="This is a BUG in FLEUR, please report")
      endif
      ALLOCATE (ipiv(mat%matsize1))

      if (mat%l_real) THEN
         ALLOCATE (work_r(mat%matsize1))
         call dgetrf(mat%matsize1, mat%matsize1, mat%data_r, size(mat%data_r, 1), ipiv, info)
         if (info .ne. 0) call judft_error("Failed to invert matrix: dpotrf failed.")
         call dgetri(mat%matsize1, mat%data_r, size(mat%data_r, 1), ipiv, work_r, size(work_r), info)
         if (info .ne. 0) call judft_error("Failed to invert matrix: dpotrf failed.")
      else
         ALLOCATE (work_c(mat%matsize1))
         call zgetrf(mat%matsize1, mat%matsize1, mat%data_c, size(mat%data_c, 1), ipiv, info)
         if (info .ne. 0) call judft_error("Failed to invert matrix: dpotrf failed.")
         call zgetri(mat%matsize1, mat%data_c, size(mat%data_c, 1), ipiv, work_c, size(work_c), info)
         if (info .ne. 0) call judft_error("Failed to invert matrix: dpotrf failed.")
      end if
      call timestop("invert matrix")
   end subroutine t_mat_inverse

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

   subroutine t_mat_save_npy(mat, filename)
      use m_judft
      implicit NONE
      class(t_mat), intent(in) :: mat
      character(len=*)         :: filename

      if (mat%l_real) then
         !call save_npy(filename, mat%data_r)
      else
         !call save_npy(filename, mat%data_c)
      endif
   end subroutine t_mat_save_npy

   function t_mat_unsymmetry(mat) result(unsymmetry)
      implicit none
      class(t_mat), intent(in) :: mat
      real                     :: unsymmetry
      integer                  :: n

      unsymmetry = 0.0

      if (mat%matsize1 /= mat%matsize2) then
         call judft_error("Rectangular matricies can't be symmetric")
      else
         n = mat%matsize1
         if (mat%l_real) THEN
            unsymmetry = maxval(mat%data_r(:n, :n) - transpose(mat%data_r(:n, :n)))
         else
            unsymmetry = maxval(abs(mat%data_c(:n, :n) - conjg(transpose(mat%data_c(:n, :n)))))
         endif
      endif
   end function t_mat_unsymmetry

END MODULE m_types_mat
