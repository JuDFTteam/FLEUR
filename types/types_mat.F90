MODULE m_types_mat
#include"cpp_double.h"
   USE m_judft
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
      REAL, ALLOCATABLE CPP_MANAGED    :: data_r(:, :)
      COMPLEX, ALLOCATABLE CPP_MANAGED :: data_c(:, :)
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
   END type t_mat
   PUBLIC t_mat
CONTAINS
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
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      CLASS(t_mat), INTENT(INOUT)   :: mat
      integer, intent(in)           :: root, comm

      integer :: ierr, full_shape(2), me

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

      if(mat%l_real) then
         call MPI_bcast(mat%data_r, product(full_shape), MPI_DOUBLE_PRECISION, root, comm, ierr)
      else
         call MPI_bcast(mat%data_c, product(full_shape), MPI_DOUBLE_COMPLEX, root, comm, ierr)
      endif
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
      TYPE(t_mat), INTENT(INOUT)   :: vec

      INTEGER:: lwork, info
      REAL, ALLOCATABLE:: work(:)
      INTEGER, allocatable::ipiv(:)

      IF ((mat%l_real .NEQV. vec%l_real) .OR. (mat%matsize1 .NE. mat%matsize2) &
          .OR. (mat%matsize1 .NE. vec%matsize1)) &
         CALL judft_error("Invalid matices in t_mat_lproblem")
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
   END SUBROUTINE t_mat_lproblem

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
      use m_constants
      CLASS(t_mat), INTENT(INOUT)            :: mat1
      CLASS(t_mat), INTENT(IN)               :: mat2
      CLASS(t_mat), INTENT(INOUT), OPTIONAL    :: res
      character(len=1), intent(in), optional :: transA, transB

      integer           :: m,n,k, lda, ldb, ldc
      character(len=1)  :: transA_i, transB_i
      type(t_mat)       :: tmp

      call timestart("t_mat_multiply")

      if(mat1%matsize1 == -1 .and. mat1%matsize2 == -1) call judft_error("mat1 not initialized")
      if(mat2%matsize1 == -1 .and. mat2%matsize2 == -1) call judft_error("mat2 not initialized")

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
               call res%free()
            else
               if(res%l_real) then
                  if(any(shape(res%data_r) < [m,n])) then
                     call res%free()
                  else
                     res%data_r = 0.0
                     res%matsize1 = m
                     res%matsize2 = n
                  endif
               else
                  if(any(shape(res%data_c) < [m,n])) then
                     call res%free()
                  else
                     res%data_c = cmplx_0
                     res%matsize1 = m
                     res%matsize2 = n
                  endif
               endif
            endif
         endif
         if(.not. res%allocated()) call res%alloc(mat1%l_real, m,n)

         ldc = merge(size(res%data_r, dim=1), size(res%data_c, dim=1), mat2%l_real)
         if(ldc < max(1,m)) call judft_error("problem with ldc")

         IF (mat1%l_real) THEN
            call dgemm(transA_i,transB_i,m,n,k, 1.0, mat1%data_r, lda, mat2%data_r, ldb, 0.0, res%data_r, ldc)
         ELSE
            call zgemm(transA_i,transB_i,m,n,k,cmplx_1, mat1%data_c, lda, mat2%data_c, ldb, cmplx_0,res%data_c, ldc)
         ENDIF
      else
         if (mat1%matsize1  /= mat1%matsize2 .or. mat2%matsize2 /= mat2%matsize1)&
            CALL judft_error("Cannot multiply matrices inplace because of non-matching dimensions", hint="This is a BUG in FLEUR, please report")

         call tmp%alloc(mat1%l_real, n,n)
         ldc = merge(size(tmp%data_r, dim=1), size(tmp%data_c, dim=1), tmp%l_real)
         if(ldc < max(1,m)) call judft_error("problem with ldc")

         if (mat1%l_real) THEN
            call dgemm(transA_i,transB_i,n,n,n, 1.0, mat1%data_r, lda, mat2%data_r, ldb, 0.0, tmp%data_r, ldc)
         ELSE
            call zgemm(transA_i,transB_i,n,n,n,cmplx_1, mat1%data_c, lda, mat2%data_c, ldb, cmplx_0, tmp%data_c, ldc)
         ENDIF
         call mat1%copy(tmp,1,1)
         call tmp%free()
      end IF
      call timestop("t_mat_multiply")
   end SUBROUTINE t_mat_multiply

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

      INTEGER:: i1, i2

      call timestart("t_mat_copy")

      select type (mat1)
      type is(t_mat)

      class default
         call judft_error("you can only copy a t_mat to a t_mat")
      end select

      i1 = mat%matsize1 - n1 + 1  !space available for first dimension
      i2 = mat%matsize2 - n2 + 1
      i1 = MIN(i1, mat1%matsize1)
      i2 = MIN(i2, mat1%matsize2)
      IF (mat%l_real) THEN
         call dlacpy("N", i1, i2, mat1%data_r, size(mat1%data_r, 1),  mat%data_r(n1,n2), size(mat%data_r,1) )
      ELSE
         call zlacpy("N", i1, i2, mat1%data_c, size(mat1%data_c, 1),  mat%data_c(n1,n2), size(mat%data_c,1) )
      END IF

      call timestop("t_mat_copy")
   END SUBROUTINE t_mat_copy

   SUBROUTINE t_mat_clear(mat)
#ifdef _OPENACC
    use openacc
#endif
      IMPLICIT NONE
      CLASS(t_mat), INTENT(INOUT):: mat
      INTEGER :: i

      IF (mat%l_real) THEN
         call CPP_LAPACK_slaset("A",mat%matsize1,mat%matsize2,0.0,0.0,mat%data_r,mat%matsize1)
      ELSE
         call CPP_LAPACK_claset("A",mat%matsize1,mat%matsize2,cmplx(0.0,0.0),cmplx(0.0,0.0),mat%data_c,mat%matsize1)
      ENDIF
#ifdef _OPENACC
      IF (mat%l_real) THEN
        if (acc_is_present(mat%data_r)) Then
          !$acc kernels present(mat%data_r)
          mat%data_r=0.0
          !$acc end kernels
        endif
      ELSE
        if (acc_is_present(mat%data_c)) Then
          !$acc kernels present(mat%data_c)
          mat%data_c=0.0
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
         call save_npy(filename, mat%data_r)
      else
         call save_npy(filename, mat%data_c)
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
