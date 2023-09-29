!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mpimat
   USE m_judft
   USE m_types_mat
   USE m_constants
#ifdef CPP_MPI
   USE mpi
#endif
   IMPLICIT NONE
   PRIVATE
   INTEGER, PARAMETER    :: DEFAULT_BLOCKSIZE = 64
   INTEGER, PARAMETER   :: dlen_ = 9

   !<This data-type extends the basic t_mat for distributed matrices.
   !<
   !<It stores the additional mpi_communicator and sets up a blacs grid for the matrix.
   !<This can be used to perform scalapack calls on the matrix with little additional input.
   !<The copy procedure is overwritten from t_mat to enable also redistribution of the matrix.
   TYPE t_blacsdata
      INTEGER:: no_use
      INTEGER:: mpi_com                          !> mpi-communiator over which matrix is distributed
      INTEGER:: blacs_desc(dlen_)                !> blacs descriptor
      !> 1: =1
      !> 2: context
      !> 3,4: global matrix size
      !> 5,6: block sizes
      !> 7,8: row/colum of grid for first row/colum of matrix
      !> 9: leading dimension of local matrix
      INTEGER:: npcol, nprow                     !> the number of columns/rows in the processor grid
      INTEGER:: mycol, myrow
   END TYPE t_blacsdata

   TYPE, EXTENDS(t_mat):: t_mpimat
      INTEGER                   :: global_size1, global_size2        !> this is the size of the full-matrix
      TYPE(t_blacsdata), POINTER :: blacsdata => null()
   CONTAINS
      PROCEDURE   :: copy => mpimat_copy     !<overwriten from t_mat, also performs redistribution
      PROCEDURE   :: move => mpimat_move     !<overwriten from t_mat, also performs redistribution
      PROCEDURE   :: free => mpimat_free     !<overwriten from t_mat, takes care of blacs-grids
      procedure   :: save_npy => mpimat_save_npy
      PROCEDURE   :: multiply => mpimat_multiply  !<overwriten from t_mat, takes care of blacs-grids
      PROCEDURE   :: init_details => mpimat_init
      PROCEDURE   :: init_template => mpimat_init_template     !<overwriten from t_mat, also calls alloc in t_mat
      PROCEDURE   :: add_transpose => mpimat_add_transpose !<overwriten from t_mat
      PROCEDURE   :: u2l => t_mpimat_u2l   ! construct full matrix if only upper triangle of hermitian matrix is given
      PROCEDURE   :: l2u => t_mpimat_l2u
      PROCEDURE   :: print_matrix
      PROCEDURE   :: from_non_dist
      procedure   :: to_non_dist
      PROCEDURE   :: transpose => mpimat_transpose
      procedure   :: print_type => mpimat_print_type
      PROCEDURE   :: linear_problem => t_mpimat_lproblem
      PROCEDURE   :: is_column_cyclic
      FINAL :: finalize, finalize_1d, finalize_2d, finalize_3d
   END TYPE t_mpimat

   PUBLIC t_mpimat, mingeselle

CONTAINS
   SUBROUTINE t_mpimat_lproblem(mat, vec)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(IN)   :: mat
      class(t_mat), INTENT(INOUT)   :: vec

      integer :: ipiv(mat%global_size1), info
#ifdef CPP_SCALAPACK
      call timestart("mpimat_lproblem")
      if (mat%l_real .neqv. vec%l_real) call judft_error("mat and vec need to be same kind")

      select type (vec)
      class is (t_mat)
         call judft_error("lproblem can only be solved if vec and mat are the same class")
      class is (t_mpimat)
         if (mat%l_real) then
            !call pdgesv(n,               nrhs,             a,         ia,ja,desca,                    ipiv
            call pdgesv(mat%global_size1, vec%global_size2, mat%data_r, 1, 1, mat%blacsdata%blacs_desc, ipiv, &
                        !  b,ib,jb,descb,info)
                        vec%data_r, 1, 1, vec%blacsdata%blacs_desc, info)
            if (info /= 0) call judft_error("Error in pdgesv for lproblem")
         else
            !call pzgesv(n,               nrhs,             a,          ia,ja,desca,                    ipiv,
            call pzgesv(mat%global_size1, vec%global_size2, mat%data_c, 1, 1, mat%blacsdata%blacs_desc, ipiv, &
                        !           b,         ib,jb, descb,info)
                        vec%data_c, 1, 1, vec%blacsdata%blacs_desc, info)

            if (info /= 0) call judft_error("Error in pzgesv for lproblem: "//int2str(info))
         end if
      end select
      call timestop("mpimat_lproblem")
#else
      call judft_error("no scala")
#endif
   end subroutine t_mpimat_lproblem

   subroutine mpimat_print_type(mat)
      implicit none
      CLASS(t_mpimat), INTENT(IN)     :: mat

      write (*, *) "type -> t_ mpimat"
   end subroutine mpimat_print_type

   SUBROUTINE mpimat_multiply(mat1, mat2, res, transA, transB)
      use m_judft
      CLASS(t_mpimat), INTENT(INOUT)     :: mat1
      CLASS(t_mat), INTENT(IN)           :: mat2
      CLASS(t_mat), INTENT(INOUT), OPTIONAL :: res
      character(len=1), intent(in), optional :: transA, transB

#ifdef CPP_SCALAPACK
      TYPE(t_mpimat)::m, r
      character(len=1)  :: transA_i, transB_i

      transA_i = "N"
      if (present(transA)) transA_i = transA
      transB_i = "N"
      if (present(transB)) transB_i = transB

      IF (.NOT. PRESENT(res)) CALL judft_error("BUG: in mpicase the multiply requires the optional result argument")
      SELECT TYPE (mat2)
      TYPE IS (t_mpimat)
         SELECT TYPE (res)
         TYPE is (t_mpimat)
            CALL m%init(mat1, mat2%global_size1, mat2%global_size2)
            CALL m%copy(mat2, 1, 1)
            CALL r%init(mat1, res%global_size1, res%global_size2)
            IF (mat1%l_real) THEN
               CALL pdgemm(transA_i, transB_i, mat1%global_size1, m%global_size2, mat1%global_size2, 1.0, &
                           mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, &
                           m%data_r, 1, 1, m%blacsdata%blacs_desc, 0.0, &
                           r%data_r, 1, 1, r%blacsdata%blacs_desc)
            ELSE
               CALL pzgemm(transA_i, transB_i, mat1%global_size1, m%global_size2, mat1%global_size2, cmplx_1, &
                           mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, &
                           m%data_c, 1, 1, m%blacsdata%blacs_desc, cmplx_0, &
                           r%data_c, 1, 1, r%blacsdata%blacs_desc)
            END IF
            CALL res%copy(r, 1, 1)
            CALL r%free()
            CALL m%free()
         CLASS default
            CALL judft_error("BUG in mpimat%multiply: res needs to be t_mpimat")
         END SELECT
      CLASS default
         CALL judft_error("BUG in mpimat%multiply: mat2 needs to be t_mpimat")
      END SELECT
#endif
   END SUBROUTINE mpimat_multiply

   subroutine mpimat_transpose(mat1, res)
      CLASS(t_mpimat), INTENT(INOUT) ::mat1
      CLASS(t_mat), INTENT(OUT), OPTIONAL ::res
      real, allocatable :: rd(:, :)
      complex, allocatable :: cd(:, :)
#ifdef CPP_SCALAPACK
      if (present(res)) Then
         select type (res)
         type is (t_mpimat)
            res%blacsdata = mat1%blacsdata
            res%matsize1 = mat1%matsize1
            res%matsize2 = mat1%matsize2
            res%global_size1 = mat1%global_size1
            res%global_size2 = mat1%global_size2
            res%l_real = mat1%l_real
         class default
            call judft_error("BUG in t_mpimat%transpose, wrong matrix type")
         end select
      END IF
      IF (mat1%l_real) THEN
         allocate (rd(size(mat1%data_r, 1), size(mat1%data_r, 2)))
         call pdtran(mat1%global_size1, mat1%global_size2, 1.0, mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, 0.0, rd, 1, 1, mat1%blacsdata%blacs_desc)
         if (present(res)) Then
            call move_alloc(rd, res%data_r)
         else
            call move_alloc(rd, mat1%data_r)
         end if
      ELSE
         allocate (cd(size(mat1%data_c, 1), size(mat1%data_c, 2)))
         call pztranc(mat1%global_size1, mat1%global_size2, cmplx(1.0, 0.0), mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, cmplx(0.0, 0.0), cd, 1, 1, mat1%blacsdata%blacs_desc)
         if (present(res)) Then
            call move_alloc(cd, res%data_c)
         else
            call move_alloc(cd, mat1%data_c)
         end if
      END IF
#else
      call judft_error("Not compiled with MPI")
#endif
   END subroutine

   SUBROUTINE print_matrix(mat, fileno)
#ifdef CPP_SCALAPACK
      USE mpi
#endif
      CLASS(t_mpimat), INTENT(INOUT) ::mat
      INTEGER:: fileno

#ifdef CPP_SCALAPACK
      INTEGER, EXTERNAL:: indxl2g
      CHARACTER(len=10)::filename
      INTEGER :: irank, isize, i, j, npr, npc, r, c, tmp, err, status(MPI_STATUS_SIZE)

      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, irank, err)
      CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, isize, err)

      tmp = 0

      IF (irank > 0) CALL MPI_RECV(tmp, 1, MPI_INTEGER, irank - 1, 0, mat%blacsdata%mpi_com, status, err) !lock
      WRITE (filename, "(a,i0)") "out.", fileno
      OPEN (fileno, file=filename, access='append')

      CALL blacs_gridinfo(mat%blacsdata%blacs_desc(2), npr, npc, r, c)
      DO i = 1, mat%matsize1
         DO j = 1, mat%matsize2
            IF (mat%l_real) THEN
               WRITE (fileno, "(5(i0,1x),2(f10.5,1x))") irank, i, j, indxl2g(i, mat%blacsdata%blacs_desc(5), r, 0, npr), &
                  indxl2g(j, mat%blacsdata%blacs_desc(6), c, 0, npc), mat%data_r(i, j)
            ELSE
               WRITE (fileno, "(5(i0,1x),2(f10.5,1x))") irank, i, j, indxl2g(i, mat%blacsdata%blacs_desc(5), r, 0, npr), &
                  indxl2g(j, mat%blacsdata%blacs_desc(6), c, 0, npc), mat%data_c(i, j)
            END IF
         END DO
      END DO
      CLOSE (fileno)
      IF (irank + 1 < isize) CALL MPI_SEND(tmp, 1, MPI_INTEGER, irank + 1, 0, mat%blacsdata%mpi_com, err)

#endif
   END SUBROUTINE print_matrix

   subroutine t_mpimat_l2u(mat)
#ifdef CPP_SCALAPACK
      USE mpi
#endif
      implicit none
      CLASS(t_mpimat), INTENT(INOUT) ::mat

      INTEGER :: i, j, i_glob, j_glob, myid, err, np
      COMPLEX, ALLOCATABLE:: tmp_c(:, :)
      REAL, ALLOCATABLE   :: tmp_r(:, :)
#ifdef CPP_SCALAPACK
      INTEGER, EXTERNAL    :: numroc, indxl2g  !SCALAPACK functions

      call timestart("t_mpimat_l2u")

      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, myid, err)
      CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, np, err)
      !Set lower part of matrix to zero

      DO i = 1, mat%matsize1
         DO j = 1, mat%matsize2
            ! Get global column corresponding to i and number of local rows up to
            ! and including the diagonal, these are unchanged in A
            i_glob = indxl2g(i, mat%blacsdata%blacs_desc(5), mat%blacsdata%myrow, 0, mat%blacsdata%nprow)
            j_glob = indxl2g(j, mat%blacsdata%blacs_desc(6), mat%blacsdata%mycol, 0, mat%blacsdata%npcol)

            IF (i_glob < j_glob) THEN
               IF (mat%l_real) THEN
                  mat%data_r(i, j) = 0.0
               ELSE
                  mat%data_c(i, j) = 0.0
               END IF
            elseif (i_glob == j_glob) THEN
               IF (mat%l_real) THEN
                  mat%data_r(i, j) = mat%data_r(i, j)*0.5
               ELSE
                  mat%data_c(i, j) = mat%data_c(i, j)*0.5
               END IF
            END IF
         END DO
      END DO

      IF (mat%l_real) THEN
         ALLOCATE (tmp_r(mat%matsize1, mat%matsize2))
         tmp_r = mat%data_r
      ELSE
         ALLOCATE (tmp_c(mat%matsize1, mat%matsize2))
         tmp_c = mat%data_c
      END IF
      CALL MPI_BARRIER(mat%blacsdata%mpi_com, i)
      IF (mat%l_real) THEN
         CALL pdgeadd('t', mat%global_size1, mat%global_size2, 1.0, tmp_r, 1, 1, mat%blacsdata%blacs_desc, 1.0, mat%data_r, 1, 1, mat%blacsdata%blacs_desc)
      ELSE
         CALL pzgeadd('c', mat%global_size1, mat%global_size2, CMPLX(1.0, 0.0), tmp_c, 1, 1, mat%blacsdata%blacs_desc, CMPLX(1.0, 0.0), mat%data_c, 1, 1, mat%blacsdata%blacs_desc)
      END IF
#endif
      call timestop("t_mpimat_l2u")
   end subroutine t_mpimat_l2u

   SUBROUTINE t_mpimat_u2l(mat)
#ifdef CPP_SCALAPACK
      USE mpi
#endif
      implicit none
      CLASS(t_mpimat), INTENT(INOUT) ::mat

      INTEGER :: i, j, i_glob, j_glob, myid, err, np
      COMPLEX, ALLOCATABLE:: tmp_c(:, :)
      REAL, ALLOCATABLE   :: tmp_r(:, :)
#ifdef CPP_SCALAPACK
      INTEGER, EXTERNAL    :: numroc, indxl2g  !SCALAPACK functions

      call timestart("t_mpimat_u2l")
      if (all(mat%blacsdata%blacs_desc(5:6) == [mat%global_size1, 1])) then
         !copy upper triangle to lower triangle
         call mingeselle(mat, mat)
      else
         CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, myid, err)
         CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, np, err)
         !Set lower part of matrix to zero

         DO i = 1, mat%matsize1
            DO j = 1, mat%matsize2
               ! Get global column corresponding to i and number of local rows up to
               ! and including the diagonal, these are unchanged in A
               i_glob = indxl2g(i, mat%blacsdata%blacs_desc(5), mat%blacsdata%myrow, 0, mat%blacsdata%nprow)
               j_glob = indxl2g(j, mat%blacsdata%blacs_desc(6), mat%blacsdata%mycol, 0, mat%blacsdata%npcol)

               IF (i_glob > j_glob) THEN
                  IF (mat%l_real) THEN
                     mat%data_r(i, j) = 0.0
                  ELSE
                     mat%data_c(i, j) = 0.0
                  END IF
               END IF
               IF (i_glob == j_glob) THEN
                  IF (mat%l_real) THEN
                     mat%data_r(i, j) = mat%data_r(i, j)/2.0
                  ELSE
                     mat%data_c(i, j) = mat%data_c(i, j)/2.0
                  END IF
               END IF
            END DO
         END DO

         IF (mat%l_real) THEN
            ALLOCATE (tmp_r(mat%matsize1, mat%matsize2))
            tmp_r = mat%data_r
         ELSE
            ALLOCATE (tmp_c(mat%matsize1, mat%matsize2))
            tmp_c = mat%data_c
         END IF
         CALL MPI_BARRIER(mat%blacsdata%mpi_com, i)
         IF (mat%l_real) THEN
            CALL pdgeadd('t', mat%global_size1, mat%global_size2, 1.0, tmp_r, 1, 1, mat%blacsdata%blacs_desc, 1.0, mat%data_r, 1, 1, mat%blacsdata%blacs_desc)
         ELSE
            CALL pzgeadd('c', mat%global_size1, mat%global_size2, CMPLX(1.0, 0.0), tmp_c, 1, 1, mat%blacsdata%blacs_desc, CMPLX(1.0, 0.0), mat%data_c, 1, 1, mat%blacsdata%blacs_desc)
         END IF
      end if !mingeselle
#endif
      call timestop("t_mpimat_u2l")
   END SUBROUTINE t_mpimat_u2l

   SUBROUTINE mpimat_add_transpose(mat, mat1)
      CLASS(t_mpimat), INTENT(INOUT) ::mat
      CLASS(t_mat), INTENT(INOUT) ::mat1

      INTEGER:: i, ii, n_size, n_rank

      call timestart("mpimat_add_transpose")
      SELECT TYPE (mat1)
      TYPE IS (t_mpimat)
#ifdef CPP_MPI
         CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, n_rank, i)
         CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, n_size, i)
#endif
         !Set lower part of matrix to zero...
         ii = 0
         DO i = n_rank + 1, MIN(mat%global_size1, mat%global_size2), n_size
            ii = ii + 1
            IF (mat%l_real) THEN
               mat%data_r(i + 1:, ii) = 0.0
               mat1%data_r(i + 1:, ii) = 0.0
               mat1%data_r(i, ii) = 0.0
            ELSE
               mat%data_c(i + 1:, ii) = 0.0
               mat1%data_c(i + 1:, ii) = 0.0
               mat1%data_c(i, ii) = 0.0
            END IF
         END DO
         IF (mat%l_real) THEN
#ifdef CPP_SCALAPACK

            CALL pdgeadd('t', mat1%global_size1, mat1%global_size2, 1.0, mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, 1.0, mat%data_r, 1, 1, mat%blacsdata%blacs_desc)
         ELSE
            CALL pzgeadd('c', mat1%global_size1, mat1%global_size2, CMPLX(1.0, 0.0), mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, CMPLX(1.0, 0.0), mat%data_c, 1, 1, mat1%blacsdata%blacs_desc)
#endif
         END IF
         !Now multiply the diagonal of the matrix by 1/2

         !ii=0
         !DO i=n_rank+1,MIN(mat%global_size1,mat%global_size2),n_size
         !   ii=ii+1
         !   IF (mat%l_real) THEN
         !      mat%data_r(i,ii)=mat%data_r(i,ii)/2
         !   ELSE
         !      mat%data_c(i,ii)=mat%data_c(i,ii)/2
         !   END IF
         !ENDDO
      CLASS default
         CALL judft_error("Inconsistent types in t_mpimat_add_transpose")
      END SELECT

      call timestop("mpimat_add_transpose")
   END SUBROUTINE mpimat_add_transpose

   SUBROUTINE mpimat_copy(mat, mat1, n1, n2)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)::mat
      CLASS(t_mat), INTENT(IN)      ::mat1
      INTEGER, INTENT(IN) ::n1, n2
      INTEGER :: irank, err

      call timestart("mpimat_copy")

      select type (mat1)
      type is (t_mpimat)

      class default
         call judft_error("you can only copy a t_mpimat to a t_mpimat")
      end select

#ifdef CPP_SCALAPACK
      SELECT TYPE (mat1)
      TYPE IS (t_mpimat)
         if (mat1%is_column_cyclic().and..not.mat%is_column_cyclic()) THEN
            call cyclic_column_to_2Dblock_cyclic(mat1,mat,n1,n2)
         else
            IF (mat%l_real) THEN
               CALL pdgemr2d(Mat1%global_size1, mat1%global_size2, mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, mat%data_r, n1, n2, mat%blacsdata%blacs_desc, mat1%blacsdata%blacs_desc(2))
            ELSE
               CALL pzgemr2d(mat1%global_size1, mat1%global_size2, mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, mat%data_c, n1, n2, mat%blacsdata%blacs_desc, mat1%blacsdata%blacs_desc(2))
            END IF
         endif   
      CLASS DEFAULT
         CALL judft_error("Wrong datatype in copy")
      END SELECT
#else
       call judft_error("Distributed matrix without SCALAPACK",calledby="mpimat_copy")      
#endif

      call timestop("mpimat_copy")
   END SUBROUTINE mpimat_copy

   SUBROUTINE from_non_dist(mat, mat1)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)::mat
      TYPE(t_mat), INTENT(IN)       ::mat1

      INTEGER:: blacs_desc(9), irank, ierr, umap(1, 1), np
#ifdef CPP_SCALAPACK
      blacs_desc = (/1, -1, mat1%matsize1, mat1%matsize2, mat1%matsize1, mat1%matsize2, 0, 0, mat1%matsize1/)

      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, irank, ierr)
      umap(1, 1) = 0
      CALL BLACS_GET(mat%blacsdata%blacs_desc(2), 10, blacs_desc(2))
      CALL BLACS_GRIDMAP(blacs_desc(2), umap, 1, 1, 1)
      IF (mat%l_real) THEN
         CALL pdgemr2d(Mat1%matsize1, mat1%matsize2, mat1%data_r, 1, 1, blacs_desc, mat%data_r, 1, 1, mat%blacsdata%blacs_desc, mat%blacsdata%blacs_desc(2))
      ELSE
         CALL pzgemr2d(mat1%matsize1, mat1%matsize2, mat1%data_c, 1, 1, blacs_desc, mat%data_c, 1, 1, mat%blacsdata%blacs_desc, mat%blacsdata%blacs_desc(2))
      END IF
#endif
   END SUBROUTINE from_non_dist

   subroutine to_non_dist(mat_in, mat_out)
      implicit none
      CLASS(t_mpimat), INTENT(IN)::mat_in
      TYPE(t_mat), INTENT(INOUT)       ::mat_out

      INTEGER:: blacs_desc(9), irank, ierr, umap(1, 1), np

#ifdef CPP_SCALAPACK
      blacs_desc = [1, -1, mat_out%matsize1, mat_out%matsize2, mat_out%matsize1, mat_out%matsize2, 0, 0, mat_out%matsize1]

      CALL MPI_COMM_RANK(mat_in%blacsdata%mpi_com, irank, ierr)
      umap(1, 1) = 0
      CALL BLACS_GET(mat_in%blacsdata%blacs_desc(2), 10, blacs_desc(2))
      CALL BLACS_GRIDMAP(blacs_desc(2), umap, 1, 1, 1)
      IF (mat_in%l_real) THEN
         !call pdgemr2d(m,                  n,                   a,           ia,ja,desca
         call pdgemr2d(mat_in%global_size1, mat_in%global_size2, mat_in%data_r, 1, 1, mat_in%blacsdata%blacs_desc, &
                       !             b,            ib,jb,descb,      ictxt)
                       mat_out%data_r, 1, 1, blacs_desc, mat_in%blacsdata%blacs_desc(2))
      ELSE
         !call pzgemr2d(m,                  n,                   a,           ia,ja,desca
         call pzgemr2d(mat_in%global_size1, mat_in%global_size2, mat_in%data_c, 1, 1, mat_in%blacsdata%blacs_desc, &
                       !             b,            ib,jb,descb,      ictxt)
                       mat_out%data_c, 1, 1, blacs_desc, mat_in%blacsdata%blacs_desc(2))
      END IF
#endif
   end subroutine to_non_dist

   subroutine mpimat_save_npy(mat, filename)
      use m_judft
      implicit NONE
      CLASS(t_mpimat), INTENT(IN)::mat
      character(len=*)         :: filename
      type(t_mat) :: tmp
      integer :: ierr, irank
#ifdef CPP_MPI
      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, irank, ierr)

      if (irank == 0) then
         call tmp%alloc(mat%l_real, mat%global_size1, mat%global_size2)
      else
         call tmp%alloc(mat%l_real,1,1)
      endif 

      call mat%to_non_dist(tmp)

      if (irank == 0) then
         call tmp%save_npy(filename)
      end if
      call tmp%free()
#endif
   end subroutine mpimat_save_npy

   SUBROUTINE mpimat_move(mat, mat1)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)::mat
      CLASS(t_mat), INTENT(INOUT)   ::mat1
      CALL mat%copy(mat1, 1, 1)
   END SUBROUTINE mpimat_move

   SUBROUTINE finalize(mat)
      IMPLICIT NONE
      TYPE(t_mpimat), INTENT(INOUT) :: mat
      CALL mpimat_free(mat)
   END SUBROUTINE finalize

   SUBROUTINE finalize_1d(mat)
      IMPLICIT NONE

      TYPE(t_mpimat), INTENT(INOUT) :: mat(:)
      INTEGER                      :: i
      DO i = 1, size(mat)
         CALL mpimat_free(mat(i))
      END DO
   END SUBROUTINE finalize_1d

   SUBROUTINE finalize_2d(mat)
      IMPLICIT NONE

      TYPE(t_mpimat), INTENT(INOUT) :: mat(:, :)
      INTEGER                      :: i, j

      DO i = 1, size(mat, dim=1)
         DO j = 1, size(mat, dim=2)
            CALL mpimat_free(mat(i, j))
         END DO
      END DO
   END SUBROUTINE finalize_2d

   SUBROUTINE finalize_3d(mat)
      IMPLICIT NONE

      TYPE(t_mpimat), INTENT(INOUT) :: mat(:, :, :)
      INTEGER                      :: i, j, k

      DO i = 1, size(mat, dim=1)
         DO j = 1, size(mat, dim=2)
            DO k = 1, size(mat, dim=3)
               CALL mpimat_free(mat(i, j, k))
            END DO
         END DO
      END DO
   END SUBROUTINE finalize_3d

   SUBROUTINE mpimat_free(mat)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT) :: mat
      INTEGER :: ierr
      IF (ALLOCATED(mat%data_r)) DEALLOCATE (mat%data_r)
      IF (ALLOCATED(mat%data_c)) DEALLOCATE (mat%data_c)
      IF (ASSOCIATED(mat%blacsdata)) THEN
         IF (mat%blacsdata%no_use > 1) THEN
            mat%blacsdata%no_use = mat%blacsdata%no_use - 1
            mat%blacsdata => null()
         ELSE
#ifdef CPP_SCALAPACK
            if (mat%blacsdata%blacs_desc(2) /= -1) THEN
               CALL BLACS_GRIDEXIT(mat%blacsdata%blacs_desc(2), ierr)
               DEALLOCATE (mat%blacsdata)
            endif   
#endif
         END IF
      END IF
   END SUBROUTINE mpimat_free

   !>Initialization of the distributed matrix.
  !!
  !! The argument l_2d controls the kind of distribution used:
  !!  - TRUE: the matrix is a Scalapack BLOCK-CYCLIC distribution
  !!  - FALSE: the matrix is distributed in a one-dimensional column cyclic distribution with blocksize 1
  !! as used in the parallel matrix setup of FLEUR
   SUBROUTINE mpimat_init(mat, l_real, matsize1, matsize2, mpi_subcom, l_2d, nb_x, nb_y, mat_name)
#ifdef CPP_MPI
      use mpi
#endif
      IMPLICIT NONE
      CLASS(t_mpimat)                      :: mat
      INTEGER, INTENT(IN), OPTIONAL        :: matsize1, matsize2, mpi_subcom
      LOGICAL, INTENT(IN), OPTIONAL        :: l_real, l_2d
      INTEGER, INTENT(IN), OPTIONAL        :: nb_y, nb_x
      character(len=*), intent(in), optional :: mat_name

#ifdef CPP_SCALAPACK
      INTEGER::nbx, nby, irank, ierr
      CALL mpi_comm_rank(MPI_COMM_WORLD, irank, ierr)

      call timestart("mpimat_init")
      ALLOCATE (mat%blacsdata, stat=ierr)
      if (mpi_subcom == MPI_COMM_NULL) Then
         mat%blacsdata%blacs_desc(2) = -1
         mat%global_size1 = matsize1
         mat%global_size2 = matsize2
         mat%matsize1 = 0
         mat%matsize2 = 0
      else
         nbx = DEFAULT_BLOCKSIZE; nby = DEFAULT_BLOCKSIZE
         IF (PRESENT(nb_x)) nbx = nb_x
         IF (PRESENT(nb_y)) nby = nb_y
         IF (.NOT. (PRESENT(matsize1) .AND. PRESENT(matsize2) .AND. PRESENT(mpi_subcom) .AND. PRESENT(l_real) .AND. PRESENT(l_2d))) &
            CALL judft_error("Optional arguments must be present in mpimat_init")
         mat%global_size1 = matsize1
         mat%global_size2 = matsize2
         mat%blacsdata%no_use = 1
      end if
      CALL priv_create_blacsgrid(mpi_subcom, l_2d, matsize1, matsize2, nbx, nby, &
                                 mat%blacsdata, mat%matsize1, mat%matsize2)

      mat%blacsdata%mpi_com = mpi_subcom
      CALL mat%alloc(l_real) !Attention,sizes determined in call to priv_create_blacsgrid
      !check if this matrix is actually distributed over MPI_COMM_SELF
      !IF (mpi_subcom == MPI_COMM_SELF) THEN
      !   CALL MPI_COMM_RANK(mpi_subcom, irank, ierr)
      !   IF (irank > 0) mat%blacsdata%blacs_desc(2) = -1
      !END IF

      call timestop("mpimat_init")
#else
    CALL juDFT_ERROR("No parallel matrix setup without SCALAPACK")
#endif
   END SUBROUTINE mpimat_init

   SUBROUTINE mpimat_init_template(mat, templ, global_size1, global_size2, mat_name)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)  :: mat
      CLASS(t_mat), INTENT(IN)        :: templ
      INTEGER, INTENT(IN), OPTIONAL    :: global_size1, global_size2
      character(len=*), intent(in), optional :: mat_name

      INTEGER::numroc
      EXTERNAL::numroc

      SELECT TYPE (templ)
      TYPE IS (t_mpimat)
         mat%l_real = templ%l_real
         IF (PRESENT(global_size1) .AND. PRESENT(global_size2)) THEN
            ALLOCATE (mat%blacsdata)
            mat%blacsdata = templ%blacsdata
            mat%blacsdata%no_use = mat%blacsdata%no_use + 1
            mat%blacsdata%blacs_desc(3) = global_size1
            mat%blacsdata%blacs_desc(4) = global_size2
            mat%global_size1 = global_size1
            mat%global_size2 = global_size2
#ifdef CPP_SCALAPACK
            mat%matsize1 = NUMROC(global_size1, mat%blacsdata%blacs_desc(5), mat%blacsdata%myrow, mat%blacsdata%blacs_desc(7), mat%blacsdata%nprow)
            mat%matsize2 = NUMROC(global_size2, mat%blacsdata%blacs_desc(6), mat%blacsdata%mycol, mat%blacsdata%blacs_desc(8), mat%blacsdata%npcol)
#endif
         ELSE
            mat%matsize1 = templ%matsize1
            mat%matsize2 = templ%matsize2
            mat%global_size1 = templ%global_size1
            mat%global_size2 = templ%global_size2
            mat%blacsdata => templ%blacsdata
            mat%blacsdata%no_use = mat%blacsdata%no_use + 1
         END IF
         CALL mat%alloc()

      CLASS default
         CALL judft_error("Mixed initialization in t_mpimat not possible(BUG)")
      END SELECT
   END SUBROUTINE mpimat_init_template

   SUBROUTINE priv_create_blacsgrid(mpi_subcom, l_2d, m1, m2, nbc, nbr, blacsdata, local_size1, local_size2)
#ifdef CPP_SCALAPACK
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mpi_subcom
      INTEGER, INTENT(IN) :: m1, m2
      INTEGER, INTENT(INOUT)::nbc, nbr
      LOGICAL, INTENT(IN) :: l_2d
      type(t_blacsdata), INTENT(OUT)::blacsdata
      INTEGER, INTENT(OUT):: local_size1, local_size2

#ifdef CPP_SCALAPACK
      INTEGER     :: myrowssca, mycolssca
      INTEGER     :: iamblacs, npblacs, np, myid, mycol, myrow
      INTEGER     :: nprow2, npcol2
      INTEGER     :: k, i, j, ictextblacs
      INTEGER     :: ierr

      INTEGER, ALLOCATABLE :: iblacsnums(:), ihelp(:), iusermap(:, :)

      EXTERNAL descinit, blacs_get
      EXTERNAL blacs_pinfo, blacs_gridinit

      !Determine rank and no of processors
      if (mpi_subcom /= MPI_COMM_NULL) THEN
         CALL MPI_COMM_RANK(mpi_subcom, myid, ierr)
         CALL MPI_COMM_SIZE(mpi_subcom, np, ierr)

         ! compute processor grid, as square as possible
         ! If not square with more rows than columns
         IF (l_2d) THEN
            distloop: DO j = INT(SQRT(REAL(np))), 1, -1
               IF ((np/j)*j == np) THEN
                  blacsdata%npcol = np/j
                  blacsdata%nprow = j
                  EXIT distloop
               END IF
            END DO distloop
         ELSE
            nbc = 1
            nbr = MAX(m1, m2)
            blacsdata%npcol = np
            blacsdata%nprow = 1
         END IF

         ALLOCATE (iblacsnums(np), ihelp(np), iusermap(blacsdata%nprow, blacsdata%npcol))

         !   A blacsdata%nprow*blacsdata%npcol processor grid will be created
         !   Row and column index myrow, mycol of this processor in the grid
         !   and distribution of A and B in ScaLAPACK
         !   The local processor will get myrowssca rows and mycolssca columns
         !   of A and B
         !

         myrow = myid/blacsdata%npcol  ! my row number in the BLACS blacsdata%nprow*blacsdata%npcol grid
         mycol = myid - (myid/blacsdata%npcol)*blacsdata%npcol  ! my column number in the BLACS blacsdata%nprow*blacsdata%npcol grid
         !
         !  Now allocate Asca to put the elements of Achi or receivebuffer to
         !
         myrowssca = (m1 - 1)/(nbr*blacsdata%nprow)*nbr + MIN(MAX(m1 - (m1 - 1)/(nbr*blacsdata%nprow)*nbr*blacsdata%nprow - nbr*myrow, 0), nbr)
         !     Number of rows the local process gets in ScaLAPACK distribution
         mycolssca = (m2 - 1)/(nbc*blacsdata%npcol)*nbc + MIN(MAX(m2 - (m2 - 1)/(nbc*blacsdata%npcol)*nbc*blacsdata%npcol - nbc*mycol, 0), nbc)

         !Get BLACS ranks for all MPI ranks
         CALL BLACS_PINFO(iamblacs, npblacs)  ! iamblacs = local process rank (e.g. myid)
         ! npblacs  = number of available processes
         iblacsnums = -2
         ihelp = -2
         ihelp(myid + 1) = iamblacs ! Get the Blacs id corresponding to the MPI id
         if (mpi_subcom /= MPI_COMM_NULL) &
            CALL MPI_ALLREDUCE(ihelp, iblacsnums, np, MPI_INTEGER, MPI_MAX, mpi_subcom, ierr)

         IF (ierr /= 0) call judft_error('Error in allreduce for BLACS nums')

         !     iblacsnums(i) is the BLACS-process number of MPI-process i-1
         k = 1
         DO i = 1, blacsdata%nprow
            DO j = 1, blacsdata%npcol
               iusermap(i, j) = iblacsnums(k)
               k = k + 1
            END DO
         END DO
      else
         CALL BLACS_PINFO(iamblacs, npblacs)
         iusermap = reshape([iamblacs], [1, 1])
         blacsdata%nprow = 1; blacsdata%npcol = 1; blacsdata%blacs_desc(2) = -1
      end if

!#ifdef CPP_BLACSDEFAULT
      !Get the Blacs default context
      CALL BLACS_GET(0, 0, ictextblacs)
!#else
!    ictextblacs=mpi_subcom
!#endif
      ! Create the Grid
      CALL BLACS_GRIDMAP(ictextblacs, iusermap, size(iusermap, 1), blacsdata%nprow, blacsdata%npcol)
      !     Now control, whether the BLACS grid is the one we wanted
      CALL BLACS_GRIDINFO(ictextblacs, nprow2, npcol2, blacsdata%myrow, blacsdata%mycol)
      IF (nprow2 /= blacsdata%nprow) THEN
         WRITE (oUnit, *) 'Wrong number of rows in BLACS grid'
         WRITE (oUnit, *) 'nprow=', blacsdata%nprow, ' nprow2=', nprow2
         call judft_error('Wrong number of rows in BLACS grid')
      END IF
      IF (npcol2 /= blacsdata%npcol) THEN
         WRITE (oUnit, *) 'Wrong number of columns in BLACS grid'
         WRITE (oUnit, *) 'npcol=', blacsdata%npcol, ' npcol2=', npcol2
         call judft_error('Wrong number of columns in BLACS grid')

      END IF
      !Create the descriptors
      if (mpi_subcom == MPI_COMM_NULL) then
         blacsdata%blacs_desc(2) = -1
         local_size1 = 0
         local_size2 = 0
      else
         IF (myrowssca==0.or.mycolssca==0) THEN
            CALL juDFT_warn("With your chosen eigenvalue parallelization scheme some MPI processes have no share of the Hamiltonian. Please check your parallelization.",hint="Either reduce the number of MPI processes or increase the k-point parallelization.")
         END IF
         CALL descinit(blacsdata%blacs_desc, m1, m2, nbr, nbc, 0, 0, ictextblacs, myrowssca, ierr)
         IF (ierr /= 0) call judft_error('Creation of BLACS descriptor failed')
         local_size1 = myrowssca
         local_size2 = mycolssca
      end if

#endif
   END SUBROUTINE priv_create_blacsgrid

   SUBROUTINE mingeselle(mat_in, mat_out)
      !---------------------------------------------------------------------+
      !                                                                     |
      ! Transfers the spin-down/spin-up part , upper triangle of the        |
      ! mat_in to the lower triangle of the mat_out.                        |
      !                                                                     |
      !      --------------------                                           |
      !      |      |  mat_out  |                                           |
      !      |      |           | nv1                                       |
      !      --------------------                                           |
      !      |      |           |                                           |
      !      |mat_in|           |                                           |
      !      |      |           |                                           |
      !      |      |           | nv2                                       |
      !      --------------------                                           |
      !                                                                     |
      ! For eigenvector-parallelization this needs some communication       |
      ! between the nodes, since this part is created 'column-wise'         |
      ! but needed row-wise.                                                |
      !                                                                     |
      !     n_send(i): number of elements to send to pe #i                  |
      !     n_recv(i): number of elements to receive from pe #i             |
      !     ns_tot,nr_tot : total number of elements to send/receive        |
      !     cs_el,cr_el: send and receive elements                          |
      !     in_pos(xy,n,i): where to put in the n'th element sent by pe #i  |
      !                                                                     |
      !  ToDo: chop the send messages to the smaller chunks                 |
      !                                                                     |
      !  Based on the old mingeselle.                                       |
      !                                    Okt. 2020                        |
      !                                    U.Alekseeva                      |
      !---------------------------------------------------------------------+

      CLASS(t_mat), INTENT(INOUT) ::mat_in
      CLASS(t_mat), INTENT(INOUT) ::mat_out
      ! ..
      ! .. Local Scalarsxi
      INTEGER n_rank, n_size
      INTEGER ki, kj, kjj, ns_tot, nr_tot, n_p, n_help, i, ii
      INTEGER ns_max, nr_max, np_s, np_r
      INTEGER inext, ifront, req_s, req_r
      INTEGER SUB_COMM
      ! ..
      ! .. Local Arrays
      INTEGER ierr
      INTEGER, ALLOCATABLE :: c_help_size(:, :)
      INTEGER, ALLOCATABLE :: n_send(:), nsr(:)
      INTEGER, ALLOCATABLE :: n_recv(:), n_r(:)
      INTEGER, ALLOCATABLE :: in_pos(:, :, :)
      COMPLEX, ALLOCATABLE :: cs_el(:, :), cr_el(:), b_b(:), c_help(:, :)
      LOGICAL, ALLOCATABLE :: nsl(:)

#ifdef CPP_MPI
      INTEGER stt(MPI_STATUS_SIZE)

      call timestart("mingeselle")
      IF (mat_out%l_real .OR. mat_in%l_real) CALL juDFT_error("Matrices should be complex", calledby="mingeselle")

      SELECT TYPE (mat_in)
      TYPE IS (t_mpimat)
         SELECT TYPE (mat_out)
         TYPE IS (t_mpimat)
            IF ((mat_in%global_size1 .NE. mat_out%global_size2) .OR. (mat_in%global_size2 .NE. mat_out%global_size1)) THEN
               CALL juDFT_error("The matrices do not match", calledby="mingeselle")
            END IF
            CALL MPI_COMM_RANK(mat_out%blacsdata%mpi_com, n_rank, ierr)
            CALL MPI_COMM_SIZE(mat_out%blacsdata%mpi_com, n_size, ierr)
            SUB_COMM = mat_out%blacsdata%mpi_com

            ALLOCATE (n_send(0:n_size - 1), n_recv(0:n_size - 1), n_r(0:n_size - 1), nsr(0:n_size - 1))
            ALLOCATE (c_help_size(2, 0:n_size - 1), nsl(0:n_size - 1))
            ns_tot = 0
            nr_tot = 0
            n_send = 0
            n_r = 0
            c_help_size = 0

            ! determine number of elements to send to other pe's
            ! and calculate the dimensions of c_helpi
            ! rows of c_help correspond to columns of mat_in and vice versa

            DO ki = 1, mat_in%matsize2
               kjj = n_rank + 1 + (ki - 1)*n_size      ! global column index
               nsr = 0
               nsl = .FALSE.
               DO kj = 1, min(kjj - 1, mat_in%matsize1)
                  ns_tot = ns_tot + 1
                  n_p = MOD(kj - 1, n_size)
                  n_send(n_p) = n_send(n_p) + 1
                  nsr(n_p) = nsr(n_p) + 1
                  nsl(n_p) = .TRUE.
               END DO
               DO n_p = 0, n_size - 1
                  IF (c_help_size(2, n_p) < nsr(n_p)) c_help_size(2, n_p) = nsr(n_p)
                  IF (nsl(n_p)) c_help_size(1, n_p) = c_help_size(1, n_p) + 1
               END DO
            END DO
            !print*, "send", n_rank, ns_tot, n_send

            ! determine number of elements to receive from other pe's

            DO ki = 1, mat_out%matsize2
               kjj = n_rank + 1 + (ki - 1)*n_size      ! global column index
               DO kj = kjj + 1, mat_out%matsize1
                  nr_tot = nr_tot + 1
                  n_p = MOD(kj - 1, n_size)
                  n_r(n_p) = n_r(n_p) + 1
               END DO
            END DO
            !print*, "recv", n_rank, nr_tot, n_r

            ! determine the maximal number of s/r-counts and allocate s/r-arrays

            ns_max = 0
            nr_max = 0
            DO n_p = 0, n_size - 1
               ns_max = MAX(ns_max, n_send(n_p))
               nr_max = MAX(nr_max, n_r(n_p))
            END DO
            !      WRITE (*,*) ns_max ,nr_max  , n_size, n_rank
            ALLOCATE (cs_el(ns_max, 0:n_size - 1), cr_el(nr_max))
            ALLOCATE (in_pos(2, nr_max, 0:n_size - 1))

            ! for every send destination:
            ! put the elements of the mat_in into the c_help,
            ! resorting them on the way: rows <-> columns
            ! then put them in the send buffers

            ALLOCATE (c_help(mat_in%matsize2, ceiling(real(mat_in%matsize1)/n_size)))
            c_help = cmplx(0.0, 0.0)
            DO n_p = 0, n_size - 1

               IF (c_help_size(2, n_p) > size(c_help, 2)) CALL juDFT_error("allocated c_help is too small", calledby="mingeselle")
               IF (c_help_size(1, n_p) > size(c_help, 1)) CALL juDFT_error("allocated c_help is too small", calledby="mingeselle")
               !print*, "c_help_size",n_rank, n_p,c_help_size(:,n_p)

               DO ki = 1, c_help_size(1, n_p)
                  DO kj = 1, min(ki, c_help_size(2, n_p))
                     kjj = (kj - 1)*n_size + n_p + 1     ! #row of the element in mat_in
                     IF (n_rank - 1 < n_p) THEN
                        c_help(ki, kj) = mat_in%data_c(kjj, ki + 1)
                     ELSE
                        c_help(ki, kj) = mat_in%data_c(kjj, ki)
                     END IF
                  END DO
               END DO

               n_help = 0
               DO kj = 1, c_help_size(2, n_p)
                  DO ki = kj, c_help_size(1, n_p)
                     n_help = n_help + 1
                     cs_el(n_help, n_p) = CONJG(c_help(ki, kj))
                  END DO
               END DO
               IF (n_help .NE. n_send(n_p)) CALL juDFT_error("Number of elements to send is wrong", calledby="mingeselle")

            END DO
            DEALLOCATE (c_help)

            ! now we look where to put in the received elements

            n_recv = 0
            DO ki = 1, mat_out%matsize2
               kjj = n_rank + 1 + (ki - 1)*n_size      ! global column index
               DO kj = kjj + 1, mat_out%matsize1
                  n_p = MOD(kj - 1, n_size)
                  n_recv(n_p) = n_recv(n_p) + 1
                  in_pos(1, n_recv(n_p), n_p) = kj
                  in_pos(2, n_recv(n_p), n_p) = ki
               END DO
            END DO
            DO n_p = 0, n_size - 1
               IF (n_recv(n_p) /= n_r(n_p)) CALL juDFT_error("n_recv.NE.n_s", calledby="mingeselle")
            END DO

            ! Mandaliet, mandaliet, min geselle kumme niet

            ifront = ibefore(n_size, n_rank)
            inext = iafter(n_size, n_rank)
            DO n_p = 0, n_size - 1

               ! determine pe's to send to and to receive from

               np_s = MOD(inext + n_p, n_size)
               np_r = MOD(ifront - n_p, n_size)
               IF (np_r .LT. 0) np_r = np_r + n_size

               ! send section: local rows i with mod(i-1,np) = np_s will be sent to proc np_s

               IF (np_s .NE. n_rank) THEN
                  CALL MPI_ISEND(cs_el(1, np_s), n_send(np_s), MPI_DOUBLE_COMPLEX, &
                                 np_s, n_rank, SUB_COMM, req_s, ierr)
               END IF

               ! receive section : local rows i  with mod(i-1,np) = np_r will be received from np_r
               ! ... skipped, if update matrix from local data:

               IF (np_r .NE. n_rank) THEN
                  CALL MPI_IRECV(cr_el, n_recv(np_r), MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, np_r, SUB_COMM, req_r, ierr)
                  CALL MPI_WAIT(req_s, stt, ierr)
                  CALL MPI_WAIT(req_r, stt, ierr)
                  DO ki = 1, n_recv(np_r)
                     mat_out%data_c(in_pos(1, ki, np_r), in_pos(2, ki, np_r)) = cr_el(ki)
                  END DO
               ELSE
                  DO ki = 1, n_recv(np_r)
                     mat_out%data_c(in_pos(1, ki, np_r), in_pos(2, ki, np_r)) = cs_el(ki, np_s)
                  END DO
               END IF
            END DO

            DEALLOCATE (cs_el, cr_el, in_pos)

         CLASS DEFAULT
            call judft_error("Wrong type (1) in mingeselle")
         END SELECT
      CLASS DEFAULT
         call judft_error("Wrong type (2) in mingeselle")
      END SELECT

      call timestop("mingeselle")
#endif
   END SUBROUTINE mingeselle
   !
   !-------------------------------------------------------------
   !
   INTEGER FUNCTION ibefore(np, p)
      !
      ! Determine (in a ring structure) which is the front process
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np  !  number of processes
      INTEGER, INTENT(IN) :: p   !  current processes

      IF (p > 0) THEN
         ibefore = p - 1
      ELSE
         ibefore = np - 1
      END IF

   END FUNCTION ibefore
   !
   !-------------------------------------------------------------
   !
   INTEGER FUNCTION iafter(np, p)
      !
      ! Determine (in a ring structure) which is the next process
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np  !  number of processes
      INTEGER, INTENT(IN) :: p   !  current processes

      IF (p < np - 1) THEN
         iafter = p + 1
      ELSE
         iafter = 0
      END IF

   END FUNCTION iafter


   subroutine cyclic_column_to_2Dblock_cyclic(mat,mat2d,s1,s2)
      implicit none 
      class(t_mpimat),intent(in)   ::mat
      class(t_mpimat),intent(inout)::mat2d
      integer,intent(in),optional  ::s1,s2

      real,allocatable::vecr(:)
      complex,allocatable::vecc(:)
      integer:: my_proc,num_proc,ierr
      integer:: nprow,npcol,myrow,mycol
      integer:: shift1,shift2
      integer:: n1,n2,blockindex,n_col,n_row

      shift1=0;shift2=0
      if (present(s1)) shift1=s1-1
      if (present(s2)) shift2=s2-1
      
      if (mat%l_real) THEN
         allocate(vecr(mat%global_size1))
      else 
         allocate(vecc(mat%global_size1)) 
      endif
#ifdef CPP_MPI      
      !process ranks for cyclic column dist
      call MPI_COMM_RANK(mat%blacsdata%mpi_com,my_proc,ierr)
      call MPI_COMM_SIZE(mat%blacsdata%mpi_com,num_proc,ierr)
      !info for 2dblock cyclic dist
      call blacs_gridinfo(mat2d%blacsdata%blacs_desc,nprow,npcol,myrow,mycol)

#ifdef __NEW_CODE      
      !create processor map blacs_row,blacs_col->mpi_rank
      allocate(pmap(nprow,npcol))
      pmap=-1
      pmap(myrow,mycol)=mp_proc
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,pmap,size(pmap),MPI_INTEGER,MPI_MAX,mat%blacsdata%mpi_com,ierr)
      if (any(pmap<0)) call judft_bug("Bug1:types_mpimat")

      DO n2=1,mat%global_size2,num_proc
         sglobal_col=n2+my_proc
         slocal_col=(n2-1)/num_proc+1
         !Which 2d column contains the data?
         blockindex=(sglobal_col+shift1-1)/mat2d%blacsdata%blacs_desc(6)
         rec_p_col=mod(blockindex,npcol)
         !now loop over column
         DO n1=1,mat%global_size1
            blockindex=(n1+shift1-1)/mat2d%blacsdata%blacs_desc(5)
            rec_p_row=mod(blockindex,nprow)
            rec_r_index(rec_p_row)=rec_r_index(rec_p_row)+1
            vecr(rec_r_index(rec_p_row),rec_p_row)=mat%data_r(n1,slocal_col)
         ENDDO
         !Send data to all columns of the processor grid involved
         DO rec_p_col=0,npcol-1
            call MPI_ISEND(vecr(1,rec_p_row),rec_r_index(rec_p_row),MPI_DOUBLE,pmap(rec_p_col,rec_p_col),sGlobal_col,mat%blacsdata%mpi_com,request(rec_p_col),ierr)   
         ENDDO
         !now each processor might have data from more columns in 2D dist
         DO n2_2d=n2,n2+num_proc
            rglobal_col=n2_2d+shift1
            !Which 2d column contains the data?
            blockindex=(rglobal_col-1)/mat2d%blacsdata%blacs_desc(6)
            rec_p_col=mod(blockindex,npcol)
            if (mycol.ne.rec_p_col) cycle !this process does not contain data
            !Where is the column comming from?
            s_col=mod(n2_2d-1,num_proc)
            !Which is the first element we store the data in?
            blockindex=(shift1-1)/mat2d%blacsdata%blacs_desc(5)
            if (myrow==0) n_row=shift1-blockindex*mat2d%blacsdata%blacs_desc(5) !first block might be incomplete
            n_row=n_row+blockindex/nprow*mat2d%blacsdata%blacs_desc(5)
            !Get the data
            CALL MPI_RECV(mat2%data_r(n_row:,rec_col),mat2%matsize2-n_row,MPI_DOUBLE,s_col,mat%blacsdata%mpi_com,rglobal_col,ierr)
         ENDDO
      ENDDO      
#endif      
      !Now we loop over columns and BC them
      DO n2=1,mat%global_size2   
         if (mat%l_real) THEN
            if (my_proc==mod(n2-1,num_proc)) vecr=mat%data_r(:,(n2-1)/num_proc+1) !This process owns the column
            call MPI_BCAST(vecr,size(vecr),MPI_DOUBLE,mod(n2-1,num_proc),mat%blacsdata%mpi_com,ierr)   
         else
            if (my_proc==mod(n2-1,num_proc)) vecc=mat%data_c(:,(n2-1)/num_proc+1)
            call MPI_BCAST(vecc,size(vecc),MPI_DOUBLE_COMPLEX,mod(n2-1,num_proc),mat%blacsdata%mpi_com,ierr)   
         endif    
         !Which 2d column contains the data?
         blockindex=(n2+shift2-1)/mat2d%blacsdata%blacs_desc(6)
         if (mycol.ne.mod(blockindex,npcol)) cycle !This process contains no data  
         n_col=(n2+shift2)-blockindex*mat2d%blacsdata%blacs_desc(6)+ &
               blockindex/npcol*mat2d%blacsdata%blacs_desc(6)
         DO n1=1,mat%global_size1
            blockindex=(n1+shift1-1)/mat2d%blacsdata%blacs_desc(5)
            if (myrow.ne.mod(blockindex,nprow)) cycle !No data here
            n_row=(n1+shift1)-blockindex*mat2d%blacsdata%blacs_desc(5)+ &
                    blockindex/nprow*mat2d%blacsdata%blacs_desc(5)
            if (mat%l_real) then 
               mat2d%data_r(n_row,n_col)=vecr(n1)
            else   
               mat2d%data_c(n_row,n_col)=vecc(n1)
            
            end if 
         enddo
      ENDDO   
#endif      
   end subroutine

   logical function is_column_cyclic(mat)
      implicit none 
      class(t_mpimat),intent(in)   ::mat

      is_column_cyclic=(mat%blacsdata%blacs_desc(5)==mat%global_size1.and.mat%blacsdata%blacs_desc(6)==1)
   end function

END MODULE m_types_mpimat
