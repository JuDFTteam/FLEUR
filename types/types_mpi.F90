!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mpi
   TYPE t_mpi
      !k-point parallelism
      INTEGER :: mpi_comm !< replaces MPI_COMM_WORLD
      INTEGER :: irank    !< rank of task in mpi_comm
      INTEGER :: isize    !< no of tasks in mpi_comm
      INTEGER, ALLOCATABLE :: k_list(:)
      INTEGER, ALLOCATABLE :: coulomb_owner(:)
      !Eigenvalue parallelism
      INTEGER :: sub_comm !< Sub-Communicator for eigenvalue parallelization (all PE working on same k-point)
      INTEGER :: n_rank   !< rank in sub_comm
      INTEGER :: n_size   !< PE per kpoint, i.e. "isize" for eigenvalue parallelization
      INTEGER, ALLOCATABLE :: ev_list(:)
      !Communicator for PE on same node
      INTEGER :: mpi_comm_same_node
      logical :: l_set_root_comm = .false. ! only create root comm once
      logical :: l_mpi_multithreaded = .false.
      integer :: root_comm ! communicator between all n_rank = 0
      !Communicator for diagonalization
      INTEGER :: diag_sub_comm
      LOGICAL :: pe_diag=.true.
   CONTAINS
      procedure :: set_errhandler    => t_mpi_set_errhandler
      procedure :: is_root => mpi_is_root
      procedure :: set_root_comm => t_mpi_set_root_comm
   END TYPE t_mpi

   INTERFACE juDFT_win_create
      MODULE PROCEDURE  juDFT_win_create_real, juDFT_win_create_cmplx, juDFT_win_create_int, &
                        juDFT_win_create_real_3D, juDFT_win_create_cmplx_3D
   END INTERFACE juDFT_win_create

   PRIVATE
   PUBLIC :: juDFT_win_create, judft_comm_split, judft_comm_split_type, t_mpi, set_root_comm, calcIndexBounds
contains
   subroutine t_mpi_set_root_comm(fmpi)
      implicit none
      class(t_mpi), intent(inout) :: fmpi

      if(.not. fmpi%l_set_root_comm ) then
            call judft_comm_split(fmpi%mpi_comm, fmpi%n_rank, 0, fmpi%root_comm)
            fmpi%l_set_root_comm = .True.
      endif
   end subroutine t_mpi_set_root_comm

   function mpi_is_root(mpi) result(is_root)
      implicit none
      class(t_mpi), intent(in) :: mpi
      logical :: is_root
      is_root = mpi%irank == 0
   end function mpi_is_root

   subroutine juDFT_win_create_real(base, size, disp_unit, info, comm, win)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      real, POINTER, ASYNCHRONOUS, intent(inout) :: base(:)
      integer, intent(in)      :: disp_unit, info, comm
      integer, intent(inout)   :: win

#ifdef CPP_MPI
      INTEGER(KIND=MPI_ADDRESS_KIND) :: SIZE
      integer                  :: err, err_handler

      call timestart("MPI_WIN_CREATE")
      CALL MPI_WIN_CREATE(base, size, disp_unit, info, comm, win, err)
      if(err /= 0) call judft_error("Can't create MPI_Win for real_data_ptr")
      call timestop("MPI_WIN_CREATE")

      call timestart("MPI_Win_create_errhandler")
      call MPI_Win_create_errhandler(judft_mpi_error_handler, err_handler, err)
      if(err /= 0) call judft_error("Can't create Error handler")
      call timestop("MPI_Win_create_errhandler")

      call timestart("MPI_WIN_SET_ERRHANDLER")
      CALL MPI_WIN_SET_ERRHANDLER(win, err_handler, err)
      if(err /= 0) call judft_error("Can't assign Error handler to Win")
      call timestop("MPI_WIN_SET_ERRHANDLER")
#else
   INTEGER :: SIZE
#endif
   end subroutine juDFT_win_create_real

   subroutine juDFT_win_create_real_3D(base, size, disp_unit, info, comm, win)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      real, POINTER, ASYNCHRONOUS, intent(inout) :: base(:,:,:)
      integer, intent(in)      :: disp_unit, info, comm
      integer, intent(inout)   :: win

#ifdef CPP_MPI
      INTEGER(KIND=MPI_ADDRESS_KIND) :: SIZE
      integer                  :: err, err_handler

      call timestart("MPI_WIN_CREATE")
      CALL MPI_WIN_CREATE(base, size, disp_unit, info, comm, win, err)
      if(err /= 0) call judft_error("Can't create MPI_Win for real_data_ptr")
      call timestop("MPI_WIN_CREATE")

      call timestart("MPI_Win_create_errhandler")
      call MPI_Win_create_errhandler(judft_mpi_error_handler, err_handler, err)
      if(err /= 0) call judft_error("Can't create Error handler")
      call timestop("MPI_Win_create_errhandler")

      call timestart("MPI_WIN_SET_ERRHANDLER")
      CALL MPI_WIN_SET_ERRHANDLER(win, err_handler, err)
      if(err /= 0) call judft_error("Can't assign Error handler to Win")
      call timestop("MPI_WIN_SET_ERRHANDLER")
#else
   INTEGER :: SIZE
#endif
   end subroutine juDFT_win_create_real_3D

   subroutine juDFT_win_create_cmplx(base, size, disp_unit, info, comm, win)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      complex, POINTER, ASYNCHRONOUS, intent(inout):: base(:)
      integer, intent(in)      :: disp_unit, info, comm
      integer, intent(inout)   :: win

#ifdef CPP_MPI
      INTEGER(KIND=MPI_ADDRESS_KIND) :: SIZE
      integer                  :: err, err_handler

      CALL MPI_WIN_CREATE(base, size, disp_unit, info, comm, win, err)
      if(err /= 0) call judft_error("Can't create MPI_Win for cmplx_data_ptr")

      call MPI_Win_create_errhandler(judft_mpi_error_handler, err_handler, err)
      if(err /= 0) call judft_error("Can't create Error handler")

      CALL MPI_WIN_SET_ERRHANDLER(win, err_handler, err)
      if(err /= 0) call judft_error("Can't assign Error handler to Win")
#else
   INTEGER :: SIZE
#endif
   end subroutine juDFT_win_create_cmplx

   subroutine juDFT_win_create_cmplx_3D(base, size, disp_unit, info, comm, win)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      complex, POINTER, ASYNCHRONOUS, intent(inout):: base(:,:,:)
      integer, intent(in)      :: disp_unit, info, comm
      integer, intent(inout)   :: win

#ifdef CPP_MPI
      INTEGER(KIND=MPI_ADDRESS_KIND) :: SIZE
      integer                  :: err, err_handler

      CALL MPI_WIN_CREATE(base, size, disp_unit, info, comm, win, err)
      if(err /= 0) call judft_error("Can't create MPI_Win for cmplx_data_ptr")

      call MPI_Win_create_errhandler(judft_mpi_error_handler, err_handler, err)
      if(err /= 0) call judft_error("Can't create Error handler")

      CALL MPI_WIN_SET_ERRHANDLER(win, err_handler, err)
      if(err /= 0) call judft_error("Can't assign Error handler to Win")
#else
   INTEGER :: SIZE
#endif
   end subroutine juDFT_win_create_cmplx_3D

   subroutine juDFT_win_create_int(base, size, disp_unit, info, comm, win)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      integer, POINTER, ASYNCHRONOUS, intent(inout) :: base(:)
      integer, intent(in)      :: disp_unit, info, comm
      integer, intent(inout)   :: win

#ifdef CPP_MPI
      INTEGER(KIND=MPI_ADDRESS_KIND) :: SIZE
      integer                  :: err, err_handler

      CALL MPI_WIN_CREATE(base, size, disp_unit, info, comm, win, err)
      if(err /= 0) call judft_error("Can't create MPI_Win for cmplx_data_ptr")

      call MPI_Win_create_errhandler(judft_mpi_error_handler, err_handler, err)
      if(err /= 0) call judft_error("Can't create Error handler")

      CALL MPI_WIN_SET_ERRHANDLER(win, err_handler, err)
      if(err /= 0) call judft_error("Can't assign Error handler to Win")
#else
   INTEGER :: SIZE
#endif
   end subroutine juDFT_win_create_int

   subroutine judft_comm_split(comm, color, key, new_comm)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      integer, intent(in)    :: comm, color, key
      integer, intent(inout) :: new_comm
#ifdef CPP_MPI
      integer                :: ierr, err_handler

      CALL MPI_COMM_SPLIT(comm,color,key,new_comm,ierr)
      if(ierr /= 0) call judft_error("Can't split comm")

      call MPI_Comm_create_errhandler(judft_mpi_error_handler, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't create Error handler")

      call MPI_Comm_Set_Errhandler(new_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to new_comm")
#endif
   end subroutine judft_comm_split

   subroutine judft_comm_split_type(comm, split_type, key, info, new_comm)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      integer, intent(in)    :: comm, split_type, key, info
      integer, intent(inout) :: new_comm
      integer                :: ierr, err_handler

#ifdef CPP_MPI
      call MPI_comm_split_type(comm, split_type, key, info, new_comm, ierr)
      if(ierr /= 0) call judft_error("Can't split comm")

      call MPI_Comm_create_errhandler(judft_mpi_error_handler, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't create Error handler")

      call MPI_Comm_Set_Errhandler(new_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to new_comm")
#endif
   end subroutine judft_comm_split_type

   subroutine t_mpi_set_errhandler(self)
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      class(t_mpi), intent(in) :: self

#ifdef CPP_MPI
      integer                  :: err_handler, ierr

      call MPI_Comm_create_errhandler(judft_mpi_error_handler, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't create Error handler")

      call MPI_Comm_Set_Errhandler(MPI_COMM_WORLD, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to MPI_COMM_WORLD")

      call MPI_Comm_Set_Errhandler(self%mpi_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to self%mpi_comm")

      call MPI_Comm_Set_Errhandler(self%sub_comm, err_handler, ierr)
      if(ierr /= 0) call judft_error("Can't assign Error handler to self%sub_comm")
#endif
   end subroutine t_mpi_set_errhandler

   subroutine judft_mpi_error_handler(comm, error_code)
#ifdef CPP_MPI
      use mpi
#endif
      use m_judft
      implicit none
      integer  :: comm, error_code
      integer             :: str_len, ierr
      character(len=3000) :: error_str

#ifdef CPP_MPI
      call MPI_ERROR_STRING(error_code, error_str, str_len, ierr)
      call judft_error("MPI failed with Error_code = " // int2str(error_code) // new_line("A") // &
                       error_str(1:str_len))
#endif
   end subroutine judft_mpi_error_handler

   SUBROUTINE calcIndexBounds(fmpi,firstIndexOverall, lastIndexOverall, firstIndexRank, lastIndexRank)

      IMPLICIT NONE

      TYPE(t_mpi), INTENT(IN)        :: fmpi
      INTEGER, INTENT(IN)            :: firstIndexOverall, lastIndexOverall
      INTEGER, INTENT(OUT)           :: firstIndexRank, lastIndexRank

      INTEGER :: chunkSize, leftoverSize, length

      length = lastIndexOverall - firstIndexOverall + 1
      chunkSize = length / fmpi%isize
      leftoverSize = MODULO(length, fmpi%isize)
      IF (fmpi%irank < leftoverSize) THEN
         firstIndexRank = fmpi%irank*(chunkSize + 1) + firstIndexOverall
         lastIndexRank = (fmpi%irank + 1)*(chunkSize + 1) + firstIndexOverall - 1
      ELSE
         firstIndexRank = leftoverSize*(chunkSize + 1) + firstIndexOverall + (fmpi%irank - leftoverSize)*chunkSize
         lastIndexRank = (firstIndexRank + chunkSize) - 1
      ENDIF
   END SUBROUTINE calcIndexBounds

END MODULE m_types_mpi
