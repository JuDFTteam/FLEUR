! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
!--------------------------------------------------------------------------------
module m_elpa
   use m_judft
   use m_constants
   use m_types_solver
   implicit none
   private
   type, extends(t_solver):: t_solver_elpa
   contains
      procedure        :: solve_gev => elpa_gev  !solver for generalized eigenvalue problem
   end type

   public get_solver_elpa

contains

   function get_solver_elpa() result(solver)
      type(t_solver_elpa), pointer::solver
      allocate (solver)
      solver%name = "elpa"
#ifdef CPP_ELPA
      solver%available = .true.
#else
      solver%available = .false.
#endif
      solver%parallel = .true.
      solver%serial = .false.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .true.
      solver%use_sp = .false.
   end function

   subroutine elpa_gev(self, hmat, smat, ne, eig, zmat)
      use m_types_mat
      use m_judft

      implicit none
      class(t_solver_elpa)                  :: self
      class(t_mat), intent(INOUT) :: hmat, smat
      integer, intent(INOUT) :: ne
      class(t_mat), allocatable, intent(OUT)   :: zmat
      real, intent(OUT)   :: eig(:)
#ifdef CPP_ELPA
      integer            :: nev, lwork, liwork, n
      integer            :: info
      integer           :: num, np, myid
      integer           :: err
      integer           :: i
      real, allocatable      :: eig2(:)
      type(t_mpimat)        :: ev_dist
      integer               :: kernel
      class(elpa_t), pointer :: elpa_obj
      logical :: firstcall = .true.

      call timestart("ELPA GEV")
      select type (hmat)
      type IS (t_mpimat)
         select type (smat)
         type IS (t_mpimat)
            call MPI_BARRIER(hmat%blacsdata%mpi_com, err)
            call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, err)
            call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, err)
            if (firstcall) then
               err = elpa_init(20180525)
               firstcall = .false.
            end if
            elpa_obj => elpa_allocate()

            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
            if (err .ne. 0) call juDFT_error('Failed to allocated "eig2"', calledby='elpa')

            call ev_dist%init(hmat)! Eigenvectors
            if (err .ne. 0) call juDFT_error('Failed to allocated "ev_dist"', calledby='elpa')

            ! Blocking factor
            if (hmat%blacsdata%blacs_desc(5) .ne. hmat%blacsdata%blacs_desc(6)) &
               call judft_error("Different block sizes for rows/columns not supported")
            call elpa_obj%set("na", hmat%global_size1, err)
            call elpa_obj%set("nev", ne, err)
            call elpa_obj%set("local_nrows", hmat%matsize1, err)
            call elpa_obj%set("local_ncols", hmat%matsize2, err)
            call elpa_obj%set("nblk", hmat%blacsdata%blacs_desc(5), err)
            call elpa_obj%set("mpi_comm_parent", hmat%blacsdata%mpi_com, err)
            call elpa_obj%set("process_row", hmat%blacsdata%myrow, err)
            call elpa_obj%set("process_col", hmat%blacsdata%mycol, err)
            call elpa_obj%set("blacs_context", hmat%blacsdata%blacs_desc(2), err)
            call elpa_obj%set("timings", 1, err)
            err = elpa_obj%setup()

#if defined(CPP_GPU)||defined(_OPENACC)
            call elpa_obj%set("gpu_hermitian_multiply", 1, err)
            !call elpa_obj%set("cannon_for_generalized",0,err)
            call elpa_obj%set("nvidia-gpu", 1, err)
            call elpa_obj%setup_gpu()
            print *, "ELPA for GPU"
            if (myid .eq. 0) call elpa_obj%store_settings("save_to_disk.txt", err)
#else
            call elpa_obj%set("solver", ELPA_SOLVER_2STAGE)
#endif

            call hmat%u2l()
            call smat%u2l()
            call elpa_obj%timer_start("ELPA")
            if (hmat%l_real) then
               call elpa_obj%generalized_eigenvectors(hmat%data_r, smat%data_r, eig2, ev_dist%data_r, .false., err)
            else
               call elpa_obj%generalized_eigenvectors(hmat%data_c, smat%data_c, eig2, ev_dist%data_c, .false., err)
            end if
            call elpa_obj%timer_stop("ELPA")
            if (myid .eq. 0) call elpa_obj%print_times("ELPA")
            call MPI_BARRIER(hmat%blacsdata%mpi_com, err)
            call elpa_deallocate(elpa_obj)
            !CALL elpa_uninit()
            ! END of ELPA stuff
            !
            !     Each process has all eigenvalues in output
            eig(:ne) = eig2(:ne)
            deallocate (eig2)
            !
            !
            !     Redistribute eigenvectors  from ScaLAPACK distribution to each process, i.e. for
            !     process i these are eigenvectors i+1, np+i+1, 2*np+i+1...
            !     Only num=num2/np eigenvectors per process
            !
            num = ne
            ne = 0
            do i = myid + 1, num, np
               ne = ne + 1
            end do
            !
            !
            allocate (t_mpimat::ev)
            call ev%init(hmat%l_real, hmat%global_size1, hmat%global_size1, hmat%blacsdata%mpi_com, .false.)
            call ev%copy(ev_dist, 1, 1)
         class DEFAULT
            call judft_error("Wrong type (1) in scalapack")
         end select
      class DEFAULT
         call judft_error("Wrong type (2) in scalapack")
      end select

      call timestop("ELPA GEV")
#endif

   end subroutine elpa_gev

end module m_elpa
