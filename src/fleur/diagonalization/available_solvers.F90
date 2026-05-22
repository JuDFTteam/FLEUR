!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_available_solvers
   use m_types_solver
   use m_lapack
   use m_writeout
   use m_dummy_diag
   use m_scalapack
   use m_chase
   use m_elsi
   use m_magma
   use m_cuda_diag
   use m_elpa
   use m_nvlamath
   implicit none
   private
   integer :: first_real_solver = 4, num_solvers = 11
   public:: t_solver, parallel_solver_available, select_solver, print_solver, list_solvers

contains

   subroutine assign_solver(i, all_solvers)
      integer, INTENT(in) :: i
      CLASS(t_solver), allocatable, INTENT(OUT) :: all_solvers
      select case (i)
      case (1)
         all_solvers = get_solver_stop()
      case (2)
         all_solvers = get_solver_dummy()
      case (3)
         all_solvers = get_solver_debugout()
      case (4)
         all_solvers = get_solver_lapack()
      case (5)
         all_solvers = get_solver_scalapack()
      case (6)
         all_solvers = get_solver_chase()
      case (7)
         all_solvers = get_solver_elsi()
      case (8)
         all_solvers = get_solver_magma()
      case (9)
         all_solvers = get_solver_cuda()
      case (10)
         all_solvers = get_solver_elpa()
      case (11)
         all_solvers = get_solver_nvlamath()
      case default
         call judft_bug("Illegal request for solver")
      end select
   end subroutine

   function select_by_name(solver_name)
      character(len=*), intent(in):: solver_name
      CLASS(t_solver), allocatable :: select_by_name

      integer:: i
      DO i = 1, num_solvers
         call assign_solver(i, select_by_name)
         if (trim(solver_name) == trim(select_by_name%name)) return
      end do
      call judft_error("No Solver/transform could be selected:"//solver_name)
   end function

   logical function parallel_solver_available()
      integer ::i
      class(t_solver), allocatable::s

      parallel_solver_available = .false.
      !make an explit loop here
      do i = 1, num_solvers
         call assign_solver(i, s)
         parallel_solver_available = parallel_solver_available .or. (s%available .and. s%parallel)
      end do
   end function parallel_solver_available

   subroutine select_solver(parallel, gpu, single_precision, diag_solver, diag_transform)
      use m_juDFT
      logical, intent(IN)           :: parallel
      logical, intent(in), optional  :: single_precision, gpu
      class(t_solver), INTENT(OUT), allocatable  :: diag_solver, diag_transform

      logical :: use_single_precision, use_gpu, generalized, fit
      character(len=30):: diag_name, trans_name
      integer :: i

      use_single_precision = .false.
      if (present(single_precision)) use_single_precision = single_precision
      use_gpu = .false.
#ifdef _OPENACC
      use_gpu = .true.
#endif
      if (present(gpu)) use_gpu = gpu

      diag_name = trim(juDFT_string_for_argument("-diag"))
      if (len_trim(diag_name) .gt. 0) then
         !solver was specified on command line
         if (index(diag_name, "+") .gt. 0) then
            ! trensformation + standard solver
            trans_name = diag_name(:index(diag_name, "+") - 1)
            diag_name = diag_name(index(diag_name, "+") + 1:)
            diag_transform = select_by_name(trans_name)
         end if
         !check if "-sp" was given
         if (index(diag_name, "-") .gt. 0) then
            use_single_precision = trim(diag_name(index(diag_name, "-") + 1:)) .eq. "sp"
            diag_name = diag_name(:index(diag_name, "-") - 1)
         end if
         !select solver from name
         diag_solver = select_by_name(diag_name)
      else
         !defaults
         do i = first_real_solver, num_solvers
            call assign_solver(i, diag_solver)
            fit = diag_solver%available
            !Check if solver fits the requirements
            if (use_gpu) fit = fit .and. diag_solver%gpu
            if (parallel) then
               fit = fit .and. diag_solver%parallel
            else
               fit = fit .and. diag_solver%serial
            end if
            if (use_single_precision) fit = fit .and. diag_solver%single_precision
            if (fit) exit
         end do
         if (.not. fit) call judft_error("No default solver found.")
      end if

      diag_solver%use_sp = use_single_precision

      !Check if a default tansformation must be selected as well
      if ((.not. diag_solver%generalized .or. diag_solver%use_sp) .and. &
          .not. allocated(diag_transform)) &
         then
         do i = first_real_solver, num_solvers
            call assign_solver(i, diag_transform)
            fit = diag_transform%available .and. diag_transform%transform
            !Check if solver fits the requirements
            if (use_gpu) fit = fit .and. diag_transform%gpu
            if (parallel) then
               fit = fit .and. diag_transform%parallel
            else
               fit = fit .and. diag_transform%serial
            end if
            if (fit) exit
         end do
         if (.not. fit) call judft_error("No transform found")
      end if

      ! check if selected solvers are OK
      generalized = .not. allocated(diag_transform)
      if (.not. allocated(diag_transform)) then
         ! we use a generalized solver
         if (.not. diag_solver%generalized) then
            print *, diag_solver%name, diag_solver%generalized
            call judft_error("No generalized solver available")
         end if
         if (parallel) then
            if (.not. (diag_solver%parallel)) &
               call judft_error("No parallel solver available for your problem", &
                                hint="You might have selected the wrong solver with the -diag option")
         else
            if (.not. (diag_solver%serial)) &
               call judft_error("No serial solver available for your problem", &
                                hint="You might have selected the wrong solver with the -diag option")
         end if
         if (use_gpu) then
            if (.not. (diag_solver%gpu)) &
               call judft_warn("No GPU solver available for your problem", &
                               hint="You might have selected the wrong solver with the -diag option")
         end if
      else
         !we use a standard solver+transform
         if (.not. allocated(diag_transform) .and. diag_solver%standard) &
            call judft_error("No standard solver available or missing transform")
         if (parallel) then
            if (.not. (diag_solver%parallel .and. diag_transform%parallel)) &
               call judft_error("No parallel solver available for your problem", &
                                hint="You might have selected the wrong solver with the -diag option")
         else
            if (.not. (diag_solver%serial .and. diag_transform%serial)) &
               call judft_error("No serial solver available for your problem", &
                                hint="You might have selected the wrong solver with the -diag option")
         end if
         if (use_gpu) then
            if (.not. (diag_solver%gpu .and. diag_transform%gpu)) &
               call judft_warn("No GPU solver available for your problem", &
                               hint="You might have selected the wrong solver with the -diag option")
         end if
      end if
   end subroutine select_solver

   function print_solver(parallel)
      logical, intent(IN):: parallel
      character(len=30):: print_solver
      class(t_solver), allocatable:: s, t

      call select_solver(parallel, diag_solver=s, diag_transform=t)
      print_solver = trim(s%name)
      if (s%use_sp) print_solver = trim(print_solver)//"-sp"
      if (allocated(t)) print_solver = t%name//'+'//trim(print_solver)
   end function

   subroutine list_solvers()
      integer:: i
      class(t_solver), allocatable::s

      write (*, '(a)') "Hints on choosing the diagonalization method (see the docu for more details):"
      write (*, '(a)') "  The `-diag` option takes a string as an argument that:"
      write (*, '(a)') "    Either simply specifies the solver used for the generalized eigenvalue problem (NAME_GEV):"
      write (*, '(a)') "      '-diag NAME_GEV'"
      write (*, '(a)') &
         "    Or the string combines (with '+') a solver of the standard problem (NAME_STD) with a transformation (NAME_TRANS):"
      write (*, '(a)') "      '-diag NAME_TRANS+NAME_STD'"
      write (*, '(a)') ""
      write (*, '(a)') "List of solvers/transforms:"
      write (*, '(a)') "   Name   available  serial  parallel  GEV  STD STD-SP  Transform  GPU: "
      do i = 1, num_solvers
         call assign_solver(i, s)
         write(*,'(a10,4x,l,7x,l,6x,l,6x,l,3x,l,2x,l,9x,l,6x,l)') s%name,s%available,s%serial,s%parallel,s%generalized,s%standard,s%single_precision,s%transform,s%gpu
      end do
      write (*, *)
   end subroutine

end module m_available_solvers
