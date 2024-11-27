!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
   implicit none
   private
   type:: solverlist
      class(t_solver), pointer::s
   end type
   type(solverlist), allocatable :: all_solvers(:)
   integer :: first_real_solver = 4
   public:: t_solver, parallel_solver_available, select_solver, print_solver, list_solvers

contains

   subroutine init()
      if (allocated(all_solvers)) return
      allocate (all_solvers(10))
      all_solvers(1)%s => get_solver_stop()
      all_solvers(2)%s => get_solver_dummy()
      all_solvers(3)%s => get_solver_debugout()!These solvers are not "real solvers and should not be selected automatically
      all_solvers(4)%s => get_solver_lapack()
      all_solvers(5)%s => get_solver_scalapack()
      all_solvers(6)%s => get_solver_chase()
      all_solvers(7)%s => get_solver_elsi()
      all_solvers(8)%s => get_solver_magma()
      all_solvers(9)%s => get_solver_cuda()
      all_solvers(10)%s => get_solver_elpa()
      !t_solver_elpa()/)

   end subroutine

   logical function parallel_solver_available()
      integer ::i
      call init()
      parallel_solver_available = .false.
      !make an explit loop here
      do i = 1, size(all_solvers)
         parallel_solver_available = parallel_solver_available .or. (all_solvers(i)%s%available .and. all_solvers(i)%s%parallel)
      end do
   end function parallel_solver_available

   function select_solver(parallel, gpu, single_precision) result(diag_solver)
      use m_juDFT
      logical, intent(IN)           :: parallel
      logical, intent(in), optional  :: single_precision, gpu
      class(t_solver), pointer  :: diag_solver

      logical :: use_single_precision, use_gpu, generalized, fit
      character(len=30):: name, trans
      integer :: i
      class(t_solver), pointer :: use_transform
      call init()
      diag_solver => null()
      use_transform => null()

      use_single_precision = .false.
      if (present(single_precision)) use_single_precision = single_precision
      use_gpu = .false.
#ifdef _OPENACC
      use_gpu = .true.
#endif
      if (present(gpu)) use_gpu = gpu

      name = trim(juDFT_string_for_argument("-diag"))

      if (len(name) .gt. 0) then
         !solver was specified on command line
         if (index("+", name) .gt. 0) then
            ! trensformation + standard solver
            trans = name(:index("+", name))
            name = name(index("+", name) + 1:)
            do i = 1, size(all_solvers)
               if (all_solvers(i)%s%name .eq. trans) use_transform => all_solvers(i)%s
            end do
            if (.not. associated(use_transform)) call judft_error("Transformation not found: "//trans)
         end if
         !check if "-sp" was given
         if (index("-", name) .gt. 0) then
            use_single_precision = name(index("-", name) + 1:) .eq. "sp"
            name = name(:index("-", name))
         end if
         !select solver from name
         do i = 1, size(all_solvers)
            if (all_solvers(i)%s%name .eq. name) diag_solver => all_solvers(i)%s
         end do
         !check if name was found
         if (.not. associated(diag_solver)) call judft_error("Solver not found: "//name)
      else
         !defaults
         do i = first_real_solver, size(all_solvers)
            fit = all_solvers(i)%s%available
            !Check if solver fits the requirements
            if (use_gpu) fit = fit .and. all_solvers(i)%s%gpu
            if (parallel) then
               fit = fit .and. all_solvers(i)%s%parallel
            else
               fit = fit .and. all_solvers(i)%s%serial
            end if
            if (use_single_precision) fit = fit .and. all_solvers(i)%s%single_precision
            if (fit) exit
         end do
         if (i .le. size(all_solvers)) diag_solver => all_solvers(i)%s
         diag_solver%use_sp = use_single_precision
      end if

      !Check if a default tansformation must be selected as well
      if (.not. diag_solver%generalized .or. diag_solver%use_sp .and. &
          .not. associated(use_transform)) &
         then
         do i = first_real_solver, size(all_solvers)
            fit = all_solvers(i)%s%available .and. all_solvers(i)%s%transform
            !Check if solver fits the requirements
            if (use_gpu) fit = fit .and. all_solvers(i)%s%gpu
            if (parallel) then
               fit = fit .and. all_solvers(i)%s%parallel
            else
               fit = fit .and. all_solvers(i)%s%serial
            end if
            if (fit) exit
         end do
         if (i .le. size(all_solvers)) use_transform => all_solvers(i)%s
      end if

      ! check if selected solvers are OK
      generalized = .not. associated(use_transform)
      if (.not. associated(use_transform)) then
         ! we use a generalized solver
         if (.not. diag_solver%generalized) call judft_error("No generalized solver available")
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
         if (.not. associated(use_transform) .or. diag_solver%standard) &
            call judft_error("No standard solver available or missing transform")
         if (parallel) then
            if (.not. (diag_solver%parallel .and. use_transform%parallel)) &
               call judft_error("No parallel solver available for your problem", &
                                hint="You might have selected the wrong solver with the -diag option")
         else
            if (.not. (diag_solver%serial .and. use_transform%serial)) &
               call judft_error("No serial solver available for your problem", &
                                hint="You might have selected the wrong solver with the -diag option")
         end if
         if (use_gpu) then
            if (.not. (diag_solver%gpu .and. use_transform%gpu)) &
               call judft_warn("No GPU solver available for your problem", &
                               hint="You might have selected the wrong solver with the -diag option")
         end if
         diag_solver%transformer => use_transform
      end if
   end function select_solver

   function print_solver(parallel)
      logical, intent(IN):: parallel
      character(len=30):: print_solver
      class(t_solver), allocatable:: s

      s = select_solver(parallel)
      print_solver = s%name
      if (s%use_sp) print_solver = print_solver//"-sp"
      if (associated(s%transformer)) print_solver = s%transformer%name//'+'//print_solver
   end function

   subroutine list_solvers()
      integer:: i

      call init()

      write (*, '(a)') "Hints on choosing the diagonalization method (see the docu for more details):"
      write (*, '(a)') "  The `-diag` option takes a string as an argument that:"
      write (*, '(a)') "    Either simply specifies the solver used for the generalized eigenvalue problem (NAME_GEV):"
      write (*, '(a)') "      '-diag NAME_GEV'"
      write (*, '(a)') &
         "    Or the string combines (with '+') a solver of the standard problem (NAME_STD) with a transformation (NAME_TRANS):"
      write (*, '(a)') "      '-diag NAME_TRANS+NAME_STD'"
      write (*, '(a)') ""
      write (*, '(a)') "List of available solvers/transforms:"
      write (*, '(a)', ADVANCE="no") "  GEV-Solvers: "
      do i = 1, size(all_solvers)
         if (all_solvers(i)%s%available .and. all_solvers(i)%s%generalized) &
            write (*, '(a)', ADVANCE="no") all_solvers(i)%s%name
      end do
      write (*, *)
      write (*, '(a)', ADVANCE="no") "  STD-Solvers: "
      do i = 1, size(all_solvers)
         if (all_solvers(i)%s%available .and. all_solvers(i)%s%standard) then
            if (all_solvers(i)%s%single_precision) then
               write (*, '(a,a)', ADVANCE="no") all_solvers(i)%s%name, trim(all_solvers(i)%s%name)//"_sp"
            else
               write (*, '(a)', ADVANCE="no") all_solvers(i)%s%name
            end if
         end if
      end do

      write (*, *)
      write (*, '(a)', ADVANCE="no") "  Transforms : "
      do i = 1, size(all_solvers)
         if (all_solvers(i)%s%available .and. all_solvers(i)%s%transform) &
            write (*, '(a)', ADVANCE="no") all_solvers(i)%s%name
      end do
      write (*, *)
   end subroutine

end module m_available_solvers
