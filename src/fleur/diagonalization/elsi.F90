!-------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_elsi
   use m_types_solver

   type, extends(t_solver)::t_solver_elsi
   contains
      procedure        :: solve_gev => elsi_GEV
   end type
   public :: get_solver_elsi

contains

   function get_solver_elsi() result(solver)
      type(t_solver_elsi), pointer::solver
      allocate (solver)
      solver%name = "elsi"
#ifdef CPP_ELSI
      solver%available = .true.
#else
      solver%available = .false.
#endif
      solver%parallel = .true.
      solver%serial = .true.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .false.
   end function

   subroutine elsi_gev(self, hmat, smat, ne, eig, zmat)

      use m_juDFT
      use m_types_mpimat
      use m_types_mat
#ifdef CPP_ELSI
      use elsi
      use mpi
#endif
      implicit none
      class(t_solver_elsi)           ::self
      class(t_mat), intent(INOUT)    :: hmat, smat
      class(t_mat), allocatable, intent(OUT)::zmat
      real, intent(out)              :: eig(:)
      integer, intent(INOUT)         :: ne

      integer, parameter             :: solver = 1 !Use ELPA

#ifdef CPP_ELSI
      !...  Local variables
      !
      integer           :: blk, nev, myid, np, ierr, i
      integer, parameter :: BLACS_DENSE = 0
      type(elsi_handle) :: eh
      class(t_mat), allocatable :: evec
      real, allocatable         :: eig_tmp(:)

      call timestart("ELSI")

      call hmat%u2l()
      call smat%u2l()

      select type (hmat)
      type IS (t_mpimat)
         select type (smat)
         type IS (t_mpimat)
            if (hmat%blacsdata%blacs_desc(5) == hmat%blacsdata%blacs_desc(6)) then
               blk = hmat%blacsdata%blacs_desc(5)
            else
               call judft_error("BUG: in ELSI the row/column blocksize must be equal")
            end if
            call elsi_init(eh, solver, 1, BLACS_DENSE, hmat%global_size1, 1.0*ne, ne)
            call elsi_set_mpi(eh, hmat%blacsdata%mpi_com)
            call elsi_set_blacs(eh, hmat%blacsdata%blacs_desc(2), blk)
            allocate (t_mpimat::evec)
            call evec%init(hmat)
            allocate (eig_tmp(hmat%global_size1))
         class DEFAULT
            call judft_error("BUG: Inconsistent matrixes in ELSI call")
         end select
      type IS (t_mat)
         select type (smat)
         type IS (t_mat)
            call elsi_init(eh, solver, 0, BLACS_DENSE, hmat%matsize1, 1.*ne, ne)
            allocate (t_mat::evec)
            call evec%init(hmat)
            allocate (eig_tmp(hmat%matsize1))

         class DEFAULT
            call judft_error("BUG: Inconsistent matrixes in ELSI call")
         end select
      end select
      call timestart("elsi_ev")
#ifdef _OPENACC
      call elsi_set_elpa_gpu(eh, 1)
#endif
      !Now perform diagonalization
      call elsi_set_output(eh, 3)
      call elsi_set_output_unit(eh, 7)
      call elsi_set_illcond_check(eh, 0)
      call elsi_reinit(eh)
      if (hmat%l_real) then
         call elsi_ev_real(eh, hmat%data_r, smat%data_r, eig_tmp, evec%data_r)
      else
         call elsi_ev_complex(eh, hmat%data_c, smat%data_c, eig_tmp, evec%data_c)
      end if
      call timestop("elsi_ev")
      call timestart("redistribute")
      !Copy data into correct data structures
      eig = eig_tmp(:size(eig))
      select type (evec)
      type IS (t_mat)
         allocate (t_mat::zmat)
         call zamt%init(evec%l_real, evec%matsize1, ne)
         if (evec%l_real) then
            zmat%data_r = evec%data_r(:, :ne)
         else
            zmat%data_c = evec%data_c(:, :ne)
         end if
      type IS (t_mpimat)
         allocate (t_mpimat::zmat)
         call zmat%init(evec%l_real, evec%global_size1, evec%global_size1, evec%blacsdata%mpi_com, .false.)
         call zmat%copy(evec, 1, 1)
         !determine ev assigned to this rank
         nev = ne
         ne = 0
         call MPI_COMM_RANK(evec%blacsdata%mpi_com, myid, ierr)
         call MPI_COMM_SIZE(evec%blacsdata%mpi_com, np, ierr)
         do i = myid + 1, nev, np
            ne = ne + 1
            !   eig(ne)=eigenvalues(i)
         end do

      end select
      call evec%free()
      call timestop("redistribute")
      call elsi_finalize(eh)
      call timestop("ELSI")

#endif
   end subroutine elsi_gev
end module m_elsi
