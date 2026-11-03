!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
program diag_test
   use m_judft
   use m_types_mat
   use m_types_mpimat
   use m_eigen_diag
   use m_io_matrix
   use m_hdf_tools
#ifdef CPP_MPI
   use mpi
#endif
   implicit none

   integer            :: matsize, ne, mode, fid
   integer            :: err, isize
   character(len=50)  :: filename
   class(t_mat), allocatable :: hmat, smat, ev
   real, allocatable   :: eig(:)
   logical            :: l_exist, l_real
   real               :: t1, t2
#ifdef CPP_MPI
   call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, isize, err)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, err)
   if (isize > 1) then
      allocate (t_mpimat::hmat)
      allocate (t_mpimat::smat)
      select type (hmat)
      type is (t_mpimat)
         call hmat%set_mpi_com(MPI_COMM_WORLD)
      end select
      select type (smat)
      type is (t_mpimat)
         call smat%set_mpi_com(MPI_COMM_WORLD)
      end select
   end if
#endif
   if (.not. allocated(hmat)) then
      allocate (t_mat::hmat)
      allocate (t_mat::smat)
   end if

   ! get filename
   filename = judft_string_for_argument("-file")
   inquire (file=trim(filename)//".hdf", exist=l_exist)
   if (.not. l_exist) call judft_error("File specified does not exist")
   call hdf_init()

   !l_real,matsize is actually only needed if file is created
   fid = open_matrix(l_real, matsize, 2, 2, trim(filename))

   call read_matrix(hmat, 1, fid)
   call read_matrix(smat, 2, fid)
   select type (hmat)
   type is (t_mpimat)
      select type (smat)
      type is (t_mpimat)
         call smat%share_blacsgrid(hmat)!make sure we use same blacs-grids
      end select
      ne = 0.15*hmat%global_size1
      allocate (eig(hmat%global_size1))
   class default
      ne = 0.15*hmat%matsize1
      allocate (eig(hmat%matsize1))
   end select

   call cpu_time(t1)
   mode = 0
   call eigen_diag(mode, hmat, smat, ne, eig, ev)
   call cpu_time(t2)
   print *, "No of eigenvalues:", ne
   print *, eig(:ne)

   print *, "Time used:", t2 - t1

   call close_matrix(fid)

end program diag_test

