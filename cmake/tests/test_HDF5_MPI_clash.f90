!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Some parallel HDF5 builds export the MPI_F08 interface through the hdf5 module.
! Together with 'use mpi' this clashes as soon as MPI_IN_PLACE (and similar
! constants) is referenced, e.g. NVHPC: "mpi_fortran_in_place cannot be a common
! block and a subprogram". The same scope layout as in FLEUR is used here.
      program test
      implicit none
      call reduce()
      contains
      subroutine reduce()
      use hdf5
      use mpi
      integer :: ierr
      real    :: x(2)
      x = 1.0
      call MPI_INIT(ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, x, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)
      end subroutine reduce
      end program
