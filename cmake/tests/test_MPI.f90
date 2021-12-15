      program test
      USE mpi
      integer:: rank,ierr
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      end program
