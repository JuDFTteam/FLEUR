subroutine fortran_check_mpi() bind(C, name="fortran_check_mpi")
    use MPI 
    logical :: flag
    INTEGER :: mySTAT(MPI_STATUS_SIZE),ierr
    call mpi_iprobe(MPI_ANY_SOURCE,1336, MPI_COMM_WORLD,flag,mystat,ierr)
  end subroutine