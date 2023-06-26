recursive subroutine fortran_check_mpi() bind(C, name="fortran_check_mpi")
    use MPI
    logical :: flag
    INTEGER :: ierr
    call mpi_iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD,flag,MPI_STATUS_IGNORE,ierr)
  end subroutine
