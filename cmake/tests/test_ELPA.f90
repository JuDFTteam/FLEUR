      program test  !Basically test if the elpa module is available and elpa_init function exists
      use elpa
      integer:: ierr,mpi_subcom, myrowblacs, mycolblacs
      integer:: mpi_comm_rows,mpi_comm_cols
      ierr=elpa_init(20171201)
      end program
