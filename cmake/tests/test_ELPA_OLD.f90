      program test
      use elpa1
      integer:: ierr,mpi_subcom, myrowblacs, mycolblacs
      integer:: mpi_comm_rows,mpi_comm_cols
      call get_elpa_row_col_comms(mpi_subcom, myrowblacs, mycolblacs,        &
     &     mpi_comm_rows, mpi_comm_cols)
      end
