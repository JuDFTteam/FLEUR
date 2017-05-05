      program test
      use elpa1
      integer:: ierr,mpi_subcom, myrowblacs, mycolblacs
      integer:: mpi_comm_rows,mpi_comm_cols,m,nb,mycolssca
      logical :: ok
      real :: bsca(10,10)
      ok=CHOLESKY_real (m,bsca,SIZE(bsca,1),nb,mycolssca,mpi_comm_rows,mpi_comm_cols,.false.)      
      end
