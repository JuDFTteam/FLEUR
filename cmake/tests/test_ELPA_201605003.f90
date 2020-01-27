      program test
      use elpa1
      integer:: ierr,mpi_subcom, myrowblacs, mycolblacs
      integer:: mpi_comm_rows,mpi_comm_cols,m,nb,mycolssca,myrowssca
      logical :: ok
      real :: bsca(10,10),asca(10,10),eigvec(10,10)
      ok= mult_at_b_real('U', 'L',m,m,bsca,myrowssca,asca,SIZE(asca,1),nb, mpi_comm_rows, mpi_comm_cols,eigvec,myrowssca)
      end
