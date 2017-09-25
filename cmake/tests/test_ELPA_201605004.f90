      program test
      use elpa1
      integer:: ierr,mpi_subcom, myrowblacs, mycolblacs
      integer:: mpi_comm_rows,mpi_comm_cols,m,nb,myrowsssca,mycolssca
      logical :: ok
      real :: bsca(10,10),asca(10,10),eigvec(10)
      ok= elpa_mult_at_b_real('U', 'L',m, m,bsca,myrowssca,mycolssca,asca,SIZE(asca,1),SIZE(asca,2),nb,&
         mpi_comm_rows, mpi_comm_cols,eigvec,myrowssca,mycolssca)
      !ok=CHOLESKY_real (m,bsca,SIZE(bsca,1),nb,mycolssca,mpi_comm_rows,mpi_comm_cols,.false.)      
      end
