program matrixredist
    use m_types_mpimat
    use mpi
    use m_juDFT
    implicit none 

    type(t_mpimat):: cyclic1,cyclic2,rw 


    integer:: i,ierr
    integer,PARAMETER:: N=40000
    
    !init MPI
    CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,i,ierr)

    !create the matrices
    call timestart("init matrices")
    call cyclic1%init(.true.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.true.)
    call cyclic2%init(.true.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.true.)
    call rw%init(.true.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.false.)
    call random_number(rw%data_r)
    call timestop("init matrices")
    !Do the redistribution
    call timestart("redist matrix1")
    call cyclic1%copy(rw,1,1)
    call timestop("redist matrix1")
    !Same without ScaLAPACK
    use_pdgemr2d=.false.
    call timestart("redist matrix2")
    call cyclic2%copy(rw,1,1)
    call timestop("redist matrix2")

    print *,maxval(abs(cyclic1%data_r-cyclic2%data_r))

    print *,"Done"
    call judft_end("Done",0)
    !Done
end    
