program matrixredist
    use m_types_mpimat
    use mpi
    implicit none 

    types(t_mpimat):: cyclic,rw 


    integer:: i,ierr
    integer,PARAMETER:: N=1000
    
    !init MPI
    CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,i,ierr)

    !create the matrices
    call timestart("init matrices")
    call cyclic%init(.false.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.true.)
    call rw%init(.false.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.false.)
    call timestop("init matrices")
    !Do the redistribution
    call timestart("redist matrix")
    call cyclic.copy(rw,1,1)
    call timestop("redist matrix")
    !Done
    call judft_end("Done")
end    