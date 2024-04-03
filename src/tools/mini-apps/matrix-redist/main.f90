program matrixredist
    use m_types_mpimat
    use mpi
    use m_juDFT
    implicit none 

    type(t_mpimat):: cyclic1,cyclic2,rw 


    integer:: i,ierr,ii,j,irank,isize
    integer,PARAMETER:: N=10000,NN=10000
    
    !init MPI
    CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,i,ierr)

    !create the matrices
    call timestart("init matrices")
    call cyclic1%init(.true.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.true.)
    call cyclic2%init(.true.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.true.)
    call rw%init(.true.,N,N,mpi_subcom=MPI_COMM_WORLD,l_2d=.false.)
    !initialize the matrix with pattern
    call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,isize,ierr)
    !i=0
    !DO ii=irank+1,n,isize
    !    i=i+1
    !    DO j=1,N
    !        rw%data_r(j,i)=j*NN+ii
    !    ENDDO
    !ENDDO    
    call random_number(rw%data_r)
    call timestop("init matrices")
    !Do the redistribution
    call timestart("redist matrix1")
    use_pdgemr2d=.true.
    call cyclic1%copy(rw,1,1)
    call timestop("redist matrix1")
    !Same without ScaLAPACK
    use_pdgemr2d=.false.
    call timestart("redist matrix2")
    call cyclic2%copy(rw,1,1)
    call timestop("redist matrix2")

    print *,"Error:",maxval(abs(cyclic1%data_r-cyclic2%data_r))>1E-10

    
    call judft_end("Done",0)
    !Done
    contains
    function to_coord(x)
        real :: x
        integer::to_coord(2)

        to_coord(1)=x/NN
        to_coord(2)=x-to_coord(1)*NN
    END function
end    
