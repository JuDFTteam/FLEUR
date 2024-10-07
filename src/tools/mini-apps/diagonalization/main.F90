program testdiag 
    use m_juDFT
    use m_types_mat
    use m_eigen_diag
#ifdef CPP_MPI
    use m_types_mpimat
#endif        
    implicit none
    
    character(len=10)         :: arg
    class(t_mat),allocatable  :: H,S,evec
    real,allocatable    :: eig(:)
    integer             :: n,ne
    integer             :: solver=0

#ifdef CPP_MPI
    block
        integer:: i,ierr
        CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,i,ierr)
    end block
#endif    
    !determine matrix size
    IF (judft_was_argument("-N")) THEN
        arg=juDFT_string_for_argument("-N")
        read(arg,*) N
    ELSE    
        N=10000
    ENDIF    

    ne=N/10 !No of eigenvectors to calculate
    allocate(eig(ne))

    call timestart("Matrix setup")
    call setup_matrices(N,H,S)
    call timestop("Matrix setup")
    call timestart("Diagonalization")
    call eigen_diag(solver,h,s,ne,eig,evec)
    call timestop("Diagonalization")

    call judft_end("done",0)
    CONTAINS
        subroutine setup_matrices(n,H,S)
        integer,intent(in)      :: n
        class(t_mat),ALLOCATABLE,INTENT(OUT):: H,S
        class(t_mat),Allocatable            :: tmp
        real,allocatable::rand(:,:,:)

        integer:: isize,ierr
#ifdef CPP_MPI
        call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)
        if (isize>1)then
            ALLOCATE (t_mpimat::s, h ,tmp)
            CALL s%init(.false., N,N, MPI_COMM_WORLD, .true.)
        else
            ALLOCATE (t_mat::s, h ,tmp)
            CALL s%init(.false., N,N)
        endif    
#else
        ALLOCATE (t_mat::s, h,tmp)
        CALL s%init(.false., N,N)
#endif    
        CALL h%init(s)
        call tmp%init(s)

        !Now fill with random numbers
        allocate(rand(2,size(S%data_c,1),size(s%data_c,2)))
        call random_seed()
        call random_number(rand)
        h%data_c=cmplx(rand(1,:,:),rand(2,:,:))
        tmp%data_c=cmplx(rand(1,:,:),rand(2,:,:))
        !make s positive definite (might actually only be semi-definite)
        call tmp%multiply(h,s,"N","C")

        !Make a random H
        call random_number(rand)
        h%data_c=cmplx(rand(1,:,:),rand(2,:,:))
        end subroutine
end