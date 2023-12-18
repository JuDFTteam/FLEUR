!--------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_test_performance
    !$ use omp_lib
    implicit none
    private
    public test_performance
contains
    subroutine check_multithreading(dim,time)
        INTEGER,INTENT(IN):: dim
        integer :: n,n_threads 
        real,intent(out):: time(:)
        real:: vector(dim*dim),sum

#ifdef _OPENMP
        call random_number(vector)
        !Do a small OMP loop to make sure OMP-lib is initialized
        sum=0
        !$OMP PARALLEL DO private(n)shared(vector)default(none)reduction(+:sum)
        DO n=1,size(vector)
            sum=sum+vector(n)
        ENDDO
        !$OMP END PARALLEL DO

        DO n_threads=1,size(time)
            time(n_threads)=omp_get_wtime()
            sum=0.0
            !$OMP PARALLEL DO num_threads(n_threads) private(n)shared(vector)default(none)reduction(+:sum)
            DO n=1,size(vector)
                sum=sum+sin(vector(n))
            ENDDO
            !$OMP end parallel do
            time(n_threads)=omp_get_wtime()-time(n_threads)
        ENDDO        
#endif
    end SUBROUTINE      
    
    subroutine check_blas_multithreading(N,time)
        integer,intent(in):: N
        real :: a(N,N),c(N,N)
        real,intent(out):: time(:)
        integer :: n_threads

#ifdef _OPENMP        
        call random_number(a)
        !Call a dgemm to initialize the BLAS
        call dgemm("N","N",N,N,N,1.0,a,N,a,N,0.0,c,N)

        DO n_threads=1,size(time)
            call omp_set_num_threads(n_threads)
            time(n_threads)=omp_get_wtime()
            call dgemm("N","N",N,N,N,1.0,a,N,a,N,0.0,c,N)
            time(n_threads)=omp_get_wtime()-time(n_threads)
        ENDDO        
#endif
    end SUBROUTINE      

    subroutine check_diag_multithreading(N,time)
        use m_eigen_diag
        use m_types_mat
        use m_types_mpimat
#ifdef CPP_MPI
        use mpi 
#endif        
        REAL,INTENT(OUT) :: time(:)
        INTEGER,INTENT(IN)  :: N
      
        INTEGER          :: ne,n_threads,solver=0,i,j
        CLASS(t_mat), ALLOCATABLE :: ev         ! eigenvectors
        REAL :: eig(N) 
        class(t_mat),allocatable :: hmat,smat,h,s
#ifdef _OPENMP        
#ifdef CPP_MPI     
        integer :: ierr,isize

        call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

        if (isize>1) THEN
            allocate(t_mpimat:: hmat,smat,h,s)
            CALL hmat%init(.true., N,N, MPI_COMM_WORLD, .true.)
        ELSE 
            allocate(t_mat:: hmat,smat,h,s)
            CALL hmat%init(.true., N,N)
        ENDIF           
#else
        allocate(t_mat:: hmat,smat,h,s)
        CALL hmat%init(.true., N,N)
#endif        
        CALL smat%init(hmat)
        CALL h%init(hmat)
        CALL s%init(hmat)

        call random_number(hmat%data_r(:,:))
        !make smat positive (semi)definite
        call hmat%multiply(hmat,smat,"T","N")
       

        ne=0.1*N
     

        ! 1 thread
        DO n_threads=1,size(time)
            h%data_r=hmat%data_r
            s%data_r=smat%data_r
            call omp_set_num_threads(n_threads)
            time(n_threads)=omp_get_wtime()
            call eigen_diag(solver,h,s,ne,eig,ev)
            time(n_threads)=omp_get_wtime()-time(n_threads)
        ENDDO    
        
#endif
    end subroutine

    subroutine test_performance()
        use m_types_lapw
        INTEGER :: dim
        real,allocatable::time(:,:)
        INTEGER :: n

        dim=lapw_dim_nvd
#ifdef _OPENMP
        if (omp_get_max_threads()==1) THEN
            write(*,*) "Only a single OMP thread available, not OMP scaling tested"
            RETURN
        ENDIF    
        write(*,*) "Testing OMP scaling"
        allocate(time(omp_get_max_threads(),3))
        write(*,*) "Testing BLAS scaling"
        call check_multithreading(dim,time(:,1))
        write(*,*) "Testing Eigenvalue solver scaling"
        call check_blas_multithreading(dim,time(:,2))
        call check_diag_multithreading(dim,time(:,3))
        write(*,*) "THREADS|      OpenMP         |        BLAS         |      Eigenvalue"
        write(*,*) "       |   time      scaling |   time      scaling |   time      scaling "
        time=time*1000
        DO n=1,size(time,1)
            write(*,"(i3,5x,3('|',f8.2,3x,f8.2,2x))") n,time(n,1),time(1,1)/time(n,1),time(n,2),time(1,2)/time(n,2),time(n,3),time(1,3)/time(n,3)
        ENDDO
        if (any(time(2,:)/time(1,:)<1.5)) then
            write(*,*) "WARNING"
            write(*,*) "No proper OMP scaling detected"
            write(*,*) "You should check if you:"
            write(*,*) "    - use a multithreaded solver/BLAS"
            write(*,*) "    - use a proper distribution of threads to your hardware"
        else
            write(*,*) "OpenMP scaling looks OK"
        endif  
#else
        write(*,*) "NO OMP-Version"
#endif        
    end subroutine    
end module
