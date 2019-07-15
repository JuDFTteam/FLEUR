!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_setupMPI
  use m_juDFT
  IMPLICIT NONE

CONTAINS
  SUBROUTINE setupMPI(nkpt,neigd,mpi)
!$  use omp_lib
    USE m_types  
    USE m_available_solvers,ONLY:parallel_solver_available
    INTEGER,INTENT(in)           :: nkpt,neigd
    TYPE(t_mpi),INTENT(inout)    :: mpi

    INTEGER :: omp=-1,i

    !$ omp=omp_get_max_threads()
    if (mpi%irank==0) THEN
       !print INFO on parallelization
       WRITE(*,*) "--------------------------------------------------------"
#ifdef CPP_MPI
       write(*,*) "Number of MPI-tasks:  ",mpi%isize
       CALL add_usage_data("MPI-PE",mpi%isize)     
#else
       CALL add_usage_data("MPI-PE",1)     
#endif
       IF (omp==-1) THEN
          write(*,*) "No OpenMP version of FLEUR."
          CALL add_usage_data("OMP",0)
       ELSE
          WRITE(*,*) "Number of OMP-threads:",omp
          CALL add_usage_data("OMP",omp)
       ENDIF
    endif
    IF (mpi%isize==1) THEN
       !give some info on available parallelisation
       CALL priv_dist_info(nkpt)
       mpi%n_rank=0
       mpi%n_size=1
       mpi%sub_comm=mpi%mpi_comm
       IF (ALLOCATED(mpi%k_list)) DEALLOCATE(mpi%k_List,mpi%ev_list)
       ALLOCATE(mpi%k_list(nkpt))
       mpi%k_list=[(i,i=1,nkpt)]
       ALLOCATE(mpi%ev_list(neigd))
       mpi%ev_list=[(i,i=1,neigd)]
    END IF
#ifdef CPP_MPI
    !Distribute the work
    CALL priv_distribute_k(nkpt,mpi)

    !Now check if parallelization is possible
    IF (mpi%n_size>1.AND..NOT.parallel_solver_available()) &
         CALL juDFT_error("MPI parallelization failed",hint="You have to either compile FLEUR with a parallel diagonalization library (ELPA,SCALAPACK...) or you have to run such that the No of kpoints can be distributed on the PEs")       
#endif

    !generate the MPI communicators
    CALL priv_create_comm(nkpt,neigd,mpi)
    if (mpi%irank==0) WRITE(*,*) "--------------------------------------------------------"

  END SUBROUTINE setupMPI


  SUBROUTINE priv_distribute_k(nkpt,mpi)
    use m_types
    implicit none
    INTEGER,INTENT(in)      :: nkpt
    TYPE(t_mpi),INTENT(inout)    :: mpi

    !-------------------------------------------------------------------------------------------
    !
    ! Distribute the k-point / eigenvector  parallelisation so, that
    ! all pe's have aproximately equal load. Maximize for k-point 
    ! parallelisation. The naming conventions are as follows:
    !
    ! groups             1               2               3             4      (n_groups = 4) 
    !                 /     \         /     \          /   \         /   \
    ! k-points:      1       2       3       4       5       6      7     8     (nkpts = 8)
    !               /|\     /|\     /|\     /|\     /|\     /|\    /|\   /|\
    ! irank        0 1 2   3 4 5   1 2 3   4 5 6   0 1 2   3 4 5  1 2 3  4 5 6  (mpi%isize = 6)
    !
    ! n_rank       0 1 2   0 1 2   0 1 2   0 1 2   0 1 2   0 1 2  0 1 2  0 1 2  (mpi%n_size = 3)
    !
    !
    ! In the above example, 6 pe's should work on 8 k-points and distribute
    ! their load in a way, that 3 pe's work on each k-points, so 2 k-points
    ! are done in parellel (n_members=2) and each processor runs a loop over
    ! 4 k-points (mpi%n_groups = 4).
    ! n_rank and n_size are the equivalents of irank and isize. The former
    ! belong to the communicator SUB_COMM, the latter to MPI_COMM.
    !
    !          G.B. `99
    !
    !-------------------------------------------------------------------------------------------
    INTEGER:: n_members,n_size_min,nk
    CHARACTER(len=1000)::txt

    n_members = MIN(nkpt,mpi%isize)
    IF (judft_was_argument("-n_min_size")) THEN
       txt=judft_string_for_argument("-n_min_size")
       READ(txt,*) n_size_min
       WRITE(*,*) "Trying to use ",n_size_min," PE per kpt"
       n_members = MIN(n_members , CEILING(REAL(mpi%isize)/n_size_min) ) 
    ENDIF
    DO  
       IF ((MOD(mpi%isize,n_members) == 0).AND.(MOD(nkpt,n_members) == 0) ) EXIT
       n_members = n_members - 1
    ENDDO
    ALLOCATE(mpi%k_list(nkpt/n_members))
    mpi%k_list=[(nk, nk=nkpt/n_members,nkpt,n_members )]

    !mpi%n_groups = nkpt/n_members
    mpi%n_size   = mpi%isize/n_members
    !mpi%n_stride = n_members
    IF (mpi%irank == 0) THEN
       WRITE(*,*) 'k-points in parallel: ',n_members
       WRITE(*,*) "pe's per k-point:     ",mpi%n_size
       WRITE(*,*) '# of k-point loops:   ',nkpt/n_members
    ENDIF
  END SUBROUTINE priv_distribute_k

  SUBROUTINE priv_create_comm(nkpt,neigd,mpi)
    use m_types
    implicit none
    INTEGER,INTENT(in)      :: nkpt,neigd
    TYPE(t_mpi),INTENT(inout)    :: mpi
#ifdef CPP_MPI
    INTEGER :: n_members,n,i,ierr,sub_group,world_group,n_start
    INTEGER :: i_mygroup(mpi%n_size)
    LOGICAL :: compact ! Deside how to distribute k-points

    compact = .true.
    n_members = mpi%isize/mpi%n_size
    
    ! now, we make the groups
    
    
    IF (compact) THEN

        ! This will distribute sub ranks in a compact manner.
        ! For example, if nkpt = 8 and mpi%isize = 6:
        
        !  -----------------------------------
        ! |  0  |  1  |  2  |  3  |  4  |  5  |    mpi%irank
        !  -----------------------------------
        ! |  0  |  1  |  3  |  0  |  1  |  2  |    mpi%n_rank
        !  -----------------------------------
        ! |        1        |        2        |    k - points
        ! |        3        |        4        |
        ! |        5        |        6        |
        ! |        7        |        8        |
        !  -----------------------------------
    
        n_start = INT(mpi%irank/mpi%n_size) + 1
        i_mygroup(1) = (n_start-1) * mpi%n_size
        do i = 2, mpi%n_size
           i_mygroup(i) = i_mygroup(i-1) + 1
        enddo

    ELSE

        ! This will distribute sub ranks in a spread manner.
        ! For example, if nkpt = 8 and mpi%isize = 6:
    
        !  -----------------------------------
        ! |  0  |  1  |  2  |  3  |  4  |  5  |    mpi%irank
        !  -----------------------------------
        ! |  0  |  1  |  3  |  0  |  1  |  2  |    mpi%n_rank
        !  -----------------------------------
        ! |  1  |  2  |  1  |  2  |  1  |  2  |    k - points
        ! |  3  |  4  |  3  |  4  |  3  |  4  |
        ! |  5  |  6  |  5  |  6  |  5  |  6  |  
        ! |  7  |  8  |  7  |  8  |  7  |  8  |
        !  -----------------------------------

        n_start = MOD(mpi%irank,n_members) + 1
        !!      n_start = INT(irank/n_size) * n_size
        n = 0
        DO i = n_start,mpi%isize,n_members
        !!      DO i = n_start+1,n_start+n_size
           n = n+1
           i_mygroup(n) = i-1
        ENDDO

    ENDIF ! compact

    CALL MPI_COMM_GROUP (mpi%MPI_COMM,WORLD_GROUP,ierr)
    CALL MPI_GROUP_INCL (WORLD_GROUP,mpi%n_size,i_mygroup,SUB_GROUP,ierr)
    CALL MPI_COMM_CREATE (mpi%MPI_COMM,SUB_GROUP,mpi%SUB_COMM,ierr)
    !write (*,"(a,i0,100i4)") "MPI:",mpi%sub_comm,mpi%irank,mpi%n_groups,mpi%n_size,n,i_mygroup

    CALL MPI_COMM_RANK (mpi%SUB_COMM,mpi%n_rank,ierr)
    ALLOCATE(mpi%ev_list(neigd/mpi%n_size+1))
    mpi%ev_list=[(i,i=mpi%irank+1,neigd,mpi%n_size)]
    
#endif
  END SUBROUTINE priv_create_comm

  SUBROUTINE priv_dist_info(nkpt)
    USE m_available_solvers,ONLY:parallel_solver_available
    IMPLICIT NONE
    INTEGER,INTENT(in)           :: nkpt

    INTEGER:: n,k_only,pe_k_only(nkpt)
    
#ifdef CPP_MPI
    !Create a list of PE that will lead to k-point parallelization only
    k_only=0
    DO n=1,nkpt
       IF (MOD(nkpt,n)==0) THEN
          k_only=k_only+1
          pe_k_only(k_only)=n
       ENDIF
    END DO
    WRITE(*,*) "Most efficient MPI parallelization for:"
    WRITE(*,*) pe_k_only(:k_only)
    !check if eigenvalue parallelization is possible
    IF (parallel_solver_available()) WRITE(*,*) "Additional eigenvalue parallelization possible"
#endif
  END SUBROUTINE priv_dist_info



END MODULE m_setupMPI


