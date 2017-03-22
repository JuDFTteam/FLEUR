!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_setupMPI
  use m_juDFT
  IMPLICIT NONE

CONTAINS
  SUBROUTINE setupMPI(nkpt,mpi)
    USE m_types
    USE m_eigen_diag,ONLY:parallel_solver_available
    INTEGER,INTENT(in)           :: nkpt
    TYPE(t_mpi),INTENT(inout)    :: mpi

    IF (mpi%isize==1) THEN
       !give some info on available parallelisation
       CALL priv_dist_info(nkpt)
       mpi%n_start=1
       mpi%n_stride=1
       mpi%n_rank=0
       mpi%n_size=1
       mpi%n_groups=1
       mpi%sub_comm=mpi%mpi_comm
    END IF
#ifdef CPP_MPI
    !Distribute the work
    CALL priv_distribute_k(nkpt,mpi)

    !Now check is parallelization is possible
    IF (mpi%n_size>1.AND..NOT.parallel_solver_available()) &
         CALL juDFT_error("MPI parallelization failed",hint="You have to either compile FLEUR with a parallel diagonalization library (ELPA,SCALAPACK...) or you have to run such that the No of kpoints can be distributed on the PEs")       
#endif

    !generate the MPI communicators
    CALL priv_create_comm(nkpt,mpi)

  END SUBROUTINE setupMPI


  SUBROUTINE priv_distribute_k(nkpt,mpi)
    use m_types
    implicit none
    INTEGER,INTENT(in)      :: nkpt
    TYPE(t_mpi),INTENT(inout)    :: mpi

    !------------------------------------------------------------------------
    !
    ! Distribute the k-point / eigenvector  parallelisation so, that
    ! all pe's have aproximately equal load. Maximize for k-point 
    ! parallelisation. The naming conventions are as follows:
    !
    ! groups             1               2          (n_groups = 2) 
    !                 /     \         /     \
    ! k-points:      1       2       3       4      (nkpts = 4)
    !               /|\     /|\     /|\     /|\    
    ! irank        01 2   34 5   01 2   34 5    (isize = 6)
    !
    ! n_rank       01 2   01 2   01 2   01 2    (n_size = 3)
    !
    ! nrec         12 3   45 6   78 9  1011 12  ...rec. no. on eig-file
    !              * *     * *     * *     *  *
    !
    ! In the above example, 6 pe's should work on 4 k-points and distribute
    ! their load in a way, that 3 pe's work on each k-points, so 2 k-points
    ! are done in parellel (n_members=2) and there are 2 groups of k-points.
    ! n_rank and n_size are the equivalents of irank and isize. The former
    ! belong to the communicator SUB_COMM, the latter to MPI_COMM.
    !
    !          G.B. `99
    !
    !------------------------------------------------------------------------
    INTEGER:: n_members,n_size_min
    CHARACTER(len=20)::txt

    n_members = MIN(nkpt,mpi%isize)
    IF (judft_was_argument("-n_size_min")) THEN
       txt=judft_string_for_argument("-n_size_min")
       READ(txt,*) n_size_min
       WRITE(*,*) "Trying to use ",n_size_min," PE per kpt"
       n_members = MIN(n_members , CEILING(REAL(mpi%isize)/n_size_min) ) 
    ENDIF
    DO  
       IF ((MOD(mpi%isize,n_members) == 0).AND.(MOD(nkpt,n_members) == 0) ) EXIT
       n_members = n_members - 1
    ENDDO
    mpi%n_groups = nkpt/n_members
    mpi%n_size   = mpi%isize/n_members
    mpi%n_stride = n_members
    IF (mpi%irank == 0) THEN
       WRITE(*,*) 'k-points in parallel: ',n_members
       WRITE(*,*) "pe's per k-point:     ",mpi%n_size
       WRITE(*,*) '# of k-point loops:   ',mpi%n_groups
    ENDIF
  END SUBROUTINE priv_distribute_k

  SUBROUTINE priv_create_comm(nkpt,mpi)
    use m_types
    implicit none
    INTEGER,INTENT(in)      :: nkpt
    TYPE(t_mpi),INTENT(inout)    :: mpi
#ifdef CPP_MPI
    INTEGER:: n_members,n,i,ierr,sub_group,world_group
    INTEGER:: i_mygroup(mpi%n_size)


    n_members = nkpt/mpi%n_groups
    !
    ! now, we make the groups
    !
    mpi%n_start = MOD(mpi%irank,n_members) + 1
    !!      n_start = INT(irank/n_size) * n_size
    n = 0
    DO i = mpi%n_start,mpi%isize,n_members
       !!      DO i = n_start+1,n_start+n_size
       n = n+1
       i_mygroup(n) = i-1
    ENDDO


    CALL MPI_COMM_GROUP (mpi%MPI_COMM,WORLD_GROUP,ierr)
    CALL MPI_GROUP_INCL (WORLD_GROUP,mpi%n_size,i_mygroup,SUB_GROUP,ierr)
    CALL MPI_COMM_CREATE (mpi%MPI_COMM,SUB_GROUP,mpi%SUB_COMM,ierr)
    write (*,"(a,i0,100i4)") "MPI:",mpi%sub_comm,mpi%irank,mpi%n_groups,mpi%n_size,n,i_mygroup

    CALL MPI_COMM_RANK (mpi%SUB_COMM,mpi%n_rank,ierr)
#endif
  END SUBROUTINE priv_create_comm

  SUBROUTINE priv_dist_info(nkpt)
    USE m_eigen_diag,ONLY:parallel_solver_available
    IMPLICIT NONE
    INTEGER,INTENT(in)        :: nkpt

    INTEGER:: n,k_only,pe_k_only(nkpt)
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
  END SUBROUTINE priv_dist_info



END MODULE m_setupMPI


