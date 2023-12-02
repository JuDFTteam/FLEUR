!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_setupMPI
#ifdef CPP_MPI
  use mpi
#endif
  use m_juDFT
  IMPLICIT NONE

CONTAINS
  SUBROUTINE setupMPI(nkpt,neigd,nbasfcn,fmpi)
!$  use omp_lib

    use m_omp_checker
    USE m_types
    USE m_available_solvers,ONLY:parallel_solver_available
    INTEGER,INTENT(in)           :: nkpt,neigd,nbasfcn
    TYPE(t_mpi),INTENT(inout)    :: fmpi

    INTEGER :: omp=-1,i,isize,localrank,gpus,ii, me, nk
    logical :: finished
#ifdef CPP_MPI
    CALL juDFT_COMM_SPLIT_TYPE(fmpi%mpi_comm,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,fmpi%mpi_comm_same_node)
#endif
    call omp_checker()
    !$ omp=omp_get_max_threads()
    if (fmpi%irank==0) THEN
       !print INFO on parallelization
       WRITE(*,*) "--------------------------------------------------------"
#ifdef CPP_MPI
       write(*,*) "Number of MPI-tasks:  ",fmpi%isize
       CALL MPI_COMM_SIZE(fmpi%mpi_comm_same_node,isize,i)
       write(*,*) "Number of PE/node  :  ",isize
       CALL add_usage_data("MPI-PE",fmpi%isize)
#else
       CALL add_usage_data("MPI-PE",1)
#endif
       IF (omp==-1) THEN
          write(*,*) "No OpenMP version of FLEUR."
          CALL add_usage_data("OMP",0)
       ELSE
          WRITE(*,*) "Number of OMP-threads:",omp
          IF(omp.EQ.1.AND.fmpi%isize.GE.6.AND.&
             ABS(NINT(REAL(nkpt)/REAL(fmpi%isize))*fmpi%isize-nkpt).GT.1.0e-7) THEN
             WRITE(*,*) ''
             WRITE(*,*) '========================================'
             WRITE(*,*) 'WARNING:'
             WRITE(*,*) 'You are making use of multiple MPI processes but no OpenMP parallelization.'
             WRITE(*,*) 'The chosen parallelization scheme is also no pure k-point parallelization.'
             WRITE(*,*) 'The performance of your calculation may benefit by also employing some'
             WRITE(*,*) 'OpenMP parallelization.'
             WRITE(*,*) ''
             WRITE(*,*) 'Fleur will proceed with the calculation.'
             WRITE(*,*) '========================================'
             WRITE(*,*) ''
          END IF

          CALL add_usage_data("OMP",omp)
       ENDIF
    endif
    call priv_distribute_gpu(fmpi)
    IF (fmpi%isize==1) THEN
       !give some info on available parallelisation
       CALL priv_dist_info(nkpt)
       fmpi%n_rank=0
       fmpi%n_size=1
       fmpi%sub_comm=fmpi%mpi_comm
       fmpi%diag_sub_comm=fmpi%mpi_comm
       IF (ALLOCATED(fmpi%k_list)) DEALLOCATE(fmpi%k_List)
       IF (ALLOCATED(fmpi%ev_list)) DEALLOCATE(fmpi%ev_list)
       ALLOCATE(fmpi%ev_list(neigd))
       ALLOCATE(fmpi%k_list(nkpt))
       fmpi%k_list=[(i,i=1,nkpt)]
       fmpi%coulomb_owner=[(0,i=1,nkpt)]
       fmpi%ev_list=[(i,i=1,neigd)]
       WRITE(*,*) "--------------------------------------------------------"
       RETURN
    END IF
#ifdef CPP_MPI
    !Distribute the work
    CALL priv_distribute_k(nkpt,nbasfcn,fmpi)
    !generate the MPI communicators
    CALL priv_create_comm(nkpt,neigd,fmpi)
    !Now check if parallelization is possible
    IF (fmpi%n_size>1) THEN
      if (judft_was_argument("-serial_diag")) THEN
        call priv_redist_for_diag(fmpi)
      else if (.NOT.parallel_solver_available()) then
        call juDFT_error("MPI parallelization failed",hint="You have to either compile FLEUR with a parallel diagonalization library (ELPA,SCALAPACK...) or you have to run such that the No of kpoints can be distributed on the PEs")
      endif
    endif
#endif


    ALLOCATE(fmpi%k_list(SIZE([(i, i=INT(fmpi%irank/fmpi%n_size)+1,nkpt,fmpi%isize/fmpi%n_size )])))
    ! this corresponds to the compact = .true. switch in priv_create_comm
    fmpi%k_list=[(i, i=INT(fmpi%irank/fmpi%n_size)+1,nkpt,fmpi%isize/fmpi%n_size )]

    ! create an array with the owners of the correct coulomb matrix
    allocate(fmpi%coulomb_owner(nkpt), source=-1)
    do nk =1,nkpt
      me = 0
      finished = .False.
      do while(.not. finished)
         if(any(nk == [(i, i=INT(me/fmpi%n_size)+1,nkpt,fmpi%isize/fmpi%n_size)] )) then
            fmpi%coulomb_owner(nk) = me
            finished = .True.
         endif
         me = me + 1
         if(me > fmpi%isize .and. .not. finished) then
            call judft_error("somehow i cant lokate this k-point")
         endif
      enddo
   enddo


    call fmpi%set_errhandler()
    if (fmpi%irank==0) WRITE(*,*) "--------------------------------------------------------"

  END SUBROUTINE setupMPI


  SUBROUTINE priv_distribute_k(nkpt,nbasfcn,fmpi)
    use m_types
    implicit none
    INTEGER,INTENT(in)      :: nkpt, nbasfcn
    TYPE(t_mpi),INTENT(inout)    :: fmpi

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
    ! irank        0 1 2   3 4 5   1 2 3   4 5 6   0 1 2   3 4 5  1 2 3  4 5 6  (fmpi%isize = 6)
    !
    ! n_rank       0 1 2   0 1 2   0 1 2   0 1 2   0 1 2   0 1 2  0 1 2  0 1 2  (fmpi%n_size = 3)
    !
    !
    ! In the above example, 6 pe's should work on 8 k-points and distribute
    ! their load in a way, that 3 pe's work on each k-points, so 2 k-points
    ! are done in parellel (n_members=2) and each processor runs a loop over
    ! 4 k-points (fmpi%n_groups = 4).
    ! n_rank and n_size are the equivalents of irank and isize. The former
    ! belong to the communicator SUB_COMM, the latter to MPI_COMM.
    !
    !          G.B. `99
    !
    !-------------------------------------------------------------------------------------------
    INTEGER:: n_members,n_size_min,nk
    CHARACTER(len=1000)::txt

    IF (judft_was_argument("-pe_per_kpt")) THEN
      txt=judft_string_for_argument("-pe_per_kpt")
      READ(txt,*) fmpi%n_size
      WRITE(*,*) "Using exactly ",fmpi%n_size," PE per kpt"
      if (mod(fmpi%isize,fmpi%n_size).ne.0) call judft_error("Parallelization error",&
           hint="If you specify -pe_per_kpt, the total number of processes needs to be a multiple.")   
      n_members=fmpi%isize/fmpi%n_size     
    ELSE
      n_members = MIN(nkpt,fmpi%isize)
      IF (judft_was_argument("-min_pe_per_kpt")) THEN
         txt=judft_string_for_argument("-min_pe_per_kpt")
         READ(txt,*) n_size_min
         WRITE(*,*) "Trying to use ",n_size_min," PE per kpt"
         n_members = MIN(n_members , CEILING(REAL(fmpi%isize)/n_size_min) )
      ENDIF
      DO
         IF ((MOD(fmpi%isize,n_members) == 0).AND.(MOD(nkpt,n_members) == 0) ) EXIT
         n_members = n_members - 1
      ENDDO
    endif  

    !fmpi%n_groups = nkpt/n_members
    fmpi%n_size   = fmpi%isize/n_members
    !fmpi%n_stride = n_members
    IF (fmpi%irank == 0) THEN
       WRITE(*,*) 'k-points in parallel: ',n_members
       WRITE(*,*) "pe's per k-point:     ",fmpi%n_size
       WRITE(*,*) '# of k-point loops:   ',nkpt/n_members
       if (mod(nkpt,n_members).ne.0) then
         Write(*,*) 'your k-point parallelism is not fully load-balanced'
       endif

       IF((REAL(nbasfcn) / REAL(fmpi%n_size)).LE.20) THEN
          WRITE(*,*) ''
          WRITE(*,*) '========================================'
          WRITE(*,*) 'WARNING:'
          WRITE(*,*) 'The chosen parallelization scheme implies very few LAPW basis functions'
          WRITE(*,*) 'per MPI process. This may lead to poor performance and other problems.'
          IF((nkpt/n_members).GE.4*fmpi%isize) THEN
             WRITE(*,*) ''
             WRITE(*,*) 'You may want to adjust the number of MPI processes such that the k-point'
             WRITE(*,*) 'parallelization is increased.'
             WRITE(*,*) ''
          END IF
          WRITE(*,*) 'Fleur will proceed with the calculation.'
          WRITE(*,*) '========================================'
          WRITE(*,*) ''
       END IF

    ENDIF
  END SUBROUTINE priv_distribute_k

  SUBROUTINE priv_create_comm(nkpt,neigd,fmpi)
    use m_types
    implicit none
    INTEGER,INTENT(in)      :: nkpt,neigd
    TYPE(t_mpi),INTENT(inout)    :: fmpi
#ifdef CPP_MPI
    INTEGER :: n_members,n,i,ierr,sub_group,world_group,n_start
    INTEGER :: i_mygroup(fmpi%n_size)
    LOGICAL :: compact ! Deside how to distribute k-points

    compact = .true.
    n_members = fmpi%isize/fmpi%n_size

    ! now, we make the groups


    IF (compact) THEN

        ! This will distribute sub ranks in a compact manner.
        ! For example, if nkpt = 8 and fmpi%isize = 6:

        !  -----------------------------------
        ! |  0  |  1  |  2  |  3  |  4  |  5  |    fmpi%irank
        !  -----------------------------------
        ! |  0  |  1  |  3  |  0  |  1  |  2  |    fmpi%n_rank
        !  -----------------------------------
        ! |        1        |        2        |    k - points
        ! |        3        |        4        |
        ! |        5        |        6        |
        ! |        7        |        8        |
        !  -----------------------------------

        n_start = INT(fmpi%irank/fmpi%n_size) + 1
        i_mygroup(1) = (n_start-1) * fmpi%n_size
        do i = 2, fmpi%n_size
           i_mygroup(i) = i_mygroup(i-1) + 1
        enddo

    ELSE

        ! This will distribute sub ranks in a spread manner.
        ! For example, if nkpt = 8 and fmpi%isize = 6:

        !  -----------------------------------
        ! |  0  |  1  |  2  |  3  |  4  |  5  |    fmpi%irank
        !  -----------------------------------
        ! |  0  |  1  |  3  |  0  |  1  |  2  |    fmpi%n_rank
        !  -----------------------------------
        ! |  1  |  2  |  1  |  2  |  1  |  2  |    k - points
        ! |  3  |  4  |  3  |  4  |  3  |  4  |
        ! |  5  |  6  |  5  |  6  |  5  |  6  |
        ! |  7  |  8  |  7  |  8  |  7  |  8  |
        !  -----------------------------------

        n_start = MOD(fmpi%irank,n_members) + 1
        !!      n_start = INT(irank/n_size) * n_size
        n = 0
        DO i = n_start,fmpi%isize,n_members
        !!      DO i = n_start+1,n_start+n_size
           n = n+1
           i_mygroup(n) = i-1
        ENDDO

    ENDIF ! compact

    CALL MPI_COMM_GROUP (fmpi%MPI_COMM,WORLD_GROUP,ierr)
    CALL MPI_GROUP_INCL (WORLD_GROUP,fmpi%n_size,i_mygroup,SUB_GROUP,ierr)
    CALL MPI_COMM_CREATE (fmpi%MPI_COMM,SUB_GROUP,fmpi%SUB_COMM,ierr)
    !write (*,"(a,i0,100i4)") "MPI:",fmpi%sub_comm,fmpi%irank,fmpi%n_groups,fmpi%n_size,n,i_mygroup

    CALL MPI_COMM_RANK (fmpi%SUB_COMM,fmpi%n_rank,ierr)
    ALLOCATE(fmpi%ev_list(neigd/fmpi%n_size+1))
    fmpi%ev_list=[(i,i=fmpi%n_rank+1,neigd,fmpi%n_size)]
    fmpi%diag_sub_comm=fmpi%sub_comm !default both are the same...
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

  subroutine priv_redist_for_diag(fmpi)
    use m_types_mpi
    type(t_mpi),intent(inout):: fmpi
#ifdef CPP_MPI
    IF (fmpi%n_rank==0) THEN
      fmpi%diag_sub_comm=MPI_COMM_SELF
      fmpi%pe_diag=.true. !actually the default
    ELSE
      fmpi%diag_sub_comm=MPI_COMM_NULL
      fmpi%pe_diag=.false.
    ENDIF
#endif
    end

   subroutine priv_distribute_gpu(fmpi)
#ifdef _OPENACC
   use openacc
#endif    
    use m_types_mpi
    type(t_mpi),intent(in):: fmpi
    INTEGER :: i, isize, gpus,localrank

#ifdef _OPENACC
call timestart("Distribute GPUs")
    gpus=acc_get_num_devices(acc_device_nvidia)
#ifdef CPP_MPI
    if (fmpi%irank==0) write(*,*) "Number of GPU per node/MPI:",gpus
    CALL MPI_COMM_SIZE(fmpi%mpi_comm_same_node,isize,i)
    if (isize>1) THEN
      if (fmpi%irank==0) write(*,*) "Number of MPI/PE per node:",isize
      if (gpus>isize) call judft_warn("You use more GPU/node as MPI-PEs/node running. This will underutilize the GPUs")
      CALL MPI_COMM_RANK(fmpi%mpi_comm_same_node,localrank,i)
      if (gpus>1) THEN 
         call acc_set_device_num(mod(localrank,gpus),acc_device_nvidia)
         write(*,*) "Assigning PE:",fmpi%irank," to local GPU:",mod(localrank,gpus)
      endif   
    ENDIF
#else
    write(*,*) "Number of GPU    :",gpus
#endif
   call timestop("Distribute GPUs")
#endif
    end subroutine

END MODULE m_setupMPI
