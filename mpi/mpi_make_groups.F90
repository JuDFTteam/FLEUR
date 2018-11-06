!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpimakegroups
  use m_juDFT
CONTAINS
  SUBROUTINE mpi_make_groups(&
       mpi,dimension,kpts, input,atoms,noco,&
       mlotot,mlolotot,&
       n_start,n_groups,n,matsize,ne, n_rank,n_size,SUB_COMM)
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
! The results (containing a subset of ew & ev's) are written as separate
! records on the file 'eig' with the vacuum energy parameter of the
! marked (*) records set to 999.9 to indicate that this is only one part
! of a k-point's results. This marker is needed in routines fermie and 
! cdnval.
!          G.B. `99
!
!------------------------------------------------------------------------
    USE m_types
    IMPLICIT NONE
    
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: dimension
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER, INTENT (IN)  :: mlotot
    INTEGER, INTENT (IN)  :: mlolotot    
    INTEGER, INTENT (OUT) :: n_start,n_groups,n,SUB_COMM
    INTEGER, INTENT (OUT) :: matsize,n_rank,n_size,ne

      INTEGER i,n_members
      INTEGER, ALLOCATABLE :: i_mygroup(:)

      INCLUDE 'mpif.h'
      INTEGER WORLD_GROUP,SUB_GROUP
      INTEGER ierr(3)
      LOGICAL l_cm
!
! first determine the number of groups of k-points to process
!
      n_groups = 0
      IF (kpts%nkpt.GT.mpi%isize) THEN           ! if there are more k-points than PEs

        IF (mod(kpts%nkpt,mpi%isize).EQ.0) THEN  ! maybe kpts%nkpt is a multiple of mpi%isize
          n_groups = kpts%nkpt/mpi%isize
          n_size = 1
        ELSE                            ! or an integer fraction of mpi%isize fits
          DO i=2,mpi%isize
            IF (mod(mpi%isize,i).EQ.0) THEN
              n_size = i
              n_groups = kpts%nkpt * i/mpi%isize
              IF (mod(kpts%nkpt,mpi%isize/i).EQ.0) GOTO 990
            ENDIF
          ENDDO
          n_groups = kpts%nkpt               ! or use all PE's per k-point
          n_size = mpi%isize
        ENDIF

      ELSEIF (kpts%nkpt.LT.mpi%isize) THEN       ! if there are more PEs than k-points

        IF (mod(mpi%isize,kpts%nkpt).EQ.0) THEN  ! maybe mpi%isize is a multiple of kpts%nkpt 
           n_groups = 1
           n_size = mpi%isize/kpts%nkpt
        ELSE                            ! or an integer fraction of kpts%nkpt fits
          DO i=kpts%nkpt-1,2,-1
            IF (mod(kpts%nkpt,i).EQ.0) THEN
               n_groups = kpts%nkpt/i
               n_size = mpi%isize/i
               IF (mod(mpi%isize,i).EQ.0) GOTO 990
            ENDIF
          ENDDO
          n_groups = kpts%nkpt               ! or use all PE's per k-point
          n_size = mpi%isize
        ENDIF

      ELSE
!
! if there are as many pe's as k-points (isize = nkpt), use one PE per
! kpoint (n_size = 1) and finish in a single k-loop (n_groups = 1)
!
        n_groups = 1
        n_size = 1
      ENDIF



 990  IF (n_groups==0)  CALL juDFT_error("mpi_make_groups:1",calledby ="mpi_make_groups")
      n_members = kpts%nkpt/n_groups 
!
! check different algorithm
!
      CALL check_memory(DIMENSION,input,atoms, mlotot,mlolotot,noco,kpts,mpi, n_size)

      write(*,*) n_size
      n_members = MIN(kpts%nkpt,mpi%isize)
      n_members = MIN(n_members , CEILING(REAL(mpi%isize)/n_size) ) + 1
      
      l_cm = .false.
      DO WHILE (.not.l_cm) 
        n_members = n_members - 1
        IF ((mod(mpi%isize,n_members) == 0).AND.&
     &      (mod(kpts%nkpt,n_members) == 0) ) THEN
           l_cm = .true.
        ENDIF
      ENDDO
      n_groups = kpts%nkpt/n_members
      n_size = mpi%isize/n_members
      IF (mpi%irank == 0) THEN
        write(*,*) 'k-points in parallel: ',n_members
        write(*,*) "pe's per k-point:     ",n_size
        write(*,*) '# of k-point loops:   ',n_groups
      ENDIF
!
! now, we make the groups
!
      n_start = mod(mpi%irank,n_members) + 1
!!      n_start = INT(irank/n_size) * n_size
      ALLOCATE ( i_mygroup(n_size) )
      n = 0
      DO i = n_start,mpi%isize,n_members
!!      DO i = n_start+1,n_start+n_size
        n = n+1
        i_mygroup(n) = i-1
      ENDDO

!      write (*,*) irank,n_groups,n_start,i_mygroup

      CALL MPI_COMM_GROUP (mpi%MPI_COMM,WORLD_GROUP,ierr)
      CALL MPI_GROUP_INCL (WORLD_GROUP,n_size,i_mygroup, SUB_GROUP,ierr)
      CALL MPI_COMM_CREATE (mpi%MPI_COMM,SUB_GROUP,SUB_COMM,ierr)

      CALL MPI_COMM_RANK (SUB_COMM,n_rank,ierr)
!
! determine number of columns per group
!
      n=0
      DO  i=1+n_rank, dimension%nbasfcn, n_size
        n=n+1
      ENDDO
      IF (n_size.GT.1) THEN
        matsize = dimension%nbasfcn * n
      ELSE
        matsize = (dimension%nbasfcn+1)*dimension%nbasfcn/2
      ENDIF
!
        ne = dimension%neigd
      ne = max(ne,5)

      DEALLOCATE (i_mygroup)

      END SUBROUTINE mpi_make_groups

!----------------------------------------------------------------------
      SUBROUTINE check_memory(DIMENSION,input,atoms, mlotot,mlolotot, noco,kpts,mpi, n_size)

!
! check the free and the (approximate) required memory ;
! determine minimal n_size to fit into the memory (hopefully).
!
        USE m_types
      IMPLICIT NONE
      type(t_mpi),INTENT(IN)         :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_kpts),INTENT(IN)        :: kpts
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_noco),INTENT(IN)        :: noco
      TYPE(t_input),INTENT(IN)       :: input

      INTEGER, INTENT (IN)  :: mlotot,mlolotot    
      INTEGER, INTENT (OUT) :: n_size

      INTEGER  err, mb
      INTEGER*8 mem, matsz, m_h
      REAL, ALLOCATABLE :: test(:)

      n_size = CEILING( real(mpi%isize)/min(kpts%nkpt,mpi%isize) )
      
 10   CONTINUE
!
! some basic arrays allocated in eigen()
!
      mem = (dimension%lmplmd+1)*atoms%ntype*4                       ! tlmplm%tuu,tlmplm%tdd etc.
      mem = mem + (dimension%lmd+1)*(2*atoms%llod+1)*max(mlotot,1)*2 ! tlmplm%tuulo ...
      mem = mem + (2*atoms%llod+1)**2 * max(mlolotot,1)    ! tlmplm%tuloulo
      IF (noco%l_noco) mem = mem * 2                      ! both spins
      mem = mem + 49*atoms%n_u*input%jspins*2                      ! lda+U, *2 for complex
      mem = mem+INT((dimension%nbasfcn*2+(dimension%lmd+1)*atoms%ntype)*0.5)+1 ! tlmplm%ind, *0.5 for integer

      matsz = dimension%nbasfcn * CEILING(REAL(dimension%nbasfcn)/n_size) ! size of a, b
#ifdef CPP_INVERSION
      mem = mem + 2 * matsz                          ! real case
#else
      mem = mem + 2 * matsz * 2                      ! complec case
#endif
!
! now the arrays in hssphn()
!
      m_h = dimension%nvd*(dimension%lmd+1)*4  + dimension%nvd*8 + atoms%nlod            ! ar, ai ..., cph, rph, vk, gk
      m_h = m_h + 2 * (2*atoms%llod+1)**2 * atoms%nlod * 3 * 2   ! alo,blo,clo
      IF (noco%l_ss) m_h = m_h * 2
      m_h = m_h + dimension%nvd*(5+atoms%lmaxd)                      ! axr, ... plegend
      IF (noco%l_ss.OR.noco%l_constr.OR.(noco%l_noco.AND.noco%l_soc)) THEN
        m_h = m_h + dimension%nvd*(atoms%lmaxd+1)*atoms%ntype*2*2          ! fj,gj
      ELSE 
        m_h = m_h + dimension%nvd*(atoms%lmaxd+1)*atoms%ntype*2            ! fj,gj
      ENDIF
      IF (noco%l_noco.AND.noco%l_soc) THEN
        m_h = m_h + dimension%nvd*(atoms%lmaxd+4)
      ENDIF
      IF (noco%l_constr) THEN
        m_h = m_h + (atoms%lmaxd+1)*atoms%ntype
      ENDIF
      IF (noco%l_noco.AND.(.NOT.noco%l_ss)) THEN
        matsz = (dimension%nvd+mlotot) * CEILING(REAL(dimension%nvd+mlotot)/n_size)
        m_h = m_h + matsz * 2 * 2                    ! aahlp,bbhlp
      ENDIF
!
! see, whether it fits
!
      mb = (mem+m_h)*8/(1024)**2
      ALLOCATE ( test(mem+m_h) , stat = err)
      WRITE(*,*) mb,'Mbytes needed  in hssphn!',err,mem
      IF ( err /= 0 ) THEN
        n_size = n_size * 2
        IF (n_size > mpi%isize) THEN
          mb = (mem+m_h)*8/(1024)**2
          WRITE(*,*) mb,'Mbytes needed  in hssphn!'
          CALL juDFT_error("mpi_make_groups: memory too small!",calledby ="mpi_make_groups")
        ENDIF 
        GOTO 10
      ENDIF
      DEALLOCATE (test)
!
! now, allocate z and jump into chani
!      
      matsz = dimension%nbasfcn * CEILING(REAL(dimension%nbasfcn)/n_size)   ! size of z
#ifdef CPP_INVERSION
      mem = mem + matsz                             ! real case
#else
      mem = mem + matsz * 2                         ! complex case
#endif
      mem = mem + matsz * 2 * 3                     ! asca,bsca,eigvec
      mem = mem + dimension%nbasfcn 
#ifdef CPP_INVERSION
      mem = mem + matsz                             ! work
#else
      mem = mem + matsz * 2                         ! work, rwork neglected
#endif
!
! see, whether it fits
!
      mb = (mem)*8/(1024)**2
      ALLOCATE ( test(mem) , stat = err)
      WRITE(*,*) mb,'Mbytes needed  in chani !',err,mem
      IF ( err /= 0 ) THEN
        n_size = n_size * 2
        IF (n_size > mpi%isize) THEN
          mb = (mem)*8/(1024)**2
          WRITE(*,*) mb,'Mbytes needed  in chani !'
          CALL juDFT_error("mpi_make_groups: memory too small!",calledby&
     &         ="mpi_make_groups")
        ENDIF
        GOTO 10
      ENDIF
      DEALLOCATE (test)
      
      END SUBROUTINE check_memory
!----------------------------------------------------------------------

      END MODULE m_mpimakegroups
