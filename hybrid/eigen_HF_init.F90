
module m_eigen_hf_init
!
!     preparations for HF and hybrid functional calculation
!
contains
  subroutine eigen_hf_init(hybrid,kpts,atoms,input,dimension,hybdat,irank2,isize2,l_real)
    USE m_types
    USE m_read_core
    USE m_util
    USE m_io_hybrid
  implicit none
  TYPE(t_hybrid),INTENT(INOUT)     :: hybrid
  TYPE(t_kpts),INTENT(IN)          :: kpts
  TYPE(t_atoms),INTENT(IN)         :: atoms
  TYPE(t_input),INTENT(IN)         :: input
  TYPE(t_dimension),INTENT(IN)     :: dimension
  INTEGER,INTENT(OUT)              :: irank2(:),isize2(:)
  TYPE(t_hybdat),INTENT(OUT)       :: hybdat
  LOGICAL,INTENT(IN)               :: l_real
  
 

  
  LOGICAL,SAVE :: init_vex=.true.
  INTEGER,SAVE :: nohf_it
  integer:: itype,ieq,l,m,i,nk,l1,l2,m1,m2,ok

  
   IF( .NOT. hybrid%l_hybrid ) THEN
        hybrid%l_calhf  = .false.
      ELSE IF( (all(hybrid%ddist .lt. 1E-5) .or. init_vex .or. nohf_it >= 50 ) .and. input%imix .gt. 10) THEN
        hybrid%l_calhf  = .true.
        init_vex = .false.
        nohf_it  = 0
      ELSE IF( input%imix .lt. 10 ) THEN
        hybrid%l_calhf  = .true.
        init_vex = .true.
        nohf_it  = 0
      ELSE
        hybrid%l_calhf = .false.
        nohf_it = nohf_it + 1
      END IF

      IF( hybrid%l_calhf ) THEN
         
        !initialize hybdat%gridf for radial integration
        CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

        !Alloc variables
        ALLOCATE(hybdat%lmaxc(atoms%ntype))
        ALLOCATE(hybdat%ne_eig(kpts%nkpt),hybdat%nbands(kpts%nkpt))
        ALLOCATE(hybdat%nobd(kpts%nkptf))
        ALLOCATE(hybdat%bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype))
        ALLOCATE(hybdat%bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype))
        ALLOCATE(hybdat%bas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype))
        ALLOCATE(hybdat%drbas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype))
        
        !sym%tau = oneD%ods%tau

        ! preparations for core states
        CALL core_init(dimension,input,atoms, hybdat%lmaxcd,hybdat%maxindxc)
        ALLOCATE( hybdat%nindxc(0:hybdat%lmaxcd,atoms%ntype), stat = ok )
        IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation hybdat%nindxc'
        ALLOCATE( hybdat%core1(atoms%jmtd,hybdat%maxindxc,0:hybdat%lmaxcd,atoms%ntype), stat=ok )
        IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation core1'
        ALLOCATE( hybdat%core2(atoms%jmtd,hybdat%maxindxc,0:hybdat%lmaxcd,atoms%ntype), stat=ok )
        IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation core2'
        ALLOCATE( hybdat%eig_c(hybdat%maxindxc,0:hybdat%lmaxcd,atoms%ntype), stat=ok      )
        IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation hybdat%eig_c'
        hybdat%nindxc = 0 ; hybdat%core1 = 0 ; hybdat%core2 = 0 ; hybdat%eig_c = 0

      
        ALLOCATE( hybdat%nbasm(kpts%nkptf) ,stat=ok)
        IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation hybdat%nbasm'
        DO nk = 1,kpts%nkptf
          hybdat%nbasm(nk) = hybrid%nbasp + hybrid%ngptm(nk)
        END DO

        ! pre-calculate gaunt coefficients

        hybdat%maxfac   = max(2*atoms%lmaxd+hybrid%maxlcutm1+1,2*hybdat%lmaxcd+2*atoms%lmaxd+1)
        ALLOCATE ( hybdat%fac( 0:hybdat%maxfac),hybdat%sfac( 0:hybdat%maxfac),stat=ok)
        IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation fac,hybdat%sfac'
        hybdat%fac(0)   = 1
        hybdat%sfac(0)  = 1
        DO i=1,hybdat%maxfac
          hybdat%fac(i)    = hybdat%fac(i-1)*i            ! hybdat%fac(i)    = i!
          hybdat%sfac(i)   = hybdat%sfac(i-1)*sqrt(i*1d0) ! hybdat%sfac(i)   = sqrt(i!)
        END DO


        ALLOCATE ( hybdat%gauntarr( 2, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:hybrid%maxlcutm1, -atoms%lmaxd:atoms%lmaxd ,-hybrid%maxlcutm1:hybrid%maxlcutm1 ),stat=ok)
        IF( ok .ne. 0 ) STOP 'eigen: failure allocation hybdat%gauntarr'
        hybdat%gauntarr = 0
        DO l2 = 0,atoms%lmaxd
          DO l1 = 0,atoms%lmaxd
            DO l = abs(l1-l2),min(l1+l2,hybrid%maxlcutm1)
              DO m = -l,l
                DO m1 = -l1,l1
                  m2 = m1 + m ! Gaunt condition -m1+m2-m = 0
                  IF(abs(m2).le.l2) hybdat%gauntarr(1,l1,l2,l,m1,m) = gaunt(l1,l2,l,m1,m2,m,hybdat%maxfac,hybdat%fac,hybdat%sfac)
                  m2 = m1 - m ! switch role of l2-index
                  IF(abs(m2).le.l2) hybdat%gauntarr(2,l1,l2,l,m1,m) = gaunt(l2,l1,l,m2,m1,m,hybdat%maxfac,hybdat%fac,hybdat%sfac)
                END DO
              END DO
            END DO
          END DO
        END DO

      END IF ! hybrid%l_calhf

#ifdef CPP_NEVER   
!#       ifdef CPP_MPI
          IF ( hybrid%l_calhf ) THEN
            ALLOCATE( nkpt_EIBZ(kpts%nkpt) )
            DO nk = 1,kpts%nkpt
              nkpt_EIBZ(nk) = symm_hf_nkpt_EIBZ(kpts%nkptf,kpts%nkpt3,nk,kpts%bkf,sym%nsym, sym%nop,sym%mrot,sym%invtab)
            END DO
            comm=init_work_dist(mpi%mpi_comm,mpi%irank,mpi%isize,kpts%nkpt,nkpt_EIBZ)
            skip_kpt = ( comm == MPI_COMM_NULL )
            DO nk = 1,kpts%nkpt
              ! determine rank in new communicator and size of the communicator
              IF ( skip_kpt(nk) ) THEN
                irank2(nk) = MPI_PROC_NULL
              ELSE
                CALL MPI_COMM_RANK(comm(nk),irank2(nk),ierr(1))
                CALL MPI_COMM_SIZE(comm(nk),isize2(nk),ierr(1))
              END IF
            END DO
          ELSE
            irank2   = 0
            isize2   = 1
            skip_kpt = .false.
          END IF
#       else
          irank2   = 0
          isize2   = 1
          !skip_kpt = .false.
#       endif
        end subroutine eigen_hf_init
      end module m_eigen_hf_init
