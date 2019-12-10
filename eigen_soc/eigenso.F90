MODULE m_eigenso
  !
  !*********************************************************************
  !     sets ur and solves the spin-orbit eigenvalue problem in the
  !     second variation procedure.
  !
  !     way: takes e.v. and e.f. from previous scalar-rel. calc.
  !     makes spin-orbit matrix elements solves e.v. and put it on 'eig'
  !
  !     Tree:  eigenso-|- readPotential
  !                    |- spnorb  : sets up s-o parameters 
  !                    |    |- soinit - sorad  : radial part
  !                    |    |- sgml            : diagonal angular parts
  !                    |    |- anglso          : non-diagonal -"-
  !                    |
  !                    |- alineso : sets up and solves e.v. problem
  !                         |- hsohelp
  !                         |- hsoham
  !
  !**********************************************************************
  !
CONTAINS
  SUBROUTINE eigenso(eig_id,mpi,DIMENSION,stars,vacuum,atoms,sphhar,&
                     obsolete,sym,cell,noco,input,kpts,oneD,vTot,enpara,results,hub1)

    USE m_types
    USE m_eig66_io, ONLY : read_eig,write_eig
    USE m_spnorb 
    USE m_alineso
    USE m_judft
#ifdef CPP_MPI
    USE m_mpi_bc_pot
#endif
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_dimension),INTENT(IN)  :: DIMENSION
    TYPE(t_oneD),INTENT(IN)       :: oneD
    TYPE(t_obsolete),INTENT(IN)   :: obsolete
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_vacuum),INTENT(IN)     :: vacuum
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_kpts),INTENT(IN)       :: kpts
    TYPE(t_sphhar),INTENT(IN)     :: sphhar
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_potden),INTENT(IN)     :: vTot
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_results),INTENT(INOUT) :: results
    TYPE(t_hub1ham),OPTIONAL,INTENT(INOUT) :: hub1

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif

    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id       
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER i,j,nk,nk_i,jspin,n ,l
    ! INTEGER n_loc,n_plus,i_plus,
    INTEGER nsz,nmat,n_stride
    LOGICAL l_socvec   !,l_all
    INTEGER wannierspin
    TYPE(t_usdus)        :: usdus
    !     ..
    !     .. Local Arrays..
    CHARACTER*3 chntype

    TYPE(t_rsoc) :: rsoc
    INTEGER, ALLOCATABLE :: neigBuffer(:,:)
    REAL,    ALLOCATABLE :: eig_so(:), eigBuffer(:,:,:)
    COMPLEX, ALLOCATABLE :: zso(:,:,:)

    TYPE(t_mat)::zmat
    TYPE(t_lapw)::lapw

    INTEGER :: ierr
    
    !  ..

    INQUIRE (4649,opened=l_socvec)

    ! To be consistent with angles should be redefined here!
    !noco%theta= -noco%theta
    !noco%phi=   noco%phi+pi_const
    ! now the definition of rotation matrices
    ! is equivalent to the def in the noco-routines

    ALLOCATE(  usdus%us(0:atoms%lmaxd,atoms%ntype,input%jspins), usdus%dus(0:atoms%lmaxd,atoms%ntype,input%jspins),&
         usdus%uds(0:atoms%lmaxd,atoms%ntype,input%jspins),usdus%duds(0:atoms%lmaxd,atoms%ntype,input%jspins),&
         usdus%ddn(0:atoms%lmaxd,atoms%ntype,input%jspins),&
         usdus%ulos(atoms%nlod,atoms%ntype,input%jspins),usdus%dulos(atoms%nlod,atoms%ntype,input%jspins),&
         usdus%uulon(atoms%nlod,atoms%ntype,input%jspins),usdus%dulon(atoms%nlod,atoms%ntype,input%jspins))
   
    IF (input%l_wann.OR.l_socvec) THEN
       wannierspin = 2
    ELSE
       wannierspin = input%jspins
    ENDIF

    !
    !---> set up and solve the eigenvalue problem
    ! --->    radial k-idp s-o matrix elements calc. and storage
    !
#if defined(CPP_MPI)
    !RMA synchronization
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
    CALL timestart("eigenso: spnorb")
    !  ..

    !Get spin-orbit coupling matrix elements
    CALL spnorb( atoms,noco,input,mpi, enpara,vTot%mt,usdus,rsoc,.TRUE.,hub1)
    !


    ALLOCATE (eig_so(2*DIMENSION%neigd))
    ALLOCATE (eigBuffer(2*DIMENSION%neigd,kpts%nkpt,wannierspin))
    ALLOCATE (neigBuffer(kpts%nkpt,wannierspin))
    results%eig = 1.0e300
    eigBuffer = 1.0e300
    results%neig = 0
    neigBuffer = 0
    rsoc%soangl(:,:,:,:,:,:) = CONJG(rsoc%soangl(:,:,:,:,:,:))
    CALL timestop("eigenso: spnorb")
    !
    !--->    loop over k-points: each can be a separate task
    DO nk_i=1,SIZE(mpi%k_list)
        nk=mpi%k_list(nk_i)
     !DO nk = mpi%n_start,n_end,n_stride
       CALL lapw%init(input,noco, kpts,atoms,sym,nk,cell,.FALSE., mpi)
       ALLOCATE( zso(lapw%nv(1)+atoms%nlotot,2*DIMENSION%neigd,wannierspin))
       zso(:,:,:) = CMPLX(0.0,0.0)
       CALL timestart("eigenso: alineso")
       CALL alineso(eig_id,lapw, mpi,DIMENSION,atoms,sym,kpts,&
            input,noco,cell,oneD,nk,usdus,rsoc,nsz,nmat, eig_so,zso)
       CALL timestop("eigenso: alineso")
       IF (mpi%irank.EQ.0) THEN
          WRITE (6,FMT=8010) nk,nsz
          WRITE (6,FMT=8020) (eig_so(i),i=1,nsz)
       ENDIF
8010   FORMAT (1x,/,/,' #k=',i6,':',/,&
            ' the',i4,' SOC eigenvalues are:')
8020   FORMAT (5x,5f12.6)

       IF (mpi%n_rank==0) THEN
          IF (input%eonly) THEN
             CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz))
             STOP 'jspin is undefined here (eigenso - eonly branch)'
             eigBuffer(:nsz,nk,jspin) = eig_so(:nsz)
             neigBuffer(nk,jspin) = nsz
          ELSE
             CALL zmat%alloc(.FALSE.,SIZE(zso,1),nsz)
             DO jspin = 1,wannierspin
                CALL timestart("eigenso: write_eig")  
                zmat%data_c=zso(:,:nsz,jspin)
                CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz),zmat=zmat)
                eigBuffer(:nsz,nk,jspin) = eig_so(:nsz)
                neigBuffer(nk,jspin) = nsz
                CALL timestop("eigenso: write_eig")
             ENDDO
          ENDIF ! (input%eonly) ELSE
       ENDIF ! n_rank == 0
       DEALLOCATE (zso)
    ENDDO ! DO nk 

#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(neigBuffer,results%neig,kpts%nkpt*input%jspins,MPI_INTEGER,MPI_SUM,mpi%mpi_comm,ierr)
    CALL MPI_ALLREDUCE(eigBuffer(:2*dimension%neigd,:,1:input%jspins),&
                     results%eig(:2*dimension%neigd,:,1:input%jspins),&
        2*dimension%neigd*kpts%nkpt*input%jspins,MPI_DOUBLE_PRECISION,MPI_MIN,mpi%mpi_comm,ierr)
                       
                       
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#else
    results%neig(:,1:input%jspins) = neigBuffer(:,1:input%jspins)
    results%eig(:2*dimension%neigd,:,1:input%jspins) = eigBuffer(:2*dimension%neigd,:,1:input%jspins)
!    results%neig(:,:) = neigBuffer(:,:)
!    results%eig(:2*dimension%neigd,:,:) = eigBuffer(:2*dimension%neigd,:,:)
#endif

    RETURN
  END SUBROUTINE eigenso
END MODULE m_eigenso
