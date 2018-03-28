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
                     obsolete,sym,cell,noco,input,kpts,oneD,vTot,enpara)

    USE m_eig66_io, ONLY : read_eig,write_eig
    USE m_spnorb 
    USE m_alineso
    USE m_types
    USE m_judft
#ifdef CPP_MPI
    USE m_mpi_bc_pot
#endif
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_obsolete),INTENT(IN)  :: obsolete
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_potden),INTENT(IN)    :: vTot
    TYPE(t_enpara),INTENT(IN)    :: enpara
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id       
    !     ..
    !     ..
    !     .. Local Scalars ..
    INTEGER i,j,nk,jspin,n ,l
    INTEGER n_loc,n_plus,i_plus,n_end,nsz,nmat
    LOGICAL l_socvec   !,l_all
    INTEGER wannierspin
    TYPE(t_usdus):: usdus
    !     ..
    !     .. Local Arrays..
    CHARACTER*3 chntype

    TYPE(t_rsoc) :: rsoc
    REAL,    ALLOCATABLE :: eig_so(:) 
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

    ALLOCATE(  usdus%us(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd), usdus%dus(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),&
         usdus%uds(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),usdus%duds(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),&
         usdus%ddn(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd),&
         usdus%ulos(atoms%nlod,atoms%ntype,DIMENSION%jspd),usdus%dulos(atoms%nlod,atoms%ntype,DIMENSION%jspd),&
         usdus%uulon(atoms%nlod,atoms%ntype,DIMENSION%jspd),usdus%dulon(atoms%nlod,atoms%ntype,DIMENSION%jspd))
   
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

    IF (SIZE(noco%theta)>1) CALL judft_warn("only first SOC-angle used in second variation")
    !Get spin-orbit coupling matrix elements
    CALL spnorb( atoms,noco,input,mpi, enpara,vTot%mt,usdus,rsoc,.TRUE.)
    !


    ALLOCATE( eig_so(2*DIMENSION%neigd) )
    rsoc%soangl(:,:,:,:,:,:,1) = CONJG(rsoc%soangl(:,:,:,:,:,:,1))
    CALL timestop("eigenso: spnorb")
    !
    !--->    loop over k-points: each can be a separate task
    !
    n_loc = INT(kpts%nkpt/mpi%isize)
    n_plus = kpts%nkpt - mpi%isize*n_loc
    i_plus = -1
    IF (mpi%irank.LT.n_plus) i_plus = 0
    n_end = (mpi%irank+1)+(n_loc+i_plus)*mpi%isize
    !
    !--->  start loop k-pts
    !
    DO  nk = mpi%irank+1,n_end,mpi%isize
       CALL lapw%init(input,noco, kpts,atoms,sym,nk,cell,.FALSE., mpi)
       ALLOCATE( zso(lapw%nv(1)+atoms%nlotot,2*DIMENSION%neigd,wannierspin))
       zso(:,:,:) = CMPLX(0.0,0.0)
       CALL timestart("eigenso: alineso")
       CALL alineso(eig_id,lapw, mpi,DIMENSION,atoms,sym,kpts,&
            input,noco,cell,oneD,nk,usdus,rsoc,nsz,nmat, eig_so,zso)
       CALL timestop("eigenso: alineso")
       IF (mpi%irank.EQ.0) THEN
          WRITE (16,FMT=8010) nk,nsz
          WRITE (16,FMT=8020) (eig_so(i),i=1,nsz)
       ENDIF
8010   FORMAT (1x,/,/,' #k=',i6,':',/,&
            ' the',i4,' SOC eigenvalues are:')
8020   FORMAT (5x,5f12.6)

       IF (input%eonly) THEN
          CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz))

       ELSE          
          CALL zmat%alloc(.FALSE.,SIZE(zso,1),nsz)
          DO jspin = 1,wannierspin
             CALL timestart("eigenso: write_eig")  
             zmat%data_c=zso(:,:nsz,jspin)
             CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz),zmat=zmat)

             CALL timestop("eigenso: write_eig")  
          ENDDO
       ENDIF ! (input%eonly) ELSE
       deallocate(zso)
    ENDDO ! DO nk 
    RETURN
  END SUBROUTINE eigenso
END MODULE m_eigenso
