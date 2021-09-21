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

#ifdef CPP_MPI
    use mpi 
#endif
CONTAINS
  SUBROUTINE eigenso(eig_id,fmpi,stars,sphhar,nococonv,vTot,enpara,results,hub1inp,hub1data,fi)

    USE m_types
    USE m_constants
    USE m_eig66_io, ONLY : read_eig,write_eig
    USE m_spnorb
    USE m_alineso
    USE m_judft
    USE m_unfold_band_kpts
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)        :: fmpi
    type(t_fleurinput), intent(in) :: fi
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_sphhar),INTENT(IN)     :: sphhar
    TYPE(t_potden),INTENT(IN)     :: vTot
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_results),INTENT(INOUT) :: results
    TYPE(t_hub1inp),OPTIONAL,INTENT(IN) :: hub1inp
    TYPE(t_hub1data),OPTIONAL,INTENT(INOUT) :: hub1data

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

    COMPLEX              :: unfoldingBuffer(SIZE(results%unfolding_weights,1),fi%kpts%nkpt,fi%input%jspins) ! needed for unfolding bandstructure fmpi case

    REAL,    ALLOCATABLE :: eig_so(:), eigBuffer(:,:,:)
    COMPLEX, ALLOCATABLE :: zso(:,:,:)

    TYPE(t_mat)::zmat
    TYPE(t_lapw)::lapw

    INTEGER :: ierr, jsp

    !  ..

    INQUIRE (4649,opened=l_socvec)

    ! To be consistent with angles should be redefined here!
    !noco%theta= -noco%theta
    !noco%phi=   noco%phi+pi_const
    ! now the definition of rotation matrices
    ! is equivalent to the def in the noco-routines

    ALLOCATE(  usdus%us(0:fi%atoms%lmaxd,fi%atoms%ntype,fi%input%jspins), usdus%dus(0:fi%atoms%lmaxd,fi%atoms%ntype,fi%input%jspins),&
         usdus%uds(0:fi%atoms%lmaxd,fi%atoms%ntype,fi%input%jspins),usdus%duds(0:fi%atoms%lmaxd,fi%atoms%ntype,fi%input%jspins),&
         usdus%ddn(0:fi%atoms%lmaxd,fi%atoms%ntype,fi%input%jspins),&
         usdus%ulos(fi%atoms%nlod,fi%atoms%ntype,fi%input%jspins),usdus%dulos(fi%atoms%nlod,fi%atoms%ntype,fi%input%jspins),&
         usdus%uulon(fi%atoms%nlod,fi%atoms%ntype,fi%input%jspins),usdus%dulon(fi%atoms%nlod,fi%atoms%ntype,fi%input%jspins))

    IF (fi%input%l_wann.OR.l_socvec) THEN
       wannierspin = 2
    ELSE
       wannierspin = fi%input%jspins
    ENDIF

    !
    !---> set up and solve the eigenvalue problem
    ! --->    radial k-idp s-o matrix elements calc. and storage
    !
#if defined(CPP_MPI)
    !RMA synchronization
    CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
    CALL timestart("eigenso: spnorb")
    !  ..

    !Get spin-orbit coupling matrix elements
    CALL spnorb( fi%atoms,fi%noco,nococonv,fi%input,fmpi, enpara,vTot%mt,usdus,rsoc,.TRUE.,hub1inp,hub1data)
    !


    ALLOCATE (eig_so(2*fi%input%neig))
    ALLOCATE (eigBuffer(2*fi%input%neig,fi%kpts%nkpt,wannierspin))
    ALLOCATE (neigBuffer(fi%kpts%nkpt,wannierspin))
    results%eig = 1.0e300
    eigBuffer = 1.0e300
    unfoldingBuffer = CMPLX(0.0,0.0)
    results%neig = 0
    neigBuffer = 0
    rsoc%soangl(:,:,:,:,:,:) = CONJG(rsoc%soangl(:,:,:,:,:,:))
    CALL timestop("eigenso: spnorb")
    !
    !--->    loop over k-points: each can be a separate task
    DO nk_i=1,SIZE(fmpi%k_list)
        nk=fmpi%k_list(nk_i)
     !DO nk = fmpi%n_start,n_end,n_stride
       CALL lapw%init(fi%input,fi%noco, nococonv,fi%kpts,fi%atoms,fi%sym,nk,fi%cell,.FALSE., fmpi)
       ALLOCATE( zso(lapw%nv(1)+fi%atoms%nlotot,2*fi%input%neig,wannierspin))
       zso(:,:,:) = CMPLX(0.0,0.0)

       CALL timestart("eigenso: alineso")
       CALL alineso(eig_id,lapw, fmpi,fi%atoms,fi%sym,fi%kpts,&
       fi%input,fi%noco,fi%cell,fi%oneD,nk,usdus,rsoc,nsz,nmat, eig_so,zso)
       CALL timestop("eigenso: alineso")
       IF (fmpi%irank.EQ.0) THEN
          WRITE (oUnit,FMT=8010) nk,nsz
          WRITE (oUnit,FMT=8020) (eig_so(i),i=1,nsz)
       ENDIF
8010   FORMAT (1x,/,/,' #k=',i6,':',/,' the',i4,' SOC eigenvalues are:')
8020   FORMAT (5x,5f12.6)

       IF (fmpi%n_rank==0) THEN
          IF (fi%input%eonly) THEN
             CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz))
             STOP 'jspin is undefined here (eigenso - eonly branch)'
             eigBuffer(:nsz,nk,jspin) = eig_so(:nsz)
             neigBuffer(nk,jspin) = nsz
          ELSE
             CALL zmat%alloc(.FALSE.,SIZE(zso,1),nsz)
             DO jspin = 1,wannierspin
                CALL timestart("eigenso: write_eig")

                call timestart("cpy zmat")
                zmat%data_c=zso(:,:nsz,jspin)
                call timestop("cpy zmat")
                
                CALL write_eig(eig_id, nk,jspin,neig=nsz,neig_total=nsz, eig=eig_so(:nsz),zmat=zmat)

                call timestart("cpy buffers")
                eigBuffer(:nsz,nk,jspin) = eig_so(:nsz)
                neigBuffer(nk,jspin) = nsz
                call timestop("cpy buffers")

                CALL timestop("eigenso: write_eig")
             ENDDO
          ENDIF ! (input%eonly) ELSE
       ENDIF ! n_rank == 0
      IF (fi%banddos%unfoldband) THEN
        !IF(modulo (fi%kpts%nkpt,fmpi%n_size).NE.0) call !juDFT_error("number fi%kpts needs to be multiple of number fmpi threads", &
        !                hint=errmsg, calledby="eigenso.F90")
        !write(*,*) 'unfodling for SOC - remember to use useOlap=F'
        jsp=1  
        CALL calculate_plot_w_n(fi%banddos,fi%cell,fi%kpts,zMat,lapw,nk,jsp,eig_so(:nsz),results,fi%input,fi%atoms,unfoldingBuffer,fmpi,fi%noco%l_soc,zso=zso)
        IF (fi%input%jspins==2) THEN
          jsp=2
          CALL calculate_plot_w_n(fi%banddos,fi%cell,fi%kpts,zMat,lapw,nk,jsp,eig_so(:nsz),results,fi%input,fi%atoms,unfoldingBuffer,fmpi,fi%noco%l_soc,zso=zso)
        ENDIF
       END IF
      DEALLOCATE (zso)
    ENDDO ! DO nk

#ifdef CPP_MPI
    IF (fi%banddos%unfoldband) THEN
        results%unfolding_weights = CMPLX(0.0,0.0)
        CALL MPI_ALLREDUCE(unfoldingBuffer,results%unfolding_weights,SIZE(results%unfolding_weights,1)*SIZE(results%unfolding_weights,2)*SIZE(results%unfolding_weights,3),MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
    END IF
    CALL MPI_ALLREDUCE(neigBuffer,results%neig,fi%kpts%nkpt*wannierspin,MPI_INTEGER,MPI_SUM,fmpi%mpi_comm,ierr)
    CALL MPI_ALLREDUCE(eigBuffer(:2*fi%input%neig,:,:),results%eig(:2*fi%input%neig,:,:),&
                       2*fi%input%neig*fi%kpts%nkpt*wannierspin,MPI_DOUBLE_PRECISION,MPI_MIN,fmpi%mpi_comm,ierr)
    CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#else
    results%unfolding_weights(:,:,:) = unfoldingBuffer(:,:,:)
    results%neig(:,:) = neigBuffer(:,:)
    results%eig(:2*fi%input%neig,:,:) = eigBuffer(:2*fi%input%neig,:,:)
#endif

    RETURN
  END SUBROUTINE eigenso
END MODULE m_eigenso
