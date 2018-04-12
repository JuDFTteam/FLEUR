MODULE m_cdnval
  use m_juDFT
CONTAINS
  SUBROUTINE cdnval(eig_id, mpi,kpts,jspin,sliceplot,noco, input,banddos,cell,atoms,enpara,stars,&
       vacuum,dimension,sphhar,sym,obsolete,vTot,oneD,coreSpecInput,den,results,&
       qvac,qvlay,qa21, chmom,clmom)
    !
    !     ***********************************************************
    !         this subroutin is a modified version of cdnval.F.
    !         it calculates a layer charge distribution and an orbital
    !         composition for each state in a film.
    !         this information is written in the  'ek_orbcomp'  file
    !                                    Yu.Koroteev  01.12.2003
    !     ***********************************************************
    !     flapw7 valence density generator
    !                                         c.l.fu
    !     zelec used to calculate ne - 6.12.95 r.pentcheva
    !
    !     changed subroutine to allow parallel writing of vacdos&dosinp
    !     used temporary direct access file 84,tmp_dos to store data used
    !     in cdninf
    !     call of cdninf changed, sympsi is called from cdnval now
    !     look for 'ifdef CPP_MPI' -blocks!               d.wortmann 6.5.99
    !
    !******** ABBREVIATIONS ************************************************
    !     nbands   : number of bands in the energy window
    !     noccbd   : number of occupied bands
    !     slice    : when set to .true. the charge density of a enery range
    !                (slice) or a single state is calculated
    !     e1s,e2s  : (used if slice) energy range for the slice. if both
    !                are set to 0.0 the charge density of the band nr. nnne
    !                is calculated
    !     pallst   : (used if slice) if set to .true. bands above the
    !                Fermi-Energy are taken into account
    !     nnne     : (used if slice) number of the band of which the charge
    !                density is calculated if e1s and e2s are 0.0
    !     kk       : (used if slice) if set to 0 all k-points are used to
    !                calculate the charge density of the slice, otherwise
    !                only k-points kk is taken into account
    !     nslibd   : number of bands in slice
    !     ener     : band energy averaged over all bands and k-points,
    !                wheighted with the l-like charge of each atom type
    !     sqal     : l-like charge of each atom type. sum over all k-points
    !                and bands
    !***********************************************************************
    !
    USE m_constants
    USE m_eig66_io,ONLY: write_dos
    USE m_radfun
    USE m_radflo
    USE m_rhomt
    USE m_rhonmt
    USE m_rhomtlo
    USE m_rhonmtlo
    USE m_mcdinit
    USE m_sympsi
    USE m_eparas      ! energy parameters and partial charges
    USE m_qal21       ! off-diagonal part of partial charges
    USE m_abcof
    USE m_nmat        ! calculate density matrix for LDA + U
    USE m_vacden
    USE m_nstm3
    USE m_pwden
    USE m_forcea8
    USE m_forcea12
    USE m_forcea21
    USE m_checkdopall
    USE m_int21       ! integrate (spin) off-diagonal radial functions
    USE m_int21lo     ! -"- for u_lo
    USE m_rhomt21     ! calculate (spin) off-diagonal MT-density coeff's
    USE m_rhonmt21    ! -"-                       non-MT-density coeff's
    USE m_cdnmt       ! calculate the density and orbital moments etc.
    USE m_orbmom      ! coeffd for orbital moments
    USE m_qmtsl       ! These subroutines divide the input%film into vacuum%layers
    USE m_qintsl      ! (slabs) and intergate the DOS in these vacuum%layers
    USE m_orbcomp     ! calculate corbital composition (like p_x,p_y,p_z)
    USE m_Ekwritesl   ! and write to file.
    USE m_abcrot2
    USE m_doswrite
    USE m_cdnread, ONLY : cdn_read0, cdn_read
    USE m_corespec, only : l_cs    ! calculation of core spectra (EELS)
    USE m_corespec_io, only : corespec_init
    USE m_corespec_eval, only : corespec_gaunt,corespec_rme,corespec_dos,corespec_ddscs
#ifdef CPP_MPI
    USE m_mpi_col_den ! collect density data from parallel nodes
#endif
    USE m_types
    USE m_xmlOutput
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT)    :: results
    TYPE(t_mpi),INTENT(IN)           :: mpi
    TYPE(t_dimension),INTENT(IN)     :: dimension
    TYPE(t_oneD),INTENT(IN)          :: oneD
    TYPE(t_enpara),INTENT(INOUT)     :: enpara
    TYPE(t_obsolete),INTENT(IN)      :: obsolete
    TYPE(t_banddos),INTENT(IN)       :: banddos
    TYPE(t_sliceplot),INTENT(IN)     :: sliceplot
    TYPE(t_input),INTENT(IN)         :: input
    TYPE(t_vacuum),INTENT(IN)        :: vacuum
    TYPE(t_noco),INTENT(IN)          :: noco
    TYPE(t_sym),INTENT(IN)           :: sym
    TYPE(t_stars),INTENT(IN)         :: stars
    TYPE(t_cell),INTENT(IN)          :: cell
    TYPE(t_kpts),INTENT(IN)          :: kpts
    TYPE(t_sphhar),INTENT(IN)        :: sphhar
    TYPE(t_atoms),INTENT(IN)         :: atoms
    TYPE(t_coreSpecInput),INTENT(IN) :: coreSpecInput
    TYPE(t_potden),INTENT(IN)        :: vTot
    TYPE(t_potden),INTENT(INOUT)     :: den

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id,jspin

    !     .. Array Arguments ..
    COMPLEX, INTENT(INOUT) :: qa21(atoms%ntype)
    REAL, INTENT   (OUT) :: chmom(atoms%ntype,dimension%jspd),clmom(3,atoms%ntype,dimension%jspd)
    REAL, INTENT (INOUT) :: qvac(dimension%neigd,2,kpts%nkpt,dimension%jspd)
    REAL, INTENT (INOUT) :: qvlay(dimension%neigd,vacuum%layerd,2,kpts%nkpt,dimension%jspd)

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    logical :: mpi_flag, mpi_status
#endif
    !     .. Local Scalars ..
    TYPE(t_lapw):: lapw
    INTEGER :: llpd
    REAL wronk,emcd_lo,emcd_up
    INTEGER i,ie,iv,ivac,j,k,l,n,ilo,isp,&
         nbands,noded,nodeu,noccbd,nslibd,na,&
         ikpt,jsp_start,jsp_end,ispin
    INTEGER  skip_t,skip_tt
    INTEGER n_size,i_rec,n_rank ,ncored,n_start,n_end,noccbd_l,nbasfcn
    LOGICAL l_fmpl,l_mcd,l_evp,l_orbcomprot,l_real
    !     ...Local Arrays ..
    INTEGER n_bands(0:dimension%neigd),ncore(atoms%ntype)
    REAL    e_mcd(atoms%ntype,input%jspins,dimension%nstd)
    REAL    eig(dimension%neigd)
    REAL    vz0(2)
    REAL    uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
    REAL    ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)

    !orbcomp
    REAL,    ALLOCATABLE :: orbcomp(:,:,:),qmtp(:,:)

    REAL,    ALLOCATABLE :: qis(:,:,:)
    !-dw
    INTEGER, ALLOCATABLE :: gvac1d(:),gvac2d(:)
    INTEGER, ALLOCATABLE :: jsym(:),ksym(:)

    REAL,    ALLOCATABLE :: we(:)

    ! radial functions
    REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:)

    REAL,    ALLOCATABLE :: sqlo(:,:,:)
    REAL,    ALLOCATABLE :: qal(:,:,:,:),sqal(:,:,:),ener(:,:,:)
    REAL,    ALLOCATABLE :: svac(:,:),pvac(:,:),mcd(:,:,:)
    REAL,    ALLOCATABLE :: enerlo(:,:,:),qmat(:,:,:,:)

    COMPLEX, ALLOCATABLE :: qstars(:,:,:,:),m_mcd(:,:,:,:)

    TYPE (t_orb)              :: orb
    TYPE (t_denCoeffs)        :: denCoeffs
    TYPE (t_denCoeffsOffdiag) :: denCoeffsOffdiag
    TYPE (t_force)            :: force
    TYPE (t_slab)             :: slab
    TYPE (t_eigVecCoeffs)     :: eigVecCoeffs

    TYPE (t_usdus)             :: usdus
    TYPE (t_zMat)              :: zMat
    INTEGER :: nkpt_extended

    l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco)
    zmat%l_real=sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco)
    !     ..
    !     ..
    llpd=(atoms%lmaxd*(atoms%lmaxd+3))/2
    !---> l_fmpl is meant as a switch to to a plot of the full magnet.
    !---> density without the atomic sphere approximation for the magnet.
    !---> density. It is not completely implemented (lo's missing).
    l_fmpl = .false.
    IF (noco%l_mperp) THEN
       !--->    when the off-diag. part of the desinsity matrix, i.e. m_x and
       !--->    m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
       !--->    is called only once. therefore, several spin loops have been
       !--->    added. if l_mperp = F, these loops run only from jspin - jspin.
       jsp_start = 1
       jsp_end   = 2
       ALLOCATE ( qmat(0:3,atoms%ntype,dimension%neigd,4) )
    ELSE
       jsp_start = jspin
       jsp_end   = jspin
       ALLOCATE (qmat(1,1,1,1) )
    ENDIF
    !---> if l_mperp = F, these variables are only needed for one spin
    !---> at a time, otherwise for both spins:
    ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )      ! Deallocation before mpi_col_den
    ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
    ALLOCATE ( jsym(dimension%neigd),ksym(dimension%neigd) )
    ALLOCATE ( gvac1d(dimension%nv2d),gvac2d(dimension%nv2d) )
    ALLOCATE ( qal(0:3,atoms%ntype,dimension%neigd,jsp_start:jsp_end) )
    ALLOCATE ( sqal(0:3,atoms%ntype,jsp_start:jsp_end) )
    ALLOCATE ( ener(0:3,atoms%ntype,jsp_start:jsp_end) )
    ALLOCATE (   sqlo(atoms%nlod,atoms%ntype,jsp_start:jsp_end) )
    ALLOCATE ( enerlo(atoms%nlod,atoms%ntype,jsp_start:jsp_end) )
    ALLOCATE ( svac(2,jsp_start:jsp_end) )
    ALLOCATE ( pvac(2,jsp_start:jsp_end) )
    ALLOCATE ( qstars(vacuum%nstars,dimension%neigd,vacuum%layerd,2) )

    ! --> Initializations

    CALL usdus%init(atoms,input%jspins)
    CALL denCoeffs%init(atoms,sphhar,jsp_start,jsp_end)
    CALL denCoeffsOffdiag%init(atoms,noco,sphhar,l_fmpl)
    CALL force%init1(input,atoms)

    IF ((l_fmpl).AND.(.not.noco%l_mperp)) CALL juDFT_error("for fmpl set noco%l_mperp = T!" ,calledby ="cdnval")

    svac(:,:) = 0.0 ; pvac(:,:) = 0.0
    sqal(:,:,:) = 0.0 ; ener(:,:,:) = 0.0

    CALL orb%init(atoms,noco,jsp_start,jsp_end)

    INQUIRE (file='mcd_inp',exist=l_mcd)
    IF (l_mcd) THEN
       OPEN (23,file='mcd_inp',STATUS='old',FORM='formatted')
       READ (23,*) emcd_lo,emcd_up
       CLOSE (23)
       ALLOCATE ( m_mcd(dimension%nstd,(3+1)**2,3*atoms%ntype,2) )
       ALLOCATE ( mcd(3*atoms%ntype,dimension%nstd,dimension%neigd) )
       IF (.not.banddos%dos) WRITE (*,*) 'For mcd-spectra set banddos%dos=T!'
    ELSE
       ALLOCATE ( m_mcd(1,1,1,1),mcd(1,1,1) )
    ENDIF

! calculation of core spectra (EELS) initializations -start-
    CALL corespec_init(input,atoms,coreSpecInput)
    IF(l_cs.AND.(mpi%isize.NE.1)) CALL juDFT_error('EELS + MPI not implemented', calledby = 'cdnval')
    IF(l_cs.AND.jspin.EQ.1) CALL corespec_gaunt()
! calculation of core spectra (EELS) initializations -end-

  
    IF (mpi%irank==0) THEN
       WRITE (6,FMT=8000) jspin
       WRITE (16,FMT=8000) jspin
       CALL openXMLElementPoly('mtCharges',(/'spin'/),(/jspin/))
    END IF
8000 FORMAT (/,/,10x,'valence density: spin=',i2)

    CALL cdn_read0(eig_id,mpi%irank,mpi%isize,jspin,dimension%jspd,&
                   noco%l_noco,n_bands,n_size)
#ifdef CPP_MPI
    ! Sinchronizes the RMA operations
    CALL MPI_BARRIER(mpi%mpi_comm,ie) 
#endif

    ALLOCATE ( qis(dimension%neigd,kpts%nkpt,dimension%jspd))

    skip_tt = dot_product(enpara%skiplo(:atoms%ntype,jspin),atoms%neq(:atoms%ntype))
    IF (noco%l_soc.OR.noco%l_noco)  skip_tt = 2 * skip_tt
    !-lo
    !---> set up l-dependent m.t. wavefunctions
    na = 1
    ncored = 0

    ALLOCATE ( flo(atoms%jmtd,2,atoms%nlod,dimension%jspd) )
    DO n = 1,atoms%ntype
       IF (input%cdinf.AND.mpi%irank==0) WRITE (6,FMT=8001) n
       DO  l = 0,atoms%lmax(n)
          DO ispin =jsp_start,jsp_end
             CALL radfun(l,n,ispin,enpara%el0(l,n,ispin),vTot%mt(1,0,n,ispin),atoms,&
                         f(1,1,l,ispin),g(1,1,l,ispin),usdus,nodeu,noded,wronk)
             IF (input%cdinf.AND.mpi%irank==0) WRITE (6,FMT=8002) l,&
                  enpara%el0(l,n,ispin),usdus%us(l,n,ispin),usdus%dus(l,n,ispin),nodeu,&
                  usdus%uds(l,n,ispin),usdus%duds(l,n,ispin),noded,usdus%ddn(l,n,ispin),&
                  wronk
          END DO
          IF (noco%l_mperp) THEN
             CALL int_21(f,g,atoms,n,l,denCoeffsOffdiag)
          END IF
       END DO
       IF (l_mcd) THEN
          CALL mcd_init(atoms,input,dimension,&
                        vTot%mt(:,0,:,:),g,f,emcd_up,emcd_lo,n,jspin,&
                        ncore,e_mcd,m_mcd)
          ncored = max(ncore(n),ncored)
       END IF

       IF(l_cs) CALL corespec_rme(atoms,input,n,dimension%nstd,&
                                  input%jspins,jspin,results%ef,&
                                  dimension%msh,vTot%mt(:,0,:,:),f,g)

       !
       !--->   generate the extra wavefunctions for the local orbitals,
       !--->   if there are any.
       !
       IF ( atoms%nlo(n) > 0 ) THEN
          DO ispin = jsp_start,jsp_end
             CALL radflo(atoms,n,ispin, enpara%ello0(1,1,ispin),vTot%mt(:,0,n,ispin), f(1,1,0,ispin),&
                         g(1,1,0,ispin),mpi, usdus, uuilon,duilon,ulouilopn, flo(:,:,:,ispin))
          END DO
       END IF

       DO ilo = 1, atoms%nlo(n)
          IF (noco%l_mperp) THEN
             CALL int_21lo(f,g,atoms,n,flo,ilo,denCoeffsOffdiag)
          END IF
       END DO

       na = na + atoms%neq(n)
    END DO
    DEALLOCATE (flo)
8001 FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',/,&
         t32,'radial function',t79,'energy derivative',/,t3,'l',t8,&
         'energy',t26,'value',t39,'derivative',t53,'nodes',t68,&
         'value',t81,'derivative',t95,'nodes',t107,'norm',t119,&
         'wronskian')
8002 FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

    IF (input%film) vz0(:) = vTot%vacz(vacuum%nmz,:,jspin)

    CALL slab%init(banddos,dimension,atoms,cell)

    IF ((banddos%ndir.EQ.-3).AND.banddos%dos) THEN
       IF (oneD%odi%d1)  CALL juDFT_error("layer-resolved feature does not work with 1D",calledby ="cdnval")

       ALLOCATE ( orbcomp(dimension%neigd,23,atoms%nat) )
       ALLOCATE ( qmtp(dimension%neigd,atoms%nat) )
       IF (.NOT.input%film) qvac(:,:,:,jspin) = 0.0
    ELSE
       ALLOCATE(orbcomp(1,1,1),qmtp(1,1))
    END IF

    !-->   loop over k-points: each can be a separate task
    IF (kpts%nkpt < mpi%isize) THEN
       l_evp = .true.
       IF (l_mcd) THEN
          mcd(:,:,:) = 0.0
       ENDIF
       ener(:,:,:) = 0.0
       sqal(:,:,:) = 0.0
       qal(:,:,:,:) = 0.0
       enerlo(:,:,:) = 0.0
       sqlo(:,:,:) = 0.0
    ELSE
       l_evp = .false.
    END IF
    ALLOCATE ( we(dimension%neigd) )
    i_rec = 0 ; n_rank = 0

    ! For k-point paralelization: 
    ! the number of iterations is adjusted to the number of MPI processes to synchronize RMA operations
    if (l_evp) then
       nkpt_extended = kpts%nkpt
    else
       nkpt_extended = (kpts%nkpt / mpi%isize + 1) * mpi%isize
    endif
    DO ikpt = 1,nkpt_extended
       i_rec = i_rec + 1
       IF ((mod(i_rec-1,mpi%isize).EQ.mpi%irank).OR.l_evp) THEN
        IF ( ikpt < kpts%nkpt + 1) THEN
          !-t3e
          we=0.0
          !--->    determine number of occupied bands and set weights (we)
          noccbd = 0
          DO i = 1,dimension%neigd ! nbands
             we(i) = results%w_iks(n_bands(n_rank)+i,ikpt,jspin)
             IF (noco%l_noco) we(i) = results%w_iks(i,ikpt,1)
             IF ((we(i).GE.1.e-8).OR.input%pallst) THEN
                noccbd = noccbd + 1
             ELSE
                we(i)=0.0
             END IF
          END DO
          ! uncomment this so that cdinf plots works for all states
          ! noccbd = neigd

          !
          ! -> Gu test: distribute ev's among the processors...
          !
          CALL lapw%init(input,noco, kpts,atoms,sym,ikpt,cell,.false., mpi)
          skip_t = skip_tt
          IF (l_evp.AND.(mpi%isize.GT.1)) THEN
             IF (banddos%dos) THEN
                noccbd_l = CEILING( real(n_bands(1)) / mpi%isize )
                n_start = mpi%irank*noccbd_l + 1
                n_end   = min( (mpi%irank+1)*noccbd_l , n_bands(1) )
             ELSE
                noccbd_l = CEILING( real(noccbd) / mpi%isize )
                n_start = mpi%irank*noccbd_l + 1
                n_end   = min( (mpi%irank+1)*noccbd_l , noccbd )
             END IF
             noccbd = n_end - n_start + 1
             IF (noccbd<1) THEN
                noccbd=0
             ELSE
                we(1:noccbd) = we(n_start:n_end)
             END IF
             IF (n_start > skip_tt) THEN
                skip_t  = 0
             END IF
             IF (n_end <= skip_tt) THEN
                skip_t  = noccbd
             END IF
             IF ((n_start <= skip_tt).AND.(n_end > skip_tt)) THEN
                skip_t  = mod(skip_tt,noccbd)
             END IF
          ELSE
             n_start = 1
             IF (banddos%dos) THEN
                noccbd_l = n_bands(1)
                n_end    = n_bands(1)
                noccbd   = n_bands(1)
             ELSE
                noccbd_l = noccbd
                n_end    = noccbd
             END IF
          END IF

          nbasfcn = MERGE(zMat%nbasfcn+lapw%nv(2),lapw%nv(1),noco%l_noco)+atoms%nlotot
          CALL zMat%init(l_real,nbasfcn,dimension%neigd)
          CALL cdn_read(eig_id,dimension%nvd,dimension%jspd,mpi%irank,mpi%isize,&
                        ikpt,jspin,zmat%nbasfcn,noco%l_ss,noco%l_noco,&
                        noccbd,n_start,n_end,nbands,eig,zMat)
#ifdef CPP_MPI
          ! Sinchronizes the RMA operations
          CALL MPI_BARRIER(mpi%mpi_comm,ie)
#endif
          !IF (l_evp.AND.(isize.GT.1)) THEN
          !  eig(1:noccbd) = eig(n_start:n_end)
          !ENDIF
          !
          IF (vacuum%nstm.EQ.3.AND.input%film) THEN
             CALL nstm3(sym,atoms,vacuum,stars,ikpt,lapw%nv(jspin),&
                        input,jspin,kpts,&
                        cell,kpts%wtkpt(ikpt),lapw%k1(:,jspin),lapw%k2(:,jspin),&
                        enpara%evac0(1,jspin),vTot%vacz(:,:,jspin),vz0,&
                        gvac1d,gvac2d)
          END IF

          IF (noccbd.EQ.0) GO TO 199
          !
          !--->    if slice, only a certain bands are taken into account
          !--->    in order to do this the coresponding eigenvalues, eigenvectors
          !--->    and weights have to be copied to the beginning of the arrays
          !--->    eig, z and we and the number of occupied bands (noccbd) has to
          !--->    changed
          IF (sliceplot%slice) THEN
             IF (mpi%irank==0) WRITE (16,FMT=*) 'NNNE',sliceplot%nnne
             IF (mpi%irank==0) WRITE (16,FMT=*) 'sliceplot%kk',sliceplot%kk
             nslibd = 0
             IF (input%pallst) we(:nbands) = kpts%wtkpt(ikpt)
             IF (sliceplot%kk.EQ.0) THEN
                IF (mpi%irank==0) THEN
                   WRITE (16,FMT='(a)') 'ALL K-POINTS ARE TAKEN IN SLICE'
                   WRITE (16,FMT='(a,i2)') ' sliceplot%slice: k-point nr.',ikpt
                END IF
                DO i = 1,nbands
                   IF (eig(i).GE.sliceplot%e1s .AND. eig(i).LE.sliceplot%e2s) THEN
                      nslibd = nslibd + 1
                      eig(nslibd) = eig(i)
                      we(nslibd) = we(i)
                      if (zmat%l_real) THEN
                         zMat%z_r(:,nslibd) = zMat%z_r(:,i)
                      else
                         zMat%z_c(:,nslibd) = zMat%z_c(:,i)
                      endif
                   END IF
                END DO
                IF (mpi%irank==0) WRITE (16,'(a,i3)') ' eigenvalues in sliceplot%slice:',nslibd
             ELSE IF (sliceplot%kk.EQ.ikpt) THEN
                IF (mpi%irank==0) WRITE (16,FMT='(a,i2)') ' sliceplot%slice: k-point nr.',ikpt
                IF ((sliceplot%e1s.EQ.0.0) .AND. (sliceplot%e2s.EQ.0.0)) THEN
                   IF (mpi%irank==0) WRITE (16,FMT='(a,i5,f10.5)') 'slice: eigenvalue nr.',&
                        sliceplot%nnne,eig(sliceplot%nnne)
                   nslibd = nslibd + 1
                   eig(nslibd) = eig(sliceplot%nnne)
                   we(nslibd) = we(sliceplot%nnne)
                   if (zmat%l_real) Then
                      zMat%z_r(:,nslibd) = zMat%z_r(:,sliceplot%nnne)
                   else
                      zMat%z_c(:,nslibd) = zMat%z_c(:,sliceplot%nnne)
                   endif
                ELSE
                   DO i = 1,nbands
                      IF (eig(i).GE.sliceplot%e1s .AND. eig(i).LE.sliceplot%e2s) THEN
                         nslibd = nslibd + 1
                         eig(nslibd) = eig(i)
                         we(nslibd) = we(i)
                         if (zmat%l_real) THEN
                            zMat%z_r(:,nslibd) = zMat%z_r(:,i)
                         else
                            zMat%z_c(:,nslibd) = zMat%z_c(:,i)
                         endif
                      END IF
                   END DO
                   IF (mpi%irank==0) WRITE (16,FMT='(a,i3)')' eigenvalues in sliceplot%slice:',nslibd
                END IF
             END IF
             noccbd = nslibd
             IF (nslibd.EQ.0) GO TO 199 !200
          END IF ! sliceplot%slice

          !--->    in normal iterations the charge density of the unoccupied
          !--->    does not need to be calculated (in pwden, vacden and abcof)
          IF (banddos%dos.AND. .NOT.(l_evp.AND.(mpi%isize.GT.1)) ) THEN
             noccbd=nbands
          END IF
          !     ----> add in spin-doubling factor
          we(:noccbd) = 2.0 * we(:noccbd) / input%jspins

          !---> pk non-collinear
          !--->    valence density in the interstitial and vacuum region
          !--->    has to be called only once (if jspin=1) in the non-collinear
          !--->    case
          !     ----> valence density in the interstitial region
          IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
             CALL timestart("cdnval: pwden")
             CALL pwden(stars,kpts,banddos,oneD,input,mpi,noco,cell,atoms,sym,ikpt,&
                        jspin,lapw,noccbd,we,eig,den,qis,results,force%f_b8,zMat)
             CALL timestop("cdnval: pwden")
          END IF
          !+new
          !--->    charge of each valence state in this k-point of the SBZ
          !--->    in the layer interstitial region of the film
          !
          IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
             IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
                CALL q_int_sl(jspin,stars,atoms,sym,cell,noccbd,lapw,slab,oneD,zMat)
             END IF
          END IF
          !-new c
          !--->    valence density in the vacuum region
          IF (input%film) THEN
             IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
                CALL timestart("cdnval: vacden")
                CALL vacden(vacuum,dimension,stars,oneD, kpts,input, cell,atoms,noco,banddos,&
                            gvac1d,gvac2d, we,ikpt,jspin,vTot%vacz(:,:,jspin),vz0, noccbd,lapw, enpara%evac0,eig,&
                            den,qvac,qvlay, qstars,zMat)
                CALL timestop("cdnval: vacden")
             END IF
             !--->       perform Brillouin zone integration and summation over the
             !--->       bands in order to determine the vacuum energy parameters.
             DO ispin = jsp_start,jsp_end
                DO ivac = 1,vacuum%nvac
                   pvac(ivac,ispin)=pvac(ivac,ispin)+dot_product(eig(:noccbd)*qvac(:noccbd,ivac,ikpt,ispin),we(:noccbd))
                   svac(ivac,ispin)=svac(ivac,ispin)+dot_product(qvac(:noccbd,ivac,ikpt,ispin),we(:noccbd))
                END DO
             END DO
          END IF

          !--->    valence density in the atomic spheres
          CALL eigVecCoeffs%init(dimension,atoms,noco,jspin,noccbd)

          DO ispin = jsp_start,jsp_end
             IF (input%l_f) THEN
                CALL timestart("cdnval: to_pulay")
                CALL force%init2(noccbd,input,atoms)
                CALL abcof(input,atoms,sym, cell,lapw,noccbd,usdus, noco,ispin,oneD,&
                           eigVecCoeffs%acof(:,0:,:,ispin),eigVecCoeffs%bcof(:,0:,:,ispin),&
                           eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat,eig,force)
                CALL timestop("cdnval: to_pulay")
             ELSE
                CALL timestart("cdnval: abcof")
                CALL abcof(input,atoms,sym, cell,lapw,noccbd,usdus, noco,ispin,oneD,&
                           eigVecCoeffs%acof(:,0:,:,ispin),eigVecCoeffs%bcof(:,0:,:,ispin),&
                           eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat)
                CALL timestop("cdnval: abcof")

             END IF

             IF (atoms%n_u.GT.0) THEN
                CALL n_mat(atoms,sym,noccbd,usdus,ispin,we,eigVecCoeffs,den%mmpMat(:,:,:,jspin))
             END IF

             !--->       perform Brillouin zone integration and summation over the
             !--->       bands in order to determine the energy parameters for each
             !--->       atom and angular momentum
             IF (.not.sliceplot%slice) THEN
                CALL eparas(ispin,atoms,noccbd,mpi,ikpt,noccbd,we,eig,&
                            skip_t,l_evp,eigVecCoeffs,usdus,&
                            ncore,l_mcd,m_mcd,enerlo(1,1,ispin),sqlo(1,1,ispin),&
                            ener(0,1,ispin),sqal(0,1,ispin),&
                            qal(0:,:,:,ispin),mcd)

                IF (noco%l_mperp.AND.(ispin == jsp_end)) THEN
                   CALL qal_21(atoms,input,noccbd,we,noco,eigVecCoeffs,denCoeffsOffdiag,qal,qmat)
                END IF
             END IF
             !
             !+new
             !--->    layer charge of each valence state in this k-point of the SBZ
             !--->    from the mt-sphere region of the film
             !
             IF (banddos%dos.AND.(banddos%ndir.EQ.-3))  THEN
                CALL q_mt_sl(ispin,atoms,noccbd,ikpt,noccbd,skip_t,noccbd,eigVecCoeffs,usdus,slab)

                INQUIRE (file='orbcomprot',exist=l_orbcomprot)
                IF (l_orbcomprot) THEN                           ! rotate ab-coeffs
                   CALL abcrot2(atoms,noccbd,eigVecCoeffs,ispin)
                END IF

                CALL orb_comp(ispin,noccbd,atoms,noccbd,usdus,eigVecCoeffs,orbcomp,qmtp)
             END IF
             !-new
             !--->          set up coefficients for the spherical and
             CALL timestart("cdnval: rhomt")
             CALL rhomt(atoms,we,noccbd,eigVecCoeffs,denCoeffs,ispin)
             CALL timestop("cdnval: rhomt")

             IF (noco%l_soc) CALL orbmom(atoms,noccbd,we,ispin,eigVecCoeffs,orb)
             !--->          non-spherical m.t. density
             CALL timestart("cdnval: rhonmt")
             CALL rhonmt(atoms,sphhar,we,noccbd,sym,eigVecCoeffs,denCoeffs,ispin)
             CALL timestop("cdnval: rhonmt")

             !--->          set up coefficients of the local orbitals and the
             !--->          flapw - lo cross terms for the spherical and
             !--->          non-spherical mt density
             CALL timestart("cdnval: rho(n)mtlo")
             CALL rhomtlo(atoms,noccbd,we,eigVecCoeffs,denCoeffs,ispin)

             CALL rhonmtlo(atoms,sphhar,noccbd,we,eigVecCoeffs,denCoeffs,ispin)
             CALL timestop("cdnval: rho(n)mtlo")

             IF (input%l_f) THEN
                CALL timestart("cdnval: force_a12/21")
                IF (.not.input%l_useapw) THEN
                   CALL force_a12(atoms,noccbd,sym,dimension,cell,oneD,&
                                  we,ispin,noccbd,usdus,eigVecCoeffs,force,results)
                ENDIF
                CALL force_a21(input,atoms,dimension,noccbd,sym,oneD,cell,we,ispin,&
                               enpara%el0(0:,:,ispin),noccbd,eig,usdus,eigVecCoeffs,force,results)
                CALL timestop("cdnval: force_a12/21")
             END IF

             IF(l_cs) THEN
                CALL corespec_dos(atoms,usdus,ispin,dimension%lmd,kpts%nkpt,ikpt,dimension%neigd,&
                                  noccbd,results%ef,banddos%sig_dos,eig,we,eigVecCoeffs)
             END IF
          END DO !--->    end loop over ispin

          IF (noco%l_mperp) THEN
             CALL rhomt21(atoms,we,noccbd,eigVecCoeffs,denCoeffsOffdiag)
             IF (l_fmpl) THEN
                CALL rhonmt21(atoms,sphhar,we,noccbd,sym,eigVecCoeffs,denCoeffsOffdiag)
             END IF
          END IF

199       CONTINUE
          IF ((banddos%dos .OR. banddos%vacdos .OR. input%cdinf)  ) THEN
             CALL timestart("cdnval: write_info")
             !--->    calculate charge distribution of each state (l-character ...)
             !--->    and write the information to the files dosinp and vacdos
             !--->    for dos and bandstructure plots

             !--dw    parallel writing of vacdos,dosinp....
             !        write data to direct access file first, write to formated file later by PE 0 only!
             !--dw    since z is no longer an argument of cdninf sympsi has to be called here!

             IF (banddos%ndir.GT.0) THEN
                CALL sympsi(lapw%bkpt,lapw%nv(jspin),lapw%k1(:,jspin),lapw%k2(:,jspin),&
                            lapw%k3(:,jspin),sym,dimension,nbands,cell,eig,noco, ksym,jsym,zMat)
             END IF

             !--dw   now write k-point data to tmp_dos
             CALL write_dos(eig_id,ikpt,jspin,qal(:,:,:,jspin),qvac(:,:,ikpt,jspin),qis(:,ikpt,jspin),&
                            qvlay(:,:,:,ikpt,jspin),qstars,ksym,jsym,mcd,slab%qintsl,&
                            slab%qmtsl(:,:),qmtp(:,:),orbcomp)

             CALL timestop("cdnval: write_info")
             !-new_sl
          END IF

          !--->  end of loop over PE's
        ELSE !(ikpt < nkpt + 1)
#ifdef CPP_MPI
          ! Synchronizes the RMA operations
          CALL MPI_BARRIER(mpi%mpi_comm,ie)
#endif
        END IF
       END IF ! --> end "IF ((mod(i_rec-1,mpi%isize).EQ.mpi%irank).OR.l_evp) THEN"
    END DO !---> end of k-point loop
    DEALLOCATE (we,f,g)
    !+t3e
#ifdef CPP_MPI
    CALL timestart("cdnval: mpi_col_den")
    DO ispin = jsp_start,jsp_end
       CALL mpi_col_den(mpi,sphhar,atoms,oneD,stars,vacuum,&
                        input,noco,l_fmpl,ispin,llpd, den%vacxy(1,1,1,ispin),&
                        den%vacz(1,1,ispin),den%pw(1,ispin), ener(0,1,ispin),sqal(0,1,ispin),&
                        results,svac(1,ispin),pvac(1,ispin),denCoeffs,&
                        sqlo(1,1,ispin),enerlo(1,1,ispin),orb,&
                        denCoeffsOffdiag,den,den%mmpMat(:,:,:,jspin))
    END DO
    CALL timestop("cdnval: mpi_col_den")
#endif

    IF(l_cs) CALL corespec_ddscs(jspin,input%jspins)

    IF (((jspin.eq.input%jspins).OR.noco%l_mperp) .AND. (banddos%dos.or.banddos%vacdos.or.input%cdinf) ) THEN
       CALL timestart("cdnval: dos")
       IF (mpi%irank==0) THEN
          CALL doswrite(eig_id,dimension,kpts,atoms,vacuum,input,banddos,&
                        sliceplot,noco,sym,cell,l_mcd,ncored,ncore,e_mcd,&
                        results%ef,results%bandgap,slab%nsld,oneD)
          IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
             CALL Ek_write_sl(eig_id,dimension,kpts,atoms,vacuum,input,jspin,sym,cell,slab)
          END IF
       END IF
#ifdef CPP_MPI                
       CALL MPI_BARRIER(mpi%mpi_comm,ie)
#endif
       CALL timestop("cdnval: dos")
    END IF

    IF (mpi%irank==0) THEN
       CALL cdnmt(dimension%jspd,atoms,sphhar,llpd,&
                  noco,l_fmpl,jsp_start,jsp_end,&
                  enpara%el0,enpara%ello0,vTot%mt(:,0,:,:),denCoeffs,&
                  usdus,orb,denCoeffsOffdiag,&
                  chmom,clmom,qa21,den%mt)

       DO ispin = jsp_start,jsp_end
          IF (.NOT.sliceplot%slice) THEN
             DO n=1,atoms%ntype
                enpara%el1(0:3,n,ispin)=ener(0:3,n,ispin)/sqal(0:3,n,ispin)
                IF (atoms%nlo(n)>0) enpara%ello1(:atoms%nlo(n),n,ispin)=enerlo(:atoms%nlo(n),n,ispin)/sqlo(:atoms%nlo(n),n,ispin)
             END DO
             IF (input%film) enpara%evac1(:vacuum%nvac,ispin)=pvac(:vacuum%nvac,ispin)/svac(:vacuum%nvac,ispin)
          END IF

          !--->      check continuity of charge density
          IF (input%cdinf) THEN
             CALL timestart("cdnval: cdninf-stuff")
             WRITE (6,FMT=8210) ispin
8210         FORMAT (/,5x,'check continuity of cdn for spin=',i2)
             CALL checkDOPAll(input,dimension,sphhar,stars,atoms,sym,vacuum,oneD,&
                              cell,den,ispin)
             CALL timestop("cdnval: cdninf-stuff")
          END IF
          !+for
          !--->      forces of equ. A8 of Yu et al.
          IF ((input%l_f)) THEN
             CALL timestart("cdnval: force_a8")
             CALL force_a8(input,atoms,sphhar,ispin,vTot%mt(:,:,:,ispin),den%mt,force,results)
             CALL timestop("cdnval: force_a8")
          END IF
          !-for
       END DO ! end of loop ispin = jsp_start,jsp_end
       CALL closeXMLElement('mtCharges')

       IF(vacuum%nvac.EQ.1) THEN
          den%vacz(:,2,:) = den%vacz(:,1,:)
          IF (sym%invs) THEN
             den%vacxy(:,:,2,:) = CONJG(den%vacxy(:,:,1,:))
          ELSE
             den%vacxy(:,:,2,:) = den%vacxy(:,:,1,:)
          END IF
       END IF

    END IF ! end of (mpi%irank==0)
    !+t3e
    !Note: no deallocation anymore, we rely on Fortran08 :-)

    IF ((jsp_end.EQ.input%jspins)) THEN
       IF ((banddos%dos.OR.banddos%vacdos).AND.(banddos%ndir/=-2))  CALL juDFT_end("DOS OK",mpi%irank)
       IF (vacuum%nstm.EQ.3)  CALL juDFT_end("VACWAVE OK",mpi%irank)
    END IF

  END SUBROUTINE cdnval
END MODULE m_cdnval
