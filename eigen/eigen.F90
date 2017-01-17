!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen
  use m_juDFT
CONTAINS
  SUBROUTINE eigen(mpi,stars,sphhar,atoms,obsolete,xcpot,&
       sym,kpts,dimension, vacuum, input, cell, enpara_in,banddos, noco,jij, oneD,hybrid,&
       it,eig_id,results)
    !*********************************************************************
    !     sets up and solves the eigenvalue problem for a basis of lapws.
    !
    ! nv,   nvd     ... actual length & dimension of EV without LO's
    ! nmat, nbasfcn                                   including LO's
    !        g. bihlmayer '96
    !**********************************************************************
    USE m_constants, ONLY : pi_const,sfp_const
    USE m_types
    USE m_lodpot
    USE m_tlmplm
    USE m_tlmplm_store
    USE m_apws
    USE m_hsmt
    USE m_hsint
    USE m_hsvac
    USE m_od_hsvac
    USE m_usetup
    USE m_loddop
    USE m_eigen_diag
#ifdef CPP_NOTIMPLEMENTED
    USE m_symm_hf,  ONLY : symm_hf_nkpt_EIBZ
    USE m_gen_bz
    USE m_gen_wavf
    USE m_hsfock
    USE m_read_core
    USE m_subvxc
    USE m_gweig
    USE m_gw_qsgw
    USE m_checkolap
#endif
    USE m_hsefunctional
    USE m_hybridmix    , ONLY: amix_pbe0,amix_hf
    USE m_util
    USE m_icorrkeys
    USE m_eig66_io, ONLY : open_eig, write_eig, close_eig,read_eig
    USE m_xmlOutput

    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: dimension
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_hybrid),INTENT(IN)    :: hybrid
    TYPE(t_enpara),INTENT(INOUT) :: enpara_in
    TYPE(t_obsolete),INTENT(IN)  :: obsolete
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_banddos),INTENT(IN)   :: banddos
    TYPE(t_jij),INTENT(IN)       :: jij
    TYPE(t_sym),INTENT(IN)       :: sym  
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(INOUT)  :: atoms !in u_setup n_u might be modified

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN) :: it
    INTEGER,INTENT(INOUT):: eig_id
    !     ..
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..
    INTEGER jsp,nk,nred,ne_all,n_u_in,ne_found
    INTEGER iter,ne,matsize  ,nrec,lh0
    INTEGER nspins,isp,l,i,j,err,gwc
    INTEGER mlotot,mlolotot,mlot_d,mlolot_d,nlot_d
    LOGICAL l_wu,lcal_qsgw,l_file,l_real,l_zref
    REAL evac_sv(dimension%jspd)
    INTEGER ::eig_id_hf=-1
    INTEGER :: nu=8
    
    !     ..
    !     .. Local Arrays ..
    INTEGER, PARAMETER :: lmaxb=3
    INTEGER, ALLOCATABLE :: matind(:,:),kveclo(:)
    INTEGER, ALLOCATABLE :: nv2(:)
    REAL,    ALLOCATABLE :: bkpt(:)
    REAL,    ALLOCATABLE :: eig(:)

    COMPLEX, ALLOCATABLE :: vpw(:,:),vzxy(:,:,:,:)
    COMPLEX, ALLOCATABLE :: vpwtot(:,:)
    REAL,    ALLOCATABLE :: vz(:,:,:),vr(:,:,:,:)
    REAL,    ALLOCATABLE :: vrtot(:,:,:,:)

    COMPLEX, ALLOCATABLE :: vs_mmp(:,:,:,:)
    TYPE(t_tlmplm)  :: td
    TYPE(t_usdus)   :: ud
    TYPE(t_lapw)    :: lapw
    Type(t_enpara)  :: enpara
    TYPE(t_zMat)    :: zMat
    TYPE(t_hamOvlp) :: hamOvlp
    !
    INTEGER fh,nn,n
    INTEGER ierr(3)

    !
    !     .. variables for HF or hybrid functional calculation ..
    !
    !      - scalar -
#ifdef CPP_NEVER
    INTEGER, INTENT(IN)     ::  maxlcutm,maxindxm,maxbasm
    INTEGER, INTENT(IN)     ::  maxindxp
    INTEGER, INTENT(IN)     ::  bands
    !     - arrays -
    INTEGER, INTENT(IN)     ::  nindxm(0:maxlcutm,atoms%ntype)
    INTEGER, INTENT(IN)     ::  lcutm(atoms%ntype)
    REAL   , INTENT(IN)     ::  basm(atoms%jmtd,maxindxm,0:maxlcutm,atoms%ntype)
#endif
    !     - local scalar -
    INTEGER                 ::  itype,ispin,isym,iisym
    INTEGER                 ::  indx,ic
    INTEGER                 ::  ll,lm,l1,l2
    INTEGER                 ::  lmaxcd
    INTEGER                 ::  maxindxc,mnobd
    INTEGER                 ::  maxfac
    INTEGER                 ::  maxbands
    LOGICAL                 ::  l_hybrid
    !     - local arrays -
#ifdef CPP_NEVER
    INTEGER                 ::  nobd(kpts%nkptf)
    INTEGER                 ::  lmaxc(atoms%ntype)
    INTEGER                 ::  g(3)
    INTEGER                 ::  nindxp(0:maxlcutm,atoms%ntype)
    INTEGER , ALLOCATABLE   ::  nkpt_EIBZ(:)
    INTEGER , ALLOCATABLE   ::  nindxc(:,:)
    INTEGER , ALLOCATABLE   ::  kveclo_eig(:,:)
    INTEGER , ALLOCATABLE   ::  nbasm(:)
    INTEGER                 ::  comm(kpts%nkpt),irank2(kpts%nkpt),isize2(kpts%nkpt)
    REAL                    ::  el_eig(0:atoms%lmaxd,atoms%ntype), ello_eig(atoms%nlod,atoms%ntype),rarr(3)
    REAL                    ::  bas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
    REAL                    ::  drbas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
    REAL,    ALLOCATABLE    ::  eig_c(:,:,:)
    REAL,    ALLOCATABLE    ::  core1(:,:,:,:),core2(:,:,:,:)
    REAL,    ALLOCATABLE    ::  gauntarr(:,:,:,:,:,:)
    REAL,    ALLOCATABLE    ::  sfac(:),fac(:)
    REAL,    ALLOCATABLE    ::  prodm(:,:,:,:)
    TYPE(PRODTYPE),ALLOCATABLE :: prod(:,:,:)
#endif
    INTEGER                 ::  ne_eig(kpts%nkpt),nbands(kpts%nkpt)
    REAL,    ALLOCATABLE    ::  eig_irr(:,:),vr0(:,:,:)
    REAL                    ::  bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
    REAL                    ::  bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)

#ifdef CPP_MPI
    INTEGER   :: sndreqd,sndreq(mpi%isize*kpts%nkpt)
#endif
    !
    !
    ! --> Allocate
    !
    ALLOCATE ( ud%uloulopn(atoms%nlod,atoms%nlod,atoms%ntype,dimension%jspd),nv2(dimension%jspd) )
    ALLOCATE ( ud%ddn(0:atoms%lmaxd,atoms%ntype,dimension%jspd),eig(dimension%neigd),bkpt(3) )
    ALLOCATE ( ud%us(0:atoms%lmaxd,atoms%ntype,dimension%jspd),ud%uds(0:atoms%lmaxd,atoms%ntype,dimension%jspd) )
    ALLOCATE ( ud%dus(0:atoms%lmaxd,atoms%ntype,dimension%jspd),ud%duds(0:atoms%lmaxd,atoms%ntype,dimension%jspd))
    ALLOCATE ( ud%ulos(atoms%nlod,atoms%ntype,dimension%jspd),ud%dulos(atoms%nlod,atoms%ntype,dimension%jspd) )
    ALLOCATE ( ud%uulon(atoms%nlod,atoms%ntype,dimension%jspd),ud%dulon(atoms%nlod,atoms%ntype,dimension%jspd) )
   ! ALLOCATE ( enpara%ello(atoms%nlod,atoms%ntype,dimension%jspd) )
   ! ALLOCATE ( enpara%el(0:atoms%lmaxd,atoms%ntype,dimension%jspd),enpara%evac(2,dimension%jspd) )
    ALLOCATE ( lapw%k1(dimension%nvd,dimension%jspd),lapw%k2(dimension%nvd,dimension%jspd),lapw%k3(dimension%nvd,dimension%jspd),lapw%rk(dimension%nvd,dimension%jspd) )
    !
    ! --> some parameters first
    !
    !     determine the total number of lo's : nlotot
    !
    mlotot = 0 ; mlolotot = 0
    DO nn = 1, atoms%ntype
       mlotot = mlotot + atoms%nlo(nn)
       mlolotot = mlolotot + atoms%nlo(nn)*(atoms%nlo(nn)+1)/2
    ENDDO
    nlot_d = atoms%nlotot !max(atoms%nlotot,1)
    ALLOCATE ( kveclo(nlot_d) )
    !     ..
    nbands     = 0
    bas1 = 0 ; bas2 = 0
    l_hybrid   = (&
         xcpot%icorr == icorr_pbe0 .OR.&
         xcpot%icorr == icorr_hse  .OR.&
         xcpot%icorr == icorr_vhse .OR.&
         xcpot%icorr == icorr_hf   .OR.&
         xcpot%icorr == icorr_exx)
    l_real=sym%invs.and..not.noco%l_noco
    if (noco%l_soc.and.l_real.and.l_hybrid ) THEN
       CALL juDFT_error('hybrid functional + SOC + inv.symmetry is not tested', calledby='eigen')
    END IF

    !
    !  if gw = 1 or 2, we are in the first or second run of a GW  calculation
    !  if gw = 1 we just proceed as normal (one round),
    !  if gw = 2 it's the second run: write out the eigenfunctions and
    !  the matrix elements with the xc-potential (needs two rounds)
    !  if gw = 3 energy-independet hermitian self-energy is read in from file
    !            spex.qsgw, transformed to the APW basis, and SCF is performed
    !
    gwc = 1
    fh = 0
    !
    ! look, if WU diagonalisation
    !
    IF (it.LT.input%isec1) THEN
       IF (mpi%irank.eq.0) WRITE (6,FMT=8110) it,input%isec1
8110   FORMAT (' IT=',i4,'  ISEC1=',i4,' standard diagonalization')
       l_wu = .false.
    ELSE
       IF (mpi%irank.eq.0) WRITE (6,FMT=8120) it,input%isec1
8120   FORMAT (' IT=',i4,'  ISEC1=',i4,' reduced diagonalization')
       l_wu = .true.
    END IF
    !
    ! load potential from file pottot (=unit 8)
    !
    ALLOCATE ( vpw(stars%n3d,dimension%jspd),vzxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd) )
    ALLOCATE ( vz(vacuum%nmzd,2,4), vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd) )
    ALLOCATE ( vr0(atoms%jmtd,atoms%ntype,dimension%jspd) ) ; vr0 = 0
    OPEN (nu,file='pottot',form='unformatted',status='old')
    IF (input%gw.eq.2) THEN
       ALLOCATE ( vpwtot(stars%n3d,dimension%jspd), vrtot(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd) )
       IF ( mpi%irank == 0 ) WRITE(6,'(A/A/A/A)')&
            &  'Info: vxc matrix elements for GW will be calculated in gw_vxc',&
            &  'Info: and stored in "vxc", the values obtained from the',&
            &  'Info: original implementation are saved to "vxc.old".'
    ENDIF
999 CONTINUE
    CALL loddop(stars,vacuum,atoms,sphhar, input,sym, nu, iter,vr,vpw,vz,vzxy)
    CLOSE(nu)
    IF (mpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/),(/it,iter/),&
                                                    reshape((/19,13,5,5/),(/2,2/)))

    !
    ! some modifications for gw-calculations
    !
    IF (input%gw.eq.2.and.gwc.eq.1) THEN
       vrtot(:,:,:,:)  = vr  ! store potential for subroutine gw_vxc
       vpwtot(:,:) = vpw !
    ENDIF

    IF (gwc==1) THEN
       vr0(:,:,:) = vr(:,0,:,:)
       lh0 = 1
    ELSE IF (gwc==2) THEN
       lh0 = 0                         ! for a input%gw-calculation, we
                                       ! now evaluate matrix elements
       DO jsp = 1,input%jspins               ! with the coulomb potential
          DO nn = 1,atoms%ntype                ! but with explicit kinetic energy
             DO j = 1,atoms%jri(nn)
                vr(j,0,nn,jsp) = vr(j,0,nn,jsp)-vr0(j,nn,jsp)*sfp_const/atoms%rmsh(j,nn)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    INQUIRE(file='fleur.qsgw',EXIST=lcal_qsgw)
    lcal_qsgw = .not. lcal_qsgw

    !
    ! set energy parameters (normally to that, what we read in)
    !
    IF (gwc /= 2) THEN
       CALL lodpot(mpi,atoms,sphhar,obsolete,vacuum,&
            input, vr,vz, enpara_in, enpara)
    ENDIF
    !
   


    !---> set up and solve the eigenvalue problem
    !---> loop over energy windows

!check if z-reflection trick can be used

    l_zref=(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).and..not.noco%l_noco) 


#if ( defined(CPP_MPI))
    IF (mpi%n_size > 1) l_zref = .false.
    IF ( hybrid%l_calhf ) THEN
       call judft_error("BUG parallelization in HF case must be fixed")
       !n_start  = 1
       !n_stride = 1
    END IF
#endif
    !Count number of matrix columns on this PE
    n=0
    DO i=1+mpi%n_rank,dimension%nbasfcn,mpi%n_size
       n=n+1
    enddo
    IF (mpi%n_size>1) THEN
       matsize = dimension%nbasfcn * n
    ELSE
       matsize = (dimension%nbasfcn+1)*dimension%nbasfcn/2
    ENDIF
    ne = max(5,dimension%neigd)

    if (l_hybrid.or.hybrid%l_calhf) THEN
       eig_id_hf=eig_id
    endif
    eig_id=open_eig(&
         mpi%mpi_comm,dimension%nbasfcn,dimension%neigd,kpts%nkpt,dimension%jspd,atoms%lmaxd,&
         atoms%nlod,atoms%ntype,atoms%nlotot,noco%l_noco,.true.,l_real,noco%l_soc,.false.,mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=dimension%nstd,nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.or.input%cdinf,l_mcd=banddos%l_mcd,l_orb=banddos%l_orb)

    IF (l_real) THEN
       ALLOCATE ( hamOvlp%a_r(matsize), stat = err )
    ELSE
       ALLOCATE ( hamOvlp%a_c(matsize), stat = err )
    endif
    IF (err.NE.0) THEN
       WRITE (*,*) 'eigen: an error occured during allocation of'
       WRITE (*,*) 'the Hamilton Matrix: ',err,'  size: ',matsize
       CALL juDFT_error("eigen: Error during allocation of Hamilton" //"matrix",calledby ="eigen")
    ENDIF
    if (l_real) THEN
       ALLOCATE ( hamOvlp%b_r(matsize), stat = err )
    else
       ALLOCATE ( hamOvlp%b_c(matsize), stat = err )
    endif

    IF (err.NE.0) THEN
       WRITE (*,*) 'eigen: an error occured during allocation of'
       WRITE (*,*) 'the overlap Matrix: ',err,'  size: ',matsize
       CALL juDFT_error("eigen: Error during allocation of overlap " //"matrix",calledby ="eigen")
    ENDIF

    hamOvlp%l_real = l_real
    hamOvlp%matsize = matsize

    ALLOCATE (  matind(dimension%nbasfcn,2) )
    !
    !--->    loop over spins
    nspins = input%jspins
    IF (noco%l_noco) nspins = 1
    !
    !        Append information about file eig to gwa
    IF(input%gw.eq.2.and.gwc.eq.1) THEN
       IF ( mpi%irank == 0 ) THEN
          OPEN(15,file='gwa',status='old',form='unformatted')
          READ(15)
          READ(15)
          READ(15)
          WRITE(15) mpi%n_start,mpi%n_stride,mpi%n_rank,mpi%n_size,dimension%nvd,&
               &                 dimension%nbasfcn,atoms%nlotot
          CLOSE(15)
       END IF
    ENDIF
    !  ..
    !  LDA+U
    n_u_in=atoms%n_u
    IF ((atoms%n_u.GT.0)) THEN
       ALLOCATE( vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins) )
       CALL u_setup(sym,atoms,lmaxb,sphhar,input, enpara%el0(0:,:,:),vr,mpi, vs_mmp,results)
    ELSE
       ALLOCATE( vs_mmp(-lmaxb:-lmaxb,-lmaxb:-lmaxb,1,2) )
    ENDIF
    !
    !--->    loop over k-points: each can be a separate task

    DO jsp = 1,nspins
       !+do

       !-do

       !
       !--->       set up k-point independent t(l'm',lm) matrices
       !
       CALL timestart("tlmplm")
       err=0
       j = 1 ; IF (noco%l_noco) j = 2
       ALLOCATE(td%tuu(0:dimension%lmplmd,atoms%ntype,j),stat=err)
       ALLOCATE(td%tud(0:dimension%lmplmd,atoms%ntype,j),stat=err)
       ALLOCATE(td%tdd(0:dimension%lmplmd,atoms%ntype,j),stat=err)
       ALLOCATE(td%tdu(0:dimension%lmplmd,atoms%ntype,j),stat=err)
       mlot_d = max(mlotot,1) ; mlolot_d = max(mlolotot,1)
       ALLOCATE(td%tdulo(0:dimension%lmd,-atoms%llod:atoms%llod,mlot_d,j),stat=err)
       ALLOCATE(td%tuulo(0:dimension%lmd,-atoms%llod:atoms%llod,mlot_d,j),stat=err)
       ALLOCATE(td%tuloulo(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,mlolot_d,j), stat=err)
       ALLOCATE(td%ind(0:dimension%lmd,0:dimension%lmd,atoms%ntype,j),stat=err )
       IF (err.NE.0) THEN
          WRITE (*,*) 'eigen: an error occured during allocation of'
          WRITE (*,*) 'the tlmplm%tuu, tlmplm%tdd etc.: ',err,'  size: ',mlotot
          CALL juDFT_error("eigen: Error during allocation of tlmplm, tdd  etc.",calledby ="eigen")
       ENDIF
       CALL tlmplm(sphhar,atoms,dimension,enpara, jsp,1,mpi, vr(1,0,1,jsp),gwc,lh0,input, td,ud)
       IF (input%l_f) call write_tlmplm(td,vs_mmp,atoms%n_u>0,1,jsp,input%jspins)
       CALL timestop("tlmplm")

       !---> pk non-collinear
       !--->       call tlmplm again for the second spin direction in
       !--->       each MT, because the t-matrices are needed for both
       !--->       spins at once in hsmt
       IF (noco%l_noco) THEN
          isp = 2
          CALL timestart("tlmplm")
          CALL tlmplm(sphhar,atoms,dimension,enpara,isp,isp,mpi, vr(1,0,1,isp),gwc,lh0,input, td,ud)
          IF (input%l_f) call write_tlmplm(td,vs_mmp,atoms%n_u>0,2,2,input%jspins)
          CALL timestop("tlmplm")
       ENDIF
       !

#ifdef CPP_MPI
       ! check that all sending operations are completed
       IF ( hybrid%l_calhf ) CALL MPI_WAITALL(sndreqd,sndreq,MPI_STATUSES_IGNORE,ierr)
#endif

       k_loop:DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride
#if defined(CPP_MPI)&&defined(CPP_NEVER)
          IF ( hybrid%l_calhf ) THEN
             ! jump to next k-point if this process is not present in communicator
             IF ( comm(nk) == MPI_COMM_NULL ) CYCLE
             ! allocate buffer for communication of the results
             IF ( irank2(nk) /= 0 ) CALL work_dist_reserve_buffer( nbands(nk) )
          END IF
#endif

          nrec =  kpts%nkpt*(jsp-1) + nk
          nrec = mpi%n_size*(nrec-1) + mpi%n_rank + 1
          !
          !--->         set up lapw list
          !
          call timestart("Setup of LAPW")
          lapw%rk = 0 ; lapw%k1 = 0 ; lapw%k2 = 0 ; lapw%k3 = 0
          CALL apws(dimension,input,noco, kpts,nk,cell,l_zref, mpi%n_size,jsp, bkpt,lapw,matind,nred)

          call timestop("Setup of LAPW")
          IF (noco%l_noco) THEN
             !--->         the file potmat contains the 2x2 matrix-potential in
             !--->         the interstitial region and the vacuum
             OPEN (25,FILE='potmat',FORM='unformatted', STATUS='old')
          ENDIF
          !
          !--->         set up interstitial hamiltonian and overlap matrices
          !
          call timestart("Interstitial Hamiltonian&Overlap")
          CALL hsint(input,noco,jij,stars, vpw(:,jsp),lapw,jsp, mpi%n_size,mpi%n_rank,kpts%bk(:,nk),cell,atoms,l_real,hamOvlp)

          call timestop("Interstitial Hamiltonian&Overlap")
          !
          !--->         update with sphere terms
          !
          IF (.not.l_wu) THEN
             call timestart("MT Hamiltonian&Overlap")
             CALL hsmt(dimension,atoms,sphhar,sym,enpara, mpi%SUB_COMM,mpi%n_size,mpi%n_rank,jsp,input,mpi,&
                  lmaxb,gwc, noco,cell, lapw, bkpt,vr, vs_mmp, oneD,ud, kveclo,td,l_real,hamOvlp)
             call timestop("MT Hamiltonian&Overlap")
          ENDIF
          !
#ifdef CPP_NOTIMPLEMENTED
          IF( l_hybrid ) THEN

             CALL hsfock(nk,atoms,lcutm,obsolete,lapw, dimension,kpts,jsp,input,hybrid,maxbasm,&
                  maxindxp,maxlcutm,maxindxm,nindxm, basm,bas1,bas2,bas1_MT,drbas1_MT,ne_eig,eig_irr,&
                  mpi%n_size,sym,cell, noco,noco,oneD, nbasp,nbasm, results,results,it,nbands(nk),maxbands,nobd,&
                  mnobd,xcpot, core1,core2,nindxc,maxindxc,lmaxc, lmaxcd, kveclo_eig,maxfac,fac,sfac,gauntarr,&
                  nindxp,prod,prodm,gwc, mpi,irank2(nk),isize2(nk),comm(nk), a)

             IF ( irank2(nk) /= 0 ) CYCLE

             IF( hybrid%l_subvxc ) THEN
                CALL subvxc(lapw,kpts(:,nk),obsolete,dimension, input,jsp,atoms, hybrid,matsize,enpara%el0,enpara%ello0,&
                     sym, nlot_d,kveclo, cell,sphhar, stars,stars, xcpot,mpi, irank2(nk),vacuum,&
                     oneD, vr(:,:,:,jsp),vpw(:,jsp), a)
             END IF

          END IF ! l_hybrid
#endif
          !
          !--->         update with vacuum terms
          !
          call timestart("Vacuum Hamiltonian&Overlap")
          IF (input%film .AND. .NOT.oneD%odi%d1) THEN
             CALL hsvac(vacuum,stars,dimension, atoms, jsp,input,vzxy(1,1,1,jsp),vz,enpara%evac0,cell, &
                  bkpt,lapw,sym, noco,jij, mpi%n_size,mpi%n_rank,nv2,l_real,hamOvlp)
          ELSEIF (oneD%odi%d1) THEN
             CALL od_hsvac(vacuum,stars,dimension, oneD,atoms, jsp,input,vzxy(1,1,1,jsp),vz, &
                  enpara%evac0,cell, bkpt,lapw, oneD%odi%M,oneD%odi%mb,oneD%odi%m_cyl,oneD%odi%n2d, &
                  mpi%n_size,mpi%n_rank,sym,noco,jij,nv2,l_real,hamOvlp)
          END IF
          call timestop("Vacuum Hamiltonian&Overlap")

#ifdef CPP_NOTIMPLEMENTED
          IF ( input%gw.eq.3.or.(input%gw.eq.2.and.gwc.eq.1.and..not.lcal_qsgw)) THEN

             CALL gw_qsgw ( lcal_qsgw, b,cell,sym,atoms,&
                  jsp,dimension,lapw, nk,kpts, matsize,oneD%tau,noco, a )


          END IF

          IF (gwc==2) THEN
             CALL gw_eig(eig_id,nk,kpts,atoms,dimension,neigd,sym,&
                  kveclo,cell, ud%us(0,1,jsp),ud%dus(0,1,jsp),ud%uds(0,1,jsp),&
                  ud%duds(0,1,jsp),ud%ddn(0,1,jsp),ud%ulos(1,1,jsp),ud%uulon(1,1,jsp),ud%dulon(1,1,jsp),&
                  ud%dulos(1,1,jsp),nrec,noco,jsp,matsize,a,sphhar,stars,stars,&
                  vrtot(1,0,1,jsp),vpwtot,vr,vpw,vs_mmp(-lmaxb,-lmaxb,1,jsp),lmaxb,oneD)
             CYCLE k_loop
          ENDIF
#endif
          IF (noco%l_noco) CLOSE (25)

          !write overlap matrix b to direct access file olap
          inquire(file='olap',exist=l_file)
          if (l_file) THEN
             if (l_real) THEN
                OPEN(88,file='olap',form='unformatted',access='direct', recl=matsize*8)
                WRITE(88,rec=nrec) hamOvlp%b_r
                CLOSE(88)
             else
                OPEN(88,file='olap',form='unformatted',access='direct', recl=matsize*16)
                WRITE(88,rec=nrec) hamOvlp%b_c
                CLOSE(88)
             endif
          endif

       
          CALL eigen_diag(jsp,eig_id,it,atoms,dimension,matsize,mpi, mpi%n_rank,mpi%n_size,ne,nk,lapw,input,&
               nred,mpi%sub_comm, sym,l_zref,matind,kveclo, noco,cell,bkpt,enpara%el0,jij,l_wu,&
               oneD,td,ud, eig,ne_found,hamOvlp,zMat)
          
          !
          !--->         output results
          !
          CALL timestart("EV output")
          ne_all=ne_found
#if defined(CPP_MPI)
          !Collect number of all eigenvalues
          CALL MPI_ALLREDUCE(ne_found,ne_all,1,MPI_INTEGER,MPI_SUM, mpi%sub_comm,ierr)
#endif
          !jij%eig_l = 0.0 ! need not be used, if hdf-file is present
          if (.not.l_real) THEN
             IF (.not.jij%l_J) THEN
                zMat%z_c(:lapw%nmat,:ne_found) = conjg(zMat%z_c(:lapw%nmat,:ne_found))
             ELSE
                zMat%z_c(:lapw%nmat,:ne_found) = cmplx(0.0,0.0)
             ENDIF
          endif
          if (l_real) THEN
             CALL write_eig(eig_id, nk,jsp,ne_found,ne_all,lapw%nv(jsp),lapw%nmat,&
                  lapw%k1(:lapw%nv(jsp),jsp),lapw%k2 (:lapw%nv(jsp),jsp),lapw%k3(:lapw%nv(jsp),jsp),&
                  bkpt, kpts%wtkpt(nk),eig(:ne_found),enpara%el0(0:,:,jsp), enpara%ello0(:,:,jsp),enpara%evac0(:,jsp),&
                  atoms%nlotot,kveclo,mpi%n_size,mpi%n_rank,z=zMat%z_r(:,:ne_found))
          else
             CALL write_eig(eig_id, nk,jsp,ne_found,ne_all,lapw%nv(jsp),lapw%nmat,&
                  lapw%k1(:lapw%nv(jsp),jsp),lapw%k2 (:lapw%nv(jsp),jsp),lapw%k3(:lapw%nv(jsp),jsp),&
                  bkpt, kpts%wtkpt(nk),eig(:ne_found),enpara%el0(0:,:,jsp), enpara%ello0(:,:,jsp),enpara%evac0(:,jsp),&
                  atoms%nlotot,kveclo,mpi%n_size,mpi%n_rank,z=zMat%z_c(:,:ne_found))
          endif
          IF (noco%l_noco) THEN
             CALL write_eig(eig_id, nk,2,ne_found,ne_all,lapw%nv(2),lapw%nmat,&
                  lapw%k1(:lapw%nv(2),2),lapw%k2 (:lapw%nv(2),2),lapw%k3(:lapw%nv(2),2),&
                  bkpt, kpts%wtkpt(nk),eig(:ne_found),enpara%el0(0:,:,2), enpara%ello0(:,:,2),enpara%evac0(:,2),&
                  atoms%nlotot,kveclo)
          ENDIF

#if defined(CPP_MPI)&&defined(CPP_NEVER)
          IF ( hybrid%l_calhf ) THEN
             IF ( isize2(nk) == 1 ) THEN
                WRITE(*,'(a,i6,a,i6,a)') 'HF: kpt ', nk, ' was done by rank ', mpi%irank, '.'
             ELSE
                WRITE(*,'(a,i6,a,i6,a,i6,a)')&
                     'HF: kpt ', nk, ' was done by rank ', mpi%irank, ' and ', isize2(nk)-1, ' more.'
             END IF
             !                ELSE
             !                  WRITE(*,'(a,i6,a,i6,a)')        '    kpt ', nk, ' was done by rank ', irank, '.'
          END IF
          !#             else
          !                WRITE (*,*) 'pe: ',irank,' wrote ',nrec
#             endif
          CALL timestop("EV output")
          !#ifdef CPP_MPI
          if (l_real) THEN
             DEALLOCATE ( zMat%z_r )
          else
             DEALLOCATE ( zMat%z_c )
endif
          !
       END DO  k_loop

       DEALLOCATE (td%tuu,td%tud,td%tdu,td%tdd)
       DEALLOCATE (td%ind,td%tuulo,td%tdulo)
       DEALLOCATE (td%tuloulo)
#ifdef CPP_NEVER
       IF ( hybrid%l_calhf ) THEN
          DEALLOCATE ( eig_irr,kveclo_eig )
       END IF
#endif
    END DO ! spin loop ends
    DEALLOCATE( vs_mmp )
    DEALLOCATE (matind)
    if (l_real) THEN
       deallocate(hamOvlp%a_r,hamOvlp%b_r)
    else
       deallocate(hamOvlp%a_c,hamOvlp%b_c)
    endif
#ifdef CPP_NEVER
    IF( hybrid%l_calhf ) THEN
       DEALLOCATE( fac,sfac,gauntarr )
       DEALLOCATE( nindxc,core1,core2,nbasm,eig_c )
    END IF
#endif
#if defined(CPP_MPI)&&defined(CPP_NEVER)
    IF ( hybrid%l_calhf ) DEALLOCATE (nkpt_EIBZ)
#endif

    IF ( input%gw.eq.2.AND.(gwc==1) )  THEN        ! go for another round
       OPEN (nu,file='potcoul',form='unformatted',status='old')
       !
       !       Generate input file abcoeff for subsequent GW calculation
       !       28.10.2003 Arno Schindlmayr
       !
       IF ( mpi%irank == 0 ) THEN
          WRITE(6,'(A)') 'Info: Write out vxc for GW and vxc.old.'
          WRITE(6,'(A)') 'Info: Write out abcoeff for GW.'
          WRITE(6,'(A)') 'Info: Write out radfun for gw_vxc and GW.'
       END IF
       OPEN (12,file='vxc.old',form='formatted',status='unknown') ! contains vxc from gw_eig
       OPEN (13,file='vxc',form='formatted',status='unknown')     ! contains vxc from gw_vxc
       OPEN (1013,file='vxcfull',form='unformatted',status='unknown')
       INQUIRE(file='fleur.qsgw',exist=l_file)
       IF(l_file) THEN
          WRITE(6,'(A)') 'Info: Write file qsgw for GW.'
          OPEN(1014,file='qsgw',form='unformatted')
       ENDIF
       OPEN (15,file='abcoeff',form='unformatted',status='unknown', action='write')
       OPEN (14,file='radfun',form='unformatted',status='unknown')
       WRITE(14) atoms%jri(1:atoms%ntype)
       OPEN (16,file='latharm',form='unformatted',status='unknown')
       WRITE(16) sphhar%nlhd,sphhar%memd
       l = 0
       DO i = 1,atoms%ntype
          j = atoms%ntypsy(sum(atoms%neq(:i-1))+1)
          WRITE(16) sphhar%nlh(j),sphhar%llh(:sphhar%nlh(j),j),sphhar%nmem(:sphhar%nlh(j),j),&
               sphhar%mlh(:sphhar%memd,:sphhar%nlh(j),j),sphhar%clnu(:sphhar%memd,:sphhar%nlh(j),j)
          DO j = 1,atoms%neq(i)
             l = l + 1
             IF(atoms%invsat(l).EQ.2) THEN
                WRITE(16) -atoms%ngopr(sym%invsatnr(l))
             ELSE
                WRITE(16)  atoms%ngopr(l)
             ENDIF
          ENDDO
       ENDDO
       CLOSE (16)
       gwc=2
       GOTO 999
    ELSE IF ( input%gw.eq.2.AND.(gwc==2) )  THEN
       CLOSE (12)
       CLOSE (13)
       CLOSE (1013)
       CLOSE (14)
       CLOSE (15)
       IF(.NOT.noco%l_soc)  THEN
          INQUIRE(1014,opened=l_file)
          IF(l_file) CLOSE(1014)
          INQUIRE(667,opened=l_file)
          IF(l_file) CLOSE(667)
          CALL juDFT_end("GW finished",mpi%irank)
       ENDIF
    ENDIF

    !     hf: write out radial potential vr0
    IF (l_hybrid.or.hybrid%l_calhf) THEN
       open(unit=120,file='vr0',form='unformatted')
       DO isp=1,dimension%jspd
          DO nn=1,atoms%ntype
             DO i=1,atoms%jmtd
                WRITE(120) vr0(i,nn,isp)
             END DO
          END DO
       END DO
       CLOSE(120)
    ENDIF

    DEALLOCATE ( vpw,vzxy,vz,vr,vr0 )

#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
    if (l_hybrid.or.hybrid%l_calhf) CALL close_eig(eig_id_hf)
    atoms%n_u=n_u_in


    IF( input%jspins .EQ. 1 .AND. l_hybrid ) THEN
       results%te_hfex%valence = 2*results%te_hfex%valence
       results%te_hfex%core    = 2*results%te_hfex%core
    END IF
    enpara_in%epara_min = minval(enpara%el0)
    enpara_in%epara_min = min(minval(enpara%ello0),enpara_in%epara_min)
!    enpara_in=enpara
  END SUBROUTINE eigen
END MODULE m_eigen
