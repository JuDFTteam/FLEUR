MODULE m_cdnval
  use m_juDFT
CONTAINS
  SUBROUTINE cdnval(eig_id, mpi,kpts,jspin,sliceplot,noco, input,banddos,cell,atoms,enpara,stars,&
       vacuum,dimension, sphhar, sym,obsolete, igq_fft,vr, vz, oneD, n_mmp,results, qpw,rhtxy,&
       rho,rht,cdom,cdomvz,cdomvxy,qa21, chmom,clmom)
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
    USE m_eig66_io,ONLY: write_dos
    USE m_radfun
    USE m_radflo
    USE m_rhomt
    USE m_rhonmt
    USE m_rhomtlo
    USE m_rhonmtlo
    USE m_mcdinit
    USE m_sphpts
    USE m_points
    USE m_sympsi
    USE m_enpara, ONLY : w_enpara,mix_enpara
    USE m_eparas      ! energy parameters and partial charges
    USE m_qal21       ! off-diagonal part of partial charges
    USE m_abcof
    USE m_topulay
    USE m_nmat        ! calculate density matrix for LDA + U
    USE m_vacden
    USE m_nstm3
    USE m_pwden
    USE m_forcea8
    USE m_forcea12
    USE m_forcea21
    USE m_checkdop    ! check continuity of density on MT radius R
    USE m_int21       ! integrate (spin) off-diagonal radial functions
    USE m_int21lo     ! -"- for u_lo
    USE m_rhomt21     ! calculate (spin) off-diagonal MT-density coeff's
    USE m_rhonmt21    ! -"-                       non-MT-density coeff's
    USE m_cdnmt       ! calculate the density and orbital moments etc.
    USE m_orbmom      ! coeffd for orbital moments
    USE m_qmtsl       ! These subroutines divide the input%film into vacuum%layers
    USE m_qintsl      ! (slabs) and intergate the DOS in these vacuum%layers
    USE m_slabdim     ! (mt + interstitial)
    USE m_slabgeom    ! (written by Yu.Koroteev, 2003/2004)
    USE m_orbcomp     ! calculate corbital composition (like p_x,p_y,p_z)
    USE m_Ekwritesl   ! and write to file.
    USE m_abcrot2
    USE m_doswrite
    USE m_cylpts
    USE m_cdnread, ONLY : cdn_read0, cdn_read
#ifdef CPP_MPI
    USE m_mpi_col_den ! collect density data from parallel nodes
#endif
    USE m_types
    USE m_xmlOutput
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_mpi),INTENT(IN)   :: mpi
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_enpara),INTENT(INOUT)   :: enpara
    TYPE(t_obsolete),INTENT(IN)   :: obsolete
    TYPE(t_banddos),INTENT(IN)   :: banddos
    TYPE(t_sliceplot),INTENT(IN)   :: sliceplot
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_kpts),INTENT(IN)   :: kpts
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: eig_id,jspin

    !     .. Array Arguments ..
    COMPLEX, INTENT(INOUT) :: qpw(stars%n3d,dimension%jspd)
    COMPLEX, INTENT(INOUT) :: rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd)
    COMPLEX, INTENT(INOUT) :: cdom(stars%n3d)
    COMPLEX, INTENT(INOUT) :: cdomvz(vacuum%nmzd,2)
    COMPLEX, INTENT(INOUT) :: cdomvxy(vacuum%nmzxyd,oneD%odi%n2d-1,2)
    COMPLEX, INTENT(INOUT) :: qa21(atoms%ntypd)
    INTEGER, INTENT (IN) :: igq_fft(0:stars%kq1d*stars%kq2d*stars%kq3d-1)
    REAL, INTENT    (IN) :: vz(vacuum%nmzd,2)
    REAL, INTENT    (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL, INTENT   (OUT) :: chmom(atoms%ntypd,dimension%jspd),clmom(3,atoms%ntypd,dimension%jspd)
    REAL, INTENT (INOUT) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL, INTENT (INOUT) :: rht(vacuum%nmzd,2,dimension%jspd)
    COMPLEX, INTENT(INOUT) :: n_mmp(-3:3,-3:3,atoms%n_u)

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    !     .. Local Scalars ..
    TYPE(t_lapw):: lapw
    INTEGER :: llpd
    REAL wk,wronk,sign,emcd_lo,emcd_up
    INTEGER i,ie,iv,ivac,j,k,l,l1,lh ,n,ilo,isp,nat,&
         nbands,noded,nodeu,noccbd,nslibd,na,&
         ikpt,npd ,jsp_start,jsp_end,ispin
    INTEGER  skip_t,skip_tt
    INTEGER n_size,i_rec,n_rank ,ncored,n_start,n_end,noccbd_l
    COMPLEX,parameter:: czero=(0.0,0.0)
    LOGICAL l_fmpl,l_mcd,l_evp,l_orbcomprot
    !     ...Local Arrays ..
    INTEGER n_bands(0:dimension%neigd),ncore(atoms%ntypd)
    REAL    cartk(3),bkpt(3),xp(3,dimension%nspd),e_mcd(atoms%ntypd,input%jspins,dimension%nstd)
    REAL    ello(atoms%nlod,atoms%ntypd,dimension%jspd),evac(2,dimension%jspd)
    REAL    epar(0:atoms%lmaxd,atoms%ntypd,dimension%jspd),evdu(2,dimension%jspd)
    REAL    eig(dimension%neigd)
    REAL    vz0(2)
    REAL    uuilon(atoms%nlod,atoms%ntypd),duilon(atoms%nlod,atoms%ntypd)
    REAL    ulouilopn(atoms%nlod,atoms%nlod,atoms%ntypd)

    INTEGER, PARAMETER :: n2max_nstm3=13

    INTEGER nsld,nsl
    !
    INTEGER, ALLOCATABLE :: nmtsl(:,:),nslat(:,:)
    REAL,    ALLOCATABLE :: zsl(:,:),volsl(:)
    REAL,    ALLOCATABLE :: volintsl(:)
    REAL,    ALLOCATABLE :: qintsl(:,:),qmtsl(:,:)
    REAL,    ALLOCATABLE :: orbcomp(:,:,:),qmtp(:,:)
    REAL,    ALLOCATABLE :: qis(:,:,:),qvac(:,:,:,:),qvlay(:,:,:,:,:)
    !-new_sl
    !-dw
    INTEGER, ALLOCATABLE :: gvac1d(:),gvac2d(:) ,kveclo(:)
    INTEGER, ALLOCATABLE :: jsym(:),ksym(:)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
    COMPLEX, ALLOCATABLE :: z(:,:)
#else
    REAL,    ALLOCATABLE :: z(:,:)
#endif
    REAL,    ALLOCATABLE :: aclo(:,:,:),acnmt(:,:,:,:,:)
    REAL,    ALLOCATABLE :: bclo(:,:,:),bcnmt(:,:,:,:,:)
    REAL,    ALLOCATABLE :: cclo(:,:,:,:),ccnmt(:,:,:,:,:),we(:)
    REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:)
    REAL,    ALLOCATABLE :: uloulopn(:,:,:,:),uloulopn21(:,:,:)
    REAL,    ALLOCATABLE :: uu(:,:,:),dd(:,:,:),du(:,:,:)
    REAL,    ALLOCATABLE :: uunmt(:,:,:,:),ddnmt(:,:,:,:)
    REAL,    ALLOCATABLE :: dunmt(:,:,:,:),udnmt(:,:,:,:),sqlo(:,:,:)
    REAL,    ALLOCATABLE :: qal(:,:,:,:),sqal(:,:,:),ener(:,:,:)
    REAL,    ALLOCATABLE :: svac(:,:),pvac(:,:),mcd(:,:,:)
    REAL,    ALLOCATABLE :: enerlo(:,:,:),qmat(:,:,:,:)
    COMPLEX, ALLOCATABLE :: acof(:,:,:,:),bcof(:,:,:,:),ccof(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: acoflo(:,:,:,:),bcoflo(:,:,:,:)
    COMPLEX, ALLOCATABLE :: cveccof(:,:,:,:,:),f_a12(:,:)
    COMPLEX, ALLOCATABLE :: e1cof(:,:,:),e2cof(:,:,:),f_a21(:,:)
    COMPLEX, ALLOCATABLE :: f_b4(:,:),f_b8(:,:)
    COMPLEX, ALLOCATABLE :: aveccof(:,:,:,:),bveccof(:,:,:,:)
    COMPLEX, ALLOCATABLE :: uloulop21(:,:,:)
    COMPLEX, ALLOCATABLE :: uunmt21(:,:,:),ddnmt21(:,:,:)
    COMPLEX, ALLOCATABLE :: dunmt21(:,:,:),udnmt21(:,:,:)
    COMPLEX, ALLOCATABLE :: qstars(:,:,:,:),m_mcd(:,:,:,:)
    TYPE (t_orb),  ALLOCATABLE :: orb(:,:,:,:)
    TYPE (t_orbl), ALLOCATABLE :: orbl(:,:,:,:)
    TYPE (t_orblo),ALLOCATABLE :: orblo(:,:,:,:,:)
    TYPE (t_mt21), ALLOCATABLE :: mt21(:,:)
    TYPE (t_lo21), ALLOCATABLE :: lo21(:,:)
    TYPE (t_usdus):: usdus
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
       ALLOCATE ( mt21(0:atoms%lmaxd,atoms%ntypd),lo21(atoms%nlod,atoms%ntypd) )  ! Deallocation at end of subroutine
       ALLOCATE ( uloulopn21(atoms%nlod,atoms%nlod,atoms%ntypd) )
       ALLOCATE ( uloulop21(atoms%nlod,atoms%nlod,atoms%ntypd) )
       ALLOCATE ( qmat(0:3,atoms%ntypd,dimension%neigd,4) )
       IF (l_fmpl) THEN
          ALLOCATE ( uunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd) )
          ALLOCATE ( ddnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd) )
          ALLOCATE ( dunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd) )
          ALLOCATE ( udnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd) )
       ELSE
          ALLOCATE ( uunmt21(1,1,1),ddnmt21(1,1,1) )
          ALLOCATE ( dunmt21(1,1,1),udnmt21(1,1,1) )
       ENDIF
    ELSE
       jsp_start = jspin
       jsp_end   = jspin
       ALLOCATE ( mt21(1,1),lo21(1,1),uunmt21(1,1,1) )
       ALLOCATE ( ddnmt21(1,1,1),dunmt21(1,1,1),udnmt21(1,1,1) )
       ALLOCATE ( uloulopn21(1,1,1),uloulop21(1,1,1),qmat(1,1,1,1) )
    ENDIF
    !
    !---> if l_mperp = F, these variables are only needed for one spin
    !---> at a time, otherwise for both spins:
    !
    ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )      ! Deallocation before mpi_col_den
    ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
    ALLOCATE (   usdus%us(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE (  usdus%uds(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE (  usdus%dus(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( usdus%duds(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( usdus%ddn(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( lapw%k1(dimension%nvd,dimension%jspd),lapw%k2(dimension%nvd,dimension%jspd),lapw%k3(dimension%nvd,dimension%jspd) )
    ALLOCATE ( jsym(dimension%neigd),ksym(dimension%neigd) )
    ALLOCATE ( gvac1d(dimension%nv2d),gvac2d(dimension%nv2d) )
    ALLOCATE (  usdus%ulos(atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( usdus%dulos(atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( usdus%uulon(atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( usdus%dulon(atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( usdus%uloulopn(atoms%nlod,atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( uu(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( dd(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( du(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( uunmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( ddnmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( dunmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( udnmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( qal(0:3,atoms%ntypd,dimension%neigd,jsp_start:jsp_end) )
    ALLOCATE ( sqal(0:3,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( ener(0:3,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE (   sqlo(atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( enerlo(atoms%nlod,atoms%ntypd,jsp_start:jsp_end) )
    ALLOCATE ( svac(2,jsp_start:jsp_end) )
    ALLOCATE ( pvac(2,jsp_start:jsp_end) )
    ALLOCATE ( qstars(vacuum%nstars,dimension%neigd,vacuum%layerd,2) )
    !
    ! --> Initializations
    !
    uu(:,:,:) = 0.0 ; dd(:,:,:) = 0.0 ; du(:,:,:) = 0.0
    IF (noco%l_mperp) THEN
       mt21(:,:)%uu = czero ; mt21(:,:)%ud = czero
       mt21(:,:)%du = czero ; mt21(:,:)%dd = czero
       lo21(:,:)%uulo = czero ; lo21(:,:)%ulou = czero
       lo21(:,:)%dulo = czero ; lo21(:,:)%ulod = czero
       uloulop21(:,:,:) = czero
    ENDIF
    uunmt(:,:,:,:) = 0.0 ; ddnmt(:,:,:,:) = 0.0
    udnmt(:,:,:,:) = 0.0 ; dunmt(:,:,:,:) = 0.0
    IF (l_fmpl) THEN
       IF (.not.noco%l_mperp)  CALL juDFT_error("for fmpl set noco%l_mperp = T!" ,calledby ="cdnval")
       uunmt21(:,:,:) = czero ; ddnmt21(:,:,:) = czero
       udnmt21(:,:,:) = czero ; dunmt21(:,:,:) = czero
    ENDIF
    svac(:,:) = 0.0 ; pvac(:,:) = 0.0
    sqal(:,:,:) = 0.0 ; ener(:,:,:) = 0.0
    !+soc
    IF (noco%l_soc) THEN
       ALLOCATE ( orb(0:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
       ALLOCATE ( orbl(atoms%nlod,-atoms%llod:atoms%llod,atoms%ntypd,jsp_start:jsp_end)     )
       ALLOCATE ( orblo(atoms%nlod,atoms%nlod,-atoms%llod:atoms%llod,atoms%ntypd,jsp_start:jsp_end))
       orb(:,:,:,:)%uu = 0.0 ; orb(:,:,:,:)%dd = 0.0
       orb(:,:,:,:)%uum = czero ; orb(:,:,:,:)%uup = czero
       orb(:,:,:,:)%ddm = czero ; orb(:,:,:,:)%ddp = czero
       orbl(:,:,:,:)%uulo = 0.0 ; orbl(:,:,:,:)%dulo = 0.0
       orbl(:,:,:,:)%uulom = czero ; orbl(:,:,:,:)%uulop = czero
       orbl(:,:,:,:)%dulom = czero ; orbl(:,:,:,:)%dulop = czero
       orblo(:,:,:,:,:)%z = 0.0
       orblo(:,:,:,:,:)%p = czero ; orblo(:,:,:,:,:)%m = czero
    ELSE
       ALLOCATE ( orb(0:0,-atoms%lmaxd:-atoms%lmaxd,1,jsp_start:jsp_end) )
       ALLOCATE ( orbl(1,-atoms%llod:-atoms%llod,1,jsp_start:jsp_end) )
       ALLOCATE ( orblo(1,1,-atoms%llod:-atoms%llod,1,jsp_start:jsp_end) )
    ENDIF
    !+for
    IF (input%l_f) THEN
       ALLOCATE ( f_a12(3,atoms%ntypd),f_a21(3,atoms%ntypd) )
       ALLOCATE ( f_b4(3,atoms%ntypd),f_b8(3,atoms%ntypd) )
       f_b4(:,:) = czero  ; f_a12(:,:) = czero
       f_b8(:,:) = czero  ; f_a21(:,:) = czero
    ELSE
       ALLOCATE ( f_b8(1,1) )
    ENDIF
    !
    INQUIRE (file='mcd_inp',exist=l_mcd)
    IF (l_mcd) THEN
       OPEN (23,file='mcd_inp',STATUS='old',FORM='formatted')
       READ (23,*) emcd_lo,emcd_up
       CLOSE (23)
       ALLOCATE ( m_mcd(dimension%nstd,(3+1)**2,3*atoms%ntypd,2) )
       ALLOCATE ( mcd(3*atoms%ntypd,dimension%nstd,dimension%neigd) )
       IF (.not.banddos%dos) WRITE (*,*) 'For mcd-spectra set banddos%dos=T!'
    ELSE
       ALLOCATE ( m_mcd(1,1,1,1),mcd(1,1,1) )
    ENDIF

    ALLOCATE ( kveclo(atoms%nlotot) )

    IF (mpi%irank==0) THEN
       WRITE (6,FMT=8000) jspin
       WRITE (16,FMT=8000) jspin
       CALL openXMLElementPoly('mtCharges',(/'spin'/),(/jspin/))
    END IF
8000 FORMAT (/,/,10x,'valence density: spin=',i2)

    CALL cdn_read0(&
         eig_id,&
         mpi%irank,mpi%isize,jspin,dimension%jspd,&
         noco%l_noco,&
         ello,evac,epar,bkpt,wk,n_bands,n_size)!keep

    !+lo
    !---> if local orbitals are used, the eigenvector has a higher
    !---> dimension then nvd
    ALLOCATE ( aclo(atoms%nlod,atoms%ntypd,jsp_start:jsp_end), &
         bclo(atoms%nlod,atoms%ntypd,jsp_start:jsp_end),&
         cclo(atoms%nlod,atoms%nlod,atoms%ntypd,jsp_start:jsp_end),&
         acnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end), &
         bcnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end), &
         ccnmt(atoms%nlod,atoms%nlod,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end) )
    aclo(:,:,:) = 0.0 ; bclo(:,:,:) = 0.0 ; ccnmt(:,:,:,:,:) = 0.0
    acnmt(:,:,:,:,:)=0.0 ; bcnmt(:,:,:,:,:)=0.0 ; cclo(:,:,:,:)=0.0

    ALLOCATE ( qis(dimension%neigd,kpts%nkptd,dimension%jspd), &
         qvac(dimension%neigd,2,kpts%nkptd,dimension%jspd), &
         qvlay(dimension%neigd,vacuum%layerd,2,kpts%nkptd,dimension%jspd) )
    qvac(:,:,:,:)=0.0 ;  qvlay(:,:,:,:,:)=0.0

    skip_tt = dot_product(enpara%skiplo(:atoms%ntype,jspin),atoms%neq(:atoms%ntype))
    IF (noco%l_soc.OR.noco%l_noco)  skip_tt = 2 * skip_tt
    !-lo
    !---> set up l-dependent m.t. wavefunctions
    na = 1
    ncored = 0

    ALLOCATE ( flo(atoms%jmtd,2,atoms%nlod,dimension%jspd) )
    DO  n = 1,atoms%ntype
       IF (input%cdinf.AND.mpi%irank==0) WRITE (6,FMT=8001) n
       DO  l = 0,atoms%lmax(n)
          DO ispin = jsp_start,jsp_end
             CALL radfun(&
                  l,n,ispin,epar(l,n,ispin),vr(1,0,n,ispin),atoms,&
                  f(1,1,l,ispin),g(1,1,l,ispin),usdus,&
                  nodeu,noded,wronk)
             IF (input%cdinf.AND.mpi%irank==0) WRITE (6,FMT=8002) l,&
                  epar(l,n,ispin),usdus%us(l,n,ispin),usdus%dus(l,n,ispin),nodeu,&
                  usdus%uds(l,n,ispin),usdus%duds(l,n,ispin),noded,usdus%ddn(l,n,ispin),&
                  wronk
          END DO
          IF (noco%l_mperp) THEN
             CALL int_21(&
                  f,g,atoms,n,l,&
                  mt21(l,n)%uun,mt21(l,n)%udn,&
                  mt21(l,n)%dun,mt21(l,n)%ddn)
          END IF
       END DO
       IF (l_mcd) THEN
          CALL mcd_init(&
               atoms,input,dimension,&
               vr(:,0,:,:),g,f,emcd_up,emcd_lo,n,jspin,&
               ncore,e_mcd,m_mcd)
          ncored = max(ncore(n),ncored)
       END IF
       !
       !--->   generate the extra wavefunctions for the local orbitals,
       !--->   if there are any.
       !
       IF ( atoms%nlo(n) > 0 ) THEN
          DO ispin = jsp_start,jsp_end
             CALL radflo(atoms,n,ispin, ello(1,1,ispin),vr(:,0,n,ispin), f(1,1,0,ispin),&
                  g(1,1,0,ispin),mpi, usdus, uuilon,duilon,ulouilopn, flo(:,:,:,ispin))
          END DO
       END IF

       DO ilo = 1, atoms%nlo(n)
          IF (noco%l_mperp) THEN
             CALL int_21lo(f,g,atoms,n, flo,ilo,&
                  lo21(ilo,n)%uulon,lo21(ilo,n)%dulon,&
                  lo21(ilo,n)%uloun,lo21(ilo,n)%ulodn,&
                  uloulopn21(1,1,n))
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

    IF (input%film) vz0(:) = vz(vacuum%nmz,:)
    nsld=1
    !+q_sl
    IF ((banddos%ndir.EQ.-3).AND.banddos%dos) THEN
       IF (oneD%odi%d1)  CALL juDFT_error("layer-resolved feature does not work with 1D",calledby ="cdnval")
       CALL slab_dim(atoms, nsld)
       ALLOCATE ( nmtsl(atoms%ntypd,nsld),nslat(atoms%natd,nsld) )
       ALLOCATE ( zsl(2,nsld),volsl(nsld) )
       ALLOCATE ( volintsl(nsld) )
       CALL slabgeom(&
            atoms,cell,nsld,&
            nsl,zsl,nmtsl,nslat,volsl,volintsl)

       ALLOCATE ( qintsl(nsld,dimension%neigd))
       ALLOCATE ( qmtsl(nsld,dimension%neigd))
       ALLOCATE ( orbcomp(dimension%neigd,23,atoms%natd) )
       ALLOCATE ( qmtp(dimension%neigd,atoms%natd) )
       IF (.NOT.input%film) qvac(:,:,:,jspin) = 0.0
    ELSE
       ALLOCATE(nmtsl(1,1),nslat(1,1),zsl(1,1),volsl(1),volintsl(1))
       ALLOCATE(qintsl(1,1),qmtsl(1,1),orbcomp(1,1,1),qmtp(1,1))
    END IF
    !-q_sl
    !
    !-->   loop over k-points: each can be a separate task
    !
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
    DO ikpt = 1,kpts%nkpt
       i_rec = i_rec + 1
       IF ((mod(i_rec-1,mpi%isize).EQ.mpi%irank).OR.l_evp) THEN
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
          IF (.NOT.ALLOCATED(z)) ALLOCATE (z(dimension%nbasfcn,dimension%neigd))
          z = 0
          CALL cdn_read(&
               eig_id,dimension%nvd,dimension%jspd,mpi%irank,mpi%isize,&
               ikpt,jspin,dimension%nbasfcn,noco%l_ss,noco%l_noco,&
               noccbd,n_start,n_end,&
               lapw%nmat,lapw%nv,ello,evdu,epar,kveclo,&
               lapw%k1,lapw%k2,lapw%k3,bkpt,wk,nbands,eig,z)
          !IF (l_evp.AND.(isize.GT.1)) THEN
          !  eig(1:noccbd) = eig(n_start:n_end)
          !ENDIF
          !
          IF (vacuum%nstm.EQ.3.AND.input%film) THEN
             CALL nstm3(&
                  sym,atoms,vacuum,stars,ikpt,lapw%nv(jspin),&
                  input,jspin,kpts,&
                  cell,wk,lapw%k1(:,jspin),lapw%k2(:,jspin),&
                  evac(1,jspin),vz,vz0,&
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
             IF (input%pallst) we(:nbands) = wk
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
                      z(:,nslibd) = z(:,i)
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
                   z(:,nslibd) = z(:,sliceplot%nnne)
                ELSE
                   DO i = 1,nbands
                      IF (eig(i).GE.sliceplot%e1s .AND. eig(i).LE.sliceplot%e2s) THEN
                         nslibd = nslibd + 1
                         eig(nslibd) = eig(i)
                         we(nslibd) = we(i)
                         z(:,nslibd) = z(:,i)
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
          we(:noccbd) = 2.*we(:noccbd)/input%jspins

          !---> pk non-collinear
          !--->    valence density in the interstitial and vacuum region
          !--->    has to be called only once (if jspin=1) in the non-collinear
          !--->    case
          !     ----> valence density in the interstitial region
          IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
             CALL timestart("cdnval: pwden")
             CALL pwden(&
                  stars,kpts,banddos,oneD,&
                  input,mpi,noco,cell,atoms,sym,ikpt,&
                  jspin,lapw,noccbd,&
                  igq_fft,we,z,&
                  eig,bkpt,&
                  qpw,cdom,qis,results%force,f_b8)
             CALL timestop("cdnval: pwden")
          END IF
          !+new
          !--->    charge of each valence state in this k-point of the SBZ
          !--->    in the layer interstitial region of the film
          !
          IF (banddos%dos.AND.(banddos%ndir.EQ.-3))  THEN
             IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
                CALL q_int_sl(&
                     jspin,stars,atoms,sym,&
                     volsl,volintsl,&
                     cell,&
                     z,noccbd,lapw,&
                     nsl,zsl,nmtsl,oneD,&
                     qintsl(:,:))
                !
             END IF
          END IF
          !-new c
          !--->    valence density in the vacuum region
          IF (input%film) THEN
             IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
                CALL timestart("cdnval: vacden")
                CALL vacden(&
                     vacuum,dimension,stars,oneD,&
                     kpts,input,&
                     cell,atoms,noco,banddos,&
                     gvac1d,gvac2d,&
                     we,ikpt,jspin,vz,vz0,&
                     noccbd,z,bkpt,lapw,&
                     evac,eig,&
                     rhtxy,rht,qvac,qvlay,&
                     qstars,cdomvz,cdomvxy)
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
          !--->    construct a(tilta) and b(tilta)
          IF (noco%l_mperp) THEN
             ALLOCATE ( acof(noccbd,0:dimension%lmd,atoms%natd,dimension%jspd),&
                                ! Deallocated before call to sympsi
                  bcof(noccbd,0:dimension%lmd,atoms%natd,dimension%jspd),                &
                  ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%natd,dimension%jspd) )
          ELSE
             ALLOCATE ( acof(noccbd,0:dimension%lmd,atoms%natd,jspin:jspin),&
                  bcof(noccbd,0:dimension%lmd,atoms%natd,jspin:jspin),&
                  ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%natd,jspin:jspin) )
          END IF

          DO ispin = jsp_start,jsp_end
             IF (input%l_f) THEN
                CALL timestart("cdnval: to_pulay")
                ALLOCATE (e1cof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd),&
                                ! Deallocated after call to force_a21
                     e2cof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd),&
                     acoflo(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%natd),&
                     bcoflo(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%natd),&
                     aveccof(3,noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd),&
                     bveccof(3,noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd),&
                     cveccof(3,-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%natd) )

                CALL to_pulay(input,atoms,noccbd,sym, lapw, noco,cell,bkpt, z,noccbd,eig,usdus,&
                     kveclo,ispin,oneD, acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                     e1cof,e2cof,aveccof,bveccof, ccof(-atoms%llod,1,1,1,ispin),acoflo,bcoflo,cveccof)
                CALL timestop("cdnval: to_pulay")

             ELSE
                CALL timestart("cdnval: abcof")
                CALL abcof(input,atoms,noccbd,sym, cell, bkpt,lapw,noccbd,z, usdus, noco,ispin,kveclo,oneD,&
                     acof(:,0:,:,ispin),bcof(:,0:,:,ispin),ccof(-atoms%llod:,:,:,:,ispin))
                CALL timestop("cdnval: abcof")

             END IF

             IF (atoms%n_u.GT.0) THEN
                CALL n_mat(atoms,sym,noccbd,usdus,ispin,we, acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                     ccof(-atoms%llod:,:,:,:,ispin), n_mmp)
             END IF
             !
             !--->       perform Brillouin zone integration and summation over the
             !--->       bands in order to determine the energy parameters for each
             !--->       atom and angular momentum
             !
             IF (.not.sliceplot%slice) THEN
                CALL eparas(ispin,atoms,noccbd,mpi,ikpt,noccbd,we,eig,ccof,&
                     skip_t,l_evp,acof(:,0:,:,ispin),bcof(:,0:,:,ispin),usdus,&
                     ncore,l_mcd,m_mcd,&
                     enerlo(1,1,ispin),sqlo(1,1,ispin),&
                     ener(0,1,ispin),sqal(0,1,ispin),&
                     qal(0:,:,:,ispin),mcd)

                IF (noco%l_mperp.AND.(ispin == jsp_end)) THEN
                   CALL qal_21(atoms, input,noccbd,we,ccof,&
                        noco,acof,bcof,mt21,lo21,uloulopn21,&
                        qal,qmat)
                END IF
             END IF
             !
             !+new
             !--->    layer charge of each valence state in this k-point of the SBZ
             !--->    from the mt-sphere region of the film
             !
             IF (banddos%dos.AND.(banddos%ndir.EQ.-3))  THEN
                CALL q_mt_sl(ispin, atoms,noccbd,nsld, ikpt,noccbd,ccof(-atoms%llod,1,1,1,ispin),&
                     skip_t,noccbd, acof(:,0:,:,ispin),bcof(:,0:,:,ispin),usdus,&
                     nmtsl,nsl, qmtsl(:,:))

                INQUIRE (file='orbcomprot',exist=l_orbcomprot)
                IF (l_orbcomprot) THEN                           ! rotate ab-coeffs
                   CALL abcrot2(atoms, noccbd,&
                        acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                        ccof(-atoms%llod:,:,:,:,ispin))
                END IF

                CALL orb_comp(ispin,noccbd,atoms,noccbd,usdus,acof(1:,0:,1:,ispin),bcof(1:,0:,1:,ispin),&
                     ccof(-atoms%llod:,1:,1:,1:,ispin), orbcomp, qmtp)
             END IF
             !-new
             !--->          set up coefficients for the spherical and
             CALL timestart("cdnval: rhomt")
             CALL rhomt(atoms,we,noccbd, acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                  uu(0:,:,ispin),dd(0:,:,ispin),du(0:,:,ispin))
             CALL timestop("cdnval: rhomt")
             !+soc
             IF (noco%l_soc) THEN
                CALL orbmom(atoms,noccbd, we,acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                     ccof(-atoms%llod:,:,:,:,ispin), orb(0:,-atoms%lmaxd:,:,ispin),orbl(:,-atoms%llod:,:,ispin),&
                     orblo(:,:,-atoms%llod:,:,ispin) )
             END IF
             !     -soc
             !--->          non-spherical m.t. density
             CALL timestart("cdnval: rhonmt")
             CALL rhonmt(atoms,sphhar, we,noccbd,sym, acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                  uunmt(0:,:,:,ispin),ddnmt(0:,:,:,ispin), udnmt(0:,:,:,ispin),dunmt(0:,:,:,ispin))
             CALL timestop("cdnval: rhonmt")

             !--->          set up coefficients of the local orbitals and the
             !--->          flapw - lo cross terms for the spherical and
             !--->          non-spherical mt density
             CALL timestart("cdnval: rho(n)mtlo")
             CALL rhomtlo(atoms,&
                  noccbd,we,acof(:,0:,:,ispin),bcof(:,0:,:,ispin),&
                  ccof(-atoms%llod:,:,:,:,ispin),&
                  aclo(:,:,ispin),bclo(:,:,ispin),cclo(:,:,:,ispin))
             !
             CALL rhonmtlo(&
                  atoms,sphhar,&
                  noccbd,we,acof(:,0:,:,ispin),&
                  bcof(:,0:,:,ispin),ccof(-atoms%llod:,:,:,:,ispin),&
                  acnmt(0:,:,:,:,ispin),bcnmt(0:,:,:,:,ispin),&
                  ccnmt(:,:,:,:,ispin))
             CALL timestop("cdnval: rho(n)mtlo")

             IF (input%l_f) THEN
                CALL timestart("cdnval: force_a12/21")
                IF (.not.input%l_useapw) THEN
                   CALL force_a12(atoms,noccbd,sym, dimension,cell,oneD,&
                        we,ispin,noccbd,usdus,acof(:,0:,:,ispin),&
                        bcof(:,0:,:,ispin),e1cof,e2cof, acoflo,bcoflo, results,f_a12)
                ENDIF
                CALL force_a21(input,atoms,dimension,noccbd,sym,&
                     oneD,cell,we,ispin,epar(0:,:,ispin),noccbd,eig,usdus,acof(:,0:,:,ispin),&
                     bcof(:,0:,:,ispin),ccof(-atoms%llod:,:,:,:,ispin), aveccof,bveccof,cveccof,&
                     results,f_a21,f_b4)

                DEALLOCATE (e1cof,e2cof,aveccof,bveccof)
                DEALLOCATE (acoflo,bcoflo,cveccof)
                CALL timestop("cdnval: force_a12/21")
             END IF
          END DO !--->    end loop over ispin

          IF (noco%l_mperp) THEN
             CALL rhomt21(atoms, we,noccbd,acof,bcof, ccof,&
                  mt21,lo21,uloulop21)
             IF (l_fmpl) THEN
                CALL rhonmt21(atoms,llpd,sphhar, we,noccbd,sym, acof,bcof,&
                     uunmt21,ddnmt21,udnmt21,dunmt21)
             END IF
          END IF

          DEALLOCATE (acof,bcof,ccof)
          !
199       CONTINUE
          IF ((banddos%dos .OR. banddos%vacdos .OR. input%cdinf)  ) THEN
             CALL timestart("cdnval: write_info")
             !
             !--->    calculate charge distribution of each state (l-character ...)
             !--->    and write the information to the files dosinp and vacdos
             !--->    for dos and bandstructure plots
             !

             !--dw    parallel writing of vacdos,dosinp....
             !        write data to direct access file first, write to formated file later by PE 0 only!
             !--dw    since z is no longer an argument of cdninf sympsi has to be called here!
             !
             cartk=matmul(bkpt,cell%bmat)
             IF (banddos%ndir.GT.0) THEN
                CALL sympsi(bkpt,lapw%nv(jspin),lapw%k1(:,jspin),lapw%k2(:,jspin),&
                     lapw%k3(:,jspin),sym,dimension,nbands,cell, z,eig,noco, ksym,jsym)
             END IF
             !
             !--dw   now write k-point data to tmp_dos
             !
             CALL write_dos(eig_id,ikpt,jspin,qal(:,:,:,jspin),qvac(:,:,ikpt,jspin),qis(:,ikpt,jspin),&
                  qvlay(:,:,:,ikpt,jspin),qstars,ksym,jsym,mcd,qintsl,&
                  qmtsl(:,:),qmtp(:,:),orbcomp)

             CALL timestop("cdnval: write_info")
             !-new_sl
          END IF

          !--->  end of loop over PE's
          DEALLOCATE (z)
       END IF ! --> end "IF ((mod(i_rec-1,mpi%isize).EQ.mpi%irank).OR.l_evp) THEN"
    END DO !---> end of k-point loop
    DEALLOCATE (we,f,g,usdus%us,usdus%dus,usdus%duds,usdus%uds,usdus%ddn)
    !+t3e
#ifdef CPP_MPI
    CALL timestart("cdnval: mpi_col_den")
    DO ispin = jsp_start,jsp_end
       CALL mpi_col_den(mpi,sphhar,atoms,oneD,stars,vacuum,&
            input,noco,l_fmpl,ispin,llpd, rhtxy(1,1,1,ispin),&
            rht(1,1,ispin),qpw(1,ispin), ener(0,1,ispin),sqal(0,1,ispin),&
            results,svac(1,ispin),pvac(1,ispin),uu(0,1,ispin),&
            dd(0,1,ispin),du(0,1,ispin),uunmt(0,1,1,ispin),ddnmt(0,1,1,ispin),&
            udnmt(0,1,1,ispin),dunmt(0,1,1,ispin),sqlo(1,1,ispin),&
            aclo(1,1,ispin),bclo(1,1,ispin),cclo(1,1,1,ispin),&
            acnmt(0,1,1,1,ispin),bcnmt(0,1,1,1,ispin),&
            ccnmt(1,1,1,1,ispin),enerlo(1,1,ispin),&
            orb(0,-atoms%lmaxd,1,ispin),orbl(1,-atoms%llod,1,ispin),&
            orblo(1,1,-atoms%llod,1,ispin),mt21,lo21,uloulop21,&
            uunmt21,ddnmt21,udnmt21,dunmt21,cdom,cdomvz,cdomvxy,n_mmp)
    END DO
    CALL timestop("cdnval: mpi_col_den")
#endif
    IF (((jspin.eq.input%jspins).OR.noco%l_mperp) .AND. (banddos%dos.or.banddos%vacdos.or.input%cdinf) ) THEN
       CALL timestart("cdnval: dos")
       IF (mpi%irank==0) THEN
          CALL doswrite(&
               eig_id,dimension,kpts,atoms,vacuum,&
               input,banddos,&
               sliceplot,noco,sym,&
               cell,&
               l_mcd,ncored,ncore,e_mcd,&
               results%ef,nsld,oneD)
          IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
             CALL Ek_write_sl(&
                  eig_id,dimension,kpts,atoms,vacuum,&
                  nsld,input,jspin,&
                  sym,cell,&
                  nsl,nslat)
          END IF
       END IF
#ifdef CPP_MPI                
       CALL MPI_BARRIER(mpi%mpi_comm,ie)
#endif
       CALL timestop("cdnval: dos")
    END IF

    IF (mpi%irank==0) THEN
       CALL cdnmt(&
            dimension%jspd,atoms,sphhar,llpd,&
            noco,l_fmpl,jsp_start,jsp_end,&
            epar,ello,vr(:,0,:,:),uu,du,dd,uunmt,udnmt,dunmt,ddnmt,&
            usdus,usdus%uloulopn,aclo,bclo,cclo,acnmt,bcnmt,ccnmt,&
            orb,orbl,orblo,mt21,lo21,uloulopn21,uloulop21,&
            uunmt21,ddnmt21,udnmt21,dunmt21,&
            chmom,clmom,&
            qa21,rho)

       DO ispin = jsp_start,jsp_end
          WRITE (6,*) 'Energy Parameters for spin:',ispin
          IF (.not.sliceplot%slice) THEN
             CALL mix_enpara(&
                  ispin,atoms,vacuum,obsolete,input,&
                  enpara,&
                  vr(:,0,:,:),vz,pvac(1,ispin),&
                  svac(1,ispin),&
                  ener(0,1,ispin),sqal(0,1,ispin),&
                  enerlo(1,1,ispin),&
                  sqlo(1,1,ispin))
             CALL w_enpara(&
                  atoms,jspin,input%film,&
                  enpara,16)
          END IF

          !--->      check continuity of charge density
          IF (input%cdinf) THEN
             CALL timestart("cdnval: cdninf-stuff")

             WRITE (6,FMT=8210) ispin
8210         FORMAT (/,5x,'check continuity of cdn for spin=',i2)
             IF (input%film .AND. .NOT.oneD%odi%d1) THEN
                !--->             vacuum boundaries
                npd = min(dimension%nspd,25)
                CALL points(xp,npd)
                DO ivac = 1,vacuum%nvac
                   sign = 3. - 2.*ivac
                   DO j = 1,npd
                      xp(3,j) = sign*cell%z1/cell%amat(3,3)
                   END DO
                   CALL checkdop(&
                        xp,npd,0,0,ivac,1,ispin,.true.,dimension,atoms,&
                        sphhar,stars,sym,&
                        vacuum,cell,oneD,&
                        qpw,rho,rhtxy,rht)
                END DO
             ELSE IF (oneD%odi%d1) THEN
                !-odim
                npd = min(dimension%nspd,25)
                CALL cylpts(xp,npd,cell%z1)
                CALL checkdop(&
                     xp,npd,0,0,ivac,1,ispin,.true.,dimension,atoms,&
                     sphhar,stars,sym,&
                     vacuum,cell,oneD,&
                     qpw,rho,rhtxy,rht)
                !+odim
             END IF
             !--->          m.t. boundaries
             nat = 1
             DO n = 1, atoms%ntype
                CALL sphpts(xp,dimension%nspd,atoms%rmt(n),atoms%pos(1,atoms%nat))
                CALL checkdop(&
                     xp,dimension%nspd,n,nat,0,-1,ispin,.true.,&
                     dimension,atoms,sphhar,stars,sym,&
                     vacuum,cell,oneD,&
                     qpw,rho,rhtxy,rht)
                nat = nat + atoms%neq(n)
             END DO
             CALL timestop("cdnval: cdninf-stuff")

          END IF
          !+for
          !--->      forces of equ. A8 of Yu et al.
          IF ((input%l_f)) THEN
             CALL timestart("cdnval: force_a8")
             CALL force_a8(input,atoms,sphhar, ispin, vr,rho,&
                  f_a12,f_a21,f_b4,f_b8,results%force)
             CALL timestop("cdnval: force_a8")
          END IF
          !-for
       END DO ! end of loop ispin = jsp_start,jsp_end
       CALL closeXMLElement('mtCharges')
    END IF ! end of (mpi%irank==0)
    !+t3e
    !Note: no deallocation anymore, we rely on Fortran08 :-)

    IF ((jsp_end.EQ.input%jspins)) THEN
       IF ((banddos%dos.OR.banddos%vacdos).AND.(banddos%ndir/=-2))  CALL juDFT_end("DOS OK",mpi%irank)
       IF (vacuum%nstm.EQ.3)  CALL juDFT_end("VACWAVE OK",mpi%irank)
    END IF
  END SUBROUTINE cdnval
END MODULE m_cdnval
