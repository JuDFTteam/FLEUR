MODULE m_calc_hybrid
  USE m_judft
CONTAINS
  SUBROUTINE calc_hybrid(hybrid,kpts,atoms,input,DIMENSION,mpi,noco,cell,vacuum,oneD,banddos,results,sym,xcpot,v,it  )
    USE m_types
    USE m_mixedbasis
    USE m_coulombmatrix
    USE m_hf_init
    USE m_hf_setup
    USE m_hsfock
    USE m_eig66_io
    USE m_apws
    USE m_io_hybrid
    USE m_icorrkeys
    IMPLICIT NONE
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_hybrid),INTENT(INOUT) :: hybrid
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_banddos),INTENT(IN)   :: banddos
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_sym),INTENT(IN)       :: sym  
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_atoms),INTENT(IN)     :: atoms!in u_setup n_u might be modified
    TYPE(t_potden),INTENT(IN)    :: v

    INTEGER,INTENT(IN) :: it

    !Local variables
    INTEGER        :: eig_id,jsp,nk,nred
    TYPE(t_hybdat) :: hybdat
    type(t_lapw)   :: lapw
    LOGICAL        :: init_vex=.TRUE. !In first call we have to init v_nonlocal
    INTEGER,SAVE   :: nohf_it
    INTEGER        ::  comm(kpts%nkpt),irank2(kpts%nkpt),isize2(kpts%nkpt)
    LOGICAL        :: l_restart=.FALSE.,l_zref
    INTEGER, ALLOCATABLE :: matind(:,:)
    REAL,    ALLOCATABLE    ::  eig_irr(:,:)
    real               :: bkpt(3)
    
    CALL open_hybrid_io1(DIMENSION,sym%invs)
    IF (kpts%nkptf==0) CALL judft_error("kpoint-set of full BZ not available",hint="to generate kpts in the full BZ you should specify a k-mesh in inp.xml")
    
    !Check if new non-local potential shall be generated

    hybrid%l_subvxc = ( hybrid%l_hybrid.AND.xcpot%icorr /= icorr_exx )
    IF (.NOT.ALLOCATED(v%pw)) THEN
       hybrid%l_calhf=.FALSE.
       !INQUIRE(file="v_x.mat",exist=hybrid%l_addhf)
       hybrid%l_addhf=.false.
       hybrid%l_subvxc = ( hybrid%l_subvxc .AND. hybrid%l_addhf)
       RETURN
    ENDIF
 
    hybrid%l_addhf=.true.
    !In first iteration allocate some memory
    IF (init_vex) THEN
       ALLOCATE(hybrid%ne_eig(kpts%nkpt),hybrid%nbands(kpts%nkpt),hybrid%nobd(kpts%nkptf))
       ALLOCATE( hybrid%nbasm(kpts%nkptf))
       init_vex=.false.
    ENDIF
    hybrid%l_calhf  = .TRUE.
    
    
    IF( .NOT. ALLOCATED(results%w_iks) )&
         ALLOCATE ( results%w_iks(DIMENSION%neigd2,kpts%nkpt,DIMENSION%jspd) )
    
    IF (.NOT.hybrid%l_calhf) RETURN !use existing non-local potential

    !check if z-reflection trick can be used

    l_zref=(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco) 
    
    ALLOCATE (  matind(DIMENSION%nbasfcn,2) )

    
    CALL timestart("generation of non-local HF potential")
    CALL timestart("Preparation for Hybrid functionals")
    CALL juDFT_WARN ("Hybrid functionals not working in this version")      
      
    eig_id=open_eig(&
         mpi%mpi_comm,dimension%nbasfcn,dimension%neigd,kpts%nkpt,dimension%jspd,atoms%lmaxd,atoms%nlod,atoms%ntype,atoms%nlotot&
         ,noco%l_noco,.FALSE.,sym%invs.AND..NOT.noco%l_noco,noco%l_soc,.FALSE.,&
         mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
         nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
         l_orb=banddos%l_orb)
    
    !construct the mixed-basis
    CALL timestart("generation of mixed basis")
    CALL mixedbasis(atoms,kpts, dimension,input,cell,sym,xcpot,hybrid, eig_id,mpi,v,l_restart)
    CALL timestop("generation of mixed basis")

         
    CALL open_hybrid_io2(hybrid,DIMENSION,atoms,sym%invs)

    CALL timestart("generation of coulomb matrix")
    CALL coulombmatrix(mpi,atoms,kpts,cell,sym,hybrid,xcpot,l_restart)
    CALL timestop("generation of coulomb matrix")
    
    CALL hf_init(hybrid,kpts,atoms,input,DIMENSION,hybdat,irank2,isize2,sym%invs)
    CALL timestop("Preparation for Hybrid functionals")
    
    CALL timestart("Calculation of non-local potential")
    DO jsp = 1,input%jspins
       CALL HF_setup(hybrid,input,sym,kpts,dimension,atoms,mpi,noco,cell,oneD,results,jsp,eig_id,&
            hybdat,irank2,it,sym%invs,v%mt(:,0,:,:),eig_irr)  
       DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride

          CALL apws(DIMENSION,input,noco, kpts,nk,cell,l_zref, mpi%n_size,jsp, bkpt,lapw,matind,nred)
  
          CALL hsfock(nk,atoms,hybrid,lapw,DIMENSION,kpts,jsp,input,hybdat,eig_irr,&
               sym,cell,noco,results,it,MAXVAL(hybrid%nobd),xcpot,&
               mpi,irank2(nk),isize2(nk),comm(nk))
       END DO
    END DO
    CALL timestop("Calculation of non-local potential")
    CALL timestop("generation of non-local HF potential")
    CALL close_eig(eig_id)
  end subroutine calc_hybrid
END MODULE m_calc_hybrid
