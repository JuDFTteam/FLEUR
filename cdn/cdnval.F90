!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnval

USE m_juDFT
#ifdef CPP_MPI
use mpi
#endif

CONTAINS

SUBROUTINE cdnval(eig_id, fmpi,kpts,jspin,noco,nococonv,input,banddos,cell,atoms,enpara,stars,&
                  vacuum,sphhar,sym,vTot ,cdnvalJob,den,regCharges,dos,vacdos,results,&
                  moments,gfinp,hub1inp,hub1data,coreSpecInput,mcd,slab,orbcomp,jDOS,greensfImagPart)

   !************************************************************************************
   !     This is the FLEUR valence density generator
   !******** ABBREVIATIONS *************************************************************
   !     noccbd   : number of occupied bands
   !     pallst   : if set to .true. bands above the Fermi-Energy are taken into account
   !     ener     : band energy averaged over all bands and k-points,
   !                wheighted with the l-like charge of each atom type
   !     sqal     : l-like charge of each atom type. sum over all k-points and bands
   !************************************************************************************

   USE m_types
   USE m_constants
   USE m_eig66_io
   USE m_genMTBasis
   USE m_calcDenCoeffs
   USE m_mcdinit
   USE m_sympsi
   USE m_eparas      ! energy parameters and partial charges
   USE m_qal21       ! off-diagonal part of partial charges
   USE m_abcof
   USE m_nmat        ! calculate density matrix for LDA + U
   USE m_nmat21
   USE m_vacden
   USE m_pwden
   USE m_forcea8
   USE m_force_sf ! Klueppelberg (force level 3)
   USE m_checkdopall
   USE m_greensfBZint
   USE m_greensfCalcImagPart
   USE m_local_hamiltonian
   USE m_greensfCalcScalarProducts
   USE m_cdnmt       ! calculate the density and orbital moments etc.
   USE m_orbmom      ! coeffd for orbital moments
   USE m_qmtsl       ! These subroutines divide the input%film into banddos%layers
   USE m_qintsl      ! (slabs) and intergate the DOS in these banddos%layers
   USE m_orbcomp     ! calculate orbital composition (like p_x,p_y,p_z)
   USE m_jDOS
   USE m_corespec, only : l_cs    ! calculation of core spectra (EELS)
   USE m_corespec_io, only : corespec_init
   USE m_corespec_eval, only : corespec_gaunt,corespec_rme,corespec_dos,corespec_ddscs
   USE m_xmlOutput
   USE m_types_dos
   USE m_types_mcd
   USE m_types_slab
   USE m_types_jDOS
   USE m_types_vacDOS
   USE m_types_orbcomp
#ifdef CPP_MPI
   USE m_mpi_col_den ! collect density data from parallel nodes
#endif
   USE m_dfpt_rhomt
   USE m_dfpt_rhonmt

   IMPLICIT NONE

   TYPE(t_results),       INTENT(INOUT) :: results
   TYPE(t_mpi),           INTENT(IN)    :: fmpi

   TYPE(t_enpara),        INTENT(IN)    :: enpara
   TYPE(t_banddos),       INTENT(IN)    :: banddos
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_vacuum),        INTENT(IN)    :: vacuum
   TYPE(t_noco),          INTENT(IN)    :: noco
   TYPE(t_nococonv),      INTENT(IN)    :: nococonv
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_stars),         INTENT(IN)    :: stars
   TYPE(t_cell),          INTENT(IN)    :: cell
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_gfinp),         INTENT(IN)    :: gfinp
   TYPE(t_hub1inp),       INTENT(IN)    :: hub1inp
   TYPE(t_potden),        INTENT(IN)    :: vTot
   TYPE(t_cdnvalJob),     INTENT(IN)    :: cdnvalJob
   TYPE(t_potden),        INTENT(INOUT) :: den
   TYPE(t_regionCharges), INTENT(INOUT) :: regCharges
   TYPE(t_dos),           INTENT(INOUT) :: dos
   TYPE(t_vacdos),        INTENT(INOUT) :: vacdos
   TYPE(t_moments),       INTENT(INOUT) :: moments
   TYPE(t_hub1data),       OPTIONAL, INTENT(INOUT) :: hub1data
   TYPE(t_coreSpecInput),  OPTIONAL, INTENT(IN)    :: coreSpecInput
   TYPE(t_mcd),            OPTIONAL, INTENT(INOUT) :: mcd
   TYPE(t_slab),           OPTIONAL, INTENT(INOUT) :: slab
   TYPE(t_orbcomp),        OPTIONAL, INTENT(INOUT) :: orbcomp
   TYPE(t_jDOS),           OPTIONAL, INTENT(INOUT) :: jDOS
   TYPE(t_greensfImagPart),OPTIONAL, INTENT(INOUT) :: greensfImagPart

   ! Scalar Arguments
   INTEGER,               INTENT(IN)    :: eig_id, jspin

   ! Local Scalars
   INTEGER :: ikpt,ikpt_i,jsp_start,jsp_end,ispin,jsp,max_length_k_list,nk
   INTEGER :: iErr,nbands,noccbd,iType
   INTEGER :: skip_t,skip_tt,nbasfcn
   LOGICAL :: l_real, l_corespec, l_empty

   ! Local Arrays
   REAL,ALLOCATABLE :: we(:),eig(:)
   REAL :: bkpt(3)
   INTEGER,ALLOCATABLE :: ev_list(:)
   REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:) ! radial functions

   TYPE (t_lapw)              :: lapw
   TYPE (t_orb)               :: orb
   TYPE (t_denCoeffs)         :: denCoeffs
   TYPE (t_denCoeffsOffdiag)  :: denCoeffsOffdiag
   TYPE (t_force)             :: force
   TYPE (t_eigVecCoeffs)      :: eigVecCoeffs
   TYPE (t_usdus)             :: usdus
   TYPE (t_mat)               :: zMat
   TYPE (t_gVacMap)           :: gVacMap
   TYPE (t_tlmplm)            :: tlmplm
   TYPE (t_greensfBZintCoeffs):: greensfBZintCoeffs
   TYPE(t_scalarGF), ALLOCATABLE :: scalarGF(:)

   CALL timestart("cdnval")

   call timestart("init")
   l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco).AND.atoms%n_hia==0

   ! Klueppelberg (force level 3)
   IF (input%l_f.AND.(input%f_level.GE.3)) THEN
      CALL init_sf(sym,cell,atoms)
   END IF

   IF (noco%l_mperp.OR.banddos%l_jDOS) THEN
      ! when the off-diag. part of the density matrix, i.e. m_x and
      ! m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
      ! is called only once. therefore, several spin loops have been
      ! added. if l_mperp = F, these loops run only from jspin - jspin.
      jsp_start = 1
      jsp_end   = 2
   ELSE
      jsp_start = jspin
      jsp_end   = jspin
   END IF

   !Do we need to consider the unoccupied states
   l_empty = banddos%dos.or.banddos%band
   IF(gfinp%n>0 .AND. PRESENT(greensfImagPart)) THEN
      l_empty = l_empty.OR.greensfImagPart%l_calc
   ENDIF

   ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)) ! Deallocation before mpi_col_den
   ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins))
   ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins))

   ! Initializations
   CALL usdus%init(atoms,input%jspins)
   CALL denCoeffs%init(atoms,sphhar,jsp_start,jsp_end)
   ! The last entry in denCoeffsOffdiag%init is l_fmpl. It is meant as a switch to a plot of the full magnet.
   ! density without the atomic sphere approximation for the magnet. density.
   CALL denCoeffsOffdiag%init(atoms,noco,sphhar,banddos%l_jDOS,any(noco%l_unrestrictMT).OR.noco%l_mperp)
   CALL force%init1(input,atoms)
   CALL orb%init(atoms,noco,jsp_start,jsp_end)

   !Greens function always considers the empty states
   IF(gfinp%n>0 .AND. PRESENT(greensfImagPart)) THEN
      IF(greensfImagPart%l_calc) THEN
         CALL greensfBZintCoeffs%init(gfinp,atoms,noco,SIZE(cdnvalJob%ev_list))
         CALL greensfCalcScalarProducts(gfinp,atoms,input,enpara,noco,sphhar,vTot,fmpi,hub1data=hub1data,&
                                        scalarProducts=scalarGF)
      ENDIF
   ENDIF


   IF (denCoeffsOffdiag%l_fmpl.AND.(.NOT.noco%l_mperp)) CALL juDFT_error("for fmpl set noco%l_mperp = T!" ,calledby ="cdnval")
   IF (banddos%l_mcd.AND..NOT.PRESENT(mcd)) CALL juDFT_error("mcd is missing",calledby ="cdnval")

   ! calculation of core spectra (EELS) initializations -start-
   l_coreSpec = .FALSE.
   IF (PRESENT(coreSpecInput)) THEN
      CALL corespec_init(input,atoms,coreSpecInput)
      IF(l_cs.AND.(fmpi%isize.NE.1)) CALL juDFT_error('EELS + fmpi not implemented', calledby = 'cdnval')
      IF(l_cs.AND.jspin.EQ.1) CALL corespec_gaunt()
      l_coreSpec = l_cs
   END IF
   ! calculation of core spectra (EELS) initializations -end-

   IF (fmpi%irank==0) THEN
      WRITE (oUnit,FMT=8000) jspin
      CALL openXMLElementPoly('mtCharges',(/'spin'/),(/jspin/))
   END IF
8000 FORMAT (/,/,10x,'valence density: spin=',i2)

   DO iType = 1, atoms%ntype
      DO ispin = 1, input%jspins
         CALL genMTBasis(atoms,enpara,vTot,fmpi,iType,ispin,usdus,f(:,:,0:,ispin),g(:,:,0:,ispin),flo(:,:,:,ispin),&
                         hub1data=hub1data)
      END DO
      IF (noco%l_mperp.OR.banddos%l_jDOS) CALL denCoeffsOffdiag%addRadFunScalarProducts(atoms,f,g,flo,iType)
      IF (banddos%l_mcd) CALL mcd_init(atoms,banddos,input,vTot%mt(:,0,:,:),g,f,mcd,iType,jspin)
      IF (l_coreSpec) CALL corespec_rme(atoms,input,iType,29,input%jspins,jspin,results%ef,&
                                        atoms%msh,vTot%mt(:,0,:,:),f,g)
   END DO
   DEALLOCATE (f,g,flo)

   skip_tt = dot_product(enpara%skiplo(:atoms%ntype,jspin),atoms%neq(:atoms%ntype))
   IF (noco%l_soc.OR.noco%l_noco) skip_tt = 2 * skip_tt

   jsp = MERGE(1,jspin,noco%l_noco)
   call timestop("init")

   max_length_k_list=size(cdnvalJob%k_list)
#ifdef CPP_MPI   
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,max_length_k_list,1,MPI_INTEGER,MPI_MAX,fmpi%mpi_comm,ierr)
#endif
   DO ikpt_i = 1,size(cdnvalJob%k_list)
      ikpt=cdnvalJob%k_list(ikpt_i)
      bkpt=kpts%bk(:,ikpt)

      CALL lapw%init(input,noco,nococonv, kpts,atoms,sym,ikpt,cell, fmpi)
      skip_t = skip_tt
      ev_list=cdnvaljob%compact_ev_list(ikpt_i,l_empty)
      noccbd = SIZE(ev_list)
      we  = cdnvalJob%weights(ev_list,ikpt)
      eig = results%eig(ev_list,ikpt,jsp)

      IF (cdnvalJob%l_evp) THEN
         IF (minval(ev_list) > skip_tt) skip_t = 0
         IF (maxval(ev_list) <= skip_tt) skip_t = noccbd
         IF ((minval(ev_list) <= skip_tt).AND.(maxval(ev_list) > skip_tt)) skip_t = mod(skip_tt,noccbd)
      END IF

      nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
      CALL zMat%init(l_real,nbasfcn,noccbd)
      CALL read_eig(eig_id,ikpt,jsp,list=ev_list,neig=nbands,zmat=zMat)
#ifdef CPP_MPI
      CALL MPI_BARRIER(fmpi%mpi_comm,iErr) ! Synchronizes the RMA operations
#endif

      IF (noccbd.LE.0) CYCLE ! Note: This jump has to be after the MPI_BARRIER is called

      ! valence density in the atomic spheres
      CALL eigVecCoeffs%init(input,atoms,jspin,noccbd,noco%l_mperp.OR.banddos%l_jDOS)

      DO ispin = jsp_start, jsp_end
         IF (input%l_f) CALL force%init2(noccbd,input,atoms)
         CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,noco,nococonv,ispin,&
                    eigVecCoeffs%abcof(:,0:,0,:,ispin),eigVecCoeffs%abcof(:,0:,1,:,ispin),&
                    eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat,eig,force)

         IF (atoms%n_u+atoms%n_opc.GT.0) CALL n_mat(atoms,sym,noccbd,usdus,ispin,we,eigVecCoeffs,den%mmpMat(:,:,:,ispin))
         IF (atoms%n_u.GT.0.AND.noco%l_mperp.AND.(ispin==jsp_end)) THEN
            call timestart("n_mat21")
            CALL n_mat21(atoms,sym,noccbd,we,denCoeffsOffdiag,eigVecCoeffs,den%mmpMat(:,:,:,3))
            call timestop("n_mat21")

         ENDIF

         ! perform Brillouin zone integration and summation over the
         ! bands in order to determine the energy parameters for each atom and angular momentum
         call timestart("eparas")
         CALL eparas(ispin,atoms,banddos,noccbd,ev_list,fmpi,ikpt,noccbd,we,eig,&
                     skip_t,cdnvalJob%l_evp,eigVecCoeffs,usdus,regCharges,dos,mcd)

         call timestop("eparas")
         IF (noco%l_mperp.AND.(ispin==jsp_end)) then
           call timestart("qal_21")
           CALL qal_21(atoms,banddos,input,noccbd,ev_list,nococonv,eigVecCoeffs,denCoeffsOffdiag,ikpt,dos)
           call timestop("qal_21")
         endif

         ! layer charge of each valence state in this k-point of the SBZ from the mt-sphere region of the film
         IF (PRESENT(slab).and.banddos%l_slab) CALL q_mt_sl(ispin,atoms,sym,noccbd,ev_list,ikpt,noccbd,skip_t,noccbd,eigVecCoeffs,usdus,slab)

         IF(banddos%l_orb) THEN

           IF (PRESENT(orbcomp)) CALL orb_comp(banddos,ispin,ikpt,noccbd,ev_list,atoms,noccbd,usdus,eigVecCoeffs,orbcomp)
         ENDIF
         !Decomposition into total angular momentum states
         IF(banddos%dos.AND.banddos%l_jDOS) THEN
            IF(PRESENT(jDOS).AND.ispin==jsp_end) THEN
               CALL jDOS_comp(ikpt,noccbd,ev_list,we,atoms,banddos,input,usdus,&
                              denCoeffsOffdiag,eigVecCoeffs,jDOS)
            ENDIF
         ENDIF
         CALL dfpt_rhomt(atoms,we,we,noccbd,ispin,ispin,[0.0,0.0,0.0],.FALSE.,eigVecCoeffs,eigVecCoeffs,denCoeffs)
         CALL dfpt_rhonmt(atoms,sphhar,we,we,noccbd,ispin,ispin,[0.0,0.0,0.0],.FALSE.,.TRUE.,sym,eigVecCoeffs,eigVecCoeffs,denCoeffs)
         CALL dfpt_rhomtlo(atoms,noccbd,we,we,ispin,ispin,[0.0,0.0,0.0],.FALSE.,eigVecCoeffs,eigVecCoeffs,denCoeffs)
         CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we,eigVecCoeffs,eigVecCoeffs,denCoeffs,ispin,ispin,.FALSE.,[0.0,0.0,0.0])

         IF (noco%l_soc) CALL orbmom(atoms,noccbd,we,ispin,eigVecCoeffs,orb)
         IF (input%l_f) THEN
           call local_ham(sphhar,atoms,sym,noco,nococonv,enpara,fmpi,vtot,vtot,den,input,hub1inp,hub1data,tlmplm,usdus,0.0)
           CALL force%addContribsA21A12(input,atoms,sym,cell ,enpara,&
           usdus,tlmplm,vtot,eigVecCoeffs,noccbd,ispin,eig,we,results,jsp_start,jspin,nbasfcn,zMat,lapw,sphhar,lapw%gvec(1,:,:),lapw%gvec(2,:,:),lapw%gvec(3,:,:),bkpt)
         ENDIF
         IF(l_coreSpec) CALL corespec_dos(atoms,usdus,ispin,atoms%lmaxd*(atoms%lmaxd+2),kpts%nkpt,ikpt,input%neig,&
                                          noccbd,results%ef,banddos%sig_dos,eig,we,eigVecCoeffs)
      END DO ! end loop over ispin
      IF (noco%l_mperp) then
        call timestart("denCoeffsOffdiag%calcCoefficients")
        CALL dfpt_rhomt(atoms,we,we,noccbd,2,1,[0.0,0.0,0.0],.FALSE.,eigVecCoeffs,eigVecCoeffs,denCoeffs)
        CALL dfpt_rhonmt(atoms,sphhar,we,we,noccbd,2,1,[0.0,0.0,0.0],.FALSE.,.FALSE.,sym,eigVecCoeffs,eigVecCoeffs,denCoeffs)
        CALL dfpt_rhomtlo(atoms,noccbd,we,we,2,1,[0.0,0.0,0.0],.FALSE.,eigVecCoeffs,eigVecCoeffs,denCoeffs)
        CALL dfpt_rhonmtlo(atoms,sphhar,sym,noccbd,we,we,eigVecCoeffs,eigVecCoeffs,denCoeffs,2,1,.FALSE.,[0.0,0.0,0.0])
        call timestop("denCoeffsOffdiag%calcCoefficients")
      endif

      IF(gfinp%n>0 .AND. PRESENT(greensfImagPart)) THEN
         IF(greensfImagPart%l_calc) THEN
            do ispin = MERGE(1,jsp_start,gfinp%l_mperp),MERGE(3,jsp_end,gfinp%l_mperp)
               CALL greensfBZint(ikpt,noccbd,ispin,gfinp,sym,atoms,noco,nococonv,input,kpts,&
                                 scalarGF,eigVecCoeffs,greensfBZintCoeffs)
               CALL greensfCalcImagPart_single_kpt(ikpt,ikpt_i,ev_list,ispin,gfinp,atoms,input,kpts,noco,fmpi,&
                                 results,greensfBZintCoeffs,greensfImagPart)
            enddo
         ENDIF
      ENDIF

      CALL gVacMap%init(sym,atoms,vacuum,stars,lapw,input,cell,kpts,enpara,vTot,ikpt,jspin)

      ! valence density in the interstitial and vacuum region has to be called only once (if jspin=1) in the non-collinear case
      IF (.NOT.((jspin.EQ.2).AND.noco%l_noco)) THEN
         ! valence density in the interstitial region
         CALL pwden(stars,kpts,banddos ,input,fmpi,noco,nococonv,cell,atoms,sym,ikpt,&
                    jspin,lapw,noccbd,ev_list,we,eig,den,results,force%f_b8,zMat,dos)
         ! charge of each valence state in this k-point of the SBZ in the layer interstitial region of the film
         IF (PRESENT(slab).AND.banddos%l_slab) CALL q_int_sl(jspin,ikpt,stars,atoms,sym,cell,noccbd,ev_list,lapw,slab ,zMat)
         ! valence density in the vacuum region
         IF (input%film) THEN
            CALL vacden(vacuum,stars,input,cell,atoms,noco,nococonv,banddos,&
                        we,ikpt,jspin,REAL(vTot%vac(:,1,:,:)),noccbd,ev_list,lapw,enpara%evac,den,zMat,vacdos,dos)
         END IF
      END IF
      IF (input%film) CALL regCharges%sumBandsVac(vacuum,vacdos,noccbd,ikpt,jsp_start,jsp_end,eig,we)

      IF (.FALSE..AND.(banddos%dos.OR.banddos%vacdos.OR.input%cdinf)) THEN
         ! since z is no longer an argument of cdninf sympsi has to be called here!
         CALL sympsi(lapw,jspin,sym,noccbd,cell,eig,noco,dos%jsym(:,ikpt,jspin),zMat)
      END IF
   END DO ! end of k-point loop

#ifdef CPP_MPI
   !print *,"Remaining Barriers:",size(cdnvalJob%k_list)+1,max_length_k_list
   DO nk=size(cdnvalJob%k_list)+1,max_length_k_list
      CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
   ENDDO
   DO ispin = jsp_start,jsp_end
      CALL mpi_col_den(fmpi,sphhar,atoms ,stars,vacuum,input,noco,ispin,dos,vacdos,&
                       results,denCoeffs,orb,denCoeffsOffdiag,den,regCharges,mcd,slab,orbcomp,jDOS)
   END DO
#endif

   IF(gfinp%n>0 .AND. PRESENT(greensfImagPart)) THEN
      IF(greensfImagPart%l_calc) THEN
         call timestart("Green's function: Imag Part collect")
         do ispin = MERGE(1,jsp_start,gfinp%l_mperp),MERGE(3,jsp_end,gfinp%l_mperp)
            CALL greensfImagPart%collect(ispin,fmpi%mpi_comm)
         enddo
         call timestop("Green's function: Imag Part collect")
      ENDIF
   ENDIF

   IF (fmpi%irank==0) THEN
      CALL timestart("cdnmt")
      CALL cdnmt(input%jspins,input,atoms,sym,sphhar,noco,jsp_start,jsp_end,enpara,banddos,&
                 vTot%mt(:,0,:,:),denCoeffs,usdus,orb,denCoeffsOffdiag,den%mt,hub1inp,moments,jDOS,hub1data=hub1data)
      CALL timestop("cdnmt")
      IF (l_coreSpec) CALL corespec_ddscs(jspin,input%jspins)
      DO ispin = jsp_start,jsp_end
         IF (input%cdinf) THEN
            WRITE (oUnit,FMT=8210) ispin
8210        FORMAT (/,5x,'check continuity of cdn for spin=',i2)
            CALL checkDOPAll(input,sphhar,stars,atoms,sym,vacuum ,cell,den,ispin)
         END IF
         IF (input%l_f) CALL force_a8(input,atoms,sym,sphhar,ispin,vTot%mt(:,:,:,ispin),den%mt,force,fmpi,results)
      END DO
      CALL closeXMLElement('mtCharges')
   END IF

   CALL timestop("cdnval")

END SUBROUTINE cdnval

END MODULE m_cdnval
