!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_metagga
   PUBLIC  :: calc_EnergyDen
   PRIVATE :: calc_EnergyDen_auxillary_weights

   type t_RS_potden
      REAL, ALLOCATABLE :: is(:,:), mt(:,:)
   end type t_RS_potden

CONTAINS
   SUBROUTINE calc_kinEnergyDen(EnergyDen_rs, vTot_rs, den_rs, kinEnergyDen_RS)
      USE m_juDFT_stop
      IMPLICIT NONE
      REAL, INTENT(in)                 :: den_RS(:,:), EnergyDen_RS(:,:), vTot_RS(:,:)
      REAL, INTENT(inout), allocatable :: kinEnergyDen_RS(:,:)
#ifdef CPP_LIBXC
      REAL, PARAMETER                  :: eps = 1e-15

      !implicit allocation
      kinEnergyDen_RS = EnergyDen_RS - vTot_RS * den_RS

      if(any(kinEnergyDen_RS < eps)) then
         write (6,*) "         lowest kinetic energy density cutoff = ", minval(kinEnergyDen_RS)
         kinEnergyDen_RS = max(kinEnergyDen_RS, eps)
      endif

#else
      CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
   END SUBROUTINE calc_kinEnergyDen

   SUBROUTINE calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars, &
                             vacuum, DIMENSION, sphhar, sym, vTot, oneD, results, EnergyDen)
      ! calculates the energy density
      ! EnergyDen = \sum_i n_i(r) \varepsilon_i
      ! where n_i(r) is the one-particle density
      ! and \varepsilon_i are the eigenenergies

      USE m_types_setup
      USE m_types_potden
      USE m_types_kpts
      USE m_types_mpi
      USE m_types_enpara
      USE m_types_misc
      USE m_types_regionCharges
      USE m_types_dos
      USE m_types_cdnval
      USE m_cdnval

      IMPLICIT NONE

      INTEGER,           INTENT(in)           :: eig_id
      TYPE(t_mpi),       INTENT(in)           :: mpi
      TYPE(t_kpts),      INTENT(in)           :: kpts
      TYPE(t_noco),      INTENT(in)           :: noco
      TYPE(t_input),     INTENT(in)           :: input
      TYPE(t_banddos),   INTENT(in)           :: banddos
      TYPE(t_cell),      INTENT(in)           :: cell
      TYPE(t_atoms),     INTENT(in)           :: atoms
      TYPE(t_enpara),    INTENT(in)           :: enpara
      TYPE(t_stars),     INTENT(in)           :: stars
      TYPE(t_vacuum),    INTENT(in)           :: vacuum
      TYPE(t_dimension), INTENT(in)           :: DIMENSION
      TYPE(t_sphhar),    INTENT(in)           :: sphhar
      TYPE(t_sym),       INTENT(in)           :: sym
      TYPE(t_potden),    INTENT(in)           :: vTot
      TYPE(t_oneD),      INTENT(in)           :: oneD
      TYPE(t_results),   INTENT(in)           :: results
      TYPE(t_potden),    INTENT(inout)        :: EnergyDen

      ! local
      INTEGER                         :: jspin

      TYPE(t_regionCharges)           :: regCharges
      TYPE(t_dos)                     :: dos
      TYPE(t_moments)                 :: moments
      TYPE(t_results)                 :: tmp_results
      TYPE(t_cdnvalJob)               :: cdnvalJob
      TYPE(t_potden)                  :: aux_den, real_den

      CALL regCharges%init(input, atoms)
      CALL dos%init(input,        atoms, DIMENSION, kpts, vacuum)
      CALL moments%init(input,    atoms)
      tmp_results = results

      DO jspin = 1,input%jspins
         CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin)

         ! replace brillouin weights with auxillary weights
         CALL calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, cdnvalJob%weights)

         CALL cdnval(eig_id, mpi, kpts, jspin, noco, input, banddos, cell, atoms, &
                     enpara, stars, vacuum, DIMENSION, sphhar, sym, vTot, oneD, cdnvalJob, &
                     EnergyDen, regCharges, dos, tmp_results, moments)
      ENDDO

   END SUBROUTINE calc_EnergyDen

   SUBROUTINE calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, f_ik)
      USE m_types_kpts
      USE m_eig66_io
      IMPLICIT NONE
      ! calculates new (auxillary-)weights as
      ! f_iks = w_iks * E_iks
      !, where  f_iks are the new (auxillary-)weights
      ! w_iks are the weights used in brillouin zone integration
      ! E_iks are the eigen energies

      INTEGER,      INTENT(in)        :: eig_id
      INTEGER,      INTENT(in)        :: jspin
      TYPE(t_kpts), INTENT(in)        :: kpts
      REAL,         INTENT(inout)     :: f_ik(:,:) ! f_ik(band_idx, kpt_idx)

      ! local vars
      REAL                       :: w_i(SIZE(f_ik,dim=1)), e_i(SIZE(f_ik,dim=1))
      INTEGER                    :: ikpt

      DO ikpt = 1,kpts%nkpt
         CALL read_eig(eig_id,ikpt,jspin, eig=e_i, w_iks=w_i)

         f_ik(:,ikpt) = e_i * w_i
      ENDDO
   END SUBROUTINE calc_EnergyDen_auxillary_weights

   SUBROUTINE calc_kinED_pw(eig_id, mpi,kpts,jspin,noco,input,banddos,cell,atoms,enpara,stars,&
                     vacuum,dimension,sphhar,sym,vTot,oneD,cdnvalJob,kinED,regCharges,dos,results,&
                     moments,coreSpecInput,mcd,slab,orbcomp)

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
      USE m_eig66_io
      USE m_genMTBasis
      USE m_calcDenCoeffs
      USE m_mcdinit
      USE m_sympsi
      USE m_eparas      ! energy parameters and partial charges
      USE m_qal21       ! off-diagonal part of partial charges
      USE m_abcof
      USE m_nmat        ! calculate density matrix for LDA + U
      USE m_vacden
      USE m_pwden
      USE m_forcea8
      USE m_checkdopall
      USE m_cdnmt       ! calculate the density and orbital moments etc.
      USE m_orbmom      ! coeffd for orbital moments
      USE m_qmtsl       ! These subroutines divide the input%film into vacuum%layers
      USE m_qintsl      ! (slabs) and intergate the DOS in these vacuum%layers
      USE m_orbcomp     ! calculate orbital composition (like p_x,p_y,p_z)
      USE m_abcrot2
      USE m_corespec, only : l_cs    ! calculation of core spectra (EELS)
      USE m_corespec_io, only : corespec_init
      USE m_corespec_eval, only : corespec_gaunt,corespec_rme,corespec_dos,corespec_ddscs
      USE m_xmlOutput
#ifdef CPP_MPI
      USE m_mpi_col_den ! collect density data from parallel nodes
#endif

      IMPLICIT NONE

      TYPE(t_results),       INTENT(INOUT) :: results
      TYPE(t_mpi),           INTENT(IN)    :: mpi
      TYPE(t_dimension),     INTENT(IN)    :: dimension
      TYPE(t_oneD),          INTENT(IN)    :: oneD
      TYPE(t_enpara),        INTENT(IN)    :: enpara
      TYPE(t_banddos),       INTENT(IN)    :: banddos
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_vacuum),        INTENT(IN)    :: vacuum
      TYPE(t_noco),          INTENT(IN)    :: noco
      TYPE(t_sym),           INTENT(IN)    :: sym
      TYPE(t_stars),         INTENT(IN)    :: stars
      TYPE(t_cell),          INTENT(IN)    :: cell
      TYPE(t_kpts),          INTENT(IN)    :: kpts
      TYPE(t_sphhar),        INTENT(IN)    :: sphhar
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_potden),        INTENT(IN)    :: vTot
      TYPE(t_cdnvalJob),     INTENT(IN)    :: cdnvalJob
      TYPE(t_potden),        INTENT(INOUT) :: kinED
      TYPE(t_regionCharges), INTENT(INOUT) :: regCharges
      TYPE(t_dos),           INTENT(INOUT) :: dos
      TYPE(t_moments),       INTENT(INOUT) :: moments
      TYPE(t_coreSpecInput), OPTIONAL, INTENT(IN)    :: coreSpecInput
      TYPE(t_mcd),           OPTIONAL, INTENT(INOUT) :: mcd
      TYPE(t_slab),          OPTIONAL, INTENT(INOUT) :: slab
      TYPE(t_orbcomp),       OPTIONAL, INTENT(INOUT) :: orbcomp

      ! Scalar Arguments
      INTEGER,               INTENT(IN)    :: eig_id, jspin

#ifdef CPP_MPI
      INCLUDE 'mpif.h'
#endif

      ! Local Scalars
      INTEGER :: ikpt,jsp_start,jsp_end,ispin,jsp
      INTEGER :: iErr,nbands,noccbd,iType
      INTEGER :: skip_t,skip_tt,nStart,nEnd,nbasfcn
      LOGICAL :: l_orbcomprot, l_real, l_dosNdir, l_corespec

      ! Local Arrays
      REAL,    ALLOCATABLE :: we(:)
      REAL,    ALLOCATABLE :: eig(:)
      REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:) ! radial functions

      TYPE (t_lapw)             :: lapw
      TYPE (t_orb)              :: orb
      TYPE (t_denCoeffs)        :: denCoeffs
      TYPE (t_denCoeffsOffdiag) :: denCoeffsOffdiag
      TYPE (t_force)            :: force
      TYPE (t_eigVecCoeffs)     :: eigVecCoeffs
      TYPE (t_usdus)            :: usdus
      TYPE (t_mat)              :: zMat, zPrime
      TYPE (t_potden)           :: kinED_comp(3)
      TYPE (t_gVacMap)          :: gVacMap

      CALL timestart("cdnval")

      l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco)
      l_dosNdir = banddos%dos.AND.(banddos%ndir.EQ.-3)

      IF (noco%l_mperp) THEN
         ! when the off-diag. part of the desinsity matrix, i.e. m_x and
         ! m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
         ! is called only once. therefore, several spin loops have been
         ! added. if l_mperp = F, these loops run only from jspin - jspin.
         jsp_start = 1
         jsp_end   = 2
      ELSE
         jsp_start = jspin
         jsp_end   = jspin
      END IF

      ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end)) ! Deallocation before mpi_col_den
      ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end))
      ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins))

      ! Initializations
      CALL usdus%init(atoms,input%jspins)
      CALL denCoeffs%init(atoms,sphhar,jsp_start,jsp_end)
      ! The last entry in denCoeffsOffdiag%init is l_fmpl. It is meant as a switch to a plot of the full magnet.
      ! density without the atomic sphere approximation for the magnet. density. It is not completely implemented (lo's missing).
      CALL denCoeffsOffdiag%init(atoms,noco,sphhar,noco%l_mtnocopot)
      CALL force%init1(input,atoms)
      CALL orb%init(atoms,noco,jsp_start,jsp_end)

      IF (denCoeffsOffdiag%l_fmpl.AND.(.NOT.noco%l_mperp)) CALL juDFT_error("for fmpl set noco%l_mperp = T!" ,calledby ="cdnval")
      IF (l_dosNdir.AND.oneD%odi%d1) CALL juDFT_error("layer-resolved feature does not work with 1D",calledby ="cdnval")
      IF (banddos%l_mcd.AND..NOT.PRESENT(mcd)) CALL juDFT_error("mcd is missing",calledby ="cdnval")

      ! calculation of core spectra (EELS) initializations -start-
      l_coreSpec = .FALSE.
      IF (PRESENT(coreSpecInput)) THEN
         CALL corespec_init(input,atoms,coreSpecInput)
         IF(l_cs.AND.(mpi%isize.NE.1)) CALL juDFT_error('EELS + MPI not implemented', calledby = 'cdnval')
         IF(l_cs.AND.jspin.EQ.1) CALL corespec_gaunt()
         l_coreSpec = l_cs
      END IF
      ! calculation of core spectra (EELS) initializations -end-

      IF (mpi%irank==0) THEN
         WRITE (6,FMT=8000) jspin
         WRITE (16,FMT=8000) jspin
         CALL openXMLElementPoly('mtCharges',(/'spin'/),(/jspin/))
      END IF
8000  FORMAT (/,/,10x,'valence density: spin=',i2)

      DO iType = 1, atoms%ntype
         DO ispin = jsp_start, jsp_end
            CALL genMTBasis(atoms,enpara,vTot,mpi,iType,ispin,usdus,f(:,:,0:,ispin),g(:,:,0:,ispin),flo(:,:,:,ispin))
         END DO
         IF (noco%l_mperp) CALL denCoeffsOffdiag%addRadFunScalarProducts(atoms,f,g,flo,iType)
         IF (banddos%l_mcd) CALL mcd_init(atoms,input,dimension,vTot%mt(:,0,:,:),g,f,mcd,iType,jspin)
         IF (l_coreSpec) CALL corespec_rme(atoms,input,iType,dimension%nstd,input%jspins,jspin,results%ef,&
                                           dimension%msh,vTot%mt(:,0,:,:),f,g)
      END DO
      DEALLOCATE (f,g,flo)

      skip_tt = dot_product(enpara%skiplo(:atoms%ntype,jspin),atoms%neq(:atoms%ntype))
      IF (noco%l_soc.OR.noco%l_noco) skip_tt = 2 * skip_tt
      ALLOCATE (we(MAXVAL(cdnvalJob%noccbd(:))))
      ALLOCATE (eig(MAXVAL(cdnvalJob%noccbd(:))))
      jsp = MERGE(1,jspin,noco%l_noco)

      DO ikpt = cdnvalJob%ikptStart, cdnvalJob%nkptExtended, cdnvalJob%ikptIncrement

         IF (ikpt.GT.kpts%nkpt) THEN
#ifdef CPP_MPI
            CALL MPI_BARRIER(mpi%mpi_comm,iErr) ! Synchronizes the RMA operations
#endif
            EXIT
         END IF

         CALL lapw%init(input,noco, kpts,atoms,sym,ikpt,cell,.false., mpi)
         skip_t = skip_tt
         noccbd = cdnvalJob%noccbd(ikpt)
         nStart = cdnvalJob%nStart(ikpt)
         nEnd = cdnvalJob%nEnd(ikpt)
         we(1:noccbd) = cdnvalJob%weights(1:noccbd,ikpt)
         eig(1:noccbd) = results%eig(nStart:nEnd,ikpt,jsp)

         IF (cdnvalJob%l_evp) THEN
            IF (nStart > skip_tt) skip_t = 0
            IF (nEnd <= skip_tt) skip_t = noccbd
            IF ((nStart <= skip_tt).AND.(nEnd > skip_tt)) skip_t = mod(skip_tt,noccbd)
         END IF

         nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
         CALL zMat%init(l_real,nbasfcn,noccbd)
         CALL read_eig(eig_id,ikpt,jsp,n_start=nStart,n_end=nEnd,neig=nbands,zmat=zMat)

         call set_zPrime(zMat, kpts%bk(:,ikpt), lapw, zPrime)
#ifdef CPP_MPI
         CALL MPI_BARRIER(mpi%mpi_comm,iErr) ! Synchronizes the RMA operations
#endif

         IF (noccbd.LE.0) CYCLE ! Note: This jump has to be after the MPI_BARRIER is called

         CALL gVacMap%init(dimension,sym,atoms,vacuum,stars,lapw,input,cell,kpts,enpara,vTot,ikpt,jspin)

         ! valence density in the interstitial and vacuum region has to be called only once (if jspin=1) in the non-collinear case
         IF (.NOT.((jspin.EQ.2).AND.noco%l_noco)) THEN
            ! valence density in the interstitial region
            CALL pwden(stars,kpts,banddos,oneD,input,mpi,noco,cell,atoms,sym,ikpt,&
                       jspin,lapw,noccbd,we,eig,kinED,results,force%f_b8,zMat,dos)
            ! charge of each valence state in this k-point of the SBZ in the layer interstitial region of the film
            IF (l_dosNdir.AND.PRESENT(slab)) CALL q_int_sl(jspin,ikpt,stars,atoms,sym,cell,noccbd,lapw,slab,oneD,zMat)
            ! valence density in the vacuum region
            IF (input%film) THEN
               CALL vacden(vacuum,dimension,stars,oneD, kpts,input,sym,cell,atoms,noco,banddos,&
                           gVacMap,we,ikpt,jspin,vTot%vacz(:,:,jspin),noccbd,lapw,enpara%evac,eig,kinED,zMat,dos)
            END IF
         END IF
         IF (input%film) CALL regCharges%sumBandsVac(vacuum,dos,noccbd,ikpt,jsp_start,jsp_end,eig,we)

         ! valence density in the atomic spheres
         CALL eigVecCoeffs%init(input,DIMENSION,atoms,noco,jspin,noccbd)
         DO ispin = jsp_start, jsp_end
            IF (input%l_f) CALL force%init2(noccbd,input,atoms)
            CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,noco,ispin,oneD,&
                       eigVecCoeffs%acof(:,0:,:,ispin),eigVecCoeffs%bcof(:,0:,:,ispin),&
                       eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat,eig,force)
            IF (atoms%n_u.GT.0) CALL n_mat(atoms,sym,noccbd,usdus,ispin,we,eigVecCoeffs,kinED%mmpMat(:,:,:,jspin))

            ! perform Brillouin zone integration and summation over the
            ! bands in order to determine the energy parameters for each atom and angular momentum
            CALL eparas(ispin,atoms,noccbd,mpi,ikpt,noccbd,we,eig,&
                        skip_t,cdnvalJob%l_evp,eigVecCoeffs,usdus,regCharges,dos,banddos%l_mcd,mcd)

            IF (noco%l_mperp.AND.(ispin==jsp_end)) THEN
               CALL qal_21(dimension,atoms,input,noccbd,noco,eigVecCoeffs,denCoeffsOffdiag,ikpt,dos)
            ENDIF

            ! layer charge of each valence state in this k-point of the SBZ from the mt-sphere region of the film
            IF (l_dosNdir) THEN
               IF (PRESENT(slab)) CALL q_mt_sl(ispin,atoms,noccbd,ikpt,noccbd,skip_t,noccbd,eigVecCoeffs,usdus,slab)

               INQUIRE (file='orbcomprot',exist=l_orbcomprot)
               IF (l_orbcomprot) CALL abcrot2(atoms,noccbd,eigVecCoeffs,ispin) ! rotate ab-coeffs

               IF (PRESENT(orbcomp)) CALL orb_comp(ispin,ikpt,noccbd,atoms,noccbd,usdus,eigVecCoeffs,orbcomp)
            END IF

            CALL calcDenCoeffs(atoms,sphhar,sym,we,noccbd,eigVecCoeffs,ispin,denCoeffs)

            IF (noco%l_soc) CALL orbmom(atoms,noccbd,we,ispin,eigVecCoeffs,orb)
            IF (input%l_f) CALL force%addContribsA21A12(input,atoms,dimension,sym,cell,oneD,enpara,&
                                                        usdus,eigVecCoeffs,noccbd,ispin,eig,we,results)
            IF(l_coreSpec) CALL corespec_dos(atoms,usdus,ispin,dimension%lmd,kpts%nkpt,ikpt,dimension%neigd,&
                                             noccbd,results%ef,banddos%sig_dos,eig,we,eigVecCoeffs)
         END DO ! end loop over ispin
         IF (noco%l_mperp) CALL denCoeffsOffdiag%calcCoefficients(atoms,sphhar,sym,eigVecCoeffs,we,noccbd)

         IF ((banddos%dos.OR.banddos%vacdos.OR.input%cdinf).AND.(banddos%ndir.GT.0)) THEN
            ! since z is no longer an argument of cdninf sympsi has to be called here!
            CALL sympsi(lapw,jspin,sym,dimension,nbands,cell,eig,noco,dos%ksym(:,ikpt,jspin),dos%jsym(:,ikpt,jspin),zMat)
         END IF
      END DO ! end of k-point loop

#ifdef CPP_MPI
      DO ispin = jsp_start,jsp_end
         CALL mpi_col_den(mpi,sphhar,atoms,oneD,stars,vacuum,input,noco,ispin,regCharges,dos,&
                          results,denCoeffs,orb,denCoeffsOffdiag,kinED,kinED%mmpMat(:,:,:,jspin),mcd,slab,orbcomp)
      END DO
#endif

      IF (mpi%irank==0) THEN
         CALL cdnmt(input%jspins,atoms,sphhar,noco,jsp_start,jsp_end,&
                    enpara,vTot%mt(:,0,:,:),denCoeffs,usdus,orb,denCoeffsOffdiag,moments,kinED%mt)
         IF (l_coreSpec) CALL corespec_ddscs(jspin,input%jspins)
         DO ispin = jsp_start,jsp_end
            IF (input%cdinf) THEN
               WRITE (6,FMT=8210) ispin
8210           FORMAT (/,5x,'check continuity of cdn for spin=',i2)
               CALL checkDOPAll(input,dimension,sphhar,stars,atoms,sym,vacuum,oneD,cell,kinED,ispin)
            END IF
            IF (input%l_f) CALL force_a8(input,atoms,sphhar,ispin,vTot%mt(:,:,:,ispin),kinED%mt,force,results)
         END DO
         CALL closeXMLElement('mtCharges')
      END IF

      CALL timestop("cdnval")

      write (*,*) "done calc_kinED_pw"
   END SUBROUTINE calc_kinED_pw

   subroutine set_zPrime(zMat, kpt, lapw, zPrime)
      USE m_types
      implicit none
      TYPE (t_mat), intent(in) :: zMat
      REAL, intent(in)         :: kpt(3) 
      TYPE(t_lapw), intent(in) :: lapw
      TYPE (t_mat)             :: zPrime

      call zPrime%init(zMat)

      if(zPrime%l_real) then
         write (*,*) "zMat is real"
      else
         write (*,*) "zMat is complex"
      endif


   end subroutine set_zPrime
END MODULE m_metagga
