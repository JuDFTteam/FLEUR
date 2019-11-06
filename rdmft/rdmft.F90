!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rdmft

CONTAINS

SUBROUTINE rdmft(eig_id,mpi,input,kpts,banddos,sliceplot,cell,atoms,enpara,stars,vacuum,dimension,&
                 sphhar,sym,field,vTot,vCoul,oneD,noco,xcpot,mpbasis,hybrid,results,coreSpecInput,archiveType,outDen)

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_eig66_io
#ifndef CPP_OLDINTEL
   USE m_cdnval
   USE m_cdngen
   USE m_cdn_io
   USE m_cdncore
   USE m_qfix
   USE m_vgen_coulomb
   USE m_convol
   USE m_intnv

   USE m_mixedbasis
   USE m_coulombmatrix
   USE m_hf_init
   USE m_hf_setup
   USE m_io_hybrid
   USE m_symm_hf
   USE m_exchange_valence_hf
   USE m_exchange_core
   USE m_symmetrizeh
   USE m_bfgs_b2

#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif
#endif

   IMPLICIT NONE

   TYPE(t_mpi),           INTENT(IN)    :: mpi
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   TYPE(t_banddos),       INTENT(IN)    :: banddos
   TYPE(t_sliceplot),     INTENT(IN)    :: sliceplot
   TYPE(t_cell),          INTENT(IN)    :: cell
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_enpara),        INTENT(INOUT) :: enpara
   TYPE(t_stars),         INTENT(IN)    :: stars
   TYPE(t_vacuum),        INTENT(IN)    :: vacuum
   TYPE(t_dimension),     INTENT(IN)    :: dimension
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_field),         INTENT(INOUT) :: field
   TYPE(t_potden),        INTENT(INOUT) :: vTot
   TYPE(t_potden),        INTENT(INOUT) :: vCoul
   TYPE(t_oneD),          INTENT(IN)    :: oneD
   TYPE(t_noco),          INTENT(INOUT) :: noco
   TYPE(t_xcpot_inbuild), INTENT(INOUT) :: xcpot
   TYPE(t_mpbasis),       intent(inout) :: mpbasis
   TYPE(t_hybrid),        INTENT(INOUT) :: hybrid
   TYPE(t_results),       INTENT(INOUT) :: results
   TYPE(t_coreSpecInput), INTENT(IN)    :: coreSpecInput
   TYPE(t_potden),        INTENT(INOUT) :: outDen

   INTEGER,               INTENT(IN)    :: eig_id
   INTEGER,               INTENT(IN)    :: archiveType
   TYPE(t_potden)                       :: EnergyDen

#ifndef CPP_OLDINTEL
   TYPE(t_cdnvalJob)                    :: cdnvalJob
   TYPE(t_potden)                       :: singleStateDen, overallDen, overallVCoul, vTotTemp
   TYPE(t_regionCharges)                :: regCharges
   TYPE(t_dos)                          :: dos
   TYPE(t_moments)                      :: moments
   TYPE(t_mat)                          :: exMat, zMat, olap, trafo, invtrafo, tmpMat, exMatLAPW
   TYPE(t_lapw)                         :: lapw
   TYPE(t_hybdat)                       :: hybdat
   INTEGER                              :: ikpt,ikpt_i,iband_i, iBand, jkpt, jBand, iAtom, i, na, itype, lh, j
   INTEGER                              :: jspin, jspmax, jsp, isp, ispin, nbasfcn, nbands
   INTEGER                              :: nsymop, nkpt_EIBZ, ikptf, iterHF, mnobd
   INTEGER                              :: iState, iStep, numStates, maxHistoryLength, numRelevantStates
   REAL                                 :: fix, potDenInt, fermiEnergyTemp, spinDegenFac
   REAL                                 :: rdmftFunctionalValue, occStateI, gradSum
   REAL                                 :: exchangeTerm, lagrangeMultiplier, equalityCriterion
   REAL                                 :: mixParam, rdmftEnergy
   REAL                                 :: sumOcc, tempOcc, addCharge, subCharge, addChargeWeight, subChargeWeight
   REAL, PARAMETER                      :: degenEps = 0.00001
   REAL, PARAMETER                      :: convCrit = 1.0e-6
   REAL, PARAMETER                      :: minOcc = 1.0e-8
   LOGICAL                              :: converged, l_qfix, l_zref
   CHARACTER(LEN=20)                    :: filename

   INTEGER                              :: nsest(dimension%neigd) ! probably too large
   INTEGER                              :: indx_sest(dimension%neigd,dimension%neigd) ! probably too large
   INTEGER                              :: rrot(3,3,sym%nsym)
   INTEGER                              :: psym(sym%nsym) ! Note: psym is only filled up to index nsymop
   INTEGER                              :: lowestState(kpts%nkpt,input%jspins)
   INTEGER                              :: highestState(kpts%nkpt,input%jspins)
   INTEGER                              :: neigTemp(kpts%nkpt,input%jspins)

   REAL                                 :: wl_iks(dimension%neigd,kpts%nkptf)

   REAL, ALLOCATABLE                    :: overallVCoulSSDen(:,:,:)
   REAL, ALLOCATABLE                    :: vTotSSDen(:,:,:)
   REAL, ALLOCATABLE                    :: dEdOcc(:,:,:)

   REAL, ALLOCATABLE                    :: exDiag(:,:,:)
   REAL, ALLOCATABLE                    :: eig_irr(:,:)

   REAL, ALLOCATABLE                    :: gradient(:), lastGradient(:)
   REAL, ALLOCATABLE                    :: parameters(:), lastParameters(:)
   REAL, ALLOCATABLE                    :: minConstraints(:)
   REAL, ALLOCATABLE                    :: maxConstraints(:)
   REAL, ALLOCATABLE                    :: gradientCorrections(:,:)
   REAL, ALLOCATABLE                    :: paramCorrections(:,:)
   REAL, ALLOCATABLE                    :: equalityLinCombi(:)

   REAL, ALLOCATABLE                    :: occupationVec(:)

   INTEGER, ALLOCATABLE                 :: parent(:)
   INTEGER, ALLOCATABLE                 :: pointer_EIBZ(:)
   INTEGER, ALLOCATABLE                 :: n_q(:)

   LOGICAL, ALLOCATABLE                 :: enabledConstraints(:)

#endif

#ifndef CPP_OLDINTEL

   WRITE(*,*) 'entered RDMFT subroutine'

   ! General initializations
   mixParam = 0.0001
   lagrangeMultiplier = 0.1 !results%ef
   spinDegenFac = 2.0 / input%jspins ! This factor is used to compensate the missing second spin in non-spinpolarized calculations

   neigTemp(:,:) = results%neig(:,:)

   ! Determine which states have to be considered
   lowestState(:,:) = -1
   highestState(:,:) = -1
   DO jspin = 1, input%jspins
      jsp = MERGE(1,jspin,noco%l_noco)
      DO ikpt = 1, kpts%nkpt
         ! determine lowest state
         DO iBand = 1, results%neig(ikpt,jsp)
            IF((results%w_iks(iBand,ikpt,jspin) / kpts%wtkpt(ikpt)).LT.(1.0-input%rdmftOccEps)) THEN
               lowestState(ikpt,jspin) = iBand
               EXIT
            END IF
         END DO
         lowestState(ikpt,jspin) = lowestState(ikpt,jspin) - input%rdmftStatesBelow
         DO iBand = lowestState(ikpt,jspin)-1, 1, -1
            IF((results%eig(iBand+1,ikpt,jsp) - results%eig(iBand,ikpt,jsp)).GT.degenEps) THEN
               EXIT
            END IF
            lowestState(ikpt,jspin) = iBand
         END DO
         lowestState(ikpt,jspin) = MAX(lowestState(ikpt,jspin),1)

         ! determine highest state
         DO iBand = results%neig(ikpt,jsp)-1, 1, -1
            IF((results%w_iks(iBand,ikpt,jspin) / kpts%wtkpt(ikpt)).GT.(0.0+input%rdmftOccEps)) THEN
               highestState(ikpt,jspin) = iBand
               EXIT
            END IF
         END DO
         highestState(ikpt,jspin) = highestState(ikpt,jspin) + input%rdmftStatesAbove
         IF((results%neig(ikpt,jsp)-1).LT.highestState(ikpt,jspin)) THEN
            WRITE(6,*) 'Error: Not enough states calculated:'
            WRITE(6,*) 'ikpt, jsp: ', ikpt, jsp
            WRITE(6,*) 'highestState(ikpt,jspin): ', highestState(ikpt,jspin)
            WRITE(6,*) 'results%neig(ikpt,jsp): ', results%neig(ikpt,jsp)
            CALL juDFT_error('Not enough states calculated', calledby = 'rdmft')
         END IF
         DO iBand = highestState(ikpt,jspin)+1, results%neig(ikpt,jsp)
            IF((results%eig(iBand,ikpt,jsp) - results%eig(iBand-1,ikpt,jsp)).GT.degenEps) THEN
               EXIT
            END IF
            highestState(ikpt,jspin) = iBand
         END DO
         IF(highestState(ikpt,jspin).EQ.results%neig(ikpt,jsp)) THEN
            WRITE(6,*) 'Error: Highest state is degenerate:'
            WRITE(6,*) 'ikpt, jsp: ', ikpt, jsp
            WRITE(6,*) 'highestState(ikpt,jspin): ', highestState(ikpt,jspin)
            WRITE(6,*) 'results%eig(highestState(ikpt,jspin),ikpt,jsp): ', results%eig(highestState(ikpt,jspin),ikpt,jsp)
            WRITE(6,*) 'results%eig(highestState(ikpt,jspin)-1,ikpt,jsp): ', results%eig(highestState(ikpt,jspin)-1,ikpt,jsp)
            CALL juDFT_error('Highest state is degenerate', calledby = 'rdmft')
         END IF
      END DO
   END DO

   ! Move occupations of relevant states well into allowed region
   numRelevantStates = SUM(highestState(:,:)) - SUM(lowestState(:,:)) + input%jspins*kpts%nkpt
   ALLOCATE(occupationVec(numRelevantStates))
   occupationVec(:) = 0.0
   sumOcc = 0.0
   tempOcc = 0.0
   addCharge = 0.0
   subCharge = 0.0
   addChargeWeight = 0.0
   subChargeWeight = 0.0
   iState = 0
   DO jspin = 1, input%jspins
      jsp = MERGE(1,jspin,noco%l_noco)
      DO ikpt = 1, kpts%nkpt
         DO iBand = lowestState(ikpt,jspin), highestState(ikpt,jspin)
            iState = iState + 1
            occupationVec(iState) = results%w_iks(iBand,ikpt,jsp) / (kpts%wtkpt(ikpt))
            sumOcc = sumOcc + results%w_iks(iBand,ikpt,jsp)
            IF(occupationVec(iState).LT.0.01) THEN
               addCharge = addCharge + (0.01-occupationVec(iState))*kpts%wtkpt(ikpt)
               addChargeWeight = addChargeWeight + kpts%wtkpt(ikpt)
            END IF
            IF(occupationVec(iState).GT.0.99) THEN
               subCharge = subCharge + (occupationVec(iState)-0.99)*kpts%wtkpt(ikpt)
               subChargeWeight = subChargeWeight + kpts%wtkpt(ikpt)
            END IF
         END DO
      END DO
   END DO
   iState = 0
   DO jspin = 1, input%jspins
      jsp = MERGE(1,jspin,noco%l_noco)
      DO ikpt = 1, kpts%nkpt
         DO iBand = lowestState(ikpt,jspin), highestState(ikpt,jspin)
            iState = iState + 1
            IF(occupationVec(iState).LT.0.01) THEN
               occupationVec(iState) = occupationVec(iState) + 0.5*(subCharge+addCharge)*(kpts%wtkpt(ikpt)/addChargeWeight)
            END IF
            IF(occupationVec(iState).GT.0.99) THEN
               occupationVec(iState) = occupationVec(iState) - 0.5*(subCharge+addCharge)*(kpts%wtkpt(ikpt)/subChargeWeight)
            END IF
            results%w_iks(iBand,ikpt,jsp) = occupationVec(iState) * kpts%wtkpt(ikpt)
            tempOcc = tempOcc + occupationVec(iState) * kpts%wtkpt(ikpt)
         END DO
      END DO
   END DO
   DEALLOCATE(occupationVec)

   ! Some more initializations

   results%neig(:,:) = highestState(:,:) + 1

   ALLOCATE(overallVCoulSSDen(MAXVAL(results%neig(1:kpts%nkpt,1:input%jspins)),kpts%nkpt,input%jspins))
   ALLOCATE(vTotSSDen(MAXVAL(results%neig(1:kpts%nkpt,1:input%jspins)),kpts%nkpt,input%jspins))
   ALLOCATE(dEdOcc(MAXVAL(results%neig(1:kpts%nkpt,1:input%jspins)),kpts%nkpt,input%jspins))

   CALL regCharges%init(input,atoms)
   CALL dos%init(input,atoms,dimension,kpts,vacuum)
   CALL moments%init(input,atoms)
   CALL overallDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
   CALL overallVCoul%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTCOUL)
   IF (ALLOCATED(vTot%pw_w)) DEALLOCATE (vTot%pw_w)
   ALLOCATE(vTot%pw_w(SIZE(overallDen%pw,1),input%jspins))
   DO jspin = 1, input%jspins
      CALL convol(stars,vTot%pw_w(:,jspin),vTot%pw(:,jspin),stars%ufft)
   END DO

   CALL vTotTemp%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTTOT)
   vTotTemp = vTot

   DO jsp = 1,SIZE(vTot%mt,4)
      DO iAtom = 1,atoms%ntype
         vTotTemp%mt(:atoms%jri(iAtom),0,iAtom,jsp)  = sfp_const * vTot%mt(:atoms%jri(iAtom),0,iAtom,jsp) / atoms%rmsh(:atoms%jri(iAtom),iAtom)
      END DO
   END DO

   vCoul%pw_w = CMPLX(0.0,0.0)
   DO jspin = 1, input%jspins
      CALL convol(stars,vCoul%pw_w(:,jspin),vCoul%pw(:,jspin),stars%ufft)
   END DO

   vTotSSDen = 0.0

   ! Calculate all single state densities
   CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin)

   numStates = 0
   DO jspin = 1, input%jspins
      jsp = MERGE(1,jspin,noco%l_noco)
      DO ikpt_i = 1, SIZE(mpi%k_list)
         ikpt= mpi%k_list(ikpt_i)
         DO iBand_i = 1,size(cdnvalJOB%ev_list)
            iband=mpi%ev_list(iband_i)
            IF (iband>highestState(ikpt,jsp)) CYCLE
            numStates = numStates + 1
            ! Construct cdnvalJob object for this state
            ! (Reasonable parallelization is not yet done - should be placed over the loops enclosing this section)
            cdnvalJob%k_list=[ikpt]
            cdnvalJob%ev_list=[iband]
            cdnvalJob%weights(iBand,ikpt) = spinDegenFac

            ! Call cdnval to construct density
            WRITE(*,*) 'Note: some optional flags may have to be reset in rdmft before the cdnval call'
            WRITE(*,*) 'This is not yet implemented!'
            CALL singleStateDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
            CALL cdnval(eig_id,mpi,kpts,jsp,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                        sphhar,sym,vTot,oneD,cdnvalJob,singleStateDen,regCharges,dos,results,moments)

            ! Store the density on disc (These are probably way too many densities to keep them in memory)
            filename = ''
            WRITE(filename,'(a,i1.1,a,i4.4,a,i5.5)') 'cdn-', jsp, '-', ikpt, '-', iBand
            IF (mpi%irank.EQ.0) THEN
               CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                                 0,-1.0,0.0,.FALSE.,singleStateDen,TRIM(ADJUSTL(filename)))
            END IF
#ifdef CPP_MPI
            CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,singleStateDen)
#endif
            ! For each state calculate Integral over KS effective potential times single state density
            potDenInt = 0.0
            CALL int_nv(jsp,stars,vacuum,atoms,sphhar,cell,sym,input,oneD,vTotTemp,singleStateDen,potDenInt)
            vTotSSDen(iBand,ikpt,jsp) = potDenInt
         END DO
      END DO
   END DO

   ! Initializations for exchange contributions

   WRITE(*,*) 'RDMFT: HF initializations start'

   IF(ALLOCATED(hybrid%ne_eig)) DEALLOCATE(hybrid%ne_eig)
   IF(ALLOCATED(hybrid%nbands)) DEALLOCATE(hybrid%nbands)
   IF(ALLOCATED(hybrid%nobd)) DEALLOCATE(hybrid%nobd)
   IF(ALLOCATED(hybrid%nbasm)) DEALLOCATE(hybrid%nbasm)
   IF(ALLOCATED(hybrid%div_vv)) DEALLOCATE(hybrid%div_vv)
   ALLOCATE(hybrid%ne_eig(kpts%nkpt),hybrid%nbands(kpts%nkpt),hybrid%nobd(kpts%nkptf,input%jspins))
   ALLOCATE(hybrid%nbasm(kpts%nkptf))
   ALLOCATE(hybrid%div_vv(DIMENSION%neigd,kpts%nkpt,input%jspins))
   l_zref = (sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco)
   iterHF = 0
   hybrid%l_calhf = .TRUE.

!   CALL open_hybrid_io1(DIMENSION,sym%invs)

   CALL mixedbasis(atoms,kpts,input,cell,xcpot,mpbasis,hybrid,enpara,mpi,vTot)

   CALL open_hybrid_io2(mpbasis, hybrid,DIMENSION,atoms,sym%invs)

   CALL coulombmatrix(mpi,atoms,kpts,cell,sym,mpbasis,hybrid,xcpot)

   CALL hf_init(hybrid,atoms,input,DIMENSION,hybdat)

   WRITE(*,*) 'RDMFT: HF initializations end'

   ALLOCATE(parent(kpts%nkptf))
   ALLOCATE(exDiag(dimension%neigd,ikpt,input%jspins))
   ALLOCATE(lastGradient(numStates+1))
   ALLOCATE(lastParameters(numStates+1))
   lastGradient = 0.0
   lastParameters = 0.0
   maxHistoryLength = 17!7
   ALLOCATE(gradientCorrections(numStates+1,maxHistoryLength))
   ALLOCATE(paramCorrections(numStates+1,maxHistoryLength))
   gradientCorrections = 0.0
   paramCorrections = 0.0
   istep = 0

   ! Occupation number optimization loop

   converged = .FALSE.
   DO WHILE (.NOT.converged)

      WRITE(*,*) 'RDMFT: convergence loop start'

      DO jspin = 1, input%jspins
         DO ikpt = 1,kpts%nkpt
            WRITE(*,*) 'jspin, ikpt: ', jspin, ikpt
            WRITE(*,'(8f10.5)') results%w_iks(1:10,ikpt,jspin)
         END DO
      END DO

      ! Calculate overall density with current occupation numbers (don't forget core electron density)
      CALL overallDen%resetPotDen()
      jspmax = input%jspins
      IF (noco%l_mperp) jspmax = 1
      DO jspin = 1,jspmax
         CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin)
         CALL cdnval(eig_id,mpi,kpts,jsp,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                     sphhar,sym,vTot,oneD,cdnvalJob,overallDen,regCharges,dos,results,moments)
      END DO

      CALL cdncore(mpi,dimension,oneD,input,vacuum,noco,sym,&
                   stars,cell,sphhar,atoms,vTot,overallDen,moments,results)
      IF (mpi%irank.EQ.0) THEN
         CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,overallDen,noco%l_noco,.TRUE.,.true.,fix)
      END IF
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,overallDen)
#endif

      ! Calculate Coulomb potential for overall density (+including external potential)
      CALL overallDen%sum_both_spin()!workden)
      CALL overallVCoul%resetPotDen()
      ALLOCATE(overallVCoul%pw_w(size(overallVCoul%pw,1),size(overallVCoul%pw,2)))
      overallVCoul%pw_w(:,:) = 0.0
      CALL vgen_coulomb(1,mpi,DIMENSION,oneD,input,field,vacuum,sym,stars,cell,sphhar,atoms,overallDen,overallVCoul)
      CALL convol(stars,overallVCoul%pw_w(:,1),overallVCoul%pw(:,1),stars%ufft)   ! Is there a problem with a second spin?!
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,overallVCoul)
#endif

      overallVCoulSSDen = 0.0
      DO jspin = 1, input%jspins
         jsp = MERGE(1,jspin,noco%l_noco)
         DO ikpt = 1, kpts%nkpt
            DO iBand = 1, highestState(ikpt,jsp)
               ! Read the single-state density from disc
               filename = ''
               WRITE(filename,'(a,i1.1,a,i4.4,a,i5.5)') 'cdn-', jsp, '-', ikpt, '-', iBand
               IF (mpi%irank.EQ.0) THEN
                  CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,&
                                   CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,singleStateDen,TRIM(ADJUSTL(filename)))
                  CALL singleStateDen%sum_both_spin()!workden)
               END IF
#ifdef CPP_MPI
               CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,singleStateDen)
#endif

               ! For each state calculate integral over Coulomb potential times single state density
               potDenInt = 0.0
               CALL int_nv(1,stars,vacuum,atoms,sphhar,cell,sym,input,oneD,vCoul,singleStateDen,potDenInt) ! Is there a problem with a second spin?!
               overallVCoulSSDen(iBand,ikpt,jsp) = potDenInt
            END DO
         END DO
      END DO

      ! Construct exchange matrix diagonal

      exDiag = 0.0
      DO jspin = 1, input%jspins
         jsp = MERGE(1,jspin,noco%l_noco)
         ! remove weights(wtkpt) in w_iks
         wl_iks = 0.0
         DO ikptf=1,kpts%nkptf
            ikpt = kpts%bkp(ikptf)
            DO iBand=1, results%neig(ikpt,jsp)
               wl_iks(iBand,ikptf) = results%w_iks(iBand,ikpt,jspin) / (kpts%wtkpt(ikpt))!*kpts%nkptf) Last term to be included after applying functional
               wl_iks(iBand,ikptf) = sqrt(wl_iks(iBand,ikptf)) ! For Müller functional
               wl_iks(iBand,ikptf) = wl_iks(iBand,ikptf) / kpts%nkptf ! This is the last part of two lines above
            END DO
         END DO

         IF(ALLOCATED(eig_irr)) DEALLOCATE (eig_irr)
         IF(ALLOCATED(hybdat%kveclo_eig)) DEALLOCATE (hybdat%kveclo_eig)
         IF(ALLOCATED(hybdat%pntgptd)) DEALLOCATE (hybdat%pntgptd)
         IF(ALLOCATED(hybdat%pntgpt)) DEALLOCATE (hybdat%pntgpt)
         IF(ALLOCATED(hybdat%prodm)) DEALLOCATE (hybdat%prodm)

         call hybdat%prod%free()

         IF(ALLOCATED(hybdat%nindxp1)) DEALLOCATE (hybdat%nindxp1)

         results%neig(:,:) = neigTemp(:,:)

         CALL HF_setup(mpbasis,hybrid,input,sym,kpts,dimension,atoms,mpi,noco,&
                       cell,oneD,results,jspin,enpara,eig_id,&
                       hybdat,sym%invs,vTot%mt(:,0,:,:),eig_irr)

         results%neig(:,:) = highestState(:,:) + 1

         mnobd = MAXVAL(hybrid%nobd(:,jsp))

         DO ikpt = 1,kpts%nkpt

            CALL lapw%init(input,noco,kpts,atoms,sym,ikpt,cell,l_zref)

            parent = 0
            CALL symm_hf_init(sym,kpts,ikpt,nsymop,rrot,psym)
            CALL symm_hf(kpts,ikpt,sym,dimension,hybdat,eig_irr,atoms,hybrid,cell,lapw,jspin,&
                         rrot,nsymop,psym,nkpt_EIBZ,n_q,parent,pointer_EIBZ,nsest,indx_sest)

            exMat%l_real=sym%invs
            CALL exchange_valence_hf(ikpt,kpts,nkpt_EIBZ, sym,atoms,mpbasis,hybrid,cell,dimension,input,jspin,hybdat,mnobd,lapw,&
                                     eig_irr,results,pointer_EIBZ,n_q,wl_iks,xcpot,noco,nsest,indx_sest,&
                                     mpi,exMat)
            CALL exchange_vccv1(ikpt,atoms,hybrid,hybdat,dimension,jspin,lapw,nsymop,nsest,indx_sest,mpi,1.0,results,exMat)

            !Start of workaround for increased functionality of symmetrizeh (call it))

            nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)

            CALL olap%alloc(sym%invs,nbasfcn)
            CALL read_olap(olap, kpts%nkpt*(jspin-1)+ikpt)
            IF (olap%l_real) THEN
               DO i = 1, nbasfcn
                  DO j = 1, i
                     olap%data_r(i,j) = olap%data_r(j,i)
                  END DO
               END DO
            ELSE
               DO i = 1, nbasfcn
                  DO j = 1, i
                     olap%data_c(i,j) = CONJG(olap%data_c(j,i))
                  END DO
               END DO
               olap%data_c = conjg(olap%data_c)
            END IF

            CALL zMat%init(olap%l_real,nbasfcn,dimension%neigd)

            CALL read_eig(eig_id,ikpt,jspin,list=[(i,i=1,hybrid%nbands(ikpt))],neig=nbands,zmat=zMat)

!            CALL read_z(zMat,kpts%nkpt*(jspin-1)+ikpt)
            zMat%matsize2 = hybrid%nbands(ikpt) ! reduce "visible matsize" for the following computations

            CALL olap%multiply(zMat,trafo)

            CALL invtrafo%alloc(olap%l_real,hybrid%nbands(ikpt),nbasfcn)
            CALL trafo%TRANSPOSE(invtrafo)
            IF(.NOT.invtrafo%l_real) invtrafo%data_c = CONJG(invtrafo%data_c)

            DO i = 1, hybrid%nbands(ikpt)
               DO j = 1, i-1
                  IF (exMat%l_real) THEN
                     exMat%data_r(i,j)=exMat%data_r(j,i)
                  ELSE
                     exMat%data_c(i,j)=conjg(exMat%data_c(j,i))
                  END IF
               END DO
            END DO

            CALL exMat%multiply(invtrafo,tmpMat)
            CALL trafo%multiply(tmpMat,exMatLAPW)

            CALL symmetrizeh(atoms,kpts%bkf(:,ikpt),jspin,lapw,sym,hybdat%kveclo_eig,cell,nsymop,psym,exMatLAPW)

            IF (.NOT.exMatLAPW%l_real) exMatLAPW%data_c=conjg(exMatLAPW%data_c)
            zMat%matsize1=MIN(zMat%matsize1,exMatLAPW%matsize2)

            CALL exMatLAPW%multiply(zMat,tmpMat)

            DO iBand = 1, highestState(ikpt,jsp)
               IF (zMat%l_real) THEN
                  exDiag(iBand,ikpt,jspin) = dot_product(zMat%data_r(:zMat%matsize1,iband),tmpMat%data_r(:,iband))
               ELSE
                  exDiag(iBand,ikpt,jspin) = REAL(dot_product(zMat%data_c(:zMat%matsize1,iband),tmpMat%data_c(:,iband)))
               END IF
            END DO

            !End of workaround for increased functionality of symmetrizeh (call it))

!            DO iBand = 1, highestState(ikpt,jsp)
!               IF (exMat%l_real) THEN
!                  exDiag(iBand,ikpt,jspin) = exMat%data_r(iBand,iBand)
!               ELSE
!                  exDiag(iBand,ikpt,jspin) = REAL(exMat%data_c(iBand,iBand))
!               END IF
!            END DO
         END DO
      END DO

      ! Calculate total energy derivative with respect to occupations (dEdOcc)

      gradSum = 0.0
      DO ispin = 1, input%jspins
         isp = MERGE(1,ispin,noco%l_noco)
!         CALL cdnvalJob%init(mpi,input,kpts,noco,results,isp,banddos=banddos)
         DO ikpt = 1, kpts%nkpt
            DO iBand = 1, highestState(ikpt,isp)
               occStateI = results%w_iks(iBand,ikpt,isp) / (kpts%wtkpt(ikpt))!*kpts%nkptf)
               occStateI = MAX(occStateI,minOcc)
!               IF(occStateI.LT.1.0e-7) occStateI = 5.0e-4 ! This is preliminary. I have to discuss what do do here.
!               occStateI = cdnvalJob%weights(iBand,ikpt)
               rdmftFunctionalValue = 0.5*0.5*SQRT(1.0/occStateI) ! for Müller functional derivative

               exchangeTerm = - rdmftFunctionalValue * exDiag(iBand,ikpt,isp) * kpts%wtkpt(ikpt) * spinDegenFac !*kpts%nkptf

               dEdOcc(iBand,ikpt,isp) = +((spinDegenFac * results%eig(iBand,ikpt,isp)) - vTotSSDen(iBand,ikpt,isp) + &
                                              overallVCoulSSDen(iBand,ikpt,isp) + exchangeTerm)

               WRITE(*,*) 'ENERGY GRADIENT CONTRIBUTIONS'
               WRITE(*,*) 'ispin, ikpt, iBand', ispin, ikpt, iBand
               WRITE(*,*) 'results%eig(iBand,ikpt,isp)', results%eig(iBand,ikpt,isp)
               WRITE(*,*) 'vTotSSDen(iBand,ikpt,isp)', vTotSSDen(iBand,ikpt,isp)
               WRITE(*,*) 'overallVCoulSSDen(iBand,ikpt,isp)', overallVCoulSSDen(iBand,ikpt,isp)
               WRITE(*,*) 'exchangeTerm', exchangeTerm
               WRITE(*,*) 'exDiag(iBand,ikpt,isp)', exDiag(iBand,ikpt,isp)
               WRITE(*,*) 'rdmftFunctionalValue', rdmftFunctionalValue


               gradSum = gradSum + dEdOcc(iBand,ikpt,isp) ! * results%w_iks(iBand,ikpt,isp)
            END DO
         END DO
      END DO
      lagrangeMultiplier = -gradSum / numStates !(input%zelec/(2.0/REAL(input%jspins)))

   WRITE(*,*) 'lagrangeMultiplier: ', lagrangeMultiplier

      ! Optimize occupation numbers

      ALLOCATE (gradient(numStates+1))
      ALLOCATE (parameters(numStates+1))
      ALLOCATE (equalityLinCombi(numStates+1))
      ALLOCATE (minConstraints(numStates+1))
      ALLOCATE (maxConstraints(numStates+1))
      ALLOCATE (enabledConstraints(numStates+1))
      enabledConstraints = .FALSE.
      enabledConstraints(1:numStates) = .TRUE.
      minConstraints(:) = 0.0
      maxConstraints(:) = 1.0
      equalityLinCombi = 0.0
      gradient = 0.0
      parameters = 0.0
      iState = 0
      DO ispin = 1, input%jspins
         isp = MERGE(1,ispin,noco%l_noco)
         DO ikpt = 1, kpts%nkpt
            DO iBand = lowestState(ikpt,isp), highestState(ikpt,isp)
               iState = iState + 1
               occStateI = results%w_iks(iBand,ikpt,isp) / kpts%wtkpt(ikpt)
               occStateI = MAX(occStateI,minOcc)
               equalityLinCombi(iState) = kpts%wtkpt(ikpt)
               gradient(iState) = dEdOcc(iBand,ikpt,isp) + lagrangeMultiplier
               gradient(numStates+1) = gradient(numStates+1) + occStateI * kpts%wtkpt(ikpt)
               parameters(iState) = occStateI
            END DO
         END DO
      END DO
      equalityCriterion = input%zelec/(2.0/REAL(input%jspins))
      gradient(numStates+1) = gradient(numStates+1) - equalityCriterion ! This should actually always be 0.0
      parameters(numStates+1) = lagrangeMultiplier

      mixParam = 0.01 / MAXVAL(ABS(gradient(:numStates)))
      WRITE(*,*) 'mixParam: ', mixParam

      CALL bfgs_b2(numStates+1,gradient,lastGradient,minConstraints,maxConstraints,enabledConstraints,parameters,&
                   lastParameters,equalityLinCombi,equalityCriterion,maxHistoryLength,paramCorrections,&
                   gradientCorrections,iStep,mixParam,converged,convCrit)

      iState = 0
      DO ispin = 1, input%jspins
         isp = MERGE(1,ispin,noco%l_noco)
         DO ikpt = 1, kpts%nkpt
            DO iBand = lowestState(ikpt,isp), highestState(ikpt,isp)
               iState = iState + 1
               results%w_iks(iBand,ikpt,isp) = MERGE(parameters(iState) * kpts%wtkpt(ikpt),0.0,parameters(iState).GT.minOcc)
            END DO
         END DO
      END DO

      DEALLOCATE (enabledConstraints,maxConstraints,minConstraints)
      DEALLOCATE (parameters,gradient,equalityLinCombi)


   END DO ! WHILE (.NOT.converged)

   WRITE(*,*) 'RDMFT: convergence loop end'

   hybrid%l_calhf = .FALSE.

   ! Calculate final overall density

   !I think we need most of cdngen at this place so I just use cdngen
   CALL outDen%resetPotDen()
   CALL cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,DIMENSION,kpts,atoms,sphhar,stars,sym,&
               enpara,cell,noco,vTot,results,oneD,coreSpecInput,archiveType,xcpot,outDen, EnergyDen)

   ! Calculate RDMFT energy
   rdmftEnergy = 0.0
   DO ispin = 1, input%jspins
      isp = MERGE(1,ispin,noco%l_noco)
!      CALL cdnvalJob%init(mpi,input,kpts,noco,results,isp,banddos=banddos)
      DO ikpt = 1, kpts%nkpt
         DO iBand = 1, highestState(ikpt,isp)
            occStateI = results%w_iks(iBand,ikpt,isp) / (kpts%wtkpt(ikpt))!*kpts%nkptf)
            rdmftFunctionalValue = 0.5*SQRT(occStateI) ! for Müller functional

            exchangeTerm = -rdmftFunctionalValue * exDiag(iBand,ikpt,isp) * kpts%wtkpt(ikpt) * spinDegenFac !*kpts%nkptf

            rdmftEnergy = rdmftEnergy + exchangeTerm + &
                          occStateI * ((spinDegenFac*results%eig(iBand,ikpt,isp)) - vTotSSDen(iBand,ikpt,isp) + &
                                       overallVCoulSSDen(iBand,ikpt,isp))
         END DO
      END DO
   END DO

   results%neig(:,:) = neigTemp(:,:)

   WRITE(6,'(a,f20.10,a)') 'RDMFT energy: ', rdmftEnergy, ' Htr'

#endif
END SUBROUTINE rdmft

END MODULE m_rdmft
