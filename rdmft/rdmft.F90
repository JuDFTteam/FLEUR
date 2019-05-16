!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rdmft

CONTAINS

SUBROUTINE rdmft(eig_id,mpi,input,kpts,banddos,sliceplot,cell,atoms,enpara,stars,vacuum,dimension,&
                 sphhar,sym,field,vTot,oneD,noco,xcpot,hybrid,results,coreSpecInput,archiveType,outDen)

   USE m_types
   USE m_juDFT
   USE m_constants
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
   TYPE(t_oneD),          INTENT(IN)    :: oneD
   TYPE(t_noco),          INTENT(INOUT) :: noco
   TYPE(t_xcpot_inbuild), INTENT(INOUT) :: xcpot
   TYPE(t_hybrid),        INTENT(INOUT) :: hybrid
   TYPE(t_results),       INTENT(INOUT) :: results
   TYPE(t_coreSpecInput), INTENT(IN)    :: coreSpecInput
   TYPE(t_potden),        INTENT(INOUT) :: outDen

   INTEGER,               INTENT(IN)    :: eig_id
   INTEGER,               INTENT(IN)    :: archiveType
   TYPE(t_potden)                       :: EnergyDen

#ifndef CPP_OLDINTEL
   TYPE(t_cdnvalJob)                    :: cdnvalJob
   TYPE(t_potden)                       :: singleStateDen, overallDen, overallVCoul
   TYPE(t_regionCharges)                :: regCharges
   TYPE(t_dos)                          :: dos
   TYPE(t_moments)                      :: moments
   TYPE(t_mat)                          :: exMat
   TYPE(t_lapw)                         :: lapw
   TYPE(t_hybdat)                       :: hybdat
   INTEGER                              :: ikpt, iBand, jkpt, jBand
   INTEGER                              :: jspin, jspmax, jsp, isp, ispin
   INTEGER                              :: nsymop, nkpt_EIBZ, ikptf, iterHF, mnobd
   INTEGER                              :: iState, iStep, numStates, maxHistoryLength, numRelevantStates
   REAL                                 :: fix, potDenInt, fermiEnergyTemp
   REAL                                 :: rdmftFunctionalValue, occStateI
   REAL                                 :: exchangeTerm, lagrangeMultiplier, equalityCriterion
   REAL                                 :: mixParam, convCrit, rdmftEnergy
   REAL                                 :: sumOcc, tempOcc, addCharge, subCharge, addChargeWeight, subChargeWeight
   REAL, PARAMETER                      :: degenEps = 0.00001
   LOGICAL                              :: converged, l_qfix, l_restart, l_zref
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
   ! General initializations
   mixParam = 0.0001
   convCrit = 1.0e-8
   lagrangeMultiplier = 0.1 !results%ef

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
            WRITE(*,*) 'occ: ', iState, occupationVec(iState)
         END DO
      END DO
   END DO
   DEALLOCATE(occupationVec)
   WRITE(*,*) 'sumOcc, tempOcc: ', sumOcc, tempOcc

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
   ALLOCATE(vTot%pw_w(SIZE(overallDen%pw,1),1))
   DO jspin = 1, input%jspins
      CALL convol(stars,vTot%pw_w(:,jspin),vTot%pw(:,jspin),stars%ufft)
   END DO

   vTotSSDen = 0.0

   ! Calculate all single state densities
   cdnvalJob%l_evp = .FALSE.
   cdnvalJob%nkptExtended = kpts%nkpt
   ALLOCATE(cdnvalJob%noccbd(kpts%nkpt))
   ALLOCATE(cdnvalJob%nStart(kpts%nkpt))
   ALLOCATE(cdnvalJob%nEnd(kpts%nkpt))
   ALLOCATE(cdnvalJob%weights(1,kpts%nkpt))
   numStates = 0
   DO jspin = 1, input%jspins
      jsp = MERGE(1,jspin,noco%l_noco)
      DO ikpt = 1, kpts%nkpt
         DO iBand = 1, highestState(ikpt,jsp)
            numStates = numStates + 1
            ! Construct cdnvalJob object for this state
            ! (Reasonable parallelization is not yet done - should be placed over the loops enclosing this section)
            cdnvalJob%ikptStart = ikpt
            cdnvalJob%ikptIncrement = kpts%nkpt
            IF(mpi%irank.EQ.0) THEN
               cdnvalJob%noccbd = 1
               cdnvalJob%nStart = iBand
               cdnvalJob%nEnd = iBand
               cdnvalJob%weights = 0.0
               cdnvalJob%weights(1,ikpt) = 1.0
            ELSE
               cdnvalJob%noccbd = 0
               cdnvalJob%nStart = 1
               cdnvalJob%nEnd = 0
               cdnvalJob%weights = 0.0
            END IF

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
            CALL int_nv(jsp,stars,vacuum,atoms,sphhar,cell,sym,input,oneD,vTot,singleStateDen,potDenInt)
            vTotSSDen(iBand,ikpt,jsp) = potDenInt
         END DO
      END DO
   END DO

   WRITE(*,*) 'Point A reached'

   ! Initializations for exchange contributions

   IF(ALLOCATED(hybrid%ne_eig)) DEALLOCATE(hybrid%ne_eig)
   IF(ALLOCATED(hybrid%nbands)) DEALLOCATE(hybrid%nbands)
   IF(ALLOCATED(hybrid%nobd)) DEALLOCATE(hybrid%nobd)
   IF(ALLOCATED(hybrid%nbasm)) DEALLOCATE(hybrid%nbasm)
   IF(ALLOCATED(hybrid%div_vv)) DEALLOCATE(hybrid%div_vv)
   ALLOCATE(hybrid%ne_eig(kpts%nkpt),hybrid%nbands(kpts%nkpt),hybrid%nobd(kpts%nkptf))
   ALLOCATE(hybrid%nbasm(kpts%nkptf))
   ALLOCATE(hybrid%div_vv(DIMENSION%neigd,kpts%nkpt,input%jspins))
   l_restart = .FALSE.
   l_zref = (sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco)
   iterHF = 0
   hybrid%l_calhf = .TRUE.

   CALL open_hybrid_io1(DIMENSION,sym%invs)

   CALL mixedbasis(atoms,kpts,dimension,input,cell,sym,xcpot,hybrid,enpara,mpi,vTot,l_restart)

   WRITE(*,*) 'Point B reached'

   CALL open_hybrid_io2(hybrid,DIMENSION,atoms,sym%invs)

   CALL coulombmatrix(mpi,atoms,kpts,cell,sym,hybrid,xcpot,l_restart)

   WRITE(*,*) 'Point C reached'

   CALL hf_init(hybrid,kpts,atoms,input,DIMENSION,hybdat,sym%invs)

   WRITE(*,*) 'Point D reached'

   ALLOCATE(parent(kpts%nkptf))
   ALLOCATE(exDiag(dimension%neigd,ikpt,input%jspins))
   ALLOCATE(lastGradient(numStates+1))
   ALLOCATE(lastParameters(numStates+1))
   lastGradient = 0.0
   lastParameters = 0.0
   maxHistoryLength = 5!7
   ALLOCATE(gradientCorrections(numStates+1,maxHistoryLength))
   ALLOCATE(paramCorrections(numStates+1,maxHistoryLength))
   gradientCorrections = 0.0
   paramCorrections = 0.0
   istep = 0

   ! Occupation number optimization loop

   converged = .FALSE.
   DO WHILE (.NOT.converged)

   WRITE(*,*) 'Point E reached'

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
         CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin,banddos=banddos)
         CALL cdnval(eig_id,mpi,kpts,jsp,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                     sphhar,sym,vTot,oneD,cdnvalJob,overallDen,regCharges,dos,results,moments)
      END DO

      WRITE(*,*) 'overallDen%pw(1,1)', overallDen%pw(1,1)

      CALL cdncore(mpi,dimension,oneD,input,vacuum,noco,sym,&
                   stars,cell,sphhar,atoms,vTot,overallDen,moments,results)
      IF (mpi%irank.EQ.0) THEN
         CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,overallDen,noco%l_noco,.TRUE.,.true.,fix)
      END IF
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,overallDen)
#endif

   WRITE(*,*) 'Point E-1 reached'

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

   WRITE(*,*) 'Point E-2 reached'

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
               CALL int_nv(1,stars,vacuum,atoms,sphhar,cell,sym,input,oneD,overallVCoul,singleStateDen,potDenInt) ! Is there a problem with a second spin?!
               overallVCoulSSDen(iBand,ikpt,jsp) = potDenInt
            END DO
         END DO
      END DO

   WRITE(*,*) 'Point E-3 reached'

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
         IF(ALLOCATED(hybdat%prod)) DEALLOCATE (hybdat%prod)
         IF(ALLOCATED(hybdat%nindxp1)) DEALLOCATE (hybdat%nindxp1)

         CALL HF_setup(hybrid,input,sym,kpts,dimension,atoms,mpi,noco,cell,oneD,results,jspin,enpara,eig_id,&
                       hybdat,iterHF,sym%invs,vTot%mt(:,0,:,:),eig_irr)

         mnobd = MAXVAL(hybrid%nobd)

         DO ikpt = 1,kpts%nkpt

            CALL lapw%init(input,noco,kpts,atoms,sym,ikpt,cell,l_zref)

            parent = 0
            CALL symm_hf_init(sym,kpts,ikpt,nsymop,rrot,psym)
            CALL symm_hf(kpts,ikpt,sym,dimension,hybdat,eig_irr,atoms,hybrid,cell,lapw,jspin,mpi,&
                         rrot,nsymop,psym,nkpt_EIBZ,n_q,parent,pointer_EIBZ,nsest,indx_sest)

            exMat%l_real=sym%invs
            CALL exchange_valence_hf(ikpt,kpts,nkpt_EIBZ, sym,atoms,hybrid,cell,dimension,input,jspin,hybdat,mnobd,lapw,&
                                     eig_irr,results,parent,pointer_EIBZ,n_q,wl_iks,iterHF,xcpot,noco,nsest,indx_sest,&
                                     mpi,exMat)
            CALL exchange_vccv1(ikpt,atoms,hybrid,hybdat,dimension,jspin,lapw,nsymop,nsest,indx_sest,mpi,1.0,results,exMat)
            DO iBand = 1, highestState(ikpt,jsp)
               IF (exMat%l_real) THEN
                  exDiag(iBand,ikpt,jspin) = exMat%data_r(iBand,iBand)
               ELSE
                  exDiag(iBand,ikpt,jspin) = REAL(exMat%data_c(iBand,iBand))
               END IF
            END DO
         END DO
      END DO

   WRITE(*,*) 'Point E-4 reached'

      ! Calculate total energy derivative with respect to occupations (dEdOcc)

      DO ispin = 1, input%jspins
         isp = MERGE(1,ispin,noco%l_noco)
!         CALL cdnvalJob%init(mpi,input,kpts,noco,results,isp,banddos=banddos)
         DO ikpt = 1, kpts%nkpt
            DO iBand = 1, highestState(ikpt,isp)
               occStateI = results%w_iks(iBand,ikpt,isp) / (kpts%wtkpt(ikpt))!*kpts%nkptf)
               IF(occStateI.LT.1.0e-7) occStateI = 5.0e-4 ! This is preliminary. I have to discuss what do do here.
!               occStateI = cdnvalJob%weights(iBand,ikpt)
               rdmftFunctionalValue = 0.5*SQRT(1.0/occStateI) ! for Müller functional derivative

               exchangeTerm = - rdmftFunctionalValue * exDiag(iBand,ikpt,isp) * kpts%wtkpt(ikpt)!*kpts%nkptf

               dEdOcc(iBand,ikpt,isp) = +(results%eig(iBand,ikpt,isp) - vTotSSDen(iBand,ikpt,isp) + &
                                              overallVCoulSSDen(iBand,ikpt,isp) + exchangeTerm) + &
                                              lagrangeMultiplier ! lagrangeMultiplier for charge conservation
            END DO
         END DO
      END DO

   WRITE(*,*) 'Point E-5 reached'

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
               equalityLinCombi(iState) = kpts%wtkpt(ikpt)
               gradient(iState) = dEdOcc(iBand,ikpt,isp)
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
               results%w_iks(iBand,ikpt,isp) = parameters(iState) * kpts%wtkpt(ikpt)
            END DO
         END DO
      END DO

      DEALLOCATE (enabledConstraints,maxConstraints,minConstraints)
      DEALLOCATE (parameters,gradient,equalityLinCombi)


   END DO ! WHILE (.NOT.converged)

   hybrid%l_calhf = .FALSE.

   WRITE(*,*) 'Point F reached'

   ! Calculate final overall density

   WRITE(*,*) 'Point G reached'

   !I think we need most of cdngen at this place so I just use cdngen
   CALL outDen%resetPotDen()
   CALL cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,DIMENSION,kpts,atoms,sphhar,stars,sym,&
               enpara,cell,noco,vTot,results,oneD,coreSpecInput,archiveType,xcpot,outDen, EnergyDen)

   WRITE(*,*) 'Point H reached'

   ! Calculate RDMFT energy
   rdmftEnergy = 0.0
   DO ispin = 1, input%jspins
      isp = MERGE(1,ispin,noco%l_noco)
!      CALL cdnvalJob%init(mpi,input,kpts,noco,results,isp,banddos=banddos)
      DO ikpt = 1, kpts%nkpt
         DO iBand = 1, highestState(ikpt,isp)
            occStateI = results%w_iks(iBand,ikpt,isp) / (kpts%wtkpt(ikpt))!*kpts%nkptf)
            rdmftFunctionalValue = SQRT(occStateI) ! for Müller functional

            exchangeTerm = -rdmftFunctionalValue * exDiag(iBand,ikpt,isp) * kpts%wtkpt(ikpt)!*kpts%nkptf

            rdmftEnergy = rdmftEnergy + exchangeTerm + occStateI * (results%eig(iBand,ikpt,isp) - vTotSSDen(iBand,ikpt,isp) + &
                                                                    overallVCoulSSDen(iBand,ikpt,isp))
         END DO
      END DO
   END DO

   results%neig(:,:) = neigTemp(:,:)

   WRITE(6,'(a,f20.10,a)') 'RDMFT energy: ', rdmftEnergy, ' Htr'

#endif
END SUBROUTINE rdmft

END MODULE m_rdmft
