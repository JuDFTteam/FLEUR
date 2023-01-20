!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rdmft

CONTAINS

SUBROUTINE rdmft(eig_id,fmpi,fi,enpara,stars,&
                 sphhar,vTot,vCoul,nococonv,xcpot,mpdata,hybdat,&
                 results,archiveType,outDen)
   use m_types_vacdos
   use m_work_package
   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_intgr, ONLY : intgr3
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
   USE m_xmlOutput
   USE m_types_dos
   use m_calc_cmt

#endif

   IMPLICIT NONE

   TYPE(t_mpi),           INTENT(IN)    :: fmpi
   type(t_fleurinput), intent(in)    :: fi
   TYPE(t_enpara),        INTENT(INOUT) :: enpara
   TYPE(t_stars),         INTENT(IN)    :: stars
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_potden),        INTENT(INOUT) :: vTot
   TYPE(t_potden),        INTENT(INOUT) :: vCoul
   TYPE(t_nococonv),      INTENT(INOUT) :: nococonv
   TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
   TYPE(t_mpdata),        intent(inout) :: mpdata
   TYPE(t_hybdat),        INTENT(INOUT) :: hybdat
   TYPE(t_results),       INTENT(INOUT) :: results
   TYPE(t_potden),        INTENT(INOUT) :: outDen

   INTEGER,               INTENT(IN)    :: eig_id
   INTEGER,               INTENT(IN)    :: archiveType
   TYPE(t_potden)                       :: EnergyDen

#ifndef CPP_OLDINTEL
   TYPE(t_cdnvalJob)                    :: cdnvalJob
   TYPE(t_potden)                       :: singleStateDen, overallDen, overallVCoul, vTotTemp
   TYPE(t_regionCharges)                :: regCharges
   TYPE(t_dos)                          :: dos
   TYPE(t_vacdos)                       :: vacdos
   TYPE(t_moments)                      :: moments
   TYPE(t_mat)                          :: exMat, zMat, olap, trafo, invtrafo, tmpMat, exMatLAPW
   TYPE(t_lapw)                         :: lapw
   type(t_work_package)                 :: work_pack
   INTEGER                              :: ikpt, ikpt_i, iBand, jkpt, jBand, iAtom, na, itype, lh, iGrid
   INTEGER                              :: jspin, jspmax, jsp, isp, ispin, nbasfcn, nbands
   INTEGER                              :: nsymop, ikptf, iterHF, ierr
   INTEGER                              :: iState, jState, iStep, numStates, numRelevantStates, convIter
   INTEGER                              :: maxHistoryLength
   INTEGER                              :: lastGroupEnd, currentGroupEnd
   REAL                                 :: fix, potDenInt, fermiEnergyTemp, tempDistance, spinDegenFac
   REAL                                 :: rdmftFunctionalValue, occStateI, gradSum
   REAL                                 :: exchangeTerm, lagrangeMultiplier, equalityCriterion
   REAL                                 :: mixParam, rdmftEnergy, occSum
   REAL                                 :: sumOcc, addCharge, subCharge, addChargeWeight, subChargeWeight
   REAL                                 :: rhs, totz, theta, temp
   REAL                                 :: averageParam, averageGrad
   REAL, PARAMETER                      :: degenEps = 0.00001
   REAL, PARAMETER                      :: convCrit = 5.0e-6, occMixParam = 0.5
   REAL, PARAMETER                      :: minOcc = 1.0e-13, minOccB = 1.0e-5
   LOGICAL                              :: converged, l_qfix, l_restart
   CHARACTER(LEN=20)                    :: filename
   CHARACTER(LEN=20)                    :: attributes(3)

   INTEGER                              :: nsest(fi%input%neig) ! probably too large
   INTEGER                              :: rrot(3,3,fi%sym%nsym)
   INTEGER                              :: psym(fi%sym%nsym) ! Note: psym is only filled up to index nsymop
   INTEGER                              :: lowestState(fi%kpts%nkpt,fi%input%jspins)
   INTEGER                              :: highestState(fi%kpts%nkpt,fi%input%jspins)
   INTEGER                              :: neigTemp(fi%kpts%nkpt,fi%input%jspins)

   REAL                                 :: wl_iks(fi%input%neig,fi%kpts%nkptf)

   REAL                                 :: vmd(fi%atoms%ntype), zintn_r(fi%atoms%ntype), dpj(fi%atoms%jmtd), mt(fi%atoms%jmtd,fi%atoms%ntype)

   REAL, ALLOCATABLE                    :: overallVCoulSSDen(:,:,:)
   REAL, ALLOCATABLE                    :: vTotSSDen(:,:,:)
   REAL, ALLOCATABLE                    :: dEdOcc(:,:,:)

   REAL, ALLOCATABLE                    :: zintn_rSSDen(:,:,:)
   REAL, ALLOCATABLE                    :: vmdSSDen(:,:,:)


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

   INTEGER, ALLOCATABLE                 :: paramGroup(:)
   INTEGER, ALLOCATABLE                 :: indx_sest(:,:)
   INTEGER, ALLOCATABLE                 :: parent(:)
   INTEGER, ALLOCATABLE                 :: pointer_EIBZ(:)
   INTEGER, ALLOCATABLE                 :: n_q(:)
   complex, allocatable                 :: cmt_nk(:,:,:)

   LOGICAL, ALLOCATABLE                 :: enabledConstraints(:)
   type(t_hybmpi)    :: glob_mpi

   complex :: c_phase(fi%input%neig)

#endif

#ifndef CPP_OLDINTEL

   WRITE(*,*) 'entered RDMFT subroutine'

   ! General initializations
   mixParam = 0.001
   lagrangeMultiplier = 0.1 !results%ef
   spinDegenFac = 2.0 / fi%input%jspins ! This factor is used to compensate the missing second spin in non-spinpolarized calculations

   neigTemp(:,:) = results%neig(:,:)

   ! Determine which states have to be considered
   lowestState(:,:) = -1
   highestState(:,:) = -1
   DO jspin = 1, fi%input%jspins
      jsp = MERGE(1,jspin,fi%noco%l_noco)
      DO ikpt = 1, fi%kpts%nkpt
         ! determine lowest state
         DO iBand = 1, results%neig(ikpt,jsp)
            IF((results%w_iks(iBand,ikpt,jspin) / fi%kpts%wtkpt(ikpt)).LT.(1.0-fi%input%rdmftOccEps)) THEN
               lowestState(ikpt,jspin) = iBand
               EXIT
            END IF
         END DO
         lowestState(ikpt,jspin) = lowestState(ikpt,jspin) - fi%input%rdmftStatesBelow
         DO iBand = lowestState(ikpt,jspin)-1, 1, -1
            IF((results%eig(iBand+1,ikpt,jsp) - results%eig(iBand,ikpt,jsp)).GT.degenEps) THEN
               EXIT
            END IF
            lowestState(ikpt,jspin) = iBand
         END DO
         lowestState(ikpt,jspin) = MAX(lowestState(ikpt,jspin),1)

         ! determine highest state
         DO iBand = results%neig(ikpt,jsp)-1, 1, -1
            IF((results%w_iks(iBand,ikpt,jspin) / fi%kpts%wtkpt(ikpt)).GT.(0.0+fi%input%rdmftOccEps)) THEN
               highestState(ikpt,jspin) = iBand
               EXIT
            END IF
         END DO
         highestState(ikpt,jspin) = highestState(ikpt,jspin) + fi%input%rdmftStatesAbove
         IF((results%neig(ikpt,jsp)-1).LT.highestState(ikpt,jspin)) THEN
            WRITE(oUnit,*) 'Error: Not enough states calculated:'
            WRITE(oUnit,*) 'ikpt, jsp: ', ikpt, jsp
            WRITE(oUnit,*) 'highestState(ikpt,jspin): ', highestState(ikpt,jspin)
            WRITE(oUnit,*) 'results%neig(ikpt,jsp): ', results%neig(ikpt,jsp)
            CALL juDFT_error('Not enough states calculated', calledby = 'rdmft')
         END IF
         DO iBand = highestState(ikpt,jspin)+1, results%neig(ikpt,jsp)
            IF((results%eig(iBand,ikpt,jsp) - results%eig(iBand-1,ikpt,jsp)).GT.degenEps) THEN
               EXIT
            END IF
            highestState(ikpt,jspin) = iBand
         END DO
         IF(highestState(ikpt,jspin).EQ.results%neig(ikpt,jsp)) THEN
            WRITE(oUnit,*) 'Error: Highest state is degenerate:'
            WRITE(oUnit,*) 'ikpt, jsp: ', ikpt, jsp
            WRITE(oUnit,*) 'highestState(ikpt,jspin): ', highestState(ikpt,jspin)
            WRITE(oUnit,*) 'results%eig(highestState(ikpt,jspin),ikpt,jsp): ', results%eig(highestState(ikpt,jspin),ikpt,jsp)
            WRITE(oUnit,*) 'results%eig(highestState(ikpt,jspin)-1,ikpt,jsp): ', results%eig(highestState(ikpt,jspin)-1,ikpt,jsp)
            CALL juDFT_error('Highest state is degenerate', calledby = 'rdmft')
         END IF
      END DO
   END DO

   IF (ANY(results%w_iksRDMFT(:,:,:).NE.0.0)) THEN
      results%w_iks(:,:,:) = results%w_iksRDMFT(:,:,:)
   END IF
!   results%w_iksRDMFT(:,:,:) = results%w_iks(:,:,:)

   ! Move occupations of relevant states well into allowed region
   numRelevantStates = SUM(highestState(:,:)) - SUM(lowestState(:,:)) + fi%input%jspins*fi%kpts%nkpt
   ALLOCATE(occupationVec(numRelevantStates))
   occupationVec(:) = 0.0
   sumOcc = 0.0
   addCharge = 0.0
   subCharge = 0.0
   addChargeWeight = 0.0
   subChargeWeight = 0.0
   iState = 0
   DO jspin = 1, fi%input%jspins
      jsp = MERGE(1,jspin,fi%noco%l_noco)
      DO ikpt = 1, fi%kpts%nkpt
         DO iBand = lowestState(ikpt,jspin), highestState(ikpt,jspin)
            iState = iState + 1
            occupationVec(iState) = results%w_iks(iBand,ikpt,jsp) / (fi%kpts%wtkpt(ikpt))
            sumOcc = sumOcc + results%w_iks(iBand,ikpt,jsp)
            IF(occupationVec(iState).LT.minOccB) THEN
               addCharge = addCharge + (minOccB-occupationVec(iState))*fi%kpts%wtkpt(ikpt)
               addChargeWeight = addChargeWeight + fi%kpts%wtkpt(ikpt)
            END IF
            IF(occupationVec(iState).GT.(1.0-minOccB)) THEN
               subCharge = subCharge + (occupationVec(iState)-(1.0-minOccB))*fi%kpts%wtkpt(ikpt)
               subChargeWeight = subChargeWeight + fi%kpts%wtkpt(ikpt)
            END IF
         END DO
      END DO
   END DO
   iState = 0
   DO jspin = 1, fi%input%jspins
      jsp = MERGE(1,jspin,fi%noco%l_noco)
      DO ikpt = 1, fi%kpts%nkpt
         DO iBand = lowestState(ikpt,jspin), highestState(ikpt,jspin)
            iState = iState + 1
            IF(occupationVec(iState).LT.minOccB) THEN
               occupationVec(iState) = occupationVec(iState) + 0.5*(subCharge+addCharge)*(fi%kpts%wtkpt(ikpt)/addChargeWeight)
            END IF
            IF(occupationVec(iState).GT.(1.0-minOccB)) THEN
               occupationVec(iState) = occupationVec(iState) - 0.5*(subCharge+addCharge)*(fi%kpts%wtkpt(ikpt)/subChargeWeight)
            END IF
!            results%w_iks(iBand,ikpt,jsp) = occupationVec(iState) * fi%kpts%wtkpt(ikpt)
         END DO
      END DO
   END DO
   DEALLOCATE(occupationVec)

   ! Some more initializations

   results%neig(:,:) = highestState(:,:) + 1

   ALLOCATE(overallVCoulSSDen(MAXVAL(results%neig(1:fi%kpts%nkpt,1:fi%input%jspins)),fi%kpts%nkpt,fi%input%jspins))
   ALLOCATE(vTotSSDen(MAXVAL(results%neig(1:fi%kpts%nkpt,1:fi%input%jspins)),fi%kpts%nkpt,fi%input%jspins))
   ALLOCATE(dEdOcc(MAXVAL(results%neig(1:fi%kpts%nkpt,1:fi%input%jspins)),fi%kpts%nkpt,fi%input%jspins))

   ALLOCATE(zintn_rSSDen(MAXVAL(results%neig(1:fi%kpts%nkpt,1:fi%input%jspins)),fi%kpts%nkpt,fi%input%jspins))
   ALLOCATE(vmdSSDen(MAXVAL(results%neig(1:fi%kpts%nkpt,1:fi%input%jspins)),fi%kpts%nkpt,fi%input%jspins))

   zintn_rSSDen(:,:,:) = 0.0
   vmdSSDen(:,:,:) = 0.0

   CALL regCharges%init(fi%input,fi%atoms)
   CALL dos%init(fi%input,fi%atoms,fi%kpts,fi%banddos,results%eig)
   CALL vacdos%init(fi%input,fi%atoms,fi%kpts,fi%banddos,results%eig)
   CALL moments%init(fmpi,fi%input,sphhar,fi%atoms)
   CALL overallDen%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,POTDEN_TYPE_DEN)
   CALL overallVCoul%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,POTDEN_TYPE_POTCOUL)
   IF (ALLOCATED(vTot%pw_w)) DEALLOCATE (vTot%pw_w)
   ALLOCATE(vTot%pw_w(SIZE(overallDen%pw,1),fi%input%jspins))
   DO jspin = 1, fi%input%jspins
      CALL convol(stars,vTot%pw_w(:,jspin),vTot%pw(:,jspin))
   END DO

   CALL vTotTemp%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,POTDEN_TYPE_POTTOT)
   vTotTemp = vTot

   DO jsp = 1,SIZE(vTot%mt,4)
      DO iAtom = 1,fi%atoms%ntype
         vTotTemp%mt(:fi%atoms%jri(iAtom),0,iAtom,jsp)  = sfp_const * vTot%mt(:fi%atoms%jri(iAtom),0,iAtom,jsp) / fi%atoms%rmsh(:fi%atoms%jri(iAtom),iAtom)
      END DO
   END DO

   vCoul%pw_w = CMPLX(0.0,0.0)
   DO jspin = 1, fi%input%jspins
      CALL convol(stars,vCoul%pw_w(:,jspin),vCoul%pw(:,jspin))
   END DO

   vTotSSDen = 0.0

   ! Calculate all single state densities

   numStates = 0
   DO jspin = 1, fi%input%jspins
      jsp = MERGE(1,jspin,fi%noco%l_noco)

      CALL cdnvalJob%init(fmpi,fi%input,fi%kpts,fi%noco,results,jsp)

      DO ikpt_i = 1, SIZE(fmpi%k_list)
         ikpt= fmpi%k_list(ikpt_i)
         DO iBand = 1, highestState(ikpt,jsp)
            numStates = numStates + 1
            ! Construct cdnvalJob object for this state
            ! (Reasonable parallelization is not yet done - should be placed over the loops enclosing this section)

            cdnvalJob%k_list=[ikpt]
!            cdnvalJob%ev_list=[iBand]
            cdnvalJob%weights(:,:) = 0.0
            cdnvalJob%weights(iBand,ikpt) = spinDegenFac

            ! Call cdnval to construct density
            WRITE(*,*) 'Note: some optional flags may have to be reset in rdmft before the cdnval call'
            WRITE(*,*) 'This is not yet implemented!'
            CALL singleStateDen%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,POTDEN_TYPE_DEN)
            CALL cdnval(eig_id,fmpi,fi%kpts,jsp,fi%noco,nococonv,fi%input,fi%banddos,fi%cell,fi%atoms,enpara,stars,fi%vacuum,&
                        sphhar,fi%sym,vTot, cdnvalJob,singleStateDen,regCharges,dos,vacdos,results,moments,&
                        fi%gfinp,fi%hub1inp)

            ! Store the density on disc (These are probably way too many densities to keep them in memory)
            filename = ''
            WRITE(filename,'(a,i1.1,a,i4.4,a,i5.5)') 'cdn-', jsp, '-', ikpt, '-', iBand
            IF (fmpi%irank.EQ.0) THEN
               CALL writeDensity(stars,fi%noco,fi%vacuum,fi%atoms,fi%cell,sphhar,fi%input,fi%sym, CDN_ARCHIVE_TYPE_CDN_const,CDN_input_DEN_const,&
                                 0,-1.0,0.0,-1.0,-1.0,.FALSE.,singleStateDen,inFilename=TRIM(ADJUSTL(filename)))
            END IF
            CALL singleStateDen%distribute(fmpi%mpi_comm)

            ! For each state calculate Integral over KS effective potential times single state density
            potDenInt = 0.0
            CALL int_nv(jsp,stars,fi%vacuum,fi%atoms,sphhar,fi%cell,fi%sym,fi%input,vTotTemp,singleStateDen,potDenInt)
            vTotSSDen(iBand,ikpt,jsp) = potDenInt

            mt(:,:) = 0.0
            DO iType = 1, fi%atoms%ntype
               DO iGrid = 1, fi%atoms%jri(iType)
                  mt(iGrid,iType) = singleStateDen%mt(iGrid,0,iType,jsp)
               END DO

               DO iGrid = 1, fi%atoms%jri(iType)
                  dpj(iGrid) = mt(iGrid,iType)/fi%atoms%rmsh(iGrid,iType)
               END DO

               CALL intgr3(dpj,fi%atoms%rmsh(1,iType),fi%atoms%dx(iType),fi%atoms%jri(iType),rhs)

               zintn_r(iType) = fi%atoms%neq(iType)*fi%atoms%zatom(iType)*sfp_const*rhs/2.0
               zintn_rSSDen(iBand,ikpt,jsp) = zintn_rSSDen(iBand,ikpt,jsp) + zintn_r(iType)

               CALL intgr3(mt(1,iType),fi%atoms%rmsh(1,iType),fi%atoms%dx(iType),fi%atoms%jri(iType),totz)

!               vmd(iType) = fi%atoms%rmt(iType)*vCoul%mt(fi%atoms%jri(iType),0,iType,1)/sfp_const + fi%atoms%zatom(iType) - totz*sfp_const
               vmd(iType) = -totz*sfp_const
               vmd(iType) = -fi%atoms%neq(iType)*fi%atoms%zatom(iType)*vmd(iType)/ (2.0*fi%atoms%rmt(iType))

               vmdSSDen(iBand,ikpt,jsp) = vmdSSDen(iBand,ikpt,jsp) + vmd(iType)
            END DO

         END DO
      END DO
   END DO

   ! Initializations for exchange contributions

   WRITE(*,*) 'RDMFT: HF initializations start'

   IF(ALLOCATED(hybdat%nbands)) DEALLOCATE(hybdat%nbands)
   IF(ALLOCATED(hybdat%nobd)) DEALLOCATE(hybdat%nobd)
   IF(ALLOCATED(hybdat%nbasm)) DEALLOCATE(hybdat%nbasm)
   IF(ALLOCATED(hybdat%div_vv)) DEALLOCATE(hybdat%div_vv)
   ALLOCATE(hybdat%nbasm(fi%kpts%nkptf))
   ALLOCATE(hybdat%div_vv(fi%input%neig,fi%kpts%nkpt,fi%input%jspins))

   iterHF = 0
   hybdat%l_calhf = .TRUE.

   CALL mixedbasis(fi%atoms,fi%kpts,fi%input,fi%cell,xcpot,fi%mpinp,mpdata,fi%hybinp, hybdat,enpara,fmpi,vTot, iterHF)

   !allocate coulomb matrix
   IF (.NOT.ALLOCATED(hybdat%coul)) ALLOCATE(hybdat%coul(fi%kpts%nkpt))
   DO ikpt = 1, fi%kpts%nkpt
      CALL hybdat%coul(ikpt)%alloc(fi, mpdata%num_radbasfn, mpdata%n_g, ikpt)
   END DO

   CALL glob_mpi%copy_mpi(fmpi)
   call work_pack%init(fi, hybdat, mpdata, glob_mpi, jsp, glob_mpi%rank, glob_mpi%size)

   CALL coulombmatrix(fmpi, fi, mpdata, hybdat, xcpot)

   DO ikpt = 1, fi%kpts%nkpt
      CALL hybdat%coul(ikpt)%mpi_bcast(fi, fmpi%mpi_comm, 0)
   END DO

   CALL hf_init(mpdata,fi,hybdat)

   WRITE(*,*) 'RDMFT: HF initializations end'

   maxHistoryLength = 7*numStates
   maxHistoryLength = 5*numStates
   maxHistoryLength = 23

   ALLOCATE(parent(fi%kpts%nkptf))
   ALLOCATE(exDiag(fi%input%neig,ikpt,fi%input%jspins))
   ALLOCATE(lastGradient(numStates+1))
   ALLOCATE(lastParameters(numStates+1))
   lastGradient = 0.0
   lastParameters = 0.0
   ALLOCATE(gradientCorrections(numStates+1,maxHistoryLength))
   ALLOCATE(paramCorrections(numStates+1,maxHistoryLength))
   gradientCorrections = 0.0
   paramCorrections = 0.0
   istep = 0

   ! Occupation number optimization loop

   convIter = 0

   converged = .FALSE.
   DO WHILE (.NOT.converged)

      WRITE(*,*) 'RDMFT: convergence loop start'
      convIter = convIter + 1
      WRITE(*,'(a,i7)') 'convIter: ', convIter

      DO jspin = 1, fi%input%jspins
         DO ikpt = 1,fi%kpts%nkpt
            WRITE(*,*) 'jspin, ikpt: ', jspin, ikpt
            WRITE(*,'(10f12.7)') results%w_iks(1:10,ikpt,jspin)
         END DO
      END DO

      ! Calculate overall density with current occupation numbers (don't forget core electron density)
      CALL overallDen%resetPotDen()
      jspmax = fi%input%jspins
      IF (fi%noco%l_mperp) jspmax = 1

      DO jspin = 1,jspmax
         CALL cdnvalJob%init(fmpi,fi%input,fi%kpts,fi%noco,results,jspin)
         CALL cdnval(eig_id,fmpi,fi%kpts,jspin,fi%noco,nococonv,fi%input,fi%banddos,fi%cell,fi%atoms,enpara,stars,fi%vacuum,&
                     sphhar,fi%sym,vTot, cdnvalJob,overallDen,regCharges,dos,vacdos,results,moments,&
                     fi%gfinp,fi%hub1inp)
      END DO

      CALL cdncore(fmpi, fi%input,fi%vacuum,fi%noco,nococonv,fi%sym,&
                   stars,fi%cell,sphhar,fi%atoms,vTot,overallDen,moments,results)
      IF (fmpi%irank.EQ.0) THEN
         CALL qfix(fmpi,stars,nococonv,fi%atoms,fi%sym,fi%vacuum,sphhar,fi%input,fi%cell, overallDen,&
                   fi%noco%l_noco,.TRUE.,l_par=.FALSE.,force_fix=.TRUE.,fix=fix)
      END IF
      CALL overallDen%distribute(fmpi%mpi_comm)

      ! Calculate Coulomb potential for overall density (+including external potential)
      CALL overallDen%sum_both_spin()!workden)
      CALL overallVCoul%resetPotDen()
      ALLOCATE(overallVCoul%pw_w(size(overallVCoul%pw,1),size(overallVCoul%pw,2)))
      overallVCoul%pw_w(:,:) = 0.0
      CALL vgen_coulomb(1,fmpi, fi%input,fi%field,fi%vacuum,fi%sym,stars,fi%cell,sphhar,fi%atoms,.FALSE.,overallDen,overallVCoul)
      CALL convol(stars,overallVCoul%pw_w(:,1),overallVCoul%pw(:,1))   ! Is there a problem with a second spin?!
      CALL overallVCoul%distribute(fmpi%mpi_comm)

      overallVCoulSSDen = 0.0
      DO jspin = 1, fi%input%jspins
         jsp = MERGE(1,jspin,fi%noco%l_noco)
         DO ikpt = 1, fi%kpts%nkpt
            DO iBand = 1, highestState(ikpt,jsp)
               ! Read the single-state density from disc
               filename = ''
               WRITE(filename,'(a,i1.1,a,i4.4,a,i5.5)') 'cdn-', jsp, '-', ikpt, '-', iBand
               IF (fmpi%irank.EQ.0) THEN
                  CALL readDensity(stars,fi%noco,fi%vacuum,fi%atoms,fi%cell,sphhar,fi%input,fi%sym, CDN_ARCHIVE_TYPE_CDN_const,&
                                   CDN_input_DEN_const,0,fermiEnergyTemp,tempDistance,l_qfix,singleStateDen,inFilename=TRIM(ADJUSTL(filename)))
                  CALL singleStateDen%sum_both_spin()!workden)
               END IF
               CALL singleStateDen%distribute(fmpi%mpi_comm)

               ! For each state calculate integral over Coulomb potential times single state density
               potDenInt = 0.0
               CALL int_nv(1,stars,fi%vacuum,fi%atoms,sphhar,fi%cell,fi%sym,fi%input, vCoul,singleStateDen,potDenInt) ! Is there a problem with a second spin?!
               overallVCoulSSDen(iBand,ikpt,jsp) = potDenInt
            END DO
         END DO
      END DO

      ! Construct exchange matrix diagonal

      exDiag = 0.0
      DO jspin = 1, fi%input%jspins
         jsp = MERGE(1,jspin,fi%noco%l_noco)
         ! remove weights(wtkpt) in w_iks
         wl_iks = 0.0
         DO ikptf = 1, fi%kpts%nkptf
            ikpt = fi%kpts%bkp(ikptf)
            DO iBand=1, results%neig(ikpt,jsp)
               wl_iks(iBand,ikptf) = results%w_iks(iBand,ikpt,jspin) / (fi%kpts%wtkpt(ikpt))!*fi%kpts%nkptf) Last term to be included after applying functional
               wl_iks(iBand,ikptf) = sqrt(wl_iks(iBand,ikptf)) ! For Müller functional
               wl_iks(iBand,ikptf) = wl_iks(iBand,ikptf) / fi%kpts%nkptf ! This is the last part of two lines above
            END DO
         END DO

         IF(ALLOCATED(eig_irr)) DEALLOCATE (eig_irr)
         IF(ALLOCATED(hybdat%pntgptd)) DEALLOCATE (hybdat%pntgptd)
         IF(ALLOCATED(hybdat%pntgpt)) DEALLOCATE (hybdat%pntgpt)
         IF(ALLOCATED(hybdat%prodm)) DEALLOCATE (hybdat%prodm)
         IF(ALLOCATED(hybdat%nindxp1)) DEALLOCATE (hybdat%nindxp1)

         results%neig(:,:) = neigTemp(:,:)

         CALL HF_setup(mpdata,fi,fmpi,nococonv,results,jspin,enpara,&
                       hybdat,vTot%mt(:,0,:,:),eig_irr)

         results%neig(:,:) = highestState(:,:) + 1

         DO ikpt = 1, fi%kpts%nkpt

            CALL lapw%init(fi%input,fi%noco,nococonv,fi%kpts,fi%atoms,fi%sym,ikpt,fi%cell)

            nbasfcn = 0
            IF(fi%noco%l_noco) then
               nbasfcn = lapw%nv(1) + lapw%nv(2) + 2*fi%atoms%nlotot
            ELSE
               nbasfcn = lapw%nv(1) + fi%atoms%nlotot
            END IF

            parent = 0
            CALL zMat%init(fi%sym%invs,nbasfcn,fi%input%neig)

            if(ikpt /= fi%kpts%bkp(ikpt)) call juDFT_error("We should be reading the parent z-mat here!")
            call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv,  fi%input, ikpt, jsp, zMat, c_phase=c_phase)
            allocate(cmt_nk(hybdat%nbands(ikpt,jsp), hybdat%maxlmindx, fi%atoms%nat), stat=ierr)
            if(ierr  /= 0) call judft_error("can't allocate cmt_nk")
            call calc_cmt(fi%atoms, fi%cell, fi%input, fi%noco, nococonv, fi%hybinp, hybdat, mpdata, fi%kpts, &
                        fi%sym,   zMat, jsp, ikpt, c_phase, cmt_nk)

            ALLOCATE (indx_sest(hybdat%nbands(ikpt,jsp), hybdat%nbands(ikpt,jsp)))
            indx_sest = 0

            call symm_hf_init(fi,ikpt,nsymop,rrot,psym)
            call symm_hf(fi,ikpt,hybdat,results,work_pack%k_packs(ikpt)%submpi, eig_irr,mpdata, c_phase,&
                         rrot,nsymop,psym,n_q,parent,nsest,indx_sest, jsp)

            exMat%l_real=fi%sym%invs
            CALL exchange_valence_hf(work_pack%k_packs(ikpt),fi,fmpi,zMat, mpdata,jspin,hybdat,lapw,&
                                     eig_irr,results,n_q,wl_iks,xcpot,nococonv,stars,nsest,indx_sest,&
                                     cmt_nk, exMat)
            deallocate(cmt_nk)
            CALL exchange_vccv1(ikpt,fi, mpdata,hybdat,jspin,lapw,glob_mpi,nsymop,nsest,indx_sest,&
                                1.0,results,cmt_nk,exMat)

            DEALLOCATE(indx_sest)

            !Start of workaround for increased functionality of fi%symmetrizeh (call it))

            nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

            CALL read_eig(hybdat%eig_id,ikpt,jspin, smat=olap)

            zMat%matsize2 = hybdat%nbands(ikpt,jsp) ! reduce "visible matsize" for the following computations

            CALL olap%multiply(zMat,trafo)

            CALL invtrafo%alloc(olap%l_real,hybdat%nbands(ikpt,jsp),nbasfcn)
            CALL trafo%TRANSPOSE(invtrafo)

            DO iBand = 1, hybdat%nbands(ikpt,jsp)
               DO jBand = 1, iBand-1
                  IF (exMat%l_real) THEN
                     exMat%data_r(iBand,jBand)=exMat%data_r(jBand,iBand)
                  ELSE
                     exMat%data_c(iBand,jBand)=conjg(exMat%data_c(jBand,iBand))
                  END IF
               END DO
            END DO

            CALL exMat%multiply(invtrafo,tmpMat)
            CALL trafo%multiply(tmpMat,exMatLAPW)

            call symmetrizeh(fi%atoms,fi%kpts%bkf(:,ikpt),jspin,lapw,fi%sym,fi%cell,nsymop,psym,exMatLAPW)

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

            !End of workaround for increased functionality of fi%symmetrizeh (call it))

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

      occSum = 0.0
      gradSum = 0.0
      DO ispin = 1, fi%input%jspins
         isp = MERGE(1,ispin,fi%noco%l_noco)
!         CALL cdnvalJob%init(fmpi,fi%input,fi%kpts,fi%noco,results,isp,fi%banddos=fi%banddos)
         DO ikpt = 1, fi%kpts%nkpt
            DO iBand = 1, highestState(ikpt,isp)
               occStateI = results%w_iks(iBand,ikpt,isp) / (fi%kpts%wtkpt(ikpt))!*fi%kpts%nkptf)
               occStateI = MAX(occStateI,minOcc)

               occStateI = MIN(occStateI,1.0-minOcc)

!               IF(occStateI.LT.1.0e-7) occStateI = 5.0e-4 ! This is preliminary. I have to discuss what do do here.
!               occStateI = cdnvalJob%weights(iBand,ikpt)
               rdmftFunctionalValue = 0.5*0.5*SQRT(1.0/occStateI) ! for Müller functional derivative


               !!! Test start
!                  occStateI = MIN(occStateI,1.0-minOcc)
!                  rdmftFunctionalValue = 0.5 * 0.5*pi_const*COS(0.5*pi_const*occStateI)

!                  rdmftFunctionalValue = ASIN(SQRT(occStateI)) * 2.0 / pi_const
!                  rdmftFunctionalValue = (SIN(rdmftFunctionalValue*pi_const/2.0)) * (COS(rdmftFunctionalValue*pi_const/2.0))**2.0 * pi_const
               !!! Test 2:
!                  occStateI = MIN(occStateI,1.0-minOcc)
!                  rdmftFunctionalValue = ASIN(SQRT(occStateI)) * 2.0 / pi_const
!                  rdmftFunctionalValue = COS(rdmftFunctionalValue*pi_const/2.0) * pi_const / 4.0
               !!! Test end

               occSum = occSum + results%w_iks(iBand,ikpt,isp)

               exchangeTerm = - rdmftFunctionalValue * exDiag(iBand,ikpt,isp) * fi%kpts%wtkpt(ikpt) * spinDegenFac !*fi%kpts%nkptf

               dEdOcc(iBand,ikpt,isp) = +((spinDegenFac * results%eig(iBand,ikpt,isp)) - vTotSSDen(iBand,ikpt,isp) + &
                                          overallVCoulSSDen(iBand,ikpt,isp) - &
                                          zintn_rSSDen(iBand,ikpt,isp) + vmdSSDen(iBand,ikpt,isp) + exchangeTerm )

               theta = ASIN(SQRT(occStateI))! * 2.0 /  pi_const
               dEdOcc(iBand,ikpt,isp) = 2.0 * sin(theta) * cos(theta) * dEdOcc(iBand,ikpt,isp)
!               dEdOcc(iBand,ikpt,isp) = dEdOcc(iBand,ikpt,isp) + 2.0 * COS(theta) * exDiag(iBand,ikpt,isp) * fi%kpts%wtkpt(ikpt) * spinDegenFac

               WRITE(*,*) 'ENERGY GRADIENT CONTRIBUTIONS'
               WRITE(*,*) 'ispin, ikpt, iBand, weight', ispin, ikpt, iBand, results%w_iks(iBand,ikpt,isp)
               WRITE(*,*) 'results%eig(iBand,ikpt,isp)', results%eig(iBand,ikpt,isp)
               WRITE(*,*) 'vTotSSDen(iBand,ikpt,isp)', vTotSSDen(iBand,ikpt,isp)
               WRITE(*,*) 'overallVCoulSSDen(iBand,ikpt,isp)', overallVCoulSSDen(iBand,ikpt,isp)
               WRITE(*,*) 'zintn_rSSDen(iBand,ikpt,isp)', zintn_rSSDen(iBand,ikpt,isp)
               WRITE(*,*) 'vmdSSDen(iBand,ikpt,isp)', vmdSSDen(iBand,ikpt,isp)
               WRITE(*,*) 'exchangeTerm', exchangeTerm
               WRITE(*,*) 'exDiag(iBand,ikpt,isp)', exDiag(iBand,ikpt,isp)
               WRITE(*,*) 'rdmftFunctionalValue', rdmftFunctionalValue

               gradSum = gradSum + dEdOcc(iBand,ikpt,isp) !* results%w_iks(iBand,ikpt,isp))
            END DO
         END DO
      END DO
      lagrangeMultiplier = -gradSum / numStates ! occSum  !(fi%input%zelec/(2.0/REAL(fi%input%jspins)))

   WRITE(*,*) 'lagrangeMultiplier: ', lagrangeMultiplier

      ! Optimize occupation numbers

      ALLOCATE (gradient(numStates+1))
      ALLOCATE (parameters(numStates+1))
      ALLOCATE (paramGroup(numStates+1))
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
      DO ispin = 1, fi%input%jspins
         isp = MERGE(1,ispin,fi%noco%l_noco)
         DO ikpt = 1, fi%kpts%nkpt
            DO iBand = lowestState(ikpt,isp), highestState(ikpt,isp)
               iState = iState + 1
               occStateI = results%w_iks(iBand,ikpt,isp) / fi%kpts%wtkpt(ikpt)

               occStateI = MAX(occStateI,minOcc)
               occStateI = MIN(occStateI,1.0-minOcc)

               theta = ASIN(SQRT(occStateI))! * 2.0 /  pi_const

               WRITE(7865,'(i7,4f15.10)') iState, occStateI, theta, sin(theta), cos(theta)

!               occStateI = MAX(occStateI,minOcc)
               equalityLinCombi(iState) = fi%kpts%wtkpt(ikpt)

!               dEdOcc(iBand,ikpt,isp) = dEdOcc(iBand,ikpt,isp) + lagrangeMultiplier

!               dEdOcc(iBand,ikpt,isp) = 2.0 * sin(theta) * cos(theta) * dEdOcc(iBand,ikpt,isp)

               gradient(iState) = dEdOcc(iBand,ikpt,isp) + lagrangeMultiplier

               gradient(numStates+1) = gradient(numStates+1) + occStateI * fi%kpts%wtkpt(ikpt)
!               parameters(iState) = theta !occStateI
               parameters(iState) = occStateI
               IF(iBand.EQ.1) THEN
                  IF(iState.EQ.1) THEN
                     paramGroup(iState) = 1
                  ELSE
                     paramGroup(iState) = paramGroup(iState-1) + 1
                  END IF
               ELSE
                  IF ((results%eig(iBand,ikpt,jsp) - results%eig(iBand-1,ikpt,jsp)).GT.degenEps) THEN
                     paramGroup(iState) = paramGroup(iState-1) + 1
                  ELSE
                     paramGroup(iState) = paramGroup(iState-1)
                  END IF
               END IF
            END DO
         END DO
      END DO
      equalityCriterion = fi%input%zelec/(2.0/REAL(fi%input%jspins))
      gradient(numStates+1) = gradient(numStates+1) - equalityCriterion ! This should actually always be 0.0
      parameters(numStates+1) = lagrangeMultiplier
      paramGroup(numStates+1) = paramGroup(numStates) + 1

      currentGroupEnd = 0
      lastGroupEnd = 0
      DO iState = 2, numStates + 1
         IF (paramGroup(iState).NE.paramGroup(iState-1)) THEN
            currentGroupEnd = iState - 1
            averageParam = 0.0
            averageGrad = 0.0
            DO jState = lastGroupEnd + 1, currentGroupEnd
               averageParam = averageParam + parameters(jState)
               averageGrad = averageGrad + gradient(jState)
            END DO
            averageParam = averageParam / (currentGroupEnd - lastGroupEnd)
            averageGrad = averageGrad / (currentGroupEnd - lastGroupEnd)
            DO jState = lastGroupEnd + 1, currentGroupEnd
               temp = ABS(parameters(jState) - averageParam) / (ABS(averageParam) + degenEps)
               IF (temp.GT.degenEps) THEN
                  WRITE(*,'(a,i0,a,f15.10,a,f15.10)') 'iState: ', jState, 'average: ', averageParam, 'parameter: ', parameters(jState)
                  CALL juDFT_error('parameter difference in paramGroup is very large (rdmft-1)')
               END IF
               temp = ABS(gradient(jState) - averageGrad) / (ABS(averageGrad) + degenEps)
               IF (temp.GT.degenEps) THEN
                  WRITE(*,'(a,i0,a,f15.10,a,f15.10)') 'iState: ', jState, 'average: ', averageGrad, 'gradient: ', gradient(jState)
                  CALL juDFT_error('Gradient difference in paramGroup is very large (rdmft-1)')
               END IF
               parameters(jState) = averageParam
               gradient(jState) = averageGrad
            END DO
            lastGroupEnd = currentGroupEnd
         END IF
      END DO

      WRITE(*,*) 'gradient(numStates+1): ', gradient(numStates+1)

!      mixParam = 0.01 / MAXVAL(ABS(gradient(:numStates)))
!      mixParam = MIN(0.0002,mixParam)
      mixParam = 0.001
      WRITE(*,*) 'mixParam: ', mixParam

      CALL bfgs_b2(numStates+1,gradient,lastGradient,minConstraints,maxConstraints,enabledConstraints,parameters,&
                   lastParameters,equalityLinCombi,equalityCriterion,maxHistoryLength,paramCorrections,&
                   gradientCorrections,iStep,mixParam,converged,convCrit)


      currentGroupEnd = 0
      lastGroupEnd = 0
      DO iState = 2, numStates + 1
         IF (paramGroup(iState).NE.paramGroup(iState-1)) THEN
            currentGroupEnd = iState - 1
            averageParam = 0.0
            averageGrad = 0.0
            DO jState = lastGroupEnd + 1, currentGroupEnd
               averageParam = averageParam + parameters(jState)
               averageGrad = averageGrad + gradient(jState)
            END DO
            averageParam = averageParam / (currentGroupEnd - lastGroupEnd)
            averageGrad = averageGrad / (currentGroupEnd - lastGroupEnd)
            DO jState = lastGroupEnd + 1, currentGroupEnd
               temp = ABS(parameters(jState) - averageParam) / (ABS(averageParam) + degenEps)
               IF (temp.GT.degenEps) THEN
                  WRITE(*,'(a,i0,a,f15.10,a,f15.10)') 'iState: ', jState, 'average: ', averageParam, 'parameter: ', parameters(jState)
                  CALL juDFT_error('parameter difference in paramGroup is very large (rdmft-2)')
               END IF
               temp = ABS(gradient(jState) - averageGrad) / (ABS(averageGrad) + degenEps)
               IF (temp.GT.degenEps) THEN
                  WRITE(*,'(a,i0,a,f15.10,a,f15.10)') 'iState: ', jState, 'average: ', averageGrad, 'gradient: ', gradient(jState)
                  CALL juDFT_error ('Gradient difference in paramGroup is very large (rdmft-2)')
               END IF
               parameters(jState) = averageParam
               gradient(jState) = averageGrad
            END DO
            lastGroupEnd = currentGroupEnd
         END IF
      END DO


      WRITE(3555,*) 'Occupation numbers:'
      iState = 0
      DO ispin = 1, fi%input%jspins
         isp = MERGE(1,ispin,fi%noco%l_noco)
         DO ikpt = 1, fi%kpts%nkpt
            DO iBand = lowestState(ikpt,isp), highestState(ikpt,isp)
               iState = iState + 1
!               parameters(iState) = (SIN(parameters(iState)*0.5*pi_const))**2.0
               WRITE(3555,'(3i7,f15.10)') iBand, ikpt, isp, parameters(iState)
               results%w_iks(iBand,ikpt,isp) = parameters(iState) * fi%kpts%wtkpt(ikpt)
!               results%w_iks(iBand,ikpt,isp) = MERGE(parameters(iState) * fi%kpts%wtkpt(ikpt),0.0,parameters(iState).GT.minOcc)
            END DO
         END DO
      END DO

      WRITE(3555,'(a,f15.10)') 'total occupation: ', SUM(parameters(:numStates))

      DEALLOCATE (enabledConstraints,maxConstraints,minConstraints)
      DEALLOCATE (parameters,gradient,paramGroup,equalityLinCombi)


   END DO ! WHILE (.NOT.converged)

   ! Only mix a part of the newly determined occupation numbers into the occupation numbers from the
   ! previous iteration.
!   WRITE(*,*) 'Test: mixing in only a part of newly determined occupation numbers!'
!   DO ispin = 1, fi%input%jspins
!      isp = MERGE(1,ispin,fi%noco%l_noco)
!      DO ikpt = 1, fi%kpts%nkpt
!         DO iBand = lowestState(ikpt,isp), highestState(ikpt,isp)
!            results%w_iks(iBand,ikpt,isp) = (1.0-occMixParam) * results%w_iksRDMFT(iBand,ikpt,isp) + &
!                                            occMixParam * results%w_iks(iBand,ikpt,isp)
!         END DO
!      END DO
!   END DO


   WRITE(2503,*) 'convIter: ', convIter

   WRITE(*,*) 'RDMFT: convergence loop end'

   hybdat%l_calhf = .FALSE.

   ! Calculate final overall density

   !I think we need most of cdngen at this place so I just use cdngen
   CALL outDen%resetPotDen()
   CALL cdngen(eig_id,fmpi,fi%input,fi%banddos,fi%sliceplot,fi%vacuum,fi%kpts,fi%atoms,sphhar,stars,fi%sym,fi%gfinp,fi%hub1inp,&
               enpara,fi%cell,fi%noco,nococonv,vTot,results, fi%corespecinput,archiveType,xcpot,outDen, EnergyDen)

   ! Calculate RDMFT energy
   rdmftEnergy = 0.0
   DO ispin = 1, fi%input%jspins
      isp = MERGE(1,ispin,fi%noco%l_noco)
!      CALL cdnvalJob%init(fmpi,fi%input,fi%kpts,fi%noco,results,isp,fi%banddos=fi%banddos)
      DO ikpt = 1, fi%kpts%nkpt
         DO iBand = 1, highestState(ikpt,isp)
            occStateI = results%w_iks(iBand,ikpt,isp) / (fi%kpts%wtkpt(ikpt))!*fi%kpts%nkptf)
            rdmftFunctionalValue = 0.5*SQRT(occStateI) ! for Müller functional

            exchangeTerm = -rdmftFunctionalValue * exDiag(iBand,ikpt,isp) * fi%kpts%wtkpt(ikpt) * spinDegenFac !*fi%kpts%nkptf

            rdmftEnergy = rdmftEnergy + exchangeTerm + &
                          occStateI * ((spinDegenFac*results%eig(iBand,ikpt,isp)) - vTotSSDen(iBand,ikpt,isp) + &
                                       0.5*overallVCoulSSDen(iBand,ikpt,isp))

               WRITE(2505,*) 'ENERGY CONTRIBUTIONS'
               WRITE(2505,*) 'ispin, ikpt, iBand, weight', ispin, ikpt, iBand, results%w_iks(iBand,ikpt,isp)
               WRITE(2505,*) 'results%eig(iBand,ikpt,isp)', occStateI * spinDegenFac * results%eig(iBand,ikpt,isp)
               WRITE(2505,*) 'vTotSSDen(iBand,ikpt,isp)', - occStateI * vTotSSDen(iBand,ikpt,isp)
               WRITE(2505,*) 'overallVCoulSSDen(iBand,ikpt,isp)', occStateI * overallVCoulSSDen(iBand,ikpt,isp)
               WRITE(2505,*) 'exchangeTerm', exchangeTerm
               WRITE(2505,*) 'exDiag(iBand,ikpt,isp)', exDiag(iBand,ikpt,isp)
               WRITE(2505,*) 'rdmftFunctionalValue', rdmftFunctionalValue

         END DO
      END DO
   END DO

   results%w_iksRDMFT(:,:,:) = results%w_iks(:,:,:)
   results%neig(:,:) = neigTemp(:,:)


   ! Madelung term (taken from totale):

   mt=0.0
   DO iType = 1, fi%atoms%ntype
      DO iGrid = 1, fi%atoms%jri(iType)
         mt(iGrid,iType) = outDen%mt(iGrid,0,iType,1) + outDen%mt(iGrid,0,iType,fi%input%jspins)
      END DO
   END DO
   IF (fi%input%jspins.EQ.1) mt=mt/2.0 !we just added the same value twice

   DO iType = 1, fi%atoms%ntype
      DO iGrid = 1,fi%atoms%jri(iType)
         dpj(iGrid) = mt(iGrid,iType)/fi%atoms%rmsh(iGrid,iType)
      END DO
      CALL intgr3(dpj,fi%atoms%rmsh(1,iType),fi%atoms%dx(iType),fi%atoms%jri(iType),rhs)

      zintn_r(iType) = fi%atoms%neq(iType)*fi%atoms%zatom(iType)*sfp_const*rhs/2.
      rdmftEnergy = rdmftEnergy - zintn_r(iType)

      CALL intgr3(mt(1,iType),fi%atoms%rmsh(1,iType),fi%atoms%dx(iType),fi%atoms%jri(iType),totz)

      vmd(iType) = fi%atoms%rmt(iType)*vCoul%mt(fi%atoms%jri(iType),0,iType,1)/sfp_const + fi%atoms%zatom(iType) - totz*sfp_const
      vmd(iType) = -fi%atoms%neq(iType)*fi%atoms%zatom(iType)*vmd(iType)/ (2.0*fi%atoms%rmt(iType))

      rdmftEnergy = rdmftEnergy + vmd(iType)

      WRITE(2505,*) '======================================='
      WRITE(2505,*) 'iType: ', iType
      WRITE(2505,*) 'zintn_r(iType): ', zintn_r(iType)
      WRITE(2505,*) 'vmd(iType): ', vmd(iType)
      WRITE(2505,*) '======================================='

   END DO

   ! Output
   WRITE(oUnit,'(a,f20.10,a)') 'RDMFT energy: ', rdmftEnergy, ' Htr'
   CALL openXMLElementPoly('rdmft',(/'energy'/),(/rdmftEnergy/))
   DO ispin = 1, fi%input%jspins
      isp = MERGE(1,ispin,fi%noco%l_noco)
      DO ikpt = 1, fi%kpts%nkpt
         CALL openXMLElementPoly('occupations',(/'spin  ','kpoint'/),(/ispin,ikpt/))
         DO iBand = lowestState(ikpt,isp), highestState(ikpt,isp)
            occStateI = results%w_iks(iBand,ikpt,isp) / (fi%kpts%wtkpt(ikpt))!*fi%kpts%nkptf)
            WRITE(attributes(1),'(i0)') iBand
            WRITE(attributes(2),'(f18.10)') results%eig(iBand,ikpt,isp)
            WRITE(attributes(3),'(f18.10)') occStateI
            CALL writeXMLElement('state',(/'index     ','energy    ','occupation'/),attributes)
         END DO
         CALL closeXMLElement('occupations')
      END DO
   END DO
   CALL closeXMLElement('rdmft')

   WRITE(2505,*) '======================================='
   WRITE(2505,*) 'convIter: ', convIter
   WRITE(2505,'(a,f20.10,a)') 'RDMFT energy: ', rdmftEnergy, ' Htr'
   WRITE(2505,*) '======================================='
   WRITE(2505,*) '======================================='
   WRITE(2505,*) '======================================='

#endif
END SUBROUTINE rdmft

END MODULE m_rdmft
