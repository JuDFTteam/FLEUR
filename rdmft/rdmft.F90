!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rdmft

CONTAINS

SUBROUTINE rdmft(eig_id,mpi,fi,enpara,stars,&
                 sphhar,vTot,vCoul,nococonv,xcpot,mpdata,hybdat,&
                 results,archiveType,outDen)

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
   USE m_io_hybinp
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
   TYPE(t_moments)                      :: moments
   TYPE(t_mat)                          :: exMat, zMat, olap, trafo, invtrafo, tmpMat, exMatLAPW
   TYPE(t_lapw)                         :: lapw
   INTEGER                              :: ikpt, ikpt_i, iBand, jkpt, jBand, iAtom, i, na, itype, lh, j
   INTEGER                              :: jspin, jspmax, jsp, isp, ispin, nbasfcn, nbands
   INTEGER                              :: nsymop, nkpt_EIBZ, ikptf, iterHF
   INTEGER                              :: iState, iStep, numStates, numRelevantStates, convIter
   INTEGER                              :: maxHistoryLength
   REAL                                 :: fix, potDenInt, fermiEnergyTemp, spinDegenFac
   REAL                                 :: rdmftFunctionalValue, occStateI, gradSum
   REAL                                 :: exchangeTerm, lagrangeMultiplier, equalityCriterion
   REAL                                 :: mixParam, rdmftEnergy, occSum
   REAL                                 :: sumOcc, addCharge, subCharge, addChargeWeight, subChargeWeight
   REAL                                 :: rhs, totz, theta
   REAL, PARAMETER                      :: degenEps = 0.00001
   REAL, PARAMETER                      :: convCrit = 5.0e-6
   REAL, PARAMETER                      :: minOcc = 1.0e-13
   LOGICAL                              :: converged, l_qfix, l_restart, l_zref
   CHARACTER(LEN=20)                    :: filename

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

   INTEGER, ALLOCATABLE                 :: indx_sest(:,:)
   INTEGER, ALLOCATABLE                 :: parent(:)
   INTEGER, ALLOCATABLE                 :: pointer_EIBZ(:)
   INTEGER, ALLOCATABLE                 :: n_q(:)

   LOGICAL, ALLOCATABLE                 :: enabledConstraints(:)
   type(t_hybmpi)    :: hybmpi

   complex :: c_phase(fi%input%neig)

#endif

#ifndef CPP_OLDINTEL

   WRITE(*,*) 'entered RDMFT subroutine'

   ! General initializations
   mixParam = 0.0001
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

   IF (ANY(results%w_iksRDMFT(:,:,:).NE.0.0)) THEN
      results%w_iks(:,:,:) = results%w_iksRDMFT(:,:,:)
   END IF

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
            IF(occupationVec(iState).LT.0.0001) THEN
               addCharge = addCharge + (0.0001-occupationVec(iState))*fi%kpts%wtkpt(ikpt)
               addChargeWeight = addChargeWeight + fi%kpts%wtkpt(ikpt)
            END IF
            IF(occupationVec(iState).GT.0.9999) THEN
               subCharge = subCharge + (occupationVec(iState)-0.9999)*fi%kpts%wtkpt(ikpt)
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
            IF(occupationVec(iState).LT.0.0001) THEN
               occupationVec(iState) = occupationVec(iState) + 0.5*(subCharge+addCharge)*(fi%kpts%wtkpt(ikpt)/addChargeWeight)
            END IF
            IF(occupationVec(iState).GT.0.9999) THEN
               occupationVec(iState) = occupationVec(iState) - 0.5*(subCharge+addCharge)*(fi%kpts%wtkpt(ikpt)/subChargeWeight)
            END IF
            results%w_iks(iBand,ikpt,jsp) = occupationVec(iState) * fi%kpts%wtkpt(ikpt)
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
   CALL dos%init(fi%input,fi%atoms,fi%kpts,fi%vacuum)
   CALL moments%init(mpi,fi%input,sphhar,fi%atoms)
   CALL overallDen%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,POTDEN_TYPE_DEN)
   CALL overallVCoul%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,POTDEN_TYPE_POTCOUL)
   IF (ALLOCATED(vTot%pw_w)) DEALLOCATE (vTot%pw_w)
   ALLOCATE(vTot%pw_w(SIZE(overallDen%pw,1),fi%input%jspins))
   DO jspin = 1, fi%input%jspins
      CALL convol(stars,vTot%pw_w(:,jspin),vTot%pw(:,jspin),stars%ufft)
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
      CALL convol(stars,vCoul%pw_w(:,jspin),vCoul%pw(:,jspin),stars%ufft)
   END DO

   vTotSSDen = 0.0

   ! Calculate all single state densities

   numStates = 0
   DO jspin = 1, fi%input%jspins
      jsp = MERGE(1,jspin,fi%noco%l_noco)

      CALL cdnvalJob%init(mpi,fi%input,fi%kpts,fi%noco,results,jsp)

      DO ikpt_i = 1, SIZE(mpi%k_list)
         ikpt= mpi%k_list(ikpt_i)
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
            CALL cdnval(eig_id,mpi,fi%kpts,jsp,fi%noco,nococonv,fi%input,fi%banddos,fi%cell,fi%atoms,enpara,stars,fi%vacuum,&
                        sphhar,fi%sym,vTot,fi%oned,cdnvalJob,singleStateDen,regCharges,dos,results,moments,&
                        fi%gfinp,fi%hub1inp)

            ! Store the density on disc (These are probably way too many densities to keep them in memory)
            filename = ''
            WRITE(filename,'(a,i1.1,a,i4.4,a,i5.5)') 'cdn-', jsp, '-', ikpt, '-', iBand
            IF (mpi%irank.EQ.0) THEN
               CALL writeDensity(stars,fi%noco,fi%vacuum,fi%atoms,fi%cell,sphhar,fi%input,fi%sym,fi%oned,CDN_ARCHIVE_TYPE_CDN_const,CDN_input_DEN_const,&
                                 0,-1.0,0.0,.FALSE.,singleStateDen,TRIM(ADJUSTL(filename)))
            END IF
#ifdef CPP_MPI
            CALL mpi_bc_potden(mpi,stars,sphhar,fi%atoms,fi%input,fi%vacuum,fi%oned,fi%noco,singleStateDen)
#endif
            ! For each state calculate Integral over KS effective potential times single state density
            potDenInt = 0.0
            CALL int_nv(jsp,stars,fi%vacuum,fi%atoms,sphhar,fi%cell,fi%sym,fi%input,fi%oned,vTotTemp,singleStateDen,potDenInt)
            vTotSSDen(iBand,ikpt,jsp) = potDenInt

            mt(:,:) = 0.0
            DO iType = 1, fi%atoms%ntype
               DO i = 1, fi%atoms%jri(iType)
                  mt(i,iType) = singleStateDen%mt(i,0,iType,jsp)
               END DO

               DO j = 1,fi%atoms%jri(iType)
                  dpj(j) = mt(j,iType)/fi%atoms%rmsh(j,iType)
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

   IF(ALLOCATED(hybdat%ne_eig)) DEALLOCATE(hybdat%ne_eig)
   IF(ALLOCATED(hybdat%nbands)) DEALLOCATE(hybdat%nbands)
   IF(ALLOCATED(hybdat%nobd)) DEALLOCATE(hybdat%nobd)
   IF(ALLOCATED(hybdat%nbasm)) DEALLOCATE(hybdat%nbasm)
   IF(ALLOCATED(hybdat%div_vv)) DEALLOCATE(hybdat%div_vv)
   ALLOCATE(hybdat%ne_eig(fi%kpts%nkpt),hybdat%nbands(fi%kpts%nkpt),hybdat%nobd(fi%kpts%nkptf,fi%input%jspins))
   ALLOCATE(hybdat%nbasm(fi%kpts%nkptf))
   ALLOCATE(hybdat%div_vv(fi%input%neig,fi%kpts%nkpt,fi%input%jspins))

   l_zref = (fi%sym%zrfs.AND.(SUM(ABS(fi%kpts%bk(3,:fi%kpts%nkpt))).LT.1e-9).AND..NOT.fi%noco%l_noco)
   iterHF = 0
   hybdat%l_calhf = .TRUE.

!   CALL open_fi%hybinp_io1(fi%sym%invs)

   CALL mixedbasis(fi%atoms,fi%kpts,fi%input,fi%cell,xcpot,fi%mpinp,mpdata,fi%hybinp, hybdat,enpara,mpi,vTot, iterHF)

   !CALL open_hybinp_io2(mpdata, fi%hybinp,hybdat,fi%input,fi%atoms,fi%sym%invs)

   if(mpi%irank == 0) CALL coulombmatrix(mpi, fi, mpdata, hybdat, xcpot)
   call hybmpi%copy_mpi(mpi)
   do i =1,fi%kpts%nkpt
      call hybdat%coul(i)%mpi_ibc(fi, hybmpi)
   enddo

   CALL hf_init(eig_id,mpdata,fi,hybdat)

   WRITE(*,*) 'RDMFT: HF initializations end'

   maxHistoryLength = 7*numStates
   maxHistoryLength = 5*numStates

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
            WRITE(*,'(10f11.6)') results%w_iks(1:10,ikpt,jspin)
         END DO
      END DO

      ! Calculate overall density with current occupation numbers (don't forget core electron density)
      CALL overallDen%resetPotDen()
      jspmax = fi%input%jspins
      IF (fi%noco%l_mperp) jspmax = 1

      DO jspin = 1,jspmax
         CALL cdnvalJob%init(mpi,fi%input,fi%kpts,fi%noco,results,jspin)
         CALL cdnval(eig_id,mpi,fi%kpts,jspin,fi%noco,nococonv,fi%input,fi%banddos,fi%cell,fi%atoms,enpara,stars,fi%vacuum,&
                     sphhar,fi%sym,vTot,fi%oned,cdnvalJob,overallDen,regCharges,dos,results,moments,&
                     fi%gfinp,fi%hub1inp)
      END DO

      CALL cdncore(mpi,fi%oned,fi%input,fi%vacuum,fi%noco,nococonv,fi%sym,&
                   stars,fi%cell,sphhar,fi%atoms,vTot,overallDen,moments,results)
      IF (mpi%irank.EQ.0) THEN
         CALL qfix(mpi,stars,fi%atoms,fi%sym,fi%vacuum,sphhar,fi%input,fi%cell,fi%oned,overallDen,fi%noco%l_noco,.TRUE.,.true.,fix)
      END IF
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,fi%atoms,fi%input,fi%vacuum,fi%oned,fi%noco,overallDen)
#endif

      ! Calculate Coulomb potential for overall density (+including external potential)
      CALL overallDen%sum_both_spin()!workden)
      CALL overallVCoul%resetPotDen()
      ALLOCATE(overallVCoul%pw_w(size(overallVCoul%pw,1),size(overallVCoul%pw,2)))
      overallVCoul%pw_w(:,:) = 0.0
      CALL vgen_coulomb(1,mpi,fi%oned,fi%input,fi%field,fi%vacuum,fi%sym,stars,fi%cell,sphhar,fi%atoms,.FALSE.,overallDen,overallVCoul)
      CALL convol(stars,overallVCoul%pw_w(:,1),overallVCoul%pw(:,1),stars%ufft)   ! Is there a problem with a second spin?!
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,fi%atoms,fi%input,fi%vacuum,fi%oned,fi%noco,overallVCoul)
#endif

      overallVCoulSSDen = 0.0
      DO jspin = 1, fi%input%jspins
         jsp = MERGE(1,jspin,fi%noco%l_noco)
         DO ikpt = 1, fi%kpts%nkpt
            DO iBand = 1, highestState(ikpt,jsp)
               ! Read the single-state density from disc
               filename = ''
               WRITE(filename,'(a,i1.1,a,i4.4,a,i5.5)') 'cdn-', jsp, '-', ikpt, '-', iBand
               IF (mpi%irank.EQ.0) THEN
                  CALL readDensity(stars,fi%noco,fi%vacuum,fi%atoms,fi%cell,sphhar,fi%input,fi%sym,fi%oned,CDN_ARCHIVE_TYPE_CDN_const,&
                                   CDN_input_DEN_const,0,fermiEnergyTemp,l_qfix,singleStateDen,TRIM(ADJUSTL(filename)))
                  CALL singleStateDen%sum_both_spin()!workden)
               END IF
#ifdef CPP_MPI
               CALL mpi_bc_potden(mpi,stars,sphhar,fi%atoms,fi%input,fi%vacuum,fi%oned,fi%noco,singleStateDen)
#endif

               ! For each state calculate integral over Coulomb potential times single state density
               potDenInt = 0.0
               CALL int_nv(1,stars,fi%vacuum,fi%atoms,sphhar,fi%cell,fi%sym,fi%input,fi%oned,vCoul,singleStateDen,potDenInt) ! Is there a problem with a second spin?!
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
         DO ikptf=1,fi%kpts%nkptf
            ikpt = fi%kpts%bkp(ikptf)
            DO iBand=1, results%neig(ikpt,jsp)
               wl_iks(iBand,ikptf) = results%w_iks(iBand,ikpt,jspin) / (fi%kpts%wtkpt(ikpt))!*fi%kpts%nkptf) Last term to be included after applying functional
               wl_iks(iBand,ikptf) = sqrt(wl_iks(iBand,ikptf)) ! For Müller functional
               wl_iks(iBand,ikptf) = wl_iks(iBand,ikptf) / fi%kpts%nkptf ! This is the last part of two lines above
            END DO
         END DO

         IF(ALLOCATED(eig_irr)) DEALLOCATE (eig_irr)
         IF(ALLOCATED(hybdat%kveclo_eig)) DEALLOCATE (hybdat%kveclo_eig)
         IF(ALLOCATED(hybdat%pntgptd)) DEALLOCATE (hybdat%pntgptd)
         IF(ALLOCATED(hybdat%pntgpt)) DEALLOCATE (hybdat%pntgpt)
         IF(ALLOCATED(hybdat%prodm)) DEALLOCATE (hybdat%prodm)
         IF(ALLOCATED(hybdat%nindxp1)) DEALLOCATE (hybdat%nindxp1)

         results%neig(:,:) = neigTemp(:,:)

         CALL HF_setup(mpdata,fi%hybinp,fi%input,fi%sym,fi%kpts,fi%atoms,mpi,fi%noco,nococonv,&
                       fi%cell,fi%oned,results,jspin,enpara,&
                       hybdat,fi%sym%invs,vTot%mt(:,0,:,:),eig_irr)

         results%neig(:,:) = highestState(:,:) + 1

         DO ikpt = 1,fi%kpts%nkpt

            CALL lapw%init(fi%input,fi%noco,nococonv,fi%kpts,fi%atoms,fi%sym,ikpt,fi%cell,l_zref)

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

            ALLOCATE (indx_sest(hybdat%nbands(ikpt), hybdat%nbands(ikpt)))
            indx_sest = 0

            call symm_hf_init(fi%sym,fi%kpts,ikpt,nsymop,rrot,psym)
            call symm_hf(fi%kpts,ikpt,fi%sym,hybdat,eig_irr,fi%input,fi%atoms,mpdata,fi%hybinp,fi%cell,lapw,&
                         fi%noco,nococonv, fi%oned, zMat, c_phase,jspin,&
                         rrot,nsymop,psym,nkpt_EIBZ,n_q,parent,pointer_EIBZ,nsest,indx_sest)

            exMat%l_real=fi%sym%invs
            CALL exchange_valence_hf(ikpt,fi,zMat, c_phase,nkpt_EIBZ,mpdata,jspin,hybdat,lapw,&
                                     eig_irr,results,pointer_EIBZ,n_q,wl_iks,xcpot,nococonv,nsest,indx_sest,&
                                     mpi,exMat)
            CALL exchange_vccv1(ikpt,fi%input,fi%atoms,fi%cell, fi%kpts, fi%sym, fi%noco,nococonv, fi%oned,&
                                mpdata,fi%hybinp,hybdat,jspin,lapw,nsymop,nsest,indx_sest,mpi,&
                                1.0,results,exMat)

            DEALLOCATE(indx_sest)

            !Start of workaround for increased functionality of fi%symmetrizeh (call it))

            nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

            CALL olap%alloc(fi%sym%invs,nbasfcn)
            CALL read_olap(olap, fi%kpts%nkpt*(jspin-1)+ikpt, nbasfcn)

            zMat%matsize2 = hybdat%nbands(ikpt) ! reduce "visible matsize" for the following computations

            CALL olap%multiply(zMat,trafo)

            CALL invtrafo%alloc(olap%l_real,hybdat%nbands(ikpt),nbasfcn)
            CALL trafo%TRANSPOSE(invtrafo)

            DO i = 1, hybdat%nbands(ikpt)
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

            call symmetrizeh(fi%atoms,fi%kpts%bkf(:,ikpt),jspin,lapw,fi%sym,hybdat%kveclo_eig,fi%cell,nsymop,psym,exMatLAPW)

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
!         CALL cdnvalJob%init(mpi,fi%input,fi%kpts,fi%noco,results,isp,fi%banddos=fi%banddos)
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
            END DO
         END DO
      END DO
      equalityCriterion = fi%input%zelec/(2.0/REAL(fi%input%jspins))
      gradient(numStates+1) = gradient(numStates+1) - equalityCriterion ! This should actually always be 0.0
      parameters(numStates+1) = lagrangeMultiplier

      WRITE(*,*) 'gradient(numStates+1): ', gradient(numStates+1)

      mixParam = 0.01 / MAXVAL(ABS(gradient(:numStates)))
!      mixParam = MIN(0.0002,mixParam)
      WRITE(*,*) 'mixParam: ', mixParam

      CALL bfgs_b2(numStates+1,gradient,lastGradient,minConstraints,maxConstraints,enabledConstraints,parameters,&
                   lastParameters,equalityLinCombi,equalityCriterion,maxHistoryLength,paramCorrections,&
                   gradientCorrections,iStep,mixParam,converged,convCrit)

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
      DEALLOCATE (parameters,gradient,equalityLinCombi)


   END DO ! WHILE (.NOT.converged)

   WRITE(2503,*) 'convIter: ', convIter

   WRITE(*,*) 'RDMFT: convergence loop end'

   hybdat%l_calhf = .FALSE.

   ! Calculate final overall density

   !I think we need most of cdngen at this place so I just use cdngen
   CALL outDen%resetPotDen()
   CALL cdngen(eig_id,mpi,fi%input,fi%banddos,fi%sliceplot,fi%vacuum,fi%kpts,fi%atoms,sphhar,stars,fi%sym,fi%gfinp,fi%hub1inp,&
               enpara,fi%cell,fi%noco,nococonv,vTot,results,fi%oned,fi%corespecinput,archiveType,xcpot,outDen, EnergyDen)

   ! Calculate RDMFT energy
   rdmftEnergy = 0.0
   DO ispin = 1, fi%input%jspins
      isp = MERGE(1,ispin,fi%noco%l_noco)
!      CALL cdnvalJob%init(mpi,fi%input,fi%kpts,fi%noco,results,isp,fi%banddos=fi%banddos)
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
      DO i = 1, fi%atoms%jri(iType)
         mt(i,iType) = outDen%mt(i,0,iType,1) + outDen%mt(i,0,iType,fi%input%jspins)
      END DO
   END DO
   IF (fi%input%jspins.EQ.1) mt=mt/2.0 !we just added the same value twice

   DO iType = 1, fi%atoms%ntype
      DO j = 1,fi%atoms%jri(iType)
         dpj(j) = mt(j,iType)/fi%atoms%rmsh(j,iType)
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


   WRITE(6,'(a,f20.10,a)') 'RDMFT energy: ', rdmftEnergy, ' Htr'

   WRITE(2505,*) '======================================='
   WRITE(2505,*) 'convIter: ', convIter
   WRITE(2505,'(a,f20.10,a)') 'RDMFT energy: ', rdmftEnergy, ' Htr'
   WRITE(2505,*) '======================================='
   WRITE(2505,*) '======================================='
   WRITE(2505,*) '======================================='

#endif
END SUBROUTINE rdmft

END MODULE m_rdmft
