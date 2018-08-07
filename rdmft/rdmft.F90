!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rdmft

CONTAINS

SUBROUTINE rdmft(eig_id,mpi,input,kpts,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                 sphhar,sym,field,vTot,oneD,noco,results)

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_cdnval
   USE m_cdn_io
   USE m_cdncore
   USE m_qfix
   USE m_vgen_coulomb
   USE m_convol
   USE m_intnv
#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif

   IMPLICIT NONE

   TYPE(t_mpi),           INTENT(IN)    :: mpi
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   TYPE(t_banddos),       INTENT(IN)    :: banddos
   TYPE(t_cell),          INTENT(IN)    :: cell
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_enpara),        INTENT(IN)    :: enpara
   TYPE(t_stars),         INTENT(IN)    :: stars
   TYPE(t_vacuum),        INTENT(IN)    :: vacuum
   TYPE(t_dimension),     INTENT(IN)    :: dimension
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_field),         INTENT(INOUT) :: field
   TYPE(t_potden),        INTENT(INOUT) :: vTot
   TYPE(t_oneD),          INTENT(IN)    :: oneD
   TYPE(t_noco),          INTENT(IN)    :: noco
   TYPE(t_results),       INTENT(INOUT) :: results

   INTEGER,               INTENT(IN)    :: eig_id

   TYPE(t_cdnvalJob)                    :: cdnvalJob
   TYPE(t_potden)                       :: singleStateDen, overallDen, overallVCoul
   TYPE(t_regionCharges)                :: regCharges
   TYPE(t_dos)                          :: dos
   TYPE(t_moments)                      :: moments
   INTEGER                              :: jspin, ikpt, iBand, jsp, jspmax
   REAL                                 :: fix, potDenInt, fermiEnergyTemp
   LOGICAL                              :: converged, l_qfix
   CHARACTER(LEN=20)                    :: filename

   REAL, ALLOCATABLE                    :: overallVCoulSSDen(:,:,:)
   REAL, ALLOCATABLE                    :: vTotSSDen(:,:,:)

   CALL juDFT_error('rdmft not yet implemented!', calledby = 'rdmft')

   ! General initializations
   ALLOCATE(overallVCoulSSDen(MAXVAL(results%neig(1:kpts%nkpt,1:input%jspins)),kpts%nkpt,input%jspins))
   ALLOCATE(vTotSSDen(MAXVAL(results%neig(1:kpts%nkpt,1:input%jspins)),kpts%nkpt,input%jspins))

   CALL regCharges%init(input,atoms)
   CALL dos%init(input,atoms,dimension,kpts,vacuum)
   CALL moments%init(input,atoms)
   CALL overallDen%init(stars,atoms,sphhar,vacuum,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)
   CALL overallVCoul%init(stars,atoms,sphhar,vacuum,input%jspins,noco%l_noco,POTDEN_TYPE_POTCOUL)
   ALLOCATE(overallVCoul%pw_w(SIZE(overallDen%pw,1),1))
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
   DO jspin = 1, input%jspins
      jsp = MERGE(1,jspin,noco%l_noco)
      DO ikpt = 1, kpts%nkpt
         DO iBand = 1, results%neig(ikpt,jsp)
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
            CALL singleStateDen%init(stars,atoms,sphhar,vacuum,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)
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

   ! Construct exchange matrix in the basis of eigenstates
   ! TODO!!!!!

   converged = .FALSE.
   DO WHILE (.NOT.converged)

      ! Calculate overall density with current occupation numbers (don't forget core electron density)
      jspmax = input%jspins
      IF (noco%l_mperp) jspmax = 1
      DO jspin = 1,jspmax
         CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin,banddos=banddos)
         CALL cdnval(eig_id,mpi,kpts,jsp,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                     sphhar,sym,vTot,oneD,cdnvalJob,overallDen,regCharges,dos,results,moments)
      END DO
      CALL cdncore(mpi,dimension,oneD,input,vacuum,noco,sym,&
                   stars,cell,sphhar,atoms,vTot,overallDen,moments,results)
      IF (mpi%irank.EQ.0) THEN
         CALL qfix(stars,atoms,sym,vacuum,sphhar,input,cell,oneD,overallDen,noco%l_noco,.TRUE.,.true.,fix)
      END IF
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,overallDen)
#endif

      ! Calculate Coulomb potential for overall density (+including external potential)
      CALL overallDen%sum_both_spin()!workden)
      CALL overallVCoul%resetPotDen()
      CALL vgen_coulomb( 1, mpi, DIMENSION, oneD, input, field, vacuum, sym, stars, cell, &
           sphhar, atoms, overallDen, overallVCoul )
      CALL convol(stars, overallVCoul%pw_w(:,1), overallVCoul%pw(:,1), stars%ufft)   ! Is there a problem with a second spin?!
#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,overallVCoul)
#endif

      overallVCoulSSDen = 0.0
      DO jspin = 1, input%jspins
         jsp = MERGE(1,jspin,noco%l_noco)
         DO ikpt = 1, kpts%nkpt
            DO iBand = 1, results%neig(ikpt,jsp)
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

      ! Optimize occupation numbers

      ! Check convergence of occupation numbers and set "converged" flag

   END DO ! WHILE (.NOT.converged)

   ! Calculate final overall density

   ! Calculate total energy

END SUBROUTINE rdmft

END MODULE m_rdmft
