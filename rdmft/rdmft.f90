!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rdmft

CONTAINS

SUBROUTINE rdmft(eig_id,mpi,input,kpts,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
                 sphhar,sym,vTot,oneD,noco,results)

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_cdnval
   USE m_cdn_io

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
   TYPE(t_potden),        INTENT(IN)    :: vTot
   TYPE(t_oneD),          INTENT(IN)    :: oneD
   TYPE(t_noco),          INTENT(IN)    :: noco
   TYPE(t_results),       INTENT(INOUT) :: results

   INTEGER,               INTENT(IN)    :: eig_id

   TYPE(t_cdnvalJob)                    :: cdnvalJob
   TYPE(t_potden)                       :: singleStateDen
   TYPE(t_regionCharges)                :: regCharges
   TYPE(t_dos)                          :: dos
   TYPE(t_moments)                      :: moments
   INTEGER                              :: jspin, ikpt, iBand, jsp
   LOGICAL                              :: converged
   CHARACTER(LEN=20)                    :: filename

   CALL juDFT_error('rdmft not yet implemented!', calledby = 'rdmft')

   converged = .FALSE.

   ! Calculate all single state densities
   CALL regCharges%init(input,atoms)
   CALL dos%init(input,atoms,dimension,kpts,vacuum)
   CALL moments%init(input,atoms)
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
            CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN_const,CDN_INPUT_DEN_const,&
                              0,-1.0,0.0,.FALSE.,singleStateDen,TRIM(ADJUSTL(filename)))
         END DO
      END DO
   END DO

   DO WHILE (.NOT.converged)

      ! Calculate overall density with current occupation numbers (don't forget core electron density)

      ! Calculate Coulomb potential for overall density (+including external potential)

      ! For all states calculate integral over Coulomb potential times single state density

      ! For all states calculate Integral over other potential contributions times single state density

      ! Construct exchange matrix in the basis of eigenstates

      ! Optimize occupation numbers

      ! Check convergence of occupation numbers and set "converged" flag

   END DO ! WHILE (.NOT.converged)

   ! Calculate final overall density

   ! Calculate total energy

END SUBROUTINE rdmft

END MODULE m_rdmft
