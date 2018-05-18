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
   USE m_cdnval

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

   INTEGER                              :: jspin, ikpt, iBand
   LOGICAL                              :: converged

   CALL juDFT_error('rdmft not yet implemented!', calledby = 'rdmft')

   converged = .FALSE.

   ! Calculate all single state densities
   DO jspin = 1, input%jspins
      DO ikpt = 1, kpts%nkpt
         DO iBand = 1, results%neig(ikpt,jspin)
            ! Construct cdnvalJob object for this state

            ! Call cdnval to construct density
!            CALL cdnval(eig_id,mpi,kpts,jspin,noco,input,banddos,cell,atoms,enpara,stars,vacuum,dimension,&
!                        sphhar,sym,vTot,oneD,cdnvalJob,outDen,regCharges,dos,results,moments,orbcomp)

            ! Store the density on disc (These are probably way too many densities to keep them in memory)
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
