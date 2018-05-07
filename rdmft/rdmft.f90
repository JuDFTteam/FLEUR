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

   CALL juDFT_error('rdmft not yet implemented!', calledby = 'rdmft')

END SUBROUTINE rdmft

END MODULE m_rdmft
