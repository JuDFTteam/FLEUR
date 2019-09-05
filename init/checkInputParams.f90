!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_checkInputParams

CONTAINS

SUBROUTINE checkInputParams(mpi,input,dimension,atoms,noco,xcpot,oneD)

   USE m_juDFT
   USE m_types

   TYPE(t_mpi),           INTENT(IN)    :: mpi
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_dimension),     INTENT(IN)    :: dimension
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_noco),          INTENT(IN)    :: noco
   CLASS(t_xcpot),        INTENT(IN)    :: xcpot
   TYPE(t_oneD),          INTENT(IN)    :: oneD

   IF(mpi%irank.NE.0) RETURN

   IF(noco%l_noco.AND.noco%l_soc.AND.ANY(atoms%nlo(:).GT.0)) THEN
      CALL juDFT_warn('In l_noco+l_soc SOC contributions to LOs are not taken into account!',&
                      calledby = 'checkInputParams',hint='If you know what you do deactivate this stop.')
   END IF

   IF (((xcpot%is_hybrid().OR.input%l_rdmft)).AND.((input%film.OR.oneD%odi%d1))) THEN
      CALL juDFT_error("2D film and 1D calculations not implemented for HF/EXX/PBE0/HSE", &
                       calledby ="fleur", hint="Use a supercell or a different functional")
   END IF

END SUBROUTINE checkInputParams


END MODULE m_checkInputParams
