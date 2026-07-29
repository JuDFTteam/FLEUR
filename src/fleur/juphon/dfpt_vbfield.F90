!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vbfield
   USE m_juDFT
CONTAINS
  SUBROUTINE dfpt_vbfield(input,stars,noco,atoms,vTot,vTotIm)
    !This subroutine calculates the Zeeman field perturbation
    USE m_types
    USE m_constants
    USE m_rotMMPmat
    
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)::input
    TYPE(t_noco),INTENT(IN) ::noco
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_atoms),INTENT(IN)::atoms
    TYPE(t_potden),INTENT(INOUT)::vTot,vTotIm

    INTEGER :: iType
    REAL    :: bsign
    LOGICAL :: l_afm

    l_afm = .FALSE.

    IF (input%jspins.NE.2) CALL judft_error("B-fields can only be used in spin-polarized calculations")
    !IF (noco%l_noco) CALL judft_error("B-fields not implemented in noco case")

    vTot%pw(:,:)     = 0.0
    vTot%mt(:,:,:,:) = 0.0
    vTotIm%mt(:,:,:,:) = 0.0

    IF (l_afm) THEN
       !No IR contribution; alternate sign between sublattices in MT
       DO iType = 1, atoms%ntype
          bsign = MERGE(-1.0, 1.0, MOD(iType,2) == 1)
          vTot%mt(:atoms%jri(iType),0,iType,1) = -bsign * sfp_const/2.0
          vTot%mt(:atoms%jri(iType),0,iType,2) =  bsign * sfp_const/2.0
       END DO
    ELSE
       !Uniform field in IR and all MT
       vTot%pw(1,1) = -1/2.
       vTot%pw(1,2) = +1/2.
       DO iType = 1, atoms%ntype
          vTot%mt(:atoms%jri(iType),0,iType,1) = -sfp_const/2.0
          vTot%mt(:atoms%jri(iType),0,iType,2) = +sfp_const/2.0
       END DO
    END IF



  END SUBROUTINE dfpt_vbfield
END MODULE m_dfpt_vbfield