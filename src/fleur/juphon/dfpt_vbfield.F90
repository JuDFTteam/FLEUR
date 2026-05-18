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
 
    real:: bohr

    IF (input%jspins.NE.2) CALL judft_error("B-fields can only be used in spin-polarized calculations")
    !IF (noco%l_noco) CALL judft_error("B-fields not implemented in noco case")
    
    !Interstitial
    !vTot%pw_w(:,1)=vTot%pw_w(:,1)-(1/2)*stars%ustep(:)
    !vTot%pw_w(:,2)=vTot%pw_w(:,2)+(1/2)*stars%ustep(:)
    vTot%pw(:,:)=0.0
    vTot%pw(1,1)=-1/2. 
    vTot%pw(1,2)=+1/2. 

    !MT-spheres
    DO iType = 1, atoms%ntype
       vTot%mt(:atoms%jri(iType),0,iType,1) = vTot%mt(:atoms%jri(iType),0,iType,1) - sfp_const/2.0
       vTot%mt(:atoms%jri(iType),0,iType,2) = vTot%mt(:atoms%jri(iType),0,iType,2) + sfp_const/2.0 

    END DO

    !Set MT imag to zero
    vTotIm%mt(:,:,:,:)=0.0

  END SUBROUTINE dfpt_vbfield
END MODULE m_dfpt_vbfield