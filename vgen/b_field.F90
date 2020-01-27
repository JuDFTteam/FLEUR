!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_bfield
   USE m_juDFT
CONTAINS
  SUBROUTINE bfield(input,noco,atoms,field,vTot)
    !This subroutine adds a Zeeman-field to the potential
    !field%b_field is the field applied everywhere
    !field%b_field_mt is the field specific to the MT-sphere of a single atom type
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)::input
    TYPE(t_noco),INTENT(IN) ::noco
    TYPE(t_atoms),INTENT(IN)::atoms
    TYPE(t_field),INTENT(IN)::field
    TYPE(t_potden),INTENT(INOUT)::vtot

    INTEGER :: n

    
    IF (.NOT.field%l_b_field) RETURN !no B-field specified

    IF (input%jspins.NE.2) CALL judft_error("B-fields can only be used in spin-polarized calculations")
    IF (noco%l_noco) CALL judft_error("B-fields not implemented in noco case")
    
    !Interstitial
    vTot%pw(1,1)=vTot%pw(1,1)-field%b_field/2.
    vTot%pw(1,2)=vTot%pw(1,2)+field%b_field/2.

    !MT-spheres
    DO n=1,atoms%ntype
       vTot%mt(:,0,n,1)=vTot%mt(:,0,n,1)-(field%b_field+field%b_field_mt(n))/2.
       vTot%mt(:,0,n,2)=vTot%mt(:,0,n,2)+(field%b_field+field%b_field_mt(n))/2.
    ENDDO

    !Vacuum
    vTot%vacz(:,:,1)=vTot%vacz(:,:,1)-field%b_field/2.
    vTot%vacz(:,:,2)=vTot%vacz(:,:,2)+field%b_field/2.
  END SUBROUTINE bfield
END MODULE m_bfield
