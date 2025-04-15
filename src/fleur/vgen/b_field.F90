!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_bfield
   USE m_juDFT
CONTAINS
  SUBROUTINE bfield(input,stars,noco,atoms,field,vTot,nococonv)
    !This subroutine adds a Zeeman-field to the potential
    !field%b_field is the field applied everywhere
    !field%b_field_mt is the field specific to the MT-sphere of a single atom type
    USE m_types
    USE m_constants
    USE m_rotMMPmat
    
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)::input
    TYPE(t_noco),INTENT(IN) ::noco
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_atoms),INTENT(IN)::atoms
    TYPE(t_field),INTENT(IN)::field
    TYPE(t_potden),INTENT(INOUT)::vtot
    TYPE(t_nococonv):: nococonv
    INTEGER :: n
    COMPLEX :: bComplexTemp
    REAL :: bExternal(4)

    
    IF (.NOT.field%l_b_field) RETURN !no B-field specified

    IF (input%jspins.NE.2) CALL judft_error("B-fields can only be used in spin-polarized calculations")
    !IF (noco%l_noco) CALL judft_error("B-fields not implemented in noco case")
    
    !Interstitial
    vTot%pw_w(:,1)=vTot%pw_w(:,1)-((field%b_field/2.0)*stars%ustep(:)*stars%nstr(:))
    vTot%pw_w(:,2)=vTot%pw_w(:,2)+((field%b_field/2.0)*stars%ustep(:)*stars%nstr(:))

    !MT-spheres
    DO n=1,atoms%ntype
       bExternal(1) = -field%b_field / 2
       bExternal(2) = field%b_field / 2

       IF (noco%l_noco) THEN
          CALL nococonv%rotdenmat(nococonv%alph(n), nococonv%beta(n), bExternal(1), bExternal(2), bComplexTemp, .false.)
          bExternal(3) = REAL(bComplexTemp)
          bExternal(4) = AIMAG(bComplexTemp)
       END IF

       vTot%mt(:atoms%jri(n),0,n,1) = vTot%mt(:atoms%jri(n),0,n,1) + sfp_const * (bExternal(1) * 2.0 - field%b_field_mt(n)) / 2.0
       vTot%mt(:atoms%jri(n),0,n,2) = vTot%mt(:atoms%jri(n),0,n,2) + sfp_const * (bExternal(2) * 2.0 + field%b_field_mt(n)) / 2.0

       IF (noco%l_noco) THEN
          vTot%mt(:atoms%jri(n),0,n,3) = vTot%mt(:atoms%jri(n),0,n,3) + sfp_const * (bExternal(3))
          vTot%mt(:atoms%jri(n),0,n,4) = vTot%mt(:atoms%jri(n),0,n,4) + sfp_const * (bExternal(4))
       END IF
    ENDDO

    !Vacuum
    IF (input%film) THEN
       vTot%vac(:,1,:,1)=vTot%vac(:,1,:,1)-field%b_field/2.
       vTot%vac(:,1,:,2)=vTot%vac(:,1,:,2)+field%b_field/2.
    END IF
  END SUBROUTINE bfield
END MODULE m_bfield
