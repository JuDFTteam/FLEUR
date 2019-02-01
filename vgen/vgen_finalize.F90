!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vgen_finalize
  USE m_juDFT
CONTAINS
  SUBROUTINE vgen_finalize(atoms,stars,vacuum,sym,noco,input,sphhar,vTot,vCoul,denRot)
    !     ***********************************************************
    !     FLAPW potential generator                           *
    !     ***********************************************************
    !     some rescaling is done here
    !     ***********************************************************
    !     in noco case vmatgen is called to generate 2x2 int-potential
    !     **********************************************************
    USE m_constants
    USE m_vmatgen
    USE m_types
    USE m_rotate_mt_den_tofrom_local
    IMPLICIT NONE
    TYPE(t_vacuum),INTENT(IN)       :: vacuum
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_stars),INTENT(IN)        :: stars
    TYPE(t_atoms),INTENT(IN)        :: atoms 
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_sphhar),INTENT(IN)       :: sphhar
    TYPE(t_potden),INTENT(INOUT)    :: vTot,vCoul,denRot
    !     ..
    !     .. Local Scalars ..
    INTEGER i,js,n

    !           ---> store v(l=0) component as r*v(l=0)/sqrt(4pi)
    
    DO js = 1,SIZE(vtot%mt,4)
       DO n = 1,atoms%ntype
          vTot%mt(:atoms%jri(n),0,n,js)  = atoms%rmsh(:atoms%jri(n),n)*vTot%mt(:atoms%jri(n),0,n,js)/sfp_const
       ENDDO
    ENDDO     ! js =1,input%jspins
    
    ! Rescale vTot%pw_w with number of stars
    IF (.NOT.noco%l_noco) THEN
       DO js=1,SIZE(vtot%pw_w,2)
          DO i=1,stars%ng3
             vTot%pw_w(i,js)=vtot%pw_w(i,js) / stars%nstr(i)
          END DO
       END DO
    ELSEIF(noco%l_noco) THEN
       CALL vmatgen(stars,atoms,vacuum,sym,input,denRot,vTot)
       IF (noco%l_mtnocoPot) CALL rotate_mt_den_from_local(atoms,sphhar,sym,denRot,vtot)
    ENDIF

    write (*,*) "Set vTot to zero in vgen_finalize()"
    vTot%pw_w = 0.0
    vTot%pw   = 0.0

    ! Rescale vCoul%pw_w with number of stars
    DO js = 1, SIZE(vCoul%pw_w,2)
       DO i = 1, stars%ng3
          vcoul%pw_w(i,js) = vcoul%pw_w(i,js) / stars%nstr(i)  !this normalization is needed for gw
       END DO
    END DO

    !Copy first vacuum into second vacuum if this was not calculated before 
    IF (vacuum%nvac==1) THEN
       vTot%vacz(:,2,:)  = vTot%vacz(:,1,:)
       IF (sym%invs) THEN
          vTot%vacxy(:,:,2,:)  = CMPLX(vTot%vacxy(:,:,1,:))
       ELSE
          vTot%vacxy(:,:,2,:)  = vTot%vacxy(:,:,1,:)
       ENDIF
    ENDIF
 
  END SUBROUTINE vgen_finalize
END MODULE m_vgen_finalize
