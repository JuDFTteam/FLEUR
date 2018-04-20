!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_potden

  !> Data type for the density or the potential
   TYPE t_potden
     INTEGER             :: iter  
     INTEGER             :: potdenType
     COMPLEX,ALLOCATABLE :: pw(:,:),pw_w(:,:)
     REAL,ALLOCATABLE    :: mt(:,:,:,:)
     REAL,ALLOCATABLE    :: vacz(:,:,:)
     COMPLEX,ALLOCATABLE :: vacxy(:,:,:,:)
     !For angles of density/potential in noco case
     REAL,ALLOCATABLE  :: theta_pw(:)
     REAL,ALLOCATABLE  :: phi_pw(:)
     REAL,ALLOCATABLE  :: theta_vacz(:,:)
     REAL,ALLOCATABLE  :: phi_vacz(:,:)
     REAL,ALLOCATABLE  :: theta_vacxy(:,:,:)
     REAL,ALLOCATABLE  :: phi_vacxy(:,:,:)
     

     ! For density matrix and associated potential matrix
     COMPLEX, ALLOCATABLE :: mmpMat(:,:,:,:)

     !this type contains two init routines that should be used to allocate
     !memory. You can either specify the datatypes or give the dimensions as integers
     !See implementation below!
   CONTAINS
     PROCEDURE :: init_potden_types
     PROCEDURE :: init_potden_simple
     PROCEDURE :: resetpotden
     GENERIC   :: init=>init_potden_types,init_potden_simple
     PROCEDURE :: copy_both_spin
     PROCEDURE :: sum_both_spin
  END TYPE t_potden

CONTAINS

  SUBROUTINE sum_both_spin(this,that)
    IMPLICIT NONE
    CLASS(t_potden),INTENT(INOUT)   :: this
    TYPE(t_potden),INTENT(INOUT),OPTIONAL :: that

    IF (PRESENT(that)) THEN
       IF (SIZE(this%pw,2)>1) THEN
          that%mt(:,0:,:,1)=this%mt(:,0:,:,1)+this%mt(:,0:,:,2)
          that%pw(:,1)=this%pw(:,1)+this%pw(:,2)
          that%vacz(:,:,1)=this%vacz(:,:,1)+this%vacz(:,:,2)
          that%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)+this%vacxy(:,:,:,2)
          IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,1)=this%pw_w(:,1)+this%pw_w(:,2)
       ELSE
          that%mt(:,0:,:,1)=this%mt(:,0:,:,1)
          that%pw(:,1)=this%pw(:,1)
          that%vacz(:,:,1)=this%vacz(:,:,1)
          that%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)
          IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,1)=this%pw_w(:,1)
       ENDIF
    ELSE
       IF (SIZE(this%pw,2)>1) THEN
          this%mt(:,0:,:,1)=this%mt(:,0:,:,1)+this%mt(:,0:,:,2)
          this%pw(:,1)=this%pw(:,1)+this%pw(:,2)
          this%vacz(:,:,1)=this%vacz(:,:,1)+this%vacz(:,:,2)
          this%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)+this%vacxy(:,:,:,2)
          IF (ALLOCATED(this%pw_w)) this%pw_w(:,1)=this%pw_w(:,1)+this%pw_w(:,2)
       ENDIF
    END IF
  END SUBROUTINE sum_both_spin
    
  SUBROUTINE copy_both_spin(this,that)
    IMPLICIT NONE
    CLASS(t_potden),INTENT(IN)   :: this
    TYPE(t_potden),INTENT(INOUT) :: that

    that%mt(:,0:,:,1)=this%mt(:,0:,:,1)
    that%pw(:,1)=this%pw(:,1)
    that%vacz(:,:,1)=this%vacz(:,:,1)
    that%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)
    IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,1)=this%pw_w(:,1)
    
    IF (SIZE(that%mt,4)==2) THEN
       that%mt(:,0:,:,2)=this%mt(:,0:,:,1)
       that%pw(:,2)=this%pw(:,1)
       that%vacz(:,:,2)=this%vacz(:,:,1)
       that%vacxy(:,:,:,2)=this%vacxy(:,:,:,1)
       IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,2)=this%pw_w(:,1)
    END IF
  END SUBROUTINE copy_both_spin
  
  SUBROUTINE init_potden_types(pd,stars,atoms,sphhar,vacuum,jspins,nocoExtraDim,potden_type)
    USE m_judft
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT):: pd 
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_sphhar),INTENT(IN):: sphhar
    TYPE(t_vacuum),INTENT(IN):: vacuum
    INTEGER,INTENT(IN)       :: jspins, potden_type
    LOGICAL,INTENT(IN)       :: nocoExtraDim

    CALL init_potden_simple(pd,stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,&
                            atoms%n_u,jspins,nocoExtraDim,potden_type,&
                            vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
  END SUBROUTINE init_potden_types

  SUBROUTINE init_potden_simple(pd,ng3,jmtd,nlhd,ntype,n_u,jspins,nocoExtraDim,potden_type,nmzd,nmzxyd,n2d)
    USE m_constants
    USE m_judft
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT) :: pd
    INTEGER,INTENT(IN)          :: ng3,jmtd,nlhd,ntype,n_u,jspins,potden_type
    LOGICAL,INTENT(IN)          :: nocoExtraDim
    INTEGER,INTENT(IN)          :: nmzd,nmzxyd,n2d

    INTEGER:: err(4)

    err=0
    pd%iter=0
    pd%potdenType=potden_type
    IF(ALLOCATED(pd%pw)) DEALLOCATE (pd%pw)
    IF(ALLOCATED(pd%mt)) DEALLOCATE (pd%mt)
    IF(ALLOCATED(pd%vacz)) DEALLOCATE (pd%vacz)
    IF(ALLOCATED(pd%vacxy)) DEALLOCATE (pd%vacxy)
    IF(ALLOCATED(pd%mmpMat)) DEALLOCATE (pd%mmpMat)
    ALLOCATE (pd%pw(ng3,MERGE(3,jspins,nocoExtraDim)),stat=err(1))
    ALLOCATE (pd%mt(jmtd,0:nlhd,ntype,jspins),stat=err(2))
    ALLOCATE (pd%vacz(nmzd,2,MERGE(4,jspins,nocoExtraDim)),stat=err(3))
    ALLOCATE (pd%vacxy(nmzxyd,n2d-1,2,MERGE(3,jspins,nocoExtraDim)),stat=err(4))

    ALLOCATE (pd%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,n_u),jspins))

    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating potential or density")
    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    pd%vacz=0.0
    pd%vacxy=CMPLX(0.0,0.0)
    pd%mmpMat = CMPLX(0.0,0.0)
  END SUBROUTINE init_potden_simple

  SUBROUTINE resetPotDen(pd)

    IMPLICIT NONE

    CLASS(t_potden),INTENT(INOUT) :: pd

    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    pd%vacz=0.0
    pd%vacxy=CMPLX(0.0,0.0)
    pd%mmpMat = CMPLX(0.0,0.0)
  END SUBROUTINE resetPotDen
END MODULE m_types_potden
