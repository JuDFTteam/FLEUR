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
     COMPLEX,ALLOCATABLE :: pw(:,:)
     REAL,ALLOCATABLE    :: mt(:,:,:,:)
     REAL,ALLOCATABLE    :: vacz(:,:,:)
     COMPLEX,ALLOCATABLE :: vacxy(:,:,:,:)
     ! For density only (noco case)
     COMPLEX, ALLOCATABLE :: cdom(:)
     COMPLEX, ALLOCATABLE :: cdomvz(:,:)
     COMPLEX, ALLOCATABLE :: cdomvxy(:,:,:)
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
  END TYPE t_potden

CONTAINS
  SUBROUTINE init_potden_types(pd,stars,atoms,sphhar,vacuum,noco,oneD,jsp,nocoExtraDim,potden_type)
    USE m_judft
    USE m_types_misc
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT):: pd
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_sphhar),INTENT(IN):: sphhar
    TYPE(t_vacuum),INTENT(IN):: vacuum
    TYPE(t_noco),INTENT(IN)  :: noco
    TYPE(t_oneD),INTENT(IN)  :: oneD
    INTEGER,INTENT(IN)       :: jsp, potden_type
    LOGICAL,INTENT(IN)       :: nocoExtraDim

    CALL init_potden_simple(pd,stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,&
                            atoms%n_u,noco%l_noco,jsp,nocoExtraDim,potden_type,&
                            vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
  END SUBROUTINE init_potden_types

  SUBROUTINE init_potden_simple(pd,ng3,jmtd,nlhd,ntype,n_u,l_noco,jsp,nocoExtraDim,potden_type,nmzd,nmzxyd,n2d)
    USE m_constants
    USE m_judft
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT) :: pd
    INTEGER,INTENT(IN)          :: ng3,jmtd,nlhd,ntype,n_u,jsp,potden_type
    LOGICAL,INTENT(IN)          :: l_noco,nocoExtraDim
    INTEGER,INTENT(IN)          :: nmzd,nmzxyd,n2d

    INTEGER:: err(4)

    err=0
    pd%iter=0
    pd%potdenType=potden_type
    IF(ALLOCATED(pd%pw)) DEALLOCATE (pd%pw)
    IF(ALLOCATED(pd%mt)) DEALLOCATE (pd%mt)
    IF(ALLOCATED(pd%vacz)) DEALLOCATE (pd%vacz)
    IF(ALLOCATED(pd%vacxy)) DEALLOCATE (pd%vacxy)
    IF(ALLOCATED(pd%cdom)) DEALLOCATE (pd%cdom)
    IF(ALLOCATED(pd%cdomvz)) DEALLOCATE (pd%cdomvz)
    IF(ALLOCATED(pd%cdomvxy)) DEALLOCATE (pd%cdomvxy)
    IF(ALLOCATED(pd%mmpMat)) DEALLOCATE (pd%mmpMat)
    ALLOCATE (pd%pw(ng3,jsp),stat=err(1))
    ALLOCATE (pd%mt(jmtd,0:nlhd,ntype,jsp),stat=err(2))
    ALLOCATE (pd%vacz(nmzd,2,MERGE(4,jsp,nocoExtraDim)),stat=err(3))
    ALLOCATE (pd%vacxy(nmzxyd,n2d-1,2,jsp),stat=err(4))
    IF (l_noco) THEN
       ALLOCATE (pd%cdom(ng3))
       ALLOCATE (pd%cdomvz(nmzd,2))
       ALLOCATE (pd%cdomvxy(nmzxyd,n2d-1,2))
    ELSE
       ALLOCATE (pd%cdom(1))
       ALLOCATE (pd%cdomvz(1,1),pd%cdomvxy(1,1,1))
    END IF

    ALLOCATE (pd%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,n_u),jsp))

    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating potential or density")
    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    pd%vacz=0.0
    pd%vacxy=CMPLX(0.0,0.0)
    pd%cdom = CMPLX(0.0,0.0)
    pd%cdomvz = CMPLX(0.0,0.0)
    pd%cdomvxy = CMPLX(0.0,0.0)
    pd%mmpMat = CMPLX(0.0,0.0)
  END SUBROUTINE init_potden_simple

  SUBROUTINE resetPotDen(pd)

    IMPLICIT NONE

    CLASS(t_potden),INTENT(INOUT) :: pd

    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    pd%vacz=0.0
    pd%vacxy=CMPLX(0.0,0.0)
    pd%cdom = CMPLX(0.0,0.0)
    pd%cdomvz = CMPLX(0.0,0.0)
    pd%cdomvxy = CMPLX(0.0,0.0)
    pd%mmpMat = CMPLX(0.0,0.0)
  END SUBROUTINE resetPotDen
END MODULE m_types_potden
