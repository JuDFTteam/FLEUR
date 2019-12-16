!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_force

IMPLICIT NONE

PRIVATE

   TYPE t_force
      COMPLEX, ALLOCATABLE :: f_a12(:,:)
      COMPLEX, ALLOCATABLE :: f_a21(:,:)
      COMPLEX, ALLOCATABLE :: f_b4(:,:)
      COMPLEX, ALLOCATABLE :: f_b8(:,:)

      COMPLEX, ALLOCATABLE :: e1cof(:,:,:)
      COMPLEX, ALLOCATABLE :: e2cof(:,:,:)
      COMPLEX, ALLOCATABLE :: aveccof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: bveccof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: cveccof(:,:,:,:,:)

      COMPLEX, ALLOCATABLE :: acoflo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: bcoflo(:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init1 => force_init1
         PROCEDURE,PASS :: init2 => force_init2
         PROCEDURE      :: addContribsA21A12
   END TYPE t_force

PUBLIC t_force

CONTAINS

SUBROUTINE force_init1(thisForce,input,atoms)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_force),     INTENT(INOUT) :: thisForce
   TYPE(t_input),      INTENT(IN)    :: input
   TYPE(t_atoms),      INTENT(IN)    :: atoms

   IF (input%l_f) THEN
      ALLOCATE (thisForce%f_a12(3,atoms%ntype))
      ALLOCATE (thisForce%f_a21(3,atoms%ntype))
      ALLOCATE (thisForce%f_b4(3,atoms%ntype))
      ALLOCATE (thisForce%f_b8(3,atoms%ntype))
   ELSE
      ALLOCATE (thisForce%f_a12(1,1))
      ALLOCATE (thisForce%f_a21(1,1))
      ALLOCATE (thisForce%f_b4(1,1))
      ALLOCATE (thisForce%f_b8(1,1))
   END IF

   thisForce%f_a12 = CMPLX(0.0,0.0)
   thisForce%f_a21 = CMPLX(0.0,0.0)
   thisForce%f_b4 = CMPLX(0.0,0.0)
   thisForce%f_b8 = CMPLX(0.0,0.0)

END SUBROUTINE force_init1

SUBROUTINE force_init2(thisForce,noccbd,input,atoms)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_force),     INTENT(INOUT) :: thisForce
   TYPE(t_input),      INTENT(IN)    :: input
   TYPE(t_atoms),      INTENT(IN)    :: atoms
   INTEGER,            INTENT(IN)    :: noccbd

   IF (ALLOCATED(thisForce%e1cof)) DEALLOCATE(thisForce%e1cof)
   IF (ALLOCATED(thisForce%e2cof)) DEALLOCATE(thisForce%e2cof)
   IF (ALLOCATED(thisForce%acoflo)) DEALLOCATE(thisForce%acoflo)
   IF (ALLOCATED(thisForce%bcoflo)) DEALLOCATE(thisForce%bcoflo)
   IF (ALLOCATED(thisForce%aveccof)) DEALLOCATE(thisForce%aveccof)
   IF (ALLOCATED(thisForce%bveccof)) DEALLOCATE(thisForce%bveccof)
   IF (ALLOCATED(thisForce%cveccof)) DEALLOCATE(thisForce%cveccof)

   IF (input%l_f) THEN
      ALLOCATE (thisForce%e1cof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%e2cof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%acoflo(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))
      ALLOCATE (thisForce%bcoflo(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))
      ALLOCATE (thisForce%aveccof(3,noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%bveccof(3,noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%cveccof(3,-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))
   ELSE
      ALLOCATE (thisForce%e1cof(1,1,1))
      ALLOCATE (thisForce%e2cof(1,1,1))
      ALLOCATE (thisForce%acoflo(1,1,1,1))
      ALLOCATE (thisForce%bcoflo(1,1,1,1))
      ALLOCATE (thisForce%aveccof(1,1,1,1))
      ALLOCATE (thisForce%bveccof(1,1,1,1))
      ALLOCATE (thisForce%cveccof(1,1,1,1,1))
   END IF

   thisForce%e1cof = CMPLX(0.0,0.0)
   thisForce%e2cof = CMPLX(0.0,0.0)
   thisForce%acoflo = CMPLX(0.0,0.0)
   thisForce%bcoflo = CMPLX(0.0,0.0)
   thisForce%aveccof = CMPLX(0.0,0.0)
   thisForce%bveccof = CMPLX(0.0,0.0)
   thisForce%cveccof = CMPLX(0.0,0.0)

END SUBROUTINE force_init2

SUBROUTINE addContribsA21A12(thisForce,input,atoms,sym,cell,oneD,enpara,&
                             usdus,eigVecCoeffs,noccbd,ispin,eig,we,results)

   USE m_types_setup
   USE m_types_usdus
   USE m_types_enpara
   USE m_types_cdnval, ONLY: t_eigVecCoeffs
   USE m_types_misc
   USE m_forcea12
   USE m_forcea21

   IMPLICIT NONE

   CLASS(t_force),       INTENT(INOUT) :: thisForce
   TYPE(t_input),        INTENT(IN)    :: input
   TYPE(t_atoms),        INTENT(IN)    :: atoms
   
   TYPE(t_sym),          INTENT(IN)    :: sym
   TYPE(t_cell),         INTENT(IN)    :: cell
   TYPE(t_oneD),         INTENT(IN)    :: oneD
   TYPE(t_enpara),       INTENT(IN)    :: enpara
   TYPE(t_usdus),        INTENT(IN)    :: usdus
   TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
   TYPE(t_results),      INTENT(INOUT) :: results

   INTEGER,              INTENT(IN)    :: noccbd
   INTEGER,              INTENT(IN)    :: ispin

   REAL,                 INTENT(IN)    :: eig(noccbd)
   REAL,                 INTENT(IN)    :: we(noccbd)

   IF (.NOT.input%l_useapw) THEN
      CALL force_a12(atoms,noccbd,sym,cell,oneD,&
                     we,ispin,noccbd,usdus,eigVecCoeffs,thisForce%acoflo,thisForce%bcoflo,&
                     thisForce%e1cof,thisForce%e2cof,thisForce%f_a12,results)
   END IF
   CALL force_a21(input,atoms,sym,oneD,cell,we,ispin,&
                  enpara%el0(0:,:,ispin),noccbd,eig,usdus,eigVecCoeffs,&
                  thisForce%aveccof,thisForce%bveccof,thisForce%cveccof,&
                  thisForce%f_a21,thisForce%f_b4,results)

END SUBROUTINE addContribsA21A12

END MODULE m_types_force
