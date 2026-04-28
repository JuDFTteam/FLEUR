!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_addContribsA21A12

    IMPLICIT NONE
  
    PRIVATE
    PUBLIC :: addContribsA21A12
  contains

SUBROUTINE addContribsA21A12(thisForce,input,atoms,sym,cell ,enpara,&
    usdus,tlmplm,vtot,abc,noccbd,ispin,eig,we,results,jsp_start,jspin,nbasfcn,zMat,lapw,sphhar,k1,k2,k3,bkpt,itype)
 use m_types_force
 USE m_types_setup
 USE m_types_lapw
 USE m_types_mat
 USE m_types_sphhar
 USE m_types_usdus
 USE m_types_tlmplm
 USE m_types_enpara
 USE m_types_abc
 USE m_types_misc
 USE m_types_potden
 USE m_forcea12
 USE m_forcea21
 USE m_force_a12_lv2

 IMPLICIT NONE

 type(t_force),       INTENT(INOUT) :: thisForce
 TYPE(t_input),        INTENT(IN)    :: input
 TYPE(t_atoms),        INTENT(IN)    :: atoms
 TYPE(t_sym),          INTENT(IN)    :: sym
 TYPE(t_cell),         INTENT(IN)    :: cell
  
 TYPE(t_enpara),       INTENT(IN)    :: enpara
 TYPE(t_usdus),        INTENT(IN)    :: usdus
 TYPE(t_tlmplm),       INTENT(IN)    :: tlmplm
 TYPE(t_potden),       INTENT(IN)    :: vtot
 TYPE(t_abc), INTENT(IN)             :: abc
 TYPE(t_results),      INTENT(INOUT) :: results
 TYPE(t_lapw), INTENT(IN)            :: lapw
 TYPE(t_mat), INTENT(IN)             :: zMat
 TYPE(t_sphhar), INTENT(IN)          :: sphhar
 INTEGER,INTENT(IN)                  :: itype

 INTEGER,              INTENT(IN)    :: noccbd,k1(:,:),k2(:,:),k3(:,:)
 INTEGER,              INTENT(IN)    :: ispin,jsp_start,jspin,nbasfcn

 REAL,                 INTENT(IN)    :: eig(noccbd),bkpt(3)
 REAL,                 INTENT(IN)    :: we(noccbd)


 IF (.NOT.input%l_useapw) THEN

    IF (input%f_level.LT.2) THEN
       CALL force_a12(atoms,noccbd,sym,cell ,&
            we,ispin,noccbd,usdus,abc,thisForce%acoflo,thisForce%bcoflo,&
            thisForce%e1cof,thisForce%e2cof,thisForce%f_a12,results,itype)
    ELSE ! Klueppelberg (force level 2)
       IF (ispin.eq.jsp_start) THEN ! since we use IS rep, this part needs to be calculated only once
          CALL force_a12_lv2(jspin,input%jspins,noccbd,input%neig,atoms%ntype,atoms%ntype,atoms%nat,nbasfcn, &! or ispin?
               sym%nop,lapw%dim_nvd(),atoms%lmaxd,cell%omtil,lapw%nv,atoms%neq,k1,k2,k3, &
               sym%invarind,sym%invarop,sym%invtab,sym%mrot,sym%ngopr,cell%amat,cell%bmat,eig,atoms%rmt,atoms%taual,we,bkpt,zMat,thisForce%f_a12,results%force,itype )
       END IF
    END IF
 END IF
 CALL force_a21(input,atoms,sym ,cell,we,ispin,&
      enpara%el0(0:,:,ispin),noccbd,eig,usdus,tlmplm,vtot,abc,&
      thisForce%aveccof,thisForce%bveccof,thisForce%cveccof,&
      thisForce%f_a21,thisForce%f_b4,results,itype)

END SUBROUTINE addContribsA21A12
end module