!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_hamil
  use m_juDFT
  implicit none
CONTAINS

  subroutine hsmt_hamil(sym,atoms,ispin,cell,lapw,td,gk,vk,fj,gj,hamovlp)
    !Calculate hamiltonian matrix
#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_hsmt_ab
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    TYPE(t_hamovlp),INTENT(INOUT) :: hamovlp
    
    INTEGER:: n,nn,na,ab_offset
    COMPLEX,ALLOCATABLE:: ab(:,:)
    complex XX(2*atoms%lmaxd*(atoms%lmaxd+2)+2)
    ALLOCATE(ab(lapw%nv(1),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

   

    !Initialize The Hamiltonian with negatively shifted overlap
    hamovlp%h_c=-1.*td%e_shift*hamovlp%s_c
    
    DO n=1,atoms%ntype
       DO nn = 1,atoms%neq(n)
          na = SUM(atoms%neq(:n-1))+nn
          IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
             
             CALL hsmt_ab(sym,atoms,ispin,n,na,cell,lapw,gk,vk,fj,gj,ab,ab_offset)
             
             ab=matmul(ab,td%h_loc(:,:,n,ispin))
        
             CALL ZHERK("U","N",lapw%nv(ispin),2*ab_offset,cmplx(1.,0),ab(:,:),size(ab,1),cmplx(1.0,0.0),HamOvlp%h_c,size(HamOvlp%h_c,1))             
          ENDIF
       end DO
    end DO
    
  end subroutine hsmt_hamil

END MODULE m_hsmt_hamil
