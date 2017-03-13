!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_hamil
  USE m_juDFT
  IMPLICIT NONE
CONTAINS

  SUBROUTINE hsmt_hamil(sym,atoms,ispin,cell,lapw,td,gk,vk,fj,gj,hamovlp)
    !Calculate hamiltonian matrix
#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_hsmt_ab
#ifdef _OPENACC
    USE cublas 
    USE cudafor
    USE openacc
#endif
#ifdef __PGI
!    USE nvtx
#endif
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
    
    INTEGER:: n,nn,na,ab_offset(NO_STREAMS)
    COMPLEX,ALLOCATABLE:: ab(:,:)
  
#ifdef _OPENACC
    TYPE(cublasHandle) :: cublas_handle
    INTEGER            :: istat,n1,n2,n3,n4,istream
    INTEGER,PARAMETER  :: NO_STREAM=2
    COMPLEX, device, ALLOCATABLE :: hloc_tmp(:,:),hc_tmp(:,:) 
    istat = cublasCreate(cublas_handle)
    ALLOCATE(hloc_tmp(SIZE(td%h_loc,1),SIZE(td%h_loc,2)))
    ALLOCATE(hc_tmp(SIZE(hamovlp%h_c,1),SIZE(hamovlp%h_c,2)))
#else
    INTEGER,PARAMETER  :: NO_STREAMS=1
    !CALL nvtxStartRange("hsmt_hamil",1)
    n_ab1=lapw%nv(1)
    n_ab2=2*atoms%lmaxd*(atoms%lmaxd+2)+2

    ALLOCATE(ab(n_ab1,n_ab2,NO_STREAMS))
    
    !Initialize The Hamiltonian with negatively shifted overlap
    hamovlp%h_c=-1.*td%e_shift*hamovlp%s_c
#ifdef _OPENACC
    hc_tmp=hamovlp%h_c
#endif
    !$acc data create(ab)
    DO n=1,atoms%ntype
#ifdef _OPENACC
       hloc_tmp=td%h_loc(:,:,n,ispin)! stream missing -> wait
       !$acc wait
       DO nn = 1,atoms%neq(n)
          na = SUM(atoms%neq(:n-1))+nn
          IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
             istream=mod(na,NO_STREAMS)
             ! istream in handle!
             CALL hsmt_ab(sym,atoms,ispin,n,na,cell,lapw,gk,vk,fj,gj,ab(:,:,istream+1),ab_offset(istream+1))
             !$acc update device (ab(:,:,istream+1)) async(istream)
#ifdef _OPENACC             
             !$acc host_data use_device(ab)
             istat = cublaszgemm_v2(cublas_handle(istream+1),CUBLAS_OP_N,CUBLAS_OP_N,n_ab1,n_ab2,n_ab2,CMPLX(1.0,0.0),ab(:,:,istream+1),n_ab1,hloc_tmp,SIZE(td%h_loc,1),CMPLX(0.,0.),ab,n_ab1)          !$acc wait(1) async(0) !  
             istat = cublaszherk_v2(cublas_handle(0),CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_N,lapw%nv(ispin),2*ab_offset(istream+1),1.,ab(:,:,istream+1),n_ab1,1.0,hc_tmp,SIZE(hamovlp%h_c,1)) 
             !$acc end host_data
#else             
             CALL zgemm("N","N",SIZE(ab,1),SIZE(ab,2),SIZE(ab,2),CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,ispin),SIZE(td%h_loc,1),CMPLX(0.,0.),ab,SIZE(ab,1))
             CALL ZHERK("U","N",lapw%nv(ispin),2*ab_offset(1),CMPLX(1.,0),ab,SIZE(ab,1),CMPLX(1.0,0.0),HamOvlp%h_c,SIZE(HamOvlp%h_c,1))             
#endif
             
          ENDIF
       END DO
    END DO
    !$acc end data
    !CALL nvtxendrange
#ifdef _OPENACC
    hamovlp%h_c=hc_tmp
#endif
  END SUBROUTINE hsmt_hamil
  
END MODULE m_hsmt_hamil
