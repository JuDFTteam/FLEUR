!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_nonsph
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC hsmt_nonsph

CONTAINS
  SUBROUTINE hsmt_nonsph(n,fmpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,td,fjgj,hmat,set0)
    USE m_hsmt_fjgj
    USE m_types
    USE m_hsmt_ab
#ifdef _OPENACC
    USE cublas
#define CPP_zgemm cublaszgemm
#define CPP_zherk cublaszherk
#define CPP_data_c data_c
#else
#define CPP_zgemm zgemm
#define CPP_zherk zherk
#define CPP_data_c hmat%data_c
#endif
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_fjgj),INTENT(IN)       :: fjgj
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)          :: n,isp,jsp,iintsp,jintsp
    COMPLEX,INTENT(IN)            :: chi
    LOGICAL                       :: set0  !if true, initialize the hmat matrix with zeros
    !     .. Array Arguments ..
    CLASS(t_mat),INTENT(INOUT)     ::hmat


    INTEGER:: nn,na,ab_size,l,ll,m,size_ab_select
    INTEGER:: size_data_c,size_ab,size_ab2 !these data-dimensions are not available so easily in openacc, hence we store them
    COMPLEX,ALLOCATABLE:: ab1(:,:),ab_select(:,:)
    COMPLEX,ALLOCATABLE:: abCoeffs(:,:),ab2(:,:),h_loc(:,:),data_c(:,:)
    real :: rchi
    complex :: cchi
    CALL timestart("non-spherical setup")

    size_ab=maxval(lapw%nv)
    if (fmpi%n_size==1) Then
       size_ab_select=size_ab
    ELSE
      size_ab_select=lapw%num_local_cols(iintsp)
    ENDIF
    ALLOCATE(ab_select(size_ab_select,2*atoms%lmaxd*(atoms%lmaxd+2)+2))
    ALLOCATE(abCoeffs(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)),ab1(size_ab,2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (iintsp.NE.jintsp) THEN
       ALLOCATE(ab2(lapw%nv(iintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
       size_ab2=lapw%nv(iintsp)
    else
       allocate(ab2(1,1))
       size_ab2=1
    endif
#ifndef _OPENACC
    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       !$OMP PARALLEL DO DEFAULT(shared)
       DO l = 1, size(hmat%data_c,2)
          hmat%data_c(:,l)=0.0
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    size_data_c=size(hmat%data_c,1)
#else
    if (hmat%l_real) THEN
       allocate(data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       size_data_c=size(data_c,1)
    else
       allocate(data_c(SIZE(hmat%data_c,1),SIZE(hmat%data_c,2)))
       size_data_c=size(data_c,1)
    endif
#endif
    allocate(h_loc(SIZE(td%h_loc_nonsph,1),SIZE(td%h_loc_nonsph,1)))
    h_loc=td%h_loc_nonsph(0:,0:,n,isp,jsp)
#ifdef _OPENACC
    !$acc enter data create(ab2,ab1,abCoeffs,data_c,ab_select)copyin(h_loc)
    !$acc kernels present(data_c) default(none)
    data_c(:,:)=0.0
    !$acc end kernels
#endif
    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((sym%invsat(na)==0) .OR. (sym%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(sym%invsat(na)==0))
          cchi=MERGE(chi,chi*2,(sym%invsat(na)==0))
          CALL hsmt_ab(sym,atoms,noco,nococonv,jsp,jintsp,n,na,cell,lapw,fjgj,abCoeffs,ab_size,.TRUE.)
          !!$acc update device(ab)
          !$acc host_data use_device(abCoeffs,ab1,h_loc)
          CALL CPP_zgemm("T","N",lapw%nv(jintsp),ab_size,ab_size,cmplx(1.0,0.0),abCoeffs,SIZE(abCoeffs,1),&
                     h_loc,size(td%h_loc_nonsph,1),cmplx(0.,0.),ab1,size_ab)
          !$acc end host_data
          !ab1=MATMUL(TRANSPOSE(abCoeffs(:ab_size,:lapw%nv(iintsp))),td%h_loc(:ab_size,:ab_size,n,isp))
          !OK now of these ab1 coeffs only a part is needed in case of MPI parallelism
          !$acc kernels default(none) present(ab_select,ab1)copyin(fmpi)
          if (fmpi%n_size>1)Then
            ab_select(:,:)=ab1(fmpi%n_rank+1:lapw%nv(jintsp):fmpi%n_size,:)
          ELSE
            ab_select(:,:)=ab1(:,:) !All of ab1 needed
          ENDIF
          !$acc end kernels
          IF (iintsp==jintsp) THEN
             IF (isp==jsp) THEN
#ifdef _OPENACC
               !$acc kernels default(none) present(ab1)
               ab1(:,:)=conjg(ab1(:,:))
               !$acc end kernels
#else
               !$omp parallel do default(shared)
               DO ll = 1, SIZE(ab1,2)
                  ab1(:,ll)=conjg(ab1(:,ll))
               ENDDO
               !$omp end parallel do
#endif
               IF (fmpi%n_size==1) THEN !use z-herk trick on single PE
                 !$acc host_data use_device(data_c,ab1)
                 IF (set0 .and. nn == 1) THEN
                   !CPP_data_c = CMPLX(0.0,0.0) 
                   CALL CPP_zherk("U","N",lapw%nv(iintsp),ab_size,Rchi,ab1,size_ab,0.0,CPP_data_c,size_data_c)
                 ELSE               
                   CALL CPP_zherk("U","N",lapw%nv(iintsp),ab_size,Rchi,ab1,size_ab,1.0,CPP_data_c,size_data_c)
                 ENDIF
                 !$acc end host_data
               ELSE
                 !$acc host_data use_device(data_c,ab1,ab_select)
                 IF (set0 .and. nn == 1) THEN
                   !CPP_data_c = CMPLX(0.0,0.0) 
                   CALL CPP_zgemm("N","T",lapw%nv(iintsp),size_ab_select,ab_size,cchi,ab1,size_ab,ab_select,lapw%num_local_cols(iintsp),CMPLX(0.0,0.0),CPP_data_c,size_data_c)
                 ELSE
                   CALL CPP_zgemm("N","T",lapw%nv(iintsp),size_ab_select,ab_size,cchi,ab1,size_ab,ab_select,lapw%num_local_cols(iintsp),CMPLX(1.0,0.0),CPP_data_c,size_data_c)
                 ENDIF
                 !$acc end host_data
               ENDIF
             ELSE !This is the case of a local off-diagonal contribution.
                !It is not Hermitian, so we need to USE zgemm CALL
                CALL hsmt_ab(sym,atoms,noco,nococonv,isp,iintsp,n,na,cell,lapw,fjgj,abCoeffs,ab_size,.TRUE.)
                !!$acc update device(abCoeffs)
                !$acc kernels default(none) present(abCoeffs)
                abCoeffs(:,:)=conjg(abCoeffs(:,:))
                !$acc end kernels
                !$acc host_data use_device(abCoeffs,data_c,ab1,ab_select)
                IF (set0 .and. nn == 1) THEN
                   !CPP_data_c = CMPLX(0.0,0.0) 
                   CALL CPP_zgemm("T","T",lapw%nv(iintsp),size_ab_select,ab_size,chi,abCoeffs,SIZE(abCoeffs,1),&
                       ab_select,size_ab_select,CMPLX(0.0,0.0),CPP_data_c,SIZE_data_c)
                ELSE
                   CALL CPP_zgemm("T","T",lapw%nv(iintsp),size_ab_select,ab_size,chi,abCoeffs,SIZE(abCoeffs,1),&
                       ab_select,size_ab_select,CMPLX(1.0,0.0),CPP_data_c,SIZE_data_c)
                ENDIF
                !$acc end host_data
             ENDIF
          ELSE  !here the l_ss off-diagonal part starts
             !Second set of abCoeffs is needed
             CALL hsmt_ab(sym,atoms,noco,nococonv,isp,iintsp,n,na,cell,lapw,fjgj,abCoeffs,ab_size,.TRUE.)
             if (isp==jsp) Then
                !!$acc update device (abCoeffs)
                !$acc host_data use_device(abCoeffs,h_loc,ab2)
                CALL CPP_zgemm("T","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),abCoeffs,SIZE(abCoeffs,1),&
                     h_loc,size(td%h_loc_nonsph,1),CMPLX(0.,0.),ab2,size_ab2)
                !$acc end host_data
                !$acc kernels  default(none) present(ab2)
                ab2(:,:)=conjg(ab2(:,:))
                !$acc end kernels
                !Multiply for Hamiltonian
                !$acc host_data use_device(ab2,ab1,data_c,ab_select)
                CALL CPP_zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,chi,ab2,size_ab2,&
                     ab_select,size_ab_select,CMPLX(1.0,0.0),CPP_data_c,size_data_c)
                !$acc end host_data
             ELSE
                !$acc kernels default(none) present(abCoeffs)
                abCoeffs(:,:)=conjg(abCoeffs(:,:))
                !$acc end kernels
                !$acc host_data use_device(abCoeffs,ab1,data_c,ab_select)
                IF (set0 .and. nn == 1) THEN
                   !CPP_data_c = CMPLX(0.0,0.0) 
                   CALL CPP_zgemm("T","T",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,cchi,abCoeffs,SIZE(abCoeffs,1),&
                       ab_select,size_ab_select,CMPLX(0.0,0.0),CPP_data_c,SIZE_data_c)
                ELSE  
                   CALL CPP_zgemm("T","T",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,cchi,abCoeffs,SIZE(abCoeffs,1),&
                       ab_select,size_ab_select,CMPLX(1.0,0.0),CPP_data_c,SIZE_data_c)
                ENDIF
                !$acc end host_data
             ENDIF
          ENDIF

       END IF
    END DO
#ifdef _OPENACC
       if (hmat%l_real) THEN
          !$acc kernels present(hmat,hmat%data_r,data_c) default(none)
          hmat%data_r=hmat%data_r+real(data_c)
          !$acc end kernels
       else
          !$acc kernels present(hmat,hmat%data_c,data_c) default(none)
          hmat%data_c=hmat%data_c+data_c
          !$acc end kernels
       endif
#else
    IF (hmat%l_real) THEN
       !$OMP PARALLEL DO DEFAULT(shared)
       DO l = 1, size(hmat%data_c,2)
          hmat%data_r(:,l) = hmat%data_r(:,l) + REAL(hmat%data_c(:,l))
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
#endif
       !$acc exit data delete(ab2,ab1,abCoeffs,data_c,ab_select,h_loc)

    CALL timestop("non-spherical setup")
  END SUBROUTINE hsmt_nonsph

END MODULE m_hsmt_nonsph
