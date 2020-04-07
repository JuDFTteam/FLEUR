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
  INTERFACE priv_noMPI
    module procedure priv_noMPI_cpu
#ifdef CPP_GPU
    module procedure priv_noMPI_gpu
#endif
  END INTERFACE
CONTAINS
  SUBROUTINE hsmt_nonsph(n,mpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,td,fjgj,hmat)
    USE m_hsmt_fjgj
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
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
    !     .. Array Arguments ..
    CLASS(t_mat),INTENT(INOUT)     ::hmat

#if defined CPP_GPU
    COMPLEX,ALLOCATABLE,DEVICE :: h_loc_dev(:,:)
#endif
    CALL timestart("non-spherical setup")
    IF (mpi%n_size==1) THEN
#if defined CPP_GPU
    ALLOCATE(h_loc_dev(size(td%h_loc,1),size(td%h_loc,2)))
    h_loc_dev(1:,1:) = CONJG(td%h_loc(0:,0:,n,isp,jsp))

       CALL priv_noMPI(n,mpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,h_loc_dev,fjgj,hmat)
#else
       CALL priv_noMPI(n,mpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,td,fjgj,hmat)
#endif
    ELSE
       CALL priv_MPI(n,mpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,td,fjgj,hmat)
    ENDIF
    CALL timestop("non-spherical setup")
  END SUBROUTINE hsmt_nonsph


  SUBROUTINE priv_noMPI_cpu(n,mpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,td,fjgj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_hsmt_fjgj
#ifdef _OPENACC
    USE cublas
#endif
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_nococonv),INTENT(IN) :: nococonv
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    TYPE(t_fjgj),INTENT(IN)     :: fjgj
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,jsp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    CLASS(t_mat),INTENT(INOUT)::hmat


    INTEGER:: nn,na,ab_size,l,ll,m,size_td
    INTEGER:: size_data_c,size_ab,size_ab1,size_ab2 !these data-dimensions are not available so easily in openacc, hence we store them
    COMPLEX,ALLOCATABLE:: ab(:,:),ab1(:,:),ab2(:,:),h_loc(:,:),data_c(:,:)
    real :: rchi

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2),ab1(lapw%nv(jintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
    size_ab=maxval(lapw%nv);size_ab1=lapw%nv(jintsp)

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
       hmat%data_c=0.0
    ENDIF
#else
    if (hmat%l_real) THEN
       allocate(data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       size_data_c=size(data_c,1)
    else
       allocate(data_c(SIZE(hmat%data_c,1),SIZE(hmat%data_c,2)))
       size_data_c=size(data_c,1)
    endif
    size_td=SIZE(td%h_loc,1)          
    allocate(h_loc(size_td,size_td))
    h_loc=td%h_loc(0:,0:,n,isp,jsp)
    print *,"hsmt_nonsp",size_td,shape(td%h_loc)
#endif
    !$acc data create(ab2,ab1,ab,data_c)copyin(h_loc)
    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((sym%invsat(na)==0) .OR. (sym%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(sym%invsat(na)==0))

          CALL hsmt_ab(sym,atoms,noco,nococonv,jsp,jintsp,n,na,cell,lapw,fjgj,ab,ab_size,.TRUE.)
#ifdef _OPENACC
          !$acc update device(ab)
          !$acc host_data use_device(ab,ab1,h_loc) 
          CALL cublaszgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,cmplx(1.0,0.0),ab,size_ab,&
                     h_loc,size(td%h_loc,1),cmplx(0.,0.),ab1,size_ab1)
          !$acc end host_data
#else
          call zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,cmplx(1.0,0.0),ab,size_ab,&
                     td%h_loc(0:,0:,n,isp,jsp),size(td%h_loc,1),cmplx(0.,0.),ab1,size_ab1)
#endif
          !ab1=MATMUL(ab(:lapw%nv(iintsp),:ab_size),td%h_loc(:ab_size,:ab_size,n,isp))
          IF (iintsp==jintsp) THEN
             IF (isp==jsp) THEN
#ifdef _OPENACC
                !$acc wait
                !$acc kernels present(ab1)
                ab1=conjg(ab1)
                !$acc end kernels
                !$acc host_data use_device(data_c,ab1) 
                CALL cublasZHERK("U","N",lapw%nv(iintsp),ab_size,Rchi,ab1,size_ab1,1.0,data_c,size_data_c)
                !$acc end host_data
#else
                call ZHERK("U","N",lapw%nv(iintsp),ab_size,Rchi,conjg(ab1),size_ab1,1.0,hmat%data_c,size(hmat%data_c,1))
#endif
             ELSE !This is the case of a local off-diagonal contribution.
                  !It is not Hermitian, so we need to USE zgemm CALL
                CALL hsmt_ab(sym,atoms,noco,nococonv,isp,iintsp,n,na,cell,lapw,fjgj,ab,ab_size,.TRUE.)
#ifdef _OPENACC
                !$acc update device(ab)
                !$acc kernels present(ab1)
                ab=conjg(ab)
                !$acc end kernels
                !$acc host_data use_device(ab,data_c,ab1) 
                CALL cublaszgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,ab,size_ab,&
                     ab1,size_ab1,CMPLX(1.0,0.0),data_c,SIZE_data_c)
                !$acc end host_data
#else
                CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,CONJG(ab),size_ab,&
                     ab1,SIZE_ab1,CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
#endif
             ENDIF
          ELSE  !here the l_ss off-diagonal part starts
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,nococonv,isp,iintsp,n,na,cell,lapw,fjgj,ab,ab_size,.TRUE.)
             if (isp==jsp) Then
#ifdef _OPENACC
                !$acc update device (ab)
                !$acc host_data use_device(ab,h_loc,ab2)
                CALL cublaszgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,size_ab,&
                     h_loc,size(td%h_loc,1),CMPLX(0.,0.),ab2,size_ab2,1))
                !$acc end host_data
                !$acc kernels present(ab2)
                ab2=conjg(ab2)
                !$acc end kernels
                !Multiply for Hamiltonian
                !$acc host_data use_device(ab2,ab1,data_c)
                CALL cublaszgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,ab2,size_ab2,&
                     ab1,size_ab1,CMPLX(1.0,0.0),data_c,size_data_c)
                !$acc end host_data
#else
               CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),&
               td%h_loc(0:,0:,n,isp,jsp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab2,SIZE(ab2,1))
               !Multiply for Hamiltonian
               CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,conjg(ab2),SIZE(ab2,1),&
               ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
#endif
             ELSE
#ifdef _OPENACC
                !$acc kernels present(ab)
                ab=conjg(ab)
                !$acc end kernels
                !$acc host_data use_device(ab,ab1,data_c)
                CALL cublaszgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,ab,size_ab,&
                     ab1,size_ab1,CMPLX(1.0,0.0),data_c,SIZE_data_c)
                !$acc end host_data
#else
               CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,CONJG(ab),SIZE(ab,1),&
               ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
#endif
             ENDIF
          ENDIF
         
         END IF
       end do
#ifdef _OPENACC      
       if (hmat%l_real) THEN
          !$acc kernels present(hmat,hmat%data_r)
          hmat%data_r=hmat%data_r+real(data_c)
          !$acc end kernels
       else
          !$acc kernels present(hmat,hmat%data_r)
          hmat%data_r=hmat%data_r+real(data_c)
          !$acc end kernels
       endif
#else
    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+REAL(data_c)
    ENDIF
#endif
       !$acc end data

 END SUBROUTINE priv_noMPI_cpu


  SUBROUTINE priv_MPI(n,mpi,sym,atoms,isp,jsp,iintsp,jintsp,chi,noco,nococonv,cell,lapw,td,fjgj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_nococonv),INTENT(IN) :: nococonv
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    TYPE(t_fjgj),INTENT(IN)     :: fjgj
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,jsp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    CLASS(t_mat),INTENT(INOUT)::hmat


    INTEGER:: nn,na,ab_size,l,ll,m,i,ii
    COMPLEX,ALLOCATABLE:: ab(:,:),ab1(:,:),ab_select(:,:)
    complex :: cchi

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2),ab1(lapw%nv(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2),ab_select(lapw%num_local_cols(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))

    !IF (iintsp.NE.jintsp) ALLOCATE(ab_select1(lapw%num_local_cols(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))

    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       hmat%data_c=0.0
    ENDIF

    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((sym%invsat(na)==0) .OR. (sym%invsat(na)==1)) THEN
          cchi=MERGE(chi,chi*2,(sym%invsat(na)==0))

          CALL hsmt_ab(sym,atoms,noco,nococonv,jsp,jintsp,n,na,cell,lapw,fjgj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian

          CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp,jsp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          !Cut out of ab1 only the needed elements here
          ab_select=ab1(mpi%n_rank+1:lapw%nv(jintsp):mpi%n_size,:)
          IF (iintsp==jintsp) THEN
            if (isp/=jsp) then
              CALL hsmt_ab(sym,atoms,noco,nococonv,isp,iintsp,n,na,cell,lapw,fjgj,ab,ab_size,.TRUE.)
              CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,cchi,CONJG(ab),SIZE(ab,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
            ELSE
              CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,cchi,CONJG(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
            ENDIF
          ELSE
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,nococonv,isp,iintsp,n,na,cell,lapw,fjgj,ab,ab_size,.TRUE.)
             if (isp==jsp) THEN
               CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp,jsp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
               !Multiply for Hamiltonian
               CALL zgemm("N","t",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,chi,conjg(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(jintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
             ELSE
               CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,cchi,CONJG(ab),SIZE(ab,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
             endif
          ENDIF
       ENDIF
    END DO
    !delete lower part of matrix
    !i=0
    !DO ii=mpi%n_rank+1,lapw%nv(iintsp),mpi%n_size
    !   i=i+1
    !   hmat%data_c(ii+1:,i)=0.0
    !ENDDO
    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+hmat%data_c
    ENDIF

  END SUBROUTINE priv_MPI


END MODULE m_hsmt_nonsph
