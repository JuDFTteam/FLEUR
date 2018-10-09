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
  SUBROUTINE hsmt_nonsph(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)          :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(IN)            :: chi
    !     .. Array Arguments ..
    REAL,INTENT(IN)               :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)     ::hmat
    CALL timestart("non-spherical setup")
    IF (mpi%n_size==1) THEN
#if defined (_CUDA)
       CALL priv_noMPI_gpu(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
#else
       CALL priv_noMPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
#endif
    ELSE
       CALL priv_MPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    ENDIF
    CALL timestop("non-spherical setup")
  END SUBROUTINE hsmt_nonsph

#if defined (_CUDA)
  SUBROUTINE priv_noMPI_gpu(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
  !   cublas: required to use generic BLAS interface
  !   cudafor: required to use CUDA runtime API routines 
  !   nvtx: profiling
    USE cublas   
    USE cudafor
    USE nvtx

    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)::hmat

    
    INTEGER:: nn,na,ab_size,l,ll,m
    real :: rchi
    COMPLEX,ALLOCATABLE,DEVICE :: c_dev(:,:), ab1_dev(:,:), ab_dev(:,:), ab2_dev(:,:)
    COMPLEX,ALLOCATABLE,DEVICE :: h_loc_dev(:,:)
    REAL,   ALLOCATABLE,DEVICE :: fj_dev(:,:,:), gj_dev(:,:,:)
    integer :: i, j, istat
    call nvtxStartRange("hsmt_nonsph",1)    

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2),ab1(lapw%nv(jintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
#ifdef _CUDA
    ALLOCATE(h_loc_dev(size(td%h_loc,1),size(td%h_loc,2)))
    ALLOCATE(ab1_dev(lapw%nv(jintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
    ALLOCATE(ab_dev(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
    h_loc_dev(1:,1:) = CONJG(td%h_loc(0:,0:,n,isp)) !WORKAROUND, var_dev=CONJG(var_dev) does not work 
    ALLOCATE(fj_dev(MAXVAL(lapw%nv),atoms%lmaxd+1,MERGE(2,1,noco%l_noco)))
    ALLOCATE(gj_dev(MAXVAL(lapw%nv),atoms%lmaxd+1,MERGE(2,1,noco%l_noco)))
    fj_dev(1:,1:,1:)= fj(1:,0:,1:)
    gj_dev(1:,1:,1:)= gj(1:,0:,1:)
    !note that basically all matrices in the GPU version are conjugates of their
    !cpu counterparts

    IF (iintsp.NE.jintsp) ALLOCATE(ab2_dev(lapw%nv(iintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       hmat%data_c=0.0
    ENDIF
    ALLOCATE(c_dev(SIZE(hmat%data_c,1),SIZE(hmat%data_c,2)))
    c_dev = hmat%data_c
    
    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))

          CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj_dev,gj_dev,ab_dev,ab_size,.TRUE.)

          !Calculate Hamiltonian
          CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab_dev,SIZE(ab_dev,1),&
                     h_loc_dev,SIZE(h_loc_dev,1),CMPLX(0.,0.),ab1_dev,SIZE(ab1_dev,1))
          !ab1=MATMUL(ab(:lapw%nv(iintsp),:ab_size),td%h_loc(:ab_size,:ab_size,n,isp))
          IF (iintsp==jintsp) THEN
             call nvtxStartRange("zherk",3)
             CALL ZHERK("U","N",lapw%nv(iintsp),ab_size,Rchi,ab1_dev,SIZE(ab1_dev,1),1.0,c_dev,SIZE(c_dev,1))
             istat = cudaDeviceSynchronize() 
             call nvtxEndRange()    
          ELSE  !here the l_ss off-diagonal part starts
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj_dev,gj_dev,ab_dev,ab_size,.TRUE.)
             CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab_dev,SIZE(ab_dev,1),&
                        h_loc_dev,SIZE(td%h_loc,1),CMPLX(0.,0.),ab2_dev,SIZE(ab2_dev,1))
             !Multiply for Hamiltonian
             !$cuf kernel do<<<*,256>>>
             do i = 1,size(ab1_dev,2)
               do j = 1,size(ab1_dev,1)
                  ab1_dev(j,i) = conjg(ab1_dev(j,i))
               enddo
             enddo

             CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,ab2_dev,SIZE(ab2_dev,1),&
                        ab1_dev,SIZE(ab1_dev,1),CMPLX(1.0,0.0),c_dev,SIZE(c_dev,1))
          ENDIF
       ENDIF
    END DO
    hmat%data_c = c_dev
    
    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+REAL(hmat%data_c)
    ENDIF

    call nvtxEndRange
 END SUBROUTINE priv_noMPI_gpu
#endif

  SUBROUTINE priv_noMPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm

    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)::hmat

    
    INTEGER:: nn,na,ab_size,l,ll,m
    COMPLEX,ALLOCATABLE:: ab(:,:),ab1(:,:),ab2(:,:)
    real :: rchi

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2),ab1(lapw%nv(jintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (iintsp.NE.jintsp) ALLOCATE(ab2(lapw%nv(iintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (hmat%l_real) THEN
       IF (ANY(SHAPE(hmat%data_c)/=SHAPE(hmat%data_r))) THEN
          DEALLOCATE(hmat%data_c)
          ALLOCATE(hmat%data_c(SIZE(hmat%data_r,1),SIZE(hmat%data_r,2)))
       ENDIF
       hmat%data_c=0.0
    ENDIF
    
    DO nn = 1,atoms%neq(n)
       na = SUM(atoms%neq(:n-1))+nn
       IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))

          CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian
          CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          !ab1=MATMUL(ab(:lapw%nv(iintsp),:ab_size),td%h_loc(:ab_size,:ab_size,n,isp))
          IF (iintsp==jintsp) THEN
             CALL ZHERK("U","N",lapw%nv(iintsp),ab_size,Rchi,CONJG(ab1),SIZE(ab1,1),1.0,hmat%data_c,SIZE(hmat%data_c,1))
          ELSE  !here the l_ss off-diagonal part starts
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
             CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab2,SIZE(ab2,1))
             !Multiply for Hamiltonian
             CALL zgemm("N","T",lapw%nv(iintsp),lapw%nv(jintsp),ab_size,chi,conjg(ab2),SIZE(ab2,1),ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
          ENDIF
       ENDIF
    END DO
    
    IF (hmat%l_real) THEN
       hmat%data_r=hmat%data_r+REAL(hmat%data_c)
    ENDIF

 END SUBROUTINE priv_noMPI


  SUBROUTINE priv_MPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
!Calculate overlap matrix
    USE m_hsmt_ab
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_tlmplm),INTENT(IN)   :: td
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
    COMPLEX,INTENT(in)   :: chi
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: fj(:,0:,:),gj(:,0:,:)
    CLASS(t_mat),INTENT(INOUT)::hmat

    
    INTEGER:: nn,na,ab_size,l,ll,m,i,ii
    COMPLEX,ALLOCATABLE:: ab(:,:),ab1(:,:),ab_select(:,:)
    real :: rchi

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
       IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
          rchi=MERGE(REAL(chi),REAL(chi)*2,(atoms%invsat(na)==0))
          
          CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian
        
          CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          !Cut out of ab1 only the needed elements here
          ab_select=ab1(mpi%n_rank+1:lapw%nv(jintsp):mpi%n_size,:)
          IF (iintsp==jintsp) THEN
             CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,CMPLX(rchi,0.0),CONJG(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
          ELSE
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
             CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
             !Multiply for Hamiltonian
             CALL zgemm("N","t",lapw%nv(iintsp),lapw%num_local_cols(jintsp),ab_size,chi,conjg(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(jintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))   
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
