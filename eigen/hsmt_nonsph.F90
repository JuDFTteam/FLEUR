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
       CALL priv_noMPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    ELSE
       CALL priv_MPI(n,mpi,sym,atoms,isp,iintsp,jintsp,chi,noco,cell,lapw,td,fj,gj,hmat)
    ENDIF
    CALL timestop("non-spherical setup")
  END SUBROUTINE hsmt_nonsph

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

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lmaxd*(atoms%lmaxd+2)+2),ab1(lapw%nv(iintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

    IF (iintsp.NE.jintsp) ALLOCATE(ab2(lapw%nv(jintsp),2*atoms%lmaxd*(atoms%lmaxd+2)+2))

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
          
          CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian
          CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          !ab1=MATMUL(ab(:lapw%nv(iintsp),:ab_size),td%h_loc(:ab_size,:ab_size,n,isp))
          IF (iintsp==jintsp) THEN
             CALL ZHERK("U","N",lapw%nv(iintsp),ab_size,Rchi,CONJG(ab1),SIZE(ab1,1),1.0,hmat%data_c,SIZE(hmat%data_c,1))
          ELSE  !here the l_ss off-diagonal part starts
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
             CALL zgemm("N","N",lapw%nv(jintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab2,SIZE(ab2,1))
             !Multiply for Hamiltonian
             CALL zgemm("N","T",lapw%nv(jintsp),lapw%nv(iintsp),ab_size,chi,conjg(ab2),SIZE(ab2,1),ab1,SIZE(ab1,1),CMPLX(1.0,0.0),hmat%data_c,SIZE(hmat%data_c,1))
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

    ALLOCATE(ab(MAXVAL(lapw%nv),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2),ab1(lapw%nv(iintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2),ab_select(lapw%num_local_cols(jintsp),2*atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+2))

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
          
          CALL hsmt_ab(sym,atoms,noco,isp,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
          !Calculate Hamiltonian
        
          CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(0:,0:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
          !Cut out of ab1 only the needed elements here
          ab_select=ab1(mpi%n_rank+1:lapw%nv(iintsp):mpi%n_size,:)
          IF (iintsp==jintsp) THEN
             CALL zgemm("N","T",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,CMPLX(rchi,0.0),CONJG(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))
          ELSE
             !Second set of ab is needed
             CALL hsmt_ab(sym,atoms,noco,isp,jintsp,n,na,cell,lapw,fj,gj,ab,ab_size,.TRUE.)
             CALL zgemm("N","N",lapw%nv(iintsp),ab_size,ab_size,CMPLX(1.0,0.0),ab,SIZE(ab,1),td%h_loc(:,:,n,isp),SIZE(td%h_loc,1),CMPLX(0.,0.),ab1,SIZE(ab1,1))
             !Multiply for Hamiltonian
             CALL zgemm("N","t",lapw%nv(iintsp),lapw%num_local_cols(iintsp),ab_size,chi,conjg(ab1),SIZE(ab1,1),ab_select,lapw%num_local_cols(iintsp),CMPLX(1.,0.0),hmat%data_c,SIZE(hmat%data_c,1))   
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
