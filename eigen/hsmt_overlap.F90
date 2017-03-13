!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_overlap
  use m_juDFT
  implicit none
CONTAINS

  subroutine hsmt_overlap(sym,atoms,ispin,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
!Calculate overlap matrix
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_hsmt_ab
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    TYPE(t_hamovlp),INTENT(INOUT) :: hamovlp

    INTEGER:: I(2),n,nn

    call timestart("analytic")
!    call hsmt_overlap_analytic(sym,atoms,ispin,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
    call timestop("analytic")
    hamovlp%h_c=hamovlp%s_c
    hamovlp%s_c=0.0
    call timestart("with zherk")
    call hsmt_overlap_zherk(sym,atoms,ispin,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
    call timestop("with zherk")
    print *,"Diff Overlap:",maxval(abs(hamovlp%h_c-hamovlp%s_c))
    print *,"At:",maxloc(abs(hamovlp%h_c-hamovlp%s_c))
    print *,"h_c:",maxval(abs(hamovlp%h_c))
    print *,"s_c:",maxval(abs(hamovlp%s_c))

    DO n=1,3
       DO nn=1,3
          print *,hamovlp%h_c(n,nn),hamovlp%s_c(n,nn)
       enddo
    enddo

    hamovlp%h_c=0.0
  end subroutine hsmt_overlap
 
 subroutine hsmt_overlap_zherk(sym,atoms,ispin,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
!Calculate overlap matrix
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_hsmt_ab
#ifdef _OPENACC
    USE cublas 
    USE cudafor
#endif

    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    TYPE(t_hamovlp),INTENT(INOUT) :: hamovlp
    
    INTEGER:: n,nn,na,ab_offset,l,m,lm
    COMPLEX,ALLOCATABLE:: ab(:,:)
#ifdef _OPENACC
    INTEGER :: istat
    COMPLEX,DEVICE,allocatable:: s_c(:,:)
    TYPE(cublasHandle) :: cublas_handle
    istat = cublasCreate(cublas_handle)
    ALLOCATE(s_c(size(hamovlp%s_c,1),size(hamovlp%s_c,2)))
    
    s_c=hamovlp%s_c
#endif

    ALLOCATE(ab(lapw%nv(1),2*atoms%lmaxd*(atoms%lmaxd+2)+2))
    
    
    !$acc data create(ab)
    DO n=1,atoms%ntype
       DO nn = 1,atoms%neq(n)
          na = SUM(atoms%neq(:n-1))+nn
          IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
             CALL hsmt_ab(sym,atoms,ispin,n,na,cell,lapw,gk,vk,fj,gj,ab,ab_offset)
             DO l=0,atoms%lmax(n)
                DO m=-l,l
                   lm=l*(l+1)+m
                   ab(:,ab_offset+1+lm)=sqrt(usdus%ddn(l,n,ispin))*ab(:,ab_offset+1+lm)
                ENDDO
             ENDDO
#ifdef _OPENACC
             !$acc update device(ab)
             !$acc host_data use_device(ab)
             istat = cublaszherk_v2(cublas_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_N,lapw%nv(ispin),ab_offset,1.,ab,size(ab,1),1.0,s_c,size(hamovlp%s_c,1)) 
             istat = cublaszherk_v2(cublas_handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_N,lapw%nv(ispin),ab_offset,1.,ab(1,ab_offset+1),size(ab,1),1.0,s_c,size(hamovlp%s_c,1)) 
             !$acc end host_data
#else
             CALL ZHERK("U","N",lapw%nv(ispin),ab_offset,cmplx(1.,0),ab(1,1),size(ab,1),cmplx(1.0,0.0),HamOvlp%s_c,size(HamOvlp%s_c,1))
             CALL ZHERK("U","N",lapw%nv(ispin),ab_offset,cmplx(1.,0),ab(1,ab_offset+1),size(ab,1),cmplx(1.0,0.0),HamOvlp%s_c,size(HamOvlp%s_c,1))
#endif
          ENDIF
       end DO
    end DO
    !$acc end data
#ifdef _OPENACC
    HamOvlp%s_c=s_c
#endif
  end subroutine hsmt_overlap_zherk

  subroutine hsmt_overlap_analytic(sym,atoms,ispin,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
    !Calculate overlap matrix
    USE m_constants
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    TYPE(t_hamovlp),INTENT(INOUT) :: hamovlp
    
    !     .. Local Scalars ..
    REAL ski(3),s
    INTEGER ki,kj,l,n,na,nn
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL, ALLOCATABLE :: plegend(:,:)
    COMPLEX, ALLOCATABLE :: phase(:,:)
    INTEGER,PARAMETER::ab_dim=1

  
    !     ..
    DO l = 0,atoms%lmaxd
       fleg1(l) = real(l+l+1)/real(l+1)
       fleg2(l) = real(l)/real(l+1)
       fl2p1(l) = real(l+l+1)/fpi_const
    END DO
    !
    !!$OMP PARALLEL default(NONE) &
    !!$OMP SHARED(sym,cell,atoms,lapw,usdus,ispin,gk,vk,fj,gj,hamovlp,fleg1,fl2p1,fleg2)&
    !!$OMP PRIVATE(ski,ki,kj,l,n,na,nn,plegend,phase)
    

    ALLOCATE(phase(maxval(lapw%nv),ab_dim))
    ALLOCATE(plegend(maxval(lapw%nv),0:atoms%lmaxd))
    
    !!$OMP DO  SCHEDULE(DYNAMIC,1)
    !$acc data  copyin(hamovlp) copy(hamovlp%s_c)& 
    !$acc& copyin(usdus,usdus%ddn,atoms,atoms%lmax,atoms%taual,atoms%neq,gk,lapw,lapw%nv,lapw%k1,lapw%k2,lapw%k3,fleg1,fleg2,fl2p1,fj,gj)
    !$acc parallel loop gang private(ski,plegend,phase) num_gangs(3000)
    DO  ki =  1, lapw%nv(ispin)
       ski = (/lapw%k1(ki,ispin),lapw%k2(ki,ispin),lapw%k3(ki,ispin)/) 
       !--->       legendre polynomials
       plegend=0.0
       plegend(:,0)=1.0
       DO kj = 1,ki
!          plegend(kj,1) = dot_product(gk(kj,:,ispin),gk(ki,:,ispin))
          plegend(kj,1) = gk(kj,1,ispin)*gk(ki,1,ispin)+ gk(kj,2,ispin)*gk(ki,2,ispin)+ gk(kj,3,ispin)*gk(ki,3,ispin)
       END DO 
       DO l = 1,atoms%lmaxd - 1
          plegend(:,l+1) = fleg1(l)*plegend(:,1)*plegend(:,l) - fleg2(l)*plegend(:,l-1)
       END DO
       !include factor fl2p1 already here
       DO l=0,atoms%lmaxd
          plegend(:ki,l)=plegend(:ki,l)*fl2p1(l)
       ENDDO
       !--->       loop over equivalent atoms
       DO n = 1,atoms%ntype
          phase=0.0
          DO nn = 1,atoms%neq(n)
             na = SUM(atoms%neq(:n-1))+nn
             !--->             set up phase factors
             DO kj = 1,ki
                s=(ski(1)-lapw%k1(kj,ispin))*atoms%taual(1,na)
                s=s+(ski(2)-lapw%k2(kj,ispin))*atoms%taual(2,na)
                s=s+(ski(3)-lapw%k3(kj,ispin))*atoms%taual(3,na)
                phase(kj,1)=phase(kj,1)+exp(cmplx(0.,1.*tpi_const)*s)
!                phase(kj,1)=phase(kj,1)+exp(cmplx(0.,1.*tpi_const)*dot_product(ski-(/lapw%k1(kj,ispin),lapw%k2(kj,ispin),lapw%k3(kj,ispin)/),atoms%taual(:,na)))
             END DO
          END DO
          
          !--->          update overlap and l-diagonal hamiltonian matrix
          !$acc loop collapse(2) vector    
          DO  l = 0,atoms%lmax(n)
             DO kj = 1,ki
                hamovlp%s_c(kj,ki)= hamovlp%s_c(kj,ki)+ phase(kj,1)*(plegend(kj,l)*( fj(ki,l,n,ispin)*fj(kj,l,n,ispin) + gj(ki,l,n,ispin)*usdus%ddn(l,n,ispin)*gj(kj,l,n,ispin) ))
             END DO
             !--->          end loop over l
          enddo
          !--->       end loop over atom types (ntype)
       enddo
      !--->    end loop over ki
    enddo
    !$acc end parallel
    !$acc wait
    !$acc end data
    !!$OMP END DO
    DEALLOCATE(phase,plegend) 
    !!$OMP END PARALLEL


  end subroutine hsmt_overlap_analytic

#if 1==2 
!new version from Markus
 subroutine hsmt_overlap_analytic(sym,atoms,ispin,cell,lapw,usdus,gk,vk,fj,gj,hamovlp)
    !Calculate overlap matrix
    USE m_constants
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    TYPE(t_hamovlp),INTENT(INOUT) :: hamovlp
    
    !     .. Local Scalars ..
    !REAL ski(3)
    real :: s
    INTEGER ki,kj,l,n,na,nn,idx1,idx2
    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
    
    REAL, ALLOCATABLE :: plegend(:,:,:),ski(:,:)
    COMPLEX, ALLOCATABLE :: phase(:,:,:)
    INTEGER,PARAMETER::ab_dim=1

  
    !     ..
    DO l = 0,atoms%lmaxd
       fleg1(l) = real(l+l+1)/real(l+1)
       fleg2(l) = real(l)/real(l+1)
       fl2p1(l) = real(l+l+1)/fpi_const
    END DO
    !
    !!$OMP PARALLEL default(NONE) &
    !!$OMP SHARED(sym,cell,atoms,lapw,usdus,ispin,gk,vk,fj,gj,hamovlp,fleg1,fl2p1,fleg2)&
    !!$OMP PRIVATE(ski,ki,kj,l,n,na,nn,plegend,phase)
    

    ALLOCATE(phase(maxval(lapw%nv),ab_dim,maxval(lapw%nv)))
    ALLOCATE(plegend(maxval(lapw%nv),0:atoms%lmaxd,maxval(lapw%nv))) ! Markus: think about using the max of
    ! ab_dim, atoms%lmaxd here and above, i.e. make phase and plegend same
    ! size. THis will allow to fuse the three initialization loops below
    ALLOCATE(ski(3,maxval(lapw%nv)))
    
    !!$OMP DO  SCHEDULE(DYNAMIC,1)
    !$acc data  copyin(hamovlp) copy(hamovlp%s_c)& 
    !$acc& copyin(usdus,usdus%ddn,atoms,atoms%lmax,atoms%taual,atoms%neq,gk,lapw,lapw%nv,lapw%k1,lapw%k2,lapw%k3,fleg1,fleg2,fl2p1,fj,gj)
    !$acc parallel loop gang create(ski,plegend,phase) num_gangs(3000)
    DO  ki =  1, lapw%nv(ispin)
       !$acc loop vector collapse(2)
       ! initialization loops
       do idx1=1,size(phase,1) 
          do idx2=1,size(phase,2)
             phase(idx1,idx2,ki)=0.0
          enddo
       enddo
       !$acc loop vector collapse(2)
       do idx1=1,size(plegend,1)
          do idx2=1,size(plegend,2)
             plegend(idx1,idx2,ki)=0.0
          enddo
       enddo
       !$acc loop vector
       do idx1=1,size(plegend,1)
          plegend(idx1,0,ki)=1.0
       enddo
       
       ski(:,ki) = (/lapw%k1(ki,ispin),lapw%k2(ki,ispin),lapw%k3(ki,ispin)/) 
       !--->       legendre polynomials
!       plegend(:,:,ki)=0.0
!       plegend(:,0,ki)=1.0

       !$acc loop vector
       DO kj = 1,ki
!          plegend(kj,1) = dot_product(gk(kj,:,ispin),gk(ki,:,ispin))
          plegend(kj,1,ki) = gk(kj,1,ispin)*gk(ki,1,ispin)+ gk(kj,2,ispin)*gk(ki,2,ispin)+ gk(kj,3,ispin)*gk(ki,3,ispin)
       END DO
       !$acc loop seq
       DO l = 1,atoms%lmaxd - 1
          plegend(:,l+1,ki) = fleg1(l)*plegend(:,1,ki)*plegend(:,l,ki) - fleg2(l)*plegend(:,l-1,ki)
       END DO
       !include factor fl2p1 already here
       !$acc loop vector
       DO l=0,atoms%lmaxd
          plegend(:ki,l,ki)=plegend(:ki,l,ki)*fl2p1(l)
       ENDDO
       !--->       loop over equivalent atoms
       !$acc loop collapse(force:3) ! Markus: force: is a PGI extension of OpenACC. It won't be needed
                                    ! if the prefix sum below is pre-calculated
       DO n = 1,atoms%ntype
!          phase(:,:,ki)=0.0
          DO nn = 1,atoms%neq(n)
             na = SUM(atoms%neq(:n-1))+nn ! Markus: think about pre-calculating this prefix sum operation
                                          ! if n is large.
                                          ! If precalculated, the force:3 can be reduc ed to 3 above
             !--->             set up phase factors
             DO kj = 1,ki
                s=(ski(1)-lapw%k1(kj,ispin))*atoms%taual(1,na)
                s=s+(ski(2)-lapw%k2(kj,ispin))*atoms%taual(2,na)
                s=s+(ski(3)-lapw%k3(kj,ispin))*atoms%taual(3,na)
!                phase(kj,1,ki)=phase(kj,1,ki)+exp(cmplx(0.,1.*tpi_const)*s)
                phase(kj,1,ki)=phase(kj,1,ki)+exp(cmplx(0.d0,tpi_const)*s) ! Markus: is tpi_const double already? Because s is, correct? -> reduce nr. of conversions
!                phase(kj,1)=phase(kj,1)+exp(cmplx(0.,1.*tpi_const)*dot_product(ski-(/lapw%k1(kj,ispin),lapw%k2(kj,ispin),lapw%k3(kj,ispin)/),atoms%taual(:,na)))
             END DO
          END DO
          
          !--->          update overlap and l-diagonal hamiltonian matrix
          !$acc loop collapse(2) vector    
          DO  l = 0,atoms%lmax(n)
             DO kj = 1,ki
                hamovlp%s_c(kj,ki)= hamovlp%s_c(kj,ki)+ phase(kj,1,ki)*(plegend(kj,l,ki)*( fj(ki,l,n,ispin)*fj(kj,l,n,ispin) + gj(ki,l,n,ispin)*usdus%ddn(l,n,ispin)*gj(kj,l,n,ispin) ))
             END DO
             !--->          end loop over l
          enddo
          !--->       end loop over atom types (ntype)
       enddo
      !--->    end loop over ki
    enddo
    !$acc end parallel
    !$acc wait
    !$acc end data
    !!$OMP END DO
    DEALLOCATE(phase,plegend) 
    !!$OMP END PARALLEL


  end subroutine hsmt_overlap_analytic
#endif

END MODULE m_hsmt_overlap
