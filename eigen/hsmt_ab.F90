!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_ab
  use m_juDFT
  implicit none

  INTERFACE hsmt_ab
    module procedure hsmt_ab_cpu
#ifdef CPP_GPU
    module procedure hsmt_ab_gpu
#endif
  END INTERFACE


CONTAINS

#ifdef CPP_GPU

  SUBROUTINE hsmt_ab_gpu(sym,atoms,noco,ispin,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,l_nonsph,abclo,alo1,blo1,clo1)
!Calculate overlap matrix
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_apws
    USE cudafor
    USE nvtx
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_noco),INTENT(IN)     :: noco
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin,n,na,iintsp
    LOGICAL,INTENT(IN)   :: l_nonsph
    INTEGER,INTENT(OUT)  :: ab_size
    !     ..
    !     .. Array Arguments ..
    REAL, DEVICE, INTENT(IN)       :: fj(:,:,:),gj(:,:,:)
    COMPLEX,DEVICE, INTENT (OUT) :: ab(:,:)
    !Optional arguments if abc coef for LOs are needed
    COMPLEX, INTENT(INOUT),OPTIONAL:: abclo(:,-atoms%llod:,:,:)
    REAL,INTENT(IN),OPTIONAL:: alo1(:),blo1(:),clo1(:)
    
    INTEGER:: np,k,l,ll1,m,lmax,nkvec,lo,lm,invsfct
    REAL   :: th,v(3),bmrot(3,3),vmult(3)
    COMPLEX,ALLOCATABLE :: ylm(:,:)
    COMPLEX,ALLOCATABLE :: c_ph(:,:)
    REAL,   ALLOCATABLE :: gkrot(:,:)
    COMPLEX:: term

    REAL,   ALLOCATABLE,DEVICE :: fj_dev(:,:,:), gj_dev(:,:,:)
    COMPLEX,ALLOCATABLE,DEVICE :: c_ph_dev(:,:)
    COMPLEX,ALLOCATABLE,DEVICE :: ylm_dev(:,:)
    REAL,   ALLOCATABLE,DEVICE :: gkrot_dev(:,:)
    INTEGER :: istat

 
    call nvtxStartRange("hsmt_ab",2)    
    lmax=MERGE(atoms%lnonsph(n),atoms%lmax(n),l_nonsph)

    ALLOCATE(c_ph_dev(lapw%nv(1),MERGE(2,1,noco%l_ss)))
    ALLOCATE(ylm_dev((lmax+1)**2,lapw%nv(1)))
    ALLOCATE(gkrot_dev(3,lapw%nv(1)))

    ALLOCATE(ylm((lmax+1)**2,lapw%nv(1)))
    ALLOCATE(c_ph(lapw%nv(1),MERGE(2,1,noco%l_ss)))
    ALLOCATE(gkrot(3,lapw%nv(1)))

    
    ab_size=lmax*(lmax+2)+1
    ab=0.0
    
    np = sym%invtab(atoms%ngopr(na))
    !--->          set up phase factors
    CALL lapw%phase_factors(iintsp,atoms%taual(:,na),noco%qss,c_ph(:,iintsp))
    c_ph_dev=c_ph   
 
    IF (np==1) THEN
       gkrot(:, 1:lapw%nv(iintsp)) = lapw%gk(:, 1:lapw%nv(iintsp),iintsp)
    ELSE
       bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
       DO k = 1,lapw%nv(iintsp)
          !-->  apply the rotation that brings this atom into the
          !-->  representative (this is the definition of ngopr(na)
          !-->  and transform to cartesian coordinates
          v(:) = lapw%vk(:,k,iintsp)
          gkrot(:,k) = MATMUL(TRANSPOSE(bmrot),v)
       END DO
    END IF

    gkrot_dev = gkrot 


    !-->  synthesize the complex conjugates of a and b

    !!$cuf kernel do <<<*,256>>>
    !DO k = 1,lapw%nv(1)
    !   !-->    generate spherical harmonics
    !   CALL ylm4_dev(lmax,gkrot_dev(:,k),ylm_dev(:,k))
    !ENDDO 

    DO k = 1,lapw%nv(1)
       call ylm4(lmax,gkrot(:,k),ylm(:,k))
    ENDDO
    ylm_dev = ylm

    call nvtxStartRange("hsmt_cuf",5)    
    !$cuf kernel do <<<*,256>>>
    DO k = 1,lapw%nv(1)
       !-->    generate spherical harmonics
       !CALL ylm4_dev(lmax,gkrot_dev(:,k),ylm_dev(:,k))
       DO l = 0,lmax
          ll1 = l* (l+1)
          DO m = -l,l               
             ab(k,ll1+m+1)         = CONJG(fj(k,l+1,iintsp)*c_ph_dev(k,iintsp)*ylm_dev(ll1+m+1,k)) 
             ab(k,ll1+m+1+ab_size) = CONJG(gj(k,l+1,iintsp)*c_ph_dev(k,iintsp)*ylm_dev(ll1+m+1,k)) 
          END DO
       END DO
    ENDDO !k-loop
    istat = cudaDeviceSynchronize() 
    call nvtxEndRange

    IF (PRESENT(abclo)) THEN
       print*, "Ooooops, TODO in hsmt_ab"
       !DO k = 1,lapw%nv(1)
       !   !determine also the abc coeffs for LOs
       !   invsfct=MERGE(1,2,atoms%invsat(na).EQ.0)
       !   term = fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)*c_ph(k,iintsp)
       !   DO lo = 1,atoms%nlo(n)
       !      l = atoms%llo(lo,n)
       !      DO nkvec=1,invsfct*(2*l+1)
       !         IF (lapw%kvec(nkvec,lo,na)==k) THEN !This k-vector is used in LO
       !            ll1 = l*(l+1) + 1
       !            DO m = -l,l
       !               lm = ll1 + m
       !               abclo(1,m,nkvec,lo) = term*ylm(k,lm)*alo1(lo)
       !               abclo(2,m,nkvec,lo) = term*ylm(k,lm)*blo1(lo)
       !               abclo(3,m,nkvec,lo) = term*ylm(k,lm)*clo1(lo)
       !            END DO
       !         END IF
       !      ENDDO
       !   ENDDO
       !ENDDO
    ENDIF
       
    ab_size=ab_size*2
    
    call nvtxEndRange
  END SUBROUTINE hsmt_ab_gpu
#endif

  SUBROUTINE hsmt_ab_cpu(sym,atoms,noco,ispin,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,l_nonsph,abclo,alo1,blo1,clo1)
!Calculate overlap matrix
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    USE m_apws
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_noco),INTENT(IN)     :: noco
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin,n,na,iintsp
    LOGICAL,INTENT(IN)   :: l_nonsph
    INTEGER,INTENT(OUT)  :: ab_size
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN)       :: fj(:,0:,:),gj(:,0:,:)
    COMPLEX, INTENT (OUT) :: ab(:,:)
    !Optional arguments if abc coef for LOs are needed
    COMPLEX, INTENT(INOUT),OPTIONAL:: abclo(:,-atoms%llod:,:,:)
    REAL,INTENT(IN),OPTIONAL:: alo1(:),blo1(:),clo1(:)
    
    INTEGER:: np,k,l,ll1,m,lmax,nkvec,lo,lm,invsfct
    COMPLEX:: term
    REAL   :: th,v(3),bmrot(3,3),vmult(3)
    COMPLEX :: ylm((atoms%lmaxd+1)**2)
    COMPLEX,ALLOCATABLE:: c_ph(:,:)
    REAL,ALLOCATABLE   :: gkrot(:,:)
    LOGICAL :: l_apw
   
    ALLOCATE(c_ph(lapw%nv(1),MERGE(2,1,noco%l_ss)))
    ALLOCATE(gkrot(3,lapw%nv(1)))

    lmax=MERGE(atoms%lnonsph(n),atoms%lmax(n),l_nonsph)
    
    ab_size=lmax*(lmax+2)+1
    l_apw=ALL(gj==0.0)
    ab=0.0
    
    np = sym%invtab(atoms%ngopr(na))
    !--->          set up phase factors
    CALL lapw%phase_factors(iintsp,atoms%taual(:,na),noco%qss,c_ph(:,iintsp))
    
    IF (np==1) THEN
       gkrot(:, 1:lapw%nv(iintsp)) = lapw%gk(:, 1:lapw%nv(iintsp),iintsp)
    ELSE
       bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
       DO k = 1,lapw%nv(iintsp)
          !-->  apply the rotation that brings this atom into the
          !-->  representative (this is the definition of ngopr(na)
          !-->  and transform to cartesian coordinates
          v(:) = lapw%vk(:,k,iintsp)
          gkrot(:,k) = MATMUL(TRANSPOSE(bmrot),v)
       END DO
    END IF
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP& SHARED(lapw,gkrot,lmax,c_ph,iintsp,ab,fj,gj,abclo,cell,atoms) &
    !$OMP& SHARED(alo1,blo1,clo1,ab_size,na,n) &
    !$OMP& PRIVATE(k,vmult,ylm,l,ll1,m,lm,term,invsfct,lo,nkvec)
    DO k = 1,lapw%nv(1)
       !-->    generate spherical harmonics
       vmult(:) =  gkrot(:,k)
       CALL ylm4(lmax,vmult,ylm)
       !-->  synthesize the complex conjugates of a and b
       DO l = 0,lmax
          ll1 = l* (l+1)
          DO m = -l,l               
             term = c_ph(k,iintsp)*ylm(ll1+m+1)
             ab(k,ll1+m+1)         = fj(k,l,iintsp)*term
             ab(k,ll1+m+1+ab_size) = gj(k,l,iintsp)*term
          END DO
       END DO
       IF (PRESENT(abclo)) THEN
          !determine also the abc coeffs for LOs
          invsfct=MERGE(1,2,atoms%invsat(na).EQ.0)
          term = fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)*c_ph(k,iintsp)
          DO lo = 1,atoms%nlo(n)
             l = atoms%llo(lo,n)
             DO nkvec=1,invsfct*(2*l+1)
                IF (lapw%kvec(nkvec,lo,na)==k) THEN !This k-vector is used in LO
                   ll1 = l*(l+1) + 1
                   DO m = -l,l
                      lm = ll1 + m
                      abclo(1,m,nkvec,lo) = term*ylm(lm)*alo1(lo)
                      abclo(2,m,nkvec,lo) = term*ylm(lm)*blo1(lo)
                      abclo(3,m,nkvec,lo) = term*ylm(lm)*clo1(lo)
                   END DO
                END IF
             ENDDO
          ENDDO
       ENDIF
       
    ENDDO !k-loop
    !$OMP END PARALLEL DO
    IF (.NOT.l_apw) ab_size=ab_size*2
    
  END SUBROUTINE hsmt_ab_cpu
END MODULE m_hsmt_ab
