!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_ab
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_ab(sym,atoms,noco,ispin,iintsp,n,na,cell,lapw,fj,gj,ab,ab_size,l_nonsph)
!Calculate overlap matrix
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
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
    
    INTEGER:: np,k,l,ll1,m,lmax
    complex:: term
    real   :: th,v(3),bmrot(3,3),vmult(3)
    complex,allocatable:: c_ph(:,:),ylm(:)
    real,allocatable   :: gkrot(:,:)
    LOGICAL :: l_apw
    
    ALLOCATE(c_ph(lapw%nv(1),1),ylm((atoms%lmaxd+1)**2))
    ALLOCATE(gkrot(3,lapw%nv(1)))

    lmax=MERGE(atoms%lnonsph(n),atoms%lmax(n),l_nonsph)
    
    ab_size=lmax*(lmax+2)+1
    l_apw=ALL(gj==0.0)
    ab=0.0
    
    np = sym%invtab(atoms%ngopr(na))
    !--->          set up phase factors
    DO k = 1,lapw%nv(iintsp)
       th= DOT_PRODUCT(lapw%gvec(:,k,iintsp)+(iintsp-1.5)*noco%qss,atoms%taual(:,na))
       c_ph(k,iintsp) = CMPLX(COS(tpi_const*th),-SIN(tpi_const*th))
    END DO
    
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
    DO k = 1,lapw%nv(1)
       !-->    generate spherical harmonics
       vmult(:) =  gkrot(:,k)
       CALL ylm4(lmax,vmult,ylm)
       !-->  synthesize the complex conjugates of a and b
       DO l = 0,lmax
          ll1 = l* (l+1)
          DO m = -l,l               
             term = c_ph(k,iintsp)*ylm(ll1+m+1)
             ab(k,ll1+m+1)         = fj(k,l,ispin)*term
             ab(k,ll1+m+1+ab_size) = gj(k,l,ispin)*term
          END DO
       END DO
    ENDDO !k-loop
    IF (.NOT.l_apw) ab_size=ab_size*2
    
  END SUBROUTINE hsmt_ab
END MODULE m_hsmt_ab
