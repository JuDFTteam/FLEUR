!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_ab
  use m_juDFT
  implicit none
CONTAINS
  subroutine hsmt_ab(sym,atoms,ispin,n,na,cell,lapw,gk,vk,fj,gj,ab,ab_offset)
!Calculate overlap matrix
#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin,n,na
    INTEGER,INTENT(OUT)  :: ab_offset
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
    COMPLEX, INTENT (OUT) :: ab(:,:)
    
    INTEGER:: np,k,l,ll1,m
    complex:: term
    real   :: th,v(3),bmrot(3,3),vmult(3)
    complex,allocatable:: c_ph(:,:),ylm(:)
    real,allocatable   :: gkrot(:,:)
    
    ALLOCATE(c_ph(lapw%nv(1),1),ylm((atoms%lmaxd+1)**2))
    allocate(gkrot(lapw%nv(1),3))


     ab_offset=atoms%lmax(n)*(atoms%lmax(n)+2)+1
     ab=0.0

     np = sym%invtab(atoms%ngopr(na))
     !--->          set up phase factors
     DO k = 1,lapw%nv(ispin)
        !th= DOT_PRODUCT((/lapw%k1(k,ispin),lapw%k2(k,ispin),lapw%k3(k,ispin)/)+(ispin-1.5)*noco%qss,atoms%taual(:,na))
        th= DOT_PRODUCT((/lapw%k1(k,ispin),lapw%k2(k,ispin),lapw%k3(k,ispin)/),atoms%taual(:,na))
        c_ph(k,ispin) = CMPLX(COS(tpi_const*th),-SIN(tpi_const*th))
     END DO
     
     IF (np==1) THEN
        gkrot( 1:lapw%nv(ispin),:) = gk( 1:lapw%nv(ispin),:,ispin)
     ELSE
        bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
        DO k = 1,lapw%nv(ispin)
           !-->  apply the rotation that brings this atom into the
           !-->  representative (this is the definition of ngopr(na)
           !-->  and transform to cartesian coordinates
           v(:) = vk(k,:,ispin)
           gkrot(k,:) = MATMUL(TRANSPOSE(bmrot),v)
        END DO
     END IF
     DO k = 1,lapw%nv(1)
        !-->    generate spherical harmonics
        vmult(:) =  gkrot(k,:)
        CALL ylm4(atoms%lmax(n),vmult,ylm)
        !-->  synthesize the complex conjugates of a and b
        DO l = 0,atoms%lmax(n)
           ll1 = l* (l+1)
           DO m = -l,l               
              term = c_ph(k,ispin)*ylm(ll1+m+1)
              ab(k,ll1+m+1)           = fj(k,l,n,ispin)*term
              ab(k,ll1+m+1+ab_offset) = gj(k,l,n,ispin)*term
           END DO
        END DO
     ENDDO !k-loop
 
end subroutine hsmt_ab
END MODULE m_hsmt_ab
