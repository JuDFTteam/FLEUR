!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_ab
  use m_juDFT
  implicit none

CONTAINS

  SUBROUTINE hsmt_ab(sym,atoms,noco,nococonv,ispin,iintsp,n,na,cell,lapw,fjgj,abCoeffs,ab_size,l_nonsph,abclo,alo1,blo1,clo1)
!Calculate overlap matrix, CPU vesion
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm 
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_nococonv),INTENT(IN) :: nococonv
    TYPE(t_fjgj),INTENT(IN)     :: fjgj
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ispin,n,na,iintsp
    LOGICAL,INTENT(IN)   :: l_nonsph
    INTEGER,INTENT(OUT)  :: ab_size
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (OUT) :: abCoeffs(:,:)
    !Optional arguments if abc coef for LOs are needed
    COMPLEX, INTENT(INOUT),OPTIONAL:: abclo(:,-atoms%llod:,:,:)
    REAL,INTENT(IN),OPTIONAL:: alo1(:),blo1(:),clo1(:)

    INTEGER:: np,k,l,ll1,m,lmax,nkvec,lo,lm,invsfct
    COMPLEX:: term
    REAL   :: th,v(3),bmrot(3,3),vmult(3)
    COMPLEX :: ylm((atoms%lmaxd+1)**2)
    COMPLEX,ALLOCATABLE:: c_ph(:,:)
    REAL,ALLOCATABLE   :: gkrot(:,:)
    LOGICAL :: l_apw, l_pres_abclo

    ALLOCATE(c_ph(maxval(lapw%nv),MERGE(2,1,noco%l_ss.or.noco%l_mtNocoPot)))
    ALLOCATE(gkrot(3,maxval(lapw%nv)))

    lmax=MERGE(atoms%lnonsph(n),atoms%lmax(n),l_nonsph)

    ab_size=lmax*(lmax+2)+1
    l_apw=ALL(fjgj%gj==0.0)
    abCoeffs=0.0

    np = sym%invtab(sym%ngopr(na))
    !--->          set up phase factors
    CALL lapw%phase_factors(iintsp,atoms%taual(:,na),nococonv%qss,c_ph(:,iintsp))

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
    l_pres_abclo = PRESENT(abclo)
#ifndef _OPENACC
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP& SHARED(lapw,gkrot,lmax,c_ph,iintsp,abCoeffs,fjgj,abclo,cell,atoms,sym) &
    !$OMP& SHARED(alo1,blo1,clo1,ab_size,na,n,ispin,l_pres_abclo) &
    !$OMP& PRIVATE(k,vmult,ylm,l,ll1,m,lm,term,invsfct,lo,nkvec)
#else
    !$acc kernels present(abCoeffs)
    abCoeffs(:,:)=0.0
    !$acc end kernels
#endif
    
    !$acc parallel loop present(fjgj%fj,fjgj%gj,abCoeffs) private(vmult,k,ylm,lm,invsfct,lo,nkvec) 
    DO k = 1,lapw%nv(iintsp)
       !-->    generate spherical harmonics
       vmult(:) =  gkrot(:,k)
       CALL ylm4(lmax,vmult,ylm)
       !-->  synthesize the complex conjugates of a and b
       !$acc  loop vector private(l,ll1,m,term)
       DO l = 0,lmax 
          ll1 = l* (l+1)
          DO m = -l,l
             term = c_ph(k,iintsp)*ylm(ll1+m+1)
             abCoeffs(ll1+m+1,k)         = fjgj%fj(l,k,ispin,iintsp)*term
             abCoeffs(ll1+m+1+ab_size,k) = fjgj%gj(l,k,ispin,iintsp)*term
          END DO
       END DO
       !$acc end loop 
       IF (l_pres_abclo) THEN
          !determine also the abc coeffs for LOs
          invsfct=MERGE(1,2,sym%invsat(na).EQ.0)
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
    !$acc end parallel loop
#ifndef _OPENACC
    !$OMP END PARALLEL DO
#endif
    IF (.NOT.l_apw) ab_size=ab_size*2

  END SUBROUTINE hsmt_ab
END MODULE m_hsmt_ab
