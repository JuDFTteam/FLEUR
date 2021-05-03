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
    COMPLEX, INTENT (INOUT) :: abCoeffs(:,:)
    !Optional arguments if abc coef for LOs are needed
    COMPLEX, INTENT(INOUT),OPTIONAL:: abclo(:,:,:,:)
    REAL,INTENT(IN),OPTIONAL:: alo1(:),blo1(:),clo1(:)

    INTEGER :: np,k,l,ll1,m,lmax,nkvec,lo,lm,invsfct,lmMin,lmMax,ll
    COMPLEX :: term,tempA,tempB
    REAL    :: v(3),bmrot(3,3),gkrot(3)
    COMPLEX :: ylm((atoms%lmaxd+1)**2),facA((atoms%lmaxd+1)**2),facB((atoms%lmaxd+1)**2)
    COMPLEX :: c_ph(maxval(lapw%nv),MERGE(2,1,noco%l_ss.or.any(noco%l_unrestrictMT).or.any(noco%l_spinoffd_ldau)))
    LOGICAL :: l_apw,l_abclo

    l_abclo=present(abclo)
    lmax = MERGE(atoms%lnonsph(n),atoms%lmax(n),l_nonsph)
    ab_size = lmax*(lmax+2)+1
    ! replace APW+lo check (may actually be a broken trick) by something simpler
!    l_apw=ALL(fjgj%gj==0.0)
    l_apw = .FALSE.

    ! We skip the initialization for speed
!    abCoeffs=0.0

    np = sym%invtab(sym%ngopr(na))
    CALL lapw%phase_factors(iintsp,atoms%taual(:,na),nococonv%qss,c_ph(:,iintsp))
    bmrot = transpose(MATMUL(1.*sym%mrot(:,:,np),cell%bmat))

#ifndef _OPENACC
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP& SHARED(lapw,lmax,c_ph,iintsp,abCoeffs,fjgj,abclo,cell,atoms,sym) &
    !$OMP& SHARED(l_abclo,alo1,blo1,clo1,ab_size,na,n,ispin,bmrot) &
    !$OMP& PRIVATE(k,ylm,l,ll1,m,lm,term,invsfct,lo,nkvec,facA,facB,v) &
    !$OMP& PRIVATE(gkrot,lmMin,lmMax,tempA,tempB)
#else
    !$acc kernels present(abCoeffs) default(none)
    abCoeffs(:,:)=0.0
    !$acc end kernels
#endif


    !$acc data copyin(atoms,atoms%llo,atoms%llod,atoms%nlo,cell,cell%omtil,atoms%rmt) if (l_abclo)
    !$acc parallel loop present(fjgj,fjgj%fj,fjgj%gj,abCoeffs) create(ylm, facA, facB) &
    !$acc copyin(lmax,lapw,lapw%nv,lapw%vk,lapw%kvec,bmrot,c_ph, sym, sym%invsat,l_abclo) &
    !$acc present(abclo,alo1,blo1,clo1)&
    !$acc private(gkrot,k,v,ylm,l,lm,invsfct,lo,facA,facB,term,invsfct,tempA,tempB,lmMin,lmMax,ll)  default(none)
    DO k = 1,lapw%nv(iintsp)
       !-->  apply the rotation that brings this atom into the
       !-->  representative (this is the definition of ngopr(na)
       !-->  and transform to cartesian coordinates
       v=lapw%vk(:,k,iintsp)
       !gkrot(:) = MATMUL(bmrot,v)
       gkrot(1) = bmrot(1,1)*v(1)+bmrot(1,2)*v(2)+bmrot(1,3)*v(3)
       gkrot(2) = bmrot(2,1)*v(1)+bmrot(2,2)*v(2)+bmrot(2,3)*v(3)
       gkrot(3) = bmrot(3,1)*v(1)+bmrot(3,2)*v(2)+bmrot(3,3)*v(3)

       !-->    generate spherical harmonics
       CALL ylm4(lmax,gkrot,ylm)
       !-->  synthesize the complex conjugates of a and b
       !$acc  loop vector private(l,tempA,tempB,lmMin,lmMax)
       DO l = 0,lmax
          tempA = fjgj%fj(k,l,ispin,iintsp)*c_ph(k,iintsp)
          tempB = fjgj%gj(k,l,ispin,iintsp)*c_ph(k,iintsp)
          lmMin = l*(l+1) + 1 - l
          lmMax = l*(l+1) + 1 + l
          facA(lmMin:lmMax) = tempA
          facB(lmMin:lmMax) = tempB
       END DO
       !$acc end loop

       abCoeffs(:ab_size,k)            = facA(:ab_size)*ylm(:ab_size)
       abCoeffs(ab_size+1:2*ab_size,k) = facB(:ab_size)*ylm(:ab_size)
       
       IF (l_abclo) THEN
          !determine also the abc coeffs for LOs
          invsfct=MERGE(1,2,sym%invsat(na).EQ.0)
          term = fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)*c_ph(k,iintsp)
          !!$acc loop vector private(lo,l,nkvec,ll1,m,lm)
          DO lo = 1,atoms%nlo(n)
             l = atoms%llo(lo,n)
             DO nkvec=1,invsfct*(2*l+1)
                IF (lapw%kvec(nkvec,lo,na)==k) THEN !This k-vector is used in LO
                   ll1 = l*(l+1) + 1
                   DO m = -l,l
                      lm = ll1 + m
                      abclo(1,m+atoms%llod+1,nkvec,lo) = term*ylm(lm)*alo1(lo)
                      abclo(2,m+atoms%llod+1,nkvec,lo) = term*ylm(lm)*blo1(lo)
                      abclo(3,m+atoms%llod+1,nkvec,lo) = term*ylm(lm)*clo1(lo)
                   END DO
                END IF
             ENDDO
          ENDDO
          !!$acc end loop
       ENDIF

    ENDDO !k-loop
    !$acc end parallel loop
    !$acc end data

#ifndef _OPENACC
    !$OMP END PARALLEL DO
#endif

    IF (.NOT.l_apw) ab_size=ab_size*2
  END SUBROUTINE hsmt_ab
END MODULE m_hsmt_ab
