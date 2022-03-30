!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP not_used
#else
#define CPP_OMP $OMP
#endif
MODULE m_hlomat
  IMPLICIT NONE
  !***********************************************************************
  ! updates the hamiltonian  matrix with the contributions from the local
  ! orbitals.
  ! p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE hlomat(input,atoms,fmpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,isp,jsp,&
       ntyp,na,fjgj,alo1,blo1,clo1, iintsp,jintsp,chi,hmat)
    !
#include"cpp_double.h"
    USE m_hsmt_ab
    USE m_types
!    USE m_types_mpimat
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_mpi),INTENT(IN)    :: fmpi
    TYPE(t_usdus),INTENT(IN)  :: ud
    TYPE(t_tlmplm),INTENT(IN) :: tlmplm
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_fjgj),INTENT(IN)   :: fjgj


    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,ntyp
    INTEGER, INTENT (IN) :: jsp,isp !spin for usdus and tlmplm
    INTEGER, INTENT (IN) :: jintsp,iintsp
    COMPLEX, INTENT (IN) :: chi
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT (IN) :: alo1(:,:),blo1(:,:),clo1(:,:)

    CLASS(t_mat),INTENT (INOUT) :: hmat
    !     ..
    !     .. Local Scalars ..
    COMPLEX axx,bxx,cxx,dtd,dtu,tdulo,tulod,tulou,tuloulo,utd,utu, tuulo
    INTEGER im,in,invsfct,l,lm,lmp,lo,lolo,lolop,lop,lp,i
    INTEGER mp,nkvec,nkvecp,lmplm,loplo,kp,m,mlo,mlolo
    INTEGER locol,lorow,ii,ij,n,k,ab_size,s
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: abCoeffs(:,:,:),ax(:),bx(:),cx(:)
    COMPLEX,ALLOCATABLE  :: abclo(:,:,:,:,:)
    !     ..


    !-->              synthesize the complex conjugates of a and b
    ALLOCATE(abCoeffs(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapw%nv),2))
    ALLOCATE(ax(MAXVAL(lapw%nv)),bx(MAXVAL(lapw%nv)),cx(MAXVAL(lapw%nv)))
    ALLOCATE(abclo(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod,2))

    !$acc data create(abcoeffs,abclo)
    !$acc data copyin(alo1,blo1,clo1)
    CALL hsmt_ab(sym,atoms,noco,nococonv,jsp,iintsp,ntyp,na,cell,lapw,fjgj,abCoeffs(:,:,1),ab_size,.TRUE.,abclo(:,:,:,:,1),alo1(:,isp),blo1(:,isp),clo1(:,isp))
    IF (isp==jsp.AND.iintsp==jintsp) THEN
       !$acc kernels present(abcoeffs)
       abcoeffs(:,:,2)=abcoeffs(:,:,1)
       abclo(:,:,:,:,2)=abclo(:,:,:,:,1)
       !$acc end kernels
#ifndef _OPENACC
       CALL CPP_BLAS_ccopy(SIZE(abCoeffs,1)*SIZE(abCoeffs,2),abCoeffs(:,:,1),1,abCoeffs(:,:,2),1)
       CALL CPP_BLAS_ccopy(SIZE(abclo,1)*SIZE(abclo,2)*SIZE(abclo,3)*SIZE(abclo,4),abclo(:,:,:,:,1),1,abclo(:,:,:,:,2),1)
#endif
    ELSE
       CALL hsmt_ab(sym,atoms,noco,nococonv,isp,jintsp,ntyp,na,cell,lapw,fjgj,abCoeffs(:,:,2),ab_size,.TRUE.,abclo(:,:,:,:,2),alo1(:,jsp),blo1(:,jsp),clo1(:,jsp))
    ENDIF
    !$acc end data


    mlo=0;mlolo=0
    DO m=1,ntyp-1
       mlo=mlo+atoms%nlo(m)
       mlolo=mlolo+atoms%nlo(m)*(atoms%nlo(m)+1)/2
    ENDDO

    IF ((sym%invsat(na) == 0) .OR. (sym%invsat(na) == 1)) THEN
       !--->    if this atom is the first of two atoms related by inversion,
       !--->    the contributions to the overlap matrix of both atoms are added
       !--->    at once. where it is made use of the fact, that the sum of
       !--->    these contributions is twice the real part of the contribution
       !--->    of each atom. note, that in this case there are twice as many
       !--->    (2*(2*l+1)) k-vectors (compare abccoflo and comments there).
       IF (sym%invsat(na) == 0) invsfct = 1
       IF (sym%invsat(na) == 1) invsfct = 2
       !
       !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abcoeffs,abclo) &
       !$acc & copyin(atoms,lapw,tlmplm,tlmplm%tulou,tlmplm%tulod,tlmplm%h_loc(:,:,ntyp,jsp,isp),lapw%nv(:),tlmplm%tdulo(:,:,:,jsp,isp),tlmplm%tuloulo(:,:,:,jsp,isp),atoms%rmt(ntyp))&
       !$acc & copyin(lapw%index_lo(:,na),tlmplm%h_loc2,tlmplm%tuulo(:,:,:,jsp,isp),atoms%llo(:,ntyp),atoms%nlo(ntyp),atoms%lnonsph(ntyp))&
       !$acc & copyin(ud,ud%us(:,ntyp,isp),ud%uds(:,ntyp,isp),ud%dus(:,ntyp,isp),ud%dulos(:,ntyp,isp),ud%duds(:,ntyp,isp))&
       !$acc & copyin(input, input%l_useapw, fmpi, fmpi%n_size, fmpi%n_rank)&
       !$acc & create(ax,bx,cx)&
       !$acc & default(none)
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          !--->       calculate the hamiltonian matrix elements with the regular
          !--->       flapw basis-functions
          DO m = -l,l
             lm = l* (l+1) + m
             DO kp = 1,lapw%nv(iintsp)
                ax(kp) = CMPLX(0.0,0.0)
                bx(kp) = CMPLX(0.0,0.0)
                cx(kp) = CMPLX(0.0,0.0)
             END DO
             !CALL timestart("hlomat11")
             !CPP_OMP PARALLEL DO DEFAULT(none) &
             !CPP_OMP& SHARED(ax,bx,cx,ntyp,isp,jsp,m,lm,lo,mlo) &
             !CPP_OMP& SHARED(lapw,abCoeffs,ab_size,iintsp) &
             !CPP_OMP& SHARED(atoms,tlmplm) &
             !CPP_OMP  PRIVATE(lp,mp,lmp,s)
             DO kp = 1,lapw%nv(iintsp)
                DO lp = 0,atoms%lnonsph(ntyp)
                   DO mp = -lp,lp
                      lmp = lp* (lp+1) + mp
                      s=tlmplm%h_loc2(ntyp)
                      ax(kp) = ax(kp) + CONJG(abCoeffs(lmp,kp,1))          *tlmplm%h_loc(lmp,lm,ntyp,jsp,isp)
                      ax(kp) = ax(kp) + CONJG(abCoeffs(ab_size/2+lmp,kp,1))*tlmplm%h_loc(s+lmp,lm,ntyp,jsp,isp)
                      bx(kp) = bx(kp) + CONJG(abCoeffs(lmp,kp,1))          *tlmplm%h_loc(lmp,s+lm,ntyp,jsp,isp)
                      bx(kp) = bx(kp) + CONJG(abCoeffs(ab_size/2+lmp,kp,1))*tlmplm%h_loc(s+lmp,s+lm,ntyp,jsp,isp)
                      cx(kp) = cx(kp) + CONJG(abCoeffs(lmp,kp,1))          *tlmplm%tuulo(lmp,m,lo+mlo,jsp,isp)
                      cx(kp) = cx(kp) + CONJG(abCoeffs(ab_size/2+lmp,kp,1))*tlmplm%tdulo(lmp,m,lo+mlo,jsp,isp)
                   END DO
                END DO
             END DO
             !CPP_OMP END PARALLEL DO
             !CALL timestop("hlomat11")
             !+t3e
             DO nkvec = 1,invsfct* (2*l+1)
                locol= lapw%nv(jintsp)+lapw%index_lo(lo,na)+nkvec !this is the column of the matrix
                IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN !only this MPI rank calculates this column
                   locol=(locol-1)/fmpi%n_size+1 !this is the column in local storage
                   IF (hmat%l_real) THEN
                      DO kp = 1,lapw%nv(iintsp)
                         hmat%data_r(kp,locol) = hmat%data_r(kp,locol) + chi*invsfct * (&
                              REAL(abclo(1,m,nkvec,lo,2))* REAL(ax(kp)) -&
                              AIMAG(abclo(1,m,nkvec,lo,2))*AIMAG(ax(kp)) +&
                              REAL(abclo(2,m,nkvec,lo,2))* REAL(bx(kp)) -&
                              AIMAG(abclo(2,m,nkvec,lo,2))*AIMAG(bx(kp)) +&
                              REAL(abclo(3,m,nkvec,lo,2))* REAL(cx(kp)) -&
                              AIMAG(abclo(3,m,nkvec,lo,2))*AIMAG(cx(kp)) )
                         IF (input%l_useapw) THEN
                            !---> APWlo
                            hmat%data_r(kp,locol) = hmat%data_r(kp,locol) + 0.25 * atoms%rmt(ntyp)**2 * chi*invsfct * (&
                                 (abCoeffs(lm,kp,1)* ud%us(l,ntyp,isp)+&
                                 abCoeffs(ab_size/2+lm,kp,1)*ud%uds(l,ntyp,isp))*&
                                 (CONJG(abclo(1,m,nkvec,lo,1))*  ud%dus(l,ntyp,isp)&
                                 +CONJG(abclo(2,m,nkvec,lo,1))* ud%duds(l,ntyp,isp)&
                                 +CONJG(abclo(3,m,nkvec,lo,1))*ud%dulos(lo,ntyp,isp) ))
                         ENDIF
                      ENDDO
                   ELSE
                      DO kp = 1,lapw%nv(iintsp)
                         hmat%data_c(kp,locol) = hmat%data_c(kp,locol) + CONJG(chi)*invsfct * (&
                              abclo(1,m,nkvec,lo,2) *  ax(kp)  +&
                              abclo(2,m,nkvec,lo,2) *  bx(kp)  +&
                              abclo(3,m,nkvec,lo,2) *  cx(kp)  )
                         IF (input%l_useapw) THEN
                            !---> APWlo
                            hmat%data_c(kp,locol)=hmat%data_c(kp,locol) + 0.25 * atoms%rmt(ntyp)**2 * CONJG(chi)*invsfct*(&
                                 (CONJG(abCoeffs(lm,kp,1))* ud%us(l,ntyp,isp)+&
                                 CONJG(abCoeffs(ab_size/2+lm,kp,1))*ud%uds(l,ntyp,isp))*&
                                 (abclo(1,m,nkvec,lo,2)*  ud%dus(l,ntyp,isp)&
                                 +abclo(2,m,nkvec,lo,2)* ud%duds(l,ntyp,isp)&
                                 +abclo(3,m,nkvec,lo,2)*ud%dulos(lo,ntyp,isp) ))
                         ENDIF
                      ENDDO
                   ENDIF
                   !--->             jump to the last matrixelement of the current row
                ENDIF
             END DO
          END DO
          !--->       calculate the hamiltonian matrix elements with other
          !--->       local orbitals at the same atom and with itself
          DO nkvec = 1,invsfct* (2*l+1)
             locol= lapw%nv(jintsp)+lapw%index_lo(lo,na)+nkvec !this is the column of the matrix
             IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN !only this MPI rank calculates this column
                locol=(locol-1)/fmpi%n_size+1 !this is the column in local storage
                !--->          calculate the hamiltonian matrix elements with other
                !--->          local orbitals at the same atom, if they have the same l
                DO lop = 1, MERGE(lo-1,atoms%nlo(ntyp),iintsp==jintsp)
                   IF (lop==lo) CYCLE
                   lp = atoms%llo(lop,ntyp)
                   DO nkvecp = 1,invsfct* (2*lp+1)
                      lorow=lapw%nv(iintsp)+lapw%index_lo(lop,na)+nkvecp
                      DO m = -l,l
                         lm = l* (l+1) + m
                         DO mp = -lp,lp
                            lmp = lp* (lp+1) + mp
                            s=tlmplm%h_loc2(ntyp)
                            !--->                   note, that utu,dtu... are the t-matrices and
                            !--->                   not their complex conjugates as in hssphn
                            !--->                   and that a,b,alo... are the complex
                            !--->                   conjugates of the a,b...-coefficients
                            utu=tlmplm%h_loc(lmp,lm,ntyp,jsp,isp)
                            dtu=tlmplm%h_loc(lmp+s,lm,ntyp,jsp,isp)
                            utd=tlmplm%h_loc(lmp,lm+s,ntyp,jsp,isp)
                            dtd=tlmplm%h_loc(lmp+s,lm+s,ntyp,jsp,isp)
                            tuulo = tlmplm%tuulo(lmp,m,lo+mlo,jsp,isp)
                            tdulo = tlmplm%tdulo(lmp,m,lo+mlo,jsp,isp)
                            tulou=tlmplm%tulou(lm,mp,lop+mlo,jsp,isp)
                            tulod=tlmplm%tulod(lm,mp,lop+mlo,jsp,isp)
                            !--->                         note that lo > lop
                            IF (lo>lop) THEN
                               lolop = ((lo-1)*lo)/2 + lop
                               tuloulo = tlmplm%tuloulo (m,mp,lolop+mlolo,jsp,isp)
                            ELSE
                               lolop = ((lop-1)*lop)/2 + lo
                               tuloulo = tlmplm%tuloulo (mp,m,lolop+mlolo,jsp,isp)
                            ENDIF
                            axx=abclo(1,m,nkvec,lo,2)*utu +&
                                 abclo(2,m,nkvec,lo,2)*utd +&
                                 abclo(3,m,nkvec,lo,2)*tuulo
                            bxx=abclo(1,m,nkvec,lo,2)*dtu +&
                                 abclo(2,m,nkvec,lo,2)*dtd +&
                                 abclo(3,m,nkvec,lo,2)*tdulo
                            cxx = &
                                 abclo(1,m,nkvec,lo,2)*tulou +&
                                 abclo(2,m,nkvec,lo,2)*tulod +&
                                 abclo(3,m,nkvec,lo,2)*tuloulo
                            IF (hmat%l_real) THEN
                               hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) + chi*invsfct * (&
                                    REAL(abclo(1,mp,nkvecp,lop,1))* REAL(axx) +&
                                    AIMAG(abclo(1,mp,nkvecp,lop,1))*AIMAG(axx) +&
                                    REAL(abclo(2,mp,nkvecp,lop,1))* REAL(bxx) +&
                                    AIMAG(abclo(2,mp,nkvecp,lop,1))*AIMAG(bxx) +&
                                    REAL(abclo(3,mp,nkvecp,lop,1))* REAL(cxx) +&
                                    AIMAG(abclo(3,mp,nkvecp,lop,1))*AIMAG(cxx) )
                            ELSE
                               hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) + CONJG(chi)*invsfct * (&
                                    CONJG(abclo(1,mp,nkvecp,lop,1)) * axx +&
                                    CONJG(abclo(2,mp,nkvecp,lop,1)) * bxx +&
                                    CONJG(abclo(3,mp,nkvecp,lop,1)) * cxx )
                            ENDIF
                         END DO
                      END DO
                   END DO
                END DO
                !--->          calculate the hamiltonian matrix elements of one local
                !--->          orbital with itself
                lop=lo
                DO nkvecp = 1,MERGE(nkvec,invsfct* (2*l+1),iintsp==jintsp)
                   lorow=lapw%nv(iintsp)+lapw%index_lo(lop,na)+nkvecp
                   DO m = -l,l
                      lm = l* (l+1) + m
                      DO mp = -l,l
                         lmp = l* (l+1) + mp
                         s=tlmplm%h_loc2(ntyp)
                         utu=tlmplm%h_loc(lmp,lm,ntyp,jsp,isp)
                         dtu=tlmplm%h_loc(lmp+s,lm,ntyp,jsp,isp)
                         utd=tlmplm%h_loc(lmp,lm+s,ntyp,jsp,isp)
                         dtd=tlmplm%h_loc(lmp+s,lm+s,ntyp,jsp,isp)
                         tuulo = tlmplm%tuulo(lmp,m,lo+mlo,jsp,isp)
                         tdulo = tlmplm%tdulo(lmp,m,lo+mlo,jsp,isp)
                         tulou = tlmplm%tulou(lm,mp,lo+mlo,jsp,isp)
                         tulod = tlmplm%tulod(lm,mp,lo+mlo,jsp,isp)
                         lolo = ((lo-1)*lo)/2 + lo
                         tuloulo =tlmplm%tuloulo(m,mp,lolo+mlolo,jsp,isp)
                         axx = abclo(1,m,nkvec,lo,2)*utu +&
                              abclo(2,m,nkvec,lo,2)*utd +&
                              abclo(3,m,nkvec,lo,2)*tuulo
                         bxx = abclo(1,m,nkvec,lo,2)*dtu +&
                              abclo(2,m,nkvec,lo,2)*dtd +&
                              abclo(3,m,nkvec,lo,2)*tdulo
                         cxx = abclo(1,m,nkvec,lo,2)*tulou +&
                              abclo(2,m,nkvec,lo,2)*tulod +&
                              abclo(3,m,nkvec,lo,2)*tuloulo
                         IF (hmat%l_real) THEN
                            hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) + chi*invsfct* (&
                                 REAL(abclo(1,mp,nkvecp,lo,1))* REAL(axx) +&
                                 AIMAG(abclo(1,mp,nkvecp,lo,1))*AIMAG(axx) +&
                                 REAL(abclo(2,mp,nkvecp,lo,1))* REAL(bxx) +&
                                 AIMAG(abclo(2,mp,nkvecp,lo,1))*AIMAG(bxx) +&
                                 REAL(abclo(3,mp,nkvecp,lo,1))* REAL(cxx) +&
                                 AIMAG(abclo(3,mp,nkvecp,lo,1))*AIMAG(cxx) )
                         ELSE
                            hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) + CONJG(chi)*invsfct* (&
                                 CONJG(abclo(1,mp,nkvecp,lo,1))*axx +&
                                 CONJG(abclo(2,mp,nkvecp,lo,1))*bxx +&
                                 CONJG(abclo(3,mp,nkvecp,lo,1))*cxx )
                         ENDIF
                      END DO
                   END DO
                END DO
             ENDIF !If this lo to be calculated by fmpi rank
          END DO
       END DO ! end of lo = 1,atoms%nlo loop
       !$acc end kernels
    END IF
    !$acc end data
  END SUBROUTINE hlomat
END MODULE m_hlomat
