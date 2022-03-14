!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_tlo
  USE m_juDFT
  !***********************************************************************
  !     sets up the extra t-matrix elements due to the local orbitals.
  !     only non=zero elements are calculated
  !
  !     p.kurz jul. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE tlo(atoms,sym,sphhar,jspin1,jspin2,jsp,ntyp,enpara,lh0,input,vr,&
       na,flo,f,g,usdus, tlmplm, one, v1imag)
    !
    !*************** ABBREVIATIONS *****************************************
    ! tuulo      : t-matrix element of the lo and the apw radial fuction
    ! tdulo      : t-matrix element of the lo and the energy derivativ of
    !              the apw radial fuction
    ! tuloulo    : t-matrix element of two los
    !c***********************************************************************
    !
    USE m_intgr, ONLY : intgr3
    USE m_gaunt, ONLY: gaunt1
    USE m_types
    USE m_constants
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_usdus),INTENT(IN)    :: usdus
    TYPE(t_tlmplm),INTENT(INOUT):: tlmplm
    TYPE(t_enpara),INTENT(IN)   :: enpara
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin1,jspin2,jsp,ntyp ,lh0,na
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd)
    REAL,    INTENT (IN) :: f(:,:,0:,:),g(:,:,0:,:) !(atoms%jmtd,2,0:atoms%lmaxd,spins)
    REAL,    INTENT (IN) :: flo(:,:,:,:)!(atoms%jmtd,2,atoms%nlod,spins)

    COMPLEX, INTENT(IN) :: one

    TYPE(t_potden), OPTIONAL, INTENT(IN) :: v1imag

    !     ..
    !     .. Local Scalars ..
    COMPLEX cil
    INTEGER i,l,lh,lm ,lmin,lmp,lo,lop,loplo,lp,lpmax,lpmax0,lpmin,lpmin0,lpp ,mem,mp,mpp,m,lmx,mlo,mlolo
    !     ..
    !     .. Local Arrays ..
    REAL x(atoms%jmtd),ulovulo(atoms%nlod*(atoms%nlod+1)/2,lh0:sphhar%nlhd)
    REAL uvulo(atoms%nlod,0:atoms%lmaxd,lh0:sphhar%nlhd),dvulo(atoms%nlod,0:atoms%lmaxd,lh0:sphhar%nlhd)
    !     ..

    DO lo = 1,atoms%nlo(ntyp)
       l = atoms%llo(lo,ntyp)
       DO lp = 0,atoms%lmax(ntyp)
          lmin = ABS(lp-l)
          !               lmin = lp - l
          lmx = lp + l
          DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
             lpp = sphhar%llh(lh,sym%ntypsy(na))
             IF ((MOD(l+lp+lpp,2).EQ.1) .OR. (lpp.LT.lmin) .OR.&
                  (lpp.GT.lmx)) THEN
                uvulo(lo,lp,lh) = 0.0
                dvulo(lo,lp,lh) = 0.0
             ELSE
                DO i = 1,atoms%jri(ntyp)
                   x(i) = (f(i,1,lp,jspin1)*flo(i,1,lo,jspin2)+ f(i,2,lp,jspin1)*flo(i,2,lo,jspin2))*vr(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),uvulo(lo,lp,lh))
                DO i = 1,atoms%jri(ntyp)
                   x(i) = (g(i,1,lp,jspin1)*flo(i,1,lo,jspin2)+ g(i,2,lp,jspin1)*flo(i,2,lo,jspin2))*vr(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),dvulo(lo,lp,lh))
             END IF
          END DO
       END DO
    END DO
    loplo = 0
    DO lop = 1,atoms%nlo(ntyp)
       lp = atoms%llo(lop,ntyp)
       DO lo = 1,lop
          l = atoms%llo(lo,ntyp)
          loplo = loplo + 1
          IF (loplo>SIZE(ulovulo,1))  CALL juDFT_error("loplo too large!!!" ,calledby ="tlo")
          DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
             lpp = sphhar%llh(lh,sym%ntypsy(na))
             lmin = ABS(lp - l)
             lmx = lp + l
             IF ((MOD(l+lp+lpp,2).EQ.1).OR.(lpp.LT.lmin).OR.(lpp.GT.lmx)) THEN
                ulovulo(loplo,lh) = 0.0
             ELSE
                DO i = 1,atoms%jri(ntyp)
                   x(i) = (flo(i,1,lop,jspin1)*flo(i,1,lo,jspin2)+flo(i,2,lop,jspin1)*flo(i,2,lo,jspin2))*vr(i,lh)
                END DO
                CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),ulovulo(loplo,lh))
             END IF
          END DO
       END DO
    END DO
    !---> generate the different t matrices
    !---> but first initialize them ( done in eigen )
    !
    !---> generate the t-matrices. for optimal performance consider only
    !---> those combinations of l,l',l'',m,m',m'' that satisfy the three
    !---> conditions for non-zero gaunt-coeff. i.e.
    !---> |l - l''| <= l' <= l + l'' (triangular condition)
    !---> m' = m + m'' and l + l' + l'' even
    !---> loop over the local orbitals
    mlo=SUM(atoms%nlo(:ntyp-1))
    DO lo = 1,atoms%nlo(ntyp)
       l = atoms%llo(lo,ntyp)
       DO m = -l,l
          lm=l*(l+1)+m
          !--->       loop over the lattice harmonics
          DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
             lpp = sphhar%llh(lh,sym%ntypsy(na))
             lpmin0 = ABS(l-lpp)
             lpmax0 = l + lpp
             !--->          check that lpmax is smaller than the max l of the
             !--->          wavefunction expansion at this atom
             lpmax = MIN(lpmax0,atoms%lmax(ntyp))
             !--->          make sure that l + l'' + lpmax is even
             lpmax = lpmax - MOD(l+lpp+lpmax,2)
             DO mem = 1,sphhar%nmem(lh,sym%ntypsy(na))
                mpp = sphhar%mlh(mem,lh,sym%ntypsy(na))
                mp = m + mpp
                lpmin = MAX(lpmin0,ABS(mp))
                !--->             make sure that l + l'' + lpmin is even
                lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                !--->             loop over l'
                DO lp = lpmin,lpmax,2
                   lmp = lp* (lp+1) + mp
                   cil = ((ImagUnit** (l-lp))*sphhar%clnu(mem,lh,sym%ntypsy(na)))* gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                   tlmplm%tuulo(lmp,m,lo+mlo,jspin1,jspin2) = &
                        tlmplm%tuulo(lmp,m,lo+mlo,jspin1,jspin2) + one*cil*uvulo(lo,lp,lh)
                   tlmplm%tdulo(lmp,m,lo+mlo,jspin1,jspin2) = &
                        tlmplm%tdulo(lmp,m,lo+mlo,jspin1,jspin2) + one*cil*dvulo(lo,lp,lh)
                  !cil = conjg(((ImagUnit** (l-lp))*sphhar%clnu(mem,lh,sym%ntypsy(na))))* gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                  tlmplm%tulou(lmp,m,lo+mlo,jspin1,jspin2) = &
                        tlmplm%tulou(lmp,m,lo+mlo,jspin1,jspin2) + one*conjg(cil*uvulo(lo,lp,lh))
                   tlmplm%tulod(lmp,m,lo+mlo,jspin1,jspin2) = &
                        tlmplm%tulod(lmp,m,lo+mlo,jspin1,jspin2) + one*conjg(cil*dvulo(lo,lp,lh))
                END DO
             END DO
          END DO
       END DO
    END DO
    !---> generate the t-matrix including two local orbitals for lo' >= lo
    !---> loop over lo'
    mlolo=DOT_PRODUCT(atoms%nlo(:ntyp-1),atoms%nlo(:ntyp-1)+1)/2
    DO lop = 1,atoms%nlo(ntyp)
       lp = atoms%llo(lop,ntyp)
       DO mp = -lp,lp
          !--->       loop over the lattice harmonics
          DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
             lpp = sphhar%llh(lh,sym%ntypsy(na))
             DO mem = 1,sphhar%nmem(lh,sym%ntypsy(na))
                mpp = sphhar%mlh(mem,lh,sym%ntypsy(na))
                m = mp - mpp
                !--->             loop over lo
                DO lo = 1,lop
                   l = atoms%llo(lo,ntyp)
                   loplo = ((lop-1)*lop)/2 + lo
                   IF ((ABS(l-lpp).LE.lp) .AND. (lp.LE. (l+lpp)) .AND.&
                        (MOD(l+lp+lpp,2).EQ.0) .AND. (ABS(m).LE.l)) THEN
                      cil = ((ImagUnit** (l-lp))*sphhar%clnu(mem,lh,sym%ntypsy(na)))* gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                      tlmplm%tuloulo(mp,m,loplo+mlolo,jspin1,jspin2) = tlmplm%tuloulo(mp,m,loplo+mlolo,jspin1,jspin2) + one*conjg(cil*ulovulo(loplo,lh))
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO
    !---> add the diagonal terms from the muffin-tin hamiltonian. these
    !---> terms have to be made hermitian. if second variation is switched
    !---> on, the t-matrices contain only the contributions from the
    !---> non-spherical hamiltonian.
    IF (.NOT.input%secvar.AND.jspin1==jspin2) THEN
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          DO m = -l,l
             lm = l* (l+1) + m
             tlmplm%tuulo(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tuulo(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%uulon(lo,ntyp,jspin1) *&
                  ( enpara%el0(l,ntyp,jspin1)+enpara%ello0(lo,ntyp,jspin1) )
             tlmplm%tdulo(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tdulo(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%dulon(lo,ntyp,jspin1) *&
                  ( enpara%el0(l,ntyp,jspin1)+enpara%ello0(lo,ntyp,jspin1) ) + 0.5 * usdus%uulon(lo,ntyp,jspin1)
            tlmplm%tulou(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tulou(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%uulon(lo,ntyp,jspin1) *&
                  ( enpara%el0(l,ntyp,jspin1)+enpara%ello0(lo,ntyp,jspin1) )
            tlmplm%tulod(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tulod(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%dulon(lo,ntyp,jspin1) *&
                 ( enpara%el0(l,ntyp,jspin1)+enpara%ello0(lo,ntyp,jspin1) ) + 0.5 * usdus%uulon(lo,ntyp,jspin1)
            IF (atoms%ulo_der(lo,ntyp).GE.1) THEN
                tlmplm%tuulo(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tuulo(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%uuilon(lo,ntyp,jspin1)
                tlmplm%tdulo(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tdulo(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%duilon(lo,ntyp,jspin1)
                tlmplm%tulou(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tulou(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%uuilon(lo,ntyp,jspin1)
                tlmplm%tulod(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tulod(lm,m,lo+mlo,jspin1,jspin2) + 0.5 * usdus%duilon(lo,ntyp,jspin1)
             ENDIF
             !+apw_lo
             IF (atoms%l_dulo(lo,ntyp)) THEN
                tlmplm%tuulo(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tuulo(lm,m,lo+mlo,jspin1,jspin2) + 0.5
                tlmplm%tdulo(lm,m,lo+mlo,jspin1,jspin2) = 0.0
                tlmplm%tulou(lm,m,lo+mlo,jspin1,jspin2) = tlmplm%tulou(lm,m,lo+mlo,jspin1,jspin2) + 0.5
                tlmplm%tulod(lm,m,lo+mlo,jspin1,jspin2) = 0.0
             ENDIF
             !+apw_lo
          END DO
       END DO
       DO lop = 1,atoms%nlo(ntyp)
          lp = atoms%llo(lop,ntyp)
          DO lo = atoms%lo1l(lp,ntyp),lop
             loplo = ((lop-1)*lop)/2 + lo
             DO m = -lp,lp
                tlmplm%tuloulo(m,m,loplo+mlolo,jspin1,jspin2) = tlmplm%tuloulo(m,m,loplo+mlolo,jspin1,jspin2) + 0.5* (enpara%ello0(lop,ntyp,jspin1)+&
                     enpara%ello0(lo,ntyp,jspin1))* usdus%uloulopn(lop,lo,ntyp,jspin1) + 0.5* (usdus%ulouilopn(lop,lo,ntyp,jspin1) +&
                     usdus%ulouilopn(lo,lop,ntyp,jspin1))
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE tlo
END MODULE m_tlo
