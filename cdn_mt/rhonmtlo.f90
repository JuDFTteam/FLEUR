!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rhonmtlo
  !
  !***********************************************************************
  ! This subroutine is the equivalent of rhomt for the local orbital
  ! contributions to the charge.
  ! acnmt, bcnmt, ccnmt are the equivalents of uunmt, ddnmt, udnmt dunmt
  ! in rhonmt
  ! p.kurz sept. 1996
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE rhonmtlo(atoms,sphhar, ne,we,acof,bcof,ccof, acnmt,bcnmt,ccnmt)
    USE m_gaunt,ONLY:gaunt1
    USE m_types
    IMPLICIT NONE
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne   
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(:)!(nobd)
    COMPLEX, INTENT (IN) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%natd)
    REAL,    INTENT (INOUT) :: acnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntypd)
    REAL,    INTENT (INOUT) :: bcnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntypd)
    REAL,    INTENT (INOUT) :: ccnmt(atoms%nlod,atoms%nlod,sphhar%nlhd,atoms%ntypd)
    !     ..
    !     .. Local Scalars ..
    COMPLEX ci,cmv,fact,cf1
    INTEGER i,jmem,l,lh,lmp,lo,lop,lp,lpmax,lpmax0,lpmin,lpmin0,m,lpp ,mp,mpp,na,neqat0,nn,ntyp
    !     ..
    !     ..

    ci = CMPLX(0.0,1.0)

    !---> for optimal performance consider only
    !---> those combinations of l,l',l'',m,m',m'' that satisfy the three
    !---> conditions for non-zero gaunt-coeff. i.e.
    !---> |l - l''| <= l' <= l + l'' (triangular condition)
    !---> m' + m'' = m and l + l' + l'' even

    neqat0 = 0
    DO ntyp = 1,atoms%ntype
       !--->    loop over the lattice harmonics
       DO lh = 1,sphhar%nlh(atoms%ntypsy(neqat0+1))
          lpp = sphhar%llh(lh,atoms%ntypsy(neqat0+1))
          DO jmem = 1,sphhar%nmem(lh,atoms%ntypsy(neqat0+1))
             mpp = sphhar%mlh(jmem,lh,atoms%ntypsy(neqat0+1))
             cmv = CONJG(sphhar%clnu(jmem,lh,atoms%ntypsy(neqat0+1)))
             DO lo = 1,atoms%nlo(ntyp)
                l = atoms%llo(lo,ntyp)
                lpmin0 = ABS(l-lpp)
                lpmax0 = l + lpp
                !--->             check that lpmax is smaller than the max l of the
                !--->             wavefunction expansion at this atom
                lpmax = MIN(lpmax0,atoms%lmax(ntyp))
                !--->             make sure that l + l'' + lpmax is even
                lpmax = lpmax - MOD(l+lpp+lpmax,2)
                DO m = -l,l

                   !--->                add flapw - local orbital cross-terms

                   !--->                add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
                   !--->                note that gaunt1(l,lp,lpp,m,mp,mpp) computes the
                   !--->                integral of conjg(y(l,m))*y(lp,mp)*y(lpp,mpp),
                   !--->                however, since the gaunt coef. are real, this is
                   !--->                the same as int. y(l,m)*conjg(y(lp,mp)*y(lpp,mpp))
                   mp = m - mpp
                   lpmin = MAX(lpmin0,ABS(mp))
                   !--->                make sure that l + l'' + lpmin is even
                   lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                   !--->                loop over l'
                   DO lp = lpmin,lpmax,2
                      lmp = lp* (lp+1) + mp
                      fact = cmv* (ci** (l-lp))*gaunt1(l,lp,lpp,m,mp,mpp,atoms%lmaxd)
                      na = neqat0
                      DO nn = 1,atoms%neq(ntyp)
                         na = na + 1
                         DO i = 1,ne
                            cf1 = fact *  ccof(m,i,lo,na)
                            acnmt(lp,lo,lh,ntyp) =acnmt(lp,lo,lh,ntyp) + we(i) * REAL(cf1 * CONJG(acof(i,lmp,na)) )
                            bcnmt(lp,lo,lh,ntyp) =bcnmt(lp,lo,lh,ntyp) + we(i) * REAL(cf1 * CONJG(bcof(i,lmp,na)) )
                         END DO
                      END DO
                   END DO

                   !--->                add terms containing gaunt1(lp,l,lpp,mp,m,mpp)
                   mp = m + mpp
                   lpmin = MAX(lpmin0,ABS(mp))
                   !--->                make sure that l + l'' + lpmin is even
                   lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                   !--->                loop over l'
                   DO lp = lpmin,lpmax,2
                      lmp = lp* (lp+1) + mp
                      fact = cmv* (ci** (lp-l))*gaunt1(lp,l,lpp,mp,m,mpp,atoms%lmaxd)
                      na = neqat0
                      DO nn = 1,atoms%neq(ntyp)
                         na = na + 1
                         DO i = 1,ne
                            cf1 = fact * CONJG(ccof(m,i,lo,na))
                            acnmt(lp,lo,lh,ntyp) = acnmt(lp,lo,lh,ntyp) + we(i) * REAL(cf1 * acof(i,lmp,na) )
                            bcnmt(lp,lo,lh,ntyp) = bcnmt(lp,lo,lh,ntyp) + we(i) * REAL(cf1 * bcof(i,lmp,na) )
                         END DO
                      END DO
                   END DO

                   !--->                add local orbital - local orbital terms
                   DO lop = 1,atoms%nlo(ntyp)
                      lp = atoms%llo(lop,ntyp)

                      !--->                   add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
                      mp = m - mpp
                      IF ((ABS(l-lpp).LE.lp) .AND.(lp.LE. (l+lpp)) .AND.(MOD(l+lp+lpp,2).EQ.0) .AND.(ABS(mp).LE.lp)) THEN
                         fact = cmv* (ci** (l-lp))*gaunt1(l,lp,lpp,m,mp,mpp,atoms%lmaxd)
                         na = neqat0
                         DO nn = 1,atoms%neq(ntyp)
                            na = na + 1
                            DO i = 1,ne
                               ccnmt(lop,lo,lh,ntyp) =&
                                    ccnmt(lop,lo,lh,ntyp) + we(i) * REAL(fact * CONJG(ccof(mp,i,lop,na))*ccof(m ,i,lo ,na))
                            END DO
                         END DO
                      END IF

                   END DO
                END DO
             END DO
          END DO
       END DO
       neqat0 = neqat0 + atoms%neq(ntyp)
    END DO

  END SUBROUTINE rhonmtlo
END MODULE m_rhonmtlo
