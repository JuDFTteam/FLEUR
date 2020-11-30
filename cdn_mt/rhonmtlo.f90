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
  SUBROUTINE rhonmtlo(atoms,sphhar,sym,ne,we,eigVecCoeffs,denCoeffs,ispin)
    USE m_gaunt,ONLY:gaunt1
    USE m_types
    use m_constants

    IMPLICIT NONE

    TYPE(t_sphhar),       INTENT(IN)    :: sphhar
    TYPE(t_atoms),        INTENT(IN)    :: atoms
    TYPE(t_sym),          INTENT(IN)    :: sym
    TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
    TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

    INTEGER, INTENT (IN) :: ne, ispin

    REAL,    INTENT (IN) :: we(:)!(nobd)

    !     .. Local Scalars ..
    COMPLEX cmv,fact,cf1
    INTEGER i,jmem,l,lh,lmp,lo,lop,lp,lpmax,lpmax0,lpmin,lpmin0,m,lpp ,mp,mpp,na,neqat0,nn,ntyp
    !     ..
    !     ..

    !---> for optimal performance consider only
    !---> those combinations of l,l',l'',m,m',m'' that satisfy the three
    !---> conditions for non-zero gaunt-coeff. i.e.
    !---> |l - l''| <= l' <= l + l'' (triangular condition)
    !---> m' + m'' = m and l + l' + l'' even

    neqat0 = 0
    DO ntyp = 1,atoms%ntype
       !--->    loop over the lattice harmonics
       DO lh = 1,sphhar%nlh(sym%ntypsy(neqat0+1))
          lpp = sphhar%llh(lh,sym%ntypsy(neqat0+1))
          DO jmem = 1,sphhar%nmem(lh,sym%ntypsy(neqat0+1))
             mpp = sphhar%mlh(jmem,lh,sym%ntypsy(neqat0+1))
             cmv = CONJG(sphhar%clnu(jmem,lh,sym%ntypsy(neqat0+1)))
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
                      fact = cmv* (ImagUnit** (l-lp))*gaunt1(l,lp,lpp,m,mp,mpp,atoms%lmaxd)
                      na = neqat0
                      DO nn = 1,atoms%neq(ntyp)
                         na = na + 1
                         DO i = 1,ne
                            cf1 = fact *  eigVecCoeffs%ccof(m,i,lo,na,ispin)
                            denCoeffs%acnmt(lp,lo,lh,ntyp,ispin) = denCoeffs%acnmt(lp,lo,lh,ntyp,ispin) +&
                                                                   we(i) * REAL(cf1 * CONJG(eigVecCoeffs%acof(i,lmp,na,ispin)) )
                            denCoeffs%bcnmt(lp,lo,lh,ntyp,ispin) = denCoeffs%bcnmt(lp,lo,lh,ntyp,ispin) +&
                                                                   we(i) * REAL(cf1 * CONJG(eigVecCoeffs%bcof(i,lmp,na,ispin)) )
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
                      fact = cmv* (ImagUnit** (lp-l))*gaunt1(lp,l,lpp,mp,m,mpp,atoms%lmaxd)
                      na = neqat0
                      DO nn = 1,atoms%neq(ntyp)
                         na = na + 1
                         DO i = 1,ne
                            cf1 = fact * CONJG(eigVecCoeffs%ccof(m,i,lo,na,ispin))
                            denCoeffs%acnmt(lp,lo,lh,ntyp,ispin) = denCoeffs%acnmt(lp,lo,lh,ntyp,ispin) +&
                                                                   we(i) * REAL(cf1 * eigVecCoeffs%acof(i,lmp,na,ispin) )
                            denCoeffs%bcnmt(lp,lo,lh,ntyp,ispin) = denCoeffs%bcnmt(lp,lo,lh,ntyp,ispin) +&
                                                                   we(i) * REAL(cf1 * eigVecCoeffs%bcof(i,lmp,na,ispin) )
                         END DO
                      END DO
                   END DO

                   !--->                add local orbital - local orbital terms
                   DO lop = 1,atoms%nlo(ntyp)
                      lp = atoms%llo(lop,ntyp)

                      !--->                   add terms containing gaunt1(l,lp,lpp,m,mp,mpp)
                      mp = m - mpp
                      IF ((ABS(l-lpp).LE.lp) .AND.(lp.LE. (l+lpp)) .AND.(MOD(l+lp+lpp,2).EQ.0) .AND.(ABS(mp).LE.lp)) THEN
                         fact = cmv* (ImagUnit** (l-lp))*gaunt1(l,lp,lpp,m,mp,mpp,atoms%lmaxd)
                         na = neqat0
                         DO nn = 1,atoms%neq(ntyp)
                            na = na + 1
                            DO i = 1,ne
                               denCoeffs%ccnmt(lop,lo,lh,ntyp,ispin) =&
                                  denCoeffs%ccnmt(lop,lo,lh,ntyp,ispin) +&
                                  we(i) * REAL(fact * CONJG(eigVecCoeffs%ccof(mp,i,lop,na,ispin))*eigVecCoeffs%ccof(m,i,lo,na,ispin))
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
