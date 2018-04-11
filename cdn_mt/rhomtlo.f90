!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rhomtlo
  !
  !***********************************************************************
  ! This subroutine is the equivalent of rhomt for the local orbital
  ! contributions to the charge.
  ! aclo,bclo,cclo are the equivalents of uu,ud,dd in rhomt
  ! p.kurz sept. 1996
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE rhomtlo(atoms,ne,we,eigVecCoeffs,denCoeffs,ispin)

    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),        INTENT(IN)    :: atoms
    TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
    TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

    INTEGER, INTENT (IN) :: ne, ispin

    REAL,    INTENT (IN) :: we(:)!(nobd)

    INTEGER i,l,lm,lo,lop ,natom,nn,ntyp,m

    natom = 0
    !---> loop over atoms
    DO ntyp = 1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          !--->       loop over the local orbitals
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             !--->       contribution of cross terms flapw - local orbitals
             DO m = -l,l
                lm = l* (l+1) + m
                DO i = 1,ne
                   denCoeffs%aclo(lo,ntyp,ispin) = denCoeffs%aclo(lo,ntyp,ispin) + we(i)*2*&
                        real(conjg(eigVecCoeffs%acof(i,lm,natom,ispin))*eigVecCoeffs%ccof(m,i,lo,natom,ispin))
                   denCoeffs%bclo(lo,ntyp,ispin) = denCoeffs%bclo(lo,ntyp,ispin) + we(i)*2*&
                        real(conjg(eigVecCoeffs%bcof(i,lm,natom,ispin))*eigVecCoeffs%ccof(m,i,lo,natom,ispin))
                END DO
             END DO
             !--->       contribution of local orbital - local orbital terms
             !--->       loop over lo'
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1,ne
                         denCoeffs%cclo(lop,lo,ntyp,ispin) = denCoeffs%cclo(lop,lo,ntyp,ispin) + we(i)*&
                              real(conjg(eigVecCoeffs%ccof(m,i,lop,natom,ispin))*eigVecCoeffs%ccof(m,i,lo,natom,ispin))
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
    END DO


  END SUBROUTINE rhomtlo
END MODULE m_rhomtlo
