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
  SUBROUTINE rhomtlo(atoms, ne,we,acof,bcof,ccof, aclo,bclo,cclo)

    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne 
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(:)!(nobd)
    COMPLEX, INTENT (IN) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:,:)!(-atoms%llod:llod,nobd,atoms%nlod,atoms%natd)
    REAL,    INTENT (INOUT):: aclo(atoms%nlod,atoms%ntypd),bclo(atoms%nlod,atoms%ntypd)
    REAL,    INTENT (INOUT):: cclo(atoms%nlod,atoms%nlod,atoms%ntypd)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lm,lo,lop ,natom,nn,ntyp,m
    !     ..


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
                   aclo(lo,ntyp) = aclo(lo,ntyp) + we(i)*2*&
                        real(conjg(acof(i,lm,natom))*ccof(m,i,lo,natom))
                   bclo(lo,ntyp) = bclo(lo,ntyp) + we(i)*2*&
                        real(conjg(bcof(i,lm,natom))*ccof(m,i,lo,natom))
                END DO
             END DO
             !--->       contribution of local orbital - local orbital terms
             !--->       loop over lo'
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1,ne
                         cclo(lop,lo,ntyp) = cclo(lop,lo,ntyp) + we(i)*&
                              real(conjg(ccof(m,i,lop,natom))*ccof(m,i,lo ,natom))
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
    END DO


  END SUBROUTINE rhomtlo
END MODULE m_rhomtlo
