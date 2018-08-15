!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abclocdn
  USE m_juDFT
  !*********************************************************************
  ! Calculates the (upper case) A, B and C coefficients for the local
  ! orbitals. The difference to abccoflo is, that a summation over the
  ! Gs ist performed. The A, B and C coeff. are set up for each eigen-
  ! state.
  ! Philipp Kurz 99/04
  !*********************************************************************
  !*************** ABBREVIATIONS ***************************************
  ! nkvec   : stores the number of G-vectors that have been found and
  !           accepted during the construction of the local orbitals.
  ! kvec    : k-vector used in hssphn to attach the local orbital 'lo'
  !           of atom 'na' to it.
  !*********************************************************************
CONTAINS
  SUBROUTINE abclocdn(atoms,sym,noco,lapw,cell,ccchi,iintsp,phase,ylm,&
       ntyp,na,k,nkvec,lo,ne,alo1,blo1,clo1,acof,bcof,ccof,zMat,l_force,fgp,force)

    USE m_types
    USE m_constants

    IMPLICIT NONE

    TYPE(t_noco),  INTENT(IN) :: noco
    TYPE(t_sym),   INTENT(IN) :: sym
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_lapw),  INTENT(IN) :: lapw
    TYPE(t_cell),  INTENT(IN) :: cell
    TYPE(t_mat),   INTENT(IN) :: zMat
    TYPE(t_force), OPTIONAL, INTENT(INOUT) :: force

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: iintsp
    INTEGER, INTENT (IN) :: k,na,ne,ntyp,nkvec,lo
    COMPLEX, INTENT (IN) :: phase
    LOGICAL, INTENT (IN) :: l_force

    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: alo1(:),blo1(:),clo1(:)
    COMPLEX, INTENT (IN) :: ylm( (atoms%lmaxd+1)**2 )
    COMPLEX, INTENT (IN) :: ccchi(2)
    COMPLEX, INTENT (INOUT) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: ccof(-atoms%llod:,:,:,:)!(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%nat)
    REAL,    OPTIONAL, INTENT (IN)    :: fgp(3)

    !     .. Local Scalars ..
    COMPLEX ctmp,term1
    INTEGER i,j,l,ll1,lm,nbasf,m,na2,lmp
    !     ..
    !     ..
    term1 = 2 * tpi_const/SQRT(cell%omtil) * ((atoms%rmt(ntyp)**2)/2) * phase
    !---> the whole program is in hartree units, therefore 1/wronskian is
    !---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
    !---> and c coefficients, is included in the t-matrices. thus, it does
    !---> not show up in the formula above.
    l = atoms%llo(lo,ntyp)
    ll1 = l* (l+1)
    nbasf=lapw%nv(iintsp)+lapw%index_lo(lo,na)+nkvec
    DO i = 1,ne
       DO m = -l,l
          lm = ll1 + m
          !+gu_con
          IF (noco%l_noco) THEN
             IF (noco%l_ss) THEN
                ctmp = term1*CONJG(ylm(ll1+m+1))*ccchi(iintsp)*zMat%data_c((iintsp-1)*(lapw%nv(1)+atoms%nlotot)+nbasf,i)
             ELSE
                ctmp = term1*CONJG(ylm(ll1+m+1))*( ccchi(1)*zMat%data_c(nbasf,i)+ccchi(2)*zMat%data_c(lapw%nv(1)+atoms%nlotot+nbasf,i) )
             ENDIF
          ELSE
             IF (zMat%l_real) THEN
                ctmp = zMat%data_r(nbasf,i)*term1*CONJG(ylm(ll1+m+1))
             ELSE
                ctmp = zMat%data_c(nbasf,i)*term1*CONJG(ylm(ll1+m+1))
             ENDIF
          ENDIF
          acof(i,lm,na) = acof(i,lm,na) + ctmp*alo1(lo)
          bcof(i,lm,na) = bcof(i,lm,na) + ctmp*blo1(lo)
          ccof(m,i,lo,na) = ccof(m,i,lo,na) + ctmp*clo1(lo)
          IF (atoms%invsat(ntyp)==1.AND.noco%l_soc.AND.sym%invs) THEN
             ctmp = zMat%data_c(nbasf,i)*CONJG(term1)*ylm(ll1+m+1)*(-1)**(l-m)
             na2 = sym%invsatnr(na)
             lmp = ll1 - m
             acof(i,lmp,na2) = acof(i,lmp,na2) +ctmp*alo1(lo)
             bcof(i,lmp,na2) = bcof(i,lmp,na2) +ctmp*blo1(lo)
             ccof(-m,i,lo,na2) = ccof(-m,i,lo,na2) +ctmp*clo1(lo)
          ENDIF
          IF (l_force) THEN
             force%acoflo(m,i,lo,na) = force%acoflo(m,i,lo,na) + ctmp*alo1(lo)
             force%bcoflo(m,i,lo,na) = force%bcoflo(m,i,lo,na) + ctmp*blo1(lo)
             DO j = 1,3
                force%aveccof(j,i,lm,na)   = force%aveccof(j,i,lm,na)   + fgp(j)*ctmp*alo1(lo)
                force%bveccof(j,i,lm,na)   = force%bveccof(j,i,lm,na)   + fgp(j)*ctmp*blo1(lo)
                force%cveccof(j,m,i,lo,na) = force%cveccof(j,m,i,lo,na) + fgp(j)*ctmp*clo1(lo)
             END DO
          END IF
       END DO
    END DO 
  
  END SUBROUTINE abclocdn
END MODULE m_abclocdn
