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
    COMPLEX, INTENT (INOUT) :: acof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT (INOUT) :: bcof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT (INOUT) :: ccof(-atoms%llod:,:,:,:)!(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%nat)
    REAL,    OPTIONAL, INTENT (IN)    :: fgp(3)

    !     .. Local Scalars ..
    COMPLEX ctmp,term1,work(ne)
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
    if (noco%l_noco) Then
      if (noco%l_ss) THEN
        work = ccchi(iintsp)*zMat%data_c((iintsp-1)*(lapw%nv(1)+atoms%nlotot)+nbasf,:ne)
      else
        work= ccchi(1)*zMat%data_c(nbasf,:ne)+ccchi(2)*zMat%data_c(lapw%nv(1)+atoms%nlotot+nbasf,:ne)
      ENDIF
    ELSE
      if (zmat%l_real) Then
          work=zmat%data_r(nbasf,:ne)
        else
          work=zmat%data_c(nbasf,:ne)
        endif
    endif

    !!$acc kernels default(none) present(acof,bcof,ccof,alo1,blo1,clo1,ccchi,ylm)create(ctmp) &
    !!$acc copyin(work,na,term1,l,ne,ll1,sym,sym%invsat,noco)
    !!$acc loop seq private(i,m,lm,ctmp,na2,lmp)
    DO i = 1,ne
      !!$acc loop seq
      DO m = -l,l
          lm = ll1 + m
          ctmp=term1*conjg(ylm(ll1+m+1))*work(i)
          acof(i,lm,na) = acof(i,lm,na) + ctmp*alo1(lo)
          bcof(i,lm,na) = bcof(i,lm,na) + ctmp*blo1(lo)
          ccof(m,i,lo,na) = ccof(m,i,lo,na) + ctmp*clo1(lo)
          IF (sym%invsat(na)==1.AND.noco%l_soc.AND.sym%invs) THEN
             ctmp = work(i)*CONJG(term1)*ylm(ll1+m+1)*(-1)**(l-m)
             na2 = sym%invsatnr(na)
             lmp = ll1 - m
             acof(i,lmp,na2) = acof(i,lmp,na2) +ctmp*alo1(lo)
             bcof(i,lmp,na2) = bcof(i,lmp,na2) +ctmp*blo1(lo)
             ccof(-m,i,lo,na2) = ccof(-m,i,lo,na2) +ctmp*clo1(lo)
          ENDIF
        END DO
        !!$acc end loop
    END DO
    !!$acc end loop
    !!$acc end kernels

    IF (l_force) THEN
      DO i = 1,ne
        DO m = -l,l
          lm = ll1 + m
          ctmp=term1*conjg(ylm(ll1+m+1))*work(i)
          force%acoflo(m,i,lo,na) = force%acoflo(m,i,lo,na) + ctmp*alo1(lo)
          force%bcoflo(m,i,lo,na) = force%bcoflo(m,i,lo,na) + ctmp*blo1(lo)
          DO j = 1,3
            force%aveccof(j,i,lm,na)   = force%aveccof(j,i,lm,na)   + fgp(j)*ctmp*alo1(lo)
            force%bveccof(j,i,lm,na)   = force%bveccof(j,i,lm,na)   + fgp(j)*ctmp*blo1(lo)
            force%cveccof(j,m,i,lo,na) = force%cveccof(j,m,i,lo,na) + fgp(j)*ctmp*clo1(lo)
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE abclocdn
END MODULE m_abclocdn
