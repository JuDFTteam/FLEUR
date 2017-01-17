MODULE m_rhomt
CONTAINS
  SUBROUTINE rhomt(atoms, we,ne,acof,bcof,uu,dd,du)
    !     ***************************************************************
    !     perform the sum over m (for each l) and bands to set up the
    !     coefficient of spherical charge densities in subroutine
    !     cdnval                                   c.l.fu
    !     ***************************************************************

    USE m_types
    IMPLICIT NONE

    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..  
    INTEGER, INTENT (IN) :: ne 
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: acof(:,0:,:)!(nobd,0:lmaxd* (lmaxd+2),natd)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:)
    REAL,    INTENT (IN) :: we(:)!(nobd)
    REAL, INTENT (INOUT) :: dd(0:atoms%lmaxd,atoms%ntype)
    REAL, INTENT (INOUT) :: du(0:atoms%lmaxd,atoms%ntype)
    REAL, INTENT (INOUT) :: uu(0:atoms%lmaxd,atoms%ntype)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lm ,n,na,natom,m
    !     ..
    natom = 0
    DO n = 1,atoms%ntype
       DO na = 1,atoms%neq(n)
          natom = natom + 1
          DO l = 0,atoms%lmax(n)
             !     -----> sum over m
             DO m = -l,l
                lm = l* (l+1) + m
                !     -----> sum over occupied bands
                DO i = 1,ne
                   uu(l,n) = uu(l,n) + we(i)* REAL(acof(i,lm,natom)*CONJG(acof(i,lm,natom)))
                   dd(l,n) = dd(l,n) + we(i)* REAL(bcof(i,lm,natom)*CONJG(bcof(i,lm,natom)))
                   du(l,n) = du(l,n) + we(i)* REAL(acof(i,lm,natom)*CONJG(bcof(i,lm,natom)))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE rhomt
END MODULE m_rhomt
