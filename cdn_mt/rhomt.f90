MODULE m_rhomt
CONTAINS
  SUBROUTINE rhomt(atoms,we,ne,eigVecCoeffs,denCoeffs,ispin)
    !     ***************************************************************
    !     perform the sum over m (for each l) and bands to set up the
    !     coefficient of spherical charge densities in subroutine
    !     cdnval                                   c.l.fu
    !     ***************************************************************

    USE m_types
    IMPLICIT NONE

    INTEGER,              INTENT(IN)    :: ne, ispin
    TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
    REAL,                 INTENT(IN)    :: we(:)!(nobd)
    TYPE(t_atoms),        INTENT(IN)    :: atoms
    TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

    INTEGER i,l,lm ,n,na,natom,m

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
                   denCoeffs%uu(l,n,ispin) = denCoeffs%uu(l,n,ispin) +&
                      we(i) * REAL(eigVecCoeffs%acof(i,lm,natom,ispin)*CONJG(eigVecCoeffs%acof(i,lm,natom,ispin)))
                   denCoeffs%dd(l,n,ispin) = denCoeffs%dd(l,n,ispin) +&
                      we(i) * REAL(eigVecCoeffs%bcof(i,lm,natom,ispin)*CONJG(eigVecCoeffs%bcof(i,lm,natom,ispin)))
                   denCoeffs%du(l,n,ispin) = denCoeffs%du(l,n,ispin) +&
                      we(i) * REAL(eigVecCoeffs%acof(i,lm,natom,ispin)*CONJG(eigVecCoeffs%bcof(i,lm,natom,ispin)))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE rhomt
END MODULE m_rhomt
