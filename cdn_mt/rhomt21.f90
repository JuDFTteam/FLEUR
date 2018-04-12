MODULE m_rhomt21
  !     ***************************************************************
  !     perform the sum over m (for each l) and bands to set up the
  !     coefficient of spherical charge densities in subroutine 
  !     cdnval                                 
  !     for offdiagonal matrix-elements in case of noncollinear magnetism 
  !     FF
  !     ***************************************************************
CONTAINS
  SUBROUTINE rhomt21(atoms,we,ne,eigVecCoeffs,denCoeffsOffdiag)

    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)               :: atoms
    TYPE(t_eigVecCoeffs),INTENT(IN)        :: eigVecCoeffs
    TYPE(t_denCoeffsOffdiag),INTENT(INOUT) :: denCoeffsOffdiag

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne 

    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(:)!(nobd)

    !     .. Local Scalars ..
    INTEGER   i,l,lm ,itype,na,natom,lo,lop,m
    natom = 0
    DO itype = 1,atoms%ntype
       DO na = 1,atoms%neq(itype)
          natom = natom + 1
          !
          !--->       normal u, du contribution
          !
          DO l = 0,atoms%lmax(itype)
             DO m = -l,l
                lm = l* (l+1) + m
                !--->           sum over occupied bands
                DO i = 1,ne
                   denCoeffsOffdiag%uu21(l,itype) = denCoeffsOffdiag%uu21(l,itype) +&
                                                    we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,2))*eigVecCoeffs%acof(i,lm,natom,1)
                   denCoeffsOffdiag%ud21(l,itype) = denCoeffsOffdiag%ud21(l,itype) +&
                                                    we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,2))*eigVecCoeffs%bcof(i,lm,natom,1)
                   denCoeffsOffdiag%du21(l,itype) = denCoeffsOffdiag%du21(l,itype) +&
                                                    we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,2))*eigVecCoeffs%acof(i,lm,natom,1)
                   denCoeffsOffdiag%dd21(l,itype) = denCoeffsOffdiag%dd21(l,itype) +&
                                                    we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,2))*eigVecCoeffs%bcof(i,lm,natom,1)
                ENDDO ! i = 1,ne
             ENDDO   ! m = -l,l
          ENDDO     ! l
          !
          !--->       loop over the local orbitals
          !
          DO lo = 1,atoms%nlo(itype)
             l = atoms%llo(lo,itype)
             !--->         contribution of cross terms flapw - local orbitals
             DO m = -l,l
                lm = l* (l+1) + m
                DO i = 1,ne
                   denCoeffsOffdiag%uulo21(lo,itype) = denCoeffsOffdiag%uulo21(lo,itype) +&
                                                       we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,2))*eigVecCoeffs%ccof(m,i,lo,natom,1)
                   denCoeffsOffdiag%dulo21(lo,itype) = denCoeffsOffdiag%dulo21(lo,itype) +&
                                                       we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,2))*eigVecCoeffs%ccof(m,i,lo,natom,1)
                   denCoeffsOffdiag%ulou21(lo,itype) = denCoeffsOffdiag%ulou21(lo,itype) +&
                                                       we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,1))*eigVecCoeffs%ccof(m,i,lo,natom,2)
                   denCoeffsOffdiag%ulod21(lo,itype) = denCoeffsOffdiag%ulod21(lo,itype) +&
                                                       we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,1))*eigVecCoeffs%ccof(m,i,lo,natom,2)
                ENDDO
             ENDDO
             !--->         contribution of local orbital - local orbital terms
             !--->         loop over lo'
             DO lop = 1,atoms%nlo(itype)
                IF (atoms%llo(lop,itype).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1,ne
                         denCoeffsOffdiag%uloulop21(lop,lo,itype) = denCoeffsOffdiag%uloulop21(lop,lo,itype)+&
                              we(i)*CONJG(eigVecCoeffs%ccof(m,i,lop,natom,2))*eigVecCoeffs%ccof(m,i,lo, natom,1)
                      ENDDO ! i = 1,ne
                   ENDDO   ! m = -l,l
                ENDIF
             ENDDO     ! lop
          ENDDO       ! lo

       ENDDO          ! na
    ENDDO             ! itype

  END SUBROUTINE rhomt21
END MODULE m_rhomt21
