MODULE m_rhomt21
  !     ***************************************************************
  !     perform the sum over m (for each l) and bands to set up the
  !     coefficient of spherical charge densities in subroutine 
  !     cdnval                                 
  !     for offdiagonal matrix-elements in case of noncollinear magnetism 
  !     FF
  !     ***************************************************************
CONTAINS
  SUBROUTINE rhomt21(atoms,we,ne,eigVecCoeffs,uu21,ud21,du21,dd21,uulo21,dulo21,ulou21,ulod21,uloulop21)

    USE m_types_setup
    USE m_types_cdnval

    IMPLICIT NONE

    TYPE(t_atoms),       INTENT(IN)    :: atoms
    TYPE(t_eigVecCoeffs),INTENT(IN)    :: eigVecCoeffs

    !     .. Scalar Arguments ..
    INTEGER,             INTENT(IN)    :: ne 

    !     .. Array Arguments ..
    REAL,                INTENT(IN)    :: we(:)!(nobd)
    COMPLEX,             INTENT(INOUT) :: uu21(atoms%lmaxd,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: ud21(atoms%lmaxd,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: du21(atoms%lmaxd,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: dd21(atoms%lmaxd,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: uulo21(atoms%nlod,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: dulo21(atoms%nlod,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: ulou21(atoms%nlod,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: ulod21(atoms%nlod,atoms%ntype)
    COMPLEX,             INTENT(INOUT) :: uloulop21(atoms%nlod,atoms%nlod,atoms%ntype)

    !     .. Local Scalars ..
    INTEGER i,l,lm,itype,na,natom,lo,lop,m
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
                   uu21(l,itype) = uu21(l,itype) + we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,2))*eigVecCoeffs%acof(i,lm,natom,1)
                   ud21(l,itype) = ud21(l,itype) + we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,2))*eigVecCoeffs%bcof(i,lm,natom,1)
                   du21(l,itype) = du21(l,itype) + we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,2))*eigVecCoeffs%acof(i,lm,natom,1)
                   dd21(l,itype) = dd21(l,itype) + we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,2))*eigVecCoeffs%bcof(i,lm,natom,1)
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
                   uulo21(lo,itype) = uulo21(lo,itype) + we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,2))*eigVecCoeffs%ccof(m,i,lo,natom,1)
                   dulo21(lo,itype) = dulo21(lo,itype) + we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,2))*eigVecCoeffs%ccof(m,i,lo,natom,1)
                   ulou21(lo,itype) = ulou21(lo,itype) + we(i)* CONJG(eigVecCoeffs%acof(i,lm,natom,1))*eigVecCoeffs%ccof(m,i,lo,natom,2)
                   ulod21(lo,itype) = ulod21(lo,itype) + we(i)* CONJG(eigVecCoeffs%bcof(i,lm,natom,1))*eigVecCoeffs%ccof(m,i,lo,natom,2)
                ENDDO
             ENDDO
             !--->         contribution of local orbital - local orbital terms
             !--->         loop over lo'
             DO lop = 1,atoms%nlo(itype)
                IF (atoms%llo(lop,itype).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1,ne
                         uloulop21(lop,lo,itype) = uloulop21(lop,lo,itype)+&
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
