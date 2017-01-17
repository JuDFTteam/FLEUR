MODULE m_rhomt21
  !     ***************************************************************
  !     perform the sum over m (for each l) and bands to set up the
  !     coefficient of spherical charge densities in subroutine 
  !     cdnval                                 
  !     for offdiagonal matrix-elements in case of noncollinear magnetism 
  !     FF
  !     ***************************************************************
CONTAINS
  SUBROUTINE rhomt21(atoms, we,ne,acof,bcof, ccof,mt21,lo21,uloulop21)

    USE m_types, ONLY : t_mt21,t_lo21
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne 
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: acof(:,0:,:,:)!(nobd,0:lmaxd* (lmaxd+2),natd,jspd)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:,:)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:,:,:) !(-llod:llod,nobd,nlod,natd,jspd)
    REAL,    INTENT (IN) :: we(:)!(nobd)
    TYPE (t_mt21), INTENT (INOUT) :: mt21(0:atoms%lmaxd,atoms%ntype)
    TYPE (t_lo21), INTENT (INOUT) :: lo21(atoms%nlod,atoms%ntype)
    COMPLEX,       INTENT (INOUT) :: uloulop21(atoms%nlod,atoms%nlod,atoms%ntype)
    !     ..
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
                   mt21(l,itype)%uu = mt21(l,itype)%uu + we(i)* CONJG(acof(i,lm,natom,2))*acof(i,lm,natom,1)
                   mt21(l,itype)%ud = mt21(l,itype)%ud + we(i)* CONJG(acof(i,lm,natom,2))*bcof(i,lm,natom,1)
                   mt21(l,itype)%du = mt21(l,itype)%du + we(i)* CONJG(bcof(i,lm,natom,2))*acof(i,lm,natom,1)
                   mt21(l,itype)%dd = mt21(l,itype)%dd + we(i)* CONJG(bcof(i,lm,natom,2))*bcof(i,lm,natom,1)
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
                   lo21(lo,itype)%uulo = lo21(lo,itype)%uulo + we(i)* CONJG(acof(i,lm,natom,2))*ccof(m,i,lo,natom,1)
                   lo21(lo,itype)%dulo = lo21(lo,itype)%dulo + we(i)* CONJG(bcof(i,lm,natom,2))*ccof(m,i,lo,natom,1)
                   lo21(lo,itype)%ulou = lo21(lo,itype)%ulou + we(i)* CONJG(acof(i,lm,natom,1))*ccof(m,i,lo,natom,2)
                   lo21(lo,itype)%ulod = lo21(lo,itype)%ulod + we(i)* CONJG(bcof(i,lm,natom,1))*ccof(m,i,lo,natom,2)
                ENDDO
             ENDDO
             !--->         contribution of local orbital - local orbital terms
             !--->         loop over lo'
             DO lop = 1,atoms%nlo(itype)
                IF (atoms%llo(lop,itype).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1,ne
                         uloulop21(lop,lo,itype) = uloulop21(lop,lo,itype)+&
                              we(i)*CONJG(ccof(m,i,lop,natom,2))*ccof(m,i,lo, natom,1)
                      ENDDO ! i = 1,ne
                   ENDDO   ! m = -l,l
                ENDIF
             ENDDO     ! lop
          ENDDO       ! lo

       ENDDO          ! na
    ENDDO             ! itype

  END SUBROUTINE rhomt21
END MODULE m_rhomt21
