MODULE m_orbmom
  !     ***************************************************************
  !     perform the sum over m (for each l) and bands to set up the
  !     coefficient of spherical contribution to orbital moment.
  !     all quantities are in the local spin-frame
  !     ***************************************************************

CONTAINS
  SUBROUTINE orbmom(atoms,ne,we,acof,bcof, ccof, orb,orbl,orblo)

    !USE m_types, ONLY : t_orb,t_orbl,t_orblo
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: acof(:,0:,:) !(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:) !(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:,:)!(-atoms%llod:llod,nobd,atoms%nlod,atoms%nat)
    REAL,    INTENT (IN) :: we(:)!(nobd)
    TYPE (t_orb),  INTENT (INOUT) :: orb(0:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%ntype)
    TYPE (t_orbl), INTENT (INOUT) :: orbl(atoms%nlod,-atoms%llod:atoms%llod,atoms%ntype)
    TYPE (t_orblo),INTENT (INOUT) :: orblo(atoms%nlod,atoms%nlod,-atoms%llod:atoms%llod,atoms%ntype)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lm ,n,na,natom,ilo,ilop,m
    COMPLEX,PARAMETER:: czero= CMPLX(0.0,0.0)

    natom = 0
    DO n = 1,atoms%ntype
       DO na = 1,atoms%neq(n)
          natom = natom + 1

          DO  l = 0,atoms%lmax(n)
             !     -----> sum over m
             DO  m = -l,l
                lm = l* (l+1) + m
                !     -----> sum over occupied bands
                DO  i = 1,ne
                   ! coeff. for lz ->
                   orb(l,m,n)%uu = orb(l,m,n)%uu + we(i)*acof(i,lm,natom)* CONJG(acof(i,lm,natom))
                   orb(l,m,n)%dd = orb(l,m,n)%dd + we(i)*bcof(i,lm,natom)* CONJG(bcof(i,lm,natom))
                   ! coeff. for l+ <M'|l+|M> with respect to M ->
                   IF (m.NE.l) THEN
                      orb(l,m,n)%uup = orb(l,m,n)%uup + we(i)*acof(i,lm,natom)* CONJG(acof(i,lm+1,natom))
                      orb(l,m,n)%ddp = orb(l,m,n)%ddp + we(i)*bcof(i,lm,natom)* CONJG(bcof(i,lm+1,natom))
                   ELSE
                      orb(l,m,n)%uup = czero
                      orb(l,m,n)%ddp = czero
                   ENDIF
                   ! coeff. for l- <M'|l-|M> with respect to M ->
                   IF (m.NE.-l) THEN
                      orb(l,m,n)%uum = orb(l,m,n)%uum + we(i)*acof(i,lm,natom)* CONJG(acof(i,lm-1,natom))
                      orb(l,m,n)%ddm = orb(l,m,n)%ddm + we(i)*bcof(i,lm,natom)* CONJG(bcof(i,lm-1,natom))
                   ELSE
                      orb(l,m,n)%uum = czero
                      orb(l,m,n)%ddm = czero
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !
          ! --> Local Orbital contribution: u,lo part
          !
          DO ilo = 1, atoms%nlo(n)
             l = atoms%llo(ilo,n)
             DO m = -l, l
                lm = l* (l+1) + m
                DO i = 1,ne
                   orbl(ilo,m,n)%uulo = orbl(ilo,m,n)%uulo + we(i) * (&
                        acof(i,lm,natom)* CONJG(ccof(m,i,ilo,natom)) +&
                        ccof(m,i,ilo,natom)* CONJG(acof(i,lm,natom)) )
                   orbl(ilo,m,n)%dulo = orbl(ilo,m,n)%dulo + we(i) * (&
                        bcof(i,lm,natom)* CONJG(ccof(m,i,ilo,natom)) +&
                        ccof(m,i,ilo,natom)* CONJG(bcof(i,lm,natom)) )
                   IF (m.NE.l) THEN
                      orbl(ilo,m,n)%uulop = orbl(ilo,m,n)%uulop + we(i) *(&
                           acof(i,lm,natom)* CONJG(ccof(m+1,i,ilo,natom))+&
                           ccof(m,i,ilo,natom)* CONJG(acof(i,lm+1,natom)))
                      orbl(ilo,m,n)%dulop = orbl(ilo,m,n)%dulop + we(i) *(&
                           bcof(i,lm,natom)* CONJG(ccof(m+1,i,ilo,natom))+&
                           ccof(m,i,ilo,natom)* CONJG(bcof(i,lm+1,natom)))
                   ELSE
                      orbl(ilo,m,n)%uulop = czero
                      orbl(ilo,m,n)%dulop = czero
                   ENDIF
                   IF (m.NE.-l) THEN
                      orbl(ilo,m,n)%uulom = orbl(ilo,m,n)%uulom + we(i) *(&
                           acof(i,lm,natom)* CONJG(ccof(m-1,i,ilo,natom))+&
                           ccof(m,i,ilo,natom)* CONJG(acof(i,lm-1,natom)))
                      orbl(ilo,m,n)%dulom = orbl(ilo,m,n)%dulom + we(i) *(&
                           bcof(i,lm,natom)* CONJG(ccof(m-1,i,ilo,natom))+&
                           ccof(m,i,ilo,natom)* CONJG(bcof(i,lm-1,natom)))
                   ELSE
                      orbl(ilo,m,n)%uulom = czero
                      orbl(ilo,m,n)%dulom = czero
                   ENDIF
                ENDDO  ! sum over eigenstates (i)
             ENDDO    ! loop over m
             !
             ! --> lo,lo' part           
             !
             DO ilop = 1, atoms%nlo(n)
                IF (atoms%llo(ilop,n).EQ.l) THEN
                   DO m = -l, l
                      DO i = 1,ne
                         orblo(ilo,ilop,m,n)%z = orblo(ilo,ilop,m,n)%z +&
                              we(i) *   ccof(m,i,ilo, natom) * CONJG( ccof(m,i,ilop,natom) ) 
                         IF (m.NE.l) THEN
                            orblo(ilo,ilop,m,n)%p = orblo(ilo,ilop,m,n)%p +&
                                 we(i) *  ccof(m,  i,ilo, natom) * CONJG( ccof(m+1,i,ilop,natom) ) 
                         ELSE
                            orblo(ilo,ilop,m,n)%p = czero
                         ENDIF
                         IF (m.NE.-l) THEN
                            orblo(ilo,ilop,m,n)%m = orblo(ilo,ilop,m,n)%m +&
                                 we(i) *  ccof(m,  i,ilo, natom) * CONJG( ccof(m-1,i,ilop,natom) )  
                         ELSE
                            orblo(ilo,ilop,m,n)%m = czero
                         ENDIF
                      ENDDO  ! sum over eigenstates (i)
                   ENDDO    ! loop over m
                ENDIF
             ENDDO      ! loop over lo's (ilop)

          ENDDO      ! loop over lo's (ilo)

       ENDDO ! sum over equiv atoms (na)
    ENDDO    ! loop over atom types (n)

    RETURN
  END SUBROUTINE orbmom
END MODULE m_orbmom
