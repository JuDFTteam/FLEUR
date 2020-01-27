MODULE m_orbmom
  !     ***************************************************************
  !     perform the sum over m (for each l) and bands to set up the
  !     coefficient of spherical contribution to orbital moment.
  !     all quantities are in the local spin-frame
  !     ***************************************************************

CONTAINS
  SUBROUTINE orbmom(atoms,ne,we,ispin,eigVecCoeffs,orb)

    !USE m_types, ONLY : t_orb,t_orbl,t_orblo
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),        INTENT(IN) :: atoms
    TYPE(t_eigVecCoeffs), INTENT(IN) :: eigVecCoeffs
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne, ispin
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(:)!(nobd)
    TYPE (t_orb), INTENT (INOUT) :: orb

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
                   orb%uu(l,m,n,ispin) = orb%uu(l,m,n,ispin) + we(i)*eigVecCoeffs%acof(i,lm,natom,ispin)*&
                                                               CONJG(eigVecCoeffs%acof(i,lm,natom,ispin))
                   orb%dd(l,m,n,ispin) = orb%dd(l,m,n,ispin) + we(i)*eigVecCoeffs%bcof(i,lm,natom,ispin)*&
                                                               CONJG(eigVecCoeffs%bcof(i,lm,natom,ispin))
                   ! coeff. for l+ <M'|l+|M> with respect to M ->
                   IF (m.NE.l) THEN
                      orb%uup(l,m,n,ispin) = orb%uup(l,m,n,ispin) + we(i)*eigVecCoeffs%acof(i,lm,natom,ispin)*&
                                                                    CONJG(eigVecCoeffs%acof(i,lm+1,natom,ispin))
                      orb%ddp(l,m,n,ispin) = orb%ddp(l,m,n,ispin) + we(i)*eigVecCoeffs%bcof(i,lm,natom,ispin)*&
                                                                    CONJG(eigVecCoeffs%bcof(i,lm+1,natom,ispin))
                   ELSE
                      orb%uup(l,m,n,ispin) = czero
                      orb%ddp(l,m,n,ispin) = czero
                   ENDIF
                   ! coeff. for l- <M'|l-|M> with respect to M ->
                   IF (m.NE.-l) THEN
                      orb%uum(l,m,n,ispin) = orb%uum(l,m,n,ispin) + we(i)*eigVecCoeffs%acof(i,lm,natom,ispin)*&
                                                                    CONJG(eigVecCoeffs%acof(i,lm-1,natom,ispin))
                      orb%ddm(l,m,n,ispin) = orb%ddm(l,m,n,ispin) + we(i)*eigVecCoeffs%bcof(i,lm,natom,ispin)*&
                                                                    CONJG(eigVecCoeffs%bcof(i,lm-1,natom,ispin))
                   ELSE
                      orb%uum(l,m,n,ispin) = czero
                      orb%ddm(l,m,n,ispin) = czero
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
                   orb%uulo(ilo,m,n,ispin) = orb%uulo(ilo,m,n,ispin) + we(i) * (&
                        eigVecCoeffs%acof(i,lm,natom,ispin)* CONJG(eigVecCoeffs%ccof(m,i,ilo,natom,ispin)) +&
                        eigVecCoeffs%ccof(m,i,ilo,natom,ispin)* CONJG(eigVecCoeffs%acof(i,lm,natom,ispin)) )
                   orb%dulo(ilo,m,n,ispin) = orb%dulo(ilo,m,n,ispin) + we(i) * (&
                        eigVecCoeffs%bcof(i,lm,natom,ispin)* CONJG(eigVecCoeffs%ccof(m,i,ilo,natom,ispin)) +&
                        eigVecCoeffs%ccof(m,i,ilo,natom,ispin)* CONJG(eigVecCoeffs%bcof(i,lm,natom,ispin)) )
                   IF (m.NE.l) THEN
                      orb%uulop(ilo,m,n,ispin) = orb%uulop(ilo,m,n,ispin) + we(i) *(&
                           eigVecCoeffs%acof(i,lm,natom,ispin)* CONJG(eigVecCoeffs%ccof(m+1,i,ilo,natom,ispin))+&
                           eigVecCoeffs%ccof(m,i,ilo,natom,ispin)* CONJG(eigVecCoeffs%acof(i,lm+1,natom,ispin)))
                      orb%dulop(ilo,m,n,ispin) = orb%dulop(ilo,m,n,ispin) + we(i) *(&
                           eigVecCoeffs%bcof(i,lm,natom,ispin)* CONJG(eigVecCoeffs%ccof(m+1,i,ilo,natom,ispin))+&
                           eigVecCoeffs%ccof(m,i,ilo,natom,ispin)* CONJG(eigVecCoeffs%bcof(i,lm+1,natom,ispin)))
                   ELSE
                      orb%uulop(ilo,m,n,ispin) = czero
                      orb%dulop(ilo,m,n,ispin) = czero
                   ENDIF
                   IF (m.NE.-l) THEN
                      orb%uulom(ilo,m,n,ispin) = orb%uulom(ilo,m,n,ispin) + we(i) *(&
                           eigVecCoeffs%acof(i,lm,natom,ispin)* CONJG(eigVecCoeffs%ccof(m-1,i,ilo,natom,ispin))+&
                           eigVecCoeffs%ccof(m,i,ilo,natom,ispin)* CONJG(eigVecCoeffs%acof(i,lm-1,natom,ispin)))
                      orb%dulom(ilo,m,n,ispin) = orb%dulom(ilo,m,n,ispin) + we(i) *(&
                           eigVecCoeffs%bcof(i,lm,natom,ispin)* CONJG(eigVecCoeffs%ccof(m-1,i,ilo,natom,ispin))+&
                           eigVecCoeffs%ccof(m,i,ilo,natom,ispin)* CONJG(eigVecCoeffs%bcof(i,lm-1,natom,ispin)))
                   ELSE
                      orb%uulom(ilo,m,n,ispin) = czero
                      orb%dulom(ilo,m,n,ispin) = czero
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
                         orb%z(ilo,ilop,m,n,ispin) = orb%z(ilo,ilop,m,n,ispin) +&
                              we(i) *   eigVecCoeffs%ccof(m,i,ilo,natom,ispin) * CONJG( eigVecCoeffs%ccof(m,i,ilop,natom,ispin) ) 
                         IF (m.NE.l) THEN
                            orb%p(ilo,ilop,m,n,ispin) = orb%p(ilo,ilop,m,n,ispin) +&
                                 we(i) *  eigVecCoeffs%ccof(m,i,ilo,natom,ispin) * CONJG( eigVecCoeffs%ccof(m+1,i,ilop,natom,ispin) ) 
                         ELSE
                            orb%p(ilo,ilop,m,n,ispin) = czero
                         ENDIF
                         IF (m.NE.-l) THEN
                            orb%m(ilo,ilop,m,n,ispin) = orb%m(ilo,ilop,m,n,ispin) +&
                                 we(i) *  eigVecCoeffs%ccof(m,i,ilo,natom,ispin) * CONJG( eigVecCoeffs%ccof(m-1,i,ilop,natom,ispin) )  
                         ELSE
                            orb%m(ilo,ilop,m,n,ispin) = czero
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
