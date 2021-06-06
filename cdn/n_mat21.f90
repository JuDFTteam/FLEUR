MODULE m_nmat21
   !------------------------------------------------------------
   !This subroutine calculates the density matrix n^{s}_{m,m'}
   !for a given atom 'n' and l-quantum number 'l'. The l's for
   !all atoms are stored in lda_u(), if lda_u()<0, no +U is used.
   !For details see Eq.(12) of Shick et al. PRB 60, 10765 (1999)
   !Part of the LDA+U package                   G.B., Oct. 2000
   !------------------------------------------------------------

   USE m_types
   USE m_constants
   USE m_symMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE n_mat21(atoms,sym,ne,we,denCoeffsOffdiag,eigVecCoeffs,den)

      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffdiag
      INTEGER,                   INTENT(IN)     :: ne
      REAL,                      INTENT(IN)     :: we(:)!(input%neig)
      TYPE(t_potden),            INTENT(INOUT)  :: den

      INTEGER iBand,l,m,lp,mp,n,natom
      INTEGER ilo,ilop,ll1,lmp,lm,i_u
      COMPLEX c_0

      COMPLEX uu(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX dd(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX du(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX ud(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      !
      ! calculate n_mat:
      !
      DO i_u = 1,atoms%n_u
         l = atoms%lda_u(i_u)%l
         n = atoms%lda_u(i_u)%atomType
         ll1 = (l+1)*l
         DO natom = SUM(atoms%neq(:n-1))+1, SUM(atoms%neq(:n))
            uu = cmplx_0
            dd = cmplx_0
            du = cmplx_0
            ud = cmplx_0
            !
            !  prepare n_mat in local frame (in noco-calculations this depends
            !                                also on alpha(n) and beta(n) )
            !
            DO m = -l,l
               lm = ll1+m
               DO mp = -l,l
                  lmp = ll1+mp
                  DO iBand = 1,ne
                     uu(m,mp) = uu(m,mp) + we(iBand) * conjg(eigVecCoeffs%acof(iBand,lmp,natom,2))*eigVecCoeffs%acof(iBand,lm,natom,1)
                     dd(m,mp) = dd(m,mp) + we(iBand) * conjg(eigVecCoeffs%bcof(iBand,lmp,natom,2))*eigVecCoeffs%bcof(iBand,lm,natom,1)
                     ud(m,mp) = ud(m,mp) + we(iBand) * conjg(eigVecCoeffs%acof(iBand,lmp,natom,2))*eigVecCoeffs%bcof(iBand,lm,natom,1)
                     du(m,mp) = du(m,mp) + we(iBand) * conjg(eigVecCoeffs%bcof(iBand,lmp,natom,2))*eigVecCoeffs%acof(iBand,lm,natom,1)
                  ENDDO
                  uu(m,mp) = CMPLX(-REAL(uu(m,mp)), AIMAG(uu(m,mp)))
                  dd(m,mp) = CMPLX(-REAL(dd(m,mp)), AIMAG(dd(m,mp)))
                  du(m,mp) = CMPLX(-REAL(du(m,mp)), AIMAG(du(m,mp)))
                  ud(m,mp) = CMPLX(-REAL(ud(m,mp)), AIMAG(ud(m,mp)))
               ENDDO
            ENDDO
            !
            !  add local orbital contribution (if there is one) (untested so far)
            !
            ! DO ilo = 1, atoms%nlo(n)
            !    IF (atoms%llo(ilo,n).EQ.l) THEN
            !       DO m = -l,l
            !          lm = ll1+m
            !          DO mp = -l,l
            !             lmp = ll1+mp
            !             c_0 = cmplx_0
            !             DO i = 1,ne
            !                c_0 = c_0 +  we(i) * ( &
            !                            conjg(eigVecCoeffs%acof(i,lmp,natom,2))*eigVecCoeffs%ccof(m,i,ilo,natom,1) * denCoeffsOffdiag%uulo21n(l,n) &
            !                          + conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,2))*eigVecCoeffs%acof(i,lm,natom,1) * denCoeffsOffdiag%ulou21n(l,n) &
            !                          + conjg(eigVecCoeffs%bcof(i,lmp,natom,2))*eigVecCoeffs%ccof(m,i,ilo,natom,1) * denCoeffsOffdiag%dulo21n(l,n) &
            !                          + conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,2))*eigVecCoeffs%bcof(i,lm,natom,1) * denCoeffsOffdiag%ulod21n(l,n))
            !             ENDDO
            !             DO ilop = 1, atoms%nlo(n)
            !                IF (atoms%llo(ilop,n).EQ.l) THEN
            !                   DO i = 1,ne
            !                      c_0 = c_0 +  we(i) * denCoeffsOffdiag%uloulop21n(ilo,ilop,n) &
            !                                 *conjg(eigVecCoeffs%ccof(mp,i,ilop,natom,2)) *eigVecCoeffs%ccof(m ,i,ilo ,natom,1)
            !                   ENDDO
            !                ENDIF
            !             ENDDO
            !             n_tmp(m,mp) = n_tmp(m,mp) + CMPLX(-REAL(c_0), AIMAG(c_0))
            !          ENDDO
            !       ENDDO
            !    ENDIF
            ! ENDDO
            !
            !  n_mmp should be rotated by D_mm' ; compare force_a21
            !
            !Note: This can be done only if the correct magnetic symmetries are
            !present. This is not the case at the moment (Jan 2020)
            den%mmpMat_uu(:,:,i_u,3) = den%mmpMat_uu(:,:,i_u,3) + symMMPmat(uu,sym,natom,l,phase=.TRUE.) * 1.0/atoms%neq(n)
            den%mmpMat_dd(:,:,i_u,3) = den%mmpMat_dd(:,:,i_u,3) + symMMPmat(dd,sym,natom,l,phase=.TRUE.) * 1.0/atoms%neq(n)
            den%mmpMat_du(:,:,i_u,3) = den%mmpMat_du(:,:,i_u,3) + symMMPmat(du,sym,natom,l,phase=.TRUE.) * 1.0/atoms%neq(n)
            den%mmpMat_ud(:,:,i_u,3) = den%mmpMat_ud(:,:,i_u,3) + symMMPmat(du,sym,natom,l,phase=.TRUE.) * 1.0/atoms%neq(n)
         ENDDO ! sum  over equivalent
      ENDDO     ! loop over u parameters

   END SUBROUTINE n_mat21
END MODULE m_nmat21
