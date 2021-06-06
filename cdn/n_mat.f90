!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_nmat
   !     ************************************************************
   !     This subroutine calculates the density matrix n^{s}_{m,m'}
   !     for a given atom 'n' and l-quantum number 'l'. The l's for
   !     all atoms are stored in lda_u(), if lda_u()<0, no +U is used.
   !     For details see Eq.(12) of Shick et al. PRB 60, 10765 (1999)
   !     Part of the LDA+U package                   G.B., Oct. 2000
   !     Extension to multiple U per atom type by G.M. 2017
   !     ************************************************************
   CONTAINS
   SUBROUTINE n_mat(atoms,sym,ne,usdus,jspin,we,eigVecCoeffs,den)

      USE m_types
      USE m_constants
      USE m_symMMPmat

      IMPLICIT NONE

      TYPE(t_usdus),       INTENT(IN)     :: usdus
      TYPE(t_sym),         INTENT(IN)     :: sym
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_eigVecCoeffs),INTENT(IN)     :: eigVecCoeffs
      INTEGER,             INTENT(IN)     :: ne,jspin
      REAL,                INTENT(IN)     :: we(:)!(input%neig)
      TYPE(t_potden),      INTENT(INOUT)  :: den

      INTEGER iBand,l,m,lp,mp,n,natom,i_u
      INTEGER ilo,ilop,ll1,lmp,lm

      COMPLEX uu(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX dd(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX du(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      !
      ! calculate n_mat:
      !
      DO i_u = 1,atoms%n_u
         n = atoms%lda_u(i_u)%atomType
         l = atoms%lda_u(i_u)%l
         ll1 = (l+1)*l
         DO natom = SUM(atoms%neq(:n-1))+1, SUM(atoms%neq(:n))
            uu = cmplx_0
            dd = cmplx_0
            du = cmplx_0
            !
            !  prepare n_mat in local frame (in noco-calculations this depends
            !                                also on alpha(n) and beta(n) )
            !
            DO m = -l,l
               lm = ll1+m
               DO mp = -l,l
                  lmp = ll1+mp
                  DO iBand = 1,ne
                     uu(m,mp) = uu(m,mp) + we(iBand) * conjg(eigVecCoeffs%acof(iBand,lmp,natom,jspin))*eigVecCoeffs%acof(iBand,lm,natom,jspin)
                     dd(m,mp) = dd(m,mp) + we(iBand) * conjg(eigVecCoeffs%bcof(iBand,lmp,natom,jspin))*eigVecCoeffs%bcof(iBand,lm,natom,jspin)
                     du(m,mp) = du(m,mp) + we(iBand) * conjg(eigVecCoeffs%bcof(iBand,lmp,natom,jspin))*eigVecCoeffs%acof(iBand,lm,natom,jspin)
                  ENDDO
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
            !                c_0 = c_0 +  we(i) * ( usdus%uulon(ilo,n,jspin) * (&
            !                            conjg(eigVecCoeffs%acof(i,lmp,natom,jspin))*eigVecCoeffs%ccof(m,i,ilo,natom,jspin) &
            !                          + conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,jspin))*eigVecCoeffs%acof(i,lm,natom,jspin) )&
            !                          + usdus%dulon(ilo,n,jspin) * (&
            !                            conjg(eigVecCoeffs%bcof(i,lmp,natom,jspin))*eigVecCoeffs%ccof(m,i,ilo,natom,jspin) &
            !                          + conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,jspin))*eigVecCoeffs%bcof(i,lm,natom,jspin)))
            !             ENDDO
            !             DO ilop = 1, atoms%nlo(n)
            !                IF (atoms%llo(ilop,n).EQ.l) THEN
            !                   DO i = 1,ne
            !                      c_0 = c_0 +  we(i) * usdus%uloulopn(ilo,ilop,n,jspin) *&
            !                                  conjg(eigVecCoeffs%ccof(mp,i,ilop,natom,jspin)) *eigVecCoeffs%ccof(m,i,ilo,natom,jspin)
            !                   ENDDO
            !                ENDIF
            !             ENDDO
            !             n_tmp(m,mp) = n_tmp(m,mp) + c_0
            !          ENDDO
            !       ENDDO
            !    ENDIF
            ! ENDDO
            !
            !  n_mmp should be rotated by D_mm' ; compare force_a21
            !
            den%mmpMat_uu(:,:,i_u,jspin) = den%mmpMat_uu(:,:,i_u,jspin) + symMMPmat(uu,sym,natom,l) * 1.0/atoms%neq(n)
            den%mmpMat_dd(:,:,i_u,jspin) = den%mmpMat_dd(:,:,i_u,jspin) + symMMPmat(dd,sym,natom,l) * 1.0/atoms%neq(n)
            den%mmpMat_du(:,:,i_u,jspin) = den%mmpMat_du(:,:,i_u,jspin) + symMMPmat(du,sym,natom,l) * 1.0/atoms%neq(n)
            den%mmpMat_ud(:,:,i_u,jspin) = den%mmpMat_ud(:,:,i_u,jspin) + symMMPmat(du,sym,natom,l) * 1.0/atoms%neq(n)
         ENDDO ! sum  over equivalent atoms
      END DO !loop over u parameters

   END SUBROUTINE n_mat
END MODULE m_nmat
