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
   SUBROUTINE n_mat(atoms,sym,ne,usdus,jspin,we,eigVecCoeffs,n_mmp)

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
      COMPLEX,             INTENT(INOUT)  :: n_mmp(-lmaxU_const:,-lmaxU_const:,:)

      INTEGER i,l,m,lp,mp,n,natom,i_u
      INTEGER ilo,ilop,ll1,lmp,lm
      COMPLEX c_0

      COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      !
      ! calculate n_mat:
      !
      DO i_u = 1,atoms%n_u
         n = atoms%lda_u(i_u)%atomType
         l = atoms%lda_u(i_u)%l
         ll1 = (l+1)*l
         DO natom = SUM(atoms%neq(:n-1))+1, SUM(atoms%neq(:n))
            n_tmp = cmplx_0
            !
            !  prepare n_mat in local frame (in noco-calculations this depends
            !                                also on alpha(n) and beta(n) )
            !
            DO m = -l,l
               lm = ll1+m
               DO mp = -l,l
                  lmp = ll1+mp
                  c_0 = cmplx_0
                  DO i = 1,ne
                     c_0 = c_0 +  we(i) * ( usdus%ddn(l,n,jspin) *&
                                 conjg(eigVecCoeffs%bcof(i,lmp,natom,jspin))*eigVecCoeffs%bcof(i,lm,natom,jspin) &
                               + conjg(eigVecCoeffs%acof(i,lmp,natom,jspin))*eigVecCoeffs%acof(i,lm,natom,jspin) )
                  ENDDO
                  n_tmp(m,mp) = c_0
               ENDDO
            ENDDO
            !
            !  add local orbital contribution (if there is one) (untested so far)
            !
            DO ilo = 1, atoms%nlo(n)
               IF (atoms%llo(ilo,n).EQ.l) THEN
                  DO m = -l,l
                     lm = ll1+m
                     DO mp = -l,l
                        lmp = ll1+mp
                        c_0 = cmplx_0
                        DO i = 1,ne
                           c_0 = c_0 +  we(i) * ( usdus%uulon(ilo,n,jspin) * (&
                                       conjg(eigVecCoeffs%acof(i,lmp,natom,jspin))*eigVecCoeffs%ccof(m,i,ilo,natom,jspin) &
                                     + conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,jspin))*eigVecCoeffs%acof(i,lm,natom,jspin) )&
                                     + usdus%dulon(ilo,n,jspin) * (&
                                       conjg(eigVecCoeffs%bcof(i,lmp,natom,jspin))*eigVecCoeffs%ccof(m,i,ilo,natom,jspin) &
                                     + conjg(eigVecCoeffs%ccof(mp,i,ilo,natom,jspin))*eigVecCoeffs%bcof(i,lm,natom,jspin)))
                        ENDDO
                        DO ilop = 1, atoms%nlo(n)
                           IF (atoms%llo(ilop,n).EQ.l) THEN
                              DO i = 1,ne
                                 c_0 = c_0 +  we(i) * usdus%uloulopn(ilo,ilop,n,jspin) *&
                                             conjg(eigVecCoeffs%ccof(mp,i,ilop,natom,jspin)) *eigVecCoeffs%ccof(m,i,ilo,natom,jspin)
                              ENDDO
                           ENDIF
                        ENDDO
                        n_tmp(m,mp) = n_tmp(m,mp) + c_0
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            !
            !  n_mmp should be rotated by D_mm' ; compare force_a21
            !
            n_mmp(:,:,i_u) = n_mmp(:,:,i_u) + conjg(symMMPmat(n_tmp,sym,natom,l)) * 1.0/atoms%neq(n)
         ENDDO ! sum  over equivalent atoms
      END DO !loop over u parameters

   END SUBROUTINE n_mat
END MODULE m_nmat
