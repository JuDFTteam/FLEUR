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
   SUBROUTINE n_mat(atoms,radfun,sym,ne,we,abc,abc1,n_mmp,ntype,jsp,jsp1)
      USE m_types_radfun
      USE m_types_abc
      USE m_types
      USE m_constants
      USE m_symMMPmat

      IMPLICIT NONE
      TYPE(t_sym),         INTENT(IN)     :: sym
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_radfun),      INTENT(IN)     :: radfun
      TYPE(t_abc)         ,INTENT(IN)     :: abc,abc1
      INTEGER,             INTENT(IN)     :: ne,ntype,jsp,jsp1
      REAL,                INTENT(IN)     :: we(:)!(input%neig)
      COMPLEX,             INTENT(INOUT)  :: n_mmp(-lmaxU_const:,-lmaxU_const:,:)

      INTEGER i,l,m,lp,mp,n,natom,i_u,i_denmat
      INTEGER ll1,lmp,lm,j,jj
      COMPLEX c_0

      COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      !
      ! calculate n_mat:
      !PRINT *,'Hello'
      !WRITE (*,*) 'Hello'
      DO i_u = 1,atoms%n_u+atoms%n_opc
         if(i_u>atoms%n_u) then
            i_denmat = i_u + atoms%n_hia
            n = atoms%lda_opc(i_u-atoms%n_u)%atomType
            l = atoms%lda_opc(i_u-atoms%n_u)%l
         else
            i_denmat = i_u
            n = atoms%lda_u(i_u)%atomType
            l = atoms%lda_u(i_u)%l
         endif
         if (n/=ntype) cycle !Only for atom types we currently have abc coefficients for
         ll1 = (l+1)*l
         DO natom = 1, atoms%neq(n) 
            n_tmp = cmplx_0
            !
            !  prepare n_mat in local frame (in noco-calculations this depends
            !                                also on alpha(n) and beta(n) )
            !
            DO m = -l,l
               lm = ll1+m
               DO mp = -l,l
                  lmp = ll1+mp
                  c_0 = 0.0
                  DO i = 1,ne
                     DO j=1,size(abc%cof,3)
                        DO jj=1,size(abc1%cof,3)
                        c_0 = c_0 +  we(i) *  conjg(abc%cof(i,lmp,j,natom))*abc%cof(i,lm,jj,natom)*radfun%integral(j, jj, l, jsp, jsp1)
                        ENDDO
                     ENDDO
                  ENDDO      
                  n_tmp(m,mp) = c_0
               ENDDO
            ENDDO
            !
            !  n_mmp should be rotated by D_mm' ; compare force_a21
            !
            n_mmp(:,:,i_denmat) = n_mmp(:,:,i_denmat) + symMMPmat(conjg(n_tmp),sym,natom+atoms%firstatom(n)-1,l) * 1.0/atoms%neq(n)
         ENDDO ! sum  over equivalent atoms
      END DO !loop over u parameters

   END SUBROUTINE n_mat
END MODULE m_nmat
