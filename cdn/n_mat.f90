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
  SUBROUTINE n_mat(atoms,sym, ne,usdus,jspin,we,eigVecCoeffs,n_mmp)
    !

    USE m_types
    USE m_constants
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(IN)        :: usdus
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_eigVecCoeffs),INTENT(IN) :: eigVecCoeffs
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne,jspin 
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(:)!(input%neig)
    COMPLEX, INTENT (INOUT) :: n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_u)
    !     ..
    !     .. Local Scalars ..
    COMPLEX c_0
    INTEGER i,j,k,l ,mp,n,it,is,isi,natom,natomTemp,n_ldau,lp,m,i_u
    INTEGER ilo,ilop,ll1,nn,lmp,lm
    REAL fac
    !     ..
    !     .. Local Arrays ..
    COMPLEX n_tmp(-3:3,-3:3),nr_tmp(-3:3,-3:3),d_tmp(-3:3,-3:3)
    COMPLEX n1_tmp(-3:3,-3:3)
    !     ..
    !
    ! calculate n_mat:
    !
    natom = 0
    i_u = 1
    DO n = 1,atoms%ntype
       DO WHILE (i_u.LE.atoms%n_u)
          IF (atoms%lda_u(i_u)%atomType.GT.n) EXIT
          natomTemp = natom
          n_tmp(:,:) = cmplx(0.0,0.0)
          l = atoms%lda_u(i_u)%l
          ll1 = (l+1)*l 
          DO nn = 1, atoms%neq(n)
             natomTemp = natomTemp + 1
             !
             !  prepare n_mat in local frame (in noco-calculations this depends 
             !                                also on alpha(n) and beta(n) )
             !
             DO m = -l,l
                lm = ll1+m
                DO mp = -l,l
                   lmp = ll1+mp
                   c_0 = cmplx(0.0,0.0)
                   DO i = 1,ne
                      c_0 = c_0 +  we(i) * ( usdus%ddn(l,n,jspin) *&
                           conjg(eigVecCoeffs%bcof(i,lmp,natomTemp,jspin))*eigVecCoeffs%bcof(i,lm,natomTemp,jspin) +&
                           conjg(eigVecCoeffs%acof(i,lmp,natomTemp,jspin))*eigVecCoeffs%acof(i,lm,natomTemp,jspin) )
                   ENDDO
                   n_tmp(m,mp) = c_0 
                ENDDO
             ENDDO
             !
             !  add local orbrbital contribution (if there is one) (untested so far)
             !
             DO ilo = 1, atoms%nlo(n)
                IF (atoms%llo(ilo,n).EQ.l) THEN

                   DO m = -l,l
                      lm = ll1+m
                      DO mp = -l,l
                         lmp = ll1+mp
                         c_0 = cmplx(0.0,0.0)
                         DO i = 1,ne
                            c_0 = c_0 +  we(i) * (  usdus%uulon(ilo,n,jspin) * (&
                                 conjg(eigVecCoeffs%acof(i,lmp,natomTemp,jspin))*eigVecCoeffs%ccof(m,i,ilo,natomTemp,jspin) +&
                                 conjg(eigVecCoeffs%ccof(mp,i,ilo,natomTemp,jspin))*eigVecCoeffs%acof(i,lm,natomTemp,jspin) )&
                                 + usdus%dulon(ilo,n,jspin) * (&
                                 conjg(eigVecCoeffs%bcof(i,lmp,natomTemp,jspin))*eigVecCoeffs%ccof(m,i,ilo,natomTemp,jspin) +&
                                 conjg(eigVecCoeffs%ccof(mp,i,ilo,natomTemp,jspin))*eigVecCoeffs%bcof(i,lm,natomTemp,jspin)))
                         ENDDO
                         DO ilop = 1, atoms%nlo(n)
                            IF (atoms%llo(ilop,n).EQ.l) THEN
                               DO i = 1,ne
                                  c_0 = c_0 +  we(i) * usdus%uloulopn(ilo,ilop,n,jspin) *&
                                       conjg(eigVecCoeffs%ccof(mp,i,ilop,natomTemp,jspin)) *eigVecCoeffs%ccof(m,i,ilo,natomTemp,jspin)
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
             DO it = 1, sym%invarind(natomTemp)

                fac = 1.0  /  ( sym%invarind(natomTemp) * atoms%neq(n) )
                is = sym%invarop(natomTemp,it)
                isi = sym%invtab(is)
                d_tmp(:,:) = cmplx(0.0,0.0)
                DO m = -l,l
                   DO mp = -l,l
                      d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                   ENDDO
                ENDDO
                nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                n1_tmp =  matmul( nr_tmp, d_tmp )
                DO m = -l,l
                   DO mp = -l,l
                      n_mmp(m,mp,i_u) = n_mmp(m,mp,i_u) + conjg(n1_tmp(m,mp)) * fac
                   ENDDO
                ENDDO

             ENDDO

          ENDDO ! sum  over equivalent atoms
          i_u = i_u + 1
       END DO
       natom = natom + atoms%neq(n)
    ENDDO     ! loop over atom types

    !     do m=-l,l
    !      write(*,'(14f12.6)') (n_mmp(m,mp),mp=-l,l)
    !     enddo
    !
    RETURN
  END SUBROUTINE n_mat
END MODULE m_nmat
