!-------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!-------------------------------------------------------------------------------

MODULE m_dfpt_rhomt
   !! Module adapted from rhomt21, which in turn was adapted from rhomt, and
   !! from rhomtlo for the lo part.
   !!
   !! It contains subroutines for the non-LO and LO spherical contributions
   !! to the density coefficients used to construct \(\rho_{0}^{\alpha}(r)\).
CONTAINS
   SUBROUTINE dfpt_rhomt(atoms,we,we1,ne,ilSpinPr,ilSpin,qpoint,l_dfpt,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
      !! Subroutine to construct all spherical MT density coefficients (for a
      !! density perturbation) without LOs in one routine. The spin input dictates,
      !! which element is gonna be built.
      !! The coefficients are of the form:
      !! \begin{aligned}
      !! d_{l,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha} &= \sum_{\nu\boldsymbol{k}}\sum_{m}\\
      !! &* A_{l,m,\lambda'}^{\sigma_{\alpha}',\nu\boldsymbol{k}*}\\
      !! &* (2\tilde{f}_{\nu\boldsymbol{k}}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{q},j,\beta~(1)}+\tilde{f}_{\nu\boldsymbol{k}\boldsymbol{q}}^{j,\beta~(1)}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}})
      !! \end{aligned}
      !! The k-point loop is performed outside this routine. In contrast to older
      !! routines, the arrays uu etc. and uu21 etc. are merged into one
      !! spin and u-order dependent array.
      !!
      !! \(\sigma_{\alpha}(')\): local spin indices \(\rightarrow\) ilSpinPr, ilSpin
      !!
      !! \(\lambda(')\): order of the radial basis function (0: u, 1: d)
      !!
      !! \(\nu\boldsymbol{k}\): State index (k-point and number of state)
      !!
      !! \(\boldsymbol{q},j,\beta\): Perturbation index; q-point, direction and atom
      !!
      !! \(\tilde{f}_{\nu\boldsymbol{k}}\): (Smeared) occupation number [perturbed for \(\tilde{f}^{(1)}\)]
      !!
      !! \(A\): Summed matching coefficients and eigenvectors [perturbed for \(A^{(1)}\)]

      USE m_types_atoms
      USE m_types_cdnval

      IMPLICIT NONE

      TYPE(t_atoms),       INTENT(IN)    :: atoms

      INTEGER,              INTENT(IN)    :: ne
      INTEGER,              INTENT(IN)    :: ilSpinPr  !! \(\sigma_{\alpha}^{'}\)
      INTEGER,              INTENT(IN)    :: ilSpin    !! \(\sigma_{\alpha}\)
      REAL,                 INTENT(IN)    :: we(ne)    !! \(\tilde{f}_{\nu\boldsymbol{k}}\)
      REAL,                 INTENT(IN)    :: we1(ne)   !! \(\tilde{f}_{\nu\boldsymbol{k}\boldsymbol{q}}^{j,\beta~(1)}\)
      REAL,                 INTENT(IN)    :: qpoint(3) !! \(\boldsymbol{q}\)
      LOGICAL,              INTENT(IN)    :: l_dfpt

      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs  !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}}\)
      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs1 !! \(A_{l,m,\lambda}^{\sigma_{\alpha},\nu\boldsymbol{k}\boldsymbol{q},j,\beta~(1)}\)

      TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs !! \(d_{l,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}\)

      INTEGER i,l,lm,itype,na,natom,lo,lop,m

      COMPLEX :: temp

      natom = 0
      DO itype = 1,atoms%ntype
         DO na = 1,atoms%neq(itype)
            natom = natom + 1
            DO l = 0,atoms%lmax(itype)
               DO m = -l,l
                  lm = l* (l+1) + m
                  DO i = 1,ne
                     ! uu/du
                     temp = we(i) * eigVecCoeffs1%abcof(i,lm,0,natom,ilSpin) ! If not DFPT, this is the base case for rhomt(21)
                     IF (l_dfpt) THEN
                        temp = temp * 2.0
                        IF (norm2(qpoint)<=1e-8) THEN
                           temp = temp + we1(i) * eigVecCoeffs%abcof(i,lm,0,natom,ilSpin)
                        END IF
                     END IF
                     denCoeffs%mt_coeff(l,itype,0,0,ilSpinPr, ilSpin) = denCoeffs%mt_coeff(l,itype,0,0,ilSpinPr, ilSpin) &
                                                                    & + we(i) * CONJG(eigVecCoeffs%abcof(i,lm,0,natom,ilSpinPr)) * temp
                     denCoeffs%mt_coeff(l,itype,1,0,ilSpinPr, ilSpin) = denCoeffs%mt_coeff(l,itype,1,0,ilSpinPr, ilSpin) &
                                                                    & + we(i) * CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ilSpinPr)) * temp
                     ! ud/dd
                     temp = we(i) * eigVecCoeffs1%abcof(i,lm,1,natom,ilSpin)
                     IF (l_dfpt) THEN
                        temp = temp * 2.0
                        IF (norm2(qpoint)<=1e-8) THEN
                           temp = temp + we1(i) * eigVecCoeffs%abcof(i,lm,1,natom,ilSpin)
                        END IF
                     END IF
                     denCoeffs%mt_coeff(l,itype,0,1,ilSpinPr, ilSpin) = denCoeffs%mt_coeff(l,itype,0,1,ilSpinPr, ilSpin) &
                                                                    & + we(i) * CONJG(eigVecCoeffs%abcof(i,lm,0,natom,ilSpinPr)) * temp
                     denCoeffs%mt_coeff(l,itype,1,1,ilSpinPr, ilSpin) = denCoeffs%mt_coeff(l,itype,1,1,ilSpinPr, ilSpin) &
                                                                    & + we(i) * CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ilSpinPr)) * temp
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE dfpt_rhomt

   SUBROUTINE dfpt_rhomtlo(atoms,ne,we,we1,ilSpinPr,ilSpin,qpoint,l_dfpt,eigVecCoeffs,eigVecCoeffs1,denCoeffs)
      !! This is a complementary routine to the one above for \(\lambda(')\ge 2\),
      !! i.e. mixed or pure LO contributions.

      USE m_types_atoms
      USE m_types_cdnval

      IMPLICIT NONE

      TYPE(t_atoms),        INTENT(IN)    :: atoms

      INTEGER, INTENT(IN) :: ne, ilSpinPr, ilSpin

      REAL,    INTENT(IN) :: we(:),we1(:)!(nobd)
      REAL,    INTENT(IN) :: qpoint(3)
      LOGICAL, INTENT(IN) :: l_dfpt

      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs, eigVecCoeffs1
      TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

      INTEGER i,l,lm,lo,lop ,natom,nn,ntyp,m

      COMPLEX :: temp

      natom = 0
      DO ntyp = 1,atoms%ntype
         DO nn = 1,atoms%neq(ntyp)
            natom = natom + 1
            DO lo = 1,atoms%nlo(ntyp)
               l = atoms%llo(lo,ntyp)
               DO m = -l,l
                  lm = l* (l+1) + m
                  DO i = 1,ne
                     temp = we(i) * eigVecCoeffs1%ccof(m,i,lo,natom,ilSpin) ! If not DFPT, this is the base case for rhomt(21)
                     IF (l_dfpt) THEN
                        temp = temp * 2.0
                        IF (norm2(qpoint)<=1e-8) THEN
                           temp = temp + we1(i) * eigVecCoeffs%ccof(m,i,lo,natom,ilSpin)
                        END IF
                     END IF
                     denCoeffs%mt_ulo_coeff(lo,ntyp,0,ilSpinPr,ilSpin) = denCoeffs%mt_ulo_coeff(lo,ntyp,0,ilSpinPr,ilSpin) &
                                                                     & + we(i) * CONJG(eigVecCoeffs%abcof(i,lm,0,natom,ilSpinPr)) * temp
                     denCoeffs%mt_ulo_coeff(lo,ntyp,1,ilSpinPr,ilSpin) = denCoeffs%mt_ulo_coeff(lo,ntyp,1,ilSpinPr,ilSpin) &
                                                                     & + we(i) * CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ilSpinPr)) * temp
                     temp = we(i) * eigVecCoeffs1%abcof(i,lm,0,natom,ilSpin)
                     IF (l_dfpt) THEN
                        temp = temp * 2.0
                        IF (norm2(qpoint)<=1e-8) THEN
                           temp = temp + we1(i) * eigVecCoeffs%abcof(i,lm,0,natom,ilSpin)
                        END IF
                     END IF
                     denCoeffs%mt_lou_coeff(lo,ntyp,0,ilSpinPr,ilSpin) = denCoeffs%mt_lou_coeff(lo,ntyp,0,ilSpinPr,ilSpin) &
                                                                     & + we(i) * CONJG(eigVecCoeffs%ccof(m,i,lo,natom,ilSpinPr)) * temp
                     temp = we(i) * eigVecCoeffs1%abcof(i,lm,1,natom,ilSpin)
                     IF (l_dfpt) THEN
                        temp = temp * 2.0
                        IF (norm2(qpoint)<=1e-8) THEN
                           temp = temp + we1(i) * eigVecCoeffs%abcof(i,lm,1,natom,ilSpin)
                        END IF
                     END IF
                     denCoeffs%mt_lou_coeff(lo,ntyp,1,ilSpinPr,ilSpin) = denCoeffs%mt_lou_coeff(lo,ntyp,1,ilSpinPr,ilSpin) &
                                                                     & + we(i) * CONJG(eigVecCoeffs%ccof(m,i,lo,natom,ilSpinPr)) * temp
                  END DO
               END DO
               DO lop = 1,atoms%nlo(ntyp)
                  IF (atoms%llo(lop,ntyp).EQ.l) THEN
                     DO m = -l,l
                        DO i = 1,ne
                           temp = we(i) * eigVecCoeffs1%ccof(m,i,lo,natom,ilSpin) ! If not DFPT, this is the base case for rhomt(21)
                           IF (l_dfpt) THEN
                              temp = temp * 2.0
                              IF (norm2(qpoint)<=1e-8) THEN
                                 temp = temp + we1(i) * eigVecCoeffs%ccof(m,i,lo,natom,ilSpin)
                              END IF
                           END IF
                           denCoeffs%mt_lolo_coeff(lop,lo,ntyp,ilSpinPr,ilSpin) = denCoeffs%mt_lolo_coeff(lop,lo,ntyp,ilSpinPr,ilSpin) &
                                                                              & + we(i) * CONJG(eigVecCoeffs%ccof(m,i,lop,natom,ilSpinPr)) * temp
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE dfpt_rhomtlo
END MODULE m_dfpt_rhomt
