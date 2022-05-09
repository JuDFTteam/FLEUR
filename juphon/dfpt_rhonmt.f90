!-------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!-------------------------------------------------------------------------------

MODULE m_dfpt_rhonmt
   !> Module adapted from rhonmt21, which in turn was adapted from rhonmt, and
   !> from rhonmtlo for the lo part.
   USE m_gaunt,ONLY:gaunt1
   USE m_types_setup
   USE m_types_cdnval
   USE m_constants

   IMPLICIT NONE

CONTAINS

   !! Subroutine to construct all MT density coefficients for a density perturbation
   !! in one routine. The spin input dictates, which is gonna be performed.
   !! The coefficients are of the form:
   !! \begin{align*}
   !! d_{l',l,L,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha} &= \sum_{\nu\bm{k}}\sum_{m\mu(L)}\\
   !! &*  c_{L,\mu}^{*}G_{l,l''(L),l'}^{m,m''(\mu),m-m''(\mu)}A_{l',m-m''(\mu),\lambda'}^{\sigma_{\alpha}',\nu\bm{k}*}\\
   !! &* (2\tilde{f}_{\nu\bm{k}}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\bm{k}\bm{q},j,\beta~(1)}+\tilde{f}_{\nu\bm{k}\bm{q}}^{j,\beta~(1)}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\bm{k}})
   !! \end{align}
   !! The k-point loop is performed outside this routine. In contrast to older
   !! routines, the arrays uunmt etc. and uunmt21 etc. are merged into one
   !! spin and u-order dependent array.
   !!
   !! \(\sigma_{\alpha}(')\): local spin indices \(\rightarrow\) ilSpinPr, ilSpin
   !! \(\lambda(')\): order of the radial basis function (0: u, 1: d)
   !! \(L\): Lattice harmonic index
   !! \(\mu(L)\): Lattice harmonic member index
   !! \(c_{L,\mu}\): Lattice harmonic coefficients
   !! \(\nu\bm{k}\): State index (k-point and number of state)
   !! \(\bm{q},j,\beta\): Perturbation index; q-point, direction and atom
   !! \(\tilde{f}_{\nu\bm{k}}\): Occupation number (smeared)
   !! \(A_{...}^{...}\): Summed matching coefficients and eigenvectors [perturbed for A^{(1)}]
   SUBROUTINE dfpt_rhonmt(atoms,sphhar,we,we1,ne,ilSpinPr,ilSpin,qpoint,l_dfpt,l_less_effort,sym,acof,a1cof,denCoeffs)

      TYPE(t_sym),          INTENT(IN)    :: sym
      TYPE(t_sphhar),       INTENT(IN)    :: sphhar
      TYPE(t_atoms),        INTENT(IN)    :: atoms

      INTEGER,              INTENT(IN)    :: ne, ilSpinPr, ilSpin
      REAL,                 INTENT(IN)    :: we(:),we1(:)!(nobd)
      REAL,                 INTENT(IN)    :: qpoint(3)
      LOGICAL,              INTENT(IN)    :: l_dfpt, l_less_effort

      COMPLEX, INTENT(IN)    :: acof(:,:,0:,:,:)!(nu,lm,iOrd,iAtom,ilSpin)
      COMPLEX, INTENT(IN)    :: a1cof(:,:,0:,:,:)!(nu,lm,iOrd,iAtom,ilSpin)

      TYPE(t_denCoeffs), INTENT(INOUT) :: denCoeffs

      COMPLEX coef, cil, coef1
      COMPLEX :: temp(ne)

!#include"cpp_double.h"
!      COMPLEX CPP_BLAS_cdotc
!      EXTERNAL CPP_BLAS_cdotc

      INTEGER jmem,l,lh,llp,llpmax,lm,lmp,lp,lv,m, mp,mv,na,natom,nn,ns,nt,lphi,lplow

      DO ns=1,sym%nsymt
         !$OMP parallel do default(none) &
         !$OMP private(lh,lp,l,lv,mp,m,mv,lm,lmp,llp,llpmax,lphi,lplow) &
         !$OMP private(cil,jmem,coef1,coef,temp,na,nt,nn,natom) &
         !$OMP shared(sym,we,we1,ne,ns,atoms,sphhar,acof,a1cof,denCoeffs,ilSpinPr,ilSpin,l_dfpt,l_less_effort,qpoint) &
         !$OMP collapse(2)
         DO lh = 1, sphhar%nlh(ns)
            DO l = 0, atoms%lmaxd
               lv = sphhar%llh(lh,ns)
               DO jmem = 1, sphhar%nmem(lh,ns)
                  mv = sphhar%mlh(jmem,lh,ns)
                  m_loop: DO m = -l,l
                     lm= l*(l+1) + m
                     mp = m - mv

                     !maximum value of lp
                     lphi  = l + lv
                     !---> check that lphi is smaller than the max l of the
                     !---> wavefunction expansion
                     lphi = MIN(lphi,atoms%lmaxd)
                     !--->  make sure that l + l'' + lphi is even
                     lphi = lphi - MOD(l+lv+lphi,2)

                     lplow = abs(l-lv)
                     lplow = MAX(lplow,ABS(mp))
                     !---> make sure that l + l'' + lplow is even
                     lplow = lplow + MOD(ABS(lphi-lplow),2)

                     IF (lplow.GT.lphi) CYCLE m_loop

                     DO lp = lplow, lphi,2
                        cil = ImagUnit**(lp-l)
                        coef1 = cil * sphhar%clnu(jmem,lh,ns)
                        lmp = lp*(lp+1) + mp
                        IF (lmp>lm.AND.l_less_effort) CYCLE m_loop

                        coef=  CONJG(coef1 * gaunt1(l,lv,lp,m,mv,mp,atoms%lmaxd))
                        IF (ABS(coef) .LT. 1e-12 ) CYCLE
                        natom= 0
                        DO nn=1,atoms%ntype
                           llp= lp*(atoms%lmax(nn)+1)+l+1
                           llpmax = (atoms%lmax(nn)+1)**2
                           IF(llp.GT.llpmax) CYCLE
                           nt= natom
                           DO na= 1,atoms%neq(nn)
                              nt= nt+1
                              IF (sym%ntypsy(nt)==ns) THEN
                                 ! uu/du
                                 temp(:) = coef * we(:) * a1cof(:,lm,0,nt,ilSpin) ! If not DFPT, this is the base case for rhonmt(21)
                                 IF (lmp/=lm.AND.l_less_effort) temp(:) = temp(:) * 2.0
                                 IF (l_dfpt) THEN
                                    temp(:) = temp(:) * 2.0
                                    IF (norm2(qpoint)<=1e-8) THEN
                                       temp(:) = temp(:) + coef * we1(:) * acof(:,lm,0,nt,ilSpin)
                                    END IF
                                 END IF
                                 denCoeffs%nmt_coeff(llp,lh,nn,0,0,ilSpinPr,ilSpin) = denCoeffs%nmt_coeff(llp,lh,nn,0,0,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,0,nt,ilSpinPr),temp(:ne))
                                 denCoeffs%nmt_coeff(llp,lh,nn,1,0,ilSpinPr,ilSpin) = denCoeffs%nmt_coeff(llp,lh,nn,1,0,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,1,nt,ilSpinPr),temp(:ne))

                                 ! dd/ud
                                 temp(:) = coef * we(:) * a1cof(:,lm,1,nt,ilSpin)
                                 IF (l_dfpt) THEN
                                    temp(:) = temp(:) * 2.0
                                    IF (norm2(qpoint)<=1e-8) THEN
                                       temp(:) = temp(:) + coef * we1(:) * acof(:,lm,1,nt,ilSpin)
                                    END IF
                                 END IF
                                 denCoeffs%nmt_coeff(llp,lh,nn,1,1,ilSpinPr,ilSpin) = denCoeffs%nmt_coeff(llp,lh,nn,1,1,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,1,nt,ilSpinPr),temp(:ne))
                                 denCoeffs%nmt_coeff(llp,lh,nn,0,1,ilSpinPr,ilSpin) = denCoeffs%nmt_coeff(llp,lh,nn,0,1,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,0,nt,ilSpinPr),temp(:ne))
                              ENDIF ! (sym%ntypsy(nt)==ns)
                           ENDDO ! na
                           natom= natom + atoms%neq(nn)
                        END DO ! nn
                     END DO
                  END DO m_loop ! m
               END DO ! jmem
            END DO ! l
         END DO ! lh
         !$OMP end parallel do
      END DO ! ns
   END SUBROUTINE dfpt_rhonmt

   SUBROUTINE dfpt_rhonmtlo(atoms,sphhar,sym,ne,we,we1,acof,a1cof,locof,lo1cof,denCoeffs,ilSpinPr,ilSpin,l_dfpt,qpoint)
      USE m_gaunt,ONLY:gaunt1
      USE m_types
      use m_constants

      IMPLICIT NONE

      TYPE(t_sphhar),       INTENT(IN)    :: sphhar
      TYPE(t_atoms),        INTENT(IN)    :: atoms
      TYPE(t_sym),          INTENT(IN)    :: sym

      INTEGER, INTENT (IN) :: ne, ilSpinPr, ilSpin

      REAL,    INTENT (IN) :: we(:),we1(:)!(nobd)

      REAL,                 INTENT(IN)    :: qpoint(3)
      LOGICAL,              INTENT(IN)    :: l_dfpt

      COMPLEX, INTENT(IN)    :: acof(:,:,0:,:,:)!(nu,lm,iOrd,iAtom,ilSpin)
      COMPLEX, INTENT(IN)    :: a1cof(:,:,0:,:,:)!(nu,lm,iOrd,iAtom,ilSpin)
      COMPLEX, INTENT(IN)    :: locof(:,:,:,:,:)!(m,nu,ilo,iAtom,ilSpin)
      COMPLEX, INTENT(IN)    :: lo1cof(:,:,:,:,:)!(m,nu,ilo,iAtom,ilSpin)

      TYPE(t_denCoeffs), INTENT(INOUT) :: denCoeffs

      COMPLEX :: cmv,fact,cf1, cf2
      INTEGER :: i,jmem,l,lh,lmp,lo,lop,lp,lpmax,lpmax0,lpmin,lpmin0,m,lpp ,mp,mpp,na,neqat0,nn,ntyp

      neqat0 = 0
      DO ntyp = 1,atoms%ntype
         DO lh = 1,sphhar%nlh(sym%ntypsy(neqat0+1))
            lpp = sphhar%llh(lh,sym%ntypsy(neqat0+1))
            DO jmem = 1,sphhar%nmem(lh,sym%ntypsy(neqat0+1))
               mpp = sphhar%mlh(jmem,lh,sym%ntypsy(neqat0+1))
               cmv = CONJG(sphhar%clnu(jmem,lh,sym%ntypsy(neqat0+1)))
               DO lo = 1,atoms%nlo(ntyp)
                  l = atoms%llo(lo,ntyp)
                  lpmin0 = ABS(l-lpp)
                  lpmax0 = l + lpp

                  lpmax = MIN(lpmax0,atoms%lmax(ntyp))
                  lpmax = lpmax - MOD(l+lpp+lpmax,2)
                  DO m = -l,l
                     mp = m - mpp
                     lpmin = MAX(lpmin0,ABS(mp))
                     lpmin = lpmin + MOD(l+lpp+lpmin,2)
                     DO lp = lpmin,lpmax,2
                        lmp = lp* (lp+1) + mp
                        fact = cmv* (ImagUnit** (l-lp))*gaunt1(l,lp,lpp,m,mp,mpp,atoms%lmaxd)
                        na = neqat0
                        DO nn = 1,atoms%neq(ntyp)
                           na = na + 1
                           DO i = 1,ne
                              cf1 = fact * we(i) * lo1cof(m,i,lo,na,ilSpin)! If not DFPT, this is the base case for rhonmtlo
                              IF (l_dfpt) THEN
                                 cf1 = cf1 * 2.0
                                 IF (norm2(qpoint)<=1e-8) THEN
                                    cf1 = cf1 + fact * we1(i) * locof(m,i,lo,na,ilSpin)
                                 END IF
                              END IF
                              denCoeffs%nmt_ulo_coeff(lp,lo,lh,ntyp,0,ilSpinPr,ilSpin) = denCoeffs%nmt_ulo_coeff(lp,lo,lh,ntyp,0,ilSpinPr,ilSpin) &
                                                                                     & + CONJG(acof(i,lmp,0,na,ilSpinPr)) * cf1
                              denCoeffs%nmt_ulo_coeff(lp,lo,lh,ntyp,1,ilSpinPr,ilSpin) = denCoeffs%nmt_ulo_coeff(lp,lo,lh,ntyp,1,ilSpinPr,ilSpin) &
                                                                                     & + CONJG(acof(i,lmp,1,na,ilSpinPr)) * cf1
                           END DO
                        END DO
                     END DO

                     mp = m + mpp
                     lpmin = MAX(lpmin0,ABS(mp))
                     lpmin = lpmin + MOD(l+lpp+lpmin,2)
                     DO lp = lpmin,lpmax,2
                        lmp = lp* (lp+1) + mp
                        fact = cmv* (ImagUnit** (lp-l))*gaunt1(lp,l,lpp,mp,m,mpp,atoms%lmaxd)
                        na = neqat0
                        DO nn = 1,atoms%neq(ntyp)
                           na = na + 1
                           DO i = 1,ne
                              cf1 = fact * we(i) * a1cof(i,lmp,0,na,ilSpin)! If not DFPT, this is the base case for rhonmtlo
                              cf2 = fact * we(i) * a1cof(i,lmp,1,na,ilSpin)
                              IF (l_dfpt) THEN
                                 cf1 = cf1 * 2.0
                                 cf2 = cf2 * 2.0
                                 IF (norm2(qpoint)<=1e-8) THEN
                                    cf1 = cf1 + fact * we1(i) * acof(i,lmp,0,na,ilSpin)
                                    cf2 = cf2 + fact * we1(i) * acof(i,lmp,1,na,ilSpin)
                                 END IF
                              END IF
                              denCoeffs%nmt_lou_coeff(lp,lo,lh,ntyp,0,ilSpinPr,ilSpin) = denCoeffs%nmt_lou_coeff(lp,lo,lh,ntyp,0,ilSpinPr,ilSpin) &
                                                                                     & + CONJG(locof(m,i,lo,na,ilSpinPr)) * cf1
                              denCoeffs%nmt_lou_coeff(lp,lo,lh,ntyp,1,ilSpinPr,ilSpin) = denCoeffs%nmt_lou_coeff(lp,lo,lh,ntyp,1,ilSpinPr,ilSpin) &
                                                                                     & + CONJG(locof(m,i,lo,na,ilSpinPr)) * cf2
                           END DO
                        END DO
                     END DO

                     DO lop = 1,atoms%nlo(ntyp)
                        lp = atoms%llo(lop,ntyp)
                        mp = m - mpp
                        IF ((ABS(l-lpp).LE.lp) .AND.(lp.LE. (l+lpp)) .AND.(MOD(l+lp+lpp,2).EQ.0) .AND.(ABS(mp).LE.lp)) THEN
                           fact = cmv* (ImagUnit** (l-lp))*gaunt1(l,lp,lpp,m,mp,mpp,atoms%lmaxd)
                           na = neqat0
                           DO nn = 1,atoms%neq(ntyp)
                              na = na + 1
                              DO i = 1,ne
                                 cf1 = fact * we(i) * lo1cof(m,i,lo,na,ilSpin)! If not DFPT, this is the base case for rhonmtlo
                                 IF (l_dfpt) THEN
                                    cf1 = cf1 * 2.0
                                    IF (norm2(qpoint)<=1e-8) THEN
                                       cf1 = cf1 + fact * we1(i) * locof(m,i,lo,na,ilSpin)
                                    END IF
                                 END IF
                                 denCoeffs%nmt_lolo_coeff(lop,lo,lh,ntyp,ilSpinPr,ilSpin) = denCoeffs%nmt_lolo_coeff(lop,lo,lh,ntyp,ilSpinPr,ilSpin) &
                                                                                        & + CONJG(locof(mp,i,lop,na,ilSpinPr)) * cf1
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
         neqat0 = neqat0 + atoms%neq(ntyp)
      END DO
   END SUBROUTINE
END MODULE m_dfpt_rhonmt
