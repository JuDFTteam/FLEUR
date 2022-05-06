!-------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!-------------------------------------------------------------------------------

MODULE m_dfpt_rhonmt
   !> Module adapted from rhonmt21, which in turn was adapted from rhonmt.
   USE m_gaunt,ONLY:gaunt1
   USE m_types_setup
   USE m_types_cdnval
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE dfpt_rhonmt(atoms,sphhar,we,we1,ne,ilSpinPr,ilSpin,sym,acof,a1cof,nmt_coeff)
      !> Subroutine to constuct all MT density coefficients for a density perturbation
      !! in one routine. The spin input dictates, which is gonna be performed.
      !! The coefficients are of the form
      !! \f{eqnarray*}{
      !! d_{l',l,L,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha} &=& \sum_{\nu\bm{k}}\sum_{m\mu(L)}\\
      !! &*&  c_{L,\mu}^{*}G_{l,l''(L),l'}^{m,m''(\mu),m-m''(\mu)}A_{l',m-m''(\mu),\lambda'}^{\sigma_{\alpha}',\nu\bm{k}}\\
      !! &*& (2\tilde{f}_{\nu\bm{k}}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\bm{k}\bm{q},j,\beta~(1)}+\tilde{f}_{\nu\bm{k}\bm{q}}^{j,\beta~(1)}A_{l,m,\lambda}^{\sigma_{\alpha},\nu\bm{k}}) \\
      !! \f}
      !! The k-point loop is performed outside this routine. In contrast to older
      !! routines, the arrays uunmt etc. and uunmt21 etc. are merged into one
      !! spin and u-order dependent array.
      !!
      !! \f$\sigma_{\alpha}(')\f$: local spin indices \f$\rightarrow\f$ ilSpinPr, ilSpin
      !! \f$\lambda(')\f$: order of the radial basis function (0: u, 1: d)
      !! \f$L\f$: Lattice harmonic index
      !! \f$\mu(L)\f$: Lattice harmonic member index
      !! \f$c_{L,\mu}\f$: Lattice harmonic coefficients
      !! \f$\nu\bm{k}\f$: State index (k-point and number of state)
      !! \f$\bm{q},j,\beta\f$: Perturbation index; q-point, direction and atom
      !! \f$\tilde{f}_{\nu\bm{k}}\f$: Occupation number (smeared)
      !! \f$A_{...}^{...}\f$: Summed matching coefficients and eigenvectors [perturbed for A^{(1)}]

      TYPE(t_sym),          INTENT(IN)    :: sym
      TYPE(t_sphhar),       INTENT(IN)    :: sphhar
      TYPE(t_atoms),        INTENT(IN)    :: atoms

      INTEGER,              INTENT(IN)    :: ne, ilSpinPr, ilSpin
      REAL,                 INTENT(IN)    :: we(:),we1(:)!(nobd)

      COMPLEX, INTENT(IN)    :: acof(:,:,:,:,:)!(nu,lmp,iOrd,iAtom,ilSpin)
      COMPLEX, INTENT(IN)    :: a1cof(:,:,:,:,:)!(nu,lmp,iOrd,iAtom,ilSpin)
      COMPLEX, INTENT(INOUT) :: nmt_coeff(:,:,:,:,:,:,:)!llp,lh,iAtom,iOrdPr,iOrd,ilSpinPr,ilSpin

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
         !$OMP shared(sym,we,we1,ne,ns,atoms,sphhar,acof,a1cof,nmt_coeff,ilSpinPr,ilSpin) &
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
                                 temp(:) = coef * 2 * we(:) * a1cof(:,lm,0,nt,ilSpin)
                                 nmt_coeff(llp,lh,nn,0,0,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,0,0,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,0,nt,ilSpinPr),temp(:ne))
                                 nmt_coeff(llp,lh,nn,1,0,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,1,0,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,1,nt,ilSpinPr),temp(:ne))

                                 temp(:) = coef * we1(:) * acof(:,lm,0,nt,ilSpin)
                                 nmt_coeff(llp,lh,nn,0,0,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,0,0,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,0,nt,ilSpinPr),temp(:ne))
                                 nmt_coeff(llp,lh,nn,1,0,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,1,0,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,1,nt,ilSpinPr),temp(:ne))
                                 ! dd/ud
                                 temp(:) = coef * 2 * we(:) * a1cof(:,lm,1,nt,ilSpin)
                                 nmt_coeff(llp,lh,nn,1,1,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,1,1,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,1,nt,ilSpinPr),temp(:ne))
                                 nmt_coeff(llp,lh,nn,0,1,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,0,1,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,0,nt,ilSpinPr),temp(:ne))

                                 temp(:) = coef * we1(:) * acof(:,lm,1,nt,ilSpin)
                                 nmt_coeff(llp,lh,nn,1,1,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,1,1,ilSpinPr,ilSpin) &
                                                                        & + dot_product(acof(:ne,lmp,1,nt,ilSpinPr),temp(:ne))
                                 nmt_coeff(llp,lh,nn,0,1,ilSpinPr,ilSpin) = nmt_coeff(llp,lh,nn,0,1,ilSpinPr,ilSpin) &
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
END MODULE m_dfpt_rhonmt
