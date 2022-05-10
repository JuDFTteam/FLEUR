!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rhonmt21
   !     *************************************************************
   !     subroutine sets up the coefficients of the spin (up,down)
   !     part of the non-spherical muffin-tin density.
   !                                                 pk`00 ff`01 gb`02
   !     Added parallelization and reworked for the efficient use with FFN.
   !                                                 R. Hilgers July '20
   !     *************************************************************
   USE m_gaunt,ONLY:gaunt1
   USE m_types_setup
   USE m_types_cdnval
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE rhonmt21(atoms,sphhar,we,ne,sym,eigVecCoeffs,uunmt21,udnmt21,dunmt21,ddnmt21)


      TYPE(t_sym),          INTENT(IN)    :: sym
      TYPE(t_sphhar),       INTENT(IN)    :: sphhar
      TYPE(t_atoms),        INTENT(IN)    :: atoms
      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs

      !     .. Scalar Arguments ..
      INTEGER,              INTENT(IN)    :: ne

      !     .. Array Arguments ..
      REAL,                 INTENT(IN)    :: we(:)!(nobd)
      COMPLEX,              INTENT(INOUT) :: uunmt21(:,:,:)!((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
      COMPLEX,              INTENT(INOUT) :: udnmt21(:,:,:)!((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
      COMPLEX,              INTENT(INOUT) :: dunmt21(:,:,:)!((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
      COMPLEX,              INTENT(INOUT) :: ddnmt21(:,:,:)!((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)

      !     .. Local Scalars ..
      COMPLEX coef, cil, coef1
      COMPLEX :: temp(ne)

#include"cpp_double.h"
      COMPLEX CPP_BLAS_cdotc
      EXTERNAL CPP_BLAS_cdotc

      INTEGER jmem,l,lh,llp,llpmax,lm,lmp,lp,lv,m, mp,mv,na,natom,nb,nn,ns,nt,lphi,lplow

      DO ns=1,sym%nsymt
         !$OMP parallel do default(none) &
         !$OMP private(lh,lp,l,lv,mp,m,mv,lm,lmp,llp,llpmax,lphi,lplow) &
         !$OMP private(cil,jmem,coef1,coef,temp,na,nt,nn,natom) &
         !$OMP shared(sym,we,ne,ns,uunmt21,udnmt21,dunmt21,ddnmt21,atoms,sphhar,eigVecCoeffs) &
         !$OMP collapse(2)
         DO lh = 1,sphhar%nlh(ns)
            DO l = 0,atoms%lmaxd
               lv = sphhar%llh(lh,ns)
               DO jmem = 1,sphhar%nmem(lh,ns)
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
                                 temp(:) = coef * we(:) * eigVecCoeffs%abcof(:,lm,0,nt,1)
                                 !uunmt21(llp,lh,nn) = uunmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%acof(:,lmp,nt,2),1,temp,1)
                                 !dunmt21(llp,lh,nn) = dunmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%bcof(:,lmp,nt,2),1,temp,1)

                                 uunmt21(llp,lh,nn) = uunmt21(llp,lh,nn) + dot_product(eigVecCoeffs%abcof(:ne,lmp,0,nt,2),temp(:ne))
                                 dunmt21(llp,lh,nn) = dunmt21(llp,lh,nn) + dot_product(eigVecCoeffs%abcof(:ne,lmp,1,nt,2),temp(:ne))

                                 temp(:) = coef * we(:) * eigVecCoeffs%abcof(:,lm,1,nt,1)
                                 !udnmt21(llp,lh,nn) = udnmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%acof(:,lmp,nt,2),1,temp,1)
                                 !ddnmt21(llp,lh,nn) = ddnmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%bcof(:,lmp,nt,2),1,temp,1)

                                 udnmt21(llp,lh,nn) = udnmt21(llp,lh,nn) + dot_product(eigVecCoeffs%abcof(:ne,lmp,0,nt,2),temp(:ne))
                                 ddnmt21(llp,lh,nn) = ddnmt21(llp,lh,nn) + dot_product(eigVecCoeffs%abcof(:ne,lmp,1,nt,2),temp(:ne))
                              ENDIF ! (sym%ntypsy(nt)==ns)
                           ENDDO ! na
                           natom= natom + atoms%neq(nn)
                        ENDDO ! nn
                     ENDDO
                  ENDDO m_loop ! m
               ENDDO ! jmem
            ENDDO ! l
         ENDDO ! lh
         !$OMP end parallel do
      ENDDO ! ns

   END SUBROUTINE rhonmt21
END MODULE m_rhonmt21
