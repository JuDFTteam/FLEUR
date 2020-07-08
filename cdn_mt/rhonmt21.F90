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
      COMPLEX, PARAMETER :: mi = (0.0,-1.0)
      COMPLEX :: temp(ne)

#include"cpp_double.h"
      COMPLEX CPP_BLAS_cdotc
      EXTERNAL CPP_BLAS_cdotc

      INTEGER jmem,l,lh,llp,llpmax,lm,lmp,lp,lv,m, mp,mv,na,natom,nb,nn,ns,nt!,lplow0,lphi,lplow,lcond

      DO ns=1,sym%nsymt
         !$OMP parallel do default(none) &
         !$OMP private(lh,lp,l,lv,cil,llp,jmem,coef1,mp,lmp,m,lm,coef,mv,temp,na,nt,nn,natom,llpmax) &
         !$OMP shared(sym,we,ne,ns,uunmt21,udnmt21,dunmt21,ddnmt21,atoms,sphhar,eigVecCoeffs) &
         !$OMP collapse(3)
         DO lh = 1,sphhar%nlh(ns)
            DO lp = 0,atoms%lmaxd
               DO l = 0,atoms%lmaxd
                  lv = sphhar%llh(lh,ns)
                  IF ( MOD(lv+l+lp,2) .NE. 0 ) CYCLE
                  cil = mi**(l-lp)
                  DO jmem = 1,sphhar%nmem(lh,ns)
                     mv = sphhar%mlh(jmem,lh,ns)
                     coef1 = cil * sphhar%clnu(jmem,lh,ns) 
                     mp_loop: DO mp = -lp,lp
                        lmp = lp*(lp+1) + mp
                        m_loop: DO m = -l,l
                           coef=  CONJG(coef1 * gaunt1(l,lv,lp,m,mv,mp,atoms%lmaxd))
                           IF (ABS(coef) .LT. 1e-12 ) CYCLE
                           lm= l*(l+1) + m
                           natom= 0
                           DO nn=1,atoms%ntype
                              llp= lp*(atoms%lmax(nn)+1)+l+1
                              llpmax = (atoms%lmax(nn)+1)**2
                              IF(llp.GT.llpmax) CYCLE
                              nt= natom
                              DO na= 1,atoms%neq(nn)
                                 nt= nt+1
                                 IF (sym%ntypsy(nt)==ns) THEN
                                    temp(:) = coef * we(:) * eigVecCoeffs%acof(:,lm,nt,1)
                                    uunmt21(llp,lh,nn) = uunmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%acof(:,lmp,nt,2),1,temp,1)
                                    dunmt21(llp,lh,nn) = dunmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%bcof(:,lmp,nt,2),1,temp,1)
                                    temp(:) = coef * we(:) * eigVecCoeffs%bcof(:,lm,nt,1)
                                    udnmt21(llp,lh,nn) = udnmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%acof(:,lmp,nt,2),1,temp,1)
                                    ddnmt21(llp,lh,nn) = ddnmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%bcof(:,lmp,nt,2),1,temp,1)
                                 ENDIF ! (sym%ntypsy(nt)==ns)
                              ENDDO ! na
                              natom= natom + atoms%neq(nn)
                           ENDDO ! nn
                        ENDDO m_loop ! m
                     ENDDO  mp_loop
                  ENDDO ! jmem
               ENDDO ! lp
            ENDDO ! l
         ENDDO ! lh
         !$OMP end parallel do
      ENDDO ! ns

   END SUBROUTINE rhonmt21
END MODULE m_rhonmt21
