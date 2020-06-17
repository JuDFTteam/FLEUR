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
CONTAINS
  SUBROUTINE rhonmt21(atoms,sphhar,we,ne,sym,eigVecCoeffs,uunmt21,udnmt21,dunmt21,ddnmt21)
#include"cpp_double.h"
    USE m_gaunt,ONLY:gaunt1
    USE m_types_setup
    USE m_types_cdnval

    IMPLICIT NONE

    TYPE(t_sym),          INTENT(IN)    :: sym
    TYPE(t_sphhar),       INTENT(IN)    :: sphhar
    TYPE(t_atoms),        INTENT(IN)    :: atoms
    TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs

    !     .. Scalar Arguments ..
    INTEGER,              INTENT(IN)    :: ne   

    !     .. Array Arguments ..
    REAL,                 INTENT(IN)    :: we(:)!(nobd)
    COMPLEX,              INTENT(INOUT) :: uunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
    COMPLEX,              INTENT(INOUT) :: udnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
    COMPLEX,              INTENT(INOUT) :: dunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
    COMPLEX,              INTENT(INOUT) :: ddnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype)
              
    !     .. Local Scalars .. 
    COMPLEX coef, cil, coef1
    COMPLEX, PARAMETER :: mi = (0.0,-1.0)
    COMPLEX CPP_BLAS_cdotc
    EXTERNAL CPP_BLAS_cdotc
 
    INTEGER jmem,l,lh,llp,lm,lmp,lp,lv,m, mp,mv,na,natom,nb,nn,ns,nt!,lplow0,lphi,lplow,lcond
    !     ..
    !
    DO ns=1,sym%nsymt
       natom= 0
       DO nn=1,atoms%ntype
          nt= natom
          DO na= 1,atoms%neq(nn)
             nt= nt+1
             IF (sym%ntypsy(nt)==ns) THEN
                !$OMP PARALLEL DO PRIVATE(lh,lp,l,lv,cil,llp,jmem,coef1,mp,lmp,m,lm,coef,mv) &
                !$OMP DEFAULT(none) &
                !$OMP SHARED(we,ne,na,nt,nn,ns,uunmt21,udnmt21,dunmt21,ddnmt21,atoms,sphhar,eigVecCoeffs) &
                !$OMP collapse(3)
                DO lh = 1,sphhar%nlh(ns)
                   DO lp = 0,atoms%lmax(nn)
                      DO l = 0,atoms%lmax(nn)
                         lv = sphhar%llh(lh,ns)
                         IF ( MOD(lv+l+lp,2) .EQ. 0 ) THEN
                            cil = mi**(l-lp)
                            llp= lp*(atoms%lmax(nn)+1)+l+1
                            DO jmem = 1,sphhar%nmem(lh,ns)
                               mv = sphhar%mlh(jmem,lh,ns)
                               coef1 = cil * sphhar%clnu(jmem,lh,ns) 
                               DO mp = -lp,lp
                                  lmp = lp*(lp+1) + mp
                                  m_loop: DO m = -l,l
                                  IF(.NOT.(l.GE.m.AND.lp.GE.mp.AND.lv.GE.mv)) CYCLE m_loop
                                     coef=  CONJG(coef1 * gaunt1(l,lv,lp,m,mv,mp,atoms%lmaxd))
                                     IF (ABS(coef) .GE. 0 ) THEN
                                        lm= l*(l+1) + m
                                        uunmt21(llp,lh,nn) = uunmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%acof(:,lmp,nt,2),1, we(:) * coef * eigVecCoeffs%acof(:,lm,nt,1),1)
                                        udnmt21(llp,lh,nn) = udnmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%acof(:,lmp,nt,2),1,we(:) * coef * eigVecCoeffs%bcof(:,lm,nt,1),1)
                                        dunmt21(llp,lh,nn) = dunmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%bcof(:,lmp,nt,2),1, we(:) * coef * eigVecCoeffs%acof(:,lm,nt,1),1)
                                        ddnmt21(llp,lh,nn) = ddnmt21(llp,lh,nn) + CPP_BLAS_cdotc(ne,eigVecCoeffs%bcof(:,lmp,nt,2),1,we(:) * coef * eigVecCoeffs%bcof(:,lm,nt,1),1)
                                     ENDIF ! (coef >= 0)
                                  ENDDO m_loop ! m
                               ENDDO  ! mp
                            ENDDO ! jmem
                         ENDIF ! ( MOD(lv+l+lp),2) == 0 )
                      ENDDO ! lp
                   ENDDO ! l
                ENDDO ! lh
                !$OMP END PARALLEL DO
             ENDIF ! (sym%ntypsy(nt)==ns)
          ENDDO ! na
          natom= natom + atoms%neq(nn)
       ENDDO ! nn

    ENDDO ! ns

    RETURN

  END SUBROUTINE rhonmt21
END MODULE m_rhonmt21
