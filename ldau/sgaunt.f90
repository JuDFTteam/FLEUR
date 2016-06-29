!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_sgaunt
CONTAINS
  SUBROUTINE sgaunt(&
       lmaxw,lmmax,lmm,&
       c)
    !************************************************************
    !*   Calculation of the Gaunt coefficients C(L2M2,L1M1,LM)  *
    !*                                                          *
    !*    l2m2                    /         *                   *
    !*   C       =C(l2m2,l1m1,lm)=\dr*Y(r)*Y(r)*Y(r)            *
    !*    lm,l1m1                 /    lm   l1m1 l2m2           *
    !*                                                          *
    !*    and C.ne.0 when l2=/l1-l/,/l1-l/+2,...,l1+l,m2=m1-m   *
    !*    Y(lm) is a complex spherical garmonic with a phase    *
    !*    after Condon and Shortley                             *
    !* Written by S.Yu.Savrasov (P.N.Lebedev Physical Institute)*
    !************************************************************
    !*    called by umtx() ; part of the LDA+U package          *
    !*                                          G.B., Oct. 2000 *
    !************************************************************

    USE m_clebsch
    USE m_constants,ONLY:pi_const
    IMPLICIT NONE
    INTEGER, INTENT (IN)  :: lmaxw,lmmax,lmm

    REAL,    INTENT (OUT) :: c(0:2*lmaxw+1,lmmax,lmmax)
    INTEGER l1,m1,l ,l1m1,l2,lm,m2,ll2,m
    REAL    aj,bj,am,bm,cj,cm,dl1,dl2,dl3,a1,a2

    DO l1 = 0,lmm
       DO m1 = -l1,l1
          DO l = 0,lmm
             DO m = -l,l
                l1m1 = l1*(l1+1)+m1+1
                lm   = l*(l+1)+m+1
                DO l2 = ABS(l1-l),l1+l,2
                   ll2 = l2/2
                   m2  = m1-m
                   IF (ABS(m2).LE.l2) THEN     !!! selection rule
                      aj = REAL(l)
                      bj = REAL(l2)
                      am = REAL(m)
                      bm = REAL(m2)
                      cj = REAL(l1)
                      cm = REAL(m1)
                      a1 = clebsch(aj,bj,0.0,0.0,cj,0.0)  !!! Clebsch-Gordan coefficients
                      a2 = clebsch(aj,bj,am,bm,cj,cm)     !!! Clebsch-Gordan coefficients
                      dl1 = REAL(2*l +1)
                      dl2 = REAL(2*l2+1)
                      dl3 = REAL(2*l1+1)
                      c(ll2,l1m1,lm)=a1*a2*SQRT(dl1*dl2/dl3/4.0/pi_const)
                   ELSEIF (ABS(m2).GT.l2)THEN
                      c(ll2,l1m1,lm)=0.0
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE sgaunt
END MODULE m_sgaunt
