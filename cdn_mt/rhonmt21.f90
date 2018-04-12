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
  !     *************************************************************
CONTAINS
  SUBROUTINE rhonmt21(atoms,sphhar,we,ne,sym,eigVecCoeffs,denCoeffsOffdiag)

    USE m_gaunt,ONLY:gaunt1
    USE m_types

    IMPLICIT NONE

    TYPE(t_sym),INTENT(IN)                 :: sym
    TYPE(t_sphhar),INTENT(IN)              :: sphhar
    TYPE(t_atoms),INTENT(IN)               :: atoms
    TYPE(t_eigVecCoeffs),INTENT(IN)        :: eigVecCoeffs
    TYPE(t_denCoeffsOffdiag),INTENT(INOUT) :: denCoeffsOffdiag

    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN) :: ne   

    !     .. Array Arguments ..
    REAL,    INTENT(IN) :: we(:)!(nobd)

    !     .. Local Scalars ..
    COMPLEX coef, cconst, cil, coef1
    COMPLEX, PARAMETER :: mi = (0.0,-1.0)
    INTEGER jmem,l,lh,llp,lm,lmp,lp,lv,m, mp,mv,na,natom,nb,nn,ns,nt
    !     ..
    !
    DO ns=1,sym%nsymt
       natom= 0
       DO nn=1,atoms%ntype
          nt= natom
          DO na= 1,atoms%neq(nn)
             nt= nt+1
             IF (atoms%ntypsy(nt)==ns) THEN

                DO lh = 1,sphhar%nlh(ns)
                   lv = sphhar%llh(lh,ns)
                   DO lp = 0,atoms%lmax(nn)
                      DO l = 0,atoms%lmax(nn)

                         IF ( MOD(lv+l+lp,2) == 0 ) THEN
                            cil = mi**(l-lp)
                            llp= lp*(atoms%lmax(nn)+1)+l+1

                            DO jmem = 1,sphhar%nmem(lh,ns)
                               mv = sphhar%mlh(jmem,lh,ns)
                               coef1 = cil * sphhar%clnu(jmem,lh,ns) 
                               DO mp = -lp,lp
                                  lmp = lp*(lp+1) + mp
                                  DO m = -l,l
                                     lm= l*(l+1) + m
                                     coef=  CONJG( coef1 *gaunt1(l,lv,lp,m,mv,mp,atoms%lmaxd) )

                                     IF (ABS(coef) >= 0 ) THEN
                                        DO nb = 1,ne
                                           cconst= we(nb) * coef
                                           denCoeffsOffdiag%uunmt21(llp,lh,nn) = denCoeffsOffdiag%uunmt21(llp,lh,nn)+ &
                                                cconst * eigVecCoeffs%acof(nb,lm,nt,1)*CONJG(eigVecCoeffs%acof(nb,lmp,nt,2))
                                           denCoeffsOffdiag%udnmt21(llp,lh,nn) = denCoeffsOffdiag%udnmt21(llp,lh,nn)+&
                                                cconst * eigVecCoeffs%bcof(nb,lm,nt,1)*CONJG(eigVecCoeffs%acof(nb,lmp,nt,2))
                                           denCoeffsOffdiag%dunmt21(llp,lh,nn) = denCoeffsOffdiag%dunmt21(llp,lh,nn)+&
                                                cconst * eigVecCoeffs%acof(nb,lm,nt,1)*CONJG(eigVecCoeffs%bcof(nb,lmp,nt,2))
                                           denCoeffsOffdiag%ddnmt21(llp,lh,nn) = denCoeffsOffdiag%ddnmt21(llp,lh,nn)+&
                                                cconst * eigVecCoeffs%bcof(nb,lm,nt,1)*CONJG(eigVecCoeffs%bcof(nb,lmp,nt,2))
                                        ENDDO ! nb
                                     ENDIF ! (coef >= 0)

                                  ENDDO ! mp
                               ENDDO ! m
                            ENDDO ! jmem

                         ENDIF ! ( MOD(lv+l+lp),2) == 0 )

                      ENDDO ! lp
                   ENDDO ! l
                ENDDO ! lh

             ENDIF ! (atoms%ntypsy(nt)==ns)
          ENDDO ! na
          natom= natom + atoms%neq(nn)
       ENDDO ! nn

    ENDDO ! ns

    RETURN

  END SUBROUTINE rhonmt21
END MODULE m_rhonmt21
