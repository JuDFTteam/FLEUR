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
  SUBROUTINE rhonmt21(atoms,llpd,sphhar, we,ne,sym,&
       acof,bcof, uunmt21,ddnmt21,udnmt21,dunmt21)
    USE m_gaunt,ONLY:gaunt1
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN) :: llpd    
    INTEGER,INTENT(IN) :: ne   
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT(IN) :: acof(:,0:,:,:)!(nobd,0:lmaxd* (lmaxd+2),natd,jspd)
    COMPLEX, INTENT(IN) :: bcof(:,0:,:,:)
    REAL,    INTENT(IN) :: we(:)!(nobd)
    COMPLEX, INTENT (INOUT) :: ddnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype  )
    COMPLEX, INTENT (INOUT) :: dunmt21((atoms%lmaxd+1)**2 ,sphhar%nlhd,atoms%ntype )
    COMPLEX, INTENT (INOUT) :: udnmt21((atoms%lmaxd+1)**2 ,sphhar%nlhd,atoms%ntype )
    COMPLEX, INTENT (INOUT) :: uunmt21((atoms%lmaxd+1)**2 ,sphhar%nlhd,atoms%ntype )
    !     ..
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
                                           uunmt21(llp,lh,nn) = uunmt21(llp,lh,nn)+ &
                                                cconst * acof(nb,lm,nt,1)*CONJG(acof(nb,lmp,nt,2))
                                           udnmt21(llp,lh,nn) = udnmt21(llp,lh,nn)+&
                                                cconst * bcof(nb,lm,nt,1)*CONJG(acof(nb,lmp,nt,2))
                                           dunmt21(llp,lh,nn) = dunmt21(llp,lh,nn)+&
                                                cconst * acof(nb,lm,nt,1)*CONJG(bcof(nb,lmp,nt,2))
                                           ddnmt21(llp,lh,nn) = ddnmt21(llp,lh,nn)+&
                                                cconst * bcof(nb,lm,nt,1)*CONJG(bcof(nb,lmp,nt,2))
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
