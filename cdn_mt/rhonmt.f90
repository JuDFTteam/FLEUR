MODULE m_rhonmt
CONTAINS
  SUBROUTINE rhonmt(atoms,sphhar, we,ne,sym, acof,bcof,&
       uunmt,ddnmt,udnmt,dunmt)
    !     *************************************************************
    !     subroutine sets up the coefficients of non-sphereical
    !     muffin-tin density                          c.l.fu
    !     *************************************************************
    USE m_gaunt,ONLY:gaunt1
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN) :: ne  
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT(IN) :: acof(:,0:,:)!(nobd,0:lmaxd* (lmaxd+2),natd)
    COMPLEX, INTENT(IN) :: bcof(:,0:,:)
    REAL,    INTENT(IN) :: we(:)!(nobd)
    REAL,INTENT(INOUT) :: ddnmt(0:,:,:)!(0:(lmaxd* (lmaxd+3))/2,nlhd,ntypd)
    REAL,INTENT(INOUT) :: dunmt(0:,:,:)
    REAL,INTENT(INOUT) :: udnmt(0:,:,:)
    REAL,INTENT(INOUT) :: uunmt(0:,:,:)
    !     ..
    !     .. Local Scalars ..
    COMPLEX cconst,cil,cmv,ci
    REAL coef
    INTEGER  :: jmem,l,lcond,lh,llp,llpmax,lm,lmp,lp,lphi,lplow,lplow0,lv
    INTEGER  :: mp,mv,na,natom,nb,nn,ns,nt,m
    !     ..
    !     ..

    ci = cmplx(0.0,1.0)
    !
    !Initialize private variables in gaunt module before parallel region
    !$      coef = gaunt1(0,0,0,0,0,0,atoms%lmaxd)

    DO ns = 1,sym%nsymt
       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP&  PRIVATE(lv,jmem,mv,cmv,l,lm,mp,m,llpmax,nt,na,nb,lplow0)&
       !$OMP&  PRIVATE(lphi,lplow,lcond,lp,cil,lmp,llp,coef,cconst&
       !$OMP& ,natom,nn) 
       DO  lh = 1,sphhar%nlh(ns)
          lv = sphhar%llh(lh,ns)
          DO  jmem = 1,sphhar%nmem(lh,ns)
             mv = sphhar%mlh(jmem,lh,ns)
             cmv = conjg(sphhar%clnu(jmem,lh,ns))
             DO  l = 0,atoms%lmaxd
                m_loop: DO m = -l,l
                   lm = l* (l+1) + m
                   mp = m - mv
                   !     -----> set up the lower and upper limit of lp in such a way that
                   !     -----> lp+l+lv is even, lp<=l, and (lp,l,lv) satisfies the
                   !     -----> triangular relation
                   lplow0 = iabs(l-lv)
                   lphi = l - mod(lv,2)
                   lplow = max(lplow0,iabs(mp))
                   lcond = iabs(lphi-lplow)
                   lplow = lplow + mod(lcond,2)
                   IF (lplow.GT.lphi) CYCLE m_loop
                   DO  lp = lplow,lphi,2
                      cil = ci** (l-lp)
                      lmp = lp* (lp+1) + mp
                      IF (lmp.GT.lm) CYCLE m_loop
                      llp = (l* (l+1))/2 + lp
                      !     -----> gaunt's coefficient
                      coef = 2.*gaunt1(l,lv,lp,m,mv,mp,atoms%lmaxd)
                      IF (lmp.EQ.lm) coef = coef/2.
                      cconst = coef* (cil*cmv)
                      natom = 0
                      DO  nn = 1,atoms%ntype
                         llpmax = (atoms%lmax(nn)* (atoms%lmax(nn)+3))/2
                         IF (llp.LE.llpmax) THEN 
                            nt = natom
                            DO  na = 1,atoms%neq(nn)
                               nt = nt + 1
                               IF (atoms%ntypsy(nt).EQ.ns) THEN
                                  DO  nb = 1,ne
                                     uunmt(llp,lh,nn) = uunmt(llp,lh,nn)&
                                          +we(nb)*real(cconst*acof(nb,lm,nt)*conjg(acof(nb,lmp,nt)))
                                     ddnmt(llp,lh,nn) = ddnmt(llp,lh,nn) +&
                                          we(nb)*real(cconst*bcof(nb,lm,nt)*conjg(bcof(nb,lmp,nt)))
                                     udnmt(llp,lh,nn) = udnmt(llp,lh,nn) +&
                                          we(nb)*real(cconst*acof(nb,lm,nt)*conjg(bcof(nb,lmp,nt)))
                                     dunmt(llp,lh,nn) = dunmt(llp,lh,nn) +&
                                          we(nb)*real(cconst*bcof(nb,lm,nt)*conjg(acof(nb,lmp,nt)))
                                  ENDDO
                               ENDIF
                            ENDDO
                         ENDIF
                         natom = natom + atoms%neq(nn)
                      ENDDO
                   ENDDO
                ENDDO m_loop
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDDO

  END SUBROUTINE rhonmt
END MODULE m_rhonmt
