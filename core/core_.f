      MODULE m_core
c...............................................................core
      CONTAINS
      SUBROUTINE core(
     >                mrad,vt,bt,zz,stval,dx,nlshell,nqntab,lqntab,jtop,
     X                ectab,
     <                rhochr,rhospn)
      USE m_felim
      USE m_findlim
      USE m_cnodes
      USE m_coredir
      USE m_ccsdnt
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad,jtop,nlshell
      REAL,    INTENT (IN) :: dx,stval,zz
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lqntab(15),nqntab(15)
      REAL,    INTENT (IN) :: bt(mrad),vt(mrad)
      REAL,    INTENT (OUT) :: rhochr(mrad),rhospn(mrad)
      REAL,    INTENT (INOUT) :: ectab(100)
C     ..
C     .. Local Scalars ..
      REAL bhf,bsum,btot,dvstep,ec,elim,qcor,spcor,tolvar,vz,
     +     xmj,ec_sv
      INTEGER i,ic,ic1,ic2,ie,iflag,ii,ilshell,imin,ir,ish,istart,iter,
     +        itermax,iv,j,jv,kap1,kap2,l,lll,muem05,n,ncor,nmatch,node,
     +        nqn,nrmax,nsh,nsol,nval,nvar,nzero,s,t
      LOGICAL ferro
C     ..
C     .. Local Arrays ..
      REAL bb(mrad),bhff(15),dedv(4,4),dv(4),dvde(4,4),err(4),errnew(4),
     +     fc(2,2,mrad),fck(2,2,mrad),gc(2,2,mrad),gck(2,2,mrad),
     +     piw(2,2),pow(2,2),qiw(2,2),qow(2,2),rc(mrad),rc2(mrad),
     +     rchr(mrad),rspn(mrad),var(4),varnew(4),vv(mrad)
      INTEGER kap(2)
      CHARACTER txtl(0:3)
      CHARACTER*3 txtk(4)
C     ..
C     .. External Functions ..
      REAL rsimp
      EXTERNAL rsimp
C     ..
C     .. External Subroutines ..
      EXTERNAL cfnorm,coreerr,corehff,nwrfst,rinvgj
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,iabs,int,sign
C     ..
C     .. Data statements ..

      DATA txtl/'s','p','d','f'/
      DATA txtk/'1/2','3/2','5/2','7/2'/
!      DATA nqntab/1,2,2,3,3,3,4,4,4,4,5,5,5,6,6/
!      DATA lqntab/0,0,1,0,1,2,0,1,2,3,0,1,2,0,1/
      DATA tolvar/1.0E-12/
      DATA dvstep/0.010/
      DATA itermax/150/
C     ..
c
      nrmax = mrad
      DO n = 1,nrmax
         rc(n) = exp(stval+ (n-1)*dx)
         rc2(n) = rc(n)**2
      END DO
c
      DO n = 1,nrmax
         rhochr(n) = 0.00
         rhospn(n) = 0.00
      END DO
c
      bsum = 0.0
      DO n = 1,jtop
         vv(n) = vt(n)
         bb(n) = bt(n)
         bsum = bsum + abs(bb(n))
      END DO
C

      ncor = 0
      DO 10 ilshell = 1,nlshell
         ncor = ncor + 2* (2*lqntab(ilshell)+1)
   10 CONTINUE
      nval = int(zz) - ncor

      IF (bsum.GT.1.0E-8) THEN
         ferro = .true.
      ELSE
         ferro = .false.
         WRITE (6,FMT=
     +'(''  PARAMAGNETIC CASE'',/,
     +   ''  *****************'')')
      END IF

      WRITE (6,FMT=
     +'(/,''  ATOMIC NUMBER  : '',F6.2,/,
     +    ''  CORE ELECTRONS : '',I5,/,
     +    ''  VAL. ELECTRONS : '',I5)') zz,ncor,nval
      WRITE (6,FMT=
     +'(  /,
     +    ''  MESH    RC(1)  : '',F12.6,//10x,
     +    ''MUE  KAP ITER    ENERGY '',
     +    ''(Hart)'', " HF (KG) ")') rc(1)


C                   ---------------------------------------
C--->               INITIALIZE QUANTUM NUMBERS  NQN  AND  L
C                   ---------------------------------------
      btot = 0.0
      ic = 0
C ---> nl  - loop
      DO 190 ilshell = 1,nlshell
         nqn = nqntab(ilshell)
         l = lqntab(ilshell)
         nsh = 2* (2*l+1)
         bhff(ilshell) = 0.0
         ish = 0
C ---> \mu  - loop
         DO 180 muem05 = -l - 1,+l
            xmj = muem05 + 0.5
            kap1 = -l - 1
            kap2 = l
            kap(1) = kap1
            kap(2) = kap2
c
            lll = l* (l+1)
            IF (abs(xmj).GT.l) THEN
               nsol = 1
            ELSE
               nsol = 2
            END IF
c
            nvar = 2*nsol
c
c----> s - loop : solutions for each (nl,\mu)
c
            DO 170 s = 1,nsol
c
               DO i = 1,2
                  DO j = 1,2
                     DO n = 1,mrad
                        gc(i,j,n) = 0.0
                        gck(i,j,n) = 0.0
                        fc(i,j,n) = 0.0
                        fck(i,j,n) = 0.0
                     END DO
                  END DO
               END DO
c
               ic = ic + 1
               ish = ish + 1
               t = 3 - s
               ec = ectab(ic)
               write(*,*) ic,ec
c+gu
               IF (ic > 1) THEN
                 IF (ectab(ic-1) > 0.1) THEN
                   vv(:) = vv(:) + 2*ec_sv
                   bb(:) = bb(:) + 2*ec_sv
                 ENDIF
               ENDIF
               IF (ec.GT.0.1) THEN !  GOTO 170
                  vv(:) = vv(:) - 2*ectab(ic)
                  bb(:) = bb(:) - 2*ectab(ic)
                  ec_sv = ectab(ic)
                  ec = -ectab(ic)
               ELSE
                 ec_sv = 0.0
               ENDIF
               istart = 1
c
               IF (ish.GT.1) GO TO 30
c
c energy limits : elim and emax
c
               CALL felim(
     >                    mrad,lll,zz,nqn,vv,rc,
     <                    elim)
               write(*,*) 'elim',elim
   20          IF (ec.LE.elim) ec = elim*0.70
c
c turning point and physical infinity : nmatch and nzero
c
               CALL findlim(
     >                      mrad,lll,ec,vv,rc,
     <                      nmatch,nzero)
               write(*,*) 'nmatch,nzero',nmatch,nzero
   30          CONTINUE
c
c nodes for large component : if IFLAG=0 number of nodes is OK
c
               iflag = 0
               CALL cnodes(mrad,iflag,s,ec,l,xmj,nqn,vv,bb,rc,dx,
     +                     nmatch,nzero,gc,fc,pow,qow,piw,qiw,node)
               IF (iflag.NE.0) GO TO 20
c
c setup of Newton-Raphson algorithm
c
               CALL nwrfst(mrad,nsol,s,t,nmatch,nzero,ferro,ec,rc,
     +                     pow,piw,gc,err,var,dv,varnew,errnew)
c
               CALL coreerr(err,var,s,nsol,pow,qow,piw,qiw)
               DO 40 iv = 1,nvar
                  dv(iv) = var(iv)
   40          CONTINUE
c
               iter = 0
C---> iterations start
   50          iter = iter + 1
C                         ----------------------------------
C                         CHECK WHETHER NUMBER OF NODES O.K.
C                         ----------------------------------
               IF (iter.GT.1) THEN
                  node = 0
                  DO 60 n = 2, (nmatch-1)
                     IF (gc(s,s,n)*gc(s,s,n-1).LT.0.0) node = node + 1
   60             CONTINUE
c
                  IF (node.NE. (nqn-l-1)) THEN
                     IF (node.GT. (nqn-l-1)) THEN
                        ec = 1.2*ec
                        write (*,*) 'up',node,ec
                     ELSE
                        ec = 0.8*ec
                        write (*,*) 'do',node,ec
                     END IF
                     istart = istart + 1
                     IF (istart.LT.20) GO TO 20
                  END IF
               END IF
c                    - atomic energy search -
               DO 90 iv = 2,nvar
                  DO 70 jv = 1,nvar
                     varnew(jv) = var(jv)
   70             CONTINUE
                  varnew(iv) = var(iv) + dv(iv)*dvstep
c
                  IF (abs(var(iv)).EQ.0.00) THEN
                     IF (ferro) THEN
                        WRITE (6,FMT='(A,I3,A)') ' VAR ',iv,
     +                    ' = 0 ??????!!!!!'
                     END IF
                  ELSE
                     IF (abs(dv(iv)/var(iv)).LT.
     +                   tolvar) varnew(iv) = var(iv)*
     +                   (1.00+sign(dvstep*tolvar,dv(iv)))
                  END IF
c
                  CALL coreerr(errnew,varnew,s,nsol,pow,qow,piw,qiw)
c
                  DO 80 ie = 1,nvar
                     IF ((errnew(ie)-err(ie)).EQ.0.00) THEN
                        dedv(ie,iv) = 0.00
                        IF ((ie.EQ.iv) .AND. .NOT.ferro) THEN
                           dedv(ie,iv) = 1.00
                        END IF
                     ELSE
                        dedv(ie,iv) = (errnew(ie)-err(ie))/
     +                                (varnew(iv)-var(iv))
                     END IF
   80             CONTINUE
   90          CONTINUE
c
               DO 100 jv = 1,nvar
                  varnew(jv) = var(jv)
  100          CONTINUE
               varnew(1) = var(1) + dv(1)*dvstep
               IF (abs(dv(1)/var(1)).LT.tolvar) varnew(1) = var(1)*
     +             (1.00+sign(dvstep*tolvar,dv(1)))
c
               CALL coredir(mrad,varnew(1),l,xmj,1,vv,bb,rc,dx,
     +                      nmatch,nzero,gc,fc,pow,qow,piw,qiw)
               CALL coredir(mrad,varnew(1),l,xmj,2,vv,bb,rc,dx,
     +                      nmatch,nzero,gc,fc,pow,qow,piw,qiw)
c
               CALL coreerr(errnew,varnew,s,nsol,pow,qow,piw,qiw)
c
               DO 110 ie = 1,nvar
                  dedv(ie,1) = (errnew(ie)-err(ie))/ (varnew(1)-var(1))
  110          CONTINUE
c
               CALL rinvgj(dvde,dedv,4,nvar)
c
               DO 130 iv = 1,nvar
                  dv(iv) = 0.00
                  DO 120 ie = 1,nvar
                     dv(iv) = dv(iv) + dvde(iv,ie)*err(ie)
  120             CONTINUE
                  var(iv) = var(iv) - dv(iv)
  130          CONTINUE
c
               IF (var(1).GT.0.00) THEN
                  WRITE (6,FMT=*) ' WARNING FROM <CORE> E=',var(1)
                  var(1) = -0.20
               END IF
c
               CALL coredir(mrad,var(1),l,xmj,1,vv,bb,rc,dx,
     +                      nmatch,nzero,gc,fc,pow,qow,piw,qiw)
               CALL coredir(mrad,var(1),l,xmj,2,vv,bb,rc,dx,
     +                      nmatch,nzero,gc,fc,pow,qow,piw,qiw)
c
               CALL coreerr(err,var,s,nsol,pow,qow,piw,qiw)
c
               ec = var(1)
c
C---> check relative change in parameters
C---> parameters 3 and 4 = 0 for paramagnetic case !
c
               IF (iter.LT.itermax) THEN
                  IF ((nsol.EQ.2) .AND. .NOT.ferro) THEN
                     DO iv = 3,4
                        err(iv) = 0.00
                        errnew(iv) = 0.00
                        var(iv) = 0.00
                        varnew(iv) = 0.00
                        dv(iv) = 0.00
                     END DO
                  END IF
                  DO 140 iv = 1,nvar
                     IF (abs(var(iv)).EQ.0.00) THEN
                        IF (ferro) THEN
                           WRITE (6,FMT='(A,I3,A)') ' VAR ',iv,
     +                       ' = 0 ??????!!!!!'
                        END IF
                     ELSE
                        IF (abs(dv(iv)/var(iv)).GT.tolvar) GO TO 50
                     END IF
  140             CONTINUE
               ELSE
                  WRITE (6,FMT=
     +'('' ITERATION NOT CONVERGED AFTER'',I3,'' STEPS !'',/,
     + '' PARAMETERS:'',4E18.10,/,'' LAST CORR.:'',4E18.10,/,
     +  '' LAST ERROR:'',4E18.10)') itermax, (var(iv),iv=1,4),
     +              (dv(iv),iv=1,4), (err(ie),ie=1,4)
               END IF
c
C---> end of iterations
c
               CALL cfnorm(mrad,s,t,nsol,nmatch,jtop,var,gc,fc,rc,rc2,
     +                     dx,gck,fck)
C                --------------------------
C--->             CALCULATE  CHARGE DENSITY
C                --------------------------
               CALL ccsdnt(mrad,s,jtop,nsol,l,xmj,kap1,kap2,gck,fck,
     +                     rc2,rchr,rspn)
               DO ir = 1,jtop
                  rhochr(ir) = rhochr(ir) + rchr(ir)
                  rhospn(ir) = rhospn(ir) + rspn(ir)
               END DO
c                -----------------
c--->             HYPERFINE FIELD
c                -----------------
               CALL corehff(mrad,kap1,kap2,xmj,s,nsol,bhf,gck,fck,
     +                      rc,dx,jtop)
               bhff(ilshell) = bhff(ilshell) + bhf
c                 - final info. -
               ectab(ic) = ec + 2*ec_sv
c
               IF (ish.LT.nsh) THEN
                  WRITE (6,FMT=8000) nqn,txtl(l),txtk(iabs(kap(s))),
     +              (2*muem05+1),kap(s),iter,ec/2.
               ELSE
                  WRITE (6,FMT=8000) nqn,txtl(l),txtk(iabs(kap(s))),
     +              (2*muem05+1),kap(s),iter,ec/2.

C              ----------------------------
C--->          CHECK CONSISTENCY OF RESULTS
C              ----------------------------
                  IF (l.NE.0) THEN
                     ic1 = ic - nsh + 1
                     ic2 = ic
                     IF (ectab(ic2).GE.ectab(ic1)) THEN
                        imin = ic1
                        vz = +1.00
                     ELSE
                        imin = ic2
                        vz = -1.00
                     END IF
                     iflag = 0
                     ii = 1
                     DO 150 i = ic1 + 1,ic2,2
                        IF (vz* (ectab(i)-ectab(i-ii)).LT.0.0) iflag = 1
                        ii = 2
  150                CONTINUE
                     IF (ectab(ic1+2).GT.ectab(imin)) iflag = 1
                     DO 160 i = ic1 + 4,ic2 - 1,2
                        IF (ectab(i).GT.ectab(imin)) iflag = 1
                        IF (vz* (ectab(i)-ectab(i-ii)).GT.0.0) iflag = 1
  160                CONTINUE
c
                     IF (ferro .AND. (iflag.EQ.1)) WRITE (6,FMT=
     +'('' >>> CHECK DATA '',                               '' E(KAP,MJ)
     + SHOULD BE MONOTONOUS AND '',                         '' E(+L,MJ)
     +< E(-L-1,MJ) '',//)')
                  END IF
               END IF
c
  170       CONTINUE
c----> s - loop end
  180    CONTINUE
C----> \mu  - loop end
         WRITE (6,FMT='(I4,A1,20X,F16.3)') nqn,txtl(l),bhff(ilshell)
         btot = bhff(ilshell) + btot
  190 CONTINUE
C----> nl  - loop end
      WRITE (6,FMT='(''HFTOT:'',20X,F16.3)') btot
c
      qcor = rsimp(mrad,rhochr,rc,jtop,dx)
      WRITE (6,FMT=8020) 'charge',qcor
      spcor = rsimp(mrad,rhospn,rc,jtop,dx)
      WRITE (6,FMT=8020) 'spin',spcor
c
C----> that's all
c
c
 8000 FORMAT (I4,A1,A3,I3,'/2',2I4,1X,F14.7,F16.3,:,F16.3,/)
 8010 FORMAT (/,'  NQN=',I2,'  L=',I2,'  KAP=',I2,'  MJ=',I2,
     +       '/2    IC=',I3,'  ISH=',I2,/,' E(',A5,')=',F15.5,/,
     +       ' NMATCH  =',I5,'    R=',F10.5,/,' NZERO   =',I5,'    R=',
     +       F10.5,/,' NODES   =',I5,'  RAT=',E11.4)
 8020 FORMAT (' integrated ',A,':',F12.8)
c
      RETURN
      END SUBROUTINE core
      END MODULE m_core
