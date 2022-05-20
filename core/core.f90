MODULE m_core

   CONTAINS
   
   SUBROUTINE core(mrad,vt,bt,zz,stval,dx,nlshell,nqntab,lqntab,jtop,ectab,rhochr,rhospn)

      USE m_constants
      USE m_felim
      USE m_findlim
      USE m_cnodes
      USE m_coredir
      USE m_ccsdnt

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mrad,jtop,nlshell
      REAL,    INTENT (IN) :: dx,stval,zz
      INTEGER, INTENT (IN) :: lqntab(15),nqntab(15)
      REAL,    INTENT (IN) :: bt(mrad),vt(mrad)
      REAL,    INTENT (OUT) :: rhochr(mrad),rhospn(mrad)
      REAL,    INTENT (INOUT) :: ectab(100)

      REAL    :: bhf,bsum,btot,dvstep,ec,elim,qcor,spcor,tolvar,vz,xmj,ec_sv
      INTEGER :: i,ic,ic1,ic2,ie,iflag,ii,ilshell,imin,ir,ish,istart,iter
      INTEGER :: itermax,iv,j,jv,kap1,kap2,l,lll,muem05,n,ncor,nmatch,node
      INTEGER :: nqn,nrmax,nsh,nsol,nval,nvar,nzero,s,t
      LOGICAL :: ferro

      REAL :: bb(mrad),bhff(15),dedv(4,4),dv(4),dvde(4,4),err(4),errnew(4)
      REAL :: fc(2,2,mrad),fck(2,2,mrad),gc(2,2,mrad),gck(2,2,mrad)
      REAL :: piw(2,2),pow(2,2),qiw(2,2),qow(2,2),rc(mrad),rc2(mrad)
      REAL :: rchr(mrad),rspn(mrad),var(4),varnew(4),vv(mrad)
      INTEGER kap(2)
      CHARACTER txtl(0:3)
      CHARACTER*3 txtk(4)

      REAL rsimp
      EXTERNAL rsimp

      EXTERNAL cfnorm,coreerr,corehff,nwrfst,rinvgj

      INTRINSIC abs,exp,iabs,int,sign

      DATA txtl/'s','p','d','f'/
      DATA txtk/'1/2','3/2','5/2','7/2'/
!      DATA nqntab/1,2,2,3,3,3,4,4,4,4,5,5,5,6,6/
!      DATA lqntab/0,0,1,0,1,2,0,1,2,3,0,1,2,0,1/
      DATA tolvar/1.0E-12/
      DATA dvstep/0.010/
      DATA itermax/150/

      nrmax = mrad
      DO n = 1,nrmax
         rc(n) = exp(stval+ (n-1)*dx)
         rc2(n) = rc(n)**2
      END DO

      DO n = 1,nrmax
         rhochr(n) = 0.00
         rhospn(n) = 0.00
      END DO

      bsum = 0.0
      DO n = 1,jtop
         vv(n) = vt(n)
         bb(n) = bt(n)
         bsum = bsum + abs(bb(n))
      END DO

      ncor = 0
      DO ilshell = 1,nlshell
         ncor = ncor + 2* (2*lqntab(ilshell)+1)
      END DO
      nval = int(zz) - ncor

      IF (bsum.GT.1.0E-8) THEN
         ferro = .true.
      ELSE
         ferro = .false.
         WRITE (oUnit,FMT='(''  PARAMAGNETIC CASE'',/,''  *****************'')')
      END IF

      WRITE (oUnit,FMT='(/,''  ATOMIC NUMBER  : '',F6.2,/,''  CORE ELECTRONS : '',I5,/,''  VAL. ELECTRONS : '',I5)') zz,ncor,nval
      WRITE (oUnit,FMT='(  /,''  MESH    RC(1)  : '',F12.6,//10x,''MUE  KAP ITER    ENERGY '',''(Hart)'', " HF (KG) ")') rc(1)


      ! ---------------------------------------
      ! INITIALIZE QUANTUM NUMBERS  NQN  AND  L
      ! ---------------------------------------
      
      btot = 0.0
      ic = 0
      
      DO ilshell = 1,nlshell ! nl  - loop
         nqn = nqntab(ilshell)
         l = lqntab(ilshell)
         nsh = 2* (2*l+1)
         bhff(ilshell) = 0.0
         ish = 0
         
         DO muem05 = -l - 1,+l ! \mu  - loop
            xmj = muem05 + 0.5
            kap1 = -l - 1
            kap2 = l
            kap(1) = kap1
            kap(2) = kap2

            lll = l* (l+1)
            IF (abs(xmj).GT.l) THEN
               nsol = 1
            ELSE
               nsol = 2
            END IF

            nvar = 2*nsol

            DO s = 1,nsol ! s - loop : solutions for each (nl,\mu)
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

               ic = ic + 1
               ish = ish + 1
               t = 3 - s
               ec = ectab(ic)
               write(*,*) ic,ec

               IF (ic > 1) THEN
                 IF (ectab(ic-1) > 0.1) THEN
                   vv(:) = vv(:) + 2*ec_sv
                   bb(:) = bb(:) + 2*ec_sv
                 ENDIF
               ENDIF
               IF (ec.GT.0.1) THEN
                  vv(:) = vv(:) - 2*ectab(ic)
                  bb(:) = bb(:) - 2*ectab(ic)
                  ec_sv = ectab(ic)
                  ec = -ectab(ic)
               ELSE
                 ec_sv = 0.0
               ENDIF
               istart = 1

               IF (ish.GT.1) GO TO 30

               ! energy limits : elim and emax
               CALL felim(mrad,lll,zz,nqn,vv,rc,elim)
               write(*,*) 'elim',elim
               
   20          IF (ec.LE.elim) ec = elim*0.70

               ! turning point and physical infinity : nmatch and nzero
               CALL findlim(mrad,lll,ec,vv,rc,nmatch,nzero)
               write(*,*) 'nmatch,nzero',nmatch,nzero
   30          CONTINUE

               ! nodes for large component : if IFLAG=0 number of nodes is OK
               iflag = 0
               CALL cnodes(mrad,iflag,s,ec,l,xmj,nqn,vv,bb,rc,dx,nmatch,nzero,gc,fc,pow,qow,piw,qiw,node)
               IF (iflag.NE.0) GO TO 20

               ! setup of Newton-Raphson algorithm
               CALL nwrfst(mrad,nsol,s,t,nmatch,nzero,ferro,ec,rc,pow,piw,gc,err,var,dv,varnew,errnew)

               CALL coreerr(err,var,s,nsol,pow,qow,piw,qiw)
               DO iv = 1,nvar
                  dv(iv) = var(iv)
               END DO

               iter = 0
               ! iterations start
   50          iter = iter + 1
   
               !----------------------------------
               ! CHECK WHETHER NUMBER OF NODES O.K.
               !----------------------------------
               IF (iter.GT.1) THEN
                  node = 0
                  DO n = 2, (nmatch-1)
                     IF (gc(s,s,n)*gc(s,s,n-1).LT.0.0) node = node + 1
                  END DO

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
               
               ! - atomic energy search -
               DO iv = 2,nvar
                  DO jv = 1,nvar
                     varnew(jv) = var(jv)
                  END DO
                  varnew(iv) = var(iv) + dv(iv)*dvstep

                  IF (abs(var(iv)).EQ.0.00) THEN
                     IF (ferro) THEN
                        WRITE (oUnit,FMT='(A,I3,A)') ' VAR ',iv,' = 0 ??????!!!!!'
                     END IF
                  ELSE
                     IF (abs(dv(iv)/var(iv)).LT.tolvar) varnew(iv) = var(iv)*(1.00+sign(dvstep*tolvar,dv(iv)))
                  END IF

                  CALL coreerr(errnew,varnew,s,nsol,pow,qow,piw,qiw)

                  DO ie = 1,nvar
                     IF ((errnew(ie)-err(ie)).EQ.0.00) THEN
                        dedv(ie,iv) = 0.00
                        IF ((ie.EQ.iv) .AND. .NOT.ferro) THEN
                           dedv(ie,iv) = 1.00
                        END IF
                     ELSE
                        dedv(ie,iv) = (errnew(ie)-err(ie)) / (varnew(iv)-var(iv))
                     END IF
                  END DO
               END DO

               DO jv = 1,nvar
                  varnew(jv) = var(jv)
               END DO
               varnew(1) = var(1) + dv(1)*dvstep
               IF (abs(dv(1)/var(1)).LT.tolvar) varnew(1) = var(1) * (1.00+sign(dvstep*tolvar,dv(1)))

               CALL coredir(mrad,varnew(1),l,xmj,1,vv,bb,rc,dx,nmatch,nzero,gc,fc,pow,qow,piw,qiw)
               CALL coredir(mrad,varnew(1),l,xmj,2,vv,bb,rc,dx,nmatch,nzero,gc,fc,pow,qow,piw,qiw)

               CALL coreerr(errnew,varnew,s,nsol,pow,qow,piw,qiw)

               DO ie = 1,nvar
                  dedv(ie,1) = (errnew(ie)-err(ie))/ (varnew(1)-var(1))
               END DO

               CALL rinvgj(dvde,dedv,4,nvar)

               DO iv = 1,nvar
                  dv(iv) = 0.00
                  DO ie = 1,nvar
                     dv(iv) = dv(iv) + dvde(iv,ie)*err(ie)
                  END DO
                  var(iv) = var(iv) - dv(iv)
               END DO

               IF (var(1).GT.0.00) THEN
                  WRITE (oUnit,FMT=*) ' WARNING FROM <CORE> E=',var(1)
                  var(1) = -0.20
               END IF

               CALL coredir(mrad,var(1),l,xmj,1,vv,bb,rc,dx,nmatch,nzero,gc,fc,pow,qow,piw,qiw)
               CALL coredir(mrad,var(1),l,xmj,2,vv,bb,rc,dx,nmatch,nzero,gc,fc,pow,qow,piw,qiw)

               CALL coreerr(err,var,s,nsol,pow,qow,piw,qiw)

               ec = var(1)

               ! check relative change in parameters
               ! parameters 3 and 4 = 0 for paramagnetic case !
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
                  DO iv = 1,nvar
                     IF (abs(var(iv)).EQ.0.00) THEN
                        IF (ferro) THEN
                           WRITE (oUnit,FMT='(A,I3,A)') ' VAR ',iv,' = 0 ??????!!!!!'
                        END IF
                     ELSE
                        IF (abs(dv(iv)/var(iv)).GT.tolvar) GO TO 50
                     END IF
                  END DO
               ELSE
                  WRITE (oUnit, FMT='('' ITERATION NOT CONVERGED AFTER'',I3,'' STEPS !'',/,'' PARAMETERS:'',4E18.10,/,'' LAST CORR.:'',4E18.10,/,'' LAST ERROR:'',4E18.10)') &
                     itermax, (var(iv),iv=1,4),(dv(iv),iv=1,4), (err(ie),ie=1,4)
               END IF

               ! end of iterations
               CALL cfnorm(mrad,s,t,nsol,nmatch,jtop,var,gc,fc,rc,rc2,dx,gck,fck)
               
               !--------------------------
               ! CALCULATE  CHARGE DENSITY
               !--------------------------
               CALL ccsdnt(mrad,s,jtop,nsol,l,xmj,kap1,kap2,gck,fck,rc2,rchr,rspn)
               DO ir = 1,jtop
                  rhochr(ir) = rhochr(ir) + rchr(ir)
                  rhospn(ir) = rhospn(ir) + rspn(ir)
               END DO
               
               !-----------------
               ! HYPERFINE FIELD
               !-----------------
               CALL corehff(mrad,kap1,kap2,xmj,s,nsol,bhf,gck,fck,rc,dx,jtop)
               bhff(ilshell) = bhff(ilshell) + bhf
               
               ! final info
               ectab(ic) = ec + 2*ec_sv

               IF (ish.LT.nsh) THEN
                  WRITE (oUnit,FMT=8000) nqn,txtl(l),txtk(iabs(kap(s))),(2*muem05+1),kap(s),iter,ec/2.
               ELSE
                  WRITE (oUnit,FMT=8000) nqn,txtl(l),txtk(iabs(kap(s))),(2*muem05+1),kap(s),iter,ec/2.

                  !----------------------------
                  ! CHECK CONSISTENCY OF RESULTS
                  !----------------------------
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
                     DO i = ic1 + 1,ic2,2
                        IF (vz* (ectab(i)-ectab(i-ii)).LT.0.0) iflag = 1
                        ii = 2
                     END DO
                     IF (ectab(ic1+2).GT.ectab(imin)) iflag = 1
                     DO i = ic1 + 4,ic2 - 1,2
                        IF (ectab(i).GT.ectab(imin)) iflag = 1
                        IF (vz* (ectab(i)-ectab(i-ii)).GT.0.0) iflag = 1
                     END DO

                     IF (ferro .AND. (iflag.EQ.1)) THEN
                        WRITE (oUnit,FMT='('' >>> CHECK DATA '',                               '' E(KAP,MJ) SHOULD BE MONOTONOUS AND '',                         '' E(+L,MJ) < E(-L-1,MJ) '',//)')
                     END IF
                  END IF
               END IF

            END DO ! s - loop end
            
         END DO ! \mu  - loop end
         
         WRITE (oUnit,FMT='(I4,A1,20X,F16.3)') nqn,txtl(l),bhff(ilshell)
         btot = bhff(ilshell) + btot
         
      END DO ! nl  - loop end
      
      WRITE (oUnit,FMT='(''HFTOT:'',20X,F16.3)') btot

      qcor = rsimp(mrad,rhochr,rc,jtop,dx)
      WRITE (oUnit,FMT=8020) 'charge',qcor
      spcor = rsimp(mrad,rhospn,rc,jtop,dx)
      WRITE (oUnit,FMT=8020) 'spin',spcor

      ! that's all


 8000 FORMAT (I4,A1,A3,I3,'/2',2I4,1X,F14.7,F16.3,:,F16.3,/)
 8010 FORMAT (/,'  NQN=',I2,'  L=',I2,'  KAP=',I2,'  MJ=',I2,'/2    IC=',I3,'  ISH=',I2,/,' E(',A5,')=',F15.5,/,' NMATCH  =',I5,'    R=',F10.5,/,' NZERO   =',I5,'    R=',F10.5,/,' NODES   =',I5,'  RAT=',E11.4)
 8020 FORMAT (' integrated ',A,':',F12.8)

   END SUBROUTINE core
END MODULE m_core
