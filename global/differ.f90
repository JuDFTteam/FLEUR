!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_differ
  USE m_juDFT
  !c
  !cdiffer   subroutine for the differential equations
  !c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
  !c-----x this version modified to accept energies up to emax=2     x----
  !c-----x                              d.d.koelling   7/7/77      x----x
  !c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
CONTAINS
  SUBROUTINE differ(fn,fl,fj,c,z,h,rnot,rn,d,msh,vr,e,a,b,ierr)
    USE m_inwint
    USE m_outint
    IMPLICIT NONE

    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)  :: msh
    INTEGER,INTENT(OUT) :: ierr
    REAL,   INTENT(IN)  :: fj,fl,fn,c,z,h,rnot,rn,d
    REAL,   INTENT(IN)  :: vr(msh)
    REAL,   INTENT(OUT) :: a(msh),b(msh)
    REAL,   INTENT(INOUT) :: e
    !     .. Local Scalars ..
    REAL cis,cs,de,del,dg,eabsv,emax,emin,fkap,g,h3,&
         qcoef,qqqq,r,ra,rb,rg,rj,s,w,wmin
    INTEGER k,ki,kj,nodes,nqnt,ntimes
    LOGICAL dbl
    !     .. Local Arrays ..
    REAL a0(5),b0(5)
    !     .. Intrinsic Functions ..
    INTRINSIC abs,exp,min,sqrt

    a = 0.0
    b = 0.0
    nqnt = fn - fl - 0.99e0
    del = 5.e-7
    emax = 2
    emin = -z*z/fn**2 - 10.0
    IF ((e.GE.emax) .OR. (e.LE.emin)) e = 0.5e0* (emax+emin)
    s      = 2.0e0* (fj-fl)
    cs     = c*s
    cis    = 1.0e0/cs
    fkap   = fj + 0.5e0
    g      = SQRT(fkap**2- (cis*z)**2)
    h3     = h/3.0e0
    ntimes_loop:DO ntimes=1,200
       CALL outint( msh,e,fkap,cs,cis,s,vr,z,rn,rnot,h,d,a0,b0,a,b,ki,nodes)
       IF (nqnt-nodes.GT.0) THEN
          IF (e.GT.emin) emin = e
          e = 0.5e0* (e+emax)
          IF ((emax-emin).GT.del) CYCLE ntimes_loop
          WRITE (6,FMT=8020) nodes,fn,fl,fj,emin,e,emax
          WRITE (6,FMT=8030) vr
          WRITE (6,FMT=8030) a
          CALL juDFT_error("differ 2: problems with solving dirac equation"&
               ,calledby ="differ")
       ENDIF
       IF (nqnt-nodes<0) THEN
          IF (e.LT.emax) emax = e
          e = 0.5e0* (e+emin)
          IF ((emax-emin).GT.del) CYCLE ntimes_loop
          WRITE (6,FMT=8010) nodes,fn,fl,fj,emin,e,emax
          WRITE (6,FMT=8030) vr
          WRITE (6,FMT=8030) a
          CALL juDFT_error("differ 1: problems with solving dirac equation"&
               ,calledby ="differ")
       ENDIF
       !**** correct number of nodes
       ra = a(ki)
       rb = b(ki)
       CALL inwint(e,fl,ki,fkap,cs,cis,s,z,h,d,rn,rnot,msh,vr,a,b,kj)
       ra = ra/a(ki)
       rb = rb/b(ki)
       a(ki:kj) = a(ki:kj)*ra
       b(ki:kj) = ra*b(ki:kj)
       dg = EXP(h*g)
       rg = rnot**g
       DO  k = 1,kj
          a(k) = a(k)*rg
          b(k) = rg*b(k)
          rg   = rg*dg
       ENDDO
       eabsv = ABS(e)
       qcoef = SQRT(eabsv+eabsv)
       rj = rnot*d** (kj-1)
       r = rj
       rg = a(kj)
       dg = b(kj)
       qqqq = MIN(ABS(rg),ABS(dg))
       IF (qqqq.LT.1.0e-25) THEN
          IF (kj+1.LE.SIZE(a)) THEN
             a(kj+1:)=0.0
             b(kj+1:)=0.0
          ENDIF
       ELSE
          wmin = (1.e-35)/qqqq
          DO k=kj+1,SIZE(a)
             r = r*d
             w = EXP(qcoef* (rj-r))
             IF (w.LT.wmin) THEN
                a(k)=0.0
                b(k)=0.0
                EXIT
             ENDIF
             a(k) = w*rg
             b(k) = w*dg
          ENDDO
       ENDIF
       r = rn
       w = r* (a(msh)**2+b(msh)**2)
       r = r + r
       rj = 1.0e0/d
       dbl = .FALSE.
       DO k=msh-1,2,-1
          r = r*rj
          dbl = .NOT. dbl
          rg = r* (a(k)**2+b(k)**2)
          w = w + rg
          IF (dbl) w = w + rg
       ENDDO
       w = h3* (w+rnot* (a(1)**2+b(1)**2))
       de = cs*a(ki)*b(ki)* (ra-rb)/ (ra*w)
       IF (de.GT.0.0e0) emin = e
       IF (de.LT.0.0e0) emax = e
       IF (ABS(de).LT.del) THEN
          e = e+de
          a =  a/SQRT(w)
          b = -b/SQRT(w)
          ierr = 0
          RETURN
       END IF
       e = 0.5e0* (emax+emin)
    END DO ntimes_loop
    WRITE (6,FMT=8000)
    WRITE (6,FMT=8030) fn,fl,fj,emin,emax
    WRITE (6,FMT=8030) e,de
    WRITE (6,FMT=8030) ra,rb,w,a(ki),b(ki)
    WRITE (6,FMT=8000)
    WRITE (6,FMT=8030) vr
    WRITE (6,FMT=8030) a
    ierr = 1
8000 FORMAT (/,/,/,/,10x,' too many tries required')
8010 FORMAT (/,/,/,/,10x,' too many nodes.',i5,3f4.1,3e15.7)
8020 FORMAT (/,/,/,/,10x,' too few nodes. ',i5,3f4.1,3e15.7)
8030 FORMAT (10x,5e14.4)
  END SUBROUTINE differ
END MODULE m_differ
