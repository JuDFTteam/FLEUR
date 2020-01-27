      MODULE m_differ
      use m_juDFT
c
cdiffer   subroutine for the differential equations
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
c-----x this version modified to accept energies up to emax=2     x----
c-----x                              d.d.koelling   7/7/77      x----x
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
      CONTAINS
      SUBROUTINE differ(
     >                  fn,fl,fj,c,z,h,rnot,rn,d,msh,vr,
     X                  e,
     <                  a,b,ierr)
      USE m_inwint
      USE m_outint

      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER,INTENT(IN)  :: msh
      INTEGER,INTENT(OUT) :: ierr
      REAL,   INTENT(IN)  :: fj,fl,fn,c,z,h,rnot,rn,d
      REAL,   INTENT(IN)  :: vr(msh)
      REAL,   INTENT(OUT) :: a(msh),b(msh)
      REAL,   INTENT(INOUT) :: e
C     ..
C     .. Local Scalars ..
      REAL cis,cs,de,del,dg,eabsv,emax,emin,fkap,g,h3,
     +     qcoef,qqqq,r,ra,rb,rg,rj,s,w,wmin
      INTEGER k,ki,kj,n,nodes,nqnt,ntimes
      LOGICAL dbl
C     ..
C     .. Local Arrays ..
      REAL a0(5),b0(5)
C     ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,min,sqrt
C     ..
      a = 0.0
      b = 0.0
      ierr = 0
      nqnt = fn - fl - 0.99e0
      n = msh
      del = 5.e-7
      emax = 2
      emin = -z*z/fn**2 - 10.0
      IF ((e.GE.emax) .OR. (e.LE.emin)) e = 0.5e0* (emax+emin)
      s      = 2.0e0* (fj-fl)
      cs     = c*s
      cis    = 1.0e0/cs
      fkap   = fj + 0.5e0
      g      = sqrt(fkap**2- (cis*z)**2)
      h3     = h/3.0e0
      ntimes = 0
   10 ntimes = ntimes + 1
      IF (ntimes.GT.200) GO TO 140
      CALL outint(
     >            msh,e,fkap,cs,cis,s,vr,z,rn,rnot,h,d,
     <            a0,b0,a,b,
     <            ki,nodes)
      IF (nqnt-nodes.EQ.0) GOTO 30
      IF (nqnt-nodes.GT.0) GOTO 130
c**** too many nodes
   20 IF (e.LT.emax) emax = e
      e = 0.5e0* (e+emin)
      IF ((emax-emin).GT.del) GO TO 10
      WRITE (6,FMT=8010) nodes,fn,fl,fj,emin,e,emax
      WRITE (6,FMT=8030) vr
      WRITE (6,FMT=8030) a
      CALL juDFT_error("differ 1: problems with solving dirac equation"
     +     ,calledby ="differ")
c**** correct number of nodes
   30 ra = a(ki)
      rb = b(ki)
      CALL inwint(
     >            e,fl,ki,fkap,cs,cis,s,z,h,d,rn,rnot,msh,vr,
     X            a,b,
     <            kj)
      ra = ra/a(ki)
      rb = rb/b(ki)
      DO 40 k = ki,kj
         a(k) = a(k)*ra
         b(k) = ra*b(k)
   40 CONTINUE
      dg = exp(h*g)
      rg = rnot**g
      DO 50 k = 1,kj
         a(k) = a(k)*rg
         b(k) = rg*b(k)
         rg   = rg*dg
   50 CONTINUE
      eabsv = abs(e)
      qcoef = sqrt(eabsv+eabsv)
      rj = rnot*d** (kj-1)
      r = rj
      rg = a(kj)
      dg = b(kj)
      qqqq = min(abs(rg),abs(dg))
      IF (qqqq.LT.1.0e-25) GO TO 70
      wmin = (1.e-35)/qqqq
      k = kj
   60 k = k + 1
      IF (k.GT.n) GO TO 90
      r = r*d
      w = exp(qcoef* (rj-r))
      IF (w.LT.wmin) GO TO 80
      a(k) = w*rg
      b(k) = w*dg
      GO TO 60
   70 k = kj + 1
   80 IF (k.GT.n) GO TO 90
      a(k) = 0.0e0
      b(k) = 0.0e0
      k = k + 1
      GO TO 80
   90 r = rn
      w = r* (a(n)**2+b(n)**2)
      k = n
      r = r + r
      rj = 1.0e0/d
      dbl = .false.
  100 k = k - 1
      r = r*rj
      dbl = .NOT. dbl
      rg = r* (a(k)**2+b(k)**2)
      w = w + rg
      IF (dbl) w = w + rg
      IF (k.GT.2) GO TO 100
      w = h3* (w+rnot* (a(1)**2+b(1)**2))
      de = cs*a(ki)*b(ki)* (ra-rb)/ (ra*w)
      IF (de.GT.0.0e0) emin = e
      IF (de.LT.0.0e0) emax = e
      e = e + de
      IF (abs(de).LT.del) GO TO 110
      IF ((e.GE.emax) .OR. (e.LE.emin)) e = 0.5e0* (emax+emin)
      GO TO 10
  110 w = 1.0e0/sqrt(w)
      !
      ! for a consistent definition of the small component in
      ! the Dirac and scalar-relativistic approximation (SRA)
      ! we have to multiply the small component with -1;
      ! then the small component of the Dirac equation
      ! also fulfills the SRA for s-states
      !                                               M.B. (July, 2012)
      DO 120 k = 1,n
         a(k) =  a(k)*w
         b(k) = -b(k)*w
  120 CONTINUE
      ierr = 0
      RETURN
c**** too few nodes
  130 IF (e.GT.emin) emin = e
      e = 0.5e0* (e+emax)
      IF ((emax-emin).GT.del) GO TO 10
      WRITE (6,FMT=8020) nodes,fn,fl,fj,emin,e,emax
      WRITE (6,FMT=8030) vr
      WRITE (6,FMT=8030) a
      CALL juDFT_error("differ 2: problems with solving dirac equation"
     +     ,calledby ="differ")
  140 WRITE (6,FMT=8000)
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
