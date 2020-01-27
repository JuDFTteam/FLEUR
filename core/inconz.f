      SUBROUTINE inconz(e,l,xmj,kap1,kap2,vv,bb,rr,xx1,xx2)

c..........................................................inconz
c initial point for outcome regular solution of dirac eq.
c order kap1=-L-1, kap2=L
c
      USE m_constants, ONLY : c_light
      IMPLICIT NONE
C     .. Scalar Arguments ..
      REAL bb,e,rr,vv,xmj
      INTEGER kap1,kap2,l
C     ..
C     .. Array Arguments ..
      REAL xx1(4),xx2(4)
C     ..
C     .. Local Scalars ..
      REAL aa11,aa12,aa21,aa22,bb1,bb2,bc0,bqq,cc,cg1,cg2,cg4,cg5,cg8,
     +     cgo,csq,det,emvpp,emvqq,rpwgpm,tz,vc0
      INTEGER i,j,m,mps,nsol
C     ..
C     .. Local Arrays ..
      REAL cgd(2),cgmd(2),gam(2),kap(2),pc(2,2,0:1),qc(2,2,0:1),wp(2,2),
     +     wq(2,2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,real,int,sqrt
C     ..
      cc = c_light(2.0)
      csq = cc*cc
c
C     EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
C VV=VV(1)
      tz = real(int(-vv*rr))
      vc0 = vv - (-tz)/rr
C BB=BB(1)
      bc0 = bb
C
C    CALCULATE G-COEFFICIENTS OF B-FIELD
C
c      KAP1 = - L - 1
c      KAP2 = + L
      cg1 = -xmj/ (kap1+0.5)
      cg5 = -xmj/ (-kap1+0.5)
      cgd(1) = cg1
      cgmd(1) = cg5
      kap(1) = real(kap1)
      gam(1) = sqrt(kap(1)**2- (tz/cc)**2)
      IF (abs(xmj).GE.l) THEN
         cg2 = 0.00
         cg4 = 0.00
         cg8 = 0.00
         nsol = 1
         cgd(2) = 0.00
         cgo = 0.00
         cgmd(2) = 0.00
         gam(2) = 0.00
         kap(2) = 0.00
      ELSE
         cg2 = -sqrt(1.0- (xmj/ (kap1+0.50))**2)
         cg4 = -xmj/ (kap2+0.50)
         cg8 = -xmj/ (-kap2+0.50)
         nsol = 2
         cgd(2) = cg4
         cgo = cg2
         cgmd(2) = cg8
         kap(2) = real(kap2)
         gam(2) = sqrt(kap(2)**2- (tz/cc)**2)
      END IF
C
      DO 10 j = 1,nsol
         i = 3 - j
         pc(j,j,0) = sqrt(abs(kap(j)-gam(j)))
         qc(j,j,0) = (kap(j)+gam(j))* (csq/tz)*pc(j,j,0)
         pc(i,j,0) = 0.0
         qc(i,j,0) = 0.0
   10 CONTINUE
C  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS
      mps = 1
      aa12 = -tz/csq
      aa21 = tz
      emvqq = (e-vc0+csq)/csq
      emvpp = -e + vc0
      bqq = bc0/csq
      DO 40 j = 1,nsol
         DO 30 m = 1,mps
            DO 20 i = 1,nsol
               bb1 = (emvqq+bqq*cgmd(i))*qc(i,j,m-1)
               bb2 = (emvpp+bc0*cgd(i))*pc(i,j,m-1) +
     +               bc0*cgo*pc(3-i,j,m-1)
               aa11 = gam(j) + m + kap(i)
               aa22 = gam(j) + m - kap(i)
               det = aa11*aa22 - aa12*aa21
               pc(i,j,m) = (bb1*aa22-aa12*bb2)/det
               qc(i,j,m) = (aa11*bb2-bb1*aa21)/det
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
C
C  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
C  FOR THE FIRST - MESH - POINT
c      RR= RC(1)
      DO 80 j = 1,nsol
         rpwgpm = rr** (gam(j))
         DO 50 i = 1,nsol
            wp(i,j) = pc(i,j,0)*rpwgpm
            wq(i,j) = qc(i,j,0)*rpwgpm
   50    CONTINUE
         DO 70 m = 1,mps
            rpwgpm = rpwgpm*rr
            DO 60 i = 1,nsol
               wp(i,j) = wp(i,j) + pc(i,j,m)*rpwgpm
               wq(i,j) = wq(i,j) + qc(i,j,m)*rpwgpm
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
C---> First point solutions construction
      IF (nsol.EQ.2) THEN
         xx1(1) = wp(1,1)
         xx1(2) = wq(1,1)
         xx1(3) = wp(2,1)
         xx1(4) = wq(2,1)
c
         xx2(1) = wp(1,2)
         xx2(2) = wq(1,2)
         xx2(3) = wp(2,2)
         xx2(4) = wq(2,2)
      ELSE
         xx1(1) = wp(1,1)
         xx1(2) = wq(1,1)
c not needed
         xx1(3) = 0.0
         xx1(4) = 0.0
         xx2(1) = 0.0
         xx2(2) = 0.0
         xx2(3) = 0.0
         xx2(4) = 0.0
      END IF
      RETURN
      END
