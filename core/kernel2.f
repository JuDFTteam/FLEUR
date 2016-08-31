      MODULE m_kernel2
C.............................................................kernel2
c solution of 4 coupled equations (logariphmic mesh) for (l,\mu)
c 1->k1, 2->k2
c d/dxP1=[-k1]P1+exp(x)[e-v+2mc^2]/c^2Q1+exp(x)B/c^2\sigma_z[-1,-1]Q1
c d/dxQ1=[k1]Q1-exp(x)[e-v]P1+exp(x)B[\sigma_z[1,1]P1+\sigma_z[1,2]P2]
c d/dxP2=[-k2]P2+exp(x)[e-v+2mc^2]/c^2Q2+exp(x)B/c^2\sigma_z[-1,-1]Q2
c d/dxQ2=[k2]Q2-exp(x)[e-v]P2+exp(x)B[\sigma_z[2,1]P1+\sigma_z[2,2]P2]
c notations:
c Pi=exp(x)*gi
c Qi=exp(x)*cfi
c \gamma=sqrt(k^2-(Z/c)^2)
c v=0.5[V(1/2)+V(-1/2)]
c b=0.5[V(1/2)-V(-1/2)]
c Rydberg units: in charge
c Hartree units: com.
c NSOL= 2 -  4 equations
c ------------                                     a. shick KFA 1996
      CONTAINS
      SUBROUTINE kernel2(mrad,nsol,xmj,k1,k2,xx1,xx2,e,v,b,ri,dx,
     +                   nmatch,nstart,dp,dq,wp,wq)


      USE m_constants, ONLY : c_light
      USE m_diff
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad
      REAL dx,e,xmj
      INTEGER k1,k2,nmatch,nsol,nstart
C     ..
C     .. Array Arguments ..
      REAL b(mrad),dp(2,2,mrad),dq(2,2,mrad),ri(mrad),v(mrad),
     +     wp(2,2,mrad),wq(2,2,mrad),xx1(4),xx2(4)
C     ..
C     .. Local Scalars ..
      REAL abz,bh,bmn,bmnp1,bn,bnp1,c1,c2,c3,cc,csq,dp11,dp12,
     +     dp13,dp14,dp21,dp22,dp23,dp24,dq11,dq12,dq13,dq14,dq21,dq22,
     +     dq23,dq24,dxd8,expdxh,gk,gk1,gk1k,gkk1,gmk,gmk1,p1c,p2c,
     +     q1c,q2c,r,vh,vmn,vmnp1,vn,vnp1,xk,xk1,xm
      INTEGER i,iny,ir,j,jri,n,nrk
C     ..
C     .. Local Arrays ..
      REAL bm(mrad),p1(mrad),p1s(mrad),p2(mrad),p2s(mrad),q1(mrad),
     +     q1s(mrad),q2(mrad),q2s(mrad),vm(mrad)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,real,sqrt
C     ..
C     .. Data statements ..
      DATA abz/0.5/
C     ..
      cc = c_light(2.0)
c
      DO ir = 1,mrad
         p1(ir) = 0.0
         q1(ir) = 0.0
         p1s(ir) = 0.0
         q1s(ir) = 0.0
c
         p2(ir) = 0.0
         q2(ir) = 0.0
         p2s(ir) = 0.0
         q2s(ir) = 0.0
c
         DO i = 1,2
            DO j = 1,2
c
               wp(i,j,ir) = 0.0
               wq(i,j,ir) = 0.0
c
               dp(i,j,ir) = 0.0
               dq(i,j,ir) = 0.0
            END DO
         END DO
      END DO
C
      csq = cc*cc
      expdxh = exp(dx/2.0)
      dxd8 = dx/8.0
      jri = nmatch
      nrk = jri
c-derivatives of potential
      CALL diff(
     >          mrad,v,dx,jri,
     <          vm)
      CALL diff(
     >          mrad,b,dx,jri,
     <          bm)
c- begin to solve coupled Dirac equation in "magnetic" field
      xk = real(k1)
      xk1 = real(k2)
      xm = xmj
      gk = -xm/ (xk+abz)
      gk1 = -xm/ (xk1+abz)
      gkk1 = -sqrt(1.0-gk*gk)
      gk1k = gkk1
      gmk = xm/ (xk-abz)
      gmk1 = xm/ (xk1-abz)
c
c     NB: V=R**2*POTENTIAL
c     NB: B=R**2*FIELD
c
      DO 30 iny = 1,nsol
C                                                                       WAB02140
         r = ri(nstart)
         vn = v(nstart)
         vmn = vm(nstart)
         bn = b(nstart)
         bmn = bm(nstart)
         IF (iny.EQ.1) THEN
            p1(nstart) = xx1(1)
            q1(nstart) = xx1(2)
            p2(nstart) = xx1(3)
            q2(nstart) = xx1(4)
         ELSE
            p1(nstart) = xx2(1)
            q1(nstart) = xx2(2)
            p2(nstart) = xx2(3)
            q2(nstart) = xx2(4)
         END IF
c---->f.po solution of Dirac eq.
         c1 = vn/r**2 - e
         c2 = 1.0 - c1/csq
c      C2=2.0-C1/CSQ !H
         c3 = bn/r/r
         p1s(nstart) = -xk*p1(nstart) + r* (c2+c3/csq*gmk)*q1(nstart)
         q1s(nstart) = xk*q1(nstart) + r*
     +                 ((c1+c3*gk)*p1(nstart)+c3*gkk1*p2(nstart))
C
         p2s(nstart) = -xk1*p2(nstart) + r* (c2+c3/csq*gmk1)*q2(nstart)
         q2s(nstart) = xk1*q2(nstart) +
     +                 r* ((c1+c3*gk1)*p2(nstart)+c3*gk1k*p1(nstart))
C   ******************************************************************
C   *                                                                *
C   *     THE RUNGE-KUTTA METHOD IS USED IN THE <NKR> FIRST STEPS    *
C   *                                                                *
C   ******************************************************************
C
         n = nstart
c-> y1
   10    p1c = p1(n)
         q1c = q1(n)
         p2c = p2(n)
         q2c = q2(n)
c-> k1=f(x1,y1)
         dp11 = dx* (-xk*p1c+r* (c2+c3/csq*gmk)*q1c)
         dq11 = dx* (xk*q1c+r* ((c1+c3*gk)*p1c+c3*gkk1*p2c))
         dp21 = dx* (-xk1*p2c+r* (c2+c3/csq*gmk1)*q2c)
         dq21 = dx* (xk1*q2c+r* ((c1+c3*gk1)*p2c+c3*gk1k*p1c))
c-> y11=y1+0.5*k1
         p1c = p1c + 0.5*dp11
         q1c = q1c + 0.5*dq11
         p2c = p2c + 0.5*dp21
         q2c = q2c + 0.5*dq21
c-> x1'=x1+0.5h
         r = r*expdxh
c-> interpolation of V and B for x1'
         vnp1 = v(n+1)
         vmnp1 = vm(n+1)
         bnp1 = b(n+1)
         bmnp1 = bm(n+1)
         vh = (vn+vnp1)*0.5 + (vmn-vmnp1)*dxd8
         bh = (bn+bnp1)*0.5 + (bmn-bmnp1)*dxd8
c-> k2=f(x1',y11)
         c1 = vh/r/r - e
         c2 = 1.0 - c1/csq
c      C2=2.0-C1/CSQ !H
         c3 = bh/r/r
C                                                                       WAB02850
         dp12 = dx* (-xk*p1c+r* (c2+c3/csq*gmk)*q1c)
         dq12 = dx* (xk*q1c+r* ((c1+c3*gk)*p1c+c3*gkk1*p2c))
         dp22 = dx* (-xk1*p2c+r* (c2+c3/csq*gmk1)*q2c)
         dq22 = dx* (xk1*q2c+r* ((c1+c3*gk1)*p2c+c3*gk1k*p1c))
c-> y12=y1+0.5*k2 (y12=y11-0.5*k1+0.5*k2)
         p1c = p1c + 0.5* (dp12-dp11)
         q1c = q1c + 0.5* (dq12-dq11)
         p2c = p2c + 0.5* (dp22-dp21)
         q2c = q2c + 0.5* (dq22-dq21)
c-> k3=f(x1',y12)
         dp13 = dx* (-xk*p1c+r* (c2+c3/csq*gmk)*q1c)
         dq13 = dx* (xk*q1c+r* ((c1+c3*gk)*p1c+c3*gkk1*p2c))
         dp23 = dx* (-xk1*p2c+r* (c2+c3/csq*gmk1)*q2c)
         dq23 = dx* (xk1*q2c+r* ((c1+c3*gk1)*p2c+c3*gk1k*p1c))
c-> y13=y1+k3 (y13=y12-0.5*k2+k3)
         p1c = p1c + dp13 - 0.5*dp12
         q1c = q1c - 0.5*dq12 + dq13
         p2c = p2c + dp23 - 0.5*dp22
         q2c = q2c - 0.5*dq22 + dq23
c-> x2=x1+h
         n = n + 1
         r = ri(n)
c-> k4=f(x2,y13)
         c1 = vnp1/r/r - e
         c2 = 1.0 - c1/csq
c      C2=2.0-C1/CSQ !H
         c3 = bnp1/r/r
c
         dp14 = dx* (-xk*p1c+r* (c2+c3/csq*gmk)*q1c)
         dq14 = dx* (xk*q1c+r* ((c1+c3*gk)*p1c+c3*gkk1*p2c))
         dp24 = dx* (-xk1*p2c+r* (c2+c3/csq*gmk1)*q2c)
         dq24 = dx* (xk1*q2c+r* ((c1+c3*gk1)*p2c+c3*gk1k*p1c))
c-> y2=y1+1/6*(k1+2*k2+2*k3+k4)
         p1(n) = p1(n-1) + (dp11+2.0* (dp12+dp13)+dp14)/6.0
         q1(n) = q1(n-1) + (dq11+2.0* (dq12+dq13)+dq14)/6.0
         p2(n) = p2(n-1) + (dp21+2.0* (dp22+dp23)+dp24)/6.0
         q2(n) = q2(n-1) + (dq21+2.0* (dq22+dq23)+dq24)/6.0
c-> f(x2,y2)
         p1s(n) = -xk*p1(n) + r* (c2+c3/csq*gmk)*q1(n)
         q1s(n) = xk*q1(n) + r* ((c1+c3*gk)*p1(n)+c3*gkk1*p2(n))
         p2s(n) = -xk1*p2(n) + r* (c2+c3/csq*gmk1)*q2(n)
         q2s(n) = xk1*q2(n) + r* ((c1+c3*gk1)*p2(n)+c3*gk1k*p1(n))
c-> redefinition of V ind B (2->1)
         vn = vnp1
         bn = bnp1
         vmn = vmnp1
         bmn = bmnp1
C                                                                       WAB03340
         IF (n-nrk.LT.0) GOTO 10
c---->Milne's method off now!
c
c storage of wave functions and their derivatives.
         DO ir = 1,jri
            wp(1,iny,ir) = p1(ir)
            wp(2,iny,ir) = p2(ir)
            wq(1,iny,ir) = q1(ir)
            wq(2,iny,ir) = q2(ir)
c derivatives
            dp(1,iny,ir) = p1s(ir)
            dp(2,iny,ir) = p2s(ir)
            dq(1,iny,ir) = q1s(ir)
            dq(2,iny,ir) = q2s(ir)
         END DO
C***********************************************************            WAB05210
   30 CONTINUE
c      write(6,*) wp
      
      END SUBROUTINE kernel2
      END MODULE m_kernel2
