      MODULE m_coredir
      CONTAINS
      SUBROUTINE coredir(mrad,e,l,xmj,iway,vv,bb,rc,dx,nmatch,nzero,
     +                   gc,fc,pow,qow,piw,qiw)

c.........................................................coredir
c   solution of dirac equation for atomic problem
c   full relativistic spin-polarized case
c   Ry units: in charge
c_______________________________________________  a. shick KFA 1996

      USE m_constants, ONLY : c_light
      USE m_crtail
      USE m_kernel1
      USE m_kernel2

      IMPLICIT NONE
c
C     .. Parameters ..
      INTEGER, INTENT (IN) :: mrad
C     ..
C     .. Scalar Arguments ..
      REAL dx,e,xmj
      INTEGER iway,l,nmatch,nzero
C     ..
C     .. Array Arguments ..
      REAL bb(mrad),fc(2,2,mrad),gc(2,2,mrad),piw(2,2),pow(2,2),
     +     qiw(2,2),qow(2,2),rc(mrad),vv(mrad)
C     ..
C     .. Local Scalars ..
      REAL cc,csq,dx1
      INTEGER i,ir,irv,j,kap1,kap2,n,nn,nsol,nstart
C     ..
C     .. Local Arrays ..
      REAL b(mrad),dp(2,2,mrad),dq(2,2,mrad),ra(mrad),v(mrad),
     +     wp(2,2,mrad),wq(2,2,mrad),xx1(4),xx2(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL inconi,inconz
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
      cc = c_light(2.0)
      csq = cc*cc
C
      kap1 = -l - 1
      kap2 = +l
c
      nsol = 2
      IF (abs(xmj).GE.l) nsol = 1
C
      IF (iway.EQ.2) GO TO 60
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C             OUTWARD INTEGRATION
C
      DO i = 1,2
         DO j = 1,2
            pow(j,i) = 0.0
            qow(j,i) = 0.0
         END DO
      END DO
c     potential redefinition
      DO ir = 1,nmatch
         v(ir) = vv(ir)*rc(ir)*rc(ir)
         b(ir) = bb(ir)*rc(ir)*rc(ir)
         ra(ir) = rc(ir)
      END DO
c initial condition
      CALL inconz(e,l,xmj,kap1,kap2,vv(1),bb(1),rc(1),xx1,xx2)
c dirac equation solution: two cases 1)nsol=2; 2)nsol=1.
      nstart = 1
      IF (nsol.EQ.2) THEN
         CALL kernel2(mrad,nsol,xmj,kap1,kap2,xx1,xx2,e,v,b,ra,dx,
     +                nmatch,nstart,dp,dq,wp,wq)
      ELSE
         CALL kernel1(mrad,xmj,kap1,xx1,e,v,b,ra,dx,nmatch,nstart,dp,
     +                dq,wp,wq)
      END IF
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
      DO 30 n = 1,nmatch
         DO 20 j = 1,nsol
            DO 10 i = 1,nsol
               gc(i,j,n) = wp(i,j,n)/rc(n)
               fc(i,j,n) = wq(i,j,n)/ (rc(n)*cc)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
c
      DO 50 j = 1,nsol
         DO 40 i = 1,nsol
            pow(i,j) = wp(i,j,nmatch)
            qow(i,j) = wq(i,j,nmatch)
   40    CONTINUE
   50 CONTINUE
      RETURN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C             INWARD INTEGRATION
C
   60 CONTINUE
      DO i = 1,2
         DO j = 1,2
            piw(j,i) = 0.0
            qiw(j,i) = 0.0
         END DO
      END DO
c initial condition
      CALL inconi(l,xmj,e,csq,rc(nzero),xx1,xx2)
c     redefinition of order & potential redefinition
      DO ir = nzero,nmatch,-1
         irv = nzero - ir + 1
         v(irv) = vv(ir)*rc(ir)*rc(ir)
         b(irv) = bb(ir)*rc(ir)*rc(ir)
         ra(irv) = rc(ir)
      END DO
      dx1 = -dx
c dirac equation solution: two cases 1)nsol=2; 2)nsol=1.
      nstart = 1
      IF (nsol.EQ.2) THEN
         CALL kernel2(mrad,nsol,xmj,kap1,kap2,xx1,xx2,e,v,b,ra,dx1,
     +                nzero-nmatch+1,nstart,dp,dq,wp,wq)
      ELSE
         CALL kernel1(mrad,xmj,kap1,xx1,e,v,b,ra,dx1,nzero-nmatch+1,
     +                nstart,dp,dq,wp,wq)
      END IF
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
      DO 90 nn = 1,nzero - nmatch + 1
         n = nzero - nn + 1
         DO 80 j = 1,nsol
            DO 70 i = 1,nsol
               gc(i,j,n) = wp(i,j,nn)/ra(nn)
               fc(i,j,n) = wq(i,j,nn)/ (ra(nn)*cc)
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE

c exponential tail

      CALL crtail(
     >            mrad,e,rc,nsol,nzero,csq,
     X            gc,fc)

      DO  j = 1,nsol
         DO  i = 1,nsol
            piw(i,j) = wp(i,j,nzero-nmatch+1)
            qiw(i,j) = wq(i,j,nzero-nmatch+1)
         ENDDO
      ENDDO
c
      END SUBROUTINE coredir
      END MODULE m_coredir
