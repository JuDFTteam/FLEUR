      MODULE m_rcerf
      use m_juDFT
c*********************************************************************
c     calculates  real( erf(x+iy) ) for z=x+iy in the first quadrant.
c             m. weinert   may 1987
c     new declaration part
c                       s. bl"ugel, IFF, Nov.97

      PRIVATE
      PUBLIC rcerf
      CONTAINS
c*********************************************************************
      REAL FUNCTION rcerf(x,y)
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      REAL    x,y
C     ..
C     .. Local Scalars ..
      COMPLEX z
      REAL wr,wi
c
      z = cmplx(x,y)
c
c--->    calculate w(z')=exp(-z'**2)*erfc(-iz') for -iz'=z
c
      CALL wofz(y,x,wr,wi)
c
c--->    erf(x+iy)=1-exp(-z**2)*conjg( w(y+ix) )
c
      RCERF = 1.0 - real( exp( -z*z ) * cmplx( wr,-wi ) )

      END FUNCTION rcerf

c*********************************************************************
      SUBROUTINE wofz(x,y,wr,wi)
c     calculates  w(z) = exp(-z**2) erfc(-iz)   for z = x + iy in the
c     first quadrant of the complex plane. based on acm algorithm 363
c     by w. gautschi, comm. acm 12, 635 (1969). the accuracy is about
c     10**-10.
c                       m. weinert    may 1987
c     new declaration part
c                       s. bl"ugel, IFF, Nov.97
c*********************************************************************
      USE m_constants
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      REAL    x,y,wr,wi
C     ..
C     .. Local Scalars ..
      INTEGER icap,nu,n
      REAL    c,h,h2,plam,r1,r2,s,s1,s2,tsqpi,t1,t2
      LOGICAL bol

c--->    2/sqrt(pi) required for the normalization
      tsqpi = 2.0/sqrt(pi_const) 
c
c--->    stop if not in first quadrant 
      IF ( x<0.0 .OR. y<0 )  CALL juDFT_error("wofz",calledby ="rcerf")
c
      IF ( y.LT.4.29 .AND. x.LT.5.33 ) THEN
        s=(1.-y/4.29)*sqrt(1.-x*x/28.41)
        h=1.6*s
        h2=2*h
        icap=6+23*s
        nu=9+21*s
        plam=h2**icap
        bol=.true.
      ELSE
        h=0.0
        icap=0
        nu=8
        bol=.false.
      END IF
c
      r1=0.0
      r2=0.0
      s1=0.0
      s2=0.0
      DO n = nu,0,-1
         t1=y+h+(n+1)*r1
         t2=x-(n+1)*r2
         c=0.5/(t1*t1+t2*t2)
         r1=c*t1
         r2=c*t2
         IF ( bol .AND. n.LE.icap ) then
           t1=plam+s1
           s1=r1*t1-r2*s2
           s2=r2*t1+r1*s2
           plam=plam/h2
         END IF
      ENDDO
      IF ( bol ) THEN
        wr=tsqpi*s1
        wi=tsqpi*s2
      ELSE
        wr=tsqpi*r1
        wi=tsqpi*r2
      END IF
      IF ( y.EQ.0.0 ) wr=exp(-x*x)

      END SUBROUTINE wofz

      END MODULE m_rcerf
