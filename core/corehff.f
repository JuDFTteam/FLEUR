c............................................................corehff
      SUBROUTINE corehff(mrad,kap1,kap2,xmj,s,nsol,bhf,gck,fck,
     +                   rc,dx,jtop)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE THE RELATIVISTIC HYPERFINE FIELDS FOR THE            *
C   *                  CURRENT  CORE STATE S                           *
C   *                                                                  *
C   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
C   *                                                                  *
C   ********************************************************************

      IMPLICIT NONE
C CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
C                                      ELECTRON CHARGE     IN ESU
C                                      BOHR-RADIUS         IN CM
c
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
C        THE FACTOR  I  HAS BEEN OMITTED
C
C     .. Parameters ..
      REAL e0
      PARAMETER (e0=1.6021892e-19*2.997930e+09)
      REAL a0
      PARAMETER (a0=0.52917706e-08)
      REAL cautog
      PARAMETER (cautog=e0/ (a0*a0))
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad
      REAL bhf,dx,xmj
      INTEGER jtop,kap1,kap2,nsol,s
C     ..
C     .. Array Arguments ..
      REAL fck(2,2,mrad),gck(2,2,mrad),rc(mrad)
C     ..
C     .. Local Scalars ..
      INTEGER n
C     ..
C     .. Local Arrays ..
      REAL ame(2,2),rint(mrad)
C     ..
C     .. External Functions ..
      REAL rsimp
      EXTERNAL rsimp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real,sqrt
C     ..
      ame(1,1) = 4.0*kap1*xmj/ (4.0*kap1*kap1-1.0)
      IF (nsol.EQ.2) THEN
         ame(2,2) = 4.0*kap2*xmj/ (4.0*kap2*kap2-1.0)
         ame(2,1) = sqrt(0.25- (xmj/real(kap1-kap2))**2)
         ame(1,2) = ame(2,1)
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (nsol.EQ.1) THEN
         DO n = 1,jtop
            rint(n) = (gck(1,s,n)*fck(1,s,n)+fck(1,s,n)*gck(1,s,n))*
     +                ame(1,1)
         END DO
      ELSE
         DO n = 1,jtop
            rint(n) = (gck(1,s,n)*fck(1,s,n)+fck(1,s,n)*gck(1,s,n))*
     +                ame(1,1) + (gck(2,s,n)*fck(2,s,n)+
     +                fck(2,s,n)*gck(2,s,n))*ame(2,2) +
     +                (gck(2,s,n)*fck(1,s,n)+fck(2,s,n)*gck(1,s,n))*
     +                ame(2,1) + (gck(1,s,n)*fck(2,s,n)+
     +                fck(1,s,n)*gck(2,s,n))*ame(1,2)
         END DO
      END IF
      bhf = -cautog*rsimp(mrad,rint,rc,jtop,dx)*0.001
c      write(6,'(''hf='',e14.7)') BHF
      RETURN
      END
