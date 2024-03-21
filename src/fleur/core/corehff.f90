MODULE m_corehff

   CONTAINS

   SUBROUTINE corehff(mrad,kap1,kap2,xmj,s,nsol,bhf,gck,fck,rc,dx,jtop)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE THE RELATIVISTIC HYPERFINE FIELDS FOR THE            *
!   *                  CURRENT  CORE STATE S                           *
!   *                                                                  *
!   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
!   *                                                                  *
!   ********************************************************************

      USE m_constants
      USE m_rsimp

      IMPLICIT NONE
! CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!                                      ELECTRON CHARGE     IN ESU
!                                      BOHR-RADIUS         IN CM
!
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
!        THE FACTOR  I  HAS BEEN OMITTED

      INTEGER, INTENT (IN) :: mrad
      INTEGER, INTENT (IN) :: jtop,kap1,kap2,nsol,s
      REAL, INTENT (OUT)   :: bhf
      REAL, INTENT (IN)    :: dx,xmj
      REAL, INTENT (IN)    :: fck(2,2,mrad),gck(2,2,mrad),rc(mrad)

      INTEGER n

      REAL e0, a0, cautog
      REAL ame(2,2),rint(mrad)
      
      a0 = bohr_to_angstrom_const * 1.0e-8
      e0 = 1.6021892e-19 * 2.997930e+09
      cautog = e0 / (a0*a0)

      ame(1,1) = 4.0*kap1*xmj/ (4.0*kap1*kap1-1.0)
      IF (nsol.EQ.2) THEN
         ame(2,2) = 4.0*kap2*xmj/ (4.0*kap2*kap2-1.0)
         ame(2,1) = sqrt(0.25- (xmj/real(kap1-kap2))**2)
         ame(1,2) = ame(2,1)
      END IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (nsol.EQ.1) THEN
         DO n = 1,jtop
            rint(n) = (gck(1,s,n)*fck(1,s,n)+fck(1,s,n)*gck(1,s,n)) * ame(1,1)
         END DO
      ELSE
         DO n = 1,jtop
            rint(n) = (gck(1,s,n)*fck(1,s,n)+fck(1,s,n)*gck(1,s,n)) * ame(1,1) + &
                      (gck(2,s,n)*fck(2,s,n)+fck(2,s,n)*gck(2,s,n)) * ame(2,2) + &
                      (gck(2,s,n)*fck(1,s,n)+fck(2,s,n)*gck(1,s,n)) * ame(2,1) + &
                      (gck(1,s,n)*fck(2,s,n)+fck(1,s,n)*gck(2,s,n)) * ame(1,2)
         END DO
      END IF
      bhf = -cautog*rsimp(mrad,rint,rc,jtop,dx)*0.001
!      write(oUnit,'(''hf='',e14.7)') BHF

   END SUBROUTINE corehff
      
END MODULE m_corehff
