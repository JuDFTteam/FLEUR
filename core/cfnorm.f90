MODULE m_cfnorm

   CONTAINS

   SUBROUTINE cfnorm(mrad,is,it,nsol,nmatch,jtop,var,gc,fc,rc,rc2,dx,gck,fck)
!...........................................................cfnorm
! wavefunctions normalization
!

      USE m_rsimp

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mrad
      REAL dx
      INTEGER is,it,jtop,nmatch,nsol

      REAL fc(2,2,mrad),fck(2,2,mrad),gc(2,2,mrad),gck(2,2,mrad),rc(mrad),rc2(mrad),var(4)

      REAL xnorm
      INTEGER i,j,k,n

      REAL rint(mrad)


      DO n = 1,mrad
         DO k = 1,nsol
            gck(k,is,n) = 0.0
            fck(k,is,n) = 0.0
         END DO
      END DO
!                         ---------------------------------
!--->                     NORMALIZE WAVEFUNCTIONS ACCORDING
!                               TO MATCHING CONDITIONS
!                         ---------------------------------
!                                    INWARD - SOLUTION
      DO n = nmatch,jtop
         DO j = 1,nsol
            DO i = 1,nsol
               gc(i,j,n) = gc(i,j,n)*var(2)
               fc(i,j,n) = fc(i,j,n)*var(2)
            END DO
         END DO
      END DO

      IF (nsol.EQ.2) THEN
!                                   OUTWARD - SOLUTION
         DO n = 1, (nmatch-1)
            DO i = 1,nsol
               gc(i,it,n) = gc(i,it,n)*var(3)
               fc(i,it,n) = fc(i,it,n)*var(3)
            END DO
         END DO
!                                    INWARD - SOLUTION
         DO n = nmatch,jtop
            DO i = 1,nsol
               gc(i,it,n) = gc(i,it,n)*var(4)
               fc(i,it,n) = fc(i,it,n)*var(4)
            END DO
         END DO
      END IF
!                                    SUM FOR EACH KAPPA
      DO n = 1,jtop
         DO k = 1,nsol
            gck(k,is,n) = 0.00
            fck(k,is,n) = 0.00
         END DO
      END DO
      DO n = 1,jtop
         DO k = 1,nsol
            DO j = 1,nsol
               gck(k,is,n) = gck(k,is,n) + gc(k,j,n)
               fck(k,is,n) = fck(k,is,n) + fc(k,j,n)
            END DO
         END DO
      END DO
!                       -----------------------------------
!                       CALCULATE  NORM  AND NORMALIZE TO 1
!                       -----------------------------------
      DO n = 1,jtop
         rint(n) = 0.00
         DO k = 1,nsol
            rint(n) = rint(n) + rc2(n)* (gck(k,is,n)**2+fck(k,is,n)**2)
         END DO
      END DO
      xnorm = rsimp(mrad,rint,rc,jtop,dx)
      xnorm = 1.00/sqrt(xnorm)
      DO n = 1,jtop
         DO k = 1,nsol
            gck(k,is,n) = gck(k,is,n)*xnorm
            fck(k,is,n) = fck(k,is,n)*xnorm
         END DO
      END DO

   END SUBROUTINE cfnorm      
END MODULE m_cfnorm
