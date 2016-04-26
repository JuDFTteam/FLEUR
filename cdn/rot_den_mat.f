      MODULE m_rotdenmat
      use m_juDFT
      CONTAINS
      SUBROUTINE rot_den_mat(
     >                       alph,beta,
     X                       rho11,rho22,rho21)
c***********************************************************************
c This subroutine rotates the direction of the magnetization of the 
c density matrix by multiplying with the unitary 2x2 spin rotation
c matrix. --> U*rho*U^dagger
c Philipp Kurz 2000-02-03
c***********************************************************************

      IMPLICIT NONE

C     .. Scalar Arguments ..
      REAL, INTENT    (IN) :: alph,beta
      REAL, INTENT    (INOUT) :: rho11
      REAL, INTENT    (INOUT) :: rho22
      COMPLEX, INTENT (INOUT) :: rho21
C     ..
C     .. Local Scalars ..
      INTEGER ispin
      REAL eps
      COMPLEX ci
C     ..
C     .. Local Arrays ..
      COMPLEX u2(2,2),rho(2,2),rhoh(2,2)
C     ..
      
      eps = 1.0e-10
      ci = cmplx(0.0,1.0)

c---> set up the unitary 2x2 spin rotation matrix U^(2)
      u2(1,1) =  exp(-ci*alph/2)*cos(beta/2)
      u2(1,2) = -exp(-ci*alph/2)*sin(beta/2)
      u2(2,1) =  exp( ci*alph/2)*sin(beta/2)
      u2(2,2) =  exp( ci*alph/2)*cos(beta/2)
      
      rho(1,1) = cmplx(rho11,0.0)
      rho(2,2) = cmplx(rho22,0.0)
      rho(2,1) =       rho21
      rho(1,2) = conjg(rho21)

c---> first calculate U*rho
      rhoh(1,1) = u2(1,1)*rho(1,1) + u2(1,2)*rho(2,1)
      rhoh(1,2) = u2(1,1)*rho(1,2) + u2(1,2)*rho(2,2)
      rhoh(2,1) = u2(2,1)*rho(1,1) + u2(2,2)*rho(2,1)
      rhoh(2,2) = u2(2,1)*rho(1,2) + u2(2,2)*rho(2,2)
c---> now calculate (U*rho)*U^dagger
      rho(1,1) = rhoh(1,1)*conjg(u2(1,1))
     +         + rhoh(1,2)*conjg(u2(1,2))
      rho(1,2) = rhoh(1,1)*conjg(u2(2,1))
     +         + rhoh(1,2)*conjg(u2(2,2))
      rho(2,1) = rhoh(2,1)*conjg(u2(1,1))
     +         + rhoh(2,2)*conjg(u2(1,2))
      rho(2,2) = rhoh(2,1)*conjg(u2(2,1))
     +         + rhoh(2,2)*conjg(u2(2,2))

c---> check wether the diagonal elements of the rotated density
c---> are real.
      DO ispin = 1,2
         IF (aimag(rho(ispin,ispin)).GT.eps) THEN
            WRITE(16,8000)
            WRITE( 6,8000)
            CALL juDFT_error("rotation of mag. failed",calledby
     +           ="rot_den_mat")
         ENDIF
      ENDDO
 8000 FORMAT('After the rotation of the density matrix in the'/
     +       'muffin-tin sphere one diagonal element of the'/
     +       '(hermitian) density matrix is not real. That means'/
     +       'that the density matrix was probably damaged.')

      rho11 = real(rho(1,1))
      rho22 = real(rho(2,2))
      rho21 =      rho(2,1)

      END SUBROUTINE rot_den_mat
      END MODULE m_rotdenmat

