      MODULE m_gaussp
!**************************************************************
!     generates gaussian points to exactly integrate spherical
!     harmonics up to lmax, i.e., (lm|l'm') for l,l'<=lmax
!     number of points = (2*lmax+1)*(lmax+1 + mod(lmax+1,2))
!**************************************************************
      CONTAINS 
      SUBROUTINE gaussp(
     >                  lmax,
     <                  vgauss,wt)

      USE m_grule
      USE m_constants
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: lmax
      REAL,    INTENT (OUT) :: vgauss(3,*),wt(*) ! points/weights

      INTEGER :: ngpt,nphi,i,j,k
      REAL    :: delphi,phi,rxy
      REAL    :: xx(lmax/2+1),w(lmax/2+1)

!   determine the number of points cos(theta); ngpt always even
      ngpt= lmax+1 + mod(lmax+1,2)
      CALL grule(            ! outputs ngpt/2 points
     >           ngpt,
     <           xx,w)

!  in phi, use nyquist frequency, i.e.,  2*(lmax+1)
      nphi = 2*lmax+1
      delphi = 2.0*pi_const/nphi

      j = 0
      DO i = 1, ngpt/2
         rxy=sqrt(1.0-xx(i)*xx(i))
         DO k=1,nphi
            phi=k*delphi
            j=j+1
            vgauss(1,j) = rxy*cos(phi)
            vgauss(2,j) = rxy*sin(phi)
            vgauss(3,j) = xx(i)
            wt(j) = w(i)*delphi
            j=j+1
            vgauss(1,j) = vgauss(1,j-1)
            vgauss(2,j) = vgauss(2,j-1)
            vgauss(3,j) = -xx(i)
            wt(j) = w(i)*delphi
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE gaussp
      END MODULE m_gaussp
