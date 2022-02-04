MODULE m_lorentzian_smooth

   USE m_constants

   IMPLICIT NONE

   !PARAMETER FOR LORENTZIAN SMOOTHING
   REAL,    PARAMETER :: cut              = 1e-8

   INTERFACE lorentzian_smooth
      PROCEDURE lorentzian_smooth_r, lorentzian_smooth_c
   END INTERFACE

   CONTAINS

   SUBROUTINE lorentzian_smooth_r(e,f,sigma,n)

      INTEGER, INTENT(IN)    :: n
      REAL,    INTENT(INOUT) :: f(:)
      REAL,    INTENT(IN)    :: sigma , e(:)

      COMPLEX, ALLOCATABLE :: f_c(:)

      f_c = f
      CALL lorentzian_smooth_c(e,f_c,sigma,n)
      f = REAL(f_c)

   END SUBROUTINE lorentzian_smooth_r

   !This is essentially smooth out of m_smooth but with a lorentzian distribution
   SUBROUTINE lorentzian_smooth_c(e,f,sigma,n)

      INTEGER, INTENT(IN)    :: n
      COMPLEX, INTENT(INOUT) :: f(:)
      REAL,    INTENT(IN)    :: sigma
      REAL,    INTENT(IN)    :: e(:)

      REAL :: dx
      COMPLEX :: f0(n)
      REAL :: ee(n)
      INTEGER :: ie , je , j1 , j2 , numPoints

      f0 = f
      f = 0.0
      ee = 0.0
      dx = e(2)-e(1)

      DO ie =1,  n

         ee(ie) = 1/pi_const * sigma/dx * 1.0/((ie-1)**2+(sigma/dx)**2)

         IF ( ee(ie).LT.cut ) EXIT
      ENDDO
      numPoints = ie - 1

      DO ie = 1, n

         j1 = ie - numPoints + 1
         j1 = MERGE(1,j1,j1.LT.1)
         j2 = ie + numPoints - 1
         j2 = MERGE(n,j2,j2.GT.n)
         DO je = j1 , j2
            f(ie) = f(ie) + ee(IABS(je-ie)+1)*f0(je)
         ENDDO

      ENDDO

   END SUBROUTINE lorentzian_smooth_c

END MODULE m_lorentzian_smooth