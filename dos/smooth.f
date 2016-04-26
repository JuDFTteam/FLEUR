      MODULE m_smooth
!
!     the function f(x), defined on a linear mesh,
!     is smoothened by a gaussian
!
      CONTAINS
      SUBROUTINE smooth(e,
     X                  f,
     >                  sigma,n)

      USE m_constants, ONLY: tpi_const
      IMPLICIT NONE

!    Arguments
      INTEGER, INTENT(IN)    :: n
      REAL,    INTENT(INOUT) :: f(n)
      REAL,    INTENT(IN)    :: sigma , e(n)

!    Locals
      REAL :: c , c2 , dx  , f0(n)
      INTEGER :: i , j , j1 , j2 , m1, m

      REAL, ALLOCATABLE :: ee(:)
 
      dx = e(2) - e(1)
      c = dx/(sigma*sqrt(tpi_const))
      c2 = -0.5*(dx/sigma)**2

      m = NINT(sqrt(log(1.0e-8/c)/c2))+1
      ALLOCATE ( ee(m) )
      DO i = 1, m
         ee(i) = c * exp(c2*(i-1)**2)
         IF ( ee(i).LT.1.E-8 ) EXIT
      ENDDO
      m1=i-1
      f0 = f
      f = 0.
     
      DO i = 1 , N
         j1 = i - m1 + 1
         IF ( j1.LT.1 ) j1 = 1
         j2 = i + m1 - 1
         IF ( j2.GT.N ) j2 = N
         DO j = j1 , j2
            f(i) = f(i) + ee(IABS(j-i)+1)*f0(j)
         ENDDO
      ENDDO
      DEALLOCATE ( ee )

      END SUBROUTINE smooth
      END MODULE m_smooth
