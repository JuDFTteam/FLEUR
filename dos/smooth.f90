MODULE m_smooth
!
!     the function f(x), defined on a linear mesh,
!     is smoothened by a gaussian
!

   INTERFACE smooth
      PROCEDURE smooth_r, smooth_c
   END INTERFACE

   CONTAINS

   SUBROUTINE smooth_r(e,f,sigma,n)

      INTEGER, INTENT(IN)    :: n
      REAL,    INTENT(INOUT) :: f(:)
      REAL,    INTENT(IN)    :: sigma , e(:)

      COMPLEX, ALLOCATABLE :: f_c(:)

      ALLOCATE(f_c(n))
      f_c = f
      CALL smooth_c(e,f_c,sigma,n)
      f = REAL(f_c)

   END SUBROUTINE smooth_r

   SUBROUTINE smooth_c(e,f,sigma,n)

   USE m_constants, ONLY: tpi_const
   IMPLICIT NONE

!    Arguments
   INTEGER, INTENT(IN)    :: n
   COMPLEX, INTENT(INOUT) :: f(:)
   REAL,    INTENT(IN)    :: sigma , e(:)

!    Locals
   REAL :: c , c2 , dx
   INTEGER :: i , j , j1 , j2 , m1, m
   COMPLEX :: f0(n)

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

   END SUBROUTINE smooth_c
END MODULE m_smooth
