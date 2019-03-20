MODULE m_kkintgr

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_kkintgr
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  Performs the Kramer-Kronig-Transformation to obtain the Green's function 
   !>  in the complex plane from the imaginary part calculated on the real axis
   !
   !------------------------------------------------------------------------------
   CONTAINS

   SUBROUTINE kkintgr_real(nz,e,ne,sigma,del,bot,im,g)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated

      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      USE m_types
      USE m_constants
      USE m_intgr

      IMPLICIT NONE

      COMPLEX,             INTENT(OUT)    :: g(:)  !green's function
      COMPLEX,             INTENT(IN)     :: e(:)  !energy contour, on which to calculate G(z)
      REAL,                INTENT(IN)     :: im(:) !previously calculated imaginary part on equidistant mesh at E+i*sigma

      INTEGER,             INTENT(IN)     :: nz
      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: sigma       
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: bot

      !Local declarations:
      INTEGER iz, j, i_sing
      COMPLEX  ez
      REAL, ALLOCATABLE :: integrand(:)
      REAL result

      ALLOCATE (integrand(ne))


      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(ne,sigma,del,nz,bot) &
      !$OMP SHARED(g,im,e) &
      !$OMP PRIVATE(j,ez,integrand,result,i_sing)

      !$OMP DO
      DO iz = 1, nz

         !We calculate the integral by removing the singularity and treating it analytically
         ez = REAL(e(iz)) - bot

         i_sing = INT(ez/del) + 1

         DO j = 1, i_sing -1

            !Formula for the integrand if the contour is on the same imaginary part as before
            integrand(j) = (im(j)-im(i_sing))/((j-1)*del-ez)

         ENDDO

         !at j = iz we estimate the integrand with the centered three-point formula for the derivative of im(e) 
         !(at the edges we use backward/forward formulas)

         IF(i_sing.EQ.1) THEN 
            integrand(i_sing) = (im(2)-im(1))/(del) 
         ELSE IF(i_sing.EQ.ne) THEN
            integrand(i_sing) = (im(ne)-im(ne-1))/del
         ELSE
            integrand(i_sing) = (im(i_sing+1)-im(i_sing-1))/(2*del) 
         END IF

         DO j = i_sing+1, ne

            !Formula for the integrand if the contour is on the same imaginary part as before
            integrand(j) = (im(j)-im(i_sing))/((j-1)*del-ez)

         ENDDO

         CALL intgz0(integrand(:), del,ne,result,.false.)

         !The singularity is treated analytically

         g(iz) = -1/pi_const * result  + ImagUnit * im(i_sing)
         !(result - im(i_sing) * LOG(ABS(ez/REAL(ez-(ne-1)*del)))) + ImagUnit * im(i_sing)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL


   ENDSUBROUTINE kkintgr_real

   SUBROUTINE kkintgr_complex(nz,e,ne,sigma,del,bot,im,g)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated

      !G(z) = -1/pi * int_bot^top dE' 1/(z-E') * Im(G(E'+i*delta))

      !The imaginary part is defined on a equidistant energy mesh with spacing del

      USE m_types
      USE m_constants
      USE m_intgr

      IMPLICIT NONE

      COMPLEX,             INTENT(OUT)    :: g(:)  !green's function
      COMPLEX,             INTENT(IN)     :: e(:)  !energy contour, on which to calculate G(z)
      REAL,                INTENT(IN)     :: im(:) !previously calculated imaginary part on equidistant mesh at E+i*sigma

      INTEGER,             INTENT(IN)     :: nz
      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: sigma       
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: bot

      INTEGER iz, j
      COMPLEX ez
      COMPLEX,ALLOCATABLE :: integrand(:)
      REAL re, imag

      ALLOCATE (integrand(ne))

      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(ne,del,nz,bot) &
      !$OMP SHARED(g,e,im) &
      !$OMP PRIVATE(j,ez,integrand,re,imag)

      !$OMP DO

      DO iz = 1, nz

         ez = e(iz) - bot 

         DO j = 1, ne

            integrand(j) = 1/(ez-(j-1)*del) * im(j)

         ENDDO
         !We don't use the normal integration routines given in m_intgr because
         !the functions are very peaky and these are handled better by a simple trapezian method
         CALL trapz(REAL(integrand(:)),del,ne,re)
         CALL trapz(AIMAG(integrand(:)),del,ne,imag)


         g(iz) = -1/pi_const*CMPLX(re,imag)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL


   END SUBROUTINE kkintgr_complex

   !trapezoidal rule (more efficient than simpsons scheme in e.g intgz0
   ! for our peak-like functions)
   SUBROUTINE trapz(y,h,n,z)

      IMPLICIT NONE

      REAL,          INTENT(IN)     :: y(n)
      REAL,          INTENT(OUT)    :: z

      INTEGER,       INTENT(IN)     :: n
      REAL,          INTENT(IN)     :: h


      INTEGER i

      z = y(1)
      DO i = 2, n-1

         z = z + 2*y(i)

      ENDDO
      z = z + y(n)

      z = z*h/2.0

      END SUBROUTINE trapz

ENDMODULE m_kkintgr 