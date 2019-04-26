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

   USE m_constants
   USE m_juDFT

   INTEGER, PARAMETER :: method_maclaurin = 1
   INTEGER, PARAMETER :: method_deriv     = 2

   INTEGER, PARAMETER :: rectangular      = 1
   INTEGER, PARAMETER :: semicircle       = 2

   CONTAINS 

   SUBROUTINE kkintgr(im,eb,del,ne,g,ez,l_conjg,shape,nz,method)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated
      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      !The dominant source of error for this routine is a insufficiently dense energy mesh on the real axis
      !TODO: Some way to estimate the error (maybe search for the sharpest peak and estimate from width)


      IMPLICIT NONE

      !Information about the integrand
      REAL,          INTENT(IN)  :: im(ne)      !part of the green's function on the real axis 
      REAL,          INTENT(IN)  :: eb          !Bottom energy cutoff
      REAL,          INTENT(IN)  :: del         !Energy step on the real axis 
      INTEGER,       INTENT(IN)  :: ne          !Number of energy points on the real axis

      !Information about the complex energy contour
      COMPLEX,       INTENT(OUT) :: g(nz)       !Green's function on the complex plane
      COMPLEX,       INTENT(IN)  :: ez(nz)      !Complex energy contour
      LOGICAL,       INTENT(IN)  :: l_conjg     !Switch determines wether we calculate g on the complex conjugate of the contour ez
      INTEGER,       INTENT(IN)  :: shape       !Determines wether we have a rectangular (1) or a semicircle contour(2)
      INTEGER,       INTENT(IN)  :: nz          !Number of energy points on the complex contour
  
      !Information about the method 
      INTEGER,       INTENT(IN)  :: method      !Integer associated with the method to be used (definitions above)

      REAL ::  im_calc(ne)  !Array where the smoothed version of im is stored

      INTEGER  iz,n1,n2
      REAL     sigma,re_n1,re_n2


      IF(shape.EQ.rectangular) THEN
         !If we have a rectangular contour we only need  to smooth the imaginary part once
         im_calc = im
         sigma = AIMAG(ez(1)) 
         IF(sigma.NE.0.0) THEN
            CALL lorentzian_smooth(del,im_calc,sigma,ne)
         ENDIF
         IF(l_conjg) im_calc = -im_calc
      ENDIF
      
      g = 0.0
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(nz,ne,method,shape,del,eb,l_conjg) &
      !$OMP SHARED(g,ez,im) &
      !$OMP PRIVATE(iz,n1,n2,sigma,re_n1,re_n2) &
      !$OMP PRIVATE(im_calc)

      !$OMP DO
      DO iz = 1, nz

         SELECT CASE(shape)

         CASE (rectangular)
            IF(nz.EQ.ne) THEN
               g(iz) = re_ire(im_calc,ne,iz,method) + ImagUnit * im_calc(iz)
            ELSE IF(nz.LT.ne) THEN
               !Next point to the left
               n1 = INT((REAL(ez(iz))-eb)/del) +1
               !next point to the right
               n2 = n1 + 1

               !Here we perform the integrations
               re_n2 = re_ire(im_calc,ne,n2,method)
               re_n1 = re_ire(im_calc,ne,n1,method)

               !Real Part (only split up to make it readable)
               g(iz) = (re_n2-re_n1)/del * (REAL(ez(iz))-(n1-1)*del-eb) + re_n1
               !Imaginary Part
               g(iz) = g(iz) + ImagUnit *( (im_calc(n2)-im_calc(n1))/del * (REAL(ez(iz))-(n1-1)*del-eb) + im_calc(n1) )
            ELSE
               CALL juDFT_error("Complex Grid finer than on the real axis",calledby="kkintgr")
            ENDIF
         CASE (semicircle)
            !For a semicircle we need to smooth to the correct imaginary part first
            im_calc = im
            sigma = AIMAG(ez(iz)) 
            IF(sigma.NE.0.0) THEN
               CALL lorentzian_smooth(del,im_calc,sigma,ne)
            ENDIF
            IF (l_conjg) im_calc = -im_calc
            !Next point to the left
            n1 = INT((REAL(ez(iz))-eb)/del) +1
            !next point to the right
            n2 = n1 + 1

            !Here we perform the integrations
            re_n2 = re_ire(im_calc,ne,n2,method)
            re_n1 = re_ire(im_calc,ne,n1,method)

            !Real Part (only split up to make it readable)
            g(iz) = (re_n2-re_n1)/del * (REAL(ez(iz))-(n1-1)*del-eb) + re_n1
            !Imaginary Part
            g(iz) = g(iz) + ImagUnit *( (im_calc(n2)-im_calc(n1))/del * (REAL(ez(iz))-(n1-1)*del-eb) + im_calc(n1) )
         CASE default
            CALL juDFT_error("Not a valid shape for the energy contour",calledby="kkintgr")
         END SELECT

      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

   END SUBROUTINE kkintgr

   REAL FUNCTION re_ire(im,ne,ire,method)

      IMPLICIT NONE 

      REAL,    INTENT(IN)  :: im(ne) !Imaginary part
      INTEGER, INTENT(IN)  :: ne     !Dimension of the energy grid
      INTEGER, INTENT(IN)  :: ire    !Position where to calculate the real part
      INTEGER, INTENT(IN)  :: method !Method to be used
      INTEGER i,j
      REAL    y

      re_ire = 0.0
      SELECT CASE(method)
               
      CASE (method_maclaurin)
         !Calculate the real part on the same energy points as the imaginary part
         !regardless of the contour
         !If i is odd skip the odd points and the other way around and use the trapezian method
         DO j = 1, INT(ne/2.0)
            IF(MOD(ire,2).EQ.0) THEN
               i = 2*j-1
            ELSE
               i = 2*j 
            ENDIF
            y = - 1/pi_const * 2.0 * im(i)/REAL(ire-i) 
            IF(j.EQ.1.OR.j.EQ.INT(ne/2.0)) y = y/2.0
            re_ire = re_ire + y
         ENDDO

      CASE (method_deriv)
         !Remove the singularity and treat it analytically
         DO j = 1, ne
            IF(j-ire.NE.0) THEN
               y = -1/pi_const * (im(j)-im(ire))/REAL(ire-j)
            ELSE
               IF(ire.EQ.1) THEN
                  y = -1/pi_const * (im(2)-im(1))
               ELSE IF(ire.EQ.ne) THEN
                  y = -1/pi_const * (im(ne)-im(ne-1))
               ELSE
                  y = -1/pi_const * (im(ire+1)-im(ire-1))/2.0
               ENDIF
            ENDIF
            IF(j.EQ.1.OR.j.EQ.INT(ne/2.0)) y = y/2.0
            re_ire = re_ire + y
         ENDDO
      CASE default
         CALL juDFT_error("No valid method for KK-integration chosen",calledby="kkintgr")
      END SELECT
   END FUNCTION re_ire

   !This is essentially smooth out of m_smooth but with a lorentzian distribution
   SUBROUTINE lorentzian_smooth(dx,f,sigma,n)

      !USE m_constants, ONLY: pi_const
      IMPLICIT NONE

!    Arguments
      INTEGER, INTENT(IN)    :: n
      REAL,    INTENT(INOUT) :: f(n)
      REAL,    INTENT(IN)    :: sigma 
      REAL,    INTENT(IN)    :: dx

!    Locals
      REAL :: c , f0(n)
      INTEGER :: i , j , j1 , j2 , m1, m

      REAL, ALLOCATABLE :: ee(:)
 
      c = dx/(pi_const) * sigma

      m = NINT(sqrt(c/1e-14-sigma**2))+1
      ALLOCATE ( ee(m) )
      DO i = 1, m
         ee(i) = c * 1.0/(sigma**2+(i-1)**2*dx**2)
         IF ( ee(i).LT.1.E-14 ) EXIT
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

      END SUBROUTINE lorentzian_smooth

      !General Purpose trapezian method integration (not used in kkintgr)
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

END MODULE m_kkintgr