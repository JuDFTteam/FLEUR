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

   !PARAMETER FOR LORENTZIAN SMOOTHING
   REAL,    PARAMETER :: cut              = 1e-8

   CONTAINS 

   SUBROUTINE kkintgr(im,eb,del,ne,g,ez,l_conjg,shape,nz,method)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated
      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      !The dominant source of error for this routine is a insufficiently dense energy mesh on the real axis
      !TODO: Some way to estimate the error (maybe search for the sharpest peak and estimate from width)

      USE m_smooth

      IMPLICIT NONE

      !Information about the integrand
      REAL,          INTENT(IN)  :: im(ne)      !Imaginary part of the green's function on the real axis 
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

      REAL ::  im_calc(ne),e(ne)  !Array where the smoothed version of im is stored

      INTEGER  iz,n1,n2,i,count
      REAL     sigma,re_n1,re_n2,im_n1,im_n2

      DO i = 1, ne
         e(i) = (i-1) * del + eb
      ENDDO
      
      CALL timestart("kkintgr: integration")
      g = 0.0
      sigma = 0.0
      !!$OMP PARALLEL DEFAULT(none) &
      !!$OMP SHARED(nz,ne,method,shape,del,eb,l_conjg) &
      !!$OMP SHARED(g,ez,im,e) &
      !!$OMP PRIVATE(iz,n1,n2,sigma,re_n1,re_n2,im_n1,im_n2,im_calc) 
   !
      !!$OMP DO
      DO iz = 1, nz
         IF(method.EQ.3) THEN
            g(iz) = g_circle(im,ne,MERGE(conjg(ez(iz)),ez(iz),l_conjg),del,eb)
         ELSE
            IF(AIMAG(ez(iz)).NE.0.0.AND.AIMAG(ez(iz)).NE.sigma) THEN
               !Sigma is changed, so we need to smooth here
               im_calc = im !Get the original version
               sigma = AIMAG(ez(iz)) 
               CALL smooth(e,im_calc,sigma,ne)   
            ENDIF
            !Next point to the left
            n1 = INT((REAL(ez(iz))-eb)/del) +1
            !next point to the right
            n2 = n1 + 1
            !Here we perform the Kramers-kronig-Integration
            re_n2 = re_ire(im_calc,ne,n2,method)
            re_n1 = re_ire(im_calc,ne,n1,method)
            !Interpolate to the energy ez(iz)
            !Real Part 
            g(iz) = (re_n2-re_n1)/del * (REAL(ez(iz))-(n1-1)*del-eb) + re_n1
            !Imaginary Part (0 outside of the energy range)
            im_n1 = MERGE(im_calc(n1),0.0,(n1.LE.ne).AND.(n1.GE.1)) 
            im_n2 = MERGE(im_calc(n2),0.0,(n2.LE.ne).AND.(n2.GE.1)) 
            g(iz) = g(iz) + ImagUnit *( (im_n2-im_n1)/del * (REAL(ez(iz))-(n1-1)*del-eb) + im_n1 )
            IF(ISNAN(AIMAG(g(iz))).OR.ISNAN(REAL(g(iz)))) THEN
               CALL juDFT_error("Kkintgr failed",calledby="kkintgr")
            ENDIF
            IF(l_conjg) g(iz) = conjg(g(iz))
         ENDIF
      ENDDO
      !!$OMP END DO
      !!$OMP END PARALLEL
      CALL timestop("kkintgr: integration")

   END SUBROUTINE kkintgr

   COMPLEX FUNCTION g_circle(im,ne,z,del,eb)

      IMPLICIT NONE

      REAL,    INTENT(IN) :: im(ne)
      INTEGER, INTENT(IN) :: ne 
      COMPLEX, INTENT(IN) :: z
      REAL,    INTENT(IN) :: del 
      REAL,    INTENT(IN) :: eb

      COMPLEX :: integrand(ne)
      INTEGER :: i

      integrand = 0.0
      DO i = 1, ne
         integrand(i) = 1.0/(z-(i-1)*del-eb) * im(i)
      ENDDO

      g_circle = -1/pi_const *( trapz(REAL(integrand(:)),del,ne) &
                        + ImagUnit * trapz(AIMAG(integrand(:)),del,ne))
   END FUNCTION g_circle

   REAL FUNCTION re_ire(im,ne,ire,method)

      IMPLICIT NONE 

      REAL,    INTENT(IN)  :: im(ne) !Imaginary part
      INTEGER, INTENT(IN)  :: ne     !Dimension of the energy grid
      INTEGER, INTENT(IN)  :: ire    !Position where to calculate the real part
      INTEGER, INTENT(IN)  :: method !Method to be used
      INTEGER i,j
      REAL    y,im_ire

      re_ire = 0.0
      im_ire = MERGE(im(ire),0.0,(ire.LE.ne).AND.(ire.GE.1)) 
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
               y = -1/pi_const * (im(j)-im_ire)/REAL(ire-j)
            ELSE
               IF(ire.EQ.1) THEN
                  y = -1/pi_const * (im(2)-im(1))
               ELSE IF(ire.EQ.ne) THEN
                  y = -1/pi_const * (im(ne)-im(ne-1))
               ELSE IF((ire.LE.ne).AND.(ire.GE.1)) THEN
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

      IMPLICIT NONE

!    Arguments
      INTEGER, INTENT(IN)    :: n
      REAL,    INTENT(INOUT) :: f(n)
      REAL,    INTENT(IN)    :: sigma 
      REAL,    INTENT(IN)    :: dx

!    Locals
      REAL :: c , f0(n)
      INTEGER :: i , j , j1 , j2 , m1, m

      INTEGER ie,je
      REAL fac,sum,energy
      

      !LORENTZIAN SMOOTHING COMPARE https://heasarc.nasa.gov/xanadu/xspec/models/glsmooth.html
      f0 = f
      f = 0.0
      DO ie = 1, n 

         energy = (ie-0.5)*dx

         fac = 1.0
         je = ie
         sum = 0.0
         DO WHILE(fac.GT.cut.AND.je.GE.1) 
            fac = 1/pi_const * ( atan(2*((je)*dx-energy)/sigma) &
                                 -atan(2*((je-1)*dx-energy)/sigma) )

            f(je) = f(je) + fac * f0(ie)
            sum = sum + fac
            je = je-1
         ENDDO

         fac = 1.0
         je = ie+1
         DO WHILE(fac.GT.cut.AND.je.LE.n) 
            fac = 1/pi_const * ( atan(2*((je)*dx-energy)/sigma) &
                                 -atan(2*((je-1)*dx-energy)/sigma) )

            f(je) = f(je) + fac * f0(ie)
            sum = sum + fac
            je = je+1
         ENDDO
      ENDDO
      IF(ANY(ISNAN(f(:)))) CALL juDFT_error("Smoothing failed", calledby="lorentzian_smooth")


   END SUBROUTINE lorentzian_smooth

   !General Purpose trapezian method integration (not used in kkintgr)
   REAL FUNCTION trapz(y,h,n)

      IMPLICIT NONE
   
      REAL,          INTENT(IN)     :: y(n)
   
      INTEGER,       INTENT(IN)     :: n
      REAL,          INTENT(IN)     :: h
   
   
      INTEGER i
   
      trapz = y(1)
      DO i = 2, n-1
   
         trapz = trapz + 2*y(i)
   
      ENDDO
      trapz = trapz + y(n)
   
      trapz = trapz*h/2.0
   
   END FUNCTION trapz

END MODULE m_kkintgr