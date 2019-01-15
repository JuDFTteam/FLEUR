MODULE m_kkintgr

CONTAINS

   SUBROUTINE kkintgr_real(n_cut,sigma,del,im,g,bot)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated

      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      USE m_types
      USE m_constants
      USE m_intgr

      IMPLICIT NONE

      COMPLEX,             INTENT(OUT)    :: g(:)  !green's function
      !COMPLEX,             INTENT(IN)     :: e(:)  !energy contour, on which to calculate G(z)
      REAL,                INTENT(IN)     :: im(:) !previously calculated imaginary part on equidistant mesh at E+i*sigma

      !INTEGER,             INTENT(IN)     :: n_z
      INTEGER,             INTENT(IN)     :: n_cut
      REAL,                INTENT(IN)     :: sigma       
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: bot

      !Local declarations:
      INTEGER iz, j
      COMPLEX  ez
      REAL, ALLOCATABLE :: integrand(:)
      REAL result

      ALLOCATE (integrand(n_cut))


      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(n_cut,sigma,del) &
      !$OMP SHARED(g,im) &
      !$OMP PRIVATE(j,ez,integrand,result)

      !$OMP DO
      DO iz = 1, n_cut

         !We calculate the integral by removing the singularity and treating it analytically

         DO j = 1, iz -1

            !Formula for the integrand if the contour is on the same imaginary part as before
            integrand(j) = (im(j)-im(iz))/((j-iz) * del)

         ENDDO

         !at j = iz we estimate the integrand with the centered three-point formula for the derivative of im(e) 
         !(at the edges we use backward/forward formulas)

         IF(iz.EQ.1) THEN 
            integrand(iz) = (im(2)-im(1))/(del) 
         ELSE IF(iz.EQ.n_cut) THEN
            integrand(iz) = (im(n_cut)-im(n_cut-1))/del
         ELSE
            integrand(iz) = (im(iz+1)-im(iz-1))/(2*del) 
         END IF

         DO j = iz+1, n_cut

            !Formula for the integrand if the contour is on the same imaginary part as before
            integrand(j) = (im(j)-im(iz))/((j-iz) * del)

         ENDDO

         CALL intgz0(integrand(:), del,n_cut,result,.false.)

         !The singularity is treated analytically

         g(iz) = -1/pi_const * (result - im(iz) * LOG(ABS(REAL(iz)/REAL(iz-1 - n_cut)))) + ImagUnit * im(iz)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL


   ENDSUBROUTINE kkintgr_real

   SUBROUTINE kkintgr_complex(nz,e,n_cut,sigma,del,im,g,bot)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated

      !G(z) = -1/pi * int_bot^top dE' 1/(z-E') * Im(G(E'+i*delta))

      USE m_types
      USE m_constants
      USE m_intgr

      IMPLICIT NONE

      COMPLEX,             INTENT(OUT)    :: g(:)  !green's function
      COMPLEX,             INTENT(IN)     :: e(:)  !energy contour, on which to calculate G(z)
      REAL,                INTENT(IN)     :: im(:) !previously calculated imaginary part on equidistant mesh at E+i*sigma

      INTEGER,             INTENT(IN)     :: nz
      INTEGER,             INTENT(IN)     :: n_cut
      REAL,                INTENT(IN)     :: sigma       
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: bot

      INTEGER i_sing, iz, j
      COMPLEX ez



      DO iz = 1, nz

         ez = e(iz) - bot

         i_sing = INT(ez/del) + 1

         g(:) = CMPLX(0.0,0.0) !TEMPORARY
      ENDDO






   END SUBROUTINE kkintgr_complex


   SUBROUTINE energy_contour(mode,g,ef,l_ef)

      ! calculates the energy contour where the greens function is calculated
      ! mode determines the kind of contour between e_bot and the fermi energy (if l_ef = .true.) 
      ! mode = 1 gives a equidistant contour with imaginary part g%sigma with g%nz points

      ! mode = 2 gives a half circle with 2**g%nz points

      USE m_types
      USE m_constants
      USE m_juDFT

      INTEGER,          INTENT(IN)     :: mode
      TYPE(t_greensf),  INTENT(INOUT)  :: g
      REAL,             INTENT(IN)     :: ef
      LOGICAL,          INTENT(IN)     :: l_ef


      INTEGER i, j, iz

      REAL e1, e2, del
      REAL psi(4), wpsi(4), r, xr, xm, c, s, a, b



      IF(mode.EQ.1) THEN

         e1 = g%e_bot
         e2 = g%e_top

         IF(l_ef) e2 = ef

         del = (e2-e1)/REAL(g%nz-1)

         DO i = 1, g%nz

            g%e(i) = (i-1)*del + e1 + ImagUnit * g%sigma

         ENDDO

         g%de(:) = del

      ELSE IF(mode.EQ.2) THEN

         np = 2**(g%nz-2)

         e1 = g%e_bot
         e2 = g%e_top

         IF(l_ef) e2 = ef


         !Radius
         r  = (e2-e1)*0.5
         xr = (e2+e1)*0.5

         a = 0.43056815579702629
         b = 0.16999052179242813

         psi(1) =    a/np
         psi(2) =    b/np
         psi(3) =   -b/np
         psi(4) =   -a/np

         a = 0.17392742256872693
         b = 0.32607257743127307

         wpsi(1) =   a/np
         wpsi(2) =   b/np
         wpsi(3) =   b/np
         wpsi(4) =   a/np

         iz = 1

         DO i = 1, np

            xm = (np-i+0.5)/np

            DO j = 1, 4

               c = cos((psi(j)+xm)**2*pi_const)
               s = sin((psi(j)+xm)**2*pi_const)

               g%e(iz) = CMPLX(xr+r*c, r*s*0.25)

               g%de(iz) = pi_const * CMPLX((psi(j)+xm)*r*s*wpsi(j)*2.0,&
                                          -(psi(j)+xm)*r*c*wpsi(j)*0.5)

               iz = iz+1

            ENDDO

         ENDDO

         g%nz = 2**(g%nz)

      ELSE

         CALL juDFT_error("Invalid mode for energy contour in Green's function calculation", calledby="energy_contour")

      END IF



   END SUBROUTINE energy_contour

ENDMODULE m_kkintgr 