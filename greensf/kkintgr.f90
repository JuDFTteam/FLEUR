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
   ! TODO: Look at FFT for Transformation
   !       How to do changing imaginary parts
   !------------------------------------------------------------------------------
   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   INTEGER, PARAMETER :: method_maclaurin = 1
   INTEGER, PARAMETER :: method_deriv     = 2
   INTEGER, PARAMETER :: method_direct    = 3
   INTEGER, PARAMETER :: method_fft       = 4

   CHARACTER(len=10), PARAMETER :: smooth_method = 'lorentzian' !(lorentzian or gaussian)

   CONTAINS

   SUBROUTINE kkintgr(im,eMesh,ez,l_conjg,g,method)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated
      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      !The dominant source of error for this routine is a insufficiently dense energy mesh on the real axis
      !TODO: Some way to estimate the error (maybe search for the sharpest peak and estimate from width)
      USE m_smooth
      USE m_lorentzian_smooth

      COMPLEX,       INTENT(IN)     :: im(:)       !Imaginary part of the green's function on the real axis
      REAL,          INTENT(IN)     :: eMesh(:)    !Energy grid on the real axis
      COMPLEX,       INTENT(IN)     :: ez(:)       !Complex energy contour
      LOGICAL,       INTENT(IN)     :: l_conjg     !Switch determines wether we calculate g on the complex conjugate of the contour ez
      COMPLEX,       INTENT(INOUT)  :: g(:)        !Green's function on the complex plane
      INTEGER,       INTENT(IN)     :: method      !Integer associated with the method to be used (definitions above)

      INTEGER  :: iz,izp,n1,n2,ne,nz
      INTEGER  :: ismooth,nsmooth
      REAL     :: eb,del
      COMPLEX  :: re_n1,re_n2,im_n1,im_n2
      INTEGER, ALLOCATABLE :: smoothInd(:)
      REAL,    ALLOCATABLE :: sigma(:)
      COMPLEX, ALLOCATABLE :: smoothed(:,:)

      nz  = SIZE(ez)
      ne  = SIZE(eMesh)
      eb  = eMesh(1)
      del = eMesh(2) - eMesh(1)

      ALLOCATE(smoothInd(nz),source=0)
      ALLOCATE(sigma(nz),source=0.0)

      IF(method.NE.method_direct) THEN
         CALL timestart("kkintgr: smoothing")
         !Smooth the imaginary part beforehand
         !Determine how many unique values there are for the imaginary part
         nsmooth = 0
         outer: DO iz = 1, nz
            DO izp = 1, iz-1
               IF(ABS(AIMAG(ez(izp))-AIMAG(ez(iz))).LT.1e-12) THEN
                  smoothInd(iz) = smoothInd(izp)
                  CYCLE outer
               ENDIF
            ENDDO
            nsmooth = nsmooth + 1
            smoothInd(iz) = nsmooth
            sigma(nsmooth) = AIMAG(ez(iz))
         ENDDO outer
         ALLOCATE(smoothed(ne,nsmooth), source=cmplx_0)
         !$OMP parallel do default(none) &
         !$OMP shared(nsmooth,smoothed,sigma,ne,eMesh,im) &
         !$OMP private(ismooth)
         DO ismooth = 1, nsmooth
            smoothed(:,ismooth) = im(:ne)
            IF(ABS(sigma(ismooth)).LT.1e-12) CYCLE
            SELECT CASE (TRIM(ADJUSTL(smooth_method)))
            CASE('lorentzian')
               CALL lorentzian_smooth(eMesh,smoothed(:,ismooth),sigma(ismooth),ne)
            CASE('gaussian')
               CALL smooth(eMesh,smoothed(:,ismooth),sigma(ismooth),ne)
            CASE DEFAULT
               CALL juDFT_error("No valid smooth_method set",&
                                hint="This is a bug in FLEUR, please report",&
                                calledby="kkintgr")
            END SELECT
         ENDDO
         !$OMP end parallel do
         CALL timestop("kkintgr: smoothing")
      ENDIF


      CALL timestart("kkintgr: integration")
      !$OMP parallel do default(none) &
      !$OMP shared(nz,ne,method,del,eb,l_conjg) &
      !$OMP shared(g,ez,eMesh,im,smoothed,smoothInd) &
      !$OMP private(iz,n1,n2,re_n1,re_n2,im_n1,im_n2)
      DO iz = 1, nz
         SELECT CASE(method)

         CASE(method_direct)
            g(iz) = kk_direct(im,eMesh,MERGE(conjg(ez(iz)),ez(iz),l_conjg))
         CASE(method_maclaurin, method_deriv)
            !Use the previously smoothed version and interpolate after
            !Next point to the left
            n1 = INT((REAL(ez(iz))-eb)/del) +1
            !next point to the right
            n2 = n1 + 1
            !Here we perform the Kramers-kronig-Integration
            re_n2 = kk_num(smoothed(:,smoothInd(iz)),ne,n2,method)
            re_n1 = kk_num(smoothed(:,smoothInd(iz)),ne,n1,method)
            !Interpolate to the energy ez(iz)
            !Real Part
            g(iz) = (re_n2-re_n1)/del * (REAL(ez(iz))-(n1-1)*del-eb) + re_n1

            !Imaginary Part (0 outside of the energy range)
            IF(n1.LE.ne.AND.n1.GE.1) THEN
               im_n1 = smoothed(n1,smoothInd(iz))
            ELSE
               im_n1 = 0.0
            ENDIF
            IF(n2.LE.ne.AND.n2.GE.1) THEN
               im_n2 = smoothed(n2,smoothInd(iz))
            ELSE
               im_n2 = 0.0
            ENDIF

            g(iz) = g(iz) + ImagUnit *( (im_n2-im_n1)/del * (REAL(ez(iz))-(n1-1)*del-eb) + im_n1 )

            IF(l_conjg) g(iz) = conjg(g(iz))

         CASE(method_fft)
            CALL juDFT_error("Not implemented yet", calledby="kkintgr")
         CASE DEFAULT
            CALL juDFT_error("Not a valid integration method", calledby="kkintgr")
         END SELECT
      ENDDO
      !$OMP end parallel do
      CALL timestop("kkintgr: integration")

   END SUBROUTINE kkintgr

   PURE COMPLEX FUNCTION kk_direct(im,eMesh,z)

      USE m_trapz

      COMPLEX, INTENT(IN) :: im(:)
      REAL,    INTENT(IN) :: eMesh(:)
      COMPLEX, INTENT(IN) :: z

      COMPLEX :: integrand(SIZE(eMesh))

      integrand = 1.0/(z-eMesh) * im
      kk_direct = -1/pi_const *trapz(integrand,eMesh(2)-eMesh(1),SIZE(eMesh))
   END FUNCTION kk_direct

   PURE COMPLEX FUNCTION kk_num(im,ne,ire,method)

      COMPLEX, INTENT(IN)  :: im(:) !Imaginary part
      INTEGER, INTENT(IN)  :: ne     !Dimension of the energy grid
      INTEGER, INTENT(IN)  :: ire    !Position where to calculate the real part
      INTEGER, INTENT(IN)  :: method !Method to be used
      INTEGER i,j
      COMPLEX y
      COMPLEX im_ire

      kk_num = cmplx_0
      IF(ire.LE.ne.AND.ire.GE.1) THEN
         im_ire = im(ire)
      ELSE
         im_ire = cmplx_0
      ENDIF

      SELECT CASE(method)

      CASE (method_maclaurin)
         !Calculate the real part on the same energy points as the imaginary part
         !regardless of the contour
         !If i is odd skip the odd points and the other way around and use the trapezian method
         DO i = MERGE(1,2,MOD(ire,2)==0), ne, 2
            y = - 1/pi_const * 2.0 * im(i)/REAL(ire-i)
            IF(i.EQ.1 .OR. i.EQ.2 .OR.&
               i.EQ.ne .OR. i.EQ.ne-1) y = y/2.0
            kk_num = kk_num + y
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
               ELSE IF((ire.LT.ne).AND.(ire.GT.1)) THEN
                  y = -1/pi_const * (im(ire+1)-im(ire-1))/2.0
               ENDIF
            ENDIF
            IF(j.EQ.1 .OR. j.EQ.ne) y = y/2.0
            kk_num = kk_num + y
         ENDDO
      CASE default
      END SELECT

   END FUNCTION kk_num

END MODULE m_kkintgr
