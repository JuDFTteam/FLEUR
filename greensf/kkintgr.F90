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
   USE m_types_mat
   USE m_smooth
   USE m_lorentzian_smooth

#include"cpp_double.h"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER :: method_maclaurin = 1
   INTEGER, PARAMETER :: method_deriv     = 2
   INTEGER, PARAMETER :: method_direct    = 3
   INTEGER, PARAMETER :: method_fft       = 4

   CHARACTER(len=10), PARAMETER :: smooth_method = 'lorentzian' !(lorentzian or gaussian)
   INTEGER, PARAMETER :: int_method(3) = [method_direct,method_direct,method_maclaurin]

   TYPE(t_mat),SAVE, ALLOCATABLE :: integration_weights(:)
   INTEGER,SAVE, ALLOCATABLE :: methods(:)
   REAL,SAVE, ALLOCATABLE :: sigma(:)
   REAL,SAVE, ALLOCATABLE :: energy_grid(:)
   REAL, SAVE :: del = -1.0

   INTERFACE kkintgr
      PROCEDURE :: kkintgr_single, kkintgr_mmpmat
      PROCEDURE :: kkintgr_mmpmat_extra1
   END INTERFACE

   PUBLIC kkintgr, kkintgr_init, kkintgr_free

   CONTAINS

   SUBROUTINE kkintgr_init(eMesh, contour, iContour, nContour, shape)

      REAL,          INTENT(IN)     :: eMesh(:)    !Energy grid on the real axis
      COMPLEX,       INTENT(IN)     :: contour(:)  !Complex energy contour
      INTEGER,       INTENT(IN)     :: iContour,nContour,shape

      INTEGER :: iz,izp,i,n1,n2
      REAL :: y

      IF(.NOT.ALLOCATED(integration_weights)) THEN
         ALLOCATE(integration_weights(nContour))
         ALLOCATE(methods(nContour), source=0)
         ALLOCATE(sigma(nContour), source=0.0)
         energy_grid = eMesh
      ENDIF

      IF(integration_weights(iContour)%allocated()) RETURN

      IF(del < 0.0) del = eMesh(2)-eMesh(1)
      IF(ABS(del-(eMesh(2)-eMesh(1)))>1e-12) CALL juDFT_error("Inconsistent energy grids", calledby="kkintgr_init")
      methods(iContour) = int_method(shape)

      CALL integration_weights(iContour)%init(.false.,SIZE(eMesh),SIZE(contour))

      IF(int_method(shape) == method_direct) THEN
         DO iz = 1, SIZE(contour)
            integration_weights(iContour)%data_c(:,iz) = 1.0/(contour(iz)-eMesh)
            integration_weights(iContour)%data_c(1,iz) = integration_weights(iContour)%data_c(1,iz)/2.0
            integration_weights(iContour)%data_c(SIZE(eMesh),iz) = integration_weights(iContour)%data_c(SIZE(eMesh),iz)/2.0
         ENDDO
      ELSE IF(int_method(shape) == method_maclaurin) THEN
         
         IF(ANY(ABS(AIMAG(contour)-AIMAG(contour(1)))>1e-12)) THEN
            CALL juDFT_error("Unsuitable contour for integration", calledby="kkintgr_init")
         ENDIF

         sigma(iContour) = AIMAG(contour(1))

         CALL integration_weights(iContour)%init(.false.,SIZE(eMesh),2*SIZE(contour))

         DO iz = 1, SIZE(contour)
            !Next point to the left of the point
            n1 = INT((REAL(contour(iz))-eMesh(1))/del) + 1

            integration_weights(iContour)%data_c(:,iz) = cmplx_0
            !Calculate the real part on the same energy points as the imaginary part
            !regardless of the contour
            !If i is odd skip the odd points and the other way around and use the trapezian method
            DO i = MERGE(1,2,MOD(n1,2)==0), SIZE(eMesh), 2
               y = -1/pi_const * 2.0/REAL(n1-i)
               IF(i.EQ.1 .OR. i.EQ.2 .OR.&
                  i.EQ.SIZE(eMesh) .OR. i.EQ.SIZE(eMesh)-1) y = y/2.0
               integration_weights(iContour)%data_c(i,iz) = y
            ENDDO
            IF(n1>=1.AND.n1<= SIZE(eMesh)) THEN
               integration_weights(iContour)%data_c(n1,iz) = ImagUnit
            ENDIF

            integration_weights(iContour)%data_c(:,iz) = integration_weights(iContour)%data_c(:,iz) * (1.0-(REAL(contour(iz))-(n1-1)*del-eMesh(1))/del)

         ENDDO

         DO izp = SIZE(contour)+1, 2*SIZE(contour)
            !Next point to the right of the point
            iz=izp-SIZE(contour)
            n2 = INT((REAL(contour(iz))-eMesh(1))/del) + 2

            integration_weights(iContour)%data_c(:,izp) = cmplx_0
            !Calculate the real part on the same energy points as the imaginary part
            !regardless of the contour
            !If i is odd skip the odd points and the other way around and use the trapezian method
            DO i = MERGE(1,2,MOD(n2,2)==0), SIZE(eMesh), 2
               y = -1/pi_const * 2.0/REAL(n2-i)
               IF(i.EQ.1 .OR. i.EQ.2 .OR.&
                  i.EQ.SIZE(eMesh) .OR. i.EQ.SIZE(eMesh)-1) y = y/2.0
               integration_weights(iContour)%data_c(i,izp) = y
            ENDDO
            IF(n2>=1.AND.n2<=SIZE(eMesh)) THEN
               integration_weights(iContour)%data_c(n2,izp) = ImagUnit
            ENDIF
            integration_weights(iContour)%data_c(:,izp) = integration_weights(iContour)%data_c(:,izp) * (REAL(contour(iz))-(n2-2)*del-eMesh(1))/del
         ENDDO

      ELSE
         CALL juDFT_error("Not a supported method", calledby="kkintgr_init")
      ENDIF

   END SUBROUTINE kkintgr_init

   SUBROUTINE kkintgr_free()

      INTEGER :: iContour

      IF(.NOT.ALLOCATED(integration_weights)) RETURN

      DO iContour = 1, SIZE(integration_weights)
         IF(integration_weights(iContour)%allocated()) CALL integration_weights(iContour)%free()
      ENDDO

      DEALLOCATE(integration_weights,methods,energy_grid,sigma)

   END SUBROUTINE kkintgr_free


   SUBROUTINE kkintgr_single(im,l_conjg,g,iContour)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated
      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      !The dominant source of error for this routine is a insufficiently dense energy mesh on the real axis
      !TODO: Some way to estimate the error (maybe search for the sharpest peak and estimate from width)

      COMPLEX,       INTENT(IN)     :: im(:)       !Imaginary part of the green's function on the real axis
      LOGICAL,       INTENT(IN)     :: l_conjg     !Switch determines wether we calculate g on the complex conjugate of the contour ez
      COMPLEX,       INTENT(INOUT)  :: g(:)        !Green's function on the complex plane
      INTEGER,       INTENT(IN)     :: iContour    !Which contour is used

      INTEGER  :: ne,nz
      COMPLEX, ALLOCATABLE :: smoothed(:)
      CHARACTER(len=1) :: transA

      IF(.NOT.integration_weights(iContour)%allocated()) THEN
         CALL juDFT_error("Integration weights not initialized",&
                          hint="This is a bug in FLEUR, please report",&
                          calledby="kkintgr")
      ENDIF

      nz  = integration_weights(iContour)%matsize2
      ne  = integration_weights(iContour)%matsize1

      IF(methods(iContour)==method_direct) THEN
         transA = 'T'
         IF(l_conjg) transA = 'C'

         CALL CPP_BLAS_cgemm(transA,'N',&
                           nz,1,ne,&
                           CMPLX(-del/pi_const,0.0),&
                           integration_weights(iContour)%data_c,ne,&
                           im,ne,&
                           cmplx_0,&
                           g,nz)
      ELSE IF(methods(iContour)==method_maclaurin) THEN

         nz  = INT(integration_weights(iContour)%matsize2/2)

         smoothed = im
         IF(ABS(sigma(iContour)).GT.1e-12) THEN
            CALL timestart('kkintgr: smoothing')
            SELECT CASE (TRIM(ADJUSTL(smooth_method)))
            CASE('lorentzian')
               CALL lorentzian_smooth(energy_grid,smoothed,sigma(iContour),ne)
            CASE('gaussian')
               CALL smooth(energy_grid,smoothed,sigma(iContour),ne)
            CASE DEFAULT
               CALL juDFT_error("No valid smooth_method set",&
                              hint="This is a bug in FLEUR, please report",&
                              calledby="kkintgr")
            END SELECT

            CALL timestop('kkintgr: smoothing')
         ENDIF

         CALL CPP_BLAS_cgemm('T','N',&
                             nz,1,ne,&
                             cmplx_1,&
                             integration_weights(iContour)%data_c(:,:nz),ne,&
                             smoothed,ne,&
                             cmplx_0,&
                             g,nz)         
         CALL CPP_BLAS_cgemm('T','N',&
                             nz,1,ne,&
                             cmplx_1,&
                             integration_weights(iContour)%data_c(:,nz+1:),ne,&
                             smoothed,ne,&
                             cmplx_1,&
                             g,nz)

         IF(l_conjg) g = conjg(g)

      ENDIF


   END SUBROUTINE kkintgr_single

   SUBROUTINE kkintgr_mmpmat(im,l_conjg,g,iContour)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated
      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      !The dominant source of error for this routine is a insufficiently dense energy mesh on the real axis
      !TODO: Some way to estimate the error (maybe search for the sharpest peak and estimate from width)
      USE m_smooth
      USE m_lorentzian_smooth

      COMPLEX,       INTENT(IN)     :: im(:,-lmaxU_const:,-lmaxU_const:)       !Imaginary part of the green's function on the real axis
      LOGICAL,       INTENT(IN)     :: l_conjg     !Switch determines wether we calculate g on the complex conjugate of the contour ez
      COMPLEX,       INTENT(INOUT)  :: g(:,-lmaxU_const:,-lmaxU_const:)        !Green's function on the complex plane
      INTEGER,       INTENT(IN)     :: iContour      !Which contour is used

      INTEGER  :: ne,nz,m,mp
      COMPLEX, ALLOCATABLE :: smoothed(:,:,:)
      CHARACTER(len=1) :: transA

      IF(.NOT.integration_weights(iContour)%allocated()) THEN
         CALL juDFT_error("Integration weights not initialized",&
                          hint="This is a bug in FLEUR, please report",&
                          calledby="kkintgr")
      ENDIF

      nz  = integration_weights(iContour)%matsize2
      ne  = integration_weights(iContour)%matsize1

      transA = 'T'
      IF(l_conjg) transA = 'C'

      IF(methods(iContour)==method_direct) THEN
         CALL CPP_BLAS_cgemm(transA,'N',&
                           nz,(2*lmaxU_const+1)**2,ne,&
                           CMPLX(-del/pi_const,0.0),&
                           integration_weights(iContour)%data_c,ne,&
                           im,ne,&
                           cmplx_0,&
                           g,nz)
      ELSE IF(methods(iContour)==method_maclaurin) THEN

         nz  = INT(integration_weights(iContour)%matsize2/2)

         ALLOCATE(smoothed, source=im)
         IF(ABS(sigma(iContour)).GT.1e-12) THEN
            CALL timestart('kkintgr: smoothing')
            SELECT CASE (TRIM(ADJUSTL(smooth_method)))
            CASE('lorentzian')
               !$OMP parallel do default(none) collapse(2) &
               !$OMP shared(energy_grid, smoothed,sigma,ne,iContour) &
               !$OMP private(m,mp)
               DO mp = -lmaxU_const, lmaxU_const
                  DO m = -lmaxU_const, lmaxU_const
                     CALL lorentzian_smooth(energy_grid,smoothed(:,m,mp),sigma(iContour),ne)
                  ENDDO
               ENDDO
               !$OMP end parallel do
            CASE('gaussian')
               !$OMP parallel do default(none) collapse(2) &
               !$OMP shared(energy_grid, smoothed,sigma,ne,iContour) &
               !$OMP private(m,mp)
               DO mp = -lmaxU_const, lmaxU_const
                  DO m = -lmaxU_const, lmaxU_const
                     CALL smooth(energy_grid,smoothed(:,m,mp),sigma(iContour),ne)
                  ENDDO
               ENDDO
               !$OMP end parallel do
            CASE DEFAULT
               CALL juDFT_error("No valid smooth_method set",&
                              hint="This is a bug in FLEUR, please report",&
                              calledby="kkintgr")
            END SELECT

            CALL timestop('kkintgr: smoothing')
         ENDIF

         CALL CPP_BLAS_cgemm('T','N',&
                              nz,(2*lmaxU_const+1)**2,ne,&
                              cmplx_1,&
                              integration_weights(iContour)%data_c(:,:nz),ne,&
                              smoothed,ne,&
                              cmplx_0,&
                              g,nz)         
         CALL CPP_BLAS_cgemm('T','N',&
                              nz,(2*lmaxU_const+1)**2,ne,&
                              cmplx_1,&
                              integration_weights(iContour)%data_c(:,nz+1:),ne,&
                              smoothed,ne,&
                              cmplx_1,&
                              g,nz)

         IF(l_conjg) g = conjg(g)

      ENDIF


   END SUBROUTINE kkintgr_mmpmat

   SUBROUTINE kkintgr_mmpmat_extra1(im,l_conjg,g,iContour)

      !calculates the Kramer Kronig Transformation on the same contour where the imaginary part was calculated
      !Re(G(E+i * delta)) = -1/pi * int_bot^top dE' P(1/(E-E')) * Im(G(E'+i*delta))

      !The dominant source of error for this routine is a insufficiently dense energy mesh on the real axis
      !TODO: Some way to estimate the error (maybe search for the sharpest peak and estimate from width)
      USE m_smooth
      USE m_lorentzian_smooth

      COMPLEX,       INTENT(IN)     :: im(:,-lmaxU_const:,-lmaxU_const:,:)       !Imaginary part of the green's function on the real axis
      LOGICAL,       INTENT(IN)     :: l_conjg     !Switch determines wether we calculate g on the complex conjugate of the contour ez
      COMPLEX,       INTENT(INOUT)  :: g(:,-lmaxU_const:,-lmaxU_const:,:)        !Green's function on the complex plane
      INTEGER,       INTENT(IN)     :: iContour      !Which contour is used

      INTEGER  :: ne,nz,m,mp,i
      COMPLEX, ALLOCATABLE :: smoothed(:,:,:,:)
      CHARACTER(len=1) :: transA

      IF(.NOT.integration_weights(iContour)%allocated()) THEN
         CALL juDFT_error("Integration weights not initialized",&
                          hint="This is a bug in FLEUR, please report",&
                          calledby="kkintgr")
      ENDIF

      nz  = integration_weights(iContour)%matsize2
      ne  = integration_weights(iContour)%matsize1

      IF(methods(iContour)==method_direct) THEN
         CALL CPP_BLAS_cgemm(transA,'N',&
                           nz,(2*lmaxU_const+1)**2*SIZE(im,4),ne,&
                           CMPLX(-del/pi_const,0.0),&
                           integration_weights(iContour)%data_c,ne,&
                           im,ne,&
                           cmplx_0,&
                           g,nz)
      ELSE IF(methods(iContour)==method_maclaurin) THEN

         nz  = INT(integration_weights(iContour)%matsize2/2)

         ALLOCATE(smoothed, source=im)
         IF(ABS(sigma(iContour)).GT.1e-12) THEN
            CALL timestart('kkintgr: smoothing')
            SELECT CASE (TRIM(ADJUSTL(smooth_method)))
            CASE('lorentzian')
               !$OMP parallel do default(none) collapse(3) &
               !$OMP shared(energy_grid, smoothed,sigma,ne,iContour) &
               !$OMP private(m,mp,i)
               DO i = 1, SIZE(smoothed,4)
                  DO mp = -lmaxU_const, lmaxU_const
                     DO m = -lmaxU_const, lmaxU_const
                        CALL lorentzian_smooth(energy_grid,smoothed(:,m,mp,i),sigma(iContour),ne)
                     ENDDO
                  ENDDO
               ENDDO
               !$OMP end parallel do
            CASE('gaussian')
               !$OMP parallel do default(none) collapse(3) &
               !$OMP shared(energy_grid, smoothed,sigma,ne,iContour) &
               !$OMP private(m,mp)
               DO i = 1, SIZE(smoothed,4)
                  DO mp = -lmaxU_const, lmaxU_const
                     DO m = -lmaxU_const, lmaxU_const
                        CALL smooth(energy_grid,smoothed(:,m,mp,i),sigma(iContour),ne)
                     ENDDO
                  ENDDO
               ENDDO
               !$OMP end parallel do
            CASE DEFAULT
               CALL juDFT_error("No valid smooth_method set",&
                              hint="This is a bug in FLEUR, please report",&
                              calledby="kkintgr")
            END SELECT

            CALL timestop('kkintgr: smoothing')
         ENDIF

         CALL CPP_BLAS_cgemm('T','N',&
                              nz,(2*lmaxU_const+1)**2*SIZE(smoothed,4),ne,&
                              cmplx_1,&
                              integration_weights(iContour)%data_c(:,:nz),ne,&
                              smoothed,ne,&
                              cmplx_0,&
                              g,nz)         
         CALL CPP_BLAS_cgemm('T','N',&
                              nz,(2*lmaxU_const+1)**2*SIZE(smoothed,4),ne,&
                              cmplx_1,&
                              integration_weights(iContour)%data_c(:,nz+1:),ne,&
                              smoothed,ne,&
                              cmplx_1,&
                              g,nz)

         IF(l_conjg) g = conjg(g)

      ENDIF

   END SUBROUTINE kkintgr_mmpmat_extra1

END MODULE m_kkintgr