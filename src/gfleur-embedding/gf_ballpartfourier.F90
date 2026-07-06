!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_ballpartfourier
      use m_juDFT
      IMPLICIT NONE
!**********************************************************             
!     Calculate the Fourier transform of part of a sphere.              
!     Frank Freimuth, April 2007                                        
!                                                                       
!     NEEDS ADJUSTMENT FOR DIFFERENT COMPILERS                          
!     see functions at end of file                                      
!                                                                       
!**********************************************************             
      REAL :: radius,kpara,korthogo 
      CONTAINS 
      SUBROUTINE ballpartfou(radius_in,zcoor,kpara_in,                  &
     &                          korthogo_in,value)                      
      USE m_quadrature 
      USE m_constants, ONLY: pi_const 
      IMPLICIT NONE 
      REAL,INTENT(IN) :: radius_in 
      REAL,INTENT(IN) :: zcoor 
      REAL,INTENT(IN) :: kpara_in 
      REAL,INTENT(IN) :: korthogo_in 
      COMPLEX,INTENT(OUT) :: value 
                                                                        
      REAL :: value_r 
      REAL :: value_i 
      REAL :: pi 
                                                                        
      REAL    :: epsabs,epsrel,abserr,sinR,cosR,sinz,cosz,kp1,kp2,kp3 
      INTEGER :: neval,qerr,limit,lenw,last 
      REAL,ALLOCATABLE :: work(:) 
      INTEGER,ALLOCATABLE :: iwork(:) 
                                                                        
      pi = pi_const 
      radius = radius_in 
      kpara = kpara_in 
      korthogo = korthogo_in 
                                                                        
                                        !analytic formula for korthogo =
      IF(ABS(korthogo) < 1.0E-10) THEN 
         sinR = SIN(kpara*Radius) 
         cosR = COS(kpara*Radius) 
         sinz = SIN(kpara*zcoor) 
         cosz = COS(kpara*zcoor) 
         kp1 = 1.0/kpara 
         kp2 = kp1*kp1 
         kp3 = kp2*kp1 
         value = pi*( 2.0*kp3*CMPLX(sinR,-cosR)+                        &
     &        2.0*kp3*CMPLX(-sinz,cosz)+                                &
     &        zcoor*zcoor*kp1*CMPLX(sinz,-cosz)+                        &
     &        radius**2  *kp1*CMPLX(-sinz,cosz)+                        &
     &        2*radius*kp2*CMPLX(-cosR,-sinR)+                          &
     &        2*zcoor* kp2*CMPLX( cosz ,sinz)  )                        
                                                                        
                                                                        
         RETURN 
      ENDIF 
                                                                        
      epsabs = 1.0E-12 
      epsrel = 1.0E-12 
      limit = 5 
      lenw = 4*limit 
      ALLOCATE(iwork(limit)) 
      ALLOCATE(work(lenw)) 
      CALL dqag(integrandr,zcoor,radius,epsabs,epsrel,6,value_r,abserr, &
     &     neval,qerr,limit,lenw,last,iwork,work)                       
                                                                        
      IF(qerr /= 0)THEN 
         PRINT*,"error in dqag" 
         PRINT*,"qerr =",qerr 
          CALL juDFT_error("dqag",calledby="gf_ballpartfourier.F90")
      ENDIF 
      CALL dqag(integrandi,zcoor,radius,epsabs,epsrel,6,value_i,abserr, &
     &     neval,qerr,limit,lenw,last,iwork,work)                       
                                                                        
      IF(qerr /= 0)THEN 
         PRINT*,"error in dqag" 
         PRINT*,"qerr =",qerr 
          CALL juDFT_error("dqag",calledby="gf_ballpartfourier.F90")
      ENDIF 
                                                                        
                                                                        
      value = CMPLX(value_r,value_i)*2*pi/korthogo 
                                                                        
                                                                        
      END SUBROUTINE ballpartfou 
                                                                        
                                                                        
!two functions that will be integrated follow                           
! both need bessel-functions that need some compiler dependend call     
                                                                        
      FUNCTION integrandr(x) RESULT(integrand) 
      USE m_gf_bessel 
                                                                        
      IMPLICIT NONE 
      REAL,INTENT(IN) :: x 
                                                                        
      REAL :: integrand 
      REAL :: r2mz2 
                                                                        
                                                                        
      r2mz2 = radius*radius-x*x 
      r2mz2 = SQRT(r2mz2) 
                                                                        
      integrand = COS(kpara*x)*r2mz2*gf_bessel1(r2mz2*korthogo) 
                                                                        
                                                                        
      RETURN 
      END FUNCTION integrandr 
                                                                        
      FUNCTION integrandi(x) RESULT(integrand) 
      USE m_gf_bessel 
                                                                        
      IMPLICIT NONE 
      REAL,INTENT(IN) :: x 
                                                                        
      REAL :: integrand 
      REAL :: r2mz2 
                                                                        
      r2mz2 = radius*radius-x*x 
      r2mz2 = SQRT(r2mz2) 
                                                                        
      integrand = SIN(kpara*x)*r2mz2*gf_bessel1(r2mz2*korthogo) 
                                                                        
                                                                        
      RETURN 
      END FUNCTION integrandi 
                                                                        
                                                                        
      END                                           
