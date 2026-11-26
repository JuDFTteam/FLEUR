!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_curvypartfourier 
      use m_juDFT
      IMPLICIT NONE
!**********************************************************             
!     Calculate the Fourier transform of the surface of                 
!     part of a sphere.                                                 
!     Frank Freimuth, October 2007                                      
!**********************************************************             
      REAL radius,kpara,korthogo,zcoor 
      CONTAINS 
      SUBROUTINE gf_curvypartfou(                                       &
     &              radius_in,circleradius,kpara_in,                    &
     &                          korthogo_in,zcoor_in,                   &
     &                          value)                                  
      USE m_quadrature 
      IMPLICIT NONE 
      REAL,INTENT(IN)     :: radius_in 
      REAL,INTENT(IN)     :: circleradius 
      REAL,INTENT(IN)     :: kpara_in 
      REAL,INTENT(IN)     :: korthogo_in 
      REAL,INTENT(IN)     :: zcoor_in 
      COMPLEX,INTENT(OUT) :: value 
                                                                        
      REAL value_r 
      REAL value_i 
                                                                        
      REAL epsabs,epsrel,abserr,sinR,cosR,sinz,cosz,kp1,kp2,kp3 
      INTEGER neval,qerr,limit,lenw,last 
      REAL,ALLOCATABLE    :: work(:) 
      INTEGER,ALLOCATABLE :: iwork(:) 
                                                                        
      radius=radius_in 
      kpara=kpara_in 
      korthogo=korthogo_in 
      zcoor=zcoor_in 
                                                                        
                                                                        
!      if(abs(kpara) .lt. 1.0e-10) then  !analytic formula for kpara=0  
!       return                                                          
!      endif                                                            
                                                                        
      epsabs=1.0E-12 
      epsrel=1.0E-12 
      limit=5 
      lenw=4*limit 
      ALLOCATE(iwork(limit)) 
      ALLOCATE(work(lenw)) 
      CALL dqag(integrandr,0.0,circleradius,epsabs,epsrel,6,value_r,    &
     &          abserr,neval,qerr,limit,lenw,last,iwork,work)           
                                                                        
      IF(qerr/=0)THEN 
         PRINT*,"error in dqag" 
         PRINT*,"qerr=",qerr 
          CALL juDFT_error("dqag",calledby="gf_curvypartfourier.F90")
      ENDIF 
!      print*,"neval=",neval                                            
!      print*,"abserr=",abserr                                          
!      print*,"last=",last                                              
      CALL dqag(integrandi,0.0,circleradius,epsabs,epsrel,6,value_i,    &
     &          abserr,neval,qerr,limit,lenw,last,iwork,work)           
      DEALLOCATE(iwork,work) 
      IF(qerr/=0)THEN 
         PRINT*,"error in dqag" 
         PRINT*,"qerr=",qerr 
          CALL juDFT_error("dqag",calledby="gf_curvypartfourier.F90")
      ENDIF 
!      print*,"neval=",neval                                            
!      print*,"abserr=",abserr                                          
!      print*,"last=",last                                              
                                                                        
      value=cmplx(value_r,value_i) 
                                                                        
      END SUBROUTINE gf_curvypartfou 
                                                                        
      FUNCTION integrandr(x) RESULT(integrand) 
      USE m_gf_bessel 
                                                                        
      IMPLICIT NONE 
      REAL,INTENT(IN) :: x 
                                                                        
      REAL :: integrand 
      REAL :: r2mz2 
      r2mz2=radius*radius-x*x 
      r2mz2=sqrt(r2mz2) 
                                                                        
                                                                        
      integrand = COS(korthogo*r2mz2)*x*gf_bessel0(x*kpara)*SQRT(1+(x   &
     &     /r2mz2)**2)                                                  
      RETURN 
      END FUNCTION integrandr 
                                                                        
      FUNCTION integrandi(x) RESULT(integrand) 
      USE m_gf_bessel 
      IMPLICIT NONE 
      REAL,INTENT(IN)::x 
                                                                        
      REAL :: integrand 
      REAL :: r2mz2 
                                                                        
      r2mz2=radius*radius-x*x 
      r2mz2=sqrt(r2mz2) 
                                                                        
      integrand = SIN(korthogo*r2mz2)*x*gf_bessel0(x*kpara)             &
     &     *SQRT(1+(x/r2mz2)**2)                                        
      RETURN 
      END FUNCTION integrandi 
                                                                        
      END                                           
