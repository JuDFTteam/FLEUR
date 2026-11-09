!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_ballpartfourier
      use m_juDFT
      use m_grule, only: grule
      use m_gf_bessel, only: gf_bessel1
      use m_quadrature, only: dqag
      use m_constants, only: pi_const
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: ballpartfou, ballpart_step
      PUBLIC :: n_ballpart_calls, n_ballpart_quad, n_ballpart_neval
!**********************************************************             
!     Calculate the Fourier transform of part of a sphere.              
!     Frank Freimuth, April 2007                                        
!                                                                       
!     NEEDS ADJUSTMENT FOR DIFFERENT COMPILERS                          
!     see functions at end of file                                      
!                                                                       
!**********************************************************             
      REAL :: radius,kpara,korthogo 
      !counters for diagnostics of the step-function setup
      INTEGER :: n_ballpart_calls = 0, n_ballpart_quad = 0, n_ballpart_neval = 0
      CONTAINS 
      SUBROUTINE ballpartfou(radius_in,zcoor,kpara_in,                  &
     &                          korthogo_in,value)                      
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
      n_ballpart_calls = n_ballpart_calls + 1
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
      n_ballpart_quad = n_ballpart_quad + 1
      n_ballpart_neval = n_ballpart_neval + neval
                                                                        
      IF(qerr /= 0)THEN 
         PRINT*,"error in dqag" 
         PRINT*,"qerr =",qerr 
          CALL juDFT_error("dqag",calledby="gf_ballpartfourier.F90")
      ENDIF 
      CALL dqag(integrandi,zcoor,radius,epsabs,epsrel,6,value_i,abserr, &
     &     neval,qerr,limit,lenw,last,iwork,work)                       
      n_ballpart_quad = n_ballpart_quad + 1
      n_ballpart_neval = n_ballpart_neval + neval
                                                                        
      IF(qerr /= 0)THEN 
         PRINT*,"error in dqag" 
         PRINT*,"qerr =",qerr 
          CALL juDFT_error("dqag",calledby="gf_ballpartfourier.F90")
      ENDIF 
                                                                        
                                                                        
      value = CMPLX(value_r,value_i)*2*pi/korthogo 
                                                                        
                                                                        
      END SUBROUTINE ballpartfou 

      SUBROUTINE ballpart_step(omtil,bmat,rmt,zcoor,taual,mx1,mx2,mx3,step)
!**********************************************************
!     Add the Fourier coefficients (k /= 0) of the part of a
!     sphere beyond a plane to step(k1,k2,k3).
!     zcoor is the signed distance of the plane from the
!     sphere centre taual (lattice coordinates).
!     With x = R cos(t) the integral over the cap is smooth
!     and a fixed Gauss-Legendre rule serves all k at once.
!**********************************************************
      IMPLICIT NONE
      REAL,INTENT(IN)       :: omtil,bmat(3,3),rmt,zcoor,taual(3)
      INTEGER,INTENT(IN)    :: mx1,mx2,mx3
      COMPLEX,INTENT(INOUT) :: step(-mx1:mx1,-mx2:mx2,-mx3:mx3)

      INTEGER :: k1,k2,k3,j,nq,ng
      REAL    :: z0,tpi,theta0,kmax,kp,kz,rg(3)
      LOGICAL :: l_ortho
      COMPLEX :: v,pre
      REAL,ALLOCATABLE    :: xg(:),wg(:),xq(:),sq(:),wq(:),a(:)
      COMPLEX,ALLOCATABLE :: e(:,:),ph1(:),ph2(:),ph3(:),row(:)

      tpi = 2.0*pi_const
      z0 = ABS(zcoor)
      IF (z0 >= rmt) RETURN
      theta0 = ACOS(z0/rmt)

      !bound on |k|, number of quadrature points
      kmax = mx1*SQRT(SUM(bmat(1,:)**2)) + mx2*SQRT(SUM(bmat(2,:)**2)) &
     &     + mx3*SQRT(SUM(bmat(3,:)**2))
      nq = CEILING(1.5*kmax*rmt) + 30
      ng = (nq+1)/2
      ALLOCATE(xg(ng),wg(ng),xq(nq),wq(nq),sq(nq),a(nq))
      CALL grule(nq,xg,wg)
      DO j = 1,ng
         xq(j) = xg(j);  xq(nq+1-j) = -xg(j)
         wq(j) = wg(j);  wq(nq+1-j) = wg(j)
      ENDDO
      !theta in [0,theta0]; integrand weight rho**2 (rho = R sin(theta))
      DO j = 1,nq
         xq(j) = 0.5*theta0*(1.0+xq(j))
         sq(j) = rmt*SIN(xq(j))
         wq(j) = 0.5*theta0*wq(j)*sq(j)*sq(j)
         xq(j) = rmt*COS(xq(j))
      ENDDO

      !phases of the sphere position
      ALLOCATE(ph1(-mx1:mx1),ph2(-mx2:mx2),ph3(-mx3:mx3))
      DO k1 = -mx1,mx1
         ph1(k1) = EXP(CMPLX(0.0,-tpi*k1*taual(1)))
      ENDDO
      DO k2 = -mx2,mx2
         ph2(k2) = EXP(CMPLX(0.0,-tpi*k2*taual(2)))
      ENDDO
      DO k3 = -mx3,mx3
         ph3(k3) = EXP(CMPLX(0.0,-tpi*k3*taual(3)))
      ENDDO

      !k_z depends on k3 only and k_perp on (k1,k2) only for a z-axis normal to the plane
      l_ortho = MAXVAL(ABS((/bmat(1,3),bmat(2,3),bmat(3,1),bmat(3,2)/))) < 1.0E-12

      IF (l_ortho) THEN
         ALLOCATE(e(nq,-mx3:mx3),row(-mx3:mx3))
         DO k3 = -mx3,mx3
            kz = k3*bmat(3,3)
            DO j = 1,nq
               e(j,k3) = EXP(CMPLX(0.0,kz*xq(j)))
            ENDDO
         ENDDO
         DO k2 = -mx2,mx2
            DO k1 = -mx1,mx1
               rg(1:2) = k1*bmat(1,1:2) + k2*bmat(2,1:2)
               kp = SQRT(rg(1)**2+rg(2)**2)
               IF (kp < 1.0E-10) THEN
                  DO k3 = -mx3,mx3
                     IF (k3 == 0) THEN
                        row(k3) = 0.0
                        CYCLE
                     ENDIF
                     CALL ballpartfou(rmt,z0,k3*bmat(3,3),0.0,v)
                     row(k3) = v
                  ENDDO
               ELSE
                  DO j = 1,nq
                     a(j) = wq(j)*gf_bessel1(kp*sq(j))*2.0*pi_const/kp
                  ENDDO
                  row = MATMUL(CMPLX(a,0.0),e)
               ENDIF
               pre = ph1(k1)*ph2(k2)/omtil
               IF (zcoor >= 0) THEN
                  DO k3 = -mx3,mx3
                     IF (k1 == 0 .AND. k2 == 0 .AND. k3 == 0) CYCLE
                     step(k1,k2,k3) = step(k1,k2,k3) + CONJG(row(k3))*pre*ph3(k3)
                  ENDDO
               ELSE
                  DO k3 = -mx3,mx3
                     IF (k1 == 0 .AND. k2 == 0 .AND. k3 == 0) CYCLE
                     step(k1,k2,k3) = step(k1,k2,k3) + row(k3)*pre*ph3(k3)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ELSE
         !general cell: every point has its own (k_z,k_perp)
         DO k3 = -mx3,mx3
            DO k2 = -mx2,mx2
               DO k1 = -mx1,mx1
                  IF (k1 == 0 .AND. k2 == 0 .AND. k3 == 0) CYCLE
                  rg = k1*bmat(1,:) + k2*bmat(2,:) + k3*bmat(3,:)
                  kp = SQRT(rg(1)**2+rg(2)**2)
                  IF (kp < 1.0E-10) THEN
                     CALL ballpartfou(rmt,z0,rg(3),0.0,v)
                  ELSE
                     v = 0.0
                     DO j = 1,nq
                        v = v + wq(j)*gf_bessel1(kp*sq(j))*EXP(CMPLX(0.0,rg(3)*xq(j)))
                     ENDDO
                     v = v*2.0*pi_const/kp
                  ENDIF
                  IF (zcoor >= 0) v = CONJG(v)
                  step(k1,k2,k3) = step(k1,k2,k3) + v*ph1(k1)*ph2(k2)*ph3(k3)/omtil
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      END SUBROUTINE ballpart_step

                                                                        
                                                                        
!two functions that will be integrated follow                           
! both need bessel-functions that need some compiler dependend call     
                                                                        
      FUNCTION integrandr(x) RESULT(integrand) 
                                                                        
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
