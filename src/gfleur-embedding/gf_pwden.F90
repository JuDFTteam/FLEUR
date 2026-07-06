!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_pwden 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_pwden(                                              &
     &     lapw,stars,jspin,omtil,G,                                    &
     &     qpw)                                                         
!*********************************************************************  
!     Calculates the interstitial charge density from the               
!     Green FUNCTION. A DOUBLE sum over the G-vectors must be performed 
!                                                                       
!     Daniel Wortmann                                                   
!*********************************************************************  
      USE m_gf_types 
      USE m_constants,ONLY:pimach 
      IMPLICIT NONE 
!                                                                       
!     .. Scalar Arguments ..                                            
      TYPE(t_lapw),INTENT(IN) :: lapw 
      INTEGER, INTENT (IN)    :: jspin 
      REAL,    INTENT (IN)    :: omtil 
!     ..                                                                
!     .. Array Arguments ..                                             
      TYPE(t_stars),INTENT(IN):: stars 
      COMPLEX, INTENT (IN)    :: G(:,:) 
                                                                        
      COMPLEX, INTENT (INOUT) :: qpw(:)
!     Locals                                                            
      REAL    :: pi 
      INTEGER :: n1,n2,i1,i2,i3,in 
      COMPLEX :: dg,phase 
      COMPLEX,ALLOCATABLE::rg(:,:,:) 
                                                                        
      pi=pimach() 
      ALLOCATE(rg(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,            &
     &     -stars%mx3:stars%mx3))                                       
                                                                        
      rg=CMPLX(0.0,0.0) 
                                ! Calculate the reduced G-function      
      DO n1=1,lapw%nv_sphere(jspin)
!         IF (ABS(lapw%k%k3(n1,jspin))>stars%mx3) CYCLE                 
         DO n2=1,lapw%nv_sphere(jspin)
!            IF (ABS(lapw%k%k3(n2,jspin))>stars%mx3) CYCLE              
            i1 = lapw%k%k1(n2,jspin) - lapw%k%k1(n1,jspin) 
            IF (iabs(i1)>stars%mx1) CYCLE 
            i2 = lapw%k%k2(n2,jspin) - lapw%k%k2(n1,jspin) 
            IF (iabs(i2)>stars%mx2) CYCLE 
            i3 = lapw%k%k3(n2,jspin) - lapw%k%k3(n1,jspin) 
            IF (iabs(i3)>stars%mx3) CYCLE 
            rg(i1,i2,i3)=rg(i1,i2,i3)+g(n2,n1) 
         ENDDO 
      ENDDO 
                                ! Now calculate the density for all star
      DO i1=-stars%mx1,stars%mx1 
         DO i2=-stars%mx2,stars%mx2 
            DO i3=-stars%mx3,stars%mx3 
               in=stars%ig(i1,i2,i3) 
               IF (in==0) CYCLE 
               phase = stars%rgphs(i1,i2,i3)/ (stars%nstr(in)*omtil)/pi 
                                                            !factor 1/Pi
               dg=CMPLX(0.0,0.5)*(rg(i1,i2,i3)-CONJG(rg(-i1,-i2,-i3))) 
               qpw(in)=qpw(in)+                             &
     &              phase*dg                                            
               ! At this point one could also calculate the coef for -i1
               ! might be faster??? But doublecounting has to be avoided
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
                                                                        
      RETURN 
      END SUBROUTINE 
      END                                           
