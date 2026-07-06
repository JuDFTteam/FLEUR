!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_bessel 
      use m_juDFT
          IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Calculate bessel functions                                       
!                 Daniel Wortmann, (08-02-19)                           
!-----------------------------------------------                        
      CONTAINS 
                                                                        
      !<-- F: gf_bessel0(x)                                             

      FUNCTION gf_bessel0(x) 
!-----------------------------------------------                        
!Cylinder Bessel function of order zero                                 
!             (last modified: 08-04-15) D. Wortmann                     
!-----------------------------------------------                        
#ifdef INTEL_COMPILER                                                   
      USE ifport 
      IMPLICIT NONE 
      REAL   ,INTENT(IN)     ::  x 
      REAL                   :: gf_bessel0 
      gf_bessel0 = dbesj0(x) 
#else                                                                   
            IMPLICIT NONE 
!     ..Arguments ..                                                    
      REAL,    INTENT  (IN) :: x 
      REAL gf_bessel0 
                                                                        
!     .. Parameters ..                                                  
      REAL,    PARAMETER :: zero = 0.0 
      INTEGER, PARAMETER :: m1 = 0 
!     ..Locals ..                                                       
      INTEGER :: m,i,mass 
      REAL :: quot 
      REAL, ALLOCATABLE :: aux(:) 
!     ..                                                                
                                                                        
                                                                        
      IF (x<zero)  CALL juDFT_error("gf_bessel0",calledby="gf_bessel.F90")
                                                                        
      IF (x==zero .AND. m1==0) THEN 
         gf_bessel0 = 1. 
         RETURN 
      END IF 
      IF (x==zero .AND. m1/=0) THEN 
         gf_bessel0 = 0. 
         RETURN 
      END IF 
                                                                        
      mass = INT( m1 + 50 + x ) 
      ALLOCATE ( aux(0:mass) ) 
      aux(mass) = 0.0 
      aux(mass-1) = 1.0E-22 
                                                                        
      DO i=mass-2,0,-1 
        aux(i) = 2*(i+1)*aux(i+1)/x - aux(i+2) 
      END DO 
                                                                        
      quot = aux(0) 
                                                                        
      DO i=1,INT( mass/2. ) 
        quot = quot + 2*aux(2*i) 
      END DO 
                                                                        
      gf_bessel0 = aux(m1)/quot 
                                                                        
      IF (x==zero .AND. m1==0) THEN 
         gf_bessel0 = 1. 
         RETURN 
      END IF 
      IF (x==zero .AND. m1/=0) THEN 
         gf_bessel0 = 0. 
         RETURN 
      END IF 
                                                                        
      DEALLOCATE ( aux ) 
                                                                        
      RETURN 
#endif                                                                  
                                                                        
      END FUNCTION 

      !>                                                                
                                                                        
      !<-- F: gf_bessel1(x)                                             
      FUNCTION gf_bessel1(x) 
!-----------------------------------------------                        
!Cylinder Bessel function of first order                                
!             (last modified: 08-02-19) D. Wortmann                     
!-----------------------------------------------                        
#ifdef INTEL_COMPILER                                                   
      USE ifport 
      IMPLICIT NONE 
      REAL   ,INTENT(IN)     ::  x 
      REAL                   :: gf_bessel1 
      gf_bessel1 = dbesj0(x) 
#else                                                                   
      IMPLICIT NONE 
!     ..                                                                
!     ..Arguments ..                                                    
                                                                        
      REAL,    INTENT  (IN) :: x 
      REAL gf_bessel1 
!                                                                       
!     .. Parameters ..                                                  
      REAL,    PARAMETER :: zero = 0.0 
      INTEGER,PARAMETER  :: m1=1 
!     ..Locals ..                                                       
      INTEGER :: i,mass 
      REAL :: quot 
      REAL, ALLOCATABLE :: aux(:) 
!     ..                                                                
                                                                        
                                                                        
      IF (x<zero)  CALL juDFT_error("cylbes2",calledby="gf_bessel.F90")
                                                                        
      IF (x==zero .AND. m1==0) THEN 
         gf_bessel1 = 1. 
         RETURN 
      END IF 
      IF (x==zero .AND. m1/=0) THEN 
         gf_bessel1 = 0. 
         RETURN 
      END IF 
                                                                        
      mass = INT( m1 + 50 + x ) 
      ALLOCATE ( aux(0:mass) ) 
      aux(mass) = 0.0 
      aux(mass-1) = 1.0E-22 
                                                                        
      DO i=mass-2,0,-1 
        aux(i) = 2*(i+1)*aux(i+1)/x - aux(i+2) 
      END DO 
                                                                        
      quot = aux(0) 
                                                                        
      DO i=1,INT( mass/2. ) 
        quot = quot + 2*aux(2*i) 
      END DO 
                                                                        
      gf_bessel1 = aux(m1)/quot 
                                                                        
      DEALLOCATE ( aux ) 
#endif                                                                  
      END FUNCTION 
      !>                                                                
                                                                        
      END                                           
