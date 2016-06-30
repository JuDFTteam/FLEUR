!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_od_cylbes
      use m_juDFT
      CONTAINS
      SUBROUTINE od_cylbes(
     >     m1,x,
     <     fJ) 

      IMPLICIT NONE
!     ..
!     ..Arguments ..
      INTEGER, INTENT  (IN) :: m1
      REAL,    INTENT  (IN) :: x
      REAL,    INTENT (OUT) :: fJ
!
!     .. Parameters ..
      REAL,    PARAMETER :: zero = 0.0
!     ..Locals ..
      INTEGER :: m,i,mass
      REAL :: quot
      REAL, ALLOCATABLE :: aux(:)
!     ..

      IF (x.LT.zero)  CALL juDFT_error("cylbes2",calledby="od_cylbes")

      IF (x.EQ.zero .AND. m1.EQ.0) THEN
         fJ = 1.
         RETURN
      END IF
      IF (x.EQ.zero .AND. m1.NE.0) THEN
         fJ = 0.
         RETURN
      END IF

      if (m1.lt.0) then
         m = -m1
      else
         m = m1
      end if
        
      mass = INT( m + 50 + x )
      ALLOCATE ( aux(0:mass) )       
      aux(mass) = 0.0
      aux(mass-1) = 1.0e-22   
     
      DO i=mass-2,0,-1
        aux(i) = 2*(i+1)*aux(i+1)/x - aux(i+2)
      END DO

      quot = aux(0)

      DO i=1,INT( mass/2. )
        quot = quot + 2*aux(2*i)
      END DO  

      IF (m1.LT.0) THEN
         fJ = ((-1)**m)*aux(m)/quot
      ELSE
         fJ = aux(m)/quot
      END IF
     
      DEALLOCATE ( aux )
      
      RETURN
      END SUBROUTINE od_cylbes
      END MODULE m_od_cylbes
    
 
