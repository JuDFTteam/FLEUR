      MODULE m_angle
      CONTAINS
      REAL ELEMENTAL FUNCTION angle(x,y)
      
c----------------------------------------------
c     calculates an angle of a vector
c     given by rectangular coordinates (x,y)
c-----------------------------------------------

      USE m_constants
      IMPLICIT NONE

      REAL, INTENT (IN) :: x,y

      IF (x.NE.0.0) THEN
         IF (y.NE.0.0) THEN
            IF (y.GT.0.0) THEN
               IF (x.LT.0.0) angle = pi_const - ATAN(y/abs(x))
               IF (x.GT.0.0) angle = ATAN(y/x)
            ELSE
               IF (x.GT.0.0) angle = - ATAN(abs(y)/x)
               IF (x.LT.0.0) angle = - pi_const 
     +                                  + ATAN(abs(y)/abs(x))
            END IF
         ELSE
            IF (x.LT.0.0) angle = -pi_const
            IF (x.GT.0.0) angle = 0.0
         END IF
      ELSE
         IF (y.LT.0.0) angle = -pi_const/2.
         IF (y.EQ.0.0) angle = 0.0 
         IF (y.GT.0.0) angle = pi_const/2.
      END IF    
      
      END FUNCTION angle
      END MODULE m_angle
