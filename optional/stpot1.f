      MODULE m_stpot1
      CONTAINS
      SUBROUTINE stpot1(
     >                  msh,n,z,rad,
     <                  vr1)
c     **********************************************************
c     create a starting potential
c     **********************************************************

      IMPLICIT NONE
C     ..
      INTEGER,INTENT (IN) :: n,msh
      REAL,INTENT (IN)    :: z,rad(msh)
      REAL,INTENT (OUT)   :: vr1(msh)
C     ..
C     .. Local Scalars ..
      INTEGER i,i1
      REAL    d1,h1
C     ..
      IF (z.LT.1.0) THEN
        DO i = 1,n
           vr1(i) = 1.e-3
        ENDDO
        RETURN
      ENDIF
C
      d1 = 0.9e0
      h1 = d1* (z-1)**0.4e0
      i1 = 0.75e0*n
      DO i = 1,n
         vr1(i) = -1.0
      ENDDO
      DO i = 1,i1
         vr1(i) = vr1(i)* ((z-1)/ (h1* (exp(rad(i)/d1)-1)+1)+1)
      ENDDO
      RETURN
      END SUBROUTINE
      END
