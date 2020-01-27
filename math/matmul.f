      MODULE m_matmul
!
! makes a matrix-matrix multiplication for 3x3 matrices A and B
! whereby      A       B
!            integer  integer  ... matmul1
!            real     real     ... matmul2
!            integer  real     ... matmul3
! additionally, integer rotation matrices + nonsymmorphic 
! translations can be multiplied in matmul4
!
      CONTAINS
!--------------------------------------------------------
      SUBROUTINE matmul1 (ma,mb,mc)

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: ma(3,3)
      INTEGER, INTENT (IN) :: mb(3,3)
      INTEGER, INTENT (OUT):: mc(3,3)

      INTEGER :: i,j,k

      INTEGER x
      DO i=1,3
         DO k=1,3
         x = 0
            DO j=1,3
               x = x + ma(i,j)*mb(j,k)
            ENDDO
         mc(i,k) = x
         END DO
      END DO

      END SUBROUTINE matmul1
!--------------------------------------------------------
      SUBROUTINE matmul2 (a,b,c)

      IMPLICIT NONE
      REAL,    INTENT (IN) :: a(3,3)
      REAL,    INTENT (IN) :: b(3,3)
      REAL,    INTENT (OUT):: c(3,3)

      INTEGER :: i,j,k
      REAL x
      DO i=1,3
         DO k=1,3
         x=0.e0
            DO j=1,3
               x = x + a(i,j)*b(j,k)
            ENDDO
         c(i,k) = x
         END DO
      END DO

      END SUBROUTINE matmul2
!--------------------------------------------------------
      SUBROUTINE matmul3 (ma,b,c)

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: ma(3,3)
      REAL,    INTENT (IN) ::  b(3,3)
      REAL,    INTENT (OUT)::  c(3,3)

      INTEGER :: i,j,k
      REAL x
      DO i=1,3
         DO k=1,3
         x=0.e0
            DO j=1,3
               x = x + ma(i,j)*b(j,k)
            ENDDO
         c(i,k) = x
         END DO
      END DO

      END SUBROUTINE matmul3
!--------------------------------------------------------
      SUBROUTINE matmul3r (ma,b,c)

      REAL, INTENT (IN) :: ma(3,3)
      REAL, INTENT (IN) ::  b(3,3)
      REAL, INTENT (OUT)::  c(3,3)

      INTEGER :: i,j,k
      REAL x
      DO i=1,3
         DO k=1,3
         x=0.e0
            DO j=1,3
               x = x + ma(i,j)*b(j,k)
            ENDDO
         c(i,k) = x
         END DO
      END DO
      RETURN
      END SUBROUTINE matmul3r
!--------------------------------------------------------
      SUBROUTINE matmul4 (ma,ta,mb,tb,mc,tc)

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: ma(3,3),mb(3,3)
      REAL,    INTENT (IN) :: ta(3),tb(3)
      INTEGER, INTENT (OUT):: mc(3,3)
      REAL,    INTENT (OUT):: tc(3)

      INTEGER :: i,j,k
      REAL x,xa(4,4),xb(4,4),xc(4,4)

      xa(:,:) = 0.0 ; xa(4,4) = 1.0 ; xb(:,:) = xa(:,:)
      xa(1:3,1:3) = real(ma(1:3,1:3)) ; xa(1:3,4) = ta(:) 
      xb(1:3,1:3) = real(mb(1:3,1:3)) ; xb(1:3,4) = tb(:)
      DO i=1,4
         DO k=1,4
         x=0.e0
            DO j=1,4
               x = x + xa(i,j)*xb(j,k)
            ENDDO
         xc(i,k) = x
         END DO
      END DO
      mc(1:3,1:3) = nint(xc(1:3,1:3)) ; tc(:) = xc(1:3,4) 
      DO i = 1,3
        IF (tc(i).GT.1.0) THEN
          tc(i) = tc(i) - int(tc(i))
        ELSEIF (tc(i).LT.0.0) THEN
          tc(i) = tc(i) + int(tc(i)) + 1
        ENDIF
      ENDDO

      END SUBROUTINE matmul4

      END MODULE m_matmul
