      MODULE m_convndim
      use m_juDFT
      CONTAINS
      SUBROUTINE convn_dim(
     >                     gmaxr,
     <                     ncvd)
c     ***********************************************************
c     determines the optimum values for the convergence parameter
c     for each atom type using the criterion discussed in
c     m. weinert, j. math. phys. 22, 2433 (1981).  each sphere
c     and l component may have different values.  (psqpw changed
c     to allow this option).
c          m. weinert july 1982
c     ***********************************************************
      IMPLICIT NONE

      REAL gmaxr
      INTEGER ncvd

      REAL z0,z(17)
      INTEGER i,n1

c     .. data statements ..
      DATA z/6.9e0,8.1e0,9.3e0,10.5e0,11.6e0,12.7e0,13.9e0,15.0e0,
     +     16.1e0,17.2e0,18.3e0,19.4e0,20.5e0,21.6e0,22.7e0,23.7e0,
     +     24.8e0/,z0/5.7e0/
c     ..

      IF (gmaxr.LT.z0) THEN
         WRITE (6,'('' gmax.r too small:'',f10.5)') gmaxr
          CALL juDFT_error("convn",calledby="convn_dim")
      END IF

      IF (gmaxr.GT.z(17)) THEN
         n1 = 0.9e0* (gmaxr-z(17))
         ncvd = 18 + n1
      ELSE
         DO i = 1,17
            IF (gmaxr.LE.z(i)) THEN
               ncvd = i
            END IF
         END DO
      END IF

      END SUBROUTINE convn_dim
      END MODULE m_convndim
