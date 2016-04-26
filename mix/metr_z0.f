      MODULE m_metrz0
      use m_juDFT
c     *****************************************************
c     calculates weights for a 7-point simpson integration
c     in analogy to intgz0
c                         r.pentcheva,24.06.96,kfa
c     work array whelp eliminated
c                         s.bluegel,  BNL 9.Aug.96
c     *****************************************************
      CONTAINS
      SUBROUTINE metr_z0(
     >                   nz,
     <                   wght)

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nz
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (OUT) :: wght(nz)
C     ..
C     .. Local Scalars ..
      INTEGER i,iz,j,iz0,jz,nsteps
      INTEGER, PARAMETER :: nz7  = 7 , nz6 = 6
      REAL,    PARAMETER :: h0 = 140. 
C     ..
C     .. Local Arrays ..
!
! lagrangian integration coefficients (simpson 7 point rule: error  h**9)
!
      INTEGER, DIMENSION(7),   PARAMETER :: ih =
     +                                       (/41,216,27,272,27,216,41/)
      REAL,    DIMENSION(7,5), PARAMETER :: a = RESHAPE(
     +  (/19087.,65112.,-46461., 37504.,-20211., 6312.,-863.,
     +     -863.,25128., 46989.,-16256.,  7299.,-2088., 271.,
     +      271.,-2760., 30819., 37504., -6771., 1608.,-191.,
     +     -191., 1608., -6771., 37504., 30819.,-2760., 271.,
     +      271.,-2088.,  7299.,-16256., 46989.,25128.,-863./),(/7,5/))
C     ..
      wght(:)=0.0
      nsteps = (nz-1)/nz6
      iz0 = nz-nz6*nsteps
      DO iz = 1, iz0-1
c
c---> for iz-points, 1<iz<iz0<7, use lagrange interpolation, error: h**9
         DO jz = 1, nz7
            wght(nz+1-iz) = wght(nz+1-iz)+ a(jz,iz)/432.0
         ENDDO
      ENDDO
c---> weights for simpson integration
      iz0 = nz - iz0 + 2
      DO j=1,nsteps
         DO i=1,nz7
            wght(iz0-i) = wght(iz0-i) + ih(i)
         ENDDO
         iz0 = iz0 - nz6
      ENDDO
      IF (iz0-1.ne.1)  CALL juDFT_error("iz0><nz0",calledby="metr_z0")
      wght(:)=wght(:)/h0

      END SUBROUTINE metr_z0
      END MODULE m_metrz0

