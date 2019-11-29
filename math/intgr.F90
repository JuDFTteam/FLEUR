MODULE m_intgr

  !**********************************************************************
  ! intgr[0-3]:
  !    Integrators of a function y(jri) on a logarithmic mesh with jri 
  !    mesh points. The output is either a scalar [z] or a field [z(jri)].
  !    Either the first meshpoint [r0,r1] and the increment [h] or the 
  !    array of meshpoints [rmsh(jri)] is supplied:
  !
  ! intgz0 & intgz1 :
  !     (in)definite integrators on linear mesh. tail=.true. will include
  !     a tail correction assuming that the function is a simple
  !     decaying exponential between the first mesh point and infinity.
  !     y contains the nmz function values tabulated at a spacing of h.
  !
  ! intgrt:
  !    Uses the basic trapezoidal rule of integration. Indefinite.
  !                                                 - A. Neukirchen, '19
  !
  !            integrator:      ---- input ----      output
  !    intgr0:   definite       y r0        h jri  | z
  !    intgr1: indefinite       y r1        h jri  | z(jri)
  !    intgr2: indefinite       y rmsh(jri) h jri  | z(jri)
  !    intgr3:   definite       y rmsh(jri) h jri  | z
  !    intgz0:   definite       y tail      n nmz  ! z
  !    intgz1: indefinite       y tail      n nmz  ! z(nmz)
  !    intgrt: indefinite       y x         jri    | z(jri)
  !    
  !                                                            m. weinert
  !**********************************************************************

#include"cpp_double.h"

  IMPLICIT NONE

  !INTRINSIC exp,log
  INTERFACE 
    REAL FUNCTION CPP_BLAS_sdot( n, f1, is1, f2, is2 )
      INTEGER, INTENT (IN) :: n, is1, is2
      REAL,    INTENT (IN) :: f1(n), f2(n)
    END FUNCTION 
  END INTERFACE

!  interface intgz1
!    module procedure intgz1Real, intgz1Complex
!  end interface

  interface intgz1Reverse
    module procedure intgz1RealReverse, intgz1ComplexReverse
  end interface
  

  INTEGER, PARAMETER, PRIVATE :: nr  = 7 , nr1 = 6
  REAL,    PARAMETER, PRIVATE :: h0 = 140. , zero = 0.0e0
  
  ! lagrangian integration coefficients (simpson 7 point rule: error  h**9)
  
  INTEGER, DIMENSION(7),   PARAMETER, PRIVATE :: ih = (/41,216,27,272,27,216,41/)

  REAL,    DIMENSION(7,5), PARAMETER, PRIVATE :: a = RESHAPE( &
       (/19087.,65112.,-46461., 37504.,-20211., 6312.,-863., &
          -863.,25128., 46989.,-16256.,  7299.,-2088., 271., &
           271.,-2760., 30819., 37504., -6771., 1608.,-191., &
          -191., 1608., -6771., 37504., 30819.,-2760., 271., &
           271.,-2088.,  7299.,-16256., 46989.,25128.,-863./),(/7,5/))


  CONTAINS

  !**********************************************************************
  SUBROUTINE intgr0( y, r0, h, jri, z )
    !**********************************************************************
    !     ..
    !     .. Arguments ..
    INTEGER,    INTENT (IN)  :: jri
    REAL,       INTENT (IN)  :: h, r0
    REAL,       INTENT (IN)  :: y(jri)
    REAL,       INTENT (OUT) :: z
    !     ..
    !     .. Locals ..
    INTEGER :: m, n0, nsteps
    REAL    :: dr, r(7)
    INTEGER :: i, j
    REAL    :: alpha, yr(7)

    !
    !--->    integral from 0 to r1 approximated by leading term in power
    !--->    series expansion of y(r)
    !
    z = zero
    IF (y(1)*y(2).GT.zero) THEN
      alpha = 1.0 + log(y(2)/y(1))/h
      IF (alpha.GT.zero) z = r0*y(1)/alpha
    ENDIF
    !
    !--->    determine steps and starting point for simpson
    !
    nsteps = (jri-1)/nr1
    n0 = jri - nr1*nsteps
    dr = exp(h)
    r(1) = r0
    DO i = 2,7
      r(i) = dr*r(i-1)
    ENDDO 
    !
    !--->    lagrange integration for points 1<j<n0, error: h**9
    !
    IF (n0.GT.1) THEN
      DO i = 1,7
        yr(i) = r(i)*y(i)
      ENDDO
      DO j = 1,n0 - 1
        z = z + h*CPP_BLAS_sdot(7,a(1,j),1,yr,1)/60480.
      ENDDO
    ENDIF
    r(1) = r(n0)
    !
    !--->    simpson integration
    !
    DO m = 1,nsteps
      DO i = 2,nr
        r(i) = dr*r(i-1)
      ENDDO
      DO i = 1,nr
        yr(i) = h*ih(i)*r(i)/h0
      ENDDO
      z = z + CPP_BLAS_sdot(nr,yr,1,y(n0),1)
      n0 = n0 + nr1
      r(1) = r(nr)
    ENDDO

    RETURN
  END SUBROUTINE intgr0

  !**********************************************************************
  SUBROUTINE intgr1( y, r1, h, jri, z )
    !**********************************************************************
    !     .. 
    !     .. Arguments ..
    INTEGER, INTENT (IN) :: jri
    REAL,    INTENT (IN) :: h,r1
    REAL,    INTENT (IN) :: y(jri)
    REAL,    INTENT (OUT):: z(jri)
    !     ..
    !     .. Locals ..
    REAL    :: dr, rr, r(7)
    INTEGER :: i, j
    REAL    :: alpha, yr(7)
    !
    !--->    integral from 0 to r1 approximated by leading term in power
    !--->    series expansion of y(r)
    !
    z(1) = zero
    IF (y(1)*y(2).GT.zero) THEN
      alpha = 1.0 + log(y(2)/y(1))/h
      IF (alpha.GT.zero) z(1) = r1*y(1)/alpha
    ENDIF
    !
    !--->    lagrange integration for points 1<j<nr, error: h**9
    !
    dr = exp(h)
    rr = r1
    DO i = 1,7
      yr(i) = rr*y(i)
      r(i) = rr
      rr = dr*rr
    ENDDO
    DO j = 1,nr - 2
      z(j+1) = z(j) + h*CPP_BLAS_sdot(7,a(1,j),1,yr,1)/60480.
    ENDDO
    !
    !--->    simpson integration, j>nr-1
    !
    DO i = 1,nr
      r(i) = h*ih(i)*r(i)/h0
    ENDDO
    DO j = nr,jri
      z(j) = z(j-nr1) + CPP_BLAS_sdot(nr,r,1,y(j-nr1),1)
      DO i = 1,7
        r(i) = dr*r(i)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE intgr1

  !**********************************************************************
  SUBROUTINE intgr2( y, rmsh, h, jri, z )
    !**********************************************************************
    !     ..
    !     .. Arguments ..
    INTEGER, INTENT (IN) :: jri
    REAL,    INTENT (IN) :: h
    REAL,    INTENT (IN) :: rmsh(jri), y(jri)
    REAL,    INTENT (OUT):: z(jri)
    !     ..
    !     .. Locals ..
    REAL    :: dr, r(7)
    INTEGER :: i, j
    REAL    :: alpha, yr(7)
    !
    !--->    integral from 0 to r1 approximated by leading term in power
    !--->    series expansion of y(r)
    !
    z(1) = zero
    IF (y(1)*y(2).GT.zero) THEN 
      alpha = 1.0 + log(y(2)/y(1))/h
      IF (alpha.GT.zero) z(1) = rmsh(1)*y(1)/alpha
    ENDIF
    !
    !--->    lagrange integration for points 1<j<nr, error: h**9
    !
    dr = exp(h)
    DO i = 1,7
      r(i) = rmsh(i)
      yr(i) = rmsh(i)*y(i)
    ENDDO
    DO j = 1,nr - 2
      z(j+1) = z(j) + h*CPP_BLAS_sdot(7,a(1,j),1,yr,1)/60480.
    ENDDO
    !
    !--->    simpson integration, j>nr-1
    !
    DO i = 1,nr
      r(i) = h*ih(i)*r(i)/h0
    ENDDO
    DO j = nr,jri
      z(j) = z(j-nr1) + CPP_BLAS_sdot(nr,r,1,y(j-nr1),1)
      DO i = 1,7
        r(i) = dr*r(i)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE intgr2

  !**********************************************************************
  SUBROUTINE intgr3( y, r, h, jri, z )
    !**********************************************************************
    !     ..
    !     .. Arguments ..
    INTEGER, INTENT (IN) :: jri
    REAL,    INTENT (IN) :: h
    REAL,    INTENT (IN) :: r(jri)
    REAL,    INTENT (IN) :: y(jri)
    REAL,    INTENT (OUT):: z
    !     ..
    !     .. Locals ..
    INTEGER :: m, n0, nsteps
    REAL    :: tiny, yr(nr), h1, z1, ih1(nr)
    INTEGER :: i, j
    REAL    :: alpha
    !
    !--->    integral from 0 to r1 approximated by leading term in power
    !--->    series expansion of y(r)
    !
    !      DO i=1,jri
    !        IF (abs(y(i)).LT.tiny) y(i) = tiny
    !      ENDDO
    !
    z = zero
    IF (y(1)*y(2).GT.zero) THEN
      alpha = 1.0 + log(y(2)/y(1))/h
      IF (alpha.GT.zero) z = r(1)*y(1)/alpha
    ENDIF
    !
    !--->    determine steps and starting point for simpson
    !
    nsteps = (jri-1)/nr1
    n0 = jri - nr1*nsteps
    !
    !--->    lagrange integration for points 1<j<n0, error: h**9
    !
    IF (n0.GT.1) THEN
      DO i = 1,7
        yr(i) = r(i)*y(i)
      ENDDO
      z1 = 0.
      DO j = 1,n0 - 1
        z1 = z1 + CPP_BLAS_sdot(7,a(1,j),1,yr,1)
      ENDDO
      z = z + z1 * h / 60480.
    END IF
    !
    !--->    simpson integration
    !
    h1 = h / h0
    DO i = 1,nr
      ih1(i) = h1 * ih(i)
    ENDDO
    DO m = 1,nsteps
      DO i = 1,nr
        yr(i) = ih1(i)*r(i+n0-1)
      ENDDO
      z = z + CPP_BLAS_sdot(nr,yr,1,y(n0),1)
      n0 = n0 + nr1
    ENDDO

    RETURN
  END SUBROUTINE intgr3

  !**********************************************************************
  SUBROUTINE intgz0( y, h, nmz, z, tail )
    !**********************************************************************
    !     ..
    !     .. Arguments ..
    INTEGER, INTENT (IN)  :: nmz
    REAL,    INTENT (IN)  :: h
    REAL,    INTENT (IN)  :: y(nmz)
    LOGICAL, INTENT (IN)  :: tail
    REAL,    INTENT (OUT) :: z
    !     ..
    !     .. Locals ..
    REAL    :: yl, ys
    INTEGER :: m, n0, nsteps
    INTEGER :: i, j
    REAL    :: alpha, yr(7)
    !
    !--->    integral from minus infinity to the first mesh point assuming
    !--->    exponential decay in this region. this contribution is limited
    !--->    to alpha>0.1 which corresponds to square of wavefunctions
    !--->    of energy > 0.00125 a.u. below the vacuum level
    !
    z = zero
    IF (tail) THEN
      IF (y(1)*y(2).GT.zero) THEN
        alpha = log(y(2)/y(1))/h
        IF (alpha.GT.0.1) z = y(1)/alpha
      ENDIF
    ENDIF
    !
    !--->    determine steps and starting point for simpson
    !
    nsteps = (nmz-1)/nr1
    n0 = nmz - nr1*nsteps
    !
    !--->    lagrange integration for points 1<j<n0, error: h**9
    !
    yl = zero
    IF (n0.GT.1) THEN
      DO j = 1, n0 - 1
        yl = yl + CPP_BLAS_sdot(7,a(1,j),1,y,1)
      ENDDO
      yl = h*yl/60480.
    END IF
    !
    !--->    simpson integration
    !
    ys = zero
    n0 = n0 - 1
    DO m = 1,nsteps
      DO i = 1,nr
        ys = ys + ih(i)*y(n0+i)
      ENDDO
      n0 = n0 + nr1
    ENDDO
    ys = h*ys/h0
    z = z + yl + ys

    RETURN
  END SUBROUTINE intgz0

  !**********************************************************************
  SUBROUTINE intgz1( y, h, nmz, z, tail )
    !**********************************************************************
    !     .. 
    !     .. Arguments ..
    INTEGER, INTENT (IN)  :: nmz
    LOGICAL, INTENT (IN)  :: tail
    REAL,    INTENT (IN)  :: h
    REAL,    INTENT (IN)  :: y(nmz)
    REAL,    INTENT (OUT) :: z(nmz)
    !     ..
    !     .. Locals ..
    REAL    :: yl, ys
    REAL, PARAMETER :: eps = 1.e-38
    INTEGER :: i, j
    REAL    :: alpha, yr(7)
    !
    !--->    integral from minus infinity to the first mesh point assuming
    !--->    exponential decay in this region. this contribution is limited
    !--->    to alpha>0.1 which corresponds to square of wavefunctions
    !--->    of energy > 0.00125 a.u. below the vacuum level
    !
    z(1) = zero
    IF (tail) THEN
      IF (abs(y(1)).GT.eps) THEN
        IF (y(1)*y(2).GT.zero) THEN
          alpha = log(y(2)/y(1))/h
          IF (alpha.GT.0.1) z(1) = y(1)/alpha
        ENDIF
      ENDIF
    ENDIF
    !
    !--->    lagrange integration for points 1<j<n0, error: h**9
    !
    DO j = 1,nr - 2
      yl = 0
      yl = yl + CPP_BLAS_sdot(7,a(1,j),1,y,1)
      z(j+1) = z(j) + h*yl/60480.
    ENDDO
    !
    !--->    simpson integration
    !
    DO j = nr,nmz
      ys = zero
      DO i = 1,nr
        ys = ys + ih(i)*y(j-nr+i)
      ENDDO
      z(j) = z(j-nr1) + h*ys/h0
    ENDDO

    RETURN
  END SUBROUTINE intgz1



  subroutine intgz1Complex( y, h, nmz, z, tail )

    integer, intent(in)  :: nmz
    logical, intent(in)  :: tail
    real,    intent(in)  :: h
    complex, intent(in)  :: y(nmz)
    complex, intent(out) :: z(nmz)

    real                 :: zr(nmz), zi(nmz)

    call intgz1(  real(y), h, nmz, zr, tail )
    call intgz1( aimag(y), h, nmz, zi, tail )
    z = cmplx( zr, zi )

  end subroutine intgz1Complex



  subroutine intgz1RealReverse( y, h, nmz, z, tail )

    integer, intent(in)  :: nmz
    logical, intent(in)  :: tail
    real,    intent(in)  :: h
    real,    intent(in)  :: y(nmz)
    real,    intent(out) :: z(nmz)

    real                 :: y_reverse(nmz), z_reverse(nmz)
    integer              :: i

    do i = 1, nmz
      y_reverse(i) = y(nmz+1-i)
    end do
    call intgz1( y_reverse, h, nmz, z_reverse, tail )
    do i = 1, nmz
      z(i) = z_reverse(nmz+1-i)
    end do    

  end subroutine intgz1RealReverse



  subroutine intgz1ComplexReverse( y, h, nmz, z, tail )

    integer, intent(in)  :: nmz
    logical, intent(in)  :: tail
    real,    intent(in)  :: h
    complex, intent(in)  :: y(nmz)
    complex, intent(out) :: z(nmz)

    complex              :: y_reverse(nmz), z_reverse(nmz)
    integer              :: i

    do i = 1, nmz
      y_reverse(i) = y(nmz+1-i)
    end do
    call intgz1Complex( y_reverse, h, nmz, z_reverse, tail )
    do i = 1, nmz
      z(i) = z_reverse(nmz+1-i)
    end do    

  end subroutine intgz1ComplexReverse

   SUBROUTINE intgrt(y,x,jri,z)
      INTEGER,    INTENT (IN)  :: jri
      REAL,       INTENT (IN)  :: x(jri), y(jri)
      REAL,       INTENT (OUT) :: z(jri)

      INTEGER                  :: i

      z=0.0

      DO i=2, jri
         z(i:)=z(i:)+(y(i-1)+y(i))*(x(i)-x(i-1))/2.0
      END DO

   END SUBROUTINE intgrt

END MODULE m_intgr
