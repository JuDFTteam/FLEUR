!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_differentiate
CONTAINS
  REAL FUNCTION difcub(x,f,xi)
    !     **********************************************************
    !     differentiate the function f, given at the
    !     points x0,x1,x2,x3 at the point xi by lagrange
    !     interpolation for polynomial of 3rd order
    !     r.p.
    !     ***********************************************************
    IMPLICIT NONE
    !     .. Scalar Arguments ..
    REAL,INTENT(IN):: xi
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN):: f(0:3),x(0:3)
    !     ..
    difcub = ((xi-x(1))* (xi-x(2))+ (xi-x(1))* (xi-x(3))+&
         (xi-x(2))* (xi-x(3)))*f(0)/ ((x(0)-x(1))* (x(0)-x(2))*&
         (x(0)-x(3))) + ((xi-x(0))* (xi-x(2))+&
         (xi-x(0))* (xi-x(3))+ (xi-x(2))* (xi-x(3)))*f(1)/&
         ((x(1)-x(0))* (x(1)-x(2))* (x(1)-x(3))) +&
         ((xi-x(0))* (xi-x(1))+ (xi-x(0))* (xi-x(3))+&
         (xi-x(1))* (xi-x(3)))*f(2)/ ((x(2)-x(0))* (x(2)-x(1))*&
         (x(2)-x(3))) + ((xi-x(0))* (xi-x(1))+&
         (xi-x(0))* (xi-x(2))+ (xi-x(1))* (xi-x(2)))*f(3)/&
         ((x(3)-x(0))* (x(3)-x(1))* (x(3)-x(2)))
    RETURN
  END FUNCTION difcub
  SUBROUTINE diff3(&
       f,dx,&
       df)
    !********************************************************************
    !     differetiation via 3-points
    !********************************************************************

    IMPLICIT NONE

    !     .. Scalar Arguments ..
    REAL,    INTENT (IN) :: dx
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT (IN)  ::  f(:)
    REAL, INTENT (OUT) :: df(:)
    !     ..
    !     .. Local Scalars ..
    INTEGER i,jri
    REAL tdx_i
    !     ..
    jri=size(f)
    tdx_i = 1./(2.*dx)
    !
    !---> first point
    df(1) = -tdx_i * (-3.*f(1)+4.*f(2)-f(3))
    !
    !---> central point formula in charge
    DO i = 2,jri - 1
       df(i) = tdx_i * (f(i+1)-f(i-1))
    END DO
    !
    !---> last point
    df(jri) = tdx_i * (3.*f(jri)-4.*f(jri-1)+f(jri-2))
    !
    RETURN
  END SUBROUTINE diff3

END MODULE m_differentiate
