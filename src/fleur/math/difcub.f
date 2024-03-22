!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_difcub

      CONTAINS

      REAL FUNCTION difcub(x,f,xi)
c     **********************************************************
c     differentiate the function f, given at the
c     points x0,x1,x2,x3 at the point xi by lagrange
c     interpolation for polynomial of 3rd order
c     r.p.
c     ***********************************************************
      IMPLICIT NONE
C     .. Scalar Arguments ..
      REAL xi
C     ..
C     .. Array Arguments ..
      REAL f(0:3),x(0:3)
C     ..
      difcub = ((xi-x(1))* (xi-x(2))+ (xi-x(1))* (xi-x(3))+
     +         (xi-x(2))* (xi-x(3)))*f(0)/ ((x(0)-x(1))* (x(0)-x(2))*
     +         (x(0)-x(3))) + ((xi-x(0))* (xi-x(2))+
     +         (xi-x(0))* (xi-x(3))+ (xi-x(2))* (xi-x(3)))*f(1)/
     +         ((x(1)-x(0))* (x(1)-x(2))* (x(1)-x(3))) +
     +         ((xi-x(0))* (xi-x(1))+ (xi-x(0))* (xi-x(3))+
     +         (xi-x(1))* (xi-x(3)))*f(2)/ ((x(2)-x(0))* (x(2)-x(1))*
     +         (x(2)-x(3))) + ((xi-x(0))* (xi-x(1))+
     +         (xi-x(0))* (xi-x(2))+ (xi-x(1))* (xi-x(2)))*f(3)/
     +         ((x(3)-x(0))* (x(3)-x(1))* (x(3)-x(2)))
      RETURN
      END FUNCTION difcub

      END MODULE m_difcub
