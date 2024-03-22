!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_boxdim
CONTAINS
   SUBROUTINE boxdim( &
      bmat, &
      arltv1, arltv2, arltv3)
!*********************************************************************
!      This subroutine determines the maximum number L, M and N
!      nrgvl, nrgvm, nrgvn, of reciprocal lattice vectors G(L,M,N)
!      along the directions G(1), G(2), G(3), respectively, for which

!                     | G(L,M,N) | <= GMAX.

!      This equation defines a "sphere" and the nrgvl,m,n define
!      the DIMension of the BOX in which the sphere is placed.

!      In reality the G(i)'s, do not form a carteasian coordinate
!      system. Therefore the "sphere" is not a sphere, but an
!      ellipsoid. For this ellipsoid the largest L, M and N component
!      is determined as arltv1, arltv2, arltv3 to construct the boxes and

!                    nrgvl  = int( gmax/arltv1 ) + 1
!                    nrgvm  = int( gmax/arltv2 ) + 1
!                    nrgvn  = int( gmax/arltv3 ) + 1

!      G(i,xyz) is stored in bmat(i,xyz)

!      routine by s.bluegel from carpar-program

!                         S. Bl"ugel, IFF, 13. Nov. 97
!               tested by S. Heinze , IFF,
!*********************************************************************
      USE m_juDFT
! SE m_constants

!     .. Parameters ..
      IMPLICIT NONE

!     .. Scalar Arguments ..
      REAL, INTENT(OUT) :: arltv1, arltv2, arltv3
!     ..
!     .. Array Arguments ..
      REAL, INTENT(IN)  :: bmat(3, 3)
!     ..
!     .. Local Scalars ..
      INTEGER :: ixyz, j, k
      REAL ::    denom, eps, one, zero
!     ..
!     .. Local Arrays ..
      REAL :: det(3, 3), rr(3, 3)
!     ..
      DATA one, eps, zero/1.0, 1e-10, 0.0/

!--->  build up quadratic form for ellipsoid

      DO j = 1, 3
         DO k = 1, 3
            rr(k, j) = zero
            DO ixyz = 1, 3
               rr(k, j) = rr(k, j) + bmat(k, ixyz)*bmat(j, ixyz)
            ENDDO
         ENDDO
      ENDDO

!---> build determinants for Cramer's rule

      det(1, 1) = rr(2, 2)*rr(3, 3) - rr(3, 2)*rr(2, 3)
      det(1, 2) = rr(2, 1)*rr(3, 3) - rr(3, 1)*rr(2, 3)
      det(1, 3) = rr(2, 1)*rr(3, 2) - rr(3, 1)*rr(2, 2)
      det(2, 1) = rr(1, 2)*rr(3, 3) - rr(3, 2)*rr(1, 3)
      det(2, 2) = rr(1, 1)*rr(3, 3) - rr(3, 1)*rr(1, 3)
      det(2, 3) = rr(1, 1)*rr(3, 2) - rr(3, 1)*rr(1, 2)
      det(3, 1) = rr(1, 2)*rr(2, 3) - rr(2, 2)*rr(1, 3)
      det(3, 2) = rr(1, 1)*rr(2, 3) - rr(2, 1)*rr(1, 3)
      det(3, 3) = rr(1, 1)*rr(2, 2) - rr(2, 1)*rr(1, 2)

!---> check on the zeros of some determinants

      DO j = 1, 3
         IF (det(j, j) < eps) THEN
            CALL juDFT_error("boxdim: determinant. Problem with det(" &
                             //int2str(j) // ","//int2str(j) // ")", calledby="boxdim")
         END IF
      ENDDO

!---> scale determinants

      DO k = 1, 3
         denom = one/det(k, k)
         DO j = 1, 3
            det(k, j) = denom*det(k, j)
         ENDDO
         det(k, k) = one/denom
      ENDDO

!---> calculate the maximum l, m, n components of the ellipsoid

      arltv1 = sqrt(rr(1, 1) + rr(2, 2)*det(1, 2)*det(1, 2) &
                    + rr(3, 3)*det(1, 3)*det(1, 3) &
                    - (rr(1, 2) + rr(1, 2))*det(1, 2) &
                    + (rr(1, 3) + rr(1, 3))*det(1, 3) &
                    - (rr(2, 3) + rr(2, 3))*det(1, 2)*det(1, 3))
      arltv2 = sqrt(rr(2, 2) + rr(1, 1)*det(2, 1)*det(2, 1) &
                    + rr(3, 3)*det(2, 3)*det(2, 3) &
                    - (rr(1, 2) + rr(1, 2))*det(2, 1) &
                    + (rr(1, 3) + rr(1, 3))*det(2, 1)*det(2, 3) &
                    - (rr(2, 3) + rr(2, 3))*det(2, 3))
      arltv3 = sqrt(rr(3, 3) + rr(1, 1)*det(3, 1)*det(3, 1) &
                    + rr(2, 2)*det(3, 2)*det(3, 2) &
                    - (rr(1, 2) + rr(1, 2))*det(3, 1)*det(3, 2) &
                    + (rr(1, 3) + rr(1, 3))*det(3, 1) &
                    - (rr(2, 3) + rr(2, 3))*det(3, 2))

      RETURN
   end SUBROUTINE boxdim
end module m_boxdim
