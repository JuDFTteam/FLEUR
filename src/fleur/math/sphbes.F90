!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_sphbes
   use m_juDFT
   IMPLICIT NONE

   private

   !Module contains three subroutines:
   ! sphbes: calculates spherical Bessel functions
   ! dsphbs: calculates the derivative of spherical Bessel functions
   ! d_sphbes: calculates spherical Bessel functions and their derivatives for a vector of arguments
   public :: sphbes, dsphbs, d_sphbes

CONTAINS
   SUBROUTINE sphbes(lmax, x, fj)
!c********************************************************************
!c calculate spherical Bessel functions and derivatives of sqrt(e)*r
!c                                   P. Marksteiner and E. Badralexe
!c********************************************************************
!     ..
!     .. Arguments ..
      INTEGER, INTENT(IN) :: lmax
      REAL, INTENT(IN) :: x
      REAL, INTENT(OUT) :: fj(0:lmax)
!
!     .. Parameters ..
      REAL, PARAMETER :: small = 1.0e-03, zero = 0.0
!     ..
!     .. Locals ..
      INTEGER i, l, min, n
      REAL fac, quot, xinv, xx
      REAL :: aux(0:int(lmax + 10 + x))
      !REAL, ALLOCATABLE :: aux(:)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs, cos, sin

      IF (lmax == 0) THEN
         IF (x .GE. small) THEN
            fj(0) = sin(x)/x
         ELSE
            fj(0) = 1 - x*x/6.*(1.-x*x/20.*(1.-x*x/42.))
         END IF
         RETURN
      END IF

      xx = x*x
      IF (x .GE. small) THEN
         xinv = 1./x
         fj(0) = sin(x)*xinv
         fj(1) = (fj(0) - cos(x))*xinv
      ELSE
         fj(0) = 1 - xx/6.*(1.-xx/20.*(1.-xx/42.))
         fj(1) = (x/3.)*(1.-xx/10.*(1.-xx/28.))
         fac = xx/15.
         DO l = 2, lmax
            fj(l) = fac*(1.-xx/(4*l + 6))
            fac = x*fac/(2*l + 3)
         END DO
         RETURN
      END IF
      IF (lmax .LT. x) THEN

         DO l = 2, lmax
            fj(l) = (2*l - 1)*xinv*fj(l - 1) - fj(l - 2)
         END DO

      ELSE IF (lmax .GE. 2) THEN
         n = INT(lmax + 10 + x)
         !  ALLOCATE( aux(0:n) )
!
! downward recursion from arbitrary starting values
!
         aux(n) = 0.
         aux(n - 1) = 1.
         DO i = n - 1, 1, -1
            aux(i - 1) = (2*i + 1)*xinv*aux(i) - aux(i + 1)
         END DO
!
! normalize with j0 or j1, whichever is larger
!
         min = 0
         IF (abs(fj(0)) .LT. abs(fj(1))) min = 1
         quot = fj(min)/aux(min)
         DO l = 2, lmax
            fj(l) = aux(l)*quot
         END DO
         ! DEALLOCATE( aux )
      END IF

      RETURN
   END SUBROUTINE sphbes

   SUBROUTINE dsphbs(lmax, x, fj, dfj)
!********************************************************************
!     calculates the derivative of the spherical bessel functions
!     dfj(l) = d jl(x)/dx
!     for l=0,lmax and argument x
!     note that the spherical bessel functions fj(l), l=0,lmax are
!     needed (call sphbes to generate these)
!                                            m. weinert
!********************************************************************
!     ..
!     .. Arguments ..
      INTEGER, INTENT(IN) :: lmax
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) ::  fj(0:lmax)
      REAL, INTENT(OUT) :: dfj(0:lmax)
!
!     .. Parameters ..
      REAL, PARAMETER :: xlim = 1.0e-04
!     ..
!     .. Local Scalars ..
      REAL fac, x2
      INTEGER l

      dfj(0) = -fj(1)
!--->    small x limit
      IF (x .LT. xlim) THEN
         x2 = 0.5*x*x
         fac = 1./3.
         DO l = 1, lmax
            dfj(l) = fac*(l - x2*(l + 2)/(2*l + 3))
            fac = x*fac/(2*l + 3)
         END DO
      ELSE
!--->    obtain dfj using recurrence relationship
         DO l = 1, lmax
            dfj(l) = fj(l - 1) - (l + 1)*fj(l)/x
         END DO
      END IF

      RETURN
   END SUBROUTINE dsphbs

   SUBROUTINE d_sphbes(lmax, xvec, fj, dfj)
      USE m_juDFT
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN) :: lmax             ! Maximum order of spherical Bessel functions
      REAL, INTENT(IN) :: xvec(:)             ! Vector of arguments for the spherical Bessel functions
      REAL, Allocatable, INTENT(OUT) :: fj(:, :)      ! 2D array to store spherical Bessel functions for all orders and arguments
      REAL, Allocatable, INTENT(OUT) :: dfj(:, :)     ! 2D array to store derivatives of spherical Bessel functions

      ! Parameters
      REAL, PARAMETER :: small = 1.0e-03, xlim = 1.0e-04

      ! Local variables
      INTEGER :: i, l, min, n, nx
      REAL :: fac, quot, x, x2
      REAL :: aux(0:INT(lmax + 10 + MAXVAL(xvec))) ! Temporary array for downward recursion

      if (allocated(fj)) deallocate (fj)
      if (allocated(dfj)) deallocate (dfj)
      allocate (fj(0:lmax, SIZE(xvec)))
      allocate (dfj(0:lmax, SIZE(xvec)))
      !$acc enter data create(fj, dfj)
      ! Get the size of the input vector
      nx = SIZE(xvec)

      !CPP_OMP PARALLEL DO default(none) shared(lmax, xvec, fj, dfj, small, xlim)&
      !CPP_OMP & private(i, l, min, n, aux, fac, quot, x)
      ! Loop over each value in the input vector
      !$acc kernels copyin(lmax, xvec, small, xlim) &
      !$acc create(aux)present(fj,dfj) &
      !$acc loop independent private(i, l, min, n, fac, quot, x)
      DO i = 1, nx
         x = xvec(i)

         ! Handle special case for lmax = 0
         IF (lmax == 0) THEN
            IF (x >= small) THEN
               fj(0, i) = SIN(x)/x
            ELSE
               fj(0, i) = 1.0 - x*x/6.0*(1.0 - x*x/20.0*(1.0 - x*x/42.0))
            END IF
            dfj(0, i) = -fj(1, i)
            CYCLE
         END IF

         ! Case for x >= small
         IF (x >= small) THEN
            fj(0, i) = SIN(x)/x
            fj(1, i) = (fj(0, i) - COS(x))/x

            ! Case for x < small (series expansion)
         ELSE
            fj(0, i) = 1.0 - x*x/6.0*(1.0 - x*x/20.0*(1.0 - x*x/42.0))
            fj(1, i) = (x/3.0)*(1.0 - x*x/10.0*(1.0 - x*x/28.0))
            fac = x*x/15.0

            DO l = 2, lmax
               fj(l, i) = fac*(1.0 - x*x/(4*l + 6))
               fac = x*fac/(2*l + 3)
            END DO
            CYCLE
         END IF

         if (lmax < x) then
            ! Use recurrence relation for l >= 2
            DO l = 2, lmax
               fj(l, i) = (2*l - 1)/x*fj(l - 1, i) - fj(l - 2, i)
            END DO
         ELSE
            ! Handle case where lmax >= x using downward recursion
            IF (lmax >= 2) THEN
               n = INT(lmax + 10 + x)

               ! Initialize downward recursion
               aux(n) = 0.0
               aux(n - 1) = 1.0

               DO min = n - 1, 1, -1
                  aux(min - 1) = (2*min + 1)/x*aux(min) - aux(min + 1)
               END DO

               ! Normalize with j0 or j1, whichever is larger
               min = 0
               IF (ABS(fj(0, i)) < ABS(fj(1, i))) min = 1
               quot = fj(min, i)/aux(min)

               DO l = 2, lmax
                  fj(l, i) = aux(l)*quot
               END DO
            END IF
         END IF
      END DO
      !CPP_OMP END PARALLEL DO
      !$acc end kernels
      
      ! Compute derivatives of spherical Bessel functions
      !CPP_OMP PARALLEL DO default(none) shared(lmax, fj, xvec, dfj, xlim) &
      !CPP_OMP & private(i, l, x, fac)
      !$acc kernels copyin(lmax, xvec,xlim)present(fj,dfj)
      !$acc loop independent private(i, l, x, fac)
      DO i = 1, nx
         x = xvec(i)
         IF (x < xlim) THEN
            fac = 1.0/3.0
            dfj(0, i) = -fj(1, i)
            DO l = 1, lmax
               dfj(l, i) = fac*(l - 0.5*x*x*(l + 2)/(2*l + 3))
               fac = x*fac/(2*l + 3)
            END DO
         ELSE
            dfj(0, i) = -fj(1, i)
            DO l = 1, lmax
               dfj(l, i) = fj(l - 1, i) - (l + 1)*fj(l, i)/x
            END DO
         END IF
      enddo   
      !CPP_OMP END PARALLEL DO
      !$acc end kernels
   END SUBROUTINE d_sphbes
END MODULE m_sphbes
