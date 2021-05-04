MODULE m_ylm
   use iso_c_binding
   IMPLICIT NONE
   PRIVATE
   PUBLIC ylm4, ylm4_batched
!************************************************************
!     generate the spherical harmonics for the vector v
!     using a stable upward recursion in l.  (see notes
!     by m. weinert.)
!          m.weinert   january 1982
!     modified by R. Podloucky (added in ynorm); July 1989
!     cleaned up    mw 1995
!
!     modified to make use of f90 constructs. note that
!     the normalization is an internal subroutine and hence
!     can only be called from here. also, no need to dimension
!     arrays for ynorm, done dynamically.          mw 1999
!
!     GPU version added
!             U.Alekseeva      Oktober 2018
!************************************************************
CONTAINS
   SUBROUTINE ylm4(lmax, v, ylm)
!      USE m_juDFT
      INTEGER, VALUE, INTENT(IN)   :: lmax
      REAL, INTENT(IN), target     :: v(3)
      COMPLEX, INTENT(OUT), target :: ylm((lmax + 1)**2)

      real, pointer    :: v_2d(:,:)
      complex, pointer :: ylm_2d(:,:)

      call c_f_pointer(c_loc(v),   v_2d,   [3,1])
      call c_f_pointer(c_loc(ylm), ylm_2d, [(lmax + 1)**2, 1])
      call ylm4_batched(lmax, v_2d, ylm_2d)
   END SUBROUTINE ylm4

   subroutine ylm4_batched(lmax, v, ylm)
      implicit none
      INTEGER, VALUE, INTENT(IN) :: lmax
      REAL, INTENT(IN) :: v(:,:)
      COMPLEX, INTENT(OUT):: ylm(:,:)

      INTEGER  :: l, lm0, m, i
      REAL     :: fac, x, y, z, xy, r, rxy, cth, sth, cph, sph, cph2
      REAL     :: c(0:max(lmax,1)), s(0:max(lmax,1))
      REAL     :: p(0:lmax, 0:lmax)
      COMPLEX  :: ylms
      REAL     :: ynorm_dev((lmax + 1)**2)
      REAL     :: a, cd
      real, parameter :: fpi = 4.0*3.1415926535897932, small = 1.0e-12

      !--->    calculate norm
      DO l = 0, lmax
         lm0 = l*(l + 1) + 1
         cd = 1.0
         a = sqrt((2*l + 1)/fpi)
         ynorm_dev(lm0) = a
         DO m = 1, l
            cd = cd/((l + 1 - m)*(l + m))
            ynorm_dev(lm0 + m) = a*sqrt(cd)
            ynorm_dev(lm0 - m) = ((-1.0)**m)*ynorm_dev(lm0 + m)
         ENDDO
      ENDDO


!--->    calculate sin and cos of theta and phi
      !$OMP parallel do default(none) &
      !$OMP private(i, x, y, z, xy, r, rxy, cth, sth, fac, m, l, p, c ,s, ylms, lm0, cph, sph, cph2)&
      !$OMP shared(ylm, ynorm_dev, lmax, v)
      do i = 1,size(v,2)
         x = v(1,i)
         y = v(2,i)
         z = v(3,i)
         xy = x*x + y*y
         r = sqrt(xy + z*z)
         rxy = sqrt(xy)

         IF (r .GT. small) THEN
            cth = z/r
            sth = rxy/r
         ELSE
            sth = 0.0
            cth = 1.0
         ENDIF
         IF (rxy .GT. small) THEN
            cph = x/rxy
            sph = y/rxy
         ELSE
            cph = 1.0
            sph = 0.0
         ENDIF

   !---> generate associated legendre functions for m.ge.0
         fac = 1.0
   !---> loop over m values
         DO m = 0, lmax - 1
            fac = -(m + m - 1)*fac
            p(m, m) = fac
            p(m + 1, m) = (m + m + 1)*cth*fac
   !--->    recurse upward in l
            DO l = m + 2, lmax
               p(l, m) = ((l + l - 1)*cth*p(l - 1, m) - (l + m - 1)*p(l - 2, m))/(l - m)
            ENDDO
            fac = fac*sth
         ENDDO
         p(lmax, lmax) = -(lmax + lmax - 1)*fac

   !--->    determine sin and cos of phi
         s(0) = 0.0
         s(1) = sph
         c(0) = 1.0
         c(1) = cph
         cph2 = cph + cph
         DO m = 2, lmax
            s(m) = cph2*s(m - 1) - s(m - 2)
            c(m) = cph2*c(m - 1) - c(m - 2)
         ENDDO

   !--->    multiply in the normalization factors
         DO l = 0, lmax
            ylm(l*(l + 1) + 1, i) = ynorm_dev(l*(l + 1) + 1)*cmplx(p(l, 0), 0.0)
         ENDDO
         DO m = 1, lmax
            DO l = m, lmax
               lm0 = l*(l + 1) + 1
               ylms = p(l, m)*cmplx(c(m), s(m))
               ylm(lm0 + m, i) = ynorm_dev(lm0 + m)*ylms
               ylm(lm0 - m, i) = conjg(ylms)*ynorm_dev(lm0 - m)
            ENDDO
         ENDDO
      enddo
      !$OMP end parallel do
   end subroutine ylm4_batched
END MODULE m_ylm
