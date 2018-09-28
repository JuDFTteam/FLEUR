      MODULE m_ylm
        IMPLICIT NONE
        private
!---> check whether  or not normalizations are needed
        REAL, ALLOCATABLE, SAVE  :: ynorm(:)
        INTEGER,           SAVE  :: lmaxd = -1  ! initial value
        public ylm4,ylmnorm_init
      CONTAINS 
      SUBROUTINE ylm4(lmax,v,ylm)
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
!************************************************************
      USE m_juDFT
      INTEGER, INTENT (IN) :: lmax
      REAL,    INTENT (IN) :: v(3)
      COMPLEX, INTENT (OUT):: ylm( (lmax+1)**2 )

      REAL,    PARAMETER   :: small = 1.0e-12

      INTEGER l,lm0,m
      REAL    fac,x,y,z,xy,r,rxy,cth,sth,cph,sph,cph2
      REAL    p(0:lmax,0:lmax),c(0:lmax),s(0:lmax)
      COMPLEX ylms


      IF (lmax.GT.lmaxd) THEN
         CALL juDFT_error("lmaxd too small in ylm4")
      ENDIF

!--->    calculate sin and cos of theta and phi
      x = v(1)
      y = v(2)
      z = v(3)
      xy  = x*x + y*y
      r   = sqrt(xy + z*z)
      rxy = sqrt(xy)

      IF (r.GT.small) THEN
         cth = z/r
         sth = rxy/r
      ELSE
         sth = 0.0
         cth = 1.0
      ENDIF
      IF(rxy.GT.small) THEN
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
         fac = -(m+m-1)*fac
         p(m,m) = fac
         p(m+1,m) = (m+m+1)*cth*fac
!--->    recurse upward in l
         DO l = m+2, lmax
            p(l,m)=((l+l-1)*cth*p(l-1,m)-(l+m-1)*p(l-2,m))/(l-m)
         ENDDO
         fac = fac*sth
      ENDDO
      p(lmax,lmax) = -(lmax+lmax-1)*fac

!--->    determine sin and cos of phi
      s(0) = 0.0
      s(1) = sph
      c(0) = 1.0
      c(1) = cph
      cph2 = cph+cph
      DO m = 2, lmax
         s(m) = cph2*s(m-1)-s(m-2)
         c(m) = cph2*c(m-1)-c(m-2)
      ENDDO

!--->    multiply in the normalization factors
      DO l=0,lmax
         ylm(l*(l+1)+1) = ynorm(l*(l+1)+1)*cmplx(p(l,0),0.0)
      ENDDO
      DO m = 1, lmax
         DO l = m, lmax
            lm0 = l*(l+1)+1
            ylms = p(l,m)*cmplx(c(m),s(m))
            ylm(lm0+m) = ynorm(lm0+m)*ylms
            ylm(lm0-m) = conjg(ylms)*ynorm(lm0-m)
         ENDDO
      ENDDO

      RETURN

      END SUBROUTINE ylm4

      SUBROUTINE ylmnorm_init(lmax)
!********************************************************************
!     normalization constants for ylm (internal subroutine has access
!     to lmax and ynorm from above)
!********************************************************************
!$    USE OMP_LIB
      use m_juDFT
      USE m_constants, ONLY : pimach
      INTEGER,INTENT(IN):: lmax
      INTEGER l,lm0,m
      REAL    a,cd,fpi
      IF ( allocated(ynorm) ) DEALLOCATE(ynorm)
      ALLOCATE ( ynorm( (lmax+1)**2 ) )   ! allocate array
      lmaxd = lmax

      fpi = 4.0*pimach()
!$    if (omp_in_parallel()) call juDFT_error(  &
!$          "ylmnorm should not called in parallel",calledby="ylm4.f")
      DO l=0,lmax
         lm0 = l*(l+1) + 1
         cd=1.0
         a = sqrt( (2*l+1)/fpi )
         ynorm(lm0) = a
         DO m=1,l
            cd=cd/( (l+1-m)*(l+m) )
            ynorm(lm0+m)  = a*sqrt(cd)
            ynorm(lm0-m) = ( (-1.0)**m )*ynorm(lm0+m)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ylmnorm_init

      END MODULE m_ylm
