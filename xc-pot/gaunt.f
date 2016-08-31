      MODULE m_gaunt
c*********************************************************************
c     Modified module to include old gaunt2 subroutine
c     the private arrays are allocated and computed in the first call to gaunt1
c                                            Daniel Wortmann
c*********************************************************************
      PRIVATE
      INTEGER,SAVE         :: lmaxdp
      REAL,SAVE,ALLOCATABLE::w(:),yr(:,:)
      PUBLIC gaunt1,gaunt2
      CONTAINS
      REAL FUNCTION gaunt1(lp,l,ls,mp,m,ms,lmaxd)
c*********************************************************************
c     gaunt computes the integral of conjg(y(lp,mp))*y(l,m)*y(ls,ms)
c     for lp+l+ls .lt. 2*ngntd
c     using gaussian quadrature as given by
c     m. abramowitz and i.a. stegun, handbook of mathematical functions,
c     nbs applied mathematics series 55 (1968), pages 887 and 916
c     m. weinert and e. wimmer
c     northwestern university march 1980
c     modified to use calculated points and weights
c     to make it dynamic.   (m.w.  jan. 1982)
c*********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER,INTENT(IN):: l,lp,ls,m,mp,ms,lmaxd
C     ..
C     .. Local Scalars ..
      REAL zero
      INTEGER i,il,ilp,ils,n
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC mod
C     ..
C     .. Data statements ..
      DATA zero/0.0e0/
C     ..
      n= (3*lmaxd)/4+1
c     
      !check if this is first call to subroutine
      IF (.not.ALLOCATED(YR)) CALL gaunt2(lmaxd)
      !check if the previous call of the subroutine was with the same lmaxd
      IF( lmaxd .ne. lmaxdp ) THEN
        DEALLOCATE(yr,w)
        CALL gaunt2(lmaxd)
      END IF
  
      gaunt1 = zero
      IF (mp.NE. (m+ms)) RETURN
      IF (MOD((l+lp+ls),2).EQ.1) RETURN
      IF ((l+lp-ls).LT.0) RETURN
      IF ((l-lp+ls).LT.0) RETURN
      IF ((lp-l+ls).LT.0) RETURN
      il = l* (l+1) + m + 1
      ilp = lp* (lp+1) + mp + 1
      ils = ls* (ls+1) + ms + 1
      DO i = 1,n
         gaunt1 = gaunt1 + w(i)*yr(i,ilp)*yr(i,il)*yr(i,ils)
      END DO
      RETURN
      END FUNCTION

c     private subroutine for initializing the private arrays!
      SUBROUTINE gaunt2(
     >                  lmaxd)
c**********************************************************************
c     sets up values needed for gaunt1
c        m. weinert  january 1982
c**********************************************************************
      USE m_constants, ONLY : pimach
      USE m_grule
!$    USE omp_lib
      IMPLICIT NONE
      
      INTEGER, INTENT (IN)  :: lmaxd
C     ..
C     .. Local Scalars ..
      REAL a,cd,cth,fac,fpi,rf,sgm,sth,t
      INTEGER k,l,lm,lomax,m,nn
      INTEGER n,lmax1d
C     ..
C     .. Local Arrays ..
      REAL p(0:lmaxd+1,0:lmaxd+1),x((3*lmaxd)/4+1)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
      if (allocated(w)) return
      !$ if (omp_in_parallel()) STOP "BUG IN GAUNT!!"
      ALLOCATE(w((3*lmaxd)/4+1),yr((3*lmaxd)/4+1,(lmaxd+1)**2))

      n = (3*lmaxd)/4+1
      lmaxdp = lmaxd
      lmax1d = lmaxd+1
c
      fpi = 4.0 * pimach()
      rf = fpi** (1./3.)
      lomax = lmax1d - 1
c--->    obtain gauss-legendre points and weights
      nn = 2*n
      CALL grule(nn,x,w)
c--->    generate associated legendre functions for m.ge.0
      DO  k = 1,n
         cth = x(k)
         sth = sqrt(1.0-cth*cth)
         fac = 1.0
c--->    loop over m values
         DO  m = 0,lomax
            fac = - (2*m-1)*fac
            p(m,m) = fac
            p(m+1,m) = (2*m+1)*cth*fac
c--->    recurse upward in l
            DO  l = m + 2,lomax
               p(l,m) = ((2*l-1)*cth*p(l-1,m)- (l+m-1)*p(l-2,m))/ (l-m)
            ENDDO
            fac = fac*sth
         ENDDO
c--->    multiply in the normalization factors
         DO  l = 0,lomax
            a = rf*sqrt((2*l+1)/fpi)
            cd = 1
            lm = l* (l+1) + 1
            yr(k,lm) = a*p(l,0)
            sgm = -1.
            DO  m = 1,l
               t = (l+1-m)* (l+m)
               cd = cd/t
               yr(k,lm+m) = a*sqrt(cd)*p(l,m)
               yr(k,lm-m) = sgm*a*sqrt(cd)*p(l,m)
               sgm = -sgm
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE

      END MODULE
