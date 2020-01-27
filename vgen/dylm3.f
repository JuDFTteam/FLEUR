      MODULE m_dylm
      use m_juDFT
c.....------------------------------------------------------------------
c     preparation of dylmt1(=d(ylm)/dtheta),
c     dylmt2(=d(dylmt1)/dtheta),
c     dylmf1, dylmf2 are for fai.
c     dylmtf=d(dylmt1)/df
c     t.a. june, 1996.
c.....------------------------------------------------------------------
      CONTAINS
      SUBROUTINE dylm3(
     >                 lmaxd,lmax,v,ylm,
     <                 dylmt1,dylmt2,dylmf1,dylmf2,dylmtf)
c
      use m_constants
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: lmaxd,lmax
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: v(3)
      COMPLEX, INTENT (IN) ::    ylm( (lmaxd+1)**2 )
      COMPLEX, INTENT (OUT):: dylmt1( (lmaxd+1)**2 )
      COMPLEX, INTENT (OUT):: dylmt2( (lmaxd+1)**2 )
      COMPLEX, INTENT (OUT):: dylmf1( (lmaxd+1)**2 )
      COMPLEX, INTENT (OUT):: dylmf2( (lmaxd+1)**2 )
      COMPLEX, INTENT (OUT):: dylmtf( (lmaxd+1)**2 )
C     ..
C     .. Local Scalars ..
      INTEGER lm1,lm,lm2,lm1m,lmm1m,lmm,lmm1,lmm2
      INTEGER l,m,ll1,llm
      COMPLEX em1f,em2f,ep1f,ep2f
      REAL cph,rxy,small,sph,x,xy,y
C     ..
C     .. Data Statements ..
      DATA small/1.0e-12/
c.....------------------------------------------------------------------
c     ..
      IF (lmax.GT.lmaxd) THEN
        WRITE (6,*) lmax,lmaxd
         CALL juDFT_error("lmax.GT.lmaxd",calledby="dylm3")
      ENDIF

c--->    calculate sin and cos of phi
      x = v(1)
      y = v(2)

      xy = x*x + y*y
      rxy = sqrt(xy)
c
      IF (rxy.gt.small) THEN
         cph = x/rxy
         sph = y/rxy
      ELSE
         cph = 1.e0
         sph = 0.e0
      ENDIF

      ep1f=cmplx(cph,sph)
      em1f=conjg(ep1f)
      ep2f=ep1f*ep1f
      em2f=em1f*em1f
c
      DO 21 l=0,lmax
         ll1 = l*(l+1) 

         DO m=-l,l
            llm = ll1 + m + 1
            dylmt1(llm) = cmplx(0.0,0.0)
            dylmt2(llm) = cmplx(0.0,0.0)
         ENDDO

         DO 23 m=-l,l
            llm = ll1 + m + 1

            lmm1m = l - m - 1
            lmm   = l - m
            lmm1  = l - m + 1
            lmm2  = l - m + 2
            lm1m  = l + m - 1
            lm    = l + m
            lm1   = l + m + 1
            lm2   = l + m + 2

            dylmt2(llm)=dylmt2(llm) -
     +                      (lmm*lm1+lmm1*lm)/4.e0*ylm(llm)

            IF (m+2.le.l) dylmt2(llm)=dylmt2(llm) +
     +           sqrt(real(lmm1m*lmm*lm1*lm2))/4*ylm(llm+2)*em2f

            IF (m+1.le.l) dylmt1(llm)=dylmt1(llm) +
     +           sqrt(real(lmm*lm1))/2*ylm(llm+1)*em1f

            IF (m-1.ge.-l) dylmt1(llm)=dylmt1(llm) -
     +           sqrt(real(lm*lmm1))/2*ylm(llm-1)*ep1f

            IF (m-2.ge.-l) dylmt2(llm)=dylmt2(llm) +
     +           sqrt(real(lmm1*lmm2*lm1m*lm))/4*ylm(llm-2)*ep2f

   23    ENDDO

         DO m=-l,l
            llm = ll1 + m + 1
            dylmf1(llm) = ImagUnit * m *    ylm(llm)
            dylmf2(llm) = -m * m *    ylm(llm)
            dylmtf(llm) = ImagUnit * m * dylmt1(llm)
         ENDDO

   21 ENDDO

      RETURN
      END SUBROUTINE dylm3
      END MODULE m_dylm
