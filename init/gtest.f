      MODULE m_gtest
      use m_juDFT
!*********************************************************************
!     test the gaussian integration mesh for this l
!     1) integrate each ylm, 2) get normalization of each ylm,
!     3) (ylm(lmax,m=lmax | ylm) for all lm
!*********************************************************************
      CONTAINS
      SUBROUTINE gtest(
     >                 lmax,ngpts,vgauss,wt)

      USE m_ylm
      IMPLICIT NONE
!
      INTEGER, INTENT (IN) :: lmax,ngpts
      REAL,    INTENT (IN) :: vgauss(3,*),wt(*) !   gaussian integration points

      INTEGER :: lm,lmmax,nn
      LOGICAL :: lerr
      COMPLEX :: temp
      REAL    :: v(3)
      COMPLEX, DIMENSION((lmax+1)**2) :: ylm,s1,s2,s3

      REAL, PARAMETER :: del = 1.e-10

      lmmax = (lmax+1)**2
      DO lm = 1, lmmax
         s1(lm)=cmplx(0.0,0.0)
         s2(lm)=cmplx(0.0,0.0)
         s3(lm)=cmplx(0.0,0.0)
      ENDDO

!---> loop over points
      DO nn = 1, ngpts
         v(:) = vgauss(:,nn)
         CALL ylm4(
     >             lmax,v,
     <             ylm)
         DO lm = 1, lmmax
            s1(lm) = s1(lm) + wt(nn)*ylm(lm)
            s2(lm) = s2(lm) + wt(nn)*conjg(ylm(lm))*ylm(lm)
            s3(lm) = s3(lm) + wt(nn)*conjg(ylm(lmmax))*ylm(lm)
         ENDDO
      ENDDO

!---> check for non-zero elements
      lerr = .false.
      DO lm=1,lmmax
         IF ((abs(s1(lm)).GT.del).AND.(lm.GT.1)) THEN
            WRITE (6,1000) lm,s1(lm)
            lerr = .true.
         ENDIF
 1000    FORMAT (20x,' integration of s.h. lm=',i3,1p,2e14.6)

         temp = s2(lm) - cmplx(1.0,0.0)
         IF (abs(temp).GT.del) THEN
            WRITE (6,1001) lm,temp
            lerr = .true.
         ENDIF
 1001    FORMAT (20x,' normalization error lm=',i3,1p,2e14.6)

         IF (abs(s3(lm)).GT.del) THEN
            IF(lm.ne.lmmax) THEN
              WRITE (6,1002) lm,s3(lm)
              lerr = .true.
            ELSEIF( abs(abs(s3(lm))-1.0).GT.del ) THEN
              WRITE (6,1002) lm,s3(lm)
              lerr = .true.
            ENDIF
         ENDIF
 1002    FORMAT (20x,' (ylm(lmax,lmax) | ylm)    lm=',i3,1p,2e14.6)
      ENDDO

      IF (lerr) THEN
        WRITE (6,'(/,'' error in gaussian points'')')
         CALL juDFT_error("Error in gaussian points",calledby="gtest")
      ENDIF

      RETURN
      END SUBROUTINE gtest
      END MODULE m_gtest
