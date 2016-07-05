MODULE m_fft2d
CONTAINS
  SUBROUTINE fft2d(&
       &                 stars,&
       &                 afft2,bfft2,&
       &                 fg,fgi,fgxy,&
       &                 stride,isn,&
       &                 gfxy )

    !*************************************************************
    !*                                                           *
    !* interface for fg2(star) -- FFT --> gfft (r)     (isn=+1)  *
    !*            or gfft (r)  -- FFT --> fg2(star)    (isn=-1)  *
    !*                                                           *
    !* dimension of gfft2 is (3*stars%k1d x 3*stars%k2d)                     *
    !* afft/bfft contain the real/imaginary part of gfft         *
    !* stars%igfft2(i,1) is the pointer from the G-sphere to stars     *
    !* stars%igfft2(i,2) is the pointer from the G-sphere to fft-grid  *
    !* stars%pgfft2(i)   contains the phases of the G-vectors of sph.  *
    !*                                                           *
    !*************************************************************
#include"cpp_double.h"
    USE m_cfft
    USE m_types
    IMPLICIT NONE
    TYPE(t_stars),INTENT(IN) :: stars
    INTEGER, INTENT (IN) :: isn,stride
    REAL                 :: fg,fgi

    REAL,   INTENT (INOUT):: afft2(0:9*stars%k1d*stars%k2d-1),bfft2(0:9*stars%k1d*stars%k2d-1)
    COMPLEX               :: fgxy(stride,stars%ng2-1)
    REAL,OPTIONAL,INTENT(IN) :: gfxy(0:) !factor to calculate the derivates, i.e. g_x

    !... local variables

    INTEGER i,ifftd2
    REAL  scale
    COMPLEX fg2(stars%n2d)

    ifftd2=9*stars%k1d*stars%k2d
    !
    IF (isn.GT.0) THEN
       !
       !  ---> put stars onto the fft-grid 
       !
       fg2(1) = CMPLX(fg,fgi)
       CALL CPP_BLAS_ccopy(stars%ng2-1,fgxy,stride,fg2(2),1)
       !fg2(2:)=fgxy(1,:)
       afft2=0.0
       bfft2=0.0
       IF (PRESENT(gfxy)) THEN
          DO i=0,stars%kimax2
             afft2(stars%igfft2(i,2))=REAL(fg2(stars%igfft2(i,1))*stars%pgfft2(i))*gfxy(i)
             bfft2(stars%igfft2(i,2))=AIMAG(fg2(stars%igfft2(i,1))*stars%pgfft2(i))*gfxy(i)
          ENDDO
       ELSE 
          DO i=0,stars%kimax2
             afft2(stars%igfft2(i,2))=REAL(fg2(stars%igfft2(i,1))*stars%pgfft2(i))
             bfft2(stars%igfft2(i,2))=AIMAG(fg2(stars%igfft2(i,1))*stars%pgfft2(i))
          ENDDO
       ENDIF
    ENDIF

    !---> now do the fft (isn=+1 : G -> r ; isn=-1 : r -> G)

    CALL cfft(afft2,bfft2,ifftd2,3*stars%k1d,3*stars%k1d,isn)
    CALL cfft(afft2,bfft2,ifftd2,3*stars%k2d,ifftd2,isn)

    IF (isn.LT.0) THEN
       !
       !  ---> collect stars from the fft-grid
       !
       DO i=1,stars%ng2
          fg2(i) = CMPLX(0.0,0.0)
       ENDDO
       scale=1.0/ifftd2
       DO i=0,stars%kimax2
          fg2(stars%igfft2(i,1))=fg2(stars%igfft2(i,1))+ CONJG( stars%pgfft2(i) ) * &
               &                 CMPLX(afft2(stars%igfft2(i,2)),bfft2(stars%igfft2(i,2)))
       ENDDO
       fg=scale*REAL(fg2(1))/stars%nstr2(1)
       fgi=scale*AIMAG(fg2(1))/stars%nstr2(1)
       DO i=2,stars%ng2
          fgxy(1,i-1)=scale*fg2(i)/stars%nstr2(i)
       ENDDO
    ENDIF

  END SUBROUTINE fft2d
END MODULE m_fft2d
