      MODULE m_fft3d
      CONTAINS
      SUBROUTINE fft3d(&
     &                 afft,bfft,fg3,&
     &                 stars,isn,scaled)

!************************************************************
!*                                                          *
!* interface for fg3(star) -- FFT --> (a,b)fft (r) (isn=+1) *
!*         or (a,b)fft (r) -- FFT --> fg3(star)    (isn=-1) *
!*                                                          *
!* dimension of (a,b)fft is (3*k1d x 3*k2d x 3*k3d)         *
!* afft and bfft contain the real/imaginary part of the FFT *
!* igfft(i,1) is the pointer from the G-sphere to stars     *
!* igfft(i,2) is the pointer from the G-sphere to fft-grid  *
!* pgfft(i)   contains the phases of the G-vectors of sph.  *
!*                                                          *
!************************************************************
      USE m_cfft
      USE m_types
      IMPLICIT NONE
    
      INTEGER, INTENT (IN) :: isn
      TYPE(t_stars),INTENT(IN):: stars
      REAL,    INTENT (INOUT) :: afft(0:27*stars%k1d*stars%k2d*stars%k3d-1)
      REAL,    INTENT (INOUT) :: bfft(0:27*stars%k1d*stars%k2d*stars%k3d-1)
      COMPLEX                 :: fg3(stars%ng3)
      LOGICAL,INTENT(IN),OPTIONAL :: scaled !< determines if coefficients are scaled by stars%nstr

      INTEGER i,ifftd
      REAL scale

      ifftd=27*stars%k1d*stars%k2d*stars%k3d
     
      IF (isn.GT.0) THEN
!
!  ---> put stars onto the fft-grid 
!
        afft=0.0
        bfft=0.0
        DO i=0,stars%kimax
          afft(stars%igfft(i,2))=real(fg3(stars%igfft(i,1))*&
     &                                   stars%pgfft(i))
          bfft(stars%igfft(i,2))=aimag(fg3(stars%igfft(i,1))*&
     &                                   stars%pgfft(i))
        ENDDO
      ENDIF

!---> now do the fft (isn=+1 : G -> r ; isn=-1 : r -> G)

      CALL cfft(afft,bfft,ifftd,3*stars%k1d,3*stars%k1d,isn)
      CALL cfft(afft,bfft,ifftd,3*stars%k2d,9*stars%k1d*stars%k2d,isn)
      CALL cfft(afft,bfft,ifftd,3*stars%k3d,ifftd,isn)

      IF (isn.LT.0) THEN
!
!  ---> collect stars from the fft-grid
!
        DO i=1,stars%ng3
          fg3(i) = cmplx(0.0,0.0)
        ENDDO
        DO i=0,stars%kimax
          fg3(stars%igfft(i,1))=fg3(stars%igfft(i,1))+  stars%pgfft(i)*&
     &                cmplx(afft(stars%igfft(i,2)),bfft(stars%igfft(i,2)))
        ENDDO
        scale=1.0/ifftd
        IF (PRESENT(scaled)) THEN
           IF (scaled) THEN
               fg3=scale*fg3/stars%nstr
           ELSE
               fg3=scale*fg3    
           ENDIF
        ELSE
           fg3=scale*fg3/stars%nstr
        ENDIF
      ENDIF

      END SUBROUTINE fft3d
      END MODULE m_fft3d
