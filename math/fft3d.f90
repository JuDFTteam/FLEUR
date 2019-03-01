MODULE m_fft3d
CONTAINS
   SUBROUTINE fft3d(&
  &                 afft, bfft, fg3,&
  &                 stars, isn, scaled)

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
      USE m_types
      USE m_fft_interface
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: isn
      TYPE(t_stars), INTENT(IN):: stars
      REAL, INTENT(INOUT) :: afft(0:27*stars%mx1*stars%mx2*stars%mx3 - 1)
      REAL, INTENT(INOUT) :: bfft(0:27*stars%mx1*stars%mx2*stars%mx3 - 1)
      COMPLEX                 :: fg3(stars%ng3)
      LOGICAL, INTENT(IN), OPTIONAL :: scaled !< determines if coefficients are scaled by stars%nstr

      INTEGER i, ifftd
      REAL scale
      COMPLEX ctmp
      LOGICAL forw
      INTEGER length_zfft(3)
      complex :: zfft(0:27*stars%mx1*stars%mx2*stars%mx3 - 1)

      ifftd = 27*stars%mx1*stars%mx2*stars%mx3

      IF (isn > 0) THEN
!
!  ---> put stars onto the fft-grid
!
         afft = 0.0
         bfft = 0.0
         DO i = 0, stars%kimax
            ctmp = fg3(stars%igfft(i, 1))*stars%pgfft(i)
            afft(stars%igfft(i, 2)) = real(ctmp)
            bfft(stars%igfft(i, 2)) = aimag(ctmp)
         ENDDO
      ENDIF

!---> now do the fft (isn=+1 : G -> r ; isn=-1 : r -> G)

      zfft = cmplx(afft, bfft)
      if (isn == -1) then
         forw = .true.
      else
         forw = .false.
      end if
      length_zfft(1) = 3*stars%mx1
      length_zfft(2) = 3*stars%mx2
      length_zfft(3) = 3*stars%mx3
      call fft_interface(3, length_zfft, zfft, forw)

      afft = real(zfft)
      bfft = aimag(zfft)

      IF (isn < 0) THEN
!
!  ---> collect stars from the fft-grid
!
         DO i = 1, stars%ng3
            fg3(i) = cmplx(0.0, 0.0)
         ENDDO
         DO i = 0, stars%kimax
            fg3(stars%igfft(i, 1)) = fg3(stars%igfft(i, 1)) + CONJG(stars%pgfft(i))* &
       &                 zfft(stars%igfft(i, 2))
         ENDDO
         scale = 1.0/ifftd
         IF (PRESENT(scaled)) THEN
            IF (scaled) THEN
               fg3 = scale*fg3/stars%nstr
            ELSE
               fg3 = scale*fg3
            ENDIF
         ELSE
            fg3 = scale*fg3/stars%nstr
         ENDIF
      ENDIF

   END SUBROUTINE fft3d
END MODULE m_fft3d
