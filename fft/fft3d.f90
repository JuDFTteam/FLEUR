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
      USE m_types_fftGrid
      IMPLICIT none

      INTEGER, INTENT(IN)       :: isn
      TYPE(t_stars), INTENT(IN) :: stars
      REAL, INTENT(INOUT)       :: afft(:)
      REAL, INTENT(INOUT)       :: bfft(:)
      COMPLEX                   :: fg3(stars%ng3) !Sometimes we call this with an intent(in) variable if isn>0 (This is somewhat unsafe)
      LOGICAL, INTENT(IN), OPTIONAL :: scaled ! < determines if coefficients are scaled by stars%nstr

      TYPE(t_fftgrid) :: fftgrid
      call fftgrid%init((/stars%mx1,stars%mx2,stars%mx3/))

      IF (isn > 0) THEN
          call fftgrid%putFieldOnGrid(stars,fg3)
      ELSE 
         fftgrid%grid=cmplx(afft,bfft)
      ENDIF       
     
      call fftgrid%perform_fft(forward=(isn<0))
!---> now do the fft (isn=+1 : G -> r ; isn=-1 : r -> G)

      if (isn >0) THEN
        afft = real(fftgrid%grid)
        bfft = aimag(fftgrid%grid)
      else
         call fftgrid%takeFieldFromGrid(stars,fg3)
         IF (PRESENT(scaled)) THEN
            IF (scaled) fg3 = fg3/stars%nstr
         ELSE
            fg3 = fg3/stars%nstr
         ENDIF
      ENDIF   

   END SUBROUTINE fft3d
END MODULE m_fft3d
