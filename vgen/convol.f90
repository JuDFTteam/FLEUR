MODULE m_convol
CONTAINS
   SUBROUTINE convol(stars, fg3, ag3)

   !************************************************************
   !*                                                          *
   !* calculate f(G) = \sum_G' U(G - G') a(G')                 *
   !*                                                          *
   !* U is already given on the real space mesh as U(r)        *
   !*                                                          *
   !*       ag3(star) -- FFT --> gfft(r,1)                     *
   !*                            gfft(r,1)=gfft(r,1) * U (r)   *
   !*       fg3(star) <- FFT --- gfft(r,1)                     *
   !*                                                          *
   !* dimension of gfft is                                     *
   !* (3*stars%mx1 x 3*stars%mx2 x 3*stars%mx3)                *
   !*                                                          *
   !************************************************************
      USE m_types_fftGrid
      USE m_juDFT
      USE m_types_stars

      IMPLICIT NONE

      TYPE(t_stars), INTENT(IN) :: stars

      COMPLEX, INTENT(INOUT)          :: fg3(:)
      COMPLEX, INTENT(IN),   OPTIONAL :: ag3(:)

      TYPE(t_fftgrid) :: fftgrid

      CALL fftgrid%init((/3*stars%mx1,3*stars%mx2,3*stars%mx3/))

      IF (SIZE(fftgrid%grid)/=SIZE(stars%ufft)) CALL judft_error("Bug in t_stars%convol")
      IF (PRESENT(ag3)) THEN
        CALL fftgrid%putFieldOnGrid(stars,ag3)
      ELSE
        !In place version
        CALL fftgrid%putFieldOnGrid(stars,fg3)
      END IF
      CALL fftgrid%perform_fft(forward=.false.)

      fftgrid%grid = fftgrid%grid*stars%ufft

      CALL fftgrid%perform_fft(forward=.true.)
      CALL fftgrid%takeFieldFromGrid(stars,fg3)
      fg3 = fg3*stars%nstr

   END SUBROUTINE convol

   SUBROUTINE dfpt_convol(stars, starsq, pw, pwq, pww)
      USE m_types_fftGrid
      USE m_juDFT
      USE m_types_stars

      IMPLICIT NONE

      TYPE(t_stars), INTENT(IN) :: stars, starsq

      COMPLEX, INTENT(IN)    :: pw(:), pwq(:)
      COMPLEX, INTENT(INOUT) :: pww(:)

      TYPE(t_fftgrid) :: fftgrid, fftgridq

      CALL fftgrid%init((/3*stars%mx1,3*stars%mx2,3*stars%mx3/))
      CALL fftgridq%init((/3*starsq%mx1,3*starsq%mx2,3*starsq%mx3/))

      IF (SIZE(fftgrid%grid)/=SIZE(stars%ufft)) CALL judft_error("Size mismatch in dfpt_convol (1)!")
      IF (SIZE(fftgridq%grid)/=SIZE(starsq%ufft1)) CALL judft_error("Size mismatch in dfpt_convol (2)!")
      IF (SIZE(fftgridq%grid)/=SIZE(fftgrid%grid)) CALL judft_error("Size mismatch in dfpt_convol (3)!")

      CALL fftgrid%putFieldOnGrid(stars,pw)
      CALL fftgridq%putFieldOnGrid(starsq,pwq)
      CALL fftgrid%perform_fft(forward=.false.)
      CALL fftgridq%perform_fft(forward=.false.)

      fftgrid%grid = fftgrid%grid*starsq%ufft1 + fftgridq%grid*stars%ufft

      CALL fftgrid%perform_fft(forward=.true.)
      CALL fftgrid%takeFieldFromGrid(stars,pww)
      pww = pww*stars%nstr

   END SUBROUTINE dfpt_convol

END MODULE m_convol
