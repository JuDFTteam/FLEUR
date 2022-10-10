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

   SUBROUTINE dfpt_convol_direct(stars, starsq, pw, pwq, pwwq)
      ! TODO: Should probably be replaced by a "finer" function with full
      !       G-grid for ustep(1)
      USE m_types_fftGrid
      USE m_juDFT
      USE m_types_stars

      IMPLICIT NONE

      TYPE(t_stars), INTENT(IN) :: stars, starsq

      COMPLEX, INTENT(IN)    :: pw(:), pwq(:)
      COMPLEX, INTENT(INOUT) :: pwwq(:)

      TYPE(t_fftgrid) :: fftgrid, fftgridq

      CALL fftgrid%init((/3*stars%mx1,3*stars%mx2,3*stars%mx3/))
      CALL fftgridq%init((/3*starsq%mx1,3*starsq%mx2,3*starsq%mx3/))

      IF (SIZE(fftgridq%grid)/=SIZE(fftgrid%grid)) CALL judft_error("Size mismatch in dfpt_convol (3)!")

      CALL fftgrid%putFieldOnGrid(stars,pw)
      CALL fftgridq%putFieldOnGrid(starsq,pwq)
      CALL fftgrid%perform_fft(forward=.false.)
      CALL fftgridq%perform_fft(forward=.false.)

      fftgrid%grid = fftgrid%grid*fftgridq%grid

      CALL fftgrid%perform_fft(forward=.true.)
      CALL fftgrid%takeFieldFromGrid(starsq,pwwq)
      pwwq = pwwq*starsq%nstr

   END SUBROUTINE dfpt_convol_direct

   SUBROUTINE dfpt_convol_big(resultstar, stars, starsfull, pw, pwfull, pww)
      USE m_types_fftGrid
      USE m_juDFT
      USE m_types_stars

      IMPLICIT NONE

      INTEGER,       INTENT(IN) :: resultstar
      TYPE(t_stars), INTENT(IN) :: stars, starsfull

      COMPLEX, INTENT(IN)    :: pw(:), pwfull(:)
      COMPLEX, INTENT(INOUT) :: pww(:)

      TYPE(t_fftgrid) :: fftgrid, fftgridfin

      CALL fftgrid%init((/3*stars%mx1,3*stars%mx2,3*stars%mx3/))

      IF (SIZE(pwfull)/=SIZE(fftgrid%grid)) CALL judft_error("Size mismatch in dfpt_convol (4)!")

      CALL fftgrid%putFieldOnGrid(stars,pw)
      CALL fftgrid%perform_fft(forward=.false.)

      IF (resultstar==1) THEN
         fftgrid%grid = fftgrid%grid*pwfull
         CALL fftgrid%perform_fft(forward=.true.)
         CALL fftgrid%takeFieldFromGrid(stars,pww)
         pww = pww*stars%nstr
      ELSE IF (resultstar==2) THEN
         CALL fftgridfin%init((/3*starsfull%mx1,3*starsfull%mx2,3*starsfull%mx3/))
         fftgridfin%grid = fftgrid%grid*pwfull
         CALL fftgridfin%perform_fft(forward=.true.)
         CALL fftgridfin%takeFieldFromGrid(starsfull,pww)
         pww = pww*starsfull%nstr
      ELSE
         CALL judft_error("Impossible star index!")
      END IF

   END SUBROUTINE dfpt_convol_big

END MODULE m_convol
