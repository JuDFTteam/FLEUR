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
   
         TYPE(t_stars), INTENT(IN)     :: stars
   
         COMPLEX, INTENT(IN),OPTIONAL   :: ag3(:)
         COMPLEX, INTENT(INOUT)         :: fg3(:)
        
         TYPE(t_fftgrid) :: fftgrid
         
         call fftgrid%init((/3*stars%mx1,3*stars%mx2,3*stars%mx3/))

         if (size(fftgrid%grid).ne.size(stars%ufft)) call judft_error("Bug in t_stars%convol")
         if (present(ag3)) THEN
           call fftgrid%putFieldOnGrid(stars,ag3)
         else
           !In place version
           call fftgrid%putFieldOnGrid(stars,fg3)
         endif
         call fftgrid%perform_fft(forward=.false.)
      
         fftgrid%grid = fftgrid%grid*stars%ufft
         
         call fftgrid%perform_fft(forward=.true.)
         call fftgrid%takeFieldFromGrid(stars,fg3)
         fg3 = fg3*stars%nstr
       
      END SUBROUTINE convol

END MODULE m_convol
