MODULE m_fft2d
CONTAINS
   SUBROUTINE fft2d(stars,afft2,bfft2,fg,isn,firstderiv,secondderiv,cell )
      !!
      !!
    !*************************************************************
    !*                                                           *
    !* interface for fg2(star) -- FFT --> gfft (r)     (isn=+1)  *
    !*            or gfft (r)  -- FFT --> fg2(star)    (isn=-1)  *
    !*                                                           *
    !* dimension of gfft2 is (3*stars%mx1 x 3*stars%mx2)                     *
    !* afft/bfft contain the real/imaginary part of gfft         *
    !* stars%igfft2(i,1) is the pointer from the G-sphere to stars     *
    !* stars%igfft2(i,2) is the pointer from the G-sphere to fft-grid  *
    !* stars%pgfft2(i)   contains the phases of the G-vectors of sph.  *
    !*                                                           *
    !*************************************************************
      USE m_types_fftgrid
      USE m_types

      IMPLICIT NONE

      TYPE(t_stars), INTENT(IN) :: stars
      
      INTEGER, INTENT (IN)   :: isn
      REAL                   :: afft2(0:9*stars%mx1*stars%mx2-1),bfft2(0:9*stars%mx1*stars%mx2-1)
      COMPLEX                :: fg(:)
    
      REAL,         OPTIONAL, INTENT(IN) :: firstderiv(3),secondderiv(3)
      TYPE(t_cell), OPTIONAL, INTENT(IN) :: cell
    
      TYPE(t_fftgrid) :: grid
      
      INTEGER :: i
   
      CALL grid%init([3*stars%mx1,3*stars%mx2,1])

      IF (isn>0) THEN
       	!  ---> put stars onto the fft-grid
       	CALL grid%putFieldOnGrid(stars,fg,cell,firstderiv=firstderiv,secondderiv=secondderiv,l_2d=.TRUE.)
    	ELSE
       	grid%grid=cmplx(afft2,bfft2)
    	END IF

    	CALL grid%perform_fft(forward=(isn<0))
    
    	IF (isn>0) THEN
        	afft2 =  REAL(grid%grid)
        	bfft2 = AIMAG(grid%grid)
      ELSE
         CALL grid%takeFieldFromGrid(stars,fg,l_2d=.TRUE.)
         !Scaling by stars%nstr is already done in previous call
         !IF (PRESENT(scaled)) THEN
         !   IF (.not.scaled) fg3 = fg3*stars%nstr
         !ENDIF
      END IF   
    	IF (isn<0) THEN
       	!  ---> collect stars from the fft-grid
			CALL grid%takeFieldFromGrid(stars, fg, l_2d=.TRUE.)
    	END IF
  END SUBROUTINE fft2d
END MODULE m_fft2d
