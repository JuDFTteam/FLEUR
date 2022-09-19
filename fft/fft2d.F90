MODULE m_fft2d
CONTAINS
  SUBROUTINE fft2d(&
       &                 stars,&
       &                 afft2,bfft2,&
       &                 fg,fgi,fgxy,&
       &                 isn,&
       &                 firstderiv,secondderiv,cell )

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
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_cell),INTENT(IN),OPTIONAL:: cell
    INTEGER, INTENT (IN) :: isn
    REAL                 :: fg,fgi

    REAL                  :: afft2(0:9*stars%mx1*stars%mx2-1),bfft2(0:9*stars%mx1*stars%mx2-1)
    COMPLEX               :: fgxy(:)
    REAL,OPTIONAL,INTENT(IN):: firstderiv(3),secondderiv(3)
    
    !... local variables

   TYPE(t_fftgrid) :: grid
   INTEGER i
   COMPLEX fg2(stars%ng2)

   
   call grid%init([3*stars%mx1,3*stars%mx2,1])

    IF (isn>0) THEN
       !  ---> put stars onto the fft-grid
       fg2(1) = CMPLX(fg,fgi)
       fg2(2:)=fgxy(:)

       call grid%putFieldOnGrid(stars,fg2,cell,firstderiv=firstderiv,secondderiv=secondderiv,l_2d=.true.)
    else
       grid%grid=cmplx(afft2,bfft2)
    endif

    call grid%perform_fft(forward=(isn<0))
    
    if (isn >0) THEN
        afft2 = real(grid%grid)
        bfft2 = aimag(grid%grid)
      else
         call grid%takeFieldFromGrid(stars,fg2,l_2d=.true.)
         !Scaling by stars%nstr is already done in previous call
         !IF (PRESENT(scaled)) THEN
         !   IF (.not.scaled) fg3 = fg3*stars%nstr
         !ENDIF
      ENDIF   
    IF (isn.LT.0) THEN
       !
       !  ---> collect stars from the fft-grid
       !
       call grid%takeFieldFromGrid(stars, fg2, l_2d=.true.)
       fg=REAL(fg2(1))
       fgi=AIMAG(fg2(1))
       DO i=2,stars%ng2
          fgxy(i-1)=fg2(i)
       ENDDO
    ENDIF

  END SUBROUTINE fft2d
END MODULE m_fft2d
