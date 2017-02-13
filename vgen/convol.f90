      MODULE m_convol
      CONTAINS
      SUBROUTINE convol(&
     &                  stars, &
     &                  fg3,&
     &                  ag3&
     &                   )

!************************************************************
!*                                                          *
!* calculate f(G) = \sum_G' U(G' - G) a(G')                 *
!*                                                          *
!*       ag3(star) -- FFT --> gfft(r,1)                     *
!*                            gfft(r,1)=gfft(r,1) * U (r)   *
!*       fg3(star) <- FFT --- gfft(r,1)                     *
!*                                                          *
!* dimension of gfft is (3*stars%mx1 x 3*stars%mx2 x 3*stars%mx3)             *
!*                                                          *
!************************************************************
      USE m_types
      USE m_fft3d
      IMPLICIT NONE

      TYPE(t_stars),INTENT(IN) :: stars
      COMPLEX, INTENT (IN)     :: ag3(stars%ng3)
      COMPLEX, INTENT (OUT)    :: fg3(stars%ng3)

      INTEGER i,ifftd
      REAL, ALLOCATABLE :: gfft(:,:)

      ifftd=27*stars%mx1*stars%mx2*stars%mx3

      ALLOCATE (gfft(0:ifftd-1,2))

      CALL fft3d(&
     &           gfft(0,1),gfft(0,2),&
     &           ag3,&
     &           stars,+1) 

      DO i=0,ifftd-1
        gfft(i,:)=gfft(i,:)*stars%ufft(i)
      ENDDO

      CALL fft3d(&
     &           gfft(0,1),gfft(0,2),&
     &           fg3,&
     &           stars,-1) 

      fg3(:stars%ng3)=fg3(:stars%ng3)*stars%nstr(:stars%ng3)

      DEALLOCATE (gfft)

      END SUBROUTINE convol
      END MODULE m_convol
