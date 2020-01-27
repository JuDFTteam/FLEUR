MODULE m_fft3dxc
CONTAINS
   SUBROUTINE fft3dxc( &
      afft, bfft, fg3, &
      k1d, k2d, k3d, ng3, kimax, isn, &
      igfft1, igfft2, pgfft, nstr)

!***********************************************************
!                                                          *
! interface for fg3(star) -- fft --> (a,b)fft (r) (isn=+1) *
!         or (a,b)fft (r) -- fft --> fg3(star)    (isn=-1) *
!                                                          *
! dimension of (a,b)fft is (k1d x k2d x k3d)               *
! afft and bfft contain the real/imaginary part of the fft *
! igfft1(i)  is the pointer from the g-sphere to stars     *
! igfft2(i)  is the pointer from the g-sphere to fft-grid  *
! pgfft(i)   contains the phases of the g-vectors of sph.  *
!                                                          *
!***********************************************************

      USE m_fft_interface
      IMPLICIT NONE

      INTEGER :: k1d, k2d, k3d, ng3, kimax, isn
      INTEGER :: igfft1(0:k1d*k2d*k3d-1), igfft2(0:k1d*k2d*k3d-1)
      INTEGER :: nstr(ng3)
      COMPLEX pgfft(0:k1d*k2d*k3d-1)
      REAL ::    afft(0:k1d*k2d*k3d-1), bfft(0:k1d*k2d*k3d-1)
      COMPLEX :: fg3(ng3)

      INTEGER :: i, ifftd
      REAL :: scale, zero
      COMPLEX :: ctmp

      LOGICAL :: forw
      INTEGER :: length_zfft(3)
      complex :: zfft(0:k1d*k2d*k3d-1)

      ifftd = k1d*k2d*k3d
      zero = 0.0

      IF (isn > 0) THEN

         !  ---> put stars onto the fft-grid

         afft = 0.0
         bfft = 0.0
         DO i = 0, kimax - 1
            ctmp = fg3(igfft1(i))*pgfft(i)
            afft(igfft2(i)) = real(ctmp)
            bfft(igfft2(i)) = aimag(ctmp)
         ENDDO
      ENDIF

!---> now do the fft (isn=+1 : g -> r ; isn=-1 : r -> g)

      zfft = cmplx(afft, bfft)
      if (isn == -1) then
         forw = .TRUE.
      else
         forw = .FALSE.
      end if
      length_zfft(1) = k1d
      length_zfft(2) = k2d
      length_zfft(3) = k3d
      call fft_interface(3, length_zfft, zfft, forw, igfft2(0:kimax-1))
      afft = real(zfft)
      bfft = aimag(zfft)

      IF (isn < 0) THEN

         !  ---> collect stars from the fft-grid

         DO i = 1, ng3
            fg3(i) = cmplx(0.0, 0.0)
         ENDDO
         scale = 1.0/ifftd
         DO i = 0, kimax - 1
            fg3(igfft1(i)) = fg3(igfft1(i)) + CONJG(pgfft(i))* &
                             zfft(igfft2(i))
         ENDDO
         DO i = 1, ng3
            fg3(i) = scale*fg3(i)/nstr(i)
         ENDDO
      ENDIF

   END SUBROUTINE fft3dxc
END MODULE m_fft3dxc
