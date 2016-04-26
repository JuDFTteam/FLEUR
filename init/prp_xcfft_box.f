      SUBROUTINE prp_xcfft_box(gmaxxc,bmat,kxc1d,kxc2d,kxc3d)
c*********************************************************************
c     Determine the dimensions kxc1d, kxc2d, kxc3d
c     of the dimension of the xc-energy and potential fft-box
c     needed for the fast calculation in interstitial region
c
c     dimensions kxc(i)d for xc-potential/energy FFT in interstitial
c     added.
c        s. bluegel, IFF, Oct. 97
c
c     subroutine boxdim added to subroutine prp_xcfft
c        s.bluegel, IFF, 18.Nov.97
c*********************************************************************
c
      USE m_ifft, ONLY : ifft235
      USE m_boxdim
      IMPLICIT NONE

      INTEGER kxc1d,kxc2d,kxc3d
      REAL bmat(3,3),gmaxxc

      REAL arltv1,arltv2,arltv3
      INTEGER mk1,mk2,mk3,iofile,ksfft
C ..

C ... Intrinsic Functions
      INTRINSIC int
C
C------->          ABBREVIATIONS
C
C   ksfft=(0,1) : KEY OF SELECTING FFT-PRDOGRAM AND RADIX-TYPE
C                      0  PROGRAM, RADIX-2 ONLY
C                      1  PROGRAM, RADIX-2, RADIX-3,RADIX-5
c   iofile      : device number for in and output
c   gmax        : cut-off wavevector for charge density
c   gmaxxc      : cut-off wavevector for xc-potential and energy
c                 in interstitial
c   arltv(i)    : length of reciprical lattice vector along
c                 direction (i)
c
C---> Determine rkmax box of size mk1, mk2, mk3,
c     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < gmaxxc
C
      CALL boxdim(bmat,arltv1,arltv2,arltv3)

c     (add 1+1 due to integer rounding, strange k_vector in BZ)
      mk1 = int(gmaxxc/arltv1) + 2
      mk2 = int(gmaxxc/arltv2) + 2
      mk3 = int(gmaxxc/arltv3) + 2
c
      kxc1d = 2*mk1
      kxc2d = 2*mk2
      kxc3d = 2*mk3
c
c---> fft's are usually fastest for low primes
c     (restrict kxcid to: kxcid=  (2**P) * (3**Q) * (5**R)
c
      iofile = 6
      ksfft = 1
      WRITE (6,*) 'minimum: kxc1d,kxc2d,kxc3d',kxc1d,kxc2d,kxc3d
      kxc1d = ifft235(iofile,ksfft,kxc1d,2.0)
      kxc2d = ifft235(iofile,ksfft,kxc2d,2.0)
      kxc3d = ifft235(iofile,ksfft,kxc3d,2.0)
      WRITE (6,*) 'ifft235: kxc1d,kxc2d,kxc3d',kxc1d,kxc2d,kxc3d

      RETURN
      END
