!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_prpxcfft
   use m_juDFT
CONTAINS
   SUBROUTINE prp_xcfft(mpi, stars, input, cell, xcpot)
!*********************************************************************
!     this subroutine prepares the necessary variables and checks
!     to calculate the xc-potential and energy in the interstitial
!     (subroutine visxc(g).f) by fast fourier transform (fft).
!     in certain cases gmaxxc is slightly reajusted to perform
!     a quick fft based on powers   (2**p) * (3**q) * (5**r) only.
!
!     dimensions kd(i)d for charge denisity box are checked.
!        s. bluegel, iff , oct. 97
!     subroutine boxdim added to subroutine to treat non-othogonal
!     lattice systems
!        s.bluegel, IFF, 18.Nov.97
!*********************************************************************

      USE m_types
      USE m_constants
      USE m_ifft, ONLY: ifft235
      USE m_boxdim

      IMPLICIT NONE

      type(t_mpi), intent(in)       :: mpi
      TYPE(t_stars), INTENT(INOUT)   :: stars
      TYPE(t_input), INTENT(IN)      :: input
      TYPE(t_cell), INTENT(IN)       :: cell
      CLASS(t_xcpot), INTENT(INOUT)  :: xcpot

!
!---> local variables
!
      INTEGER ksfft, mxc1, mxc2, mxc3, istr, iofile
      REAL arltv1, arltv2, arltv3, gmxxc_new
!
!---> intrinsic functions
!
      INTRINSIC int, min
!
!------->          background
!
!        determine the limits  for the xc energy/pot fft-box.
!        since the xc functional has a non-linear dependence
!        of the charge density in many cases a gmax to determine
!        the xc-potential or energy of gmax=2*rkmax=gmaxq is
!        in many cases not sufficent. Therefore, we use a larger
!        gmax called gmaxxc to evaluate the xc-potential and energy
!        in the interstitial region. This leads to a larger xc-fft-box.
!        the box dimension is  determined by :
!        m(i) >= gmaxxc*a(i)/(2*pi), where a(i) is the
!        magnitude of the ith basis vector. (this is only true
!        for <a(i)|a(j)> = 0 for i .ne. j., otherwise see boxdim)
!        remember: gmax used through out the flapw code must be
!        larger than or equal  gmax >= gmaxxc >= gmaxp*rkmax.
!
!
!------->          abbreviations
!
!   ksfft=(0,1) : key of selecting fft-prdogram and radix-type
!                      0  program, radix-2 only
!                      1  program, radix-2, radix-3,radix-5
!   iofile      : device number for in and output
!   rkmax       : cut-off for |g+k|
!   gmax        : actually used largest g-vector used througout
!               : the calculation
!   gmaxxc      : cut-off wavevector for fast fourier transform of
!                 xc-potential and energy
!                 gmaxxc usually gmaxxc >= 2*rkmax
!   kxc1,2,3d   : dimensions of the xc-pot/energy fft-box
!                 ( in positive fft domain )
!   kxc1,2,3_fft: actual size of the xc-pot/energy fft-box
!                 ( in positive fft domain )
!   nxc3_fft     : number of stars in the xc-pot/energy sphere
!   kmxxc_fft    : number of g-vectors in the xc-pot/energy sphere
!   arltv(i)    : length of reciprical lattice vector along
!                 direction (i)
!
      IF(.NOT. xcpot%needs_grad()) xcpot%gmaxxc = stars%gmax
      if(mpi%irank == 0) WRITE(oUnit, '('' gmaxxc should be: 2*kmax <= gmaxxc <= gmax '')')
      IF(abs(xcpot%gmaxxc - stars%gmax) .le. 10.0**(-6)) THEN
         if(mpi%irank == 0) WRITE(oUnit, '('' concerning memory, you may want to choose'',&
        &              '' a smaller value for gmax'')')
      END IF
      IF(xcpot%gmaxxc .LE. 10.0**(-6)) THEN
         if(mpi%irank == 0) WRITE(oUnit, '(" gmaxxc=0 : gmaxxc=gmax choosen as default",&
     &              " value")')
         if(mpi%irank == 0) WRITE(oUnit, '(" concerning memory, you may want to choose",&
     &              " a smaller value for gmax")')
         xcpot%gmaxxc = stars%gmax
      END IF
      IF(xcpot%gmaxxc .LE. 2*input%rkmax) THEN
         if(mpi%irank == 0) WRITE(oUnit, '('' concerning accuracy and total energy'',&
     &              '' convergence, you may want'',/,&
     &              '' to choose a larger gmaxxc '')')
      END IF
      if(mpi%irank == 0) write(oUnit, '('' gmaxxc ='',f10.6)') xcpot%gmaxxc
!
!---> Determine dimensions of fft-box of size mxc1, mxc2, mxc3,
!     for which |G(mxc1,mxc2,mxc3)| < Gmaxxc
!
      CALL boxdim(&
     &            cell%bmat,&
     &            arltv1, arltv2, arltv3)!keep
!
      mxc1 = int(xcpot%gmaxxc/arltv1) + 1
      mxc2 = int(xcpot%gmaxxc/arltv2) + 1
      mxc3 = int(xcpot%gmaxxc/arltv3) + 1
!
!----> fft done in positive domain
!      (remember number of grid points not 2*m+1 because of cyclic
!       boundary condition)
!
      mxc1 = mxc1 + mxc1
      mxc2 = mxc2 + mxc2
      mxc3 = mxc3 + mxc3
!
!---> fft's are usually fastest for low primes
!     (restrict kqid to: kwid=  (2**p) * (3**q) * (5**r)
!
      iofile = oUnit
      ksfft = 1
      stars%kxc1_fft = ifft235(ksfft, mxc1, 2.0)
      stars%kxc2_fft = ifft235(ksfft, mxc2, 2.0)
      stars%kxc3_fft = ifft235(ksfft, mxc3, 2.0)
!
!---> for mxc = 2**p, fft is very fast. if mq very close to 2**p
!     choose this, even 2**p < mxc . therefore:
!
      gmxxc_new = min(real(stars%kxc1_fft)*arltv1, &
     &                 real(stars%kxc2_fft)*arltv2, real(stars%kxc3_fft)*arltv3)
      gmxxc_new = 0.5*gmxxc_new

      IF(gmxxc_new .LT. xcpot%gmaxxc) THEN
         WRITE(oUnit, '('' gmaxxc recalculated '')')
         WRITE(oUnit, 2100) xcpot%gmaxxc, gmxxc_new, gmxxc_new*gmxxc_new
         xcpot%gmaxxc = gmxxc_new
      ENDIF
!
!-----> gmax => gmaxxc
!       (otherwise too few elements in arrays defined in strng1)
!
      IF(gmxxc_new .GT. stars%gmax) THEN
         if(mpi%irank == 0) THEN
            WRITE(oUnit, '('' gmax must be at least gmxxc_new'')')
            WRITE(oUnit, '('' increase gmax , or reduce gmaxxc'')')
            WRITE(oUnit, '('' gmxxc_new ='',f10.3,''  gmax ='',f10.3)') &
        &                gmxxc_new, stars%gmax
!cc          CALL juDFT_error("gmxxc_new.gt.gmax",calledby="prp_xcfft")
         endif
      ENDIF
!
!------> check dimensions of fft chargedensity box used in pwden.f
!
      IF(stars%kxc1_fft .GT. stars%kxc1_fft .OR. stars%kxc2_fft .gt. stars%kxc2_fft .OR. &
    &                            stars%kxc3_fft .gt. stars%kxc3_fft) THEN
         if(mpi%irank == 0) THEN
            WRITE(oUnit, '('' box dim. for fft too small'')')
            WRITE(oUnit, 2110) stars%kxc1_fft, stars%kxc2_fft, stars%kxc3_fft, stars%kxc1_fft, stars%kxc2_fft, stars%kxc3_fft
            CALL juDFT_error("mxc[1,2,3]d>kxc[1,2,3]d ", calledby&
       &         ="prp_xcfft")
         endif
      ENDIF
2110  FORMAT(' stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft ', 6i5)
!
!-----> how many stars are in charge density sphere?
!       assume stars are ordered according to length
!
      stars%nxc3_fft = 0
      DO istr = 1, stars%ng3
         IF(stars%sk3(istr) .LE. xcpot%gmaxxc) THEN
            stars%nxc3_fft = istr
         ENDIf
      ENDDO
!
      IF(stars%nxc3_fft .EQ. 0) THEN
         if(mpi%irank == 0) THEN
            WRITE(oUnit, '('' presumably ng3 too small '')')
            WRITE(oUnit, '('' sk3max, gmaxxc '', 2f10.6)')&
        &                stars%sk3(stars%ng3), xcpot%gmaxxc
            CALL juDFT_error("nxc3_fft==0", calledby="prp_xcfft")
         endif
      ENDIF
!
      IF(stars%nxc3_fft .GT. stars%ng3) THEN
         if(mpi%irank == 0) THEN
            WRITE(oUnit, '('' nxc3_fft > n3d '')')
            WRITE(oUnit, '('' nxc3_fft, n3d '',2i10)') stars%nxc3_fft, stars%ng3
         endif
         CALL juDFT_error("nxc3_fft>n3d ", calledby="prp_xcfft")
      ENDIF
!
!-----> check that all nxc3_fft stars fit into the xc-pot/energy fft-box
!
      DO istr = 1, stars%nxc3_fft
         IF((2*stars%kv3(1, istr) .gt. stars%kxc1_fft) .OR.&
      &       (2*stars%kv3(2, istr) .gt. stars%kxc2_fft) .OR.&
      &       (2*stars%kv3(3, istr) .gt. stars%kxc3_fft)) THEN
            if(mpi%irank == 0) THEN
               WRITE(oUnit, '('' not all nxc3_fft stars in xc-pot/eng fft box'')')
               WRITE(oUnit, '('' inconsistency in def.s see also strgn1'')')
               WRITE(oUnit, '('' kxc1_fft,kxc2_fft,kxc3_fft,kv1,kv2,kv3 '',6i5)')&
          &                 stars%kxc1_fft, stars%kxc2_fft, stars%kxc3_fft, 2*stars%kv3(1, istr),&
          &                 2*stars%kv3(2, istr), 2*stars%kv3(3, istr)
            endif
            CALL juDFT_error("2*stars([1,2,3],istr)>mxc[1,2,3]d", calledby&
       &         ="prp_xcfft")
         ENDIF
      ENDDO
!
!-----> how many g-vectors belong to these nxc3_fft stars
!
      stars%kmxxc_fft = 0
      DO istr = 1, stars%nxc3_fft
         stars%kmxxc_fft = stars%kmxxc_fft + stars%nstr(istr)
      ENDDO

      IF(stars%kmxxc_fft .gt. stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft) then
         if(mpi%irank == 0) THEN
            WRITE(oUnit, '('' array dimensions in later subroutines too'',&
        &             '' small'')')
         endif
      ENDIF

2100  format(/, 1x, 'old gmaxxc  =', f10.5, '(a.u.)**(-1) ==>  new gmaxxc'&
      &      , '  =', f10.5, '(a.u.)**(-1) ', /, t38, '==>  new e_cut(xc)   =',&
      &            f10.5, ' ry')

   END SUBROUTINE prp_xcfft
END MODULE m_prpxcfft
