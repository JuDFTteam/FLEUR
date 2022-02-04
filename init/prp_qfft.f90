!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_prpqfft
      use m_juDFT
      CONTAINS
      SUBROUTINE prp_qfft(l_write,stars,cell,noco,input)
!*********************************************************************
!     This subroutine prepares the necessary variables and checks
!     to calculate the plane wave chargedensity in the interstitial
!     (subroutine pwden.f) by fast fourier transform (FFT).
!     In certain cases rkmax is slightly reajusted to perform
!     a quick FFT based on powers   (2**P) * (3**Q) * (5**R) only.
!
!     dimensions kd(i)d for charge denisity box are checked.
!        s. bluegel, JRCAT, Feb. 97
!
!     subroutine boxdim added to subroutine to treat non-orthogonal
!     lattice systems.
!        s.bluegel, IFF, 18.Nov.97
!*********************************************************************

      USE m_types
      USE m_constants
      USE m_ifft, ONLY : ifft235
      USE m_boxdim

      IMPLICIT NONE

!---> Arguments
      !
      LOGICAL,INTENT(IN)           :: l_write
      TYPE(t_stars),INTENT(INOUT)  :: stars
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_input),INTENT(INOUT)  :: input  !rkmax might be modified
!
!
!---> local variables
!
      INTEGER ksfft,mq1,mq2,mq3,istr,iofile,ng2_fft,kmxq2_fft
      REAL    arltv1,arltv2,arltv3,rknew
!
!---> intrinsic functions
!    
      INTRINSIC int,min
!
!---> data statement
!
      REAL gmaxp
      DATA gmaxp /2.0/
!
!------->          BACKGROUND       
!
!        Determine the limits  for the charge density fft-box
!        The charge density is the  square of the wavefunction.
!        Since the largest G-vector of the wavefunction is
!        given by |G + k| < rkmax, the largest G -vector gmaxq
!        contributing the charge-density is exactly :
!        gmaxq = gmaxp*rkmax,  with gmaxp =2. 
!        Therefore the box dimension must be  determined by :
!        m(i) >= gmaxp * rkmax*a(i)/(2*pi), where a(i) is the 
!        magnitude of the ith basis vector. (this is only true
!        for <a(i)|a(j)> = 0 for i .ne. j., otherwise see boxdim)
!        REMEMBER: gmax used through out the FLAPW code can be
!        larger than gmaxp*rkmax. 
!
!
!------->          ABBREVIATIONS
!
!   ksfft=(0,1) : KEY OF SELECTING FFT-PRDOGRAM AND RADIX-TYPE
!                      0  PROGRAM, RADIX-2 ONLY
!                      1  PROGRAM, RADIX-2, RADIX-3,RADIX-5
!   iofile      : device number for in and output
!   gmax        : actually used gmax troughtout the FLAPW code
!   gmaxq       : cut-off wavevector for charge density
!   rkmax       : cut-off for |g+k|
!   gmaxp       : gmaxp = gmaxq/rkmax, ideal: gmaxp=2
!   kq1,2,3d    : dimensions of the charge density fft-box
!                 ( in positive fft domain )
!   mq1,2,3     : actual size of the charge density fft-box
!                 ( in positive fft domain )
!   nq3_fft     : number of stars in the charge density sphere
!   kmxq_fft    : number of g-vectors in the charge density sphere
!   arltv(i)    : length of reciprical lattice vector along
!                 direction (i)
!
!---> Determine dimensions of fft-box of size mq1, mq2, mq3,
!     for which |G(mq1,mq2,mq3)| < Gmaxp*Kmax
!
      CALL boxdim(&
     &            cell%bmat,&
     &            arltv1,arltv2,arltv3)
!
      mq1 = int( gmaxp*input%rkmax/arltv1 ) + 1
      mq2 = int( gmaxp*input%rkmax/arltv2 ) + 1
      mq3 = int( gmaxp*input%rkmax/arltv3 ) + 1

!---> add + 1 in spin spiral calculation, to make sure that all G's are 
!---> still within the FFT-box after being shifted by the spin spiral
!---> q-vector.
      IF (noco%l_ss) THEN
         mq1 = mq1 + 1
         mq2 = mq2 + 1
         mq3 = mq3 + 1
      ENDIF         
!
!----> fft done in positive domain 
!      (Remember number of grid points not 2*m+1 because of cyclic 
!       boundary condition)
!
      mq1 = mq1 + mq1
      mq2 = mq2 + mq2
      mq3 = mq3 + mq3
!
!---> fft's are usually fastest for low primes
!     (restrict kqid to: kwid=  (2**P) * (3**Q) * (5**R)
!
      iofile = oUnit
      ksfft  = 1
      stars%kq1_fft = ifft235(ksfft,mq1,gmaxp)
      stars%kq2_fft = ifft235(ksfft,mq2,gmaxp)
      stars%kq3_fft = ifft235(ksfft,mq3,gmaxp)
!+gb
!      kq1_fft = mq1
!      kq2_fft = mq2
!      kq3_fft = mq3
!
!---> For mq = 2**P, FFT is very fast. If mq very close to 2**P
!     choose this, even 2**P < mq . Therefore:
!     factor 0.5 added by S.B. 21.Aug.97 
!
      rknew = min( real( stars%kq1_fft )*arltv1, real( stars%kq2_fft )*arltv2, &
     &             real( stars%kq3_fft )*arltv3)
      rknew = 0.5 * rknew / gmaxp
      IF (rknew.LT.input%rkmax) THEN
         IF (l_write) WRITE (oUnit,'('' rkmax and true gmax recalculated '')')
         IF (l_write) WRITE (oUnit,2100) input%rkmax, rknew, rknew*rknew
         IF (l_write) WRITE (oUnit,2200) gmaxp*rknew, gmaxp*rknew*gmaxp*rknew
         input%rkmax = rknew
      ENDIF
!
!-----> gmax => gmaxp*rkmax 
!       (otherwise too few elements in arrays defined in strng1)
!
      IF (2*input%rkmax.GT.stars%gmax) THEN
         WRITE (oUnit,'('' gmax must be at least 2*rkmax'')')
         WRITE (oUnit,'('' increase gmax , or reduce rkmax'')')
         WRITE (oUnit,'('' rkmax ='',f10.3,''  gmax ='',f10.3)') input%rkmax,stars%gmax
         CALL juDFT_error("rkmax,gmax",calledby ="prp_qfft",hint&
     &        ="gmax must be at least 2*rkmax")
      ENDIF
!
!     mq1 = kq1_fft
!     mq2 = kq2_fft
!     mq3 = kq3_fft
!

!
!-----> how many stars are in charge density sphere?
!       assume stars are ordered according to length
!
!---> 3d stars
      stars%ng3_fft = 0
      DO istr = 1 , stars%ng3
         IF ( stars%sk3(istr).LE.gmaxp*input%rkmax ) THEN
            stars%ng3_fft = istr
         ENDIF
      ENDDO
!---> 2d stars
      ng2_fft = 0
      DO istr = 1 , stars%ng2
         IF ( stars%sk2(istr).LE.gmaxp*input%rkmax ) THEN
            ng2_fft = istr
         ENDIF
      ENDDO
!
      IF ( stars%ng3_fft.EQ.0 ) THEN
         WRITE(oUnit,'('' presumably ng3 too small '')')
         WRITE(oUnit,'('' sk3max, gmaxp*rkmax '', 2f10.6)') &
     &                stars%sk3(stars%ng3),gmaxp*input%rkmax
         CALL juDFT_error("presumably ng3 too small","prp_qfft")
      ENDIF
!
      IF ( stars%ng3_fft.GT.stars%ng3 ) THEN
         WRITE(oUnit,'('' nq3_fft > n3d '')')
         WRITE(oUnit,'('' nq3_fft, n3d '',2i10)') stars%ng3_fft, stars%ng3
          CALL juDFT_error("nq3_fft > n3d",calledby="prp_qfft")
      ENDIF
!
!-----> check that all nq3_fft stars fit into the charge density FFT-box
!
      DO istr = 1 , stars%ng3_fft
        IF ( ( 2*stars%kv3(1,istr).GT.stars%kq1_fft ) .OR.&
     &       ( 2*stars%kv3(2,istr).GT.stars%kq2_fft ) .OR.&
     &       ( 2*stars%kv3(3,istr).GT.stars%kq3_fft ) ) THEN
          WRITE (oUnit,'('' not all nq3_fft stars in chg. den. FFT box'')')
          WRITE (oUnit,'('' inconsistency in def.s see also strgn1'')')
          WRITE (oUnit,'('' mq1d,mq2d,mq3d,kv1,kv2,kv3 '',6i5)')&
     &                  stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,2*stars%kv3(1,istr),2*stars%kv3(2,istr),&
     &                                 2*stars%kv3(3,istr)
          CALL juDFT_error("not all nq3_fft stars in chg. den. FFT box",&
     &         calledby ="prp_qfft")
        ENDIF
      ENDDO
!
!-----> how many G-vectors belong to these nq3_fft stars                    
!
!---> 3d vectors
      stars%kmxq_fft = 0
      DO istr = 1 , stars%ng3_fft
         stars%kmxq_fft = stars%kmxq_fft + stars%nstr(istr)
      ENDDO
      IF ( stars%kmxq_fft .GT. stars%kq1_fft*stars%kq2_fft*stars%kq3_fft ) THEN
         IF (l_write) WRITE (oUnit,'('' array dimensions in later subroutines too'',&
     &             '' small'',2i10)') stars%kmxq_fft,stars%kq1_fft*stars%kq2_fft*stars%kq3_fft
      ENDIF
!---> 2d vectors
      kmxq2_fft = 0
      DO istr = 1 , ng2_fft
         kmxq2_fft = kmxq2_fft + stars%nstr2(istr)
      ENDDO
      IF ( kmxq2_fft .GT. stars%kq1_fft*stars%kq2_fft ) THEN
         IF (l_write) WRITE (oUnit,'('' array dimensions in later subroutines too'',&
     &             '' small'',2i10)') kmxq2_fft,stars%kq1_fft*stars%kq2_fft
      ENDIF

!
 2100 FORMAT (/,1x,'old rkmax   =',f10.5, '(a.u.)**(-1) ==>  new rkmax '&
     &      ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new E_cut(wf)   =',&
     &            f10.5,' Ry')
 2200 FORMAT (/,1x,'true gmax   =',f10.5, '(a.u.)**(-1)',&
     &       /,t38,'==>  new E_cut(chg)  =', f10.5,' Ry')

      END SUBROUTINE prp_qfft
      END MODULE m_prpqfft
