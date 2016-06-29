!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_prpqfft
      use m_juDFT
      CONTAINS
      SUBROUTINE prp_qfft(&
     &                    stars,&
     &                    cell,noco,&
     &                    input)
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
!
      USE m_ifft, ONLY : ifft235
      USE m_boxdim
      USE m_types
      IMPLICIT NONE
!
!---> Arguments
!
      TYPE(t_stars),INTENT(INOUT)  :: stars
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_input),INTENT(INOUT)  :: input  !<rkmax might be modified
!
!
!---> local variables
!
      INTEGER ksfft,mq1,mq2,mq3,istr,iofile
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
      iofile = 6
      ksfft  = 1
      stars%kq1_fft = ifft235 (iofile,ksfft,mq1,gmaxp)
      stars%kq2_fft = ifft235 (iofile,ksfft,mq2,gmaxp)
      stars%kq3_fft = ifft235 (iofile,ksfft,mq3,gmaxp)
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
         WRITE (6,'('' rkmax and true gmax recalculated '')')
         WRITE (6,2100) input%rkmax, rknew, rknew*rknew
         WRITE (6,2200) gmaxp*rknew, gmaxp*rknew*gmaxp*rknew
         WRITE (16,'('' rkmax and true gmax recalculated '')')
         WRITE (16,2100) input%rkmax, rknew, rknew*rknew
         WRITE (16,2200) gmaxp*rknew, gmaxp*rknew*gmaxp*rknew
         input%rkmax = rknew
      ENDIF
!
!-----> gmax => gmaxp*rkmax 
!       (otherwise too few elements in arrays defined in strng1)
!
      IF (2*input%rkmax.GT.stars%gmax) THEN
         WRITE (6,'('' gmax must be at least 2*rkmax'')')
         WRITE (6,'('' increase gmax , or reduce rkmax'')')
         WRITE (6,'('' rkmax ='',f10.3,''  gmax ='',f10.3)') input%rkmax,stars%gmax
         CALL juDFT_error("rkmax,gmax",calledby ="prp_qfft",hint&
     &        ="gmax must be at least 2*rkmax")
      ENDIF
!
!     mq1 = kq1_fft
!     mq2 = kq2_fft
!     mq3 = kq3_fft
!
!------> check dimensions of FFT chargedensity box used in pwden.f
!
       IF ( stars%kq1_fft.gt.stars%kq1d .OR. stars%kq2_fft.gt.stars%kq2d .OR. stars%kq3_fft.gt.stars%kq3d) THEN
          WRITE ( 6,'('' box dim. for FFT too small'')')
          WRITE ( 6,'('' kq1_fft,kq2_fft,kq3_fft,kq1d,kq2d,kq3d '',6i5)')&
     &                   stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,stars%kq1d,stars%kq2d,stars%kq3d
          WRITE (16,'('' box dim. for FFT too small'')')
          WRITE (16,'('' mq1d,mq2d,mq3d,kq1d,kq2d,kq3d '',6i5)')&
     &                 stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,stars%kq1d,stars%kq2d,stars%kq3d
          CALL juDFT_error("box dim. for FFT too small",calledby&
     &         ="prp_qfft")
       ENDIF
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
      stars%ng2_fft = 0
      DO istr = 1 , stars%ng2
         IF ( stars%sk2(istr).LE.gmaxp*input%rkmax ) THEN
            stars%ng2_fft = istr
         ENDIF
      ENDDO
!
      IF ( stars%ng3_fft.EQ.0 ) THEN
         WRITE(6,'('' presumably ng3 too small '')')
         WRITE(6,'('' sk3max, gmaxp*rkmax '', 2f10.6)') &
     &                stars%sk3(stars%ng3),gmaxp*input%rkmax
         WRITE(16,'('' presumably ng3 too small '')')
         WRITE(16,'('' sk3max, gmaxp*rkmax '', 2f10.6)')&
     &                stars%sk3(stars%ng3),gmaxp*input%rkmax
         CALL juDFT_error("presumably ng3 too small","prp_qfft")
      ENDIF
!
      IF ( stars%ng3_fft.GT.stars%n3d ) THEN
         WRITE(6,'('' nq3_fft > n3d '')')
         WRITE(6,'('' nq3_fft, n3d '',2i10)') stars%ng3_fft, stars%n3d
         WRITE(16,'('' nq3_fft > n3d '')')
         WRITE(16,'('' nq3_fft, n3d '',2i10)') stars%ng3_fft, stars%n3d
          CALL juDFT_error("nq3_fft > n3d",calledby="prp_qfft")
      ENDIF
!
!-----> check that all nq3_fft stars fit into the charge density FFT-box
!
      DO istr = 1 , stars%ng3_fft
        IF ( ( 2*stars%kv3(1,istr).GT.stars%kq1_fft ) .OR.&
     &       ( 2*stars%kv3(2,istr).GT.stars%kq2_fft ) .OR.&
     &       ( 2*stars%kv3(3,istr).GT.stars%kq3_fft ) ) THEN
          WRITE (6,'('' not all nq3_fft stars in chg. den. FFT box'')')
          WRITE (6,'('' inconsistency in def.s see also strgn1'')')
          WRITE (6,'('' mq1d,mq2d,mq3d,kv1,kv2,kv3 '',6i5)')&
     &                  stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,2*stars%kv3(1,istr),2*stars%kv3(2,istr),&
     &                                 2*stars%kv3(3,istr)
          WRITE (16,'('' not all nq3_fft stars in chg. den. FFT box'')')
          WRITE (16,'('' inconsistency in def.s see also strgn1'')')
          WRITE (16,'('' mq1d,mq2d,mq3d,kv1,kv2,kv3 '',6i5)')&
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
      IF ( stars%kmxq_fft .GT. stars%kq1d*stars%kq2d*stars%kq3d ) THEN
         WRITE (6,'('' array dimensions in later subroutines too'',&
     &             '' small'',2i10)') stars%kmxq_fft,stars%kq1d*stars%kq2d*stars%kq3d
      ENDIF
!---> 2d vectors
      stars%kmxq2_fft = 0
      DO istr = 1 , stars%ng2_fft
         stars%kmxq2_fft = stars%kmxq2_fft + stars%nstr2(istr)
      ENDDO
      IF ( stars%kmxq2_fft .GT. stars%kq1d*stars%kq2d ) THEN
         WRITE (6,'('' array dimensions in later subroutines too'',&
     &             '' small'',2i10)') stars%kmxq2_fft,stars%kq1d*stars%kq2d
      ENDIF

!
 2100 FORMAT (/,1x,'old rkmax   =',f10.5, '(a.u.)**(-1) ==>  new rkmax '&
     &      ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new E_cut(wf)   =',&
     &            f10.5,' Ry')
 2200 FORMAT (/,1x,'true gmax   =',f10.5, '(a.u.)**(-1)',&
     &       /,t38,'==>  new E_cut(chg)  =', f10.5,' Ry')

      END SUBROUTINE prp_qfft
      END MODULE m_prpqfft
