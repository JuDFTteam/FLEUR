MODULE m_prpxcfftmap
  USE m_juDFT
  !*********************************************************************
  !     this subroutine prepares the pointer which identifies a
  !     threedimensional g-vector in the positive domain of the
  !     xc (=charge density) fft-box in order to map a 3-d g-vector
  !     onto stars in case of the backtransform for fft of the
  !     charge density. correspondes  to igfft(*,2)
  !     Further it sets up the x,y, and z component of the 3-dimensional
  !     g-vector in the original domain of all g-vectors used for fft.
  !     it is independent of spin and k-point.
  !     pointer is built up when ever the chargedensity is calculated
  !     in order to save memory
  !
  !        s. bluegel, IFF, Aug. 97   
  !*********************************************************************
CONTAINS
  SUBROUTINE prp_xcfft_map(&
       &                         stars,sym,&
       &                         cell,&
       &                         igxc_fft,gxc_fft)
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_cell),INTENT(IN)  :: cell
    !
    !
    INTEGER,INTENT(OUT) ::    igxc_fft(0:stars%kxc1d*stars%kxc2d*stars%kxc3d-1)
    REAL   ,INTENT(OUT) ::    gxc_fft(0:stars%kxc1d*stars%kxc2d*stars%kxc3d-1,3)
    !
    !---> local variables
    !
    LOGICAL NEW
    INTEGER istr,iop,iopm1,il,im,in,kidx,iv1d,ifftq1,ifftq2
    INTEGER nop_local,norm,kr(3,sym%nop)

    !------->          abbreviations
    !
    !     kxc1d  : dimension of the charge density fft box in the pos. domain
    !     kxc2d  : defined in dimens.f program (subroutine apws).1,2,3 indic
    !     kxc3d  ; a_1, a_2, a_3 directions.
    !     kq(i) : i=1,2,3 actual length of the fft-box for which fft is done
    !     nstr  : number of members (arms) of reciprocal lattice (g) vector
    !             of each star.
    !     nxc3_fft: number of stars in the  charge density  fft-box
    !     kmxxc_fft: number of g-vectors forming the nxc3_fft stars in the
    !               charge density sphere
    !     gxc_fft : contains x,y,z components of g-vectors contributing to FFT.
    !
    !-----> prepare pointer which identifies a threedimensional g-vector
    !       in the positive domain of the charge density fft-box.
    !       correspondes  to igfft(*,2)
    !
    kidx    = 0
    ifftq1  = stars%kxc1_fft
    ifftq2  = stars%kxc1_fft*stars%kxc2_fft
    !
    DO istr = 1 , stars%nxc3_fft
       !
       nop_local=sym%nop
       IF (stars%kv3(3,istr).EQ.0) nop_local=sym%nop2
       !
       DO iop = 1,nop_local
          kr(1,iop) = stars%kv3(1,istr)*sym%mrot(1,1,iop)&
               &                + stars%kv3(2,istr)*sym%mrot(2,1,iop)&
               &                + stars%kv3(3,istr)*sym%mrot(3,1,iop)
          kr(2,iop) = stars%kv3(1,istr)*sym%mrot(1,2,iop)&
               &                + stars%kv3(2,istr)*sym%mrot(2,2,iop)&
               &                + stars%kv3(3,istr)*sym%mrot(3,2,iop)
          kr(3,iop) = stars%kv3(1,istr)*sym%mrot(1,3,iop)&
               &                + stars%kv3(2,istr)*sym%mrot(2,3,iop)&
               &                + stars%kv3(3,istr)*sym%mrot(3,3,iop)
       ENDDO
       !
       DO iop = 1 , nop_local
          NEW=.TRUE.
          DO iopm1 = 1 , iop - 1
             norm=(kr(1,iop)-kr(1,iopm1))**2 +&
                  &             (kr(2,iop)-kr(2,iopm1))**2 +&
                  &             (kr(3,iop)-kr(3,iopm1))**2
             IF(norm.EQ.0) NEW=.FALSE.
          ENDDO

          IF (NEW) THEN
             il=kr(1,iop)
             im=kr(2,iop)
             in=kr(3,iop)
             gxc_fft(kidx,1) = cell%bmat(1,1)*il+cell%bmat(2,1)*im+cell%bmat(3,1)*in 
             gxc_fft(kidx,2) = cell%bmat(1,2)*il+cell%bmat(2,2)*im+cell%bmat(3,2)*in
             gxc_fft(kidx,3) = cell%bmat(1,3)*il+cell%bmat(2,3)*im+cell%bmat(3,3)*in
             IF (il.LT.0) il=il+stars%kxc1_fft
             IF (im.LT.0) im=im+stars%kxc2_fft
             IF (in.LT.0) in=in+stars%kxc3_fft
             iv1d = in*ifftq2 + im*ifftq1 + il
             igxc_fft(kidx)=iv1d
             kidx=kidx+1
          ENDIF
       ENDDO
    ENDDO
    !
    IF (kidx .NE. stars%kmxxc_fft) THEN
       WRITE (6,'('' something wrong with stars%kmxxc_fft or nxc3_fft'')')
       WRITE (6,'('' stars%kmxxc_fft, acutal kidx '',2i5)')&
            &                stars%kmxxc_fft, kidx
       CALL juDFT_error("kidx /= stars",calledby ="prp_xcfft_map"&
            &        ,hint ="something wrong with kmxxc_fft or nxc3_fft")
    ENDIF
    !
  END SUBROUTINE prp_xcfft_map
END MODULE m_prpxcfftmap
