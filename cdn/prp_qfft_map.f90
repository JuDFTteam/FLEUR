!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_prpqfftmap
  use m_juDFT
CONTAINS
  SUBROUTINE prp_qfft_map(stars,sym,input)
    !*********************************************************************
    !     This subroutine prepares the pointer which identifies a 
    !     threedimensional g-vector in the positive domain of the 
    !     charge density fft-box in order to map a 3-d g-vector
    !     onto stars in case of the backtransform for fft of the 
    !     charge density. correspondes  to igfft(*,2)     
    !     it is independent of spin and k-point. 
    !     pointer is built up when ever the chargedensity is calculated
    !     in order to save memory
    !
    !        s. bluegel, JRCAT, Feb. 97
    !*********************************************************************
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(INOUT):: stars
    !
    !
    !---> local variables
    !
    LOGICAL new
    INTEGER istr,iop,iopm1,il,im,in,kid2x,kidx,iv1d,ifftq1,ifftq2
    INTEGER norm,kr(3,sym%nop),nop_local

    !------->          ABBREVIATIONS
    !
    !     kq1d  : dimension of the charge density FFT box in the pos. domain
    !     kq2d  : defined in dimens.f program (subroutine apws).1,2,3 indicate
    !     kq3d  ; a_1, a_2, a_3 directions.
    !     kq(i) : i=1,2,3 actual length of the fft-box for which FFT is done.
    !     nstr  : number of members (arms) of reciprocal lattice (g) vector
    !             of each star.
    !     nq3_fft: number of stars in the  charge density  FFT-box
    !     kmxq_fft: number of g-vectors forming the nq3_fft stars in the
    !               charge density sphere
    !
    !-----> prepare pointer which identifies a threedimensional g-vector
    !       in the positive domain of the charge density fft-box.
    !       correspondes  to igfft(*,2)     
    !
    kidx    = 0
    kid2x   = 0
    ifftq1  = stars%kq1_fft
    ifftq2  = stars%kq1_fft*stars%kq2_fft
    !
    DO istr = 1 , stars%ng3_fft
       !
       nop_local=sym%nop
       IF (stars%kv3(3,istr).EQ.0) nop_local=sym%nop2
       !
       DO iop = 1,nop_local
          kr(1,iop) = stars%kv3(1,istr)*sym%mrot(1,1,iop) &
               + stars%kv3(2,istr)*sym%mrot(2,1,iop)+ stars%kv3(3,istr)*sym%mrot(3,1,iop)
          kr(2,iop) = stars%kv3(1,istr)*sym%mrot(1,2,iop) &
               + stars%kv3(2,istr)*sym%mrot(2,2,iop)+ stars%kv3(3,istr)*sym%mrot(3,2,iop)
          kr(3,iop) = stars%kv3(1,istr)*sym%mrot(1,3,iop) &
               + stars%kv3(2,istr)*sym%mrot(2,3,iop) + stars%kv3(3,istr)*sym%mrot(3,3,iop)
       ENDDO
       !
       DO iop = 1 , nop_local
          new=.true.
          DO iopm1 = 1 , iop - 1
             norm=(kr(1,iop)-kr(1,iopm1))**2 +&
                  (kr(2,iop)-kr(2,iopm1))**2 +(kr(3,iop)-kr(3,iopm1))**2
             IF (norm.EQ.0) new=.false.
          ENDDO

          IF (new) THEN
             il=kr(1,iop)
             im=kr(2,iop)
             in=kr(3,iop)
             if(il.lt.0) il=il+stars%kq1_fft
             if(im.lt.0) im=im+stars%kq2_fft
             if(in.lt.0) in=in+stars%kq3_fft
             iv1d = in*ifftq2 + im*ifftq1 + il
             stars%igq_fft(kidx)=iv1d 
             kidx=kidx+1
             IF (input%film.AND.(stars%kv3(3,istr).EQ.0)) THEN
                iv1d = im*ifftq1 + il
                stars%igq2_fft(kid2x)=iv1d 
                kid2x=kid2x+1
             ENDIF
          ENDIF
       ENDDO
       !
    ENDDO
    !
    IF (kidx .NE. stars%kmxq_fft) THEN
       WRITE (6,'('' something wrong with stars%kmxq_fft or nq3_fft'')')
       WRITE (6,'('' stars%kmxq_fft, acutal kidx '',2i5)') &
            &                stars%kmxq_fft, kidx
       CALL juDFT_error("something wrong with stars or nq3_fft"&
            &        ,calledby ="prp_qfft_map")
    ENDIF

  END SUBROUTINE prp_qfft_map
END MODULE m_prpqfftmap
