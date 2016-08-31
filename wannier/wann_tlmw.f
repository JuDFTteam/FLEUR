!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_tlmw
      use m_juDFT
      contains
      subroutine wann_tlmw(
     >               nwfs,lwf,mrwf,
     >               l_nocosoc,spi,jspin,
     <               tlmwf)
c*********************************************************************
c     Express the angular part of the trial orbital
c     in terms of spherical harmonics.
c     Y. Mokrousov
c*********************************************************************
c     Modifications for soc and noco.
c     Tidied-up version. 
c     Frank Freimuth
c*********************************************************************

      implicit none

      integer, intent (in)  :: nwfs
      integer, intent (in)  :: lwf(nwfs),mrwf(nwfs)
      logical, intent (in)  :: l_nocosoc
      integer, intent (in)  :: spi(nwfs),jspin
      complex, intent (out) :: tlmwf(0:3,-3:3,nwfs)

      integer   :: nwf,l,lr,mr,m,i,ii,j,jj
      complex   :: tlm(0:3,-3:3,1:7)
      complex   :: ic
      intrinsic :: cmplx,sqrt,cos,sin

      ic = cmplx(0.,1.)

c..in this part the 'primary' matrices for the basic non-hybridized
c..states are constructed, the coefficients are written to the tlm 
c..these coefficients arise due to the fact that the tws comprise the
c..real harmonics, while the wave functions operate in terms of normal 
c..ones

      tlm(:,:,:) = cmplx(0.,0.)

      do l = 0,3

         if (l.eq.0) then
c..s-state
           tlm(0,0,1) = cmplx(1.,0.)

         elseif (l.eq.1) then
c..pz-state
           tlm(1,0,1) = cmplx(1.,0.)
c..px-state
           tlm(1,-1,2) =  cmplx(1.,0.)/(sqrt(2.))
           tlm(1,1,2)  = -cmplx(1.,0.)/(sqrt(2.))
c..py-state
           tlm(1,-1,3) = ic/(sqrt(2.))              
           tlm(1,1,3)  = ic/(sqrt(2.))

         elseif (l.eq.2) then
c..dz2-state
           tlm(2,0,1) = cmplx(1.,0.)
c..dxz-state
           tlm(2,1,2)  = -cmplx(1.,0.)/(sqrt(2.))
           tlm(2,-1,2) = cmplx(1.,0.)/(sqrt(2.))
c..dyz-state
           tlm(2,1,3)  = ic/(sqrt(2.))
           tlm(2,-1,3) = ic/(sqrt(2.))
c..dx2-y2-state
           tlm(2,2,4)  =  cmplx(1.,0.)/(sqrt(2.))
           tlm(2,-2,4) =  cmplx(1.,0.)/(sqrt(2.))
c..dxy-state
           tlm(2,2,5)  = -ic/(sqrt(2.))
           tlm(2,-2,5) =  ic/(sqrt(2.))

         elseif (l.eq.3)then
c..
           tlm(3,0,1) = cmplx(1.,0.)
c..
           tlm(3,1,2)  = -cmplx(1.,0.)/(sqrt(2.))
           tlm(3,-1,2) = cmplx(1.,0.)/(sqrt(2.))
c..
           tlm(3,1,3)  = ic/(sqrt(2.))
           tlm(3,-1,3) = ic/(sqrt(2.))
c..
           tlm(3,2,4)  =  cmplx(1.,0.)/(sqrt(2.))
           tlm(3,-2,4) =  cmplx(1.,0.)/(sqrt(2.))
c..
           tlm(3,2,5)  = -ic/(sqrt(2.))
           tlm(3,-2,5) =  ic/(sqrt(2.))
c..
           tlm(3,3,6)  = -cmplx(1.,0.)/(sqrt(2.))
           tlm(3,-3,6) =  cmplx(1.,0.)/(sqrt(2.))
c..
           tlm(3,3,7)  =  ic/(sqrt(2.))
           tlm(3,-3,7) =  ic/(sqrt(2.))
         else
            CALL juDFT_error("no tlm for this l",calledby ="wann_tlmw")
         endif

      enddo

c..now we are ready for more complex hybridized states, which
c..correspond to the negative values of lwf

      tlmwf(:,:,:) = cmplx(0.,0.)

      do nwf = 1,nwfs
         
         if( ((3-2*jspin).ne.spi(nwf)).and.l_nocosoc ) cycle

         lr = lwf(nwf)
         mr = mrwf(nwf)

         if (lr.ge.0) then
               tlmwf(lr,:,nwf) = tlm(lr,:,mr)  
         elseif (lr.eq.-1) then

            if (mr.eq.1) then
c..sp-1
              tlmwf(0,:,nwf) =  (1./(sqrt(2.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) =  (1./(sqrt(2.)))*tlm(1,:,2)
            elseif (mr.eq.2) then
c..sp-2
              tlmwf(0,:,nwf) =  (1./(sqrt(2.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(2.)))*tlm(1,:,2)
            endif 

         elseif (lr.eq.-2) then
              
            if (mr.eq.1) then
c..sp2-1
              tlmwf(0,:,nwf) =  (1./(sqrt(3.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(6.)))*tlm(1,:,2) +
     +                          (1./(sqrt(2.)))*tlm(1,:,3)
            elseif (mr.eq.2) then
c..sp2-2
              tlmwf(0,:,nwf) =  (1./(sqrt(3.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(6.)))*tlm(1,:,2)
     +                         -(1./(sqrt(2.)))*tlm(1,:,3)
            elseif (mr.eq.3) then
c..sp2-3
              tlmwf(0,:,nwf) =  (1./(sqrt(3.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) =  (2./(sqrt(6.)))*tlm(1,:,2)
            endif

         elseif (lr.eq.-3) then

            if (mr.eq.1) then
c..sp3-1
              tlmwf(0,:,nwf) = 0.5*tlm(0,:,1) 
              tlmwf(1,:,nwf) = 0.5*( tlm(1,:,2)+tlm(1,:,3)+tlm(1,:,1) ) 
            elseif (mr.eq.2) then
c..sp3-2
              tlmwf(0,:,nwf) = 0.5*tlm(0,:,1) 
              tlmwf(1,:,nwf) = 0.5*( tlm(1,:,2)-tlm(1,:,3)-tlm(1,:,1) ) 
            elseif (mr.eq.3) then
c..sp3-4
              tlmwf(0,:,nwf) = 0.5*tlm(0,:,1) 
              tlmwf(1,:,nwf) = 0.5*(-tlm(1,:,2)+tlm(1,:,3)-tlm(1,:,1) ) 
            elseif (mr.eq.4) then
c..sp3-4
              tlmwf(0,:,nwf) = 0.5*tlm(0,:,1) 
              tlmwf(1,:,nwf) = 0.5*(-tlm(1,:,2)-tlm(1,:,3)+tlm(1,:,1) ) 
            endif  

         elseif (lr.eq.-4) then

            if (mr.eq.1) then
c..sp3d-1
              tlmwf(0,:,nwf) =  (1./(sqrt(3.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(6.)))*tlm(1,:,2) 
     &                         +(1./(sqrt(2.)))*tlm(1,:,3)
            elseif (mr.eq.2) then
c..sp3d-2
              tlmwf(0,:,nwf) =  (1./(sqrt(3.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(6.)))*tlm(1,:,2)
     &                         -(1./(sqrt(2.)))*tlm(1,:,3)            
            elseif (mr.eq.3) then 
c..sp3d-3
              tlmwf(0,:,nwf) =  (1./(sqrt(3.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) =  (2./(sqrt(6.)))*tlm(1,:,2)
            elseif (mr.eq.4) then
c..sp3d-4
              tlmwf(1,:,nwf) =  (1./(sqrt(2.)))*tlm(1,:,1)
              tlmwf(2,:,nwf) =  (1./(sqrt(2.)))*tlm(2,:,1)
            elseif (mr.eq.5) then
c..sp3d-5
              tlmwf(1,:,nwf) = -(1./(sqrt(2.)))*tlm(1,:,1)
              tlmwf(2,:,nwf) =  (1./(sqrt(2.)))*tlm(2,:,1)
            endif 

         elseif (lr.eq.-5) then

            if (mr.eq.1) then
c..sp3d2-1
              tlmwf(0,:,nwf) =  (1./(sqrt(6.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(2.)))*tlm(1,:,2)
              tlmwf(2,:,nwf) = -(1./(sqrt(12.)))*tlm(2,:,1)
     &                         +0.5*tlm(2,:,4) 
            elseif (mr.eq.2) then
c..sp3d2-2
              tlmwf(0,:,nwf) =  (1./(sqrt(6.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) =  (1./(sqrt(2.)))*tlm(1,:,2)
              tlmwf(2,:,nwf) = -(1./(sqrt(12.)))*tlm(2,:,1)
     &                         +0.5*tlm(2,:,4) 
            elseif (mr.eq.3) then
c..sp3d2-3
              tlmwf(0,:,nwf) =  (1./(sqrt(6.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(2.)))*tlm(1,:,3)
              tlmwf(2,:,nwf) = -(1./(sqrt(12.)))*tlm(2,:,1)
     &                         -0.5*tlm(2,:,4) 
            elseif (mr.eq.4) then
c..sp3d2-4
              tlmwf(0,:,nwf) =  (1./(sqrt(6.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) =  (1./(sqrt(2.)))*tlm(1,:,3)
              tlmwf(2,:,nwf) = -(1./(sqrt(12.)))*tlm(2,:,1)
     &                         -0.5*tlm(2,:,4) 
            elseif (mr.eq.5) then
c..sp3d2-5
              tlmwf(0,:,nwf) =  (1./(sqrt(6.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) = -(1./(sqrt(2.)))*tlm(1,:,1)
              tlmwf(2,:,nwf) =  (1./(sqrt(3.)))*tlm(2,:,1)
            elseif (mr.eq.6) then
c..sp3d2-5
              tlmwf(0,:,nwf) =  (1./(sqrt(6.)))*tlm(0,:,1)
              tlmwf(1,:,nwf) =  (1./(sqrt(2.)))*tlm(1,:,1)
              tlmwf(2,:,nwf) =  (1./(sqrt(3.)))*tlm(2,:,1)
            endif

         else
            CALL juDFT_error("no tlmw for this lr",calledby ="wann_tlmw")
         endif

      enddo

      end subroutine wann_tlmw
      end module m_wann_tlmw
