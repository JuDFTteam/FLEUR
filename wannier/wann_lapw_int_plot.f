      MODULE m_wann_lapw_int_plot

c***************************************************************
c         plot wannierfunction in interstitial
c         Frank Freimuth, November 2006
c***************************************************************

      CONTAINS
      SUBROUTINE wann_lapw_int_plot(point,bmat,unigrid,wannint,
     <                              xdnout)
      use m_cotra, only: cotra1
      use m_constants, only: pimach
      implicit none
      real,intent(in)::point(3)
      real,intent(in)::bmat(3,3)
      integer,intent(in)::unigrid(4)
      complex,intent(in)::wannint(-unigrid(4):unigrid(4),
     ,                            -unigrid(4):unigrid(4),
     ,                            -unigrid(4):unigrid(4))
      complex,intent(out)::xdnout

      real xiiu(3),tpi,arg
      complex ci,factor1,factor2,factor3
      complex factor1p,factor2p,factor3p,factor
      integer pw1,pw2,pw3

      tpi=2*pimach()
      call cotra1(point,xiiu,bmat)



      xdnout=cmplx(0.0,0.0)

      if(.false.)then

      arg=(tpi*xiiu(1))/unigrid(1)
      factor1=cmplx(cos(arg),sin(arg))
      arg=(tpi*xiiu(2))/unigrid(2)
      factor2=cmplx(cos(arg),sin(arg))
      arg=(tpi*xiiu(3))/unigrid(3)
      factor3=cmplx(cos(arg),sin(arg))

      do pw3=-unigrid(4),unigrid(4)
       factor3p=(factor3)**pw3
       do pw2=-unigrid(4),unigrid(4)
        factor2p=(factor2)**pw2
        do pw1=-unigrid(4),unigrid(4)
         factor1p=(factor1)**pw1

          factor=factor1p*factor2p*factor3p
          xdnout=xdnout+factor*wannint(pw1,pw2,pw3)

        enddo
       enddo
      enddo 

      else

      do pw3=-unigrid(4),unigrid(4)
       do pw2=-unigrid(4),unigrid(4)
        do pw1=-unigrid(4),unigrid(4)
          arg=tpi*(pw1*xiiu(1)/unigrid(1)+pw2*xiiu(2)/unigrid(2)+
     +                  pw3*xiiu(3)/unigrid(3)) 
          factor=cmplx(cos(arg),sin(arg))
          xdnout=xdnout+factor*wannint(pw1,pw2,pw3)

        enddo
       enddo
      enddo 


      endif   
      END SUBROUTINE wann_lapw_int_plot
      END MODULE m_wann_lapw_int_plot
