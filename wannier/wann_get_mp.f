      module m_wann_get_mp
      use m_juDFT
      contains
      subroutine wann_get_mp(
     >               nkpts,kpoints,
     <               num)
c**************************************
c     Determine the structure of the
c     Monkhorst-Pack mesh.
c     
c     Frank Freimuth
c**************************************
      implicit none
      integer,intent(in)  :: nkpts
      real,intent(in)     :: kpoints(:,:)
      integer,intent(out) :: num(:)
      
      integer :: dim,iter
      real    :: maxi,mini,increm,compare

      IF(SIZE(kpoints,1)/=3)      CALL juDFT_error("wann_get_mp: 1"
     +     ,calledby ="wann_get_mp")
      IF(SIZE(kpoints,2)/=nkpts)  CALL juDFT_error("wann_get_mp: 2"
     +     ,calledby ="wann_get_mp")
      IF(SIZE(num,1)/=3)          CALL juDFT_error("wann_get_mp: 3"
     +     ,calledby ="wann_get_mp")

      do dim=1,3
         maxi=maxval(kpoints(dim,:))
         mini=minval(kpoints(dim,:))
         if(mini==maxi)then
            num(dim)=1
         else   
            increm=maxi-mini
            do iter=1,nkpts
               compare=maxi-kpoints(dim,iter)
               if(abs(compare).lt.1e-6)cycle
               if(compare.lt.increm) then
                  increm=compare
               endif   
            enddo
            num(dim)=(maxi-mini)/increm+1.01
         endif   
      enddo
      write(6,*)"wann_get_mp: determination of mp-grid parameters:"
      write(6,*)"mp_1=",num(1),"mp_2=",num(2),"mp_3=",num(3)
      IF(num(1)*num(2)*num(3)/=nkpts)  CALL juDFT_error
     +     ("mysterious kpoints",calledby ="wann_get_mp")

      end subroutine wann_get_mp
      end module m_wann_get_mp
