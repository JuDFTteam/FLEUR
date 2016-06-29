!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_readcenters
      use m_juDFT
      contains
      subroutine wann_readcenters(
     >               nwfs,jspin,
     <               centers,spreads)
c******************************************************
c     Read the centers and spreads (optional) from the
c     file "WF"jspin".1"
c     Frank Freimuth
c******************************************************
      implicit none
      integer,intent(in)        :: nwfs,jspin
      real,intent(out)          :: centers(3,nwfs)
      real,intent(out),optional :: spreads(nwfs)

      logical                   :: l_readspread
      character(len=3)          :: spin12(2)
      data  spin12/'WF1','WF2'/
      integer                   :: line,ios,linefinal
      integer                   :: wanind,wanindtmp
      character(len=30)         :: task

      line=0
      linefinal=0
      l_readspread=.false.
      if(present(spreads)) l_readspread=.true.
      open(100,file=spin12(jspin)//'.wout')
      do while(.true.) 
        read(100,'(a)',iostat=ios)task
        if(ios.ne.0)exit
        line=line+1
        if(index(task,"Final State").ne.0)linefinal=line
      enddo
      IF(linefinal==0) CALL juDFT_error("Final State not found",calledby
     +     ="wann_readcenters")
      rewind(100)
      line=0
      do while(.true.)
        read(100,'(a)',iostat=ios)task
        if(ios.ne.0)exit
        line=line+1
        if(line.eq.linefinal)exit
      enddo
      wanind=0
      do while(.true.)
         read(100,'(a)',iostat=ios)task
         if(ios.ne.0)exit
         if(index(task,"Sum of centres and spreads").ne.0)exit
         wanind=wanind+1
         backspace(100)
         IF(wanind>nwfs) CALL juDFT_error("wanind.gt.nwfs",calledby
     +        ="wann_readcenters")
         if(l_readspread)then
           read(100,fmt=1000)wanindtmp,centers(:,wanind),spreads(wanind)
         else
           read(100,fmt=2000)wanindtmp,centers(:,wanind)
         endif
         IF(wanindtmp/=wanind) CALL juDFT_error("wanindtmp",calledby
     +        ="wann_readcenters")
      enddo
      close(100)
c      do wanind=1,wanindtmp
c         print*,centres(:,wanind),spreads(wanind)
c      enddo   
1000  format(22x,i5,3x,f10.6,1x,f10.6,1x,f10.6,2x,f15.8)
2000  format(22x,i5,3x,f10.6,1x,f10.6,1x,f10.6)      
      end subroutine wann_readcenters
      end module m_wann_readcenters
