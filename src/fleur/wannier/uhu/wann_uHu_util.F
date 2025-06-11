      module m_wann_uHu_util
      implicit none
      contains

      subroutine array_split(nkpts,nnodes,counts,displs)
      implicit none
      integer, intent(in) :: nkpts,nnodes
      integer, intent(inout) :: counts(0:nnodes-1),displs(0:nnodes-1)
      integer :: ratio,remainder,i

      counts = 0; displs = 0

      ratio = nkpts / nnodes
      remainder = mod(nkpts,nnodes)
      do i=0,nnodes-1
       if(i<remainder) then
        counts(i) = ratio+1
        displs(i) = i*(ratio+1)
       else
        counts(i) = ratio
        displs(i) = remainder*(ratio+1) + (i-remainder)*ratio
       endif
      enddo

      end subroutine array_split


      integer function kpt_on_node(nkpt,nnodes,counts,displs)
      implicit none
      integer,intent(in) :: nkpt,nnodes
      integer,intent(in) :: counts(0:nnodes-1),displs(0:nnodes-1)
      integer :: i

      if(nnodes.eq.1) then
       kpt_on_node=0
       return
      endif

      if(nkpt.le.counts(0)) then
       kpt_on_node=0
       return
      elseif(nkpt.gt.displs(nnodes-1)) then
       kpt_on_node = nnodes-1
       return
      endif

      do i=2,nnodes-1
       kpt_on_node=i-1
       if(nkpt.le.displs(i)) return
      enddo

      end function kpt_on_node

      end module m_wann_uHu_util
