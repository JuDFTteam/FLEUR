      module m_wann_write_mmnk2
      contains
      subroutine wann_write_mmnk2(
     >               l_p0,fullnkpts,nntot_q,wann,
     >               nbnd,bpt_q,gb_q,isize,irank,
     >               fname,mmnk_q,l_unformatted)
c**********************************************************

c**********************************************************
      use m_types
      use m_fleurenv
      implicit none
      logical, intent(in)     :: l_p0,l_unformatted
      integer, intent(in)     :: fullnkpts
      integer, intent(in)     :: nntot_q
      type(t_wann),intent(in) :: wann

      integer, intent(in)     :: nbnd
      integer, intent(in)     :: bpt_q(:)
      integer, intent(in)     :: gb_q(:,:)

      integer, intent(in)     :: isize,irank
      
      CHARACTER(len=30), INTENT(IN) :: fname
      complex, intent(in)  :: mmnk_q(:,:,:,:)

      integer          :: ikpt,i,j
      integer          :: ikpt_b
      character(len=3) :: spin12(2)
      integer          :: cpu_index
      data   spin12/'WF1' , 'WF2'/

#ifdef CPP_MPI
      include 'mpif.h'
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#include "cpp_double.h"
#endif

#ifdef CPP_MPI
c******************************************************
c     Collect contributions to the mmnk matrix from the
c     various processors.
c******************************************************
      if(isize.ne.1)then
      do ikpt=1,fullnkpts
       if(l_p0)then
        do cpu_index=1,isize-1
         if(mod(ikpt-1,isize).eq.cpu_index)then
          do ikpt_b=1,nntot_q !nearest neighbors
             call MPI_RECV(
     &               mmnk_q(1:nbnd,1:nbnd,ikpt_b,ikpt),nbnd*nbnd,
     &               CPP_MPI_COMPLEX,cpu_index,5*fullnkpts,
     &               MPI_COMM_WORLD,stt,ierr)

          enddo !nearest neighbors
         endif !processors
        enddo !cpu_index
       else
        if(mod(ikpt-1,isize).eq.irank)then
         do ikpt_b=1,nntot_q !loop over nearest neighbors
            call MPI_SEND(
     &              mmnk_q(1:nbnd,1:nbnd,ikpt_b,ikpt),
     &              nbnd*nbnd,CPP_MPI_COMPLEX,0,5*fullnkpts,
     &              MPI_COMM_WORLD,ierr)
         enddo !loop over nearest neighbors
        endif !processors
       endif ! l_p0
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      enddo !ikpt
      endif !isize
#endif


c******************************************************
c     Write mmnk matrix to file.
c******************************************************
      if (l_p0) then
        if(.not.l_unformatted) then
         open(305,file=trim(fname))
         write(305,*)'Overlaps between parameter points'
         write(305,'(3i5)')nbnd,fullnkpts,nntot_q
         do ikpt=1,fullnkpts
            do ikpt_b=1,nntot_q
               write(305,'(2i5,3x,3i4)')ikpt,bpt_q(ikpt_b),
     >                                  gb_q(1:3,ikpt_b)
            do i=1,nbnd
               do j=1,nbnd
                  write(305,'(2f24.18)')
     >             real(mmnk_q(j,i,ikpt_b,ikpt)),
     >            -aimag(mmnk_q(j,i,ikpt_b,ikpt))
               enddo
            enddo
            enddo
         enddo
         close(305)
        else
         open(305,file=trim(fname),form='unformatted')
         write(305)nbnd,fullnkpts,nntot_q
         write(305)bpt_q,gb_q
         write(305)conjg(mmnk_q)     
         close(305)
        endif
      endif !l_p0

      end subroutine wann_write_mmnk2
      end module m_wann_write_mmnk2
