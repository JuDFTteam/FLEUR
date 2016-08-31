!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_write_mmnk
      contains
      subroutine wann_write_mmnk(
     >               mpi_comm,jspin2,l_p0,fullnkpts,nntot,wann,
     >               maptopair,pair_to_do,nbnd,bpt,gb,
     >               isize,irank,
     <               mmnk)
c**********************************************************
c     MPI-Version: Collect the contributions to the matrix
c     M^{k,b}_{mn} from the various processors.
c
c     Symmetry: Compose the M^{k,b}_{mn} matrix for the
c     full BZ from the pieces of the IBZ.
c
c     Write the matrix M^{k,b}_{mn} to file WF1.mmn/WF2.mmn
c     Frank Freimuth
c**********************************************************
      use m_types
      use m_juDFT
      implicit none
      integer, intent(in)     :: jspin2,mpi_comm
      logical, intent(in)     :: l_p0
      integer, intent(in)     :: fullnkpts
      integer, intent(in)     :: nntot
      type(t_wann),intent(in) :: wann

      integer, intent(in)     :: maptopair(:,:,:) !maptopair(3,fullnkpts,nntot)
      integer, intent(in)     :: pair_to_do(:,:)  !pair_to_do(fullnkpts,nntot)
      integer, intent(in)     :: nbnd
      integer, intent(in)     :: bpt(:,:)
      integer, intent(in)     :: gb(:,:,:)

      integer, intent(in)     :: isize,irank

      complex, intent(inout)  :: mmnk(:,:,:,:) !mmnk(nbnd,nbnd,nntot,fullnkpts)

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
          do ikpt_b=1,nntot !nearest neighbors
           if(pair_to_do(ikpt,ikpt_b).ne.0)then
             call MPI_RECV(
     &               mmnk(1:nbnd,1:nbnd,ikpt_b,ikpt),nbnd*nbnd,
     &               CPP_MPI_COMPLEX,cpu_index,5*fullnkpts+
     &               pair_to_do(ikpt,ikpt_b),mpi_comm,stt,ierr)

           endif !pairtodo
          enddo !nearest neighbors
         endif !processors
        enddo !cpu_index
       else
        if(mod(ikpt-1,isize).eq.irank)then
         do ikpt_b=1,nntot !loop over nearest neighbors
          if(pair_to_do(ikpt,ikpt_b).ne.0)then
            call MPI_SEND(
     &              mmnk(1:nbnd,1:nbnd,ikpt_b,ikpt),
     &              nbnd*nbnd,CPP_MPI_COMPLEX,0,5*fullnkpts+
     &              pair_to_do(ikpt,ikpt_b),mpi_comm,ierr)
          endif !pairtodo
         enddo !loop over nearest neighbors
        endif !processors
       endif ! l_p0
       call MPI_BARRIER(mpi_comm,ierr)
      enddo !ikpt
      endif !isize
#endif

c****************************************************
c     Symmetry: complete the mmnk matrix.
c****************************************************
      if(l_p0)then 
       do ikpt=1,fullnkpts
        do ikpt_b=1,nntot
         if(pair_to_do(ikpt,ikpt_b).eq.0)then
          if(maptopair(3,ikpt,ikpt_b).eq.1)then !conjugation selection
             mmnk(:,:,ikpt_b,ikpt)=conjg(transpose(mmnk(:,:,
     &          maptopair(2,ikpt,ikpt_b),maptopair(1,ikpt,ikpt_b))))
          elseif(maptopair(3,ikpt,ikpt_b).eq.2)then !rotation
             mmnk(:,:,ikpt_b,ikpt)=mmnk(:,:,maptopair
     &         (2,ikpt,ikpt_b),maptopair(1,ikpt,ikpt_b))
          elseif(maptopair(3,ikpt,ikpt_b).eq.3)then !rotation&reflection
             mmnk(:,:,ikpt_b,ikpt)=conjg( mmnk(:,:,maptopair
     &         (2,ikpt,ikpt_b),maptopair(1,ikpt,ikpt_b)) )
          else !something wrong
             call juDFT_error('maptopair')
          endif!maptopair
         endif!pairtodo
        enddo!ikpt_b
       enddo!ikpt
      endif

c******************************************************
c     Write mmnk matrix to file.
c******************************************************
      if (l_p0) then
       open (305,file=spin12(jspin2)//'.mmn')
       write (305,*) 'Overlaps of the wavefunct. the k- and b-points'
       write (305,'(3i5)') nbnd,fullnkpts,nntot
       do ikpt = 1,fullnkpts
        do ikpt_b = 1,nntot
         write (305,'(2i5,3x,3i4)') ikpt,bpt(ikpt_b,ikpt),
     &                                  gb(1:3,ikpt_b,ikpt)
         do i = 1,nbnd
          do j = 1,nbnd
c           write (305,'(2f18.12)')
           write (305,'(2f24.18)')
     &         real(mmnk(j,i,ikpt_b,ikpt)),-aimag(mmnk(j,i,ikpt_b,ikpt))
          enddo
         enddo
        enddo
       enddo !ikpt
       close (305)
      endif !l_p0

      end subroutine wann_write_mmnk
      end module m_wann_write_mmnk
