!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_write_amn
      contains
      subroutine wann_write_amn(
     >               mpi_comm,l_p0,filename,title,
     >               nbnd,fullnkpts,nwfs,
     >               irank,isize,l_freeformat,
     <               amn,l_unformatted)
c**********************************************************
c     This subroutine is used to write several matrices to 
c     files: WF1.mmn, WF1.amn, etc. The corresponding 
c     filename has to be provided as input. To be concrete
c     all explanations given in the following refer to
c     WF1.amn/WF2.amn.
c
c     MPI-Version: Collect the contributions to the matrix
c     A^{k}_{mn} from the various processors.
c
c     Write the matrix A^{k}_{mn} to file WF1.amn/WF2.amn
c
c     Frank Freimuth
c**********************************************************
      implicit none
      integer, intent(in)    :: mpi_comm
      logical, intent(in)    :: l_p0,l_unformatted
      character, intent(in)  :: filename*(*)
      character, intent(in)  :: title*(*)

      integer, intent(in)    :: nbnd
      integer, intent(in)    :: fullnkpts
      integer, intent(in)    :: nwfs

      integer, intent(in)    :: irank,isize
      logical, intent(in)    :: l_freeformat

      complex, intent(inout) :: amn(:,:,:)

      integer :: ikpt,nwf,i
      integer :: cpu_index
#ifdef CPP_MPI
      include 'mpif.h'
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#include "cpp_double.h"
#endif

#ifdef CPP_MPI
c******************************************************
c     Collect contributions to the amn matrix from the
c     various processors.
c******************************************************
      if(isize.ne.1)then
       do ikpt=1,fullnkpts
        if(l_p0)then
         do cpu_index=1,isize-1
          if(mod(ikpt-1,isize).eq.cpu_index)then
           call MPI_RECV(
     &             amn(1:nbnd,1:nwfs,ikpt),nbnd*nwfs,
     &             CPP_MPI_COMPLEX,cpu_index,
     &             ikpt,mpi_comm,stt,ierr)
          endif !processors
         enddo !cpu_index
        else
         if(mod(ikpt-1,isize).eq.irank)then
           call MPI_SEND(
     &             amn(1:nbnd,1:nwfs,ikpt),nbnd*nwfs,
     &             CPP_MPI_COMPLEX,0,
     &             ikpt,mpi_comm,ierr)
         endif !processors
        endif ! l_p0
        call MPI_BARRIER(mpi_comm,ierr)
       enddo !ikpt 
      endif !isize
#endif

      if(l_p0)then
       if(l_unformatted)then
          open (305,file=filename,form='unformatted')
       elseif(l_freeformat)then
          open(305,file=filename,recl=1000)
       else
          open (305,file=filename)
       endif
       
       if(l_unformatted)then
        write(305)nbnd,fullnkpts,nwfs
        write(305)amn
       else
        write (305,*)title
        write (305,'(i5,i7,i5)') nbnd,fullnkpts,nwfs
        if(l_freeformat)then
         do ikpt = 1,fullnkpts
          do nwf = 1,nwfs
           do i = 1,nbnd
            write (305,*) i,nwf,ikpt,
     &              real(amn(i,nwf,ikpt)),aimag(amn(i,nwf,ikpt))
           enddo !i
          enddo !nwf
         enddo !ikpt
        else
         do ikpt = 1,fullnkpts
          do nwf = 1,nwfs
           do i = 1,nbnd
            write (305,'(i5,i5,i7,3x,2f18.12)') i,nwf,ikpt,
     &              real(amn(i,nwf,ikpt)),aimag(amn(i,nwf,ikpt))
           enddo !i
          enddo !nwf
         enddo !ikpt
        endif 
       endif
       close(305)
      endif

      end subroutine wann_write_amn
      end module m_wann_write_amn
