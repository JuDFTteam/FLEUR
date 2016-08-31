!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_write_hsomtx
      contains
      subroutine wann_write_hsomtx(
     >               mpi_comm,l_p0,filename,title,
     >               spin1,spin2,nbnd,nwfs,fullnkpts,
     >               irank,isize,
     <               hsomtx)
c*************************************************************
c     This subroutine is used to write several matrices to 
c     files. The corresponding 
c     filename has to be provided as input. 
c
c     MPI-Version: Collect the contributions to the matrix
c     from the various processors.
c
c     Frank Freimuth
c*************************************************************
      implicit none
      integer, intent(in)    :: mpi_comm
      logical, intent(in)    :: l_p0
      character, intent(in)  :: filename*(*)
      character, intent(in)  :: title*(*)

      integer, intent(in)    :: spin1
      integer, intent(in)    :: spin2
      integer, intent(in)    :: nbnd
      integer, intent(in)    :: nwfs
      integer, intent(in)    :: fullnkpts

      integer, intent(in)    :: irank,isize

      complex, intent(inout) :: hsomtx(:,:,:,:,:)

      integer :: ikpt,i,j,ii,jj
      integer :: cpu_index
#ifdef CPP_MPI
      include 'mpif.h'
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#include "cpp_double.h"
#endif

#ifdef CPP_MPI
c**********************************************************
c     Collect contributions to the hsomtx matrix from the
c     various processors.
c**********************************************************
      if(isize.ne.1)then
       do ikpt=1,fullnkpts
        if(l_p0)then
         do cpu_index=1,isize-1
          if(mod(ikpt-1,isize).eq.cpu_index)then
           call MPI_RECV(
     &             hsomtx(1:spin1,1:spin2,1:nbnd,
     &             1:nwfs,ikpt),spin1*spin2*nbnd*nwfs,
     &             CPP_MPI_COMPLEX,cpu_index,
     &             ikpt,mpi_comm,stt,ierr)
          endif !processors
         enddo !cpu_index
        else
         if(mod(ikpt-1,isize).eq.irank)then
           call MPI_SEND(
     &             hsomtx(1:spin1,1:spin2,1:nbnd,
     &             1:nwfs,ikpt),spin1*spin2*nbnd*nwfs,
     &             CPP_MPI_COMPLEX,0,
     &             ikpt,mpi_comm,ierr)
         endif !processors
        endif ! l_p0
        call MPI_BARRIER(mpi_comm,ierr)
       enddo !ikpt 
      endif !isize
#endif

      write(6,*)"wann_write_hsomtx"

      if(l_p0)then
       open (305,file=filename)
c       write (305,*)title
c       write (305,'(3i5)') nbnd,nbnd,fullnkpts
       do ikpt=1,fullnkpts
        do i = 1,nwfs
         do j = 1,nbnd
          do ii=1,spin1
           do jj=1,spin2
            write (305,'(5i5,3x,2f18.12)') jj,ii,j,i,ikpt,
     &              real(hsomtx(jj,ii,j,i,ikpt)),
     &              aimag(hsomtx(jj,ii,j,i,ikpt))
           enddo !jj 
          enddo !ii
         enddo !j
        enddo !i
       enddo !ikpt
       close(305)
      endif

      end subroutine wann_write_hsomtx
      end module m_wann_write_hsomtx
