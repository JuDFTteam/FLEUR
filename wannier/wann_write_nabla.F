!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_write_nabla
      contains
      subroutine wann_write_nabla(
     >               mpi_comm,l_p0,filename,title,
     >               nbnd,fullnkpts,nwfs,
     >               irank,isize,
     <               nablamat)
c*************************************************************
c     This subroutine is used to write several matrices to 
c     files: WF1.nabl, WF1.surfcurr, etc. The corresponding 
c     filename has to be provided as input. To be concrete
c     all explanations given in the following refer to
c     WF1.nabl/WF2.nabl.
c
c     MPI-Version: Collect the contributions to the matrix
c     grad^{k}_{mn} from the various processors.
c
c     Write the matrix grad^{k}_{mn} to file WF1.nabl/WF2.nabl
c
c     Frank Freimuth
c*************************************************************
      implicit none
      integer, intent(in)    :: mpi_comm
      logical, intent(in)    :: l_p0
      character, intent(in)  :: filename*(*)
      character, intent(in)  :: title*(*)

      integer, intent(in)    :: nbnd
      integer, intent(in)    :: fullnkpts
      integer, intent(in)    :: nwfs

      integer, intent(in)    :: irank,isize

      complex, intent(inout) :: nablamat(:,:,:,:)

      integer :: ikpt,i,j,k
      integer :: cpu_index
#ifdef CPP_MPI
      include 'mpif.h'
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#include "cpp_double.h"
#endif

#ifdef CPP_MPI
c**********************************************************
c     Collect contributions to the nablamat matrix from the
c     various processors.
c**********************************************************
      if(isize.ne.1)then
       do ikpt=1,fullnkpts
        if(l_p0)then
         do cpu_index=1,isize-1
          if(mod(ikpt-1,isize).eq.cpu_index)then
           call MPI_RECV(
     &             nablamat(1:3,1:nbnd,1:nbnd,ikpt),nbnd*nbnd*3,
     &             CPP_MPI_COMPLEX,cpu_index,
     &             ikpt,mpi_comm,stt,ierr)
          endif !processors
         enddo !cpu_index
        else
         if(mod(ikpt-1,isize).eq.irank)then
           call MPI_SEND(
     &             nablamat(1:3,1:nbnd,1:nbnd,ikpt),nbnd*nbnd*3,
     &             CPP_MPI_COMPLEX,0,
     &             ikpt,mpi_comm,ierr)
         endif !processors
        endif ! l_p0
        call MPI_BARRIER(mpi_comm,ierr)
       enddo !ikpt 
      endif !isize
#endif

      write(*,*)"wn: fullnkpts=",fullnkpts
      write(6,*)"wn: fullnkpts=",fullnkpts

      if(l_p0)then
       open (305,file=filename)
       write (305,*)title
       write (305,'(3i5)') nbnd,nbnd,fullnkpts
       do ikpt=1,fullnkpts
       do i = 1,nbnd
        do j = 1,nbnd
         do k = 1,3  
          write (305,'(3i5,3x,2f18.12)') i,j,ikpt,
     &              real(nablamat(k,j,i,ikpt)),
     &              aimag(nablamat(k,j,i,ikpt))
         enddo !k
        enddo !j
       enddo !i
       enddo !ikpt
       close(305)
      endif

      end subroutine wann_write_nabla
      end module m_wann_write_nabla
