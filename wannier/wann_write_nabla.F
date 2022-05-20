!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_write_nabla
            use m_juDFT
      contains
      subroutine wann_write_nabla(
     >               fmpi_comm,l_p0,filename,title,
     >               nbnd,fullnkpts,nwfs,
     >               irank,isize,l_unformatted,
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

      USE m_constants
#ifdef CPP_MPI
      USE mpi
#endif

      implicit none

      integer, intent(in)    :: fmpi_comm
      logical, intent(in)    :: l_p0
      character, intent(in)  :: filename*(*)
      character, intent(in)  :: title*(*)

      integer, intent(in)    :: nbnd
      integer, intent(in)    :: fullnkpts
      integer, intent(in)    :: nwfs

      integer, intent(in)    :: irank,isize
      logical, intent(in)    :: l_unformatted

      complex, intent(inout) :: nablamat(:,:,:,:)

      integer :: ikpt,i,j,k
      integer :: cpu_index
#ifdef CPP_MPI
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#endif

#ifdef CPP_MPI
c**********************************************************
c     Collect contributions to the nablamat matrix from the
c     various processors.
c**********************************************************
      call timestart("wann_write_nabla")

      if(isize.ne.1)then
       do ikpt=1,fullnkpts
        if(l_p0)then
         do cpu_index=1,isize-1
          if(mod(ikpt-1,isize).eq.cpu_index)then
           call MPI_RECV(
     &             nablamat(1:3,1:nbnd,1:nbnd,ikpt),nbnd*nbnd*3,
     &             MPI_DOUBLE_COMPLEX,cpu_index,
     &             ikpt,fmpi_comm,stt,ierr(1))
          endif !processors
         enddo !cpu_index
        else
         if(mod(ikpt-1,isize).eq.irank)then
           call MPI_SEND(
     &             nablamat(1:3,1:nbnd,1:nbnd,ikpt),nbnd*nbnd*3,
     &             MPI_DOUBLE_COMPLEX,0,
     &             ikpt,fmpi_comm,ierr(1))
         endif !processors
        endif ! l_p0
        call MPI_BARRIER(fmpi_comm,ierr(1))
       enddo !ikpt 
      endif !isize
#endif

      write(*,*)"wn: fullnkpts=",fullnkpts
      write(oUnit,*)"wn: fullnkpts=",fullnkpts

      if(l_p0)then
       if(l_unformatted)then
        open(305,file=trim(filename)//'_unf',form='unformatted')
        write(305)nbnd,nbnd,fullnkpts
        write(305)nablamat(1:3,1:nbnd,1:nbnd,1:fullnkpts)
       else !l_unformatted   
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
       endif !l_unformatted 
       close(305)
      endif

      call timestop("wann_write_nabla")
      end subroutine wann_write_nabla
      end module m_wann_write_nabla
