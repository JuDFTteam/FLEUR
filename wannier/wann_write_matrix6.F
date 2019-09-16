!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_write_matrix6
      contains
      subroutine wann_write_matrix6(
     >               mpi_comm,l_p0,filename,title,
     >               spin1,spin2,num_bands1,num_bands2,
     >               num_dims,fullnkpts,
     >               irank,isize,l_unformatted,
     <               matrix6)
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
      integer, intent(in)    :: num_bands1
      integer, intent(in)    :: num_bands2

      integer, intent(in)    :: num_dims
      integer, intent(in)    :: fullnkpts

      integer, intent(in)    :: irank,isize
      logical, intent(in)    :: l_unformatted

      complex, intent(inout) :: matrix6(:,:,:,:,:,:)

      integer :: ikpt,i,j,ii,jj
      integer :: cpu_index,dir
#ifdef CPP_MPI
      include 'mpif.h'
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#include "cpp_double.h"
#endif

#ifdef CPP_MPI
c**********************************************************
c     Collect contributions to the matrix6 matrix from the
c     various processors.
c**********************************************************
      if(isize.ne.1)then
       do ikpt=1,fullnkpts
        if(l_p0)then
         do cpu_index=1,isize-1
          if(mod(ikpt-1,isize).eq.cpu_index)then
           call MPI_RECV(
     &             matrix6(1:spin1,1:spin2,1:num_bands1,
     &             1:num_bands2,1:num_dims,ikpt),
     &             spin1*spin2*num_bands1*num_bands2*num_dims,
     &             CPP_MPI_COMPLEX,cpu_index,
     &             ikpt,mpi_comm,stt,ierr)
          endif !processors
         enddo !cpu_index
        else
         if(mod(ikpt-1,isize).eq.irank)then
           call MPI_SEND(
     &             matrix6(1:spin1,1:spin2,1:num_bands1,
     &             1:num_bands2,1:num_dims,ikpt),
     &             spin1*spin2*num_bands1*num_bands2*num_dims,
     &             CPP_MPI_COMPLEX,0,
     &             ikpt,mpi_comm,ierr)
         endif !processors
        endif ! l_p0
        call MPI_BARRIER(mpi_comm,ierr)
       enddo !ikpt 
      endif !isize
#endif

      write(6,*)"wann_write_matrix6"

      if(l_p0)then
       if(l_unformatted)then
	open(305,file=trim(filename)//'_unf',form='unformatted')
        write(305)spin2,spin1,num_bands1,num_bands2,
     &          num_dims,fullnkpts
        do ikpt=1,fullnkpts
	 write(305)matrix6(:,:,:,:,:,ikpt)
        enddo
       else
        open (305,file=filename)
       write (305,*)title
       write (305,'(6i5)') spin2,spin1,num_bands1,num_bands2,
     &          num_dims,fullnkpts
        do ikpt=1,fullnkpts
         do dir=1,num_dims  
          do i = 1,num_bands2
           do j = 1,num_bands1
            do ii=1,spin1
             do jj=1,spin2
              write (305,'(6i5,3x,2f18.12)') jj,ii,j,i,dir,ikpt,
     &              real(matrix6(jj,ii,j,i,dir,ikpt)),
     &              aimag(matrix6(jj,ii,j,i,dir,ikpt))
             enddo !jj 
            enddo !ii
           enddo !j
          enddo !i
         enddo !dir 
        enddo !ikpt
       endif 
       close(305)
      endif

      end subroutine wann_write_matrix6
      end module m_wann_write_matrix6
