c**************************c
c   write out uHu matrix   c
c**************************c
      module m_wann_write_uHu
      contains
      subroutine wann_write_uHu(
     >               jspin2,l_p0,fullnkpts,nntot,nntot2,wann,
     >               nbnd,bpt,gb,isize,irank,fending,ftype,
     <               uHu_in,nkpt_loc,counts,displs,nnodes,
     >               l_unformatted,l_symcc,l_check)
      use m_types
      use m_wann_uHu_symcheck

      implicit none
      integer, intent(in)     :: jspin2
      logical, intent(in)     :: l_p0,l_unformatted,l_symcc,l_check
      integer, intent(in)     :: fullnkpts,nkpt_loc
      integer, intent(in)     :: nntot,nntot2,nnodes
      type(t_wann),intent(in) :: wann

      integer, intent(in)     :: nbnd
      integer, intent(in)     :: bpt(nntot,fullnkpts)
      integer, intent(in)     :: gb(3,nntot,fullnkpts)
      integer, intent(in)     :: counts(0:nnodes-1),displs(0:nnodes-1)

      integer, intent(in)     :: isize,irank
      
      CHARACTER(len=12), INTENT(IN) :: fending  !for file ending
      CHARACTER(len=*),  INTENT(IN) :: ftype
      complex, intent(inout)  :: uHu_in(nbnd,nbnd,nntot2,nntot,nkpt_loc)

      complex, allocatable :: uHu(:,:,:,:,:)
      integer          :: ikpt,i,j,length
      integer          :: ikpt_b,ikpt_b2
      character(len=3) :: spin12(2)
      integer          :: cpu_index
      character(len=60)             :: header
      data   spin12/'WF1' , 'WF2'/

#ifdef CPP_MPI
      include 'mpif.h'
      integer :: ierr(3)
      integer :: stt(MPI_STATUS_SIZE)
#include "cpp_double.h"
#endif

      if(isize.gt.1) then
      if(l_p0) allocate(uHu(nbnd,nbnd,nntot2,nntot,fullnkpts))
#ifdef CPP_MPI
c******************************************************
c     Collect contributions to the mmnk matrix from the
c     various processors.
c******************************************************
      length = nbnd*nbnd*nntot2*nntot
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_GATHERV(
     >       uHu_in,length*nkpt_loc,MPI_DOUBLE_COMPLEX,
     >       uHu,length*counts,length*displs,MPI_DOUBLE_COMPLEX,
     >       0,MPI_COMM_WORLD,ierr)
#else
c      uHu = uHu_in
#endif
      endif

      header='Elements uHu at k+b1 and k+b2'

c******************************************************
c     Write mmnk matrix to file.
c******************************************************
      if (l_p0) then
       write(*,*)'symmetry-complete uHu: ',l_symcc
 
       if(l_symcc.and.(nntot.ne.nntot2)) stop 'wann_write_uHu'

       if(.false. .and. l_symcc) then
       ! exploit symmetry to complete matrix
       do ikpt = 1,fullnkpts
        do ikpt_b = 1,nntot
         do ikpt_b2 = 1,ikpt_b-1
          do i=1,nbnd
           do j=1,nbnd
            if(isize.gt.1) then
            uHu(j,i,ikpt_b,ikpt_b2,ikpt)
     >       = conjg(uHu(i,j,ikpt_b2,ikpt_b,ikpt))
            else
            uHu_in(j,i,ikpt_b,ikpt_b2,ikpt)
     >       = conjg(uHu_in(i,j,ikpt_b2,ikpt_b,ikpt))
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
       endif

       if(.not.l_unformatted) then
        open (305,file=spin12(jspin2)//trim(fending)//'.uHu'
     >                               //trim(ftype))
        write (305,*) 'Elements uHu at k+b1 and k+b2'
        write (305,'(4i5)') nbnd,fullnkpts,nntot,nntot2
        write (305,*)"   "
        write (305,*)"   "
        do ikpt = 1,fullnkpts
         do ikpt_b = 1,nntot
          do ikpt_b2 = 1,nntot2
    !         write(305,'(i6,i6,i6)')ikpt,ikpt_b,ikpt_b2
           do i = 1,nbnd
            do j = 1,nbnd
             if(isize.gt.1) then
             write (305,'(2f24.18)')
     &          real(uHu(j,i,ikpt_b2,ikpt_b,ikpt)),
!     &         -aimag(uHu(j,i,ikpt_b2,ikpt_b,ikpt))
     &         aimag(uHu(j,i,ikpt_b2,ikpt_b,ikpt))
             else
             write (305,'(2f24.18)')
     &          real(uHu_in(j,i,ikpt_b2,ikpt_b,ikpt)),
!     &         -aimag(uHu_in(j,i,ikpt_b2,ikpt_b,ikpt))
     &         aimag(uHu_in(j,i,ikpt_b2,ikpt_b,ikpt))
             endif
            enddo
           enddo
          enddo
         enddo
        enddo !ikpt
        close (305)
       else
        open (305,file=spin12(jspin2)//trim(fending)//'.uHu'
     >                               //trim(ftype),form='unformatted')
        write (305) header
        write (305) nbnd,fullnkpts,nntot !,nntot2
!        write (305) bpt,gb
        do ikpt = 1,fullnkpts
         do ikpt_b = 1,nntot
          do ikpt_b2 = 1,nntot2
           if(isize.gt.1) then
!        write (305) conjg(uHu)
             write (305) uHu(:,:,ikpt_b2,ikpt_b,ikpt)
           else
!        write (305) conjg(uHu_in)
             write (305) uHu_in(:,:,ikpt_b2,ikpt_b,ikpt)
           endif
           
          enddo
         enddo
        enddo
        
        close(305)
       endif
       

       if((trim(ftype).ne.'_kq').and.l_check) then
        write(*,*)'perform symcheck...'
        if(isize.gt.1) then
           CALL wann_uHu_symcheck(uHu,nbnd,nntot,nntot2,fullnkpts)
        else
           CALL wann_uHu_symcheck(uHu_in,nbnd,nntot,nntot2,fullnkpts)
        endif
       endif

      endif !l_p0

      if(allocated(uHu)) deallocate( uHu )

#ifdef CPP_MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      end subroutine wann_write_uHu
      end module m_wann_write_uHu
