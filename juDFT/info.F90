!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_judft_info
  private
  integer:: info_index
  character(len=50),allocatable:: messages(:)
  character(len=20),allocatable:: groups(:)
  public judft_info,judft_write_infos
contains
  subroutine judft_info(message,group)
    implicit none
    character(len=*),intent(in)::message,group

    integer:: irank=0,ierr
#ifdef CPP_MPI
    include 'mpif.h'
    LOGICAL:: l_mpi
    CALL mpi_initialized(l_mpi,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
#else
    logical :: l_mpi=.false.
#endif

    if (l_mpi) THEN
       write(*,9000) irank,group,message 
    else
       write(*,8000) group,message
    endif
    
    !Store info in arrays
    info_index=info_index+1
    !perhaps (more) storage is needed?
    if (.not.allocated(messages)) call priv_reallocate()
    if (info_index>size(messages)) call priv_reallocate()
    
    messages(info_index)=message
    groups(info_index)=group
9000 FORMAT("PE:",i3,a20," WARNING: ",a50)
8000 FORMAT(a20,' WARNING: ',a50)
  end subroutine judft_info

  subroutine priv_reallocate()
    integer:: n
    character(len=50),allocatable:: old_m(:)
    character(len=20),allocatable:: old_g(:)
    if (.not.allocated(messages)) THEN
       info_index=0
       allocate(messages(10),groups(10))
       return
    endif
    
    n=size(messages)
    allocate(old_m(n),old_g(n))
    old_m=messages
    old_g=groups
    deallocate(messages,groups)
    allocate(messages(n+10),groups(n+10))
    messages(:n)=old_m
    groups(:n)=old_g
  end subroutine priv_reallocate

  subroutine judft_write_infos()
    integer::n
    
    write(6,*) "The following WARNING messages have been issued:"
    DO n=1,info_index
       write(6,8000) groups(n),messages(n)
    enddo
8000 FORMAT(a20,' WARNING: ',a50)
  end subroutine judft_write_infos
end module m_judft_info
