!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mixing_history
  USE m_types_mixvector
  IMPLICIT NONE
  PRIVATE
  INTEGER:: iter_stored=0
  TYPE(t_mixvector),ALLOCATABLE::sm_store(:),fsm_store(:)
  PUBLIC :: mixing_history,mixing_history_reset,mixing_history_store
  PUBLIC :: mixing_history_open,mixing_history_close,mixing_history_limit
CONTAINS
  
  SUBROUTINE mixing_history_open(mpi,maxiter)
    USE m_types,ONLY:t_mpi
    INTEGER,INTENT(IN)    :: maxiter
    TYPE(t_mpi),INTENT(in):: mpi
    
    CHARACTER(len=20):: filename
    LOGICAL          :: l_fileexist
    INTEGER          :: n


    IF (iter_stored>0) RETURN ! History in memory found, no need to do IO
    IF (mpi%isize>1) THEN
       WRITE(filename,'(a,i0)') "mixing_history.",mpi%irank
    ELSE
       filename="mixing_history"
    ENDIF
    INQUIRE(file=filename,exist=l_fileexist)
    IF (.NOT.l_fileexist) RETURN !No previous data
    OPEN(888,file=filename,status='old',form='unformatted')
    READ(888) iter_stored
    IF (.NOT.ALLOCATED(sm_store)) ALLOCATE(sm_store(maxiter),fsm_store(maxiter))
    DO n=1,MIN(iter_stored,maxiter)
       READ(888) sm_store(n)
       READ(888) fsm_store(n)
    ENDDO
    CLOSE(888)
  END SUBROUTINE mixing_history_open

  SUBROUTINE mixing_history_close(mpi)
    USE m_types,ONLY:t_mpi
    TYPE(t_mpi),INTENT(in):: mpi
    
    CHARACTER(len=20):: filename
    INTEGER          :: n


    IF (iter_stored==0) RETURN ! Nothing found to be stored
    IF (mpi%isize>1) THEN
       WRITE(filename,'(a,i0)') "mixing_history.",mpi%irank
    ELSE
       filename="mixing_history"
    ENDIF
    OPEN(888,file=filename,form='unformatted',status='replace')
    WRITE(888) iter_stored
    DO n=1,iter_stored
       WRITE(888) sm_store(n)
       WRITE(888) fsm_store(n)
    ENDDO
    CLOSE(888)
    DEALLOCATE(sm_store,fsm_store)
    iter_stored=0
  END SUBROUTINE mixing_history_close
    
  
  SUBROUTINE mixing_history(imix,maxiter,inden,outden,sm,fsm,it)
    USE m_types
    implicit none
    INTEGER,INTENT(in)::imix,maxiter
    type(t_potden),intent(inout)::inden,outden
    type(t_mixvector),ALLOCATABLE::sm(:),fsm(:)
    INTEGER,INTENT(out)::it

    INTEGER:: n

    if (.not.allocated(sm_store)) THEN
       allocate(sm_store(maxiter),fsm_store(maxiter))
    endif
    IF (iter_stored+1==maxiter.AND.imix<8) iter_stored=0 !This is a broyden method which has to 
                                                            !be reset as soon as maxiter is reached
    it=iter_stored+1
    allocate(sm(it),fsm(it))
    CALL sm(it)%alloc()
    CALL fsm(it)%alloc()
    CALL sm(it)%from_density(inDen)
    CALL fsm(it)%from_density(outDen)
    !store the difference fsm - sm in fsm
    fsm(it) = fsm(it) - sm(it)
    do n=1,it-1 !Copy from storage
       sm(n)=sm_store(n)
       fsm(n)=fsm_store(n)
    ENDDO
    if(iter_stored<maxiter) THEN
       iter_stored=iter_stored+1
       sm_store(iter_stored)=sm(iter_stored)
       fsm_store(iter_stored)=fsm(iter_stored)
    else
       sm_store(:maxiter)=sm(2:maxiter+1)
       fsm_store(:maxiter)=fsm(2:maxiter+1)
    endif
  end subroutine mixing_history

  SUBROUTINE mixing_history_reset(mpi)
    USE m_types,ONLY:t_mpi
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(in)::mpi
    iter_stored=0
    PRINT *, "Reset of history"
    IF (mpi%irank==0) CALL system('rm mixing_history*')
  END SUBROUTINE mixing_history_reset

  subroutine mixing_history_limit(len)
    IMPLICIT NONE
    INTEGER,INTENT(in)::len

    if (iter_stored>len) then
       fsm_store(:len)=fsm_store(iter_stored-len+1:iter_stored)
       sm_store(:len)=sm_store(iter_stored-len+1:iter_stored)
       iter_stored=len
    end if
  end subroutine mixing_history_limit
  
  SUBROUTINE mixing_history_store(fsm)
    IMPLICIT NONE
    TYPE(t_mixvector),INTENT(IN)::fsm
    IF (iter_stored>0) fsm_store(iter_stored)=fsm
  END SUBROUTINE mixing_history_store
end MODULE m_mixing_history
