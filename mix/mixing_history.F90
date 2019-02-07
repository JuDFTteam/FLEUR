!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mixing_history
  use m_types_mixvector
  implicit none
  integer:: iter_stored=0
  type(t_mixvector),allocatable::sm_store(:),fsm_store(:)
contains
  
  SUBROUTINE mixing_history(imix,maxiter,inden,outden,sm,fsm,it)
    use m_types
    implicit none
    INTEGER,INTENT(in)::imix,maxiter
    type(t_potden),intent(in)::inden,outden
    type(t_mixvector),ALLOCATABLE::sm(:),fsm(:)
    INTEGER,INTENT(out)::it

    INTEGER:: n

    if (.not.allocated(sm_store)) THEN
       allocate(sm_store(maxiter),fsm_store(maxiter))
    endif
    IF (iter_stored+1==maxiter.AND.imix.NE.9) iter_stored=0 !This is a broyden method which has to 
                                                            !be reset as soon as maxiter is reached
    it=MIN(iter_stored+1,maxiter)
    allocate(sm(it),fsm(it))
    CALL sm(it)%alloc()
    CALL fsm(it)%alloc()
    CALL sm(it)%from_density(inDen)
    CALL fsm(it)%from_density(outDen)
    !store the difference fsm - sm in fsm
    fsm(it) = fsm(it) - sm(it)
    do n=it-1,1,-1 !Copy from storage
       sm(n)=sm_store(n)
       fsm(n)=fsm_store(n)
    ENDDO
    if(iter_stored<maxiter) THEN
       iter_stored=iter_stored+1
       sm_store(:iter_stored)=sm(:iter_stored)
       fsm_store(:iter_stored)=fsm(:iter_stored)
    else
       sm_store(:maxiter-1)=sm(2:maxiter)
       fsm_store(:maxiter-1)=fsm(2:maxiter)
    endif
  end subroutine mixing_history

  SUBROUTINE mixing_history_reset()
    IMPLICIT NONE
    iter_stored=0
  END SUBROUTINE mixing_history_reset

end MODULE m_mixing_history
