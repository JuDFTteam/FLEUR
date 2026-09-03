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
  PUBLIC :: dfpt_mixing_history_reset
CONTAINS

  SUBROUTINE mixing_history_open(mpi,maxiter,basename)
    USE m_types,ONLY:t_mpi
    INTEGER,INTENT(IN)    :: maxiter
    TYPE(t_mpi),INTENT(in):: mpi

    CHARACTER(len=20), OPTIONAL, INTENT(IN) :: basename

   CHARACTER(len=80):: filename
    LOGICAL          :: l_fileexist
   LOGICAL          :: l_has_history
   INTEGER          :: n, expected_history_files, existing_history_files


    IF (iter_stored>0) RETURN ! History in memory found, no need to do IO
    CALL mixing_history_file_count(mpi,expected_history_files,existing_history_files,l_has_history,basename)
    IF (l_has_history .AND. existing_history_files /= expected_history_files) THEN
       IF (mpi%irank==0) THEN
          WRITE(*,'(a,i0,a,i0,a)') 'Reset of history: found ', existing_history_files, &
                                   ' mixing history file(s) for ', mpi%isize, ' MPI task(s).'
       END IF
       CALL mixing_history_reset(mpi,basename)
       RETURN
    END IF

    IF (mpi%isize>1) THEN
       IF (.NOT.PRESENT(basename)) THEN
         WRITE(filename,'(a,i0)') "mixing_history.",mpi%irank
       ELSE
          WRITE(filename,'(a,i0)') TRIM(basename)//"_mixing_history.",mpi%irank
       END IF
    ELSE
       IF (.NOT.PRESENT(basename)) THEN
          filename="mixing_history"
       ELSE
          filename=TRIM(basename)//"mixing_history"
       END IF
    ENDIF
    INQUIRE(file=filename,exist=l_fileexist)
    IF (.NOT.l_fileexist) RETURN !No previous data
! I comment out this extra code path for the PGI compiler. It seems to be
! not needed.
!#ifdef CPP_GPU
!    PRINT *,"Warning PGI compiler does not support reading of history"
!#else
    OPEN(888,file=filename,status='old',form='unformatted')
    READ(888) iter_stored
    IF (.NOT.ALLOCATED(sm_store)) ALLOCATE(sm_store(maxiter),fsm_store(maxiter))
    DO n=1,MIN(iter_stored,maxiter)
       call sm_store(n)%read_unformatted(888)
       call fsm_store(n)%read_unformatted(888)
    ENDDO
    CLOSE(888,status='DELETE') !The same history should not be read twice!
!#endif
  END SUBROUTINE mixing_history_open

  SUBROUTINE mixing_history_close(mpi,basename)
    USE m_types,ONLY:t_mpi
    TYPE(t_mpi),INTENT(in):: mpi

    CHARACTER(len=20), OPTIONAL, INTENT(IN) :: basename

   CHARACTER(len=80):: filename
    INTEGER          :: n


    IF (iter_stored==0) RETURN ! Nothing found to be stored
    IF (mpi%isize>1) THEN
       IF (.NOT.PRESENT(basename)) THEN
          WRITE(filename,'(a,i0)') "mixing_history.",mpi%irank
       ELSE
          WRITE(filename,'(a,i0)') TRIM(basename)//"_mixing_history.",mpi%irank
       END IF
    ELSE
       IF (.NOT.PRESENT(basename)) THEN
          filename="mixing_history"
       ELSE
          filename=TRIM(basename)//"_mixing_history"
       END IF
    ENDIF
    IF (.NOT.PRESENT(basename)) THEN
      OPEN(888,file=filename,form='unformatted',status='replace')
    ELSE
      OPEN(888,file=filename,form='unformatted',status='unknown')
    END IF
    WRITE(888) iter_stored
    DO n=1,iter_stored
      call sm_store(n)%write_unformatted(888)
      call fsm_store(n)%write_unformatted(888)
    ENDDO
    CLOSE(888)
    DEALLOCATE(sm_store,fsm_store)
    iter_stored=0
  END SUBROUTINE mixing_history_close


  SUBROUTINE mixing_history(imix,maxiter,inden,outden,sm,fsm,it,nmzxyd)
    USE m_types
    implicit none
    INTEGER,INTENT(in)::imix,maxiter
    type(t_potden),intent(inout)::inden,outden
    type(t_mixvector),ALLOCATABLE::sm(:),fsm(:)
    INTEGER,INTENT(out)::it
    INTEGER,INTENT(IN) :: nmzxyd


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
    IF (.NOT.ALLOCATED(inDen%mtIm)) THEN
      CALL sm(it)%from_density(inDen,nmzxyd)
      CALL fsm(it)%from_density(outDen,nmzxyd)
    ELSE
      CALL sm(it)%from_density(inDen,nmzxyd)
      CALL fsm(it)%from_density(outDen,nmzxyd)
    END IF
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

  SUBROUTINE mixing_history_reset(mpi,basename)
    USE m_types,ONLY:t_mpi
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(in)::mpi
    CHARACTER(len=20), OPTIONAL, INTENT(IN) :: basename
    iter_stored=0
    IF (ALLOCATED(sm_store)) DEALLOCATE(sm_store,fsm_store)
    IF (mpi%irank==0) PRINT *, "Reset of history"
    IF (mpi%irank==0) THEN
       IF (.NOT.PRESENT(basename)) THEN
          CALL system('rm -f mixing_history*')
       ELSE
          CALL system('rm -f '//TRIM(basename)//'_mixing_history* '//TRIM(basename)//'mixing_history*')
       END IF
    END IF
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

  SUBROUTINE dfpt_mixing_history_reset()
    IMPLICIT NONE
    iter_stored=0
    IF (ALLOCATED(sm_store)) DEALLOCATE(sm_store,fsm_store)
END SUBROUTINE dfpt_mixing_history_reset

  SUBROUTINE mixing_history_file_count(mpi,expected_files,existing_files,l_has_history,basename)
    USE m_types,ONLY:t_mpi
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN) :: mpi
    INTEGER,INTENT(OUT) :: expected_files, existing_files
    LOGICAL,INTENT(OUT) :: l_has_history
    CHARACTER(len=20), OPTIONAL, INTENT(IN) :: basename

    CHARACTER(len=80) :: filename
    LOGICAL :: l_fileexist
    INTEGER :: rank_count

    expected_files = MERGE(mpi%isize,1,mpi%isize>1)
    existing_files = 0
    l_has_history = .FALSE.

    IF (.NOT.PRESENT(basename)) THEN
       filename = 'mixing_history'
       INQUIRE(file=filename,exist=l_fileexist)
       IF (l_fileexist) existing_files = existing_files + 1
    ELSE
       filename = TRIM(basename)//'mixing_history'
       INQUIRE(file=filename,exist=l_fileexist)
       IF (l_fileexist) existing_files = existing_files + 1
       filename = TRIM(basename)//'_mixing_history'
       INQUIRE(file=filename,exist=l_fileexist)
       IF (l_fileexist) existing_files = existing_files + 1
    END IF

    rank_count = 0
    DO
       IF (.NOT.PRESENT(basename)) THEN
          WRITE(filename,'(a,i0)') 'mixing_history.', rank_count
       ELSE
          WRITE(filename,'(a,i0)') TRIM(basename)//'_mixing_history.', rank_count
       END IF
       INQUIRE(file=filename,exist=l_fileexist)
       IF (.NOT.l_fileexist) EXIT
       existing_files = existing_files + 1
       rank_count = rank_count + 1
    END DO

    l_has_history = existing_files > 0
  END SUBROUTINE mixing_history_file_count

end MODULE m_mixing_history
