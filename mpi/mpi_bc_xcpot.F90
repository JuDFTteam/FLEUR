!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mpi_bc_xcpot
CONTAINS
  SUBROUTINE mpi_bc_xcpot(xcpot,mpi)
    USE m_types
    IMPLICIT NONE
    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT):: xcpot
    TYPE(t_mpi),INTENT(IN)                  :: mpi

#ifdef CPP_MPI
    LOGICAL         :: l_relcor
    CHARACTER(len=4):: namex
    INTEGER         :: ierr,n
    INCLUDE 'mpif.h'

    IF (mpi%isize==1) RETURN !nothing to be done with only one PE
    !First determine type on pe0
    IF (mpi%irank==0) THEN
       SELECT TYPE(xcpot)
       TYPE IS (t_xcpot_inbuild)
          n=1
       CLASS DEFAULT
          CALL judft_error("Type could not be determined in mpi_bc_xcpot")
       END SELECT
    END IF
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF (mpi%irank.NE.0) THEN
       IF (ALLOCATED(xcpot)) DEALLOCATE(xcpot)
       !Now we know the types and can allocate on the other PE type dependend
       SELECT CASE(n)
       CASE(1)
          ALLOCATE(t_xcpot_inbuild::xcpot)
       CASE DEFAULT
          CALL judft_error("Type bcast failed in mpi_bc_xcpot")
       END SELECT
    END IF
    !Now we can do the the type dependend bc
    SELECT TYPE(xcpot)
    TYPE IS (t_xcpot_inbuild)
       IF (mpi%irank==0) THEN
          namex=xcpot%get_name()
          l_relcor=xcpot%DATA%krla==1
          n=SIZE(xcpot%lda_atom)
       ENDIF
       CALL MPI_BCAST(namex,4,MPI_CHARACTER,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(l_relcor,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
       CALL MPI_BCAST(n,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
       IF (mpi%irank.NE.0)  CALL xcpot%init(namex,l_relcor,n)
       CALL MPI_BCAST(xcpot%lda_atom,n,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    END SELECT
#endif
  END SUBROUTINE mpi_bc_xcpot
END MODULE m_mpi_bc_xcpot
