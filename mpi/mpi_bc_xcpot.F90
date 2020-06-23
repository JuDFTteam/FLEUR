!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mpi_bc_xcpot
  use m_judft
#ifdef CPP_MPI
  use mpi 
#endif
CONTAINS
  SUBROUTINE mpi_bc_xcpot(xcpot,fmpi)

    USE m_types
    USE m_types_xcpot_libxc
    IMPLICIT NONE
    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT):: xcpot
    TYPE(t_mpi),INTENT(IN)                  :: fmpi

#ifdef CPP_MPI
    LOGICAL           :: l_relcor
    CHARACTER(len=100):: namex
    INTEGER           :: ierr,n,i(5)

    IF (fmpi%isize==1) RETURN !nothing to be done with only one PE
    !First determine type on pe0
    IF (fmpi%irank==0) THEN
       SELECT TYPE(xcpot)
       TYPE IS (t_xcpot_inbuild)
          n=1
       TYPE IS (t_xcpot_libxc)
          n=2
       CLASS DEFAULT
          CALL judft_error("Type could not be determined in mpi_bc_xcpot")
       END SELECT
    END IF
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)
    IF (fmpi%irank /= 0) THEN
       IF (ALLOCATED(xcpot)) DEALLOCATE(xcpot)
       !Now we know the types and can allocate on the other PE type dependend
       SELECT CASE(n)
       CASE(1)
          ALLOCATE(t_xcpot_inbuild::xcpot)
       CASE(2)
          ALLOCATE(t_xcpot_libxc::xcpot)
       CASE DEFAULT
          CALL judft_error("Type bcast failed in mpi_bc_xcpot")
       END SELECT
    END IF
    !Now we can do the the type dependend bc
    call xcpot%mpi_bc(fmpi%mpi_comm,0)
#endif
  END SUBROUTINE mpi_bc_xcpot
END MODULE m_mpi_bc_xcpot
