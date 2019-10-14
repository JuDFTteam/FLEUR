!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_bc_tool
  USE m_judft
  IMPLICIT NONE
  PRIVATE
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif
  !This interface is used to broadcast data. On the recieving PE the data-array is first allocated to
  !have the same shape as the one on irank
  INTERFACE mpi_bc
     MODULE PROCEDURE  mpi_bc_int,mpi_bc_int1,mpi_bc_int2,mpi_bc_int3,mpi_bc_int4,mpi_bc_int5
     MODULE PROCEDURE  mpi_bc_real_3,mpi_bc_real_fixed2,mpi_bc_real,mpi_bc_real1,mpi_bc_real2,mpi_bc_real3,mpi_bc_real4,mpi_bc_real5
     MODULE PROCEDURE  mpi_bc_complex,mpi_bc_complex1,mpi_bc_complex2,mpi_bc_complex3,mpi_bc_complex4,mpi_bc_complex5
     MODULE PROCEDURE  mpi_bc_logical,mpi_bc_logical1,mpi_bc_logical2
     MODULE PROCEDURE  mpi_bc_character_fixed1
  END INTERFACE mpi_bc
  PUBLIC :: mpi_bc
CONTAINS
  SUBROUTINE mpi_bc_logical(i,irank,mpi_comm)
    IMPLICIT NONE
    LOGICAL,INTENT(INOUT):: i
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0

#ifdef CPP_MPI
    CALL MPI_BCAST(i,1,MPI_LOGICAL,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_logical

  SUBROUTINE mpi_bc_logical1(i,irank,mpi_comm)
    IMPLICIT NONE
    LOGICAL,ALLOCATABLE,INTENT(INOUT) :: i(:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(1),iup(1),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,1,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,1,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_LOGICAL,irank,mpi_comm,ierr)

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
#endif
  END SUBROUTINE mpi_bc_logical1

  SUBROUTINE mpi_bc_logical2(i,irank,mpi_comm)
    IMPLICIT NONE
    LOGICAL,ALLOCATABLE,INTENT(INOUT) :: i(:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(2),iup(2),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,2,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,2,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1),ilow(2):iup(2)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_LOGICAL,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_logical2

  SUBROUTINE mpi_bc_int(i,irank,mpi_comm)
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: i
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0

#ifdef CPP_MPI
    CALL MPI_BCAST(i,1,MPI_INTEGER,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_int

  SUBROUTINE mpi_bc_int1(i,irank,mpi_comm)
    IMPLICIT NONE
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: i(:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(1),iup(1),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,1,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,1,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,irank,mpi_comm,ierr)

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
#endif
  END SUBROUTINE mpi_bc_int1

  SUBROUTINE mpi_bc_int2(i,irank,mpi_comm)
    IMPLICIT NONE
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: i(:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(2),iup(2),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,2,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,2,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1),ilow(2):iup(2)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_int2



  SUBROUTINE mpi_bc_int3(i,irank,mpi_comm)
    IMPLICIT NONE
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: i(:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(3),iup(3),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,3,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,3,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_int3

  SUBROUTINE mpi_bc_int4(i,irank,mpi_comm)
    IMPLICIT NONE
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: i(:,:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(4),iup(4),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,4,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,4,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (allocated(i).and.myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3),ilow(4):iup(4)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_int4

  SUBROUTINE mpi_bc_int5(i,irank,mpi_comm)
    IMPLICIT NONE
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: i(:,:,:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(5),iup(5),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (allocated(i).and.myrank==irank) THEN
       ilow=LBOUND(i)
       iup=UBOUND(i)
    END IF
    CALL MPI_BCAST(ilow,5,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,5,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(i)) DEALLOCATE(i)
       ALLOCATE(i(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3),ilow(4):iup(4),ilow(5):iup(5)))
    ENDIF

    CALL MPI_BCAST(i,SIZE(i),MPI_INTEGER,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_int5

  !
  ! now the same for reals
  !

  SUBROUTINE mpi_bc_real_fixed2(irank,mpi_comm,r2)!Special routine for non-allocatable 2d arrays
    IMPLICIT NONE
    REAL,INTENT(INOUT)   :: r2(:,:)
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0
#ifdef CPP_MPI

    CALL MPI_BCAST(r2,size(r2),MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real_fixed2

  SUBROUTINE mpi_bc_real_3(irank,mpi_comm,r)!Special routine for non-allocatable 3x3 arrays
    IMPLICIT NONE
    REAL,INTENT(INOUT)   :: r(3)
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0
#ifdef CPP_MPI

    CALL MPI_BCAST(r,3,MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real_3

  SUBROUTINE mpi_bc_real(r,irank,mpi_comm)
    IMPLICIT NONE
    REAL,INTENT(INOUT)   :: r
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0
#ifdef CPP_MPI

    CALL MPI_BCAST(r,1,MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real


  SUBROUTINE mpi_bc_real1(r,irank,mpi_comm)
    IMPLICIT NONE
    REAL   ,ALLOCATABLE,INTENT(INOUT) :: r(:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(1),iup(1),myrank
    ilow=0;iup=0
#ifdef CPP_MPI
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(r).and.myrank==irank) THEN
       ilow=LBOUND(r)
       iup=UBOUND(r)
    END IF
    CALL MPI_BCAST(ilow,1,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,1,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(r)) DEALLOCATE(r)
       ALLOCATE(r(ilow(1):iup(1)))
    ENDIF

    CALL MPI_BCAST(r,SIZE(r),MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real1

  SUBROUTINE mpi_bc_real2(r,irank,mpi_comm)
    IMPLICIT NONE
    REAL   ,ALLOCATABLE,INTENT(INOUT) :: r(:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(2),iup(2),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(r).and.myrank==irank) THEN
       ilow=LBOUND(r)
       iup=UBOUND(r)
    END IF
    CALL MPI_BCAST(ilow,2,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,2,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(r)) DEALLOCATE(r)
       ALLOCATE(r(ilow(1):iup(1),ilow(2):iup(2)))
    ENDIF

    CALL MPI_BCAST(r,SIZE(r),MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real2

  SUBROUTINE mpi_bc_real3(r,irank,mpi_comm)
    IMPLICIT NONE
    REAL   ,ALLOCATABLE,INTENT(INOUT) :: r(:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(3),iup(3),myrank

#ifdef CPP_MPI
    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(r).and.myrank==irank) THEN
       ilow=LBOUND(r)
       iup=UBOUND(r)
    END IF
    CALL MPI_BCAST(ilow,3,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,3,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(r)) DEALLOCATE(r)
       ALLOCATE(r(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3)))
    ENDIF

    CALL MPI_BCAST(r,SIZE(r),MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real3

  SUBROUTINE mpi_bc_real4(r,irank,mpi_comm)
    IMPLICIT NONE
    REAL   ,ALLOCATABLE,INTENT(INOUT) :: r(:,:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(4),iup(4),myrank
#ifdef CPP_MPI

    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(r).and.myrank==irank) THEN
       ilow=LBOUND(r)
       iup=UBOUND(r)
    END IF
    CALL MPI_BCAST(ilow,4,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,4,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(r)) DEALLOCATE(r)
       ALLOCATE(r(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3),ilow(4):iup(4)))
    ENDIF

    CALL MPI_BCAST(r,SIZE(r),MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)

#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real4

  SUBROUTINE mpi_bc_real5(r,irank,mpi_comm)
    IMPLICIT NONE
    REAL   ,ALLOCATABLE,INTENT(INOUT) :: r(:,:,:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(5),iup(5),myrank
#ifdef CPP_MPI

    ilow=0;iup=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(r).and.myrank==irank) THEN
       ilow=LBOUND(r)
       iup=UBOUND(r)
    END IF
    CALL MPI_BCAST(ilow,5,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,5,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(r)) DEALLOCATE(r)
       ALLOCATE(r(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3),ilow(4):iup(4),ilow(5):iup(5)))
    ENDIF

    CALL MPI_BCAST(r,SIZE(r),MPI_DOUBLE_PRECISION,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_real5

  !
  ! And Complex!!
  !

  SUBROUTINE mpi_bc_complex(c,irank,mpi_comm)
    IMPLICIT NONE
    COMPLEX,INTENT(INOUT)   :: c
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0
#ifdef CPP_MPI

    CALL MPI_BCAST(c,1,MPI_DOUBLE_COMPLEX,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_complex

  SUBROUTINE mpi_bc_complex1(c,irank,mpi_comm)
    IMPLICIT NONE
    COMPLEX,ALLOCATABLE,INTENT(INOUT) :: c(:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(1),iup(1),myrank

#ifdef CPP_MPI
    iup=0;ilow=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(c).and.myrank==irank) THEN
       ilow=LBOUND(c)
       iup=UBOUND(c)
    END IF
    CALL MPI_BCAST(ilow,1,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,1,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(c)) DEALLOCATE(c)
       ALLOCATE(c(ilow(1):iup(1)))
    ENDIF

    CALL MPI_BCAST(c,SIZE(c),MPI_DOUBLE_COMPLEX,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_complex1

  SUBROUTINE mpi_bc_complex2(c,irank,mpi_comm)
    IMPLICIT NONE
    COMPLEX,ALLOCATABLE,INTENT(INOUT) :: c(:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(2),iup(2),myrank
#ifdef CPP_MPI

    iup=0;ilow=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(c).and.myrank==irank) THEN
       ilow=LBOUND(c)
       iup=UBOUND(c)
    END IF
    CALL MPI_BCAST(ilow,2,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,2,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(c)) DEALLOCATE(c)
       ALLOCATE(c(ilow(1):iup(1),ilow(2):iup(2)))
    ENDIF

    CALL MPI_BCAST(c,SIZE(c),MPI_DOUBLE_COMPLEX,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_complex2

  SUBROUTINE mpi_bc_complex3(c,irank,mpi_comm)
    IMPLICIT NONE
    COMPLEX,ALLOCATABLE,INTENT(INOUT) :: c(:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(3),iup(3),myrank
#ifdef CPP_MPI

    iup=0;ilow=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(c).and.myrank==irank) THEN
       ilow=LBOUND(c)
       iup=UBOUND(c)
    END IF
    CALL MPI_BCAST(ilow,3,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,3,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(c)) DEALLOCATE(c)
       ALLOCATE(c(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3)))
    ENDIF

    CALL MPI_BCAST(c,SIZE(c),MPI_DOUBLE_COMPLEX,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_complex3

  SUBROUTINE mpi_bc_complex4(c,irank,mpi_comm)
    IMPLICIT NONE
    COMPLEX,ALLOCATABLE,INTENT(INOUT) :: c(:,:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(4),iup(4),myrank
#ifdef CPP_MPI

    iup=0;ilow=0

    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(c).and.myrank==irank) THEN
       ilow=LBOUND(c)
       iup=UBOUND(c)
    END IF
    CALL MPI_BCAST(ilow,4,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,4,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(c)) DEALLOCATE(c)
       ALLOCATE(c(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3),ilow(4):iup(4)))
    ENDIF

    CALL MPI_BCAST(c,SIZE(c),MPI_DOUBLE_COMPLEX,irank,mpi_comm,ierr)
#endif

    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_complex4

  SUBROUTINE mpi_bc_complex5(c,irank,mpi_comm)
    IMPLICIT NONE
    COMPLEX,ALLOCATABLE,INTENT(INOUT) :: c(:,:,:,:,:)
    INTEGER,INTENT(IN)                :: irank,mpi_comm

    INTEGER:: ierr=0,ilow(5),iup(5),myrank
#ifdef CPP_MPI

    iup=0;ilow=0
    CALL MPI_COMM_RANK(mpi_comm,myrank,ierr)
    IF (ALLOCATED(c).and.myrank==irank) THEN
       ilow=LBOUND(c)
       iup=UBOUND(c)
    END IF
    CALL MPI_BCAST(ilow,5,MPI_INTEGER,0,mpi_comm,ierr)
    CALL MPI_BCAST(iup,5,MPI_INTEGER,0,mpi_comm,ierr)
    If (All(Ilow==0).And.All(Iup==0)) then
 print *,"Empty mpi_bc" ; return
 endif
    IF (myrank.NE.irank) THEN
       IF (ALLOCATED(c)) DEALLOCATE(c)
       ALLOCATE(c(ilow(1):iup(1),ilow(2):iup(2),ilow(3):iup(3),ilow(4):iup(4),ilow(5):iup(5)))
    ENDIF

    CALL MPI_BCAST(c,SIZE(c),MPI_DOUBLE_COMPLEX,irank,mpi_comm,ierr)

#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_complex5


  SUBROUTINE mpi_bc_character_fixed1(irank,mpi_comm,c)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(INOUT)   :: c
    INTEGER,INTENT(IN)   :: mpi_comm,irank

    INTEGER:: ierr=0
#ifdef CPP_MPI

    CALL MPI_BCAST(c,len(c),MPI_CHARACTER,irank,mpi_comm,ierr)
#endif
    IF (ierr.NE.0) CALL judft_error("MPI_BCAST failed")
  END SUBROUTINE mpi_bc_character_fixed1


END MODULE m_mpi_bc_tool
