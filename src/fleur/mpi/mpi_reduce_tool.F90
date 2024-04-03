!--------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mpi_reduce_tool

   USE m_juDFT
#ifdef CPP_MPI
   USE mpi
#endif

   IMPLICIT NONE
   PRIVATE
  
   INTERFACE mpi_sum_reduce
      MODULE PROCEDURE mpi_sum_reduce_int1, mpi_sum_reduce_int2, mpi_sum_reduce_int3
      MODULE PROCEDURE mpi_sum_reduce_real1, mpi_sum_reduce_real2, mpi_sum_reduce_real3
      MODULE PROCEDURE mpi_sum_reduce_complex1, mpi_sum_reduce_complex2, mpi_sum_reduce_complex3
   END INTERFACE mpi_sum_reduce

   INTERFACE mpi_lor_reduce
      MODULE PROCEDURE mpi_lor_reduce_bool1, mpi_lor_reduce_bool2, mpi_lor_reduce_bool3
   END INTERFACE mpi_lor_reduce
   
   PUBLIC :: mpi_sum_reduce, mpi_lor_reduce
   
   CONTAINS
   
   ! INTEGER SUBROUTINES:
   
   SUBROUTINE mpi_sum_reduce_int1(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: sourceArray(:)
      INTEGER, INTENT(INOUT) :: targetArray(:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_INTEGER, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:) = sourceArray(:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_int1
   
   SUBROUTINE mpi_sum_reduce_int2(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: sourceArray(:,:)
      INTEGER, INTENT(INOUT) :: targetArray(:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_INTEGER, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:,:) = sourceArray(:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_int2
   
   SUBROUTINE mpi_sum_reduce_int3(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: sourceArray(:,:,:)
      INTEGER, INTENT(INOUT) :: targetArray(:,:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_INTEGER, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:,:,:) = sourceArray(:,:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_int3
   
   ! REAL SUBROUTINES:
   
   SUBROUTINE mpi_sum_reduce_real1(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      REAL,    INTENT(IN)    :: sourceArray(:)
      REAL,    INTENT(INOUT) :: targetArray(:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:) = sourceArray(:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_real1

   SUBROUTINE mpi_sum_reduce_real2(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      REAL,    INTENT(IN)    :: sourceArray(:,:)
      REAL,    INTENT(INOUT) :: targetArray(:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:,:) = sourceArray(:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_real2

   SUBROUTINE mpi_sum_reduce_real3(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      REAL,    INTENT(IN)    :: sourceArray(:,:,:)
      REAL,    INTENT(INOUT) :: targetArray(:,:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:,:,:) = sourceArray(:,:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_real3

   ! COMPLEX SUBROUTINES:

   SUBROUTINE mpi_sum_reduce_complex1(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      COMPLEX, INTENT(IN)    :: sourceArray(:)
      COMPLEX, INTENT(INOUT) :: targetArray(:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:) = sourceArray(:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_complex1

   SUBROUTINE mpi_sum_reduce_complex2(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      COMPLEX, INTENT(IN)    :: sourceArray(:,:)
      COMPLEX, INTENT(INOUT) :: targetArray(:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:,:) = sourceArray(:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_complex2

   SUBROUTINE mpi_sum_reduce_complex3(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      COMPLEX, INTENT(IN)    :: sourceArray(:,:,:)
      COMPLEX, INTENT(INOUT) :: targetArray(:,:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, mpi_comm, ierr)
#else
      targetArray(:,:,:) = sourceArray(:,:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_sum_reduce_complex3

   ! LOGICAL SUBROUTINES :

   SUBROUTINE mpi_lor_reduce_bool1(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: sourceArray(:)
      LOGICAL, INTENT(INOUT) :: targetArray(:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_LOGICAL, MPI_LOR, 0, mpi_comm, ierr)
#else
      targetArray(:) = sourceArray(:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_lor_reduce_bool1

   SUBROUTINE mpi_lor_reduce_bool2(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: sourceArray(:,:)
      LOGICAL, INTENT(INOUT) :: targetArray(:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_LOGICAL, MPI_LOR, 0, mpi_comm, ierr)
#else
      targetArray(:,:) = sourceArray(:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_lor_reduce_bool2

   SUBROUTINE mpi_lor_reduce_bool3(sourceArray, targetArray, mpi_comm)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: sourceArray(:,:,:)
      LOGICAL, INTENT(INOUT) :: targetArray(:,:,:)
      INTEGER, INTENT(IN)    :: mpi_comm

      INTEGER :: ierr=0
      INTEGER :: length

      length = SIZE(sourceArray)
      IF(length.NE.SIZE(targetArray)) CALL judft_error("MPI_REDUCE failed: Array size mismatch.")

#ifdef CPP_MPI
      CALL MPI_REDUCE(sourceArray, targetArray, length, MPI_LOGICAL, MPI_LOR, 0, mpi_comm, ierr)
#else
      targetArray(:,:,:) = sourceArray(:,:,:)
#endif
      IF (ierr.NE.0) CALL judft_error("MPI_REDUCE failed")
   END SUBROUTINE mpi_lor_reduce_bool3


END MODULE m_mpi_reduce_tool
