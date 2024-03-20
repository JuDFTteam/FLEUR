!--------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_parallelLoop

   TYPE t_parallelLoop
      INTEGER :: overallMinIndex
      INTEGER :: overallMaxIndex
      INTEGER :: bunchMinIndex
      INTEGER :: bunchMaxIndex
   CONTAINS
      PROCEDURE :: t_parallelLoop_init
      GENERIC   :: init => t_parallelLoop_init
   END TYPE t_parallelLoop

   PUBLIC :: t_parallelLoop

CONTAINS

   SUBROUTINE t_parallelLoop_init(this, iBunch, numBunches, overallMinIndex, overallMaxIndex)

      IMPLICIT NONE

      CLASS(t_parallelLoop), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)                 :: iBunch, numBunches
      INTEGER, INTENT(IN)                 :: overallMinIndex, overallMaxIndex

      INTEGER :: chunkSize, leftoverSize, length

      this%overallMinIndex = overallMinIndex
      this%overallMaxIndex = overallMaxIndex

      length = overallMaxIndex - overallMinIndex + 1
      chunkSize = length / numBunches
      leftoverSize = MODULO(length, numBunches)
      IF (iBunch < leftoverSize) THEN
         this%bunchMinIndex = iBunch*(chunkSize + 1) + overallMinIndex
         this%bunchMaxIndex = (iBunch + 1)*(chunkSize + 1) + overallMinIndex - 1
      ELSE
         this%bunchMinIndex = leftoverSize*(chunkSize + 1) + overallMinIndex + (iBunch - leftoverSize)*chunkSize
         this%bunchMaxIndex = (this%bunchMinIndex + chunkSize) - 1
      ENDIF
   END SUBROUTINE t_parallelLoop_init

   SUBROUTINE calcNumberComputationBunches(minIndex, maxIndex, maxBunchSize, numBunches)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: minIndex
      INTEGER, INTENT(IN)  :: maxIndex
      INTEGER, INTENT(IN)  :: maxBunchSize
      INTEGER, INTENT(OUT) :: numBunches

      REAL :: length

      length = maxIndex - minIndex + 1
      numBunches = CEILING(length/(REAL(maxBunchSize)))

   END SUBROUTINE calcNumberComputationBunches

   INTEGER FUNCTION getNumberOfThreads()
      !$ use omp_lib
      IMPLICIT NONE
      INTEGER :: numThreads
      numThreads = 1
      !$omp parallel shared(numThreads)
      !$ numThreads = omp_get_num_threads()
      !$omp end parallel
      getNumberOfThreads = numThreads
   END FUNCTION getNumberOfThreads

END MODULE m_types_parallelLoop
