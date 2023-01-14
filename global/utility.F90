!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_utility
  USE m_juDFT
   IMPLICIT NONE

   CONTAINS

     
   SUBROUTINE getComputerArchitectures(architectures, numArchitectures)
      IMPLICIT NONE
      INTEGER         , INTENT(OUT) :: numArchitectures
      CHARACTER(LEN=*), INTENT(OUT) :: architectures(11)
      numArchitectures = 0
      architectures = ''
#ifdef CPP_AIX
      numArchitectures = numArchitectures + 1
      architectures(numArchitectures) = 'AIX'
#endif
   END SUBROUTINE getComputerArchitectures

   SUBROUTINE getPrecision(precisionString)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(OUT) :: precisionString
      
      REAL    :: dummy
      INTEGER :: realLength
      
      dummy = 0.0
      realLength = STORAGE_SIZE(dummy)
      
      IF(realLength.EQ.64) THEN
         precisionString = 'DOUBLE'        
      ELSE IF (realLength.EQ.32) THEN
!         precisionString = 'SINGLE'
         CALL juDFT_error("You compiled with single precision, this is most probably wrong!",calledby ="dimens")
      ELSE
         WRITE(*,*) 'real length: ', realLength
         CALL juDFT_error("You compiled with unknown precision, this is most probably wrong!",calledby ="dimens")
      END IF

   END SUBROUTINE getPrecision

   SUBROUTINE getTargetStructureProperties(specifiers, numSpecifiers)
      IMPLICIT NONE
      INTEGER         , INTENT(OUT) :: numSpecifiers
      CHARACTER(LEN=*), INTENT(OUT) :: specifiers(11)
      numSpecifiers = 0
      specifiers = ''
   END SUBROUTINE getTargetStructureProperties

   SUBROUTINE getAdditionalCompilationFlags(flags, numFlags)
      IMPLICIT NONE
      INTEGER         , INTENT(OUT) :: numFlags
      CHARACTER(LEN=*), INTENT(OUT) :: flags(11)
      numFlags = 0
      flags = ''
#ifdef CPP_MPI
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_MPI'
#endif
#ifdef CPP_HDF
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_HDF'
#endif
#ifdef CPP_WANN
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_WANN'
#endif
#ifdef CPP_NOSPMVEC
      numFlags = numFlags + 1
      flags(numFlags) = '+NOSPMVEC'
#endif
#ifdef CPP_IRAPPROX
      numFlags = numFlags + 1
      flags(numFlags) = '+IRAPPROX'
#endif
   END SUBROUTINE getAdditionalCompilationFlags

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

   SUBROUTINE calcComputationBunchBounds(numBunches, iBunch, firstIndexOverall, lastIndexOverall, firstIndexBunch, lastIndexBunch)

      IMPLICIT NONE

      INTEGER, INTENT(IN)            :: numBunches, iBunch
      INTEGER, INTENT(IN)            :: firstIndexOverall, lastIndexOverall
      INTEGER, INTENT(OUT)           :: firstIndexBunch, lastIndexBunch

      INTEGER :: chunkSize, leftoverSize, length

      length = lastIndexOverall - firstIndexOverall + 1
      chunkSize = length / numBunches
      leftoverSize = MODULO(length, numBunches)
      IF (iBunch < leftoverSize) THEN
         firstIndexBunch = iBunch*(chunkSize + 1) + firstIndexOverall
         lastIndexBunch = (iBunch + 1)*(chunkSize + 1) + firstIndexOverall - 1
      ELSE
         firstIndexBunch = leftoverSize*(chunkSize + 1) + firstIndexOverall + (iBunch - leftoverSize)*chunkSize
         lastIndexBunch = (firstIndexBunch + chunkSize) - 1
      ENDIF

   END SUBROUTINE calcComputationBunchBounds

END MODULE m_utility
