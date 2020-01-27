!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_utility

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
#ifdef CPP_DOUBLE                
      precisionString = 'DOUBLE'        
#else
      precisionString = 'SINGLE'
      CALL juDFT_warn("You compiled with single precision, this is most probably wrong!",calledby ="dimens")
#endif
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

END MODULE m_utility
