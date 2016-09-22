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
#ifdef CPP_APC
      numArchitectures = numArchitectures + 1
      architectures(numArchitectures) = 'APC'
#endif
#ifdef CPP_DEC
      numArchitectures = numArchitectures + 1
      architectures(numArchitectures) = 'DEC'
#endif
#ifdef CPP_AIX
      numArchitectures = numArchitectures + 1
      architectures(numArchitectures) = 'AIX'
#endif
#ifdef CPP_T90
      numArchitectures = numArchitectures + 1
#ifdef CPP_MPI
      architectures(numArchitectures) = 'T3E'
#else
      architectures(numArchitectures) = 'T90'
#endif
#endif
   END SUBROUTINE getComputerArchitectures

   SUBROUTINE getPrecision(precisionString)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(OUT) :: precisionString
#ifdef CPP_DOUBLE                
      precisionString = 'DOUBLE'        
#else
      precisionString = 'SINGLE'
#ifndef CPP_T90
      CALL juDFT_error(" define CPP_DOUBLE on non-Cray architectures!",calledby ="dimens")
#endif
#endif
   END SUBROUTINE getPrecision

   SUBROUTINE getTargetStructureProperties(specifiers, numSpecifiers)
      IMPLICIT NONE
      INTEGER         , INTENT(OUT) :: numSpecifiers
      CHARACTER(LEN=*), INTENT(OUT) :: specifiers(11)
      numSpecifiers = 0
      specifiers = ''
      numSpecifiers = numSpecifiers + 1
#ifdef CPP_INVERSION
      specifiers(numSpecifiers) = "InversionSymmetry"
#else
      specifiers(numSpecifiers) = "noInversionSymmetry"
#endif
      numSpecifiers = numSpecifiers + 1
#ifdef CPP_SOC
      specifiers(numSpecifiers) = "SpinOrbitCoupling"
#else
      specifiers(numSpecifiers) = "noSpinOrbitCoupling"
#endif
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
#ifdef CPP_APW
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_APW'
#endif
#ifdef CPP_CORE
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_CORE'
#endif
#ifdef CPP_HTML
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_HTML'
#endif
#ifdef CPP_HDF
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_HDF'
#endif
#ifdef CPP_WANN
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_WANN'
#endif
#ifdef CPP_600
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_600'
#endif
#ifdef CPP_GF
      numFlags = numFlags + 1
      flags(numFlags) = 'CPP_GF'
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
