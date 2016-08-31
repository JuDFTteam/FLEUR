!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_hybridmix

      IMPLICIT NONE

      REAL, PARAMETER       ::  amix_pbe0 = 0.25
      REAL, PARAMETER       ::  amix_hse  = 0.25
      REAL, PARAMETER       ::  omega_hse = 0.11
      REAL, PARAMETER       ::  amix_hf   = 1.00

      CONTAINS

      ! functions for variable HSE functional

      ! if a value for x is given, aMix is overwritten
      ! return the current value of aMix
      REAL FUNCTION aMix_VHSE(x) RESULT (res)

        REAL, INTENT(IN), OPTIONAL :: x
        REAL, SAVE :: aMix = aMix_HSE

        IF ( PRESENT(x) ) aMix = x
        res = aMix

      END FUNCTION aMix_VHSE

      ! if a value for x is given, omega is overwritten
      ! return the current value of omega
      REAL FUNCTION omega_VHSE(x) RESULT (res)

        REAL, INTENT(IN), OPTIONAL :: x
        REAL, SAVE :: omega = omega_HSE

        IF ( PRESENT(x) ) omega = x
        res = omega

      END FUNCTION omega_VHSE 

      END MODULE m_hybridmix
