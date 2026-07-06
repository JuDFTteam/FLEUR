!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_bessel
      use m_juDFT
      IMPLICIT NONE
!-----------------------------------------------
! DESC:Calculate cylindrical Bessel functions
!      In the port to current FLEUR the hand-written downward-recurrence
!      implementation (and a buggy Intel ifport branch calling dbesj0
!      for first order) was replaced by the F2008 intrinsics.
!                 Daniel Wortmann, (08-02-19)
!-----------------------------------------------
      CONTAINS

      !<-- F: gf_bessel0(x)
      FUNCTION gf_bessel0(x)
!-----------------------------------------------
!Cylinder Bessel function of order zero
!-----------------------------------------------
      REAL, INTENT(IN) :: x
      REAL             :: gf_bessel0
      IF (x < 0.0) CALL juDFT_error("gf_bessel0: negative argument",     &
     &                              calledby="gf_bessel")
      gf_bessel0 = BESSEL_J0(x)
      END FUNCTION
      !>

      !<-- F: gf_bessel1(x)
      FUNCTION gf_bessel1(x)
!-----------------------------------------------
!Cylinder Bessel function of first order
!-----------------------------------------------
      REAL, INTENT(IN) :: x
      REAL             :: gf_bessel1
      IF (x < 0.0) CALL juDFT_error("gf_bessel1: negative argument",     &
     &                              calledby="gf_bessel")
      gf_bessel1 = BESSEL_J1(x)
      END FUNCTION
      !>

      END
