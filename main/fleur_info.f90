!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This module realizes the Fleur info mode: It prints out some information
!!! about the charge density file and then ends the program.
!!!
!!!                             GM'17
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_fleur_info

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE fleur_info()

      USE m_juDFT
      USE m_cdn_io

      IMPLICIT NONE

      LOGICAL       :: l_exist

      IF (.NOT.juDFT_was_argument("-info")) RETURN

      WRITE(*,*) 'Fleur info mode'
      WRITE(*,*) '================================================='
      WRITE(*,*) ''
      CALL printDensityFileInfo()
      WRITE(*,*) '================================================='
      CALL juDFT_error("Fleur info output completed")
   END SUBROUTINE fleur_info

END MODULE m_fleur_info
