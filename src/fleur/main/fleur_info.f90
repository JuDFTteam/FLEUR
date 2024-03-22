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

   SUBROUTINE fleur_info(kpts)

      USE m_juDFT
      USE m_cdn_io
      USE m_setupMPI
      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_kpts), INTENT(IN)     :: kpts

      LOGICAL       :: l_exist

      WRITE(*,*) ''
      WRITE(*,'(a)') ' ========== k-point set info =========='
      WRITE(*,'(2a)') ' Selected k-point list: ', TRIM(ADJUSTL(kpts%kptsName))
      WRITE(*,'(2a)') ' k-point list type: ', TRIM(ADJUSTL(kptsKindString_consts(kpts%kptsKind)))
      IF(kpts%kptsKind.EQ.KPTS_KIND_MESH) WRITE(*,'(a,i0,a,i0,a,i0)') ' ', kpts%nkpt3(1), ' x ', kpts%nkpt3(2), ' x ', kpts%nkpt3(3)
      WRITE(*,'(a,i0)') ' Number of k points: ', kpts%nkpt
      WRITE(*,*) ''

      IF (.NOT.juDFT_was_argument("-info")) RETURN

      WRITE(*,*) 'Fleur info mode'
      WRITE(*,*) '================================================='
      WRITE(*,*) ''
      CALL priv_dist_info(kpts%nkpt)
      WRITE(*,*) ''
      CALL printDensityFileInfo()
      WRITE(*,*) '================================================='
      CALL juDFT_end("Fleur info output completed")
   END SUBROUTINE fleur_info

END MODULE m_fleur_info
