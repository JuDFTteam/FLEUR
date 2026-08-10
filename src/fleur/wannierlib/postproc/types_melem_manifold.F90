!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_melem_manifold
   !> The band window plus what the wannierisation makes of it: how many Wannier functions
   !> come out of those bands, and the energy window that selected them.
   !>
   !> It EXTENDS t_melem_window rather than repeating it, so a consumer that only needs the
   !> selection takes the parent and one that needs the Wannier side takes this. That is
   !> what keeps num_wann and the disentanglement window out of the operator layer, which
   !> has no use for either: the coarse pass is handed the window and never learns what the
   !> bands were selected for.

   USE m_types_melem_window
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_melem_window) :: t_melem_manifold
      INTEGER :: num_wann    = -1      !> Wannier functions coming out
      REAL    :: dis_win_min = 0.0     !> lower edge of the energy window
      REAL    :: dis_win_max = 0.0     !> upper edge
   CONTAINS
      PROCEDURE :: init => melem_manifold_init
   END TYPE t_melem_manifold

   PUBLIC :: t_melem_manifold

CONTAINS

   SUBROUTINE melem_manifold_init(this, num_bands, num_wann, dis_win_min, dis_win_max, &
                                  min_band, max_band)
      CLASS(t_melem_manifold), INTENT(OUT) :: this
      INTEGER,                 INTENT(IN)  :: num_bands, num_wann
      REAL,                    INTENT(IN)  :: dis_win_min, dis_win_max
      INTEGER,                 INTENT(IN)  :: min_band, max_band

      CALL this%init_window(num_bands, min_band, max_band)

      this%num_wann    = num_wann
      this%dis_win_min = dis_win_min
      this%dis_win_max = dis_win_max

      !> A window given the wrong way round selects no band at all, and the interpolation
      !> would return zeros rather than fail.
      IF (num_wann < 1 .OR. num_bands < num_wann) &
         CALL judft_error("t_melem_manifold: there must be at least one Wannier function and &
                          &no more of them than bands", calledby="melem_manifold_init")
      IF (dis_win_max < dis_win_min) &
         CALL judft_error("t_melem_manifold: the energy window is inverted", &
                          calledby="melem_manifold_init")
   END SUBROUTINE melem_manifold_init

END MODULE m_types_melem_manifold
