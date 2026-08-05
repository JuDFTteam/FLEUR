!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_melem_manifold
   !> Which Bloch bands the matrix elements are built on: where the window sits in the
   !> eig file, how many bands that is, how many Wannier functions come out of them, and
   !> the energy window that selected them.
   !>
   !> The position and the size are the same fact stated twice, and the code that reads
   !> the eig file needs the position while the code that interpolates needs the size.
   !> Keeping both here is what lets init check that they agree.

   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE t_melem_manifold
      INTEGER :: num_bands   = -1      !> Bloch bands entering the disentanglement
      INTEGER :: num_wann    = -1      !> Wannier functions coming out
      INTEGER :: min_band    = -1      !> first band of the window, counted in the eig file
      INTEGER :: max_band    = -1      !> last one
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

      this%num_bands   = num_bands
      this%num_wann    = num_wann
      this%min_band    = min_band
      this%max_band    = max_band
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
      !> The window is read from three attributes of which any one may be left out and
      !> derived from the other two, so all three given can disagree.
      IF (min_band < 1 .OR. max_band < min_band) &
         CALL judft_error("t_melem_manifold: the band window is empty or starts before the &
                          &first band", calledby="melem_manifold_init")
      IF (max_band - min_band + 1 /= num_bands) &
         CALL judft_error("t_melem_manifold: the band window does not hold num_bands bands", &
                          calledby="melem_manifold_init")
   END SUBROUTINE melem_manifold_init

END MODULE m_types_melem_manifold
