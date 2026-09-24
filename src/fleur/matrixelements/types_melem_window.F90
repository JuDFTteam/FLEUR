!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_melem_window
   !> Which Bloch bands the matrix elements are built on: where the window sits in the
   !> eig file and how many bands that is.
   !>
   !> The position and the size are the same fact stated twice, and the two are needed by
   !> different code: reading the eig file wants the position, allocating a matrix over
   !> the bands wants the size. Keeping both is what lets init check that they agree.
   !>
   !> This is the whole of what the operator layer needs to know about the selection. What
   !> the bands are selected FOR -- how many Wannier functions come out of them, which
   !> energy window picked them -- is a fact about the wannierisation and lives with it, in
   !> t_melem_manifold, which extends this type.

   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE t_melem_window
      INTEGER :: num_bands = -1      !> how many bands the window holds
      INTEGER :: min_band  = -1      !> first band of the window, counted in the eig file
      INTEGER :: max_band  = -1      !> last one
   CONTAINS
      PROCEDURE :: init_window
   END TYPE t_melem_window

   PUBLIC :: t_melem_window

CONTAINS

   SUBROUTINE init_window(this, num_bands, min_band, max_band)
      CLASS(t_melem_window), INTENT(INOUT) :: this
      INTEGER,               INTENT(IN)    :: num_bands, min_band, max_band

      this%num_bands = num_bands
      this%min_band  = min_band
      this%max_band  = max_band

      !> The window is read from three attributes of which any one may be left out and
      !> derived from the other two, so all three given can disagree.
      IF (min_band < 1 .OR. max_band < min_band) &
         CALL judft_error("t_melem_window: the band window is empty or starts before the &
                          &first band", calledby="init_window")
      IF (max_band - min_band + 1 /= num_bands) &
         CALL judft_error("t_melem_window: the band window does not hold num_bands bands", &
                          calledby="init_window")
   END SUBROUTINE init_window

END MODULE m_types_melem_window
