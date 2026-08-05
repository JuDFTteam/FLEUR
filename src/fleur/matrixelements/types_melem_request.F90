!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_melem_request
   !> Which matrix elements were asked for.
   !>
   !> There are two independent lists, because the two things one can ask for are not the
   !> same: an operator can be exported as O(R) in real space, or interpolated onto an
   !> output domain, or both, and the two lists spell some operators differently. The
   !> three flags summarise both lists for the coarse-mesh pass, which has to produce a
   !> Bloch matrix before either list can be served and does not care which one asked.

   USE m_judft

   IMPLICIT NONE
   PRIVATE

   TYPE t_melem_request
      !> Coarse-mesh Bloch matrices to produce. True if either list needs them.
      LOGICAL :: l_spin   = .FALSE.
      LOGICAL :: l_orbmom = .FALSE.
      LOGICAL :: l_socop  = .FALSE.
      !> Real-space O(R) export.
      LOGICAL :: l_operators_r = .FALSE.
      INTEGER :: n_op_r = 0
      CHARACTER(LEN=20), ALLOCATABLE :: op_r_name(:)
      !> Band interpolation. op_total selects the site-summed projection.
      INTEGER :: n_ops = 0
      CHARACTER(LEN=20), ALLOCATABLE :: op_name(:)
      INTEGER,           ALLOCATABLE :: op_total(:)
   CONTAINS
      PROCEDURE :: init => melem_request_init
   END TYPE t_melem_request

   PUBLIC :: t_melem_request

CONTAINS

   SUBROUTINE melem_request_init(this, l_spin, l_orbmom, l_socop, &
                                 l_operators_r, op_r_name, op_name, op_total)
      CLASS(t_melem_request), INTENT(OUT) :: this
      LOGICAL,          INTENT(IN) :: l_spin, l_orbmom, l_socop, l_operators_r
      CHARACTER(LEN=*), INTENT(IN) :: op_r_name(:)   !> the O(R) list, possibly empty
      CHARACTER(LEN=*), INTENT(IN) :: op_name(:)     !> the interpolation list, possibly empty
      INTEGER,          INTENT(IN) :: op_total(:)

      this%l_spin   = l_spin
      this%l_orbmom = l_orbmom
      this%l_socop  = l_socop

      !> The counts are the extents of the lists rather than a second statement about
      !> them, so a list and its length cannot come apart.
      this%l_operators_r = l_operators_r
      this%n_op_r        = SIZE(op_r_name)
      this%op_r_name     = op_r_name
      this%n_ops         = SIZE(op_name)
      this%op_name       = op_name
      this%op_total      = op_total

      IF (SIZE(op_total) /= SIZE(op_name)) &
         CALL judft_error("t_melem_request: every interpolated operator needs its own total flag", &
                          calledby="melem_request_init")
   END SUBROUTINE melem_request_init

END MODULE m_types_melem_request
