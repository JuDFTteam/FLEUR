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

   !> The operator names this layer knows, one list per thing that can be asked for. They
   !> are here, next to the record that carries them, because every consumer of the record
   !> already has it in scope -- and because a name that is spelled here and nowhere else
   !> should be a loud failure rather than an operator that quietly does not happen.
   !>
   !> ADDING AN OPERATOR: put its name in the right list below, then give it a branch in
   !> each place that dispatches on that list. Miss one and the run stops there naming it,
   !> instead of finishing with the operator silently absent.
   !>   interpolation : melem_run (compute) and melem_domains (name the output files)
   !>   operators_r   : melem_operators_r
   CHARACTER(LEN=20), PARAMETER, PUBLIC :: MELEM_OP_INTERP(8) = [ &
      CHARACTER(LEN=20) :: 'hamiltonian', 'spin', 'orbital', 'soc', &
                           'velocity', 'spinCurrent', 'orbitalCurrent', 'eigenstates']
   CHARACTER(LEN=20), PARAMETER, PUBLIC :: MELEM_OP_R(5) = [ &
      CHARACTER(LEN=20) :: 'hamiltonian', 'position', 'spin', 'orbital', 'spin_orbit']

   PUBLIC :: t_melem_request, melem_op_known

CONTAINS

   SUBROUTINE melem_request_init(this, l_spin, l_orbmom, l_socop, &
                                 l_operators_r, op_r_name, op_name, op_total)
      CLASS(t_melem_request), INTENT(OUT) :: this
      LOGICAL,          INTENT(IN) :: l_spin, l_orbmom, l_socop, l_operators_r
      CHARACTER(LEN=*), INTENT(IN) :: op_r_name(:)   !> the O(R) list, possibly empty
      CHARACTER(LEN=*), INTENT(IN) :: op_name(:)     !> the interpolation list, possibly empty
      INTEGER,          INTENT(IN) :: op_total(:)

      INTEGER :: iop

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

      !> A name nobody downstream recognises used to cost a line in the output and an
      !> operator that never appeared. It stops the run here instead, while it is still
      !> obvious that what is wrong is the input.
      DO iop = 1, this%n_ops
         IF (.NOT.melem_op_known(this%op_name(iop), MELEM_OP_INTERP)) &
            CALL judft_error('t_melem_request: "'//TRIM(this%op_name(iop))// &
                             '" is not an operator this layer can interpolate', &
                             hint=op_list_hint(MELEM_OP_INTERP), calledby="melem_request_init")
      END DO
      DO iop = 1, this%n_op_r
         IF (.NOT.melem_op_known(this%op_r_name(iop), MELEM_OP_R)) &
            CALL judft_error('t_melem_request: "'//TRIM(this%op_r_name(iop))// &
                             '" is not an operator this layer can export in real space', &
                             hint=op_list_hint(MELEM_OP_R), calledby="melem_request_init")
      END DO
   END SUBROUTINE melem_request_init

   !> Whether a name appears in one of the lists above.
   LOGICAL FUNCTION melem_op_known(name, list)
      CHARACTER(LEN=*), INTENT(IN) :: name
      CHARACTER(LEN=*), INTENT(IN) :: list(:)
      INTEGER :: i
      melem_op_known = .FALSE.
      DO i = 1, SIZE(list)
         IF (TRIM(name) == TRIM(list(i))) THEN
            melem_op_known = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION melem_op_known

   !> The accepted names on one line, so the error says what to write instead.
   FUNCTION op_list_hint(list) RESULT(txt)
      CHARACTER(LEN=*), INTENT(IN) :: list(:)
      CHARACTER(LEN=200) :: txt
      INTEGER :: i
      txt = 'accepted: '
      DO i = 1, SIZE(list)
         txt = TRIM(txt)//' '//TRIM(list(i))
      END DO
   END FUNCTION op_list_hint

END MODULE m_types_melem_request
