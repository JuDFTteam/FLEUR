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
   USE m_types_melem_optable, ONLY: WANNIERLIB_INTERP, WANNIERLIB_OPR, &
                                    melem_exposed_find, melem_exposed_names

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
      PROCEDURE :: init     => melem_request_init
      PROCEDURE :: has_op_r => melem_request_has_op_r
      PROCEDURE :: has_op   => melem_request_has_op
      PROCEDURE :: needs_op => melem_request_needs_op
   END TYPE t_melem_request

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
         IF (melem_exposed_find(this%op_name(iop), WANNIERLIB_INTERP) == 0) &
            CALL judft_error('t_melem_request: "'//TRIM(this%op_name(iop))// &
                             '" is not an operator this layer can interpolate', &
                             hint=melem_exposed_names(WANNIERLIB_INTERP), calledby="melem_request_init")
      END DO
      DO iop = 1, this%n_op_r
         IF (melem_exposed_find(this%op_r_name(iop), WANNIERLIB_OPR) == 0) &
            CALL judft_error('t_melem_request: "'//TRIM(this%op_r_name(iop))// &
                             '" is not an operator this layer can export in real space', &
                             hint=melem_exposed_names(WANNIERLIB_OPR), calledby="melem_request_init")
      END DO
   END SUBROUTINE melem_request_init

   !> Whether the real-space list asked for this operator. Callers used to be handed a
   !> flag per operator, which meant the same fact was stored twice and could disagree.
   LOGICAL FUNCTION melem_request_has_op_r(this, name)
      CLASS(t_melem_request), INTENT(IN) :: this
      CHARACTER(LEN=*),       INTENT(IN) :: name
      melem_request_has_op_r = .FALSE.
      IF (.NOT.this%l_operators_r)        RETURN
      IF (.NOT.ALLOCATED(this%op_r_name)) RETURN
      melem_request_has_op_r = melem_op_known(name, this%op_r_name)
   END FUNCTION melem_request_has_op_r

   !> Whether the INTERPOLATION list asked for this operator. The three summary flags
   !> cannot answer it: both lists set them, so l_spin is true for a real-space export as
   !> much as for an interpolated one -- and the two are not producible in the same spin
   !> configurations, which is the whole reason the distinction matters.
   LOGICAL FUNCTION melem_request_has_op(this, name)
      CLASS(t_melem_request), INTENT(IN) :: this
      CHARACTER(LEN=*),       INTENT(IN) :: name
      melem_request_has_op = .FALSE.
      IF (.NOT.ALLOCATED(this%op_name)) RETURN
      melem_request_has_op = melem_op_known(name, this%op_name)
   END FUNCTION melem_request_has_op

   !> Whether either list asks for something that needs this CATALOGUE entry built. The
   !> two lists spell the same operator differently, and the currents ask for one without
   !> ever naming it, so the answer is a lookup through the exposure tables rather than a
   !> string match on what the user wrote.
   LOGICAL FUNCTION melem_request_needs_op(this, opname, interp_only) RESULT(l_needed)
      CLASS(t_melem_request), INTENT(IN) :: this
      CHARACTER(LEN=*),       INTENT(IN) :: opname
      LOGICAL, OPTIONAL,      INTENT(IN) :: interp_only
      INTEGER :: i, k
      LOGICAL :: l_interp_only
      l_needed = .FALSE.
      l_interp_only = .FALSE.
      IF (PRESENT(interp_only)) l_interp_only = interp_only
      IF (ALLOCATED(this%op_name)) THEN
         DO i = 1, this%n_ops
            k = melem_exposed_find(this%op_name(i), WANNIERLIB_INTERP)
            IF (k > 0) THEN
               IF (TRIM(WANNIERLIB_INTERP(k)%operator) == TRIM(opname)) THEN
                  l_needed = .TRUE.
                  RETURN
               END IF
            END IF
         END DO
      END IF
      IF (l_interp_only) RETURN
      IF (this%l_operators_r .AND. ALLOCATED(this%op_r_name)) THEN
         DO i = 1, this%n_op_r
            k = melem_exposed_find(this%op_r_name(i), WANNIERLIB_OPR)
            IF (k > 0) THEN
               IF (TRIM(WANNIERLIB_OPR(k)%operator) == TRIM(opname)) THEN
                  l_needed = .TRUE.
                  RETURN
               END IF
            END IF
         END DO
      END IF
   END FUNCTION melem_request_needs_op

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


END MODULE m_types_melem_request
