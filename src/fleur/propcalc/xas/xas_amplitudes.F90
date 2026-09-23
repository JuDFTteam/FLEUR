!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_amplitudes
   USE m_juDFT, ONLY: juDFT_error
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: t_xas_transition_amplitudes
   PUBLIC :: xas_emission_from_absorption

   TYPE t_xas_transition_amplitudes
      INTEGER :: n_mj = 0
      COMPLEX, ALLOCATABLE :: m(:)
      ! Maps m(:) one-to-one to physical doubled core magnetic sublevels 2*m_j.
      INTEGER, ALLOCATABLE :: twice_mj(:)
      REAL :: absM2 = 0.0
   CONTAINS
      PROCEDURE :: set_from_matrix_row => xas_transition_amplitudes_set_from_matrix_row
      PROCEDURE :: set_twice_mj => xas_transition_amplitudes_set_twice_mj
      PROCEDURE :: compute_absM2 => xas_transition_amplitudes_compute_absM2
   END TYPE t_xas_transition_amplitudes

CONTAINS

   PURE FUNCTION xas_emission_from_absorption(absorption) RESULT(emission)
      ! Central convention used by the band-to-core emission wrappers:
      ! callers pass the physical outgoing polarization unchanged.  The
      ! Hermitian conjugation of <band|D_epsilon|core> is performed here once.
      COMPLEX, INTENT(IN) :: absorption(:, :)
      COMPLEX :: emission(SIZE(absorption, 1), SIZE(absorption, 2))

      emission = CONJG(absorption)
   END FUNCTION xas_emission_from_absorption

   SUBROUTINE xas_transition_amplitudes_set_from_matrix_row(this, matrix, band, l_skip_tiny, twice_mj)
      CLASS(t_xas_transition_amplitudes), INTENT(INOUT) :: this
      COMPLEX,                            INTENT(IN)    :: matrix(:, :)
      INTEGER,                            INTENT(IN)    :: band
      LOGICAL, OPTIONAL,                  INTENT(IN)    :: l_skip_tiny
      INTEGER, OPTIONAL,                  INTENT(IN)    :: twice_mj(:)

      IF (band < 1 .OR. band > SIZE(matrix, 1)) THEN
         CALL juDFT_error("Invalid band in XAS transition amplitudes", calledby="m_xas_amplitudes")
      END IF

      this%n_mj = SIZE(matrix, 2)
      IF (ALLOCATED(this%m)) THEN
         IF (SIZE(this%m) /= this%n_mj) DEALLOCATE(this%m)
      END IF
      IF (.NOT. ALLOCATED(this%m)) ALLOCATE(this%m(this%n_mj))

      this%m = matrix(band, :)
      IF (PRESENT(twice_mj)) CALL this%set_twice_mj(twice_mj)
      CALL this%compute_absM2(l_skip_tiny)
   END SUBROUTINE xas_transition_amplitudes_set_from_matrix_row

   SUBROUTINE xas_transition_amplitudes_set_twice_mj(this, twice_mj)
      CLASS(t_xas_transition_amplitudes), INTENT(INOUT) :: this
      INTEGER,                            INTENT(IN)    :: twice_mj(:)

      IF (this%n_mj /= 0 .AND. SIZE(twice_mj) /= this%n_mj) THEN
         CALL juDFT_error("XAS twice_mj mapping does not match amplitude dimension", calledby="m_xas_amplitudes")
      END IF
      this%n_mj = SIZE(twice_mj)
      IF (ALLOCATED(this%twice_mj)) THEN
         IF (SIZE(this%twice_mj) /= this%n_mj) DEALLOCATE(this%twice_mj)
      END IF
      IF (.NOT. ALLOCATED(this%twice_mj)) ALLOCATE(this%twice_mj(this%n_mj))
      this%twice_mj = twice_mj
   END SUBROUTINE xas_transition_amplitudes_set_twice_mj

   SUBROUTINE xas_transition_amplitudes_compute_absM2(this, l_skip_tiny)
      CLASS(t_xas_transition_amplitudes), INTENT(INOUT) :: this
      LOGICAL, OPTIONAL,                  INTENT(IN)    :: l_skip_tiny

      INTEGER :: i_mj
      REAL    :: abs_m2, matrix_strength
      LOGICAL :: l_skip

      IF (.NOT. ALLOCATED(this%m)) THEN
         CALL juDFT_error("XAS transition amplitudes are not allocated", calledby="m_xas_amplitudes")
      END IF
      this%n_mj = SIZE(this%m)
      IF (ALLOCATED(this%twice_mj)) THEN
         IF (SIZE(this%twice_mj) /= this%n_mj) THEN
            CALL juDFT_error("XAS twice_mj mapping does not match allocated amplitudes", calledby="m_xas_amplitudes")
         END IF
      END IF

      l_skip = .FALSE.
      IF (PRESENT(l_skip_tiny)) l_skip = l_skip_tiny
      IF (l_skip) THEN
         abs_m2 = 0.0
         DO i_mj = 1, this%n_mj
            matrix_strength = ABS(this%m(i_mj))**2
            IF (matrix_strength < TINY(matrix_strength)) CYCLE
            abs_m2 = abs_m2 + matrix_strength
         END DO
         this%absM2 = abs_m2
      ELSE
         this%absM2 = SUM(ABS(this%m)**2)
      END IF
   END SUBROUTINE xas_transition_amplitudes_compute_absM2

END MODULE m_xas_amplitudes
