!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_amplitudes
   USE m_juDFT, ONLY: juDFT_error
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: t_xas_transition_amplitudes

   TYPE t_xas_transition_amplitudes
      INTEGER :: n_mj = 0
      COMPLEX, ALLOCATABLE :: m(:)
      REAL :: absM2 = 0.0
   CONTAINS
      PROCEDURE :: set_from_matrix_row => xas_transition_amplitudes_set_from_matrix_row
      PROCEDURE :: compute_absM2 => xas_transition_amplitudes_compute_absM2
   END TYPE t_xas_transition_amplitudes

CONTAINS

   SUBROUTINE xas_transition_amplitudes_set_from_matrix_row(this, matrix, band, l_skip_tiny)
      CLASS(t_xas_transition_amplitudes), INTENT(INOUT) :: this
      COMPLEX,                            INTENT(IN)    :: matrix(:, :)
      INTEGER,                            INTENT(IN)    :: band
      LOGICAL, OPTIONAL,                  INTENT(IN)    :: l_skip_tiny

      IF (band < 1 .OR. band > SIZE(matrix, 1)) THEN
         CALL juDFT_error("Invalid band in XAS transition amplitudes", calledby="m_xas_amplitudes")
      END IF

      this%n_mj = SIZE(matrix, 2)
      IF (ALLOCATED(this%m)) THEN
         IF (SIZE(this%m) /= this%n_mj) DEALLOCATE(this%m)
      END IF
      IF (.NOT. ALLOCATED(this%m)) ALLOCATE(this%m(this%n_mj))

      this%m = matrix(band, :)
      CALL this%compute_absM2(l_skip_tiny)
   END SUBROUTINE xas_transition_amplitudes_set_from_matrix_row

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
