!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_spectrum
   USE m_juDFT, ONLY: juDFT_error
   IMPLICIT NONE
   PRIVATE

   REAL, PARAMETER :: xas_occ_tol = 1.0e-10
   REAL, PARAMETER :: xas_sqrt_two_pi = 2.5066282746310005024

   PUBLIC :: xas_accumulate_matrix_spectrum
   PUBLIC :: xas_gaussian_broadening

CONTAINS

   SUBROUTINE xas_accumulate_matrix_spectrum(energy_grid, eig_band, occ_band, wk, epsilon_c, matrix, eta, intensity, &
                                             broadening_type)
      ! Accumulates
      !   I(omega) += sum_band sum_mj wk*(1-occ_band(band))*|M(band,mj)|^2
      !               * g_eta(omega - (eig_band(band)-epsilon_c)).
      !
      ! occ_band(:) must contain the normalized band occupation f_nk in [0,1],
      ! without any spin-degeneracy factor. It must not be raw results%w_iks.
      ! In a FLEUR driver use, for a fixed ikpt and jsp,
      !   occ_band(:) = results%w_iks(:,ikpt,jsp) / kpts%wtkpt(ikpt)
      ! and pass
      !   wk = kpts%wtkpt(ikpt).
      !
      ! Unit convention: eig_band, epsilon_c, energy_grid, and eta must all be
      ! in the same energy unit. In FLEUR this is typically Hartree.
      REAL,              INTENT(IN)    :: energy_grid(:)
      REAL,              INTENT(IN)    :: eig_band(:)
      REAL,              INTENT(IN)    :: occ_band(:)
      REAL,              INTENT(IN)    :: wk
      REAL,              INTENT(IN)    :: epsilon_c
      COMPLEX,           INTENT(IN)    :: matrix(:, :)
      REAL,              INTENT(IN)    :: eta
      REAL,              INTENT(INOUT) :: intensity(:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: broadening_type

      INTEGER :: i_grid, i_band, i_mj
      REAL    :: transition_energy, occ_factor, strength, matrix_strength, gaussian, delta

      CALL xas_check_spectrum_inputs(energy_grid, eig_band, occ_band, matrix, eta, intensity, broadening_type)

      DO i_band = 1, SIZE(eig_band)
         occ_factor = 1.0 - occ_band(i_band)
         IF (occ_factor <= xas_occ_tol) CYCLE

         strength = 0.0
         DO i_mj = 1, SIZE(matrix, 2)
            matrix_strength = ABS(matrix(i_band, i_mj))**2
            IF (matrix_strength < TINY(matrix_strength)) CYCLE
            strength = strength + matrix_strength
         END DO
         IF (strength < TINY(strength)) CYCLE

         transition_energy = eig_band(i_band) - epsilon_c
         strength = wk*occ_factor*strength
         IF (strength < TINY(strength)) CYCLE
         DO i_grid = 1, SIZE(energy_grid)
            delta = energy_grid(i_grid) - transition_energy
            gaussian = xas_gaussian_broadening(delta, eta)
            IF (gaussian == 0.0) CYCLE
            intensity(i_grid) = intensity(i_grid) + strength*gaussian
         END DO
      END DO
   END SUBROUTINE xas_accumulate_matrix_spectrum

   REAL FUNCTION xas_gaussian_broadening(x, eta) RESULT(value)
      REAL, INTENT(IN) :: x, eta
      REAL             :: exponent

      IF (eta <= 0.0) THEN
         CALL juDFT_error("eta must be positive in xas_gaussian_broadening", calledby="m_xas_spectrum")
      END IF

      exponent = -0.5*(x/eta)**2
      IF (exponent < -700.0) THEN
         value = 0.0
      ELSE
         value = EXP(exponent)/(eta*xas_sqrt_two_pi)
      END IF
   END FUNCTION xas_gaussian_broadening

   SUBROUTINE xas_check_spectrum_inputs(energy_grid, eig_band, occ_band, matrix, eta, intensity, broadening_type)
      REAL,              INTENT(IN) :: energy_grid(:)
      REAL,              INTENT(IN) :: eig_band(:)
      REAL,              INTENT(IN) :: occ_band(:)
      COMPLEX,           INTENT(IN) :: matrix(:, :)
      REAL,              INTENT(IN) :: eta
      REAL,              INTENT(IN) :: intensity(:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: broadening_type

      IF (eta <= 0.0) THEN
         CALL juDFT_error("eta must be positive in xas_accumulate_matrix_spectrum", calledby="m_xas_spectrum")
      END IF
      IF (SIZE(eig_band) /= SIZE(occ_band)) THEN
         CALL juDFT_error("eig_band and occ_band sizes differ in xas_accumulate_matrix_spectrum", calledby="m_xas_spectrum")
      END IF
      IF (SIZE(matrix, 1) < SIZE(eig_band)) THEN
         CALL juDFT_error("matrix band dimension is too small in xas_accumulate_matrix_spectrum", calledby="m_xas_spectrum")
      END IF
      IF (SIZE(intensity) /= SIZE(energy_grid)) THEN
         CALL juDFT_error("intensity and energy_grid sizes differ in xas_accumulate_matrix_spectrum", calledby="m_xas_spectrum")
      END IF
      IF (PRESENT(broadening_type)) THEN
         SELECT CASE (TRIM(ADJUSTL(broadening_type)))
         CASE ("gaussian", "Gaussian", "GAUSSIAN", "")
         CASE DEFAULT
            CALL juDFT_error("Only Gaussian broadening is implemented for XAS", calledby="m_xas_spectrum")
         END SELECT
      END IF
   END SUBROUTINE xas_check_spectrum_inputs

END MODULE m_xas_spectrum
