!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_io
   USE m_juDFT, ONLY: juDFT_error
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: xas_write_spectrum_text

CONTAINS

   SUBROUTINE xas_write_spectrum_text(filename, energy_grid, intensity, energy_unit)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      REAL,             INTENT(IN) :: energy_grid(:)
      REAL,             INTENT(IN) :: intensity(:)
      CHARACTER(LEN=*), INTENT(IN) :: energy_unit

      INTEGER :: io_unit, i_grid, io_status
      CHARACTER(LEN=256) :: io_message

      IF (SIZE(energy_grid) /= SIZE(intensity)) THEN
         CALL juDFT_error("energy_grid and intensity sizes differ in xas_write_spectrum_text", calledby="m_xas_io")
      END IF

      OPEN(NEWUNIT=io_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write", FORM="formatted", &
           IOSTAT=io_status, IOMSG=io_message)
      IF (io_status /= 0) THEN
         CALL juDFT_error("Cannot open XAS spectrum text file: "//TRIM(io_message), calledby="m_xas_io")
      END IF

      WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# XAS spectrum"
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# energy_unit: "//TRIM(energy_unit)
      IF (io_status == 0) WRITE(io_unit, '(a,i0)', IOSTAT=io_status, IOMSG=io_message) "# n_grid: ", SIZE(energy_grid)
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# intensity: arbitrary units"
      IF (io_status == 0) WRITE(io_unit, '(a)', IOSTAT=io_status, IOMSG=io_message) "# columns: energy intensity"

      DO i_grid = 1, SIZE(energy_grid)
         IF (io_status /= 0) EXIT
         WRITE(io_unit, '(2es24.16e3)', IOSTAT=io_status, IOMSG=io_message) energy_grid(i_grid), intensity(i_grid)
      END DO

      CLOSE(io_unit)

      IF (io_status /= 0) THEN
         CALL juDFT_error("Cannot write XAS spectrum text file: "//TRIM(io_message), calledby="m_xas_io")
      END IF
   END SUBROUTINE xas_write_spectrum_text

END MODULE m_xas_io
