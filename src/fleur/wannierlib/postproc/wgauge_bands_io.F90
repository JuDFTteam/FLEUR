!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>  The single place that knows the on-disk layout of the interpolated band files, the
!>  way m_wgauge_io is the single place that knows the layout of the real-space operators.
!>
!>  One row per point of the output domain:
!>
!>     kdist   [ E_1  c_1 .. c_ncomp ]  [ E_2  c_1 .. c_ncomp ]  ...
!>
!>  the energy always in f14.8 and the components in whatever the quantity needs -- a
!>  velocity and a Berry curvature do not share a magnitude. A file with no components,
!>  which is the plain band structure, is the same row with ncomp = 0.
MODULE m_wgauge_bands_io
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wgauge_bands_open, wgauge_bands_row

CONTAINS

   !> Open <fname>.dat and put its one-line header in. The '.dat' is added here so that no
   !> caller has to remember it.
   SUBROUTINE wgauge_bands_open(iu, fname, header)
      INTEGER, INTENT(OUT) :: iu
      CHARACTER(LEN=*), INTENT(IN) :: fname, header
      INTEGER :: ierr

      OPEN(newunit=iu, file=TRIM(fname)//'.dat', status='replace', iostat=ierr)
      IF (ierr /= 0) CALL juDFT_error('cannot write '//TRIM(fname)//'.dat', &
                                      calledby='wgauge_bands_open')
      WRITE(iu, '(a)') TRIM(header)
   END SUBROUTINE wgauge_bands_open

   !> One row. evals is already in the unit the file advertises; comp is (ncomp, num_wann)
   !> and absent for a file that carries the energies alone.
   SUBROUTINE wgauge_bands_row(iu, kdist, evals, comp, fmt_comp)
      INTEGER, INTENT(IN) :: iu
      REAL, INTENT(IN) :: kdist
      REAL, INTENT(IN) :: evals(:)
      REAL, INTENT(IN), OPTIONAL :: comp(:, :)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt_comp
      INTEGER :: m, a

      IF (PRESENT(comp) .NEQV. PRESENT(fmt_comp)) CALL juDFT_error( &
         'wgauge_bands_row needs a format for the components it is given', &
         calledby='wgauge_bands_row')

      WRITE(iu, '(f12.6)', advance='no') kdist
      DO m = 1, SIZE(evals)
         WRITE(iu, '(2x,f14.8)', advance='no') evals(m)
         IF (PRESENT(comp)) THEN
            DO a = 1, SIZE(comp, 1)
               WRITE(iu, '('//TRIM(fmt_comp)//')', advance='no') comp(a, m)
            END DO
         END IF
      END DO
      WRITE(iu, '(a)') ''
   END SUBROUTINE wgauge_bands_row

END MODULE m_wgauge_bands_io
