!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The interpolated band structure, written out.
!>
!>  It takes a domain whose bands are already built -- H_W(k), the transform to real space
!>  and the diagonalization all belong to m_wgauge_hamk, which every interpolation driver
!>  shares -- so what is decided here is only the two units the bands are written in.
MODULE m_wgauge_interpolate_ham
  USE m_juDFT
  USE m_wgauge_bands_io, ONLY: wgauge_bands_open, wgauge_bands_row
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_wgauge_hamk, ONLY : t_wgauge_hgauge
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wgauge_interpolate_ham
CONTAINS

  SUBROUTINE wgauge_interpolate_ham(hg, out1, out2, irank)
    TYPE(t_wgauge_hgauge), INTENT(IN) :: hg   !> the domain's bands, already diagonalized
    !> The names the files take, decided by the caller from the exposure table plus the
    !> domain suffix. Unallocated off rank 0, which never reaches them.
    CHARACTER(LEN=*), INTENT(IN) :: out1
    CHARACTER(LEN=*), INTENT(IN) :: out2
    INTEGER, INTENT(IN) :: irank              ! MPI rank (write only on master)

    INTEGER :: ip, np, iu, iuev

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)
    CALL timestart('wgauge_interpolate_ham')
    np = SIZE(hg%kdist)

    CALL wgauge_bands_open(iu,   out1, '# kdist   E_1..E_numwann  (Htr)')
    CALL wgauge_bands_open(iuev, out2, '# kdist   E_1..E_numwann  (eV, absolute)')
    DO ip = 1, np
      CALL wgauge_bands_row(iu,   hg%kdist(ip), hg%evals(:, ip))
      CALL wgauge_bands_row(iuev, hg%kdist(ip), hartree_to_ev_const * hg%evals(:, ip))
    END DO
    CLOSE(iu); CLOSE(iuev)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib interpolation: wrote '//TRIM(out1)//'.dat (', np, ' k-points)'
    CALL timestop('wgauge_interpolate_ham')
  END SUBROUTINE wgauge_interpolate_ham

END MODULE m_wgauge_interpolate_ham
