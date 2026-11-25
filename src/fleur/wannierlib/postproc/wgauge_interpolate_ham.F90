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
  USE m_constants, ONLY : oUnit
  USE m_wgauge_hamk, ONLY : t_wgauge_hgauge
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wgauge_interpolate_ham
CONTAINS

  SUBROUTINE wgauge_interpolate_ham(hg, irank)
    TYPE(t_wgauge_hgauge), INTENT(IN) :: hg   !> the domain's bands, already diagonalized
    INTEGER, INTENT(IN) :: irank              ! MPI rank (write only on master)

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)
    !> Nothing to write. The band energies are a property of the output domain and not of
    !> any one quantity measured on it, so they go in with the domain's group, in
    !> electronvolts, whether or not this operator was asked for. Asking for it is still
    !> how a run says it wants the bands and nothing else, which is why the name stays and
    !> why this says so instead of passing in silence.
    WRITE(oUnit,'(a,i0,a)') 'wannierlib interpolation: the band energies of this domain (', &
      SIZE(hg%kdist), ' k-points) are written with the domain itself'
  END SUBROUTINE wgauge_interpolate_ham

END MODULE m_wgauge_interpolate_ham
