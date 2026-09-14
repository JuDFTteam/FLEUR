!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmkb_vac
  USE m_types
  USE m_wann_mmkb_vac
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_mmkb_vac(vacchi, noco, atoms, cell, vacuum, stars, lapw, lapw_b, zMat, zMat_b, qss, vz, vz_b, &
                                 evac, evac_b, jspin, jspin_b, nslibd, nslibd_b, gb, mmn)
    COMPLEX, INTENT(IN) :: vacchi
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_vacuum), INTENT(IN) :: vacuum
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_lapw), INTENT(IN) :: lapw_b
    TYPE(t_mat), INTENT(IN) :: zMat
    TYPE(t_mat), INTENT(IN) :: zMat_b
    REAL, INTENT(IN) :: qss(3)
    REAL, INTENT(IN) :: vz(vacuum%nmzd, 2)
    REAL, INTENT(IN) :: vz_b(vacuum%nmzd, 2)
    REAL, INTENT(IN) :: evac(2)
    REAL, INTENT(IN) :: evac_b(2)
    INTEGER, INTENT(IN) :: jspin
    INTEGER, INTENT(IN) :: jspin_b
    INTEGER, INTENT(IN) :: nslibd
    INTEGER, INTENT(IN) :: nslibd_b
    INTEGER, INTENT(IN) :: gb(3)
    COMPLEX, INTENT(INOUT) :: mmn(:, :)

    CALL wann_mmkb_vac(vacchi, noco%l_noco, atoms%nlotot, qss, SIZE(mmn, 1), cell%z1, vacuum%nmzd, lapw%dim_nv2d(), &
                       stars%mx1, stars%mx2, stars%mx3, stars%ng3, vacuum%nvac, stars%ig, vacuum%nmz, vacuum%delz, stars%ig2, &
                       cell%area, cell%bmat, cell%bbmat, evac, evac_b, lapw%bkpt, lapw_b%bkpt, vz, vz_b, nslibd, nslibd_b, &
                       jspin, jspin_b, lapw%k1, lapw%k2, lapw%k3, lapw_b%k1, lapw_b%k2, lapw_b%k3, SIZE(lapw%k1, 2), lapw%dim_nvd(), &
                       zMat%matsize1, zMat%matsize2, zMat, zMat_b, lapw%nv, lapw_b%nv, cell%omtil, gb, mmn)
  END SUBROUTINE wannierlib_mmkb_vac

END MODULE m_wannierlib_mmkb_vac
