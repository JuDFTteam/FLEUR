!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmk0_vac
  USE m_types
  USE m_wann_mmk0_vac
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_mmk0_vac(noco, atoms, cell, vacuum, stars, enpara, lapw, zMat, qss, vz, jspin, nslibd, mmn)
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_vacuum), INTENT(IN) :: vacuum
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_mat), INTENT(IN) :: zMat
    REAL, INTENT(IN) :: qss(3)
    REAL, INTENT(IN) :: vz(vacuum%nmzd, 2)
    INTEGER, INTENT(IN) :: jspin
    INTEGER, INTENT(IN) :: nslibd
    COMPLEX, INTENT(INOUT) :: mmn(:, :)

    CALL wann_mmk0_vac(noco%l_noco, atoms%nlotot, qss, cell%z1, vacuum%nmzd, lapw%dim_nv2d(), stars%mx1, stars%mx2, stars%mx3, &
                       stars%ng3, vacuum%nvac, stars%ig, vacuum%nmz, vacuum%delz, stars%ig2, cell%area, cell%bmat, cell%bbmat, &
                       enpara%evac0(:, jspin), lapw%bkpt, vz, nslibd, jspin, lapw%k1, lapw%k2, lapw%k3, SIZE(lapw%k1, 2), &
                       lapw%dim_nvd(), zMat%matsize1, zMat%matsize2, zMat, lapw%nv, cell%omtil, mmn)
  END SUBROUTINE wannierlib_mmk0_vac

END MODULE m_wannierlib_mmk0_vac
