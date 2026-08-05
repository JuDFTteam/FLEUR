!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmnkb
  USE m_juDFT
  USE m_wannierlib_mmkb_int
  USE m_melem_get_z
  USE m_wannierlib_mmkb_sph
  USE m_types
  USE m_types_abc
  USE m_types_spinor_layout, ONLY: t_spinor_layout
  USE m_types_atoms
  USE m_types_cell
  USE m_types_input
  USE m_types_kpts
  USE m_types_noco
  USE m_types_nococonv
  USE m_types_radfun
  USE m_types_sym
  USE m_types_usdus
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_mmnkb(min_band, max_band, num_bands, nntot, nk, kpts, nnkp, gkpb, kdiff, ujug, atoms, cell, input, sym, noco, nococonv, usdus, &
                              radfun, abc, jspin, jspin_rad, eig_id, stars, lapw, zMat, mmn, nk_local)
    INTEGER, INTENT(IN) :: min_band, max_band   !> forwarded to get_z for the neighbour k
    INTEGER, INTENT(IN) :: num_bands
    INTEGER, INTENT(IN) :: nntot
    INTEGER, INTENT(IN) :: nk
    TYPE(t_kpts), INTENT(IN) :: kpts
    INTEGER, INTENT(IN) :: nnkp(:, :)
    INTEGER, INTENT(IN) :: gkpb(:, :, :)
    REAL, INTENT(IN) :: kdiff(:, :)
    COMPLEX, INTENT(IN) :: ujug(:, :, :, :, :, :)
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    TYPE(t_usdus), INTENT(IN) :: usdus
    TYPE(t_radfun), INTENT(IN) :: radfun(atoms%ntype)
    TYPE(t_abc), INTENT(IN) :: abc(:)
    INTEGER, INTENT(IN) :: jspin       ! spin fisico (record del eig)
    INTEGER, INTENT(IN) :: jspin_rad   ! indice radial (=1 si jspins=1)
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_mat), INTENT(IN) :: zMat
    COMPLEX, ALLOCATABLE, INTENT(INOUT) :: mmn(:, :, :, :)
    INTEGER, INTENT(IN) :: nk_local

    TYPE(t_mat) :: zMat_b
    TYPE(t_abc) :: abc_b(atoms%ntype)
    TYPE(t_lapw) :: lapw_b
    INTEGER :: kk, nk_b, itype
    LOGICAL :: l_real_wann
    TYPE(t_spinor_layout) :: layout, layout_b

    IF (.NOT.ALLOCATED(mmn)) THEN
      IF ((num_bands > 0) .AND. (kpts%nkpt > 0) .AND. (nntot > 0)) THEN
        ALLOCATE(mmn(num_bands, num_bands, nntot, kpts%nkpt))
        mmn = CMPLX(0.0, 0.0)
      END IF
    END IF
  
    l_real_wann = input%l_real .AND. .NOT. noco%l_soc
    CALL layout%init(input, noco, lapw, atoms)
    DO kk = 1, nntot
      nk_b = nnkp(nk, kk)
      CALL melem_get_z(min_band, max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, nk_b, jspin, l_real_wann, lapw_b, zMat_b)
      CALL layout_b%init(input, noco, lapw_b, atoms)
      DO itype = 1, atoms%ntype
         CALL abc_b(itype)%init(input, atoms, num_bands, itype)
         CALL abc_b(itype)%calc_abc(input, atoms, sym, cell, lapw_b, num_bands, usdus, noco, nococonv, jspin_rad, itype, zMat_b)
      END DO

      CALL wannierlib_mmnkb_int(stars, lapw, lapw_b, jspin_rad, jspin_rad, zMat, zMat_b, gkpb(:, nk, kk), mmn(:, :, kk, nk_local), &
                                ioff=layout%row_offset(jspin), ioff_b=layout_b%row_offset(jspin))
      CALL wannierlib_mmkb_sph(atoms, abc, abc_b, kpts%bkf(:, nnkp(nk, kk)), gkpb(:, nk, kk), kpts%bkf(:, nk), ujug, kdiff, nntot, mmn(:, :, kk, nk_local))
    END DO

    
    
  END SUBROUTINE wannierlib_mmnkb

  

END MODULE m_wannierlib_mmnkb
