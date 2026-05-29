!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmnkb
  USE m_juDFT
  USE m_wannierlib_mmkb_int
  USE m_wannierlib_get_z
  USE m_wannierlib_mmkb_sph
  USE m_wannierlib_mmk0_sph
  USE m_types
  USE m_types_abc
  USE m_types_atoms
  USE m_types_cell
  USE m_types_input
  USE m_types_kpts
  USE m_types_noco
  USE m_types_nococonv
  USE m_types_radfun
  USE m_types_sym
  USE m_types_usdus
  USE m_types_wannierlib
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_mmnkb(this, num_bands, nntot, nk, kpts, nnkp, gkpb, kdiff, ujug, atoms, cell, input, sym, noco, nococonv, usdus, &
                              radfun, abc, jspin, eig_id, stars, lapw, zMat, mmn)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
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
    INTEGER, INTENT(IN) :: jspin
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_mat), INTENT(IN) :: zMat
    COMPLEX, ALLOCATABLE, INTENT(INOUT) :: mmn(:, :, :, :)

    TYPE(t_mat) :: zMat_b
    TYPE(t_abc) :: abc_b(atoms%ntype)
    TYPE(t_lapw) :: lapw_b
    INTEGER :: kk, nk_b, itype

    IF (.NOT.ALLOCATED(mmn)) THEN
      IF ((num_bands > 0) .AND. (kpts%nkpt > 0) .AND. (nntot > 0)) THEN
        ALLOCATE(mmn(num_bands, num_bands, nntot, kpts%nkpt))
        mmn = CMPLX(0.0, 0.0)
      END IF
    END IF
  
    DO kk = 1, nntot
      nk_b = nnkp(kk, nk)
      CALL wannierlib_get_z(this, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, nk_b, jspin, input%l_real, lapw_b, zMat_b)
      DO itype = 1, atoms%ntype
         CALL abc_b(itype)%init(input, atoms, radfun(itype)%n_r, num_bands, itype)
         CALL abc_b(itype)%calc_abc(input, atoms, sym, cell, lapw_b, num_bands, usdus, noco, nococonv, jspin, itype, zMat_b)
      END DO

      CALL wannierlib_mmkb_sph(atoms, abc_b, abc, kpts%bk(:, nnkp(kk, nk)), gkpb(:, kk, nk), kpts%bk(:, nk), ujug, kdiff, nntot, mmn(:, :, kk, nk))
      CALL wannierlib_mmnkb_int(stars, lapw, lapw_b, jspin, jspin, zMat, zMat_b, gkpb(:, kk, nk), mmn(:, :, kk, nk))
    END DO

    
    
  END SUBROUTINE wannierlib_mmnkb

  

END MODULE m_wannierlib_mmnkb
