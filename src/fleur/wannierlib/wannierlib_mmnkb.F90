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
  USE m_melem_mmkb_int
  USE m_matrix_element_factory, ONLY: matrix_element_states
  USE m_melem_mmkb_sph
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
  USE m_types_enpara
  USE m_types_potden
  USE m_types_mpi
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_mmnkb(min_band, max_band, num_bands, nntot, nk, kpts, nnkp, gkpb, kdiff, ujug, atoms, cell, input, sym, noco, nococonv, usdus, &
                              radfun, abc, jspin, jspin_rad, eig_id, stars, lapw, zMat, mmn, nk_local, &
                              enpara, vtot, fmpi)
    INTEGER, INTENT(IN) :: min_band, max_band   !> the band window asked for at the neighbour k
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
    TYPE(t_enpara), INTENT(IN) :: enpara   !> the factory generates the radial functions itself
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_mpi), INTENT(IN) :: fmpi

    TYPE(t_mat), POINTER :: zMat_b(:)   !> points into the factory cache, one entry per record
    TYPE(t_abc), POINTER :: abc_b(:, :) !> (2,ntype), likewise
    TYPE(t_lapw) :: lapw_b
    INTEGER :: kk, nk_b, irec
    INTEGER, ALLOCATABLE :: ev_list(:)
    TYPE(t_spinor_layout) :: layout, layout_b

    IF (.NOT.ALLOCATED(mmn)) THEN
      IF ((num_bands > 0) .AND. (kpts%nkpt > 0) .AND. (nntot > 0)) THEN
        ALLOCATE(mmn(num_bands, num_bands, nntot, kpts%nkpt))
        mmn = CMPLX(0.0, 0.0)
      END IF
    END IF
  
    ev_list = [(irec, irec = min_band, max_band)]
    !> Non-collinearly the whole spinor is one record; otherwise each spin channel is its
    !> own, and this pass reaches its block by row offset rather than by stacking them.
    irec = MERGE(1, jspin, noco%l_noco)
    CALL layout%init(input, noco, lapw, atoms)
    DO kk = 1, nntot
      nk_b = nnkp(nk, kk)
      !> The neighbour's basis is needed here as well as by the factory, so it is built
      !> here and handed over. The states themselves stay in the factory cache, which holds
      !> more than one k-point: asking for this neighbour does not discard the k it belongs
      !> to, nor the neighbour before it.
      CALL lapw_b%init(input, noco, nococonv, kpts, atoms, sym, nk_b, cell)
      CALL matrix_element_states(eig_id, nk_b, input, atoms, sym, cell, noco, nococonv, &
                                 enpara, lapw_b, vtot, fmpi, zMat_b, abc_b, ev_list=ev_list, &
                                 l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), kpts=kpts)
      CALL layout_b%init(input, noco, lapw_b, atoms)

      CALL melem_mmkb_int(stars, lapw, lapw_b, jspin_rad, jspin_rad, zMat, zMat_b(irec), gkpb(:, nk, kk), mmn(:, :, kk, nk_local), &
                                ioff=layout%row_offset(jspin), ioff_b=layout_b%row_offset(jspin))
      CALL melem_mmkb_sph(atoms, abc, abc_b(jspin, :), kpts%bkf(:, nnkp(nk, kk)), gkpb(:, nk, kk), kpts%bkf(:, nk), ujug, kdiff, nntot, mmn(:, :, kk, nk_local))
    END DO

    
    
  END SUBROUTINE wannierlib_mmnkb

  

END MODULE m_wannierlib_mmnkb
