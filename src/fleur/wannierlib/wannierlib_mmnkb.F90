!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmnkb
  USE m_juDFT
  USE m_melem_overlap, ONLY: melem_overlap_states
  USE m_types_melem_vacabc, ONLY: t_melem_vacabc
  USE m_matrix_element_factory, ONLY: matrix_element_states
  USE m_types
  USE m_types_abc
  USE m_types_spinor_layout, ONLY: t_spinor_layout
  USE m_types_atoms
  USE m_types_cell
  USE m_types_input
  USE m_types_kpts
  USE m_types_noco
  USE m_types_nococonv
  USE m_types_sym
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_types_melem_bmesh, ONLY: t_melem_bmesh
  USE m_types_enpara
  USE m_types_potden
  USE m_types_mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_mmnkb
CONTAINS

  SUBROUTINE wannierlib_mmnkb(manifold, bmesh, nk, kpts, ujug, atoms, cell, input, sym, noco, nococonv, &
                              abc, jspin, jspin_rad, eig_id, stars, lapw, zMat, mmn, nk_local, &
                              enpara, vtot, fmpi, vacuum)
    TYPE(t_melem_manifold), INTENT(IN) :: manifold   !> the band window, and how wide it is
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh   !> which k is the b-th neighbour, and by which G
    INTEGER, INTENT(IN) :: nk
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: ujug(:, :, :, :, :, :)
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    TYPE(t_abc), INTENT(IN) :: abc(:)
    INTEGER, INTENT(IN) :: jspin       ! physical spin (the eig record)
    INTEGER, INTENT(IN) :: jspin_rad   ! radial index (=1 when jspins=1)
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_mat), INTENT(IN) :: zMat
    COMPLEX, ALLOCATABLE, INTENT(INOUT) :: mmn(:, :, :, :)
    INTEGER, INTENT(IN) :: nk_local
    TYPE(t_enpara), INTENT(IN) :: enpara   !> the factory generates the radial functions itself
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_mpi), INTENT(IN) :: fmpi
    !> Only a film uses it: the expansions stay unbuilt otherwise and the overlap skips them.
    TYPE(t_vacuum), INTENT(IN) :: vacuum

    TYPE(t_mat), POINTER :: zMat_b(:)   !> points into the factory cache, one entry per record
    TYPE(t_abc), POINTER :: abc_b(:, :) !> (2,ntype), likewise
    TYPE(t_lapw) :: lapw_b
    INTEGER :: kk, nk_b, irec
    INTEGER, ALLOCATABLE :: ev_list(:)
    TYPE(t_spinor_layout) :: layout, layout_b
    TYPE(t_melem_vacabc) :: vac, vac_b

    IF (.NOT.ALLOCATED(mmn)) THEN
      IF ((manifold%num_bands > 0) .AND. (kpts%nkpt > 0) .AND. (bmesh%nntot > 0)) THEN
        ALLOCATE(mmn(manifold%num_bands, manifold%num_bands, bmesh%nntot, kpts%nkpt))
        mmn = CMPLX(0.0, 0.0)
      END IF
    END IF
  
    ev_list = [(irec, irec = manifold%min_band, manifold%max_band)]
    !> Non-collinearly the whole spinor is one record; otherwise each spin channel is its
    !> own, and this pass reaches its block by row offset rather than by stacking them.
    irec = MERGE(1, jspin, noco%l_noco)
    CALL layout%init(input, noco, lapw, atoms)
    !> This k is the bra of every neighbour below, so its expansion is built once.
    IF (input%film) CALL vac%calc(vacuum, cell, enpara, vtot, lapw, jspin_rad, zMat, &
                                  manifold%num_bands, ioff=layout%row_offset(jspin))
    DO kk = 1, bmesh%nntot
      nk_b = bmesh%nnlist(nk, kk)
      !> The neighbour's basis is needed here as well as by the factory, so it is built
      !> here and handed over. The states themselves stay in the factory cache, which holds
      !> more than one k-point: asking for this neighbour does not discard the k it belongs
      !> to, nor the neighbour before it.
      CALL lapw_b%init(input, noco, nococonv, kpts, atoms, sym, nk_b, cell)
      CALL matrix_element_states(eig_id, nk_b, input, atoms, sym, cell, noco, nococonv, &
                                 enpara, lapw_b, vtot, fmpi, zMat_b, abc_b, ev_list=ev_list, &
                                 l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), kpts=kpts)
      CALL layout_b%init(input, noco, lapw_b, atoms)
      IF (input%film) CALL vac_b%calc(vacuum, cell, enpara, vtot, lapw_b, jspin_rad, &
                                      zMat_b(irec), manifold%num_bands, &
                                      ioff=layout_b%row_offset(jspin))

      CALL melem_overlap_states(stars, atoms, lapw, lapw_b, zMat, zMat_b(irec), &
                                abc, abc_b(jspin, :), jspin_rad, jspin_rad, &
                                kpts%bkf(:, nk), kpts%bkf(:, bmesh%nnlist(nk, kk)), &
                                bmesh%gkpb(:, nk, kk), ujug, bmesh%kdiff, bmesh%nntot, &
                                ioff_a=layout%row_offset(jspin), &
                                ioff_b=layout_b%row_offset(jspin), &
                                ovl=mmn(:, :, kk, nk_local), vac_a=vac, vac_b=vac_b)
    END DO

    
    
  END SUBROUTINE wannierlib_mmnkb

  

END MODULE m_wannierlib_mmnkb
