!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_melem_get_z
  !> The eigenvectors of one k-point, read and laid out.
  !>
  !> It answers the same question for every consumer -- give me the states of this k in the
  !> shape the contraction expects -- and it needs nothing of whatever feature is asking:
  !> the band window arrives as two integers. Unlike a cache keyed on one k-point it can
  !> also be asked for a k outside the irreducible zone, which it serves by rotating the
  !> states of the parent, and for two different k-points at once.
  USE m_juDFT
  USE m_eig66_io, ONLY: read_eig
  USE m_types_atoms
  USE m_types_cell
  USE m_types_input
  USE m_types_kpts
  USE m_types_mat
  USE m_types_noco
  USE m_types_nococonv
  USE m_types_sym
  USE m_types_lapw
  USE m_types_spinor_layout, ONLY: t_spinor_layout
  IMPLICIT NONE
CONTAINS

  !> min_band/max_band select the band window; they are plain integers so that this
  !> routine, and everything that only forwards them, needs nothing of the Wannier input.
  SUBROUTINE melem_get_z(min_band, max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, nk, jspin, l_real, lapw, zMat)
    INTEGER, INTENT(IN) :: min_band, max_band
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: nk
    INTEGER, INTENT(IN) :: jspin
    LOGICAL, INTENT(IN) :: l_real
    TYPE(t_lapw), INTENT(OUT) :: lapw
    TYPE(t_mat), INTENT(OUT) :: zMat

    INTEGER :: num_selected, i, nbasfcn
    TYPE(t_spinor_layout) :: layout
    INTEGER, ALLOCATABLE :: ev_list(:)

    IF (min_band < 1 .OR. max_band < min_band) THEN
      CALL juDFT_error("Invalid band window in melem_get_z", calledby="melem_get_z")
    END IF

    num_selected = max_band - min_band + 1
    ALLOCATE(ev_list(num_selected))
    ev_list = [(i, i=min_band, max_band)]

    CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell)

    ! noco: el espinor 2N esta en UN solo record (jspin=1). Para noco leemos siempre el
    ! record 1 y devolvemos el 2N COMPLETO; abcof/mmkb_int extraen el spin via offsets.
    CALL layout%init(input, noco, lapw, atoms)
    nbasfcn = layout%nbasfcn(lapw, atoms)
    CALL zMat%init(l_real, nbasfcn, num_selected)
    CALL read_eig(eig_id, nk, MERGE(1, jspin, noco%l_noco), list=ev_list, zmat=zMat,&
          kpts=kpts,input=input,noco=noco,nococonv=nococonv,sym=sym,atoms=atoms,cell=cell)

    DEALLOCATE(ev_list)
  END SUBROUTINE melem_get_z

END MODULE m_melem_get_z