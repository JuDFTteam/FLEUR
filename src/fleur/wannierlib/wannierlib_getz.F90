!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_get_z
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
  USE m_types_wannierlib
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_get_z(this, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, nk, jspin, l_real, lapw, zMat)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
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
    INTEGER, ALLOCATABLE :: ev_list(:)

    IF (this%min_band < 1 .OR. this%max_band < this%min_band) THEN
      CALL juDFT_error("Invalid band window in wannierlib_get_z", calledby="wannierlib_get_z")
    END IF

    num_selected = this%max_band - this%min_band + 1
    ALLOCATE(ev_list(num_selected))
    ev_list = [(i, i=this%min_band, this%max_band)]

    CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, nk, cell)

    nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
    CALL zMat%init(l_real, nbasfcn, num_selected)
    CALL read_eig(eig_id, nk, jspin, list=ev_list, zmat=zMat)

    DEALLOCATE(ev_list)
  END SUBROUTINE wannierlib_get_z

END MODULE m_wannierlib_get_z