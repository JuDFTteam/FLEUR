!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_w90_adapter
  USE m_constants, ONLY : oUnit, namat_const
  USE m_juDFT
  USE m_xmlOutput
  USE m_types_atoms
  USE m_types_cell
  USE m_types_kpts
  USE m_types_mpi
  USE m_types_wannierlib
#ifdef CPP_WANNLIB_API
  USE w90_library, ONLY : lib_common_type, w90_set_comm, w90_set_option, w90_input_setopt, &
                          w90_get_nn, w90_get_nnkp, w90_get_gkpb, w90_set_eigval, &
                          w90_set_u_opt, w90_set_m_local, w90_set_u_matrix, &
                          w90_disentangle, w90_project_overlap, w90_wannierise, &
                          w90_get_centres, w90_get_spreads
#endif
  IMPLICIT NONE
#ifdef CPP_WANNLIB_API
  TYPE(lib_common_type), SAVE, PUBLIC :: wannierlib_w90main
#endif
CONTAINS

  SUBROUTINE init_w90(this, atoms, cell, kpts, fmpi, spinors, nntot_w90, nnkp, gkpb, distk)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_mpi), INTENT(IN) :: fmpi
    LOGICAL, INTENT(IN) :: spinors
    INTEGER, INTENT(OUT) :: nntot_w90
    INTEGER, ALLOCATABLE, INTENT(OUT) :: nnkp(:, :)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: gkpb(:, :, :)
    INTEGER, INTENT(IN) :: distk(:)

    INTEGER :: nn
    INTEGER :: na, itype
    INTEGER :: ierr
    INTEGER :: mp_grid(3)
    LOGICAL :: gamma_only
    CHARACTER(LEN=2) :: atom_symbol
    CHARACTER(LEN=32) :: seedname
    CHARACTER(LEN=128), ALLOCATABLE :: atoms_frac(:)

    IF (.NOT.this%l_wannierize) RETURN

#ifndef CPP_WANNLIB_API
    CALL juDFT_error('wannierlib requires Wannier90 module API (w90_library). Rebuild with CPP_WANNLIB_API and matching w90_library.mod.', &
                     calledby='init_w90')
#else
    

      CALL w90_set_comm(wannierlib_w90main, fmpi%mpi_comm)

    mp_grid = kpts%nkpt3  !TODO: is this correct????
    gamma_only = .FALSE.  !TODO: determine from kpts if this is a gamma-only calculation
 

    CALL w90_set_option(wannierlib_w90main, 'num_bands', this%num_bands)
    CALL w90_set_option(wannierlib_w90main, 'num_kpts', kpts%nkptf)
    CALL w90_set_option(wannierlib_w90main, 'num_wann', this%num_wann)
    CALL w90_set_option(wannierlib_w90main, 'num_atoms', atoms%nat)
    CALL w90_set_option(wannierlib_w90main, 'mp_grid', mp_grid)
    CALL w90_set_option(wannierlib_w90main, 'kpoints', kpts%bkf)
    CALL w90_set_option(wannierlib_w90main, 'unit_cell_cart', cell%amat)
    CALL w90_set_option(wannierlib_w90main, 'distk', distk)

    ALLOCATE(atoms_frac(atoms%nat), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating atoms_frac buffer', calledby='init_w90')
    DO na = 1, atoms%nat
      itype = atoms%itype(na)
      atom_symbol = ADJUSTL(namat_const(atoms%nz(itype)))
      WRITE(atoms_frac(na), '(a,1x,3(f24.16,1x))') TRIM(atom_symbol), atoms%taual(1, na), atoms%taual(2, na), atoms%taual(3, na)
    END DO
    !CALL w90_set_option(wannierlib_w90main, 'atoms_frac', atoms_frac)

    CALL w90_set_option(wannierlib_w90main, 'gamma_only', gamma_only)
    CALL w90_set_option(wannierlib_w90main, 'spinors', spinors)
    CALL w90_set_option(wannierlib_w90main, 'dump_inputs', .FALSE.)
    IF (this%dis_num_iter > 0) CALL w90_set_option(wannierlib_w90main, 'num_iter', this%dis_num_iter)
    CALL w90_set_option(wannierlib_w90main, 'dis_win_min', this%dis_win_min)
    CALL w90_set_option(wannierlib_w90main, 'dis_win_max', this%dis_win_max)
    CALL w90_set_option(wannierlib_w90main, 'dis_froz_min', this%dis_froz_min)
    CALL w90_set_option(wannierlib_w90main, 'dis_froz_max', this%dis_froz_max)
    IF (this%dis_mix_ratio > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_mix_ratio', this%dis_mix_ratio)
    IF (this%dis_conv_tol > 0.0) CALL w90_set_option(wannierlib_w90main, 'conv_tol', this%dis_conv_tol)

    seedname = 'fleur_wlib_internal'
    CALL w90_input_setopt(wannierlib_w90main, seedname, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_input_setopt failed in wannierlib adapter', calledby='init_w90')

    CALL w90_get_nn(wannierlib_w90main, nn, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_get_nn failed in wannierlib adapter', calledby='init_w90')
    nntot_w90 = nn
    IF (nntot_w90 <= 0) CALL juDFT_error('wannierlib invalid nntot from w90_get_nn', calledby='init_w90')

    ALLOCATE(nnkp(kpts%nkptf, nntot_w90), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating nnkp', calledby='init_w90')
    CALL w90_get_nnkp(wannierlib_w90main, nnkp, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_get_nnkp failed in wannierlib adapter', calledby='init_w90')

    ALLOCATE(gkpb(3, kpts%nkptf, nntot_w90), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating gkpb', calledby='init_w90')
    CALL w90_get_gkpb(wannierlib_w90main, gkpb, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_get_gkpb failed in wannierlib adapter', calledby='init_w90')
#endif
  END SUBROUTINE init_w90

  SUBROUTINE run_w90(this, mmn, amn, eig)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    COMPLEX, TARGET, INTENT(IN) :: mmn(:, :, :, :)
    COMPLEX, TARGET, INTENT(IN) :: amn(:, :, :)
    REAL, TARGET, INTENT(IN) :: eig(:, :)

    INTEGER :: ierr,num_kpts
    COMPLEX, ALLOCATABLE :: u_matrix(:, :, :), mmn_local(:, :, :, :), amn_local(:, :, :)


    IF (.NOT.this%l_wannierize) RETURN
    num_kpts = SIZE(amn, 3)

#ifndef CPP_WANNLIB_API
    CALL juDFT_error('wannierlib requires Wannier90 module API (w90_library). Rebuild with CPP_WANNLIB_API and matching w90_library.mod.', &
                     calledby='run_w90')
#else


    ALLOCATE(u_matrix(this%num_wann, this%num_wann, num_kpts), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating U matrix', calledby='run_w90')
    u_matrix = CMPLX(0.0, 0.0)
  ALLOCATE(mmn_local(SIZE(mmn,1), SIZE(mmn,2), SIZE(mmn,3), SIZE(mmn,4)), stat=ierr)
  IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating mmn buffer', calledby='run_w90')
  ALLOCATE(amn_local(SIZE(amn,1), SIZE(amn,2), SIZE(amn,3)), stat=ierr)
  IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating amn buffer', calledby='run_w90')
  mmn_local = mmn
  amn_local = amn
    !TODO u_matrix should be stored (in results)  
    CALL w90_set_eigval(wannierlib_w90main, eig)
  CALL w90_set_m_local(wannierlib_w90main, mmn_local)
  CALL w90_set_u_opt(wannierlib_w90main, amn_local)
    CALL w90_set_u_matrix(wannierlib_w90main, u_matrix)

    IF (this%num_bands > this%num_wann) THEN
      CALL w90_disentangle(wannierlib_w90main, oUnit, oUnit, ierr)
      IF (ierr /= 0) CALL juDFT_error('w90_disentangle failed in wannierlib adapter', calledby='run_w90')
    END IF

    CALL w90_project_overlap(wannierlib_w90main, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_project_overlap failed in wannierlib adapter', calledby='run_w90')

    CALL w90_wannierise(wannierlib_w90main, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_wannierise failed in wannierlib adapter', calledby='run_w90')
#endif
  END SUBROUTINE run_w90

  SUBROUTINE report_w90(this)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this

    INTEGER :: ierr, iw
    REAL, ALLOCATABLE :: wann_centres(:, :), wann_spreads(:)
    CHARACTER(LEN=30) :: attributes(5)
    CHARACTER(LEN=6) :: attr_names(5)

    IF (.NOT.this%l_wannierize) RETURN

#ifndef CPP_WANNLIB_API
    CALL juDFT_error('wannierlib requires Wannier90 module API (w90_library). Rebuild with CPP_WANNLIB_API and matching w90_library.mod.', &
                     calledby='report_w90')
#else
    IF (this%num_wann <= 0) THEN
      CALL juDFT_error('wannierlib invalid num_wann for report', calledby='report_w90')
    END IF

    ALLOCATE(wann_centres(3, this%num_wann), wann_spreads(this%num_wann), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating centres/spreads buffers', calledby='report_w90')
    CALL w90_get_centres(wannierlib_w90main, wann_centres)
    CALL w90_get_spreads(wannierlib_w90main, wann_spreads)

    WRITE(oUnit,*)
    WRITE(oUnit,'(a)') 'Wannier90 centres and spreads'
    WRITE(oUnit,'(a)') '  iwf        x                  y                  z               spread'
    DO iw = 1, this%num_wann
      WRITE(oUnit,'(i5,4(2x,f18.10))') iw, wann_centres(1, iw), wann_centres(2, iw), wann_centres(3, iw), wann_spreads(iw)
    END DO

    CALL openXMLElementNoAttributes('wannierlibReport')
    attr_names = (/ 'index ', 'x     ', 'y     ', 'z     ', 'spread' /)
    DO iw = 1, this%num_wann
      attributes = ''
      WRITE(attributes(1),'(i0)') iw
      WRITE(attributes(2),'(f20.12)') wann_centres(1, iw)
      WRITE(attributes(3),'(f20.12)') wann_centres(2, iw)
      WRITE(attributes(4),'(f20.12)') wann_centres(3, iw)
      WRITE(attributes(5),'(f20.12)') wann_spreads(iw)
      CALL writeXMLElement('wannierFunction', attr_names, attributes)
    END DO
    CALL closeXMLElement('wannierlibReport')
#endif
  END SUBROUTINE report_w90

END MODULE m_wannierlib_w90_adapter
