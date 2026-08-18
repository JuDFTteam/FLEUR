!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_w90_adapter
  USE m_constants, ONLY : oUnit, namat_const, hartree_to_ev_const, bohr_to_angstrom_const
  USE m_juDFT
  USE m_xmlOutput
  USE m_types_atoms
  USE m_types_cell
  USE m_types_kpts
  USE m_types_mpi
  USE m_types_wannierlib
  USE m_types_melem_bmesh
#ifdef CPP_WANNLIB_API
  USE w90_library, ONLY : lib_common_type, w90_set_comm, w90_set_option, w90_input_setopt, &
                          w90_get_nn, w90_get_nnkp, w90_get_gkpb, w90_set_eigval, &
                          w90_set_u_opt, w90_set_m_local, w90_set_u_matrix, &
                          w90_disentangle, w90_project_overlap, w90_wannierise, &
                          w90_get_centres, w90_get_spreads
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_w90, run_w90, wannierlib_get_bmesh, report_w90
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
    
    call timestart('init_w90')
    CALL w90_set_comm(wannierlib_w90main, fmpi%mpi_comm)

    mp_grid = kpts%nkpt3  !TODO: is this correct????
    gamma_only = .FALSE.  !TODO: determine from kpts if this is a gamma-only calculation
 

    CALL w90_set_option(wannierlib_w90main, 'num_bands', this%num_bands)
    CALL w90_set_option(wannierlib_w90main, 'num_kpts', kpts%nkptf)
    CALL w90_set_option(wannierlib_w90main, 'num_wann', this%num_wann)
    CALL w90_set_option(wannierlib_w90main, 'num_atoms', atoms%nat)
    CALL w90_set_option(wannierlib_w90main, 'mp_grid', mp_grid)
    CALL w90_set_option(wannierlib_w90main, 'kpoints', kpts%bkf)
    CALL w90_set_option(wannierlib_w90main, 'unit_cell_cart', bohr_to_angstrom_const*cell%amat)
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
    IF (this%dis_num_iter > 0) CALL w90_set_option(wannierlib_w90main, 'dis_num_iter', this%dis_num_iter)  ! disentanglement (XML numIter)
    IF (this%num_iter > 0)     CALL w90_set_option(wannierlib_w90main, 'num_iter', this%num_iter)          ! MLWF/wannierise (XML wannIter)
    CALL w90_set_option(wannierlib_w90main, 'dis_win_min', hartree_to_ev_const*this%dis_win_min)
    CALL w90_set_option(wannierlib_w90main, 'dis_win_max', hartree_to_ev_const*this%dis_win_max)
    CALL w90_set_option(wannierlib_w90main, 'dis_froz_min', hartree_to_ev_const*this%dis_froz_min)
    CALL w90_set_option(wannierlib_w90main, 'dis_froz_max', hartree_to_ev_const*this%dis_froz_max)
    IF (this%dis_mix_ratio > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_mix_ratio', this%dis_mix_ratio)
    IF (this%dis_conv_tol > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_conv_tol', this%dis_conv_tol)  ! disentanglement (XML disConvTol)
    IF (this%conv_tol > 0.0)     CALL w90_set_option(wannierlib_w90main, 'conv_tol', this%conv_tol)          ! MLWF/wannierise (XML wannConvTol)
    !> Only set when asked: Wannier90's own default is .FALSE., so staying silent keeps
    !> every existing run byte-identical.
    IF (this%precond) CALL w90_set_option(wannierlib_w90main, 'precond', .TRUE.)                            ! MLWF/wannierise (XML precond)

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
    call timestop('init_w90')
#endif
  END SUBROUTINE init_w90

  !> Run the Wannier90 library: hand over the eigenvalues, the overlaps M and the projections A,
  !> then disentangle -> project -> wannierise. Returns the resulting gauge factors, which is all
  !> any downstream consumer needs:
  !>    u_matrix  the MLWF gauge U(k)
  !>    u_opt     the disentanglement matrix (the projections after w90 has processed them)
  !> The full Wannier gauge is V(k) = u_opt(k) . u_matrix(k). Evaluating operators in that gauge
  !> is NOT this module's job: the adapter only translates to and from the library.
  SUBROUTINE run_w90(this, cell, kpts, mmn, amn, eig, irank, u_matrix, u_opt)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, TARGET, INTENT(IN) :: mmn(:, :, :, :)
    COMPLEX, TARGET, INTENT(IN) :: amn(:, :, :)
    REAL, TARGET, INTENT(IN) :: eig(:, :)
    INTEGER, INTENT(IN) :: irank
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: u_matrix(:, :, :)   ! (nw,nw,nk) MLWF gauge
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: u_opt(:, :, :)      ! (nb,nw,nk) disentangled projections

    INTEGER :: ierr, num_kpts
    COMPLEX, ALLOCATABLE :: mmn_local(:, :, :, :)
    REAL, ALLOCATABLE, TARGET :: eig_ev(:, :)         ! eigenvalues in eV for the w90 library (banner honesty)

    ! allocate the INTENT(OUT) results before any early return, so a caller can never pass
    ! unallocated arrays on to the matrix-element layer
    ALLOCATE(u_matrix(0, 0, 0), u_opt(0, 0, 0))
    IF (.NOT.this%l_wannierize) RETURN
    DEALLOCATE(u_matrix, u_opt)
    num_kpts = SIZE(amn, 3)

#ifndef CPP_WANNLIB_API
    CALL juDFT_error('wannierlib requires Wannier90 module API (w90_library). Rebuild with CPP_WANNLIB_API and matching w90_library.mod.', &
                     calledby='run_w90')
#else
    call timestart('run_w90')

    ALLOCATE(u_matrix(this%num_wann, this%num_wann, num_kpts), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating U matrix', calledby='run_w90')
    u_matrix = CMPLX(0.0, 0.0)
    ALLOCATE(mmn_local(SIZE(mmn,1), SIZE(mmn,2), SIZE(mmn,3), SIZE(mmn,4)), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating mmn buffer', calledby='run_w90')
    ALLOCATE(u_opt(SIZE(amn,1), SIZE(amn,2), SIZE(amn,3)), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating amn buffer', calledby='run_w90')
    mmn_local = mmn
    u_opt = amn
    !TODO u_matrix should be stored (in results)
    ALLOCATE(eig_ev(SIZE(eig,1), SIZE(eig,2)), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating eig_ev buffer', calledby='run_w90')
    eig_ev = hartree_to_ev_const * eig   ! w90 library expects eV; FLEUR eig are Hartree
    CALL w90_set_eigval(wannierlib_w90main, eig_ev)
    CALL w90_set_m_local(wannierlib_w90main, mmn_local)
    CALL w90_set_u_opt(wannierlib_w90main, u_opt)
    CALL w90_set_u_matrix(wannierlib_w90main, u_matrix)

    IF (this%num_bands > this%num_wann) THEN
      CALL w90_disentangle(wannierlib_w90main, oUnit, oUnit, ierr)
      IF (ierr /= 0) CALL juDFT_error('w90_disentangle failed in wannierlib adapter', calledby='run_w90')
    END IF

    CALL w90_project_overlap(wannierlib_w90main, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_project_overlap failed in wannierlib adapter', calledby='run_w90')

    CALL w90_wannierise(wannierlib_w90main, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_wannierise failed in wannierlib adapter', calledby='run_w90')

    call timestop('run_w90')
#endif
  END SUBROUTINE run_w90

  !> Export the coarse-mesh b-shell / neighbour information that the Wannier90 kmesh setup
  !> produced, as a plain t_melem_bmesh. This is the ONE piece of Wannier90 state the
  !> matrix-element layer needs (for the position / Berry-connection operator), and handing it
  !> over like this keeps wannierlib_w90main private to this module and keeps m_melem_* free of
  !> the Wannier90 library. Call after run_w90 (the centres are only meaningful once wannierised).
  SUBROUTINE wannierlib_get_bmesh(this, kpts, bmesh)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    !> INOUT and not OUT: the topology is already in here, put there before the
    !> wannierisation ran, and OUT would default-initialise it away. What this adds is the
    !> part only Wannier90 knows -- the shell weights and the b vectors it chose.
    TYPE(t_melem_bmesh), INTENT(INOUT) :: bmesh
#ifdef CPP_WANNLIB_API
    IF (ALLOCATED(bmesh%wb)) DEALLOCATE(bmesh%wb)
    IF (ALLOCATED(bmesh%bk)) DEALLOCATE(bmesh%bk)
    IF (ALLOCATED(bmesh%centres)) DEALLOCATE(bmesh%centres)
    IF (bmesh%nntot > 0 .AND. wannierlib_w90main%kmesh_info%nntot /= bmesh%nntot) &
      CALL juDFT_error("wannierlib_get_bmesh: Wannier90 counts a different number of "// &
                       "neighbours than the topology already set", calledby="wannierlib_get_bmesh")
    bmesh%nntot = wannierlib_w90main%kmesh_info%nntot
    IF (bmesh%nntot < 1) RETURN
    ! w90 shapes: nnlist(num_kpts, nntot), wb(nntot), bk(3, nntot, num_kpts)
    ALLOCATE(bmesh%wb(bmesh%nntot))
    ALLOCATE(bmesh%bk(3, bmesh%nntot, kpts%nkptf))
    !> The neighbour list may already be here, put there before the wannierisation ran and
    !> read since by whoever built the overlaps. Overwriting it would hand a later spin
    !> channel a different list from the one its overlaps were computed with, so it is only
    !> filled when nobody set it -- and checked against this source when somebody did.
    IF (.NOT.ALLOCATED(bmesh%nnlist)) THEN
      ALLOCATE(bmesh%nnlist(kpts%nkptf, bmesh%nntot))
      bmesh%nnlist = wannierlib_w90main%kmesh_info%nnlist(1:kpts%nkptf, 1:bmesh%nntot)
    ELSE IF (ANY(bmesh%nnlist /= wannierlib_w90main%kmesh_info%nnlist(1:kpts%nkptf, 1:bmesh%nntot))) THEN
      CALL juDFT_error("wannierlib_get_bmesh: the neighbour list already set does not match "// &
                       "the one Wannier90 reports", calledby="wannierlib_get_bmesh")
    END IF
    bmesh%wb     = wannierlib_w90main%kmesh_info%wb(1:bmesh%nntot)
    bmesh%bk     = wannierlib_w90main%kmesh_info%bk(:, 1:bmesh%nntot, 1:kpts%nkptf)
    IF (this%num_wann > 0) THEN
      ALLOCATE(bmesh%centres(3, this%num_wann))
      CALL w90_get_centres(wannierlib_w90main, bmesh%centres)
    END IF
#else
    bmesh%nntot = 0   ! no Wannier90 library -> no b-mesh; the position operator is skipped
#endif
  END SUBROUTINE wannierlib_get_bmesh

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
