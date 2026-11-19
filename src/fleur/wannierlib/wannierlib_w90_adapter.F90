!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_w90_adapter
  USE m_constants, ONLY : oUnit, namat_const, hartree_to_ev_const, bohr_to_angstrom_const
  USE m_juDFT
  USE m_judft_string, ONLY : int2str
  USE m_xmlOutput
  USE m_types_atoms
  USE m_types_cell
  USE m_types_kpts
  USE m_types_mpi
  USE m_types_wannierlib
  USE m_types_wgauge_bmesh
#ifdef CPP_WANNLIB_API
  USE w90_library, ONLY : lib_common_type, w90_set_comm, w90_set_option, w90_input_setopt, &
                          w90_get_nn, w90_get_nnkp, w90_get_gkpb, w90_set_eigval, &
                          w90_set_u_opt, w90_set_m_local, w90_set_u_matrix, &
                          w90_disentangle, w90_project_overlap, w90_wannierise, w90_create_kmesh, &
                          w90_get_centres, w90_get_spreads, w90_transport
  USE w90_library_extra, ONLY : write_chkpt
#endif
  IMPLICIT NONE

  !> How many neighbour shells Wannier90 may look through for a set of b vectors that
  !> satisfies the completeness relation. The classic path (wann_wan90prep.F) writes 200 and
  !> nobody in the group ever tuned it, because in a bulk cell 200 is already generous.
  !>
  !> A FILM is what makes this a real parameter. The out-of-plane reciprocal vector has the
  !> fixed length 2pi/dtilda while the in-plane ones shrink as the mesh is refined, so the one
  !> shell that carries any z information sinks further down a list ordered by length with
  !> every refinement. Run out of list before reaching it and Wannier90 rejects shell after
  !> shell by SVD and stops -- measured on Fe monolayer, dtilda = 9.68 bohr: 200 shells reach
  !> 1.0885 Ang^-1 and 2pi/dtilda is 1.227 Ang^-1, so 48x48x1 dies just short of it.
  !>
  !> 1000 covers the meshes this is used with (through 64x64x1). It costs a search that grows
  !> linearly and one nnshell(num_kpts, search_shells) table -- 16 MB at 4096 k-points, which
  !> is nothing next to what a run of that size already holds. Wannier90 sets no upper bound
  !> of its own; it dimensions its arrays from this value.
  INTEGER, PARAMETER :: WANNIERLIB_SEARCH_SHELLS = 1000

  PRIVATE
  PUBLIC :: init_w90, run_w90, wannierlib_get_bmesh, report_w90, WANNIERLIB_SEARCH_SHELLS
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

    IF (.NOT.this%l_wannierize) RETURN

#ifndef CPP_WANNLIB_API
    CALL juDFT_error('wannierlib requires Wannier90 module API (w90_library). Rebuild with CPP_WANNLIB_API and matching w90_library.mod.', &
                     calledby='init_w90')
#else
    
    call timestart('init_w90')
    CALL set_w90_options(this, atoms, cell, kpts, fmpi, spinors, distk)
    CALL get_w90_topology(kpts, nntot_w90, nnkp, gkpb)
    call timestop('init_w90')
#endif
  END SUBROUTINE init_w90


#ifdef CPP_WANNLIB_API
  !> Hand Wannier90 everything the input states, and commit it. Nothing is decided here:
  !> what is set is what <wannierlib> carries, converted where the two differ on units --
  !> Hartree here and electronvolt there, Bohr here and Angstrom there.
  SUBROUTINE set_w90_options(this, atoms, cell, kpts, fmpi, spinors, distk)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_mpi), INTENT(IN) :: fmpi
    LOGICAL, INTENT(IN) :: spinors
    INTEGER, INTENT(IN) :: distk(:)

    INTEGER :: na, itype, ierr
    INTEGER :: mp_grid(3)
    LOGICAL :: gamma_only
    CHARACTER(LEN=2) :: atom_symbol
    CHARACTER(LEN=32) :: seedname
    CHARACTER(LEN=128), ALLOCATABLE :: atoms_frac(:)

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
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating atoms_frac buffer', calledby='set_w90_options')
    DO na = 1, atoms%nat
      itype = atoms%itype(na)
      atom_symbol = ADJUSTL(namat_const(atoms%nz(itype)))
      WRITE(atoms_frac(na), '(a,1x,3(f24.16,1x))') TRIM(atom_symbol), atoms%taual(1, na), atoms%taual(2, na), atoms%taual(3, na)
    END DO
    !CALL w90_set_option(wannierlib_w90main, 'atoms_frac', atoms_frac)

    CALL w90_set_option(wannierlib_w90main, 'gamma_only', gamma_only)
    CALL w90_set_option(wannierlib_w90main, 'spinors', spinors)
    CALL w90_set_option(wannierlib_w90main, 'dump_inputs', .FALSE.)
    !> >= 0 and not > 0: num_iter = 0 is a value the schema documents -- no minimisation,
    !> the projection gauge kept -- and the sentinel for silence is -1. Guarding on > 0 let
    !> a requested 0 fall through to Wannier90's own default of 100 iterations, so an input
    !> asking for no minimisation quietly got a hundred steps of it.
    IF (this%dis_num_iter >= 0) CALL w90_set_option(wannierlib_w90main, 'dis_num_iter', this%dis_num_iter)  ! disentanglement (XML numIter)
    IF (this%num_iter >= 0)     CALL w90_set_option(wannierlib_w90main, 'num_iter', this%num_iter)          ! MLWF/wannierise (XML wannIter)
    !> How many neighbour shells w90 may try before it gives up on the B1 condition. Its own
    !> default is 36, which is not enough for a film: the mesh is dense in plane and a single
    !> point in z, so shell after shell is rejected by the SVD test and 36 runs out at
    !> 14x14x1. The classic path has written 200 into the .win since it was added
    !> (wann_wan90prep.F:387), and 200 is what the group uses; passing nothing left the
    !> library mode on w90's default and made it fail where the classic would not.
    CALL w90_set_option(wannierlib_w90main, 'search_shells', WANNIERLIB_SEARCH_SHELLS)
    !> How close two |b| have to be to count as one shell, and how exactly the completeness
    !> relation has to hold. Silence leaves Wannier90 on its own 1e-6.
    IF (this%kmesh_tol >= 0.0) CALL w90_set_option(wannierlib_w90main, 'kmesh_tol', this%kmesh_tol)
    CALL w90_set_option(wannierlib_w90main, 'dis_win_min', hartree_to_ev_const*this%dis_win_min)
    CALL w90_set_option(wannierlib_w90main, 'dis_win_max', hartree_to_ev_const*this%dis_win_max)
    CALL w90_set_option(wannierlib_w90main, 'dis_froz_min', hartree_to_ev_const*this%dis_froz_min)
    CALL w90_set_option(wannierlib_w90main, 'dis_froz_max', hartree_to_ev_const*this%dis_froz_max)
    IF (this%dis_mix_ratio > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_mix_ratio', this%dis_mix_ratio)
    IF (this%dis_conv_tol > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_conv_tol', this%dis_conv_tol)  ! disentanglement (XML disConvTol)
    IF (this%conv_tol > 0.0)     CALL w90_set_option(wannierlib_w90main, 'conv_tol', this%conv_tol)          ! MLWF/wannierise (XML wannConvTol)
    !> Without a window the MLWF iteration never tests the tolerance it was just given: the
    !> Wannier90 reader forces conv_window to -1 and the check in wannierise.F90 asks for
    !> above one, so conv_tol alone buys nothing there and the run costs its whole num_iter.
    !> The disentanglement is not in that position -- its window defaults to 3 and a negative
    !> one is an input error. Both are said only when asked, so silence keeps the behaviour
    !> every stored reference was made with.
    IF (this%dis_conv_window > 0) CALL w90_set_option(wannierlib_w90main, 'dis_conv_window', this%dis_conv_window)
    IF (this%conv_window > 0)     CALL w90_set_option(wannierlib_w90main, 'conv_window', this%conv_window)
    !> Projectability disentanglement: set only when asked, so silence leaves Wannier90 on
    !> its own energy-window default. The two thresholds go with it -- W90 validates both
    !> to [0,1] and aborts on anything else.
    IF (this%dis_froz_proj) THEN
       CALL w90_set_option(wannierlib_w90main, 'dis_froz_proj', .TRUE.)
       CALL w90_set_option(wannierlib_w90main, 'dis_proj_min', this%dis_proj_min)
       CALL w90_set_option(wannierlib_w90main, 'dis_proj_max', this%dis_proj_max)
    END IF
    IF (this%precond) CALL w90_set_option(wannierlib_w90main, 'precond', .TRUE.)                            ! MLWF/wannierise (XML precond)
    !> The order of the finite-difference stencil Wannier90 builds its b-shells with.
    !> Its default is 1, the single-shell formula of Marzari and Vanderbilt, whose error
    !> in the k derivative is O(b^2). What feels it are the quantities DEFINED as a
    !> derivative -- A, and above all C and F, which are second derivatives and so carry
    !> that error squared. Behind the sentinel, like everything optional: silence tells
    !> Wannier90 nothing and leaves it on its own default.
    IF (this%fd_order > 0) CALL w90_set_option(wannierlib_w90main, 'higher_order_n', this%fd_order)

    !> The keywords are Wannier90's own, spelled as its .win reader spells them:
    !> w90_set_option only appends the pair to a settings store, and what interprets it
    !> later is the very parser that reads a .win. A name it does not know is an error
    !> there, not a silent no-op, so a typo stops the run.
    IF (this%num_cg_steps >= 0)    CALL w90_set_option(wannierlib_w90main, 'num_cg_steps', this%num_cg_steps)
    IF (this%trial_step > 0.0)     CALL w90_set_option(wannierlib_w90main, 'trial_step', this%trial_step)
    IF (this%fixed_step > 0.0)     CALL w90_set_option(wannierlib_w90main, 'fixed_step', this%fixed_step)
    IF (this%guiding_centres)      CALL w90_set_option(wannierlib_w90main, 'guiding_centres', .TRUE.)
    IF (this%num_guide_cycles >= 0) CALL w90_set_option(wannierlib_w90main, 'num_guide_cycles', this%num_guide_cycles)
    IF (this%num_no_guide_iter >= 0) CALL w90_set_option(wannierlib_w90main, 'num_no_guide_iter', this%num_no_guide_iter)
    IF (this%conv_noise_amp > 0.0) CALL w90_set_option(wannierlib_w90main, 'conv_noise_amp', this%conv_noise_amp)
    IF (this%conv_noise_num >= 0)  CALL w90_set_option(wannierlib_w90main, 'conv_noise_num', this%conv_noise_num)

    !> <transport>: Wannier90's ballistic conductance. one_dim_axis is NOT optional once the
    !> stage is on -- w90_input_setopt below refuses the input without it -- so it is set from
    !> the same block that switches the stage on and cannot be left behind. Only the bulk mode
    !> is offered: it asks for a periodic direction and nothing else, while lcr wants a
    !> lead-conductor-lead cell that no case here has, and an option that cannot be exercised
    !> is a promise the schema cannot keep.
    IF (this%l_transport) THEN
       CALL w90_set_option(wannierlib_w90main, 'transport', .TRUE.)
       CALL w90_set_option(wannierlib_w90main, 'transport_mode', 'bulk')
       SELECT CASE (this%tran_axis)
       CASE (1)
          CALL w90_set_option(wannierlib_w90main, 'one_dim_axis', 'x')
       CASE (2)
          CALL w90_set_option(wannierlib_w90main, 'one_dim_axis', 'y')
       CASE (3)
          CALL w90_set_option(wannierlib_w90main, 'one_dim_axis', 'z')
       CASE DEFAULT
          CALL juDFT_error('wannierlib: <transport> needs axis="x", "y" or "z"', &
                           calledby='set_w90_options')
       END SELECT
       !> The transport stage translates the Wannier centres into one cell before it cuts
       !> the principal layers, and by default it takes the centre of that cell from the
       !> ATOMIC POSITIONS. Wannier90 never gets them from us and cannot: atoms_frac is a
       !> .win BLOCK, and w90_set_option has no overload for an array of strings -- only a
       !> scalar one -- so there is no call that would deliver them. Its own interface
       !> carries a commented-out w90_set_option_b1d, so the gap is known upstream.
       !>
       !> With no atoms the average is a sum of nothing over a count of nothing, and
       !> hamiltonian.F90 divides by atom_data%num_atoms to get NaN centres and a
       !> transmission that is identically zero -- a curve that looks computed and is not.
       !> Naming the centre explicitly is Wannier90's own way out: its reader sets
       !> automatic_translation = .FALSE. as soon as the keyword appears, and the atomic
       !> average is never reached. Half a cell puts the centres in [0,1), the home cell.
       CALL w90_set_option(wannierlib_w90main, 'translation_centre_frac', [0.5, 0.5, 0.5])
       !> Hartree here and eV there, like every other energy that crosses this boundary.
       IF (this%tran_win_min /= this%tran_win_max) THEN
          CALL w90_set_option(wannierlib_w90main, 'tran_win_min', hartree_to_ev_const*this%tran_win_min)
          CALL w90_set_option(wannierlib_w90main, 'tran_win_max', hartree_to_ev_const*this%tran_win_max)
       END IF
       IF (this%tran_energy_step > 0.0) &
          CALL w90_set_option(wannierlib_w90main, 'tran_energy_step', hartree_to_ev_const*this%tran_energy_step)
       !> Hartree and Bohr on this side, electron volts and Angstrom on the other, as for
       !> the windows above and for the cell itself.
       IF (this%tran_hr_cutoff >= 0.0) &
          CALL w90_set_option(wannierlib_w90main, 'hr_cutoff', hartree_to_ev_const*this%tran_hr_cutoff)
       IF (this%tran_dist_cutoff >= 0.0) &
          CALL w90_set_option(wannierlib_w90main, 'dist_cutoff', bohr_to_angstrom_const*this%tran_dist_cutoff)
    END IF

    !> The stem the checkpoint and its witness are named after. It is the same stem the
    !> .amn, .mmn, .eig and .win carry, so that the five files a standalone wannier90.x
    !> reads sit side by side under one name and none of them has to be renamed first.
    seedname = 'WF1'
    CALL w90_input_setopt(wannierlib_w90main, seedname, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_input_setopt failed in wannierlib adapter', calledby='set_w90_options')
  END SUBROUTINE set_w90_options

  !> The neighbour topology Wannier90 chose, read back out. It is available only once the
  !> options are committed, because that is when the k-mesh exists.
  SUBROUTINE get_w90_topology(kpts, nntot_w90, nnkp, gkpb)
    TYPE(t_kpts), INTENT(IN) :: kpts
    INTEGER, INTENT(OUT) :: nntot_w90
    INTEGER, ALLOCATABLE, INTENT(OUT) :: nnkp(:, :)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: gkpb(:, :, :)

    INTEGER :: nn, ierr

    CALL w90_get_nn(wannierlib_w90main, nn, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_get_nn failed in wannierlib adapter', calledby='get_w90_topology')
    nntot_w90 = nn
    IF (nntot_w90 <= 0) CALL juDFT_error('wannierlib invalid nntot from w90_get_nn', calledby='get_w90_topology')

    ALLOCATE(nnkp(kpts%nkptf, nntot_w90), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating nnkp', calledby='get_w90_topology')
    CALL w90_get_nnkp(wannierlib_w90main, nnkp, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_get_nnkp failed in wannierlib adapter', calledby='get_w90_topology')

    ALLOCATE(gkpb(3, kpts%nkptf, nntot_w90), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating gkpb', calledby='get_w90_topology')
    CALL w90_get_gkpb(wannierlib_w90main, gkpb, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_get_gkpb failed in wannierlib adapter', calledby='get_w90_topology')
  END SUBROUTINE get_w90_topology
#endif

  !> Run the Wannier90 library: hand over the eigenvalues, the overlaps M and the projections A,
  !> then disentangle -> project -> wannierise. Returns the resulting gauge factors, which is all
  !> any downstream consumer needs:
  !>    u_matrix  the MLWF gauge U(k)
  !>    u_opt     the disentanglement matrix (the projections after w90 has processed them)
  !> The full Wannier gauge is V(k) = u_opt(k) . u_matrix(k). Evaluating operators in that gauge
  !> is NOT this module's job: the adapter only translates to and from the library.
  SUBROUTINE run_w90(this, cell, kpts, mmn, amn, eig, jspin, n_channels, &
                     u_matrix, u_opt, mw_chk)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, TARGET, INTENT(IN) :: mmn(:, :, :, :)
    COMPLEX, TARGET, INTENT(IN) :: amn(:, :, :)
    REAL, TARGET, INTENT(IN) :: eig(:, :)
    !> Which channel this call is wannierising, and how many there are. Only the checkpoint
    !> name depends on them, and only when there are two.
    INTEGER, INTENT(IN) :: jspin, n_channels
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: u_matrix(:, :, :)   ! (nw,nw,nk) MLWF gauge
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: u_opt(:, :, :)      ! (nb,nw,nk) disentangled projections
    !> The overlaps as Wannier90 leaves them: the num_wann corner of the buffer it was
    !> given, which by then holds V^dagger(k) M V(k_b) in the final gauge and not the
    !> Bloch overlaps that went in. It is the same block the checkpoint carries, so the
    !> two can be compared against each other while both exist.
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mw_chk(:, :, :, :)  ! (nw,nw,nntot,nk_loc)

    INTEGER :: ierr, num_kpts
    CHARACTER(LEN=32) :: seedname
    COMPLEX, ALLOCATABLE :: mmn_local(:, :, :, :)
    REAL, ALLOCATABLE, TARGET :: eig_ev(:, :)         ! eigenvalues in eV for the w90 library (banner honesty)

    ! allocate the INTENT(OUT) results before any early return, so a caller can never pass
    ! unallocated arrays on to the matrix-element layer
    ALLOCATE(u_matrix(0, 0, 0), u_opt(0, 0, 0), mw_chk(0, 0, 0, 0))
    IF (.NOT.this%l_wannierize) RETURN
    DEALLOCATE(u_matrix, u_opt, mw_chk)
    num_kpts = SIZE(amn, 3)

    !> One checkpoint per spin channel. Collinearly the two channels are two separate
    !> wannierisations with two gauges, and both used to write the one name the seedname
    !> carries, so the second silently overwrote the first and only one channel was ever
    !> handed out. The name is set here, before any checkpoint is written, because
    !> set_w90_options runs once for the whole calculation.
    !>
    !> A single channel is WF1, the same stem the one set of Bloch files carries; two
    !> channels are WF1 and WF2, one per spin, as those files already are.
    IF (n_channels > 1) THEN
      WRITE (seedname, '(a,i0)') 'WF', jspin
      !> Assigned, not set through w90_input_setopt: that call PROCESSES the option store
      !> and empties it on the way out ("must be empty here; emptied by w90_input_setopt
      !> call", library_interface.F90:458), so a second call fails. The field is the one
      !> write_chkpt takes to build the file name -- it passes common_data%seedname at
      !> library_extra.F90:393 -- so this is the whole of what has to change, and it
      !> changes nothing else about w90's state.
      wannierlib_w90main%seedname = seedname
    END IF

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
      !> Optional, and tried first rather than instead: it fills u_opt with a subspace that
      !> holds the same number of states from each spin channel at every k. Wannier90 keeps
      !> the num_wann largest eigenvalues globally, which at an occasional k takes one state
      !> too many from one channel -- and because the Wannier index is shared by every k, a
      !> single such point makes the spin label inconsistent everywhere and forces the gauge
      !> to mix. Whenever the premise does not hold the routine writes why and returns
      !> .FALSE., and the ordinary disentanglement runs as before.
      CALL w90_disentangle(wannierlib_w90main, oUnit, oUnit, ierr)
      IF (ierr /= 0) CALL juDFT_error('w90_disentangle failed in wannierlib adapter', calledby='run_w90')
    END IF

    CALL w90_project_overlap(wannierlib_w90main, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_project_overlap failed in wannierlib adapter', calledby='run_w90')

    !> The point a reader resumes the minimisation from, so it is written here rather than
    !> straight after the disentanglement: this is the state w90_wannierise is about to be
    !> handed. Every rank calls it -- the writer gathers the distributed matrices itself.
    IF (this%export%w90) CALL wannierlib_write_chkpt('postdis')

    CALL w90_wannierise(wannierlib_w90main, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_wannierise failed in wannierlib adapter', calledby='run_w90')

    !> Again at the end, to the same file: the minimised gauge replaces the projected one,
    !> which is what wannier90.x leaves behind and what its restart expects to find.
    IF (this%export%w90) CALL wannierlib_write_chkpt('postwann')

    !> Transport last, because it reads the gauge: tran_main builds its own real-space
    !> Hamiltonian from u_matrix, u_matrix_opt and eigval, which are exactly what the calls
    !> above have already set, so nothing further is handed to it. That H(R) is Wannier90's
    !> own, with its own Wigner-Seitz convention, and is NOT the one <operatorsRs> writes --
    !> two objects of the same name that must not be read into each other.
    !>
    !> Every rank calls it although tran_main is serial: w90_transport runs it on rank 0 and
    !> then synchronises the error state across the communicator, so a rank that skipped the
    !> call would leave that synchronisation half done.
    IF (this%l_transport) THEN
      CALL w90_transport(wannierlib_w90main, oUnit, oUnit, ierr)
      IF (ierr /= 0) CALL juDFT_error('w90_transport failed in wannierlib adapter', calledby='run_w90')
    END IF

    !> Taken here and not beside the checkpoint write, because what the caller is owed is
    !> what Wannier90 ended with: every stage above rotates this buffer, and a copy made
    !> earlier would be the gauge of a stage that has since been superseded.
    ALLOCATE(mw_chk(this%num_wann, this%num_wann, SIZE(mmn_local, 3), SIZE(mmn_local, 4)), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating the overlap handback', &
                                    calledby='run_w90')
    mw_chk = mmn_local(:this%num_wann, :this%num_wann, :, :)

    call timestop('run_w90')
#endif
  END SUBROUTINE run_w90

#ifdef CPP_WANNLIB_API
  !> Ask Wannier90 to write its own checkpoint. The file and its layout are Wannier90's;
  !> what FLEUR decides is when to ask. 'label' is the position it records -- 'postdis'
  !> before the minimisation, 'postwann' after -- and both go to the same <seedname>.chk,
  !> which is how wannier90.x behaves and how its restart reads the position back.
  !>
  !> Every rank calls it: write_chkpt gathers the distributed overlaps itself, so a rank
  !> that skipped the call would leave that gather half done.
  SUBROUTINE wannierlib_write_chkpt(label)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER :: ierr
    CALL write_chkpt(wannierlib_w90main, label, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('write_chkpt failed in wannierlib adapter', &
                                    calledby='wannierlib_write_chkpt')
  END SUBROUTINE wannierlib_write_chkpt
#endif


  !> Export the coarse-mesh b-shell / neighbour information that the Wannier90 kmesh setup
  !> produced, as a plain t_wgauge_bmesh. This is the ONE piece of Wannier90 state the
  !> matrix-element layer needs (for the position / Berry-connection operator), and handing it
  !> over like this keeps wannierlib_w90main private to this module and keeps m_melem_* free of
  !> the Wannier90 library. Call after run_w90 (the centres are only meaningful once wannierised).
  SUBROUTINE wannierlib_get_bmesh(this, kpts, bmesh)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    !> INOUT and not OUT: the topology is already in here, put there before the
    !> wannierisation ran, and OUT would default-initialise it away. What this adds is the
    !> part only Wannier90 knows -- the shell weights and the b vectors it chose.
    TYPE(t_wgauge_bmesh), INTENT(INOUT) :: bmesh
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
      CALL gauge_centres(bmesh%centres)
    END IF
#else
    bmesh%nntot = 0   ! no Wannier90 library -> no b-mesh; the position operator is skipped
#endif
  END SUBROUTINE wannierlib_get_bmesh

  !> The centres and the spreads of the gauge Wannier90 produced. The spreads are optional
  !> because one caller only places the position operator and has no use for them.
  SUBROUTINE gauge_centres(centres, spreads)
    REAL, INTENT(OUT) :: centres(:, :)
    REAL, INTENT(OUT), OPTIONAL :: spreads(:)
#ifdef CPP_WANNLIB_API
    CALL w90_get_centres(wannierlib_w90main, centres)
    IF (PRESENT(spreads)) CALL w90_get_spreads(wannierlib_w90main, spreads)
#endif
  END SUBROUTINE gauge_centres

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
    CALL gauge_centres(wann_centres, wann_spreads)

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
