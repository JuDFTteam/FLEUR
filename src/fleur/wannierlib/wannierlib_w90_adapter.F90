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
  USE m_wannierlib_disentangle_spin, ONLY: wannierlib_disentangle_spin
#ifdef CPP_WANNLIB_API
  USE w90_library, ONLY : lib_common_type, w90_set_comm, w90_set_option, w90_input_setopt, &
                          w90_get_nn, w90_get_nnkp, w90_get_gkpb, w90_set_eigval, &
                          w90_set_u_opt, w90_set_m_local, w90_set_u_matrix, &
                          w90_disentangle, w90_project_overlap, w90_wannierise, w90_create_kmesh, &
                          w90_get_centres, w90_get_spreads
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
  PUBLIC :: init_w90, run_w90, wannierlib_get_bmesh, report_w90
#ifdef CPP_WANNLIB_API
  TYPE(lib_common_type), SAVE, PUBLIC :: wannierlib_w90main
  !> One Wannier90 instance per spin channel, half as wide, for the spin-balanced mode. They
  !> are separate objects and not one reused twice because the width is fixed at setup: the
  !> k-mesh, the b-shells and every array Wannier90 dimensions from num_wann are built inside
  !> w90_input_setopt. Three live instances are legal -- the 4.x encapsulation refactor moved
  !> all state into common_data, and there is not one module-level SAVE variable left in
  !> library_interface, wannierise, overlap or kmesh.
  TYPE(lib_common_type), SAVE :: wannierlib_w90chan(2)
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
    INTEGER :: ich
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
 

    !> In the spin-balanced mode FLEUR does the disentanglement itself and hands over a
    !> problem that is already num_wann wide, so Wannier90 is told so from the start: it then
    !> takes the supported nb == nw path and never needs its own disentanglement. The energy
    !> windows are deliberately not passed in that mode -- they were already applied on this
    !> side, and repeating them on a manifold that no longer has an energy axis is meaningless.
    IF (this%dis_spin_balanced) THEN
      CALL w90_set_option(wannierlib_w90main, 'num_bands', this%num_wann)
    ELSE
      CALL w90_set_option(wannierlib_w90main, 'num_bands', this%num_bands)
    END IF
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
    !> How many neighbour shells w90 may try before it gives up on the B1 condition. Its own
    !> default is 36, which is not enough for a film: the mesh is dense in plane and a single
    !> point in z, so shell after shell is rejected by the SVD test and 36 runs out at
    !> 14x14x1. The classic path has written 200 into the .win since it was added
    !> (wann_wan90prep.F:387), and 200 is what the group uses; passing nothing left the
    !> library mode on w90's default and made it fail where the classic would not.
    CALL w90_set_option(wannierlib_w90main, 'search_shells', WANNIERLIB_SEARCH_SHELLS)
    IF (.NOT. this%dis_spin_balanced) &
      CALL w90_set_option(wannierlib_w90main, 'dis_win_min', hartree_to_ev_const*this%dis_win_min)
    IF (.NOT. this%dis_spin_balanced) &
      CALL w90_set_option(wannierlib_w90main, 'dis_win_max', hartree_to_ev_const*this%dis_win_max)
    IF (.NOT. this%dis_spin_balanced) &
      CALL w90_set_option(wannierlib_w90main, 'dis_froz_min', hartree_to_ev_const*this%dis_froz_min)
    IF (.NOT. this%dis_spin_balanced) &
      CALL w90_set_option(wannierlib_w90main, 'dis_froz_max', hartree_to_ev_const*this%dis_froz_max)
    IF (this%dis_mix_ratio > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_mix_ratio', this%dis_mix_ratio)
    IF (this%dis_conv_tol > 0.0) CALL w90_set_option(wannierlib_w90main, 'dis_conv_tol', this%dis_conv_tol)  ! disentanglement (XML disConvTol)
    IF (this%conv_tol > 0.0)     CALL w90_set_option(wannierlib_w90main, 'conv_tol', this%conv_tol)          ! MLWF/wannierise (XML wannConvTol)
    !> Only set when asked: Wannier90's own default is .FALSE., so staying silent keeps
    !> every existing run byte-identical.
    !> Projectability disentanglement: only set when asked, so silence leaves W90 on its own
    !> energy-window default and every existing run byte-identical. The two thresholds go with
    !> it -- W90 validates both to [0,1] and would abort on anything else.
    IF (this%dis_froz_proj) THEN
       CALL w90_set_option(wannierlib_w90main, 'dis_froz_proj', .TRUE.)
       CALL w90_set_option(wannierlib_w90main, 'dis_proj_min', this%dis_proj_min)
       CALL w90_set_option(wannierlib_w90main, 'dis_proj_max', this%dis_proj_max)
    END IF
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

    !> The per-channel instances are built here, where the cell, the mesh and the atoms are
    !> in scope, and used in run_w90. Only when the mode is on: an unused instance would
    !> still print its banner and build its own k-mesh.
    IF (this%dis_spin_balanced) THEN
      DO ich = 1, 2
        CALL setup_channel(wannierlib_w90chan(ich), this, atoms, cell, kpts, fmpi, spinors, &
                           distk, ich)
      END DO
    END IF
    call timestop('init_w90')
#endif
  END SUBROUTINE init_w90

!> Guarded whole: unlike the routines around it, this one names lib_common_type in its
!> SIGNATURE, and a dummy-argument declaration cannot be handled by the #ifndef-and-error
!> pattern the rest of the file uses. Without the module API of Wannier90 -- the CI image
!> still carries the 1.2 library, where w90_library.mod does not exist -- the type is not
!> declared and the file does not compile. Only reached from the guarded branch of init_w90.
#ifdef CPP_WANNLIB_API
  !> One Wannier90 instance for a single spin channel: num_wann/2 functions, and num_bands
  !> equal to it because the subspace was already chosen on this side. The option list is
  !> deliberately shorter than the main one -- there is no disentanglement here, so the energy
  !> windows and every dis_* setting would be describing a stage that does not run. What it
  !> does share is the geometry and the mesh, which is what fixes the b-shells: the two
  !> channels must end up with the same neighbours and weights as the main instance, since the
  !> overlaps handed to them were built with those.
  SUBROUTINE setup_channel(w90, this, atoms, cell, kpts, fmpi, spinors, distk, ich)
    TYPE(lib_common_type), INTENT(INOUT) :: w90
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_mpi), INTENT(IN) :: fmpi
    LOGICAL, INTENT(IN) :: spinors
    INTEGER, INTENT(IN) :: distk(:)
    INTEGER, INTENT(IN) :: ich

    INTEGER :: ierr, nper
    CHARACTER(LEN=32) :: seedname

    nper = this%num_wann/2
    CALL w90_set_comm(w90, fmpi%mpi_comm)
    CALL w90_set_option(w90, 'num_bands', nper)
    CALL w90_set_option(w90, 'num_wann', nper)
    CALL w90_set_option(w90, 'num_kpts', kpts%nkptf)
    CALL w90_set_option(w90, 'num_atoms', atoms%nat)
    CALL w90_set_option(w90, 'mp_grid', kpts%nkpt3)
    CALL w90_set_option(w90, 'kpoints', kpts%bkf)
    CALL w90_set_option(w90, 'unit_cell_cart', bohr_to_angstrom_const*cell%amat)
    CALL w90_set_option(w90, 'distk', distk)
    CALL w90_set_option(w90, 'gamma_only', .FALSE.)
    CALL w90_set_option(w90, 'spinors', spinors)
    CALL w90_set_option(w90, 'dump_inputs', .FALSE.)
    CALL w90_set_option(w90, 'search_shells', WANNIERLIB_SEARCH_SHELLS)
    IF (this%num_iter > 0) CALL w90_set_option(w90, 'num_iter', this%num_iter)
    IF (this%conv_tol > 0.0) CALL w90_set_option(w90, 'conv_tol', this%conv_tol)
    IF (this%precond) CALL w90_set_option(w90, 'precond', .TRUE.)

    WRITE(seedname, '(a,i0)') 'fleur_wlib_channel', ich
    CALL w90_input_setopt(w90, seedname, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_input_setopt failed for a spin channel', &
                                    calledby='setup_channel')

    !> The b-shells are NOT built by w90_input_setopt. On the main instance they get built as
    !> a side effect of w90_get_nn, which calls this when setup_complete is still false; here
    !> it is called for what it does. Without it the instance reaches the wannierisation with
    !> no b vectors and no weights, the initial spread comes out as exactly zero, and ZHEEV
    !> fails inside wann_main with info=8.
    CALL w90_create_kmesh(w90, oUnit, oUnit, ierr)
    IF (ierr /= 0) CALL juDFT_error('w90_create_kmesh failed for a spin channel', &
                                    calledby='setup_channel')
  END SUBROUTINE setup_channel
#endif

  !> Which Wannier columns belong to each spin channel. Recomputed from proj_spin rather than
  !> handed down from the disentanglement, so the places that need it cannot drift apart, and
  !> read as an index list rather than as a halving: the projections do come out in blocks
  !> today, but nothing in this routine depends on that staying true.
  SUBROUTINE channel_columns(this, idx)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: idx(:, :)       !> (num_wann/2, 2)
    INTEGER :: iw, ich, n(2)
    n = 0
    DO iw = 1, this%num_wann
      ich = MERGE(1, 2, this%proj_spin(iw) > 0)
      n(ich) = n(ich) + 1
      IF (n(ich) > SIZE(idx, 1)) CALL juDFT_error( &
        'wannierlib: more projections in one spin channel than half the Wannier functions', &
        calledby='channel_columns')
      idx(n(ich), ich) = iw
    END DO
    IF (ANY(n /= SIZE(idx, 1))) CALL juDFT_error( &
      'wannierlib: the two spin channels do not hold the same number of projections', &
      calledby='channel_columns')
  END SUBROUTINE channel_columns

  !> Run the Wannier90 library: hand over the eigenvalues, the overlaps M and the projections A,
  !> then disentangle -> project -> wannierise. Returns the resulting gauge factors, which is all
  !> any downstream consumer needs:
  !>    u_matrix  the MLWF gauge U(k)
  !>    u_opt     the disentanglement matrix (the projections after w90 has processed them)
  !> The full Wannier gauge is V(k) = u_opt(k) . u_matrix(k). Evaluating operators in that gauge
  !> is NOT this module's job: the adapter only translates to and from the library.
  SUBROUTINE run_w90(this, cell, kpts, mmn, amn, eig, irank, u_matrix, u_opt, s0)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, TARGET, INTENT(IN) :: mmn(:, :, :, :)
    COMPLEX, TARGET, INTENT(IN) :: amn(:, :, :)
    REAL, TARGET, INTENT(IN) :: eig(:, :)
    INTEGER, INTENT(IN) :: irank
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: u_matrix(:, :, :)   ! (nw,nw,nk) MLWF gauge
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: u_opt(:, :, :)      ! (nb,nw,nk) disentangled projections
    !> The spin operator in the Bloch basis, (nb,nb,3,nk). Only the spin-balanced mode uses
    !> it, and only to label each band; absent means that mode has to fall back.
    COMPLEX, INTENT(IN), OPTIONAL :: s0(:, :, :, :)

    INTEGER :: ierr, num_kpts
    COMPLEX, ALLOCATABLE :: mmn_local(:, :, :, :)
    REAL, ALLOCATABLE, TARGET :: eig_ev(:, :)         ! eigenvalues in eV for the w90 library (banner honesty)
    LOGICAL :: l_balanced
    INTEGER :: nw, nnt, ik, nn, jk, iw
    INTEGER :: nper, ich, i, j
    INTEGER, ALLOCATABLE :: idx(:, :)
    COMPLEX, ALLOCATABLE, TARGET :: mred(:, :, :, :), ared(:, :, :)
    REAL, ALLOCATABLE, TARGET :: ered(:, :)
    !> TARGET, and kept alive until the wannierisation of their channel has returned:
    !> w90_set_* stores pointers to these, it does not copy them.
    COMPLEX, ALLOCATABLE, TARGET :: mch(:, :, :, :), ach(:, :, :), uch(:, :, :)
    REAL, ALLOCATABLE, TARGET :: ech(:, :)

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
    IF (this%dis_spin_balanced) THEN
      !> FLEUR chooses the subspace, and what Wannier90 receives is already num_wann wide.
      IF (SIZE(mmn_local, 4) /= num_kpts) CALL juDFT_error( &
         "wannierlib: the spin-balanced disentanglement needs every k-point on one rank", &
         hint="run this case on one MPI rank, or drop spinBalanced", calledby="run_w90")
      CALL wannierlib_disentangle_spin( &
        hartree_to_ev_const*this%dis_win_min, hartree_to_ev_const*this%dis_win_max, &
        hartree_to_ev_const*this%dis_froz_min, hartree_to_ev_const*this%dis_froz_max, &
        wannierlib_w90main%kmesh_info%nntot, wannierlib_w90main%kmesh_info%nnlist, &
        wannierlib_w90main%kmesh_info%wb, eig_ev, mmn_local, amn, this%proj_spin, &
        u_opt, l_balanced, s0)
      IF (.NOT. l_balanced) CALL juDFT_error( &
         "wannierlib: the spin-balanced disentanglement does not apply to this case", &
         hint="the reason is in the output just above; drop spinBalanced to use the "// &
              "ordinary one", calledby="run_w90")

      nw = this%num_wann
      nnt = wannierlib_w90main%kmesh_info%nntot
      ALLOCATE(mred(nw, nw, nnt, num_kpts), ared(nw, nw, num_kpts), ered(nw, num_kpts))
      DO ik = 1, num_kpts
        DO nn = 1, nnt
          jk = wannierlib_w90main%kmesh_info%nnlist(ik, nn)
          mred(:, :, nn, ik) = MATMUL(CONJG(TRANSPOSE(u_opt(:, :, ik))), &
                                      MATMUL(mmn_local(:, :, nn, ik), u_opt(:, :, jk)))
        END DO
        ared(:, :, ik) = MATMUL(CONJG(TRANSPOSE(u_opt(:, :, ik))), amn(:, :, ik))
        DO iw = 1, nw
          ered(iw, ik) = REAL(DOT_PRODUCT(u_opt(:, iw, ik), eig_ev(:, ik)*u_opt(:, iw, ik)))
        END DO
      END DO
      !> Nothing is handed to the main instance here: in this mode it does not wannierise.
      !> The reduced problem goes to the per-channel instances further down.
    ELSE
      CALL w90_set_eigval(wannierlib_w90main, eig_ev)
      CALL w90_set_m_local(wannierlib_w90main, mmn_local)
      CALL w90_set_u_opt(wannierlib_w90main, u_opt)
      CALL w90_set_u_matrix(wannierlib_w90main, u_matrix)
    END IF

    IF (this%num_bands > this%num_wann) THEN
      !> Optional, and tried first rather than instead: it fills u_opt with a subspace that
      !> holds the same number of states from each spin channel at every k. Wannier90 keeps
      !> the num_wann largest eigenvalues globally, which at an occasional k takes one state
      !> too many from one channel -- and because the Wannier index is shared by every k, a
      !> single such point makes the spin label inconsistent everywhere and forces the gauge
      !> to mix. Whenever the premise does not hold the routine writes why and returns
      !> .FALSE., and the ordinary disentanglement runs as before.
      !> Not in the balanced mode: there Wannier90 was set up with num_bands = num_wann
      !> and has nothing left to disentangle.
      IF (.NOT. this%dis_spin_balanced) THEN
        CALL w90_disentangle(wannierlib_w90main, oUnit, oUnit, ierr)
        IF (ierr /= 0) CALL juDFT_error('w90_disentangle failed in wannierlib adapter', calledby='run_w90')
      END IF
    END IF

    IF (this%dis_spin_balanced) THEN
      !> One wannierisation per spin channel. Restricting the minimisation itself is not an
      !> option -- Wannier90 is read, never patched -- but it does not need to be: M_red is
      !> block diagonal to machine precision (measured, 0.0000 % of cross weight on three
      !> magnetisation axes), so the two channels are two independent problems of num_wann/2
      !> functions and the mixing simply has nowhere to happen.
      !>
      !> WHY THIS AND NOT JUST A BLOCKED SUBSPACE. u_opt already comes out blocked by
      !> construction, and so do A_red and M_red; the minimisation therefore starts inside the
      !> block manifold and, in exact arithmetic, the gradient cannot take it out. Wannier90
      !> leaves it anyway: measured on bcc Fe, the final u_matrix carries 21.6 % of its weight
      !> in the cross blocks with the moment along x and 22.3 % along y, and the Wannier
      !> functions stop being spin eigenstates (|<s>| drops to 0.57 on average, 0.01 at worst).
      !> With the moment along z it stays blocked to 0.0000 % and the functions are exactly
      !> pure -- which is why the case looked solved for as long as only z was tested.
      !>
      !> The cost is not a cost. Minimising inside the manifold reaches Omega_tot 21.421 along
      !> x against 21.599 for the unconstrained run: the blocked optimum is the BETTER one
      !> there. An earlier estimate of +0.717 came from projecting a mixed solution onto the
      !> blocks, which is not the same thing as minimising within them and is only an upper
      !> bound.
      nper = this%num_wann/2
      ALLOCATE(idx(nper, 2))
      CALL channel_columns(this, idx)
      ALLOCATE(mch(nper, nper, nnt, num_kpts), ach(nper, nper, num_kpts), &
               ech(nper, num_kpts), uch(nper, nper, num_kpts))
      DO ich = 1, 2
        DO ik = 1, num_kpts
          DO nn = 1, nnt
            DO j = 1, nper
              DO i = 1, nper
                mch(i, j, nn, ik) = mred(idx(i, ich), idx(j, ich), nn, ik)
              END DO
            END DO
          END DO
          DO j = 1, nper
            DO i = 1, nper
              ach(i, j, ik) = ared(idx(i, ich), idx(j, ich), ik)
            END DO
          END DO
          DO i = 1, nper
            ech(i, ik) = ered(idx(i, ich), ik)
          END DO
        END DO
        uch = CMPLX(0.0, 0.0)
        CALL w90_set_eigval(wannierlib_w90chan(ich), ech)
        CALL w90_set_m_local(wannierlib_w90chan(ich), mch)
        CALL w90_set_u_opt(wannierlib_w90chan(ich), ach)
        CALL w90_set_u_matrix(wannierlib_w90chan(ich), uch)
        CALL w90_project_overlap(wannierlib_w90chan(ich), oUnit, oUnit, ierr)
        IF (ierr /= 0) CALL juDFT_error('w90_project_overlap failed for a spin channel', &
                                        calledby='run_w90')
        CALL w90_wannierise(wannierlib_w90chan(ich), oUnit, oUnit, ierr)
        IF (ierr /= 0) CALL juDFT_error('w90_wannierise failed for a spin channel', &
                                        calledby='run_w90')
        !> Straight back into the block it came from. What stays outside the two blocks is
        !> the zero u_matrix was allocated with, and that is not a convention: there is
        !> nothing to put there, the channels do not talk to each other.
        DO ik = 1, num_kpts
          DO j = 1, nper
            DO i = 1, nper
              u_matrix(idx(i, ich), idx(j, ich), ik) = uch(i, j, ik)
            END DO
          END DO
        END DO
      END DO
      WRITE (oUnit, '(a,i0,a)') 'wannierlib: wannierised each spin channel on its own, ', &
        nper, ' Wannier functions in each'
    ELSE
      CALL w90_project_overlap(wannierlib_w90main, oUnit, oUnit, ierr)
      IF (ierr /= 0) CALL juDFT_error('w90_project_overlap failed in wannierlib adapter', calledby='run_w90')

      CALL w90_wannierise(wannierlib_w90main, oUnit, oUnit, ierr)
      IF (ierr /= 0) CALL juDFT_error('w90_wannierise failed in wannierlib adapter', calledby='run_w90')
    END IF

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
      CALL gauge_centres(this, bmesh%centres)
    END IF
#else
    bmesh%nntot = 0   ! no Wannier90 library -> no b-mesh; the position operator is skipped
#endif
  END SUBROUTINE wannierlib_get_bmesh

  !> The centres and spreads of whichever run actually produced the gauge. In the balanced
  !> mode that is the two per-channel instances, laid back out in Wannier index order -- not
  !> the main one, which never wannierises there and would hand back whatever its arrays
  !> happened to be initialised to.
  SUBROUTINE gauge_centres(this, centres, spreads)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    REAL, INTENT(OUT) :: centres(:, :)
    REAL, INTENT(OUT), OPTIONAL :: spreads(:)
#ifdef CPP_WANNLIB_API
    INTEGER :: ich, i, nper
    INTEGER, ALLOCATABLE :: idx(:, :)
    REAL, ALLOCATABLE :: c(:, :), s(:)

    IF (.NOT. this%dis_spin_balanced) THEN
      CALL w90_get_centres(wannierlib_w90main, centres)
      IF (PRESENT(spreads)) CALL w90_get_spreads(wannierlib_w90main, spreads)
      RETURN
    END IF

    nper = this%num_wann/2
    ALLOCATE(idx(nper, 2), c(3, nper), s(nper))
    CALL channel_columns(this, idx)
    DO ich = 1, 2
      CALL w90_get_centres(wannierlib_w90chan(ich), c)
      CALL w90_get_spreads(wannierlib_w90chan(ich), s)
      DO i = 1, nper
        centres(:, idx(i, ich)) = c(:, i)
        IF (PRESENT(spreads)) spreads(idx(i, ich)) = s(i)
      END DO
    END DO
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
    CALL gauge_centres(this, wann_centres, wann_spreads)

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
