!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_w90_adapter
  USE m_constants, ONLY : oUnit, namat_const, hartree_to_ev_const, bohr_to_angstrom_const
  USE m_juDFT
  USE m_wannierlib_ft, ONLY : wannierlib_ft_to_real, wannierlib_ws_vectors, wannierlib_ft_to_real_reduce
  USE m_xmlOutput
  USE m_types_atoms
  USE m_types_cell
  USE m_types_kpts
  USE m_types_mpi
  USE m_types_wannierlib
  USE m_wannierlib_interpolate
  USE m_wannierlib_interpolate_op
  USE m_wannierlib_interpolate_velocity
  USE m_wannierlib_interpolate_current
  USE m_wannierlib_interpolate_eigenstates
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

  SUBROUTINE run_w90(this, cell, kpts, mmn, amn, eig, irank, &
                     s0_loc, l0_loc, soc0_loc, soc4_loc, s0pa_loc, distk, mpi_comm, wf_channel, spin_suffix, l0col_loc, v_out)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, TARGET, INTENT(IN) :: mmn(:, :, :, :)
    COMPLEX, TARGET, INTENT(IN) :: amn(:, :, :)
    REAL, TARGET, INTENT(IN) :: eig(:, :)
    INTEGER, INTENT(IN) :: irank
    ! per-rank coarse operator slices (ascending global-k order = gk_loc); all interpolation
    ! and current drivers consume these via the distributed FT-reduce (no full-mesh gather).
    COMPLEX, INTENT(IN) :: s0_loc(:, :, :, :), l0_loc(:, :, :, :, :), soc4_loc(:, :, :, :) ! spin / orbital / 2x2 SOC blocks
    COMPLEX, INTENT(IN) :: soc0_loc(:, :, :, :)       ! (nb,nb,1,nk_loc) per-rank SOC slice (soc interpolation reduce)
    COMPLEX, INTENT(IN) :: s0pa_loc(:, :, :, :, :)    ! (nb,nb,3,nat,nk_loc) per-rank per-atom MT spin slice (per-atom interpolation reduce)
    INTEGER, INTENT(IN) :: distk(:), mpi_comm
    INTEGER, INTENT(IN), OPTIONAL :: wf_channel            ! collinear jspins=2 spin channel (1/2) for the operators_r WFn seedname; default 1
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spin_suffix   ! '_spin1'/'_spin2' for collinear jspins=2; empty otherwise
    COMPLEX, INTENT(IN), OPTIONAL :: l0col_loc(:, :, :, :)  ! (nb,nb,3,nk_loc) collinear per-channel Bloch orbital L (operators_r)
    COMPLEX, INTENT(OUT), OPTIONAL :: v_out(:, :, :)        ! (num_bands,num_wann,nkptf) gauge V=u_opt.u_matrix (collinear combined spin)

    INTEGER :: ierr,num_kpts,iop,nbnd_c,nat_c,nkc,wf_ch
    LOGICAL :: l_collinear
    INTEGER :: idom, ndom, nkl_c, jkl, aw_nrpts
    INTEGER, ALLOCATABLE :: gk_loc(:)                 ! global-k indices of this rank's coarse slice (interp reduce)
    INTEGER, ALLOCATABLE :: aw_irvec(:, :), aw_ndegen(:)  ! Wigner-Seitz R-mesh of the reduced Berry connection
    CHARACTER(LEN=8) :: dkind(3), dsuf(3)
    CHARACTER(LEN=16) :: ssfx
    LOGICAL :: lex
    COMPLEX, ALLOCATABLE :: u_matrix(:, :, :), mmn_local(:, :, :, :), amn_local(:, :, :)
    COMPLEX, ALLOCATABLE :: aw_r(:, :, :, :)          ! (nw,nw,nrpts,3) reduced Wannier Berry connection A^(W)_alpha(R)
    COMPLEX, ALLOCATABLE :: hamk_r(:, :, :)           ! (nw,nw,nk) Wannier-gauge Hamiltonian (for H(R) export)
    LOGICAL :: l_hr_done, l_ar_done                   ! tight-binding H(R)/A(R) written once (domain-independent)
    REAL, ALLOCATABLE, TARGET :: eig_ev(:, :)         ! eigenvalues in eV for the w90 library (banner honesty)


    IF (.NOT.this%l_wannierize) RETURN
    num_kpts = SIZE(amn, 3)
    ssfx = ''
    IF (PRESENT(spin_suffix)) ssfx = TRIM(spin_suffix)
    wf_ch = 1
    IF (PRESENT(wf_channel)) wf_ch = wf_channel
    l_collinear = (LEN_TRIM(ssfx) > 0)   ! collinear jspins=2 -> per-channel operators_r (WF1/WF2)

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
  ALLOCATE(amn_local(SIZE(amn,1), SIZE(amn,2), SIZE(amn,3)), stat=ierr)
  IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating amn buffer', calledby='run_w90')
  mmn_local = mmn
  amn_local = amn
    !TODO u_matrix should be stored (in results)  
    ALLOCATE(eig_ev(SIZE(eig,1), SIZE(eig,2)), stat=ierr)
    IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating eig_ev buffer', calledby='run_w90')
    eig_ev = hartree_to_ev_const * eig   ! w90 library expects eV; FLEUR eig are Hartree
    CALL w90_set_eigval(wannierlib_w90main, eig_ev)
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

    ! export the full gauge V(k) = u_opt(k) . u_matrix(k) (disentangled + MLWF), needed by the
    ! collinear combined spin operator which rotates the cross-spin overlap with both channels' V.
    IF (PRESENT(v_out)) THEN
      DO iop = 1, num_kpts
        v_out(:, :, iop) = MATMUL(amn_local(:, :, iop), u_matrix(:, :, iop))
      END DO
    END IF

    ! global-k indices owned by this rank, ascending order -> matches the per-rank coarse slices
    ! s0_loc/l0_loc/soc0_loc/s0pa_loc (built in the same distk order); used by the distributed
    ! operator interpolation reduce below.
    nkl_c = COUNT(distk == irank); ALLOCATE(gk_loc(nkl_c)); jkl = 0
    DO iop = 1, SIZE(distk)
      IF (distk(iop) == irank) THEN; jkl = jkl + 1; gk_loc(jkl) = iop; END IF
    END DO

    ! ---- opt-in output domains (<path>/<plane>/<grid>); none declared -> legacy single pass ----
    ! order matters: generated domains (plane/grid) overwrite kpts_interpol and are renamed;
    ! the unsuffixed path/legacy domain runs LAST so its base-named output is not clobbered
    ! and it restores the user's original kpts_interpol before interpolating.
    ndom = 0
    IF (this%l_dom_plane) THEN; ndom=ndom+1; dkind(ndom)='plane'; dsuf(ndom)='_plane'; END IF
    IF (this%l_dom_grid)  THEN; ndom=ndom+1; dkind(ndom)='grid';  dsuf(ndom)='_grid';  END IF
    IF (this%l_dom_path)  THEN; ndom=ndom+1; dkind(ndom)='path';  dsuf(ndom)='';       END IF
    IF (ndom == 0) THEN; ndom=1; dkind(1)='legacy'; dsuf(1)=''; END IF
    ! back up a user-provided kpts_interpol that a generated (plane/grid) domain would overwrite
    IF (irank==0 .AND. (this%l_dom_plane .OR. this%l_dom_grid)) THEN
      INQUIRE(file='kpts_interpol', exist=lex)
      IF (lex) CALL wl_shell('cp -f kpts_interpol .kpts_interpol_userbak')
    END IF

    l_hr_done = .FALSE.; l_ar_done = .FALSE.   ! (legacy flags; operators_r now owns H(R)/A(R))
    ! real-space operator export (Fourier step 3, standalone format) -- once, before interpolation
    CALL wannierlib_write_operators_r(this, cell, kpts, eig, u_matrix, amn_local, &
                                      s0_loc, l0_loc, soc4_loc, distk, mpi_comm, mmn, irank, wf_ch, l_collinear, l0col_loc)
    DO idom = 1, ndom
    IF (irank==0) CALL wannierlib_write_domain_kpts(this, TRIM(dkind(idom)))

    ! Wannier-gauge interpolation: dispatch by looping over the requested operator list.
    ! Each operator supplies its own per-rank Bloch slice on the coarse mesh (s0/l0/soc0_loc);
    ! steps (2)-(5) are the shared generic driver. u_matrix = MLWF gauge, amn_local = u_opt.
    DO iop = 1, this%n_ops
      SELECT CASE (TRIM(this%op_name(iop)))
      CASE ('hamiltonian')
        CALL wannierlib_interpolate(this, cell, kpts, eig, u_matrix, amn_local, irank)
      CASE ('spin')
        ! total spin (MT-sum + interstitial): via the generic operator driver (3 comps)
        IF (this%op_total(iop) == 1) &
          CALL wannierlib_interpolate_operator(this, cell, kpts, eig, u_matrix, amn_local, &
                                               s0_loc, gk_loc, 3, 'bands_wann_spin.dat', irank, mpi_comm)
        ! per-atom (site-resolved) muffin-tin spin moment: 3*nat components in one file
        IF (this%op_peratom(iop) == 1) THEN
          nbnd_c = SIZE(s0pa_loc, 1); nat_c = SIZE(s0pa_loc, 4); nkc = SIZE(s0pa_loc, 5)
          CALL wannierlib_interpolate_operator(this, cell, kpts, eig, u_matrix, amn_local, &
                                               RESHAPE(s0pa_loc, (/nbnd_c, nbnd_c, 3*nat_c, nkc/)), &
                                               gk_loc, 3*nat_c, 'bands_wann_spin_peratom.dat', irank, mpi_comm)
        END IF
      CASE ('orbital')
        nbnd_c = SIZE(l0_loc, 1); nat_c = SIZE(l0_loc, 4); nkc = SIZE(l0_loc, 5)
        ! total (site-summed) orbital moment
        IF (this%op_total(iop) == 1) &
          CALL wannierlib_interpolate_operator(this, cell, kpts, eig, u_matrix, amn_local, &
                                               SUM(l0_loc, DIM=4), gk_loc, 3, 'bands_wann_orbmom.dat', irank, mpi_comm)
        ! per-atom (site-resolved): flatten (comp,atom) -> 3*nat components in one file
        IF (this%op_peratom(iop) == 1) &
          CALL wannierlib_interpolate_operator(this, cell, kpts, eig, u_matrix, amn_local, &
                                               RESHAPE(l0_loc, (/nbnd_c, nbnd_c, 3*nat_c, nkc/)), &
                                               gk_loc, 3*nat_c, 'bands_wann_orbmom_peratom.dat', irank, mpi_comm)
      CASE ('soc')
        CALL wannierlib_interpolate_operator(this, cell, kpts, eig, u_matrix, amn_local, &
                                             soc0_loc, gk_loc, 1, 'bands_wann_soc.dat', irank, mpi_comm)
      CASE ('velocity')
        ! Wannier Berry connection A^(W)_alpha(R): built distributed from the local overlaps and
        ! reduced (collective, all ranks); the centre check (rank 0) calibrates conj/sign. Built
        ! once and reused across output domains. amn_local = u_opt.
        IF (.NOT.ALLOCATED(aw_r)) THEN
          CALL wannierlib_build_berry_aw_r(this, cell, kpts, mmn, gk_loc, amn_local, u_matrix, mpi_comm, &
                                           aw_r, aw_irvec, aw_ndegen, aw_nrpts)
          IF (irank == 0) CALL wannierlib_check_berry_centres(this, aw_r, aw_irvec, aw_nrpts)
        END IF
        CALL wannierlib_interpolate_velocity(this, cell, kpts, eig, u_matrix, amn_local, &
                                             aw_r, aw_irvec, aw_ndegen, aw_nrpts, irank)
      CASE ('spinCurrent')
        ! operator part distributed like the generic driver: local Bloch slice + gk_loc + reduce
        CALL wannierlib_interpolate_current(this, cell, kpts, eig, u_matrix, amn_local, &
                                            s0_loc, gk_loc, 'bands_wann_spincurrent.dat', irank, mpi_comm)
      CASE ('orbitalCurrent')
        CALL wannierlib_interpolate_current(this, cell, kpts, eig, u_matrix, amn_local, &
                                            SUM(l0_loc, DIM=4), gk_loc, 'bands_wann_orbcurrent.dat', irank, mpi_comm)
      CASE ('eigenstates')
        ! Wannier-Hamiltonian eigenvectors C(k') (the H-gauge rotation U^(H)), written as a matrix
        CALL wannierlib_interpolate_eigenstates(this, cell, kpts, eig, u_matrix, amn_local, irank)
      CASE DEFAULT
        IF (irank == 0) WRITE(oUnit,'(a)') 'wannierlib: operator "'//TRIM(this%op_name(iop))//&
                                           '" not yet implemented -> skipped'
      END SELECT
    END DO
    ! rename this domain's outputs (plane/grid -> _plane/_grid; path/legacy: no suffix)
    IF (irank==0 .AND. LEN_TRIM(TRIM(dsuf(idom))//TRIM(ssfx))>0) &
      CALL wannierlib_rename_domain_outputs(this, TRIM(dsuf(idom))//TRIM(ssfx))
    END DO   ! idom

    ! restore the user's original kpts_interpol if we overwrote it
    IF (irank==0 .AND. (this%l_dom_plane .OR. this%l_dom_grid)) THEN
      INQUIRE(file='.kpts_interpol_userbak', exist=lex)
      IF (lex) CALL wl_shell('mv -f .kpts_interpol_userbak kpts_interpol')
    END IF

    call timestop('run_w90')
#endif
  END SUBROUTINE run_w90

  ! Write kpts_interpol for a generated output domain (plane/grid). For an explicit
  ! <path file="..">, copy that file to kpts_interpol; legacy/default path uses the
  ! existing kpts_interpol as-is. Rank-0 file I/O only.
  SUBROUTINE wannierlib_write_domain_kpts(this, kind)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: kind
    INTEGER :: i, j, k, iu, np
    REAL :: t1, t2, kf(3)
    LOGICAL :: lex2
    np = 0
    SELECT CASE (TRIM(kind))
    CASE ('path', 'legacy')
      ! restore the user's original kpts_interpol if a generated (plane/grid) domain overwrote it
      INQUIRE(file='.kpts_interpol_userbak', exist=lex2)
      IF (lex2) CALL wl_shell('cp -f .kpts_interpol_userbak kpts_interpol')
      ! <path listName=".."> : write the (optionally subdivided) named list as kpts_interpol.
      ! Takes precedence over @file; path_kpts was filled from kpts.xml in read_xml.
      IF (TRIM(kind) == 'path' .AND. this%path_np > 0) THEN
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') this%path_np
        DO i = 1, this%path_np
          WRITE(iu,'(3(f18.12,1x))') this%path_kpts(:, i)
        END DO
        CLOSE(iu)
        RETURN
      END IF
      ! explicit <path file="..">: use that file as the k-list (default 'kpts_interpol' -> no-op)
      IF (TRIM(kind) == 'path' .AND. TRIM(this%path_file) /= 'kpts_interpol') THEN
        INQUIRE(file=TRIM(this%path_file), exist=lex2)
        IF (.NOT. lex2) CALL juDFT_error('wannierlib: <path>/@file "'//TRIM(this%path_file)//'" not found', &
                                         calledby='wannierlib_write_domain_kpts')
        CALL wl_shell('cp -f '//TRIM(this%path_file)//' kpts_interpol')
      END IF
      RETURN
    CASE ('plane')
      IF (this%plane_np > 0) THEN     ! <plane listName=".."/>: use the named list as-is
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') this%plane_np
        DO i = 1, this%plane_np
          WRITE(iu,'(3(f18.12,1x))') this%plane_kpts(:, i)
        END DO
        CLOSE(iu)
      ELSE
        np = this%plane_n1 * this%plane_n2
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') np
        DO i = 0, this%plane_n1 - 1
          t1 = REAL(i) / REAL(MAX(1, this%plane_n1 - 1))
          DO j = 0, this%plane_n2 - 1
            t2 = REAL(j) / REAL(MAX(1, this%plane_n2 - 1))
            kf = this%plane_origin + t1 * this%plane_v1 + t2 * this%plane_v2
            WRITE(iu,'(3(f18.12,1x))') kf
          END DO
        END DO
        CLOSE(iu)
      END IF
    CASE ('grid')
      IF (this%grid_np > 0) THEN      ! <grid listName=".."/>: use the named list as-is
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') this%grid_np
        DO i = 1, this%grid_np
          WRITE(iu,'(3(f18.12,1x))') this%grid_kpts(:, i)
        END DO
        CLOSE(iu)
      ELSE
        np = this%grid_mesh(1) * this%grid_mesh(2) * this%grid_mesh(3)
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') np
        DO i = 0, this%grid_mesh(1) - 1
        DO j = 0, this%grid_mesh(2) - 1
        DO k = 0, this%grid_mesh(3) - 1
          kf(1) = (REAL(i) + this%grid_shift(1)) / REAL(this%grid_mesh(1))
          kf(2) = (REAL(j) + this%grid_shift(2)) / REAL(this%grid_mesh(2))
          kf(3) = (REAL(k) + this%grid_shift(3)) / REAL(this%grid_mesh(3))
          WRITE(iu,'(3(f18.12,1x))') kf
        END DO
        END DO
        END DO
        CLOSE(iu)
      END IF
    END SELECT
  END SUBROUTINE wannierlib_write_domain_kpts

  ! Rename this domain's operator output files bands_wann_<x>.dat -> bands_wann_<x><suffix>.dat
  SUBROUTINE wannierlib_rename_domain_outputs(this, suffix)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    INTEGER :: iop
    DO iop = 1, this%n_ops
      SELECT CASE (TRIM(this%op_name(iop)))
      CASE ('hamiltonian')
        CALL ren('bands_wann_interpol', suffix)
        CALL ren('bands_wann_interpol_ev', suffix)
      CASE ('spin')
        IF (this%op_total(iop) == 1)   CALL ren('bands_wann_spin', suffix)
        IF (this%op_peratom(iop) == 1) CALL ren('bands_wann_spin_peratom', suffix)
      CASE ('orbital')
        IF (this%op_total(iop) == 1)   CALL ren('bands_wann_orbmom', suffix)
        IF (this%op_peratom(iop) == 1) CALL ren('bands_wann_orbmom_peratom', suffix)
      CASE ('soc')
        CALL ren('bands_wann_soc', suffix)
      CASE ('velocity')
        CALL ren('bands_wann_velocity', suffix)
        CALL ren('bands_wann_berrycurv', suffix)
      CASE ('spinCurrent')
        CALL ren('bands_wann_spincurrent', suffix)
      CASE ('orbitalCurrent')
        CALL ren('bands_wann_orbcurrent', suffix)
      CASE ('eigenstates')
        CALL ren('bands_wann_eigenstates', suffix)
      END SELECT
    END DO
  CONTAINS
    SUBROUTINE ren(base, suf)
      CHARACTER(LEN=*), INTENT(IN) :: base, suf
      LOGICAL :: lexr
      INQUIRE(file=TRIM(base)//'.dat', exist=lexr)
      IF (lexr) CALL wl_shell('mv -f '//TRIM(base)//'.dat '//TRIM(base)//TRIM(suf)//'.dat')
    END SUBROUTINE ren
  END SUBROUTINE wannierlib_rename_domain_outputs

  ! Run a shell command (synchronous) and abort with a clear message if it fails,
  ! so a failed cp/mv in the domain file-shuffling never passes silently.
  SUBROUTINE wl_shell(cmd)
    CHARACTER(LEN=*), INTENT(IN) :: cmd
    INTEGER :: cs, es
    cs = 0; es = 0
    CALL EXECUTE_COMMAND_LINE(cmd, wait=.TRUE., cmdstat=cs, exitstat=es)
    IF (cs /= 0 .OR. es /= 0) CALL juDFT_error('wannierlib: shell command failed: '//TRIM(cmd), &
                                               calledby='wl_shell')
  END SUBROUTINE wl_shell

  ! Build the Wannier-gauge Berry connection in real space A^(W)_alpha(R), distributed:
  ! each rank forms A^(W)_alpha(k) = i sum_b w_b b_alpha (M^(W,b)(k) - delta) on its OWN k-slice
  ! (M^(W,b)(k) = V(k)^dagger M^(0,b)(k) V(k_b), V = u_opt.u_matrix, k_b = nnlist(b,k)) from the
  ! local overlaps mmn_loc (global k = gk_loc), then reduces coarse -> R with the distributed
  ! FT-reduce (collective). The full-mesh overlaps are never gathered. A^(W)(R) is exactly the
  ! position operator A(R) = <0n|r|Rm> and the R->k' interpolant used by the velocity/Berry curvature.
  ! kmesh (wb/bk/nnlist) from wannierlib_w90main. Collective over mpi_comm; result valid on all ranks.
  SUBROUTINE wannierlib_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, mpi_comm, &
                                         aw_r, irvec, ndegen, nrpts)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)    ! (nb,nb,nntot,nk_loc) this rank's overlap slice
    INTEGER, INTENT(IN) :: gk_loc(:)              ! (nk_loc) global k index of each slice entry
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (nb,nw,nk) full mesh
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (nw,nw,nk) full mesh
    INTEGER, INTENT(IN) :: mpi_comm
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: aw_r(:, :, :, :)   ! (nw,nw,nrpts,3) reduced Berry connection / A(R)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
#ifdef CPP_WANNLIB_API
    INTEGER :: nb, nw, nk, nnt, k, kb, nn, a, i, kl, nkl
    REAL :: wb, b(3)
    COMPLEX, ALLOCATABLE :: Vk(:, :), Vkb(:, :), Mw(:, :), tmp(:, :), aw_loc(:, :, :, :), a1(:, :, :)
    nb = this%num_bands; nw = this%num_wann; nk = kpts%nkptf
    nnt = wannierlib_w90main%kmesh_info%nntot
    nkl = SIZE(gk_loc)
    ALLOCATE(aw_loc(nw, nw, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(Vk(nb, nw), Vkb(nb, nw), Mw(nw, nw), tmp(nb, nw))
    DO kl = 1, nkl
      k = gk_loc(kl)
      Vk = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
      DO nn = 1, nnt
        kb = wannierlib_w90main%kmesh_info%nnlist(k, nn)   ! w90 shape: nnlist(num_kpts, nntot)
        IF (kb < 1 .OR. kb > nk) CALL juDFT_error('wannierlib: nnlist neighbour index out of range', &
                                                  calledby='wannierlib_build_berry_aw_r')
        wb = wannierlib_w90main%kmesh_info%wb(nn)
        b  = wannierlib_w90main%kmesh_info%bk(:, nn, k)
        Vkb = MATMUL(u_opt(:, :, kb), u_matrix(:, :, kb))
        tmp = MATMUL(mmn_loc(:, :, nn, kl), Vkb)            ! (nb x nw)
        Mw  = MATMUL(CONJG(TRANSPOSE(Vk)), tmp)             ! (nw x nw) = M^(W,b)(k)
        DO i = 1, nw
          Mw(i, i) = Mw(i, i) - CMPLX(1.0, 0.0)            ! subtract delta on the diagonal
        END DO
        DO a = 1, 3
          aw_loc(:, :, a, kl) = aw_loc(:, :, a, kl) + CMPLX(0.0, wb*b(a)) * Mw
        END DO
      END DO
    END DO
    DEALLOCATE(Vk, Vkb, Mw, tmp)
    DO a = 1, 3
      CALL wannierlib_ft_to_real_reduce(cell, kpts, aw_loc(:, :, a, :), gk_loc, mpi_comm, a1, irvec, ndegen, nrpts)
      IF (a == 1) ALLOCATE(aw_r(nw, nw, nrpts, 3))
      aw_r(:, :, :, a) = a1; DEALLOCATE(a1)
    END DO
    DEALLOCATE(aw_loc)
#else
    nrpts = 0
    ALLOCATE(aw_r(1, 1, 1, 3), irvec(3, 1), ndegen(1))
#endif
  END SUBROUTINE wannierlib_build_berry_aw_r

  ! Validation: the diagonal of A^(W)_alpha at R=0 is the Wannier centre <r_alpha>_n
  ! (Marzari-Vanderbilt). Since A^(W)_alpha(R=0) = (1/Nk) sum_k A^(W)_alpha(k), it is exactly
  ! the R=0 entry of the reduced aw_r. Compare to w90_get_centres to calibrate the conj/sign
  ! convention of the overlaps. Writes berry_centre_check.dat (rank 0).
  SUBROUTINE wannierlib_check_berry_centres(this, aw_r, irvec, nrpts)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)     ! (nw,nw,nrpts,3)
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts
#ifdef CPP_WANNLIB_API
    INTEGER :: nw, a, n, iu, irpt, irpt0
    COMPLEX :: aR0
    REAL, ALLOCATABLE :: wc(:, :)
    nw = this%num_wann
    irpt0 = 0
    DO irpt = 1, nrpts
      IF (ALL(irvec(:, irpt) == 0)) THEN; irpt0 = irpt; EXIT; END IF
    END DO
    IF (irpt0 == 0) RETURN   ! no R=0 vector (should not happen) -> skip the check
    ALLOCATE(wc(3, nw))
    CALL w90_get_centres(wannierlib_w90main, wc)
    OPEN(newunit=iu, file='berry_centre_check.dat', status='replace')
    WRITE(iu,'(a)') '# n  alpha    Re[A_nn(R=0)]        -Re[A_nn(R=0)]       w90_centre(Bohr)'
    DO n = 1, nw
      DO a = 1, 3
        aR0 = aw_r(n, n, irpt0, a)
        WRITE(iu,'(2i4,3(2x,f18.10))') n, a, REAL(aR0), -REAL(aR0), wc(a, n)
      END DO
    END DO
    CLOSE(iu)
    DEALLOCATE(wc)
    WRITE(oUnit,'(a)') 'wannierlib: wrote berry_centre_check.dat (A^(W) R=0 diag vs w90 centres)'
#endif
  END SUBROUTINE wannierlib_check_berry_centres

  ! Build the Wannier-gauge Hamiltonian ham_k = U^dag diag(eigval2) U (same as m_wannierlib_interpolate).
  SUBROUTINE wannierlib_build_hamk(this, eig, u_matrix, u_opt, ham_k)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    REAL,    INTENT(IN) :: eig(:, :)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :), u_opt(:, :, :)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: ham_k(:, :, :)
    INTEGER :: nw, nb, nk, k, i, j, m, counter
    LOGICAL :: have_dis
    REAL, ALLOCATABLE :: eigval2(:, :), eigval_opt(:)
    nw = this%num_wann; nb = this%num_bands; nk = SIZE(u_matrix, 3); have_dis = (nb > nw)
    ALLOCATE(eigval2(nw, nk), source=0.0)
    IF (have_dis) THEN
      ALLOCATE(eigval_opt(nb))
      DO k = 1, nk
        counter = 0; eigval_opt = 0.0
        DO j = 1, nb
          IF (eig(j, k) >= this%dis_win_min .AND. eig(j, k) <= this%dis_win_max) THEN
            counter = counter + 1; eigval_opt(counter) = eig(j, k)
          END IF
        END DO
        DO m = 1, nw
          DO i = 1, counter
            eigval2(m, k) = eigval2(m, k) + eigval_opt(i) * ABS(u_opt(i, m, k))**2
          END DO
        END DO
      END DO
      DEALLOCATE(eigval_opt)
    ELSE
      eigval2(1:nw, :) = eig(1:nw, :)
    END IF
    ALLOCATE(ham_k(nw, nw, nk), source=CMPLX(0.0, 0.0))
    DO k = 1, nk
      DO j = 1, nw
        DO i = 1, nw
          DO m = 1, nw
            ham_k(i, j, k) = ham_k(i, j, k) + eigval2(m, k) * CONJG(u_matrix(m, i, k)) * u_matrix(m, j, k)
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(eigval2)
  END SUBROUTINE wannierlib_build_hamk

  ! Write H(R) in Wannier90 seedname_hr.dat format (energies in eV). Rank-0 only.
  SUBROUTINE wannierlib_write_hr(this, cell, kpts, ham_k, wfpref)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: ham_k(:, :, :)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref   ! seedname prefix 'WF1'/'WF2' (collinear jspins=2 channel); default 'WF1'
    COMPLEX, ALLOCATABLE :: hr(:, :, :)
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    INTEGER :: nrpts, nw, irpt, i, j, iu, c
    CHARACTER(LEN=64) :: fn
    CALL wannierlib_ft_to_real(cell, ham_k, kpts, hr, irvec, ndegen, nrpts)
    nw = this%num_wann
    fn = 'WF1_hr'
    IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_hr'
    OPEN(newunit=iu, file=TRIM(fn)//'.dat', status='replace')
    WRITE(iu,'(a)') ' written by FLEUR wannierlib : H(R) in eV, W90 hr format'
    WRITE(iu,'(i12)') nw
    WRITE(iu,'(i12)') nrpts
    c = 0
    DO irpt = 1, nrpts
      WRITE(iu,'(i5)',advance='no') ndegen(irpt); c = c + 1
      IF (MOD(c,15) == 0) WRITE(iu,'(a)') ''
    END DO
    IF (MOD(c,15) /= 0) WRITE(iu,'(a)') ''
    DO irpt = 1, nrpts
      DO j = 1, nw
        DO i = 1, nw
          WRITE(iu,'(5i5,2f12.6)') irvec(:,irpt), i, j, &
              hartree_to_ev_const*REAL(hr(i,j,irpt)), hartree_to_ev_const*AIMAG(hr(i,j,irpt))
        END DO
      END DO
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (H(R), eV)'
    DEALLOCATE(hr, irvec, ndegen)
  END SUBROUTINE wannierlib_write_hr

  ! Write A(R) = <0n|r_alpha|Rm> in Wannier90 seedname_r.dat format (positions in Angstrom). Rank-0.
  ! aw_r is the already-reduced Berry connection A^(W)(R) (= A(R)), so no Fourier transform here.
  SUBROUTINE wannierlib_write_ar(this, aw_r, irvec, nrpts, wfpref)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)          ! (nw,nw,nrpts,3) reduced A(R)
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref   ! seedname prefix 'WF1'/'WF2' (collinear jspins=2 channel); default 'WF1'
    INTEGER :: nw, irpt, i, j, a, iu
    CHARACTER(LEN=64) :: fn
    REAL, PARAMETER :: bohr2ang = 0.5291772109
    nw = this%num_wann
    fn = 'WF1_r'
    IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_r'
    OPEN(newunit=iu, file=TRIM(fn)//'.dat', status='replace')
    WRITE(iu,'(a)') ' written by FLEUR wannierlib : A(R)=<0n|r|Rm> in Ang, W90 r format'
    WRITE(iu,'(i12)') nw
    WRITE(iu,'(i12)') nrpts
    DO irpt = 1, nrpts
      DO j = 1, nw
        DO i = 1, nw
          WRITE(iu,'(5i5,6f12.6)') irvec(:,irpt), i, j, &
            (bohr2ang*REAL(aw_r(i,j,irpt,a)), bohr2ang*AIMAG(aw_r(i,j,irpt,a)), a=1,3)
        END DO
      END DO
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (A(R), Ang)'
  END SUBROUTINE wannierlib_write_ar

  ! ---------------------------------------------------------------------------
  ! <operators_r>: real-space Wannier operator matrices O(R) (Fourier step 3),
  ! in Wannier90/FLEUR standalone format for external transport post-processing.
  ! No band interpolation. Rank 0 writes; spin/orbital/spin_orbit use a distributed FT-reduce
  ! over MPI ranks. Files: WF1_hr.dat (H, eV), WF1_r.dat (position, Ang), rspauli.1 (spin),
  ! anglmomrs.1 (orbital), rssocmat.1 (SOC), wig_vectors.
  ! ---------------------------------------------------------------------------
  SUBROUTINE wannierlib_write_operators_r(this, cell, kpts, eig, u_matrix, u_opt, &
                                          s0_loc, l0_loc, soc4_loc, distk, mpi_comm, mmn_loc, irank, wf_channel, l_collinear, l0col_loc)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)         ! (nw,nw,nk) MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)            ! (nb,nw,nk) disentangled (amn_local)
    COMPLEX, INTENT(IN) :: s0_loc(:, :, :, :), l0_loc(:, :, :, :, :), soc4_loc(:, :, :, :) ! per-rank coarse slices
    COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)       ! (nb,nb,nntot,nk_loc) this rank's overlap slice (position/Berry)
    INTEGER, INTENT(IN) :: distk(:), mpi_comm, irank
    INTEGER, INTENT(IN) :: wf_channel               ! collinear jspins=2 spin channel (1/2); 1 in the spinor case
    LOGICAL, INTENT(IN) :: l_collinear              ! collinear jspins=2 (spin channels separable, per-channel WFn files)
    COMPLEX, INTENT(IN), OPTIONAL :: l0col_loc(:, :, :, :)  ! (nb,nb,3,nk_loc) collinear per-channel Bloch orbital L
    INTEGER :: iop, nb, ik, j, nkl, aw_nrpts
    LOGICAL :: l_wig_done
    CHARACTER(LEN=8) :: wfpref                       ! 'WF1'/'WF2' seedname prefix for H(R)/position
    INTEGER, ALLOCATABLE :: gk_loc(:), aw_irvec(:, :), aw_ndegen(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), aw_r(:, :, :, :), o0l(:, :, :, :), vloc(:, :, :)
    IF (.NOT. this%l_operators_r .OR. this%n_op_r < 1) RETURN   ! all ranks (reduce is collective)
    CALL timestart('wannierlib_write_operators_r')
    nb = this%num_bands
    nkl = COUNT(distk == irank); ALLOCATE(gk_loc(nkl)); j = 0
    DO ik = 1, SIZE(distk)
      IF (distk(ik) == irank) THEN; j = j + 1; gk_loc(j) = ik; END IF
    END DO
    ! precompute the gauge V(gk)=u_opt.u_matrix for this rank's k-slice ONCE (reused by spin/orbital/soc)
    ALLOCATE(vloc(this%num_bands, this%num_wann, MAX(1, nkl)), source=CMPLX(0.0,0.0))
    DO j = 1, nkl
      vloc(:, :, j) = MATMUL(u_opt(:, :, gk_loc(j)), u_matrix(:, :, gk_loc(j)))
    END DO
    ! collinear jspins=2: per-channel seedname WF1/WF2 for H(R)/position; spinor case keeps WF1.
    WRITE(wfpref, '(a,i0)') 'WF', wf_channel
    IF (irank == 0) THEN; l_wig_done = .FALSE.; CALL wannierlib_write_wig_once(cell, kpts, l_wig_done); END IF
    DO iop = 1, this%n_op_r
      SELECT CASE (TRIM(this%op_r_name(iop)))
      CASE ('hamiltonian')   ! cheap, no coarse array -> rank 0 serial
        IF (irank == 0) THEN
          CALL wannierlib_build_hamk(this, eig, u_matrix, u_opt, ham_k)
          CALL wannierlib_write_hr(this, cell, kpts, ham_k, TRIM(wfpref))
          DEALLOCATE(ham_k)
        END IF
      CASE ('position')      ! Berry connection A^(W)(R)=A(R): distributed reduce over the local overlaps
        CALL wannierlib_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, mpi_comm, &
                                         aw_r, aw_irvec, aw_ndegen, aw_nrpts)
        IF (irank == 0) CALL wannierlib_write_ar(this, aw_r, aw_irvec, aw_nrpts, TRIM(wfpref))
        IF (ALLOCATED(aw_r)) DEALLOCATE(aw_r, aw_irvec, aw_ndegen)
      CASE ('spin')          ! distributed reduce over the coarse spin slice
        ! collinear: the spin operator is combined (2N) across both channels, so it is written once
        ! after both wannierisations by wannierlib_rspauli_collinear (main), not per channel here.
        IF (l_collinear) THEN
          CONTINUE
        ELSE
          CALL wannierlib_op_rs_distributed(this, cell, kpts, vloc, s0_loc, gk_loc, 3, mpi_comm, irank, .FALSE., 'rspauli.1')
        END IF
      CASE ('orbital')
        IF (l_collinear) THEN
          ! per-channel Bloch L built in main's mmn loop (single spin) -> reduce -> anglmomrs.{1,2}
          IF (PRESENT(l0col_loc)) THEN
            CALL wannierlib_op_rs_distributed(this, cell, kpts, vloc, l0col_loc, gk_loc, 3, mpi_comm, irank, .FALSE., &
                                              'anglmomrs.'//ACHAR(48+wf_channel))
          ELSE IF (irank == 0) THEN
            WRITE(oUnit,'(a)') 'wannierlib operators_r: collinear orbital slice missing -> skipped'
          END IF
        ELSE
          ALLOCATE(o0l(nb, nb, 3, SIZE(l0_loc, 5))); o0l = SUM(l0_loc, DIM=4)
          CALL wannierlib_op_rs_distributed(this, cell, kpts, vloc, o0l, gk_loc, 3, mpi_comm, irank, .FALSE., 'anglmomrs.1')
          DEALLOCATE(o0l)
        END IF
      CASE ('spin_orbit')
        IF (l_collinear) THEN
          IF (irank == 0) WRITE(oUnit,'(a)') 'wannierlib operators_r: spin_orbit has no collinear (no-SOC) meaning -> skipped'
        ELSE
          CALL wannierlib_op_rs_distributed(this, cell, kpts, vloc, soc4_loc, gk_loc, 4, mpi_comm, irank, .TRUE., 'rssocmat.1')
        END IF
      CASE DEFAULT
        IF (irank == 0) WRITE(oUnit,'(a)') 'wannierlib operators_r: unknown operator "'//TRIM(this%op_r_name(iop))//'" -> skipped'
      END SELECT
    END DO
    DEALLOCATE(gk_loc, vloc)
    CALL timestop('wannierlib_write_operators_r')
  END SUBROUTINE wannierlib_write_operators_r

  ! Distributed real-space export of a coarse operator (spin/orbital/soc): each rank builds
  ! ow_loc = V(gk)^dagger o0_loc V(gk) for its k-slice, FT-reduces to O(R), rank 0 writes.
  ! is_soc=.TRUE. -> rssocmat.1 format (R i j jj ii); else rspauli/anglmomrs (R i j comp).
  SUBROUTINE wannierlib_op_rs_distributed(this, cell, kpts, vloc, o0_loc, gk_loc, ncomp, mpi_comm, irank, is_soc, fname)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: vloc(:, :, :)                      ! (nb,nw,>=nk_loc) precomputed gauge V(gk)
    COMPLEX, INTENT(IN) :: o0_loc(:, :, :, :)                 ! (nb,nb,ncomp,>=nk_loc) local Bloch op
    INTEGER, INTENT(IN) :: gk_loc(:), ncomp, mpi_comm, irank
    LOGICAL, INTENT(IN) :: is_soc
    CHARACTER(LEN=*), INTENT(IN) :: fname
    COMPLEX, ALLOCATABLE :: ow_loc(:, :, :, :), or_(:, :, :, :), o1(:, :, :), tmp(:, :)
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    INTEGER :: nb, nw, nkl, kl, a, nrpts, irpt, i, j, kk, ii, jj, c, iu
    nb = SIZE(o0_loc, 1); nw = this%num_wann; nkl = SIZE(gk_loc)
    ALLOCATE(ow_loc(nw, nw, ncomp, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(tmp(nb, nw))
    DO kl = 1, nkl
      DO a = 1, ncomp
        tmp = MATMUL(o0_loc(:, :, a, kl), vloc(:, :, kl))
        ow_loc(:, :, a, kl) = MATMUL(CONJG(TRANSPOSE(vloc(:, :, kl))), tmp)
      END DO
    END DO
    DEALLOCATE(tmp)
    DO a = 1, ncomp
      CALL wannierlib_ft_to_real_reduce(cell, kpts, ow_loc(:, :, a, :), gk_loc, mpi_comm, o1, irvec, ndegen, nrpts)
      IF (a == 1) ALLOCATE(or_(nw, nw, nrpts, ncomp))
      or_(:, :, :, a) = o1; DEALLOCATE(o1)
    END DO
    DEALLOCATE(ow_loc)
    IF (irank == 0) THEN
      OPEN(newunit=iu, file=TRIM(fname), status='replace')
      IF (is_soc) THEN
        DO irpt = 1, nrpts; DO i = 1, nw; DO j = 1, nw; DO ii = 1, 2; DO jj = 1, 2
          c = (ii-1)*2 + jj
          WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
            irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, jj, ii, REAL(or_(i,j,irpt,c)), AIMAG(or_(i,j,irpt,c))
        END DO; END DO; END DO; END DO; END DO
      ELSE
        DO irpt = 1, nrpts; DO j = 1, nw; DO i = 1, nw; DO kk = 1, ncomp
          WRITE(iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
            irvec(1,irpt), irvec(2,irpt), irvec(3,irpt), i, j, kk, REAL(or_(i,j,irpt,kk)), AIMAG(or_(i,j,irpt,kk))
        END DO; END DO; END DO; END DO
      END IF
      CLOSE(iu)
      WRITE(oUnit, '(a,i0,a)') 'wannierlib: wrote '//TRIM(fname)//' (', nrpts, ' R-vectors, distributed FT)'
    END IF
    IF (ALLOCATED(or_)) DEALLOCATE(or_)
    IF (ALLOCATED(irvec)) DEALLOCATE(irvec, ndegen)
  END SUBROUTINE wannierlib_op_rs_distributed

  ! wig_vectors (once): idx R1 R2 R3 ndegen (list-directed; = orbitrans rpts).
  SUBROUTINE wannierlib_write_wig_once(cell, kpts, done)
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    LOGICAL, INTENT(INOUT) :: done
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    INTEGER :: nrpts, ii, iu
    IF (done) RETURN
    CALL wannierlib_ws_vectors(cell, kpts%nkpt3, irvec, ndegen, nrpts)
    OPEN(newunit=iu, file='wig_vectors', status='replace')
    DO ii = 1, nrpts
      WRITE(iu, *) ii, irvec(1, ii), irvec(2, ii), irvec(3, ii), ndegen(ii)
    END DO
    CLOSE(iu)
    WRITE(oUnit, '(a,i0,a)') 'wannierlib: wrote wig_vectors (', nrpts, ' R-vectors)'
    DEALLOCATE(irvec, ndegen)
    done = .TRUE.
  END SUBROUTINE wannierlib_write_wig_once

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
