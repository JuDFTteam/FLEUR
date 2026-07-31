!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Driver for the matrix elements of physical operators in the Wannier representation.
!>
!>  This module owns the whole operator pipeline, so that the wannierization side
!>  (m_wannierlib_main, m_wannierlib_w90_adapter) only has to hand over its results and
!>  never has to grow another operator argument:
!>
!>    t_melem_coarse%init   allocate the per-rank Bloch-basis operator slices and decide
!>                          which of them the requested <operator>/<operators_r> lists need
!>    t_melem_coarse%calc   ONE shared coarse k-pass (read_eig + abc per k) feeding every
!>                          Bloch-basis provider (spin / orbital / SOC). Needs only the
!>                          ab-initio eigenstates, so it runs BEFORE the wannierization.
!>    melem_run             everything that needs the Wannier gauge V = u_opt.u_matrix:
!>                          the real-space O(R) export, and the band interpolation of every
!>                          requested operator over every requested output domain.
!>    melem_rspauli_collinear   the collinear (jspins=2) combined 2N spin operator, which can
!>                          only be assembled once BOTH spin channels have been wannierised.
!>
!>  The Bloch-basis providers themselves live in m_melem_spin / m_melem_orbmom /
!>  m_melem_socmat, the k <-> R Fourier core in m_melem_ft, the real-space writers in
!>  m_melem_operators_r and the per-operator band interpolation in the m_melem_interpolate_*
!>  modules. Nothing here knows about the Wannier90 library: the only piece of W90 state the
!>  operators need (the b-shell weights of the position/Berry operator, and the reference
!>  Wannier centres) arrives as a plain t_melem_bmesh.
MODULE m_melem_driver
   USE m_juDFT
   USE m_constants, ONLY: ImagUnit, oUnit
   USE m_types_atoms
   USE m_types_cell
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_sym
   USE m_types_enpara
   USE m_types_mpi
   USE m_types_potden
   USE m_types_stars
   USE m_types_usdus
   USE m_types_mat
   USE m_types_radfun
   USE m_types_abc
   USE m_types_wannierlib
   USE m_types_melem_bmesh
   USE m_melem_spin, ONLY: melem_spin_peratom, melem_spin_bloch, melem_spin_mt_block
   USE m_melem_orbmom, ONLY: melem_orbmom_bloch, melem_orbmom_bloch_collinear
   USE m_melem_socmat, ONLY: melem_socmat_bloch
   USE m_melem_ft, ONLY: melem_ft_to_real_reduce
   USE m_melem_domains, ONLY: melem_write_domain_kpts, melem_rename_domain_outputs, melem_shell
   USE m_melem_operators_r, ONLY: melem_write_operators_r, melem_build_berry_aw_r, melem_check_berry_centres
   USE m_melem_interpolate_ham, ONLY: melem_interpolate_ham
   USE m_melem_interpolate_op, ONLY: melem_interpolate_operator
   USE m_melem_interpolate_velocity, ONLY: melem_interpolate_velocity
   USE m_melem_interpolate_current, ONLY: melem_interpolate_current
   USE m_melem_interpolate_eigenstates, ONLY: melem_interpolate_eigenstates
   USE m_wannierlib_get_z, ONLY: wannierlib_get_z
   USE m_wannierlib_mmkb_int, ONLY: wannierlib_mmnkb_int
   IMPLICIT NONE
   PRIVATE

   !> The Bloch-basis coarse-mesh operator matrices O^0(k), per-rank slices only: the full
   !> coarse mesh is never materialized. Slice entries are stored in ascending global-k order,
   !> matching the gk_loc convention of the distributed FT-reduce.
   TYPE :: t_melem_coarse
      COMPLEX, ALLOCATABLE :: s0(:, :, :, :)      !< (nb,nb,3,nk_loc)      spin
      COMPLEX, ALLOCATABLE :: l0(:, :, :, :, :)   !< (nb,nb,3,nat,nk_loc)  orbital L per atom
      COMPLEX, ALLOCATABLE :: soc0(:, :, :, :)    !< (nb,nb,1,nk_loc)      SOC
      COMPLEX, ALLOCATABLE :: soc4(:, :, :, :)    !< (nb,nb,4,nk_loc)      2x2 SOC spinor blocks
      COMPLEX, ALLOCATABLE :: s0pa(:, :, :, :, :) !< (nb,nb,3,nat,nk_loc)  per-atom MT spin
      !> collinear jspins=2 only: per-channel Bloch orbital L, filled from the wannierization's
      !> own mmn k-loop (see add_collinear_orbital) because the spinor coarse pass does not run.
      COMPLEX, ALLOCATABLE :: l0col(:, :, :, :)   !< (nb,nb,3,nk_loc)
      !> collinear jspins=2 only: the gauge V of each spin channel, needed by the combined 2N
      !> spin operator. Filled by melem_run, consumed by melem_rspauli_collinear.
      COMPLEX, ALLOCATABLE :: v_ch(:, :, :, :)    !< (nb,nw,nkptf,2)
      LOGICAL :: l_col_orb = .FALSE.   !< collinear jspins=2 AND orbital requested in <operators_r>
      LOGICAL :: l_col_spin = .FALSE.  !< collinear jspins=2 AND spin requested in <operators_r>
      !> .TRUE. only when the spinor coarse slices were really allocated (an operator is requested
      !> AND we have spinor wavefunctions). Gates %calc, so it can never write into the stubs.
      LOGICAL :: l_active = .FALSE.
   CONTAINS
      PROCEDURE :: init => melem_coarse_init
      PROCEDURE :: calc => melem_coarse_calc
      PROCEDURE :: alloc_collinear_orbital => melem_coarse_alloc_collinear_orbital
      PROCEDURE :: free => melem_coarse_free
   END TYPE t_melem_coarse

   PUBLIC :: t_melem_coarse, melem_run, melem_rspauli_collinear, melem_orbmom_bloch_collinear

CONTAINS

   !> Allocate the per-rank Bloch slices for the operators the input actually asks for, and
   !> derive the collinear special cases. Stub (1,1,1,1) allocations elsewhere keep the
   !> downstream assumed-shape dummies valid without branching at every call site.
   SUBROUTINE melem_coarse_init(this, wann, atoms, input, kpts, fmpi, distk, l_spinors)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wann
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)
      LOGICAL, INTENT(IN) :: l_spinors   !< noco%l_noco .OR. noco%l_soc

      INTEGER :: nkc_loc, iop

      ! Operator Bloch matrices on the coarse mesh: the k-loop is DISTRIBUTED over ranks (each
      ! its distk slice -> parallel get_z I/O) into per-rank local arrays. Every consumer works
      ! on those slices plus a distributed FT-reduce, so the full mesh is never assembled.
      this%l_active = (wann%l_spin .OR. wann%l_orbmom .OR. wann%l_socop) .AND. l_spinors
      IF (this%l_active) THEN
         nkc_loc = MAX(1, COUNT(distk == fmpi%irank))
         ALLOCATE (this%s0(wann%num_bands, wann%num_bands, 3, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%l0(wann%num_bands, wann%num_bands, 3, atoms%nat, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc0(wann%num_bands, wann%num_bands, 1, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc4(wann%num_bands, wann%num_bands, 4, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%s0pa(wann%num_bands, wann%num_bands, 3, atoms%nat, nkc_loc), source=cmplx(0.0, 0.0))
      ELSE
         ALLOCATE (this%s0(1, 1, 1, 1)); ALLOCATE (this%l0(1, 1, 1, 1, 1)); ALLOCATE (this%soc4(1, 1, 1, 1))
         ALLOCATE (this%soc0(1, 1, 1, 1)); ALLOCATE (this%s0pa(1, 1, 1, 1, 1))
      END IF

      ! collinear jspins=2 (no SOC/noco): the coarse spin/orbital slices above are spinor-only
      ! (stubs), so the per-channel orbital operators_r builds its own Bloch L in the mmn k-loop
      ! (reusing that loop's abc), and the combined 2N spin operator (rspauli.1) is assembled
      ! after both channels wannierise from their gauges v_ch + the cross-spin overlap.
      this%l_col_orb = .FALSE.; this%l_col_spin = .FALSE.
      IF (input%jspins == 2 .AND. .NOT. l_spinors .AND. wann%l_operators_r) THEN
         DO iop = 1, wann%n_op_r
            IF (TRIM(wann%op_r_name(iop)) == 'orbital') this%l_col_orb = .TRUE.
            IF (TRIM(wann%op_r_name(iop)) == 'spin') this%l_col_spin = .TRUE.
         END DO
      END IF
      IF (this%l_col_spin) THEN
         ALLOCATE (this%v_ch(wann%num_bands, wann%num_wann, kpts%nkptf, 2), source=cmplx(0.0, 0.0))
      ELSE
         ALLOCATE (this%v_ch(1, 1, 1, 1))
      END IF
      ! l0col is (re)sized per spin channel by alloc_collinear_orbital; allocate a stub now so it
      ! is never passed unallocated to the OPTIONAL dummy of melem_write_operators_r.
      IF (.NOT. ALLOCATED(this%l0col)) ALLOCATE (this%l0col(1, 1, 1, 1))
   END SUBROUTINE melem_coarse_init

   !> Per-spin-channel (re)allocation of the collinear orbital slice. Called once per channel
   !> by the wannierization before its mmn k-loop starts feeding add_collinear_orbital.
   SUBROUTINE melem_coarse_alloc_collinear_orbital(this, wann, nk_local)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wann
      INTEGER, INTENT(IN) :: nk_local

      IF (ALLOCATED(this%l0col)) DEALLOCATE (this%l0col)
      IF (this%l_col_orb) THEN
         ALLOCATE (this%l0col(wann%num_bands, wann%num_bands, 3, MAX(1, nk_local)), source=cmplx(0.0, 0.0))
      ELSE
         ALLOCATE (this%l0col(1, 1, 1, 1))
      END IF
   END SUBROUTINE melem_coarse_alloc_collinear_orbital

   !> Build the requested Bloch-basis operator matrices on this rank's coarse-mesh k-slice:
   !> spin S0(k), per-atom spin, orbital L0(k) and/or SOC, sharing the (get_z + calc_abc) work
   !> in ONE k-pass. Needs only the ab-initio spinor eigenstates (no gauge), so it runs before
   !> the wannierization.
   SUBROUTINE melem_coarse_calc(this, wann, atoms, input, sym, cell, noco, nococonv, kpts, &
                                stars, usdus, radfun, enpara, fmpi, vtot, eig_id, l_real_wann, distk)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wann
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_usdus), INTENT(INOUT) :: usdus     ! INOUT: spnorb (SOC) fills it
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_potden), INTENT(IN) :: vtot
      INTEGER, INTENT(IN) :: eig_id
      LOGICAL, INTENT(IN) :: l_real_wann
      INTEGER, INTENT(IN) :: distk(:)   ! rank owner of each global k (distributes the loop)

      TYPE(t_abc), ALLOCATABLE :: abc_s(:, :)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      TYPE(t_mat) :: zc(2)   ! the two spinor components when get_z does not stack them
      INTEGER :: ikpt, itype, isp, il, jspin_rad

      IF (.NOT. this%l_active) RETURN   ! nothing requested, or no spinor wavefunctions -> slices are stubs

      ALLOCATE (abc_s(2, atoms%ntype))
      il = 0
      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE   ! this rank only computes its own k-slice
         il = il + 1                            ! local (ascending global-k) index for the reduce
         ! Load this k-point's eigenvector(s).
         !   l_noco=T          : get_z returns the whole 2N spinor from record 1.
         !   l_soc=T, l_noco=F : get_z returns only N rows and the two spinor components live
         !     in records 1 and 2, so read both and stack them into the 2N layout the
         !     interstitial expects. Reading record 1 alone leaves the spin-down half unread:
         !     the muffin-tin counts the up block twice and the interstitial addresses a down
         !     block that is not there (non-magnetic Pt then sums to <sigma_z> = +N/2, not 0).
         IF (noco%l_noco) THEN
            CALL wannierlib_get_z(wann, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                  ikpt, 1, l_real_wann, lapw, zMat)
         ELSE
            CALL wannierlib_get_z(wann, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                  ikpt, 1, l_real_wann, lapw, zc(1))
            CALL wannierlib_get_z(wann, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                  ikpt, 2, l_real_wann, lapw, zc(2))
            CALL melem_stack_spinor(zc(1), zc(2), zMat)
         END IF
         DO isp = 1, 2
            ! The index handed to calc_abc must belong to the zMat it is given and must never
            ! exceed input%jspins -- calc_abc rejects that since the develop SOC refactor, and
            ! the radial arrays it indexes are only that wide anyway.
            jspin_rad = MERGE(1, isp, input%jspins == 1)
            DO itype = 1, atoms%ntype
               CALL abc_s(isp, itype)%init(input, atoms, wann%num_bands, itype)
               IF (noco%l_noco) THEN
                  CALL abc_s(isp, itype)%calc_abc(input, atoms, sym, cell, lapw, wann%num_bands, usdus, &
                                                  noco, nococonv, jspin_rad, itype, zMat)
               ELSE
                  CALL abc_s(isp, itype)%calc_abc(input, atoms, sym, cell, lapw, wann%num_bands, usdus, &
                                                  noco, nococonv, jspin_rad, itype, zc(isp))
               END IF
            END DO
         END DO
         IF (wann%l_spin) CALL melem_spin_peratom(atoms, abc_s, radfun, nococonv, this%s0pa(:, :, :, :, il))
         IF (wann%l_spin) CALL melem_spin_bloch(atoms, abc_s, radfun, nococonv, stars, lapw, zMat, &
                                                wann%num_bands, ikpt, this%s0(:, :, :, il), ikpt <= 3)
         IF (wann%l_orbmom) CALL melem_orbmom_bloch(atoms, abc_s, radfun, this%l0(:, :, :, :, il))
         IF (wann%l_socop) CALL melem_socmat_bloch(atoms, noco, nococonv, input, fmpi, enpara, vtot, &
                                                  usdus, abc_s, wann%num_bands, this%soc0(:, :, :, il), &
                                                  this%soc4(:, :, :, il))
      END DO

      DEALLOCATE (abc_s)
   END SUBROUTINE melem_coarse_calc

   !> ===================== the driver =====================
   !>
   !> Everything that needs the Wannier gauge. Called by the wannierization once per
   !> wannierised spin channel, with the gauge factors it just obtained:
   !>
   !>   (1) fill this channel's gauge V = u_opt.u_matrix (collinear combined spin needs it later)
   !>   (2) real-space export of the requested <operators_r> -> O(R) files
   !>   (3) for each requested output domain (<path>/<plane>/<grid>, or the legacy single pass):
   !>       interpolate every requested <operator> onto that domain's k-set and write its
   !>       bands_wann_*.dat, then suffix the files so domains/channels do not overwrite.
   !>
   !> Collective over fmpi%mpi_comm (the FT-reduces are); the file IO is rank 0.
   SUBROUTINE melem_run(wann, cell, kpts, eig, u_matrix, u_opt, coarse, mmn, bmesh, distk, fmpi, &
                        wf_channel, spin_suffix)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wann
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      REAL, INTENT(IN) :: eig(:, :)                    !< (nb,nk)
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)         !< (nw,nw,nk) MLWF gauge
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)            !< (nb,nw,nk) disentangled
      TYPE(t_melem_coarse), INTENT(INOUT) :: coarse    !< INOUT: v_ch of this channel is filled
      COMPLEX, INTENT(IN) :: mmn(:, :, :, :)           !< (nb,nb,nntot,nk_loc) this rank's overlap slice
      TYPE(t_melem_bmesh), INTENT(IN) :: bmesh         !< b-shell weights (position/velocity operators)
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN), OPTIONAL :: wf_channel            !< collinear spin channel (1/2); default 1
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spin_suffix  !< '_spin1'/'_spin2' collinear; empty otherwise

      INTEGER :: iop, k, nbnd_c, nat_c, nkc, wf_ch, irank, mpi_comm
      INTEGER :: idom, ndom, nkl_c, jkl, aw_nrpts
      INTEGER, ALLOCATABLE :: gk_loc(:), aw_irvec(:, :), aw_ndegen(:)
      COMPLEX, ALLOCATABLE :: aw_r(:, :, :, :)         ! (nw,nw,nrpts,3) Wannier Berry connection A^(W)(R)
      CHARACTER(LEN=8) :: dkind(3), dsuf(3)
      CHARACTER(LEN=16) :: ssfx
      LOGICAL :: lex, l_collinear

      IF (.NOT. wann%l_wannierize) RETURN
      irank = fmpi%irank; mpi_comm = fmpi%mpi_comm
      ssfx = ''
      IF (PRESENT(spin_suffix)) ssfx = TRIM(spin_suffix)
      wf_ch = 1
      IF (PRESENT(wf_channel)) wf_ch = wf_channel
      l_collinear = (LEN_TRIM(ssfx) > 0)   ! collinear jspins=2 -> per-channel operators_r (WF1/WF2)

      CALL timestart('melem_run')

      ! (1) the full gauge V(k) = u_opt(k) . u_matrix(k) (disentangled + MLWF) of this channel,
      ! needed by the collinear combined spin operator, which rotates the cross-spin overlap
      ! with BOTH channels' V once they are both available.
      IF (coarse%l_col_spin) THEN
         DO k = 1, SIZE(u_opt, 3)
            coarse%v_ch(:, :, k, wf_ch) = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
         END DO
      END IF

      ! global-k indices owned by this rank, ascending order -> matches the per-rank coarse
      ! slices in `coarse` (built in the same distk order); used by the distributed reduces.
      nkl_c = COUNT(distk == irank); ALLOCATE (gk_loc(nkl_c)); jkl = 0
      DO iop = 1, SIZE(distk)
         IF (distk(iop) == irank) THEN; jkl = jkl + 1; gk_loc(jkl) = iop; END IF
      END DO

      ! ---- opt-in output domains (<path>/<plane>/<grid>); none declared -> legacy single pass ----
      ! order matters: generated domains (plane/grid) overwrite kpts_interpol and are renamed;
      ! the unsuffixed path/legacy domain runs LAST so its base-named output is not clobbered
      ! and it restores the user's original kpts_interpol before interpolating.
      ndom = 0
      IF (wann%l_dom_plane) THEN; ndom = ndom + 1; dkind(ndom) = 'plane'; dsuf(ndom) = '_plane'; END IF
      IF (wann%l_dom_grid) THEN; ndom = ndom + 1; dkind(ndom) = 'grid'; dsuf(ndom) = '_grid'; END IF
      IF (wann%l_dom_path) THEN; ndom = ndom + 1; dkind(ndom) = 'path'; dsuf(ndom) = ''; END IF
      IF (ndom == 0) THEN; ndom = 1; dkind(1) = 'legacy'; dsuf(1) = ''; END IF
      ! back up a user-provided kpts_interpol that a generated (plane/grid) domain would overwrite
      IF (irank == 0 .AND. (wann%l_dom_plane .OR. wann%l_dom_grid)) THEN
         INQUIRE (file='kpts_interpol', exist=lex)
         IF (lex) CALL melem_shell('cp -f kpts_interpol .kpts_interpol_userbak')
      END IF

      ! (2) real-space operator export (Fourier step 3, standalone format) -- once, before interpolation
      CALL melem_write_operators_r(wann, cell, kpts, eig, u_matrix, u_opt, &
                                   coarse%s0, coarse%l0, coarse%soc4, bmesh, distk, mpi_comm, mmn, &
                                   irank, wf_ch, l_collinear, coarse%l0col)

      ! (3) Wannier-gauge interpolation: dispatch by looping over the requested operator list.
      ! Each operator supplies its own per-rank Bloch slice on the coarse mesh (coarse%s0/l0/soc0);
      ! the remaining steps are the shared generic driver m_melem_interpolate_op.
      DO idom = 1, ndom
         IF (irank == 0) CALL melem_write_domain_kpts(wann, TRIM(dkind(idom)))

         DO iop = 1, wann%n_ops
            SELECT CASE (TRIM(wann%op_name(iop)))
            CASE ('hamiltonian')
               CALL melem_interpolate_ham(wann, cell, kpts, eig, u_matrix, u_opt, irank)
            CASE ('spin')
               ! total spin (MT-sum + interstitial): via the generic operator driver (3 comps)
               IF (wann%op_total(iop) == 1) &
                  CALL melem_interpolate_operator(wann, cell, kpts, eig, u_matrix, u_opt, &
                                                  coarse%s0, gk_loc, 3, 'bands_wann_spin.dat', irank, mpi_comm)
               ! per-atom (site-resolved) muffin-tin spin moment: 3*nat components in one file
               IF (wann%op_peratom(iop) == 1) THEN
                  nbnd_c = SIZE(coarse%s0pa, 1); nat_c = SIZE(coarse%s0pa, 4); nkc = SIZE(coarse%s0pa, 5)
                  CALL melem_interpolate_operator(wann, cell, kpts, eig, u_matrix, u_opt, &
                                                  RESHAPE(coarse%s0pa, (/nbnd_c, nbnd_c, 3*nat_c, nkc/)), &
                                                  gk_loc, 3*nat_c, 'bands_wann_spin_peratom.dat', irank, mpi_comm)
               END IF
            CASE ('orbital')
               nbnd_c = SIZE(coarse%l0, 1); nat_c = SIZE(coarse%l0, 4); nkc = SIZE(coarse%l0, 5)
               ! total (site-summed) orbital moment
               IF (wann%op_total(iop) == 1) &
                  CALL melem_interpolate_operator(wann, cell, kpts, eig, u_matrix, u_opt, &
                                                  SUM(coarse%l0, DIM=4), gk_loc, 3, 'bands_wann_orbmom.dat', irank, mpi_comm)
               ! per-atom (site-resolved): flatten (comp,atom) -> 3*nat components in one file
               IF (wann%op_peratom(iop) == 1) &
                  CALL melem_interpolate_operator(wann, cell, kpts, eig, u_matrix, u_opt, &
                                                  RESHAPE(coarse%l0, (/nbnd_c, nbnd_c, 3*nat_c, nkc/)), &
                                                  gk_loc, 3*nat_c, 'bands_wann_orbmom_peratom.dat', irank, mpi_comm)
            CASE ('soc')
               CALL melem_interpolate_operator(wann, cell, kpts, eig, u_matrix, u_opt, &
                                               coarse%soc0, gk_loc, 1, 'bands_wann_soc.dat', irank, mpi_comm)
            CASE ('velocity')
               ! Wannier Berry connection A^(W)_alpha(R): built distributed from the local overlaps
               ! and reduced (collective, all ranks); the centre check (rank 0) calibrates conj/sign.
               ! Built once and reused across output domains.
               IF (.NOT. ALLOCATED(aw_r)) THEN
                  CALL melem_build_berry_aw_r(wann, cell, kpts, mmn, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                              aw_r, aw_irvec, aw_ndegen, aw_nrpts)
                  IF (irank == 0) CALL melem_check_berry_centres(wann, aw_r, aw_irvec, aw_nrpts, bmesh)
               END IF
               CALL melem_interpolate_velocity(wann, cell, kpts, eig, u_matrix, u_opt, &
                                               aw_r, aw_irvec, aw_ndegen, aw_nrpts, irank)
            CASE ('spinCurrent')
               ! operator part distributed like the generic driver: local Bloch slice + gk_loc + reduce
               CALL melem_interpolate_current(wann, cell, kpts, eig, u_matrix, u_opt, &
                                              coarse%s0, gk_loc, 'bands_wann_spincurrent.dat', irank, mpi_comm)
            CASE ('orbitalCurrent')
               CALL melem_interpolate_current(wann, cell, kpts, eig, u_matrix, u_opt, &
                                              SUM(coarse%l0, DIM=4), gk_loc, 'bands_wann_orbcurrent.dat', irank, mpi_comm)
            CASE ('eigenstates')
               ! Wannier-Hamiltonian eigenvectors C(k') (the H-gauge rotation U^(H)), as a matrix
               CALL melem_interpolate_eigenstates(wann, cell, kpts, eig, u_matrix, u_opt, irank)
            CASE DEFAULT
               IF (irank == 0) WRITE (oUnit, '(a)') 'wannierlib: operator "'//TRIM(wann%op_name(iop))// &
                  '" not yet implemented -> skipped'
            END SELECT
         END DO
         ! rename this domain's outputs (plane/grid -> _plane/_grid; path/legacy: no suffix)
         IF (irank == 0 .AND. LEN_TRIM(TRIM(dsuf(idom))//TRIM(ssfx)) > 0) &
            CALL melem_rename_domain_outputs(wann, TRIM(dsuf(idom))//TRIM(ssfx))
      END DO   ! idom

      ! restore the user's original kpts_interpol if we overwrote it
      IF (irank == 0 .AND. (wann%l_dom_plane .OR. wann%l_dom_grid)) THEN
         INQUIRE (file='.kpts_interpol_userbak', exist=lex)
         IF (lex) CALL melem_shell('mv -f .kpts_interpol_userbak kpts_interpol')
      END IF

      IF (ALLOCATED(aw_r)) DEALLOCATE (aw_r, aw_irvec, aw_ndegen)
      DEALLOCATE (gk_loc)
      CALL timestop('melem_run')
   END SUBROUTINE melem_run

   !> Collinear (jspins=2, no SOC/noco) combined 2N spin operator -> rspauli.1.
   !> Over the joint {up-WF, down-WF} space: sigma_z is block-diagonal (+/- identity, WFs are
   !> orthonormal per channel); sigma_x/sigma_y come from the cross-spin overlap
   !> o_ud(k) = <psi_k^up | psi_k^dn> (interstitial via wannierlib_mmnkb_int with b=0, muffin-tin
   !> via melem_spin_mt_block fed with BOTH channels' abc), rotated to the WF gauge with each
   !> channel's V = u_opt.u_matrix (coarse%v_ch), then FT-reduced (distributed). Rank 0 writes
   !> rspauli.1 (2N, standalone R i j comp Re Im format). Reuses our own abc machinery -- no
   !> updown.mmn0 needed.
   SUBROUTINE melem_rspauli_collinear(wann, atoms, input, sym, cell, noco, nococonv, kpts, &
                                      stars, usdus, radfun, eig_id, l_real_wann, distk, fmpi, coarse)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wann
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_usdus), INTENT(IN) :: usdus
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER, INTENT(IN) :: eig_id
      LOGICAL, INTENT(IN) :: l_real_wann
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_melem_coarse), INTENT(IN) :: coarse   !< carries v_ch (num_bands,num_wann,nkptf,2)

      TYPE(t_lapw) :: lapw_u, lapw_d
      TYPE(t_mat)  :: zMat_u, zMat_d
      TYPE(t_abc), ALLOCATABLE :: abc_both(:, :)   ! (ntype,2): 1=up, 2=dn
      INTEGER :: nb, nw, n2, nkl, kl, gk, itype, ikpt, iu, irpt, i, j, kk, nrpts, gb(3)
      INTEGER, ALLOCATABLE :: gk_loc(:), irvec(:, :), ndegen(:)
      COMPLEX, ALLOCATABLE :: o_uu(:, :), o_dd(:, :), o_ud(:, :), o_du(:, :), Xk(:, :), tmp(:, :)
      COMPLEX, ALLOCATABLE :: sig_loc(:, :, :, :), s1(:, :, :), sr(:, :, :, :)

      IF (.NOT. coarse%l_col_spin) RETURN

      nb = wann%num_bands; nw = wann%num_wann; n2 = 2*nw; gb = 0
      nkl = COUNT(distk == fmpi%irank); ALLOCATE (gk_loc(nkl)); j = 0
      DO ikpt = 1, SIZE(distk)
         IF (distk(ikpt) == fmpi%irank) THEN; j = j + 1; gk_loc(j) = ikpt; END IF
      END DO
      ALLOCATE (abc_both(atoms%ntype, 2))
      ALLOCATE (o_uu(nb, nb), o_dd(nb, nb), o_ud(nb, nb), o_du(nb, nb), Xk(nw, nw), tmp(nb, nw))
      ALLOCATE (sig_loc(n2, n2, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))

      DO kl = 1, nkl
         gk = gk_loc(kl)
         ! up + down eigenvectors + local coefficients at this k (same basis, different eigenvectors)
         CALL wannierlib_get_z(wann, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, gk, 1, l_real_wann, lapw_u, zMat_u)
         DO itype = 1, atoms%ntype
            CALL abc_both(itype, 1)%init(input, atoms, nb, itype)
            CALL abc_both(itype, 1)%calc_abc(input, atoms, sym, cell, lapw_u, nb, usdus, noco, nococonv, 1, itype, zMat_u)
         END DO
         CALL wannierlib_get_z(wann, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, gk, 2, l_real_wann, lapw_d, zMat_d)
         DO itype = 1, atoms%ntype
            CALL abc_both(itype, 2)%init(input, atoms, nb, itype)
            CALL abc_both(itype, 2)%calc_abc(input, atoms, sym, cell, lapw_d, nb, usdus, noco, nococonv, MERGE(1, 2, input%jspins == 1), itype, zMat_d)
         END DO
         ! cross-spin overlap o_ud = <up|dn>: interstitial (b=0) + muffin-tin (both channels' abc)
         o_uu = CMPLX(0.0, 0.0); o_dd = CMPLX(0.0, 0.0); o_ud = CMPLX(0.0, 0.0); o_du = CMPLX(0.0, 0.0)
         CALL wannierlib_mmnkb_int(stars, lapw_u, lapw_d, 1, 1, zMat_u, zMat_d, gb, o_ud, ioff=0, ioff_b=0)
         CALL melem_spin_mt_block(atoms, abc_both, radfun, o_uu, o_dd, o_ud, o_du)
         ! rotate to the WF gauge: X = V_up^dagger o_ud V_dn
         tmp = MATMUL(o_ud, coarse%v_ch(:, :, gk, 2))
         Xk = MATMUL(CONJG(TRANSPOSE(coarse%v_ch(:, :, gk, 1))), tmp)
         ! assemble the 2N Pauli in the WF gauge (sigma_z block-diagonal +/- I by orthonormality)
         DO i = 1, nw
            sig_loc(i, i, 3, kl) = CMPLX(1.0, 0.0)
            sig_loc(nw + i, nw + i, 3, kl) = CMPLX(-1.0, 0.0)
         END DO
         sig_loc(1:nw, nw + 1:n2, 1, kl) = Xk
         sig_loc(nw + 1:n2, 1:nw, 1, kl) = CONJG(TRANSPOSE(Xk))
         sig_loc(1:nw, nw + 1:n2, 2, kl) = -ImagUnit*Xk
         sig_loc(nw + 1:n2, 1:nw, 2, kl) = ImagUnit*CONJG(TRANSPOSE(Xk))
      END DO
      DEALLOCATE (o_uu, o_dd, o_ud, o_du, Xk, tmp, abc_both)

      ! FT-reduce each of the 3 components (collective), rank 0 writes rspauli.1
      DO kk = 1, 3
         CALL melem_ft_to_real_reduce(cell, kpts, sig_loc(:, :, kk, :), gk_loc, fmpi%mpi_comm, s1, irvec, ndegen, nrpts)
         IF (kk == 1) ALLOCATE (sr(n2, n2, nrpts, 3))
         sr(:, :, :, kk) = s1; DEALLOCATE (s1)
      END DO
      IF (fmpi%irank == 0) THEN
         OPEN (newunit=iu, file='rspauli.1', status='replace')
         DO irpt = 1, nrpts
            DO j = 1, n2
               DO i = 1, n2
                  DO kk = 1, 3
                     WRITE (iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
                        irvec(1, irpt), irvec(2, irpt), irvec(3, irpt), i, j, kk, REAL(sr(i, j, irpt, kk)), AIMAG(sr(i, j, irpt, kk))
                  END DO
               END DO
            END DO
         END DO
         CLOSE (iu)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote rspauli.1 (combined 2N collinear spin, ', nrpts, ' R-vectors, distributed FT)'
      END IF
      DEALLOCATE (sig_loc, gk_loc)
      IF (ALLOCATED(sr)) DEALLOCATE (sr)
      IF (ALLOCATED(irvec)) DEALLOCATE (irvec, ndegen)
   END SUBROUTINE melem_rspauli_collinear

   SUBROUTINE melem_coarse_free(this)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this

      IF (ALLOCATED(this%s0)) DEALLOCATE (this%s0)
      IF (ALLOCATED(this%l0)) DEALLOCATE (this%l0)
      IF (ALLOCATED(this%soc0)) DEALLOCATE (this%soc0)
      IF (ALLOCATED(this%soc4)) DEALLOCATE (this%soc4)
      IF (ALLOCATED(this%s0pa)) DEALLOCATE (this%s0pa)
      IF (ALLOCATED(this%l0col)) DEALLOCATE (this%l0col)
      IF (ALLOCATED(this%v_ch)) DEALLOCATE (this%v_ch)
      this%l_col_orb = .FALSE.; this%l_col_spin = .FALSE.
   END SUBROUTINE melem_coarse_free

   !> Stack the two spin components of a collinear-basis spinor (l_soc=T, l_noco=F, one eig
   !> record each) into the single 2N matrix the l_noco path produces, so every consumer --
   !> notably the interstitial in melem_spin, which addresses the down block at row offset
   !> nv(1)+nlotot -- sees one layout. The shared LAPW basis is checked, not assumed.
   SUBROUTINE melem_stack_spinor(zup, zdn, zspinor)
      TYPE(t_mat), INTENT(IN)  :: zup, zdn
      TYPE(t_mat), INTENT(OUT) :: zspinor
      INTEGER :: n, nb
      IF (zup%matsize1 /= zdn%matsize1 .OR. zup%matsize2 /= zdn%matsize2) &
         CALL juDFT_error('melem_stack_spinor: the two spin records differ in shape', &
                          calledby='melem_stack_spinor')
      IF (zup%l_real .OR. zdn%l_real) &
         CALL juDFT_error('melem_stack_spinor: a SOC spinor cannot be stored in a real matrix', &
                          calledby='melem_stack_spinor')
      n = zup%matsize1; nb = zup%matsize2
      CALL zspinor%init(.FALSE., 2*n, nb)
      zspinor%data_c(1:n,     1:nb) = zup%data_c(1:n, 1:nb)
      zspinor%data_c(n+1:2*n, 1:nb) = zdn%data_c(1:n, 1:nb)
   END SUBROUTINE melem_stack_spinor

END MODULE m_melem_driver
