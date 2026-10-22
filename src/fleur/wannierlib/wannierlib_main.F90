!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!>  Driver of the library-mode wannierization: build the projections A_mn(k) and the
!>  neighbour overlaps M_mn(k,b) that Wannier90 needs, drive the Wannier90 library, and
!>  report the resulting centres and spreads.
!>
!>  Matrix elements of physical operators (spin, orbital moment, SOC, position, velocity,
!>  currents) are NOT built here: they live in src/fleur/matrixelements/.
!>  This module only hands that layer what it cannot obtain on its
!>  own -- the k-point distribution, the per-k abc coefficients of the collinear channels,
!>  the Wannier gauge, and the b-mesh -- so that adding an operator never touches this file.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_main
   USE m_juDFT
   USE m_matrix_element_factory, ONLY: matrix_element_factory_reset, &
                                       matrix_element_release_anchor, matrix_element_radial
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_wannierlib_request, ONLY: wannierlib_build_request
   USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
   USE m_types_wgauge_domains, ONLY: t_wgauge_domains
   USE m_wannierlib_build_amn_mmn, ONLY: wannierlib_build_amn_mmn, wannierlib_reduce_amn, &
                                      wannierlib_gather_mmn
   USE m_wannierlib_plot, ONLY: wannierlib_plot_wf
  USE m_wannierlib_cf, ONLY: wannierlib_uiu, wannierlib_uhu
   USE m_wannierlib_w90_adapter
   USE m_melem_coarse, ONLY: t_melem_coarse
   USE m_wgauge_run, ONLY: wgauge_run
   USE m_wannierlib_export_sym, ONLY: wannierlib_export_sym
   USE m_wannierlib_sym_kstar, ONLY: t_wl_kstar, wl_build_kstar
   USE m_wannierlib_dband, ONLY: wl_build_dband
   USE m_wannierlib_export_orbitrans, ONLY: wannierlib_export_orbitrans
   USE m_wannierlib_export_gauge, ONLY: wannierlib_export_gauge, wannierlib_store_gauge
USE m_wannierlib_mmnkb, ONLY: wannierlib_kdiff
USE m_wannierlib_band_window, ONLY: wannierlib_default_windows, wannierlib_create_eig
   USE m_wannierlib_export_bloch, ONLY: wannierlib_write_eig, wannierlib_write_s0, &
                                     wannierlib_write_mmn, wannierlib_write_amn
   USE m_wgauge_spin_collinear, ONLY: wgauge_rspauli_collinear, wgauge_anglmom_collinear, wgauge_soc_collinear
   USE m_types_wgauge_bmesh, ONLY: t_wgauge_bmesh
   USE m_constants, ONLY: oUnit, hartree_to_ev_const
   USE m_types_atoms
   USE m_types_cell
   USE m_types_vacuum
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_sym
   USE m_types_enpara
   USE m_types_mpi
   USE m_types_potden
   USE m_types_mat
   USE m_types_radfun
   USE m_types_abc
   USE m_types_wannierlib
   USE m_types_misc
   USE m_types_stars
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_main
CONTAINS

   SUBROUTINE wannierlib_main(this, atoms, cell, input, kpts, sym, noco, nococonv, stars, enpara, fmpi, vtot, results, eig_id, vacuum)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_potden), INTENT(IN) :: vtot
      TYPE(t_results), INTENT(INOUT) :: results
      INTEGER, INTENT(IN) :: eig_id

      INTEGER :: nntot_w90, jspin
      COMPLEX, ALLOCATABLE :: amn(:, :, :)
      COMPLEX, ALLOCATABLE :: mmn(:, :, :, :)
      COMPLEX, ALLOCATABLE :: mmn_full(:, :, :, :)
      !> A copy of the settings, because the energy windows may be left out of the input and
      !> are filled in from the bands here. It cannot be done in place: this comes from the
      !> fleurinput container, which reaches the caller as INTENT(IN).
      TYPE(t_wannierlib_wannierize) :: wl
      REAL, ALLOCATABLE :: kdiff(:, :)
      INTEGER, ALLOCATABLE :: nnkp(:, :), gkpb(:, :, :)
      INTEGER, ALLOCATABLE :: distk(:)
      real, allocatable :: eig(:, :)
      COMPLEX, ALLOCATABLE :: u_matrix(:, :, :), u_opt(:, :, :)   ! Wannier gauge factors from w90
      !> The gauge of each spin channel, kept across the channel loop: the combined 2N spin
      !> operator needs BOTH, while u_matrix and u_opt are released at the end of every
      !> iteration. It lives here rather than in the coarse object, which holds the Bloch
      !> matrices from before any gauge exists.
      COMPLEX, ALLOCATABLE :: v_ch(:, :, :, :)   ! (num_bands, num_wann, nkptf, 2)
      TYPE(t_melem_coarse) :: melem   ! the operator (matrix-element) side
      TYPE(t_wgauge_bmesh) :: bmesh    ! b-shell weights handed to the operator side
      TYPE(t_radfun), POINTER :: radfun(:)  ! likewise; the factory owns them
      COMPLEX, ALLOCATABLE :: f0_loc(:, :, :, :, :)  ! (nw,nw,3,3,nk_loc) geometric tensor
      COMPLEX, ALLOCATABLE :: c0_loc(:, :, :, :, :)  ! (nw,nw,3,3,nk_loc) the same one with H inside
      LOGICAL :: l_wannierlib_spinors
      TYPE(t_melem_request) :: request
      TYPE(t_wgauge_manifold) :: manifold
      TYPE(t_wgauge_domains) :: domains
      TYPE(t_wl_kstar) :: kstar        ! k-point stars, for the .dmn and the symmetrizer
      COMPLEX, ALLOCATABLE :: dband(:, :, :, :)   ! the group on the bands, when it can be had
      LOGICAL :: l_dband, l_realph, l_tr
      CHARACTER(LEN=7) :: amn_file
      CHARACTER(LEN=3) :: spin12(2)
      CHARACTER(LEN=6) :: spin_sfx

      IF (.NOT. this%l_wannierize) RETURN


      l_wannierlib_spinors = noco%l_noco .OR. noco%l_soc
      ! A.1: input%l_real stays TRUE with inversion even under SOC (n_denmat=0), and
      ! reading a complex spinor into a real buffer kills the imaginary part of the MMN.
      spin12 = (/'WF1', 'WF2'/)

      !> The factory cache is module state, and what it holds belongs to whoever ran last.
      !> Nothing guarantees they cleared it: secvar_soc does so on entry and on exit, this
      !> did so only on exit, and the two stay out of each other's way only because they
      !> are sequential phases of fleur.F90. Clearing it on entry costs one rebuild of the
      !> radial functions and removes the assumption.
      CALL matrix_element_factory_reset()

      !> The radial functions come from the factory, which keeps them for the operators
      !> anyway. Generating a second set here produced the same numbers twice.
      CALL matrix_element_radial(atoms, input, enpara, fmpi, vtot, radfun)

      ! Settled before anything else: it distributes the coarse-operator k-loop too.
      CALL wannierlib_distribute_k(kpts, fmpi, distk)

      ! Hand the k-distribution to the matrix-element layer and let it build whatever Bloch-basis
      ! operators the input asks for, in one shared coarse k-pass. Needs only the ab-initio
      ! eigenstates (no gauge), so it runs before the wannierization. All of it is a no-op when
      ! no operator is requested.
      CALL wannierlib_build_request(this, input, l_wannierlib_spinors, request)

      wl = this
      CALL wannierlib_default_windows(wl, results, kpts, input, l_wannierlib_spinors)

      CALL manifold%init(wl%num_bands, wl%num_wann, wl%dis_win_min, wl%dis_win_max, &
                         wl%min_band, wl%max_band)
      CALL domains%init(this%n_domains, this%dom_kset, this%dom_suffix)

      CALL melem%init(request, manifold, atoms, input, kpts, fmpi, distk, l_wannierlib_spinors)
      CALL melem%calc(request, manifold, atoms, input, sym, cell, noco, nococonv, kpts, &
                      stars, enpara, fmpi, vtot, eig_id, distk, vacuum=vacuum)

      !> wl, not this: wannierlib_main takes the input object as INTENT(IN), so the windows
      !> the derivation fills in live only in the copy. Handing w90 the original passes it
      !> whatever the input stated and zeros for the rest, which is a window of [0,0] and a
      !> disentanglement over an empty subspace. manifold%init and run_w90 already take wl;
      !> this call was the one left behind.
      CALL init_w90(wl, atoms, cell, kpts, fmpi, l_wannierlib_spinors, nntot_w90, nnkp, gkpb, distk)
      CALL wannierlib_kdiff(kpts%nkptf, nntot_w90, kpts%bkf, nnkp, gkpb, kdiff)
      !> The neighbour topology is settled here, before anything is wannierised, and it is
      !> the same for both spin channels. The shell weights are added later, per channel.
      CALL bmesh%set_neighbours(nntot_w90, nnkp, gkpb, kdiff)

      DO jspin = 1, MERGE(2, 1, input%jspins == 2 .AND. (.NOT. l_wannierlib_spinors))

         CALL wannierlib_build_amn_mmn(this, manifold, bmesh, atoms, cell, input, kpts, sym, &
                                       noco, nococonv, stars, enpara, fmpi, vtot, eig_id, &
                                       radfun, distk, kdiff, nntot_w90, jspin, &
                                       l_wannierlib_spinors, amn, mmn, vacuum)

         ! amn was filled only on each rank's distk slice (zeros elsewhere) -> sum to the full set
         CALL wannierlib_reduce_amn(fmpi, amn)

         mmn = conjg(mmn)

         call wannierlib_create_eig(this, results, kpts, MERGE(1, jspin, l_wannierlib_spinors), eig)

         !> The three files Wannier90 was handed, so a standalone wannier90.x can repeat the
         !> run: .amn, .mmn and .eig. They are written together and after create_eig because
         !> without the eigenvalues the other two cannot be used -- w90 stops at
         !> "No .eig file found. Needed for disentanglement".
         !> amn is complete on every rank after the reduce above; mmn is one k-slice per rank
         !> and has to be gathered first. Only rank 0 writes, and neither writer collects
         !> anything itself: what they are handed is already whole.
         !> The flag alone decides. This used to write on one rank whatever was asked, which
         !> left the .mmn of a large mesh -- hundreds of megabytes of text, the reason the
         !> flag is opt-in at all -- on the disk of every serial run that never wanted it.
         !> wannierberri asks for them too: WannierBerri does no DFT, so its wannierisation
         !> consumes these three and takes the symmetry table as an argument beside them.
         !> What wannier90="T" writes is a subset of what wannierberri="T" writes.
         IF (this%export%w90 .OR. this%export%symmetrizer) THEN
            CALL wannierlib_gather_mmn(fmpi, distk, kpts%nkptf, mmn, mmn_full)
            IF (fmpi%irank == 0) THEN
               amn_file = spin12(jspin)//'.amn'
               CALL wannierlib_write_amn(amn, amn_file, "Testing amn")

               CALL wannierlib_write_mmn(this, mmn_full, kpts, nnkp, gkpb, jspin)
               CALL wannierlib_write_eig(eig, spin12(jspin)//'.eig')
            END IF
            DEALLOCATE (mmn_full)
         END IF

         !> The k-point stars under the space group. The stars and the matrices are both
         !> indexed by a star representative and an operation, so the partition is built ONCE
         !> here rather than inside the writer: a table derived twice could disagree, and a
         !> consumer handed a star table that does not match the matrices beside it
         !> symmetrises at the wrong k -- it converges, the centres come out symmetric, and
         !> the result is wrong.
         IF (this%export%symmetrizer) THEN
            !> Time reversal goes in when it is a symmetry of the Hamiltonian, which is to
            !> say without magnetism. It is not a choice the user should have to make: a
            !> magnet with it would be wrong and the smaller irreducible set would hide it,
            !> and a non-magnet without it is merely a bigger set. Measured on bcc Fe with
            !> spin-orbit: T by itself misses by 7.2e-03 while the operations that preserve
            !> the moment close at 2.9e-12.
            l_tr = (input%jspins == 1) .AND. (.NOT. noco%l_noco)
            CALL wl_build_kstar(kpts, sym, fmpi, kstar, l_tr)

            !> d_band asks for the state at one k expressed at another, and two things
            !> stand in the way of that rotation. One is the basis itself -- a
            !> centrosymmetric crystal is solved in REAL eigenvectors, and a non-symmorphic
            !> operation multiplies each coefficient by exp(-2.pi.i.(k+G).tau), which a real
            !> array cannot hold. Where tau vanishes the phase is 1 and a real basis is no
            !> obstacle, so the test is on tau and not on the lattice. The other is
            !> non-collinear magnetism, where the rotation is written but not yet right.
            !> The plumbing is in place and measured: the stacked record is carried whole
            !> because calc_abc builds a local spin direction out of BOTH halves of it, the
            !> spin rotation goes into the eigenvector before the coefficients are rebuilt
            !> rather than into the overlap afterwards, each direction reads its own radial
            !> set, and the rule carries the frame factor calc_abc's own line dictates. On
            !> Mn3Ir, whose four moments sit at 120 degrees and out of plane, the identity
            !> comes out at 2.7e-14. What is still red there is the rotation of the other
            !> operations: no leading band set closes, and | D D^+ - 1 | is 9.97e-01. The
            !> identity compares a state with itself, so a clean identity places the fault
            !> in the per-operation rotation and not in the rebuild.
            !>
            !> Fe cannot settle it either way. Its moments lie along z, which is a COLLINEAR
            !> magnet driven through the non-collinear machinery: every spin frame is real
            !> and diagonal, so both placements of a conjugation agree and the identity is
            !> green whatever the frames do. Until a real non-collinear magnet is at the
            !> rounding the matrices are refused, because the rule OVERWRITES the orders
            !> carrying a local orbital and a wrong one there is invisible downstream.
            !> Spin-orbit on its own is green: one potential serves both components, one
            !> frame the crystal, and the 2x2 is the one wl_su2 builds.
            !>
            !> The decision is taken here instead of inside the export, so that a run which
            !> cannot have d_band still gets the stars and the group; the file carries
            !> has_d_band and the reason goes to the log. The real-basis half of the test
            !> reads input%l_real rather than a rule rebuilt from the symmetry: it is the
            !> very flag eigen_hssetup hands to t_mat, so it cannot disagree with what the
            !> eigenvectors turned out to be. Rebuilding it from sym%invs was tried and came
            !> out FALSE for diamond, whose eigenvectors are real.
            l_realph = input%l_real .AND. ANY(ABS(sym%tau(:, :sym%nop)) > 1e-12)
            l_dband = (.NOT. l_realph) .AND. (.NOT. noco%l_noco)
            IF (l_dband) THEN
               CALL wl_build_dband(kstar, manifold, atoms, cell, input, kpts, sym, noco, &
                                   nococonv, enpara, vtot, stars, fmpi, radfun, eig_id, &
                                   jspin, l_wannierlib_spinors, dband)
            ELSE IF (fmpi%irank == 0) THEN
               WRITE (oUnit, '(/,a)') &
                  "Wannier90/WannierBerri symmetry export: d_band NOT written"
               IF (noco%l_noco) &
                  WRITE (oUnit, '(a)') "  the magnetisation is non-collinear"
               IF (l_realph) &
                  WRITE (oUnit, '(a)') &
                  "  the eigenvectors are real and an operation has tau /= 0"
            END IF
         END IF

         !> The stars and the group, in FLEUR's own HDF. A reader on the consumer's side
         !> turns them into that consumer's symmetrizer, so a change in its internal layout
         !> does not reach Fortran. Written here, beside the .amn/.mmn/.eig of the same
         !> channel, because the symmetry tables and the matrices have to share one k-point
         !> ORDER: if the two ever part company every operation is applied at the wrong k and
         !> nothing complains -- it converges, the centres come out symmetric, and the bands
         !> are wrong.
         IF (this%export%symmetrizer) THEN
            CALL wannierlib_export_sym(kstar, manifold, wl, atoms, cell, kpts, sym, &
                                       fmpi, jspin, l_wannierlib_spinors, dband)
            IF (ALLOCATED(dband)) DEALLOCATE (dband)
         END IF
         ! collinear jspins=2 (no SOC/noco): the two spin channels wannierise separately;
         ! tag each channel's interpolation outputs so spin 2 does not overwrite spin 1.
         spin_sfx = ''
         IF (input%jspins == 2 .AND. .NOT. l_wannierlib_spinors) WRITE(spin_sfx, '(a,i0)') '_spin', jspin

         ! wannierise this channel -> the gauge factors u_matrix (MLWF) and u_opt (disentangled)
         !> Published before run_w90 so it is the operator as built, with no gauge on it.
         IF (request%has_op_r('spin') .AND. ALLOCATED(melem%s0) .AND. this%export%bloch) &
            CALL wannierlib_write_s0(fmpi, distk, kpts%nkptf, melem%s0, &
                                     spin12(jspin)//'_s0.dat')
         CALL run_w90(wl, cell, kpts, mmn, amn, eig, fmpi%irank, u_matrix, u_opt)

         !> The gauge goes to results before anything else touches it, so whoever
         !> interpolates with it downstream gets the factors this run produced.
         CALL wannierlib_store_gauge(fmpi, distk, kpts, input, wl%num_bands, wl%num_wann, &
                                     jspin, u_matrix, u_opt, results)

         !> Every real-space operator that is assembled from both channels at once needs
         !> this gauge, and the three of them have to be named here: u_opt and u_matrix are
         !> released at the end of this iteration, so one that is not asked for now cannot
         !> be served later.
         CALL wannierlib_keep_gauge(melem%n_channels, &
                                    request%has_op_r('spin') .OR. &
                                    request%has_op_r('orbital') .OR. &
                                    request%has_op_r('spinOrbit'), manifold, &
                                    kpts%nkptf, jspin, u_opt, u_matrix, v_ch)

         !> Right here and not from v_ch: that one is only assembled on the collinear
         !> two-channel path, while u_opt and u_matrix are in hand on both, so writing from
         !> them covers the spinor route as well.
         IF (this%export%gauge) &
            CALL wannierlib_export_gauge(kpts, fmpi, jspin, u_opt, u_matrix)

         ! Draw the Wannier functions, if anybody asked. Here and not earlier because the
         ! gauge is what makes them: it costs a second pass over this rank's k-points to
         ! read the states back, and what comes out is num_wann grids rather than the
         ! whole basis tabulated at every k.
         IF (this%l_plot_wf) &
            CALL wannierlib_plot_wf(this, atoms, cell, input, sym, stars, vacuum, noco, nococonv, &
                                    enpara, vtot, kpts, fmpi, eig_id, radfun, distk, &
                                    jspin, l_wannierlib_spinors, u_matrix, u_opt)

         ! everything operator-related happens in the matrix-element layer, which only needs
         ! the gauge, the overlaps and the b-mesh. The exception is right below: an operator
         ! whose input needs the wavefunctions AFTER the gauge is known can be built neither
         ! there nor in the post-processing, and only those touch this file.
         ! Kept BEFORE report_w90 so the ordering of the operator messages in `out` is unchanged.
         CALL wannierlib_get_bmesh(this, kpts, bmesh)
         !> F wants the gauge at two neighbours at once, so it can be built neither with the
         !> coarse matrices, which come before any gauge, nor in the post-processing, which
         !> cannot reach the wavefunctions. It goes here, and only when it was asked for.
         IF (request%has_op_r('fmn')) &
            CALL wannierlib_uiu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                                radfun, jspin, l_wannierlib_spinors, eig_id, stars, enpara, &
                                vtot, vacuum, fmpi, distk, u_matrix, u_opt, f0_loc)
         !> C asks for what F asks -- the gauge at two neighbours at once -- and for the
         !> eigenvalues of the bra on top, so it lives here for the same reason.
         IF (request%has_op_r('cmn')) &
            CALL wannierlib_uhu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                                radfun, jspin, l_wannierlib_spinors, eig_id, stars, enpara, &
                                vtot, vacuum, fmpi, distk, u_matrix, u_opt, c0_loc)
         CALL wgauge_run(request, manifold, domains, cell, kpts, eig, u_matrix, u_opt, melem, f0_loc, c0_loc, &
                        mmn, bmesh, distk, fmpi, &
                        wf_channel=jspin, spin_suffix=TRIM(spin_sfx))

         !> After wgauge_run so the message order in `out` keeps the operators together,
         !> and inside the channel loop because a collinear two-channel run lays its
         !> centres down one channel at a time.
         IF (this%export%orbitrans) &
            CALL wannierlib_export_orbitrans(this, atoms, cell, bmesh, fmpi, jspin, melem%n_channels)

         !> Rank 0 owns `out` and `out.xml`, and the centres are plain data on every
         !> rank, so the report is written whatever the number of ranks.
         IF (fmpi%irank == 0) CALL report_w90(this)

         IF (ALLOCATED(amn)) DEALLOCATE (amn)
         IF (ALLOCATED(mmn)) DEALLOCATE (mmn)
         IF (ALLOCATED(u_matrix)) DEALLOCATE (u_matrix)
         IF (ALLOCATED(u_opt)) DEALLOCATE (u_opt)
         IF (ALLOCATED(f0_loc)) DEALLOCATE (f0_loc)
         IF (ALLOCATED(c0_loc)) DEALLOCATE (c0_loc)
         IF (ALLOCATED(eig)) DEALLOCATE (eig)

      END DO

      ! collinear combined 2N spin operator rspauli.1: only assemblable once BOTH channels have
      ! been wannierised, since it rotates the cross-spin overlap with both gauges.
      IF (melem%n_channels == 2 .AND. request%has_op_r('spin')) &
         CALL wgauge_rspauli_collinear(this%num_wann, melem%x0, v_ch, cell, kpts, distk, fmpi)
      ! Same reason for the orbital moment: block-diagonal, but one block per gauge.
      IF (melem%n_channels == 2 .AND. request%has_op_r('orbital')) &
         CALL wgauge_anglmom_collinear(this%num_wann, melem%l0, v_ch, cell, kpts, distk, fmpi)
      ! And the spin-orbit coupling, which here is an operator on the basis rather than a
      ! term in the eigenproblem that produced it.
      IF (melem%n_channels == 2 .AND. request%has_op_r('spinOrbit')) &
         CALL wgauge_soc_collinear(this%num_wann, melem%soc4, v_ch, cell, kpts, distk, fmpi)
      IF (ALLOCATED(v_ch)) DEALLOCATE (v_ch)

      !> Freed here and not per channel: the neighbour topology in it was set before the
      !> spin loop and every channel reads it.
      CALL bmesh%free()
      CALL matrix_element_release_anchor()
      CALL melem%free()
      !> Give back what the factory cached while walking the k-points. It holds the states
      !> and their coefficients for several k at a time, which is the largest thing this
      !> pass leaves behind, and nothing after it reads them.
      CALL matrix_element_factory_reset()
      IF (ALLOCATED(distk)) DEALLOCATE(distk)

   END SUBROUTINE wannierlib_main

   !> Which rank owns each global k-point.
   !>
   !> Prefers the map fmpi already carries, so this pass is distributed exactly like the one
   !> that filled it. Without such a map the k-points are cut into contiguous blocks, and that
   !> is a fallback and not a policy: it ignores whatever the rest of the run decided, and
   !> exists so that a missing map does not stop the calculation.
   SUBROUTINE wannierlib_distribute_k(kpts, fmpi, distk)
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_mpi),  INTENT(IN) :: fmpi
      INTEGER, ALLOCATABLE, INTENT(OUT) :: distk(:)   !> (nkptf) owning rank of each k

      INTEGER :: ikpt, nk_per_rank, ierr

      ALLOCATE(distk(kpts%nkptf), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating distk', &
                                      calledby='wannierlib_distribute_k')
      IF (ALLOCATED(fmpi%coulomb_owner) .AND. SIZE(fmpi%coulomb_owner) == kpts%nkptf) THEN
         distk = fmpi%coulomb_owner
      ELSE
         !> Block size, not the count this rank ends up with -- the two differ on the last rank.
         nk_per_rank = kpts%nkptf/MAX(1, fmpi%isize)
         IF (MOD(kpts%nkptf, MAX(1, fmpi%isize)) > 0) nk_per_rank = nk_per_rank + 1
         DO ikpt = 1, kpts%nkptf
            distk(ikpt) = MIN((ikpt - 1)/MAX(1, nk_per_rank), MAX(0, fmpi%isize - 1))
         END DO
      END IF
   END SUBROUTINE wannierlib_distribute_k

   !> Keep this spin channel's Wannier gauge V = u_opt u_matrix while both factors are alive.
   !>
   !> The operators that are assembled from both channels at once need it, and none of them
   !> can be built until BOTH channels have wannierised, while u_opt and u_matrix are
   !> released at the end of every iteration. Anything that wants a gauge later has to have
   !> kept it here.
   !>
   !> A no-op unless there are two channels and a caller asked for it, so it can be called
   !> unconditionally. Which operators want one is not visible from here, so the caller
   !> decides, and one whose request is left out of that decision finds nothing.
   SUBROUTINE wannierlib_keep_gauge(n_channels, l_want_gauge, manifold, nkptf, jspin, &
                                    u_opt, u_matrix, v_ch)
      INTEGER, INTENT(IN) :: n_channels
      LOGICAL, INTENT(IN) :: l_want_gauge
      TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
      INTEGER, INTENT(IN) :: nkptf, jspin
      COMPLEX, INTENT(IN) :: u_opt(:, :, :), u_matrix(:, :, :)
      COMPLEX, ALLOCATABLE, INTENT(INOUT) :: v_ch(:, :, :, :)   !> (nb, nw, nkptf, 2)

      INTEGER :: ik

      IF (n_channels /= 2 .OR. .NOT. l_want_gauge) RETURN
      IF (.NOT. ALLOCATED(v_ch)) &
         ALLOCATE (v_ch(manifold%num_bands, manifold%num_wann, nkptf, 2), source=CMPLX(0.0, 0.0))
      DO ik = 1, SIZE(u_opt, 3)
         v_ch(:, :, ik, jspin) = MATMUL(u_opt(:, :, ik), u_matrix(:, :, ik))
      END DO
   END SUBROUTINE wannierlib_keep_gauge









END MODULE m_wannierlib_main
