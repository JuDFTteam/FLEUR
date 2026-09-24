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
   USE m_types, ONLY: t_stars, t_results
   USE m_matrix_element_factory, ONLY: matrix_element_factory_reset, &
                                       matrix_element_release_anchor, matrix_element_radial
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_optable, ONLY: WANNIERLIB_INTERP, WANNIERLIB_OPR, &
                                    melem_exposed_find, melem_exposed_names
   USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
   USE m_types_wgauge_domains, ONLY: t_wgauge_domains
   USE m_wannierlib_build_amn_mmn, ONLY: wannierlib_build_amn_mmn, wannierlib_reduce_amn, &
                                      wannierlib_gather_mmn
   USE m_wannierlib_plot, ONLY: wannierlib_plot_wf
  USE m_wannierlib_uiu, ONLY: wannierlib_uiu
  USE m_wannierlib_uhu, ONLY: wannierlib_uhu
   USE m_wannierlib_w90_adapter
   USE m_melem_coarse, ONLY: t_melem_coarse
   USE m_wgauge_run, ONLY: wgauge_run
   USE m_wannierlib_export_basis, ONLY: wannierlib_export_basis
   USE m_wannierlib_export_gauge, ONLY: wannierlib_export_gauge, wannierlib_store_gauge
USE m_wannierlib_mmnkb, ONLY: wannierlib_kdiff
USE m_wannierlib_band_window, ONLY: wannierlib_default_windows, wannierlib_create_eig
   USE m_wannierlib_export_bloch, ONLY: wannierlib_write_eig, wannierlib_write_s0, &
                                     wannierlib_write_mmn
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
   USE m_types_usdus
   USE m_types_mat
   USE m_types_radfun
   USE m_types_abc
   USE m_types_wannierlib
   use m_wann_write_amn
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
      TYPE(t_usdus), POINTER :: usdus       ! into the factory cache
      TYPE(t_radfun), POINTER :: radfun(:)  ! likewise; the factory owns them
      COMPLEX, ALLOCATABLE :: f0_loc(:, :, :, :, :)  ! (nw,nw,3,3,nk_loc) geometric tensor
      COMPLEX, ALLOCATABLE :: c0_loc(:, :, :, :, :)  ! (nw,nw,3,3,nk_loc) the same one with H inside
      LOGICAL :: l_wannierlib_spinors
      TYPE(t_melem_request) :: request
      !> The catalogue entry each requested name needs built, one per name.
      CHARACTER(LEN=20), ALLOCATABLE :: op_needs(:), op_r_needs(:)
      INTEGER :: iop, krow
      TYPE(t_wgauge_manifold) :: manifold
      TYPE(t_wgauge_domains) :: domains
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
      CALL matrix_element_radial(atoms, input, enpara, fmpi, vtot, radfun, usdus)

      ! Settled before anything else: it distributes the coarse-operator k-loop too.
      CALL wannierlib_distribute_k(kpts, fmpi, distk)

      ! Hand the k-distribution to the matrix-element layer and let it build whatever Bloch-basis
      ! operators the input asks for, in one shared coarse k-pass. Needs only the ab-initio
      ! eigenstates (no gauge), so it runs before the wannierization. All of it is a no-op when
      ! no operator is requested.
      ! what the matrix-element layer is being asked for, and on which bands
      !> The two lists are resolved against the exposure tables HERE, because those tables
      !> say how the wannierization spells things and that is not a fact about the layer
      !> being asked. What it receives is the catalogue entry each name needs built. A name
      !> no table carries stops the run with the accepted names, rather than leaving an
      !> operator silently absent from the output.
      ALLOCATE (op_needs(SIZE(this%ops)), op_r_needs(SIZE(this%op_r_name)))
      DO iop = 1, SIZE(this%ops)
         krow = melem_exposed_find(this%ops(iop)%name, WANNIERLIB_INTERP)
         IF (krow == 0) CALL juDFT_error('wannierlib: "'//TRIM(this%ops(iop)%name)// &
            '" is not an operator that can be interpolated', &
            hint=melem_exposed_names(WANNIERLIB_INTERP), calledby='wannierlib_main')
         op_needs(iop) = WANNIERLIB_INTERP(krow)%operator
      END DO
      DO iop = 1, SIZE(this%op_r_name)
         krow = melem_exposed_find(this%op_r_name(iop), WANNIERLIB_OPR)
         IF (krow == 0) CALL juDFT_error('wannierlib: "'//TRIM(this%op_r_name(iop))// &
            '" has no real-space export', &
            hint=melem_exposed_names(WANNIERLIB_OPR), calledby='wannierlib_main')
         op_r_needs(iop) = WANNIERLIB_OPR(krow)%operator
      END DO

      CALL request%init(this%l_spin, this%l_orbmom, this%l_socop, this%l_operators_r, &
                        this%op_r_name, this%ops%name, this%ops%total, &
                        op_r_needs, op_needs, this%l_ws_distance)
      !> On a film, refuse by name rather than as a whole. The pair overlap has its vacuum
      !> half now, so the wannierization itself works and with it everything built from mmn:
      !> hamiltonian, eigenstates, position, velocity, bmn. Three quantities still stop at
      !> the interstitial and would come out short without saying so.
      IF (input%film) THEN
         !> wannierlib_uiu and wannierlib_uhu ask melem_overlap_states for the pair overlap
         !> without the vacuum expansions, and C also needs a momentum that has no vacuum
         !> part at all.
         IF (request%has_op_r('fmn') .OR. request%has_op_r('cmn')) CALL juDFT_error( &
            'wannierlib: fmn and cmn have no vacuum contribution and cannot be asked for on a film', &
            hint='everything built from the neighbour overlaps does work on a film; these two '// &
                 'need the vacuum halves of the pair overlap and of the momentum', &
            calledby='wannierlib_main')
         !> The spin operator sums the muffin tins and the interstitial. In a film there is a
         !> third region carrying spin density and this layer does not reach it. The orbital
         !> moment and the spin-orbit operator are muffin-tin quantities by construction, not
         !> by omission, so they are unaffected.
         IF (request%has_op('spin') .OR. request%has_op_r('spin')) CALL juDFT_error( &
            'wannierlib: the spin operator has no vacuum contribution and cannot be asked for on a film', &
            hint='it sums the muffin tins and the interstitial only; orbital and spin_orbit are '// &
                 'muffin-tin quantities and do work on a film', &
            calledby='wannierlib_main')
      END IF

      wl = this
      CALL wannierlib_default_windows(wl, results, kpts, input, l_wannierlib_spinors)

      CALL manifold%init(wl%num_bands, wl%num_wann, wl%dis_win_min, wl%dis_win_max, &
                         wl%min_band, wl%max_band)
      CALL domains%init(this%n_domains, this%dom_kset, this%dom_suffix)

      !> Two channels are two eigenproblems, wannierised one after the other, so there is no
      !> single spin matrix over them to interpolate: within a channel sigma_z is +/-1 by
      !> orthonormality, and the transverse part lives in the cross-spin block, which needs
      !> BOTH gauges and therefore both wannierizations. It exists only as the combined 2N
      !> operator in real space. Caught here rather than in the coarse pass, which would
      !> otherwise leave a stub-sized slice to reach the interpolation driver: the summary
      !> flags are set by the real-space list too, so they cannot tell the two lists apart.
      IF (input%jspins == 2 .AND. .NOT.l_wannierlib_spinors) THEN
         IF (request%needs_op('spin', interp_only=.TRUE.)) CALL judft_error( &
            "the spin operator cannot be interpolated when the two spin channels are "// &
            "wannierised separately", &
            hint="ask for it in <operators_r>: there both channels are combined into the "// &
                 "2N rspauli.1, which is the only form this operator has here", &
            calledby="wannierlib_main")
         IF (request%needs_op('spin_orbit', interp_only=.TRUE.)) CALL judft_error( &
            "the spin-orbit operator cannot be interpolated when the two spin channels are "// &
            "wannierised separately", &
            hint="ask for it in <operators_r>: there both channels are combined into the "// &
                 "2N rssocmat.1, which is the only form this operator has here", &
            calledby="wannierlib_main")
      END IF

      CALL melem%init(request, manifold, atoms, input, kpts, fmpi, distk, l_wannierlib_spinors)
      CALL melem%calc(request, manifold, atoms, input, sym, cell, noco, nococonv, kpts, &
                      stars, enpara, fmpi, vtot, eig_id, distk)

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
                                       radfun, usdus, distk, kdiff, nntot_w90, jspin, &
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
         !> and has to be gathered first. Only rank 0 writes, and wann_write_amn is told
         !> isize=1 so it does not try to collect a distribution it does not have.
         IF (this%export%w90 .OR. fmpi%isize == 1) THEN
            CALL wannierlib_gather_mmn(fmpi, distk, kpts%nkptf, mmn, mmn_full)
            IF (fmpi%irank == 0) THEN
               amn_file = spin12(jspin)//'.amn'
               call wann_write_amn(fmpi%mpi_comm, .true., amn_file, "Testing amn", &
                                   this%num_bands, kpts%nkptf, this%num_wann, &
                                   0, 1, .false., .false., &
                                   amn, .false.)

               CALL wannierlib_write_mmn(this, mmn_full, kpts, nnkp, gkpb, jspin)
               CALL wannierlib_write_eig(eig, spin12(jspin)//'.eig')
            END IF
            DEALLOCATE (mmn_full)
         END IF

         !> The interstitial wave functions, for irrep and through it the site-symmetric
         !> wannierisation in WannierBerri. Written here so that it lands next to the
         !> .amn/.mmn/.eig of the same channel and shares their k-point ORDER -- if the two
         !> orders ever part company, every symmetry operation is applied at the wrong k and
         !> the failure is silent: it converges, the centres come out symmetric, and the
         !> bands are wrong.
         IF (this%export%basis) &
            CALL wannierlib_export_basis(this, manifold, atoms, cell, input, kpts, sym, &
                                         noco, nococonv, enpara, vtot, fmpi, eig_id, jspin)
         ! collinear jspins=2 (no SOC/noco): the two spin channels wannierise separately;
         ! tag each channel's interpolation outputs so spin 2 does not overwrite spin 1.
         spin_sfx = ''
         IF (input%jspins == 2 .AND. .NOT. l_wannierlib_spinors) WRITE(spin_sfx, '(a,i0)') '_spin', jspin

         ! wannierise this channel -> the gauge factors u_matrix (MLWF) and u_opt (disentangled)
         !> The spin operator goes down only when it was actually asked for: it is allocated
         !> whenever melem is active but filled only on request, and a zero array would make
         !> the band labelling refuse every case instead of falling back to the orbital
         !> criterion, which is the right answer when there is no spin operator to use.
         IF (request%has_op_r('spin') .AND. ALLOCATED(melem%s0)) THEN
            !> Published before run_w90 so it is the operator as built, with no gauge on it.
            IF (this%export%bloch) &
               CALL wannierlib_write_s0(fmpi, distk, kpts%nkptf, melem%s0, &
                                        spin12(jspin)//'_s0.dat')
            CALL run_w90(wl, cell, kpts, mmn, amn, eig, fmpi%irank, u_matrix, u_opt, &
                         s0=melem%s0)
         ELSE
            CALL run_w90(wl, cell, kpts, mmn, amn, eig, fmpi%irank, u_matrix, u_opt)
         END IF

         !> The gauge goes to results before anything else touches it, so whoever
         !> interpolates with it downstream gets the factors this run produced.
         CALL wannierlib_store_gauge(fmpi, distk, kpts, input, wl%num_bands, wl%num_wann, &
                                     jspin, u_matrix, u_opt, results)

         CALL wannierlib_keep_gauge(melem%n_channels, request%has_op_r('spin'), manifold, &
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
                                vtot, fmpi, distk, u_matrix, u_opt, f0_loc)
         !> C asks for what F asks -- the gauge at two neighbours at once -- and for the
         !> eigenvalues of the bra on top, so it lives here for the same reason.
         IF (request%has_op_r('cmn')) &
            CALL wannierlib_uhu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                                radfun, jspin, l_wannierlib_spinors, eig_id, stars, enpara, &
                                vtot, fmpi, distk, u_matrix, u_opt, c0_loc)
         CALL wgauge_run(request, manifold, domains, cell, kpts, eig, u_matrix, u_opt, melem, f0_loc, c0_loc, &
                        mmn, bmesh, distk, fmpi, &
                        wf_channel=jspin, spin_suffix=TRIM(spin_sfx))

         if (fmpi%isize == 1) CALL report_w90(this)

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
      IF (melem%n_channels == 2 .AND. request%has_op_r('spin_orbit')) &
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
   !> Only the collinear two-channel spin operator needs it, and only it can say so: the
   !> combined 2N matrix cannot be assembled until BOTH channels have wannierised, while
   !> u_opt and u_matrix are released at the end of every iteration. Anything that wants a
   !> gauge later has to have kept it here.
   !>
   !> A no-op unless there are two channels and the spin operator was asked for, so the
   !> caller can call it unconditionally.
   SUBROUTINE wannierlib_keep_gauge(n_channels, l_want_spin, manifold, nkptf, jspin, &
                                    u_opt, u_matrix, v_ch)
      INTEGER, INTENT(IN) :: n_channels
      LOGICAL, INTENT(IN) :: l_want_spin
      TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
      INTEGER, INTENT(IN) :: nkptf, jspin
      COMPLEX, INTENT(IN) :: u_opt(:, :, :), u_matrix(:, :, :)
      COMPLEX, ALLOCATABLE, INTENT(INOUT) :: v_ch(:, :, :, :)   !> (nb, nw, nkptf, 2)

      INTEGER :: ik

      IF (n_channels /= 2 .OR. .NOT. l_want_spin) RETURN
      IF (.NOT. ALLOCATED(v_ch)) &
         ALLOCATE (v_ch(manifold%num_bands, manifold%num_wann, nkptf, 2), source=CMPLX(0.0, 0.0))
      DO ik = 1, SIZE(u_opt, 3)
         v_ch(:, :, ik, jspin) = MATMUL(u_opt(:, :, ik), u_matrix(:, :, ik))
      END DO
   END SUBROUTINE wannierlib_keep_gauge









END MODULE m_wannierlib_main
