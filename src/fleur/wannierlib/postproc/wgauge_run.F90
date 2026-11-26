!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Everything that needs the Wannier gauge V = u_opt . u_matrix: the real-space O(R)
!>  export, and the band interpolation of every requested operator over every requested
!>  output domain.
!>
!>  It is driven by a list of operator names and a list of output domains, which is why it
!>  belongs with the matrix elements and not with the wannierization: only the interpolation
!>  step is specific to the Wannier gauge -- the export formats and the output domains are
!>  the same service for any caller.
MODULE m_wgauge_run
   USE m_juDFT
   USE m_constants, ONLY: oUnit, hartree_to_ev_const
   USE m_types_cell
   USE m_types_kpts
   USE m_types_mpi
   USE m_types_wgauge_bmesh
   USE m_melem_coarse, ONLY: t_melem_coarse
   USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_wgauge_domains, ONLY: t_wgauge_domains
   USE m_wgauge_operators_r, ONLY: wgauge_write_operators_r
   USE m_wgauge_coeff_b, ONLY: wgauge_build_bk
   USE m_wgauge_coeff_a, ONLY: wgauge_build_berry_aw_r, wgauge_check_berry_centres
   USE m_wgauge_hamk, ONLY: t_wgauge_hgauge, wgauge_build_hamk, wgauge_hgauge_domain
   USE m_wgauge_bands_io, ONLY: WL_BANDS_HID, t_wgauge_bands_file, &
                                wgauge_bands_file_open, wgauge_bands_file_close, &
                                wgauge_bands_domain, wgauge_bands_domain_close
   USE m_wgauge_interpolate_ham, ONLY: wgauge_interpolate_ham
   USE m_wgauge_interpolate_op, ONLY: wgauge_interpolate_operator
   USE m_wgauge_interpolate_velocity, ONLY: wgauge_interpolate_velocity
   USE m_wgauge_ft, ONLY: wgauge_mdrs_set, wgauge_mdrs_clear, wgauge_ft_to_real
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wgauge_run

CONTAINS

   SUBROUTINE wgauge_run(request, manifold, domains, cell, kpts, eig, u_matrix, u_opt, coarse, f0_loc, c0_loc, &
                        mmn, bmesh, distk, fmpi, wf_channel, spin_suffix, b0k_out, gk_loc_out)
      TYPE(t_melem_request), INTENT(IN) :: request
      TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
      TYPE(t_wgauge_domains), INTENT(IN) :: domains
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      REAL, INTENT(IN) :: eig(:, :)                    !< (nb,nk)
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)         !< (nw,nw,nk) MLWF gauge
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)            !< (nb,nw,nk) disentangled
      TYPE(t_melem_coarse), INTENT(IN) :: coarse
      !> (nw,nw,3,3,nk_loc) the geometric tensor. Unallocated unless it was asked for: it is
      !> the one operator that cannot be built from what this routine is handed.
      COMPLEX, ALLOCATABLE, INTENT(IN) :: f0_loc(:, :, :, :, :)
      !> The same one, with the Hamiltonian inside.
      COMPLEX, ALLOCATABLE, INTENT(IN) :: c0_loc(:, :, :, :, :)
      COMPLEX, INTENT(IN) :: mmn(:, :, :, :)           !< (nb,nb,nntot,nk_loc) this rank's overlap slice
      TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh         !< b-shell weights (position/velocity operators)
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN), OPTIONAL :: wf_channel            !< collinear spin channel (1/2); default 1
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spin_suffix  !< '_spin1'/'_spin2' collinear; empty otherwise
      !> B(k) and this rank's global-k list, handed up for the <operatorsKs> file. The file
      !> is written by the driver because the writer lives in export/, which this directory
      !> does not reach into: what belongs here is building B(k), not writing it.
      COMPLEX, ALLOCATABLE, INTENT(OUT), OPTIONAL :: b0k_out(:, :, :, :)
      INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL :: gk_loc_out(:)

      INTEGER :: iop, wf_ch, irank, mpi_comm
      INTEGER :: idom, ndom, nkl_c, jkl, aw_nrpts
      !> The one interpolation file of this spin channel, and the group of the domain being
      !> written into it. Rank 0 only: no other rank opens either.
      TYPE(t_wgauge_bands_file) :: ifh
      INTEGER(WL_BANDS_HID) :: dgid
      CHARACTER(LEN=40) :: dgrp
      INTEGER, ALLOCATABLE :: gk_loc(:), aw_irvec(:, :), aw_ndegen(:)
      !> The three-component spin slice of one collinear channel, assembled from the
      !> channel's own block so the generic driver takes it like any other operator.
      COMPLEX, ALLOCATABLE :: s_ch(:, :, :, :)
      COMPLEX, ALLOCATABLE :: aw_r(:, :, :, :)         ! (nw,nw,nrpts,3) Wannier Berry connection A^(W)(R)
      COMPLEX, ALLOCATABLE :: b0k(:, :, :, :)          ! (nw,nw,3,nk_loc) B(k) for the k-space export
      !> H_W(k) and H(R): the Hamiltonian side of every interpolation, built once here.
      COMPLEX, ALLOCATABLE :: ham_k(:, :, :), ham_r(:, :, :)
      INTEGER, ALLOCATABLE :: h_irvec(:, :), h_ndegen(:)
      INTEGER :: h_nrpts
      TYPE(t_wgauge_hgauge) :: hg                      ! the current domain's bands and C(k')
      CHARACTER(LEN=16) :: ssfx
      CHARACTER(LEN=80) :: dsfx
      INTEGER :: np_dom
      REAL, ALLOCATABLE :: kfrac(:, :)
      LOGICAL :: l_collinear

      irank = fmpi%irank; mpi_comm = fmpi%mpi_comm
      ssfx = ''
      IF (PRESENT(spin_suffix)) ssfx = TRIM(spin_suffix)
      wf_ch = 1
      IF (PRESENT(wf_channel)) wf_ch = wf_channel
      l_collinear = (LEN_TRIM(ssfx) > 0)   ! collinear jspins=2 -> per-channel operatorsRs (WF1/WF2)

      CALL timestart('wgauge_run')

      ! global-k indices owned by this rank, ascending order -> matches the per-rank coarse
      ! slices in `coarse` (built in the same distk order); used by the distributed reduces.
      nkl_c = COUNT(distk == irank); ALLOCATE (gk_loc(nkl_c)); jkl = 0
      DO iop = 1, SIZE(distk)
         IF (distk(iop) == irank) THEN; jkl = jkl + 1; gk_loc(jkl) = iop; END IF
      END DO
      IF (PRESENT(gk_loc_out)) gk_loc_out = gk_loc

      ! ---- output domains: one per <domain>, in the order they were declared ----
      !> A suffixed domain's files are renamed after it runs, so an unsuffixed domain
      !> declared before one that is suffixed would have its base-named output taken. The
      !> reader rejects a repeated suffix, which leaves at most one unsuffixed domain, and
      !> declaring it last is the user's business rather than something to reorder here.
      ndom = domains%n
      !> Only a contradiction is an error: operators to interpolate but nowhere to do it.
      !> This pass runs on every wannierisation, so most calculations arrive here with
      !> neither, and that is not a request that failed -- it is no request at all.
      IF (request%n_ops > 0 .AND. ndom == 0) CALL juDFT_error( &
         'wannierlib: <interpolation> asks for operators but declares no output domain', &
         hint='add a <domain> with a listName naming a kPointList', &
         calledby='wgauge_run')

      ! (2) real-space operator export (Fourier step 3, standalone format) -- once, before interpolation
      CALL wgauge_write_operators_r(manifold, request, cell, kpts, eig, u_matrix, u_opt, &
                                   coarse%s0, coarse%l0, coarse%soc4, coarse%t0, coarse%hxc0, f0_loc, c0_loc, bmesh, distk, mpi_comm, mmn, &
                                   irank, wf_ch, l_collinear)

      !> (2b) <operatorsKs> can also ask for B, C and F. They are not in the Bloch basis
      !> the other half of that block promises: contracting over PAIRS of neighbours needs
      !> a band index that means the same state at k+b1 and at k+b2, and the Bloch one does
      !> not, so these three carry the Wannier gauge and go to a file of their own. The
      !> request is independent of <operatorsRs> -- either block may ask for a tensor and
      !> neither obliges the other -- so this runs whether or not the pass above did.
      !>
      !> C and F arrive already built, so only B is made here; the three travel up to the
      !> driver, which writes them.
      IF (request%has_op_k('bmn')) &
         CALL wgauge_build_bk(manifold, kpts, eig, u_matrix, u_opt, mmn, gk_loc, bmesh, &
                             b0k, irank)
      IF (PRESENT(b0k_out) .AND. ALLOCATED(b0k)) CALL MOVE_ALLOC(b0k, b0k_out)
      IF (ALLOCATED(b0k)) DEALLOCATE (b0k)

      !> MDRS, switched on once for every operator and every output domain rather than per
      !> call: interpolating H with it and an operator without it would put the two in
      !> different gauges, and the projection of one on the other is what gets written.
      !>
      !> Deliberately AFTER the O(R) export: what wgauge_write_operators_r writes is H(R) in
      !> the plain Wannier90 convention, which is what its readers expect -- Wannier90 itself
      !> stores _hr.dat that way and applies the replica average when it reads it back.
      IF (request%l_ws_distance) THEN
         IF (.NOT. ALLOCATED(bmesh%centres)) CALL juDFT_error( &
            'wannierlib: useWsDistance needs the Wannier centres and none were reported', &
            hint='it requires the Wannier90 library build (CPP_WANNLIB_API)', &
            calledby='wgauge_run')
         CALL wgauge_mdrs_set(cell, kpts%nkpt3, bmesh%centres, irank)
      END IF

      !> H_W(k) and its transform to real space. Neither depends on the output domain, and
      !> H(R) is what the band, eigenstate, operator and velocity drivers all read, so it is
      !> built once: the convention the Hamiltonian is assembled in then has a single place
      !> to reach, and the four drivers cannot drift apart.
      !>
      !> After the MDRS switch, never before it: every transform back to k reads those
      !> tables, so anything interpolated ahead of the switch would be in the other gauge.
      !>
      !> Rank 0 only, which is where the full U(k) lives and where every driver interpolates.
      !> ndom and not n_ops: what needs H(R) is the loop over domains below, which runs
      !> whether or not any quantity was asked for on them. An <interpolation> with a domain
      !> and no <operator> used to reach wgauge_hgauge_domain with ham_r unallocated, and
      !> that dummy is an assumed-shape array rather than an allocatable, so it was read
      !> before it existed. Asking for no operator is how a run says it wants the bands and
      !> nothing else, so it has to work.
      IF (irank == 0 .AND. ndom > 0) THEN
         CALL wgauge_build_hamk(manifold, eig, u_matrix, u_opt, ham_k)
         CALL wgauge_ft_to_real(cell, ham_k, kpts, ham_r, h_irvec, h_ndegen, h_nrpts)
         DEALLOCATE (ham_k)
      END IF

      ! (3) Wannier-gauge interpolation: dispatch by looping over the requested operator list.
      ! Each operator supplies its own per-rank Bloch slice on the coarse mesh (coarse%s0/l0/soc0);
      ! the remaining steps are the shared generic driver m_wgauge_interpolate_op.
      !> One file for the whole pass, holding every domain and every quantity asked for on
      !> it. Opened before the loop, because a domain is a group inside it.
      IF (irank == 0) CALL wgauge_bands_file_open(wf_ch, kpts%nkpt, ifh)

      DO idom = 1, ndom
         !> The domain's k-points and the tail its files carry. Both are rank-0 material:
         !> kfrac stays unallocated elsewhere and every driver returns before touching it,
         !> while the suffix is broadcast because it is small and the name is built here.
         IF (ALLOCATED(kfrac)) DEALLOCATE(kfrac)
         IF (irank == 0) THEN
            np_dom = domains%kset(idom)%nkpt
            ALLOCATE(kfrac(3, np_dom))
            kfrac = domains%kset(idom)%bk(:, 1:np_dom)
         END IF
         !> H(k') and its eigenvectors for this domain, once for the operators that follow.
         IF (irank == 0) CALL wgauge_hgauge_domain(cell, kfrac, ham_r, h_irvec, h_ndegen, h_nrpts, hg)

         dsfx = TRIM(domains%suffix(idom))//TRIM(ssfx)
         !> The domain's own group: its abscissa, its band energies and where its
         !> high-symmetry points sit along it. Written once here, so that no quantity below
         !> repeats them and whoever plots the bands derives nothing from kpts.xml.
         !> Named for the domain and not for the spin channel: the file is already
         !> WF<channel>_interp.hdf, so carrying the channel in the group as well would put
         !> it in the path twice. The domain that declares no suffix is the base one, and a
         !> group cannot be called nothing, so it is the one that takes 'default'.
         dgrp = TRIM(domains%suffix(idom))
         IF (dgrp(1:1) == '_') dgrp = dgrp(2:)
         IF (LEN_TRIM(dgrp) == 0) dgrp = 'default'
         IF (irank == 0) CALL wgauge_bands_domain(ifh, TRIM(dgrp), TRIM(domains%suffix(idom)), &
                                                  domains%kset(idom), &
                                                 hg%kdist, hg%evals, dgid)

         DO iop = 1, request%n_ops
            SELECT CASE (TRIM(request%op_name(iop)))
            CASE ('hamiltonian')
               CALL wgauge_interpolate_ham(hg, irank)
            CASE ('spin')
               ! total spin (MT-sum + interstitial): via the generic operator driver (3 comps)
               ! bound=1: <sigma> of a normalised spinor cannot exceed 1, so the driver can
               ! say so when the interpolation overshoots. Only spin has a bound this simple.
               IF (l_collinear) THEN
                  !> This channel's own block in z, and zero in x and y. The zeros are the
                  !> answer and not an omission: the transverse expectation value of a
                  !> spin-pure eigenstate is zero, and collinearly every eigenstate is one.
                  IF (.NOT. ALLOCATED(coarse%sz0)) CALL juDFT_error( &
                     'wannierlib: the collinear spin channel block was never built', &
                     hint='melem_coarse allocates sz0 when spin is requested; it was not', &
                     calledby='wgauge_run')
                  ALLOCATE (s_ch(SIZE(coarse%sz0, 1), SIZE(coarse%sz0, 2), 3, SIZE(coarse%sz0, 4)), &
                            source=CMPLX(0.0, 0.0))
                  s_ch(:, :, 3, :) = coarse%sz0(:, :, wf_ch, :)
                  CALL wgauge_interpolate_operator(manifold, cell, kpts, u_matrix, u_opt, &
                                                  s_ch, gk_loc, 3, kfrac, hg, dgid, 'spin', 'hbar/2', irank, mpi_comm, &
                                                  bound=1.0)
                  DEALLOCATE (s_ch)
               ELSE
                  CALL wgauge_interpolate_operator(manifold, cell, kpts, u_matrix, u_opt, &
                                                  coarse%s0, gk_loc, 3, kfrac, hg, dgid, 'spin', 'hbar/2', irank, mpi_comm, &
                                                  bound=1.0)
               END IF
            CASE ('orbital')
               ! total (site-summed) orbital moment
               CALL wgauge_interpolate_operator(manifold, cell, kpts, u_matrix, u_opt, &
                                               SUM(coarse%l0(:, :, :, :, wf_ch, :), DIM=4), gk_loc, 3, kfrac, &
                                               hg, dgid, 'orbital', 'hbar', irank, mpi_comm)
            CASE ('spinOrbit')
               CALL wgauge_interpolate_operator(manifold, cell, kpts, u_matrix, u_opt, &
                                               coarse%soc0, gk_loc, 1, kfrac, hg, dgid, 'spinOrbit', 'eV', irank, mpi_comm, &
                                               scale=hartree_to_ev_const)
            CASE ('exchangeTorque')
               ! T = -(s x m), three Cartesian components like the spin, and like it
               ! only defined with spinors -- the collinear case is refused upstream.
               CALL wgauge_interpolate_operator(manifold, cell, kpts, u_matrix, u_opt, &
                                               coarse%t0, gk_loc, 3, kfrac, hg, dgid, 'exchangeTorque', 'eV', irank, mpi_comm, &
                                               scale=hartree_to_ev_const)
            CASE ('velocity')
               ! Wannier Berry connection A^(W)_alpha(R): built distributed from the local overlaps
               ! and reduced (collective, all ranks); the centre check (rank 0) calibrates conj/sign.
               ! Built once and reused across output domains.
               IF (.NOT. ALLOCATED(aw_r)) THEN
                  CALL wgauge_build_berry_aw_r(manifold, cell, kpts, mmn, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                              aw_r, aw_irvec, aw_ndegen, aw_nrpts)
                  IF (irank == 0) CALL wgauge_check_berry_centres(manifold, aw_r, aw_irvec, aw_nrpts, bmesh)
               END IF
               CALL wgauge_interpolate_velocity(cell, ham_r, h_irvec, h_ndegen, h_nrpts, &
                                               aw_r, aw_irvec, aw_ndegen, aw_nrpts, kfrac, hg, &
                                               dgid, irank)
            CASE DEFAULT
               !> The name is in WANNIERLIB_INTERP or it would not have got past the request,
               !> so what is missing is the branch here, not the operator.
               CALL judft_bug('wgauge_run: "'//TRIM(request%op_name(iop))// &
                  '" is an accepted operator with no branch in this pass')
            END SELECT
         END DO
         IF (irank == 0) CALL wgauge_bands_domain_close(dgid)
      END DO   ! idom
      IF (irank == 0) THEN
         CALL wgauge_bands_file_close(ifh)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib interpolation: wrote '//TRIM(ifh%name)// &
            ' (', ndom, ' domain(s))'
      END IF


      CALL wgauge_mdrs_clear()   ! the state must not outlive the wannierization it describes
      IF (ALLOCATED(aw_r)) DEALLOCATE (aw_r, aw_irvec, aw_ndegen)
      IF (ALLOCATED(ham_r)) DEALLOCATE (ham_r, h_irvec, h_ndegen)
      DEALLOCATE (gk_loc)
      CALL timestop('wgauge_run')
   END SUBROUTINE wgauge_run

END MODULE m_wgauge_run
