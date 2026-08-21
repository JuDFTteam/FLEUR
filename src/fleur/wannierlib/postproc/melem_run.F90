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
MODULE m_melem_run
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_cell
   USE m_types_kpts
   USE m_types_mpi
   USE m_types_melem_bmesh
   USE m_melem_coarse, ONLY: t_melem_coarse
   USE m_types_melem_manifold, ONLY: t_melem_manifold
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_domains, ONLY: t_melem_domains
   USE m_melem_domains, ONLY: melem_write_domain_kpts, melem_rename_domain_outputs, melem_shell
   USE m_melem_operators_r, ONLY: melem_write_operators_r
   USE m_melem_coeff_a, ONLY: melem_build_berry_aw_r, melem_check_berry_centres
   USE m_melem_interpolate_ham, ONLY: melem_interpolate_ham
   USE m_melem_interpolate_op, ONLY: melem_interpolate_operator
   USE m_melem_interpolate_velocity, ONLY: melem_interpolate_velocity
   USE m_melem_interpolate_eigenstates, ONLY: melem_interpolate_eigenstates
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: melem_run

CONTAINS

   SUBROUTINE melem_run(request, manifold, domains, cell, kpts, eig, u_matrix, u_opt, coarse, f0_loc, c0_loc, &
                        mmn, bmesh, distk, fmpi, wf_channel, spin_suffix)
      TYPE(t_melem_request), INTENT(IN) :: request
      TYPE(t_melem_manifold), INTENT(IN) :: manifold
      TYPE(t_melem_domains), INTENT(IN) :: domains
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
      TYPE(t_melem_bmesh), INTENT(IN) :: bmesh         !< b-shell weights (position/velocity operators)
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN), OPTIONAL :: wf_channel            !< collinear spin channel (1/2); default 1
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spin_suffix  !< '_spin1'/'_spin2' collinear; empty otherwise

      INTEGER :: iop, wf_ch, irank, mpi_comm
      INTEGER :: idom, ndom, nkl_c, jkl, aw_nrpts
      INTEGER, ALLOCATABLE :: gk_loc(:), aw_irvec(:, :), aw_ndegen(:)
      COMPLEX, ALLOCATABLE :: aw_r(:, :, :, :)         ! (nw,nw,nrpts,3) Wannier Berry connection A^(W)(R)
      CHARACTER(LEN=16) :: ssfx
      LOGICAL :: l_collinear

      irank = fmpi%irank; mpi_comm = fmpi%mpi_comm
      ssfx = ''
      IF (PRESENT(spin_suffix)) ssfx = TRIM(spin_suffix)
      wf_ch = 1
      IF (PRESENT(wf_channel)) wf_ch = wf_channel
      l_collinear = (LEN_TRIM(ssfx) > 0)   ! collinear jspins=2 -> per-channel operators_r (WF1/WF2)

      CALL timestart('melem_run')

      ! global-k indices owned by this rank, ascending order -> matches the per-rank coarse
      ! slices in `coarse` (built in the same distk order); used by the distributed reduces.
      nkl_c = COUNT(distk == irank); ALLOCATE (gk_loc(nkl_c)); jkl = 0
      DO iop = 1, SIZE(distk)
         IF (distk(iop) == irank) THEN; jkl = jkl + 1; gk_loc(jkl) = iop; END IF
      END DO

      ! ---- output domains: one per <domain>, in the order they were declared ----
      !> A suffixed domain's files are renamed after it runs, so an unsuffixed domain
      !> declared before one that is suffixed would have its base-named output taken. The
      !> reader rejects a repeated suffix, which leaves at most one unsuffixed domain, and
      !> declaring it last is the user's business rather than something to reorder here.
      ndom = domains%n
      !> Only a contradiction is an error: operators to interpolate but nowhere to do it.
      !> This pass runs on every wannierisation, so most calculations arrive here with
      !> neither, and that is not a request that failed -- it is no request at all.
      !>
      !> Asking for operators without saying where used to fall through to a hand-written
      !> k-point file, and without one every driver reported "skipped": the run produced
      !> nothing and said so only in the output file.
      IF (request%n_ops > 0 .AND. ndom == 0) CALL juDFT_error( &
         'wannierlib: <interpolation> asks for operators but declares no output domain', &
         hint='add a <domain> with a listName naming a kPointList', &
         calledby='melem_run')

      ! (2) real-space operator export (Fourier step 3, standalone format) -- once, before interpolation
      CALL melem_write_operators_r(manifold, request, cell, kpts, eig, u_matrix, u_opt, &
                                   coarse%s0, coarse%l0, coarse%soc4, f0_loc, c0_loc, bmesh, distk, mpi_comm, mmn, &
                                   irank, wf_ch, l_collinear)

      ! (3) Wannier-gauge interpolation: dispatch by looping over the requested operator list.
      ! Each operator supplies its own per-rank Bloch slice on the coarse mesh (coarse%s0/l0/soc0);
      ! the remaining steps are the shared generic driver m_melem_interpolate_op.
      DO idom = 1, ndom
         IF (irank == 0) CALL melem_write_domain_kpts(domains, idom)

         DO iop = 1, request%n_ops
            SELECT CASE (TRIM(request%op_name(iop)))
            CASE ('hamiltonian')
               CALL melem_interpolate_ham(manifold, cell, kpts, eig, u_matrix, u_opt, irank)
            CASE ('spin')
               ! total spin (MT-sum + interstitial): via the generic operator driver (3 comps)
               IF (request%op_total(iop) == 1) &
                  CALL melem_interpolate_operator(manifold, cell, kpts, eig, u_matrix, u_opt, &
                                                  coarse%s0, gk_loc, 3, 'bands_wann_spin.dat', irank, mpi_comm)
            CASE ('orbital')
               ! total (site-summed) orbital moment
               IF (request%op_total(iop) == 1) &
                  CALL melem_interpolate_operator(manifold, cell, kpts, eig, u_matrix, u_opt, &
                                                  SUM(coarse%l0(:, :, :, :, wf_ch, :), DIM=4), gk_loc, 3, 'bands_wann_orbmom.dat', irank, mpi_comm)
            CASE ('soc')
               CALL melem_interpolate_operator(manifold, cell, kpts, eig, u_matrix, u_opt, &
                                               coarse%soc0, gk_loc, 1, 'bands_wann_soc.dat', irank, mpi_comm)
            CASE ('velocity')
               ! Wannier Berry connection A^(W)_alpha(R): built distributed from the local overlaps
               ! and reduced (collective, all ranks); the centre check (rank 0) calibrates conj/sign.
               ! Built once and reused across output domains.
               IF (.NOT. ALLOCATED(aw_r)) THEN
                  CALL melem_build_berry_aw_r(manifold, cell, kpts, mmn, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                              aw_r, aw_irvec, aw_ndegen, aw_nrpts)
                  IF (irank == 0) CALL melem_check_berry_centres(manifold, aw_r, aw_irvec, aw_nrpts, bmesh)
               END IF
               CALL melem_interpolate_velocity(manifold, cell, kpts, eig, u_matrix, u_opt, &
                                               aw_r, aw_irvec, aw_ndegen, aw_nrpts, irank)
            CASE ('eigenstates')
               ! Wannier-Hamiltonian eigenvectors C(k') (the H-gauge rotation U^(H)), as a matrix
               CALL melem_interpolate_eigenstates(manifold, cell, kpts, eig, u_matrix, u_opt, irank)
            CASE DEFAULT
               !> The name is in WANNIERLIB_INTERP or it would not have got past the request,
               !> so what is missing is the branch here, not the operator.
               CALL judft_bug('melem_run: "'//TRIM(request%op_name(iop))// &
                  '" is an accepted operator with no branch in this pass')
            END SELECT
         END DO
         !> Rename this domain's outputs. Nested rather than a single .AND.: Fortran does
         !> not promise to stop evaluating at the first false, and domains%suffix exists
         !> only on rank 0.
         IF (irank == 0) THEN
            IF (LEN_TRIM(TRIM(domains%suffix(idom))//TRIM(ssfx)) > 0) &
               CALL melem_rename_domain_outputs(request, TRIM(domains%suffix(idom))//TRIM(ssfx))
         END IF
      END DO   ! idom


      IF (ALLOCATED(aw_r)) DEALLOCATE (aw_r, aw_irvec, aw_ndegen)
      DEALLOCATE (gk_loc)
      CALL timestop('melem_run')
   END SUBROUTINE melem_run

END MODULE m_melem_run
