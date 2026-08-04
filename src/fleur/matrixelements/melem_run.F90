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
   USE m_types_wannierlib
   USE m_types_melem_bmesh
   USE m_melem_coarse, ONLY: t_melem_coarse
   USE m_melem_domains, ONLY: melem_write_domain_kpts, melem_rename_domain_outputs, melem_shell
   USE m_melem_operators_r, ONLY: melem_write_operators_r, melem_build_berry_aw_r, melem_check_berry_centres
   USE m_melem_interpolate_ham, ONLY: melem_interpolate_ham
   USE m_melem_interpolate_op, ONLY: melem_interpolate_operator
   USE m_melem_interpolate_velocity, ONLY: melem_interpolate_velocity
   USE m_melem_interpolate_current, ONLY: melem_interpolate_current
   USE m_melem_interpolate_eigenstates, ONLY: melem_interpolate_eigenstates
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: melem_run

CONTAINS

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

      INTEGER :: iop, k, wf_ch, irank, mpi_comm
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
            CASE ('orbital')
               ! total (site-summed) orbital moment
               IF (wann%op_total(iop) == 1) &
                  CALL melem_interpolate_operator(wann, cell, kpts, eig, u_matrix, u_opt, &
                                                  SUM(coarse%l0, DIM=4), gk_loc, 3, 'bands_wann_orbmom.dat', irank, mpi_comm)
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

END MODULE m_melem_run
