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
   USE m_types, ONLY: t_stars, t_results
   USE m_matrix_element_factory, ONLY: matrix_element_factory_reset, matrix_element_states, &
                                       matrix_element_release_anchor, matrix_element_radial
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_manifold, ONLY: t_melem_manifold
   USE m_types_melem_domains, ONLY: t_melem_domains
   USE m_wannierlib_amn
   USE m_wannierlib_mmnkb
   USE m_melem_ujugaunt
  USE m_wannierlib_uiu, ONLY: wannierlib_uiu
  USE m_wannierlib_uhu, ONLY: wannierlib_uhu
   USE m_wannierlib_w90_adapter
   USE m_melem_coarse, ONLY: t_melem_coarse
   USE m_melem_run, ONLY: melem_run
   USE m_melem_spin_collinear, ONLY: melem_rspauli_collinear
   USE m_types_melem_bmesh, ONLY: t_melem_bmesh
   USE m_constants, ONLY: oUnit
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
   USE m_types_usdus
   USE m_types_mat
   USE m_types_radfun
   USE m_types_spinor_layout, ONLY: radial_slot
   USE m_types_abc
   USE m_types_wannierlib

   use m_wann_write_amn
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_main
CONTAINS

   SUBROUTINE wannierlib_main(this, atoms, cell, input, kpts, sym, noco, nococonv, stars, enpara, fmpi, vtot, results, eig_id)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_potden), INTENT(IN) :: vtot
      TYPE(t_results), INTENT(IN) :: results
      INTEGER, INTENT(IN) :: eig_id

      INTEGER :: ikpt, nntot_w90, ierr, jspin, jspin_comp, irec, ib
      INTEGER :: jspin_rad_done   ! radial set the cached ujug was built for; -1 = none yet
      INTEGER, ALLOCATABLE :: ev_list(:)
      COMPLEX, ALLOCATABLE :: amn(:, :, :)
      COMPLEX, ALLOCATABLE :: mmn(:, :, :, :)
      COMPLEX, ALLOCATABLE :: ujug(:, :, :, :, :, :)
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
      INTEGER :: ik_gauge
      TYPE(t_melem_coarse) :: melem   ! the operator (matrix-element) side
      TYPE(t_melem_bmesh) :: bmesh    ! b-shell weights handed to the operator side
      TYPE(t_usdus), POINTER :: usdus       ! into the factory cache
      TYPE(t_lapw) :: lapw
      TYPE(t_mat), POINTER :: zmat_p(:)     ! into the factory cache
      TYPE(t_radfun), POINTER :: radfun(:)  ! likewise; the factory owns them
      COMPLEX, ALLOCATABLE :: f0_loc(:, :, :, :, :)  ! (nw,nw,3,3,nk_loc) geometric tensor
      COMPLEX, ALLOCATABLE :: c0_loc(:, :, :, :, :)  ! (nw,nw,3,3,nk_loc) el mismo con H dentro
      TYPE(t_abc), POINTER :: abc_p(:, :)   ! (2,ntype), likewise
      LOGICAL :: l_wannierlib_spinors
      TYPE(t_melem_request) :: request
      TYPE(t_melem_manifold) :: manifold
      TYPE(t_melem_domains) :: domains
      LOGICAL :: l_nocosoc
      INTEGER :: jspin_rad
      CHARACTER(LEN=7) :: amn_file
      CHARACTER(LEN=3) :: spin12(2)
      INTEGER :: ik_local, nk_local
      CHARACTER(LEN=6) :: spin_sfx

      IF (.NOT. this%l_wannierize) RETURN

      l_wannierlib_spinors = noco%l_noco .OR. noco%l_soc
      l_nocosoc = noco%l_noco .AND. (.NOT. noco%l_soc)
      ! A.1: input%l_real stays TRUE with inversion even under SOC (n_denmat=0), and
      ! reading a complex spinor into a real buffer kills the imaginary part of the MMN.
      spin12 = (/'WF1', 'WF2'/)

      !> The radial functions come from the factory, which keeps them for the operators
      !> anyway. Generating a second set here produced the same numbers twice.
      CALL matrix_element_radial(atoms, input, enpara, fmpi, vtot, radfun, usdus)

      ! distk: which rank owns each global k (moved up: distributes the coarse-operator k-loop)
      ALLOCATE(distk(kpts%nkptf), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating distk', calledby='wannierlib_main')
      IF (ALLOCATED(fmpi%coulomb_owner) .AND. SIZE(fmpi%coulomb_owner) == kpts%nkptf) THEN
         distk = fmpi%coulomb_owner
      ELSE
         ! Fallback to contiguous block distribution if no global owner map is available.
         nk_local = kpts%nkptf/MAX(1, fmpi%isize)
         IF (MOD(kpts%nkptf, MAX(1, fmpi%isize)) > 0) nk_local = nk_local + 1
         DO ikpt = 1, kpts%nkptf
            distk(ikpt) = MIN((ikpt - 1)/MAX(1, nk_local), MAX(0, fmpi%isize - 1))
         END DO
      END IF

      ! Hand the k-distribution to the matrix-element layer and let it build whatever Bloch-basis
      ! operators the input asks for, in one shared coarse k-pass. Needs only the ab-initio
      ! eigenstates (no gauge), so it runs before the wannierization. All of it is a no-op when
      ! no operator is requested.
      ! what the matrix-element layer is being asked for, and on which bands
      CALL request%init(this%l_spin, this%l_orbmom, this%l_socop, this%l_operators_r, &
                        this%op_r_name, this%op_name, this%op_total)
      CALL manifold%init(this%num_bands, this%num_wann, this%dis_win_min, this%dis_win_max, &
                         this%min_band, this%max_band)
      CALL domains%init(this%l_dom_path, this%l_dom_plane, this%l_dom_grid, this%path_file, &
                        this%path_kset, this%plane_kset, this%grid_kset)

      CALL melem%init(request, manifold, atoms, input, kpts, fmpi, distk, l_wannierlib_spinors)
      CALL melem%calc(request, manifold, atoms, input, sym, cell, noco, nococonv, kpts, &
                      stars, enpara, fmpi, vtot, eig_id, distk)

      CALL init_w90(this, atoms, cell, kpts, fmpi, l_wannierlib_spinors, nntot_w90, nnkp, gkpb, distk)
      CALL wannierlib_kdiff(kpts%nkptf, nntot_w90, kpts%bkf, nnkp, gkpb, kdiff)
      !> The neighbour topology is settled here, before anything is wannierised, and it is
      !> the same for both spin channels. The shell weights are added later, per channel.
      CALL bmesh%set_neighbours(nntot_w90, nnkp, gkpb, kdiff)

      DO jspin = 1, MERGE(2, 1, input%jspins == 2 .AND. (.NOT. l_wannierlib_spinors))

         ! calculate the  matrices for all k-points
         ALLOCATE (amn(this%num_bands, this%num_wann, kpts%nkptf), stat=ierr, source=cmplx(0.0, 0.0))
         IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating amn buffer', calledby='wannierlib_main')

         nk_local = COUNT(distk == fmpi%irank)
         ALLOCATE(mmn(this%num_bands, this%num_bands, nntot_w90, nk_local), stat=ierr, source=cmplx(0.0, 0.0))
         IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating local mmn buffer', calledby='wannierlib_main')

         jspin_rad_done = -1
         DO jspin_comp = MERGE(1, jspin, l_wannierlib_spinors), MERGE(2, jspin, l_wannierlib_spinors)
            ! jspin_comp = the eig record (spinor up/down); jspin_rad = the radial index.
            jspin_rad = radial_slot(radfun, jspin_comp)
            !> These depend on the radial set and on nothing else in this loop. With a single
            !> potential both spinor components read the same set, so the second pass would
            !> rebuild an identical array -- and it is the largest thing built here, lm by lm
            !> by radial pair by type by neighbour.
            IF (jspin_rad /= jspin_rad_done) THEN
               CALL melem_ujugaunt(atoms, cell, nntot_w90, kdiff, radfun, radfun, jspin_rad, jspin_rad, .FALSE., 1, ujug)
               jspin_rad_done = jspin_rad
            END IF

            ev_list = [(ib, ib = this%min_band, this%max_band)]
            ik_local = 0
            DO ikpt = 1, kpts%nkptf
               IF (distk(ikpt) /= fmpi%irank) CYCLE   ! each rank computes only its k-slice -> parallel eigenvector I/O
               !> This k stays put while its neighbours come and go: mmnkb asks the factory
               !> for one per b, and without the anchor this one would be the oldest by the
               !> third of them and be overwritten under the pointers held here.
               CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)
               CALL matrix_element_states(eig_id, ikpt, input, atoms, sym, cell, noco, &
                                          nococonv, enpara, lapw, vtot, fmpi, zmat_p, abc_p, &
                                          ev_list=ev_list, &
                                          l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                          kpts=kpts, l_anchor=.TRUE.)
               !> Non-collinearly the whole spinor is one record; otherwise each channel is
               !> its own and the spin block is reached by row offset further down.
               irec = MERGE(1, jspin_comp, noco%l_noco)

               CALL wannierlib_amn(this, atoms, kpts, ikpt, usdus, radfun, abc_p(jspin_comp, :), l_nocosoc, jspin_comp, jspin_rad, amn(:, :, ikpt))

               ik_local = ik_local + 1
               CALL wannierlib_mmnkb(manifold, bmesh, ikpt, kpts, &
                                     ujug, atoms, cell, input, sym, noco, nococonv, &
                                     abc_p(jspin_comp, :), jspin_comp, jspin_rad, eig_id, stars, lapw, &
                                     zmat_p(irec), mmn, ik_local, enpara, vtot, fmpi)
            END DO

         END DO
         IF (ALLOCATED(ujug)) DEALLOCATE (ujug)

         ! amn was filled only on each rank's distk slice (zeros elsewhere) -> sum to the full set
         CALL wannierlib_reduce_amn(fmpi, amn)

         mmn = conjg(mmn)

         IF (fmpi%isize == 1) THEN
            amn_file = spin12(jspin)//'.amn'
            call wann_write_amn(fmpi%mpi_comm, .true., amn_file, "Testing amn", &
                                this%num_bands, kpts%nkptf, this%num_wann, &
                                0, 1, .false., .false., &
                                amn, .false.)

            CALL wannierlib_write_mmn(this, mmn, kpts, nnkp, gkpb, jspin)
         END IF

         call wannierlib_create_eig(this, results, kpts, MERGE(1, jspin, l_wannierlib_spinors), eig)
         ! collinear jspins=2 (no SOC/noco): the two spin channels wannierise separately;
         ! tag each channel's interpolation outputs so spin 2 does not overwrite spin 1.
         spin_sfx = ''
         IF (input%jspins == 2 .AND. .NOT. l_wannierlib_spinors) WRITE(spin_sfx, '(a,i0)') '_spin', jspin

         ! wannierise this channel -> the gauge factors u_matrix (MLWF) and u_opt (disentangled)
         CALL run_w90(this, cell, kpts, mmn, amn, eig, fmpi%irank, u_matrix, u_opt)

         !> Keep this channel's gauge while its two factors are still alive.
         IF (melem%n_channels == 2 .AND. request%has_op_r('spin')) THEN
            IF (.NOT. ALLOCATED(v_ch)) &
               ALLOCATE (v_ch(manifold%num_bands, manifold%num_wann, kpts%nkptf, 2), &
                         source=CMPLX(0.0, 0.0))
            DO ik_gauge = 1, SIZE(u_opt, 3)
               v_ch(:, :, ik_gauge, jspin) = MATMUL(u_opt(:, :, ik_gauge), u_matrix(:, :, ik_gauge))
            END DO
         END IF

         ! everything operator-related happens in the matrix-element layer, which only needs the
         ! gauge, the overlaps and the b-mesh. Adding an operator does not touch this file.
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
         CALL melem_run(request, manifold, domains, cell, kpts, eig, u_matrix, u_opt, melem, f0_loc, c0_loc, &
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
         CALL melem_rspauli_collinear(this%num_wann, melem%x0, v_ch, cell, kpts, distk, fmpi)
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


   subroutine wannierlib_create_eig(this, results, kpts, jspin, eig)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_results), INTENT(IN) :: results
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: jspin
      REAL, ALLOCATABLE, INTENT(OUT) :: eig(:, :)

      INTEGER :: nkpt
      allocate (eig(this%num_bands, kpts%nkptf))

      DO nkpt = 1, kpts%nkptf
         eig(:this%num_bands, nkpt) = results%eig(this%min_band:this%max_band, kpts%bkp(nkpt), jspin)
      END DO

   END SUBROUTINE wannierlib_create_eig

   SUBROUTINE wannierlib_kdiff(num_kpts, nntot, bk, nnkp, gkpb, kdiff)
      INTEGER, INTENT(IN) :: num_kpts
      INTEGER, INTENT(IN) :: nntot
      REAL, INTENT(IN) :: bk(:, :)
      INTEGER, INTENT(IN) :: nnkp(:, :)
      INTEGER, INTENT(IN) :: gkpb(:, :, :)
      REAL, ALLOCATABLE, INTENT(OUT) :: kdiff(:, :)

      INTEGER :: k, kk, ikpt, kd
      REAL :: kdiffvec(3)

      ALLOCATE (kdiff(3, nntot))
      kdiff = 0.0

      kd = 1
      DO k = 1, num_kpts
         k_loop: DO kk = 1, nntot
            kdiffvec = bk(:, nnkp(k, kk)) + REAL(gkpb(:, k, kk)) - bk(:, k)
            DO ikpt = 1, kd - 1
               IF (ALL(ABS(kdiff(:, ikpt) - kdiffvec) <= 1.0e-4)) CYCLE k_loop
            END DO

            IF (kd > nntot) THEN
               CALL juDFT_error("problem in wannierlib_kdiff", calledby="wannierlib_kdiff")
            end if
            kdiff(:, kd) = kdiffvec
            kd = kd + 1
         END DO k_loop
      END DO
   END SUBROUTINE wannierlib_kdiff

   ! Complete the distributed amn: each rank filled only its distk k-slice (zeros elsewhere),
   ! so an MPI_ALLREDUCE(SUM) reassembles the full amn on every rank (needed by w90_set_u_opt).
   SUBROUTINE wannierlib_reduce_amn(fmpi, amn)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      COMPLEX, INTENT(INOUT) :: amn(:, :, :)
#ifdef CPP_MPI
      INTEGER :: ierr
      IF (fmpi%isize > 1) &
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, amn, SIZE(amn), MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
#endif
   END SUBROUTINE wannierlib_reduce_amn

   subroutine wannierlib_write_mmn(this, mmnk, kpts, nnkp, gkpb, jspin, fending)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      COMPLEX, INTENT(IN) :: mmnk(:, :, :, :)
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: nnkp(:, :)
      INTEGER, INTENT(IN) :: gkpb(:, :, :)
      INTEGER, INTENT(IN), OPTIONAL :: jspin
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fending

      INTEGER :: ikpt, ikpt_b, i, j, nntot, nbnd, jspin2
      CHARACTER(LEN=3) :: spin12(2)
      CHARACTER(LEN=12) :: file_ending

      spin12 = (/'WF1', 'WF2'/)
      jspin2 = 1
      IF (PRESENT(jspin)) jspin2 = jspin
      IF (PRESENT(fending)) THEN
         file_ending = TRIM(fending)
      ELSE
         file_ending = ''
      END IF

      nntot = SIZE(mmnk, 3)
      nbnd = SIZE(mmnk, 1)

      IF (SIZE(mmnk, 2) /= nbnd) THEN
         CALL juDFT_error('wannierlib_write_mmn: mmnk must be square in band indices', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (SIZE(mmnk, 4) /= kpts%nkptf) THEN
         CALL juDFT_error('wannierlib_write_mmn: k-point dimension mismatch', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (SIZE(nnkp, 1) /= kpts%nkptf .OR. SIZE(nnkp, 2) /= nntot) THEN
         CALL juDFT_error('wannierlib_write_mmn: nnkp dimension mismatch', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (SIZE(gkpb, 1) /= 3 .OR. SIZE(gkpb, 2) /= kpts%nkptf .OR. SIZE(gkpb, 3) /= nntot) THEN
         CALL juDFT_error('wannierlib_write_mmn: gkpb dimension mismatch', &
                          calledby='wannierlib_write_mmn')
      END IF
      IF (jspin2 < 1 .OR. jspin2 > 2) THEN
         CALL juDFT_error('wannierlib_write_mmn: jspin must be 1 or 2', &
                          calledby='wannierlib_write_mmn')
      END IF

      open (305, file=spin12(jspin2)//TRIM(file_ending)//'.mmn')
      write (305, *) 'Overlaps of the wavefunct. the k- and b-points'
      write (305, '(3i5)') nbnd, kpts%nkptf, nntot
      DO ikpt = 1, kpts%nkptf
         DO ikpt_b = 1, nntot
            write (305, '(2i5,3x,3i4)') ikpt, nnkp(ikpt, ikpt_b), gkpb(1:3, ikpt, ikpt_b)
            DO i = 1, nbnd
               DO j = 1, nbnd
                  write (305, '(2f24.18)') real(mmnk(j, i, ikpt_b, ikpt)), &
                     aimag(mmnk(j, i, ikpt_b, ikpt))
               END DO
            END DO
         END DO
      END DO
      close (305)
   END SUBROUTINE wannierlib_write_mmn

END MODULE m_wannierlib_main
