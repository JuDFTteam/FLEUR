!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_main
   USE m_types, ONLY: t_stars, t_results
   USE m_wannierlib_get_z
   USE m_wannierlib_amn
   USE m_wannierlib_mmnkb
   USE m_wannierlib_ujugaunt
   USE m_wannierlib_w90_adapter
   USE m_wannierlib_spin_melem, ONLY: wannierlib_spin_peratom, wannierlib_spin_bloch
   USE m_wannierlib_orbmom_melem
   USE m_wannierlib_socmat_melem
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
   USE m_types_abc
   USE m_types_wannierlib

   use m_wann_write_amn
   IMPLICIT NONE
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

      INTEGER :: ikpt, itype, nntot_w90, ierr, jspin, jspin_comp, iop
      COMPLEX, ALLOCATABLE :: amn(:, :, :)
      COMPLEX, ALLOCATABLE :: mmn(:, :, :, :)
      COMPLEX, ALLOCATABLE :: mmn_full(:, :, :, :)   ! (num_bands,num_bands,nntot,nkptf) full overlaps on rank 0 (Berry/interband velocity)
      LOGICAL :: l_need_mmn_full                     ! velocity/current requested -> gather mmn to rank 0
      LOGICAL :: l_need_gather                        ! spin/orbital current requested -> gather full coarse arrays to rank 0
      COMPLEX, ALLOCATABLE :: ujug(:, :, :, :, :, :)
      REAL, ALLOCATABLE :: kdiff(:, :)
      INTEGER, ALLOCATABLE :: nnkp(:, :), gkpb(:, :, :)
      INTEGER, ALLOCATABLE :: distk(:)
      real, allocatable :: eig(:, :)
      COMPLEX, ALLOCATABLE :: s0_coarse(:, :, :, :)   ! (num_bands,num_bands,3,nkptf) Bloch spin, rank 0
      COMPLEX, ALLOCATABLE :: l0_coarse(:, :, :, :, :) ! (num_bands,num_bands,3,nat,nkptf) Bloch L per atom, rank 0
      COMPLEX, ALLOCATABLE :: soc0_coarse(:, :, :, :) ! (num_bands,num_bands,1,nkptf) Bloch SOC, rank 0
      COMPLEX, ALLOCATABLE :: soc4_coarse(:, :, :, :) ! (num_bands,num_bands,4,nkptf) 2x2 SOC blocks (rssocmat.1)
      COMPLEX, ALLOCATABLE :: s0pa_coarse(:, :, :, :, :) ! (num_bands,num_bands,3,nat,nkptf) per-atom MT spin, rank 0
      COMPLEX, ALLOCATABLE :: s0_loc(:,:,:,:), l0_loc(:,:,:,:,:), soc0_loc(:,:,:,:), soc4_loc(:,:,:,:), s0pa_loc(:,:,:,:,:) ! per-rank slices
      INTEGER :: nkc_loc   ! # coarse k-points owned by this rank
      TYPE(t_usdus) :: usdus
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      TYPE(t_radfun) :: radfun(atoms%ntype)
      TYPE(t_abc) :: abc(atoms%ntype)
      LOGICAL :: l_wannierlib_spinors
      LOGICAL :: l_nocosoc
      LOGICAL :: l_real_wann
      INTEGER :: jspin_rad
      CHARACTER(LEN=7) :: amn_file
      CHARACTER(LEN=3) :: spin12(2)
      INTEGER :: ik_local, nk_local
      CHARACTER(LEN=6) :: spin_sfx

      IF (.NOT. this%l_wannierize) RETURN

      l_wannierlib_spinors = noco%l_noco .OR. noco%l_soc
      l_nocosoc = noco%l_noco .AND. (.NOT. noco%l_soc)
      ! A.1: input%l_real queda TRUE con inversion aunque haya SOC (n_denmat=0);
      ! leer un espinor complejo en buffer real mata la parte imaginaria del MMN.
      l_real_wann = input%l_real .AND. .NOT. noco%l_soc
      spin12 = (/'WF1', 'WF2'/)

      !Setup of data structures for amn and mmn calculation for all k-points
      CALL usdus%init(atoms, input%jspins)
      DO itype = 1, atoms%ntype
         CALL radfun(itype)%generate_radial_functions(atoms, input, enpara, fmpi, vtot, itype, usdus=usdus)
      END DO

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

      ! Operator Bloch matrices on the coarse mesh: k-loop DISTRIBUTED over ranks (each its
      ! distk slice -> parallel get_z I/O), into per-rank local arrays, then MPI_Gatherv to
      ! full-on-rank-0 (interpolation consumes the full set). Byte-identical to the serial path.
      IF ((this%l_spin .OR. this%l_orbmom .OR. this%l_socop) .AND. l_wannierlib_spinors) THEN
         nkc_loc = MAX(1, COUNT(distk == fmpi%irank))
         ALLOCATE(s0_loc(this%num_bands, this%num_bands, 3, nkc_loc), source=cmplx(0.0,0.0))
         ALLOCATE(l0_loc(this%num_bands, this%num_bands, 3, atoms%nat, nkc_loc), source=cmplx(0.0,0.0))
         ALLOCATE(soc0_loc(this%num_bands, this%num_bands, 1, nkc_loc), source=cmplx(0.0,0.0))
         ALLOCATE(soc4_loc(this%num_bands, this%num_bands, 4, nkc_loc), source=cmplx(0.0,0.0))
         ALLOCATE(s0pa_loc(this%num_bands, this%num_bands, 3, atoms%nat, nkc_loc), source=cmplx(0.0,0.0))
         CALL wannierlib_operator_coarse(this, atoms, input, sym, cell, noco, nococonv, kpts, &
                                         stars, usdus, radfun, enpara, fmpi, vtot, eig_id, l_real_wann, &
                                         distk, s0_loc, l0_loc, soc0_loc, soc4_loc, s0pa_loc)
         ! Fase 3a: the interpolation operators (spin/orbital/soc, total & per-atom) now consume
         ! the per-rank LOCAL coarse slices + a distributed FT-reduce inside run_w90, so they no
         ! longer need the full-mesh arrays materialized on rank 0. Only the spin/orbital CURRENTS
         ! still read the full s0_coarse/l0_coarse (Fase 3b will migrate them) -> gather on demand.
         l_need_gather = .FALSE.
         DO iop = 1, this%n_ops
            SELECT CASE (TRIM(this%op_name(iop)))
            CASE ('spinCurrent', 'orbitalCurrent'); l_need_gather = .TRUE.
            END SELECT
         END DO
         IF (l_need_gather) THEN
            IF (fmpi%irank == 0) THEN
               ALLOCATE(s0_coarse(this%num_bands, this%num_bands, 3, kpts%nkptf))
               ALLOCATE(l0_coarse(this%num_bands, this%num_bands, 3, atoms%nat, kpts%nkptf))
            ELSE
               ALLOCATE(s0_coarse(1,1,1,1)); ALLOCATE(l0_coarse(1,1,1,1,1))
            END IF
            CALL wannierlib_gather_coarse(fmpi, distk, kpts%nkptf, this%num_bands*this%num_bands*3, s0_loc, s0_coarse)
            CALL wannierlib_gather_coarse(fmpi, distk, kpts%nkptf, this%num_bands*this%num_bands*3*atoms%nat, l0_loc, l0_coarse)
         ELSE
            ALLOCATE(s0_coarse(1,1,1,1)); ALLOCATE(l0_coarse(1,1,1,1,1))
         END IF
         ! soc0/soc4/s0pa full-mesh arrays are no longer consumed (interp uses local slices) -> stubs
         ALLOCATE(soc0_coarse(1,1,1,1)); ALLOCATE(soc4_coarse(1,1,1,1)); ALLOCATE(s0pa_coarse(1,1,1,1,1))
         ! keep ALL per-rank locals (s0_loc/l0_loc/soc0_loc/soc4_loc/s0pa_loc) alive for the run_w90 reduce
      ELSE
         ALLOCATE(s0_coarse(1, 1, 1, 1)); ALLOCATE(l0_coarse(1, 1, 1, 1, 1)); ALLOCATE(soc0_coarse(1, 1, 1, 1)); ALLOCATE(soc4_coarse(1, 1, 1, 1))
         ALLOCATE(s0pa_coarse(1, 1, 1, 1, 1))
         ALLOCATE(s0_loc(1,1,1,1)); ALLOCATE(l0_loc(1,1,1,1,1)); ALLOCATE(soc4_loc(1,1,1,1))
         ALLOCATE(soc0_loc(1,1,1,1)); ALLOCATE(s0pa_loc(1,1,1,1,1))
      END IF

      CALL init_w90(this, atoms, cell, kpts, fmpi, l_wannierlib_spinors, nntot_w90, nnkp, gkpb, distk)
      CALL wannierlib_kdiff(kpts%nkptf, nntot_w90, kpts%bkf, nnkp, gkpb, kdiff)

      DO jspin = 1, MERGE(2, 1, input%jspins == 2 .AND. (.NOT. l_wannierlib_spinors))

         ! calculate the  matrices for all k-points
         ALLOCATE (amn(this%num_bands, this%num_wann, kpts%nkptf), stat=ierr, source=cmplx(0.0, 0.0))
         IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating amn buffer', calledby='wannierlib_main')

         nk_local = COUNT(distk == fmpi%irank)
         ALLOCATE(mmn(this%num_bands, this%num_bands, nntot_w90, nk_local), stat=ierr, source=cmplx(0.0, 0.0))
         IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating local mmn buffer', calledby='wannierlib_main')

         DO jspin_comp = MERGE(1, jspin, l_wannierlib_spinors), MERGE(2, jspin, l_wannierlib_spinors)
            ! jspin_comp = record del eig (spinor up/down). jspin_rad = indice radial:
            ! con jspins=1 solo existe 1 set de radiales -> usar 1 para ambas componentes.
            jspin_rad = MERGE(1, jspin_comp, input%jspins == 1)
            CALL wannierlib_ujugaunt(atoms, cell, nntot_w90, kdiff, radfun, radfun, jspin_rad, jspin_rad, .FALSE., 1, ujug)

            ik_local = 0
            DO ikpt = 1, kpts%nkptf
               IF (distk(ikpt) /= fmpi%irank) CYCLE   ! each rank computes only its k-slice -> parallel eigenvector I/O
               CALL wannierlib_get_z(this, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                     ikpt, jspin_comp, l_real_wann, lapw, zMat)

               DO itype = 1, atoms%ntype
                  CALL abc(itype)%init(input, atoms, radfun(itype)%n_r, this%num_bands, itype)
                  CALL abc(itype)%calc_abc(input, atoms, sym, cell, lapw, this%num_bands, usdus, &
                                           noco, nococonv, jspin_rad, itype, zMat)
               END DO

               CALL wannierlib_amn(this, atoms, kpts, ikpt, usdus, radfun, abc, l_nocosoc, jspin_comp, jspin_rad, amn(:, :, ikpt))

               ik_local = ik_local + 1
               CALL wannierlib_mmnkb(this, this%num_bands, nntot_w90, ikpt, kpts, nnkp, gkpb, kdiff, &
                                     ujug, atoms, cell, input, sym, noco, nococonv, usdus, &
                                     radfun, abc, jspin_comp, jspin_rad, eig_id, stars, lapw, zMat, mmn, ik_local)
            END DO

            IF (ALLOCATED(ujug)) DEALLOCATE (ujug)
         END DO

         ! amn was filled only on each rank's distk slice (zeros elsewhere) -> sum to the full set
         CALL wannierlib_reduce_amn(fmpi, amn)

         mmn = conjg(mmn)

         ! --- gather the distributed (per-rank) mmn into a full-nkptf buffer on rank 0,
         !     needed by the Berry connection for the interband velocity / currents. ---
         l_need_mmn_full = .FALSE.
         DO iop = 1, this%n_ops
            SELECT CASE (TRIM(this%op_name(iop)))
            CASE ('velocity', 'spinCurrent', 'orbitalCurrent', 'position_r'); l_need_mmn_full = .TRUE.
            END SELECT
         END DO
         DO iop = 1, this%n_op_r
            IF (TRIM(this%op_r_name(iop)) == 'position') l_need_mmn_full = .TRUE.
         END DO
         IF (l_need_mmn_full) THEN
            CALL wannierlib_gather_mmn(fmpi, distk, kpts%nkptf, this%num_bands, nntot_w90, mmn, mmn_full)
            ! validation hook: dump the gathered overlaps (rank 0) so a parallel run can be
            ! compared byte-for-byte against a serial reference .mmn (WF*_gathered vs WF*).
            IF (fmpi%irank == 0) CALL wannierlib_write_mmn(this, mmn_full, kpts, nnkp, gkpb, jspin, '_gathered')
         ELSE
            ALLOCATE(mmn_full(1, 1, 1, 1))
         END IF

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
         CALL run_w90(this, cell, kpts, mmn, amn, eig, fmpi%irank, s0_coarse, l0_coarse, soc0_coarse, soc4_coarse, s0pa_coarse, mmn_full, &
                      s0_loc, l0_loc, soc0_loc, soc4_loc, s0pa_loc, distk, fmpi%mpi_comm, &
                      spin_suffix=TRIM(spin_sfx))
         if (fmpi%isize == 1) CALL report_w90(this)

         IF (ALLOCATED(amn)) DEALLOCATE (amn)
         IF (ALLOCATED(mmn)) DEALLOCATE (mmn)
         IF (ALLOCATED(mmn_full)) DEALLOCATE (mmn_full)
         IF (ALLOCATED(eig)) DEALLOCATE (eig)

      END DO

      IF (ALLOCATED(distk)) DEALLOCATE(distk)
      IF (ALLOCATED(s0_coarse)) DEALLOCATE(s0_coarse)
      IF (ALLOCATED(l0_coarse)) DEALLOCATE(l0_coarse)
      IF (ALLOCATED(soc0_coarse)) DEALLOCATE(soc0_coarse)
      IF (ALLOCATED(soc4_coarse)) DEALLOCATE(soc4_coarse)
      IF (ALLOCATED(s0_loc)) DEALLOCATE(s0_loc)
      IF (ALLOCATED(l0_loc)) DEALLOCATE(l0_loc)
      IF (ALLOCATED(soc0_loc)) DEALLOCATE(soc0_loc)
      IF (ALLOCATED(soc4_loc)) DEALLOCATE(soc4_loc)
      IF (ALLOCATED(s0pa_loc)) DEALLOCATE(s0pa_loc)
      IF (ALLOCATED(s0pa_coarse)) DEALLOCATE(s0pa_coarse)

   END SUBROUTINE wannierlib_main

   !> Build the requested Bloch-basis operator matrices on the FULL coarse mesh
   !> (rank 0): spin S0(k) (if l_spin) and/or orbital L0(k) (if l_orbmom), sharing
   !> the (get_z + calc_abc) work in one k-pass. Needs only the ab-initio spinor
   !> eigenstates (no U), so it runs before the wannierization.
   SUBROUTINE wannierlib_operator_coarse(this, atoms, input, sym, cell, noco, nococonv, kpts, &
                                         stars, usdus, radfun, enpara, fmpi, vtot, eig_id, l_real_wann, &
                                         distk, s0_coarse, l0_coarse, soc0_coarse, soc4_coarse, s0pa_coarse)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
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
      INTEGER, INTENT(IN) :: distk(:)   ! rank owner of each global k (distribute the loop)
      COMPLEX, INTENT(OUT) :: s0_coarse(:, :, :, :)   ! (num_bands,num_bands,3,nkptf) spin
      COMPLEX, INTENT(OUT) :: l0_coarse(:, :, :, :, :) ! (num_bands,num_bands,3,nat,nkptf) orbital L per atom
      COMPLEX, INTENT(OUT) :: soc0_coarse(:, :, :, :) ! (num_bands,num_bands,1,nkptf) SOC
      COMPLEX, INTENT(OUT) :: soc4_coarse(:, :, :, :) ! (num_bands,num_bands,4,nkptf) 2x2 SOC blocks
      COMPLEX, INTENT(OUT) :: s0pa_coarse(:, :, :, :, :) ! (num_bands,num_bands,3,nat,nkptf) per-atom MT spin

      TYPE(t_abc), ALLOCATABLE :: abc_s(:, :)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      INTEGER :: ikpt, itype, isp, il

      ALLOCATE(abc_s(atoms%ntype, 2))
      il = 0
      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE   ! this rank only computes its own k-slice
         il = il + 1                            ! local (ascending global-k) index for gather
         ! noco: get_z returns the FULL spinor; both local components via calc_abc(jspin=1,2)
         CALL wannierlib_get_z(this, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                               ikpt, 1, l_real_wann, lapw, zMat)
         DO isp = 1, 2
            DO itype = 1, atoms%ntype
               CALL abc_s(itype, isp)%init(input, atoms, radfun(itype)%n_r, this%num_bands, itype)
               CALL abc_s(itype, isp)%calc_abc(input, atoms, sym, cell, lapw, this%num_bands, usdus, &
                                               noco, nococonv, isp, itype, zMat)
            END DO
         END DO
         IF (this%l_spin) CALL wannierlib_spin_peratom(atoms, abc_s, radfun, nococonv, s0pa_coarse(:, :, :, :, il))
         IF (this%l_spin) CALL wannierlib_spin_bloch(atoms, abc_s, radfun, nococonv, stars, lapw, zMat, &
                                    this%num_bands, ikpt, s0_coarse(:, :, :, il), ikpt <= 3)
         IF (this%l_orbmom) CALL wannierlib_orbmom_bloch(atoms, abc_s, radfun, l0_coarse(:, :, :, :, il))
         IF (this%l_socop) CALL wannierlib_socmat_bloch(atoms, noco, nococonv, input, fmpi, enpara, vtot, &
                                    usdus, abc_s, this%num_bands, soc0_coarse(:, :, :, il), soc4_coarse(:, :, :, il))
      END DO

      DEALLOCATE(abc_s)
   END SUBROUTINE wannierlib_operator_coarse

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

   ! Gather the k-distributed overlap matrix mmn (each rank owns COUNT(distk==irank)
   ! k-points, in ascending global-k order) into a full-nkptf buffer on rank 0. Robust
   ! for ANY distk map: MPI_Gatherv of per-rank chunks + reorder to global-k order.
   ! mmn_full is (nb,nb,nntot,nkptf) on rank 0 and (1,1,1,1) elsewhere. Serial: a copy.
   SUBROUTINE wannierlib_gather_mmn(fmpi, distk, nkptf, nb, nntot, mmn_loc, mmn_full)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:), nkptf, nb, nntot
      COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)                 ! (nb,nb,nntot,nk_local)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: mmn_full(:, :, :, :)
      INTEGER :: bs, r, ik, ierr, nk_local
#ifdef CPP_MPI
      INTEGER, ALLOCATABLE :: recvcnt(:), displ(:), locidx(:)
      COMPLEX, ALLOCATABLE :: recvbuf(:)
#endif
      bs = nb*nb*nntot
      nk_local = SIZE(mmn_loc, 4)

      IF (fmpi%isize == 1) THEN                     ! serial: local IS the full set
         ALLOCATE(mmn_full(nb, nb, nntot, nkptf))
         mmn_full = mmn_loc
         RETURN
      END IF

#ifdef CPP_MPI
      ALLOCATE(recvcnt(0:fmpi%isize-1), displ(0:fmpi%isize-1))
      DO r = 0, fmpi%isize-1
         recvcnt(r) = bs * COUNT(distk == r)
      END DO
      displ(0) = 0
      DO r = 1, fmpi%isize-1
         displ(r) = displ(r-1) + recvcnt(r-1)
      END DO
      IF (fmpi%irank == 0) THEN
         ALLOCATE(recvbuf(bs*nkptf))
      ELSE
         ALLOCATE(recvbuf(1))
      END IF
      CALL MPI_Gatherv(mmn_loc, bs*nk_local, MPI_DOUBLE_COMPLEX, &
                       recvbuf, recvcnt, displ, MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr)
      IF (fmpi%irank == 0) THEN
         ALLOCATE(mmn_full(nb, nb, nntot, nkptf))
         ALLOCATE(locidx(0:fmpi%isize-1)); locidx = 0
         DO ik = 1, nkptf                           ! place each block at its global-k slot
            r = distk(ik)
            mmn_full(:, :, :, ik) = RESHAPE(recvbuf(displ(r)+locidx(r)*bs+1 : displ(r)+locidx(r)*bs+bs), &
                                            (/ nb, nb, nntot /))
            locidx(r) = locidx(r) + 1
         END DO
         DEALLOCATE(locidx)
      ELSE
         ALLOCATE(mmn_full(1, 1, 1, 1))
      END IF
      DEALLOCATE(recvcnt, displ, recvbuf)
#else
      ALLOCATE(mmn_full(nb, nb, nntot, nkptf))
      mmn_full = mmn_loc
#endif
   END SUBROUTINE wannierlib_gather_mmn

   ! Gather a per-rank coarse-operator slice arr_loc(bs, nk_local) (ascending global-k order)
   ! into arr_full(bs, nkptf) on rank 0. bs = product of all leading dims (nb*nb*ncomp[*nat]).
   ! Mirrors wannierlib_gather_mmn but generic in bs; arrays passed by sequence association.
   SUBROUTINE wannierlib_gather_coarse(fmpi, distk, nkptf, bs, arr_loc, arr_full)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:), nkptf, bs
      COMPLEX, INTENT(IN)  :: arr_loc(bs, *)    ! (bs, nk_local)
      COMPLEX, INTENT(OUT) :: arr_full(bs, *)   ! (bs, nkptf) on rank 0; stub elsewhere (untouched)
      INTEGER :: r, ik, ierr, nk_local
#ifdef CPP_MPI
      INTEGER, ALLOCATABLE :: recvcnt(:), displ(:), locidx(:)
      COMPLEX, ALLOCATABLE :: recvbuf(:)
#endif
      nk_local = COUNT(distk == fmpi%irank)
      IF (fmpi%isize == 1) THEN
         arr_full(:, 1:nkptf) = arr_loc(:, 1:nkptf)
         RETURN
      END IF
#ifdef CPP_MPI
      ALLOCATE(recvcnt(0:fmpi%isize-1), displ(0:fmpi%isize-1))
      DO r = 0, fmpi%isize-1
         recvcnt(r) = bs * COUNT(distk == r)
      END DO
      displ(0) = 0
      DO r = 1, fmpi%isize-1
         displ(r) = displ(r-1) + recvcnt(r-1)
      END DO
      IF (fmpi%irank == 0) THEN
         ALLOCATE(recvbuf(bs*nkptf))
      ELSE
         ALLOCATE(recvbuf(1))
      END IF
      CALL MPI_Gatherv(arr_loc, bs*nk_local, MPI_DOUBLE_COMPLEX, &
                       recvbuf, recvcnt, displ, MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr)
      IF (fmpi%irank == 0) THEN
         ALLOCATE(locidx(0:fmpi%isize-1)); locidx = 0
         DO ik = 1, nkptf
            r = distk(ik)
            arr_full(:, ik) = recvbuf(displ(r)+locidx(r)*bs+1 : displ(r)+locidx(r)*bs+bs)
            locidx(r) = locidx(r) + 1
         END DO
         DEALLOCATE(locidx)
      END IF
      DEALLOCATE(recvbuf, recvcnt, displ)
#endif
   END SUBROUTINE wannierlib_gather_coarse

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
