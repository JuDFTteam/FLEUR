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

      INTEGER :: ikpt, itype, nntot_w90, ierr, jspin, jspin_comp
      COMPLEX, ALLOCATABLE :: amn(:, :, :)
      COMPLEX, ALLOCATABLE :: mmn(:, :, :, :)
      COMPLEX, ALLOCATABLE :: ujug(:, :, :, :, :, :)
      REAL, ALLOCATABLE :: kdiff(:, :)
      INTEGER, ALLOCATABLE :: nnkp(:, :), gkpb(:, :, :)
      INTEGER, ALLOCATABLE :: distk(:)
      real, allocatable :: eig(:, :)
      TYPE(t_usdus) :: usdus
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      TYPE(t_radfun) :: radfun(atoms%ntype)
      TYPE(t_abc) :: abc(atoms%ntype)
      LOGICAL :: l_wannierlib_spinors
      LOGICAL :: l_nocosoc
      CHARACTER(LEN=7) :: amn_file
      CHARACTER(LEN=3) :: spin12(2)
      INTEGER :: ik_local, nk_local

      IF (.NOT. this%l_wannierize) RETURN

      l_wannierlib_spinors = noco%l_noco .OR. noco%l_soc
      l_nocosoc = noco%l_noco .AND. (.NOT. noco%l_soc)
      spin12 = (/'WF1', 'WF2'/)

      !Setup of data structures for amn and mmn calculation for all k-points
      CALL usdus%init(atoms, input%jspins)
      DO itype = 1, atoms%ntype
         CALL radfun(itype)%generate_radial_functions(atoms, input, enpara, fmpi, vtot, itype, usdus=usdus)
      END DO

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
            CALL wannierlib_ujugaunt(atoms, cell, nntot_w90, kdiff, radfun, radfun, jspin_comp, jspin_comp, .FALSE., 1, ujug)

            ik_local = 0
            DO ikpt = 1, kpts%nkptf
               CALL wannierlib_get_z(this, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                     ikpt, jspin_comp, input%l_real, lapw, zMat)

               DO itype = 1, atoms%ntype
                  CALL abc(itype)%init(input, atoms, radfun(itype)%n_r, this%num_bands, itype)
                  CALL abc(itype)%calc_abc(input, atoms, sym, cell, lapw, this%num_bands, usdus, &
                                           noco, nococonv, jspin_comp, itype, zMat)
               END DO

               CALL wannierlib_amn(this, atoms, kpts, ikpt, usdus, radfun, abc, l_nocosoc, jspin_comp, amn(:, :, ikpt))

               IF (distk(ikpt) == fmpi%irank) THEN
                  ik_local = ik_local + 1
                  CALL wannierlib_mmnkb(this, this%num_bands, nntot_w90, ikpt, kpts, nnkp, gkpb, kdiff, &
                                        ujug, atoms, cell, input, sym, noco, nococonv, usdus, &
                                        radfun, abc, jspin_comp, eig_id, stars, lapw, zMat, mmn, ik_local)
               END IF
            END DO

            IF (ALLOCATED(ujug)) DEALLOCATE (ujug)
         END DO

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
         CALL run_w90(this, mmn, amn, eig)
         if (fmpi%isize == 1) CALL report_w90(this)

         IF (ALLOCATED(amn)) DEALLOCATE (amn)
         IF (ALLOCATED(mmn)) DEALLOCATE (mmn)
         IF (ALLOCATED(eig)) DEALLOCATE (eig)

      END DO

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
