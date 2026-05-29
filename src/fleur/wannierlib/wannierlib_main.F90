!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
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
   IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_main(this, atoms, cell, input, kpts, sym, noco, nococonv, stars, enpara, fmpi, vtot, results,  l_nocosoc, jspin, eig_id)
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
      LOGICAL, INTENT(IN) :: l_nocosoc
      INTEGER, INTENT(IN) :: jspin
      INTEGER, INTENT(IN) :: eig_id

      INTEGER :: ikpt, itype, nntot_w90
      COMPLEX, ALLOCATABLE :: amn(:, :, :)
      COMPLEX, ALLOCATABLE :: mmn(:, :, :, :)
      COMPLEX, ALLOCATABLE :: ujug(:, :, :, :, :, :)
      REAL, ALLOCATABLE :: kdiff(:, :)
      INTEGER, ALLOCATABLE :: nnkp(:, :), gkpb(:, :, :)
      real, allocatable :: eig(:, :)
      TYPE(t_usdus) :: usdus
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat
      TYPE(t_radfun) :: radfun(atoms%ntype)
      TYPE(t_abc) :: abc(atoms%ntype)

      IF (.NOT. this%l_wannierize) RETURN

      CALL init_w90(this, kpts, nntot_w90, nnkp, gkpb)

      !Setup of data structures for amn and mmn calculation for all k-points
      CALL usdus%init(atoms, input%jspins)
      DO itype = 1, atoms%ntype
         CALL radfun(itype)%generate_radial_functions(atoms, input, enpara, fmpi, vtot, itype, usdus=usdus)
      END DO

      CALL wannierlib_kdiff(kpts%nkpt, nntot_w90, kpts%bk, nnkp, gkpb, kdiff)
      CALL wannierlib_ujugaunt(atoms, cell, nntot_w90, kdiff, radfun, radfun, jspin, jspin, .FALSE., 1, ujug)

      ! calculate the  matrices for all k-points
      ALLOCATE (amn(this%num_bands, this%num_wann, kpts%nkptf), source=cmplx(0.0, 0.0))
      allocate (mmn(this%num_bands, this%num_bands, nntot_w90, kpts%nkptf), source=cmplx(0.0, 0.0))

      DO ikpt = 1, kpts%nkptf
         CALL wannierlib_get_z(this, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, ikpt, jspin, input%l_real, lapw, zMat)

         DO itype = 1, atoms%ntype
            CALL abc(itype)%init(input, atoms, radfun(itype)%n_r, this%num_bands, itype)
            CALL abc(itype)%calc_abc(input, atoms, sym, cell, lapw, this%num_bands, usdus, noco, nococonv, jspin, itype, zMat)
         END DO

         CALL wannierlib_amn(this, atoms, cell, input, kpts, ikpt, usdus, radfun, abc, l_nocosoc, jspin, amn(:, :, ikpt))
      CALL wannierlib_mmnkb(this, this%num_bands, nntot_w90, ikpt, kpts, nnkp, gkpb, kdiff, ujug, atoms, cell, input, sym, noco, nococonv, usdus, &
                               radfun, abc, jspin, eig_id, stars, lapw, zMat, mmn)
      END DO

      call wannierlib_create_eig(this, results, kpts, jspin, eig)
      CALL run_w90(this, mmn, amn, eig)
      CALL report_w90(this)
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
            kdiffvec = bk(:, nnkp(kk, k)) + REAL(gkpb(:, kk, k)) - bk(:, k)
            DO ikpt = 1, kd - 1
               IF (ALL(ABS(kdiff(:, ikpt) - kdiffvec) <= 1.0e-4)) CYCLE k_loop
            END DO

            IF (kd > nntot) CALL juDFT_error("problem in wannierlib_kdiff", calledby="wannierlib_kdiff")
            kdiff(:, kd) = kdiffvec
            kd = kd + 1
         END DO k_loop
      END DO
   END SUBROUTINE wannierlib_kdiff

END MODULE m_wannierlib_main
