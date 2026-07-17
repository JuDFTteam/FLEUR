!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_abcoeff_store
   !! Optional storage for the LAPW matching coefficients (abCoeffs) computed in
   !! hsmt_ab. The coefficients are kept in "slots" indexed by
   !! (k-point, igSpin, ilSpin, atom). On a later call for the same indices the
   !! stored array can be moved back (via move_alloc, no copy) so that hsmt_ab can
   !! skip the recomputation entirely.
   !!
   !! Usage by a driver:
   !!   1. CALL abcoeff_store_init(input,noco,kpts,atoms)   ! prepare slots, switch storage on
   !!   2. ... hsmt_ab uses abcoeff_store_alloc to allocate-or-retrieve abCoeffs ...
   !!   3. CALL abcoeff_store_save(abCoeffs,nk,igSpin,ilSpin,na)  ! move data back after usage
   !!   4. CALL abcoeff_store_free()                         ! empty the slots when data is stale
   !!
   !! NOTE: the cached coefficients depend on fjgj/lapw/geometry, which change
   !! between SCF iterations. It is the caller's responsibility to free (or re-init)
   !! the storage whenever those inputs change. The storage is global module state
   !! and is not safe for concurrent hsmt_ab calls from different threads.

   USE m_juDFT

   IMPLICIT NONE
   PRIVATE

   !> Logical that determines whether the storage is used by hsmt_ab.
   LOGICAL, PUBLIC, SAVE :: l_use_abcoeff_store = .FALSE.

   !> One slot holds the matching coefficients of a single (k,igSpin,ilSpin,atom).
   TYPE :: t_abcoeff_slot
      COMPLEX, ALLOCATABLE :: ab(:, :)
   END TYPE t_abcoeff_slot

   !> Slots indexed by (nk, igSpin, ilSpin, na).
   TYPE(t_abcoeff_slot), ALLOCATABLE, SAVE :: storage(:, :, :, :)

   PUBLIC :: abcoeff_store_init
   PUBLIC :: abcoeff_store_alloc
   PUBLIC :: abcoeff_store_save
   PUBLIC :: abcoeff_store_free

CONTAINS

   SUBROUTINE abcoeff_store_init(input, noco, kpts, atoms)
      !! Prepare the storage slots and switch the storage on. The spin dimensions
      !! are sized to the actual spin structure of the calculation (mirroring the
      !! loop bounds in hsmt and the c_ph dimension in hsmt_ab), not hardwired to 2.
      USE m_types_input
      USE m_types_noco
      USE m_types_kpts
      USE m_types_atoms

      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_noco),  INTENT(IN) :: noco
      TYPE(t_kpts),  INTENT(IN) :: kpts
      TYPE(t_atoms), INTENT(IN) :: atoms

      INTEGER :: n_igSpin, n_ilSpin, ierr

      ! Only enable the storage if the user requested it on the command line
      ! (-abcoeff_store, registered in add_fleur_arguments). Return without
      ! allocating any slots or switching the storage on otherwise.
      IF (.NOT.juDFT_was_argument("-abcoeff_store")) RETURN

      ! local spin in the atom: 1..2 for noco, otherwise 1..jspins (see hsmt.F90)
      n_ilSpin = MERGE(2, input%jspins, noco%l_noco)
      ! global/interstitial spin: doubles only when both spin sets are needed
      ! (same condition as the c_ph dimension in hsmt_ab)
      n_igSpin = MERGE(2, 1, noco%l_ss .OR. ANY(noco%l_unrestrictMT) .OR. ANY(noco%l_spinoffd_ldau))

      ! Re-init is allowed: drop any previous slots first.
      IF (ALLOCATED(storage)) CALL abcoeff_store_free()

      ALLOCATE (storage(kpts%nkpt, n_igSpin, n_ilSpin, atoms%nat), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error("Couldn't allocate abCoeff storage")

      l_use_abcoeff_store = .TRUE.
   END SUBROUTINE abcoeff_store_init

   LOGICAL FUNCTION abcoeff_store_alloc(abCoeffs, dim1, dim2, nk, igSpin, ilSpin, na) RESULT(l_reused)
      !! Allocate-or-retrieve the matching coefficient array. Replaces the bare
      !! ALLOCATE in hsmt_ab.
      !!
      !! Returns .TRUE.  if stored coefficients with matching dimensions were moved
      !!                 into abCoeffs -> the caller may skip the computation.
      !! Returns .FALSE. if only a fresh array was allocated -> the caller must
      !!                 (re)compute the coefficients.
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: abCoeffs(:, :)
      INTEGER,              INTENT(IN)  :: dim1, dim2, nk, igSpin, ilSpin, na

      INTEGER :: ierr

      l_reused = .FALSE.

      IF (l_use_abcoeff_store .AND. ALLOCATED(storage)) THEN
         IF (nk     >= 1 .AND. nk     <= SIZE(storage, 1) .AND. &
             igSpin >= 1 .AND. igSpin <= SIZE(storage, 2) .AND. &
             ilSpin >= 1 .AND. ilSpin <= SIZE(storage, 3) .AND. &
             na     >= 1 .AND. na     <= SIZE(storage, 4)) THEN
            IF (ALLOCATED(storage(nk, igSpin, ilSpin, na)%ab)) THEN
               IF (SIZE(storage(nk, igSpin, ilSpin, na)%ab, 1) == dim1 .AND. &
                   SIZE(storage(nk, igSpin, ilSpin, na)%ab, 2) == dim2) THEN
                  ! Correctly dimensioned data is available: hand it back.
                  CALL move_alloc(storage(nk, igSpin, ilSpin, na)%ab, abCoeffs)
                  l_reused = .TRUE.
                  RETURN
               ELSE
                  ! Stored data has the wrong dimensions: discard it.
                  DEALLOCATE (storage(nk, igSpin, ilSpin, na)%ab)
               END IF
            END IF
         END IF
      END IF

      ! No (valid) stored data: allocate a fresh array.
      ALLOCATE (abCoeffs(dim1, dim2), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error("Couldn't allocate abCoeffs")
   END FUNCTION abcoeff_store_alloc

   SUBROUTINE abcoeff_store_save(abCoeffs, nk, igSpin, ilSpin, na)
      !! Move abCoeffs back into its slot for later reuse. No-op (and abCoeffs is
      !! left untouched) when the storage is off, not initialised, or the indices
      !! are out of bounds.
      COMPLEX, ALLOCATABLE, INTENT(INOUT) :: abCoeffs(:, :)
      INTEGER,              INTENT(IN)    :: nk, igSpin, ilSpin, na

      IF (.NOT. (l_use_abcoeff_store .AND. ALLOCATED(storage))) RETURN
      IF (.NOT. ALLOCATED(abCoeffs)) RETURN
      IF (nk     >= 1 .AND. nk     <= SIZE(storage, 1) .AND. &
          igSpin >= 1 .AND. igSpin <= SIZE(storage, 2) .AND. &
          ilSpin >= 1 .AND. ilSpin <= SIZE(storage, 3) .AND. &
          na     >= 1 .AND. na     <= SIZE(storage, 4)) THEN
         CALL move_alloc(abCoeffs, storage(nk, igSpin, ilSpin, na)%ab)
      END IF
   END SUBROUTINE abcoeff_store_save

   SUBROUTINE abcoeff_store_free()
      !! Empty the storage: deallocate each slot's allocatable %ab component. The
      !! slot array itself and the l_use_abcoeff_store flag are intentionally left
      !! unchanged.
      INTEGER :: i, j, k, l
      IF (ALLOCATED(storage)) THEN
          DO i = 1, SIZE(storage, 1)
              DO j = 1, SIZE(storage, 2)
                  DO k = 1, SIZE(storage, 3)
                      DO l = 1, SIZE(storage, 4)
                          IF (ALLOCATED(storage(i, j, k, l)%ab)) THEN
                              DEALLOCATE (storage(i, j, k, l)%ab)
                          END IF
                      END DO
                  END DO
              END DO
          END DO
      END IF
   END SUBROUTINE abcoeff_store_free

END MODULE m_abcoeff_store
