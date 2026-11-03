!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_momentum
   USE m_constants, ONLY: tpi_const
   USE m_juDFT, ONLY: juDFT_error
   USE m_types_kpts, ONLY: t_kpts
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: rixs_build_kq_map
   PUBLIC :: rixs_fold_rlu
   PUBLIC :: rixs_site_phase
   PUBLIC :: rixs_phased_site_amplitude

CONTAINS

   PURE FUNCTION rixs_fold_rlu(q_full_rlu) RESULT(q_reduced_rlu)
      ! Fold a full reciprocal-lattice vector only for electronic k-point
      ! mapping. Never use this result in the absorber-site photon phase.
      REAL, INTENT(IN) :: q_full_rlu(3)
      REAL :: q_reduced_rlu(3)

      q_reduced_rlu = q_full_rlu
      WHERE (ABS(q_reduced_rlu - ANINT(q_reduced_rlu)) < 1.0e-8) q_reduced_rlu = ANINT(q_reduced_rlu)
      q_reduced_rlu = q_reduced_rlu - FLOOR(q_reduced_rlu)
   END FUNCTION rixs_fold_rlu

   PURE COMPLEX FUNCTION rixs_site_phase(q_full_rlu, tau_fractional) RESULT(phase)
      ! q_full_rlu is in the reciprocal basis dual to cell%amat and
      ! tau_fractional is a dimensionless direct-cell coordinate.  In FLEUR
      ! atoms%taual has this representation; atoms%pos is Cartesian and is not
      ! a valid argument here.
      REAL, INTENT(IN) :: q_full_rlu(3), tau_fractional(3)

      phase = EXP(CMPLX(0.0, tpi_const*DOT_PRODUCT(q_full_rlu, tau_fractional)))
   END FUNCTION rixs_site_phase

   PURE COMPLEX FUNCTION rixs_phased_site_amplitude(q_full_rlu, tau_fractional, local_amplitude) RESULT(amplitude)
      REAL,    INTENT(IN) :: q_full_rlu(3), tau_fractional(3)
      COMPLEX, INTENT(IN) :: local_amplitude

      amplitude = rixs_site_phase(q_full_rlu, tau_fractional)*local_amplitude
   END FUNCTION rixs_phased_site_amplitude

   SUBROUTINE rixs_build_kq_map(kpts, q_full_rlu, kq_index, reciprocal_shift)
      TYPE(t_kpts), INTENT(IN) :: kpts
      REAL,         INTENT(IN) :: q_full_rlu(3)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: kq_index(:), reciprocal_shift(:, :)

      INTEGER :: ikpt, jkpt
      INTEGER, ALLOCATABLE :: preimage_count(:)
      REAL :: target_unfolded(3), target_folded(3), candidate(3)
      CHARACTER(LEN=256) :: message

      IF (kpts%nkpt /= kpts%nkptf) THEN
         CALL juDFT_error("Finite-Q RIXS k+q mapping requires nkpt == nkptf.", calledby="m_rixs_momentum")
      END IF
      IF (.NOT. ALLOCATED(kpts%bk) .OR. kpts%nkpt < 1) THEN
         CALL juDFT_error("Finite-Q RIXS k+q mapping requires an allocated nonempty k-point list.", &
                          calledby="m_rixs_momentum")
      END IF

      ALLOCATE(kq_index(kpts%nkpt), reciprocal_shift(3, kpts%nkpt), preimage_count(kpts%nkpt), SOURCE=0)
      DO ikpt = 1, kpts%nkpt
         target_unfolded = kpts%bk(:, ikpt) + q_full_rlu
         target_folded = rixs_fold_rlu(target_unfolded)
         DO jkpt = 1, kpts%nkpt
            candidate = rixs_fold_rlu(kpts%bk(:, jkpt))
            IF (ALL(ABS(target_folded - candidate) < 1.0e-6)) THEN
               IF (kq_index(ikpt) /= 0) THEN
                  WRITE(message, '(a,i0,a,3f14.8)') "Finite-Q RIXS k+q map is nonunique for k-point ", ikpt, &
                                                    " and folded target ", target_folded
                  CALL juDFT_error(TRIM(message), calledby="m_rixs_momentum")
               END IF
               kq_index(ikpt) = jkpt
            END IF
         END DO
         IF (kq_index(ikpt) == 0) THEN
            WRITE(message, '(a,i0,a,3f14.8,a,3f14.8)') "Finite-Q RIXS found no k+q image for k-point ", ikpt, &
               "; full Q = ", q_full_rlu, "; folded target = ", target_folded
            CALL juDFT_error(TRIM(message), calledby="m_rixs_momentum")
         END IF
         ! Report the reciprocal lattice vector against the mesh representative
         ! that is actually used for the intermediate state, so that the
         ! diagnostic satisfies k_v + Q = k_n + reciprocal_shift for the k_v and
         ! k_n columns written next to it.  kpts%bk representatives need not lie
         ! in [0,1), so this differs in general from the shift measured against
         ! the folded target.  This value is diagnostic only.
         reciprocal_shift(:, ikpt) = NINT(target_unfolded - kpts%bk(:, kq_index(ikpt)))
         preimage_count(kq_index(ikpt)) = preimage_count(kq_index(ikpt)) + 1
      END DO

      IF (ANY(preimage_count /= 1)) THEN
         CALL juDFT_error("Finite-Q RIXS k+q map is not a one-to-one permutation of the full k mesh.", &
                          calledby="m_rixs_momentum")
      END IF
      DO ikpt = 1, kpts%nkpt
         IF (ABS(kpts%wtkpt(ikpt) - kpts%wtkpt(kq_index(ikpt))) > 1.0e-10) THEN
            CALL juDFT_error("Finite-Q RIXS k+q permutation does not preserve full-mesh k-point weights.", &
                             calledby="m_rixs_momentum")
         END IF
      END DO
   END SUBROUTINE rixs_build_kq_map

END MODULE m_rixs_momentum
