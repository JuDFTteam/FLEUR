!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_d_spin_frame_transform
   USE m_juDFT, ONLY: juDFT_error
   USE m_local_d_spin_analysis, ONLY: local_d_complex_to_real_unitary, &
                                      local_d_orb_xy, local_d_orb_yz, local_d_orb_xz, &
                                      local_d_orb_x2y2, local_d_orb_z2
   USE m_local_d_spin_projector_core, ONLY: local_d_l, local_d_n_orbitals, local_d_n_spins
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: local_d_orbital_global_to_local
   PUBLIC :: spinor_local_to_global
   PUBLIC :: local_d_spin_native_to_structural
   PUBLIC :: transform_local_d_spin_density
   PUBLIC :: transform_local_d_spin_density_with_matrix

CONTAINS

   SUBROUTINE local_d_orbital_global_to_local(local_to_global, orbital_global_to_local)
      ! The real l=2 basis is represented by orthonormal symmetric-traceless
      ! Cartesian tensors Q_a. For x_global=R*x_local, the same scalar orbital
      ! has Q_local=R^T Q_global R. The returned complex matrix maps spherical-
      ! harmonic coefficients as c_local=D_global_to_local*c_global.
      REAL, INTENT(IN) :: local_to_global(3, 3)
      COMPLEX, INTENT(OUT) :: orbital_global_to_local(-local_d_l:local_d_l, &
                                                       -local_d_l:local_d_l)

      REAL :: tensors(3, 3, local_d_n_orbitals)
      REAL :: rotated_tensor(3, 3), real_transform(local_d_n_orbitals, local_d_n_orbitals)
      COMPLEX :: complex_to_real(-local_d_l:local_d_l, local_d_n_orbitals)
      INTEGER :: local_orbital, global_orbital

      CALL validate_rotation(local_to_global, "structural local-to-global rotation")
      CALL initialize_real_d_tensors(tensors)
      DO global_orbital = 1, local_d_n_orbitals
         rotated_tensor = MATMUL(TRANSPOSE(local_to_global), &
                                 MATMUL(tensors(:, :, global_orbital), local_to_global))
         DO local_orbital = 1, local_d_n_orbitals
            real_transform(local_orbital, global_orbital) = &
               SUM(tensors(:, :, local_orbital)*rotated_tensor)
         END DO
      END DO
      CALL local_d_complex_to_real_unitary(complex_to_real)
      orbital_global_to_local = MATMUL(complex_to_real, &
         MATMUL(CMPLX(real_transform, 0.0), CONJG(TRANSPOSE(complex_to_real))))
   END SUBROUTINE local_d_orbital_global_to_local

   SUBROUTINE spinor_local_to_global(local_to_global, spin_local_to_global)
      ! Return the SU(2) active rotation U satisfying
      ! U (sigma.v_local) U^dagger = sigma.(R v_local).
      ! Spin basis order is up,down. The deterministic quaternion sign chosen
      ! here is conventional: U and -U produce identical density transforms.
      REAL, INTENT(IN) :: local_to_global(3, 3)
      COMPLEX, INTENT(OUT) :: spin_local_to_global(local_d_n_spins, local_d_n_spins)

      REAL :: quaternion(4), norm_value

      CALL validate_rotation(local_to_global, "structural local-to-global rotation")
      CALL rotation_to_quaternion(local_to_global, quaternion)
      norm_value = SQRT(DOT_PRODUCT(quaternion, quaternion))
      quaternion = quaternion/norm_value
      CALL canonicalize_quaternion_sign(quaternion)
      spin_local_to_global(1, 1) = CMPLX(quaternion(1), -quaternion(4))
      spin_local_to_global(1, 2) = CMPLX(-quaternion(3), -quaternion(2))
      spin_local_to_global(2, 1) = CMPLX(quaternion(3), -quaternion(2))
      spin_local_to_global(2, 2) = CMPLX(quaternion(1), quaternion(4))
   END SUBROUTINE spinor_local_to_global

   SUBROUTINE local_d_spin_native_to_structural(local_to_global, native_mt_to_global, &
                                                orbital_global_to_local, spin_native_to_local, &
                                                combined_native_to_local)
      ! FLEUR nococonv%umat is native-MT -> global. The spin map into the
      ! structural frame is U_struct(local->global)^dagger*nococonv%umat.
      ! The optional combined matrix is ordered (m=-2..2, spin=1,2).
      REAL, INTENT(IN) :: local_to_global(3, 3)
      COMPLEX, INTENT(IN) :: native_mt_to_global(local_d_n_spins, local_d_n_spins)
      COMPLEX, INTENT(OUT) :: orbital_global_to_local(-local_d_l:local_d_l, &
                                                       -local_d_l:local_d_l)
      COMPLEX, INTENT(OUT) :: spin_native_to_local(local_d_n_spins, local_d_n_spins)
      COMPLEX, OPTIONAL, INTENT(OUT) :: combined_native_to_local(2*local_d_n_orbitals, &
                                                                 2*local_d_n_orbitals)

      COMPLEX :: spin_local_global(local_d_n_spins, local_d_n_spins)
      INTEGER :: m_local, m_native, spin_local, spin_native, row, column

      CALL validate_unitary(native_mt_to_global, "native-MT-to-global spin transformation")
      CALL local_d_orbital_global_to_local(local_to_global, orbital_global_to_local)
      CALL spinor_local_to_global(local_to_global, spin_local_global)
      spin_native_to_local = MATMUL(CONJG(TRANSPOSE(spin_local_global)), native_mt_to_global)
      IF (PRESENT(combined_native_to_local)) THEN
         DO m_local = -local_d_l, local_d_l
            DO spin_local = 1, local_d_n_spins
               row = (m_local + local_d_l)*local_d_n_spins + spin_local
               DO m_native = -local_d_l, local_d_l
                  DO spin_native = 1, local_d_n_spins
                     column = (m_native + local_d_l)*local_d_n_spins + spin_native
                     combined_native_to_local(row, column) = &
                        orbital_global_to_local(m_local, m_native)*spin_native_to_local(spin_local, spin_native)
                  END DO
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE local_d_spin_native_to_structural

   SUBROUTINE transform_local_d_spin_density(rho_native, local_to_global, native_mt_to_global, rho_structural, &
                                             orbital_global_to_local, spin_native_to_local)
      ! rho_structural=T rho_native T^dagger with
      ! T=D_global_to_local (x) [U_struct^dagger U_native_to_global].
      COMPLEX, INTENT(IN) :: rho_native(-local_d_l:, :, -local_d_l:, :)
      REAL, INTENT(IN) :: local_to_global(3, 3)
      COMPLEX, INTENT(IN) :: native_mt_to_global(local_d_n_spins, local_d_n_spins)
      COMPLEX, INTENT(OUT) :: rho_structural(-local_d_l:local_d_l, local_d_n_spins, &
                                             -local_d_l:local_d_l, local_d_n_spins)
      COMPLEX, OPTIONAL, INTENT(OUT) :: orbital_global_to_local(-local_d_l:local_d_l, &
                                                                -local_d_l:local_d_l)
      COMPLEX, OPTIONAL, INTENT(OUT) :: spin_native_to_local(local_d_n_spins, local_d_n_spins)

      COMPLEX :: orbital_transform(-local_d_l:local_d_l, -local_d_l:local_d_l)
      COMPLEX :: spin_transform(local_d_n_spins, local_d_n_spins)
      COMPLEX :: combined(2*local_d_n_orbitals, 2*local_d_n_orbitals)

      IF (SIZE(rho_native, 1) /= local_d_n_orbitals .OR. SIZE(rho_native, 2) /= local_d_n_spins .OR. &
          SIZE(rho_native, 3) /= local_d_n_orbitals .OR. SIZE(rho_native, 4) /= local_d_n_spins) THEN
         CALL juDFT_error("Native d-spin density must have dimensions (5,2,5,2)", &
                          calledby="m_local_d_spin_frame_transform")
      END IF
      CALL local_d_spin_native_to_structural(local_to_global, native_mt_to_global, &
                                             orbital_transform, spin_transform, combined)
      CALL transform_local_d_spin_density_with_matrix(rho_native,combined,rho_structural)
      IF (PRESENT(orbital_global_to_local)) orbital_global_to_local = orbital_transform
      IF (PRESENT(spin_native_to_local)) spin_native_to_local = spin_transform
   END SUBROUTINE transform_local_d_spin_density

   SUBROUTINE transform_local_d_spin_density_with_matrix(rho_native, combined_native_to_local, rho_structural)
      ! Apply a previously constructed common orbital-spin frame transform.
      ! This entry point permits a site-level transform to be cached and reused
      ! for many bands without changing the established transformation algebra.
      COMPLEX, INTENT(IN) :: rho_native(-local_d_l:, :, -local_d_l:, :)
      COMPLEX, INTENT(IN) :: combined_native_to_local(2*local_d_n_orbitals,2*local_d_n_orbitals)
      COMPLEX, INTENT(OUT) :: rho_structural(-local_d_l:local_d_l, local_d_n_spins, &
                                             -local_d_l:local_d_l, local_d_n_spins)

      COMPLEX :: rho_flat(2*local_d_n_orbitals, 2*local_d_n_orbitals)
      COMPLEX :: transformed(2*local_d_n_orbitals, 2*local_d_n_orbitals)
      INTEGER :: m, mp, ispin, jspin, row, column

      IF (SIZE(rho_native, 1) /= local_d_n_orbitals .OR. SIZE(rho_native, 2) /= local_d_n_spins .OR. &
          SIZE(rho_native, 3) /= local_d_n_orbitals .OR. SIZE(rho_native, 4) /= local_d_n_spins) THEN
         CALL juDFT_error("Native d-spin density must have dimensions (5,2,5,2)", &
                          calledby="m_local_d_spin_frame_transform")
      END IF
      DO m = -local_d_l, local_d_l
         DO ispin = 1, local_d_n_spins
            row = (m + local_d_l)*local_d_n_spins + ispin
            DO mp = -local_d_l, local_d_l
               DO jspin = 1, local_d_n_spins
                  column = (mp + local_d_l)*local_d_n_spins + jspin
                  rho_flat(row, column) = rho_native(m, ispin, mp, jspin)
               END DO
            END DO
         END DO
      END DO
      transformed = MATMUL(combined_native_to_local, &
                           MATMUL(rho_flat,CONJG(TRANSPOSE(combined_native_to_local))))
      DO m = -local_d_l, local_d_l
         DO ispin = 1, local_d_n_spins
            row = (m + local_d_l)*local_d_n_spins + ispin
            DO mp = -local_d_l, local_d_l
               DO jspin = 1, local_d_n_spins
                  column = (mp + local_d_l)*local_d_n_spins + jspin
                  rho_structural(m, ispin, mp, jspin) = transformed(row, column)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE transform_local_d_spin_density_with_matrix

   SUBROUTINE initialize_real_d_tensors(tensors)
      REAL, INTENT(OUT) :: tensors(3, 3, local_d_n_orbitals)
      REAL :: inverse_sqrt_two, inverse_sqrt_six

      inverse_sqrt_two = 1.0/SQRT(2.0)
      inverse_sqrt_six = 1.0/SQRT(6.0)
      tensors = 0.0
      tensors(1, 2, local_d_orb_xy) = inverse_sqrt_two
      tensors(2, 1, local_d_orb_xy) = inverse_sqrt_two
      tensors(2, 3, local_d_orb_yz) = inverse_sqrt_two
      tensors(3, 2, local_d_orb_yz) = inverse_sqrt_two
      tensors(1, 3, local_d_orb_xz) = inverse_sqrt_two
      tensors(3, 1, local_d_orb_xz) = inverse_sqrt_two
      tensors(1, 1, local_d_orb_x2y2) = inverse_sqrt_two
      tensors(2, 2, local_d_orb_x2y2) = -inverse_sqrt_two
      tensors(1, 1, local_d_orb_z2) = -inverse_sqrt_six
      tensors(2, 2, local_d_orb_z2) = -inverse_sqrt_six
      tensors(3, 3, local_d_orb_z2) = 2.0*inverse_sqrt_six
   END SUBROUTINE initialize_real_d_tensors

   SUBROUTINE rotation_to_quaternion(rotation, quaternion)
      REAL, INTENT(IN) :: rotation(3, 3)
      REAL, INTENT(OUT) :: quaternion(4)
      REAL :: scale, trace_value

      trace_value = rotation(1,1) + rotation(2,2) + rotation(3,3)
      IF (trace_value > 0.0) THEN
         scale = 2.0*SQRT(1.0 + trace_value)
         quaternion(1) = 0.25*scale
         quaternion(2) = (rotation(3,2) - rotation(2,3))/scale
         quaternion(3) = (rotation(1,3) - rotation(3,1))/scale
         quaternion(4) = (rotation(2,1) - rotation(1,2))/scale
      ELSE IF (rotation(1,1) > rotation(2,2) .AND. rotation(1,1) > rotation(3,3)) THEN
         scale = 2.0*SQRT(1.0 + rotation(1,1) - rotation(2,2) - rotation(3,3))
         quaternion(1) = (rotation(3,2) - rotation(2,3))/scale
         quaternion(2) = 0.25*scale
         quaternion(3) = (rotation(1,2) + rotation(2,1))/scale
         quaternion(4) = (rotation(1,3) + rotation(3,1))/scale
      ELSE IF (rotation(2,2) > rotation(3,3)) THEN
         scale = 2.0*SQRT(1.0 + rotation(2,2) - rotation(1,1) - rotation(3,3))
         quaternion(1) = (rotation(1,3) - rotation(3,1))/scale
         quaternion(2) = (rotation(1,2) + rotation(2,1))/scale
         quaternion(3) = 0.25*scale
         quaternion(4) = (rotation(2,3) + rotation(3,2))/scale
      ELSE
         scale = 2.0*SQRT(1.0 + rotation(3,3) - rotation(1,1) - rotation(2,2))
         quaternion(1) = (rotation(2,1) - rotation(1,2))/scale
         quaternion(2) = (rotation(1,3) + rotation(3,1))/scale
         quaternion(3) = (rotation(2,3) + rotation(3,2))/scale
         quaternion(4) = 0.25*scale
      END IF
   END SUBROUTINE rotation_to_quaternion

   SUBROUTINE canonicalize_quaternion_sign(quaternion)
      REAL, INTENT(INOUT) :: quaternion(4)
      INTEGER :: i

      DO i = 1, 4
         IF (ABS(quaternion(i)) > 100.0*EPSILON(quaternion(i))) THEN
            IF (quaternion(i) < 0.0) quaternion = -quaternion
            RETURN
         END IF
      END DO
   END SUBROUTINE canonicalize_quaternion_sign

   SUBROUTINE validate_rotation(rotation, description)
      REAL, INTENT(IN) :: rotation(3, 3)
      CHARACTER(LEN=*), INTENT(IN) :: description
      REAL :: determinant

      determinant = rotation(1,1)*(rotation(2,2)*rotation(3,3)-rotation(2,3)*rotation(3,2)) &
                  - rotation(1,2)*(rotation(2,1)*rotation(3,3)-rotation(2,3)*rotation(3,1)) &
                  + rotation(1,3)*(rotation(2,1)*rotation(3,2)-rotation(2,2)*rotation(3,1))
      IF (MAXVAL(ABS(MATMUL(TRANSPOSE(rotation), rotation) - identity_3())) > 1.0e-8 .OR. &
          ABS(determinant - 1.0) > 1.0e-8) THEN
         CALL juDFT_error(TRIM(description)//" must be a proper orthogonal matrix", &
                          calledby="m_local_d_spin_frame_transform")
      END IF
   END SUBROUTINE validate_rotation

   SUBROUTINE validate_unitary(unitary, description)
      COMPLEX, INTENT(IN) :: unitary(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: description
      COMPLEX :: identity(SIZE(unitary, 2), SIZE(unitary, 2))
      INTEGER :: i

      IF (SIZE(unitary, 1) /= SIZE(unitary, 2)) THEN
         CALL juDFT_error(TRIM(description)//" must be square", calledby="m_local_d_spin_frame_transform")
      END IF
      identity = CMPLX(0.0, 0.0)
      DO i = 1, SIZE(identity, 1)
         identity(i, i) = CMPLX(1.0, 0.0)
      END DO
      IF (MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(unitary)), unitary) - identity)) > 1.0e-8) THEN
         CALL juDFT_error(TRIM(description)//" must be unitary", calledby="m_local_d_spin_frame_transform")
      END IF
   END SUBROUTINE validate_unitary

   FUNCTION identity_3() RESULT(identity)
      REAL :: identity(3, 3)
      INTEGER :: i

      identity = 0.0
      DO i = 1, 3
         identity(i, i) = 1.0
      END DO
   END FUNCTION identity_3

END MODULE m_local_d_spin_frame_transform
