!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_angular
   USE m_juDFT, ONLY: juDFT_error
   USE m_clebsch, ONLY: clebsch
   USE m_constants, ONLY: pi_const
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: xas_sigma_down = -1
   INTEGER, PARAMETER, PUBLIC :: xas_sigma_up = 1
   INTEGER, PARAMETER, PUBLIC :: xas_edge_l2_twice_j = 1
   INTEGER, PARAMETER, PUBLIC :: xas_edge_l3_twice_j = 3

   INTEGER, PARAMETER :: xas_gaunt_lmax = 12
   REAL, PARAMETER    :: xas_sqrt_half = 0.70710678118654752440
   REAL, PARAMETER    :: xas_check_tol = 1.0e-8
   COMPLEX, PARAMETER :: ImagUnit = CMPLX(0.0, 1.0)

   PUBLIC :: xas_cartesian_to_spherical
   PUBLIC :: xas_dipole_angular_coeff
   PUBLIC :: xas_dipole_angular_table
   PUBLIC :: xas_angular_selfcheck
   PUBLIC :: xas_print_angular_sumrule

CONTAINS

   SUBROUTINE xas_cartesian_to_spherical(eps_cart, eps_sph)
      ! Converts Cartesian components to spherical components in the convention
      ! eps_{+1}=-(eps_x+i eps_y)/sqrt(2), eps_0=eps_z,
      ! eps_{-1}= (eps_x-i eps_y)/sqrt(2).
      COMPLEX, INTENT(IN)  :: eps_cart(3)
      COMPLEX, INTENT(OUT) :: eps_sph(-1:1)

      eps_sph(-1) = xas_sqrt_half*(eps_cart(1) - ImagUnit*eps_cart(2))
      eps_sph(0) = eps_cart(3)
      eps_sph(1) = -xas_sqrt_half*(eps_cart(1) + ImagUnit*eps_cart(2))
   END SUBROUTINE xas_cartesian_to_spherical

   COMPLEX FUNCTION xas_dipole_angular_coeff(l, m, lc, twice_j, twice_mj, sigma, eps_sph) RESULT(coeff)
      ! Computes sum_q sum_mc eps_q <lc mc, 1/2 sigma | j mj>
      !                     * int dOmega Y_lm^* Y_1q Y_lc,mc.
      !
      ! Clebsch convention: <lc mc, 1/2 ms | j mj>.
      ! The half-integer quantum numbers are represented exactly by doubled
      ! integers: sigma=2 m_s = +/-1, twice_j=2j, twice_mj=2m_j.
      INTEGER, INTENT(IN) :: l, m, lc, twice_j, twice_mj, sigma
      COMPLEX, INTENT(IN) :: eps_sph(-1:1)

      INTEGER :: mc, q, lmax_needed
      REAL    :: cg, gnt_complex
      COMPLEX :: contribution

      IF (sigma /= xas_sigma_down .AND. sigma /= xas_sigma_up) THEN
         CALL juDFT_error("sigma must be +/-1 in xas_dipole_angular_coeff", calledby="m_xas_angular")
      END IF
      IF (ABS(m) > l) THEN
         coeff = CMPLX(0.0, 0.0)
         RETURN
      END IF
      IF (ABS(twice_mj) > twice_j) THEN
         coeff = CMPLX(0.0, 0.0)
         RETURN
      END IF

      lmax_needed = MAX(l, lc, 1)
      IF (lmax_needed > xas_gaunt_lmax) THEN
         CALL juDFT_error("Increase xas_gaunt_lmax in m_xas_angular", calledby="m_xas_angular")
      END IF

      coeff = CMPLX(0.0, 0.0)
      DO q = -1, 1
         IF (ABS(eps_sph(q)) < xas_check_tol) CYCLE
         DO mc = -lc, lc
            cg = clebsch(REAL(lc), 0.5, REAL(mc), 0.5*REAL(sigma), 0.5*REAL(twice_j), 0.5*REAL(twice_mj))
            IF (ABS(cg) < xas_check_tol) CYCLE
            gnt_complex = xas_complex_gaunt(l, m, 1, q, lc, mc)
            IF (ABS(gnt_complex) < xas_check_tol) CYCLE
            contribution = eps_sph(q)*cg*gnt_complex
            coeff = coeff + contribution
         END DO
      END DO
   END FUNCTION xas_dipole_angular_coeff

   REAL FUNCTION xas_complex_gaunt(l1, m1, l2, m2, l3, m3) RESULT(gnt)
      ! Complex spherical-harmonic Gaunt coefficient
      !
      !   int dOmega conjg(Y_l1,m1) Y_l2,m2 Y_l3,m3
      !
      ! evaluated through Wigner 3j symbols expressed with FLEUR's Clebsch-
      ! Gordan routine. This is the convention matched to the complex
      ! spherical polarization components used by XAS.
      INTEGER, INTENT(IN) :: l1, m1, l2, m2, l3, m3

      REAL :: tj000, tjm

      gnt = 0.0
      IF (ABS(m1) > l1 .OR. ABS(m2) > l2 .OR. ABS(m3) > l3) RETURN
      IF (-m1 + m2 + m3 /= 0) RETURN
      IF (MOD(l1 + l2 + l3, 2) == 1) RETURN
      IF ((l1 + l2 - l3) < 0) RETURN
      IF ((l1 - l2 + l3) < 0) RETURN
      IF ((-l1 + l2 + l3) < 0) RETURN

      tj000 = xas_wigner_3j_from_clebsch(l1, l2, l3, 0, 0, 0)
      IF (ABS(tj000) < xas_check_tol) RETURN
      tjm = xas_wigner_3j_from_clebsch(l1, l2, l3, -m1, m2, m3)
      IF (ABS(tjm) < xas_check_tol) RETURN

      gnt = xas_phase(m1)*SQRT(REAL((2*l1 + 1)*(2*l2 + 1)*(2*l3 + 1))/(4.0*pi_const))*tj000*tjm
   END FUNCTION xas_complex_gaunt

   REAL FUNCTION xas_wigner_3j_from_clebsch(j1, j2, j3, m1, m2, m3) RESULT(tj)
      INTEGER, INTENT(IN) :: j1, j2, j3, m1, m2, m3

      tj = 0.0
      IF (m1 + m2 + m3 /= 0) RETURN
      tj = xas_phase(j1 - j2 - m3)*clebsch(REAL(j1), REAL(j2), REAL(m1), REAL(m2), REAL(j3), REAL(-m3)) &
           /SQRT(REAL(2*j3 + 1))
   END FUNCTION xas_wigner_3j_from_clebsch

   INTEGER FUNCTION xas_phase(n) RESULT(phase)
      INTEGER, INTENT(IN) :: n

      IF (MOD(ABS(n), 2) == 0) THEN
         phase = 1
      ELSE
         phase = -1
      END IF
   END FUNCTION xas_phase

   SUBROUTINE xas_dipole_angular_table(lmax, lc, twice_j, twice_mj, eps_sph, coeff)
      ! coeff(i_sigma,lm+1), with i_sigma=1 for sigma=-1 and
      ! i_sigma=2 for sigma=+1.
      INTEGER, INTENT(IN)  :: lmax, lc, twice_j, twice_mj
      COMPLEX, INTENT(IN)  :: eps_sph(-1:1)
      COMPLEX, INTENT(OUT) :: coeff(:, :)

      INTEGER :: l, m, lm, sigma

      coeff = CMPLX(0.0, 0.0)
      DO sigma = -1, 1, 2
         DO l = 0, lmax
            DO m = -l, l
               lm = l*(l + 1) + m
               coeff((sigma + 3)/2, lm + 1) = xas_dipole_angular_coeff(l, m, lc, twice_j, twice_mj, sigma, eps_sph)
            END DO
         END DO
      END DO
   END SUBROUTINE xas_dipole_angular_table

   LOGICAL FUNCTION xas_angular_selfcheck(verbose) RESULT(ok)
      LOGICAL, OPTIONAL, INTENT(IN) :: verbose

      COMPLEX :: eps_cart(3), eps_sph(-1:1), c
      INTEGER :: m
      LOGICAL :: l_verbose

      l_verbose = .FALSE.
      IF (PRESENT(verbose)) l_verbose = verbose

      ok = .TRUE.

      eps_cart = CMPLX(0.0, 0.0)
      eps_cart(3) = CMPLX(1.0, 0.0)
      CALL xas_cartesian_to_spherical(eps_cart, eps_sph)
      ok = ok .AND. ABS(eps_sph(0) - CMPLX(1.0, 0.0)) < xas_check_tol
      ok = ok .AND. ABS(eps_sph(-1)) < xas_check_tol
      ok = ok .AND. ABS(eps_sph(1)) < xas_check_tol

      eps_cart = CMPLX(0.0, 0.0)
      eps_cart(1) = CMPLX(xas_sqrt_half, 0.0)
      eps_cart(2) = ImagUnit*xas_sqrt_half
      CALL xas_cartesian_to_spherical(eps_cart, eps_sph)
      ok = ok .AND. ABS(eps_sph(-1) - CMPLX(1.0, 0.0)) < xas_check_tol
      ok = ok .AND. ABS(eps_sph(0)) < xas_check_tol
      ok = ok .AND. ABS(eps_sph(1)) < xas_check_tol

      eps_cart = CMPLX(0.0, 0.0)
      eps_cart(1) = CMPLX(xas_sqrt_half, 0.0)
      eps_cart(2) = -ImagUnit*xas_sqrt_half
      CALL xas_cartesian_to_spherical(eps_cart, eps_sph)
      ok = ok .AND. ABS(eps_sph(1) + CMPLX(1.0, 0.0)) < xas_check_tol
      ok = ok .AND. ABS(eps_sph(0)) < xas_check_tol
      ok = ok .AND. ABS(eps_sph(-1)) < xas_check_tol

      eps_cart = CMPLX(0.0, 0.0)
      eps_cart(3) = CMPLX(1.0, 0.0)
      CALL xas_cartesian_to_spherical(eps_cart, eps_sph)
      DO m = -1, 1
         c = xas_dipole_angular_coeff(1, m, 1, xas_edge_l3_twice_j, 1, xas_sigma_up, eps_sph)
         ok = ok .AND. ABS(c) < xas_check_tol
      END DO
      c = xas_dipole_angular_coeff(0, 0, 1, xas_edge_l3_twice_j, 1, xas_sigma_up, eps_sph)
      ok = ok .AND. ABS(c) > xas_check_tol
      c = xas_dipole_angular_coeff(2, 0, 1, xas_edge_l3_twice_j, 1, xas_sigma_up, eps_sph)
      ok = ok .AND. ABS(c) > xas_check_tol

      IF (l_verbose) THEN
         IF (ok) THEN
            WRITE (*, '(a)') "m_xas_angular self-check passed."
         ELSE
            WRITE (*, '(a)') "m_xas_angular self-check failed."
         END IF
      END IF
   END FUNCTION xas_angular_selfcheck

   SUBROUTINE xas_print_angular_sumrule(lc, twice_j, unit)
      ! Diagnostic-only pure angular sum rule for dipole XAS. This uses the
      ! production angular coefficient C but no bands, abc coefficients,
      ! radial integrals, occupations, or k points.
      INTEGER, INTENT(IN) :: lc, twice_j
      INTEGER, INTENT(IN) :: unit

      COMPLEX :: eps_cart(3), eps_sph(-1:1), coeff
      REAL    :: strength(2, 3), strength_mj(2, 3, twice_j + 1)
      REAL    :: avg, rel_xy, rel_xz, rel_yz
      INTEGER :: i_l, n_l_values, l_values(2), i_pol, i_mj, twice_mj, sigma, l, m

      n_l_values = 0
      IF (lc - 1 >= 0) THEN
         n_l_values = n_l_values + 1
         l_values(n_l_values) = lc - 1
      END IF
      n_l_values = n_l_values + 1
      l_values(n_l_values) = lc + 1
      strength = 0.0
      strength_mj = 0.0

      DO i_pol = 1, 3
         eps_cart = CMPLX(0.0, 0.0)
         eps_cart(i_pol) = CMPLX(1.0, 0.0)
         CALL xas_cartesian_to_spherical(eps_cart, eps_sph)
         DO i_l = 1, n_l_values
            l = l_values(i_l)
            DO i_mj = 1, twice_j + 1
               twice_mj = -twice_j + 2*(i_mj - 1)
               DO sigma = -1, 1, 2
                  DO m = -l, l
                     coeff = xas_dipole_angular_coeff(l, m, lc, twice_j, twice_mj, sigma, eps_sph)
                     strength(i_l, i_pol) = strength(i_l, i_pol) + ABS(coeff)**2
                     strength_mj(i_l, i_pol, i_mj) = strength_mj(i_l, i_pol, i_mj) + ABS(coeff)**2
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL xas_write_sumrule_line(unit, "XAS DEBUG ANGULAR SUMRULE header: lc twice_j", REAL(lc), REAL(twice_j))
      DO i_l = 1, n_l_values
         CALL xas_angular_relative(strength(i_l, :), avg, rel_xy, rel_xz, rel_yz)
         CALL xas_write_sumrule_strength(unit, "XAS DEBUG ANGULAR SUMRULE total", l_values(i_l), 0, &
                                         strength(i_l, :), avg, rel_xy, rel_xz, rel_yz)
         DO i_mj = 1, twice_j + 1
            twice_mj = -twice_j + 2*(i_mj - 1)
            CALL xas_angular_relative(strength_mj(i_l, :, i_mj), avg, rel_xy, rel_xz, rel_yz)
            CALL xas_write_sumrule_strength(unit, "XAS DEBUG ANGULAR SUMRULE mj", l_values(i_l), twice_mj, &
                                            strength_mj(i_l, :, i_mj), avg, rel_xy, rel_xz, rel_yz)
         END DO
      END DO
   END SUBROUTINE xas_print_angular_sumrule

   SUBROUTINE xas_angular_relative(strength, avg, rel_xy, rel_xz, rel_yz)
      REAL, INTENT(IN)  :: strength(3)
      REAL, INTENT(OUT) :: avg, rel_xy, rel_xz, rel_yz

      avg = SUM(strength)/3.0
      rel_xy = 0.0
      rel_xz = 0.0
      rel_yz = 0.0
      IF (ABS(avg) > TINY(avg)) THEN
         rel_xy = (strength(1) - strength(2))/avg
         rel_xz = (strength(1) - strength(3))/avg
         rel_yz = (strength(2) - strength(3))/avg
      END IF
   END SUBROUTINE xas_angular_relative

   SUBROUTINE xas_write_sumrule_strength(unit, label, l, twice_mj, strength, avg, rel_xy, rel_xz, rel_yz)
      INTEGER,          INTENT(IN) :: unit, l, twice_mj
      CHARACTER(LEN=*), INTENT(IN) :: label
      REAL,             INTENT(IN) :: strength(3), avg, rel_xy, rel_xz, rel_yz

      WRITE(*, '(a,a,i0,a,i0,a,3es18.10,a,es18.10,a,3es12.4)') TRIM(label), " l=", l, &
         " twice_mj=", twice_mj, " Ax Ay Az=", strength, " avg=", avg, " rel_xy rel_xz rel_yz=", &
         rel_xy, rel_xz, rel_yz
      IF (unit > 0) THEN
         WRITE(unit, '(a,a,i0,a,i0,a,3es18.10,a,es18.10,a,3es12.4)') TRIM(label), " l=", l, &
            " twice_mj=", twice_mj, " Ax Ay Az=", strength, " avg=", avg, " rel_xy rel_xz rel_yz=", &
            rel_xy, rel_xz, rel_yz
      END IF
   END SUBROUTINE xas_write_sumrule_strength

   SUBROUTINE xas_write_sumrule_line(unit, label, value1, value2)
      INTEGER,          INTENT(IN) :: unit
      CHARACTER(LEN=*), INTENT(IN) :: label
      REAL,             INTENT(IN) :: value1, value2

      WRITE(*, '(a,2es18.10)') TRIM(label)//" ", value1, value2
      IF (unit > 0) WRITE(unit, '(a,2es18.10)') TRIM(label)//" ", value1, value2
   END SUBROUTINE xas_write_sumrule_line

END MODULE m_xas_angular
