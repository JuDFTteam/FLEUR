!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_radial
   USE m_intgr, ONLY: intgr3
   USE m_juDFT, ONLY: juDFT_error
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_radfun, ONLY: t_radfun
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: xas_radial_dipole_integrals
   PUBLIC :: xas_radial_print_largest

CONTAINS

   SUBROUTINE xas_radial_dipole_integrals(atoms, itype, radfun, p_core, radial_int)
      ! Computes int_0^RMT rho U_l(rho) P_core(rho) d rho for all l,
      ! spin, and radial-function channels stored in t_radfun.
      !
      ! U_l and P_core are FLEUR reduced radial functions. This corresponds
      ! to int rho^3 u_l(rho) R_core(rho) d rho with U_l=rho*u_l and
      ! P_core=rho*R_core.
      !
      ! radial_int(order,l,spin):
      !   order=1 is u_l, order=2 is dot u_l, order>=3 are LO channels
      !   in the same order as radfun%r/order and radfun%n_r(l).
      TYPE(t_atoms),  INTENT(IN)  :: atoms
      INTEGER,        INTENT(IN)  :: itype
      TYPE(t_radfun), INTENT(IN)  :: radfun
      REAL,           INTENT(IN)  :: p_core(:)
      REAL,           INTENT(OUT) :: radial_int(:, 0:, :)

      INTEGER :: ispin, l, order, jri, max_order
      REAL    :: integrand(atoms%jri(itype))

      IF (.NOT. ALLOCATED(radfun%n_r)) THEN
         CALL juDFT_error("radfun%n_r is not allocated in xas_radial_dipole_integrals", calledby="m_xas_radial")
      END IF
      IF (.NOT. ALLOCATED(radfun%r)) THEN
         CALL juDFT_error("radfun%r is not allocated in xas_radial_dipole_integrals", calledby="m_xas_radial")
      END IF

      jri = atoms%jri(itype)
      IF (SIZE(p_core) < jri) THEN
         CALL juDFT_error("p_core is shorter than the MT radial mesh in xas_radial_dipole_integrals", calledby="m_xas_radial")
      END IF
      IF (UBOUND(radial_int, 2) < atoms%lmax(itype)) THEN
         CALL juDFT_error("radial_int l dimension is too small in xas_radial_dipole_integrals", calledby="m_xas_radial")
      END IF
      IF (SIZE(radial_int, 3) < SIZE(radfun%r, 5)) THEN
         CALL juDFT_error("radial_int spin dimension is too small in xas_radial_dipole_integrals", calledby="m_xas_radial")
      END IF
      max_order = MAXVAL(radfun%n_r(0:atoms%lmax(itype)))
      IF (SIZE(radial_int, 1) < max_order) THEN
         CALL juDFT_error("radial_int radial-order dimension is too small in xas_radial_dipole_integrals", calledby="m_xas_radial")
      END IF

      radial_int = 0.0
      DO ispin = 1, SIZE(radfun%r, 5)
         DO l = 0, atoms%lmax(itype)
            DO order = 1, radfun%n_r(l)
               integrand(1:jri) = atoms%rmsh(1:jri, itype) * radfun%r(1:jri, 1, order, l, ispin) * p_core(1:jri)
               CALL intgr3(integrand, atoms%rmsh(1:jri, itype), atoms%dx(itype), jri, radial_int(order, l, ispin))
            END DO
         END DO
      END DO
   END SUBROUTINE xas_radial_dipole_integrals

   SUBROUTINE xas_radial_print_largest(radial_int, n_print, out_unit)
      REAL,    INTENT(IN)           :: radial_int(:, 0:, :)
      INTEGER, INTENT(IN), OPTIONAL :: n_print
      INTEGER, INTENT(IN), OPTIONAL :: out_unit

      INTEGER :: iprint, nout, unit, order, l, ispin
      INTEGER :: best_order, best_l, best_spin
      LOGICAL :: used(SIZE(radial_int, 1), LBOUND(radial_int, 2):UBOUND(radial_int, 2), SIZE(radial_int, 3))
      REAL    :: best_abs

      nout = 10
      IF (PRESENT(n_print)) nout = n_print
      unit = 6
      IF (PRESENT(out_unit)) unit = out_unit

      used = .FALSE.
      WRITE (unit, '(a)') "Largest XAS radial dipole integrals:"
      DO iprint = 1, MIN(nout, SIZE(radial_int))
         best_abs = -1.0
         best_order = 0
         best_l = 0
         best_spin = 0
         DO ispin = 1, SIZE(radial_int, 3)
            DO l = LBOUND(radial_int, 2), UBOUND(radial_int, 2)
               DO order = 1, SIZE(radial_int, 1)
                  IF (used(order, l, ispin)) CYCLE
                  IF (ABS(radial_int(order, l, ispin)) > best_abs) THEN
                     best_abs = ABS(radial_int(order, l, ispin))
                     best_order = order
                     best_l = l
                     best_spin = ispin
                  END IF
               END DO
            END DO
         END DO
         IF (best_abs < 0.0) EXIT
         used(best_order, best_l, best_spin) = .TRUE.
         WRITE (unit, '(a,i4,a,i4,a,i4,a,es16.8,a,es16.8)') "  order=", best_order, " l=", best_l, &
            " spin=", best_spin, " value=", radial_int(best_order, best_l, best_spin), &
            " abs=", best_abs
      END DO
   END SUBROUTINE xas_radial_print_largest

END MODULE m_xas_radial
