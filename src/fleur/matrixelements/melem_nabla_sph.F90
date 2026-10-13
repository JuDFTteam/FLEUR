!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The muffin-tin half of <u_a| e^{-i b.r} p |u_b>, from the table the gradient was
!>  tabulated into. Same shape as the plain overlap: a phase per atom and a matrix product
!>  per pair of radial slots, with one sum more because the gradient has three components.
!>
!>  The table holds the gradient in spherical components. Both the turn to cartesian ones
!>  and the -i that makes the gradient into the momentum happen here, so that what leaves
!>  is the same quantity the interstitial half returns and the two can simply be added.
MODULE m_melem_nabla_sph
   USE m_juDFT
   USE m_types_abc
   USE m_types_atoms
   USE m_constants, ONLY: tpi_const
   USE m_types_sym
   USE m_types_cell
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_nabla_sph
CONTAINS

   !> Accumulates into pnk, which the interstitial half has usually written into first.
   !>
   !> kdiff has to hold the difference the two sides span, exactly as for the overlap: the
   !> table is indexed by it and a missing entry is a stop rather than a guess.
   SUBROUTINE melem_nabla_sph(atoms, sym, cell, abc, abc_b, bbpt, gb, bkpt, nujug, kdiff, nntot, pnk)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_abc), INTENT(IN) :: abc(:)
      TYPE(t_abc), INTENT(IN) :: abc_b(:)
      REAL, INTENT(IN) :: bbpt(3)
      INTEGER, INTENT(IN) :: gb(3)
      REAL, INTENT(IN) :: bkpt(3)
      COMPLEX, INTENT(IN) :: nujug(0:, 0:, :, :, :, :, -1:)
      REAL, INTENT(IN) :: kdiff(:, :)
      INTEGER, INTENT(IN) :: nntot
      COMPLEX, INTENT(INOUT) :: pnk(:, :, :)   !> (nbnd, nbnd, 3), cartesian

      INTEGER :: n, nn, na, lm, r1, r2, nene, a, mu, nbnd
      REAL :: rph, cph, th
      REAL :: bpt(3)
      COMPLEX :: phasefac
      COMPLEX, ALLOCATABLE :: sph(:, :, :), loc(:, :, :)
      COMPLEX :: hlp(3, -1:1)
      COMPLEX, PARAMETER :: img = CMPLX(0.0, 1.0)
      REAL :: s2

      CALL timestart("melem_nabla_sph")

      nbnd = SIZE(pnk, 1)
      bpt(:) = bbpt(:) + gb(:) - bkpt(:)
      DO nene = 1, nntot
         IF (ALL(ABS(bpt(:) - kdiff(:, nene)) < 1e-4)) EXIT
      END DO
      IF (nene == nntot + 1) CALL juDFT_error("cannot find matching nearest neighbor k", &
                                              calledby="melem_nabla_sph")

      ALLOCATE (sph(nbnd, nbnd, -1:1), loc(nbnd, nbnd, 3))

      !> Spherical to cartesian, and the -i of the momentum. Hoisted out of the loop because
      !> it does not depend on the atom; what does depend on it is the rotation below.
      s2 = 1.0/SQRT(2.0)
      hlp = CMPLX(0.0, 0.0)
      hlp(1, -1) = CMPLX(s2, 0.0);   hlp(1, 1) = CMPLX(-s2, 0.0)
      hlp(2, -1) = -img*s2;          hlp(2, 1) = -img*s2
      hlp(3, 0) = CMPLX(1.0, 0.0)

      DO n = 1, atoms%ntype
         DO nn = 1, atoms%neq(n)
            na = atoms%firstAtom(n) + nn - 1
            th = tpi_const*DOT_PRODUCT(bpt, atoms%taual(:, na))
            rph = 2.0*tpi_const*COS(th)
            cph = 2.0*tpi_const*SIN(th)
            phasefac = CMPLX(rph, cph)

            !> Per atom now, not accumulated over them: the three components that come
            !> out of the contraction are a VECTOR in THIS atom's local frame, and vectors
            !> in different frames cannot be added. Summing first and converting once was
            !> the bug: with equivalent atoms related by inversion the two contributions
            !> cancelled instead of adding, and the muffin-tin term vanished.
            sph = CMPLX(0.0, 0.0)
            lm = atoms%lmax(n)*(atoms%lmax(n) + 2)
            DO mu = -1, 1
               DO r1 = 1, MAXVAL(abc(n)%n_r)
                  DO r2 = 1, MAXVAL(abc_b(n)%n_r)
                     ! Same contraction as the overlap, one per spherical component. The ket
                     ! side is conjugated here, which is the convention the assembled matrix
                     ! is undone by further up.
                     sph(:, :, mu) = sph(:, :, mu) + phasefac* &
                                     MATMUL(abc(n)%cof(:, 0:lm, r1, nn), &
                                            TRANSPOSE(MATMUL(CONJG(abc_b(n)%cof(:, 0:lm, r2, nn)), &
                                                             nujug(0:lm, 0:lm, r1, r2, n, nene, mu))))
                  END DO
               END DO
            END DO
            !> The closing is conjugated because the contraction above is, and the -i of
            !> the momentum becomes +i for the same reason. Done HERE, in the local frame,
            !> because what follows has to act on cartesian components.
            DO a = 1, 3
               loc(:, :, a) = img*(CONJG(hlp(a, -1))*sph(:, :, -1) &
                                   + CONJG(hlp(a, 0))*sph(:, :, 0) &
                                   + CONJG(hlp(a, 1))*sph(:, :, 1))
            END DO

            !> The components add straight in because the coefficients arrive expanded
            !> about the global frame, so the contraction is already global. Were they
            !> expanded about each atom's own frame this would be wrong: p is a vector and
            !> its components mix under that atom's operation, so each one would have to be
            !> turned before being summed with the others.
            DO a = 1, 3
               pnk(:, :, a) = pnk(:, :, a) + loc(:, :, a)
            END DO
         END DO
      END DO

      DEALLOCATE (sph, loc)
      CALL timestop("melem_nabla_sph")
   END SUBROUTINE melem_nabla_sph

END MODULE m_melem_nabla_sph
