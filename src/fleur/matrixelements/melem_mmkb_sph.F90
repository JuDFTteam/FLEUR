!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The muffin-tin half of <u_a| e^{-i b.r} |u_b>: a phase per atom and a matrix product
!>  per pair of radial slots, against the table m_melem_ujugaunt tabulated beforehand.
!>
!>  The b it works at is found in kdiff BY VALUE rather than by neighbour index: kdiff is
!>  deduplicated, so its slot order and the neighbour order coincide only by accident.
!>
!>  Accumulates into the same mmn the interstitial half fills.
MODULE m_melem_mmkb_sph
   USE m_juDFT
   USE m_types_abc
   USE m_types_atoms
   USE m_constants, ONLY: tpi_const
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_mmkb_sph
CONTAINS

   SUBROUTINE melem_mmkb_sph(atoms, abc, abc_b, bbpt, gb, bkpt, &
                                  ujug, &
                                  kdiff, nntot, mmn)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_abc), INTENT(IN) :: abc(:)
      TYPE(t_abc), INTENT(IN) :: abc_b(:)
      REAL, INTENT(IN) :: bbpt(3)
      INTEGER, INTENT(IN) :: gb(3)
      REAL, INTENT(IN) :: bkpt(3)
      COMPLEX, INTENT(IN) :: ujug(0:, 0:, :, :, :, :)
      REAL, INTENT(IN) :: kdiff(:, :)
      INTEGER, INTENT(IN) :: nntot
      COMPLEX, INTENT(INOUT) :: mmn(:, :)

      INTEGER :: n, nn, na, lm, r1, r2
      INTEGER :: nene
      REAL :: rph, cph, th
      REAL :: bpt(3)
      COMPLEX :: phasefac

      CALL timestart("melem_mmkb_sph")

      bpt(:) = bbpt(:) + gb(:) - bkpt(:)
      DO nene = 1, nntot
         IF (ALL(ABS(bpt(:) - kdiff(:, nene)) < 1e-4)) EXIT
      END DO
      IF (nene == nntot + 1) CALL juDFT_error("cannot find matching nearest neighbor k", calledby="melem_mmkb_sph")

      DO n = 1, atoms%ntype
         DO nn = 1, atoms%neq(n)
            na = atoms%firstAtom(n) + nn - 1
            th = tpi_const*dot_product(bpt, atoms%taual(:, na))
            rph = 2.0*tpi_const*COS(th)
            cph = 2.0*tpi_const*SIN(th)
            phasefac = CMPLX(rph, cph)

            lm = atoms%lmax(n)*(atoms%lmax(n) + 2)
            DO r1 = 1, maxval(abc(n)%n_r)
               DO r2 = 1, maxval(abc_b(n)%n_r)
                  ! M_MT = phasefac * A . (conjg(B).ujug)^T. The ket side is conjugated
                  ! here and the assembled overlap is conjugated once more by whoever
                  ! collects it; both halves of that convention are needed together.
                  mmn = mmn + phasefac*matmul(abc(n)%cof(:, 0:lm, r1, nn), &
                        transpose(matmul(conjg(abc_b(n)%cof(:, 0:lm, r2, nn)), ujug(0:lm, 0:lm, r1, r2, n, nene))))
               enddo
            END DO

         END DO
      END DO

      CALL timestop("melem_mmkb_sph")
   END SUBROUTINE melem_mmkb_sph

END MODULE m_melem_mmkb_sph
