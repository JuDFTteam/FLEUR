!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmkb_sph
   USE m_juDFT
   USE m_types_abc
   USE m_types_atoms
   USE m_constants, ONLY: tpi_const
   IMPLICIT NONE
CONTAINS

   SUBROUTINE wannierlib_mmkb_sph(atoms,  abc, abc_b, bbpt, gb, bkpt, &
                                  ujug, &
                                  kdiff, nntot, mmn)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_abc), INTENT(IN) :: abc(:)
      TYPE(t_abc), INTENT(IN) :: abc_b(:)
      REAL, INTENT(IN) :: bbpt(3)
      INTEGER, INTENT(IN) :: gb(3)
      REAL, INTENT(IN) :: bkpt(3)
      COMPLEX, INTENT(IN) :: ujug(:, :, 0:, 0:, :, :)
      REAL, INTENT(IN) :: kdiff(:, :)
      INTEGER, INTENT(IN) :: nntot
      COMPLEX, INTENT(INOUT) :: mmn(:, :)

      INTEGER :: i, j, n, nn, na, l, lp, m, mp, lm, lmp, ll, llp, r1, r2
      INTEGER :: nene
      REAL :: rph, cph, th
      REAL :: bpt(3)
      COMPLEX :: phasefac

      CALL timestart("wannierlib_mmkb_sph")

     
      bpt(:) = bbpt(:) + gb(:) - bkpt(:)
      DO nene = 1, nntot
         IF (ALL(ABS(bpt(:) - kdiff(:, nene)) < 1e-4)) EXIT
      END DO
      IF (nene == nntot + 1) CALL juDFT_error("cannot find matching nearest neighbor k", calledby="wannierlib_mmkb_sph")

      DO n = 1, atoms%ntype
         DO nn = 1, atoms%neq(n)
            na = atoms%firstAtom(n) + nn-1
            th=tpi_const*dot_product(bpt, atoms%taual(:, na))
            rph = 2.0*tpi_const*COS(th)
            cph = 2.0*tpi_const*SIN(th)
            phasefac = CMPLX(rph, cph)

            DO l = 0, atoms%lmax(n)
               ll = l*(l + 1)
               DO m = -l, l
                  lm = ll + m
                  DO lp = 0, atoms%lmax(n)
                     llp = lp*(lp + 1)
                     DO mp = -lp, lp
                        lmp = llp + mp
                        DO r1 = 1, abc(n)%n_r(l)
                           DO r2 = 1, abc_b(n)%n_r(lp)

                              DO j = 1, size(mmn,2)
                                 DO i = 1, size(mmn,1)
                                    mmn(i, j) = mmn(i, j) + phasefac* &
                                        abc(n)%cof(i, lm, r1, nn)*CONJG(abc_b(n)%cof(j, lmp, r2, nn))*ujug(r1, r2, lmp, lm, n, nene)
                                 END DO
                              END DO
                           end do
                        end do
                     END DO
                  END DO
               END DO
            END DO

         END DO
      END DO

      CALL timestop("wannierlib_mmkb_sph")
   END SUBROUTINE wannierlib_mmkb_sph

END MODULE m_wannierlib_mmkb_sph
