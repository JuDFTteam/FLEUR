!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_mmk0_sph
   USE m_juDFT
   USE m_types_abc
   USE m_types_atoms
   USE m_types_radfun
   IMPLICIT NONE
CONTAINS

   SUBROUTINE wannierlib_mmk0_sph(atoms, noccbd, abc, radfun, mmn)
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: noccbd
      TYPE(t_abc), INTENT(IN) :: abc(:)
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      COMPLEX, INTENT(INOUT) :: mmn(:, :)

      INTEGER :: i, j, l, m, natom, ntyp, lm, ll1, n_r, n_r2
      COMPLEX :: sumab

      CALL timestart("wannierlib_mmk0_sph")

      DO j = 1, size(mmn, 2)
         DO i = 1, size(mmn, 1)
            DO ntyp = 1, atoms%ntype
               DO l = 0, atoms%lmax(ntyp)
                  ll1 = l*(l + 1)
                  DO m = -l, l
                     lm = ll1 + m
                     DO natom = 1, atoms%neq(ntyp)
                        sumab = CMPLX(0.0, 0.0)
                        DO n_r = 1, abc(ntyp)%n_r(l)
                           DO n_r2 = 1, abc(ntyp)%n_r(l)
                              sumab = sumab + abc(ntyp)%cof(i, lm, n_r, natom)* &
                                      CONJG(abc(ntyp)%cof(j, lm, n_r2, natom))* &
                                      radfun(ntyp)%integral(n_r, n_r2, l, 1, 1)
                           END DO
                        END DO
                        mmn(i, j) = mmn(i, j) + sumab
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL timestop("wannierlib_mmk0_sph")
   END SUBROUTINE wannierlib_mmk0_sph

END MODULE m_wannierlib_mmk0_sph
