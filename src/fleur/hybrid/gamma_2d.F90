!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>q -> 0 treatment for films: weight of the 2 pi/(A q) head missing from the k-mesh sum.
!>The finite g = 0 remainder, -2 pi |z-z'|, is part of the Coulomb matrix itself.
MODULE m_gamma_2d
   USE m_juDFT
   USE m_constants
   USE m_types_cell
   USE m_types_kpts

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: divergence_2d, divergence_2d_core

CONTAINS

   SUBROUTINE divergence_2d(cell, kpts, divergence)
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      REAL, INTENT(OUT)        :: divergence

      REAL :: nkpt3(3)

      nkpt3 = kpts%calcNkpt3()
      divergence = divergence_2d_core(cell%omtil/cell%amat(3, 3), cell%bmat(1, :), &
                                      cell%bmat(2, :), NINT(nkpt3(1)), NINT(nkpt3(2)))
   END SUBROUTINE divergence_2d

   !>D = (A/(2 pi)^2) int_BZ d^2q/q - (1/N_k) sum_{q/=0} 1/q.  The BZ integral is done exactly
   !>in polar form, int R(phi) dphi, with R the distance to the Wigner-Seitz boundary.
   REAL FUNCTION divergence_2d_core(area, b1, b2, n1, n2)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: area, b1(3), b2(3)
      INTEGER, INTENT(IN) :: n1, n2

      INTEGER, PARAMETER :: NPHI = 20000
      INTEGER :: ip, i, j, di, dj, nk2
      REAL    :: phi, u(2), rr, gg(2), gn2, gdu, acc_i, acc_s, q(2), qn, qbest

      CALL timestart("divergence_2d")

      acc_i = 0.0
      DO ip = 1, NPHI
         phi = tpi_const*(ip - 0.5)/NPHI
         u = [COS(phi), SIN(phi)]
         rr = HUGE(1.0)
         DO i = -2, 2
            DO j = -2, 2
               IF (i == 0 .AND. j == 0) CYCLE
               gg = [i*b1(1) + j*b2(1), i*b1(2) + j*b2(2)]
               gn2 = gg(1)**2 + gg(2)**2
               gdu = gg(1)*u(1) + gg(2)*u(2)
               IF (gdu > 1e-12) rr = MIN(rr, 0.5*gn2/gdu)
            END DO
         END DO
         acc_i = acc_i + rr
      END DO
      acc_i = acc_i*tpi_const/NPHI
      acc_i = area/(tpi_const*tpi_const)*acc_i

      nk2 = MAX(n1*n2, 1)
      acc_s = 0.0
      DO i = 0, n1 - 1
         DO j = 0, n2 - 1
            IF (i == 0 .AND. j == 0) CYCLE
            qbest = HUGE(1.0)
            DO di = -1, 1
               DO dj = -1, 1
                  q = [(REAL(i)/n1 + di)*b1(1) + (REAL(j)/n2 + dj)*b2(1), &
                       (REAL(i)/n1 + di)*b1(2) + (REAL(j)/n2 + dj)*b2(2)]
                  qn = SQRT(q(1)**2 + q(2)**2)
                  IF (qn < qbest) qbest = qn
               END DO
            END DO
            IF (qbest > 1e-12) acc_s = acc_s + 1.0/qbest
         END DO
      END DO
      acc_s = acc_s/nk2

      divergence_2d_core = acc_i - acc_s
      CALL timestop("divergence_2d")
   END FUNCTION divergence_2d_core

END MODULE m_gamma_2d
