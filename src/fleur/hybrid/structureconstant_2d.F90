!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>2D lattice sums S_lm = sum_T exp(2 pi i k.T) conjg(Y_lm(T+r0)) |T+r0|^(-l-1) for films,
!>same layout as structureconstant.  l >= L_SPLIT: direct sum; l < L_SPLIT: Ewald-type method.
!>At Gamma the divergent l = 0 head is omitted (m_gamma_2d); only its finite remainder is kept.
MODULE m_structureconstant_2d
   USE m_juDFT
   USE m_constants
   USE m_ylm
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_structconst_2d_parry, ONLY: parry_structconst_2d
   USE m_structconst_2d_weinert, ONLY: weinert_structconst_2d
   USE m_types_atoms
   USE m_types_cell
   USE m_types_hybinp
   USE m_types_kpts
   USE m_types_mpi
   IMPLICIT NONE
   PRIVATE

   ! low-l method: Parry (production), Weinert pseudocharge or direct sum (cross-checks)
   INTEGER, PARAMETER, PUBLIC :: STRUCTCONST_2D_PARRY = 1
   INTEGER, PARAMETER, PUBLIC :: STRUCTCONST_2D_WEINERT = 2
   INTEGER, PARAMETER, PUBLIC :: STRUCTCONST_2D_BRUTE = 3

   INTEGER, PARAMETER, PUBLIC :: structconst_2d_method = STRUCTCONST_2D_PARRY

   INTEGER, PARAMETER, PUBLIC :: L_SPLIT = 8

   INTEGER, PARAMETER :: NSHELL_REF = 400

   PUBLIC :: structconst_2d_direct, structconst_2d_lowl
   PUBLIC :: structureconstant_2d

   REAL, PARAMETER, PUBLIC :: CONVTOL_2D = 1e-18

CONTAINS

   SUBROUTINE structconst_2d_direct(lmax, qrel, r0, amat, nshell, conv_tol, s)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: lmax, nshell
      REAL, INTENT(IN)     :: qrel(3), r0(3), amat(3, 3), conv_tol
      COMPLEX, INTENT(OUT) :: s((lmax + 1)**2)

      INTEGER :: n1, n2, l, m, lm, ish, maxl
      INTEGER :: conv(0:lmax)
      REAL    :: rv(3), a, ph, apow, shell_max(0:lmax)
      COMPLEX :: y((lmax + 1)**2), cexp, cdum

      s = cmplx_0
      conv = HUGE(1)
      maxl = lmax

      DO ish = 0, nshell
         shell_max = 0.0
         DO n1 = -ish, ish
            DO n2 = -ish, ish
               IF (MAX(ABS(n1), ABS(n2)) /= ish) CYCLE

               rv = n1*amat(:, 1) + n2*amat(:, 2) + r0
               a = NORM2(rv)
               IF (a < 1e-12) CYCLE

               ph = qrel(1)*n1 + qrel(2)*n2
               cexp = EXP(ImagUnit*tpi_const*ph)

               CALL ylm4(maxl, rv, y)
               y = CONJG(y)

               DO l = 0, maxl
                  IF (ish > conv(l)) CYCLE
                  apow = a**(-l - 1)
                  cdum = cexp*apow
                  shell_max(l) = MAX(shell_max(l), apow)
                  lm = l**2
                  DO m = -l, l
                     lm = lm + 1
                     s(lm) = s(lm) + cdum*y(lm)
                  END DO
               END DO
            END DO
         END DO

         IF (ish >= 1) THEN
            DO l = 0, lmax
               IF (conv(l) == HUGE(1) .AND. shell_max(l) < conv_tol) conv(l) = ish
            END DO
            DO WHILE (maxl > 0)
               IF (conv(maxl) == HUGE(1)) EXIT
               maxl = maxl - 1
            END DO
            IF (ALL(conv /= HUGE(1))) EXIT
         END IF
      END DO
   END SUBROUTINE structconst_2d_direct

   SUBROUTINE structconst_2d_lowl(lmax, qrel, r0, amat, bmat, conv_tol, s)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: lmax
      REAL, INTENT(IN)     :: qrel(3), r0(3), amat(3, 3), bmat(3, 3), conv_tol
      COMPLEX, INTENT(OUT) :: s((lmax + 1)**2)

      SELECT CASE (structconst_2d_method)
      CASE (STRUCTCONST_2D_PARRY)
         CALL parry_structconst_2d(lmax, qrel, r0, amat, bmat, -1.0, s)
      CASE (STRUCTCONST_2D_BRUTE)
         CALL structconst_2d_direct(lmax, qrel, r0, amat, NSHELL_REF, conv_tol, s)
      CASE (STRUCTCONST_2D_WEINERT)
         CALL weinert_structconst_2d(lmax, qrel, r0, amat, bmat, -1, s)
      CASE DEFAULT
         CALL juDFT_error("unknown structconst_2d_method", calledby="structconst_2d_lowl")
      END SELECT
   END SUBROUTINE structconst_2d_lowl

   SUBROUTINE structureconstant_2d(structconst, cell, hybinp, atoms, kpts, fmpi)
      IMPLICIT NONE
      COMPLEX, INTENT(INOUT)     :: structconst(:, :, :, :)
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_mpi), INTENT(IN)    :: fmpi

      INTEGER :: lmax, ic1, ic2, ikpt, l, lm, i
      REAL    :: r0(3), qrel(3)
      INTEGER :: ierr, buf_sz, root
      COMPLEX :: s((2*hybinp%lexp + 1)**2), s_low((2*hybinp%lexp + 1)**2)

      CALL timestart("calc struc_const. 2D")
      lmax = 2*hybinp%lexp
      structconst = cmplx_0

      DO ic2 = 1 + fmpi%irank, atoms%nat, fmpi%isize
         !$OMP PARALLEL DO default(none) schedule(dynamic) collapse(2) &
         !$OMP private(ic1, ikpt, r0, qrel, s, s_low, l, lm) &
         !$OMP shared(atoms, cell, kpts, structconst, lmax, ic2)
         DO ic1 = 1, atoms%nat
            DO ikpt = 1, kpts%nkpt
               IF (ic2 == 1 .OR. ic1 /= ic2) THEN
                  r0 = MATMUL(cell%amat, atoms%taual(:, ic2) - atoms%taual(:, ic1))
                  qrel = kpts%bk(:, ikpt)

                  CALL structconst_2d_direct(lmax, qrel, r0, cell%amat, 200, CONVTOL_2D, s)
                  CALL structconst_2d_lowl(L_SPLIT - 1, qrel, r0, cell%amat, cell%bmat, &
                                           CONVTOL_2D, s_low)

                  DO l = L_SPLIT, lmax
                     DO lm = l**2 + 1, (l + 1)**2
                        structconst(lm, ic1, ic2, ikpt) = s(lm)
                     END DO
                  END DO
                  DO l = 0, MIN(L_SPLIT - 1, lmax)
                     DO lm = l**2 + 1, (l + 1)**2
                        structconst(lm, ic1, ic2, ikpt) = s_low(lm)
                     END DO
                  END DO
               END IF
            END DO
         END DO
         !$OMP END PARALLEL DO
      END DO

#ifdef CPP_MPI
      CALL timestart("bcast struc_const. 2D")
      buf_sz = SIZE(structconst, 1)*SIZE(structconst, 2)
      DO ic2 = 1, atoms%nat
         root = MOD(ic2 - 1, fmpi%isize)
         DO ikpt = 1, kpts%nkpt
            CALL MPI_Bcast(structconst(1, 1, ic2, ikpt), buf_sz, MPI_DOUBLE_COMPLEX, &
                           root, fmpi%mpi_comm, ierr)
         END DO
      END DO
      CALL timestop("bcast struc_const. 2D")
#endif

      DO i = 2, atoms%nat
         structconst(:, i, i, :) = structconst(:, 1, 1, :)
      END DO

      IF (fmpi%irank == 0) WRITE (oUnit, '(/,A)') "### subroutine: structureconstant_2d ###"

      CALL timestop("calc struc_const. 2D")
   END SUBROUTINE structureconstant_2d

END MODULE m_structureconstant_2d
