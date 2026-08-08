!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Invariants that any matrix-element provider has to satisfy, checked on its result
!>  at one k-point. They know nothing about which operator it is, which is what makes
!>  them worth having: a new provider gets them for free.
!>
!>  Three checks, and each one catches a mistake we have actually made or nearly made:
!>
!>    finite    NaN or Inf in the result. Reading eigenvector storage that was never
!>              written gives harmless zeros in serial but stale window memory under
!>              MPI-RMA, and that arrives as garbage rather than as a crash.
!>    non-zero  a result that is identically zero. An operator whose slice was never
!>              filled, or which was gated behind a request nobody made, looks like a
!>              small number rather than like the absence of a calculation.
!>    hermitian |O - O^dagger| / max|O|. This is the one that pays for the module: get
!>              the band indices the wrong way round, or put CONJG on the wrong factor,
!>              and the matrix stops being Hermitian immediately -- while every other
!>              symptom of that mistake is a plausible-looking number.
!>
!>  Generalised from m_melem_spin's melem_spin_sumrule, which does the same thing for
!>  the spin operator alone and has been catching things there since it was written.
!>  What is NOT here is anything operator-specific: a sum rule, a bound on the
!>  diagonal, a comparison against a reference. Those belong with the operator.
!>
!>  It WARNS and does not abort. A tolerance that is right for one operator on one
!>  mesh is not obviously right for the next, and a check meant to help someone
!>  develop should not be able to kill a production run. Pass l_ok to a caller that
!>  wants to decide for itself.
MODULE m_melem_check
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_mat
   USE m_types_mpimat
   USE m_types_matelements
   IMPLICIT NONE
   PRIVATE

   !> Relative, and compared against max|O|, so it does not depend on the scale of the
   !> operator: a spin matrix is of order one and a Hamiltonian is not.
   REAL, PARAMETER :: MELEM_CHECK_TOL = 1.0e-10

   PUBLIC :: melem_check_provider, melem_check_matrix, MELEM_CHECK_TOL

CONTAINS

   !> The checks on a plain (nb,nb,ncomp) Bloch matrix, for a caller that holds one
   !> rather than a provider object -- an assembled coarse slice, say.
   SUBROUTINE melem_check_matrix(o, name, ik, tol, l_hermitian, l_ok)
      COMPLEX,           INTENT(IN)  :: o(:, :, :)   !> (nb,nb,ncomp)
      CHARACTER(LEN=*),  INTENT(IN)  :: name
      INTEGER,           INTENT(IN)  :: ik
      REAL,    OPTIONAL, INTENT(IN)  :: tol
      LOGICAL, OPTIONAL, INTENT(IN)  :: l_hermitian  !> default .TRUE.
      LOGICAL, OPTIONAL, INTENT(OUT) :: l_ok

      INTEGER :: nb, nc, i, j, a, n_bad
      REAL    :: t, amax, dev, rel
      LOGICAL :: lherm, ok
      COMPLEX :: z

      nb = SIZE(o, 1); nc = SIZE(o, 3)
      t = MELEM_CHECK_TOL; IF (PRESENT(tol)) t = tol
      lherm = .TRUE.;      IF (PRESENT(l_hermitian)) lherm = l_hermitian
      ok = .TRUE.

      IF (SIZE(o, 2) /= nb) CALL judft_bug("melem_check_matrix: the matrix is not square")

      WRITE(oUnit, '(a,i0,a)') 'wannierlib provider check ['//TRIM(name)//'], k = ', ik, ':'

      !> Written as "not (|x| <= HUGE)" and not as "|x| > HUGE" on purpose: every
      !> comparison with a NaN is false, so the negated form is the one that catches it.
      n_bad = 0
      DO a = 1, nc
         DO j = 1, nb
            DO i = 1, nb
               z = o(i, j, a)
               IF (.NOT. (ABS(REAL(z)) <= HUGE(1.0) .AND. ABS(AIMAG(z)) <= HUGE(1.0))) &
                  n_bad = n_bad + 1
            END DO
         END DO
      END DO

      IF (n_bad > 0) THEN
         !> Everything below would propagate the NaN into its own answer, so stop here.
         WRITE(oUnit, '(a,i0,a)') '   FAIL  finite     : ', n_bad, ' entries are NaN or Inf'
         WRITE(oUnit, '(a)')      '         a likely cause is reading states that were never stored'
         ok = .FALSE.
      ELSE
         WRITE(oUnit, '(a)') '   ok    finite'
         amax = MAXVAL(ABS(o))
         IF (amax <= 0.0) THEN
            WRITE(oUnit, '(a)') '   FAIL  non-zero   : every entry is zero -- was the result ever filled?'
            ok = .FALSE.
         ELSE
            WRITE(oUnit, '(a,es11.3)') '   ok    non-zero   : max|O| = ', amax
            IF (lherm) THEN
               dev = 0.0
               DO a = 1, nc
                  DO j = 1, nb
                     DO i = 1, nb
                        dev = MAX(dev, ABS(o(i, j, a) - CONJG(o(j, i, a))))
                     END DO
                  END DO
               END DO
               rel = dev/amax
               IF (rel > t) THEN
                  WRITE(oUnit, '(a,es11.3)') '   FAIL  hermitian  : max|O-O^dagger|/max|O| = ', rel
                  WRITE(oUnit, '(a)') '         the usual causes are a transposed band index (i,j), CONJG on the'
                  WRITE(oUnit, '(a)') '         wrong one of the two coefficients, and cross-spin blocks that are'
                  WRITE(oUnit, '(a)') '         not each other adjoint'
                  ok = .FALSE.
               ELSE
                  WRITE(oUnit, '(a,es11.3)') '   ok    hermitian  : ', rel
               END IF
            END IF
         END IF
      END IF

      IF (.NOT. ok) WRITE(oUnit, '(a)') &
         '   WARNING: this result does not satisfy the generic invariants'
      IF (PRESENT(l_ok)) l_ok = ok
   END SUBROUTINE melem_check_matrix

   !> The same checks on a provider object, which is the form a new operator has at hand.
   SUBROUTINE melem_check_provider(op, name, ik, tol, l_ok)
      CLASS(t_matelements), INTENT(IN)  :: op
      CHARACTER(LEN=*),     INTENT(IN)  :: name
      INTEGER,              INTENT(IN)  :: ik
      REAL,    OPTIONAL,    INTENT(IN)  :: tol
      LOGICAL, OPTIONAL,    INTENT(OUT) :: l_ok

      INTEGER :: nsp, nb, i, j
      REAL    :: t
      LOGICAL :: ok, ok_all
      COMPLEX, ALLOCATABLE :: z(:, :, :)

      t = MELEM_CHECK_TOL; IF (PRESENT(tol)) t = tol
      ok_all = .TRUE.

      !> The two stores are exclusive AS A RESULT: an operator that fills comp still has
      !> mat allocated, only never written, so checking mat would report an all-zero
      !> matrix for a perfectly good operator. comp wins whenever it exists.
      IF (ALLOCATED(op%comp)) THEN
         CALL melem_check_matrix(op%comp, name, ik, t, .TRUE., ok_all)
         IF (PRESENT(l_ok)) l_ok = ok_all
         RETURN
      END IF

      IF (.NOT. ALLOCATED(op%mat)) THEN
         WRITE(oUnit, '(a)') 'wannierlib provider check ['//TRIM(name)// &
                             ']: nothing computed yet -- skipped'
         IF (PRESENT(l_ok)) l_ok = .TRUE.
         RETURN
      END IF

      !> A distributed block holds only this rank's columns, so i and j are LOCAL indices
      !> and O(i,j) does not sit opposite O(j,i). The check would not be weaker here, it
      !> would be meaningless, so it is refused rather than reported.
      SELECT TYPE (blk => op%mat(1, 1))
      TYPE IS (t_mpimat)
         WRITE(oUnit, '(a)') 'wannierlib provider check ['//TRIM(name)// &
                             ']: the blocks are distributed -- skipped (these checks are local)'
         IF (PRESENT(l_ok)) l_ok = .TRUE.
         RETURN
      END SELECT

      !> The blocks are assembled into the single 2N spinor matrix and checked ONCE, which
      !> is not a shortcut but the only correct way to normalise. Per block it is wrong
      !> twice over: the cross-spin blocks are not Hermitian on their own -- they are each
      !> other's adjoint -- and, worse, dividing their deviation by their own magnitude
      !> divides roundoff by roundoff whenever the operator is nearly spin-diagonal, which
      !> is exactly the case for a magnet quantised along z. On the 2N matrix both
      !> statements are the same one, and the scale is the scale of the whole operator.
      nsp = SIZE(op%mat, 1)
      nb  = op%mat(1, 1)%matsize1
      ALLOCATE(z(nsp*nb, nsp*nb, 1))
      DO i = 1, nsp
         DO j = 1, nsp
            z((i-1)*nb + 1:i*nb, (j-1)*nb + 1:j*nb, 1) = block_complex(op%mat(i, j))
         END DO
      END DO
      CALL melem_check_matrix(z, name, ik, t, .TRUE., ok)
      ok_all = ok_all .AND. ok

      DEALLOCATE(z)
      IF (PRESENT(l_ok)) l_ok = ok_all
   END SUBROUTINE melem_check_provider

   !> A block as complex whatever it is stored as, so the checks are written once.
   FUNCTION block_complex(m) RESULT(z)
      CLASS(t_mat), INTENT(IN) :: m
      COMPLEX, ALLOCATABLE :: z(:, :)
      IF (m%l_real) THEN
         z = CMPLX(m%data_r, 0.0)
      ELSE
         z = m%data_c
      END IF
   END FUNCTION block_complex

END MODULE m_melem_check
