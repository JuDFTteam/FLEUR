!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_kptgen_hybrid

   USE m_juDFT

CONTAINS

! this programm generates an aequdistant kpoint set including the
! Gamma point; it is reduced to IBZ and written in kpts (M.B.)

! odified for types D.W.

   SUBROUTINE kptgen_hybrid(film, grid, cell, sym, kpts, l_soc, l_onlyIdentitySym)

      USE m_types_cell
      USE m_types_sym
      USE m_types_kpts
      USE m_constants

      IMPLICIT NONE

      LOGICAL, INTENT(IN)           :: film
      INTEGER, INTENT(IN)           :: grid(3)
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_kpts), INTENT(INOUT)   :: kpts
! - scalars -
      LOGICAL, INTENT(IN)   ::  l_soc
      LOGICAL, INTENT(IN)   ::  l_OnlyIdentitySym
! - local scalars -
      INTEGER                  ::  i, j, k, nkpt
      INTEGER                  ::  ikpt, ikpt0, nkpti
      INTEGER                  ::  nsym, nop
! - local arrays -
      INTEGER, ALLOCATABLE   ::  rot(:, :, :), rrot(:, :, :)
      INTEGER, ALLOCATABLE   ::  neqkpt(:)
      INTEGER, ALLOCATABLE   ::  pkpt(:, :, :), kptp(:), symkpt(:), iarr(:), &
                                iarr2(:)
      REAL, ALLOCATABLE      ::  rtau(:, :)
      REAL, ALLOCATABLE      ::  bk(:, :), bkhlp(:, :)
      REAL, ALLOCATABLE      ::  rarr(:)
      LOGICAL               ::  ldum

      nkpt = grid(1)*grid(2)*grid(3)
      ALLOCATE(bk(3, nkpt), bkhlp(3, nkpt))

      ikpt = 0
      DO i = 0, grid(1) - 1
         DO j = 0, grid(2) - 1
            DO k = 0, grid(3) - 1
               ikpt = ikpt + 1
               bk(:, ikpt) = [1.0*i/grid(1), &
                              1.0*j/grid(2), &
                              1.0*k/grid(3)]
            END DO
         END DO
      END DO

      IF(ikpt /= nkpt) call judft_error( 'failure: number of k-points')

      IF(sym%invs .OR. l_soc) THEN
         nsym = sym%nop
      ELSE
         nsym = 2*sym%nop
      END IF
      nop = sym%nop
      IF(l_OnlyIdentitySym) THEN
         nop = 1
         nsym = 1
      END IF

      ALLOCATE(rot(3, 3, nsym), rtau(3, nsym))

      DO i = 1, nop
         rot(:, :, i) = sym%mrot(:, :, i)
         rtau(:, i) = sym%tau(:, i)
      END DO

      DO i = nop + 1, nsym
         rot(:, :, i) = rot(:, :, i - nop)
         rtau(:, i) = rtau(:, i - nop)
      END DO

      IF(any(rot(:, :, 1) - reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/),(/3, 3/)) /= 0)) &
         call judft_error( 'kptgen: First symmetry operation is not the identity.')

      ALLOCATE(rrot(3, 3, nsym))

      DO i = 1, nop
         rrot(:, :, i) = transpose(rot(:, :, sym%invtab(i)))
      END DO

      DO i = nop + 1, nsym
         rrot(:, :, i) = -rrot(:, :, i - nop)
      END DO

      ALLOCATE(kptp(nkpt), symkpt(nkpt), rarr(3), iarr2(3), iarr(nkpt))
      ALLOCATE(pkpt(grid(1) + 1, grid(2) + 1, grid(3) + 1))
      pkpt = 0
      DO ikpt = 1, nkpt
         iarr2 = nint(bk(:, ikpt)*grid) + 1
         pkpt(iarr2(1), iarr2(2), iarr2(3)) = ikpt
      END DO

      pkpt(grid(1) + 1, :, :) = pkpt(1, :, :)
      pkpt(:, grid(2) + 1, :) = pkpt(:, 1, :)
      pkpt(:, :, grid(3) + 1) = pkpt(:, :, 1)

      IF(any(pkpt == 0)) THEN
         CALL juDFT_error('kptgen: Definition of pkpt-pointer failed.', &
                          calledby='kptgen_hybrid')
      END IF
      iarr = 1
      ldum = .FALSE.
      DO i = 1, nkpt
         IF(iarr(i) == 0) CYCLE
         kptp(i) = i
         symkpt(i) = 1
         DO k = 2, nsym
            rarr = matmul(rrot(:, :, k), bk(:, i))*grid
            iarr2 = nint(rarr)
            IF(any(abs(iarr2 - rarr) > 1e-10)) THEN
               WRITE(oUnit, '(A,I3,A)') 'kptgen: Symmetry operation', k, &
                  ' incompatible with k-point set.'
               ldum = .TRUE.
            END IF
            iarr2 = modulo(iarr2, grid) + 1
            IF(any(iarr2 > grid)) &
               call judft_error( 'kptgen: pointer indices exceed pointer dimensions.')
            j = pkpt(iarr2(1), iarr2(2), iarr2(3))
            IF(j == 0) call judft_error( 'kptgen: k-point index is zero (bug?)')
            IF(iarr(j) == 0 .OR. j == i) CYCLE
            iarr(j) = 0
            kptp(j) = i
            symkpt(j) = k
         END DO
      END DO
      IF(ldum) &
         call judft_error( 'kptgen: Some symmetry operations are incompatible &
         with k-point set.')
      i = 0
      DO ikpt = 1, nkpt
         IF(iarr(ikpt) == 1) THEN
            i = i + 1
            iarr(ikpt) = i
         END IF
      END DO
      nkpti = i
      DO ikpt = 1, nkpt
         IF(iarr(ikpt) == 0) THEN
            i = i + 1
            iarr(ikpt) = i
         END IF
      END DO
      bk(:, iarr) = bk
      kptp = iarr(kptp)
      kptp(iarr) = kptp
      symkpt(iarr) = symkpt
      DO i = 1, grid(1) + 1
         DO j = 1, grid(2) + 1
            DO k = 1, grid(3) + 1
               pkpt(i, j, k) = iarr(pkpt(i, j, k))
            END DO
         END DO
      END DO

      ALLOCATE(neqkpt(nkpti))
      neqkpt = 0
      DO ikpt0 = 1, nkpti
         DO ikpt = 1, nkpt
            IF(kptp(ikpt) == ikpt0) neqkpt(ikpt0) = neqkpt(ikpt0) + 1
         END DO
      END DO

!     Do not do any IO, but store in kpts
      kpts%nkpt = nkpti
      if(allocated(kpts%bk)) deallocate(kpts%bk)
      if(allocated(kpts%wtkpt)) deallocate(kpts%wtkpt)
      ALLOCATE(kpts%bk(3, kpts%nkpt), kpts%wtkpt(kpts%nkpt))

      DO ikpt = 1, nkpti
         kpts%bk(:, ikpt) = bk(:, ikpt)
         kpts%wtkpt(ikpt) = neqkpt(ikpt)
      END DO

   CONTAINS

! Returns least common multiple of the integers iarr(1:n).
      FUNCTION kgv(iarr, n)
         IMPLICIT NONE
         INTEGER                 :: kgv
         INTEGER, INTENT(IN)  :: n, iarr(n)
         LOGICAL              :: lprim(2:maxval(iarr))
         INTEGER, ALLOCATABLE :: prim(:), expo(:)
         INTEGER                 :: nprim, marr
         INTEGER                 :: i, j, ia, k
! Determine prime numbers
         marr = maxval(iarr)
         lprim = .TRUE.
         DO i = 2, marr
            j = 2
            DO WHILE(i*j <= marr)
               lprim(i*j) = .FALSE.
               j = j + 1
            END DO
         END DO
         nprim = count(lprim)
         ALLOCATE(prim(nprim), expo(nprim))
         j = 0
         DO i = 2, marr
            IF(lprim(i)) THEN
               j = j + 1
               prim(j) = i
            END IF
         END DO
! Determine least common multiple
         expo = 0
         DO i = 1, n
            ia = iarr(i)
            IF(ia == 0) CYCLE
            DO j = 1, nprim
               k = 0
               DO WHILE(ia/prim(j)*prim(j) == ia)
                  k = k + 1
                  ia = ia/prim(j)
               END DO
               expo(j) = max(expo(j), k)
            END DO
         END DO
         kgv = 1
         DO j = 1, nprim
            kgv = kgv*prim(j)**expo(j)
         END DO
         DEALLOCATE(prim, expo)
      END FUNCTION kgv

!     function modulo1 maps kpoint into first BZ
      FUNCTION modulo1(kpoint, nkpt, a, b, c)

         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: nkpt, a, b, c
         REAL, INTENT(IN)    :: kpoint(3)
         REAL                   :: modulo1(3)
         INTEGER                :: help(3), nkpt3(3)

         nkpt3 = (/a, b, c/)
         modulo1 = kpoint*nkpt3
         help = nint(modulo1)
         IF(any(abs(help - modulo1) > 1e-8)) THEN
            modulo1 = kpoint*nkpt3
            WRITE(*, *) modulo1
            help = nint(modulo1)
            WRITE(*, *) help
            WRITE(oUnit, '(A,F5.3,2('','',F5.3),A)') 'modulo1: argument (', &
               kpoint, ') is not an element of the k-point set.'
            CALL juDFT_error( &
               'modulo1: argument not an element of k-point set.', &
               calledby='kptgen_hybrid:modulo1')
         END IF
         modulo1 = modulo(help, nkpt3)*1.0/nkpt3

      END FUNCTION modulo1

   END SUBROUTINE kptgen_hybrid

END MODULE m_kptgen_hybrid
