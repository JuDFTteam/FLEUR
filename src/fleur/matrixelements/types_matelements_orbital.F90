!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  One Cartesian component of the orbital angular momentum at one atom, in the
!>  Bloch basis:
!>
!>      L_alpha,mn = < psi_m | L_alpha | psi_n >_MT ,   alpha = x,y,z
!>
!>  L acts on the spatial part alone, so it is the identity in spin space: the
!>  cross-spin blocks vanish and the expectation value over a spinor is the sum of
!>  its two components. The result is therefore a single matrix and not a 2x2
!>  block, both spinor flags stay .FALSE., and the sum over the spinor components
!>  is taken here.
!>
!>  L_x and L_y are both built from L_+ and L_-, so an instance per component
!>  repeats the radial and angular sums rather than sharing them.
!>
!>  Evaluated in the complex spherical harmonics, where L_z is diagonal (= m) and
!>  L_+ / L_- connect m -> m+1 / m -> m-1 with sqrt((l-+m)(l+-m+1)). The radial
!>  part is the plain overlap, since L does not act on it. Interstitial and vacuum
!>  are left out: the orbital moment is a muffin-tin quantity.
MODULE m_types_matelements_orbital
   USE m_types_matelements
   USE m_types_mat
   USE m_types_abc
   USE m_types_radfun
   USE m_types_usdus
   USE m_types_atoms
   USE m_constants, ONLY: ImagUnit
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matelements) :: t_matelements_orbital
      TYPE(t_atoms), POINTER :: atoms => NULL()
      INTEGER :: compo = 0     !> which component: 1 = L_x, 2 = L_y, 3 = L_z
      INTEGER :: ntyp  = 0     !> atom type
      INTEGER :: iat   = 0     !> which of the equivalent atoms of that type
   CONTAINS
      PROCEDURE :: init                 => init
      PROCEDURE :: calc_matrix_elements => calc_matrix_elements
   END TYPE t_matelements_orbital

   PUBLIC :: t_matelements_orbital

CONTAINS

   !> Bind the object to one component of one atom. The site is not summed over
   !> here: an instance covers one site, and the caller adds them up if it wants the
   !> total. For L that is a plain sum, since it is spin-diagonal and its trace is
   !> invariant under the local-to-global frame rotation.
   SUBROUTINE init(this, atoms, compo, ntyp, iat)
      CLASS(t_matelements_orbital), INTENT(INOUT) :: this
      TYPE(t_atoms), TARGET, INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: compo, ntyp, iat

      IF (compo < 1 .OR. compo > 3) &
         CALL judft_bug("init: the component index must be 1, 2 or 3")
      IF (ntyp < 1 .OR. ntyp > atoms%ntype) &
         CALL judft_bug("init: the atom type is out of range")
      IF (iat < 1 .OR. iat > atoms%neq(ntyp)) &
         CALL judft_bug("init: the equivalent atom is out of range")

      this%spinoroperator = .FALSE.
      this%spinorwavefcts = .FALSE.

      this%atoms => atoms
      this%compo = compo
      this%ntyp  = ntyp
      this%iat   = iat
   END SUBROUTINE init

   SUBROUTINE calc_matrix_elements(this, zmat, abc, radfun, usdus)
      CLASS(t_matelements_orbital), INTENT(INOUT) :: this
      TYPE(t_mat),    INTENT(IN) :: zmat(:)   !> unused, L works on the abc coefficients only
      TYPE(t_abc),    INTENT(IN) :: abc(:, :) !> (2 spin, ntype) local-frame coefficients
      TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype)
      TYPE(t_usdus),  INTENT(IN) :: usdus     !> unused, the radial integrals are in radfun

      INTEGER :: nb, i, j, l, ll1, mm, lm, n_r, n_r2, s, sr
      REAL    :: lplus, lminus, w
      COMPLEX :: cz, cp, cm     ! L_z, L_+ and L_- for one (i,j) of this atom

      IF (.NOT.ALLOCATED(this%mat)) &
         CALL judft_bug("calc_matrix_elements: the result matrix is not allocated")
      IF (SIZE(this%mat, 1) /= 1 .OR. SIZE(this%mat, 2) /= 1) &
         CALL judft_bug("calc_matrix_elements: L is spin-diagonal, so the result is a single block")
      IF (SIZE(abc, 1) /= 2 .OR. SIZE(abc, 2) /= this%atoms%ntype) &
         CALL judft_bug("calc_matrix_elements: the abc coefficients must have shape (2,ntype)")

      nb = SIZE(abc(1, 1)%cof, 1)
      IF (this%mat(1, 1)%matsize1 /= nb) &
         CALL judft_bug("calc_matrix_elements: the matrix size does not match the abc coefficients")

      !> radfun%integral is allocated (.,.,.,jspins,jspins), so with a single radial
      !> set both spinor components have to read slot 1; indexing 2 there ran past
      !> the array. The bound is read from the array itself.
      sr = MERGE(1, 2, SIZE(radfun(1)%integral, 4) < 2)

      DO j = 1, nb                     ! ket band
         DO i = 1, nb                  ! bra band
            cz = CMPLX(0.0, 0.0); cp = CMPLX(0.0, 0.0); cm = CMPLX(0.0, 0.0)
            DO l = 0, this%atoms%lmax(this%ntyp)
               ll1 = l*(l + 1)
               DO mm = -l, l
                  lm = ll1 + mm
                  lplus  = SQRT(REAL((l - mm)*(l + mm + 1)))   ! <m+1|L+|m>
                  lminus = SQRT(REAL((l + mm)*(l - mm + 1)))   ! <m-1|L-|m>
                  DO s = 1, 2          ! L is spin-diagonal: sum both spinor components
                     DO n_r = 1, abc(s, this%ntyp)%n_r(l)
                        DO n_r2 = 1, abc(s, this%ntyp)%n_r(l)
                           w = radfun(this%ntyp)%integral(n_r, n_r2, l, MIN(s, sr), MIN(s, sr))
                           cz = cz + abc(s, this%ntyp)%cof(i, lm, n_r, this%iat) &
                                   * CONJG(abc(s, this%ntyp)%cof(j, lm, n_r2, this%iat))*REAL(mm)*w
                           IF (mm < l) cp = cp + abc(s, this%ntyp)%cof(i, lm, n_r, this%iat) &
                                   * CONJG(abc(s, this%ntyp)%cof(j, lm + 1, n_r2, this%iat))*lplus*w
                           IF (mm > -l) cm = cm + abc(s, this%ntyp)%cof(i, lm, n_r, this%iat) &
                                   * CONJG(abc(s, this%ntyp)%cof(j, lm - 1, n_r2, this%iat))*lminus*w
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            !> Assigned and not accumulated: one instance covers one component of one
            !> atom, so each element gets exactly one contribution. Adding onto a
            !> cleared matrix would also turn a -0.0 result into +0.0, and the sign of
            !> zero is visible in the exported file.
            SELECT CASE (this%compo)
            CASE (1)
               this%mat(1, 1)%data_c(i, j) = 0.5 * (cp + cm)               ! L_x = (L+ + L-)/2
            CASE (2)
               this%mat(1, 1)%data_c(i, j) = -0.5 * ImagUnit * (cp - cm)   ! L_y = (L+ - L-)/(2i)
            CASE (3)
               this%mat(1, 1)%data_c(i, j) = cz                            ! L_z
            END SELECT
         END DO
      END DO
   END SUBROUTINE calc_matrix_elements

END MODULE m_types_matelements_orbital
