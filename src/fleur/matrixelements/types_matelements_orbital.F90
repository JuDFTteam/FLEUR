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
!>  block, and both spinor flags stay .FALSE.
!>
!>  The coefficients handed in are either a spinor, whose two components are summed
!>  here, or a single collinear channel, which is one component and stands alone.
!>  Which of the two is declared at init: a radial slot is given for a channel and
!>  withheld for a spinor, whose components read the slots of their own potentials.
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
   USE m_types_spinor_layout, ONLY: radial_slot
   USE m_types_usdus
   USE m_types_atoms
   USE m_constants, ONLY: ImagUnit
   USE m_judft
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matelements) :: t_matelements_orbital
      TYPE(t_atoms), POINTER :: atoms => NULL()
      INTEGER :: ntyp  = 0     !> atom type
      INTEGER :: iat   = 0     !> which of the equivalent atoms of that type
      INTEGER :: jspin_rad = 0 !> radial slot of the single channel served; 0 = spinor
   CONTAINS
      PROCEDURE :: init                 => init
      PROCEDURE :: calc_matrix_elements => calc_matrix_elements
   END TYPE t_matelements_orbital

   PUBLIC :: t_matelements_orbital

CONTAINS

   !> Bind the object to one atom. All three Cartesian components come out of one pass,
   !> since L_x and L_y are both built from the same L_+ and L_-. The site is not summed over
   !> here: an instance covers one site, and the caller adds them up if it wants the
   !> total. For L that is a plain sum, since it is spin-diagonal and its trace is
   !> invariant under the local-to-global frame rotation.
   SUBROUTINE init(this, atoms, ntyp, iat, jspin_rad)
      CLASS(t_matelements_orbital), INTENT(INOUT) :: this
      TYPE(t_atoms), TARGET, INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: ntyp, iat
      INTEGER, INTENT(IN), OPTIONAL :: jspin_rad   !> present: the coefficients are one channel

      IF (ntyp < 1 .OR. ntyp > atoms%ntype) &
         CALL judft_bug("init: the atom type is out of range")
      IF (iat < 1 .OR. iat > atoms%neq(ntyp)) &
         CALL judft_bug("init: the equivalent atom is out of range")

      this%jspin_rad = 0
      IF (PRESENT(jspin_rad)) THEN
         IF (jspin_rad < 1) CALL judft_bug("init: the radial slot is out of range")
         this%jspin_rad = jspin_rad
      END IF

      this%spinoroperator = .FALSE.
      this%spinorwavefcts = .FALSE.

      this%atoms => atoms
      this%ntyp  = ntyp
      this%iat   = iat
   END SUBROUTINE init

   SUBROUTINE calc_matrix_elements(this, zmat, abc, radfun, usdus)
      CLASS(t_matelements_orbital), INTENT(INOUT) :: this
      TYPE(t_mat),    INTENT(IN) :: zmat(:)   !> unused, L works on the abc coefficients only
      TYPE(t_abc),    INTENT(IN) :: abc(:, :) !> (2 spin, ntype) local-frame coefficients
      TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype)
      TYPE(t_usdus),  INTENT(IN) :: usdus     !> unused, the radial integrals are in radfun

      INTEGER :: nb, i, j, l, ll1, mm, lm, n_r, n_r2, s, sr, n_comp, slot(2)
      REAL    :: lplus, lminus, w
      COMPLEX :: cz, cp, cm     ! L_z, L_+ and L_- for one (i,j) of this atom

      IF (.NOT.ALLOCATED(this%mat)) &
         CALL judft_bug("calc_matrix_elements: the result matrix is not allocated")
      IF (SIZE(this%mat, 1) /= 1 .OR. SIZE(this%mat, 2) /= 1) &
         CALL judft_bug("calc_matrix_elements: L is spin-diagonal, so the result is a single block")
      n_comp = SIZE(abc, 1)
      IF (n_comp /= 1 .AND. n_comp /= 2) &
         CALL judft_bug("calc_matrix_elements: the coefficients are either one channel or a spinor")
      IF (SIZE(abc, 2) /= this%atoms%ntype) &
         CALL judft_bug("calc_matrix_elements: the coefficients must be indexed by atom type")
      IF ((this%jspin_rad /= 0) .NEQV. (n_comp == 1)) &
         CALL judft_bug("calc_matrix_elements: a radial slot belongs to a channel and not to a spinor")

      nb = SIZE(abc(1, 1)%cof, 1)
      IF (this%mat(1, 1)%matsize1 /= nb) &
         CALL judft_bug("calc_matrix_elements: the matrix size does not match the abc coefficients")

      !> radfun%integral is allocated (.,.,.,jspins,jspins), so with a single radial
      !> set both spinor components have to read slot 1; indexing 2 there ran past
      !> the array. The bound is read from the array itself.
      sr = radial_slot(radfun, 2)
      IF (n_comp == 1) THEN
         slot(1) = this%jspin_rad          ! the channel names its own potential
      ELSE
         slot(1) = 1; slot(2) = sr         ! a spinor reads both, clamped to what exists
      END IF

      DO j = 1, nb                     ! ket band
         DO i = 1, nb                  ! bra band
            cz = CMPLX(0.0, 0.0); cp = CMPLX(0.0, 0.0); cm = CMPLX(0.0, 0.0)
            DO l = 0, this%atoms%lmax(this%ntyp)
               ll1 = l*(l + 1)
               DO mm = -l, l
                  lm = ll1 + mm
                  lplus  = SQRT(REAL((l - mm)*(l + mm + 1)))   ! <m+1|L+|m>
                  lminus = SQRT(REAL((l + mm)*(l - mm + 1)))   ! <m-1|L-|m>
                  DO s = 1, n_comp     ! L is spin-diagonal: a spinor sums, a channel stands alone
                     DO n_r = 1, abc(s, this%ntyp)%n_r(l)
                        DO n_r2 = 1, abc(s, this%ntyp)%n_r(l)
                           w = radfun(this%ntyp)%integral(n_r, n_r2, l, slot(s), slot(s))
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
            this%comp(i, j, 1, 1, 1) = 0.5 * (cp + cm)               ! L_x = (L+ + L-)/2
            this%comp(i, j, 1, 1, 2) = -0.5 * ImagUnit * (cp - cm)   ! L_y = (L+ - L-)/(2i)
            this%comp(i, j, 1, 1, 3) = cz                            ! L_z
         END DO
      END DO
   END SUBROUTINE calc_matrix_elements

END MODULE m_types_matelements_orbital
