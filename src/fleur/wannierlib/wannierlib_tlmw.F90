!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The angular half of a projection: the trial orbital the user asked for, written in the
!>  basis of complex spherical harmonics the wavefunction is expanded in.
!>
!>  tlm holds the real harmonics of Wannier90's convention -- p_x = (Y_1,-1 - Y_1,1)/sqrt(2),
!>  p_y = i(Y_1,-1 + Y_1,1)/sqrt(2), and so on -- and tlmwf then picks, per Wannier function,
!>  the combination its (l, mr) names.
!>
!>  A NEGATIVE l is a hybrid, and that is where the coefficients stop being one row of tlm:
!>  sp3 mixes l=0 and l=1 with 1/sqrt(3) and 1/sqrt(6), sp3d2 reaches l=2 as well. Those are
!>  the guesses that localise on a bond rather than on an atom.
!>
!>  m_wannierlib_rad_twd carries the same table for the spin-orbit case, where the orbital is
!>  named by (j, m_j) and the coefficients are Clebsch-Gordan instead.
MODULE m_wannierlib_tlmw
  USE m_juDFT
  USE m_types_wannierlib
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_tlmw
CONTAINS

  SUBROUTINE wannierlib_tlmw(wannierlib, nwfs, l_spinors, jspin, tlmwf)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: wannierlib
    INTEGER, INTENT(IN) :: nwfs
    !> True whenever the run carries spinors (noco OR soc): the column guard below
    !> keeps each spinor component to the projections of its own channel, and without
    !> it both components fill all num_wann columns and amn comes out rank num_wann/2.
    !>
    !> Careful with the name l_nocosoc, which this routine deliberately does NOT take.
    !> It means opposite things in the two halves of the tree: (noco AND NOT soc) in
    !> wannierlib_main, but (noco OR soc) in wann_optional/wann_postproc -- and the
    !> latter is what FLEUR v26 calls l_nocosoc in wann_tlmw, the routine this one
    !> mirrors. So l_spinors here is the SAME condition as v26's guard, not a
    !> departure from it.
    LOGICAL, INTENT(IN) :: l_spinors
    INTEGER, INTENT(IN) :: jspin
    COMPLEX, INTENT(OUT) :: tlmwf(0:3, -3:3, nwfs)

    INTEGER :: nwf, lr, mr
    COMPLEX :: tlm(0:3, -3:3, 1:7)
    COMPLEX :: ic

    IF (nwfs > SIZE(wannierlib%proj)) THEN
      CALL juDFT_error('wannierlib_tlmw: nwfs exceeds configured projections', calledby='wannierlib_tlmw')
    END IF

    CALL timestart('wannierlib_tlmw')
    ic = CMPLX(0.0, 1.0)
    tlm = CMPLX(0.0, 0.0)

    !> The third index of tlm is Wannier90's mr, which counts the real harmonics of
    !> one l in a fixed order of its own: for l=1 that order is z, x, y; for l=2 it is
    !> z^2, xz, yz, x^2-y^2, xy; for l=3 it is z^3, xz^2, yz^2, z(x^2-y^2), xyz,
    !> x(x^2-3y^2), y(3x^2-y^2). Each row holds that one real harmonic written in the
    !> complex Y_lm, which is why the table is complex, and why a plain (l, mr)
    !> projection needs nothing more than the row.
    !>
    !> mr arrives already checked: projection_m_count in types_wannierlib refuses one
    !> outside the count of its l while the projection list is expanded, so no guard
    !> is repeated here and an all-zero row cannot be reached from the input.
    tlm(0, 0, 1) = CMPLX(1.0, 0.0)

    tlm(1, 0, 1) = CMPLX(1.0, 0.0)
    tlm(1, -1, 2) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(1, 1, 2) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(1, -1, 3) = ic / SQRT(2.0)
    tlm(1, 1, 3) = ic / SQRT(2.0)

    tlm(2, 0, 1) = CMPLX(1.0, 0.0)
    tlm(2, 1, 2) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, -1, 2) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, 1, 3) = ic / SQRT(2.0)
    tlm(2, -1, 3) = ic / SQRT(2.0)
    tlm(2, 2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, -2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, 2, 5) = -ic / SQRT(2.0)
    tlm(2, -2, 5) = ic / SQRT(2.0)

    tlm(3, 0, 1) = CMPLX(1.0, 0.0)
    tlm(3, 1, 2) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, -1, 2) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, 1, 3) = ic / SQRT(2.0)
    tlm(3, -1, 3) = ic / SQRT(2.0)
    tlm(3, 2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, -2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, 2, 5) = -ic / SQRT(2.0)
    tlm(3, -2, 5) = ic / SQRT(2.0)
    tlm(3, 3, 6) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, -3, 6) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, 3, 7) = ic / SQRT(2.0)
    tlm(3, -3, 7) = ic / SQRT(2.0)

    tlmwf = CMPLX(0.0, 0.0)

    DO nwf = 1, nwfs
      IF (l_spinors .AND. ((3 - 2 * jspin) /= wannierlib%proj(nwf)%spin)) CYCLE

      lr = wannierlib%proj(nwf)%l
      mr = wannierlib%proj(nwf)%m

      !> A negative l is a hybrid, and projection_m_count gives it 2, 3, 4, 5 or 6
      !> members: -1 sp, -2 sp2, -3 sp3, -4 sp3d, -5 sp3d2. They are the guesses that
      !> span more than one l, so each is written l by l instead of copied as a row.
      !>
      !> -4 repeating -2's first three is deliberate, not duplication: sp3d IS sp2
      !> plus two combinations that reach l=2, so its mr=1..3 are the sp2 ones.
      IF (lr >= 0) THEN
        tlmwf(lr, :, nwf) = tlm(lr, :, mr)
      ELSE IF (lr == -1) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(2.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 2)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(2.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 2)
        END IF
      ELSE IF (lr == -2) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) + (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) - (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (2.0 / SQRT(6.0)) * tlm(1, :, 2)
        END IF
      ELSE IF (lr == -3) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (tlm(1, :, 2) + tlm(1, :, 3) + tlm(1, :, 1))
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (tlm(1, :, 2) - tlm(1, :, 3) - tlm(1, :, 1))
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (-tlm(1, :, 2) + tlm(1, :, 3) - tlm(1, :, 1))
        ELSE IF (mr == 4) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (-tlm(1, :, 2) - tlm(1, :, 3) + tlm(1, :, 1))
        END IF
      ELSE IF (lr == -4) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) + (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) - (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (2.0 / SQRT(6.0)) * tlm(1, :, 2)
        ELSE IF (mr == 4) THEN
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(2.0)) * tlm(2, :, 1)
        ELSE IF (mr == 5) THEN
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(2.0)) * tlm(2, :, 1)
        END IF
      ELSE IF (lr == -5) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 2)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) + 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 2)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) + 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 3)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) - 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 4) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 3)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) - 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 5) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(3.0)) * tlm(2, :, 1)
        ELSE IF (mr == 6) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(3.0)) * tlm(2, :, 1)
        END IF
      ELSE
        CALL juDFT_error('no tlmw for this lr', calledby='wannierlib_tlmw')
      END IF
    END DO

    CALL timestop('wannierlib_tlmw')
  END SUBROUTINE wannierlib_tlmw

END MODULE m_wannierlib_tlmw
