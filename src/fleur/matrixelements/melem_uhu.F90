!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The muffin-tin part of <u_{k+b1}|H_k|u_{k+b2}>, assembled from the three uHu passes and
!>  contracted against the matching coefficients of the two neighbours.
!>
!>  This is the direct route to the C of the modern theory. m_wannierlib_cf reaches the same
!>  element through an identity that resolves H against the bra and never applies it; the two
!>  are alternatives, and comparing them is the only check either one has.
!>
!>  WHAT IS LEFT OUT, and none of it is small: the interstitial, the vacuum, the non-spherical
!>  part of the potential, and the local orbitals. What comes back is a part of the element
!>  and not the element, which is why nothing here writes a uHu file. The local orbitals are
!>  absent rather than wrong: their slots stay at zero, because the spherical Hamiltonian
!>  needs their own linearisation energy, ello0, and the radial pass is given el0.
!>
!>  TWO b VECTORS, EACH MEASURED FROM THE PARENT POINT k, and they enter in three different
!>  ways: the modulus of each dresses its own side with j_L(|b| r), the direction of each
!>  enters through its own Y_lm, and only their difference reaches the atom-position phase.
!>  That difference is also what the overlap table is indexed by, which is what lets the phase
!>  and the coefficient contraction be taken from the overlap's own routine untouched.
MODULE m_melem_uhu
   USE m_juDFT
   USE m_constants, ONLY: fpi_const
   USE m_types_atoms
   USE m_types_cell
   USE m_types_abc
   USE m_types_enpara
   USE m_types_radfun
   USE m_melem_bgeom, ONLY: melem_blen, melem_bcart
   USE m_melem_uhu_dress, ONLY: melem_uhu_dress
   USE m_melem_uhu_gaunt, ONLY: melem_uhu_radpass, melem_uhu_angpass
   USE m_melem_mmkb_sph, ONLY: melem_mmkb_sph
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_uhu_states
CONTAINS

   !> Accumulates onto uhu, the way the overlap accumulates onto its own matrix.
   SUBROUTINE melem_uhu_states(atoms, cell, enpara, radfun, jspin_rad, abc_a, abc_b, &
                               b1frac, b2frac, bkpt_a, bkpt_b, gpar, asym_max, uhu)
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER, INTENT(IN)        :: jspin_rad   !> radial set, and the potential with it
      TYPE(t_abc), INTENT(IN)    :: abc_a(:)    !> the bra, at k+b1
      TYPE(t_abc), INTENT(IN)    :: abc_b(:)    !> the ket, at k+b2
      REAL, INTENT(IN)           :: b1frac(3)   !> b1 measured from k, internal coordinates
      REAL, INTENT(IN)           :: b2frac(3)   !> b2 measured from k, internal coordinates
      REAL, INTENT(IN)           :: bkpt_a(3), bkpt_b(3)
      INTEGER, INTENT(IN)        :: gpar(3)     !> difference of the two folding vectors
      REAL, INTENT(OUT)          :: asym_max    !> largest R_MT surface term met on the way
      COMPLEX, INTENT(INOUT)     :: uhu(:, :)

      !> u and udot. A local orbital would need ello0, so its slot is left out of the pass
      !> rather than filled with the wrong energy.
      INTEGER, PARAMETER :: nslot_apw = 2

      COMPLEX, ALLOCATABLE :: uhug(:, :, :, :, :, :)
      REAL, ALLOCATABLE :: uj_a(:, :, :, :, :), duj_a(:, :, :, :, :)
      REAL, ALLOCATABLE :: uj_b(:, :, :, :, :), duj_b(:, :, :, :, :)
      REAL, ALLOCATABLE :: uhur(:, :, :, :, :, :, :)
      REAL :: rk_a, rk_b, bc_a(3), bc_b(3), kdiff(3, 1)
      INTEGER :: n, lwn, lmd, total_nr, nr

      CALL timestart("melem_uhu_states")
      lmd = atoms%lmaxd*(atoms%lmaxd + 2)
      total_nr = 0
      DO n = 1, atoms%ntype
         total_nr = MAX(total_nr, MAXVAL(radfun(n)%n_r(:)))
      END DO
      nr = MIN(nslot_apw, total_nr)
      asym_max = 0.0

      rk_a = melem_blen(cell, b1frac)
      rk_b = melem_blen(cell, b2frac)
      bc_a = melem_bcart(cell, b1frac)
      bc_b = melem_bcart(cell, b2frac)

      ALLOCATE (uj_a(atoms%jmtd, 2, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd), &
                duj_a(atoms%jmtd, 2, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd), &
                uj_b(atoms%jmtd, 2, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd), &
                duj_b(atoms%jmtd, 2, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd))
      ALLOCATE (uhur(total_nr, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:atoms%lmaxd, &
                     0:atoms%lmaxd, 0:atoms%lmaxd))
      ALLOCATE (uhug(0:lmd, 0:lmd, total_nr, total_nr, atoms%ntype, 1), &
                source=CMPLX(0.0, 0.0))

      DO n = 1, atoms%ntype
         lwn = atoms%lmax(n)
         CALL melem_uhu_dress(atoms, n, rk_a, radfun(n), jspin_rad, lwn, uj_a, duj_a)
         CALL melem_uhu_dress(atoms, n, rk_b, radfun(n), jspin_rad, lwn, uj_b, duj_b)

         !> Handed the APW slots alone, which is what keeps the local orbitals out: the pass
         !> takes the extent of its slot loop from what it is given.
         uhur = 0.0
         CALL melem_uhu_radpass(atoms, n, lwn, enpara%el0(:, n, jspin_rad), &
                                enpara%el0(:, n, jspin_rad), enpara%vr(:, n, jspin_rad), &
                                uj_a(:, :, 1:nr, :, :), duj_a(:, :, 1:nr, :, :), &
                                uj_b(:, :, 1:nr, :, :), duj_b(:, :, 1:nr, :, :), &
                                uhur, asym_max)
         CALL melem_uhu_angpass(atoms, n, lwn, bc_a, bc_b, uhur, 1, uhug)
      END DO

      !> ONE MORE 4pi THAN AN OVERLAP CARRIES. m_melem_mmkb_sph supplies the single factor
      !> that expanding one plane wave brings, because that is all an overlap has; uHu
      !> expands one on each side, and the second factor belongs to the table.
      uhug = fpi_const*uhug

      !> One pair per table, so the lookup by value has a single entry to find and the
      !> difference it is keyed on is the one the phase needs.
      kdiff(:, 1) = bkpt_b + gpar - bkpt_a
      CALL melem_mmkb_sph(atoms, abc_a, abc_b, bkpt_b, gpar, bkpt_a, uhug, kdiff, 1, uhu)

      DEALLOCATE (uj_a, duj_a, uj_b, duj_b, uhur, uhug)
      CALL timestop("melem_uhu_states")
   END SUBROUTINE melem_uhu_states

END MODULE m_melem_uhu
