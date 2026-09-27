!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>  The geometry of a neighbour vector b, in the two forms the tabulated matrix elements
!>  need it: its length, which selects the radial table, and its Cartesian direction, which
!>  selects the spherical harmonics.
!>
!>  THE ORDER OF THE CONTRACTION. `cell%bmat` holds the reciprocal vectors as ROWS, so
!>
!>      fractional -> Cartesian     MATMUL(v, bmat)        contracts the vector index
!>      Cartesian  -> fractional    MATMUL(bmat, r)/2pi    contracts the row index
!>
!>  and both orders are correct, for opposite directions. They are also interchangeable
!>  whenever `amat` is symmetric, which every cubic cell is, so a cubic case cannot tell a
!>  mistake here from a right answer. That is what the single place is for.
MODULE m_melem_bgeom
   USE m_types_cell
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_blen, melem_bcart

CONTAINS

   !> |b|, through the reciprocal metric. Callers group vectors on the exact bit pattern
   !> of this number, so it has to be computed the same way everywhere or two tables that
   !> should share a modulus stop sharing it.
   PURE REAL FUNCTION melem_blen(cell, bfrac) RESULT(rk)
      TYPE(t_cell), INTENT(IN) :: cell
      REAL, INTENT(IN) :: bfrac(3)
      rk = SQRT(DOT_PRODUCT(bfrac, MATMUL(cell%bbmat, bfrac)))
   END FUNCTION melem_blen

   !> b in Cartesian coordinates, which is what ylm4 expects.
   PURE FUNCTION melem_bcart(cell, bfrac) RESULT(bcart)
      TYPE(t_cell), INTENT(IN) :: cell
      REAL, INTENT(IN) :: bfrac(3)
      REAL :: bcart(3)
      bcart = MATMUL(bfrac, cell%bmat)
   END FUNCTION melem_bcart

END MODULE m_melem_bgeom
