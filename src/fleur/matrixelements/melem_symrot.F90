!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The rotation that takes one atom's muffin-tin frame to the global one, in cartesian
!>  coordinates.
!>
!>  WHY ANY OF THIS. abcof expands the muffin-tin coefficients in the LOCAL frame of each
!>  atom: it rotates k+G by TRANSPOSE(mrot(:,:,invtab(ngopr(na)))) before expanding, so that
!>  every atom of a type shares one set of radial functions and one table. For a SCALAR
!>  operator the frame cancels in the product A* A and there is nothing to do, which is why
!>  the pair overlap closes to machine zero with equivalent atoms. For a VECTOR operator it
!>  does not: its components live in that local frame, and summing sites without bringing
!>  them to a common frame adds vectors expressed in different axes.
!>
!>  AND THE CLASS MATTERS AS MUCH AS THE MATRIX, because the momentum and the orbital moment
!>  do not transform the same way:
!>
!>      p' = R p                 (polar)
!>      L' = det(R) R L          (axial)
!>
!>  Under inversion det(R) = -1, and that factor is exactly what leaves L untouched while it
!>  flips p. Applying R to both alike would break L, which is correct today. Hence the two
!>  entry points rather than one.
MODULE m_melem_symrot
   USE m_types_sym
   USE m_types_cell
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_rot_polar, melem_rot_axial

CONTAINS

   !>  R for a polar vector (the momentum): local frame of atom na -> global, cartesian.
   FUNCTION melem_rot_polar(sym, cell, na) RESULT(r)
      TYPE(t_sym),  INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER,      INTENT(IN) :: na
      REAL :: r(3, 3)

      !> mrot lives in LATTICE coordinates and amat carries the lattice vectors as COLUMNS,
      !> so a cartesian vector transforms with amat . mrot . amat^-1.
      !>
      !> invtab(ngopr), NOT ngopr, and this was MEASURED rather than derived. abcof and
      !> types_abc build the coefficients by going global -> local with
      !> TRANSPOSE(mrot(:,:,invtab(ngopr))), so coming back is the inverse operation, which
      !> is what invtab names. The matrix is amat and not bmat: k+G is a reciprocal vector
      !> and converts with bmat, but what rotates here are the cartesian components of a
      !> gradient, built on real-space spherical harmonics.
      !>
      !> The four spellings -- ngopr or invtab(ngopr), transposed or not -- were swept on a
      !> case of THREE atoms related by a proper rotation, and only this one closes:
      !>
      !>     mrot(ngopr)             3.77e-1  3.46e-1  1.62e-6
      !>     mrot(invtab(ngopr))     2.17e-6  1.88e-6  1.62e-6   <-- this one
      !>     TRANSPOSE(mrot(ngopr))  1.25e-1  2.76e-1  1.62e-6
      !>     TRANSPOSE(mrot(invtab)) 2.53e-1  6.61e-1  1.62e-6
      !>
      !> The sweep NEEDED three atoms. Under inversion ngopr == invtab and R = -I is its own
      !> inverse, so all four agree and the choice is invisible -- a wrong one passes every
      !> test a two-atom cell can offer. The z column above says the same from the other
      !> side: C3 turns about z, every spelling leaves that component alone, and all four
      !> agree there to the last digit.
      r = MATMUL(cell%amat, MATMUL(REAL(sym%mrot(:, :, sym%invtab(sym%ngopr(na)))), &
                                   inverse3(cell%amat)))
   END FUNCTION melem_rot_polar

   !>  det(R) R for an axial vector (the orbital moment).
   FUNCTION melem_rot_axial(sym, cell, na) RESULT(r)
      TYPE(t_sym),  INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER,      INTENT(IN) :: na
      REAL :: r(3, 3)

      r = melem_rot_polar(sym, cell, na)
      r = det3(r)*r
   END FUNCTION melem_rot_axial

   PURE REAL FUNCTION det3(a)
      REAL, INTENT(IN) :: a(3, 3)
      det3 = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
           - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
           + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
   END FUNCTION det3

   PURE FUNCTION inverse3(a) RESULT(b)
      REAL, INTENT(IN) :: a(3, 3)
      REAL :: b(3, 3), d
      d = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
      b(1,1) =  (a(2,2)*a(3,3) - a(2,3)*a(3,2))/d
      b(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))/d
      b(1,3) =  (a(1,2)*a(2,3) - a(1,3)*a(2,2))/d
      b(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))/d
      b(2,2) =  (a(1,1)*a(3,3) - a(1,3)*a(3,1))/d
      b(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))/d
      b(3,1) =  (a(2,1)*a(3,2) - a(2,2)*a(3,1))/d
      b(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))/d
      b(3,3) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1))/d
   END FUNCTION inverse3

END MODULE m_melem_symrot
