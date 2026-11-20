!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_spin_sym
   !! Action of the crystal symmetry operations on spin (axial vector) quantities.
   USE m_constants
   USE m_juDFT
   USE m_types_cell
   USE m_types_sym
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: axial_rot_cart, inplane_spin_sym_factor

CONTAINS

   FUNCTION axial_rot_cart(sym, cell, iop) RESULT(p)
      !! Cartesian matrix det(R)*R of operation iop, i.e. its action on an axial vector.
      TYPE(t_sym),  INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER,      INTENT(IN) :: iop
      REAL :: p(3,3)

      INTEGER :: m(3,3), det

      m = sym%mrot(:, :, iop)
      det = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
          - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
          + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))
      ! cell%bmat = tpi_const * amat^{-1}
      p = det*MATMUL(cell%amat, MATMUL(REAL(m), cell%bmat/tpi_const))
   END FUNCTION axial_rot_cart

   REAL FUNCTION inplane_spin_sym_factor(sym, cell)
      !! Factor by which the star (lattice harmonic) coefficients of the in-plane
      !! magnetization (m_x,m_y) have to be scaled when spin and lattice are coupled (SOC).
      !! The scalar star projection averages over all operations without rotating the spin;
      !! the correct projection averages det(R)*R over the group. Its in-plane block is the
      !! identity if no operation acts on (m_x,m_y) (only E and I) and vanishes otherwise.
      TYPE(t_sym),  INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell

      REAL, PARAMETER :: tol = 1.0e-6
      INTEGER :: iop
      REAL    :: p(3,3), pavg(2,2), unit2(2,2)

      pavg = 0.0
      DO iop = 1, sym%nop
         p = axial_rot_cart(sym, cell, iop)
         pavg = pavg + p(1:2, 1:2)
      END DO
      pavg = pavg/sym%nop

      unit2 = RESHAPE([1.0, 0.0, 0.0, 1.0], [2,2])
      IF (MAXVAL(ABS(pavg)) < tol) THEN
         inplane_spin_sym_factor = 0.0
      ELSE IF (MAXVAL(ABS(pavg - unit2)) < tol) THEN
         inplane_spin_sym_factor = 1.0
      ELSE
         inplane_spin_sym_factor = 1.0
         CALL juDFT_warn("Symmetry operations incompatible with a noncollinear SOC calculation", &
                          hint="Each operation must leave the z axis invariant as an axial vector. "// &
                               "Recreate the input with inpgen while l_soc=T is set, or use '-nosym'.", &
                          calledby="inplane_spin_sym_factor")
      END IF
   END FUNCTION inplane_spin_sym_factor

END MODULE m_spin_sym
