!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Thin d-shell adapter for the Nordstrom/Bultmark multipole transformation.
!
! rho_physical is the physical ket-bra density
!
!   rho(m,s,m',s') = C(m,s) CONJG(C(m',s')),
!
! with m=-2,...,+2 and spin order s=1 (+1/2), s=2 (-1/2).  The caller chooses
! the frame: this adapter performs no orbital or spin rotation.  For the usual
! coupled-tensor interpretation, both sectors should already use one common
! coordinate frame.  A DFT+U mmpMat is not valid direct input to this routine.
!--------------------------------------------------------------------------------
MODULE m_local_d_spin_multipoles
   USE m_spin_orbital_multipoles, ONLY: spin_orbital_density_to_multipoles
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: d_l = 2

   PUBLIC :: local_d_spin_density_to_multipoles

CONTAINS

   PURE SUBROUTINE local_d_spin_density_to_multipoles(rho_physical, w)
      COMPLEX, INTENT(IN) :: rho_physical(-d_l:, 1:, -d_l:, 1:)
      COMPLEX, INTENT(OUT) :: w(0:, 0:, 0:, -(2*d_l+1):)

      IF (ANY(SHAPE(rho_physical) /= [5, 2, 5, 2])) THEN
         ERROR STOP "local_d_spin_density_to_multipoles: rho must have shape (5,2,5,2)"
      END IF
      IF (ANY(SHAPE(w) /= [5, 2, 6, 11])) THEN
         ERROR STOP "local_d_spin_density_to_multipoles: w must have shape (5,2,6,11)"
      END IF

      CALL spin_orbital_density_to_multipoles(d_l, rho_physical, w)
   END SUBROUTINE local_d_spin_density_to_multipoles

END MODULE m_local_d_spin_multipoles
