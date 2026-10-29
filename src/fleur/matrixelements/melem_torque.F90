!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  What acts on the four exchange-field blocks once they exist: the Cartesian
!>  components of the exchange torque, and the check that they are perpendicular to
!>  the magnetisation. The blocks come from t_matelements_torque, which is the
!>  provider; this module is not one, in the same way m_melem_spin is not.
!>
!>  Freimuth, Mokrousov et al., PRB 90, 174423 (2014):
!>
!>      Omega^xc(r) = [V_min(r) - V_maj(r)] / (2 mu_B)      T(r) = -mu_B sigma x Omega^xc(r)
!>
!>  The provider already carries B_xc = (V_up - V_dn)/2 inside its blocks, so
!>  Omega = B_xc * m_hat and the mu_B cancels:  T_alpha = -(sigma_beta B_xc) x m_hat.
!>  With the moment along z that leaves T = (-<sigma_y B>, +<sigma_x B>, 0), which is
!>  why the classic route only ever writes the spin off-diagonal block as perpmag.
MODULE m_melem_torque
   USE m_judft
   USE m_constants, ONLY: oUnit
   USE m_types_atoms
   USE m_types_nococonv
   USE m_melem_spin, ONLY: melem_pauli_from_blocks
   USE m_types_kpts
   USE m_types_mpi
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_torque_axis, melem_torque_from_blocks, melem_torque_perp_check
CONTAINS

   !> The magnetisation direction shared by every site, in the global frame.
   !>
   !> The cross product is site-local: Omega points along the moment of the site the
   !> exchange field belongs to. Once the provider has summed the sites into one pair
   !> of blocks that information is gone, so the components below are only meaningful
   !> when all the sites agree. l_common says whether they do, and the caller is
   !> expected to refuse rather than to pick one axis and carry on.
   SUBROUTINE melem_torque_axis(atoms, nococonv, axis, l_common)
      TYPE(t_atoms),    INTENT(IN)  :: atoms
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      REAL,             INTENT(OUT) :: axis(3)
      LOGICAL,          INTENT(OUT) :: l_common

      REAL, PARAMETER :: tol = 1.0e-10
      INTEGER :: n
      REAL    :: a, b, v(3)

      a = nococonv%alph(1); b = nococonv%beta(1)
      axis = [SIN(b)*COS(a), SIN(b)*SIN(a), COS(b)]
      l_common = .TRUE.
      DO n = 2, atoms%ntype
         a = nococonv%alph(n); b = nococonv%beta(n)
         v = [SIN(b)*COS(a), SIN(b)*SIN(a), COS(b)]
         IF (MAXVAL(ABS(v - axis)) > tol) l_common = .FALSE.
      END DO
   END SUBROUTINE melem_torque_axis


   !> The three Cartesian components at one k from the four exchange-field blocks.
   !> The blocks must already be in the <i|O|j> convention and in the global frame,
   !> which is what the provider hands out.
   SUBROUTINE melem_torque_from_blocks(o_uu, o_dd, o_ud, o_du, axis, t0)
      COMPLEX, INTENT(IN)  :: o_uu(:,:), o_dd(:,:), o_ud(:,:), o_du(:,:)
      REAL,    INTENT(IN)  :: axis(3)     !> m_hat, global frame
      COMPLEX, INTENT(OUT) :: t0(:,:,:)   !> (nb,nb,3) 1=x 2=y 3=z

      COMPLEX, ALLOCATABLE :: s0(:,:,:)
      INTEGER :: nb

      nb = SIZE(o_uu, 1)
      IF (SIZE(t0,3) /= 3) CALL judft_bug("melem_torque_from_blocks: t0 needs three components")
      ALLOCATE(s0(nb, nb, 3))

      !> <sigma_alpha B_xc>: the same assembly the spin operator does, on blocks that
      !> carry the exchange field as their weight instead of a plain overlap.
      CALL melem_pauli_from_blocks(o_uu, o_dd, o_ud, o_du, s0)

      !> T = -(s x m_hat). With m_hat = z this is (-s_y, +s_x, 0).
      t0(:,:,1) = -(s0(:,:,2)*axis(3) - s0(:,:,3)*axis(2))
      t0(:,:,2) = -(s0(:,:,3)*axis(1) - s0(:,:,1)*axis(3))
      t0(:,:,3) = -(s0(:,:,1)*axis(2) - s0(:,:,2)*axis(1))

      DEALLOCATE(s0)
   END SUBROUTINE melem_torque_from_blocks


   !> T is a cross product with m_hat, so its projection on m_hat vanishes identically.
   !> That makes it a check with no reference attached: it fails on a sign error in the
   !> Pauli assembly, on blocks taken in the conjugate convention, and on an axis that
   !> does not match the one the blocks were rotated into.
   !>
   !> It cannot see an overall sign or scale, since both survive the projection. Those
   !> are fixed by comparing against the stored perpmag, not here.
   SUBROUTINE melem_torque_perp_check(t0, axis, ik, tol, l_ok)
      COMPLEX, INTENT(IN)  :: t0(:,:,:)
      REAL,    INTENT(IN)  :: axis(3)
      INTEGER, INTENT(IN)  :: ik
      REAL,    INTENT(IN)  :: tol
      LOGICAL, INTENT(OUT) :: l_ok

      COMPLEX, ALLOCATABLE :: proj(:,:)
      REAL :: rel, scale
      INTEGER :: nb

      nb = SIZE(t0, 1)
      ALLOCATE(proj(nb, nb))
      proj = t0(:,:,1)*axis(1) + t0(:,:,2)*axis(2) + t0(:,:,3)*axis(3)

      scale = MAX(MAXVAL(ABS(t0)), TINY(1.0))
      rel   = MAXVAL(ABS(proj)) / scale
      l_ok  = rel < tol

      WRITE(oUnit, '(a,i0,a,es10.3,a,a)') &
         'wannierlib exchange-torque perpendicularity, k = ', ik, &
         ': max |T.m| / max|T| = ', rel, '   ', MERGE('OK  ', 'FAIL', l_ok)
      IF (.NOT. l_ok) WRITE(oUnit, '(a)') &
         '  the torque is not perpendicular to the moment: the blocks, the Pauli '// &
         'assembly or the axis disagree'

      DEALLOCATE(proj)
   END SUBROUTINE melem_torque_perp_check

END MODULE m_melem_torque
