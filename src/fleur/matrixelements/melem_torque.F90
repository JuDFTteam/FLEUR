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
   USE m_wannierlib_export_paths, ONLY: wannierlib_export_file
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_torque_axis, melem_torque_from_blocks, melem_torque_perp_check
   PUBLIC :: melem_torque_write_updown
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

   !> The cross-spin exchange-field block over the whole mesh, in the Bloch basis and in
   !> the format the classic route writes as updown.perpmag: a title, nbnd/nkpt/nbnd, and
   !> one free-format record per element.
   !>
   !> Written before any gauge on purpose. updown.perpmag is a Bloch-basis quantity -- it
   !> is what wann_updown produces and what its consumers read -- so rotating it first
   !> would make it a different object that happens to share a name.
   !>
   !> Index order follows wann_perpmag.f: the FIRST index is the bra and carries spin up,
   !> the second is the ket and carries spin down.
   SUBROUTINE melem_torque_write_updown(tx0, kpts, distk, fmpi)
      USE m_constants, ONLY: oUnit
#ifdef CPP_MPI
      USE mpi
#endif
      !> (nbnd, nbnd, nk_loc) on this rank, in ascending global-k order.
      COMPLEX, INTENT(IN) :: tx0(:, :, :)
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi

      INTEGER :: nb, nkl, kl, gk, i, j, k, iu, ierr
      COMPLEX, ALLOCATABLE :: full(:, :, :)

      nb = SIZE(tx0, 1)
      nkl = COUNT(distk == fmpi%irank)
      IF (SIZE(tx0, 3) < nkl) CALL juDFT_error( &
         'melem_torque_write_updown: fewer blocks than k-points on this rank', &
         calledby='melem_torque_write_updown')

      !> Every rank holds whole matrices for its own k and zeros elsewhere, so the mesh is
      !> assembled with one sum rather than a gather with displacements.
      ALLOCATE (full(nb, nb, kpts%nkptf), source=CMPLX(0.0, 0.0))
      kl = 0
      DO i = 1, SIZE(distk)
         IF (distk(i) /= fmpi%irank) CYCLE
         kl = kl + 1
         full(:, :, i) = tx0(:, :, kl)
      END DO
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, full, SIZE(full), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                         fmpi%mpi_comm, ierr)
#endif

      IF (fmpi%irank == 0) THEN
         OPEN (newunit=iu, file=wannierlib_export_file('updown.perpmag'), status='replace', &
               recl=1000)
         WRITE (iu, *) 'Overlaps of the wavefunct. with Beff'
         WRITE (iu, '(i5,i7,i5)') nb, kpts%nkptf, nb
         DO k = 1, kpts%nkptf
            DO j = 1, nb
               DO i = 1, nb
                  WRITE (iu, *) i, j, k, REAL(full(i, j, k)), AIMAG(full(i, j, k))
               END DO
            END DO
         END DO
         CLOSE (iu)
         WRITE (oUnit, '(a,i0,a,i0,a)') &
            'wannierlib: wrote updown.perpmag (exchange torque, Bloch basis, ', nb, &
            ' bands, ', kpts%nkptf, ' k-points)'
      END IF
      DEALLOCATE (full)
   END SUBROUTINE melem_torque_write_updown

END MODULE m_melem_torque
