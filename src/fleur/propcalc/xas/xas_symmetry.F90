!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_symmetry
   USE m_constants, ONLY: tpi_const
   USE m_dwigner, ONLY: d_wigner
   USE m_grp_k, ONLY: euler
   USE m_juDFT, ONLY: juDFT_error
   USE m_types_abc, ONLY: t_abc
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_cell, ONLY: t_cell
   USE m_types_kpts, ONLY: t_kpts
   USE m_types_nococonv, ONLY: t_nococonv
   USE m_types_sym, ONLY: t_sym
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: xas_count_star_members
   PUBLIC :: xas_star_member_weight
   PUBLIC :: xas_star_operation
   PUBLIC :: xas_cart_rotation_from_sym
   PUBLIC :: xas_rotate_lab_polarization_for_parent
   PUBLIC :: xas_rotate_abc_star_member
   PUBLIC :: xas_rotate_abc_star_member_spinor
   PUBLIC :: xas_su2_from_sym
   PUBLIC :: xas_local_spin_transform

CONTAINS

   SUBROUTINE xas_count_star_members(kpts, ikpt, nstar)
      TYPE(t_kpts), INTENT(IN)  :: kpts
      INTEGER,      INTENT(IN)  :: ikpt
      INTEGER,      INTENT(OUT) :: nstar

      IF (ikpt < 1 .OR. ikpt > kpts%nkpt) THEN
         CALL juDFT_error("Invalid parent k-point index in xas_count_star_members", calledby="m_xas_symmetry")
      END IF
      IF (.NOT. ALLOCATED(kpts%bkp)) THEN
         CALL juDFT_error("kpts%bkp is not allocated in xas_count_star_members", calledby="m_xas_symmetry")
      END IF

      nstar = COUNT(kpts%bkp(:) == ikpt)
   END SUBROUTINE xas_count_star_members

   SUBROUTINE xas_star_member_weight(kpts, ikpt, wk_star)
      TYPE(t_kpts), INTENT(IN)  :: kpts
      INTEGER,      INTENT(IN)  :: ikpt
      REAL,         INTENT(OUT) :: wk_star

      INTEGER :: nstar

      CALL xas_count_star_members(kpts, ikpt, nstar)
      IF (nstar <= 0) THEN
         CALL juDFT_error("No full-zone star members for parent k-point in xas_star_member_weight", &
                          calledby="m_xas_symmetry")
      END IF

      wk_star = kpts%wtkpt(ikpt)/REAL(nstar)
   END SUBROUTINE xas_star_member_weight

   SUBROUTINE xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
      TYPE(t_sym), INTENT(IN)  :: sym
      INTEGER,     INTENT(IN)  :: bksym
      INTEGER,     INTENT(OUT) :: spatial_iop
      LOGICAL,     INTENT(OUT) :: l_time_reversal

      IF (bksym <= sym%nop) THEN
         spatial_iop = bksym
         l_time_reversal = .FALSE.
      ELSE
         spatial_iop = bksym - sym%nop
         l_time_reversal = .TRUE.
      END IF

      IF (spatial_iop < 1 .OR. spatial_iop > sym%nop) THEN
         CALL juDFT_error("Invalid spatial symmetry index in xas_star_operation", calledby="m_xas_symmetry")
      END IF
   END SUBROUTINE xas_star_operation

   SUBROUTINE xas_cart_rotation_from_sym(sym, cell, spatial_iop, r_cart)
      TYPE(t_sym),  INTENT(IN)  :: sym
      TYPE(t_cell), INTENT(IN)  :: cell
      INTEGER,      INTENT(IN)  :: spatial_iop
      REAL,         INTENT(OUT) :: r_cart(3, 3)

      IF (spatial_iop < 1 .OR. spatial_iop > sym%nop) THEN
         CALL juDFT_error("Invalid spatial symmetry index in xas_cart_rotation_from_sym", &
                          calledby="m_xas_symmetry")
      END IF

      ! Same direct-lattice to Cartesian vector rotation convention as rotate_forces:
      ! R_cart = A * M_internal * A^{-1}, with A^{-1} = bmat / 2pi in FLEUR.
      r_cart = MATMUL(cell%amat, MATMUL(REAL(sym%mrot(:, :, spatial_iop)), cell%bmat/tpi_const))
   END SUBROUTINE xas_cart_rotation_from_sym

   SUBROUTINE xas_rotate_lab_polarization_for_parent(sym, cell, bksym, eps_lab_cart, eps_parent_cart)
      TYPE(t_sym),  INTENT(IN)  :: sym
      TYPE(t_cell), INTENT(IN)  :: cell
      INTEGER,      INTENT(IN)  :: bksym
      COMPLEX,      INTENT(IN)  :: eps_lab_cart(3)
      COMPLEX,      INTENT(OUT) :: eps_parent_cart(3)

      REAL    :: r_cart(3, 3)
      COMPLEX :: eps_source(3)
      INTEGER :: spatial_iop
      LOGICAL :: l_time_reversal

      CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
      CALL xas_cart_rotation_from_sym(sym, cell, spatial_iop, r_cart)

      eps_source = eps_lab_cart
      IF (l_time_reversal) eps_source = CONJG(eps_source)

      eps_parent_cart = MATMUL(CMPLX(TRANSPOSE(r_cart), 0.0), eps_source)
   END SUBROUTINE xas_rotate_lab_polarization_for_parent

   SUBROUTINE xas_su2_from_sym(sym, cell, spatial_iop, su_global)
      TYPE(t_sym),  INTENT(IN)  :: sym
      TYPE(t_cell), INTENT(IN)  :: cell
      INTEGER,      INTENT(IN)  :: spatial_iop
      COMPLEX,      INTENT(OUT) :: su_global(2, 2)

      REAL :: alpha, beta, gamma

      IF (spatial_iop < 1 .OR. spatial_iop > sym%nop) THEN
         CALL juDFT_error("Invalid spatial symmetry index in xas_su2_from_sym", calledby="m_xas_symmetry")
      END IF

      ! Match the spin-1/2 rotation convention used by grp_k/sympsi.
      ! The Euler helper takes the proper part of the spatial operation.
      ! su_global maps parent global spinor components to the symmetry-
      ! related star-member global spinor components for this operation.
      CALL euler(spatial_iop, sym, cell, alpha, beta, gamma)
      su_global(1, 1) = COS(beta/2.0)*EXP(CMPLX(0.0, -(alpha + gamma)/2.0))
      su_global(1, 2) = -SIN(beta/2.0)*EXP(CMPLX(0.0, -(alpha - gamma)/2.0))
      su_global(2, 1) = SIN(beta/2.0)*EXP(CMPLX(0.0,  (alpha - gamma)/2.0))
      su_global(2, 2) = COS(beta/2.0)*EXP(CMPLX(0.0,  (alpha + gamma)/2.0))
   END SUBROUTINE xas_su2_from_sym

   SUBROUTINE xas_local_spin_transform(atoms, nococonv, parent_atom, mapped_atom, su_global, su_local)
      TYPE(t_atoms),    INTENT(IN)  :: atoms
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      INTEGER,          INTENT(IN)  :: parent_atom, mapped_atom
      COMPLEX,          INTENT(IN)  :: su_global(2, 2)
      COMPLEX,          INTENT(OUT) :: su_local(2, 2)

      COMPLEX :: u_parent(2, 2), u_mapped(2, 2)

      IF (parent_atom < 1 .OR. parent_atom > atoms%nat) THEN
         CALL juDFT_error("Invalid parent atom index in xas_local_spin_transform", calledby="m_xas_symmetry")
      END IF
      IF (mapped_atom < 1 .OR. mapped_atom > atoms%nat) THEN
         CALL juDFT_error("Invalid mapped atom index in xas_local_spin_transform", calledby="m_xas_symmetry")
      END IF

      ! abc%calc_abc forms local-spin coefficients as
      !   abc_local(sigma) = sum_s conjg(U_global_local(s,sigma)) * abc_global(s),
      ! i.e. abc_local = U^\dagger abc_global. Therefore the parent-local to
      ! mapped-local spin transformation is U_mapped^\dagger S_global U_parent.
      u_parent = nococonv%umat(atoms%itype(parent_atom))
      u_mapped = nococonv%umat(atoms%itype(mapped_atom))
      su_local = MATMUL(CONJG(TRANSPOSE(u_mapped)), MATMUL(su_global, u_parent))
   END SUBROUTINE xas_local_spin_transform

   SUBROUTINE xas_rotate_abc_star_member(abc_parent, atoms, sym, cell, itype, bksym, lmax, abc_star)
      ! Transform MT expansion coefficients from the parent k point to one
      ! full-zone star member. The m-block convention follows waveftrafo_symm:
      !
      !   A_star(l,m) = sum_mprime A_parent(l,mprime) * D_l(mprime,m,iop)
      !
      ! For time-reversal-added operations the parent coefficients are
      ! conjugated before the Wigner-D multiplication. Non-symmorphic phase
      ! factors are not included here; for an incoherent local single-absorber
      ! intensity they are common to all l channels of one atom and cancel in
      ! |M|^2. They must be restored before using this helper for coherent
      ! multi-atom amplitudes.
      !
      ! The current hardwired XAS validation sums atoms explicitly. This helper
      ! keeps the t_abc local-atom indexing of the selected atom type and checks
      ! that the symmetry operation maps selected atoms inside the same atom
      ! type. General production code should use the mapped atom index to place
      ! the rotated coefficients in the mapped absorber channel.
      TYPE(t_abc),   INTENT(IN)  :: abc_parent
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_sym),   INTENT(IN)  :: sym
      TYPE(t_cell),  INTENT(IN)  :: cell
      INTEGER,       INTENT(IN)  :: itype, bksym, lmax
      TYPE(t_abc),   INTENT(OUT) :: abc_star

      INTEGER :: spatial_iop, l, iOrd, band, iAtom_l, iAtom, mapped_atom, mapped_l
      INTEGER :: lm1, lm2
      INTEGER :: mrot_one(3, 3, 1)
      LOGICAL :: l_time_reversal
      COMPLEX :: source(2*lmax + 1)
      COMPLEX, ALLOCATABLE :: d_wgn(:, :, :, :)

      IF (.NOT. ALLOCATED(abc_parent%cof)) THEN
         CALL juDFT_error("abc_parent%cof is not allocated in xas_rotate_abc_star_member", &
                          calledby="m_xas_symmetry")
      END IF
      CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
      abc_star = abc_parent
      IF (lmax >= 1) THEN
         ALLOCATE(d_wgn(-lmax:lmax, -lmax:lmax, 1:lmax, 1))
         mrot_one(:, :, 1) = sym%mrot(:, :, spatial_iop)
         CALL d_wigner(1, mrot_one, cell%bmat, lmax, d_wgn)
      END IF

      DO iAtom_l = LBOUND(abc_parent%cof, 4), UBOUND(abc_parent%cof, 4)
         iAtom = atoms%firstAtom(itype) + iAtom_l - 1
         IF (ALLOCATED(sym%mapped_atom)) THEN
            mapped_atom = sym%mapped_atom(spatial_iop, iAtom)
            IF (mapped_atom < atoms%firstAtom(itype) .OR. mapped_atom >= atoms%firstAtom(itype) + atoms%neq(itype)) THEN
               CALL juDFT_error("XAS star mode 3 currently requires mapped absorber atoms to stay in the same atom type", &
                                calledby="m_xas_symmetry")
            END IF
            mapped_l = mapped_atom - atoms%firstAtom(itype) + 1
            IF (mapped_l /= iAtom_l) THEN
               CALL juDFT_error("XAS star mode 3 currently assumes local absorbing atom maps to itself", &
                                calledby="m_xas_symmetry")
            END IF
         END IF

         IF (l_time_reversal) THEN
            abc_star%cof(:, 0, :, iAtom_l) = CONJG(abc_parent%cof(:, 0, :, iAtom_l))
         END IF

         DO l = 1, lmax
            lm1 = l*l
            lm2 = l*(l + 2)
            IF (lm2 > UBOUND(abc_parent%cof, 2)) CYCLE
            DO iOrd = 1, MIN(abc_parent%n_r(l), SIZE(abc_parent%cof, 3))
               DO band = 1, SIZE(abc_parent%cof, 1)
                  source(1:2*l + 1) = abc_parent%cof(band, lm1:lm2, iOrd, iAtom_l)
                  IF (l_time_reversal) source(1:2*l + 1) = CONJG(source(1:2*l + 1))
                  abc_star%cof(band, lm1:lm2, iOrd, iAtom_l) = &
                     MATMUL(source(1:2*l + 1), d_wgn(-l:l, -l:l, l, 1))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE xas_rotate_abc_star_member

   SUBROUTINE xas_rotate_abc_star_member_spinor(abc_parent_spin, atoms, sym, cell, nococonv, itype, bksym, lmax, abc_star_spin)
      ! Spinor star transform for unitary spatial operations only:
      !
      !   A_star(m,s) = sum_mprime,sum_sprime
      !      A_parent(mprime,sprime) * D_l(mprime,m;iop) * S_local(s,sprime)
      !
      ! The orbital D_l ordering is the same as xas_rotate_abc_star_member.
      ! S_local = U_mapped^\dagger S_global U_parent follows abc%calc_abc,
      ! where local coefficients are U^\dagger times global spinors.
      ! Time-reversal-related star members are antiunitary and are deliberately
      ! guarded until the spin flip/conjugation convention is implemented.
      TYPE(t_abc),       INTENT(IN)  :: abc_parent_spin(:)
      TYPE(t_atoms),     INTENT(IN)  :: atoms
      TYPE(t_sym),       INTENT(IN)  :: sym
      TYPE(t_cell),      INTENT(IN)  :: cell
      TYPE(t_nococonv),  INTENT(IN)  :: nococonv
      INTEGER,           INTENT(IN)  :: itype, bksym, lmax
      TYPE(t_abc),       INTENT(OUT) :: abc_star_spin(:)

      INTEGER :: spatial_iop, l, iOrd, band, iAtom_l, iAtom, mapped_atom, mapped_l
      INTEGER :: lm1, lm2, s, sp
      INTEGER :: mrot_one(3, 3, 1)
      LOGICAL :: l_time_reversal
      COMPLEX :: su_global(2, 2), su_local(2, 2)
      COMPLEX :: source(2, 2*lmax + 1), orbital_rot(2, 2*lmax + 1)
      COMPLEX, ALLOCATABLE :: d_wgn(:, :, :, :)

      IF (SIZE(abc_parent_spin) < 2 .OR. SIZE(abc_star_spin) < 2) THEN
         CALL juDFT_error("Spinor XAS star rotation requires two local-spin t_abc objects", calledby="m_xas_symmetry")
      END IF
      IF (.NOT. ALLOCATED(abc_parent_spin(1)%cof) .OR. .NOT. ALLOCATED(abc_parent_spin(2)%cof)) THEN
         CALL juDFT_error("abc_parent_spin%cof is not allocated in xas_rotate_abc_star_member_spinor", &
                          calledby="m_xas_symmetry")
      END IF

      CALL xas_star_operation(sym, bksym, spatial_iop, l_time_reversal)
      IF (l_time_reversal) THEN
         CALL juDFT_error("XAS spinor time-reversal star handling is not implemented yet", calledby="m_xas_symmetry")
      END IF

      abc_star_spin(1) = abc_parent_spin(1)
      abc_star_spin(2) = abc_parent_spin(2)
      abc_star_spin(1)%cof = CMPLX(0.0, 0.0)
      abc_star_spin(2)%cof = CMPLX(0.0, 0.0)

      IF (lmax >= 1) THEN
         ALLOCATE(d_wgn(-lmax:lmax, -lmax:lmax, 1:lmax, 1))
         mrot_one(:, :, 1) = sym%mrot(:, :, spatial_iop)
         CALL d_wigner(1, mrot_one, cell%bmat, lmax, d_wgn)
      END IF
      CALL xas_su2_from_sym(sym, cell, spatial_iop, su_global)

      DO iAtom_l = LBOUND(abc_parent_spin(1)%cof, 4), UBOUND(abc_parent_spin(1)%cof, 4)
         iAtom = atoms%firstAtom(itype) + iAtom_l - 1
         mapped_atom = iAtom
         IF (ALLOCATED(sym%mapped_atom)) THEN
            mapped_atom = sym%mapped_atom(spatial_iop, iAtom)
            IF (mapped_atom < atoms%firstAtom(itype) .OR. mapped_atom >= atoms%firstAtom(itype) + atoms%neq(itype)) THEN
               CALL juDFT_error("XAS spinor star rotation currently requires mapped absorber atoms to stay in the same atom type", &
                                calledby="m_xas_symmetry")
            END IF
            mapped_l = mapped_atom - atoms%firstAtom(itype) + 1
            IF (mapped_l /= iAtom_l) THEN
               CALL juDFT_error("XAS spinor star rotation currently assumes local absorbing atom maps to itself", &
                                calledby="m_xas_symmetry")
            END IF
         END IF
         CALL xas_local_spin_transform(atoms, nococonv, iAtom, mapped_atom, su_global, su_local)

         DO l = 0, lmax
            lm1 = l*l
            lm2 = l*(l + 2)
            IF (lm2 > UBOUND(abc_parent_spin(1)%cof, 2)) CYCLE
            DO iOrd = 1, MIN(abc_parent_spin(1)%n_r(l), SIZE(abc_parent_spin(1)%cof, 3))
               DO band = 1, SIZE(abc_parent_spin(1)%cof, 1)
                  source = CMPLX(0.0, 0.0)
                  source(1, 1:2*l + 1) = abc_parent_spin(1)%cof(band, lm1:lm2, iOrd, iAtom_l)
                  source(2, 1:2*l + 1) = abc_parent_spin(2)%cof(band, lm1:lm2, iOrd, iAtom_l)
                  IF (l == 0) THEN
                     orbital_rot(:, 1) = source(:, 1)
                  ELSE
                     DO sp = 1, 2
                        orbital_rot(sp, 1:2*l + 1) = MATMUL(source(sp, 1:2*l + 1), d_wgn(-l:l, -l:l, l, 1))
                     END DO
                  END IF
                  DO s = 1, 2
                     abc_star_spin(s)%cof(band, lm1:lm2, iOrd, iAtom_l) = &
                        su_local(s, 1)*orbital_rot(1, 1:2*l + 1) + &
                        su_local(s, 2)*orbital_rot(2, 1:2*l + 1)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE xas_rotate_abc_star_member_spinor

END MODULE m_xas_symmetry
