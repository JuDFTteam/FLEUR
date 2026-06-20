!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_symmetry
   USE m_constants, ONLY: tpi_const
   USE m_dwigner, ONLY: d_wigner
   USE m_grp_k, ONLY: euler
   USE m_inv3, ONLY: inv3
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
   PUBLIC :: xas_print_symmetry_rotation_diagnostics

   ! Temporary atom-mapping diagnostic for XAS star reconstruction. Enable
   ! manually when checking how sym%mapped_atom maps global/type-local indices.
   LOGICAL, PARAMETER :: xas_debug_atom_mapping = .FALSE.

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
      ! i.e. abc_local = U^\dagger abc_global. For expansion coefficients the
      ! star operation uses the inverse spinor rotation, S_R^\dagger. Therefore
      ! the parent-local to mapped-local transform is
      !   U_mapped^\dagger S_R^\dagger U_parent.
      u_parent = nococonv%umat(atoms%itype(parent_atom))
      u_mapped = nococonv%umat(atoms%itype(mapped_atom))
      su_local = MATMUL(CONJG(TRANSPOSE(u_mapped)), MATMUL(CONJG(TRANSPOSE(su_global)), u_parent))
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
      ! The hardwired XAS validation sums atoms explicitly. Symmetry operations
      ! may permute equivalent atoms of the selected type, so the rotated
      ! coefficients are written into the mapped type-local atom slot. Mappings
      ! to genuinely different atom types remain guarded.
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
      abc_star%cof = CMPLX(0.0, 0.0)
      IF (lmax >= 1) THEN
         ALLOCATE(d_wgn(-lmax:lmax, -lmax:lmax, 1:lmax, 1))
         mrot_one(:, :, 1) = sym%mrot(:, :, spatial_iop)
         CALL d_wigner(1, mrot_one, cell%bmat, lmax, d_wgn)
      END IF

      DO iAtom_l = LBOUND(abc_parent%cof, 4), UBOUND(abc_parent%cof, 4)
         iAtom = atoms%firstAtom(itype) + iAtom_l - 1
         mapped_atom = iAtom
         IF (ALLOCATED(sym%mapped_atom)) mapped_atom = sym%mapped_atom(spatial_iop, iAtom)
         IF (mapped_atom < 1 .OR. mapped_atom > atoms%nat) THEN
            mapped_atom = xas_find_mapped_absorber(atoms, sym, spatial_iop, itype, iAtom)
         END IF
         CALL xas_check_mapped_absorber(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, &
                                        "XAS star mode 3")

         IF (l_time_reversal) THEN
            abc_star%cof(:, 0, :, mapped_l) = CONJG(abc_parent%cof(:, 0, :, iAtom_l))
         ELSE
            abc_star%cof(:, 0, :, mapped_l) = abc_parent%cof(:, 0, :, iAtom_l)
         END IF

         DO l = 1, lmax
            lm1 = l*l
            lm2 = l*(l + 2)
            IF (lm2 > UBOUND(abc_parent%cof, 2)) CYCLE
            DO iOrd = 1, MIN(abc_parent%n_r(l), SIZE(abc_parent%cof, 3))
               DO band = 1, SIZE(abc_parent%cof, 1)
                  source(1:2*l + 1) = abc_parent%cof(band, lm1:lm2, iOrd, iAtom_l)
                  IF (l_time_reversal) source(1:2*l + 1) = CONJG(source(1:2*l + 1))
                  abc_star%cof(band, lm1:lm2, iOrd, mapped_l) = &
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
      ! S_local = U_mapped^\dagger S_global^\dagger U_parent follows
      ! abc%calc_abc, where local coefficients are U^\dagger times global
      ! spinors. The adjoint is required because these are expansion
      ! coefficients transformed to the star-member basis.
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
         IF (ALLOCATED(sym%mapped_atom)) mapped_atom = sym%mapped_atom(spatial_iop, iAtom)
         IF (mapped_atom < 1 .OR. mapped_atom > atoms%nat) THEN
            mapped_atom = xas_find_mapped_absorber(atoms, sym, spatial_iop, itype, iAtom)
         END IF
         CALL xas_check_mapped_absorber(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, &
                                        "XAS spinor star rotation")
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
                     abc_star_spin(s)%cof(band, lm1:lm2, iOrd, mapped_l) = &
                        su_local(s, 1)*orbital_rot(1, 1:2*l + 1) + &
                        su_local(s, 2)*orbital_rot(2, 1:2*l + 1)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE xas_rotate_abc_star_member_spinor

   SUBROUTINE xas_check_mapped_absorber(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, context)
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_sym),   INTENT(IN)  :: sym
      INTEGER,       INTENT(IN)  :: spatial_iop, itype, iAtom_l, mapped_atom
      INTEGER,       INTENT(OUT) :: mapped_l
      CHARACTER(*),  INTENT(IN)  :: context

      INTEGER :: iAtom

      iAtom = atoms%firstAtom(itype) + iAtom_l - 1
      mapped_l = -1
      IF (mapped_atom >= 1 .AND. mapped_atom <= atoms%nat) THEN
         mapped_l = mapped_atom - atoms%firstAtom(atoms%itype(mapped_atom)) + 1
      END IF

      IF (mapped_atom < 1 .OR. mapped_atom > atoms%nat) THEN
         IF (xas_debug_atom_mapping) THEN
            CALL xas_print_atom_mapping_diagnostic(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, context)
         END IF
         CALL juDFT_error(TRIM(context)//" found invalid mapped absorber atom index", calledby="m_xas_symmetry")
      END IF
      IF (atoms%itype(mapped_atom) /= itype) THEN
         IF (xas_debug_atom_mapping) THEN
            CALL xas_print_atom_mapping_diagnostic(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, context)
         END IF
         CALL juDFT_error(TRIM(context)//" currently requires mapped absorber atoms to stay in the same atom type", &
                          calledby="m_xas_symmetry")
      END IF
      mapped_l = mapped_atom - atoms%firstAtom(itype) + 1
      IF (mapped_l < 1 .OR. mapped_l > atoms%neq(itype)) THEN
         IF (xas_debug_atom_mapping) THEN
            CALL xas_print_atom_mapping_diagnostic(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, context)
         END IF
         CALL juDFT_error(TRIM(context)//" found inconsistent mapped type-local atom index", calledby="m_xas_symmetry")
      END IF

      IF (xas_debug_atom_mapping) THEN
         CALL xas_print_atom_mapping_diagnostic(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, context)
      END IF

      IF (iAtom < 1) CALL juDFT_error("Invalid absorber atom index in xas_check_mapped_absorber", &
                                      calledby="m_xas_symmetry")
   END SUBROUTINE xas_check_mapped_absorber

   INTEGER FUNCTION xas_find_mapped_absorber(atoms, sym, spatial_iop, itype, iAtom) RESULT(mapped_atom)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym),   INTENT(IN) :: sym
      INTEGER,       INTENT(IN) :: spatial_iop, itype, iAtom

      INTEGER :: candidate
      REAL    :: mapped_tau_raw(3), residual(3)
      REAL, PARAMETER :: eps = 1.0e-6

      mapped_atom = 0
      mapped_tau_raw = MATMUL(REAL(sym%mrot(:, :, spatial_iop)), atoms%taual(:, iAtom)) + sym%tau(:, spatial_iop)
      DO candidate = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
         residual = mapped_tau_raw - atoms%taual(:, candidate)
         residual = residual - REAL(NINT(residual))
         IF (MAXVAL(ABS(residual)) < eps) THEN
            mapped_atom = candidate
            EXIT
         END IF
      END DO
   END FUNCTION xas_find_mapped_absorber

   SUBROUTINE xas_print_atom_mapping_diagnostic(atoms, sym, spatial_iop, itype, iAtom_l, mapped_atom, mapped_l, context)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym),   INTENT(IN) :: sym
      INTEGER,       INTENT(IN) :: spatial_iop, itype, iAtom_l, mapped_atom, mapped_l
      CHARACTER(*),  INTENT(IN) :: context

      INTEGER :: iAtom, mapped_type, mapped_z
      REAL    :: mapped_tau_raw(3), lattice_shift(3)

      iAtom = atoms%firstAtom(itype) + iAtom_l - 1
      mapped_type = -1
      mapped_z = -1
      mapped_tau_raw = MATMUL(REAL(sym%mrot(:, :, spatial_iop)), atoms%taual(:, iAtom)) + sym%tau(:, spatial_iop)
      lattice_shift = 0.0
      IF (mapped_atom >= 1 .AND. mapped_atom <= atoms%nat) THEN
         mapped_type = atoms%itype(mapped_atom)
         mapped_z = atoms%nz(mapped_type)
         lattice_shift = mapped_tau_raw - atoms%taual(:, mapped_atom)
         lattice_shift = lattice_shift - REAL(NINT(lattice_shift))
      END IF

      WRITE(*, '(a)') "XAS DEBUG atom mapping diagnostic: "//TRIM(context)
      WRITE(*, '(a,i0,a,i0,a,i0,a,i0)') "  original_global=", iAtom, " original_type=", itype, &
                                         " original_local=", iAtom_l, " original_Z=", atoms%nz(itype)
      WRITE(*, '(a,i0,a,i0)') "  spatial_iop=", spatial_iop, " mapped_global=", mapped_atom
      WRITE(*, '(a,i0,a,i0,a,i0)') "  mapped_type=", mapped_type, " mapped_local=", mapped_l, &
                                   " mapped_Z=", mapped_z
      WRITE(*, '(a,3es16.8)') "  tau_original=", atoms%taual(:, iAtom)
      WRITE(*, '(a,3es16.8)') "  tau_after_sym=", mapped_tau_raw
      IF (mapped_atom >= 1 .AND. mapped_atom <= atoms%nat) THEN
         WRITE(*, '(a,3es16.8)') "  tau_mapped=", atoms%taual(:, mapped_atom)
         WRITE(*, '(a,3es16.8)') "  residual_after_lattice_translation=", lattice_shift
      END IF
   END SUBROUTINE xas_print_atom_mapping_diagnostic

   SUBROUTINE xas_print_symmetry_rotation_diagnostics(sym, cell, unit)
      TYPE(t_sym),  INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER,      INTENT(IN) :: unit

      INTEGER :: iop, irow, det_int
      REAL    :: r_cart(3, 3), r_check(3, 3), r_wigner(3, 3), bmat_inv(3, 3)
      REAL    :: alpha, beta, gamma, det_real, dt, ortho_norm
      LOGICAL :: l_axis_exchange

      CALL inv3(cell%bmat, bmat_inv, dt)
      CALL xas_debug_write_line(unit, "XAS DEBUG SYM ROT: operation raw_mrot R_cart det orthonorm proper/improper Euler")
      DO iop = 1, sym%nop
         det_int = xas_det3_int(sym%mrot(:, :, iop))
         CALL xas_cart_rotation_from_sym(sym, cell, iop, r_cart)
         det_real = xas_det3_real(r_cart)
         r_check = MATMUL(TRANSPOSE(r_cart), r_cart)
         r_check(1, 1) = r_check(1, 1) - 1.0
         r_check(2, 2) = r_check(2, 2) - 1.0
         r_check(3, 3) = r_check(3, 3) - 1.0
         ortho_norm = MAXVAL(ABS(r_check))
         CALL euler(iop, sym, cell, alpha, beta, gamma)

         ! This is the proper Cartesian matrix whose inverse is used inside
         ! m_dwigner when raw sym%mrot and cell%bmat are passed to d_wigner.
         r_wigner = REAL(det_int)*MATMUL(bmat_inv, MATMUL(REAL(sym%mrot(:, :, iop)), cell%bmat))

         l_axis_exchange = xas_is_axis_permutation(r_cart)
         WRITE(*, '(a,i0,a,i0,a,es12.4,a,l1,a,3es14.6)') &
            "XAS DEBUG SYM ROT iop=", iop, " det_int=", det_int, " det_cart=", det_real, &
            " axis_perm=", l_axis_exchange, " euler=", alpha, beta, gamma
         IF (unit > 0) WRITE(unit, '(a,i0,a,i0,a,es12.4,a,l1,a,3es14.6)') &
            "XAS DEBUG SYM ROT iop=", iop, " det_int=", det_int, " det_cart=", det_real, &
            " axis_perm=", l_axis_exchange, " euler=", alpha, beta, gamma
         WRITE(*, '(a,es12.4)') "XAS DEBUG SYM ROT orthonorm=", ortho_norm
         IF (unit > 0) WRITE(unit, '(a,es12.4)') "XAS DEBUG SYM ROT orthonorm=", ortho_norm
         DO irow = 1, 3
            WRITE(*, '(a,i0,a,3i5,a,3es14.6,a,3es14.6)') "XAS DEBUG SYM ROT row ", irow, &
               " raw=", sym%mrot(irow, :, iop), " R_cart=", r_cart(irow, :), " R_wigner_proper=", r_wigner(irow, :)
            IF (unit > 0) WRITE(unit, '(a,i0,a,3i5,a,3es14.6,a,3es14.6)') "XAS DEBUG SYM ROT row ", irow, &
               " raw=", sym%mrot(irow, :, iop), " R_cart=", r_cart(irow, :), " R_wigner_proper=", r_wigner(irow, :)
         END DO
      END DO
   END SUBROUTINE xas_print_symmetry_rotation_diagnostics

   INTEGER FUNCTION xas_det3_int(mat) RESULT(det)
      INTEGER, INTENT(IN) :: mat(3, 3)

      det = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(3, 2)*mat(2, 3)) + &
            mat(1, 2)*(mat(2, 3)*mat(3, 1) - mat(2, 1)*mat(3, 3)) + &
            mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))
   END FUNCTION xas_det3_int

   REAL FUNCTION xas_det3_real(mat) RESULT(det)
      REAL, INTENT(IN) :: mat(3, 3)

      det = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(3, 2)*mat(2, 3)) + &
            mat(1, 2)*(mat(2, 3)*mat(3, 1) - mat(2, 1)*mat(3, 3)) + &
            mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))
   END FUNCTION xas_det3_real

   LOGICAL FUNCTION xas_is_axis_permutation(mat) RESULT(l_perm)
      REAL, INTENT(IN) :: mat(3, 3)

      INTEGER :: irow, icol
      REAL, PARAMETER :: eps = 1.0e-6

      l_perm = .TRUE.
      DO irow = 1, 3
         IF (COUNT(ABS(ABS(mat(irow, :)) - 1.0) < eps) /= 1) l_perm = .FALSE.
         IF (COUNT(ABS(mat(irow, :)) > eps) /= 1) l_perm = .FALSE.
      END DO
      DO icol = 1, 3
         IF (COUNT(ABS(ABS(mat(:, icol)) - 1.0) < eps) /= 1) l_perm = .FALSE.
         IF (COUNT(ABS(mat(:, icol)) > eps) /= 1) l_perm = .FALSE.
      END DO
   END FUNCTION xas_is_axis_permutation

   SUBROUTINE xas_debug_write_line(unit, text)
      INTEGER,          INTENT(IN) :: unit
      CHARACTER(LEN=*), INTENT(IN) :: text

      WRITE(*, '(a)') TRIM(text)
      IF (unit > 0) WRITE(unit, '(a)') TRIM(text)
   END SUBROUTINE xas_debug_write_line

END MODULE m_xas_symmetry
