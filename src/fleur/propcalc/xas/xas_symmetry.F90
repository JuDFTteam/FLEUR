!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_symmetry
   USE m_constants, ONLY: tpi_const
   USE m_dwigner, ONLY: d_wigner
   USE m_juDFT, ONLY: juDFT_error
   USE m_types_abc, ONLY: t_abc
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_cell, ONLY: t_cell
   USE m_types_kpts, ONLY: t_kpts
   USE m_types_sym, ONLY: t_sym
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: xas_count_star_members
   PUBLIC :: xas_star_member_weight
   PUBLIC :: xas_star_operation
   PUBLIC :: xas_cart_rotation_from_sym
   PUBLIC :: xas_rotate_lab_polarization_for_parent
   PUBLIC :: xas_rotate_abc_star_member

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

END MODULE m_xas_symmetry
