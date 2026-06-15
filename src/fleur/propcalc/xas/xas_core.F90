!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xas_core
   USE m_constants, ONLY: c_light
   USE m_differ, ONLY: differ
   USE m_intgr, ONLY: intgr3
   USE m_juDFT, ONLY: juDFT_error
   USE m_types_atoms, ONLY: t_atoms
   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: t_xas_core_state
      CHARACTER(LEN=8) :: edge_name = ""
      INTEGER          :: n = 0
      INTEGER          :: lc = -1
      INTEGER          :: kappa = 0
      INTEGER          :: twice_j = 0
      REAL             :: energy = 0.0
      REAL             :: occupation = 0.0
      INTEGER          :: nr = 0
      REAL             :: norm = 0.0
      INTEGER, ALLOCATABLE :: twice_mj(:)
      REAL,    ALLOCATABLE :: p_core(:)
      REAL,    ALLOCATABLE :: q_core(:)
   END TYPE t_xas_core_state

   PUBLIC :: xas_extract_core_states
   PUBLIC :: xas_print_core_states

CONTAINS

   SUBROUTINE xas_extract_core_states(atoms, itype, edge_name, v_mt, core_states)
      TYPE(t_atoms),                   INTENT(IN)  :: atoms
      INTEGER,                         INTENT(IN)  :: itype
      CHARACTER(LEN=*),                INTENT(IN)  :: edge_name
      REAL,                            INTENT(IN)  :: v_mt(:)
      TYPE(t_xas_core_state), ALLOCATABLE, INTENT(OUT) :: core_states(:)

      INTEGER :: edge_lc, edge_n, edge_twice_j, edge_kappa
      INTEGER :: nst, i_state, n_selected, i_selected
      INTEGER :: nr, msh, ierr, i_mj
      INTEGER :: nprnc(maxval(atoms%econf%num_states))
      INTEGER :: kappa(maxval(atoms%econf%num_states))
      REAL    :: occ(maxval(atoms%econf%num_states), 1)
      REAL    :: vrd(atoms%msh), a(atoms%msh), b(atoms%msh), norm_integrand(atoms%jri(itype))
      REAL    :: d, rn, energy, fj

      CALL xas_edge_quantum_numbers(edge_name, edge_n, edge_lc, edge_twice_j, edge_kappa)

      nr = atoms%jri(itype)
      msh = atoms%msh
      IF (SIZE(v_mt) < nr) THEN
         CALL juDFT_error("v_mt is shorter than the MT radial mesh in xas_extract_core_states", calledby="m_xas_core")
      END IF
      IF (msh < nr) THEN
         CALL juDFT_error("atoms%msh is smaller than atoms%jri in xas_extract_core_states", calledby="m_xas_core")
      END IF

      CALL atoms%econf(itype)%get_core(nst, nprnc, kappa, occ)
      n_selected = 0
      DO i_state = 1, nst
         IF (nprnc(i_state) /= edge_n) CYCLE
         IF (kappa(i_state) /= edge_kappa) CYCLE
         IF (xas_l_from_kappa(kappa(i_state)) /= edge_lc) CYCLE
         IF (xas_twice_j_from_kappa(kappa(i_state)) /= edge_twice_j) CYCLE
         n_selected = n_selected + 1
      END DO

      ALLOCATE(core_states(n_selected))
      IF (n_selected == 0) RETURN

      vrd(1:nr) = v_mt(1:nr)
      IF (msh > nr) THEN
         vrd(nr + 1:msh) = 0.0
         DO i_state = nr + 1, msh
            vrd(i_state) = vrd(nr) + vrd(nr) / REAL(nr - msh) * REAL(i_state - nr)
         END DO
      END IF

      d = EXP(atoms%dx(itype))
      rn = atoms%rmsh(1, itype)*(d**(msh - 1))
      i_selected = 0
      DO i_state = 1, nst
         IF (nprnc(i_state) /= edge_n) CYCLE
         IF (kappa(i_state) /= edge_kappa) CYCLE
         IF (xas_l_from_kappa(kappa(i_state)) /= edge_lc) CYCLE
         IF (xas_twice_j_from_kappa(kappa(i_state)) /= edge_twice_j) CYCLE

         i_selected = i_selected + 1
         fj = 0.5*REAL(edge_twice_j)
         energy = -2.0*(atoms%zatom(itype)/REAL(edge_n + edge_lc))**2
         CALL differ(REAL(edge_n), REAL(edge_lc), fj, c_light(1.0), atoms%zatom(itype), atoms%dx(itype), &
            atoms%rmsh(1, itype), rn, d, msh, vrd, energy, a, b, ierr)
         IF (ierr /= 0) THEN
            CALL juDFT_error("differ failed in xas_extract_core_states", calledby="m_xas_core")
         END IF

         core_states(i_selected)%edge_name = ADJUSTL(edge_name)
         core_states(i_selected)%n = edge_n
         core_states(i_selected)%lc = edge_lc
         core_states(i_selected)%kappa = kappa(i_state)
         core_states(i_selected)%twice_j = edge_twice_j
         core_states(i_selected)%energy = energy
         core_states(i_selected)%occupation = occ(i_state, 1)
         core_states(i_selected)%nr = nr

         ALLOCATE(core_states(i_selected)%twice_mj(edge_twice_j + 1))
         DO i_mj = 1, edge_twice_j + 1
            core_states(i_selected)%twice_mj(i_mj) = -edge_twice_j + 2*(i_mj - 1)
         END DO

         ALLOCATE(core_states(i_selected)%p_core(nr), core_states(i_selected)%q_core(nr))
         core_states(i_selected)%p_core(:) = a(1:nr)
         core_states(i_selected)%q_core(:) = b(1:nr)
         norm_integrand(:) = a(1:nr)**2 + b(1:nr)**2
         CALL intgr3(norm_integrand, atoms%rmsh(1:nr, itype), atoms%dx(itype), nr, core_states(i_selected)%norm)
      END DO
   END SUBROUTINE xas_extract_core_states

   SUBROUTINE xas_print_core_states(core_states, out_unit)
      TYPE(t_xas_core_state), INTENT(IN) :: core_states(:)
      INTEGER, OPTIONAL,      INTENT(IN) :: out_unit

      INTEGER :: i_state, unit

      unit = 6
      IF (PRESENT(out_unit)) unit = out_unit

      WRITE (unit, '(a)') "XAS core states:"
      DO i_state = 1, SIZE(core_states)
         WRITE (unit, '(a,a,a,i3,a,i3,a,i4,a,i4,a,es18.10,a,f8.4,a,i6,a,es14.6)') &
            "  edge=", TRIM(core_states(i_state)%edge_name), &
            " n=", core_states(i_state)%n, &
            " lc=", core_states(i_state)%lc, &
            " kappa=", core_states(i_state)%kappa, &
            " twice_j=", core_states(i_state)%twice_j, &
            " energy=", core_states(i_state)%energy, &
            " occupation=", core_states(i_state)%occupation, &
            " nr=", core_states(i_state)%nr, &
            " norm=", core_states(i_state)%norm
      END DO
   END SUBROUTINE xas_print_core_states

   SUBROUTINE xas_edge_quantum_numbers(edge_name, n, lc, twice_j, kappa)
      CHARACTER(LEN=*), INTENT(IN)  :: edge_name
      INTEGER,          INTENT(OUT) :: n, lc, twice_j, kappa

      CHARACTER(LEN=:), ALLOCATABLE :: edge

      edge = ADJUSTL(edge_name)
      SELECT CASE (TRIM(edge))
      CASE ("L3", "l3", "2p3/2", "2P3/2")
         n = 2
         lc = 1
         twice_j = 3
         kappa = -2
      CASE ("L2", "l2", "2p1/2", "2P1/2")
         n = 2
         lc = 1
         twice_j = 1
         kappa = 1
      CASE DEFAULT
         CALL juDFT_error("Unsupported XAS edge in xas_edge_quantum_numbers", calledby="m_xas_core")
      END SELECT
   END SUBROUTINE xas_edge_quantum_numbers

   INTEGER FUNCTION xas_l_from_kappa(kappa) RESULT(l)
      INTEGER, INTENT(IN) :: kappa

      IF (kappa > 0) THEN
         l = kappa
      ELSE
         l = -kappa - 1
      END IF
   END FUNCTION xas_l_from_kappa

   INTEGER FUNCTION xas_twice_j_from_kappa(kappa) RESULT(twice_j)
      INTEGER, INTENT(IN) :: kappa

      twice_j = 2*ABS(kappa) - 1
   END FUNCTION xas_twice_j_from_kappa

END MODULE m_xas_core
