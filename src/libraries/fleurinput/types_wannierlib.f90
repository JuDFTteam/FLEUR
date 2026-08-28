!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_types_wannierlib
  USE m_types_melem_optable, ONLY: WANNIERLIB_INTERP, WANNIERLIB_OPR, melem_exposed_find
  USE m_types_kpts
  USE m_juDFT
  USE m_types_atoms
  USE m_types_noco
  USE m_types_fleurinput_base
  USE m_constants,ONLY: ounit
  IMPLICIT NONE
  PRIVATE

  TYPE, EXTENDS(t_fleurinput_base) :: t_wannierlib_wannierize
    LOGICAL :: l_wannierize = .FALSE.
    ! Convenience flags DERIVED from the operator list ops(:) below (set in read_xml).
    ! They gate the coarse-mesh provider calls; the actual dispatch loops over ops(:).
    LOGICAL :: l_interpolation = .FALSE.   ! an <operator name="hamiltonian"> is requested
    LOGICAL :: l_spin = .FALSE.            ! an <operator name="spin"> is requested
    LOGICAL :: l_orbmom = .FALSE.          ! an <operator name="orbital"> is requested
    LOGICAL :: l_socop = .FALSE.           ! an <operator name="soc"> is requested
    LOGICAL :: l_operators_r = .FALSE.     ! an <operators_r> block (real-space O(R) export) is present
    !> Opt-in: put the Wannier functions themselves on a real-space grid and write them
    !> as XSF. Off by default because it costs a second pass over the k-points, reading
    !> the states back once the gauge is known; nothing else in the run needs it.
    !> XML @plotWF.
    LOGICAL :: l_plot_wf = .FALSE.

    ! --- opt-in output domains (<domain> children of <interpolation>) ---
    ! One <domain> per set of k-points the operators are evaluated on, repeatable and
    ! without kinds: @listName (required) names a kPointList of kpts.xml, @suffix tails the
    ! output names, @npts subdivides each segment of the list. Operators declared with no
    ! domain are an error rather than a silent no-op.
    !> One entry per <domain>. The count and the suffixes ARE broadcast -- the first sets
    !> the domain loop bound on every rank, the second is read inside that loop, where the
    !> output name is built. The k-sets are not: only rank 0 reaches the writers. Reading
    !> dom_kset off rank 0 is undefined behaviour, and a two-rank run can pass with it.
    INTEGER :: n_domains = 0
    TYPE(t_kpts), ALLOCATABLE :: dom_kset(:)
    CHARACTER(LEN=64), ALLOCATABLE :: dom_suffix(:)

    !> Minimum distance replica selection (W90 use_ws_distance); XML @useWsDistance on
    !> <interpolation>. Off by default; Wannier90's own default is on, and the difference is
    !> deliberate -- with it off this reproduces every band this code has written so far.
    !>
    !> The plain Fourier sum phases the hopping <0m|H|Rn> as e^{ikR}, as if it joined cell
    !> ORIGINS. It does not: m sits at tau_m and n at tau_n, so the real separation is
    !> R + tau_n - tau_m. Assigning by R alone is discontinuous at the Wigner-Seitz boundary,
    !> H(R) then decays slowly, and a slowly decaying Fourier series rings. Averaging over
    !> the tied minimum-distance replicas restores the real distance.
    !>
    !> It earns its keep where the centres sit far from the origin -- bond-centred, or pushed
    !> into the vacuum of a film -- and where the cell is anisotropic.
    LOGICAL :: l_ws_distance = .FALSE.

    ! --- operator table: one entry per <operator name=".."> child of <interpolation> ---
    INTEGER :: n_ops = 0
    CHARACTER(LEN=20), ALLOCATABLE :: op_name(:)    ! operator identifier (hamiltonian/spin/orbital/soc/...)
    CHARACTER(LEN=32), ALLOCATABLE :: op_comp(:)    ! requested components (Phase 2); '' = all of the operator rank
    INTEGER, ALLOCATABLE :: op_total(:)             ! 1 = write summed-over-atoms projection (default)

    ! --- operators_r table: one entry per <operators_r>/<operator name=".."> child ---
    INTEGER :: n_op_r = 0
    CHARACTER(LEN=20), ALLOCATABLE :: op_r_name(:)  ! hamiltonian/position/spin/orbital(/spin_orbit)

    INTEGER :: num_wann = 0
    INTEGER :: num_bands = 0
    INTEGER :: min_band = -1
    INTEGER :: max_band = -1

    REAL :: dis_win_min = 0.0
    REAL :: dis_win_max = 0.0
    REAL :: dis_froz_min = 0.0
    REAL :: dis_froz_max = 0.0
    INTEGER :: dis_num_iter = 0
    INTEGER :: num_iter = 0      ! MLWF/wannierise iterations (W90 num_iter); XML @wannNumIter
    !> Projectability disentanglement (W90 dis_froz_proj / dis_proj_min / dis_proj_max).
    !> Off by default, which is also Wannier90's default -- and W90 says why in its own
    !> source: "upon reading AMN we do not know where it comes from". We do: ours are the
    !> atomic projections this code builds, so the criterion is meaningful here.
    !>
    !> It replaces the ENERGY frozen window by one on p_i(k) = sum_j |A_ij(k)|^2, the share
    !> of a Bloch state that lies in the span of the trial orbitals. Below proj_min a state
    !> is discarded outright -- even from the middle of the outer window, which an energy cut
    !> cannot do -- and at or above proj_max it is frozen.
    !>
    !> It needs a full-rank amn: with the spinor guard broken the matrix came out rank
    !> num_wann/2 and p ran up to 1.98, which W90 rejects outright. Fixed in 38b46f7e4.
    LOGICAL :: dis_froz_proj = .FALSE.
    REAL :: dis_proj_min = 0.01
    REAL :: dis_proj_max = 0.95
    REAL :: dis_mix_ratio = 0.0
    REAL :: dis_conv_tol = 0.0
    REAL :: conv_tol = 0.0       ! MLWF/wannierise convergence (W90 conv_tol); XML @wannConvTol
    !> Preconditioned gradient in the spread minimisation (W90 precond); XML @precond.
    !> Off by default, which is also Wannier90's default.
    !>
    !> It changes the descent inside a valley, not the functional, so it cannot move the
    !> minimum -- and measurement confirms it does not: on bcc Fe with the quantisation
    !> axis along x, y and z it leaves Omega_I identical to every digit and Omega_total
    !> within 0.02%, and the axis whose localisation was broken came out slightly worse.
    !> It is exposed because the option exists and someone will want it, not because it
    !> fixed anything here.
    LOGICAL :: precond = .FALSE.

    INTEGER, ALLOCATABLE :: proj_ntype(:)
    INTEGER, ALLOCATABLE :: proj_atom(:)
    CHARACTER(LEN=20), ALLOCATABLE :: proj_species(:)
    INTEGER, ALLOCATABLE :: proj_l(:)
    INTEGER, ALLOCATABLE :: proj_m(:)
    INTEGER, ALLOCATABLE :: proj_spin(:)
    INTEGER, ALLOCATABLE :: proj_rwf(:)
    REAL, ALLOCATABLE :: proj_alpha(:)
    REAL, ALLOCATABLE :: proj_beta(:)
    REAL, ALLOCATABLE :: proj_gamma(:)
    REAL, ALLOCATABLE :: proj_zona(:)
    REAL, ALLOCATABLE :: proj_regio(:)
    REAL, ALLOCATABLE :: proj_j(:)
    REAL, ALLOCATABLE :: proj_mj(:)
    REAL, ALLOCATABLE :: proj_weight(:)
    REAL, ALLOCATABLE :: proj_shift(:, :)
  CONTAINS
    PROCEDURE :: init => init_wannierlib
    PROCEDURE :: read_xml => read_xml_wannierlib
    PROCEDURE :: mpi_bc => mpi_bc_wannierlib
  END TYPE t_wannierlib_wannierize

  PUBLIC :: t_wannierlib_wannierize

CONTAINS

  SUBROUTINE init_wannierlib(this, atoms, noco)
    CLASS(t_wannierlib_wannierize), INTENT(INOUT) :: this
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_noco), INTENT(IN) :: noco

    INTEGER :: i, j, nold, nnew, mcount, mval, mrepeat, spin_mult, ispin, itype, type_mult, nn, na
    INTEGER, ALLOCATABLE :: new_proj_ntype(:), new_proj_atom(:), new_proj_l(:), new_proj_m(:), new_proj_spin(:), new_proj_rwf(:)
    CHARACTER(LEN=20), ALLOCATABLE :: new_proj_species(:)
    REAL, ALLOCATABLE :: new_proj_alpha(:), new_proj_beta(:), new_proj_gamma(:), new_proj_zona(:), new_proj_regio(:)
    REAL, ALLOCATABLE :: new_proj_j(:), new_proj_mj(:), new_proj_weight(:), new_proj_shift(:, :)
    LOGICAL :: expand_spin

    IF (.NOT.ALLOCATED(this%proj_l)) RETURN !No Wannier projections configured, nothing to do
    IF (.NOT.ALLOCATED(atoms%speciesName)) CALL juDFT_error('wannierlib: atoms%speciesName not allocated in init', calledby='init_wannierlib')
    
    nold = SIZE(this%proj_l)

    IF (.NOT.ALLOCATED(this%proj_ntype)) THEN
      ALLOCATE(this%proj_ntype(nold))
      this%proj_ntype = 0
    END IF
    IF (.NOT.ALLOCATED(this%proj_atom)) THEN
      ALLOCATE(this%proj_atom(nold))
      this%proj_atom = 1
    END IF

    expand_spin = noco%l_noco .OR. noco%l_soc

    nnew = 0
    DO i = 1, nold
      mcount = projection_m_count(this%proj_l(i))

      IF (this%proj_m(i) == 0) THEN
        IF (mcount <= 0) THEN
          CALL juDFT_error('wannierlib: unsupported projection l for m expansion', calledby='init_wannierlib')
        END IF
        mrepeat = mcount
      ELSE
        IF ((mcount > 0) .AND. ((this%proj_m(i) < 1) .OR. (this%proj_m(i) > mcount))) THEN
          CALL juDFT_error('wannierlib: projection m out of range for l', calledby='init_wannierlib')
        END IF
        mrepeat = 1
      END IF

      spin_mult = 1
      IF (expand_spin .AND. (this%proj_spin(i) == 0)) spin_mult = 2

      type_mult = 0
      DO itype = 1, atoms%ntype
        IF (TRIM(atoms%speciesName(itype)) == TRIM(this%proj_species(i))) type_mult = type_mult + atoms%neq(itype)
      END DO
      IF (type_mult <= 0) THEN
        CALL juDFT_error('wannierlib: no matching atom type for projection species name', calledby='init_wannierlib')
      END IF

      nnew = nnew + mrepeat * spin_mult * type_mult
    END DO

    ALLOCATE(new_proj_ntype(nnew), new_proj_atom(nnew), new_proj_species(nnew), new_proj_l(nnew), new_proj_m(nnew), new_proj_spin(nnew), new_proj_rwf(nnew))
    ALLOCATE(new_proj_alpha(nnew), new_proj_beta(nnew), new_proj_gamma(nnew), new_proj_zona(nnew), new_proj_regio(nnew))
    ALLOCATE(new_proj_j(nnew), new_proj_mj(nnew), new_proj_weight(nnew), new_proj_shift(3, nnew))

    j = 0
    DO i = 1, nold
      mcount = projection_m_count(this%proj_l(i))

      IF ((this%proj_m(i) == 0) .AND. (mcount > 0)) THEN
        mrepeat = mcount
      ELSE
        mrepeat = 1
      END IF

      spin_mult = 1
      IF (expand_spin .AND. (this%proj_spin(i) == 0)) spin_mult = 2

      DO itype = 1, atoms%ntype
        IF (TRIM(atoms%speciesName(itype)) /= TRIM(this%proj_species(i))) CYCLE
        DO nn = 1, atoms%neq(itype)
          na = atoms%firstAtom(itype) + nn - 1
          DO mval = 1, mrepeat
          DO ispin = 1, spin_mult
            j = j + 1
            new_proj_ntype(j) = itype
            new_proj_atom(j) = nn
            new_proj_species(j) = this%proj_species(i)
            new_proj_l(j) = this%proj_l(i)
            ! BUGFIX: when the user requested m=0 ("expand all m"), the auto-generated
            ! harmonic index must be the 1-based mval -- even when only one m exists
            ! (l=0/s: mrepeat=1). The old test (mrepeat>1) left s with the literal
            ! proj_m=0, which is out of range for tlm(:,:,1:7) -> uninitialised read.
            IF (this%proj_m(i) == 0) THEN
              new_proj_m(j) = mval
            ELSE
              new_proj_m(j) = this%proj_m(i)
            END IF
            IF (spin_mult == 2) THEN
              IF (ispin == 1) THEN
                new_proj_spin(j) = 1
              ELSE
                new_proj_spin(j) = -1
              END IF
            ELSE
              new_proj_spin(j) = this%proj_spin(i)
            END IF
            new_proj_rwf(j) = this%proj_rwf(i)
            new_proj_alpha(j) = this%proj_alpha(i)
            new_proj_beta(j) = this%proj_beta(i)
            new_proj_gamma(j) = this%proj_gamma(i)
            new_proj_zona(j) = this%proj_zona(i)
            new_proj_regio(j) = this%proj_regio(i)
            new_proj_j(j) = this%proj_j(i)
            new_proj_mj(j) = this%proj_mj(i)
            new_proj_weight(j) = this%proj_weight(i)
            new_proj_shift(:, j) = this%proj_shift(:, i) 
          END DO
          END DO
        END DO
      END DO
    END DO

    CALL move_alloc(new_proj_ntype, this%proj_ntype)
    CALL move_alloc(new_proj_atom, this%proj_atom)
    CALL move_alloc(new_proj_species, this%proj_species)
    CALL move_alloc(new_proj_l, this%proj_l)
    CALL move_alloc(new_proj_m, this%proj_m)
    CALL move_alloc(new_proj_spin, this%proj_spin)
    CALL move_alloc(new_proj_rwf, this%proj_rwf)
    CALL move_alloc(new_proj_alpha, this%proj_alpha)
    CALL move_alloc(new_proj_beta, this%proj_beta)
    CALL move_alloc(new_proj_gamma, this%proj_gamma)
    CALL move_alloc(new_proj_zona, this%proj_zona)
    CALL move_alloc(new_proj_regio, this%proj_regio)
    CALL move_alloc(new_proj_j, this%proj_j)
    CALL move_alloc(new_proj_mj, this%proj_mj)
    CALL move_alloc(new_proj_weight, this%proj_weight)
    CALL move_alloc(new_proj_shift, this%proj_shift)

    this%num_wann = nnew

    IF (this%min_band == 0) THEN
      IF (noco%l_noco .OR. noco%l_soc) THEN
        this%min_band = 2 * atoms%nlotot + 1
      ELSE
        this%min_band = atoms%nlotot + 1
      END IF
    END IF

    IF ((this%max_band == 0) .AND. (this%num_bands > 0)) THEN
      this%max_band = this%min_band + this%num_bands - 1
    END IF

    IF (this%num_bands == 0) THEN
      IF ((this%min_band > 0) .AND. (this%max_band > 0)) THEN
        this%num_bands = this%max_band - this%min_band + 1
      END IF
    END IF

    IF ((this%min_band > 0) .AND. (this%max_band > 0) .AND. (this%num_bands > 0)) THEN
      IF (this%max_band /= this%min_band + this%num_bands - 1) THEN
        CALL juDFT_error('wannierlib: inconsistent numBands and [minBand,maxBand]', calledby='init_wannierlib')
      END IF
    END IF

    WRITE(oUnit, '(A)') 'wannierlib: expanded projection list'
    WRITE(oUnit, '(A)') ' idx  species               type atom   l   m spin rwf      alpha       beta      gamma       zona      regio          j         mj    weight      shift_x      shift_y      shift_z'
    WRITE(oUnit, '(A)') '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    DO i = 1, this%num_wann
      WRITE(oUnit, '(I4,1X,A20,1X,I4,1X,I4,1X,I3,1X,I3,1X,I4,1X,I3,1X,&
     &             F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,&
     &             F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5)') &
           i, TRIM(this%proj_species(i)), this%proj_ntype(i), this%proj_atom(i), this%proj_l(i), this%proj_m(i), this%proj_spin(i), this%proj_rwf(i), &
           this%proj_alpha(i), this%proj_beta(i), this%proj_gamma(i), this%proj_zona(i), this%proj_regio(i), this%proj_j(i), this%proj_mj(i), this%proj_weight(i), &
           this%proj_shift(1, i), this%proj_shift(2, i), this%proj_shift(3, i)
    END DO
    WRITE(oUnit, '(A)') '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  END SUBROUTINE init_wannierlib

  INTEGER FUNCTION projection_m_count(l)
    INTEGER, INTENT(IN) :: l

    SELECT CASE (l)
    CASE (0)
      projection_m_count = 1
    CASE (1)
      projection_m_count = 3
    CASE (2)
      projection_m_count = 5
    CASE (3)
      projection_m_count = 7
    CASE (-1)
      projection_m_count = 2
    CASE (-2)
      projection_m_count = 3
    CASE (-3)
      projection_m_count = 4
    CASE (-4)
      projection_m_count = 5
    CASE (-5)
      projection_m_count = 6
    CASE DEFAULT
      projection_m_count = 0
    END SELECT
  END FUNCTION projection_m_count

  SUBROUTINE mpi_bc_wannierlib(this, mpi_comm, irank)
    USE m_mpi_bc_tool
    CLASS(t_wannierlib_wannierize), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: mpi_comm
    INTEGER, INTENT(IN), OPTIONAL :: irank
    INTEGER :: rank

    rank = 0
    IF (PRESENT(irank)) rank = irank

    CALL mpi_bc(this%l_wannierize, rank, mpi_comm)
    CALL mpi_bc(this%l_interpolation, rank, mpi_comm)
    CALL mpi_bc(this%l_ws_distance, rank, mpi_comm)
    CALL mpi_bc(this%l_spin, rank, mpi_comm)
    CALL mpi_bc(this%l_orbmom, rank, mpi_comm)
    CALL mpi_bc(this%l_socop, rank, mpi_comm)
    ! Only the domain COUNT is broadcast: it sets the loop bound on every rank. The k-sets
    ! and suffixes stay on rank 0, which is the only one that writes the k-set file and
    ! renames the outputs -- so on any other rank dom_kset is unallocated, and reading or
    ! validating it there is a bug, not a missing broadcast.
    CALL mpi_bc(this%n_domains, rank, mpi_comm)
    !> The suffixes ARE broadcast, unlike the k-sets: they are a handful of short strings,
    !> and the output name is built inside the domain loop, which every rank turns. Leaving
    !> them on rank 0 segfaulted rank 7 -- and not rank 1, so a two-rank suite passed.
    CALL mpi_bc(this%dom_suffix, rank, mpi_comm)
    CALL mpi_bc(this%n_ops, rank, mpi_comm)
    CALL mpi_bc(this%op_name, rank, mpi_comm)
    CALL mpi_bc(this%op_comp, rank, mpi_comm)
    CALL mpi_bc(this%op_total, rank, mpi_comm)
    CALL mpi_bc(this%l_operators_r, rank, mpi_comm)
    CALL mpi_bc(this%l_plot_wf, rank, mpi_comm)
    CALL mpi_bc(this%n_op_r, rank, mpi_comm)
    CALL mpi_bc(this%op_r_name, rank, mpi_comm)
    CALL mpi_bc(this%num_wann, rank, mpi_comm)
    CALL mpi_bc(this%num_bands, rank, mpi_comm)
    CALL mpi_bc(this%min_band, rank, mpi_comm)
    CALL mpi_bc(this%max_band, rank, mpi_comm)
    CALL mpi_bc(this%dis_win_min, rank, mpi_comm)
    CALL mpi_bc(this%dis_win_max, rank, mpi_comm)
    CALL mpi_bc(this%dis_froz_min, rank, mpi_comm)
    CALL mpi_bc(this%dis_froz_max, rank, mpi_comm)
    CALL mpi_bc(this%dis_num_iter, rank, mpi_comm)
    CALL mpi_bc(this%num_iter, rank, mpi_comm)
    CALL mpi_bc(this%dis_froz_proj, rank, mpi_comm)
    CALL mpi_bc(this%dis_proj_min, rank, mpi_comm)
    CALL mpi_bc(this%dis_proj_max, rank, mpi_comm)
    CALL mpi_bc(this%dis_mix_ratio, rank, mpi_comm)
    CALL mpi_bc(this%dis_conv_tol, rank, mpi_comm)
    CALL mpi_bc(this%conv_tol, rank, mpi_comm)
    CALL mpi_bc(this%precond, rank, mpi_comm)
    CALL mpi_bc(this%proj_ntype, rank, mpi_comm)
    CALL mpi_bc(this%proj_atom, rank, mpi_comm)
    CALL mpi_bc(this%proj_species, rank, mpi_comm)
    CALL mpi_bc(this%proj_l, rank, mpi_comm)
    CALL mpi_bc(this%proj_m, rank, mpi_comm)
    CALL mpi_bc(this%proj_spin, rank, mpi_comm)
    CALL mpi_bc(this%proj_rwf, rank, mpi_comm)
    CALL mpi_bc(this%proj_alpha, rank, mpi_comm)
    CALL mpi_bc(this%proj_beta, rank, mpi_comm)
    CALL mpi_bc(this%proj_gamma, rank, mpi_comm)
    CALL mpi_bc(this%proj_zona, rank, mpi_comm)
    CALL mpi_bc(this%proj_regio, rank, mpi_comm)
    CALL mpi_bc(this%proj_j, rank, mpi_comm)
    CALL mpi_bc(this%proj_mj, rank, mpi_comm)
    CALL mpi_bc(this%proj_weight, rank, mpi_comm)
    CALL mpi_bc(this%proj_shift, rank, mpi_comm)
  END SUBROUTINE mpi_bc_wannierlib

  SUBROUTINE read_xml_wannierlib(this, xml)
    USE m_types_xml
    USE m_constants
    CLASS(t_wannierlib_wannierize), INTENT(INOUT) :: this
    TYPE(t_xml), INTENT(INOUT) :: xml

    CHARACTER(LEN=200) :: xPathA
    CHARACTER(LEN=200) :: xPathP
    CHARACTER(LEN=32) :: spin_label
    CHARACTER(LEN=64) :: dom_name
    INTEGER :: idom, jdom, dom_npts
    CHARACTER(LEN=255) :: sbuf
    INTEGER :: numberNodes, ios
    INTEGER :: nSpecies, iType, nProjType, iProj, ip, nProjTotal, iRow
    LOGICAL :: has_wannierize
    CHARACTER(LEN=20) :: species_name

    xPathA = '/fleurInput/output/wannierlib'
    numberNodes = xml%getNumberOfNodes(xPathA)
    IF (numberNodes /= 1) RETURN

    this%l_wannierize = .TRUE.

    has_wannierize = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@wannierize') == 1
    IF (has_wannierize) THEN
      this%l_wannierize = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@wannierize'))
    END IF
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@plotWF') == 1) &
      this%l_plot_wf = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plotWF'))

    ! Operator selection is now the <interpolation>/<operator> list (parsed below),
    ! not flat @interpolation/@spinoperator/@orbmom/@socop attributes.

    xPathA = '/fleurInput/output/wannierlib/bands'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      this%num_bands = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numBands'))

      this%min_band = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minBand'))
      this%max_band = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxBand'))
    END IF

      
    xPathA = '/fleurInput/output/wannierlib/disentanglement'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      this%dis_win_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disWinMin'))
      this%dis_win_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disWinMax'))
      this%dis_froz_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disFrozMin'))
      this%dis_froz_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disFrozMax'))
      this%dis_froz_proj = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disFrozProj'))
      this%dis_proj_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disProjMin'))
      this%dis_proj_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disProjMax'))
      this%dis_num_iter = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numIter'))
      this%dis_mix_ratio = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mixRatio'))
      this%dis_conv_tol = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convTol'))
    END IF

    ! --- wannierization (MLWF iteration controls); defaults to the disentanglement values ---
    this%num_iter = this%dis_num_iter
    this%conv_tol = this%dis_conv_tol
    xPathA = '/fleurInput/output/wannierlib/wannierization'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      this%num_iter = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numIter'))
      this%conv_tol = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convTol'))
      this%precond = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@precond'))
    END IF

    ! --- interpolation domain + operator list ---
    xPathA = '/fleurInput/output/wannierlib/interpolation'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      ! --- output domains: one <domain> per set of k-points, repeatable ---
      this%n_domains = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/domain')
      IF (this%n_domains > 0) THEN
        ALLOCATE(this%dom_kset(this%n_domains))
        ALLOCATE(this%dom_suffix(this%n_domains))
        this%dom_suffix = ''
        DO idom = 1, this%n_domains
          WRITE(xPathP, '(a,i0,a)') TRIM(ADJUSTL(xPathA))//'/domain[', idom, ']'
          ! @listName is use="required" in the schema, so absence cannot reach here; an
          ! EMPTY name can, and would send read_kpts_by_name looking for a nameless list.
          dom_name = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@listName'))
          IF (LEN_TRIM(dom_name) == 0) &
            CALL juDFT_error('wannierlib: <domain> has an empty listName', &
                             calledby='read_xml_wannierlib')
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@suffix') == 1) &
            this%dom_suffix(idom) = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@suffix'))
          dom_npts = 0
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@npts') == 1) &
            dom_npts = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@npts'))
          ! Pulled now, while the XML is still open.
          CALL read_domain_kset_wannierlib(this%dom_kset(idom), xml, dom_name, dom_npts)
        END DO
        !> Two domains sharing a suffix would write over each other, and the unsuffixed one
        !> is the base name, so there is room for exactly one of those.
        DO idom = 1, this%n_domains - 1
          DO jdom = idom + 1, this%n_domains
            IF (TRIM(this%dom_suffix(idom)) == TRIM(this%dom_suffix(jdom))) &
              CALL juDFT_error('wannierlib: two <domain> entries share the suffix "'// &
                               TRIM(this%dom_suffix(idom))//'"', calledby='read_xml_wannierlib')
          END DO
        END DO
      END IF

      this%l_ws_distance = evaluateFirstBoolOnly( &
        xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@useWsDistance'))
      this%n_ops = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF

    ALLOCATE(this%op_name(this%n_ops), this%op_comp(this%n_ops), this%op_total(this%n_ops))
    DO iProj = 1, this%n_ops
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/output/wannierlib/interpolation/operator[', iProj, ']'
      this%op_name(iProj) = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@name'))
      this%op_comp(iProj) = ''
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@comp') == 1) &
         this%op_comp(iProj) = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@comp'))
      this%op_total(iProj) = 1
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@total') == 1) THEN
        IF (.NOT. evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@total'))) this%op_total(iProj) = 0
      END IF

      ! derive convenience flags that gate the coarse-mesh provider calls.
      ! Class-C operators auto-request their prerequisite coarse matrices:
      !   velocity is built from H alone, so it needs no coarse provider.
      !> Which coarse matrix a name needs built for it is a column of the exposure table,
      !> so a name added there is picked up here without being remembered twice. What is
      !> left below maps the three CATALOGUE entries onto the three flags, and that mapping
      !> only grows if the catalogue itself does.
      iRow = melem_exposed_find(this%op_name(iProj), WANNIERLIB_INTERP)
      IF (iRow > 0) THEN
        IF (TRIM(WANNIERLIB_INTERP(iRow)%name) == 'hamiltonian') this%l_interpolation = .TRUE.
        SELECT CASE (TRIM(WANNIERLIB_INTERP(iRow)%operator))
        CASE ('spin');       this%l_spin = .TRUE.
        CASE ('orbital');    this%l_orbmom = .TRUE.
        CASE ('spin_orbit'); this%l_socop = .TRUE.
        END SELECT
      END IF
    END DO

    ! --- operators_r: real-space operator matrices O(R) (Fourier step 3, no interpolation).
    !     Groups <operator name=".."> children written to standalone-format files. ---
    xPathA = '/fleurInput/output/wannierlib/operators_r'
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) == 1) THEN
      this%l_operators_r = .TRUE.
      this%n_op_r = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF
    ALLOCATE(this%op_r_name(this%n_op_r))
    DO iProj = 1, this%n_op_r
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/output/wannierlib/operators_r/operator[', iProj, ']'
      this%op_r_name(iProj) = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@name'))
      iRow = melem_exposed_find(this%op_r_name(iProj), WANNIERLIB_OPR)
      IF (iRow > 0) THEN
        SELECT CASE (TRIM(WANNIERLIB_OPR(iRow)%operator))
        CASE ('spin');       this%l_spin = .TRUE.
        CASE ('orbital');    this%l_orbmom = .TRUE.
        CASE ('spin_orbit'); this%l_socop = .TRUE.
        END SELECT
      END IF
    END DO

    nSpecies = xml%getNumberOfNodes('/fleurInput/atomSpecies/species')
    nProjTotal = 0
    DO iType = 1, nSpecies
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/atomSpecies/species[', iType, ']/wannierproj'
      nProjTotal = nProjTotal + xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP)))
    END DO

    IF (nProjTotal > 0) THEN
      ALLOCATE(this%proj_ntype(nProjTotal), this%proj_atom(nProjTotal), this%proj_species(nProjTotal), this%proj_l(nProjTotal), this%proj_m(nProjTotal), this%proj_spin(nProjTotal))
      ALLOCATE(this%proj_rwf(nProjTotal), this%proj_alpha(nProjTotal), this%proj_beta(nProjTotal), this%proj_gamma(nProjTotal))
      ALLOCATE(this%proj_zona(nProjTotal), this%proj_regio(nProjTotal), this%proj_j(nProjTotal), this%proj_mj(nProjTotal))
      ALLOCATE(this%proj_weight(nProjTotal), this%proj_shift(3, nProjTotal))
      ip = 0
      DO iType = 1, nSpecies
        WRITE(xPathA, '(A,I0,A)') '/fleurInput/atomSpecies/species[', iType, ']/@name'
        IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) == 1) THEN
          species_name = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))))
        ELSE
          CALL juDFT_error('wannierlib: species name missing while reading projections', calledby='read_xml_wannierlib')
        END IF

        WRITE(xPathP, '(A,I0,A)') '/fleurInput/atomSpecies/species[', iType, ']/wannierproj'
        nProjType = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP)))
        DO iProj = 1, nProjType
          ip = ip + 1
          WRITE(xPathP, '(A,I0,A,I0,A)') '/fleurInput/atomSpecies/species[', iType, ']/wannierproj[', iProj, ']'
       
          this%proj_species(ip) = species_name
          this%proj_l(ip) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@l'))
          this%proj_m(ip) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@m'))
          this%proj_rwf(ip) = 0
          this%proj_alpha(ip) = 0.0
          this%proj_beta(ip) = 0.0
          this%proj_gamma(ip) = 0.0
          this%proj_zona(ip) = 0.0
          this%proj_regio(ip) = 1.0
          this%proj_j(ip) = -1.0
          this%proj_mj(ip) = 0.0
          this%proj_weight(ip) = 1.0
          this%proj_shift(:, ip) = 0.0

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@rwf') == 1) THEN
            this%proj_rwf(ip) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@rwf'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@zona') == 1) THEN
            this%proj_zona(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@zona'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@regio') == 1) THEN
            this%proj_regio(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@regio'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@alpha') == 1) THEN
            this%proj_alpha(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@alpha'))
          ELSE IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@phi') == 1) THEN
            this%proj_alpha(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@phi'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@beta') == 1) THEN
            this%proj_beta(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@beta'))
          ELSE IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@theta') == 1) THEN
            this%proj_beta(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@theta'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@gamma') == 1) THEN
            this%proj_gamma(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@gamma'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@j') == 1) THEN
            this%proj_j(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@j'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@mj') == 1) THEN
            this%proj_mj(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@mj'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@weight') == 1) THEN
            this%proj_weight(ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@weight'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shiftX') == 1) THEN
            this%proj_shift(1, ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shiftX'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shiftY') == 1) THEN
            this%proj_shift(2, ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shiftY'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shiftZ') == 1) THEN
            this%proj_shift(3, ip) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shiftZ'))
          END IF

          spin_label = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@spin'))
          SELECT CASE (spin_label(1:1))
          CASE ('u', 'U')
            this%proj_spin(ip) = 1
          CASE ('d', 'D')
            this%proj_spin(ip) = -1
          CASE DEFAULT
            this%proj_spin(ip) = 0
          END SELECT
        END DO
      END DO
    END IF

  END SUBROUTINE read_xml_wannierlib

  !> One domain's k-set: the named kPointList, optionally with each segment subdivided into
  !> npts pieces. Read with the standard t_kpts routine, which follows the xi:include, so a
  !> list may live in kpts.xml or in any file included from inp.xml.
  !>
  !> Subdividing only makes sense for a list whose order is a path; on a plane or a mesh the
  !> caller leaves npts at its default and the list is taken as it stands.
  !>
  !> Runs wherever read_xml runs, but the result is consumed only on rank 0, so it is never
  !> broadcast -- see the note on dom_kset.
  SUBROUTINE read_domain_kset_wannierlib(kset, xml, listname, npts)
    USE m_types_xml
    TYPE(t_kpts), INTENT(OUT) :: kset
    TYPE(t_xml), INTENT(INOUT) :: xml
    CHARACTER(LEN=*), INTENT(IN) :: listname
    INTEGER, INTENT(IN) :: npts

    TYPE(t_kpts) :: raw_kset
    INTEGER :: nraw, iseg, isub, np, fac
    REAL, ALLOCATABLE :: raw(:, :)
    REAL :: t

    IF (.NOT. raw_kset%read_kpts_by_name(TRIM(xml%filename_add_xml)//"inp.xml", TRIM(listname))) &
      CALL juDFT_error('wannierlib: <domain>/@listName "'//TRIM(listname)// &
                       '" not found in kPointLists', calledby='read_domain_kset_wannierlib')
    nraw = raw_kset%nkpt
    IF (nraw < 1) CALL juDFT_error('wannierlib: <domain>/@listName "'//TRIM(listname)// &
                                   '" is empty', calledby='read_domain_kset_wannierlib')
    IF (nraw < 2 .AND. npts > 1) &
      CALL juDFT_error('wannierlib: @npts needs a list with at least two k-points to subdivide', &
                       calledby='read_domain_kset_wannierlib')
    ALLOCATE(raw(3, nraw))
    raw = raw_kset%bk(:, 1:nraw)

    fac = npts
    IF (fac < 1) fac = 1
    IF (fac == 1) THEN
      np = nraw
      ALLOCATE(kset%bk(3, np))
      kset%bk = raw
    ELSE
      ALLOCATE(kset%bk(3, (nraw - 1) * fac + 1))
      np = 0
      DO iseg = 1, nraw - 1
        DO isub = 0, fac - 1
          np = np + 1
          t = REAL(isub) / REAL(fac)
          kset%bk(:, np) = (1.0 - t) * raw(:, iseg) + t * raw(:, iseg + 1)
        END DO
      END DO
      np = np + 1
      kset%bk(:, np) = raw(:, nraw)
    END IF
    kset%nkpt = np
    DEALLOCATE(raw)
  END SUBROUTINE read_domain_kset_wannierlib


END MODULE m_types_wannierlib
