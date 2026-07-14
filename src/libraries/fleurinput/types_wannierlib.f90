!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_types_wannierlib
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

    ! --- opt-in output domains (<path>/<plane>/<grid> children of <interpolation>) ---
    ! Each declared domain interpolates the operators on its own k-set and writes
    ! bands_wann_*[_plane|_grid].dat. If NONE is declared, the legacy single-pass
    ! behaviour (read kpts_interpol, unsuffixed output) is used -> byte-identical.
    LOGICAL :: l_dom_path  = .FALSE.       ! <path/>  1D band path (from path_file)
    LOGICAL :: l_dom_plane = .FALSE.       ! <plane/> 2D plane origin + i*v1 + j*v2
    LOGICAL :: l_dom_grid  = .FALSE.       ! <grid/>  3D uniform full-BZ mesh
    CHARACTER(LEN=64) :: path_file = 'kpts_interpol'
    REAL    :: plane_origin(3) = 0.0, plane_v1(3) = 0.0, plane_v2(3) = 0.0
    INTEGER :: plane_n1 = 0, plane_n2 = 0
    INTEGER :: grid_mesh(3) = 0
    REAL    :: grid_shift(3) = 0.0
    ! <path listName=".." npts="F"/>: reuse an existing named kPointList from kpts.xml
    ! (e.g. the band path "path-2") as the interpolation k-set, subdividing each segment
    ! into F pieces (npts<=1 -> use the list points as-is). listName takes precedence over
    ! @file. The k-points are extracted from the XML in read_xml (into path_kpts) and are
    ! consumed ONLY on rank 0 in write_domain_kpts -> not broadcast (like plane/grid).
    CHARACTER(LEN=64) :: path_listname = ''
    INTEGER :: path_npts = 0
    INTEGER :: path_np = 0
    REAL, ALLOCATABLE :: path_kpts(:, :)
    ! <plane listName=".."/> / <grid listName=".."/>: reuse a named kPointList (as-is)
    ! from kpts.xml as the plane/grid k-set instead of generating it inline. Filled in
    ! read_xml (into plane_kpts/grid_kpts), consumed on rank 0 in write_domain_kpts.
    CHARACTER(LEN=64) :: plane_listname = ''
    INTEGER :: plane_np = 0
    REAL, ALLOCATABLE :: plane_kpts(:, :)
    CHARACTER(LEN=64) :: grid_listname = ''
    INTEGER :: grid_np = 0
    REAL, ALLOCATABLE :: grid_kpts(:, :)

    ! --- operator table: one entry per <operator name=".."> child of <interpolation> ---
    INTEGER :: n_ops = 0
    CHARACTER(LEN=20), ALLOCATABLE :: op_name(:)    ! operator identifier (hamiltonian/spin/orbital/soc/...)
    CHARACTER(LEN=32), ALLOCATABLE :: op_comp(:)    ! requested components (Phase 2); '' = all of the operator rank
    INTEGER, ALLOCATABLE :: op_total(:)             ! 1 = write summed-over-atoms projection (default)
    INTEGER, ALLOCATABLE :: op_peratom(:)           ! 1 = also write per-atom (site-resolved) projection (Phase 2)

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
    REAL :: dis_mix_ratio = 0.0
    REAL :: dis_conv_tol = 0.0
    REAL :: conv_tol = 0.0       ! MLWF/wannierise convergence (W90 conv_tol); XML @wannConvTol

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
    CALL mpi_bc(this%l_spin, rank, mpi_comm)
    CALL mpi_bc(this%l_orbmom, rank, mpi_comm)
    CALL mpi_bc(this%l_socop, rank, mpi_comm)
    ! Output-domain flags: only these must agree on all ranks (they set the domain-loop
    ! count ndom). The plane/grid params (fixed-size arrays) and path_file are consumed
    ! ONLY on rank 0 (which writes kpts_interpol and does all domain file I/O), so they
    ! are not broadcast -- and mpi_bc has no specific for explicit-shape arrays anyway.
    CALL mpi_bc(this%l_dom_path, rank, mpi_comm)
    CALL mpi_bc(this%l_dom_plane, rank, mpi_comm)
    CALL mpi_bc(this%l_dom_grid, rank, mpi_comm)
    CALL mpi_bc(this%n_ops, rank, mpi_comm)
    CALL mpi_bc(this%op_name, rank, mpi_comm)
    CALL mpi_bc(this%op_comp, rank, mpi_comm)
    CALL mpi_bc(this%op_total, rank, mpi_comm)
    CALL mpi_bc(this%op_peratom, rank, mpi_comm)
    CALL mpi_bc(this%l_operators_r, rank, mpi_comm)
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
    CALL mpi_bc(this%dis_mix_ratio, rank, mpi_comm)
    CALL mpi_bc(this%dis_conv_tol, rank, mpi_comm)
    CALL mpi_bc(this%conv_tol, rank, mpi_comm)
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
    CHARACTER(LEN=255) :: sbuf
    INTEGER :: numberNodes, ios
    INTEGER :: nSpecies, iType, nProjType, iProj, ip, nProjTotal
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
    END IF

    ! --- interpolation domain + operator list ---
    xPathA = '/fleurInput/output/wannierlib/interpolation'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      ! --- opt-in output domains: <path>/<plane>/<grid> children of <interpolation> ---
      xPathP = TRIM(ADJUSTL(xPathA))//'/path'
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))) == 1) THEN
        this%l_dom_path = .TRUE.
        IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@file') == 1) &
          this%path_file = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@file'))
        IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@listName') == 1) &
          this%path_listname = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@listName'))
        IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@npts') == 1) &
          this%path_npts = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@npts'))
        ! <path listName=".."> : pull the named list's k-points now (XML is available here)
        IF (LEN_TRIM(this%path_listname) > 0) CALL read_named_kpath_wannierlib(this, xml)
      END IF
      xPathP = TRIM(ADJUSTL(xPathA))//'/plane'
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))) == 1) THEN
        this%l_dom_plane = .TRUE.
        IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@listName') == 1) THEN
          ! reuse a named kPointList from kpts.xml as-is
          this%plane_listname = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@listName'))
          CALL read_named_klist_wannierlib(xml, this%plane_listname, this%plane_kpts, this%plane_np)
        ELSE
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@origin') == 1) THEN
            sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@origin')
            READ(sbuf, *, iostat=ios) this%plane_origin
            IF (ios /= 0) CALL juDFT_error('wannierlib: <plane>/@origin needs 3 reals', calledby='read_xml_wannierlib')
          END IF
          sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@v1')
          READ(sbuf, *, iostat=ios) this%plane_v1
          IF (ios /= 0) CALL juDFT_error('wannierlib: <plane>/@v1 needs 3 reals', calledby='read_xml_wannierlib')
          sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@v2')
          READ(sbuf, *, iostat=ios) this%plane_v2
          IF (ios /= 0) CALL juDFT_error('wannierlib: <plane>/@v2 needs 3 reals', calledby='read_xml_wannierlib')
          this%plane_n1 = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@n1'))
          this%plane_n2 = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@n2'))
        END IF
      END IF
      xPathP = TRIM(ADJUSTL(xPathA))//'/grid'
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))) == 1) THEN
        this%l_dom_grid = .TRUE.
        IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@listName') == 1) THEN
          this%grid_listname = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@listName'))
          CALL read_named_klist_wannierlib(xml, this%grid_listname, this%grid_kpts, this%grid_np)
        ELSE
          sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@mesh')
          READ(sbuf, *, iostat=ios) this%grid_mesh
          IF (ios /= 0) CALL juDFT_error('wannierlib: <grid>/@mesh needs 3 integers', calledby='read_xml_wannierlib')
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shift') == 1) THEN
            sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shift')
            READ(sbuf, *, iostat=ios) this%grid_shift
            IF (ios /= 0) CALL juDFT_error('wannierlib: <grid>/@shift needs 3 reals', calledby='read_xml_wannierlib')
          END IF
        END IF
      END IF

      this%n_ops = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF

    ALLOCATE(this%op_name(this%n_ops), this%op_comp(this%n_ops), this%op_total(this%n_ops), this%op_peratom(this%n_ops))
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
      this%op_peratom(iProj) = 0
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@perAtom') == 1) THEN
        IF (evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@perAtom'))) this%op_peratom(iProj) = 1
      END IF

      ! derive convenience flags that gate the coarse-mesh provider calls.
      ! Class-C operators auto-request their prerequisite coarse matrices:
      !   spinCurrent = {v,sigma} needs the spin coarse matrix; orbitalCurrent needs L.
      !   velocity is built from H alone, so it needs no coarse provider.
      SELECT CASE (TRIM(this%op_name(iProj)))
      CASE ('hamiltonian');    this%l_interpolation = .TRUE.
      CASE ('spin');           this%l_spin = .TRUE.
      CASE ('orbital');        this%l_orbmom = .TRUE.
      CASE ('soc');            this%l_socop = .TRUE.
      CASE ('spinCurrent');    this%l_spin = .TRUE.
      CASE ('orbitalCurrent'); this%l_orbmom = .TRUE.
      END SELECT
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
      SELECT CASE (TRIM(this%op_r_name(iProj)))
      CASE ('spin');       this%l_spin = .TRUE.
      CASE ('orbital');    this%l_orbmom = .TRUE.
      CASE ('spin_orbit'); this%l_socop = .TRUE.
      END SELECT
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

  ! Extract a named kPointList (e.g. the band path "path-2") from the XML tree and store it
  ! as the <path> interpolation k-set, subdividing each segment into path_npts pieces. Uses
  ! the same XML idiom as t_kpts%read_xml (evaluatefirst consumes the kx,ky,kz tokens, so it
  ! also handles fraction coordinates like "11.00/24.00"). Runs wherever read_xml runs; the
  ! result (path_kpts) is consumed only on rank 0, so it is not broadcast.
  SUBROUTINE read_named_kpath_wannierlib(this, xml)
    USE m_types_xml
    USE m_calculator, ONLY: evaluatefirst
    CLASS(t_wannierlib_wannierize), INTENT(INOUT) :: this
    TYPE(t_xml), INTENT(INOUT) :: xml

    CHARACTER(LEN=255) :: path, path2, str
    INTEGER :: nlist, i, idx, nraw, iseg, isub, np, fac
    REAL, ALLOCATABLE :: raw(:, :)
    REAL :: t

    ! resolve listName -> kPointList index (same scan as t_kpts%read_xml)
    nlist = xml%getNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList')
    idx = 0
    DO i = 1, nlist
      WRITE (path, '(a,i0,a)') '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', i, ']'
      IF (TRIM(ADJUSTL(this%path_listname)) == TRIM(ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(path))//'/@name')))) THEN
        idx = i
        EXIT
      END IF
    END DO
    IF (idx == 0) CALL juDFT_error('wannierlib: <path>/@listName "'//TRIM(this%path_listname)// &
                                   '" not found in kPointLists', calledby='read_named_kpath_wannierlib')

    WRITE (path, '(a,i0,a)') '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', idx, ']'
    nraw = xml%getNumberOfNodes(TRIM(ADJUSTL(path))//'/kPoint')
    IF (nraw < 2) CALL juDFT_error('wannierlib: <path>/@listName needs a list with >= 2 k-points', &
                                   calledby='read_named_kpath_wannierlib')
    ALLOCATE(raw(3, nraw))
    DO i = 1, nraw
      WRITE (path2, '(a,a,i0,a)') TRIM(ADJUSTL(path)), '/kPoint[', i, ']'
      str = xml%getAttributeValue(TRIM(ADJUSTL(path2)), .TRUE.)
      raw(1, i) = evaluatefirst(str)
      raw(2, i) = evaluatefirst(str)
      raw(3, i) = evaluatefirst(str)
    END DO

    fac = this%path_npts
    IF (fac < 1) fac = 1
    IF (fac == 1) THEN
      ! use the list points as-is
      np = nraw
      ALLOCATE(this%path_kpts(3, np))
      this%path_kpts = raw
    ELSE
      ! subdivide each segment into fac pieces: (nraw-1)*fac + 1 points total
      ALLOCATE(this%path_kpts(3, (nraw - 1) * fac + 1))
      np = 0
      DO iseg = 1, nraw - 1
        DO isub = 0, fac - 1
          np = np + 1
          t = REAL(isub) / REAL(fac)
          this%path_kpts(:, np) = (1.0 - t) * raw(:, iseg) + t * raw(:, iseg + 1)
        END DO
      END DO
      np = np + 1
      this%path_kpts(:, np) = raw(:, nraw)   ! closing endpoint
    END IF
    this%path_np = np
    DEALLOCATE(raw)
  END SUBROUTINE read_named_kpath_wannierlib

  SUBROUTINE read_named_klist_wannierlib(xml, listname, arr, np)
    ! Load a named kPointList from kpts.xml as-is (no subdivision) into arr(3,np).
    ! Used by <plane listName=".."/> and <grid listName=".."/>.
    USE m_types_xml
    USE m_calculator, ONLY: evaluatefirst
    TYPE(t_xml), INTENT(INOUT)       :: xml
    CHARACTER(LEN=*), INTENT(IN)     :: listname
    REAL, ALLOCATABLE, INTENT(OUT)   :: arr(:, :)
    INTEGER, INTENT(OUT)             :: np

    CHARACTER(LEN=255) :: path, path2, str
    INTEGER :: nlist, i, idx

    nlist = xml%getNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList')
    idx = 0
    DO i = 1, nlist
      WRITE (path, '(a,i0,a)') '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', i, ']'
      IF (TRIM(ADJUSTL(listname)) == TRIM(ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(path))//'/@name')))) THEN
        idx = i
        EXIT
      END IF
    END DO
    IF (idx == 0) CALL juDFT_error('wannierlib: listName "'//TRIM(ADJUSTL(listname))// &
                                   '" not found in kPointLists', calledby='read_named_klist_wannierlib')

    WRITE (path, '(a,i0,a)') '/fleurInput/cell/bzIntegration/kPointLists/kPointList[', idx, ']'
    np = xml%getNumberOfNodes(TRIM(ADJUSTL(path))//'/kPoint')
    IF (np < 1) CALL juDFT_error('wannierlib: listName "'//TRIM(ADJUSTL(listname))// &
                                 '" has no k-points', calledby='read_named_klist_wannierlib')
    ALLOCATE(arr(3, np))
    DO i = 1, np
      WRITE (path2, '(a,a,i0,a)') TRIM(ADJUSTL(path)), '/kPoint[', i, ']'
      str = xml%getAttributeValue(TRIM(ADJUSTL(path2)), .TRUE.)
      arr(1, i) = evaluatefirst(str)
      arr(2, i) = evaluatefirst(str)
      arr(3, i) = evaluatefirst(str)
    END DO
  END SUBROUTINE read_named_klist_wannierlib

END MODULE m_types_wannierlib
