!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_types_wannierlib
  USE m_types_melem_optable, ONLY: WANNIERLIB_INTERP, WANNIERLIB_OPR, WANNIERLIB_OPK, &
                                  melem_exposed_find
  USE m_types_kpts
  USE m_juDFT
  USE m_types_atoms
  USE m_types_noco
  USE m_types_fleurinput_base
  USE m_constants,ONLY: ounit
  USE m_mpi_bc_tool
  USE m_types_xml
  USE m_calculator, ONLY: evaluateList
  USE m_constants
  IMPLICIT NONE
  PRIVATE

  !> One initial guess: the site it sits on, the orbital it is, and how it is oriented.
  !>
  !> SEQUENCE with fixed-size members only, like t_utype in t_atoms, because the list is
  !> broadcast field by field and an array of a type with allocatable components cannot be.
  !> The win is at the other end: init expands one entry per species into one per atom, and
  !> with a record that copies whole, a field added here cannot be silently dropped there.
  TYPE t_wannierlib_proj
    SEQUENCE
    INTEGER :: ntype = 0
    INTEGER :: atom = 1
    INTEGER :: l = 0
    INTEGER :: m = 0
    INTEGER :: spin = 0
    INTEGER :: rwf = 0
    REAL :: alpha = 0.0
    REAL :: beta = 0.0
    REAL :: gamma = 0.0
    REAL :: zona = 0.0
    REAL :: regio = 1.0
    REAL :: j = -1.0
    REAL :: mj = 0.0
    REAL :: weight = 1.0
    REAL :: shift(3) = 0.0
    CHARACTER(LEN=20) :: species = ''
  END TYPE t_wannierlib_proj

  !> One <operator> of <interpolation>: which operator to project onto the
  !> interpolated bands. Declaring it is what asks for it.
  TYPE t_wannierlib_op
    SEQUENCE
    CHARACTER(LEN=20) :: name = ''   ! hamiltonian/spin/orbital/spinOrbit/velocity/...
  END TYPE t_wannierlib_op

  !> The <export> block: copy out matrices the run already holds, for someone outside.
  !>
  !> PROVISIONAL. Nothing inside FLEUR reads any of these back. They exist so a run can be
  !> cross-checked against a standalone wannier90.x, and so tools that consume Wannier data
  !> -- WannierBerri, irrep -- can be fed without recomputing. Grouped rather than left
  !> as loose flags precisely because it is provisional: one object is one seam, whether
  !> the block is later promoted or dropped.
  TYPE t_wannierlib_export
    SEQUENCE
    !> The five files a standalone wannier90.x needs, as one set: the .win that declares
    !> the run, the three matrices, and the .chk Wannier90 writes itself in its own format.
    LOGICAL :: w90 = .FALSE.     ! <export wannier90="T"/>:  .win/.amn/.mmn/.eig/.chk
    !> Not a consumer but a fact about the run: the b-shells and the neighbour topology,
    !> fixed before any gauge exists and shared by everything built from them.
    LOGICAL :: bmesh = .FALSE.   ! <export bmesh="T"/>:          WF<spin>_bmesh.hdf
    LOGICAL :: orbitrans = .FALSE. ! <export orbitrans="T"/>:    input/, renamed, + four structural
    LOGICAL :: juwip = .FALSE.   ! <export juwip="T"/>:         w90files/, names unchanged
    !> The field keeps the content as its name and the attribute names the consumer, as
    !> w90 does for wannier90.
    LOGICAL :: basis = .FALSE.        ! <export irrep="T"/>:   WF<spin>_basis.hdf
    !> The one flag here that names no consumer: .xsf and .vts are read by XCrySDen and
    !> by ParaView, and by more than those two, so the file says what it holds instead.
    LOGICAL :: plot = .FALSE.    ! <export plotWF="T"/>:          WF<n>_<nnn>.xsf and .vts
    !> Cost is one wann_real per grid point per band per k-point, so this number is the
    !> whole price of the plot: 40^3 is 500x the work of 8^3. Left at 40^3, what it was
    !> before it could be asked for.
    INTEGER :: plot_grid(3) = [40, 40, 40]  ! <export plotGrid="40 40 40"/>
    !> How many cells the region spans. plot_region used to hardwire one, which cuts any
    !> function wider than the cell: in graphene (a = 2.46 A) a sigma function still holds
    !> half its maximum at the cell wall. The grid is scaled with it, so the resolution is
    !> the one plotGrid asked for and the cost grows with the product of the three.
    INTEGER :: plot_cells(3) = [1, 1, 1]    ! <export plotCells="1 1 1"/>
  END TYPE t_wannierlib_export

  TYPE, EXTENDS(t_fleurinput_base) :: t_wannierlib_wannierize
    LOGICAL :: l_wannierize = .FALSE.
    ! Convenience flags DERIVED from the operator list ops(:) below (set in read_xml).
    ! They gate the coarse-mesh provider calls; the actual dispatch loops over ops(:).
    LOGICAL :: l_interpolation = .FALSE.   ! an <operator name="hamiltonian"> is requested
    LOGICAL :: l_spin = .FALSE.            ! an <operator name="spin"> is requested
    LOGICAL :: l_orbmom = .FALSE.          ! an <operator name="orbital"> is requested
    LOGICAL :: l_socop = .FALSE.           ! an <operator name="spinOrbit"> is requested
    LOGICAL :: l_operators_r = .FALSE.     ! an <operatorsRs> block (real-space O(R) export) is present
    LOGICAL :: l_operators_k = .FALSE.     ! an <operatorsKs> block (Bloch basis, no gauge) is present
    !> The <export> block: one boolean per artefact, all off by default. They are grouped
    !> under one element because they answer the same question -- hand the run's own
    !> matrices to someone outside -- and differ only in who reads them. What <operatorsRs>
    !> and <interpolation> write is NOT here: those compute something first and writing it
    !> is the last step, while these copy out what the run already holds.
    TYPE(t_wannierlib_export) :: export
    ! --- opt-in output domains (<domain> children of <interpolation>) ---
    ! One <domain> per set of k-points the operators are evaluated on, repeatable and
    ! without kinds: @listName (required) names a kPointList of kpts.xml, @suffix tails the
    ! output names, @refinementPoints subdivides each segment of the list. Operators declared with no
    ! domain are an error rather than a silent no-op.
    !> One entry per <domain>. The count and the suffixes ARE broadcast -- the first sets
    !> the domain loop bound on every rank, the second is read inside that loop, where the
    !> output name is built. The k-sets are not: only rank 0 reaches the writers. Reading
    !> dom_kset off rank 0 is undefined behaviour, and a two-rank run can pass with it.
    INTEGER :: n_domains = 0
    TYPE(t_kpts), ALLOCATABLE :: dom_kset(:)
    CHARACTER(LEN=64), ALLOCATABLE :: dom_suffix(:)

    !> The path <export orbitrans="T"/> falls back to when the run declares no <domain>.
    !> Read here rather than in the export so the XML is opened once, and broadcast no more
    !> than dom_kset is: only rank 0 reaches the writer, and l_kpath stays false elsewhere.
    LOGICAL :: l_kpath = .FALSE.
    CHARACTER(LEN=64) :: kpath_name = ''
    TYPE(t_kpts) :: kpath_kset

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
    !> Always allocated, possibly to size zero, so SIZE() is the count and the two cannot
    !> come apart. One entry per <operator> of <interpolation>.
    TYPE(t_wannierlib_op), ALLOCATABLE :: ops(:)

    ! --- operatorsRs table: one entry per <operatorsRs>/<operator name=".."> child ---
    !> One name per <operator> of <operatorsRs>. A list of names is all this is, so it
    !> stays a list of names rather than a record with one field in it.
    CHARACTER(LEN=20), ALLOCATABLE :: op_r_name(:)

    ! --- operatorsKs table: one entry per <operatorsKs>/<operator name=".."> child.
    !     Kept apart from the list above for the same reason the request keeps them
    !     apart: neither one implies the other. ---
    CHARACTER(LEN=20), ALLOCATABLE :: op_k_name(:)

    INTEGER :: num_wann = 0
    INTEGER :: num_bands = 0
    INTEGER :: min_band = -1
    INTEGER :: max_band = -1

    !> The energy windows, in Hartree and absolute. Each may be left out of the input, in
    !> which case it is filled in from the bands themselves once they are known -- see
    !> wannierlib_default_windows.
    REAL :: dis_win_min = 0.0
    REAL :: dis_win_max = 0.0
    REAL :: dis_froz_min = 0.0
    REAL :: dis_froz_max = 0.0
    !> -1 and not 0 for "the input did not say", because for both counts ZERO IS A VALUE:
    !> the schema types each numIter as nonNegativeInteger, and 0 means do not minimise --
    !> for the wannierise count keep the gauge the projections give, for the disentanglement
    !> keep the subspace the projections and the frozen window define. With 0 as the
    !> sentinel neither of the two could be asked for.
    INTEGER :: dis_num_iter = -1
    INTEGER :: num_iter = -1     ! MLWF/wannierise iterations (W90 num_iter); XML @wannNumIter
    !> <transport>: Wannier90's ballistic conductance stage. The element's presence is the
    !> switch, so there is no separate flag to fall out of step with it. The axis is held as
    !> 1/2/3 rather than a character because every other field here is numeric and a scalar
    !> string would be the only one of its kind crossing mpi_bc; 0 means the input named no
    !> axis, which the adapter refuses rather than guessing.
    LOGICAL :: l_transport = .FALSE.
    INTEGER :: tran_axis = 0            ! 1=x, 2=y, 3=z; W90 one_dim_axis
    REAL    :: tran_win_min = 0.0       ! Hartree, absolute; min == max means "not given"
    REAL    :: tran_win_max = 0.0
    REAL    :: tran_energy_step = 0.0   ! Hartree; 0 means "not given"
    !> The two cutoffs that decide the principal layer: Wannier90 zeroes every hopping below
    !> hrCutoff and every one reaching beyond distCutoff, and the largest |R| left along the
    !> axis is the layer. Negative is the sentinel here too, because zero is a value for
    !> both: it is Wannier90's own default for hr_cutoff, and for dist_cutoff it would cut
    !> every hopping there is.
    REAL    :: tran_hr_cutoff = -1.0    ! Hartree; < 0 means "not given"
    REAL    :: tran_dist_cutoff = -1.0  ! Bohr;    < 0 means "not given"
    !> Projectability disentanglement (W90 dis_froz_proj / dis_proj_min / dis_proj_max).
    !> Off by default, as in Wannier90. The criterion is meaningful here because the
    !> projections are the atomic ones this code builds, which a reader of an arbitrary
    !> AMN file cannot assume.
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
    !> Wannier90 only tests convergence when its window is above one, and ITS TWO WINDOWS
    !> DO NOT START FROM THE SAME PLACE. For the MLWF iteration the reader forces -1 before
    !> reading the keyword (wannier90_readwrite.F90), so silence there means it never tests
    !> and runs the whole iteration count out. The disentanglement is not in that position:
    !> its own default is 3 and its reader rejects a negative one, so it does test.
    !> Zero here means we say nothing and leave each of them where it already was.
    INTEGER :: dis_conv_window = 0   ! disentanglement (W90 dis_conv_window); XML @convWindow
    INTEGER :: conv_window = 0       ! MLWF/wannierise (W90 conv_window);     XML @convWindow
    !> Preconditioned gradient in the spread minimisation (W90 precond); XML @precond.
    !> Off by default, which is also Wannier90's default.
    !>
    !> It changes the descent inside a valley, not the functional, so it cannot move the
    !> minimum.
    LOGICAL :: precond = .FALSE.
    !> The rest of Wannier90's minimisation controls. Each one is inert unless its
    !> attribute is present, so silence leaves Wannier90 on its own default and an older
    !> inp.xml keeps its numbers. None of them is stored in the checkpoint either, so a
    !> wannier90.x resuming from one takes them from the .win and not from here.
    INTEGER :: num_cg_steps = -1          ! W90 num_cg_steps;      XML @numCgSteps
    REAL    :: trial_step = -1.0          ! W90 trial_step;        XML @trialStep
    REAL    :: fixed_step = -1.0          ! W90 fixed_step;        XML @fixedStep
    LOGICAL :: guiding_centres = .FALSE.  ! W90 guiding_centres;   XML @guidingCentres
    INTEGER :: num_guide_cycles = -1      ! W90 num_guide_cycles;  XML @numGuideCycles
    INTEGER :: num_no_guide_iter = -1     ! W90 num_no_guide_iter; XML @numNoGuideIter
    REAL    :: conv_noise_amp = -1.0      ! W90 conv_noise_amp;    XML @convNoiseAmp
    INTEGER :: conv_noise_num = -1        ! W90 conv_noise_num;    XML @convNoiseNum
    !> Order of the finite-difference stencil for the k derivatives (W90 higher_order_n);
    !> XML @fdOrder. Zero is our sentinel for "say nothing", so Wannier90 keeps its own
    !> default of 1 -- the single-shell formula with an O(b^2) error. The quantities that
    !> feel it are the ones defined as k derivatives: A, and above all C and F, which are
    !> second derivatives and therefore carry the error squared.
    INTEGER :: fd_order = 0
    !> Tolerance Wannier90 groups the b-vectors into shells with, and checks the completeness
    !> relation sum_b w_b b_a b_b = delta_ab with -- Eq. (B1) of Marzari and Vanderbilt. A
    !> wannierisation runs on a uniform full-zone mesh, so when the shells come out wrong the
    !> lattice is the first suspect and not the mesh: a slightly distorted cell splits one
    !> shell into several, or fails (B1), and Wannier90 stops there with nothing left to
    !> turn. Nothing here checks the mesh, though -- mp_grid is taken from kpts%nkpt3 and the
    !> k-list from whichever list was selected -- so a shifted list would reach Wannier90 as
    !> a Monkhorst-Pack grid it is not, and would land in the same place. Negative is the
    !> sentinel for "say nothing", since zero is a value -- under it no two distances would
    !> ever count as equal and every shell would hold one vector.
    REAL :: kmesh_tol = -1.0

    LOGICAL :: l_intp = .FALSE.         ! Do interpolation
    REAL, ALLOCATABLE :: kpts_fine(:,:) ! kPoints to interpolate on 

    !> One entry per <wannierproj> as read, and one per Wannier function after init has
    !> expanded them over the atoms of each species.
    TYPE(t_wannierlib_proj), ALLOCATABLE :: proj(:)
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

    INTEGER :: i, j, nold, nnew, mcount, mval, mrepeat, spin_mult, itype, type_mult, nn, na
    INTEGER :: ipass, npass, spin_here
    TYPE(t_wannierlib_proj), ALLOCATABLE :: new_proj(:)
    LOGICAL :: expand_spin

    IF (.NOT.ALLOCATED(this%proj)) RETURN !No Wannier projections configured, nothing to do
    IF (.NOT.ALLOCATED(atoms%speciesName)) CALL juDFT_error('wannierlib: atoms%speciesName not allocated in init', calledby='init_wannierlib')
    
    nold = SIZE(this%proj)

    expand_spin = noco%l_noco .OR. noco%l_soc

    nnew = 0
    DO i = 1, nold
      mcount = projection_m_count(this%proj(i)%l)

      IF (this%proj(i)%m == 0) THEN
        IF (mcount <= 0) THEN
          CALL juDFT_error('wannierlib: unsupported projection l for m expansion', calledby='init_wannierlib')
        END IF
        mrepeat = mcount
      ELSE
        IF ((mcount > 0) .AND. ((this%proj(i)%m < 1) .OR. (this%proj(i)%m > mcount))) THEN
          CALL juDFT_error('wannierlib: projection m out of range for l', calledby='init_wannierlib')
        END IF
        mrepeat = 1
      END IF

      spin_mult = 1
      IF (expand_spin .AND. (this%proj(i)%spin == 0)) spin_mult = 2

      type_mult = 0
      DO itype = 1, atoms%ntype
        IF (TRIM(atoms%speciesName(itype)) == TRIM(this%proj(i)%species)) type_mult = type_mult + atoms%neq(itype)
      END DO
      IF (type_mult <= 0) THEN
        CALL juDFT_error('wannierlib: no matching atom type for projection species name', calledby='init_wannierlib')
      END IF

      nnew = nnew + mrepeat * spin_mult * type_mult
    END DO

    ALLOCATE(new_proj(nnew))

    !> Spin is the OUTER loop, not the inner one. The Wannier functions therefore come
    !> out ordered IN BLOCKS -- one channel first, then the other -- rather than interleaved
    !> (up, dn, up, dn, ...) as before. That gives the Wannier index the SAME convention the
    !> collinear 2N route uses, where FLEUR builds the basis as [channel 1, channel 2].
    !> Carrying two conventions in one code is a silent trap: slicing the first half of the
    !> index as "channel 1" is right in one route and scrambles the quadrants in the other,
    !> and nothing reports it.
    !>
    !> The SET of projections does not change, only its order, and that permutation reaches
    !> every operator derived from it alike: it is a similarity transformation, so no
    !> physical quantity moves. What does change byte for byte is the content of the output
    !> files of the spinor route.
    !>
    !> What this does NOT mean: with spin-orbit coupling spin is not a property of a Wannier
    !> function at all, so the blocks here order the TRIAL ORBITALS, not the character of the
    !> functions that come out.
    npass = MERGE(2, 1, expand_spin)
    j = 0
    DO ipass = 1, npass
    DO i = 1, nold
      mcount = projection_m_count(this%proj(i)%l)

      IF ((this%proj(i)%m == 0) .AND. (mcount > 0)) THEN
        mrepeat = mcount
      ELSE
        mrepeat = 1
      END IF

      spin_mult = 1
      IF (expand_spin .AND. (this%proj(i)%spin == 0)) spin_mult = 2

      IF (spin_mult == 2) THEN
        spin_here = MERGE(1, -1, ipass == 1)
      ELSE
        !> A projection whose spin the user fixed is emitted ONCE, on the pass of its
        !> own spin, so that it lands in the block it belongs to.
        spin_here = this%proj(i)%spin
        IF (npass > 1) THEN
          IF ((ipass == 1) .AND. (spin_here < 0)) CYCLE
          IF ((ipass == 2) .AND. (spin_here >= 0)) CYCLE
        END IF
      END IF

      DO itype = 1, atoms%ntype
        IF (TRIM(atoms%speciesName(itype)) /= TRIM(this%proj(i)%species)) CYCLE
        DO nn = 1, atoms%neq(itype)
          na = atoms%firstAtom(itype) + nn - 1
          DO mval = 1, mrepeat
            j = j + 1
            !> The whole record first, then only what the expansion changes. Copying field
            !> by field is how a field added to the type gets forgotten here.
            new_proj(j) = this%proj(i)
            new_proj(j)%ntype = itype
            new_proj(j)%atom = nn
            new_proj(j)%spin = spin_here
            ! With m=0 ("expand all m") the harmonic index is the 1-based mval, also when
            ! only one m exists (l=0/s: mrepeat=1): the literal m=0 is out of range for
            ! tlm(:,:,1:7) and would be read uninitialised.
            IF (this%proj(i)%m == 0) new_proj(j)%m = mval
          END DO
        END DO
      END DO
    END DO
    END DO

    CALL move_alloc(new_proj, this%proj)

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
           i, TRIM(this%proj(i)%species), this%proj(i)%ntype, this%proj(i)%atom, this%proj(i)%l, this%proj(i)%m, this%proj(i)%spin, this%proj(i)%rwf, &
           this%proj(i)%alpha, this%proj(i)%beta, this%proj(i)%gamma, this%proj(i)%zona, this%proj(i)%regio, this%proj(i)%j, this%proj(i)%mj, this%proj(i)%weight, &
           this%proj(i)%shift(1), this%proj(i)%shift(2), this%proj(i)%shift(3)
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
    CLASS(t_wannierlib_wannierize), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: mpi_comm
    INTEGER, INTENT(IN), OPTIONAL :: irank
    INTEGER :: rank, myrank, ierr, iproj, nproj, iop, nops

    rank = 0
    IF (PRESENT(irank)) rank = irank
    !> Needed by every block below that reallocates an array of records on the receivers.
    !> It is read in more than one place, so it is set once here and not next to the first
    !> use -- which is how it ended up undefined in the second.
    myrank = rank
#ifdef CPP_MPI
    CALL MPI_COMM_RANK(mpi_comm, myrank, ierr)
#endif

    CALL mpi_bc(this%l_wannierize, rank, mpi_comm)
    CALL mpi_bc(this%l_interpolation, rank, mpi_comm)
    CALL mpi_bc(this%l_ws_distance, rank, mpi_comm)
    CALL mpi_bc(this%l_spin, rank, mpi_comm)
    CALL mpi_bc(this%l_orbmom, rank, mpi_comm)
    CALL mpi_bc(this%l_socop, rank, mpi_comm)
    ! The domain COUNT is broadcast: it sets the loop bound on every rank. The k-SETS are
    ! not, because only rank 0 reaches the writers -- so on any other rank dom_kset is
    ! unallocated, and reading or validating it there is a bug, not a missing broadcast.
    ! The suffixes ARE broadcast; see below.
    CALL mpi_bc(this%n_domains, rank, mpi_comm)
    !> The suffixes ARE broadcast, unlike the k-sets: they are a handful of short strings,
    !> and the output name is built inside the domain loop, which every rank turns. Leaving
    !> them on rank 0 segfaulted rank 7 -- and not rank 1, so a two-rank suite passed.
    CALL mpi_bc(this%dom_suffix, rank, mpi_comm)
    nops = 0
    IF (ALLOCATED(this%ops)) nops = SIZE(this%ops)
    CALL mpi_bc(nops, rank, mpi_comm)
#ifdef CPP_MPI
    IF (myrank /= rank) THEN
      IF (ALLOCATED(this%ops)) DEALLOCATE(this%ops)
      ALLOCATE(this%ops(nops))
    END IF
#endif
    DO iop = 1, nops
      CALL mpi_bc(rank, mpi_comm, this%ops(iop)%name)
    END DO
    CALL mpi_bc(this%l_operators_r, rank, mpi_comm)
    CALL mpi_bc(this%export%w90, rank, mpi_comm)
    CALL mpi_bc(this%export%orbitrans, rank, mpi_comm)
    CALL mpi_bc(this%export%juwip, rank, mpi_comm)
    CALL mpi_bc(this%export%plot, rank, mpi_comm)
    CALL mpi_bc(rank, mpi_comm, this%export%plot_grid)
    CALL mpi_bc(rank, mpi_comm, this%export%plot_cells)
    CALL mpi_bc(this%export%basis, rank, mpi_comm)
    CALL mpi_bc(this%export%bmesh, rank, mpi_comm)
    CALL mpi_bc(this%l_transport, rank, mpi_comm)
    CALL mpi_bc(this%tran_axis, rank, mpi_comm)
    CALL mpi_bc(this%tran_win_min, rank, mpi_comm)
    CALL mpi_bc(this%tran_win_max, rank, mpi_comm)
    CALL mpi_bc(this%tran_energy_step, rank, mpi_comm)
    CALL mpi_bc(this%tran_hr_cutoff, rank, mpi_comm)
    CALL mpi_bc(this%tran_dist_cutoff, rank, mpi_comm)
    CALL mpi_bc(this%op_r_name, rank, mpi_comm)
    CALL mpi_bc(this%l_operators_k, rank, mpi_comm)
    CALL mpi_bc(this%op_k_name, rank, mpi_comm)
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
    CALL mpi_bc(this%dis_conv_window, rank, mpi_comm)
    CALL mpi_bc(this%conv_window, rank, mpi_comm)
    CALL mpi_bc(this%precond, rank, mpi_comm)
    CALL mpi_bc(this%num_cg_steps, rank, mpi_comm)
    CALL mpi_bc(this%trial_step, rank, mpi_comm)
    CALL mpi_bc(this%fixed_step, rank, mpi_comm)
    CALL mpi_bc(this%guiding_centres, rank, mpi_comm)
    CALL mpi_bc(this%num_guide_cycles, rank, mpi_comm)
    CALL mpi_bc(this%num_no_guide_iter, rank, mpi_comm)
    CALL mpi_bc(this%conv_noise_amp, rank, mpi_comm)
    CALL mpi_bc(this%conv_noise_num, rank, mpi_comm)
    CALL mpi_bc(this%fd_order, rank, mpi_comm)
    CALL mpi_bc(this%kmesh_tol, rank, mpi_comm)
    !> An array of a derived type does not go through the generic mpi_bc, which sizes the
    !> receivers from the message itself. The count travels first, then each record field by
    !> field -- the shape t_atoms uses for lda_u.
    nproj = 0
    IF (ALLOCATED(this%proj)) nproj = SIZE(this%proj)
    CALL mpi_bc(nproj, rank, mpi_comm)
#ifdef CPP_MPI
    IF (myrank /= rank) THEN
      IF (ALLOCATED(this%proj)) DEALLOCATE(this%proj)
      IF (nproj > 0) ALLOCATE(this%proj(nproj))
    END IF
#endif
    DO iproj = 1, nproj
      CALL mpi_bc(this%proj(iproj)%ntype, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%atom, rank, mpi_comm)
      !> Fixed-size members take the other argument order of the generic: the
      !> allocatable forms are mpi_bc(x, rank, comm), these are mpi_bc(rank, comm, x).
      CALL mpi_bc(rank, mpi_comm, this%proj(iproj)%species)
      CALL mpi_bc(this%proj(iproj)%l, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%m, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%spin, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%rwf, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%alpha, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%beta, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%gamma, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%zona, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%regio, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%j, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%mj, rank, mpi_comm)
      CALL mpi_bc(this%proj(iproj)%weight, rank, mpi_comm)
      CALL mpi_bc(rank, mpi_comm, this%proj(iproj)%shift)
    END DO
    CALL mpi_bc(this%kpts_fine, rank, mpi_comm)
    CALL mpi_bc(this%l_intp, rank, mpi_comm)

  END SUBROUTINE mpi_bc_wannierlib

  SUBROUTINE read_xml_wannierlib(this, xml)
    CLASS(t_wannierlib_wannierize), INTENT(INOUT) :: this
    TYPE(t_xml), INTENT(INOUT) :: xml

    CHARACTER(LEN=200) :: xPathA
    CHARACTER(LEN=200) :: xPathP
    CHARACTER(LEN=64) :: dom_name
    INTEGER :: idom, jdom, dom_npts
    CHARACTER(LEN=255) :: sbuf
    INTEGER :: numberNodes, ios
    INTEGER :: nSpecies, iType, nProjType, iProj, ip, nProjTotal, iRow
    INTEGER :: nOpsRead = 0, nOprRead = 0, nOpkRead = 0
    LOGICAL :: has_wannierize
    CHARACTER(LEN=20) :: species_name
    TYPE(t_kpts) :: kpts_temp
    CHARACTER(len=40) :: kptsPath
    REAL, ALLOCATABLE :: gridlist(:)


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
  
    xPathA = '/fleurInput/output/wannierlib/interpolation'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      this%l_intp = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_intp'))

      IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@kptsPath') == 1) THEN
        kptsPath = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kptsPath')))
        ! Initialize kpts_temp with the number of k-points from the kpts.xml file
        kpts_temp%nkpt = xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/kPointLists/kPointList[@name="'//TRIM(kptsPath)//'"]/kPoint')
        IF (kpts_temp%nkpt > 0) THEN
          ALLOCATE(kpts_temp%bk(3, kpts_temp%nkpt))
          ALLOCATE(kpts_temp%wtkpt(kpts_temp%nkpt))
          IF (kpts_temp%read_kpts_by_name(trim(xml%filename_add_xml)//"inp.xml", kptsPath)) THEN
            IF (ALLOCATED(this%kpts_fine)) DEALLOCATE(this%kpts_fine)
              ALLOCATE(this%kpts_fine(3, kpts_temp%nkpt))
              this%kpts_fine = kpts_temp%bk
            END IF
        END IF
      END IF
    END IF
  

      
    xPathA = '/fleurInput/output/wannierlib/disentanglement'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      !> Each window is optional. "Not given" is carried by the value itself rather than by
      !> a flag: an absent window keeps its 0.0 default, and dis_win_min == dis_win_max is
      !> what the adapter already reads as "no window". Adding fields to this type is not
      !> free -- it forces most of the tree to recompile, and ifx 2025 hits an internal
      !> compiler error in dfpt_interpolation when it does.
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@disWinMin') == 1) &
        this%dis_win_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disWinMin'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@disWinMax') == 1) &
        this%dis_win_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disWinMax'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@disFrozMin') == 1) &
        this%dis_froz_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disFrozMin'))
      this%dis_froz_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disFrozMax'))
      this%dis_froz_proj = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disFrozProj'))
      this%dis_proj_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disProjMin'))
      this%dis_proj_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disProjMax'))
      this%dis_num_iter = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numIter'))
      this%dis_mix_ratio = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mixRatio'))
      this%dis_conv_tol = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convTol'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@convWindow') == 1) &
        this%dis_conv_window = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convWindow'))
    END IF

    ! --- wannierization (MLWF iteration controls); defaults to the disentanglement values ---
    this%num_iter = this%dis_num_iter
    this%conv_tol = this%dis_conv_tol
    this%conv_window = this%dis_conv_window
    xPathA = '/fleurInput/output/wannierlib/wannierization'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      this%num_iter = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numIter'))
      this%conv_tol = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convTol'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@convWindow') == 1) &
        this%conv_window = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convWindow'))
      this%precond = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@precond'))
      !> Like convWindow: read only when the attribute is present, so that silence
      !> leaves Wannier90 on its own default and an older inp.xml keeps its numbers.
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@fdOrder') == 1) &
        this%fd_order = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@fdOrder'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@kmeshTol') == 1) &
        this%kmesh_tol = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@kmeshTol'))

      !> Same rule for all of these: read only when present.
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@numCgSteps') == 1) &
        this%num_cg_steps = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numCgSteps'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@trialStep') == 1) &
        this%trial_step = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@trialStep'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@fixedStep') == 1) &
        this%fixed_step = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@fixedStep'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@guidingCentres') == 1) &
        this%guiding_centres = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@guidingCentres'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@numGuideCycles') == 1) &
        this%num_guide_cycles = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numGuideCycles'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@numNoGuideIter') == 1) &
        this%num_no_guide_iter = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numNoGuideIter'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@convNoiseAmp') == 1) &
        this%conv_noise_amp = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convNoiseAmp'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@convNoiseNum') == 1) &
        this%conv_noise_num = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@convNoiseNum'))
    END IF

    ! --- transport (Landauer); the element being there is what switches it on ---
    xPathA = '/fleurInput/output/wannierlib/transport'
    IF (xml%getNumberOfNodes(xPathA) == 1) THEN
      this%l_transport = .TRUE.
      !> The schema restricts this to x, y or z, so the default arm cannot be reached from a
      !> validated input. It is left to fall through to 0 all the same, and the adapter says
      !> so: an axis that never arrived must not become a silent guess at the direction.
      SELECT CASE (TRIM(ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@axis'))))
      CASE ('x')
        this%tran_axis = 1
      CASE ('y')
        this%tran_axis = 2
      CASE ('z')
        this%tran_axis = 3
      END SELECT
      !> Both ends of the window, or neither. With exactly one given the other keeps its
      !> default of zero, the pair is still unequal, and it is passed on as a window running
      !> from eight electron volts down to nothing: Wannier90 takes the floor of a negative
      !> number of steps, never enters its energy loop, and writes two files holding their
      !> header and no data while reporting success. Equal ends already mean "say nothing",
      !> so only the half-given pair has to be refused.
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@energyMin') /= &
          xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@energyMax')) &
        CALL juDFT_error('wannierlib: <transport> needs energyMin and energyMax together', &
                         calledby='read_xml_wannierlib')
      !> Same convention as the disentanglement windows: absent means equal, equal means
      !> "say nothing", and silence leaves Wannier90 on its own range.
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@energyMin') == 1) &
        this%tran_win_min = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyMin'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@energyMax') == 1) &
        this%tran_win_max = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyMax'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@energyStep') == 1) &
        this%tran_energy_step = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyStep'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@hrCutoff') == 1) &
        this%tran_hr_cutoff = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@hrCutoff'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@distCutoff') == 1) &
        this%tran_dist_cutoff = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@distCutoff'))
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
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@refinementPoints') == 1) &
            dom_npts = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@refinementPoints'))
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
      nOpsRead = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF

    ALLOCATE(this%ops(nOpsRead))
    DO iProj = 1, nOpsRead
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/output/wannierlib/interpolation/operator[', iProj, ']'
      this%ops(iProj)%name = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@name'))

      ! derive convenience flags that gate the coarse-mesh provider calls.
      ! Class-C operators auto-request their prerequisite coarse matrices:
      !   velocity is built from H alone, so it needs no coarse provider.
      !> Which coarse matrix a name needs built for it is a column of the exposure table,
      !> so a name added there is picked up here without being remembered twice. What is
      !> left below maps the three CATALOGUE entries onto the three flags, and that mapping
      !> only grows if the catalogue itself does.
      iRow = melem_exposed_find(this%ops(iProj)%name, WANNIERLIB_INTERP)
      IF (iRow > 0) THEN
        IF (TRIM(WANNIERLIB_INTERP(iRow)%name) == 'hamiltonian') this%l_interpolation = .TRUE.
        SELECT CASE (TRIM(WANNIERLIB_INTERP(iRow)%operator))
        CASE ('spin');       this%l_spin = .TRUE.
        CASE ('orbital');    this%l_orbmom = .TRUE.
        CASE ('spinOrbit'); this%l_socop = .TRUE.
        END SELECT
      END IF
    END DO

    ! --- export: hand the run's own matrices to a reader outside. One boolean per
    !     artefact rather than one format enumeration, because they are not alternatives:
    !     a run can want the .amn/.mmn to be replayed by wannier90.x AND the interstitial
    !     for irrep AND the gauge, and asking for one must not turn the others off. ---
    xPathA = '/fleurInput/output/wannierlib/export'
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) == 1) THEN
      this%export%w90 = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@wannier90'))
      this%export%orbitrans = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@orbitrans'))
      this%export%juwip = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@juwip'))
      this%export%plot = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@plotWF'))
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@plotGrid') == 1) THEN
        !> Through a variable and not straight from the call: evaluateList takes its
        !> string as INTENT(INOUT), as types_sliceplot does for the same attribute.
        sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plotGrid')
        CALL evaluateList(gridlist, sbuf)
        IF (SIZE(gridlist) /= 3) CALL juDFT_error( &
          'wannierlib: <export plotGrid=".."/> takes three numbers', &
          calledby='read_xml_wannierlib')
        !> Two points per direction is the least that spans the cell at all: the step is
        !> v/(n-1), so n=1 divides by zero and n<1 is not a grid.
        IF (ANY(NINT(gridlist) < 2)) CALL juDFT_error( &
          'wannierlib: <export plotGrid=".."/> needs at least two points per direction', &
          calledby='read_xml_wannierlib')
        this%export%plot_grid = NINT(gridlist)
      END IF
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@plotCells') == 1) THEN
        sbuf = xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plotCells')
        CALL evaluateList(gridlist, sbuf)
        IF (SIZE(gridlist) /= 3) CALL juDFT_error( &
          'wannierlib: <export plotCells=".."/> takes three numbers', &
          calledby='read_xml_wannierlib')
        !> One cell is the whole cell; zero or less is not a region.
        IF (ANY(NINT(gridlist) < 1)) CALL juDFT_error( &
          'wannierlib: <export plotCells=".."/> needs at least one cell per direction', &
          calledby='read_xml_wannierlib')
        this%export%plot_cells = NINT(gridlist)
      END IF
      this%export%basis = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@irrep'))
      this%export%bmesh = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@bmesh'))
      !> wannierberri="T" is a superset of the two flags above rather than a third artefact:
      !> that wannierization consumes the .amn/.mmn/.eig/.chk and takes its symmetry table
      !> from what irrep reads out of the interstitial, so one flag asks for everything the
      !> route needs and nobody has to remember to light two.
      IF (evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@wannierberri'))) THEN
        this%export%w90 = .TRUE.
        this%export%basis = .TRUE.
      END IF
    END IF

    !> No <domain>, but somebody downstream still has to be told which k to resolve on.
    !> Tied to the orbitrans flag and not made general because orbitrans is the only reader
    !> that opens input/kpts today: a run exporting for anyone else gains no file and pays
    !> no time, and the day a second consumer wants a path it asks here by name.
    !>
    !> The path is not invented. It is the first list the input already carries with
    !> type="path", which inpgen writes for the higher-symmetry Bravais lattices and for
    !> films but not for a monoclinic or triclinic cell, so the default names something the
    !> caller can read off kpts.xml, and a run with no path list keeps today's silence.
    IF (this%export%orbitrans .AND. this%n_domains < 1) THEN
      xPathA = '/fleurInput/cell/bzIntegration/kPointLists/kPointList[@type="path"]'
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) > 0) THEN
        dom_name = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'[1]/@name'))
        IF (LEN_TRIM(dom_name) > 0) THEN
          CALL read_domain_kset_wannierlib(this%kpath_kset, xml, dom_name, 0)
          this%kpath_name = dom_name
          this%l_kpath = .TRUE.
        END IF
      END IF
    END IF

    ! --- operatorsKs: what the run holds in k, before the Fourier sum to R. The four
    !     Bloch-basis matrices come straight out of the coarse pass and need no gauge, so
    !     asking only for those costs no wannierisation; bmn, cmn and fmn are pair
    !     contractions and do carry the gauge, so they go to a file of their own. ---
    xPathA = '/fleurInput/output/wannierlib/operatorsKs'
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) == 1) THEN
      this%l_operators_k = .TRUE.
      nOpkRead = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF
    ALLOCATE(this%op_k_name(nOpkRead))
    DO iProj = 1, nOpkRead
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/output/wannierlib/operatorsKs/operator[', iProj, ']'
      this%op_k_name(iProj) = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@name'))
      iRow = melem_exposed_find(this%op_k_name(iProj), WANNIERLIB_OPK)
      IF (iRow > 0) THEN
        SELECT CASE (TRIM(WANNIERLIB_OPK(iRow)%operator))
        CASE ('spin');      this%l_spin = .TRUE.
        CASE ('orbital');   this%l_orbmom = .TRUE.
        CASE ('spinOrbit'); this%l_socop = .TRUE.
        END SELECT
      END IF
    END DO

    ! --- operatorsRs: real-space operator matrices O(R) (Fourier step 3, no interpolation).
    !     Groups <operator name=".."> children written to standalone-format files. ---
    xPathA = '/fleurInput/output/wannierlib/operatorsRs'
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) == 1) THEN
      this%l_operators_r = .TRUE.
      nOprRead = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF
    ALLOCATE(this%op_r_name(nOprRead))
    DO iProj = 1, nOprRead
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/output/wannierlib/operatorsRs/operator[', iProj, ']'
      this%op_r_name(iProj) = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@name'))
      iRow = melem_exposed_find(this%op_r_name(iProj), WANNIERLIB_OPR)
      IF (iRow > 0) THEN
        SELECT CASE (TRIM(WANNIERLIB_OPR(iRow)%operator))
        CASE ('spin');       this%l_spin = .TRUE.
        CASE ('orbital');    this%l_orbmom = .TRUE.
        CASE ('spinOrbit'); this%l_socop = .TRUE.
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
      !> The defaults live in t_wannierlib_proj, so an entry the input leaves out is
      !> already what it should be and does not need setting here.
      ALLOCATE(this%proj(nProjTotal))
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
       
          this%proj(ip)%species = species_name
          this%proj(ip)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@l'))
          this%proj(ip)%m = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@m'))
          !> The spin a projection belongs to. It was declared in the schema, allocated and
          !> broadcast, but never read: proj_spin reached the expansion uninitialised and a
          !> projection could not be pinned to one channel. That left the (l, j, m_j) branch
          !> unreachable from an input file -- with spinors every projection is doubled, both
          !> copies carry the same (l, j, m_j), and the j-resolved builder ignores the spin
          !> index, so the pair came out as two identical columns and A(k) lost rank.
          !> 0 means "either", which is what makes the expansion split it in two.
          !> The only place @spin is read. A second, unguarded read used to follow the
          !> rest of the attributes and overwrite this one by its first character, which
          !> turned the accepted "1", "+1" and "-1" back into 0 -- a projection asked for
          !> on one channel came out on both. It also asked for the attribute whether or
          !> not it was there, and an absent xPath goes through judft_warn, which is fatal
          !> unless -warn_only: a <wannierproj> written without spin= aborted the run.
          this%proj(ip)%spin = 0
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@spin') == 1) THEN
            sbuf = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@spin'))
            SELECT CASE (TRIM(sbuf))
            CASE (''); this%proj(ip)%spin = 0
            CASE ('up', '1', '+1'); this%proj(ip)%spin = 1
            CASE ('down', '-1'); this%proj(ip)%spin = -1
            CASE DEFAULT
              CALL juDFT_error('wannierlib: <wannierproj> spin="'//TRIM(sbuf)// &
                               '" is not one of up/down (or empty for both)', &
                               calledby='read_xml_wannierlib')
            END SELECT
          END IF
          this%proj(ip)%rwf = 0
          this%proj(ip)%alpha = 0.0
          this%proj(ip)%beta = 0.0
          this%proj(ip)%gamma = 0.0
          this%proj(ip)%zona = 0.0
          this%proj(ip)%regio = 1.0
          this%proj(ip)%j = -1.0
          this%proj(ip)%mj = 0.0
          this%proj(ip)%weight = 1.0
          this%proj(ip)%shift(:) = 0.0

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@rwf') == 1) THEN
            this%proj(ip)%rwf = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@rwf'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@zona') == 1) THEN
            this%proj(ip)%zona = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@zona'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@regio') == 1) THEN
            this%proj(ip)%regio = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@regio'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@alpha') == 1) THEN
            this%proj(ip)%alpha = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@alpha'))
          ELSE IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@phi') == 1) THEN
            this%proj(ip)%alpha = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@phi'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@beta') == 1) THEN
            this%proj(ip)%beta = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@beta'))
          ELSE IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@theta') == 1) THEN
            this%proj(ip)%beta = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@theta'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@gamma') == 1) THEN
            this%proj(ip)%gamma = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@gamma'))
          END IF

          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@j') == 1) THEN
            this%proj(ip)%j = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@j'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@mj') == 1) THEN
            this%proj(ip)%mj = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@mj'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@weight') == 1) THEN
            this%proj(ip)%weight = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@weight'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shiftX') == 1) THEN
            this%proj(ip)%shift(1) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shiftX'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shiftY') == 1) THEN
            this%proj(ip)%shift(2) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shiftY'))
          END IF
          IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@shiftZ') == 1) THEN
            this%proj(ip)%shift(3) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@shiftZ'))
          END IF
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
    TYPE(t_kpts), INTENT(OUT) :: kset
    TYPE(t_xml), INTENT(INOUT) :: xml
    CHARACTER(LEN=*), INTENT(IN) :: listname
    INTEGER, INTENT(IN) :: npts

    TYPE(t_kpts) :: raw_kset
    INTEGER :: nraw, iseg, isub, np, fac, isp
    REAL, ALLOCATABLE :: raw(:, :)
    REAL :: t

    !> l_labels: the interpolated band files carry the high-symmetry points, so whoever
    !> plots them does not have to work them out from kpts.xml a second time.
    IF (.NOT. raw_kset%read_kpts_by_name(TRIM(xml%filename_add_xml)//"inp.xml", &
                                         TRIM(listname), l_labels=.TRUE.)) &
      CALL juDFT_error('wannierlib: <domain>/@listName "'//TRIM(listname)// &
                       '" not found in kPointLists', calledby='read_domain_kset_wannierlib')
    nraw = raw_kset%nkpt
    IF (nraw < 1) CALL juDFT_error('wannierlib: <domain>/@listName "'//TRIM(listname)// &
                                   '" is empty', calledby='read_domain_kset_wannierlib')
    IF (nraw < 2 .AND. npts > 1) &
      CALL juDFT_error('wannierlib: @refinementPoints needs a list with at least two k-points to subdivide', &
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

    !> The high-symmetry labels the list carries, kept across the subdivision: a special
    !> point at raw index i lands at (i-1)*fac + 1, and the last one at np.
    kset%numSpecialPoints = raw_kset%numSpecialPoints
    IF (raw_kset%numSpecialPoints > 0) THEN
      ALLOCATE(kset%specialPointIndices(raw_kset%numSpecialPoints))
      ALLOCATE(kset%specialPointNames(raw_kset%numSpecialPoints))
      ALLOCATE(kset%specialPoints(3, raw_kset%numSpecialPoints))
      kset%specialPointNames = raw_kset%specialPointNames
      kset%specialPoints     = raw_kset%specialPoints
      DO isp = 1, raw_kset%numSpecialPoints
        kset%specialPointIndices(isp) = (raw_kset%specialPointIndices(isp) - 1) * fac + 1
      END DO
    END IF

    DEALLOCATE(raw)
  END SUBROUTINE read_domain_kset_wannierlib


END MODULE m_types_wannierlib
