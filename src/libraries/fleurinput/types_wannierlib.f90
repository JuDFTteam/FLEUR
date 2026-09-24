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

  !> One <operator> of <interpolation>: which operator, which components of it, and
  !> whether the site-summed projection is wanted.
  TYPE t_wannierlib_op
    SEQUENCE
    CHARACTER(LEN=20) :: name = ''   ! hamiltonian/spin/orbital/spin_orbit/velocity/...
    CHARACTER(LEN=32) :: comp = ''   ! requested components; '' = all of the operator rank
    INTEGER :: total = 1             ! 1 = write the summed-over-atoms projection
  END TYPE t_wannierlib_op

  !> The <export> block: copy out matrices the run already holds, for someone outside.
  !>
  !> PROVISIONAL. Nothing inside FLEUR reads any of these back. They exist so a run can be
  !> cross-checked against a standalone wannier90.x, and so tools that consume Wannier data
  !> -- WannierBerri, irrep -- can be fed without recomputing. Grouped rather than left as
  !> four loose flags precisely because it is provisional: one object is one seam, whether
  !> the block is later promoted or dropped.
  TYPE t_wannierlib_export
    SEQUENCE
    LOGICAL :: w90 = .FALSE.     ! <export wannier90="T"/>:      .amn/.mmn/.eig
    LOGICAL :: basis = .FALSE.   ! <export wannierberri="T"/>:   WF<n>_basis.hdf
    LOGICAL :: gauge = .FALSE.   ! <export gauge="T"/>:          WF<n>_gauge.hdf, u_opt and u_mlwf
    LOGICAL :: bloch = .FALSE.   ! <export blochOperators="T"/>: WF<n>_s0.dat
  END TYPE t_wannierlib_export

  TYPE, EXTENDS(t_fleurinput_base) :: t_wannierlib_wannierize
    LOGICAL :: l_wannierize = .FALSE.
    ! Convenience flags DERIVED from the operator list ops(:) below (set in read_xml).
    ! They gate the coarse-mesh provider calls; the actual dispatch loops over ops(:).
    LOGICAL :: l_interpolation = .FALSE.   ! an <operator name="hamiltonian"> is requested
    LOGICAL :: l_spin = .FALSE.            ! an <operator name="spin"> is requested
    LOGICAL :: l_orbmom = .FALSE.          ! an <operator name="orbital"> is requested
    LOGICAL :: l_socop = .FALSE.           ! an <operator name="spin_orbit"> is requested
    LOGICAL :: l_operators_r = .FALSE.     ! an <operators_r> block (real-space O(R) export) is present
    !> The <export> block: one boolean per artefact, all off by default. They are grouped
    !> under one element because they answer the same question -- hand the run's own
    !> matrices to someone outside -- and differ only in who reads them. What <operators_r>
    !> and <interpolation> write is NOT here: those compute something first and writing it
    !> is the last step, while these copy out what the run already holds.
    TYPE(t_wannierlib_export) :: export
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
    !> Always allocated, possibly to size zero, so SIZE() is the count and the two cannot
    !> come apart. One entry per <operator> of <interpolation>.
    TYPE(t_wannierlib_op), ALLOCATABLE :: ops(:)

    ! --- operators_r table: one entry per <operators_r>/<operator name=".."> child ---
    !> One name per <operator> of <operators_r>. A list of names is all this is, so it
    !> stays a list of names rather than a record with one field in it.
    CHARACTER(LEN=20), ALLOCATABLE :: op_r_name(:)

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
    !> Keep the same number of states per spin channel at every k. Only meaningful when
    !> spin is a good quantum number; the routine checks and refuses otherwise.
    LOGICAL :: dis_spin_balanced = .FALSE.
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
            ! BUGFIX: when the user requested m=0 ("expand all m"), the auto-generated
            ! harmonic index must be the 1-based mval -- even when only one m exists
            ! (l=0/s: mrepeat=1). The old test (mrepeat>1) left s with the literal
            ! m=0, which is out of range for tlm(:,:,1:7) -> uninitialised read.
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
    USE m_mpi_bc_tool
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
      CALL mpi_bc(rank, mpi_comm, this%ops(iop)%comp)
      CALL mpi_bc(this%ops(iop)%total, rank, mpi_comm)
    END DO
    CALL mpi_bc(this%l_operators_r, rank, mpi_comm)
    CALL mpi_bc(this%export%w90, rank, mpi_comm)
    CALL mpi_bc(this%export%basis, rank, mpi_comm)
    CALL mpi_bc(this%export%gauge, rank, mpi_comm)
    CALL mpi_bc(this%export%bloch, rank, mpi_comm)
    CALL mpi_bc(this%l_plot_wf, rank, mpi_comm)
    CALL mpi_bc(this%op_r_name, rank, mpi_comm)
    CALL mpi_bc(this%num_wann, rank, mpi_comm)
    CALL mpi_bc(this%num_bands, rank, mpi_comm)
    CALL mpi_bc(this%min_band, rank, mpi_comm)
    CALL mpi_bc(this%max_band, rank, mpi_comm)
    CALL mpi_bc(this%dis_spin_balanced, rank, mpi_comm)
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
    USE m_types_xml
    USE m_constants
    USE m_types_kpts
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
    INTEGER :: nOpsRead = 0, nOprRead = 0
    LOGICAL :: has_wannierize
    CHARACTER(LEN=20) :: species_name
    TYPE(t_kpts) :: kpts_temp
    CHARACTER(len=40) :: kptsPath


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
      this%dis_spin_balanced = evaluateFirstBoolOnly( &
        xml%getAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spinBalanced'))
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
      nOpsRead = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF

    ALLOCATE(this%ops(nOpsRead))
    DO iProj = 1, nOpsRead
      WRITE(xPathP, '(A,I0,A)') '/fleurInput/output/wannierlib/interpolation/operator[', iProj, ']'
      this%ops(iProj)%name = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@name'))
      this%ops(iProj)%comp = ''
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@comp') == 1) &
         this%ops(iProj)%comp = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@comp'))
      this%ops(iProj)%total = 1
      IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathP))//'/@total') == 1) THEN
        IF (.NOT. evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@total'))) this%ops(iProj)%total = 0
      END IF

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
        CASE ('spin_orbit'); this%l_socop = .TRUE.
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
      this%export%basis = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@wannierberri'))
      this%export%gauge = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@gauge'))
      this%export%bloch = evaluateFirstBoolOnly(xml%getAttributeValue( &
        TRIM(ADJUSTL(xPathA))//'/@blochOperators'))
    END IF

    ! --- operators_r: real-space operator matrices O(R) (Fourier step 3, no interpolation).
    !     Groups <operator name=".."> children written to standalone-format files. ---
    xPathA = '/fleurInput/output/wannierlib/operators_r'
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))) == 1) THEN
      this%l_operators_r = .TRUE.
      nOprRead = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/operator')
    END IF
    ALLOCATE(this%op_r_name(nOprRead))
    DO iProj = 1, nOprRead
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

          spin_label = ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathP))//'/@spin'))
          SELECT CASE (spin_label(1:1))
          CASE ('u', 'U')
            this%proj(ip)%spin = 1
          CASE ('d', 'D')
            this%proj(ip)%spin = -1
          CASE DEFAULT
            this%proj(ip)%spin = 0
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
