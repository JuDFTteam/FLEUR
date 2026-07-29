!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Abstract base type of all Wannier matrix-element operators O.
!>
!>  Every operator O (spin, orbital moment, SOC, ...) that is transported into the
!>  Wannier representation goes through the same pipeline; this module owns all the
!>  operator-INDEPENDENT stages, a concrete operator only implements the deferred
!>  Bloch-basis builder calc_bloch (one k at a time):
!>
!>    (1) calc_coarse   : k-loop over THIS RANK's coarse-mesh k-points (distk slice):
!>                        read the eigenstates (read_eig + band-window selection via
!>                        m_wannierlib_get_z), compute the abc coefficients
!>                        (t_abc%calc_abc, both spins in spinor mode), and hand both
!>                        to the DEFERRED calc_bloch -> O^0(k) local slice o0_loc.
!>    (2) to_wannier    : gauge rotation O_W = V^dagger O^0 V with V = u_opt.u_matrix
!>                        on the local k-slice, then the distributed FT-reduce
!>                        (m_melem_ft) k -> R  ->  o_r(nw,nw,nrpts,ncomp).
!>                        [COLLECTIVE over mpi_comm]
!>    (3) interpolate   : R -> arbitrary k' list (Wannier interpolation of O).
!>    (4) write_bloch / write_realspace : plain-text IO of O^0(k) resp. O(R).
!>                        !! The file FORMATS ARE PRELIMINARY PLACEHOLDERS and will
!>                        !! be specified later; read_* are stubs until then.
!>
!>  Shared FLEUR objects reach the deferred builder through the context bundle
!>  t_mel_context (POINTER components, filled once by the driver).
!>
!>  Planned extensions (NOT part of this module yet):
!>   * position/Berry connection: needs the neighbour overlaps M(k,b) + b-mesh
!>     weights instead of a single-k Bloch matrix -> an intermediate abstract type
!>     (e.g. t_matrixelement_overlap) overriding calc_coarse with an mmn-style loop.
!>   * hamiltonian: bypasses the Bloch pass (built from eig + U on the full mesh);
!>     a derived type overriding calc_coarse (no-op) and to_wannier.
!>   * orbital / soc derived types: wrappers of melem_orbmom_bloch(_collinear)
!>     and melem_socmat_bloch (the latter mutates ctx%usdus and consumes
!>     ctx%enpara/vtot/fmpi, which the context already carries).
!>   * band-projected output (zheev on H(k'), <O>_n = [C^dagger O C]_nn, the
!>     bands_wann_*.dat files): composed by the future driver from the interpolate()
!>     results of a Hamiltonian instance and the operator instances.
!>   * write_bloch currently reassembles the full coarse mesh via an ALLREDUCE
!>     (O(nb^2*ncomp*nkptf) buffer); switch to a gatherv/streaming writer once the
!>     plain-text format is frozen.
MODULE m_types_matrixelement
   USE m_juDFT
   USE m_constants, ONLY : oUnit
   USE m_types_atoms
   USE m_types_cell
   USE m_types_input
   USE m_types_sym
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_stars
   USE m_types_kpts
   USE m_types_enpara
   USE m_types_potden
   USE m_types_mpi
   USE m_types_usdus
   USE m_types_radfun
   USE m_types_lapw
   USE m_types_mat
   USE m_types_abc
   USE m_types_wannierlib
   IMPLICIT NONE
   PRIVATE

   !> Context bundle: POINTERs to the shared FLEUR objects the per-k Bloch builders
   !> need. Built ONCE by the driver via %init.
   !> IMPORTANT: all actuals passed to %init MUST have the TARGET attribute in the
   !> driver -- the pointers stay associated for the lifetime of the context.
   TYPE :: t_mel_context
      TYPE(t_atoms),    POINTER :: atoms => NULL()
      TYPE(t_cell),     POINTER :: cell => NULL()
      TYPE(t_input),    POINTER :: input => NULL()
      TYPE(t_sym),      POINTER :: sym => NULL()
      TYPE(t_noco),     POINTER :: noco => NULL()
      TYPE(t_nococonv), POINTER :: nococonv => NULL()
      TYPE(t_stars),    POINTER :: stars => NULL()
      TYPE(t_kpts),     POINTER :: kpts => NULL()
      TYPE(t_enpara),   POINTER :: enpara => NULL()
      TYPE(t_potden),   POINTER :: vtot => NULL()
      TYPE(t_mpi),      POINTER :: fmpi => NULL()
      TYPE(t_usdus),    POINTER :: usdus => NULL()     ! mutable target: the SOC builder fills its spnorb parts
      TYPE(t_radfun),   POINTER :: radfun(:) => NULL() ! (ntype), radials already generated
      TYPE(t_wannierlib_wannierize), POINTER :: wann => NULL()
      INTEGER :: eig_id = -1
      ! run-mode switches derived once in init (same rules as wannierlib_main):
      LOGICAL :: l_spinors = .FALSE.   ! noco%l_noco .OR. noco%l_soc -> full 2N spinor in record 1
      LOGICAL :: l_real_wann = .FALSE. ! input%l_real .AND. .NOT. noco%l_soc (a complex spinor read
                                       ! into a real buffer would lose the imaginary part)
      INTEGER :: jspin = 1             ! collinear jspins=2 channel (1/2); 1 in spinor mode
   CONTAINS
      PROCEDURE :: init => mel_context_init
   END TYPE t_mel_context

   TYPE, ABSTRACT :: t_matrixelement
      ! -- identification (fixed by the derived type's init BEFORE init_melem) --
      CHARACTER(LEN=20) :: name = 'unnamed' ! operator id; default prefix of the IO files
      INTEGER :: ncomp = 0                  ! # components alpha (spin: 3, soc: 1 or 4, per-atom: 3*nat flattened)
      INTEGER :: nsites = 1                 ! informational: >1 -> ncomp is a (ncomp/nsites) x nsites flattening
      ! -- dimensions (copied from t_wannierlib_wannierize in init_melem) --
      INTEGER :: nb = 0                     ! num_bands = size of the selected band window
      INTEGER :: nw = 0                     ! num_wann  (Wannier dimension, from to_wannier on)
      INTEGER :: min_band = -1, max_band = -1
      ! -- MPI / coarse-mesh distribution bookkeeping --
      INTEGER :: mpi_comm = 0, irank = 0, isize = 1
      INTEGER :: nk_loc = 0                 ! # coarse k-points owned by this rank
      INTEGER, ALLOCATABLE :: distk(:)      ! (nkptf) owner rank of each global coarse k
      INTEGER, ALLOCATABLE :: gk_loc(:)     ! (nk_loc) global k of this rank's slice, ascending
      ! -- Bloch-basis O^0: per-rank slice only, the full mesh is NEVER materialized --
      COMPLEX, ALLOCATABLE :: o0_loc(:, :, :, :)  ! (nb, nb, ncomp, MAX(1,nk_loc))
      ! -- real-space Wannier O(R): replicated on all ranks after to_wannier --
      COMPLEX, ALLOCATABLE :: o_r(:, :, :, :)     ! (nw, nw, nrpts, ncomp)
      INTEGER, ALLOCATABLE :: irvec(:, :)         ! (3, nrpts) Wigner-Seitz R-vectors (W90 convention)
      INTEGER, ALLOCATABLE :: ndegen(:)           ! (nrpts)
      INTEGER :: nrpts = 0
   CONTAINS
      PROCEDURE(calc_bloch_iface), DEFERRED :: calc_bloch ! the ONE thing an operator implements
      ! NOTE: the base initializer is deliberately NOT named "init" -- Fortran has no
      ! super-call, so each derived type declares its own init (any argument list)
      ! and chains CALL this%init_melem(ctx, name, ncomp [, nsites]). Keep it that way.
      PROCEDURE :: init_melem
      PROCEDURE :: calc_coarse     ! template method: k-loop -> read_eig -> abc -> calc_bloch
      PROCEDURE :: to_wannier      ! gauge rotation + FT-reduce k -> R          [COLLECTIVE]
      PROCEDURE :: interpolate     ! R -> arbitrary k' list                     [local]
      PROCEDURE :: write_bloch     ! plain-text O^0(k), PRELIMINARY format      [COLLECTIVE]
      PROCEDURE :: write_realspace ! plain-text O(R),   PRELIMINARY format      [rank 0]
      PROCEDURE :: read_bloch      ! stub until the format is frozen
      PROCEDURE :: read_realspace  ! stub until the format is frozen
      PROCEDURE :: free
   END TYPE t_matrixelement

   !> Polymorphic wrapper so several operators can share ONE coarse k-pass
   !> (one read_eig + one abc set per k), see matrixelement_calc_coarse_all.
   TYPE :: t_matrixelement_list
      CLASS(t_matrixelement), POINTER :: p => NULL()
   END TYPE t_matrixelement_list

   ABSTRACT INTERFACE
      !> Build this operator's Bloch-basis matrices at ONE coarse k (global index ik).
      !> lapw/zMat are the per-k workspace from wannierlib_get_z (band window already
      !> selected, zMat holds the full spinor in spinor mode). abc(itype, isp) are the
      !> matching coefficients; convention: in spinor mode (ctx%l_spinors) BOTH spin
      !> slots isp=1,2 are filled (calc_abc jspin=1,2 of the one spinor zMat); in the
      !> collinear jspins=2 case only slot ctx%jspin is filled. o0_k(nb,nb,ncomp)
      !> aliases this operator's o0_loc(:,:,:,ik_local) and arrives pre-zeroed.
      SUBROUTINE calc_bloch_iface(this, ctx, ik, lapw, zMat, abc, o0_k)
         IMPORT t_matrixelement, t_mel_context, t_lapw, t_mat, t_abc
         CLASS(t_matrixelement), INTENT(INOUT) :: this
         TYPE(t_mel_context),    INTENT(IN)    :: ctx  ! usdus target mutable through the pointer (SOC)
         INTEGER,                INTENT(IN)    :: ik
         TYPE(t_lapw),           INTENT(IN)    :: lapw
         TYPE(t_mat),            INTENT(IN)    :: zMat
         TYPE(t_abc),            INTENT(IN)    :: abc(:, :)     ! (ntype, 2)
         COMPLEX,                INTENT(INOUT) :: o0_k(:, :, :) ! (nb, nb, ncomp)
      END SUBROUTINE calc_bloch_iface
   END INTERFACE

   PUBLIC :: t_matrixelement, t_matrixelement_list, t_mel_context
   PUBLIC :: matrixelement_calc_coarse_all

CONTAINS

   !> Associate the context pointers and derive the run-mode switches.
   !> ALL actuals must be TARGETs in the caller (see type documentation).
   SUBROUTINE mel_context_init(this, atoms, cell, input, sym, noco, nococonv, stars, kpts, &
                               enpara, vtot, fmpi, usdus, radfun, wann, eig_id, jspin)
      CLASS(t_mel_context), INTENT(OUT) :: this
      TYPE(t_atoms),    TARGET, INTENT(IN)    :: atoms
      TYPE(t_cell),     TARGET, INTENT(IN)    :: cell
      TYPE(t_input),    TARGET, INTENT(IN)    :: input
      TYPE(t_sym),      TARGET, INTENT(IN)    :: sym
      TYPE(t_noco),     TARGET, INTENT(IN)    :: noco
      TYPE(t_nococonv), TARGET, INTENT(IN)    :: nococonv
      TYPE(t_stars),    TARGET, INTENT(IN)    :: stars
      TYPE(t_kpts),     TARGET, INTENT(IN)    :: kpts
      TYPE(t_enpara),   TARGET, INTENT(IN)    :: enpara
      TYPE(t_potden),   TARGET, INTENT(IN)    :: vtot
      TYPE(t_mpi),      TARGET, INTENT(IN)    :: fmpi
      TYPE(t_usdus),    TARGET, INTENT(INOUT) :: usdus     ! SOC builder fills spnorb parts
      TYPE(t_radfun),   TARGET, INTENT(IN)    :: radfun(:) ! (ntype), radials already generated
      TYPE(t_wannierlib_wannierize), TARGET, INTENT(IN) :: wann
      INTEGER, INTENT(IN)           :: eig_id
      INTEGER, INTENT(IN), OPTIONAL :: jspin   ! collinear jspins=2 channel; default 1

      this%atoms => atoms; this%cell => cell; this%input => input; this%sym => sym
      this%noco => noco; this%nococonv => nococonv; this%stars => stars; this%kpts => kpts
      this%enpara => enpara; this%vtot => vtot; this%fmpi => fmpi; this%usdus => usdus
      this%radfun => radfun; this%wann => wann
      this%eig_id = eig_id

      this%l_spinors   = noco%l_noco .OR. noco%l_soc
      this%l_real_wann = input%l_real .AND. .NOT. noco%l_soc
      this%jspin = 1
      IF (PRESENT(jspin)) this%jspin = jspin
      IF (this%jspin < 1 .OR. this%jspin > 2) &
         CALL juDFT_error('t_mel_context: jspin must be 1 or 2', calledby='mel_context_init')
      IF (this%l_spinors .AND. this%jspin /= 1) &
         CALL juDFT_error('t_mel_context: jspin has no meaning in spinor mode (use 1)', &
                          calledby='mel_context_init')
   END SUBROUTINE mel_context_init

   !> Base initializer, chained from the derived type's own init AFTER it fixed
   !> name/ncomp/nsites: copies the dimensions from ctx%wann, sets up the MPI
   !> k-distribution (distk owner map + this rank's gk_loc slice) and allocates
   !> the local Bloch slice o0_loc. All ranks, no communication (but ctx must be
   !> consistent across ranks so every rank builds the same distk).
   SUBROUTINE init_melem(this, ctx, name, ncomp, nsites)
      CLASS(t_matrixelement), INTENT(INOUT) :: this
      TYPE(t_mel_context),    INTENT(IN)    :: ctx
      CHARACTER(LEN=*),       INTENT(IN)    :: name
      INTEGER,                INTENT(IN)    :: ncomp
      INTEGER, INTENT(IN), OPTIONAL         :: nsites

      INTEGER :: ik, nk_blk, ierr

      IF (ncomp < 1) CALL juDFT_error('t_matrixelement: ncomp must be >= 1', calledby='init_melem')
      this%name = name
      this%ncomp = ncomp
      this%nsites = 1
      IF (PRESENT(nsites)) this%nsites = nsites

      IF (ctx%wann%min_band < 1 .OR. ctx%wann%max_band < ctx%wann%min_band) &
         CALL juDFT_error('t_matrixelement: invalid band window in wannierlib input', calledby='init_melem')
      this%nb = ctx%wann%num_bands
      this%nw = ctx%wann%num_wann
      this%min_band = ctx%wann%min_band
      this%max_band = ctx%wann%max_band

      this%mpi_comm = ctx%fmpi%mpi_comm
      this%irank = ctx%fmpi%irank
      this%isize = ctx%fmpi%isize

      ! distk: which rank owns each global coarse k (same rule as wannierlib_main)
      IF (ALLOCATED(this%distk)) DEALLOCATE(this%distk)
      ALLOCATE(this%distk(ctx%kpts%nkptf), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error('t_matrixelement: failed allocating distk', calledby='init_melem')
      IF (ALLOCATED(ctx%fmpi%coulomb_owner) .AND. SIZE(ctx%fmpi%coulomb_owner) == ctx%kpts%nkptf) THEN
         this%distk = ctx%fmpi%coulomb_owner
      ELSE
         ! fallback to a contiguous block distribution if no global owner map is available
         nk_blk = ctx%kpts%nkptf/MAX(1, this%isize)
         IF (MOD(ctx%kpts%nkptf, MAX(1, this%isize)) > 0) nk_blk = nk_blk + 1
         DO ik = 1, ctx%kpts%nkptf
            this%distk(ik) = MIN((ik - 1)/MAX(1, nk_blk), MAX(0, this%isize - 1))
         END DO
      END IF

      ! gk_loc: this rank's global-k indices, ascending (= storage order of o0_loc)
      this%nk_loc = COUNT(this%distk == this%irank)
      IF (ALLOCATED(this%gk_loc)) DEALLOCATE(this%gk_loc)
      ALLOCATE(this%gk_loc(this%nk_loc))
      this%nk_loc = 0
      DO ik = 1, ctx%kpts%nkptf
         IF (this%distk(ik) /= this%irank) CYCLE
         this%nk_loc = this%nk_loc + 1
         this%gk_loc(this%nk_loc) = ik
      END DO

      IF (ALLOCATED(this%o0_loc)) DEALLOCATE(this%o0_loc)
      ALLOCATE(this%o0_loc(this%nb, this%nb, this%ncomp, MAX(1, this%nk_loc)), &
               stat=ierr, source=CMPLX(0.0, 0.0))
      IF (ierr /= 0) CALL juDFT_error('t_matrixelement: failed allocating o0_loc', calledby='init_melem')
   END SUBROUTINE init_melem

   !> Template method for a single operator: runs the shared coarse k-loop with
   !> just this operator. To share one k-pass (one read_eig + one abc set per k)
   !> between several operators, call matrixelement_calc_coarse_all directly.
   SUBROUTINE calc_coarse(this, ctx)
      CLASS(t_matrixelement), TARGET, INTENT(INOUT) :: this
      TYPE(t_mel_context),            INTENT(IN)    :: ctx

      TYPE(t_matrixelement_list) :: ops(1)

      ops(1)%p => this
      CALL matrixelement_calc_coarse_all(ops, ctx)
   END SUBROUTINE calc_coarse

   !> The shared coarse-mesh k-loop (template body): for each k owned by this rank,
   !> read the eigenstates (band window selected), build the abc coefficients, and
   !> invoke every operator's deferred calc_bloch on the same workspace.
   !> All ranks, no MPI calls (each rank computes only its distk slice).
   SUBROUTINE matrixelement_calc_coarse_all(ops, ctx)
      USE m_wannierlib_get_z
      TYPE(t_matrixelement_list), INTENT(IN) :: ops(:)
      TYPE(t_mel_context),        INTENT(IN) :: ctx

      TYPE(t_lapw) :: lapw
      TYPE(t_mat)  :: zMat
      TYPE(t_abc), ALLOCATABLE :: abc(:, :)
      INTEGER :: iop, ik, il, isp, isp1, isp2, itype, jspin_rad

      IF (SIZE(ops) < 1) RETURN
      DO iop = 1, SIZE(ops)
         IF (.NOT. ASSOCIATED(ops(iop)%p)) &
            CALL juDFT_error('matrixelement_calc_coarse_all: unassociated operator', &
                             calledby='matrixelement_calc_coarse_all')
         IF (.NOT. ALLOCATED(ops(iop)%p%o0_loc)) &
            CALL juDFT_error('matrixelement_calc_coarse_all: operator "'//TRIM(ops(iop)%p%name)// &
                             '" not initialized (call init first)', calledby='matrixelement_calc_coarse_all')
         IF (ops(iop)%p%nb /= ctx%wann%num_bands) &
            CALL juDFT_error('matrixelement_calc_coarse_all: band-window mismatch for operator "'// &
                             TRIM(ops(iop)%p%name)//'"', calledby='matrixelement_calc_coarse_all')
         IF (ANY(ops(iop)%p%distk /= ops(1)%p%distk)) &
            CALL juDFT_error('matrixelement_calc_coarse_all: operators use different k-distributions', &
                             calledby='matrixelement_calc_coarse_all')
      END DO

      CALL timestart('melem coarse k-loop')

      ! spin slots to fill (see calc_bloch_iface): spinor mode -> both; collinear -> ctx%jspin only
      IF (ctx%l_spinors) THEN
         isp1 = 1; isp2 = 2
      ELSE
         isp1 = ctx%jspin; isp2 = ctx%jspin
      END IF

      ALLOCATE(abc(ctx%atoms%ntype, 2))
      il = 0
      DO ik = 1, ctx%kpts%nkptf
         IF (ops(1)%p%distk(ik) /= ops(1)%p%irank) CYCLE   ! this rank computes only its own k-slice
         il = il + 1

         ! spinor mode: the full 2N spinor sits in record 1; collinear: record ctx%jspin
         CALL wannierlib_get_z(ctx%wann, ctx%eig_id, ctx%input, ctx%atoms, ctx%noco, ctx%nococonv, &
                               ctx%kpts, ctx%sym, ctx%cell, ik, MERGE(1, ctx%jspin, ctx%l_spinors), &
                               ctx%l_real_wann, lapw, zMat)

         DO isp = isp1, isp2
            ! radial spin index: with jspins=1 only one set of radials exists -> use 1
            jspin_rad = MERGE(1, isp, ctx%input%jspins == 1)
            DO itype = 1, ctx%atoms%ntype
               CALL abc(itype, isp)%init(ctx%input, ctx%atoms, ctx%wann%num_bands, itype)
               CALL abc(itype, isp)%calc_abc(ctx%input, ctx%atoms, ctx%sym, ctx%cell, lapw, &
                                             ctx%wann%num_bands, ctx%usdus, ctx%noco, ctx%nococonv, &
                                             jspin_rad, itype, zMat)
            END DO
         END DO

         DO iop = 1, SIZE(ops)
            ops(iop)%p%o0_loc(:, :, :, il) = CMPLX(0.0, 0.0)
            CALL ops(iop)%p%calc_bloch(ctx, ik, lapw, zMat, abc, ops(iop)%p%o0_loc(:, :, :, il))
         END DO
      END DO
      DEALLOCATE(abc)

      CALL timestop('melem coarse k-loop')
   END SUBROUTINE matrixelement_calc_coarse_all

   !> Rotate the local Bloch slice into the Wannier gauge, O_W = V^dagger O^0 V with
   !> V(k) = u_opt(k).u_matrix(k), and Fourier transform to real space via the
   !> distributed FT-reduce -> o_r/irvec/ndegen/nrpts, replicated on all ranks.
   !> COLLECTIVE over mpi_comm (contains an MPI_ALLREDUCE); u_matrix/u_opt must be
   !> the full-mesh matrices, replicated on every rank.
   SUBROUTINE to_wannier(this, ctx, u_matrix, u_opt)
      USE m_melem_ft, ONLY : melem_ft_to_real_reduce
      CLASS(t_matrixelement), INTENT(INOUT) :: this
      TYPE(t_mel_context),    INTENT(IN)    :: ctx
      COMPLEX, INTENT(IN) :: u_matrix(:, :, :)  ! (nw, nw, nkptf) MLWF gauge, full mesh
      COMPLEX, INTENT(IN) :: u_opt(:, :, :)     ! (nb, nw, nkptf) disentanglement, full mesh

      COMPLEX, ALLOCATABLE :: vloc(:, :, :), ow_loc(:, :, :, :), tmp(:, :), o1(:, :, :)
      INTEGER :: kl, a, gk

      IF (.NOT. ALLOCATED(this%o0_loc)) &
         CALL juDFT_error('t_matrixelement%to_wannier: no Bloch data (call calc_coarse first)', &
                          calledby='to_wannier')
      CALL timestart('melem to_wannier')

      ALLOCATE(vloc(this%nb, this%nw, MAX(1, this%nk_loc)), source=CMPLX(0.0, 0.0))
      DO kl = 1, this%nk_loc
         gk = this%gk_loc(kl)
         vloc(:, :, kl) = MATMUL(u_opt(:, :, gk), u_matrix(:, :, gk))
      END DO

      ALLOCATE(ow_loc(this%nw, this%nw, this%ncomp, MAX(1, this%nk_loc)), source=CMPLX(0.0, 0.0))
      ALLOCATE(tmp(this%nb, this%nw))
      DO kl = 1, this%nk_loc
         DO a = 1, this%ncomp
            tmp = MATMUL(this%o0_loc(:, :, a, kl), vloc(:, :, kl))
            ow_loc(:, :, a, kl) = MATMUL(CONJG(TRANSPOSE(vloc(:, :, kl))), tmp)
         END DO
      END DO
      DEALLOCATE(tmp, vloc)

      IF (ALLOCATED(this%o_r)) DEALLOCATE(this%o_r)
      DO a = 1, this%ncomp
         CALL melem_ft_to_real_reduce(ctx%cell, ctx%kpts, ow_loc(:, :, a, :), this%gk_loc, &
                                           this%mpi_comm, o1, this%irvec, this%ndegen, this%nrpts)
         IF (a == 1) ALLOCATE(this%o_r(this%nw, this%nw, this%nrpts, this%ncomp))
         this%o_r(:, :, :, a) = o1
         DEALLOCATE(o1)
      END DO
      DEALLOCATE(ow_loc)

      CALL timestop('melem to_wannier')
   END SUBROUTINE to_wannier

   !> Wannier-interpolate the matrix element onto an arbitrary fractional k-list:
   !> O(k') = sum_R e^{+i 2pi k'.R} / ndegen(R) * O(R), per component. Local (o_r is
   !> replicated after to_wannier), usable on any rank.
   SUBROUTINE interpolate(this, kfrac, o_interp)
      USE m_melem_ft, ONLY : melem_ft_rtok
      CLASS(t_matrixelement), INTENT(IN)  :: this
      REAL,                   INTENT(IN)  :: kfrac(:, :)           ! (3, np) fractional k'
      COMPLEX, ALLOCATABLE,   INTENT(OUT) :: o_interp(:, :, :, :)  ! (nw, nw, ncomp, np)

      COMPLEX, ALLOCATABLE :: o_one(:, :, :)
      INTEGER :: a

      IF (.NOT. ALLOCATED(this%o_r)) &
         CALL juDFT_error('t_matrixelement%interpolate: no O(R) data (call to_wannier first)', &
                          calledby='interpolate')
      CALL timestart('melem interpolate')
      ALLOCATE(o_interp(this%nw, this%nw, this%ncomp, SIZE(kfrac, 2)))
      DO a = 1, this%ncomp
         CALL melem_ft_rtok(this%o_r(:, :, :, a), this%irvec, this%ndegen, this%nrpts, &
                                 kfrac, o_one)
         o_interp(:, :, a, :) = o_one
         DEALLOCATE(o_one)
      END DO
      CALL timestop('melem interpolate')
   END SUBROUTINE interpolate

   !> Write the Bloch-basis O^0(k) on the full coarse mesh to a plain-text file.
   !> COLLECTIVE: the full mesh is reassembled from the per-rank slices with a
   !> zero-padded ALLREDUCE (the wannierlib_reduce_amn pattern), rank 0 writes.
   !> !! FORMAT PLACEHOLDER -- the final format will be specified later. !!
   SUBROUTINE write_bloch(this, ctx, fname)
#ifdef CPP_MPI
      USE mpi
#endif
      CLASS(t_matrixelement), INTENT(IN) :: this
      TYPE(t_mel_context),    INTENT(IN) :: ctx
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fname   ! default TRIM(name)//'_o0k.dat'

      COMPLEX, ALLOCATABLE :: o0_full(:, :, :, :)
      CHARACTER(LEN=64) :: fn
      INTEGER :: kl, ik, a, i, j, iu, ierr

      IF (.NOT. ALLOCATED(this%o0_loc)) &
         CALL juDFT_error('t_matrixelement%write_bloch: no Bloch data (call calc_coarse first)', &
                          calledby='write_bloch')

      ! NOTE: O(nb^2*ncomp*nkptf) buffer -- to be replaced by a gatherv/streaming
      ! writer once the plain-text format is frozen.
      ALLOCATE(o0_full(this%nb, this%nb, this%ncomp, ctx%kpts%nkptf), source=CMPLX(0.0, 0.0))
      DO kl = 1, this%nk_loc
         o0_full(:, :, :, this%gk_loc(kl)) = this%o0_loc(:, :, :, kl)
      END DO
#ifdef CPP_MPI
      IF (this%isize > 1) &
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, o0_full, SIZE(o0_full), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                            this%mpi_comm, ierr)
#endif

      IF (this%irank == 0) THEN
         fn = TRIM(this%name)//'_o0k.dat'
         IF (PRESENT(fname)) fn = fname
         OPEN(newunit=iu, file=TRIM(fn), status='replace')
         WRITE(iu, '(a)') '# FLEUR wannierlib matrix element, Bloch basis O^0(k) -- PRELIMINARY FORMAT'
         WRITE(iu, '(a,1x,5i8)') TRIM(this%name), this%nb, this%ncomp, ctx%kpts%nkptf, &
                                 this%min_band, this%max_band
         DO ik = 1, ctx%kpts%nkptf
            DO a = 1, this%ncomp
               WRITE(iu, '(2i8)') ik, a
               DO j = 1, this%nb
                  DO i = 1, this%nb
                     WRITE(iu, '(2i5,2f24.18)') i, j, REAL(o0_full(i, j, a, ik)), AIMAG(o0_full(i, j, a, ik))
                  END DO
               END DO
            END DO
         END DO
         CLOSE(iu)
         WRITE(oUnit, '(a)') 'wannierlib melem: wrote '//TRIM(fn)//' (Bloch O^0(k), preliminary format)'
      END IF
      DEALLOCATE(o0_full)
   END SUBROUTINE write_bloch

   !> Write the real-space Wannier representation O(R) to a plain-text file.
   !> Rank 0 only (early return elsewhere -> safe to call from all ranks).
   !> !! FORMAT PLACEHOLDER -- the final format will be specified later. !!
   !> Operator-specific legacy formats (WFn_hr.dat, rspauli.1 spinor-block layout, ...)
   !> are obtained by derived types OVERRIDING this binding.
   SUBROUTINE write_realspace(this, fname)
      CLASS(t_matrixelement), INTENT(IN) :: this
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fname   ! default TRIM(name)//'_or.dat'

      CHARACTER(LEN=64) :: fn
      INTEGER :: irpt, a, i, j, iu, c

      IF (this%irank /= 0) RETURN
      IF (.NOT. ALLOCATED(this%o_r)) &
         CALL juDFT_error('t_matrixelement%write_realspace: no O(R) data (call to_wannier first)', &
                          calledby='write_realspace')

      fn = TRIM(this%name)//'_or.dat'
      IF (PRESENT(fname)) fn = fname
      OPEN(newunit=iu, file=TRIM(fn), status='replace')
      WRITE(iu, '(a)') '# FLEUR wannierlib matrix element, real space O(R) -- PRELIMINARY FORMAT'
      WRITE(iu, '(a,1x,4i8)') TRIM(this%name), this%nw, this%nrpts, this%ncomp, this%nsites
      c = 0
      DO irpt = 1, this%nrpts
         WRITE(iu, '(i5)', advance='no') this%ndegen(irpt); c = c + 1
         IF (MOD(c, 15) == 0) WRITE(iu, '(a)') ''
      END DO
      IF (MOD(c, 15) /= 0) WRITE(iu, '(a)') ''
      DO irpt = 1, this%nrpts
         DO a = 1, this%ncomp
            DO j = 1, this%nw
               DO i = 1, this%nw
                  WRITE(iu, '(6i5,2f20.8)') this%irvec(:, irpt), i, j, a, &
                     REAL(this%o_r(i, j, irpt, a)), AIMAG(this%o_r(i, j, irpt, a))
               END DO
            END DO
         END DO
      END DO
      CLOSE(iu)
      WRITE(oUnit, '(a,i0,a)') 'wannierlib melem: wrote '//TRIM(fn)//' (O(R), ', this%nrpts, &
                               ' R-vectors, preliminary format)'
   END SUBROUTINE write_realspace

   !> Stub: reading is defined together with the final file format.
   SUBROUTINE read_bloch(this, ctx, fname)
      CLASS(t_matrixelement), INTENT(INOUT) :: this
      TYPE(t_mel_context),    INTENT(IN)    :: ctx
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fname

      CALL juDFT_error('t_matrixelement%read_bloch: IO format not yet frozen', calledby='read_bloch')
   END SUBROUTINE read_bloch

   !> Stub: reading is defined together with the final file format.
   SUBROUTINE read_realspace(this, fname)
      CLASS(t_matrixelement), INTENT(INOUT) :: this
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fname

      CALL juDFT_error('t_matrixelement%read_realspace: IO format not yet frozen', calledby='read_realspace')
   END SUBROUTINE read_realspace

   SUBROUTINE free(this)
      CLASS(t_matrixelement), INTENT(INOUT) :: this

      IF (ALLOCATED(this%distk)) DEALLOCATE(this%distk)
      IF (ALLOCATED(this%gk_loc)) DEALLOCATE(this%gk_loc)
      IF (ALLOCATED(this%o0_loc)) DEALLOCATE(this%o0_loc)
      IF (ALLOCATED(this%o_r)) DEALLOCATE(this%o_r)
      IF (ALLOCATED(this%irvec)) DEALLOCATE(this%irvec)
      IF (ALLOCATED(this%ndegen)) DEALLOCATE(this%ndegen)
      this%nrpts = 0; this%nk_loc = 0
      this%nb = 0; this%nw = 0; this%ncomp = 0; this%nsites = 1
      this%min_band = -1; this%max_band = -1
   END SUBROUTINE free

END MODULE m_types_matrixelement
