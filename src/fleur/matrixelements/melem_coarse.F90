!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The coarse-mesh pass: the Bloch-basis matrices O^0(k) of every requested operator, on
!>  the ab-initio eigenstates and BEFORE any Wannier gauge exists.
!>
!>    %init   allocate the per-rank slices and decide which of them the requested
!>            <operator>/<operators_r> lists actually need
!>    %calc   ONE k-pass that builds every requested operator through
!>            matrix_element_factory, which reads the states and their coefficients and
!>            keeps them, so one k is read once however many operators ask for it
!>
!>  Everything that needs the gauge lives in m_wgauge_run; the collinear combined spin
!>  operator, which needs both channels wannierised, in m_wgauge_spin_collinear.
MODULE m_melem_coarse
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_atoms
   USE m_types_cell
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_sym
   USE m_types_enpara
   USE m_types_mpi
   USE m_types_potden
   USE m_melem_symrot, ONLY: melem_rot_axial
   USE m_types_stars
   USE m_types_mat
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_window, ONLY: t_melem_window
   USE m_melem_spin, ONLY: melem_pauli_from_blocks, melem_spin_sumrule
   USE m_types_matelements_spin, ONLY: t_matelements_spin
   USE m_types_vacuum
   USE m_types_matelements_soc, ONLY: t_matelements_soc
   USE m_types_rsoc, ONLY: t_rsoc
   USE m_types_matelements_orbital, ONLY: t_matelements_orbital
   USE m_matrix_element_factory, ONLY: matrix_element_factory
   USE m_melem_check, ONLY: melem_check_provider, melem_check_matrix
   IMPLICIT NONE
   PRIVATE

   !> The Bloch-basis coarse-mesh operator matrices O^0(k), per-rank slices only: the full
   !> coarse mesh is never materialized. Slice entries are stored in ascending global-k order,
   !> matching the gk_loc convention of the distributed FT-reduce.
   TYPE :: t_melem_coarse
      COMPLEX, ALLOCATABLE :: s0(:, :, :, :)      !< (nb,nb,3,nk_loc)      spin
      !> Orbital L per atom and per channel. One channel when the states are spinors, two
      !> when they are separate eigenproblems -- the same quantity either way, so it is one
      !> array and not two.
      COMPLEX, ALLOCATABLE :: l0(:, :, :, :, :, :) !< (nb,nb,3,nat,channel,nk_loc)
      COMPLEX, ALLOCATABLE :: soc0(:, :, :, :)    !< (nb,nb,1,nk_loc)      SOC
      COMPLEX, ALLOCATABLE :: soc4(:, :, :, :)    !< (nb,nb,4,nk_loc)      2x2 SOC spinor blocks
      !> collinear jspins=2 only: the cross-spin overlap <up|dn> per k, in the Bloch basis.
      !> It is the ingredient of the combined 2N spin operator, and it needs no gauge, so it
      !> is built here with the other coarse matrices even though the operator itself cannot
      !> be assembled until both channels are wannierised.
      COMPLEX, ALLOCATABLE :: x0(:, :, :)         !< (nb,nb,nk_loc)
      !> How many independent eigenproblems the states come from: two when jspins=2 without
      !> spinors, one otherwise. It is a fact about the calculation, where the flags it
      !> replaces were a copy of what the request already said, kept next to it and able to
      !> disagree.
      INTEGER :: n_channels = 1
      !> .TRUE. only when the spinor coarse slices were really allocated (an operator is requested
      !> AND we have spinor wavefunctions). Gates %calc, so it can never write into the stubs.
      LOGICAL :: l_active = .FALSE.
   CONTAINS
      PROCEDURE :: init => melem_coarse_init
      PROCEDURE :: calc => melem_coarse_calc
      PROCEDURE :: free => melem_coarse_free
   END TYPE t_melem_coarse

   PUBLIC :: t_melem_coarse

CONTAINS

   SUBROUTINE melem_coarse_init(this, request, window, atoms, input, kpts, fmpi, distk, l_spinors)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_melem_request), INTENT(IN) :: request
      CLASS(t_melem_window), INTENT(IN) :: window
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)
      LOGICAL, INTENT(IN) :: l_spinors   !< noco%l_noco .OR. noco%l_soc

      INTEGER :: nkc_loc
      LOGICAL :: l_ch_orb, l_ch_spin, l_ch_soc

      ! Operator Bloch matrices on the coarse mesh: the k-loop is DISTRIBUTED over ranks (each
      ! its distk slice, so the reads the factory does for it are parallel too) into per-rank
      ! local arrays. Every consumer works on those slices plus a distributed FT-reduce, so the
      ! full mesh is never assembled.
      this%n_channels = MERGE(2, 1, input%jspins == 2 .AND. .NOT. l_spinors)
      this%l_active = (request%l_spin .OR. request%l_orbmom .OR. request%l_socop) .AND. l_spinors
      !> L is spin-diagonal, so each channel has its own and either list is reason enough
      !> to build it. The cross-spin overlap is only ever wanted by the real-space export:
      !> the interpolated spin operator does not exist here, see the guard below.
      l_ch_orb  = this%n_channels == 2 .AND. request%needs_op('orbital')
      l_ch_spin = this%n_channels == 2 .AND. request%has_op_r('spin')
      !> Two channels do NOT mean no spin-orbit coupling. Keeping it out of the eigenproblem
      !> and adding it afterwards as an operator in the Wannier basis is a construction of
      !> its own, and the one that keeps the Wannier functions spin-pure: with SOC in the
      !> first variation the states are spinors, nothing in the spread functional penalises
      !> a spin character that varies with k, and the exported S(R) and L(R) stop being
      !> localised. Like the spin, it needs both gauges, so it is a real-space export only.
      l_ch_soc  = this%n_channels == 2 .AND. request%has_op_r('spin_orbit')
      nkc_loc = MAX(1, COUNT(distk == fmpi%irank))

      IF (this%l_active) THEN
         ALLOCATE (this%s0(window%num_bands, window%num_bands, 3, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc0(window%num_bands, window%num_bands, 1, nkc_loc), source=cmplx(0.0, 0.0))
      END IF
      !> soc0 is the sum over the four blocks, which only means anything for a spinor; the
      !> blocks themselves are what both paths export, so they are allocated for both.
      IF (this%l_active .OR. l_ch_soc) THEN
         ALLOCATE (this%soc4(window%num_bands, window%num_bands, 4, nkc_loc), source=cmplx(0.0, 0.0))
      END IF
      !> Otherwise they are not allocated at all. A (1,1,1,1) stub is a VALID array of
      !> the wrong shape: indexing it as if it had the right one does not fail, it
      !> returns whatever is there. Unallocated there is nothing to index wrongly, and
      !> whoever passes it on has to declare it ALLOCATABLE and decide what to do when
      !> it is absent.

      !> L is the one slice both paths fill, so both fill the same array. n_channels is 1 for
      !> spinors and 2 for separate channels, and the two cases never overlap.
      IF (this%l_active .OR. l_ch_orb) THEN
         ALLOCATE (this%l0(window%num_bands, window%num_bands, 3, atoms%nat, &
                           this%n_channels, nkc_loc), source=cmplx(0.0, 0.0))
      END IF

      ! collinear jspins=2 (no SOC/noco): the slices above are spinor-only and stay stubs, so
      ! this case has a pass of its own. What it needs differs in kind: L per channel rather
      ! than one matrix over a spinor, and for the combined 2N spin operator (rspauli.1) only
      ! the cross-spin overlap, since that operator cannot be assembled until both channels
      ! have wannierised and both gauges exist.
      IF (.NOT. this%l_active) THEN
         IF (request%l_socop .AND. .NOT. l_ch_soc) CALL judft_error( &
            "melem_coarse: the spin-orbit operator was requested without spin-orbit coupling", &
            hint="remove the operator, or switch on l_soc", calledby="melem_coarse_init")
         IF (request%l_orbmom .AND. .NOT. l_ch_orb) CALL judft_error( &
            "melem_coarse: the orbital operator has no producer in this spin configuration", &
            hint="it needs spinors (l_soc or l_noco), or jspins=2", &
            calledby="melem_coarse_init")
         IF (request%l_spin .AND. .NOT. l_ch_spin) CALL judft_error( &
            "melem_coarse: the spin operator has no producer in this spin configuration", &
            hint="it needs spinors (l_soc or l_noco), or jspins=2 with an <operators_r> block", &
            calledby="melem_coarse_init")
      END IF

      IF (l_ch_spin) THEN
         ALLOCATE (this%x0(window%num_bands, window%num_bands, MAX(1, COUNT(distk == fmpi%irank))), &
                   source=cmplx(0.0, 0.0))
      END IF
   END SUBROUTINE melem_coarse_init

   SUBROUTINE melem_coarse_calc(this, request, window, atoms, input, sym, cell, noco, nococonv, kpts, &
                                stars, enpara, fmpi, vtot, eig_id, distk, vacuum)
      !> One pass over this rank's k-slice, building every requested operator through the
      !> factory. Spin channels that wannierise separately are a loop, not a second pass:
      !> what differs between one spinor and two channels is which index the result carries
      !> and, for the spin operator, which part of it anyone will use.
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_melem_request), INTENT(IN) :: request
      CLASS(t_melem_window), INTENT(IN) :: window
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_potden), INTENT(IN) :: vtot
      INTEGER, INTENT(IN) :: eig_id
      INTEGER, INTENT(IN) :: distk(:)   ! rank owner of each global k (distributes the loop)
      !> Only a film has one, and only the spin operator needs it here: its four blocks are
      !> overlaps over the whole cell, and on a film that is three regions and not two.
      TYPE(t_vacuum), INTENT(IN), OPTIONAL :: vacuum

      TYPE(t_lapw) :: lapw
      TYPE(t_matelements_spin) :: spinop
      TYPE(t_matelements_soc) :: socop
      TYPE(t_matelements_orbital), ALLOCATABLE :: orbop(:, :)   ! (nat, channel)
      TYPE(t_rsoc) :: rsoc
      INTEGER, ALLOCATABLE :: ev_list(:)
      LOGICAL :: l_spinor_records, l_do_spin, l_do_orb, l_do_soc
      INTEGER :: ikpt, itype, il, na, iatom, ib, ch, a

      !> With two channels only the real-space list can be served, since the interpolated
      !> forms need a gauge that does not exist yet; with one, either list can. Asking the
      !> right question per case is what keeps this loop from having to know which is which.
      !> WHICH SYMMETRY OPERATION TAKES EACH ATOM ONTO ITS REPRESENTATIVE, and of which
      !> class, reported once per run.
      !>
      !> abcof builds the muffin-tin coefficients in the LOCAL frame of each atom: it
      !> rotates k+G by transpose(mrot(:,:,invtab(ngopr(na)))) before expanding. For a
      !> SCALAR operator that cancels in the product A* A and there is nothing to do. For
      !> a VECTOR one it does not: its components live in that local frame and have to be
      !> carried to the global one before the sites are summed.
      !>
      !> The class matters as much as the matrix, because p and L do NOT transform alike:
      !>     p' = R.p              (polar)
      !>     L' = det(R).R.L       (axial)
      !> Under inversion det(R) = -1, and that factor is exactly what leaves L alone while
      !> flipping p. A fix that applied R to both alike would break L, which works today.
      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, '(/,a)') 'melem symmetry: how each atom maps onto its representative'
         WRITE (oUnit, '(a)')   '   atom type neq  ngopr invtab invsat  det(R_cart)   R_cart (rows)'
         DO na = 1, atoms%nat
            iatom = 0
            DO itype = 1, atoms%ntype
               IF (na >= atoms%firstAtom(itype) .AND. &
                   na < atoms%firstAtom(itype) + atoms%neq(itype)) iatom = itype
            END DO
            CALL melem_report_sym(sym, cell, atoms, na, iatom)
         END DO
      END IF

      IF (this%n_channels == 2) THEN
         l_do_spin = request%has_op_r('spin')
         l_do_orb  = request%has_op_r('orbital') .OR. request%has_op('orbital')
         l_do_soc  = request%has_op_r('spin_orbit')
         IF (.NOT. (l_do_spin .OR. l_do_orb .OR. l_do_soc)) RETURN
      ELSE
         IF (.NOT. this%l_active) RETURN   ! nothing requested, or no spinors -> slices are stubs
         l_do_spin = request%l_spin
         l_do_orb  = request%l_orbmom
         l_do_soc  = request%l_socop
      END IF

      !> The relativistic radial SOC integrals and the L.S angular matrix depend on the
      !> potential and the quantisation axis, not on k, so they are built once here. The
      !> angular part is evaluated on the axis the calculation is quantised along.
      IF (l_do_soc) THEN
         !> The SOC operator distributes its column band index over the eigenvector
         !> sub-communicator, while this pass gives every rank whole matrices for its own
         !> k-points. With n_size > 1 it would fill only part of each column block.
         IF (fmpi%n_size /= 1) CALL judft_error( &
            "melem_coarse_calc: the SOC operator needs whole matrices per rank", &
            hint="run k-parallel (n_size = 1); eigenvector parallelism is not supported here", &
            calledby="melem_coarse_calc")
         CALL rsoc%init(atoms)
         CALL rsoc%rad_matrix(atoms, noco, nococonv, input, fmpi, enpara, vtot)
         CALL rsoc%angles(atoms, fmpi, nococonv%theta, nococonv%phi)
      END IF

      !> Set up once, outside the k loop: what an instance binds to -- a site, and a channel
      !> when there is more than one -- is the same at every k. The k dependence arrives with
      !> the coefficients, and init_mat clears the result and reuses the allocation.
      IF (l_do_orb) THEN
         ALLOCATE (orbop(atoms%nat, this%n_channels))
         DO ch = 1, this%n_channels
            na = 0
            DO itype = 1, atoms%ntype
               DO iatom = 1, atoms%neq(itype)
                  na = na + 1
                  !> melem_rot_axial, not melem_rot_polar: L is axial and carries the
                  !> det(R). For a representative atom this comes back as the identity and
                  !> init then leaves the arithmetic alone.
                  IF (this%n_channels == 2) THEN
                     CALL orbop(na, ch)%init(atoms, itype, iatom, channel=ch, &
                                             rot=melem_rot_axial(sym, cell, na))
                  ELSE
                     CALL orbop(na, ch)%init(atoms, itype, iatom, &
                                             rot=melem_rot_axial(sym, cell, na))
                  END IF
               END DO
            END DO
         END DO
      END IF

      !> One band window for every operator at every k, in the form the factory selects with.
      ev_list = [(ib, ib = window%min_band, window%max_band)]

      !> Records 1 and 2 hold the two halves of one spinor whenever the eigenvectors are not
      !> already stored as whole 2N spinors, which is the l_soc=T, l_noco=F case. Two channels
      !> are two states, never halves of one, so they are never stacked.
      l_spinor_records = this%n_channels == 1 .AND. .NOT. noco%l_noco

      il = 0
      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE   ! this rank only computes its own k-slice
         il = il + 1                            ! local (ascending global-k) index for the reduce
         !> The basis at this k. It is built here rather than inside the factory because the
         !> operators need it at their own init; the states themselves the factory reads.
         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)

         IF (l_do_spin) THEN
            !> The operator keeps the four spin blocks either way. Over a spinor the three
            !> Pauli components follow from them and are what the export wants; over two
            !> channels only the cross-spin block is wanted, since the combined 2N operator
            !> cannot be assembled until both channels have wannierised.
            IF (input%film) THEN
               IF (.NOT.PRESENT(vacuum)) CALL judft_error( &
                  'melem_coarse: the spin operator on a film needs the vacuum', &
                  hint='pass vacuum down to melem_coarse_calc', calledby='melem_coarse_calc')
               CALL spinop%init(atoms, stars, lapw, nococonv, input, noco, &
                                vacuum=vacuum, cell=cell, enpara=enpara, vtot=vtot)
            ELSE
               CALL spinop%init(atoms, stars, lapw, nococonv, input, noco)
            END IF
            CALL matrix_element_factory(spinop, eig_id, ikpt, input, atoms, sym, cell, noco, &
                                        nococonv, enpara, lapw, vtot, fmpi, ev_list=ev_list, &
                                        l_both_spinors=l_spinor_records, kpts=kpts)
            !> The generic invariants, once per run. On k = 1 and not on the first LOCAL
            !> k, so the rank that owns it reports and the others stay quiet.
            IF (ikpt == 1) CALL melem_check_provider(spinop, 'spin', ikpt)
            IF (this%n_channels == 1) THEN
               CALL melem_pauli_from_blocks(spinop%mat(1,1)%data_c, spinop%mat(2,2)%data_c, &
                                            spinop%mat(1,2)%data_c, spinop%mat(2,1)%data_c, &
                                            this%s0(:, :, :, il))
               !> The provider check cleared the four spin blocks; this one clears the Pauli
               !> combination built out of them, which is what leaves this module. A
               !> transposed block index survives the first check and fails here.
               IF (ikpt == 1) CALL melem_check_matrix(this%s0(:, :, :, il), 'spin (assembled)', ikpt)
               IF (ikpt <= 3) CALL melem_spin_sumrule(this%s0(:, :, :, il), &
                                                      spinop%mat(1,1)%data_c, spinop%mat(2,2)%data_c, &
                                                      ikpt, tol=1.0e-3)
            ELSE
               this%x0(:, :, il) = spinop%mat(1, 2)%data_c
               !> Not Hermitian on its own -- the cross-spin block is the adjoint of the
               !> other one, not of itself -- so only finite and non-zero are asked of it.
               IF (ikpt == 1) CALL melem_check_matrix(this%x0(:, :, il:il), 'spin cross-block', &
                                                     ikpt, l_hermitian=.FALSE.)
            END IF
         END IF

         IF (l_do_orb) THEN
            !> L is spin-diagonal, so a channel has its own and neither needs the other.
            !> Each instance covers one site and one channel, and whoever wants a total sums
            !> the sites -- wgauge_run.F90 does SUM(l0, DIM=4).
            !>
            !> That sum DOES need the local-to-global rotation, and the note that used to
            !> stand here said it does not. The claim holds for inversion only: L is axial,
            !> L' = det(R) . R . L, and for R = -I that is the identity. Under a proper
            !> rotation the components mix and the site sum is wrong without it.
            DO ch = 1, this%n_channels
               DO na = 1, atoms%nat
                  CALL matrix_element_factory(orbop(na, ch), eig_id, ikpt, input, atoms, sym, &
                                              cell, noco, nococonv, enpara, lapw, vtot, fmpi, &
                                              ev_list=ev_list, l_both_spinors=l_spinor_records, &
                                              kpts=kpts)
                  this%l0(:, :, 1:3, na, ch, il) = orbop(na, ch)%comp(:, :, 1:3)
                  IF (ikpt == 1 .AND. ch == 1) THEN
                     CALL melem_check_provider(orbop(na, ch), 'orbital', ikpt)
                     CALL melem_check_matrix(this%l0(:, :, 1:3, na, ch, il), &
                                             'orbital (assembled)', ikpt)
                  END IF
                  !> How |L|^2 is shared between the three components, SITE BY SITE.
                  !>
                  !> Hermiticity is no use here: the L of each site in its own local frame
                  !> is already hermitian, and a real orthogonal rotation preserves that.
                  !> The momentum was only detectable because its anti-hermitian part is a
                  !> surface term that has to cancel AGAINST ANOTHER REGION, and L is pure
                  !> muffin-tin.
                  !>
                  !> What a lost rotation does move is the SHARE: a rotation mixes the
                  !> components and conserves their sum. So the same crystal described with
                  !> the equivalent atoms grouped (neq=2) and split (neq=1,1) must give the
                  !> SAME total and, if the rotation is missing, a DIFFERENT share.
                  IF (ikpt == 1 .AND. ch == 1) THEN
                     WRITE (oUnit, '(a,i3,a,3es13.5,a,es13.5)') &
                        'melem orbital |L_a|_F  site  ', na, ' :', &
                        (SQRT(SUM(ABS(this%l0(:, :, a, na, ch, il))**2)), a=1, 3), &
                        '   total ', SQRT(SUM(ABS(this%l0(:, :, 1:3, na, ch, il))**2))
                  END IF
               END DO
               !> The SITE SUM, which is the quantity the consumers actually take and the
               !> only one here that a missing rotation can move.
               !>
               !> The per-site norms above cannot see it. A C3 about z mixes x with y, and
               !> in the test cell |L_x| and |L_y| agree to five digits, so each norm comes
               !> back unchanged whether the frames were rotated or not -- the metric is
               !> nearly invariant under the very operation it is meant to probe. The sum
               !> over sites has no such symmetry: if sites 2 and 3 are left in their local
               !> frames, their contributions add along the wrong axes.
               IF (ikpt == 1 .AND. ch == 1) THEN
                  WRITE (oUnit, '(a,3es13.5,a,es13.5)') &
                     'melem orbital |SUM_sites L_a|_F     :', &
                     (SQRT(SUM(ABS(SUM(this%l0(:, :, a, :, ch, il), DIM=3))**2)), a=1, 3), &
                     '   total ', &
                     SQRT(SUM(ABS(SUM(this%l0(:, :, 1:3, :, ch, il), DIM=4))**2))
               END IF
            END DO
         END IF

         IF (l_do_soc) THEN
            !> The operator keeps the four spin blocks. A spinor wavefunction has both
            !> components, so its expectation value of a spinor operator is the sum of all
            !> four; the blocks themselves are what the real-space export carries.
            CALL socop%init(atoms, noco, input, sym, cell, enpara, lapw, vtot, rsoc, fmpi, nococonv)
            CALL matrix_element_factory(socop, eig_id, ikpt, input, atoms, sym, cell, noco, &
                                        nococonv, enpara, lapw, vtot, fmpi, ev_list=ev_list, &
                                        l_both_spinors=l_spinor_records, kpts=kpts)
            IF (ikpt == 1) CALL melem_check_provider(socop, 'spin_orbit', ikpt)
            !> The blocks arrive as conj(<i|H_so|j>) and are restored to the convention
            !> the type documents here, where the real-space export consumes them.
            this%soc4(:, :, 1, il) = CONJG(socop%mat(1, 1)%data_c)
            this%soc4(:, :, 2, il) = CONJG(socop%mat(1, 2)%data_c)
            this%soc4(:, :, 3, il) = CONJG(socop%mat(2, 1)%data_c)
            this%soc4(:, :, 4, il) = CONJG(socop%mat(2, 2)%data_c)
            !> Only for a spinor: the sum over the blocks is the expectation value of a
            !> spinor operator, and with two channels there is no spinor to take it over.
            IF (this%l_active) this%soc0(:, :, 1, il) = CONJG(socop%mat(1, 1)%data_c + socop%mat(1, 2)%data_c &
                                   + socop%mat(2, 1)%data_c + socop%mat(2, 2)%data_c)
            !> The four blocks and their sum at one k, for the grouped-against-split
            !> control.
            !>
            !> SOC is the one operator here whose coefficients arrive in the LOCAL frame
            !> while its angular table soangl is built in the GLOBAL one and carries no
            !> atom index. secvar_soc avoids the mismatch by forcing sym%ngopr = 1 before
            !> the coefficients are made; this path hands every operator the same sym, and
            !> the other operators now need the local frame. Every SOC case validated so
            !> far has a single atom, where that rotation is the identity and the question
            !> cannot arise.
            !>
            !> Reading this: the DIAGONAL spin blocks set the noise floor and the CROSS
            !> blocks carry the signal. A local-frame rotation about z inserts exp(-i m a)
            !> into the lm coefficients; L_z S_z connects m to m and the phases cancel,
            !> while L_+ S_- and L_- S_+ connect m to m+-1 and leave exp(-+i a) behind.
            !> Tilting the spin axis off the rotation axis is not an option anyway: the
            !> rotation then leaves the symmetry group and the atoms stop being equivalent.
            IF (ikpt == 1) THEN
               WRITE (oUnit, '(a,4es13.5)') &
                  'melem soc |H_so block|_F 11 12 21 22 :', &
                  (SQRT(SUM(ABS(this%soc4(:, :, a, il))**2)), a=1, 4)
               IF (this%l_active) WRITE (oUnit, '(a,es13.5)') &
                  'melem soc |H_so summed|_F            :', &
                  SQRT(SUM(ABS(this%soc0(:, :, 1, il))**2))
            END IF
            IF (ikpt == 1 .AND. this%l_active) THEN
               !> The sum over the four blocks is Hermitian even though two of them are not:
               !> 12 + 21 is, being a matrix plus its adjoint.
               CALL melem_check_matrix(this%soc0(:, :, 1:1, il), 'spin_orbit (assembled)', ikpt)
               CALL melem_check_matrix(this%soc4(:, :, :, il), 'spin_orbit blocks', ikpt, &
                                       l_hermitian=.FALSE.)
            END IF
         END IF
      END DO

      IF (ALLOCATED(orbop)) DEALLOCATE (orbop)
      IF (ALLOCATED(ev_list)) DEALLOCATE (ev_list)
   END SUBROUTINE melem_coarse_calc

   SUBROUTINE melem_coarse_free(this)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this

      IF (ALLOCATED(this%s0)) DEALLOCATE (this%s0)
      IF (ALLOCATED(this%l0)) DEALLOCATE (this%l0)
      IF (ALLOCATED(this%soc0)) DEALLOCATE (this%soc0)
      IF (ALLOCATED(this%soc4)) DEALLOCATE (this%soc4)
      IF (ALLOCATED(this%x0)) DEALLOCATE (this%x0)
      this%n_channels = 1
      !> l_active gates both the allocations above and the calc pass, so leaving it set
      !> would describe slices that are no longer there.
      this%l_active = .FALSE.
   END SUBROUTINE melem_coarse_free

   !>  One line per atom: the operation that takes it onto its representative, in cartesian
   !>  form, with its determinant. det = +1 is a proper rotation, which mixes the components
   !>  of p and of L alike; det = -1 is improper, and there p changes sign relative to L.
   !>
   !>  invsat is reported too because the local-to-global rotation of the spin-orbit
   !>  coefficients takes a different branch for an atom reached through inversion
   !>  (invsat == 2): it borrows its partner's operation and picks up (-1)**l.
   SUBROUTINE melem_report_sym(sym, cell, atoms, na, itype)
      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER,       INTENT(IN) :: na, itype
      INTEGER :: nap, inap, i, isat
      REAL    :: rl(3,3), rc(3,3), ainv(3,3), det

      nap  = sym%ngopr(na)
      inap = sym%invtab(nap)
      isat = -1
      IF (ALLOCATED(sym%invsat)) isat = sym%invsat(na)
      rl = REAL(sym%mrot(:, :, nap))
      !> mrot esta en coordenadas de RED. Un vector cartesiano se transforma con
      !> amat . mrot . amat^-1, y amat lleva los vectores de red como COLUMNAS.
      ainv = inverse3(cell%amat)
      rc = MATMUL(cell%amat, MATMUL(rl, ainv))
      det = rc(1,1)*(rc(2,2)*rc(3,3) - rc(2,3)*rc(3,2)) &
          - rc(1,2)*(rc(2,1)*rc(3,3) - rc(2,3)*rc(3,1)) &
          + rc(1,3)*(rc(2,1)*rc(3,2) - rc(2,2)*rc(3,1))
      WRITE (oUnit, '(3i5,3i7,f12.4,a,3f8.4)') na, itype, atoms%neq(itype), nap, inap, isat, det, &
         '   ', (rc(1, i), i = 1, 3)
      WRITE (oUnit, '(51x,a,3f8.4)') '   ', (rc(2, i), i = 1, 3)
      WRITE (oUnit, '(51x,a,3f8.4)') '   ', (rc(3, i), i = 1, 3)
   END SUBROUTINE melem_report_sym

   PURE FUNCTION inverse3(a) RESULT(b)
      REAL, INTENT(IN) :: a(3,3)
      REAL :: b(3,3), d
      d = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
      b(1,1) =  (a(2,2)*a(3,3) - a(2,3)*a(3,2))/d
      b(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))/d
      b(1,3) =  (a(1,2)*a(2,3) - a(1,3)*a(2,2))/d
      b(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))/d
      b(2,2) =  (a(1,1)*a(3,3) - a(1,3)*a(3,1))/d
      b(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))/d
      b(3,1) =  (a(2,1)*a(3,2) - a(2,2)*a(3,1))/d
      b(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))/d
      b(3,3) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1))/d
   END FUNCTION inverse3

END MODULE m_melem_coarse
