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
!>  Everything that needs the gauge lives in m_melem_run; the collinear combined spin
!>  operator, which needs both channels wannierised, in m_melem_spin_collinear.
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
   USE m_types_stars
   USE m_types_mat
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_manifold, ONLY: t_melem_manifold
   USE m_melem_spin, ONLY: melem_pauli_from_blocks, melem_spin_sumrule
   USE m_types_matelements_spin, ONLY: t_matelements_spin
   USE m_types_matelements_soc, ONLY: t_matelements_soc
   USE m_types_rsoc, ONLY: t_rsoc
   USE m_types_matelements_orbital, ONLY: t_matelements_orbital
   USE m_matrix_element_factory, ONLY: matrix_element_factory
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
      !> collinear jspins=2 only: the gauge V of each spin channel, needed by the combined 2N
      !> spin operator.
      COMPLEX, ALLOCATABLE :: v_ch(:, :, :, :)    !< (nb,nw,nkptf,2)
      !> How many spin channels wannierise separately: two when jspins=2 without spinors,
      !> one otherwise. It is a fact about the calculation, where the flags it replaces were
      !> a copy of what the request already said, kept next to it and able to disagree.
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

   SUBROUTINE melem_coarse_init(this, request, manifold, atoms, input, kpts, fmpi, distk, l_spinors)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_melem_request), INTENT(IN) :: request
      TYPE(t_melem_manifold), INTENT(IN) :: manifold
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)
      LOGICAL, INTENT(IN) :: l_spinors   !< noco%l_noco .OR. noco%l_soc

      INTEGER :: nkc_loc
      LOGICAL :: l_ch_orb, l_ch_spin

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
      nkc_loc = MAX(1, COUNT(distk == fmpi%irank))

      IF (this%l_active) THEN
         ALLOCATE (this%s0(manifold%num_bands, manifold%num_bands, 3, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc0(manifold%num_bands, manifold%num_bands, 1, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc4(manifold%num_bands, manifold%num_bands, 4, nkc_loc), source=cmplx(0.0, 0.0))
      END IF
      !> Y si no, no se reservan. Un tapon (1,1,1,1) es un array VALIDO con la forma
      !> equivocada: indexarlo como si tuviera la buena no falla, solo devuelve lo que haya.
      !> Sin reservar no hay nada que indexar mal, y quien lo pase tiene que declararlo
      !> ALLOCATABLE y decidir que hace cuando no esta.

      !> L is the one slice both paths fill, so both fill the same array. n_channels is 1 for
      !> spinors and 2 for separate channels, and the two cases never overlap.
      IF (this%l_active .OR. l_ch_orb) THEN
         ALLOCATE (this%l0(manifold%num_bands, manifold%num_bands, 3, atoms%nat, &
                           this%n_channels, nkc_loc), source=cmplx(0.0, 0.0))
      END IF

      ! collinear jspins=2 (no SOC/noco): the slices above are spinor-only and stay stubs, so
      ! this case has a pass of its own. What it needs differs in kind: L per channel rather
      ! than one matrix over a spinor, and for the combined 2N spin operator (rspauli.1) only
      ! the cross-spin overlap, since that operator cannot be assembled until both channels
      ! have wannierised and their gauges v_ch exist.
      !> An operator nobody will build must not reach the export: the slices stay at their
      !> stub size, the export reads them anyway, and what comes out is small enough to pass
      !> for numerical noise instead of for the absence of a calculation.
      !> Two channels are two eigenproblems, wannierised one after the other, so there is
      !> no single spin matrix over them to interpolate: within a channel sigma_z is +/-1
      !> by orthonormality, and the transverse part lives in the cross-spin block, which
      !> needs BOTH gauges and therefore both wannierizations. It exists only as the
      !> combined 2N operator in real space. Said here because the summary flag is set by
      !> the real-space list too: without this the request passed and the stub-sized slice
      !> reached the interpolation driver, where the shapes do not conform.
      IF (this%n_channels == 2 .AND. request%needs_op('spin', interp_only=.TRUE.)) CALL judft_error( &
         "melem_coarse: the spin operator cannot be interpolated when the two spin "// &
         "channels are wannierised separately", &
         hint="ask for it in <operators_r>: there both channels are combined into the "// &
              "2N rspauli.1, which is the only form this operator has here", &
         calledby="melem_coarse_init")

      IF (.NOT. this%l_active) THEN
         IF (request%l_socop) CALL judft_error( &
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
         ALLOCATE (this%v_ch(manifold%num_bands, manifold%num_wann, kpts%nkptf, 2), source=cmplx(0.0, 0.0))
         ALLOCATE (this%x0(manifold%num_bands, manifold%num_bands, MAX(1, COUNT(distk == fmpi%irank))), &
                   source=cmplx(0.0, 0.0))
      END IF
   END SUBROUTINE melem_coarse_init

   SUBROUTINE melem_coarse_calc(this, request, manifold, atoms, input, sym, cell, noco, nococonv, kpts, &
                                stars, enpara, fmpi, vtot, eig_id, distk)
      !> One pass over this rank's k-slice, building every requested operator through the
      !> factory. Spin channels that wannierise separately are a loop, not a second pass:
      !> what differs between one spinor and two channels is which index the result carries
      !> and, for the spin operator, which part of it anyone will use.
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_melem_request), INTENT(IN) :: request
      TYPE(t_melem_manifold), INTENT(IN) :: manifold
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

      TYPE(t_lapw) :: lapw
      TYPE(t_matelements_spin) :: spinop
      TYPE(t_matelements_soc) :: socop
      TYPE(t_matelements_orbital), ALLOCATABLE :: orbop(:, :)   ! (nat, channel)
      TYPE(t_rsoc) :: rsoc
      INTEGER, ALLOCATABLE :: ev_list(:)
      LOGICAL :: l_spinor_records, l_do_spin, l_do_orb
      INTEGER :: ikpt, itype, il, na, iatom, ib, ch

      !> With two channels only the real-space list can be served, since the interpolated
      !> forms need a gauge that does not exist yet; with one, either list can. Asking the
      !> right question per case is what keeps this loop from having to know which is which.
      IF (this%n_channels == 2) THEN
         l_do_spin = request%has_op_r('spin')
         l_do_orb  = request%has_op_r('orbital') .OR. request%has_op('orbital') .OR. &
                     request%has_op('orbitalCurrent')
         IF (.NOT. (l_do_spin .OR. l_do_orb)) RETURN
      ELSE
         IF (.NOT. this%l_active) RETURN   ! nothing requested, or no spinors -> slices are stubs
         l_do_spin = request%l_spin
         l_do_orb  = request%l_orbmom
      END IF

      !> The relativistic radial SOC integrals and the L.S angular matrix depend on the
      !> potential and the quantisation axis, not on k, so they are built once here. The
      !> angular part is evaluated on the axis the calculation is quantised along.
      IF (request%l_socop) THEN
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
                  IF (this%n_channels == 2) THEN
                     CALL orbop(na, ch)%init(atoms, itype, iatom, channel=ch)
                  ELSE
                     CALL orbop(na, ch)%init(atoms, itype, iatom)
                  END IF
               END DO
            END DO
         END DO
      END IF

      !> One band window for every operator at every k, in the form the factory selects with.
      ev_list = [(ib, ib = manifold%min_band, manifold%max_band)]

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
            CALL spinop%init(atoms, stars, lapw, nococonv, input, noco)
            CALL matrix_element_factory(spinop, eig_id, ikpt, input, atoms, sym, cell, noco, &
                                        nococonv, enpara, lapw, vtot, fmpi, ev_list=ev_list, &
                                        l_both_spinors=l_spinor_records, kpts=kpts)
            IF (this%n_channels == 1) THEN
               CALL melem_pauli_from_blocks(spinop%mat(1,1)%data_c, spinop%mat(2,2)%data_c, &
                                            spinop%mat(1,2)%data_c, spinop%mat(2,1)%data_c, &
                                            this%s0(:, :, :, il))
               IF (ikpt <= 3) CALL melem_spin_sumrule(this%s0(:, :, :, il), &
                                                      spinop%mat(1,1)%data_c, spinop%mat(2,2)%data_c, &
                                                      ikpt, tol=1.0e-3)
            ELSE
               this%x0(:, :, il) = spinop%mat(1, 2)%data_c
            END IF
         END IF

         IF (l_do_orb) THEN
            !> L is spin-diagonal, so a channel has its own and neither needs the other. Each
            !> instance covers one site and one channel; whoever wants a total sums the sites,
            !> and that sum needs no local-to-global rotation, L being spin-diagonal with a
            !> frame-invariant trace.
            DO ch = 1, this%n_channels
               DO na = 1, atoms%nat
                  CALL matrix_element_factory(orbop(na, ch), eig_id, ikpt, input, atoms, sym, &
                                              cell, noco, nococonv, enpara, lapw, vtot, fmpi, &
                                              ev_list=ev_list, l_both_spinors=l_spinor_records, &
                                              kpts=kpts)
                  this%l0(:, :, 1:3, na, ch, il) = orbop(na, ch)%comp(:, :, 1:3)
               END DO
            END DO
         END IF

         IF (request%l_socop) THEN
            !> The operator keeps the four spin blocks. A spinor wavefunction has both
            !> components, so its expectation value of a spinor operator is the sum of all
            !> four; the blocks themselves are what the real-space export carries.
            CALL socop%init(atoms, noco, input, sym, cell, enpara, lapw, vtot, rsoc, fmpi, nococonv)
            CALL matrix_element_factory(socop, eig_id, ikpt, input, atoms, sym, cell, noco, &
                                        nococonv, enpara, lapw, vtot, fmpi, ev_list=ev_list, &
                                        l_both_spinors=l_spinor_records, kpts=kpts)
            this%soc4(:, :, 1, il) = socop%mat(1, 1)%data_c
            this%soc4(:, :, 2, il) = socop%mat(1, 2)%data_c
            this%soc4(:, :, 3, il) = socop%mat(2, 1)%data_c
            this%soc4(:, :, 4, il) = socop%mat(2, 2)%data_c
            this%soc0(:, :, 1, il) = socop%mat(1, 1)%data_c + socop%mat(1, 2)%data_c &
                                   + socop%mat(2, 1)%data_c + socop%mat(2, 2)%data_c
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
      IF (ALLOCATED(this%v_ch)) DEALLOCATE (this%v_ch)
      IF (ALLOCATED(this%x0)) DEALLOCATE (this%x0)
      this%n_channels = 1
   END SUBROUTINE melem_coarse_free

END MODULE m_melem_coarse
