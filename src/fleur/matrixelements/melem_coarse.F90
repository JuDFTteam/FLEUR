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
!>    %calc   ONE shared k-pass -- read_eig + abc per k -- feeding the spin, orbital and
!>            spin-orbit providers, so the eigenvectors are read once for all of them
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
   USE m_types_usdus
   USE m_types_mat
   USE m_types_radfun
   USE m_types_spinor_layout, ONLY: radial_slot, melem_stack_spinor
   USE m_types_abc
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_manifold, ONLY: t_melem_manifold
   USE m_melem_spin, ONLY: melem_pauli_from_blocks, melem_spin_sumrule, melem_spin_mt_block
   USE m_melem_orbmom, ONLY: melem_orbmom_bloch_collinear
   USE m_melem_overlap, ONLY: melem_overlap_interstitial
   USE m_types_matelements_spin, ONLY: t_matelements_spin
   USE m_types_matelements_soc, ONLY: t_matelements_soc
   USE m_types_rsoc, ONLY: t_rsoc
   USE m_types_matelements_orbital, ONLY: t_matelements_orbital
   USE m_melem_get_z, ONLY: melem_get_z
   IMPLICIT NONE
   PRIVATE

   !> The Bloch-basis coarse-mesh operator matrices O^0(k), per-rank slices only: the full
   !> coarse mesh is never materialized. Slice entries are stored in ascending global-k order,
   !> matching the gk_loc convention of the distributed FT-reduce.
   TYPE :: t_melem_coarse
      COMPLEX, ALLOCATABLE :: s0(:, :, :, :)      !< (nb,nb,3,nk_loc)      spin
      COMPLEX, ALLOCATABLE :: l0(:, :, :, :, :)   !< (nb,nb,3,nat,nk_loc)  orbital L per atom
      COMPLEX, ALLOCATABLE :: soc0(:, :, :, :)    !< (nb,nb,1,nk_loc)      SOC
      COMPLEX, ALLOCATABLE :: soc4(:, :, :, :)    !< (nb,nb,4,nk_loc)      2x2 SOC spinor blocks
      !> collinear jspins=2 only: per-channel Bloch orbital L, filled from the wannierization's
      !> own mmn k-loop (see add_collinear_orbital) because the spinor coarse pass does not run.
      COMPLEX, ALLOCATABLE :: l0col(:, :, :, :, :) !< (nb,nb,3,channel,nk_loc)
      !> collinear jspins=2 only: the cross-spin overlap <up|dn> per k, in the Bloch basis.
      !> It is the ingredient of the combined 2N spin operator, and it needs no gauge, so it
      !> is built here with the other coarse matrices even though the operator itself cannot
      !> be assembled until both channels are wannierised.
      COMPLEX, ALLOCATABLE :: x0(:, :, :)         !< (nb,nb,nk_loc)
      !> collinear jspins=2 only: the gauge V of each spin channel, needed by the combined 2N
      !> spin operator.
      COMPLEX, ALLOCATABLE :: v_ch(:, :, :, :)    !< (nb,nw,nkptf,2)
      LOGICAL :: l_col_orb = .FALSE.   !< collinear jspins=2 AND orbital requested in <operators_r>
      LOGICAL :: l_col_spin = .FALSE.  !< collinear jspins=2 AND spin requested in <operators_r>
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

      INTEGER :: nkc_loc, iop

      ! Operator Bloch matrices on the coarse mesh: the k-loop is DISTRIBUTED over ranks (each
      ! its distk slice -> parallel get_z I/O) into per-rank local arrays. Every consumer works
      ! on those slices plus a distributed FT-reduce, so the full mesh is never assembled.
      this%l_active = (request%l_spin .OR. request%l_orbmom .OR. request%l_socop) .AND. l_spinors
      IF (this%l_active) THEN
         nkc_loc = MAX(1, COUNT(distk == fmpi%irank))
         ALLOCATE (this%s0(manifold%num_bands, manifold%num_bands, 3, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%l0(manifold%num_bands, manifold%num_bands, 3, atoms%nat, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc0(manifold%num_bands, manifold%num_bands, 1, nkc_loc), source=cmplx(0.0, 0.0))
         ALLOCATE (this%soc4(manifold%num_bands, manifold%num_bands, 4, nkc_loc), source=cmplx(0.0, 0.0))
      ELSE
         ALLOCATE (this%s0(1, 1, 1, 1)); ALLOCATE (this%l0(1, 1, 1, 1, 1)); ALLOCATE (this%soc4(1, 1, 1, 1))
         ALLOCATE (this%soc0(1, 1, 1, 1))
      END IF

      ! collinear jspins=2 (no SOC/noco): the coarse spin/orbital slices above are spinor-only
      ! (stubs), so the per-channel orbital operators_r builds its own Bloch L in the mmn k-loop
      ! (reusing that loop's abc), and the combined 2N spin operator (rspauli.1) is assembled
      ! after both channels wannierise from their gauges v_ch + the cross-spin overlap.
      this%l_col_orb = .FALSE.; this%l_col_spin = .FALSE.
      IF (input%jspins == 2 .AND. .NOT. l_spinors .AND. request%l_operators_r) THEN
         DO iop = 1, request%n_op_r
            IF (TRIM(request%op_r_name(iop)) == 'orbital') this%l_col_orb = .TRUE.
            IF (TRIM(request%op_r_name(iop)) == 'spin') this%l_col_spin = .TRUE.
         END DO
      END IF
      !> An operator nobody will build must not reach the export: the slices stay at their
      !> stub size, the export reads them anyway, and what comes out is small enough to pass
      !> for numerical noise instead of for the absence of a calculation.
      IF (.NOT. this%l_active) THEN
         IF (request%l_socop) CALL judft_error( &
            "melem_coarse: the spin-orbit operator was requested without spin-orbit coupling", &
            hint="remove the operator, or switch on l_soc", calledby="melem_coarse_init")
         IF (request%l_orbmom .AND. .NOT. this%l_col_orb) CALL judft_error( &
            "melem_coarse: the orbital operator has no producer in this spin configuration", &
            hint="it needs spinors (l_soc or l_noco), or jspins=2 with an <operators_r> block", &
            calledby="melem_coarse_init")
         IF (request%l_spin .AND. .NOT. this%l_col_spin) CALL judft_error( &
            "melem_coarse: the spin operator has no producer in this spin configuration", &
            hint="it needs spinors (l_soc or l_noco), or jspins=2 with an <operators_r> block", &
            calledby="melem_coarse_init")
      END IF

      IF (this%l_col_spin) THEN
         ALLOCATE (this%v_ch(manifold%num_bands, manifold%num_wann, kpts%nkptf, 2), source=cmplx(0.0, 0.0))
         ALLOCATE (this%x0(manifold%num_bands, manifold%num_bands, MAX(1, COUNT(distk == fmpi%irank))), &
                   source=cmplx(0.0, 0.0))
      ELSE
         ALLOCATE (this%v_ch(1, 1, 1, 1)); ALLOCATE (this%x0(1, 1, 1))
      END IF
      ! Both channels are held at once: they are produced in one k-pass before either
      ! wannierization runs, and consumed one channel at a time afterwards.
      IF (this%l_col_orb) THEN
         ALLOCATE (this%l0col(manifold%num_bands, manifold%num_bands, 3, 2, &
                              MAX(1, COUNT(distk == fmpi%irank))), source=cmplx(0.0, 0.0))
      ELSE
         ALLOCATE (this%l0col(1, 1, 1, 1, 1))
      END IF
   END SUBROUTINE melem_coarse_init

   SUBROUTINE melem_coarse_calc(this, request, manifold, atoms, input, sym, cell, noco, nococonv, kpts, &
                                stars, usdus, radfun, enpara, fmpi, vtot, eig_id, l_real_wann, distk)
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
      TYPE(t_usdus), INTENT(IN) :: usdus
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_potden), INTENT(IN) :: vtot
      INTEGER, INTENT(IN) :: eig_id
      LOGICAL, INTENT(IN) :: l_real_wann
      INTEGER, INTENT(IN) :: distk(:)   ! rank owner of each global k (distributes the loop)

      TYPE(t_abc), ALLOCATABLE :: abc_s(:, :)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat) :: zMat(1)
      TYPE(t_matelements_spin) :: spinop
      TYPE(t_matelements_soc) :: socop
      !> One instance per (Cartesian component, atom): L_x and L_y are different
      !> operators, and the site resolution is an output of its own.
      TYPE(t_matelements_orbital), ALLOCATABLE :: orbop(:)
      TYPE(t_rsoc) :: rsoc
      TYPE(t_mat) :: zc(2)   ! the two spinor components when get_z does not stack them
      INTEGER :: ikpt, itype, isp, il, jspin_rad, ic, na, iatom

      !> Row 5 -- two spin channels, no spinors -- has no coarse operator slices, but the
      !> cross-spin overlap its combined spin operator needs is still a k-by-k Bloch matrix.
      IF (this%l_col_spin .OR. this%l_col_orb) &
         CALL melem_coarse_collinear(this, manifold, atoms, input, sym, cell, noco, nococonv, &
                                      kpts, stars, usdus, radfun, eig_id, l_real_wann, distk, fmpi)

      IF (.NOT. this%l_active) RETURN   ! nothing requested, or no spinor wavefunctions -> slices are stubs

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

      !> Set up once, outside the k loop: what an instance binds to -- a component and
      !> a site -- is the same at every k, since L has no interstitial part and so
      !> never needs the plane-wave set of a given k. The k dependence of the matrix
      !> elements arrives with the abc coefficients, once per k; init_mat clears the
      !> result matrices there and reuses the allocation.
      IF (request%l_orbmom) THEN
         ALLOCATE (orbop(atoms%nat))
         na = 0
         DO itype = 1, atoms%ntype
            DO iatom = 1, atoms%neq(itype)
               na = na + 1
               CALL orbop(na)%init(atoms, itype, iatom)
            END DO
         END DO
      END IF

      ALLOCATE (abc_s(2, atoms%ntype))
      il = 0
      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE   ! this rank only computes its own k-slice
         il = il + 1                            ! local (ascending global-k) index for the reduce
         ! Load this k-point's eigenvector(s).
         !   l_noco=T          : get_z returns the whole 2N spinor from record 1.
         !   l_soc=T, l_noco=F : get_z returns only N rows and the two spinor components live
         !     in records 1 and 2, so read both and stack them into the 2N layout the
         !     interstitial expects. Reading record 1 alone leaves the spin-down half unread:
         !     the muffin-tin counts the up block twice and the interstitial addresses a down
         !     block that is not there (non-magnetic Pt then sums to <sigma_z> = +N/2, not 0).
         IF (noco%l_noco) THEN
            CALL melem_get_z(manifold%min_band, manifold%max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                  ikpt, 1, l_real_wann, lapw, zMat(1))
         ELSE
            CALL melem_get_z(manifold%min_band, manifold%max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                  ikpt, 1, l_real_wann, lapw, zc(1))
            CALL melem_get_z(manifold%min_band, manifold%max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, &
                                  ikpt, 2, l_real_wann, lapw, zc(2))
            CALL melem_stack_spinor(zc(1), zc(2), zMat(1))
         END IF
         DO isp = 1, 2
            ! The index handed to calc_abc must belong to the zMat it is given and must
            ! never exceed the width of the radial arrays it indexes.
            jspin_rad = radial_slot(radfun, isp)
            DO itype = 1, atoms%ntype
               CALL abc_s(isp, itype)%init(input, atoms, manifold%num_bands, itype)
               IF (noco%l_noco) THEN
                  CALL abc_s(isp, itype)%calc_abc(input, atoms, sym, cell, lapw, manifold%num_bands, usdus, &
                                                  noco, nococonv, jspin_rad, itype, zMat(1))
               ELSE
                  CALL abc_s(isp, itype)%calc_abc(input, atoms, sym, cell, lapw, manifold%num_bands, usdus, &
                                                  noco, nococonv, jspin_rad, itype, zc(isp))
               END IF
            END DO
         END DO
         IF (request%l_spin) THEN
            !The operator keeps the four spin blocks; the three Pauli components follow
            !from them, so only the blocks are computed here.
            CALL spinop%init(atoms, stars, lapw, nococonv, input, noco)
            CALL spinop%init_mat(manifold%num_bands)
            CALL spinop%calc_matrix_elements(zMat, abc_s, radfun, usdus)
            CALL melem_pauli_from_blocks(spinop%mat(1,1)%data_c, spinop%mat(2,2)%data_c, &
                                         spinop%mat(1,2)%data_c, spinop%mat(2,1)%data_c, &
                                         this%s0(:, :, :, il))
            IF (ikpt <= 3) CALL melem_spin_sumrule(this%s0(:, :, :, il), &
                                                   spinop%mat(1,1)%data_c, spinop%mat(2,2)%data_c, &
                                                   ikpt, tol=1.0e-3)
         END IF
         IF (request%l_orbmom) THEN
            !The site-summed total is a plain sum over the last index, because L needs
            !no local-to-global rotation: it is spin-diagonal and its trace is
            !frame-invariant.
            DO na = 1, atoms%nat
               CALL orbop(na)%init_mat(manifold%num_bands, n_alpha=3)
               CALL orbop(na)%calc_matrix_elements(zMat, abc_s, radfun, usdus)
               this%l0(:, :, 1:3, na, il) = orbop(na)%comp(:, :, 1, 1, 1:3)
            END DO
         END IF
         IF (request%l_socop) THEN
            !The operator keeps the four spin blocks. A spinor wavefunction has both
            !components, so its expectation value of a spinor operator is the sum of all
            !four; the blocks themselves are what the real-space export carries.
            CALL socop%init(atoms, noco, input, sym, cell, enpara, lapw, vtot, rsoc, fmpi, nococonv)
            CALL socop%init_mat(manifold%num_bands)
            CALL socop%calc_matrix_elements(zMat, abc_s, radfun, usdus)
            this%soc4(:, :, 1, il) = socop%mat(1, 1)%data_c
            this%soc4(:, :, 2, il) = socop%mat(1, 2)%data_c
            this%soc4(:, :, 3, il) = socop%mat(2, 1)%data_c
            this%soc4(:, :, 4, il) = socop%mat(2, 2)%data_c
            this%soc0(:, :, 1, il) = socop%mat(1, 1)%data_c + socop%mat(1, 2)%data_c &
                                   + socop%mat(2, 1)%data_c + socop%mat(2, 2)%data_c
         END IF
      END DO

      DEALLOCATE (abc_s)
   END SUBROUTINE melem_coarse_calc

   SUBROUTINE melem_coarse_collinear(this, manifold, atoms, input, sym, cell, noco, nococonv, &
                                      kpts, stars, usdus, radfun, eig_id, l_real_wann, distk, fmpi)
      !> What a collinear jspins=2 calculation needs from the eigenvectors, per k-point of
      !> this rank's slice: interstitial part plus the muffin-tin contraction of the two
      !> channels' matching coefficients. The two channels are separate eigenproblems with
      !> the same basis, so both records are read at every k.
      CLASS(t_melem_coarse), INTENT(INOUT) :: this
      TYPE(t_melem_manifold), INTENT(IN) :: manifold
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_usdus), INTENT(IN) :: usdus
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER, INTENT(IN) :: eig_id
      LOGICAL, INTENT(IN) :: l_real_wann
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi

      TYPE(t_lapw) :: lapw_u, lapw_d
      TYPE(t_mat)  :: zMat_u, zMat_d
      TYPE(t_abc), ALLOCATABLE :: abc_both(:, :)   ! (2,ntype): 1=up, 2=dn
      COMPLEX, ALLOCATABLE :: o_uu(:, :), o_dd(:, :), o_ud(:, :), o_du(:, :)
      INTEGER :: nb, ikpt, il, itype, ch

      nb = manifold%num_bands
      ALLOCATE (abc_both(2, atoms%ntype))
      ALLOCATE (o_uu(nb, nb), o_dd(nb, nb), o_ud(nb, nb), o_du(nb, nb))

      il = 0
      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE
         il = il + 1
         CALL melem_get_z(manifold%min_band, manifold%max_band, eig_id, input, atoms, noco, &
                          nococonv, kpts, sym, cell, ikpt, 1, l_real_wann, lapw_u, zMat_u)
         DO itype = 1, atoms%ntype
            CALL abc_both(1, itype)%init(input, atoms, nb, itype)
            CALL abc_both(1, itype)%calc_abc(input, atoms, sym, cell, lapw_u, nb, usdus, noco, &
                                             nococonv, 1, itype, zMat_u)
         END DO
         CALL melem_get_z(manifold%min_band, manifold%max_band, eig_id, input, atoms, noco, &
                          nococonv, kpts, sym, cell, ikpt, 2, l_real_wann, lapw_d, zMat_d)
         DO itype = 1, atoms%ntype
            CALL abc_both(2, itype)%init(input, atoms, nb, itype)
            CALL abc_both(2, itype)%calc_abc(input, atoms, sym, cell, lapw_d, nb, usdus, noco, &
                                             nococonv, MERGE(1, 2, input%jspins == 1), itype, zMat_d)
         END DO
         IF (this%l_col_spin) THEN
            o_uu = CMPLX(0.0, 0.0); o_dd = CMPLX(0.0, 0.0); o_ud = CMPLX(0.0, 0.0); o_du = CMPLX(0.0, 0.0)
            CALL melem_overlap_interstitial(stars, lapw_u, lapw_d, zMat_u, zMat_d, 0, 0, o_ud)
            CALL melem_spin_mt_block(atoms, abc_both, radfun, o_uu, o_dd, o_ud, o_du)
            this%x0(:, :, il) = o_ud
         END IF
         !> L is spin-diagonal, so each channel has its own and neither needs the other.
         IF (this%l_col_orb) THEN
            DO ch = 1, 2
               CALL melem_orbmom_bloch_collinear(atoms, abc_both(ch, :), radfun, &
                                                 MERGE(1, ch, input%jspins == 1), &
                                                 this%l0col(:, :, :, ch, il))
            END DO
         END IF
      END DO

      DEALLOCATE (abc_both, o_uu, o_dd, o_ud, o_du)
   END SUBROUTINE melem_coarse_collinear

   SUBROUTINE melem_coarse_free(this)
      CLASS(t_melem_coarse), INTENT(INOUT) :: this

      IF (ALLOCATED(this%s0)) DEALLOCATE (this%s0)
      IF (ALLOCATED(this%l0)) DEALLOCATE (this%l0)
      IF (ALLOCATED(this%soc0)) DEALLOCATE (this%soc0)
      IF (ALLOCATED(this%soc4)) DEALLOCATE (this%soc4)
      IF (ALLOCATED(this%l0col)) DEALLOCATE (this%l0col)
      IF (ALLOCATED(this%v_ch)) DEALLOCATE (this%v_ch)
      IF (ALLOCATED(this%x0)) DEALLOCATE (this%x0)
      this%l_col_orb = .FALSE.; this%l_col_spin = .FALSE.
   END SUBROUTINE melem_coarse_free

END MODULE m_melem_coarse
