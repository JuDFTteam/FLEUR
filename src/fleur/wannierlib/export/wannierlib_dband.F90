!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The representation of the space group on the BANDS,
!>
!>     d_band(m, n, iop, ir) = < psi_(m, R k) | R | psi_(n, k) >,   k = the star representative,
!>
!> which is what an external symmetriser needs beside the representation on the Wannier
!> functions: the gauge it rotates is U(k), carrying one band index and one Wannier index,
!> so a rotation of U needs a matrix on each side and only D_wann can be rebuilt from the
!> projection table.
!>
!> HOW THE ROTATED STATE IS BUILT, AND WHY NOT BY ROTATING IT. The operation {R|tau} sends
!> a state at k to one at Rk whose coefficient at the rotated vector R(k+G) is the
!> coefficient at k+G times exp(-2 pi i R(k+G).tau). Applying that to a t_lapw in place --
!> which is what the old wann_kptsrotate did -- rotates the integer labels k1/k2/k3 and
!> leaves gvec, rk, gk and vk behind, and those Cartesian arrays are exactly what abcof
!> reads to build the muffin-tin coefficients: the interstitial would be rotated and the
!> spheres would not, which no norm notices.
!>
!> So the state is not rotated, it is RE-INDEXED. The basis of the image point is built the
!> way FLEUR builds any basis, and the rotated coefficients are written into that basis:
!> rotation maps the G-set of k onto the G-set of Rk, so the two lists hold the same vectors
!> in a different order and the map between them is a permutation. Everything downstream
!> then sees one ordinary basis with one ordinary coefficient vector.
!>
!> WHAT IT REFUSES TO DO, rather than do silently. Local orbitals: their rows of the
!> coefficient vector belong to an ATOM, so the operation permutes them among the atoms of
!> an orbit and the plane-wave permutation above does not reach them. Spinors: a spatial
!> rotation turns the spin as well, through the SU(2) matrix of the rotation, and that
!> factor multiplies the whole thing. Both are named in the error, because a d_band that
!> quietly ignored either would still be unitary and still be block-diagonal.
!>
!> TIME REVERSAL costs one conjugation and no new machinery. T R sends k+G to -R(k+G) and
!> conjugates, so the same permutation runs with the rotation negated and the coefficients
!> conjugated. The translation factor comes out UNCHANGED: conjugating turns
!> exp(-2 pi i R(k+G).tau) into its inverse, and R(k+G) = -(k' + G') for an antiunitary
!> operation, so the two sign changes cancel and one expression serves both halves.
MODULE m_wannierlib_dband
   USE m_juDFT
   USE m_constants, ONLY: oUnit, tpi_const, ImagUnit
   USE m_wannierlib_sym_kstar, ONLY: t_wl_kstar
   USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
   USE m_types_spinor_layout, ONLY: t_spinor_layout, radial_slot
   USE m_matrix_element_factory, ONLY: matrix_element_states, matrix_element_abc_of
   USE m_melem_overlap, ONLY: melem_overlap_states
   USE m_melem_ujugaunt, ONLY: melem_ujugaunt
   USE m_types_abc
   USE m_types_radfun
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wl_build_dband

   !> The 2x2 identity, for the places where the spin halves do not mix.
   COMPLEX, PARAMETER :: EYE2(2, 2) = RESHAPE([CMPLX(1.0, 0.0), CMPLX(0.0, 0.0), &
                                               CMPLX(0.0, 0.0), CMPLX(1.0, 0.0)], [2, 2])

   !> The radial integrals of the b = 0 kernel, one table per LOCAL spin direction.
   !>
   !> One table was enough while one potential served both components, which is a
   !> spin-orbit run. A non-collinear one has two potentials and two radial sets, and a
   !> component contracted against the wrong one is a silent error of the size of the
   !> exchange splitting. Keyed on the direction rather than rebuilt inside the k loop,
   !> because b = 0 leaves the integrals independent of k.
   TYPE t_ujug_slot
      COMPLEX, ALLOCATABLE :: t(:, :, :, :, :, :)
   END TYPE t_ujug_slot

CONTAINS

   !> Every rank computes the whole table. The eigenvectors of an arbitrary k are what this
   !> needs, every rank asks for them in the same order, and the table is small next to the
   !> mmn; distributing it would buy a collective and a reason for the ranks to disagree.
   SUBROUTINE wl_build_dband(star, manifold, atoms, cell, input, kpts, sym, noco, nococonv, &
                             enpara, vtot, stars, fmpi, radfun, eig_id, jspin, l_spinors, &
                             dband)
      USE m_types_atoms; USE m_types_cell; USE m_types_input; USE m_types_kpts
      USE m_types_sym; USE m_types_noco; USE m_types_nococonv; USE m_types_enpara
      USE m_types_potden; USE m_types_stars; USE m_types_mpi; USE m_types_lapw
      USE m_types_mat
      USE m_dwigner, ONLY: d_wigner
      TYPE(t_wl_kstar), INTENT(IN) :: star
      TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_potden), INTENT(IN) :: vtot
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER, INTENT(IN) :: eig_id, jspin
      LOGICAL, INTENT(IN) :: l_spinors
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: dband(:, :, :, :)   !> (nb,nb,nsym,nkirr)

      TYPE(t_lapw) :: lapw_irr, lapw_img
      TYPE(t_mat) :: z_irr(2), z_rot(2)
      TYPE(t_mat), POINTER :: zmat_p(:)
      TYPE(t_abc), POINTER :: abc_p(:, :)
      TYPE(t_abc), ALLOCATABLE :: abc_rot(:, :), abc_irr(:, :), abc_one(:)
      COMPLEX, ALLOCATABLE :: dwgn(:, :, :, :)
      TYPE(t_spinor_layout) :: layout, layout_irr
      COMPLEX, ALLOCATABLE :: oblk(:, :, :, :)
      TYPE(t_ujug_slot) :: ujug(2)
      INTEGER, ALLOCATABLE :: ev_list(:), perm(:)
      REAL :: kdiff0(3, 1), rrot(3, 3), kimg(3)
      INTEGER :: ir, iop, jop, kirr, jk, nb, ginv(3, 3), g0(3), j, rsign
      INTEGER :: js, js2, jrec, krec, n_comp, jsp_loc(2), n_rad, ja, jb, i0, j0, nh
      COMPLEX :: usp(2, 2), amt(2, 2), fspin(2, 2, atoms%ntype)
      TYPE(t_mat) :: zperm
      LOGICAL :: l_anti, l_stack, l_two_recs
      REAL :: dev_u, dev_a, gap_lo

      nb = manifold%num_bands
      !> Two components to carry when the states are spinors, one otherwise. Everything
      !> below is written for both: the orbital half of the rotation acts on each component
      !> alone and the spin half is the 2x2 that mixes them, which is why the two can be
      !> kept apart all the way to the overlap.
      n_comp = MERGE(2, 1, l_spinors)
      !> A spinor arrives one of two ways, and the difference decides where the spin
      !> rotation goes. After a spin-orbit second variation the eig file holds TWO records,
      !> one per component, and the components are independent all the way to the overlap:
      !> the 2x2 can wait until the end. A non-collinear Hamiltonian was 2N x 2N and ONE
      !> record holds both halves, and there they are not independent -- calc_abc builds
      !> the coefficients of a local spin direction out of BOTH halves, weighted by the
      !> atom's own frame. So the 2x2 has to be applied to the eigenvector, before the
      !> coefficients are rebuilt, and the overlap is then a plain sum over the directions.
      l_stack = l_spinors .AND. noco%l_noco
      l_two_recs = l_spinors .AND. .NOT. noco%l_noco
      !> The local spin direction each component carries. It selects the radial set and,
      !> in the non-collinear case, the frame. A spinor's component IS that direction
      !> whichever way the spinor arrived: after a collinear second variation the two
      !> components are the two channels, and with jspins=2 those have their own potentials
      !> and so their own radial sets. radial_slot clamps, so a spinor over a single
      !> potential still reads slot 1 and the jspins=1 case is unchanged. Only a scalar
      !> state takes the channel -- which is what the coefficients a few lines below have
      !> always done, and contracting them with another direction's radial set is what
      !> cost WannFeBccSOC and WannFeAFMColSOC four parts in a hundred.
      DO js = 1, n_comp
         jsp_loc(js) = MERGE(js, radial_slot(radfun, MERGE(js, jspin, l_spinors)), l_stack)
      END DO
      ALLOCATE (dband(nb, nb, star%nsym, star%nkirr), source=CMPLX(0.0, 0.0))
      !> <psi^s(Rk)|rot psi^s'(k)>, the overlap with NO spin rotation in it. The spin
      !> rotation enters the overlap linearly, so the table a reader wants is one 2x2
      !> contraction of these four blocks, and the convention sits in that one line
      !> instead of being spread through the rotation.
      ALLOCATE (oblk(nb, nb, 2, 2))
      ALLOCATE (abc_rot(2, atoms%ntype), abc_irr(2, atoms%ntype))
      ev_list = [(j, j=manifold%min_band, manifold%max_band)]

      !> b = 0: the two sides sit at the same k, so the muffin-tin table is the one entry
      !> the neighbour list never holds, built the way the M(k,k) invariant builds it.
      !> Wigner matrices for the whole group, over the full lmax rather than the lmaxU the
      !> sym carries: these rotate the muffin-tin coefficients of every l the basis has.
      ALLOCATE (dwgn(-atoms%lmaxd:atoms%lmaxd, -atoms%lmaxd:atoms%lmaxd, &
                     0:atoms%lmaxd, sym%nop))
      CALL d_wigner(sym%nop, sym%mrot, cell%bmat, atoms%lmaxd, dwgn(:, :, 1:, :))
      dwgn(:, :, 0, :) = CMPLX(1.0, 0.0)
      dev_u = 0.0
      dev_a = 0.0
      gap_lo = 0.0

      kdiff0 = 0.0
      n_rad = radial_slot(radfun, 2)
      DO ja = 1, n_rad
         CALL melem_ujugaunt(atoms, cell, 1, kdiff0, radfun, radfun, ja, ja, .FALSE., 1, &
                             ujug(ja)%t)
      END DO

      CALL timestart("wl_build_dband")
      DO ir = 1, star%nkirr
         kirr = star%kirr(ir)
         CALL lapw_irr%init(input, noco, nococonv, kpts, atoms, sym, kirr, cell)
         CALL matrix_element_states(eig_id, kirr, input, atoms, sym, cell, noco, nococonv, &
                                    enpara, lapw_irr, vtot, fmpi, zmat_p, abc_p, &
                                    ev_list=ev_list, l_both_spinors=l_two_recs, kpts=kpts)
         !> Copied, not pointed at: the slot it lives in is keyed on k and the image points
         !> below go through the same cache, which holds three. The rule below reads the
         !> representative's muffin-tin coefficients for the same reason, while the cache
         !> has already moved on to the image point.
         !> Which record is which: the two records of a spinor are its two components,
         !> while two collinear channels are unrelated states and only this one is wanted.
         CALL layout_irr%init(input, noco, lapw_irr, atoms, l_both_spinors=l_two_recs)
         !> A stacked record is copied once and whole; the coefficients come per local
         !> direction either way, because the factory builds a set per component.
         IF (l_stack) z_irr(1) = zmat_p(1)
         DO js = 1, n_comp
            jrec = MERGE(js, jspin, l_two_recs)
            IF (.NOT. l_stack) z_irr(js) = zmat_p(jrec)
            abc_irr(js, :) = abc_p(MERGE(js, jspin, l_spinors), :)
         END DO

         DO iop = 1, star%nsym
            jk = star%kirr2k(ir, iop)
            !> The image point's own basis, in FLEUR's own order. The rotated state is
            !> written into it, so nothing downstream sees a rotated anything.
            CALL lapw_img%init(input, noco, nococonv, kpts, atoms, sym, jk, cell)
            CALL matrix_element_states(eig_id, jk, input, atoms, sym, cell, noco, nococonv, &
                                       enpara, lapw_img, vtot, fmpi, zmat_p, abc_p, &
                                       ev_list=ev_list, l_both_spinors=l_two_recs, kpts=kpts)
            CALL layout%init(input, noco, lapw_img, atoms, l_both_spinors=l_two_recs)

            !> A k vector transforms under the inverse transpose, which is how the stars
            !> were built: operation iop here is the same operation iop there. Past nop
            !> the operation is that same rotation with time reversal, so the sign flips
            !> and the coefficients get conjugated below.
            l_anti = star%antiunit(iop)
            jop = star%op(iop)
            rsign = MERGE(-1, 1, l_anti)
            ginv = sym%mrot(:, :, sym%invtab(jop))
            rrot = rsign*TRANSPOSE(REAL(ginv))
            !> AGAINST THE BASIS'S OWN k, NOT kpts%bkf. The two disagree: bkf holds the
            !> mesh in [0,1) while lapw%init builds the basis about the same point written
            !> in [-1/2,1/2), so for a k such as 2/3 they differ by a lattice vector. The
            !> G labels below belong to the basis, so folding against bkf puts that lattice
            !> vector into g0 and the permutation lands one shell away. Every mesh of the
            !> test set is 2x2x2, where the two conventions agree at 0 and at 1/2 and the
            !> mistake cannot show; an odd mesh is what brings it out -- measured on a 3x3x3
            !> mesh, where 58 of 148 vectors had no image.
            kimg = MATMUL(rrot, lapw_irr%bkpt)
            !> Rotating lands on the image point up to a reciprocal lattice vector. Folding
            !> it away here leaves both sides at one k and spares the overlap a gb. The
            !> label moves WITH g0, not against it: the physical vector k'+G' rewritten
            !> about bkf(jk) is bkf(jk) + (G' + g0). The other sign puts the vector 2*g0
            !> away and wl_gmap stops the run, which is how this was found.
            g0 = NINT(kimg - lapw_img%bkpt)
            !> The fold has to be exact: a rotated k that is not the image point up to a
            !> lattice vector means the stars and the basis disagree about which k this is.
            IF (ANY(ABS(kimg - lapw_img%bkpt - g0) > 1e-9)) CALL juDFT_error( &
               "wannierlib: a rotated k-point is not its image up to a lattice vector", &
               calledby="wl_build_dband")

            CALL wl_gmap(lapw_irr, lapw_img, ginv, g0, rsign, perm)
            !> The orbital half, component by component. One permutation of reciprocal
            !> vectors serves both: a spinor's two components are expanded in the SAME
            !> basis, so what differs between them is the coefficients and nothing else.
            !> The lattice operation drags the spin only where spin-orbit coupling ties
            !> the two together. Without it they are independent groups and the SU(2) of
            !> the rotation does not belong here: sympsi asks grp_k for those matrices
            !> under the same condition. Time reversal flips the spin regardless, so the
            !> antiunitary factor stays. Measured on hcp Co: putting the rotation in where
            !> it does not belong costs 9.9e-01 with the moments at beta = pi/2, where a
            !> spurious rotation mixes the two halves of the spinor the most, and 6.4e-02
            !> with them along z, where it is only a phase on a pure state.
            IF (l_spinors .AND. (noco%l_soc .OR. (.NOT. noco%l_noco))) THEN
               CALL wl_su2(sym, cell, jop, l_anti, usp)
            ELSE IF (l_spinors .AND. l_anti) THEN
               usp = CMPLX(0.0, 0.0)
               usp(1, 2) = CMPLX(-1.0, 0.0)
               usp(2, 1) = CMPLX(1.0, 0.0)
            ELSE
               usp = EYE2
            END IF

            IF (l_stack) THEN
               !> Permute first, into the half each component belongs to, and mix the
               !> halves afterwards: z_rot^s = sum_s' U(s,s') (perm z)^s'. Two passes and
               !> not one, because every component of the mixture reads every component of
               !> the permuted copy.
               nh = layout%row_dn
               CALL zperm%init(.FALSE., z_irr(1)%matsize1, nb)
               CALL zperm%clear()
               CALL z_rot(1)%init(.FALSE., z_irr(1)%matsize1, nb)
               CALL z_rot(1)%clear()
               DO js = 1, 2
                  CALL wl_rotate_coeffs(lapw_irr, lapw_img, perm, sym%tau(:, jop), l_anti, &
                                        z_irr(1), zperm, layout_irr%row_offset(js), &
                                        layout%row_offset(js))
               END DO
               DO js = 1, 2
                  i0 = layout%row_offset(js)
                  DO js2 = 1, 2
                     j0 = layout%row_offset(js2)
                     z_rot(1)%data_c(i0 + 1:i0 + nh, :) = z_rot(1)%data_c(i0 + 1:i0 + nh, :) &
                        + usp(js, js2)*zperm%data_c(j0 + 1:j0 + nh, :)
                  END DO
               END DO
               CALL zperm%free()
               DO js = 1, n_comp
                  CALL matrix_element_abc_of(input, atoms, sym, cell, noco, nococonv, enpara, &
                                             fmpi, vtot, lapw_img, z_rot(1), nb, jsp_loc(js), &
                                             abc_one)
                  abc_rot(js, :) = abc_one
                  DEALLOCATE (abc_one)
               END DO
            ELSE
               DO js = 1, n_comp
                  CALL z_rot(js)%init(z_irr(js)%l_real, z_irr(js)%matsize1, &
                                      z_irr(js)%matsize2)
                  CALL z_rot(js)%clear()
                  CALL wl_rotate_coeffs(lapw_irr, lapw_img, perm, sym%tau(:, jop), l_anti, &
                                        z_irr(js), z_rot(js))
                  CALL matrix_element_abc_of(input, atoms, sym, cell, noco, nococonv, enpara, &
                                             fmpi, vtot, lapw_img, z_rot(js), nb, jsp_loc(js), &
                                             abc_one)
                  abc_rot(js, :) = abc_one
                  DEALLOCATE (abc_one)
               END DO
            END IF

            !> The frame each type's coefficients are expanded about, folded into the spin
            !> half of the rule. Nothing mixes unless the record was stacked.
            !>
            !> calc_abc builds the local direction j as a^loc_j = sum_s A(s,j) z^s with
            !> A = conjg(nococonv%umat), that is a^loc = transpose(A) z. A is unitary, so
            !> transpose(A) conjg(A) = 1 and z = conjg(A) a^loc. A unitary operation sends
            !> z -> U z and therefore a^loc -> transpose(A) U conjg(A) a^loc; an antiunitary
            !> one sends z -> U conjg(z) = U A conjg(a^loc), leaving transpose(A) U A on the
            !> conjugated coefficients. The unitary factor is the identity for U = 1, which
            !> is what the identity gate asks of it.
            fspin = CMPLX(0.0, 0.0)
            DO j = 1, atoms%ntype
               IF (l_stack) THEN
                  amt = TRANSPOSE(CONJG(nococonv%umat(j)))
                  IF (l_anti) THEN
                     fspin(:, :, j) = MATMUL(MATMUL(amt, usp), TRANSPOSE(amt))
                  ELSE
                     fspin(:, :, j) = MATMUL(MATMUL(amt, usp), CONJG(TRANSPOSE(amt)))
                  END IF
               ELSE
                  fspin(:, :, j) = EYE2
               END IF
            END DO
            CALL wl_rotate_mt(atoms, sym, jop, star%amap(:, iop), l_anti, kimg, dwgn, &
                              n_comp, fspin, abc_irr, abc_rot, dev_u, dev_a, gap_lo)
            !> The spin rotation is inside the eigenvector where the record was stacked, so
            !> what is left to contract with is the identity.
            IF (l_stack) usp = EYE2

            !> The blocks the contraction actually reads, and no others: a diagonal 2x2
            !> leaves half of them out, and a stacked record leaves the off-diagonal ones
            !> out because its rotation already mixed the halves. The representative's side
            !> is the record the cache holds; the rotated side is one matrix per component
            !> unless it is stacked, where both sides are reached by offset.
            oblk = CMPLX(0.0, 0.0)
            DO js = 1, n_comp
               DO js2 = 1, n_comp
                  IF (ABS(usp(js, js2)) < 1e-14) CYCLE
                  jrec = MERGE(1, MERGE(js, jspin, l_two_recs), l_stack)
                  krec = MERGE(1, js2, l_stack)
                  ja = radial_slot(radfun, jsp_loc(js))
                  jb = radial_slot(radfun, jsp_loc(js2))
                  IF (ja /= jb) CALL juDFT_error( &
                     "wannierlib: the b=0 radial table holds no mixed pair of spin sets", &
                     hint="two components reading different potentials cannot overlap here", &
                     calledby="wl_build_dband")
                  CALL melem_overlap_states(stars, atoms, lapw_img, lapw_img, &
                                            zmat_p(jrec), z_rot(krec), &
                                            abc_p(MERGE(js, jspin, l_spinors), :), &
                                            abc_rot(js2, :), ja, jb, &
                                            lapw_img%bkpt, lapw_img%bkpt, [0, 0, 0], &
                                            ujug(ja)%t, kdiff0, 1, &
                                            layout%row_offset(js), &
                                            MERGE(layout%row_offset(js2), 0, l_stack), &
                                            oblk(:, :, js, js2))
               END DO
            END DO
            !> CONJUGATED. The overlap kernel is the mmn's, which puts its conjugation on
            !> the second argument, so what comes back is the conjugate of <a|R|b>. Which of
            !> the two the file should hold is not a matter of taste: the relation a reader
            !> applies is A(Rk) = d_band . A(k) . D_wann^+, and of the SIXTEEN placements of
            !> a conjugate, a transpose and an adjoint over the three matrices exactly one
            !> satisfies it -- measured on Si, 1.4e-12 against 1.4e+00 for the other fifteen.
            !> That one is this conjugate. Unitarity, the identity, the block structure and
            !> the group product table are all blind to it: the conjugate of a representation
            !> is a representation.
            dband(:, :, iop, ir) = CMPLX(0.0, 0.0)
            DO js = 1, n_comp
               DO js2 = 1, n_comp
                  dband(:, :, iop, ir) = dband(:, :, iop, ir) &
                                         + usp(js, js2)*CONJG(oblk(:, :, js, js2))
               END DO
            END DO
            DO js = 1, MERGE(1, n_comp, l_stack)
               CALL z_rot(js)%free()
            END DO
         END DO
      END DO
      CALL timestop("wl_build_dband")

      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, '(/,a)') &
            "Wannier90/WannierBerri symmetry export: muffin-tin rotation rule"
         WRITE (oUnit, '(a,es10.3)') "  vs the rebuilt A and B, unitary ops     ", dev_u
         IF (star%nsym > star%nop) THEN
            WRITE (oUnit, '(a,es10.3)') &
               "  vs the rebuilt A and B, antiunitary ops ", dev_a
         END IF
         IF (atoms%nlotot > 0) WRITE (oUnit, '(a,es10.3)') &
            "  how far the rebuild was where overridden ", gap_lo
      END IF
      IF (fmpi%irank == 0) CALL wl_report_dband(star, sym, dband)
      DEALLOCATE (oblk)
   END SUBROUTINE wl_build_dband

   !> Where each plane wave of the representative's basis lands in the image point's basis.
   !>
   !> Built as a lookup over the integer labels rather than by searching the list: the two
   !> hold the same vectors and nothing says in the same order, so a search is quadratic in
   !> a list of thousands and this is one pass each way.
   !>
   !> It has to come out a permutation, and both halves of that are checked: a vector with
   !> no image means the rotation left the basis, and two vectors sharing an image means the
   !> map is not invertible. Either way the coefficient vector would be built from part of
   !> the state, with a norm that looks right because what is missing is a few of the
   !> shortest vectors.
   SUBROUTINE wl_gmap(lapw_a, lapw_b, ginv, g0, rsign, perm)
      USE m_types_lapw
      TYPE(t_lapw), INTENT(IN) :: lapw_a, lapw_b
      INTEGER, INTENT(IN) :: ginv(3, 3), g0(3), rsign
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: perm(:)

      INTEGER, ALLOCATABLE :: box(:, :, :), seen(:)
      INTEGER :: lo(3), hi(3), g(3), i, nva, nvb

      nva = lapw_a%nv(1)
      nvb = lapw_b%nv(1)
      IF (nva /= nvb) CALL juDFT_error( &
         "wannierlib: a k-point and its symmetry image have different basis sizes", &
         calledby="wl_gmap")

      lo(1) = MINVAL(lapw_b%k1(:nvb, 1)); hi(1) = MAXVAL(lapw_b%k1(:nvb, 1))
      lo(2) = MINVAL(lapw_b%k2(:nvb, 1)); hi(2) = MAXVAL(lapw_b%k2(:nvb, 1))
      lo(3) = MINVAL(lapw_b%k3(:nvb, 1)); hi(3) = MAXVAL(lapw_b%k3(:nvb, 1))
      ALLOCATE (box(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), source=0)
      DO i = 1, nvb
         box(lapw_b%k1(i, 1), lapw_b%k2(i, 1), lapw_b%k3(i, 1)) = i
      END DO

      IF (ALLOCATED(perm)) DEALLOCATE (perm)
      ALLOCATE (perm(nva), source=0)
      DO i = 1, nva
         g = rsign*MATMUL(TRANSPOSE(ginv), &
                          (/lapw_a%k1(i, 1), lapw_a%k2(i, 1), lapw_a%k3(i, 1)/)) + g0
         IF (ANY(g < lo) .OR. ANY(g > hi)) CALL juDFT_error( &
            "wannierlib: a rotated reciprocal vector is outside the image basis", &
            calledby="wl_gmap")
         perm(i) = box(g(1), g(2), g(3))
         IF (perm(i) == 0) CALL juDFT_error( &
            "wannierlib: a rotated reciprocal vector is not in the image basis", &
            calledby="wl_gmap")
      END DO
      DEALLOCATE (box)

      ALLOCATE (seen(nva), source=0)
      DO i = 1, nva
         IF (seen(perm(i)) /= 0) CALL juDFT_error( &
            "wannierlib: two reciprocal vectors rotate onto the same one", &
            calledby="wl_gmap")
         seen(perm(i)) = i
      END DO
      DEALLOCATE (seen)
   END SUBROUTINE wl_gmap

   !> The rotated coefficient vector, in the image point's basis.
   !>
   !> The non-symmorphic factor is the phase the rotated plane wave picks up from the
   !> translation, exp(-2 pi i (k' + G').tau) with k' + G' the vector AFTER folding, which is
   !> the same vector before it. A symmorphic group leaves it at one, so a case whose tau is
   !> all zero does not exercise this line at all.
   SUBROUTINE wl_rotate_coeffs(lapw_a, lapw_b, perm, tau, l_anti, z_a, z_b, ioff_a, ioff_b)
      USE m_types_lapw; USE m_types_mat
      TYPE(t_lapw), INTENT(IN) :: lapw_a, lapw_b
      INTEGER, INTENT(IN) :: perm(:)
      REAL, INTENT(IN) :: tau(3)
      LOGICAL, INTENT(IN) :: l_anti
      TYPE(t_mat), INTENT(IN) :: z_a
      TYPE(t_mat), INTENT(INOUT) :: z_b
      !> Row at which this component's block begins on each side. Zero unless the spinor
      !> arrives as one matrix of both halves, where the permutation of reciprocal vectors
      !> is the same for the two and only the block they live in differs.
      INTEGER, OPTIONAL, INTENT(IN) :: ioff_a, ioff_b

      INTEGER :: i, j, ia, ib
      REAL :: arg
      COMPLEX :: ph

      ia = 0; ib = 0
      IF (PRESENT(ioff_a)) ia = ioff_a
      IF (PRESENT(ioff_b)) ib = ioff_b

      !> Unreachable from wannierlib_main, which refuses d_band for this combination and
      !> writes the stars alone. It stays as the contract of this routine: a caller that
      !> hands over a real z together with a non-symmorphic operation is asking for a
      !> complex phase on a real array, and silently dropping it would look like a gauge.
      IF (z_a%l_real .AND. ANY(ABS(tau) > 1e-12)) CALL juDFT_error( &
         "wannierlib: a real eigenvector cannot carry the non-symmorphic phase", &
         hint="d_band is refused for this combination; see wannierlib_main", &
         calledby="wl_rotate_coeffs")

      DO i = 1, lapw_a%nv(1)
         j = perm(i)
         arg = -tpi_const*DOT_PRODUCT(lapw_b%bkpt + &
                                      (/REAL(lapw_b%k1(j, 1)), REAL(lapw_b%k2(j, 1)), &
                                        REAL(lapw_b%k3(j, 1))/), tau)
         IF (z_a%l_real) THEN
            z_b%data_r(ib + j, :) = z_a%data_r(ia + i, :)   ! conjugating real data is a no-op
         ELSE
            ph = EXP(ImagUnit*arg)
            IF (l_anti) THEN
               z_b%data_c(ib + j, :) = ph*CONJG(z_a%data_c(ia + i, :))
            ELSE
               z_b%data_c(ib + j, :) = ph*z_a%data_c(ia + i, :)
            END IF
         END IF
      END DO
   END SUBROUTINE wl_rotate_coeffs

   !> The muffin-tin coefficients of the rotated state, by the rotation rule instead of from
   !> the rotated eigenvector.
   !>
   !> cof(nu, lm, iOrd, iAtom) holds the A coefficient at iOrd = 1, the B at 2 and the local
   !> orbitals from 3 on, and a rotation acts on the lm index alone -- the radial order is a
   !> spectator. One rule therefore covers all of them, which is the whole point: the A and B
   !> orders can be rebuilt from the rotated eigenvector and the local orbitals cannot. Their
   !> rows of the eigenvector are written in a basis the image point does not share, because
   !> lapw%kvec picks its 2l+1 reciprocal vectors per orbital by a greedy independence test
   !> at each k, and the image point's choice is not the image of this one.
   !>
   !> So the rule is applied to every order, COMPARED against the rebuilt coefficients in
   !> every l that has no local orbital, and used in every l that has one -- the whole
   !> channel, not just the orbital's own order, because the rebuilt A and B of that l are
   !> fed by the same rows. That comparison is the control, and it runs on
   !> every case whether or not it has local orbitals: it is what says the atom map and the
   !> phase are right rather than merely plausible. Unitary and antiunitary operations are
   !> reported apart, because the antiunitary half conjugates the source and whether the
   !> conjugation commutes with the rebuild is a question, not an assumption.
   !>
   !> Map and phase come from where mapatom built the map: R.taual(a) + tau = taual(a') + L
   !> with L a lattice vector, and L is what the Bloch factor of the image point sees.
   SUBROUTINE wl_rotate_mt(atoms, sym, jop, amap, l_anti, kimg, dwgn, ncomp, fspin, &
                           src, dst, dev_u, dev_a, gap)
      USE m_types_atoms; USE m_types_sym
      USE m_constants, ONLY: tpi_const, ImagUnit
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      INTEGER, INTENT(IN) :: jop, amap(:)
      LOGICAL, INTENT(IN) :: l_anti
      REAL, INTENT(IN) :: kimg(3)
      COMPLEX, INTENT(IN) :: dwgn(-atoms%lmaxd:, -atoms%lmaxd:, 0:, :)
      !> How many local spin directions travel together, and the 2x2 that mixes them for
      !> each type. They travel together because where the magnetisation is non-collinear
      !> they do not rotate independently: the coefficients of an atom live in that atom's
      !> own spin frame, so the rotation seen from inside the frame is the global one
      !> conjugated by the frame, and that mixes the directions. The factor is not guessed
      !> -- calc_abc builds a local direction as a^loc = A^T g with A = conjg(umat) of the
      !> type, so A^T U conjg(A) is what the rule has to carry, and A^T U A on the
      !> antiunitary half, where the conjugation of the state is already in the branch
      !> below. The identity for every type when nothing mixes.
      INTEGER, INTENT(IN) :: ncomp
      COMPLEX, INTENT(IN) :: fspin(:, :, :)
      TYPE(t_abc), INTENT(IN) :: src(:, :)
      TYPE(t_abc), INTENT(INOUT) :: dst(:, :)
      REAL, INTENT(INOUT) :: dev_u, dev_a, gap

      INTEGER :: n, na, nb, ia, ja, l, j, lm1, lm2, nord, m, is, is2, nuc, nlm, nu
      LOGICAL :: l_lo
      REAL :: glat(3)
      COMPLEX :: ph
      COMPLEX, ALLOCATABLE :: blk(:, :), tflip(:, :), rot(:, :, :)

      !> Sized once for the largest type and used in slices, rather than allocated inside
      !> the loop over l and the radial orders: the shapes only shrink from here.
      nu = 0
      DO n = 1, atoms%ntype
         nu = MAX(nu, SIZE(src(1, n)%cof, 1))
      END DO
      ALLOCATE (rot(nu, 2*atoms%lmaxd + 1, ncomp), blk(nu, 2*atoms%lmaxd + 1))
      ALLOCATE (tflip(-atoms%lmaxd:atoms%lmaxd, -atoms%lmaxd:atoms%lmaxd))

      DO n = 1, atoms%ntype
         nuc = SIZE(src(1, n)%cof, 1)
         DO na = 1, atoms%neq(n)
            ia = atoms%firstAtom(n) + na - 1
            ja = amap(ia)
            IF (ja < 1) CALL juDFT_error("wannierlib: an atom is mapped nowhere", &
                                         calledby="wl_rotate_mt")
            nb = ja - atoms%firstAtom(n) + 1
            !> A symmetry operation cannot take an atom out of its own type, and if the map
            !> said so the coefficients would be rotated against a different radial set.
            IF (nb < 1 .OR. nb > atoms%neq(n)) CALL juDFT_error( &
               "wannierlib: the atom map leaves the type", calledby="wl_rotate_mt")

            glat = MATMUL(REAL(sym%mrot(:, :, jop)), atoms%taual(:, ia)) &
                   + sym%tau(:, jop) - atoms%taual(:, ja)
            IF (ANY(ABS(glat - ANINT(glat)) > 1e-6)) CALL juDFT_error( &
               "wannierlib: the atom map is not a lattice jump", calledby="wl_rotate_mt")
            ph = EXP(-tpi_const*ImagUnit*DOT_PRODUCT(kimg, ANINT(glat)))

            DO l = 0, atoms%lmax(n)
               lm1 = l**2
               lm2 = l*(l + 2)
               nlm = lm2 - lm1 + 1
               nord = MIN(src(1, n)%n_r(l), SIZE(src(1, n)%cof, 3))
               !> Whether this l can be controlled at all. A local orbital does not merely
               !> add an order of its own: abclocdn feeds its coefficient back into the A
               !> and the B of the same l through alo1 and blo1, so where one sits the
               !> rebuilt orders are built from the zeroed rows as well and are no
               !> reference. There the rule takes the whole channel; everywhere else it is
               !> checked and nothing is written.
               l_lo = nord > 2
               IF (l_anti) THEN
                  !> The conjugation itself, in the lm basis. conj(Y_lm) = (-1)^m Y_l,-m,
                  !> so conjugating the state is not conjugating its coefficients: m flips
                  !> and the sign rides along.
                  tflip(-l:l, -l:l) = CMPLX(0.0, 0.0)
                  DO m = -l, l
                     tflip(-m, m) = CMPLX((-1.0)**m, 0.0)
                  END DO
               END IF
               DO j = 1, nord
                  DO is = 1, ncomp
                     IF (l_anti) THEN
                        !> The (-1)^l is the i**l that abcof leaves in these coefficients,
                        !> turned into (-i)**l by the conjugation -- the same factor
                        !> rot_to_unrotated carries as ifac**l. Measured on Pt, which has
                        !> time reversal: 2.4e-15 this way, and 1.703 for every one of the
                        !> seven other placements of the flip, of the two signs and of the
                        !> conjugate of D. The l = 0 block passes under all of them, which is
                        !> what says the atom map and the phase are not the question here.
                        rot(:nuc, :nlm, is) = (-1.0)**l &
                           *MATMUL(MATMUL(CONJG(src(is, n)%cof(:, lm1:lm2, j, na)), &
                                          CONJG(dwgn(-l:l, -l:l, l, jop))), tflip(-l:l, -l:l))
                     ELSE
                        !> d_wgn goes in as it is, neither transposed nor conjugated. Measured
                        !> against the rebuilt coefficients on wurtzite AlN, where a proper
                        !> operation has a symmetric block and only an improper one tells the
                        !> four placements apart: D passes everywhere, while the transpose and
                        !> the conjugate each fail on one half of the group. t_abc%rotate and
                        !> rot_to_unrotated conjugate because they go the other way, from the
                        !> local frame to the global one; and the parity of an improper
                        !> operation is already inside d_wgn, so it is not added here.
                        rot(:nuc, :nlm, is) = MATMUL(src(is, n)%cof(:, lm1:lm2, j, na), &
                                                     dwgn(-l:l, -l:l, l, jop))
                     END IF
                  END DO
                  DO is = 1, ncomp
                     blk(:nuc, :nlm) = CMPLX(0.0, 0.0)
                     DO is2 = 1, ncomp
                        blk(:nuc, :nlm) = blk(:nuc, :nlm) &
                                          + fspin(is, is2, n)*rot(:nuc, :nlm, is2)
                     END DO
                     blk(:nuc, :nlm) = ph*blk(:nuc, :nlm)
                     IF (l_lo) THEN
                        !> How far the rebuild was where the rule overrode it. Reported, not
                        !> asserted: it is the evidence that the override did something. A
                        !> number down at the rounding would say the local orbital rows never
                        !> reached these coefficients, and then this whole branch is a no-op
                        !> dressed as a fix.
                        gap = MAX(gap, MAXVAL(ABS(blk(:nuc, :nlm) &
                                                  - dst(is, n)%cof(:, lm1:lm2, j, nb))))
                        dst(is, n)%cof(:, lm1:lm2, j, nb) = blk(:nuc, :nlm)
                     ELSE IF (l_anti) THEN
                        !> The control: this channel's orders are also in dst, rebuilt from
                        !> the rotated eigenvector by the routine that owns that job.
                        dev_a = MAX(dev_a, MAXVAL(ABS(blk(:nuc, :nlm) &
                                                      - dst(is, n)%cof(:, lm1:lm2, j, nb))))
                     ELSE
                        dev_u = MAX(dev_u, MAXVAL(ABS(blk(:nuc, :nlm) &
                                                      - dst(is, n)%cof(:, lm1:lm2, j, nb))))
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE wl_rotate_mt

   !> The 2x2 the spin half of a symmetry operation contributes to the band table.
   !>
   !> It is the INVERSE rotation that belongs here, d^+, because the orbital half is the
   !> inverse too: d_wigner takes its angles from the inverse matrix while the euler of
   !> m_grp_k takes them from the direct one, and the two halves only agree once this one
   !> is adjointed. The two choices differ by r -> 1/r on the relative phase of the two
   !> components, which is REAL for alpha = 0 and pi, so the convention is invisible
   !> everywhere except on a rotation of 90 degrees -- and there it is worth 0.9 of
   !> unitarity. Measured on Pt with spin-orbit, where the four C4z and their four
   !> partners under time reversal were the only operations failing, and a sweep of r
   !> against the unitarity of the lowest Kramers doublet put the minimum at
   !> exp(-i(alpha+gamma)) in all sixteen.
   !>
   !> An antiunitary operation carries -i sigma_y on top of that, which SWAPS the two
   !> components: what it contributes is built from the OFF-DIAGONAL overlap blocks, and
   !> no phase on the diagonal pair can stand in for it.
   !>
   !> xas_su2_from_sym builds a 2x2 from the same euler angles and has done since June, but
   !> it is NOT this matrix: it is the one before the adjoint above, and it has no
   !> antiunitary branch. The two therefore cannot stand in for each other, and which of
   !> them a caller needs depends on how it pairs the spin half with d_wigner's orbital
   !> one. That pairing is what was measured here; what XAS needs is for XAS to say.
   SUBROUTINE wl_su2(sym, cell, jop, l_anti, u)
      USE m_types_sym; USE m_types_cell
      USE m_grp_k, ONLY: euler
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER, INTENT(IN) :: jop
      LOGICAL, INTENT(IN) :: l_anti
      COMPLEX, INTENT(OUT) :: u(2, 2)

      REAL :: alpha, beta, gamma
      COMPLEX :: d(2, 2), sy(2, 2)

      CALL euler(jop, sym, cell, alpha, beta, gamma)
      d(1, 1) = COS(beta/2.0)*EXP(CMPLX(0.0, -(alpha + gamma)/2.0))
      d(1, 2) = -SIN(beta/2.0)*EXP(CMPLX(0.0, -(alpha - gamma)/2.0))
      d(2, 1) = SIN(beta/2.0)*EXP(CMPLX(0.0, (alpha - gamma)/2.0))
      d(2, 2) = COS(beta/2.0)*EXP(CMPLX(0.0, (alpha + gamma)/2.0))
      u = TRANSPOSE(CONJG(d))

      IF (l_anti) THEN
         sy = CMPLX(0.0, 0.0)
         sy(1, 2) = CMPLX(-1.0, 0.0)
         sy(2, 1) = CMPLX(1.0, 0.0)
         u = MATMUL(u, sy)
      END IF
   END SUBROUTINE wl_su2

   !> The multiplication table of the exported group, in this export's own slot indices.
   !>
   !> prod(i,j) is the slot carrying the rotation and translation of i applied after j, or
   !> zero where that product left the group -- which it may, because only the operations
   !> FLEUR could map are exported and a subgroup of a group need not be closed once some of
   !> its members are struck out. Antiunitary composed with antiunitary is unitary, so the
   !> product is looked for among the slots of the resulting kind rather than all of them.
   !> What uses the table is the check that D(i.j) = D(i)D(j), which is the one statement
   !> about d_band that a wrong but self-consistent set of matrices cannot satisfy.
   !>
   !> sym%multab is FLEUR's own table and mapatom builds it outside the force branch, so it
   !> is always there. It still cannot answer for these slots: it is indexed by sym's
   !> numbering and knows nothing of the antiunitary half. It could serve as the lookup this
   !> loop does by hand, at the price of depending on its argument order -- {R_j}{R_i} and
   !> not the other way round -- and that trade has not been taken.
   SUBROUTINE wl_group_table(star, sym, prod)
      USE m_types_sym
      TYPE(t_wl_kstar), INTENT(IN) :: star
      TYPE(t_sym), INTENT(IN) :: sym
      INTEGER, ALLOCATABLE, INTENT(OUT) :: prod(:, :)

      INTEGER :: i, j, k, iop, jop, kop, mm(3, 3)
      REAL :: tt(3), dt(3)
      LOGICAL :: anti

      ALLOCATE (prod(star%nsym, star%nsym))
      prod = 0
      DO i = 1, star%nsym
         iop = star%op(i)
         DO j = 1, star%nsym
            jop = star%op(j)
            mm = MATMUL(sym%mrot(:, :, iop), sym%mrot(:, :, jop))
            tt = MATMUL(REAL(sym%mrot(:, :, iop)), sym%tau(:, jop)) + sym%tau(:, iop)
            anti = star%antiunit(i) .NEQV. star%antiunit(j)
            !> The product of two slots is a slot of the same kind, and with the two
            !> halves holding different rotations it has to be looked for among those.
            DO k = 1, star%nsym
               IF (star%antiunit(k) .NEQV. anti) CYCLE
               kop = star%op(k)
               IF (.NOT. ALL(sym%mrot(:, :, kop) == mm)) CYCLE
               dt = tt - sym%tau(:, kop)
               !> The translation is only defined modulo a lattice vector.
               IF (ANY(ABS(dt - ANINT(dt)) > 1e-6)) CYCLE
               prod(i, j) = k
               EXIT
            END DO
         END DO
      END DO
   END SUBROUTINE wl_group_table

   !> The properties that hold whatever the convention, so that a run says whether the
   !> table is worth writing without a reader on the other side: unitarity, that the
   !> identity came out the identity, and the group product table.
   SUBROUTINE wl_report_dband(star, sym, dband)
      USE m_types_sym
      TYPE(t_wl_kstar), INTENT(IN) :: star
      TYPE(t_sym), INTENT(IN) :: sym
      COMPLEX, INTENT(IN) :: dband(:, :, :, :)

      INTEGER :: ir, iop, nb, i, j, k, nc, ncut, npair
      INTEGER, ALLOCATABLE :: prod(:, :)
      REAL :: dev_u, dev_all, dev, dev_p, dev_id
      COMPLEX :: ph
      COMPLEX, ALLOCATABLE :: p(:, :)

      nb = SIZE(dband, 1)
      ALLOCATE (p(nb, nb))

      !> Over the WHOLE window this number is ~1 whenever the top multiplet is cut by the
      !> band window, and it says nothing about the table: a rotated state of the last band
      !> then has weight on a band that was not retained, so those rows cannot be unitary.
      !> What is worth reporting is the largest leading set that IS closed, which is also
      !> the set a consumer has to restrict to. Swept here rather than taken from the
      !> eigenvalues, so that the number is a property of this table alone.
      ncut = 0
      dev_all = 0.0
      DO nc = nb, 1, -1
         dev = 0.0
         DO ir = 1, star%nkirr
            DO iop = 1, star%nsym
               p(:nc, :nc) = MATMUL(dband(:nc, :nc, iop, ir), &
                                    TRANSPOSE(CONJG(dband(:nc, :nc, iop, ir))))
               DO i = 1, nc
                  p(i, i) = p(i, i) - CMPLX(1.0, 0.0)
               END DO
               dev = MAX(dev, MAXVAL(ABS(p(:nc, :nc))))
            END DO
         END DO
         IF (nc == nb) dev_all = dev
         IF (dev < 1e-6) THEN
            ncut = nc; dev_u = dev; EXIT
         END IF
      END DO

      dev_id = 0.0
      DO ir = 1, star%nkirr
         p = dband(:, :, 1, ir)
         DO i = 1, nb
            p(i, i) = p(i, i) - CMPLX(1.0, 0.0)
         END DO
         dev_id = MAX(dev_id, MAXVAL(ABS(p)))
      END DO

      !> D(i.j) = D(i) D(j) up to a global phase -- a non-symmorphic operation carries one,
      !> and it cancels out of everything a reader computes -- with the second factor
      !> conjugated when i is antiunitary, which is where the antilinearity lands. Pairs are
      !> taken inside the little group of the representative, the operations that leave it
      !> where it is, because only there is the product a matrix of the same star.
      !>
      !> This is the one number here that sees the (m,n) order of the table: unitarity, the
      !> identity and the block structure are all satisfied by its transpose. Reading d_band
      !> the other way round lands at 1e-01 and worse.
      npair = 0
      dev_p = 0.0
      IF (ncut > 0) THEN
         CALL wl_group_table(star, sym, prod)
         DO ir = 1, star%nkirr
            DO i = 1, star%nsym
               IF (star%kirr2k(ir, i) /= star%kirr(ir)) CYCLE
               DO j = 1, star%nsym
                  IF (star%kirr2k(ir, j) /= star%kirr(ir)) CYCLE
                  k = prod(i, j)
                  IF (k == 0) CYCLE
                  IF (star%kirr2k(ir, k) /= star%kirr(ir)) CYCLE
                  npair = npair + 1
                  IF (star%antiunit(i)) THEN
                     p(:ncut, :ncut) = MATMUL(dband(:ncut, :ncut, i, ir), &
                                              CONJG(dband(:ncut, :ncut, j, ir)))
                  ELSE
                     p(:ncut, :ncut) = MATMUL(dband(:ncut, :ncut, i, ir), &
                                              dband(:ncut, :ncut, j, ir))
                  END IF
                  ph = SUM(CONJG(dband(:ncut, :ncut, k, ir))*p(:ncut, :ncut))
                  IF (ABS(ph) > 1e-12) THEN
                     ph = ph/ABS(ph)
                  ELSE
                     ph = CMPLX(1.0, 0.0)
                  END IF
                  dev_p = MAX(dev_p, MAXVAL(ABS(p(:ncut, :ncut) &
                                                - ph*dband(:ncut, :ncut, k, ir))))
               END DO
            END DO
         END DO
         DEALLOCATE (prod)
      END IF

      WRITE (oUnit, '(/,a)') "Wannier90/WannierBerri symmetry export: d_band"
      WRITE (oUnit, '(a,i0,a,i0,a,i0,a)') "  table ", nb, "x", nb, " for ", &
         star%nkirr*star%nsym, " (star, operation) pairs"
      IF (ncut == 0) THEN
         WRITE (oUnit, '(a,es10.3)') "  NO leading band set is closed; | D D^+ - 1 | = ", dev_all
      ELSE
         WRITE (oUnit, '(a,i0,a,i0,a,es10.3)') "  closed up to band ", ncut, " of ", nb, &
            ", where | D D^+ - 1 | = ", dev_u
         IF (ncut < nb) WRITE (oUnit, '(a,es10.3,a)') &
            "  over the whole window it is ", dev_all, ", which is the cut and not the table"
      END IF
      WRITE (oUnit, '(a,es10.3)') "  | D(identity) - 1 | max ", dev_id
      WRITE (oUnit, '(a,i0)') "  little-group pairs checked = ", npair
      IF (npair > 0) WRITE (oUnit, '(a,es10.3)') &
         "  | D(i.j) - D(i)D(j) | max = ", dev_p
      DEALLOCATE (p)
   END SUBROUTINE wl_report_dband

END MODULE m_wannierlib_dband
