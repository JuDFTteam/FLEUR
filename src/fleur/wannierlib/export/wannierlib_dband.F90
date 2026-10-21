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
      TYPE(t_mat) :: z_irr, z_rot
      TYPE(t_mat), POINTER :: zmat_p(:)
      TYPE(t_abc), POINTER :: abc_p(:, :)
      TYPE(t_abc), ALLOCATABLE :: abc_rot(:)
      TYPE(t_spinor_layout) :: layout
      COMPLEX, ALLOCATABLE :: ujug0(:, :, :, :, :, :), ovl(:, :)
      INTEGER, ALLOCATABLE :: ev_list(:), perm(:)
      REAL :: kdiff0(3, 1), rrot(3, 3), kimg(3)
      INTEGER :: ir, iop, jop, kirr, jk, nb, ioff, ginv(3, 3), g0(3), j, jsp_rad, rsign
      LOGICAL :: l_anti

      !> Both limits are of the rotation, not of this table: see the module header.
      IF (atoms%nlotot > 0) CALL juDFT_error( &
         "wannierlib: d_band is not implemented for a basis with local orbitals", &
         hint="the symmetry export still writes the stars and the group without it", &
         calledby="wl_build_dband")
      IF (l_spinors) CALL juDFT_error( &
         "wannierlib: d_band is not implemented for spinor states", &
         hint="the spin half of the rotation is missing; a collinear run without SOC works", &
         calledby="wl_build_dband")

      nb = manifold%num_bands
      !> Which radial set this component reads: a single potential serves both.
      jsp_rad = radial_slot(radfun, jspin)
      ALLOCATE (dband(nb, nb, star%nsym, star%nkirr), source=CMPLX(0.0, 0.0))
      ALLOCATE (ovl(nb, nb))
      ev_list = [(j, j=manifold%min_band, manifold%max_band)]

      !> b = 0: the two sides sit at the same k, so the muffin-tin table is the one entry
      !> the neighbour list never holds, built the way the M(k,k) invariant builds it.
      kdiff0 = 0.0
      CALL melem_ujugaunt(atoms, cell, 1, kdiff0, radfun, radfun, jsp_rad, jsp_rad, .FALSE., 1, ujug0)

      CALL timestart("wl_build_dband")
      DO ir = 1, star%nkirr
         kirr = star%kirr(ir)
         CALL lapw_irr%init(input, noco, nococonv, kpts, atoms, sym, kirr, cell)
         CALL matrix_element_states(eig_id, kirr, input, atoms, sym, cell, noco, nococonv, &
                                    enpara, lapw_irr, vtot, fmpi, zmat_p, abc_p, &
                                    ev_list=ev_list, l_both_spinors=.FALSE., kpts=kpts)
         !> Copied, not pointed at: the slot it lives in is keyed on k and the image points
         !> below go through the same cache, which holds three.
         z_irr = zmat_p(MERGE(1, jspin, noco%l_noco))

         DO iop = 1, star%nsym
            jk = star%kirr2k(ir, iop)
            !> The image point's own basis, in FLEUR's own order. The rotated state is
            !> written into it, so nothing downstream sees a rotated anything.
            CALL lapw_img%init(input, noco, nococonv, kpts, atoms, sym, jk, cell)
            CALL matrix_element_states(eig_id, jk, input, atoms, sym, cell, noco, nococonv, &
                                       enpara, lapw_img, vtot, fmpi, zmat_p, abc_p, &
                                       ev_list=ev_list, l_both_spinors=.FALSE., kpts=kpts)
            CALL layout%init(input, noco, lapw_img, atoms)
            ioff = layout%row_offset(jspin)

            !> A k vector transforms under the inverse transpose, which is how the stars
            !> were built: operation iop here is the same operation iop there. Past nop
            !> the operation is that same rotation with time reversal, so the sign flips
            !> and the coefficients get conjugated below.
            l_anti = star%antiunit(iop)
            jop = MERGE(iop - star%nop, iop, l_anti)
            rsign = MERGE(-1, 1, l_anti)
            ginv = sym%mrot(:, :, sym%invtab(jop))
            rrot = rsign*TRANSPOSE(REAL(ginv))
            kimg = MATMUL(rrot, kpts%bkf(:, kirr))
            !> Rotating lands on the image point up to a reciprocal lattice vector. Folding
            !> it away here leaves both sides at one k and spares the overlap a gb. The
            !> label moves WITH g0, not against it: the physical vector k'+G' rewritten
            !> about bkf(jk) is bkf(jk) + (G' + g0). The other sign puts the vector 2*g0
            !> away and wl_gmap stops the run, which is how this was found.
            g0 = NINT(kimg - kpts%bkf(:, jk))

            CALL wl_gmap(lapw_irr, lapw_img, ginv, g0, rsign, perm)
            CALL z_rot%init(z_irr%l_real, z_irr%matsize1, z_irr%matsize2)
            CALL z_rot%clear()
            CALL wl_rotate_coeffs(lapw_irr, lapw_img, perm, sym%tau(:, jop), l_anti, &
                                  z_irr, z_rot)

            CALL matrix_element_abc_of(input, atoms, sym, cell, noco, nococonv, enpara, fmpi, &
                                       vtot, lapw_img, z_rot, nb, jsp_rad, abc_rot)

            ovl = CMPLX(0.0, 0.0)
            CALL melem_overlap_states(stars, atoms, lapw_img, lapw_img, &
                                      zmat_p(MERGE(1, jspin, noco%l_noco)), z_rot, &
                                      abc_p(jspin, :), abc_rot, jsp_rad, jsp_rad, &
                                      lapw_img%bkpt, lapw_img%bkpt, [0, 0, 0], &
                                      ujug0, kdiff0, 1, ioff, ioff, ovl)
            !> CONJUGATED. The overlap kernel is the mmn's, which puts its conjugation on
            !> the second argument, so what comes back is the conjugate of <a|R|b>. Which of
            !> the two the file should hold is not a matter of taste: the relation a reader
            !> applies is A(Rk) = d_band . A(k) . D_wann^+, and of the SIXTEEN placements of
            !> a conjugate, a transpose and an adjoint over the three matrices exactly one
            !> satisfies it -- measured on Si, 1.4e-12 against 1.4e+00 for the other fifteen.
            !> That one is this conjugate. Unitarity, the identity, the block structure and
            !> the group product table are all blind to it: the conjugate of a representation
            !> is a representation.
            dband(:, :, iop, ir) = CONJG(ovl)
            DEALLOCATE (abc_rot)
            CALL z_rot%free()
         END DO
      END DO
      CALL timestop("wl_build_dband")

      IF (fmpi%irank == 0) CALL wl_report_dband(star, sym, dband)
      DEALLOCATE (ujug0, ovl)
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
   SUBROUTINE wl_rotate_coeffs(lapw_a, lapw_b, perm, tau, l_anti, z_a, z_b)
      USE m_types_lapw; USE m_types_mat
      TYPE(t_lapw), INTENT(IN) :: lapw_a, lapw_b
      INTEGER, INTENT(IN) :: perm(:)
      REAL, INTENT(IN) :: tau(3)
      LOGICAL, INTENT(IN) :: l_anti
      TYPE(t_mat), INTENT(IN) :: z_a
      TYPE(t_mat), INTENT(INOUT) :: z_b

      INTEGER :: i, j
      REAL :: arg
      COMPLEX :: ph

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
            z_b%data_r(j, :) = z_a%data_r(i, :)     ! conjugating real data is a no-op
         ELSE
            ph = EXP(ImagUnit*arg)
            IF (l_anti) THEN
               z_b%data_c(j, :) = ph*CONJG(z_a%data_c(i, :))
            ELSE
               z_b%data_c(j, :) = ph*z_a%data_c(i, :)
            END IF
         END IF
      END DO
   END SUBROUTINE wl_rotate_coeffs

   !> The two properties that hold whatever the convention, so that a run says whether the
   !> table is worth writing without a reader on the other side: unitarity, and that the
   !> identity came out the identity.
   !> prod(i, j) = the operation that is i composed with j, or 0 when it is not in the set.
   !>
   !> An element is (R, t, a) with a the time-reversal flag, and
   !> (R1,t1,a1).(R2,t2,a2) = (R1 R2, R1 t2 + t1, a1 xor a2): T commutes with the spatial
   !> part and does not act on the coordinates. The flag takes part in the match, so a
   !> unitary operation never pairs with its antiunitary partner. FLEUR keeps no
   !> multiplication table, and a search costs nothing at these sizes.
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
         iop = MERGE(i - star%nop, i, star%antiunit(i))
         DO j = 1, star%nsym
            jop = MERGE(j - star%nop, j, star%antiunit(j))
            mm = MATMUL(sym%mrot(:, :, iop), sym%mrot(:, :, jop))
            tt = MATMUL(REAL(sym%mrot(:, :, iop)), sym%tau(:, jop)) + sym%tau(:, iop)
            anti = star%antiunit(i) .NEQV. star%antiunit(j)
            DO kop = 1, star%nop
               IF (.NOT. ALL(sym%mrot(:, :, kop) == mm)) CYCLE
               dt = tt - sym%tau(:, kop)
               !> The translation is only defined modulo a lattice vector.
               IF (ANY(ABS(dt - ANINT(dt)) > 1e-6)) CYCLE
               prod(i, j) = MERGE(kop + star%nop, kop, anti)
               EXIT
            END DO
         END DO
      END DO
   END SUBROUTINE wl_group_table

   SUBROUTINE wl_report_dband(star, sym, dband)
      USE m_types_sym
      TYPE(t_wl_kstar), INTENT(IN) :: star
      TYPE(t_sym), INTENT(IN) :: sym
      COMPLEX, INTENT(IN) :: dband(:, :, :, :)

      INTEGER :: ir, iop, nb, i, j, k, nc, ncut, npair
      INTEGER, ALLOCATABLE :: prod(:, :)
      REAL :: dev_u, dev_1, dev, dev_p
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
         IF (dev < 1e-6) THEN
            ncut = nc; dev_u = dev; EXIT
         END IF
         IF (nc == nb) dev_1 = dev        ! guardado para el informe de la ventana entera
      END DO

      dev = dev_1
      dev_1 = 0.0
      DO ir = 1, star%nkirr
         p = dband(:, :, 1, ir)
         DO i = 1, nb
            p(i, i) = p(i, i) - CMPLX(1.0, 0.0)
         END DO
         dev_1 = MAX(dev_1, MAXVAL(ABS(p)))
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
         WRITE (oUnit, '(a,es10.3)') "  NO leading band set is closed; | D D^+ - 1 | = ", dev
      ELSE
         WRITE (oUnit, '(a,i0,a,i0,a,es10.3)') "  closed up to band ", ncut, " of ", nb, &
            ", where | D D^+ - 1 | = ", dev_u
         IF (ncut < nb) WRITE (oUnit, '(a,es10.3,a)') &
            "  over the whole window it is ", dev, ", which is the cut and not the table"
      END IF
      WRITE (oUnit, '(a,es10.3)') "  | D(identity) - 1 | max ", dev_1
      WRITE (oUnit, '(a,i0)') "  little-group pairs checked = ", npair
      IF (npair > 0) WRITE (oUnit, '(a,es10.3)') &
         "  | D(i.j) - D(i)D(j) | max = ", dev_p
      DEALLOCATE (p)
   END SUBROUTINE wl_report_dband

END MODULE m_wannierlib_dband
