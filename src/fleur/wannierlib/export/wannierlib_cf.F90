!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The two pair contractions of the modern theory, which differ only in what sits between
!>  the neighbours:
!>
!>      F_ab(k) = sum_{b1,b2} w_b1 b1_a w_b2 b2_b <u_{k+b1}      |u_{k+b2}>
!>      C_ab(k) = sum_{b1,b2} w_b1 b1_a w_b2 b2_b <u_{k+b1}|H_k  |u_{k+b2}>
!>
!>  Both are Eq. (84) of Lopez, Vanderbilt, Thonhauser and Souza, PRB 85, 014435 (2012),
!>  with the finite-difference derivative of their Eq. (85) on each side; F is the same
!>  contraction with the Hamiltonian left out, and is the quantum geometric tensor.
!>
!>  The Hamiltonian is never applied. H_k differs from the one the bra is an eigenstate of
!>  by an operator linear in the momentum,
!>
!>      H_k - H_{k+b1} = -b1.(p + k) - b1^2/2 ,
!>
!>  so the matrix element follows from the eigenvalue at k+b1, the pair overlap and the
!>  momentum between the two neighbours,
!>
!>      <u_{k+b1}|H_k|u_{k+b2}> = E_{k+b1} M - b1.Pi - (b1^2/2) M ,
!>
!>  which is what lets C be built without ever writing a uHu file.
!>
!>  TWO SETS OF b VECTORS APPEAR BELOW AND MIXING THEM IS SILENT. The identity comes from
!>  H_q = (p+q)^2/2 + V in Hartree, so the b multiplying Pi is in atomic units and goes
!>  through cell%bmat. The b carrying the weights of Eq. (84) is the one Wannier90 supplies,
!>  in inverse Angstrom, and is what puts the result in the units its consumer expects.
!>
!>  Both need the wavefunctions, so neither can run in the post-processing, and both need
!>  the gauge at two different k-points, so neither can run with the coarse matrices. They
!>  therefore sit between the two, and what they pass on is the contracted
!>  (num_wann, num_wann, 3, 3) block per k rather than the pair overlaps, which are far
!>  larger than the states they came from.
MODULE m_wannierlib_cf
  USE m_juDFT
  USE m_melem_overlap, ONLY: melem_overlap_states
  USE m_melem_ujugaunt, ONLY: melem_ujugaunt
  USE m_melem_nablaujugaunt, ONLY: melem_nablaujugaunt
  USE m_melem_nabla_int, ONLY: melem_nabla_int
  USE m_melem_nabla_sph, ONLY: melem_nabla_sph
  USE m_melem_nabla_vac, ONLY: melem_nabla_vac
  USE m_melem_uhu, ONLY: melem_uhu_states
  USE m_matrix_element_factory, ONLY: matrix_element_states, matrix_element_release_anchor
  USE m_eig66_io, ONLY: read_eig
  USE m_constants, ONLY: oUnit
  USE m_types
  USE m_types_abc
  USE m_types_radfun
  USE m_types_spinor_layout, ONLY: t_spinor_layout, radial_slot
  USE m_types_atoms
  USE m_types_cell
  USE m_types_input
  USE m_types_kpts
  USE m_types_noco
  USE m_types_nococonv
  USE m_types_sym
  USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
  USE m_types_wgauge_bmesh, ONLY: t_wgauge_bmesh
  USE m_types_enpara
  USE m_types_potden
  USE m_types_vacuum
  USE m_types_melem_vacabc, ONLY: t_melem_vacabc
  USE m_types_mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_uiu, wannierlib_uhu
  !> C at (1,b2) and at (b2,1), summed over the spin components, at one k-point. One row
  !> and one column of the pair table, which is what the adjoint check needs.
  COMPLEX, ALLOCATABLE :: pair_row(:, :, :), pair_col(:, :, :)

CONTAINS


  !> F on this rank's k-slice, in the Wannier gauge.
  SUBROUTINE wannierlib_uiu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                            radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, &
                            vacuum, fmpi, distk, u_matrix, u_opt, f0)
    TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
    TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    TYPE(t_radfun), INTENT(IN) :: radfun(:)
    INTEGER, INTENT(IN) :: jspin          !> the channel being wannierised
    LOGICAL, INTENT(IN) :: l_spinors      !> the states carry both spin components
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_vacuum), INTENT(IN) :: vacuum
    TYPE(t_mpi), INTENT(IN) :: fmpi
    INTEGER, INTENT(IN) :: distk(:)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)   !> (nw,nw,nk)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)      !> (nb,nw,nk)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: f0(:, :, :, :, :)  !> (nw,nw,3,3,nk_loc)

    CALL cf_contract(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, radfun, &
                     jspin, l_spinors, eig_id, stars, enpara, vtot, vacuum, fmpi, distk, &
                     u_matrix, u_opt, .FALSE., f0)
  END SUBROUTINE wannierlib_uiu

  !> C on this rank's k-slice, in the Wannier gauge.
  SUBROUTINE wannierlib_uhu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                            radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, &
                            vacuum, fmpi, distk, u_matrix, u_opt, c0)
    TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
    TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    TYPE(t_radfun), INTENT(IN) :: radfun(:)
    INTEGER, INTENT(IN) :: jspin
    LOGICAL, INTENT(IN) :: l_spinors
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_vacuum), INTENT(IN) :: vacuum
    TYPE(t_mpi), INTENT(IN) :: fmpi
    INTEGER, INTENT(IN) :: distk(:)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: c0(:, :, :, :, :)

    CALL cf_contract(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, radfun, &
                     jspin, l_spinors, eig_id, stars, enpara, vtot, vacuum, fmpi, distk, &
                     u_matrix, u_opt, .TRUE., c0)
  END SUBROUTINE wannierlib_uhu

  !> The contraction both share. l_ham decides whether the Hamiltonian sits between the
  !> neighbours, which is the only difference between C and F.
  !>
  !> One radial table per radial set, shared by every pair: the muffin-tin half of an overlap
  !> depends on the two sides only through b2 - b1, and there are far fewer distinct
  !> differences than pairs. Both spin components accumulate into the same result when the
  !> states are spinors.
  SUBROUTINE cf_contract(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                         radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, vacuum, &
                         fmpi, distk, u_matrix, u_opt, l_ham, o0)
    TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
    TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    TYPE(t_radfun), INTENT(IN) :: radfun(:)
    INTEGER, INTENT(IN) :: jspin
    LOGICAL, INTENT(IN) :: l_spinors
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_vacuum), INTENT(IN) :: vacuum
    TYPE(t_mpi), INTENT(IN) :: fmpi
    INTEGER, INTENT(IN) :: distk(:)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)
    LOGICAL, INTENT(IN) :: l_ham          !> C when true, F when false
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: o0(:, :, :, :, :)

    COMPLEX, ALLOCATABLE :: vgauge(:, :, :), ujp(:, :, :, :, :, :), nujp(:, :, :, :, :, :, :)
    REAL, ALLOCATABLE :: kdp(:, :)
    INTEGER :: npair, nb, nw, nkl, k, kl, jcomp, jrad, lo, hi, jeig, kdiag, i
    LOGICAL :: l_uhu, l_here
    REAL :: diag(26)
    REAL :: dpair, nrmp, wa
    INTEGER :: ia, ic, ik

    CALL timestart(MERGE("wannierlib_uhu", "wannierlib_uiu", l_ham))
    nb = manifold%num_bands
    nw = manifold%num_wann
    nkl = COUNT(distk == fmpi%irank)
    ALLOCATE (o0(nw, nw, 3, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    IF (nb <= 0 .OR. nw <= 0 .OR. bmesh%nntot <= 0 .OR. nkl <= 0) THEN
      CALL timestop(MERGE("wannierlib_uhu", "wannierlib_uiu", l_ham))
      RETURN
    END IF
    IF (.NOT. ALLOCATED(bmesh%wb) .OR. .NOT. ALLOCATED(bmesh%bk)) CALL juDFT_error( &
      'wannierlib: '//MERGE('C', 'F', l_ham)//' needs the b-shell weights, which the '// &
      'wannierisation produces', calledby='cf_contract')

    !> Off unless asked for. The direct uHu costs orders of magnitude more than the
    !> identity it checks, so it runs at a single k-point and reports; nothing downstream
    !> changes, which is what makes the two coexist in one run.
    l_uhu = l_ham .AND. juDFT_was_argument("-uhu_direct")
    diag = 0.0
    kdiag = 0

    ALLOCATE (vgauge(nb, nw, kpts%nkptf))
    DO k = 1, kpts%nkptf
      vgauge(:, :, k) = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
    END DO

    CALL bmesh%pair_diffs(kpts%bkf, kdp, npair)

    IF (l_ham .AND. fmpi%irank == 0) &
      ALLOCATE (pair_row(nb, nb, bmesh%nntot), pair_col(nb, nb, bmesh%nntot), &
                source=CMPLX(0.0, 0.0))

    lo = MERGE(1, jspin, l_spinors)
    hi = MERGE(2, jspin, l_spinors)
    jeig = MERGE(1, jspin, l_spinors)
    DO jcomp = lo, hi
      jrad = radial_slot(radfun, jcomp)
      CALL melem_ujugaunt(atoms, cell, npair, kdp, radfun, radfun, jrad, jrad, .FALSE., 1, ujp)
      !> Allocated either way so the core takes it unconditionally; only C fills it.
      IF (l_ham) THEN
        CALL melem_nablaujugaunt(atoms, cell, enpara, npair, kdp, radfun, radfun, jrad, jrad, nujp)
      ELSE
        ALLOCATE (nujp(1, 1, 1, 1, 1, 1, 1), source=CMPLX(0.0, 0.0))
      END IF
      kl = 0
      DO k = 1, kpts%nkptf
        IF (distk(k) /= fmpi%irank) CYCLE
        kl = kl + 1
        l_here = l_uhu .AND. fmpi%irank == 0 .AND. kl == 1
        IF (l_here) kdiag = k
        CALL cf_one_k(manifold, bmesh, k, kl, kpts, ujp, nujp, kdp, npair, atoms, cell, &
                      input, sym, noco, nococonv, jcomp, jrad, jeig, eig_id, stars, &
                      enpara, vtot, vacuum, fmpi, l_ham, vgauge, radfun, l_here, diag, o0)
      END DO
      DEALLOCATE (ujp, nujp)
    END DO

    !> The pair loop keeps one neighbour anchored while it walks the others; nothing after
    !> this holds a state by pointer, so the slot goes back to the pool.
    CALL matrix_element_release_anchor()
    IF (l_uhu .AND. fmpi%irank == 0 .AND. diag(2) > 0.0) CALL cf_report(kdiag, diag)
    !> The hermiticity of the momentum costs nothing and is measured whenever there is a
    !> momentum: it is the only exact test of the sign of the vacuum half, and keeping quiet
    !> about it would be a return to the state where a wrong sign went unreported.
    IF (l_ham .AND. fmpi%irank == 0 .AND. MAXVAL(diag(20:22)) > 0.0) THEN
      WRITE (oUnit, '(/,a)') 'wannierlib: the momentum at b1 = b2, where it must be hermitian'
      WRITE (oUnit, '(a,3es12.4)') '   interstitial only,     x y z  ', &
         (diag(16 + i)/MAX(diag(19 + i), 1.0e-30), i=1, 3)
      WRITE (oUnit, '(a,3es12.4)') '   MT+int only,           x y z  ', &
         (diag(10 + i)/MAX(diag(13 + i), 1.0e-30), i=1, 3)
      WRITE (oUnit, '(a,3es12.4)') '   with the vacuum,       x y z  ', &
         (diag(4 + i)/MAX(diag(7 + i), 1.0e-30), i=1, 3)
      WRITE (oUnit, '(a)') '   All three regions at once. The vacuum half alone is NOT hermitian'
      WRITE (oUnit, '(a)') '   -- restricting p to a region leaves a boundary term -- so this'
      WRITE (oUnit, '(a)') '   only closes if the vacuum enters with the right sign. x and y'
      WRITE (oUnit, '(a)') '   test the ket wavevector, z tests the derivative and its sign.'
    END IF
    IF (ALLOCATED(pair_row)) THEN
      !> C(1,b) and C(b,1) are the same matrix element read the other way round, so one has
      !> to be the adjoint of the other. The pairs with b1 /= b2 are what C is built from,
      !> and this is what constrains them.
      !>
      !> IT DOES NOT CLOSE, AND THAT IS NOT A FAULT. Resolving H against the eigenvalue of
      !> the bra is exact only for an exact eigenfunction, and a variational solution is one
      !> only inside its own basis. Writing 1 = P + Q for the basis at k+b1, what is left
      !> over is <u(k+b1)| H Q |u(k+b2)>, the Hamiltonian reaching the part of the ket
      !> outside the bra's basis. It falls as the basis grows, and at b1 = b2 the two bases
      !> coincide, Q|u> = 0, and the term vanishes identically.
      !>
      !> What it CANNOT see: a correction whose coefficient is symmetric in b1 and b2. The
      !> gauge term enters C as (b1.b2) M, and a symmetric scalar cancels between the two
      !> orders, so removing it altogether leaves this number unchanged.
      WRITE (oUnit, '(/,a)') 'wannierlib: C at (1,b2) against C at (b2,1), one the adjoint of the other'
      WRITE (oUnit, '(a)')   '        b2     |dC|          |C|        relative'
      dpair = 0.0
      DO i = 1, bmesh%nntot
        nrmp = MAXVAL(ABS(pair_row(:, :, i)))
        wa = MAXVAL(ABS(pair_row(:, :, i) - CONJG(TRANSPOSE(pair_col(:, :, i)))))
        WRITE (oUnit, '(i10,3es14.5)') i, wa, nrmp, wa/MAX(nrmp, 1.0e-30)
        dpair = MAX(dpair, wa/MAX(nrmp, 1.0e-30))
      END DO
      !> the worst PER PAIR, not max|dC|/max|C|: taking the two maxima separately divides a
      !> large numerator from one pair by a large denominator from another and understates
      WRITE (oUnit, '(a,es12.4)') '   worst per pair                  ', dpair
      DEALLOCATE (pair_row, pair_col)
    END IF
    !> HOW MUCH OF THE NON-ADJOINT PART OF C REACHES O. The departure above is in matrix
    !> elements; what matters is what survives the weighted sum over pairs and bands, and
    !> the two differ by orders of magnitude.
    !>
    !> Symmetrising C over the pair swap -- ham(b1,b2) -> [ham(b1,b2) + ham(b2,b1)+]/2 --
    !> is EXACTLY the same as O(a,c) -> [O(a,c) + O(c,a)+]/2, because relabelling b1 <-> b2
    !> in the weighted sum swaps the cartesian indices. So this costs no storage: it is the
    !> departure of the assembled tensor from that symmetry.
    !>
    !> C IS NOT SYMMETRISED. The departure is reported and nothing else: it is the only
    !> observable handle on the error above, and imposing the symmetry would remove it while
    !> leaving the part that shares it untouched.
    IF (l_ham .AND. fmpi%irank == 0 .AND. SIZE(o0, 5) > 0) THEN
      wa = 0.0; nrmp = 0.0
      DO ik = 1, SIZE(o0, 5)
        DO ia = 1, 3
          DO ic = 1, 3
            wa = MAX(wa, MAXVAL(ABS(o0(:, :, ia, ic, ik) &
                                    - CONJG(TRANSPOSE(o0(:, :, ic, ia, ik))))))
            nrmp = MAX(nrmp, MAXVAL(ABS(o0(:, :, ia, ic, ik))))
          END DO
        END DO
      END DO
      WRITE (oUnit, '(/,a)') 'wannierlib: how much of the non-adjoint part of C reaches O'
      WRITE (oUnit, '(a,2es13.5)') '   max|O(a,c) - O(c,a)+|   max|O|   ', wa, nrmp
      WRITE (oUnit, '(a,es13.5)')  '   relative                        ', wa/MAX(nrmp, 1.0e-30)
      WRITE (oUnit, '(a)') '   C itself is left alone. This is what imposing the symmetry'
      WRITE (oUnit, '(a)') '   would change, and so the consequence of the departure above.'
    END IF
    DEALLOCATE (vgauge, kdp)
    CALL timestop(MERGE("wannierlib_uhu", "wannierlib_uiu", l_ham))
  END SUBROUTINE cf_contract

  !> One k-point's contribution, accumulated into its slice of o0.
  !>
  !> Every ordered pair of neighbours is visited, the diagonal b1 = b2 included: there the
  !> overlap is the identity in the ab initio gauge and the term reduces to the plain gauge
  !> product, which is what makes it worth computing rather than assuming.
  !>
  !> The overlap and the momentum arrive conjugated, and both are undone here, per pair:
  !> the identity below mixes them with real eigenvalues and real b components, so it is
  !> only meaningful once all three are in the same convention.
  !>
  !> The momentum the nabla modules return is that of the full wavefunction, p_psi. The
  !> identity asks for p_u + k on the periodic part, and since psi = e^{i(k+b2).r} u the
  !> two differ by exactly -b2, which is subtracted here.
  !>
  !> b1 stays anchored in the factory while its partners are fetched, or the states it holds
  !> by pointer would be overwritten by the third k asked for.
  SUBROUTINE cf_one_k(manifold, bmesh, nk, nk_local, kpts, ujug_pair, nujug_pair, &
                      kdiff_pair, npair, atoms, cell, input, sym, noco, nococonv, jspin, &
                      jspin_rad, jeig, eig_id, stars, enpara, vtot, vacuum, fmpi, l_ham, &
                      vgauge, radfun, l_diag, diag, o0)
    TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
    TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
    INTEGER, INTENT(IN) :: nk, nk_local
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: ujug_pair(:, :, :, :, :, :)
    COMPLEX, INTENT(IN) :: nujug_pair(:, :, :, :, :, :, :)
    REAL, INTENT(IN) :: kdiff_pair(:, :)
    INTEGER, INTENT(IN) :: npair
    TYPE(t_atoms), INTENT(IN) :: atoms
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_input), INTENT(IN) :: input
    TYPE(t_sym), INTENT(IN) :: sym
    TYPE(t_noco), INTENT(IN) :: noco
    TYPE(t_nococonv), INTENT(IN) :: nococonv
    INTEGER, INTENT(IN) :: jspin       !> the eig record
    INTEGER, INTENT(IN) :: jspin_rad   !> the radial index
    INTEGER, INTENT(IN) :: jeig        !> the eigenvalue record
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_vacuum), INTENT(IN) :: vacuum
    TYPE(t_mpi), INTENT(IN) :: fmpi
    LOGICAL, INTENT(IN) :: l_ham
    COMPLEX, INTENT(IN) :: vgauge(:, :, :)
    TYPE(t_radfun), INTENT(IN) :: radfun(:)
    LOGICAL, INTENT(IN) :: l_diag     !> also build C the direct way, at this k
    REAL, INTENT(INOUT) :: diag(:)    !> the two routes projected onto each other
    COMPLEX, INTENT(INOUT) :: o0(:, :, :, :, :)

    TYPE(t_mat), POINTER :: zmat_a(:), zmat_b(:)
    TYPE(t_abc), POINTER :: abc_a(:, :), abc_b(:, :)
    TYPE(t_lapw) :: lapw_a, lapw_b
    TYPE(t_spinor_layout) :: layout_a, layout_b
    !> The vacuum half of the pair overlap. Empty (nv2 = 0) when there is no vacuum to
    !> expand, and melem_overlap_states skips it on that alone, so they are handed over
    !> unconditionally and the caller never branches on the geometry.
    TYPE(t_melem_vacabc) :: vac_a, vac_b
    LOGICAL :: l_vac
    COMPLEX, ALLOCATABLE :: ovl(:, :), pnk(:, :, :), ham(:, :), tmp(:, :), mw(:, :)
    COMPLEX, ALLOCATABLE :: hud(:, :)
    REAL, ALLOCATABLE :: ea(:)
    REAL :: asym
    INTEGER, ALLOCATABLE :: ev_list(:)
    INTEGER :: b1, b2, ka, kb, gpar(3), irec, nb, nw, i, a, c
    REAL :: b1c(3), b2c(3), wa, wc

    nb = manifold%num_bands
    nw = manifold%num_wann
    ev_list = [(i, i = manifold%min_band, manifold%max_band)]
    irec = MERGE(1, jspin, noco%l_noco)
    ALLOCATE (ovl(nb, nb), tmp(nb, nw), mw(nw, nw))
    !> Both halves or neither. C is built from the identity
    !> <u_{k+b1}|H_k|u_{k+b2}> = E.M - b1.P - (b1^2/2).M, which replaces H by the eigenvalue
    !> of the bra and therefore holds only for the integral over the WHOLE cell: it does not
    !> split by region. An M that carried the vacuum next to a P that did not would be wrong
    !> by exactly b1.P_vac, and silently. Now that melem_nabla_vac exists the two go together,
    !> and the one flag is what keeps them that way.
    l_vac = input%film
    IF (l_ham) ALLOCATE (pnk(nb, nb, 3), ham(nb, nb), ea(input%neig))
    IF (l_diag) ALLOCATE (hud(nb, nb))

    DO b1 = 1, bmesh%nntot
      ka = bmesh%nnlist(nk, b1)
      CALL lapw_a%init(input, noco, nococonv, kpts, atoms, sym, ka, cell)
      CALL matrix_element_states(eig_id, ka, input, atoms, sym, cell, noco, nococonv, &
                                 enpara, lapw_a, vtot, fmpi, zmat_a, abc_a, ev_list=ev_list, &
                                 l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                 kpts=kpts, l_anchor=.TRUE.)
      CALL layout_a%init(input, noco, lapw_a, atoms)
      IF (l_vac) CALL vac_a%calc(vacuum, cell, enpara, vtot, lapw_a, jspin_rad, &
                                 zmat_a(irec), nb, ioff=layout_a%row_offset(jspin))

      IF (l_ham) THEN
        !> The eigenvalues of the bra, which is the side the Hamiltonian is resolved on.
        !> Read from the eigenvector file at the parent point: results%eig is not filled
        !> per k at the moment this runs.
        ea = 0.0
        CALL read_eig(eig_id, kpts%bkp(ka), jeig, eig=ea)
        !> In atomic units, because the identity is. NOT bmesh%bk, which is Wannier90's
        !> and in inverse Angstrom.
        b1c = MATMUL(bmesh%shell_vector(kpts%bkf, nk, b1), cell%bmat)
      END IF

      DO b2 = 1, bmesh%nntot
        kb = bmesh%nnlist(nk, b2)
        !> Each side was folded back into the mesh by a reciprocal lattice vector of its own,
        !> and only their difference acts on the stored states.
        gpar = bmesh%gkpb(:, nk, b2) - bmesh%gkpb(:, nk, b1)
        CALL lapw_b%init(input, noco, nococonv, kpts, atoms, sym, kb, cell)
        CALL matrix_element_states(eig_id, kb, input, atoms, sym, cell, noco, nococonv, &
                                   enpara, lapw_b, vtot, fmpi, zmat_b, abc_b, ev_list=ev_list, &
                                   l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                   kpts=kpts)
        CALL layout_b%init(input, noco, lapw_b, atoms)
        IF (l_vac) CALL vac_b%calc(vacuum, cell, enpara, vtot, lapw_b, jspin_rad, &
                                   zmat_b(irec), nb, ioff=layout_b%row_offset(jspin))

        IF (l_ham) b2c = MATMUL(bmesh%shell_vector(kpts%bkf, nk, b2), cell%bmat)

        ovl = CMPLX(0.0, 0.0)
        CALL melem_overlap_states(stars, atoms, lapw_a, lapw_b, zmat_a(irec), zmat_b(irec), &
                                  abc_a(jspin, :), abc_b(jspin, :), jspin_rad, jspin_rad, &
                                  kpts%bkf(:, ka), kpts%bkf(:, kb), gpar, &
                                  ujug_pair, kdiff_pair, npair, &
                                  ioff_a=layout_a%row_offset(jspin), &
                                  ioff_b=layout_b%row_offset(jspin), ovl=ovl, &
                                  vac_a=vac_a, vac_b=vac_b)

        IF (l_ham) THEN
          !> Conjugated here and folded into the MATMUL on the other branch: the two are
          !> not interchangeable, the sum order differs and the last digits with it.
          ovl = CONJG(ovl)
          pnk = CMPLX(0.0, 0.0)
          CALL melem_nabla_int(stars, cell, lapw_a, lapw_b, jspin_rad, jspin_rad, &
                               zmat_a(irec), zmat_b(irec), kpts%bkf(:, kb), gpar, pnk, &
                               ioff=layout_a%row_offset(jspin), &
                               ioff_b=layout_b%row_offset(jspin))
          !> The interstitial alone, before the muffin tin joins it.
          IF (b1 == b2) THEN
            DO a = 1, 3
              diag(16 + a) = MAX(diag(16 + a), &
                                 MAXVAL(ABS(pnk(:, :, a) - CONJG(TRANSPOSE(pnk(:, :, a))))))
              diag(19 + a) = MAX(diag(19 + a), MAXVAL(ABS(pnk(:, :, a))))
            END DO
          END IF

          CALL melem_nabla_sph(atoms, sym, cell, abc_a(jspin, :), abc_b(jspin, :), &
                               kpts%bkf(:, kb), gpar, kpts%bkf(:, ka), &
                               nujug_pair, kdiff_pair, npair, pnk)
          !> MT + interstitial. Second of the three snapshots, and all three are taken HERE,
          !> before the CONJG and before the -b2c*ovl correction: taking them at different
          !> points of the pipeline made them disagree even in bulk, where there is no vacuum
          !> to add.
          IF (b1 == b2) THEN
            DO a = 1, 3
              diag(10 + a) = MAX(diag(10 + a), &
                                 MAXVAL(ABS(pnk(:, :, a) - CONJG(TRANSPOSE(pnk(:, :, a))))))
              diag(13 + a) = MAX(diag(13 + a), MAXVAL(ABS(pnk(:, :, a))))
            END DO
          END IF

          !> Before the CONJG, like the two halves above it: melem_nabla_vac is written in
          !> the same unconjugated-bra convention they are.
          IF (l_vac) CALL melem_nabla_vac(vac_a, vac_b, gpar, pnk)

          !> And all three regions. What comparing the three rows separates: the interstitial
          !> is anti-hermitian on its own -- that part of it is a surface term -- and in a cell
          !> periodic in x and y the only surfaces it sees are the spheres, so in the plane the
          !> muffin tin has to cancel it, vacuum or no vacuum. If x and y do not close already
          !> in the second row, what is missing is not in the third region.
          IF (b1 == b2) THEN
            DO a = 1, 3
              diag(4 + a) = MAX(diag(4 + a), &
                                MAXVAL(ABS(pnk(:, :, a) - CONJG(TRANSPOSE(pnk(:, :, a))))))
              diag(7 + a) = MAX(diag(7 + a), MAXVAL(ABS(pnk(:, :, a))))
            END DO
          END IF
          pnk = CONJG(pnk)

          !> p_psi -> p_u + k. One term, and it depends on the pair.
          DO a = 1, 3
            pnk(:, :, a) = pnk(:, :, a) - CMPLX(b2c(a), 0.0)*ovl
          END DO

          !> The identity. The eigenvalue multiplies from the left because the bra is what
          !> it belongs to; the quadratic term is a scalar and does not care.
          ham = -CMPLX(0.5*DOT_PRODUCT(b1c, b1c), 0.0)*ovl
          DO a = 1, 3
            ham = ham - CMPLX(b1c(a), 0.0)*pnk(:, :, a)
          END DO
          DO i = 1, nb
            ham(i, :) = ham(i, :) + ea(manifold%min_band + i - 1)*ovl(i, :)
          END DO

          !> The same element by the other route. Conjugated to reach the convention ham is
          !> in, because it comes back the way melem_overlap_states hands an overlap over
          !> and that one is conjugated above. Muffin-tin sphere and spherical potential
          !> only, so what this settles is the sign, not the digits.
          IF (l_diag) THEN
            hud = CMPLX(0.0, 0.0)
            CALL melem_uhu_states(atoms, cell, enpara, radfun, jspin_rad, &
                                  abc_a(jspin, :), abc_b(jspin, :), &
                                  bmesh%shell_vector(kpts%bkf, nk, b1), &
                                  bmesh%shell_vector(kpts%bkf, nk, b2), &
                                  kpts%bkf(:, ka), kpts%bkf(:, kb), gpar, asym, hud)
            hud = CONJG(hud)
            diag(1) = diag(1) + REAL(SUM(CONJG(ham)*hud))
            diag(2) = diag(2) + SUM(ABS(ham)**2)
            diag(3) = diag(3) + SUM(ABS(hud)**2)
            diag(4) = MAX(diag(4), asym)
            !> The DIAGONAL pair kept apart from the rest.
            !>
            !> At b1 = b2 the identity is EXACT -- bra and ket come from the same LAPW basis,
            !> so resolving H against the bra is legitimate and the residue is 1e-16. The
            !> direct uHu therefore has a machine-precision target there, one that assumes
            !> nothing about b1 /= b2, which is the case under question. What is left over is
            !> exactly what the direct route is still missing: the interstitial, the
            !> non-spherical potential and the local orbitals.
            !>
            !> Mixed over all pairs, as diag(1..3) does it, the number cannot be read: the
            !> 0.863 measured so far averages a case with a known answer together with one
            !> whose answer is the open question.
            IF (b1 == b2) THEN
              diag(23) = diag(23) + REAL(SUM(CONJG(ham)*hud))
              diag(24) = diag(24) + SUM(ABS(ham)**2)
              diag(25) = diag(25) + SUM(ABS(hud)**2)
              diag(26) = MAX(diag(26), asym)
            END IF
          END IF

          !> ACCUMULATED over the spin components, not assigned: for a spinor the overlap
          !> of one component alone is not the identity, and the relation checked below
          !> holds for the whole state and not component by component.
          IF (l_ham .AND. fmpi%irank == 0 .AND. nk_local == 1) THEN
            IF (b1 == 1) pair_row(:, :, b2) = pair_row(:, :, b2) + ham
            IF (b2 == 1) pair_col(:, :, b1) = pair_col(:, :, b1) + ham
          END IF

          tmp = MATMUL(ham, vgauge(:, :, kb))
        ELSE
          tmp = MATMUL(CONJG(ovl), vgauge(:, :, kb))
        END IF
        mw = MATMUL(CONJG(TRANSPOSE(vgauge(:, :, ka))), tmp)

        !> The weights of Eq. (84), which are Wannier90's and in its units.
        DO a = 1, 3
          wa = bmesh%wb(b1)*bmesh%bk(a, b1, nk)
          DO c = 1, 3
            wc = bmesh%wb(b2)*bmesh%bk(c, b2, nk)
            o0(:, :, a, c, nk_local) = o0(:, :, a, c, nk_local) + (wa*wc)*mw
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (ovl, tmp, mw)
    IF (ALLOCATED(pnk)) DEALLOCATE (pnk, ham, ea)
    IF (ALLOCATED(hud)) DEALLOCATE (hud)
  END SUBROUTINE cf_one_k

  !> What the two routes say about each other, in three numbers.
  !>
  !> The projection is the one that carries a verdict: it is C's sign, and nothing that the
  !> direct route leaves out can turn it over, only shrink it. The ratio of magnitudes is
  !> below one by whatever the interstitial, the non-spherical potential and the local
  !> orbitals are worth at this k, so it is a scale and not an error.
  SUBROUTINE cf_report(nk, diag)
    INTEGER, INTENT(IN) :: nk
    REAL, INTENT(IN) :: diag(:)

    WRITE (oUnit, '(/,a)') 'wannierlib: C by the identity against C by a direct uHu'
    WRITE (oUnit, '(a,i8)')     '   k-point                        ', nk
    WRITE (oUnit, '(a,f12.6)')  '   Re <C_id|C_uHu> / <C_id|C_id>  ', diag(1)/diag(2)
    WRITE (oUnit, '(a,f12.6)')  '   |C_uHu| / |C_id|               ', SQRT(diag(3)/diag(2))
    WRITE (oUnit, '(a,es12.4)') '   largest R_MT surface term      ', diag(4)
    WRITE (oUnit, '(a)') '   Muffin-tin sphere and spherical potential only: the sign of the'
    WRITE (oUnit, '(a)') '   projection is the result, its magnitude is not.'
    !> The diagonal pairs on their own; see the call site for why they are worth separating.
    IF (diag(24) > 0.0) THEN
      WRITE (oUnit, '(/,a)')      '   only the DIAGONAL pairs b1 = b2, where the identity is exact'
      WRITE (oUnit, '(a,f12.6)')  '   Re <C_id|C_uHu> / <C_id|C_id>  ', diag(23)/diag(24)
      WRITE (oUnit, '(a,f12.6)')  '   |C_uHu| / |C_id|               ', SQRT(diag(25)/diag(24))
      WRITE (oUnit, '(a,es12.4)') '   largest R_MT surface term      ', diag(26)
      WRITE (oUnit, '(a)') '   Here the identity carries no error of its own, so the shortfall'
      WRITE (oUnit, '(a)') '   is what the direct route still leaves out: interstitial,'
      WRITE (oUnit, '(a)') '   non-spherical potential and local orbitals.'
    END IF
  END SUBROUTINE cf_report

END MODULE m_wannierlib_cf
