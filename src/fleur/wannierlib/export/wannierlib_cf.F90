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
  USE m_matrix_element_factory, ONLY: matrix_element_states, matrix_element_release_anchor
  USE m_eig66_io, ONLY: read_eig
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
  USE m_types_mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_uiu, wannierlib_uhu
CONTAINS

  !> F on this rank's k-slice, in the Wannier gauge.
  SUBROUTINE wannierlib_uiu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                            radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, &
                            distk, u_matrix, u_opt, f0)
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
    TYPE(t_mpi), INTENT(IN) :: fmpi
    INTEGER, INTENT(IN) :: distk(:)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)   !> (nw,nw,nk)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)      !> (nb,nw,nk)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: f0(:, :, :, :, :)  !> (nw,nw,3,3,nk_loc)

    CALL cf_contract(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, radfun, &
                     jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, distk, u_matrix, &
                     u_opt, .FALSE., f0)
  END SUBROUTINE wannierlib_uiu

  !> C on this rank's k-slice, in the Wannier gauge.
  SUBROUTINE wannierlib_uhu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                            radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, &
                            distk, u_matrix, u_opt, c0)
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
    TYPE(t_mpi), INTENT(IN) :: fmpi
    INTEGER, INTENT(IN) :: distk(:)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: c0(:, :, :, :, :)

    CALL cf_contract(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, radfun, &
                     jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, distk, u_matrix, &
                     u_opt, .TRUE., c0)
  END SUBROUTINE wannierlib_uhu

  !> The contraction both share. l_ham decides whether the Hamiltonian sits between the
  !> neighbours, which is the only difference between C and F.
  !>
  !> One radial table per radial set, shared by every pair: the muffin-tin half of an overlap
  !> depends on the two sides only through b2 - b1, and there are far fewer distinct
  !> differences than pairs. Both spin components accumulate into the same result when the
  !> states are spinors.
  SUBROUTINE cf_contract(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                         radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, &
                         distk, u_matrix, u_opt, l_ham, o0)
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
    TYPE(t_mpi), INTENT(IN) :: fmpi
    INTEGER, INTENT(IN) :: distk(:)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)
    LOGICAL, INTENT(IN) :: l_ham          !> C when true, F when false
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: o0(:, :, :, :, :)

    COMPLEX, ALLOCATABLE :: vgauge(:, :, :), ujp(:, :, :, :, :, :), nujp(:, :, :, :, :, :, :)
    REAL, ALLOCATABLE :: kdp(:, :)
    INTEGER :: npair, nb, nw, nkl, k, kl, jcomp, jrad, lo, hi, jeig

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

    ALLOCATE (vgauge(nb, nw, kpts%nkptf))
    DO k = 1, kpts%nkptf
      vgauge(:, :, k) = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
    END DO

    CALL bmesh%pair_diffs(kpts%bkf, kdp, npair)

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
        CALL cf_one_k(manifold, bmesh, k, kl, kpts, ujp, nujp, kdp, npair, atoms, cell, &
                      input, sym, noco, nococonv, jcomp, jrad, jeig, eig_id, stars, &
                      enpara, vtot, fmpi, l_ham, vgauge, o0)
      END DO
      DEALLOCATE (ujp, nujp)
    END DO

    !> The pair loop keeps one neighbour anchored while it walks the others; nothing after
    !> this holds a state by pointer, so the slot goes back to the pool.
    CALL matrix_element_release_anchor()
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
                      jspin_rad, jeig, eig_id, stars, enpara, vtot, fmpi, l_ham, vgauge, o0)
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
    TYPE(t_mpi), INTENT(IN) :: fmpi
    LOGICAL, INTENT(IN) :: l_ham
    COMPLEX, INTENT(IN) :: vgauge(:, :, :)
    COMPLEX, INTENT(INOUT) :: o0(:, :, :, :, :)

    TYPE(t_mat), POINTER :: zmat_a(:), zmat_b(:)
    TYPE(t_abc), POINTER :: abc_a(:, :), abc_b(:, :)
    TYPE(t_lapw) :: lapw_a, lapw_b
    TYPE(t_spinor_layout) :: layout_a, layout_b
    COMPLEX, ALLOCATABLE :: ovl(:, :), pnk(:, :, :), ham(:, :), tmp(:, :), mw(:, :)
    REAL, ALLOCATABLE :: ea(:)
    INTEGER, ALLOCATABLE :: ev_list(:)
    INTEGER :: b1, b2, ka, kb, gpar(3), irec, nb, nw, i, a, c
    REAL :: b1c(3), b2c(3), wa, wc

    nb = manifold%num_bands
    nw = manifold%num_wann
    ev_list = [(i, i = manifold%min_band, manifold%max_band)]
    irec = MERGE(1, jspin, noco%l_noco)
    ALLOCATE (ovl(nb, nb), tmp(nb, nw), mw(nw, nw))
    IF (l_ham) ALLOCATE (pnk(nb, nb, 3), ham(nb, nb), ea(input%neig))

    DO b1 = 1, bmesh%nntot
      ka = bmesh%nnlist(nk, b1)
      CALL lapw_a%init(input, noco, nococonv, kpts, atoms, sym, ka, cell)
      CALL matrix_element_states(eig_id, ka, input, atoms, sym, cell, noco, nococonv, &
                                 enpara, lapw_a, vtot, fmpi, zmat_a, abc_a, ev_list=ev_list, &
                                 l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                 kpts=kpts, l_anchor=.TRUE.)
      CALL layout_a%init(input, noco, lapw_a, atoms)

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

        IF (l_ham) b2c = MATMUL(bmesh%shell_vector(kpts%bkf, nk, b2), cell%bmat)

        ovl = CMPLX(0.0, 0.0)
        CALL melem_overlap_states(stars, atoms, lapw_a, lapw_b, zmat_a(irec), zmat_b(irec), &
                                  abc_a(jspin, :), abc_b(jspin, :), jspin_rad, jspin_rad, &
                                  kpts%bkf(:, ka), kpts%bkf(:, kb), gpar, &
                                  ujug_pair, kdiff_pair, npair, &
                                  ioff_a=layout_a%row_offset(jspin), &
                                  ioff_b=layout_b%row_offset(jspin), ovl=ovl)

        IF (l_ham) THEN
          !> Conjugated here and folded into the MATMUL on the other branch: the two are
          !> not interchangeable, the sum order differs and the last digits with it.
          ovl = CONJG(ovl)
          pnk = CMPLX(0.0, 0.0)
          CALL melem_nabla_int(stars, cell, lapw_a, lapw_b, jspin_rad, jspin_rad, &
                               zmat_a(irec), zmat_b(irec), kpts%bkf(:, kb), gpar, pnk, &
                               ioff=layout_a%row_offset(jspin), &
                               ioff_b=layout_b%row_offset(jspin))
          CALL melem_nabla_sph(atoms, abc_a(jspin, :), abc_b(jspin, :), &
                               kpts%bkf(:, kb), gpar, kpts%bkf(:, ka), &
                               nujug_pair, kdiff_pair, npair, pnk)
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
  END SUBROUTINE cf_one_k

END MODULE m_wannierlib_cf
