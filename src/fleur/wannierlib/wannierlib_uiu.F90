!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The pair overlaps <u_{k+b1}|u_{k+b2}> contracted into
!>
!>      F_ab(k) = sum_{b1,b2} w_b1 b1_a w_b2 b2_b <u_{k+b1}|u_{k+b2}>
!>
!>  which is the quantum geometric tensor <d_a u_n|d_b u_m> once the finite-difference
!>  derivative |d_a u_k> = sum_b w_b b_a |u_{k+b}> is put on both sides. That derivative is
!>  Eq. (85) of Lopez, Vanderbilt, Thonhauser and Souza, PRB 85, 014435 (2012), and the
!>  contraction is the structure of their Eq. (84) with the Hamiltonian left out.
!>
!>  Two things decide where this lives. It needs the wavefunctions, so it cannot run in the
!>  post-processing, which is handed slices and no way to reach the factory. And it needs
!>  the gauge at two different k-points, V(k+b1) on the left and V(k+b2) on the right, so it
!>  cannot run with the coarse matrices either, which are built before any gauge exists. It
!>  therefore sits between the two, and what it passes on is the contracted
!>  (num_wann, num_wann, 3, 3) block per k rather than the pair overlaps, which are far
!>  larger than the states they came from.
MODULE m_wannierlib_uiu
  USE m_juDFT
  USE m_melem_overlap, ONLY: melem_overlap_states
  USE m_melem_ujugaunt, ONLY: melem_ujugaunt
  USE m_matrix_element_factory, ONLY: matrix_element_states, matrix_element_release_anchor
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
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_types_melem_bmesh, ONLY: t_melem_bmesh
  USE m_types_enpara
  USE m_types_potden
  USE m_types_mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_uiu
CONTAINS

  !> F on this rank's k-slice, in the Wannier gauge.
  !>
  !> One radial table per radial set, shared by every pair: the muffin-tin half of an overlap
  !> depends on the two sides only through b2 - b1, and there are far fewer distinct
  !> differences than pairs. Both spin components accumulate into the same F when the states
  !> are spinors, which is the same rule the neighbour overlaps follow.
  SUBROUTINE wannierlib_uiu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                            radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, &
                            distk, u_matrix, u_opt, f0)
    TYPE(t_melem_manifold), INTENT(IN) :: manifold
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh
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

    COMPLEX, ALLOCATABLE :: ujp(:, :, :, :, :, :), vgauge(:, :, :)
    REAL,    ALLOCATABLE :: kdp(:, :)
    INTEGER :: npair, nb, nw, nkl, k, kl, jcomp, jrad, lo, hi

    nb = manifold%num_bands
    nw = manifold%num_wann
    nkl = COUNT(distk == fmpi%irank)
    ALLOCATE (f0(nw, nw, 3, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    IF (nb <= 0 .OR. nw <= 0 .OR. bmesh%nntot <= 0 .OR. nkl <= 0) RETURN
    IF (.NOT. ALLOCATED(bmesh%wb) .OR. .NOT. ALLOCATED(bmesh%bk)) CALL juDFT_error( &
      'wannierlib: F needs the b-shell weights, which the wannierisation produces', &
      calledby='wannierlib_uiu')

    ALLOCATE (vgauge(nb, nw, kpts%nkptf))
    DO k = 1, kpts%nkptf
      vgauge(:, :, k) = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
    END DO

    CALL bmesh%pair_diffs(kpts%bkf, kdp, npair)

    lo = MERGE(1, jspin, l_spinors)
    hi = MERGE(2, jspin, l_spinors)
    DO jcomp = lo, hi
      jrad = radial_slot(radfun, jcomp)
      CALL melem_ujugaunt(atoms, cell, npair, kdp, radfun, radfun, jrad, jrad, .FALSE., 1, ujp)
      kl = 0
      DO k = 1, kpts%nkptf
        IF (distk(k) /= fmpi%irank) CYCLE
        kl = kl + 1
        CALL uiu_one_k(manifold, bmesh, k, kl, kpts, ujp, kdp, npair, atoms, cell, input, &
                       sym, noco, nococonv, jcomp, jrad, eig_id, stars, enpara, vtot, fmpi, &
                       vgauge, f0)
      END DO
      DEALLOCATE (ujp)
    END DO

    !> The pair loop keeps one neighbour anchored while it walks the others; nothing after
    !> this holds a state by pointer, so the slot goes back to the pool.
    CALL matrix_element_release_anchor()
    DEALLOCATE (vgauge, kdp)
  END SUBROUTINE wannierlib_uiu

  !> One k-point's contribution, accumulated into its slice of f0.
  !>
  !> Every ordered pair of neighbours is visited, the diagonal b1 = b2 included: there the
  !> overlap is the identity in the ab initio gauge and the term reduces to the plain gauge
  !> product, which is what makes it worth computing rather than assuming.
  !>
  !> The overlap arrives conjugated, which is the convention melem_mmkb_sph documents and
  !> which the neighbour overlaps undo in one sweep after their own loop; there is no such
  !> sweep here, so it is undone per pair.
  !>
  !> b1 stays anchored in the factory while its partners are fetched, or the states it holds
  !> by pointer would be overwritten by the third k asked for.
  SUBROUTINE uiu_one_k(manifold, bmesh, nk, nk_local, kpts, ujug_pair, kdiff_pair, npair, &
                       atoms, cell, input, sym, noco, nococonv, jspin, jspin_rad, &
                       eig_id, stars, enpara, vtot, fmpi, vgauge, f0)
    TYPE(t_melem_manifold), INTENT(IN) :: manifold
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh
    INTEGER, INTENT(IN) :: nk, nk_local
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: ujug_pair(:, :, :, :, :, :)
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
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_mpi), INTENT(IN) :: fmpi
    COMPLEX, INTENT(IN) :: vgauge(:, :, :)
    COMPLEX, INTENT(INOUT) :: f0(:, :, :, :, :)

    TYPE(t_mat), POINTER :: zmat_a(:), zmat_b(:)
    TYPE(t_abc), POINTER :: abc_a(:, :), abc_b(:, :)
    TYPE(t_lapw) :: lapw_a, lapw_b
    TYPE(t_spinor_layout) :: layout_a, layout_b
    COMPLEX, ALLOCATABLE :: ovl(:, :), tmp(:, :), mw(:, :)
    INTEGER, ALLOCATABLE :: ev_list(:)
    INTEGER :: b1, b2, ka, kb, gpar(3), irec, nb, nw, i, a, c
    REAL :: wa, wc

    nb = manifold%num_bands
    nw = manifold%num_wann
    ev_list = [(i, i = manifold%min_band, manifold%max_band)]
    irec = MERGE(1, jspin, noco%l_noco)
    ALLOCATE (ovl(nb, nb), tmp(nb, nw), mw(nw, nw))

    DO b1 = 1, bmesh%nntot
      ka = bmesh%nnlist(nk, b1)
      CALL lapw_a%init(input, noco, nococonv, kpts, atoms, sym, ka, cell)
      CALL matrix_element_states(eig_id, ka, input, atoms, sym, cell, noco, nococonv, &
                                 enpara, lapw_a, vtot, fmpi, zmat_a, abc_a, ev_list=ev_list, &
                                 l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                 kpts=kpts, l_anchor=.TRUE.)
      CALL layout_a%init(input, noco, lapw_a, atoms)

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

        ovl = CMPLX(0.0, 0.0)
        CALL melem_overlap_states(stars, atoms, lapw_a, lapw_b, zmat_a(irec), zmat_b(irec), &
                                  abc_a(jspin, :), abc_b(jspin, :), jspin_rad, jspin_rad, &
                                  kpts%bkf(:, ka), kpts%bkf(:, kb), gpar, &
                                  ujug_pair, kdiff_pair, npair, &
                                  ioff_a=layout_a%row_offset(jspin), &
                                  ioff_b=layout_b%row_offset(jspin), ovl=ovl)

        tmp = MATMUL(CONJG(ovl), vgauge(:, :, kb))
        mw  = MATMUL(CONJG(TRANSPOSE(vgauge(:, :, ka))), tmp)

        DO a = 1, 3
          wa = bmesh%wb(b1)*bmesh%bk(a, b1, nk)
          DO c = 1, 3
            wc = bmesh%wb(b2)*bmesh%bk(c, b2, nk)
            f0(:, :, a, c, nk_local) = f0(:, :, a, c, nk_local) + (wa*wc)*mw
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (ovl, tmp, mw)
  END SUBROUTINE uiu_one_k

END MODULE m_wannierlib_uiu
