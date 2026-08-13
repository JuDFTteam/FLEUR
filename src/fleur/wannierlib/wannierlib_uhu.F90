!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The Hamiltonian between pairs of neighbours, contracted into
!>
!>      C_ab(k) = sum_{b1,b2} w_b1 b1_a w_b2 b2_b <u_{k+b1}|H_k|u_{k+b2}>
!>
!>  which is <d_a u|H|d_b u> once the finite-difference derivative of Eq. (85) is put on
!>  both sides. Eq. (84) of Lopez, Vanderbilt, Thonhauser and Souza, PRB 85, 014435 (2012);
!>  the same contraction F carries, with the Hamiltonian left in.
!>
!>  The Hamiltonian is never applied. H_k differs from the one the bra is an eigenstate
!>  of by an operator linear in the momentum,
!>
!>      H_k - H_{k+b1} = -b1.(p + k) - b1^2/2
!>
!>  so the whole matrix element follows from three things that are already available: the
!>  eigenvalue at k+b1, the pair overlap, and the momentum between the two neighbours,
!>
!>      <u_{k+b1}|H_k|u_{k+b2}> = E_{k+b1} M - b1.Pi - (b1^2/2) M .
!>
!>  That is what lets this run without ever writing a uHu file, which for a fine mesh is
!>  larger than everything else the run produces together.
!>
!>  Two sets of b vectors appear below and mixing them is silent. The identity comes from
!>  H_q = (p+q)^2/2 + V in Hartree, so the b multiplying Pi is in atomic units and goes
!>  through cell%bmat. The b carrying the weights of Eq. (84) is the one Wannier90 supplies,
!>  in inverse Angstrom, and is what puts C in the units its consumer expects.
MODULE m_wannierlib_uhu
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
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_types_melem_bmesh, ONLY: t_melem_bmesh
  USE m_types_enpara
  USE m_types_potden
  USE m_types_mpi
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_uhu
CONTAINS

  !> C on this rank's k-slice, in the Wannier gauge.
  !>
  !> Two radial tables per radial set, both indexed by the distinct b2 - b1: the plain
  !> overlap needs one and the gradient needs its own, and they are built together because
  !> they are swept over the same vectors.
  SUBROUTINE wannierlib_uhu(manifold, bmesh, kpts, atoms, cell, input, sym, noco, nococonv, &
                            radfun, jspin, l_spinors, eig_id, stars, enpara, vtot, fmpi, &
                            distk, u_matrix, u_opt, c0)
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
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: c0(:, :, :, :, :)  !> (nw,nw,3,3,nk_loc)

    COMPLEX, ALLOCATABLE :: vgauge(:, :, :), ujp(:, :, :, :, :, :), nujp(:, :, :, :, :, :, :)
    REAL, ALLOCATABLE :: kdp(:, :)
    INTEGER :: npair, nb, nw, nkl, k, kl, jcomp, jrad, lo, hi, jeig

    CALL timestart("wannierlib_uhu")

    nb = manifold%num_bands
    nw = manifold%num_wann
    nkl = COUNT(distk == fmpi%irank)
    ALLOCATE (c0(nw, nw, 3, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    IF (nb <= 0 .OR. nw <= 0 .OR. bmesh%nntot <= 0 .OR. nkl <= 0) RETURN
    IF (.NOT. ALLOCATED(bmesh%wb) .OR. .NOT. ALLOCATED(bmesh%bk)) CALL juDFT_error( &
      'wannierlib: C needs the b-shell weights, which the wannierisation produces', &
      calledby='wannierlib_uhu')

    ALLOCATE (vgauge(nb, nw, kpts%nkptf))
    DO k = 1, kpts%nkptf
      vgauge(:, :, k) = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
    END DO

    CALL uhu_diffs(bmesh, kpts%bkf, kdp, npair)

    lo = MERGE(1, jspin, l_spinors)
    hi = MERGE(2, jspin, l_spinors)
    !> The eigenvalue record: one for spinors, the spin channel otherwise. The eigenvalues
    !> belong to the state, not to the spinor component it is being split into.
    jeig = MERGE(1, jspin, l_spinors)
    DO jcomp = lo, hi
      jrad = radial_slot(radfun, jcomp)
      CALL melem_ujugaunt(atoms, cell, npair, kdp, radfun, radfun, jrad, jrad, .FALSE., 1, ujp)
      CALL melem_nablaujugaunt(atoms, cell, enpara, npair, kdp, radfun, radfun, jrad, jrad, nujp)
      kl = 0
      DO k = 1, kpts%nkptf
        IF (distk(k) /= fmpi%irank) CYCLE
        kl = kl + 1
        CALL uhu_one_k(manifold, bmesh, k, kl, kpts, ujp, nujp, kdp, npair, atoms, cell, &
                       input, sym, noco, nococonv, jcomp, jrad, jeig, eig_id, stars, &
                       enpara, vtot, fmpi, vgauge, c0)
      END DO
      DEALLOCATE (ujp, nujp)
    END DO

    CALL matrix_element_release_anchor()
    DEALLOCATE (vgauge, kdp)
    CALL timestop("wannierlib_uhu")
  END SUBROUTINE wannierlib_uhu

  !> The distinct b2 - b1 vectors, which are what both radial tables are indexed by.
  !> Identical in contract to the one the pair overlaps use: deduplicated by value, in
  !> internal coordinates, and swept over every k rather than assumed uniform.
  SUBROUTINE uhu_diffs(bmesh, bkf, kdiff_pair, npair)
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh
    REAL, INTENT(IN) :: bkf(:, :)
    REAL, ALLOCATABLE, INTENT(OUT) :: kdiff_pair(:, :)
    INTEGER, INTENT(OUT) :: npair

    REAL :: d(3)
    INTEGER :: k, b1, b2, i
    LOGICAL :: seen

    ALLOCATE (kdiff_pair(3, bmesh%nntot**2))
    kdiff_pair = 0.0
    npair = 0

    DO k = 1, SIZE(bmesh%nnlist, 1)
      DO b1 = 1, bmesh%nntot
        DO b2 = 1, bmesh%nntot
          d = bmesh%shell_vector(bkf, k, b2) - bmesh%shell_vector(bkf, k, b1)
          seen = .FALSE.
          DO i = 1, npair
            IF (ALL(ABS(kdiff_pair(:, i) - d) <= 1.0e-4)) THEN
              seen = .TRUE.
              EXIT
            END IF
          END DO
          IF (seen) CYCLE
          IF (npair == SIZE(kdiff_pair, 2)) CALL juDFT_error( &
            'wannierlib: more distinct b2-b1 vectors than pairs of neighbours', &
            calledby='uhu_diffs')
          npair = npair + 1
          kdiff_pair(:, npair) = d
        END DO
      END DO
    END DO
  END SUBROUTINE uhu_diffs

  !> One k-point's contribution, accumulated into its slice of c0.
  !>
  !> Both the overlap and the momentum arrive in the conjugated convention that
  !> melem_mmkb_sph documents, and both are undone here, per pair, before the identity is
  !> assembled: the identity mixes them with real eigenvalues and real b components, so it
  !> is only meaningful once the two are in the same convention as each other and as those.
  !>
  !> The momentum the two modules return is that of the full wavefunction, p_psi. What the
  !> identity asks for is p_u + k on the periodic part, and since psi = e^{i(k+b2).r} u the
  !> two differ by exactly -b2. That is the one term of this assembly which lives in no
  !> module, and it is subtracted here.
  SUBROUTINE uhu_one_k(manifold, bmesh, nk, nk_local, kpts, ujug_pair, nujug_pair, &
                       kdiff_pair, npair, atoms, cell, input, sym, noco, nococonv, jspin, &
                       jspin_rad, jeig, eig_id, stars, enpara, vtot, fmpi, vgauge, c0)
    TYPE(t_melem_manifold), INTENT(IN) :: manifold
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh
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
    INTEGER, INTENT(IN) :: jspin
    INTEGER, INTENT(IN) :: jspin_rad
    INTEGER, INTENT(IN) :: jeig      !> the eigenvalue record
    INTEGER, INTENT(IN) :: eig_id
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_enpara), INTENT(IN) :: enpara
    TYPE(t_potden), INTENT(IN) :: vtot
    TYPE(t_mpi), INTENT(IN) :: fmpi
    COMPLEX, INTENT(IN) :: vgauge(:, :, :)
    COMPLEX, INTENT(INOUT) :: c0(:, :, :, :, :)

    TYPE(t_mat), POINTER :: zmat_a(:), zmat_b(:)
    TYPE(t_abc), POINTER :: abc_a(:, :), abc_b(:, :)
    TYPE(t_lapw) :: lapw_a, lapw_b
    TYPE(t_spinor_layout) :: layout_a, layout_b
    COMPLEX, ALLOCATABLE :: ovl(:, :), pnk(:, :, :), ham(:, :), tmp(:, :), hw(:, :)
    REAL, ALLOCATABLE :: ea(:)
    INTEGER, ALLOCATABLE :: ev_list(:)
    INTEGER :: b1, b2, ka, kb, gpar(3), irec, nb, nw, i, a, c
    REAL :: b1c(3), b2c(3), wa, wc

    nb = manifold%num_bands
    nw = manifold%num_wann
    ev_list = [(i, i = manifold%min_band, manifold%max_band)]
    irec = MERGE(1, jspin, noco%l_noco)
    ALLOCATE (ovl(nb, nb), pnk(nb, nb, 3), ham(nb, nb), tmp(nb, nw), hw(nw, nw))

    DO b1 = 1, bmesh%nntot
      ka = bmesh%nnlist(nk, b1)
      CALL lapw_a%init(input, noco, nococonv, kpts, atoms, sym, ka, cell)
      CALL matrix_element_states(eig_id, ka, input, atoms, sym, cell, noco, nococonv, &
                                 enpara, lapw_a, vtot, fmpi, zmat_a, abc_a, ev_list=ev_list, &
                                 l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                 kpts=kpts, l_anchor=.TRUE.)
      CALL layout_a%init(input, noco, lapw_a, atoms)

      !> The eigenvalues of the bra, which is the side the Hamiltonian is resolved on. Read
      !> from the eigenvector file at the parent point: results%eig is not filled per k at
      !> the moment this runs.
      IF (.NOT. ALLOCATED(ea)) ALLOCATE (ea(input%neig))
      ea = 0.0
      CALL read_eig(eig_id, kpts%bkp(ka), jeig, eig=ea)

      !> In atomic units, because the identity is. NOT bmesh%bk, which is Wannier90's and
      !> in inverse Angstrom.
      b1c = MATMUL(bmesh%shell_vector(kpts%bkf, nk, b1), cell%bmat)

      DO b2 = 1, bmesh%nntot
        kb = bmesh%nnlist(nk, b2)
        gpar = bmesh%gkpb(:, nk, b2) - bmesh%gkpb(:, nk, b1)
        CALL lapw_b%init(input, noco, nococonv, kpts, atoms, sym, kb, cell)
        CALL matrix_element_states(eig_id, kb, input, atoms, sym, cell, noco, nococonv, &
                                   enpara, lapw_b, vtot, fmpi, zmat_b, abc_b, ev_list=ev_list, &
                                   l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                   kpts=kpts)
        CALL layout_b%init(input, noco, lapw_b, atoms)

        b2c = MATMUL(bmesh%shell_vector(kpts%bkf, nk, b2), cell%bmat)

        ovl = CMPLX(0.0, 0.0)
        CALL melem_overlap_states(stars, atoms, lapw_a, lapw_b, zmat_a(irec), zmat_b(irec), &
                                  abc_a(jspin, :), abc_b(jspin, :), jspin_rad, jspin_rad, &
                                  kpts%bkf(:, ka), kpts%bkf(:, kb), gpar, &
                                  ujug_pair, kdiff_pair, npair, &
                                  ioff_a=layout_a%row_offset(jspin), &
                                  ioff_b=layout_b%row_offset(jspin), ovl=ovl)

        pnk = CMPLX(0.0, 0.0)
        CALL melem_nabla_int(stars, cell, lapw_a, lapw_b, jspin_rad, jspin_rad, &
                             zmat_a(irec), zmat_b(irec), kpts%bkf(:, kb), gpar, pnk, &
                             ioff=layout_a%row_offset(jspin), &
                             ioff_b=layout_b%row_offset(jspin))
        CALL melem_nabla_sph(atoms, abc_a(jspin, :), abc_b(jspin, :), &
                             kpts%bkf(:, kb), gpar, kpts%bkf(:, ka), &
                             nujug_pair, kdiff_pair, npair, pnk)

        !> Out of the conjugated convention, both of them, before they are mixed.
        ovl = CONJG(ovl)
        pnk = CONJG(pnk)

        !> p_psi -> p_u + k. One term, and it depends on the pair.
        DO a = 1, 3
          pnk(:, :, a) = pnk(:, :, a) - CMPLX(b2c(a), 0.0)*ovl
        END DO

        !> The identity. The eigenvalue multiplies from the left because the bra is what it
        !> belongs to; the quadratic term is a scalar and does not care.
        ham = -CMPLX(0.5*DOT_PRODUCT(b1c, b1c), 0.0)*ovl
        DO a = 1, 3
          ham = ham - CMPLX(b1c(a), 0.0)*pnk(:, :, a)
        END DO
        DO i = 1, nb
          ham(i, :) = ham(i, :) + ea(manifold%min_band + i - 1)*ovl(i, :)
        END DO

        tmp = MATMUL(ham, vgauge(:, :, kb))
        hw  = MATMUL(CONJG(TRANSPOSE(vgauge(:, :, ka))), tmp)

        !> And the weights of Eq. (84), which are Wannier90's and in its units.
        DO a = 1, 3
          wa = bmesh%wb(b1)*bmesh%bk(a, b1, nk)
          DO c = 1, 3
            wc = bmesh%wb(b2)*bmesh%bk(c, b2, nk)
            c0(:, :, a, c, nk_local) = c0(:, :, a, c, nk_local) + (wa*wc)*hw
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (ovl, pnk, ham, tmp, hw)
    IF (ALLOCATED(ea)) DEALLOCATE (ea)
  END SUBROUTINE uhu_one_k

END MODULE m_wannierlib_uhu
