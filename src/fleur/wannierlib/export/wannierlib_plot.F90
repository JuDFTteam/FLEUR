!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Opt-in: draw the Wannier functions, by summing them here rather than by exporting
!>  the ingredients for somebody else to sum.
!>
!>      w_n(r) = (1/N) sum_k e^{ik.r} sum_m V_mn(k) u_mk(r),   V(k) = u_opt(k) . u_matrix(k)
!>
!>  The usual route writes one UNK file per k-point -- every band of that k tabulated on
!>  a real-space grid -- and lets Wannier90 read them back and do the sum. That is the
!>  summand travelling instead of the sum: gigabytes on disk for a mesh of any size, and
!>  Wannier90 then has to hold the lot. Both gauge matrices are already in memory when
!>  this is called, and so are the states, so the sum is done where the pieces are. What
!>  survives is num_wann grids instead of num_bands * num_kpts of them.
!>
!>  Nothing here is new physics. `wann_real` puts one Bloch state on one grid point --
!>  muffin-tin by spherical harmonics, interstitial by plane waves -- `wann_plot_vac`
!>  does the vacuum of a film, `m_xsf_io` writes the format XCrySDen and VESTA read, and
!>  a .vts file goes beside it for ParaView, which does not read that one.
!>  The same combination is what `wann_plot_um_dat` (FF, 2006) does in the classic path;
!>  this is that algorithm with the gauge taken from memory instead of from a .chk file,
!>  and the plotting region taken from inp.xml instead of from a side-channel plot_inp.
!>
!>  The repacking is spelled out because it is the one thing that would silently produce
!>  garbage if it were wrong. Modern FLEUR carries the matching coefficients in `t_abc`
!>  and the radial functions in `t_radfun`, both indexed by the radial-function slot `j`
!>  WITHIN an l channel:
!>
!>      abc(ntyp)%cof(band, lm, j, atom_local)      j = 1 -> u, 2 -> udot, 3.. -> LOs of that l
!>      radfun(ntyp)%r(ir, 1:2, j, l, spin)         1:2 = large/small component
!>
!>  while `wann_real` predates that unification and wants the three families apart, with
!>  the local orbitals indexed by the type's GLOBAL lo counter and the atom by its GLOBAL
!>  index. Mapping a global `ilo` back to its slot means counting how many local orbitals
!>  of the same l came before it.
MODULE m_wannierlib_plot
   USE m_juDFT
   USE m_constants
   USE m_types_abc
   USE m_types_radfun
   USE m_types_wannierlib
   USE m_types_spinor_layout, ONLY: radial_slot, t_spinor_layout
   USE m_matrix_element_factory, ONLY: matrix_element_states
   USE m_wann_real
   USE m_wann_plot_vac
   USE m_wann_2dvacabcof
   USE m_xsf_io
#ifdef CPP_MPI
   USE mpi
#endif
   USE m_types_atoms
   USE m_types_cell
   USE m_types_enpara
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_mat
   USE m_types_mpi
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_potden
   USE m_types_stars
   USE m_types_sym
   USE m_types_vacuum
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_plot_wf

CONTAINS

   !> One pass over this rank's k-points, accumulating every Wannier function at once.
   !> Called AFTER Wannier90 has run, because the gauge is what it returns; the states are
   !> read a second time, which costs an eigenvector read per k and no diagonalisation.
   SUBROUTINE wannierlib_plot_wf(this, atoms, cell, input, sym, stars, vacuum, noco, nococonv, &
                                 enpara, vtot, kpts, fmpi, eig_id, radfun, distk, &
                                 jspin, l_wannierlib_spinors, u_matrix, u_opt)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_atoms),     INTENT(IN) :: atoms
      TYPE(t_cell),      INTENT(IN) :: cell
      TYPE(t_input),     INTENT(IN) :: input
      TYPE(t_sym),       INTENT(IN) :: sym
      TYPE(t_stars),     INTENT(IN) :: stars
      TYPE(t_vacuum),    INTENT(IN) :: vacuum
      TYPE(t_noco),      INTENT(IN) :: noco
      TYPE(t_nococonv),  INTENT(IN) :: nococonv
      TYPE(t_enpara),    INTENT(IN) :: enpara
      TYPE(t_potden),    INTENT(IN) :: vtot
      TYPE(t_kpts),      INTENT(IN) :: kpts
      TYPE(t_mpi),       INTENT(IN) :: fmpi
      INTEGER,           INTENT(IN) :: eig_id
      TYPE(t_radfun),    INTENT(IN) :: radfun(:)
      INTEGER,           INTENT(IN) :: distk(:)        !> (nkptf) owning rank of each k
      INTEGER,           INTENT(IN) :: jspin
      LOGICAL,           INTENT(IN) :: l_wannierlib_spinors
      COMPLEX,           INTENT(IN) :: u_matrix(:, :, :)   !> (nw,nw,nk)
      COMPLEX,           INTENT(IN) :: u_opt(:, :, :)      !> (nb,nw,nk)

      REAL    :: vec1(3), vec2(3), vec3(3), zero(3), point(3), pt(3), poinint(3), phas
      INTEGER :: grid(3), mesh, posi, ix, iy, iz, i1, i2, i3, ii3
      LOGICAL :: twodim
      INTEGER :: ikpt, nbn, nb, nw, ib, ierr, irec, addnoco, jspin_rad, nv2d
      INTEGER :: jspin_comp, icomp, ncomp, kidx
      INTEGER :: nt, nq, na, ivac, jvac
      REAL    :: s
      COMPLEX :: xdnout, factor
      INTEGER, ALLOCATABLE :: ev_list(:)
      COMPLEX, ALLOCATABLE :: wf(:, :, :), vgauge(:, :)
      COMPLEX, ALLOCATABLE :: acof(:, :, :), bcof(:, :, :), ccof(:, :, :, :)
      REAL,    ALLOCATABLE :: ff(:, :, :, :), gg(:, :, :, :), flo(:, :, :, :)
      REAL,    ALLOCATABLE :: vz(:, :, :), u(:, :, :), ue(:, :, :)
      COMPLEX, ALLOCATABLE :: ac(:, :, :), bc(:, :, :)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat), POINTER :: zmat_p(:)
      TYPE(t_abc), POINTER :: abc_p(:, :)
      TYPE(t_spinor_layout) :: slay

      IF (.NOT. this%export%plot) RETURN

      CALL timestart('wannierlib_plot_wf')

      nb = this%num_bands
      nw = this%num_wann
      !> A spinor is drawn component by component and the two are combined only at the
      !> end, because neither half is a function on its own: what can be plotted is the
      !> norm of the pair. A scalar state has one component and the loop below runs once.
      ncomp = MERGE(2, 1, l_wannierlib_spinors)

      CALL plot_region(cell, this%export%plot_grid, this%export%plot_cells, &
                       twodim, grid, vec1, vec2, vec3, zero)
      mesh = grid(1)*grid(2)*grid(3)

      ALLOCATE (wf(nw, mesh, ncomp), source=CMPLX(0.0, 0.0), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error('wannierlib_plot_wf failed allocating the grid', &
                                      calledby='wannierlib_plot_wf')
      ALLOCATE (vgauge(nb, nw))
      !> The vacuum potential each spin component sees. One per component and not one
      !> per channel: a spinor whose two halves were expanded in the same vacuum potential
      !> is a film with no exchange splitting outside the slab, which is not this one.
      !> The fourth extent is clamped, because a single potential serves both components.
      ALLOCATE (vz(vacuum%nmzd, 2, ncomp), source=0.0)
      IF (input%film .AND. ALLOCATED(vtot%vac)) THEN
         DO icomp = 1, ncomp
            jspin_comp = MERGE(icomp, jspin, l_wannierlib_spinors)
            vz(:, 1:MIN(2, SIZE(vtot%vac, 3)), icomp) = &
               REAL(vtot%vac(:, 1, 1:MIN(2, SIZE(vtot%vac, 3)), MIN(jspin_comp, SIZE(vtot%vac, 4))))
         END DO
      END IF

      ev_list = [(ib, ib=this%min_band, this%max_band)]

      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE

         !> The full Wannier gauge, as the adapter documents it: disentangle, then rotate.
         !> Without disentanglement u_opt is the identity and this is u_matrix alone.
         vgauge = MATMUL(u_opt(:, :, ikpt), u_matrix(:, :, ikpt))

         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)
         nv2d = lapw%dim_nv2d()   ! only defined once this k's basis is set up
         !> l_both_spinors is what makes the factory hand back the two halves of one
         !> second-variation spinor; without it the records are independent channels.
         CALL matrix_element_states(eig_id, ikpt, input, atoms, sym, cell, noco, &
                                    nococonv, enpara, lapw, vtot, fmpi, zmat_p, abc_p, &
                                    ev_list=ev_list, &
                                    l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                    kpts=kpts)
         !> Where a spin block begins inside one record: zero when the two components are
         !> two records, row_dn when one record already holds the whole spinor. Asked of
         !> t_spinor_layout rather than re-derived from the flags here.
         CALL slay%init(input, noco, lapw, atoms)

         !> jspin_comp is the spin DIRECTION this pass draws. For a scalar state it is the
         !> channel being wannierised and the loop runs once; for a spinor it is the
         !> component, and the two are summed onto separate grids because neither half is
         !> a function on its own.
         DO jspin_comp = MERGE(1, jspin, l_wannierlib_spinors), &
                         MERGE(2, jspin, l_wannierlib_spinors)
            icomp     = MERGE(jspin_comp, 1, l_wannierlib_spinors)
            irec      = MERGE(1, jspin_comp, noco%l_noco)
            addnoco   = slay%row_offset(jspin_comp)
            jspin_rad = radial_slot(radfun, jspin_comp)
            !> lapw keeps one set of G vectors per channel and only a jspins=2 run has a
            !> second one, so a spinor over a single potential reads the first set twice.
            kidx      = MERGE(jspin_comp, 1, input%jspins == 2)

            CALL repack_coefficients(atoms, sym, noco, radfun, abc_p(jspin_comp, :), jspin_rad, nb, &
                                     ff, gg, flo, acof, bcof, ccof)

            IF (input%film) THEN
               IF (ALLOCATED(ac)) DEALLOCATE (ac, bc, u, ue)
               ALLOCATE (ac(nv2d, nb, 2), bc(nv2d, nb, 2), u(vacuum%nmzd, nv2d, vacuum%nvac), &
                         ue(vacuum%nmzd, nv2d, vacuum%nvac))
               CALL wann_2dvacabcof(nv2d, nb, vacuum%nvac, vacuum%nmzd, vacuum%nmz, cell%omtil, &
                                    vz(:, :, icomp), &
                                    lapw%nv(kidx), lapw%bkpt, cell%z1, lapw%dim_nvd(), &
                                    lapw%k1(:, kidx), lapw%k2(:, kidx), lapw%k3(:, kidx), &
                                    enpara%evac0(:, kidx), cell%bbmat, vacuum%delz, cell%bmat, &
                                    lapw%dim_nbasfcn(), input%neig, zmat_p(irec), &
                                    ac, bc, u, ue, addnoco, noco%l_ss, nococonv%qss, jspin_comp)
            END IF

            ii3 = MERGE(0, 3, input%film)

            DO nbn = 1, nb
               DO iz = 0, grid(3) - 1
                  DO iy = 0, grid(2) - 1
                     xloop: DO ix = 0, grid(1) - 1
                        posi = ix + 1 + iy*grid(1) + iz*grid(1)*grid(2)
                        point = zero + vec1*(REAL(ix)/REAL(grid(1) - 1)) + vec2*(REAL(iy)/REAL(grid(2) - 1))
                        IF (.NOT. twodim) point = point + vec3*(REAL(iz)/REAL(grid(3) - 1))

                        !> e^{ik.r}: wann_real returns the periodic part, the Bloch phase is here.
                        poinint = MATMUL(cell%bmat, point)/tpi_const
                        phas = tpi_const*DOT_PRODUCT(lapw%bkpt, poinint)
                        factor = CMPLX(COS(phas), SIN(phas))

                        !> inside a muffin tin? the sphere may be the one of a neighbouring cell
                        DO i1 = -3, 3
                           DO i2 = -3, 3
                              DO i3 = -ii3, ii3
                                 pt = point + MATMUL(cell%amat, [REAL(i1), REAL(i2), REAL(i3)])
                                 na = 0
                                 DO nt = 1, atoms%ntype
                                    DO nq = 1, atoms%neq(nt)
                                       na = na + 1
                                       s = NORM2(atoms%pos(:, na) - pt)
                                       IF (s < atoms%rmsh(atoms%jri(nt), nt)) THEN
                                          CALL eval_state(pt, nt, na, 0, 1)
                                          wf(:, posi, icomp) = wf(:, posi, icomp) + xdnout*vgauge(nbn, :)*factor
                                          CYCLE xloop
                                       END IF
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO

                        !> vacuum of a film
                        IF (input%film .AND. ABS(point(3)) >= cell%z1) THEN
                           ivac = MERGE(2, 1, point(3) < 0.0)
                           jvac = MERGE(1, ivac, vacuum%nvac == 1)
                           CALL wann_plot_vac(point, cell%z1, vacuum%nmzd, nv2d, stars%ng3, &
                                              vacuum%nvac, vacuum%nmz, vacuum%delz, cell%bmat, &
                                              cell%bbmat, enpara%evac0(:, kidx), lapw%bkpt, vz(:, :, icomp), jspin_comp, &
                                              lapw%k1(:, kidx), lapw%k2(:, kidx), lapw%k3(:, kidx), &
                                              lapw%dim_nvd(), lapw%dim_nbasfcn(), input%neig, &
                                              lapw%nv(kidx), cell%omtil, nb, &
                                              ac(:, nbn, ivac), bc(:, nbn, ivac), &
                                              u(:, :, jvac), ue(:, :, jvac), xdnout)
                           wf(:, posi, icomp) = wf(:, posi, icomp) + xdnout*vgauge(nbn, :)*factor
                           CYCLE xloop
                        END IF

                        !> interstitial
                        CALL eval_state(point, 0, 0, 0, 2)
                        wf(:, posi, icomp) = wf(:, posi, icomp) + xdnout*vgauge(nbn, :)*factor
                     END DO xloop
                  END DO
               END DO
            END DO
         END DO
      END DO

      IF (ALLOCATED(ac)) DEALLOCATE (ac, bc, u, ue)

#ifdef CPP_MPI
      CALL collect_wf(wf, fmpi)
#endif
      wf = wf/REAL(kpts%nkptf)

      IF (fmpi%irank == 0) CALL write_plots(wf, atoms, cell, input, twodim, grid, &
                                          vec1, vec2, vec3, zero, jspin)

      CALL timestop('wannierlib_plot_wf')

   CONTAINS

      !> The `wann_real` argument list, in one place so the loop above stays readable.
      SUBROUTINE eval_state(p, n, nat, iv, iflag)
         REAL,    INTENT(IN) :: p(3)
         INTEGER, INTENT(IN) :: n, nat, iv, iflag
         REAL :: pp(3)
         pp = p   ! wann_real takes it INTENT(INOUT)
         CALL wann_real(pp, n, nat, iv, iflag, lapw%bkpt, nbn, &
                        stars%ng3, vacuum%nmzxyd, stars%ng2, MAXVAL(sym%ntypsy), &
                        atoms%lmaxd, atoms%jmtd, atoms%nat, atoms%ntype, vacuum%nmzd, &
                        sym%nop, sym%nop2, sym%mrot, sym%tau, sym%invtab, &
                        stars%ng3, vacuum%nvac, sym%invs, cell%z1, vacuum%delz, &
                        vacuum%nmz, vacuum%nmzxy, stars%ng2, &
                        atoms%lmax, atoms%rmsh, atoms%jri, atoms%pos, sym%ngopr, sym%ntypsy, &
                        lapw%dim_nvd(), cell%omtil, cell%amat, cell%bmat, &
                        atoms%nlod, atoms%llod, atoms%nlo, atoms%llo, &
                        ff, gg, flo, acof(nbn, :, :), bcof(nbn, :, :), ccof(:, nbn, :, :), &
                        zmat_p(irec), lapw%nv(kidx), &
                        lapw%k1(:, kidx), lapw%k2(:, kidx), lapw%k3(:, kidx), &
                        atoms%lmaxd*(atoms%lmaxd + 2), lapw%dim_nbasfcn(), &
                        noco%l_ss, nococonv%qss, jspin_comp, addnoco, xdnout)
      END SUBROUTINE eval_state

   END SUBROUTINE wannierlib_plot_wf

   !> Where to draw: one unit cell, centred on the origin, on a fixed grid.
   !>
   !> Deliberately not taken from <plotting>/<plot>. That block exists and carries exactly
   !> the right fields, but it is written by default -- inpgen and AiiDA both emit a 2D
   !> density-plot region nobody asked for -- so borrowing it would silently draw a flat
   !> slice of somebody else's plot instead of the Wannier function. A region that is
   !> always the same is one less thing to get wrong; whether it is big enough is measured
   !> and reported rather than left to be assumed.
   SUBROUTINE plot_region(cell, want_grid, want_cells, twodim, grid, vec1, vec2, vec3, zero)
      TYPE(t_cell), INTENT(IN)  :: cell
      INTEGER, INTENT(IN)  :: want_grid(3)
      INTEGER, INTENT(IN)  :: want_cells(3)
      LOGICAL, INTENT(OUT) :: twodim
      INTEGER, INTENT(OUT) :: grid(3)
      REAL,    INTENT(OUT) :: vec1(3), vec2(3), vec3(3), zero(3)

      twodim = .FALSE.
      !> The region is what the input asks for, not one cell. A Wannier function wider than
      !> the cell is cut by the walls, and the loop below already evaluates anywhere: the
      !> Bloch phase comes from the real coordinate and the muffin-tin search spans +-3
      !> cells. The grid scales with the region, so the resolution is the one asked for.
      grid   = want_grid*want_cells
      vec1   = cell%amat(:, 1)*want_cells(1)
      vec2   = cell%amat(:, 2)*want_cells(2)
      vec3   = cell%amat(:, 3)*want_cells(3)
      zero   = -0.5*(vec1 + vec2 + vec3)
   END SUBROUTINE plot_region

   !> t_abc / t_radfun -> the separate acof/bcof/ccof and ff/gg/flo that wann_real wants.
   SUBROUTINE repack_coefficients(atoms, sym, noco, radfun, abc, jspin_rad, nbands, &
                                  ff, gg, flo, acof, bcof, ccof)
      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_sym),    INTENT(IN) :: sym
      TYPE(t_noco),   INTENT(IN) :: noco
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      TYPE(t_abc),    INTENT(IN) :: abc(:)
      INTEGER,        INTENT(IN) :: jspin_rad, nbands
      REAL,    ALLOCATABLE, INTENT(INOUT) :: ff(:, :, :, :), gg(:, :, :, :), flo(:, :, :, :)
      COMPLEX, ALLOCATABLE, INTENT(INOUT) :: acof(:, :, :), bcof(:, :, :), ccof(:, :, :, :)

      INTEGER :: lmd, nlod, llod, ntyp, na, nat_local, l, m, lm, ilo, ir

      lmd  = atoms%lmaxd*(atoms%lmaxd + 2)
      nlod = MAX(1, atoms%nlod)
      llod = MAX(1, atoms%llod)

      IF (.NOT. ALLOCATED(acof)) THEN
         ALLOCATE (acof(nbands, 0:lmd, atoms%nat), bcof(nbands, 0:lmd, atoms%nat))
         ALLOCATE (ccof(-llod:llod, nbands, nlod, atoms%nat))
         ALLOCATE (ff(atoms%ntype, atoms%jmtd, 2, 0:atoms%lmaxd))
         ALLOCATE (gg(atoms%ntype, atoms%jmtd, 2, 0:atoms%lmaxd))
         ALLOCATE (flo(atoms%ntype, atoms%jmtd, 2, nlod))
      END IF
      acof = CMPLX(0.0, 0.0); bcof = CMPLX(0.0, 0.0); ccof = CMPLX(0.0, 0.0)
      ff = 0.0; gg = 0.0; flo = 0.0

      !> --- radial functions: slot 1 is u, slot 2 is udot, the rest are the LOs ---
      DO ntyp = 1, atoms%ntype
         DO l = 0, atoms%lmax(ntyp)
            DO ir = 1, atoms%jri(ntyp)
               ff(ntyp, ir, 1:2, l) = radfun(ntyp)%r(ir, 1:2, 1, l, jspin_rad)
               IF (abc(ntyp)%n_r(l) >= 2) gg(ntyp, ir, 1:2, l) = radfun(ntyp)%r(ir, 1:2, 2, l, jspin_rad)
            END DO
         END DO
         DO ilo = 1, atoms%nlo(ntyp)
            l = atoms%llo(ilo, ntyp)
            DO ir = 1, atoms%jri(ntyp)
               flo(ntyp, ir, 1:2, ilo) = radfun(ntyp)%r(ir, 1:2, atoms%slot_of_lo(ilo, ntyp), l, jspin_rad)
            END DO
         END DO
      END DO

      !> --- matching coefficients: same slot convention, atom index made global ---
      na = 0
      DO ntyp = 1, atoms%ntype
         DO nat_local = 1, atoms%neq(ntyp)
            na = na + 1
            DO l = 0, atoms%lmax(ntyp)
               DO m = -l, l
                  lm = l*(l + 1) + m
                  acof(:, lm, na) = abc(ntyp)%cof(1:nbands, lm, 1, nat_local)
                  IF (abc(ntyp)%n_r(l) >= 2) bcof(:, lm, na) = abc(ntyp)%cof(1:nbands, lm, 2, nat_local)
               END DO
            END DO
            DO ilo = 1, atoms%nlo(ntyp)
               l = atoms%llo(ilo, ntyp)
                  DO m = -l, l
                  lm = l*(l + 1) + m
                  ccof(m, :, ilo, na) = abc(ntyp)%cof(1:nbands, lm, atoms%slot_of_lo(ilo, ntyp), nat_local)
               END DO
            END DO
         END DO
      END DO

      !> An atom that abcof reached through its inversion partner carries coefficients in
      !> the rotated frame; the plot needs the global one. The guard is the producer's own
      !> condition (t_abc: any(invsat==2) .and. .not.l_soc) -- with SOC the partners were
      !> built directly and flipping them here would be the error, not the fix.
      !> Nothing turns the coefficients here: they arrive expanded about the global frame
      !> already, which is the frame the real-space grid below is laid out in. An expansion
      !> about each atom's own frame would have to be undone first, and for an atom reached
      !> through its inversion partner that means undoing the partner's rotation and the
      !> inversion with it.
   END SUBROUTINE repack_coefficients

#ifdef CPP_MPI
   !> Each rank summed its own share of the k-mesh; the whole sum lives on rank 0.
   SUBROUTINE collect_wf(wf, fmpi)
      COMPLEX,     INTENT(INOUT) :: wf(:, :, :)
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      INTEGER :: ierr
      IF (fmpi%isize == 1) RETURN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, wf, SIZE(wf), MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
   END SUBROUTINE collect_wf
#endif

   !> One XSF per Wannier function, two data blocks in it.
   !>
   !> For a scalar state they are the real part -- which is what a two-sign isosurface plot
   !> is made of -- and the modulus. The imaginary part is not written but its weight is
   !> reported, because that number is the check on the phase fix below and a file nobody
   !> opens is not.
   !>
   !> A spinor has neither. Each component carries its own phase, so there is no global one
   !> to divide out and no real part to contour; dividing both by the phase of one of them
   !> would make a picture that depends on which was picked. What is drawn instead is the
   !> norm of the pair, sqrt(|w_up|^2 + |w_dn|^2), and beside it the spin density
   !> |w_up|^2 - |w_dn|^2 -- both unchanged by either phase. What is reported is then the
   !> spin polarisation of the function, which says whether spin-orbit coupling mixed it or
   !> only tilted it, and a value far from +-1 in a weakly coupled magnet is the signal
   !> that the two components were not the two halves of the same state.
   SUBROUTINE write_plots(wf, atoms, cell, input, twodim, grid, vec1, vec2, vec3, zero, jspin)
      COMPLEX,       INTENT(INOUT) :: wf(:, :, :)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_input), INTENT(IN) :: input
      LOGICAL,       INTENT(IN) :: twodim
      INTEGER,       INTENT(IN) :: grid(3), jspin
      REAL,          INTENT(IN) :: vec1(3), vec2(3), vec3(3), zero(3)

      INTEGER :: nw, mesh, ncomp, nplo, ipos, imax, iunit, ix, iy, iz
      REAL    :: total, edge, wup, wdn
      REAL    :: reweight, imweight, report
      COMPLEX :: phase
      LOGICAL :: l_spinor
      REAL, ALLOCATABLE :: fld1(:), fld2(:), mag(:)
      CHARACTER(LEN=8)  :: nm1, nm2
      CHARACTER(LEN=40) :: fname

      nw    = SIZE(wf, 1)
      mesh  = SIZE(wf, 2)
      ncomp = SIZE(wf, 3)
      l_spinor = ncomp == 2
      ALLOCATE (fld1(mesh), fld2(mesh), mag(mesh))
      nm1 = MERGE('norm    ', 'real    ', l_spinor)
      nm2 = MERGE('spin    ', 'modulus ', l_spinor)

      WRITE (oUnit, '(/,a)') 'Wannier functions on a real-space grid (plotWF)'
      WRITE (oUnit, '(a,3i5)') '  grid           ', grid
      WRITE (oUnit, '(a,i5)')  '  functions      ', nw
      IF (l_spinor) THEN
         WRITE (oUnit, '(a)')  '  spinor: each file holds the norm of the pair and the spin density'
         WRITE (oUnit, '(a)')  '   WF    spin polarisation   weight on the box edge'
         WRITE (oUnit, '(a)')  '                             (+1 is pure up, -1 pure down; a large'
         WRITE (oUnit, '(a)')  '                              edge weight means the cell clips it)'
      ELSE
         WRITE (oUnit, '(a)')  '   WF    max|Im|/max|Re|   weight on the box edge'
         WRITE (oUnit, '(a)')  '                             (both should be small; a large edge'
         WRITE (oUnit, '(a)')  '                              weight means the cell clips the function)'
      END IF

      DO nplo = 1, nw
         IF (l_spinor) THEN
            !> No phase is taken off: the two fields below do not have one.
            DO ipos = 1, mesh
               wup = ABS(wf(nplo, ipos, 1))**2
               wdn = ABS(wf(nplo, ipos, 2))**2
               fld1(ipos) = SQRT(wup + wdn)
               fld2(ipos) = wup - wdn
            END DO
            mag = fld1
            wup = SUM(ABS(wf(nplo, :, 1))**2)
            wdn = SUM(ABS(wf(nplo, :, 2))**2)
            report = 0.0
            IF (wup + wdn > 0.0) report = (wup - wdn)/(wup + wdn)
         ELSE
            !> A Wannier function is real up to one global phase; take it off, using the point
            !> where the function is largest so the phase is well determined there.
            imax = 1
            DO ipos = 1, mesh
               IF (ABS(wf(nplo, ipos, 1)) > ABS(wf(nplo, imax, 1))) imax = ipos
            END DO
            phase = wf(nplo, imax, 1)/ABS(wf(nplo, imax, 1))
            wf(nplo, :, 1) = wf(nplo, :, 1)/phase

            reweight = MAXVAL(ABS(REAL(wf(nplo, :, 1))))
            imweight = MAXVAL(ABS(AIMAG(wf(nplo, :, 1))))
            IF (reweight > 0.0) imweight = imweight/reweight
            fld1 = REAL(wf(nplo, :, 1))
            fld2 = ABS(wf(nplo, :, 1))
            mag  = fld2
            report = imweight
         END IF

         !> How much of the function sits on the outermost layer of the box. If this is not
         !> small the unit cell is too tight for this Wannier function and the picture is cut.
         total = 0.0; edge = 0.0
         DO iz = 0, grid(3) - 1
            DO iy = 0, grid(2) - 1
               DO ix = 0, grid(1) - 1
                  ipos = ix + 1 + iy*grid(1) + iz*grid(1)*grid(2)
                  total = total + mag(ipos)**2
                  IF (ix == 0 .OR. ix == grid(1) - 1 .OR. iy == 0 .OR. iy == grid(2) - 1 .OR. &
                      (.NOT. twodim .AND. (iz == 0 .OR. iz == grid(3) - 1))) &
                     edge = edge + mag(ipos)**2
               END DO
            END DO
         END DO
         IF (total > 0.0) edge = edge/total
         WRITE (oUnit, '(i5,3x,es12.4,6x,es12.4)') nplo, report, edge

         WRITE (fname, '(a,i1,a,i3.3,a)') 'WF', jspin, '_', nplo, '.xsf'
         OPEN (NEWUNIT=iunit, file=TRIM(fname))
         CALL xsf_WRITE_atoms(iunit, atoms, input%film, cell%amat)
         CALL xsf_WRITE_header(iunit, twodim, TRIM(fname), vec1, vec2, vec3, zero, grid)
         DO ipos = 1, mesh
            WRITE (iunit, *) fld1(ipos)
         END DO
         CALL xsf_WRITE_newblock(iunit, twodim, vec1, vec2, vec3, zero, grid)
         DO ipos = 1, mesh
            WRITE (iunit, *) fld2(ipos)
         END DO
         CALL xsf_WRITE_endblock(iunit, twodim)
         CLOSE (iunit)

         !> The same function again for ParaView, which does not read XSF. Same grid and
         !> the same two fields, built once above, so the two files cannot disagree about
         !> what they are showing.
         WRITE (fname, '(a,i1,a,i3.3,a)') 'WF', jspin, '_', nplo, '.vts'
         CALL write_vts(TRIM(fname), fld1, nm1, fld2, nm2, mag, l_spinor, &
                        twodim, grid, vec1, vec2, vec3, zero)
      END DO

      WRITE (oUnit, '(a,i0,a,i0,a)') '  wrote ', nw, ' .xsf and ', nw, ' .vts files '// &
         '(real part and modulus in each)'
   END SUBROUTINE write_plots


   !> One Wannier function for ParaView, as a VTK structured grid (.vts).
   !>
   !> STRUCTURED GRID AND NOT IMAGE DATA (.vti), although the grid is regular: an image
   !> dataset stores one spacing per CARTESIAN axis, and the three vectors spanned here are
   !> the Bravais vectors, which meet at 60 degrees in an fcc or bcc cell. Written that way
   !> a Wannier function comes out sheared into a box it does not live in -- a picture that
   !> looks perfectly fine and is wrong. So every point carries its own cartesian coordinate.
   !> This was measured the hard way in August 2026 on FeO2 before it was understood.
   !>
   !> Angstrom, like the XSF beside it, so the two agree with each other and with the atom
   !> positions anyone loads next to them.
   SUBROUTINE write_vts(fname, fld1, nm1, fld2, nm2, mag, l_spinor, &
                        twodim, grid, vec1, vec2, vec3, zero)
      CHARACTER(LEN=*), INTENT(IN) :: fname, nm1, nm2
      !> The two fields to write and their names, plus the magnitude the isosurface levels
      !> are taken from -- the modulus of a scalar state, the norm of a spinor.
      REAL,    INTENT(IN) :: fld1(:), fld2(:), mag(:)
      LOGICAL, INTENT(IN) :: l_spinor
      LOGICAL, INTENT(IN) :: twodim
      INTEGER, INTENT(IN) :: grid(3)
      REAL,    INTENT(IN) :: vec1(3), vec2(3), vec3(3), zero(3)

      INTEGER :: iunit, ix, iy, iz, nz, mesh, ipos
      REAL    :: p(3)
      CHARACTER(LEN=40) :: ext

      !> A 2D region spans vec1 and vec2 only, and then the third dimension is a single
      !> layer rather than grid(3) of them, which is also what wf was filled with.
      nz   = MERGE(1, grid(3), twodim)
      mesh = grid(1)*grid(2)*nz
      WRITE (ext, '(6(i0,1x))') 0, grid(1) - 1, 0, grid(2) - 1, 0, nz - 1

      OPEN (NEWUNIT=iunit, file=fname)
      WRITE (iunit, '(a)') '<?xml version="1.0"?>'
      !> The level to contour at, carried in the file rather than left in `out`: these are
      !> copied to a laptop to make figures and the run directory stays behind. Chosen by
      !> ENCLOSED CHARGE and not as a fraction of the maximum -- 0.45*max encloses under
      !> half of the function and draws a stunted lobe, which is how the first pictures of
      !> these came out.
      IF (l_spinor) THEN
         WRITE (iunit, '(a,3(a,f10.6),a)') '<!-- isosurface by enclosed charge: ', &
            '  70% ', iso_by_charge(mag, 0.70), &
            '  85% ', iso_by_charge(mag, 0.85), &
            '  95% ', iso_by_charge(mag, 0.95), '   (contour the norm at c) -->'
      ELSE
         WRITE (iunit, '(a,3(a,f10.6),a)') '<!-- isosurface by enclosed charge: ', &
            '  70% ', iso_by_charge(mag, 0.70), &
            '  85% ', iso_by_charge(mag, 0.85), &
            '  95% ', iso_by_charge(mag, 0.95), '   (contour the real part at +c and -c) -->'
      END IF
      WRITE (iunit, '(a)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
      WRITE (iunit, '(a)') '  <StructuredGrid WholeExtent="'//TRIM(ext)//'">'
      WRITE (iunit, '(a)') '    <Piece Extent="'//TRIM(ext)//'">'
      WRITE (iunit, '(a)') '      <Points>'
      WRITE (iunit, '(a)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      !> x fastest, then y, then z: the order wf was filled in, and the order VTK reads.
      !> The last point along each direction repeats the first cell image, so the vectors
      !> span the whole cell and the step is v_i/(n_i - 1).
      DO iz = 0, nz - 1
         DO iy = 0, grid(2) - 1
            DO ix = 0, grid(1) - 1
               p = zero + vec1*(REAL(ix)/REAL(grid(1) - 1)) + vec2*(REAL(iy)/REAL(grid(2) - 1))
               IF (.NOT. twodim) p = p + vec3*(REAL(iz)/REAL(grid(3) - 1))
               WRITE (iunit, '(3(f14.7,1x))') p*bohr_to_angstrom_const
            END DO
         END DO
      END DO
      WRITE (iunit, '(a)') '        </DataArray>'
      WRITE (iunit, '(a)') '      </Points>'
      WRITE (iunit, '(a)') '      <PointData Scalars="'//TRIM(nm1)//'">'
      WRITE (iunit, '(a)') '        <DataArray type="Float32" Name="'//TRIM(nm1)//'" format="ascii">'
      DO ipos = 1, mesh
         WRITE (iunit, '(es16.7)') fld1(ipos)
      END DO
      WRITE (iunit, '(a)') '        </DataArray>'
      WRITE (iunit, '(a)') '        <DataArray type="Float32" Name="'//TRIM(nm2)//'" format="ascii">'
      DO ipos = 1, mesh
         WRITE (iunit, '(es16.7)') fld2(ipos)
      END DO
      WRITE (iunit, '(a)') '        </DataArray>'
      WRITE (iunit, '(a)') '      </PointData>'
      WRITE (iunit, '(a)') '    </Piece>'
      WRITE (iunit, '(a)') '  </StructuredGrid>'
      WRITE (iunit, '(a)') '</VTKFile>'
      CLOSE (iunit)
   END SUBROUTINE write_vts


   !> The |w| at which the region |w| > c holds `frac` of the integral of |w|^2.
   !>
   !> By bisection rather than by sorting: the enclosed weight falls monotonically with c,
   !> so forty halvings pin it down to one part in 10^12 of the peak, and nothing has to be
   !> ordered. The volume element cancels, being the same at every point of the grid.
   FUNCTION iso_by_charge(mag, frac) RESULT(c)
      REAL, INTENT(IN) :: mag(:)
      REAL, INTENT(IN) :: frac
      REAL :: c, lo, hi, total, inside
      INTEGER :: it

      total = SUM(mag**2)
      c = 0.0
      IF (total <= 0.0) RETURN
      lo = 0.0; hi = MAXVAL(mag)
      DO it = 1, 40
         c = 0.5*(lo + hi)
         inside = SUM(mag**2, MASK=mag > c)
         IF (inside > frac*total) THEN
            lo = c      ! still enclosing too much: the level has to come up
         ELSE
            hi = c
         END IF
      END DO
   END FUNCTION iso_by_charge



END MODULE m_wannierlib_plot
