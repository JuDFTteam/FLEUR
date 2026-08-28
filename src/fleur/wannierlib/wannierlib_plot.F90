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
!>  does the vacuum of a film, and `m_xsf_io` writes the format XCrySDen and VESTA read.
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
   USE m_types
   USE m_types_abc
   USE m_types_radfun
   USE m_types_wannierlib
   USE m_types_spinor_layout, ONLY: radial_slot
   USE m_matrix_element_factory, ONLY: matrix_element_states
   USE m_wann_real
   USE m_wann_abinv
   USE m_wann_plot_vac
   USE m_wann_2dvacabcof
   USE m_xsf_io
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
      INTEGER :: nt, nq, na, ivac, jvac
      REAL    :: s
      COMPLEX :: xdnout, factor
      INTEGER, ALLOCATABLE :: ev_list(:)
      COMPLEX, ALLOCATABLE :: wf(:, :), vgauge(:, :)
      COMPLEX, ALLOCATABLE :: acof(:, :, :), bcof(:, :, :), ccof(:, :, :, :)
      REAL,    ALLOCATABLE :: ff(:, :, :, :), gg(:, :, :, :), flo(:, :, :, :)
      REAL,    ALLOCATABLE :: vz(:, :), u(:, :, :), ue(:, :, :)
      COMPLEX, ALLOCATABLE :: ac(:, :, :), bc(:, :, :)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat), POINTER :: zmat_p(:)
      TYPE(t_abc), POINTER :: abc_p(:, :)

      IF (.NOT. this%l_plot_wf) RETURN
      !> The spinor case would need both components summed and only their norm drawn, and
      !> the phase convention below stops meaning anything. Refused rather than guessed.
      IF (l_wannierlib_spinors) CALL juDFT_error( &
         'wannierlib plotWF is not implemented for spinor Wannier functions (noco+SOC)', &
         calledby='wannierlib_plot_wf')

      CALL timestart('wannierlib_plot_wf')

      nb = this%num_bands
      nw = this%num_wann
      jspin_rad = radial_slot(radfun, jspin)

      CALL plot_region(cell, twodim, grid, vec1, vec2, vec3, zero)
      mesh = grid(1)*grid(2)*grid(3)

      ALLOCATE (wf(nw, mesh), source=CMPLX(0.0, 0.0), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error('wannierlib_plot_wf failed allocating the grid', &
                                      calledby='wannierlib_plot_wf')
      ALLOCATE (vgauge(nb, nw))
      ALLOCATE (vz(vacuum%nmzd, 2), source=0.0)
      IF (input%film .AND. ALLOCATED(vtot%vac)) &
         vz(:, 1:MIN(2, SIZE(vtot%vac, 3))) = REAL(vtot%vac(:, 1, 1:MIN(2, SIZE(vtot%vac, 3)), jspin))

      ev_list = [(ib, ib=this%min_band, this%max_band)]

      DO ikpt = 1, kpts%nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE

         !> The full Wannier gauge, as the adapter documents it: disentangle, then rotate.
         !> Without disentanglement u_opt is the identity and this is u_matrix alone.
         vgauge = MATMUL(u_opt(:, :, ikpt), u_matrix(:, :, ikpt))

         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)
         nv2d = lapw%dim_nv2d()   ! only defined once this k's basis is set up
         CALL matrix_element_states(eig_id, ikpt, input, atoms, sym, cell, noco, &
                                    nococonv, enpara, lapw, vtot, fmpi, zmat_p, abc_p, &
                                    ev_list=ev_list, kpts=kpts)
         irec = MERGE(1, jspin, noco%l_noco)
         addnoco = 0
         IF (noco%l_noco .AND. jspin == 2) addnoco = lapw%nv(1) + atoms%nlotot

         CALL repack_coefficients(atoms, sym, noco, radfun, abc_p(jspin, :), jspin_rad, nb, &
                                  ff, gg, flo, acof, bcof, ccof)

         IF (input%film) THEN
            IF (ALLOCATED(ac)) DEALLOCATE (ac, bc, u, ue)
            ALLOCATE (ac(nv2d, nb, 2), bc(nv2d, nb, 2), u(vacuum%nmzd, nv2d, vacuum%nvac), &
                      ue(vacuum%nmzd, nv2d, vacuum%nvac))
            CALL wann_2dvacabcof(nv2d, nb, vacuum%nvac, vacuum%nmzd, vacuum%nmz, cell%omtil, vz, &
                                 lapw%nv(jspin), lapw%bkpt, cell%z1, lapw%dim_nvd(), &
                                 lapw%k1(:, jspin), lapw%k2(:, jspin), lapw%k3(:, jspin), &
                                 enpara%evac0(:, jspin), cell%bbmat, vacuum%delz, cell%bmat, &
                                 lapw%dim_nbasfcn(), input%neig, zmat_p(irec), &
                                 ac, bc, u, ue, addnoco, noco%l_ss, nococonv%qss, jspin)
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
                                       wf(:, posi) = wf(:, posi) + xdnout*vgauge(nbn, :)*factor
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
                                           cell%bbmat, enpara%evac0(:, jspin), lapw%bkpt, vz, jspin, &
                                           lapw%k1(:, jspin), lapw%k2(:, jspin), lapw%k3(:, jspin), &
                                           lapw%dim_nvd(), lapw%dim_nbasfcn(), input%neig, &
                                           lapw%nv(jspin), cell%omtil, nb, &
                                           ac(:, nbn, ivac), bc(:, nbn, ivac), &
                                           u(:, :, jvac), ue(:, :, jvac), xdnout)
                        wf(:, posi) = wf(:, posi) + xdnout*vgauge(nbn, :)*factor
                        CYCLE xloop
                     END IF

                     !> interstitial
                     CALL eval_state(point, 0, 0, 0, 2)
                     wf(:, posi) = wf(:, posi) + xdnout*vgauge(nbn, :)*factor
                  END DO xloop
               END DO
            END DO
         END DO
      END DO

      IF (ALLOCATED(ac)) DEALLOCATE (ac, bc, u, ue)

#ifdef CPP_MPI
      CALL collect_wf(wf, fmpi)
#endif
      wf = wf/REAL(kpts%nkptf)

      IF (fmpi%irank == 0) CALL write_xsf(wf, atoms, cell, input, twodim, grid, &
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
                        zmat_p(irec), lapw%nv(jspin), &
                        lapw%k1(:, jspin), lapw%k2(:, jspin), lapw%k3(:, jspin), &
                        atoms%lmaxd*(atoms%lmaxd + 2), lapw%dim_nbasfcn(), &
                        noco%l_ss, nococonv%qss, jspin, addnoco, xdnout)
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
   SUBROUTINE plot_region(cell, twodim, grid, vec1, vec2, vec3, zero)
      TYPE(t_cell), INTENT(IN)  :: cell
      LOGICAL, INTENT(OUT) :: twodim
      INTEGER, INTENT(OUT) :: grid(3)
      REAL,    INTENT(OUT) :: vec1(3), vec2(3), vec3(3), zero(3)

      twodim = .FALSE.
      grid   = [40, 40, 40]
      vec1   = cell%amat(:, 1)
      vec2   = cell%amat(:, 2)
      vec3   = cell%amat(:, 3)
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

      INTEGER :: lmd, nlod, llod, ntyp, na, nat_local, l, m, lm, ilo, ir, nseen

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
         !> A global `ilo` is the nseen-th local orbital of its own l, so its slot is 2+nseen.
         DO ilo = 1, atoms%nlo(ntyp)
            l = atoms%llo(ilo, ntyp)
            nseen = COUNT(atoms%llo(1:ilo, ntyp) == l)
            DO ir = 1, atoms%jri(ntyp)
               flo(ntyp, ir, 1:2, ilo) = radfun(ntyp)%r(ir, 1:2, 2 + nseen, l, jspin_rad)
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
               nseen = COUNT(atoms%llo(1:ilo, ntyp) == l)
               DO m = -l, l
                  lm = l*(l + 1) + m
                  ccof(m, :, ilo, na) = abc(ntyp)%cof(1:nbands, lm, 2 + nseen, nat_local)
               END DO
            END DO
         END DO
      END DO

      !> An atom that abcof reached through its inversion partner carries coefficients in
      !> the rotated frame; the plot needs the global one. The guard is the producer's own
      !> condition (t_abc: any(invsat==2) .and. .not.l_soc) -- with SOC the partners were
      !> built directly and flipping them here would be the error, not the fix.
      IF (ANY(sym%invsat == 2) .AND. .NOT. noco%l_soc) CALL wann_abinv(atoms, sym, acof, bcof, ccof)
   END SUBROUTINE repack_coefficients

#ifdef CPP_MPI
   !> Each rank summed its own share of the k-mesh; the whole sum lives on rank 0.
   SUBROUTINE collect_wf(wf, fmpi)
      USE mpi
      COMPLEX,     INTENT(INOUT) :: wf(:, :)
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      INTEGER :: ierr
      IF (fmpi%isize == 1) RETURN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, wf, SIZE(wf), MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
   END SUBROUTINE collect_wf
#endif

   !> One XSF per Wannier function, two data blocks in it: the real part -- which is what
   !> a two-sign isosurface plot is made of -- and the modulus. The imaginary part is not
   !> written but its weight is reported, because that number is the check on the phase
   !> fix below and a file nobody opens is not.
   SUBROUTINE write_xsf(wf, atoms, cell, input, twodim, grid, vec1, vec2, vec3, zero, jspin)
      COMPLEX,       INTENT(INOUT) :: wf(:, :)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_input), INTENT(IN) :: input
      LOGICAL,       INTENT(IN) :: twodim
      INTEGER,       INTENT(IN) :: grid(3), jspin
      REAL,          INTENT(IN) :: vec1(3), vec2(3), vec3(3), zero(3)

      INTEGER :: nw, mesh, nplo, ipos, imax, iunit, ix, iy, iz
      REAL    :: total, edge
      REAL    :: reweight, imweight
      COMPLEX :: phase
      CHARACTER(LEN=40) :: fname

      nw   = SIZE(wf, 1)
      mesh = SIZE(wf, 2)

      WRITE (oUnit, '(/,a)') 'Wannier functions on a real-space grid (plotWF)'
      WRITE (oUnit, '(a,3i5)') '  grid           ', grid
      WRITE (oUnit, '(a,i5)')  '  functions      ', nw
      WRITE (oUnit, '(a)')     '   WF    max|Im|/max|Re|   weight on the box edge'
      WRITE (oUnit, '(a)')     '                             (both should be small; a large edge'
      WRITE (oUnit, '(a)')     '                              weight means the cell clips the function)'

      DO nplo = 1, nw
         !> A Wannier function is real up to one global phase; take it off, using the point
         !> where the function is largest so the phase is well determined there.
         imax = 1
         DO ipos = 1, mesh
            IF (ABS(wf(nplo, ipos)) > ABS(wf(nplo, imax))) imax = ipos
         END DO
         phase = wf(nplo, imax)/ABS(wf(nplo, imax))
         wf(nplo, :) = wf(nplo, :)/phase

         reweight = MAXVAL(ABS(REAL(wf(nplo, :))))
         imweight = MAXVAL(ABS(AIMAG(wf(nplo, :))))
         IF (reweight > 0.0) imweight = imweight/reweight

         !> How much of the function sits on the outermost layer of the box. If this is not
         !> small the unit cell is too tight for this Wannier function and the picture is cut.
         total = 0.0; edge = 0.0
         DO iz = 0, grid(3) - 1
            DO iy = 0, grid(2) - 1
               DO ix = 0, grid(1) - 1
                  ipos = ix + 1 + iy*grid(1) + iz*grid(1)*grid(2)
                  total = total + ABS(wf(nplo, ipos))**2
                  IF (ix == 0 .OR. ix == grid(1) - 1 .OR. iy == 0 .OR. iy == grid(2) - 1 .OR. &
                      (.NOT. twodim .AND. (iz == 0 .OR. iz == grid(3) - 1))) &
                     edge = edge + ABS(wf(nplo, ipos))**2
               END DO
            END DO
         END DO
         IF (total > 0.0) edge = edge/total
         WRITE (oUnit, '(i5,3x,es12.4,6x,es12.4)') nplo, imweight, edge

         WRITE (fname, '(a,i1,a,i3.3,a)') 'WF', jspin, '_', nplo, '.xsf'
         OPEN (NEWUNIT=iunit, file=TRIM(fname))
         CALL xsf_WRITE_atoms(iunit, atoms, input%film, cell%amat)
         CALL xsf_WRITE_header(iunit, twodim, TRIM(fname), vec1, vec2, vec3, zero, grid)
         DO ipos = 1, mesh
            WRITE (iunit, *) REAL(wf(nplo, ipos))
         END DO
         CALL xsf_WRITE_newblock(iunit, twodim, vec1, vec2, vec3, zero, grid)
         DO ipos = 1, mesh
            WRITE (iunit, *) ABS(wf(nplo, ipos))
         END DO
         CALL xsf_WRITE_endblock(iunit, twodim)
         CLOSE (iunit)
      END DO

      WRITE (oUnit, '(a,i0,a)') '  wrote ', nw, ' XSF files (block A: real part, block B: modulus)'
   END SUBROUTINE write_xsf

END MODULE m_wannierlib_plot
