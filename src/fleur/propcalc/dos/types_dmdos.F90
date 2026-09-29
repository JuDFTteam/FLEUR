!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dmdos
   !! Band- and k-resolved local density matrix of the DOS atoms.
   !! For eigenstate nu at k, block (l,l') and spin block (s,s') the stored matrix is
   !!    rho(m,mp) = sum_ij cof_s(nu,l'mp,i) conjg(cof_s'(nu,lm,j)) <u_i^{l's}|u_j^{ls'}>
   !! which is the n_mmp convention: rho(m,mp) = <l'mp s|rho|lm s'>.
   !! The spin blocks are 11 (and 22) collinear, and 11,22,21,12 in the non-collinear case.
   !! Symmetrized matrices are averaged over the site symmetry and the equivalent atoms,
   !! otherwise every selected atom is stored in the frame of its muffin-tin functions.
#ifdef CPP_MPI
   use mpi
#endif
   use m_judft
   use m_constants
   use m_types_eigdos
   implicit none
   PRIVATE
   public t_dmdos

   TYPE, extends(t_eigdos):: t_dmdos
      INTEGER :: nsb = 0            !number of spin blocks
      LOGICAL :: l_sym = .FALSE.
      LOGICAL :: l_global_frame = .FALSE.
      INTEGER, ALLOCATABLE :: lpairs(:, :)     !(2,npair) blocks (l,l'), l<=l'
      INTEGER, ALLOCATABLE :: elem_atom(:)     !atom of each element (first selected atom of the type if symmetrized)
      INTEGER, ALLOCATABLE :: elem_type(:)
      INTEGER, ALLOCATABLE :: elem_neq(:)
      REAL, ALLOCATABLE    :: elem_euler(:, :) !(3,nelem) orbital frame applied in postprocessing
      COMPLEX, ALLOCATABLE :: rho(:, :, :, :, :, :, :) !(-3:3,-3:3,nsb,neig,nkpt,npair,nelem)
      COMPLEX, ALLOCATABLE :: dwgn(:, :, :, :) !(-3:3,-3:3,0:3,nop) Wigner matrices of the symmetry operations
      COMPLEX, ALLOCATABLE :: sphase(:)        !(nop) spin off-diagonal phase of the operations
      CHARACTER(len=20), ALLOCATABLE :: weight_names(:)
      INTEGER, ALLOCATABLE :: weight_index(:, :) !(3,nweights) element, pair, component
   CONTAINS
      PROCEDURE, PASS :: init => dmdos_init
      PROCEDURE :: calc_dm
      PROCEDURE :: collect
      PROCEDURE :: get_weight_eig
      PROCEDURE :: get_num_weights
      PROCEDURE :: get_weight_name
      PROCEDURE :: sym_weights
      PROCEDURE :: postprocessing
      PROCEDURE :: write_extra
   END TYPE t_dmdos

CONTAINS

   SUBROUTINE dmdos_init(this, input, atoms, kpts, banddos, noco, sym, cell, eig)
      USE m_types_input
      USE m_types_atoms
      USE m_types_kpts
      USE m_types_banddos
      USE m_types_noco
      USE m_types_sym
      USE m_types_cell
      USE m_dwigner
      USE m_angles
      CLASS(t_dmdos), INTENT(INOUT) :: this
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_kpts), INTENT(IN)      :: kpts
      TYPE(t_banddos), INTENT(IN)   :: banddos
      TYPE(t_noco), INTENT(IN)      :: noco
      TYPE(t_sym), INTENT(IN)       :: sym
      TYPE(t_cell), INTENT(IN)      :: cell
      REAL, INTENT(IN)              :: eig(:, :, :)

      character, parameter :: spdf(0:3) = ["s", "p", "d", "f"]
      character(len=3), parameter :: comp_name(4) = ["   ", "mx ", "my ", "mz "]
      INTEGER :: n, p, l, iop, ncomp, ic, nw
      COMPLEX :: d(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, lmaxU_const, sym%nop)
      TYPE(t_sym) :: sym_phase

      this%l_initialized = .TRUE.
      this%name_of_dos = "DensityMatrix"
      this%nsb = merge(4, input%jspins, noco%l_noco)
      this%l_global_frame = noco%l_noco .AND. banddos%global_frame
      SELECT CASE (banddos%dm_sym)
      CASE (1)
         this%l_sym = .TRUE.
      CASE (2)
         this%l_sym = .FALSE.
      CASE DEFAULT
         this%l_sym = kpts%kptsKind /= KPTS_KIND_PATH
      END SELECT
      this%lpairs = banddos%dm_lpairs

      IF (this%l_sym) THEN
         this%elem_type = banddos%dos_typelist
         allocate (this%elem_atom(size(this%elem_type)))
         DO n = 1, size(this%elem_type)
            this%elem_atom(n) = minval(banddos%dos_atomlist, mask=atoms%itype(banddos%dos_atomlist) == this%elem_type(n))
         END DO
      ELSE
         this%elem_atom = banddos%dos_atomlist
         this%elem_type = atoms%itype(banddos%dos_atomlist)
      END IF
      this%elem_neq = atoms%neq(this%elem_type)
      allocate (this%elem_euler(3, size(this%elem_atom)), source=0.0)

      allocate (this%rho(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, this%nsb, input%neig, kpts%nkpt, &
                         size(this%lpairs, 2), size(this%elem_atom)), source=cmplx_0)
      IF (noco%l_noco) THEN
         this%eig = eig(:, :, 1:1)
      ELSE
         this%eig = eig
      END IF

      IF (this%l_sym) THEN
         call d_wigner(sym%nop, sym%mrot, cell%bmat, lmaxU_const, d)
         allocate (this%dwgn(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, 0:lmaxU_const, sym%nop), source=cmplx_0)
         DO iop = 1, sym%nop
            this%dwgn(0, 0, 0, iop) = cmplx_1
            DO l = 1, lmaxU_const
               this%dwgn(-l:l, -l:l, l, iop) = d(-l:l, -l:l, l, iop)
            END DO
         END DO
         allocate (this%sphase(sym%nop), source=cmplx_1)
         IF (noco%l_noco) THEN
            sym_phase = sym
            IF (allocated(sym_phase%phase)) deallocate (sym_phase%phase)
            allocate (sym_phase%phase(sym%nop))
            call angles(sym_phase)
            DO iop = 1, sym%nop
               this%sphase(iop) = exp(ImagUnit*sym_phase%phase(sym%invtab(iop)))
            END DO
         END IF
      END IF

      !scalar weights: trace per spin (collinear) or charge and magnetization (non-collinear)
      ncomp = merge(4, 1, noco%l_noco)
      nw = count(this%lpairs(1, :) == this%lpairs(2, :))*size(this%elem_atom)*ncomp
      allocate (this%weight_names(nw), this%weight_index(3, nw))
      nw = 0
      DO n = 1, size(this%elem_atom)
         DO p = 1, size(this%lpairs, 2)
            l = this%lpairs(1, p)
            IF (l /= this%lpairs(2, p)) CYCLE
            DO ic = 1, ncomp
               nw = nw + 1
               write (this%weight_names(nw), "(a,i0,a,a)") "DM:", this%elem_atom(n), spdf(l), trim(comp_name(ic))
               this%weight_index(:, nw) = [n, p, ic]
            END DO
         END DO
      END DO
   END SUBROUTINE dmdos_init

   SUBROUTINE calc_dm(this, abc_s, abc_sp, radfun, atoms, sym, ev_list, itype, ikpt, ispin, ispinpr)
      !! Adds the contribution of atom type itype for the spin pair (ispin,ispinpr) with ispinpr<=ispin.
      !! For ispin/=ispinpr both off-diagonal blocks 21 and 12 are computed.
      USE m_types_atoms
      USE m_types_sym
      USE m_types_abc
      USE m_types_radfun
      USE m_intgr
      CLASS(t_dmdos), INTENT(INOUT) :: this
      TYPE(t_abc), INTENT(IN)       :: abc_s, abc_sp
      TYPE(t_radfun), INTENT(IN)    :: radfun
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_sym), INTENT(IN)       :: sym
      INTEGER, INTENT(IN)           :: ev_list(:)
      INTEGER, INTENT(IN)           :: itype, ikpt, ispin, ispinpr

      INTEGER :: natom, na, n, p, l, lp, ib, s, sp, b, i, m, mp, j, jj, nr, nrp, jri
      REAL    :: ovl(size(radfun%integral, 1), size(radfun%integral, 1))
      REAL, ALLOCATABLE :: rf(:)
      COMPLEX :: mat(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const), bvec(size(abc_s%cof, 3))

      if (.not. this%l_initialized) return
      jri = atoms%jri(itype)
      allocate (rf(jri))

      DO natom = 1, atoms%neq(itype)
         na = atoms%firstAtom(itype) + natom - 1
         IF (this%l_sym) THEN
            n = findloc(this%elem_type, itype, dim=1)
         ELSE
            n = findloc(this%elem_atom, na, dim=1)
         END IF
         IF (n == 0) CYCLE
         DO p = 1, size(this%lpairs, 2)
            l = this%lpairs(1, p); lp = this%lpairs(2, p)
            IF (lp > atoms%lmax(itype)) CYCLE
            nr = radfun%n_r(l); nrp = radfun%n_r(lp)
            DO ib = 1, merge(1, 2, ispin == ispinpr)
               !column (l') belongs to spin s, row (l) to spin sp
               IF (ib == 1) THEN
                  s = ispin; sp = ispinpr
               ELSE
                  s = ispinpr; sp = ispin
               END IF
               b = block_index(this, s, sp)
               IF (l == lp) THEN
                  ovl(:nrp, :nr) = radfun%integral(:nrp, :nr, l, s, sp)
               ELSE
                  DO j = 1, nrp
                     DO jj = 1, nr
                        rf = radfun%r(:jri, 1, j, lp, s)*radfun%r(:jri, 1, jj, l, sp) + &
                             radfun%r(:jri, 2, j, lp, s)*radfun%r(:jri, 2, jj, l, sp)
                        CALL intgr0(rf, atoms%rmsh(1, itype), atoms%dx(itype), jri, ovl(j, jj))
                     END DO
                  END DO
               END IF
               DO i = 1, size(ev_list)
                  mat = cmplx_0
                  DO m = -l, l
                     IF (ib == 1) THEN
                        bvec(:nrp) = matmul(ovl(:nrp, :nr), conjg(abc_sp%cof(i, l*(l + 1) + m, :nr, natom)))
                     ELSE
                        bvec(:nrp) = matmul(ovl(:nrp, :nr), conjg(abc_s%cof(i, l*(l + 1) + m, :nr, natom)))
                     END IF
                     DO mp = -lp, lp
                        IF (ib == 1) THEN
                           mat(m, mp) = sum(abc_s%cof(i, lp*(lp + 1) + mp, :nrp, natom)*bvec(:nrp))
                        ELSE
                           mat(m, mp) = sum(abc_sp%cof(i, lp*(lp + 1) + mp, :nrp, natom)*bvec(:nrp))
                        END IF
                     END DO
                  END DO
                  IF (this%l_sym) THEN
                     this%rho(:, :, b, ev_list(i), ikpt, p, n) = this%rho(:, :, b, ev_list(i), ikpt, p, n) + &
                        site_symmetrize(this, mat, sym, na, l, lp, s /= sp)/atoms%neq(itype)
                  ELSE
                     this%rho(:, :, b, ev_list(i), ikpt, p, n) = mat
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE calc_dm

   INTEGER FUNCTION block_index(this, s, sp)
      CLASS(t_dmdos), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: s, sp
      IF (s == sp) THEN
         block_index = s
      ELSE IF (s == 2) THEN
         block_index = 3
      ELSE
         block_index = 4
      END IF
   END FUNCTION block_index

   FUNCTION site_symmetrize(this, mat, sym, na, l, lp, l_offdiag) RESULT(mat_sym)
      !! Average over the operations that map atom na onto itself, as symMMPmat (also for l=0 blocks)
      USE m_types_sym
      CLASS(t_dmdos), INTENT(IN) :: this
      COMPLEX, INTENT(IN)        :: mat(-lmaxU_const:, -lmaxU_const:)
      TYPE(t_sym), INTENT(IN)    :: sym
      INTEGER, INTENT(IN)        :: na, l, lp
      LOGICAL, INTENT(IN)        :: l_offdiag
      COMPLEX :: mat_sym(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const)

      INTEGER :: iop, nops
      COMPLEX :: ph

      mat_sym = cmplx_0
      nops = 0
      DO iop = 1, sym%nop
         IF (sym%mapped_atom(iop, na) /= na) CYCLE
         nops = nops + 1
         ph = cmplx_1
         IF (l_offdiag) ph = this%sphase(iop)
         mat_sym = mat_sym + ph*matmul(conjg(transpose(this%dwgn(:, :, l, iop))), matmul(mat, this%dwgn(:, :, lp, iop)))
      END DO
      mat_sym = mat_sym/nops
   END FUNCTION site_symmetrize

   SUBROUTINE collect(this, fmpi, jspin)
      !! Sum the data of all ranks on rank 0. In a non-collinear run all blocks are collected once (jspin=1),
      !! otherwise the block of the current spin.
      USE m_types_mpi
      CLASS(t_dmdos), INTENT(INOUT) :: this
      TYPE(t_mpi), INTENT(IN)       :: fmpi
      INTEGER, INTENT(IN)           :: jspin
#ifdef CPP_MPI
      INTEGER, PARAMETER :: max_buf = 4000000
      INTEGER :: bs, be, n, p, k0, k1, nk, ierr
      COMPLEX, ALLOCATABLE :: sbuf(:, :, :, :, :), rbuf(:, :, :, :, :)

      if (.not. this%l_initialized) return
      IF (this%nsb == 4) THEN
         IF (jspin /= 1) RETURN
         bs = 1; be = 4
      ELSE
         bs = jspin; be = jspin
      END IF
      nk = max(1, max_buf/(size(this%rho, 1)*size(this%rho, 2)*(be - bs + 1)*size(this%rho, 4)))
      DO n = 1, size(this%rho, 7)
         DO p = 1, size(this%rho, 6)
            DO k0 = 1, size(this%rho, 5), nk
               k1 = min(k0 + nk - 1, size(this%rho, 5))
               sbuf = this%rho(:, :, bs:be, :, k0:k1, p, n)
               allocate (rbuf, mold=sbuf)
               CALL MPI_REDUCE(sbuf, rbuf, size(sbuf), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr)
               IF (fmpi%irank == 0) this%rho(:, :, bs:be, :, k0:k1, p, n) = rbuf
               deallocate (rbuf)
            END DO
         END DO
      END DO
#endif
   END SUBROUTINE collect

   SUBROUTINE sym_weights(this)
      !! Average the matrices over degenerate states
      CLASS(t_dmdos), INTENT(INOUT) :: this
      INTEGER :: b, ikpt, i, j, s, ne

      ne = min(size(this%eig, 1), size(this%rho, 4))
      DO b = 1, this%nsb
         s = merge(1, b, this%nsb == 4)
         DO ikpt = 1, size(this%rho, 5)
            i = 1
            DO WHILE (i < ne)
               j = 0
               DO WHILE (i + j + 1 <= ne)
                  IF (abs(this%eig(i, ikpt, s) - this%eig(i + j + 1, ikpt, s)) >= 1E-5) EXIT
                  j = j + 1
               END DO
               IF (j > 0) THEN
                  this%rho(:, :, b, i, ikpt, :, :) = sum(this%rho(:, :, b, i:i + j, ikpt, :, :), dim=3)/(j + 1)
                  this%rho(:, :, b, i + 1:i + j, ikpt, :, :) = spread(this%rho(:, :, b, i, ikpt, :, :), 3, j)
               END IF
               i = i + j + 1
            END DO
         END DO
      END DO
   END SUBROUTINE sym_weights

   SUBROUTINE postprocessing(this, noco, nococonv, banddos, alldos, ef)
      !! Rotates into the requested orbital and spin frames and checks the traces against the l-resolved DOS
      USE m_types_noco
      USE m_types_nococonv
      USE m_types_banddos
      USE m_types_dos
      USE m_dwigner
      CLASS(t_dmdos), INTENT(INOUT) :: this
      TYPE(t_noco), INTENT(IN)      :: noco
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      TYPE(t_banddos), INTENT(IN)   :: banddos
      CLASS(t_eigdos_list), INTENT(IN), OPTIONAL :: alldos(:)
      REAL, INTENT(IN), OPTIONAL    :: ef

      INTEGER :: n, na, itype, p, l, lp, b, i, ikpt, m, mp
      REAL    :: amx(3, 3, 1), imx(3, 3), angle(3), co, si
      COMPLEX :: dmat(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, 0:lmaxU_const)
      COMPLEX :: d(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, lmaxU_const, 1)
      COMPLEX :: u(2, 2), s2(2, 2)
      TYPE(t_dos), POINTER :: dos

      WRITE (oUnit, "(/,a,f10.3,a)") " DensityMatrix DOS: storage of the matrices ", &
         16.0*size(this%rho, kind=8)/1024.0**3, " GB"

      !orbital frame
      DO n = 1, size(this%elem_atom)
         na = this%elem_atom(n); itype = this%elem_type(n)
         IF (banddos%align_to_spin(na)) THEN
            angle = [0.0, nococonv%beta(itype), 0.0]
            IF (abs(angle(2)) > 1.0e-6) angle(3) = pi_const/2.0 - nococonv%alph(itype)
         ELSE
            angle = [banddos%alpha(na), banddos%beta(na), banddos%gamma(na)]
         END IF
         this%elem_euler(:, n) = angle
         IF (all(abs(angle) < 1.0e-12)) CYCLE
         CALL euler(angle(1), angle(2), angle(3), amx)
         imx = 0.0; imx(1, 1) = 1.0; imx(2, 2) = 1.0; imx(3, 3) = 1.0
         CALL d_wigner(1, amx, imx, lmaxU_const, d)
         dmat = cmplx_0
         dmat(0, 0, 0) = cmplx_1
         DO l = 1, lmaxU_const
            dmat(-l:l, -l:l, l) = d(-l:l, -l:l, l, 1)
         END DO
         !the coefficients transform with conjg(D), hence rho -> D_l rho D_l'^+
         DO p = 1, size(this%lpairs, 2)
            l = this%lpairs(1, p); lp = this%lpairs(2, p)
            DO ikpt = 1, size(this%rho, 5)
               DO i = 1, size(this%rho, 4)
                  DO b = 1, this%nsb
                     this%rho(:, :, b, i, ikpt, p, n) = matmul(dmat(:, :, l), &
                        matmul(this%rho(:, :, b, i, ikpt, p, n), conjg(transpose(dmat(:, :, lp)))))
                  END DO
               END DO
            END DO
         END DO
      END DO

      !spin frame: same transformation as rotdenmat(toGlobal) of t_dos
      IF (this%l_global_frame) THEN
         DO n = 1, size(this%elem_atom)
            itype = this%elem_type(n)
            co = cos(0.5*nococonv%beta(itype)); si = sin(0.5*nococonv%beta(itype))
            u(1, :) = [co, -si]
            u(2, :) = [si, co]
            u(1, :) = u(1, :)*exp(-0.5*ImagUnit*nococonv%alph(itype))
            u(2, :) = u(2, :)*exp(0.5*ImagUnit*nococonv%alph(itype))
            DO p = 1, size(this%lpairs, 2)
               DO ikpt = 1, size(this%rho, 5)
                  DO i = 1, size(this%rho, 4)
                     DO mp = -lmaxU_const, lmaxU_const
                        DO m = -lmaxU_const, lmaxU_const
                           s2(1, 1) = this%rho(m, mp, 1, i, ikpt, p, n)
                           s2(2, 2) = this%rho(m, mp, 2, i, ikpt, p, n)
                           s2(2, 1) = this%rho(m, mp, 3, i, ikpt, p, n)
                           s2(1, 2) = this%rho(m, mp, 4, i, ikpt, p, n)
                           s2 = matmul(u, matmul(s2, conjg(transpose(u))))
                           this%rho(m, mp, 1, i, ikpt, p, n) = s2(1, 1)
                           this%rho(m, mp, 2, i, ikpt, p, n) = s2(2, 2)
                           this%rho(m, mp, 3, i, ikpt, p, n) = s2(2, 1)
                           this%rho(m, mp, 4, i, ikpt, p, n) = s2(1, 2)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF

      !consistency check against the l-resolved weights of t_dos
      IF (.NOT. present(alldos)) RETURN
      NULLIFY (dos)
      DO n = 1, size(alldos)
         IF (.not. associated(alldos(n)%p)) CYCLE
         SELECT TYPE (dd => alldos(n)%p)
         TYPE IS (t_dos)
            dos => dd
            EXIT
         END SELECT
      END DO
      IF (.NOT. associated(dos)) RETURN
      WRITE (oUnit, "(a,es12.3)") " DensityMatrix DOS: max. deviation of the traces from the l-resolved DOS weights:", &
         trace_deviation(this, dos, banddos)
   END SUBROUTINE postprocessing

   REAL FUNCTION trace_deviation(this, dos, banddos)
      !! Largest difference between Tr rho_ll and the corresponding MT:<type><l> weight. Unsymmetrized data is
      !! compared as average over the equivalent atoms, the spin off-diagonal part only without symmetrization.
      USE m_types_dos
      USE m_types_banddos
      CLASS(t_dmdos), INTENT(IN)  :: this
      TYPE(t_dos), INTENT(IN)     :: dos
      TYPE(t_banddos), INTENT(IN) :: banddos

      INTEGER :: it, itype, n_dos, p, l, b, m, n, cnt, nb, ne, nk
      REAL, ALLOCATABLE :: tr(:, :, :)

      trace_deviation = 0.0
      ne = min(size(this%rho, 4), size(dos%qal, 3))
      nk = size(this%rho, 5)
      nb = merge(this%nsb, min(this%nsb, 2), .NOT. this%l_sym)
      allocate (tr(ne, nk, 4))
      DO it = 1, size(banddos%dos_typelist)
         itype = banddos%dos_typelist(it)
         n_dos = banddos%map_atomtype(itype)
         DO p = 1, size(this%lpairs, 2)
            l = this%lpairs(1, p)
            IF (l /= this%lpairs(2, p)) CYCLE
            tr = 0.0
            cnt = 0
            DO n = 1, size(this%elem_atom)
               IF (this%elem_type(n) /= itype) CYCLE
               cnt = cnt + 1
               DO b = 1, nb
                  DO m = -l, l
                     IF (b == 4) THEN
                        tr(:, :, 4) = tr(:, :, 4) + aimag(this%rho(m, m, 3, :ne, :, p, n))
                     ELSE
                        tr(:, :, b) = tr(:, :, b) + real(this%rho(m, m, b, :ne, :, p, n))
                     END IF
                  END DO
               END DO
            END DO
            IF (cnt == 0) CYCLE
            IF (.NOT. this%l_sym .AND. cnt /= this%elem_neq(findloc(this%elem_type, itype, dim=1))) CYCLE
            DO b = 1, nb
               trace_deviation = max(trace_deviation, maxval(abs(tr(:, :, b)/cnt - dos%qal(l, n_dos, :ne, :nk, b))))
            END DO
         END DO
      END DO
   END FUNCTION trace_deviation

   INTEGER FUNCTION get_num_weights(this)
      CLASS(t_dmdos), INTENT(IN) :: this
      get_num_weights = 0
      IF (allocated(this%weight_names)) get_num_weights = size(this%weight_names)
   END FUNCTION get_num_weights

   CHARACTER(len=20) FUNCTION get_weight_name(this, id)
      CLASS(t_dmdos), INTENT(IN) :: this
      INTEGER, INTENT(IN)        :: id
      IF (id > get_num_weights(this)) CALL judft_error("Not enough weight names in t_dmdos")
      get_weight_name = this%weight_names(id)
   END FUNCTION get_weight_name

   FUNCTION get_weight_eig(this, id)
      !! Traces of the diagonal blocks: per spin (collinear) or charge, mx, my, mz (non-collinear)
      CLASS(t_dmdos), INTENT(IN) :: this
      INTEGER, INTENT(IN)        :: id
      REAL, ALLOCATABLE          :: get_weight_eig(:, :, :)

      INTEGER :: n, p, ic, l, m, b
      COMPLEX, ALLOCATABLE :: tr(:, :, :)

      n = this%weight_index(1, id); p = this%weight_index(2, id); ic = this%weight_index(3, id)
      l = this%lpairs(1, p)
      allocate (tr(size(this%rho, 4), size(this%rho, 5), this%nsb), source=cmplx_0)
      DO b = 1, this%nsb
         DO m = -l, l
            tr(:, :, b) = tr(:, :, b) + this%rho(m, m, b, :, :, p, n)
         END DO
      END DO
      IF (this%nsb < 4) THEN
         get_weight_eig = real(tr)
      ELSE
         allocate (get_weight_eig(size(tr, 1), size(tr, 2), 1))
         SELECT CASE (ic)
         CASE (1)
            get_weight_eig(:, :, 1) = real(tr(:, :, 1) + tr(:, :, 2))
         CASE (2)
            get_weight_eig(:, :, 1) = 2.0*real(tr(:, :, 3))
         CASE (3)
            get_weight_eig(:, :, 1) = 2.0*aimag(tr(:, :, 3))
         CASE (4)
            get_weight_eig(:, :, 1) = real(tr(:, :, 1) - tr(:, :, 2))
         CASE DEFAULT
            CALL judft_bug("t_dmdos: invalid weight component")
         END SELECT
      END IF
   END FUNCTION get_weight_eig

   SUBROUTINE write_extra(this, hdf_id)
#ifdef CPP_HDF
      USE HDF5
      USE m_banddos_io
#endif
      CLASS(t_dmdos), INTENT(IN) :: this
#ifdef CPP_HDF
      INTEGER(HID_T), INTENT(IN) :: hdf_id
#else
      INTEGER, INTENT(IN)        :: hdf_id
#endif
#ifdef CPP_HDF
      CALL writeDensityMatrixData(hdf_id, this%name_of_dos, this%lpairs, this%elem_atom, this%elem_type, this%l_sym, &
                                  this%l_global_frame, this%elem_euler, this%rho, this%eig)
#endif
   END SUBROUTINE write_extra

END MODULE m_types_dmdos
