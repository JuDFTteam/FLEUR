!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfprep_layer
   !Construction of the film setup of one layer cut out of the bulk
   !supercell: atoms, symmetry (the subgroup of the bulk operations
   !compatible with the film), atom types and the common k-point set.
   USE m_juDFT
   USE m_constants
   USE m_types_atoms
   USE m_types_cell
   USE m_types_sym
   USE m_types_input
   USE m_types_vacuum
   USE m_types_kpts
   USE m_types_enparaXML
   USE m_gfprep_inp, ONLY: t_layergeo
   USE m_gfprep_bulk, ONLY: t_gfprep_struct
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_gfprep_layer, gfprep_make_layer, gfprep_kpts

   TYPE t_gfprep_layer
      TYPE(t_layergeo)  :: geo
      !film setup as written to the layer inp.xml
      TYPE(t_atoms)     :: atoms
      TYPE(t_enparaXML) :: enpara
      TYPE(t_sym)       :: sym
      TYPE(t_cell)      :: cell
      TYPE(t_vacuum)    :: vacuum
      TYPE(t_input)     :: input
      !layer type -> bulk type, layer atom -> bulk atom
      INTEGER, ALLOCATABLE :: typmap(:), atmap(:)
   END TYPE t_gfprep_layer

   REAL, PARAMETER :: eps = 1.0E-5

CONTAINS

   SUBROUTINE gfprep_make_layer(b, geo, gmax_scale, lay)
      TYPE(t_gfprep_struct), INTENT(IN)   :: b
      TYPE(t_layergeo), INTENT(IN)        :: geo
      REAL, INTENT(IN)                    :: gmax_scale  !Gmax(layer)/Gmax(bulk)
      TYPE(t_gfprep_layer), INTENT(INOUT) :: lay

      !candidate atoms: in-plane internal coordinates, z relative to the
      !layer center (cartesian), bulk atom and bulk type
      REAL, ALLOCATABLE    :: xy(:, :), z(:)
      INTEGER, ALLOCATABLE :: nb(:), tb(:), opmap(:, :), itype(:), order(:)
      REAL, ALLOCATABLE    :: trans(:, :)
      INTEGER :: na, i3, nat, n, iop, nop, a, t, ntype, i, k, it
      REAL    :: cz, zr, rmt

      lay%geo = geo
      cz = b%fi%cell%amat(3, 3)

      !atoms whose MT spheres lie completely within |z|<=d/2
      ALLOCATE (xy(2, 5*b%fi%atoms%nat), z(5*b%fi%atoms%nat), nb(5*b%fi%atoms%nat), tb(5*b%fi%atoms%nat))
      nat = 0
      DO na = 1, b%fi%atoms%nat
         n = b%fi%atoms%itype(na)
         rmt = b%fi%atoms%rmt(n)
         DO i3 = FLOOR((geo%zc - geo%d)/cz) - 1, CEILING((geo%zc + geo%d)/cz) + 1
            zr = b%fi%atoms%pos(3, na) + i3*cz - geo%zc
            IF (ABS(zr) + rmt > 0.5*geo%d) CYCLE
            nat = nat + 1
            IF (nat > SIZE(z)) CALL priv_grow(xy, z, nb, tb)
            xy(:, nat) = b%fi%atoms%taual(1:2, na)
            z(nat) = zr
            nb(nat) = na
            tb(nat) = n
         END DO
      END DO
      IF (nat == 0) CALL juDFT_error("No atoms in layer", calledby="gfprep_make_layer")

      !bulk operations, combined with the pure translations of the
      !supercell, that leave the layer invariant (one per rotation)
      CALL priv_pure_translations(b%fi%atoms, trans)
      ALLOCATE (lay%sym%mrot(3, 3, b%fi%sym%nop), lay%sym%tau(3, b%fi%sym%nop), opmap(nat, b%fi%sym%nop))
      nop = 0
      DO iop = 1, b%fi%sym%nop
         DO it = 1, SIZE(trans, 2)
            IF (priv_layer_op(b%fi%sym%mrot(:, :, iop), b%fi%sym%tau(:, iop) + trans(:, it), geo%zc/cz, &
                              xy(:, :nat), z(:nat), tb(:nat), opmap(:, nop + 1))) THEN
               nop = nop + 1
               lay%sym%mrot(:, :, nop) = b%fi%sym%mrot(:, :, iop)
               lay%sym%tau(1:2, nop) = MODULO(b%fi%sym%tau(1:2, iop) + trans(1:2, it) + eps, 1.0) - eps
               lay%sym%tau(3, nop) = 0.0
               EXIT
            END IF
         END DO
      END DO
      lay%sym%nop = nop
      lay%sym%nop2 = COUNT(lay%sym%mrot(3, 3, :nop) == 1)
      lay%sym%mrot = lay%sym%mrot(:, :, :nop)
      lay%sym%tau = lay%sym%tau(:, :nop)

      !atom types: orbits of the layer operations
      ALLOCATE (itype(nat), order(nat))
      itype = 0
      ntype = 0
      k = 0
      DO a = 1, nat
         IF (itype(a) /= 0) CYCLE
         ntype = ntype + 1
         k = k + 1
         order(k) = a
         itype(a) = ntype
         DO iop = 1, nop
            i = opmap(a, iop)
            IF (itype(i) == 0) THEN
               itype(i) = ntype
               k = k + 1
               order(k) = i
            ELSE IF (itype(i) /= ntype) THEN
               CALL juDFT_error("BUG: inconsistent atom orbits", calledby="gfprep_make_layer")
            END IF
         END DO
      END DO

      CALL priv_fill_atoms(b, lay, ntype, itype(order), xy(:, order), z(order), nb(order), tb(order))

      !cell, vacuum and input of the film
      lay%cell%amat = 0.0
      lay%cell%amat(1:2, 1:2) = b%fi%cell%amat(1:2, 1:2)
      lay%cell%amat(3, 3) = geo%dt
      lay%vacuum%dvac = geo%d
      lay%vacuum%nvac = 2
      lay%input = b%fi%input
      lay%input%film = .TRUE.
      lay%input%gmax = gmax_scale*lay%input%gmax
      lay%input%zelec = 0.0
      DO t = 1, ntype
         lay%input%zelec = lay%input%zelec + lay%atoms%neq(t)*lay%atoms%econf(t)%valence_electrons
      END DO

      WRITE (oUnit, "(a,i0,a,i0,a,i0,a)") "gfprep: layer with ", nat, " atoms, ", ntype, &
         " types and ", nop, " symmetry operations"
   END SUBROUTINE gfprep_make_layer

   SUBROUTINE priv_pure_translations(atoms, trans)
      !internal translations (including 0) mapping the bulk supercell
      !onto itself, e.g. (0,0,1/2) for a cell doubled along z
      TYPE(t_atoms), INTENT(IN)      :: atoms
      REAL, ALLOCATABLE, INTENT(OUT) :: trans(:, :)

      REAL    :: t(3), d(3), cand(3, atoms%nat)
      INTEGER :: na, nb, nc, n
      LOGICAL :: l_ok, l_found

      n = 0
      DO na = 1, atoms%nat
         IF (atoms%itype(na) /= atoms%itype(1)) CYCLE
         t = atoms%taual(:, na) - atoms%taual(:, 1)
         t = t - NINT(t)
         l_ok = .TRUE.
         DO nb = 1, atoms%nat
            l_found = .FALSE.
            DO nc = 1, atoms%nat
               IF (atoms%itype(nc) /= atoms%itype(nb)) CYCLE
               d = atoms%taual(:, nb) + t - atoms%taual(:, nc)
               IF (ALL(ABS(d - NINT(d)) < eps)) l_found = .TRUE.
            END DO
            l_ok = l_ok .AND. l_found
         END DO
         IF (.NOT. l_ok) CYCLE
         n = n + 1
         cand(:, n) = t
      END DO
      trans = cand(:, :n)
   END SUBROUTINE priv_pure_translations

   LOGICAL FUNCTION priv_layer_op(mrot, tau, zc_int, xy, z, tb, map)
      !Checks whether the bulk operation {mrot|tau} is a symmetry of the
      !layer centered at zc_int (internal bulk coordinate). map returns
      !the image of each atom.
      INTEGER, INTENT(IN)  :: mrot(3, 3), tb(:)
      REAL, INTENT(IN)     :: tau(3), zc_int, xy(:, :), z(:)
      INTEGER, INTENT(OUT) :: map(:)

      INTEGER :: a, b
      REAL    :: t3, xyn(2), zn, dxy(2)

      priv_layer_op = .FALSE.
      map = 0
      IF (ANY(mrot(1:2, 3) /= 0) .OR. ANY(mrot(3, 1:2) /= 0)) RETURN
      !z-translation in the frame of the layer center
      t3 = tau(3) + (mrot(3, 3) - 1)*zc_int
      IF (ABS(t3 - NINT(t3)) > eps) RETURN
      DO a = 1, SIZE(z)
         xyn = MATMUL(REAL(mrot(1:2, 1:2)), xy(:, a)) + tau(1:2)
         zn = mrot(3, 3)*z(a)
         DO b = 1, SIZE(z)
            IF (tb(b) /= tb(a)) CYCLE
            IF (ABS(z(b) - zn) > eps) CYCLE
            dxy = xyn - xy(:, b)
            IF (ALL(ABS(dxy - NINT(dxy)) < eps)) THEN
               map(a) = b
               EXIT
            END IF
         END DO
         IF (map(a) == 0) RETURN
      END DO
      priv_layer_op = .TRUE.
   END FUNCTION priv_layer_op

   SUBROUTINE priv_fill_atoms(b, lay, ntype, itype, xy, z, nb, tb)
      !Sets up the atoms (and the energy parameter quantum numbers) of the
      !layer; the species data of each layer type is copied from its bulk type
      TYPE(t_gfprep_struct), INTENT(IN)   :: b
      TYPE(t_gfprep_layer), INTENT(INOUT) :: lay
      INTEGER, INTENT(IN) :: ntype, itype(:), nb(:), tb(:)
      REAL, INTENT(IN)    :: xy(:, :), z(:)

      INTEGER :: t, n, a, i, nu, nhia, jsp
      TYPE(t_atoms) :: atb

      atb = b%fi%atoms
      ASSOCIATE (at => lay%atoms)
         at%ntype = ntype
         at%nat = SIZE(z)
         at%nlod = atb%nlod
         at%llod = atb%llod
         ALLOCATE (lay%typmap(ntype), lay%atmap(at%nat))
         DO t = 1, ntype
            lay%typmap(t) = tb(FINDLOC(itype, t, 1))
         END DO
         lay%atmap = nb
         ALLOCATE (at%neq(ntype))
         DO t = 1, ntype
            at%neq(t) = COUNT(itype == t)
         END DO
         ASSOCIATE (m => lay%typmap)
            at%nz = atb%nz(m)
            at%zatom = atb%zatom(m)
            at%rmt = atb%rmt(m)
            at%jri = atb%jri(m)
            at%dx = atb%dx(m)
            at%lmax = atb%lmax(m)
            at%lnonsph = atb%lnonsph(m)
            at%lapw_l = atb%lapw_l(m)
            at%econf = atb%econf(m)
            at%speciesName = atb%speciesName(m)
            at%l_geo = atb%l_geo(m)
            at%relax = atb%relax(:, m)
            at%bmu = atb%bmu(m)
            at%flipSpinPhi = atb%flipSpinPhi(m)
            at%flipSpinTheta = atb%flipSpinTheta(m)
            at%flipSpinScale = atb%flipSpinScale(m)
            at%nlo = atb%nlo(m)
            at%llo = atb%llo(:, m)
            at%ulo_der = atb%ulo_der(:, m)
         END ASSOCIATE
         at%lmaxd = MAXVAL(at%lmax)

         !LDA+U and Hubbard 1: first all U entries, then all HIA entries
         ALLOCATE (at%lda_u(MAX(1, 4*ntype)))
         nu = 0
         DO t = 1, ntype
            DO i = 1, atb%n_u
               IF (atb%lda_u(i)%atomType /= lay%typmap(t)) CYCLE
               nu = nu + 1
               at%lda_u(nu) = atb%lda_u(i)
               at%lda_u(nu)%atomType = t
            END DO
         END DO
         nhia = 0
         DO t = 1, ntype
            DO i = atb%n_u + 1, atb%n_u + atb%n_hia
               IF (atb%lda_u(i)%atomType /= lay%typmap(t)) CYCLE
               nhia = nhia + 1
               at%lda_u(nu + nhia) = atb%lda_u(i)
               at%lda_u(nu + nhia)%atomType = t
            END DO
         END DO
         at%n_u = nu
         at%n_hia = nhia

         !positions: taual(3) is internal with respect to dTilda
         ALLOCATE (at%taual(3, at%nat), at%pos(3, at%nat), at%label(at%nat))
         DO a = 1, at%nat
            at%taual(1:2, a) = xy(:, a)
            at%taual(3, a) = z(a)/lay%geo%dt
            at%label(a) = atb%label(nb(a))
         END DO
         at%pos = 0.0
      END ASSOCIATE

      !quantum numbers of the energy parameters
      n = SIZE(b%fi%enparaXML%qn_el, 3)
      ALLOCATE (lay%enpara%qn_el(0:3, ntype, n), lay%enpara%qn_ello(MAX(1, atb%nlod), ntype, n))
      lay%enpara%qn_ello = 0
      DO jsp = 1, n
         lay%enpara%qn_el(:, :, jsp) = b%fi%enparaXML%qn_el(:, lay%typmap, jsp)
         IF (atb%nlod > 0) lay%enpara%qn_ello(:, :, jsp) = b%fi%enparaXML%qn_ello(:, lay%typmap, jsp)
      END DO
      lay%enpara%evac0 = eVac0Default_const
   END SUBROUTINE priv_fill_atoms

   SUBROUTINE gfprep_kpts(lays, kgrid, l_timerev, kpts)
      !The 2D k-point set shared by all layers: a Gamma-centered grid
      !reduced with the in-plane point operations common to all layers
      TYPE(t_gfprep_layer), INTENT(IN) :: lays(:)
      INTEGER, INTENT(IN)              :: kgrid(2)
      LOGICAL, INTENT(IN)              :: l_timerev
      TYPE(t_kpts), INTENT(OUT)        :: kpts

      INTEGER, ALLOCATABLE :: rot(:, :, :), owner(:, :)
      REAL, ALLOCATABLE    :: bk(:, :), wt(:)
      INTEGER :: nrot, iop, l, i1, i2, j1, j2, ir, nk, is
      REAL    :: k(2), q(2)
      LOGICAL :: l_common

      IF (ANY(kgrid < 1)) CALL juDFT_error("Invalid 2D k-point grid", &
                                           hint="set kgrid in &GFPREP", calledby="gfprep_kpts")
      !in-plane rotations common to all layers
      ALLOCATE (rot(2, 2, lays(1)%sym%nop))
      nrot = 0
      DO iop = 1, lays(1)%sym%nop
         l_common = .TRUE.
         DO l = 2, SIZE(lays)
            l_common = l_common .AND. priv_has_rot(lays(l)%sym, lays(1)%sym%mrot(1:2, 1:2, iop))
         END DO
         IF (.NOT. l_common) CYCLE
         IF (ANY([(ALL(rot(:, :, ir) == lays(1)%sym%mrot(1:2, 1:2, iop)), ir=1, nrot)])) CYCLE
         nrot = nrot + 1
         rot(:, :, nrot) = lays(1)%sym%mrot(1:2, 1:2, iop)
      END DO

      ALLOCATE (owner(0:kgrid(1) - 1, 0:kgrid(2) - 1), bk(2, kgrid(1)*kgrid(2)), wt(kgrid(1)*kgrid(2)))
      owner = 0
      nk = 0
      DO i2 = 0, kgrid(2) - 1
         DO i1 = 0, kgrid(1) - 1
            IF (owner(i1, i2) /= 0) CYCLE
            nk = nk + 1
            k = [REAL(i1)/kgrid(1), REAL(i2)/kgrid(2)]
            bk(:, nk) = k - NINT(k)
            wt(nk) = 0.0
            DO ir = 1, nrot
               DO is = 1, MERGE(2, 1, l_timerev)
                  q = MATMUL(k, REAL(rot(:, :, ir)))*(3 - 2*is)
                  IF (ANY(ABS(q*kgrid - NINT(q*kgrid)) > eps)) CALL juDFT_error( &
                     "The k-point grid is not compatible with the symmetry", &
                     hint="choose kgrid with equal entries", calledby="gfprep_kpts")
                  j1 = MODULO(NINT(q(1)*kgrid(1)), kgrid(1))
                  j2 = MODULO(NINT(q(2)*kgrid(2)), kgrid(2))
                  IF (owner(j1, j2) == 0) THEN
                     owner(j1, j2) = nk
                     wt(nk) = wt(nk) + 1.0
                  END IF
               END DO
            END DO
         END DO
      END DO

      kpts%kptsName = "default"
      kpts%kptsKind = KPTS_KIND_UNSPECIFIED
      kpts%nkpt = nk
      ALLOCATE (kpts%bk(3, nk), kpts%wtkpt(nk))
      kpts%bk(1:2, :) = bk(:, :nk)
      kpts%bk(3, :) = 0.0
      kpts%wtkpt = wt(:nk)/SUM(wt(:nk))
      kpts%nkpt3 = [kgrid(1), kgrid(2), 1]
      WRITE (oUnit, "(a,i0,a,i0,a,i0,a,i0,a)") "gfprep: ", kgrid(1), "x", kgrid(2), " k-grid, ", &
         nk, " k-points (", nrot, " common in-plane rotations)"
   END SUBROUTINE gfprep_kpts

   LOGICAL FUNCTION priv_has_rot(sym, r)
      TYPE(t_sym), INTENT(IN) :: sym
      INTEGER, INTENT(IN)     :: r(2, 2)
      INTEGER :: iop
      priv_has_rot = .FALSE.
      DO iop = 1, sym%nop
         IF (ALL(sym%mrot(1:2, 1:2, iop) == r)) priv_has_rot = .TRUE.
      END DO
   END FUNCTION priv_has_rot

   SUBROUTINE priv_grow(xy, z, nb, tb)
      REAL, ALLOCATABLE, INTENT(INOUT)    :: xy(:, :), z(:)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: nb(:), tb(:)
      REAL, ALLOCATABLE    :: xy2(:, :), z2(:)
      INTEGER, ALLOCATABLE :: i2(:)
      INTEGER :: n
      n = SIZE(z)
      ALLOCATE (xy2(2, 2*n), z2(2*n))
      xy2(:, :n) = xy; z2(:n) = z
      CALL MOVE_ALLOC(xy2, xy); CALL MOVE_ALLOC(z2, z)
      ALLOCATE (i2(2*n)); i2(:n) = nb; CALL MOVE_ALLOC(i2, nb)
      ALLOCATE (i2(2*n)); i2(:n) = tb; CALL MOVE_ALLOC(i2, tb)
   END SUBROUTINE priv_grow

END MODULE m_gfprep_layer
