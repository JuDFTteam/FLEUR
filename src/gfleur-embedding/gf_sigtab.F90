!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_sigtab
   !Scattering amplitudes of a periodic embedding potential on a plane,
   !   sigma_u(p) = <p+u|Sigma|p>,
   !p an in-plane momentum (fractional coordinates of the substrate
   !reciprocal lattice, extended zone), u a substrate reciprocal vector.
   !The amplitudes are taken from an embedding-potential file (t_embfile).
   !On a regular Gamma-centered k-mesh the stored k-points are unfolded with
   !the symmetry operations and time reversal, and sigma_u(p) is interpolated
   !in p (4x4 Lagrange, fallback 2x2). For a k-list that is not a mesh only
   !exact lookups are possible.
   USE hdf5
   USE m_hdf_tools
   USE m_juDFT
   USE m_constants, ONLY: oUnit, pi_const
   USE m_gf_io2dmat, ONLY: gf_io2dmat_kmesh
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_embfile, t_sigtab, gf_kcanon

   INTEGER, PARAMETER :: SIGTAB_GRID = 1, SIGTAB_LIST = 2
   PUBLIC :: SIGTAB_GRID, SIGTAB_LIST
   INTEGER, PARAMETER :: INTERP_EXACT = 0, INTERP_LINEAR = 1, INTERP_CUBIC = 3
   PUBLIC :: INTERP_EXACT, INTERP_LINEAR, INTERP_CUBIC

   !Reader of an embedding-potential file (layout of gf_io2dmat, mode E)
   TYPE t_embfile
      CHARACTER(LEN=200)   :: name = ""
      INTEGER(HID_T)       :: fid = -1, varid = -1
      INTEGER              :: ng = 0, spins = 0, ne = 0, nk = 0
      COMPLEX, ALLOCATABLE :: en(:)
      REAL, ALLOCATABLE    :: kpts(:, :)
      INTEGER, ALLOCATABLE :: gvecs(:, :)
      INTEGER              :: gbox(2) = 0
      INTEGER, ALLOCATABLE :: gmap(:, :)
      INTEGER              :: nsym = 0
      INTEGER, ALLOCATABLE :: symop(:, :, :)
      REAL, ALLOCATABLE    :: symtau(:, :)
      LOGICAL :: l_amat = .FALSE., l_rkmax = .FALSE., l_zplane = .FALSE., l_symtau = .FALSE.
      LOGICAL :: l_timerev = .FALSE., l_noco = .FALSE., l_soc = .FALSE.
      REAL    :: amat2d(2, 2) = 0.0, rkmax = 0.0, zplane = -1.0
   CONTAINS
      PROCEDURE :: open => embfile_open
      PROCEDURE :: read => embfile_read
      PROCEDURE :: close => embfile_close
      PROCEDURE :: find_energy => embfile_find_energy
      PROCEDURE :: find_k => embfile_find_k
      PROCEDURE :: gindex => embfile_gindex
      PROCEDURE :: insphere => embfile_insphere
   END TYPE t_embfile

   !Table of sigma_u(p) for one energy and spin channel
   TYPE t_sigtab
      INTEGER :: mode = 0
      !substrate frame: real-space lattice (columns), reciprocal lattice (rows)
      REAL    :: asub(2, 2) = 0.0, bsub(2, 2) = 0.0
      INTEGER :: nmesh(2) = 0
      INTEGER :: nk = 0
      !table k-points (fractional, canonical in [-1/2,1/2); list mode: sorted by k(1))
      REAL, ALLOCATABLE    :: k(:, :)
      !grid mode: mesh index (0:N1-1,0:N2-1) -> table k-point
      INTEGER, ALLOCATABLE :: kidx(:, :)
      !list mode: table k-point -> file k-point and its integer shift k_file=k+kshift
      INTEGER, ALLOCATABLE :: kfile(:), kshift(:, :)
      !Umklapp list and its index map
      INTEGER :: nu = 0, umax = 0
      INTEGER, ALLOCATABLE :: u(:, :), umap(:, :)
      !1: scalar per spin, 2: 2x2 spinor blocks (noco file, or a collinear file
      !for a noco ad-layer: diagonal blocks)
      INTEGER :: ns = 1
      !the data of the file are collinear (spins 1 or 2)
      LOGICAL :: l_collsrc = .TRUE.
      !s(ns*ns,ig,ik,iu) = sigma_u(k+g), g the file vector ig
      COMPLEX, ALLOCATABLE :: s(:, :, :, :)
      LOGICAL, ALLOCATABLE :: filled(:, :, :)
      !in-plane operations as maps of fractional reciprocal row vectors f -> f wop,
      !the phase exp(-2 pi i u'.tau) and time reversal
      INTEGER :: nop = 1
      INTEGER, ALLOCATABLE :: wop(:, :, :)
      REAL, ALLOCATABLE    :: tau(:, :)
      LOGICAL :: l_timerev = .FALSE.
      INTEGER :: interp = INTERP_CUBIC
      REAL    :: eps_k = 1.0E-5
      !vacuum extrapolation outside the cutoff
      LOGICAL :: l_vacextrap = .FALSE.
      REAL    :: vvac = 0.0
      COMPLEX :: energy = (0.0, 0.0)
      !diagnostics (reset with reset_stat)
      INTEGER :: n_exact = 0, n_cubic = 0, n_linear = 0, n_extrap = 0, n_fail = 0
      !interr: max |cubic-linear| over max |sigma| of the evaluations
      REAL    :: interr = 0.0, valmax = 0.0, unfold_dev = 0.0, unfold_scale = 0.0
      INTEGER :: n_skipped = 0
   CONTAINS
      PROCEDURE :: init => sigtab_init
      PROCEDURE :: fill => sigtab_fill
      PROCEDURE :: eval => sigtab_eval
      PROCEDURE :: uindex => sigtab_uindex
      PROCEDURE :: reset_stat => sigtab_reset_stat
      PROCEDURE :: memory => sigtab_memory
   END TYPE t_sigtab

CONTAINS

   !-------------------------------------------------------------------------
   ! t_embfile
   !-------------------------------------------------------------------------
   SUBROUTINE embfile_open(f, name)
      CLASS(t_embfile), INTENT(INOUT) :: f
      CHARACTER(LEN=*), INTENT(IN)    :: name

      INTEGER(HID_T)   :: did, fspace
      INTEGER(HSIZE_T) :: dims(7), maxdims(7)
      INTEGER          :: hdferr, n, iflag, idims(2)
      LOGICAL          :: l_exist

      f%name = name
      INQUIRE (FILE=TRIM(name)//".hdf", EXIST=l_exist)
      IF (.NOT. l_exist) CALL juDFT_error("Embedding potential file not found: "//TRIM(name)//".hdf", &
                                          calledby="gf_sigtab")
      CALL io_hdfopen(TRIM(name)//".hdf", H5F_ACC_RDONLY_F, f%fid, hdferr)
      CALL io_dopen(f%fid, "2Dmat", f%varid, hdferr)
      CALL h5dget_space_f(f%varid, fspace, hdferr)
      CALL h5sget_simple_extent_dims_f(fspace, dims, maxdims, hdferr)
      CALL h5sclose_f(fspace, hdferr)
      f%ng = INT(dims(1)); f%spins = INT(dims(4)); f%ne = INT(dims(5)); f%nk = INT(dims(6))
      IF (dims(2) /= dims(1)) CALL juDFT_error("2Dmat of "//TRIM(name)//" is not square", calledby="gf_sigtab")

      ALLOCATE (f%en(f%ne))
      CALL io_dopen(f%fid, "energies", did, hdferr)
      CALL io_read(did, (/-1, 1/), (/1, f%ne/), "en", f%en)
      CALL io_dclose(did, hdferr)

      CALL io_dopen(f%fid, "gvecs", did, hdferr)
      CALL io_datadim(did, idims)
      IF (idims(2) /= f%ng) CALL juDFT_error("gvecs and 2Dmat of "//TRIM(name)//" differ in size", &
                                                 hint="T-matrix files cannot be used here", calledby="gf_sigtab")
      ALLOCATE (f%gvecs(2, f%ng))
      CALL io_read(did, (/1, 1/), (/2, f%ng/), "gvecs", f%gvecs)
      CALL io_dclose(did, hdferr)
      f%gbox = MAXVAL(ABS(f%gvecs), DIM=2)
      ALLOCATE (f%gmap(-f%gbox(1):f%gbox(1), -f%gbox(2):f%gbox(2)))
      f%gmap = 0
      DO n = 1, f%ng
         f%gmap(f%gvecs(1, n), f%gvecs(2, n)) = n
      END DO

      ALLOCATE (f%kpts(2, f%nk))
      CALL io_dopen(f%fid, "kpts", did, hdferr)
      CALL io_read(did, (/1, 1/), (/2, f%nk/), "kpts", f%kpts)
      f%nsym = 0
      IF (io_attexists(did, "nsym")) CALL io_read_att(did, "nsym", f%nsym)
      ALLOCATE (f%symop(2, 2, MAX(f%nsym, 1)), f%symtau(2, MAX(f%nsym, 1)))
      f%symop = 0; f%symtau = 0.0
      f%symop(1, 1, 1) = 1; f%symop(2, 2, 1) = 1
      IF (f%nsym > 0) CALL io_read_att(did, "symop", f%symop)
      f%l_symtau = io_attexists(did, "symtau")
      IF (f%l_symtau .AND. f%nsym > 0) CALL io_read_att(did, "symtau", f%symtau)
      f%l_amat = io_attexists(did, "amat2d")
      IF (f%l_amat) CALL io_read_att(did, "amat2d", f%amat2d)
      f%l_rkmax = io_attexists(did, "rkmax")
      IF (f%l_rkmax) CALL io_read_att(did, "rkmax", f%rkmax)
      f%l_zplane = io_attexists(did, "zplane")
      IF (f%l_zplane) CALL io_read_att(did, "zplane", f%zplane)
      IF (io_attexists(did, "l_timerev")) CALL io_read_att(did, "l_timerev", f%l_timerev)
      IF (io_attexists(did, "l_noco")) CALL io_read_att(did, "l_noco", f%l_noco)
      IF (io_attexists(did, "l_soc")) CALL io_read_att(did, "l_soc", f%l_soc)
      CALL io_dclose(did, hdferr)
      iflag = 0
      IF (.NOT. f%l_amat) iflag = 1
      WRITE (oUnit, "(a,a,a,i0,a,i0,a,i0,a,i0,a,i0)") "embedding file ", TRIM(name), ".hdf: ng=", f%ng, &
         " spins=", f%spins, " energies=", f%ne, " k-points=", f%nk, " nsym=", f%nsym
      IF (iflag == 1) WRITE (oUnit, "(a)") "   (no cell metadata amat2d in the file)"
   END SUBROUTINE embfile_open

   SUBROUTINE embfile_read(f, ie, ispin, ik, mat)
      !the matrix of energy ie, spin index ispin (noco: 1..4 = 11,22,12,21), k-point ik
      CLASS(t_embfile), INTENT(IN) :: f
      INTEGER, INTENT(IN)          :: ie, ispin, ik
      COMPLEX, INTENT(OUT)         :: mat(:, :)

      REAL, ALLOCATABLE :: a(:, :, :)

      IF (ispin > f%spins) CALL juDFT_error("spin index not in "//TRIM(f%name), calledby="gf_sigtab")
      ALLOCATE (a(f%ng, f%ng, 2))
      CALL io_read(f%varid, (/1, 1, 1, ispin, ie, ik/), (/f%ng, f%ng, 1, 1, 1, 1/), "A", a(:, :, 1))
      CALL io_read(f%varid, (/1, 1, 2, ispin, ie, ik/), (/f%ng, f%ng, 1, 1, 1, 1/), "A", a(:, :, 2))
      mat = CMPLX(a(:, :, 1), a(:, :, 2))
   END SUBROUTINE embfile_read

   SUBROUTINE embfile_close(f)
      CLASS(t_embfile), INTENT(INOUT) :: f
      INTEGER :: hdferr
      IF (f%fid < 0) RETURN
      CALL io_dclose(f%varid, hdferr)
      CALL io_hdfclose(f%fid, hdferr)
      f%fid = -1
   END SUBROUTINE embfile_close

   INTEGER FUNCTION embfile_find_energy(f, z, eps) RESULT(ie)
      CLASS(t_embfile), INTENT(IN) :: f
      COMPLEX, INTENT(IN)          :: z
      REAL, INTENT(IN)             :: eps
      INTEGER :: n
      ie = 0
      DO n = 1, f%ne
         IF (ABS(REAL(f%en(n)) - REAL(z)) < eps .AND. ABS(AIMAG(f%en(n)) - AIMAG(z)) < eps) ie = n
      END DO
   END FUNCTION embfile_find_energy

   INTEGER FUNCTION embfile_find_k(f, k, eps, nshift) RESULT(ik)
      !file k-point equal to k modulo a reciprocal vector: kpts(:,ik) = k + nshift
      CLASS(t_embfile), INTENT(IN) :: f
      REAL, INTENT(IN)             :: k(2), eps
      INTEGER, INTENT(OUT)         :: nshift(2)
      INTEGER :: n
      REAL    :: d(2)
      ik = 0
      nshift = 0
      DO n = 1, f%nk
         d = f%kpts(:, n) - k
         IF (ALL(ABS(d - NINT(d)) < eps)) THEN
            ik = n
            nshift = NINT(d)
            RETURN
         END IF
      END DO
   END FUNCTION embfile_find_k

   INTEGER FUNCTION embfile_gindex(f, g) RESULT(ig)
      CLASS(t_embfile), INTENT(IN) :: f
      INTEGER, INTENT(IN)          :: g(2)
      ig = 0
      IF (ANY(ABS(g) > f%gbox)) RETURN
      ig = f%gmap(g(1), g(2))
   END FUNCTION embfile_gindex

   SUBROUTINE embfile_insphere(f, ik, mat, bsub, inside)
      !G-vectors of the basis of k-point ik: |k+g|<=rkmax if the file has
      !rkmax, else the nonzero diagonal of mat
      CLASS(t_embfile), INTENT(IN) :: f
      INTEGER, INTENT(IN)          :: ik
      COMPLEX, INTENT(IN)          :: mat(:, :)
      REAL, INTENT(IN)             :: bsub(2, 2)
      LOGICAL, INTENT(OUT)         :: inside(:)
      INTEGER :: n
      REAL    :: q(2)
      DO n = 1, f%ng
         IF (f%l_rkmax) THEN
            q = MATMUL(f%kpts(:, ik) + f%gvecs(:, n), bsub)
            inside(n) = (DOT_PRODUCT(q, q) <= f%rkmax**2*(1.0 + 1.0E-10))
         ELSE
            inside(n) = (ABS(mat(n, n)) > 0.0)
         END IF
      END DO
   END SUBROUTINE embfile_insphere

   !-------------------------------------------------------------------------
   ! t_sigtab
   !-------------------------------------------------------------------------
   SUBROUTINE sigtab_init(tab, f, ulist, interp, asub, l_unfold, l_spinor)
      !k-structure of the table and the Umklapp list (closed under the
      !operations); asub: substrate lattice (columns, substrate frame);
      !l_spinor: 2x2 blocks also from a collinear file
      CLASS(t_sigtab), INTENT(INOUT) :: tab
      TYPE(t_embfile), INTENT(IN)    :: f
      INTEGER, INTENT(IN)            :: ulist(:, :)
      INTEGER, INTENT(IN)            :: interp
      REAL, INTENT(IN)               :: asub(2, 2)
      LOGICAL, INTENT(IN)            :: l_unfold
      LOGICAL, INTENT(IN), OPTIONAL  :: l_spinor

      INTEGER, ALLOCATABLE :: utmp(:, :), cover(:, :)
      INTEGER :: n, iop, nu, i, j, m(2), it, ntr, ik, nfound
      REAL    :: rc(2, 2), w(2, 2), kc(2), kk(2)
      LOGICAL :: l_new, l_mesh

      tab%asub = asub
      tab%bsub = 2.0*pi_const*priv_inv2(asub)
      tab%interp = interp
      tab%l_collsrc = (f%spins <= 2)
      tab%ns = 1
      IF (.NOT. tab%l_collsrc) tab%ns = 2
      IF (PRESENT(l_spinor)) THEN
         IF (l_spinor) tab%ns = 2
      END IF

      !operations as integer maps of fractional reciprocal row vectors:
      !Cartesian R = A mrot A^-1, f' = f B R^T B^-1
      tab%nop = 1
      !without SOC the operations act on collinear spin components unchanged;
      !a noco file needs its full mesh
      IF (l_unfold .AND. tab%l_collsrc) tab%nop = MAX(f%nsym, 1)
      ALLOCATE (tab%wop(2, 2, tab%nop), tab%tau(2, tab%nop))
      DO iop = 1, tab%nop
         rc = MATMUL(MATMUL(asub, REAL(f%symop(:, :, iop))), priv_inv2(asub))
         w = MATMUL(MATMUL(tab%bsub, TRANSPOSE(rc)), priv_inv2(tab%bsub))
         IF (ANY(ABS(w - NINT(w)) > 1.0E-5)) CALL juDFT_error( &
            "Symmetry operation of "//TRIM(f%name)//" does not map the lattice onto itself", &
            hint="check amat2d/symop of the substrate file", calledby="gf_sigtab")
         tab%wop(:, :, iop) = NINT(w)
         tab%tau(:, iop) = f%symtau(:, iop)
      END DO
      tab%l_timerev = l_unfold .AND. f%l_timerev .AND. tab%l_collsrc
      IF (l_unfold .AND. f%nsym > 1 .AND. .NOT. f%l_symtau) &
         CALL juDFT_warn("The substrate file has no symtau: fractional translations are taken as zero", &
                         calledby="gf_sigtab")

      !Umklapp list closed under the operations
      ALLOCATE (utmp(2, SIZE(ulist, 2)*tab%nop + 1))
      nu = 0
      DO n = 1, SIZE(ulist, 2)
         DO iop = 1, tab%nop
            m = MATMUL(ulist(:, n), tab%wop(:, :, iop))
            l_new = .TRUE.
            DO i = 1, nu
               IF (ALL(utmp(:, i) == m)) l_new = .FALSE.
            END DO
            IF (l_new) THEN
               nu = nu + 1
               utmp(:, nu) = m
            END IF
         END DO
      END DO
      tab%nu = nu
      tab%u = utmp(:, :nu)
      tab%umax = MAXVAL(ABS(tab%u))
      ALLOCATE (tab%umap(-tab%umax:tab%umax, -tab%umax:tab%umax))
      tab%umap = 0
      DO n = 1, nu
         tab%umap(tab%u(1, n), tab%u(2, n)) = n
      END DO

      !regular mesh covered by the stored k-points and their images?
      CALL gf_io2dmat_kmesh(f%kpts, tab%nmesh, l_mesh)
      l_mesh = ALL(tab%nmesh > 0)
      ntr = 1
      IF (tab%l_timerev) ntr = 2
      IF (l_mesh) THEN
         ALLOCATE (cover(0:tab%nmesh(1) - 1, 0:tab%nmesh(2) - 1))
         cover = 0
         DO ik = 1, f%nk
            DO iop = 1, tab%nop
               DO it = 1, ntr
                  kk = MATMUL(f%kpts(:, ik), REAL(tab%wop(:, :, iop)))*(3 - 2*it)
                  m = MODULO(NINT(kk*tab%nmesh), tab%nmesh)
                  cover(m(1), m(2)) = 1
               END DO
            END DO
         END DO
         l_mesh = ALL(cover == 1)
         IF (.NOT. l_mesh) WRITE (oUnit, "(a,i0,a,i0,a)") "sigtab: the k-points do not cover the ", &
            tab%nmesh(1), "x", tab%nmesh(2), " mesh, exact lookups only"
      END IF

      IF (l_mesh) THEN
         tab%mode = SIGTAB_GRID
         tab%nk = tab%nmesh(1)*tab%nmesh(2)
         ALLOCATE (tab%k(2, tab%nk), tab%kidx(0:tab%nmesh(1) - 1, 0:tab%nmesh(2) - 1))
         n = 0
         DO j = 0, tab%nmesh(2) - 1
            DO i = 0, tab%nmesh(1) - 1
               n = n + 1
               tab%k(:, n) = (/REAL(MODULO(i + tab%nmesh(1)/2, tab%nmesh(1)) - tab%nmesh(1)/2)/tab%nmesh(1), &
                               REAL(MODULO(j + tab%nmesh(2)/2, tab%nmesh(2)) - tab%nmesh(2)/2)/tab%nmesh(2)/)
               tab%kidx(i, j) = n
            END DO
         END DO
      ELSE
         !list mode: canonical k, sorted by the first component, no unfolding
         tab%mode = SIGTAB_LIST
         tab%nop = 1
         tab%l_timerev = .FALSE.
         tab%nk = f%nk
         ALLOCATE (tab%k(2, f%nk), tab%kfile(f%nk), tab%kshift(2, f%nk))
         DO ik = 1, f%nk
            kc = gf_kcanon(f%kpts(:, ik), tab%eps_k)
            tab%k(:, ik) = kc
            tab%kfile(ik) = ik
            tab%kshift(:, ik) = NINT(f%kpts(:, ik) - kc)
         END DO
         CALL priv_sortk(tab%k, tab%kfile, tab%kshift)
         nfound = 0
         DO ik = 2, f%nk
            IF (ALL(ABS(tab%k(:, ik) - tab%k(:, ik - 1)) < tab%eps_k)) nfound = nfound + 1
         END DO
         IF (nfound > 0) CALL juDFT_warn("The substrate k-list contains equivalent k-points", calledby="gf_sigtab")
      END IF
      ALLOCATE (tab%s(tab%ns**2, f%ng, tab%nk, tab%nu), tab%filled(f%ng, tab%nk, tab%nu))
      WRITE (oUnit, "(a,a,a,i0,a,i0,a,i0,a,f10.2,a)") "sigtab: ", MERGE("mesh", "list", tab%mode == SIGTAB_GRID), &
         " mode, k-points: ", tab%nk, ", Umklapp vectors: ", tab%nu, ", operations: ", tab%nop, &
         ", memory: ", tab%memory(), " MB"
   END SUBROUTINE sigtab_init

   REAL FUNCTION sigtab_memory(tab) RESULT(mb)
      CLASS(t_sigtab), INTENT(IN) :: tab
      mb = 0.0
      IF (ALLOCATED(tab%s)) mb = (16.0*SIZE(tab%s) + 4.0*SIZE(tab%filled))/1024.0**2
   END FUNCTION sigtab_memory

   SUBROUTINE sigtab_fill(tab, f, ie, ispin)
      !fill the table from energy ie, spin ispin (ignored for a noco file)
      CLASS(t_sigtab), INTENT(INOUT) :: tab
      TYPE(t_embfile), INTENT(IN)    :: f
      INTEGER, INTENT(IN)            :: ie, ispin

      COMPLEX, ALLOCATABLE :: mat(:, :, :)
      LOGICAL, ALLOCATABLE :: inside(:)
      INTEGER :: ik, ic, ig1, ig2, iu, iu2, iop, it, ntr, kt, ig, ip(2), g(2), u(2), up(2), n
      INTEGER, PARAMETER :: comp(4) = (/1, 4, 3, 2/)
      REAL    :: p(2), pp(2), x(2)
      COMPLEX :: v(4), ph

      tab%s = CMPLX(0.0, 0.0)
      tab%filled = .FALSE.
      tab%unfold_dev = 0.0
      tab%unfold_scale = 0.0
      tab%n_skipped = 0
      ALLOCATE (mat(f%ng, f%ng, tab%ns**2), inside(f%ng))
      ntr = 1
      IF (tab%l_timerev) ntr = 2
      DO ik = 1, f%nk
         IF (tab%ns == 1) THEN
            CALL f%read(ie, MIN(ispin, f%spins), ik, mat(:, :, 1))
         ELSE IF (tab%l_collsrc) THEN
            !collinear file in 2x2 blocks: diagonal 11 and 22
            CALL f%read(ie, 1, ik, mat(:, :, 1))
            CALL f%read(ie, MIN(2, f%spins), ik, mat(:, :, 4))
            mat(:, :, 2:3) = CMPLX(0.0, 0.0)
         ELSE
            !2x2 blocks in column order 11,21,12,22 from the spins 11,22,12,21
            DO ic = 1, 4
               CALL f%read(ie, comp(ic), ik, mat(:, :, ic))
            END DO
         END IF
         CALL f%insphere(ik, mat(:, :, 1), tab%bsub, inside)
         DO ig2 = 1, f%ng
            IF (.NOT. inside(ig2)) CYCLE
            DO ig1 = 1, f%ng
               IF (.NOT. inside(ig1)) CYCLE
               u = f%gvecs(:, ig1) - f%gvecs(:, ig2)
               iu = tab%uindex(u)
               IF (iu == 0) CYCLE
               p = f%kpts(:, ik) + f%gvecs(:, ig2)
               v(:tab%ns**2) = mat(ig1, ig2, :)
               DO iop = 1, tab%nop
                  pp = MATMUL(p, REAL(tab%wop(:, :, iop)))
                  up = MATMUL(u, tab%wop(:, :, iop))
                  ph = EXP(CMPLX(0.0, -2.0*pi_const*DOT_PRODUCT(REAL(up), tab%tau(:, iop))))
                  DO it = 1, ntr
                     !time reversal: sigma_u(p) = sigma_u(-p-u)
                     IF (it == 2) pp = -pp - up
                     iu2 = tab%uindex(up)
                     IF (iu2 == 0) CYCLE
                     IF (tab%mode == SIGTAB_GRID) THEN
                        x = pp*tab%nmesh
                        ip = NINT(x)
                        IF (ANY(ABS(x - ip) > 1.0E-4)) CALL juDFT_error("k-point not on the mesh", &
                                                                        calledby="gf_sigtab")
                        CALL priv_meshsplit(tab, ip, kt, g)
                     ELSE
                        kt = priv_findk(tab, gf_kcanon(pp, tab%eps_k))
                        IF (kt == 0) CALL juDFT_error("internal error in sigtab list mode", calledby="gf_sigtab")
                        g = NINT(pp - tab%k(:, kt))
                     END IF
                     ig = f%gindex(g)
                     IF (ig == 0) THEN
                        tab%n_skipped = tab%n_skipped + 1
                        CYCLE
                     END IF
                     IF (tab%filled(ig, kt, iu2)) THEN
                        tab%unfold_dev = MAX(tab%unfold_dev, MAXVAL(ABS(tab%s(:, ig, kt, iu2) - ph*v(:tab%ns**2))))
                     END IF
                     tab%s(:, ig, kt, iu2) = ph*v(:tab%ns**2)
                     tab%filled(ig, kt, iu2) = .TRUE.
                     tab%unfold_scale = MAX(tab%unfold_scale, MAXVAL(ABS(v(:tab%ns**2))))
                  END DO
               END DO
            END DO
         END DO
      END DO
      n = COUNT(tab%filled)
      WRITE (oUnit, "(a,i0,a,i0,a,es10.2,a,es10.2,a)") "sigtab: filled ", n, " of ", SIZE(tab%filled), &
         " entries, max deviation of symmetry images ", tab%unfold_dev, " (scale ", tab%unfold_scale, ")"
      IF (tab%unfold_scale > 0.0 .AND. tab%unfold_dev > 1.0E-6*tab%unfold_scale) &
         CALL juDFT_warn("Symmetry images of the substrate embedding potential differ", &
                         hint="check symop/symtau of the substrate file or use a full k-mesh (gfprep kfull=t)", &
                         calledby="gf_sigtab")
   END SUBROUTINE sigtab_fill

   INTEGER FUNCTION sigtab_uindex(tab, u) RESULT(iu)
      CLASS(t_sigtab), INTENT(IN) :: tab
      INTEGER, INTENT(IN)         :: u(2)
      iu = 0
      IF (ANY(ABS(u) > tab%umax)) RETURN
      iu = tab%umap(u(1), u(2))
   END FUNCTION sigtab_uindex

   SUBROUTINE sigtab_reset_stat(tab)
      CLASS(t_sigtab), INTENT(INOUT) :: tab
      tab%n_exact = 0; tab%n_cubic = 0; tab%n_linear = 0; tab%n_extrap = 0; tab%n_fail = 0
      tab%interr = 0.0
      tab%valmax = 0.0
   END SUBROUTINE sigtab_reset_stat

   SUBROUTINE sigtab_eval(tab, f, p, iu, val, ok)
      !sigma_u(p), p fractional (substrate), iu index of u in the table
      CLASS(t_sigtab), INTENT(INOUT) :: tab
      TYPE(t_embfile), INTENT(IN)    :: f
      REAL, INTENT(IN)               :: p(2)
      INTEGER, INTENT(IN)            :: iu
      COMPLEX, INTENT(OUT)           :: val(:)
      LOGICAL, INTENT(OUT)           :: ok

      REAL    :: x(2), wx(4), wy(4), lx(2), ly(2), q(2)
      INTEGER :: i0(2), i, j, ip(2)
      COMPLEX :: c(SIZE(val), 4, 4), lin(SIZE(val)), arg
      LOGICAL :: okc(4, 4)

      val = CMPLX(0.0, 0.0)
      ok = .FALSE.
      IF (tab%mode == SIGTAB_GRID) THEN
         x = p*tab%nmesh
         IF (ALL(ABS(x - NINT(x)) < 1.0E-6)) THEN
            CALL priv_lookup(tab, f, NINT(x), iu, val, ok)
            IF (ok) tab%n_exact = tab%n_exact + 1
         ELSE IF (tab%interp /= INTERP_EXACT) THEN
            !4x4 stencil around x
            i0 = FLOOR(x) - 1
            DO j = 1, 4
               DO i = 1, 4
                  CALL priv_lookup(tab, f, i0 + (/i - 1, j - 1/), iu, c(:, i, j), okc(i, j))
               END DO
            END DO
            CALL priv_lagrange4(x(1) - i0(1), wx)
            CALL priv_lagrange4(x(2) - i0(2), wy)
            lx = (/1.0 - (x(1) - i0(1) - 1), x(1) - i0(1) - 1/)
            ly = (/1.0 - (x(2) - i0(2) - 1), x(2) - i0(2) - 1/)
            IF (ALL(okc(2:3, 2:3))) THEN
               lin = 0.0
               DO j = 1, 2
                  DO i = 1, 2
                     lin = lin + lx(i)*ly(j)*c(:, i + 1, j + 1)
                  END DO
               END DO
            END IF
            IF (tab%interp == INTERP_CUBIC .AND. ALL(okc)) THEN
               DO j = 1, 4
                  DO i = 1, 4
                     val = val + wx(i)*wy(j)*c(:, i, j)
                  END DO
               END DO
               ok = .TRUE.
               tab%n_cubic = tab%n_cubic + 1
               tab%interr = MAX(tab%interr, MAXVAL(ABS(val - lin)))
            ELSE IF (ALL(okc(2:3, 2:3))) THEN
               val = lin
               ok = .TRUE.
               tab%n_linear = tab%n_linear + 1
            END IF
         END IF
      ELSE
         i = priv_findk(tab, gf_kcanon(p, tab%eps_k))
         IF (i > 0) THEN
            ip = NINT(p - tab%k(:, i))
            CALL priv_lookup_list(tab, f, i, ip, iu, val, ok)
            IF (ok) tab%n_exact = tab%n_exact + 1
         END IF
      END IF
      IF (ok) THEN
         tab%valmax = MAX(tab%valmax, MAXVAL(ABS(val)))
         RETURN
      END IF
      IF (tab%l_vacextrap) THEN
         !free-electron vacuum beyond the cutoff: diagonal only
         val = CMPLX(0.0, 0.0)
         IF (ALL(tab%u(:, iu) == 0)) THEN
            q = MATMUL(p, tab%bsub)
            arg = DOT_PRODUCT(q, q) + 2.0*(tab%vvac - tab%energy)
            IF (ABS(AIMAG(tab%energy)) < 1.0E-12) arg = CMPLX(REAL(arg), -1.0E-12)
            val(1) = -0.5*SQRT(arg)
            IF (tab%ns == 2) val(4) = val(1)
         END IF
         ok = .TRUE.
         tab%n_extrap = tab%n_extrap + 1
      ELSE
         tab%n_fail = tab%n_fail + 1
      END IF
   END SUBROUTINE sigtab_eval

   SUBROUTINE priv_lookup(tab, f, ip, iu, val, ok)
      !grid mode: value at the mesh point ip (extended zone, units 1/N)
      TYPE(t_sigtab), INTENT(IN)  :: tab
      TYPE(t_embfile), INTENT(IN) :: f
      INTEGER, INTENT(IN)         :: ip(2), iu
      COMPLEX, INTENT(OUT)        :: val(:)
      LOGICAL, INTENT(OUT)        :: ok
      INTEGER :: g(2), ig, kt
      CALL priv_meshsplit(tab, ip, kt, g)
      ig = f%gindex(g)
      ok = .FALSE.
      val = CMPLX(0.0, 0.0)
      IF (ig == 0) RETURN
      ok = tab%filled(ig, kt, iu)
      IF (ok) val = tab%s(:, ig, kt, iu)
   END SUBROUTINE priv_lookup

   SUBROUTINE priv_meshsplit(tab, ip, kt, g)
      !mesh point ip (units 1/N) = k + g with k in [-1/2,1/2): the G-list of a
      !file covers |k+g|<=rkmax for these k only
      TYPE(t_sigtab), INTENT(IN) :: tab
      INTEGER, INTENT(IN)        :: ip(2)
      INTEGER, INTENT(OUT)       :: kt, g(2)
      INTEGER :: c(2), m(2)
      c = MODULO(ip + tab%nmesh/2, tab%nmesh) - tab%nmesh/2
      g = (ip - c)/tab%nmesh
      m = MODULO(c, tab%nmesh)
      kt = tab%kidx(m(1), m(2))
   END SUBROUTINE priv_meshsplit

   SUBROUTINE priv_lookup_list(tab, f, kt, g, iu, val, ok)
      TYPE(t_sigtab), INTENT(IN)  :: tab
      TYPE(t_embfile), INTENT(IN) :: f
      INTEGER, INTENT(IN)         :: kt, g(2), iu
      COMPLEX, INTENT(OUT)        :: val(:)
      LOGICAL, INTENT(OUT)        :: ok
      INTEGER :: ig
      ig = f%gindex(g)
      ok = .FALSE.
      val = CMPLX(0.0, 0.0)
      IF (ig == 0) RETURN
      ok = tab%filled(ig, kt, iu)
      IF (ok) val = tab%s(:, ig, kt, iu)
   END SUBROUTINE priv_lookup_list

   SUBROUTINE priv_lagrange4(t, w)
      !weights of the cubic Lagrange interpolation on the nodes 0,1,2,3
      REAL, INTENT(IN)  :: t
      REAL, INTENT(OUT) :: w(4)
      w(1) = -(t - 1.0)*(t - 2.0)*(t - 3.0)/6.0
      w(2) = t*(t - 2.0)*(t - 3.0)/2.0
      w(3) = -t*(t - 1.0)*(t - 3.0)/2.0
      w(4) = t*(t - 1.0)*(t - 2.0)/6.0
   END SUBROUTINE priv_lagrange4

   FUNCTION gf_kcanon(k, eps) RESULT(kc)
      !k modulo reciprocal vectors in [-1/2,1/2), values within eps of +1/2 map to -1/2
      REAL, INTENT(IN) :: k(2), eps
      REAL             :: kc(2)
      kc = k - FLOOR(k + 0.5 + eps)
   END FUNCTION gf_kcanon

   INTEGER FUNCTION priv_findk(tab, kc) RESULT(ik)
      !list mode: binary search in k(1), then scan
      TYPE(t_sigtab), INTENT(IN) :: tab
      REAL, INTENT(IN)           :: kc(2)
      INTEGER :: lo, hi, mid, n
      ik = 0
      lo = 1; hi = tab%nk
      DO WHILE (lo < hi)
         mid = (lo + hi)/2
         IF (tab%k(1, mid) < kc(1) - tab%eps_k) THEN
            lo = mid + 1
         ELSE
            hi = mid
         END IF
      END DO
      DO n = lo, tab%nk
         IF (tab%k(1, n) > kc(1) + tab%eps_k) EXIT
         IF (ALL(ABS(tab%k(:, n) - kc) < tab%eps_k)) THEN
            ik = n
            RETURN
         END IF
      END DO
   END FUNCTION priv_findk

   SUBROUTINE priv_sortk(k, kf, ks)
      !insertion sort by k(1)
      REAL, INTENT(INOUT)    :: k(:, :)
      INTEGER, INTENT(INOUT) :: kf(:), ks(:, :)
      INTEGER :: i, j, tf, ts(2)
      REAL    :: t(2)
      DO i = 2, SIZE(k, 2)
         t = k(:, i); tf = kf(i); ts = ks(:, i)
         j = i - 1
         DO WHILE (j >= 1)
            IF (k(1, j) <= t(1)) EXIT
            k(:, j + 1) = k(:, j); kf(j + 1) = kf(j); ks(:, j + 1) = ks(:, j)
            j = j - 1
         END DO
         k(:, j + 1) = t; kf(j + 1) = tf; ks(:, j + 1) = ts
      END DO
   END SUBROUTINE priv_sortk

   FUNCTION priv_inv2(a) RESULT(b)
      REAL, INTENT(IN) :: a(2, 2)
      REAL             :: b(2, 2)
      REAL :: det
      det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
      IF (ABS(det) < 1.0E-12) CALL juDFT_error("singular 2x2 matrix", calledby="gf_sigtab")
      b(1, 1) = a(2, 2)/det; b(2, 2) = a(1, 1)/det
      b(1, 2) = -a(1, 2)/det; b(2, 1) = -a(2, 1)/det
   END FUNCTION priv_inv2

END MODULE m_gf_sigtab
