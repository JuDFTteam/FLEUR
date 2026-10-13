!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfprep_inp
   !Input of gfprep (&GFPREP namelist in gfprep_inp) and the layer geometry:
   !detection of the atomic planes of the bulk supercell, choice of the
   !embedding planes and the d, dt, c parameters of each layer.
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_atoms
   USE m_types_cell
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: t_gfprep_inp, t_layergeo, gfprep_read_inp, gfprep_geometry

   INTEGER, PARAMETER :: maxdim = 1000

   TYPE t_gfprep_inp
      !directory of the bulk calculation (with trailing /)
      CHARACTER(LEN=200) :: bulkdir = "./"
      !prefix of the layer input files: <prefix><i>_inp.xml
      CHARACTER(LEN=20)  :: prefix = "layer"
      !atoms closer than plane_tol in z belong to the same atomic plane
      REAL               :: plane_tol = 0.2
      !auto mode: atomic planes per layer and range of planes
      INTEGER            :: planes_per_layer = 1
      INTEGER            :: first_plane = 0, last_plane = 0
      !explicit embedding planes (cartesian z in a.u.)
      REAL, ALLOCATABLE  :: zplanes(:)
      !dt = d + dt_addon
      REAL               :: dt_addon = 1.5
      !constant shift of the potential, 999 means -V_int(G=0)
      REAL               :: vshift = 0.0
      LOGICAL, ALLOCATABLE :: spinflip(:)
      !2D k-point grid
      INTEGER            :: kgrid(2) = 0
      !Gmax of the layers relative to the bulk; the window transform of the
      !interstitial needs more plane waves than the bulk
      REAL               :: gmax_scale = 1.5
   END TYPE t_gfprep_inp

   TYPE t_layergeo
      !embedding planes, physical height, center
      REAL :: z1, z2, c, zc
      !width of the region containing the MT spheres and the film period
      REAL :: d, dt
   END TYPE t_layergeo

CONTAINS

   SUBROUTINE gfprep_read_inp(inp)
      !Reads the optional file gfprep_inp. All entries have defaults.
      TYPE(t_gfprep_inp), INTENT(OUT) :: inp

      CHARACTER(LEN=200) :: bulkdir
      CHARACTER(LEN=20)  :: prefix
      REAL    :: plane_tol, dt_addon, vshift, gmax_scale
      INTEGER :: planes_per_layer, first_plane, last_plane, kgrid(2)
      REAL    :: zplanes(maxdim)
      LOGICAL :: spinflip(maxdim)
      NAMELIST /GFPREP/ bulkdir, prefix, plane_tol, planes_per_layer, first_plane, &
         last_plane, zplanes, dt_addon, vshift, spinflip, kgrid, gmax_scale

      LOGICAL :: l_exist
      INTEGER :: n, ios

      bulkdir = inp%bulkdir; prefix = inp%prefix
      plane_tol = inp%plane_tol; dt_addon = inp%dt_addon; vshift = inp%vshift
      planes_per_layer = inp%planes_per_layer
      first_plane = inp%first_plane; last_plane = inp%last_plane
      kgrid = inp%kgrid
      gmax_scale = inp%gmax_scale
      zplanes = HUGE(1.0)
      spinflip = .FALSE.

      INQUIRE (FILE="gfprep_inp", EXIST=l_exist)
      IF (l_exist) THEN
         OPEN (99, FILE="gfprep_inp", STATUS="old", ACTION="read")
         READ (99, NML=GFPREP, IOSTAT=ios)
         CLOSE (99)
         IF (ios /= 0) CALL juDFT_error("Error reading &GFPREP in gfprep_inp", calledby="gfprep_read_inp")
      ELSE
         WRITE (*, *) "No gfprep_inp found, using defaults"
      END IF

      inp%bulkdir = ADJUSTL(bulkdir)
      n = LEN_TRIM(inp%bulkdir)
      IF (n == 0) THEN
         inp%bulkdir = "./"
      ELSE IF (inp%bulkdir(n:n) /= "/") THEN
         inp%bulkdir = TRIM(inp%bulkdir)//"/"
      END IF
      inp%prefix = prefix
      inp%plane_tol = plane_tol; inp%dt_addon = dt_addon; inp%vshift = vshift
      inp%planes_per_layer = planes_per_layer
      inp%first_plane = first_plane; inp%last_plane = last_plane
      inp%kgrid = kgrid
      inp%gmax_scale = gmax_scale
      IF (inp%gmax_scale < 1.0) CALL juDFT_error("gmax_scale<1", calledby="gfprep_read_inp")
      IF (inp%dt_addon <= 0.0) CALL juDFT_error("dt_addon must be positive", calledby="gfprep_read_inp")
      n = COUNT(zplanes < HUGE(1.0))
      IF (n > 0) THEN
         IF (n < 2) CALL juDFT_error("zplanes needs at least two planes", calledby="gfprep_read_inp")
         inp%zplanes = zplanes(:n)
         IF (ANY(inp%zplanes(2:) <= inp%zplanes(:n - 1))) then
           CALL juDFT_error("zplanes must be given in ascending order", calledby="gfprep_read_inp")
         end if
      END IF
      inp%spinflip = spinflip
      IF (inp%planes_per_layer < 1) CALL juDFT_error("planes_per_layer<1", calledby="gfprep_read_inp")
   END SUBROUTINE gfprep_read_inp

   SUBROUTINE gfprep_geometry(inp, atoms, cell, geo)
      !Determines the embedding planes and the geometry of all layers and
      !documents the choice in gfprep_planes.txt.
      TYPE(t_gfprep_inp), INTENT(INOUT) :: inp
      TYPE(t_atoms), INTENT(IN)         :: atoms
      TYPE(t_cell), INTENT(IN)          :: cell
      TYPE(t_layergeo), ALLOCATABLE, INTENT(OUT) :: geo(:)

      REAL, ALLOCATABLE :: zplane(:), zemb(:)
      INTEGER, ALLOCATABLE :: nplane(:)
      INTEGER :: np, nl, i, p0, p1
      REAL    :: cz

      cz = cell%amat(3, 3)
      CALL priv_atomic_planes(atoms, cz, inp%plane_tol, zplane, nplane)
      np = SIZE(zplane)

      IF (ALLOCATED(inp%zplanes)) THEN
         zemb = inp%zplanes
      ELSE
         !auto mode: embedding planes between the atomic planes
         IF (inp%first_plane == 0) inp%first_plane = 1
         IF (inp%last_plane == 0) inp%last_plane = inp%first_plane + np - 1
         IF (inp%last_plane < inp%first_plane) CALL juDFT_error("last_plane<first_plane", calledby="gfprep_geometry")
         IF (MOD(inp%last_plane - inp%first_plane + 1, inp%planes_per_layer) /= 0) then
           CALL juDFT_error("Number of planes is not a multiple of planes_per_layer", calledby="gfprep_geometry")
         end if
         nl = (inp%last_plane - inp%first_plane + 1)/inp%planes_per_layer
         ALLOCATE (zemb(nl + 1))
         DO i = 0, nl
            p0 = inp%first_plane + i*inp%planes_per_layer
            zemb(i + 1) = 0.5*(priv_zplane(p0 - 1) + priv_zplane(p0))
         END DO
      END IF

      nl = SIZE(zemb) - 1
      ALLOCATE (geo(nl))
      DO i = 1, nl
         CALL priv_layer_geometry(atoms, cz, zemb(i), zemb(i + 1), inp%dt_addon, geo(i))
      END DO
      IF (ABS(SUM(geo(:)%c)/cz - NINT(SUM(geo(:)%c)/cz)) > 1.0E-6) THEN
         WRITE (*, "(a)") "gfprep: note, the layers do not cover an integer number of bulk periods"
         WRITE (oUnit, "(a)") "gfprep: note, the layers do not cover an integer number of bulk periods"
      END IF

      !document the setup
      OPEN (98, FILE="gfprep_planes.txt", STATUS="replace", ACTION="write")
      WRITE (98, "(a)") "# atomic planes of the bulk supercell (one period)"
      WRITE (98, "(a)") "# plane     z(a.u.)  atoms"
      DO i = 1, np
         WRITE (98, "(i7,f12.6,i7)") i, zplane(i), nplane(i)
      END DO
      WRITE (98, "(a)") "# layers: z1, z2 = embedding planes, c = z2-z1, d, dt"
      WRITE (98, "(a)") "# layer         z1          z2           c           d          dt"
      DO i = 1, nl
         WRITE (98, "(i7,5f12.6)") i, geo(i)%z1, geo(i)%z2, geo(i)%c, geo(i)%d, geo(i)%dt
      END DO
      WRITE (98, "(a)") "# to change the embedding planes set zplanes=... in &GFPREP"
      CLOSE (98)
      WRITE (oUnit, "(a,i0,a)") "gfprep: ", np, " atomic planes per bulk period"
      DO i = 1, nl
         WRITE (oUnit, "(a,i4,a,2f12.6,a,3f12.6)") "layer", i, " planes:", geo(i)%z1, geo(i)%z2, &
            " c,d,dt:", geo(i)%c, geo(i)%d, geo(i)%dt
      END DO

   CONTAINS
      REAL FUNCTION priv_zplane(p)
         !z of atomic plane p, continued periodically beyond one period
         INTEGER, INTENT(IN) :: p
         INTEGER :: ip, iper
         iper = FLOOR(REAL(p - 1)/np)
         ip = p - iper*np
         priv_zplane = zplane(ip) + iper*cz
      END FUNCTION priv_zplane
   END SUBROUTINE gfprep_geometry

   SUBROUTINE priv_atomic_planes(atoms, cz, tol, zplane, nplane)
      !Groups the atoms of one bulk period into planes along z
      TYPE(t_atoms), INTENT(IN)         :: atoms
      REAL, INTENT(IN)                  :: cz, tol
      REAL, ALLOCATABLE, INTENT(OUT)    :: zplane(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nplane(:)

      REAL    :: zs(atoms%nat), zsum(atoms%nat)
      INTEGER :: cnt(atoms%nat)
      INTEGER :: na, i, j, np

      zs = MODULO(atoms%pos(3, :atoms%nat), cz)
      !insertion sort, nat is small
      DO i = 2, atoms%nat
         DO j = i, 2, -1
            IF (zs(j) >= zs(j - 1)) EXIT
            zs(j - 1:j) = [zs(j), zs(j - 1)]
         END DO
      END DO
      np = 1; zsum(1) = zs(1); cnt(1) = 1
      DO na = 2, atoms%nat
         IF (zs(na) - zs(na - 1) < tol) THEN
            zsum(np) = zsum(np) + zs(na); cnt(np) = cnt(np) + 1
         ELSE
            np = np + 1; zsum(np) = zs(na); cnt(np) = 1
         END IF
      END DO
      !merge a plane split by the periodic boundary
      IF (np > 1 .AND. zs(1) + cz - zs(atoms%nat) < tol) THEN
         zsum(1) = zsum(1) + zsum(np) - cnt(np)*cz
         cnt(1) = cnt(1) + cnt(np)
         np = np - 1
      END IF
      !a merged first plane may lie slightly below 0, the order is kept
      ALLOCATE (zplane(np), nplane(np))
      zplane = zsum(:np)/cnt(:np)
      nplane = cnt(:np)
   END SUBROUTINE priv_atomic_planes

   SUBROUTINE priv_layer_geometry(atoms, cz, z1, z2, dt_addon, geo)
      !Geometry of the layer between the planes z1<z2: d is the smallest
      !symmetric slab around the center containing the MT spheres of all
      !atoms between the planes (at least c)
      TYPE(t_atoms), INTENT(IN)     :: atoms
      REAL, INTENT(IN)              :: cz, z1, z2, dt_addon
      TYPE(t_layergeo), INTENT(OUT) :: geo

      INTEGER :: na, i3
      REAL    :: z, ext
      LOGICAL :: l_found
      REAL, PARAMETER :: eps = 1.0E-6

      geo%z1 = z1; geo%z2 = z2
      geo%c = z2 - z1
      geo%zc = 0.5*(z1 + z2)
      IF (geo%c <= 0.0) CALL juDFT_error("Embedding planes not ascending", calledby="gfprep_geometry")
      ext = 0.5*geo%c
      l_found = .FALSE.
      DO na = 1, atoms%nat
         DO i3 = FLOOR(z1/cz) - 2, CEILING(z2/cz) + 2
            z = atoms%pos(3, na) + i3*cz
            IF (ABS(z - z1) < eps .OR. ABS(z - z2) < eps) then
              CALL juDFT_error("An atom lies on an embedding plane", calledby="gfprep_geometry")
            end if
            IF (z > z1 .AND. z < z2) THEN
               l_found = .TRUE.
               ext = MAX(ext, ABS(z - geo%zc) + atoms%rmt(atoms%itype(na)))
            END IF
         END DO
      END DO
      IF (.NOT. l_found) CALL juDFT_error("No atoms between two embedding planes", calledby="gfprep_geometry")
      geo%d = 2.0*ext + 1.0E-4
      geo%dt = geo%d + dt_addon
   END SUBROUTINE priv_layer_geometry

END MODULE m_gfprep_inp
