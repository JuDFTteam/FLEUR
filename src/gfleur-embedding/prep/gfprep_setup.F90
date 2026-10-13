!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gfprep_setup
   !Writes the geometry information of a layer to gf_setup.hdf:
   !  layer-i/stepinfo: MT spheres within the layer and MT spheres of
   !                    outer atoms cut by the embedding planes
   !  layer-i/embdesc:  descriptors of the left and right embedding plane
   USE hdf5
   USE m_juDFT
   USE m_constants
   USE m_hdf_tools
   USE m_types_atoms
   USE m_types_cell
   USE m_gfprep_inp, ONLY: t_layergeo
   USE m_gf_embdesc, ONLY: t_embdesc, gf_writedesc_setup
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gfprep_write_setup

CONTAINS

   SUBROUTINE gfprep_write_setup(ilay, geo, atoms, cell)
      INTEGER, INTENT(IN)          :: ilay
      TYPE(t_layergeo), INTENT(IN) :: geo
      TYPE(t_atoms), INTENT(IN)    :: atoms   !of the bulk
      TYPE(t_cell), INTENT(IN)     :: cell    !of the bulk

      REAL, ALLOCATABLE    :: newpos(:, :), newrmt(:), cutpos(:, :), cutrmt(:)
      REAL, ALLOCATABLE    :: apos(:, :, :, :), armt(:, :, :)
      INTEGER, ALLOCATABLE :: az(:, :, :)
      TYPE(t_embdesc)      :: desc_left, desc_right
      INTEGER :: na, n, i3, nnew, ncut, nmax, nin(2), nout(2), nall, iz
      REAL    :: cz, zr, rmt, c, xy(2)
      INTEGER(HID_T) :: fid, gid, sid
      INTEGER :: hdferr
      LOGICAL :: l_exist

      cz = cell%amat(3, 3)
      c = geo%c
      nmax = 5*atoms%nat*(CEILING(geo%d/cz) + 4)
      ALLOCATE (newpos(3, nmax), newrmt(nmax), cutpos(3, nmax), cutrmt(nmax))
      ALLOCATE (apos(3, nmax, 3, 2), armt(nmax, 3, 2), az(nmax, 3, 2))
      apos = 0.0; armt = 0.0; az = 0
      nnew = 0; ncut = 0; nin = 0; nout = 0; nall = 0

      DO na = 1, atoms%nat
         n = atoms%itype(na)
         rmt = atoms%rmt(n)
         iz = NINT(atoms%zatom(n))
         xy = MATMUL(cell%amat(1:2, 1:2), atoms%taual(1:2, na))
         DO i3 = FLOOR((geo%zc - geo%d)/cz) - 1, CEILING((geo%zc + geo%d)/cz) + 1
            zr = atoms%pos(3, na) + i3*cz - geo%zc
            IF (ABS(zr - 0.5*c) < rmt .OR. ABS(zr + 0.5*c) < rmt) THEN
               !the MT sphere is cut by an embedding plane
               IF (ABS(zr) + rmt > 0.5*geo%d) THEN
                  ncut = ncut + 1
                  cutpos(:, ncut) = [atoms%taual(1:2, na), zr/geo%dt]
                  cutrmt(ncut) = rmt
               END IF
               IF (ABS(zr) <= 0.5*c) THEN
                  !atom inside the embedding region
                  IF (ABS(zr + 0.5*c) <= rmt) CALL priv_add(nin(1), 1, 1, zr + 0.5*c)
                  IF (ABS(zr - 0.5*c) <= rmt) CALL priv_add(nin(2), 1, 2, 0.5*c - zr)
               ELSE
                  !atom outside the embedding region
                  IF (zr <= 0.0) CALL priv_add(nout(1), 2, 1, ABS(zr + 0.5*c))
                  IF (zr >= 0.0) CALL priv_add(nout(2), 2, 2, ABS(0.5*c - zr))
               END IF
            END IF
            IF (ABS(zr) + rmt <= 0.5*geo%d) THEN
               !atom in the slab
               nnew = nnew + 1
               newpos(:, nnew) = [atoms%taual(1:2, na), zr/geo%dt]
               newrmt(nnew) = rmt
               nall = nall + 1
               apos(1:2, nall, 3, :) = SPREAD(xy, 2, 2)
               apos(3, nall, 3, 2) = zr + 0.5*c
               apos(3, nall, 3, 1) = 0.5*c - zr
               armt(nall, 3, :) = rmt
               az(nall, 3, :) = iz
            END IF
         END DO
      END DO
      WRITE (oUnit, "(a,i0,a,i0,a,i0,a)") "gfprep: layer ", ilay, ": ", nnew, &
         " MT spheres in the slab, ", ncut, " cut outer spheres"

      !stepinfo
      INQUIRE (FILE="gf_setup.hdf", EXIST=l_exist)
      IF (l_exist) THEN
         CALL h5fopen_f("gf_setup.hdf", H5F_ACC_RDWR_F, fid, hdferr)
      ELSE
         CALL h5fcreate_f("gf_setup.hdf", H5F_ACC_TRUNC_F, fid, hdferr)
      END IF
      CALL io_check("Could not open gf_setup.hdf", hdferr)
      CALL h5gcreate_f(fid, io_layername(ilay), gid, hdferr)
      CALL h5gcreate_f(gid, "stepinfo", sid, hdferr)
      CALL io_check("Could not create stepinfo in gf_setup.hdf", hdferr)
      CALL io_write_att(sid, "c", c)
      CALL io_write_att(sid, "ncutat", ncut)
      CALL io_write_att(sid, "nnewat", nnew)
      CALL io_write_att(sid, "newpos", newpos(:, :nnew))
      CALL io_write_att(sid, "newrmt", newrmt(:nnew))
      IF (ncut > 0) THEN
         CALL io_write_att(sid, "cutpos", cutpos(:, :ncut))
         CALL io_write_att(sid, "cutrmt", cutrmt(:ncut))
      END IF
      CALL h5gclose_f(sid, hdferr)
      CALL h5gclose_f(gid, hdferr)
      CALL h5fclose_f(fid, hdferr)

      !embedding-plane descriptors
      CALL priv_desc(desc_left, 1, nin(1), nout(1))
      CALL priv_desc(desc_right, 2, nin(2), nout(2))
      CALL gf_writedesc_setup(ilay, desc_right, desc_left)

   CONTAINS
      SUBROUTINE priv_add(cnt, kind, side, zdist)
         INTEGER, INTENT(INOUT) :: cnt
         INTEGER, INTENT(IN)    :: kind, side
         REAL, INTENT(IN)       :: zdist
         cnt = cnt + 1
         apos(1:2, cnt, kind, side) = xy
         apos(3, cnt, kind, side) = zdist
         armt(cnt, kind, side) = rmt
         az(cnt, kind, side) = iz
      END SUBROUTINE priv_add

      SUBROUTINE priv_desc(desc, side, ni, no)
         TYPE(t_embdesc), INTENT(OUT) :: desc
         INTEGER, INTENT(IN)          :: side, ni, no
         INTEGER :: m
         desc%valid = .TRUE.
         desc%cut_atoms_in = ni
         desc%cut_atoms_out = no
         desc%all_atoms = nall
         desc%ucell = cell%amat(1:2, 1:2)
         desc%dist_aux = 0.5*(geo%d - c)
         desc%dvec = [0.0, 0.0, c]
         m = MAX(ni, no, nall)
         ALLOCATE (desc%atoms_pos(3, m, 3), desc%atoms_rmt(m, 3), desc%atoms_z(m, 3))
         desc%atoms_pos = apos(:, :m, :, side)
         desc%atoms_rmt = armt(:m, :, side)
         desc%atoms_z = az(:m, :, side)
      END SUBROUTINE priv_desc
   END SUBROUTINE gfprep_write_setup

END MODULE m_gfprep_setup
