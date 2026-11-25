!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>  The single place that knows how the interpolated bands sit inside WF<spin>_interp.hdf,
!>  the way m_wgauge_io is the single place that knows the layout of the real-space operators.
!>
!>  One group per output domain, carrying what the domain alone decides, and one dataset per
!>  quantity interpolated on it:
!>
!>     /<suffix>/   kdist(np)   energy(num_wann,np)   kpts/...
!>                  velocity(3,num_wann,np)   spin(3,num_wann,np)   ...
!>
!>  The abscissa and the energies were a column of every .dat file this replaces, written
!>  again for each quantity. Here the domain owns them once and a quantity carries only its
!>  own components, so adding a quantity adds a dataset and nothing else.
!>
!>  Electronvolts and angstroem throughout, which is what the real-space operator files
!>  already use and what every reader of them expects. Fleur works in hartree and bohr, so
!>  the conversion happens here, at the one place that writes: a caller hands over what it
!>  computed and does not multiply first. Each dataset still records its unit as an
!>  attribute, which is why the second copy of the energies in electronvolts is gone rather
!>  than merely moved -- there is no longer a header to disagree with.
!>
!>  The labels of the high-symmetry points go into a kpts group laid out like the one in
!>  banddos.hdf, so a reader that already draws a Fleur band structure meets no second
!>  convention here.
!>
!>  Whether a domain traces a line, covers a plane or fills the zone is copied from the list
!>  the user named and not decided here: t_wgauge_domains deliberately does not know it, and
!>  the k-point list already carries it as kptsKind.
!>
!>  Needs an HDF5 build. Without one there is nowhere to write any of this, so every entry
!>  point says so instead of writing less than it was asked for.
MODULE m_wgauge_bands_io
   USE m_juDFT
   USE m_types_kpts
   USE m_constants, ONLY: hartree_to_ev_const, bohr_to_angstrom_const
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools
   USE m_wannierlib_hdf_util, ONLY: wr_r1, wr_r2, wr_i1, wl_str1, wl_grp, wl_grp_close, &
                                    wl_ds3_open, wl_ds3_slab, wl_ds_close, &
                                    wl_hdf_create, wl_hdf_root
#endif
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: WL_BANDS_HID, t_wgauge_bands_file, t_wgauge_bands_ds
   PUBLIC :: wgauge_bands_file_open, wgauge_bands_file_close
   PUBLIC :: wgauge_bands_domain, wgauge_bands_domain_close
   PUBLIC :: wgauge_bands_open, wgauge_bands_comp, wgauge_bands_close

   !> An HDF5 object id in a kind that exists whether or not this was built with HDF5. The
   !> callers hold one between calls, so declaring them INTEGER(HID_T) would put a
   !> declaration inside a CPP_HDF branch and its use outside. That breaks the build without
   !> the feature, which is the build no test here runs.
   INTEGER, PARAMETER :: WL_BANDS_HID = SELECTED_INT_KIND(18)

   !> The layout this writes. A reader checks it before anything else, so a file written by
   !> an older Fleur is refused rather than misread.
   INTEGER, PARAMETER :: WL_INTERP_VERSION = 1

   !> The one file of a spin channel: WF<jspin>_interp.hdf, and its root group, which is
   !> what the domains hang from. Both ids travel together because closing one without the
   !> other leaves the file unreadable, and the name is kept so the log can say what it wrote.
   TYPE t_wgauge_bands_file
      INTEGER(WL_BANDS_HID) :: fid = 0
      INTEGER(WL_BANDS_HID) :: gid = 0
      CHARACTER(LEN=80) :: name = ''
   END TYPE t_wgauge_bands_file

   !> A dataset this module has open, and the name to name in a message when a write on it
   !> fails. The handle alone identifies it to HDF5; the name is for whoever reads the error.
   TYPE t_wgauge_bands_ds
      INTEGER(WL_BANDS_HID) :: id = 0
      CHARACTER(LEN=40) :: name = ''
      INTEGER :: ncomp = 0
      INTEGER :: nw = 0
      INTEGER :: np = 0
   END TYPE t_wgauge_bands_ds

CONTAINS

   !> Make WF<jspin>_interp.hdf and leave it open at its root group.
   !>
   !> Created through the same roster and the same root as the other files this directory
   !> writes, so it carries their layout version, their k-point count and their provenance
   !> stamp without restating any of it. The file is opened here and not by the caller
   !> because that is what keeps every HDF5 type inside this module: a caller holds an id,
   !> never an id of a type that only exists in an HDF5 build.
   SUBROUTINE wgauge_bands_file_open(jspin, nkpt, fh)
      INTEGER, INTENT(IN) :: jspin, nkpt
      TYPE(t_wgauge_bands_file), INTENT(OUT) :: fh
#ifdef CPP_HDF
      INTEGER(HID_T) :: f, g
      CHARACTER(LEN=80) :: fname

      CALL wl_hdf_create('interp', jspin, f, fname)
      CALL wl_hdf_root(f, WL_INTERP_VERSION, nkpt, g)
      fh%fid = INT(f, WL_BANDS_HID)
      fh%gid = INT(g, WL_BANDS_HID)
      fh%name = fname
#else
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_file_open')
#endif
   END SUBROUTINE wgauge_bands_file_open

   SUBROUTINE wgauge_bands_file_close(fh)
      TYPE(t_wgauge_bands_file), INTENT(INOUT) :: fh
#ifdef CPP_HDF
      INTEGER :: e

      CALL wl_grp_close(INT(fh%gid, HID_T))
      CALL h5fclose_f(INT(fh%fid, HID_T), e)
      fh%fid = 0
      fh%gid = 0
#else
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_file_close')
#endif
   END SUBROUTINE wgauge_bands_file_close

   !> Make the group of one output domain and fill in everything that belongs to the domain
   !> and not to any quantity measured on it: the abscissa, the labels of the high-symmetry
   !> points, and the band energies.
   !>
   !> kdist arrives in bohr^-1 and evals in hartree, the units Fleur computes them in, and
   !> both are converted on the way out.
   !>
   !> A domain with no labelled points writes no labels, which is what a plane or a mesh
   !> should do, and nSpecialPoints then reads 0 rather than being absent.
   SUBROUTINE wgauge_bands_domain(fh, name, suffix, kset, kdist, evals, gid)
      TYPE(t_wgauge_bands_file), INTENT(IN) :: fh
      !> The group this domain takes, and the suffix the input asked for. They are not the
      !> same string: a group cannot be called nothing and cannot sensibly begin with an
      !> underscore, so the name is sanitised while the suffix is kept verbatim. Without
      !> both, suffix="_plane" and suffix="plane" would be indistinguishable in the output
      !> and nothing could say which <domain> a group came from.
      CHARACTER(LEN=*), INTENT(IN) :: name, suffix
      TYPE(t_kpts), INTENT(IN) :: kset
      REAL, INTENT(IN) :: kdist(:), evals(:, :)
      INTEGER(WL_BANDS_HID), INTENT(OUT) :: gid
#ifdef CPP_HDF
      INTEGER(HID_T) :: dgid, kgid
      INTEGER :: np, nw, isp, ip

      np = SIZE(kdist)
      nw = SIZE(evals, 1)
      IF (SIZE(evals, 2) /= np) CALL juDFT_error( &
         'the domain has one number of points and its energies another', &
         calledby='wgauge_bands_domain')

      CALL wl_grp(INT(fh%gid, HID_T), TRIM(name), dgid)
      !> Written only when there is one: its absence is how the base domain says so.
      IF (LEN_TRIM(suffix) > 0) CALL io_write_attchar0(dgid, 'suffix', TRIM(suffix))
      CALL io_write_attint0(dgid, 'npoints', np)
      CALL io_write_attint0(dgid, 'numWann', nw)
      CALL io_write_attint0(dgid, 'kptsKind', kset%kptsKind)
      CALL io_write_attchar0(dgid, 'listName', TRIM(kset%kptsName))

      CALL wr_r1(dgid, 'kdist', (/np/), kdist/bohr_to_angstrom_const, units='1/Angstrom')
      CALL wr_r2(dgid, 'energy', (/nw, np/), hartree_to_ev_const*evals, units='eV')

      CALL wl_grp(dgid, 'kpts', kgid)
      CALL io_write_attint0(kgid, 'nkpt', np)
      CALL io_write_attint0(kgid, 'nSpecialPoints', kset%numSpecialPoints)
      IF (kset%numSpecialPoints > 0) THEN
         DO isp = 1, kset%numSpecialPoints
            ip = kset%specialPointIndices(isp)
            IF (ip < 1 .OR. ip > np) CALL juDFT_error( &
               'a special point falls outside the domain it labels', &
               calledby='wgauge_bands_domain')
         END DO
         CALL wl_str1(kgid, 'specialPointLabels', kset%numSpecialPoints, &
                      LEN(kset%specialPointNames(1)), kset%specialPointNames)
         CALL wr_i1(kgid, 'specialPointIndices', (/kset%numSpecialPoints/), &
                    kset%specialPointIndices)
      END IF
      CALL wl_grp_close(kgid)

      gid = INT(dgid, WL_BANDS_HID)
#else
      gid = 0
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_domain')
#endif
   END SUBROUTINE wgauge_bands_domain

   SUBROUTINE wgauge_bands_domain_close(gid)
      INTEGER(WL_BANDS_HID), INTENT(IN) :: gid
#ifdef CPP_HDF
      CALL wl_grp_close(INT(gid, HID_T))
#else
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_domain_close')
#endif
   END SUBROUTINE wgauge_bands_domain_close

   !> One dataset for one interpolated quantity, shaped (ncomp, num_wann, npoints) and opened
   !> at its full extent. The components are formed inside the loop over the domain and never
   !> exist as one array, so they arrive a k-point at a time.
   !>
   !> units says what the numbers are in, and it is recorded rather than implied by the name:
   !> a velocity and a Berry curvature are both three components per band.
   SUBROUTINE wgauge_bands_open(gid, name, ncomp, nw, np, units, ds)
      INTEGER(WL_BANDS_HID), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name, units
      INTEGER, INTENT(IN) :: ncomp, nw, np
      TYPE(t_wgauge_bands_ds), INTENT(OUT) :: ds
#ifdef CPP_HDF
      INTEGER(HID_T) :: d

      IF (ncomp < 1) CALL juDFT_error('a quantity with no components has no dataset', &
                                      calledby='wgauge_bands_open')
      CALL wl_ds3_open(INT(gid, HID_T), name, (/ncomp, nw, np/), d)
      CALL io_write_attchar0(d, 'units', TRIM(units))
      ds%id = INT(d, WL_BANDS_HID)
      ds%name = name
      ds%ncomp = ncomp
      ds%nw = nw
      ds%np = np
#else
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_open')
#endif
   END SUBROUTINE wgauge_bands_open

   !> The components of every band at one point of the domain, comp being (ncomp, num_wann).
   !>
   !> The extents are checked against what the dataset was opened with. A slab of the wrong
   !> shape, or a point outside the domain, is a write HDF5 would take without complaint and
   !> a reader would find as a band of zeros much later.
   SUBROUTINE wgauge_bands_comp(ds, ip, comp)
      TYPE(t_wgauge_bands_ds), INTENT(IN) :: ds
      INTEGER, INTENT(IN) :: ip
      REAL, INTENT(IN) :: comp(:, :)
#ifdef CPP_HDF
      IF (SIZE(comp, 1) /= ds%ncomp .OR. SIZE(comp, 2) /= ds%nw) CALL juDFT_error( &
         'a slab of the wrong shape for dataset '//TRIM(ds%name), &
         calledby='wgauge_bands_comp')
      IF (ip < 1 .OR. ip > ds%np) CALL juDFT_error( &
         'a point outside the domain of dataset '//TRIM(ds%name), &
         calledby='wgauge_bands_comp')
      CALL wl_ds3_slab(INT(ds%id, HID_T), TRIM(ds%name), (/1, 1, ip/), &
                       (/ds%ncomp, ds%nw, 1/), RESHAPE(comp, (/ds%ncomp, ds%nw, 1/)))
#else
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_comp')
#endif
   END SUBROUTINE wgauge_bands_comp

   SUBROUTINE wgauge_bands_close(ds)
      TYPE(t_wgauge_bands_ds), INTENT(INOUT) :: ds
#ifdef CPP_HDF
      CALL wl_ds_close(INT(ds%id, HID_T))
      ds%id = 0
#else
      CALL juDFT_error('interpolation output needs an HDF5 build', &
                       calledby='wgauge_bands_close')
#endif
   END SUBROUTINE wgauge_bands_close

END MODULE m_wgauge_bands_io
