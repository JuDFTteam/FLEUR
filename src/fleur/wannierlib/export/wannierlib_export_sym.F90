!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Writes what a Wannier post-processing code needs in order to symmetrise a wannierisation:
!> the k-point stars under the space group, and the group that produced them.
!>
!> WHY A FILE OF OUR OWN AND NOT THE CONSUMER'S. WannierBerri reads its symmetrizer from a
!> .npz whose keys are its own internal layout, and Wannier90 reads a .dmn; both belong to
!> somebody else and both have moved -- WannierBerri dropped .dmn support entirely. What is
!> written here is FLEUR's own HDF, versioned at the root, and a short reader on the
!> consumer's side turns it into whatever that consumer wants. A format change then costs a
!> python file rather than a Fortran one, and FLEUR never has to track a layout it does not
!> own.
!>
!> NOTHING IS CONVERTED HERE. Indices are 1-based because that is how FLEUR holds them, the
!> rotations are sym%mrot as stored -- the fractional-basis matrices acting on real-space
!> coordinates, NOT the reciprocal-space rotation the stars were built with -- and the
!> k-points keep kpts%bkf's own (3,nk) layout. A reader wanting 0-based indices, the inverse
!> transpose or an (nk,3) array does that itself. Converting on both sides is how an
!> off-by-one survives: it then has to be wrong in two places to show up, and it never is.
!>
!> WHAT IS DELIBERATELY ABSENT.
!>   mp_grid     t_kpts documents nkpt3 as "not supposed to be reliable", and a consumer can
!>               recover the mesh from the k-point list, which is written.
!>   D_wann      not written, and it does not have to be: the consumer builds the
!>               representation of the group on the Wannier functions from the projections,
!>               which /proj declares. What only FLEUR knows is the table, not the matrices.
!>
!> WHAT IS WRITTEN, AND THE VERSION THAT SAYS SO. The root carries a layout version: 2 for
!> the file without band representations, 3 for the file with them. A reader decides from
!> that one integer instead of probing for a group.
!>   /kstar       kirr, k2kirr, k2kirr_sym, weight, kirr2k -- the stars, and which operation
!>                carries each k of the mesh to the representative of its own.
!>   /spacegroup  mrot, tau, invtab, antiunitary -- the operations the stars were built
!>                from. NOT all of sym; see the index rule below.
!>   /proj        what every Wannier function was projected from, which is what lets a
!>                reader build D_wann for itself.
!>   /d_band      re and im of (num_bands, num_bands, nsym, nkirr). Version 3 only: a run
!>                whose eigenvectors are real, or whose atom map came out empty, writes
!>                version 2 and the group is simply not there.
!>   /meta        kpoints, amat, bmat, taual.
!>
!> THE OPERATION INDEX IS THIS FILE'S OWN, not sym's. /spacegroup holds only the operations
!> that survived the atom map, in the order they are written, and every index in /kstar and
!> every third index of /d_band counts in THAT order. A reader that uses one of them to
!> index sym reads a different operation, and nothing in the file would contradict it.
!> /antiunitary says which kind each slot is: those marked 1 act on the conjugated state.
!>
!> The lattice IS written, although inp.xml carries it too and a second copy can disagree
!> with the first. The reason is that a symmetry consumer has to use the cell the operations
!> were derived in, and rebuilding it from inp.xml means applying <bulkLattice scale> and
!> trusting that nothing transformed the cell in between -- inpgen does move the origin. One
!> 3x3 written beside the operations removes that whole class of silent failure.
MODULE m_wannierlib_export_sym
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_wannierlib_sym_kstar, ONLY: t_wl_kstar
   USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
   USE m_eulerrot, ONLY: eulerrot
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools
   USE m_wannierlib_hdf_util, ONLY: wr_r1, wr_r2, wr_r3, wr_r4, wr_i1, wr_i2, wr_i3, &
                                    wl_hdf_create, wl_hdf_root
#endif
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_sym
CONTAINS

   !> Write WF<spin>_sym.hdf: the stars, the group they were built from, what each Wannier
   !> function was projected from and, where it could be computed, the group's action on the
   !> bands. Rank 0 writes and the other ranks return at once. The module header above is the
   !> layout, the version rule and the list of what is left out on purpose.
   SUBROUTINE wannierlib_export_sym(star, manifold, wannierize, atoms, cell, kpts, sym, &
                                    fmpi, jspin, l_spinors, dband)
      USE m_types_cell; USE m_types_kpts; USE m_types_sym; USE m_types_mpi
      USE m_types_atoms, ONLY: t_atoms
      USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
      TYPE(t_wl_kstar), INTENT(IN) :: star
      TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wannierize
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: jspin
      LOGICAL, INTENT(IN) :: l_spinors
      !> The group on the bands. Optional and allowed to arrive unallocated, because a run
      !> with local orbitals or with spinor states cannot build it yet and still wants the
      !> rest of this file; has_d_band below says which of the two the reader is holding.
      COMPLEX, ALLOCATABLE, INTENT(IN), OPTIONAL :: dband(:, :, :, :)

#ifdef CPP_HDF
      CHARACTER(LEN=32) :: filename
      INTEGER(HID_T) :: fid, gid
      INTEGER :: err, nw, iw, na
      LOGICAL :: l_dband
      INTEGER, ALLOCATABLE :: iatom(:), mr(:, :, :), iv(:), au(:)
      REAL, ALLOCATABLE :: ta(:, :)
      INTEGER :: iop, jop
      REAL, ALLOCATABLE :: pos(:, :), shift(:, :), euler(:, :), amx(:, :, :)

      IF (fmpi%irank /= 0) RETURN
      !> wl_build_kstar either filled this or stopped the run, so an empty table here means
      !> the export was reached without it -- a wiring mistake, not a bad input.
      IF (star%nkirr == 0) CALL juDFT_error( &
         "wannierlib: the symmetry export was reached with no k-point stars built", &
         calledby="wannierlib_export_sym")

      l_dband = .FALSE.
      IF (PRESENT(dband)) THEN
         IF (ALLOCATED(dband)) l_dband = .TRUE.
      END IF
      IF (l_dband) THEN
         IF (SIZE(dband, 1) /= manifold%num_bands .OR. SIZE(dband, 2) /= manifold%num_bands &
             .OR. SIZE(dband, 3) /= star%nsym .OR. SIZE(dband, 4) /= star%nkirr) &
            CALL juDFT_error("wannierlib: d_band does not have the shape of this k-set", &
                             calledby="wannierlib_export_sym")
      END IF

      CALL wl_hdf_create('sym', jspin, fid, filename)

      !> The stars. kirr and kirr2k index the FULL mesh (1..nk); k2kirr is a star number
      !> (1..nkirr); k2kirr_sym is an operation index into mrot below.
      CALL h5gcreate_f(fid, 'kstar', gid, err)
      CALL wr_i1(gid, 'kirr', (/star%nkirr/), star%kirr(:star%nkirr))
      CALL wr_i1(gid, 'k2kirr', (/star%nk/), star%k2kirr(:star%nk))
      CALL wr_i1(gid, 'k2kirr_sym', (/star%nk/), star%k2kirr_sym(:star%nk))
      CALL wr_i1(gid, 'weight', (/star%nkirr/), star%weight(:star%nkirr))
      CALL wr_i2(gid, 'kirr2k', (/star%nkirr, star%nsym/), star%kirr2k(:star%nkirr, :star%nsym))
      CALL h5gclose_f(gid, err)

      !> The group, as t_sym holds it: mrot in the fractional basis, tau the fractional
      !> translation, invtab the index of each operation's inverse.
      !>
      !> Written out to nsym rather than to nop. With time reversal the second half is the
      !> same rotation and the same translation with antiunitary set, and a reader indexing
      !> mrot(iop) for any operation the star table names should find it there instead of
      !> having to know the halving rule -- a reader that got that rule wrong would apply
      !> the right rotation without the conjugation, which is a symmetrisation that
      !> converges. invtab follows the group: the inverse of T R is T R^-1.
      ALLOCATE (mr(3, 3, star%nsym), ta(3, star%nsym), iv(star%nsym), au(star%nsym))
      DO iop = 1, star%nsym
         jop = star%op(iop)
         mr(:, :, iop) = sym%mrot(:, :, jop)
         ta(:, iop) = sym%tau(:, jop)
         iv(iop) = star%invtab(iop)
         au(iop) = MERGE(1, 0, star%antiunit(iop))
      END DO
      CALL h5gcreate_f(fid, 'spacegroup', gid, err)
      CALL wr_i3(gid, 'mrot', (/3, 3, star%nsym/), mr)
      CALL wr_r2(gid, 'tau', (/3, star%nsym/), ta)
      CALL wr_i1(gid, 'invtab', (/star%nsym/), iv)
      CALL wr_i1(gid, 'antiunitary', (/star%nsym/), au)
      CALL h5gclose_f(gid, err)
      DEALLOCATE (mr, ta, iv, au)

      !> ONE ENTRY PER WANNIER FUNCTION, in the order init left wannierize%proj in, which IS
      !> the column order of the .amn written beside this file. That order is the one thing
      !> here a consumer cannot derive: init expands the declared projections as
      !> spin -> projection -> atom type -> atom -> m, so a species holding several
      !> equivalent atoms interleaves the atoms with the harmonics, while a consumer
      !> rebuilding the same set from a Wyckoff orbit numbers it its own way. The two differ
      !> by a permutation, and a permuted D_wann converges and gives wrong bands, so the
      !> table is written out rather than described.
      !>
      !> l < 0 is a Wannier90 hybrid (-1 sp, -2 sp2, -3 sp3, -4 sp3d, -5 sp3d2) and m is
      !> Wannier90's mr over REAL harmonics, not the magnetic quantum number.
      !>
      !> pos is the centre the .amn actually used: the atom plus the projection's shift,
      !> which enters A_mn as exp(-2 pi i k.shift) and so moves the trial orbital off the
      !> nucleus. shift is written beside it because a reader comparing pos against inp.xml
      !> otherwise has no way to tell a shifted projection from a moved atom.
      !>
      !> amx is built with eulerrot, the same routine wannierlib_amn rotates the angular
      !> coefficients with, so the Euler convention is applied once, here, and never
      !> re-derived from alpha/beta/gamma on the other side. The angles go out too, for a
      !> reader that wants to check this matrix rather than trust it.
      nw = manifold%num_wann
      IF (.NOT. ALLOCATED(wannierize%proj)) CALL juDFT_error( &
         "wannierlib: the symmetry export found no projection table", &
         hint="add <wannierproj .../> children to a <species>", &
         calledby="wannierlib_export_sym")
      IF (SIZE(wannierize%proj) < nw) CALL juDFT_error( &
         "wannierlib: fewer projections than Wannier functions", &
         hint="the projection table is what fixes the column order of the .amn", &
         calledby="wannierlib_export_sym")

      ALLOCATE (iatom(nw), pos(3, nw), shift(3, nw), euler(3, nw), amx(3, 3, nw))
      DO iw = 1, nw
         !> proj%atom counts within its type, which is what the augmentation coefficients
         !> are indexed by; the global atom is what a consumer needs to place it.
         na = atoms%firstAtom(wannierize%proj(iw)%ntype) + wannierize%proj(iw)%atom - 1
         iatom(iw) = na
         shift(:, iw) = wannierize%proj(iw)%shift(:)
         pos(:, iw) = atoms%taual(:, na) + shift(:, iw)
         euler(:, iw) = (/wannierize%proj(iw)%alpha, wannierize%proj(iw)%beta, &
                          wannierize%proj(iw)%gamma/)
      END DO
      CALL eulerrot(nw, wannierize%proj(:nw)%alpha, wannierize%proj(:nw)%beta, &
                    wannierize%proj(:nw)%gamma, amx)

      CALL h5gcreate_f(fid, 'proj', gid, err)
      CALL wr_i1(gid, 'l', (/nw/), wannierize%proj(:nw)%l)
      CALL wr_i1(gid, 'm', (/nw/), wannierize%proj(:nw)%m)
      CALL wr_i1(gid, 'spin', (/nw/), wannierize%proj(:nw)%spin)
      CALL wr_i1(gid, 'species', (/nw/), wannierize%proj(:nw)%ntype)
      CALL wr_i1(gid, 'atom', (/nw/), iatom)
      CALL wr_i1(gid, 'rwf', (/nw/), wannierize%proj(:nw)%rwf)
      CALL wr_r2(gid, 'pos', (/3, nw/), pos)
      CALL wr_r2(gid, 'shift', (/3, nw/), shift)
      CALL wr_r2(gid, 'euler', (/3, nw/), euler)
      CALL wr_r3(gid, 'amx', (/3, 3, nw/), amx)
      CALL wr_r1(gid, 'weight', (/nw/), wannierize%proj(:nw)%weight)
      CALL wr_r1(gid, 'zona', (/nw/), wannierize%proj(:nw)%zona)
      !> j < 0 marks a projection given as (l, m): the j-resolved branch of the trial
      !> orbitals is the one that reads j and mj, and it only runs for spinor states.
      CALL wr_r1(gid, 'j', (/nw/), wannierize%proj(:nw)%j)
      CALL wr_r1(gid, 'mj', (/nw/), wannierize%proj(:nw)%mj)
      CALL h5gclose_f(gid, err)
      DEALLOCATE (iatom, pos, shift, euler, amx)

      !> Real and imaginary parts as two datasets rather than one compound: a compound type
      !> has to be declared on both sides and agreed on, and this is read by a script.
      IF (l_dband) THEN
         CALL h5gcreate_f(fid, 'd_band', gid, err)
         CALL wr_r4(gid, 're', (/manifold%num_bands, manifold%num_bands, star%nsym, &
                                 star%nkirr/), REAL(dband))
         CALL wr_r4(gid, 'im', (/manifold%num_bands, manifold%num_bands, star%nsym, &
                                 star%nkirr/), AIMAG(dband))
         CALL h5gclose_f(gid, err)
      END IF

      CALL h5gcreate_f(fid, 'meta', gid, err)
      CALL wr_r2(gid, 'kpoints', (/3, star%nk/), kpts%bkf(:, :star%nk))
      CALL wr_r2(gid, 'amat', (/3, 3/), cell%amat)
      CALL wr_r2(gid, 'bmat', (/3, 3/), cell%bmat)
      CALL wr_r2(gid, 'taual', (/3, atoms%nat/), atoms%taual(:, :atoms%nat))
      CALL h5gclose_f(gid, err)

      !> version and nkpt come from wl_hdf_root, which is the part of the contract every
      !> WF<spin>_<what>.hdf shares. time_reversal states that the stars were built under the
      !> unitary group alone: a consumer with an antiunitary slot fills it with .FALSE.
      CALL wl_hdf_root(fid, MERGE(3, 2, l_dband), star%nk, gid)
      CALL io_write_attint0(gid, 'nkirr', star%nkirr)
      CALL io_write_attint0(gid, 'nsym', star%nsym)
      CALL io_write_attint0(gid, 'num_bands', manifold%num_bands)
      CALL io_write_attint0(gid, 'num_wann', manifold%num_wann)
      CALL io_write_attint0(gid, 'jspin', jspin)
      CALL io_write_attint0(gid, 'nop', star%nop)
      CALL io_write_attint0(gid, 'time_reversal', MERGE(1, 0, star%nsym > star%nop))
      CALL io_write_attint0(gid, 'has_d_band', MERGE(1, 0, l_dband))
      !> spinor says the states are two-component (noco OR soc), which is what decides
      !> whether a consumer builds the representation over spinors or over two independent
      !> scalar channels. FLEUR knows it; deducing it on the other side from the +1/-1 of
      !> /proj/spin works but is a guess, and a projection whose spin the user fixed in a
      !> collinear run carries the same marks.
      CALL io_write_attint0(gid, 'spinor', MERGE(1, 0, l_spinors))
      CALL h5gclose_f(gid, err)
      CALL h5fclose_f(fid, err)

      WRITE (oUnit, '(a,3(i0,a))') 'wannierlib: wrote '//TRIM(filename)//' (', &
         star%nkirr, ' stars of ', star%nk, ' k-points, ', nw, &
         ' projections'//MERGE(', with d_band)', '; no d_band)  ', l_dband)
#else
      CALL juDFT_error("wannierlib: the symmetry export needs FLEUR built with HDF5", &
                       calledby="wannierlib_export_sym")
#endif
   END SUBROUTINE wannierlib_export_sym

END MODULE m_wannierlib_export_sym
