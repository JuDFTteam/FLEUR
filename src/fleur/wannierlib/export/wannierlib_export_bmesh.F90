!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The finite-difference b-mesh: which k-point neighbours which, and with what weight.
!>
!> WHY THIS IS ITS OWN FILE. The b-mesh is a property of the WANNIERISATION, not of any one
!> consumer: it is settled inside init_w90, before a gauge exists, before the operators are
!> built and before anything is exported, and it is the same for every quantity any of those
!> produce. It travelled inside the operator export first, which put a general fact behind
!> one particular request -- a run that wants the topology and no operators had to ask for
!> operators. It is published once, here, and the other exports do not repeat it.
!>
!> WHAT IT IS FOR. Any derivative in k is a finite difference over these shells:
!>
!>    df/dk_a = sum_b w_b b_a [ f(k+b) - f(k) ]
!>
!> and that sum cannot be taken without knowing WHICH stored k-point is k+b. On a regular
!> mesh a reader could search for it, but at the zone boundary k+b leaves the zone and comes
!> back as another point plus a reciprocal lattice vector; gkpb is that vector, and without
!> it the match is ambiguous exactly where it matters.
!>
!> The weights satisfy sum_b w_b b_a b_c = delta_ac, which is the statement that makes the
!> difference above a derivative. It is true in the units written here and false in any
!> other, so it doubles as a check that the file was read the way it was written.
!>
!> LAYOUT. HDF reverses the dimension order, so a reader sees every array with its indices
!> back to front: nnlist is written (nk, nntot) and arrives (nntot, nk), gkpb is written
!> (3, nk, nntot) and arrives (nntot, nk, 3), bk is written (3, nntot, nk) and arrives
!> (nk, nntot, 3). Index origin is Fortran's: nnlist holds 1-based k-point numbers.
!>
!> bk carries a k index because it needs one. Wannier90 chooses the shells per k-point,
!> so neighbour slot b does NOT name the same vector at every k, and on a mesh coarse
!> enough that +b and -b reach the same stored point it demonstrably does not. A reader
!> that collapses bk over k gets a shell list that is right at the first k and wrong
!> after it, with every norm and every weight still in place.
MODULE m_wannierlib_export_bmesh
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_mpi, ONLY: t_mpi
   USE m_types_kpts, ONLY: t_kpts
   USE m_types_wgauge_bmesh, ONLY: t_wgauge_bmesh
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_bmesh
CONTAINS

   !> Write WF<jspin>_bmesh.hdf. Rank 0 only: nothing here is distributed, every rank holds
   !> the same mesh, and the weights come from Wannier90 which every rank called.
   SUBROUTINE wannierlib_export_bmesh(bmesh, kpts, jspin, fmpi)
#ifdef CPP_HDF
      USE hdf5
      USE m_wannierlib_hdf_util
#endif
      TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
      TYPE(t_kpts),         INTENT(IN) :: kpts
      INTEGER,              INTENT(IN) :: jspin
      TYPE(t_mpi),          INTENT(IN) :: fmpi
#ifdef CPP_HDF
      INTEGER(HID_T)    :: fid, gid
      CHARACTER(LEN=80) :: filename
      INTEGER           :: e, nk

      IF (fmpi%irank /= 0) RETURN
      nk = kpts%nkptf

      !> An empty mesh means Wannier90 never ran its shell search, so there is nothing to
      !> write and saying so beats leaving a file that looks complete.
      IF (bmesh%nntot < 1) CALL juDFT_error( &
         "wannierlib: <export bmesh='T'/> was asked for but no b-shells exist; "// &
         "Wannier90 has not chosen them", calledby="wannierlib_export_bmesh")
      CALL assert_shape('nnlist', ALLOCATED(bmesh%nnlist))
      CALL assert_shape('gkpb',   ALLOCATED(bmesh%gkpb))
      CALL assert_shape('wb',     ALLOCATED(bmesh%wb))
      CALL assert_shape('bk',     ALLOCATED(bmesh%bk))

      CALL wl_hdf_create('bmesh', jspin, fid, filename)
      CALL wl_hdf_root(fid, 1, nk, gid)

      CALL wr_i1(gid, 'nntot', [1], [bmesh%nntot])
      !> Which stored k-point is the neighbour, 1-based, and the reciprocal lattice vector
      !> that brings it back into the zone. The pair is what makes k+b unambiguous.
      CALL wr_i2(gid, 'nnlist', [nk, bmesh%nntot], bmesh%nnlist(1:nk, 1:bmesh%nntot))
      CALL wr_i3(gid, 'gkpb', [3, nk, bmesh%nntot], bmesh%gkpb(:, 1:nk, 1:bmesh%nntot))
      CALL wr_r1(gid, 'wb', [bmesh%nntot], bmesh%wb(1:bmesh%nntot))
      CALL wr_r3(gid, 'bk', [3, bmesh%nntot, nk], bmesh%bk(:, 1:bmesh%nntot, 1:nk))
      !> The k list itself, so the file stands alone: nnlist indexes into it.
      CALL wr_r2(gid, 'kpts', [3, nk], kpts%bkf(:, 1:nk))

      CALL h5gclose_f(gid, e)
      CALL h5fclose_f(fid, e)
      WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: wrote '//TRIM(filename)//' (', nk, &
         ' k-points, ', bmesh%nntot, ' neighbours each)'
#else
      CALL juDFT_error("wannierlib: the b-mesh export needs FLEUR built with HDF5", &
                       calledby="wannierlib_export_bmesh")
#endif
   END SUBROUTINE wannierlib_export_bmesh

#ifdef CPP_HDF
   !> A piece of the mesh that should be here and is not. Each one is filled in a different
   !> place -- the topology before the wannierisation, the weights after Wannier90 chose
   !> them -- so naming the missing one names the step that did not happen.
   SUBROUTINE assert_shape(name, l_there)
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL,          INTENT(IN) :: l_there
      IF (l_there) RETURN
      CALL juDFT_error("wannierlib: the b-mesh has no "//name//"; it is incomplete at the "// &
                       "point it is written", calledby="wannierlib_export_bmesh")
   END SUBROUTINE assert_shape
#endif

END MODULE m_wannierlib_export_bmesh
