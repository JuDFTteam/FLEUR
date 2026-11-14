!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The operator matrices in the BLOCH basis, O(k), as the coarse pass built them and before
!> any gauge exists.
!>
!> WHY THIS IS NOT <operatorsRs> UNDER ANOTHER FILE NAME. O(R) is the end of a chain:
!> wannierise the manifold, rotate every matrix into that gauge, Fourier the result. A
!> consumer that means to build its own Wannier functions -- or none at all -- would have
!> to pay for a gauge it is going to discard and then undo it. These matrices come out one
!> step after the eigenvectors and commit to nothing: the row index is FLEUR's own band
!> index, the last index is the full-zone k list, and what to do next is the reader's.
!>
!> WHAT TRAVELS WITH THEM. The finite-difference b-shells, because they are the one piece
!> of the Wannier90 setup a reader cannot rebuild from inp.xml and the one that is already
!> fixed here: the shell search runs inside init_w90, so wb and bk exist before the first
!> iteration of anything. Without them no derivative in k can be taken at all.
!>
!> WHAT DOES NOT. No centres -- they are a property of a gauge, and the point of this file
!> is that there is not one yet. No Hamiltonian, which in this basis is the eigenvalue list
!> that out.xml and the .eig already carry. No position operator, which at a single k is
!> not a matrix but a finite difference over the shells above.
!>
!> LAYOUT. One rank-5 dataset per operator, (2, nb, nb, ncomp, nkpt) in Fortran order, real
!> part then imaginary. It is written one component at a time: the whole array and the
!> buffer it is gathered through do not fit side by side for an orbital moment on many
!> atoms. HDF reverses the dimension order, so a reader sees (nkpt, ncomp, nb, nb, 2) AND
!> GETS EVERY MATRIX TRANSPOSED -- the third axis is the column and the fourth the row.
!> Undoing that is not optional and nothing on the reader's side will report it: the
!> transpose of a Hermitian matrix is Hermitian, keeps its norm and its eigenvalues, and
!> quietly answers a different question. Transposed back, and rotated with the gauge in
!> the checkpoint, these reproduce rspauli.1 and anglmomrs.1 at R = 0 on bcc Fe to 7e-9,
!> which is the precision those files are printed with. For the orbital moment ncomp runs
!> over (cartesian, atom, channel), fastest index first, and dims carries nat and
!> n_channels so a reader can undo that too.
MODULE m_wannierlib_export_ks
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_mpi, ONLY: t_mpi
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_ks
CONTAINS

   !> Write WF<jspin>_operators_ks.hdf. Every rank takes part -- the matrices are spread
   !> over the ranks by k and have to be gathered -- but only rank 0 holds the file open.
   SUBROUTINE wannierlib_export_ks(l_spin, l_orbital, l_spinorbit, l_torque, &
                                   s0, l0, soc4, t0, x0, n_channels, nat, &
                                   wb, bk, distk, nkptf, nb, jspin, fmpi)
#ifdef CPP_HDF
      USE hdf5
      USE m_wannierlib_hdf_util
#endif
      LOGICAL,              INTENT(IN) :: l_spin, l_orbital, l_spinorbit, l_torque
      COMPLEX, ALLOCATABLE, INTENT(IN) :: s0(:, :, :, :)          !< (nb,nb,3,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: l0(:, :, :, :, :, :)    !< (nb,nb,3,nat,channel,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: soc4(:, :, :, :)        !< (nb,nb,4,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: t0(:, :, :, :)          !< (nb,nb,3,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: x0(:, :, :)             !< (nb,nb,nk_loc)
      INTEGER,              INTENT(IN) :: n_channels, nat
      REAL,    ALLOCATABLE, INTENT(IN) :: wb(:)                   !< (nntot) shell weights
      REAL,    ALLOCATABLE, INTENT(IN) :: bk(:, :, :)             !< the b vectors themselves
      INTEGER,              INTENT(IN) :: distk(:)                !< (nkptf) owning rank of each k
      INTEGER,              INTENT(IN) :: nkptf, nb, jspin
      TYPE(t_mpi),          INTENT(IN) :: fmpi
#ifdef CPP_HDF
      INTEGER(HID_T)    :: fid, gid, did
      CHARACTER(LEN=80) :: filename
      INTEGER           :: ic, iat, ich, icomp, ncomp, e
      LOGICAL           :: l_root

      l_root = (fmpi%irank == 0)

      !> The band count travels as an argument because the caller knows it as the width of
      !> the manifold, not as the shape of an array that may not be allocated. Where both
      !> are available they have to agree, and a mismatch means the manifold moved between
      !> the operator pass and here, which no later check would catch.
      !> Asking for an operator this run does not build has to stop here. These matrices
      !> are assembled in a spinor basis only, so a collinear run reaches this point with
      !> nothing to write, and a file that comes back quietly missing what was asked for
      !> is worse than no file: the name is in inp.xml, so the reader believes it is there.
      CALL assert_asked('spin',           l_spin,      ALLOCATED(s0))
      CALL assert_asked('orbital',        l_orbital,   ALLOCATED(l0))
      CALL assert_asked('spinOrbit',      l_spinorbit, ALLOCATED(soc4))
      CALL assert_asked('exchangeTorque', l_torque,    ALLOCATED(t0))

      CALL assert_nb('spin',           s0,   nb)
      CALL assert_nb('spinOrbit',      soc4, nb)
      CALL assert_nb('exchangeTorque', t0,   nb)
      IF (ALLOCATED(l0)) THEN
         IF (SIZE(l0, 1) /= nb .OR. SIZE(l0, 2) /= nb) &
            CALL juDFT_error("wannierlib: the orbital moment has a different band count "// &
                             "than the manifold", calledby="wannierlib_export_ks")
         IF (SIZE(l0, 4) /= nat) &
            CALL juDFT_error("wannierlib: the orbital moment has a different atom count "// &
                             "than the structure", calledby="wannierlib_export_ks")
      END IF

      IF (l_root) THEN
         CALL wl_hdf_create('operators_ks', jspin, fid, filename)
         CALL wl_hdf_root(fid, 1, nkptf, gid)
         CALL wr_i1(gid, 'dims', [4], [nb, nkptf, nat, n_channels])
         IF (ALLOCATED(wb)) CALL wr_r1(gid, 'bmesh_wb', [SIZE(wb)], wb)
         IF (ALLOCATED(bk)) CALL wr_r3(gid, 'bmesh_bk', SHAPE(bk), bk)
      END IF

      did = -1
      IF (l_spin .AND. ALLOCATED(s0)) THEN
         ncomp = SIZE(s0, 3)
         IF (l_root) CALL wl_ds5_open(gid, 'spin', [2, nb, nb, ncomp, nkptf], did)
         DO ic = 1, ncomp
            CALL put_component(did, 'spin', ic, l_root, fmpi, distk, nkptf, nb, s0(:, :, ic, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      !> The cross-spin block of a collinear jspins=2 run. It has no component axis of its
      !> own and is not an operator the user can name: it goes out whenever the spin
      !> matrices do and the run has it, because without it those matrices are one quadrant
      !> of a 2x2 and a reader cannot tell which.
      IF (l_spin .AND. ALLOCATED(s0) .AND. ALLOCATED(x0)) THEN
         IF (l_root) CALL wl_ds5_open(gid, 'spin_cross', [2, nb, nb, 1, nkptf], did)
         CALL put_component(did, 'spin_cross', 1, l_root, fmpi, distk, nkptf, nb, x0)
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_orbital .AND. ALLOCATED(l0)) THEN
         ncomp = SIZE(l0, 3)*SIZE(l0, 4)*SIZE(l0, 5)
         IF (l_root) CALL wl_ds5_open(gid, 'orbital', [2, nb, nb, ncomp, nkptf], did)
         icomp = 0
         DO ich = 1, SIZE(l0, 5)
            DO iat = 1, SIZE(l0, 4)
               DO ic = 1, SIZE(l0, 3)
                  icomp = icomp + 1
                  CALL put_component(did, 'orbital', icomp, l_root, fmpi, distk, nkptf, nb, &
                                     l0(:, :, ic, iat, ich, :))
               END DO
            END DO
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_spinorbit .AND. ALLOCATED(soc4)) THEN
         ncomp = SIZE(soc4, 3)
         IF (l_root) CALL wl_ds5_open(gid, 'spinOrbit', [2, nb, nb, ncomp, nkptf], did)
         DO ic = 1, ncomp
            CALL put_component(did, 'spinOrbit', ic, l_root, fmpi, distk, nkptf, nb, soc4(:, :, ic, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_torque .AND. ALLOCATED(t0)) THEN
         ncomp = SIZE(t0, 3)
         IF (l_root) CALL wl_ds5_open(gid, 'exchangeTorque', [2, nb, nb, ncomp, nkptf], did)
         DO ic = 1, ncomp
            CALL put_component(did, 'exchangeTorque', ic, l_root, fmpi, distk, nkptf, nb, t0(:, :, ic, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_root) THEN
         CALL h5gclose_f(gid, e)
         CALL h5fclose_f(fid, e)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: wrote '//TRIM(filename)//' (', nb, &
            ' bands, ', nkptf, ' k-points, no gauge applied)'
      END IF
#else
      CALL juDFT_error("wannierlib: the Bloch-basis operator export needs FLEUR built "// &
                       "with HDF5", calledby="wannierlib_export_ks")
#endif
   END SUBROUTINE wannierlib_export_ks

#ifdef CPP_HDF
   !> Gather one component over the ranks that own the k-points and write it as one slab of
   !> an already open dataset. Every rank has to reach the reduction; only rank 0 writes.
   SUBROUTINE put_component(did, name, icomp, l_root, fmpi, distk, nkptf, nb, loc)
      USE hdf5
      USE m_wannierlib_hdf_util, ONLY: wl_ds5_slab
#ifdef CPP_MPI
      USE mpi
#endif
      INTEGER(HID_T),   INTENT(IN) :: did
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER,          INTENT(IN) :: icomp, nkptf, nb
      LOGICAL,          INTENT(IN) :: l_root
      TYPE(t_mpi),      INTENT(IN) :: fmpi
      INTEGER,          INTENT(IN) :: distk(:)
      COMPLEX,          INTENT(IN) :: loc(:, :, :)     !< (nb,nb,nk_loc)
      COMPLEX, ALLOCATABLE :: full(:, :, :)
      REAL,    ALLOCATABLE :: buf(:, :, :, :, :)
      INTEGER :: ikpt, ik_local, ierr

      ALLOCATE (full(nb, nb, nkptf), source=CMPLX(0.0, 0.0))
      ik_local = 0
      DO ikpt = 1, nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE
         ik_local = ik_local + 1
         IF (ik_local > SIZE(loc, 3)) &
            CALL juDFT_error("wannierlib: more k-points assigned to this rank than the "// &
                             "operator holds", calledby="put_component")
         full(:, :, ikpt) = loc(:, :, ik_local)
      END DO
#ifdef CPP_MPI
      IF (fmpi%isize > 1) CALL MPI_ALLREDUCE(MPI_IN_PLACE, full, SIZE(full), &
                                             MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
#endif
      IF (.NOT. l_root) RETURN
      ALLOCATE (buf(2, nb, nb, 1, nkptf))
      buf(1, :, :, 1, :) = REAL(full)
      buf(2, :, :, 1, :) = AIMAG(full)
      CALL wl_ds5_slab(did, name, [1, 1, 1, icomp, 1], [2, nb, nb, 1, nkptf], buf)
   END SUBROUTINE put_component

   !> An operator named in <operatorsKs> that the run never built. The message says what
   !> is missing rather than how to get it, because the condition is one fact -- the basis
   !> is not spinorial -- reached from several directions.
   SUBROUTINE assert_asked(name, l_asked, l_built)
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL,          INTENT(IN) :: l_asked, l_built
      IF (.NOT. l_asked) RETURN
      IF (l_built)       RETURN
      CALL juDFT_error('wannierlib: <operatorsKs> asked for "'//name//'", which this run '// &
                       'does not build; the Bloch-basis operators need a spinor basis '// &
                       '(l_soc or l_noco)', calledby="wannierlib_export_ks")
   END SUBROUTINE assert_asked

   !> The band count of an operator has to be the manifold's, on both indices.
   SUBROUTINE assert_nb(name, arr, nb)
      CHARACTER(LEN=*),     INTENT(IN) :: name
      COMPLEX, ALLOCATABLE, INTENT(IN) :: arr(:, :, :, :)
      INTEGER,              INTENT(IN) :: nb
      IF (.NOT. ALLOCATED(arr)) RETURN
      IF (SIZE(arr, 1) == nb .AND. SIZE(arr, 2) == nb) RETURN
      CALL juDFT_error("wannierlib: "//name//" has a different band count than the "// &
                       "manifold", calledby="wannierlib_export_ks")
   END SUBROUTINE assert_nb
#endif

END MODULE m_wannierlib_export_ks
