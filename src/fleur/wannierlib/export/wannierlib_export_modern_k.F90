!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> B(k), C(k) and F(k): the three tensors of the modern theory in the Wannier gauge, one
!> step before the Fourier sum that writes WF<n>_bmn.dat, _cmn.dat and _fmn.dat.
!>
!> Asked for by KS="T", which collects what the run holds in k whether or not a
!> gauge is in it, and written to a file of its own rather than beside the Bloch-basis
!> matrices, because these three DO carry a gauge and nothing in a dataset's shape says so.
!>
!> WHY A GAUGE HAS TO BE IN THEM. The contraction over pairs of neighbours b1, b2 needs a
!> band index that means the same state at every k, and the Bloch one does not: band m at
!> k+b1 and band m at k+b2 are different states. Fixing a gauge is what makes the sum over
!> b well defined, and the gauge that enters is the one at the NEIGHBOURS, not at k --
!> V(k+b1)^dagger M V(k+b2) has no V(k) to factor out -- so there is no gauge-free B(k),
!> C(k) or F(k) of shape (nw, nw, ..) to write instead. What does exist without a gauge is
!> the pair matrices themselves, <u_{k+b1}|u_{k+b2}> and the same with H in the middle --
!> the uIu and uHu of the classic route -- and they are nntot^2 larger: on bcc Fe with 36
!> bands, 512 k-points and nntot = 12 that is 1.66 GB against the 56 MB written here.
!>
!> WHAT IT IS FOR. The Fourier sum is the only step left between this file and the .dat
!> files, so a consumer can take the three tensors onto its own R set, or its own
!> convention, without FLEUR running again. It also makes that last step checkable from
!> outside: B(R) rebuilt outside FLEUR from the Bloch-basis overlaps reproduces
!> WF1_bmn.dat to 6e-08, and this file carries the other two down the same path.
!>
!> AND IT CARRIES THE K-POINTS ITSELF. O(R) = (1/N) sum_k exp(-2pi i k.R) O(k) needs the k
!> the tensors were built on, and a consumer that has to go and find them somewhere else
!> can pick up a list in another order and get a wrong answer with no error: the sum is
!> defined for any list of the right length. They were in WF<n>_bmesh.hdf, which is a
!> separate opt-in flag, so a file could be written that nothing on its own could use.
!> Fractional, in the same order the datasets are indexed in -- that order is the contract,
!> not the values. amat goes with them, in Bohr as FLEUR holds it (NOT the Angstrom of the
!> operators above), for a consumer that wants its R set in Cartesian.
!>
!> UNITS. The ones the .dat files are printed in, so the pair is consistent and the
!> Fourier sum between them is a plain phase sum with nothing converted on the way:
!> B in eV*Angstrom, C in eV*Angstrom^2, F in Angstrom^2. B and C arrive in Hartree and
!> are converted HERE, beside each other, because a layout's unit contract belongs to
!> whoever writes it -- which is what the cart2/cart2e pair does on the .dat side. F
!> carries no energy and is passed through. Getting this wrong is invisible in F and in
!> the shape of C: it came out of the first run a clean factor of 27.2114 small.
!>
!> LAYOUT, version 2 -- version 1 had no kpts and no amat. One rank-5 dataset per tensor,
!> (2, nw, nw, ncomp, nkpt) in Fortran order, real part then imaginary, written one
!> component at a time. HDF reverses the dimension order,
!> so a reader sees (nkpt, ncomp, nw, nw, 2) AND GETS EVERY MATRIX TRANSPOSED; undoing
!> that is not optional and nothing on the reader's side reports it. The nine components
!> of C and F run with the second Cartesian index fastest, c = 3*(a-1) + b, which is the
!> order the cart2 layout writes on disk.
MODULE m_wannierlib_export_modern_k
   USE m_juDFT
   USE m_constants, ONLY: oUnit, hartree_to_ev_const
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_modern_k
CONTAINS

   !> Write WF<jspin>_modern_k.hdf. Collective: each rank owns a k-slice and every tensor
   !> has to be gathered before rank 0 can write it.
   SUBROUTINE wannierlib_export_modern_k(b0_loc, c0_loc, f0_loc, gk_loc, kpts, cell, nw, &
                                         jspin, mpicm, irank, l_want_b, l_want_c, l_want_f)
      USE m_types_kpts, ONLY: t_kpts
      USE m_types_cell, ONLY: t_cell
#ifdef CPP_HDF
      USE hdf5
      USE m_wannierlib_hdf_util
#endif
      COMPLEX, ALLOCATABLE, INTENT(IN) :: b0_loc(:, :, :, :)      !< (nw,nw,3,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: c0_loc(:, :, :, :, :)   !< (nw,nw,3,3,nk_loc)
      COMPLEX, ALLOCATABLE, INTENT(IN) :: f0_loc(:, :, :, :, :)   !< (nw,nw,3,3,nk_loc)
      INTEGER,              INTENT(IN) :: gk_loc(:)               !< global k of each slice entry
      TYPE(t_kpts),         INTENT(IN) :: kpts                    !< the full-zone list the tensors are on
      TYPE(t_cell),         INTENT(IN) :: cell
      INTEGER,              INTENT(IN) :: nw, jspin, mpicm, irank
      !> Which of the three KS="T" named. Not the same question as whether the
      !> slice exists: RS="T" may have built C or F for its own .dat file, and a
      !> tensor nobody asked for in k has no business in this one.
      LOGICAL,              INTENT(IN) :: l_want_b, l_want_c, l_want_f
#ifdef CPP_HDF
      INTEGER(HID_T)    :: fid, gid, did
      CHARACTER(LEN=80) :: filename
      INTEGER           :: a, b, ic, e, nkptf
      LOGICAL           :: l_root, l_b, l_c, l_f

      nkptf = kpts%nkptf
      l_root = (irank == 0)
      l_b = l_want_b .AND. ALLOCATED(b0_loc)
      l_c = l_want_c .AND. ALLOCATED(c0_loc)
      l_f = l_want_f .AND. ALLOCATED(f0_loc)

      !> A name that was asked for and whose slice is missing: the b-shells are the only
      !> way that happens, and B says so where it is built. Saying it again here beats a
      !> file with a dataset silently absent from it.
      IF (.NOT. (l_b .OR. l_c .OR. l_f)) THEN
         IF (l_root) WRITE (oUnit, '(a)') 'wannierlib: KS="T" asked for B(k), '// &
            'C(k) or F(k) and none of them could be built -> nothing written'
         RETURN
      END IF

      IF (l_root) THEN
         CALL wl_hdf_create('modern_k', jspin, fid, filename)
         CALL wl_hdf_root(fid, 2, nkptf, gid)
         CALL wr_i1(gid, 'dims', [2], [nw, nkptf])
         CALL wr_r2(gid, 'kpts', [3, nkptf], kpts%bkf(:, 1:nkptf))
         CALL wr_r2(gid, 'amat', [3, 3], cell%amat)
      END IF

      did = -1
      IF (l_b) THEN
         IF (l_root) CALL wl_ds5_open(gid, 'bmn', [2, nw, nw, 3, nkptf], did)
         DO a = 1, 3
            CALL put_k(did, 'bmn', l_root, gk_loc, nkptf, nw, mpicm, a, &
                       hartree_to_ev_const*b0_loc(:, :, a, :))
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_c) THEN
         IF (l_root) CALL wl_ds5_open(gid, 'cmn', [2, nw, nw, 9, nkptf], did)
         ic = 0
         DO a = 1, 3
            DO b = 1, 3
               ic = ic + 1
               CALL put_k(did, 'cmn', l_root, gk_loc, nkptf, nw, mpicm, ic, &
                          hartree_to_ev_const*c0_loc(:, :, a, b, :))
            END DO
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_f) THEN
         IF (l_root) CALL wl_ds5_open(gid, 'fmn', [2, nw, nw, 9, nkptf], did)
         ic = 0
         DO a = 1, 3
            DO b = 1, 3
               ic = ic + 1
               CALL put_k(did, 'fmn', l_root, gk_loc, nkptf, nw, mpicm, ic, f0_loc(:, :, a, b, :))
            END DO
         END DO
         IF (l_root) CALL wl_ds_close(did)
      END IF

      IF (l_root) THEN
         CALL h5gclose_f(gid, e)
         CALL h5fclose_f(fid, e)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: wrote '//TRIM(filename)//' (', nw, &
            ' Wannier functions, ', nkptf, ' k-points and the list they are on, Wannier gauge)'
      END IF
#else
      CALL juDFT_error("wannierlib: bmn/cmn/fmn with KS='T' need FLEUR built with HDF5", &
                       calledby="wannierlib_export_modern_k")
#endif
   END SUBROUTINE wannierlib_export_modern_k

#ifdef CPP_HDF
   !> One component, gathered from the k-slices and written as a slab. The gather is a sum
   !> onto a zeroed array because each k belongs to exactly one rank, which is the idiom
   !> the rest of postproc reduces with.
   SUBROUTINE put_k(did, name, l_root, gk_loc, nkptf, nw, mpicm, icomp, loc)
      USE hdf5
      USE m_wannierlib_hdf_util, ONLY: wl_ds5_slab
#ifdef CPP_MPI
      USE mpi
#endif
      INTEGER(HID_T),   INTENT(IN) :: did
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL,        INTENT(IN) :: l_root
      INTEGER,        INTENT(IN) :: gk_loc(:), nkptf, nw, mpicm, icomp
      COMPLEX,        INTENT(IN) :: loc(:, :, :)        !< (nw,nw,nk_loc)
      COMPLEX, ALLOCATABLE :: full(:, :, :)
      REAL,    ALLOCATABLE :: buf(:, :, :, :, :)
      INTEGER :: kl, ierr

      ALLOCATE (full(nw, nw, nkptf), source=CMPLX(0.0, 0.0))
      DO kl = 1, SIZE(gk_loc)
         full(:, :, gk_loc(kl)) = loc(:, :, kl)
      END DO
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, full, SIZE(full), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                         mpicm, ierr)
#endif
      IF (l_root) THEN
         ALLOCATE (buf(2, nw, nw, 1, nkptf))
         buf(1, :, :, 1, :) = REAL(full)
         buf(2, :, :, 1, :) = AIMAG(full)
         CALL wl_ds5_slab(did, name, [1, 1, 1, icomp, 1], [2, nw, nw, 1, nkptf], buf)
         DEALLOCATE (buf)
      END IF
      DEALLOCATE (full)
   END SUBROUTINE put_k
#endif

END MODULE m_wannierlib_export_modern_k
