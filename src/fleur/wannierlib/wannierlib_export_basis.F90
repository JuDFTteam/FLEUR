!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Writes the two things that no other FLEUR output carries, for external symmetry and
!> Wannier post-processing (irrep, WannierBerri): the plane-wave index list of the basis at
!> each k-point, and the eigenvector coefficients on it.
!>
!> WHY ONLY THESE TWO. Everything else a consumer needs FLEUR already publishes, so putting
!> it here as well would be a second copy that can disagree with the first:
!>    cell, positions, atomic numbers, rkmax, spin angles   inp.xml
!>    k-point list                                          kpts.xml
!>    eigenvalues, Fermi energy                             out.xml
!> The two exceptions are not available anywhere. The G vectors are rebuilt by lapw%init on
!> every call -- a box scan filtered by |k+G| < rkmax and then SORTED -- and never stored,
!> and the rows of an eigenvector are indexed by that ordering, so a reader that regenerated
!> the list itself would have to reproduce the sort exactly or pair every coefficient with
!> the wrong G. And the eigenvectors themselves are gone: close_eig removes the eig file
!> when the run ends.
!>
!> NOTHING IS INTERPRETED HERE. The records are written as they were read, with all rows --
!> the local-orbital ones included -- and the integers that say where each block starts are
!> written beside them. Deciding which rows are the interstitial, converting units, building
!> a moment direction out of two angles: all of that belongs to the reader. The one thing
!> that cannot be left to the reader is the LAYOUT of what is written, because getting it
!> wrong produces wrong symmetry characters and no error, so the file states it rather than
!> expecting the reader to know it.
MODULE m_wannierlib_export_basis
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_export_basis
CONTAINS

   SUBROUTINE wannierlib_export_basis(this, manifold, atoms, cell, input, kpts, sym, noco, &
                                      nococonv, enpara, vtot, fmpi, eig_id, jspin)
      USE m_types_wannierlib
      USE m_types_melem_manifold, ONLY: t_melem_manifold
      USE m_types_atoms; USE m_types_cell; USE m_types_input; USE m_types_kpts
      USE m_types_sym; USE m_types_noco; USE m_types_nococonv; USE m_types_enpara
      USE m_types_potden; USE m_types_mpi; USE m_types_lapw; USE m_types_mat
      USE m_types_abc
      USE m_types_spinor_layout, ONLY: t_spinor_layout
      USE m_matrix_element_factory, ONLY: matrix_element_states
#ifdef CPP_HDF
      USE hdf5
      USE m_hdf_tools
#endif
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_melem_manifold), INTENT(IN) :: manifold
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_potden), INTENT(IN) :: vtot
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: eig_id, jspin

#ifdef CPP_HDF
      TYPE(t_lapw) :: lapw
      TYPE(t_spinor_layout) :: layout
      TYPE(t_mat), POINTER :: zmat(:)
      TYPE(t_abc), POINTER :: abc(:, :)
      CHARACTER(LEN=32) :: filename, kptname
      LOGICAL :: l_ex
      INTEGER :: ik, ib, ig, ir, nb, nrow, ng
      INTEGER, ALLOCATABLE :: ev_list(:)
      REAL, ALLOCATABLE :: buf(:, :, :, :)
      INTEGER(HID_T) :: fid, gid
      INTEGER :: n(4), err

      IF (fmpi%irank /= 0) RETURN
      !> A k of the full mesh whose eigenvectors live at a symmetry-related parent would have
      !> to be rotated into place, and that rotation is what the consumer is being asked to
      !> determine. Wannierisation wants the full-zone mesh anyway (inpgen -noKsym).
      IF (kpts%nkpt /= kpts%nkptf) CALL juDFT_error( &
         "wannierlib: the basis export needs the full-zone k-mesh", &
         hint="generate it with inpgen -noKsym", calledby="wannierlib_export_basis")
      !> A film has a third region, the vacuum, which no plane-wave coefficient describes.
      IF (input%film) CALL juDFT_error( &
         "wannierlib: the basis export is not defined for a film", &
         calledby="wannierlib_export_basis")

      nb = manifold%num_bands
      ev_list = [(ib, ib=manifold%min_band, manifold%max_band)]
      WRITE (filename, '(a,i0,a)') 'WF', jspin, '_basis.hdf'
      INQUIRE (FILE=TRIM(filename), EXIST=l_ex)
      IF (l_ex) CALL system('rm '//TRIM(filename))
      CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, fid, err, H5P_DEFAULT_F, H5P_DEFAULT_F)

      DO ik = 1, kpts%nkptf
         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ik, cell)
         CALL matrix_element_states(eig_id, ik, input, atoms, sym, cell, noco, nococonv, &
                                    enpara, lapw, vtot, fmpi, zmat, abc, ev_list=ev_list, &
                                    l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                    kpts=kpts)
         !> l_both_spinors has to be stated with the SAME value the states were read with:
         !> after a second variation the two records are halves of one spinor, and without
         !> saying so the layout reports two independent spin channels instead.
         CALL layout%init(input, noco, lapw, atoms, &
                          l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco))
         ng = MAXVAL(lapw%nv(:input%jspins))
         nrow = zmat(1)%matsize1

         WRITE (kptname, '(a,i0)') '/kpt_', ik
         CALL h5gcreate_f(fid, TRIM(kptname), gid, err)
         !> The row layout, so that the reader takes it from the file instead of knowing it.
         CALL io_write_attint1(gid, 'nv', lapw%nv(:input%jspins))
         CALL io_write_attint0(gid, 'nlotot', atoms%nlotot)
         CALL io_write_attint0(gid, 'nrec', layout%nrec)
         CALL io_write_attint0(gid, 'row_dn', layout%row_dn)
         CALL io_write_attint0(gid, 'nrow', nrow)
         CALL io_write_attint0(gid, 'nband', nb)
         n(:3) = (/3, ng, input%jspins/)
         CALL wr_i3(gid, 'gvec', n(:3), lapw%gvec(:, :ng, :input%jspins))

         !> Every row as it was read, local orbitals included: which of them are the
         !> interstitial is the reader's decision, and it has the integers above to make it.
         ALLOCATE (buf(2, nrow, nb, layout%nrec))
         DO ir = 1, layout%nrec
            DO ib = 1, nb
               DO ig = 1, nrow
                  IF (zmat(ir)%l_real) THEN
                     buf(1, ig, ib, ir) = zmat(ir)%data_r(ig, ib); buf(2, ig, ib, ir) = 0.0
                  ELSE
                     buf(1, ig, ib, ir) = REAL(zmat(ir)%data_c(ig, ib))
                     buf(2, ig, ib, ir) = AIMAG(zmat(ir)%data_c(ig, ib))
                  END IF
               END DO
            END DO
         END DO
         n = (/2, nrow, nb, layout%nrec/)
         CALL wr_r4(gid, 'z', n, buf)
         DEALLOCATE (buf)
         CALL h5gclose_f(gid, err)
      END DO
      CALL h5gopen_f(fid, '/', gid, err)
      CALL io_write_attint0(gid, 'version', 3)
      CALL io_write_attint0(gid, 'nkpt', kpts%nkptf)
      CALL h5gclose_f(gid, err)
      CALL h5fclose_f(fid, err)
      WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote '//TRIM(filename)//' (', kpts%nkptf, &
         ' k-points; cell, k-list, eigenvalues and spin angles are in inp.xml/kpts.xml/out.xml)'
#else
      CALL juDFT_error("wannierlib: the basis export needs FLEUR built with HDF5", &
                       calledby="wannierlib_export_basis")
#endif
   END SUBROUTINE wannierlib_export_basis

#ifdef CPP_HDF
   !> Dataspace, dataset, write, close -- the same four steps every time, collected so they
   !> are not repeated at each call site.
   SUBROUTINE wr_r4(gid, name, n, dat)
      USE hdf5
      USE m_hdf_tools
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:)
      REAL, INTENT(IN) :: dat(:, :, :, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(4, INT(n(:4), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_real4(did, (/1, 1, 1, 1/), n(:4), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_r4

   SUBROUTINE wr_i3(gid, name, n, dat)
      USE hdf5
      USE m_hdf_tools
      INTEGER(HID_T), INTENT(IN) :: gid
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: n(:), dat(:, :, :)
      INTEGER(HID_T) :: sid, did
      INTEGER :: e
      CALL h5screate_simple_f(3, INT(n(:3), HSIZE_T), sid, e)
      CALL h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, e)
      CALL h5sclose_f(sid, e)
      CALL io_write_integer3(did, (/1, 1, 1/), n(:3), name, dat)
      CALL h5dclose_f(did, e)
   END SUBROUTINE wr_i3
#endif

END MODULE m_wannierlib_export_basis
