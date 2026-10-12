!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_dfpt_dynmat_io
#ifdef CPP_HDF
   use hdf5
   use m_hdf_tools
#endif
   use m_juDFT
   use m_constants
   use m_types_atoms
   use m_types_cell
   use m_types_kpts

   implicit none
   private
   public :: init_dynmat_file, write_dynmat_q, read_dynmat_q, write_dynmat_gamma, read_dynmat_gamma

   integer, parameter :: DYNMAT_FILE_VERSION = 1

#ifdef CPP_HDF
   interface split_cmplx
      module procedure split_cmplx1, split_cmplx2
   end interface split_cmplx
#endif

contains

   subroutine init_dynmat_file(filename, atoms, cell, qpts, l_truncate)
      !! Creates a dynmat HDF file with its header; an existing file is kept (unless l_truncate) and checked.

      character(len=*), intent(in) :: filename
      type(t_atoms),    intent(in) :: atoms
      type(t_cell),     intent(in) :: cell
      type(t_kpts),     intent(in) :: qpts
      logical,          intent(in) :: l_truncate

#ifdef CPP_HDF
      logical        :: l_exist
      integer        :: hdferr
      integer(HID_T) :: fid

      inquire(file=trim(filename), exist=l_exist)
      if (l_exist .and. .not. l_truncate) then
         call open_checked(filename, 3 * atoms%nat, H5F_ACC_RDONLY_F, fid)
      else
         call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, fid, hdferr)
         if (hdferr /= 0) call judft_error("Could not create "//trim(filename), calledby="init_dynmat_file")
         call write_header(fid, atoms, cell, qpts)
      end if
      call io_hdfclose(fid)
#else
      call judft_error("DFPT dynmat IO requires FLEUR compiled with HDF5.", calledby="init_dynmat_file")
#endif
   end subroutine init_dynmat_file

   subroutine write_dynmat_q(filename, iQ, qvec, weight, eigVals, freqs, dynMat, eigVecs, label)
      !! Writes all data of one q-point into /qpoints/q_<iQ>, replacing an existing group.
      !! freqs are in 1/cm; the THz values are derived from them.

      character(len=*),           intent(in) :: filename
      integer,                    intent(in) :: iQ
      real,                       intent(in) :: qvec(3)
      real,                       intent(in) :: weight
      real,                       intent(in) :: eigVals(:)
      complex,                    intent(in) :: freqs(:)
      complex,          optional, intent(in) :: dynMat(:,:)
      complex,          optional, intent(in) :: eigVecs(:,:)
      character(len=*), optional, intent(in) :: label

#ifdef CPP_HDF
      integer(HID_T)                :: fid, gid
      character(len=:), allocatable :: gname

      call open_checked(filename, size(eigVals), H5F_ACC_RDWR_F, fid)
      gname = '/qpoints/q_'//int2str(iQ)
      if (io_groupexists(fid, gname)) call io_gdelete(fid, gname)

      call io_gcreate(fid, gname, gid)
      call io_write_attint0(gid, 'iQ', iQ)
      call io_write_attreal1(gid, 'qvec', qvec)
      call io_write_attreal0(gid, 'weight', weight)
      if (present(label)) then
         if (len_trim(label) > 0) call io_write_attchar0(gid, 'label', trim(label))
      end if
      call io_write_var(gid, 'eigenvalues', eigVals)
      call io_write_var(gid, 'frequencies_cm', split_cmplx(freqs))
      call io_write_var(gid, 'frequencies_THz', split_cmplx(freqs / 33.0))
      if (present(dynMat)) call io_write_var(gid, 'dynmat', split_cmplx(dynMat))
      if (present(eigVecs)) call io_write_var(gid, 'eigenvectors', split_cmplx(eigVecs))
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      call judft_error("DFPT dynmat IO requires FLEUR compiled with HDF5.", calledby="write_dynmat_q")
#endif
   end subroutine write_dynmat_q

   subroutine read_dynmat_q(filename, iQ, dynMat, qvec)
      !! Reads the dynamical matrix of q-point iQ from a dynmat HDF file.

      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: iQ
      complex,          intent(out) :: dynMat(:,:)
      real, optional,   intent(out) :: qvec(3)

#ifdef CPP_HDF
      integer(HID_T)                :: fid, gid
      character(len=:), allocatable :: gname
      real, allocatable             :: buf(:,:,:)

      call open_checked(filename, size(dynMat, 1), H5F_ACC_RDONLY_F, fid)
      gname = '/qpoints/q_'//int2str(iQ)
      if (.not. io_groupexists(fid, gname)) call judft_error("q-point "//int2str(iQ)//" missing in "//trim(filename), calledby="read_dynmat_q")
      call io_gopen(fid, gname, gid)
      if (.not. io_dataexists(gid, 'dynmat')) call judft_error("No dynmat stored for q-point "//int2str(iQ)//" in "//trim(filename), calledby="read_dynmat_q")
      allocate(buf(2, size(dynMat, 1), size(dynMat, 2)))
      call io_read_var(gid, 'dynmat', buf)
      dynMat = cmplx(buf(1,:,:), buf(2,:,:))
      if (present(qvec)) call io_read_attreal1(gid, 'qvec', qvec)
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      dynMat = cmplx(0.0, 0.0)
      if (present(qvec)) qvec = 0.0
      call judft_error("DFPT dynmat IO requires FLEUR compiled with HDF5.", calledby="read_dynmat_q")
#endif
   end subroutine read_dynmat_q

   subroutine write_dynmat_gamma(filename, dynMat0)
      !! Stores the uncorrected Gamma dynamical matrix used for the acoustic sum rule.

      character(len=*), intent(in) :: filename
      complex,          intent(in) :: dynMat0(:,:)

#ifdef CPP_HDF
      integer(HID_T) :: fid, gid

      call open_checked(filename, size(dynMat0, 1), H5F_ACC_RDWR_F, fid)
      if (io_groupexists(fid, '/sumrule')) call io_gdelete(fid, '/sumrule')
      call io_gcreate(fid, '/sumrule', gid)
      call io_write_var(gid, 'dynmat_gamma', split_cmplx(dynMat0))
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      call judft_error("DFPT dynmat IO requires FLEUR compiled with HDF5.", calledby="write_dynmat_gamma")
#endif
   end subroutine write_dynmat_gamma

   subroutine read_dynmat_gamma(filename, dynMat0)
      !! Reads the uncorrected Gamma dynamical matrix for the acoustic sum rule.

      character(len=*), intent(in)  :: filename
      complex,          intent(out) :: dynMat0(:,:)

#ifdef CPP_HDF
      integer(HID_T)    :: fid, gid
      real, allocatable :: buf(:,:,:)

      call open_checked(filename, size(dynMat0, 1), H5F_ACC_RDONLY_F, fid)
      if (.not. io_groupexists(fid, '/sumrule')) then
         call judft_error("Sum rule needs the Gamma dynamical matrix, missing in "//trim(filename)//". Calculate Gamma first.", calledby="read_dynmat_gamma")
      end if
      call io_gopen(fid, '/sumrule', gid)
      allocate(buf(2, size(dynMat0, 1), size(dynMat0, 2)))
      call io_read_var(gid, 'dynmat_gamma', buf)
      dynMat0 = cmplx(buf(1,:,:), buf(2,:,:))
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      dynMat0 = cmplx(0.0, 0.0)
      call judft_error("DFPT dynmat IO requires FLEUR compiled with HDF5.", calledby="read_dynmat_gamma")
#endif
   end subroutine read_dynmat_gamma

#ifdef CPP_HDF
   subroutine open_checked(filename, ndim, access, fid)
      !! Opens an existing dynmat HDF file and checks that it matches the dimension ndim = 3*nat.

      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: ndim
      integer,          intent(in)  :: access
      integer(HID_T),   intent(out) :: fid

      logical        :: l_exist
      integer        :: hdferr, nat_file
      integer(HID_T) :: gid

      inquire(file=trim(filename), exist=l_exist)
      if (.not. l_exist) call judft_error(trim(filename)//" not found.", calledby="open_checked")
      call io_hdfopen(trim(filename), access, fid, hdferr)
      if (hdferr /= 0) call judft_error("Could not open "//trim(filename), calledby="open_checked")
      call io_gopen(fid, '/general', gid)
      call io_read_attint0(gid, 'nAtoms', nat_file)
      call io_gclose(gid)
      if (3 * nat_file /= ndim) call judft_error(trim(filename)//" belongs to a cell with a different number of atoms.", calledby="open_checked")
   end subroutine open_checked

   subroutine write_header(fid, atoms, cell, qpts)
      !! Writes version, cell, atom and q-set information of a new dynmat file.

      integer(HID_T), intent(in) :: fid
      type(t_atoms),  intent(in) :: atoms
      type(t_cell),   intent(in) :: cell
      type(t_kpts),   intent(in) :: qpts

      integer(HID_T) :: gid
      integer        :: iatom
      integer        :: atomicNumbers(atoms%nat)
      real           :: masses(atoms%nat)

      call io_gcreate(fid, '/meta', gid)
      call io_write_attint0(gid, 'version', DYNMAT_FILE_VERSION)
      call io_gclose(gid)

      do iatom = 1, atoms%nat
         atomicNumbers(iatom) = atoms%nz(atoms%itype(iatom))
         masses(iatom) = atomicMasses_const(atomicNumbers(iatom)) * massInElectronMasses
      end do

      call io_gcreate(fid, '/general', gid)
      call io_write_attint0(gid, 'nAtoms', atoms%nat)
      call io_write_attint0(gid, 'nTypes', atoms%ntype)
      call io_write_var(gid, 'amat', cell%amat)
      call io_write_var(gid, 'bmat', cell%bmat)
      call io_write_var(gid, 'positions', atoms%taual(:, :atoms%nat))
      call io_write_var(gid, 'atomicNumbers', atomicNumbers)
      call io_write_var(gid, 'masses', masses)
      call io_gclose(gid)

      call io_gcreate(fid, '/qpoints', gid)
      call io_write_attint0(gid, 'nq', qpts%nkpt)
      call io_write_attchar0(gid, 'kind', trim(adjustl(kptsKindString_consts(qpts%kptsKind))))
      if (qpts%numSpecialPoints > 0 .and. allocated(qpts%specialPointIndices)) then
         call io_write_var(gid, 'specialPointIndices', qpts%specialPointIndices(:qpts%numSpecialPoints))
      end if
      call io_gclose(gid)
   end subroutine write_header

   function split_cmplx1(c) result(r)
      complex, intent(in) :: c(:)
      real                :: r(2, size(c))

      r(1,:) = real(c)
      r(2,:) = aimag(c)
   end function split_cmplx1

   function split_cmplx2(c) result(r)
      complex, intent(in) :: c(:,:)
      real                :: r(2, size(c, 1), size(c, 2))

      r(1,:,:) = real(c)
      r(2,:,:) = aimag(c)
   end function split_cmplx2
#endif

end module m_dfpt_dynmat_io
