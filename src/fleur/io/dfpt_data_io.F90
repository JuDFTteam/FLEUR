!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_dfpt_data_io
   !! DFPT response tensors in dfpt_data.hdf, arrays in Fortran order:
   !!   /dielectricTensor        @qlim   tensor(3,3)       (field direction, field direction)
   !!   /dielectricTensorStatic          tensor(3,3)
   !!   /bornEffectiveCharges    @nAtoms charges(3,3,nat)  (displacement direction, field direction, atom)
#ifdef CPP_HDF
   use hdf5
   use m_hdf_tools
#endif
   use m_juDFT
   use m_types_atoms

   implicit none
   private
   public :: write_dfpt_dielten, read_dfpt_dielten, write_dfpt_dielten_static, write_dfpt_bec, read_dfpt_bec

   character(len=*), parameter :: DATA_FILE = "dfpt_data.hdf"
   character(len=*), parameter :: GRP_DIELTEN = "/dielectricTensor"
   character(len=*), parameter :: GRP_DIELTEN_STATIC = "/dielectricTensorStatic"
   character(len=*), parameter :: GRP_BEC = "/bornEffectiveCharges"

contains

   subroutine write_dfpt_dielten(dielten, qlim)
      !! High-frequency dielectric tensor.

      real, intent(in) :: dielten(3,3)
      real, intent(in) :: qlim

#ifdef CPP_HDF
      integer(HID_T) :: fid, gid

      call open_data_file(fid)
      call new_group(fid, GRP_DIELTEN, gid)
      call io_write_attreal0(gid, 'qlim', qlim)
      call io_write_var(gid, 'tensor', dielten)
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      call judft_error("DFPT data IO requires FLEUR compiled with HDF5.", calledby="write_dfpt_dielten")
#endif
   end subroutine write_dfpt_dielten

   subroutine read_dfpt_dielten(dielten)
      real, intent(out) :: dielten(3,3)

#ifdef CPP_HDF
      integer(HID_T) :: fid, gid

      call open_group(GRP_DIELTEN, "the E-field run (l_efield)", fid, gid)
      call io_read_var(gid, 'tensor', dielten)
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      dielten = 0.0
      call judft_error("DFPT data IO requires FLEUR compiled with HDF5.", calledby="read_dfpt_dielten")
#endif
   end subroutine read_dfpt_dielten

   subroutine write_dfpt_dielten_static(dielten)
      real, intent(in) :: dielten(3,3)

#ifdef CPP_HDF
      integer(HID_T) :: fid, gid

      call open_data_file(fid)
      call new_group(fid, GRP_DIELTEN_STATIC, gid)
      call io_write_var(gid, 'tensor', dielten)
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      call judft_error("DFPT data IO requires FLEUR compiled with HDF5.", calledby="write_dfpt_dielten_static")
#endif
   end subroutine write_dfpt_dielten_static

   subroutine write_dfpt_bec(atoms, bec)
      !! bec(iType, displacement direction, field direction) is stored per atom.

      type(t_atoms), intent(in) :: atoms
      real,          intent(in) :: bec(:,:,:)

#ifdef CPP_HDF
      integer(HID_T) :: fid, gid
      integer        :: iAtom
      real           :: charges(3,3,atoms%nat)

      do iAtom = 1, atoms%nat
         charges(:,:,iAtom) = bec(atoms%itype(iAtom),:,:)
      end do
      call open_data_file(fid)
      call new_group(fid, GRP_BEC, gid)
      call io_write_attint0(gid, 'nAtoms', atoms%nat)
      call io_write_var(gid, 'charges', charges)
      call io_gclose(gid)
      call io_hdfclose(fid)
#else
      call judft_error("DFPT data IO requires FLEUR compiled with HDF5.", calledby="write_dfpt_bec")
#endif
   end subroutine write_dfpt_bec

   subroutine read_dfpt_bec(atoms, bec)
      !! Fills bec(iType, displacement direction, field direction) from the per-atom data.

      type(t_atoms), intent(in)  :: atoms
      real,          intent(out) :: bec(:,:,:)

#ifdef CPP_HDF
      integer(HID_T) :: fid, gid
      integer        :: iAtom, nat_file
      real           :: charges(3,3,atoms%nat)

      bec = 0.0
      if (atoms%ntype /= atoms%nat) call judft_error("Born effective charges need one atom per atom type.", calledby="read_dfpt_bec")
      call open_group(GRP_BEC, "the Born effective charge run (l_borneffcharge)", fid, gid)
      call io_read_attint0(gid, 'nAtoms', nat_file)
      if (nat_file /= atoms%nat) call judft_error(DATA_FILE//" belongs to a cell with a different number of atoms.", calledby="read_dfpt_bec")
      call io_read_var(gid, 'charges', charges)
      call io_gclose(gid)
      call io_hdfclose(fid)
      do iAtom = 1, atoms%nat
         bec(atoms%itype(iAtom),:,:) = charges(:,:,iAtom)
      end do
#else
      bec = 0.0
      call judft_error("DFPT data IO requires FLEUR compiled with HDF5.", calledby="read_dfpt_bec")
#endif
   end subroutine read_dfpt_bec

#ifdef CPP_HDF
   subroutine open_data_file(fid)
      !! Opens dfpt_data.hdf for writing, creating it if missing.

      integer(HID_T), intent(out) :: fid

      logical :: l_exist
      integer :: hdferr

      inquire(file=DATA_FILE, exist=l_exist)
      if (l_exist) then
         call io_hdfopen(DATA_FILE, H5F_ACC_RDWR_F, fid, hdferr)
      else
         call h5fcreate_f(DATA_FILE, H5F_ACC_TRUNC_F, fid, hdferr)
      end if
      if (hdferr /= 0) call judft_error("Could not open "//DATA_FILE, calledby="open_data_file")
   end subroutine open_data_file

   subroutine new_group(fid, gname, gid)
      !! Creates gname, replacing an existing group of that name.

      integer(HID_T),   intent(in)  :: fid
      character(len=*), intent(in)  :: gname
      integer(HID_T),   intent(out) :: gid

      if (io_groupexists(fid, gname)) call io_gdelete(fid, gname)
      call io_gcreate(fid, gname, gid)
   end subroutine new_group

   subroutine open_group(gname, producer, fid, gid)
      !! Opens dfpt_data.hdf read-only and the group gname; errors name the run that writes it.

      character(len=*), intent(in)  :: gname, producer
      integer(HID_T),   intent(out) :: fid, gid

      logical :: l_exist
      integer :: hdferr

      inquire(file=DATA_FILE, exist=l_exist)
      if (.not. l_exist) call judft_error(DATA_FILE//" not found; it is written by "//producer//".", calledby="open_group")
      call io_hdfopen(DATA_FILE, H5F_ACC_RDONLY_F, fid, hdferr)
      if (hdferr /= 0) call judft_error("Could not open "//DATA_FILE, calledby="open_group")
      if (.not. io_groupexists(fid, gname)) then
         call judft_error(gname//" missing in "//DATA_FILE//"; it is written by "//producer//".", calledby="open_group")
      end if
      call io_gopen(fid, gname, gid)
   end subroutine open_group
#endif

end module m_dfpt_data_io
