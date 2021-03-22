module m_store_load_hybrid
#ifdef CPP_HDF
   USE hdf5
#endif
   use m_juDFT
   use m_types
   use m_mpi_bc_tool
   use m_juDFT
   use m_types_mpimat

   character(len=*), parameter :: hybstore_fname = "hybrid.h5"
   public store_hybrid_data, load_hybrid_data
#ifdef CPP_HDF
   private open_file, open_datasetr, write_int_2d, close_dataset, close_file
#endif
contains
   subroutine load_hybrid_data(fi, fmpi, hybdat, mpdata)
      use m_constants
      implicit none
      type(t_fleurinput), intent(in)     :: fi
      type(t_mpi), intent(in)            :: fmpi
      type(t_hybdat), intent(inout)      :: hybdat
      type(t_mpdata), intent(inout)      :: mpdata

      logical :: l_exist
      integer, allocatable :: dims(:)
      character(len=:), allocatable :: dset_name
      integer                       :: ierr, nk, jsp
      real, allocatable             :: tmp(:, :)

#ifdef CPP_HDF
      integer(HID_T)   :: dset_id
      INTEGER(HID_T)   :: file_id

      call timestart("load_hybrid_data")

      if (fmpi%is_root()) INQUIRE (file='hybrid.h5', exist=l_exist)
      call mpi_bc(l_exist, 0, fmpi%mpi_comm)

      if (.not. allocated(hybdat%v_x)) then
         IF (fmpi%n_size == 1) THEN
            ALLOCATE (t_mat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
         ELSE
            ALLOCATE (t_mpimat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
         END IF
      end if

      if (l_exist) then
         if (fmpi%is_root()) then
            call timestart("read part")
            file_id = open_file()

            dset_id = open_dataset(file_id, "nbands")
            if (.not. allocated(hybdat%nbands)) allocate (hybdat%nbands(fi%kpts%nkptf, fi%input%jspins))
            call read_int_2d(dset_id, hybdat%nbands)
            call close_dataset(dset_id)

            dset_id = open_dataset(file_id, "nobd")
            if (.not. allocated(hybdat%nobd)) allocate (hybdat%nobd(fi%kpts%nkptf, fi%input%jspins))
            call read_int_2d(dset_id, hybdat%nobd)
            call close_dataset(dset_id)

            do jsp = 1, fi%input%jspins
               do nk = 1, fi%kpts%nkpt
                  if (fi%sym%invs) then
                     dset_name = "vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                     dset_id = open_dataset(file_id, dset_name)
                     dims = get_dims(dset_id)
                     call hybdat%v_x(nk, jsp)%alloc(fi%sym%invs, dims(1), dims(2))
                     call read_dbl_2d(dset_id, hybdat%v_x(nk, jsp)%data_r)
                     call close_dataset(dset_id)
                  else
                     dset_name = "r_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                     dset_id = open_dataset(file_id, dset_name)

                     ! get dimensions and alloc space
                     dims = get_dims(dset_id)
                     call hybdat%v_x(nk, jsp)%alloc(fi%sym%invs, dims(1), dims(2))
                     allocate (tmp(dims(1), dims(2)), stat=ierr)
                     if (ierr /= 0) call juDFT_error("can't alloc tmp")

                     ! get real part
                     call read_dbl_2d(dset_id, tmp)
                     hybdat%v_x(nk, jsp)%data_c = tmp
                     call close_dataset(dset_id)

                     ! get complex part
                     dset_name = "c_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                     dset_id = open_dataset(file_id, dset_name)
                     call read_dbl_2d(dset_id, tmp)
                     hybdat%v_x(nk, jsp)%data_c = hybdat%v_x(nk, jsp)%data_c + ImagUnit*tmp
                     call close_dataset(dset_id)
                     deallocate (tmp)
                  end if
               end do
            end do
            call timestop("read part")
         end if

         call timestart("bcast part")
         call mpi_bc(hybdat%nbands, 0, fmpi%mpi_comm)
         call mpi_bc(hybdat%nobd, 0, fmpi%mpi_comm)
         do jsp = 1, fi%input%jspins
            do nk = 1, fi%kpts%nkpt
               call hybdat%v_x(nk, jsp)%bcast(0, fmpi%mpi_comm)
            end do
         end do
         call timestop("bcast part")

         call mpdata%set_num_radfun_per_l(fi%atoms)
         call hybdat%set_maxlmindx(fi%atoms, mpdata%num_radfun_per_l)

         hybdat%l_addhf = .True.
         hybdat%l_subvxc = .True.
      end if

      call timestop("load_hybrid_data")
#endif
   end subroutine load_hybrid_data

   subroutine store_hybrid_data(fi, hybdat)
      implicit none
      type(t_fleurinput), intent(in)     :: fi
      type(t_hybdat), intent(in)         :: hybdat

      integer                       :: error, nk, jsp
      character(len=:), allocatable :: dset_name
#ifdef CPP_HDF
      integer(HID_T)   :: dset_id
      INTEGER(HID_T)   :: file_id

      call timestart("store_hybrid_data")
      file_id = open_file()

      dset_id = open_dataset(file_id, "nbands", [fi%kpts%nkptf, fi%input%jspins], H5T_NATIVE_INTEGER)
      call write_int_2d(dset_id, hybdat%nbands)
      call close_dataset(dset_id)

      dset_id = open_dataset(file_id, "nobd", [fi%kpts%nkptf, fi%input%jspins], H5T_NATIVE_INTEGER)
      call write_int_2d(dset_id, hybdat%nobd)
      call close_dataset(dset_id)

      dset_id = open_dataset(file_id, "bkf", [3, fi%kpts%nkptf], H5T_NATIVE_DOUBLE)
      call write_dbl_2d(dset_id, fi%kpts%bkf)
      call close_dataset(dset_id)

      dset_id = open_dataset(file_id, "bkp", [fi%kpts%nkptf, 1], H5T_NATIVE_INTEGER)
      call write_int_1d(dset_id, fi%kpts%bkp)
      call close_dataset(dset_id)

      ! hdf5 only knows reals
      do jsp = 1, fi%input%jspins
         do nk = 1, fi%kpts%nkpt
            if (fi%sym%invs) then
               dset_name = "vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
               dset_id = open_dataset(file_id, dset_name, shape(hybdat%v_x(nk, jsp)%data_r), H5T_NATIVE_DOUBLE)
               call write_dbl_2d(dset_id, hybdat%v_x(nk, jsp)%data_r)
               call close_dataset(dset_id)
            else
               dset_name = "r_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
               dset_id = open_dataset(file_id, dset_name, shape(hybdat%v_x(nk, jsp)%data_c), H5T_NATIVE_DOUBLE)
               call write_dbl_2d(dset_id, real(hybdat%v_x(nk, jsp)%data_c))
               call close_dataset(dset_id)

               dset_name = "c_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
               dset_id = open_dataset(file_id, dset_name, shape(hybdat%v_x(nk, jsp)%data_c), H5T_NATIVE_DOUBLE)
               call write_dbl_2d(dset_id, aimag(hybdat%v_x(nk, jsp)%data_c))
               call close_dataset(dset_id)
            end if
         end do
      end do

      call close_file(file_id)
      call timestop("store_hybrid_data")
#endif
   end subroutine store_hybrid_data

#ifdef CPP_HDF
   function open_file() result(file_id)
      implicit none
      integer(HID_T) :: file_id

      logical :: file_exist
      integer :: error

      INQUIRE (file='hybrid.h5', exist=file_exist)

      if (file_exist) then
         CALL h5fopen_f(hybstore_fname, H5F_ACC_RDWR_F, file_id, error)
         if (error /= 0) call juDFT_error("cant't open hdf5 file")
      else
         CALL h5fcreate_f(hybstore_fname, H5F_ACC_TRUNC_F, file_id, error)
         if (error /= 0) call juDFT_error("cant't create hdf5 file")
      end if
   end function open_file

   subroutine close_file(file_id)
      implicit none
      integer(HID_T), intent(in) :: file_id

      integer :: ierr

      CALL h5fclose_f(file_id, ierr)
      if (ierr /= 0) call juDFT_error("can't close hdf5 file")
   end subroutine close_file

   function open_dataset(file_id, dsetname, in_dims, type_id) result(dset_id)
      implicit NONE
      integer(HID_T), intent(in)           :: file_id
      character(len=*), intent(in)         :: dsetname
      integer, optional, intent(in)        :: in_dims(:)
      integer(HID_T), intent(in), optional :: type_id

      INTEGER(HID_T) :: dset_id
      integer :: ierr
      logical :: dset_exists
      integer(HID_T)   :: dspace_id
      INTEGER(HSIZE_T), allocatable :: dims(:)

      call h5lexists_f(file_id, dsetname, dset_exists, ierr)
      if (ierr /= 0) call juDFT_error("Can't check if dataset exists")

      if (dset_exists) then
         call h5dopen_f(file_id, dsetname, dset_id, ierr)
      else
         if (present(in_dims)) then
            allocate (dims(size(in_dims)))
            dims = in_dims
         else
            call juDFT_error("dims needed for file creation")
         end if

         CALL h5screate_simple_f(2, dims, dspace_id, ierr)
         if (ierr /= 0) call juDFT_error("can't create dataspace")

         if (.not. present(type_id)) call juDFT_error("type_id needed for dataset creation")
         CALL h5dcreate_f(file_id, dsetname, type_id, dspace_id, dset_id, ierr)
         if (ierr /= 0) call juDFT_error("creating data set failed")

         CALL h5sclose_f(dspace_id, ierr)
         if (ierr /= 0) call juDFT_error("can't close dataspace")
      end if
   end function open_dataset

   subroutine close_dataset(dset_id)
      implicit none
      INTEGER(HID_T), intent(in)         :: dset_id

      integer :: ierr

      CALL h5dclose_f(dset_id, ierr)
   end subroutine close_dataset

   subroutine write_int_1d(dset_id, mtx)
      implicit none
      INTEGER(HID_T), intent(in)         :: dset_id
      integer, intent(in)                :: mtx(:)

      INTEGER(HSIZE_T), DIMENSION(2)     :: data_dims
      integer                            :: ierr

      data_dims(1) = size(mtx)
      data_dims(2) = 1

      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mtx, data_dims, ierr)
      if (ierr /= 0) call juDFT_error("can't write int 1d")
   end subroutine write_int_1d

   subroutine write_int_2d(dset_id, mtx)
      implicit none
      INTEGER(HID_T), intent(in)         :: dset_id
      integer, intent(in)                :: mtx(:, :)

      INTEGER(HSIZE_T), DIMENSION(2)     :: data_dims
      integer                            :: ierr

      data_dims = shape(mtx)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mtx, data_dims, ierr)
      if (ierr /= 0) call juDFT_error("can't write int 2d")
   end subroutine write_int_2d

   subroutine read_int_2d(dset_id, mtx)
      implicit none
      INTEGER(HID_T), intent(in)         :: dset_id
      integer, intent(inout)             :: mtx(:, :)

      INTEGER(HSIZE_T), DIMENSION(2)     :: data_dims
      integer                            :: ierr

      data_dims = shape(mtx)
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, mtx, data_dims, ierr)
      if (ierr /= 0) call juDFT_error("can't read int 2d")
   end subroutine read_int_2d

   subroutine write_dbl_2d(dset_id, mtx)
      implicit none
      INTEGER(HID_T), intent(in)         :: dset_id
      real, intent(in)                   :: mtx(:, :)

      INTEGER(HSIZE_T), DIMENSION(2)     :: data_dims
      integer                            :: ierr

      data_dims = shape(mtx)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mtx, data_dims, ierr)
      if (ierr /= 0) call juDFT_error("can't write 2d real mtx")
   end subroutine write_dbl_2d

   subroutine read_dbl_2d(dset_id, mtx)
      implicit none
      INTEGER(HID_T), intent(in)     :: dset_id
      real, intent(inout)            :: mtx(:, :)

      INTEGER(HSIZE_T), DIMENSION(2)     :: data_dims
      integer                            :: ierr

      data_dims = shape(mtx)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mtx, data_dims, ierr)
      if (ierr /= 0) call juDFT_error("can't read int 2d")
   end subroutine read_dbl_2d

   function get_dims(dset_id) result(dims)
      implicit none
      INTEGER(HID_T), intent(in)         :: dset_id
      integer, allocatable               :: dims(:)

      integer(HSIZE_T), allocatable      :: hdims(:), hmaxdims(:)
      integer(HID_T) :: dataspace_id
      integer        :: ierr, ndims

      call h5dget_space_f(dset_id, dataspace_id, ierr)
      if (ierr /= 0) call juDFT_error("can't get dataspace")

      call h5sget_simple_extent_ndims_f(dataspace_id, ndims, ierr)
      if (ierr /= 0) call juDFT_error("can't get ndims")

      allocate (hdims(ndims), hmaxdims(ndims))

      call h5sget_simple_extent_dims_f(dataspace_id, hdims, hmaxdims, ierr)
      if (ierr /= ndims) call juDFT_error("can't get dims")

      dims = hdims

      call h5sclose_f(dataspace_id, ierr)
      if (ierr /= 0) call juDFT_error("can't close dataspace")
   end function get_dims

#endif
end module m_store_load_hybrid
