module m_store_load_hybrid
#ifdef CPP_HDF
   USE hdf5
#endif
   use m_juDFT
   use m_types
   use m_mpi_bc_tool
   use m_juDFT
   use m_types_mpimat
   use m_distrib_vx

   character(len=*), parameter :: hybstore_fname = "hybrid.h5"
   public store_hybrid_data, load_hybrid_data
#ifdef CPP_HDF
   private open_file, open_datasetr, write_int_2d, close_dataset, close_file
#endif
contains
   subroutine load_hybrid_data(fi, fmpi, hybdat, mpdata)
      use m_constants
      use m_mixing_history
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
      type(t_mat)                   :: vx_tmp

#ifdef CPP_HDF
      integer(HID_T)   :: dset_id
      INTEGER(HID_T)   :: file_id


      if(.not. hybdat%l_subvxc) then
         call timestart("load_hybrid_data")
         if (fmpi%is_root()) INQUIRE (file='hybrid.h5', exist=l_exist)
         call mpi_bc(l_exist, 0, fmpi%mpi_comm)


         if(l_exist .and. (.not. allocated(hybdat%v_x))) then
#ifdef CPP_MPI
            call MPI_Barrier(fmpi%mpi_comm, ierr)
#endif
            call mixing_history_reset(fmpi)

            IF (fmpi%n_size == 1) THEN
               ALLOCATE (t_mat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
            ELSE
               ALLOCATE (t_mpimat::hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))
            END IF

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
               call timestop("read part")
            end if

            do jsp = 1, fi%input%jspins
               do nk = 1, fi%kpts%nkpt
                  if(fmpi%is_root()) then
                     call timestart("read part")
                     if (fi%sym%invs) then
                        dset_name = "vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                        dset_id = open_dataset(file_id, dset_name)
                        dims = get_dims(dset_id)
                        call vx_tmp%alloc(fi%sym%invs, dims(1), dims(2))
                        call read_dbl_2d(dset_id, vx_tmp%data_r)
                        call close_dataset(dset_id)
                     else
                        dset_name = "r_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                        dset_id = open_dataset(file_id, dset_name)

                        ! get dimensions and alloc space
                        dims = get_dims(dset_id)
                        call vx_tmp%alloc(fi%sym%invs, dims(1), dims(2))
                        allocate (tmp(dims(1), dims(2)), stat=ierr)
                        if (ierr /= 0) call juDFT_error("can't alloc tmp")

                        ! get real part
                        call read_dbl_2d(dset_id, tmp)
                        vx_tmp%data_c = tmp
                        call close_dataset(dset_id)

                        ! get complex part
                        dset_name = "c_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                        dset_id = open_dataset(file_id, dset_name)
                        call read_dbl_2d(dset_id, tmp)
                        vx_tmp%data_c = vx_tmp%data_c + ImagUnit*tmp
                        call close_dataset(dset_id)
                        deallocate (tmp)
                     end if
                     call timestop("read part")
                  endif

                  call mpi_bc(dims, 0, fmpi%mpi_comm)
                  call distrib_single_vx(fi, fmpi, jsp, nk, 0, vx_tmp, hybdat, dims=dims)
                  call vx_tmp%free()
               end do
            end do

            call timestart("bcast part")
            call mpi_bc(hybdat%nbands, 0, fmpi%mpi_comm)
            call mpi_bc(hybdat%nobd, 0, fmpi%mpi_comm)      
            call timestop("bcast part")

            call mpdata%set_num_radfun_per_l(fi%atoms)
            call hybdat%set_maxlmindx(fi%atoms, mpdata%num_radfun_per_l)

            hybdat%l_addhf = .True.
            hybdat%l_subvxc = .True.
         end if
         call timestop("load_hybrid_data")
      endif
#endif
   end subroutine load_hybrid_data

   subroutine store_hybrid_data(fi, fmpi, hybdat)
      implicit none
      type(t_fleurinput), intent(in)     :: fi
      type(t_mpi), intent(in)            :: fmpi
      type(t_hybdat), intent(in)         :: hybdat

      integer                       :: error, nk, jsp
      character(len=:), allocatable :: dset_name
      type(t_mat) :: vx_tmp
#ifdef CPP_HDF
      integer(HID_T)   :: dset_id
      INTEGER(HID_T)   :: file_id

      call timestart("store_hybrid_data")

      if(fmpi%irank == 0) then
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
      endif

      ! hdf5 only knows reals
      do jsp = 1, fi%input%jspins
         do nk = 1, fi%kpts%nkpt
            call collect_vx(fi, fmpi, hybdat, nk, jsp, vx_tmp)
            if(fmpi%irank == 0 ) then
               if (fi%sym%invs) then
                  dset_name = "vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                  dset_id = open_dataset(file_id, dset_name, shape(vx_tmp%data_r), H5T_NATIVE_DOUBLE)
                  call write_dbl_2d(dset_id, vx_tmp%data_r)
                  call close_dataset(dset_id)
               else
                  dset_name = "r_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                  dset_id = open_dataset(file_id, dset_name, shape(vx_tmp%data_c), H5T_NATIVE_DOUBLE)
                  call write_dbl_2d(dset_id, real(vx_tmp%data_c))
                  call close_dataset(dset_id)

                  dset_name = "c_vx_nk="//int2str(nk)//"_jsp="//int2str(jsp)
                  dset_id = open_dataset(file_id, dset_name, shape(vx_tmp%data_c), H5T_NATIVE_DOUBLE)
                  call write_dbl_2d(dset_id, aimag(vx_tmp%data_c))
                  call close_dataset(dset_id)
               end if
               call vx_tmp%free()
            endif
         end do
      end do

      if(fmpi%irank == 0) call close_file(file_id)
      call timestop("store_hybrid_data")
#endif
   end subroutine store_hybrid_data

   subroutine collect_vx(fi, fmpi, hybdat, nk, jsp, vx_tmp)
      use m_glob_tofrom_loc
      implicit none 
      type(t_fleurinput), intent(in)     :: fi
      type(t_mpi), intent(in)            :: fmpi
      type(t_hybdat), intent(in)         :: hybdat
      integer, intent(in)                :: nk, jsp 
      type(t_mat), intent(inout)         :: vx_tmp
      
      integer, parameter :: recver = 0 ! HDF is node on global root
      integer :: sender, ierr, buff(2), i, pe_i, i_loc
      logical :: l_mpimat

      ! find out and bcast what kind of matrix we are using
      sender = merge(fmpi%irank, -1, any(fmpi%k_list == nk) .and. fmpi%n_rank == 0)
#ifdef CPP_MPI
      call MPI_Allreduce(MPI_IN_PLACE, sender, 1, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
#endif

      select type(vx_origin => hybdat%v_x(nk,jsp)) 
      class is(t_mpimat) 
         if(sender == fmpi%irank) then
            buff = [vx_origin%global_size1, vx_origin%global_size2]
         endif
      class is(t_mat)
         if(sender == fmpi%irank) then
            buff = [vx_origin%matsize1, vx_origin%matsize2]
         endif
      end select
#ifdef CPP_MPI
      call MPI_Bcast(buff, 2, MPI_INTEGER, sender, fmpi%mpi_comm, ierr)
#endif

      if(fmpi%irank == recver) then
         call vx_tmp%init(fi%sym%invs, buff(1), buff(2))
      endif

      do i = 1, buff(2)
         call glob_to_loc(fmpi, i, pe_i, i_loc)
         sender = merge(fmpi%irank, -1, pe_i == fmpi%n_rank .and. any(fmpi%k_list == nk))
#ifdef CPP_MPI
         call MPI_Allreduce(MPI_IN_PLACE, sender, 1, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
#endif

         if(sender == recver .and. fmpi%irank == recver) then
            if(vx_tmp%l_real) then
               vx_tmp%data_r(:,i) = hybdat%v_x(nk,jsp)%data_r(:,i_loc)
            else
               vx_tmp%data_c(:,i) = hybdat%v_x(nk,jsp)%data_c(:,i_loc)
            endif
#ifdef CPP_MPI
         elseif(sender == fmpi%irank) then
            if(fi%sym%invs) then
               call MPI_Send(hybdat%v_x(nk,jsp)%data_r(:,i_loc), buff(1), MPI_DOUBLE_PRECISION, recver, 100+i, fmpi%mpi_comm,ierr)
            else
               call MPI_Send(hybdat%v_x(nk,jsp)%data_c(:,i_loc), buff(1), MPI_DOUBLE_COMPLEX, recver, 100+i, fmpi%mpi_comm,ierr)
            endif
         elseif(fmpi%irank == recver) then
            if(vx_tmp%l_real) then
               call MPI_Recv(vx_tmp%data_r(:,i), buff(1), MPI_DOUBLE_PRECISION, sender, 100+i, fmpi%mpi_comm, MPI_STATUS_IGNORE, ierr)
            else
               call MPI_Recv(vx_tmp%data_c(:,i), buff(1), MPI_DOUBLE_COMPLEX, sender, 100+i, fmpi%mpi_comm, MPI_STATUS_IGNORE, ierr)
            endif
#endif
         endif
      enddo
   end subroutine collect_vx

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
