module m_types_coul
   use m_types_mat
   use m_mtir_size
   use m_types_fleurinput
   use m_judft
#ifdef CPP_MPI
   use mpi
#endif
   implicit none
   private
   type t_coul
      REAL, ALLOCATABLE      :: mt1_r(:, :, :, :)
      REAL, ALLOCATABLE      :: mt2_r(:, :, :, :), mt3_r(:, :, :)
      COMPLEX, ALLOCATABLE   :: mt1_c(:, :, :, :)
      COMPLEX, ALLOCATABLE   :: mt2_c(:, :, :, :), mt3_c(:, :, :)
      type(t_mat)            :: mtir
#ifdef CPP_MPI
      integer                :: comm = MPI_COMM_NULL ! communicator for this coulomb matrix
#else
      integer                :: comm = -1
#endif
      logical                :: l_participate = .False. ! am i somehow involved with this coulomb mtx
   contains
      procedure :: alloc => t_coul_alloc
      procedure :: mini_alloc => t_coul_mini_alloc
      procedure :: free => t_coul_free
      procedure :: mpi_bcast => t_coul_mpi_bc
   end type t_coul
   public t_coul
contains

   subroutine t_coul_mpi_bc(coul, fi, communicator, root)
      use m_types_fleurinput
      use m_types_hybmpi
      use m_judft

      implicit none
      class(t_coul)                  :: coul
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: root, communicator
#ifdef CPP_MPI
      integer :: ierr

      call timestart("Bcast coulomb_mtx")
      call timestart("Bcast small stuff")
      if (fi%sym%invs) THEN
         call MPI_Bcast(coul%mt1_r, size(coul%mt1_r), MPI_DOUBLE_PRECISION, root, communicator, ierr)
         call MPI_Bcast(coul%mt2_r, size(coul%mt2_r), MPI_DOUBLE_PRECISION, root, communicator, ierr)
         call MPI_Bcast(coul%mt3_r, size(coul%mt3_r), MPI_DOUBLE_PRECISION, root, communicator, ierr)
      else
         call MPI_Bcast(coul%mt1_c, size(coul%mt1_c), MPI_DOUBLE_COMPLEX, root, communicator, ierr)
         call MPI_Bcast(coul%mt2_c, size(coul%mt2_c), MPI_DOUBLE_COMPLEX, root, communicator, ierr)
         call MPI_Bcast(coul%mt3_c, size(coul%mt3_c), MPI_DOUBLE_COMPLEX, root, communicator, ierr)
      endif
      call timestop("Bcast small stuff")

      call coul%mtir%bcast(root, communicator)
      call timestop("Bcast coulomb_mtx")
#endif
   end subroutine t_coul_mpi_bc

   subroutine t_coul_free(coul)
      implicit none
      class(t_coul), intent(inout) :: coul
      integer :: ierr

      if (allocated(coul%mt1_r)) deallocate (coul%mt1_r)
      if (allocated(coul%mt1_c)) deallocate (coul%mt1_c)
      if (allocated(coul%mt2_r)) deallocate (coul%mt2_r)
      if (allocated(coul%mt3_r)) deallocate (coul%mt3_r)
      if (allocated(coul%mt2_c)) deallocate (coul%mt2_c)
      if (allocated(coul%mt3_c)) deallocate (coul%mt3_c)
      call coul%mtir%free()

#ifdef CPP_MPI
      if (coul%comm /= MPI_COMM_NULL) call MPI_comm_free(coul%comm, ierr)
#endif
   end subroutine t_coul_free

   subroutine t_coul_alloc(coul, fi, num_radbasfn, n_g, ikpt, l_print)
      implicit NONE
      class(t_coul), intent(inout) :: coul
      type(t_fleurinput), intent(in)    :: fi
      integer, intent(in) :: num_radbasfn(:, :), n_g(:), ikpt
      logical, intent(in), optional :: l_print
      integer :: info, isize

      isize = mtir_size(fi, n_g, ikpt)

      if (present(l_print)) then
         if (l_print) then
            write (*, *) "Coulomb dimensions:"
            write (*, *) "real:", fi%sym%invs
            write (*, *) "mt1 -> ["//int2str(maxval(num_radbasfn) - 1)// &
               ", "//int2str(maxval(num_radbasfn) - 1)// &
               ", "//int2str(maxval(fi%hybinp%lcutm1) + 1)// &
               ", "//int2str(fi%atoms%ntype)//"]"
            write (*, *) "mt2 -> ["//int2str(maxval(num_radbasfn) - 1)// &
               ", "//int2str(2*maxval(fi%hybinp%lcutm1) + 1)// &
               ", "//int2str(maxval(fi%hybinp%lcutm1) + 2)// &
               ", "//int2str(fi%atoms%nat)//"]"
            write (*, *) "mt3 -> ["//int2str(maxval(num_radbasfn) - 1)// &
               ", "//int2str(fi%atoms%nat)// &
               ", "//int2str(fi%atoms%nat)//"]"
            write (*, *) "mtir-> ["//int2str(isize)// &
               ", "//int2str(isize)//"]"
         endif
      endif

      if (fi%sym%invs) THEN
         if (.not. allocated(coul%mt1_r)) then
            allocate (coul%mt1_r(maxval(num_radbasfn) - 1, &
                                 maxval(num_radbasfn) - 1, &
                                 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=info)
            if (info /= 0) call judft_error("Can't allocate coul%mt1_r")
         endif
         if (.not. allocated(coul%mt2_r)) then
            allocate (coul%mt2_r(maxval(num_radbasfn) - 1, &
                                 -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                                 0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info)
            if (info /= 0) call judft_error("Can't allocate coul%mt2_r")
         endif

         if (.not. allocated(coul%mt3_r)) then
            allocate (coul%mt3_r(maxval(num_radbasfn) - 1, fi%atoms%nat, fi%atoms%nat), stat=info)
            if (info /= 0) call judft_error("Can't allocate coul%mt3_r")
         endif
         if (.not. coul%mtir%allocated()) call coul%mtir%alloc(.True., isize, isize)
      else
         if (.not. allocated(coul%mt1_c)) then
            allocate (coul%mt1_c(maxval(num_radbasfn) - 1, &
                                 maxval(num_radbasfn) - 1, &
                                 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=info)
            if (info /= 0) call judft_error("Can't allocate coul%mt1_c")
         endif
         if (.not. allocated(coul%mt2_c)) then
            allocate (coul%mt2_c(maxval(num_radbasfn) - 1, &
                                 -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                                 0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info)
            if (info /= 0) call judft_error("Can't allocate coul%mt2_c")
         endif

         if (.not. allocated(coul%mt3_c)) then
            allocate (coul%mt3_c(maxval(num_radbasfn) - 1, fi%atoms%nat, fi%atoms%nat), stat=info)
            if (info /= 0) call judft_error("Can't allocate coul%mt3_c")
         endif
         if (.not. coul%mtir%allocated()) call coul%mtir%alloc(.False., isize, isize)
      endif
   end subroutine t_coul_alloc

   subroutine t_coul_mini_alloc(coul, fi)
      implicit NONE
      type(t_fleurinput), intent(in)  :: fi
      class(t_coul), intent(inout)    :: coul

      if (.not. allocated(coul%mt1_r)) then
         allocate (coul%mt1_r(1,1,1,1))
      endif
      if (.not. allocated(coul%mt2_r)) then
         allocate (coul%mt2_r(1, 1, 1, 1))
      endif

      if (.not. allocated(coul%mt3_r)) then
         allocate (coul%mt3_r(1, 1, 1))
      endif
      if (.not. allocated(coul%mt1_c)) then
         allocate (coul%mt1_c(1, 1, 1, 1))
      endif
      if (.not. allocated(coul%mt2_c)) then
         allocate (coul%mt2_c(1, 1, 1, 1))
      endif

      if (.not. allocated(coul%mt3_c)) then
         allocate (coul%mt3_c(1, 1, 1))
      endif

      call coul%mtir%alloc(fi%sym%invs, 1, 1)
   end subroutine t_coul_mini_alloc
end module m_types_coul
