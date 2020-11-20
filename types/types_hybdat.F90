MODULE m_types_hybdat
   use m_types_usdus
   use m_types_mat
#ifdef CPP_MPI
   use mpi 
#endif
   IMPLICIT NONE

   type t_mtir_block
      integer, allocatable :: pe(:)   ! list where each coulomb is stored
#ifdef CPP_MPI
      integer(MPI_ADDRESS_KIND), allocatable :: slot(:) ! local slot where it's stored
#else 
      integer, allocatable :: slot(:) ! local slot where it's stored
#endif
      logical              :: l_real  ! real or complex storage

      real, pointer    :: r(:,:,:) ! (:,:,i_slot)
      complex, pointer :: c(:,:,:) ! (:,:,i_slot)

      !mpi stuff 
      integer :: handle
   contains 
      procedure :: alloc => t_mtir_block_alloc
      procedure :: init  => t_mtir_block_init
      procedure :: read  => t_mtir_block_read
   end type t_mtir_block

   type t_coul
      REAL, ALLOCATABLE      :: mt1_r(:, :, :, :)
      REAL, ALLOCATABLE      :: mt2_r(:, :, :, :),   mt3_r(:, :, :)
      REAL, ALLOCATABLE      :: mtir_r(:, :)
      COMPLEX, ALLOCATABLE   :: mt1_c(:, :, :, :)
      COMPLEX, ALLOCATABLE   :: mt2_c(:, :, :, :),   mt3_c(:, :, :)
      COMPLEX, ALLOCATABLE   :: mtir_c(:, :)
   contains 
      procedure :: init     => t_coul_init
      procedure :: alloc    => t_coul_alloc
      procedure :: free     => t_coul_free
      procedure :: mpi_bc  => t_coul_mpi_bc
      procedure :: size_MB  => t_coul_size_MB
   end type t_coul

   TYPE t_hybdat
      LOGICAL                :: l_subvxc = .false.
      LOGICAL                :: l_calhf = .false.
      LOGICAL                :: l_addhf = .false.
      LOGICAL                :: l_print_iob_splitting = .True.
      INTEGER                :: lmaxcd, maxindxc
      INTEGER                :: maxfac
      REAL, ALLOCATABLE      :: gridf(:, :)
      INTEGER, ALLOCATABLE   :: nindxc(:, :)
      INTEGER, ALLOCATABLE   :: lmaxc(:)
      REAL, ALLOCATABLE      :: core1(:, :, :, :), core2(:, :, :, :)
      REAL, ALLOCATABLE      :: eig_c(:, :, :)
      !INTEGER, ALLOCATABLE   :: kveclo_eig(:, :)
      REAL, ALLOCATABLE      :: sfac(:), fac(:)
      REAL, ALLOCATABLE      :: gauntarr(:, :, :, :, :, :)
      REAL, ALLOCATABLE      :: bas1(:, :, :, :), bas2(:, :, :, :)
      REAL, ALLOCATABLE      :: bas1_MT(:, :, :), drbas1_MT(:, :, :)
      REAL, ALLOCATABLE      :: prodm(:, :, :, :)
      REAL, ALLOCATABLE      :: div_vv(:, :, :)
      INTEGER, ALLOCATABLE   :: pntgptd(:)
      INTEGER, ALLOCATABLE   :: pntgpt(:, :, :, :)
      INTEGER, ALLOCATABLE   :: nindxp1(:, :)
      INTEGER, ALLOCATABLE   :: ne_eig(:)
      INTEGER, ALLOCATABLE   :: nbands(:)
      INTEGER, ALLOCATABLE   :: nobd(:, :)
      INTEGER                :: maxlmindx = -1
      COMPLEX, ALLOCATABLE   :: stepfunc(:, :, :)
      INTEGER                :: nbasp = -1
      INTEGER                :: eig_id = -1
      INTEGER, ALLOCATABLE   :: nbasm(:)

      ! coulomb matrix stuff
      type(t_coul), allocatable     :: coul(:)
      type(t_mtir_block)              :: mtir

      type(t_usdus)             :: usdus
      type(t_mat), allocatable  :: v_x(:,:) ! nkpt, jsp
   contains
      procedure :: set_stepfunction => set_stepfunction
      procedure :: free       => free_hybdat
      procedure :: allocate   => allocate_hybdat
      procedure :: set_states => set_states_hybdat
   END TYPE t_hybdat

contains
   subroutine t_mtir_block_read(mtir, fi, n_g, ik, out_mtx)
      use m_types_fleurinput
      use m_types_mat 
      use m_types_mpi
      use m_judft
      implicit none
      class(t_mtir_block), intent(in)   :: mtir
      type(t_fleurinput), intent(in)    :: fi
      integer, intent(in)               :: ik, n_g(:)
      type(t_mat), intent(inout)        :: out_mtx 

      integer :: isize, irank, ierr, max_size, i

      call timestart("t_mtir_block_read")
      isize = mtir_size(fi, n_g, ik)
      max_size = 0
      do i = 1,fi%kpts%nkpt
         max_size = max(max_size, mtir_size(fi, n_g, i))
      enddo

#ifdef CPP_MPI 
      call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
#else 
      irank = 0
#endif
      call out_mtx%init(mtir%l_real, max_size, max_size)
      ! check if I have it locally
      if(mtir%pe(ik) == irank) then 
         if(mtir%l_real) then 
            call dlacpy("A", isize, isize, mtir%r(1,1,mtir%slot(ik)+1), size(mtir%r, 1), out_mtx%data_r, max_size)
         else
            call zlacpy("A", isize, isize, mtir%c(1,1,mtir%slot(ik)+1), size(mtir%c, 1), out_mtx%data_c, max_size)
         endif
      else
#ifdef CPP_MPI 
         CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, mtir%pe(ik), 0, mtir%handle, ierr)

         if(mtir%l_real) then 
            CALL MPI_GET(out_mtx%data_r, max_size**2, MPI_DOUBLE_PRECISION, mtir%pe(ik), mtir%slot(ik),&
                         max_size**2, MPI_DOUBLE_PRECISION, mtir%handle, ierr)
         else
            CALL MPI_GET(out_mtx%data_c, max_size**2, MPI_DOUBLE_COMPLEX, mtir%pe(ik), mtir%slot(ik),&
                         max_size**2, MPI_DOUBLE_COMPLEX, mtir%handle, ierr)
         endif 

         CALL MPI_WIN_UNLOCK(mtir%pe(ik), mtir%handle, ierr)
#else 
         call judft_error("without MPI everything should be local (but it isn't i guess)")
#endif
      endif

      call timestop("t_mtir_block_read")
   end subroutine t_mtir_block_read

   subroutine t_mtir_block_init(mtir, fi, fmpi)
      use m_types_fleurinput
      use m_types_mpi
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none
      class(t_mtir_block), intent(inout)  :: mtir
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpi), INTENT(IN)           :: fmpi

      integer :: i, cnt_slot, ierr
      integer, allocatable :: tmp_slot(:)

      call timestart("t_mtir_block_init")

#ifdef CPP_MPI
      mtir%pe   = -1 
      allocate(tmp_slot(fi%kpts%nkpt), source=-1)

      cnt_slot = 0
      do i = 1,fi%kpts%nkpt 
         if(any(i == fmpi%k_list .and. fmpi%n_rank == 0)) then 
            mtir%pe(i)   = fmpi%irank
            tmp_slot(i) = cnt_slot 
            cnt_slot = cnt_slot + 1
         endif 
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, mtir%pe,   fi%kpts%nkpt, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, tmp_slot, fi%kpts%nkpt, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
      mtir%slot = tmp_slot
#else 
       mtir%pe = 0 
       mtir%slot = [(i, i=1,fi%kpts%nkpt)]
#endif
      call timestop("t_mtir_block_init")
   end subroutine t_mtir_block_init

   subroutine t_mtir_block_alloc(mtir, fi, fmpi, n_g, my_n_k)
      use m_types_mpi
      use m_types_fleurinput
      use m_judft
      implicit none 
      class(t_mtir_block), intent(inout)  :: mtir
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      integer, intent(in)               :: n_g(:), my_n_k

      integer :: ikpt, max_coul_size, ierr, slot_size, irank, type_size
      integer(kind=MPI_ADDRESS_KIND) :: win_size

      call timestart("t_mtir_block_alloc")

      mtir%l_real = fi%sym%invs 

      max_coul_size = 0
      do ikpt = 1,fi%kpts%nkpt
         max_coul_size = max(max_coul_size, mtir_size(fi, n_g, ikpt))
      enddo

#ifdef CPP_MPI 
      call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
      if(mtir%l_real) then
         CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, type_size, ierr)
      else
         CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX, type_size, ierr)
      endif
      slot_size = type_size * max_coul_size**2
      win_size  = slot_size * my_n_k
#else 
      irank = 0
#endif

      write (*,*) "max_coul_size on create", max_coul_size

      if(mtir%l_real) then        
         if(.not. associated(mtir%r)) then
            allocate(mtir%r(max_coul_size, max_coul_size, my_n_k), stat=ierr)
            if(ierr /= 0) call juDFT_error("can't allocate mtir%r of size: " // &
                                              int2str(max_coul_size) // "^2 x " // int2str(my_n_k))

#ifdef CPP_MPI 
            if(fmpi%isize > 1) then
               write (*,*) "[" // int2str(irank) //"]: slot_size, win_size", slot_size, win_size
               call judft_win_create_real(mtir%r, win_size, slot_size, &
                                    MPI_INFO_NULL, fmpi%mpi_comm, mtir%handle)
            endif
#endif
         endif

      else
         if(.not. associated(mtir%c)) then
            allocate(mtir%c(max_coul_size, max_coul_size, my_n_k), stat=ierr)
            if(ierr /= 0) call juDFT_error("can't allocate mtir%c of size: " // &
                                             int2str(max_coul_size) // "^2 x " // int2str(my_n_k)) 

#ifdef CPP_MPI 
            if(fmpi%isize > 1) then
               write (*,*) "[" // int2str(irank) //"]: slot_size, win_size", slot_size, win_size
               call judft_win_create_cmplx(mtir%c, win_size, slot_size, &
                                    MPI_INFO_NULL, fmpi%mpi_comm, mtir%handle)
            endif
#endif
         endif
      endif

      if(.not. allocated(mtir%pe)) then
         allocate(mtir%pe(fi%kpts%nkpt), stat=ierr)
         if(ierr /= 0) call judft_error("can't alloc mtir%pe")
      endif 

      if(.not. allocated(mtir%slot)) then
         allocate(mtir%slot(fi%kpts%nkpt), stat=ierr)
         if(ierr /= 0) call judft_error("can't alloc mtir%slot")
      endif

      call timestop("t_mtir_block_alloc")
   end subroutine t_mtir_block_alloc

   subroutine set_states_hybdat(hybdat, fi, results, jsp)
      use m_judft
      use m_types_misc
      use m_types_fleurinput
      implicit none 
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      type(t_results), intent(in)    :: results
      integer, intent(in)            :: jsp

      integer :: nk, i

      DO nk = 1, fi%kpts%nkpt, 1
         hybdat%ne_eig(nk) = results%neig(nk, jsp)
         hybdat%nobd(nk,jsp) = COUNT(results%w_iks(:hybdat%ne_eig(nk), nk, jsp) > 0.0)
         hybdat%nbands(nk) = fi%hybinp%bands1
      END do

      DO nk = 1, fi%kpts%nkpt
         IF (hybdat%nobd(nk,jsp) > hybdat%nbands(nk)) THEN
            hybdat%nbands(nk) = hybdat%nobd(nk,jsp)
         END IF
      ENDDO
      
      ! spread hybdat%nobd from IBZ to whole BZ
      DO nk = 1, fi%kpts%nkptf
         i = fi%kpts%bkp(nk)
         hybdat%nbands(nk) = hybdat%nbands(i)
         hybdat%nobd(nk,jsp) = hybdat%nobd(i,jsp)
      END DO
   end subroutine set_states_hybdat

   function t_coul_size_MB(coul) result(size_MB)
      implicit none 
      class(t_coul), intent(in) :: coul
      real  :: size_MB 

      size_MB = 0
      
      ! real parts
      if(allocated(coul%mt1_r))   size_MB = size_MB + 8 * 1e-6 * size(coul%mt1_r) 
      if(allocated(coul%mt2_r))   size_MB = size_MB + 8 * 1e-6 * size(coul%mt2_r) 
      if(allocated(coul%mt3_r))   size_MB = size_MB + 8 * 1e-6 * size(coul%mt3_r) 
      if(allocated(coul%mtir_r))  size_MB = size_MB + 8 * 1e-6 * size(coul%mtir_r) 

      ! complex parts
      if(allocated(coul%mt1_c))   size_MB = size_MB + 16 * 1e-6 * size(coul%mt1_c) 
      if(allocated(coul%mt2_c))   size_MB = size_MB + 16 * 1e-6 * size(coul%mt2_r) 
      if(allocated(coul%mt3_c))   size_MB = size_MB + 16 * 1e-6 * size(coul%mt3_r) 
      if(allocated(coul%mtir_c))  size_MB = size_MB + 16 * 1e-6 * size(coul%mtir_r) 
   end function

   subroutine t_coul_mpi_bc(coul, fi, communicator, root)
      use m_types_fleurinput
      use m_types_hybmpi
      use m_judft
#ifdef CPP_MPI
      use mpi
#endif
      implicit none 
      class(t_coul)                  :: coul
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: root, communicator
#ifdef CPP_MPI
      integer :: ierr

      if (fi%sym%invs) THEN
         call MPI_Bcast(coul%mt1_r, size(coul%mt1_r), MPI_DOUBLE_PRECISION, root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mt1_r failed")

         call MPI_Bcast(coul%mt2_r,   size(coul%mt2_r),   MPI_DOUBLE_PRECISION, root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mt2_r failed")

         call MPI_Bcast(coul%mt3_r,   size(coul%mt3_r),   MPI_DOUBLE_PRECISION, root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mt3_r failed")

         call MPI_Bcast(coul%mtir_r,  size(coul%mtir_r),  MPI_DOUBLE_PRECISION, root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mtir_r failed")
      else 
         call MPI_Bcast(coul%mt1_c, size(coul%mt1_c),     MPI_DOUBLE_COMPLEX,  root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mt1_c failed")

         call MPI_Bcast(coul%mt2_c,   size(coul%mt2_c),   MPI_DOUBLE_COMPLEX , root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mt2_r failed")

         call MPI_Bcast(coul%mt3_c,   size(coul%mt3_c),   MPI_DOUBLE_COMPLEX , root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mt3_r failed")

         call MPI_Bcast(coul%mtir_c,  size(coul%mtir_c),  MPI_DOUBLE_COMPLEX , root, communicator, ierr)
         if(ierr /= 0) call judft_error("MPI_Bcast of coul%mtir_r failed")
      endif
#endif
   end subroutine t_coul_mpi_bc

   subroutine t_coul_free(coul)
      implicit none 
      class(t_coul), intent(inout) :: coul 

      if(allocated(coul%mt1_r)) deallocate(coul%mt1_r)
      if(allocated(coul%mt1_c)) deallocate(coul%mt1_c)
      if(allocated(coul%mt2_r)) deallocate(coul%mt2_r)
      if(allocated(coul%mt3_r)) deallocate(coul%mt3_r)
      if(allocated(coul%mtir_r)) deallocate(coul%mtir_r)
      if(allocated(coul%mt2_c)) deallocate(coul%mt2_c)
      if(allocated(coul%mt3_c)) deallocate(coul%mt3_c)
      if(allocated(coul%mtir_c)) deallocate(coul%mtir_c)
   end subroutine t_coul_free

   subroutine t_coul_alloc(coul, fi, num_radbasfn, n_g, ikpt, l_print)
      use m_types_fleurinput
      use m_judft
      implicit NONE 
      class(t_coul), intent(inout) :: coul
      type(t_fleurinput), intent(in)    :: fi
      integer, intent(in) :: num_radbasfn(:, :), n_g(:), ikpt
      logical, intent(in), optional :: l_print
      integer :: info, isize, l, itype

      isize = mtir_size(fi, n_g, ikpt)

      if(present(l_print)) then 
         if(l_print) then 
            write (*,*) "Coulomb dimensions:"
            write (*,*) "real:", fi%sym%invs
            write (*,*) "mt1 -> [" // int2str(maxval(num_radbasfn) - 1) //&
                              ", "  // int2str(maxval(num_radbasfn) - 1) //&
                              ", "  // int2str(maxval(fi%hybinp%lcutm1)+1)// &
                              ", "  // int2str(fi%atoms%ntype) //  "]"
            write (*,*) "mt2 -> [" // int2str(maxval(num_radbasfn) - 1) // &
                              ", " // int2str(2*maxval(fi%hybinp%lcutm1)+1) // &
                              ", " // int2str(maxval(fi%hybinp%lcutm1)+2) // &
                              ", " // int2str(fi%atoms%nat) // "]"
            write (*,*) "mt3 -> [" // int2str(maxval(num_radbasfn) - 1) // &
                              ", " // int2str(fi%atoms%nat) // &
                              ", " // int2str(fi%atoms%nat) // "]"
            write (*,*) "mtir-> [" // int2str(isize) // &
                              ", " // int2str(isize) // "]"
         endif
      endif
 
      if (fi%sym%invs) THEN      
         if(.not. allocated(coul%mt1_r)) then 
            allocate(coul%mt1_r(maxval(num_radbasfn) - 1,&
                              maxval(num_radbasfn) - 1,&
                              0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt1_r")
         endif
         if(.not. allocated(coul%mt2_r)) then
            allocate(coul%mt2_r(maxval(num_radbasfn) - 1,&
                                 -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1),&
                                 0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt2_r")
         endif

         if(.not. allocated(coul%mt3_r)) then
             allocate(coul%mt3_r(maxval(num_radbasfn) - 1,fi%atoms%nat, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt3_r")
         endif 

         if(.not. allocated(coul%mtir_r)) then
            allocate(coul%mtir_r(isize,isize), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mtir_r")
         endif 
      else      
         if(.not. allocated(coul%mt1_c)) then 
            allocate(coul%mt1_c(maxval(num_radbasfn) - 1,&
                              maxval(num_radbasfn) - 1,&
                              0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt1_c")
         endif
         if(.not. allocated(coul%mt2_c)) then
            allocate(coul%mt2_c(maxval(num_radbasfn) - 1,&
                                 -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                                 0:maxval(fi%hybinp%lcutm1) + 1, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt2_c")
         endif 

         if(.not. allocated(coul%mt3_c)) then
            allocate(coul%mt3_c(maxval(num_radbasfn) - 1, fi%atoms%nat, fi%atoms%nat), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mt3_c")
         endif 

         if(.not. allocated(coul%mtir_c)) then
            allocate(coul%mtir_c(isize, isize), stat=info)
            if(info /= 0) call judft_error("Can't allocate coul%mtir_c")
         endif 
      endif
   end subroutine t_coul_alloc

   subroutine t_coul_init(coul)
      implicit none 
      class(t_coul), intent(inout) :: coul
      
      if(allocated(coul%mt1_r)) coul%mt1_r = 0
      if(allocated(coul%mt1_c)) coul%mt1_c = 0

      if(allocated(coul%mt2_r))   coul%mt2_r = 0
      if(allocated(coul%mt3_r))   coul%mt3_r = 0
      
      if(allocated(coul%mt2_c))   coul%mt2_c = 0
      if(allocated(coul%mt3_c))   coul%mt3_c = 0
   end subroutine t_coul_init

   subroutine allocate_hybdat(hybdat, fi, num_radfun_per_l)
      use m_types_fleurinput
      use m_judft
      implicit none
      class(t_hybdat), intent(inout) :: hybdat
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: num_radfun_per_l(:, :)
      integer                        :: ok(12)

      ok = -1
      allocate(hybdat%lmaxc(fi%atoms%ntype), &
               stat=ok(1), source=0)
      allocate(hybdat%bas1(fi%atoms%jmtd, maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(2), source=0.0)
      allocate(hybdat%bas2(fi%atoms%jmtd, maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(3), source=0.0)
      allocate(hybdat%bas1_MT(maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(4), source=0.0)
      allocate(hybdat%drbas1_MT(maxval(num_radfun_per_l), 0:fi%atoms%lmaxd, fi%atoms%ntype), &
               stat=ok(5), source=0.0)

      ! core allocs
      allocate(hybdat%nindxc(0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(6), source=0)
      allocate(hybdat%core1(fi%atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(7), source=0.0)
      allocate(hybdat%core2(fi%atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(8), source=0.0)
      allocate(hybdat%eig_c(hybdat%maxindxc, 0:hybdat%lmaxcd, fi%atoms%ntype), &
               stat=ok(9), source=0.0)

      allocate(hybdat%fac(0:hybdat%maxfac), stat=ok(10), source=0.0)
      allocate(hybdat%sfac(0:hybdat%maxfac), stat=ok(11), source=0.0)

      ALLOCATE(hybdat%gauntarr(2, 0:fi%atoms%lmaxd, 0:fi%atoms%lmaxd, 0:maxval(fi%hybinp%lcutm1), &
                               -fi%atoms%lmaxd:fi%atoms%lmaxd, -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1)), &
               stat=ok(12), source=0.0)

      if(any(ok /= 0)) then
         write(*, *) "allocation of hybdat failed. Error in array no.:"
         write(*, *) maxloc(abs(ok))
         call juDFT_error("allocation of hybdat failed. Error in array no is the outfile")
      endif
   end subroutine allocate_hybdat

   subroutine free_hybdat(hybdat)
      implicit none
      class(t_hybdat), intent(inout) :: hybdat

      if(allocated(hybdat%lmaxc)) deallocate(hybdat%lmaxc)
      if(allocated(hybdat%bas1)) deallocate(hybdat%bas1)
      if(allocated(hybdat%bas2)) deallocate(hybdat%bas2)
      if(allocated(hybdat%bas1_MT)) deallocate(hybdat%bas1_MT)
      if(allocated(hybdat%drbas1_MT)) deallocate(hybdat%drbas1_MT)
      if(allocated(hybdat%nindxc)) deallocate(hybdat%nindxc)
      if(allocated(hybdat%core1)) deallocate(hybdat%core1)
      if(allocated(hybdat%core2)) deallocate(hybdat%core2)
      if(allocated(hybdat%eig_c)) deallocate(hybdat%eig_c)
      if(allocated(hybdat%fac)) deallocate(hybdat%fac)
      if(allocated(hybdat%sfac)) deallocate(hybdat%sfac)
      if(allocated(hybdat%gauntarr)) deallocate(hybdat%gauntarr)
   end subroutine free_hybdat

   subroutine set_stepfunction(hybdat, cell, atoms, g, svol)
      use m_types_cell
      use m_types_atoms
      use m_judft
      implicit none
      class(t_hybdat), INTENT(INOUT) :: hybdat
      type(t_cell), INTENT(in)    :: cell
      type(t_atoms), INTENT(in)    :: atoms
      integer, INTENT(in)    :: g(3)
      real, INTENT(in)    :: svol
      integer :: i, j, k, ok

      if(.not. allocated(hybdat%stepfunc)) then
         call timestart("setup stepfunction")
         ALLOCATE(hybdat%stepfunc(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
         IF(ok /= 0) then
            call juDFT_error('wavefproducts_inv5: error allocation stepfunc')
         endif

         DO i = -g(1), g(1)
            DO j = -g(2), g(2)
               DO k = -g(3), g(3)
                  hybdat%stepfunc(i, j, k) = stepfunction(cell, atoms,[i, j, k])/svol
               END DO
            END DO
         END DO
         call timestop("setup stepfunction")
      endif

   end subroutine set_stepfunction

   !private subroutine
   FUNCTION stepfunction(cell, atoms, g)
      USE m_types_cell
      USE m_types_atoms
      USE m_constants
      IMPLICIT NONE

      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN) :: g(3)
      COMPLEX             :: stepfunction  !Is real in inversion case
      REAL                :: gnorm, gnorm3, r, fgr
      INTEGER             :: itype, ieq, icent

      gnorm = gptnorm(g, cell%bmat)
      gnorm3 = gnorm**3
      IF(abs(gnorm) < 1e-12) THEN
         stepfunction = 1
         DO itype = 1, atoms%ntype
            stepfunction = stepfunction - atoms%neq(itype)*atoms%volmts(itype)/cell%omtil
         END DO
      ELSE
         stepfunction = 0
         icent = 0
         DO itype = 1, atoms%ntype
            r = gnorm*atoms%rmt(itype)
            fgr = fpi_const*(sin(r) - r*cos(r))/gnorm3/cell%omtil
            DO ieq = 1, atoms%neq(itype)
               icent = icent + 1
               stepfunction = stepfunction - fgr*exp(-cmplx(0., tpi_const*dot_product(atoms%taual(:, icent), g)))
            ENDDO
         ENDDO
      ENDIF

   END FUNCTION stepfunction

   PURE FUNCTION gptnorm(gpt, bmat)
      IMPLICIT NONE
      REAL                :: gptnorm
      INTEGER, INTENT(IN)  :: gpt(3)
      REAL, INTENT(IN)     :: bmat(3, 3)

      gptnorm = norm2(matmul(gpt(:), bmat(:, :)))

   END FUNCTION gptnorm

   function mtir_size(fi, n_g, ikpt) result(isize)
      use m_types_fleurinput
      implicit none 
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: n_g(:), ikpt

      integer :: isize, itype, l

      isize = 0
      do itype = 1, fi%atoms%ntype
         do l = 0, fi%hybinp%lcutm1(itype)
            isize = isize + (2*l + 1)* fi%atoms%neq(itype)
         enddo 
      enddo

      isize = isize + n_g(ikpt)
   end function mtir_size
END MODULE m_types_hybdat
