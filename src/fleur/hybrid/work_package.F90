module m_work_package
   use m_types
   use m_distribute_mpi
   use m_divide_most_evenly
   use m_mtir_size
#ifdef _OPENACC 
   use openacc
   use iso_c_binding
#endif
   implicit none
   private
   type,public:: t_band_package  
      integer :: start_idx, psize, rank, size
   contains 
      procedure :: init => t_band_package_init
   end type t_band_package

   type,public:: t_q_package 
      integer :: rank, size, ptr
      type(t_hybmpi) :: submpi
      type(t_band_package), allocatable :: band_packs(:)
   contains
      procedure :: init => t_q_package_init
      procedure :: free => t_q_package_free
   end type t_q_package 

   type,public:: t_qwps
      type(t_q_package), allocatable :: q_packs 
   end type t_qwps

   type,public:: t_k_package
      integer :: nk, rank, size
      type(t_hybmpi) :: submpi
      type(t_q_package), allocatable :: q_packs(:)
   contains
      procedure :: init  => t_k_package_init 
      procedure :: print => t_k_package_print
      procedure :: free  => t_k_package_free
   end type t_k_package 

   type,public:: t_work_package 
      integer :: rank, size, n_kpacks, max_kpacks
      type(t_k_package), allocatable :: k_packs(:)
      type(t_hybmpi) :: submpi
   contains
      procedure :: init  => t_work_package_init 
      procedure :: print => t_work_package_print
      procedure :: owner_nk => t_work_package_owner_nk
      procedure :: has_nk => t_work_package_has_nk
      procedure :: free => t_work_package_free
   end type t_work_package

contains
   subroutine t_work_package_free(work_pack)
      implicit none 
      class(t_work_package), intent(inout) :: work_pack 
      integer :: i

      if(allocated(work_pack%k_packs)) then 
         do i = 1, size(work_pack%k_packs)
            call work_pack%k_packs(i)%free() 
         enddo
         deallocate(work_pack%k_packs)
      endif
   end subroutine t_work_package_free 

   subroutine t_k_package_free(k_pack)
      implicit none 
      class(t_k_package), intent(inout) :: k_pack 
      integer :: i

      if(allocated(k_pack%q_packs)) then 
         do i = 1, size(k_pack%q_packs)
            call k_pack%q_packs(i)%free()
         enddo 
         deallocate(k_pack%q_packs)
      endif
   end subroutine t_k_package_free 

   subroutine t_q_package_free(q_pack) 
      implicit none 
      class(t_q_package), intent(inout) :: q_pack 

      if(allocated(q_pack%band_packs)) deallocate(q_pack%band_packs)
   end subroutine t_q_package_free

   subroutine t_work_package_init(work_pack, fi, hybdat, mpdata, wp_mpi, jsp, rank, size) 
      implicit none 
      class(t_work_package), intent(inout) :: work_pack
      type(t_fleurinput), intent(in)       :: fi
      type(t_hybdat), intent(in)           :: hybdat 
      type(t_mpdata), intent(in)           :: mpdata
      type(t_hybmpi), intent(in)           :: wp_mpi
      integer, intent(in)                  :: rank, size, jsp

      call timestart("t_work_package_init")
      work_pack%rank    = rank
      work_pack%size    = size
      work_pack%submpi  = wp_mpi

      call split_into_work_packages(work_pack, fi, hybdat, mpdata, jsp)

      call timestop("t_work_package_init")
   end subroutine t_work_package_init

   subroutine t_k_package_init(k_pack, fi, hybdat, mpdata, k_wide_mpi, jsp, nk)
      implicit none 
      class(t_k_package), intent(inout) :: k_pack
      type(t_fleurinput), intent(in)    :: fi
      type(t_hybdat), intent(in)        :: hybdat
      type(t_mpdata), intent(in)        :: mpdata
      type(t_hybmpi), intent(in)        :: k_wide_mpi
      type(t_hybmpi)                    :: q_wide_mpi

      integer, intent(in)  :: nk, jsp
      integer              :: iq, jq, loc_num_qs, i, cnt, n_groups, idx, q_rank, w_cnt
      integer, allocatable :: loc_qs(:)

      n_groups = min(k_wide_mpi%size, fi%kpts%EIBZ(nk)%nkpt)
      allocate(loc_qs(n_groups), source=0)
      do w_cnt = 1, n_groups 
         do i = w_cnt, fi%kpts%EIBZ(nk)%nkpt, n_groups 
            loc_qs(w_cnt) = loc_qs(w_cnt) + 1 
         enddo
      enddo

      call distribute_mpi(loc_qs, k_wide_mpi, q_wide_mpi, q_rank)

      k_pack%submpi = k_wide_mpi
      k_pack%nk = nk 
      
      allocate(k_pack%q_packs(loc_qs(q_rank+1)))
      cnt = 0
      do iq = q_rank+1,fi%kpts%EIBZ(nk)%nkpt, n_groups
         cnt = cnt + 1
         jq = fi%kpts%EIBZ(nk)%pointer(iq)
         call k_pack%q_packs(cnt)%init(fi, hybdat, mpdata, q_wide_mpi, jsp, nk, iq, jq)
      enddo
   end subroutine t_k_package_init

   subroutine t_q_package_init(q_pack, fi, hybdat, mpdata, q_wide_mpi, jsp, nk, rank, ptr)
      implicit none 
      class(t_q_package), intent(inout) :: q_pack 
      type(t_fleurinput), intent(in)    :: fi
      type(t_hybdat), intent(in)        :: hybdat
      type(t_mpdata), intent(in)        :: mpdata
      type(t_hybmpi), intent(in)        :: q_wide_mpi
      integer, intent(in)               :: rank, ptr, jsp, nk

      integer              :: target_psize
      integer              :: n_parts, ikqpt, i
      integer, allocatable :: start_idx(:), psize(:)

      q_pack%submpi = q_wide_mpi
      q_pack%rank   = rank 
      q_pack%size   = fi%kpts%EIBZ(nk)%nkpt
      q_pack%ptr    = ptr

      ikqpt = fi%kpts%get_nk(fi%kpts%to_first_bz(fi%kpts%bkf(:,nk) + fi%kpts%bkf(:,ptr)))
      n_parts = calc_n_parts(fi, hybdat, mpdata%n_g, q_pack, ikqpt, jsp)
      
      allocate(start_idx(n_parts), psize(n_parts))
      allocate(q_pack%band_packs(n_parts))

      call divide_most_evenly(hybdat%nobd(ikqpt, jsp), n_parts, start_idx, psize)

      do i = 1, n_parts
         call q_pack%band_packs(i)%init(start_idx(i), psize(i), i, n_parts)
      enddo
   end subroutine t_q_package_init

   subroutine t_band_package_init(band_pack, start_idx, psize, rank, size)
      implicit none 
      class(t_band_package), intent(inout) :: band_pack 
      integer, intent(in)                  :: rank, size, start_idx, psize

      band_pack%start_idx = start_idx
      band_pack%psize     = psize 
      band_pack%rank      = rank
      band_pack%size      = size
   end subroutine t_band_package_init

   subroutine t_work_package_print(work_pack)
      implicit none
      class(t_work_package), intent(inout) :: work_pack
      integer :: i 

      write (*,*) "WP (" // int2str(work_pack%rank) // "/" // int2str(work_pack%size) // ") has: "
      do i = 1,size(work_pack%k_packs)
         call work_pack%k_packs(i)%print()
      enddo
   end subroutine t_work_package_print 

   subroutine t_k_package_print(k_pack)
      implicit none 
      class(t_k_package), intent(in) :: k_pack 

      write (*,*) "kpoint: "
      write (*,*) "nk = ", k_pack%nk
   end subroutine t_k_package_print

   subroutine split_into_work_packages(work_pack, fi, hybdat, mpdata, jsp)
#ifdef CPP_MPI
      use mpi 
#endif
      implicit none 
      class(t_work_package), intent(inout) :: work_pack
      type(t_fleurinput), intent(in)       :: fi
      type(t_hybdat), intent(in)           :: hybdat
      type(t_mpdata), intent(in)           :: mpdata
      integer, intent(in)                  :: jsp
      integer :: k_cnt, i, ierr
      
      if(work_pack%rank < modulo(fi%kpts%nkpt, work_pack%size)) then
         work_pack%n_kpacks = ceiling(1.0*fi%kpts%nkpt / work_pack%size)
      else 
         work_pack%n_kpacks = floor(1.0*fi%kpts%nkpt / work_pack%size)
      endif
      allocate(work_pack%k_packs(work_pack%n_kpacks))

#ifdef CPP_MPI
      call MPI_AllReduce(work_pack%n_kpacks, work_pack%max_kpacks, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else    
      work_pack%max_kpacks = work_pack%n_kpacks
#endif
      if(work_pack%n_kpacks /= work_pack%max_kpacks) then
         call judft_warn("Your parallization is not efficient. Make sure that nkpts%pe == 0 or nkpts <= pe")
      endif 

      
      ! get my k-list
      k_cnt = 1
      do i = work_pack%rank+1, fi%kpts%nkpt, work_pack%size
         work_pack%k_packs(k_cnt)%rank = k_cnt -1
         work_pack%k_packs(k_cnt)%size = work_pack%n_kpacks

         call work_pack%k_packs(k_cnt)%init(fi, hybdat, mpdata, work_pack%submpi, jsp, i)
         k_cnt = k_cnt + 1
      enddo
   end subroutine split_into_work_packages


   function t_work_package_owner_nk(work_pack, nk) result(owner) 
      use m_types_hybmpi
      implicit none 
      class(t_work_package), intent(in) :: work_pack
      integer, intent(in)               :: nk
      integer                           :: owner

      owner = modulo(nk-1, work_pack%size)
   end function t_work_package_owner_nk

   function t_work_package_has_nk(work_pack, nk) result(has_nk) 
      implicit none 
      class(t_work_package), intent(in) :: work_pack
      integer, intent(in)               :: nk
      logical                           :: has_nk
      integer :: i 

      has_nk = .false.
      do i = 1, work_pack%n_kpacks 
         if (work_pack%k_packs(i)%nk == nk) then
            has_nk = .True.
            exit
         endif
      enddo
   end function t_work_package_has_nk

   function calc_n_parts(fi, hybdat, n_g, q_pack, ikqpt, jsp) result(n_parts)
      implicit none 
      type(t_fleurinput), intent(in) :: fi
      type(t_hybdat), intent(in)     :: hybdat
      integer, intent(in)            :: n_g(:), ikqpt, jsp
      class(t_q_package), intent(in) :: q_pack 
      
      integer :: n_parts, me, ierr, ikpt

      integer(8), parameter :: i8_one = 1
      integer(8)            :: coulomb_size, exch_size, indx_size, nsest_size, target_size, rc_factor 
      integer(8)            :: cprod_size, spmm_peak, max_peak
      integer(8)            :: max_nbasm, max_nbands, psize

      rc_factor  = merge(8, 16, fi%sym%invs)
      max_nbasm  = maxval(hybdat%nbasm)
      max_nbands = maxval(hybdat%nbands)

      target_size = target_memsize(fi, hybdat)
      coulomb_size = 0.0
      do ikpt = 1,fi%kpts%nkpt
         coulomb_size = max(int(mtir_size(fi, n_g, ikpt),kind=8)**2, coulomb_size)
      enddo
      ! size in byte
      coulomb_size = rc_factor * coulomb_size
      exch_size    = rc_factor * maxval(i8_one*hybdat%nbands)**2
      indx_size    = 4 *         maxval(i8_one*hybdat%nbands)**2
      nsest_size   = 4 *         maxval(i8_one*hybdat%nbands)

      psize = maxval(hybdat%nobd)
      do while(psize > 1)
         cprod_size = max_nbasm * max_nbands * psize * rc_factor

         spmm_peak = 2*cprod_size + coulomb_size + exch_size + indx_size + nsest_size

         max_peak = maxval([spmm_peak])

         if(max_peak <= target_size) then 
            exit 
         endif
         psize = psize - 1
      enddo

      n_parts = ceiling(1.0*maxval(hybdat%nobd)/psize)
      do while(mod(n_parts, q_pack%submpi%size) /= 0)
         n_parts = n_parts + 1
      enddo

      if(n_parts > hybdat%nobd(ikqpt, jsp)) then 
         write (*,*) "too many parts... reducing to nobd"
         n_parts = hybdat%nobd(ikqpt, jsp)
      endif
#ifdef CPP_MPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
#else
      me = 0
#endif
      !if(me == 0) write (*,*) "psize: " // int2str(psize) // " max_peak: " // int2str(max_peak) // " nparts: " // int2str(n_parts)
   end function calc_n_parts

   integer(8) function target_memsize(fi, hybdat)

      implicit none 
      type(t_fleurinput), intent(in) :: fi
      type(t_hybdat), intent(in)     :: hybdat

#ifdef _OPENACC    
      integer           :: ikpt
      integer(C_SIZE_T) :: gpu_mem
      real              :: coulomb_size, exch_size

      gpu_mem = acc_get_property(0,acc_device_current, acc_property_free_memory)
      target_memsize = int(0.75*gpu_mem, kind=8)
#else
      target_memsize = int(15e9, kind=8) ! 15 Gb
#endif
   end function target_memsize
end module m_work_package
