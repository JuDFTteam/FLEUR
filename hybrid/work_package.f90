module m_work_package
   use m_types
   implicit none

   type t_band_package  
      integer :: start_idx, psize, rank, size
   contains 
      procedure :: init => t_band_package_init
   end type t_band_package

   type t_q_package 
      integer :: rank, size, ptr
      type(t_band_package), allocatable :: band_packs(:)
   contains
      procedure :: init => t_q_package_init
      procedure :: free => t_q_package_free
   end type t_q_package 

   type t_k_package
      integer :: nk, rank, size
      type(t_q_package), allocatable :: q_packs(:)
   contains
      procedure :: init  => t_k_package_init 
      procedure :: print => t_k_package_print
      procedure :: free  => t_k_package_free
   end type t_k_package 

   type t_work_package 
      integer :: rank, size
      type(t_k_package), allocatable :: k_packs(:)
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

   subroutine t_work_package_init(work_pack, fi, hybdat, jsp, rank, size) 
      implicit none 
      class(t_work_package), intent(inout) :: work_pack
      type(t_fleurinput), intent(in)       :: fi
      type(t_hybdat), intent(in)           :: hybdat 
      integer, intent(in)                  :: rank, size, jsp

      work_pack%rank = rank
      work_pack%size = size
      call split_into_work_packages(work_pack, fi, hybdat, jsp)

   end subroutine 

   subroutine t_k_package_init(k_pack, fi, hybdat, jsp, nk)
      implicit none 
      class(t_k_package), intent(inout) :: k_pack
      type(t_fleurinput), intent(in)    :: fi
      type(t_hybdat), intent(in)        :: hybdat
      integer, intent(in) :: nk, jsp
      integer             :: iq, jq

      k_pack%nk = nk
      allocate(k_pack%q_packs(fi%kpts%EIBZ(nk)%nkpt)) 
      do iq = 1,fi%kpts%EIBZ(nk)%nkpt
         jq = fi%kpts%EIBZ(nk)%pointer(iq)
         call k_pack%q_packs(iq)%init(fi, hybdat, jsp, nk, iq, jq)
      enddo
   end subroutine t_k_package_init

   subroutine t_q_package_init(q_pack, fi, hybdat, jsp, nk, rank, ptr)
      implicit none 
      class(t_q_package), intent(inout) :: q_pack 
      type(t_fleurinput), intent(in)    :: fi
      type(t_hybdat), intent(in)        :: hybdat
      integer, intent(in)               :: rank, ptr, jsp, nk

      real                 :: target_psize
      integer              :: n_parts, ikqpt, i
      integer, allocatable :: start_idx(:), psize(:)

      q_pack%rank  = rank 
      q_pack%size  = fi%kpts%EIBZ(nk)%nkpt
      q_pack%ptr = ptr

   ! arrays should be less than 5 gb
      if(fi%sym%invs) then
         target_psize = 5e9/( 8.0 * maxval(hybdat%nbasm) * MIN(fi%hybinp%bands1, fi%input%neig)) 
      else
         target_psize = 5e9/(16.0 * maxval(hybdat%nbasm) * MIN(fi%hybinp%bands1, fi%input%neig)) 
      endif

      ikqpt = fi%kpts%get_nk(fi%kpts%to_first_bz(fi%kpts%bkf(:,nk) + fi%kpts%bkf(:,ptr)))

      n_parts = ceiling(hybdat%nobd(ikqpt, jsp)/target_psize)
      allocate(start_idx(n_parts), psize(n_parts))
      allocate(q_pack%band_packs(n_parts))

      call split_band_loop(hybdat%nobd(ikqpt, jsp), n_parts, start_idx, psize)

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

   subroutine split_into_work_packages(work_pack, fi, hybdat, jsp)
      implicit none 
      class(t_work_package), intent(inout) :: work_pack
      type(t_fleurinput), intent(in)       :: fi
      type(t_hybdat), intent(in)           :: hybdat
      integer, intent(in)                  :: jsp
      integer :: my_num_ks, k_cnt, i 
      
      if(work_pack%rank < modulo(fi%kpts%nkpt, work_pack%size)) then
         my_num_ks = ceiling(1.0*fi%kpts%nkpt / work_pack%size)
      else 
         my_num_ks = floor(1.0*fi%kpts%nkpt / work_pack%size)
      endif

      allocate(work_pack%k_packs(my_num_ks))
      
      ! get my k-list
      k_cnt = 1
      do i = work_pack%rank+1 ,fi%kpts%nkpt ,work_pack%size
         work_pack%k_packs(k_cnt)%rank = k_cnt -1
         work_pack%k_packs(k_cnt)%size = my_num_ks

         call work_pack%k_packs(k_cnt)%init(fi, hybdat, jsp, i)
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
      do i = 1, work_pack%k_packs(1)%size 
         if (work_pack%k_packs(i)%nk == nk) then
            has_nk = .True.
            exit
         endif
      enddo
   end function t_work_package_has_nk

   subroutine split_band_loop(n_total, n_parts, start_idx, psize)
      use m_types
      implicit none
      integer, intent(in)                 :: n_total, n_parts
      integer, allocatable, intent(inout) :: start_idx(:), psize(:)

      integer             :: i, big_size, small_size, end_idx

      if(allocated(start_idx)) deallocate(start_idx)
      if(allocated(psize)) deallocate(psize)
      allocate(start_idx(n_parts), psize(n_parts))

      small_size = floor((1.0*n_total)/n_parts)
      big_size = small_size +1

      end_idx = 0
      do i = 1,n_parts
         psize(i) = merge(big_size, small_size,i <= mod(n_total, n_parts))

         start_idx(i) = end_idx + 1
         end_idx = start_idx(i) + psize(i) - 1
      enddo
   end subroutine split_band_loop
end module m_work_package
