module m_work_package
   use m_types
   implicit none

   type t_band_package  
      integer :: start_idx, end_idx
   end type t_band_package

   type t_q_package 
      integer :: n_q 
      integer :: comm_q_group
   end type t_q_package 

   type t_k_package
      integer :: nk, rank, size
      type(t_q_package), allocatable :: q_packs(:)
   contains
      procedure :: init  => t_k_package_init 
      procedure :: print => t_k_package_print
   end type t_k_package 

   type t_work_package 
      integer :: rank, size
      type(t_k_package), allocatable :: k_packs(:)
   contains
      procedure :: init  => t_work_package_init 
      procedure :: print => t_work_package_print
      procedure :: owner_nk => t_work_package_owner_nk
      procedure :: has_nk => t_work_package_has_nk
   end type t_work_package





contains
   subroutine t_work_package_init(work_pack, fi, rank, size) 
      implicit none 
      class(t_work_package), intent(inout) :: work_pack
      type(t_fleurinput), intent(in)       :: fi
      integer, intent(in)                  :: rank, size

      work_pack%rank = rank
      work_pack%size = size
      call split_into_work_packages(work_pack, fi)

   end subroutine 

   subroutine t_k_package_init(k_pack, fi, nk)
      implicit none 
      class(t_k_package), intent(inout) :: k_pack
      type(t_fleurinput), intent(in)       :: fi
      integer, intent(in) :: nk
      integer             :: nkpt_eibz

      k_pack%nk = nk
      allocate(k_pack%q_packs(fi%kpts%nkpt_EIBZ(nk))) 
   end subroutine 

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

   subroutine split_into_work_packages(work_pack, fi)
      implicit none 
      class(t_work_package), intent(inout) :: work_pack
      type(t_fleurinput), intent(in)       :: fi
      integer :: my_num_ks, k_cnt, i

      if(fi%kpts%nkpt < work_pack%size) call judft_error("not enough k-points for work_packages")


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

         call work_pack%k_packs(k_cnt)%init(fi, i)
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
end module m_work_package
