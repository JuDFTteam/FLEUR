module m_omp_checker
   
contains
   subroutine omp_checker()
      use omp_lib
      use m_judft
      use, intrinsic :: iso_c_binding
      implicit none

#ifdef __linux__ 
      interface 
         function findmycpu() bind(c)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int) :: findmycpu
         end function findmycpu
      end interface

      integer(kind=c_int), allocatable :: cpu(:)
      integer :: me, num_threads, mycpu, i


      !$omp parallel shared(cpu) private(me, num_threads, mycpu)
      me = omp_get_thread_num()
      num_threads = omp_get_num_threads()
      mycpu = findmycpu()

      if(me == 0) allocate(cpu(num_threads), source=-1)

      !$omp barrier

      cpu(me+1) = mycpu
      !$omp end parallel

      do i = 1,size(cpu)
         if(count(cpu(i) == cpu) /= 1) then
            call judft_warn("Multiple OMPs on one core: There are " // &
                           int2str(count(cpu(i) == cpu)) // " on cpu " // &
                           int2str(cpu))
            exit
         endif
      enddo
#endif
   end subroutine omp_checker
end module m_omp_checker
