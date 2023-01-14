module m_omp_checker

contains
   subroutine omp_checker()
      !$ use omp_lib
      use m_judft
      USE m_constants
      use, intrinsic :: iso_c_binding
      implicit none
      
#ifdef CPP_SCHED
      interface
         function findmycpu() bind(c)
            use, intrinsic :: iso_c_binding
            integer(kind=c_int) :: findmycpu
         end function findmycpu
      end interface

      integer(kind=c_int), allocatable :: cpu(:)
      integer :: me, num_threads, mycpu, i
      logical :: l_problem


      !$omp parallel shared(cpu) private(me, num_threads, mycpu) default(none)
!$    if (.false.) then
      me = 0
      num_threads = 1
!$    endif
!$    me = omp_get_thread_num()
!$    num_threads = omp_get_num_threads()
      mycpu = findmycpu()

      if(me == 0) allocate(cpu(num_threads), source=-1)

      !$omp barrier

      cpu(me+1) = mycpu
      !$omp end parallel

      l_problem = .False.
      do i = 1,size(cpu)
         if(count(cpu(i) == cpu) /= 1) then
            WRITE(*,*) "The OpenMP parallelism seems to be weird"
            WRITE(*,*) "Multiple OMPs on one core: There are " // int2str(count(cpu(i) == cpu)) // &
                       " on cpu " // int2str(cpu(i))
            WRITE(oUnit,*) "The OpenMP parallelism seems to be weird"
            WRITE(oUnit,*) "Multiple OMPs on one core: There are " // int2str(count(cpu(i) == cpu)) // &
                           " on cpu " // int2str(cpu(i))

            l_problem = .True. 
            exit
         endif
      enddo

      if(l_problem) then
         !$omp parallel default(none) private(me, mycpu)
!$       me = omp_get_thread_num()  
         mycpu = findmycpu()
         write (*,*) "me: ", me, "my cpu:", mycpu 
         !$omp end parallel
      endif

#endif

   end subroutine omp_checker

   INTEGER FUNCTION getNumberOfThreads()
      !$ use omp_lib
      IMPLICIT NONE
      INTEGER :: numThreads
      numThreads = 1
      !$omp parallel
!$    numThreads = omp_get_num_threads()
      !$omp end parallel
      getNumberOfThreads = numThreads
   END FUNCTION getNumberOfThreads

end module m_omp_checker
