MODULE m_judft_para

CONTAINS
   subroutine juDFT_check_para()
      implicit none
      call check_omp_para()
   end subroutine juDFT_check_para

   subroutine check_omp_para()
      use omp_lib
      use m_judft_string
      use m_judft_stop
      implicit none
      real(8)  :: summe, t_omp, t_seq
      integer(4)  :: rank, size, ierr
      integer  :: i, omp_threads
      integer, parameter :: loop_end = 300000000

      summe = 0.0
      !$omp parallel reduction(+: t_omp)
         omp_threads = OMP_GET_NUM_THREADS()
         t_omp = OMP_GET_WTIME()
         !$omp do reduction(+:summe)
            do i = 1, loop_end*omp_threads
               summe = summe + 1.0
            enddo
         !$omp end do
         t_omp = OMP_GET_WTIME() - t_omp
      !$omp end parallel

      t_omp = t_omp / omp_threads

      summe = summe / omp_threads

      t_seq = OMP_GET_WTIME()
      do i = 1, loop_end
         summe = summe - 1.0
      enddo
      t_seq = OMP_GET_WTIME() - t_seq

      if( abs(t_seq/t_omp -1.0) < 0.1)then
         write (*,*) "Parallelization OK"
      else
         write (*,*) "t_omp = ", t_omp
         write (*,*) "t_seq = ", t_seq
         write (*,*) "Summe = ", summe

         call juDFT_warn("OMP parallelization underperform with a parallel efficiency of " // &
                         float2str(t_seq/t_omp), hint="check if your slurm files is set properly")
      endif
   end subroutine check_omp_para

END MODULE m_judft_para
