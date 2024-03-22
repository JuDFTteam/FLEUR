module m_divide_most_evenly 
   use m_judft
contains
   subroutine divide_most_evenly(n_total, n_parts, start_idx, psize)
      implicit none
      integer, intent(in)                 :: n_total, n_parts
      integer, allocatable, intent(inout) :: start_idx(:), psize(:)

      integer             :: i, big_size, small_size, end_idx

      if(allocated(start_idx)) deallocate(start_idx)
      if(allocated(psize)) deallocate(psize)
      allocate(start_idx(n_parts), psize(n_parts))

      if(n_parts == 0) call judft_error("You need more than 0 parts")
      if(n_parts > n_total) call judft_error("You can't have more n_parts, than n_total")

      small_size = floor((1.0*n_total)/n_parts)
      big_size = small_size +1

      end_idx = 0
      do i = 1,n_parts
         psize(i) = merge(big_size, small_size,i <= mod(n_total, n_parts))
         if(psize(i) == 0) then
            write (*,*) "n_total, n_parts", n_total, n_parts
            call judft_warn("some band_packs have 0 bands")
         endif
         start_idx(i) = end_idx + 1
         end_idx = start_idx(i) + psize(i) - 1
      enddo
   end subroutine divide_most_evenly
end module m_divide_most_evenly