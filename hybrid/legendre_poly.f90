module m_legendre_poly
   use m_types
   use m_constants
contains
   subroutine legendre_poly(x, P)
      implicit none
      real, intent(in)                 :: x(:)
      type(t_mat), intent(inout)       :: P

      integer :: n_max, n

      n_max =  P%matsize2 - 1
      if(P%matsize1 /= size(x))  then 
         write (*,*) "P:", P%matsize1, P%matsize2 
         write (*,*) "x:", size(x) 
         write (*,*) "n_max+1:", n_max+1 

         call judft_error("Leg dimensions don't agree")
      endif 
 
      if (n_max >= 0) then
         if(P%l_real) then
            P%data_r(:, 1) = 1.0
         else 
            P%data_c(:, 1) = cmplx_1
         endif
      endif

      if (n_max >= 1) then
         if(P%l_real) then
            P%data_r(:, 2) = x
         else 
            P%data_c(:, 2) = x
         endif
      endif

      do n = 1, n_max - 1
         if(P%l_real) then
            P%data_r(:, n + 2) = ((2*n + 1)*x*P%data_r(:, n+1) - n*P%data_r(:, n))/(n + 1.0)
         else
            P%data_c(:, n + 2) = ((2*n + 1)*x*P%data_c(:, n+1) - n*P%data_c(:, n))/(n + 1.0)
         endif
      enddo
   end subroutine legendre_poly
end module m_legendre_poly