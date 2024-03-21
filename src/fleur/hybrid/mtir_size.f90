module m_mtir_size
   implicit none
contains

   function mtir_size(fi, n_g, ikpt) result(isize)
      use m_types_fleurinput
      implicit none
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: n_g(:), ikpt

      integer :: isize, itype, l

      isize = 0
      do itype = 1, fi%atoms%ntype
         do l = 0, fi%hybinp%lcutm1(itype)
            isize = isize + (2*l + 1)*fi%atoms%neq(itype)
         enddo
      enddo

      isize = isize + n_g(ikpt)
   end function mtir_size
end module m_mtir_size
