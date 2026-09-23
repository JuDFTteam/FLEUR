module m_mtir_size
   implicit none
contains

   !> Size of the *reduced* MT+IR layout used by t_coul%mtir.
   !>
   !> Not the same as hybdat%nbasm: this counts (2l+1) per atom and l without the
   !> radial multiplicity num_radbasfn, because the MT block has already been
   !> contracted at that point.  Do not mix the two layouts or the accessors in
   !> t_hybdat (which describe the full mixed product basis) with each other.
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
