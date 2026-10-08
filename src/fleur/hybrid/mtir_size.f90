!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_mtir_size
   USE m_vac_const, ONLY: NVAC_MPB
   USE m_types_fleurinput
   implicit none
   PRIVATE
   PUBLIC :: mtir_size
contains

   !> Size of the *reduced* MT+IR layout used by t_coul%mtir.
   !>
   !> Not the same as hybdat%nbasm: this counts (2l+1) per atom and l without the
   !> radial multiplicity num_radbasfn, because the MT block has already been
   !> contracted at that point.  Do not mix the two layouts or the accessors in
   !> t_hybdat (which describe the full mixed product basis) with each other.
   function mtir_size(fi, n_g, ikpt, n_g_vac) result(isize)
      implicit none
      type(t_fleurinput), intent(in) :: fi
      integer, intent(in)            :: n_g(:), ikpt
      integer, intent(in), optional  :: n_g_vac(:)

      integer :: isize, itype, l

      isize = 0
      do itype = 1, fi%atoms%ntype
         do l = 0, fi%hybinp%lcutm1(itype)
            isize = isize + (2*l + 1)*fi%atoms%neq(itype)
         enddo
      enddo

      isize = isize + n_g(ikpt)

      ! films: vacuum carriers per (vacuum, G||) and the two g = 0 second carriers
      if (present(n_g_vac)) then
         if (fi%input%film) then
            isize = isize + n_g_vac(ikpt)*NVAC_MPB + NVAC_MPB
         end if
      endif
   end function mtir_size
end module m_mtir_size
