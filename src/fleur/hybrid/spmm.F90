!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_spmm
   use m_types_fleurinput
   use m_types_mpdata
   use m_judft
   implicit none
   private
   public :: calc_ibasm
contains
   function calc_ibasm(fi, mpdata) result(ibasm)
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      integer :: ibasm, iatom, itype, ieq, l, m

      call timestart("calc_ibasm")
      ibasm = 0
      iatom = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, fi%hybinp%lcutm1(itype)
               DO m = -l, l
                  ibasm = ibasm + mpdata%num_radbasfn(l, itype) - 1
               END DO
            END DO
         END DO
      END DO
      call timestop("calc_ibasm")
   end function calc_ibasm
end module m_spmm
