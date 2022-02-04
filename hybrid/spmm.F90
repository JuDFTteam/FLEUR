module m_spmm
contains
   function calc_ibasm(fi, mpdata) result(ibasm)
      use m_types
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
