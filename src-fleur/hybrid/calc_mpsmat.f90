module m_calc_mpsmat
   use m_types
   use m_judft
   use m_constants
contains
   subroutine calc_mpsmat(fi, mpdata, smat)
      implicit none
      type(t_mat), intent(inout) :: smat
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata

      integer :: igpt2, igpt1

      call timestart("calc smat")
      call smat%alloc(.False., mpdata%num_gpts(),mpdata%num_gpts())
      !$OMP PARALLEL DO default(shared) schedule(guided)&
      !$OMP private(igpt2,igpt1)&
      !$OMP shared(mpdata, fi, smat)
      DO igpt2 = 1, mpdata%num_gpts()
         DO igpt1 = 1, mpdata%num_gpts() !igpt2
            smat%data_c(igpt1, igpt2) = calc_smat_elem(fi, mpdata, igpt1, igpt2)
         END DO
      END DO
      !$OMP END PARALLEL DO

      ! call smat%u2l()
      call timestop("calc smat")
   end subroutine calc_mpsmat

   function calc_smat_elem(fi, mpdata, igpt1, igpt2) result(elem)
         ! Calculate the hermitian matrix smat(i,j) = sum(a) integral(MT(a)) exp[i(Gj-Gi)r] dr
      implicit none
      type(t_fleurinput), intent(in)  :: fi
      TYPE(t_mpdata), intent(in)      :: mpdata 
      integer, intent(in)             :: igpt1, igpt2
      complex :: elem

      integer :: g(3), ic, itype
      real    :: gnorm, rdum 

      elem = 0.0
      g = mpdata%g(:, igpt2) - mpdata%g(:, igpt1)
      gnorm = gptnorm(g, fi%cell%bmat)
      IF (abs(gnorm) < 1e-12) THEN
         DO itype = 1, fi%atoms%ntype
            elem = elem + fi%atoms%neq(itype)*fpi_const*fi%atoms%rmt(itype)**3/3
         END DO
      ELSE
         do ic = 1,fi%atoms%nat
            itype = fi%atoms%itype(ic)
            rdum = fi%atoms%rmt(itype)*gnorm
            rdum = fpi_const*(SIN(rdum) - rdum*COS(rdum))/gnorm**3
            elem = elem + rdum*EXP(ImagUnit*tpi_const*dot_PRODUCT(fi%atoms%taual(:, ic), g))
         END DO
      END IF
   end function calc_smat_elem
end module m_calc_mpsmat
