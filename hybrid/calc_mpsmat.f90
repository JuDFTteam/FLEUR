module m_calc_mpsmat
contains
   subroutine calc_mpsmat(fi, mpdata, smat)
      use m_types
      use m_judft
      use m_constants
      implicit none
      type(t_mat), intent(inout) :: smat
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata

      integer :: igpt2, igpt1, ic, itype, ineq, g(3)
      real    :: gnorm, rdum

      call timestart("calc smat")
      call smat%alloc(.False., mpdata%num_gpts(),mpdata%num_gpts())
      !$OMP PARALLEL DO default(none) schedule(guided)&
      !$OMP private(ic,igpt2,igpt1, g, gnorm, itype, rdum)&
      !$OMP shared(mpdata, fi, smat)
      DO igpt2 = 1, mpdata%num_gpts()
         DO igpt1 = 1, igpt2
            g = mpdata%g(:, igpt2) - mpdata%g(:, igpt1)
            gnorm = gptnorm(g, fi%cell%bmat)
            IF (abs(gnorm) < 1e-12) THEN
               DO itype = 1, fi%atoms%ntype
                  smat%data_c(igpt1, igpt2) = smat%data_c(igpt1, igpt2) + fi%atoms%neq(itype)*fpi_const*fi%atoms%rmt(itype)**3/3
               END DO
            ELSE
               ic = 0
               do ic =1,fi%atoms%nat
                  itype = fi%atoms%itype(ic)
                  rdum = fi%atoms%rmt(itype)*gnorm
                  rdum = fpi_const*(SIN(rdum) - rdum*COS(rdum))/gnorm**3
                  smat%data_c(igpt1, igpt2) = smat%data_c(igpt1, igpt2) &
                                       + rdum*EXP(ImagUnit*tpi_const*dot_PRODUCT(fi%atoms%taual(:, ic), g))
               END DO
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

      call smat%u2l()
      call timestop("calc smat")
   end subroutine calc_mpsmat
end module m_calc_mpsmat
