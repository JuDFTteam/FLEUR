module m_gamma_double_gpt_loop
   use m_types 
   use m_juDFT
   use m_constants
   use m_glob_tofrom_loc
contains
   subroutine gamma_double_gpt_loop(fi, fmpi, hybdat, mpdata, sphbesmoment, gmat, ngptm1, pgptm1, pqnrm, coul) 
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      type(t_hybdat), intent(in)        :: hybdat
      TYPE(t_mpdata), intent(in)        :: mpdata
      real, intent(in)                  :: sphbesmoment(:, :, :), gmat(:,:)
      integer, intent(in)               :: ngptm1(:), pgptm1(:,:), pqnrm(:,:)
      type(t_mat), intent(inout)        :: coul

      real        :: rdum
      integer     :: igpt0, igpt2, igptp2, iqnrm1, iqnrm2, igpt1, igptp1
      integer     :: ix, iy, ix_loc, iatm1, iatm2, itype1, itype2, pe_ix
      REAL        :: q1(3), q2(3), qnorm, qnorm1, qnorm2, rdum1
      complex     :: cdum

      rdum = (fpi_const)**(1.5)/fi%cell%vol**2*gmat(1, 1)
      call timestart("double gpt loop")
      DO igpt0 = 1, ngptm1(1)
         igpt2 = pgptm1(igpt0, 1)
         if(igpt2 /= 1) then
            ix = hybdat%nbasp + igpt2
            call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
            if(pe_ix == fmpi%n_rank) then
               iqnrm2 = pqnrm(igpt2, 1)
               igptp2 = mpdata%gptm_ptr(igpt2, 1)
               q2 = MATMUL(mpdata%g(:, igptp2), fi%cell%bmat)
               qnorm2 = norm2(q2)

               !$OMP PARALLEL DO default(none) schedule(dynamic) &
               !$OMP shared(igpt2, hybdat, fi, pqnrm, mpdata, q2, qnorm2, igptp2) &
               !$OMP shared(coul, ix_loc, rdum, sphbesmoment, iqnrm2)&
               !$OMP private(igpt1, iy, iqnrm1, igptp1, q1, qnorm1, rdum1, iatm1) &
               !$OMP private(itype1, iatm2, itype2, cdum)
               DO igpt1 = 2, igpt2
                  iy = hybdat%nbasp + igpt1
                  iqnrm1 = pqnrm(igpt1, 1)
                  igptp1 = mpdata%gptm_ptr(igpt1, 1)
                  q1 = MATMUL(mpdata%g(:, igptp1), fi%cell%bmat)
                  qnorm1 = norm2(q1)
                  rdum1 = dot_PRODUCT(q1, q2)/(qnorm1*qnorm2)
                  do iatm1 = 1,fi%atoms%nat
                     itype1 = fi%atoms%itype(iatm1)
                     do iatm2 = 1,fi%atoms%nat 
                        itype2 = fi%atoms%itype(iatm2)
                        cdum = EXP(CMPLX(0.0, 1.0)*tpi_const* &
                                 (-dot_PRODUCT(mpdata%g(:, igptp1), fi%atoms%taual(:, iatm1)) &
                                    + dot_PRODUCT(mpdata%g(:, igptp2), fi%atoms%taual(:, iatm2))))
                        coul%data_c(iy, ix_loc) = coul%data_c(iy, ix_loc) + rdum*cdum*( &
                                          -sphbesmoment(1, itype1, iqnrm1) &
                                          *sphbesmoment(1, itype2, iqnrm2)*rdum1/3 &
                                          - sphbesmoment(0, itype1, iqnrm1) &
                                          *sphbesmoment(2, itype2, iqnrm2)/6 &
                                          - sphbesmoment(2, itype1, iqnrm1) &
                                          *sphbesmoment(0, itype2, iqnrm2)/6 &
                                          + sphbesmoment(0, itype1, iqnrm1) &
                                          *sphbesmoment(1, itype2, iqnrm2)/qnorm2/2 &
                                          + sphbesmoment(1, itype1, iqnrm1) &
                                          *sphbesmoment(0, itype2, iqnrm2)/qnorm1/2)
                     END DO
                  END DO
               END DO
               !$OMP END PARALLEL DO
            endif !pe_ix
         endif
      END DO
      call timestop("double gpt loop")
   end subroutine gamma_double_gpt_loop 
end module m_gamma_double_gpt_loop