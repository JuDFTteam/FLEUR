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
      real, intent(in)                  :: sphbesmoment(0:, :, :), gmat(:,:)
      integer, intent(in)               :: ngptm1(:), pgptm1(:,:), pqnrm(:,:)
      complex, intent(inout)            :: coul(:,:)

      real        :: rdum
      integer     :: igpt0, igpt2, igptp2, iqnrm1, iqnrm2, igpt1, igptp1
      integer     :: ix, iy, ix_loc, iatm1, iatm2, itype1, itype2, pe_ix
      REAL        :: qnorm, rdum1
      complex     :: cdum
      
      real, allocatable :: qs(:,:), qnorms(:)

      call timestart("double gpt loop")

      call setup_q_and_qnorm(fi, mpdata, qs, qnorms)

      rdum = (fpi_const)**(1.5)/fi%cell%vol**2*gmat(1, 1)

      !$OMP PARALLEL DO default(none) schedule(dynamic) &
      !$OMP shared(fmpi, pgptm1, ngptm1, hybdat, fi, pqnrm, mpdata) &
      !$OMP shared(coul, rdum, sphbesmoment, qs, qnorms)&
      !$OMP private(igpt1, igpt2, iy, iqnrm1, igptp1, rdum1, iatm1, iqnrm2) &
      !$OMP private(itype1, iatm2, itype2, cdum, ix, pe_ix, ix_loc, igptp2)
      DO igpt0 = 1, ngptm1(1)
         igpt2 = pgptm1(igpt0, 1)
         if(igpt2 /= 1) then
            ix = hybdat%nbasp + igpt2
            pe_ix = mod((ix-1), fmpi%n_size)
            ix_loc = ((ix-1)/fmpi%n_size) +1
            if(pe_ix == fmpi%n_rank) then
               iqnrm2 = pqnrm(igpt2, 1)
               igptp2 = mpdata%gptm_ptr(igpt2, 1)
               DO igpt1 = 2, igpt2
                  iy = hybdat%nbasp + igpt1
                  iqnrm1 = pqnrm(igpt1, 1)
                  igptp1 = mpdata%gptm_ptr(igpt1, 1)
                  rdum1 = dot_PRODUCT(qs(:,igptp1), qs(:,igptp2))/(qnorms(igptp1)*qnorms(igptp2))
                  do iatm1 = 1,fi%atoms%nat
                     itype1 = fi%atoms%itype(iatm1)
                     do iatm2 = 1,fi%atoms%nat 
                        itype2 = fi%atoms%itype(iatm2)
                        cdum = EXP(CMPLX(0.0, 1.0)*tpi_const* &
                                 (-dot_PRODUCT(mpdata%g(:, igptp1), fi%atoms%taual(:, iatm1)) &
                                    + dot_PRODUCT(mpdata%g(:, igptp2), fi%atoms%taual(:, iatm2))))
                        coul(iy, ix_loc) = coul(iy, ix_loc) + rdum*cdum*( &
                                          -sphbesmoment(1, itype1, iqnrm1) &
                                          *sphbesmoment(1, itype2, iqnrm2)*rdum1/3 &
                                          - sphbesmoment(0, itype1, iqnrm1) &
                                          *sphbesmoment(2, itype2, iqnrm2)/6 &
                                          - sphbesmoment(2, itype1, iqnrm1) &
                                          *sphbesmoment(0, itype2, iqnrm2)/6 &
                                          + sphbesmoment(0, itype1, iqnrm1) &
                                          *sphbesmoment(1, itype2, iqnrm2)/qnorms(igptp2)/2 &
                                          + sphbesmoment(1, itype1, iqnrm1) &
                                          *sphbesmoment(0, itype2, iqnrm2)/qnorms(igptp1)/2)
                     END DO
                  END DO
               END DO
            endif !pe_ix
         endif
      END DO
      !$OMP end parallel do
      call timestop("double gpt loop")
   end subroutine gamma_double_gpt_loop 

   subroutine setup_q_and_qnorm(fi, mpdata, qs, qnorms)
      implicit none 
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata
      real, intent(inout), allocatable  :: qs(:,:), qnorms(:)

      real, allocatable :: g(:,:), tmp(:,:)
      integer :: i, num_gs

      num_gs = size(mpdata%g,2)
      allocate(qs(3,num_gs), source=0.0)

      g = real(mpdata%g)
      call dgemm("T", "N", 3,num_gs,3, 1.0,   fi%cell%bmat, 3,   g, 3, 0.0,  qs, 3)
      qnorms = norm2(qs, dim=1)
   end subroutine setup_q_and_qnorm
end module m_gamma_double_gpt_loop