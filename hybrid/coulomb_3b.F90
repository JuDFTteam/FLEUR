module m_coulomb_3b
contains

    subroutine coul_3b_igpt1_loop(fi, hybdat, mpdata, sphbesmoment, pqnrm, carr2, carr2a, carr2b, igpt2, ikpt, ix, coulmat)
        use m_types
        use m_calc_l_m_from_lm
        implicit none 
        type(t_fleurinput), intent(in)    :: fi
        TYPE(t_hybdat), INTENT(INOUT)     :: hybdat
        TYPE(t_mpdata), intent(in)        :: mpdata
        real, intent(in)                  :: sphbesmoment(:,:,:)
        
        integer, intent(in)    :: igpt2, ikpt, ix, pqnrm(:,:)
        complex, intent(in)    :: carr2(:,:), carr2a(:,:), carr2b(:,:)
        complex, intent(inout) :: coulmat(:,:)

        integer :: igpt1, iqnrm1, iy, ic, lm, l, m, itype
        complex :: csum, cdum

        call timestart("igpt1")
        DO igpt1 = 1, igpt2
           iy = hybdat%nbasp + igpt1
           !igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
           iqnrm1 = pqnrm(igpt1, ikpt)
           csum = 0
           !$OMP PARALLEL DO default(none) &
           !$OMP private(ic, itype, lm, l, m, cdum) &
           !$OMP shared(fi, carr2b, sphbesmoment, iqnrm1, igpt1, carr2, carr2a) &
           !$OMP reduction(+: csum) &
           !$OMP collapse(2)
           do ic = 1, fi%atoms%nat
              do lm = 1, (fi%hybinp%lexp+1)**2
                 itype = fi%atoms%itype(ic)
                 call calc_l_m_from_lm(lm, l, m)
                 cdum = carr2b(ic, igpt1)*sphbesmoment(l, itype, iqnrm1)
                 csum = csum + cdum*carr2(ic, lm)*CONJG(carr2a(lm, igpt1)) ! for coulomb
              END DO
           END DO
           !$OMP end parallel do
           coulmat(iy,ix) = coulmat(iy,ix) + csum/fi%cell%vol
        END DO
        call timestop("igpt1")

    end subroutine coul_3b_igpt1_loop

end module m_coulomb_3b