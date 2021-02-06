module m_coulomb_3b
contains
    subroutine coulomb_3b(fi, fmpi, mpdata, hybdat, structconst, sphbesmoment, gmat, pqnrm, ngptm1, pgptm1, ikpt, coulmat)
        use m_types
        use m_constants
        use m_ylm
        use m_calc_l_m_from_lm
        implicit none
        type(t_fleurinput), intent(in)    :: fi
        TYPE(t_mpi), INTENT(IN)           :: fmpi
        TYPE(t_mpdata), intent(in)        :: mpdata
        TYPE(t_hybdat), INTENT(INOUT)     :: hybdat
        complex, intent(in)               :: structconst(:,:,:,:)
        real, intent(in)                  :: sphbesmoment(:, :, :), gmat(:, :)
        integer, intent(in)               :: ikpt, pqnrm(:, :), ngptm1(:), pgptm1(:, :)
        complex, intent(inout)            :: coulmat(:,:)

        integer :: igpt, igptp, iqnrm, iqnrm1, iqnrm2, ic, igpt0, igpt1, igpt2, ix, iy, igptp2
        integer :: iatom, itype, itype2, root, ierr
        integer :: lm, l, m, lm1, l1, m1, lm2, l2, m2, iat2
        REAL                 :: q(3)
        COMPLEX              :: y((fi%hybinp%lexp + 1)**2), cdum, cexp, csum
        COMPLEX, ALLOCATABLE :: carr2(:, :), carr2a(:, :), carr2b(:, :), structconst1(:, :)

                    ! group together quantities which depend only on l,m and igpt -> carr2a
        allocate (carr2a((fi%hybinp%lexp + 1)**2, maxval(mpdata%n_g)), carr2b(fi%atoms%nat, maxval(mpdata%n_g)))
        carr2a = 0; carr2b = 0
        DO igpt = 1, mpdata%n_g(ikpt)
           igptp = mpdata%gptm_ptr(igpt, ikpt)
           iqnrm = pqnrm(igpt, ikpt)
           q = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)

           call ylm4(fi%hybinp%lexp, q, y)

           y = CONJG(y)
           lm = 0
           DO l = 0, fi%hybinp%lexp
              DO M = -l, l
                 lm = lm + 1
                 carr2a(lm, igpt) = fpi_const*CMPLX(0.0, 1.0)**(l)*y(lm)
              END DO
           END DO
           DO ic = 1, fi%atoms%nat
              carr2b(ic, igpt) = EXP(-CMPLX(0.0, 1.0)*tpi_const* &
                                     dot_PRODUCT(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%atoms%taual(:, ic)))
           END DO
        END DO

        !finally we can loop over the plane waves (G: igpt1,igpt2)
        call timestart("loop over plane waves")
        allocate (carr2(fi%atoms%nat, (fi%hybinp%lexp + 1)**2), &
                  structconst1(fi%atoms%nat, (2*fi%hybinp%lexp + 1)**2), source=cmplx_0)
           
        DO igpt0 = 1+fmpi%n_rank, ngptm1(ikpt), fmpi%n_size !1,ngptm1(ikpt)
           igpt2 = pgptm1(igpt0, ikpt)
           ix = hybdat%nbasp + igpt2
           igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
           iqnrm2 = pqnrm(igpt2, ikpt)

           carr2 = 0

           call timestart("itype loops")
           do iatom = 1,fi%atoms%nat
              itype2 = fi%atoms%itype(iatom)
              cexp = CONJG(carr2b(iatom, igpt2))
              structconst1(:, :) = transpose(structconst(:, :, iatom, ikpt))
              
              !$OMP PARALLEL DO default(none) private(lm1,l1,m1,lm2,l2,m2,cdum,l,lm, iat2) &
              !$OMP shared(fi, sphbesmoment, itype2, iqnrm2, cexp, carr2a, igpt2, carr2, gmat, structconst1) 
              DO lm1 = 1, (fi%hybinp%lexp+1)**2
                 call calc_l_m_from_lm(lm1, l1, m1)
                 do lm2 = 1, (fi%hybinp%lexp+1)**2
                    call calc_l_m_from_lm(lm2, l2, m2)
                    cdum = (-1)**(l2 + m2)*sphbesmoment(l2, itype2, iqnrm2)*cexp*carr2a(lm2, igpt2)*gmat(lm1, lm2)
                    l = l1 + l2
                    lm = l**2 + l - l1 - m2 + (m1 + l1) + 1
                    do iat2 =1,fi%atoms%nat
                       carr2(iat2, lm1) = carr2(iat2,lm1) + cdum*structconst1(iat2, lm)
                    enddo
                 enddo
              enddo
              !$OMP end parallel do
           end do ! iatom

           call timestop("itype loops")

           call timestart("igpt1")
           iy = hybdat%nbasp
           DO igpt1 = 1, igpt2
              iy = iy + 1
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
        END DO !igpt0
        deallocate (carr2, carr2a, carr2b, structconst1)

#ifdef CPP_MPI
        call timestart("bcast itype&igpt1 loop")
        do igpt0 = 1, ngptm1(ikpt)
           root = mod(igpt0 - 1,fmpi%n_size)
           igpt2 = pgptm1(igpt0, ikpt)
           ix = hybdat%nbasp + igpt2
           call MPI_Bcast(coulmat(hybdat%nbasp+1,ix), igpt2, MPI_DOUBLE_COMPLEX, root, fmpi%sub_comm, ierr)
        enddo
        call timestop("bcast itype&igpt1 loop")
#endif
        call timestop("loop over plane waves")
    end subroutine coulomb_3b
end module m_coulomb_3b

