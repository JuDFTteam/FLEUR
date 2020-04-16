module m_wavefproducts_noinv
   USE m_types_hybdat

CONTAINS
   SUBROUTINE wavefproducts_noinv(fi, ik, z_k, iq, jsp, bandoi, bandof, lapw, hybdat, mpdata, nococonv, nkqpt, cprod)
      USE m_types
      use m_juDFT
      use m_constants, only: cmplx_0
      IMPLICIT NONE

      type(t_fleurinput), intent(in)  :: fi
      type(t_nococonv), intent(in)    :: nococonv
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)      :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat
      type(t_mat), intent(in)         :: z_k ! z_k is also z_k_p since ik < nkpt
      type(t_mat), intent(inout)      :: cprod

!     - scalars -
      INTEGER, INTENT(IN)        ::  ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(INOUT)     ::  nkqpt

      INTEGER              :: g_t(3), psize
      REAL                 :: kqpt(3), kqpthlp(3)
      complex              :: c_phase_k(hybdat%nbands(ik))
      complex, allocatable :: c_phase_kqpt(:)
      type(t_mat)          :: z_kqpt_p

      call timestart("wavefproducts_noinv")
      cprod%data_c = cmplx_0
      nkqpt = 0

      ! calculate nkpqt
      kqpthlp = fi%kpts%bkf(:, ik) + fi%kpts%bkf(:, iq)
      kqpt = fi%kpts%to_first_bz(kqpthlp)
      g_t = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      nkqpt = fi%kpts%get_nk(kqpt)
      allocate (c_phase_kqpt(hybdat%nbands(nkqpt)))

      IF (.not. fi%kpts%is_kpt(kqpt)) then
         call juDFT_error('wavefproducts: k-point not found')
      endif

      call wavefproducts_noinv_IS(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, nococonv, &
                                  nkqpt, z_k, c_phase_k, z_kqpt_p, c_phase_kqpt, cprod)

      call wavefproducts_noinv_MT(fi, ik, iq, bandoi, bandof, nococonv, mpdata, hybdat, &
                                  jsp, nkqpt, z_k, c_phase_k, z_kqpt_p, c_phase_kqpt, cprod)

      call timestop("wavefproducts_noinv")

   END SUBROUTINE wavefproducts_noinv

   subroutine wavefproducts_noinv_IS(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, nococonv, &
                                     nkqpt, z_k_p, c_phase_k, z_kqpt_p, c_phase_kqpt, cprod)
      use m_types
      use m_constants
      use m_wavefproducts_aux
      use m_judft
      use m_io_hybinp
      implicit NONE
      type(t_fleurinput), intent(in)  :: fi
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)      :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat
      type(t_mat), intent(in)         :: z_k_p
      type(t_mat), intent(inout)      :: z_kqpt_p, cprod

!     - scalars -
      INTEGER, INTENT(IN)      ::  ik, iq, jsp, g_t(3), bandoi, bandof
      INTEGER, INTENT(IN)      ::  nkqpt

!     - arrays -
      complex, intent(inout)    :: c_phase_k(hybdat%nbands(ik)), c_phase_kqpt(hybdat%nbands(nkqpt))

!     - local scalars -
      INTEGER                 :: ic, n1, n2, iob, iband
      INTEGER                 :: ig1, ig2, ig, psize
      INTEGER                 :: igptm, iigptm, ngpt0, nbasfcn

      COMPLEX                 ::  cdum, cdum1

      TYPE(t_lapw)            ::  lapw_nkqpt

!      - local arrays -
      INTEGER                 ::  g(3)
      INTEGER, ALLOCATABLE    ::  gpt0(:, :)
      INTEGER, ALLOCATABLE    ::  pointer(:, :, :)

      COMPLEX                 ::  carr1(bandoi:bandof)
      COMPLEX                 ::  carr(bandoi:bandof, hybdat%nbands(ik))
      TYPE(t_mat)             ::  z_kqpt
      COMPLEX, ALLOCATABLE    ::  z0(:, :), ctmp(:, :, :)

      call timestart("wavefproducts_noinv5 IR")
      !
      ! compute G's fulfilling |bk(:,nkqpt) + G| <= rkmax
      !
      CALL lapw_nkqpt%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nkqpt, fi%cell, fi%sym%zrfs)
      nbasfcn = lapw_nkqpt%hyb_num_bas_fun(fi)
      call z_kqpt%alloc(.false., nbasfcn, fi%input%neig)
      call z_kqpt_p%init(z_kqpt)

      ! read in z at k-point ik and nkqpt
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, nkqpt, jsp, z_kqpt, &
                  c_phase=c_phase_kqpt, parent_z=z_kqpt_p)

      g = maxval(abs(lapw%gvec(:, :lapw%nv(jsp), jsp)), dim=2) &
          + maxval(abs(lapw_nkqpt%gvec(:, :lapw_nkqpt%nv(jsp), jsp)), dim=2) &
          + maxval(abs(mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(iq), iq))), dim=2) + 1

      psize = bandof-bandoi+1

      call hybdat%set_stepfunction(fi%cell, fi%atoms, g, sqrt(fi%cell%omtil))

      !
      ! convolute phi(n,k) with the step function and store in cpw0
      !

      !(1) prepare list of G vectors
      call prep_list_of_gvec(lapw, mpdata, g, g_t, iq, jsp, pointer, gpt0, ngpt0)

      !(2) calculate convolution
      call timestart("calc convolution")
      call timestart("step function")
      ALLOCATE (z0(bandoi:bandof, ngpt0), source=cmplx_0)

      DO ig2 = 1, lapw_nkqpt%nv(jsp)
         if (z_kqpt%l_real) then
            carr1 = z_kqpt%data_r(ig2, bandoi:bandof)
         else
            carr1 = z_kqpt%data_c(ig2, bandoi:bandof)
         endif
         DO ig = 1, ngpt0
            g = gpt0(:, ig) - lapw_nkqpt%gvec(:, ig2, jsp)
            cdum = hybdat%stepfunc(g(1), g(2), g(3))
            DO n2 = bandoi,bandof
               z0(n2, ig) = z0(n2, ig) + carr1(n2)*cdum
            END DO
         END DO
      END DO
      call timestop("step function")

      call timestart("hybrid g")
      allocate (ctmp(bandoi:bandof, hybdat%nbands(ik), mpdata%n_g(iq)), source=(0.0, 0.0))
      if (z_k_p%l_real) then
         !$OMP PARALLEL DO default(none) &
         !$OMP private(igptm, ig1, iigptm, g, ig2, n1, n2) &
         !$OMP shared(mpdata, lapw, pointer, hybdat, ctmp, z0, z_k_p, g_t, jsp, iq, ik, bandoi, bandof) &
         !$OMP collapse(2)
         DO igptm = 1, mpdata%n_g(iq)
            DO ig1 = 1, lapw%nv(jsp)
               iigptm = mpdata%gptm_ptr(igptm, iq)
               g = lapw%gvec(:, ig1, jsp) + mpdata%g(:, iigptm) - g_t
               ig2 = pointer(g(1), g(2), g(3))
               IF (ig2 == 0) call juDFT_error('wavefproducts_noinv2: pointer undefined')

               DO n1 = 1, hybdat%nbands(ik)
                  DO n2 = bandoi,bandof
                     ctmp(n2, n1, igptm) = ctmp(n2, n1, igptm) + z_k_p%data_r(ig1, n1)*z0(n2, ig2)
                  END DO
               END DO

            END DO
         END DO
         !$OMP END PARALLEL DO
      else
         !$OMP PARALLEL DO default(none) &
         !$OMP private(igptm, ig1, iigptm, g, ig2, n1, n2) &
         !$OMP shared(mpdata, lapw, pointer, hybdat, ctmp, z0, z_k_p, g_t, jsp, iq, ik, bandoi, bandof) &
         !$OMP collapse(2)
         DO igptm = 1, mpdata%n_g(iq)
            DO ig1 = 1, lapw%nv(jsp)
               iigptm = mpdata%gptm_ptr(igptm, iq)
               g = lapw%gvec(:, ig1, jsp) + mpdata%g(:, iigptm) - g_t
               ig2 = pointer(g(1), g(2), g(3))
               IF (ig2 == 0) call juDFT_error('wavefproducts_noinv2: pointer undefined')

               DO n1 = 1, hybdat%nbands(ik)
                  DO n2 = bandoi,bandof
                     ctmp(n2, n1, igptm) = ctmp(n2, n1, igptm) + conjg(z_k_p%data_c(ig1, n1))*z0(n2, ig2)
                  END DO
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
      endif

      call timestart("copy to cprod")
      do igptm = 1, mpdata%n_g(iq)
         ic = hybdat%nbasp + igptm
         do iob = 1,psize 
            do iband = 1, hybdat%nbands(ik)
               cprod%data_c(ic, iob + (iband-1)*psize) = ctmp(iob, iband, igptm)
            enddo 
         enddo
      enddo
      call timestop("copy to cprod")

      call timestop("hybrid g")
      deallocate (z0, pointer, gpt0)
      call timestop("calc convolution")

      call timestop("wavefproducts_noinv5 IR")
   end subroutine wavefproducts_noinv_IS

   subroutine wavefproducts_noinv_MT(fi, ik, iq, bandoi, bandof, nococonv, mpdata, hybdat, jsp, ikqpt, &
                                     z_k_p, c_phase_k, z_kqpt_p, c_phase_kqpt, cprod)
      use m_types
      USE m_constants
      use m_io_hybinp
      use m_judft
      use m_wavefproducts_aux
      use m_calc_cmt
      IMPLICIT NONE
      type(t_fleurinput), intent(in)  :: fi
      type(t_nococonv), intent(in)    :: nococonv
      TYPE(t_mpdata), INTENT(IN)      :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat
      type(t_mat), intent(in)         :: z_k_p, z_kqpt_p
      type(t_mat), intent(inout)      :: cprod

      !     - scalars -
      INTEGER, INTENT(IN)     ::  ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(IN)     ::  ikqpt

      !     - arrays -
      complex, intent(in)     :: c_phase_k(hybdat%nbands(ik))
      complex, intent(in)     :: c_phase_kqpt(hybdat%nbands(ikqpt))

      !     - local scalars -
      INTEGER                 ::  ic, l, n, l1, l2, n1, n2, lm_0, lm1_0, lm2_0
      INTEGER                 ::  lm, lm1, lm2, m1, m2, i, ll, j, k
      INTEGER                 ::  itype, ieq, ic1, m, iob, iband, psize

      COMPLEX                 ::  atom_phase

      LOGICAL                 ::  offdiag

      !      - local arrays -
      INTEGER                 ::  lmstart(0:fi%atoms%lmaxd, fi%atoms%ntype)

      COMPLEX                 ::  carr(bandoi:bandof, hybdat%nbands(ik))
      COMPLEX                 ::  cmt_ikqpt(hybdat%nbands(ikqpt), hybdat%maxlmindx, fi%atoms%nat)
      COMPLEX                 ::  cmt_nk(hybdat%nbands(ik), hybdat%maxlmindx, fi%atoms%nat)

      call timestart("wavefproducts_noinv5 MT")
      psize = bandof-bandoi+1
      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      call timestart("set lmstart")
      DO itype = 1, fi%atoms%ntype
         DO l = 0, fi%atoms%lmax(itype)
            lmstart(l, itype) = sum([(mpdata%num_radfun_per_l(ll, itype)*(2*ll + 1), ll=0, l - 1)])
         END DO
      END DO
      call timestop("set lmstart")

      ! read in cmt coefficients from direct access file cmt
      call calc_cmt(fi%atoms, fi%cell, fi%input, fi%noco, nococonv, fi%hybinp, hybdat, mpdata, fi%kpts, &
                    fi%sym, fi%oneD, z_k_p, jsp, ik, c_phase_k, cmt_nk)
      call calc_cmt(fi%atoms, fi%cell, fi%input, fi%noco, nococonv, fi%hybinp, hybdat, mpdata, fi%kpts, &
                    fi%sym, fi%oneD, z_kqpt_p, jsp, ikqpt, c_phase_kqpt, cmt_ikqpt)

      call timestart("loop over l, l1, l2, n, n1, n2")
      !$OMP PARALLEL PRIVATE(m, carr, lm1, m1, m2, lm2, i,j,k, &
      !$OMP lm, n1, l1, n2, l2, offdiag, lm1_0, lm2_0, itype, ieq, &
      !$OMP ic, lm_0)
      lm_0 = 0
      ic = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            ic = ic + 1
            ic1 = 0

            atom_phase = exp(-ImagUnit*tpi_const*dot_product(fi%kpts%bkf(:, iq), fi%atoms%taual(:, ic)))

            DO l = 0, fi%hybinp%lcutm1(itype)

               DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products
                  call mpdata%set_nl(n, l, itype, n1, l1, n2, l2)

                  IF (mod(l1 + l2 + l, 2) == 0) THEN
                     offdiag = (l1 /= l2) .or. (n1 /= n2) ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                     !(leading to the same basis-function product)

                     lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                     lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                     lm = lm_0
                     !$OMP DO
                     DO m = -l, l
                        carr = 0.0

                        lm1 = lm1_0 + n1 ! go to lm index for m1=-l1
                        DO m1 = -l1, l1
                           m2 = m1 + m ! Gaunt condition -m1+m2-m=0
                           IF (abs(m2) <= l2) THEN
                              lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              IF (abs(hybdat%gauntarr(1, l1, l2, l, m1, m)) > 1e-12) THEN
                                 carr = carr + hybdat%gauntarr(1, l1, l2, l, m1, m) &
                                        *outer_prod(cmt_ikqpt(bandoi:bandof, lm2, ic), &
                                                    conjg(cmt_nk(1:hybdat%nbands(ik), lm1, ic)))
                              END IF
                           END IF

                           m2 = m1 - m ! switch role of b1 and b2
                           IF (abs(m2) <= l2 .and. offdiag) THEN
                              lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              IF (abs(hybdat%gauntarr(2, l1, l2, l, m1, m)) > 1e-12) THEN
                                 carr = carr + hybdat%gauntarr(2, l1, l2, l, m1, m) &
                                        *outer_prod(cmt_ikqpt(bandoi:bandof, lm1, ic), &
                                                    conjg(cmt_nk(1:hybdat%nbands(ik), lm2, ic)))
                              END IF
                           END IF

                           lm1 = lm1 + mpdata%num_radfun_per_l(l1, itype) ! go to lm start index for next m1-quantum number

                        END DO  !m1

                        lm = lm_0 + (m + l)*mpdata%num_radbasfn(l, itype)
                        do k = 1, hybdat%nbands(ik)
                           do j = 1, psize
                              DO i = 1, mpdata%num_radbasfn(l, itype)
                                 cprod%data_c(i + lm, j + (k-1)*psize) &
                                      = cprod%data_c(i + lm, j + (k-1)*psize) &
                                          + hybdat%prodm(i, n, l, itype)*carr(j+bandoi-1, k)*atom_phase
                              ENDDO
                           end do
                        end do
                     END DO
                     !$OMP END  DO
                  ENDIF
               END DO
               lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
            END DO
         END DO
      END DO
      !$OMP END PARALLEL
      call timestop("loop over l, l1, l2, n, n1, n2")
      call timestop("wavefproducts_noinv5 MT")
   end subroutine wavefproducts_noinv_MT
end module m_wavefproducts_noinv
