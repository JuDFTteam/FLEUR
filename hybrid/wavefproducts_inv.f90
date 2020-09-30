module m_wavefproducts_inv
   USE m_types_hybdat

   USE m_constants
   USE m_judft
   USE m_types
   USE m_types_hybinp
   USE m_util
   USE m_io_hybinp
   USE m_wrapper
   USE m_constants
   USE m_io_hybinp
   USE m_wavefproducts_aux

CONTAINS
   SUBROUTINE wavefproducts_inv(fi, ik, z_k, iq, jsp, bandoi, bandof, lapw, hybdat, mpdata, nococonv, stars, ikqpt, cmt_nk, cprod)
      IMPLICIT NONE
      type(t_fleurinput), intent(in):: fi
      TYPE(t_mpdata), intent(in)    :: mpdata
      type(t_nococonv), intent(in)  :: nococonv
      type(t_stars), intent(in)     :: stars
      type(t_mat), intent(in)       :: z_k  ! = z_k_p since ik < nkpt
      TYPE(t_lapw), INTENT(IN)      :: lapw
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      type(t_mat), intent(inout)    :: cprod

      ! - scalars -
      INTEGER, INTENT(IN)      :: jsp, ik, iq, bandoi, bandof
      INTEGER, INTENT(INOUT)   :: ikqpt
      complex, intent(in)  :: cmt_nk(:,:,:)


      ! - local scalars -
      INTEGER                 ::    g_t(3), psize
      REAL                    ::    kqpt(3), kqpthlp(3)

      type(t_mat) ::  z_kqpt_p
      complex, allocatable :: c_phase_kqpt(:)

      CALL timestart("wavefproducts_inv5")
      cprod%data_r = 0.0
      ikqpt = -1
      kqpthlp = fi%kpts%bkf(:, ik) + fi%kpts%bkf(:, iq)
      ! kqpt can lie outside the first BZ, transfer it back
      kqpt = fi%kpts%to_first_bz(kqpthlp)
      g_t = nint(kqpt - kqpthlp)

      ! determine number of kqpt
      ikqpt = fi%kpts%get_nk(kqpt)
      allocate (c_phase_kqpt(hybdat%nbands(fi%kpts%bkp(ikqpt))))
      IF (.not. fi%kpts%is_kpt(kqpt)) call juDFT_error('wavefproducts_inv5: k-point not found')

      call wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                 ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)

      call wavefproducts_inv_MT(fi, nococonv, jsp, bandoi, bandof, ik, iq, hybdat, mpdata, &
                                ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cmt_nk, cprod)

      CALL timestop("wavefproducts_inv5")

   END SUBROUTINE wavefproducts_inv

   subroutine wavefproducts_inv_IS(fi, jsp, lapw, ik, iq, bandoi, bandof, g_t, hybdat, mpdata, &
                                   nococonv, ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)

      implicit NONE
      type(t_fleurinput), intent(in):: fi
      TYPE(t_mpdata), intent(in)  :: mpdata
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      TYPE(t_lapw), INTENT(IN)      :: lapw
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      TYPE(t_mat), intent(in)       :: z_k ! = z_k_p since ik < nkpt
      TYPE(t_mat), intent(inout)    :: z_kqpt_p, cprod
      ! - scalars -
      INTEGER, INTENT(IN)      :: jsp, ik, iq, g_t(3), bandoi, bandof
      INTEGER, INTENT(IN)      :: ikqpt

      ! - arrays -
      complex, intent(inout)   :: c_phase_kqpt(hybdat%nbands(ikqpt))

      ! - local scalars -
      INTEGER                 ::    ic, ig, ig2, ig1, ok, igptm, iigptm, psize, b_idx
      INTEGER                 ::    ngpt0, n1, n2, nbasfcn, ierr, iband, iob
      REAL                    ::    rdum, rdum1
      TYPE(t_lapw)            ::    lapw_nkqpt
      type(t_mat)             ::    z_kqpt

      ! - local arrays -
      INTEGER, ALLOCATABLE    ::    pointer(:, :, :), gpt0(:, :)
      INTEGER                 ::    g(3)

      REAL                    ::    rarr1(bandoi:bandof)
      REAL, ALLOCATABLE       ::    z0(:, :), rtmp(:, :, :), rarr2(:,:)

      CALL timestart("wavefproducts_inv5 IR")
      allocate(rarr2(bandoi:bandof, hybdat%nbands(ik)), stat=ok, source=0.0)
      if(ok /= 0) call juDFT_error("Can't allocaten rarr2 in wavefproducts_inv_IS")
      psize = bandof-bandoi+1
      !
      ! compute G's fulfilling |bk(:,ikqpt) + G| <= rkmax
      !
      CALL lapw_nkqpt%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ikqpt, fi%cell, fi%sym%zrfs)
      nbasfcn = lapw_nkqpt%hyb_num_bas_fun(fi)
      call z_kqpt%alloc(.true., nbasfcn, fi%input%neig)
      call z_kqpt_p%init(z_kqpt)

      ! read in z at k-point ikqpt
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, ikqpt, jsp, z_kqpt, &
                  c_phase=c_phase_kqpt, parent_z=z_kqpt_p)

      g = maxval(abs(lapw%gvec(:, :lapw%nv(jsp), jsp)), dim=2) &
          + maxval(abs(lapw_nkqpt%gvec(:, :lapw_nkqpt%nv(jsp), jsp)), dim=2) &
          + maxval(abs(mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(iq), iq))), dim=2) + 1

      call hybdat%set_stepfunction(fi%cell, fi%atoms, g, sqrt(fi%cell%omtil))
      !
      ! convolute phi(n,k) with the step function and store in cpw0
      !

      !(1) prepare list of G vectors
      call prep_list_of_gvec(lapw, mpdata, g, g_t, iq, jsp, pointer, gpt0, ngpt0)

      !(2) calculate convolution
      call timestart("calc convolution")
      allocate (z0(bandoi:bandof, ngpt0), stat=ok, source=0.0)
      IF (ok /= 0) call juDFT_error('wavefproducts_inv5: error allocation z0')

      call timestart("step function")
      DO ig2 = 1, lapw_nkqpt%nv(jsp)
         rarr1 = z_kqpt%data_r(ig2, bandoi:bandof)
         DO ig = 1, ngpt0
            g = gpt0(:, ig) - lapw_nkqpt%gvec(:, ig2, jsp)
            rdum = REAL(hybdat%stepfunc(g(1), g(2), g(3)))
            DO n2 = bandoi,bandof
               z0(n2, ig) = z0(n2, ig) + rarr1(n2)*rdum
            END DO
         END DO
      END DO
      call timestop("step function")

      call timestart("hybrid g")
      allocate (rtmp(bandoi:bandof, hybdat%nbands(ik), mpdata%n_g(iq)), source=0.0)
      !$OMP PARALLEL DO default(none) &
      !$OMP private(igptm, ig1, iigptm, g, ig2, n1, n2) &
      !$OMP shared(mpdata, lapw, pointer, hybdat, rtmp, z0, z_k, g_t, jsp, iq, ik, bandoi, bandof) &
      !$OMP collapse(2)
      DO igptm = 1, mpdata%n_g(iq)
         DO ig1 = 1, lapw%nv(jsp)
            iigptm = mpdata%gptm_ptr(igptm, iq)
            g = lapw%gvec(:, ig1, jsp) + mpdata%g(:, iigptm) - g_t
            ig2 = pointer(g(1), g(2), g(3))

            IF (ig2 == 0) call juDFT_error('wavefproducts_inv5: pointer undefined')

            DO n1 = 1, hybdat%nbands(ik)
               DO n2 = bandoi,bandof
                  rtmp(n2, n1, igptm) = rtmp(n2, n1, igptm) + z_k%data_r(ig1, n1)*z0(n2, ig2)
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      call timestart("copy to cprod")
      do igptm = 1, mpdata%n_g(iq)
         ic = hybdat%nbasp + igptm
         do iob = 1,psize 
            b_idx = iob - 1 + bandoi
            do iband = 1, hybdat%nbands(ik)
               cprod%data_r(ic, iob + (iband-1)*psize) = rtmp(b_idx, iband, igptm)
            enddo 
         enddo
      enddo
      call timestop("copy to cprod")

      deallocate (rtmp)
      call timestop("hybrid g")
      call timestop("calc convolution")

      deallocate (z0, pointer, gpt0)
      CALL timestop("wavefproducts_inv5 IR")

   end subroutine wavefproducts_inv_IS

   subroutine wavefproducts_inv_MT(fi, nococonv, jsp, bandoi, bandof, ik, iq, hybdat, mpdata, &
                                   ikqpt, z_k_p, z_kqpt_p, c_phase_kqpt, ccmt_nk, cprod)
      use m_calc_cmt
      implicit NONE
      type(t_fleurinput), intent(in):: fi
      TYPE(t_mpdata), INTENT(IN)   :: mpdata
      type(t_nococonv), intent(in)  :: nococonv
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      type(t_mat), intent(in)       :: z_k_p
      type(t_mat), intent(in)       :: z_kqpt_p
      type(t_mat), intent(inout)    :: cprod

      ! - scalars -
      INTEGER, INTENT(IN)      :: ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(IN)      :: ikqpt

      ! - arrays -
      complex, intent(in)        :: c_phase_kqpt(hybdat%nbands(ikqpt)), ccmt_nk(:,:,:)

      ! - local scalars -
      INTEGER                 ::    i, iband, psize, iob
      INTEGER                 ::    iatom, iiatom, itype, ieq, ishift
      INTEGER                 ::    ibando, iatom1, iatom2, ioffset
      INTEGER                 ::    l, p, l1, m1, l2, m2, p1, p2, n, ok
      INTEGER                 ::    lm, lm1, lm2, lm_0, lm_00, lm1_0, lm2_0
      INTEGER                 ::    j, ll, m, lmp1, lmp2, lmp3, lmp4, lp1, lp2
      REAL                    ::    rdum, rfac, rfac1, rfac2, rdum1, rdum2
      REAL                    ::    add1, add2
      REAL                    ::    fac, fac1, fac2
      REAL                    ::    monepl1, monepl2, monepl, monepm1, monepm, moneplm, monepl1m1
      COMPLEX                 ::    cdum, cfac
      COMPLEX, PARAMETER      ::    img = (0.0, 1.0)
      LOGICAL                 ::    offdiag

      ! - local arrays -
      INTEGER                 ::    lmstart(0:fi%atoms%lmaxd, fi%atoms%ntype)

      REAL                    ::    cmt_nk(hybdat%nbands(ik), hybdat%maxlmindx, fi%atoms%nat)
      REAL, allocatable       ::    cmt_nkqpt(:,:,:)
      REAL, allocatable       ::    rarr2(:,:),rarr3(:,:,:)

      COMPLEX                 ::    cmplx_exp(fi%atoms%nat), cexp_nk(fi%atoms%nat)
      COMPLEX, ALLOCATABLE    ::    ccmt_nk2(:, :, :)
      COMPLEX, ALLOCATABLE    ::    ccmt_nkqpt(:, :, :)


      allocate(rarr2(bandoi:bandof, hybdat%nbands(ik)), stat=ok, source=0.0)
      if(ok /= 0) call juDFT_error("Can't alloc rarr2 in wavefproducts_inv_MT")
      allocate(rarr3(2,bandoi:bandof, hybdat%nbands(ik)), stat=ok, source=0.0)
      if(ok /= 0) call juDFT_error("Can't alloc rarr3 in wavefproducts_inv_MT")

      psize = bandof-bandoi+1
      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      DO itype = 1, fi%atoms%ntype
         DO l = 0, fi%atoms%lmax(itype)
            lmstart(l, itype) = sum([(mpdata%num_radfun_per_l(ll, itype)*(2*ll + 1), ll=0, l - 1)])
         END DO
      END DO

      ! read in cmt coefficient at k-point ik
      !ccmt_nkqpt(input%neig, hybdat%maxlmindx, fi%atoms%nat), &
      allocate (ccmt_nk2(hybdat%nbands(ik), hybdat%maxlmindx, fi%atoms%nat), &
                ccmt_nkqpt(psize, hybdat%maxlmindx, fi%atoms%nat), &
                source=cmplx(0.0, 0.0), stat=ok)
      IF (ok /= 0) call juDFT_error('wavefproducts_inv5: error allocation ccmt_nk/ccmt_nkqpt')
      allocate(cmt_nkqpt(bandoi:bandof, hybdat%maxlmindx, fi%atoms%nat), source=0.0, stat=ok)
      if (ok /= 0) call juDFT_error("error allocating cmt_nkqpt")

      !read in cmt coefficients at k+q point
      call calc_cmt(fi%atoms, fi%cell, fi%input, fi%noco, nococonv, fi%hybinp, hybdat, mpdata, fi%kpts, &
                    fi%sym, fi%oneD, z_kqpt_p, jsp, ikqpt, c_phase_kqpt, ccmt_nkqpt)

      iatom = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1

            cmplx_exp(iatom) = exp((-img)*tpi_const* &
                                   dot_product(fi%kpts%bkf(:, iq) + fi%kpts%bkf(:, ik), fi%atoms%taual(:, iatom)))

            cexp_nk(iatom) = exp((-img)*tpi_const* &
                                 dot_product(fi%kpts%bkf(:, ik), fi%atoms%taual(:, iatom)))
         END DO
      END DO

      rfac = 1./sqrt(2.0)
      cfac = -img/sqrt(2.0)
      iatom = 0
      DO itype = 1, fi%atoms%ntype
         DO ieq = 1, fi%atoms%neq(itype)
            iatom = iatom + 1
            ! determine the number of the inverse symmetric atom belonging to iatom
            IF (fi%sym%invsatnr(iatom) == 0) THEN
               iiatom = iatom
            ELSE
               iiatom = fi%sym%invsatnr(iatom)
            END IF
            ! the cmt coefficients at iatom and iiatom are made real in one step
            IF (iiatom < iatom) CYCLE
            lm1 = 0
            DO l = 0, fi%atoms%lmax(itype)
               DO m = -l, l
                  rdum = (-1)**(l + m)
                  DO p = 1, mpdata%num_radfun_per_l(l, itype)
                     lm1 = lm1 + 1
                     ! lm index at l,-m
                     lm2 = lm1 - 2*m*mpdata%num_radfun_per_l(l, itype)

                     IF (iatom == iiatom) THEN
                        IF (m < 0) THEN
                           cmt_nkqpt(:, lm1, iatom) = real((ccmt_nkqpt(:, lm1, iatom) + rdum*ccmt_nkqpt(:, lm2, iiatom))*cmplx_exp(iatom)*rfac)

                           cmt_nk(:, lm1, iatom) = real((ccmt_nk(:, lm1, iatom) + rdum*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*rfac)
                        ELSE IF (m > 0) THEN

                           cmt_nkqpt(:, lm1, iatom) = real((ccmt_nkqpt(:, lm1, iatom) - rdum*ccmt_nkqpt(:, lm2, iiatom))*cmplx_exp(iatom)*cfac)

                           cmt_nk(:, lm1, iatom) = real((ccmt_nk(:, lm1, iatom) - rdum*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*cfac)
                        ELSE
                           IF (mod(l, 2) == 0) THEN
                              cmt_nkqpt(:, lm1, iatom) = real(ccmt_nkqpt(:, lm1, iatom)*cmplx_exp(iatom))
                              cmt_nk(:, lm1, iatom) = real(ccmt_nk(:, lm1, iatom)*cexp_nk(iatom))
                           ELSE
                              cmt_nkqpt(:, lm1, iatom) = real(ccmt_nkqpt(:, lm1, iatom)*(-img)*cmplx_exp(iatom))
                              cmt_nk(:, lm1, iatom) = real(ccmt_nk(:, lm1, iatom)*(-img)*cexp_nk(iatom))
                           END IF
                        END IF
                     ELSE
                        cdum = rdum*cmplx_exp(iatom)*cmplx_exp(iiatom)
                        cmt_nkqpt(:, lm1, iatom) = real((ccmt_nkqpt(:, lm1, iatom) + cdum*ccmt_nkqpt(:, lm2, iiatom))*rfac)

                        cmt_nkqpt(:, lm1, iiatom) = real((ccmt_nkqpt(:, lm1, iatom) - cdum*ccmt_nkqpt(:, lm2, iiatom))*cfac)

                        cdum = rdum*cexp_nk(iatom)*cexp_nk(iiatom)
                        cmt_nk(:, lm1, iatom) = real((ccmt_nk(:, lm1, iatom) + cdum*ccmt_nk(:, lm2, iiatom))*rfac)

                        cmt_nk(:, lm1, iiatom) = real((ccmt_nk(:, lm1, iatom) - cdum*ccmt_nk(:, lm2, iiatom))*cfac)
                     END IF

                  END DO
               END DO
            END DO

         END DO
      END DO
      deallocate (ccmt_nkqpt)

      lm_0 = 0
      lm_00 = 0
      iatom1 = 0
      iiatom = 0

      DO itype = 1, fi%atoms%ntype
         ioffset = sum([((2*ll + 1)*mpdata%num_radbasfn(ll, itype), ll=0, fi%hybinp%lcutm1(itype))])
         lm_0 = lm_00
         DO ieq = 1, fi%atoms%neq(itype)
            iatom1 = iatom1 + 1
            IF (fi%sym%invsatnr(iatom1) == 0) THEN
               iatom2 = iatom1
            ELSE
               iatom2 = fi%sym%invsatnr(iatom1)
            END IF
            IF (iatom1 > iatom2) CYCLE

            IF (iatom1 /= iatom2) THEN
               call timestart("iatom1 neq iatom2")
               ! loop over l of mixed basis
               DO l = 0, fi%hybinp%lcutm1(itype)
                  ! loop over basis functions products, which belong to l
                  DO n = 1, hybdat%nindxp1(l, itype)

                     ! determine l1,p1 and l2,p2 for the basis functions, which can generate l
                     l1 = mpdata%l1(n, l, itype)
                     l2 = mpdata%l2(n, l, itype)
                     p1 = mpdata%n1(n, l, itype)
                     p2 = mpdata%n2(n, l, itype)

                     ! condition for Gaunt coefficients
                     IF (mod(l + l1 + l2, 2) /= 0) CYCLE

                     offdiag = l1 /= l2 .or. p1 /= p2 ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                     !(leading to the same basis-function product)

                     lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                     lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                     lm = lm_0
                     lp1 = lm1_0 + p1
                     lp2 = lm2_0 + p2

                     ! sum over different m of mixed basis functions with qn l
                     DO m = -l, l
                        rarr3 = 0.0

                        ! go to lm index for m1=-l1
                        lmp1 = lm1_0 + p1

                        DO m1 = -l1, l1
                           ! Gaunt condition -m1+m2-m=0
                           m2 = m1 + m
                           IF (abs(m2) <= l2) THEN
                              lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              ! precalculated Gaunt coefficient
                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                              IF (abs(rdum) > 1e-12) THEN
                                 !$OMP parallel do default(none) collapse(2) &
                                 !$OMP private(iband, ibando, rdum1, rdum2) &
                                 !$OMP shared(hybdat, bandoi, bandof, rdum, rarr3, cmt_nkqpt,cmt_nk) &
                                 !$OMP shared(iatom1, iatom2,lmp1,lmp2, ik)
                                 DO iband = 1, hybdat%nbands(ik)
                                    DO ibando = bandoi,bandof
                                       rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                       rdum2 = rdum*cmt_nk(iband, lmp1, iatom2)

                                       rarr3(1, ibando, iband) = rarr3(1, ibando, iband) &
                                                                 + rdum1*cmt_nkqpt(ibando, lmp2, iatom1) + rdum2*cmt_nkqpt(ibando, lmp2, iatom2)

                                       rarr3(2, ibando, iband) = rarr3(2, ibando, iband) &
                                                                 + rdum1*cmt_nkqpt(ibando, lmp2, iatom2) - rdum2*cmt_nkqpt(ibando, lmp2, iatom1)

                                    END DO  !ibando
                                 END DO  !iband
                                 !$OMP END parallel DO
                              END IF  ! rdum
                           END IF  ! abs(m2) .le. l2

                           m2 = m1 - m ! switch role of b1 and b2
                           IF (abs(m2) <= l2 .and. offdiag) THEN
                              lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              rdum = hybdat%gauntarr(2, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                              IF (abs(rdum) > 1e-12) THEN
                                 !$OMP parallel do default(none) collapse(2) &
                                 !$OMP private(iband, ibando, rdum1, rdum2) &
                                 !$OMP shared(hybdat, bandoi, bandof, rdum, rarr3, cmt_nkqpt,cmt_nk) &
                                 !$OMP shared(iatom1, iatom2,lmp1,lmp2, ik)
                                 DO iband = 1, hybdat%nbands(ik)
                                    DO ibando = bandoi,bandof
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                       rdum2 = rdum*cmt_nk(iband, lmp2, iatom2)
                                       rarr3(1, ibando, iband) = rarr3(1, ibando, iband) &
                                                                 + rdum1*cmt_nkqpt(ibando, lmp1, iatom1) + rdum2*cmt_nkqpt(ibando, lmp1, iatom2)

                                       rarr3(2, ibando, iband) = rarr3(2, ibando, iband) &
                                                                 + rdum1*cmt_nkqpt(ibando, lmp1, iatom2) - rdum2*cmt_nkqpt(ibando, lmp1, iatom1)
                                    END DO  !ibando
                                 END DO  !iband
                                 !$OMP END parallel DO
                              END IF  ! rdum .ne. 0
                           END IF  ! abs(m2) .le. l2 .and. offdiag

                           ! go to lmp start index for next m1-quantum number
                           lmp1 = lmp1 + mpdata%num_radfun_per_l(l1, itype)

                        END DO  !m1

                        ishift = -2*m*mpdata%num_radbasfn(l, itype)

                        ! go to lm mixed basis startindx for l and m
                        lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                        lm2 = lm + (iatom2 - 1 - iiatom)*ioffset + ishift

                        rdum = tpi_const*dot_product(fi%kpts%bkf(:, iq), fi%atoms%taual(:, iatom1))
                        rfac1 = sin(rdum)/sqrt(2.0)
                        rfac2 = cos(rdum)/sqrt(2.0)
                        !$OMP PARALLEL DO default(none) collapse(3) &
                        !$OMP private(iband, ibando, i, iob, rdum1, rdum2, add1, add2, j) &
                        !$OMP shared(cprod, hybdat, psize, lm1, lm2, l, n, itype, rarr3)&
                        !$OMP shared(bandoi,bandof,rfac1,rfac2, ik, mpdata)
                        DO iband = 1, hybdat%nbands(ik)
                           DO ibando = bandoi,bandof
                              DO i = 1, mpdata%num_radbasfn(l, itype)
                                 iob = ibando + 1 - bandoi
                                 rdum1 = rarr3(1, ibando, iband)
                                 rdum2 = rarr3(2, ibando, iband)
                                 add1 = rdum1*rfac2 + rdum2*rfac1
                                 add2 = rdum2*rfac2 - rdum1*rfac1
                                 j = lm1 + i
                                 cprod%data_r(j, iob + (iband-1)*psize) &
                                    = cprod%data_r(j, iob + (iband-1)*psize) + hybdat%prodm(i, n, l, itype)*add1
                                 j = lm2 + i
                                 cprod%data_r(j, iob + (iband-1)*psize) &
                                    = cprod%data_r(j, iob + (iband-1)*psize) + hybdat%prodm(i, n, l, itype)*add2

                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband
                        !$OMP end parallel do

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpdata%num_radbasfn(l, itype)

                     END DO  !m

                  END DO !n
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
                  IF (lm /= lm_0) call juDFT_error('wavefproducts_inv5: counting of lm-index incorrect (bug?)')
               END DO !l
               call timestop("iatom1 neq iatom2")
            ELSE !case: iatom1==iatom2
               call timestart("iatom1 eq iatom2")

               ! loop over l of mixed basis
               monepl = -1
               DO l = 0, fi%hybinp%lcutm1(itype)
                  monepl = -monepl
                  ! loop over basis functions products, which belong to l
                  DO n = 1, hybdat%nindxp1(l, itype)

                     ! determine l1,p1 and l2,p2 for the basis functions, which can generate l
                     l1 = mpdata%l1(n, l, itype)
                     l2 = mpdata%l2(n, l, itype)
                     p1 = mpdata%n1(n, l, itype)
                     p2 = mpdata%n2(n, l, itype)

                     ! condition for Gaunt coefficients
                     IF (mod(l + l1 + l2, 2) /= 0) CYCLE

                     offdiag = l1 /= l2 .or. p1 /= p2 ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                     !(leading to the same basis-function product)

                     lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                     lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                     lm = lm_0
                     lp1 = lm1_0 + p1
                     lp2 = lm2_0 + p2

                     ! calculate (-1)**l1 and (-1)**l2 (monepl = minus one power l)
                     monepl1 = (-1)**l1
                     monepl2 = (-1)**l2

                     ! sum over different m of mixed basis functions with qn l

                     !
                     !case m<0
                     !

                     monepm = -monepl
                     DO m = -l, -1
                        monepm = -monepm
                        moneplm = monepl*monepm

                        ! calculate the contributions which are identical for m>0 and m <0
                        rarr2 = 0.0
                        IF (abs(m) <= l2) THEN
                           lmp1 = lp1 + l1*mpdata%num_radfun_per_l(l1, itype)
                           IF (mod(l1, 2) == 0) THEN
                              lmp2 = lp2 + (m + l2)*mpdata%num_radfun_per_l(l2, itype)
                           ELSE
                              lmp2 = lp2 + (-m + l2)*mpdata%num_radfun_per_l(l2, itype)
                           END IF

                           rdum = hybdat%gauntarr(1, l1, l2, l, 0, m)
                           IF (abs(rdum) > 1e-12) THEN
                              DO iband = 1, hybdat%nbands(ik)
                                 rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                 IF (mod(l1, 2) /= 0) rdum1 = moneplm*rdum1
                                 DO ibando = bandoi,bandof
                                    rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, -m, m)
                              IF (abs(rdum) > 1e-12) THEN
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l1, 2) == 0) rdum1 = moneplm*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp1, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0
                           END IF  ! offdiag

                        END IF  ! abs(m) .le. l2

                        IF (abs(m) <= l1) THEN
                           IF (mod(l2, 2) == 0) THEN
                              lmp3 = lp1 + (m + l1)*mpdata%num_radfun_per_l(l1, itype)
                           ELSE
                              lmp3 = lp1 + (-m + l1)*mpdata%num_radfun_per_l(l1, itype)
                           END IF
                           lmp2 = lp2 + l2*mpdata%num_radfun_per_l(l2, itype)

                           rdum = hybdat%gauntarr(1, l1, l2, l, -m, m)
                           IF (abs(rdum) > 1e-12) THEN
                              DO iband = 1, hybdat%nbands(ik)
                                 rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                 IF (mod(l2, 2) == 0) rdum1 = moneplm*rdum1
                                 DO ibando = bandoi,bandof
                                    rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, 0, m)
                              IF (abs(rdum) > 1e-12) THEN
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l2, 2) /= 0) rdum1 = moneplm*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp3, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0
                           END IF  ! offdiag

                        END IF  ! abs(m) .le. l2

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !go to lm index for m1=-l1
                        lmp1 = lp1
                        monepm1 = -monepl1
                        DO m1 = -l1, l1
                           monepm1 = -monepm1
                           IF (m1 == 0) THEN
                              lmp1 = lmp1 + mpdata%num_radfun_per_l(l1, itype)
                              CYCLE
                           END IF
                           ! (-1)**(l1+m1)
                           monepl1m1 = monepl1*monepm1
                           m2 = m1 + m
                           IF (abs(m2) <= l2 .and. m2 /= 0) THEN
                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                              IF (abs(rdum) > 1e-12) THEN
                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (-m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                 END IF
                                 rdum = rdum/sqrt(2.0)
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(1, l2, l1, l, m2, -m)
                                 IF (abs(rdum) > 1e-12) THEN
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1
                                    ELSE
                                       lmp3 = lmp1 - 2*m1*mpdata%num_radfun_per_l(l1, itype)
                                       fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                    END IF
                                    rdum = moneplm*rdum/sqrt(2.0)
                                    DO iband = 1, hybdat%nbands(ik)
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                       IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                       DO ibando = bandoi,bandof
                                          rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband
                                 END IF  ! rdum .ne. 0
                              END IF  ! offdiag

                           END IF  ! abs(m2) .le. l2 .and. m2 .ne. 0

                           m2 = m1 - m
                           IF (abs(m2) <= l2 .and. m2 /= 0) THEN

                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, -m)
                              IF (abs(rdum) > 1e-12) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (-m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                 END IF
                                 rdum = moneplm*rdum/sqrt(2.0)
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!moneplm*rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband

                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(1, l2, l1, l, m2, m)
                                 IF (abs(rdum) > 1e-12) THEN
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    IF (sign(1, m1) + sign(1, m2) /= 0) THEN
                                       lmp3 = lmp1
                                    ELSE
                                       lmp3 = lmp1 - 2*m1*mpdata%num_radfun_per_l(l1, itype)
                                       fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                    END IF
                                    rdum = rdum/sqrt(2.0)
                                    DO iband = 1, hybdat%nbands(ik)
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                       IF (sign(1, m1) + sign(1, m2) == 0) rdum1 = fac*rdum1
                                       DO ibando = bandoi,bandof
                                          rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband
                                 END IF  ! rdum .ne. 0
                              END IF  ! offdiag

                           END IF  ! abs(m2) .le. l2 .and. m1 .ne. 0

                           !go to lmp start index for next m1-quantum number
                           lmp1 = lmp1 + mpdata%num_radfun_per_l(l1, itype)
                        END DO  ! m1

                        ! go to lm mixed basis startindx for l and m
                        lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                        DO iband = 1, hybdat%nbands(ik)
                           DO ibando = bandoi,bandof
                              iob = ibando + 1 - bandoi
                              rdum = rarr2(ibando, iband)
                              DO i = 1, mpdata%num_radbasfn(l, itype)
                                 j = lm1 + i
                                 cprod%data_r(j, iob + (iband-1)*psize) &
                                    = cprod%data_r(j, iob + (iband-1)*psize) + hybdat%prodm(i, n, l, itype)*rdum
                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpdata%num_radbasfn(l, itype)

                     END DO  ! m=-l,-1

                     !
                     !case m=0
                     !

                     m = 0
                     rarr2 = 0.0
                     lmp1 = lp1

                     monepm1 = -monepl1

                     DO m1 = -l1, l1
                        m2 = m1
                        monepm1 = -monepm1
                        IF (abs(m2) <= l2) THEN

                           IF (mod(l, 2) == 0) THEN
                              lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              !lmp3 and lmp4 are variables, which avoid an if clause in the loop
                              lmp3 = lmp2
                              lmp4 = lmp1
                           ELSE
                              lmp2 = lp2 + (-m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              !lmp3 and lmp3 are variables, which avoid an if clause in the loop
                              lmp3 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              lmp4 = lmp1 - 2*m1*mpdata%num_radfun_per_l(l1, itype)

                              fac1 = monepl1*monepm1 ! (-1)**(l1+m1)
                              fac2 = monepl2*monepm1 ! (-1)**(l2+m1)
                           END IF

                           !precalculated Gaunt coefficient
                           IF (mod(l, 2) == 0) THEN
                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                           ELSE
                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)*fac1
                           END IF
                           IF (abs(rdum) > 1e-12) THEN
                              DO iband = 1, hybdat%nbands(ik)
                                 rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                 DO ibando = bandoi,bandof
                                    rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum.ne.0

                           !change role of b1 and b2
                           IF (offdiag) THEN
                              IF (mod(l, 2) == 0) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, m)
                              ELSE
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, m)*fac2
                              END IF
                              IF (abs(rdum) > 1e-12) THEN
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp4, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum.ne.0
                           END IF  ! offdiag

                        END IF  ! abs(m2).le.l2

                        ! go to lmp start index for next m1-quantum number
                        lmp1 = lmp1 + mpdata%num_radfun_per_l(l1, itype)
                     END DO  !m1

                     ! go to lm mixed basis startindx for l and m
                     lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                     DO iband = 1, hybdat%nbands(ik)
                        DO ibando = bandoi,bandof
                           iob = ibando + 1 - bandoi
                           rdum = rarr2(ibando, iband)
                           DO i = 1, mpdata%num_radbasfn(l, itype)
                              j = lm1 + i
                              cprod%data_r(j, iob + (iband-1)*psize) &
                                 = cprod%data_r(j, iob + (iband-1)*psize) + hybdat%prodm(i, n, l, itype)*rdum
                           END DO  !i -> loop over mixed basis functions
                        END DO  !ibando
                     END DO  !iband

                     ! go to lm start index for next m-quantum number
                     lm = lm + mpdata%num_radbasfn(l, itype)

                     !
                     ! case: m>0
                     !

                     rarr2 = 0.0
                     monepm = 1
                     DO m = 1, l
                        monepm = -monepm
                        moneplm = monepl*monepm

                        ! calculate the contributions which are identical for m>0 and m <0
                        rarr2 = 0.0
                        IF (abs(m) <= l2) THEN
                           lmp1 = lp1 + l1*mpdata%num_radfun_per_l(l1, itype)
                           IF (mod(l1, 2) == 0) THEN
                              lmp2 = lp2 + (m + l2)*mpdata%num_radfun_per_l(l2, itype)
                           ELSE
                              lmp2 = lp2 + (-m + l2)*mpdata%num_radfun_per_l(l2, itype)
                           END IF

                           rdum = hybdat%gauntarr(1, l1, l2, l, 0, m)
                           IF (abs(rdum) > 1e-12) THEN
                              DO iband = 1, hybdat%nbands(ik)
                                 rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                 IF (mod(l1, 2) /= 0) rdum1 = moneplm*rdum1
                                 DO ibando = bandoi,bandof
                                    rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, -m, m)
                              IF (abs(rdum) > 1e-12) THEN
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l1, 2) == 0) rdum1 = moneplm*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp1, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0
                           END IF  ! offdiag

                        END IF  ! abs(m) .le. l2

                        IF (abs(m) <= l1) THEN
                           IF (mod(l2, 2) == 0) THEN
                              lmp3 = lp1 + (m + l1)*mpdata%num_radfun_per_l(l1, itype)
                           ELSE
                              lmp3 = lp1 + (-m + l1)*mpdata%num_radfun_per_l(l1, itype)
                           END IF
                           lmp2 = lp2 + l2*mpdata%num_radfun_per_l(l2, itype)

                           rdum = hybdat%gauntarr(1, l1, l2, l, -m, m)
                           IF (abs(rdum) > 1e-12) THEN
                              DO iband = 1, hybdat%nbands(ik)
                                 rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                 IF (mod(l2, 2) == 0) rdum1 = moneplm*rdum1
                                 DO ibando = bandoi,bandof
                                    rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, 0, m)
                              IF (abs(rdum) > 1e-12) THEN
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l2, 2) /= 0) rdum1 = moneplm*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp3, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0
                           END IF  ! offdiag

                        END IF  ! abs(m) .le. l2

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !go to lm index for m1=-l1
                        lmp1 = lp1
                        monepm1 = -monepl1
                        DO m1 = -l1, l1
                           monepm1 = -monepm1
                           IF (m1 == 0) THEN
                              lmp1 = lmp1 + mpdata%num_radfun_per_l(l1, itype)
                              CYCLE
                           END IF
                           m2 = m1 + m
                           ! (-1)**(l1+m1)
                           monepl1m1 = monepl1*monepm1
                           IF (abs(m2) <= l2 .and. m2 /= 0) THEN
                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                              IF (abs(rdum) > 1e-12) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (-m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    fac = -moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))/2
                                 END IF

                                 rdum = -moneplm*monepl1m1*rdum/sqrt(2.0)
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!-moneplm*monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband

                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, -m)
                                 IF (abs(rdum) > 1e-12) THEN
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1 - 2*m1*mpdata%num_radfun_per_l(l1, itype)
                                    ELSE
                                       lmp3 = lmp1
                                       fac = 1/2.*monepl1m1*(sign(1, m2) - sign(1, m1))
                                    END IF
                                    rdum = monepl1m1*moneplm*rdum/sqrt(2.0)
                                    DO iband = 1, hybdat%nbands(ik)
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!monepl1m1*moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                       IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                       DO ibando = bandoi,bandof
                                          rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband

                                 END IF  ! rdum

                              END IF  ! offdiag
                           END IF  ! abs(m2) .le. l2 .and. m2 .ne. 0

                           m2 = m1 - m
                           IF (abs(m2) <= l2 .and. m2 /= 0) THEN

                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, -m)
                              IF (abs(rdum) > 1e-12) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (-m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m1) - sign(1, m2))
                                 END IF
                                 rdum = monepl1m1*rdum/sqrt(2.0)
                                 DO iband = 1, hybdat%nbands(ik)
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = rdum1*fac
                                    DO ibando = bandoi,bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, m)
                                 IF (abs(rdum) > 1e-12) THEN
                                    lmp2 = lp2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1 - 2*m1*mpdata%num_radfun_per_l(l1, itype)
                                    ELSE
                                       lmp3 = lmp1
                                       fac = -monepl1m1*(sign(1, m1) - sign(1, m2))/2
                                    END IF
                                    rdum = -monepl1m1*rdum/sqrt(2.0)
                                    DO iband = 1, hybdat%nbands(ik)
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!-monepl1m1*rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                       IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                       DO ibando = bandoi,bandof
                                          rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt_nkqpt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband

                                 END IF  ! rdum .ne. 0
                              END IF  ! offdiag

                           END IF  !  abs(m2) .le. l2 .and. m2 .ne. 0

                           !go to lmp start index for next m1-quantum number
                           lmp1 = lmp1 + mpdata%num_radfun_per_l(l1, itype)
                        END DO  ! m1

                        ! multiply rarr2 by (-1)**(l+m+1)
                        rarr2(:, :) = (-1)*moneplm*rarr2(:, :)

                        ! go to lm mixed basis startindx for l and m
                        lm1 = lm + (iatom1 - 1 - iiatom)*ioffset

                        DO iband = 1, hybdat%nbands(ik)
                           DO ibando = bandoi,bandof
                              iob = ibando + 1 - bandoi
                              rdum = rarr2(ibando, iband)
                              DO i = 1, mpdata%num_radbasfn(l, itype)
                                 j = lm1 + i
                                 cprod%data_r(j, iob + (iband-1)*psize)&
                                    = cprod%data_r(j, iob + (iband-1)*psize) + hybdat%prodm(i, n, l, itype)*rdum
                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpdata%num_radbasfn(l, itype)

                     END DO  ! m=1,l

                  END DO !n
                  lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the m start index of the next l-quantum number
                  IF (lm /= lm_0) call juDFT_error('wavefproducts_inv5: counting of lm-index incorrect (bug?)')
               END DO !l

               call timestop("iatom1 eq iatom2")
            END IF  ! iatom1 .ne. iatom2

            lm_0 = lm_00
         END DO !ieq
         iiatom = iiatom + fi%atoms%neq(itype)
         lm_00 = lm_00 + fi%atoms%neq(itype)*ioffset
      END DO  !itype
   end subroutine wavefproducts_inv_MT
end module m_wavefproducts_inv
