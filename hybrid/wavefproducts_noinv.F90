module m_wavefproducts_noinv
   USE m_types_hybdat

CONTAINS
   SUBROUTINE wavefproducts_noinv(fi, ik, z_k, iq, jsp, bandoi, bandof, lapw, hybdat, mpdata, nococonv, stars, ikqpt, cmt_nk, cprod)
      USE m_types
      use m_juDFT
      use m_wavefproducts_aux
      use m_constants, only: cmplx_0
      IMPLICIT NONE

      type(t_fleurinput), intent(in)  :: fi
      type(t_nococonv), intent(in)    :: nococonv
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)      :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat
      type(t_mat), intent(in)         :: z_k ! z_k is also z_k_p since ik < nkpt
      type(t_stars), intent(in)       :: stars
      type(t_mat), intent(inout)      :: cprod

!     - scalars -
      INTEGER, INTENT(IN)        ::  ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(INOUT)     ::  ikqpt

      complex, intent(in)  :: cmt_nk(:,:,:)

      INTEGER              :: g_t(3)
      REAL                 :: kqpt(3), kqpthlp(3)
      complex, allocatable :: c_phase_kqpt(:)
      type(t_mat)          :: z_kqpt_p, cprod_tmp

      call timestart("wavefproducts_noinv")
      cprod%data_c = cmplx_0
      ikqpt = 0

      ! calculate ikqpt
      kqpthlp = fi%kpts%bkf(:, ik) + fi%kpts%bkf(:, iq)
      kqpt = fi%kpts%to_first_bz(kqpthlp)

      ! if k+q outside of first BZ put we need this shift
      g_t = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      ikqpt = fi%kpts%get_nk(kqpt)
      allocate (c_phase_kqpt(hybdat%nbands(ikqpt,jsp)))
      call cprod_tmp%init(cprod)

      IF (.not. fi%kpts%is_kpt(kqpt)) then
         call juDFT_error('wavefproducts: k-point not found')
      endif

      call wavefproducts_IS_FFT(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, stars, nococonv, &
                                  ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)

      ! call wavefproducts_noinv_IS(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, nococonv, &
      !                             ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)


      call wavefproducts_noinv_MT(fi, ik, iq, bandoi, bandof, nococonv, mpdata, hybdat, &
                                  jsp, ikqpt, z_kqpt_p, c_phase_kqpt, cmt_nk, cprod)

      call timestop("wavefproducts_noinv")

   END SUBROUTINE wavefproducts_noinv

   subroutine wavefproducts_noinv_IS(fi, ik, iq, g_t, jsp, bandoi, bandof, mpdata, hybdat, lapw, nococonv, &
                                     ikqpt, z_k, z_kqpt_p, c_phase_kqpt, cprod)
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
      type(t_mat), intent(in)         :: z_k
      type(t_mat), intent(inout)      :: z_kqpt_p, cprod

!     - scalars -
      INTEGER, INTENT(IN)      ::  ik, iq, jsp, g_t(3), bandoi, bandof
      INTEGER, INTENT(IN)      ::  ikqpt

!     - arrays -
      complex, intent(inout)    :: c_phase_kqpt(hybdat%nbands(ikqpt,jsp))

!     - local scalars -
      INTEGER                 :: ic, n1, n2, iob, iband, ok
      INTEGER                 :: ig1, ig2, ig, psize, b_idx
      INTEGER                 :: igptm, iigptm, ngpt0, nbasfcn

      COMPLEX                 ::  cdum

      TYPE(t_lapw)            ::  lapw_ikqpt

!      - local arrays -
      INTEGER                 ::  g(3)
      INTEGER, ALLOCATABLE    ::  gpt0(:, :)
      INTEGER, ALLOCATABLE    ::  pointer(:, :, :)

      COMPLEX                 ::  carr1(bandoi:bandof)
      TYPE(t_mat)             ::  z_kqpt
      COMPLEX, ALLOCATABLE    ::  z0(:, :), ctmp(:, :, :), carr(:,:)

      call timestart("wavefproducts_noinv5 IR")
      allocate(carr(bandoi:bandof, hybdat%nbands(ik,jsp)), stat=ok, source=cmplx_0)
      if(ok /= 0) call juDFT_error("Can't alloc carr in wavefproducts_noinv_IS")
      !
      ! compute G's fulfilling |bk(:,ikqpt) + G| <= rkmax
      !
      CALL lapw_ikqpt%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ikqpt, fi%cell, fi%sym%zrfs)
      nbasfcn = lapw_ikqpt%hyb_num_bas_fun(fi)
      call z_kqpt%alloc(.false., nbasfcn, fi%input%neig)
      call z_kqpt_p%init(z_kqpt)

      ! read in z at k-point ik and ikqpt
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv, fi%input, ikqpt, jsp, z_kqpt, &
                  c_phase=c_phase_kqpt, parent_z=z_kqpt_p)

      g = maxval(abs(lapw%gvec(:, :lapw%nv(jsp), jsp)), dim=2) &
          + maxval(abs(lapw_ikqpt%gvec(:, :lapw_ikqpt%nv(jsp), jsp)), dim=2) &
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

      DO ig2 = 1, lapw_ikqpt%nv(jsp)
         if (z_kqpt%l_real) then
            carr1 = z_kqpt%data_r(ig2, bandoi:bandof)
         else
            carr1 = z_kqpt%data_c(ig2, bandoi:bandof)
         endif
         DO ig = 1, ngpt0
            g = gpt0(:, ig) - lapw_ikqpt%gvec(:, ig2, jsp)
            cdum = hybdat%stepfunc(g(1), g(2), g(3))
            DO n2 = bandoi,bandof
               z0(n2, ig) = z0(n2, ig) + carr1(n2)*cdum
            END DO
         END DO
      END DO
      call timestop("step function")

      call timestart("hybrid g")
      allocate (ctmp(bandoi:bandof, hybdat%nbands(ik,jsp), mpdata%n_g(iq)), source=(0.0, 0.0))
      if (z_k%l_real) then
         !$OMP PARALLEL DO default(none) &
         !$OMP private(igptm, ig1, iigptm, g, ig2, n1, n2) &
         !$OMP shared(mpdata, lapw, pointer, hybdat, ctmp, z0, z_k, g_t, jsp, iq, ik, bandoi, bandof) &
         !$OMP collapse(2)
         DO igptm = 1, mpdata%n_g(iq)
            DO ig1 = 1, lapw%nv(jsp)
               iigptm = mpdata%gptm_ptr(igptm, iq)
               g = lapw%gvec(:, ig1, jsp) + mpdata%g(:, iigptm) - g_t
               ig2 = pointer(g(1), g(2), g(3))
               IF (ig2 == 0) call juDFT_error('wavefproducts_noinv2: pointer undefined')

               DO n1 = 1, hybdat%nbands(ik,jsp)
                  DO n2 = bandoi,bandof
                     ctmp(n2, n1, igptm) = ctmp(n2, n1, igptm) + z_k%data_r(ig1, n1)*z0(n2, ig2)
                  END DO
               END DO

            END DO
         END DO
         !$OMP END PARALLEL DO
      else
         !$OMP PARALLEL DO default(none) &
         !$OMP private(igptm, ig1, iigptm, g, ig2, n1, n2) &
         !$OMP shared(mpdata, lapw, pointer, hybdat, ctmp, z0, z_k, g_t, jsp, iq, ik, bandoi, bandof) &
         !$OMP collapse(2)
         DO igptm = 1, mpdata%n_g(iq)
            DO ig1 = 1, lapw%nv(jsp)
               iigptm = mpdata%gptm_ptr(igptm, iq)
               g = lapw%gvec(:, ig1, jsp) + mpdata%g(:, iigptm) - g_t
               ig2 = pointer(g(1), g(2), g(3))
               IF (ig2 == 0) call juDFT_error('wavefproducts_noinv2: pointer undefined')

               DO n1 = 1, hybdat%nbands(ik,jsp)
                  DO n2 = bandoi,bandof
                     ctmp(n2, n1, igptm) = ctmp(n2, n1, igptm) + conjg(z_k%data_c(ig1, n1))*z0(n2, ig2)
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
            b_idx = iob - 1 + bandoi
            do iband = 1, hybdat%nbands(ik,jsp)
               cprod%data_c(ic, iob + (iband-1)*psize) = ctmp(b_idx, iband, igptm)
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
                                     z_kqpt_p, c_phase_kqpt, cmt_nk, cprod)
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
      type(t_mat), intent(in)         :: z_kqpt_p
      type(t_mat), intent(inout)      :: cprod

      !     - scalars -
      INTEGER, INTENT(IN)     ::  ik, iq, jsp, bandoi, bandof
      INTEGER, INTENT(IN)     ::  ikqpt

      !     - arrays -
      complex, intent(in)     :: c_phase_kqpt(hybdat%nbands(ikqpt,jsp))

      complex, intent(in)    :: cmt_nk(:,:,:)

      !     - local scalars -
      INTEGER                 ::  iatm, l, n, l1, l2, n1, n2, lm_0, lm1_0, lm2_0
      INTEGER                 ::  lm, lm1, lm2, m1, m2, i, ll, j, k, ok
      INTEGER                 ::  itype, ieq, m, psize

      COMPLEX                 ::  atom_phase, cscal

      LOGICAL                 ::  offdiag

      !      - local arrays -
      INTEGER                 ::  lmstart(0:fi%atoms%lmaxd, fi%atoms%ntype)

      COMPLEX, allocatable    ::  cmt_ikqpt(:,:,:)

      call timestart("wavefproducts_noinv5 MT")
      allocate(cmt_ikqpt(bandoi:bandof, hybdat%maxlmindx, fi%atoms%nat), stat=ok, source=cmplx_0)
      if(ok /= 0) call juDFT_error("alloc cmt_ikqpt")
      
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
                    fi%sym, fi%oneD, z_kqpt_p, jsp, ikqpt, c_phase_kqpt, cmt_ikqpt)


      ! $acc enter data copyin(cprod, cprod%data_c,  hybdat, hybdat%nbands, hybdat%nindxp1, hybdat%gauntarr, hybdat%prodm) &
      ! $acc copyin(bandoi, bandof, lmstart, lm_0, mpdata, mpdata%num_radfun_per_l, mpdata%l1, mpdata%l2, mpdata%n1, mpdata%n2)&
      ! $acc copyin(cmt_ikqpt, cmt_nk, psize, ik)

      call timestart("loop over l, l1, l2, n, n1, n2")
      lm_0 = 0
      do iatm = 1,fi%atoms%nat 
         itype = fi%atoms%itype(iatm)
         atom_phase = exp(-ImagUnit*tpi_const*dot_product(fi%kpts%bkf(:, iq), fi%atoms%taual(:, iatm)))

         ! The default(shared) in the OMP part of the following loop is needed to avoid compilation issues on gfortran 7.5.
         DO l = 0, fi%hybinp%lcutm1(itype)
! #ifdef _OPENACC
!             !$acc parallel loop default(none) copyin(itype, l, iatm, atom_phase) collapse(2)&
!             !$acc private(k,j,n, n1, l1, n2, l2, offdiag, lm1_0, lm2_0, lm, m, cscal, lm1, m1, m2, lm2, i) &
!             !$acc present(cprod, cprod%data_c, hybdat, hybdat%nbands, hybdat%nindxp1, hybdat%gauntarr, hybdat%prodm)&
!             !$acc present(bandoi, bandof, lmstart, lm_0, mpdata, mpdata%num_radfun_per_l, mpdata%l1, mpdata%l2, mpdata%n1, mpdata%n2)&
!             !$acc present(cmt_ikqpt, cmt_nk, psize, ik)
! #else
            !$OMP PARALLEL DO default(shared) collapse(2) schedule(dynamic) & 
            !$OMP private(k,j,n, n1, l1, n2, l2, offdiag, lm1_0, lm2_0, lm, m, cscal, lm1, m1, m2, lm2, i)&
            !$OMP shared(hybdat, bandoi, bandof, lmstart, lm_0, mpdata, cmt_ikqpt, cmt_nk, cprod, itype, l) &
            !$OMP shared(iatm, psize, atom_phase, ik)
! #endif
            do k = 1, hybdat%nbands(ik,jsp)
               do j = bandoi, bandof 
                  DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products
                     ! don't call object funcktions in acc
                     l1 = mpdata%l1(n, l, itype) !
                     l2 = mpdata%l2(n, l, itype) ! current basis-function mpdatauct
                     n1 = mpdata%n1(n, l, itype) ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
                     n2 = mpdata%n2(n, l, itype) !
                     ! call mpdata%set_nl(n, l, itype, n1, l1, n2, l2)

                     IF (mod(l1 + l2 + l, 2) == 0) THEN
                        offdiag = (l1 /= l2) .or. (n1 /= n2) ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                        !(leading to the same basis-function product)

                        lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                        lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                        lm = lm_0
                        DO m = -l, l
                           cscal = 0.0

                           lm1 = lm1_0 + n1 ! go to lm index for m1=-l1
                           DO m1 = -l1, l1
                              m2 = m1 + m ! Gaunt condition -m1+m2-m=0
                              
                              IF (abs(m2) <= l2) THEN
                                 lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                 IF (abs(hybdat%gauntarr(1, l1, l2, l, m1, m)) > 1e-12) THEN
                                    cscal = cscal + hybdat%gauntarr(1, l1, l2, l, m1, m) &
                                                            * cmt_ikqpt(j, lm2, iatm) &
                                                               * conjg(cmt_nk(k, lm1, iatm))
                                 END IF
                              END IF

                              m2 = m1 - m ! switch role of b1 and b2
                              IF (abs(m2) <= l2 .and. offdiag) THEN
                                 lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                                 IF (abs(hybdat%gauntarr(2, l1, l2, l, m1, m)) > 1e-12) THEN
                                    cscal = cscal + hybdat%gauntarr(2, l1, l2, l, m1, m) &
                                                            * cmt_ikqpt(j, lm1, iatm) & 
                                                               * conjg(cmt_nk(k, lm2, iatm))
                                 END IF
                              END IF

                              lm1 = lm1 + mpdata%num_radfun_per_l(l1, itype) ! go to lm start index for next m1-quantum number

                           END DO  !m1

                           lm = lm_0 + (m + l)*mpdata%num_radbasfn(l, itype)
                           DO i = 1, mpdata%num_radbasfn(l, itype)
                              cprod%data_c(i + lm, (j-bandoi+1) + (k-1)*psize) &
                                 = cprod%data_c(i + lm, (j-bandoi+1) + (k-1)*psize) &
                                       + hybdat%prodm(i, n, l, itype)*cscal*atom_phase
                           ENDDO
                        END DO
                     ENDIF
                  END DO !n
               enddo  !j
            enddo !k
#ifdef _OPENACC
#else
            !$OMP END PARALLEL DO
#endif
            lm_0 = lm_0 + mpdata%num_radbasfn(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
         END DO
      END DO
      ! $acc exit data copyout(cprod, cprod%data_c) &
      ! $acc delete(hybdat, hybdat%nbands, hybdat%nindxp1, hybdat%gauntarr, hybdat%prodm) &
      ! $acc delete(bandoi, bandof, lmstart, lm_0, mpdata, mpdata%num_radfun_per_l, mpdata%l1, mpdata%l2, mpdata%n1, mpdata%n2) &
      ! $acc delete(cmt_ikqpt, cmt_nk, psize, ik)
      call timestop("loop over l, l1, l2, n, n1, n2")
      call timestop("wavefproducts_noinv5 MT")
   end subroutine wavefproducts_noinv_MT
end module m_wavefproducts_noinv
