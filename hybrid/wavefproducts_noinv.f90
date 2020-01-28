module m_wavefproducts_noinv
      USE m_types_hybdat

CONTAINS
   SUBROUTINE wavefproducts_noinv5(bandoi, bandof, ik, iq, &
                                    input, jsp, cell, atoms, mpdata, hybinp,&
                                   hybdat, kpts, lapw, sym, noco,&
                                   nkqpt, cprod)
      USE m_types
      use m_juDFT
      use m_constants, only: cmplx_0
      IMPLICIT NONE

      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_noco), INTENT(IN)        :: noco
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)     :: mpdata
      TYPE(t_hybinp), INTENT(IN)      :: hybinp
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat

!     - scalars -
      INTEGER, INTENT(IN)        ::  bandoi, bandof
      INTEGER, INTENT(IN)        ::  ik, iq, jsp
      INTEGER, INTENT(INOUT)     ::  nkqpt

!     - arrays -

      COMPLEX, INTENT(INOUT)    ::  cprod(hybdat%maxbasm1, bandoi:bandof, hybdat%nbands(ik))

      INTEGER        :: g_t(3)
      REAL           :: kqpt(3), kqpthlp(3)


      call timestart("wavefproducts_noinv5")
      cprod = cmplx_0; nkqpt = 0

      ! calculate nkpqt
      kqpthlp = kpts%bkf(:,ik) + kpts%bkf(:,iq)
      kqpt = kpts%to_first_bz(kqpthlp)
      g_t  = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      nkqpt = kpts%get_nk(kqpt)
      IF (.not. kpts%is_kpt(kqpt)) then
         call juDFT_error('wavefproducts: k-point not found')
      endif

      call wavefproducts_noinv5_IS(bandoi, bandof, ik, iq, g_t,&
                                         input, jsp, cell, atoms, mpdata, hybinp,&
                                        hybdat, kpts, lapw, sym, noco,&
                                        nkqpt, cprod)

      call wavefproducts_noinv_MT(bandoi, bandof, ik, iq, &
                                   input,atoms, mpdata, hybinp, hybdat, kpts, &
                                  nkqpt, cprod)
      call timestop("wavefproducts_noinv5")

   END SUBROUTINE wavefproducts_noinv5

   subroutine wavefproducts_noinv5_IS(bandoi, bandof, ik, iq, g_t, &
                                       input, jsp, cell, atoms, mpdata, hybinp,&
                                      hybdat, kpts, lapw, sym, noco,&
                                      nkqpt, cprod)
      use m_types
      use m_constants
      use m_wavefproducts_aux
      use m_judft
      use m_io_hybinp
      implicit NONE
      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_noco), INTENT(IN)        :: noco
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_mpdata), intent(in)  :: mpdata
      TYPE(t_hybinp), INTENT(IN)      :: hybinp
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat

!     - scalars -
      INTEGER, INTENT(IN)      ::  bandoi, bandof
      INTEGER, INTENT(IN)      ::  ik, iq, jsp, g_t(3)
      INTEGER, INTENT(IN)      ::  nkqpt

!     - arrays -

      COMPLEX, INTENT(INOUT)    ::  cprod(hybdat%maxbasm1, bandoi:bandof, hybdat%nbands(ik))

!     - local scalars -
      INTEGER                 :: ic, n1, n2
      INTEGER                 :: ig1, ig2, ig
      INTEGER                 :: igptm, iigptm, ngpt0, nbasfcn

      COMPLEX                 ::  cdum, cdum1

      TYPE(t_lapw)            ::  lapw_nkqpt

!      - local arrays -
      INTEGER                 ::  g(3)
      INTEGER, ALLOCATABLE    ::  gpt0(:,:)
      INTEGER, ALLOCATABLE    ::  pointer(:,:,:)



      COMPLEX                 ::  carr1(bandoi:bandof)
      COMPLEX                 ::  carr(bandoi:bandof, hybdat%nbands(ik))
      TYPE(t_mat)             ::  z_nk, z_kqpt
      COMPLEX, ALLOCATABLE    ::  z0(:,:)


      call timestart("wavefproducts_noinv5 IR")
      cprod = cmplx_0

      !
      ! compute G's fulfilling |bk(:,nkqpt) + G| <= rkmax
      !
      CALL lapw_nkqpt%init(input, noco, kpts, atoms, sym, nkqpt, cell, sym%zrfs)
      nbasfcn = calc_number_of_basis_functions(lapw, atoms, noco)
      call z_nk%alloc(.false., nbasfcn, input%neig)
      nbasfcn = calc_number_of_basis_functions(lapw_nkqpt, atoms, noco)
      call z_kqpt%alloc(.false., nbasfcn, input%neig)

      ! read in z at k-point ik and nkqpt
      call timestart("read_z")
      call read_z(atoms, cell, mpdata, hybdat, hybinp, kpts, sym, noco, input, lapw, ik, jsp, z_nk)
      call read_z(atoms, cell, mpdata, hybdat, hybinp, kpts, sym, noco, input, lapw, nkqpt, jsp, z_kqpt)
      call timestop("read_z")

      g = maxval(abs(lapw%gvec(:,:lapw%nv(jsp), jsp)), dim=2) &
        + maxval(abs(lapw_nkqpt%gvec(:,:lapw_nkqpt%nv(jsp), jsp)), dim=2)&
        + maxval(abs(mpdata%g(:,mpdata%gptm_ptr(:mpdata%n_g(iq), iq))), dim=2) + 1

      call hybdat%set_stepfunction(cell, atoms, g, sqrt(cell%omtil))

      !
      ! convolute phi(n,k) with the step function and store in cpw0
      !

      !(1) prepare list of G vectors
      call prep_list_of_gvec(lapw, mpdata, g, g_t, iq, jsp, pointer, gpt0, ngpt0)

      !(2) calculate convolution
      call timestart("calc convolution")
      call timestart("step function")
      ALLOCATE(z0(bandoi:bandof, ngpt0), source=cmplx_0)

      DO ig2 = 1, lapw_nkqpt%nv(jsp)
         if(z_kqpt%l_real) then
            carr1 = z_kqpt%data_r(ig2, bandoi:bandof)
         else
            carr1 = z_kqpt%data_c(ig2, bandoi:bandof)
         endif
         DO ig = 1, ngpt0
            g = gpt0(:,ig) - lapw_nkqpt%gvec(:,ig2, jsp)
            cdum = hybdat%stepfunc(g(1), g(2), g(3))
            DO n2 = bandoi, bandof
               z0(n2, ig) = z0(n2, ig) + carr1(n2)*cdum
            END DO
         END DO
      END DO
      call timestop("step function")

      call timestart("hybinp g")
      ic = hybdat%nbasp
      DO igptm = 1, mpdata%n_g(iq)
         carr = 0
         ic = ic + 1
         iigptm = mpdata%gptm_ptr(igptm, iq)

         DO ig1 = 1, lapw%nv(jsp)
            g = lapw%gvec(:,ig1, jsp) + mpdata%g(:,iigptm) - g_t
            ig2 = pointer(g(1), g(2), g(3))

            IF (ig2 == 0) call juDFT_error('wavefproducts_noinv2: pointer undefined')

            DO n1 = 1, hybdat%nbands(ik)
               if(z_nk%l_real) then
                  cdum1 = z_nk%data_r(ig1, n1)
               ELSE
                  cdum1 = conjg(z_nk%data_c(ig1, n1))
               endif
               DO n2 = bandoi, bandof
                  carr(n2, n1) = carr(n2, n1) + cdum1*z0(n2, ig2)
               END DO
            END DO

         END DO
         cprod(ic, :,:) = carr(:,:)
      END DO
      call timestop("hybinp g")
      deallocate(z0, pointer, gpt0)
      call timestop("calc convolution")

      call timestop("wavefproducts_noinv5 IR")
   end subroutine wavefproducts_noinv5_IS


   subroutine wavefproducts_noinv_MT(bandoi, bandof, ik, iq, &
                                      input,atoms, mpdata, hybinp, hybdat, kpts, &
                                     nkqpt, cprod)
      use m_types
      USE m_constants
      use m_io_hybinp
      use m_judft
      use m_wavefproducts_aux
      IMPLICIT NONE
      TYPE(t_input),INTENT(IN)         :: input
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybinp), INTENT(IN)      :: hybinp
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat

      !     - scalars -
      INTEGER, INTENT(IN)      ::  bandoi, bandof
      INTEGER, INTENT(IN)      ::  ik, iq
      INTEGER, INTENT(IN)     ::  nkqpt

      !     - arrays -

      COMPLEX, INTENT(INOUT)    ::  cprod(hybdat%maxbasm1, bandoi:bandof, hybdat%nbands(ik))

      !     - local scalars -
      INTEGER                 ::  ic, l, n, l1, l2, n1, n2, lm_0, lm1_0, lm2_0
      INTEGER                 ::  lm, lm1, lm2, m1, m2, i, ll, j, k
      INTEGER                 ::  itype, ieq, ic1, m

      COMPLEX                 ::  atom_phase

      LOGICAL                 ::  offdiag

      !      - local arrays -
      INTEGER                 ::  lmstart(0:atoms%lmaxd, atoms%ntype)

      COMPLEX                 ::  carr(bandoi:bandof, hybdat%nbands(ik))
      COMPLEX                 ::  cmt(input%neig, hybdat%maxlmindx, atoms%nat)
      COMPLEX                 ::  cmt_nk(input%neig, hybdat%maxlmindx, atoms%nat)

      call timestart("wavefproducts_noinv5 MT")
      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      call timestart("set lmstart")
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            lmstart(l, itype) = sum([(mpdata%num_radfun_per_l(ll, itype)*(2*ll+1), ll=0, l-1)])
         END DO
      END DO
      call timestop("set lmstart")

      ! read in cmt coefficients from direct access file cmt
      call timestart("read_cmt")
      call read_cmt(cmt_nk(:,:,:), ik)
      call read_cmt(cmt(:,:,:), nkqpt)
      call timestop("read_cmt")


      call timestart("loop over l, l1, l2, n, n1, n2")
      !$OMP PARALLEL PRIVATE(m, carr, lm1, m1, m2, lm2, i,j,k, &
      !$OMP lm, n1, l1, n2, l2, offdiag, lm1_0, lm2_0, itype, ieq, &
      !$OMP ic, lm_0)
      lm_0 = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            ic1 = 0

            atom_phase = exp(-ImagUnit*tpi_const*dot_product(kpts%bkf(:,iq), atoms%taual(:,ic)))

            DO l = 0, hybinp%lcutm1(itype)

               DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products
                  call mpdata%set_nl(n,l,itype, n1,l1,n2,l2)

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
                                             * outer_prod(cmt(bandoi:bandof, lm2, ic), &
                                                          conjg(cmt_nk(1:hybdat%nbands(ik), lm1, ic)))
                              END IF
                           END IF

                           m2 = m1 - m ! switch role of b1 and b2
                           IF (abs(m2) <= l2 .and. offdiag) THEN
                              lm2 = lm2_0 + n2 + (m2 + l2)*mpdata%num_radfun_per_l(l2, itype)
                              IF (abs(hybdat%gauntarr(2, l1, l2, l, m1, m)) > 1e-12) THEN
                                 carr = carr + hybdat%gauntarr(2, l1, l2, l, m1, m) &
                                             * outer_prod(cmt(bandoi:bandof, lm1, ic),&
                                                          conjg(cmt_nk(1:hybdat%nbands(ik), lm2, ic)))
                              END IF
                           END IF

                           lm1 = lm1 + mpdata%num_radfun_per_l(l1, itype) ! go to lm start index for next m1-quantum number

                        END DO  !m1

                        lm = lm_0 + (m+l) * mpdata%num_radbasfn(l,itype)
                        do k = 1,hybdat%nbands(ik)
                           do j = bandoi, bandof
                              DO i = 1, mpdata%num_radbasfn(l, itype)
                                 cprod(i+lm,j,k) = cprod(i+lm,j,k) &
                                       + hybdat%prodm(i, n, l, itype)*carr(j,k) *atom_phase
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
