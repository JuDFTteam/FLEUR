module m_wavefproducts_inv

CONTAINS
   SUBROUTINE wavefproducts_inv5(bandi, bandf, bandoi, bandof, dimension, input,&
                                 jsp, atoms, lapw, kpts, nk, iq, hybdat, hybrid,&
                                 cell, nbasm_mt, sym, noco, nkqpt, cprod)

      USE m_util, ONLY: modulo1
      USE m_types_hybrid, ONLY: gptnorm
      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_io_hybrid
      use m_wavefproducts_aux

      IMPLICIT NONE
      TYPE(t_dimension), INTENT(IN) :: dimension
      TYPE(t_hybrid), INTENT(IN)    :: hybrid
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_noco), INTENT(IN)      :: noco
      TYPE(t_sym), INTENT(IN)       :: sym
      TYPE(t_cell), INTENT(IN)      :: cell
      TYPE(t_kpts), INTENT(IN)      :: kpts
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_lapw), INTENT(IN)      :: lapw
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat

      ! - scalars -
      INTEGER, INTENT(IN)      :: bandi, bandf, bandoi, bandof
      INTEGER, INTENT(IN)      :: jsp, nk, iq
      INTEGER, INTENT(IN)      :: nbasm_mt
      INTEGER, INTENT(OUT)     :: nkqpt

      ! - arrays -
      REAL, INTENT(OUT)        ::    cprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)

      ! - local scalars -
      INTEGER                 ::    g_t(3)
      REAL                    ::    kqpt(3), kqpthlp(3)


      CALL timestart("wavefproducts_inv5")
      kqpthlp = kpts%bkf(:, nk) + kpts%bkf(:, iq)
      ! kqpt can lie outside the first BZ, transfer it back
      kqpt = kpts%to_first_bz(kqpthlp)
      g_t = nint(kqpt - kqpthlp)

      ! determine number of kqpt
      nkqpt = kpts%get_nk(kqpt)
      IF (.not. kpts%is_kpt(kqpt)) call juDFT_error('wavefproducts_inv5: k-point not found')


      ! 
      ! call wavefproducts_inv_IS_using_noinv(bandi, bandf, bandoi, bandof, dimension, input,&
      !                           jsp, atoms, lapw, kpts, nk, iq, g_t, hybdat, hybrid,&
      !                           cell, nbasm_mt, sym, noco, nkqpt, cprod)
      !

      call wavefproducts_inv_IS(bandi, bandf, bandoi, bandof, dimension, input,&
                                jsp, atoms, lapw, kpts, nk, iq, g_t, hybdat, hybrid,&
                                cell, nbasm_mt, sym, noco, nkqpt, cprod)
      ! call save_npy("cprod.npy", cprod)
      !
      ! call juDFT_error("i dont want anymore")
      call wavefproducts_inv5_MT(bandi, bandf, bandoi, bandof, dimension,&
                                atoms, kpts, nk, iq, hybdat, hybrid,&
                                sym, nkqpt, cprod)

      CALL timestop("wavefproducts_inv5")

   END SUBROUTINE wavefproducts_inv5

   subroutine wavefproducts_inv_IS_using_noinv(bandi, bandf, bandoi, bandof, dimension, input,&
                                 jsp, atoms, lapw, kpts, nk, iq, g_t, hybdat, hybrid,&
                                 cell, nbasm_mt, sym, noco, nkqpt, rprod)
     use m_types
     use m_wavefproducts_noinv
     use m_judft
     implicit NONE

     TYPE(t_dimension), INTENT(IN) :: dimension
     TYPE(t_hybrid), INTENT(IN)    :: hybrid
     TYPE(t_input), INTENT(IN)     :: input
     TYPE(t_noco), INTENT(IN)      :: noco
     TYPE(t_sym), INTENT(IN)       :: sym
     TYPE(t_cell), INTENT(IN)      :: cell
     TYPE(t_kpts), INTENT(IN)      :: kpts
     TYPE(t_atoms), INTENT(IN)     :: atoms
     TYPE(t_lapw), INTENT(IN)      :: lapw
     TYPE(t_hybdat), INTENT(INOUT) :: hybdat

     ! - scalars -
     INTEGER, INTENT(IN)  :: bandi, bandf, bandoi, bandof, jsp, nk, iq, nbasm_mt
     INTEGER, INTENT(IN)  :: nkqpt, g_t(3)

     ! - arrays -
     REAL, INTENT(OUT) :: rprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)
     COMPLEX :: cprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)

     call wavefproducts_noinv5_IS(bandi, bandf, bandoi, bandof, nk, iq, g_t,&
                                  dimension, input, jsp, cell, atoms, hybrid,&
                                  hybdat, kpts, lapw, sym, nbasm_mt, noco,&
                                  nkqpt, cprod)
    call save_npy("cprod_using_noinv.npy", cprod)
    rprod = real(cprod)
   end subroutine wavefproducts_inv_IS_using_noinv

   subroutine wavefproducts_inv_IS(bandi, bandf, bandoi, bandof, dimension, input,&
                                 jsp, atoms, lapw, kpts, nk, iq, g_t, hybdat, hybrid,&
                                 cell, nbasm_mt, sym, noco, nkqpt, cprod)
     use m_types
     use m_constants
     use m_wavefproducts_aux
     use m_judft
     use m_io_hybrid
     implicit NONE
     TYPE(t_dimension), INTENT(IN) :: dimension
     TYPE(t_hybrid), INTENT(IN)    :: hybrid
     TYPE(t_input), INTENT(IN)     :: input
     TYPE(t_noco), INTENT(IN)      :: noco
     TYPE(t_sym), INTENT(IN)       :: sym
     TYPE(t_cell), INTENT(IN)      :: cell
     TYPE(t_kpts), INTENT(IN)      :: kpts
     TYPE(t_atoms), INTENT(IN)     :: atoms
     TYPE(t_lapw), INTENT(IN)      :: lapw
     TYPE(t_hybdat), INTENT(INOUT) :: hybdat

     ! - scalars -
     INTEGER, INTENT(IN)      :: bandi, bandf, bandoi, bandof
     INTEGER, INTENT(IN)      :: jsp, nk, iq, g_t(3)
     INTEGER, INTENT(IN)      :: nbasm_mt
     INTEGER, INTENT(IN)      :: nkqpt

     ! - arrays -
     REAL, INTENT(OUT)        ::    cprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)

     ! - local scalars -
     INTEGER                 ::    ic, ig, ig2, ig1, ok, igptm, iigptm
     INTEGER                 ::    ngpt0, n1, n2, nbasfcn
     REAL                    ::    rdum, rdum1
     TYPE(t_lapw)            ::    lapw_nkqpt

     ! - local arrays -
     INTEGER, ALLOCATABLE    ::    pointer(:, :, :), gpt0(:, :)
     INTEGER                 ::    g(3)

     REAL                    ::    rarr1(bandoi:bandof)
     REAL                    ::    rarr2(bandoi:bandof, bandf - bandi + 1)
     REAL, ALLOCATABLE       ::    z0(:, :)

     TYPE(t_mat)             :: z_nk, z_kqpt

     CALL timestart("wavefproducts_inv5 IR")
     cprod = 0.0

     !
     ! compute G's fulfilling |bk(:,nkqpt) + G| <= rkmax
     !
     CALL lapw_nkqpt%init(input, noco, kpts, atoms, sym, nkqpt, cell, sym%zrfs)
     nbasfcn = calc_number_of_basis_functions(lapw, atoms, noco)
     call z_nk%alloc(.true., nbasfcn, dimension%neigd)
     nbasfcn = calc_number_of_basis_functions(lapw_nkqpt, atoms, noco)
     call z_kqpt%alloc(.true., nbasfcn, dimension%neigd)

     ! read in z at k-point nk and nkqpt
     call timestart("read_z")
     CALL read_z(z_nk, nk)
     call read_z(z_kqpt, nkqpt)
     call timestop("read_z")

     g = maxval(abs(lapw%gvec(:, :lapw%nv(jsp), jsp)), dim=2) &
    &  + maxval(abs(lapw_nkqpt%gvec(:, :lapw_nkqpt%nv(jsp), jsp)), dim=2)&
    &  + maxval(abs(hybrid%gptm(:, hybrid%pgptm(:hybrid%ngptm(iq), iq))), dim=2) + 1

     call hybdat%set_stepfunction(cell, atoms, g, sqrt(cell%omtil))
     !
     ! convolute phi(n,k) with the step function and store in cpw0
     !

     !(1) prepare list of G vectors
     call prep_list_of_gvec(lapw, hybrid, g, g_t, iq, jsp, pointer, gpt0, ngpt0)

     !(2) calculate convolution
     call timestart("calc convolution")
     ALLOCATE (z0(bandoi:bandof, ngpt0), stat=ok, source=0.0)
     IF (ok /= 0) call juDFT_error('wavefproducts_inv5: error allocation z0')

     call timestart("step function")
     DO ig2 = 1, lapw_nkqpt%nv(jsp)
        rarr1 = z_kqpt%data_r(ig2, bandoi:bandof)
        DO ig = 1, ngpt0
           g = gpt0(:, ig) - lapw_nkqpt%gvec(:, ig2, jsp)
           rdum = REAL(hybdat%stepfunc(g(1), g(2), g(3)))
           DO n2 = bandoi, bandof
              z0(n2, ig) = z0(n2, ig) + rarr1(n2)*rdum
           END DO
        END DO
     END DO
     call timestop("step function")

     call timestart("hybrid gptm")
     ic = nbasm_mt
     DO igptm = 1, hybrid%ngptm(iq)
        rarr2 = 0
        ic = ic + 1
        iigptm = hybrid%pgptm(igptm, iq)

        DO ig1 = 1, lapw%nv(jsp)
           g = lapw%gvec(:, ig1, jsp) + hybrid%gptm(:, iigptm) - g_t
           ig2 = pointer(g(1), g(2), g(3))

           IF (ig2 == 0) call juDFT_error('wavefproducts_inv5: pointer undefined')

           DO n1 = 1, bandf - bandi + 1
              rdum1 = z_nk%data_r(ig1, n1)
              DO n2 = bandoi, bandof
                 rarr2(n2, n1) = rarr2(n2, n1) + rdum1*z0(n2, ig2)
              END DO
           END DO

        END DO
        cprod(ic, :, :) = rarr2(:, :)
     END DO
     call timestop("hybrid gptm")
     call timestop("calc convolution")

     DEALLOCATE (z0, pointer, gpt0)
     CALL timestop("wavefproducts_inv5 IR")
   end subroutine wavefproducts_inv_IS

   subroutine wavefproducts_inv5_MT(bandi, bandf, bandoi, bandof, dimension,&
                                   atoms, kpts, nk, iq, hybdat, hybrid,&
                                   sym, nkqpt, cprod)
     use m_types
     use m_judft
     use m_io_hybrid
     use m_constants
     implicit NONE
     TYPE(t_dimension), INTENT(IN) :: dimension
     TYPE(t_hybrid), INTENT(IN)    :: hybrid
     TYPE(t_sym), INTENT(IN)       :: sym
     TYPE(t_kpts), INTENT(IN)      :: kpts
     TYPE(t_atoms), INTENT(IN)     :: atoms
     TYPE(t_hybdat), INTENT(INOUT) :: hybdat

     ! - scalars -
     INTEGER, INTENT(IN)      :: bandi, bandf, bandoi, bandof
     INTEGER, INTENT(IN)      :: nk, iq
     INTEGER, INTENT(IN)      :: nkqpt

     ! - arrays -
     REAL, INTENT(INOUT)        ::    cprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)

     ! - local scalars -
     INTEGER                 ::    i, iband
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
     INTEGER                 ::    lmstart(0:atoms%lmaxd, atoms%ntype)

     REAL                    ::    cmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat)
     REAL                    ::    cmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)
     REAL                    ::    rarr2(bandoi:bandof, bandf - bandi + 1)
     REAL                    ::    rarr3(2, bandoi:bandof, bandf - bandi + 1)

     COMPLEX                 ::    cmplx_exp(atoms%nat), cexp_nk(atoms%nat)
     COMPLEX, ALLOCATABLE    ::    ccmt_nk(:, :, :)
     COMPLEX, ALLOCATABLE    ::    ccmt(:, :, :)

     ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
     DO itype = 1, atoms%ntype
        DO l = 0, atoms%lmax(itype)
           lmstart(l, itype) = sum([(hybrid%nindx(ll, itype)*(2*ll + 1), ll=0, l - 1)])
        END DO
     END DO

     ! read in cmt coefficient at k-point nk
     ALLOCATE (ccmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat), &
               ccmt(dimension%neigd, hybrid%maxlmindx, atoms%nat), &
               source=cmplx(0.0, 0.0), stat=ok)
     IF (ok /= 0) call juDFT_error('wavefproducts_inv5: error allocation ccmt_nk/ccmt')

     call read_cmt(ccmt_nk, nk)
     !read in cmt coefficients at k+q point
     call read_cmt(ccmt, nkqpt)

     iatom = 0
     DO itype = 1, atoms%ntype
        DO ieq = 1, atoms%neq(itype)
           iatom = iatom + 1

           cmplx_exp(iatom) = exp((-img)*tpi_const* &
                                  dot_product(kpts%bkf(:, iq) + kpts%bkf(:, nk), atoms%taual(:, iatom)))

           cexp_nk(iatom) = exp((-img)*tpi_const* &
                                dot_product(kpts%bkf(:, nk), atoms%taual(:, iatom)))
        END DO
     END DO

     rfac = 1./sqrt(2.0)
     cfac = -img/sqrt(2.0)
     iatom = 0
     DO itype = 1, atoms%ntype
        DO ieq = 1, atoms%neq(itype)
           iatom = iatom + 1
           ! determine the number of the inverse symmetric atom belonging to iatom
           IF (sym%invsatnr(iatom) == 0) THEN
              iiatom = iatom
           ELSE
              iiatom = sym%invsatnr(iatom)
           END IF
           ! the cmt coefficients at iatom and iiatom are made real in one step
           IF (iiatom < iatom) CYCLE
           lm1 = 0
           DO l = 0, atoms%lmax(itype)
              DO m = -l, l
                 rdum = (-1)**(l + m)
                 DO p = 1, hybrid%nindx(l, itype)
                    lm1 = lm1 + 1
                    ! lm index at l,-m
                    lm2 = lm1 - 2*m*hybrid%nindx(l, itype)

                    IF (iatom == iiatom) THEN
                       IF (m < 0) THEN
                          cmt(:, lm1, iatom) = real((ccmt(:, lm1, iatom) + rdum*ccmt(:, lm2, iiatom))*cmplx_exp(iatom)*rfac)

                          cmt_nk(:, lm1, iatom) = real((ccmt_nk(:, lm1, iatom) + rdum*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*rfac)
                       ELSE IF (m > 0) THEN

                          cmt(:, lm1, iatom) = real((ccmt(:, lm1, iatom) - rdum*ccmt(:, lm2, iiatom))*cmplx_exp(iatom)*cfac)

                          cmt_nk(:, lm1, iatom) = real((ccmt_nk(:, lm1, iatom) - rdum*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*cfac)
                       ELSE
                          IF (mod(l, 2) == 0) THEN
                             cmt(:, lm1, iatom) = real(ccmt(:, lm1, iatom)*cmplx_exp(iatom))
                             cmt_nk(:, lm1, iatom) = real(ccmt_nk(:, lm1, iatom)*cexp_nk(iatom))
                          ELSE
                             cmt(:, lm1, iatom) = real(ccmt(:, lm1, iatom)*(-img)*cmplx_exp(iatom))
                             cmt_nk(:, lm1, iatom) = real(ccmt_nk(:, lm1, iatom)*(-img)*cexp_nk(iatom))
                          END IF
                       END IF
                    ELSE
                       cdum = rdum*cmplx_exp(iatom)*cmplx_exp(iiatom)
                       cmt(:, lm1, iatom) = real((ccmt(:, lm1, iatom) + cdum*ccmt(:, lm2, iiatom))*rfac)

                       cmt(:, lm1, iiatom) = real((ccmt(:, lm1, iatom) - cdum*ccmt(:, lm2, iiatom))*cfac)

                       cdum = rdum*cexp_nk(iatom)*cexp_nk(iiatom)
                       cmt_nk(:, lm1, iatom) = real((ccmt_nk(:, lm1, iatom) + cdum*ccmt_nk(:, lm2, iiatom))*rfac)

                       cmt_nk(:, lm1, iiatom) = real((ccmt_nk(:, lm1, iatom) - cdum*ccmt_nk(:, lm2, iiatom))*cfac)
                    END IF

                 END DO
              END DO
           END DO

        END DO
     END DO
     DEALLOCATE (ccmt_nk, ccmt)

     lm_0 = 0
     lm_00 = 0
     iatom1 = 0
     iiatom = 0

     DO itype = 1, atoms%ntype
        ioffset = sum([((2*ll + 1)*hybrid%nindxm1(ll, itype), ll=0, hybrid%lcutm1(itype))])
        lm_0 = lm_00
        DO ieq = 1, atoms%neq(itype)
           iatom1 = iatom1 + 1
           IF (sym%invsatnr(iatom1) == 0) THEN
              iatom2 = iatom1
           ELSE
              iatom2 = sym%invsatnr(iatom1)
           END IF
           IF (iatom1 > iatom2) CYCLE

           IF (iatom1 /= iatom2) THEN
              call timestart("iatom1 neq iatom2")
              ! loop over l of mixed basis
              DO l = 0, hybrid%lcutm1(itype)
                 ! loop over basis functions products, which belong to l
                 DO n = 1, hybdat%nindxp1(l, itype)

                    ! determine l1,p1 and l2,p2 for the basis functions, which can generate l
                    l1 = hybdat%prod(n, l, itype)%l1
                    l2 = hybdat%prod(n, l, itype)%l2
                    p1 = hybdat%prod(n, l, itype)%n1
                    p2 = hybdat%prod(n, l, itype)%n2

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
                             lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                             ! precalculated Gaunt coefficient
                             rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                             IF (abs(rdum) > 1e-12) THEN
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                   rdum2 = rdum*cmt_nk(iband, lmp1, iatom2)
                                   ! loop over occupied bands
                                   DO ibando = bandoi, bandof

                                      rarr3(1, ibando, iband) = rarr3(1, ibando, iband)&
               &                    + rdum1*cmt(ibando, lmp2, iatom1) + rdum2*cmt(ibando, lmp2, iatom2)

                                      rarr3(2, ibando, iband) = rarr3(2, ibando, iband)&
               &                    + rdum1*cmt(ibando, lmp2, iatom2) - rdum2*cmt(ibando, lmp2, iatom1)

                                   END DO  !ibando
                                END DO  !iband
                             END IF  ! rdum
                          END IF  ! abs(m2) .le. l2

                          m2 = m1 - m ! switch role of b1 and b2
                          IF (abs(m2) <= l2 .and. offdiag) THEN
                             lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                             rdum = hybdat%gauntarr(2, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                             IF (abs(rdum) > 1e-12) THEN
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                   rdum2 = rdum*cmt_nk(iband, lmp2, iatom2)
                                   ! loop over occupied bands
                                   DO ibando = bandoi, bandof
                                      rarr3(1, ibando, iband) = rarr3(1, ibando, iband)&
               &                    + rdum1*cmt(ibando, lmp1, iatom1) + rdum2*cmt(ibando, lmp1, iatom2)

                                      rarr3(2, ibando, iband) = rarr3(2, ibando, iband)&
               &                    + rdum1*cmt(ibando, lmp1, iatom2) - rdum2*cmt(ibando, lmp1, iatom1)
                                   END DO  !ibando
                                END DO  !iband
                             END IF  ! rdum .ne. 0
                          END IF  ! abs(m2) .le. l2 .and. offdiag

                          ! go to lmp start index for next m1-quantum number
                          lmp1 = lmp1 + hybrid%nindx(l1, itype)

                       END DO  !m1

                       ishift = -2*m*hybrid%nindxm1(l, itype)

                       ! go to lm mixed basis startindx for l and m
                       lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                       lm2 = lm + (iatom2 - 1 - iiatom)*ioffset + ishift

                       rdum = tpi_const*dot_product(kpts%bkf(:, iq), atoms%taual(:, iatom1))
                       rfac1 = sin(rdum)/sqrt(2.0)
                       rfac2 = cos(rdum)/sqrt(2.0)
                       DO iband = bandi, bandf
                          DO ibando = bandoi, bandof
                             rdum1 = rarr3(1, ibando, iband)
                             rdum2 = rarr3(2, ibando, iband)
                             add1 = rdum1*rfac2 + rdum2*rfac1
                             add2 = rdum2*rfac2 - rdum1*rfac1
                             DO i = 1, hybrid%nindxm1(l, itype)
                                j = lm1 + i
                                cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*add1
                                j = lm2 + i
                                cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*add2

                             END DO  !i -> loop over mixed basis functions
                          END DO  !ibando
                       END DO  !iband

                       ! go to lm start index for next m-quantum number
                       lm = lm + hybrid%nindxm1(l, itype)

                    END DO  !m

                 END DO !n
                 lm_0 = lm_0 + hybrid%nindxm1(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
                 IF (lm /= lm_0) call juDFT_error('wavefproducts_inv5: counting of lm-index incorrect (bug?)')
              END DO !l
              call timestop("iatom1 neq iatom2")
           ELSE !case: iatom1==iatom2
              call timestart("iatom1 eq iatom2")

              ! loop over l of mixed basis
              monepl = -1
              DO l = 0, hybrid%lcutm1(itype)
                 monepl = -monepl
                 ! loop over basis functions products, which belong to l
                 DO n = 1, hybdat%nindxp1(l, itype)

                    ! determine l1,p1 and l2,p2 for the basis functions, which can generate l
                    l1 = hybdat%prod(n, l, itype)%l1
                    l2 = hybdat%prod(n, l, itype)%l2
                    p1 = hybdat%prod(n, l, itype)%n1
                    p2 = hybdat%prod(n, l, itype)%n2

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
                          lmp1 = lp1 + l1*hybrid%nindx(l1, itype)
                          IF (mod(l1, 2) == 0) THEN
                             lmp2 = lp2 + (m + l2)*hybrid%nindx(l2, itype)
                          ELSE
                             lmp2 = lp2 + (-m + l2)*hybrid%nindx(l2, itype)
                          END IF

                          rdum = hybdat%gauntarr(1, l1, l2, l, 0, m)
                          IF (abs(rdum) > 1e-12) THEN
                             DO iband = bandi, bandf
                                rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                IF (mod(l1, 2) /= 0) rdum1 = moneplm*rdum1
                                DO ibando = bandoi, bandof
                                   rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                END DO  ! ibando
                             END DO  ! iband
                          END IF  ! rdum .ne. 0

                          IF (offdiag) THEN
                             rdum = hybdat%gauntarr(1, l2, l1, l, -m, m)
                             IF (abs(rdum) > 1e-12) THEN
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                   IF (mod(l1, 2) == 0) rdum1 = moneplm*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp1, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband
                             END IF  ! rdum .ne. 0
                          END IF  ! offdiag

                       END IF  ! abs(m) .le. l2

                       IF (abs(m) <= l1) THEN
                          IF (mod(l2, 2) == 0) THEN
                             lmp3 = lp1 + (m + l1)*hybrid%nindx(l1, itype)
                          ELSE
                             lmp3 = lp1 + (-m + l1)*hybrid%nindx(l1, itype)
                          END IF
                          lmp2 = lp2 + l2*hybrid%nindx(l2, itype)

                          rdum = hybdat%gauntarr(1, l1, l2, l, -m, m)
                          IF (abs(rdum) > 1e-12) THEN
                             DO iband = bandi, bandf
                                rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                IF (mod(l2, 2) == 0) rdum1 = moneplm*rdum1
                                DO ibando = bandoi, bandof
                                   rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                END DO  ! ibando
                             END DO  ! iband
                          END IF  ! rdum .ne. 0

                          IF (offdiag) THEN
                             rdum = hybdat%gauntarr(1, l2, l1, l, 0, m)
                             IF (abs(rdum) > 1e-12) THEN
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                   IF (mod(l2, 2) /= 0) rdum1 = moneplm*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                             lmp1 = lmp1 + hybrid%nindx(l1, itype)
                             CYCLE
                          END IF
                          ! (-1)**(l1+m1)
                          monepl1m1 = monepl1*monepm1
                          m2 = m1 + m
                          IF (abs(m2) <= l2 .and. m2 /= 0) THEN
                             rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                             IF (abs(rdum) > 1e-12) THEN
                                IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                ELSE
                                   lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                   fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                END IF
                                rdum = rdum/sqrt(2.0)
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                   IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband
                             END IF  ! rdum .ne. 0

                             IF (offdiag) THEN
                                rdum = hybdat%gauntarr(1, l2, l1, l, m2, -m)
                                IF (abs(rdum) > 1e-12) THEN
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                   IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                      lmp3 = lmp1
                                   ELSE
                                      lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                      fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                   END IF
                                   rdum = moneplm*rdum/sqrt(2.0)
                                   DO iband = bandi, bandf
                                      rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                      IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                      DO ibando = bandoi, bandof
                                         rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                ELSE
                                   lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                   fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                END IF
                                rdum = moneplm*rdum/sqrt(2.0)
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!moneplm*rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                   IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband

                             END IF  ! rdum .ne. 0

                             IF (offdiag) THEN
                                rdum = hybdat%gauntarr(1, l2, l1, l, m2, m)
                                IF (abs(rdum) > 1e-12) THEN
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                   IF (sign(1, m1) + sign(1, m2) /= 0) THEN
                                      lmp3 = lmp1
                                   ELSE
                                      lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                      fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                   END IF
                                   rdum = rdum/sqrt(2.0)
                                   DO iband = bandi, bandf
                                      rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                      IF (sign(1, m1) + sign(1, m2) == 0) rdum1 = fac*rdum1
                                      DO ibando = bandoi, bandof
                                         rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
                                      END DO  ! ibando
                                   END DO  ! iband
                                END IF  ! rdum .ne. 0
                             END IF  ! offdiag

                          END IF  ! abs(m2) .le. l2 .and. m1 .ne. 0

                          !go to lmp start index for next m1-quantum number
                          lmp1 = lmp1 + hybrid%nindx(l1, itype)
                       END DO  ! m1

                       ! go to lm mixed basis startindx for l and m
                       lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                       DO iband = bandi, bandf
                          DO ibando = bandoi, bandof
                             rdum = rarr2(ibando, iband)
                             DO i = 1, hybrid%nindxm1(l, itype)
                                j = lm1 + i
                                cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                             END DO  !i -> loop over mixed basis functions
                          END DO  !ibando
                       END DO  !iband

                       ! go to lm start index for next m-quantum number
                       lm = lm + hybrid%nindxm1(l, itype)

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
                             lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                             !lmp3 and lmp4 are variables, which avoid an if clause in the loop
                             lmp3 = lmp2
                             lmp4 = lmp1
                          ELSE
                             lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                             !lmp3 and lmp3 are variables, which avoid an if clause in the loop
                             lmp3 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                             lmp4 = lmp1 - 2*m1*hybrid%nindx(l1, itype)

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
                             DO iband = bandi, bandf
                                rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                DO ibando = bandoi, bandof
                                   rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
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
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp4, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband
                             END IF  ! rdum.ne.0
                          END IF  ! offdiag

                       END IF  ! abs(m2).le.l2

                       ! go to lmp start index for next m1-quantum number
                       lmp1 = lmp1 + hybrid%nindx(l1, itype)
                    END DO  !m1

                    ! go to lm mixed basis startindx for l and m
                    lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                    DO iband = bandi, bandf
                       DO ibando = bandoi, bandof
                          rdum = rarr2(ibando, iband)
                          DO i = 1, hybrid%nindxm1(l, itype)
                             j = lm1 + i
                             cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                          END DO  !i -> loop over mixed basis functions
                       END DO  !ibando
                    END DO  !iband

                    ! go to lm start index for next m-quantum number
                    lm = lm + hybrid%nindxm1(l, itype)

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
                          lmp1 = lp1 + l1*hybrid%nindx(l1, itype)
                          IF (mod(l1, 2) == 0) THEN
                             lmp2 = lp2 + (m + l2)*hybrid%nindx(l2, itype)
                          ELSE
                             lmp2 = lp2 + (-m + l2)*hybrid%nindx(l2, itype)
                          END IF

                          rdum = hybdat%gauntarr(1, l1, l2, l, 0, m)
                          IF (abs(rdum) > 1e-12) THEN
                             DO iband = bandi, bandf
                                rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                IF (mod(l1, 2) /= 0) rdum1 = moneplm*rdum1
                                DO ibando = bandoi, bandof
                                   rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                END DO  ! ibando
                             END DO  ! iband
                          END IF  ! rdum .ne. 0

                          IF (offdiag) THEN
                             rdum = hybdat%gauntarr(1, l2, l1, l, -m, m)
                             IF (abs(rdum) > 1e-12) THEN
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                   IF (mod(l1, 2) == 0) rdum1 = moneplm*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp1, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband
                             END IF  ! rdum .ne. 0
                          END IF  ! offdiag

                       END IF  ! abs(m) .le. l2

                       IF (abs(m) <= l1) THEN
                          IF (mod(l2, 2) == 0) THEN
                             lmp3 = lp1 + (m + l1)*hybrid%nindx(l1, itype)
                          ELSE
                             lmp3 = lp1 + (-m + l1)*hybrid%nindx(l1, itype)
                          END IF
                          lmp2 = lp2 + l2*hybrid%nindx(l2, itype)

                          rdum = hybdat%gauntarr(1, l1, l2, l, -m, m)
                          IF (abs(rdum) > 1e-12) THEN
                             DO iband = bandi, bandf
                                rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                IF (mod(l2, 2) == 0) rdum1 = moneplm*rdum1
                                DO ibando = bandoi, bandof
                                   rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                END DO  ! ibando
                             END DO  ! iband
                          END IF  ! rdum .ne. 0

                          IF (offdiag) THEN
                             rdum = hybdat%gauntarr(1, l2, l1, l, 0, m)
                             IF (abs(rdum) > 1e-12) THEN
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                   IF (mod(l2, 2) /= 0) rdum1 = moneplm*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                             lmp1 = lmp1 + hybrid%nindx(l1, itype)
                             CYCLE
                          END IF
                          m2 = m1 + m
                          ! (-1)**(l1+m1)
                          monepl1m1 = monepl1*monepm1
                          IF (abs(m2) <= l2 .and. m2 /= 0) THEN
                             rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                             IF (abs(rdum) > 1e-12) THEN

                                IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                   lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                ELSE
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                   fac = -moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))/2
                                END IF

                                rdum = -moneplm*monepl1m1*rdum/sqrt(2.0)
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!-moneplm*monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                   IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband

                             END IF  ! rdum .ne. 0

                             IF (offdiag) THEN
                                rdum = hybdat%gauntarr(2, l1, l2, l, m1, -m)
                                IF (abs(rdum) > 1e-12) THEN
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                   IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                      lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                   ELSE
                                      lmp3 = lmp1
                                      fac = 1/2.*monepl1m1*(sign(1, m2) - sign(1, m1))
                                   END IF
                                   rdum = monepl1m1*moneplm*rdum/sqrt(2.0)
                                   DO iband = bandi, bandf
                                      rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!monepl1m1*moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                      IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                      DO ibando = bandoi, bandof
                                         rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                                   lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                ELSE
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                   fac = 1/2.*moneplm*monepl1m1*(sign(1, m1) - sign(1, m2))
                                END IF
                                rdum = monepl1m1*rdum/sqrt(2.0)
                                DO iband = bandi, bandf
                                   rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sqrt(2.0)
                                   IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = rdum1*fac
                                   DO ibando = bandoi, bandof
                                      rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                   END DO  ! ibando
                                END DO  ! iband
                             END IF  ! rdum .ne. 0

                             IF (offdiag) THEN
                                rdum = hybdat%gauntarr(2, l1, l2, l, m1, m)
                                IF (abs(rdum) > 1e-12) THEN
                                   lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                   IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                      lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                   ELSE
                                      lmp3 = lmp1
                                      fac = -monepl1m1*(sign(1, m1) - sign(1, m2))/2
                                   END IF
                                   rdum = -monepl1m1*rdum/sqrt(2.0)
                                   DO iband = bandi, bandf
                                      rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!-monepl1m1*rdum*cmt_nk(iband,lmp2,iatom1)/sqrt(2.0)
                                      IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                      DO ibando = bandoi, bandof
                                         rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
                                      END DO  ! ibando
                                   END DO  ! iband

                                END IF  ! rdum .ne. 0
                             END IF  ! offdiag

                          END IF  !  abs(m2) .le. l2 .and. m2 .ne. 0

                          !go to lmp start index for next m1-quantum number
                          lmp1 = lmp1 + hybrid%nindx(l1, itype)
                       END DO  ! m1

                       ! multiply rarr2 by (-1)**(l+m+1)
                       rarr2(:, :) = (-1)*moneplm*rarr2(:, :)

                       ! go to lm mixed basis startindx for l and m
                       lm1 = lm + (iatom1 - 1 - iiatom)*ioffset

                       DO iband = bandi, bandf
                          DO ibando = bandoi, bandof
                             rdum = rarr2(ibando, iband)
                             DO i = 1, hybrid%nindxm1(l, itype)
                                j = lm1 + i
                                cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                             END DO  !i -> loop over mixed basis functions
                          END DO  !ibando
                       END DO  !iband

                       ! go to lm start index for next m-quantum number
                       lm = lm + hybrid%nindxm1(l, itype)

                    END DO  ! m=1,l

                 END DO !n
                 lm_0 = lm_0 + hybrid%nindxm1(l, itype)*(2*l + 1) ! go to the m start index of the next l-quantum number
                 IF (lm /= lm_0) call juDFT_error('wavefproducts_inv5: counting of lm-index incorrect (bug?)')
              END DO !l

              call timestop("iatom1 eq iatom2")
           END IF  ! iatom1 .ne. iatom2

           lm_0 = lm_00
        END DO !ieq
        iiatom = iiatom + atoms%neq(itype)
        lm_00 = lm_00 + atoms%neq(itype)*ioffset
     END DO  !itype
   end subroutine wavefproducts_inv5_MT
end module m_wavefproducts_inv
