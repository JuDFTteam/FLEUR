!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wavefproducts
   USE m_judft
   PRIVATE
   PUBLIC wavefproducts_noinv, wavefproducts_noinv5
   PUBLIC wavefproducts_inv, wavefproducts_inv5
CONTAINS

   SUBROUTINE wavefproducts_noinv(bandi, bandf, nk, iq, dimension, input, jsp,&                  !cprod,&
  &                 cell, atoms, hybrid, hybdat,&
  &                 kpts, mnobd,&
  &                 lapw, sym, noco, nbasm_mt, nkqpt, cprod)

      USE m_constants
      USE m_util, ONLY: modulo1
      USE m_wrapper
      USE m_types
      USE m_io_hybrid
      IMPLICIT NONE

      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_hybrid), INTENT(IN)      :: hybrid
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_noco), INTENT(IN)        :: noco
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat

!     - scalars -
      INTEGER, INTENT(IN)      ::  nk, iq, jsp
      INTEGER, INTENT(IN)      :: mnobd
      INTEGER, INTENT(IN)      :: nbasm_mt
      INTEGER, INTENT(IN)      ::  bandi, bandf
      INTEGER, INTENT(OUT)     ::  nkqpt

!     - arrays -

      COMPLEX, INTENT(OUT)    ::  cprod(hybrid%maxbasm1, mnobd, bandf - bandi + 1)

!     - local scalars -
      INTEGER                 ::  ic, l, n, l1, l2, n1, n2, lm_0, lm1_0, lm2_0, lm, lm1, lm2, m1, m2, i, j, ll
      INTEGER                 ::  itype, ieq, ikpt, ikpt1, ikpt2, igpt, igptp, igpt1, igpt2, iband, iband1, iband2
      INTEGER                 ::  k, ic1, ioffset, ibando
      INTEGER                 ::  q, idum, m
      INTEGER                 ::  nbasm_ir
      INTEGER                 ::  nbasmmt, nbasfcn
      INTEGER                 ::  ok
      REAL                    ::  rdum, svol, s2, pi
      REAL                    ::  mtthr = 0
      COMPLEX                 ::  cdum, cdum0
      COMPLEX, PARAMETER       ::  img = (0.0, 1.0)

      LOGICAL                 ::  offdiag
!      - local arrays -
      INTEGER                 ::  iarr(lapw%nv(jsp))
      INTEGER                 ::  g(3), ghelp(3), lmstart(0:atoms%lmaxd, atoms%ntype)
      INTEGER                 ::  gpthlp(3, dimension%nvd), nvhlp(input%jspins)
      INTEGER                 ::  gpt_nk(3, lapw%nv(jsp))
      INTEGER                 ::  g_t(3)
      INTEGER                 ::  gsum(3)
      REAL                    ::  bkpt(3)
      REAL                    ::  bkhlp(3)
      REAL                    ::  kqpt(3), kqpthlp(3)
      COMPLEX                 ::  carr(1:mnobd, bandf - bandi + 1)
!      COMPLEX                 :: chelp(maxbasm,mnobd,bandf-bandi+1,nkpt_EIBZ)
      COMPLEX                 ::  cexp
      COMPLEX                 ::  z_help(lapw%nv(jsp))
      COMPLEX                 ::  cmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      COMPLEX                 ::  cmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      COMPLEX, ALLOCATABLE     ::  cprod_ir(:, :, :)
      TYPE(t_mat)             :: z_nk, z_kqpt
      TYPE(t_lapw)            :: lapw_nkqpt
      CALL timestart("wavefproducts_noinv")

      ! preparations

      svol = sqrt(cell%omtil)
      s2 = sqrt(2.0)

      gpt_nk(1, :) = lapw%k1(:lapw%nv(jsp), jsp)
      gpt_nk(2, :) = lapw%k2(:lapw%nv(jsp), jsp)
      gpt_nk(3, :) = lapw%k3(:lapw%nv(jsp), jsp)

      !
      ! compute k+q point for given q point in EIBZ(k)
      !
      kqpthlp = kpts%bkf(:, nk) + kpts%bkf(:, iq)
      ! k+q can lie outside the first BZ, transfer
      ! it back into the 1. BZ
      kqpt = modulo1(kqpthlp, kpts%nkpt3)
      g_t(:) = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      nkqpt = 0
      DO ikpt = 1, kpts%nkptf
         IF (maxval(abs(kqpt - kpts%bkf(:, ikpt))) <= 1E-06) THEN
            nkqpt = ikpt
            EXIT
         END IF
      END DO
      IF (nkqpt == 0) STOP 'wavefproducts: k-point not found'

      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            lmstart(l, itype) = sum((/(hybrid%nindx(ll, itype)*(2*ll + 1), ll=0, l - 1)/))
         END DO
      END DO

      nbasm_ir = maxval(mpbasis%ngptm)
      ALLOCATE (cprod_ir(bandf - bandi + 1, mnobd, nbasm_ir), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts: failure allocation cprod_ir'
      cprod_ir = 0

      cprod = 0

      CALL lapw_nkqpt%init(input, noco, kpts, atoms, sym, nkqpt, cell, sym%zrfs)
      nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
      call z_nk%alloc(.false., nbasfcn, dimension%neigd)
      nbasfcn = MERGE(lapw_nkqpt%nv(1) + lapw_nkqpt%nv(2) + 2*atoms%nlotot, lapw_nkqpt%nv(1) + atoms%nlotot, noco%l_noco)
      call z_kqpt%alloc(.false., nbasfcn, dimension%neigd)

      call read_z(z_nk, kpts%nkptf*(jsp - 1) + nk)
      call read_z(z_kqpt, kpts%nkptf*(jsp - 1) + nkqpt)

      ! read in cmt coefficients from direct access file cmt
      call read_cmt(cmt_nk, nk)

      ! IR contribution

      CALL timestart("wavefproducts_noinv IR")

      DO igpt = 1, mpbasis%ngptm(iq)
         igptp = mpbasis%gptm_ptr(igpt, iq)
         ghelp = mpbasis%gptm(:, igptp) - g_t(:)
         DO i = 1, lapw%nv(jsp)
            gsum(:) = ghelp + gpt_nk(:, i)
            IF (all(abs(gsum) <= hybdat%pntgptd)) THEN
               iarr(i) = hybdat%pntgpt(gsum(1), gsum(2), gsum(3), nkqpt)
            ELSE
               iarr(i) = 0
            END IF

         END DO

         DO iband1 = 1, hybrid%nobd(nkqpt,jsp)
            where (iarr > 0)
            z_help(:) = z_kqpt%data_c(iarr(:), iband1)
            elsewhere
            z_help = 0.0
            end where

            DO iband = bandi, bandf
               cprod_ir(iband, iband1, igpt) = 1/svol*dotprod(z_nk%data_c(:lapw%nv(jsp), iband), z_help)
            END DO !iband

         END DO  !iband1

      END DO !igpt
      CALL timestop("wavefproducts_noinv IR")

      !
      ! MT contribution
      !
      call read_cmt(cmt, nkqpt)
      lm_0 = 0
      ic = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            ic1 = 0

            cexp = exp(-2*img*pi_const*dot_product(kpts%bkf(:, iq), atoms%taual(:, ic)))

            DO l = 0, hybrid%lcutm1(itype)
               DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products

                  l1 = hybdat%prod(n, l, itype)%l1 !
                  l2 = hybdat%prod(n, l, itype)%l2 ! current basis-function product
                  n1 = hybdat%prod(n, l, itype)%n1 ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
                  n2 = hybdat%prod(n, l, itype)%n2 !

                  IF (mod(l1 + l2 + l, 2) /= 0) cycle

                  offdiag = l1 /= l2 .or. n1 /= n2 ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                  !(leading to the same basis-function product)

                  lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                  lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                  lm = lm_0
                  DO m = -l, l

                     carr = 0.0

                     lm1 = lm1_0 + n1 ! go to lm index for m1=-l1
                     DO m1 = -l1, l1
                        m2 = m1 + m ! Gaunt condition -m1+m2-m=0
                        IF (abs(m2) <= l2) THEN
                           lm2 = lm2_0 + n2 + (m2 + l2)*hybrid%nindx(l2, itype)
                           rdum = hybdat%gauntarr(1, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 cdum = rdum*conjg(cmt_nk(iband, lm1, ic)) !nk
                                 DO iband1 = 1, mnobd
                                    carr(iband1, iband) = carr(iband1, iband) + cdum*cmt(iband1, lm2, ic) !ikpt

                                 END DO
                              END DO
                           END IF
                        END IF

                        m2 = m1 - m ! switch role of b1 and b2
                        IF (abs(m2) <= l2 .and. offdiag) THEN
                           lm2 = lm2_0 + n2 + (m2 + l2)*hybrid%nindx(l2, itype)
                           rdum = hybdat%gauntarr(2, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 cdum = rdum*conjg(cmt_nk(iband, lm2, ic)) !nk
                                 DO iband1 = 1, mnobd
                                    carr(iband1, iband) = carr(iband1, iband) + cdum*cmt(iband1, lm1, ic)
                                 END DO
                              END DO
                           END IF
                        END IF

                        lm1 = lm1 + hybrid%nindx(l1, itype) ! go to lm start index for next m1-quantum number

                     END DO  !m1

                     DO iband = bandi, bandf
                        DO iband1 = 1, mnobd
                           cdum = carr(iband1, iband)*cexp
                           DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                              j = lm + i
                              cprod(j, iband1, iband) = cprod(j, iband1, iband) + hybdat%prodm(i, n, l, itype)*cdum
                           END DO

                        END DO
                     END DO

                     lm = lm + mpbasis%num_rad_bas_fun(l, itype) ! go to lm start index for next m-quantum number

                  END DO

               END DO
               lm_0 = lm_0 + mpbasis%num_rad_bas_fun(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               IF (lm /= lm_0) STOP 'wavefproducts: counting of lm-index incorrect (bug?)'
            END DO
         END DO
      END DO

      CALL timestop("wavefproducts_noinv")

      ic = nbasm_mt
      DO igpt = 1, mpbasis%ngptm(iq)
         ic = ic + 1
         DO ibando = 1, mnobd
            DO iband = bandi, bandf
               cprod(ic, ibando, iband) = cprod_ir(iband, ibando, igpt)
            END DO
         END DO
      END DO

   END SUBROUTINE wavefproducts_noinv

   SUBROUTINE wavefproducts_inv(&
  &                  bandi, bandf, dimension, input, jsp, atoms,&
  &                  lapw, kpts,&
  &                  nk, iq, hybdat, mnobd, hybrid,&
  &                  parent, cell,&
  &                  nbasm_mt, sym, noco,&
  &                  nkqpt, cprod)

      USE m_util, ONLY: modulo1
      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_io_hybrid
      IMPLICIT NONE
      TYPE(t_hybdat), INTENT(IN)   :: hybdat
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_noco), INTENT(IN)   :: noco
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_lapw), INTENT(IN)   :: lapw

      ! - scalars -
      INTEGER, INTENT(IN)      ::    bandi, bandf
      INTEGER, INTENT(IN)      ::    jsp, nk, iq
      INTEGER, INTENT(IN)      ::    mnobd
      INTEGER, INTENT(IN)      ::    nbasm_mt
      INTEGER, INTENT(OUT)     ::    nkqpt

      ! - arrays -
      INTEGER, INTENT(IN)      ::    parent(kpts%nkptf)

      REAL, INTENT(OUT)       ::    cprod(hybrid%maxbasm1, mnobd, bandf - bandi + 1)

      ! - local scalars -
      INTEGER                 ::    i, ikpt, ic, iband, iband1, igpt, igptp, ibando, iatom, iiatom, itype, ieq, ishift, ioffset, iatom1, iatom2
      INTEGER                 ::    l, p, l1, m1, l2, m2, p1, p2, n, ok
      INTEGER                 ::    lm, lm1, lm2, lm_0, lm_00, lm1_0, lm2_0, lmp1, lmp2, lmp3, lmp4, lp1, lp2
      INTEGER                 ::    j, ll, m, nbasfcn
      INTEGER                 :: nbasm_ir
      REAL                    ::    svol, sr2
      REAL                    ::    rdum, rfac, rfac1, rfac2, rdum1, rdum2
      REAL                    ::    sin1, sin2, cos1, cos2, add1, add2
      REAL                    ::    fac, fac1, fac2
      REAL                    ::    monepl1, monepl2, monepl, monepm1, monepm, moneplm, monepl1m1
      COMPLEX, PARAMETER       ::    img = (0.0, 1.0)
      COMPLEX                 ::    fexp
      COMPLEX                 ::    cdum, cconst, cfac
      LOGICAL                 ::    offdiag

      ! - local arrays -
      INTEGER                 ::    iarr(lapw%nv(jsp))
      INTEGER                 ::    gpt_nk(3, lapw%nv(jsp)), ghelp(3)
      INTEGER                 ::    gsum(3)
      INTEGER                 ::    g_t(3)
      INTEGER                 ::    lmstart(0:atoms%lmaxd, atoms%ntype)
      INTEGER                 ::    lmstart2(0:maxval(hybrid%lcutm1), atoms%nat)
      REAL                    ::    kqpt(3), kqpthlp(3)

      REAL, ALLOCATABLE        ::    cprod_ir(:, :, :)

      REAL                    ::    z_help(lapw%nv(jsp))

      REAL                    ::    cmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      REAL                    ::    cmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)

      COMPLEX                 ::    ccmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      COMPLEX                 ::    ccmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)

      REAL                    ::    rarr1(1:mnobd, bandf - bandi + 1)
      REAL                    ::    rarr(2, 1:mnobd, bandf - bandi + 1)
      COMPLEX                 ::    cmthlp(dimension%neigd), cmthlp1(dimension%neigd)
      COMPLEX                 ::    cexp(atoms%nat), cexp_nk(atoms%nat)
      TYPE(t_mat)             :: z_nk, z_kqpt
      TYPE(t_lapw)            :: lapw_nkqpt

      CALL timestart("wavefproducts_inv")
      CALL timestart("wavefproducts_inv IR")
      svol = sqrt(cell%omtil)
      sr2 = sqrt(2.0)

      nbasm_ir = maxval(mpbasis%ngptm)
      ALLOCATE (cprod_ir(bandf - bandi + 1, mnobd, nbasm_ir))
      cprod_ir = 0
      gpt_nk(1, :) = lapw%k1(:lapw%nv(jsp), jsp)
      gpt_nk(2, :) = lapw%k2(:lapw%nv(jsp), jsp)
      gpt_nk(3, :) = lapw%k3(:lapw%nv(jsp), jsp)

      !
      ! compute k+q point for q (iq) in EIBZ(k)
      !

      kqpthlp = kpts%bkf(:, nk) + kpts%bkf(:, iq)
      ! kqpt can lie outside the first BZ, transfer it back
      kqpt = modulo1(kqpthlp, kpts%nkpt3)
      g_t(:) = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      nkqpt = 0
      DO ikpt = 1, kpts%nkptf
         IF (maxval(abs(kqpt - kpts%bkf(:, ikpt))) <= 1E-06) THEN
            nkqpt = ikpt
            EXIT
         END IF
      END DO
      IF (nkqpt == 0) STOP 'wavefproducts: k-point not found'

      ! read in z at current k-point nk

      CALL lapw_nkqpt%init(input, noco, kpts, atoms, sym, nkqpt, cell, sym%zrfs)
      nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
      call z_nk%alloc(.true., nbasfcn, dimension%neigd)
      nbasfcn = MERGE(lapw_nkqpt%nv(1) + lapw_nkqpt%nv(2) + 2*atoms%nlotot, lapw_nkqpt%nv(1) + atoms%nlotot, noco%l_noco)
      call z_kqpt%alloc(.true., nbasfcn, dimension%neigd)

      call read_z(z_nk, kpts%nkptf*(jsp - 1) + nk)
      call read_z(z_kqpt, kpts%nkptf*(jsp - 1) + nkqpt)

      DO igpt = 1, mpbasis%ngptm(iq)
         igptp = mpbasis%gptm_ptr(igpt, iq)
         ghelp = mpbasis%gptm(:, igptp) - g_t(:)
         DO i = 1, lapw%nv(jsp)
            gsum(:) = ghelp + gpt_nk(:, i)
            IF (all(abs(gsum) <= hybdat%pntgptd)) THEN
               iarr(i) = hybdat%pntgpt(gsum(1), gsum(2), gsum(3), nkqpt)
            ELSE
               iarr(i) = 0
            END IF

         END DO

         DO iband1 = 1, hybrid%nobd(nkqpt,jsp)
            where (iarr > 0)
            z_help(:) = z_kqpt%data_r(iarr(:), iband1)
            elsewhere
            z_help = 0.0
            end where
            DO iband = bandi, bandf
               cprod_ir(iband, iband1, igpt) = 1/svol*dotprod(z_nk%data_r(:lapw%nv(jsp), iband), z_help)

            END DO !iband

         END DO  !iband1

      END DO !igpt

      CALL timestop("wavefproducts_inv IR")

      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            lmstart(l, itype) = sum((/(hybrid%nindx(ll, itype)*(2*ll + 1), ll=0, l - 1)/))
         END DO
      END DO

      ! read in cmt coefficient at k-point nk

      CALL read_cmt(ccmt_nk, nk)

      !read in cmt coefficients at k+q point
      call read_cmt(ccmt, nkqpt)

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            cexp(iatom) = exp((-img)*tpi_const*dotprod(kpts%bkf(:, iq) + kpts%bkf(:, nk), atoms%taual(:, iatom)))

            cexp_nk(iatom) = exp((-img)*tpi_const*dotprod(kpts%bkf(:, nk), atoms%taual(:, iatom)))
         END DO
      END DO

      rfac = 1./sr2
      cfac = -img/sr2
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
                  DO p = 1, hybrid%nindx(l, itype)
                     lm1 = lm1 + 1
                     ! lm index at l,-m
                     lm2 = lm1 - 2*m*hybrid%nindx(l, itype)

                     IF (iatom == iiatom) THEN
                        IF (m < 0) THEN
                           cmt(:, lm1, iatom) = (ccmt(:, lm1, iatom) + (-1)**(l + m)*ccmt(:, lm2, iiatom))*cexp(iatom)*rfac

                           cmt_nk(:, lm1, iatom) = (ccmt_nk(:, lm1, iatom) + (-1)**(l + m)*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*rfac
                        ELSE IF (m > 0) THEN

                           cmt(:, lm1, iatom) = (ccmt(:, lm1, iatom) - (-1)**(l + m)*ccmt(:, lm2, iiatom))*cexp(iatom)*cfac

                           cmt_nk(:, lm1, iatom) = (ccmt_nk(:, lm1, iatom) - (-1)**(l + m)*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*cfac
                        ELSE
                           IF (mod(l, 2) == 0) THEN
                              cmt(:, lm1, iatom) = ccmt(:, lm1, iatom)*cexp(iatom)
                              cmt_nk(:, lm1, iatom) = ccmt_nk(:, lm1, iatom)*cexp_nk(iatom)
                           ELSE
                              cmt(:, lm1, iatom) = ccmt(:, lm1, iatom)*(-img)*cexp(iatom)
                              cmt_nk(:, lm1, iatom) = ccmt_nk(:, lm1, iatom)*(-img)*cexp_nk(iatom)
                           END IF
                        END IF
                     ELSE
                        cmt(:, lm1, iatom) = (ccmt(:, lm1, iatom) + (-1)**(l + m)*ccmt(:, lm2, iiatom))*rfac

                        cmt(:, lm1, iiatom) = (ccmt(:, lm1, iatom) - (-1)**(l + m)*ccmt(:, lm2, iiatom))*cfac

                        cmt_nk(:, lm1, iatom) = (ccmt_nk(:, lm1, iatom) + (-1)**(l + m)*ccmt_nk(:, lm2, iiatom))*rfac

                        cmt_nk(:, lm1, iiatom) = (ccmt_nk(:, lm1, iatom) - (-1)**(l + m)*ccmt_nk(:, lm2, iiatom))*cfac
                     END IF

                  END DO
               END DO
            END DO

         END DO
      END DO

      cprod = 0.0

      lm_0 = 0
      lm_00 = 0
      iatom1 = 0
      iiatom = 0

      DO itype = 1, atoms%ntype
         ioffset = sum((/((2*ll + 1)*mpbasis%num_rad_bas_fun(ll, itype), ll=0, hybrid%lcutm1(itype))/))
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
                        rarr = 0.0

                        ! go to lm index for m1=-l1
                        lmp1 = lm1_0 + p1

                        DO m1 = -l1, l1
                           ! Gaunt condition -m1+m2-m=0
                           m2 = m1 + m
                           IF (abs(m2) <= l2) THEN
                              lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                              ! precalculated Gaunt coefficient
                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, m)
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                    rdum2 = rdum*cmt_nk(iband, lmp1, iatom2)
                                    ! loop over occupied bands
                                    DO ibando = 1, mnobd!hybrid%nobd(peibz(ikpt))

                                       rarr(1, ibando, iband) = rarr(1, ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1) + rdum2*cmt(ibando, lmp2, iatom2)

                                       rarr(2, ibando, iband) = rarr(2, ibando, iband) + rdum1*cmt(ibando, lmp2, iatom2) - rdum2*cmt(ibando, lmp2, iatom1)

                                    END DO  !ibando
                                 END DO  !iband
                              END IF  ! rdum
                           END IF  ! abs(m2) .le. l2

                           m2 = m1 - m ! switch role of b1 and b2
                           IF (abs(m2) <= l2 .and. offdiag) THEN
                              lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                              rdum = hybdat%gauntarr(2, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    rdum2 = rdum*cmt_nk(iband, lmp2, iatom2)
                                    ! loop over occupied bands
                                    DO ibando = 1, mnobd!hybrid%nobd(peibz(ikpt)
                                       rarr(1, ibando, iband) = rarr(1, ibando, iband) + rdum1*cmt(ibando, lmp1, iatom1) + rdum2*cmt(ibando, lmp1, iatom2)

                                       rarr(2, ibando, iband) = rarr(2, ibando, iband) + rdum1*cmt(ibando, lmp1, iatom2) - rdum2*cmt(ibando, lmp1, iatom1)
                                    END DO  !ibando
                                 END DO  !iband
                              END IF  ! rdum .ne. 0
                           END IF  ! abs(m2) .le. l2 .and. offdiag

                           ! go to lmp start index for next m1-quantum number
                           lmp1 = lmp1 + hybrid%nindx(l1, itype)

                        END DO  !m1

                        ishift = -2*m*mpbasis%num_rad_bas_fun(l, itype)

                        ! go to lm mixed basis startindx for l and m
                        lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                        lm2 = lm + (iatom2 - 1 - iiatom)*ioffset + ishift

                        rdum = tpi_const*dotprod(kpts%bkf(:, iq), atoms%taual(:, iatom1))
                        rfac1 = sin(rdum)/sr2
                        rfac2 = cos(rdum)/sr2
                        DO iband = bandi, bandf
                           DO ibando = 1, mnobd
                              rdum1 = rarr(1, ibando, iband)
                              rdum2 = rarr(2, ibando, iband)
!                       sin1  = rdum1*rfac1
!                       cos1  = rdum1*rfac2
!                       sin2  = rdum2*rfac1
!                       cos2  = rdum2*rfac2
                              add1 = rdum1*rfac2 + rdum2*rfac1
                              add2 = rdum2*rfac2 - rdum1*rfac1
                              DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                                 j = lm1 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*add1!( cos1 + sin2 )
                                 j = lm2 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*add2!( cos2 - sin1 )

                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpbasis%num_rad_bas_fun(l, itype)

                     END DO  !m

                  END DO !n
                  lm_0 = lm_0 + mpbasis%num_rad_bas_fun(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
                  IF (lm /= lm_0) STOP 'wavefproducts: counting of lm-index incorrect (bug?)'
               END DO !l

            ELSE !case: iatom1==iatom2

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
                        rarr1 = 0.0
                        IF (abs(m) <= l2) THEN
                           lmp1 = lp1 + l1*hybrid%nindx(l1, itype)
                           IF (mod(l1, 2) == 0) THEN
                              lmp2 = lp2 + (m + l2)*hybrid%nindx(l2, itype)
                           ELSE
                              lmp2 = lp2 + (-m + l2)*hybrid%nindx(l2, itype)
                           END IF

                           rdum = hybdat%gauntarr(1, l1, l2, l, 0, m)
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                 IF (mod(l1, 2) /= 0) rdum1 = moneplm*rdum1
                                 DO ibando = 1, mnobd
                                    rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, -m, m)
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l1, 2) == 0) rdum1 = moneplm*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp1, iatom1)
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
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                 IF (mod(l2, 2) == 0) rdum1 = moneplm*rdum1
                                 DO ibando = 1, mnobd
                                    rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, 0, m)
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l2, 2) /= 0) rdum1 = moneplm*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                              IF (rdum /= 0) THEN
                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                 END IF
                                 rdum = rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(1, l2, l1, l, m2, -m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1
                                    ELSE
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                       fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                    END IF
                                    rdum = moneplm*rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sr2
                                       IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                       DO ibando = 1, mnobd
                                          rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband
                                 END IF  ! rdum .ne. 0
                              END IF  ! offdiag

                           END IF  ! abs(m2) .le. l2 .and. m2 .ne. 0

                           m2 = m1 - m
                           IF (abs(m2) <= l2 .and. m2 /= 0) THEN

                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, -m)
                              IF (rdum /= 0) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                 END IF
                                 rdum = moneplm*rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!moneplm*rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband

                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(1, l2, l1, l, m2, m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m1) + sign(1, m2) /= 0) THEN
                                       lmp3 = lmp1
                                    ELSE
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                       fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                    END IF
                                    rdum = rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!rdum*cmt_nk(iband,lmp2,iatom1)/sr2
                                       IF (sign(1, m1) + sign(1, m2) == 0) rdum1 = fac*rdum1
                                       DO ibando = 1, mnobd
                                          rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                           DO ibando = 1, mnobd
                              rdum = rarr1(ibando, iband)
                              DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                                 j = lm1 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpbasis%num_rad_bas_fun(l, itype)

                     END DO  ! m=-l,-1

                     !
                     !case m=0
                     !

                     m = 0
                     rarr1 = 0.0
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
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                 DO ibando = 1, mnobd
                                    rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
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
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp4, iatom1)
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
                        DO ibando = 1, mnobd
                           rdum = rarr1(ibando, iband)
                           DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                              j = lm1 + i
                              cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                           END DO  !i -> loop over mixed basis functions
                        END DO  !ibando
                     END DO  !iband

                     ! go to lm start index for next m-quantum number
                     lm = lm + mpbasis%num_rad_bas_fun(l, itype)

                     !
                     ! case: m>0
                     !

                     rarr1 = 0.0
                     monepm = 1
                     DO m = 1, l
                        monepm = -monepm
                        moneplm = monepl*monepm

                        ! calculate the contributions which are identical for m>0 and m <0
                        rarr1 = 0.0
                        IF (abs(m) <= l2) THEN
                           lmp1 = lp1 + l1*hybrid%nindx(l1, itype)
                           IF (mod(l1, 2) == 0) THEN
                              lmp2 = lp2 + (m + l2)*hybrid%nindx(l2, itype)
                           ELSE
                              lmp2 = lp2 + (-m + l2)*hybrid%nindx(l2, itype)
                           END IF

                           rdum = hybdat%gauntarr(1, l1, l2, l, 0, m)
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)
                                 IF (mod(l1, 2) /= 0) rdum1 = moneplm*rdum1
                                 DO ibando = 1, mnobd
                                    rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, -m, m)
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l1, 2) == 0) rdum1 = moneplm*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp1, iatom1)
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
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 rdum1 = rdum*cmt_nk(iband, lmp3, iatom1)
                                 IF (mod(l2, 2) == 0) rdum1 = moneplm*rdum1
                                 DO ibando = 1, mnobd
                                    rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                 END DO  ! ibando
                              END DO  ! iband
                           END IF  ! rdum .ne. 0

                           IF (offdiag) THEN
                              rdum = hybdat%gauntarr(1, l2, l1, l, 0, m)
                              IF (rdum /= 0) THEN
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)
                                    IF (mod(l2, 2) /= 0) rdum1 = moneplm*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
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
                              IF (rdum /= 0) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = -moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))/2
                                 END IF

                                 rdum = -moneplm*monepl1m1*rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!-moneplm*monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband

                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, -m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                    ELSE
                                       lmp3 = lmp1
                                       fac = 1/2.*monepl1m1*(sign(1, m2) - sign(1, m1))
                                    END IF
                                    rdum = monepl1m1*moneplm*rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!monepl1m1*moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sr2
                                       IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                       DO ibando = 1, mnobd
                                          rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband

                                 END IF  ! rdum

                              END IF  ! offdiag
                           END IF  ! abs(m2) .le. l2 .and. m2 .ne. 0

                           m2 = m1 - m
                           IF (abs(m2) <= l2 .and. m2 /= 0) THEN

                              rdum = hybdat%gauntarr(1, l1, l2, l, m1, -m)
                              IF (rdum /= 0) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m1) - sign(1, m2))
                                 END IF
                                 rdum = monepl1m1*rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = rdum1*fac
                                    DO ibando = 1, mnobd
                                       rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                    ELSE
                                       lmp3 = lmp1
                                       fac = -monepl1m1*(sign(1, m1) - sign(1, m2))/2
                                    END IF
                                    rdum = -monepl1m1*rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!-monepl1m1*rdum*cmt_nk(iband,lmp2,iatom1)/sr2
                                       IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                       DO ibando = 1, mnobd
                                          rarr1(ibando, iband) = rarr1(ibando, iband) + rdum1*cmt(ibando, lmp3, iatom1)
                                       END DO  ! ibando
                                    END DO  ! iband

                                 END IF  ! rdum .ne. 0
                              END IF  ! offdiag

                           END IF  !  abs(m2) .le. l2 .and. m2 .ne. 0

                           !go to lmp start index for next m1-quantum number
                           lmp1 = lmp1 + hybrid%nindx(l1, itype)
                        END DO  ! m1

                        ! multiply rarr1 by (-1)**(l+m+1)
                        rarr1(:, :) = (-1)*moneplm*rarr1(:, :)

                        ! go to lm mixed basis startindx for l and m
                        lm1 = lm + (iatom1 - 1 - iiatom)*ioffset

                        DO iband = bandi, bandf
                           DO ibando = 1, mnobd
                              rdum = rarr1(ibando, iband)
                              DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                                 j = lm1 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpbasis%num_rad_bas_fun(l, itype)

                     END DO  ! m=1,l

                  END DO !n
                  lm_0 = lm_0 + mpbasis%num_rad_bas_fun(l, itype)*(2*l + 1) ! go to the m start index of the next l-quantum number
                  IF (lm /= lm_0) STOP 'wavefproducts: counting of lm-index incorrect (bug?)'
               END DO !l

            END IF  ! iatom1 .ne. iatom2

            lm_0 = lm_00
         END DO !ieq
         iiatom = iiatom + atoms%neq(itype)
         lm_00 = lm_00 + atoms%neq(itype)*ioffset
      END DO  !itype

      ic = nbasm_mt
      DO igpt = 1, mpbasis%ngptm(iq)
         ic = ic + 1
         DO ibando = 1, mnobd
            DO iband = bandi, bandf
               cprod(ic, ibando, iband) = cprod_ir(iband, ibando, igpt)
            END DO
         END DO
      END DO

      CALL timestop("wavefproducts_inv")

   END SUBROUTINE wavefproducts_inv

   SUBROUTINE wavefproducts_inv5(&
  &                    bandi, bandf, bandoi, bandof,&
  &                    dimension, input, jsp, atoms,&
  &                    lapw, kpts,&
  &                    nk, iq, hybdat, mnobd, hybrid,&
  &                    parent, cell,&
  &                    nbasm_mt, sym,&
  &                    noco,&
  &                    nkqpt, cprod)

      USE m_util, ONLY: modulo1
      USE m_olap, ONLY: gptnorm
      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_io_hybrid
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
      INTEGER, INTENT(IN)      :: mnobd
      INTEGER, INTENT(IN)      :: nbasm_mt
      INTEGER, INTENT(OUT)     :: nkqpt

      ! - arrays -
      INTEGER, INTENT(IN)      ::    parent(kpts%nkptf)

      REAL, INTENT(OUT)        ::    cprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)

      ! - local scalars -
      INTEGER                 ::    i, ikpt, ic, iband, iband1, igpt, igptp, ig, ig2, ig1
      INTEGER                 ::    iatom, iiatom, itype, ieq, ishift
      INTEGER                 ::    ibando, iatom1, iatom2, ioffset
      INTEGER                 ::    k, l, p, l1, m1, l2, m2, p1, p2, n, ok
      INTEGER                 ::    igptm, iigptm
      INTEGER                 ::    lm, lm1, lm2, lm_0, lm_00, lm1_0, lm2_0, lmp1, lmp2, lmp3, lmp4, lp1, lp2
      INTEGER                 ::    j, ll, m
      INTEGER                 ::    nbasm_ir, ngpt0
      INTEGER                 ::    n1, n2, nbasfcn
      REAL                    ::    svol, sr2
      REAL                    ::    rdum, rfac, rfac1, rfac2, rdum1, rdum2
      REAL                    ::    sin1, sin2, cos1, cos2, add1, add2
      REAL                    ::    fac, fac1, fac2
      REAL                    ::    monepl1, monepl2, monepl, monepm1, monepm, moneplm, monepl1m1
      COMPLEX                 ::    fexp
      COMPLEX                 ::    cdum, cconst, cfac
      COMPLEX, PARAMETER       ::    img = (0.0, 1.0)
      LOGICAL                 ::    offdiag
      TYPE(t_lapw)            ::    lapw_nkqpt

      ! - local arrays -
      INTEGER                 ::    g(3), g_t(3)
      INTEGER                 ::    lmstart(0:atoms%lmaxd, atoms%ntype)
      INTEGER, ALLOCATABLE    ::    gpt0(:, :)
      INTEGER, ALLOCATABLE    ::    pointer(:, :, :)

      REAL                    ::    kqpt(3), kqpthlp(3)
      REAL                    ::    bkpt(3)
      REAL                    ::    cmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      REAL                    ::    cmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      REAL                    ::    rarr1(bandoi:bandof)
      REAL                    ::    rarr2(bandoi:bandof, bandf - bandi + 1)
      REAL                    ::    rarr3(2, bandoi:bandof, bandf - bandi + 1)
      REAL, ALLOCATABLE       ::    z0(:, :)

      COMPLEX                 ::    cexp(atoms%nat), cexp_nk(atoms%nat)
      COMPLEX, ALLOCATABLE     ::    ccmt_nk(:, :, :)
      COMPLEX, ALLOCATABLE     ::    ccmt(:, :, :)
      TYPE(t_mat)             :: z_nk, z_kqpt

      CALL timestart("wavefproducts_inv5")
      CALL timestart("wavefproducts_inv5 IR")

      cprod = 0
      svol = sqrt(cell%omtil)
      sr2 = sqrt(2.0)

      nbasm_ir = maxval(mpbasis%ngptm)

      !
      ! compute k+q point for q (iq) in EIBZ(k)
      !

      kqpthlp = kpts%bkf(:, nk) + kpts%bkf(:, iq)
      ! kqpt can lie outside the first BZ, transfer it back
      kqpt = modulo1(kqpthlp, kpts%nkpt3)
      g_t(:) = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      nkqpt = 0
      DO ikpt = 1, kpts%nkptf
         IF (maxval(abs(kqpt - kpts%bkf(:, ikpt))) <= 1E-06) THEN
            nkqpt = ikpt
            EXIT
         END IF
      END DO
      IF (nkqpt == 0) STOP 'wavefproducts_inv5: k-point not found'

      !
      ! compute G's fulfilling |bk(:,nkqpt) + G| <= rkmax
      !
      CALL lapw_nkqpt%init(input, noco, kpts, atoms, sym, nkqpt, cell, sym%zrfs)

      nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
      call z_nk%alloc(.true., nbasfcn, dimension%neigd)
      nbasfcn = MERGE(lapw_nkqpt%nv(1) + lapw_nkqpt%nv(2) + 2*atoms%nlotot, lapw_nkqpt%nv(1) + atoms%nlotot, noco%l_noco)
      call z_kqpt%alloc(.true., nbasfcn, dimension%neigd)

      ! read in z at k-point nk and nkqpt
      call timestart("read_z")
      CALL read_z(z_nk, kpts%nkptf*(jsp - 1) + nk)
      call read_z(z_kqpt, kpts%nkptf*(jsp - 1) + nkqpt)
      call timestop("read_z")

      g(1) = maxval(abs(lapw%k1(:lapw%nv(jsp), jsp))) &
     &     + maxval(abs(lapw_nkqpt%k1(:lapw_nkqpt%nv(jsp), jsp)))&
     &     + maxval(abs(mpbasis%gptm(1, mpbasis%gptm_ptr(:mpbasis%ngptm(iq), iq)))) + 1
      g(2) = maxval(abs(lapw%k2(:lapw%nv(jsp), jsp)))&
     &     + maxval(abs(lapw_nkqpt%k2(:lapw_nkqpt%nv(jsp), jsp)))&
     &     + maxval(abs(mpbasis%gptm(2, mpbasis%gptm_ptr(:mpbasis%ngptm(iq), iq)))) + 1
      g(3) = maxval(abs(lapw%k3(:lapw%nv(jsp), jsp)))&
     &     + maxval(abs(lapw_nkqpt%k3(:lapw_nkqpt%nv(jsp), jsp)))&
     &     + maxval(abs(mpbasis%gptm(3, mpbasis%gptm_ptr(:mpbasis%ngptm(iq), iq)))) + 1

      ALLOCATE (pointer(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts_inv5: error allocation pointer'
      ALLOCATE (gpt0(3, size(pointer)), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts_inv5: error allocation gpt0'

      if (.not. allocated(hybdat%stepfunc_r)) then
         call timestart("setup stepfunction")
         ALLOCATE (hybdat%stepfunc_r(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
         IF (ok /= 0) then
            call juDFT_error('wavefproducts_inv5: error allocation stepfunc_r')
         endif

         DO i = -g(1), g(1)
            DO j = -g(2), g(2)
               DO k = -g(3), g(3)
                  hybdat%stepfunc_r(i, j, k) = stepfunction(cell, atoms, (/i, j, k/))
               END DO
            END DO
         END DO
         call timestop("setup stepfunction")
      endif

      !
      ! convolute phi(n,k) with the step function and store in cpw0
      !

      !(1) prepare list of G vectors
      call timestart("prep list of Gvec")
      pointer = 0
      ic = 0
      DO ig1 = 1, lapw%nv(jsp)
         DO igptm = 1, mpbasis%ngptm(iq)
            iigptm = mpbasis%gptm_ptr(igptm, iq)
            g(1) = lapw%k1(ig1, jsp) + mpbasis%gptm(1, iigptm) - g_t(1)
            g(2) = lapw%k2(ig1, jsp) + mpbasis%gptm(2, iigptm) - g_t(2)
            g(3) = lapw%k3(ig1, jsp) + mpbasis%gptm(3, iigptm) - g_t(3)
            IF (pointer(g(1), g(2), g(3)) == 0) THEN
               ic = ic + 1
               gpt0(:, ic) = g
               pointer(g(1), g(2), g(3)) = ic
            END IF
         END DO
      END DO
      ngpt0 = ic
      call timestop("prep list of Gvec")

      !(2) calculate convolution
      call timestart("calc convolution")
      ALLOCATE (z0(bandoi:bandof, ngpt0), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts_inv5: error allocation z0'
      z0 = 0
      call timestart("step function")
      DO ig2 = 1, lapw_nkqpt%nv(jsp)
         rarr1 = z_kqpt%data_r(ig2, bandoi:bandof)
         DO ig = 1, ngpt0
            g(1) = gpt0(1, ig) - lapw_nkqpt%k1(ig2, jsp)
            g(2) = gpt0(2, ig) - lapw_nkqpt%k2(ig2, jsp)
            g(3) = gpt0(3, ig) - lapw_nkqpt%k3(ig2, jsp)
            rdum = hybdat%stepfunc_r(g(1), g(2), g(3))/svol
            DO n2 = bandoi, bandof
               z0(n2, ig) = z0(n2, ig) + rarr1(n2)*rdum
            END DO
         END DO
      END DO
      call timestop("step function")

      call timestart("hybrid gptm")
      ic = nbasm_mt
      DO igptm = 1, mpbasis%ngptm(iq)
         rarr2 = 0
         ic = ic + 1
         iigptm = mpbasis%gptm_ptr(igptm, iq)

         DO ig1 = 1, lapw%nv(jsp)
            g(1) = lapw%k1(ig1, jsp) + mpbasis%gptm(1, iigptm) - g_t(1)
            g(2) = lapw%k2(ig1, jsp) + mpbasis%gptm(2, iigptm) - g_t(2)
            g(3) = lapw%k3(ig1, jsp) + mpbasis%gptm(3, iigptm) - g_t(3)

            ig2 = pointer(g(1), g(2), g(3))

            IF (ig2 == 0) THEN
               STOP 'wavefproducts_inv5: pointer undefined'
            END IF

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

      WRITE (2005, *) 'Point B'
      DO n2 = 1, 1
         DO n1 = 1, 2
            DO ic = 1, 20
               WRITE (2010, '(3i7,f15.8)') ic, n1, n2, cprod(ic, n1, n2)
            END DO
         END DO
      END DO

      DEALLOCATE (z0, pointer, gpt0)
      CALL timestop("wavefproducts_inv5 IR")

      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            lmstart(l, itype) = sum((/(hybrid%nindx(ll, itype)*(2*ll + 1), ll=0, l - 1)/))
         END DO
      END DO

      ! read in cmt coefficient at k-point nk
      ALLOCATE (ccmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat), ccmt(dimension%neigd, hybrid%maxlmindx, atoms%nat), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts_inv5: error allocation ccmt_nk/ccmt'

      call read_cmt(ccmt_nk, nk)
      !read in cmt coefficients at k+q point
      call read_cmt(ccmt, nkqpt)

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            cexp(iatom) = exp((-img)*tpi_const*dotprod(kpts%bkf(:, iq) + kpts%bkf(:, nk), atoms%taual(:, iatom)))

            cexp_nk(iatom) = exp((-img)*tpi_const*dotprod(kpts%bkf(:, nk), atoms%taual(:, iatom)))
         END DO
      END DO

      rfac = 1./sr2
      cfac = -img/sr2
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
                           cmt(:, lm1, iatom) = (ccmt(:, lm1, iatom) + rdum*ccmt(:, lm2, iiatom))*cexp(iatom)*rfac

                           cmt_nk(:, lm1, iatom) = (ccmt_nk(:, lm1, iatom) + rdum*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*rfac
                        ELSE IF (m > 0) THEN

                           cmt(:, lm1, iatom) = (ccmt(:, lm1, iatom) - rdum*ccmt(:, lm2, iiatom))*cexp(iatom)*cfac

                           cmt_nk(:, lm1, iatom) = (ccmt_nk(:, lm1, iatom) - rdum*ccmt_nk(:, lm2, iiatom))*cexp_nk(iatom)*cfac
                        ELSE
                           IF (mod(l, 2) == 0) THEN
                              cmt(:, lm1, iatom) = ccmt(:, lm1, iatom)*cexp(iatom)
                              cmt_nk(:, lm1, iatom) = ccmt_nk(:, lm1, iatom)*cexp_nk(iatom)
                           ELSE
                              cmt(:, lm1, iatom) = ccmt(:, lm1, iatom)*(-img)*cexp(iatom)
                              cmt_nk(:, lm1, iatom) = ccmt_nk(:, lm1, iatom)*(-img)*cexp_nk(iatom)
                           END IF
                        END IF
                     ELSE
                        cdum = rdum*cexp(iatom)*cexp(iiatom)
                        cmt(:, lm1, iatom) = (ccmt(:, lm1, iatom) + cdum*ccmt(:, lm2, iiatom))*rfac

                        cmt(:, lm1, iiatom) = (ccmt(:, lm1, iatom) - cdum*ccmt(:, lm2, iiatom))*cfac

                        cdum = rdum*cexp_nk(iatom)*cexp_nk(iiatom)
                        cmt_nk(:, lm1, iatom) = (ccmt_nk(:, lm1, iatom) + cdum*ccmt_nk(:, lm2, iiatom))*rfac

                        cmt_nk(:, lm1, iiatom) = (ccmt_nk(:, lm1, iatom) - cdum*ccmt_nk(:, lm2, iiatom))*cfac
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
         ioffset = sum((/((2*ll + 1)*mpbasis%num_rad_bas_fun(ll, itype), ll=0, hybrid%lcutm1(itype))/))
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
                              IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
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

                        ishift = -2*m*mpbasis%num_rad_bas_fun(l, itype)

                        ! go to lm mixed basis startindx for l and m
                        lm1 = lm + (iatom1 - 1 - iiatom)*ioffset
                        lm2 = lm + (iatom2 - 1 - iiatom)*ioffset + ishift

                        rdum = tpi_const*dotprod(kpts%bkf(:, iq), atoms%taual(:, iatom1))
                        rfac1 = sin(rdum)/sr2
                        rfac2 = cos(rdum)/sr2
                        DO iband = bandi, bandf
                           DO ibando = bandoi, bandof
                              rdum1 = rarr3(1, ibando, iband)
                              rdum2 = rarr3(2, ibando, iband)
!                       sin1  = rdum1*rfac1
!                       cos1  = rdum1*rfac2
!                       sin2  = rdum2*rfac1
!                       cos2  = rdum2*rfac2
                              add1 = rdum1*rfac2 + rdum2*rfac1
                              add2 = rdum2*rfac2 - rdum1*rfac1
                              DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                                 j = lm1 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*add1!( cos1 + sin2 )
                                 j = lm2 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*add2!( cos2 - sin1 )

                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpbasis%num_rad_bas_fun(l, itype)

                     END DO  !m

                  END DO !n
                  lm_0 = lm_0 + mpbasis%num_rad_bas_fun(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
                  IF (lm /= lm_0) STOP 'wavefproducts_inv5: counting of lm-index incorrect (bug?)'
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
                           IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
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
                           IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                 END IF
                                 rdum = rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = bandoi, bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(1, l2, l1, l, m2, -m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1
                                    ELSE
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                       fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                    END IF
                                    rdum = moneplm*rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sr2
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
                              IF (rdum /= 0) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))
                                 END IF
                                 rdum = moneplm*rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!moneplm*rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = bandoi, bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband

                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(1, l2, l1, l, m2, m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m1) + sign(1, m2) /= 0) THEN
                                       lmp3 = lmp1
                                    ELSE
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                       fac = 1/2.*monepl1m1*(sign(1, m1) - sign(1, m2))
                                    END IF
                                    rdum = rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!rdum*cmt_nk(iband,lmp2,iatom1)/sr2
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
                              DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                                 j = lm1 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpbasis%num_rad_bas_fun(l, itype)

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
                           IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
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
                           DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                              j = lm1 + i
                              cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                           END DO  !i -> loop over mixed basis functions
                        END DO  !ibando
                     END DO  !iband

                     ! go to lm start index for next m-quantum number
                     lm = lm + mpbasis%num_rad_bas_fun(l, itype)

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
                           IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
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
                           IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN
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
                              IF (rdum /= 0) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = -moneplm*monepl1m1*(sign(1, m2) - sign(1, m1))/2
                                 END IF

                                 rdum = -moneplm*monepl1m1*rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!-moneplm*monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = fac*rdum1
                                    DO ibando = bandoi, bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband

                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, -m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                    ELSE
                                       lmp3 = lmp1
                                       fac = 1/2.*monepl1m1*(sign(1, m2) - sign(1, m1))
                                    END IF
                                    rdum = monepl1m1*moneplm*rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!monepl1m1*moneplm*rdum*cmt_nk(iband,lmp2,iatom1)/sr2
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
                              IF (rdum /= 0) THEN

                                 IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                    lmp2 = lp2 + (-m2 + l2)*hybrid%nindx(l2, itype)
                                 ELSE
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    fac = 1/2.*moneplm*monepl1m1*(sign(1, m1) - sign(1, m2))
                                 END IF
                                 rdum = monepl1m1*rdum/sr2
                                 DO iband = bandi, bandf
                                    rdum1 = rdum*cmt_nk(iband, lmp1, iatom1)!monepl1m1*rdum*cmt_nk(iband,lmp1,iatom1)/sr2
                                    IF (sign(1, m2) + sign(1, m1) == 0) rdum1 = rdum1*fac
                                    DO ibando = bandoi, bandof
                                       rarr2(ibando, iband) = rarr2(ibando, iband) + rdum1*cmt(ibando, lmp2, iatom1)
                                    END DO  ! ibando
                                 END DO  ! iband
                              END IF  ! rdum .ne. 0

                              IF (offdiag) THEN
                                 rdum = hybdat%gauntarr(2, l1, l2, l, m1, m)
                                 IF (rdum /= 0) THEN
                                    lmp2 = lp2 + (m2 + l2)*hybrid%nindx(l2, itype)
                                    IF (sign(1, m2) + sign(1, m1) /= 0) THEN
                                       lmp3 = lmp1 - 2*m1*hybrid%nindx(l1, itype)
                                    ELSE
                                       lmp3 = lmp1
                                       fac = -monepl1m1*(sign(1, m1) - sign(1, m2))/2
                                    END IF
                                    rdum = -monepl1m1*rdum/sr2
                                    DO iband = bandi, bandf
                                       rdum1 = rdum*cmt_nk(iband, lmp2, iatom1)!-monepl1m1*rdum*cmt_nk(iband,lmp2,iatom1)/sr2
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
                              DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                                 j = lm1 + i
                                 cprod(j, ibando, iband) = cprod(j, ibando, iband) + hybdat%prodm(i, n, l, itype)*rdum
                              END DO  !i -> loop over mixed basis functions
                           END DO  !ibando
                        END DO  !iband

                        ! go to lm start index for next m-quantum number
                        lm = lm + mpbasis%num_rad_bas_fun(l, itype)

                     END DO  ! m=1,l

                  END DO !n
                  lm_0 = lm_0 + mpbasis%num_rad_bas_fun(l, itype)*(2*l + 1) ! go to the m start index of the next l-quantum number
                  IF (lm /= lm_0) STOP 'wavefproducts_inv5: counting of lm-index incorrect (bug?)'
               END DO !l

               call timestop("iatom1 eq iatom2")
            END IF  ! iatom1 .ne. iatom2

            lm_0 = lm_00
         END DO !ieq
         iiatom = iiatom + atoms%neq(itype)
         lm_00 = lm_00 + atoms%neq(itype)*ioffset
      END DO  !itype
      CALL timestop("wavefproducts_inv5")

   END SUBROUTINE wavefproducts_inv5

   SUBROUTINE wavefproducts_noinv5(&
  &                      bandi, bandf, bandoi, bandof,&
  &                      nk, iq, dimension, input, jsp,&
  &                      cell, atoms, hybrid,&
  &                      hybdat,&
  &                      kpts,&
  &                      mnobd,&
  &                      lapw, sym,&
  &                      nbasm_mt,&
  &                      noco,&
  &                      nkqpt, cprod)

      USE m_constants
      USE m_util, ONLY: modulo1
      USE m_olap, ONLY: gptnorm
      USE m_trafo
      USE m_wrapper
      USE m_types
      USE m_io_hybrid
      IMPLICIT NONE
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_noco), INTENT(IN)        :: noco
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_hybrid), INTENT(IN)      :: hybrid
      TYPE(t_hybdat), INTENT(INOUT)   :: hybdat

!     - scalars -
      INTEGER, INTENT(IN)      ::  bandi, bandf, bandoi, bandof
      INTEGER, INTENT(IN)      ::  nk, iq, jsp
      INTEGER, INTENT(IN)      :: mnobd, nbasm_mt
      INTEGER, INTENT(OUT)     ::  nkqpt

!     - arrays -

      COMPLEX, INTENT(OUT)    ::  cprod(hybrid%maxbasm1, bandoi:bandof, bandf - bandi + 1)

!     - local scalars -
      INTEGER                 ::  ic, l, n, l1, l2, n1, n2, lm_0, lm1_0, lm2_0, lm, lm1, lm2, m1, m2, i, j, ll
      INTEGER                 ::  itype, ieq, ikpt, ikpt1, ikpt2, igpt, igptp, igpt1, igpt2, iband, iband1, iband2
      INTEGER                 ::  k, ic1, ioffset, ibando, ig1, ig2, ig
      INTEGER                 ::  igptm, iigptm
      INTEGER                 ::  q, idum
      INTEGER                 :: nbasm_ir, ngpt0
      INTEGER                 ::  nbasmmt, nbasfcn
      INTEGER                 ::  ok, m

      REAL                    ::  rdum, svol, s2

      COMPLEX                 ::  cdum, cdum0, cdum1
      COMPLEX                 ::  cexp
      COMPLEX, PARAMETER       ::  img = (0.0, 1.0)

      LOGICAL                 ::  offdiag
      TYPE(t_lapw)            ::    lapw_nkqpt

!      - local arrays -
      INTEGER                 ::  g(3), g_t(3)
      INTEGER                 ::  lmstart(0:atoms%lmaxd, atoms%ntype)
      INTEGER, ALLOCATABLE    ::  gpt0(:, :)
      INTEGER, ALLOCATABLE    ::  pointer(:, :, :)

      REAL                    ::  bkpt(3)
      REAL                    ::  kqpt(3), kqpthlp(3)

      COMPLEX                 ::  carr1(bandoi:bandof)
      COMPLEX                 ::  carr2(bandoi:bandof, bandf - bandi + 1)
      TYPE(t_mat)             ::  z_nk, z_kqpt
      COMPLEX                 ::  cmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      COMPLEX                 ::  cmt_nk(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      COMPLEX, ALLOCATABLE     ::  z0(:, :)

      call timestart("wavefproducts_noinv5")
      call timestart("wavefproducts_noinv5 IR")
      cprod = 0
      svol = sqrt(cell%omtil)
      s2 = sqrt(2.0)

      nbasm_ir = maxval(mpbasis%ngptm)

      !
      ! compute k+q point for given q point in EIBZ(k)
      !
      kqpthlp = kpts%bkf(:, nk) + kpts%bkf(:, iq)
      ! k+q can lie outside the first BZ, transfer
      ! it back into the 1. BZ
      kqpt = modulo1(kqpthlp, kpts%nkpt3)
      g_t(:) = nint(kqpt - kqpthlp)
      ! determine number of kqpt
      nkqpt = 0
      DO ikpt = 1, kpts%nkptf
         IF (maxval(abs(kqpt - kpts%bkf(:, ikpt))) <= 1E-06) THEN
            nkqpt = ikpt
            EXIT
         END IF
      END DO
      IF (nkqpt == 0) STOP 'wavefproducts: k-point not found'

      !
      ! compute G's fulfilling |bk(:,nkqpt) + G| <= rkmax
      !
      CALL lapw_nkqpt%init(input, noco, kpts, atoms, sym, nkqpt, cell, sym%zrfs)
      nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
      call z_nk%alloc(.false., nbasfcn, dimension%neigd)
      nbasfcn = MERGE(lapw_nkqpt%nv(1) + lapw_nkqpt%nv(2) + 2*atoms%nlotot, lapw_nkqpt%nv(1) + atoms%nlotot, noco%l_noco)
      call z_kqpt%alloc(.false., nbasfcn, dimension%neigd)

      ! read in z at k-point nk and nkqpt
      call timestart("read_z")
      call read_z(z_nk, kpts%nkptf*(jsp - 1) + nk)
      call read_z(z_kqpt, kpts%nkptf*(jsp - 1) + nkqpt)
      call timestop("read_z")

      g(1) = maxval(abs(lapw%k1(:lapw%nv(jsp), jsp))) &
     &     + maxval(abs(lapw_nkqpt%k1(:lapw_nkqpt%nv(jsp), jsp)))&
     &     + maxval(abs(mpbasis%gptm(1, mpbasis%gptm_ptr(:mpbasis%ngptm(iq), iq)))) + 1
      g(2) = maxval(abs(lapw%k2(:lapw%nv(jsp), jsp)))&
     &     + maxval(abs(lapw_nkqpt%k2(:lapw_nkqpt%nv(jsp), jsp)))&
     &     + maxval(abs(mpbasis%gptm(2, mpbasis%gptm_ptr(:mpbasis%ngptm(iq), iq)))) + 1
      g(3) = maxval(abs(lapw%k3(:lapw%nv(jsp), jsp)))&
     &     + maxval(abs(lapw_nkqpt%k3(:lapw_nkqpt%nv(jsp), jsp)))&
     &     + maxval(abs(mpbasis%gptm(3, mpbasis%gptm_ptr(:mpbasis%ngptm(iq), iq)))) + 1

      ALLOCATE (pointer(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts_noinv2: error allocation pointer'
      ALLOCATE (gpt0(3, size(pointer)), stat=ok)
      IF (ok /= 0) STOP 'wavefproducts_noinv2: error allocation gpt0'

      if (.not. allocated(hybdat%stepfunc_c)) then
         call timestart("setup stepfunc")
         ALLOCATE (hybdat%stepfunc_c(-g(1):g(1), -g(2):g(2), -g(3):g(3)), stat=ok)
         IF (ok /= 0) then
            call juDFT_error('wavefproducts_noinv2: error allocation stepfunc')
         endif
         DO i = -g(1), g(1)
            DO j = -g(2), g(2)
               DO k = -g(3), g(3)
                  hybdat%stepfunc_c(i, j, k) = stepfunction(cell, atoms, (/i, j, k/))
               END DO
            END DO
         END DO
         call timestop("setup stepfunc")
      endif

      !
      ! convolute phi(n,k) with the step function and store in cpw0
      !

      !(1) prepare list of G vectors
      call timestart("prep list of Gvec")
      pointer = 0
      ic = 0
      DO ig1 = 1, lapw%nv(jsp)
         DO igptm = 1, mpbasis%ngptm(iq)
            iigptm = mpbasis%gptm_ptr(igptm, iq)
            g(1) = lapw%k1(ig1, jsp) + mpbasis%gptm(1, iigptm) - g_t(1)
            g(2) = lapw%k2(ig1, jsp) + mpbasis%gptm(2, iigptm) - g_t(2)
            g(3) = lapw%k3(ig1, jsp) + mpbasis%gptm(3, iigptm) - g_t(3)
            IF (pointer(g(1), g(2), g(3)) == 0) THEN
               ic = ic + 1
               gpt0(:, ic) = g
               pointer(g(1), g(2), g(3)) = ic
            END IF
         END DO
      END DO
      ngpt0 = ic
      call timestop("prep list of Gvec")

      !(2) calculate convolution
      call timestart("calc convolution")
      call timestart("step function")
      ALLOCATE (z0(bandoi:bandof, ngpt0))
      z0 = 0
      DO ig2 = 1, lapw_nkqpt%nv(jsp)
         carr1 = z_kqpt%data_c(ig2, bandoi:bandof)
         DO ig = 1, ngpt0
            g(1) = gpt0(1, ig) - lapw_nkqpt%k1(ig2, jsp)
            g(2) = gpt0(2, ig) - lapw_nkqpt%k2(ig2, jsp)
            g(3) = gpt0(3, ig) - lapw_nkqpt%k3(ig2, jsp)
            cdum = hybdat%stepfunc_c(g(1), g(2), g(3))/svol
            DO n2 = bandoi, bandof
               z0(n2, ig) = z0(n2, ig) + carr1(n2)*cdum
            END DO
         END DO
      END DO
      call timestop("step function")

      call timestart("hybrid gptm")
      ic = nbasm_mt
      DO igptm = 1, mpbasis%ngptm(iq)
         carr2 = 0
         ic = ic + 1
         iigptm = mpbasis%gptm_ptr(igptm, iq)

         DO ig1 = 1, lapw%nv(jsp)
            g(1) = lapw%k1(ig1, jsp) + mpbasis%gptm(1, iigptm) - g_t(1)
            g(2) = lapw%k2(ig1, jsp) + mpbasis%gptm(2, iigptm) - g_t(2)
            g(3) = lapw%k3(ig1, jsp) + mpbasis%gptm(3, iigptm) - g_t(3)

            ig2 = pointer(g(1), g(2), g(3))

            IF (ig2 == 0) THEN
               STOP 'wavefproducts_noinv2: pointer undefined'
            END IF

            DO n1 = 1, bandf - bandi + 1
               cdum1 = conjg(z_nk%data_c(ig1, n1))
               DO n2 = bandoi, bandof
                  carr2(n2, n1) = carr2(n2, n1) + cdum1*z0(n2, ig2)
               END DO
            END DO

         END DO
         cprod(ic, :, :) = carr2(:, :)
      END DO
      call timestop("hybrid gptm")
      DEALLOCATE (z0, pointer, gpt0)
      call timestop("calc convolution")

      call timestop("wavefproducts_noinv5 IR")

!       RETURN

      !
      ! MT contribution
      !

      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            lmstart(l, itype) = sum((/(hybrid%nindx(ll, itype)*(2*ll + 1), ll=0, l - 1)/))
         END DO
      END DO

      ! read in cmt coefficients from direct access file cmt
      call read_cmt(cmt_nk(:, :, :), nk)
      call read_cmt(cmt(:, :, :), nkqpt)

      lm_0 = 0
      ic = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            ic1 = 0

            cexp = exp(-img*tpi_const*dot_product(kpts%bkf(:, iq), atoms%taual(:, ic)))

            DO l = 0, hybrid%lcutm1(itype)
               DO n = 1, hybdat%nindxp1(l, itype) ! loop over basis-function products

                  l1 = hybdat%prod(n, l, itype)%l1 !
                  l2 = hybdat%prod(n, l, itype)%l2 ! current basis-function product
                  n1 = hybdat%prod(n, l, itype)%n1 ! = bas(:,n1,l1,itype)*bas(:,n2,l2,itype) = b1*b2
                  n2 = hybdat%prod(n, l, itype)%n2 !

                  IF (mod(l1 + l2 + l, 2) /= 0) cycle

                  offdiag = l1 /= l2 .or. n1 /= n2 ! offdiag=true means that b1*b2 and b2*b1 are different combinations
                  !(leading to the same basis-function product)

                  lm1_0 = lmstart(l1, itype) ! start at correct lm index of cmt-coefficients
                  lm2_0 = lmstart(l2, itype) ! (corresponding to l1 and l2)

                  lm = lm_0
                  DO m = -l, l

                     carr2 = 0.0

                     lm1 = lm1_0 + n1 ! go to lm index for m1=-l1
                     DO m1 = -l1, l1
                        m2 = m1 + m ! Gaunt condition -m1+m2-m=0
                        IF (abs(m2) <= l2) THEN
                           lm2 = lm2_0 + n2 + (m2 + l2)*hybrid%nindx(l2, itype)
                           rdum = hybdat%gauntarr(1, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 cdum = rdum*conjg(cmt_nk(iband, lm1, ic)) !nk
                                 DO iband1 = bandoi, bandof
                                    carr2(iband1, iband) = carr2(iband1, iband) + cdum*cmt(iband1, lm2, ic) !ikpt

                                 END DO
                              END DO
                           END IF
                        END IF

                        m2 = m1 - m ! switch role of b1 and b2
                        IF (abs(m2) <= l2 .and. offdiag) THEN
                           lm2 = lm2_0 + n2 + (m2 + l2)*hybrid%nindx(l2, itype)
                           rdum = hybdat%gauntarr(2, l1, l2, l, m1, m) ! precalculated Gaunt coefficient
                           IF (rdum /= 0) THEN
                              DO iband = bandi, bandf
                                 cdum = rdum*conjg(cmt_nk(iband, lm2, ic)) !nk
                                 DO iband1 = bandoi, bandof
                                    carr2(iband1, iband) = carr2(iband1, iband) + cdum*cmt(iband1, lm1, ic)
                                 END DO
                              END DO
                           END IF
                        END IF

                        lm1 = lm1 + hybrid%nindx(l1, itype) ! go to lm start index for next m1-quantum number

                     END DO  !m1

                     DO iband = bandi, bandf
                        DO iband1 = bandoi, bandof
                           cdum = carr2(iband1, iband)*cexp
                           DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                              j = lm + i
                              cprod(j, iband1, iband) = cprod(j, iband1, iband) + hybdat%prodm(i, n, l, itype)*cdum
                           END DO

                        END DO
                     END DO

                     lm = lm + mpbasis%num_rad_bas_fun(l, itype) ! go to lm start index for next m-quantum number

                  END DO

               END DO
               lm_0 = lm_0 + mpbasis%num_rad_bas_fun(l, itype)*(2*l + 1) ! go to the lm start index of the next l-quantum number
               IF (lm /= lm_0) STOP 'wavefproducts_noinv2: counting of lm-index incorrect (bug?)'
            END DO
         END DO
      END DO

      call timestop("wavefproducts_noinv5")

   END SUBROUTINE wavefproducts_noinv5

   !private subroutine
   FUNCTION stepfunction(cell, atoms, g)
      USE m_types
      USE m_constants
      USE m_olap
      IMPLICIT NONE

      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN) :: g(3)
      COMPLEX             :: stepfunction  !Is real in inversion case
      REAL                :: gnorm, gnorm3, r, fgr
      INTEGER             :: itype, ieq, icent

      gnorm = gptnorm(g, cell%bmat)
      gnorm3 = gnorm**3
      IF (gnorm == 0) THEN
         stepfunction = 1
         DO itype = 1, atoms%ntype
            stepfunction = stepfunction - atoms%neq(itype)*atoms%volmts(itype)/cell%omtil
         END DO
      ELSE
         stepfunction = 0
         icent = 0
         DO itype = 1, atoms%ntype
            r = gnorm*atoms%rmt(itype)
            fgr = fpi_const*(sin(r) - r*cos(r))/gnorm3/cell%omtil
            DO ieq = 1, atoms%neq(itype)
               icent = icent + 1
               stepfunction = stepfunction - fgr*exp(-cmplx(0., tpi_const*dot_product(atoms%taual(:, icent), g)))
            ENDDO
         ENDDO
      ENDIF

   END FUNCTION stepfunction

END MODULE m_wavefproducts
