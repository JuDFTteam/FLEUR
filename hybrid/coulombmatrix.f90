!
!     Calculates the Coulomb matrix
!
!     v      =  < M    | v | M    >
!      k,IJ        k,I        k,J
!
!     with the mixed-basis functions M (indices I and J).
!
!     Note that
!                 *
!     v      =  v     .
!      k,JI      k,IJ
!
!     In the code: coulomb(IJ,k) = v     where only the upper triangle (I<=J) is stored.
!                                   k,IJ
!
!     The Coulomb matrix v(IJ,k) diverges at the Gamma-point. Here, we apply the decomposition
!
!              (0)        (1)   *        2-l              (0)*   (0)    (1)*        m  (1)
!     v     = v    + SUM v   * Y  (k) / k        with    v    = v   ,  v      = (-1)  v
!      k,IJ    IJ     lm  IJ    lm                        JI     IJ     JI,lm          IJ,l,-m
!
!     where a = atom index, R  = position vector, T  = Wigner-Seitz radius (scalar).
!                            a                     0
!                                    (0)
!     In the code: coulomb(IJ,1)  = v    where only the upper triangle (I<=J) is stored,
!                                    IJ
!                                    (1)
!                  coulfac(IJ,lm) = v                                    IJ,lm
!
!     For the PW contribution we have to construct plane waves within the MT spheres with the help
!     of spherical Bessel functions. The value lexp (LEXP in gwinp) is the corresponding cutoff.
!
MODULE m_coulombmatrix

CONTAINS

   SUBROUTINE coulombmatrix(mpi, fi, mpdata, hybdat, xcpot, my_k_list)
      USE m_types
      USE m_types_hybdat
      USE m_juDFT
      USE m_constants
      USE m_trafo, ONLY: symmetrize, bramat_trafo
      USE m_intgrf, ONLY: intgrf, intgrf_init
      use m_util, only: primitivef
      USE m_hsefunctional, ONLY: change_coulombmatrix
      USE m_wrapper
      USE m_io_hybinp
      use m_ylm
      use m_sphbes, only: sphbes
      use m_calc_l_m_from_lm
      use m_calc_mpsmat
      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_mpi), INTENT(IN)           :: mpi
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat
      integer, intent(in)               :: my_k_list(:)

      ! - local scalars -
      INTEGER                    :: inviop
      INTEGER                    :: nqnrm, iqnrm, iqnrm1, iqnrm2, iqnrmstart, iqnrmstep
      INTEGER                    :: itype, l, ix, iy, iy0, i, j, lm, l1, l2, m1, m2, ineq, idum, ikpt
      INTEGER                    :: lm1, lm2, itype1, itype2, ineq1, ineq2, n, n1, n2, ng
      INTEGER                    :: ic, ic1, ic2, ic3, ic4
      INTEGER                    :: igpt, igpt1, igpt2, igptp, igptp1, igptp2
      INTEGER                    :: isym, isym1, isym2, igpt0
      INTEGER                    :: ok, iatm1, iatm2
      INTEGER                    :: m, im
      INTEGER                    :: maxfac

      LOGICAL                    :: lsym

      REAL                       :: rdum, rdum1, rdum2
      REAL                       :: svol, qnorm, qnorm1, qnorm2, gnorm
      REAL                       :: fcoulfac

      COMPLEX                    :: cdum, cexp, csum

      ! - local arrays -
      INTEGER                    :: g(3)
      INTEGER, ALLOCATABLE   :: pqnrm(:, :)
      INTEGER                    :: rrot(3, 3, fi%sym%nsym), invrrot(3, 3, fi%sym%nsym)
      INTEGER, ALLOCATABLE   :: iarr(:), POINTER(:, :, :, :)!,pointer(:,:,:)
      INTEGER, ALLOCATABLE   :: nsym_gpt(:, :), sym_gpt(:, :, :)
      INTEGER                    :: nsym1(fi%kpts%nkpt + 1), sym1(fi%sym%nsym, fi%kpts%nkpt + 1)

      INTEGER, ALLOCATABLE   ::  ngptm1(:)
      INTEGER, ALLOCATABLE   ::  pgptm1(:, :)

      REAL                       :: q(3), q1(3), q2(3)
      REAL                       :: integrand(fi%atoms%jmtd), primf1(fi%atoms%jmtd), primf2(fi%atoms%jmtd)
      REAL                       :: moment(maxval(mpdata%num_radbasfn), 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), &
                                    moment2(maxval(mpdata%num_radbasfn), fi%atoms%ntype)
      REAL                       :: sphbes_var(fi%atoms%jmtd, 0:maxval(fi%hybinp%lcutm1))
      REAL                       :: sphbesmoment1(fi%atoms%jmtd, 0:maxval(fi%hybinp%lcutm1))
      REAL                       :: rarr(0:fi%hybinp%lexp + 1), rarr1(0:maxval(fi%hybinp%lcutm1))
      REAL, ALLOCATABLE   :: gmat(:, :), qnrm(:)
      REAL, ALLOCATABLE   :: sphbesmoment(:, :, :)
      REAL, ALLOCATABLE   :: sphbes0(:, :, :)
      REAL, ALLOCATABLE   :: olap(:, :, :, :), integral(:, :, :, :)
      REAL, ALLOCATABLE   :: gridf(:, :)
      REAL                       :: facA(0:MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1))
      REAL                       :: facB(0:MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1))
      REAL                       :: facC(-1:MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1))

      COMPLEX     :: structconst((2*fi%hybinp%lexp + 1)**2, fi%atoms%nat, fi%atoms%nat, fi%kpts%nkpt)             ! nw = 1
      COMPLEX     :: y((fi%hybinp%lexp + 1)**2)
      COMPLEX     :: dwgn(-maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), 0:maxval(fi%hybinp%lcutm1), fi%sym%nsym)
      COMPLEX, ALLOCATABLE   :: carr2(:, :), carr2a(:, :), carr2b(:, :)
      COMPLEX, ALLOCATABLE   :: structconst1(:, :)

      INTEGER                    :: ishift, ishift1
      INTEGER                    :: iatom, iatom1, entry_len
      INTEGER                    :: indx1, indx2, indx3, indx4
      TYPE(t_mat)                :: coul_mtmt, mat, coulmat, smat, tmp
      type(t_mat), allocatable   :: coulomb(:)

      CALL timestart("Coulomb matrix setup")
      call timestart("prep in coulomb")
      if (mpi%is_root()) write (*, *) "start of coulomb calculation"

      call mat%alloc(.True., maxval(mpdata%num_radbasfn), maxval(mpdata%num_radbasfn))

      svol = SQRT(fi%cell%vol)
      fcoulfac = fpi_const/fi%cell%vol
      maxfac = MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1)

      facA(0) = 1                    !
      facB(0) = 1                    ! Define:
      facC(-1:0) = 1                    ! facA(i)    = i!
      DO i = 1, maxfac                       ! facB(i)   = sqrt(i!)
         facA(i) = facA(i - 1)*i            ! facC(i) = (2i+1)!!
         facB(i) = facB(i - 1)*SQRT(i*1.0) !
         facC(i) = facC(i - 1)*(2*i + 1)   !
      END DO

      CALL intgrf_init(fi%atoms%ntype, fi%atoms%jmtd, fi%atoms%jri, fi%atoms%dx, fi%atoms%rmsh, gridf)

      !     Calculate the structure constant
      CALL structureconstant(structconst, fi%cell, fi%hybinp, fi%atoms, fi%kpts, mpi)

      IF (mpi%irank == 0) WRITE (oUnit, '(//A)') '### subroutine: coulombmatrix ###'

      !
      !     Matrix allocation
      !

      call timestart("coulomb allocation")
      call coul_mtmt%alloc(.False., maxval(hybdat%nbasm), maxval(hybdat%nbasm))

      allocate(coulomb(fi%kpts%nkpt))
      DO im = 1, size(my_k_list)
         ikpt = my_k_list(im)
         call coulomb(ikpt)%alloc(.False., hybdat%nbasm(ikpt), hybdat%nbasm(ikpt))
      enddo
      call timestop("coulomb allocation")

      IF (mpi%irank == 0) then
         write (oUnit,*) "Size of coulomb matrix: " //&
                            float2str(sum([(coulomb(my_k_list(i))%size_mb(), i=1,size(my_k_list))])) // " MB"
      endif

      !     Generate Symmetry:
      !     Reduce list of g-Points so that only one of each symm-equivalent is calculated
      ! calculate rotations in reciprocal space
      DO isym = 1, fi%sym%nsym
         IF (isym <= fi%sym%nop) THEN
            inviop = fi%sym%invtab(isym)
            rrot(:, :, isym) = TRANSPOSE(fi%sym%mrot(:, :, inviop))
            DO l = 0, maxval(fi%hybinp%lcutm1)
               dwgn(:, :, l, isym) = TRANSPOSE(fi%hybinp%d_wgn2(-maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), &
                                                                -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), l, isym))
            END DO
         ELSE
            inviop = isym - fi%sym%nop
            rrot(:, :, isym) = -rrot(:, :, inviop)
            dwgn(:, :, :, isym) = dwgn(:, :, :, inviop)
            DO l = 0, maxval(fi%hybinp%lcutm1)
               DO m1 = -l, l
                  DO m2 = -l, -1
                     cdum = dwgn(m1, m2, l, isym)
                     dwgn(m1, m2, l, isym) = dwgn(m1, -m2, l, isym)*(-1)**m2
                     dwgn(m1, -m2, l, isym) = cdum*(-1)**m2
                  END DO
               END DO
            END DO
         END IF
      END DO
      invrrot(:, :, :fi%sym%nop) = rrot(:, :, fi%sym%invtab)
      IF (fi%sym%nsym > fi%sym%nop) THEN
         invrrot(:, :, fi%sym%nop + 1:) = rrot(:, :, fi%sym%invtab + fi%sym%nop)
      END IF

      ! Get symmetry operations that leave bk(:,ikpt) invariant -> sym1
      nsym1 = 0
      DO ikpt = 1, fi%kpts%nkpt
         isym1 = 0
         DO isym = 1, fi%sym%nsym
            ! temporary fix until bramat_trafo is correct
            ! for systems with symmetries including translations
            IF (isym > fi%sym%nop) THEN
               isym2 = isym - fi%sym%nop
            ELSE
               isym2 = isym
            END IF
            IF (ANY(abs(fi%sym%tau(:, isym2)) > 1e-12)) CYCLE

            IF (ALL(ABS(MATMUL(rrot(:, :, isym), fi%kpts%bk(:, ikpt)) - fi%kpts%bk(:, ikpt)) < 1e-12)) THEN
               isym1 = isym1 + 1
               sym1(isym1, ikpt) = isym
            END IF
         END DO
         nsym1(ikpt) = isym1
      END DO
      ! Define reduced lists of G points -> pgptm1(:,ikpt), ikpt=1,..,nkpt
      !if(allocated(pgptm1)) deallocate(fi%hybinp%pgptm1)
      allocate (pgptm1(maxval(mpdata%n_g), fi%kpts%nkptf), source=0) !in mixedbasis
      allocate (iarr(maxval(mpdata%n_g)), source=0)
      allocate (POINTER(fi%kpts%nkpt, &
                        MINVAL(mpdata%g(1, :)) - 1:MAXVAL(mpdata%g(1, :)) + 1, &
                        MINVAL(mpdata%g(2, :)) - 1:MAXVAL(mpdata%g(2, :)) + 1, &
                        MINVAL(mpdata%g(3, :)) - 1:MAXVAL(mpdata%g(3, :)) + 1), &
                source=0)
      allocate (ngptm1, mold=mpdata%n_g)
      ngptm1 = 0

      DO ikpt = 1, fi%kpts%nkpt
         DO igpt = 1, mpdata%n_g(ikpt)
            g = mpdata%g(:, mpdata%gptm_ptr(igpt, ikpt))
            POINTER(ikpt, g(1), g(2), g(3)) = igpt
         END DO
         iarr = 0
         j = 0
         DO igpt = mpdata%n_g(ikpt), 1, -1
            IF (iarr(igpt) == 0) THEN
               j = j + 1
               pgptm1(j, ikpt) = igpt
               DO isym1 = 1, nsym1(ikpt)
                  g = MATMUL(rrot(:, :, sym1(isym1, ikpt)), mpdata%g(:, mpdata%gptm_ptr(igpt, ikpt)))
                  i = POINTER(ikpt, g(1), g(2), g(3))
                  IF (i == 0) call judft_error('coulombmatrix: zero pointer (bug?)')
                  iarr(i) = 1
               END DO
            END IF
         END DO
         ngptm1(ikpt) = j
      END DO
      deallocate (iarr)

      ! Distribute the work as equally as possible over the processes
      call timestop("prep in coulomb")

      call timestart("define gmat")
      ! Define gmat (symmetric)
      allocate (gmat((fi%hybinp%lexp + 1)**2, (fi%hybinp%lexp + 1)**2))
      gmat = 0
      lm1 = 0
      DO l1 = 0, fi%hybinp%lexp
         DO m1 = -l1, l1
            lm1 = lm1 + 1
            lm2 = 0
            lp1: DO l2 = 0, l1
               DO m2 = -l2, l2
                  lm2 = lm2 + 1
                  IF (lm2 > lm1) EXIT lp1 ! Don't cross the diagonal!
                  gmat(lm1, lm2) = facB(l1 + l2 + m2 - m1)*facB(l1 + l2 + m1 - m2)/ &
                                   (facB(l1 + m1)*facB(l1 - m1)*facB(l2 + m2)*facB(l2 - m2))/ &
                                   SQRT(1.0*(2*l1 + 1)*(2*l2 + 1)*(2*(l1 + l2) + 1))*(fpi_const)**1.5
                  gmat(lm2, lm1) = gmat(lm1, lm2)
               END DO
            END DO LP1
         END DO
      END DO
      call timestop("define gmat")

      ! Calculate moments of MT functions
      call timestart("calc moments of MT")
      DO itype = 1, fi%atoms%ntype
         DO l = 0, fi%hybinp%lcutm1(itype)
            DO i = 1, mpdata%num_radbasfn(l, itype)
               ! note that mpdata%radbasfn_mt already contains the factor rgrid
               moment(i, l, itype) = intgrf(fi%atoms%rmsh(:, itype)**(l + 1)*mpdata%radbasfn_mt(:, i, l, itype), &
                                            fi%atoms, itype, gridf)
            END DO
         END DO
         DO i = 1, mpdata%num_radbasfn(0, itype)
            moment2(i, itype) = intgrf(fi%atoms%rmsh(:, itype)**3*mpdata%radbasfn_mt(:, i, 0, itype), &
                                       fi%atoms, itype, gridf)
         END DO
      END DO
      call timestop("calc moments of MT")

      call timestart("getnorm")
      ! Look for different qnorm = |k+G|, definition of qnrm and pqnrm.
      CALL getnorm(fi%kpts, mpdata%g, mpdata%n_g, mpdata%gptm_ptr, qnrm, nqnrm, pqnrm, fi%cell)
      allocate (sphbesmoment(0:fi%hybinp%lexp, fi%atoms%ntype, nqnrm), &
                olap(maxval(mpdata%num_radbasfn), 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype, nqnrm), &
                integral(maxval(mpdata%num_radbasfn), 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype, nqnrm))
      sphbes_var = 0
      sphbesmoment = 0
      sphbesmoment1 = 0
      olap = 0
      integral = 0

      ! Calculate moments of spherical Bessel functions (for (2) and (3))              (->sphbesmoment)
      ! Calculate overlap of spherical Bessel functions with basis functions (for (2)) (->olap)
      ! Calculate overlap of sphbesmoment1(r,l)         with basis functions (for (2)) (->integral)
      ! We use           sphbes(r,l) = j_l(qr)
      ! and       sphbesmoment1(r,l) = 1/r**(l-1) * INT(0..r) r'**(l+2) * j_l(qr') dr'
      !                                + r**(l+2) * INT(r..S) r'**(1-l) * j_l(qr') dr' .

      iqnrmstart = mpi%irank + 1
      iqnrmstep = mpi%isize
      call timestop("getnorm")

      call timestart("Bessel calculation")
      !DO iqnrm = iqnrmstart, nqnrm, iqnrmstep
      do iqnrm = 1, nqnrm
         qnorm = qnrm(iqnrm)
         DO itype = 1, fi%atoms%ntype
            ng = fi%atoms%jri(itype)
            rdum = fi%atoms%rmt(itype)
            sphbes_var = 0
            sphbesmoment1 = 0
            IF (abs(qnorm) < 1e-12) THEN
               sphbesmoment(0, itype, iqnrm) = rdum**3/3
               DO i = 1, ng
                  sphbes_var(i, 0) = 1
                  sphbesmoment1(i, 0) = fi%atoms%rmsh(i, itype)**2/3 + (rdum**2 - fi%atoms%rmsh(i, itype)**2)/2
               END DO
            ELSE
               call sphbes(fi%hybinp%lexp + 1, qnorm*rdum, rarr)
               DO l = 0, fi%hybinp%lexp
                  sphbesmoment(l, itype, iqnrm) = rdum**(l + 2)*rarr(l + 1)/qnorm
               END DO
               DO i = ng, 1, -1
                  rdum = fi%atoms%rmsh(i, itype)
                  call sphbes(fi%hybinp%lcutm1(itype) + 1, qnorm*rdum, rarr)
                  DO l = 0, fi%hybinp%lcutm1(itype)
                     sphbes_var(i, l) = rarr(l)
                     IF (l /= 0) THEN; rdum1 = -rdum**(1 - l)*rarr(l - 1)
                     ELSE; rdum1 = -COS(qnorm*rdum)/qnorm
                     ENDIF
                     IF (i == ng) rarr1(l) = rdum1
                     sphbesmoment1(i, l) = (rdum**(l + 2)*rarr(l + 1)/rdum**(l + 1) &
                                            + (rarr1(l) - rdum1)*rdum**l)/qnorm
                  END DO
               END DO
            END IF
            DO l = 0, fi%hybinp%lcutm1(itype)
               DO n = 1, mpdata%num_radbasfn(l, itype)
                  ! note that mpdata%radbasfn_mt already contains one factor rgrid
                  olap(n, l, itype, iqnrm) = &
                     intgrf(fi%atoms%rmsh(:, itype)*mpdata%radbasfn_mt(:, n, l, itype)*sphbes_var(:, l), &
                            fi%atoms, itype, gridf)

                  integral(n, l, itype, iqnrm) = &
                     intgrf(fi%atoms%rmsh(:, itype)*mpdata%radbasfn_mt(:, n, l, itype)*sphbesmoment1(:, l), &
                            fi%atoms, itype, gridf)

               END DO
            END DO
         END DO
      END DO
      call timestop("Bessel calculation")

      !
      !     (1) Case < MT | v | MT >
      !


      !       (1a) r,r' in same MT
      call timestart("loop 1")
      ix = 0
      iy = 0
      iy0 = 0
      DO itype = 1, fi%atoms%ntype
         DO ineq = 1, fi%atoms%neq(itype)
            ! Here the diagonal block matrices do not depend on ineq. In (1b) they do depend on ineq, though,
            DO l = 0, fi%hybinp%lcutm1(itype)
               mat%matsize1=mpdata%num_radbasfn(l, itype)
               mat%matsize2=mpdata%num_radbasfn(l, itype)
               DO n2 = 1, mpdata%num_radbasfn(l, itype)
                  ! note that mpdata%radbasfn_mt already contains the factor rgrid
                  CALL primitivef(primf1, mpdata%radbasfn_mt(:, n2, l, itype) &
                                    *fi%atoms%rmsh(:, itype)**(l + 1), fi%atoms%rmsh, fi%atoms%dx, &
                                    fi%atoms%jri, fi%atoms%jmtd, itype, fi%atoms%ntype)
                  ! -itype is to enforce inward integration
                  CALL primitivef(primf2, mpdata%radbasfn_mt(:fi%atoms%jri(itype), n2, l, itype) &
                                    /fi%atoms%rmsh(:fi%atoms%jri(itype), itype)**l, fi%atoms%rmsh, fi%atoms%dx, &
                                    fi%atoms%jri, fi%atoms%jmtd, -itype, fi%atoms%ntype)

                  primf1(:fi%atoms%jri(itype)) = primf1(:fi%atoms%jri(itype))/fi%atoms%rmsh(:fi%atoms%jri(itype), itype)**l
                  primf2 = primf2*fi%atoms%rmsh(:, itype)**(l + 1)

                  DO n1 = 1, n2
                     integrand = mpdata%radbasfn_mt(:, n1, l, itype)*(primf1 + primf2)
                     mat%data_r(n1, n2) = fpi_const/(2*l + 1) * intgrf(integrand, fi%atoms, itype, gridf)
                  END DO
               END DO

               ! distribute mat for m=-l,l on coulomb in block-matrix form
               DO M = -l, l
                  DO n2 = 1, mpdata%num_radbasfn(l, itype)
                     ix = ix + 1
                     iy = iy0
                     DO n1 = 1, n2
                        iy = iy + 1
                        i = ix*(ix - 1)/2 + iy
                        j = n2*(n2 - 1)/2 + n1
                        coul_mtmt%data_c(iy, ix) = mat%data_r(n1, n2)
                     END DO
                  END DO
                  iy0 = ix
               END DO

            END DO
         END DO
      END DO
      call timestop("loop 1")
      
      call coul_mtmt%u2l()

      call coulmat%alloc(.False., hybdat%nbasp, hybdat%nbasp)

      DO im = 1, size(my_k_list)
         ikpt = my_k_list(im)

         ! only the first rank handles the MT-MT part
         call timestart("MT-MT part")

         ix = 0
         ic2 = 0
         DO itype2 = 1, fi%atoms%ntype
            DO ineq2 = 1, fi%atoms%neq(itype2)
               ic2 = ic2 + 1
               lm2 = 0
               DO l2 = 0, fi%hybinp%lcutm1(itype2)
                  DO m2 = -l2, l2
                     lm2 = lm2 + 1
                     DO n2 = 1, mpdata%num_radbasfn(l2, itype2)
                        ix = ix + 1

                        iy = 0
                        ic1 = 0
                        lp2: DO itype1 = 1, itype2
                           DO ineq1 = 1, fi%atoms%neq(itype1)
                              ic1 = ic1 + 1
                              lm1 = 0
                              DO l1 = 0, fi%hybinp%lcutm1(itype1)
                                 DO m1 = -l1, l1
                                    lm1 = lm1 + 1
                                    DO n1 = 1, mpdata%num_radbasfn(l1, itype1)
                                       iy = iy + 1
                                       IF (iy > ix) EXIT lp2 ! Don't cross the diagonal!
                                       rdum = (-1)**(l2 + m2)*moment(n1, l1, itype1)*moment(n2, l2, itype2)*gmat(lm1, lm2)
                                       l = l1 + l2
                                       lm = l**2 + l + m1 - m2 + 1
                                       idum = ix*(ix - 1)/2 + iy
                                       coulmat%data_c(iy, ix) = coul_mtmt%data_c(iy,ix) &
                                                            + EXP(CMPLX(0.0, 1.0)*tpi_const* &
                                                                  dot_PRODUCT(fi%kpts%bk(:, ikpt), &
                                                                              fi%atoms%taual(:, ic2) - fi%atoms%taual(:, ic1))) &
                                                            *rdum*structconst(lm, ic1, ic2, ikpt)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO lp2

                     END DO
                  END DO
               END DO
            END DO
         END DO
         
         call coulmat%u2l()
         IF (fi%sym%invs) THEN
            !symmetrize makes the Coulomb matrix real symmetric               
            CALL symmetrize(coulmat%data_c, hybdat%nbasp, hybdat%nbasp, 3, .FALSE., &
               fi%atoms, fi%hybinp%lcutm1, maxval(fi%hybinp%lcutm1), &
               mpdata%num_radbasfn, fi%sym)
         ENDIF

         call coulomb(ikpt)%copy(coulmat, 1,1)
         call timestop("MT-MT part")

      END DO

      IF (maxval(mpdata%n_g) /= 0) THEN ! skip calculation of plane-wave contribution if mixed basis does not contain plane waves

         !
         !     (2) Case < MT | v | PW >
         !

         !     (2a) r in MT, r' everywhere
         !     (2b) r,r' in same MT
         !     (2c) r,r' in different MT

         call timestart("loop over interst.")
         DO im = 1, size(my_k_list)
            ikpt = my_k_list(im)
            call loop_over_interst(fi, hybdat, mpdata, structconst, sphbesmoment, moment, moment2, &
                                   qnrm, facc, gmat, integral, olap, pqnrm, pgptm1, ngptm1, ikpt, coulomb(ikpt))

            call coulomb(ikpt)%u2l()
         END DO
         call timestop("loop over interst.")
         deallocate (olap, integral)

         !
         !     (3) Case < PW | v | PW >
         !
         !     (3a) r,r' everywhere; r everywhere, r' in MT; r in MT, r' everywhere
         ! Calculate the hermitian matrix smat(i,j) = sum(a) integral(MT(a)) exp[i(Gj-Gi)r] dr

         call calc_mpsmat(fi, mpdata, smat)

         ! Coulomb matrix, contribution (3a)
         call timestart("coulomb matrix 3a")
         DO im = 1, size(my_k_list)
            ikpt = my_k_list(im)

            DO igpt0 = 1, ngptm1(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
               ix = hybdat%nbasp + igpt2
               iy = hybdat%nbasp
               q2 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp2), fi%cell%bmat)
               rdum2 = SUM(q2**2)
               IF (abs(rdum2) > 1e-12) rdum2 = fpi_const/rdum2

               DO igpt1 = 1, igpt2
                  igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
                  iy = iy + 1
                  q1 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp1), fi%cell%bmat)
                  idum = ix*(ix - 1)/2 + iy
                  rdum1 = SUM(q1**2)
                  IF (abs(rdum1) > 1e-12) rdum1 = fpi_const/rdum1

                  IF (ikpt == 1) THEN
                     IF (igpt1 /= 1) THEN
                        coulomb(1)%data_c(iy,ix) = -smat%data_c(igptp1, igptp2)*rdum1/fi%cell%vol
                     END IF
                     IF (igpt2 /= 1) THEN
                        coulomb(1)%data_c(iy,ix) = coulomb(1)%data_c(iy,ix) - smat%data_c(igptp1, igptp2)*rdum2/fi%cell%vol
                     END IF
                  ELSE
                     coulomb(ikpt)%data_c(iy,ix) = -smat%data_c(igptp1, igptp2)*(rdum1 + rdum2)/fi%cell%vol
                  END IF
               END DO
               IF (ikpt /= 1 .OR. igpt2 /= 1) THEN                  !
                  coulomb(ikpt)%data_c(iy,ix) = coulomb(ikpt)%data_c(iy,ix) + rdum2
               END IF                                            !
            END DO
            call coulomb(ikpt)%u2l()
         END DO
         call timestop("coulomb matrix 3a")
         !     (3b) r,r' in different MT

         call timestart("coulomb matrix 3b")
         DO im = 1, size(my_k_list)
            ikpt = my_k_list(im)
            if (mpi%is_root()) write (*, *) "coulomb pw-loop nk: ("//int2str(ikpt)//"/"//int2str(fi%kpts%nkpt)//")"
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
                      structconst1(fi%atoms%nat, (2*fi%hybinp%lexp + 1)**2))
            carr2 = 0; structconst1 = 0

            DO igpt0 = 1, ngptm1(ikpt)!1,ngptm1(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               ix = hybdat%nbasp + igpt2
               igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
               iqnrm2 = pqnrm(igpt2, ikpt)
               iatom = 0
               carr2 = 0
               call timestart("itype loops")
               DO itype2 = 1, fi%atoms%ntype
                  DO ineq2 = 1, fi%atoms%neq(itype2)
                     iatom = iatom + 1
                     cexp = CONJG(carr2b(iatom, igpt2))
                     structconst1(:, :) = transpose(structconst(:, :, iatom, ikpt))
                     ! this is a nested loop over
                     ! l=1..hyb%lexp{
                     !    m=-l..l{}
                     ! }
                     !$OMP PARALLEL DO default(none) private(lm1,l1,m1,lm2,l2,m2,cdum,l,lm) &
                     !$OMP shared(fi, sphbesmoment, itype2, iqnrm2, cexp, carr2a, igpt2, carr2, gmat, structconst1) &
                     !$OMP collapse(2)
                     DO lm1 = 1, (fi%hybinp%lexp+1)**2
                        do lm2 = 1, (fi%hybinp%lexp+1)**2
                           call calc_l_m_from_lm(lm1, l1, m1)
                           call calc_l_m_from_lm(lm2, l2, m2)
                           cdum = (-1)**(l2 + m2)*sphbesmoment(l2, itype2, iqnrm2)*cexp*carr2a(lm2, igpt2)
                           l = l1 + l2
                           lm = l**2 + l - l1 - m2 + (m1 + l1) + 1
                           carr2(:, lm1) = carr2(:, lm1) + cdum*gmat(lm1, lm2)*structconst1(:, lm)
                        END DO
                     enddo
                     !$OMP end parallel do
                  END DO
               END DO
               call timestop("itype loops")

               call timestart("igpt1")
               iy = hybdat%nbasp
               DO igpt1 = 1, igpt2
                  iy = iy + 1
                  igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
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
                  coulomb(ikpt)%data_c(iy,ix) = coulomb(ikpt)%data_c(iy,ix) + csum/fi%cell%vol
               END DO
               call timestop("igpt1")
            END DO
            deallocate (carr2, carr2a, carr2b, structconst1)
            call coulomb(ikpt)%u2l() 
            call timestop("loop over plane waves")
         END DO !ikpt
         call timestop("coulomb matrix 3b")

         if(any(my_k_list == 1)) then
            !     Add corrections from higher orders in (3b) to coulomb(:,1)
            ! (1) igpt1 > 1 , igpt2 > 1  (finite G vectors)
            call timestart("add corrections from higher orders")
            rdum = (fpi_const)**(1.5)/fi%cell%vol**2*gmat(1, 1)
            !$OMP PARALLEL DO default(none) schedule(dynamic) &
            !$OMP private(igpt0, igpt2, ix, iqnrm2, igptp2, q2, qnorm2, igpt1, iy)&
            !$OMP private(iqnrm1,igptp1,q1,qnorm1,rdum1,iatm1,itype1,iatm2,itype2,cdum)&
            !$OMP shared(ngptm1, pgptm1, hybdat, mpdata, coulomb, sphbesmoment,pqnrm,fi,rdum)
            DO igpt0 = 1, ngptm1(1)
               igpt2 = pgptm1(igpt0, 1)
               IF (igpt2 == 1) CYCLE
               ix = hybdat%nbasp + igpt2
               iqnrm2 = pqnrm(igpt2, 1)
               igptp2 = mpdata%gptm_ptr(igpt2, 1)
               q2 = MATMUL(mpdata%g(:, igptp2), fi%cell%bmat)
               qnorm2 = norm2(q2)
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
                        coulomb(1)%data_c(iy, ix) = coulomb(1)%data_c(iy, ix) + rdum*cdum*( &
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
            END DO
            !$OMP END PARALLEL DO

            call coulomb(1)%u2l()   

            ! (2) igpt1 = 1 , igpt2 > 1  (first G vector vanishes, second finite)
            iy = hybdat%nbasp + 1
            DO igpt0 = 1, ngptm1(1)
               igpt2 = pgptm1(igpt0, 1); IF (igpt2 == 1) CYCLE
               ix = hybdat%nbasp + igpt2
               iqnrm2 = pqnrm(igpt2, 1)
               igptp2 = mpdata%gptm_ptr(igpt2, 1)
               qnorm2 = qnrm(iqnrm2)
               idum = ix*(ix - 1)/2 + iy
               DO itype1 = 1, fi%atoms%ntype
                  DO ineq1 = 1, fi%atoms%neq(itype1)
                     ic2 = 0
                     DO itype2 = 1, fi%atoms%ntype
                        DO ineq2 = 1, fi%atoms%neq(itype2)
                           ic2 = ic2 + 1
                           cdum = EXP(CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT(mpdata%g(:, igptp2), fi%atoms%taual(:, ic2)))
                           coulomb(1)%data_c(iy, ix) = coulomb(1)%data_c(iy, ix) &
                                             + rdum*cdum*fi%atoms%rmt(itype1)**3*( &
                                             +sphbesmoment(0, itype2, iqnrm2)/30*fi%atoms%rmt(itype1)**2 &
                                             - sphbesmoment(2, itype2, iqnrm2)/18 &
                                             + sphbesmoment(1, itype2, iqnrm2)/6/qnorm2)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            call coulomb(1)%u2l()


            ! (2) igpt1 = 1 , igpt2 = 1  (vanishing G vectors)
            iy = hybdat%nbasp + 1
            ix = hybdat%nbasp + 1
            idum = ix*(ix - 1)/2 + iy
            DO itype1 = 1, fi%atoms%ntype
               DO ineq1 = 1, fi%atoms%neq(itype1)
                  DO itype2 = 1, fi%atoms%ntype
                     DO ineq2 = 1, fi%atoms%neq(itype2)
                        coulomb(1)%data_c(iy, ix) = coulomb(1)%data_c(iy, ix) &
                                          + rdum*fi%atoms%rmt(itype1)**3*fi%atoms%rmt(itype2)**3* &
                                          (fi%atoms%rmt(itype1)**2 + fi%atoms%rmt(itype2)**2)/90
                     END DO
                  END DO
               END DO
            END DO
            call coulomb(1)%u2l()
            call timestop("add corrections from higher orders")
         endif ! 1 in my_k_list


         !     (3c) r,r' in same MT

         ! Calculate sphbesintegral
         call timestart("sphbesintegral")
         allocate (sphbes0(-1:fi%hybinp%lexp + 2, fi%atoms%ntype, nqnrm),&
              &           carr2((fi%hybinp%lexp + 1)**2, maxval(mpdata%n_g)))
         sphbes0 = 0; carr2 = 0
         DO iqnrm = 1, nqnrm
            DO itype = 1, fi%atoms%ntype
               rdum = qnrm(iqnrm)*fi%atoms%rmt(itype)
               call sphbes(fi%hybinp%lexp + 2, rdum, sphbes0(0, itype, iqnrm))
               IF (abs(rdum) > 1e-12) sphbes0(-1, itype, iqnrm) = COS(rdum)/rdum
            END DO
         END DO
         call timestop("sphbesintegral")

         call timestart("loop 2")
         DO im = 1, size(my_k_list)
            ikpt = my_k_list(im)
            call timestart("harmonics setup")
            DO igpt = 1, mpdata%n_g(ikpt)
               igptp = mpdata%gptm_ptr(igpt, ikpt)
               q = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)
               call ylm4(fi%hybinp%lexp, q, carr2(:, igpt))
            END DO
            call timestop("harmonics setup")
            call perform_double_g_loop(fi, hybdat, mpdata, sphbes0, carr2, ngptm1,pgptm1,&
                                       pqnrm,qnrm, nqnrm, ikpt, coulomb(ikpt))
            call coulomb(ikpt)%u2l()
         END DO
         call timestop("loop 2")
         deallocate (carr2)

         !
         !     Symmetry-equivalent G vectors
         !
         ! All elements are needed so send all data to all processes treating the
         ! respective k-points

         allocate (carr2(maxval(hybdat%nbasm), 2), iarr(maxval(mpdata%n_g)))
         allocate (nsym_gpt(mpdata%num_gpts(), fi%kpts%nkpt), &
                   sym_gpt(MAXVAL(nsym1), mpdata%num_gpts(), fi%kpts%nkpt))
         nsym_gpt = 0; sym_gpt = 0
         call timestart("loop 3")
         DO im = 1, size(my_k_list)
            ikpt = my_k_list(im)
            carr2 = 0; iarr = 0
            iarr(pgptm1(:ngptm1(ikpt), ikpt)) = 1
            DO igpt0 = 1, ngptm1(ikpt)
               lsym = (1 <= igpt0) .AND. (ngptm1(ikpt) >= igpt0)
               igpt2 = pgptm1(igpt0, ikpt)
               carr2(:hybdat%nbasm(ikpt),2) = coulomb(ikpt)%data_c(:hybdat%nbasm(ikpt),hybdat%nbasp + igpt2)

               IF (lsym) THEN
                  ic = 1
                  sym_gpt(ic, igpt0, ikpt) = igpt2
               END IF
               DO isym1 = 2, nsym1(ikpt)
                  isym = sym1(isym1, ikpt)
                  CALL bramat_trafo(carr2(:, 2), igpt2, ikpt, isym, .FALSE., POINTER(ikpt, :, :, :), &
                                    fi%sym, rrot(:, :, isym), invrrot(:, :, isym), mpdata, fi%hybinp, &
                                    fi%kpts, maxval(fi%hybinp%lcutm1), fi%atoms, fi%hybinp%lcutm1, &
                                    mpdata%num_radbasfn, maxval(mpdata%num_radbasfn), dwgn(:, :, :, isym), &
                                    hybdat%nbasp, hybdat%nbasm, carr2(:, 1), igpt1)
                  IF (iarr(igpt1) == 0) THEN
                     CALL bramat_trafo(carr2(:, 2), igpt2, ikpt, isym, .TRUE., POINTER(ikpt, :, :, :), &
                                       fi%sym, rrot(:, :, isym), invrrot(:, :, isym), mpdata, fi%hybinp, &
                                       fi%kpts, maxval(fi%hybinp%lcutm1), fi%atoms, fi%hybinp%lcutm1, &
                                       mpdata%num_radbasfn, maxval(mpdata%num_radbasfn), &
                                       dwgn(:, :, :, isym), hybdat%nbasp, hybdat%nbasm, carr2(:, 1), igpt1)
                     l = (hybdat%nbasp + igpt1 - 1)*(hybdat%nbasp + igpt1)/2
                     coulomb(ikpt)%data_c(:hybdat%nbasp + igpt1,hybdat%nbasp + igpt1) = carr2(:hybdat%nbasp + igpt1, 1)
                     coulomb(ikpt)%data_c(hybdat%nbasp + igpt1,:hybdat%nbasp + igpt1) = conjg(carr2(:hybdat%nbasp + igpt1, 1))

                     iarr(igpt1) = 1
                     IF (lsym) THEN
                        ic = ic + 1
                        sym_gpt(ic, igpt0, ikpt) = igpt1
                     END IF
                  END IF
               END DO
               nsym_gpt(igpt0, ikpt) = ic
            END DO ! igpt0
            call coulomb(ikpt)%u2l()
         END DO ! ikpt
         call timestop("loop 3")
         call timestart("gap 1:")
         deallocate (carr2, iarr, pgptm1)
      END IF
      deallocate (qnrm, pqnrm)

      IF (xcpot%is_name("hse") .OR. xcpot%is_name("vhse")) THEN
         !
         ! The HSE functional is realized subtracting erf/r from
         ! the normal Coulomb matrix
         !
         call judft_error("HSE is not implemented")
      ELSE
         if(any(my_k_list == 1)) then
            CALL subtract_sphaverage(fi%sym, fi%cell, fi%atoms, mpdata, &
                                    fi%hybinp, hybdat, hybdat%nbasm, gridf, coulomb(1))
         endif
      END IF
      
      ! transform Coulomb matrix to the biorthogonal set
      ! REFACTORING HINT: THIS IS DONE WTIH THE INVERSE OF OLAP
      ! IT CAN EASILY BE REWRITTEN AS A LINEAR SYSTEM
      call timestop("gap 1:")
      DO im = 1, size(my_k_list)
         ikpt = my_k_list(im)
         call apply_inverse_olaps(mpdata, fi%atoms, fi%cell, hybdat, fi%sym, fi%kpts, ikpt, coulomb(ikpt))
         call coulomb(ikpt)%u2l()
      enddo

      !call plot_coulombmatrix() -> code was shifted to plot_coulombmatrix.F90
      !
      ! rearrange coulomb matrix
      !
      if(.not. allocated(hybdat%coul)) allocate(hybdat%coul(fi%kpts%nkpt))
      call timestart("loop bla")
      DO ikpt = 1, fi%kpts%nkpt
         call hybdat%coul(ikpt)%init()
      enddo

      DO im = 1, size(my_k_list)
         ikpt = my_k_list(im)
         ! unpack coulomb into coulomb(ikpt)

         ! only one processor per k-point calculates MT convolution
         !
         ! store m-independent part of Coulomb matrix in MT spheres
         ! in coulomb_mt1(:mpdata%num_radbasfn(l,itype)-1,:mpdata%num_radbasfn(l,itype)-1,l,itype)
         !
         call timestart("m-indep. part of coulomb mtx")
         indx1 = 0
         DO itype = 1, fi%atoms%ntype
            DO ineq = 1, fi%atoms%neq(itype)
               DO l = 0, fi%hybinp%lcutm1(itype)

                  IF (ineq == 1) THEN
                     DO n = 1, mpdata%num_radbasfn(l, itype) - 1
                        if (fi%sym%invs) THEN
                           hybdat%coul(ikpt)%mt1_r(n, 1:mpdata%num_radbasfn(l, itype) - 1, l, itype) &
                              = real(coulomb(ikpt)%data_c(indx1 + n, indx1 + 1:indx1 + mpdata%num_radbasfn(l, itype) - 1))
                        else
                           hybdat%coul(ikpt)%mt1_c(n, 1:mpdata%num_radbasfn(l, itype) - 1, l, itype) &
                              = real(coulomb(ikpt)%data_c(indx1 + n, indx1 + 1:indx1 + mpdata%num_radbasfn(l, itype) - 1))
                        endif
                     END DO
                  END IF

                  indx1 = indx1 + (2*l + 1)*mpdata%num_radbasfn(l, itype)
               END DO
            END DO
         END DO
         call timestop("m-indep. part of coulomb mtx")

         !
         ! store m-dependent and atom-dependent part of Coulomb matrix in MT spheres
         ! in coulomb_mt2(:mpdata%num_radbasfn(l,itype)-1,-l:l,l,iatom)
         !
         call timestart("m-dep. part of coulomb mtx")
         indx1 = 0
         iatom = 0
         DO itype = 1, fi%atoms%ntype
            DO ineq = 1, fi%atoms%neq(itype)
               iatom = iatom + 1
               DO l = 0, fi%hybinp%lcutm1(itype)
                  DO M = -l, l
                     if (fi%sym%invs) THEN
                        hybdat%coul(ikpt)%mt2_r(:mpdata%num_radbasfn(l, itype) - 1, M, l, iatom) &
                           = real(coulomb(ikpt)%data_c(indx1 + 1:indx1 + mpdata%num_radbasfn(l, itype) - 1, indx1 + mpdata%num_radbasfn(l, itype)))
                     else
                        hybdat%coul(ikpt)%mt2_c(:mpdata%num_radbasfn(l, itype) - 1, M, l, iatom) &
                           = coulomb(ikpt)%data_c(indx1 + 1:indx1 + mpdata%num_radbasfn(l, itype) - 1, indx1 + mpdata%num_radbasfn(l, itype))
                     endif

                     indx1 = indx1 + mpdata%num_radbasfn(l, itype)

                  END DO
               END DO
            END DO
         END DO
         call timestop("m-dep. part of coulomb mtx")

         !
         ! due to the subtraction of the divergent part at the Gamma point
         ! additional contributions occur
         !
         call timestart("gamma point treatment")
         IF (ikpt == 1) THEN
            !
            ! store the contribution of the G=0 plane wave with the MT l=0 functions in
            ! coulomb_mt2(:mpdata%num_radbasfn(l=0,itype),0,maxval(fi%hybinp%lcutm1)+1,iatom)
            !
            ic = 0
            iatom = 0
            DO itype = 1, fi%atoms%ntype
               DO ineq = 1, fi%atoms%neq(itype)
                  iatom = iatom + 1
                  DO n = 1, mpdata%num_radbasfn(0, itype) - 1
                     if (fi%sym%invs) THEN
                        hybdat%coul(ikpt)%mt2_r(n, 0, maxval(fi%hybinp%lcutm1) + 1, iatom) =  real(coulomb(ikpt)%data_c(ic + n, hybdat%nbasp + 1))
                     else
                        hybdat%coul(ikpt)%mt2_c(n, 0, maxval(fi%hybinp%lcutm1) + 1, iatom) = coulomb(ikpt)%data_c(ic + n, hybdat%nbasp + 1)
                     endif
                  END DO
                  ic = ic + SUM([((2*l + 1)*mpdata%num_radbasfn(l, itype), l=0, fi%hybinp%lcutm1(itype))])
               END DO
            END DO

            !
            ! store the contributions between the MT s-like functions at atom1 and
            ! and the constant function at a different atom2
            !
            iatom = 0
            ic = 0
            DO itype = 1, fi%atoms%ntype
               ishift = SUM([((2*l + 1)*mpdata%num_radbasfn(l, itype), l=0, fi%hybinp%lcutm1(itype))])
               DO ineq = 1, fi%atoms%neq(itype)
                  iatom = iatom + 1
                  ic1 = ic + mpdata%num_radbasfn(0, itype)

                  iatom1 = 0
                  ic2 = 0
                  DO itype1 = 1, fi%atoms%ntype
                     ishift1 = SUM([((2*l1 + 1)*mpdata%num_radbasfn(l1, itype1), l1=0, fi%hybinp%lcutm1(itype1))])
                     DO ineq1 = 1, fi%atoms%neq(itype1)
                        iatom1 = iatom1 + 1
                        ic3 = ic2 + 1
                        ic4 = ic3 + mpdata%num_radbasfn(0, itype1) - 2

                        IF (fi%sym%invs) THEN
                           hybdat%coul(ikpt)%mt3_r(:mpdata%num_radbasfn(0, itype1) - 1, iatom, iatom1) = real(coulomb(ikpt)%data_c(ic1, ic3:ic4))
                        ELSE
                           hybdat%coul(ikpt)%mt3_c(:mpdata%num_radbasfn(0, itype1) - 1, iatom, iatom1) &
                              = CONJG(coulomb(ikpt)%data_c(ic1, ic3:ic4))
                        ENDIF
                        ic2 = ic2 + ishift1
                     END DO
                  END DO

                  ic = ic + ishift
               END DO
            END DO

            !test
            iatom = 0
            DO itype = 1, fi%atoms%ntype
               DO ineq = 1, fi%atoms%neq(itype)
                  iatom = iatom + 1
                  if (fi%sym%invs) THEN
                     IF (MAXVAL(ABS(hybdat%coul(ikpt)%mt2_r(:mpdata%num_radbasfn(0, itype) - 1, 0, 0, iatom) &
                                    - hybdat%coul(ikpt)%mt3_r(:mpdata%num_radbasfn(0, itype) - 1, iatom, iatom))) > 1E-08) &
                        call judft_error('coulombmatrix: coulomb_mt2 and coulomb_mt3 are inconsistent')

                  else
                     IF (MAXVAL(ABS(hybdat%coul(ikpt)%mt2_c(:mpdata%num_radbasfn(0, itype) - 1, 0, 0,iatom) &
                                    - hybdat%coul(ikpt)%mt3_c(:mpdata%num_radbasfn(0, itype) - 1, iatom,iatom))) > 1E-08) &
                        call judft_error('coulombmatrix: coulomb_mt2 and coulomb_mt3 are inconsistent')
                  endif
               END DO
            END DO
         END IF
         call timestop("gamma point treatment")

         !
         ! add the residual MT contributions, i.e. those functions with an moment,
         ! to the matrix coulomb_mtir, which is fully occupied
         !
         call timestart("residual MT contributions")
         ic = 0
         DO itype = 1, fi%atoms%ntype
            DO ineq = 1, fi%atoms%neq(itype)
               DO l = 0, fi%hybinp%lcutm1(itype)
                  DO M = -l, l
                     ic = ic + 1
                  END DO
               END DO
            END DO
         END DO

         indx1 = 0; indx2 = 0; indx3 = 0; indx4 = 0
         DO itype = 1, fi%atoms%ntype
            DO ineq = 1, fi%atoms%neq(itype)
               DO l = 0, fi%hybinp%lcutm1(itype)
                  DO M = -l, l
                     indx1 = indx1 + 1
                     indx3 = indx3 + mpdata%num_radbasfn(l, itype)

                     indx2 = 0
                     indx4 = 0
                     DO itype1 = 1, fi%atoms%ntype
                        DO ineq1 = 1, fi%atoms%neq(itype1)
                           DO l1 = 0, fi%hybinp%lcutm1(itype1)
                              DO m1 = -l1, l1
                                 indx2 = indx2 + 1
                                 indx4 = indx4 + mpdata%num_radbasfn(l1, itype1)
                                 IF (indx4 < indx3) CYCLE
                                 IF (fi%sym%invs) THEN
                                    hybdat%coul(ikpt)%mtir_r(indx1, indx2) = real(coulomb(ikpt)%data_c(indx3, indx4))
                                    hybdat%coul(ikpt)%mtir_r(indx2, indx1) = hybdat%coul(ikpt)%mtir_r(indx1, indx2)
                                 ELSE
                                    hybdat%coul(ikpt)%mtir_c(indx1, indx2) = coulomb(ikpt)%data_c(indx3, indx4)
                                    hybdat%coul(ikpt)%mtir_c(indx2, indx1) = CONJG(hybdat%coul(ikpt)%mtir_c(indx1, indx2))
                                 ENDIF
                              END DO
                           END DO
                        END DO
                     END DO

                     DO igpt = 1, mpdata%n_g(ikpt)
                        indx2 = indx2 + 1
                        IF (fi%sym%invs) THEN
                           hybdat%coul(ikpt)%mtir_r(indx1, indx2) = real(coulomb(ikpt)%data_c(indx3, hybdat%nbasp + igpt))
                           hybdat%coul(ikpt)%mtir_r(indx2, indx1) = hybdat%coul(ikpt)%mtir_r(indx1, indx2)
                        ELSE
                           hybdat%coul(ikpt)%mtir_c(indx1, indx2) = coulomb(ikpt)%data_c(indx3, hybdat%nbasp + igpt)
                           hybdat%coul(ikpt)%mtir_c(indx2, indx1) = CONJG(hybdat%coul(ikpt)%mtir_c(indx1, indx2))
                        ENDIF

                     END DO

                  END DO
               END DO
            END DO
         END DO
         call timestop("residual MT contributions")

         IF (indx1 /= ic) call judft_error('coulombmatrix: error index counting')

         !
         ! add ir part to the matrix coulomb_mtir
         !
         if (fi%sym%invs) THEN
            hybdat%coul(ikpt)%mtir_r(ic + 1:ic + mpdata%n_g(ikpt), ic + 1:ic + mpdata%n_g(ikpt)) &
               = real(coulomb(ikpt)%data_c(hybdat%nbasp + 1:hybdat%nbasm(ikpt), hybdat%nbasp + 1:hybdat%nbasm(ikpt)))
            ic2 = indx1 + mpdata%n_g(ikpt)
         else
            hybdat%coul(ikpt)%mtir_c(ic + 1:ic + mpdata%n_g(ikpt), ic + 1:ic + mpdata%n_g(ikpt)) &
               = coulomb(ikpt)%data_c(hybdat%nbasp + 1:hybdat%nbasm(ikpt), hybdat%nbasp + 1:hybdat%nbasm(ikpt))
            ic2 = indx1 + mpdata%n_g(ikpt)
         end if

         call coulomb(ikpt)%free()
      END DO ! ikpt
      call timestop("loop bla")
      CALL timestop("Coulomb matrix setup")

   END SUBROUTINE coulombmatrix

   !     Calculate body of Coulomb matrix at Gamma point: v_IJ = SUM(G) c^*_IG c_JG 4*pi/G**2 .
   !     For this we must subtract from coulomb(:,1) the spherical average of a term that comes
   !     from the fact that MT functions have k-dependent Fourier coefficients (see script).
   SUBROUTINE subtract_sphaverage(sym, cell, atoms, mpdata, hybinp, hybdat, nbasm1, gridf, coulomb)

      USE m_types
      USE m_constants
      USE m_wrapper
      USE m_trafo
      USE m_util
      use m_intgrf
      USE m_olap
      IMPLICIT NONE

      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_mpdata), intent(in)  :: mpdata
      TYPE(t_hybinp), INTENT(IN)    :: hybinp
      TYPE(t_hybdat), INTENT(IN)    :: hybdat

      INTEGER, INTENT(IN)    :: nbasm1(:)
      REAL, INTENT(IN)    :: gridf(:, :)
      type(t_mat), intent(inout) :: coulomb

      ! - local scalars -
      INTEGER               :: l, i, j, n, nn, itype, ieq, M

      ! - local arrays -
      TYPE(t_mat) :: olap
      !COMPLEX , ALLOCATABLE :: constfunc(:)  !can also be real in inversion case
      COMPLEX      :: coeff(nbasm1(1)), cderiv(nbasm1(1), -1:1), claplace(nbasm1(1))

      call timestart("subtract_sphaverage")
      CALL olap%alloc(sym%invs, mpdata%n_g(1), mpdata%n_g(1), 0.)

      n = nbasm1(1)
      nn = n*(n + 1)/2
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(1), 1)), mpdata%n_g(1), atoms, cell)

      ! Define coefficients (coeff) and their derivatives (cderiv,claplace)
      coeff = 0
      cderiv = 0
      claplace = 0
      j = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, hybinp%lcutm1(itype)
               DO M = -l, l
                  DO i = 1, mpdata%num_radbasfn(l, itype)
                     j = j + 1
                     IF (l == 0) THEN
                        coeff(j) = SQRT(fpi_const) &
                                   *intgrf(atoms%rmsh(:, itype)*mpdata%radbasfn_mt(:, i, 0, itype), &
                                           atoms, itype, gridf) &
                                   /SQRT(cell%vol)

                        claplace(j) = -SQRT(fpi_const) &
                                      *intgrf(atoms%rmsh(:, itype)**3*mpdata%radbasfn_mt(:, i, 0, itype), &
                                              atoms, itype, gridf) &
                                      /SQRT(cell%vol)

                     ELSE IF (l == 1) THEN
                        cderiv(j, M) = -SQRT(fpi_const/3)*CMPLX(0.0, 1.0) &
                                       *intgrf(atoms%rmsh(:, itype)**2*mpdata%radbasfn_mt(:, i, 1, itype), &
                                               atoms, itype, gridf) &
                                       /SQRT(cell%vol)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      IF (olap%l_real) THEN
         coeff(hybdat%nbasp + 1:n) = olap%data_r(1, 1:n - hybdat%nbasp)
      else
         coeff(hybdat%nbasp + 1:n) = olap%data_c(1, 1:n - hybdat%nbasp)
      END IF
      IF (sym%invs) THEN
         CALL symmetrize(coeff, 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(claplace, 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(cderiv(:, -1), 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(cderiv(:, 0), 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(cderiv(:, 1), 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
      ENDIF
      ! Subtract head contributions from coulomb(:nn,1) to obtain the body
      l = 0
      DO j = 1, n
         DO i = 1, j
            l = l + 1
            coulomb%data_c(i,j) = coulomb%data_c(i,j) - fpi_const/3 &
                                       *(dot_PRODUCT(cderiv(i, :), cderiv(j, :)) &
                                       + (CONJG(coeff(i))*claplace(j) &
                                          + CONJG(claplace(i))*coeff(j))/2)
         END DO
      END DO

      call coulomb%u2l()
      call timestop("subtract_sphaverage")
   END SUBROUTINE subtract_sphaverage

   !     -----------------------------------------------------------------------------------------------

   !     Calculates the structure constant
   !                                                        1               *      ^
   !     structconst(lm,ic1,ic2,k) = SUM exp(ikT) -----------------------  Y  ( T + R(ic) )
   !                                  T           | T + R(ic1) - R(ic2) |   lm
   !
   !     with T = lattice vectors
   !
   !     An Ewald summation method devised by O.K. Andersen is used for l<5
   !     (see e.g. H.L. Skriver, "The LMTO method", Springer 1984).
   !     (The real-space function G can be calculated with gfunction.f)
   !

   SUBROUTINE structureconstant(structconst, cell, hybinp, atoms, kpts, mpi)

      USE m_juDFT
      USE m_types
      USE m_constants
      USE m_rorder, ONLY: rorderp, rorderpf
      use m_ylm
      IMPLICIT NONE

      TYPE(t_mpi), INTENT(IN)   :: mpi

      TYPE(t_hybinp), INTENT(IN) :: hybinp

      TYPE(t_cell), INTENT(IN)   :: cell

      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_kpts), INTENT(IN)    :: kpts
      ! - scalars -

      ! - arrays -
      COMPLEX, INTENT(INOUT)   ::  structconst(:, :, :, :)

      ! - local scalars -
      INTEGER                   ::  i, ic1, ic2, lm, ikpt, l, ishell, nshell
      INTEGER                   ::  m
      INTEGER                   ::  nptsh, maxl

      REAL                      ::  rad, rrad, rdum
      REAL                      ::  a, a1, aa
      REAL                      ::  pref, rexp
      REAL                      ::  scale

      COMPLEX                   ::  cdum, cexp

      LOGICAL, SAVE          ::  first = .TRUE.
      ! - local arrays -
      INTEGER                   ::  conv(0:2*hybinp%lexp)
      INTEGER, ALLOCATABLE     ::  pnt(:), ptsh(:, :)

      REAL                      ::  rc(3), ra(3), k(3), ki(3), ka(3)
      REAL                      ::  convpar(0:2*hybinp%lexp), g(0:2*hybinp%lexp)
      REAL, ALLOCATABLE     ::  radsh(:)

      COMPLEX                   ::  y((2*hybinp%lexp + 1)**2)
      COMPLEX                   ::  shlp((2*hybinp%lexp + 1)**2, kpts%nkpt)
      REAL, PARAMETER           :: CONVPARAM = 1e-18
      ! Do some additional shells ( real-space and Fourier-space sum )
      INTEGER, PARAMETER        :: ADDSHELL1 = 40
      INTEGER, PARAMETER        :: ADDSHELL2 = 0

      call timestart("calc struc_const.")

      IF (mpi%irank /= 0) first = .FALSE.

      rdum = cell%vol**(1.0/3) ! define "average lattice parameter"

      ! ewaldlambda = ewaldscale
      scale = hybinp%ewaldlambda/rdum

      !       lambda = ewaldlambda / rdum

      pref = fpi_const/(scale**3*cell%vol)

      DO l = 0, 2*hybinp%lexp
         convpar(l) = CONVPARAM/scale**(l + 1)
      END DO

      IF (first) THEN
         WRITE (oUnit, '(//A)') '### subroutine: structureconstant ###'
         WRITE (oUnit, '(/A)') 'Real-space sum:'
      END IF

      !
      !     Determine cutoff radii for real-space and Fourier-space summation
      ! (1) real space
      call timestart("determine cutoff radii")
      a = 1
1     rexp = EXP(-a)
      g(0) = rexp/a*(1 + a*11/16*(1 + a*3/11*(1 + a/9)))
      g(1) = rexp/a**2*(1 + a*(1 + a/2*(1 + a*7/24*(1 + a/7))))
      g(2) = rexp/a**3*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a*3/16 &
                                                          *(1 + a/9))))))
      g(3) = rexp/a**4*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                   *(1 + a/8)))))))
      g(4) = rexp/a**5*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                   *(1 + a/7*(1 + a/8*(1 + a/10)))))))))
      g(5) = rexp/a**6*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                   *(1 + a/7*(1 + a/8*(1 + a/9*(1 + a/10))))))))))
      g(6) = rexp/a**7*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                   *(1 + a/7*(1 + a/8*(1 + a/9*(1 + a/10*(1 + a/11*(1 + a/12))))))))))))
      g(7) = rexp/a**8*(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6 &
                                                                   *(1 + a/7*(1 + a/8*(1 + a/9*(1 + a/10*(1 + a/11*(1 + a/12*(1 + a/13)))))))))))))
      DO l = 8, 2*hybinp%lexp
         g(l) = a**(-l - 1)
      END DO
      IF (ANY(g > convpar/10)) THEN ! one digit more accuracy for real-space sum
         a = a + 1
         GOTO 1
      END IF
      rad = a/scale
      call timestop("determine cutoff radii")

      ! (2) Fourier space
      call timestart("fourier space")
      a = 1
2     aa = (1 + a**2)**(-1)
      g(0) = pref*aa**4/a**2
      g(1) = pref*aa**4/a
      g(2) = pref*aa**5/3
      g(3) = pref*aa**5*a/15
      g(4) = pref*aa**6*a**2/105
      g(5) = pref*aa**6*a**3/945
      g(6) = pref*aa**7*a**4/10395
      g(7) = pref*aa**7*a**5/135135
      IF (ANY(g > convpar)) THEN
         a = a + 1
         GOTO 2
      END IF
      rrad = a*scale
      call timestop("fourier space")

      IF (first) THEN
         WRITE (oUnit, '(/A,2F10.5)') 'Cutoff radii: ', rad, rrad
         WRITE (oUnit, '(/A)') 'Real-space sum'
      END IF

      !
      !     Determine atomic shells
      call timestart("determine atomic shell")
      CALL getshells(ptsh, nptsh, radsh, nshell, rad, cell%amat, first)
      call timestop("determine atomic shell")

      allocate (pnt(nptsh))
      structconst = 0

      !
      !     Real-space sum
      !
      call timestart("realspace sum")
      DO ic2 = 1, atoms%nat
         DO ic1 = 1, atoms%nat
            IF (ic2 /= 1 .AND. ic1 == ic2) CYCLE
            rc = MATMUL(cell%amat, (atoms%taual(:, ic2) - atoms%taual(:, ic1)))
            DO i = 1, nptsh
               ra = MATMUL(cell%amat, ptsh(:, i)) + rc
               a = norm2(ra)
               radsh(i) = a
            END DO
            call timestart("rorderpf")
            CALL rorderpf(pnt, radsh, nptsh, MAX(0, INT(LOG(nptsh*0.001)/LOG(2.0))))
            call timestop("rorderpf")
            ptsh = ptsh(:, pnt)
            radsh = radsh(pnt)
            maxl = 2*hybinp%lexp
            a1 = HUGE(a1)  ! stupid initial value
            ishell = 1
            conv = HUGE(i)
            shlp = 0
            DO i = 1, nptsh
               IF (ALL(conv /= HUGE(i))) EXIT
               IF (i /= 1) THEN
                  IF (ABS(radsh(i) - radsh(i - 1)) > 1e-10) ishell = ishell + 1
               ENDIF
               ra = MATMUL(cell%amat, ptsh(:, i)) + rc
               a = scale*norm2(ra)
               IF (abs(a) < 1e-12) THEN
                  CYCLE
               ELSE IF (ABS(a - a1) > 1e-10) THEN
                  a1 = a
                  rexp = EXP(-a)
                  IF (ishell <= conv(0)) g(0) = rexp/a &
                                                *(1 + a*11/16*(1 + a*3/11*(1 + a/9)))
                  IF (ishell <= conv(1)) g(1) = rexp/a**2 &
                                                *(1 + a*(1 + a/2*(1 + a*7/24*(1 + a/7))))
                  IF (ishell <= conv(2)) g(2) = rexp/a**3 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a*3/16*(1 + a/9))))))
                  IF (ishell <= conv(3)) g(3) = rexp/a**4 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/8)))))))
                  IF (ishell <= conv(4)) g(4) = rexp/a**5 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8 &
                                                                                                               *(1 + a/10)))))))))
                  IF (ishell <= conv(5)) g(5) = rexp/a**6 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8*(1 + a/9 &
                                                                                                                        *(1 + a/10))))))))))
                  IF (ishell <= conv(6)) g(6) = rexp/a**7 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8*(1 + a/9 &
                                                                                                                        *(1 + a/10*(1 + a/11*(1 + a/12))))))))))))
                  IF (ishell <= conv(7)) g(7) = rexp/a**8 &
                                                *(1 + a*(1 + a/2*(1 + a/3*(1 + a/4*(1 + a/5*(1 + a/6*(1 + a/7*(1 + a/8*(1 + a/9 &
                                                                                                                        *(1 + a/10*(1 + a/11*(1 + a/12*(1 + a/13)))))))))))))
                  DO l = 8, maxl
                     IF (ishell <= conv(l)) g(l) = a**(-l - 1)
                  END DO
                  DO l = 0, maxl
                     IF (conv(l) == HUGE(i) .AND. g(l) < convpar(l)/10) conv(l) = ishell + ADDSHELL1
                  END DO
               END IF
               IF (ishell > conv(maxl) .AND. maxl /= 0) maxl = maxl - 1
               call ylm4(maxl, ra, y)
               y = CONJG(y)
               DO ikpt = 1, kpts%nkpt
                  rdum = kpts%bk(1, ikpt)*ptsh(1, i) + kpts%bk(2, ikpt)*ptsh(2, i) + kpts%bk(3, ikpt)*ptsh(3, i)
                  cexp = EXP(CMPLX(0.0, 1.0)*tpi_const*rdum)
                  lm = 0
                  DO l = 0, maxl
                     IF (ishell <= conv(l)) THEN
                        cdum = cexp*g(l)
                        DO M = -l, l
                           lm = lm + 1
                           shlp(lm, ikpt) = shlp(lm, ikpt) + cdum*y(lm)
                        END DO
                     ELSE
                        lm = lm + 2*l + 1
                     END IF
                  END DO
               END DO
            END DO
            structconst(:, ic1, ic2, :) = shlp
         END DO
      END DO
      call timestop("realspace sum")

      deallocate (ptsh, radsh)

      IF (first) WRITE (oUnit, '(/A)') 'Fourier-space sum'

      !
      !     Determine reciprocal shells
      !
      call timestart("determince reciproc. shell")
      CALL getshells(ptsh, nptsh, radsh, nshell, rrad, cell%bmat, first)
      call timestop("determince reciproc. shell")
      ! minimum nonzero reciprocal-shell radius (needed in routines concerning the non-local hartree-fock exchange)
      !hybinp%radshmin = radsh(2)
      !
      !     Fourier-space sum
      !
      call timestart("fourierspace sum")
      DO ikpt = 1, kpts%nkpt
         k = kpts%bk(:, ikpt)
         maxl = MIN(7, hybinp%lexp*2)
         ishell = 1
         conv = HUGE(i)
         DO i = 1, nptsh
            IF (i > 1) THEN
               IF (ABS(radsh(i) - radsh(i - 1)) > 1e-10) ishell = ishell + 1
            ENDIF
            ki = ptsh(:, i) + k - NINT(k) ! -nint(...) transforms to Wigner-Seitz cell ( i.e. -0.5 <= x,y,z < 0.5 )
            ka = MATMUL(ki, cell%bmat)
            a = norm2(ka)/scale
            aa = (1 + a**2)**(-1)
            IF (ABS(a - a1) > 1e-10) THEN
               a1 = a
               IF (abs(a) < 1e-12) THEN
                  g(0) = pref*(-4)
                  g(1) = 0
               ELSE
                  IF (ishell <= conv(0)) g(0) = pref*aa**4/a**2
                  IF (ishell <= conv(1)) g(1) = pref*aa**4/a
               END IF
               IF (ishell <= conv(2)) g(2) = pref*aa**5/3
               IF (ishell <= conv(3)) g(3) = pref*aa**5*a/15
               IF (ishell <= conv(4)) g(4) = pref*aa**6*a**2/105
               IF (ishell <= conv(5)) g(5) = pref*aa**6*a**3/945
               IF (ishell <= conv(6)) g(6) = pref*aa**7*a**4/10395
               IF (ishell <= conv(7)) g(7) = pref*aa**7*a**5/135135
               IF (ishell > 1) THEN
                  DO l = 0, 7
                     IF (conv(l) == HUGE(i) .AND. g(l) < convpar(l)) conv(l) = ishell + ADDSHELL2
                  END DO
               END IF
            END IF

            IF (ishell > conv(maxl) .AND. maxl /= 0) maxl = maxl - 1
            call ylm4(maxl, ka, y)
            IF (norm2(ka(:)) .LT. 1.0e-16) y(2:(maxl + 1)**2) = CMPLX(0.0, 0.0)
            cdum = 1.0
            lm = 0
            DO l = 0, maxl
               IF (ishell <= conv(l)) THEN
                  DO M = -l, l
                     lm = lm + 1
                     y(lm) = CONJG(y(lm))*cdum*g(l)
                  END DO
               ELSE
                  y(lm + 1:lm + 2*l + 1) = 0
                  lm = lm + 2*l + 1
               END IF
               cdum = cdum*CMPLX(0.0, 1.0)
            END DO
            DO ic2 = 1, atoms%nat
               DO ic1 = 1, atoms%nat
                  IF (ic2 /= 1 .AND. ic1 == ic2) CYCLE
                  cexp = EXP(CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT(ki, atoms%taual(:, ic1) - atoms%taual(:, ic2)))
                  DO lm = 1, (maxl + 1)**2
                     structconst(lm, ic1, ic2, ikpt) = structconst(lm, ic1, ic2, ikpt) + cexp*y(lm)
                  END DO
               END DO
            END DO
         END DO
      END DO
      call timestop("fourierspace sum")

      !
      !     Add contribution for l=0 to diagonal elements and rescale structure constants
      !
      structconst(1, 1, 1, :) = structconst(1, 1, 1, :) - 5.0/16/SQRT(fpi_const)
      DO i = 2, atoms%nat
         structconst(:, i, i, :) = structconst(:, 1, 1, :)
      END DO
      DO l = 0, 2*hybinp%lexp
         structconst(l**2 + 1:(l + 1)**2, :, :, :) = structconst(l**2 + 1:(l + 1)**2, :, :, :)*scale**(l + 1)
      END DO

      rad = (cell%vol*3/4/pi_const)**(1.0/3) ! Wigner-Seitz radius (rad is recycled)

      !     Calculate accuracy of Gamma-decomposition
      IF (ALL(abs(kpts%bk) > 1e-12)) THEN
         a = 1e30 ! ikpt = index of shortest non-zero k-point
         DO i = 2, kpts%nkpt
            rdum = SUM(MATMUL(kpts%bk(:, i), cell%bmat)**2)
            IF (rdum < a) THEN
               ikpt = i
               a = rdum
            END IF
         END DO
         rdum = norm2(MATMUL(kpts%bk(:, ikpt), cell%bmat))
         a = 0
         DO ic2 = 1, atoms%nat
            DO ic1 = 1, MAX(1, ic2 - 1)
               a = a + ABS(structconst(1, ic1, ic2, ikpt) - &
                           (structconst(1, ic1, ic2, 1) + SQRT(fpi_const)/cell%vol/rdum**2* &
                            EXP(-CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT( &
                                kpts%bk(:, ikpt), atoms%taual(:, ic2) - atoms%taual(:, ic1)))))**2
            END DO
         END DO
         a = SQRT(a/atoms%nat**2)
         aa = SQRT(SUM(ABS(structconst(1, :, :, ikpt))**2)/atoms%nat**2)
         IF (first) WRITE (oUnit, '(/A,F8.5,A,F8.5,A)') 'Accuracy of Gamma-decomposition (structureconstant):', a, ' (abs)', a/aa, ' (rel)'
      ENDIF
      deallocate (ptsh, radsh)

      first = .FALSE.

      call timestop("calc struc_const.")
   END SUBROUTINE structureconstant

   !     -----------------

   !     Determines all shells of the crystal defined by lat and vol with radii smaller than rad.
   !     The lattice points (number = nptsh) are stored in ptsh, their corresponding lengths (shell radii) in radsh.

   SUBROUTINE getshells(ptsh, nptsh, radsh, nshell, rad, lat, lwrite)

      USE m_juDFT
      USE m_constants
      USE m_rorder, ONLY: rorderpf

      IMPLICIT NONE

      LOGICAL, INTENT(IN)    :: lwrite
      INTEGER, INTENT(INOUT)   :: nptsh, nshell
      INTEGER, ALLOCATABLE   :: ptsh(:, :)
      REAL, ALLOCATABLE   :: radsh(:)
      REAL, INTENT(IN)    :: rad, lat(:, :)
      REAL                    :: r(3), rdum
      INTEGER, ALLOCATABLE   :: pnt(:)
      INTEGER                 :: n, i, ix, iy, iz, ok
      LOGICAL                 :: found
      INTEGER, ALLOCATABLE   :: ihelp(:, :)
      REAL, ALLOCATABLE   :: rhelp(:)

      allocate (ptsh(3, 100000), radsh(100000), stat=ok)
      IF (ok /= 0) call judft_error('getshells: failure allocation ptsh/radsh')

      ptsh = 0
      radsh = 0

      n = 0
      i = 0
      DO
         found = .FALSE.
         DO ix = -n, n
            DO iy = -(n - ABS(ix)), n - ABS(ix)
               iz = n - ABS(ix) - ABS(iy)
1              r = ix*lat(:, 1) + iy*lat(:, 2) + iz*lat(:, 3)
               rdum = SUM(r**2)
               IF (rdum < rad**2) THEN
                  found = .TRUE.
                  i = i + 1
                  IF (i > SIZE(radsh)) THEN
                     allocate (rhelp(SIZE(radsh)), ihelp(3, SIZE(ptsh, 2)), stat=ok)
                     IF (ok /= 0) call judft_error('getshells: failure allocation rhelp/ihelp')
                     rhelp = radsh
                     ihelp = ptsh
                     deallocate (radsh, ptsh)
                     allocate (radsh(SIZE(rhelp) + 100000), ptsh(3, SIZE(ihelp, 2) + 100000), stat=ok)
                     IF (ok /= 0) call judft_error('getshells: failure re-allocation ptsh/radsh')
                     radsh(1:SIZE(rhelp)) = rhelp
                     ptsh(:, 1:SIZE(ihelp, 2)) = ihelp
                     deallocate (rhelp, ihelp)
                  END IF
                  ptsh(:, i) = [ix, iy, iz]
                  radsh(i) = SQRT(rdum)
               END IF
               IF (iz > 0) THEN
                  iz = -iz
                  GOTO 1
               END IF
            END DO
         END DO
         IF (.NOT. found) EXIT
         n = n + 1
      END DO
      nptsh = i

      allocate (pnt(nptsh))

      !reallocate radsh ptsh
      allocate (rhelp(nptsh), ihelp(3, nptsh))
      rhelp = radsh(1:nptsh)
      ihelp = ptsh(:, 1:nptsh)
      deallocate (radsh, ptsh)
      allocate (radsh(nptsh), ptsh(3, nptsh))
      radsh = rhelp
      ptsh = ihelp
      deallocate (rhelp, ihelp)

      CALL rorderpf(pnt, radsh, nptsh, MAX(0, INT(LOG(nptsh*0.001)/LOG(2.0))))
      radsh = radsh(pnt)
      ptsh = ptsh(:, pnt)
      nshell = 1
      DO i = 2, nptsh
         IF (radsh(i) - radsh(i - 1) > 1e-10) nshell = nshell + 1
      END DO

      IF (lwrite) &
         WRITE (oUnit, '(A,F10.5,A,I7,A,I5,A)') &
         '  Sphere of radius', rad, ' contains', &
         nptsh, ' lattice points and', nshell, ' shells.'

   END SUBROUTINE getshells

   !     ---------

   !     Returns a list of (k+G) vector lengths in qnrm(1:nqnrm) and the corresponding pointer pqnrm(1:ngpt(ikpt),ikpt)
   SUBROUTINE getnorm(kpts, gpt, ngpt, pgpt, qnrm, nqnrm, pqnrm, cell)
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts

      INTEGER, INTENT(IN)  :: ngpt(:), gpt(:, :), pgpt(:, :)!(dim,kpts%nkpt)
      REAL, ALLOCATABLE :: qnrm(:), help(:)
      INTEGER, INTENT(INOUT) :: nqnrm
      INTEGER, ALLOCATABLE :: pqnrm(:, :)
      INTEGER               :: i, j, ikpt, igpt, igptp
      REAL                  :: q(3), qnorm

      allocate (qnrm(MAXVAL(ngpt)*kpts%nkpt), pqnrm(MAXVAL(ngpt), kpts%nkpt))
      i = 0
      DO ikpt = 1, kpts%nkpt
         igptloop: DO igpt = 1, ngpt(ikpt)
            igptp = pgpt(igpt, ikpt)
            IF (igptp == 0) call judft_error('getnorm: zero pointer (bug?)')
            q = MATMUL(kpts%bk(:, ikpt) + gpt(:, igptp), cell%bmat)
            qnorm = norm2(q)
            DO j = 1, i
               IF (ABS(qnrm(j) - qnorm) < 1e-12) THEN
                  pqnrm(igpt, ikpt) = j
                  CYCLE igptloop
               END IF
            END DO
            i = i + 1
            qnrm(i) = qnorm
            pqnrm(igpt, ikpt) = i
         END DO igptloop
      END DO
      nqnrm = i

      allocate (help(nqnrm))
      help(1:nqnrm) = qnrm(1:nqnrm)
      deallocate (qnrm)
      allocate (qnrm(1:nqnrm))
      qnrm = help

   END SUBROUTINE getnorm

   FUNCTION sphbessel_integral(atoms, itype, qnrm, nqnrm, iqnrm1, iqnrm2, l, hybinp, &
                               sphbes0, l_warnin, l_warnout)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN)  :: itype, nqnrm, iqnrm1, iqnrm2, l
      REAL, INTENT(IN)  :: qnrm(nqnrm), sphbes0(-1:hybinp%lexp + 2, atoms%ntype, nqnrm)
      LOGICAL, INTENT(IN), OPTIONAL  ::  l_warnin
      LOGICAL, INTENT(INOUT), OPTIONAL  ::  l_warnout
      REAL                  :: sphbessel_integral
      REAL                  :: q1, q2, dq, s, sb01, sb11, sb21, sb31, sb02, sb12
      REAL                  :: sb22, sb32, a1, a2, da, b1, b2, db, c1, c2, dc, r1, r2
      LOGICAL               :: l_warn, l_warned

      IF (PRESENT(l_warnin)) THEN
         l_warn = l_warnin
      ELSE
         l_warn = .TRUE.
      END IF
      l_warned = .FALSE.

      q1 = qnrm(iqnrm1)
      q2 = qnrm(iqnrm2)
      s = atoms%rmt(itype)
      IF (abs(q1) < 1e-12 .AND. abs(q2) < 1e-12) THEN
         IF (l > 0) THEN
            sphbessel_integral = 0
         ELSE
            sphbessel_integral = 2*s**5/15
         ENDIF
      ELSE IF (abs(q1) < 1e-12 .OR. abs(q2) < 1e-12) THEN
         IF (l > 0) THEN
            sphbessel_integral = 0
         ELSE IF (abs(q1) < 1e-12) THEN
            sphbessel_integral = s**3/(3*q2**2)*(q2*s*sphbes0(1, itype, iqnrm2) &
                                                 + sphbes0(2, itype, iqnrm2))
         ELSE
            sphbessel_integral = s**3/(3*q1**2)*(q1*s*sphbes0(1, itype, iqnrm1) &
                                                 + sphbes0(2, itype, iqnrm1))
         ENDIF
      ELSE IF (abs(q1 - q2) < 1e-12) THEN
         sphbessel_integral = s**3/(2*q1**2)*((2*l + 3)*sphbes0(l + 1, itype, iqnrm1)**2 - &
                                              (2*l + 1)*sphbes0(l, itype, iqnrm1)*sphbes0(l + 2, itype, iqnrm1))
      ELSE ! We use either if two fromulas that are stable for high and small q1/q2 respectively
         sb01 = sphbes0(l - 1, itype, iqnrm1)
         sb11 = sphbes0(l, itype, iqnrm1)
         sb21 = sphbes0(l + 1, itype, iqnrm1)
         sb31 = sphbes0(l + 2, itype, iqnrm1)
         sb02 = sphbes0(l - 1, itype, iqnrm2)
         sb12 = sphbes0(l, itype, iqnrm2)
         sb22 = sphbes0(l + 1, itype, iqnrm2)
         sb32 = sphbes0(l + 2, itype, iqnrm2)
         dq = q1**2 - q2**2
         a1 = q2/q1*sb21*sb02
         a2 = q1/q2*sb22*sb01
         da = a1 - a2
         b1 = sb31*sb12
         b2 = sb32*sb11
         db = b1 - b2
         c1 = sb21*sb22/(q1*q2)
         c2 = db/dq*(2*l + 1)/(2*l + 3)
         dc = c1 + c2
         r1 = ABS(da/a1)
         r2 = MIN(ABS(db/b1), ABS(dc/c1))
         ! Ensure numerical stability. If both formulas are not sufficiently stable, the program stops.
         IF (r1 > r2) THEN
            IF (r1 < 1e-6 .AND. l_warn) THEN
               WRITE (oUnit, '(A,E12.5,A,E12.5,A)') 'sphbessel_integral: Warning! Formula One possibly unstable. Ratios:', &
                  r1, '(', r2, ')'
               WRITE (oUnit, '(A,2F15.10,I4)') '                    Current qnorms and atom type:', q1, q2, itype
               l_warned = .TRUE.
            END IF
            sphbessel_integral = s**3/dq*da
         ELSE
            IF (r2 < 1e-6 .AND. l_warn) THEN
               WRITE (oUnit, '(A,E13.5,A,E13.5,A)') 'sphbessel_integral: Warning! Formula Two possibly unstable. Ratios:', &
                  r2, '(', r1, ')'
               WRITE (oUnit, '(A,2F15.10,I4)') '                    Current qnorms and atom type:', &
                  q1, q2, itype
               l_warned = .TRUE.
            END IF
            sphbessel_integral = s**3*dc
         END IF
      END IF

      IF (PRESENT(l_warnout)) l_warnout = l_warned

   END FUNCTION sphbessel_integral

   subroutine apply_inverse_olaps(mpdata, atoms, cell, hybdat, sym, kpts, ikpt, coulomb)
      USE m_olap, ONLY: olap_pw
      USE m_types
      use m_judft
      implicit none
      type(t_mpdata), intent(in) :: mpdata
      type(t_atoms), intent(in)  :: atoms
      type(t_cell), intent(in)   :: cell
      type(t_hybdat), intent(in) :: hybdat
      type(t_sym), intent(in)    :: sym
      type(t_kpts), intent(in)   :: kpts
      type(t_mat), intent(inout) :: coulomb
      integer, intent(in)        :: ikpt

      type(t_mat)     :: olap, coulhlp, coul_submtx
      integer         :: nbasm

      call timestart("solve olap linear eq. sys")
      nbasm = hybdat%nbasp + mpdata%n_g(ikpt)
      CALL olap%alloc(sym%invs, mpdata%n_g(ikpt), mpdata%n_g(ikpt), 0.0)
      !calculate IR overlap-matrix
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(ikpt), ikpt)), mpdata%n_g(ikpt), atoms, cell)

      ! perform O^-1 * coulhlp%data_r(hybdat%nbasp + 1:, :) = x
      ! rewritten as O * x = C

      call coul_submtx%alloc(sym%invs, mpdata%n_g(ikpt), nbasm)
      if (coul_submtx%l_real) then
         coul_submtx%data_r = real(coulomb%data_c(hybdat%nbasp + 1:, :))
      else
         coul_submtx%data_c = coulomb%data_c(hybdat%nbasp + 1:, :)
      endif

      call olap%linear_problem(coul_submtx)

      if (coul_submtx%l_real) then
         coulomb%data_c(hybdat%nbasp + 1:, :) = coul_submtx%data_r
         coul_submtx%data_r = real(transpose(coulomb%data_c(:, hybdat%nbasp + 1:)))
      else
         coulomb%data_c(hybdat%nbasp + 1:, :) = coul_submtx%data_c
         coul_submtx%data_c = conjg(transpose(coulomb%data_c(:, hybdat%nbasp + 1:)))
      endif

      ! perform  coulomb%data_r(hybdat%nbasp + 1:, :) * O^-1  = X
      ! rewritten as O^T * x^T = C^T

      ! reload O, since the solver destroys it.
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(ikpt), ikpt)), mpdata%n_g(ikpt), atoms, cell)
      ! Notice O = O^T since it's symmetric
      call olap%linear_problem(coul_submtx)

      if (coul_submtx%l_real) then
         coulomb%data_c(:, hybdat%nbasp + 1:) = transpose(coul_submtx%data_r)
      else
         coulomb%data_c(:, hybdat%nbasp + 1:) = conjg(transpose(coul_submtx%data_c))
      endif

      call coul_submtx%free()
      call olap%free()
      call timestop("solve olap linear eq. sys")
   end subroutine apply_inverse_olaps

   subroutine loop_over_interst(fi, hybdat, mpdata, structconst, sphbesmoment, moment, moment2, &
                                qnrm, facc, gmat, integral, olap, pqnrm, pgptm1, ngptm1, ikpt, coulmat)
      use m_types
      use m_juDFT
      use m_ylm, only: ylm4
      use m_constants, only: fpi_const, tpi_const
      USE m_trafo, ONLY: symmetrize
      use m_calc_l_m_from_lm
      implicit none

      type(t_fleurinput), intent(in)    :: fi
      type(t_hybdat), intent(in)        :: hybdat
      type(t_mpdata), intent(in)        :: mpdata
      REAL, intent(in)                  :: sphbesmoment(0:, :, :), qnrm(:), facC(-1:), gmat(:, :), moment(:, 0:, :), moment2(:, :)
      real, intent(in)                  :: integral(:, 0:, :, :), olap(:, 0:, :, :)
      integer, intent(in)               :: ikpt, ngptm1(:), pqnrm(:, :), pgptm1(:, :)
      complex, intent(in)               :: structconst(:, :, :, :)
      type(t_mat), intent(inout)        :: coulmat

      integer  :: igpt0, igpt, igptp, iqnrm, niter
      integer  :: ix, iy, ic, itype, lm, l, m, itype1, ic1, l1, m1, lm1
      integer  :: l2, m2, lm2, n, i, idum, iatm, j_type, j_l, iy_start, j_m, j_lm
      real     :: q(3), qnorm, svol
      COMPLEX  :: y((fi%hybinp%lexp + 1)**2), y1((fi%hybinp%lexp + 1)**2), y2((fi%hybinp%lexp + 1)**2)
      complex  :: csum, csumf(9), cdum, cexp
      integer, allocatable :: lm_arr(:), ic_arr(:)

      coulmat%data_c(:hybdat%nbasp,hybdat%nbasp+1:) = 0
      svol = SQRT(fi%cell%vol)
      ! start to loop over interstitial plane waves
      DO igpt0 = 1, ngptm1(ikpt) !1,ngptm1(ikpt)
         igpt = pgptm1(igpt0, ikpt)
         igptp = mpdata%gptm_ptr(igpt, ikpt)
         ix = hybdat%nbasp + igpt
         q = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)
         qnorm = norm2(q)
         iqnrm = pqnrm(igpt, ikpt)
         IF (ABS(qnrm(iqnrm) - qnorm) > 1e-12) then
            call judft_error('coulombmatrix: qnorm does not equal corresponding & element in qnrm (bug?)') ! We shouldn't st op here!
         endif

         call ylm4(2, MATMUL(fi%kpts%bk(:, fi%kpts%nkpt), fi%cell%bmat), y1)
         call ylm4(2, MATMUL(mpdata%g(:, igptp), fi%cell%bmat), y2)
         call ylm4(fi%hybinp%lexp, q, y)
         y1 = CONJG(y1); y2 = CONJG(y2); y = CONJG(y)

         ! this unrolls the do ic=1,atoms%nat{do lm=1,..{}} 
         call collapse_ic_and_lm_loop(fi%atoms, fi%hybinp%lcutm1, niter, ic_arr, lm_arr)

         !$OMP PARALLEL DO default(none) &
         !$OMP private(ic, lm, itype, l, m, csum, csumf, ic1, itype1, cexp, lm1, l2, cdum, m2, lm2, iy) &
         !$OMP private(j_m, j_type, iy_start, l1, m1) &
         !$OMP shared(ic_arr, lm_arr, fi, mpdata, olap, qnorm, moment, integral, hybdat, coulmat, svol) &
         !$OMP shared(moment2, ix, igpt, facc, structconst, y, y1, y2, gmat, iqnrm, sphbesmoment, ikpt) &
         !$OMP shared(igptp, niter) &
         !$OMP schedule(dynamic)
         do i = 1,niter 
            ic = ic_arr(i)
            lm = lm_arr(i)

            itype = fi%atoms%itype(ic)
            call calc_l_m_from_lm(lm, l, m)

            ! calculate sum over lm and centers for (2c) -> csum, csumf
            csum = 0
            csumf = 0
            do ic1 = 1, fi%atoms%nat
               itype1 = fi%atoms%itype(ic1)
               cexp = fpi_const*EXP(CMPLX(0.0, 1.0)*tpi_const &
                                    *(dot_PRODUCT(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%atoms%taual(:, ic1)) &
                                       - dot_PRODUCT(fi%kpts%bk(:, ikpt), fi%atoms%taual(:, ic))))

               do lm1 = 1, (fi%hybinp%lexp+1)**2
                  call calc_l_m_from_lm(lm1, l1, m1)
                  l2 = l + l1 ! for structconst
                  cdum = sphbesmoment(l1, itype1, iqnrm)*CMPLX(0.0, 1.0)**(l1)*cexp
                  m2 = M - m1              ! for structconst
                  lm2 = l2**2 + l2 + m2 + 1 !
                  csum = csum - (-1)**(m1 + l1)*gmat(lm1, lm)*y(lm1)*cdum*structconst(lm2, ic, ic1, ikpt)
               END DO

               ! add contribution of (2c) to csum and csumf coming from linear and quadratic orders of Y_lm*(G) / G * j_(l+1)(GS)
               IF (ikpt == 1 .AND. l <= 2) THEN
                  cexp = EXP(CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT(mpdata%g(:, igptp), fi%atoms%taual(:, ic1))) &
                           *gmat(lm, 1)*fpi_const/fi%cell%vol
                  csumf(lm) = csumf(lm) - cexp*SQRT(fpi_const)* &
                              CMPLX(0.0, 1.0)**l*sphbesmoment(0, itype1, iqnrm)/facC(l - 1)
                  IF (l == 0) THEN
                     IF (igpt /= 1) THEN
                        csum = csum - cexp*(sphbesmoment(0, itype1, iqnrm)*fi%atoms%rmt(itype1)**2 - &
                                             sphbesmoment(2, itype1, iqnrm)*2.0/3)/10
                     ELSE
                        csum = csum - cexp*fi%atoms%rmt(itype1)**5/30
                     END IF
                  ELSE IF (l == 1) THEN
                     csum = csum + cexp*CMPLX(0.0, 1.0)*SQRT(fpi_const) &
                              *sphbesmoment(1, itype1, iqnrm)*y(lm)/3
                  END IF
               END IF
            END DO

            ! add contribution of (2a) to csumf
            IF (ikpt == 1 .AND. igpt == 1 .AND. l <= 2) THEN
               csumf(lm) = csumf(lm) + (fpi_const)**2*CMPLX(0.0, 1.0)**l/facC(l)
            END IF

            ! finally define coulomb
            cdum = (fpi_const)**2*CMPLX(0.0, 1.0)**(l)*y(lm) &
                     *EXP(CMPLX(0.0, 1.0)*tpi_const &
                        *dot_PRODUCT(mpdata%g(:, igptp), fi%atoms%taual(:, ic)))

            !calclate iy_start on the fly for OpenMP
            iy_start = 0
            do iatm = 1, ic-1
               j_type = fi%atoms%itype(iatm)
               do j_l = 0,fi%hybinp%lcutm1(j_type)
                  iy_start = iy_start + mpdata%num_radbasfn(j_l, j_type) * (2*j_l+1)
               end do
            end do
            do j_lm = 1,lm-1
               call calc_l_m_from_lm(j_lm, j_l, j_m)
               iy_start = iy_start + mpdata%num_radbasfn(j_l, itype) 
            enddo

            DO n = 1, mpdata%num_radbasfn(l, itype)
               iy = iy_start + n

               IF (ikpt == 1 .AND. igpt == 1) THEN
                  IF (l == 0) coulmat%data_c(iy, ix) = &
                     -cdum*moment2(n, itype)/6/svol         ! (2a)
                  coulmat%data_c(iy, ix) = coulmat%data_c(iy, ix) &
                                                   + (-cdum/(2*l + 1)*integral(n, l, itype, iqnrm) & ! (2b)&
                                                      + csum*moment(n, l, itype))/svol          ! (2c)
               ELSE
                  coulmat%data_c(iy, ix) = &
                     (cdum*olap(n, l, itype, iqnrm)/qnorm**2 &  ! (2a)&
                        - cdum/(2*l + 1)*integral(n, l, itype, iqnrm) & ! (2b)&
                        + csum*moment(n, l, itype))/svol          ! (2c)

               END IF
            END DO
         END DO ! collapsed atom & lm loop (ic)
         !$OMP END PARALLEL DO
      END DO

      IF (fi%sym%invs) THEN
         CALL symmetrize(coulmat%data_c(:hybdat%nbasp,hybdat%nbasp+1:), hybdat%nbasp, mpdata%n_g(ikpt), 1, .FALSE., &
                         fi%atoms, fi%hybinp%lcutm1, maxval(fi%hybinp%lcutm1), mpdata%num_radbasfn, fi%sym)
      ENDIF
   endsubroutine loop_over_interst

   subroutine perform_double_g_loop(fi, hybdat, mpdata, sphbes0, carr2, ngptm1,pgptm1,pqnrm,qnrm, nqnrm, ikpt, coulomb)
      use m_juDFT
      use m_types
      use m_constants, only: tpi_const,fpi_const
      use omp_lib
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      integer, intent(in)               :: ikpt, ngptm1(:), pqnrm(:,:),pgptm1(:, :), nqnrm
      real, intent(in)                  :: qnrm(:), sphbes0(:, :, :)
      complex, intent(in)               :: carr2(:, :)
      !complex, intent(inout)            :: coulomb(:) ! only at ikpt
      type(t_mat), intent(inout)        :: coulomb

      integer :: igpt0, igpt1, igpt2, ix, iy, igptp1, igptp2, iqnrm1, iqnrm2
      integer :: ic, itype, lm, m, idum, l, i
      real    :: q1(3), q2(3)
      complex :: y1((fi%hybinp%lexp + 1)**2), y2((fi%hybinp%lexp + 1)**2)
      COMPLEX :: cexp1(fi%atoms%ntype)
      complex :: cdum, cdum1 
      logical :: ldum
      integer, parameter :: lock_size = 100
      integer(kind=omp_lock_kind) :: lock(0:lock_size-1)

      call timestart("double g-loop")

      ! create lock for race-condition in coulomb
      do i =0,lock_size-1
         call omp_init_lock(lock(i))
      enddo

      !$OMP PARALLEL DO default(none) &
      !$OMP private(igpt0, igpt1, igpt2, ix, igptp2, iqnrm2, q2, y2, iy,igptp1, iqnrm1, q1) &
      !$OMP private(y1, ic, itype, cexp1, lm, cdum, l, cdum1, m, idum, ldum) &
      !$OMP shared(coulomb, ngptm1, ikpt, pgptm1, hybdat, mpdata, pqnrm, fi) &
      !$OMP shared(lock, nqnrm, sphbes0, qnrm, carr2) &
      !$OMP schedule(dynamic)
      DO igpt0 = 1, ngptm1(ikpt)!1,ngptm1(ikpt)
         igpt2 = pgptm1(igpt0, ikpt)
         ix = hybdat%nbasp + igpt2
         igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
         iqnrm2 = pqnrm(igpt2, ikpt)
         q2 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp2), fi%cell%bmat)
         y2 = CONJG(carr2(:, igpt2))
         DO igpt1 = 1, igpt2
            iy = hybdat%nbasp + igpt1
            igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
            iqnrm1 = pqnrm(igpt1, ikpt)
            q1 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp1), fi%cell%bmat)
            y1 = carr2(:, igpt1)
            cexp1 = 0
            do ic = 1,fi%atoms%nat 
               itype = fi%atoms%itype(ic)
               cexp1(itype) = cexp1(itype) + &
                              EXP(CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT( &
                                    (mpdata%g(:, igptp2) - mpdata%g(:, igptp1)), fi%atoms%taual(:, ic)))
            ENDDO
            lm = 0
            cdum = 0
            DO l = 0, fi%hybinp%lexp
               cdum1 = 0
               DO itype = 1, fi%atoms%ntype
                  cdum1 = cdum1 + cexp1(itype)*sphbessel_integral( &
                           fi%atoms, itype, qnrm, nqnrm, &
                           iqnrm1, iqnrm2, l, fi%hybinp, &
                           sphbes0, .False., ldum) &
                           /(2*l + 1)
               END DO
               DO M = -l, l
                  lm = lm + 1
                  cdum = cdum + cdum1*y1(lm)*y2(lm)
               ENDDO
            ENDDO
            idum = ix*(ix - 1)/2 + iy
            call omp_set_lock(lock(modulo(idum,lock_size)))
            coulomb%data_c(iy,ix) = coulomb%data_c(iy,ix) + (fpi_const)**3*cdum/fi%cell%vol
            call omp_unset_lock(lock(modulo(idum,lock_size)))
         END DO
      END DO
      !$OMP END PARALLEL DO
      do i =0,lock_size-1
         call omp_destroy_lock(lock(i))
      enddo
      call timestop("double g-loop")
   end subroutine perform_double_g_loop

   subroutine collapse_ic_and_lm_loop(atoms, lcutm1, niter, ic_arr, lm_arr)
      use m_types
      implicit none 
      type(t_atoms), intent(in) :: atoms 
      integer, intent(in)       :: lcutm1(:)
      integer, intent(out)      :: niter 
      integer, intent(inout), allocatable :: ic_arr(:),  lm_arr(:) 

      integer :: ic, lm, itype 

      if(allocated(ic_arr)) deallocate(ic_arr)
      if(allocated(lm_arr)) deallocate(lm_arr)

      niter = 0
      do ic = 1, atoms%nat
         itype = atoms%itype(ic)
         do lm = 1,(lcutm1(itype)+1)**2
            niter = niter + 1
         enddo 
      enddo

      allocate( lm_arr(niter), ic_arr(niter))
      niter = 0
      do ic = 1, atoms%nat
         itype = atoms%itype(ic)
         do lm = 1,(lcutm1(itype)+1)**2
            niter = niter + 1
            lm_arr(niter)    = lm
            ic_arr(niter)    = ic
         enddo 
      enddo
   end subroutine collapse_ic_and_lm_loop
END MODULE m_coulombmatrix
