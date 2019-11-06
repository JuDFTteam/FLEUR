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

   SUBROUTINE coulombmatrix(mpi, atoms, kpts, cell, sym, mpbasis, hybrid, xcpot)

      USE m_types
      USE m_juDFT
      USE m_constants, ONLY: pi_const
      USE m_olap, ONLY: olap_pw
      use m_types_hybrid, only: gptnorm
      USE m_trafo, ONLY: symmetrize, bramat_trafo
      USE m_intgrf, ONLY: intgrf, intgrf_init
      use m_util, only: primitivef
      USE m_hsefunctional, ONLY: change_coulombmatrix
      USE m_wrapper
      USE m_io_hybrid
      use m_ylm
      use m_sphbes, only: sphbes
      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN)     :: xcpot
      TYPE(t_mpi), INTENT(IN)       :: mpi
      TYPE(t_mpbasis), intent(in)  :: mpbasis
      TYPE(t_hybrid), INTENT(INOUT) :: hybrid
      TYPE(t_sym), INTENT(IN)       :: sym
      TYPE(t_cell), INTENT(IN)      :: cell
      TYPE(t_kpts), INTENT(IN)      :: kpts
      TYPE(t_atoms), INTENT(IN)     :: atoms


      ! - local scalars -
      INTEGER                    :: inviop
      INTEGER                    :: nqnrm, iqnrm, iqnrm1, iqnrm2, iqnrmstart, iqnrmstep
      INTEGER                    :: itype, l, ix, iy, iy0, i, j, lm, l1, l2, m1, m2, ineq, idum, ikpt, ikpt0, ikpt1
      INTEGER                    :: lm1, lm2, itype1, itype2, ineq1, ineq2, n, n1, n2, ng
      INTEGER                    :: ic, ic1, ic2, ic3, ic4
      INTEGER                    :: igpt, igpt1, igpt2, igptp, igptp1, igptp2
      INTEGER                    :: isym, isym1, isym2, igpt0
      INTEGER                    :: ok
      INTEGER                    :: m
      INTEGER                    :: ikptmin, ikptmax, nkminmax
      INTEGER                    :: maxfac

      LOGICAL                    :: lsym

      REAL                       :: rdum, rdum1, rdum2
      REAL                       :: svol, qnorm, qnorm1, qnorm2, gnorm
      REAL                       :: fcoulfac
      REAL                       :: time1, time2

      COMPLEX                    :: cdum, cdum1, cexp, csum

      ! - local arrays -
      INTEGER                    :: g(3)
      INTEGER                    :: nbasm1(kpts%nkptf)
      INTEGER, ALLOCATABLE   :: pqnrm(:, :)
      INTEGER                    :: rrot(3, 3, sym%nsym), invrrot(3, 3, sym%nsym)
      INTEGER, ALLOCATABLE   :: iarr(:), POINTER(:, :, :, :)!,pointer(:,:,:)
      INTEGER                    :: igptmin(kpts%nkpt), igptmax(kpts%nkpt)
      INTEGER, ALLOCATABLE   :: nsym_gpt(:, :), sym_gpt(:, :, :)
      INTEGER                    :: nsym1(kpts%nkpt + 1), sym1(sym%nsym, kpts%nkpt + 1)

      INTEGER, ALLOCATABLE   ::  ngptm1(:)
      INTEGER, ALLOCATABLE   ::  pgptm1(:,:)

      LOGICAL                    :: calc_mt(kpts%nkpt)

      REAL                       :: q(3), q1(3), q2(3)
      REAL                       :: integrand(atoms%jmtd), primf1(atoms%jmtd), primf2(atoms%jmtd)
      REAL                       :: mat(maxval(mpbasis%num_rad_bas_fun)*(maxval(mpbasis%num_rad_bas_fun) + 1)/2)
      REAL                       :: moment(maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), atoms%ntype), &
                                    moment2(maxval(mpbasis%num_rad_bas_fun), atoms%ntype)
      REAL                       :: sphbes_var(atoms%jmtd, 0:maxval(hybrid%lcutm1))
      REAL                       :: sphbesmoment1(atoms%jmtd, 0:maxval(hybrid%lcutm1))
      REAL                       :: rarr(0:hybrid%lexp + 1), rarr1(0:maxval(hybrid%lcutm1))
      REAL, ALLOCATABLE   :: gmat(:, :), qnrm(:)
      REAL, ALLOCATABLE   :: sphbesmoment(:, :, :)
      REAL, ALLOCATABLE   :: sphbes0(:, :, :)
      REAL, ALLOCATABLE   :: olap(:, :, :, :), integral(:, :, :, :)
      REAL, ALLOCATABLE   :: gridf(:, :)
      REAL                       :: facA(0:MAX(2*atoms%lmaxd + maxval(hybrid%lcutm1) + 1, 4*MAX(maxval(hybrid%lcutm1), hybrid%lexp) + 1))
      REAL                       :: facB(0:MAX(2*atoms%lmaxd + maxval(hybrid%lcutm1) + 1, 4*MAX(maxval(hybrid%lcutm1), hybrid%lexp) + 1))
      REAL                       :: facC(-1:MAX(2*atoms%lmaxd + maxval(hybrid%lcutm1) + 1, 4*MAX(maxval(hybrid%lcutm1), hybrid%lexp) + 1))

      COMPLEX     :: cexp1(atoms%ntype), csumf(9)
      COMPLEX     :: structconst((2*hybrid%lexp + 1)**2, atoms%nat, atoms%nat, kpts%nkpt)             ! nw = 1
      COMPLEX     :: y((hybrid%lexp + 1)**2), y1((hybrid%lexp + 1)**2), y2((hybrid%lexp + 1)**2)
      COMPLEX     :: dwgn(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1), sym%nsym)
      COMPLEX, ALLOCATABLE   :: smat(:, :)
      COMPLEX, ALLOCATABLE   :: coulmat(:, :)
      COMPLEX, ALLOCATABLE   :: carr2(:, :), carr2a(:, :), carr2b(:, :)
      COMPLEX, ALLOCATABLE   :: structconst1(:, :)
      REAL, ALLOCATABLE   :: coulomb_mt1(:, :, :, :, :)

      !REAL       , ALLOCATABLE   :: coulomb(:,:) !At the moment we always calculate a complex coulomb matrix
      REAL, ALLOCATABLE   :: coulomb_mt2_r(:, :, :, :, :), coulomb_mt3_r(:, :, :, :)
      REAL, ALLOCATABLE   :: coulomb_mtir_r(:, :, :), coulombp_mtir_r(:, :)
      COMPLEX, ALLOCATABLE   :: coulomb(:, :)
      COMPLEX, ALLOCATABLE   :: coulomb_mt2_c(:, :, :, :, :), coulomb_mt3_c(:, :, :, :)
      COMPLEX, ALLOCATABLE   :: coulomb_mtir_c(:, :, :), coulombp_mtir_c(:, :)

      INTEGER                    :: ishift, ishift1
      INTEGER                    :: iatom, iatom1
      INTEGER                    :: indx1, indx2, indx3, indx4
      LOGICAL                    :: l_warn, l_warned!.true.!.false.
      TYPE(t_mat)                :: olapm, coulhlp

      CALL timestart("Coulomb matrix setup")

      svol = SQRT(cell%vol)
      fcoulfac = 4*pi_const/cell%vol
      maxfac = MAX(2*atoms%lmaxd + maxval(hybrid%lcutm1) + 1, 4*MAX(maxval(hybrid%lcutm1), hybrid%lexp) + 1)

      facA(0) = 1                    !
      facB(0) = 1                    ! Define:
      facC(-1:0) = 1                    ! facA(i)    = i!
      DO i = 1, maxfac                       ! facB(i)   = sqrt(i!)
         facA(i) = facA(i - 1)*i            ! facC(i) = (2i+1)!!
         facB(i) = facB(i - 1)*SQRT(i*1.0) !
         facC(i) = facC(i - 1)*(2*i + 1)   !
      END DO

      CALL intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, gridf)

      nbasm1 = hybrid%nbasp + mpbasis%ngptm(:)

      !     Calculate the structure constant
      CALL structureconstant(structconst, cell, hybrid, atoms, kpts, mpi)

      IF (mpi%irank == 0) WRITE (6, '(//A)') '### subroutine: coulombmatrix ###'

      !
      !     Matrix allocation
      !

      call timestart("coulomb allocation")
      IF (ALLOCATED(coulomb)) deallocate(coulomb)

      allocate(coulomb(hybrid%maxbasm1*(hybrid%maxbasm1 + 1)/2, kpts%nkpt), stat=ok)
      IF (ok /= 0) call judft_error('coulombmatrix: failure allocation coulomb matrix')
      coulomb = 0
      call timestop("coulomb allocation")

      IF (mpi%irank == 0) WRITE (6, '(/A,F6.1," MB")') 'Size of coulomb matrix:', 16.0/1048576*SIZE(coulomb)

      !     Generate Symmetry:
      !     Reduce list of g-Points so that only one of each symm-equivalent is calculated

      IF (mpi%irank == 0) WRITE (6, '(/A)', advance='no') 'Setup for symmetry...'
      CALL cpu_TIME(time1)
      ! calculate rotations in reciprocal space
      DO isym = 1, sym%nsym
         IF (isym <= sym%nop) THEN
            inviop = sym%invtab(isym)
            rrot(:, :, isym) = TRANSPOSE(sym%mrot(:, :, inviop))
            DO l = 0, maxval(hybrid%lcutm1)
               dwgn(:, :, l, isym) = TRANSPOSE(hybrid%d_wgn2(-maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), &
                                                             -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), l, isym))
            END DO
         ELSE
            inviop = isym - sym%nop
            rrot(:, :, isym) = -rrot(:, :, inviop)
            dwgn(:, :, :, isym) = dwgn(:, :, :, inviop)
            DO l = 0, maxval(hybrid%lcutm1)
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
      invrrot(:, :, :sym%nop) = rrot(:, :, sym%invtab)
      IF (sym%nsym > sym%nop) THEN
         invrrot(:, :, sym%nop + 1:) = rrot(:, :, sym%invtab + sym%nop)
      END IF

      ! Get symmetry operations that leave bk(:,ikpt) invariant -> sym1
      nsym1 = 0
      DO ikpt = 1, kpts%nkpt
         isym1 = 0
         DO isym = 1, sym%nsym
            ! temporary fix until bramat_trafo is correct
            ! for systems with symmetries including translations
            IF (isym > sym%nop) THEN
               isym2 = isym - sym%nop
            ELSE
               isym2 = isym
            END IF
            IF (ANY(abs(sym%tau(:, isym2)) > 1e-12)) CYCLE

            IF (ALL(ABS(MATMUL(rrot(:, :, isym), kpts%bk(:, ikpt)) - kpts%bk(:, ikpt)) < 1e-12)) THEN
               isym1 = isym1 + 1
               sym1(isym1, ikpt) = isym
            END IF
         END DO
         nsym1(ikpt) = isym1
      END DO
      ! Define reduced lists of G points -> pgptm1(:,ikpt), ikpt=1,..,nkpt
      !if(allocated(pgptm1)) deallocate(hybrid%pgptm1)
      allocate(pgptm1(maxval(mpbasis%ngptm),kpts%nkptf), source=0) !in mixedbasis
      allocate(iarr(maxval(mpbasis%ngptm)), source=0)
      allocate(POINTER(kpts%nkpt,&
                      MINVAL(mpbasis%gptm(1, :)) - 1:MAXVAL(mpbasis%gptm(1, :)) + 1, &
                      MINVAL(mpbasis%gptm(2, :)) - 1:MAXVAL(mpbasis%gptm(2, :)) + 1, &
                      MINVAL(mpbasis%gptm(3, :)) - 1:MAXVAL(mpbasis%gptm(3, :)) + 1), &
                      source=0)
      allocate(ngptm1, mold=mpbasis%ngptm)
      ngptm1 = 0

      DO ikpt = 1, kpts%nkpt
         DO igpt = 1, mpbasis%ngptm(ikpt)
            g = mpbasis%gptm(:, mpbasis%gptm_ptr(igpt, ikpt))
            POINTER(ikpt, g(1), g(2), g(3)) = igpt
         END DO
         iarr = 0
         j = 0
         DO igpt = mpbasis%ngptm(ikpt), 1, -1
            IF (iarr(igpt) == 0) THEN
               j = j + 1
               pgptm1(j, ikpt) = igpt
               DO isym1 = 1, nsym1(ikpt)
                  g = MATMUL(rrot(:, :, sym1(isym1, ikpt)), mpbasis%gptm(:, mpbasis%gptm_ptr(igpt, ikpt)))
                  i = POINTER(ikpt, g(1), g(2), g(3))
                  IF (i == 0) call judft_error('coulombmatrix: zero pointer (bug?)')
                  iarr(i) = 1
               END DO
            END IF
         END DO
         ngptm1(ikpt) = j
      END DO
      deallocate(iarr)

      IF (mpi%irank == 0) WRITE (6, '(12X,A)', advance='no') 'done'
      CALL cpu_TIME(time2)
      IF (mpi%irank == 0) WRITE (6, '(2X,A,F8.2,A)') '( Timing:', time2 - time1, ' )'

      ! Distribute the work as equally as possible over the processes
      ikptmin = 1
      ikptmax = kpts%nkpt
      igptmin = 1
      igptmax = ngptm1(:kpts%nkpt)
      calc_mt = .TRUE.
      nkminmax = kpts%nkpt

      IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') 'Preparations...'
      CALL cpu_TIME(time1)

      call timestart("define gmat")
      ! Define gmat (symmetric)
      i = (hybrid%lexp + 1)**2
      allocate(gmat(i, i))
      gmat = 0
      lm1 = 0
      DO l1 = 0, hybrid%lexp
         DO m1 = -l1, l1
            lm1 = lm1 + 1
            lm2 = 0
            lp1: DO l2 = 0, l1
               DO m2 = -l2, l2
                  lm2 = lm2 + 1
                  IF (lm2 > lm1) EXIT lp1 ! Don't cross the diagonal!
                  gmat(lm1, lm2) = facB(l1 + l2 + m2 - m1)*facB(l1 + l2 + m1 - m2)/ &
                                   (facB(l1 + m1)*facB(l1 - m1)*facB(l2 + m2)*facB(l2 - m2))/ &
                                   SQRT(1.0*(2*l1 + 1)*(2*l2 + 1)*(2*(l1 + l2) + 1))*(4*pi_const)**1.5
                  gmat(lm2, lm1) = gmat(lm1, lm2)
               END DO
            END DO LP1
         END DO
      END DO
      call timestop("define gmat")

      ! Calculate moments of MT functions
      call timestart("calc moments of MT")
      DO itype = 1, atoms%ntype
         DO l = 0, hybrid%lcutm1(itype)
            DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
               ! note that hybrid%basm1 already contains the factor rgrid
               moment(i, l, itype) = intgrf(atoms%rmsh(:, itype)**(l + 1)*hybrid%basm1(:, i, l, itype), &
                                            atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
            END DO
         END DO
         DO i = 1, mpbasis%num_rad_bas_fun(0, itype)
            moment2(i, itype) = intgrf(atoms%rmsh(:, itype)**3*hybrid%basm1(:, i, 0, itype), &
                                       atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
         END DO
      END DO
      call timestop("calc moments of MT")

      call timestart("getnorm")
      ! Look for different qnorm = |k+G|, definition of qnrm and pqnrm.
      CALL getnorm(kpts, mpbasis%gptm, mpbasis%ngptm, mpbasis%gptm_ptr, qnrm, nqnrm, pqnrm, cell)
      allocate(sphbesmoment(0:hybrid%lexp, atoms%ntype, nqnrm), &
                olap(maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), atoms%ntype, nqnrm), &
                integral(maxval(mpbasis%num_rad_bas_fun), 0:maxval(hybrid%lcutm1), atoms%ntype, nqnrm))
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
      DO iqnrm = iqnrmstart, nqnrm, iqnrmstep
         qnorm = qnrm(iqnrm)
         DO itype = 1, atoms%ntype
            ng = atoms%jri(itype)
            rdum = atoms%rmt(itype)
            sphbes_var = 0
            sphbesmoment1 = 0
            IF (abs(qnorm) < 1e-12) THEN
               sphbesmoment(0, itype, iqnrm) = rdum**3/3
               DO i = 1, ng
                  sphbes_var(i, 0) = 1
                  sphbesmoment1(i, 0) = atoms%rmsh(i, itype)**2/3 + (rdum**2 - atoms%rmsh(i, itype)**2)/2
               END DO
            ELSE
               call sphbes(hybrid%lexp + 1, qnorm*rdum, rarr)
               DO l = 0, hybrid%lexp
                  sphbesmoment(l, itype, iqnrm) = rdum**(l + 2)*rarr(l + 1)/qnorm
               END DO
               DO i = ng, 1, -1
                  rdum = atoms%rmsh(i, itype)
                  call sphbes(hybrid%lcutm1(itype) + 1, qnorm*rdum, rarr)
                  DO l = 0, hybrid%lcutm1(itype)
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
            DO l = 0, hybrid%lcutm1(itype)
               DO n = 1, mpbasis%num_rad_bas_fun(l, itype)
                  ! note that hybrid%basm1 already contains one factor rgrid
                  olap(n, l, itype, iqnrm) = &
                     intgrf(atoms%rmsh(:, itype)*hybrid%basm1(:, n, l, itype)*sphbes_var(:, l), &
                            atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

                  integral(n, l, itype, iqnrm) = &
                     intgrf(atoms%rmsh(:, itype)*hybrid%basm1(:, n, l, itype)*sphbesmoment1(:, l), &
                            atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

               END DO
            END DO
         END DO
      END DO
      call timestop("Bessel calculation")

      IF (mpi%irank == 0) THEN
         WRITE (6, '(18X,A)', advance='no') 'done'
         CALL cpu_TIME(time2)
         WRITE (6, '(2X,A,F8.2,A)', advance='no') '( Timing:', time2 - time1, ' )'
         WRITE (6, *)
      END IF

      !
      !     (1) Case < MT | v | MT >
      !

      IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') '< MT | v | MT > contribution...'

      CALL cpu_TIME(time1)

      IF (ANY(calc_mt)) THEN

         !       (1a) r,r' in same MT
         call timestart("loop 1")
         ix = 0
         iy = 0
         iy0 = 0
         DO itype = 1, atoms%ntype
            DO ineq = 1, atoms%neq(itype)
               ! Here the diagonal block matrices do not depend on ineq. In (1b) they do depend on ineq, though,
               DO l = 0, hybrid%lcutm1(itype)
                  DO n2 = 1, mpbasis%num_rad_bas_fun(l, itype)
                     ! note that hybrid%basm1 already contains the factor rgrid
                     CALL primitivef(primf1, hybrid%basm1(:, n2, l, itype) &
                                     *atoms%rmsh(:, itype)**(l + 1), atoms%rmsh, atoms%dx, &
                                     atoms%jri, atoms%jmtd, itype, atoms%ntype)
                     ! -itype is to enforce inward integration
                     CALL primitivef(primf2, hybrid%basm1(:atoms%jri(itype), n2, l, itype) &
                                     /atoms%rmsh(:atoms%jri(itype), itype)**l, atoms%rmsh, atoms%dx, &
                                     atoms%jri, atoms%jmtd, -itype, atoms%ntype)

                     primf1(:atoms%jri(itype)) = primf1(:atoms%jri(itype))/atoms%rmsh(:atoms%jri(itype), itype)**l
                     primf2 = primf2*atoms%rmsh(:, itype)**(l + 1)

                     DO n1 = 1, n2
                        integrand = hybrid%basm1(:, n1, l, itype)*(primf1 + primf2)
                        !                 call intgr0( (4*pimach())/(2*l+1)*integrand,rmsh(1,itype),dx(itype),jri(itype),mat(n2*(n2-1)/2+n1) )
                        mat(n2*(n2 - 1)/2 + n1) = (4*pi_const)/(2*l + 1) &
                                                  *intgrf(integrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, &
                                                          atoms%ntype, itype, gridf)
                     END DO
                  END DO

                  ! distribute mat for m=-l,l on coulomb in block-matrix form
                  DO M = -l, l
                     DO n2 = 1, mpbasis%num_rad_bas_fun(l, itype)
                        ix = ix + 1
                        iy = iy0
                        DO n1 = 1, n2
                           iy = iy + 1
                           i = ix*(ix - 1)/2 + iy
                           j = n2*(n2 - 1)/2 + n1
                           coulomb(i, kpts%nkpt) = mat(j)
                        END DO
                     END DO
                     iy0 = ix
                  END DO

               END DO
            END DO
         END DO
         call timestop("loop 1")

         !       (1b) r,r' in different MT

         allocate(coulmat(hybrid%nbasp, hybrid%nbasp), stat=ok)
         IF (ok /= 0) call judft_error('coulombmatrix: failure allocation coulmat')
         coulmat = 0

      END IF

      DO ikpt = ikptmin, ikptmax

         ! only the first rank handles the MT-MT part
         call timestart("MT-MT part")
         IF (calc_mt(ikpt)) THEN

            ix = 0
            ic2 = 0
            DO itype2 = 1, atoms%ntype
               DO ineq2 = 1, atoms%neq(itype2)
                  ic2 = ic2 + 1
                  lm2 = 0
                  DO l2 = 0, hybrid%lcutm1(itype2)
                     DO m2 = -l2, l2
                        lm2 = lm2 + 1
                        DO n2 = 1, mpbasis%num_rad_bas_fun(l2, itype2)
                           ix = ix + 1

                           iy = 0
                           ic1 = 0
                           lp2: DO itype1 = 1, itype2
                              DO ineq1 = 1, atoms%neq(itype1)
                                 ic1 = ic1 + 1
                                 lm1 = 0
                                 DO l1 = 0, hybrid%lcutm1(itype1)
                                    DO m1 = -l1, l1
                                       lm1 = lm1 + 1
                                       DO n1 = 1, mpbasis%num_rad_bas_fun(l1, itype1)
                                          iy = iy + 1
                                          IF (iy > ix) EXIT lp2 ! Don't cross the diagonal!
                                          rdum = (-1)**(l2 + m2)*moment(n1, l1, itype1)*moment(n2, l2, itype2)*gmat(lm1, lm2)
                                          l = l1 + l2
                                          lm = l**2 + l + m1 - m2 + 1
                                          idum = ix*(ix - 1)/2 + iy
                                          coulmat(iy, ix) = coulomb(idum, kpts%nkpt) &
                                                            + EXP(CMPLX(0.0, 1.0)*2*pi_const* &
                                                                  dot_PRODUCT(kpts%bk(:, ikpt), &
                                                                              atoms%taual(:, ic2) - atoms%taual(:, ic1))) &
                                                            *rdum*structconst(lm, ic1, ic2, ikpt)
                                          coulmat(ix, iy) = CONJG(coulmat(iy, ix))
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

            IF (sym%invs) THEN
               !symmetrize makes the Coulomb matrix real symmetric
               CALL symmetrize(coulmat, hybrid%nbasp, hybrid%nbasp, 3, .FALSE., &
                               atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                               mpbasis%num_rad_bas_fun, sym)
            ENDIF

            coulomb(:hybrid%nbasp*(hybrid%nbasp + 1)/2, ikpt) = packmat(coulmat)

         END IF
         call timestop("MT-MT part")

      END DO
      IF (ANY(calc_mt)) deallocate(coulmat)

      IF (mpi%irank == 0) THEN
         WRITE (6, '(2X,A)', advance='no') 'done'
         CALL cpu_TIME(time2)
         WRITE (6, '(2X,A,F8.2,A)', advance='no') '( Timing:', time2 - time1, ' )'
         WRITE (6, *)
      END IF

      IF (maxval(mpbasis%ngptm) /= 0) THEN ! skip calculation of plane-wave contribution if mixed basis does not contain plane waves

         !
         !     (2) Case < MT | v | PW >
         !

         IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') '< MT | v | PW > contribution...'

         CALL cpu_TIME(time1)

         !     (2a) r in MT, r' everywhere
         !     (2b) r,r' in same MT
         !     (2c) r,r' in different MT

         allocate(coulmat(hybrid%nbasp, maxval(mpbasis%ngptm)), stat=ok)
         IF (ok /= 0) call judft_error('coulombmatrix: failure allocation coulmat')
         coulmat = 0

         call timestart("loop over interst.")
         DO ikpt = ikptmin, ikptmax !1,kpts%nkpt

            coulmat = 0
            ! start to loop over interstitial plane waves
            DO igpt0 = igptmin(ikpt), igptmax(ikpt) !1,ngptm1(ikpt)
               igpt = pgptm1(igpt0, ikpt)
               igptp = mpbasis%gptm_ptr(igpt, ikpt)
               ix = hybrid%nbasp + igpt
               q = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp), cell%bmat)
               qnorm = norm2(q)
               iqnrm = pqnrm(igpt, ikpt)
               IF (ABS(qnrm(iqnrm) - qnorm) > 1e-12) then
                  call judft_error('coulombmatrix: qnorm does not equal corresponding & element in qnrm (bug?)') ! We shouldn't stop here!
               endif

               call timestart("harmonics")
               call ylm4(2, MATMUL(kpts%bk(:, kpts%nkpt), cell%bmat), y1)
               call ylm4(2, MATMUL(mpbasis%gptm(:, igptp), cell%bmat), y2)
               call ylm4(hybrid%lexp, q, y)
               call timestop("harmonics")
               y1 = CONJG(y1); y2 = CONJG(y2); y = CONJG(y)

               iy = 0
               ic = 0
               DO itype = 1, atoms%ntype
                  DO ineq = 1, atoms%neq(itype)
                     ic = ic + 1
                     lm = 0
                     DO l = 0, hybrid%lcutm1(itype)
                        DO M = -l, l
                           lm = lm + 1

                           ! calculate sum over lm and centers for (2c) -> csum, csumf
                           csum = 0
                           csumf = 0
                           ic1 = 0
                           DO itype1 = 1, atoms%ntype
                              DO ineq1 = 1, atoms%neq(itype1)
                                 ic1 = ic1 + 1
                                 cexp = 4*pi_const*EXP(CMPLX(0.0, 1.0)*2*pi_const &
                                                       *(dot_PRODUCT(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp), atoms%taual(:, ic1)) &
                                                         - dot_PRODUCT(kpts%bk(:, ikpt), atoms%taual(:, ic))))

                                 lm1 = 0
                                 DO l1 = 0, hybrid%lexp
                                    l2 = l + l1 ! for structconst
                                    idum = 1
                                    cdum = sphbesmoment(l1, itype1, iqnrm)*CMPLX(0.0, 1.0)**(l1)*cexp
                                    DO m1 = -l1, l1
                                       lm1 = lm1 + 1
                                       m2 = M - m1              ! for structconst
                                       lm2 = l2**2 + l2 + m2 + 1 !
                                       csum = csum - idum*gmat(lm1, lm)*y(lm1)*cdum*structconst(lm2, ic, ic1, ikpt)
                                       idum = -idum ! factor (-1)*(l1+m1)
                                    END DO
                                 END DO

                                 ! add contribution of (2c) to csum and csumf coming from linear and quadratic orders of Y_lm*(G) / G * j_(l+1)(GS)
                                 IF (ikpt == 1 .AND. l <= 2) THEN
                                    cexp = EXP(CMPLX(0.0, 1.0)*2*pi_const*dot_PRODUCT(mpbasis%gptm(:, igptp), atoms%taual(:, ic1))) &
                                           *gmat(lm, 1)*4*pi_const/cell%vol
                                    csumf(lm) = csumf(lm) - cexp*SQRT(4*pi_const)* &
                                                CMPLX(0.0, 1.0)**l*sphbesmoment(0, itype1, iqnrm)/facC(l - 1)
                                    IF (l == 0) THEN
                                       IF (igpt /= 1) THEN
                                          csum = csum - cexp*(sphbesmoment(0, itype1, iqnrm)*atoms%rmt(itype1)**2 - &
                                                              sphbesmoment(2, itype1, iqnrm)*2.0/3)/10
                                       ELSE
                                          csum = csum - cexp*atoms%rmt(itype1)**5/30
                                       END IF
                                    ELSE IF (l == 1) THEN
                                       csum = csum + cexp*CMPLX(0.0, 1.0)*SQRT(4*pi_const) &
                                              *sphbesmoment(1, itype1, iqnrm)*y(lm)/3
                                    END IF
                                 END IF

                              END DO
                           END DO

                           ! add contribution of (2a) to csumf
                           IF (ikpt == 1 .AND. igpt == 1 .AND. l <= 2) THEN
                              csumf(lm) = csumf(lm) + (4*pi_const)**2*CMPLX(0.0, 1.0)**l/facC(l)
                           END IF

                           ! finally define coulomb
                           idum = ix*(ix - 1)/2
                           cdum = (4*pi_const)**2*CMPLX(0.0, 1.0)**(l)*y(lm) &
                                  *EXP(CMPLX(0.0, 1.0)*2*pi_const &
                                       *dot_PRODUCT(mpbasis%gptm(:, igptp), atoms%taual(:, ic)))
                           DO n = 1, mpbasis%num_rad_bas_fun(l, itype)
                              iy = iy + 1

                              IF (ikpt == 1 .AND. igpt == 1) THEN
                                 IF (l == 0) coulmat(iy, ix - hybrid%nbasp) = &
                                    -cdum*moment2(n, itype)/6/svol         ! (2a)
                                 coulmat(iy, ix - hybrid%nbasp) = coulmat(iy, ix - hybrid%nbasp) &
                                                                  + (-cdum/(2*l + 1)*integral(n, l, itype, iqnrm) & ! (2b)&
                                                                     + csum*moment(n, l, itype))/svol          ! (2c)
                              ELSE
                                 coulmat(iy, ix - hybrid%nbasp) = &
                                    (cdum*olap(n, l, itype, iqnrm)/qnorm**2 &  ! (2a)&
                                     - cdum/(2*l + 1)*integral(n, l, itype, iqnrm) & ! (2b)&
                                     + csum*moment(n, l, itype))/svol          ! (2c)

                              END IF

                           END DO

                        END DO
                     END DO
                  END DO

               END DO
            END DO

            IF (sym%invs) THEN
               CALL symmetrize(coulmat, hybrid%nbasp, mpbasis%ngptm(ikpt), 1, .FALSE., &
                               atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), mpbasis%num_rad_bas_fun, sym)
            ENDIF

            M = hybrid%nbasp*(hybrid%nbasp + 1)/2
            DO i = 1, mpbasis%ngptm(ikpt)
               DO j = 1, hybrid%nbasp + i
                  M = M + 1
                  IF (j <= hybrid%nbasp) coulomb(M, ikpt) = coulmat(j, i)
               END DO
            END DO
         END DO
         call timestop("loop over interst.")

         deallocate(coulmat, olap, integral)

         IF (mpi%irank == 0) THEN
            WRITE (6, '(2X,A)', advance='no') 'done'
            CALL cpu_TIME(time2)
            WRITE (6, '(2X,A,F8.2,A)') '( Timing:', time2 - time1, ' )'
         END IF

         !
         !     (3) Case < PW | v | PW >
         !

         IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') '< PW | v | PW > contribution...'

         CALL cpu_TIME(time1)

         !     (3a) r,r' everywhere; r everywhere, r' in MT; r in MT, r' everywhere

         CALL cpu_TIME(time1)
         ! Calculate the hermitian matrix smat(i,j) = sum(a) integral(MT(a)) exp[i(Gj-Gi)r] dr
         call timestart("calc smat")
         allocate(smat(mpbasis%num_gpts(), mpbasis%num_gpts()))
         smat = 0
         DO igpt2 = 1, mpbasis%num_gpts()
            DO igpt1 = 1, igpt2
               g = mpbasis%gptm(:, igpt2) - mpbasis%gptm(:, igpt1)
               gnorm = gptnorm(g, cell%bmat)
               IF (abs(gnorm) < 1e-12) THEN
                  DO itype = 1, atoms%ntype
                     smat(igpt1, igpt2) = smat(igpt1, igpt2) + atoms%neq(itype)*4*pi_const*atoms%rmt(itype)**3/3
                  END DO
               ELSE
                  ic = 0
                  DO itype = 1, atoms%ntype
                     rdum = atoms%rmt(itype)*gnorm
                     rdum = 4*pi_const*(SIN(rdum) - rdum*COS(rdum))/gnorm**3
                     DO ineq = 1, atoms%neq(itype)
                        ic = ic + 1
                        smat(igpt1, igpt2) = smat(igpt1, igpt2) &
                                             + rdum*EXP(CMPLX(0.0, 1.0)*2*pi_const*dot_PRODUCT(atoms%taual(:, ic), g))
                     END DO
                  END DO
               END IF
               smat(igpt2, igpt1) = CONJG(smat(igpt1, igpt2))
            END DO
         END DO
         call timestop("calc smat")

         ! Coulomb matrix, contribution (3a)
         call timestart("coulomb matrix")
         DO ikpt = ikptmin, ikptmax

            DO igpt0 = igptmin(ikpt), igptmax(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               igptp2 = mpbasis%gptm_ptr(igpt2, ikpt)
               ix = hybrid%nbasp + igpt2
               iy = hybrid%nbasp
               q2 = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp2), cell%bmat)
               rdum2 = SUM(q2**2)
               IF (abs(rdum2) > 1e-12) rdum2 = 4*pi_const/rdum2

               DO igpt1 = 1, igpt2
                  igptp1 = mpbasis%gptm_ptr(igpt1, ikpt)
                  iy = iy + 1
                  q1 = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp1), cell%bmat)
                  idum = ix*(ix - 1)/2 + iy
                  rdum1 = SUM(q1**2)
                  IF (abs(rdum1) > 1e-12) rdum1 = 4*pi_const/rdum1

                  IF (ikpt == 1) THEN
                     IF (igpt1 /= 1) THEN
                        coulomb(idum, 1) = -smat(igptp1, igptp2)*rdum1/cell%vol
                     END IF
                     IF (igpt2 /= 1) THEN
                        coulomb(idum, 1) = coulomb(idum, 1) - smat(igptp1, igptp2)*rdum2/cell%vol
                     END IF
                  ELSE
                     coulomb(idum, ikpt) = -smat(igptp1, igptp2)*(rdum1 + rdum2)/cell%vol
                  END IF
               END DO
               IF (ikpt /= 1 .OR. igpt2 /= 1) THEN                  !
                  coulomb(idum, ikpt) = coulomb(idum, ikpt) + rdum2 ! diagonal term
               END IF                                            !
            END DO

         END DO
         call timestop("coulomb matrix")
         !     (3b) r,r' in different MT

         call timestart("loop 4:")
         DO ikpt = ikptmin, ikptmax!1,kpts%nkpt

            ! group together quantities which depend only on l,m and igpt -> carr2a
            allocate(carr2a((hybrid%lexp + 1)**2, maxval(mpbasis%ngptm)), carr2b(atoms%nat, maxval(mpbasis%ngptm)))
            carr2a = 0; carr2b = 0
            DO igpt = 1, mpbasis%ngptm(ikpt)
               igptp = mpbasis%gptm_ptr(igpt, ikpt)
               iqnrm = pqnrm(igpt, ikpt)
               q = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp), cell%bmat)

               call timestart("harmonics")
               call ylm4(hybrid%lexp, q, y)
               call timestop("harmonics")

               y = CONJG(y)
               lm = 0
               DO l = 0, hybrid%lexp
                  DO M = -l, l
                     lm = lm + 1
                     carr2a(lm, igpt) = 4*pi_const*CMPLX(0.0, 1.0)**(l)*y(lm)
                  END DO
               END DO
               DO ic = 1, atoms%nat
                  carr2b(ic, igpt) = EXP(-CMPLX(0.0, 1.0)*2*pi_const* &
                                         dot_PRODUCT(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp), atoms%taual(:, ic)))
               END DO
            END DO

            !finally we can loop over the plane waves (G: igpt1,igpt2)
            call timestart("loop over plane waves")
            allocate(carr2(atoms%nat, (hybrid%lexp + 1)**2), &
                      structconst1(atoms%nat, (2*hybrid%lexp + 1)**2))
            carr2 = 0; structconst1 = 0
            DO igpt0 = igptmin(ikpt), igptmax(ikpt)!1,ngptm1(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               ix = hybrid%nbasp + igpt2
               igptp2 = mpbasis%gptm_ptr(igpt2, ikpt)
               iqnrm2 = pqnrm(igpt2, ikpt)
               ic2 = 0
               carr2 = 0
               DO itype2 = 1, atoms%ntype
                  DO ineq2 = 1, atoms%neq(itype2)
                     ic2 = ic2 + 1
                     cexp = CONJG(carr2b(ic2, igpt2))
                     lm2 = 0
                     DO ic1 = 1, atoms%nat
                        structconst1(ic1, :) = structconst(:, ic1, ic2, ikpt)
                     END DO
                     DO l2 = 0, hybrid%lexp
                        idum = 1
                        DO m2 = -l2, l2
                           lm2 = lm2 + 1
                           cdum = idum*sphbesmoment(l2, itype2, iqnrm2)*cexp*carr2a(lm2, igpt2)
                           IF (abs(cdum) > 1e-12) THEN
                              lm1 = 0
                              DO l1 = 0, hybrid%lexp
                                 l = l1 + l2
                                 M = -l1 - m2 !first loop of m1
                                 lm = l**2 + l + M
                                 DO m1 = -l1, l1
                                    lm1 = lm1 + 1
                                    lm = lm + 1
                                    cdum1 = cdum*gmat(lm1, lm2)
                                    DO ic1 = 1, atoms%nat
                                       carr2(ic1, lm1) = carr2(ic1, lm1) + cdum1*structconst1(ic1, lm)
                                    END DO
                                 END DO
                              END DO
                           END IF
                           idum = -idum !factor (-1)**(l+M)
                        END DO
                     END DO
                  END DO
               END DO

               iy = hybrid%nbasp
               DO igpt1 = 1, igpt2
                  iy = iy + 1
                  igptp1 = mpbasis%gptm_ptr(igpt1, ikpt)
                  iqnrm1 = pqnrm(igpt1, ikpt)
                  csum = 0
                  ic = 0
                  DO itype = 1, atoms%ntype
                     DO ineq = 1, atoms%neq(itype)
                        ic = ic + 1
                        cexp = carr2b(ic, igpt1)
                        lm = 0
                        DO l = 0, hybrid%lexp
                           cdum = cexp*sphbesmoment(l, itype, iqnrm1)
                           DO M = -l, l
                              lm = lm + 1
                              csum = csum + cdum*carr2(ic, lm)*CONJG(carr2a(lm, igpt1)) ! for coulomb
                           END DO
                        END DO
                     END DO
                  END DO
                  idum = ix*(ix - 1)/2 + iy
                  coulomb(idum, ikpt) = coulomb(idum, ikpt) + csum/cell%vol
               END DO
            END DO
            deallocate(carr2, carr2a, carr2b, structconst1)
            call timestop("loop over plane waves")
         END DO !ikpt
         call timestop("loop 4:")
         !     Add corrections from higher orders in (3b) to coulomb(:,1)
         ! (1) igpt1 > 1 , igpt2 > 1  (finite G vectors)
         call timestart("add corrections from higher orders")
         rdum = (4*pi_const)**(1.5)/cell%vol**2*gmat(1, 1)
         DO igpt0 = 1, ngptm1(1)
            igpt2 = pgptm1(igpt0, 1); IF (igpt2 == 1) CYCLE
            ix = hybrid%nbasp + igpt2
            iqnrm2 = pqnrm(igpt2, 1)
            igptp2 = mpbasis%gptm_ptr(igpt2, 1)
            q2 = MATMUL(mpbasis%gptm(:, igptp2), cell%bmat)
            qnorm2 = norm2(q2)
            iy = hybrid%nbasp + 1
            DO igpt1 = 2, igpt2
               iy = iy + 1
               idum = ix*(ix - 1)/2 + iy
               iqnrm1 = pqnrm(igpt1, 1)
               igptp1 = mpbasis%gptm_ptr(igpt1, 1)
               q1 = MATMUL(mpbasis%gptm(:, igptp1), cell%bmat)
               qnorm1 = norm2(q1)
               rdum1 = dot_PRODUCT(q1, q2)/(qnorm1*qnorm2)
               ic1 = 0
               DO itype1 = 1, atoms%ntype
                  DO ineq1 = 1, atoms%neq(itype1)
                     ic1 = ic1 + 1
                     ic2 = 0
                     DO itype2 = 1, atoms%ntype
                        DO ineq2 = 1, atoms%neq(itype2)
                           ic2 = ic2 + 1
                           cdum = EXP(CMPLX(0.0, 1.0)*2*pi_const* &
                                      (-dot_PRODUCT(mpbasis%gptm(:, igptp1), atoms%taual(:, ic1)) &
                                       + dot_PRODUCT(mpbasis%gptm(:, igptp2), atoms%taual(:, ic2))))
                           coulomb(idum, 1) = coulomb(idum, 1) + rdum*cdum*( &
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
            END DO
         END DO
         ! (2) igpt1 = 1 , igpt2 > 1  (first G vector vanishes, second finite)
         iy = hybrid%nbasp + 1
         DO igpt0 = 1, ngptm1(1)
            igpt2 = pgptm1(igpt0, 1); IF (igpt2 == 1) CYCLE
            ix = hybrid%nbasp + igpt2
            iqnrm2 = pqnrm(igpt2, 1)
            igptp2 = mpbasis%gptm_ptr(igpt2, 1)
            qnorm2 = qnrm(iqnrm2)
            idum = ix*(ix - 1)/2 + iy
            DO itype1 = 1, atoms%ntype
               DO ineq1 = 1, atoms%neq(itype1)
                  ic2 = 0
                  DO itype2 = 1, atoms%ntype
                     DO ineq2 = 1, atoms%neq(itype2)
                        ic2 = ic2 + 1
                        cdum = EXP(CMPLX(0.0, 1.0)*2*pi_const*dot_PRODUCT(mpbasis%gptm(:, igptp2), atoms%taual(:, ic2)))
                        coulomb(idum, 1) = coulomb(idum, 1) &
                                           + rdum*cdum*atoms%rmt(itype1)**3*( &
                                           +sphbesmoment(0, itype2, iqnrm2)/30*atoms%rmt(itype1)**2 &
                                           - sphbesmoment(2, itype2, iqnrm2)/18 &
                                           + sphbesmoment(1, itype2, iqnrm2)/6/qnorm2)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         ! (2) igpt1 = 1 , igpt2 = 1  (vanishing G vectors)
         iy = hybrid%nbasp + 1
         ix = hybrid%nbasp + 1
         idum = ix*(ix - 1)/2 + iy
         DO itype1 = 1, atoms%ntype
            DO ineq1 = 1, atoms%neq(itype1)
               DO itype2 = 1, atoms%ntype
                  DO ineq2 = 1, atoms%neq(itype2)
                     coulomb(idum, 1) = coulomb(idum, 1) &
                                        + rdum*atoms%rmt(itype1)**3*atoms%rmt(itype2)**3* &
                                        (atoms%rmt(itype1)**2 + atoms%rmt(itype2)**2)/90
                  END DO
               END DO
            END DO
         END DO
         call timestop("add corrections from higher orders")

         !     (3c) r,r' in same MT

         ! Calculate sphbesintegral
         call timestart("sphbesintegral")
         allocate(sphbes0(-1:hybrid%lexp + 2, atoms%ntype, nqnrm),&
              &           carr2((hybrid%lexp + 1)**2, maxval(mpbasis%ngptm)))
         sphbes0 = 0; carr2 = 0
         DO iqnrm = 1, nqnrm
            DO itype = 1, atoms%ntype
               rdum = qnrm(iqnrm)*atoms%rmt(itype)
               call sphbes(hybrid%lexp + 2, rdum, sphbes0(0, itype, iqnrm))
               IF (abs(rdum) > 1e-12) sphbes0(-1, itype, iqnrm) = COS(rdum)/rdum
            END DO
         END DO
         call timestop("sphbesintegral")

         l_warn = (mpi%irank == 0)
         call timestart("loop 2")
         DO ikpt = ikptmin, ikptmax!1,nkpt

            DO igpt = 1, mpbasis%ngptm(ikpt)
               igptp = mpbasis%gptm_ptr(igpt, ikpt)
               q = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp), cell%bmat)
               call timestart("harmonics")
               call ylm4(hybrid%lexp, q, carr2(:, igpt))
               call timestop("harmonics")
            END DO

            DO igpt0 = igptmin(ikpt), igptmax(ikpt)!1,ngptm1(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               ix = hybrid%nbasp + igpt2
               igptp2 = mpbasis%gptm_ptr(igpt2, ikpt)
               iqnrm2 = pqnrm(igpt2, ikpt)
               q2 = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp2), cell%bmat)
               y2 = CONJG(carr2(:, igpt2))
               iy = hybrid%nbasp
               DO igpt1 = 1, igpt2
                  iy = iy + 1
                  igptp1 = mpbasis%gptm_ptr(igpt1, ikpt)
                  iqnrm1 = pqnrm(igpt1, ikpt)
                  q1 = MATMUL(kpts%bk(:, ikpt) + mpbasis%gptm(:, igptp1), cell%bmat)
                  y1 = carr2(:, igpt1)
                  cexp1 = 0
                  ic = 0
                  DO itype = 1, atoms%ntype
                     DO ineq = 1, atoms%neq(itype)
                        ic = ic + 1
                        cexp1(itype) = cexp1(itype) + &
                                       EXP(CMPLX(0.0, 1.0)*2*pi_const*dot_PRODUCT( &
                                           (mpbasis%gptm(:, igptp2) - mpbasis%gptm(:, igptp1)), atoms%taual(:, ic)))
                     ENDDO
                  ENDDO
                  lm = 0
                  cdum = 0
                  DO l = 0, hybrid%lexp
                     cdum1 = 0
                     DO itype = 1, atoms%ntype
                        cdum1 = cdum1 + cexp1(itype)*sphbessel_integral( &
                                atoms, itype, qnrm, nqnrm, &
                                iqnrm1, iqnrm2, l, hybrid, &
                                sphbes0, l_warn, l_warned) &
                                /(2*l + 1)
                        l_warn = l_warn .AND. .NOT. l_warned ! only warn once
                     END DO
                     DO M = -l, l
                        lm = lm + 1
                        cdum = cdum + cdum1*y1(lm)*y2(lm)
                     ENDDO
                  ENDDO
                  idum = ix*(ix - 1)/2 + iy
                  coulomb(idum, ikpt) = coulomb(idum, ikpt) + (4*pi_const)**3*cdum/cell%vol
               END DO
            END DO

         END DO
         call timestop("loop 2")
         deallocate(carr2)

         IF (mpi%irank == 0) THEN
            WRITE (6, '(2X,A)', advance='no') 'done'
            CALL cpu_TIME(time2)
            WRITE (6, '(2X,A,F8.2,A)') '( Timing:', time2 - time1, ' )'
         END IF

         !
         !     Symmetry-equivalent G vectors
         !

         IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') 'Symm.-equiv. matrix elements...'
         CALL cpu_TIME(time1)
         ! All elements are needed so send all data to all processes treating the
         ! respective k-points

         allocate(carr2(hybrid%maxbasm1, 2), iarr(maxval(mpbasis%ngptm)))
         allocate(nsym_gpt(mpbasis%num_gpts(), kpts%nkpt), &
                   sym_gpt(MAXVAL(nsym1), mpbasis%num_gpts(), kpts%nkpt))
         nsym_gpt = 0; sym_gpt = 0
         call timestart("loop 3")
         DO ikpt = ikptmin, ikptmax
            carr2 = 0; iarr = 0
            iarr(pgptm1(:ngptm1(ikpt), ikpt)) = 1
            DO igpt0 = 1, ngptm1(ikpt) !igptmin(ikpt),igptmax(ikpt)
               lsym = ((igptmin(ikpt) <= igpt0) .AND. &
                       (igptmax(ikpt) >= igpt0))
               igpt2 = pgptm1(igpt0, ikpt)
               j = (hybrid%nbasp + igpt2 - 1)*(hybrid%nbasp + igpt2)/2
               i = hybrid%nbasp + igpt2
               carr2(1:i, 2) = coulomb(j + 1:j + i, ikpt)
               j = j + i
               DO i = hybrid%nbasp + igpt2 + 1, nbasm1(ikpt)
                  j = j + i - 1
                  IF (sym%invs) THEN
                     carr2(i, 2) = coulomb(j, ikpt)
                  ELSE
                     carr2(i, 2) = CONJG(coulomb(j, ikpt))
                  ENDIF
               END DO
               IF (lsym) THEN
                  ic = 1
                  sym_gpt(ic, igpt0, ikpt) = igpt2
               END IF
               DO isym1 = 2, nsym1(ikpt)
                  isym = sym1(isym1, ikpt)
                  CALL bramat_trafo( &
                     carr2(:, 1), igpt1, &
                     carr2(:, 2), igpt2, ikpt, isym, .FALSE., POINTER(ikpt, :, :, :), &
                     sym, rrot(:, :, isym), invrrot(:, :, isym), mpbasis, hybrid, &
                     kpts, maxval(hybrid%lcutm1), atoms, hybrid%lcutm1, &
                     mpbasis%num_rad_bas_fun, maxval(mpbasis%num_rad_bas_fun), dwgn(:, :, :, isym), &
                     hybrid%nbasp, nbasm1)
                  IF (iarr(igpt1) == 0) THEN
                     CALL bramat_trafo( &
                        carr2(:, 1), igpt1, &
                        carr2(:, 2), igpt2, ikpt, isym, .TRUE., POINTER(ikpt, :, :, :), &
                        sym, rrot(:, :, isym), invrrot(:, :, isym), mpbasis, hybrid, &
                        kpts, maxval(hybrid%lcutm1), atoms, hybrid%lcutm1, &
                        mpbasis%num_rad_bas_fun, maxval(mpbasis%num_rad_bas_fun), &
                        dwgn(:, :, :, isym), hybrid%nbasp, nbasm1)
                     l = (hybrid%nbasp + igpt1 - 1)*(hybrid%nbasp + igpt1)/2
                     coulomb(l + 1:l + hybrid%nbasp + igpt1, ikpt) = carr2(:hybrid%nbasp + igpt1, 1)
                     iarr(igpt1) = 1
                     IF (lsym) THEN
                        ic = ic + 1
                        sym_gpt(ic, igpt0, ikpt) = igpt1
                     END IF
                  END IF
               END DO
               nsym_gpt(igpt0, ikpt) = ic
            END DO ! igpt0
         END DO ! ikpt
         call timestop("loop 3")
         call timestart("gap 1:")
         deallocate(carr2, iarr, pgptm1)
         IF (mpi%irank == 0) THEN
            WRITE (6, '(2X,A)', advance='no') 'done'
            CALL cpu_TIME(time2)
            WRITE (6, '(2X,A,F8.2,A)') '( Timing:', time2 - time1, ' )'
         END IF
      END IF
      deallocate(qnrm, pqnrm)

      CALL cpu_TIME(time1)
      IF (xcpot%is_name("hse") .OR. xcpot%is_name("vhse")) THEN
         !
         ! The HSE functional is realized subtracting erf/r from
         ! the normal Coulomb matrix
         !
      ELSE
         IF (ikptmin == 1) CALL subtract_sphaverage(sym, cell, atoms, mpbasis, hybrid, nbasm1, gridf, coulomb)
      END IF

      ! transform Coulomb matrix to the biorthogonal set
      IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') 'Transform to biorthogonal set...'
      CALL cpu_TIME(time1)
      call timestop("gap 1:")
      call timestart("calc eigenvalues olap_pw")
      DO ikpt = ikptmin, ikptmax

         !calculate IR overlap-matrix
         CALL olapm%alloc(sym%invs, mpbasis%ngptm(ikpt), mpbasis%ngptm(ikpt), 0.0)

         CALL olap_pw(olapm, mpbasis%gptm(:, mpbasis%gptm_ptr(:mpbasis%ngptm(ikpt), ikpt)), mpbasis%ngptm(ikpt), atoms, cell)

         !         !calculate eigenvalues of olapm
         !         ALLOCATE( eval(ngptm(ikpt)),evec(ngptm(ikpt),ngptm(ikpt)) )
         !         CALL diagonalize(evec,eval,olapm)
         !
         !         !
         !         ! small eigenvalues lead to inaccuries in the inversion
         !         ! however it seems that these do not play an important role
         !         ! for the HF exchange
         !         ! thus we do not do a SingularValueDecomposition
         !         !
         !
         !         IF( any(eval .le. 1E-06 ) ) THEN
         ! !           WRITE(*,*) count( eval .le. 1E-06 )
         ! !           WRITE(*,*) 'eval .le. 1E-06'
         ! !           ALLOCATE( involapm(ngptm(ikpt),ngptm(ikpt)) )
         !           olapm = 0
         !           DO i = 1,ngptm(ikpt)
         !             IF( eval(i) .le. 1E-06) CYCLE
         !             olapm(i,i) = 1/eval(i)
         !           END DO
         !
         !           ALLOCATE( invevec(ngptm(ikpt),ngptm(ikpt)) )
         ! if (sym%invs) then
         !           invevec =         transpose(evec)
         ! else
         !           invevec = conjg( transpose(evec) )
         ! endif
         !
         !           olapm   = matmul(evec,matmul(olapm,invevec) )
         !
         !           DEALLOCATE(invevec)!,involapm)
         !         ELSE
         !calculate inverse overlap-matrix
         CALL olapm%inverse()
         !         END IF

         !unpack matrix coulomb
         CALL coulhlp%from_packed(sym%invs, nbasm1(ikpt), REAL(coulomb(:, ikpt)), coulomb(:, ikpt))

         if (olapm%l_real) THEN
            !multiply with inverse olap from right hand side
            coulhlp%data_r(:, hybrid%nbasp + 1:) = MATMUL(coulhlp%data_r(:, hybrid%nbasp + 1:), olapm%data_r)
            !multiply with inverse olap from left side
            coulhlp%data_r(hybrid%nbasp + 1:, :) = MATMUL(olapm%data_r, coulhlp%data_r(hybrid%nbasp + 1:, :))
         else
            !multiply with inverse olap from right hand side
            coulhlp%data_c(:, hybrid%nbasp + 1:) = MATMUL(coulhlp%data_c(:, hybrid%nbasp + 1:), olapm%data_c)
            !multiply with inverse olap from left side
            coulhlp%data_c(hybrid%nbasp + 1:, :) = MATMUL(olapm%data_c, coulhlp%data_c(hybrid%nbasp + 1:, :))
         end if
         coulomb(:(nbasm1(ikpt)*(nbasm1(ikpt) + 1))/2, ikpt) = coulhlp%to_packed()

      END DO
      call timestop("calc eigenvalues olap_pw")

      IF (mpi%irank == 0) THEN
         WRITE (6, '(1X,A)', advance='no') 'done'
         CALL cpu_TIME(time2)
         WRITE (6, '(2X,A,F8.2,A)') '( Timing:', time2 - time1, ' )'
      END IF

      !call plot_coulombmatrix() -> code was shifted to plot_coulombmatrix.F90

      IF (mpi%irank == 0) WRITE (6, '(A)', advance='no') 'Writing of data to file...'
      CALL cpu_TIME(time1)
      !
      ! rearrange coulomb matrix
      !

      allocate(coulomb_mt1(maxval(mpbasis%num_rad_bas_fun) - 1, maxval(mpbasis%num_rad_bas_fun) - 1, 0:maxval(hybrid%lcutm1), atoms%ntype, 1))
      ic = (maxval(hybrid%lcutm1) + 1)**2*atoms%nat
      idum = ic + maxval(mpbasis%ngptm)
      idum = (idum*(idum + 1))/2
      if (sym%invs) THEN
         allocate(coulomb_mt2_r(maxval(mpbasis%num_rad_bas_fun) - 1, -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1) + 1, atoms%nat, 1))
         allocate(coulomb_mt3_r(maxval(mpbasis%num_rad_bas_fun) - 1, atoms%nat, atoms%nat, 1))
         allocate(coulomb_mtir_r(ic + maxval(mpbasis%ngptm), ic + maxval(mpbasis%ngptm), 1))
         allocate(coulombp_mtir_r(idum, 1))
      else
         allocate(coulomb_mt2_c(maxval(mpbasis%num_rad_bas_fun) - 1, -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1), 0:maxval(hybrid%lcutm1) + 1, atoms%nat, 1))
         allocate(coulomb_mt3_c(maxval(mpbasis%num_rad_bas_fun) - 1, atoms%nat, atoms%nat, 1))
         allocate(coulomb_mtir_c(ic + maxval(mpbasis%ngptm), ic + maxval(mpbasis%ngptm), 1))
         allocate(coulombp_mtir_c(idum, 1))
      endif
      call timestart("loop bla")
      DO ikpt = ikptmin, ikptmax
         ikpt0 = 1
         ikpt1 = 1
         ! initialize arrays
         if (sym%invs) THEN
            coulomb_mt1 = 0; coulomb_mt2_r = 0
            coulomb_mt3_r = 0; coulombp_mtir_r = 0
         else
            coulomb_mt1 = 0; coulomb_mt2_c = 0
            coulomb_mt3_c = 0; coulombp_mtir_c = 0
         endif
         ! unpack coulomb into coulhlp

         call coulhlp%from_packed(sym%invs, nbasm1(ikpt), real(coulomb(:, ikpt)), coulomb(:, ikpt))

         ! only one processor per k-point calculates MT convolution
         IF (calc_mt(ikpt)) THEN

            !
            ! store m-independent part of Coulomb matrix in MT spheres
            ! in coulomb_mt1(:mpbasis%num_rad_bas_fun(l,itype)-1,:mpbasis%num_rad_bas_fun(l,itype)-1,l,itype)
            !
            call timestart("m-indep. part of coulomb mtx")
            indx1 = 0
            DO itype = 1, atoms%ntype
               DO ineq = 1, atoms%neq(itype)
                  DO l = 0, hybrid%lcutm1(itype)

                     IF (ineq == 1) THEN
                        DO n = 1, mpbasis%num_rad_bas_fun(l, itype) - 1
                           if (coulhlp%l_real) THEN
                              coulomb_mt1(n, 1:mpbasis%num_rad_bas_fun(l, itype) - 1, l, itype, ikpt0) &
                                 = coulhlp%data_r(indx1 + n, indx1 + 1:indx1 + mpbasis%num_rad_bas_fun(l, itype) - 1)
                           else
                              coulomb_mt1(n, 1:mpbasis%num_rad_bas_fun(l, itype) - 1, l, itype, ikpt0) &
                                 = coulhlp%data_c(indx1 + n, indx1 + 1:indx1 + mpbasis%num_rad_bas_fun(l, itype) - 1)
                           end if
                        END DO
                     END IF

                     indx1 = indx1 + (2*l + 1)*mpbasis%num_rad_bas_fun(l, itype)
                  END DO
               END DO
            END DO
            call timestop("m-indep. part of coulomb mtx")

            !
            ! store m-dependent and atom-dependent part of Coulomb matrix in MT spheres
            ! in coulomb_mt2(:mpbasis%num_rad_bas_fun(l,itype)-1,-l:l,l,iatom)
            !
            call timestart("m-dep. part of coulomb mtx")
            indx1 = 0
            iatom = 0
            DO itype = 1, atoms%ntype
               DO ineq = 1, atoms%neq(itype)
                  iatom = iatom + 1
                  DO l = 0, hybrid%lcutm1(itype)
                     DO M = -l, l
                        if (coulhlp%l_real) THEN
                           coulomb_mt2_r(:mpbasis%num_rad_bas_fun(l, itype) - 1, M, l, iatom, ikpt0) &
                              = coulhlp%data_r(indx1 + 1:indx1 + mpbasis%num_rad_bas_fun(l, itype) - 1, indx1 + mpbasis%num_rad_bas_fun(l, itype))
                        else
                           coulomb_mt2_c(:mpbasis%num_rad_bas_fun(l, itype) - 1, M, l, iatom, ikpt0) &
                              = coulhlp%data_c(indx1 + 1:indx1 + mpbasis%num_rad_bas_fun(l, itype) - 1, indx1 + mpbasis%num_rad_bas_fun(l, itype))
                        endif

                        indx1 = indx1 + mpbasis%num_rad_bas_fun(l, itype)

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
               ! coulomb_mt2(:mpbasis%num_rad_bas_fun(l=0,itype),0,maxval(hybrid%lcutm1)+1,iatom)
               !
               ic = 0
               iatom = 0
               DO itype = 1, atoms%ntype
                  DO ineq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     DO n = 1, mpbasis%num_rad_bas_fun(0, itype) - 1
                        if (coulhlp%l_real) THEN
                           coulomb_mt2_r(n, 0, maxval(hybrid%lcutm1) + 1, iatom, ikpt0) = coulhlp%data_r(ic + n, hybrid%nbasp + 1)
                        else
                           coulomb_mt2_c(n, 0, maxval(hybrid%lcutm1) + 1, iatom, ikpt0) = coulhlp%data_c(ic + n, hybrid%nbasp + 1)
                        endif
                     END DO
                     ic = ic + SUM([((2*l + 1)*mpbasis%num_rad_bas_fun(l, itype), l=0, hybrid%lcutm1(itype))])
                  END DO
               END DO

               !
               ! store the contributions between the MT s-like functions at atom1 and
               ! and the constant function at a different atom2
               !
               iatom = 0
               ic = 0
               DO itype = 1, atoms%ntype
                  ishift = SUM([((2*l + 1)*mpbasis%num_rad_bas_fun(l, itype), l=0, hybrid%lcutm1(itype))])
                  DO ineq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     ic1 = ic + mpbasis%num_rad_bas_fun(0, itype)

                     iatom1 = 0
                     ic2 = 0
                     DO itype1 = 1, atoms%ntype
                        ishift1 = SUM([((2*l1 + 1)*mpbasis%num_rad_bas_fun(l1, itype1), l1=0, hybrid%lcutm1(itype1))])
                        DO ineq1 = 1, atoms%neq(itype1)
                           iatom1 = iatom1 + 1
                           ic3 = ic2 + 1
                           ic4 = ic3 + mpbasis%num_rad_bas_fun(0, itype1) - 2

                           IF (sym%invs) THEN
                              coulomb_mt3_r(:mpbasis%num_rad_bas_fun(0, itype1) - 1, iatom, iatom1, ikpt0) = coulhlp%data_r(ic1, ic3:ic4)
                           ELSE
                              coulomb_mt3_c(:mpbasis%num_rad_bas_fun(0, itype1) - 1, iatom, iatom1, ikpt0) &
                                 = CONJG(coulhlp%data_c(ic1, ic3:ic4))
                           ENDIF
                           ic2 = ic2 + ishift1
                        END DO
                     END DO

                     ic = ic + ishift
                  END DO
               END DO

               !test
               iatom = 0
               DO itype = 1, atoms%ntype
                  DO ineq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     if (sym%invs) THEN
                        IF (MAXVAL(ABS(coulomb_mt2_r(:mpbasis%num_rad_bas_fun(0, itype) - 1, 0, 0, &
                                                     iatom, ikpt0) &
                                       - coulomb_mt3_r(:mpbasis%num_rad_bas_fun(0, itype) - 1, iatom, &
                                                       iatom, ikpt0))) &
                            > 1E-08) &
                           call judft_error('coulombmatrix: coulomb_mt2 and coulomb_mt3 are inconsistent')

                     else
                        IF (MAXVAL(ABS(coulomb_mt2_c(:mpbasis%num_rad_bas_fun(0, itype) - 1, 0, 0, &
                                                     iatom, ikpt0) &
                                       - coulomb_mt3_c(:mpbasis%num_rad_bas_fun(0, itype) - 1, iatom, &
                                                       iatom, ikpt0))) &
                            > 1E-08) &
                           call judft_error('coulombmatrix: coulomb_mt2 and coulomb_mt3 are inconsistent')
                     endif
                  END DO
               END DO
            END IF
            call timestop("gamma point treatment")

         END IF ! calc_mt

         !
         ! add the residual MT contributions, i.e. those functions with an moment,
         ! to the matrix coulomb_mtir, which is fully occupied
         !
         call timestart("residual MT contributions")
         ic = 0
         DO itype = 1, atoms%ntype
            DO ineq = 1, atoms%neq(itype)
               DO l = 0, hybrid%lcutm1(itype)
                  DO M = -l, l
                     ic = ic + 1
                  END DO
               END DO
            END DO
         END DO

         indx1 = 0; indx2 = 0; indx3 = 0; indx4 = 0
         DO itype = 1, atoms%ntype
            DO ineq = 1, atoms%neq(itype)
               DO l = 0, hybrid%lcutm1(itype)
                  DO M = -l, l
                     indx1 = indx1 + 1
                     indx3 = indx3 + mpbasis%num_rad_bas_fun(l, itype)

                     indx2 = 0
                     indx4 = 0
                     DO itype1 = 1, atoms%ntype
                        DO ineq1 = 1, atoms%neq(itype1)
                           DO l1 = 0, hybrid%lcutm1(itype1)
                              DO m1 = -l1, l1
                                 indx2 = indx2 + 1
                                 indx4 = indx4 + mpbasis%num_rad_bas_fun(l1, itype1)
                                 IF (indx4 < indx3) CYCLE
                                 IF (calc_mt(ikpt)) THEN
                                    IF (sym%invs) THEN
                                       coulomb_mtir_r(indx1, indx2, ikpt1) = coulhlp%data_r(indx3, indx4)
                                       coulomb_mtir_r(indx2, indx1, ikpt1) = coulomb_mtir_r(indx1, indx2, ikpt1)
                                    ELSE
                                       coulomb_mtir_c(indx1, indx2, ikpt1) = coulhlp%data_c(indx3, indx4)
                                       coulomb_mtir_c(indx2, indx1, ikpt1) = CONJG(coulomb_mtir_c(indx1, indx2, ikpt1))
                                    ENDIF
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO

                     DO igpt = 1, mpbasis%ngptm(ikpt)
                        indx2 = indx2 + 1
                        IF (sym%invs) THEN
                           coulomb_mtir_r(indx1, indx2, ikpt1) = coulhlp%data_r(indx3, hybrid%nbasp + igpt)
                           coulomb_mtir_r(indx2, indx1, ikpt1) = coulomb_mtir_r(indx1, indx2, ikpt1)
                        ELSE
                           coulomb_mtir_c(indx1, indx2, ikpt1) = coulhlp%data_c(indx3, hybrid%nbasp + igpt)
                           coulomb_mtir_c(indx2, indx1, ikpt1) = CONJG(coulomb_mtir_c(indx1, indx2, ikpt1))
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
         if (sym%invs) THEN
            coulomb_mtir_r(ic + 1:ic + mpbasis%ngptm(ikpt), ic + 1:ic + mpbasis%ngptm(ikpt), ikpt1) &
               = coulhlp%data_r(hybrid%nbasp + 1:nbasm1(ikpt), hybrid%nbasp + 1:nbasm1(ikpt))
            ic2 = indx1 + mpbasis%ngptm(ikpt)
            coulombp_mtir_r(:ic2*(ic2 + 1)/2, ikpt0) = packmat(coulomb_mtir_r(:ic2, :ic2, ikpt1))
         else
            coulomb_mtir_c(ic + 1:ic + mpbasis%ngptm(ikpt), ic + 1:ic + mpbasis%ngptm(ikpt), ikpt1) &
               = coulhlp%data_c(hybrid%nbasp + 1:nbasm1(ikpt), hybrid%nbasp + 1:nbasm1(ikpt))
            ic2 = indx1 + mpbasis%ngptm(ikpt)
            coulombp_mtir_c(:ic2*(ic2 + 1)/2, ikpt0) = packmat(coulomb_mtir_c(:ic2, :ic2, ikpt1))
         end if
         call timestart("write coulomb_spm")
         if (sym%invs) THEN
            CALL write_coulomb_spm_r(ikpt, coulomb_mt1(:, :, :, :, 1), coulomb_mt2_r(:, :, :, :, 1), &
                                     coulomb_mt3_r(:, :, :, 1), coulombp_mtir_r(:, 1))
!!$       print *,"DEBUG"
!!$       DO n1=1,SIZE(coulomb_mt1,1)
!!$          DO n2=1,SIZE(coulomb_mt1,2)
!!$             DO i=1,SIZE(coulomb_mt1,3)
!!$                DO j=1,SIZE(coulomb_mt1,4)
!!$                   WRITE(732,*) n1,n2,i-1,j,coulomb_mt2_r(n1,n2,i-1,j,1)
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
         else
            call write_coulomb_spm_c(ikpt, coulomb_mt1(:, :, :, :, 1), coulomb_mt2_c(:, :, :, :, 1), &
                                     coulomb_mt3_c(:, :, :, 1), coulombp_mtir_c(:, 1))
         endif
         call timestop("write coulomb_spm")

      END DO ! ikpt
      call timestop("loop bla")

      if (sym%invs) THEN
         deallocate(coulomb_mt1, coulomb_mt2_r, coulomb_mt3_r, coulomb_mtir_r, coulombp_mtir_r)
      else
         deallocate(coulomb_mt1, coulomb_mt2_c, coulomb_mt3_c, coulomb_mtir_c, coulombp_mtir_c)
      end if

      IF (mpi%irank == 0) THEN
         WRITE (6, '(7X,A)', advance='no') 'done'
         CALL cpu_TIME(time2)
         WRITE (6, '(2X,A,F8.2,A)') '( Timing:', time2 - time1, ' )'
      END IF

      CALL timestop("Coulomb matrix setup")

   END SUBROUTINE coulombmatrix

   !     Calculate body of Coulomb matrix at Gamma point: v_IJ = SUM(G) c^*_IG c_JG 4*pi/G**2 .
   !     For this we must subtract from coulomb(:,1) the spherical average of a term that comes
   !     from the fact that MT functions have k-dependent Fourier coefficients (see script).
   SUBROUTINE subtract_sphaverage(sym, cell, atoms, mpbasis, hybrid, nbasm1, gridf, coulomb)

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
      TYPE(t_mpbasis), intent(in)  :: mpbasis
      TYPE(t_hybrid), INTENT(IN)    :: hybrid

      INTEGER, INTENT(IN)    :: nbasm1(:)
      REAL, INTENT(IN)    :: gridf(:,:)
      COMPLEX, INTENT(INOUT) :: coulomb(:,:)

      ! - local scalars -
      INTEGER               :: l, i, j, n, nn, itype, ieq, M

      ! - local arrays -
      TYPE(t_mat) :: olap
      !COMPLEX , ALLOCATABLE :: constfunc(:)  !can also be real in inversion case
      COMPLEX      :: coeff(nbasm1(1)), cderiv(nbasm1(1), -1:1), claplace(nbasm1(1))

      CALL olap%alloc(sym%invs, mpbasis%ngptm(1), mpbasis%ngptm(1), 0.)

      n = nbasm1(1)
      nn = n*(n + 1)/2
      CALL olap_pw(olap, mpbasis%gptm(:, mpbasis%gptm_ptr(:mpbasis%ngptm(1), 1)), mpbasis%ngptm(1), atoms, cell)

      ! Define coefficients (coeff) and their derivatives (cderiv,claplace)
      coeff = 0
      cderiv = 0
      claplace = 0
      j = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            DO l = 0, hybrid%lcutm1(itype)
               DO M = -l, l
                  DO i = 1, mpbasis%num_rad_bas_fun(l, itype)
                     j = j + 1
                     IF (l == 0) THEN
                        coeff(j) = SQRT(4*pi_const) &
                                   *intgrf(atoms%rmsh(:, itype)*hybrid%basm1(:, i, 0, itype), &
                                           atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf) &
                                   /SQRT(cell%vol)

                        claplace(j) = -SQRT(4*pi_const) &
                                      *intgrf(atoms%rmsh(:, itype)**3*hybrid%basm1(:, i, 0, itype), &
                                              atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf) &
                                      /SQRT(cell%vol)

                     ELSE IF (l == 1) THEN
                        cderiv(j, M) = -SQRT(4*pi_const/3)*CMPLX(0.0, 1.0) &
                                       *intgrf(atoms%rmsh(:, itype)**2*hybrid%basm1(:, i, 1, itype), &
                                               atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf) &
                                       /SQRT(cell%vol)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      IF (olap%l_real) THEN
         coeff(hybrid%nbasp + 1:n) = olap%data_r(1, 1:n - hybrid%nbasp)
      else
         coeff(hybrid%nbasp + 1:n) = olap%data_c(1, 1:n - hybrid%nbasp)
      END IF
      IF (sym%invs) THEN
         CALL symmetrize(coeff, 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                         mpbasis%num_rad_bas_fun, sym)
         CALL symmetrize(claplace, 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                         mpbasis%num_rad_bas_fun, sym)
         CALL symmetrize(cderiv(:, -1), 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                         mpbasis%num_rad_bas_fun, sym)
         CALL symmetrize(cderiv(:, 0), 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                         mpbasis%num_rad_bas_fun, sym)
         CALL symmetrize(cderiv(:, 1), 1, nbasm1(1), 2, .FALSE., &
                         atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                         mpbasis%num_rad_bas_fun, sym)
      ENDIF
      ! Subtract head contributions from coulomb(:nn,1) to obtain the body
      l = 0
      DO j = 1, n
         DO i = 1, j
            l = l + 1
            coulomb(l, 1) = coulomb(l, 1) - 4*pi_const/3 &
                            *(dot_PRODUCT(cderiv(i, :), cderiv(j, :)) &
                              + (CONJG(coeff(i))*claplace(j) &
                                 + CONJG(claplace(i))*coeff(j))/2)
         END DO
      END DO
      coeff(hybrid%nbasp + 1) = 1.0
      coeff(hybrid%nbasp + 2:) = 0.0
      IF (sym%invs) THEN

         CALL desymmetrize(coeff, 1, nbasm1(1), 2, &
                           atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                           mpbasis%num_rad_bas_fun, sym)
         CALL symmetrize(coeff, nbasm1(1), 1, 1, .FALSE., &
                         atoms, hybrid%lcutm1, maxval(hybrid%lcutm1), &
                         mpbasis%num_rad_bas_fun, sym)
      ENDIF
      ! Explicit normalization here in order to prevent failure of the diagonalization in diagonalize_coulomb
      ! due to inaccuracies in the overlap matrix (which can make it singular).
      !constfunc = coeff / SQRT ( ( SUM(ABS(coeff(:hybrid%nbasp))**2) + dot_product ( coeff(hybrid%nbasp+1:), MATMUL(olap,coeff(hybrid%nbasp+1:)) ) ) )

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

   ! Convergence parameter
#define CONVPARAM 1e-18
   ! Do some additional shells ( real-space and Fourier-space sum )
#define ADDSHELL1 40
#define ADDSHELL2 0

   SUBROUTINE structureconstant(structconst, cell, hybrid, atoms, kpts, mpi)

      USE m_constants, ONLY: pi_const
      USE m_rorder, ONLY: rorderp, rorderpf
      USE m_types
      USE m_juDFT
      use m_ylm
      IMPLICIT NONE

      TYPE(t_mpi), INTENT(IN)   :: mpi

      TYPE(t_hybrid), INTENT(INOUT)   :: hybrid

      TYPE(t_cell), INTENT(IN)   :: cell

      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_kpts), INTENT(IN)    :: kpts
      ! - scalars -

      ! - arrays -
      COMPLEX, INTENT(INOUT)   ::  structconst(:,:,:,:)

      ! - local scalars -
      INTEGER                   ::  i, ic1, ic2, lm, ikpt, l, ishell, nshell
      INTEGER                   ::  m
      INTEGER                   ::  nptsh, maxl

      REAL                      ::  rad, rrad, rdum
      REAL                      ::  a, a1, aa
      REAL                      ::  pref, rexp
      REAL                      ::  time1, time2
      REAL                      ::  scale

      COMPLEX                   ::  cdum, cexp

      LOGICAL, SAVE          ::  first = .TRUE.
      ! - local arrays -
      INTEGER                   ::  conv(0:2*hybrid%lexp)
      INTEGER, ALLOCATABLE     ::  pnt(:), ptsh(:, :)

      REAL                      ::  rc(3), ra(3), k(3), ki(3), ka(3)
      REAL                      ::  convpar(0:2*hybrid%lexp), g(0:2*hybrid%lexp)
      REAL, ALLOCATABLE     ::  radsh(:)

      COMPLEX                   ::  y((2*hybrid%lexp + 1)**2)
      COMPLEX                   ::  shlp((2*hybrid%lexp + 1)**2, kpts%nkpt)

      call timestart("calc struc_const.")

      IF (mpi%irank /= 0) first = .FALSE.

      rdum = cell%vol**(1.0/3) ! define "average lattice parameter"

      ! ewaldlambda = ewaldscale
      scale = hybrid%ewaldlambda/rdum

      !       lambda = ewaldlambda / rdum

      pref = 4*pi_const/(scale**3*cell%vol)

      DO l = 0, 2*hybrid%lexp
         convpar(l) = CONVPARAM/scale**(l + 1)
      END DO

      IF (first) THEN
         WRITE (6, '(//A)') '### subroutine: structureconstant ###'
         WRITE (6, '(/A)') 'Real-space sum:'
      END IF

      !
      !     Determine cutoff radii for real-space and Fourier-space summation
      ! (1) real space
      call timestart("determine cutoff radii")
      a = 1
      g = 1e99
      do while(ANY(g > convpar/10))
         rexp = EXP(-a)
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
         DO l = 8, 2*hybrid%lexp
            g(l) = a**(-l - 1)
         END DO
         a = a + 1
      enddo

      rad = a/scale
      call timestop("determine cutoff radii")

      ! (2) Fourier space
      call timestart("fourier space")
      a = 1
      g(0:7) = 1e99
      do while(ANY(g > convpar))
         aa = (1 + a**2)**(-1)
         g(0) = pref*aa**4/a**2
         g(1) = pref*aa**4/a
         g(2) = pref*aa**5/3
         g(3) = pref*aa**5*a/15
         g(4) = pref*aa**6*a**2/105
         g(5) = pref*aa**6*a**3/945
         g(6) = pref*aa**7*a**4/10395
         g(7) = pref*aa**7*a**5/135135
         a = a + 1
      enddo
      rrad = a*scale
      call timestop("fourier space")

      IF (first) THEN
         WRITE (6, '(/A,2F10.5)') 'Cutoff radii: ', rad, rrad
         WRITE (6, '(/A)') 'Real-space sum'
      END IF

      !
      !     Determine atomic shells
      call timestart("determine atomic shell")
      CALL getshells(ptsh, nptsh, radsh, nshell, rad, cell%amat, first)
      call timestop("determine atomic shell")

      allocate(pnt(nptsh))
      structconst = 0

      !
      !     Real-space sum
      !
      CALL cpu_TIME(time1)

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
            maxl = 2*hybrid%lexp
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
               call timestart("harmonics")
               call ylm4(maxl, ra, y)
               call timestop("harmonics")
               y = CONJG(y)
               call timestart("kloop")
               DO ikpt = 1, kpts%nkpt
                  rdum = kpts%bk(1, ikpt)*ptsh(1, i) + kpts%bk(2, ikpt)*ptsh(2, i) + kpts%bk(3, ikpt)*ptsh(3, i)
                  cexp = EXP(CMPLX(0.0, 1.0)*2*pi_const*rdum)
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
               call timestop("kloop")
            END DO
            structconst(:, ic1, ic2, :) = shlp
         END DO
      END DO
      call timestop("realspace sum")

      deallocate(ptsh, radsh)

      CALL cpu_TIME(time2)
      IF (first) WRITE (6, '(A,F7.2)') '  Timing: ', time2 - time1
      CALL cpu_TIME(time1)

      IF (first) WRITE (6, '(/A)') 'Fourier-space sum'

      !
      !     Determine reciprocal shells
      !
      call timestart("determince reciproc. shell")
      CALL getshells(ptsh, nptsh, radsh, nshell, rrad, cell%bmat, first)
      call timestop("determince reciproc. shell")
      ! minimum nonzero reciprocal-shell radius (needed in routines concerning the non-local hartree-fock exchange)
      !hybrid%radshmin = radsh(2)
      !
      !     Fourier-space sum
      !
      call timestart("fourierspace sum")
      DO ikpt = 1, kpts%nkpt
         k = kpts%bk(:, ikpt)
         maxl = MIN(7, hybrid%lexp*2)
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
            call timestart("harmonics")
            call ylm4(maxl, ka, y)
            IF(norm2(ka(:)).LT.1.0e-16) y(2:(maxl+1)**2)=CMPLX(0.0,0.0)
            call timestop("harmonics")
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
                  cexp = EXP(CMPLX(0.0, 1.0)*2*pi_const*dot_PRODUCT(ki, atoms%taual(:, ic1) - atoms%taual(:, ic2)))
                  DO lm = 1, (maxl + 1)**2
                     structconst(lm, ic1, ic2, ikpt) = structconst(lm, ic1, ic2, ikpt) + cexp*y(lm)
                  END DO
               END DO
            END DO
         END DO
      END DO
      call timestop("fourierspace sum")

      CALL cpu_TIME(time2)
      IF (first) WRITE (6, '(A,F7.2)') '  Timing: ', time2 - time1

      !
      !     Add contribution for l=0 to diagonal elements and rescale structure constants
      !
      structconst(1, 1, 1, :) = structconst(1, 1, 1, :) - 5.0/16/SQRT(4*pi_const)
      DO i = 2, atoms%nat
         structconst(:, i, i, :) = structconst(:, 1, 1, :)
      END DO
      DO l = 0, 2*hybrid%lexp
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
                           (structconst(1, ic1, ic2, 1) + SQRT(4*pi_const)/cell%vol/rdum**2* &
                            EXP(-CMPLX(0.0, 1.0)*2*pi_const*dot_PRODUCT( &
                                kpts%bk(:, ikpt), atoms%taual(:, ic2) - atoms%taual(:, ic1)))))**2
            END DO
         END DO
         a = SQRT(a/atoms%nat**2)
         aa = SQRT(SUM(ABS(structconst(1, :, :, ikpt))**2)/atoms%nat**2)
         IF (first) WRITE (6, '(/A,F8.5,A,F8.5,A)') 'Accuracy of Gamma-decomposition (structureconstant):', a, ' (abs)', a/aa, ' (rel)'
      ENDIF
      deallocate(ptsh, radsh)

      first = .FALSE.

      call timestop("calc struc_const.")
   END SUBROUTINE structureconstant

   !     -----------------

   !     Determines all shells of the crystal defined by lat and vol with radii smaller than rad.
   !     The lattice points (number = nptsh) are stored in ptsh, their corresponding lengths (shell radii) in radsh.

   SUBROUTINE getshells(ptsh, nptsh, radsh, nshell, rad, lat, lwrite)
      USE m_rorder, ONLY: rorderpf
      USE m_juDFT
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: lwrite
      INTEGER, INTENT(OUT)   :: nptsh, nshell
      INTEGER, ALLOCATABLE   :: ptsh(:, :)
      REAL, ALLOCATABLE   :: radsh(:)
      REAL, INTENT(IN)    :: rad, lat(:,:)
      REAL                    :: r(3), rdum
      INTEGER, ALLOCATABLE   :: pnt(:)
      INTEGER                 :: n, i, ix, iy, iz, ok
      LOGICAL                 :: found
      INTEGER, ALLOCATABLE   :: ihelp(:, :)
      REAL, ALLOCATABLE   :: rhelp(:)

      allocate(ptsh(3, 100000), radsh(100000), stat=ok)
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
                     allocate(rhelp(SIZE(radsh)), ihelp(3, SIZE(ptsh, 2)), stat=ok)
                     IF (ok /= 0) call judft_error('getshells: failure allocation rhelp/ihelp')
                     rhelp = radsh
                     ihelp = ptsh
                     deallocate(radsh, ptsh)
                     allocate(radsh(SIZE(rhelp) + 100000), ptsh(3, SIZE(ihelp, 2) + 100000), stat=ok)
                     IF (ok /= 0) call judft_error('getshells: failure re-allocation ptsh/radsh')
                     radsh(1:SIZE(rhelp)) = rhelp
                     ptsh(:, 1:SIZE(ihelp, 2)) = ihelp
                     deallocate(rhelp, ihelp)
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

      allocate(pnt(nptsh))

      !reallocate radsh ptsh
      allocate(rhelp(nptsh), ihelp(3, nptsh))
      rhelp = radsh(1:nptsh)
      ihelp = ptsh(:, 1:nptsh)
      deallocate(radsh, ptsh)
      allocate(radsh(nptsh), ptsh(3, nptsh))
      radsh = rhelp
      ptsh = ihelp
      deallocate(rhelp, ihelp)

      CALL rorderpf(pnt, radsh, nptsh, MAX(0, INT(LOG(nptsh*0.001)/LOG(2.0))))
      radsh = radsh(pnt)
      ptsh = ptsh(:, pnt)
      nshell = 1
      DO i = 2, nptsh
         IF (radsh(i) - radsh(i - 1) > 1e-10) nshell = nshell + 1
      END DO

      IF (lwrite) &
         WRITE (6, '(A,F10.5,A,I7,A,I5,A)') &
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

      INTEGER, INTENT(IN)  :: ngpt(:), gpt(:,:), pgpt(:, :)!(dim,kpts%nkpt)
      REAL, ALLOCATABLE :: qnrm(:), help(:)
      INTEGER, INTENT(OUT) :: nqnrm
      INTEGER, ALLOCATABLE :: pqnrm(:, :)
      INTEGER               :: i, j, ikpt, igpt, igptp
      REAL                  :: q(3), qnorm

      allocate(qnrm(MAXVAL(ngpt)*kpts%nkpt), pqnrm(MAXVAL(ngpt), kpts%nkpt))
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

      allocate(help(nqnrm))
      help(1:nqnrm) = qnrm(1:nqnrm)
      deallocate(qnrm)
      allocate(qnrm(1:nqnrm))
      qnrm = help

   END SUBROUTINE getnorm

   FUNCTION sphbessel_integral(atoms, itype, qnrm, nqnrm, iqnrm1, iqnrm2, l, hybrid, &
                               sphbes0, l_warnin, l_warnout)

      USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN)  :: itype, nqnrm, iqnrm1, iqnrm2, l
      REAL, INTENT(IN)  :: qnrm(nqnrm), sphbes0(-1:hybrid%lexp + 2, atoms%ntype, nqnrm)
      LOGICAL, INTENT(IN), OPTIONAL  ::  l_warnin
      LOGICAL, INTENT(OUT), OPTIONAL  ::  l_warnout
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
      ELSE IF (abs(q1 -q2) < 1e-12) THEN
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
               WRITE (6, '(A,E12.5,A,E12.5,A)') 'sphbessel_integral: Warning! Formula One possibly unstable. Ratios:', &
                  r1, '(', r2, ')'
               WRITE (6, '(A,2F15.10,I4)') '                    Current qnorms and atom type:', q1, q2, itype
               l_warned = .TRUE.
            END IF
            sphbessel_integral = s**3/dq*da
         ELSE
            IF (r2 < 1e-6 .AND. l_warn) THEN
               WRITE (6, '(A,E13.5,A,E13.5,A)') 'sphbessel_integral: Warning! Formula Two possibly unstable. Ratios:', &
                  r2, '(', r1, ')'
               WRITE (6, '(A,2F15.10,I4)') '                    Current qnorms and atom type:', &
                  q1, q2, itype
               l_warned = .TRUE.
            END IF
            sphbessel_integral = s**3*dc
         END IF
      END IF

      IF (PRESENT(l_warnout)) l_warnout = l_warned

   END FUNCTION sphbessel_integral

   !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

END MODULE m_coulombmatrix
