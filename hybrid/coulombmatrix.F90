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
#ifdef CPP_MPI
   use mpi
#endif
   use m_types
   USE m_intgrf, ONLY: intgrf, intgrf_init
   use m_sphbes, only: sphbes
   use m_glob_tofrom_loc
   USE m_trafo, ONLY: symmetrize_mpimat, symmetrize, bramat_trafo
   use m_gamma_double_gpt_loop
CONTAINS

   SUBROUTINE coulombmatrix(fmpi, fi, mpdata, hybdat, xcpot)
      use m_work_package
      use m_structureconstant
      USE m_types
      USE m_types_mpimat
      use m_types_mat
      USE m_types_hybdat
      USE m_juDFT
      USE m_constants
      use m_util, only: primitivef
      USE m_hsefunctional, ONLY: change_coulombmatrix
      USE m_wrapper
      USE m_io_hybrid
      use m_ylm
      use m_calc_l_m_from_lm
      use m_calc_mpsmat
      use m_copy_coul
      use m_apply_inverse_olap
      IMPLICIT NONE

      TYPE(t_xcpot_inbuild), INTENT(IN) :: xcpot
      TYPE(t_mpi), INTENT(IN)           :: fmpi
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      ! - local scalars -
      INTEGER                    :: inviop
      INTEGER                    :: nqnrm, iqnrm, iqnrm1, iqnrm2
      INTEGER                    :: itype, l, ix, iy, iy0, i, j, lm, l1, l2, m1, m2, ineq, ikpt
      INTEGER                    :: lm1, lm2, itype1, itype2, ineq1, ineq2, n1, n2, iat2
      INTEGER                    :: ic, ic1, ic2
      INTEGER                    :: igpt, igpt1, igpt2, igptp, igptp1, igptp2
      INTEGER                    :: isym, isym1, isym2, igpt0
      INTEGER                    :: iatm1, iatm2
      INTEGER                    :: m, im
      INTEGER                    :: maxfac, ix_loc, pe, pe_ix

      LOGICAL                    :: lsym

      REAL                       :: rdum, rdum1, rdum2
      REAL                       :: svol, qnorm1, qnorm2
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

      REAL                       :: q(3), q1(3), q2(3), mtmt_term
      REAL                       :: integrand(fi%atoms%jmtd), primf1(fi%atoms%jmtd), primf2(fi%atoms%jmtd)
      REAL                       :: moment(maxval(mpdata%num_radbasfn), 0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype), &
                                    moment2(maxval(mpdata%num_radbasfn), fi%atoms%ntype)
      REAL, ALLOCATABLE   :: gmat(:, :), qnrm(:)
      REAL, ALLOCATABLE   :: sphbesmoment(:, :, :)
      REAL, ALLOCATABLE   :: sphbes0(:, :, :)
      REAL, ALLOCATABLE   :: olap(:, :, :, :), integral(:, :, :, :)
      REAL, ALLOCATABLE   :: gridf(:, :)
      REAL                       :: facA(0:MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1))
      REAL                       :: facB(0:MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1))
      REAL                       :: facC(-1:MAX(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 4*MAX(maxval(fi%hybinp%lcutm1), fi%hybinp%lexp) + 1))
      REAL    :: sphbes_var(fi%atoms%jmtd, 0:maxval(fi%hybinp%lcutm1))
      REAL    :: sphbesmoment1(fi%atoms%jmtd, 0:maxval(fi%hybinp%lcutm1))

      COMPLEX     :: y((fi%hybinp%lexp + 1)**2), smat
      COMPLEX     :: dwgn(-maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), -maxval(fi%hybinp%lcutm1):maxval(fi%hybinp%lcutm1), 0:maxval(fi%hybinp%lcutm1), fi%sym%nsym)
      COMPLEX, ALLOCATABLE   :: carr2(:, :), carr2a(:, :), carr2b(:, :)
      COMPLEX, ALLOCATABLE   :: structconst(:,:,:,:)

      INTEGER                    :: ierr
      INTEGER                    :: iatom, mtmt_idx
      TYPE(t_mat)                :: mat
      type(t_mat), allocatable   :: mtmt_repl(:)
      class(t_mat), allocatable  :: coul(:)

      CALL timestart("Coulomb matrix setup")
      call timestart("prep in coulomb")
      if (fmpi%is_root()) write (*, *) "start of coulomb calculation"

      allocate(structconst((2*fi%hybinp%lexp + 1)**2, fi%atoms%nat, fi%atoms%nat, fi%kpts%nkpt), stat=ierr)
      if(ierr /= 0) call judft_error("can't allocate structconst. error: " // int2str(ierr))

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
      CALL structureconstant(structconst, fi%cell, fi%hybinp, fi%atoms, fi%kpts, fmpi)

      IF (fmpi%irank == 0) WRITE (oUnit, '(//A)') '### subroutine: coulombmatrix ###'

      !
      !     Matrix allocation
      !

      call timestart("coulomb allocation")

      if(fmpi%n_size == 1) then
         allocate(t_mat::coul(fi%kpts%nkpt))
      else
         allocate(t_mpimat::coul(fi%kpts%nkpt))
      endif
      do ikpt = 1, fi%kpts%nkpt 
         if(any(ikpt == fmpi%k_list))then 
            call coul(ikpt)%init(.False., hybdat%nbasm(ikpt), hybdat%nbasm(ikpt), fmpi%sub_comm, .false.)
         else
            call coul(ikpt)%init(.False., 1, 1, fmpi%sub_comm, .false.)
         endif 
      enddo
      call timestop("coulomb allocation") 

      IF (fmpi%irank == 0) then
         write (oUnit,*) "Size of coulomb matrix: " //&
                            float2str(sum([(coul(fmpi%k_list(i))%size_mb(), i=1,size(fmpi%k_list))])) // " MB"
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

      call timestop("getnorm")

      call bessel_calculation(fi, fmpi, mpdata, nqnrm, gridf, qnrm, sphbesmoment, olap, integral)

      !
      !     (1) Case < MT | v | MT >
      !


      !       (1a) r,r' in same MT
      call timestart("loop 1")
      ix = 0
      iy = 0
      iy0 = 0


      call mat%alloc(.True., maxval(mpdata%num_radbasfn), maxval(mpdata%num_radbasfn))

      allocate(mtmt_repl(calc_num_mtmts(fi)))

      mtmt_idx = 0
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
                  mtmt_idx = mtmt_idx + 1
                  call mtmt_repl(mtmt_idx)%alloc(.True., mpdata%num_radbasfn(l, itype), mpdata%num_radbasfn(l, itype))

                  DO n2 = 1, mpdata%num_radbasfn(l, itype)
                     DO n1 = 1, n2
                        mtmt_repl(mtmt_idx)%data_r(n1, n2) = mat%data_r(n1, n2)
                     END DO
                  END DO
                  call mtmt_repl(mtmt_idx)%u2l()
               END DO

            END DO
         END DO
      END DO
      call mat%free()
      call timestop("loop 1")

      DO im = 1, size(fmpi%k_list)
         ikpt = fmpi%k_list(im)

         ! only the first rank handles the MT-MT part
         call timestart("MT-MT part")
         ix = 0
         ic2 = 0
         mtmt_idx = 0
         DO itype2 = 1, fi%atoms%ntype
            DO ineq2 = 1, fi%atoms%neq(itype2)
               ic2 = ic2 + 1
               lm2 = 0
               DO l2 = 0, fi%hybinp%lcutm1(itype2)
                  DO m2 = -l2, l2
                     mtmt_idx = mtmt_idx + 1
                     lm2 = lm2 + 1
                     DO n2 = 1, mpdata%num_radbasfn(l2, itype2)
                        ix = ix + 1
                        call glob_to_loc(fmpi, ix, pe, ix_loc)
                        if(fmpi%n_rank == pe) then
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
                                          rdum = (-1)**(l2 + m2)*moment(n1, l1, itype1)*moment(n2, l2, itype2)*gmat(lm1, lm2)
                                          l = l1 + l2
                                          lm = l**2 + l + m1 - m2 + 1

                                          if(itype2 /= itype1 .or. ineq2 /= ineq1 .or. l2 /= l1 .or. m2 /= m1) then 
                                             mtmt_term = 0.0
                                          else 
                                             mtmt_term = mtmt_repl(mtmt_idx)%data_r(n1, n2)
                                          endif
                                       
                                          coul(ikpt)%data_c(iy, ix_loc) = mtmt_term + EXP(CMPLX(0.0, 1.0)*tpi_const* &
                                                                     dot_PRODUCT(fi%kpts%bk(:, ikpt), &
                                                                                 fi%atoms%taual(:, ic2) - fi%atoms%taual(:, ic1))) &
                                                               *rdum*structconst(lm, ic1, ic2, ikpt)
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO lp2
                        endif ! pe = n_rank
                     END DO
                  END DO
               END DO
            END DO
         END DO

         IF (fi%sym%invs) THEN
            !symmetrize makes the Coulomb matrix real symmetric     
            CALL symmetrize_mpimat(fi, fmpi, coul(ikpt)%data_c, [1,1],[hybdat%nbasp, hybdat%nbasp], &
                                   3, .FALSE., mpdata%num_radbasfn)
         ENDIF
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
         DO im = 1, size(fmpi%k_list)
            ikpt = fmpi%k_list(im)
            call loop_over_interst(fi, hybdat, mpdata, fmpi, structconst, sphbesmoment, moment, moment2, &
                                   qnrm, facc, gmat, integral, olap, pqnrm, pgptm1, ngptm1, ikpt, coul(ikpt))

         END DO

         call timestop("loop over interst.")
         deallocate (olap, integral)

         !
         !     (3) Case < PW | v | PW >
         !
         !     (3a) r,r' everywhere; r everywhere, r' in MT; r in MT, r' everywhere

         ! Coulomb matrix, contribution (3a)
         call timestart("coulomb matrix 3a")
         DO im = 1, size(fmpi%k_list)
            ikpt = fmpi%k_list(im)

            DO igpt0 = 1, ngptm1(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
               ix = hybdat%nbasp + igpt2
               call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
               if(fmpi%n_rank == pe_ix) then
                  q2 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp2), fi%cell%bmat)
                  rdum2 = SUM(q2**2)
                  IF (abs(rdum2) > 1e-12) rdum2 = fpi_const/rdum2

                  DO igpt1 = 1, igpt2
                     igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
                     iy = hybdat%nbasp + igpt1
                     q1 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp1), fi%cell%bmat)
                     rdum1 = SUM(q1**2)
                     IF (abs(rdum1) > 1e-12) rdum1 = fpi_const/rdum1
                     smat = calc_smat_elem(fi, mpdata, igptp1, igptp2)

                     IF (ikpt == 1) THEN
                        IF (igpt1 /= 1) THEN
                           coul(1)%data_c(iy,ix_loc) = -smat*rdum1/fi%cell%vol 
                        END IF
                        IF (igpt2 /= 1) THEN
                           coul(1)%data_c(iy,ix_loc) &
                              = coul(1)%data_c(iy,ix_loc) - smat*rdum2/fi%cell%vol
                        END IF
                     ELSE
                        coul(ikpt)%data_c(iy,ix_loc) = -smat*(rdum1 + rdum2)/fi%cell%vol
                     END IF
                  END DO
                  IF (ikpt /= 1 .OR. igpt2 /= 1) THEN
                     coul(ikpt)%data_c(iy,ix_loc) = coul(ikpt)%data_c(iy,ix_loc)  + rdum2
                  END IF
               endif
            END DO
         END DO
         call timestop("coulomb matrix 3a")
         !     (3b) r,r' in different MT

         call timestart("coulomb matrix 3b")
         DO im = 1, size(fmpi%k_list)
            ikpt = fmpi%k_list(im)
            if (fmpi%is_root()) write (*, *) "coulomb pw-loop nk: ("//int2str(ikpt)//"/"//int2str(fi%kpts%nkpt)//")"
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

            !$OMP parallel default(none) private(carr2, igpt0, igpt1, igpt2, ix, iy, ix_loc, pe_ix, iatom, lm1, l1, m1)&
            !$OMP private(lm2, l2, m2, cdum, ic, lm, l, m, itype, itype2, igptp1, csum, igptp2, iqnrm1, iqnrm2, cexp)&
            !$OMP shared(fmpi, fi, ikpt, ngptm1, pgptm1, pqnrm, coul, gmat, structconst, sphbesmoment, hybdat, mpdata)&
            !$OMP shared(carr2a, carr2b)
            allocate(carr2(fi%atoms%nat, (fi%hybinp%lexp + 1)**2))               
            !$OMP do schedule(dynamic)
            DO igpt0 = 1, ngptm1(ikpt)
               igpt2 = pgptm1(igpt0, ikpt)
               ix = hybdat%nbasp + igpt2
               call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
               if(fmpi%n_rank == pe_ix) then
                  igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
                  iqnrm2 = pqnrm(igpt2, ikpt)

                  carr2 = 0

!                  call timestart("itype loops")
                  do iatom = 1,fi%atoms%nat
                     itype2 = fi%atoms%itype(iatom)
                     cexp = CONJG(carr2b(iatom, igpt2))
                     
                     DO lm1 = 1, (fi%hybinp%lexp+1)**2
                        call calc_l_m_from_lm(lm1, l1, m1)
                        do lm2 = 1, (fi%hybinp%lexp+1)**2
                           call calc_l_m_from_lm(lm2, l2, m2)
                           cdum = (-1)**(l2 + m2)*sphbesmoment(l2, itype2, iqnrm2)*cexp*carr2a(lm2, igpt2)*gmat(lm1, lm2)
                           l = l1 + l2
                           lm = l**2 + l - l1 - m2 + (m1 + l1) + 1
                           do iat2 =1,fi%atoms%nat
                              carr2(iat2, lm1) = carr2(iat2,lm1) + cdum*structconst(lm, iat2, iatom, ikpt)
                           enddo
                        enddo
                     enddo
                  end do ! iatom

!                  call timestop("itype loops")

!                  call timestart("igpt1")
                  DO igpt1 = 1, igpt2
                     iy = hybdat%nbasp + igpt1
                     igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
                     iqnrm1 = pqnrm(igpt1, ikpt)
                     csum = 0
                     do ic = 1, fi%atoms%nat
                        do lm = 1, (fi%hybinp%lexp+1)**2
                           itype = fi%atoms%itype(ic)
                           call calc_l_m_from_lm(lm, l, m)
                           cdum = carr2b(ic, igpt1)*sphbesmoment(l, itype, iqnrm1)
                           csum = csum + cdum*carr2(ic, lm)*CONJG(carr2a(lm, igpt1)) ! for coulomb
                        END DO
                     END DO
                     coul(ikpt)%data_c(iy,ix_loc) = coul(ikpt)%data_c(iy,ix_loc) + csum/fi%cell%vol
                  END DO
!                  call timestop("igpt1")
               endif ! pe_ix
            END DO !igpt0
            !$omp end do
            deallocate (carr2) 
            !$OMP end parallel
            deallocate(carr2a, carr2b)

            call timestop("loop over plane waves")
         END DO !ikpt
         call timestop("coulomb matrix 3b")

         ! check if I own the gamma point
         if(any(fmpi%k_list == 1)) then
            !     Add corrections from higher orders in (3b) to coulomb(:,1)
            ! (1) igpt1 > 1 , igpt2 > 1  (finite G vectors)
            call timestart("add corrections from higher orders")
            call gamma_double_gpt_loop(fi, fmpi, hybdat, mpdata, sphbesmoment, gmat,  ngptm1, pgptm1, pqnrm, coul(1)%data_c) 

            rdum = (fpi_const)**(1.5)/fi%cell%vol**2*gmat(1, 1)

            ! (2) igpt1 = 1 , igpt2 > 1  (first G vector vanishes, second finite)
            call timestart("igpt1=1 loop")
            iy = hybdat%nbasp + 1
            DO igpt0 = 1, ngptm1(1)
               igpt2 = pgptm1(igpt0, 1)
               IF (igpt2 /= 1) then
                  ix = hybdat%nbasp + igpt2
                  call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
                  if(fmpi%n_rank == pe_ix) then
                     iqnrm2 = pqnrm(igpt2, 1)
                     igptp2 = mpdata%gptm_ptr(igpt2, 1)
                     qnorm2 = qnrm(iqnrm2)
                     DO itype1 = 1, fi%atoms%ntype
                        DO ineq1 = 1, fi%atoms%neq(itype1)
                           ic2 = 0
                           DO itype2 = 1, fi%atoms%ntype
                              DO ineq2 = 1, fi%atoms%neq(itype2)
                                 ic2 = ic2 + 1
                                 cdum = EXP(CMPLX(0.0, 1.0)*tpi_const*dot_PRODUCT(mpdata%g(:, igptp2), fi%atoms%taual(:, ic2)))
                                 coul(1)%data_c(iy, ix_loc) = coul(1)%data_c(iy, ix_loc) &
                                                   + rdum*cdum*fi%atoms%rmt(itype1)**3*( &
                                                   +sphbesmoment(0, itype2, iqnrm2)/30*fi%atoms%rmt(itype1)**2 &
                                                   - sphbesmoment(2, itype2, iqnrm2)/18 &
                                                   + sphbesmoment(1, itype2, iqnrm2)/6/qnorm2)
                              END DO
                           END DO
                        END DO
                     END DO
                  endif !pe_ix
               endif
            END DO
            call timestop("igpt1=1 loop")

            ! (2) igpt1 = 1 , igpt2 = 1  (vanishing G vectors)
            call timestart("igpt1=igpt2=1 loop")
            iy = hybdat%nbasp + 1
            ix = hybdat%nbasp + 1
            call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
            if(pe_ix == fmpi%n_rank) then
               DO itype1 = 1, fi%atoms%ntype
                  DO ineq1 = 1, fi%atoms%neq(itype1)
                     DO itype2 = 1, fi%atoms%ntype
                        DO ineq2 = 1, fi%atoms%neq(itype2)
                           coul(1)%data_c(iy, ix_loc) = coul(1)%data_c(iy, ix_loc) &
                                             + rdum*fi%atoms%rmt(itype1)**3*fi%atoms%rmt(itype2)**3* &
                                             (fi%atoms%rmt(itype1)**2 + fi%atoms%rmt(itype2)**2)/90
                        END DO
                     END DO
                  END DO
               END DO
            endif ! pe_ix
            call timestop("igpt1=igpt2=1 loop")
            call timestop("add corrections from higher orders")
         endif

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
         DO im = 1, size(fmpi%k_list)
            ikpt = fmpi%k_list(im)
            call timestart("harmonics setup")
            DO igpt = 1, mpdata%n_g(ikpt)
               igptp = mpdata%gptm_ptr(igpt, ikpt)
               q = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)
               call ylm4(fi%hybinp%lexp, q, carr2(:, igpt))
            END DO
            call timestop("harmonics setup")
            call perform_double_g_loop(fi, hybdat, fmpi, mpdata, sphbes0, carr2, ngptm1,pgptm1,&
                                       pqnrm,qnrm, nqnrm, ikpt, coul(ikpt))
            ! this one is needed
#ifdef CPP_MPI
            call timestart("post dblgloop barrier")
            call MPI_Barrier(fmpi%sub_comm, ierr)
            call timestop("post dblgloop barrier")
#endif
            call coul(ikpt)%u2l()
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
         DO im = 1, size(fmpi%k_list)
            ikpt = fmpi%k_list(im)
            carr2 = 0; iarr = 0
            iarr(pgptm1(:ngptm1(ikpt), ikpt)) = 1
            DO igpt0 = 1, ngptm1(ikpt)
               lsym = (1 <= igpt0) .AND. (ngptm1(ikpt) >= igpt0)
               igpt2 = pgptm1(igpt0, ikpt)
               ix = hybdat%nbasp + igpt2
               call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
               if(pe_ix == fmpi%n_rank) then 
                  carr2(:hybdat%nbasm(ikpt),2) = coul(ikpt)%data_c(:hybdat%nbasm(ikpt),ix_loc)
               endif
#ifdef CPP_MPI
               call timestart("bcast carr2")
               call MPI_Bcast(carr2(1,2), hybdat%nbasm(ikpt), MPI_DOUBLE_COMPLEX, pe_ix, fmpi%sub_comm, ierr)
               call timestop("bcast carr2")
#endif

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
                     ix = hybdat%nbasp + igpt1
                     call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
                     if(pe_ix == fmpi%n_rank) then
                        coul(ikpt)%data_c(:hybdat%nbasp + igpt1,ix_loc) = carr2(:hybdat%nbasp + igpt1, 1)
                     endif 

                     do ix = 1,hybdat%nbasp + igpt1
                        call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
                        if(pe_ix == fmpi%n_rank) coul(ikpt)%data_c(hybdat%nbasp + igpt1, ix_loc) = conjg(carr2(ix, 1))
                     enddo

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
         deallocate (carr2, iarr, pgptm1)
      END IF
      deallocate (qnrm, pqnrm)

      IF (xcpot%is_name("hse") .OR. xcpot%is_name("vhse")) THEN
         !
         ! The HSE functional is realized subtracting erf/r from
         ! the normal Coulomb matrix
         !
      ELSE
         ! check for gamma
         if(any(fmpi%k_list == 1)) then
            CALL subtract_sphaverage(fi%sym, fi%cell, fi%atoms, mpdata, &
                                    fi%hybinp, hybdat, fmpi, hybdat%nbasm, gridf, coul(1))
         endif
      END IF
      
      ! transform Coulomb matrix to the biorthogonal set
      call timestop("gap 1:")
      DO im = 1, size(fmpi%k_list)
         ikpt = fmpi%k_list(im)
         call apply_inverse_olaps(mpdata, fi%atoms, fi%cell, hybdat, fmpi, fi%sym, ikpt, coul(ikpt))
         ! lower to upper, because the lower half is better in memory
#ifdef CPP_MPI
         call timestart("post inverse barrier")
         call MPI_BARRIER(fmpi%sub_comm, ierr)
         call timestop("post inverse barrier")
#endif
         call coul(ikpt)%l2u()
      enddo

      !call plot_coulombmatrix() -> code was shifted to plot_coulombmatrix.F90
      !
      ! rearrange coulomb matrix
      !
      if(.not. allocated(hybdat%coul)) allocate(hybdat%coul(fi%kpts%nkpt))

      DO im = 1, size(fmpi%k_list)
         ikpt = fmpi%k_list(im)
         ! unpack coulomb into coulomb(ikpt)
         call copy_from_dense_to_sparse(fi, fmpi, mpdata, coul, ikpt, hybdat)
      END DO ! ikpt
      CALL timestop("Coulomb matrix setup")

   END SUBROUTINE coulombmatrix

   !     Calculate body of Coulomb matrix at Gamma point: v_IJ = SUM(G) c^*_IG c_JG 4*pi/G**2 .
   !     For this we must subtract from coulomb(:,1) the spherical average of a term that comes
   !     from the fact that MT functions have k-dependent Fourier coefficients (see script).
   SUBROUTINE subtract_sphaverage(sym, cell, atoms, mpdata, hybinp, hybdat, fmpi, nbasm1, gridf, coulomb)

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
      type(t_mpi), intent(in)    :: fmpi

      INTEGER, INTENT(IN)    :: nbasm1(:)
      REAL, INTENT(IN)    :: gridf(:, :)
      class(t_mat), intent(inout) :: coulomb

      ! - local scalars -
      INTEGER               :: l, ix, iy, ix_loc, pe_ix, i, j, n, nn, itype, ieq, M, ierr

      ! - local arrays -
      TYPE(t_mat) :: olap
      !COMPLEX , ALLOCATABLE :: constfunc(:)  !can also be real in inversion case
      COMPLEX      :: coeff(1,nbasm1(1)), cderiv(-1:1, nbasm1(1)), claplace(1,nbasm1(1))

      call timestart("subtract_sphaverage")
      CALL olap%alloc(sym%invs, mpdata%n_g(1), mpdata%n_g(1), 0.)

      n = nbasm1(1)
      nn = n*(n + 1)/2
      CALL olap_pw(olap, mpdata%g(:, mpdata%gptm_ptr(:mpdata%n_g(1), 1)), mpdata%n_g(1), atoms, cell, fmpi)

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
                        coeff(1,j) = SQRT(fpi_const) &
                                   *intgrf(atoms%rmsh(:, itype)*mpdata%radbasfn_mt(:, i, 0, itype), &
                                           atoms, itype, gridf) &
                                   /SQRT(cell%vol)

                        claplace(1,j) = -SQRT(fpi_const) &
                                      *intgrf(atoms%rmsh(:, itype)**3*mpdata%radbasfn_mt(:, i, 0, itype), &
                                              atoms, itype, gridf) &
                                      /SQRT(cell%vol)

                     ELSE IF (l == 1) THEN
                        cderiv(M,j) = -SQRT(fpi_const/3)*CMPLX(0.0, 1.0) &
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
         coeff(1,hybdat%nbasp + 1:n) = olap%data_r(1, 1:n - hybdat%nbasp)
      else
         coeff(1,hybdat%nbasp + 1:n) = olap%data_c(1, 1:n - hybdat%nbasp)
      END IF
      IF (sym%invs) THEN
         CALL symmetrize(coeff, 1, nbasm1(1), 2, &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(claplace, 1, nbasm1(1), 2, &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(cderiv(-1:-1,:), 1, nbasm1(1), 2, &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(cderiv(0:0,:), 1, nbasm1(1), 2, &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
         CALL symmetrize(cderiv(1:1,:), 1, nbasm1(1), 2, &
                         atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                         mpdata%num_radbasfn, sym)
      ENDIF
      ! Subtract head contributions from coulomb(:nn,1) to obtain the body
      l = 0
      DO ix = 1, n
         call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
         if(fmpi%n_rank == pe_ix) then
            DO iy = 1, ix
               l = l + 1
               coulomb%data_c(iy,ix_loc) = coulomb%data_c(iy,ix_loc) - fpi_const/3 &
                                          *(dot_PRODUCT(cderiv(:,iy), cderiv(:,ix)) &
                                          + (CONJG(coeff(1,iy))*claplace(1,ix) &
                                             + CONJG(claplace(1,iy))*coeff(1,ix))/2)
            END DO
         endif
      END DO

      !needed bc apply inverse uses lower half 
#ifdef CPP_MPI
      call timestart("post subtr avg barrier")
      call MPI_Barrier(fmpi%sub_comm, ierr)
      call timestop("post subtr avg barrier")
#endif
      call coulomb%u2l()
      call timestop("subtract_sphaverage")
   END SUBROUTINE subtract_sphaverage



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

      allocate (qnrm(MAXVAL(ngpt)*kpts%nkpt), source=0.0)
      allocate (pqnrm(MAXVAL(ngpt), kpts%nkpt), source=0)
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

   subroutine loop_over_interst(fi, hybdat, mpdata, fmpi, structconst, sphbesmoment, moment, moment2, &
                                qnrm, facc, gmat, integral, olap, pqnrm, pgptm1, ngptm1, ikpt, coul)
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
      type(t_mpi), intent(in)           :: fmpi
      REAL, intent(in)                  :: sphbesmoment(0:, :, :), qnrm(:), facC(-1:), gmat(:, :), moment(:, 0:, :), moment2(:, :)
      real, intent(in)                  :: integral(:, 0:, :, :), olap(:, 0:, :, :)
      integer, intent(in)               :: ikpt, ngptm1(:), pqnrm(:, :), pgptm1(:, :)
      complex, intent(in)               :: structconst(:, :, :, :)
      class(t_mat), intent(inout)       :: coul

      integer  :: igpt0, igpt, igptp, iqnrm, niter
      integer  :: ix, iy, ic, itype, lm, l, m, itype1, ic1, l1, m1, lm1, loc_from
      integer  :: l2, m2, lm2, n, i, iatm, j_type, j_l, iy_start, j_m, j_lm, pe_ix, ix_loc
      real     :: q(3), qnorm, svol, tmp_vec(3)
      COMPLEX  :: y((fi%hybinp%lexp + 1)**2), y1((fi%hybinp%lexp + 1)**2), y2((fi%hybinp%lexp + 1)**2)
      complex  :: csum, csumf(9), cdum, cexp
      integer, allocatable :: lm_arr(:), ic_arr(:)


      call range_from_glob_to_loc(fmpi, hybdat%nbasp+1, loc_from)
      coul%data_c(:hybdat%nbasp,loc_from:) = 0 

      svol = SQRT(fi%cell%vol)
      ! start to loop over interstitial plane waves
      !DO igpt0 = 1, ngptm1(ikpt)
      do igpt0 = 1, ngptm1(ikpt)
         igpt = pgptm1(igpt0, ikpt)
         igptp = mpdata%gptm_ptr(igpt, ikpt)
         ix = hybdat%nbasp + igpt
         call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
         if(pe_ix == fmpi%n_rank) then 
            q = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)
            qnorm = norm2(q)
            iqnrm = pqnrm(igpt, ikpt)
            IF (ABS(qnrm(iqnrm) - qnorm) > 1e-12) then
               call judft_error('coulombmatrix: qnorm does not equal corresponding & element in qnrm (bug?)') ! We shouldn't st op here!
            endif

            tmp_vec = MATMUL(fi%kpts%bk(:, fi%kpts%nkpt), fi%cell%bmat)
            call ylm4(2, tmp_vec, y1)
            tmp_vec = MATMUL(mpdata%g(:, igptp), fi%cell%bmat)
            call ylm4(2, tmp_vec, y2)
            call ylm4(fi%hybinp%lexp, q, y)
            y1 = CONJG(y1); y2 = CONJG(y2); y = CONJG(y)

            ! this unrolls the do ic=1,atoms%nat{do lm=1,..{}} 
            call collapse_ic_and_lm_loop(fi%atoms, fi%hybinp%lcutm1, niter, ic_arr, lm_arr)

            !$OMP PARALLEL DO default(none) &
            !$OMP private(ic, lm, itype, l, m, csum, csumf, ic1, itype1, cexp, lm1, l2, cdum, m2, lm2, iy) &
            !$OMP private(j_m, j_type, iy_start, l1, m1) &
            !$OMP shared(ic_arr, lm_arr, fi, mpdata, olap, qnorm, moment, integral, hybdat, svol) &
            !$OMP shared(moment2, ix, igpt, facc, structconst, y, y1, y2, gmat, iqnrm, sphbesmoment, ikpt) &
            !$OMP shared(igptp, niter, fmpi, pe_ix, coul, ix_loc) &
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
                     IF (l == 0) coul%data_c(iy, ix_loc) = -cdum*moment2(n, itype)/6/svol 
                     coul%data_c(iy, ix_loc) = coul%data_c(iy, ix_loc) &
                                                      + (-cdum/(2*l + 1)*integral(n, l, itype, iqnrm) & ! (2b)&
                                                         + csum*moment(n, l, itype))/svol          ! (2c)
                  ELSE
                     coul%data_c(iy, ix_loc) = (cdum*olap(n, l, itype, iqnrm)/qnorm**2 &  ! (2a)&
                                                - cdum/(2*l + 1)*integral(n, l, itype, iqnrm) & ! (2b)&
                                                + csum*moment(n, l, itype))/svol          ! (2c)
                  endif
               END DO
            END DO ! collapsed atom & lm loop (ic)
            !$OMP END PARALLEL DO
         endif !pe_ix
      END DO

      IF (fi%sym%invs) THEN
         call symmetrize_mpimat(fi, fmpi, coul%data_c, [1,hybdat%nbasp+1], [hybdat%nbasp, hybdat%nbasp+mpdata%n_g(ikpt)],&
                                1, .false., mpdata%num_radbasfn)
      ENDIF
   endsubroutine loop_over_interst

   subroutine perform_double_g_loop(fi, hybdat, fmpi, mpdata, sphbes0, carr2, ngptm1,pgptm1,pqnrm,qnrm, nqnrm, ikpt, coulomb)
      use m_juDFT
      use m_constants, only: tpi_const,fpi_const
      use m_sphbessel_integral
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      type(t_mpi), intent(in)           :: fmpi
      integer, intent(in)               :: ikpt, ngptm1(:), pqnrm(:,:),pgptm1(:, :), nqnrm
      real, intent(in)                  :: qnrm(:), sphbes0(:, :, :)
      complex, intent(in)               :: carr2(:, :)
      !complex, intent(inout)            :: coulomb(:) ! only at ikpt
      class(t_mat), intent(inout)        :: coulomb

      integer :: igpt0, igpt1, igpt2, ix, iy, igptp1, igptp2, iqnrm1, iqnrm2
      integer :: ic, itype, lm, m, l, pe_ix, ix_loc
      real    :: q1(3), q2(3)
      complex :: y1((fi%hybinp%lexp + 1)**2), y2((fi%hybinp%lexp + 1)**2)
      COMPLEX :: cexp1(fi%atoms%ntype)
      complex :: cdum, cdum1 
      logical :: ldum

      call timestart("double g-loop")

      DO igpt0 = 1, ngptm1(ikpt)
         igpt2 = pgptm1(igpt0, ikpt)
         ix = hybdat%nbasp + igpt2
         call glob_to_loc(fmpi, ix, pe_ix, ix_loc)
         if(pe_ix == fmpi%n_rank) then
            igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
            iqnrm2 = pqnrm(igpt2, ikpt)
            q2 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp2), fi%cell%bmat)
            y2 = CONJG(carr2(:, igpt2))
            
            !$OMP PARALLEL DO default(none) &
            !$OMP private(igpt1, iy, igptp1, iqnrm1, q1, y1, cexp1, ic, itype, lm) &
            !$OMP private(cdum, l, cdum1, m, ldum) &
            !$OMP shared(igpt2, coulomb, hybdat, mpdata, ikpt, fi, carr2, pqnrm, igptp2)&
            !$OMP shared(qnrm, sphbes0, iqnrm2, nqnrm, y2, ix_loc)
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
               coulomb%data_c(iy,ix_loc) = coulomb%data_c(iy,ix_loc) + (fpi_const)**3*cdum/fi%cell%vol
            END DO
            !$OMP end parallel do
         endif !pe_ix
      END DO
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

   subroutine bessel_calculation(fi, fmpi, mpdata, nqnrm, gridf, qnrm, sphbesmoment, olap, integral)
      implicit NONE 
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpi), intent(in)           :: fmpi
      type(t_mpdata), intent(in)        :: mpdata
      integer, intent(in)               :: nqnrm
      real, intent(in)                  :: gridf(:,:), qnrm(:)
      real, intent(inout)               :: sphbesmoment(0:,:,:), olap(:,0:,:,:), integral(:,0:,:,:)

      integer :: iqnrm, itype, i, l, n, ng, buf_sz, root, ierr
      REAL    :: sphbes_var(fi%atoms%jmtd, 0:maxval(fi%hybinp%lcutm1))
      REAL    :: sphbesmoment1(fi%atoms%jmtd, 0:maxval(fi%hybinp%lcutm1))
      REAL    :: rarr(0:fi%hybinp%lexp + 1), rarr1(0:maxval(fi%hybinp%lcutm1))
      real    :: rdum, qnorm, rdum1

      call timestart("Bessel calculation")
      
      do iqnrm = 1+fmpi%irank, nqnrm, fmpi%isize
         qnorm = qnrm(iqnrm)
         !$OMP parallel do default(none) &
         !$OMP shared(olap, integral, sphbesmoment, fi,qnorm, iqnrm, mpdata, gridf) &
         !$OMP private(itype, rdum, sphbes_var, sphbesmoment1, ng, rarr, rarr1, rdum1, i, l)
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

#ifdef CPP_MPI 
      call timestart("bcast bessel")
      do iqnrm = 1, nqnrm
         root = mod(iqnrm - 1,fmpi%isize)
         buf_sz = size(olap,1) * size(olap,2) * size(olap,3)
         call MPI_Bcast(olap(1,0,1,iqnrm), buf_sz, MPI_DOUBLE_PRECISION, root, fmpi%mpi_comm, ierr)

         buf_sz = size(integral,1) * size(integral,2) * size(integral,3)
         call MPI_Bcast(integral(1,0,1,iqnrm), buf_sz, MPI_DOUBLE_PRECISION, root, fmpi%mpi_comm, ierr)

         buf_sz = size(sphbesmoment,1) * size(sphbesmoment,2)
         call MPI_Bcast(sphbesmoment(0,1,iqnrm), buf_sz, MPI_DOUBLE_PRECISION, root, fmpi%mpi_comm, ierr)
      enddo
      call timestop("bcast bessel")
#endif
      call timestop("Bessel calculation")
   end subroutine bessel_calculation

   function calc_num_mtmts(fi) result(num_mtmt)
      implicit none 
      type(t_fleurinput), intent(in) :: fi
      integer                        :: num_mtmt 

      integer :: itype, ineq, l, m

      num_mtmt = 0
      DO itype = 1, fi%atoms%ntype
         DO ineq = 1, fi%atoms%neq(itype)
            DO l = 0, fi%hybinp%lcutm1(itype)
               DO M = -l, l
                  num_mtmt = num_mtmt + 1
               enddo
            enddo 
         enddo 
      enddo 
   end function calc_num_mtmts
   
END MODULE m_coulombmatrix
