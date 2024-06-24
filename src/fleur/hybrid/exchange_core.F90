!     Calculates the HF exchange term
!
!                                          s          s*          s            s*
!                                       phi    (r) phi     (r) phi     (r') phi    (r')
!                         occ.             n_1k       n'k+q       n'k+q        n_2k
!     exchange(n,q)  =  - SUM  INT INT  ------------------------------------------- dr dr'
!                         k,n'                           | r - r' |
!
!                         occ                  s          s    ~        ~       s         s
!                    =  - SUM  SUM  v     < phi      | phi     M    > < M    phi     | phi      >
!                         k,n' I,J   k,IJ      n'k+q      n_1k  q,I      q,J    n_2k      n'k+q
!
!     for the different combinations of n_1 and n_2 and  n' runs over the core states.
!     ( n_1,n_2:  valence-valence, core-core,core-valence )
!
!     It is done directly without employing the mixed basis set.

MODULE m_exchange_core
   USE m_types_hybdat
#ifdef CPP_MPI 
   use mpi 
#endif
#ifdef _OPENACC
   USE cublas
   use openacc
#define CPP_zgemm cublaszgemm
#else
#define CPP_zgemm zgemm
#endif

CONTAINS
   SUBROUTINE exchange_vccv1(nk, fi, mpdata, hybdat, jsp, lapw, submpi,&
                             nsymop, nsest, indx_sest, a_ex, results, cmt, mat_ex)
      use m_juDFT
      USE m_types
      USE m_constants
      use m_wavefproducts_aux
      USE m_util
      use m_intgrf
      USE m_wrapper
      USE m_io_hybrid
      use m_calc_cmt
      IMPLICIT NONE
      type(t_fleurinput), intent(in) :: fi
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_mpdata), intent(in)     :: mpdata
      TYPE(t_lapw), INTENT(IN)       :: lapw
      type(t_hybmpi), intent(in)     :: submpi

      !     -scalars -
      INTEGER, INTENT(IN)      :: jsp
      INTEGER, INTENT(IN)      :: nk
      INTEGER, INTENT(IN)      :: nsymop
      REAL, INTENT(IN)         :: a_ex
      !     - arays -
      INTEGER, INTENT(IN)      ::  nsest(:), indx_sest(:,:)
      complex, intent(in)      :: cmt(:,:,:)

      TYPE(t_mat), INTENT(INOUT):: mat_ex
      !     - local scalars -
      INTEGER                 ::  iatom, itype, l, l1, l2
      INTEGER                 ::  ll, lm, m1, m2, p1, p2, n, n1, n2, nn2, i, j
      INTEGER                 ::  m

      REAL                    ::  rdum

      !     - local arrays -
      INTEGER, ALLOCATABLE     ::  larr(:), larr2(:)
      INTEGER, ALLOCATABLE     ::  parr(:), parr2(:)

      integer                 ::  nbasfcn, ierr, buf_sz, ld_integral, ld_carr, ld_tmp, ld_dotres
      REAL                    ::  integrand(fi%atoms%jmtd)
      REAL                    ::  primf1(fi%atoms%jmtd), primf2(fi%atoms%jmtd)
      REAL, ALLOCATABLE       ::  fprod(:, :), fprod2(:, :)

      COMPLEX, ALLOCATABLE    :: carr2(:, :), carr3(:, :), ctmp_vec(:)
      type(t_mat)             :: exchange
      complex, allocatable    :: integral(:,:), carr(:,:), tmp(:,:), dot_result(:,:)

      call timestart("exchange_vccv1")
      ! read in mt wavefunction coefficients from file cmt
      nbasfcn = calc_number_of_basis_functions(lapw, fi%atoms, fi%noco)
      
      call exchange%alloc(mat_ex%l_real, hybdat%nbands(nk,jsp), hybdat%nbands(nk,jsp))
      allocate(fprod(fi%atoms%jmtd, 5), stat=ierr)
      if(ierr /= 0) call judft_error("alloc fprod failed")
      
      allocate(larr(5), parr(5))

      iatom = 0
      rdum = 0

      call timestart("atom_loop")
      allocate(dot_result(hybdat%nbands(nk,jsp), hybdat%nbands(nk,jsp)), stat=ierr, source=cmplx_0)
      if(ierr /= 0) call judft_error("can't alloc dot_result")
      ld_dotres   = size(dot_result,1)

      
      !$acc data copyin(indx_sest, nsest, hybdat, hybdat%nbands, exchange) create(dot_result) copyout(exchange%data_r, exchange%data_c)
         if(exchange%l_real) then
            !$acc kernels present(exchange, exchange%data_r)
            exchange%data_r(:,:) = 0.0
            !$acc end kernels
         else 
            !$acc kernels present(exchange, exchange%data_c)
            exchange%data_c(:,:) = cmplx_0
            !$acc end kernels
         endif
         do iatom = 1+submpi%rank,fi%atoms%nat, submpi%size 
            itype = fi%atoms%itype(iatom)
            DO l1 = 0, hybdat%lmaxc(itype)
               DO p1 = 1, hybdat%nindxc(l1, itype)

                  DO l = 0, fi%hybinp%lcutm1(itype)

                     ! Define core-valence product functions
                     call timestart("Define core-valence prod.-func")
                     n = 0
                     DO l2 = 0, fi%atoms%lmax(itype)
                        IF (l < ABS(l1 - l2) .OR. l > l1 + l2) CYCLE

                        DO p2 = 1, mpdata%num_radfun_per_l(l2, itype)
                           n = n + 1
                           M = SIZE(fprod, 2)
                           IF (n > M) THEN
                              allocate(fprod2(fi%atoms%jmtd, M), larr2(M), parr2(M))
                              fprod2 = fprod; larr2 = larr; parr2 = parr
                              deallocate(fprod, larr, parr)
                              allocate(fprod(fi%atoms%jmtd, M + 5), larr(M + 5), parr(M + 5))
                              fprod(:, :M) = fprod2
                              larr(:M) = larr2
                              parr(:M) = parr2
                              deallocate(fprod2, larr2, parr2)
                           END IF
                           fprod(:fi%atoms%jri(itype), n) = (hybdat%core1(:fi%atoms%jri(itype), p1, l1, itype)*hybdat%bas1(:fi%atoms%jri(itype), p2, l2, itype) &
                                                            + hybdat%core2(:fi%atoms%jri(itype), p1, l1, itype)*hybdat%bas2(:fi%atoms%jri(itype), p2, l2, itype))/fi%atoms%rmsh(:fi%atoms%jri(itype), itype)
                           larr(n) = l2
                           parr(n) = p2
                        END DO
                     END DO
                     call timestop("Define core-valence prod.-func")

                     ! Evaluate radial integrals (special part of Coulomb matrix : contribution from single MT)

                     call timestart("Eval rad. integr")
                     allocate(integral(n,n), stat=ierr, source=cmplx_0)
                     if(ierr /= 0) call judft_error("can't allocate integral")
                     ld_integral = size(integral, 1)

                     allocate(carr(n,hybdat%nbands(nk,jsp)), stat=ierr, source=cmplx_0)
                     ld_carr     = size(carr,1)
                     if(ierr /= 0) call judft_error("can't allocate carr")
                     
                     allocate(tmp(n,hybdat%nbands(nk,jsp)), stat=ierr, source=cmplx_0)
                     if(ierr /= 0) call judft_error("can't allocate carr")
                     ld_tmp      = size(tmp,1)

                     allocate(carr2(n, lapw%nv(jsp)), carr3(n, lapw%nv(jsp)), ctmp_vec(n), stat=ierr)
                     if(ierr /= 0) call judft_error("can't allocate something i guess")

                     DO i = 1, n
                        CALL primitivef(primf1, fprod(:fi%atoms%jri(itype), i)*fi%atoms%rmsh(:fi%atoms%jri(itype), itype)**(l + 1),&
                                       fi%atoms%rmsh, fi%atoms%dx, fi%atoms%jri, fi%atoms%jmtd, itype, fi%atoms%ntype)
                        ! -itype is to enforce inward integration
                        CALL primitivef(primf2, fprod(:fi%atoms%jri(itype), i)/fi%atoms%rmsh(:fi%atoms%jri(itype), itype)**l,&
                                       fi%atoms%rmsh, fi%atoms%dx, fi%atoms%jri, fi%atoms%jmtd, -itype, fi%atoms%ntype)  
                        primf1(:fi%atoms%jri(itype)) = primf1(:fi%atoms%jri(itype))/fi%atoms%rmsh(:fi%atoms%jri(itype), itype)**l
                        primf2(:fi%atoms%jri(itype)) = primf2(:fi%atoms%jri(itype))*fi%atoms%rmsh(:fi%atoms%jri(itype), itype)**(l + 1)
                        DO j = 1, n
                           integrand = fprod(:, j)*(primf1 + primf2)
                           integral(i, j) = fpi_const/(2*l + 1)*intgrf(integrand, fi%atoms, itype, hybdat%gridf)
                        END DO
                     END DO
                     call timestop("Eval rad. integr")

                     !$acc enter data copyin(integral) create(tmp)
                     ! Add everything up
                     call timestart("Add everything up")
                     DO m1 = -l1, l1
                        DO M = -l, l
                           m2 = m1 + M

                           call timestart("set carr")
                           carr = 0
                           !$OMP PARALLEL DO default(none) collapse(2)&
                           !$OMP private(n1, i, ll, lm, l2)&
                           !$OMP shared(hybdat, n, m2, mpdata, carr, cmt, larr, itype, parr, iatom)&
                           !$OMP shared(m1, M, l, l1, nk, jsp)
                           DO n1 = 1, hybdat%nbands(nk,jsp)
                              DO i = 1, n
                                 ll = larr(i)
                                 IF (ABS(m2) > ll) CYCLE

                                 lm = SUM([((2*l2 + 1)*mpdata%num_radfun_per_l(l2, itype), l2=0, ll - 1)]) &
                                       + (m2 + ll)*mpdata%num_radfun_per_l(ll, itype) + parr(i)

                                 carr(i, n1) = cmt(n1, lm, iatom)*gaunt(l1, ll, l, m1, m2, M, hybdat%maxfac, hybdat%fac, hybdat%sfac)

                              END DO
                           enddo
                           !$OMP END PARALLEL DO
                           call timestop("set carr")
                           
                           call timestart("zgemms with cpy")
                           !$acc enter data copyin(carr)

                           call timestart("zgemms without cpy")
                           !$acc host_data use_device(integral, carr, tmp, dot_result)
                           call CPP_zgemm("N", "N", n, hybdat%nbands(nk,jsp), n, cmplx_1, integral, ld_integral, carr, ld_carr, cmplx_0, tmp, ld_tmp)
                           call CPP_zgemm("C", "N", hybdat%nbands(nk,jsp), hybdat%nbands(nk,jsp), n, cmplx_1, carr, ld_carr, tmp, ld_tmp, cmplx_0, dot_result, ld_dotres)
                           !$acc end host_data
                           !$acc wait
                           call timestop("zgemms without cpy")
                           !$acc exit data delete(carr)
                           call timestop("zgemms with cpy")

                           call timestart("add to exchange")
                           if(exchange%l_real) then
#ifdef _OPENACC
                              !$acc parallel loop default(none) private(n1, n2, nn2) &
                              !$acc present(exchange, exchange%data_r, dot_result, indx_sest, nsest, hybdat, hybdat%nbands)
#else
                              !$OMP PARALLEL DO default(none) schedule(dynamic, 10)&
                              !$OMP private(n1, n2, nn2)&
                              !$OMP shared(hybdat, nsest, indx_sest, exchange, dot_result, nk, jsp)
#endif
                              DO n1 = 1, hybdat%nbands(nk,jsp)
                                 DO n2 = 1, nsest(n1)!n1
                                    nn2 = indx_sest(n2, n1)
                                    if(nn2 <= n1) then
                                       exchange%data_r(nn2, n1) = exchange%data_r(nn2, n1) + real(dot_result(n1,nn2))
                                    endif
                                 END DO
                              END DO
#ifdef _OPENACC
                              !$acc end parallel loop
#else
                              !$OMP END PARALLEL DO
#endif
                           else
#ifdef _OPENACC
                              !$acc parallel loop default(none) private(n1, n2, nn2) &
                              !$acc present(exchange, exchange%data_c, dot_result, indx_sest, nsest, hybdat, hybdat%nbands)
#else
                              !$OMP PARALLEL DO default(none) schedule(dynamic, 10)&
                              !$OMP private(n1, n2, nn2)&
                              !$OMP shared(hybdat, nsest, indx_sest, exchange, dot_result, nk, jsp)
#endif
                              DO n1 = 1, hybdat%nbands(nk,jsp)
                                 DO n2 = 1, nsest(n1)!n1
                                    nn2 = indx_sest(n2, n1)
                                    if(nn2 <= n1) then
                                       exchange%data_c(nn2, n1) = exchange%data_c(nn2, n1) + dot_result(n1,nn2)
                                    endif
                                 END DO
                              END DO
#ifdef _OPENACC
                              !$acc end parallel loop
#else
                              !$OMP END PARALLEL DO
#endif
                           endif
                           call timestop("add to exchange")
                        END DO
                     END DO
                     !$acc exit data delete(integral, tmp)
                     call timestop("Add everything up")
                     deallocate(integral, carr, tmp)
                     deallocate(carr2, carr3, ctmp_vec)
                  END DO
               END DO
            END DO
         END DO
      !$acc end data


      deallocate(dot_result)

      call timestop("atom_loop")

      buf_sz = hybdat%nbands(nk,jsp)**2

#ifdef CPP_MPI
      call timestart("exchange reduce")
      if(exchange%l_real) then 
         if(submpi%rank == 0) then
            call MPI_REDUCE(MPI_IN_PLACE, exchange%data_r, buf_sz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, submpi%comm, ierr)
         else
            call MPI_REDUCE(exchange%data_r, MPI_IN_PLACE, buf_sz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, submpi%comm, ierr)
         endif
      else
         if(submpi%rank == 0) then
            call MPI_REDUCE(MPI_IN_PLACE, exchange%data_c, buf_sz, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, submpi%comm, ierr)
         else
            call MPI_REDUCE(exchange%data_c, MPI_IN_PLACE, buf_sz, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, submpi%comm, ierr)
         endif
      endif
      call timestop("exchange reduce")
#endif

      call timestart("calc te_hfex%core")
      if(submpi%rank == 0) then
         DO n1 = 1, hybdat%nobd(nk,jsp)
            if(exchange%l_real) then
               results%te_hfex%core = real(results%te_hfex%Core - a_ex*results%w_iks(n1, nk, jsp)*exchange%data_r(n1, n1))
            else 
               results%te_hfex%core = real(results%te_hfex%Core - a_ex*results%w_iks(n1, nk, jsp)*exchange%data_c(n1, n1))
            endif
         END DO
      endif 
#ifdef CPP_MPI
      call MPI_Bcast(results%te_hfex%core, 1, MPI_DOUBLE_PRECISION, 0, submpi%comm, ierr)
#endif
      call timestop("calc te_hfex%core")

      ! add the core-valence contribution to the exchange matrix mat_ex
      ! factor 1/nsymop is needed due to the symmetrization in symmetrizeh

      call timestart("copy to mat_ex")
      if(submpi%rank == 0) then
         IF (mat_ex%l_real) THEN
            ! mat_ex%data_r = mat_ex%data_r + exchange%data_r/nsymop
            call daxpy(buf_sz, 1.0/ nsymop, exchange%data_r, 1, mat_ex%data_r, 1)
         ELSE
            !mat_ex%data_c = mat_ex%data_c + CONJG(exchange%data_c)/nsymop
            call zlacgv(buf_sz, exchange%data_c, 1)
            call zaxpy(buf_sz, cmplx_1/nsymop, exchange%data_c, 1, mat_ex%data_c, 1)
         END IF
      endif
      call timestop("copy to mat_ex")
      call timestop("exchange_vccv1")
   END SUBROUTINE exchange_vccv1

   SUBROUTINE exchange_cccc(nk, atoms, hybdat, ncstd, sym, kpts, a_ex, results)

      USE m_types
      USE m_constants
      USE m_util
      use m_intgrf
      USE m_wrapper
      USE m_gaunt
      USE m_trafo
      USE m_io_hybrid
      use m_juDFT

      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(IN)   :: hybdat
      TYPE(t_results), INTENT(INOUT)   :: results
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms

      ! - scalars -
      INTEGER, INTENT(IN)    ::  nk, ncstd

      REAL, INTENT(IN)    ::  a_ex

      ! - arays -

      ! - local scalars -
      INTEGER               ::  itype, ieq, icst, icst1, icst2, iatom, iatom0
      INTEGER               ::  l1, l2, l, ll, llmax, it2
      INTEGER               ::  m1, m2, mm, m
      INTEGER               ::  n1, n2, n

      REAL                  ::  rdum, rdum1
      ! - local arrays -
      INTEGER               ::  point(hybdat%maxindxc, -hybdat%lmaxcd:hybdat%lmaxcd, 0:hybdat%lmaxcd, atoms%nat)
      REAL                  ::  rprod(atoms%jmtd), primf1(atoms%jmtd), primf2(atoms%jmtd), integrand(atoms%jmtd)
      COMPLEX               ::  exch(ncstd, ncstd)

      !       IF ( irank == 0 ) THEN
      !         WRITE(oUnit,'(//A)') '### core-core-core-core exchange ###'
      !         WRITE(oUnit,'(/A)') '        k-point       band    exchange'
      !       END IF

      ! set up point

      call timestart("exchange_cccc")
      icst = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybdat%lmaxc(itype)
               DO M = -l, l
                  DO n = 1, hybdat%nindxc(l, itype)
                     icst = icst + 1
                     point(n, M, l, iatom) = icst
                  END DO
               END DO
            END DO
         END DO
      END DO

      llmax = 2*hybdat%lmaxcd
      exch = 0

      !$OMP PARALLEL DO default(none) schedule(dynamic)&
      !$OMP PRIVATE(itype, l1,l2,l,ll, m1,m2,M, mm, rdum, n,n1,n2, rprod, primf1, primf2)&
      !$OMP PRIVATE(integrand, iatom0, iatom, rdum1, icst1, icst2)&
      !$OMP SHARED(atoms, hybdat, llmax, point, exch)
      DO itype = 1, atoms%ntype
         iatom0 = atoms%firstAtom(itype) - 1
         DO l1 = 0, hybdat%lmaxc(itype)  ! left core state
            DO l2 = 0, hybdat%lmaxc(itype)  ! right core state
               DO l = 0, hybdat%lmaxc(itype)   ! occupied core state

                  DO ll = ABS(l1 - l), l1 + l
                     IF (ll < ABS(l - l2) .OR. ll > l + l2) CYCLE
                     IF (MOD(l + l1 + ll, 2) /= 0) CYCLE
                     IF (MOD(l + l2 + ll, 2) /= 0) CYCLE

                     DO m1 = -l1, l1
                        m2 = m1
                        IF (ABS(m2) > l2) CYCLE
                        DO M = -l, l
                           mm = M - m1
                           IF (ABS(mm) > ll) CYCLE
                           rdum = fpi_const/(2*ll + 1)*gaunt1(l, ll, l1, M, mm, m1, llmax)*gaunt1(l, ll, l2, M, mm, m2, llmax)

                           DO n = 1, hybdat%nindxc(l, itype)
                              DO n2 = 1, hybdat%nindxc(l2, itype)
                                 rprod(:atoms%jri(itype)) = (hybdat%core1(:atoms%jri(itype), n, l, itype)*hybdat%core1(:atoms%jri(itype), n2, l2, itype) &
                                                             + hybdat%core2(:atoms%jri(itype), n, l, itype)*hybdat%core2(:atoms%jri(itype), n2, l2, itype))/atoms%rmsh(:atoms%jri(itype), itype)

                                 CALL primitivef(primf1, rprod(:)*atoms%rmsh(:, itype)**(ll + 1), atoms%rmsh, atoms%dx, atoms%jri, atoms%jmtd, itype, atoms%ntype)
                                 CALL primitivef(primf2, rprod(:atoms%jri(itype))/atoms%rmsh(:atoms%jri(itype), itype)**ll, atoms%rmsh, atoms%dx, atoms%jri, atoms%jmtd, -itype, atoms%ntype)  ! -itype is to enforce inward integration

                                 primf1(:atoms%jri(itype)) = primf1(:atoms%jri(itype))/atoms%rmsh(:atoms%jri(itype), itype)**ll
                                 primf2(:atoms%jri(itype)) = primf2(:atoms%jri(itype))*atoms%rmsh(:atoms%jri(itype), itype)**(ll + 1)

                                 DO n1 = 1, hybdat%nindxc(l1, itype)

                                    rprod(:atoms%jri(itype)) = (hybdat%core1(:atoms%jri(itype), n, l, itype)*hybdat%core1(:atoms%jri(itype), n1, l1, itype) &
                                                                + hybdat%core2(:atoms%jri(itype), n, l, itype)*hybdat%core2(:atoms%jri(itype), n1, l1, itype))/atoms%rmsh(:atoms%jri(itype), itype)

                                    integrand = rprod*(primf1 + primf2)

                                    rdum1 = rdum*intgrf(integrand, atoms, itype, hybdat%gridf)

                                    iatom = iatom0
                                    DO ieq = 1, atoms%neq(itype)
                                       iatom = iatom + 1
                                       icst1 = point(n1, m1, l1, iatom)
                                       icst2 = point(n2, m2, l2, iatom)
                                       ! no race-cond since iatoms are different between loops                                       
                                       exch(icst1, icst2) = exch(icst1, icst2) + rdum1
                                    END DO
                                 END DO  !n1

                              END DO  !n2
                           END DO  !n

                        END DO  !M
                     END DO  !m1

                  END DO  !ll

               END DO  !l
            END DO  !l2
         END DO  !l1
      END DO  !itype
      !$OMP END PARALLEL DO

      IF (sym%invs) THEN
         CALL symmetrize(exch, ncstd, ncstd, 3, atoms, hybdat%lmaxc, hybdat%lmaxcd, hybdat%nindxc, sym)
         IF (ANY(ABS(AIMAG(exch)) > 1E-6)) call judft_error('exchange_cccc: exch possesses significant imaginary part')
      ENDIF
      !       DO icst = 1,ncstd
      !         IF ( irank == 0 )
      !           WRITE(oUnit,'(    ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,F12.5)')bkpt,icst,REAL(exch(icst,icst))*(-hartree_to_ev_const)
      !       END DO

      ! add core exchange contributions to the te_hfex

      DO icst1 = 1, ncstd
         results%te_hfex%core = real(results%te_hfex%core - a_ex*kpts%wtkpt(nk)*exch(icst1, icst1))
      END DO

      call timestop("exchange_cccc")
   END SUBROUTINE exchange_cccc
END MODULE m_exchange_core
