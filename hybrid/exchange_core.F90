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

CONTAINS
   SUBROUTINE exchange_vccv1(nk, input,atoms, cell, kpts, sym, noco, nococonv, oneD,&
                             mpdata, hybinp, hybdat, jsp, lapw, &
                             nsymop, nsest, indx_sest, mpi, a_ex, results, mat_ex)
      use m_juDFT
      USE m_types
      USE m_constants
      use m_wavefproducts_aux
      USE m_util
      use m_intgrf
      USE m_wrapper
      USE m_io_hybinp
      use m_calc_cmt

      IMPLICIT NONE
      TYPE(t_input),INTENT(IN)::     input
      TYPE(t_hybdat), INTENT(IN)   :: hybdat
      TYPE(t_results), INTENT(INOUT)   :: results
      TYPE(t_mpi), INTENT(IN)   :: mpi
      TYPE(t_mpdata), intent(in)   :: mpdata
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_atoms), INTENT(IN)   :: atoms
      type(t_cell), intent(in)   :: cell
      type(t_kpts), intent(in)   :: kpts
      type(t_sym), intent(in)    :: sym
      type(t_noco), intent(in)   :: noco
      type(t_nococonv), intent(in):: nococonv
      type(t_oneD), intent(in)   :: oneD
      TYPE(t_lapw), INTENT(IN)   :: lapw

      !     -scalars -
      INTEGER, INTENT(IN)      :: jsp
      INTEGER, INTENT(IN)      :: nk
      INTEGER, INTENT(IN)      ::  nsymop
      REAL, INTENT(IN)         ::  a_ex
      !     - arays -
      INTEGER, INTENT(IN)      ::  nsest(:), indx_sest(:,:)


      TYPE(t_mat), INTENT(INOUT):: mat_ex
      !     - local scalars -
      INTEGER                 ::  iatom, ieq, itype, ic, l, l1, l2
      INTEGER                 ::  ll, lm, m1, m2, p1, p2, n, n1, n2, nn2, i, j
      INTEGER                 ::  m

      REAL                    ::  rdum
      REAL                    ::  sum_offdia

      !     - local arrays -
      INTEGER, ALLOCATABLE     ::  larr(:), larr2(:)
      INTEGER, ALLOCATABLE     ::  parr(:), parr2(:)

      integer                 :: nbasfcn
      REAL                    ::  integrand(atoms%jmtd)
      REAL                    ::  primf1(atoms%jmtd), primf2(atoms%jmtd)
      REAL, ALLOCATABLE       ::  fprod(:, :), fprod2(:, :)
      complex, ALLOCATABLE    ::  integral(:, :)

      COMPLEX                 ::  cmt(hybdat%nbands(nk), hybdat%maxlmindx, atoms%nat)
      COMPLEX                 ::  exchange(hybdat%nbands(nk), hybdat%nbands(nk))
      complex                 :: c_phase(hybdat%nbands(nk))
      COMPLEX, ALLOCATABLE    :: carr(:, :), carr2(:, :), carr3(:, :), ctmp_vec(:)
      type(t_mat)             :: zmat

      complex, external   :: zdotc

      call timestart("exchange_vccv1")
      ! read in mt wavefunction coefficients from file cmt
      nbasfcn = calc_number_of_basis_functions(lapw, atoms, noco)
      CALL zmat%init(sym%invs, nbasfcn, input%neig)
      if(nk /= kpts%bkp(nk)) call juDFT_error("We should be reading the parent z-mat here!")
      call read_z(atoms, cell, hybdat, kpts, sym, noco, nococonv,  input, kpts%bkp(nk), jsp, zmat, c_phase=c_phase)
      ! zmat%matsize2 = hybdat%nbands(nk)
      call calc_cmt(atoms, cell, input, noco,nococonv, hybinp, hybdat, mpdata, kpts, &
                          sym, oneD, zmat, jsp, nk, c_phase, cmt)
      call zmat%free()

      allocate(fprod(atoms%jmtd, 5), larr(5), parr(5))

      exchange = 0
      iatom = 0
      rdum = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l1 = 0, hybdat%lmaxc(itype)
               DO p1 = 1, hybdat%nindxc(l1, itype)

                  DO l = 0, hybinp%lcutm1(itype)

                     ! Define core-valence product functions
                     call timestart("Define core-valence prod.-func")
                     n = 0
                     DO l2 = 0, atoms%lmax(itype)
                        IF (l < ABS(l1 - l2) .OR. l > l1 + l2) CYCLE

                        DO p2 = 1, mpdata%num_radfun_per_l(l2, itype)
                           n = n + 1
                           M = SIZE(fprod, 2)
                           IF (n > M) THEN
                              allocate(fprod2(atoms%jmtd, M), larr2(M), parr2(M))
                              fprod2 = fprod; larr2 = larr; parr2 = parr
                              deallocate(fprod, larr, parr)
                              allocate(fprod(atoms%jmtd, M + 5), larr(M + 5), parr(M + 5))
                              fprod(:, :M) = fprod2
                              larr(:M) = larr2
                              parr(:M) = parr2
                              deallocate(fprod2, larr2, parr2)
                           END IF
                           fprod(:atoms%jri(itype), n) = (hybdat%core1(:atoms%jri(itype), p1, l1, itype)*hybdat%bas1(:atoms%jri(itype), p2, l2, itype) &
                                                          + hybdat%core2(:atoms%jri(itype), p1, l1, itype)*hybdat%bas2(:atoms%jri(itype), p2, l2, itype))/atoms%rmsh(:atoms%jri(itype), itype)
                           larr(n) = l2
                           parr(n) = p2
                        END DO
                     END DO
                     call timestop("Define core-valence prod.-func")

                     ! Evaluate radial integrals (special part of Coulomb matrix : contribution from single MT)

                     call timestart("Eval rad. integr")
                     allocate(integral(n, n), carr(n, hybdat%nbands(nk)), carr2(n, lapw%nv(jsp)), carr3(n, lapw%nv(jsp)), ctmp_vec(n))

                     DO i = 1, n
                        CALL primitivef(primf1, fprod(:atoms%jri(itype), i)*atoms%rmsh(:atoms%jri(itype), itype)**(l + 1), atoms%rmsh, atoms%dx, atoms%jri, atoms%jmtd, itype, atoms%ntype)
                        CALL primitivef(primf2, fprod(:atoms%jri(itype), i)/atoms%rmsh(:atoms%jri(itype), itype)**l, atoms%rmsh, atoms%dx, atoms%jri, atoms%jmtd, -itype, atoms%ntype)  ! -itype is to enforce inward integration

                        primf1(:atoms%jri(itype)) = primf1(:atoms%jri(itype))/atoms%rmsh(:atoms%jri(itype), itype)**l
                        primf2(:atoms%jri(itype)) = primf2(:atoms%jri(itype))*atoms%rmsh(:atoms%jri(itype), itype)**(l + 1)
                        DO j = 1, n
                           integrand = fprod(:, j)*(primf1 + primf2)
                           integral(i, j) = fpi_const/(2*l + 1)*intgrf(integrand, atoms, itype, hybdat%gridf)
                        END DO
                     END DO
                     call timestop("Eval rad. integr")

                     ! Add everything up
                     call timestart("Add everything up")
                     DO m1 = -l1, l1
                        DO M = -l, l
                           m2 = m1 + M

                           carr = 0
                           DO n1 = 1, hybdat%nbands(nk)

                              DO i = 1, n
                                 ll = larr(i)
                                 IF (ABS(m2) > ll) CYCLE

                                 lm = SUM([((2*l2 + 1)*mpdata%num_radfun_per_l(l2, itype), l2=0, ll - 1)]) &
                                      + (m2 + ll)*mpdata%num_radfun_per_l(ll, itype) + parr(i)

                                 carr(i, n1) = cmt(n1, lm, iatom)*gaunt(l1, ll, l, m1, m2, M, hybdat%maxfac, hybdat%fac, hybdat%sfac)

                              END DO
                              DO n2 = 1, nsest(n1)!n1
                                 nn2 = indx_sest(n2, n1)
                                 call zgemv("N", n, n, cmplx_1, integral, n, carr(1,nn2), 1, cmplx_0, ctmp_vec, 1)
                                 exchange(nn2, n1) = exchange(nn2, n1) + zdotc(n, carr(1,n1), 1, ctmp_vec, 1)
                              END DO
                           END DO
                        END DO
                     END DO
                     call timestop("Add everything up")

                     deallocate(integral, carr, carr2, carr3, ctmp_vec)

                  END DO
               END DO
            END DO
         END DO
      END DO

      IF (mat_ex%l_real) THEN
         IF (ANY(ABS(AIMAG(exchange)) > 1e-10)) THEN
            IF (mpi%irank == 0) WRITE (oUnit, '(A)') 'exchangeCore: Warning! Unusually large imaginary component.'
            WRITE (*, *) MAXVAL(ABS(AIMAG(exchange)))
            call judft_error('exchangeCore: Unusually large imaginary component.')
         END IF
      ENDIF

      DO n1 = 1, hybdat%nobd(nk,jsp)
         results%te_hfex%core = real(results%te_hfex%Core - a_ex*results%w_iks(n1, nk, jsp)*exchange(n1, n1))
      END DO

      ! add the core-valence contribution to the exchange matrix mat_ex
      ! factor 1/nsymop is needed due to the symmetrization in symmetrizeh

      ic = 0
      sum_offdia = 0
      IF (mat_ex%l_real) THEN
         mat_ex%data_r = real(mat_ex%data_r + exchange/nsymop)
      ELSE
         mat_ex%data_c = mat_ex%data_c + CONJG(exchange)/nsymop
      END IF
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
      USE m_io_hybinp
      use m_juDFT
      use omp_lib

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
         iatom0 = sum([(atoms%neq(it2), it2=1,itype-1)])
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
         CALL symmetrize(exch, ncstd, ncstd, 3, .FALSE., atoms, hybdat%lmaxc, hybdat%lmaxcd, hybdat%nindxc, sym)
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
