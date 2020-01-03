
MODULE m_hf_init
   !
   !     preparations for HF and hybrid functional calculation
   !
CONTAINS
   SUBROUTINE hf_init(mpdata, hybrid, atoms, input,  hybdat)
      USE m_types
      USE m_hybrid_core
      USE m_util
      use m_intgrf
      USE m_io_hybrid
      USE m_types_hybdat
      IMPLICIT NONE
      TYPE(t_mpdata), intent(inout) :: mpdata
      TYPE(t_hybrid), INTENT(INOUT)     :: hybrid
      TYPE(t_atoms), INTENT(IN)         :: atoms
      TYPE(t_input), INTENT(IN)         :: input
      TYPE(t_hybdat), INTENT(OUT)       :: hybdat

      INTEGER:: l, m, i, l1, l2, m1, m2, ok

      !initialize hybdat%gridf for radial integration
      CALL intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, hybdat%gridf)

      !Alloc variables
      allocate(hybdat%lmaxc(atoms%ntype), source=0)
      allocate(hybdat%bas1(atoms%jmtd, maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
      allocate(hybdat%bas2(atoms%jmtd, maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
      allocate(hybdat%bas1_MT(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)
      allocate(hybdat%drbas1_MT(maxval(mpdata%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype), source=0.0)

      ! preparations for core states
      CALL core_init( input, atoms, hybdat%lmaxcd, hybdat%maxindxc)
      allocate(hybdat%nindxc(0:hybdat%lmaxcd, atoms%ntype), stat=ok, source=0)
      IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%nindxc')
      allocate(hybdat%core1(atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), stat=ok, source=0.0)
      IF (ok /= 0) call judft_error('eigen_hf: failure allocation core1')
      allocate(hybdat%core2(atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), stat=ok, source=0.0)
      IF (ok /= 0) call judft_error('eigen_hf: failure allocation core2')
      allocate(hybdat%eig_c(hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), stat=ok, source=0.0)
      IF (ok /= 0) call judft_error('eigen_hf: failure allocation hybdat%eig_c')

      ! pre-calculate gaunt coefficients

      hybdat%maxfac = max(2*atoms%lmaxd + maxval(hybrid%lcutm1) + 1, 2*hybdat%lmaxcd + 2*atoms%lmaxd + 1)
      allocate(hybdat%fac(0:hybdat%maxfac), hybdat%sfac(0:hybdat%maxfac), stat=ok, source=0.0)
      IF (ok /= 0) call judft_error('eigen_hf: failure allocation fac,hybdat%sfac')
      hybdat%fac(0) = 1
      hybdat%sfac(0) = 1
      DO i = 1, hybdat%maxfac
         hybdat%fac(i) = hybdat%fac(i - 1)*i            ! hybdat%fac(i)    = i!
         hybdat%sfac(i) = hybdat%sfac(i - 1)*sqrt(i*1.0) ! hybdat%sfac(i)   = sqrt(i!)
      END DO

      ALLOCATE(hybdat%gauntarr(2, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:maxval(hybrid%lcutm1),&
                           -atoms%lmaxd:atoms%lmaxd, -maxval(hybrid%lcutm1):maxval(hybrid%lcutm1)),&
                            stat=ok, source=0.0)
      IF (ok /= 0) call judft_error('eigen: failure allocation hybdat%gauntarr')

      DO l2 = 0, atoms%lmaxd
         DO l1 = 0, atoms%lmaxd
            DO l = abs(l1 - l2), min(l1 + l2, maxval(hybrid%lcutm1))
               DO m = -l, l
                  DO m1 = -l1, l1
                     m2 = m1 + m ! Gaunt condition -m1+m2-m = 0
                     IF (abs(m2) <= l2) hybdat%gauntarr(1, l1, l2, l, m1, m) = gaunt(l1, l2, l, m1, m2, m, hybdat%maxfac, hybdat%fac, hybdat%sfac)
                     m2 = m1 - m ! switch role of l2-index
                     IF (abs(m2) <= l2) hybdat%gauntarr(2, l1, l2, l, m1, m) = gaunt(l2, l1, l, m2, m1, m, hybdat%maxfac, hybdat%fac, hybdat%sfac)
                  END DO
               END DO
            END DO
         END DO
      END DO

      !skip_kpt = .false.

   END SUBROUTINE hf_init
END MODULE m_hf_init
