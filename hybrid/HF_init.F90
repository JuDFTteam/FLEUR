
MODULE m_hf_init
   !
   !     preparations for HF and hybrid functional calculation
   !
CONTAINS
   SUBROUTINE hf_init(hybrid, atoms, input, DIMENSION, hybdat)
      USE m_types
      USE m_util
      USE m_io_hybrid
      IMPLICIT NONE
      TYPE(t_hybrid), INTENT(INOUT)     :: hybrid
      TYPE(t_atoms), INTENT(IN)         :: atoms
      TYPE(t_input), INTENT(IN)         :: input
      TYPE(t_dimension), INTENT(IN)     :: DIMENSION
      TYPE(t_hybdat), INTENT(OUT)       :: hybdat

      INTEGER:: l, m, i, l1, l2, m1, m2, ok

      !initialize hybdat%gridf for radial integration
      CALL intgrf_init(atoms%ntype, atoms%jmtd, atoms%jri, atoms%dx, atoms%rmsh, hybdat%gridf)

      !Alloc variables
      ALLOCATE (hybdat%lmaxc(atoms%ntype))
      ALLOCATE (hybdat%bas1(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype))
      ALLOCATE (hybdat%bas2(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype))
      ALLOCATE (hybdat%bas1_MT(hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype))
      ALLOCATE (hybdat%drbas1_MT(hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype))

      !sym%tau = oneD%ods%tau

      ! preparations for core states
      CALL core_init(dimension, input, atoms, hybdat%lmaxcd, hybdat%maxindxc)
      ALLOCATE (hybdat%nindxc(0:hybdat%lmaxcd, atoms%ntype), stat=ok)
      IF (ok /= 0) STOP 'eigen_hf: failure allocation hybdat%nindxc'
      ALLOCATE (hybdat%core1(atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), stat=ok)
      IF (ok /= 0) STOP 'eigen_hf: failure allocation core1'
      ALLOCATE (hybdat%core2(atoms%jmtd, hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), stat=ok)
      IF (ok /= 0) STOP 'eigen_hf: failure allocation core2'
      ALLOCATE (hybdat%eig_c(hybdat%maxindxc, 0:hybdat%lmaxcd, atoms%ntype), stat=ok)
      IF (ok /= 0) STOP 'eigen_hf: failure allocation hybdat%eig_c'
      hybdat%nindxc = 0; hybdat%core1 = 0; hybdat%core2 = 0; hybdat%eig_c = 0

      ! pre-calculate gaunt coefficients

      hybdat%maxfac = max(2*atoms%lmaxd + hybrid%maxlcutm1 + 1, 2*hybdat%lmaxcd + 2*atoms%lmaxd + 1)
      ALLOCATE (hybdat%fac(0:hybdat%maxfac), hybdat%sfac(0:hybdat%maxfac), stat=ok)
      IF (ok /= 0) STOP 'eigen_hf: failure allocation fac,hybdat%sfac'
      hybdat%fac(0) = 1
      hybdat%sfac(0) = 1
      DO i = 1, hybdat%maxfac
         hybdat%fac(i) = hybdat%fac(i - 1)*i            ! hybdat%fac(i)    = i!
         hybdat%sfac(i) = hybdat%sfac(i - 1)*sqrt(i*1.0) ! hybdat%sfac(i)   = sqrt(i!)
      END DO

      ALLOCATE (hybdat%gauntarr(2, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:hybrid%maxlcutm1, -atoms%lmaxd:atoms%lmaxd, -hybrid%maxlcutm1:hybrid%maxlcutm1), stat=ok)
      IF (ok /= 0) STOP 'eigen: failure allocation hybdat%gauntarr'
      hybdat%gauntarr = 0
      DO l2 = 0, atoms%lmaxd
         DO l1 = 0, atoms%lmaxd
            DO l = abs(l1 - l2), min(l1 + l2, hybrid%maxlcutm1)
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
