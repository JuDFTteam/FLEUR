MODULE m_hf_init
   !
   !     preparations for HF and fi%hybinp functional calculation
   !
CONTAINS
   SUBROUTINE hf_init(mpdata,fi, hybdat)
      USE m_types
      USE m_hybrid_core
      USE m_util
      use m_intgrf
      USE m_types_hybdat
      IMPLICIT NONE
      TYPE(t_mpdata), intent(inout)     :: mpdata
      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_hybdat), INTENT(INOUT)     :: hybdat

      INTEGER:: l, m, i, l1, l2, m1, m2


      !initialize hybdat%gridf for radial integration
      CALL intgrf_init(fi%atoms%ntype, fi%atoms%jmtd, fi%atoms%jri, fi%atoms%dx, fi%atoms%rmsh, hybdat%gridf)
      ! preparations for core states
      CALL core_init( fi%input, fi%atoms, hybdat%lmaxcd, hybdat%maxindxc)
      hybdat%maxfac = max(2*fi%atoms%lmaxd + maxval(fi%hybinp%lcutm1) + 1, 2*hybdat%lmaxcd + 2*fi%atoms%lmaxd + 1)

      !Alloc variables
      call hybdat%free()
      call hybdat%allocate(fi, mpdata%num_radfun_per_l)

      ! pre-calculate gaunt coefficients
      hybdat%fac(0) = 1
      hybdat%sfac(0) = 1
      DO i = 1, hybdat%maxfac
         hybdat%fac(i) = hybdat%fac(i - 1)*i            ! hybdat%fac(i)    = i!
         hybdat%sfac(i) = hybdat%sfac(i - 1)*sqrt(i*1.0) ! hybdat%sfac(i)   = sqrt(i!)
      END DO

      DO l2 = 0, fi%atoms%lmaxd
         DO l1 = 0, fi%atoms%lmaxd
            DO l = abs(l1 - l2), min(l1 + l2, maxval(fi%hybinp%lcutm1))
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
   END SUBROUTINE hf_init
END MODULE m_hf_init
