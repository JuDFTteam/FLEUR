!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_ab
   USE m_juDFT

   IMPLICIT NONE

CONTAINS

   SUBROUTINE hsmt_ab(sym,atoms,noco,nococonv,ilSpin,igSpin,n,na,cell,lapw,fjgj,abCoeffs,ab_size,l_nonsph,abclo,alo1,blo1,clo1)
      ! Construct the matching coefficients of the form
      ! a_{l,m,order}^{\mu,\bm{G}}(\bm{k}) = e^{i \bm{K}\cdot\bm{\tau}_{\mu}} *
      !                                      (Y_{l}^{m}(R^{\mu}\bm{K}))^{*} *
      !                                      fjgj
      ! with \bm{K} = \bm{k + G + \pm q/2}
      !
      ! This translates to (e.g. order = 1 == acof):
      ! fjgj = \frac{4\pi}{\sqrt{\Omega}W} * (udot_{l}(R_{\alpha}) * K * j_{l})^{'}(K R_{\alpha})
      !                                    -  udot_{l}^{'}(R_{\alpha}) * j_{l})(K R_{\alpha})
      ! abCoeffs(lm, k)                    = c_ph(k,igSpin) *
      !                                      CONJG(ylm(lm, k)) *
      !                                      fjgj%fj(k,l,ilSpin,igSpin)
      ! A necessary factor i**l is omitted here and instead calculated where
      ! the coefficients are used.
      ! Also, fjgj carry information of the global and local spins \sigma_{g}
      ! and \sigma{\alpha} through a K prefactor [\pm q in K] and the u[dot]
      ! respectively. The former also appears in the complex phase factor.
      USE m_constants, ONLY : fpi_const,tpi_const
      USE m_types
      USE m_ylm
      USE m_hsmt_fjgj

      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_lapw),     INTENT(IN)    :: lapw
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_fjgj),     INTENT(IN)    :: fjgj
    !     ..
    !     .. Scalar Arguments ..
      INTEGER,          INTENT(IN)    :: ilSpin, n, na, igSpin
      LOGICAL,          INTENT(IN)    :: l_nonsph
      INTEGER,          INTENT(OUT)   :: ab_size
    !     ..
    !     .. Array Arguments ..
      COMPLEX,          INTENT(INOUT) :: abCoeffs(:,:)
    !Optional arguments if abc coef for LOs are needed
      COMPLEX, INTENT(INOUT), OPTIONAL :: abclo(:, :, :, :)
      REAL,    INTENT(IN),    OPTIONAL :: alo1(:), blo1(:), clo1(:)

      INTEGER :: np, k, l, ll1, m, lmax, nkvec, lo, lm, invsfct, lmMin, lmMax, ierr
      COMPLEX :: term
      REAL    :: bmrot(3, 3)
      COMPLEX :: c_ph(MAXVAL(lapw%nv), MERGE(2, 1, noco%l_ss.OR.ANY(noco%l_unrestrictMT) &
                                                          & .OR.ANY(noco%l_spinoffd_ldau)))
      LOGICAL :: l_apw, l_abclo

      REAL,    ALLOCATABLE :: gkrot(:, :)
      COMPLEX, ALLOCATABLE :: ylm(:, :)

      l_abclo = PRESENT(abclo)
      lmax = MERGE(atoms%lnonsph(n), atoms%lmax(n), l_nonsph)
      ab_size = lmax*(lmax+2) + 1
      ! TODO: replace APW+lo check (may actually be a broken trick) by something simpler
      ! l_apw=ALL(fjgj%gj==0.0)
      l_apw = .FALSE.

      ! We skip the initialization for speed
      ! abCoeffs=0.0

      np = sym%invtab(sym%ngopr(na))
      CALL lapw%phase_factors(igSpin, atoms%taual(:, na), nococonv%qss, c_ph(:, igSpin))
      bmrot = TRANSPOSE(MATMUL(1.0 * sym%mrot(:, :, np), cell%bmat))

      ALLOCATE(ylm((lmax+1)**2, lapw%nv(igSpin)), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error("Couldn't allocate ylm")
      ALLOCATE(gkrot(3,lapw%nv(igSpin)), stat=ierr)
      IF (ierr /= 0) CALL juDFT_error("Couldn't allocate gkrot")

      ! Generate spherical harmonics
      ! gkrot = matmul(bmrot, lapw%vk(:,:,igSpin))
      ! These two lines should eventually move to the GPU:
      CALL dgemm("N","N", 3, lapw%nv(igSpin), 3, 1.0, bmrot, 3, lapw%vk(:,:,igSpin), 3, 0.0, gkrot, 3)
      CALL ylm4_batched(lmax,gkrot,ylm)

#ifndef _OPENACC
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP& SHARED(lapw,lmax,c_ph,igSpin,abCoeffs,fjgj,abclo,cell,atoms,sym) &
      !$OMP& SHARED(l_abclo,alo1,blo1,clo1,ab_size,na,n,ilSpin,bmrot, ylm) &
      !$OMP& PRIVATE(k,l,ll1,m,lm,term,invsfct,lo,nkvec) &
      !$OMP& PRIVATE(lmMin,lmMax)
#else
      !$acc kernels present(abCoeffs) default(none)
      abCoeffs(:,:)=0.0
      !$acc end kernels
#endif

      !$acc data copyin(atoms,atoms%llo,atoms%llod,atoms%nlo,cell,cell%omtil,atoms%rmt) if (l_abclo)
      !$acc parallel loop present(fjgj,fjgj%fj,fjgj%gj,abCoeffs) vector_length(32)&
      !$acc copyin(lmax,lapw,lapw%nv,lapw%vk,lapw%kvec,bmrot,c_ph, sym, sym%invsat,l_abclo, ylm) &
      !$acc present(abclo,alo1,blo1,clo1)&
      !$acc private(k,l,lm,invsfct,lo,term,lmMin,lmMax)  default(none)
      DO k = 1,lapw%nv(igSpin)
         !$acc  loop vector private(l,lmMin,lmMax)
         DO l = 0,lmax
            lmMin = l*(l+1) + 1 - l
            lmMax = l*(l+1) + 1 + l
            abCoeffs(lmMin:lmMax, k)                = fjgj%fj(k,l,ilSpin,igSpin)*c_ph(k,igSpin) * CONJG(ylm(lmMin:lmMax, k))
            abCoeffs(ab_size+lmMin:ab_size+lmMax,k) = fjgj%gj(k,l,ilSpin,igSpin)*c_ph(k,igSpin) * CONJG(ylm(lmMin:lmMax, k))
         END DO
         !$acc end loop

         IF (l_abclo) THEN
            ! Determine also the abc coeffs for LOs
            invsfct=MERGE(1,2,sym%invsat(na).EQ.0)
            term = fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)*c_ph(k,igSpin)
            !!$acc loop vector private(lo,l,nkvec,ll1,m,lm)
            DO lo = 1,atoms%nlo(n)
               l = atoms%llo(lo,n)
               DO nkvec=1,invsfct*(2*l+1)
                  IF (lapw%kvec(nkvec,lo,na)==k) THEN !This k-vector is used in LO
                     ll1 = l*(l+1) + 1
                     DO m = -l,l
                        lm = ll1 + m
                        abclo(1,m+atoms%llod+1,nkvec,lo) = term*CONJG(ylm(lm,k))*alo1(lo)
                        abclo(2,m+atoms%llod+1,nkvec,lo) = term*CONJG(ylm(lm,k))*blo1(lo)
                        abclo(3,m+atoms%llod+1,nkvec,lo) = term*CONJG(ylm(lm,k))*clo1(lo)
                     END DO
                  END IF
               END DO
            END DO
            !!$acc end loop
         END IF
      END DO !k-loop
      !$acc end parallel loop
      !$acc end data

#ifndef _OPENACC
      !$OMP END PARALLEL DO
#endif

      IF (.NOT.l_apw) ab_size=ab_size*2
   END SUBROUTINE hsmt_ab
END MODULE m_hsmt_ab
