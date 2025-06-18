!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP noomp
#else
#define CPP_OMP $omp
#endif
MODULE m_hsmt_ab
   !! Module to produce matching coefficients.
   USE m_juDFT
   implicit none



CONTAINS

   SUBROUTINE hsmt_ab(sym,atoms,noco,nococonv,ilSpin,igSpin,n,na,cell,lapw,fjgj,abCoeffs,ab_size,l_nonsph,abclo,alo1,blo1,clo1)
      !! Construct the matching coefficients of the form:
      !! \begin{aligned}
      !! a_{l,m,p}^{\mu,\boldsymbol{G}}(\boldsymbol{k}) = e^{i \boldsymbol{K}\cdot\boldsymbol{\tau}_{\mu}}
      !! Y_{l}^{m*}(R^{\mu}\boldsymbol{K})
      !! f_{l,p}^{\alpha}(K)
      !! \end{aligned}
      !! with \(\boldsymbol{K} = \boldsymbol{k + G \pm q/2}\).
      !!
      !! For example, for p = 0 == acof this translates to:
      !! \begin{aligned}
      !! f_{l,0}^{\alpha}(K) = \frac{4\pi}{\sqrt{\Omega}W}
      !! (\overset{.}{u}_{l}(R_{\alpha}) * K * j_{l}^{'}(K R_{\alpha})
      !! -\overset{.}{u}_{l}^{'}(R_{\alpha}) * j_{l}(K R_{\alpha}))
      !! \end{aligned}
      !! And in the actual code into:
      !!
      !! abCoeffs(lm, k) = c_ph(k,igSpin) * CONJG(ylm(lm, k)) * fjgj%fj(k,l,ilSpin,igSpin)
      !!
      !! A factor of \(i^l\) is omitted here and instead calculated where
      !! the coefficients are used.
      !! Also, \(f_{l,p}^{\alpha}(K)\) carries information about the global and
      !! local spins \(\sigma_{g}\) and \(\sigma_{\alpha}\) through the K prefactor
      !! [\(\pm q\) in K] and the \(u/\overset{.}{u}\)
      !! respectively. The former also appears in the complex phase factor.

      USE m_constants, ONLY : fpi_const,tpi_const
      USE m_types
      USE m_ylm
      USE m_hsmt_fjgj
      USE m_matmul_dgemm

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
      REAL    :: term, bmrot(3, 3)
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
      call timestart("init")
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
      !CALL dgemm("N","N", 3, lapw%nv(igSpin), 3, 1.0, bmrot, 3, lapw%vk(:,:,igSpin), 3, 0.0, gkrot, 3)
      call blas_matmul(3,lapw%nv(igSpin),3,bmrot,lapw%vk(:,:,igspin),gkrot)
      CALL ylm4_batched(lmax,gkrot,ylm)
      call timestop("init")
      call timestart("loop")
      !CPP_OMP PARALLEL DO DEFAULT(none) &
      !CPP_OMP& SHARED(lapw,lmax,c_ph,igSpin,abCoeffs,fjgj) &
      !CPP_OMP& SHARED(ab_size,ilSpin,ylm) &
      !CPP_OMP& PRIVATE(k,l,lmMin,lmMax)
#ifdef _OPENACC
      !$acc kernels present(abCoeffs) default(none)
      abCoeffs(:,:)=0.0
      !$acc end kernels
#endif
      !$acc data copyin(c_ph,ylm)
      !$acc kernels present(fjgj,fjgj%fj,fjgj%gj,abCoeffs,c_ph, ylm)&
      !$acc copyin(lmax,ab_size,igSpin,ilSpin,lapw,lapw%nv) &
      !$acc private(k,l,lmMin,lmMax)  default(none)
      DO k = 1,lapw%nv(igSpin)
         !$acc  loop vector private(l,lmMin,lmMax)
         DO l = 0,lmax
            lmMin = l*(l+1) + 1 - l
            lmMax = l*(l+1) + 1 + l
            abCoeffs(lmMin:lmMax, k)                = fjgj%fj(k,l,ilSpin,igSpin)*c_ph(k,igSpin) * CONJG(ylm(lmMin:lmMax, k))
            abCoeffs(ab_size+lmMin:ab_size+lmMax,k) = fjgj%gj(k,l,ilSpin,igSpin)*c_ph(k,igSpin) * CONJG(ylm(lmMin:lmMax, k))
         END DO
         !$acc end loop
      enddo
      !$acc end kernels
      !CPP_OMP END PARALLEL DO

      IF (l_abclo) THEN
         ! Determine also the abc coeffs for LOs
         invsfct=MERGE(1,2,sym%invsat(na).EQ.0)
         term = fpi_const/SQRT(cell%omtil)* ((atoms%rmt(n)**2)/2)
         !$acc kernels present(abclo,alo1,blo1,clo1,c_ph,ylm) default(none)
         !$acc copyin(atoms,atoms%nlo,atoms%llo,invsfct,lapw,lapw%kvec,term,igSpin) 
         !$acc loop private(lo,l)
         DO lo = 1,atoms%nlo(n)
            l = atoms%llo(lo,n)
            !$acc loop private(nkvec,ll1)
            !CPP_OMP parallel do default(none) &
            !CPP_OMP& shared(lapw,abclo,alo1,blo1,clo1,term,c_ph,igSpin,l,invsfct) &
            !CPP_OMP& private(nkvec,ll1,m,lm)
            DO nkvec=1,invsfct*(2*l+1)
               k=lapw%kvec(nkvec,lo,na)
               ll1 = l*(l+1) + 1
               !$acc loop private(m,lm)
               DO m = -l,l
                  lm = ll1 + m
                  abclo(1,m+atoms%llod+1,nkvec,lo) = term*c_ph(k,igSpin)*CONJG(ylm(lm,k))*alo1(lo)
                  abclo(2,m+atoms%llod+1,nkvec,lo) = term*c_ph(k,igSpin)*CONJG(ylm(lm,k))*blo1(lo)
                  abclo(3,m+atoms%llod+1,nkvec,lo) = term*c_ph(k,igSpin)*CONJG(ylm(lm,k))*clo1(lo)
               END DO
               !$acc end loop
            END DO
            !CPP_OMP end parallel do
            !$acc end loop
         END DO
         !$acc end loop
         !$acc end kernels
      END IF
      !$acc end data

      call timestop("loop")


      IF (.NOT.l_apw) ab_size=ab_size*2
   END SUBROUTINE hsmt_ab
END MODULE m_hsmt_ab
