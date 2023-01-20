! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_sph
   USE m_juDFT
   IMPLICIT NONE

   INTERFACE hsmt_sph
#ifdef _OPENACC
      MODULE PROCEDURE hsmt_sph_acc
#else
      MODULE PROCEDURE hsmt_sph_cpu
#endif
   END INTERFACE

CONTAINS

   SUBROUTINE hsmt_sph_acc(n,atoms,fmpi,isp,input,nococonv,igSpinPr,igSpin,chi,lapw,el,e_shift,usdus,fjgj,smat,hmat,set0,l_fullj,lapwq,fjgjq)
      USE m_constants, ONLY : fpi_const, tpi_const
      USE m_types
      USE m_hsmt_fjgj
#ifdef CPP_GPU
      USE nvtx
#endif

      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_lapw),     INTENT(IN)    :: lapw
      TYPE(t_usdus),    INTENT(IN)    :: usdus
      TYPE(t_fjgj),     INTENT(IN)    :: fjgj
      CLASS(t_mat),     INTENT(INOUT) :: smat, hmat
      LOGICAL,          INTENT(IN)    :: l_fullj, set0  !if true, initialize the smat matrix with zeros

      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwq
      TYPE(t_fjgj), OPTIONAL, INTENT(IN) :: fjgjq

      ! Scalar Arguments
      INTEGER, INTENT(IN) :: n, isp, igSpinPr, igSpin
      COMPLEX, INTENT(IN) :: chi

      ! Array Arguments
      REAL,    INTENT(IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
      REAL,    INTENT(IN) :: e_shift!(atoms%ntype,input%jspins)

      ! Local Scalars
      REAL :: tnn(3), elall, fjkiln, gjkiln, ddnln, ski(3)
      REAL :: apw_lo1, apw_lo2, w1

      INTEGER :: ikG0, ikG, ikGPr, l, nn, l3, jv, kj_off, kj_vec

      ! Local Arrays
      REAL :: fleg1(0:atoms%lmaxd), fleg2(0:atoms%lmaxd), fl2p1(0:atoms%lmaxd)
      REAL :: qssAdd(3), qssAddPr(3)
      REAL :: plegend(0:2)
      REAL :: xlegend
      REAL :: VecHelpS, VecHelpH
      REAL :: cph_re, cph_im
      REAL :: dot, fct, fct2

      COMPLEX :: cfac

      LOGICAL :: l_samelapw

      TYPE(t_lapw) :: lapwPr

      CALL timestart("spherical setup")
      l_samelapw = .FALSE.
      IF (.NOT.PRESENT(lapwq)) l_samelapw = .TRUE.
      IF (.NOT.l_samelapw) THEN
         lapwPr = lapwq
      ELSE
         lapwPr = lapw
      END IF
      !call nvtxStartRange("hsmt_sph",1)
      DO l = 0,atoms%lmaxd
         fleg1(l) = REAL(l+l+1)/REAL(l+1)
         fleg2(l) = REAL(l)/REAL(l+1)
         fl2p1(l) = REAL(l+l+1)/fpi_const
      END DO ! l

      qssAdd   = MERGE(-nococonv%qss/2, +nococonv%qss/2,   igSpin.EQ.1)
      qssAddPr = MERGE(-nococonv%qss/2, +nococonv%qss/2, igSpinPr.EQ.1)
      !$acc  data &
      !$acc&   copyin(igSpin,igSpinPr,n,fleg1,fleg2,isp,fl2p1,el,e_shift,chi,qssAdd,qssAddPr,l_fullj)&
      !$acc&   copyin(lapw,lapwPr,atoms,fmpi,input,usdus)&
      !$acc&   copyin(lapw%nv,lapw%gvec,lapw%gk,lapwPr%nv,lapwPr%gvec,lapwPr%gk,lapw%bkpt,lapwPr%bkpt)&
      !$acc&   copyin(atoms%lmax,atoms%rmt,atoms%lnonsph,atoms%firstAtom,atoms%neq,atoms%taual)&
      !$acc&   copyin(fmpi%n_size,fmpi%n_rank)&
      !$acc&   copyin(input%l_useapw)&
      !$acc&   copyin(usdus%dus,usdus%uds,usdus%us,usdus%ddn,usdus%duds)&
      !$acc&   present(fjgj)&
      !$acc&   present(hmat,smat,hmat%data_c,hmat%data_r,smat%data_r,smat%data_c)

      !$acc parallel default(none)
      !$acc loop gang
      DO ikG =  fmpi%n_rank+1, lapw%nv(igSpin), fmpi%n_size
         !$acc loop  vector independent&
         !$acc &    PRIVATE(ikGPr,ikG0,ski,plegend,tnn,vechelps,vechelph,xlegend,fjkiln,gjkiln,ddnln,elall,l3,l,fct,fct2,cph_re,cph_im,cfac,dot)
         DO  ikGPr = 1, MERGE(lapwPr%nv(igSpinPr),MIN(ikG,lapwPr%nv(igSpinPr)),l_fullj)
            ikG0 = (ikG-1)/fmpi%n_size + 1
            ski = lapw%gvec(:,ikG,igSpin) + qssAdd(:) + lapw%bkpt + lapw%qphon

            ! Update overlap and l-diagonal hamiltonian matrix
            VecHelpS = 0.0
            VecHelpH = 0.0

            ! x for legendre polynomials
            xlegend = dot_product(lapwPr%gk(1:3,ikGPr,igSpinPr),lapw%gk(1:3,ikG,igSpin))

            !$acc loop seq
            DO  l = 0,atoms%lmax(n)
               fjkiln = fjgj%fj(ikG,l,isp,igSpin)
               gjkiln = fjgj%gj(ikG,l,isp,igSpin)
               IF (input%l_useapw) THEN
                  w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
                  apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 + fjkiln * usdus%us(l,n,isp)  * usdus%dus(l,n,isp) )
                  apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 + gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
               END IF ! useapw
               ddnln = usdus%ddn(l,n,isp)
               elall = el(l,n,isp)

               IF (l<=atoms%lnonsph(n).AND..NOT.l_fullj) elall=elall-e_shift!(isp)

               ! Legendre polynomials
               l3 = modulo(l, 3)

               IF (l == 0) THEN
                  plegend(0) = 1.0
               ELSE IF (l == 1) THEN
                  plegend(1) = xlegend
               ELSE
                  plegend(l3) = fleg1(l-1)*xlegend*plegend(modulo(l-1,3)) &
                            & - fleg2(l-1)*plegend(modulo(l-2,3))
               END IF ! l

               fct  = plegend(l3) * fl2p1(l) * ( fjkiln*fjgj%fj(ikGPr,l,isp,igSpinPr) &
                                             & + gjkiln*fjgj%gj(ikGPr,l,isp,igSpinPr)*ddnln )
               IF (.NOT.l_fullj) THEN
                  fct2 = plegend(l3)*fl2p1(l) * 0.5 * ( gjkiln*fjgj%fj(ikGPr,l,isp,igSpinPr) &
                                                    & + fjkiln*fjgj%gj(ikGPr,l,isp,igSpinPr) )
               ELSE
                  fct2 = plegend(l3)*fl2p1(l) * gjkiln*fjgj%fj(ikGPr,l,isp,igSpinPr)
               END IF

               VecHelpS = VecHelpS + fct
               VecHelpH = VecHelpH + fct*elall + fct2

               IF (input%l_useapw) THEN
                  VecHelpH = VecHelpH + plegend(l3) * ( apw_lo1*fjgj%fj(ikGPr,l,isp,igSpinPr) &
                                                    & + apw_lo2*fjgj%gj(l,ikGPr,isp,igSpinPr) )
               END IF ! useapw
            END DO ! l
            !$end acc
            ! Set up phase factors
            cph_re = 0.0
            cph_im = 0.0
            DO nn = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
               tnn(1:3) = tpi_const*atoms%taual(1:3,nn)

               dot = DOT_PRODUCT(ski(1:3) - lapwPr%gvec(1:3,ikGPr,igSpinPr) - qssAddPr(1:3) - lapwPr%bkpt - lapwPr%qphon, tnn(1:3))

               cph_re = cph_re + COS(dot)
               cph_im = cph_im + SIN(dot)
               ! IF (igSpinPr.NE.igSpin) cph_im=-cph_im
            END DO ! nn

            cfac = CMPLX(cph_re,cph_im)
!            ! Prefactor: i * (k + G + qssAdd - k' - G' - qssAdd')
!            IF (l_fullj) THEN
!               pref = ImagUnit * MATMUL(ski(1:3) - lapwPr%gvec(1:3,ikGPr,igSpinPr) - qssAddPr(1:3) - lapwPr%bkpt, bmat)
!               cfac = pref(idir) * cfac
!            END IF

            IF (smat%l_real) THEN
               IF (set0) THEN
                  smat%data_r(ikGPr,ikG0) = cph_re * VecHelpS
               ELSE
                  smat%data_r(ikGPr,ikG0) = &
                  smat%data_r(ikGPr,ikG0) + cph_re * VecHelpS
               END IF
               hmat%data_r(ikGPr,ikG0) = &
               hmat%data_r(ikGPr,ikG0) + cph_re * VecHelpH
            ELSE  ! real
               IF (set0) THEN
                  smat%data_c(ikGPr,ikG0) = chi*cfac * VecHelpS
               ELSE
                  smat%data_c(ikGPr,ikG0) = &
                  smat%data_c(ikGPr,ikG0) + chi*cfac * VecHelpS
               END IF
               hmat%data_c(ikGPr,ikG0) = &
               hmat%data_c(ikGPr,ikG0) + chi*cfac * VecHelpH
            END IF ! real
         END DO ! kj_off
         !$acc end loop
      END DO ! ikG

      !$acc end loop
      !$acc end parallel
      !$acc end data
      !$acc wait
      CALL timestop("spherical setup")
      !call nvtxEndRange()
      RETURN
   END SUBROUTINE hsmt_sph_acc

   SUBROUTINE hsmt_sph_cpu(n,atoms,fmpi,isp,input,nococonv,igSpinPr,igSpin,chi,lapw,el,e_shift,usdus,fjgj,smat,hmat,set0,l_fullj,lapwq, fjgjq)
      USE m_constants, ONLY : fpi_const, tpi_const
      USE m_types
      USE m_hsmt_fjgj
#ifdef CPP_GPU
      USE nvtx
#endif

      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_lapw),     INTENT(IN)    :: lapw
      TYPE(t_usdus),    INTENT(IN)    :: usdus
      TYPE(t_fjgj),     INTENT(IN)    :: fjgj
      CLASS(t_mat),     INTENT(INOUT) :: smat, hmat
      LOGICAL,          INTENT(IN)    :: l_fullj, set0  !if true, initialize the smat matrix with zeros

      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwq
      TYPE(t_fjgj), OPTIONAL, INTENT(IN) :: fjgjq

      ! Scalar Arguments
      INTEGER, INTENT(IN) :: n, isp, igSpinPr, igSpin
      COMPLEX, INTENT(IN) :: chi

      ! Array Arguments
      REAL,    INTENT(IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
      REAL,    INTENT(IN) :: e_shift!(atoms%ntype,input%jspins)

      ! Local Scalars
      REAL :: tnn(3), elall, fjkiln, gjkiln, ddnln, ski(3)
      REAL :: apw_lo1, apw_lo2, w1

      INTEGER :: ikG0, ikG, ikGPr, l, nn, kj_end, l3, jv, kj_off, kj_vec

      LOGICAL :: l_samelapw

      TYPE(t_lapw) :: lapwPr
      TYPE(t_fjgj) :: fjgjPr

      ! Local Arrays
      REAL :: fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
      REAL :: qssAdd(3),qssAddPr(3)

      REAL, ALLOCATABLE :: plegend(:,:)
      REAL, ALLOCATABLE :: xlegend(:)
      REAL, ALLOCATABLE :: VecHelpS(:),VecHelpH(:)
      REAL, ALLOCATABLE :: cph_re(:), cph_im(:)
      REAL, ALLOCATABLE :: dot(:), fct(:), fct2(:)

      COMPLEX, ALLOCATABLE :: cfac(:)

      INTEGER, PARAMETER :: NVEC = 128

      INTEGER :: NVEC_rem  !remainder

      CALL timestart("spherical setup")
      l_samelapw = .FALSE.
      IF (.NOT.PRESENT(lapwq)) l_samelapw = .TRUE.
      IF (.NOT.l_samelapw) THEN
         lapwPr = lapwq
         fjgjPr = fjgjq
      ELSE
         lapwPr = lapw
         fjgjPr = fjgj
      END IF
      !call nvtxStartRange("hsmt_sph",1)
      DO l = 0,atoms%lmaxd
         fleg1(l) = REAL(l+l+1)/REAL(l+1)
         fleg2(l) = REAL(l)/REAL(l+1)
         fl2p1(l) = REAL(l+l+1)/fpi_const
      END DO ! l

      !$OMP     PARALLEL DEFAULT(NONE)&
      !$OMP     SHARED(lapw,lapwPr,atoms,nococonv,fmpi,input,usdus,smat,hmat)&
      !$OMP     SHARED(igSpin,igSpinPr,n,fleg1,fleg2,fjgj,fjgjPr,isp,fl2p1,el,e_shift,chi,set0,l_fullj)&
      !$OMP     PRIVATE(ikG0,ikG,ski,ikGPr,kj_off,kj_vec,plegend,xlegend,l,l3,kj_end,qssAdd,qssAddPr,fct2)&
      !$OMP     PRIVATE(cph_re,cph_im,cfac,dot,nn,tnn,fjkiln,gjkiln)&
      !$OMP     PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct)&
      !$OMP     PRIVATE(VecHelpS,VecHelpH,NVEC_rem)
      ALLOCATE(cph_re(NVEC),cph_im(NVEC),cfac(NVEC))
      ALLOCATE(dot(NVEC),fct(NVEC),fct2(NVEC))
      ALLOCATE(plegend(NVEC,0:2))
      ALLOCATE(xlegend(NVEC))
      ALLOCATE(VecHelpS(NVEC),VecHelpH(NVEC))
      qssAdd   = MERGE(-nococonv%qss/2,+nococonv%qss/2,igSpin  .EQ.1)
      qssAddPr = MERGE(-nococonv%qss/2,+nococonv%qss/2,igSpinPr.EQ.1)
      !$OMP      DO SCHEDULE(DYNAMIC,1)
      DO ikG =  fmpi%n_rank+1, lapw%nv(igSpin), fmpi%n_size
         kj_end = MERGE(lapwPr%nv(igSpinPr),min(ikG,lapwPr%nv(igSpinPr)),l_fullj)
         ikG0 = (ikG-1)/fmpi%n_size + 1
         ski = lapw%gvec(:,ikG,igSpin) + qssAdd(:) + lapw%bkpt + lapw%qphon
         DO kj_off = 1, lapwPr%nv(igSpinPr), NVEC
            NVEC_rem = NVEC
            kj_vec = kj_off - 1 + NVEC
            IF (kj_vec > kj_end) THEN
               kj_vec = kj_end
               NVEC_rem = kj_end - kj_off + 1
            END IF
            IF (NVEC_rem<0 ) exit
            ! Update overlap and l-diagonal hamiltonian matrix
            VecHelpS = 0.0
            VecHelpH = 0.0

            ! x for legendre polynomials
            DO jv = 0, NVEC_rem-1
               ikGPr = jv + kj_off
               xlegend(jv+1) =dot_product(lapwPr%gk(1:3,ikGPr,igSpinPr),lapw%gk(1:3,ikG,igSpin))
            END DO ! ikGPr

            DO  l = 0,atoms%lmax(n)
               fjkiln = fjgj%fj(ikG,l,isp,igSpin)
               gjkiln = fjgj%gj(ikG,l,isp,igSpin)

               IF (input%l_useapw) THEN
                  w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
                  apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 + fjkiln * usdus%us(l,n,isp)  * usdus%dus(l,n,isp) )
                  apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 + gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
               END IF ! useapw

               IF (l_fullj) THEN
                  !apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( usdus%us(l,n,isp)  * usdus%duds(l,n,isp) * gjkiln &
                  !                                           & + usdus%us(l,n,isp)  * usdus%dus(l,n,isp)  * fjkiln )
                  !apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( usdus%uds(l,n,isp) * usdus%dus(l,n,isp)  * fjkiln &
                  !                                           & + usdus%uds(l,n,isp) * usdus%duds(l,n,isp) * gjkiln)
                  w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
                  apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 + fjkiln * usdus%us(l,n,isp)  * usdus%dus(l,n,isp) )
                  apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 + gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
               END IF

               ddnln = usdus%ddn(l,n,isp)
               elall = el(l,n,isp)

               IF (l<=atoms%lnonsph(n).AND..NOT.l_fullj) elall=elall-e_shift!(isp)

               ! Legendre polynomials
               l3 = modulo(l, 3)
               IF (l == 0) THEN
                  plegend(:NVEC_REM,0) = 1.0
               ELSE IF (l == 1) THEN
                  plegend(:NVEC_REM,1) = xlegend(:NVEC_REM)
               ELSE
                  plegend(:NVEC_REM,l3) = fleg1(l-1)*xlegend(:NVEC_REM)*plegend(:NVEC_REM,modulo(l-1,3)) &
                                      & - fleg2(l-1)*plegend(:NVEC_REM,modulo(l-2,3))
               END IF ! l

               fct(:NVEC_REM)  = plegend(:NVEC_REM,l3) * fl2p1(l) * ( fjkiln*fjgjPr%fj(kj_off:kj_vec,l,isp,igSpinPr) &
                                                                  & + gjkiln*fjgjPr%gj(kj_off:kj_vec,l,isp,igSpinPr)*ddnln )

               !IF (.NOT.l_fullj) THEN
               IF (.TRUE.) THEN
                  fct2(:NVEC_REM) = plegend(:NVEC_REM,l3) * fl2p1(l) * 0.5 * ( gjkiln*fjgjPr%fj(kj_off:kj_vec,l,isp,igSpinPr) &
                                                                           & + fjkiln*fjgjPr%gj(kj_off:kj_vec,l,isp,igSpinPr) )
               ELSE
                  fct2(:NVEC_REM) = plegend(:NVEC_REM,l3) * fl2p1(l) * gjkiln*fjgjPr%fj(kj_off:kj_vec,l,isp,igSpinPr)
               END IF

               VecHelpS(:NVEC_REM) = VecHelpS(:NVEC_REM) + fct(:NVEC_REM)
               VecHelpH(:NVEC_REM) = VecHelpH(:NVEC_REM) + fct(:NVEC_REM)*elall + fct2(:NVEC_REM)

               !IF (input%l_useapw.OR.(l_fullj.AND.l==0)) THEN
               IF (input%l_useapw.OR.l_fullj) THEN ! correction
                  VecHelpH(:NVEC_REM) = VecHelpH(:NVEC_REM) + plegend(:NVEC_REM,l3) * ( apw_lo1*fjgjPr%fj(kj_off:kj_vec,l,isp,igSpinPr) + &
                                                                                      & apw_lo2*fjgjPr%gj(kj_off:kj_vec,l,isp,igSpinPr) )
               END IF ! useapw
            END DO ! l
            ! Set up phase factors
            cph_re = 0.0
            cph_im = 0.0
            DO nn = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
               tnn(1:3) = tpi_const*atoms%taual(1:3,nn)
               DO jv = 0, NVEC_rem-1
                  ikGPr = jv + kj_off
                  dot(jv+1) = DOT_PRODUCT(ski(1:3) - lapwPr%gvec(1:3,ikGPr,igSpinPr) - qssAddPr(1:3) - lapwPr%bkpt - lapwPr%qphon, tnn(1:3))
               END DO ! ikGPr
               cph_re(:NVEC_REM) = cph_re(:NVEC_REM) + COS(dot(:NVEC_REM))
               cph_im(:NVEC_REM) = cph_im(:NVEC_REM) + SIN(dot(:NVEC_REM))
               cfac(:NVEC_REM) = CMPLX(cph_re(:NVEC_REM),cph_im(:NVEC_REM))
               ! IF (igSpinPr.NE.igSpin) cph_im=-cph_im
            END DO ! nn

            IF (smat%l_real) THEN
               IF (set0) THEN
                  smat%data_r(kj_off:kj_vec,ikG0) = cph_re(:NVEC_REM) * VecHelpS(:NVEC_REM)
               ELSE
                  smat%data_r(kj_off:kj_vec,ikG0) = &
                  smat%data_r(kj_off:kj_vec,ikG0) + cph_re(:NVEC_REM) * VecHelpS(:NVEC_REM)
               END IF
               hmat%data_r(kj_off:kj_vec,ikG0) = &
               hmat%data_r(kj_off:kj_vec,ikG0) + cph_re(:NVEC_REM) * VecHelpH(:NVEC_REM)
            ELSE  ! real
               IF (set0) THEN
                  smat%data_c(kj_off:kj_vec,ikG0) = chi*cfac(:NVEC_REM) * VecHelpS(:NVEC_REM)
               ELSE
                  smat%data_c(kj_off:kj_vec,ikG0) = &
                  smat%data_c(kj_off:kj_vec,ikG0) + chi*cfac(:NVEC_REM) * VecHelpS(:NVEC_REM)
               END IF
               hmat%data_c(kj_off:kj_vec,ikG0) = &
               hmat%data_c(kj_off:kj_vec,ikG0) + chi*cfac(:NVEC_REM) * VecHelpH(:NVEC_REM)
            END IF ! real
         END DO ! kj_off
      END DO ! ikG
      !$OMP     END DO
      DEALLOCATE(plegend)
      DEALLOCATE(VecHelpS,VecHelpH)
      !$OMP     END PARALLEL
      CALL timestop("spherical setup")
      RETURN
   END SUBROUTINE hsmt_sph_cpu

END MODULE m_hsmt_sph
