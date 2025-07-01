!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP none
#else
#define CPP_OMP $omp
#endif

MODULE m_hsmt_sph
   USE m_juDFT
   IMPLICIT NONE
!> Module for the spherical part of the Hamiltonian and overlap matrix
!> This module contains the subroutine hsmt_sph, which calculates the spherical part of
!> the Hamiltonian and overlap matrix for the given input parameters.
   private
   public :: hsmt_sph
CONTAINS

   SUBROUTINE hsmt_sph(n,atoms,fmpi,isp,input,nococonv,igSpinPr,igSpin,chi,lapw,el,e_shift,usdus,fjgj,smat,hmat,set0,l_fullj,lapwq,fjgjq)
      USE m_constants, ONLY : fpi_const, tpi_const
      USE m_types
      USE m_hsmt_fjgj


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
      TYPE(t_fjgj) :: fjgjPr

      CALL timestart("spherical setup")
      l_samelapw = .FALSE.
      IF (.NOT.PRESENT(lapwq)) l_samelapw = .TRUE.
      IF (.NOT.l_samelapw) THEN
         lapwPr = lapwq
         fjgjPr = fjgjq 
      ELSE
         lapwPr = lapw
         fjgjPr = fjgj  !TODO do we have to transfer fjgjPr to GPU?
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
      !$acc&   present(fjgj,fjgj%fj,fjgj%gj,fjgjPr,fjgjPr%fj,fjgjPr%gj)&
      !$acc&   present(hmat,smat,hmat%data_c,hmat%data_r,smat%data_r,smat%data_c)

      !$acc parallel default(none)
      !$acc loop gang
      !CPP_OMP     PARALLEL DO SCHEDULE(dynamic,1) DEFAULT(NONE)&
      !CPP_OMP     SHARED(lapw,lapwPr,atoms,nococonv,fmpi,input,usdus,smat,hmat)&
      !CPP_OMP     SHARED(igSpin,igSpinPr,n,fleg1,fleg2,fjgj,fjgjPr,isp,fl2p1,el,e_shift,chi,set0,l_fullj,qssAdd,qssAddPr)&
      !CPP_OMP     PRIVATE(ikG0,ikG,ski,ikGPr,plegend,xlegend,l,l3,fct2)&
      !CPP_OMP     PRIVATE(cph_re,cph_im,cfac,dot,nn,tnn,fjkiln,gjkiln)&
      !CPP_OMP     PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct)&
      !CPP_OMP     PRIVATE(VecHelpS,VecHelpH)
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
               fjkiln = fjgj%fj(l,ikG,isp,igSpin)
               gjkiln = fjgj%gj(l,ikG,isp,igSpin)
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

               fct  = plegend(l3) * fl2p1(l) * ( fjkiln*fjgjPr%fj(l,ikGPr,isp,igSpinPr) &
                                             & + gjkiln*fjgjPr%gj(l,ikGPr,isp,igSpinPr)*ddnln )
               !IF (.NOT.l_fullj) THEN
               IF (.TRUE.) THEN
                  fct2 = plegend(l3)*fl2p1(l) * 0.5 * ( gjkiln*fjgjPr%fj(l,ikGPr,isp,igSpinPr) &
                                                    & + fjkiln*fjgjPr%gj(l,ikGPr,isp,igSpinPr) )
               ELSE
                  fct2 = plegend(l3)*fl2p1(l) * gjkiln*fjgjPr%fj(l,ikGPr,isp,igSpinPr)
               END IF

               VecHelpS = VecHelpS + fct
               VecHelpH = VecHelpH + fct*elall + fct2

               IF (input%l_useapw) THEN
                  VecHelpH = VecHelpH + plegend(l3) * ( apw_lo1*fjgjPr%fj(l,ikGPr,isp,igSpinPr) &
                                                    & + apw_lo2*fjgjPr%gj(l,ikGPr,isp,igSpinPr) )
               END IF ! useapw
            END DO ! l
            !$acc end loop 
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
         END DO ! ikGkjPr
         !$acc end loop
      END DO ! ikG
      !CPP_OMP  END PARALLEL DO
      !$acc end loop
      !$acc end parallel
      !$acc end data
      !$acc wait
      CALL timestop("spherical setup")
      !call nvtxEndRange()
      RETURN
   END SUBROUTINE hsmt_sph

END MODULE m_hsmt_sph
