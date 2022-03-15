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
      module procedure hsmt_sph_acc
#else
      module procedure hsmt_sph_cpu
#endif
   END INTERFACE

CONTAINS

SUBROUTINE hsmt_sph_acc(n,atoms,fmpi,isp,input,nococonv,iintsp,jintsp,chi,lapw,el,e_shift,usdus,fjgj,smat,hmat,set0)
   USE m_constants, ONLY : fpi_const,tpi_const
   USE m_types
   USE m_hsmt_fjgj

   IMPLICIT NONE
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_mpi),INTENT(IN)        :: fmpi
   TYPE(t_nococonv),INTENT(IN)   :: nococonv
   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_lapw),INTENT(IN)       :: lapw
   TYPE(t_usdus),INTENT(IN)      :: usdus
   TYPE(t_fjgj),INTENT(IN)       :: fjgj
   CLASS(t_mat),INTENT(INOUT)    :: smat,hmat
   LOGICAL,INTENT(IN)            :: set0  !if true, initialize the smat matrix with zeros
    !     ..
   !     .. Scalar Arguments ..
   INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
   COMPLEX, INTENT(IN)  :: chi
   !     ..
   !     .. Array Arguments ..
   REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
   REAL,    INTENT (IN) :: e_shift!(atoms%ntype,input%jspins)

   !     ..
   !     .. Local Scalars ..
   REAL tnn(3), elall,fjkiln,gjkiln,ddnln,ski(3)
   REAL apw_lo1,apw_lo2,w1

   INTEGER kii,ki,kj,l,nn,kj_end,l3,jv,kj_off,kj_vec,l31

   !     ..
   !     .. Local Arrays ..
   REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
   REAL qssbti(3),qssbtj(3)
   REAL :: plegend(0:2)
   REAL :: xlegend
   REAL :: VecHelpS,VecHelpH
   REAL :: cph_re, cph_im
   REAL :: dot, fct, fct2

   CALL timestart("spherical setup")
   DO l = 0,atoms%lmaxd
      fleg1(l) = REAL(l+l+1)/REAL(l+1)
      fleg2(l) = REAL(l)/REAL(l+1)
      fl2p1(l) = REAL(l+l+1)/fpi_const
   END DO ! l
   qssbti=MERGE(- nococonv%qss/2,+ nococonv%qss/2,jintsp.EQ.1)
   qssbtj=MERGE(- nococonv%qss/2,+ nococonv%qss/2,iintsp.EQ.1)
   !$acc  data &
   !$acc&   copyin(jintsp,iintsp,n,fleg1,fleg2,isp,fl2p1,el,e_shift,chi,qssbti,qssbtj)&
   !$acc&   copyin(lapw,atoms,fmpi,input,usdus)&
   !$acc&   copyin(lapw%nv,lapw%gvec,lapw%gk)&
   !$acc&   copyin(atoms%lmax,atoms%rmt,atoms%lnonsph,atoms%neq,atoms%taual)&
   !$acc&   copyin(fmpi%n_size,fmpi%n_rank)&
   !$acc&   copyin(input%l_useapw)&
   !$acc&   copyin(usdus%dus,usdus%uds,usdus%us,usdus%ddn,usdus%duds)&
   !$acc&   present(fjgj)&
   !$acc&   present(hmat,smat,hmat%data_c,hmat%data_r,smat%data_r,smat%data_c)

   !$acc parallel default(none)
   !$acc loop gang
   DO  ki =  fmpi%n_rank+1, lapw%nv(jintsp), fmpi%n_size
      !$acc loop  vector independent&
      !$acc &    PRIVATE(kj, kii,ski,plegend,tnn,vechelps,vechelph,xlegend,fjkiln,gjkiln,ddnln,elall,l3,l,fct,fct2,cph_re,cph_im,dot)
      DO  kj = 1, min(ki,lapw%nv(iintsp))
         kii=(ki-1)/fmpi%n_size+1
         ski = lapw%gvec(:,ki,jintsp) + qssbti(:)

         !--->          update overlap and l-diagonal hamiltonian matrix
         VecHelpS = 0.0
         VecHelpH = 0.0

         !--->       x for legendre polynomials
         xlegend =dot_product(lapw%gk(1:3,kj,iintsp),lapw%gk(1:3,ki,jintsp))
         !$acc loop seq
         DO  l = 0,atoms%lmax(n)
            fjkiln = fjgj%fj(ki,l,isp,jintsp)
            gjkiln = fjgj%gj(ki,l,isp,jintsp)
            IF (input%l_useapw) THEN
               w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
               apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 + fjkiln * usdus%us(l,n,isp)  * usdus%dus(l,n,isp) )
               apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 + gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
            ENDIF ! useapw
            ddnln = usdus%ddn(l,n,isp)
            elall = el(l,n,isp)
            IF (l<=atoms%lnonsph(n)) elall=elall-e_shift!(isp)

            !--->       legendre polynomials
            l3 = modulo(l, 3)
            IF (l == 0) THEN
               plegend(0) = 1.0
            ELSE IF (l == 1) THEN
               plegend(1) = xlegend
            ELSE
               plegend(l3) = fleg1(l-1)*xlegend*plegend(modulo(l-1,3)) - fleg2(l-1)*plegend(modulo(l-2,3))
            END IF ! l
            fct  = plegend(l3)*fl2p1(l)       * ( fjkiln*fjgj%fj(kj,l,isp,iintsp) + gjkiln*fjgj%gj(kj,l,isp,iintsp)*ddnln )
            fct2 = plegend(l3)*fl2p1(l) * 0.5 * ( gjkiln*fjgj%fj(kj,l,isp,iintsp) + fjkiln*fjgj%gj(kj,l,isp,iintsp) )

            VecHelpS = VecHelpS + fct
            VecHelpH = VecHelpH + fct*elall + fct2

            IF (input%l_useapw) THEN
               VecHelpH = VecHelpH + plegend(l3) * ( apw_lo1*fjgj%fj(kj,l,isp,iintsp) + apw_lo2*fjgj%gj(l,kj,isp,iintsp) )
            ENDIF ! useapw

            !--->          end loop over l
         ENDDO ! l
         !$end acc
         !--->             set up phase factors
         cph_re = 0.0
         cph_im = 0.0
         DO nn = SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n))
            tnn(1:3) = tpi_const*atoms%taual(1:3,nn)

            dot = DOT_PRODUCT(ski(1:3) - lapw%gvec(1:3,kj,iintsp) - qssbtj(1:3), tnn(1:3))

            cph_re = cph_re + COS(dot)
            cph_im = cph_im - SIN(dot)
            ! IF (iintsp.NE.jintsp) cph_im=-cph_im
         END DO ! nn

         IF (smat%l_real) THEN
            IF (set0) THEN
               smat%data_r(kj,kii) = cph_re * VecHelpS
            ELSE
               smat%data_r(kj,kii) = &
               smat%data_r(kj,kii) + cph_re * VecHelpS
            ENDIF
            hmat%data_r(kj,kii) = &
            hmat%data_r(kj,kii) + cph_re * VecHelpH
         ELSE  ! real
            IF (set0) THEN
               smat%data_c(kj,kii) = chi*cmplx(cph_re,cph_im) * VecHelpS
            ELSE
               smat%data_c(kj,kii) = &
               smat%data_c(kj,kii) + chi*cmplx(cph_re,cph_im) * VecHelpS
            ENDIF
            hmat%data_c(kj,kii) = &
            hmat%data_c(kj,kii) + chi*cmplx(cph_re,cph_im) * VecHelpH
         ENDIF ! real


      END DO ! kj_off
      !$acc end loop
      !--->    end loop over ki
   ENDDO

   !$acc end loop
   !$acc end parallel
   !$acc end data
   !$acc wait
   CALL timestop("spherical setup")
   RETURN
END SUBROUTINE hsmt_sph_acc

SUBROUTINE hsmt_sph_cpu(n,atoms,fmpi,isp,input,nococonv,iintsp,jintsp,chi,lapw,el,e_shift,usdus,fjgj,smat,hmat,set0)
   USE m_constants, ONLY : fpi_const,tpi_const
   USE m_types
   USE m_hsmt_fjgj
   IMPLICIT NONE
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_mpi),INTENT(IN)        :: fmpi
   TYPE(t_nococonv),INTENT(IN)   :: nococonv
   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_lapw),INTENT(IN)       :: lapw
   TYPE(t_usdus),INTENT(IN)      :: usdus
   TYPE(t_fjgj),INTENT(IN)       :: fjgj
   CLASS(t_mat),INTENT(INOUT)    :: smat,hmat
   LOGICAL,INTENT(IN)            :: set0  !if true, initialize the smat matrix with zeros
    !     ..
   !     .. Scalar Arguments ..
   INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
   COMPLEX, INTENT(IN)  :: chi
   !     ..
   !     .. Array Arguments ..
   REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
   REAL,    INTENT (IN) :: e_shift!(atoms%ntype,input%jspins)

   !     ..
   !     .. Local Scalars ..
   REAL tnn(3), elall,fjkiln,gjkiln,ddnln,ski(3)
   REAL apw_lo1,apw_lo2,w1

   INTEGER kii,ki,kj,l,nn,kj_end,l3,jv,kj_off,kj_vec,l31

   !     ..
   !     .. Local Arrays ..
   REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
   REAL qssbti(3),qssbtj(3)
   REAL, ALLOCATABLE :: plegend(:,:)
   REAL, ALLOCATABLE :: xlegend(:)
   REAL, ALLOCATABLE :: VecHelpS(:),VecHelpH(:)
   REAL, ALLOCATABLE :: cph_re(:), cph_im(:)
   REAL, ALLOCATABLE :: dot(:), fct(:), fct2(:)
   INTEGER, PARAMETER :: NVEC = 128
   INTEGER :: NVEC_rem  !remainder

   CALL timestart("spherical setup")
   DO l = 0,atoms%lmaxd
      fleg1(l) = REAL(l+l+1)/REAL(l+1)
      fleg2(l) = REAL(l)/REAL(l+1)
      fl2p1(l) = REAL(l+l+1)/fpi_const
   END DO ! l
  !$OMP     PARALLEL DEFAULT(NONE)&
  !$OMP     SHARED(lapw,atoms,nococonv,fmpi,input,usdus,smat,hmat)&
  !$OMP     SHARED(jintsp,iintsp,n,fleg1,fleg2,fjgj,isp,fl2p1,el,e_shift,chi,set0)&
  !$OMP     PRIVATE(kii,ki,ski,kj,kj_off,kj_vec,plegend,xlegend,l,l3,kj_end,qssbti,qssbtj,fct2)&
  !$OMP     PRIVATE(cph_re,cph_im,dot,nn,tnn,fjkiln,gjkiln)&
  !$OMP     PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct)&
  !$OMP     PRIVATE(VecHelpS,VecHelpH,NVEC_rem)
   ALLOCATE(cph_re(NVEC),cph_im(NVEC))
   ALLOCATE(dot(NVEC),fct(NVEC),fct2(NVEC))
   ALLOCATE(plegend(NVEC,0:2))
   ALLOCATE(xlegend(NVEC))
   ALLOCATE(VecHelpS(NVEC),VecHelpH(NVEC))
   qssbti=MERGE(- nococonv%qss/2,+ nococonv%qss/2,jintsp.EQ.1)
   qssbtj=MERGE(- nococonv%qss/2,+ nococonv%qss/2,iintsp.EQ.1)
   !$OMP      DO SCHEDULE(DYNAMIC,1)
   DO  ki =  fmpi%n_rank+1, lapw%nv(jintsp), fmpi%n_size
      kj_end=min(ki,lapw%nv(iintsp))
      kii=(ki-1)/fmpi%n_size+1
      ski = lapw%gvec(:,ki,jintsp) + qssbti(:)
      DO  kj_off = 1, lapw%nv(iintsp), NVEC
         NVEC_rem = NVEC
         kj_vec = kj_off - 1 + NVEC
         IF (kj_vec > kj_end) THEN
            kj_vec = kj_end
            NVEC_rem = kj_end - kj_off + 1
         ENDIF
         if (NVEC_rem<0 ) exit
         !--->          update overlap and l-diagonal hamiltonian matrix
         VecHelpS = 0.0
         VecHelpH = 0.0

         !--->       x for legendre polynomials
         DO jv = 0, NVEC_rem-1
            kj = jv + kj_off
            xlegend(jv+1) =dot_product(lapw%gk(1:3,kj,iintsp),lapw%gk(1:3,ki,jintsp))
         END DO ! kj
         DO  l = 0,atoms%lmax(n)

            fjkiln = fjgj%fj(ki,l,isp,jintsp)
            gjkiln = fjgj%gj(ki,l,isp,jintsp)

            IF (input%l_useapw) THEN
               w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
               apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 + fjkiln * usdus%us(l,n,isp)  * usdus%dus(l,n,isp) )
               apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 + gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
            ENDIF ! useapw

            ddnln = usdus%ddn(l,n,isp)
            elall = el(l,n,isp)
            IF (l<=atoms%lnonsph(n)) elall=elall-e_shift!(isp)

            !--->       legendre polynomials
            l3 = modulo(l, 3)
            IF (l == 0) THEN
               plegend(:NVEC_REM,0) = 1.0
            ELSE IF (l == 1) THEN
               plegend(:NVEC_REM,1) = xlegend(:NVEC_REM)
            ELSE
               plegend(:NVEC_REM,l3) = fleg1(l-1)*xlegend(:NVEC_REM)*plegend(:NVEC_REM,modulo(l-1,3)) - fleg2(l-1)*plegend(:NVEC_REM,modulo(l-2,3))
            END IF ! l

            fct(:NVEC_REM)  = plegend(:NVEC_REM,l3)*fl2p1(l)       * ( fjkiln*fjgj%fj(kj_off:kj_vec,l,isp,iintsp) + gjkiln*fjgj%gj(kj_off:kj_vec,l,isp,iintsp)*ddnln )
            fct2(:NVEC_REM) = plegend(:NVEC_REM,l3)*fl2p1(l) * 0.5 * ( gjkiln*fjgj%fj(kj_off:kj_vec,l,isp,iintsp) + fjkiln*fjgj%gj(kj_off:kj_vec,l,isp,iintsp) )

            VecHelpS(:NVEC_REM) = VecHelpS(:NVEC_REM) + fct(:NVEC_REM)
            VecHelpH(:NVEC_REM) = VecHelpH(:NVEC_REM) + fct(:NVEC_REM)*elall + fct2(:NVEC_REM)

            IF (input%l_useapw) THEN
               VecHelpH(:NVEC_REM) = VecHelpH(:NVEC_REM) + plegend(:NVEC_REM,l3) * ( apw_lo1*fjgj%fj(kj_off:kj_vec,l,isp,iintsp) + apw_lo2*fjgj%gj(kj_off:kj_vec,l,isp,iintsp) )
            ENDIF ! useapw

            !--->          end loop over l
         ENDDO ! l
         !--->             set up phase factors
         cph_re = 0.0
         cph_im = 0.0
         DO nn = SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n))
            tnn(1:3) = tpi_const*atoms%taual(1:3,nn)
            DO jv = 0, NVEC_rem-1
               kj = jv + kj_off
               dot(jv+1) = DOT_PRODUCT(ski(1:3) - lapw%gvec(1:3,kj,iintsp) - qssbtj(1:3), tnn(1:3))
            END DO ! kj
            cph_re(:NVEC_REM) = cph_re(:NVEC_REM) + COS(dot(:NVEC_REM))
            cph_im(:NVEC_REM) = cph_im(:NVEC_REM) - SIN(dot(:NVEC_REM))
            ! IF (iintsp.NE.jintsp) cph_im=-cph_im
         END DO ! nn

         IF (smat%l_real) THEN
            IF (set0) THEN
               smat%data_r(kj_off:kj_vec,kii) = cph_re(:NVEC_REM) * VecHelpS(:NVEC_REM)
            ELSE
               smat%data_r(kj_off:kj_vec,kii) = &
               smat%data_r(kj_off:kj_vec,kii) + cph_re(:NVEC_REM) * VecHelpS(:NVEC_REM)
            ENDIF
            hmat%data_r(kj_off:kj_vec,kii) = &
            hmat%data_r(kj_off:kj_vec,kii) + cph_re(:NVEC_REM) * VecHelpH(:NVEC_REM)
         ELSE  ! real
            IF (set0) THEN
               smat%data_c(kj_off:kj_vec,kii) = chi*cmplx(cph_re(:NVEC_REM),cph_im(:NVEC_REM)) * VecHelpS(:NVEC_REM)
            ELSE
               smat%data_c(kj_off:kj_vec,kii) = &
               smat%data_c(kj_off:kj_vec,kii) + chi*cmplx(cph_re(:NVEC_REM),cph_im(:NVEC_REM)) * VecHelpS(:NVEC_REM)
            ENDIF
            hmat%data_c(kj_off:kj_vec,kii) = &
            hmat%data_c(kj_off:kj_vec,kii) + chi*cmplx(cph_re(:NVEC_REM),cph_im(:NVEC_REM)) * VecHelpH(:NVEC_REM)
         ENDIF ! real


      END DO ! kj_off
      !--->    end loop over ki
   ENDDO
   !$OMP     END DO
   DEALLOCATE(plegend)
   DEALLOCATE(VecHelpS,VecHelpH)
   !$OMP     END PARALLEL
   CALL timestop("spherical setup")
   RETURN
END SUBROUTINE hsmt_sph_cpu



END MODULE m_hsmt_sph
