!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_sph
   USE m_juDFT
   IMPLICIT NONE

   INTERFACE hsmt_sph
      module procedure hsmt_sph_cpu
#ifdef CPP_GPU
      module procedure hsmt_sph_gpu
#endif
   END INTERFACE

CONTAINS

SUBROUTINE hsmt_sph_cpu(n,atoms,mpi,isp,input,noco,iintsp,jintsp,chi,lapw,el,e_shift,usdus,fj,gj,smat,hmat)
   USE m_constants, ONLY : fpi_const,tpi_const
   USE m_types
   IMPLICIT NONE
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_mpi),INTENT(IN)        :: mpi
   TYPE(t_noco),INTENT(IN)       :: noco
   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_lapw),INTENT(IN)       :: lapw
   TYPE(t_usdus),INTENT(IN)      :: usdus
   CLASS(t_mat),INTENT(INOUT)    :: smat,hmat
   !     ..
   !     .. Scalar Arguments ..
   INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
   COMPLEX, INTENT(IN)  :: chi
   !     ..
   !     .. Array Arguments ..
   REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
   REAL,    INTENT (IN) :: e_shift!(atoms%ntype,input%jspins)
   REAL,    INTENT (IN) :: fj(:,0:,:),gj(:,0:,:)
   !     ..
   !     .. Local Scalars ..
   REAL tnn(3), elall,fjkiln,gjkiln,ddnln,ski(3)
   REAL apw_lo1,apw_lo2,w1

   INTEGER kii,ki,kj,l,nn,kj_end,l3,jv,kj_off,kj_vec

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
!$OMP     SHARED(lapw,atoms,noco,mpi,input,usdus,smat,hmat)&
!$OMP     SHARED(jintsp,iintsp,n,fleg1,fleg2,fj,gj,isp,fl2p1,el,e_shift,chi)&
!$OMP     PRIVATE(kii,ki,ski,kj,kj_off,kj_vec,plegend,xlegend,l,l3,kj_end,qssbti,qssbtj,fct2)&
!$OMP     PRIVATE(cph_re,cph_im,dot,nn,tnn,fjkiln,gjkiln)&
!$OMP     PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct)&
!$OMP     PRIVATE(VecHelpS,VecHelpH,NVEC_rem)
   ALLOCATE(cph_re(NVEC),cph_im(NVEC))
   ALLOCATE(dot(NVEC),fct(NVEC),fct2(NVEC))
   ALLOCATE(plegend(NVEC,0:2))
   ALLOCATE(xlegend(NVEC))
   ALLOCATE(VecHelpS(NVEC),VecHelpH(NVEC))
   qssbti=MERGE(- noco%qss/2,+ noco%qss/2,jintsp.EQ.1)
   qssbtj=MERGE(- noco%qss/2,+ noco%qss/2,iintsp.EQ.1)
!$OMP      DO SCHEDULE(DYNAMIC,1)
   DO  ki =  mpi%n_rank+1, lapw%nv(jintsp), mpi%n_size
      kii=(ki-1)/mpi%n_size+1
      ski(1:3) = lapw%gvec(1:3,ki,jintsp) + qssbti(1:3)

      IF (smat%l_real) THEN
         kj_end = ki
      ELSE
         kj_end = MIN(ki,lapw%nv(iintsp))
      ENDIF
      NVEC_rem = NVEC
      DO  kj_off = 1, kj_end, NVEC
         kj_vec = kj_off - 1 + NVEC
         IF (kj_vec > kj_end) THEN
             kj_vec = kj_end
             NVEC_rem = kj_end - kj_off + 1
         ENDIF

         !--->          update overlap and l-diagonal hamiltonian matrix
         VecHelpS = 0.0
         VecHelpH = 0.0
         
         !--->       x for legendre polynomials
         DO jv = 0, NVEC_rem-1
            kj = jv + kj_off
            xlegend(jv+1) = DOT_PRODUCT(lapw%gk(1:3,kj,iintsp), lapw%gk(1:3,ki,jintsp))
         END DO ! kj
         
         DO  l = 0,atoms%lmax(n)
         
            fjkiln = fj(ki,l,jintsp)
            gjkiln = gj(ki,l,jintsp)

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
               plegend(:,0) = 1.0
            ELSE IF (l == 1) THEN
               plegend(:NVEC_REM,1) = xlegend(:NVEC_REM)
            ELSE
               plegend(:NVEC_REM,l3) = fleg1(l-1)*xlegend(:NVEC_REM)*plegend(:NVEC_REM,modulo(l-1,3)) - fleg2(l-1)*plegend(:NVEC_REM,modulo(l-2,3))
            END IF ! l

            fct(:NVEC_REM)  = plegend(:NVEC_REM,l3)*fl2p1(l)       * ( fjkiln*fj(kj_off:kj_vec,l,iintsp) + gjkiln*gj(kj_off:kj_vec,l,iintsp)*ddnln )
            fct2(:NVEC_REM) = plegend(:NVEC_REM,l3)*fl2p1(l) * 0.5 * ( gjkiln*fj(kj_off:kj_vec,l,iintsp) + fjkiln*gj(kj_off:kj_vec,l,iintsp) )

            VecHelpS(:NVEC_REM) = VecHelpS(:NVEC_REM) + fct(:NVEC_REM)
            VecHelpH(:NVEC_REM) = VecHelpH(:NVEC_REM) + fct(:NVEC_REM)*elall + fct2(:NVEC_REM)

            IF (input%l_useapw) THEN
               VecHelpH(:NVEC_REM) = VecHelpH(:NVEC_REM) + plegend(:NVEC_REM,l3) * ( apw_lo1*fj(kj_off:kj_vec,l,iintsp) + apw_lo2*gj(kj_off:kj_vec,l,iintsp) )
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
            smat%data_r(kj_off:kj_vec,kii) = &
            smat%data_r(kj_off:kj_vec,kii) + cph_re * VecHelpS(:NVEC_REM)
            hmat%data_r(kj_off:kj_vec,kii) = &
            hmat%data_r(kj_off:kj_vec,kii) + cph_re * VecHelpH(:NVEC_REM)
         ELSE  ! real
            smat%data_c(kj_off:kj_vec,kii) = &
            smat%data_c(kj_off:kj_vec,kii) + chi*cmplx(cph_re,cph_im) * VecHelpS(:NVEC_REM)
            hmat%data_c(kj_off:kj_vec,kii) = &
            hmat%data_c(kj_off:kj_vec,kii) + chi*cmplx(cph_re,cph_im) * VecHelpH(:NVEC_REM)
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


#ifdef CPP_GPU
   ATTRIBUTES(global)&
      SUBROUTINE HsmtSphGpuKernel_real(loop_size,iintsp,jintsp,nv,lmaxd,lmax,ki_start,ki_end,ki_step,nn_start,nn_end,&
                                       lnonsph,qssbti,qssbtj,gvec,gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,taual,ddn,el,e_shift,&
                                       smat_data,hmat_data,&
                                       uds,dus,us,duds,rmt)
   INTEGER, VALUE, INTENT(IN) :: loop_size
   INTEGER,VALUE,INTENT(IN) :: iintsp,jintsp,lmaxd,lmax,ki_start,ki_end,ki_step,nn_start,nn_end,lnonsph
   REAL,INTENT(IN)    :: qssbti(3),qssbtj(3)
   INTEGER,INTENT(IN) :: gvec(:,:,:),nv(2)
   REAL,INTENT(IN)    :: gk(:,:,:)
   REAL,INTENT(IN)    :: fleg1(0:lmaxd),fleg2(0:lmaxd),fl2p1(0:lmaxd)
   REAL,INTENT(IN)    :: fl2p1bt(0:lmaxd)
   REAL,INTENT(IN)    :: fj(:,0:,:),gj(:,0:,:)
   REAL,INTENT(IN)    :: taual(:,:)
   REAL,INTENT(IN)    :: ddn(0:lmaxd)
   REAL,INTENT(IN)    :: el(0:lmaxd)
   REAL,VALUE,INTENT(IN)    :: e_shift
   REAL,INTENT(INOUT) :: smat_data(:,:),hmat_data(:,:)
   !+APW
   REAL,INTENT(IN),OPTIONAL   :: uds(0:lmaxd),dus(0:lmaxd),us(0:lmaxd),duds(0:lmaxd)
   REAL,INTENT(IN),OPTIONAL   :: rmt
   !-APW

   REAL,   PARAMETER :: tpi_const=2.*3.1415926535897932
   REAL, ALLOCATABLE :: plegend(:)
   REAL cph
   REAL tnn(3), elall,fct,fct2,fjkiln,gjkiln,ddnln,ski(3)
   REAL apw_lo1,apw_lo2,apw1,w1
   INTEGER kii,ki,kj,l,nn,k
   INTEGER :: loop_start, loop_end, i

   ALLOCATE(plegend(0:lmaxd))
   plegend=0.0
   plegend(0)=1.0

   k = (blockidx%x-1)*blockdim%x + threadidx%x

   !TODO!!!
   !for seq, i.e. ki_start = 1, ki_step = 1
   loop_start = (k-1) * loop_size + 1
   loop_end = loop_start + loop_size - 1
   if (loop_end > ki_end ) loop_end = ki_end

   DO ki = loop_start,loop_end,ki_step
      !DO  ki =  ki_start,ki_end,ki_step
      DO kj = 1,ki
         kii=(ki-1)/ki_step+1
         ski = gvec(:,ki,jintsp) + qssbti
         !--->             set up phase factors
         cph = 0.0
         DO nn = nn_start,nn_end
            tnn = tpi_const*taual(:,nn)
            cph = cph + COS(DOT_PRODUCT(ski-gvec(:,kj,iintsp)-qssbtj,tnn))
            ! IF (iintsp.NE.jintsp) cph(kj)=CONJG(cph(kj))
         ENDDO
         !--->       legendre polynomials
         plegend(1) = DOT_PRODUCT(gk(:,kj,iintsp),gk(:,ki,jintsp))
         DO l = 1,lmax - 1
            plegend(l+1) = fleg1(l)*plegend(1)*plegend(l) - fleg2(l)*plegend(l-1)
         END DO
         DO  l = 0,lmax
            fjkiln = fj(ki,l,jintsp)
            gjkiln = gj(ki,l,jintsp)
            !+APW
            IF (PRESENT(uds)) THEN
               w1 = 0.5 * ( uds(l)*dus(l) + us(l)*duds(l) )
               apw_lo1 = fl2p1(l) * 0.5 * rmt**2 * ( gjkiln * w1 +&
                                                    fjkiln * us(l) * dus(l) )
               apw_lo2 = fl2p1(l) * 0.5 * rmt**2 * ( fjkiln * w1 +&
                                                    gjkiln * uds(l) * duds(l) )
            ENDIF
            !-APW
            ddnln =  ddn(l)
            elall = el(l)
            IF (l<=lnonsph) elall=elall-e_shift!(isp)
            !DO kj = 1,ki
            fct  = plegend(l)*fl2p1(l)*&
                   ( fjkiln*fj(kj,l,iintsp) + gjkiln*gj(kj,l,iintsp)*ddnln )
            fct2 = plegend(l)*fl2p1bt(l) * ( fjkiln*gj(kj,l,iintsp) + gjkiln*fj(kj,l,iintsp) )

            smat_data(kj,kii)=smat_data(kj,kii)+cph*fct
            hmat_data(kj,kii)=hmat_data(kj,kii) + cph * ( fct * elall + fct2)
            !+APW
            IF (PRESENT(uds)) THEN
               apw1 = cph * plegend(l)  * &
                      ( apw_lo1 * fj(kj,l,iintsp) + apw_lo2 * gj(kj,l,iintsp) )
               hmat_data(kj,kii)=hmat_data(kj,kii) + apw1
            ENDIF
            !-APW
            !ENDDO
            !--->          end loop over l
         ENDDO
      ENDDO
      !--->    end loop over ki
   ENDDO
   DEALLOCATE(plegend)
END SUBROUTINE HsmtSphGpuKernel_real

ATTRIBUTES(global)&
   SUBROUTINE HsmtSphGpuKernel_cmplx(loop_size,iintsp,jintsp,nv,lmaxd,lmax,ki_start,ki_end,ki_step,nn_start,nn_end,&
                                     lnonsph,chi,qssbti,qssbtj,gvec,gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,taual,ddn,el,e_shift,&
                                     smat_data,hmat_data,&
                                     uds,dus,us,duds,rmt)
INTEGER, VALUE, INTENT(IN) :: loop_size
INTEGER, VALUE, INTENT(IN) :: iintsp,jintsp,lmaxd,lmax,ki_start,ki_end,ki_step,nn_start,nn_end,lnonsph
COMPLEX, VALUE, INTENT(IN)  :: chi
REAL,INTENT(IN)    :: qssbti(3),qssbtj(3)
INTEGER,INTENT(IN) :: gvec(:,:,:),nv(2)
REAL,INTENT(IN)    :: gk(:,:,:)
REAL,INTENT(IN)    :: fleg1(0:lmaxd),fleg2(0:lmaxd),fl2p1(0:lmaxd)
REAL,INTENT(IN)    :: fl2p1bt(0:lmaxd)
REAL,INTENT(IN)    :: fj(:,0:,:),gj(:,0:,:)
REAL,INTENT(IN)    :: taual(:,:)
REAL,INTENT(IN)    :: ddn(0:lmaxd)
REAL,INTENT(IN)    :: el(0:lmaxd)
REAL, VALUE,INTENT(IN)    :: e_shift
COMPLEX,INTENT(INOUT) :: smat_data(:,:),hmat_data(:,:)
!+APW
REAL,INTENT(IN),OPTIONAL   :: uds(0:lmaxd),dus(0:lmaxd),us(0:lmaxd),duds(0:lmaxd)
REAL,INTENT(IN),OPTIONAL   :: rmt
!-APW

REAL,   PARAMETER :: tpi_const=2.*3.1415926535897932
REAL, ALLOCATABLE :: plegend(:)
COMPLEX :: cph
REAL apw_lo1,apw_lo2,w1
COMPLEX capw1
REAL tnn(3), elall,fct,fct2,fjkiln,gjkiln,ddnln,ski(3)
INTEGER kii,ki,kj,kjj,l,nn,kj_end,k
INTEGER :: loop_start, loop_end, i

ALLOCATE(plegend(0:lmaxd))
plegend=0.0
plegend(0)=1.0

k = (blockidx%x-1)*blockdim%x + threadidx%x

!TODO!!!
!for seq, i.e. ki_start = 1, ki_step = 1
loop_start = (k-1) * loop_size + 1
loop_end = loop_start + loop_size - 1
if (loop_end > ki_end ) loop_end = ki_end

DO ki = loop_start,loop_end,ki_step
   !DO  ki =  ki_start,ki_end,ki_step
   kj_end = MIN(ki,nv(iintsp))
   DO kj = 1,kj_end
      kii=(ki-1)/ki_step+1
      ski = gvec(:,ki,jintsp) + qssbti

      !--->             set up phase factors
      cph = 0.0
      DO nn = nn_start,nn_end
         tnn = tpi_const*taual(:,nn)
         cph = cph +&
               CMPLX(COS(DOT_PRODUCT(ski-gvec(:,kj,iintsp)-qssbtj,tnn)),&
                     SIN(DOT_PRODUCT(gvec(:,kj,iintsp)+qssbtj-ski,tnn)))
         ! IF (iintsp.NE.jintsp) cph(kj)=CONJG(cph(kj))
      ENDDO
      !--->       legendre polynomials
      plegend(1) = DOT_PRODUCT(gk(:,kj,iintsp),gk(:,ki,jintsp))
      DO l = 1,lmax - 1
         plegend(l+1) = fleg1(l)*plegend(1)*plegend(l) - fleg2(l)*plegend(l-1)
      END DO
      DO  l = 0,lmax
         fjkiln = fj(ki,l,jintsp)
         gjkiln = gj(ki,l,jintsp)
         !+APW
         IF (PRESENT(uds)) THEN
            w1 = 0.5 * ( uds(l)*dus(l) + us(l)*duds(l) )
            apw_lo1 = fl2p1(l) * 0.5 * rmt**2 * ( gjkiln * w1 +&
                                                 fjkiln * us(l) * dus(l) )
            apw_lo2 = fl2p1(l) * 0.5 * rmt**2 * ( fjkiln * w1 +&
                                                 gjkiln * uds(l) * duds(l) )
         ENDIF
         !-APW
         ddnln =  ddn(l)
         elall = el(l)
         IF (l<=lnonsph) elall=elall-e_shift!(isp)
         fct  = plegend(l)*fl2p1(l)*&
                ( fjkiln*fj(kj,l,iintsp) + gjkiln*gj(kj,l,iintsp)*ddnln )
         fct2 = plegend(l)*fl2p1bt(l) * ( fjkiln*gj(kj,l,iintsp) + gjkiln*fj(kj,l,iintsp) )

         smat_data(kj,kii)=smat_data(kj,kii) + chi * cph * fct
         hmat_data(kj,kii)=hmat_data(kj,kii) + chi * cph * ( fct * elall + fct2)
         !+APW
         IF (PRESENT(uds)) THEN
            capw1 = cph*plegend(l)&
                    * ( apw_lo1 * fj(kj,l,iintsp) + apw_lo2 * gj(kj,l,iintsp) )
            hmat_data(kj,kii)=hmat_data(kj,kii) + capw1
         ENDIF
         !-APW
         !--->          end loop over l
      ENDDO
   ENDDO
   !--->    end loop over ki
ENDDO
DEALLOCATE(plegend)
END SUBROUTINE HsmtSphGpuKernel_cmplx

SUBROUTINE hsmt_sph_gpu(n,atoms,mpi,isp,input,noco,iintsp,jintsp,chi,lapw,el,e_shift,usdus,fj,gj,smat,hmat)
   USE m_constants, ONLY : fpi_const,tpi_const
   USE m_types
   USE nvtx
   IMPLICIT NONE
   TYPE(t_input),INTENT(IN)      :: input
   TYPE(t_mpi),INTENT(IN)        :: mpi
   TYPE(t_noco),INTENT(IN)       :: noco
   TYPE(t_atoms),INTENT(IN)      :: atoms
   TYPE(t_lapw),INTENT(IN)       :: lapw
   TYPE(t_usdus),INTENT(IN)      :: usdus
   CLASS(t_mat),INTENT(INOUT)     :: smat,hmat
   !     ..
   !     .. Scalar Arguments ..
   INTEGER, INTENT (IN) :: n,isp,iintsp,jintsp
   COMPLEX, INTENT(IN)  :: chi
   !     ..
   !     .. Array Arguments ..
   REAL,MANAGED,INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,input%jspins)
   REAL,    INTENT (IN) :: e_shift!(atoms%ntype,input%jspins)
   REAL,MANAGED,INTENT(IN)    :: fj(:,0:,:),gj(:,0:,:)
   !     ..
   !     .. Local Scalars ..
   INTEGER l
   INTEGER :: grid, block, loop_size

   !     ..
   !     .. Local Arrays ..
   REAL,MANAGED :: fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)
   REAL,MANAGED :: fl2p1bt(0:atoms%lmaxd)
   REAL,MANAGED :: qssbti(3),qssbtj(3)
   INTEGER, DEVICE :: nv_dev(2)

   call nvtxStartRange("hsmt_sph",2)
   CALL timestart("spherical setup")
   print*, "HsmtSph_GPU"

   DO l = 0,atoms%lmaxd
      fleg1(l) = REAL(l+l+1)/REAL(l+1)
      fleg2(l) = REAL(l)/REAL(l+1)
      fl2p1(l) = REAL(l+l+1)/fpi_const
      fl2p1bt(l) = fl2p1(l)*0.5
   END DO
   qssbti=MERGE(- noco%qss/2,+ noco%qss/2,jintsp.EQ.1)
   qssbtj=MERGE(- noco%qss/2,+ noco%qss/2,iintsp.EQ.1)

   ! pretty ugly solution
   nv_dev = lapw%nv
   loop_size = 1
   block = 32   ! number of threads in a block
   grid = ceiling(real(lapw%nv(jintsp))/(loop_SIZE*block))
   !loop_size = max(lapw%nv(jintsp)/(grid*block),1)   !number of iterations performed by each thread
   !if (loop_size * grid*block < lapw%nv(jintsp)) loop_size = loop_size + 1
   IF (input%l_useapw) THEN
      !TODO!!!!
      ! APW case is not testet
      IF (smat%l_real) THEN
         CALL HsmtSphGpuKernel_real<<<grid,block>>>(loop_size,iintsp,jintsp,nv_dev,&
                                                    atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                                            lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                                                    qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                                                    usdus%ddn(:,n,isp),el(:,n,isp),e_shift,&
                                                    smat%data_r,hmat%data_r,&
                                           usdus%uds(:,n,isp),usdus%dus(:,n,isp),usdus%us(:,n,isp),usdus%duds(:,n,isp),atoms%rmt(n))
      ELSE
         CALL HsmtSphGpuKernel_cmplx<<<grid,block>>>(loop_size,iintsp,jintsp,nv_dev,&
                                                     atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                                            lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                                                   chi,qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                                                     usdus%ddn(:,n,isp),el(:,n,isp),e_shift,&
                                                     smat%data_c,hmat%data_c,&
                                           usdus%uds(:,n,isp),usdus%dus(:,n,isp),usdus%us(:,n,isp),usdus%duds(:,n,isp),atoms%rmt(n))
      ENDIF
   ELSE
      IF (smat%l_real) THEN
         CALL HsmtSphGpuKernel_real<<<grid,block>>>(loop_size,iintsp,jintsp,nv_dev,&
                                                    atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                                            lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                                                    qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                                                    usdus%ddn(:,n,isp),el(:,n,isp),e_shift,smat%data_r,hmat%data_r)
      ELSE
         CALL HsmtSphGpuKernel_cmplx<<<grid,block>>>(loop_size,iintsp,jintsp,nv_dev,&
                                                     atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                                            lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                                                   chi,qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                                                     usdus%ddn(:,n,isp),el(:,n,isp),e_shift,smat%data_c,hmat%data_c)
      ENDIF
   ENDIF
   CALL timestop("spherical setup")

   call nvtxEndRange
   RETURN
END SUBROUTINE hsmt_sph_gpu
#endif

END MODULE m_hsmt_sph
