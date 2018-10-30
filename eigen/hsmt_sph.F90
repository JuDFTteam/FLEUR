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

#ifdef CPP_GPU
  ATTRIBUTES(global) SUBROUTINE HsmtSphGpuKernel_real(grid,block,iintsp,jintsp,nv,lmaxd,lmax,ki_start,ki_end,ki_step,nn_start,nn_end,&
                       lnonsph,qssbti,qssbtj,gvec,gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,taual,ddn,el,e_shift,&
                       smat_data,hmat_data,&
                       uds,dus,us,duds,rmt)
     INTEGER, VALUE, INTENT(IN) :: grid, block
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
     REAL, ALLOCATABLE :: plegend(:,:)
     COMPLEX, ALLOCATABLE :: cph(:)
     REAL tnn(3), elall,fct,fct2,fjkiln,gjkiln,ddnln,ski(3)
     REAL apw_lo1,apw_lo2,apw1,w1
     INTEGER kii,ki,kj,l,nn,k
     INTEGER :: loop_start, loop_end, i, loop_size

     ALLOCATE(cph(MAXVAL(nv)))
     ALLOCATE(plegend(MAXVAL(nv),0:lmaxd))
     plegend=0.0
     plegend(:,0)=1.0

     k = (blockidx%x-1)*blockdim%x + threadidx%x
 
     !TODO!!!     
     !for seq, i.e. ki_start = 1, ki_step = 1
     loop_size = max(ki_end/(grid*block),1)
     if (loop_size * grid*block < ki_end) loop_size = loop_size + 1
     loop_start = (k-1) * loop_size + 1
     loop_end = loop_start + loop_size - 1
     if (loop_end > ki_end ) loop_end = ki_end

     DO ki = loop_start,loop_end,ki_step
     !DO  ki =  ki_start,ki_end,ki_step 
       kii=(ki-1)/ki_step+1
       ski = gvec(:,ki,jintsp) + qssbti
       !--->       legendre polynomials
       DO kj = 1,ki
          plegend(kj,1) = DOT_PRODUCT(gk(:,kj,iintsp),gk(:,ki,jintsp))
       END DO
       DO l = 1,lmax - 1
          plegend(:ki,l+1) = fleg1(l)*plegend(:ki,1)*plegend(:ki,l) - fleg2(l)*plegend(:ki,l-1)
       END DO
       !--->             set up phase factors
       cph = 0.0
       DO nn = nn_start,nn_end
          tnn = tpi_const*taual(:,nn)
          DO kj = 1,ki
             cph(kj) = cph(kj) +&
                  CMPLX(COS(DOT_PRODUCT(ski-gvec(:,kj,iintsp)-qssbtj,tnn)),&
                  SIN(DOT_PRODUCT(gvec(:,kj,iintsp)+qssbtj-ski,tnn)))
            ! IF (iintsp.NE.jintsp) cph(kj)=CONJG(cph(kj))
          END DO
       END DO

       !--->          update overlap and l-diagonal hamiltonian matrix
       DO  l = 0,lmax
          !+APW
          IF (PRESENT(uds)) THEN
             w1 = 0.5 * ( uds(l)*dus(l) + us(l)*duds(l) )
             apw_lo1 = fl2p1(l) * 0.5 * rmt**2 * ( gjkiln * w1 +&
                  fjkiln * us(l) * dus(l) )
             apw_lo2 = fl2p1(l) * 0.5 * rmt**2 * ( fjkiln * w1 +&
                  gjkiln * uds(l) * duds(l) )
          ENDIF
          !-APW
          fjkiln = fj(ki,l,jintsp)
          gjkiln = gj(ki,l,jintsp)
          ddnln =  ddn(l)
          elall = el(l)
          IF (l<=lnonsph) elall=elall-e_shift!(isp)
          DO kj = 1,ki
              fct  = plegend(kj,l)*fl2p1(l)*&
                     ( fjkiln*fj(kj,l,iintsp) + gjkiln*gj(kj,l,iintsp)*ddnln )
              fct2 = plegend(kj,l)*fl2p1bt(l) * ( fjkiln*gj(kj,l,iintsp) + gjkiln*fj(kj,l,iintsp) )

              smat_data(kj,kii)=smat_data(kj,kii)+REAL(cph(kj))*fct
              hmat_data(kj,kii)=hmat_data(kj,kii) + REAL(cph(kj)) * ( fct * elall + fct2)
              !+APW
              IF (PRESENT(uds)) THEN
                 apw1 = REAL(cph(kj)) * plegend(kj,l)  * &
                        ( apw_lo1 * fj(kj,l,iintsp) + apw_lo2 * gj(kj,l,iintsp) )
                 hmat_data(kj,kii)=hmat_data(kj,kii) + apw1
              ENDIF
              !-APW
          ENDDO
       !--->          end loop over l
       ENDDO
     !--->    end loop over ki
     ENDDO
     DEALLOCATE(plegend)
     DEALLOCATE(cph)
  END SUBROUTINE HsmtSphGpuKernel_real

  ATTRIBUTES(global) SUBROUTINE HsmtSphGpuKernel_cmplx(grid,block,iintsp,jintsp,nv,lmaxd,lmax,ki_start,ki_end,ki_step,nn_start,nn_end,&
                       lnonsph,chi,qssbti,qssbtj,gvec,gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,taual,ddn,el,e_shift,&
                       smat_data,hmat_data,&
                       uds,dus,us,duds,rmt)
     INTEGER, VALUE, INTENT(IN) :: grid, block
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
     REAL, ALLOCATABLE :: plegend(:,:)
     REAL, ALLOCATABLE :: VecHelpS(:),VecHelpH(:)
     COMPLEX, ALLOCATABLE :: cph(:)
     REAL apw_lo1,apw_lo2,w1
     COMPLEX capw1
     REAL tnn(3), elall,fct,fct2,fjkiln,gjkiln,ddnln,ski(3)
     INTEGER kii,ki,kj,l,nn,kj_end,k
     INTEGER :: loop_start, loop_end, i, loop_size

     ALLOCATE(cph(MAXVAL(nv)))
     ALLOCATE(plegend(MAXVAL(nv),0:lmaxd))
     plegend=0.0
     plegend(:,0)=1.0
     ALLOCATE(VecHelpS(MAXVAL(nv)),VecHelpH(MAXVAL(nv)))

     k = (blockidx%x-1)*blockdim%x + threadidx%x
 
     !TODO!!!     
     !for seq, i.e. ki_start = 1, ki_step = 1
     loop_size = max(ki_end/(grid*block),1)
     if (loop_size * grid*block < ki_end) loop_size = loop_size + 1
     loop_start = (k-1) * loop_size + 1
     loop_end = loop_start + loop_size - 1
     if (loop_end > ki_end ) loop_end = ki_end

     DO ki = loop_start,loop_end,ki_step
     !DO  ki =  ki_start,ki_end,ki_step 
       kii=(ki-1)/ki_step+1
       ski = gvec(:,ki,jintsp) + qssbti
       !--->       legendre polynomials
       DO kj = 1,ki
          plegend(kj,1) = DOT_PRODUCT(gk(:,kj,iintsp),gk(:,ki,jintsp))
       END DO
       DO l = 1,lmax - 1
          plegend(:ki,l+1) = fleg1(l)*plegend(:ki,1)*plegend(:ki,l) - fleg2(l)*plegend(:ki,l-1)
       END DO
       !--->             set up phase factors
       cph = 0.0
       DO nn = nn_start,nn_end
          tnn = tpi_const*taual(:,nn)
          DO kj = 1,ki
             cph(kj) = cph(kj) +&
                  CMPLX(COS(DOT_PRODUCT(ski-gvec(:,kj,iintsp)-qssbtj,tnn)),&
                  SIN(DOT_PRODUCT(gvec(:,kj,iintsp)+qssbtj-ski,tnn)))
            ! IF (iintsp.NE.jintsp) cph(kj)=CONJG(cph(kj))
          END DO
       END DO

       !--->          update overlap and l-diagonal hamiltonian matrix
       kj_end = MIN(ki,nv(iintsp)) 
       VecHelpS = 0.d0
       VecHelpH = 0.d0
       DO  l = 0,lmax
          !+APW
          IF (PRESENT(uds)) THEN
             w1 = 0.5 * ( uds(l)*dus(l) + us(l)*duds(l) )
             apw_lo1 = fl2p1(l) * 0.5 * rmt**2 * ( gjkiln * w1 +&
                  fjkiln * us(l) * dus(l) )
             apw_lo2 = fl2p1(l) * 0.5 * rmt**2 * ( fjkiln * w1 +&
                  gjkiln * uds(l) * duds(l) )
          ENDIF
          !-APW
          fjkiln = fj(ki,l,jintsp)
          gjkiln = gj(ki,l,jintsp)
          ddnln =  ddn(l)
          elall = el(l)
          IF (l<=lnonsph) elall=elall-e_shift!(isp)
          DO kj = 1,kj_end
              fct  = plegend(kj,l)*fl2p1(l)*&
                     ( fjkiln*fj(kj,l,iintsp) + gjkiln*gj(kj,l,iintsp)*ddnln )
              fct2 = plegend(kj,l)*fl2p1bt(l) * ( fjkiln*gj(kj,l,iintsp) + gjkiln*fj(kj,l,iintsp) )
 
              VecHelpS(kj) = VecHelpS(kj) + fct
              VecHelpH(kj) = VecHelpH(kj) + fct*elall + fct2 
              !+APW
              IF (PRESENT(uds)) THEN
                 capw1 = cph(kj)*plegend(kj,l)&
                        * ( apw_lo1 * fj(kj,l,iintsp) + apw_lo2 * gj(kj,l,iintsp) )
                 hmat_data(kj,kii)=hmat_data(kj,kii) + capw1
              ENDIF
              !-APW
          END DO
       !--->          end loop over l
       ENDDO
       smat_data(:kj_end,kii)=smat_data(:kj_end,kii) + chi*cph(:kj_end) * VecHelpS(:kj_end)
       hmat_data(:kj_end,kii)=hmat_data(:kj_end,kii) + chi*cph(:kj_end) * VecHelpH(:kj_end)
     !--->    end loop over ki
     ENDDO
     DEALLOCATE(plegend)
     DEALLOCATE(cph)
     DEALLOCATE(VecHelpS,VecHelpH)
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
    INTEGER :: grid, block

    !     ..
    !     .. Local Arrays ..
    REAL,MANAGED :: fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL,MANAGED :: fl2p1bt(0:atoms%lmaxd)
    REAL,MANAGED :: qssbti(3),qssbtj(3)
    INTEGER, DEVICE :: nv_dev(2)

    call nvtxStartRange("hsmt_sph",2)    
    CALL timestart("spherical setup")

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
    block = 256
    grid = lapw%nv(jintsp)/(block*4) + 1
    IF (input%l_useapw) THEN
       !TODO!!!!
       ! APW case is not testet
       IF (smat%l_real) THEN
          CALL HsmtSphGpuKernel_real<<<grid,block>>>(grid,block,iintsp,jintsp,nv_dev,atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                    lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                    qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                    usdus%ddn(:,n,isp),el(:,n,isp),e_shift,&
                    smat%data_r,hmat%data_r,&
                    usdus%uds(:,n,isp),usdus%dus(:,n,isp),usdus%us(:,n,isp),usdus%duds(:,n,isp),atoms%rmt(n))
       ELSE
          CALL HsmtSphGpuKernel_cmplx<<<grid,block>>>(grid,block,iintsp,jintsp,nv_dev,atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                    lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                    chi,qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                    usdus%ddn(:,n,isp),el(:,n,isp),e_shift,&
                    smat%data_c,hmat%data_c,&
                    usdus%uds(:,n,isp),usdus%dus(:,n,isp),usdus%us(:,n,isp),usdus%duds(:,n,isp),atoms%rmt(n))
       ENDIF
    ELSE
       IF (smat%l_real) THEN
          CALL HsmtSphGpuKernel_real<<<grid,block>>>(grid,block,iintsp,jintsp,nv_dev,atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
                    lapw%nv(jintsp), mpi%n_size,SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n)),atoms%lnonsph(n),&
                    qssbti,qssbtj,lapw%gvec,lapw%gk,fleg1,fleg2,fl2p1,fl2p1bt,fj,gj,atoms%taual,&
                    usdus%ddn(:,n,isp),el(:,n,isp),e_shift,smat%data_r,hmat%data_r)
       ELSE
          CALL HsmtSphGpuKernel_cmplx<<<grid,block>>>(grid,block,iintsp,jintsp,nv_dev,atoms%lmaxd,atoms%lmax(n),mpi%n_rank+1,&
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
    CLASS(t_mat),INTENT(INOUT)     :: smat,hmat
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
    REAL tnn(3), elall,fct,fjkiln,gjkiln,ddnln,ski(3)
    REAL apw_lo1,apw_lo2,apw1,w1

    COMPLEX capw1
    INTEGER kii,ki,kj,l,nn

    !     ..
    !     .. Local Arrays ..
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL fl2p1bt(0:atoms%lmaxd)
    REAL qssbti(3),qssbtj(3)
    REAL, ALLOCATABLE :: plegend(:,:)
    COMPLEX, ALLOCATABLE :: cph(:)
    LOGICAL apw(0:atoms%lmaxd)

    CALL timestart("spherical setup")

    DO l = 0,atoms%lmaxd
       fleg1(l) = REAL(l+l+1)/REAL(l+1)
       fleg2(l) = REAL(l)/REAL(l+1)
       fl2p1(l) = REAL(l+l+1)/fpi_const
       fl2p1bt(l) = fl2p1(l)*0.5
    END DO
    !$OMP PARALLEL DEFAULT(SHARED)&
    !$OMP PRIVATE(kii,ki,ski,kj,plegend,l)&
    !$OMP PRIVATE(cph,nn,tnn,fjkiln,gjkiln)&
    !$OMP PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct,apw1)&
    !$OMP PRIVATE(capw1) 
    ALLOCATE(cph(MAXVAL(lapw%nv)))
    ALLOCATE(plegend(MAXVAL(lapw%nv),0:atoms%lmaxd))
    plegend=0.0
    plegend(:,0)=1.0
    qssbti=MERGE(- noco%qss/2,+ noco%qss/2,jintsp.EQ.1)
    qssbtj=MERGE(- noco%qss/2,+ noco%qss/2,iintsp.EQ.1)
    !$OMP  DO SCHEDULE(DYNAMIC,1)
    DO  ki =  mpi%n_rank+1, lapw%nv(jintsp), mpi%n_size
       kii=(ki-1)/mpi%n_size+1
       ski = lapw%gvec(:,ki,jintsp) + qssbti
       !--->       legendre polynomials
       DO kj = 1,ki
          plegend(kj,1) = DOT_PRODUCT(lapw%gk(:,kj,iintsp),lapw%gk(:,ki,jintsp))
       END DO
       DO l = 1,atoms%lmax(n) - 1
          plegend(:ki,l+1) = fleg1(l)*plegend(:ki,1)*plegend(:ki,l) - fleg2(l)*plegend(:ki,l-1)
       END DO
       !--->             set up phase factors
       cph = 0.0
       DO nn = SUM(atoms%neq(:n-1))+1,SUM(atoms%neq(:n))
          tnn = tpi_const*atoms%taual(:,nn)
          DO kj = 1,ki
             cph(kj) = cph(kj) +&
                  CMPLX(COS(DOT_PRODUCT(ski-lapw%gvec(:,kj,iintsp)-qssbtj,tnn)),&
                  SIN(DOT_PRODUCT(lapw%gvec(:,kj,iintsp)+qssbtj-ski,tnn)))
            ! IF (iintsp.NE.jintsp) cph(kj)=CONJG(cph(kj))
          END DO
       END DO

       !--->          update overlap and l-diagonal hamiltonian matrix
       DO  l = 0,atoms%lmax(n)
          fjkiln = fj(ki,l,jintsp)
          gjkiln = gj(ki,l,jintsp)
          !
          IF (input%l_useapw) THEN
             w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + &
                  usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
             apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 +&
                  fjkiln * usdus%us(l,n,isp) * usdus%dus(l,n,isp) )
             apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 +&
                  gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
             !
          ENDIF
          ddnln =  usdus%ddn(l,n,isp)
          elall = el(l,n,isp)
          IF (l<=atoms%lnonsph(n)) elall=elall-e_shift!(isp)
          IF (smat%l_real) THEN
             DO kj = 1,ki
                fct  = plegend(kj,l)*fl2p1(l)*&
                     ( fjkiln*fj(kj,l,iintsp) + gjkiln*gj(kj,l,iintsp)*ddnln )
                smat%data_r(kj,kii)=smat%data_r(kj,kii)+REAL(cph(kj))*fct
                hmat%data_r(kj,kii)=hmat%data_r(kj,kii) + REAL(cph(kj)) * &
                     ( fct * elall + plegend(kj,l) * fl2p1bt(l) *&
                     ( fjkiln*gj(kj,l,iintsp) + gjkiln*fj(kj,l,iintsp) ) )
                !+APW
                IF (input%l_useapw) THEN
                   apw1 = REAL(cph(kj)) * plegend(kj,l)  * &
                        ( apw_lo1 * fj(kj,l,iintsp) + apw_lo2 * gj(kj,l,iintsp) )
                   hmat%data_r(kj,kii)=hmat%data_r(kj,kii) + apw1
                ENDIF
                !-APW
             ENDDO
          ELSE
             DO kj = 1,MIN(ki,lapw%nv(iintsp))
                fct  = chi*plegend(kj,l)*fl2p1(l)*&
                     ( fjkiln*fj(kj,l,iintsp) + gjkiln*gj(kj,l,iintsp)*ddnln )

                smat%data_c(kj,kii)=smat%data_c(kj,kii) + cph(kj)*fct
                hmat%data_c(kj,kii)=hmat%data_c(kj,kii) + cph(kj) * ( fct*elall &
                     + chi*plegend(kj,l)*fl2p1bt(l) * ( fjkiln*gj(kj,l,iintsp) + gjkiln*fj(kj,l,iintsp) ) )

                IF (input%l_useapw) THEN
                   capw1 = cph(kj)*plegend(kj,l)&
                        * ( apw_lo1 * fj(kj,l,iintsp) + apw_lo2 * gj(kj,l,iintsp) )
                   hmat%data_c(kj,kii)=hmat%data_c(kj,kii) + capw1
                ENDIF
             END DO

          ENDIF
       !--->          end loop over l
       ENDDO
    !--->    end loop over ki
    ENDDO
    !$OMP END DO
    DEALLOCATE(plegend)
    DEALLOCATE(cph)
    DEALLOCATE(VecHelpS,VecHelpH)
    !$OMP END PARALLEL
    CALL timestop("spherical setup")

    RETURN
  END SUBROUTINE hsmt_sph_cpu
END MODULE m_hsmt_sph
