!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abcof
   use openacc
   use omp_lib
   use cublas_v2
   use cudafor
CONTAINS

  ! The subroutine abcof calculates the A, B, and C coefficients for the
  ! eigenfunctions. Also some force contributions can be calculated.
  SUBROUTINE abcof(input,atoms,sym, cell,lapw,ne,usdus,&
                   noco,nococonv,jspin,oneD, acof,bcof,ccof,zMat,eig,force)
#ifdef _OPENACC
#ifdef __PGI
!    use cublas
#endif
#define CPP_ACC acc
#define CPP_OMP no_OMP_used
#define zgemm_acc cublaszgemm
#else
#define CPP_ACC No_acc_used
#define CPP_OMP OMP
#define zgemm_acc zgemm
#endif
    USE m_juDFT
    USE m_types
    USE m_constants
    USE m_ylm
    USE m_setabc1lo
    USE m_abclocdn
    USE m_hsmt_fjgj
    USE m_hsmt_ab

    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)             :: input
    TYPE(t_usdus),INTENT(IN)             :: usdus
    TYPE(t_lapw),INTENT(IN)              :: lapw
    TYPE(t_oneD),INTENT(IN)              :: oneD
    TYPE(t_noco),INTENT(IN)              :: noco
    TYPE(t_nococonv),INTENT(IN)          :: nococonv
    TYPE(t_sym),INTENT(IN)               :: sym
    TYPE(t_cell),INTENT(IN)              :: cell
    TYPE(t_atoms),INTENT(IN)             :: atoms
    TYPE(t_mat),INTENT(IN)               :: zMat
    TYPE(t_force),OPTIONAL,INTENT(INOUT) :: force
    


    ! scalar arguments
    INTEGER, INTENT(IN)        :: ne
    INTEGER, INTENT(IN)        :: jspin

    ! array arguments
    COMPLEX, INTENT(OUT)       :: acof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT(OUT)       :: bcof(:,0:,:)!(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
    COMPLEX, INTENT(OUT)       :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    REAL, OPTIONAL, INTENT(IN) :: eig(:)!(input%neig)

    ! Local objects
    TYPE(t_fjgj) :: fjgj

    ! Local scalars
    INTEGER :: i,iLAPW,l,ll1,lm,nap,jAtom,lmp,m,nkvec,iAtom,iType,acof_size
    INTEGER :: inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap,abSize, num_streams
    integer(kind=cuda_stream_kind) :: stream
    REAL    :: tmk, qss(3), s2h
    COMPLEX :: phase, c_1
    complex :: zero = (0.0, 0.0), one = (1.0, 0.0)
    LOGICAL :: l_force

    ! Local arrays
    REAL    :: fg(3),fgp(3),fgr(3),fk(3),fkp(3),fkr(3)
    REAL    :: alo1(atoms%nlod,input%jspins),blo1(atoms%nlod,input%jspins)
    REAL    :: clo1(atoms%nlod,input%jspins)
    COMPLEX :: ylm((atoms%lmaxd+1)**2)
    COMPLEX :: ccchi(2,2)
    REAL,    ALLOCATABLE :: realCoeffs(:,:), imagCoeffs(:,:), workTrans_r(:,:)
    REAL,    ALLOCATABLE :: fgpl(:,:)
    COMPLEX, ALLOCATABLE :: s2h_e(:,:)
    COMPLEX, ALLOCATABLE :: work_c(:,:), workTrans_c(:,:), workTrans_cf(:,:)
    COMPLEX, ALLOCATABLE :: abCoeffs(:,:)
    COMPLEX, ALLOCATABLE :: abTemp(:,:)
    COMPLEX, ALLOCATABLE :: helpMat_c(:,:), helpMat_force(:,:)


  integer(kind=4) :: ierr
  integer :: async_num
  type(cublasHandle) :: handle
  complex, allocatable :: a(:), b(:), c(:)
  complex :: alpha, beta
  integer(kind=4), parameter :: n = 1024
  integer(kind=cuda_count_kind), parameter :: cuda_stack_size = 8192*4

    CALL timestart("abcof")


    ! Initializations
    acof_size=size(acof,1)
    !$acc data copyout(acof, bcof, ccof) copyin(zero, one, zmat, zmat%data_c)
    !$acc kernels present(acof,bcof,ccof) default(none)
    acof(:,:,:)   = CMPLX(0.0,0.0)
    bcof(:,:,:)   = CMPLX(0.0,0.0)
    ccof(:,:,:,:) = CMPLX(0.0,0.0)
    !$acc end kernels
    l_force = .FALSE.

    ierr = cudaDeviceSetLimit(cudaLimitStackSize, cuda_stack_size)
    num_streams = 2

    ! loop over atoms
    !$OMP parallel default(none) num_threads(num_streams) private(fjgj, abCoeffs, abTemp, fgpl, work_c, handle, stream, ierr, itype) &
    !$OMP private(absize, async_num, iAtom, i, lm)&
    !$OMP shared(ne, acof, bcof, ccof, lapw, input, atoms, cell, sym, one, zero, jspin, usdus, iintsp, nvmax, noco, nococonv, zmat, &
    !$OMP acof_size)
    ! Allocations
    CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,jspin,noco)
    ALLOCATE(abCoeffs(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
    ALLOCATE(abTemp(SIZE(acof,1),0:2*SIZE(acof,2)-1))
    ALLOCATE(fgpl(3,MAXVAL(lapw%nv)))
    ALLOCATE (work_c(MAXVAL(lapw%nv),ne))
        
    ierr = cublascreate(handle)
    ! ierr = cudaStreamCreate(stream)
    async_num = omp_get_thread_num()
    stream = acc_get_cuda_stream(async_num)
    ierr = cublasSetStream(handle, stream)

    write (*,*) "Flag 1", async_num

    !$acc data create(abTemp,fjgj,fjgj%fj,fjgj%gj,work_c,abcoeffs) async(async_num)
    !$OMP do
    DO iAtom = 1, atoms%nat
       iType = atoms%itype(iAtom)

       CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,iType,jspin)
       !$acc update device (fjgj%fj,fjgj%gj) async(async_num)

       iintsp = 1 !Only one spin
       nvmax=lapw%nv(jspin)
       ! Filling of work array (modified zMat)
       !$acc kernels present(zmat, zmat%data_c, work_c) default(none) async(async_num)
       DO i = 1, ne
         work_c(:nvmax,i) = zmat%data_c(:nvmax,i)
       END DO
       !$acc end kernels

       write (*,*) "Flag 2", async_num, iAtom
       ! Calculation of a, b coefficients for LAPW basis functions
       CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,iintsp,iType,iAtom,cell,lapw,fjgj,abCoeffs,abSize,.FALSE., p_async_num=async_num)
       abSize = abSize / 2
!
!       ! Obtaining A, B coefficients for eigenfunctions
!
!       !$acc host_data use_device(work_c,abCoeffs,abTemp)
!       ierr= cublaszgemm_v2(handle, CUBLAS_OP_T,CUBLAS_OP_C,ne,2*abSize,nvmax,one,   work_c,MAXVAL(lapw%nv),abCoeffs,2*atoms%lmaxd*(atoms%lmaxd+2)+2,zero, abTemp, acof_size)
!       !$acc end host_data
!
!      !$acc kernels present(acof,bcof,abTemp) default(none) async(async_num)
!      DO lm = 0, absize-1
!        DO i = 1, ne
!          acof(i,lm,iAtom) = acof(i,lm,iAtom) + abTemp(i,lm)
!          bcof(i,lm,iAtom) = bcof(i,lm,iAtom) + abTemp(i,absize+lm)
!        END DO
!      END DO
!      !$acc end kernels
    END DO ! loop over atoms
    !$omp end do

    !$acc end data
    
    DEALLOCATE(abCoeffs)
    DEALLOCATE(abTemp)
    DEALLOCATE(fgpl)
    DEALLOCATE(work_c)
    !$acc wait(async_num)
    !$omp end parallel

    !$acc end data 
    CALL timestop("abcof")

  END SUBROUTINE abcof
END MODULE m_abcof
