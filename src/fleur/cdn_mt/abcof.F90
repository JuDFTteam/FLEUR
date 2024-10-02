!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abcof

CONTAINS

  ! The subroutine abcof calculates the A, B, and C coefficients for the
  ! eigenfunctions. Also some force contributions can be calculated.
  SUBROUTINE abcof(input,atoms,sym, cell,lapw,ne,usdus,&
                   noco,nococonv,jspin , acof,bcof,ccof,zMat,eig,force,nat_start,nat_stop)
#ifdef _OPENACC
    use cublas
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
    INTEGER,OPTIONAL,INTENT(IN):: nat_start,nat_stop

    ! Local objects
    TYPE(t_fjgj) :: fjgj

    ! Local scalars
    INTEGER :: i,iLAPW,l,ll1,lm,nap,jAtom,lmp,m,nkvec,iAtom,iType,acof_size,iAtom_l,jatom_l
    INTEGER :: inv_f,ie,ilo,kspin,iintsp,nintsp,nvmax,lo,inap,abSize
    REAL    :: tmk, qss(3), s2h
    COMPLEX :: phase, c_1, c_2
    LOGICAL :: l_force,l_useinversionsym

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
    COMPLEX, ALLOCATABLE :: helpMat_force(:,:)


    CALL timestart("abcof")

    ! Checks
    IF (zmat%l_real) THEN
       IF (noco%l_noco) CALL judft_error("BUG in abcof, l_noco but real?")
    ENDIF

    ! Allocations
    CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,jspin,noco)
    ALLOCATE(abCoeffs(2*atoms%lmaxd*(atoms%lmaxd+2)+2,MAXVAL(lapw%nv)))
    ALLOCATE(abTemp(SIZE(acof,1),0:2*SIZE(acof,2)-1))
    ALLOCATE(fgpl(3,MAXVAL(lapw%nv)))
    ALLOCATE (work_c(MAXVAL(lapw%nv),ne))

    ! Initializations
    acof_size=size(acof,1)
    !$acc enter data create(abTemp,fjgj,fjgj%fj,fjgj%gj,work_c,abcoeffs)
    acof(:,:,:)   = CMPLX(0.0,0.0)
    bcof(:,:,:)   = CMPLX(0.0,0.0)
    ccof(:,:,:,:) = CMPLX(0.0,0.0)
    l_force = .FALSE.
    IF(PRESENT(eig).AND.input%l_f) l_force = .TRUE.
    IF(l_force) THEN
       force%acoflo  = CMPLX(0.0,0.0)
       force%bcoflo  = CMPLX(0.0,0.0)
       force%e1cof   = CMPLX(0.0,0.0)
       force%e2cof   = CMPLX(0.0,0.0)
       force%aveccof = CMPLX(0.0,0.0)
       force%bveccof = CMPLX(0.0,0.0)
       force%cveccof = CMPLX(0.0,0.0)
       ALLOCATE(helpMat_force(ne,atoms%lmaxd*(atoms%lmaxd+2)+1))
       ALLOCATE(workTrans_cf(ne,MAXVAL(lapw%nv)))
       ALLOCATE(s2h_e(ne,MAXVAL(lapw%nv)))
    ENDIF

    !Use inversion symmetry explicitely
    l_useinversionsym=any(sym%invsat==2)!.and.(.not.noco%l_soc).and.(.not.present(nat_start))

    
    ! loop over atoms
    DO iAtom = 1,atoms%nat
       !There might be a parallelization over atoms...
      iAtom_l=iAtom
       if (present(nat_start).and.present(nat_stop)) THEN
         if (iatom<nat_start.or.iatom>nat_stop) cycle
         iAtom_l=iAtom-nat_start+1
       endif   
       if (sym%invsat(iatom)==2.and. l_useinversionsym) cycle
       iType = atoms%itype(iAtom)

       CALL timestart("fjgj coefficients")
       CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,iType,jspin)
       !!$acc update device (fjgj%fj,fjgj%gj)
       CALL timestop("fjgj coefficients")

       CALL setabc1lo(atoms,iType,usdus,jspin,alo1,blo1,clo1)

          ! generate the spinors (chi)
       IF(noco%l_noco) ccchi=conjg(nococonv%umat(itype))


       nintsp = 1
       IF (noco%l_ss) nintsp = 2
       ! loop over the interstitial spin
       DO iintsp = 1,nintsp

          nvmax=lapw%nv(jspin)
          IF (noco%l_ss) nvmax=lapw%nv(iintsp)
          qss = MERGE(-1.0,1.0,iintsp.EQ.1)*nococonv%qss/2.0
          
     
             ! Filling of work array (modified zMat)
             CALL timestart("fill work array")
             IF (noco%l_noco) THEN
                IF (noco%l_ss) THEN
                   !$acc kernels copyin(zmat,zMat%data_c,ccchi,atoms,lapw,lapw%nv) present(work_c)default(none)
                   ! the coefficients of the spin-down basis functions are
                   ! stored in the second half of the eigenvector
                   kspin = (iintsp-1)*(lapw%nv(1)+atoms%nlotot)
                   work_c(:nvmax,:) = ccchi(iintsp,jspin)*zMat%data_c(kspin+1:kspin+nvmax,:ne)
                   !$acc end kernels
                ELSE
                   ! perform sum over the two interstitial spin directions
                   ! and take into account the spin boundary conditions
                   ! (jspin counts the local spin directions inside each MT)
                   !$acc kernels copyin(atoms,zMat,zMat%data_c,ccchi,lapw) present(work_c) default(none)
                   kspin = lapw%nv(1)+atoms%nlotot
                   work_c(:nvmax,:) = ccchi(1,jspin)*zMat%data_c(:nvmax,:ne) + ccchi(2,jspin)*zMat%data_c(kspin+1:kspin+nvmax,:ne)
                   !$acc end kernels
                END IF
             ELSE
                IF (zmat%l_real) THEN
                   !$CPP_OMP PARALLEL DO default(shared) private(i)
                   !$acc kernels copyin(zmat,zMat%data_r)present(work_c)default(none)
                   DO i = 1, ne
#ifdef _OPENACC
                      work_c(:nvmax,i) = zmat%data_r(:nvmax,i)
#else
                      work_c(:nvmax,i) = 0.0
                      CALL dcopy(nvmax,zMat%data_r(:,i),1,work_c(:,i),2)
#endif
                   END DO
                   !$acc end kernels
                   !$CPP_OMP END PARALLEL DO
                ELSE
                   !$CPP_OMP PARALLEL DO default(shared) private(i)
                   !$acc kernels copyin(zMat,zMat%data_c)present(work_c) default(none)
                   DO i = 1, ne
#ifdef _OPENACC
                      work_c(:nvmax,i) = zmat%data_c(:nvmax,i)
#else
                      CALL zcopy(nvmax,zMat%data_c(:,i),1,work_c(:,i),1)
#endif
                   END DO
                   !$acc end kernels
                   !$CPP_OMP END PARALLEL DO
                END IF
             END IF

             CALL timestop("fill work array")

             ! Calculation of a, b coefficients for LAPW basis functions
             CALL timestart("hsmt_ab")
             !!$acc data copyin(fjgj,fjgj%fj,fjgj%gj) copyout(abcoeffs)
             CALL hsmt_ab(sym,atoms,noco,nococonv,jspin,iintsp,iType,iAtom,cell,lapw,fjgj,abCoeffs,abSize,.FALSE.)
             !!$acc end data
             abSize = abSize / 2
             CALL timestop("hsmt_ab")

             ! Obtaining A, B coefficients for eigenfunctions
             CALL timestart("gemm")

             ! variant with zgemm


             !$acc host_data use_device(work_c,abCoeffs,abTemp)
             CALL zgemm_acc("T","T",ne,2*abSize,nvmax,CMPLX(1.0,0.0),work_c,MAXVAL(lapw%nv),abCoeffs,2*atoms%lmaxd*(atoms%lmaxd+2)+2,CMPLX(0.0,0.0),abTemp,acof_size)
             !$acc end host_data
             !$acc update self(abTemp)
             !stop "DEBUG"
             !$OMP PARALLEL DO default(shared) private(i,lm) collapse(2)
             DO lm = 0, absize-1
                DO i = 1, ne
                   acof(i,lm,iAtom_l) = acof(i,lm,iAtom_l) + abTemp(i,lm)
                   bcof(i,lm,iAtom_l) = bcof(i,lm,iAtom_l) + abTemp(i,absize+lm)
                END DO
             END DO
             !$OMP END PARALLEL DO

             CALL timestop("gemm")

             CALL timestart("local orbitals")
             ! Treatment of local orbitals
             !!$acc data copyin(alo1,blo1,clo1,ccchi)create(ylm)
             DO lo = 1, atoms%nlo(iType)               
                DO nkvec = 1, lapw%nkvec(lo,iAtom)
                   iLAPW = lapw%kvec(nkvec,lo,iAtom)
                   fg(:) = MERGE(lapw%gvec(:,iLAPW,iintsp),lapw%gvec(:,iLAPW,jspin),noco%l_ss) + qss + lapw%qPhon
                   fk = lapw%bkpt + fg(:)
                   tmk = tpi_const * DOT_PRODUCT(fk(:),atoms%taual(:,iAtom))
                   phase = CMPLX(COS(tmk),SIN(tmk))

                    
                   nap = sym%ngopr(iAtom)
                   inap = sym%invtab(nap)
                   fkr = MATMUL(TRANSPOSE(sym%mrot(:,:,inap)),fk(:))
                   fgr = MATMUL(TRANSPOSE(sym%mrot(:,:,inap)),fg(:))
                   
                   fkp = MATMUL(fkr,cell%bmat)
                   fgp = MATMUL(fgr,cell%bmat)

                   CALL ylm4(atoms%lmax(iType),fkp,ylm)
                   !!$acc update device(ylm)
                   CALL abclocdn(atoms,noco,lapw,cell,ccchi(:,jspin),iintsp,phase,ylm,iType,iAtom,iLAPW,nkvec,&
                                 lo,ne,alo1(:,jspin),blo1(:,jspin),clo1(:,jspin),acof,bcof,ccof,zMat,l_force,fgp,force,iAtom_l)
                END DO
             END DO ! loop over LOs
             !!$acc end data
             CALL timestop("local orbitals")

            
             ! Force contributions
             IF (atoms%l_geo(iType).AND.l_force) THEN
               !$acc  update self(abcoeffs,work_c)
               CALL timestart("transpose work array")
               ! For transposing the work array an OpenMP parallelization with explicit loops is used.
               ! This solution works fastest on all compilers. Note that this section can actually be
               ! a bottleneck without parallelization if many OpenMP threads are used.
               IF (zmat%l_real) THEN
                  ALLOCATE (workTrans_r(ne,nvmax))
                  !$OMP PARALLEL DO default(shared) private(i,iLAPW) collapse(2)
                  DO i = 1,ne
                     DO iLAPW = 1, nvmax
                        workTrans_r(i,iLAPW) = work_c(iLAPW,i)
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
               ELSE
                  ALLOCATE (workTrans_c(ne,nvmax))
                  !$OMP PARALLEL DO default(shared) private(i,iLAPW) collapse(2)
                  DO i = 1,ne
                     DO iLAPW = 1, nvmax
                        workTrans_c(i,iLAPW) = work_c(iLAPW,i)
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
               ENDIF
               CALL timestop("transpose work array")
            
                CALL timestart("force contributions")
                DO iLAPW = 1,nvmax

                   fg(:) = MERGE(lapw%gvec(:,iLAPW,iintsp),lapw%gvec(:,iLAPW,jspin),noco%l_ss) + qss
                   fk = lapw%bkpt + fg(:)
                   s2h = 0.5 * DOT_PRODUCT(fk,MATMUL(cell%bbmat,fk))
                   IF (zmat%l_real) THEN
                      s2h_e(:ne,iLAPW) = CMPLX((s2h-eig(:ne)) * workTrans_r(:ne,iLAPW))
                   ELSE
                      s2h_e(:ne,iLAPW) = (s2h-eig(:ne)) * workTrans_c(:ne,iLAPW)
                   ENDIF
                    
                      nap = sym%ngopr(iAtom)
                      inap = sym%invtab(nap)
                      fgr = MATMUL(TRANSPOSE(sym%mrot(:,:,inap)),fg(:))
                   
                   fgpl(:,iLAPW) = MATMUL(fgr,cell%bmat)
                ENDDO

                workTrans_cf = CMPLX(0.0,0.0)
                helpMat_force = CMPLX(0.0,0.0)

                CALL zgemm("N","T",ne,abSize,nvmax,CMPLX(1.0,0.0),s2h_e,ne,abCoeffs,size(abcoeffs,1),CMPLX(1.0,0.0),force%e1cof(:,:,iAtom),ne)
                CALL zgemm("N","T",ne,abSize,nvmax,CMPLX(1.0,0.0),s2h_e,ne,abCoeffs(1+abSize,1),size(abcoeffs,1),CMPLX(1.0,0.0),force%e2cof(:,:,iAtom),ne)
                DO i =1,3
                   IF (zmat%l_real) THEN
                      DO iLAPW = 1,nvmax
                         workTrans_cf(:,iLAPW) = CMPLX(workTrans_r(:,iLAPW) * fgpl(i,iLAPW))
                      ENDDO
                   ELSE
                      DO iLAPW = 1,nvmax
                         workTrans_cf(:,iLAPW) = workTrans_c(:,iLAPW) * fgpl(i,iLAPW)
                      ENDDO
                   ENDIF
                   CALL zgemm("N","T",ne,abSize,nvmax,CMPLX(1.0,0.0),workTrans_cf,ne,abCoeffs,size(abCoeffs,1),CMPLX(0.0,0.0),helpMat_force,ne)
                   force%aveccof(i,:,:,iAtom) = force%aveccof(i,:,:,iAtom) + helpMat_force(:,:)
                   CALL zgemm("N","T",ne,abSize,nvmax,CMPLX(1.0,0.0),workTrans_cf,ne,abCoeffs(1+abSize,1),size(abcoeffs,1),CMPLX(0.0,0.0),helpMat_force,ne)
                   force%bveccof(i,:,:,iAtom) = force%bveccof(i,:,:,iAtom) + helpMat_force(:,:)
                ENDDO
                CALL timestop("force contributions")
             END IF
             IF (atoms%l_geo(iType).AND.l_force) THEN
                IF (zmat%l_real) THEN
                   DEALLOCATE (workTrans_r)
                ELSE
                   DEALLOCATE (workTrans_c)
                ENDIF
             END IF
       END DO ! loop over interstitial spin
    END DO ! loop over atoms
    !$acc exit data delete(abTemp,fjgj%fj,fjgj%gj,work_c,abcoeffs)
    !$acc exit data delete(fjgj)
    DEALLOCATE(work_c)
    IF(l_force) THEN
       DEALLOCATE(helpMat_force)
       DEALLOCATE(workTrans_cf)
       DEALLOCATE(s2h_e)
    ENDIF

    ! Treatment of atoms inversion symmetric to others
    IF (l_useinversionsym) THEN
       !Comment on SOC case:
       !
       !                           -p,n       (l+m)   p,n  *
       ! Usually, we exploit that A     = (-1)      (A    )  if p and -p are the positions
       !                           l,m                l,-m
       ! of two atoms related by inversion symmetry and the coefficients are considered to
       ! be in the local frame of the representative atom. This is possible, if z is real.
       ! After SOC, however, the eigenvectors z are complex and this is no longer possible
       ! so the z has to enter, not z*. This is done within the k-loop.
       !                                    -p,n       m   p,n  *
       ! When called from hsohelp, we need A     = (-1)  (A    ) because we don't have to
       !                                     l,m           l,-m                    l
       ! rotate, but in the sums in hsoham only products A*  A   enter and the (-1) cancels.
       !                                                  lm  lm
       DO iAtom = 1, atoms%nat
          iType = atoms%itype(iAtom)
          iAtom_l=iAtom
          if (present(nat_start).and.present(nat_stop)) THEN
            if (iatom<nat_start.or.iatom>nat_stop) cycle
            iAtom_l=iAtom-nat_start+1
          endif   
          IF (sym%invsat(iAtom).EQ.1) THEN
             jAtom = sym%invsatnr(iAtom)
             jatom_l=jatom
             if (present(nat_start).and.present(nat_stop)) THEN
               if (jatom<nat_start.or.jatom>nat_stop) call judft_bug("MPI distribution failed in 2nd variation SOC")
               jAtom_l=jAtom-nat_start+1
             endif
             DO ilo = 1,atoms%nlo(iType)
                l = atoms%llo(ilo,iType)
                DO m = -l,l
                   inv_f = (-1)**(m+l)
                   DO ie = 1,ne
                      ccof(m,ie,ilo,jatom_l) = inv_f * CONJG( ccof(-m,ie,ilo,iatom_l))
                      IF(l_force) THEN
                         force%acoflo(m,ie,ilo,jatom_l) = inv_f * CONJG(force%acoflo(-m,ie,ilo,iatom_l))
                         force%bcoflo(m,ie,ilo,jatom_l) = inv_f * CONJG(force%bcoflo(-m,ie,ilo,iatom_l))
                         force%cveccof(:,m,ie,ilo,jatom_l) = -inv_f * CONJG(force%cveccof(:,-m,ie,ilo,iatom_l))
                      END IF
                   END DO
                END DO
             END DO
             DO l = 0,atoms%lmax(iType)
                ll1 = l* (l+1)
                DO m =-l,l
                   lm  = ll1 + m
                   lmp = ll1 - m
                   inv_f = (-1)**(m+l)
                   acof(:ne,lm,jatom_l) = inv_f * CONJG(acof(:ne,lmp,iatom_l))
                   bcof(:ne,lm,jatom_l) = inv_f * CONJG(bcof(:ne,lmp,iatom_l))
                   IF (atoms%l_geo(iType).AND.l_force) THEN
                      force%e1cof(:ne,lm,jatom_l) = inv_f * CONJG(force%e1cof(:ne,lmp,iatom_l))
                      force%e2cof(:ne,lm,jatom_l) = inv_f * CONJG(force%e2cof(:ne,lmp,iatom_l))
                      force%aveccof(:,:ne,lm,jatom_l) = -inv_f * CONJG(force%aveccof(:,:ne,lmp,iatom_l))
                      force%bveccof(:,:ne,lm,jatom_l) = -inv_f * CONJG(force%bveccof(:,:ne,lmp,iatom_l))
                   END IF
                END DO
             END DO
             
          END IF
       END DO
    END IF

    CALL timestop("abcof")

  END SUBROUTINE abcof
END MODULE m_abcof
