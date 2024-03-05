!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_eigen

#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT

#ifdef _OPENACC
   USE cublas
#define CPP_zgemv cublaszgemv
#else
#define CPP_zgemv zgemv
#endif

   IMPLICIT NONE

CONTAINS

   SUBROUTINE dfpt_eigen(fi, sphhar, results, resultsq, results1, fmpi, enpara, nococonv, starsq, v1real, v1imag, vTot, inden, bqpt, &
                             eig_id, q_eig_id, dfpt_eig_id, iDir, iDtype, killcont, l_real, sh_den, dfpt_eig_id2)

      USE m_types
      USE m_constants
      USE m_dfpt_eigen_hssetup
      USE m_pot_io
      USE m_util
      USE m_eig66_io, ONLY : write_eig, read_eig
      USE m_xmlOutput
      USE m_types_mpimat
      USE m_dfpt_tlmplm
      USE m_local_hamiltonian

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_sphhar),     INTENT(IN)       :: sphhar
      TYPE(t_results),INTENT(INOUT):: results, resultsq, results1
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_enpara),INTENT(IN) :: enpara
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      TYPE(t_stars),INTENT(IN)     :: starsq
      TYPE(t_potden),INTENT(IN)    :: inden, v1real, v1imag, vTot
      REAL,         INTENT(IN)     :: bqpt(3)
      INTEGER,      INTENT(IN)     :: eig_id, q_eig_id, dfpt_eig_id, iDir, iDtype, killcont(6)
      LOGICAL,      INTENT(IN)     :: l_real, sh_den
      INTEGER, OPTIONAL, INTENT(IN) :: dfpt_eig_id2

      INTEGER n_size,n_rank
      INTEGER i,err,nk,jsp,nk_i,neigd2

      INTEGER              :: ierr, iNupr

      REAL :: bkpt(3), q_loop(3)

      INTEGER                   :: nu

      LOGICAL                   :: old_and_wrong

      COMPLEX                   :: wtfq

      TYPE(t_tlmplm) :: td, tdV1
      TYPE(t_potden) :: vx
      TYPE(t_hub1data) :: hub1data
      TYPE(t_usdus)             :: ud
      TYPE(t_lapw)              :: lapw, lapwq
      CLASS(t_mat), ALLOCATABLE :: zMatk, zMatq, zMat1, zMat2
      CLASS(t_mat), ALLOCATABLE :: hmat,smat

      INTEGER                   :: dealloc_stat, nbasfcnq, nbasfcn, neigk, neigq, noccbd, noccbdq, noccbdmin
      character(len=300)        :: errmsg
      INTEGER, ALLOCATABLE      :: ev_list(:), q_ev_list(:)
      COMPLEX, ALLOCATABLE      :: tempVec(:), tempMat1(:), tempMat2(:), z1H(:,:), z1S(:,:), tempMat3(:), z1H2(:,:), z1S2(:,:)
      REAL,    ALLOCATABLE      :: eigk(:), eigq(:), eigs1(:), eigBuffer(:,:,:)

      COMPLEX  zdotc
      EXTERNAL zdotc


      old_and_wrong = .FALSE.

      CALL vx%copyPotDen(vTot)
      ALLOCATE(vx%pw_w, mold=vx%pw)
      vx%pw_w = vTot%pw_w

      ! Get the (lm) matrix elements for V1 and H0
      CALL ud%init(fi%atoms,fi%input%jspins)
      CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,.FALSE.)
      CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,vTot,vx,inden,fi%input,fi%hub1inp,hub1data,td,ud,0.0,.TRUE.)
      
      ALLOCATE(eigBuffer(fi%input%neig,fi%kpts%nkpt,fi%input%jspins))
      eigBuffer = 0.0
      results1%eig = 1.0e300

      DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
         k_loop:DO nk_i = 1,size(fmpi%k_list)
            nk=fmpi%k_list(nk_i)

            ! Get the required eigenvectors and values at k for occupied bands:
            bkpt = fi%kpts%bk(:, nk)

            q_loop = bqpt

            CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)
            CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi, q_loop)

            noccbd  = COUNT(results%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
            noccbdq = COUNT(resultsq%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)

            noccbdmin = MIN(noccbdq,noccbd)

            nbasfcn  = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
            nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

            IF (fmpi%n_size == 1) THEN
               ALLOCATE (t_mat::zMatk)
               ALLOCATE (t_mat::zMatq)
            ELSE
               ALLOCATE (t_mpimat::zMatk)
               ALLOCATE (t_mpimat::zMatq)
            END IF

            ! Initialize the expansion coefficient matrices at k and k+q
            ! Then read all the stuff into it
            CALL zMatk%init(l_real,nbasfcn,noccbd)
            CALL zMatq%init(l_real,nbasfcnq,nbasfcnq)

            ALLOCATE(ev_list(noccbd))
            ev_list = (/(i, i=1,noccbd, 1)/)
            ALLOCATE(q_ev_list(nbasfcnq))
            q_ev_list = (/(i, i=1,nbasfcnq, 1)/)

            ALLOCATE(eigk(noccbd))
            ALLOCATE(eigq(nbasfcnq))
            ALLOCATE(eigs1(noccbd))

            CALL timestart("Read eigenstuff at k/k+q")
            CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=zMatk)
            CALL read_eig(q_eig_id, nk, jsp, list=q_ev_list, neig=neigq, eig=eigq, zmat=zMatq)
            CALL timestop("Read eigenstuff at k/k+q")

            ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
            CALL timestart("Setup of matrix perturbations")
            CALL dfpt_eigen_hssetup(jsp,fmpi,fi,enpara,nococonv,starsq,ud,td,tdV1,vTot,v1real,lapw,lapwq,iDir,iDtype,hmat,smat,nk,killcont)
            CALL timestop("Setup of matrix perturbations")

            IF (fmpi%n_size == 1) THEN
               ALLOCATE (t_mat::zMat1)
            ELSE
               ALLOCATE (t_mpimat::zMat1)
            END IF

            ! Initialize the expansion coefficient perturbation matrix
            CALL zMat1%init(.FALSE.,nbasfcnq,noccbd)

            ALLOCATE(z1H,mold=zMat1%data_c)
            ALLOCATE(z1S,mold=zMat1%data_c)
            z1H = CMPLX(0.0,0.0)
            z1S = CMPLX(0.0,0.0)

            ! For the dynmat calculation: initialize the auxiliary coefficients
            IF (.NOT.sh_den.AND..NOT.old_and_wrong) THEN
               IF (fmpi%n_size == 1) THEN
                  ALLOCATE (t_mat::zMat2)
               ELSE
                  ALLOCATE (t_mpimat::zMat2)
               END IF

               CALL zMat2%init(.FALSE.,nbasfcnq,noccbd)

               ALLOCATE(z1H2,mold=zMat2%data_c)
               ALLOCATE(z1S2,mold=zMat2%data_c)
               z1H2 = CMPLX(0.0,0.0)
               z1S2 = CMPLX(0.0,0.0)
            END IF

            ! Allocate auxiliary quantities
            ALLOCATE(tempVec(nbasfcnq))
            ALLOCATE(tempMat1(nbasfcnq))
            ALLOCATE(tempMat2(nbasfcnq))
            IF (.NOT.sh_den.AND..NOT.old_and_wrong) ALLOCATE(tempMat3(nbasfcnq))

            CALL timestart("Matrix multiplications")
            DO nu = 1, noccbd
               eigs1(nu) = 0.0

               ! TODO: At the moment H1 and S1 are handled seperately; this was
               !       for debugging purposes; possibly revert for efficiency.
               IF (l_real) THEN ! l_real for zMatk
                  tempVec(:nbasfcnq) = MATMUL(hmat%data_c,zMatk%data_r(:nbasfcn,nu))
               ELSE
                  CALL CPP_zgemv('N',nbasfcnq,nbasfcn,CMPLX(1.0,0.0),hmat%data_c,nbasfcnq,zMatk%data_c(:nbasfcn,nu),1,CMPLX(0.0,0.0),tempVec,1)
                  !tempVec(:nbasfcnq) = MATMUL(hmat%data_c,zMatk%data_c(:nbasfcn,nu))
               END IF

               IF (norm2(q_loop).LT.1e-8) THEN
                  IF (nbasfcnq.NE.nbasfcn) CALL juDFT_error("nbasfcnq/=nbasfcn for q=0", calledby="dfpt_eigen.F90")
                  IF (l_real) THEN
                     eigs1(nu) = REAL(DOT_PRODUCT(zMatk%data_r(:nbasfcn,nu),tempVec))
                  ELSE
                     eigs1(nu) = REAL(zdotc(nbasfcn,zMatk%data_c(:nbasfcn,nu),1,tempVec,1))
                     !eigs1(nu) = REAL(DOT_PRODUCT(zMatk%data_c(:nbasfcn,nu),tempVec))
                  END IF
               ELSE
                  eigs1(nu) = 0
               END IF

               IF (zMatq%l_real) THEN ! l_real for zMatq
                  tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
               ELSE
                  CALL CPP_zgemv('C',nbasfcnq,nbasfcnq,CMPLX(1.0,0.0),zmatq%data_c,nbasfcnq,tempvec,1,CMPLX(0.0,0.0),tempMat1,1)
                  !tempMat1(:nbasfcnq) = MATMUL(CONJG(TRANSPOSE(zMatq%data_c)),tempvec)
               END IF

               ! tempMat1 = H^{(1}_{\nu'\nu}
               DO iNupr = 1, nbasfcnq
                  IF (.NOT.sh_den.AND.old_and_wrong) THEN
                     IF (norm2(bqpt)<1e-8.AND.iNupr==nu) THEN
                        tempMat2(iNupr) = 0.0
                     ELSE IF (ABS(eigq(iNupr)-eigk(nu))<fi%juPhon%eDiffCut) THEN
                        tempMat2(iNupr) = 0.0
                     ELSE
                        tempMat2(iNupr) = 1.0/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr)
                     END IF
                  ELSE IF (.NOT.sh_den.AND..NOT.old_and_wrong) THEN
                     IF (norm2(bqpt)<1e-8.AND.iNupr==nu) THEN
                        tempMat2(iNupr) = 0.0
                        tempMat3(iNupr) = 0.0
                     ELSE IF (ABS(eigq(iNupr)-eigk(nu))<fi%juPhon%eDiffCut) THEN
                        tempMat2(iNupr) = 0.0
                        ! Additional correction term that constitutes new
                        ! coefficients:
                        tempMat3(iNupr) = 0.5 * tempMat1(iNupr)
                     ! TODO: This part of the correction had no effect whatsoever yet.
                     !       Reactivate for misbehaving materials and see if there are
                     !       changes.
                     ELSE IF (iNuPr<=noccbdmin.AND.nu<=noccbdmin) THEN
                        wtfq = resultsq%w_iks(iNupr,nk,jsp)/fi%kpts%wtkpt(nk)
                        tempMat2(iNupr) = 1.0/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr) &
                                      & *(1.0-wtfq)
                        ! Additional correction term that constitutes new
                        ! coefficients:
                        tempMat3(iNupr) = 0.5 * tempMat1(iNupr) * wtfq
                     ELSE
                        tempMat2(iNupr) = 1.0/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr)
                        tempMat3(iNupr) = 0.0
                     END IF
                  ELSE
                     IF (norm2(bqpt)<1e-8.AND.iNupr==nu) THEN
                        tempMat2(iNupr) = 0.0
                     ELSE IF (ABS(eigq(iNupr)-eigk(nu))<fi%juPhon%eDiffCut) THEN
                        tempMat2(iNupr) = 0.0
                     ELSE IF (iNuPr<=noccbdmin.AND.nu<=noccbdmin) THEN
                        wtfq = resultsq%w_iks(iNupr,nk,jsp)/fi%kpts%wtkpt(nk)
                        tempMat2(iNupr) = 1.0/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr) &
                                       & *(1.0-wtfq)
                     ELSE
                        tempMat2(iNupr) = 1.0/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr)
                     END IF
                  END IF
               END DO

               IF (zMatq%l_real) THEN
                  z1H(:nbasfcnq,nu) = -MATMUL(zMatq%data_r,tempMat2(:nbasfcnq))
                  IF (.NOT.sh_den.AND..NOT.old_and_wrong) z1H2(:nbasfcnq,nu) = -MATMUL(zMatq%data_r,tempMat3(:nbasfcnq))
               ELSE
                  CALL CPP_zgemv('N',nbasfcnq,nbasfcnq,CMPLX(-1.0,0.0),zMatq%data_c,nbasfcnq,tempMat2,1,CMPLX(0.0,0.0),z1H(:nbasfcnq,nu),1)
                  !z1H(:nbasfcnq,nu) = -MATMUL(zMatq%data_c,tempMat2(:nbasfcnq))
                  IF (.NOT.sh_den.AND..NOT.old_and_wrong) CALL CPP_zgemv('N',nbasfcnq,nbasfcnq,CMPLX(-1.0,0.0),zmatq%data_c,nbasfcnq,tempMat3,1,CMPLX(0.0,0.0),z1H2(:nbasfcnq,nu),1)
                  !IF (.NOT.sh_den.AND..NOT.old_and_wrong) z1H2(:nbasfcnq,nu) = -MATMUL(zMatq%data_c,tempMat3(:nbasfcnq))
               END IF

               IF (l_real) THEN ! l_real for zMatk
                  tempVec(:nbasfcnq) = MATMUL(smat%data_c,zMatk%data_r(:nbasfcn,nu))
               ELSE
                  CALL CPP_zgemv('N',nbasfcnq,nbasfcn,CMPLX(1.0,0.0),smat%data_c,nbasfcnq,zMatk%data_c(:nbasfcn,nu),1,CMPLX(0.0,0.0),tempVec,1)
                  !tempVec(:nbasfcnq) = MATMUL(smat%data_c,zMatk%data_c(:nbasfcn,nu))
               END IF

               IF (norm2(q_loop).LT.1e-8) THEN
                  IF (nbasfcnq.NE.nbasfcn) CALL juDFT_error("nbasfcnq/=nbasfcn for q=0", calledby="dfpt_eigen.F90")
                  IF (l_real) THEN
                     eigs1(nu) = eigs1(nu) - eigk(nu)*REAL(DOT_PRODUCT(zMatk%data_r(:nbasfcn,nu),tempVec))
                  ELSE
                     eigs1(nu) = eigs1(nu) - eigk(nu)*REAL(zdotc(nbasfcn,zMatk%data_c(:nbasfcn,nu),1,tempVec,1))
                     !eigs1(nu) = eigs1(nu) - eigk(nu)*REAL(DOT_PRODUCT(zMatk%data_c(:nbasfcn,nu),tempVec))
                  END IF
               ELSE
                  eigs1(nu) = 0
               END IF

               IF (zMatq%l_real) THEN ! l_real for zMatq
                  tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
               ELSE
                  CALL CPP_zgemv('C',nbasfcnq,nbasfcnq,CMPLX(1.0,0.0),zmatq%data_c,nbasfcnq,tempvec,1,CMPLX(0.0,0.0),tempMat1,1)
                  !tempMat1(:nbasfcnq) = MATMUL(CONJG(TRANSPOSE(zMatq%data_c)),tempvec)
               END IF
                  
               ! tempMat1 = S^{(1}_{\nu'\nu}
               DO iNupr = 1, nbasfcnq
                  IF (.NOT.sh_den.AND.old_and_wrong) THEN
                     IF (norm2(bqpt)<1e-8.AND.iNupr==nu) THEN
                        tempMat2(iNupr) = 0.0
                     ELSE IF (ABS(eigq(iNupr)-eigk(nu))<fi%juPhon%eDiffCut) THEN
                        tempMat2(iNupr) = 0.5*tempMat1(iNupr)
                     ELSE
                        tempMat2(iNupr) = -eigk(nu)/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr)
                     END IF
                  ELSE IF (.NOT.sh_den.AND..NOT.old_and_wrong) THEN
                     IF (norm2(bqpt)<1e-8.AND.iNupr==nu) THEN
                        tempMat2(iNupr) = 0.0
                        tempMat3(iNupr) = 0.0
                     ELSE IF (ABS(eigq(iNupr)-eigk(nu))<fi%juPhon%eDiffCut) THEN
                        tempMat2(iNupr) = 0.5 * tempMat1(iNupr)
                        ! Additional correction term that constitutes new
                        ! coefficients:
                        tempMat3(iNupr) = -0.5 * eigk(nu) * tempMat1(iNupr)
                     ! TODO: This part of the correction had no effect whatsoever yet.
                     !       Reactivate for misbehaving materials and see if there are
                     !       changes.
                     ELSE IF (iNuPr<=noccbdmin.AND.nu<=noccbdmin) THEN
                        wtfq = resultsq%w_iks(iNupr,nk,jsp)/fi%kpts%wtkpt(nk)
                        tempMat2(iNupr) = -eigk(nu)/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr) &
                                      & *(1.0-wtfq) &
                                      & +0.5*tempMat1(iNupr)*wtfq
                        ! Additional correction term that constitutes new
                        ! coefficients:
                        tempMat3(iNupr) = -0.5 * eigq(iNupr) * tempMat1(iNupr) * wtfq
                     ELSE
                        tempMat2(iNupr) = -eigk(nu)/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr)
                        tempMat3(iNupr) = 0.0
                     END IF
                  ELSE
                     IF (norm2(bqpt)<1e-8.AND.iNupr==nu) THEN
                        tempMat2(iNupr) = 0.0
                     ELSE IF (ABS(eigq(iNupr)-eigk(nu))<fi%juPhon%eDiffCut) THEN
                        tempMat2(iNupr) = 0.5*tempMat1(iNupr)
                     ELSE IF (iNuPr<=noccbdmin.AND.nu<=noccbdmin) THEN
                        wtfq = resultsq%w_iks(iNupr,nk,jsp)/fi%kpts%wtkpt(nk)

                        tempMat2(iNupr) = -eigk(nu)/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr) &
                                       & *(1.0-wtfq) &
                                       & +0.5*tempMat1(iNupr)*wtfq
                     ELSE
                        tempMat2(iNupr) = -eigk(nu)/(eigq(iNupr)-eigk(nu))*tempMat1(iNupr)
                     END IF
                  END IF
               END DO

               IF (zMatq%l_real) THEN
                  z1S(:nbasfcnq,nu) = -MATMUL(zMatq%data_r,tempMat2(:nbasfcnq))
                  IF (.NOT.sh_den.AND..NOT.old_and_wrong) z1S2(:nbasfcnq,nu) = -MATMUL(zMatq%data_r,tempMat3(:nbasfcnq))
               ELSE
                  CALL CPP_zgemv('N',nbasfcnq,nbasfcnq,CMPLX(-1.0,0.0),zmatq%data_c,nbasfcnq,tempMat2,1,CMPLX(0.0,0.0),z1S(:nbasfcnq,nu),1)
                  !z1S(:nbasfcnq,nu) = -MATMUL(zMatq%data_c,tempMat2(:nbasfcnq))
                  IF (.NOT.sh_den.AND..NOT.old_and_wrong) CALL CPP_zgemv('N',nbasfcnq,nbasfcnq,CMPLX(-1.0,0.0),zmatq%data_c,nbasfcnq,tempMat3,1,CMPLX(0.0,0.0),z1S2(:nbasfcnq,nu),1)
                  !IF (.NOT.sh_den.AND..NOT.old_and_wrong) z1S2(:nbasfcnq,nu) = -MATMUL(zMatq%data_c,tempMat3(:nbasfcnq))
               END IF

               zMat1%data_c(:nbasfcnq,nu) = z1H(:nbasfcnq,nu) + z1S(:nbasfcnq,nu)
               IF (.NOT.sh_den.AND..NOT.old_and_wrong) zMat2%data_c(:nbasfcnq,nu) = z1H2(:nbasfcnq,nu) + z1S2(:nbasfcnq,nu)
            END DO

            results1%neig = results%neig

            CALL timestop("Matrix multiplications")

            CALL smat%free()
            CALL hmat%free()
            DEALLOCATE(hmat,smat, stat=dealloc_stat, errmsg=errmsg)
            IF(dealloc_stat /= 0) CALL juDFT_error("Deallocation failed for hmat or smat", hint=errmsg, calledby="dfpt_eigen.F90")

!#ifdef CPP_MPI
                     !CALL MPI_BARRIER(fmpi%mpi_comm,iErr) ! Synchronizes the RMA operations
!#endif

            ! Output results
            CALL timestart("EV1 output")

            IF (fmpi%n_rank == 0) THEN
#ifdef CPP_MPI
               CALL MPI_COMM_RANK(fmpi%diag_sub_comm,n_rank,err)
               CALL MPI_COMM_SIZE(fmpi%diag_sub_comm,n_size,err)
#else
               n_rank = 0; n_size=1;
#endif
               CALL write_eig(dfpt_eig_id, nk, jsp, noccbd, noccbd, &
                              eigs1(:noccbd), n_start=n_size,n_end=n_rank,zMat=zMat1)
               eigBuffer(:noccbd,nk,jsp) = eigs1(:noccbd)
               IF (.NOT.sh_den.AND..NOT.old_and_wrong) CALL write_eig(dfpt_eig_id2, nk, jsp, noccbd, noccbd, &
                                                                      eigs1(:noccbd), n_start=n_size,n_end=n_rank,zMat=zMat2)
            ELSE
               IF (fmpi%pe_diag) CALL write_eig(dfpt_eig_id, nk, jsp, noccbd, &
                              n_start=fmpi%n_size,n_end=fmpi%n_rank,zMat=zMat1)
               IF ((.NOT.sh_den).AND.(.NOT.old_and_wrong).AND.fmpi%pe_diag) CALL write_eig(dfpt_eig_id2, nk, jsp, noccbd, &
                                                                              n_start=fmpi%n_size,n_end=fmpi%n_rank,zMat=zMat2)
            END IF

#if defined(CPP_MPI)
            ! RMA synchronization
            CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
            CALL timestop("EV1 output")

            IF (ALLOCATED(ev_list)) DEALLOCATE(ev_list)
            IF (ALLOCATED(q_ev_list)) DEALLOCATE(q_ev_list)
            IF (ALLOCATED(eigk)) DEALLOCATE(eigk)
            IF (ALLOCATED(eigq)) DEALLOCATE(eigq)
            IF (ALLOCATED(eigs1)) DEALLOCATE(eigs1)
            IF (ALLOCATED(tempVec)) DEALLOCATE(tempVec)
            IF (ALLOCATED(tempMat1)) DEALLOCATE(tempMat1)
            IF (ALLOCATED(tempMat2)) DEALLOCATE(tempMat2)
            IF (ALLOCATED(tempMat3)) DEALLOCATE(tempMat3)
            IF (ALLOCATED(z1H)) DEALLOCATE(z1H)
            IF (ALLOCATED(z1S)) DEALLOCATE(z1S)
            IF (ALLOCATED(z1H2)) DEALLOCATE(z1H2)
            IF (ALLOCATED(z1S2)) DEALLOCATE(z1S2)
            IF (ALLOCATED(zmatk)) THEN
               CALL zMatk%free()
               DEALLOCATE(zMatk)
            END IF
            IF (ALLOCATED(zmatq)) THEN
               CALL zMatq%free()
               DEALLOCATE(zMatq)
            END IF
            IF (ALLOCATED(zmat1)) THEN
               CALL zMat1%free()
               DEALLOCATE(zMat1)
            END IF
            IF (ALLOCATED(zmat2)) THEN
               CALL zMat2%free()
               DEALLOCATE(zMat2)
            END IF
         END DO  k_loop
      END DO ! spin loop ends
      neigd2 = MIN(fi%input%neig,lapw%dim_nbasfcn())
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(eigBuffer(:neigd2,:,:),results1%eig(:neigd2,:,:),neigd2*fi%kpts%nkpt*fi%input%jspins,MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#else
      results1%eig(:neigd2,:,:) = eigBuffer(:neigd2,:,:)
#endif

   END SUBROUTINE dfpt_eigen
END MODULE m_dfpt_eigen
