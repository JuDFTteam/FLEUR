!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_eigen_new

#ifdef CPP_MPI
   USE mpi
#endif
   USE m_juDFT

   IMPLICIT NONE

CONTAINS

   SUBROUTINE dfpt_eigen_new(fi, sphhar, results, resultsq, fmpi, enpara, nococonv, starsq, v1real, v1imag, vTot, inden, bqpt, &
                             eig_id, q_eig_id, dfpt_eig_id, iDir, iDtype, killcont, l_real, dfpt_tag)

      USE m_types
      USE m_constants
      USE m_dfpt_eigen_hssetup
      USE m_pot_io
      USE m_util
      USE m_eig66_io, ONLY : write_eig, read_eig
      USE m_xmlOutput
      USE m_types_mpimat
      USE m_invert_HepsS
      USE m_dfpt_tlmplm
      USE m_mt_setup
      !USE m_npy

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_sphhar),     INTENT(IN)       :: sphhar
      TYPE(t_results),INTENT(INOUT):: results, resultsq
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_enpara),INTENT(IN) :: enpara
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      TYPE(t_stars),INTENT(IN)     :: starsq
      TYPE(t_potden),INTENT(IN)    :: inden, v1real, v1imag, vTot
      REAL,         INTENT(IN)     :: bqpt(3)
      INTEGER,      INTENT(IN)     :: eig_id, q_eig_id, dfpt_eig_id, iDir, iDtype, killcont(6)
      LOGICAL,      INTENT(IN)     :: l_real
      CHARACTER(len=20), INTENT(IN) :: dfpt_tag

      INTEGER n_size,n_rank
      INTEGER i,err,nk,jsp,nk_i

      INTEGER              :: ierr

      REAL :: bkpt(3)

      INTEGER                   :: nu

      TYPE(t_tlmplm) :: td, tdV1
      TYPE(t_potden) :: vx
      TYPE(t_hub1data) :: hub1data
      TYPE(t_usdus)             :: ud
      TYPE(t_lapw)              :: lapw, lapwq
      TYPE(t_kpts)              :: kqpts ! basically kpts, but with q added onto each one.
      CLASS(t_mat), ALLOCATABLE :: zMatk, zMatq, zMat1
      CLASS(t_mat), ALLOCATABLE :: hmat,smat

      INTEGER                   :: dealloc_stat, nbasfcnq, nbasfcn, neigk, neigq, noccbd, noccbdq
      character(len=300)        :: errmsg
      INTEGER, ALLOCATABLE      :: ev_list(:), q_ev_list(:), k_selection(:)
      COMPLEX, ALLOCATABLE      :: tempVec(:), tempMat1(:), tempMat2(:)
      REAL,    ALLOCATABLE      :: eigk(:), eigq(:), eigs1(:)

      CLASS(t_mat), ALLOCATABLE :: invE(:), matE(:)

      CALL timestart("dfpt_eigen")

      ALLOCATE(k_selection(7))
      k_selection = [1,45,77,255,366,412,512] 

      CALL vx%copyPotDen(vTot)
      ALLOCATE(vx%pw_w, mold=vx%pw)
      vx%pw_w = vTot%pw_w

      call ud%init(fi%atoms,fi%input%jspins)
      kqpts = fi%kpts
      ! Modify this from kpts only in DFPT case.
      DO nk_i = 1, fi%kpts%nkpt
         kqpts%bk(:, nk_i) = kqpts%bk(:, nk_i) + bqpt
      END DO

      CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,.FALSE.)
      CALL mt_setup(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,nococonv,enpara,fi%hub1inp,hub1data,inden,vTot,vx,fmpi,td,ud,0.0,.TRUE.)

      DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
         k_loop:DO nk_i = 1,size(fmpi%k_list)
               nk=fmpi%k_list(nk_i)

               ! Get the required eigenvectors and values at k for occupied bands:
               bkpt = fi%kpts%bk(:, nk)

               CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)
               CALL lapwq%init(fi%input, fi%noco, nococonv, kqpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)

               noccbd  = COUNT(results%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
               noccbdq = COUNT(resultsq%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)

               nbasfcn  = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
               nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

               IF (fmpi%n_size == 1) THEN
                  ALLOCATE (t_mat::zMatk)
                  ALLOCATE (t_mat::zMatq)
               ELSE
                  ALLOCATE (t_mpimat::zMatk)
                  ALLOCATE (t_mpimat::zMatq)
               END IF
               CALL zMatk%init(l_real,nbasfcn,noccbd)
               CALL zMatq%init(l_real,nbasfcnq,nbasfcnq)
               ALLOCATE(ev_list(noccbd))
               ev_list = (/(i, i=1,noccbd, 1)/)
               ALLOCATE(q_ev_list(noccbdq))
               q_ev_list = (/(i, i=1,nbasfcnq, 1)/)
               ALLOCATE(eigk(noccbd))
               ALLOCATE(eigq(nbasfcnq))
               ALLOCATE(eigs1(noccbd))

               CALL timestart("Read eigenstuff at k/k+q")
               CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=zMatk)
               CALL read_eig(q_eig_id, nk, jsp, list=q_ev_list, neig=neigq, eig=eigq, zmat=zMatq)
               CALL timestop("Read eigenstuff at k/k+q")

               CALL timestart("Energy inversion")
               CALL invert_HepsS(fmpi, fi%atoms, fi%noco, fi%juPhon, lapwq, zMatq, eigq, eigk, neigq, noccbd, zMatq%l_real, invE, matE)
               CALL timestop("Energy inversion")

               ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
               CALL timestart("Setup of matrix perturbations")
               CALL dfpt_eigen_hssetup(jsp,fmpi,fi,enpara,nococonv,starsq,ud,td,tdV1,v1real,lapw,lapwq,iDir,iDtype,smat,hmat,nk,killcont)
               CALL timestop("Setup of matrix perturbations")

               IF (fmpi%n_size == 1) THEN
                  ALLOCATE (t_mat::zMat1)
               ELSE
                  ALLOCATE (t_mpimat::zMat1)
               END IF
               CALL zMat1%init(.FALSE.,nbasfcnq,noccbd)

               ALLOCATE(tempVec(nbasfcnq))
               ALLOCATE(tempMat1(nbasfcnq))
               ALLOCATE(tempMat2(neigq))

               !TODO: Optimize this with (SCA)LAPACK CALLS
               DO nu = 1, noccbd
                  IF (l_real) THEN ! l_real for zMatk
                     tempVec(:nbasfcnq) = MATMUL(hmat%data_c-eigk(nu)*smat%data_c,zMatk%data_r(:nbasfcn,nu))
                  ELSE
                     tempVec(:nbasfcnq) = MATMUL(hmat%data_c-eigk(nu)*smat%data_c,zMatk%data_c(:nbasfcn,nu))
                  END IF

                  IF (zMatq%l_real) THEN ! l_real for zMatq
                     tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                  ELSE
                     tempMat1(:nbasfcnq) = MATMUL(CONJG(TRANSPOSE(zMatq%data_c)),tempvec)
                  END IF

                  tempMat2(:neigq) = MATMUL(invE(nu)%data_r,tempMat1)

                  IF (norm2(bqpt).LT.1e-8) THEN
                     IF (nbasfcnq.NE.nbasfcn) CALL juDFT_error("nbasfcnq/=nbasfcn for q=0", calledby="dfpt_eigen.F90")
                     IF (l_real) THEN
                        eigs1 = DOT_PRODUCT(zMatk%data_r(:nbasfcn,nu),tempVec)
                     ELSE
                        eigs1 = DOT_PRODUCT(zMatk%data_c(:nbasfcn,nu),tempVec) !real(?)
                     END IF
                  ELSE
                     eigs1 = 0
                  END IF

                  IF (zMatq%l_real) THEN
                     zMat1%data_c(:nbasfcnq,nu) = -MATMUL(zMatq%data_r,tempMat2(:neigq))
                  ELSE
                     zMat1%data_c(:nbasfcnq,nu) = -MATMUL(zMatq%data_c,tempMat2(:neigq))
                  END IF

                  IF (ANY(nk==k_selection)) THEN
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_"//int2str(nu)//"_tempVec.npy",tempVec)
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_"//int2str(nu)//"_HS1band.npy",tempMat1)
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_"//int2str(nu)//"_matE.npy",matE(nu)%data_r)
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_"//int2str(nu)//"_invE.npy",invE(nu)%data_r)
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_"//int2str(nu)//"_z1band.npy",tempMat2)
                  END IF
               END DO

               IF (ANY(nk==k_selection)) THEN
                  IF (l_real) THEN ! l_real for zMatk
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_zMatk.npy",zMatk%data_r)
                  ELSE
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_zMatk.npy",zMatk%data_c)
                  END IF

                  IF (zMatq%l_real) THEN ! l_real for zMatq
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_zMatkq.npy",zMatq%data_r)
                  ELSE
                     CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_zMatkq.npy",zMatq%data_c)
                  END IF

                  CALL save_npy(TRIM(dfpt_tag)//"_"//int2str(nk)//"_z1.npy",zMat1%data_c)
               END IF

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
                  ELSE
                     IF (fmpi%pe_diag) CALL write_eig(dfpt_eig_id, nk, jsp, noccbd, &
                                    n_start=fmpi%n_size,n_end=fmpi%n_rank,zMat=zMat1)
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

          END DO  k_loop
        END DO ! spin loop ends

        CALL timestop("dfpt_eigen")

   END SUBROUTINE dfpt_eigen_new
END MODULE m_dfpt_eigen_new
