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

   IMPLICIT NONE

CONTAINS

   SUBROUTINE dfpt_eigen(fi, jsp, nk, results, fmpi, enpara, nococonv, starsq, v1, lapwq, td, tdV1, ud, zMatq, eigq, bqpt, neigq, eig_id, dfpt_eig_id, iDir, iDtype, killcont, l_real)

      USE m_types
      USE m_constants
      USE m_dfpt_eigen_hssetup
      USE m_pot_io
      USE m_util
      USE m_eig66_io, ONLY : write_eig, read_eig
      USE m_xmlOutput
      USE m_types_mpimat
      USE m_invert_HepsS
      USE m_npy

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_results),INTENT(IN):: results
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_enpara),INTENT(IN) :: enpara
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      TYPE(t_stars),INTENT(IN)     :: starsq
      TYPE(t_potden),INTENT(IN)    :: v1
      TYPE(t_lapw)  ,INTENT(IN)    :: lapwq
      TYPE(t_tlmplm) ,INTENT(IN)           :: td, tdV1
      TYPE(t_usdus)  ,INTENT(IN)           :: ud
      CLASS(t_mat), INTENT(IN)     :: zMatq
      REAL,         INTENT(IN)     :: eigq(:), bqpt(3)
      INTEGER,      INTENT(IN)     :: neigq, eig_id, dfpt_eig_id, iDir, iDtype, nk, jsp, killcont(6)
      LOGICAL,      INTENT(IN)     :: l_real

      INTEGER n_size,n_rank
      INTEGER i,err

      INTEGER              :: ierr

      REAL :: bkpt(3)

      INTEGER                   :: nu

      TYPE(t_lapw)              :: lapw
      CLASS(t_mat), ALLOCATABLE :: zMatk, zMat1
      CLASS(t_mat), ALLOCATABLE :: hmat,smat

      INTEGER                   :: dealloc_stat, nbasfcnq, nbasfcn, neigk, noccbd
      character(len=300)        :: errmsg
      INTEGER, ALLOCATABLE      :: ev_list(:)
      COMPLEX, ALLOCATABLE      :: tempVec(:), tempMat1(:), tempMat2(:)
      REAL,    ALLOCATABLE      :: eigk(:), eigs1(:)

      CLASS(t_mat), ALLOCATABLE :: invHepsS(:), invE(:)

      CALL timestart("dfpt_eigen")

      ! Get the required eigenvectors and values at k for occupied bands:
      bkpt = fi%kpts%bk(:, nk)

      CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)

      noccbd = COUNT(results%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
      nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::zMatk)
      ELSE
         ALLOCATE (t_mpimat::zMatk)
      END IF
      CALL zMatk%init(l_real,nbasfcn,noccbd)
      ALLOCATE(ev_list(noccbd))
      ev_list = (/(i, i=1,noccbd, 1)/)
      ALLOCATE(eigk(noccbd))
      ALLOCATE(eigs1(noccbd))

      CALL timestart("Read eigenstuff at k")
      CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=zMatk)
      CALL timestop("Read eigenstuff at k")
      CALL invert_HepsS(fmpi, fi%atoms, fi%noco, fi%juPhon, lapwq, zMatq, eigq, eigk, neigq, noccbd, zMatq%l_real, invHepsS, invE)

      ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
      CALL timestart("Setup of matrix perturbations")
      CALL dfpt_eigen_hssetup(jsp,fmpi,fi,enpara,nococonv,starsq,ud,td,tdV1,v1,lapw,lapwq,iDir,iDtype,smat,hmat,nk,killcont)
      CALL timestop("Setup of matrix perturbations")

      nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
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
            IF (nk==1) CALL save_npy("new_zKet_"//int2str(nk)//".npy",zMatk%data_r)
         ELSE
            tempVec(:nbasfcnq) = MATMUL(hmat%data_c-eigk(nu)*smat%data_c,zMatk%data_c(:nbasfcn,nu))
            IF (nk==1) CALL save_npy("new_zKet_"//int2str(nk)//".npy",zMatk%data_c)
         END IF

         IF (zMatq%l_real) THEN ! l_real for zMatq
            tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
            IF (nk==1) CALL save_npy("new_zBra_"//int2str(nk)//".npy",zMatq%data_r)
         ELSE
            tempMat1(:nbasfcnq) = MATMUL(CONJG(TRANSPOSE(zMatq%data_c)),tempvec)
            IF (nk==1) CALL save_npy("new_zBra_"//int2str(nk)//".npy",zMatq%data_c)
         END IF
         IF (nk==1) CALL save_npy("new_hepss1band_"//int2str(nk)//"_"//int2str(nu)//".npy", tempMat1)

         tempMat2(:neigq) = MATMUL(invE(nu)%data_r,tempMat1)
         IF (nk==1) CALL save_npy("new_z1band_"//int2str(nk)//"_"//int2str(nu)//".npy",tempMat2)

         ! TODO: Reactivate this!
         !IF (norm2(bqpt).LT.1e-8) THEN
         !   IF (nbasfcnq.NE.nbasfcn) CALL juDFT_error("nbasfcnq/=nbasfcn for q=0", calledby="dfpt_eigen.F90")
         !   IF (l_real) THEN
         !      eigs1 = DOT_PRODUCT(zMatk%data_r(:nbasfcn,nu),tempVec)
         !   ELSE
         !      eigs1 = DOT_PRODUCT(zMatk%data_c(:nbasfcn,nu),tempVec) !real(?)
         !   END IF
         !ELSE
            eigs1 = 0
         !END IF
         !IF (nk==1) CALL save_npy("new_eig1_"//int2str(nk)//".npy",eigs1)
         IF (zMatq%l_real) THEN
            zMat1%data_c(:nbasfcnq,nu) = -MATMUL(zMatq%data_r,tempMat2(:neigq))
         ELSE
            zMat1%data_c(:nbasfcnq,nu) = -MATMUL(zMatq%data_c,tempMat2(:neigq))
         END IF
      END DO

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

        IF (ALLOCATED(zmatk)) THEN
          CALL zMatk%free()
          DEALLOCATE(zMatk)
        END IF
        IF (ALLOCATED(zmat1)) THEN
          CALL zMat1%free()
          DEALLOCATE(zMat1)
        END IF
        !TODO: Is this ok?
        !IF (ALLOCATED(invHepsS(1))) THEN
          !DO nu = 1, noccbd
             !CALL invHepsS(nu)%free()
          !END DO
          !DEALLOCATE(invHepsS)
        !END IF

        CALL timestop("dfpt_eigen")

   END SUBROUTINE dfpt_eigen
END MODULE m_dfpt_eigen
