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

    SUBROUTINE dfpt_eigen()

#include"cpp_double.h"
      USE m_types
      USE m_constants
      USE m_eigen_hssetup
      USE m_pot_io
      USE m_eigen_diag
      USE m_mt_setup
      USE m_util
      USE m_eig66_io, ONLY : write_eig, read_eig
      USE m_xmlOutput
      USE m_types_mpimat
      use m_store_load_hybrid
      USE m_invert_HepsS

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_kpts), INTENT(IN) :: kqpts ! basically kpts, but with q added onto each one.
      TYPE(t_results),INTENT(INOUT):: results
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_mpi),INTENT(IN)       :: fmpi

      TYPE(t_mpdata), intent(inout):: mpdata
      TYPE(t_hybdat), INTENT(INOUT):: hybdat
      TYPE(t_enpara),INTENT(IN) :: enpara
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_potden),INTENT(IN)    :: v
      TYPE(t_lapw)  ,INTENT(IN)    :: lapwq
      CLASS(t_mat), INTENT(IN)     :: zMatq
      REAL,         INTENT(IN)     :: eigq(:)
      INTEGER,      INTENT(IN)     :: neigq

      ! Scalar Arguments
      INTEGER,INTENT(IN)    :: iter
      INTEGER,INTENT(IN)    :: eig_id

      TYPE(t_kpts) :: kqpts ! basically kpts, but with q added onto each one.

      ! Local Scalars
      INTEGER jsp,nk,nred,ne_all,ne_found,neigd2
      INTEGER ne, nk_i,n_size,n_rank
      INTEGER isp,i,j,err
      LOGICAL l_file,l_real,l_zref
      INTEGER :: solver=0
      ! Local Arrays
      INTEGER              :: ierr
      INTEGER              :: neigBuffer(fi%kpts%nkpt,fi%input%jspins)

      INTEGER, ALLOCATABLE :: nvBuffer(:,:), nvBufferTemp(:,:), nvfullBuffer(:,:), GbasVecBuffer(:, :, :, :)
      REAL,    ALLOCATABLE :: bkpt(:)
      REAL,    ALLOCATABLE :: eigBuffer(:,:,:)

      INTEGER                   :: jsp_m, i_kpt_m, i_m

      TYPE(t_tlmplm)            :: td
      TYPE(t_usdus)             :: ud
      TYPE(t_lapw)              :: lapw
      CLASS(t_mat), ALLOCATABLE :: zMat
      CLASS(t_mat), ALLOCATABLE :: hmat,smat
      CLASS(t_mat), ALLOCATABLE :: smat_unfold !used for unfolding bandstructure

      INTEGER                   :: dealloc_stat
      character(len=300)        :: errmsg
      INTEGER, ALLOCATABLE      :: ev_list(:)

      l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco).AND.atoms%n_hia==0

      neigBuffer = 0
      eigBuffer = 1.0e300
      nvBuffer = 0
      nvBufferTemp = 0

        nvBuffer(nk,jsp) = lapw%nv(jsp)

        ! Get the required eigenvectors and values at k for occupied bands:
        bkpt = fi%kpts%bk(:, nk)

        CALL lapw%init(input, noco, nococonv, fi%kpts, atoms, sym, nk, cell, .FALSE., fmpi)

        noccbd = COUNT(results%w_iks(:,:,jsp)*2.0/input%jspins>1.e-8)
        nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)

        CALL zMatk%init(l_real,nbasfcn,noccbd)
        ALLOCATE(ev_list(noccbd))
        ev_list = (/(i, i=1,noccbd, 1)/)

        CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=zMatk)

        CALL invert_HepsS(fi%atoms, fi%noco, fi%juPhon, lapwq, zMatq, eigq, eigs, neigq, noccbd, l_real, invHepsS)

        ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
        CALL timestart("Setup of H&S matrices")
        CALL dfpt_eigen_hssetup(jsp,fmpi,fi,mpdata,results,xcpot,enpara,nococonv,stars,sphhar,hybdat,ud,td,v,lapw,nk,smat,hmat,&
                              & kqpts,bqpt,dfpt_eig_id,starsq,tdV1,v1real,lapwq)
        CALL timestop("Setup of H&S matrices")

        nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*atoms%nlotot,lapwq%nv(1)+atoms%nlotot,noco%l_noco)
        CALL z1%init(.FALSE.,nbasfcnq,noccbd)

        ! PSEUDO CODE: This is what we need to calculate eventually.
        !IF (q.EQ.0) THEN
        !    eigs1 = conjg(z0)*(hmat-eigs*smat)*z0
        !END IF

        !z1%data_c = -invHepsS*(hmat-eigs*smat)*z0

        CALL smat%free()
        CALL hmat%free()
        DEALLOCATE(hmat,smat, stat=dealloc_stat, errmsg=errmsg)
        IF(dealloc_stat /= 0) CALL juDFT_error("Deallocation failed for hmat or smat", hint=errmsg, calledby="dfpt_eigen.F90")

!#ifdef CPP_MPI
            !CALL MPI_BARRIER(fmpi%mpi_comm,iErr) ! Synchronizes the RMA operations
!#endif

        ! Output results
        CALL timestart("EV output")

        IF (fmpi%n_rank == 0) THEN
            ! Only process 0 writes out the value of ne_all and the
            ! eigenvalues.
#ifdef CPP_MPI
            call MPI_COMM_RANK(fmpi%diag_sub_comm,n_rank,err)
            call MPI_COMM_SIZE(fmpi%diag_sub_comm,n_size,err)
#else
            n_rank = 0; n_size=1;
#endif

            CALL write_eig(dfpt_eig_id, nk, jsp, noccbd, noccbd, &
                         eigs1(:noccbd), n_start=n_size,n_end=n_rank,zMat=z1)
            eigBuffer(:noccbd, nk, jsp) = eigs1(:noccbd)
        ELSE
            if (fmpi%pe_diag) CALL write_eig(dfpt_eig_id, nk, jsp, noccbd, &
                           n_start=fmpi%n_size,n_end=fmpi%n_rank,zMat=z1)
        ENDIF

#if defined(CPP_MPI)
        ! RMA synchronization
        CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
        CALL timestop("EV output")

        IF (allocated(zmat)) THEN
          call zMat%free()
          deallocate(zMat)
        ENDIF

      neigd2 = MIN(fi%input%neig,lapw%dim_nbasfcn())
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(eigBuffer(:neigd2,:,:),results%eig(:neigd2,:,:),neigd2*fi%kpts%nkpt*fi%input%jspins,MPI_DOUBLE_PRECISION,MPI_MIN,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(nvBuffer(:,:),nvBufferTemp(:,:),size(nvbuffer),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#else
      results%eig(:neigd2,:,:) = eigBuffer(:neigd2,:,:)
      nvBufferTemp(:,:) = nvBuffer(:,:)
#endif

   END SUBROUTINE dfpt_eigen
END MODULE m_dfpt_eigen
