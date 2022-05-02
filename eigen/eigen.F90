!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!<@brief
!<The eigen routine sets up and solves the generalized eigenvalue problem
!More description at end of file
MODULE m_eigen
#ifdef CPP_MPI
   use mpi
#endif
   USE m_juDFT
   IMPLICIT NONE
CONTAINS
   !>The eigenvalue problem is constructed and solved in this routine. The following steps are performed:
   !> 1. Preparation: generate energy parameters, open eig-file
   !> 2. CALL to mt_setup() : this constructs the local Hamiltonian (i.e. the Hamiltonian in the \f$ u,\dot u, u_{lo} \f$ basis) LDA+U is also added here
   !> 3. within the (collinear)spin and k-point loop: CALL to eigen_hssetup() to generate the matrices, CALL to eigen_diag() to perform diagonalization
   !> 4. writing (saving) of eigenvectors
   !>
   !>@author D. Wortmann
   !
   ! Modifications done to use this with DFPT phonons:
   ! a) We need additional MT-integrals from mt_setup that cover the potential variation V1.
   ! b) The eigenvalues are to be evaluated for k+q, not k.
   ! c) Additionally, load in the occupied states for k without q.
   ! d) The work isn't done once the eigenvectors and eigenvalues are found. There is post-
   !    processing to gain the perturbed eigenvalues and eigenvectors.
   ! e) The latter are the actual output for the routine when used for DFPT. They are saved
   !    the same way as the eigenvalues before, but for a shifted eig_id.
   SUBROUTINE eigen(fi,fmpi,stars,sphhar,xcpot,&
                    enpara,nococonv,mpdata,hybdat,&
                    iter,eig_id,results,inden,v,vx,hub1data,nvfull,GbasVec_eig,bqpt,dfpt_eig_id,starsq,v1real,v1imag)

#include"cpp_double.h"
      USE m_types
      USE m_constants
      USE m_eigen_hssetup
      USE m_pot_io
      USE m_eigen_diag
      !USE m_hsefunctional
      USE m_mt_setup
      USE m_util
      !USE m_icorrkeys
      USE m_eig66_io, ONLY : write_eig, read_eig
      USE m_xmlOutput

      USE m_symmetrize_matrix
      USE m_unfold_band_kpts !used for unfolding bands
      USE m_types_mpimat
      use m_store_load_hybrid
      USE m_dfpt_tlmplm

      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      TYPE(t_results),INTENT(INOUT):: results
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_mpi),INTENT(IN)       :: fmpi

      TYPE(t_mpdata), intent(inout):: mpdata
      TYPE(t_hybdat), INTENT(INOUT):: hybdat
      TYPE(t_enpara),INTENT(INOUT) :: enpara
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_potden),INTENT(IN)    :: inden !
      TYPE(t_hub1data),INTENT(INOUT):: hub1data
      TYPE(t_potden), INTENT(IN)   :: vx
      TYPE(t_potden),INTENT(IN)    :: v

!    EXTERNAL MPI_BCAST    !only used by band_unfolding to broadcast the gvec

      ! Scalar Arguments
      INTEGER,INTENT(IN)    :: iter
      INTEGER,INTENT(IN)    :: eig_id

      INTEGER, OPTIONAL, ALLOCATABLE, INTENT(OUT) :: nvfull(:, :), GbasVec_eig(:, :, :, :)

      REAL,    OPTIONAL, INTENT(IN) :: bqpt
      INTEGER, OPTIONAL, INTENT(IN) :: dfpt_eig_id

      TYPE(t_stars),  OPTIONAL, INTENT(IN) :: starsq
      TYPE(t_potden), OPTIONAL, INTENT(IN) :: v1real, v1imag

      ! Local Scalars
      INTEGER jsp,nk,nred,ne_all,ne_found,neigd2
      INTEGER ne, nk_i,n_size,n_rank
      INTEGER isp,i,j,err
      LOGICAL l_real,l_zref
      INTEGER :: solver=0
      ! Local Arrays
      INTEGER              :: ierr
      INTEGER              :: neigBuffer(fi%kpts%nkpt,fi%input%jspins)

      COMPLEX              :: unfoldingBuffer(SIZE(results%unfolding_weights,1),fi%kpts%nkpt,fi%input%jspins) ! needed for unfolding bandstructure fmpi case

      INTEGER, ALLOCATABLE :: nvBuffer(:,:), nvBufferTemp(:,:), nvfullBuffer(:,:), GbasVecBuffer(:, :, :, :)
      REAL,    ALLOCATABLE :: bkpt(:)
      REAL,    ALLOCATABLE :: eig(:), eigBuffer(:,:,:)

      TYPE(t_tlmplm)            :: td, tdV1, tdmod
      TYPE(t_usdus)             :: ud, uddummy
      TYPE(t_lapw)              :: lapw
      CLASS(t_mat), ALLOCATABLE :: zMat
      CLASS(t_mat), ALLOCATABLE :: hmat,smat
      CLASS(t_mat), ALLOCATABLE :: smat_unfold !used for unfolding bandstructure
      TYPE(t_kpts)              :: kqpts ! basically kpts, but with q added onto each one.
      TYPE(t_hub1data)          :: hub1datadummy

      ! Variables for HF or fi%hybinp functional calculation
      INTEGER                   :: comm(fi%kpts%nkpt),irank2(fi%kpts%nkpt),isize2(fi%kpts%nkpt), dealloc_stat
      character(len=300)        :: errmsg
      real                      :: alpha_hybrid

      LOGICAL                   :: l_dfpteigen

      l_dfpteigen = PRESENT(bqpt)

      kqpts = fi%kpts
      ! Modify this from kpts only in DFPT case.
      IF (l_dfpteigen) THEN
          DO nk_i = 1, fi%kpts%nkpt
              kqpts%bk(3, nk_i) = kqpts%bk(3, nk_i) + bqpt
          END DO
      END IF

      call ud%init(fi%atoms,fi%input%jspins)
      call uddummy%init(fi%atoms,fi%input%jspins)
      ALLOCATE(eig(fi%input%neig))
      ALLOCATE(bkpt(3))
      ALLOCATE(eigBuffer(fi%input%neig,fi%kpts%nkpt,fi%input%jspins))
      ALLOCATE(nvBuffer(fi%kpts%nkpt,MERGE(1,fi%input%jspins,fi%noco%l_noco)),nvBufferTemp(fi%kpts%nkpt,MERGE(1,fi%input%jspins,fi%noco%l_noco)))

      ! check if z-reflection trick can be used
      l_zref=(fi%sym%zrfs.AND.(SUM(ABS(kqpts%bk(3,:fi%kpts%nkpt))).LT.1e-9).AND..NOT.fi%noco%l_noco)
      IF (fmpi%n_size > 1) l_zref = .FALSE.

      !IF (fmpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/),(/iter,v%iter/),&
      !                                                RESHAPE((/19,13,5,5/),(/2,2/)))

      ! Set up and solve the eigenvalue problem
      !   loop over spins
      !     set up k-point independent t(l'm',lm) matrices

      alpha_hybrid = MERGE(xcpot%get_exchange_weight(),0.0,hybdat%l_subvxc)
      CALL mt_setup(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,nococonv,enpara,fi%hub1inp,hub1data,inden,v,vx,fmpi,td,ud,alpha_hybrid)
      ! Get matrix elements of perturbed potential and modified H/S in DFPT case.
      IF (l_dfpteigen) THEN
          hub1datadummy = hub1data
          CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,v,fmpi,tdV1,v1real,v1imag)
          CALL mt_setup(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,nococonv,enpara,fi%hub1inp,hub1datadummy,inden,v,vx,fmpi,tdmod,uddummy,alpha_hybrid,.TRUE.)
      END IF

      neigBuffer = 0
      results%neig = 0
      results%eig = 1.0e300
      eigBuffer = 1.0e300
      unfoldingBuffer = CMPLX(0.0,0.0)
      nvBuffer = 0
      nvBufferTemp = 0

      IF (PRESENT(nvfull)) THEN
          ALLOCATE(nvfull(MERGE(1,fi%input%jspins,fi%noco%l_noco), fi%kpts%nkpt))
          nvfull = 0
          ALLOCATE(nvfullBuffer(MERGE(1,fi%input%jspins,fi%noco%l_noco), fi%kpts%nkpt))
          nvfullBuffer = 0
          ALLOCATE(GbasVec_eig(3, fi%input%neig, fi%kpts%nkpt, MERGE(1,fi%input%jspins,fi%noco%l_noco)))
          GbasVec_eig = 0
          ALLOCATE(GbasVecBuffer(3, fi%input%neig, fi%kpts%nkpt, MERGE(1,fi%input%jspins,fi%noco%l_noco)))
          GbasVecBuffer = 0
      END IF

      DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
         k_loop:DO nk_i = 1,size(fmpi%k_list)
            nk=fmpi%k_list(nk_i)
            ! Set up lapw list
            CALL lapw%init(fi%input,fi%noco,nococonv, kqpts, fi%atoms, fi%sym, nk, fi%cell, l_zref, fmpi)

            call timestart("Setup of H&S matrices")
            CALL eigen_hssetup(jsp,fmpi,fi,mpdata,results,vx,xcpot,enpara,nococonv,stars,sphhar,hybdat,ud,td,v,lapw,nk,smat,hmat)
            CALL timestop("Setup of H&S matrices")

            nvBuffer(nk,jsp) = lapw%nv(jsp)

            IF (PRESENT(nvfull)) THEN
                nvfullBuffer(jsp, nk) = lapw%nv(jsp) + fi%atoms%nlotot
                GbasVecBuffer(:, :lapw%nv(jsp), nk, jsp) = lapw%gvec(:, :lapw%nv(jsp), jsp)
            END IF

            ne_all=fi%input%neig
            IF(ne_all < 0) ne_all = lapw%nmat
            IF(ne_all > lapw%nmat) ne_all = lapw%nmat

            !Try to symmetrize matrix
            CALL symmetrize_matrix(fmpi,fi%noco,kqpts,nk,hmat,smat)

            IF (fi%banddos%unfoldband .AND. (.NOT. fi%noco%l_soc)) THEN
               select type(smat)
               type is (t_mat)
                  allocate(t_mat::smat_unfold)
                  select type(smat_unfold)
                  type is (t_mat)
                     smat_unfold=smat
                     smat_unfold%data_c=CONJG(smat%data_c)
                  end select
               type is (t_mpimat)
                  allocate(t_mpimat::smat_unfold)
                  select type(smat_unfold)
                  type is (t_mpimat)
                     smat_unfold=smat
                     smat_unfold%data_c=CONJG(smat%data_c)
                  end select
               end select
            END IF

            ! Solve generalized eigenvalue problem.
            !     ne_all ... number of eigenpairs searched (and found) on this node
            !                on input, overall number of eigenpairs searched,
            !                on output, local number of eigenpairs found
            !     eig ...... all eigenvalues, output
            !     zMat ..... local eigenvectors, output
            if (fmpi%pe_diag) THEN
              CALL eigen_diag(solver,hmat,smat,ne_all,eig,zMat,nk,jsp,iter)
            ELSE
              ne_all=0
            endif
            CALL smat%free()
            CALL hmat%free()
            DEALLOCATE(hmat,smat, stat=dealloc_stat, errmsg=errmsg)
            if(dealloc_stat /= 0) call juDFT_error("deallocate failed for hmat or smat",&
                                                   hint=errmsg, calledby="eigen.F90")

            ! Output results
            CALL timestart("EV output")
            ne_found=ne_all
            if (fmpi%pe_diag) THEN
#if defined(CPP_MPI)
              ! Collect number of all eigenvalues
              CALL MPI_ALLREDUCE(ne_found,ne_all,1,MPI_INTEGER,MPI_SUM, fmpi%diag_sub_comm,ierr)
              ne_all=MIN(fi%input%neig,ne_all)
#endif
            endif

            IF (fmpi%n_rank == 0) THEN
                ! Only process 0 writes out the value of ne_all and the
                ! eigenvalues.
#ifdef CPP_MPI
                call MPI_COMM_RANK(fmpi%diag_sub_comm,n_rank,err)
                call MPI_COMM_SIZE(fmpi%diag_sub_comm,n_size,err)
#else
                n_rank = 0; n_size=1;
#endif
                !IF (.NOT.l_dfpteigen) THEN
                    CALL write_eig(eig_id, nk,jsp,ne_found,ne_all,&
                               eig(:ne_all),n_start=n_size,n_end=n_rank,zMat=zMat)
                !ELSE
                !    CALL dfpt_eigen()
                !    RETURN
                !END IF
                eigBuffer(:ne_all,nk,jsp) = eig(:ne_all)
            ELSE
                !IF (.NOT.l_dfpteigen) THEN
                    if (fmpi%pe_diag) CALL write_eig(eig_id, nk,jsp,ne_found,&
                                  n_start=fmpi%n_size,n_end=fmpi%n_rank,zMat=zMat)
                !ELSE
                !    if (fmpi%pe_diag) CALL dfpt_eigen()
                !    RETURN
                !END IF
            ENDIF
            neigBuffer(nk,jsp) = ne_found
#if defined(CPP_MPI)
            ! RMA synchronization
            CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
            CALL timestop("EV output")

            IF (fi%banddos%unfoldband .AND. (.NOT. fi%noco%l_soc)) THEN
               IF(modulo (fi%kpts%nkpt,fmpi%n_size).NE.0) call juDFT_error("number fi%kpts needs to be multiple of number fmpi threads",&
                   hint=errmsg, calledby="eigen.F90")
               CALL calculate_plot_w_n(fi%banddos,fi%cell,fi%kpts,zMat,lapw,nk,jsp,eig,results,fi%input,fi%atoms,unfoldingBuffer,fmpi,fi%noco%l_soc,smat_unfold=smat_unfold)
               CALL smat_unfold%free()
               DEALLOCATE(smat_unfold, stat=dealloc_stat, errmsg=errmsg)
               if(dealloc_stat /= 0) call juDFT_error("deallocate failed for smat_unfold",&
                                                      hint=errmsg, calledby="eigen.F90")
            END IF

            IF (allocated(zmat)) THEN
              call zMat%free()
              deallocate(zMat)
            ENDIF
         END DO  k_loop
      END DO ! spin loop ends

      neigd2 = MIN(fi%input%neig,lapw%dim_nbasfcn())
#ifdef CPP_MPI
      IF (fi%banddos%unfoldband .AND. (.NOT. fi%noco%l_soc)) THEN
         results%unfolding_weights = CMPLX(0.0,0.0)
       CALL MPI_ALLREDUCE(unfoldingBuffer,results%unfolding_weights,SIZE(results%unfolding_weights,1)*SIZE(results%unfolding_weights,2)*SIZE(results%unfolding_weights,3),CPP_MPI_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
      END IF
      CALL MPI_ALLREDUCE(neigBuffer,results%neig,fi%kpts%nkpt*fi%input%jspins,MPI_INTEGER,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(eigBuffer(:neigd2,:,:),results%eig(:neigd2,:,:),neigd2*fi%kpts%nkpt*fi%input%jspins,MPI_DOUBLE_PRECISION,MPI_MIN,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(nvBuffer(:,:),nvBufferTemp(:,:),size(nvbuffer),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      IF (PRESENT(nvfull)) THEN
          CALL MPI_ALLREDUCE(nvfullBuffer(:,:),nvfull(:,:),size(nvfullBuffer),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(GbasVecBuffer(:,:,:,:),GbasVec_eig(:,:,:,:),size(GbasVecBuffer),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      END IF
      CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#else
      results%neig(:,:) = neigBuffer(:,:)
      results%eig(:neigd2,:,:) = eigBuffer(:neigd2,:,:)
      results%unfolding_weights(:,:,:) = unfoldingBuffer(:,:,:)
      nvBufferTemp(:,:) = nvBuffer(:,:)
      IF (PRESENT(nvfull)) THEN
          nvfull(:,:) = nvfullBuffer(:,:)
          GbasVec_eig(:,:,:,:) = GbasVecBuffer(:,:,:,:)
      END IF
#endif

      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,'(a)') ''
         WRITE(oUnit,'(a)') '              basis set size:'
         WRITE(oUnit,'(a)') '      jsp    ikpt     nv      LOs  overall'
         DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO nk = 1, fi%kpts%nkpt
               WRITE(oUnit,'(5i8)') jsp, nk, nvBufferTemp(nk,jsp), fi%atoms%nlotot, nvBufferTemp(nk,jsp) + fi%atoms%nlotot
            END DO
         END DO
      END IF


      IF( fi%input%jspins .EQ. 1 .AND. fi%hybinp%l_hybrid ) THEN
         results%te_hfex%valence = 2*results%te_hfex%valence
         IF(hybdat%l_calhf) results%te_hfex%core = 2*results%te_hfex%core
      END IF
      enpara%epara_min = MINVAL(enpara%el0)
      enpara%epara_min = MIN(MINVAL(enpara%ello0),enpara%epara_min)
   END SUBROUTINE eigen
END MODULE m_eigen
