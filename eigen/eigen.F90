!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!<@brief
!<The eigen routine sets up and solves the generalized eigenvalue problem
!More description at end of file
MODULE m_eigen
   USE m_juDFT
   IMPLICIT NONE
CONTAINS
   !>The eigenvalue problem is constructed and solved in this routine. The following steps are performed:
   !> 1. Preparation: generate energy parameters, open eig-file
   !> 2. CALL to mt_setup() : this constructs the local Hamiltonian (i.e. the Hamiltonian in the \f$ u,\dot u, u_{lo} \f$ basis) LDA+U is also added here
   !> 3. within the (collinear)spin and k-point loop: CALL to eigen_hssetup() to generate the matrices, CALL to eigen_diag() to perform diagonalization
   !> 4. writing (saving) of eigenvectors
   !>
   !> The matrices generated and diagonalized here are of type m_mat as defined in m_types_mat.
   !>@author D. Wortmann
   SUBROUTINE eigen(mpi,stars,sphhar,atoms,obsolete,xcpot,sym,kpts,DIMENSION,vacuum,input,&
                    cell,enpara,banddos,noco,oneD,hybrid,iter,eig_id,results,inden,v,vx)

#include"cpp_double.h"
      USE m_constants, ONLY : pi_const,sfp_const
      USE m_types
      USE m_eigen_hssetup
      USE m_pot_io
      USE m_eigen_diag
      USE m_add_vnonlocal
      USE m_subvxc
      !USE m_hsefunctional
      USE m_mt_setup
      USE m_util
      USE m_io_hybrid
      !USE m_icorrkeys
      USE m_eig66_io, ONLY : open_eig, write_eig, close_eig,read_eig
      USE m_xmlOutput
#ifdef CPP_MPI
      USE m_mpi_bc_potden
#endif
      USE m_symmetrize_matrix
      USE m_unfold_band_kpts !used for unfolding bands
      USE m_types_mpimat

      IMPLICIT NONE
      TYPE(t_results),INTENT(INOUT):: results
      CLASS(t_xcpot),INTENT(IN)    :: xcpot
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_dimension),INTENT(IN) :: DIMENSION
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_hybrid),INTENT(INOUT) :: hybrid
      TYPE(t_enpara),INTENT(INOUT) :: enpara
      TYPE(t_obsolete),INTENT(IN)  :: obsolete
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_banddos),INTENT(IN)   :: banddos
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_kpts),INTENT(INOUT)   :: kpts
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_potden),INTENT(IN)    :: inden,vx
      TYPE(t_potden),INTENT(INOUT) :: v    !u_setup will modify the potential matrix

#ifdef CPP_MPI
      INCLUDE 'mpif.h'
#endif

!    EXTERNAL MPI_BCAST    !only used by band_unfolding to broadcast the gvec

      ! Scalar Arguments
      INTEGER,INTENT(IN)    :: iter
      INTEGER,INTENT(IN)    :: eig_id

      ! Local Scalars
      INTEGER jsp,nk,nred,ne_all,ne_found
      INTEGER ne,lh0
      INTEGER isp,i,j,err
      LOGICAL l_wu,l_file,l_real,l_zref
      INTEGER :: solver=0
      ! Local Arrays
      INTEGER              :: ierr(3)
      INTEGER              :: neigBuffer(kpts%nkpt,input%jspins)

      COMPLEX              :: unfoldingBuffer(SIZE(results%unfolding_weights,1),kpts%nkpt,input%jspins) ! needed for unfolding bandstructure mpi case

      INTEGER, PARAMETER   :: lmaxb = 3
      REAL,    ALLOCATABLE :: bkpt(:)
      REAL,    ALLOCATABLE :: eig(:)
      COMPLEX, ALLOCATABLE :: vs_mmp(:,:,:,:)

      INTEGER                   :: jsp_m, i_kpt_m, i_m

      TYPE(t_tlmplm)            :: td
      TYPE(t_usdus)             :: ud
      TYPE(t_lapw)              :: lapw
      CLASS(t_mat), ALLOCATABLE :: zMat
      CLASS(t_mat), ALLOCATABLE :: hmat,smat
      CLASS(t_mat), ALLOCATABLE :: smat_unfold !used for unfolding bandstructure

      ! Variables for HF or hybrid functional calculation
      INTEGER                   :: comm(kpts%nkpt),irank2(kpts%nkpt),isize2(kpts%nkpt), dealloc_stat
      character(len=300)        :: errmsg

      call ud%init(atoms,input%jspins)
      ALLOCATE (eig(DIMENSION%neigd),bkpt(3))

      l_real=sym%invs.AND..NOT.noco%l_noco

      ! check if z-reflection trick can be used
      l_zref=(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco)
      IF (mpi%n_size > 1) l_zref = .FALSE.

#ifdef CPP_MPI
      CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,v)
#endif

      !IF (mpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/),(/iter,v%iter/),&
      !                                                RESHAPE((/19,13,5,5/),(/2,2/)))

      ! Set up and solve the eigenvalue problem
      !   loop over spins
      !     set up k-point independent t(l'm',lm) matrices
      CALL mt_setup(atoms,sym,sphhar,input,noco,enpara,inden,v,mpi,results,DIMENSION,td,ud)

      neigBuffer = 0
      results%neig = 0
      results%eig = 1.0e300
      unfoldingBuffer = CMPLX(0.0,0.0)

      DO jsp = 1,MERGE(1,input%jspins,noco%l_noco)
         k_loop:DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride

            ! Set up lapw list
            CALL lapw%init(input,noco, kpts,atoms,sym,nk,cell,l_zref, mpi)
            call timestart("Setup of H&S matrices")
            CALL eigen_hssetup(jsp,mpi,DIMENSION,hybrid,enpara,input,vacuum,noco,sym,&
                               stars,cell,sphhar,atoms,ud,td,v,lapw,l_real,smat,hmat)
            CALL timestop("Setup of H&S matrices")

            IF(hybrid%l_hybrid.OR.input%l_rdmft) THEN

               DO i = 1, hmat%matsize1
                  DO j = 1, i
                     IF (hmat%l_real) THEN
                        IF ((i.LE.5).AND.(j.LE.5)) THEN
                           WRITE(1233,'(2i7,2f15.8)') i, j, hmat%data_r(i,j), hmat%data_r(j,i)
                        END IF
                     ELSE
                        IF ((i.LE.5).AND.(j.LE.5)) THEN
                           WRITE(1233,'(2i7,4f15.8)') i, j, hmat%data_c(i,j), hmat%data_c(j,i)
                        END IF
                     ENDIF
                  END DO
               END DO

               ! Write overlap matrix smat to direct access file olap
               print *,"Wrong overlap matrix used, fix this later"
               CALL write_olap(smat,(jsp-1)*kpts%nkpt+nk) ! Note: At this moment this only works without MPI parallelization
            END IF ! hybrid%l_hybrid.OR.input%l_rdmft

            IF(hybrid%l_hybrid) THEN
               PRINT *,"TODO"
!             STOP "TODO"
               PRINT *,"BASIS:", lapw%nv(jsp), atoms%nlotot
               IF (hybrid%l_addhf) CALL add_Vnonlocal(nk,lapw,atoms,hybrid,dimension,kpts,jsp,results,xcpot,noco,hmat)

               IF(hybrid%l_subvxc) THEN
                  CALL subvxc(lapw,kpts%bk(:,nk),DIMENSION,input,jsp,v%mt(:,0,:,:),atoms,ud,hybrid,enpara%el0,enpara%ello0,&
                              sym,cell,sphhar,stars,xcpot,mpi,oneD,hmat,vx)
               END IF
            END IF ! hybrid%l_hybrid

            l_wu=.FALSE.
            ne_all=DIMENSION%neigd
            IF(ne_all < 0) ne_all = lapw%nmat
            IF(ne_all > lapw%nmat) ne_all = lapw%nmat

            !Try to symmetrize matrix
            CALL symmetrize_matrix(mpi,noco,kpts,nk,hmat,smat)

            IF (banddos%unfoldband) THEN
               select type(smat)
               type is (t_mat)
                  allocate(t_mat::smat_unfold)
                  select type(smat_unfold)
                  type is (t_mat)
                     smat_unfold=smat
                  end select
               type is (t_mpimat)
                  allocate(t_mpimat::smat_unfold)
                  select type(smat_unfold)
                  type is (t_mpimat)
                     smat_unfold=smat
                  end select
               end select
            END IF

            ! Solve generalized eigenvalue problem.
            !     ne_all ... number of eigenpairs searched (and found) on this node
            !                on input, overall number of eigenpairs searched,
            !                on output, local number of eigenpairs found
            !     eig ...... all eigenvalues, output
            !     zMat ..... local eigenvectors, output
            CALL eigen_diag(solver,hmat,smat,ne_all,eig,zMat,nk,jsp,iter)
              
            CALL smat%free()
            CALL hmat%free()
            DEALLOCATE(hmat,smat, stat=dealloc_stat, errmsg=errmsg)
            if(dealloc_stat /= 0) call juDFT_error("deallocate failed for hmat or smat",&
                                                   hint=errmsg, calledby="eigen.F90")

            ! Output results
            CALL timestart("EV output")
#if defined(CPP_MPI)
            ! Collect number of all eigenvalues
            ne_found=ne_all
            CALL MPI_ALLREDUCE(ne_found,ne_all,1,MPI_INTEGER,MPI_SUM, mpi%sub_comm,ierr)
            ne_all=MIN(DIMENSION%neigd,ne_all)
#else
            ne_found=ne_all
#endif
            IF (.NOT.zMat%l_real) THEN
               zMat%data_c(:lapw%nmat,:ne_found) = CONJG(zMat%data_c(:lapw%nmat,:ne_found))
            END IF
            IF (mpi%n_rank == 0) THEN
                ! Only process 0 writes out the value of ne_all and the
                ! eigenvalues. 
                ! Trying to use MPI_PUT for the very same slot by all processes
                ! causes problems with IntelMPI/2019
                !        Mai 2019                 U. Alekseeva      
                CALL write_eig(eig_id, nk,jsp,ne_found,ne_all,&
                           eig(:ne_all),n_start=mpi%n_size,n_end=mpi%n_rank,zMat=zMat)
            ELSE
                CALL write_eig(eig_id, nk,jsp,ne_found,&
                           n_start=mpi%n_size,n_end=mpi%n_rank,zMat=zMat)
            ENDIF
            neigBuffer(nk,jsp) = ne_found
#if defined(CPP_MPI)
            ! RMA synchronization
            CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
            CALL timestop("EV output")

            IF (banddos%unfoldband) THEN
               CALL calculate_plot_w_n(banddos,cell,kpts,smat_unfold,zMat,lapw,nk,jsp,eig,results,input,atoms,unfoldingBuffer,mpi)
               CALL smat_unfold%free()
               DEALLOCATE(smat_unfold, stat=dealloc_stat, errmsg=errmsg)
               if(dealloc_stat /= 0) call juDFT_error("deallocate failed for smat_unfold",&
                                                      hint=errmsg, calledby="eigen.F90")
            END IF

            call zMat%free()
            deallocate(zMat)
         END DO  k_loop
      END DO ! spin loop ends

#ifdef CPP_MPI
      IF (banddos%unfoldband) THEN
         results%unfolding_weights = CMPLX(0.0,0.0)
       CALL MPI_ALLREDUCE(unfoldingBuffer,results%unfolding_weights,SIZE(results%unfolding_weights,1)*SIZE(results%unfolding_weights,2)*SIZE(results%unfolding_weights,3),CPP_MPI_COMPLEX,MPI_SUM,mpi%mpi_comm,ierr)
      END IF
      CALL MPI_ALLREDUCE(neigBuffer,results%neig,kpts%nkpt*input%jspins,MPI_INTEGER,MPI_SUM,mpi%mpi_comm,ierr)
      CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#else
      results%neig(:,:) = neigBuffer(:,:)
      results%unfolding_weights(:,:,:) = unfoldingBuffer(:,:,:)
#endif

      ! Sorry for the following strange workaround to fill the results%eig array.
      ! At some point someone should have a closer look at how the eigenvalues are
      ! distributed and fill the array without using the eigenvalue-IO.
      DO jsp = 1,MERGE(1,input%jspins,noco%l_noco)
         DO nk = 1,kpts%nkpt
            CALL read_eig(eig_id,nk,jsp,results%neig(nk,jsp),results%eig(:,nk,jsp))
#ifdef CPP_MPI
            CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
         END DO
      END DO

      !IF (hybrid%l_hybrid.OR.hybrid%l_calhf) CALL close_eig(eig_id)

      IF( input%jspins .EQ. 1 .AND. hybrid%l_hybrid ) THEN
         results%te_hfex%valence = 2*results%te_hfex%valence
         IF(hybrid%l_calhf) results%te_hfex%core = 2*results%te_hfex%core
      END IF
      enpara%epara_min = MINVAL(enpara%el0)
      enpara%epara_min = MIN(MINVAL(enpara%ello0),enpara%epara_min)
   END SUBROUTINE eigen
END MODULE m_eigen

