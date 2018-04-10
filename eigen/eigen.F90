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
  SUBROUTINE eigen(mpi,stars,sphhar,atoms,obsolete,xcpot,&
       sym,kpts,DIMENSION, vacuum, input, cell, enpara,banddos, noco, oneD,hybrid,&
       it,eig_id,results,inden,v,vx)
    USE m_constants, ONLY : pi_const,sfp_const
    USE m_types
    USE m_apws
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
    USE m_mpi_bc_pot
#endif

    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
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
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(INOUT)  :: atoms!in u_setup n_u might be modified
    TYPE(t_potden),INTENT(IN)    :: inden
    TYPE(t_potden),INTENT(INOUT) :: v,vx
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN) :: it
    INTEGER,INTENT(INOUT):: eig_id
    !     ..
    !-odim
    !+odim
    !     ..
    !     .. Local Scalars ..
    INTEGER jsp,nk,nred,ne_all,ne_found
    INTEGER ne  ,lh0
    INTEGER isp,i,j,err
    LOGICAL l_wu,l_file,l_real,l_zref
    
    !     ..
    !     .. Local Arrays ..
    INTEGER, PARAMETER :: lmaxb=3
    REAL,    ALLOCATABLE :: bkpt(:)
    REAL,    ALLOCATABLE :: eig(:)

    COMPLEX, ALLOCATABLE :: vs_mmp(:,:,:,:)
    TYPE(t_tlmplm)  :: td
    TYPE(t_usdus)   :: ud
    TYPE(t_lapw)    :: lapw
    CLASS(t_Mat),ALLOCATABLE    :: zMat
    CLASS(t_mat),ALLOCATABLE    :: hmat,smat
    TYPE(T_mat)     :: olap
    !
    INTEGER nn,n
    INTEGER ierr(3)


    !     .. variables for HF or hybrid functional calculation ..
    !
    INTEGER                 ::  comm(kpts%nkpt),irank2(kpts%nkpt),isize2(kpts%nkpt)
    
    !
    !
    ! --> Allocate
    !
    call ud%init(atoms,DIMENSION%jspd)
    ALLOCATE ( eig(DIMENSION%neigd),bkpt(3) )
    
    !
    ! --> some parameters first
    
    l_real=sym%invs.AND..NOT.noco%l_noco
    !check if z-reflection trick can be used

    l_zref=(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco) 
    IF (mpi%n_size > 1) l_zref = .FALSE.
 
#ifdef CPP_MPI
    CALL mpi_bc_pot(mpi,stars,sphhar,atoms,input,vacuum,&
                    v%iter,v%mt,v%pw,v%vacz,v%vacxy)
#endif
    IF (mpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/),(/it,v%iter/),&
                                                    RESHAPE((/19,13,5,5/),(/2,2/)))
    
     eig_id=open_eig(&
          mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,kpts%nkpt,DIMENSION%jspd,atoms%lmaxd,&
         atoms%nlod,atoms%ntype,atoms%nlotot,noco%l_noco,.TRUE.,l_real,noco%l_soc,.FALSE.,&
         mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
         nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
         l_orb=banddos%l_orb)
     !
     !---> set up and solve the eigenvalue problem
     !--->    loop over spins
     !--->       set up k-point independent t(l'm',lm) matrices
     !
     CALL mt_setup(atoms,sym,sphhar,input,noco,enpara,inden,v,mpi,results,DIMENSION,td,ud)
   
    DO jsp = 1,MERGE(1,input%jspins,noco%l_noco)
       k_loop:DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride

          !
          !--->         set up lapw list
        
          CALL lapw%init(input,noco, kpts,atoms,sym,nk,cell,l_zref, mpi)
          call timestart("Setup of H&S matrices")
          CALL eigen_hssetup(jsp,mpi,DIMENSION,hybrid,enpara,input,vacuum,noco,sym,&
               stars,cell,sphhar,atoms,ud,td,v,lapw,l_real,smat,hmat)
          CALL timestop("Setup of H&S matrices")
        
          IF( hybrid%l_hybrid ) THEN
             !write overlap matrix b to direct access file olap
             print *,"Wrong overlap matrix used, fix this later"
             !call write_olap(smat,nrec)
             STOP "TODO"
             PRINT *,"BASIS:",lapw%nv(jsp),atoms%nlotot
             !if (hybrid%l_addhf) CALL add_Vnonlocal(nk,hybrid,dimension, kpts,jsp,results,xcpot,hamovlp)
             
             
             IF( hybrid%l_subvxc ) THEN
                !CALL subvxc(lapw,kpts%bk(:,nk),DIMENSION,input,jsp,v%mt(:,0,:,:),atoms,ud,hybrid,enpara%el0,enpara%ello0,&
                 !    sym, atoms%nlotot,kveclo, cell,sphhar, stars, xcpot,mpi,&
                 !    oneD,  hamovlp,vx)
             END IF

          END IF ! hybrid%l_hybrid
          l_wu=.FALSE.
          ne_all=DIMENSION%neigd
          if (allocated(zmat)) deallocate(zmat)
          CALL eigen_diag(hmat,smat,ne_all,eig,zMat)
          DEALLOCATE(hmat,smat)
          !
          !--->         output results
          !
          CALL timestart("EV output")
#if defined(CPP_MPI)
          !Collect number of all eigenvalues
          ne_found=ne_all
          CALL MPI_ALLREDUCE(ne_found,ne_all,1,MPI_INTEGER,MPI_SUM, mpi%sub_comm,ierr)
          ne_all=MIN(DIMENSION%neigd,ne_all)
#else
          ne_found=ne_all
#endif          
          IF (.NOT.l_real) THEN
                zMat%data_c(:lapw%nmat,:ne_found) = CONJG(zMat%data_c(:lapw%nmat,:ne_found))
          ENDIF
    	  CALL write_eig(eig_id, nk,jsp,ne_found,ne_all,&
                  eig(:ne_found),n_start=mpi%n_size,n_end=mpi%n_rank,zmat=zMat)
#if defined(CPP_MPI)
          !RMA synchronization
          CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif

          CALL timestop("EV output")
          !
       END DO  k_loop

    END DO ! spin loop ends

  
#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
    !IF (hybrid%l_hybrid.OR.hybrid%l_calhf) CALL close_eig(eig_id)
  

    IF( input%jspins .EQ. 1 .AND. hybrid%l_hybrid ) THEN
       results%te_hfex%valence = 2*results%te_hfex%valence
       results%te_hfex%core    = 2*results%te_hfex%core
    END IF
    enpara%epara_min = MINVAL(enpara%el0)
    enpara%epara_min = MIN(MINVAL(enpara%ello0),enpara%epara_min)
  END SUBROUTINE eigen
END MODULE m_eigen

