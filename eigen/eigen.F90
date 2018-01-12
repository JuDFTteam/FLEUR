!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen
  USE m_juDFT
CONTAINS
  SUBROUTINE eigen(mpi,stars,sphhar,atoms,obsolete,xcpot,&
       sym,kpts,DIMENSION, vacuum, input, cell, enpara_in,banddos, noco,jij, oneD,hybrid,&
       it,eig_id,results,v,vx)
    !*********************************************************************
    !     sets up and solves the eigenvalue problem for a basis of lapws.
    !
    ! nv,   nvd     ... actual length & dimension of EV without LO's
    ! nmat, nbasfcn                                   including LO's
    !        g. bihlmayer '96
    !**********************************************************************
    USE m_constants, ONLY : pi_const,sfp_const
    USE m_types
    USE m_lodpot
  
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
    TYPE(t_enpara),INTENT(INOUT) :: enpara_in
    TYPE(t_obsolete),INTENT(IN)  :: obsolete
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_banddos),INTENT(IN)   :: banddos
    TYPE(t_jij),INTENT(IN)       :: jij
    TYPE(t_sym),INTENT(IN)       :: sym  
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(INOUT)  :: atoms!in u_setup n_u might be modified
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
    INTEGER nspins,isp,i,j,err
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
    TYPE(t_enpara)  :: enpara
    TYPE(t_zMat)    :: zMat
    TYPE(t_lapwmat) :: hmat,smat
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
   
    !
    ! set energy parameters (normally to that, what we read in)
    !
    CALL lodpot(mpi,atoms,sphhar,obsolete,vacuum,&
            input, v%mt,v%vacz, enpara_in, enpara)
    !
    
     eig_id=open_eig(&
          mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,kpts%nkpt,DIMENSION%jspd,atoms%lmaxd,&
         atoms%nlod,atoms%ntype,atoms%nlotot,noco%l_noco,.TRUE.,l_real,noco%l_soc,.FALSE.,&
         mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
         nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
         l_orb=banddos%l_orb)
    !
    !---> set up and solve the eigenvalue problem
    !--->    loop over spins
    nspins = input%jspins
    IF (noco%l_noco) nspins = 1
    !
    !  ..
    !  LDA+U
  
    !
    !--->    loop over k-points: each can be a separate task

    DO jsp = 1,nspins
       !
       !--->       set up k-point independent t(l'm',lm) matrices
       !
       smat%l_real=l_real;hmat%l_real=l_real
       CALL mt_setup(jsp,atoms,sym,sphhar,input,noco,enpara,v,mpi,results,DIMENSION,td,ud)
       k_loop:DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride

          !
          !--->         set up lapw list
        
          CALL apws(DIMENSION,input,noco, kpts,atoms,sym,nk,cell,l_zref, mpi%n_size,jsp, bkpt,lapw,nred)

          CALL eigen_hssetup(jsp,mpi,DIMENSION,oned,hybrid,enpara,input,vacuum,noco,jij,sym,&
               stars,cell,kpts,sphhar,atoms,ud,td,v,bkpt,lapw,smat,hmat)

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
          lapw%nmat=lapw%nv(1)+atoms%nlotot
          IF (noco%l_noco) lapw%nmat=lapw%nmat+lapw%nv(2)+atoms%nlotot
          l_wu=.false.
          CALL eigen_diag(jsp,eig_id,it,atoms,DIMENSION,mpi, mpi%n_rank,mpi%n_size, DIMENSION%neigd,nk,lapw,input,&
               nred,mpi%sub_comm, sym,l_zref, noco,cell,bkpt,enpara%el0,jij,l_wu,&
               oneD,td,ud, eig,ne_found,smat,hmat,zMat)
          
          !
          !--->         output results
          !
          CALL timestart("EV output")
          ne_all=ne_found
#if defined(CPP_MPI)
          !Collect number of all eigenvalues
          CALL MPI_ALLREDUCE(ne_found,ne_all,1,MPI_INTEGER,MPI_SUM, mpi%sub_comm,ierr)
          ne_all=MIN(DIMENSION%neigd,ne_all)
#endif
          !jij%eig_l = 0.0 ! need not be used, if hdf-file is present
          IF (.NOT.l_real) THEN
             IF (.NOT.jij%l_J) THEN
                zMat%z_c(:lapw%nmat,:ne_found) = CONJG(zMat%z_c(:lapw%nmat,:ne_found))
             ELSE
                zMat%z_c(:lapw%nmat,:ne_found) = CMPLX(0.0,0.0)
             ENDIF
          ENDIF
	  zmat%nbands=ne_found
          CALL write_eig(eig_id, nk,jsp,ne_found,ne_all,lapw%nv(jsp),lapw%nmat,&
                  lapw%k1(:lapw%nv(jsp),jsp),lapw%k2 (:lapw%nv(jsp),jsp),lapw%k3(:lapw%nv(jsp),jsp),&
                  bkpt, kpts%wtkpt(nk),eig(:ne_found),el=enpara%el0(0:,:,jsp),ello=enpara%ello0(:,:,jsp),evac=enpara%evac0(:,jsp),&
                  nlotot=atoms%nlotot,n_start=mpi%n_size,n_end=mpi%n_rank,zmat=zMat)
          IF (noco%l_noco) THEN
             CALL write_eig(eig_id, nk,2,ne_found,ne_all,lapw%nv(2),lapw%nmat,&
                  lapw%k1(:lapw%nv(2),2),lapw%k2 (:lapw%nv(2),2),lapw%k3(:lapw%nv(2),2),&
                  bkpt, kpts%wtkpt(nk),eig(:ne_found),el=enpara%el0(0:,:,2),ello= enpara%ello0(:,:,2),evac=enpara%evac0(:,2),&
                  nlotot=atoms%nlotot)
          ENDIF
#if defined(CPP_MPI)
          !RMA synchronization
          CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif

          CALL timestop("EV output")
          !#ifdef CPP_MPI
          IF (l_real) THEN
             DEALLOCATE ( zMat%z_r )
          ELSE
             DEALLOCATE ( zMat%z_c )
          ENDIF
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
    enpara_in%epara_min = MINVAL(enpara%el0)
    enpara_in%epara_min = MIN(MINVAL(enpara%ello0),enpara_in%epara_min)
  END SUBROUTINE eigen
END MODULE m_eigen
