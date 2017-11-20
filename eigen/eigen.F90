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
    USE m_tlmplm_cholesky
    USE m_tlmplm_store
    USE m_apws
    USE m_hsmt_simple
    USE m_hs_int
    USE m_hsvac
    USE m_od_hsvac
    USE m_usetup
    USE m_pot_io
    USE m_eigen_diag
    USE m_add_vnonlocal
    USE m_subvxc
    !USE m_hsefunctional
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
    INTEGER jsp,nk,nred,ne_all,n_u_in,ne_found
    INTEGER ne  ,nrec,lh0
    INTEGER nspins,isp,i,j,err
    INTEGER mlotot,mlolotot
    LOGICAL l_wu,l_file,l_real,l_zref
    
    !     ..
    !     .. Local Arrays ..
    INTEGER, PARAMETER :: lmaxb=3
    INTEGER, ALLOCATABLE :: matind(:,:),kveclo(:)
    INTEGER, ALLOCATABLE :: nv2(:)
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
    INTEGER fh,nn,n
    INTEGER ierr(3)


    !     .. variables for HF or hybrid functional calculation ..
    !
    INTEGER                 ::  comm(kpts%nkpt),irank2(kpts%nkpt),isize2(kpts%nkpt)
    
#ifdef CPP_MPI
    INTEGER   :: sndreqd,sndreq(mpi%isize*kpts%nkpt)
#endif
    !
    !
    ! --> Allocate
    !
    call ud%init(atoms,DIMENSION%jspd)
    ALLOCATE ( nv2(DIMENSION%jspd) )
    ALLOCATE ( eig(DIMENSION%neigd),bkpt(3) )
   
    !
    ! --> some parameters first
    !
    !     determine the total number of lo's : nlotot
    !
    mlotot =sum(atoms%nlo)
  
    mlolotot = 0
    DO nn = 1, atoms%ntype
       mlolotot = mlolotot + atoms%nlo(nn)*(atoms%nlo(nn)+1)/2
    ENDDO
    ALLOCATE ( kveclo(atoms%nlotot) )
    !     ..
    
    l_real=sym%invs.AND..NOT.noco%l_noco
    IF (noco%l_soc.AND.l_real.AND.hybrid%l_hybrid ) THEN
       CALL juDFT_error('hybrid functional + SOC + inv.symmetry is not tested', calledby='eigen')
    END IF

    
    fh = 0
    !
    ! look, if WU diagonalisation
    !
    IF (it.LT.input%isec1) THEN
       IF (mpi%irank.EQ.0) WRITE (6,FMT=8110) it,input%isec1
8110   FORMAT (' IT=',i4,'  ISEC1=',i4,' standard diagonalization')
       l_wu = .FALSE.
    ELSE
       IF (mpi%irank.EQ.0) WRITE (6,FMT=8120) it,input%isec1
8120   FORMAT (' IT=',i4,'  ISEC1=',i4,' reduced diagonalization')
       l_wu = .TRUE.
    END IF
   
    !IF (mpi%irank.EQ.0) THEN
    !   CALL readPotential(stars,vacuum,atoms,sphhar,input,sym,POT_ARCHIVE_TYPE_TOT_const,&
    !                      v%iter,v%mt,v%pw,v%vacz,v%vacxy)
    !END IF
#ifdef CPP_MPI
    CALL mpi_bc_pot(mpi,stars,sphhar,atoms,input,vacuum,&
                    v%iter,v%mt,v%pw,v%vacz,v%vacxy)
#endif

999 CONTINUE
    IF (mpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/),(/it,v%iter/),&
                                                    RESHAPE((/19,13,5,5/),(/2,2/)))

   
    !
    ! set energy parameters (normally to that, what we read in)
    !
    CALL lodpot(mpi,atoms,sphhar,obsolete,vacuum,&
            input, v%mt,v%vacz, enpara_in, enpara)
    !


    !---> set up and solve the eigenvalue problem
    !---> loop over energy windows

!check if z-reflection trick can be used

    l_zref=(sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkpt))).LT.1e-9).AND..NOT.noco%l_noco) 


#if ( defined(CPP_MPI))
    IF (mpi%n_size > 1) l_zref = .FALSE.
    IF ( hybrid%l_calhf ) THEN
       CALL judft_error("BUG parallelization in HF case must be fixed")
       !n_start  = 1
       !n_stride = 1
    END IF
#endif
    ne = DIMENSION%neigd

    eig_id=open_eig(&
          mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,kpts%nkpt,DIMENSION%jspd,atoms%lmaxd,&
         atoms%nlod,atoms%ntype,atoms%nlotot,noco%l_noco,.TRUE.,l_real,noco%l_soc,.FALSE.,&
         mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
         nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
         l_orb=banddos%l_orb)

    
   
    ALLOCATE (  matind(DIMENSION%nbasfcn,2) )
    !
    !--->    loop over spins
    nspins = input%jspins
    IF (noco%l_noco) nspins = 1
    !
    !  ..
    !  LDA+U
    n_u_in=atoms%n_u
    IF ((atoms%n_u.GT.0)) THEN
       ALLOCATE( vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins) )
       CALL u_setup(sym,atoms,lmaxb,sphhar,input, enpara%el0(0:,:,:),v%mt,mpi, vs_mmp,results)
    ELSE
       ALLOCATE( vs_mmp(-lmaxb:-lmaxb,-lmaxb:-lmaxb,1,2) )
    ENDIF
    !
    !--->    loop over k-points: each can be a separate task

    DO jsp = 1,nspins
       

       !
       !--->       set up k-point independent t(l'm',lm) matrices
       !
       CALL timestart("tlmplm")
       err=0
       j = 1 ; IF (noco%l_noco) j = 2
       ALLOCATE(td%tuu(0:DIMENSION%lmplmd,atoms%ntype,j),stat=err)
       ALLOCATE(td%tud(0:DIMENSION%lmplmd,atoms%ntype,j),stat=err)
       ALLOCATE(td%tdd(0:DIMENSION%lmplmd,atoms%ntype,j),stat=err)
       ALLOCATE(td%tdu(0:DIMENSION%lmplmd,atoms%ntype,j),stat=err)
       mlotot = MAX(mlotot,1) 
       ALLOCATE(td%tdulo(0:DIMENSION%lmd,-atoms%llod:atoms%llod,mlotot,j),stat=err)
       ALLOCATE(td%tuulo(0:DIMENSION%lmd,-atoms%llod:atoms%llod,mlotot,j),stat=err)
       ALLOCATE(td%tuloulo(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,MAX(mlolotot,1),j), stat=err)
       ALLOCATE(td%ind(0:DIMENSION%lmd,0:DIMENSION%lmd,atoms%ntype,j),stat=err )
       ALLOCATE(td%h_loc(0:2*atoms%lmaxd*(atoms%lmaxd+2)+1,0:2*atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%ntype,j))
       IF (err.NE.0) THEN
          WRITE (*,*) 'eigen: an error occured during allocation of'
          WRITE (*,*) 'the tlmplm%tuu, tlmplm%tdd etc.: ',err,'  size: ',mlotot
          CALL juDFT_error("eigen: Error during allocation of tlmplm, tdd  etc.",calledby ="eigen")
       ENDIF
       lh0=1
       CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,mpi,v%mt(:,0,1,jsp),input,vs_mmp, td,ud)
       IF (input%l_f) CALL write_tlmplm(td,vs_mmp,atoms%n_u>0,1,jsp,input%jspins)
       CALL timestop("tlmplm")

       !---> pk non-collinear
       !--->       call tlmplm again for the second spin direction in
       !--->       each MT, because the t-matrices are needed for both
       !--->       spins at once in hsmt
       IF (noco%l_noco) THEN
          isp = 2
          CALL timestart("tlmplm")
          !CALL tlmplm(sphhar,atoms,DIMENSION,enpara,isp,isp,mpi, v%mt(1,0,1,isp),lh0,input, td,ud)
          IF (input%l_f) CALL write_tlmplm(td,vs_mmp,atoms%n_u>0,2,2,input%jspins)
          CALL timestop("tlmplm")
       ENDIF
       !

#ifdef CPP_MPI
       ! check that all sending operations are completed
       IF ( hybrid%l_calhf ) CALL MPI_WAITALL(sndreqd,sndreq,MPI_STATUSES_IGNORE,ierr)
#endif

       k_loop:DO nk = mpi%n_start,kpts%nkpt,mpi%n_stride

          nrec =  kpts%nkpt*(jsp-1) + nk
          !
          !--->         set up lapw list
          !
          CALL timestart("Setup of LAPW")

          CALL apws(DIMENSION,input,noco, kpts,nk,cell,l_zref, mpi%n_size,jsp, bkpt,lapw,matind,nred)
          lapw%nmat=lapw%nv(jsp)+atoms%nlotot
          CALL hmat%init(lapw,mpi,atoms%nlotot,noco%l_noco,input%jspins,l_real)
          CALL smat%init(lapw,mpi,atoms%nlotot,noco%l_noco,input%jspins,l_real)
          
          CALL timestop("Setup of LAPW")
          IF (noco%l_noco) THEN
             !--->         the file potmat contains the 2x2 matrix-potential in
             !--->         the interstitial region and the vacuum
             OPEN (25,FILE='potmat',FORM='unformatted', STATUS='old')
          ENDIF
          !
          !--->         set up interstitial hamiltonian and overlap matrices
          !
          CALL timestart("Interstitial Hamiltonian&Overlap")
          CALL hs_int(input,noco,stars,lapw,mpi,cell,jsp,kpts%bk(:,nk),v%pw,smat,hmat)

          

          CALL timestop("Interstitial Hamiltonian&Overlap")
          !
          !--->         update with sphere terms
          !
          IF (.NOT.l_wu) THEN
             CALL timestart("MT Hamiltonian&Overlap")
             CALL hsmt_simple(jsp,kpts%bk(:,nk),DIMENSION,input,sym,cell,atoms,lapw,td,noco,ud,enpara,hmat,smat)
             CALL timestop("MT Hamiltonian&Overlap")
          ENDIF
          !

         
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
          !
          !--->         update with vacuum terms
          !
          CALL timestart("Vacuum Hamiltonian&Overlap")
          IF (input%film .AND. .NOT.oneD%odi%d1) THEN
             CALL hsvac(vacuum,stars,DIMENSION, atoms, jsp,input,v%vacxy,v%vacz,enpara%evac0,cell, &
                  bkpt,lapw,sym, noco,jij, mpi%n_size,mpi%n_rank,nv2,l_real,hmat,smat)
          ELSEIF (oneD%odi%d1) THEN
        !     CALL od_hsvac(vacuum,stars,DIMENSION, oneD,atoms, jsp,input,v%vacxy(1,1,1,jsp),v%vacz, &
        !          enpara%evac0,cell, bkpt,lapw, oneD%odi%M,oneD%odi%mb,oneD%odi%m_cyl,oneD%odi%n2d, &
        !          mpi%n_size,mpi%n_rank,sym,noco,jij,nv2,l_real,hamOvlp)
          END IF
          CALL timestop("Vacuum Hamiltonian&Overlap")

          IF (noco%l_noco) CLOSE (25)

         
          
          
          CALL eigen_diag(jsp,eig_id,it,atoms,DIMENSION,mpi, mpi%n_rank,mpi%n_size,ne,nk,lapw,input,&
               nred,mpi%sub_comm, sym,l_zref,matind,kveclo, noco,cell,bkpt,enpara%el0,jij,l_wu,&
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
                  nlotot=atoms%nlotot,kveclo=kveclo,n_start=mpi%n_size,n_end=mpi%n_rank,zmat=zMat)
          IF (noco%l_noco) THEN
             CALL write_eig(eig_id, nk,2,ne_found,ne_all,lapw%nv(2),lapw%nmat,&
                  lapw%k1(:lapw%nv(2),2),lapw%k2 (:lapw%nv(2),2),lapw%k3(:lapw%nv(2),2),&
                  bkpt, kpts%wtkpt(nk),eig(:ne_found),el=enpara%el0(0:,:,2),ello= enpara%ello0(:,:,2),evac=enpara%evac0(:,2),&
                  nlotot=atoms%nlotot,kveclo=kveclo)
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

       DEALLOCATE (td%tuu,td%tud,td%tdu,td%tdd)
       DEALLOCATE (td%ind,td%tuulo,td%tdulo)
       DEALLOCATE (td%tuloulo)
    END DO ! spin loop ends
    DEALLOCATE( vs_mmp )
    DEALLOCATE (matind)
#if defined(CPP_MPI)&&defined(CPP_NEVER)
    IF ( hybrid%l_calhf ) DEALLOCATE (nkpt_EIBZ)
#endif

  
#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%MPI_COMM,ierr)
#endif
    !IF (hybrid%l_hybrid.OR.hybrid%l_calhf) CALL close_eig(eig_id)
    atoms%n_u=n_u_in


    IF( input%jspins .EQ. 1 .AND. hybrid%l_hybrid ) THEN
       results%te_hfex%valence = 2*results%te_hfex%valence
       results%te_hfex%core    = 2*results%te_hfex%core
    END IF
    enpara_in%epara_min = MINVAL(enpara%el0)
    enpara_in%epara_min = MIN(MINVAL(enpara%ello0),enpara_in%epara_min)
!    enpara_in=enpara
  END SUBROUTINE eigen
END MODULE m_eigen
