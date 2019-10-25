!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur
   IMPLICIT NONE
CONTAINS
   SUBROUTINE fleur_execute(mpi_comm)

    !     ***************************************************************
    !
    !     based on flapw7 (c.l.fu, m.weinert, e.wimmer):
    !     full potential linearized augmented plane wave method for thin
    !     films and superlattices (version 7 ---- general symmetry)
    !     symmetry part       ---  e.wimmer
    !     potential generator ---  c.l.fu,r.podloucky
    !     matrix elements     ---  m.weinert
    !     charge density      ---  c.l.fu
    !                                c.l.fu        1987
    !     2nd variation diagon.  --- r.-q. wu      1992
    !     forces a la Yu et al   --- r.podloucky   1995
    !     full relativistic core --- a.shick       1996
    !     broyden mixing         --- r.pentcheva   1996
    !     gga (pw91, pbe)        --- t.asada       1997
    !     local orbitals         --- p.kurz        1997
    !     automatic symmetry     --- w.hofer       1997
    !     core tails & start     --- r.abt         1998
    !     spin orbit coupling    --- a.shick,x.nie 1998
    !     non-colinear magnet.   --- p.kurz        1999
    !     one-dimensional        --- y.mokrousov   2002
    !     exchange parameters    --- m.lezaic      2004
    !
    !                       g.bihlmayer, s.bluegel 1999
    !     ***************************************************************
    !----------------------------------------
    ! this routine is the main PROGRAM

   USE m_types
   USE m_constants
   USE m_fleur_init
   USE m_optional
   USE m_cdn_io
   USE m_mixing_history
   USE m_qfix
   USE m_vgen
   USE m_vgen_coulomb
   USE m_writexcstuff
   USE m_vmatgen
   USE m_eigen
   USE m_eigenso
   USE m_fermie
   USE m_cdngen
   USE m_totale
   USE m_potdis
   USE m_mix
   USE m_xmlOutput
   USE m_juDFT_time
   USE m_calc_hybrid
   USE m_rdmft
   USE m_io_hybrid
   USE m_wann_optional
   USE m_wannier
   USE m_bs_comfort
   USE m_dwigner
   USE m_ylm
   USE m_metagga
   USE m_plot
   USE m_xcBfield
#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif
   USE m_eig66_io
   USE m_chase_diag
   USE m_writeBasis
   !$ USE omp_lib
   IMPLICIT NONE

   INTEGER, INTENT(IN)             :: mpi_comm

    TYPE(t_input)                   :: input
    TYPE(t_field)                   :: field, field2
    TYPE(t_dimension)               :: DIMENSION
    TYPE(t_atoms)                   :: atoms
    TYPE(t_sphhar)                  :: sphhar
    TYPE(t_cell)                    :: cell
    TYPE(t_stars)                   :: stars
    TYPE(t_sym)                     :: sym
    TYPE(t_noco)                    :: noco
    TYPE(t_vacuum)                  :: vacuum
    TYPE(t_sliceplot)               :: sliceplot
    TYPE(t_banddos)                 :: banddos
    TYPE(t_obsolete)                :: obsolete
    TYPE(t_enpara)                  :: enpara
    TYPE(t_results)                 :: results
    TYPE(t_kpts)                    :: kpts
    TYPE(t_hybrid)                  :: hybrid
    TYPE(t_oneD)                    :: oneD
    TYPE(t_mpi)                     :: mpi
    TYPE(t_coreSpecInput)           :: coreSpecInput
    TYPE(t_wann)                    :: wann
    TYPE(t_potden)                  :: vTot, vx, vCoul, vTemp, vxcForPlotting
    TYPE(t_potden)                  :: inDen, outDen, EnergyDen, dummyDen
    TYPE(t_potden), DIMENSION(3)    :: testDen

    CLASS(t_xcpot),     ALLOCATABLE :: xcpot
    CLASS(t_forcetheo), ALLOCATABLE :: forcetheo

    ! local scalars
    INTEGER :: eig_id,archiveType, num_threads
    INTEGER :: iter,iterHF,i
    LOGICAL :: l_opti,l_cont,l_qfix,l_real
    REAL    :: fix
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: ierr(2),n
#endif

    mpi%mpi_comm = mpi_comm

    CALL timestart("Initialization")
    CALL fleur_init(mpi,input,field,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,sliceplot,&
                    banddos,obsolete,enpara,xcpot,results,kpts,hybrid,oneD,coreSpecInput,wann,l_opti)
    CALL timestop("Initialization")

    IF ( ( input%preconditioning_param /= 0 ) .AND. oneD%odi%d1 ) THEN
      CALL juDFT_error('Currently no preconditioner for 1D calculations', calledby = 'fleur')
    END IF

    IF (l_opti) CALL optional(mpi,atoms,sphhar,vacuum,dimension,&
                              stars,input,sym,cell,sliceplot,obsolete,xcpot,noco,oneD)

    IF (input%l_wann.AND.(mpi%irank==0).AND.(.NOT.wann%l_bs_comf)) THEN
       IF(mpi%isize.NE.1) CALL juDFT_error('No Wannier+MPI at the moment',calledby = 'fleur')
       CALL wann_optional(input,kpts,atoms,sym,cell,oneD,noco,wann)
    END IF
  
    iter     = 0
    iterHF   = 0
    l_cont = (iter < input%itmax)
    
    IF (mpi%irank.EQ.0) CALL openXMLElementNoAttributes('scfLoop')

    ! Initialize and load inDen density (start)
    CALL inDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
    
    archiveType = CDN_ARCHIVE_TYPE_CDN1_const
    IF (noco%l_noco) archiveType = CDN_ARCHIVE_TYPE_NOCO_const
    IF(mpi%irank.EQ.0) THEN
       CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                        0,results%ef,l_qfix,inDen)
       CALL timestart("Qfix")
       CALL qfix(mpi,stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.false.,fix)
       CALL timestop("Qfix")
       CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                         0,-1.0,results%ef,.FALSE.,inDen)
    END IF
    
    IF ((sliceplot%iplot.NE.0 ).AND.(mpi%irank==0) ) THEN
       CALL makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                      noco, inDen, PLOT_INPDEN, sliceplot) 
    END IF 

    ! Initialize and load inDen density (end)

    ! Initialize potentials (start)
    CALL vTot%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTTOT)
    CALL vCoul%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTCOUL)
    CALL vx%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTCOUL)
    CALL vTemp%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTTOT)
    ! Initialize potentials (end)

    ! Open/allocate eigenvector storage (start)
    l_real=sym%invs.AND..NOT.noco%l_noco
    eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,DIMENSION%neigd,kpts%nkpt,input%jspins,&
                    noco%l_noco,.NOT.INPUT%eig66(1),l_real,noco%l_soc,INPUT%eig66(1),mpi%n_size)

#ifdef CPP_CHASE
    CALL init_chase(mpi,dimension,input,atoms,kpts,noco,sym%invs.AND..NOT.noco%l_noco)
#endif
    ! Open/allocate eigenvector storage (end)

    scfloop:DO WHILE (l_cont)

       iter = iter + 1
       IF (mpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/),&
                                                       (/iter,inden%iter/), RESHAPE((/19,13,5,5/),(/2,2/)))

!!$       !+t3e
!!$       IF (input%alpha.LT.10.0) THEN
!!$
!!$          IF (iter.GT.1) THEN
!!$             input%alpha = input%alpha - NINT(input%alpha)
!!$          END IF

       !CALL resetIterationDependentTimers()
       CALL timestart("Iteration")
       IF (mpi%irank.EQ.0) THEN
          WRITE (6,FMT=8100) iter
8100      FORMAT (/,10x,'   iter=  ',i5)
       ENDIF !mpi%irank.eq.0
       input%total = .TRUE.

#ifdef CPP_CHASE
       CALL chase_distance(results%last_distance)
#endif

#ifdef CPP_MPI
       CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,inDen)
#endif

       dimension%neigd2 = dimension%neigd
       IF (noco%l_soc) dimension%neigd2 = dimension%neigd*2

       !HF
       !$ num_threads = omp_get_max_threads()
       !$ call omp_set_num_threads(1)
       IF (hybrid%l_hybrid) THEN
          SELECT TYPE(xcpot)
          TYPE IS(t_xcpot_inbuild)
             CALL calc_hybrid(eig_id,hybrid,kpts,atoms,input,DIMENSION,mpi,noco,&
                              cell,oneD,enpara,results,sym,xcpot,vTot,iter,iterHF)
          END SELECT
          IF(hybrid%l_calhf) THEN
             call mixing_history_reset(mpi)
             iter = 0
          END IF
       ENDIF
       !RDMFT
       IF(input%l_rdmft) THEN
          CALL open_hybrid_io1(DIMENSION,sym%invs)
       END IF
       IF(.not.input%eig66(1))THEN
          CALL reset_eig(eig_id,noco%l_soc) ! This has to be placed after the calc_hybrid call but before eigen
       END IF
       !$ call omp_set_num_threads(num_threads)

       !#endif

!!$             DO pc = 1, wann%nparampts
!!$                !---> gwf
!!$                IF (wann%l_sgwf.OR.wann%l_ms) THEN
!!$                   noco%qss(:) = wann%param_vec(:,pc)
!!$                   noco%alph(:) = wann%param_alpha(:,pc)
!!$                ELSE IF (wann%l_socgwf) THEN
!!$                   IF(wann%l_dim(2)) noco%phi   = tpi_const * wann%param_vec(2,pc)
!!$                   IF(wann%l_dim(3)) noco%theta = tpi_const * wann%param_vec(3,pc)
!!$                END IF
       !---< gwf

       CALL timestart("generation of potential")
       CALL vgen(hybrid,field,input,xcpot,DIMENSION,atoms,sphhar,stars,vacuum,sym,&
                 obsolete,cell,oneD,sliceplot,mpi,results,noco,EnergyDen,inDen,vTot,vx,vCoul)
       CALL timestop("generation of potential")

       IF ((sliceplot%iplot.NE.0 ).AND.(mpi%irank==0) ) THEN            
          CALL makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                         noco, vTot, PLOT_POT_TOT, sliceplot)      
!          CALL makeplots(sym,stars,vacuum,atoms,sphhar,input,cell,oneD,noco,sliceplot,vCoul,PLOT_POT_COU)
!          CALL subPotDen(vxcForPlotting,vTot,vCoul)
!          CALL makeplots(sym,stars,vacuum,atoms,sphhar,input,cell,oneD,noco,sliceplot,vxcForPlotting,PLOT_POT_VXC
       END IF 

#ifdef CPP_MPI
       CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

       CALL forcetheo%start(vtot,mpi%irank==0)
       forcetheoloop:DO WHILE(forcetheo%next_job(iter==input%itmax,atoms,noco))

          CALL timestart("gen. of hamil. and diag. (total)")
          CALL timestart("eigen")
          vTemp = vTot
          CALL timestart("Updating energy parameters")
          CALL enpara%update(mpi,atoms,vacuum,input,vToT)
          CALL timestop("Updating energy parameters")
          IF(.not.input%eig66(1))THEN
            CALL eigen(mpi,stars,sphhar,atoms,xcpot,sym,kpts,DIMENSION,vacuum,input,&
                     cell,enpara,banddos,noco,oneD,hybrid,iter,eig_id,results,inDen,vTemp,vx)
          ENDIF             
          vTot%mmpMat = vTemp%mmpMat
!!$          eig_idList(pc) = eig_id
          CALL timestop("eigen")

          ! add all contributions to total energy
#ifdef CPP_MPI
          ! send all result of local total energies to the r
          IF (hybrid%l_hybrid.AND.hybrid%l_calhf) THEN
             IF (mpi%irank==0) THEN
                CALL MPI_Reduce(MPI_IN_PLACE,results%te_hfex%core,1,MPI_REAL8,MPI_SUM,0,mpi%mpi_comm,ierr(1))
             ELSE
                CALL MPI_Reduce(results%te_hfex%core,MPI_IN_PLACE,1,MPI_REAL8,MPI_SUM,0, mpi%mpi_comm,ierr(1))
             END IF
             IF (mpi%irank==0) THEN
                CALL MPI_Reduce(MPI_IN_PLACE,results%te_hfex%valence,1,MPI_REAL8,MPI_SUM,0,mpi%mpi_comm,ierr(1))
             ELSE
                CALL MPI_Reduce(results%te_hfex%valence,MPI_IN_PLACE,1,MPI_REAL8,MPI_SUM,0, mpi%mpi_comm,ierr(1))
             END IF
          END IF
#endif

          ! WRITE(6,fmt='(A)') 'Starting 2nd variation ...'
          IF (noco%l_soc.AND..NOT.noco%l_noco.AND..NOT.INPUT%eig66(1)) &
             CALL eigenso(eig_id,mpi,DIMENSION,stars,vacuum,atoms,sphhar,&
                          obsolete,sym,cell,noco,input,kpts, oneD,vTot,enpara,results)
          CALL timestop("gen. of hamil. and diag. (total)")

#ifdef CPP_MPI
          CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

          ! fermi level and occupancies
          IF (noco%l_soc.AND.(.NOT.noco%l_noco)) DIMENSION%neigd = 2*DIMENSION%neigd

	  IF (input%gw.GT.0) THEN
	    IF (mpi%irank.EQ.0) THEN
	       CALL writeBasis(input,noco,kpts,atoms,sym,cell,enpara,vTot,vCoul,vx,mpi,DIMENSION,&
		  	     results,eig_id,oneD,sphhar,stars,vacuum)
	    END IF
	    IF (input%gw.EQ.2) THEN
	       CALL juDFT_end("GW data written. Fleur ends.",mpi%irank)
	    END IF
	  END IF

          !IF ((mpi%irank.EQ.0)) THEN
             CALL timestart("determination of fermi energy")

             IF (noco%l_soc.AND.(.NOT.noco%l_noco)) THEN
                input%zelec = input%zelec*2
                CALL fermie(eig_id,mpi,kpts,input,noco,enpara%epara_min,cell,results)
                results%seigscv = results%seigscv/2
                results%ts = results%ts/2
                input%zelec = input%zelec/2
             ELSE
                CALL fermie(eig_id,mpi,kpts,input,noco,enpara%epara_min,cell,results)
             ENDIF
             CALL timestop("determination of fermi energy")

!!$          !+Wannier
!!$          IF(wann%l_bs_comf)THEN
!!$             IF(pc.EQ.1) THEN
!!$                OPEN(777,file='out_eig.1')
!!$                OPEN(778,file='out_eig.2')
!!$                OPEN(779,file='out_eig.1_diag')
!!$                OPEN(780,file='out_eig.2_diag')
!!$             END IF
!!$
!!$             CALL bs_comfort(eig_id,DIMENSION,input,noco,kpts%nkpt,pc)
!!$
!!$             IF(pc.EQ.wann%nparampts)THEN
!!$                CLOSE(777)
!!$                CLOSE(778)
!!$                CLOSE(779)
!!$                CLOSE(780)
!!$             END IF
!!$          END IF
!!$          !-Wannier

          !ENDIF
#ifdef CPP_MPI
          CALL MPI_BCAST(results%ef,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(results%w_iks,SIZE(results%w_iks),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
#endif

          IF (forcetheo%eval(eig_id,DIMENSION,atoms,kpts,sym,cell,noco,input,mpi,oneD,enpara,vToT,results)) THEN
             IF (noco%l_soc.AND.(.NOT.noco%l_noco)) DIMENSION%neigd=DIMENSION%neigd/2
             CYCLE forcetheoloop
          ENDIF

          
          !+Wannier functions
          IF ((input%l_wann).AND.(.NOT.wann%l_bs_comf)) THEN
             CALL wannier(DIMENSION,mpi,input,kpts,sym,atoms,stars,vacuum,sphhar,oneD,&
                  wann,noco,cell,enpara,banddos,sliceplot,vTot,results,&
                  (/eig_id/),(sym%invs).AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco),kpts%nkpt)
          END IF
          !-Wannier

          ! charge density generation
          CALL timestart("generation of new charge density (total)")
          CALL outDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
          outDen%iter = inDen%iter
          CALL cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum, &
                      dimension,kpts,atoms,sphhar,stars,sym,&
                      enpara,cell,noco,vTot,results,oneD,coreSpecInput,&
                      archiveType,xcpot,outDen,EnergyDen)
           
          IF ((sliceplot%iplot.NE.0 ).AND.(mpi%irank==0) ) THEN        
!               CDN including core charge
                CALL makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, outDen, PLOT_OUTDEN_Y_CORE, sliceplot)
!!               CDN subtracted by core charge
!                CALL makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, &
!                               cell, noco, outDen, PLOT_OUTDEN_N_CORE, sliceplot)
          END IF 

          IF (input%l_rdmft) THEN
             SELECT TYPE(xcpot)
                TYPE IS(t_xcpot_inbuild)
                   CALL rdmft(eig_id,mpi,input,kpts,banddos,sliceplot,cell,atoms,enpara,stars,vacuum,dimension,&
                              sphhar,sym,field,vTot,vCoul,oneD,noco,xcpot,hybrid,results,coreSpecInput,archiveType,outDen)
             END SELECT
          END IF

          IF (noco%l_soc.AND.(.NOT.noco%l_noco)) DIMENSION%neigd=DIMENSION%neigd/2

#ifdef CPP_MPI
          CALL MPI_BCAST(enpara%evac,SIZE(enpara%evac),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(enpara%evac0,SIZE(enpara%evac0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(enpara%el0,SIZE(enpara%el0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(enpara%ello0,SIZE(enpara%ello0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)

          IF (noco%l_noco) THEN
             DO n= 1,atoms%ntype
                IF (noco%l_relax(n)) THEN
                   CALL MPI_BCAST(noco%alph(n),1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                   CALL MPI_BCAST(noco%beta(n),1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                ENDIF
             ENDDO
             IF (noco%l_constr) THEN
                CALL MPI_BCAST(noco%b_con,SIZE(noco%b_con),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
             ENDIF
          ENDIF
#endif
          CALL timestop("generation of new charge density (total)")

             
!!$             !----> output potential and potential difference
!!$             IF (obsolete%disp) THEN
!!$                reap = .FALSE.
!!$                input%total = .FALSE.
!!$                CALL timestart("generation of potential (total)")
!!$                CALL vgen(hybrid,reap,input,xcpot,DIMENSION, atoms,sphhar,stars,vacuum,sym,&
!!$                     obsolete,cell,oneD,sliceplot,mpi, results,noco,outDen,inDenRot,vTot,vx,vCoul)
!!$                CALL timestop("generation of potential (total)")
!!$
!!$                CALL potdis(stars,vacuum,atoms,sphhar, input,cell,sym)
!!$             END IF
             
             ! total energy
             CALL timestart('determination of total energy')
             CALL totale(mpi,atoms,sphhar,stars,vacuum,DIMENSION,sym,input,noco,cell,oneD,&
                         xcpot,hybrid,vTot,vCoul,iter,inDen,results)
             CALL timestop('determination of total energy')
          IF (hybrid%l_hybrid) CALL close_eig(eig_id)

       END DO forcetheoloop

       CALL forcetheo%postprocess()

       CALL enpara%mix(mpi,atoms,vacuum,input,vTot%mt(:,0,:,:),vtot%vacz)
       field2 = field

       ! mix input and output densities
       CALL mix_charge(field2,DIMENSION,mpi,(iter==input%itmax.OR.judft_was_argument("-mix_io")),&
            stars,atoms,sphhar,vacuum,input,&
            sym,cell,noco,oneD,archiveType,xcpot,iter,inDen,outDen,results)
       
       IF ((sliceplot%iplot.NE.0 ).AND.(mpi%irank==0) ) THEN        
!               CDN including core charge
                CALL makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, outDen, PLOT_MIXDEN_Y_CORE, sliceplot)
!!               CDN subtracted by core charge
!                CALL makeplots(sym,stars,vacuum,atoms,sphhar,input,cell,oneD,noco,sliceplot,inDen,PLOT_MIXDEN_N_CORE)
!                CALL makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, &
!                               cell, noco, outDen, PLOT_OUTDEN_N_CORE, sliceplot)
          END IF 


       IF(mpi%irank == 0) THEN
         WRITE (6,FMT=8130) iter
8130     FORMAT (/,5x,'******* it=',i3,'  is completed********',/,/)
         WRITE(*,*) "Iteration:",iter," Distance:",results%last_distance
         CALL timestop("Iteration")
       END IF ! mpi%irank.EQ.0
          
#ifdef CPP_MPI
       CALL MPI_BCAST(results%last_distance,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif
       CALL priv_geo_end(mpi)

       l_cont = .TRUE.
       IF (hybrid%l_hybrid) THEN
          IF(hybrid%l_calhf) THEN
             l_cont = l_cont.AND.(iterHF < input%itmax)
             l_cont = l_cont.AND.(input%mindistance<=results%last_distance)
             CALL check_time_for_next_iteration(iterHF,l_cont)
          ELSE
             l_cont = l_cont.AND.(iter < 50) ! Security stop for non-converging nested PBE calculations
          END IF
          IF (hybrid%l_subvxc) THEN
             results%te_hfex%valence = 0
          END IF
       ELSE
          l_cont = l_cont.AND.(iter < input%itmax)
          ! MetaGGAs need a at least 2 iterations
          l_cont = l_cont.AND.((input%mindistance<=results%last_distance).OR.input%l_f & 
                               .OR. (xcpot%exc_is_MetaGGA() .and. iter == 1))
          CALL check_time_for_next_iteration(iter,l_cont)
       END IF

       !CALL writeTimesXML()

       IF (mpi%irank.EQ.0) THEN
          IF (isCurrentXMLElement("iteration")) CALL closeXMLElement('iteration')
       END IF

  !Break SCF loop if Plots were generated in ongoing run (iplot=/=0).
!!       IF(sliceplot%iplot.NE.0) THEN
!!          CALL juDFT_end("Stopped self consistency loop after plots have been generated.")
!!       END IF


    END DO scfloop ! DO WHILE (l_cont)

    ! Test: Build a field, for which the theoretical divergence etc. are known and
    ! compare with the result of the routine.

    !CALL builddivtest(stars,atoms,sphhar,vacuum,sym,cell,1,testDen)
    !CALL makeBxc(stars,atoms,sphhar,vacuum,input,noco,vTot,testDen)
    CALL dummyDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
    CALL matrixsplit(stars, atoms, sphhar, vacuum, input, noco, 1.0, inDen, dummyDen, testDen(1), testDen(2), testDen(3))
    CALL checkplotinp()
    CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, .FALSE., .FALSE., 'testDen             ', testDen(1), testDen(1), testDen(2), testDen(3))
    !CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, .FALSE., .FALSE., 'testDeny            ', testDen(2))
    !CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, .FALSE., .FALSE., 'testDenz            ', testDen(3))
    CALL sourcefree(mpi,dimension,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,testDen)

    CALL add_usage_data("Iterations",iter)

    IF (mpi%irank.EQ.0) CALL closeXMLElement('scfLoop')

    CALL close_eig(eig_id)

    CALL juDFT_end("all done",mpi%irank)
    
   CONTAINS
      SUBROUTINE priv_geo_end(mpi)
         TYPE(t_mpi),INTENT(IN)::mpi
         LOGICAL :: l_exist
         !Check if a new input was generated
         INQUIRE (file='inp_new',exist=l_exist)
         IF (l_exist) THEN
            CALL juDFT_end(" GEO new inp created ! ",mpi%irank)
         END IF
         !check for inp.xml
         INQUIRE (file='inp_new.xml',exist=l_exist)
         IF (.NOT.l_exist) RETURN
         IF (mpi%irank==0) THEN
            CALL system('mv inp.xml inp_old.xml')
            CALL system('mv inp_new.xml inp.xml')
            INQUIRE (file='qfix',exist=l_exist)
            IF (l_exist) THEN
               OPEN(2,file='qfix')
               WRITE(2,*)"F"
               CLOSE(2)
               PRINT *,"qfix set to F"
            ENDIF
            CALL mixing_history_reset(mpi)
         ENDIF
         CALL juDFT_end(" GEO new inp.xml created ! ",mpi%irank)
      END SUBROUTINE priv_geo_end
    
   END SUBROUTINE fleur_execute
END MODULE m_fleur
