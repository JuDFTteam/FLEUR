
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
    USE m_broyd_io
    USE m_qfix
    USE m_vgen
    USE m_rhodirgen
    USE m_writexcstuff
    USE m_vmatgen
    USE m_eigen
    USE m_eigenso
    USE m_fermie
    USE m_force0
    USE m_cdngen
    USE m_totale
    USE m_potdis
    USE m_mix
    USE m_xmlOutput
    USE m_juDFT_time
    USE m_calc_hybrid
    USE m_wann_optional
    USE m_wannier
    USE m_bs_comfort
    USE m_gen_map
    USE m_dwigner
    USE m_ylm
#ifdef CPP_MPI
    USE m_mpi_bc_all,  ONLY : mpi_bc_all
    USE m_mpi_bc_potden
#endif
    USE m_eig66_io,   ONLY : open_eig, close_eig
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: mpi_comm

    !     Types, these variables contain a lot of data!
    TYPE(t_input)    :: input
    TYPE(t_dimension):: DIMENSION
    TYPE(t_atoms)    :: atoms
    TYPE(t_sphhar)   :: sphhar
    TYPE(t_cell)     :: cell
    TYPE(t_stars)    :: stars
    TYPE(t_sym)      :: sym
    TYPE(t_noco)     :: noco
    TYPE(t_vacuum)   :: vacuum
    TYPE(t_sliceplot):: sliceplot
    TYPE(t_banddos)  :: banddos
    TYPE(t_obsolete) :: obsolete
    TYPE(t_enpara)   :: enpara
    TYPE(t_xcpot)    :: xcpot
    TYPE(t_results)  :: results
    TYPE(t_kpts)     :: kpts
    TYPE(t_hybrid)   :: hybrid
    TYPE(t_oneD)     :: oneD
    TYPE(t_mpi)      :: mpi
    TYPE(t_coreSpecInput) :: coreSpecInput
    TYPE(t_wann)     :: wann
    TYPE(t_potden)   :: vTot,vx,vCoul,vTemp
    TYPE(t_potden)   :: inDen, outDen, inDenRot
    CLASS(t_forcetheo),ALLOCATABLE:: forcetheo

    !     .. Local Scalars ..
    INTEGER:: eig_id, archiveType
    INTEGER:: n,it,ithf
    LOGICAL:: reap,l_opti,l_cont,l_qfix, l_wann_inp
    REAL   :: fermiEnergyTemp, fix
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER:: ierr(2)
#endif

    mpi%mpi_comm = mpi_comm
 
    CALL timestart("Initialization")
    CALL fleur_init(mpi,input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
         sliceplot,banddos,obsolete,enpara,xcpot,results,kpts,hybrid,&
         oneD,coreSpecInput,wann,l_opti)
    CALL timestop("Initialization")



    IF (l_opti) &
         CALL OPTIONAL(mpi,atoms,sphhar,vacuum,DIMENSION,&
         stars,input,sym,cell,sliceplot,obsolete,xcpot,noco,oneD)
 

    !+Wannier
    INQUIRE (file='wann_inp',exist=l_wann_inp)
    input%l_wann = input%l_wann.OR.l_wann_inp
    IF (input%l_wann.AND.(mpi%irank==0).AND.(.NOT.wann%l_bs_comf)) THEN
       IF(mpi%isize.NE.1) CALL juDFT_error('No Wannier+MPI at the moment',calledby = 'fleur')
       CALL wann_optional(input,kpts,atoms,sym,cell,oneD,noco,wann)
    END IF
    IF (wann%l_gwf) input%itmax = 1

    !-Wannier


    it     = 0
    ithf   = 0
    l_cont = ( it < input%itmax )
    
    results%last_distance = -1.0
    IF (mpi%irank.EQ.0) CALL openXMLElementNoAttributes('scfLoop')

    ! Initialize and load inDen density (start)
    CALL inDen%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)
    CALL inDenRot%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)

    archiveType = CDN_ARCHIVE_TYPE_CDN1_const
    IF (noco%l_noco) archiveType = CDN_ARCHIVE_TYPE_NOCO_const

    IF(mpi%irank.EQ.0) THEN
       CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
            0,fermiEnergyTemp,l_qfix,inDen)
       CALL timestart("Qfix")
       CALL qfix(stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.false.,fix)
       CALL timestop("Qfix")
       CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
            0,-1.0,0.0,.FALSE.,inDen)
    END IF
    ! Initialize and load inDen density (end)

    ! Initialize potentials (start)
    CALL vTot%init(stars,atoms,sphhar,vacuum,noco,oneD,DIMENSION%jspd,noco%l_noco,POTDEN_TYPE_POTTOT)
    CALL vCoul%init(stars,atoms,sphhar,vacuum,noco,oneD,DIMENSION%jspd,noco%l_noco,POTDEN_TYPE_POTCOUL)
    CALL vx%init(stars,atoms,sphhar,vacuum,noco,oneD,DIMENSION%jspd,.FALSE.,POTDEN_TYPE_POTCOUL)
    CALL vTemp%init(stars,atoms,sphhar,vacuum,noco,oneD,DIMENSION%jspd,noco%l_noco,POTDEN_TYPE_POTTOT)
    ! Initialize potentials (end)

    scfloop:DO WHILE (l_cont)

       it = it + 1
       IF (mpi%irank.EQ.0) CALL openXMLElementFormPoly('iteration',(/'numberForCurrentRun','overallNumber      '/)&
            ,(/it,inden%iter/), RESHAPE((/19,13,5,5/),(/2,2/)))
       inDenRot = inDen

!!$       !+t3e
!!$       IF (input%alpha.LT.10.0) THEN
!!$
!!$          IF (it.GT.1) THEN
!!$             input%alpha = input%alpha - NINT(input%alpha)
!!$          END IF

       CALL resetIterationDependentTimers()
       CALL timestart("Iteration")
       IF (mpi%irank.EQ.0) THEN
          !-t3e
          WRITE (6,FMT=8100) it
          WRITE (16,FMT=8100) it
8100      FORMAT (/,10x,'   it=    ',i5)


          !      ----> potential generator
          !
          !---> pk non-collinear
          !--->        reload the density matrix from file rhomat_in
          !--->        calculate spin-up and -down density for USE in the
          !--->        potential generator and store the direction of
          !--->        magnetization on file dirofmag
          IF (noco%l_noco) THEN
             CALL timestart("gen. spin-up and -down density")
             CALL rhodirgen(DIMENSION,sym,stars,atoms,sphhar,&
                  vacuum,cell,input,noco,oneD,inDenRot)
             CALL timestop("gen. spin-up and -down density")
          ENDIF
          !---> pk non-collinear

          reap=.NOT.obsolete%disp
          input%total = .TRUE.
       ENDIF !mpi%irank.eq.0
#ifdef CPP_MPI
       CALL MPI_BCAST(input%total,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif

#ifdef CPP_MPI
       CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,inDen)
       CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,inDenRot)
#endif


       IF ( noco%l_soc ) THEN
          dimension%neigd2 = dimension%neigd*2
       ELSE
          dimension%neigd2 = dimension%neigd
       END IF


       !HF
       IF (hybrid%l_hybrid) CALL  calc_hybrid(hybrid,kpts,atoms,input,DIMENSION,mpi,noco,&
            cell,vacuum,oneD,banddos,results,sym,xcpot,vTot,it)
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
       CALL vgen(hybrid,reap,input,xcpot,DIMENSION, atoms,sphhar,stars,vacuum,&
            sym,obsolete,cell, oneD,sliceplot,mpi ,results,noco,inDen,inDenRot,vTot,vx,vCoul)
       CALL timestop("generation of potential")

       IF (mpi%irank.EQ.0) THEN
          !---> pk non-collinear
          !--->          generate the four component matrix potential from spin up
          !--->          and down potentials and direction of the magnetic field
          IF (noco%l_noco) THEN
             CALL timestart("generation of potential-matrix")
             CALL vmatgen(stars, atoms,sphhar,vacuum,sym,input,oneD,inDenRot,vTot)
             CALL timestop("generation of potential-matrix")
          ENDIF
          !
       ENDIF ! mpi%irank.eq.0



#ifdef CPP_MPI
       CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

       CALL forcetheo%start()

       forcetheoloop:DO WHILE(forcetheo%next_job(it==input%itmax,noco))


          CALL timestart("generation of hamiltonian and diagonalization (total)")
          CALL timestart("eigen")
          vTemp = vTot
          CALL enpara%update(mpi,atoms,vacuum,input,vToT)
          CALL eigen(mpi,stars,sphhar,atoms,obsolete,xcpot,&
               sym,kpts,DIMENSION,vacuum,input,cell,enpara,banddos,noco,oneD,hybrid,&
               it,eig_id,results,inDenRot,vTemp,vx)
          vTot%mmpMat = vTemp%mmpMat
!!$          eig_idList(pc) = eig_id
          CALL timestop("eigen")
          !
          !                   add all contributions to total energy
          !
#ifdef CPP_MPI
          ! send all result of local total energies to the r
          IF (mpi%irank==0) THEN
             CALL MPI_Reduce(MPI_IN_PLACE,results%te_hfex%valence,&
                  1,MPI_REAL8,MPI_SUM,0,mpi%mpi_comm,ierr(1))
             CALL MPI_Reduce(MPI_IN_PLACE,results%te_hfex%core,&
                  1,MPI_REAL8,MPI_SUM,0,mpi%mpi_comm,ierr(1))
          ELSE
             CALL MPI_Reduce(results%te_hfex%valence,MPI_IN_PLACE,&
                  1,MPI_REAL8,MPI_SUM,0, mpi%mpi_comm,ierr(1))
             CALL MPI_Reduce(results%te_hfex%core,MPI_IN_PLACE,&
                  1,MPI_REAL8,MPI_SUM,0, mpi%mpi_comm,ierr(1))
          ENDIF
          !                                  END IF
#endif



          ! WRITE(6,fmt='(A)') 'Starting 2nd variation ...'
          IF (noco%l_soc.AND..NOT.noco%l_noco) &
               CALL eigenso(eig_id,mpi,DIMENSION,stars,vacuum,atoms,sphhar,&
               obsolete,sym,cell,noco,input,kpts, oneD,vTot,enpara)
          CALL timestop("generation of hamiltonian and diagonalization (total)")

#ifdef CPP_MPI
          CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

          !
          !              ----> fermi level and occupancies

          IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) DIMENSION%neigd = 2*DIMENSION%neigd
          IF( .NOT. ALLOCATED(results%w_iks) )&
               ALLOCATE ( results%w_iks(DIMENSION%neigd,kpts%nkpt,DIMENSION%jspd) )
          IF ( (mpi%irank.EQ.0)) THEN
             CALL timestart("determination of fermi energy")

             IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) THEN
                input%zelec = input%zelec*2
                CALL fermie(eig_id,mpi,kpts,obsolete,&
                     input,noco,enpara%epara_min,cell,results)
                results%seigscv = results%seigscv/2
                results%ts = results%ts/2
                input%zelec = input%zelec/2
             ELSE
                CALL fermie(eig_id,mpi,kpts,obsolete,&
                     input,noco,enpara%epara_min,cell,results)
             ENDIF
             CALL timestop("determination of fermi energy")
!!$             
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
          ENDIF
#ifdef CPP_MPI
          CALL MPI_BCAST(results%ef,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(results%w_iks,SIZE(results%w_iks),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
#endif


          IF (forcetheo%eval(eig_id,DIMENSION,atoms,kpts,sym,&
               cell,noco, input,mpi, oneD,enpara,vToT,results)) THEN
             IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) DIMENSION%neigd=DIMENSION%neigd/2
             CYCLE forcetheoloop
          ENDIF


          CALL force_0(results)!        ----> initialise force_old
          !        ----> charge density
          
!!$          !+Wannier functions
!!$          IF ((input%l_wann).AND.(.NOT.wann%l_bs_comf)) THEN
!!$             CALL wannier(DIMENSION,mpi,input,kpts,sym,atoms,stars,vacuum,sphhar,oneD,&
!!$                  wann,noco,cell,enpara,banddos,sliceplot,vTot,results,&
!!$                  eig_idList,(sym%invs).AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco),kpts%nkpt)
!!$          END IF
!!$          IF (wann%l_gwf) CALL juDFT_error("provide wann_inp if l_gwf=T", calledby = "fleur")
!!$          !-Wannier

          CALL timestart("generation of new charge density (total)")
          CALL outDen%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)
          outDen%iter = inDen%iter
          CALL cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,&
               DIMENSION,kpts,atoms,sphhar,stars,sym,obsolete,&
               enpara,cell,noco,vTot,results,oneD,coreSpecInput,&
               archiveType,outDen)

          IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) DIMENSION%neigd=DIMENSION%neigd/2
          !+t3e
#ifdef CPP_MPI
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
          IF (mpi%irank.EQ.0) THEN
             !-t3e
             
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
             
             !----> total energy
             CALL timestart('determination of total energy')
             CALL totale(atoms,sphhar,stars,vacuum,DIMENSION,&
                  sym,input,noco,cell,oneD,xcpot,hybrid,vTot,vCoul,it,inDen,results)

             CALL timestop('determination of total energy')
          ENDIF ! mpi%irank.EQ.0
          IF ( hybrid%l_hybrid ) CALL close_eig(eig_id)

       END DO forcetheoloop

       CALL forcetheo%postprocess()
       
       CALL enpara%mix(mpi,atoms,vacuum,input,vTot%mt(:,0,:,:),vtot%vacz)
       IF (mpi%irank.EQ.0) THEN
          !          ----> mix input and output densities
          CALL timestart("mixing")
          CALL mix(stars,atoms,sphhar,vacuum,input,sym,cell,noco,oneD,hybrid,archiveType,inDen,outDen,results)
          CALL timestop("mixing")
          
          WRITE (6,FMT=8130) it
          WRITE (16,FMT=8130) it
8130      FORMAT (/,5x,'******* it=',i3,'  is completed********',/,/)
          WRITE(*,*) "Iteration:",it," Distance:",results%last_distance
          CALL timestop("Iteration")
          !+t3e
       ENDIF ! mpi%irank.EQ.0
       
          
#ifdef CPP_MPI
       CALL MPI_BCAST(results%last_distance,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
       CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif
       CALL priv_geo_end(mpi)

!!$       IF ( hybrid%l_calhf ) ithf = ithf + 1
!!$    IF ( hybrid%l_subvxc ) THEN
!!$       l_cont = ( ithf < input%itmax )
!!$       results%te_hfex%core    = 0
!!$       results%te_hfex%valence = 0
!!$    ELSE
       l_cont = ( it < input%itmax )
!!$    END IF
       CALL writeTimesXML()
       CALL check_time_for_next_iteration(it,l_cont)

       l_cont = l_cont.AND.((input%mindistance<=results%last_distance).OR.input%l_f)

       IF ((mpi%irank.EQ.0).AND.(isCurrentXMLElement("iteration"))) THEN
          CALL closeXMLElement('iteration')
       END IF

    END DO scfloop ! DO WHILE (l_cont)

    IF (mpi%irank.EQ.0) CALL closeXMLElement('scfLoop')
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
         CALL resetBroydenHistory()
      ENDIF
      CALL juDFT_end(" GEO new inp.xml created ! ",mpi%irank)
    END SUBROUTINE priv_geo_end
    
  END SUBROUTINE fleur_execute
END MODULE m_fleur
