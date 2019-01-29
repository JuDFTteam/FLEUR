!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_setup
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleur_setup(mpi,job,&
       input,field,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
       sliceplot,banddos,obsolete,enpara,xcpot,results,kpts,hybrid,&
       oneD,coreSpecInput,wann,l_opti)
    USE m_types
    USE m_judft
    USE m_juDFT_init
    USE m_init_wannier_defaults
    USE m_rinpXML
    USE m_postprocessInput
    USE m_gen_map
    USE m_dwigner
    USE m_gen_bz
    USE m_ylm
    USE m_InitParallelProcesses
    USE m_xmlOutput
    USE m_constants
    USE m_winpXML
    USE m_writeOutParameters
    USE m_setupMPI
    USE m_cdn_io
    USE m_fleur_info
    USE m_checks
    USE m_writeOutHeader
    USE m_fleur_init_old
    USE m_types_xcpot_inbuild
    USE m_mpi_bc_xcpot

#ifdef CPP_MPI
    USE m_mpi_bc_all,  ONLY : mpi_bc_all
#ifndef CPP_OLDINTEL
    USE m_mpi_dist_forcetheorem
#endif
#endif
#ifdef CPP_HDF
    USE m_hdf_tools
#endif
    IMPLICIT NONE
    !     Types, these variables contain a lot of data!
    TYPE(t_mpi)    ,INTENT(INOUT):: mpi
    TYPE(t_job)      ,INTENT(OUT):: job
    TYPE(t_input)    ,INTENT(OUT):: input
    TYPE(t_field),    INTENT(OUT) :: field
    TYPE(t_dimension),INTENT(OUT):: DIMENSION
    TYPE(t_atoms)    ,INTENT(OUT):: atoms
    TYPE(t_sphhar)   ,INTENT(OUT):: sphhar
    TYPE(t_cell)     ,INTENT(OUT):: cell
    TYPE(t_stars)    ,INTENT(OUT):: stars
    TYPE(t_sym)      ,INTENT(OUT):: sym
    TYPE(t_noco)     ,INTENT(OUT):: noco
    TYPE(t_vacuum)   ,INTENT(OUT):: vacuum
    TYPE(t_sliceplot),INTENT(OUT):: sliceplot
    TYPE(t_banddos)  ,INTENT(OUT):: banddos
    TYPE(t_obsolete) ,INTENT(OUT):: obsolete 
    TYPE(t_enpara)   ,INTENT(OUT):: enpara
    CLASS(t_xcpot),ALLOCATABLE,INTENT(OUT):: xcpot
    TYPE(t_results)  ,INTENT(OUT):: results
    TYPE(t_kpts)     ,INTENT(OUT):: kpts
    TYPE(t_hybrid)   ,INTENT(OUT):: hybrid
    TYPE(t_oneD)     ,INTENT(OUT):: oneD
    TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput
    TYPE(t_wann)     ,INTENT(OUT):: wann
    CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT)::forcetheo
    LOGICAL,          INTENT(OUT):: l_opti


    INTEGER, ALLOCATABLE          :: xmlElectronStates(:,:)
    INTEGER, ALLOCATABLE          :: atomTypeSpecies(:)
    INTEGER, ALLOCATABLE          :: speciesRepAtomType(:)
    REAL, ALLOCATABLE             :: xmlCoreOccs(:,:,:)
    LOGICAL, ALLOCATABLE          :: xmlPrintCoreStates(:,:)
    CHARACTER(len=3), ALLOCATABLE :: noel(:)
    !     .. Local Scalars ..
    INTEGER    :: i,n,l,m1,m2,isym,iisym,numSpecies,pc,iAtom,iType
    COMPLEX    :: cdum
    CHARACTER(len=4)              :: namex
    CHARACTER(len=12)             :: relcor, tempNumberString
    CHARACTER(LEN=20)             :: filename
    REAL                          :: a1(3),a2(3),a3(3)
    REAL                          :: dtild, phi_add
    LOGICAL                       :: l_found, l_kpts, l_exist, l_krla


    CALL priv_setup_libraries()


    CALL check_command_line()
    CALL priv_setup_output()




    !Now read in the inp.xml file
    IF (mpi%irank==0) THEN
       CALL inp_open()
       !The setup types
       CALL cell%read_xml()
       CALL sym%read_xml()
       CALL noco%read_xml()
       CALL input%read_xml()
       CALL atoms%read_xml()
       CALL sphhar%read_xml()
       CALL stars%read_xml()
       CALL noco%read_xml()
       CALL vacuum%read_xml()
       CALL sliceplot%read_xml()
       CALL banddos%read_xml()
       CALL coreSpecInput%read_xml()
       CALL field%read_xml()
       !these features are not really ready yet...
       CALL wannier%read_xml()
       
       CALL hybrid%read_xml()!fehlt
       CALL kpts%read_xml()!?fehlt
       call enpara%read_xml() !fehlt
       
       !must be different...
       CALL xcpot%read_xml()!felht
       CALL foretheo%read_xml()!felht

       
       !The job type
       CALL job%read_xml()
       

       CALL writeOutParameters(mpi,input,sym,stars,atoms,vacuum,obsolete,kpts,&
            oneD,hybrid,cell,banddos,sliceplot,xcpot,&
            noco,dimension,enpara,sphhar)
       CALL fleur_info(kpts)
       CALL deleteDensities()

       
    ENDIF
    !check mode will only run the init-part of FLEUR
    IF (judft_was_argument("-check")) CALL judft_end("Check-mode done",mpi%irank)

    CALL oneD%defaults() !fehlt
   
 
    CALL ylmnorm_init(atoms%lmaxd)
    !
    !--> determine dimensions should be removed later on
    !
    DIMENSION%nbasfcn = DIMENSION%nvd + atoms%nat*atoms%nlod*(2*atoms%llod+1)
    DIMENSION%neigd = minNeigd
    DIMENSION%neigd = MAX(5,NINT(0.75*input%zelec))
    IF (noco%l_noco) DIMENSION%nbasfcn = 2*DIMENSION%nbasfcn

    CALL priv_setup_hybrid()

    !Finalize the MPI setup
    CALL setupMPI(kpts%nkpt,mpi)

    call collect_usage_data()
    CALL results%init(dimension,input,atoms,kpts,noco)

    IF (mpi%irank.EQ.0) THEN
       CALL setStartingDensity(noco%l_noco)
    END IF


 

  

  
    IF (input%l_inpXML) THEN            
       ALLOCATE(noel(1))
       IF (mpi%irank.EQ.0) THEN
          WRITE (6,*) 'XML code path used: Calculation parameters are stored in out.xml'
          ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
          ALLOCATE(atomTypeSpecies(1),speciesRepAtomType(1))
          ALLOCATE(xmlElectronStates(1,1),xmlPrintCoreStates(1,1))
          ALLOCATE(xmlCoreOccs(1,1,1))
          namex = '    '
          relcor = '            '
          a1 = 0.0
          a2 = 0.0
          a3 = 0.0
          CALL r_inpXML(&
               job,atoms,obsolete,vacuum,input,stars,sliceplot,banddos,DIMENSION,forcetheo,&
               cell,sym,xcpot,noco,oneD,hybrid,kpts,enpara,coreSpecInput,wann,&
               noel,namex,relcor,a1,a2,a3,dtild,xmlElectronStates,&
               xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType,&
               l_kpts)
       END IF
       CALL mpi_bc_xcpot(xcpot,mpi)
       CALL postprocessInput(mpi,job,input,field,sym,stars,atoms,vacuum,obsolete,kpts,&
            oneD,hybrid,cell,banddos,sliceplot,xcpot,forcetheo,&
            noco,dimension,enpara,sphhar,l_opti,noel,l_kpts)

       IF (mpi%irank.EQ.0) THEN
          filename = ''
          numSpecies = SIZE(speciesRepAtomType)
          CALL w_inpXML(&
               atoms,obsolete,vacuum,input,stars,sliceplot,forcetheo,banddos,&
               cell,sym,xcpot,noco,oneD,hybrid,kpts,job,kpts%nkpt3,kpts%l_gamma,&
               noel,namex,relcor,a1,a2,a3,dtild,input%comment,&
               xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
               atomTypeSpecies,speciesRepAtomType,.TRUE.,filename,&
               .TRUE.,numSpecies,enpara)

          DEALLOCATE(atomTypeSpecies,speciesRepAtomType)
          DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
       END IF

       DEALLOCATE(noel)

#ifdef CPP_MPI
       CALL initParallelProcesses(atoms,vacuum,input,stars,sliceplot,banddos,&
            DIMENSION,cell,sym,xcpot,noco,oneD,hybrid,&
            kpts,enpara,sphhar,mpi,obsolete)
#ifndef CPP_OLDINTEL
       CALL mpi_dist_forcetheorem(mpi,forcetheo)
#endif
#endif

    END IF ! end of else branch of "IF (input%l_inpXML) THEN"
    !

    
  CONTAINS

    SUBROUTINE priv_setup_hybrid()
 
      IF (    xcpot%is_hybrid() ) THEN
         
         IF (input%film .OR. oneD%odi%d1)&
              &    CALL juDFT_error("2D film and 1D calculations not implemented"&
              &                 //"for HF/EXX/PBE0/HSE", calledby ="fleur",&
            &                 hint="Use a supercell or a different functional")
         
         !             IF( ANY( atoms%l_geo  ) )&
       !                  &     CALL juDFT_error("Forces not implemented for HF/PBE0/HSE ",&
       !                  &                    calledby ="fleur")

       !calculate whole Brilloun zone
       !CALL gen_bz(kpts,sym)
       CALL gen_map(&
            &          atoms,sym,oneD,hybrid)
       !
       ! calculate d_wgn
       !
       ALLOCATE (hybrid%d_wgn2(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,sym%nsym))
       CALL d_wigner(sym%nop,sym%mrot,cell%bmat,atoms%lmaxd,hybrid%d_wgn2(:,:,1:,:sym%nop))
       hybrid%d_wgn2(:,:,0,:) = 1

       DO isym = sym%nop+1,sym%nsym
          iisym = isym - sym%nop
          DO l = 0,atoms%lmaxd
             DO m2 = -l,l
                DO m1 = -l,-1
                   cdum                  = hybrid%d_wgn2( m1,m2,l,iisym)
                   hybrid%d_wgn2( m1,m2,l,isym) = hybrid%d_wgn2(-m1,m2,l,iisym)*(-1)**m1
                   hybrid%d_wgn2(-m1,m2,l,isym) = cdum                  *(-1)**m1
                END DO
                hybrid%d_wgn2(0,m2,l,isym) = hybrid%d_wgn2(0,m2,l,iisym)
             END DO
          END DO
       END DO
    ELSE
       ALLOCATE(hybrid%map(0,0),hybrid%tvec(0,0,0),hybrid%d_wgn2(0,0,0,0))
       hybrid%l_calhf   = .FALSE.
    END IF
  END SUBROUTINE priv_setup_hybrid
    SUBROUTINE priv_setup_libraries()
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER ierr(3)
    CALL MPI_COMM_RANK (mpi%mpi_comm,mpi%irank,ierr)
    CALL MPI_COMM_SIZE (mpi%mpi_comm,mpi%isize,ierr)
#else
    mpi%irank=0 ; mpi%isize=1; mpi%mpi_comm=1
#endif
#ifdef CPP_HDF
    CALL hdf_init()
#endif
  END SUBROUTINE priv_setup_libraries
   
  SUBROUTINE priv_setup_output()
    IF (mpi%irank.EQ.0) THEN
       CALL startXMLOutput()
       IF (judft_was_argument("-info")) THEN
          CLOSE(6)
          OPEN (6,status='SCRATCH')
       ELSE
          IF (.not.judft_was_argument("-no_out")) &
               OPEN (6,file='out',form='formatted',status='unknown')
       ENDIF
       CALL writeOutHeader()
       IF (judft_was_argument("-info")) THEN
          OPEN (16,status='SCRATCH')
       ELSE
          OPEN (16,file='inf',form='formatted',status='unknown')
       ENDIF
    ENDIF
    END SUBROUTINE priv_setup_output

    SUBROUTINE collect_usage_data()
    !Collect some usage info
    CALL add_usage_data("A-Types",atoms%ntype)
    CALL add_usage_data("Atoms",atoms%nat)
    CALL add_usage_data("Real",sym%invs.AND..NOT.noco%l_noco)
    CALL add_usage_data("Spins",input%jspins)
    CALL add_usage_data("Noco",noco%l_noco)
    CALL add_usage_data("SOC",noco%l_soc)
    CALL add_usage_data("SpinSpiral",noco%l_ss)
    CALL add_usage_data("PlaneWaves",DIMENSION%nvd)
    CALL add_usage_data("LOs",atoms%nlotot)
    CALL add_usage_data("Iterations",job%itmax)
  END SUBROUTINE collect_usage_data

    
  END SUBROUTINE fleur_setup
      END MODULE m_fleur_init
