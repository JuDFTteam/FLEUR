!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_init
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleur_init(mpi,&
       input,field,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
       sliceplot,banddos,enpara,xcpot,results,kpts,hybrid,&
       oneD,coreSpecInput,wann)
    USE m_types
    USE m_fleurinput_read_xml
    USE m_fleurinput_mpi_bc
    USE m_judft
    USE m_juDFT_init
    USE m_init_wannier_defaults
    USE m_postprocessInput
    USE m_gen_map
    USE m_dwigner
    !USE m_gen_bz
    USE m_ylm
    USE m_InitParallelProcesses
    USE m_xmlOutput
    USE m_constants
    USE m_winpXML
    USE m_writeOutParameters
    USE m_setupMPI
    USE m_cdn_io
    USE m_fleur_info
    USE m_mixing_history
    USE m_checks
    USE m_writeOutHeader
    !USE m_fleur_init_old
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
    TYPE(t_enpara)   ,INTENT(OUT):: enpara
    CLASS(t_xcpot),ALLOCATABLE,INTENT(OUT):: xcpot
    TYPE(t_results)  ,INTENT(OUT):: results
    TYPE(t_kpts)     ,INTENT(OUT):: kpts
    TYPE(t_hybrid)   ,INTENT(OUT):: hybrid
    TYPE(t_oneD)     ,INTENT(OUT):: oneD
    TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput
    TYPE(t_wann)     ,INTENT(OUT):: wann
    CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT)::forcetheo

    type(t_enparaXML)::enparaXML
    TYPE(t_forcetheo_data)::forcetheo_data


    INTEGER, ALLOCATABLE          :: xmlElectronStates(:,:)
    INTEGER, ALLOCATABLE          :: atomTypeSpecies(:)
    INTEGER, ALLOCATABLE          :: speciesRepAtomType(:)
    REAL, ALLOCATABLE             :: xmlCoreOccs(:,:,:)
    LOGICAL, ALLOCATABLE          :: xmlPrintCoreStates(:,:)
    !     .. Local Scalars ..
    INTEGER    :: i,n,l,m1,m2,isym,iisym,numSpecies,pc,iAtom,iType
    COMPLEX    :: cdum
    CHARACTER(len=4)              :: namex
    CHARACTER(len=12)             :: relcor, tempNumberString
    CHARACTER(LEN=20)             :: filename
    REAL                          :: a1(3),a2(3),a3(3)
    REAL                          :: dtild, phi_add
    LOGICAL                       :: l_found, l_kpts, l_exist, l_krla

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER ierr(3)
    CALL MPI_COMM_RANK (mpi%mpi_comm,mpi%irank,ierr)
    CALL MPI_COMM_SIZE (mpi%mpi_comm,mpi%isize,ierr)
#else
    mpi%irank=0 ; mpi%isize=1; mpi%mpi_comm=1
#endif
    !determine if we use an xml-input file
    INQUIRE (file='inp.xml',exist=l_found)
    IF (.NOT.l_found) THEN
       CALL judft_error("No input file found",calledby='fleur_init',hint="To use FLEUR, you have to provide an 'inp.xml' file in the working directory")
    END IF

    CALL check_command_line()
#ifdef CPP_HDF
    CALL hdf_init()
#endif
    IF (mpi%irank.EQ.0) THEN
       CALL startFleur_XMLOutput()
       IF (judft_was_argument("-info")) THEN
          CLOSE(6)
          OPEN (6,status='SCRATCH')
       ELSE
          IF (.not.judft_was_argument("-no_out")) &
               OPEN (6,file='out',form='formatted',status='unknown')
       ENDIF
       CALL writeOutHeader()
       OPEN (16,status='SCRATCH')
    ENDIF

    ALLOCATE(t_xcpot_inbuild::xcpot)

    IF (mpi%irank.EQ.0) THEN
       CALL fleurinput_read_xml(cell,sym,atoms,input,noco,vacuum,field,&
            sliceplot,banddos,hybrid,oneD,coreSpecInput,wann,&
            xcpot,forcetheo_data,kpts,enparaXML)
    END IF

    CALL fleurinput_mpi_bc(cell,sym,atoms,input,noco,vacuum,field,&
         sliceplot,banddos,hybrid,oneD,coreSpecInput,wann,&
         xcpot,forcetheo_data,kpts,enparaXML,mpi%mpi_comm)


    CALL timestart("postprocessInput")
    CALL postprocessInput(mpi,input,field,sym,stars,atoms,vacuum,kpts,&
         oneD,hybrid,cell,banddos,sliceplot,xcpot,forcetheo,forcetheo_data,&
         noco,DIMENSION,enpara,enparaxml,sphhar,l_kpts)
    CALL timestop("postprocessInput")

    IF (mpi%irank.EQ.0) THEN
       CALL w_inpXML(&
            atoms,vacuum,input,stars,sliceplot,forcetheo,banddos,&
            cell,sym,xcpot,noco,oneD,hybrid,kpts,enpara,&
            .TRUE.,[.TRUE.,.TRUE.,.TRUE.,.TRUE.])
    END IF


    CALL ylmnorm_init(atoms%lmaxd)
    !
    !--> determine more dimensions
    !
    DIMENSION%nbasfcn = DIMENSION%nvd + atoms%nat*atoms%nlod*(2*atoms%llod+1)
    DIMENSION%lmd     = atoms%lmaxd* (atoms%lmaxd+2)
    IF (noco%l_noco) DIMENSION%nbasfcn = 2*DIMENSION%nbasfcn


    IF (mpi%irank.EQ.0) THEN
       CALL writeOutParameters(mpi,input,sym,stars,atoms,vacuum,kpts,&
            oneD,hybrid,cell,banddos,sliceplot,xcpot,&
            noco,DIMENSION,enpara,sphhar)
       CALL fleur_info(kpts)
       CALL deleteDensities()
    END IF

    !Finalize the MPI setup
    CALL setupMPI(kpts%nkpt,dimension%neigd,mpi)

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
    CALL add_usage_data("nkpt", kpts%nkpt)

#ifdef CPP_GPU
    CALL add_usage_data("gpu_per_node",1)
#else
    CALL add_usage_data("gpu_per_node",0)
#endif

    CALL results%init(DIMENSION,input,atoms,kpts,noco)

    IF (mpi%irank.EQ.0) THEN
       IF(input%gw.NE.0) CALL mixing_history_reset(mpi)
       CALL setStartingDensity(noco%l_noco)
    END IF

    !new check mode will only run the init-part of FLEUR
    IF (judft_was_argument("-check")) CALL judft_end("Check-mode done",mpi%irank)
  CONTAINS
    SUBROUTINE init_hybrid()
      IF (xcpot%is_hybrid().OR.input%l_rdmft) THEN
         IF (input%film.OR.oneD%odi%d1) THEN
            CALL juDFT_error("2D film and 1D calculations not implemented for HF/EXX/PBE0/HSE", &
                 calledby ="fleur", hint="Use a supercell or a different functional")
         END IF

         !             IF( ANY( atoms%l_geo  ) )&
         !                  &     CALL juDFT_error("Forces not implemented for HF/PBE0/HSE ",&
         !                  &                    calledby ="fleur")

         !calculate whole Brilloun zone
         !CALL gen_bz(kpts,sym)
         CALL gen_map(atoms,sym,oneD,hybrid)

         ! calculate d_wgn
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
         hybrid%l_calhf = .FALSE.
         ALLOCATE(hybrid%map(0,0),hybrid%tvec(0,0,0),hybrid%d_wgn2(0,0,0,0))
         IF(input%l_rdmft) THEN
            hybrid%l_calhf = .FALSE.
         END IF
      ENDIF
    END SUBROUTINE init_hybrid

    SUBROUTINE init_wannier()
      ! Initializations for Wannier functions (start)
      IF (mpi%irank.EQ.0) THEN
         wann%l_gwf = wann%l_ms.OR.wann%l_sgwf.OR.wann%l_socgwf

         IF(wann%l_gwf) THEN
            WRITE(*,*)'running HDWF-extension of FLEUR code'
            WRITE(*,*)'with l_sgwf =',wann%l_sgwf,' and l_socgwf =',wann%l_socgwf

            IF(wann%l_socgwf.AND. .NOT.noco%l_soc) THEN
               CALL juDFT_error("set l_soc=T if l_socgwf=T",calledby="fleur_init")
            END IF

            IF((wann%l_ms.OR.wann%l_sgwf).AND..NOT.(noco%l_noco.AND.noco%l_ss)) THEN
               CALL juDFT_error("set l_noco=l_ss=T for l_sgwf.or.l_ms",calledby="fleur_init")
            END IF

            IF((wann%l_ms.OR.wann%l_sgwf).AND.wann%l_socgwf) THEN
               CALL juDFT_error("(l_ms.or.l_sgwf).and.l_socgwf",calledby="fleur_init")
            END IF

            INQUIRE(FILE=wann%param_file,EXIST=l_exist)
            IF(.NOT.l_exist) THEN
               CALL juDFT_error("where is param_file"//TRIM(wann%param_file)//"?",calledby="fleur_init")
            END IF
            OPEN (113,file=wann%param_file,status='old')
            READ (113,*) wann%nparampts,wann%scale_param
            CLOSE(113)
         ELSE
            wann%nparampts=1
            wann%scale_param=1.0
         END IF
      END IF

      ALLOCATE (wann%param_vec(3,wann%nparampts))
      ALLOCATE (wann%param_alpha(atoms%ntype,wann%nparampts))

      IF(mpi%irank.EQ.0) THEN
         IF(wann%l_gwf) THEN
            OPEN(113,file=wann%param_file,status='old')
            READ(113,*)!header
            WRITE(6,*) 'parameter points for HDWFs generation:'
            IF(wann%l_sgwf.OR.wann%l_ms) THEN
               WRITE(6,*)'      q1       ','      q2       ','      q3'
            ELSE IF(wann%l_socgwf) THEN
               WRITE(6,*)'      --       ','     phi       ','    theta'
            END IF

            DO pc = 1, wann%nparampts
               READ(113,'(3(f14.10,1x))') wann%param_vec(1,pc), wann%param_vec(2,pc), wann%param_vec(3,pc)
               wann%param_vec(:,pc) = wann%param_vec(:,pc) / wann%scale_param
               WRITE(6,'(3(f14.10,1x))') wann%param_vec(1,pc), wann%param_vec(2,pc), wann%param_vec(3,pc)
               IF(wann%l_sgwf.OR.wann%l_ms) THEN
                  iAtom = 1
                  DO iType = 1, atoms%ntype
                     phi_add = tpi_const*(wann%param_vec(1,pc)*atoms%taual(1,iAtom) +&
                          wann%param_vec(2,pc)*atoms%taual(2,iAtom) +&
                          wann%param_vec(3,pc)*atoms%taual(3,iAtom))
                     wann%param_alpha(iType,pc) = noco%alph(iType) + phi_add
                     iAtom = iAtom + atoms%neq(iType)
                  END DO
               END IF
            END DO

            IF(ANY(wann%param_vec(1,:).NE.wann%param_vec(1,1))) wann%l_dim(1)=.TRUE.
            IF(ANY(wann%param_vec(2,:).NE.wann%param_vec(2,1))) wann%l_dim(2)=.TRUE.
            IF(ANY(wann%param_vec(3,:).NE.wann%param_vec(3,1))) wann%l_dim(3)=.TRUE.

            CLOSE(113)

            IF(wann%l_dim(1).AND.wann%l_socgwf) THEN
               CALL juDFT_error("do not specify 1st component if l_socgwf",calledby="fleur_init")
            END IF
         END IF!(wann%l_gwf)
      END IF!(mpi%irank.EQ.0)

#ifdef CPP_MPI
      CALL MPI_BCAST(wann%param_vec,3*wann%nparampts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(wann%param_alpha,atoms%ntype*wann%nparampts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(wann%l_dim,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

      ! Initializations for Wannier functions (end)
    END SUBROUTINE init_wannier
  END SUBROUTINE fleur_init
END MODULE m_fleur_init
