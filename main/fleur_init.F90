!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_init
#ifdef CPP_MPI 
   use mpi 
#endif
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleur_init(fmpi,&
       input,field,atoms, sphhar,cell,stars,sym,noco,nococonv,vacuum,forcetheo,&
       sliceplot,banddos,enpara,xcpot,results,kpts,mpinp,hybinp,&
       oneD,coreSpecInput,gfinp,hub1inp,wann)
    USE m_types
    USE m_fleurinput_read_xml
    USE m_fleurinput_mpi_bc
    USE m_types_mpinp
    USE m_judft
    USE m_juDFT_init
    USE m_init_wannier_defaults
    USE m_dwigner
    USE m_ylm
    !USE m_InitParallelProcesses
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
    USE m_prpxcfft
    USE m_make_stars
    USE m_make_sphhar
    USE m_convn
    USE m_efield
    USE m_fleurinput_postprocess
    USE m_make_forcetheo
    USE m_lapwdim
    USE m_gaunt, ONLY: gaunt_init
#ifdef CPP_MPI
    !USE m_mpi_bc_all,  ONLY : mpi_bc_all
#ifndef CPP_OLDINTEL
    USE m_mpi_dist_forcetheorem
#endif
#endif
#ifdef CPP_HDF
    USE m_hdf_tools
#endif
    IMPLICIT NONE
    !     Types, these variables contain a lot of data!
    TYPE(t_mpi)    ,INTENT(INOUT):: fmpi
    TYPE(t_input)    ,INTENT(OUT):: input
    TYPE(t_field),    INTENT(OUT) :: field

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
    TYPE(t_mpinp)    ,INTENT(OUT):: mpinp
    TYPE(t_hybinp)   ,INTENT(OUT):: hybinp
    TYPE(t_oneD)     ,INTENT(OUT):: oneD
    TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput
    TYPE(t_wann)     ,INTENT(OUT):: wann
    CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT)::forcetheo
    TYPE(t_gfinp)    ,INTENT(OUT):: gfinp
    TYPE(t_hub1inp)  ,INTENT(OUT):: hub1inp
    TYPE(t_nococonv), INTENT(OUT):: nococonv
    TYPE(t_enparaXML)::enparaXML
    TYPE(t_forcetheo_data)::forcetheo_data


    INTEGER, ALLOCATABLE          :: xmlElectronStates(:,:)
    INTEGER, ALLOCATABLE          :: atomTypeSpecies(:)
    INTEGER, ALLOCATABLE          :: speciesRepAtomType(:)
    REAL, ALLOCATABLE             :: xmlCoreOccs(:,:,:)
    LOGICAL, ALLOCATABLE          :: xmlPrintCoreStates(:,:)
    !     .. Local Scalars ..
    INTEGER    :: i,n,l,m1,m2,isym,iisym,numSpecies,pc,iAtom,iType,minneigd
    COMPLEX    :: cdum
    CHARACTER(len=4)              :: namex
    CHARACTER(len=12)             :: relcor, tempNumberString
    CHARACTER(LEN=20)             :: filename
    REAL                          :: a1(3),a2(3),a3(3)
    REAL                          :: dtild, phi_add
    LOGICAL                       :: l_found, l_kpts, l_exist, l_krla

#ifdef CPP_MPI
    INTEGER ierr(3)
    CALL MPI_COMM_RANK (fmpi%mpi_comm,fmpi%irank,ierr(1))
    CALL MPI_COMM_SIZE (fmpi%mpi_comm,fmpi%isize,ierr(1))
#else
    fmpi%irank=0 ; fmpi%isize=1; fmpi%mpi_comm=1
#endif
    CALL check_command_line()
#ifdef CPP_HDF
    CALL hdf_init()
#endif
    IF (fmpi%irank.EQ.0) THEN
       CALL startFleur_XMLOutput()
       IF (judft_was_argument("-info")) THEN
          CLOSE(oUnit)
          OPEN (oUnit,status='SCRATCH')
       ELSE
          IF (.NOT.judft_was_argument("-no_out")) &
               OPEN (oUnit,file='out',form='formatted',status='unknown')
       ENDIF
       CALL writeOutHeader()
       !this should be removed, it deletes output of old inf file
       OPEN (16,status='SCRATCH')
    ENDIF

    ALLOCATE(t_xcpot_inbuild::xcpot)
    !Only PE==0 reads the input and does basic postprocessing
    IF (fmpi%irank.EQ.0) THEN
       CALL fleurinput_read_xml(cell,sym,atoms,input,noco,vacuum,field,&
            sliceplot,banddos,mpinp,hybinp,oneD,coreSpecInput,&
            wann,xcpot,forcetheo_data,kpts,enparaXML,gfinp,hub1inp)
       CALL fleurinput_postprocess(Cell,Sym,Atoms,Input,Noco,Vacuum,&
            Banddos,Oned,Xcpot,Kpts,gfinp)
    END IF
    !Distribute input to all PE
    CALL fleurinput_mpi_bc(Cell,Sym,Atoms,Input,Noco,Vacuum,Field,&
         Sliceplot,Banddos,mpinp,hybinp,Oned,Corespecinput,Wann,&
         Xcpot,Forcetheo_data,Kpts,Enparaxml,gfinp,hub1inp,fmpi%Mpi_comm)
    !Remaining init is done using all PE
    CALL nococonv%init(noco)
    CALL nococonv%init_ss(noco,atoms)
    CALL ylmnorm_init(MAX(atoms%lmaxd, 2*hybinp%lexp))
    CALL gaunt_init(atoms%lmaxd+1)
    CALL enpara%init_enpara(atoms,input%jspins,input%film,enparaXML)
    CALL make_sphhar(fmpi%irank==0,atoms,sphhar,sym,cell,oneD)
    ! Store structure data (has to be performed before calling make_stars)
    CALL storeStructureIfNew(input,stars, atoms, cell, vacuum, oneD, sym, fmpi,sphhar,noco)
    CALL make_stars(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot,oneD,noco,fmpi)
    CALL make_forcetheo(forcetheo_data,cell,sym,atoms,forcetheo)
    CALL lapw_dim(kpts,cell,input,noco,nococonv,oneD,forcetheo,atoms)
    CALL input%init(noco,hybinp%l_hybrid,lapw_dim_nbasfcn)
    CALL oned%init(atoms) !call again, because make_stars modified it :-)
    CALL hybinp%init(atoms, cell, input, oneD, sym, xcpot)
    CALL kpts%init(cell, sym, input%film, hybinp%l_hybrid .or. input%l_rdmft)
    IF(fmpi%irank==0)CALL gfinp%init(atoms, sym, noco, cell, input)
    CALL gfinp%mpi_bc(fmpi%mpi_comm) !THis has to be rebroadcasted because there could be new gf elements after init_gfinp
    CALL prp_xcfft(fmpi,stars,input,cell,xcpot)
    CALL convn(fmpi%irank==0,atoms,stars)
    IF (fmpi%irank==0) CALL e_field(atoms,stars,sym,vacuum,cell,input,field%efield)
    IF (fmpi%isize>1) CALL field%mpi_bc(fmpi%mpi_comm,0)

    !At some point this should be enabled for noco as well
    IF (.NOT.noco%l_noco) &
         CALL transform_by_moving_atoms(fmpi,stars,atoms,vacuum, cell, sym, sphhar,input,oned,noco)

#ifndef _OPENACC    
    IF (fmpi%irank.EQ.0) THEN
       CALL w_inpXML(&
            atoms,vacuum,input,stars,sliceplot,forcetheo,banddos,&
            cell,sym,xcpot,noco,oneD,mpinp,hybinp,kpts,enpara,gfinp,&
            .TRUE.,[.TRUE.,.TRUE.,.TRUE.,.TRUE.])
    END IF
#endif
    !
    !--> determine more dimensions
    !



    IF (fmpi%irank.EQ.0) THEN
       CALL writeOutParameters(fmpi,input,sym,stars,atoms,vacuum,kpts,&
            oneD,hybinp,cell,banddos,sliceplot,xcpot,&
            noco,enpara,sphhar)
       CALL fleur_info(kpts)
       CALL deleteDensities()
    END IF

    !Finalize the fmpi setup
    CALL setupMPI(kpts%nkpt,input%neig,fmpi)

    !Collect some usage info
    CALL add_usage_data("A-Types",atoms%ntype)
    CALL add_usage_data("Atoms",atoms%nat)
    CALL add_usage_data("Real",sym%invs.AND..NOT.noco%l_noco&
         .AND..NOT.(noco%l_soc.AND.atoms%n_u+atoms%n_hia.GT.0))
    CALL add_usage_data("Spins",input%jspins)
    CALL add_usage_data("Noco",noco%l_noco)
    CALL add_usage_data("SOC",noco%l_soc)
    CALL add_usage_data("SpinSpiral",noco%l_ss)
    CALL add_usage_data("PlaneWaves",lapw_dim_nvd)
    CALL add_usage_data("LOs",atoms%nlotot)
    CALL add_usage_data("nkpt", kpts%nkpt)

#ifdef CPP_GPU
    CALL add_usage_data("gpu_per_node",1)
#else
    CALL add_usage_data("gpu_per_node",0)
#endif

    CALL results%init(input,atoms,kpts,noco)

    IF (fmpi%irank.EQ.0) THEN
       IF(input%gw.NE.0) CALL mixing_history_reset(fmpi)
       CALL setStartingDensity(noco%l_noco)
    END IF

    !new check mode will only run the init-part of FLEUR
    IF (judft_was_argument("-check")) CALL judft_end("Check-mode done",fmpi%irank)
#ifdef CPP_MPI
    CALL MPI_BARRIER(fmpi%mpi_comm,ierr(1))
#endif
  CONTAINS
    SUBROUTINE init_wannier()
      ! Initializations for Wannier functions (start)
      IF (fmpi%irank.EQ.0) THEN
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

      IF(fmpi%irank.EQ.0) THEN
         IF(wann%l_gwf) THEN
            OPEN(113,file=wann%param_file,status='old')
            READ(113,*)!header
            WRITE(oUnit,*) 'parameter points for HDWFs generation:'
            IF(wann%l_sgwf.OR.wann%l_ms) THEN
               WRITE(oUnit,*)'      q1       ','      q2       ','      q3'
            ELSE IF(wann%l_socgwf) THEN
               WRITE(oUnit,*)'      --       ','     phi       ','    theta'
            END IF

            DO pc = 1, wann%nparampts
               READ(113,'(3(f14.10,1x))') wann%param_vec(1,pc), wann%param_vec(2,pc), wann%param_vec(3,pc)
               wann%param_vec(:,pc) = wann%param_vec(:,pc) / wann%scale_param
               WRITE(oUnit,'(3(f14.10,1x))') wann%param_vec(1,pc), wann%param_vec(2,pc), wann%param_vec(3,pc)
               IF(wann%l_sgwf.OR.wann%l_ms) THEN
                  iAtom = 1
                  DO iType = 1, atoms%ntype
                     phi_add = tpi_const*(wann%param_vec(1,pc)*atoms%taual(1,iAtom) +&
                          wann%param_vec(2,pc)*atoms%taual(2,iAtom) +&
                          wann%param_vec(3,pc)*atoms%taual(3,iAtom))
                     wann%param_alpha(iType,pc) = nococonv%alph(iType) + phi_add
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
      END IF!(fmpi%irank.EQ.0)

#ifdef CPP_MPI
      CALL MPI_BCAST(wann%param_vec,3*wann%nparampts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr(1))
      CALL MPI_BCAST(wann%param_alpha,atoms%ntype*wann%nparampts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr(1))
      CALL MPI_BCAST(wann%l_dim,3,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr(1))
#endif

      ! Initializations for Wannier functions (end)
    END SUBROUTINE init_wannier
  END SUBROUTINE fleur_init
END MODULE m_fleur_init
