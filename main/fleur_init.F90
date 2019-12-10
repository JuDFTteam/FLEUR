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
             sliceplot,banddos,obsolete,enpara,xcpot,results,kpts,mpbasis,hybrid,&
             oneD,coreSpecInput,wann,hub1,l_opti)
          USE m_types
          USE m_judft
          USE m_juDFT_init
          USE m_init_wannier_defaults
          USE m_rinpXML
          USE m_postprocessInput
          USE m_gen_map
          USE m_dwigner
          USE m_gen_bz
          USE m_calc_tetra
          USE m_calc_tria
          USE m_ylm
          USE m_InitParallelProcesses
          USE m_checkInputParams
          USE m_xmlOutput
          USE m_constants
          USE m_winpXML
          USE m_writeOutParameters
          USE m_setupMPI
          USE m_cdn_io
          USE m_fleur_info
          USE m_mixing_history
          USE m_checks
          USE m_prpqfftmap
          USE m_writeOutHeader
          USE m_fleur_init_old
          USE m_types_xcpot_inbuild
          USE m_mpi_bc_xcpot
          USE m_wann_read_inp
          use m_gaunt, only: gaunt_init

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
          TYPE(t_obsolete) ,INTENT(OUT):: obsolete
          TYPE(t_enpara)   ,INTENT(OUT):: enpara
          CLASS(t_xcpot),ALLOCATABLE,INTENT(OUT):: xcpot
          TYPE(t_results)  ,INTENT(OUT):: results
          TYPE(t_kpts)     ,INTENT(OUT):: kpts
          TYPE(t_mpbasis), intent(inout):: mpbasis
          TYPE(t_hybrid)   ,INTENT(OUT):: hybrid
          TYPE(t_oneD)     ,INTENT(OUT):: oneD
          TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput
          TYPE(t_wann)     ,INTENT(OUT):: wann
          CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT)::forcetheo
          TYPE(t_hub1ham)  ,INTENT(OUT):: hub1 
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
          CHARACTER(len=12)             :: relcor
          CHARACTER(LEN=20)             :: filename
          REAL                          :: a1(3),a2(3),a3(3)
          REAL                          :: dtild, phi_add
          LOGICAL                       :: l_found, l_kpts, l_exist
          LOGICAL                       :: l_wann_inp

#ifdef CPP_MPI
          INCLUDE 'mpif.h'
          INTEGER ierr(3)
          CALL MPI_COMM_RANK (mpi%mpi_comm,mpi%irank,ierr)
          CALL MPI_COMM_SIZE (mpi%mpi_comm,mpi%isize,ierr)
#else
          mpi%irank=0 ; mpi%isize=1; mpi%mpi_comm=1
#endif
          !determine if we use an xml-input file
          INQUIRE (file='inp.xml',exist=input%l_inpXML)
          INQUIRE(file='inp',exist=l_found)
          IF (input%l_inpXML) THEN
             !xml found, we will use it, check if we also have a inp-file
             IF (l_found) CALL judft_warn("Both inp & inp.xml given.", calledby="fleur_init",hint="Please delete one of the input files")
          ELSE
             IF (.NOT.l_found) CALL judft_error("No input file found",calledby='fleur_init',hint="To use FLEUR, you have to provide either an 'inp' or an 'inp.xml' file in the working directory")
          END IF

          CALL check_command_line()
#ifdef CPP_HDF
          CALL hdf_init()
#endif
          !call juDFT_check_para()
          CALL field%init(input)
          input%eig66(1)=.FALSE.
          input%gw                = -1
          input%gw_neigd          =  0
          !-t3e
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
             !OPEN (16,status='SCRATCH')
          ENDIF

          input%l_rdmft = .FALSE.

          input%l_wann = .FALSE.
          CALL initWannierDefaults(wann)

          input%minDistance = 0.0
          input%ldauLinMix = .FALSE.
          input%ldauMixParam = 0.05
          input%ldauSpinf = 1.0
          input%ldauAdjEnpara = .FALSE.
          input%pallst = .FALSE.
          input%scaleCell = 1.0
          input%scaleA1 = 1.0
          input%scaleA2 = 1.0
          input%scaleC = 1.0
          input%forcealpha = 1.0
          input%forcemix = 2 ! BFGS is default.
          input%epsdisp = 0.00001
          input%epsforce = 0.00001
          input%numBandsKPoints = -1

          input%l_gfsphavg = .TRUE.
          input%l_gfmperp = .FALSE.
          input%l_gf = .FALSE.
          input%l_dftspinpol = .FALSE.
          input%l_resolvent = .FALSE.
          input%l_hist = .FALSE.
          input%gf_ne = 1300
          input%gf_ellow = 0.0
          input%gf_elup = 0.0
          input%gf_mode = 0
          input%gf_n1 = 5
          input%gf_n2 = 128
          input%gf_n3 = 5
          input%gf_sigma = 0.0314
          input%gf_nmatsub = 5
          input%gf_n = 128
          input%gf_alpha = 1.0
          input%gf_et = 0.0
          input%gf_eb = 0.0
          input%gfTet = .FALSE.
          input%minoccDistance = 0.01
          input%minmatDistance = 0.001
          atoms%n_hia = 0
          atoms%n_gf = 0
          atoms%n_j0 = 0

          kpts%ntet = 1
          kpts%numSpecialPoints = 1

          sliceplot%iplot=0
          sliceplot%kk = 0
          sliceplot%e1s = 0.0
          sliceplot%e2s = 0.0
          sliceplot%nnne = 0

          banddos%l_mcd = .FALSE.
          banddos%e_mcd_lo = -10.0
          banddos%e_mcd_up = 0.0

          banddos%unfoldband = .FALSE.
          banddos%s_cell_x = 1
          banddos%s_cell_y = 1
          banddos%s_cell_z = 1

          noco%l_mtNocoPot = .FALSE.

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
                CALL timestart("r_inpXML")
                CALL r_inpXML(&
                     atoms,obsolete,vacuum,input,stars,sliceplot,banddos,DIMENSION,forcetheo,field,&
                     cell,sym,xcpot,noco,oneD,mpbasis,hybrid,kpts,enpara,coreSpecInput,wann,&
                     noel,namex,relcor,a1,a2,a3,dtild,xmlElectronStates,&
                     xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType,&
                     l_kpts,hub1)
                CALL timestop("r_inpXML") 
             END IF
             CALL mpi_bc_xcpot(xcpot,mpi)
#ifdef CPP_MPI
#ifndef CPP_OLDINTEL
             CALL mpi_dist_forcetheorem(mpi,forcetheo)
#endif
#endif

             CALL timestart("postprocessInput")
             CALL postprocessInput(mpi,input,field,sym,stars,atoms,vacuum,obsolete,kpts,&
                                   oneD,mpbasis,hybrid,cell,banddos,sliceplot,xcpot,forcetheo,&
                                   noco,hub1,dimension,enpara,sphhar,l_opti,l_kpts)
             CALL timestop("postprocessInput")

             IF (mpi%irank.EQ.0) THEN
                filename = ''
                numSpecies = SIZE(speciesRepAtomType)
                CALL w_inpXML(&
                              atoms,obsolete,vacuum,input,stars,sliceplot,forcetheo,banddos,&
                              cell,sym,xcpot,noco,oneD,mpbasis,hybrid,kpts,kpts%nkpt3,kpts%l_gamma,&
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
                  DIMENSION,cell,sym,xcpot,noco,oneD,mpbasis,hybrid,&
                  kpts,enpara,sphhar,mpi,obsolete,hub1)
#endif

          ELSE ! else branch of "IF (input%l_inpXML) THEN"
             CALL fleur_init_old(mpi,&
                  input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
                  sliceplot,banddos,obsolete,enpara,xcpot,kpts,mpbasis,hybrid,&
                  oneD,coreSpecInput,l_opti)
          END IF ! end of else branch of "IF (input%l_inpXML) THEN"
          !

          IF (.NOT.mpi%irank==0) CALL enpara%init(atoms,input%jspins,.FALSE.)
                   !-odim
          oneD%odd%nq2 = oneD%odd%n2d
          oneD%odd%kimax2 = oneD%odd%nq2 - 1
          oneD%odd%nat = atoms%nat

          oneD%odi%d1 = oneD%odd%d1 ; oneD%odi%mb = oneD%odd%mb ; oneD%odi%M = oneD%odd%M
          oneD%odi%k3 = oneD%odd%k3 ; oneD%odi%chi = oneD%odd%chi ; oneD%odi%rot = oneD%odd%rot
          oneD%odi%invs = oneD%odd%invs ; oneD%odi%zrfs = oneD%odd%zrfs
          oneD%odi%n2d = oneD%odd%n2d ; oneD%odi%nq2 = oneD%odd%nq2 ; oneD%odi%nn2d = oneD%odd%nn2d
          oneD%odi%kimax2 = oneD%odd%kimax2 ; oneD%odi%m_cyl = oneD%odd%m_cyl
          oneD%odi%ig => oneD%ig1 ; oneD%odi%kv => oneD%kv1 ; oneD%odi%nst2 => oneD%nstr1

          oneD%ods%nop = oneD%odd%nop ; oneD%ods%nat = oneD%odd%nat
          oneD%ods%mrot => oneD%mrot1 ; oneD%ods%tau => oneD%tau1 ; oneD%ods%ngopr => oneD%ngopr1
          oneD%ods%invtab => oneD%invtab1 ; oneD%ods%multab => oneD%multab1

          oneD%odl%nn2d = oneD%odd%nn2d
          oneD%odl%igf => oneD%igfft1 ; oneD%odl%pgf => oneD%pgfft1

          oneD%odg%nn2d = oneD%odd%nn2d
          oneD%odg%pgfx => oneD%pgft1x ; oneD%odg%pgfy => oneD%pgft1y
          oneD%odg%pgfxx => oneD%pgft1xx ; oneD%odg%pgfyy => oneD%pgft1yy ; oneD%odg%pgfxy => oneD%pgft1xy
          !+odim
          !

#ifdef CPP_MPI
          CALL MPI_BCAST(l_opti,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(noco%l_noco,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(noco%l_soc,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(input%strho ,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(input%jspins,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(atoms%n_u,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(atoms%n_hia,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(atoms%n_gf,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(atoms%n_j0,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(atoms%lmaxd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          call MPI_BCAST( input%preconditioning_param, 1, MPI_DOUBLE_PRECISION, 0, mpi%mpi_comm, ierr )
#endif
          CALL ylmnorm_init(max(atoms%lmaxd, 2*hybrid%lexp))
          CALL gaunt_init(max(atoms%lmaxd+1, 2*hybrid%lexp))
          !
          !--> determine more dimensions
          !
          atoms%nlotot = 0
          IF(mpi%irank.EQ.0) THEN
             DO n = 1, atoms%ntype
                DO l = 1,atoms%nlo(n)
                   atoms%nlotot = atoms%nlotot + atoms%neq(n) * ( 2*atoms%llo(l,n) + 1 )
                ENDDO
             ENDDO
             DIMENSION%nbasfcn = DIMENSION%nvd + atoms%nlotot
          END IF
          DIMENSION%lmd     = atoms%lmaxd* (atoms%lmaxd+2)
          DIMENSION%lmplmd  = (DIMENSION%lmd* (DIMENSION%lmd+3))/2

          ALLOCATE (stars%igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))
          ALLOCATE (stars%igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1))
#ifdef CPP_MPI
          CALL mpi_bc_all(mpi,stars,sphhar,atoms,obsolete,sym,kpts,DIMENSION,input,field,&
                          banddos,sliceplot,vacuum,cell,enpara,noco,oneD,mpbasis,hybrid,hub1)
#endif

          ! Set up pointer for backtransformation from g-vector in positive
          ! domain of carge density fftibox into stars
          CALL prp_qfft_map(stars,sym,input,stars%igq2_fft,stars%igq_fft)

          !-t3e
          !-odim
          oneD%odd%nq2 = oneD%odd%n2d
          oneD%odd%kimax2 = oneD%odd%nq2 - 1
          oneD%odd%nat = atoms%nat

          oneD%odi%d1 = oneD%odd%d1 ; oneD%odi%mb = oneD%odd%mb ; oneD%odi%M = oneD%odd%M
          oneD%odi%k3 = oneD%odd%k3 ; oneD%odi%chi = oneD%odd%chi ; oneD%odi%rot = oneD%odd%rot
          oneD%odi%invs = oneD%odd%invs ; oneD%odi%zrfs = oneD%odd%zrfs
          oneD%odi%n2d = oneD%odd%n2d ; oneD%odi%nq2 = oneD%odd%nq2 ; oneD%odi%nn2d = oneD%odd%nn2d
          oneD%odi%kimax2 = oneD%odd%kimax2 ; oneD%odi%m_cyl = oneD%odd%m_cyl
          oneD%odi%ig => oneD%ig1 ; oneD%odi%kv => oneD%kv1 ; oneD%odi%nst2 => oneD%nstr1

          oneD%ods%nop = oneD%odd%nop ; oneD%ods%nat = oneD%odd%nat
          oneD%ods%mrot => oneD%mrot1 ; oneD%ods%tau => oneD%tau1 ; oneD%ods%ngopr => oneD%ngopr1
          oneD%ods%invtab => oneD%invtab1 ; oneD%ods%multab => oneD%multab1

          oneD%odl%nn2d = oneD%odd%nn2d
          oneD%odl%igf => oneD%igfft1 ; oneD%odl%pgf => oneD%pgfft1

          oneD%odg%nn2d = oneD%odd%nn2d
          oneD%odg%pgfx => oneD%pgft1x ; oneD%odg%pgfy => oneD%pgft1y
          oneD%odg%pgfxx => oneD%pgft1xx ; oneD%odg%pgfyy => oneD%pgft1yy ; oneD%odg%pgfxy => oneD%pgft1xy
          !+odim
          IF (noco%l_noco) DIMENSION%nbasfcn = 2*DIMENSION%nbasfcn

          IF( sym%invs .OR. noco%l_soc ) THEN
             sym%nsym = sym%nop
          ELSE
             ! combine time reversal symmetry with the spatial symmetry opera
             ! thus the symmetry operations are doubled
             sym%nsym = 2*sym%nop
          END IF

          CALL checkInputParams(mpi,input,dimension,atoms,sym,noco,xcpot,oneD,forcetheo)

          ! Initializations for Wannier functions (start)
          IF (mpi%irank.EQ.0) THEN
#ifdef CPP_WANN
             INQUIRE(FILE='plotbscomf',EXIST=wann%l_bs_comf)
             WRITE(*,*)'l_bs_comf=',wann%l_bs_comf
             WRITE(*,*) 'Logical variables for wannier functions to be read in!!'
#endif
             wann%l_gwf = wann%l_ms.OR.wann%l_sgwf.OR.wann%l_socgwf

             if(wann%l_gwf) then
                WRITE(*,*)'running HDWF-extension of FLEUR code'
                WRITE(*,*)'with l_sgwf =',wann%l_sgwf,' and l_socgwf =',wann%l_socgwf

                IF(wann%l_socgwf.AND. .NOT.noco%l_soc) THEN
                  CALL juDFT_error("set l_soc=T if l_socgwf=T",calledby="fleur_init")
                END IF

                IF((wann%l_ms.or.wann%l_sgwf).AND..NOT.(noco%l_noco.AND.noco%l_ss)) THEN
                   CALL juDFT_error("set l_noco=l_ss=T for l_sgwf.or.l_ms",calledby="fleur_init")
                END IF

                IF((wann%l_ms.or.wann%l_sgwf).and.wann%l_socgwf) THEN
                   CALL juDFT_error("(l_ms.or.l_sgwf).and.l_socgwf",calledby="fleur_init")
                END IF

                INQUIRE(FILE=wann%param_file,EXIST=l_exist)
                IF(.NOT.l_exist) THEN
                   CALL juDFT_error("where is param_file"//trim(wann%param_file)//"?",calledby="fleur_init")
                END IF
                OPEN (113,file=wann%param_file,status='old')
                READ (113,*) wann%nparampts,wann%scale_param
                CLOSE(113)
             ELSE
                wann%nparampts=1
                wann%scale_param=1.0
             END IF
          END IF

#ifdef CPP_MPI
          CALL MPI_BCAST(wann%l_bs_comf,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(wann%l_gwf,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(wann%nparampts,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(wann%scale_param,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)

          CALL MPI_BCAST(wann%l_sgwf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(wann%l_socgwf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(wann%l_ms,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

          ALLOCATE (wann%param_vec(3,wann%nparampts))
          ALLOCATE (wann%param_alpha(atoms%ntype,wann%nparampts))

          IF(mpi%irank.EQ.0) THEN
             IF(wann%l_gwf) THEN
                OPEN(113,file=wann%param_file,status='old')
                READ(113,*)!header
                write(6,*) 'parameter points for HDWFs generation:'
                IF(wann%l_sgwf.or.wann%l_ms) THEN
                   WRITE(6,*)'      q1       ','      q2       ','      q3'
                ELSE IF(wann%l_socgwf) THEN
                   WRITE(6,*)'      --       ','     phi       ','    theta'
                END IF

                DO pc = 1, wann%nparampts
                   READ(113,'(3(f14.10,1x))') wann%param_vec(1,pc), wann%param_vec(2,pc), wann%param_vec(3,pc)
                   wann%param_vec(:,pc) = wann%param_vec(:,pc) / wann%scale_param
                   WRITE(6,'(3(f14.10,1x))') wann%param_vec(1,pc), wann%param_vec(2,pc), wann%param_vec(3,pc)
                   IF(wann%l_sgwf.or.wann%l_ms) THEN
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

                IF(ANY(wann%param_vec(1,:).NE.wann%param_vec(1,1))) wann%l_dim(1)=.true.
                IF(ANY(wann%param_vec(2,:).NE.wann%param_vec(2,1))) wann%l_dim(2)=.true.
                IF(ANY(wann%param_vec(3,:).NE.wann%param_vec(3,1))) wann%l_dim(3)=.true.

                CLOSE(113)

                IF(wann%l_dim(1).and.wann%l_socgwf) THEN
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

          IF (xcpot%is_hybrid().OR.input%l_rdmft) THEN

!             IF( ANY( atoms%l_geo  ) )&
!                  &     CALL juDFT_error("Forces not implemented for HF/PBE0/HSE ",&
!                  &                    calledby ="fleur")

             !calculate whole Brilloun zone
             !CALL gen_bz(kpts,sym)
             CALL gen_map(atoms,sym,oneD,hybrid)

             ! calculate d_wgn
             ALLOCATE (hybrid%d_wgn2(-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,sym%nsym))
             hybrid%d_wgn2 =  CMPLX(0.0,0.0)
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
             IF ( banddos%dos .AND. banddos%ndir == -3) THEN
                WRITE(*,*) 'Recalculating k point grid to cover the full BZ.'
                CALL gen_bz(kpts,sym)
                kpts%nkpt = kpts%nkptf
                IF(ALLOCATED(kpts%bk)) DEALLOCATE(kpts%bk,kpts%wtkpt)
                ALLOCATE(kpts%bk(3,kpts%nkptf),kpts%wtkpt(kpts%nkptf))
                kpts%bk(:,:) = kpts%bkf(:,:)
                IF (kpts%nkpt3(1)*kpts%nkpt3(2)*kpts%nkpt3(3).NE.kpts%nkptf) THEN
                   IF(kpts%l_gamma) THEN
                      kpts%wtkpt = 1.0 / (kpts%nkptf-1)
                      DO i = 1, kpts%nkptf
                         IF(ALL(kpts%bk(:,i).EQ.0.0)) THEN
                            kpts%wtkpt(i) = 0.0
                         END IF
                      END DO
                   ELSE
                      CALL juDFT_error("nkptf does not match product of nkpt3(i).",calledby="fleur_init")
                   END IF
                ELSE
                   kpts%wtkpt = 1.0 / kpts%nkptf
                END IF
             END IF
             IF (atoms%n_gf>0) THEN
              IF(.NOT.input%tria.AND..NOT.input%l_hist.AND..NOT.banddos%band) THEN
                IF(input%film) THEN
                  CALL calc_tria(kpts,cell,input,sym)
                ELSE
                  !Calculate regular decomposition into tetrahedra (make_tetra doesnt seem to work for most meshes)
                  IF(kpts%nkptf.EQ.0) CALL gen_bz(kpts,sym)
                  CALL calc_tetra(kpts,cell,input,sym)
                ENDIF
              ENDIF
             ENDIF
             ALLOCATE(hybrid%map(0,0),hybrid%tvec(0,0,0),hybrid%d_wgn2(0,0,0,0))
             hybrid%l_calhf = .FALSE.
          END IF

          IF(input%l_rdmft) THEN
             hybrid%l_calhf = .FALSE.
          END IF

          IF (mpi%irank.EQ.0) THEN
             CALL writeOutParameters(mpi,input,sym,stars,atoms,vacuum,obsolete,kpts,&
                                     oneD,hybrid,cell,banddos,sliceplot,xcpot,&
                                     noco,dimension,enpara,sphhar)
             CALL fleur_info(kpts)
             CALL deleteDensities()
          END IF

          !Finalize the MPI setup
          CALL setupMPI(kpts%nkpt,DIMENSION%neigd*MERGE(2,1,noco%l_soc.AND..NOT.noco%l_noco),mpi)

          !Collect some usage info
          CALL add_usage_data("A-Types",atoms%ntype)
          CALL add_usage_data("Atoms",atoms%nat)
          CALL add_usage_data("Real",sym%invs.AND..NOT.noco%l_noco.AND..NOT.(noco%l_soc.AND.atoms%n_u+atoms%n_hia>0))
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

          INQUIRE (file='wann_inp',exist=l_wann_inp)
          input%l_wann = input%l_wann.OR.l_wann_inp
 
          IF(input%l_wann) THEN
            CALL wann_read_inp(DIMENSION,input,noco,mpi,wann)
          END IF

          CALL results%init(dimension,input,atoms,kpts,noco)

          IF (mpi%irank.EQ.0) THEN
             IF(input%gw.NE.0) CALL mixing_history_reset(mpi)
             CALL setStartingDensity(noco%l_noco)
          END IF
          
          !new check mode will only run the init-part of FLEUR
          IF (judft_was_argument("-check")) CALL judft_end("Check-mode done",mpi%irank)




        END SUBROUTINE fleur_init
      END MODULE m_fleur_init
