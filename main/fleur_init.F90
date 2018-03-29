!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_init
      IMPLICIT NONE
      CONTAINS
        SUBROUTINE fleur_init(mpi,&
             input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
             sliceplot,banddos,obsolete,enpara,xcpot,results,jij,kpts,hybrid,&
             oneD,coreSpecInput,wann,l_opti)
          USE m_judft
          USE m_juDFT_init
          USE m_types
          USE m_init_wannier_defaults
          USE m_rinpXML
          USE m_postprocessInput
          USE m_dimens
          USE m_inped
          USE m_setup
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
          USE m_prpqfftmap
          USE m_writeOutHeader
#ifdef CPP_MPI
          USE m_mpi_bc_all,  ONLY : mpi_bc_all
#endif
#ifdef CPP_HDF
          USE m_hdf_tools
#endif
          IMPLICIT NONE
          !     Types, these variables contain a lot of data!
          TYPE(t_mpi)    ,INTENT(INOUT):: mpi
          TYPE(t_input)    ,INTENT(OUT):: input
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
          TYPE(t_xcpot)    ,INTENT(OUT):: xcpot
          TYPE(t_results)  ,INTENT(OUT):: results
          TYPE(t_jij)      ,INTENT(OUT):: jij
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
             IF (l_found.AND..NOT.(juDFT_was_argument("-xmlInput").OR.juDFT_was_argument("-xml"))) &
                  CALL judft_warn("Both inp & inp.xml given.",calledby="fleur_init",hint="Please delete one of the input files or specify -xml to use inp.xml")
          ELSE
             IF(juDFT_was_argument("-xmlInput").OR.juDFT_was_argument("-xml")) &
                  CALL judft_error("inp.xml not found",calledby="fleur_init",hint="You gave the -xml option but provided no inp.xml file")
             IF (.NOT.l_found) CALL judft_error("No input file found",calledby='fleur_init',hint="To use FLEUR, you have to provide either an 'inp' or an 'inp.xml' file in the working directory")
          END IF



          CALL check_command_line()
#ifdef CPP_HDF
          CALL hdf_init()
#endif
          results%seigscv         = 0.0
          results%te_vcoul        = 0.0
          results%te_veff         = 0.0
          results%te_exc          = 0.0
          results%te_hfex%valence = 0.0
          results%te_hfex%core    = 0.0
          results%e_ldau          = 0.0
          results%ts              = 0.0
          input%gw                = -1
          input%gw_neigd          =  0
          !-t3e
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

          input%l_wann = .FALSE.
          CALL initWannierDefaults(wann)

          input%minDistance = 0.0
          input%ldauLinMix = .FALSE.
          input%ldauMixParam = 0.05
          input%ldauSpinf = 1.0
          input%pallst = .FALSE.
          input%scaleCell = 1.0
          input%scaleA1 = 1.0
          input%scaleA2 = 1.0
          input%scaleC = 1.0

          kpts%ntet = 1
          kpts%numSpecialPoints = 1

          sliceplot%iplot=.FALSE.
          sliceplot%kk = 0
          sliceplot%e1s = 0.0
          sliceplot%e2s = 0.0
          sliceplot%nnne = 0

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
                     atoms,obsolete,vacuum,input,stars,sliceplot,banddos,DIMENSION,forcetheo,&
                     cell,sym,xcpot,noco,Jij,oneD,hybrid,kpts,enpara,coreSpecInput,wann,&
                     noel,namex,relcor,a1,a2,a3,dtild,xmlElectronStates,&
                     xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType,&
                     l_kpts)

                ALLOCATE (results%force(3,atoms%ntype,DIMENSION%jspd))
                ALLOCATE (results%force_old(3,atoms%ntype))
                results%force(:,:,:) = 0.0
             END IF

             CALL postprocessInput(mpi,input,sym,stars,atoms,vacuum,obsolete,kpts,&
                                   oneD,hybrid,jij,cell,banddos,sliceplot,xcpot,&
                                   noco,dimension,enpara,sphhar,l_opti,noel,l_kpts)

             IF (mpi%irank.EQ.0) THEN
                filename = ''
                numSpecies = SIZE(speciesRepAtomType)
                CALL w_inpXML(&
                              atoms,obsolete,vacuum,input,stars,sliceplot,forcetheo,banddos,&
                              cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,kpts%nkpt3,kpts%l_gamma,&
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
                  DIMENSION,cell,sym,xcpot,noco,jij,oneD,hybrid,&
                  kpts,enpara,sphhar,mpi,results,obsolete)
#endif

          ELSE ! else branch of "IF (input%l_inpXML) THEN"
             ALLOCATE(t_forcetheo::forcetheo) !default no forcetheorem type
             namex = '    '
             relcor = '            '

             !--- J< 
             jij%l_wr=.TRUE.
             jij%nqptd=1
             jij%nmagn=1
             jij%mtypes=1
             jij%phnd=1
             !--- J>

             CALL dimens(mpi,input,sym,stars,atoms,sphhar,DIMENSION,vacuum,&
                         obsolete,kpts,oneD,hybrid,Jij)

             DIMENSION%nn2d= (2*stars%mx1+1)* (2*stars%mx2+1)
             DIMENSION%nn3d= (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)
             !-odim
             IF (oneD%odd%d1) THEN
                oneD%odd%k3 = stars%mx3
                oneD%odd%nn2d = (2*(oneD%odd%k3) + 1)*(2*(oneD%odd%M) + 1)
             ELSE
                oneD%odd%k3 = 0 ; oneD%odd%M =0 ; oneD%odd%nn2d = 1
                oneD%odd%mb = 0
             ENDIF
             !-odim
             ALLOCATE ( atoms%nz(atoms%ntype),atoms%relax(3,atoms%ntype),atoms%nlhtyp(atoms%ntype))
             ALLOCATE ( sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd),stars%ustep(stars%ng3) )
             ALLOCATE ( stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3),stars%ig2(stars%ng3) )
             ALLOCATE ( atoms%jri(atoms%ntype),stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3),sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd) )
             ALLOCATE (sym%mrot(3,3,sym%nop),sym%tau(3,sym%nop))
             ALLOCATE ( atoms%lmax(atoms%ntype),sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))!,sym%mrot(3,3,sym%nop) )
             ALLOCATE ( atoms%ncv(atoms%ntype),atoms%neq(atoms%ntype),atoms%ngopr(atoms%nat) )
             ALLOCATE ( sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd) )
             ALLOCATE ( stars%nstr2(stars%ng2),atoms%ntypsy(atoms%nat),stars%nstr(stars%ng3) )
             ALLOCATE ( stars%igfft(0:DIMENSION%nn3d-1,2),stars%igfft2(0:DIMENSION%nn2d-1,2),atoms%nflip(atoms%ntype) )
             ALLOCATE ( atoms%ncst(atoms%ntype) )
             ALLOCATE ( vacuum%izlay(vacuum%layerd,2) )
             ALLOCATE ( sym%invarop(atoms%nat,sym%nop),sym%invarind(atoms%nat) )
             ALLOCATE ( sym%multab(sym%nop,sym%nop),sym%invtab(sym%nop) )
             ALLOCATE ( atoms%invsat(atoms%nat),sym%invsatnr(atoms%nat) )
             ALLOCATE ( atoms%lnonsph(atoms%ntype) )
             ALLOCATE ( atoms%dx(atoms%ntype),atoms%pos(3,atoms%nat))!,sym%tau(3,sym%nop) )
             ALLOCATE ( atoms%rmsh(atoms%jmtd,atoms%ntype),atoms%rmt(atoms%ntype),stars%sk2(stars%ng2),stars%sk3(stars%ng3) )
             ALLOCATE ( stars%phi2(stars%ng2) )
             ALLOCATE ( atoms%taual(3,atoms%nat),atoms%volmts(atoms%ntype),atoms%zatom(atoms%ntype) )
             ALLOCATE ( enpara%el0(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd) )
             ALLOCATE ( enpara%evac0(2,DIMENSION%jspd),stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3)  )
             ALLOCATE ( results%force(3,atoms%ntype,DIMENSION%jspd) )
             ALLOCATE ( results%force_old(3,atoms%ntype) )
             ALLOCATE ( kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt) )
             ALLOCATE ( stars%pgfft(0:DIMENSION%nn3d-1),stars%pgfft2(0:DIMENSION%nn2d-1) )
             ALLOCATE ( stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1) )
             ALLOCATE ( atoms%bmu(atoms%ntype),atoms%vr0(atoms%ntype) )
             ALLOCATE ( enpara%lchange(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd) )
             ALLOCATE ( enpara%lchg_v(2,DIMENSION%jspd),atoms%l_geo(atoms%ntype) )
             ALLOCATE ( atoms%nlo(atoms%ntype),atoms%llo(atoms%nlod,atoms%ntype),enpara%skiplo(atoms%ntype,DIMENSION%jspd) )
             ALLOCATE ( enpara%ello0(atoms%nlod,atoms%ntype,DIMENSION%jspd),enpara%llochg(atoms%nlod,atoms%ntype,DIMENSION%jspd) )
             ALLOCATE ( atoms%lo1l(0:atoms%llod,atoms%ntype),atoms%nlol(0:atoms%llod,atoms%ntype),atoms%lapw_l(atoms%ntype) )
             ALLOCATE ( noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype),noco%l_relax(atoms%ntype) )
             ALLOCATE ( jij%alph1(atoms%ntype),jij%l_magn(atoms%ntype),jij%M(atoms%ntype) )
             ALLOCATE ( jij%magtype(atoms%ntype),jij%nmagtype(atoms%ntype) )
             ALLOCATE ( noco%b_con(2,atoms%ntype),atoms%lda_u(atoms%ntype),atoms%l_dulo(atoms%nlod,atoms%ntype) )
             ALLOCATE ( enpara%enmix(DIMENSION%jspd),sym%d_wgn(-3:3,-3:3,3,sym%nop) )
             ALLOCATE ( atoms%ulo_der(atoms%nlod,atoms%ntype) )
             ALLOCATE ( atoms%numStatesProvided(atoms%ntype))
             ALLOCATE ( kpts%ntetra(4,kpts%ntet), kpts%voltet(kpts%ntet))
             !+odim
             ALLOCATE ( oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M) )
             ALLOCATE ( oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d) )
             ALLOCATE ( oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop) )
             ALLOCATE ( oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop) )
             ALLOCATE ( oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1) )
             stars%sk2(:) = 0.0 ; stars%phi2(:) = 0.0
             !-odim

             ! HF/hybrid functionals/EXX
             ALLOCATE ( hybrid%nindx(0:atoms%lmaxd,atoms%ntype) )
           
             kpts%specificationType = 0
             atoms%numStatesProvided(:) = 0
             input%l_coreSpec = .FALSE.

             jij%M(:)             = 0.0
             jij%l_magn(:)        =.FALSE.

             atoms%vr0(:)         = 0.0
             results%force(:,:,:) = 0.0

             CALL timestart("preparation:stars,lattice harmonics,+etc")

             !+t3e
             IF (mpi%irank.EQ.0) THEN
                !-t3e
                CALL inped(atoms,obsolete,vacuum,input,banddos,xcpot,sym,&
                           cell,sliceplot,noco,&
                           stars,oneD,jij,hybrid,kpts,a1,a2,a3,namex,relcor)
                !
                IF (xcpot%is_gga()) THEN
                   ALLOCATE (stars%ft2_gfx(0:DIMENSION%nn2d-1),stars%ft2_gfy(0:DIMENSION%nn2d-1))
                   ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
                             oneD%pgft1xy(0:oneD%odd%nn2d-1),&
                             oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
                ELSE
                   ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
                   ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
                             oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
                ENDIF
                oneD%odd%nq2 = oneD%odd%n2d
                oneD%odi%nq2 = oneD%odd%nq2
                !-odim
                !+t3e
                INQUIRE(file="cdn1",exist=l_opti)
                IF (noco%l_noco) INQUIRE(file="rhomat_inp",exist=l_opti)
                l_opti=.NOT.l_opti
                IF ((sliceplot%iplot).OR.(input%strho).OR.(input%swsp).OR.&
                     &    (input%lflip).OR.(input%l_bmt)) l_opti = .TRUE.
                !

                namex=xcpot%get_name()
                l_krla = xcpot%krla.EQ.1
             END IF ! mpi%irank.eq.0

#ifdef CPP_MPI
             CALL MPI_BCAST(namex,4,MPI_CHARACTER,0,mpi%mpi_comm,ierr)
             CALL MPI_BCAST(l_krla,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif
             IF (mpi%irank.NE.0) THEN
                CALL xcpot%init(namex,l_krla)
             END IF

             CALL setup(mpi,atoms,kpts,DIMENSION,sphhar,&
                        obsolete,sym,stars,oneD,input,noco,&
                        vacuum,cell,xcpot,&
                        sliceplot,enpara,l_opti)

             IF (mpi%irank.EQ.0) THEN
                !
                stars%ng3=stars%ng3 ; stars%ng2=stars%ng2 
                !+t3e
                banddos%l_orb = .FALSE.
                banddos%orbCompAtom = 0

                ALLOCATE(xcpot%lda_atom(atoms%ntype))
                ALLOCATE(noco%socscale(atoms%ntype))
                xcpot%lda_atom(:) = .FALSE.
                noco%socscale(:) = 1.0

                IF(juDFT_was_argument("-toXML")) THEN
                   WRITE(*,*) ''
                   WRITE(*,*) 'Please note:'
                   WRITE(*,*) 'The inp to xml input conversion is experimental and'
                   WRITE(*,*) 'only made for basic inp files without sophisticated'
                   WRITE(*,*) 'parametrizations. You might have to adjust the generated'
                   WRITE(*,*) 'file by hand to really obtain an adequate input file.'
                   WRITE(*,*) 'Also the generated XML input file is not meant to be'
                   WRITE(*,*) 'beautiful.'
                   WRITE(*,*) ''
                   ALLOCATE(noel(atoms%ntype),atomTypeSpecies(atoms%ntype),speciesRepAtomType(atoms%ntype))
                   ALLOCATE(xmlElectronStates(29,atoms%ntype),xmlPrintCoreStates(29,atoms%ntype))
                   ALLOCATE(xmlCoreOccs(1,1,1),atoms%label(atoms%nat))
                   filename = 'inpConverted.xml'
                   xmlElectronStates = noState_const
                   xmlPrintCoreStates = .FALSE.
                   DO i = 1, atoms%nat
                      WRITE(atoms%label(i),'(i0)') i
                   END DO
                   DO i = 1, atoms%ntype
                      noel(i) = namat_const(atoms%nz(i))
                      atomTypeSpecies(i) = i
                      speciesRepAtomType(i) = i
                   END DO
                   numSpecies = SIZE(speciesRepAtomType)
                   ALLOCATE(atoms%speciesName(numSpecies))
                   atoms%speciesName = ''
                   DO i = 1, numSpecies
                      tempNumberString = ''
                      WRITE(tempNumberString,'(i0)') i
                      atoms%speciesName(i) = TRIM(ADJUSTL(noel(speciesRepAtomType(i)))) // '-' // TRIM(ADJUSTL(tempNumberString))
                   END DO
                   a1(:) = a1(:) / input%scaleCell
                   a2(:) = a2(:) / input%scaleCell
                   a3(:) = a3(:) / input%scaleCell
                   kpts%specificationType = 3
                   sym%symSpecType = 3
                   CALL w_inpXML(&
                                 atoms,obsolete,vacuum,input,stars,sliceplot,forcetheo,banddos,&
                                 cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,kpts%nkpt3,kpts%l_gamma,&
                                 noel,namex,relcor,a1,a2,a3,dtild,input%comment,&
                                 xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
                                 atomTypeSpecies,speciesRepAtomType,.FALSE.,filename,&
                                 .TRUE.,numSpecies,enpara)
                   DEALLOCATE(atoms%speciesName, atoms%label)
                   DEALLOCATE(noel,atomTypeSpecies,speciesRepAtomType)
                   DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
                   CALL juDFT_end("Fleur inp to XML input conversion completed.")
                END IF
             END IF ! mpi%irank.eq.0
             CALL timestop("preparation:stars,lattice harmonics,+etc")

          END IF ! end of else branch of "IF (input%l_inpXML) THEN"
          !
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
          CALL MPI_BCAST(atoms%lmaxd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
#endif
          CALL ylmnorm_init(atoms%lmaxd)
          !
          !--> determine more dimensions
          !
          DIMENSION%nbasfcn = DIMENSION%nvd + atoms%nat*atoms%nlod*(2*atoms%llod+1)
          DIMENSION%lmd     = atoms%lmaxd* (atoms%lmaxd+2)
          DIMENSION%lmplmd  = (DIMENSION%lmd* (DIMENSION%lmd+3))/2


          IF (mpi%irank.EQ.0) THEN

             !--- J< 
             jij%l_jenerg = .FALSE.
             IF (jij%l_J) THEN
                OPEN (113,file='qpts',status='old')
                INQUIRE(file='jenerg',exist=jij%l_jenerg)
                OPEN (114,file='jenerg',status='unknown')
                READ (113,*) jij%nqptd
             ENDIF
          ENDIF!(mpi%irank.EQ.0)

#ifdef CPP_MPI
          CALL MPI_BCAST(jij%nqptd,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
          CALL MPI_BCAST(jij%l_jenerg,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif
          ALLOCATE ( jij%qj(3,jij%nqptd) )
          jij%nqpt=jij%nqptd
          !+t3e
          IF (mpi%irank.EQ.0) THEN
             !-t3e
             IF (jij%l_J) THEN
                IF(jij%l_disp)THEN
                   WRITE(6,*) 'q points for the magnon spectrum calculation:'
                ELSE
                   WRITE(6,*) 'q points for the J-constants calculation:'
                ENDIF
                WRITE(6,*)'      q1       ','      q2       ','      q3'
                DO i=1,jij%nqpt
                   READ (113,3333) jij%qj(1:3,i)
                   WRITE(6,3333) jij%qj(1:3,i)
                ENDDO
3333            FORMAT(3(f14.10,1x))
             ENDIF !(jij%l_J)
             !+t3e
          ENDIF
          !-t3e
          !--- J>

          !Now check for additional input files
          IF (mpi%irank.EQ.0) THEN
             IF(.NOT.banddos%l_orb) THEN
                INQUIRE(file='orbcomp',exist=banddos%l_orb)
                IF (banddos%l_orb) THEN
                   OPEN (111,file='orbcomp',form='formatted')
                   READ (111,*) banddos%orbCompAtom
                   CLOSE (111)
                END IF
             END IF
             INQUIRE(file='mcd_inp',exist=banddos%l_mcd)
          END IF


#ifdef CPP_MPI
          CALL mpi_bc_all(&
               &           mpi,stars,sphhar,atoms,obsolete,&
               &           sym,kpts,jij,DIMENSION,input,&
               &           banddos,sliceplot,vacuum,cell,enpara,&
               &           noco,oneD,xcpot,hybrid)
          ! initialize record length of the eig file

#endif

          ! Set up pointer for backtransformation from g-vector in positive 
          ! domain of carge density fftibox into stars
          ALLOCATE (stars%igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))
          ALLOCATE (stars%igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1))
          CALL prp_qfft_map(stars,sym,input,stars%igq2_fft,stars%igq_fft)

          atoms%nlotot = 0
          DO n = 1, atoms%ntype
             DO l = 1,atoms%nlo(n)
                atoms%nlotot = atoms%nlotot + atoms%neq(n) * ( 2*atoms%llo(l,n) + 1 )
             ENDDO
          ENDDO

          jij%qn = 0.0
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
          !
          !--- J<
          IF (jij%l_J) THEN
             input%itmax = 1
             jij%phnd=2
             jij%nkpt_l = CEILING(REAL(kpts%nkpt)/mpi%isize)
             ALLOCATE ( jij%eig_l(DIMENSION%neigd+5,jij%nkpt_l) )
          ELSE
             jij%nkpt_l = 1
             ALLOCATE ( jij%eig_l(DIMENSION%neigd+5,1) )
          ENDIF
          !--- J>

          IF( sym%invs .OR. noco%l_soc ) THEN
             sym%nsym = sym%nop
          ELSE
             ! combine time reversal symmetry with the spatial symmetry opera
             ! thus the symmetry operations are doubled
             sym%nsym = 2*sym%nop
          END IF

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
             IF ( banddos%dos .AND. banddos%ndir == -3 ) THEN
                WRITE(*,*) 'Recalculating k point grid to cover the full BZ.'
                CALL gen_bz(kpts,sym)
                kpts%nkpt = kpts%nkptf
                DEALLOCATE(kpts%bk,kpts%wtkpt)
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
             ALLOCATE(hybrid%map(0,0),hybrid%tvec(0,0,0),hybrid%d_wgn2(0,0,0,0))
             hybrid%l_calhf   = .FALSE.
          END IF
 
          IF (mpi%irank.EQ.0) THEN
             CALL writeOutParameters(mpi,input,sym,stars,atoms,vacuum,obsolete,kpts,&
                                     oneD,hybrid,jij,cell,banddos,sliceplot,xcpot,&
                                     noco,dimension,enpara,sphhar)
             CALL fleur_info(kpts)
             CALL deleteDensities()
          END IF

          !Finalize the MPI setup
          CALL setupMPI(kpts%nkpt,mpi)

          IF (mpi%irank.EQ.0) THEN
             CALL setStartingDensity(noco%l_noco)
          END IF

          !new check mode will only run the init-part of FLEUR
          IF (judft_was_argument("-check")) CALL judft_end("Check-mode done",mpi%irank)

          !check for broken feature
          IF ((mpi%n_size>1).and.(ANY(atoms%nlo(:)>0)).and.(noco%l_noco)) call judft_warn("Eigenvector parallelization is broken for noco&LOs")

        END SUBROUTINE fleur_init
      END MODULE m_fleur_init
