!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_init
      IMPLICIT NONE
      CONTAINS
        SUBROUTINE fleur_init(mpi,&
             input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,&
             sliceplot,banddos,obsolete,enpara,xcpot,results,jij,kpts,hybrid,&
             oneD,l_opti)
          USE m_judft
          USE m_juDFT_init
          USE m_types
          USE m_rinpXML
          USE m_dimens
          USE m_inped
          USE m_setup
          USE m_gen_map
          USE m_dwigner
          USE m_gen_bz
          USE m_icorrkeys
          USE m_ylm
          USE m_InitParallelProcesses
          USE m_xmlOutput
          USE m_winpXML
          USE m_setupMPI
          USE m_cdn_io
          USE m_fleur_info
          USE m_checks
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
          LOGICAL,          INTENT(OUT):: l_opti


          INTEGER, ALLOCATABLE          :: xmlElectronStates(:,:)
          INTEGER, ALLOCATABLE          :: atomTypeSpecies(:)
          INTEGER, ALLOCATABLE          :: speciesRepAtomType(:)
          REAL, ALLOCATABLE             :: xmlCoreOccs(:,:,:)
          LOGICAL, ALLOCATABLE          :: xmlPrintCoreStates(:,:)
          CHARACTER(len=3), ALLOCATABLE :: noel(:)
          !     .. Local Scalars ..
          INTEGER    :: i,n,l,m1,m2,isym,iisym,numSpecies
          COMPLEX    :: cdum
          CHARACTER(len=4)              :: namex
          CHARACTER(len=12)             :: relcor
          CHARACTER(LEN=20)             :: filename
          REAL                          :: a1(3),a2(3),a3(3)
          REAL                          :: scale, dtild
          LOGICAL                       :: l_found


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


#ifdef CPP_MPI
          INCLUDE 'mpif.h'
          INTEGER ierr(3)
          CALL MPI_COMM_RANK (mpi%mpi_comm,mpi%irank,ierr)
          CALL MPI_COMM_SIZE (mpi%mpi_comm,mpi%isize,ierr)

          sliceplot%iplot=.FALSE.
#else
          mpi%irank=0 ; mpi%isize=1; mpi%mpi_comm=1
#endif
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
#if !(defined(__TOS_BGQ__)||defined(__PGI))
             !Do not open out-file on BlueGene
             OPEN (6,file='out',form='formatted',status='unknown')
#endif
             OPEN (16,file='inf',form='formatted',status='unknown')
          ENDIF

          input%minDistance = 0.0
          kpts%ntet = 1
          kpts%numSpecialPoints = 1
          IF (input%l_inpXML) THEN
             IF (mpi%irank.EQ.0) THEN
                ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
                ALLOCATE(noel(1),atomTypeSpecies(1),speciesRepAtomType(1))
                ALLOCATE(xmlElectronStates(1,1),xmlPrintCoreStates(1,1))
                ALLOCATE(xmlCoreOccs(1,1,1))
                namex = '    '
                relcor = '            '
                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                scale = 1.0
                CALL r_inpXML(&
                     atoms,obsolete,vacuum,input,stars,sliceplot,banddos,DIMENSION,&
                     cell,sym,xcpot,noco,Jij,oneD,hybrid,kpts,enpara,sphhar,l_opti,&
                     noel,namex,relcor,a1,a2,a3,scale,dtild,xmlElectronStates,&
                     xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType)

                ALLOCATE (results%force(3,atoms%ntype,DIMENSION%jspd))
                ALLOCATE (results%force_old(3,atoms%ntype))
                results%force(:,:,:) = 0.0

                filename = ''
                numSpecies = SIZE(speciesRepAtomType)
                CALL w_inpXML(&
                     &                        atoms,obsolete,vacuum,input,stars,sliceplot,banddos,&
                     &                        cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,(/1,1,1/),kpts%l_gamma,&
                     &                        noel,namex,relcor,a1,a2,a3,scale,dtild,input%comment,&
                     &                        xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
                     &                        atomTypeSpecies,speciesRepAtomType,.TRUE.,filename,&
                     &                        .TRUE.,numSpecies,enpara)
                DEALLOCATE(noel,atomTypeSpecies,speciesRepAtomType)
                DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
             END IF

#ifdef CPP_MPI
             CALL initParallelProcesses(atoms,vacuum,input,stars,sliceplot,banddos,&
                  DIMENSION,cell,sym,xcpot,noco,jij,oneD,hybrid,&
                  kpts,enpara,sphhar,mpi,results,obsolete)
#endif

          ELSE ! else branch of "IF (input%l_inpXML) THEN"

             CALL dimens(&
                  &            mpi,input,&
                  &            sym,stars,&
                  &            atoms,sphhar,&
                  &            DIMENSION,vacuum,&
                  &            obsolete,kpts,&
                  &            oneD,hybrid,Jij)


             DIMENSION%nn2d= (2*stars%mx1+1)* (2*stars%mx2+1)
             DIMENSION%nn3d= (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)
             !-odim
             IF (oneD%odd%d1) THEN
                oneD%odd%k3 = stars%mx3
                oneD%odd%nn2d = (2*(oneD%odd%k3) + 1)*(2*(oneD%odd%M) + 1)
             ELSE
                oneD%odd%k3 = 0 ; oneD%odd%M =0 ; oneD%odd%nn2d = 1
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
             ALLOCATE ( noco%soc_opt(atoms%ntype+2) )
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
             ALLOCATE ( hybrid%select1(4,atoms%ntype),hybrid%lcutm1(atoms%ntype),&
                  &           hybrid%select2(4,atoms%ntype),hybrid%lcutm2(atoms%ntype),hybrid%lcutwf(atoms%ntype) )
             ALLOCATE ( hybrid%ddist(DIMENSION%jspd) )
             hybrid%ddist     = 1.
             !

             atoms%numStatesProvided(:) = 0

             atoms%vr0(:)         = 0.0
             jij%M(:)             = 0.0
             jij%l_magn(:)        =.FALSE.
             results%force(:,:,:) = 0.0

             CALL timestart("preparation:stars,lattice harmonics,+etc")
             !--- J< 
             jij%l_wr=.TRUE.
             jij%nqptd=1
             jij%nmagn=1
             jij%mtypes=1
             jij%phnd=1
             !--- J>
             !+t3e
             IF (mpi%irank.EQ.0) THEN
                !-t3e
                CALL inped( &
                     &           atoms,obsolete,vacuum,&
                     &           input,banddos,xcpot,sym,&
                     &           cell,sliceplot,noco,&
                     &           stars,oneD,jij,hybrid,kpts)
                !
                IF (xcpot%igrd.NE.0) THEN
                   ALLOCATE (stars%ft2_gfx(0:DIMENSION%nn2d-1),stars%ft2_gfy(0:DIMENSION%nn2d-1))
                   !-odim
                   ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
                        &             oneD%pgft1xy(0:oneD%odd%nn2d-1),&
                        &             oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
                ELSE
                   ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
                   !-odim
                   ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
                        &             oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
                ENDIF
                oneD%odd%nq2 = oneD%odd%n2d
                !-odim
                !+t3e
                INQUIRE(file="cdn1",exist=l_opti)
                IF (noco%l_noco) INQUIRE(file="rhomat_inp",exist=l_opti)
                l_opti=.NOT.l_opti
                IF ((sliceplot%iplot).OR.(input%strho).OR.(input%swsp).OR.&
                     &    (input%lflip).OR.(obsolete%l_f2u).OR.(obsolete%l_u2f).OR.(input%l_bmt)) l_opti = .TRUE.
                !
                CALL setup(&
                     &     atoms,kpts,DIMENSION,sphhar,&
                     &     obsolete,sym,stars,oneD,input,noco,&
                     &     vacuum,cell,xcpot,&
                     &     sliceplot,enpara,l_opti)
                !
                stars%ng3=stars%ng3 ; stars%ng2=stars%ng2 
                !+t3e
             ENDIF ! mpi%irank.eq.0
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


#ifdef CPP_MPI
          CALL mpi_bc_all(&
               &           mpi,stars,sphhar,atoms,obsolete,&
               &           sym,kpts,jij,DIMENSION,input,&
               &           banddos,sliceplot,vacuum,cell,enpara,&
               &           noco,oneD,xcpot,hybrid)
          ! initialize record length of the eig file

#endif 
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


          IF (    (xcpot%icorr.EQ.icorr_hf ) .OR. (xcpot%icorr.EQ.icorr_pbe0)&
               &    .OR.(xcpot%icorr.EQ.icorr_exx) .OR. (xcpot%icorr.EQ.icorr_hse)&
               &    .OR.(xcpot%icorr.EQ.icorr_vhse) ) THEN
             IF (input%film .OR. oneD%odi%d1)&
                  &    CALL juDFT_error("2D film and 1D calculations not implemented"&
                  &                 //"for HF/EXX/PBE0/HSE", calledby ="fleur",&
                  &                 hint="Use a supercell or a different functional")

             IF( ANY( atoms%l_geo  ) )&
                  &     CALL juDFT_error("Forces not implemented for HF/PBE0/HSE ",&
                  &                    calledby ="fleur")

             IF (.NOT. obsolete%pot8) STOP 'Choose pot8=T'
             !calculate whole Brilloun zone
             CALL gen_bz(kpts,sym)
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
                CALL gen_bz(kpts,sym)
             END IF
             ALLOCATE(hybrid%map(0,0),hybrid%tvec(0,0,0),hybrid%d_wgn2(0,0,0,0))
             hybrid%l_calhf   = .FALSE.
          END IF

          IF (mpi%irank.EQ.0) THEN
             CALL fleur_info(kpts)
          END IF

          !Finalize the MPI setup
          CALL setupMPI(kpts%nkpt,mpi)

          IF (mpi%irank.EQ.0) THEN
             CALL setStartingDensity(noco%l_noco)
          END IF

          !new check mode will only run the init-part of FLEUR
          IF (judft_was_argument("-check")) CALL judft_end("Check-mode done",mpi%irank)


        END SUBROUTINE fleur_init
     END MODULE
