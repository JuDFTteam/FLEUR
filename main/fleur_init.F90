      MODULE m_fleur_init
      IMPLICIT NONE
      CONTAINS
        SUBROUTINE fleur_init(ivers,mpi,&
                 input,dimension,atoms,sphhar,cell,stars,sym,noco,vacuum,&
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
          USE icorrkeys
          USE m_ylm
          USE m_InitParallelProcesses
#ifdef CPP_MPI
          USE m_mpi_bc_all,  ONLY : mpi_bc_all
#endif
          IMPLICIT NONE
          !     Types, these variables contain a lot of data!
          CHARACTER(len=9),INTENT(IN)  :: ivers
          TYPE(t_mpi)    ,INTENT(INOUT):: mpi
          TYPE(t_input)    ,INTENT(OUT):: input
          TYPE(t_dimension),INTENT(OUT):: dimension
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
         

          !     .. Local Scalars ..
          INTEGER    :: i,n,l,m1,m2,isym,iisym
          COMPLEX    :: cdum
#ifdef CPP_MPI
          INCLUDE 'mpif.h'
          INTEGER ierr(3)
          CALL MPI_COMM_RANK (mpi%mpi_comm,mpi%irank,ierr)
          CALL MPI_COMM_SIZE (mpi%mpi_comm,mpi%isize,ierr)

          sliceplot%iplot=.FALSE.
#else
          mpi%irank=0 ; mpi%isize=1; mpi%mpi_comm=1
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
#ifndef  __TOS_BGQ__
             !Do not open out-file on BlueGene
             OPEN (6,file='out',form='formatted',status='unknown')
#endif
             OPEN (16,file='inf',form='formatted',status='unknown')
          ENDIF

          input%l_inpXML = .FALSE.
          kpts%numSpecialPoints = 1
          INQUIRE (file='inp.xml',exist=input%l_inpXML)
          IF(.NOT.juDFT_was_argument("-xmlInput")) THEN
             input%l_inpXML = .FALSE.
          END IF
          IF (input%l_inpXML) THEN
             IF (mpi%irank.EQ.0) THEN
                CALL r_inpXML(&
                              atoms,obsolete,vacuum,input,stars,sliceplot,banddos,dimension,&
                              cell,sym,xcpot,noco,Jij,oneD,hybrid,kpts,enpara,sphhar,l_opti)

                ALLOCATE (results%force(3,atoms%ntype,dimension%jspd))
                ALLOCATE (results%force_old(3,atoms%ntype))
                results%force(:,:,:) = 0.0

                WRITE(*,*) 'TODO: Distribute parameters and arrays to other parallel processes!'
             END IF

#ifdef CPP_MPI
             CALL initParallelProcesses(atoms,vacuum,input,stars,sliceplot,banddos,&
                                        dimension,cell,sym,xcpot,noco,jij,oneD,hybrid,&
                                        kpts,enpara,sphhar,mpi,results,obsolete)
#endif

          ELSE ! else branch of "IF (input%l_inpXML) THEN"

          CALL dimens(&
               &            mpi,ivers,input,&
               &            sym,stars,&
               &            atoms,sphhar,&
               &            dimension,vacuum,&
               &            obsolete,kpts,&
               &            oneD,hybrid,Jij)


          dimension%nn2d= (2*stars%k1d+1)* (2*stars%k2d+1)
          dimension%nn3d= (2*stars%k1d+1)* (2*stars%k2d+1)* (2*stars%k3d+1)
          !-odim
          IF (oneD%odd%d1) THEN
             oneD%odd%k3 = stars%k3d
             oneD%odd%nn2d = (2*(oneD%odd%k3) + 1)*(2*(oneD%odd%M) + 1)
          ELSE
             oneD%odd%k3 = 0 ; oneD%odd%M =0 ; oneD%odd%nn2d = 1
          ENDIF
          !-odim
          atoms%nat = atoms%natd ! This is preliminary. The value of nat changes later.
          ALLOCATE ( atoms%nz(atoms%ntypd),atoms%relax(3,atoms%ntypd),atoms%nlhtyp(atoms%ntype))
          ALLOCATE ( sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd),stars%ustep(stars%n3d) )
          ALLOCATE ( stars%ig(-stars%k1d:stars%k1d,-stars%k2d:stars%k2d,-stars%k3d:stars%k3d),stars%ig2(stars%n3d),stars%igz(stars%n3d) )
          ALLOCATE ( atoms%jri(atoms%ntypd),stars%kv2(2,stars%n2d),stars%kv3(3,stars%n3d),sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd) )
          ALLOCATE (sym%mrot(3,3,sym%nop),sym%tau(3,sym%nop))
          ALLOCATE ( atoms%lmax(atoms%ntypd),sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))!,sym%mrot(3,3,sym%nop) )
          ALLOCATE ( atoms%ncv(atoms%ntypd),atoms%neq(atoms%ntypd),atoms%ngopr(atoms%natd) )
          ALLOCATE ( sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd) )
          ALLOCATE ( stars%nstr2(stars%n2d),atoms%ntypsy(atoms%natd),stars%nstr(stars%n3d) )
          ALLOCATE ( stars%igfft(0:dimension%nn3d-1,2),stars%igfft2(0:dimension%nn2d-1,2),atoms%nflip(atoms%ntypd) )
          ALLOCATE ( atoms%ncst(atoms%ntypd) )
          ALLOCATE ( vacuum%izlay(vacuum%layerd,2) )
          ALLOCATE ( sym%invarop(atoms%natd,sym%nop),sym%invarind(atoms%natd) )
          ALLOCATE ( sym%multab(sym%nop,sym%nop),sym%invtab(sym%nop) )
          ALLOCATE ( atoms%invsat(atoms%natd),sym%invsatnr(atoms%natd) )
          ALLOCATE ( atoms%lnonsph(atoms%ntypd) )
          ALLOCATE ( atoms%dx(atoms%ntypd),atoms%pos(3,atoms%natd))!,sym%tau(3,sym%nop) )
          ALLOCATE ( atoms%rmsh(atoms%jmtd,atoms%ntypd),atoms%rmt(atoms%ntypd),stars%sk2(stars%n2d),stars%sk3(stars%n3d) )
          ALLOCATE ( stars%phi2(stars%n2d) )
          ALLOCATE ( atoms%taual(3,atoms%natd),atoms%volmts(atoms%ntypd),atoms%zatom(atoms%ntypd) )
          ALLOCATE ( enpara%el0(0:atoms%lmaxd,atoms%ntypd,dimension%jspd) )
          ALLOCATE ( enpara%evac0(2,dimension%jspd),stars%rgphs(-stars%k1d:stars%k1d,-stars%k2d:stars%k2d,-stars%k3d:stars%k3d)  )
          ALLOCATE ( results%force(3,atoms%ntypd,dimension%jspd) )
          ALLOCATE ( results%force_old(3,atoms%ntypd) )
          ALLOCATE ( kpts%bk(3,kpts%nkptd),kpts%wtkpt(kpts%nkptd) )
          ALLOCATE ( stars%pgfft(0:dimension%nn3d-1),stars%pgfft2(0:dimension%nn2d-1) )
          ALLOCATE ( stars%ufft(0:27*stars%k1d*stars%k2d*stars%k3d-1) )
          ALLOCATE ( atoms%bmu(atoms%ntypd),atoms%vr0(atoms%ntypd) )
          ALLOCATE ( enpara%lchange(0:atoms%lmaxd,atoms%ntypd,dimension%jspd) )
          ALLOCATE ( enpara%lchg_v(2,dimension%jspd),atoms%l_geo(atoms%ntypd) )
          ALLOCATE ( atoms%nlo(atoms%ntypd),atoms%llo(atoms%nlod,atoms%ntypd),enpara%skiplo(atoms%ntypd,dimension%jspd) )
          ALLOCATE ( enpara%ello0(atoms%nlod,atoms%ntypd,dimension%jspd),enpara%llochg(atoms%nlod,atoms%ntypd,dimension%jspd) )
          ALLOCATE ( atoms%lo1l(0:atoms%llod,atoms%ntypd),atoms%nlol(0:atoms%llod,atoms%ntypd),atoms%lapw_l(atoms%ntypd) )
          ALLOCATE ( noco%alph(atoms%ntypd),noco%beta(atoms%ntypd),noco%l_relax(atoms%ntypd) )
          ALLOCATE ( jij%alph1(atoms%ntypd),jij%l_magn(atoms%ntypd),jij%M(atoms%ntypd) )
          ALLOCATE ( jij%magtype(atoms%ntypd),jij%nmagtype(atoms%ntypd) )
          ALLOCATE ( noco%b_con(2,atoms%ntypd),atoms%lda_u(atoms%ntypd),atoms%l_dulo(atoms%nlod,atoms%ntypd) )
          ALLOCATE ( enpara%enmix(dimension%jspd),sym%d_wgn(-3:3,-3:3,3,sym%nop) )
          ALLOCATE ( atoms%ulo_der(atoms%nlod,atoms%ntypd) )
          ALLOCATE ( noco%soc_opt(atoms%ntypd+2) )
          ALLOCATE ( atoms%numStatesProvided(atoms%ntypd))
          !+odim
          ALLOCATE ( oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M) )
          ALLOCATE ( oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d) )
          ALLOCATE ( oneD%ngopr1(atoms%natd),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop) )
          ALLOCATE ( oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop) )
          ALLOCATE ( oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1) )
          stars%sk2(:) = 0.0 ; stars%phi2(:) = 0.0
          !-odim

          ! HF/hybrid functionals/EXX
          ALLOCATE ( hybrid%nindx(0:atoms%lmaxd,atoms%ntypd) )
          ALLOCATE ( hybrid%select1(4,atoms%ntypd),hybrid%lcutm1(atoms%ntypd),&
               &           hybrid%select2(4,atoms%ntypd),hybrid%lcutm2(atoms%ntypd),hybrid%lcutwf(atoms%ntypd) )
          ALLOCATE ( hybrid%ddist(dimension%jspd) )
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
                ALLOCATE (stars%ft2_gfx(0:dimension%nn2d-1),stars%ft2_gfy(0:dimension%nn2d-1))
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
             IF (.NOT.(input%strho.OR.obsolete%l_f2u.OR.obsolete%l_u2f.OR.sliceplot%iplot)) THEN
                IF (noco%l_noco) THEN
                   INQUIRE (file='rhomat_inp',exist=input%strho) ! if no density (rhoma
                ELSE
                   INQUIRE (file='cdn1',exist=input%strho)       ! if no density (cdn1)
                ENDIF
                input%strho = .NOT.input%strho                ! create a starting density
             ENDIF
             l_opti = .FALSE.
             IF ((sliceplot%iplot).OR.(input%strho).OR.(input%swsp).OR.&
                  &    (input%lflip).OR.(obsolete%l_f2u).OR.(obsolete%l_u2f).OR.(input%l_bmt)) l_opti = .TRUE.

             obsolete%form76 = .FALSE.
             IF (noco%l_soc.AND.obsolete%form66) THEN
                IF (.NOT.input%eonly)  CALL juDFT_error("form66 = T only with eonly = T !",calledby="fleur")
                obsolete%form66 = .FALSE.
                obsolete%form76 = .TRUE.
             ENDIF
             !
             CALL setup(&
                  &     atoms,kpts,dimension,sphhar,&
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
          oneD%odd%nat = atoms%natd

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
          dimension%nbasfcn = dimension%nvd + atoms%natd*atoms%nlod*(2*atoms%llod+1)
          dimension%lmd     = atoms%lmaxd* (atoms%lmaxd+2)
          dimension%lmplmd  = (dimension%lmd* (dimension%lmd+3))/2

          sym%l_zref=.FALSE.

          IF ((mpi%irank.EQ.0).AND.(.NOT.l_opti)) THEN
             IF (sym%zrfs.AND.(SUM(ABS(kpts%bk(3,:kpts%nkptd))).LT.1e-9)) sym%l_zref=.TRUE.
          ENDIF
          IF (mpi%irank.EQ.0) THEN
             IF (noco%l_noco) sym%l_zref = .FALSE.

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
               &           sym,kpts,jij,dimension,input,&
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
          oneD%odd%nat = atoms%natd

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
          IF (noco%l_noco) dimension%nbasfcn = 2*dimension%nbasfcn
          !
          !--- J<
          IF (jij%l_J) THEN
             input%itmax = 1
             jij%phnd=2
             jij%nkpt_l = CEILING(REAL(kpts%nkptd)/mpi%isize)
             ALLOCATE ( jij%eig_l(dimension%neigd+5,jij%nkpt_l) )
          ELSE
             jij%nkpt_l = 1
             ALLOCATE ( jij%eig_l(dimension%neigd+5,1) )
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
                IF ( obsolete%nwd /= 1 )&
                     &     STOP 'orbital decomposed DOS only implemented for 1 window!'
                CALL gen_bz(kpts,sym)
             END IF
             ALLOCATE(hybrid%map(0,0),hybrid%tvec(0,0,0),hybrid%d_wgn2(0,0,0,0))
             hybrid%l_calhf   = .FALSE.
          END IF

     END SUBROUTINE
     END MODULE
