MODULE m_cdnval
  use m_juDFT
CONTAINS
  SUBROUTINE cdnval(eig_id, mpi,kpts,jspin,sliceplot,noco, input,banddos,cell,atoms,enpara,stars,&
                    vacuum,dimension,sphhar,sym,obsolete,vTot,oneD,coreSpecInput,cdnvalKLoop,den,regCharges,results,&
                    moments)
    !
    !     ***********************************************************
    !         this subroutin is a modified version of cdnval.F.
    !         it calculates a layer charge distribution and an orbital
    !         composition for each state in a film.
    !         this information is written in the  'ek_orbcomp'  file
    !                                    Yu.Koroteev  01.12.2003
    !     ***********************************************************
    !     flapw7 valence density generator
    !                                         c.l.fu
    !     zelec used to calculate ne - 6.12.95 r.pentcheva
    !
    !     changed subroutine to allow parallel writing of vacdos&dosinp
    !     used temporary direct access file 84,tmp_dos to store data used
    !     in cdninf
    !     call of cdninf changed, sympsi is called from cdnval now
    !     look for 'ifdef CPP_MPI' -blocks!               d.wortmann 6.5.99
    !
    !******** ABBREVIATIONS ************************************************
    !     nbands   : number of bands in the energy window
    !     noccbd   : number of occupied bands
    !     slice    : when set to .true. the charge density of a enery range
    !                (slice) or a single state is calculated
    !     e1s,e2s  : (used if slice) energy range for the slice. if both
    !                are set to 0.0 the charge density of the band nr. nnne
    !                is calculated
    !     pallst   : (used if slice) if set to .true. bands above the
    !                Fermi-Energy are taken into account
    !     nnne     : (used if slice) number of the band of which the charge
    !                density is calculated if e1s and e2s are 0.0
    !     kk       : (used if slice) if set to 0 all k-points are used to
    !                calculate the charge density of the slice, otherwise
    !                only k-points kk is taken into account
    !     nslibd   : number of bands in slice
    !     ener     : band energy averaged over all bands and k-points,
    !                wheighted with the l-like charge of each atom type
    !     sqal     : l-like charge of each atom type. sum over all k-points
    !                and bands
    !***********************************************************************
    !
    USE m_constants
    USE m_eig66_io,ONLY: write_dos
    USE m_genMTBasis
    USE m_calcDenCoeffs
    USE m_mcdinit
    USE m_sympsi
    USE m_eparas      ! energy parameters and partial charges
    USE m_qal21       ! off-diagonal part of partial charges
    USE m_abcof
    USE m_nmat        ! calculate density matrix for LDA + U
    USE m_vacden
    USE m_pwden
    USE m_forcea8
    USE m_forcea12
    USE m_forcea21
    USE m_checkdopall
    USE m_int21       ! integrate (spin) off-diagonal radial functions
    USE m_int21lo     ! -"- for u_lo
    USE m_rhomt21     ! calculate (spin) off-diagonal MT-density coeff's
    USE m_rhonmt21    ! -"-                       non-MT-density coeff's
    USE m_cdnmt       ! calculate the density and orbital moments etc.
    USE m_orbmom      ! coeffd for orbital moments
    USE m_qmtsl       ! These subroutines divide the input%film into vacuum%layers
    USE m_qintsl      ! (slabs) and intergate the DOS in these vacuum%layers
    USE m_orbcomp     ! calculate corbital composition (like p_x,p_y,p_z)
    USE m_Ekwritesl   ! and write to file.
    USE m_abcrot2
    USE m_doswrite
    USE m_eig66_io, ONLY : read_eig
    USE m_corespec, only : l_cs    ! calculation of core spectra (EELS)
    USE m_corespec_io, only : corespec_init
    USE m_corespec_eval, only : corespec_gaunt,corespec_rme,corespec_dos,corespec_ddscs
#ifdef CPP_MPI
    USE m_mpi_col_den ! collect density data from parallel nodes
#endif
    USE m_types
    USE m_xmlOutput
    IMPLICIT NONE
    TYPE(t_results),       INTENT(INOUT) :: results
    TYPE(t_mpi),           INTENT(IN)    :: mpi
    TYPE(t_dimension),     INTENT(IN)    :: dimension
    TYPE(t_oneD),          INTENT(IN)    :: oneD
    TYPE(t_enpara),        INTENT(INOUT) :: enpara
    TYPE(t_obsolete),      INTENT(IN)    :: obsolete
    TYPE(t_banddos),       INTENT(IN)    :: banddos
    TYPE(t_sliceplot),     INTENT(IN)    :: sliceplot
    TYPE(t_input),         INTENT(IN)    :: input
    TYPE(t_vacuum),        INTENT(IN)    :: vacuum
    TYPE(t_noco),          INTENT(IN)    :: noco
    TYPE(t_sym),           INTENT(IN)    :: sym
    TYPE(t_stars),         INTENT(IN)    :: stars
    TYPE(t_cell),          INTENT(IN)    :: cell
    TYPE(t_kpts),          INTENT(IN)    :: kpts
    TYPE(t_sphhar),        INTENT(IN)    :: sphhar
    TYPE(t_atoms),         INTENT(IN)    :: atoms
    TYPE(t_coreSpecInput), INTENT(IN)    :: coreSpecInput
    TYPE(t_potden),        INTENT(IN)    :: vTot
    TYPE(t_cdnvalKLoop),   INTENT(IN)    :: cdnvalKLoop
    TYPE(t_potden),        INTENT(INOUT) :: den
    TYPE(t_regionCharges), INTENT(INOUT) :: regCharges
    TYPE(t_moments),       INTENT(INOUT) :: moments

    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN)              :: eig_id,jspin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif

    !     .. Local Scalars ..
    INTEGER :: ikpt,jsp_start,jsp_end,ispin,jsp
    INTEGER :: i,ie,ivac,j,k,l,n,ilo,nbands,noccbd
    INTEGER :: skip_t,skip_tt
    INTEGER :: nStart,nEnd,nbasfcn
    LOGICAL :: l_fmpl,l_evp,l_orbcomprot,l_real, l_write

    !     ...Local Arrays ..
    INTEGER, ALLOCATABLE :: jsym(:),ksym(:)
    REAL,    ALLOCATABLE :: we(:)
    REAL,    ALLOCATABLE :: eig(:)
    REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:) ! radial functions

    TYPE (t_lapw)             :: lapw
    TYPE (t_orb)              :: orb
    TYPE (t_denCoeffs)        :: denCoeffs
    TYPE (t_denCoeffsOffdiag) :: denCoeffsOffdiag
    TYPE (t_force)            :: force
    TYPE (t_slab)             :: slab
    TYPE (t_eigVecCoeffs)     :: eigVecCoeffs
    TYPE (t_mcd)              :: mcd
    TYPE (t_usdus)            :: usdus
    TYPE (t_zMat)             :: zMat
    TYPE (t_orbcomp)          :: orbcomp
    TYPE (t_gVacMap)          :: gVacMap

    l_real = sym%invs.AND.(.NOT.noco%l_soc).AND.(.NOT.noco%l_noco)

    !---> l_fmpl is meant as a switch to to a plot of the full magnet.
    !---> density without the atomic sphere approximation for the magnet.
    !---> density. It is not completely implemented (lo's missing).
    l_fmpl = .false.
    IF (noco%l_mperp) THEN
       !--->    when the off-diag. part of the desinsity matrix, i.e. m_x and
       !--->    m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
       !--->    is called only once. therefore, several spin loops have been
       !--->    added. if l_mperp = F, these loops run only from jspin - jspin.
       jsp_start = 1
       jsp_end   = 2
    ELSE
       jsp_start = jspin
       jsp_end   = jspin
    ENDIF
    !---> if l_mperp = F, these variables are only needed for one spin
    !---> at a time, otherwise for both spins:
    ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )      ! Deallocation before mpi_col_den
    ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
    ALLOCATE ( flo(atoms%jmtd,2,atoms%nlod,dimension%jspd) )
    ALLOCATE ( jsym(dimension%neigd),ksym(dimension%neigd) )

    ! --> Initializations

    CALL usdus%init(atoms,input%jspins)
    CALL denCoeffs%init(atoms,sphhar,jsp_start,jsp_end)
    CALL denCoeffsOffdiag%init(atoms,noco,sphhar,l_fmpl)
    CALL force%init1(input,atoms)
    CALL orb%init(atoms,noco,jsp_start,jsp_end)
    CALL mcd%init1(banddos,dimension,input,atoms)
    CALL slab%init(banddos,dimension,atoms,cell)
    CALL orbcomp%init(banddos,dimension,atoms)

    IF ((l_fmpl).AND.(.not.noco%l_mperp)) CALL juDFT_error("for fmpl set noco%l_mperp = T!" ,calledby ="cdnval")
    IF ((banddos%ndir.EQ.-3).AND.banddos%dos.AND.oneD%odi%d1) CALL juDFT_error("layer-resolved feature does not work with 1D",calledby ="cdnval")

! calculation of core spectra (EELS) initializations -start-
    CALL corespec_init(input,atoms,coreSpecInput)
    IF(l_cs.AND.(mpi%isize.NE.1)) CALL juDFT_error('EELS + MPI not implemented', calledby = 'cdnval')
    IF(l_cs.AND.jspin.EQ.1) CALL corespec_gaunt()
! calculation of core spectra (EELS) initializations -end-

    IF (mpi%irank==0) THEN
       WRITE (6,FMT=8000) jspin
       WRITE (16,FMT=8000) jspin
       CALL openXMLElementPoly('mtCharges',(/'spin'/),(/jspin/))
    END IF
8000 FORMAT (/,/,10x,'valence density: spin=',i2)

    l_write = input%cdinf.AND.mpi%irank==0

    DO n = 1,atoms%ntype

       DO ispin = jsp_start, jsp_end
          CALL genMTBasis(atoms,enpara,vTot,mpi,n,ispin,l_write,usdus,f(:,:,0:,ispin),g(:,:,0:,ispin),flo(:,:,:,ispin))
       END DO

       IF (noco%l_mperp) THEN
          DO l = 0,atoms%lmax(n)
             CALL int_21(f,g,atoms,n,l,denCoeffsOffdiag)
          END DO
          DO ilo = 1, atoms%nlo(n)
             CALL int_21lo(f,g,atoms,n,flo,ilo,denCoeffsOffdiag)
          END DO
       END IF

       IF (banddos%l_mcd) CALL mcd_init(atoms,input,dimension,vTot%mt(:,0,:,:),g,f,mcd,n,jspin)

       IF(l_cs) CALL corespec_rme(atoms,input,n,dimension%nstd,input%jspins,jspin,results%ef,&
                                  dimension%msh,vTot%mt(:,0,:,:),f,g)
    END DO
    DEALLOCATE (f,g,flo)

    skip_tt = dot_product(enpara%skiplo(:atoms%ntype,jspin),atoms%neq(:atoms%ntype))
    IF (noco%l_soc.OR.noco%l_noco)  skip_tt = 2 * skip_tt
    ALLOCATE (we(MAXVAL(cdnvalKLoop%noccbd(:))))
    ALLOCATE (eig(MAXVAL(cdnvalKLoop%noccbd(:))))
    jsp = MERGE(1,jspin,noco%l_noco)

    DO ikpt = cdnvalKLoop%ikptStart, cdnvalKLoop%nkptExtended, cdnvalKLoop%ikptIncrement

       IF (ikpt.GT.kpts%nkpt) THEN
#ifdef CPP_MPI
          CALL MPI_BARRIER(mpi%mpi_comm,ie) ! Synchronizes the RMA operations
#endif
          EXIT
       END IF

       CALL lapw%init(input,noco, kpts,atoms,sym,ikpt,cell,.false., mpi)
       skip_t = skip_tt
       noccbd = cdnvalKLoop%noccbd(ikpt)
       nStart = cdnvalKLoop%nStart(ikpt)
       nEnd = cdnvalKLoop%nEnd(ikpt)

       we=0.0
       IF(noccbd.GT.0) we(1:noccbd) = results%w_iks(nStart:nEnd,ikpt,jsp)
       IF ((sliceplot%slice).AND.(input%pallst)) we(:) = kpts%wtkpt(ikpt)

       IF (cdnvalKLoop%l_evp) THEN
          IF (nStart > skip_tt) skip_t = 0
          IF (nEnd <= skip_tt) skip_t = noccbd
          IF ((nStart <= skip_tt).AND.(nEnd > skip_tt)) skip_t = mod(skip_tt,noccbd)
       END IF

       nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
       CALL zMat%init(l_real,nbasfcn,noccbd)
       CALL read_eig(eig_id,ikpt,jsp,n_start=nStart,n_end=nEnd,neig=nbands,zmat=zMat)
#ifdef CPP_MPI
       CALL MPI_BARRIER(mpi%mpi_comm,ie) ! Synchronizes the RMA operations
#endif

       eig(1:noccbd) = results%eig(nStart:nEnd,ikpt,jsp)

       CALL gVacMap%init(dimension,sym,atoms,vacuum,stars,lapw,input,cell,kpts,enpara,vTot,ikpt,jspin)

       IF (noccbd.EQ.0) GO TO 199

       !     ----> add in spin-doubling factor
       we(:noccbd) = 2.0 * we(:noccbd) / input%jspins

       !--->    valence density in the interstitial and vacuum region
       !--->    has to be called only once (if jspin=1) in the non-collinear case

       !     ----> valence density in the interstitial region
       IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
          CALL timestart("cdnval: pwden")
          CALL pwden(stars,kpts,banddos,oneD,input,mpi,noco,cell,atoms,sym,ikpt,&
                     jspin,lapw,noccbd,we,eig,den,regCharges%qis,results,force%f_b8,zMat)
          CALL timestop("cdnval: pwden")
       END IF

       !--->    charge of each valence state in this k-point of the SBZ
       !--->    in the layer interstitial region of the film
       IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
          IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
             CALL q_int_sl(jspin,stars,atoms,sym,cell,noccbd,lapw,slab,oneD,zMat)
          END IF
       END IF

       !--->    valence density in the vacuum region
       IF (input%film) THEN
          IF (.NOT.((jspin.EQ.2) .AND. noco%l_noco)) THEN
             CALL timestart("cdnval: vacden")
             CALL vacden(vacuum,dimension,stars,oneD, kpts,input, cell,atoms,noco,banddos,&
                         gVacMap,we,ikpt,jspin,vTot%vacz(:,:,jspin),noccbd,lapw,enpara%evac0,eig,&
                         den,regCharges%qvac,regCharges%qvlay,regCharges%qstars,zMat)
             CALL timestop("cdnval: vacden")
          END IF
          !--->       perform Brillouin zone integration and summation over the
          !--->       bands in order to determine the vacuum energy parameters.
          DO ispin = jsp_start,jsp_end
             DO ivac = 1,vacuum%nvac
                regCharges%pvac(ivac,ispin)=regCharges%pvac(ivac,ispin)+dot_product(eig(:noccbd)*regCharges%qvac(:noccbd,ivac,ikpt,ispin),we(:noccbd))
                regCharges%svac(ivac,ispin)=regCharges%svac(ivac,ispin)+dot_product(regCharges%qvac(:noccbd,ivac,ikpt,ispin),we(:noccbd))
             END DO
          END DO
       END IF

       !--->    valence density in the atomic spheres
       CALL eigVecCoeffs%init(dimension,atoms,noco,jspin,noccbd)

       DO ispin = jsp_start,jsp_end
          IF (input%l_f) CALL force%init2(noccbd,input,atoms)

          CALL timestart("cdnval: abcof")
          CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,noco,ispin,oneD,&
                     eigVecCoeffs%acof(:,0:,:,ispin),eigVecCoeffs%bcof(:,0:,:,ispin),&
                     eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat,eig,force)
          CALL timestop("cdnval: abcof")

          IF (atoms%n_u.GT.0) CALL n_mat(atoms,sym,noccbd,usdus,ispin,we,eigVecCoeffs,den%mmpMat(:,:,:,jspin))

          !--->       perform Brillouin zone integration and summation over the
          !--->       bands in order to determine the energy parameters for each
          !--->       atom and angular momentum
          IF (.not.sliceplot%slice) THEN
             CALL eparas(ispin,atoms,noccbd,mpi,ikpt,noccbd,we,eig,&
                         skip_t,cdnvalKLoop%l_evp,eigVecCoeffs,usdus,regCharges,mcd,banddos%l_mcd)

             IF (noco%l_mperp.AND.(ispin == jsp_end)) THEN
                CALL qal_21(dimension,atoms,input,noccbd,noco,eigVecCoeffs,denCoeffsOffdiag,regCharges)
             END IF
          END IF

          !--->    layer charge of each valence state in this k-point of the SBZ
          !--->    from the mt-sphere region of the film
          IF (banddos%dos.AND.(banddos%ndir.EQ.-3))  THEN
             CALL q_mt_sl(ispin,atoms,noccbd,ikpt,noccbd,skip_t,noccbd,eigVecCoeffs,usdus,slab)

             INQUIRE (file='orbcomprot',exist=l_orbcomprot)
             IF (l_orbcomprot) CALL abcrot2(atoms,noccbd,eigVecCoeffs,ispin) ! rotate ab-coeffs

             CALL orb_comp(ispin,noccbd,atoms,noccbd,usdus,eigVecCoeffs,orbcomp)
          END IF

          CALL calcDenCoeffs(atoms,sphhar,sym,we,noccbd,eigVecCoeffs,ispin,denCoeffs)

          IF (noco%l_soc) CALL orbmom(atoms,noccbd,we,ispin,eigVecCoeffs,orb)

          IF (input%l_f) THEN
             CALL timestart("cdnval: force_a12/21")
             IF (.not.input%l_useapw) THEN
                CALL force_a12(atoms,noccbd,sym,dimension,cell,oneD,&
                               we,ispin,noccbd,usdus,eigVecCoeffs,force,results)
             ENDIF
             CALL force_a21(input,atoms,dimension,noccbd,sym,oneD,cell,we,ispin,&
                            enpara%el0(0:,:,ispin),noccbd,eig,usdus,eigVecCoeffs,force,results)
             CALL timestop("cdnval: force_a12/21")
          END IF

          IF(l_cs) THEN
             CALL corespec_dos(atoms,usdus,ispin,dimension%lmd,kpts%nkpt,ikpt,dimension%neigd,&
                               noccbd,results%ef,banddos%sig_dos,eig,we,eigVecCoeffs)
          END IF
       END DO !--->    end loop over ispin

       IF (noco%l_mperp) THEN
          CALL rhomt21(atoms,we,noccbd,eigVecCoeffs,denCoeffsOffdiag)
          IF (l_fmpl) CALL rhonmt21(atoms,sphhar,we,noccbd,sym,eigVecCoeffs,denCoeffsOffdiag)
       END IF

199    CONTINUE
       IF ((banddos%dos .OR. banddos%vacdos .OR. input%cdinf)  ) THEN
          CALL timestart("cdnval: write_info")
          !--->    calculate charge distribution of each state (l-character ...)
          !--->    and write the information to the files dosinp and vacdos
          !--->    for dos and bandstructure plots

          !--dw    since z is no longer an argument of cdninf sympsi has to be called here!
          IF (banddos%ndir.GT.0) CALL sympsi(lapw,jspin,sym,dimension,nbands,cell,eig,noco,ksym,jsym,zMat)

          CALL write_dos(eig_id,ikpt,jspin,regCharges,slab,orbcomp,ksym,jsym,mcd%mcd)

          CALL timestop("cdnval: write_info")
       END IF
    END DO !---> end of k-point loop

#ifdef CPP_MPI
    CALL timestart("cdnval: mpi_col_den")
    DO ispin = jsp_start,jsp_end
       CALL mpi_col_den(mpi,sphhar,atoms,oneD,stars,vacuum,input,noco,l_fmpl,ispin,regCharges,&
                        results,denCoeffs,orb,denCoeffsOffdiag,den,den%mmpMat(:,:,:,jspin))
    END DO
    CALL timestop("cdnval: mpi_col_den")
#endif

    IF (mpi%irank==0) THEN
       CALL cdnmt(dimension%jspd,atoms,sphhar,noco,l_fmpl,jsp_start,jsp_end,&
                  enpara,vTot%mt(:,0,:,:),denCoeffs,usdus,orb,denCoeffsOffdiag,moments,den%mt)

       IF(l_cs) CALL corespec_ddscs(jspin,input%jspins)

       IF (((jspin.eq.input%jspins).OR.noco%l_mperp) .AND. (banddos%dos.or.banddos%vacdos.or.input%cdinf) ) THEN
          CALL timestart("cdnval: dos")
          CALL doswrite(eig_id,dimension,kpts,atoms,vacuum,input,banddos,&
                        sliceplot,noco,sym,cell,mcd,results,slab%nsld,oneD)
          IF (banddos%dos.AND.(banddos%ndir.EQ.-3)) THEN
             CALL Ek_write_sl(eig_id,dimension,kpts,atoms,vacuum,input,jspin,sym,cell,slab)
          END IF
          CALL timestop("cdnval: dos")
       END IF

       DO ispin = jsp_start,jsp_end

          !--->      check continuity of charge density
          IF (input%cdinf) THEN
             CALL timestart("cdnval: cdninf-stuff")
             WRITE (6,FMT=8210) ispin
8210         FORMAT (/,5x,'check continuity of cdn for spin=',i2)
             CALL checkDOPAll(input,dimension,sphhar,stars,atoms,sym,vacuum,oneD,cell,den,ispin)
             CALL timestop("cdnval: cdninf-stuff")
          END IF

          !--->      forces of equ. A8 of Yu et al.
          IF ((input%l_f)) THEN
             CALL timestart("cdnval: force_a8")
             CALL force_a8(input,atoms,sphhar,ispin,vTot%mt(:,:,:,ispin),den%mt,force,results)
             CALL timestop("cdnval: force_a8")
          END IF

       END DO ! end of loop ispin = jsp_start,jsp_end
       CALL closeXMLElement('mtCharges')

       IF(vacuum%nvac.EQ.1) THEN
          den%vacz(:,2,:) = den%vacz(:,1,:)
          IF (sym%invs) THEN
             den%vacxy(:,:,2,:) = CONJG(den%vacxy(:,:,1,:))
          ELSE
             den%vacxy(:,:,2,:) = den%vacxy(:,:,1,:)
          END IF
       END IF

    END IF ! end of (mpi%irank==0)
#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%mpi_comm,ie) ! Synchronizes the RMA operations
#endif

    IF ((jsp_end.EQ.input%jspins)) THEN
       IF ((banddos%dos.OR.banddos%vacdos).AND.(banddos%ndir/=-2))  CALL juDFT_end("DOS OK",mpi%irank)
       IF (vacuum%nstm.EQ.3)  CALL juDFT_end("VACWAVE OK",mpi%irank)
    END IF

  END SUBROUTINE cdnval
END MODULE m_cdnval
