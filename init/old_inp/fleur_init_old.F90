!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_init_old
  IMPLICIT NONE
CONTAINS
  !> Collection of code for old-style inp-file treatment
  SUBROUTINE fleur_init_old(mpi,&
       input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
       sliceplot,banddos,obsolete,enpara,xcpot,kpts,hybrid,&
       oneD,coreSpecInput,l_opti)
    USE m_types
    USE m_judft
    USE m_dimens
    USE m_inped
    USE m_setup
    USE m_constants
    USE m_winpXML
#ifdef CPP_MPI
#ifndef CPP_OLDINTEL
    USE m_mpi_dist_forcetheorem
#endif
#endif

    IMPLICIT NONE
    !     Types, these variables contain a lot of data!
    TYPE(t_mpi)    ,INTENT(INOUT)  :: mpi
    TYPE(t_input)    ,INTENT(INOUT):: input
    TYPE(t_dimension),INTENT(OUT)  :: DIMENSION
    TYPE(t_atoms)    ,INTENT(OUT)  :: atoms
    TYPE(t_sphhar)   ,INTENT(OUT)  :: sphhar
    TYPE(t_cell)     ,INTENT(OUT)  :: cell
    TYPE(t_stars)    ,INTENT(OUT)  :: stars
    TYPE(t_sym)      ,INTENT(OUT)  :: sym
    TYPE(t_noco)     ,INTENT(OUT)  :: noco
    TYPE(t_vacuum)   ,INTENT(OUT)  :: vacuum
    TYPE(t_sliceplot),INTENT(INOUT):: sliceplot
    TYPE(t_banddos)  ,INTENT(OUT)  :: banddos
    TYPE(t_obsolete) ,INTENT(OUT)  :: obsolete 
    TYPE(t_enpara)   ,INTENT(OUT)  :: enpara
    CLASS(t_xcpot),INTENT(OUT),ALLOCATABLE  :: xcpot
    TYPE(t_kpts)     ,INTENT(INOUT):: kpts
    TYPE(t_hybrid)   ,INTENT(OUT)  :: hybrid
    TYPE(t_oneD)     ,INTENT(OUT)  :: oneD
    TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput
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
    INTEGER:: ierr
    INCLUDE 'mpif.h'
#endif

    ALLOCATE(t_forcetheo::forcetheo) !default no forcetheorem type
    ALLOCATE(t_xcpot_inbuild::xcpot)

    SELECT TYPE(xcpot)
    TYPE IS (t_xcpot_inbuild)
    namex = '    '
    relcor = '            '

    CALL dimens(mpi,input,sym,stars,atoms,sphhar,DIMENSION,vacuum,&
         obsolete,kpts,oneD,hybrid)
    stars%kimax2= (2*stars%mx1+1)* (2*stars%mx2+1)-1
    stars%kimax = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)-1
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
    ALLOCATE (atoms%corestateoccs(1,2,atoms%ntype));atoms%corestateoccs=0.0
    ALLOCATE ( sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd),stars%ustep(stars%ng3) )
    ALLOCATE ( stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3),stars%ig2(stars%ng3) )
    ALLOCATE ( atoms%jri(atoms%ntype),stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3),sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd) )
    ALLOCATE (sym%mrot(3,3,sym%nop),sym%tau(3,sym%nop))
    ALLOCATE ( atoms%lmax(atoms%ntype),sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))!,sym%mrot(3,3,sym%nop) )
    ALLOCATE ( atoms%ncv(atoms%ntype),atoms%neq(atoms%ntype),atoms%ngopr(atoms%nat) )
    ALLOCATE ( sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd) )
    ALLOCATE ( stars%nstr2(stars%ng2),atoms%ntypsy(atoms%nat),stars%nstr(stars%ng3) )
    ALLOCATE ( stars%igfft(0:stars%kimax,2),stars%igfft2(0:stars%kimax2,2))
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
    ALLOCATE ( stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3)  )
    ALLOCATE ( kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt) )
    ALLOCATE ( stars%pgfft(0:stars%kimax),stars%pgfft2(0:stars%kimax2) )
    ALLOCATE ( stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1) )
    ALLOCATE ( atoms%bmu(atoms%ntype) )
    ALLOCATE ( atoms%l_geo(atoms%ntype) )
    ALLOCATE ( atoms%nlo(atoms%ntype),atoms%llo(atoms%nlod,atoms%ntype) )
    ALLOCATE ( atoms%lo1l(0:atoms%llod,atoms%ntype),atoms%nlol(0:atoms%llod,atoms%ntype),atoms%lapw_l(atoms%ntype) )
    ALLOCATE ( noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype),noco%l_relax(atoms%ntype) )
    ALLOCATE ( noco%b_con(2,atoms%ntype),atoms%lda_u(atoms%ntype),atoms%l_dulo(atoms%nlod,atoms%ntype) )
    ALLOCATE ( atoms%gfelem(atoms%ntype),atoms%j0(atoms%ntype))
    ALLOCATE ( sym%d_wgn(-3:3,-3:3,3,sym%nop) )
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

    atoms%nlo(:) = 0
    atoms%llo(:,:) = -1
    input%eig66(1)=.FALSE.
    ! HF/hybrid functionals/EXX
    ALLOCATE ( hybrid%nindx(0:atoms%lmaxd,atoms%ntype) )

    kpts%specificationType = 0
    atoms%numStatesProvided(:) = 0
    input%l_coreSpec = .FALSE.




    IF(.NOT.juDFT_was_argument("-toXML")) THEN
       PRINT *,"--------------WARNING----------------------"
       PRINT *,"You are using the old-style FLEUR inp file."
       PRINT *,"Please be warned that not all features are"
       PRINT *,"implemented/tested in this mode. Please "
       PRINT *,"consider switching to xml input. You might"
       PRINT *,"find the -toXML command line option useful."
       PRINT *,"--------------WARNING----------------------"
    ELSE
       IF (mpi%isize>1) CALL judft_error("Do not call -toXML with more than a single PE")
    ENDIF
    CALL timestart("preparation:stars,lattice harmonics,+etc")

    !+t3e
    IF (mpi%irank.EQ.0) THEN
       !-t3e
       CALL inped(atoms,obsolete,vacuum,input,banddos,xcpot,sym,&
            cell,sliceplot,noco,&
            stars,oneD,hybrid,kpts,a1,a2,a3,namex,relcor)
       !
       IF (xcpot%needs_grad()) THEN
          ALLOCATE (stars%ft2_gfx(0:stars%kimax2),stars%ft2_gfy(0:stars%kimax2))
          ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
               oneD%pgft1xy(0:oneD%odd%nn2d-1),&
               oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
       ELSE
          ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
          ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
               oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
       ENDIF
       oneD%odd%nq2 = stars%ng2!oneD%odd%n2d
       oneD%odi%nq2 = oneD%odd%nq2
       !-odim
       !+t3e
       INQUIRE(file="cdn1",exist=l_opti)
       IF (noco%l_noco) INQUIRE(file="rhomat_inp",exist=l_opti)
       l_opti=.NOT.l_opti
       IF ((sliceplot%iplot.NE.0).OR.(input%strho).OR.(input%swsp).OR.&
            &    (input%lflip).OR.(input%l_bmt)) l_opti = .TRUE.
       !

       namex=xcpot%get_name()
       l_krla = xcpot%data%krla.EQ.1
    END IF ! mpi%irank.eq.0

#ifdef CPP_MPI
    CALL MPI_BCAST(namex,4,MPI_CHARACTER,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(l_krla,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(input%eig66(1),1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(atoms%ntype,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
#ifndef CPP_OLDINTEL
    CALL mpi_dist_forcetheorem(mpi,forcetheo)
#endif
#endif
    IF (mpi%irank.NE.0) THEN
       CALL xcpot%init(namex,l_krla,atoms%ntype)
    END IF

    CALL setup(mpi,atoms,kpts,DIMENSION,sphhar,&
         obsolete,sym,stars,oneD,input,noco,&
         vacuum,cell,xcpot,&
         sliceplot,enpara,l_opti)

    IF (mpi%irank.EQ.0) THEN
       banddos%l_orb = .FALSE.
       banddos%orbCompAtom = 0

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
          ALLOCATE(hybrid%lcutm1(atoms%ntype),hybrid%lcutwf(atoms%ntype),hybrid%select1(4,atoms%ntype))
          filename = 'inpConverted.xml'
          xmlElectronStates = noState_const
          xmlPrintCoreStates = .FALSE.
          DO i = 1, atoms%nat
             WRITE(atoms%label(i),'(i0)') i
          END DO
          DO iType = 1, atoms%ntype
             noel(iType) = namat_const(atoms%nz(iType))
             atomTypeSpecies(iType) = iType
             speciesRepAtomType(iType) = iType

             hybrid%lcutm1(iType) = 4
             hybrid%lcutwf(iType) = atoms%lmax(iType) - atoms%lmax(iType) / 10
             hybrid%select1(:,iType) = (/4, 0, 4, 2 /)
          END DO
          hybrid%gcutm1 = input%rkmax - 0.5
          hybrid%tolerance1 = 1.0e-4
          hybrid%ewaldlambda = 3
          hybrid%lexp = 16
          hybrid%bands1 = max( nint(input%zelec)*10, 60 )

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
               cell,sym,xcpot,noco,oneD,hybrid,kpts,kpts%nkpt3,kpts%l_gamma,&
               noel,namex,relcor,a1,a2,a3,cell%amat(3,3),input%comment,&
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
    END SELECT
  END SUBROUTINE fleur_init_old
END MODULE m_fleur_init_old
