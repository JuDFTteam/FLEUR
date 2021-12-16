!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mix

  !------------------------------------------------------------------------
  !  mixing of charge densities or potentials:
  !    IMIX = 0 : linear mixing
  !    IMIX = 3 : Broyden's First method
  !    IMIX = 5 : Broyden's Second method
  !    IMIX = 7 : Generalized Anderson method
  !------------------------------------------------------------------------

contains

  SUBROUTINE mix_charge( field,   fmpi, l_writehistory,&
       stars, atoms, sphhar, vacuum, input, sym, cell, noco, nococonv,&
       oneD, archiveType, xcpot, iteration, inDen, outDen, results, l_runhia, sliceplot)

    use m_juDFT
    use m_constants
    use m_cdn_io
    use m_stmix
    use m_broyden
    use m_qfix
    use m_types
    use m_umix
    use m_checkMMPmat
    USE m_kerker
    use m_pulay
    use m_a_pulay
    use m_types_mixvector
    USE m_distance
    use m_mixing_history
    use m_RelaxSpinAxisMagn
    USE m_plot
    implicit none

    type(t_oneD),      intent(in)    :: oneD
    type(t_input),     intent(in)    :: input
    type(t_vacuum),    intent(in)    :: vacuum
    type(t_noco),      intent(in)    :: noco
    type(t_nococonv),  intent(in)    :: nococonv
    TYPE(t_sym),TARGET,INTENT(in)    :: sym
    TYPE(t_stars),TARGET,INTENT(in)  :: stars
    TYPE(t_cell),TARGET,INTENT(in)   :: cell
    TYPE(t_sphhar),TARGET,INTENT(in) :: sphhar
    type(t_field),     intent(inout) :: field
    TYPE(t_sliceplot), INTENT(IN)    :: sliceplot

    type(t_mpi),       intent(in)    :: fmpi
    TYPE(t_atoms),TARGET,INTENT(in)  :: atoms
    class(t_xcpot), intent(in)       :: xcpot
    type(t_potden),    intent(inout) :: outDen
    type(t_results),   intent(inout) :: results
    type(t_potden),    intent(inout) :: inDen
    integer,           intent(in)    :: archiveType
    integer,           intent(inout) :: iteration
    LOGICAL,           INTENT(IN)    :: l_writehistory
    LOGICAL,           INTENT(IN)    :: l_runhia

    real                             :: fix
    type(t_potden)                   :: resDen, vYukawa
    TYPE(t_mixvector),ALLOCATABLE    :: sm(:), fsm(:)
    TYPE(t_mixvector)                :: fsm_mag
    LOGICAL                          :: l_densitymatrix,l_firstItU
    INTEGER                          :: it,maxiter
    INTEGER                          :: indStart_noDenmatmixing, indEnd_noDenmatmixing


    CALL timestart("Charge Density Mixing")
    l_densitymatrix=.FALSE.
    l_firstItU=.FALSE.
    !The density/potential matrices for DFT+U are split into two parts
    ! 1:atoms%n_u Are the elements for normal DFT+U
    ! atoms%n_u+1:atoms%n_u+atoms%n_hia are the elements for DFT+Hubbard 1
    ! atoms%n_u+atoms%n_hia+1:atoms%n_u+atoms%n_hia+atoms%n_opc are the elements fro DFT+OPC
    !The latter are never mixed and held constant
    indStart_noDenmatmixing = atoms%n_u + 1
    indEnd_noDenmatmixing = atoms%n_u + atoms%n_hia + atoms%n_opc

    IF (atoms%n_u>0) THEN
       l_firstItU = ALL(inDen%mmpMat(:,:,1:atoms%n_u,:)==0.0)
       l_densitymatrix=.NOT.input%ldaulinmix.AND..NOT.l_firstItU
       IF (fmpi%irank==0) CALL u_mix(input,atoms,noco,inDen%mmpMat,outDen%mmpMat)
    ENDIF

    CALL timestart("Reading of distances")
    CALL mixvector_init(fmpi%mpi_comm,l_densitymatrix,oneD,input,vacuum,noco,stars,cell,sphhar,atoms,sym)
    CALL timestart("read history")
    CALL mixing_history_open(fmpi,input%maxiter)
    CALL timestop("read history")
    maxiter=MERGE(1,input%maxiter,input%imix==0)
    CALL mixing_history(input%imix,maxiter,inden,outden,sm,fsm,it)

    CALL distance(fmpi%irank,cell%vol,input%jspins,fsm(it),inDen,outDen,results,fsm_Mag)
    CALL timestop("Reading of distances")

    ! Preconditioner for relaxation of Magnetic moments
    call precond_noco(it,vacuum,sphhar,stars,sym,oneD,cell,noco,nococonv,input,atoms,inden,outden,fsm(it))


    ! KERKER PRECONDITIONER
    IF( input%preconditioning_param /= 0 )  THEN
       CALL timestart("Preconditioner")
       CALL kerker( field,  fmpi, &
                    stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                    oneD, inDen, outDen, fsm(it) )
       !Store modified density in history
       CALL mixing_history_store(fsm(it))
       CALL timestop("Preconditioner")
    END IF

    if (atoms%n_u>0.and.fmpi%irank.ne.0.and.input%ldaulinmix) inden%mmpMat(:,:,:atoms%n_u,:)=0.0

    !mixing of the densities
    CALL timestart("Mixing")
    SELECT CASE(input%imix)
    CASE(0)
       IF (fmpi%irank==0) WRITE(oUnit, fmt='(a,f10.5,a,f10.5)' ) &
            'STRAIGHT MIXING: alpha=',input%alpha," spin-mixing=",MERGE(input%alpha*input%spinf,0.,input%jspins>1)
       CALL stmix(atoms,input,noco,fsm(it),fsm_mag,sm(it))
    CASE(3,5)
       CALL judft_error("Broyden 1/2 method not implemented")
    CASE(7)
       IF (fmpi%irank==0) WRITE(oUnit, fmt='(a,f10.5,a,i0)' ) &
            'GENERALIZED ANDERSON MIXING: alpha=',input%alpha," History-length=",it-1
       Call broyden(input%alpha,fsm,sm)
    CASE(9)
       IF (fmpi%irank==0) WRITE(oUnit, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL pulay(input%alpha,fsm,sm,0)
    CASE(11)
       IF (fmpi%irank==0) WRITE(oUnit, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'PERIODIC PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL pulay(input%alpha,fsm,sm,input%maxiter)
    CASE(13)
       IF (fmpi%irank==0) WRITE(oUnit, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'RESTARTED PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL pulay(input%alpha,fsm,sm,0)
       IF (it==input%maxiter) CALL mixing_history_limit(0) !Restarting Pulay
    CASE(15)
       IF (fmpi%irank==0) WRITE(oUnit, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'ADAPTED PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL a_pulay(input%alpha,fsm,sm)
    CASE DEFAULT
       CALL judft_error("Unknown Mixing schema")
    END SELECT
    CALL timestop("Mixing")

    CALL timestart("Postprocessing")
    !extracte mixed density
    inDen%pw=0.0;inDen%mt=0.0
    IF (ALLOCATED(inDen%vacz)) inden%vacz=0.0
    IF (ALLOCATED(inDen%vacxy)) inden%vacxy=0.0
    IF (ALLOCATED(inDen%mmpMat).AND.l_densitymatrix) inden%mmpMat(:,:,:atoms%n_u,:)=0.0
    CALL sm(it)%to_density(inDen)
    IF (atoms%n_u>0.AND.l_firstItU) THEN
       !No density matrix was present
       !but is now created...
       inden%mmpMAT(:,:,:atoms%n_u,:)=outden%mmpMat(:,:,:atoms%n_u,:)
       !Delete the history without U
       CALL mixing_history_reset(fmpi)
       CALL mixvector_reset()
    ENDIF


    IF(atoms%n_u>0) THEN
       !When the mixing of the density matrix is done together
       !with the charge density depending on the mixing scheme
       !it can become unstable
       !Check whether the mixed density matrix makes sense
       !And correct invalid elements
       CALL checkMMPmat(1,atoms%n_u,l_densitymatrix,fmpi,atoms,input,outden,inden)
    ENDIF

    IF(atoms%n_hia+atoms%n_opc>0) THEN
       !For LDA+HIA we don't use any mixing of the density matrices we just pass it on
       inDen%mmpMat(:,:,indStart_noDenmatmixing:indEnd_noDenmatmixing,:) = &
               outDen%mmpMat(:,:,indStart_noDenmatmixing:indEnd_noDenmatmixing,:)
       IF(l_runhia) THEN
          iteration = 1
          CALL mixing_history_reset(fmpi)
          CALL mixvector_reset()
       ENDIF
    ENDIF

    if(iteration == 1 .and. xcpot%vx_is_MetaGGA()) then
       CALL mixing_history_reset(fmpi)
       CALL mixvector_reset()
    endif

    if (input%periodicMixingReset > 0 .and. mod(iteration, input%periodicMixingReset) == 0) then
       call mixing_history_reset(fmpi)
       call mixvector_reset()
    endif

    call timestart("qfix")
    !fix charge of the new density
    IF (fmpi%irank==0) CALL qfix(fmpi,stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.FALSE.,.FALSE., fix)
    call timestop("qfix")



    IF(vacuum%nvac.EQ.1) THEN
       inDen%vacz(:,2,:) = inDen%vacz(:,1,:)
       IF (sym%invs) THEN
          inDen%vacxy(:,:,2,:) = CONJG(inDen%vacxy(:,:,1,:))
       ELSE
          inDen%vacxy(:,:,2,:) = inDen%vacxy(:,:,1,:)
       END IF
    END IF

    call timestart("Density output")
    !write out mixed density (but not for a plotting run)
    IF ((fmpi%irank==0).AND.(sliceplot%iplot==0)) CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
         1,results%last_distance,results%ef,results%last_mmpmatDistance,results%last_occDistance,.TRUE.,inDen)

#ifdef CPP_HDF
    IF (fmpi%irank==0.and.judft_was_argument("-last_extra")) THEN
       CALL system("rm cdn_last.hdf")
       CALL writeDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
            1,results%last_distance,results%ef,results%last_mmpmatDistance,results%last_occDistance,.TRUE.,&
            inDen,'cdn_last')
       CALL writeCoreDensity(input,atoms,inDen%mtCore,inDen%tec,inDen%qint,'cdn_last')
    END IF
#endif
    call timestop("Density output")
    inDen%iter = inDen%iter + 1

    IF (l_writehistory.AND.input%imix.NE.0) CALL mixing_history_close(fmpi)

    CALL timestop("Postprocessing")
    CALL timestop("Charge Density Mixing")
  END SUBROUTINE mix_charge

END MODULE m_mix
