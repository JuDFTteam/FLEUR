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

  SUBROUTINE mix_charge( field, DIMENSION,  mpi, l_writehistory,&
       stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
       oneD, archiveType, inDen, outDen, results ,l_runhia)

    use m_juDFT
    use m_constants
    use m_cdn_io
    use m_stmix
    use m_broyden
    use m_qfix
    use m_types
    use m_umix
    USE m_kerker
    use m_pulay
    use m_a_pulay
    use m_types_mixvector
    USE m_distance
    use m_mixing_history
    implicit none

    type(t_oneD),      intent(in)    :: oneD
    type(t_input),     intent(in)    :: input
    type(t_vacuum),    intent(in)    :: vacuum
    type(t_noco),      intent(in)    :: noco
    type(t_sym),       intent(in)    :: sym
    TYPE(t_stars),TARGET,INTENT(in)  :: stars
    TYPE(t_cell),TARGET,INTENT(in)   :: cell
    TYPE(t_sphhar),TARGET,INTENT(in) :: sphhar
    type(t_field),     intent(inout) :: field
    type(t_dimension), intent(in)    :: dimension
    type(t_mpi),       intent(in)    :: mpi
    TYPE(t_atoms),TARGET,INTENT(in)  :: atoms 
    type(t_potden),    intent(inout) :: outDen
    type(t_results),   intent(inout) :: results
    type(t_potden),    intent(inout) :: inDen
    integer,           intent(in)    :: archiveType
    LOGICAL,           INTENT(IN)    :: l_writehistory
    LOGICAL,           INTENT(IN)    :: l_runhia

    real                             :: fix
    type(t_potden)                   :: resDen, vYukawa
    TYPE(t_mixvector),ALLOCATABLE    :: sm(:), fsm(:)
    TYPE(t_mixvector)                :: fsm_mag
    LOGICAL                          :: l_densitymatrix
    INTEGER                          :: it,maxiter


    CALL timestart("Charge Density Mixing")
    l_densitymatrix=.FALSE.
    IF (atoms%n_u>0) THEN
       l_densitymatrix=.NOT.input%ldaulinmix
       IF (mpi%irank==0) CALL u_mix(input,atoms,inDen%mmpMat,outDen%mmpMat)
       IF (ALL(inDen%mmpMat==0.0)) THEN
          l_densitymatrix=.FALSE.
          inDen%mmpMat=outDen%mmpMat
          if (mpi%irank.ne.0) inden%mmpmat=0.0 
       ENDIF
    ENDIF

    IF(atoms%n_hia>0) THEN
      inDen%mmpMat(:,:,atoms%n_u+1:atoms%n_hia,:) = outDen%mmpMat(:,:,atoms%n_u+1:atoms%n_hia,:)
    ENDIF 

    CALL timestart("Reading of distances")
    CALL mixvector_init(mpi%mpi_comm,l_densitymatrix,oneD,input,vacuum,noco,sym,stars,cell,sphhar,atoms)

    CALL mixing_history_open(mpi,input%maxiter)

    maxiter=MERGE(1,input%maxiter,input%imix==0)
    CALL mixing_history(input%imix,maxiter,inden,outden,sm,fsm,it)

    CALL distance(mpi%irank,cell%vol,input%jspins,fsm(it),inDen,outDen,results,fsm_Mag)
    CALL timestop("Reading of distances")
 
    ! KERKER PRECONDITIONER
    IF( input%preconditioning_param /= 0 )  THEN 
       CALL timestart("Preconditioner")
       CALL kerker( field, DIMENSION, mpi, &
                    stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                    oneD, inDen, outDen, fsm(it) )
       !Store modified density in history
       CALL mixing_history_store(fsm(it))
       CALL timestop("Preconditioner")
    END IF
  
    CALL timestart("Mixing")
    !mixing of the densities
    SELECT CASE(input%imix)
    CASE(0)
       IF (mpi%irank==0) WRITE( 6, fmt='(a,f10.5,a,f10.5)' ) &
            'STRAIGHT MIXING: alpha=',input%alpha," spin-mixing=",MERGE(input%alpha*input%spinf,0.,input%jspins>1)
       CALL stmix(atoms,input,noco,fsm(it),fsm_mag,sm(it))
    CASE(3,5)
       CALL judft_error("Broyden 1/2 method not implemented")
    CASE(7)
       IF (mpi%irank==0) WRITE( 6, fmt='(a,f10.5,a,i0)' ) &
            'GENERALIZED ANDERSON MIXING: alpha=',input%alpha," History-length=",it-1
       Call broyden(input%alpha,fsm,sm)
    CASE(9)
       IF (mpi%irank==0) WRITE( 6, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL pulay(input%alpha,fsm,sm,0)
    CASE(11)
       IF (mpi%irank==0) WRITE( 6, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'PERIODIC PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL pulay(input%alpha,fsm,sm,input%maxiter)
    CASE(13)
       IF (mpi%irank==0) WRITE( 6, fmt='(a,f10.5,a,i0,a,i0)' ) &
            'RESTARTED PULAY MIXING: alpha=',input%alpha," History-length=",it-1,"/",input%maxiter
       CALL pulay(input%alpha,fsm,sm,0)
       IF (it==input%maxiter) CALL mixing_history_limit(0) !Restarting Pulay 
    CASE(15)
       IF (mpi%irank==0) WRITE( 6, fmt='(a,f10.5,a,i0,a,i0)' ) &
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
    IF (ALLOCATED(inDen%mmpMat).AND.l_densitymatrix) inden%mmpMat=0.0
    CALL sm(it)%to_density(inDen)
    IF (atoms%n_u>0.AND..NOT.l_densitymatrix.AND..NOT.input%ldaulinmix) THEN
       !No density matrix was present 
       !but is now created...
       inden%mmpMAT=outden%mmpMat
       CALL mixing_history_reset(mpi)
       CALL mixvector_reset()
    ENDIF

    IF (atoms%n_hia>0.AND.l_runhia) THEN
      CALL mixing_history_reset(mpi)
      CALL mixvector_reset()
   ENDIF

    !fix charge of the new density
    IF (mpi%irank==0) CALL qfix(mpi,stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.FALSE., fix)



    IF(vacuum%nvac.EQ.1) THEN
       inDen%vacz(:,2,:) = inDen%vacz(:,1,:)
       IF (sym%invs) THEN
          inDen%vacxy(:,:,2,:) = CONJG(inDen%vacxy(:,:,1,:))
       ELSE
          inDen%vacxy(:,:,2,:) = inDen%vacxy(:,:,1,:)
       END IF
    END IF


    !write out mixed density
    IF (mpi%irank==0) CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
         1,results%last_distance,results%ef,.TRUE.,inDen)

#ifdef CPP_HDF
    IF (mpi%irank==0.and.judft_was_argument("-last_extra")) THEN
       CALL system("rm cdn_last.hdf")
       CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
            1,results%last_distance,results%ef,.TRUE.,inDen,'cdn_last')

    END IF
#endif

    inDen%iter = inDen%iter + 1

    IF (l_writehistory.AND.input%imix.NE.0) CALL mixing_history_close(mpi)

    CALL timestop("Postprocessing")

    CALL timestop("Charge Density Mixing")
  END SUBROUTINE mix_charge
  
END MODULE m_mix
