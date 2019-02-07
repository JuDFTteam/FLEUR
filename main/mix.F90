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

  SUBROUTINE mix( field, DIMENSION,  mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, archiveType, inDen, outDen, results )

    use m_juDFT
    use m_constants
    use m_cdn_io
    use m_stmix
    use m_broyden
    use m_qfix
    use m_types
    use m_umix
    USE m_kerker
    use m_types_mixvector
    USE m_distance
    use m_mixing_history
    implicit none

    type(t_oneD),      intent(in)    :: oneD
    type(t_input),     intent(in)    :: input
    type(t_vacuum),    intent(in)    :: vacuum
    type(t_noco),      intent(in)    :: noco
    type(t_sym),       intent(in)    :: sym
    type(t_stars),     intent(in)    :: stars
    type(t_cell),      intent(in)    :: cell
    type(t_sphhar),    intent(in)    :: sphhar
    type(t_field),     intent(inout) :: field
    type(t_dimension), intent(in)    :: dimension
    type(t_mpi),       intent(in)    :: mpi
    type(t_atoms),     intent(inout) :: atoms !n_u is modified temporarily
    type(t_potden),    intent(inout) :: outDen
    type(t_results),   intent(inout) :: results
    type(t_potden),    intent(inout) :: inDen
    integer,           intent(in)    :: archiveType

    real                             :: fix
    type(t_potden)                   :: resDen, vYukawa
    TYPE(t_mixvector),ALLOCATABLE    :: sm(:), fsm(:)
    TYPE(t_mixvector)                :: fsm_mag
    LOGICAL                          :: l_densitymatrix
    INTEGER                          :: it,maxiter

    
    MPI0_a: IF( mpi%irank == 0 ) THEN
      !determine type of mixing:
      !imix=0:straight, imix=o broyden first, imix=5:broyden second
      !imix=:generalozed anderson mixing
      select case( input%imix )
        case( 0 )
           WRITE( 6, fmt='(a,2f10.5)' ) 'STRAIGHT MIXING',input%alpha
           IF (input%jspins.EQ.1) WRITE (6,FMT='(a,2f10.5)')&
                &    'charge density mixing parameter:',input%alpha
           IF (input%jspins.EQ.2) WRITE (6,FMT='(a,2f10.5)')&
                &    'spin density mixing parameter:',input%alpha*input%spinf

        case( 3 )
          write( 6, fmt='(a,f10.5)' ) 'BROYDEN FIRST MIXING',input%alpha
        case( 5 )
          write( 6, fmt='(a,f10.5)' ) 'BROYDEN SECOND MIXING',input%alpha
        case( 7 )
          write( 6, fmt='(a,f10.5)' ) 'ANDERSON GENERALIZED',input%alpha
        case default
          call juDFT_error( "mix: input%imix =/= 0,3,5,7 ", calledby ="mix" )
      end select

      if ( input%jspins == 2 .and. input%imix /= 0 ) then
        write( 6, '(''WARNING : for QUASI-NEWTON METHODS SPINF=1'')' )
      end if
   ENDIF MPI0_a

   l_densitymatrix=.FALSE.
   IF (atoms%n_u>0) THEN
      l_densitymatrix=.NOT.input%ldaulinmix
      CALL u_mix(input,atoms,inDen%mmpMat,outDen%mmpMat)
      IF (ALL(inDen%mmpMat==0.0)) THEN
         l_densitymatrix=.FALSE.
         inDen%mmpMat=outDen%mmpMat
      ENDIF
   ENDIF
   CALL mixvector_init(mpi%mpi_comm,l_densitymatrix,oneD,input,vacuum,noco,sym,stars,cell,sphhar,atoms)
   maxiter=merge(1,input%maxiter,input%imix==0)
   CALL mixing_history(input%imix,maxiter,inden,outden,sm,fsm,it)
  
   CALL distance(mpi%irank,cell%vol,input%jspins,fsm(it),inDen,outDen,results,fsm_Mag)
   
    ! KERKER PRECONDITIONER
    IF( input%preconditioning_param /= 0 )  call kerker(field, DIMENSION, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, inDen, outDen, fsm(it) )
    
    
    !mixing of the densities
    if(input%imix==0.or.it==1) CALL stmix(atoms,input,noco,fsm(it),fsm_mag,sm(it))
    !if(it>1.and.input%imix==9) CALL pulay(input%alpha,fsm,sm)
    if(it>1.and.(input%imix==3.or.input%imix==5.or.input%imix==7)) Call broyden(input%alpha,fsm,sm)

    !initiatlize mixed density and extract it 
    CALL sm(it)%to_density(inDen)
      
    !fix charge of the new density
    CALL qfix(mpi,stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.FALSE., fix)

   

    IF(vacuum%nvac.EQ.1) THEN
       inDen%vacz(:,2,:) = inDen%vacz(:,1,:)
       IF (sym%invs) THEN
          inDen%vacxy(:,:,2,:) = CONJG(inDen%vacxy(:,:,1,:))
       ELSE
          inDen%vacxy(:,:,2,:) = inDen%vacxy(:,:,1,:)
       END IF
    END IF
    
    IF (atoms%n_u>0.AND..NOT.l_densitymatrix.AND..NOT.input%ldaulinmix) THEN
       !No density matrix was present 
       !but is now created...
       CALL mixing_history_reset()
       CALL mixvector_reset()
    ENDIF
    !write out mixed density
    CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
         1,results%last_distance,results%ef,.TRUE.,inDen)
    
#ifdef CPP_HDF
    IF (judft_was_argument("-last_extra")) THEN
       CALL system("rm cdn_last.hdf")
       CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
            1,results%last_distance,results%ef,.TRUE.,inDen,'cdn_last')
       
    END IF
#endif

    inDen%iter = inDen%iter + 1
    


  END SUBROUTINE mix
  
END MODULE m_mix
