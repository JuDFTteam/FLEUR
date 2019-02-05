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

  SUBROUTINE mix( field, xcpot, DIMENSION, obsolete, sliceplot, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, hybrid, archiveType, inDen, outDen, results )

    use m_juDFT
    use m_constants
    use m_cdn_io
    use m_broyd_io
    use m_brysh1
    use m_stmix
    use m_broyden
    use m_broyden2
    use m_brysh2
    use m_metric
    use m_qfix
    use m_types
    use m_xmlOutput
    use m_umix
    USE m_kerker
    use m_types_mixvector
#ifdef CPP_MPI
    use m_mpi_bc_potden
#endif
    implicit none

    type(t_oneD),      intent(in)    :: oneD
    type(t_hybrid),    intent(in)    :: hybrid 
    type(t_input),     intent(in)    :: input
    type(t_vacuum),    intent(in)    :: vacuum
    type(t_noco),      intent(in)    :: noco
    type(t_sym),       intent(in)    :: sym
    type(t_stars),     intent(in)    :: stars
    type(t_cell),      intent(in)    :: cell
    type(t_sphhar),    intent(in)    :: sphhar
    type(t_field),     intent(inout) :: field
    class(t_xcpot),    intent(in)    :: xcpot
    type(t_dimension), intent(in)    :: dimension
    type(t_obsolete),  intent(in)    :: obsolete
    type(t_sliceplot), intent(in)    :: sliceplot
    type(t_mpi),       intent(in)    :: mpi
    type(t_atoms),     intent(inout) :: atoms !n_u is modified temporarily
    type(t_potden),    intent(inout) :: outDen
    type(t_results),   intent(inout) :: results
    type(t_potden),    intent(inout) :: inDen
    integer,           intent(in)    :: archiveType

    real                             :: fix
    type(t_potden)                   :: resDen, vYukawa
    TYPE(t_mixvector)                :: sm, fsm, fmMet, smMet,fsm_mag
    LOGICAL                          :: l_densitymatrix

    
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
    
   CALL sm%alloc()
   CALL fsm%alloc()
   CALL fmMet%alloc()
   CALL fsm_mag%alloc()
   !put input charge density into array sm
   !(in the spin polarized case the arrays sm and fsm consist of spin up and spin down densities)
   
   CALL sm%from_density(inDen)
   CALL fsm%from_density(outDen)
   
   !store the difference fsm - sm in fsm
   fsm = fsm - sm
   ! calculate Magnetisation-difference
   CALL fsm_mag%from_density(outden,swapspin=.true.)
   fsm_mag=fsm_mag-sm

   ! Apply metric w to fsm and store in fmMet:  w |fsm>
   fmMet=fsm%apply_metric()
   
   !CALL distance(fsm,fmMet)
   
    ! KERKER PRECONDITIONER
    IF( input%preconditioning_param /= 0 ) THEN
       call kerker(field, DIMENSION, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, inDen, outDen, fsm )
    ENDIF

    
    !mixing of the densities
    SELECT CASE(input%imix)
    CASE(0)
       CALL stmix(atoms,input,noco,fsm,fsm_mag,sm)
    CASE(9)
       CALL pulay()
    CASE default
       PRINT *,"New Broyden"
       
       !CALL broyden(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
       !     hybrid,mmap,nmaph,mapmt,mapvac2,nmap,fsm,sm)
    END SELECT

    !initiatlize mixed density and extract it 
    CALL sm%to_density(inDen)
      
    !fix charge of the new density
    CALL qfix(mpi,stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.FALSE., fix)

    IF(atoms%n_u.NE.n_u_keep) THEN
       inDen%mmpMat = outden%mmpMat
    END IF

    atoms%n_u=n_u_keep

    IF(vacuum%nvac.EQ.1) THEN
       inDen%vacz(:,2,:) = inDen%vacz(:,1,:)
       IF (sym%invs) THEN
          inDen%vacxy(:,:,2,:) = CONJG(inDen%vacxy(:,:,1,:))
       ELSE
          inDen%vacxy(:,:,2,:) = inDen%vacxy(:,:,1,:)
       END IF
    END IF
    
    IF (atoms%n_u > 0) THEN
       IF (.NOT.l_densityMatrixPresent) THEN
          inDen%mmpMat(:,:,:,:) = outDen%mmpMat(:,:,:,:)
          CALL resetBroydenHistory()
       END IF
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
    
7900 FORMAT (/,'---->    distance of charge densities for spin ',i2,'                 it=',i5,':',f13.6,' me/bohr**3')
7901 FORMAT (/,'----> HF distance of charge densities for spin ',i2,'                 it=',i5,':',f13.6,' me/bohr**3')
8000 FORMAT (/,'---->    distance of charge densities for it=',i5,':', f13.6,' me/bohr**3')
8001 FORMAT (/,'----> HF distance of charge densities for it=',i5,':', f13.6,' me/bohr**3')
8010 FORMAT (/,'---->    distance of spin densities for it=',i5,':', f13.6,' me/bohr**3')
8011 FORMAT (/,'----> HF distance of spin densities for it=',i5,':', f13.6,' me/bohr**3')
8020 FORMAT (4d25.14)
8030 FORMAT (10i10)

  END SUBROUTINE mix
  
END MODULE m_mix
