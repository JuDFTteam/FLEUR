!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_mix

  !------------------------------------------------------------------------
  !  mixing of charge densities or potentials:
  !    IMIX = 0 : linear mixing                                     
  !    IMIX = 3 : Broyden's First method                            
  !    IMIX = 5 : Broyden's Second method                           
  !    IMIX = 7 : Generalized Anderson method                       
  !------------------------------------------------------------------------

contains

  subroutine mix( field, xcpot, dimension, obsolete, sliceplot, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, hybrid, archiveType, inDen, outDen, results )

#include"cpp_double.h"

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
    use m_kerker
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

    real                             :: fix, intfac, vacfac
    integer                          :: i, imap, js, n, lh
    integer                          :: mmap, mmaph, nmaph, nmap, mapmt, mapvac, mapvac2
    integer                          :: iofl, n_u_keep
    logical                          :: l_exist, l_ldaU, l_densityMatrixPresent, l_pot
    real, allocatable                :: sm(:), fsm(:), fmMet(:), smMet(:)
    character(len=20)                :: attributes(2)
    complex                          :: n_mmpTemp(-3:3,-3:3,max(1,atoms%n_u),input%jspins)
    type(t_potden)                   :: resDen, vYukawa
    integer                          :: ierr(2)



    MPI0_a: if( mpi%irank == 0 ) then


      !determine type of mixing:
      !imix=0:straight, imix=o broyden first, imix=5:broyden second
      !imix=:generalozed anderson mixing
      select case( input%imix )
        case( 0 )
          write( 6, fmt='(a,2f10.5)' ) 'STRAIGHT MIXING',input%alpha
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

      call sm%init(oneD,input,vacuum,noco,sym,stars,cell,sphhar,atoms)
      call fsm%alloc()
      call fmMet%alloc()
      !put input charge density into array sm
      !(in the spin polarized case the arrays sm and fsm consist of spin up and spin down densities)
      call sm%from_density(inDen)
      call fsm%from_density(outDen)
      !store the difference fsm - sm in fsm
      fsm = fsm - sm
      ! Apply metric w to fsm and store in fmMet:  w |fsm>
      fmMet=fsm%apply_metric()

      call distance()

    end if MPI0_a

    ! KERKER PRECONDITIONER
    IF( input%preconditioning_param /= 0 ) THEN
       call kerker(field, DIMENSION, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, inDen, outDen, fsm )
    ENDIF
    MPI0_c: if( mpi%irank == 0 ) then
       
    !mixing of the densities
      IF (input%imix.EQ.0) THEN
         CALL stmix(atoms,input,noco, nmap,nmaph,fsm, sm)
      ELSE
         CALL broyden(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
                      hybrid,mmap,nmaph,mapmt,mapvac2,nmap,fsm,sm)

!        Replace the broyden call above by the commented metric and broyden2 calls
!        below to switch on the continuous restart of the Broyden method.
         ! Apply metric w to sm and store in smMet:  w |sm>
!         CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
!                     mmap,nmaph,mapmt,mapvac2,sm,smMet,l_pot)
!  
!         CALL broyden2(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
!                       hybrid,mmap,nmaph,mapmt,mapvac2,nmap,fsm,sm,fmMet,smMet)
      END IF

      !initiatlize mixed density and extract it with brysh2 call
      inDen%mmpMat = CMPLX(0.0,0.0)
      call sm%to_density(inDen)
      
      !fix charge of the new density
      CALL qfix(mpi,stars,atoms,sym,vacuum, sphhar,input,cell,oneD,inDen,noco%l_noco,.FALSE.,.false., fix)

      IF(atoms%n_u.NE.n_u_keep) THEN
         inDen%mmpMat = n_mmpTemp
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

    end if MPI0_c

  end subroutine mix

end module m_mix
