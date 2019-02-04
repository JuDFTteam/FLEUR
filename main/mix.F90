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
    real                             :: dist(6)
    real, allocatable                :: sm(:), fsm(:), fmMet(:), smMet(:)
    character(len=20)                :: attributes(2)
    complex                          :: n_mmpTemp(-3:3,-3:3,max(1,atoms%n_u),input%jspins)
    type(t_potden)                   :: resDen, vYukawa
    integer                          :: ierr(2)


    !External functions
    real :: CPP_BLAS_sdot
    external :: CPP_BLAS_sdot

    ! YM: I have exported 'vol' from outside, be aware
    !     IF (film) THEN
    !        vol = 2.0 * z1 * area
    !     ELSE
    !        vol = omtil
    !     ENDIF

    MPI0_a: if( mpi%irank == 0 ) then

      l_densityMatrixPresent = any( inDen%mmpMat(:,:,:,:) /= 0.0 )

      !In systems without inversions symmetry the interstitial star-
      !coefficients are complex. Thus twice as many numbers have to be
      !stored.
      intfac = 2.0
      if ( sym%invs ) intfac = 1.0

      !The corresponding is true for the coeff. of the warping vacuum
      !density depending on the two dimensional inversion.
      vacfac = 2.0
      if ( sym%invs2 ) vacfac = 1.0

      mmaph = intfac * stars%ng3 + atoms%ntype * ( sphhar%nlhd + 1 ) * atoms%jmtd + &
              vacfac * vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + vacuum%nmzd * vacuum%nvac
      mmap  = mmaph * input%jspins
      !in a non-collinear calculations extra space is needed for the
      !off-diag. part of the density matrix. these coeff. are generally
      !complex independ of invs and invs2.
      if ( noco%l_noco ) then
         mmap = mmap + 2 * stars%ng3 + 2 * vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + &
              2 * vacuum%nmzd * vacuum%nvac
         IF (noco%l_mtnocopot) mmap= mmap+ 2*atoms%ntype * ( sphhar%nlhd + 1 ) * atoms%jmtd 
      end if

      ! LDA+U (start)
      n_mmpTemp = inDen%mmpMat
      n_u_keep = atoms%n_u
      if ( atoms%n_u > 0 ) call u_mix( input, atoms, inDen%mmpMat, outDen%mmpMat )
      if ( l_densityMatrixPresent ) then
        !In an LDA+U caclulation, also the density matrix is included in the
        !supervectors (sm,fsm) if no linear mixing is performed on it.
        if ( input%ldauLinMix ) then
          atoms%n_u = 0
        else
          mmap = mmap + 7 * 7 * 2 * atoms%n_u * input%jspins ! add 7*7 complex numbers per atoms%n_u and spin
        end if
      else
        atoms%n_u = 0
      end if
      ! LDA+U (end)

      allocate( sm(mmap), fsm(mmap) )
      allocate( smMet(mmap), fmMet(mmap) )
      dist(:) = 0.0

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

      !put input charge density into array sm
      !(in the spin polarized case the arrays sm and fsm consist of spin up and spin down densities)
      call brysh1( input, stars, atoms, sphhar, noco, vacuum, sym, oneD, &
                   intfac, vacfac, inDen, nmap, nmaph, mapmt, mapvac, mapvac2, sm ) 

      !put output charge density into array fsm
      call brysh1( input, stars, atoms, sphhar, noco, vacuum, sym, oneD, &
                   intfac, vacfac, outDen, nmap, nmaph, mapmt, mapvac, mapvac2, fsm )

      !store the difference fsm - sm in fsm
      fsm(:nmap) = fsm(:nmap) - sm(:nmap)

      l_pot = .false.
      ! Apply metric w to fsm and store in fmMet:  w |fsm>
      call metric( cell, atoms, vacuum, sphhar, input, noco, stars, sym, oneD, &
                   mmap, nmaph, mapmt, mapvac2, fsm, fmMet, l_pot )

      !calculate the distance of charge densities for each spin
      IF(hybrid%l_calhf) THEN
         CALL openXMLElement('densityConvergence',(/'units  ','comment'/),(/'me/bohr^3','HF       '/))
      ELSE
         CALL openXMLElement('densityConvergence',(/'units'/),(/'me/bohr^3'/))
      END IF

      DO js = 1,input%jspins
         dist(js) = CPP_BLAS_sdot(nmaph,fsm(nmaph*(js-1)+1),1, fmMet(nmaph*(js-1)+1),1)

         attributes = ''
         WRITE(attributes(1),'(i0)') js
         WRITE(attributes(2),'(f20.10)') 1000*SQRT(ABS(dist(js)/cell%vol))
         CALL writeXMLElementForm('chargeDensity',(/'spin    ','distance'/),attributes,reshape((/4,8,1,20/),(/2,2/)))
         IF( hybrid%l_calhf ) THEN
            WRITE ( 6,FMT=7901) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
         ELSE
            WRITE ( 6,FMT=7900) js,inDen%iter,1000*SQRT(ABS(dist(js)/cell%vol))
         END IF
      END DO
      IF (noco%l_noco) dist(6) = CPP_BLAS_sdot((nmap-2*nmaph), fsm(nmaph*2+1),1,fmMet(nmaph*2+1),1)
      IF (noco%l_noco) WRITE (6,FMT=7900) 3,inDen%iter,1000*SQRT(ABS(dist(6)/cell%vol))

      !calculate the distance of total charge and spin density
      !|rho/m(o) - rho/m(i)| = |rh1(o) -rh1(i)|+ |rh2(o) -rh2(i)| +/_
      !                        +/_2<rh2(o) -rh2(i)|rh1(o) -rh1(i)>
      IF (input%jspins.EQ.2) THEN
         dist(3) = CPP_BLAS_sdot(nmaph,fsm,1,fmMet(nmaph+1),1)
         dist(4) = dist(1) + dist(2) + 2.0e0*dist(3)
         dist(5) = dist(1) + dist(2) - 2.0e0*dist(3)
         CALL writeXMLElementFormPoly('overallChargeDensity',(/'distance'/),&
                                      (/1000*SQRT(ABS(dist(4)/cell%vol))/),reshape((/10,20/),(/1,2/)))
         CALL writeXMLElementFormPoly('spinDensity',(/'distance'/),&
                                      (/1000*SQRT(ABS(dist(5)/cell%vol))/),reshape((/19,20/),(/1,2/)))
         IF( hybrid%l_calhf ) THEN
            WRITE ( 6,FMT=8001) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
            WRITE ( 6,FMT=8011) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
         ELSE
            WRITE ( 6,FMT=8000) inDen%iter,1000*SQRT(ABS(dist(4)/cell%vol))
            WRITE ( 6,FMT=8010) inDen%iter,1000*SQRT(ABS(dist(5)/cell%vol))
         END IF

         !dist/vol should always be >= 0 ,
         !but for dist=0 numerically you might obtain dist/vol < 0
         !(e.g. when calculating non-magnetic systems with jspins=2).
      END IF
      results%last_distance=maxval(1000*SQRT(ABS(dist/cell%vol)))
      DEALLOCATE (smMet,fmMet)
      CALL closeXMLElement('densityConvergence')

    end if MPI0_a

    ! KERKER PRECONDITIONER
    IF( input%preconditioning_param /= 0 ) THEN
       call kerker(field, DIMENSION, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, inDen, outDen, fsm ,mapmt,mapvac,mapvac2,nmap,nmaph  )
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

      CALL brysh2(input,stars,atoms,sphhar,noco,vacuum,sym,sm,oneD,inDen) 
      DEALLOCATE (sm,fsm)

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
