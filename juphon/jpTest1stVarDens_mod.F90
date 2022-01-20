module jpTest1stVarDens

  implicit none

  contains

  subroutine Test1stVarDens( atoms, sym, lathar, cell, input, results, stars, kpts, dimens, qpts, usdus, ngdp2km, ngdp, logUnit,   &
      & testz1Phi0ContSw, testRho1IRsw, testRho1MTsw, testRho1BasCorrSw, testPlotRho03Dsw, testGradRho0PathSw, noPtsCon, iqpt,     &
      & paPoX, paPoY, paPoZ, harSw, extSw, xcSw, gdp, GbasVec, z, gdp2Ind, rho0IR, ne, mapGbas, nv, kveclo, rbas1, rbas2, &
      & ilo2p, nobd, mlh_atom, nmem_atom, clnu_atom, rho0MT, nRadFun, iloTable, kpq2kPrVec, mapKpq2K, gdp2iLim )

#include "cppmacro.h"
    use m_types
    use m_jpIOnMixing, only : loadDensity
    use mod_juPhonUtils, only : calcGrFinLH
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_sphhar),                 intent(in) :: lathar
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_results),                intent(in) :: results
    type(t_stars),                  intent(in) :: stars
    type(t_kpts),                   intent(in) :: kpts
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: qpts
    type(t_usdus),                  intent(in) :: usdus

    ! Scalar parameter
    integer,                        intent(in) :: ngdp2km
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit
    logical,                        intent(in) :: testz1Phi0ContSw
    logical,                        intent(in) :: testRho1IRsw
    logical,                        intent(in) :: testRho1MTsw
    logical,                        intent(in) :: testRho1BasCorrSw
    logical,                        intent(in) :: testPlotRho03Dsw
    logical,                        intent(in) :: testGradRho0PathSw
    integer,                        intent(in) :: noPtsCon !number of points to check continuity
    integer,                        intent(in) :: iqpt
    real,                           intent(in) :: paPoX
    real,                           intent(in) :: paPoY
    real,                           intent(in) :: paPoZ
    logical,                        intent(in) :: harSw
    logical,                        intent(in) :: extSw
    logical,                        intent(in) :: xcSw

    ! Array parameter
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: GbasVec(:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)
    integer,                        intent(in) :: gdp2Ind( : , : , : )
    complex,                        intent(in) :: rho0IR(:, :)
    integer,                        intent(in) :: ne(:)
    integer,                        intent(in) :: mapGbas(:, :, :) !todo rename ilst to this everywhere
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kveclo(:,:)
    real,                           intent(in) :: rbas1(:,:,:,:,:)
    real,                           intent(in) :: rbas2(:,:,:,:,:)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: nobd(:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: iloTable(:, 0:,: )
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    integer,                        intent(in) :: gdp2iLim(2, 3)

    ! Scalar variables
    integer                                    :: ikpt
    integer                                    :: idir
    integer                                    :: iDatom
    integer                                    :: iDtype
    integer                                    :: iDeqat
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: imesh
    logical                                    :: densFilePresent

    ! Array variables
    complex,           allocatable             :: z1nG(:, :, : )
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: rho1IR(:, :)
    complex,           allocatable             :: rho1MT(:, :, :, :)
    character(len=22)                          :: filename



    if (.false.) then
      call benchMarkMTGrad( atoms, input, lathar, sym, cell, clnu_atom, nmem_atom, mlh_atom  )
    end if

    if ( testz1Phi0ContSw ) then
      ! Tests continuity of z1 phi0 for all displaced atoms, k-points and direction of displacements
      allocate( z1nG(dimens%nbasfcn, maxval(nobd(:, :)), 3) )
      write(logUnit, '(a)') 'Continuity check of z1'
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          do ikpt = 1, kpts%nkpt
            call ReadInz1( atoms, ikpt, iqpt, mapKpq2K(ikpt, iqpt), iDatom, nobd, nv, z1nG )
            do idir = 1, 3
              write(logUnit, '(a,i3)')    'index of k-point         : ', ikpt
              write(logUnit, '(a,1x,i2)') 'displaced atom index     : ', iDatom
              write(logUnit, '(a,2x,i1)') 'direction of displacement: ', idir
              write(logUnit, '(a)') '_____________________________________'
              call checkz1ph0Cont( atoms, cell, sym, dimens, usdus, noPtsCon, nobd(ikpt, 1), kpts%bk(:, ikpt), ikpt, GbasVec,      &
                & mapGbas, nv, z1nG(:, :, idir), nRadFun, kveclo, rbas1, rbas2, iloTable, logUnit, sum(nobd(:64, 1) ), idir )
              write(logUnit, *)
            end do ! idir
          end do ! ikpt
        end do ! iDeqat
      end do ! iDtype
      deallocate(z1nG)
    end if

    if ( testPlotRho03Dsw ) then
      ! Creates an xsf file for the 3D representation of the unperturbed density
      call test3DplotPotDens( atoms, stars, sym, cell, dimens, lathar, input, ngdp, clnu_atom, nmem_atom, mlh_atom, gdp, rho0IR,   &
        & rho0MT, logUnit )
    end if

    if ( testGradRho0PathSw ) then
      write(*, '(2x,a)') 'Performing TestGradRho0PathSw...'
      ! todo This test is also in the sternheimer module. So if it does not work we might have a look there
      ! For q equals the Gamma point, we read in the converged rho1, calculate the final rho1 and compare it to the gradient of the
      ! unperturbed density
      allocate( rho1IR( ngdp, 3), rho1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3 ) )
      rho1IR(:, :) = cmplx(0.0, 0.0)
      rho1MT(:, :, :, :) = cmplx(0.0, 0.0)
      call calcGrFinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), grRho0MT )
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          do idir = 1, 3
            write(filename, '(a10,i1,a4,i1,a4,i2)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', 1
            inquire(file=filename, exist=densFilePresent)
            if ( densFilePresent ) then
              write(*, '(a,1x,i1,1x,a,1x,i2,1x,a,a)') 'Variation of density for displaced atom', iDatom, 'and q-point', 1,      &
                'is read from ', filename
              call loadDensity( atoms, ngdp, filename, rho1IR(:, idir), rho1MT(:, :, :, idir) )
            else
              call juDFT_warn( filename//' not found! Its first density variation is set to zero.', calledby='Test1stVarDens',     &
                & hint='Perform calculation to gather required rho1.' )
            end if
            do oqn_l = 0, atoms%lmax(iDtype)! + 1
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + mqn_m + 1
                do imesh = 1, atoms%jri(iDtype)
                  ! If we are at the displaced atom we have to add the gradient of the density
                  rho1MT(imesh, lm, iDatom, idir) =  rho1MT(imesh, lm, iDatom, idir) - grRho0MT(imesh, lm, iDatom, idir)
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! idir
          call checkGradDens0s( sym, cell, input, lathar, stars, atoms, rho1IR(:, :), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom,   &
            & rho0MT, rho1MT(:, :, :, :), paPoX, paPoY, paPoZ, harSw, extSw, xcSw, rho0IR, rho0MT )
        end do ! ieqat
      end do ! dispType
    end if

    if ( testRho1IRsw .or. testRho1MTsw .or. testRho1BasCorrSw) then
      write(*, *)
      write(*, '(a)') 'Initiating density variation test(s)...'
      write(*, '(a)') '--------------------------------'
      write(logUnit, *)
    write ( logUnit, * )
      write(logUnit, '(a)') 'Density variation test(s)'
      write(logUnit, '(a)') '*************************'
      write(logUnit, *)
    else
      write ( logUnit, * )
      write(*, '(a)') '-----------------------------------'
      write(*, '(a)') 'DISABLED density variation test(s)!'
      write(*, '(a)') '-----------------------------------'
      write(logUnit, '(a)') 'DISABLED density variation tests!'
      write(logUnit, '(a)') '*********************************'
      return
    end if


    if ( testRho1IRsw ) then
      ! Uses the juPhon routine for the linear IR density variations to reproduce the IR density from Fleur
      write(*, '(2x,a)') 'Performing TestRho1IRroutines...'
      call testRho1IRroutines( cell, input, results, stars, kpts, dimens, qpts, nobd, nv, GbasVec, z, &
        & rho0IR, ne, mapGbas, logUnit )
    else
      write(*, '(2x,a)') 'DISABLED TestRho1IRroutines!'
    write(logUnit, '(a)') 'Test reproducing interstitial unperturbed density of FLEUR with juPhon double sum method routine for the linear variation of the density'
    write(logUnit, '(a)') '----------------------------------------------------------------------------------------------------------------------------------------'
      write (logUnit, '(a)')'                                                                                                                                       |__DISABLED!'
    end if

    if ( testRho1MTsw ) then
      ! Uses the juPhon routine for the linear MT density variations to reproduce the MT density from Fleur
      write(*, '(2x,a)') 'Performing TestRho1MTroutines...'
      call testRho1MTroutine( atoms, stars, dimens, sym, cell, kpts, usdus, input, results, lathar, mlh_atom, clnu_atom, nmem_atom,&
        & nv, mapGbas, GbasVec, ne, z, kveclo, rbas1, rbas2, ilo2p, nobd, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, rho0IR, rho0MT,   &
        & logUnit )
    else
      write(*, '(2x,a)') 'DISABLED TestRho1MTroutines!'

      write(logUnit, '(a)') 'Test reproducing muffin-tin unperturbed density of FLEUR with juPhon FFT method routine for the linear variation of the density'
      write(logUnit, '(a)') '-------------------------------------------------------------------------------------------------------------------------------'
      write(logUnit, '(a)') '                                                                                                                              |__DISABLED!'
    end if

    if ( testRho1BasCorrSw ) then
      write(*, '(2x,a)') 'Performing TestRho1BasCorr...'
      ! Chooses a special z1 sothat the first density variation without core terms results in zero.
      call TestRho1BasCorr( atoms, stars, sym, cell, dimens, kpts,  input, usdus, results, nv, ne, nobd, gdp, ngdp, ngdp2km,       &
        & mapGbas, GbasVec, z, kveclo, ilo2p, gdp2Ind, gdp2iLim, rho0IR, rbas1, rbas2, kpq2kPrVec, logUnit )
    else
      write(*, '(2x,a)') 'DISABLED TestRho1BasCorr!'
      write(logUnit, '(a)') 'Test of basis set correction for linear density variation'
      write(logUnit, '(a)') '---------------------------------------------------------'
      write(logUnit, '(a)') '                                                        |_DISABLED!'
    end if


  end subroutine Test1stVarDens

  subroutine ReadInz1(atoms, ikpt, iqpt, ikpq, iDatom, nobd, nv, z1nG)

    use m_types, only : t_atoms
    use mod_juPhonUtils, only : fopen, fclose
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms

    ! Scalar parameter
    integer,                       intent(in)  :: ikpt
    integer,                       intent(in)  :: iqpt
    integer,                       intent(in)  :: ikpq
    integer,                       intent(in)  :: iDatom

    ! Array parameter
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: nv(:, :)
    complex,                       intent(out) :: z1nG(:, :, :)

    ! Scalar variables
    integer                                    :: iband
    integer                                    :: idir
    integer                                    :: iBas
    logical                                    :: st

    ! Array variables
    character(len=:), allocatable              :: filename
    character(len=15)                          :: filenameTemp


    z1nG(:, :, :) = cmplx(0.0, 0.0)
    if (ikpt < 10 .and. iqpt >= 10 ) then
      write(filenameTemp, '(a,i1,a,i0,a,i2.2)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
    else if( ikpt >= 10 .and. iqpt < 10) then
      write(filenameTemp, '(a,i1,a,i2.2,a,i0)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
    else if( ikpt < 10 .and. iqpt < 10) then
      write(filenameTemp, '(a,i1,a,i2.2,a,i2.2)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
    else
      write(filenameTemp, '(a,i1,a,i0,a,i0)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
    end if
    filename = trim(filenameTemp)
    inquire (file=filename, exist=st)
    if ( st ) then
      call fopen( 1000, name=filename, status='old', action='read', form='unformatted')
      rewind(1000)
      ! ATTENTION: Do not change order of loops, otherwise wrong input!
      do iband = 1, nobd(ikpt, 1)
        do idir = 1, 3
          do iBas = 1, nv(1, ikpq) + atoms%nlotot
            read(1000) z1nG(iBas, iband, idir)
          end do ! iBas
        end do ! idir
      end do ! iband
      call fclose(1000)
    else
      call juDFT_warn( filename//' not found! Its z1 is set to zero.', calledby='readInz1', hint='Perform calculation to gather&
        & required z1.' )
    end if

  end subroutine ReadInz1

  ! Try to generate the density unperturbed from FLEUR
  subroutine testRho1IRroutines( cell, input, results, stars, kpts, dimens, qpts, nobd, nv, GbasVec, z, &
      & rho0IRFleur, ne, mapGbas, logUnit )

#include "cppmacro.h"

    use m_jpDens1stVar, only : calcRho1IRValDS, calcRho1IRValFFT
    use m_types
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_jpPotDensHelper, only : convertStar2G, genPotDensGvecs

    implicit none

    ! Type parameter
    type(t_cell),                 intent(in) :: cell
    type(t_input),                intent(in) :: input
    type(t_results),              intent(in) :: results
    type(t_stars),                intent(in) :: stars
    type(t_kpts),                 intent(in) :: kpts
    type(t_dimension),            intent(in) :: dimens
    type(t_kpts),                 intent(in) :: qpts

    ! Local scalar
    integer,                      intent(in) :: logUnit
    integer,                      intent(in) :: nobd(:, :)
    integer,                      intent(in) :: nv(:, :)

    ! Array parameter
    integer,                      intent(in) :: GbasVec(:, :)
    MCOMPLEX,                     intent(in) :: z(:, :, :, :)
    complex,                      intent(in) :: rho0IRFleur(:, :)
    integer,                      intent(in) :: ne(:)
    integer,                      intent(in) :: mapGbas(:, :, :) !todo rename ilst to this everywhere

    integer                                  :: ikpt
    integer                                  :: ngdp
    integer                                  :: ngdp2km

    complex,         allocatable             :: rho0IRfft(:)
    complex,         allocatable             :: rho0IRds(:, :)
    complex,         allocatable             :: z1Dummy(:, :, :)
    complex,         allocatable             :: rho0IRFleurG(:)
    complex,         allocatable             :: qpwF(:, :)
    integer,         allocatable             :: kpq2kPrVecDummy(:, :, :)
    integer,         allocatable             :: gdp(:, :)
    integer,         allocatable             :: gdp2Ind( : , : , : )
    integer                                  :: gdp2iLim(2, 3)

    integer                                  :: iG
    logical                                  :: passed = .true.

    !todo The G-vector set has to be built again because if gdp2Ind contains G-vectors until Gmax the Fleur density is not
    ! reproduced any more, this does only work if it goes until 2kmax, Check this further!!!
    call genPotDensGvecs(stars, cell, input, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, .true.)

    ! Try to get out unperturbed density
    allocate(rho0IRds(ngdp, 3))
    allocate(rho0IRfft(ngdp))
    allocate(rho0IRFleurG(ngdp))
    rho0IRds = 0.0
    rho0IRfft = 0.0
    rho0IRFleurG = 0.0

    allocate(qpwF(stars%n3d, 1))
    read(7800) qpwF
    call convertStar2G( qpwF(:, 1), rho0IRFleurG, stars, ngdp, gdp )

    allocate(z1Dummy(dimens%nbasfcn, maxval(nobd(:, :)), 3))
    allocate(kpq2kPrVecDummy(3, kpts%nkpt, qpts%nkpt))
    kpq2kPrVecDummy = 0
    do ikpt = 1, kpts%nkpt
      z1Dummy = cmplx(0.0, 0.0)
      z1Dummy(:, :nobd(ikpt, 1), 1) = z(:, :nobd(ikpt, 1), ikpt, 1)
      !is ne as parameter really correct here not nobd?
      call calcRho1IRValDS( cell, results, nobd(ikpt, 1), nv(1, :), ikpt, 1, ikpt, 1,&
        & GbasVec(:, :), z(:, :, ikpt, 1), z1Dummy, rho0IRds, gdp2Ind, mapGbas, gdp2iLim, kpq2kPrVecDummy )
!      call calcRho1IRValFFT( cell, stars, results, ikpt, ngdp2km, nobd(ikpt, 1), ne(ikpt), nv(1, ikpt), gdp(:, :), & ! todo gpot only to ngdp2km
!        &GbasVec(:, mapGbas(:nv(1, ikpt), ikpt, 1)), z(:, :, ikpt, 1), z1Dummy(:, :, 1), rho0IRfft(:) )
    end do

    ! we have to multiply a factor of 2 because we produce twice the density because the factor 2 stems from the product rule when deriving the density. 
    ! for the first variation of the density this factor of 2 is essential but when reproducing the density this factor of 2 is too much
    ! note also that when time-reversal symmetry is broken we have to rewrite this test we cannot exploit it any more to get the factor two.
    ! time-reversal symmetry is broken for spin-orbit-coupling in combination with magnetism
    if (.true.) then
      do iG = 1, ngdp
        write(8010, '(i5,1x,2(f15.8))') iG, rho0IRds(iG, 1)
        write(8011, '(i5,1x,2(f15.8))') iG, 2 * rho0IRFleurG(iG)
      end do
    end if
    write(logUnit, *)
    write(logUnit, '(a)') 'Test reproducing interstitial unperturbed density of FLEUR with juPhon double sum method routine for the linear variation of the density'
    write(logUnit, '(a)') '----------------------------------------------------------------------------------------------------------------------------------------'
    write(logUnit, *)

    do iG = 1, ngdp
      if ( abs(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG)) > 9e-6 )  then
        if( abs(rho0IRds(iG, 1)) < 1e-12) then
          write (logUnit, '(a)')'                                                                                                                                       |__failed'
          write (logUnit, '(a)') 'ToBeChecked   Reference   Absolute Difference'
          write (logUnit, '(5f15.8)') rho0IRds(iG, 1), 2 * rho0IRFleurG(iG), abs(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG)), (abs(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG)) / (2 * rho0IRFleurG(iG)))
          passed = .false.
        else if (abs(real(rho0IRds(iG, 1)) < 1e-12)) then
          if ( (abs(aimag(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG))) / aimag(2 * rho0IRFleurG(iG)) > 9e-4)) then
            write (logUnit, '(a)')'                                                                                                                                       |__failed'
            write (logUnit, '(a)') 'ToBeChecked   Reference   Absolute Difference   Relative Difference'
            write (logUnit, '(5f15.8)') rho0IRds(iG, 1), 2 * rho0IRFleurG(iG), abs(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG)), (abs(aimag(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG))) / (aimag(2 * rho0IRFleurG(iG))))
            passed = .false.
            call juDFT_warn('Interstitial density from double sum method is not consistent with FLEUR density.', calledby='testRho1IRroutines', hint='Check logfile.')
          end if
        else
          if ( (abs(real(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG))) / real(2 * rho0IRFleurG(iG)) > 9e-4) ) then
            write (logUnit, '(a)')'                                                                                                                                       |__failed'
            write (logUnit, '(a)') 'ToBeChecked   Reference   Absolute Difference   Relative Difference'
            write (logUnit, '(5f15.8)') rho0IRds(iG, 1), 2 * rho0IRFleurG(iG), abs(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG)), (abs(rho0IRds(iG, 1) - 2 * rho0IRFleurG(iG)) / (2 * real(rho0IRFleurG(iG))))
            passed = .false.
            call juDFT_warn('Interstitial density from double sum method is not consistent with FLEUR density.', calledby='testRho1IRroutines', hint='Check logfile.')
          end if
        end if
        end if
    end do ! iG

!    if ( any( abs(rho0IRds(:, 1) - 2 * rho0IRFleurG(:)) > 4e-5 ) ) then
!      write (logUnit, '(a)')'                                                                                                                                       |__failed'
!      call juDFT_warn('Interstitial density from double sum method is not consistent with FLEUR density.', calledby='testRho1IRroutines', hint='Fix bug.')
!    else
!      write (logUnit, '(a)')'                                                                                                                                       |__passed'
!    end if
    if (passed) then
      write (logUnit, '(a)')'                                                                                                                                       |__passed'
    end if


    !todo the FFT density method has to be activated at some point.
!    if ( any( abs(rho0IRfft(:) -  2 * rho0IRFleurG(:)) > 1e-9 ) ) then
!      call juDFT_error('Interstitial density from FFT method is not consistent with FLEUR density.', calledby='operateInput', hint='Fix bug.')
!    end if
!    write(logUnit, '(a)') 'Test reproducing interstitial unperturbed density of FLEUR with juPhon FFT method routine for the linear variation of the density'
!    write(logUnit, '(a)') '---------------------------------------------------------------------------------------------------------------------------------'
!    write(logUnit, '(a)') '                                                                                                                                |__passed'
!    write(logUnit, *)

  end subroutine testRho1IRroutines

  ! Put in acof and bcof and ccof and try to create the unperturbed density from FLEUR
  subroutine testRho1MTroutine(atoms, stars, dimens, sym, cell, kpts, usdus, input, results, lathar, mlh_atom, clnu_atom, nmem_atom, nv, ilst, GbasVec, ne, z, kveclo, rbas1, rbas2, ilo2p, nobd, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, rho0IR, rho0MTFleur, logUnit)

#include "cppmacro.h"

    use m_types, only : t_atoms, t_stars, t_dimension, t_sym, t_cell, t_kpts, t_usdus, t_input, t_results, t_noco, t_sphhar
    use m_abcof
    use m_od_types, only : od_inp, od_sym
    use m_jpDens1stVar, only : calcVarZContrRho1MT, multRadSolVzcRho1MT
    use m_juDFT_NOstopNO, only : juDFT_warn
    use m_JPConstants, only : fpi

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in)    :: atoms
    type(t_stars),                  intent(in)    :: stars
    type(t_dimension),              intent(in)    :: dimens
    type(t_sym),                    intent(in)    :: sym
    type(t_cell),                   intent(in)    :: cell
    type(t_kpts),                   intent(in)    :: kpts
    type(t_usdus),                  intent(in)    :: usdus
    type(t_input),                  intent(in)    :: input
    type(t_results),                intent(in)    :: results
    type(t_sphhar),                 intent(in)    :: lathar

    integer, intent(in) :: ngdp
    integer, intent(in) :: logUnit
    integer, intent(in) :: ngdp2km
    ! Array parameter
    integer,                        intent(in)    :: nv(:, :)
    integer,                        intent(in)    :: ilst(:, :, :)
    integer,                        intent(in)    :: GbasVec(:, :)
    integer,                        intent(in)    :: ne(:)
    MCOMPLEX,                       intent(in)    :: z(:, :, :, :)
    integer,                        intent(in)    :: kveclo(:,:)
    real,                           intent(in)    :: rbas1(:,:,:,:,:)
    real,                           intent(in)    :: rbas2(:,:,:,:,:)
    integer,                        intent(in)    :: ilo2p(:, :)
    integer,                      intent(in) :: nobd(:, :)
    real,                        intent(in)    :: rho0MTFleur(:, 0:, :, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    integer,            intent(in) :: nmem_atom(0:, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    integer,                      intent(in) :: gdp2iLim(2, 3)
    integer,         intent(in)    :: gdp2Ind( :, :, : )
    complex,                      intent(in) :: rho0IR(:, :)
    integer,           intent(in) :: gdp(:, :)

    ! Local scalar variable
    integer                                       :: nmat
    integer                                       :: ikpt
    integer                                       :: iatom
    integer                                       :: itype
    integer                                       :: ieqat
    integer                                       :: ptsym
    integer                                       :: ilh
    integer                                       :: oqn_l
    integer                                       :: mqn_m
    integer                                       :: lm
    integer                                       :: imem
    integer                                       :: imesh
    integer                      :: llpd
    integer :: ii, ll, jj

    ! Local array variables
    complex,           allocatable                :: acof(:, :, :)
    complex,           allocatable                :: bcof(:, :, :)
    complex,           allocatable                :: ccof(:, :, :, :)
    complex,           allocatable                :: rho0MTjuPhon(:, :, :, :)
    complex,           allocatable                :: rho0MTFleurSpH(:, :, :, :)
    complex,     allocatable        :: uu(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: du(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: ud(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: dd(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: aclo(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: bclo(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: cclo(:,:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: uunmt(:,:,:,:,:)  ! The allocation might be better a bit more outside
    complex,     allocatable        :: udnmt(:,:,:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: dunmt(:,:,:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: ddnmt(:,:,:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: uunmtFleur(:,:,:)  ! The allocation might be better a bit more outside
    real,     allocatable        :: udnmtFleur(:,:,:)! The allocation might be better a bit more outside
    real,     allocatable        :: dunmtFleur(:,:,:)! The allocation might be better a bit more outside
    real,     allocatable        :: ddnmtFleur(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: uunmtFleurSH(:,:,:)  ! The allocation might be better a bit more outside
    complex,     allocatable        :: udnmtFleurSH(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: dunmtFleurSH(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: ddnmtFleurSH(:,:,:)! The allocation might be better a bit more outside
    real,     allocatable        :: uunmtHelp(:,:,:)! The allocation might be better a bit more outside
    real,     allocatable        :: ddnmtHelp(:,:,:)! The allocation might be better a bit more outside
    real,     allocatable        :: udnmtHelp(:,:,:)! The allocation might be better a bit more outside
    real,     allocatable        :: dunmtHelp(:,:,:)! The allocation might be better a bit more outside
    complex,     allocatable        :: acnmt(:,:,:,:,:) ! The allocation might be better a bit more outside
    complex,     allocatable        :: bcnmt(:,:,:,:,:) ! The allocation might be better a bit more outside
    complex,     allocatable        :: ccnmt(:,:,:,:,:)  ! The allocation might be better a bit more outside
    real,    allocatable :: acnmtFleur(:, :, :, :, :)
    real,    allocatable :: bcnmtFleur(:, :, :, :, :)
    real,    allocatable :: ccnmtFleur(:, :, :, :, :)
    complex,    allocatable :: acnmtFleurSH(:, :, :, :)
    complex,    allocatable :: bcnmtFleurSH(:, :, :, :)
    complex,    allocatable :: ccnmtFleurSH(:, :, :, :)
    real,                        allocatable    :: rho0ValFleur(:, :, :, :)
    integer, allocatable :: ngoprI(:)
    complex,         allocatable             :: rho0IRds(:, :)


    ! Local type variables !todo beware maybe not take them from fleur_init might be dangerous
    type(t_noco)                                  :: noco
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    integer :: ilo, jlo, mqnLH_m, oqnLH_l


   ! write(1006, '(f15.8)') results%w_iks
   ! call testSpecialz1(atoms, stars, sym, cell, dimens, kpts,  input, usdus, results, nv, ne, nobd, gdp, ngdp, ngdp2km, ilst, GbasVec, z, kveclo, ilo2p, gdp2Ind, gdp2iLim, rho0IR, rbas1, rbas2)
   ! write(*, *) 'test finished'
   ! NOstopNO
    !to be substituted by dimens%lmd, when lm starts from 0
    allocate(rho0MTFleurSpH(atoms%jmtd, 0:(atoms%lmaxd + 1)**2, atoms%nat, input%jspins))
    allocate(rho0MTjuPhon(atoms%jmtd, 0:(atoms%lmaxd + 1 )**2, atoms%nat, 3))
    allocate( noco%alph(atoms%ntype), noco%beta(atoms%ntype) ) ! Dummy variables to run abcof
    allocate( uu( 0 : atoms%lmaxd, atoms%nat, 3 ), du( 0 : atoms%lmaxd, atoms%nat, 3 ), dd( 0 : atoms%lmaxd, atoms%nat, 3 ), ud( 0 : atoms%lmaxd, atoms%nat, 3 ) )
    allocate( aclo(atoms%nlod, atoms%nat, 3), bclo(atoms%nlod, atoms%nat, 3), cclo(atoms%nlod, atoms%nlod, atoms%nat, 3) )
    allocate( uunmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & ddnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd,  (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & udnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd,  (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & dunmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd,  (atoms%lmaxd + 1)**2, atoms%nat, 3) )
    allocate( acnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &bcnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &ccnmt(atoms%nlod, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3) )
    allocate(ngoprI(atoms%nat))
    ! This method is also in hsohelp
    ngoprI(:) = 1 ! The coefficients are multiplied with unity matrix sothat we have a global coordinate system and the coefficients are not rotated, this is equal to using abrot

    rho0MTjuPhon = cmplx(0.0, 0.0)
    uu = 0.0
    du = 0.0
    ud = 0.0
    dd = 0.0
    aclo = 0.0
    bclo = 0.0
    cclo = 0.0
    uunmt = 0.0
    ddnmt = 0.0
    udnmt = 0.0
    dunmt = 0.0
    acnmt = 0.0
    bcnmt = 0.0
    ccnmt = 0.0

    do ikpt = 1, kpts%nkpt
!Note: allocation and deallocation has to be made within the loop because the routines are fitted to changing sizes of the acof/bcof/ccof arrays
      allocate( acof(nobd(ikpt, 1), 0:dimens%lmd, atoms%nat), bcof(nobd(ikpt, 1), 0:dimens%lmd, atoms%nat), &
        &ccof(-atoms%llod:atoms%llod, nobd(ikpt, 1), atoms%nlod, atoms%nat) )
      acof = cmplx(0.0, 0.0)
      bcof = cmplx(0.0, 0.0)
      ccof = cmplx(0.0, 0.0)
      nmat = nv(1, ikpt) + atoms%nlotot
      ! We have to produce unrotated coefficients
      call abcof ( atoms%lmaxd, atoms%ntype, dimens%neigd, nobd(ikpt, 1), atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
        & dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), &
        & GbasVec(1, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), GbasVec(2, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
        & GbasVec(3, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv(:, ikpt),  nmat, nobd(ikpt, 1), z(:, :, ikpt, 1), usdus%us(:, :, 1), &
        & usdus%dus(:, :, 1), usdus%uds(:, :, 1), usdus%duds(:, :, 1), usdus%ddn(:, :, 1), atoms%invsat, sym%invsatnr, &
        & usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, noco%alph, noco%beta, noco%qss, kveclo(:, ikpt), odi, ods, &
        & acof, bcof, ccof)


      call calcVarZContrRho1MT( acof, bcof, ccof, acof, bcof, ccof, uu(:, :, 1), du(:, :, 1), dd(:, :, 1), ud(:, :, 1), aclo(:, :, 1), bclo(:, :, 1), cclo(:, :, :, 1), uunmt(:, :, :, :, 1), udnmt(:, :, :, :, 1), dunmt(:, :, :, :, 1), ddnmt(:, :, :, :, 1), acnmt(:, :, :, :, 1), &
        &bcnmt(:, :, :, :, 1), ccnmt(:, :, :, :, 1), atoms, nobd(ikpt, 1), rbas1, rbas2, 2 * results%w_iks(:nobd(ikpt, 1), ikpt, 1), ilo2p )

      deallocate(acof, bcof, ccof)
    end do

    ! Multiply radial solutions and add final contribution to rho1MT
    call multRadSolVzcRho1MT( atoms, aclo, bclo, cclo, acnmt, bcnmt, ccnmt, rbas1, rbas2, uu, du, ud, dd, uunmt, udnmt, dunmt, ddnmt, ilo2p, rho0MTjuPhon )

    allocate(rho0ValFleur(atoms%jmtd, 0: lathar%nlhd, atoms%ntype, input%jspins))
    rho0ValFleur = 0
    read(1040) rho0ValFleur
    iatom = 0
    rho0MTFleurSpH = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = oqn_l * (oqn_l + 1) + mqn_m
            do imesh = 1, atoms%jri(itype)
              rho0MTFleurSpH(imesh, lm, iatom, 1) = rho0MTFleurSpH(imesh, lm, iatom, 1) + rho0ValFleur(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom) / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do
    end do
    !todo this seems a bit strange here with the factor of sqrt fpi, is it in the zs then we have to change the other tests?

    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
  !      do imesh = 1, atoms%jri(iatom)
  !        rho0MTjuPhon(imesh, 0, iatom, 1) = rho0MTjuPhon(imesh, 0, iatom, 1) / sqrt(fpi)
  !      end do
        do lm = 0, (atoms%lmaxd + 1)**2
    !      do imesh = 1, atoms%jri(iatom)
    !        rho0MTjuPhon(imesh, lm, iatom, 1) = rho0MTjuPhon(imesh, lm, iatom, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
    !      end do
        end do
      end do
    end do
    if (.false.) then
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do lm = 0, atoms%lmax(itype) * (atoms%lmax(itype) + 1)
              do imesh = 1, atoms%jri(itype)
                write(1082, '(2i8,2f15.8)') imesh, lm, abs(2 * rho0MTFleurSpH(imesh, lm, iatom, 1) - rho0MTjuPhon(imesh, lm, iatom, 1))
                if (2 * abs(real(rho0MTFleurSpH(imesh, lm, iatom, 1))) < 1e-8 .and. 2 * abs(aimag(rho0MTFleurSpH(imesh, lm, iatom, 1))) > 1e-8) then
                write(1080, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, 2 * abs(real(rho0MTFleurSpH(imesh, lm, iatom, 1))), 2 * aimag(rho0MTFleurSpH(imesh, lm, iatom, 1))
              else if (2 * abs(real(rho0MTFleurSpH(imesh, lm, iatom, 1))) > 1e-8 .and. 2 * abs(aimag(rho0MTFleurSpH(imesh, lm, iatom, 1))) < 1e-8) then

                write(1080, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, 2 * real(rho0MTFleurSpH(imesh, lm, iatom, 1)), 2* abs(aimag(rho0MTFleurSpH(imesh, lm, iatom, 1)))
              else if (2 * abs(real(rho0MTFleurSpH(imesh, lm, iatom, 1))) < 1e-8 .and. 2* abs(aimag(rho0MTFleurSpH(imesh, lm, iatom, 1))) < 1e-8) then
                   write(1080, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, 2 * abs(real(rho0MTFleurSpH(imesh, lm, iatom, 1))), 2* abs(aimag(rho0MTFleurSpH(imesh, lm, iatom, 1)))
              else
                   write(1080, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, 2 * real(rho0MTFleurSpH(imesh, lm, iatom, 1)), 2 * aimag(rho0MTFleurSpH(imesh, lm, iatom, 1))
              end if

                if (abs(real(rho0MTjuPhon(imesh, lm, iatom, 1))) < 1e-8 .and. abs(aimag(rho0MTjuPhon(imesh, lm, iatom, 1))) > 1e-8) then
                write(1081, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, abs(real(rho0MTjuPhon(imesh, lm, iatom, 1))), aimag(rho0MTjuPhon(imesh, lm, iatom, 1))
              else if (abs(real(rho0MTjuPhon(imesh, lm, iatom, 1))) > 1e-8 .and. abs(aimag(rho0MTjuPhon(imesh, lm, iatom, 1))) < 1e-8) then

                write(1081, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, real(rho0MTjuPhon(imesh, lm, iatom, 1)), abs(aimag(rho0MTjuPhon(imesh, lm, iatom, 1)))
              else if (abs(real(rho0MTjuPhon(imesh, lm, iatom, 1))) < 1e-8 .and. abs(aimag(rho0MTjuPhon(imesh, lm, iatom, 1))) < 1e-8) then
                   write(1081, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, abs(real(rho0MTjuPhon(imesh, lm, iatom, 1))), abs(aimag(rho0MTjuPhon(imesh, lm, iatom, 1)))
              else
                   write(1081, '(3(i5,1x),2(f15.8,1x))') imesh, lm, iatom, real(rho0MTjuPhon(imesh, lm, iatom, 1)), aimag(rho0MTjuPhon(imesh, lm, iatom, 1))
              end if
              end do ! imesh
            end do ! lm
          end do ! ieqat
        end do ! itype
    end if

    write(logUnit, '(a)') 'Test reproducing muffin-tin unperturbed density of FLEUR with juPhon FFT method routine for the linear variation of the density'
    write(logUnit, '(a)') '-------------------------------------------------------------------------------------------------------------------------------'
    write(logUnit, *)

    if ( any( abs(2 * rho0MTFleurSpH(:, :, :, 1) - rho0MTjuPhon(:, :, :, 1)) > 9e-6 ) ) then
      write(logUnit, '(a)') '                                                                                                                                |__failed'
      call juDFT_warn('Muffin-tin density is not consistent with FLEUR density.', calledby='testRho1MTroutine', hint='Fix bug.')
    else
      write(logUnit, '(a)') '                                                                                                                                |__passed'
    end if

  end subroutine testRho1MTroutine

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Similiar to the CheckDensnPot Test, the continuity of the wavefunction at the MT boundary is tested.
  !>
  !> @details
  !> For noPotsCon random points at the MT sphere the continuity of the wavefunctions are tested.
  !>
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[out] kpts       : K-points type, see types.f90.
  !> @param[out] sym        : Symmetries type, see types.f90.
  !> @param[out] dimens     : Dimension type, see types.f90.
  !> @param[out] usdus      : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] noPotsCon  : Number of points for which continuity test should be performed.
  !> @param[out] GbasVec    : G-basis vectors
  !> @param[out] ilst       : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                          pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] nv         : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] ne         : Number of eigenvalues per k-point.
  !> @param[out] z          : Kohn-Sham eigenvectors.
  !> @param[out] nRadFun    : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] kveclo     : Basis G-vectors of local orbitals.
  !> @param[out] rbas1      : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2      : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] iloTable   : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkz1ph0Cont( atoms, cell,  sym, dimens, usdus, noPtsCon, nobd, bkpt, ikpt, GbasVec, ilst, nv, z1, nRadFun, kveclo, rbas1, &
    & rbas2, iloTable, logUnit, nobdSum, idir )

#include "cppmacro.h"

    use m_JPConstants
    use m_types, only : t_noco, t_atoms, t_cell, t_kpts, t_sym, t_dimension, t_usdus
    use m_abcof
    use m_cotra
    use m_ylm_old
    use m_od_types, only: od_inp, od_sym
    use m_juDFT_time, only : TimeStart, TimeNOstopNO

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in) :: atoms
    type(t_cell),              intent(in) :: cell
    type(t_sym),               intent(in) :: sym
    type(t_dimension),         intent(in) :: dimens
    type(t_usdus),             intent(in) :: usdus

    ! Scalar parameter
    integer,                   intent(in) :: noPtsCon !number of points to check continuity
    integer,                   intent(in) :: nobd
    integer,                   intent(in) :: logUnit
    integer,                   intent(in)  :: ikpt

    ! Array parameter
    real,                      intent(in) :: bkpt(:)
    integer,                   intent(in) :: GbasVec(:, :)
    integer,                   intent(in) :: ilst(:, :, :)
    integer,                   intent(in) :: nv(:, :)
    MCOMPLEX,                  intent(in) :: z1(:,:)
    integer,                   intent(in) :: nRadFun(0:, :)
    integer,                   intent(in) :: kveclo(:, :)
    real,                      intent(in) :: rbas1(:,:,0:,:,:)
    real,                      intent(in) :: rbas2(:,:,0:,:,:)
    integer,                   intent(in) :: iloTable(:, 0:,: )
    integer, intent(in) :: nobdSum
    integer, intent(in) :: idir

    ! Local type variables
    type(t_noco)                          :: noco_type
    type (od_inp)                         :: odi_type !Both types, this and the one beneath have to be taken from od_types module
    type (od_sym)                         :: ods_type

    ! Local scalar variables
    integer                               :: jspin
    integer                               :: nmat
    integer                               :: itype
    integer                               :: irandPt
    real                                  :: euclidNorm
    integer                               :: iatom
    integer                               :: l
    integer                               :: m
    integer                               :: ieqat
    integer                               :: lm
    integer                               :: lmp
    complex                               :: cdum
    integer                               :: p
    integer                               :: ilo
    integer                               :: iband
    complex                               :: sumOfMTFuncVal
    complex                               :: sumInterst
    integer                               :: iG
    integer                               :: symOpr_temp
    integer                               :: imatmulc, imatmulr
    integer                               :: l_temp
    integer                               :: lm_temp,maxlmp
    real                                  :: av1, rms1, dms1, av2, rms2, dms2
    complex                               :: expFunc

    ! Local array variables
    complex,      allocatable             :: acof(:, :, :)
    complex,      allocatable             :: bcof(:, :, :)
    complex,      allocatable             :: ccof(:, :, :, :)
    real                                  :: Gext_temp(3)
    real                                  :: rvecExt(3)
    real,         allocatable             :: randPtsCart(:, :, :)
    complex,     allocatable              :: coeffMT(:, :)
    real                                  :: GridPtsMT_Temp(3)
    real                                  :: rotatedVector(3)
    complex,     allocatable              :: ylm(:)
    real,        allocatable              :: mismatch(:)
    real                                  :: currPnt(3)
    integer, allocatable :: ngoprI(:)


    allocate( noco_type%alph(atoms%ntype), noco_type%beta(atoms%ntype) ) !Up to now those variables are only of dummy character
    allocate( randPtsCart(3, noPtsCon, atoms%ntype) )
    allocate( mismatch(atoms%nat) )
    allocate(ngoprI(atoms%nat))

    ! Generate random points
    call TimeStart('randomPoints')
    do itype = 1, atoms%ntype
      do irandPt = 1, noPtsCon
        call random_number( randPtsCart(:, irandPt, itype) )
        euclidNorm = norm2( randPtsCart(:, irandPt, itype) )
        randPtsCart(:, irandPt, itype) = randPtsCart(:, irandPt, itype) / euclidNorm * atoms%rmt(itype)
      end do
    end do
    call TimeNOstopNO('randomPoints')

    maxlmp = maxval( (/ (sum( (/ ((2*l+1)* nRadFun(l,itype), l = 0, atoms%lmax(itype)) /) ), itype=1, atoms%ntype) /) )

    write(*, *) logUnit
    write(logUnit,*) '---------------------------------------------------------------'
    write(logUnit,*) 'Check wave functions for continuity at MT boundary:'
    write(logUnit,*)
    write(logUnit,*) '   jspin    ikpt    band-averaged abs. discontinuity per atom'
    ngoprI(:) = 1
    do jspin =  1, dimens%jspd

        mismatch = 0

        nmat = nv(jspin, ikpt) + atoms%nlotot !dimension of Hamiltonian

        allocate(acof(nobd, 0:dimens%lmd, atoms%nat), bcof(nobd, 0:dimens%lmd, atoms%nat), &
          &ccof(-atoms%llod:atoms%llod, nobd, atoms%nlod, atoms%nat))
        acof = cmplx(0.0, 0.0)
        bcof = cmplx(0.0, 0.0)
        ccof = cmplx(0.0, 0.0)

        call Abcof( atoms%lmaxd, atoms%ntype, dimens%neigd, nobd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
          &dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, &
          &atoms%neq, atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, bkpt, GbasVec(1, ilst(:nv(jspin, ikpt), ikpt, jspin)), &
          &GbasVec(2, ilst(:nv(jspin, ikpt), ikpt, jspin)), GbasVec(3, ilst(: nv(jspin, ikpt), ikpt, jspin)), nv(:, ikpt), nmat, &
          &nobd, z1(:, :), usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, &
          &usdus%ulos, usdus%uulon, usdus%dulon, usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, noco_type%l_noco, &
          &noco_type%l_ss, jspin, noco_type%alph, noco_type%beta, noco_type%qss, kveclo(:, ikpt), odi_type ,ods_type, acof, bcof, ccof )


        iatom   = 0
        do itype = 1, atoms%ntype
          allocate( ylm( ( atoms%lmax(itype) +  1 )**2 ) )
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            allocate( coeffMT(nobd, maxlmp) )
            coeffMT = 0
            lmp   =   0
            lm    =  -1
            cdum  = -iu
            do l = 0, atoms%lmax(itype)
              cdum = cdum * iu
              do m = -l, l
                lm = lm + 1
                !p = 1
                lmp = lmp + 1
                coeffMT(:nobd, lmp) = cdum * acof(:nobd, lm, iatom)
                !p = 2
                lmp = lmp + 1
                coeffMT(:nobd, lmp) = cdum * bcof(:nobd, lm, iatom)
                !LOs
                do p = 3, nRadFun(l, itype)
                  ilo = iloTable(p, l, itype)
                  lmp = lmp + 1
                  coeffMT(:nobd, lmp) = cdum * ccof(m, :nobd, ilo, iatom)
                end do
              end do
            end do

            do irandPt = 1, noPtsCon

              ! Determine current point coordinates in local coordinate system
              currPnt(:) = randPtsCart(:, irandPt, itype) + atoms%pos(:, iatom)
              !Interstitial part, part coming from the local orbitals is 0 in interstitial
              sumInterst = 0
              do iG = 1, nv(jspin, ikpt)

                !todo rewrite
                call cotra3(bkpt(:)+GbasVec(:, ilst(iG, ikpt, jspin)),Gext_temp,cell%bmat)
                expFunc = exp( iu *  dot_product( Gext_temp, currPnt ) )

                do iband = 1,nobd
                  ! this is no linear run in storage but instead of ne * nv matrix vector operations and dot products we only do nv
                  ! sum up all contributions of the wavefunctions looping over the IR basis functions
                  sumInterst = sumInterst + z1(iG, iband)  * expFunc
                end do
              end do
              sumInterst = sumInterst *  cell%vol**(-0.5)

              !Evaluate MT sphere function at same random points
!              symOpr_temp = atoms%ngopr(iatom) ! symmetry operation mapping local to global coordinate system
              symOpr_temp = 1!atoms%ngopr(iatom) ! symmetry operation mapping local to global coordinate system
             ! if (symOpr_temp /= 1) then
             !   if ( atoms%invsat(iatom) == 2 ) then
             !     rvecExt = - randPtsCart(:, irandPt, itype)
             !   else
             !     ! rotate vector with fitting symmetry operation if not correct vector
             !     call cotra1(randPtsCart(:, irandPt, itype), GridPtsMT_Temp(:), cell%bmat) !allocate and deallocate randPtsGrid
             !     rotatedVector = matmul( sym%mrot(:,:,symOpr_temp), GridPtsMT_Temp )
             !     call cotra0(rotatedVector, rvecExt, cell%amat)
             !   end if
             ! else
                rvecExt = randPtsCart(:, irandPt, itype)
             ! end if
              call Ylm4( atoms%lmax(itype), rvecExt, ylm )

              ! die kleine Komponente ist eher uninteressant in diesem Zusammenhang, da nur die grosse Komponente an die Planewave gematcht ist und wir daher nur deren Stetigkeit testen kann
              sumOfMTFuncVal = 0
              lmp = 0
              do l = 0, atoms%lmax(itype)
                l_temp = 1 +  l * (l + 1)
                do m =  - l , l
                  lm_temp = l_temp +  m
                  do p = 1, nRadFun(l, itype)
                    lmp = lmp + 1
                    do iband = 1, nobd
                      !sum up all contributions of the wavefunctions (hidden in the matching coefficients coeffMT) looping over all basis
                      !functions of the MT
                      sumOfMTFuncVal = sumOfMTFuncVal + coeffMT(iband, lmp) * rbas1(atoms%jri(itype), p, l, itype, jspin) *  &
                        &ylm(lm_temp) / atoms%rmt(itype)
                    end do
                  end do
                end do
              end do

            end do ! irandPt
            ! Compare the IR value of the wavefunction with the MT value of the wavefunction
            mismatch(iatom) = mismatch(iatom) + abs( sumOfMTFuncVal - sumInterst ) /  nobdSum / noPtsCon
            deallocate(coeffMT)
          end do ! ieqat
        deallocate(ylm)
      end do ! itype
      deallocate(acof, bcof, ccof)

      !todo output should be reviewed
      write(logUnit,'(2(4x,i4),4x,3es15.8)') idir,ikpt,mismatch, mismatch/ abs(sumOfMTFuncVal) * 100
    end do  !jspin

  end subroutine checkz1ph0Cont

  ! Makes data for a 3dplot of the density
    subroutine test3DplotPotDens( atoms, stars, sym, cell, dimens, lathar, input, ngdp, clnu_atom, nmem_atom, mlh_atom, gdp, &
        &rho0IR, rho0MT, logUnit )

      use m_types
      use m_jpPlotObservables, only : plotVestPotDens0Wrap, plotVestPotDens1
      use m_jpPotDensHelper, only : convertStar2G

      implicit none

      ! Type parameter
      type(t_atoms),                              intent(in) :: atoms
      type(t_stars),                              intent(in) :: stars
      type(t_sym),                                intent(in) :: sym
      type(t_cell),                               intent(in) :: cell
      type(t_dimension),                          intent(in) :: dimens
      type(t_sphhar),                             intent(in) :: lathar
      type(t_input),                              intent(in) :: input

      integer,                                    intent(in) :: ngdp
      integer,                                    intent(in) :: logUnit

      complex,                                    intent(in) :: clnu_atom(:, 0:, :)
      integer,                                    intent(in) :: nmem_atom(0:, :)
      integer,                                    intent(in) :: mlh_atom(:, 0:, :)
      integer,                                    intent(in) :: gdp(:, :)
      complex,                                    intent(in) :: rho0IR(:, :)
      real,                                       intent(in) :: rho0MT(:, 0:, :, :)

      type(t_sliceplot)                                      :: sliceplot

      integer                                                :: itype
      integer                                                :: ieqat
      integer                                                :: iatom
      integer                                                :: ptsym
      integer                                                :: ilh
      integer                                                :: oqn_l
      integer                                                :: lm
      integer                                                :: imesh
      real                                                   :: imem
      real                                                   :: mqn_m
      real,               parameter                          :: qpt(3) = 0.0
      complex,                       allocatable             :: rho0IRpw(:, :)
      complex,                       allocatable             :: rho0MTsh(:, :, :, :)

      integer :: ii


      ! Plot unperturbed density with adapted plotdop routine
      sliceplot%plpot = .false. ! plot density and not potential
      call plotVestPotDens0Wrap( atoms, stars, sym, cell, dimens, lathar, input, sliceplot, rho0IR, rho0MT )

      ! Transform densities to plane waves and spherical harmonics
      allocate(rho0IRpw(ngdp, 1))
      allocate(rho0MTsh(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, input%jspins))

      rho0IRpw = cmplx(0.0, 0.0)
      rho0MTsh = cmplx(0.0, 0.0)

      call convertStar2G(rho0IR(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp)

 !     write(1924, *) rho0IR
 !     do ii = 1, ngdp
 !     write(1925, '(3(i5),2(es15.8))') gdp(1, ii), gdp(2, ii), gdp(3, ii), rho0IRpw(ii,1)
 !     end do
 !     NOstopNO

      rho0MTsh = 0
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ptsym = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l = lathar%llh(ilh, ptsym)
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MT(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh
        end do
      end do

      ! Call routine actually thought for variation of the density with transformed unperturbed density to test it
     ! qpt(:) = 0.0
      call plotVestPotDens1( dimens, stars, lathar, atoms, input, sym, cell,  gdp, ngdp, qpt, rho0IRpw(:, 1), rho0MTsh(:, :, :, 1), 'plot1.xsf', 1 , .TRUE., "plot_inp")

      write(logUnit, *) 
      write(logUnit, '(a)') '-----------------------------------------------------------------------------------------------------------------------------------------------------'
      write(logUnit, '(a)') 'Test reproducing the 3D unperturbed density of FLEUR with the juPhon routine for the linear variation of the density'
      write (logUnit, '(a)')'  --> Compare plot0.xsf from FLEUR (header might be in other units and differ by a factor of approximately 0.52917720859.) with plot1.xsf from juPhon'
      write (logUnit, '(a)')'  --> Find the maximal deviation executing the python script testPlotting3D.py'
      write(logUnit, '(a)') '-----------------------------------------------------------------------------------------------------------------------------------------------------'
      write(logUnit, *) 


    end subroutine test3DplotPotDens

    ! For q = 0 we assume a special z1 = -i (k + G) z0. For this choice, the density variation is supposed to be zero and the
    ! Hellmann-Feynman contributions should be almost zero except for an surface integral. This is the idea of the test.
    subroutine TestRho1BasCorr(atoms, stars, sym, cell, dimens, kpts,  input, usdus, results, nv, ne, nobd, gdp, ngdp, ngdp2km,    &
        & ilst, GbasVec, z0, kveclo, ilo2p, gdp2Ind, gdp2iLim, rho0IR, rbas1, rbas2, kpq2kPrVec, logUnit )

      use m_types
      use m_juDFT_NOstopNO, only : juDFT_warn
      use m_jPConstants, only : iu
      use m_jpDens1stVar, only : calcKdepValRho1MT, multRadSolVzcRho1MT, calcRho1IRValDS
      use m_jpPotDensHelper, only : convertStar2G

      implicit none

      ! Type parameters
      type(t_atoms),                 intent(in) :: atoms
      type(t_stars),                 intent(in) :: stars
      type(t_sym),                   intent(in) :: sym
      type(t_cell),                  intent(in) :: cell
      type(t_dimension),             intent(in) :: dimens
      type(t_kpts),                  intent(in) :: kpts
      type(t_input),                 intent(in) :: input
      type(t_usdus),                 intent(in) :: usdus
      type(t_results),               intent(in) :: results

      ! Scalar parameter
      integer,                       intent(in) :: ngdp
      integer,                       intent(in) :: logUnit
      integer,                       intent(in) :: ngdp2km
      ! Array parameters
      integer,                       intent(in) :: nv(:, :)
      integer,                       intent(in) :: ne(:)
      integer,                       intent(in) :: ilst(:, :, :)
      integer,                       intent(in) :: GbasVec(:, :)
      MCOMPLEX,                      intent(in) :: z0(:, :, :, :)
      integer,                       intent(in) :: kveclo(:,:)
      integer,                       intent(in) :: ilo2p(:, :)
      real,                          intent(in) :: rbas1(:, :, 0:, :, :)
      real,                          intent(in) :: rbas2(:, :, 0:, :, :)
      integer,                       intent(in) :: nobd(:, :)
      integer,                       intent(in) :: gdp2Ind( : , : , : )
      integer,                       intent(in) :: gdp2iLim(2, 3)
    complex,                      intent(in) :: rho0IR(:, :)
    integer,                         intent(in) :: kpq2kPrVec(:, :, :)
    complex       ,allocatable                           :: rho0IRnCt(:, :)
    integer,           intent(in) :: gdp(:, :)

      ! Local scalar variables
      integer                                   :: ikpt
      integer                                   :: iBas
      integer                                   :: ieig
      integer                                   :: idir
      integer                                   :: dispAtInd
      integer                                   :: dispType
      integer                                   :: ieqat
      integer :: iG
      logical :: passed = .true.


      ! Local array variables
      real                                      :: kExt(3)
      real                                      :: Gext(3)
      complex,           allocatable            :: z1Special(:, :, :, :)
      complex,           allocatable            :: rho1MT(:, :, :, :, :)
      complex,           allocatable            :: uu(:,:,:)
      complex,           allocatable            :: du(:,:,:)
      complex,           allocatable            :: ud(:,:,:)
      complex,           allocatable            :: dd(:,:,:)
      complex,           allocatable            :: aclo(:,:,:)
      complex,           allocatable            :: bclo(:,:,:)
      complex,           allocatable            :: cclo(:,:,:,:)
      complex,           allocatable            :: uunmt(:,:,:,:,:)
      complex,           allocatable            :: udnmt(:,:,:,:,:)
      complex,           allocatable            :: dunmt(:,:,:,:,:)
      complex,           allocatable            :: ddnmt(:,:,:,:,:)
      complex,           allocatable            :: acnmt(:,:,:,:,:)
      complex,           allocatable            :: bcnmt(:,:,:,:,:)
      complex,           allocatable            :: ccnmt(:,:,:,:,:)
    complex,           allocatable   :: rho1IRDS(:, :, :)
    complex,                      allocatable :: rho0IRG(:)
    complex,           allocatable   :: gradRhoIR(:, :)
    real :: invsAtNr
    real, allocatable:: polVec(:, :)
    complex, allocatable :: atomDirSumIR(:)
    complex, allocatable :: atomDirSumMT(:, :)

      write(logUnit, '(a)') 'Test of basis set correction for linear density variation'
      write(logUnit, '(a)') '---------------------------------------------------------'

    allocate(rho0IRG(ngdp))
    allocate(polVec(3, atoms%nat))
 !   read(7800) rho0IRncT

      allocate( rho1IRDS(ngdp, 3, atoms%nat) )
      allocate( gradRhoIR(ngdp, 3) )
      allocate( z1Special(dimens%nbasfcn, dimens%neigd, 3, kpts%nkpt) )
      ! Construct z1 = -i (k + G) z0
      z1Special = 0
!      write(*, *) 'z1Special test have to be reviewed again'
      !TODO actually the factor 1 / atoms%nat is not required, if we for example have two atoms we can also evaluate it for every 
      ! atom itself. Think about it. it is only a question of perspective
      invsAtNr = 1.! / real(atoms%nat)
      do ikpt = 1, kpts%nkpt
        kExt(:) = matmul( cell%bmat, kpts%bk(:, ikpt) )
        do idir = 1, 3
          do iBas = 1, nv(1, ikpt)
            Gext(:) = matmul( cell%bmat, GbasVec(:, ilst(iBas, ikpt, 1)) )
            do ieig = 1, nobd(ikpt, 1) !LOs mit berücksichtigen, auch mit Gs?, + or - q?
              z1Special(iBas, ieig, idir, ikpt) = z1Special(iBas, ieig, idir, ikpt) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) ) * z0(iBas, ieig, ikpt, 1)

            end do
          end do
          do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
            Gext(:) = matmul( cell%bmat, GbasVec(:, ilst(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) )
            do ieig = 1, nobd(ikpt, 1) !LOs mit berücksichtigen, auch mit Gs?, + or - q?
              z1Special(iBas, ieig, idir, ikpt) = z1Special(iBas, ieig, idir, ikpt) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) ) * z0(iBas, ieig, ikpt, 1)
            end do
          end do
        end do
      end do

      ! Calculate the density response
      allocate( uu( 0 : atoms%lmaxd, atoms%nat, 3 ), du( 0 : atoms%lmaxd, atoms%nat, 3 ), dd( 0 : atoms%lmaxd, atoms%nat, 3 ), ud( 0 : atoms%lmaxd, atoms%nat, 3))
      allocate( uunmt( 0 : atoms%lmaxd, 0: atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
        & ddnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
        & udnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
        & dunmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3) )
      allocate( acnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
        &bcnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
        &ccnmt(atoms%nlod, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3) )
      allocate( aclo(atoms%nlod, atoms%nat, 3), bclo(atoms%nlod, atoms%nat, 3), cclo(atoms%nlod, atoms%nlod, atoms%nat, 3) )

      uu = 0.0
      du = 0.0
      ud = 0.0
      dd = 0.0
      aclo = 0.0
      bclo = 0.0
      cclo = 0.0
      uunmt = 0.0
      ddnmt = 0.0
      udnmt = 0.0
      dunmt = 0.0
      acnmt = 0.0
      bcnmt = 0.0
      ccnmt = 0.0
      ! they have to be set to zero for every atom!
      dispAtInd = 0

      allocate( rho1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat))
      rho1MT = cmplx(0.0, 0.0)
      rho1IRDS = cmplx(0.0, 0.0)
      do dispType = 1, atoms%ntype
        do ieqat = 1, atoms%neq(dispType)
          dispAtInd = dispAtInd + 1
      uu = 0.0
      du = 0.0
      ud = 0.0
      dd = 0.0
      aclo = 0.0
      bclo = 0.0
      cclo = 0.0
      uunmt = 0.0
      ddnmt = 0.0
      udnmt = 0.0
      dunmt = 0.0
      acnmt = 0.0
      bcnmt = 0.0
      ccnmt = 0.0
          do ikpt = 1, kpts%nkpt
            call calcKdepValRho1MT( atoms, dimens, sym, cell, kpts, usdus, results, ikpt, ikpt, dispAtInd, nv, ilst, GbasVec, nobd(:, 1), uu, &
              &du, ud, dd, aclo, bclo, cclo, uunmt, udnmt, dunmt, ddnmt, acnmt, bcnmt, ccnmt, z0(:, :, ikpt, 1), z1Special(:, :, :, ikpt), &
              & kveclo, rbas1, rbas2, ilo2p )
            !write(1005, '(i3, i6)') ikpt, nv(1, ikpt)
        !    write(*, *) 'nobd', ikpt, nobd(ikpt, 1)
        !    write(1005, *) ikpt, 2 * results%w_iks(:, ikpt, 1)
            do idir = 1, 3
              call calcRho1IRValDS( cell, results, nobd(ikpt, 1), nv(1, :), ikpt, 1, ikpt, idir,&
                & GbasVec(:, :), z0(:, :, ikpt, 1), z1Special(:, :, :, ikpt), rho1IRDS(:, :, dispAtInd), gdp2Ind, ilst, gdp2iLim, kpq2kPrVec)
            end do
          end do ! ikpt
          call multRadSolVzcRho1MT( atoms, aclo, bclo, cclo, acnmt, bcnmt, ccnmt, rbas1, rbas2, uu, du, ud, dd, uunmt, udnmt, dunmt, &
            &ddnmt, ilo2p, rho1MT(:, :, :, :, dispAtInd) )
        end do
      end do

      ! Determine gradient of density in interstitial
     ! write(1015, *) rho0IR(1, 1)
     ! this needs to be n3d not n3g otherwise we get an error
     !todo this does also work for core terms in the interstitial, muffin-tin not known yet
     allocate(rho0IRncT(stars%n3d, 1))
     rho0IRncT(:, :) = cmplx(0.,0.)
     rewind(7800)
     read(7800) rho0IRncT
!     write(1827,'(2(f15.8))') rho0IRncT
     call convertStar2G(rho0IRncT(:, 1), rho0IRG, stars, ngdp, gdp)

     ! do iG = 1, ngdp
     !   write(1013, '(4(i8,1x), 2(f15.8))') iG, gdp(1, iG), gdp(2, iG), gdp(3, iG), rho0IRG(iG)
     ! end do

      gradRhoIR = 0.0
      do iG = 1, ngdp2km
        Gext(:) = matmul( cell%bmat, gdp(:, iG) )
        do idir = 1, 3
          gradRhoIR(iG, idir) = gradRhoIR(iG, idir) - iu * Gext(idir)  * rho0IRG(iG)
        end do
      end do


      allocate(atomDirSumMT(atoms%jmtd, (atoms%lmaxd + 1)**2))
      allocate(atomDirSumIR(ngdp))

      ! Normalization of polarization vectors
      polVec(1, :) = 1. / real(atoms%nat) * 1.
      polVec(2, :) = 0.
      polVec(3, :) = 0.

     ! do dispAtInd = 1, atoms%nat
     !   do idir = 1, 3
     !     write(1008, *) idir, dispAtInd, polVec(idir, dispAtInd)
     !   end do
     ! end do

     ! We are only interested in the density of the displaced atom, because only that cancels away
      atomDirSumMT = cmplx(0.0, 0.0)
      do dispAtInd = 1, atoms%nat
        do idir = 1, 3
          atomDirSumMT(:, :) = atomDirSumMT(:, :) + rho1MT(:, :, dispAtInd, idir, dispAtInd) * polVec(idir, dispAtInd)
        end do
      end do

      atomDirSumIR = cmplx(0.0, 0.0)
      do dispAtInd = 1, atoms%nat
        do idir = 1, 3
          do iG = 1, ngdp
            atomDirSumIR(iG) = atomDirSumIR(iG) + rho1IRDS(iG, idir, dispAtInd) * polVec(idir, dispAtInd)
          end do
        end do
      end do

     ! do iG = 1, ngdp
     !   write(1002, '(4(i4,1x), 2(f15.8))') iG, gdp(1, iG), gdp(2, iG), gdp(3, iG), atomDirSumIR(iG)
     !   write(1003, '(4(i4,1x), 2(f15.8))') iG, gdp(1, iG), gdp(2, iG), gdp(3, iG), gradRhoIR(iG, 1)
     ! end do

     ! do iG = 1, ngdp
     !   write(1010, '(4(i4,1x))') iG, gdp(1, iG), gdp(2, iG), gdp(3, iG)
     ! end do



     ! do iG = 1, ngdp2km
     !   if (abs(gradRhoIR(iG, 1) -  atomDirSumIR(iG)) < 1e-8) then
     !     write(1007, '(4(i4,1x), 4(f15.8), i5)') iG, gdp(1, iG), gdp(2, iG), gdp(3, iG), gradRhoIR(iG, 1), atomDirSumIR(iG), stars%ig(gdp(1, iG), gdp(2, iG), gdp(3, iG))
     !   end if
     ! end do


      ! Compare gradient of density with density response
      if (any(abs(gradRhoIR(:, 1) -  atomDirSumIR(:)) > 1e-8)) then
        passed = .false.
        call juDFT_warn("z1Special Test failed for IR dir 1",calledby ="testSpecialz1")
      end if


!     dispAtInd = 0
!     do dispType = 1, atoms%ntype
!       do ieqat = 1, atoms%neq(dispType)
!         dispAtInd = dispAtInd + 1
         if (any(abs(atomDirSumMT(:, :)) > 1e-8)) then
           passed = .false.
           call juDFT_warn("z1Special Test failed for MT dir 1", calledby ="testSpecialz1")
         end if
!       end do
!     end do

      ! Normalization of polarization vectors
      polVec(1, :) = 0.
      polVec(2, :) = 1. / real(atoms%nat) * 1.
      polVec(3, :) = 0.

      ! Only if we are at the displaced atom the basis correction cancels the density away
      atomDirSumMT = cmplx(0.0, 0.0)
      do dispAtInd = 1, atoms%nat
        do idir = 1, 3
          atomDirSumMT(:, :) = atomDirSumMT(:, :) + rho1MT(:, :, dispAtInd, idir, dispAtInd) * polVec(idir, dispAtInd)
        end do
      end do

      atomDirSumIR = cmplx(0.0, 0.0)
      do dispAtInd = 1, atoms%nat
        do idir = 1, 3
          atomDirSumIR(:) = atomDirSumIR(:) + rho1IRDS(:, idir, dispAtInd) * polVec(idir, dispAtInd)
        end do
      end do


      ! Compare gradient of density with density response
      if (any(abs(gradRhoIR(:, 2) -  atomDirSumIR) > 1e-8)) then
        passed = .false.
        call juDFT_warn("z1Special Test failed for IR dir 2",calledby ="testSpecialz1")
      end if


     dispAtInd = 0
     do dispType = 1, atoms%ntype
       do ieqat = 1, atoms%neq(dispType)
         dispAtInd = dispAtInd + 1
         if (any(abs(atomDirSumMT(:, :)) > 1e-8)) then
           passed = .false.
           call juDFT_warn("z1Special Test failed for MT dir 2", calledby ="testSpecialz1")
         end if
       end do
     end do

      ! Normalization of polarization vectors
      polVec(1, :) = 0.
      polVec(2, :) = 0.
      polVec(3, :) = 1. / real(atoms%nat) * 1.

      atomDirSumMT = cmplx(0.0, 0.0)
      do dispAtInd = 1, atoms%nat
        do idir = 1, 3
          atomDirSumMT(:, :) = atomDirSumMT(:, :) + rho1MT(:, :, dispAtInd, idir, dispAtInd) * polVec(idir, dispAtInd)
        end do
      end do

      atomDirSumIR = cmplx(0.0, 0.0)
      do dispAtInd = 1, atoms%nat
        do idir = 1, 3
          atomDirSumIR(:) = atomDirSumIR(:) + rho1IRDS(:, idir, dispAtInd) * polVec(idir, dispAtInd)
        end do
      end do


      ! Compare gradient of density with density response
      if (any(abs(gradRhoIR(:, 3) -  atomDirSumIR) > 1e-8)) then
        passed = .false.
        call juDFT_warn("z1Special Test failed for IR dir 3",calledby ="testSpecialz1")
      end if


     dispAtInd = 0
     do dispType = 1, atoms%ntype
       do ieqat = 1, atoms%neq(dispType)
         dispAtInd = dispAtInd + 1
         if (any(abs(atomDirSumMT(:, :)) > 1e-8)) then
           passed = .false.
           call juDFT_warn("z1Special Test failed for MT dir 3", calledby ="testSpecialz1")
         end if
       end do
     end do

    if (passed) then
      write (logUnit, '(a)')  '                                                        |__passed'
    else
      write (logUnit, '(a)')  '                                                        |__failed'
    end if

    end subroutine TestRho1BasCorr

    ! Plots the gradient of the density and compares it to rho1 and is actually recycled from the routine that compares the gradient from weinert with tht numerical gradient of the effective potential
  subroutine checkGradDens0s(symT, cellT, inputT, sphharT, starsT, atomsT, vGradCoul0IR, ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MT, gradVrjuPhon, paPoXCo, paPoYCo, paPoZCo, harSw, extSw, xcSw, vpwStar, vrFleur)

    use m_types
    use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
    use m_gaunt
    use m_ylm_old
    use m_cotra
    use m_juDFT_NOstopNO
    use mod_juPhonUtils, only : fopen, fclose, calcGrFinLH
    use m_jpPotDensHelper, only : convertStar2G

    implicit none

    type(t_sym),        intent(in) :: symT
    type(t_cell),       intent(in) :: cellT
    type(t_sphhar),     intent(in) :: sphharT
    type(t_stars),      intent(in) :: starsT
    type(t_atoms),      intent(in) :: atomsT
    type(t_input),      intent(in) :: inputT
    complex,            intent(in) :: vGradCoul0IR(:, :)
    real,               intent(in) :: vrFleur(:, 0:, :, :)
    integer,            intent(in) :: ngdp
    integer,            intent(in) :: gdp(:, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    integer,            intent(in) :: nmem_atom(0:, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    real,               intent(in) :: rho0MT(:, 0:, :, :)
    complex,            intent(in) :: vpwStar(starsT%n3d, 1)
    real,               intent(in) :: paPoXCo
    real,               intent(in) :: paPoYCo
    real,               intent(in) :: paPoZCo
    logical,            intent(in) :: harSw
    logical,            intent(in) :: extSw
    logical,            intent(in) :: xcSw

    integer                        :: iGvec
    real                           :: Gext(3)
    complex                        :: vCoulBenchmark(ngdp, 3)
    complex                        :: gradVcMT(atomsT%jmtd, (atomsT%lmaxd + 1)**2, atomsT%nat, 3)
    real                           :: sqr4pi3
    real                           :: Vc0nSym(atomsT%jmtd, (atomsT%lmaxd + 1 +  1)**2, atomsT%nat) ! todo is this the correct dimension
    integer                        :: iatom
    integer                        :: itype
    integer                        :: ineq
    integer                        :: symType
    integer                        :: lh
    integer                        :: oqn_l
    integer                        :: lm_temp
    integer                        :: mem
    integer                        :: mems
    integer                        :: mqn_m
    integer                        :: mqn_mpp
    integer                        :: lm
    integer                        :: lmMinus
    integer                        :: lmPlus
    real                           :: Vc0MTDerived(atomsT%jmtd)
    real                           :: tempGaunt1, tempGaunt2
    integer                        :: imesh
    complex                        :: vpw_G(ngdp)
    complex                        :: vpw_G_coul(ngdp)
    complex                        :: vpw_G_hart(ngdp)
    complex                        :: vpw_G_xc(ngdp)
    complex              :: vpwCoulUwStars1(starsT%n3d, 1)
    complex              :: vpwCoulUwStars2(starsT%n3d, 1)
    complex              :: vpwStarDiff(starsT%n3d, 1)
    complex                        :: vpwRfleur(3)
    complex                        :: vpwRjuPhon(3)
    complex                        :: exponentialjuPhon
    complex                        :: exponentialFLEUR
    complex                        :: exponential
    integer           :: isIR
    integer           :: wasIR
    integer          :: c1, c2, c3
    real                           :: ucpath(3)
    real                              :: ucpathA(3), ucpathB(3), ucpathC(3)
    real                           :: ucpathExt(3)
    real                      :: ucpathExtc(3), ucpathExta(3)
    real                           :: dx
    integer                        :: ii
    real                           :: x
    complex                        :: newvpwFLEUR(starsT%n3d, 1) ! make this declaration more general!
    complex                        :: newvpw_G(ngdp)
    complex                        :: onlyPotFLEUR
    complex                        :: onlypsqFLEUR
    integer                        :: ieq
    integer                        :: ext2brav(3)
    integer                        :: ext2bravLock(3)
    integer                        :: atomLock
    integer                        :: itypeLock
    logical                        :: pseudoPot = .false.
    logical                        :: realPot = .true.
    integer                   :: iterations
  real   :: vrFleurDiff(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurxc(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrTemp(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  complex,             intent(in) :: gradVrjuPhon(:, :, :, :)
  complex,      allocatable :: gradVrFleur(:, :, :, :)
  real    :: direc(3)
  real :: direcExt(3)
  complex                                     :: ylm((atomsT%lmaxd + 1)**2)
  complex                                     :: vMTfleur(3, atomsT%jmtd)
  complex                                     :: vMTjuPhon(3, atomsT%jmtd)
  integer                          :: idirec
  real                              :: av, dmx, rms
  real                              :: IRvalue(1), MTvalue(1)
  integer                             :: jspin
  integer                               :: ilh
  real        :: potFleurMT(atomsT%jmtd)
  integer    :: symAt_temp, symOpr_temp
  real      :: rvec(3), rvec_int(3)
  integer    :: imatmulr, imatmulc
  real :: rotatedVector(3)
  integer :: ilath
  real :: linCombBas
  integer :: l_temp
  integer :: imem

logical :: fBenchSw

    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw
    fBenchSw = .true.
    !vpwStar = 0
    vpwStarDiff = 0
    !vrFLEUR = 0
    vrFLEURDiff = 0

!    !todo do this with pottot and potcoul
!    if ( harSw .and. .not.extSw .and. .not.xcSw ) then ! Hartree potential
!      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else if ( .not.harSw .and. extSw .and. .not.xcSw ) then ! external potential
!      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStarDiff
!      call fclose( 1000 )
!
!      vpwStar = vpwStar - vpwStarDiff
!
!
!      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
!        read(1000) vrFLEURDiff
!      call fclose(1000)
!
!      vrFLEUR = vrFLEUR - vrFLEURDiff
!
!    else if ( .not.harSw .and. .not.extSw .and. xcSw ) then ! exchange-correlation potential
!      call fopen( 1000, name='vpw_xc', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else if ( harSw .and. extSw .and. .not.xcSw ) then ! Coulomb potential
!      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else if ( harSw .and. extSw .and. xcSw ) then ! effective potential
!      call fopen( 1000, name='vpw_eff', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!      call fopen(1000, name='v0MTFLEUR_eff', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else
!      call juDFT_warn('FLEUR benchmark for this potential combination not available. Setting benchmark curve equals zero.', &
!        & calledby='checkGradPot0s', hint='Choose Hartree, external, exchange-correlation, Coulomb or effective potential', &
!        & file='jpTestPotential_mod.f90')
!      fBenchSw = .false.
!
!    end if

    call ConvertStar2G( vpwStar(:, 1), vpw_G, starsT, ngdp, gdp )

    ! In this case we compare with minus the gradient
    do iGvec = 1, ngdp
      Gext(:) = matmul(cellT%bmat, gdp(:, iGvec))
      vCoulBenchmark(iGvec, :) = -iu * Gext(:) * vpw_G(iGvec)
    enddo
    call calcGrFinLH(atomsT, sphharT, clnu_atom, nmem_atom, mlh_atom, vrFleur(:, :, :, 1), gradVrFleur)


    !call fopen(1111, name='pathPseudoPot', status='replace', action='write', form='formatted')
    !call fopen(1222, name='potFLEUR', status='replace', action='write', form='formatted')
    call fopen(1333, name='pathGrDens', status='replace', action='write', form='formatted')
    call fopen(1444, name='pathGrDensCont', status='replace', action='write', form='formatted')

    dx = 1. / 600.
    x = 0 - dx
    wasIR = -1
    do ii = 0, 600
      x = x + dx
      ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
      iatom = 1
      isIR = 0
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          do c1 = -1, 1
            do c2 = -1, 1
              do c3 = -1, 1
                ext2brav = [c1, c2, c3]
                ucpathExt= matmul(cellT%amat, ucpath - atomsT%taual(:, iatom) - ext2brav)
                if ( norm2 ( ucpathExt ) < atomsT%rmt(itype) ) then
                  isIR = iatom ! todo atomLock is redundant!
                  atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                  ext2bravLock = ext2brav
                  itypeLock = itype
                end if
              end do
            end do
          end do
          iatom = iatom + 1
        end do
      end do

      if ( pseudoPot ) then

        vpwRfleur = cmplx(0.,0.)
        vpwRjuPhon = cmplx(0.,0.)

        do iGvec = 1, ngdp
          exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
          vpwRfleur(:)  = vpwRfleur(:)  + vCoulBenchmark(iGvec, :) * exponential
          vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :)   * exponential
        end do

      !  write (1111, '(13(es16.8E3, 2x),i2, 3(2x, i2))') x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), &
      !    & aimag(vpwRfleur(2)), real(vpwRfleur(3)), aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)),     &
      !    & real(vpwRjuPhon(2)), aimag(vpwRjuPhon(2)), real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), &
      !    & ext2bravLock(2), ext2bravLock(3)

      end if

      if ( realPot ) then ! if no pseudoPotMode
        if ( wasIR == -1 ) then ! initial point of path
          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
              !todo floating is inexact here for some reasons....
              onlyPotFLEUR = onlyPotFLEUR +1!+ !vpw_G(iGvec) * exponentialFLEUR ! pure potential
            end do

            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
        else ! wasIR is available
          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do

            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          else if (wasIR == 0 .and. isIR >= 0 ) then ! crosses MT surface into the atom

            ucpathA = ucpath ! ucpathA is in atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
            do iterations = 1, 1000
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              !  ! mtsphere found
                exit
              else
                if ( ( norm2(ucpathExtc) - atomsT%rmt( itypeLock ) ) < 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations > 1000) then
              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
            end if

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0., 0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpathC))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), wasIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock 
            direcExt = matmul(cellT%amat, direc)
            direcExt = direcExt / norm2(direcExt)
       !     call ylmnorm_init(atomsT%lmaxd + 1)
            !call ylm4(atomsT%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
       !     call ylmnorm_init(atomsT%lmaxd)

            vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atomsT%lmax(itypeLock)! + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atomsT%jri(itypeLock)
                    vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) - gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, atomLock, idirec) * ylm(lm)
                  end do
                end do
              end do
            end do

            do imesh = atomsT%jri(itypeLock), 1, -1
              write(1333, '(13(es20.10E3, 2x), i2, 3(2x, i2))')&
                &-atomsT%rmsh(imesh, itypeLock), real(vMTfleur(1, imesh)), aimag(vMTfleur(1, imesh)), &
                &real(vMTfleur(2, imesh)), aimag(vMTfleur(2, imesh)), real(vMTfleur(3, imesh)), &
                &aimag(vMTfleur(3, imesh)), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            do idirec = 1, 3
              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              else
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          else if (wasIR >= 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

            ucpathA = ucpath ! is out of atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
            do iterations = 1, 1000 ! really so much iterations needed?
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
                ! mtsphere found
                exit
              else
                if ( ( norm2( ucpathExtc ) - atomsT%rmt( itypeLock ) ) > 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations >= 1000) then
              write(*, *) 'Warning: Iteration not converged'
            end if

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            direcExt = direcExt / norm2(direcExt)
      !      call ylmnorm_init(atomsT%lmaxd + 1)
            !call ylm4(atomsT%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
      !      call ylmnorm_init(atomsT%lmaxd)

            vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atomsT%lmax(itypeLock) !+ 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atomsT%jri(itypeLock)
                    vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) - gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm,  atomLock, idirec) * ylm(lm)
                  end do
                end do
              end do
            end do
            do imesh = 1, atomsT%jri(itypeLock)
              write(1333, '(13(es20.10E3, 2x),i2, 3(2x, i2))')& 
                &atomsT%rmsh(imesh, itypeLock), real(vMTfleur(1, imesh)), aimag(vMTfleur(1, imesh)), &
                &real(vMTfleur(2, imesh)), aimag(vMTfleur(2, imesh)), real(vMTfleur(3, imesh)), &
                &aimag(vMTfleur(3, imesh)), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                &aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            !todo seems to be the real ptoential here
            !iatom = 0
            potFleurMT = 0
              !Evaluate MT sphere function at same random points
            !symAt_temp  = atomsT%ntypsy(atomLock) ! symmetry of atom iatom
            !symOpr_temp = atomsT%ngopr (atomLock) ! symmetry operation mapping local to global coordinate system
            !do imesh = 1, atomsT%jri(itypeLock)
            !  rvec = direcExt
            !  if (symOpr_temp /= 1) then
            !  !Transform into internal real space coordinates
            !  !todo might be wrong bmat
            !    call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
            !    do imatmulr = 1, 3
            !      rotatedVector(imatmulr) = 0.0
            !      do imatmulc = 1, 3
            !        rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) * rvec_int(imatmulc)
            !      end do
            !    end do
            !    call cotra0(rotatedVector, rvec, cellT%amat)
            !  end if
            !  call ylm4(atomsT%lmax(itypeLock), rvec, ylm)

            !  do ilath = 0, sphharT%nlh(symAt_temp)
            !    linCombBas = 0.0
            !    l_temp = sphharT%llh(ilath, symAt_temp) * (sphharT%llh(ilath, symAt_temp) + 1) + 1 !understand this
            !    do imem = 1, sphharT%nmem(ilath, symAt_temp)
            !      lm_temp = l_temp + sphharT%mlh(imem, ilath, symAt_temp)
            !      linCombBas = linCombBas + real(clnu_atom(imem, ilath, symAt_temp) * ylm(lm_temp))
            !    end do
            !      potFleurMT(imesh) = potFleurMT(imesh) + vrFleur(imesh, ilath, itypeLock, 1) * linCombBas
            !  end do
            !end do
            !call fopen(1020, name='potMTfleur', status='replace', action='write', form='formatted')
            !do imesh = 1, atomsT%jri(itypeLock)
            !  write(1020, '(3(es16.8E3, 2x))') atomsT%rmsh(imesh, itypeLock), potFleurMT(imesh), real(vMTfLEUR(1, imesh))
            !end do
            !call fclose(1020)


            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0., 0.)
            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpathC))
              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR) ! todo the lines down really needed?
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            do idirec = 1, 3
              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              else
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          end if! if still in MT then wait for crossing

        end if !is not first point on path
      end if ! is not in pseudopot mode
      wasIR = isIR
    end do
    !call fclose(1111)
    !call fclose(1222)
    call fclose(1333)
    call fclose(1444)
  end subroutine checkGradDens0s
  !todo unite with checkDensnPot from jpTestInput
  ! this only checks the continuity of the potential and writes it to the log file
  ! This routine also supports lmax functions instead of lmax + 1 and complex functions
subroutine checkjuPhDens1(atomsT, cellT, ngdp, fIR, fMufT, qpt, gdp, noPts2chCont, logFileUnit)

    use m_types
    use m_cotra
    use m_ylm_old
    use mod_juPhonUtils
    use m_jPConstants
    implicit none
    type(t_atoms),        intent(in) :: atomsT
    type(t_cell),         intent(in) :: cellT
    integer,              intent(in) :: noPts2chCont !number of points to check continuity
    integer,              intent(in) :: ngdp
    integer,              intent(in) :: gdp(:, :)
    integer,              intent(in)            :: logFileUnit ! why not as intent(in)?
    complex,             intent(in) :: fIR(:, :)
    complex,             intent(in) :: fMufT(:, :, :, :)
    real, intent(in) :: qpt(:)
!    real,                 intent(in)  :: rmt(:) !MT-radii of atom type itype
!    integer,              intent(in)  :: ntype ! number of atom types
!    integer,              intent(in)  :: nop !number of symmetry operations
!    integer,              intent(in)  :: nq3 ! number of stars
!    logical,              intent(in)  :: symor !are symmetry operations symorphic or noct
!    integer,              intent(in)  :: kv3(:, :) !reciprocal g-vector of star
!    integer,              intent(in)  :: mrot(:, :, :) !symmetry operations, rotation matrices
!    real,              intent(in)  :: tau(:, :) !symmetry operations, translation vectors
!    integer,              intent(in)  :: invtab(:) !lists which symmetry operations are inverse to symmetry operations with indesx isym
!    real,                 intent(in)  :: bmat(:, :) !reciprocal Bravais matrix transposed
!    real,                 intent(in)  :: amat(:, :) !Bravais matrix
!    real,                 intent(in)  :: pos(:, :) !Postions of atoms in Brillouin zone
!    integer,              intent(in)  :: nmem(:, :) !number of lattice harmonics
!    integer,              intent(in)  :: llh(:, :) !ls of  lattice harmonics
!    complex,              intent(in)  :: clnu(:, :, :) !phasefactors in linear combination of lattice harmonics
!    logical,              intent(in)  :: density ! is function density?
!    complex,              intent(in)  :: fpw(:, :)
!    real,                 intent(in)  :: fr(:, :, :, :)
!    integer,              intent(in)  :: nstr(:)
!    integer,              intent(in)  :: nat
!    integer,              intent(in)  :: ntypsy(:)
!    integer,              intent(in)  :: ngopr(:)
!    integer,              intent(in)  :: lmax(:)
!    integer,              intent(in)  :: nlh(:)
!    integer,              intent(in)  :: mlh(:, :, :)
!    integer,              intent(in)  :: jri(:)
!    integer,              intent(in)  :: lmaxd
!    integer,              intent(in)  :: neq(:)
!
!
    integer                           :: oqn_l
    integer                           :: mqn_m
    integer                           :: lm_lcontr
    integer                           :: lm
    integer                           :: irandPt
    integer                           :: ieq !loop variable
    real                              :: euclidNorm !variable to store eculidian L2 norm of vector
    real,   allocatable               :: randPtsCart(:, :, :) !cartesian dimension, number of Pts, itype, random points
    real,   allocatable               :: randPtsGrid(:, :, :) !random points on grid, coordinates, number of points
    real,   allocatable               :: randPtsMTLoc(:, :, :) !random points in local atomic coordinate (coordinates, number of points, ntype)
    integer                           :: itype !certain atom type
!    complex                           :: starAtRanPt(nq3) !value of star at certain point ! ng3 = nq3
    complex,   allocatable               :: sumfuncvalI(:, :, :) !sum of function values coming from interstitial
    integer                         :: iGvec
    integer                          :: idirec
!    real                              :: scaling !scaling factor which makes difference between density and potential or wavefunction
!    integer                           :: icoord !loop variable
    integer                           :: imatmulc, imatmulr !loop variable
    real                              :: rvec(3) !randoom point
    real                              :: rvec_int(3)
    real  :: rotatedVector(3) !rotated vector after symmetry operation
!    integer                           :: istar !loop variable
    integer                           :: iatom
!    integer                           :: symAt_temp
    integer                           :: symOpr_temp
    complex                           :: ylm((atomsT%lmaxd + 1)**2) !TODO eigentlich lmax(itype)
    complex,    allocatable              :: sumOfMTFuncVal(:, :,:)
!    integer                           :: ilath
!    real                              :: linCombBas
!    integer                           :: l_temp
!    integer                           :: imem
!    integer                           :: lm_temp
    real                              :: av(3), dms(3), rms(3)
!    integer                           :: ispin

    allocate(randPtsCart(3, noPts2chCont, atomsT%nat), randPtsGrid(3, noPts2chCont, atomsT%ntype), randPtsMTLoc(3, noPts2chCont, atomsT%ntype) )

    allocate(sumfuncvalI(3, noPts2chCont,atomsT%nat), sumOfMTFuncVal(3, noPts2chCont,atomsT%nat))
   !Generate noPts2chCont random points on the MT sphere
    randPtsCart(:, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        call sphpts(randPtsCart(:, :, iatom), noPts2chCont, atomsT%rmt(itype), atomsT%pos(:, iatom))
!      do irandPt = 1, noPts2chCont
!        call random_number(randPtsCart(:, irandPt, itype))
!        euclidNorm = norm2(randPtsCart(:, irandPt, itype))
!        randPtsCart(:, irandPt, itype) = randPtsCart(:, irandPt, itype) / euclidNorm * atomsT%rmt(itype)
      end do
    end do
    

    !Evaluate Stars for this random points at each MT sphere boundary in the unit cell
    sumfuncValI = cmplx(0.0) !todo this is probably unneccessary, because it is already set to 0 in the loop
    iatom       = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        do irandPt = 1, noPts2chCont
!          rvec = atomsT%pos(:,iatom) + randPtsCart(:,irandPt,itype)
          rvec(:) = randPtsCart(:, irandPt, iatom)
          ! transform rvec in internal, real-space coordinates --> rvec_int
          call cotra1(rvec, rvec_int, cellT%bmat) !todo do this inline from now on!!!
          ! Evaluate stars at random point
          sumfuncValI(:, irandPt,iatom) = cmplx(0., 0.) !todo somehow this does not make sense
          do iGvec = 1, ngdp !nq3 is number of stars
            do idirec = 1, 3
              sumfuncValI(idirec, irandPt,iatom) = sumfuncvalI(idirec, irandPt,iatom) + fIR(iGvec, idirec) * exp(iu * tpi * dot_product(gdp(:, iGvec) + qpt(:), rvec_int(:))) !jsp noch einfügen, warum Multiplikation with G-vectors
            end do ! idirec
          end do  !iGvec
        end do  !irandPt
      end do  !ieq
    end do  !itype

    if (all(qpt(:) < 1e-6)) then
      if (any(aimag(sumfuncValI(:, :, :)) > 1e-8)) then
        write(*, *) "Imaginary part of observable potential should be 0"
        !NOstopNO
      end if
    end if


    iatom = 0
    sumOfMTFuncVal = cmplx(0., 0.)
    do itype = 1, atomsT%ntype
      !Evaluate MT sphere function at same random points
      randPtsGrid = 0.0
      do ieq = 1, atomsT%neq(itype)
        iatom       = iatom + 1
        !symOpr_temp = 1!atomsT%ngopr(iatom) ! symmetry operation mapping local to global coordinate system
        do irandPt = 1, noPts2chCont
          !rvec = randPtsCart(:,irandPt,itype)
          ! We have to subtract the atom position to get into a local coordinate system of atom iatom
          rvec = randPtsCart(:,irandPt,iatom) - atomsT%pos(:, iatom)
         ! if (symOpr_temp /= 1) then
         !   !Transform into internal real space coordinates
         !   call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
         !   do imatmulr = 1, 3 !todo optimize this
         !     rotatedVector(imatmulr) = 0.0
         !     do imatmulc = 1, 3
         !       rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) *rvec_int(imatmulc)
         !     end do
         !   end do
         !   call cotra0(rotatedVector, rvec, cellT%amat)
         ! end if
!          call ylmnorm_init(atomsT%lmaxd + 1)
          !call ylm4(atomsT%lmax(itype) + 1, rvec, ylm)
          call ylm4(atomsT%lmax(itype) , rvec, ylm)
!          call ylmnorm_init(atomsT%lmaxd)

          do oqn_l = 0, atomsT%lmax(itype)! + 1
            lm_lcontr = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lcontr + mqn_m
              do idirec = 1, 3
                sumOfMTFuncVal(idirec, irandPt, iatom) = sumOfMTFuncVal(idirec, irandPt, iatom) + fMufT(atomsT%jri(itype), lm, iatom, idirec) * ylm(lm)
              end do ! idir
            end do
          end do
        end do
      end do
    end do

    if (all(qpt(:) < 1e-6)) then
      if (any(aimag(sumOfMTFuncVal(:, :, :)) > 1e-8)) then
        write(*, *) "Imaginary part of observable potential should be 0"
        !NOstopNO
      end if
    end if

    write ( logFileUnit, '(a)' ) "Continuity of the unperturbed potential's gradient:"
    write(logFileUnit,*)         '---------------------------------------------------'
    iatom = 0
    do itype = 1,atomsT%ntype
      do ieq = 1,atomsT%neq(itype)
        iatom = iatom + 1
        write(logFileUnit,*) '  Atom:',iatom
        write(logFileUnit,'(a)') '  coordinate(cartesian)      IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        do irandPt = 1,noPts2chCont
            write(logFileUnit,'(3f8.4,3x,3(4f12.8, 3x))') randPtsCart(:,irandPt,iatom),sumfuncValI(1, irandPt,iatom),sumOfMTFuncVal(1, irandPt,iatom), sumfuncValI(2, irandPt,iatom),sumOfMTFuncVal(2, irandPt,iatom), sumfuncValI(3, irandPt,iatom),sumOfMTFuncVal(3, irandPt,iatom)
        end do
        do idirec = 1, 3
          call fitchkNeg(real(sumfuncvalI(idirec, :,iatom)), real(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
        end do
        write(logFileUnit, *)
      !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'real average absolute value interstitial :',av(1), av(2), av(3)
        write(logFileUnit,'(2x,A,6f15.8,A/)')               'real rms(x, y, z), real dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
        av(:) = 0.
        rms(:) = 0.
        dms(:) = 0.
        do idirec = 1, 3
          call fitchkNeg(aimag(sumfuncvalI(idirec, :,iatom)), aimag(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
        end do
        write(logFileUnit, *)
      !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
        write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'aimag average absolute value interstitial :',av(1), av(2), av(3)
        write(logFileUnit,'(2x,A,6f15.8,A/)')               'aimag rms(x, y, z), dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
      end do
    end do
  end subroutine checkjuPhDens1

  ! This is a fork of fitchk to handle quantities where the interstitial varies a lot so the average is almost zero
  subroutine fitchkNeg(f1, f2, n, av, rms, dmx)

    implicit none

    integer, intent(in) :: n

    real,    intent(in) :: f1(n)
    real,    intent(in) :: f2(n)

    real,    intent(out) :: av
    real,    intent(out) :: dmx
    real,    intent(out) :: rms


    real                 :: d
    integer              :: i

    av = 0.
    rms = 0.
    dmx = 0.
    do i = 1,n
       av = av + abs(f1(i))
       d = (f1(i)-f2(i))**2
       dmx = max(d,dmx)
       rms = rms + d
    end do
    av = av/n
    if (abs(av).LT.1.e-30) then
       rms = 0.
       dmx = 0.
       return
    end if
    rms = sqrt(rms/n)/av*100.
    dmx = sqrt(dmx)/av*100.

  end subroutine fitchkNeg

  ! this is for testGoldstein Modes and irrelevant now
  subroutine benchMarkMTGrad( atoms, input, lathar, sym, cell, clnu_atom, nmem_atom, mlh_atom  )

    use m_types
    use mod_juPhonUtils, only : calcGrFinLH, calcGrR2FinLH
    use m_jpConstants, only : fpi

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell

    ! Array parameters
    complex,                        intent(in)  :: clnu_atom(:, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:, 0:, :)

    ! Scalar variables
    integer                                     :: iatom
    integer                                     :: itype
    integer                                     :: ieqat
    integer                                     :: idir
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: imesh

    ! Array variables
    real,              allocatable              :: testFunction(:,:,:,:)
    complex,           allocatable              :: grTestFunction(:, :, :, :)
    complex,           allocatable              :: r2GrTestFunction(:, :, :, :)

    allocate( testFunction(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    testFunction(:, :, :, :) = cmplx(0., 0.)


    ! Testfunction f(r) = r
    testFunction(:, :, :, :) = cmplx(0., 0.)

    do imesh = 1, atoms%jri(1)
      testFunction(imesh, 0, 1, 1) = atoms%rmsh(imesh, 1)**3 * sqrt(fpi)
    end do

    call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1),  r2GrTestFunction)

    do imesh = 1, atoms%jri(1)
      testFunction(imesh, 0, 1, 1) = testFunction(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2
    end do

    call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1), grTestFunction)

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2GrTestFunction(imesh, lm, iatom, idir) = r2GrTestFunction(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(5794, '(4(i8),2x,3(f24.9))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype), grTestFunction(imesh, lm, iatom, idir)
                write(5795, '(4(i8),2x,3(f24.9))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype), r2GrTestFunction(imesh, lm, iatom, idir)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! Testfunction f(r) = 1 / r
    testFunction(:, :, :, :) = cmplx(0., 0.)
    do imesh = 1, atoms%jri(1)
      testFunction(imesh, 0, 1, 1) = atoms%rmsh(imesh, 1) * sqrt(fpi)
    end do

    call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1),  r2GrTestFunction)

    do imesh = 1, atoms%jri(1)
     testFunction(imesh, 0, 1, 1) = testFunction(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2
    end do

    call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1), grTestFunction)

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2GrTestFunction(imesh, lm, iatom, idir) = r2GrTestFunction(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    if (.true.) then
      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  write(5797, '(4(i8),2x,3(f24.9))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype), grTestFunction(imesh, lm, iatom, idir)
                  write(5798, '(4(i8),2x,3(f24.9))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype), r2GrTestFunction(imesh, lm, iatom, idir)
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idir
    end if


    ! Testfunction f(r) = 1 / r^2
    testFunction(:, :, :, :) = cmplx(0., 0.)
    do imesh = 1, atoms%jri(1)
      testFunction(imesh, 0, 1, 1) = sqrt(fpi)
    end do
    deallocate(r2GrTestFunction, grTestFunction)

    call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1),  r2GrTestFunction)

    do imesh = 1, atoms%jri(1)
      testFunction(imesh, 0, 1, 1) = testFunction(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2
    end do

    call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1), grTestFunction)

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2GrTestFunction(imesh, lm, iatom, idir) = r2GrTestFunction(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(5802, '(4(i8),2x,3(f30.9))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype), grTestFunction(imesh, lm, iatom, idir)
                write(5803, '(4(i8),2x,3(f30.9))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype), r2GrTestFunction(imesh, lm, iatom, idir)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    write(*, *) 'benchmark r2ornotgradient finished'
  end subroutine benchMarkMTGrad

end module jpTest1stVarDens
