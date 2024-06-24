! todo Test whether with q and -q Eii2 gives the same results
! todo Plot all Eii2 and write the matrices to a file
module m_jpTestDynMatHF

  use m_types

  implicit none

  contains

  subroutine TestDynMatHF( atoms, lathar, input, stars, cell, dimens, sym, qpts, harSw, extSw, xcSw, testEii2PsDensSw, testEii2LatPeriodQSw, ngdp, paPoX, &
      & paPoY, paPoZ, memd_atom, gdp, mlh_atom, nmem_atom, clnu_atom, vExt2IR, vExt2MT, rho0IR, rho0MT, logUnit )

    use m_jp2ndOrdQuant, only : GenVext2
    use m_jpPlotObservables
    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_sphhar),                 intent(in) :: lathar
    type(t_input),                  intent(in) :: input
    type(t_stars),                  intent(in) :: stars
    type(t_cell),                   intent(in) :: cell
    type(t_dimension),              intent(in) :: dimens
    type(t_sym),                    intent(in) :: sym
    type(t_kpts),                   intent(in) :: qpts

    ! Scalar parameters
    integer, intent(in) :: logUnit
    logical,                        intent(in) :: harSw
    logical,                        intent(in) :: extSw
    logical,                        intent(in) :: xcSw
    integer,                        intent(in) :: ngdp
    real,                           intent(in) :: paPoX
    real,                           intent(in) :: paPoY
    real,                           intent(in) :: paPoZ
    integer,                        intent(in) :: memd_atom
    logical,                        intent(in) :: testEii2PsDensSw
    logical,                        intent(in) :: testEii2LatPeriodQSw

    ! Array parameters
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    real,                          intent(in)  :: rho0MT(:, 0:, :, :)
    complex,                       intent(in)  :: rho0IR(:,:)

    ! Array variables
    complex,           allocatable             :: vExt2IR(:, :, :, :)
    complex,           allocatable             :: vExt2MT(:, :, :, :)


    logical :: testgradGradTrlYlm
    logical :: testLegPolyVext2

    if ( .false. ) then
      call GenVext2(atoms, cell, dimens, ngdp, gdp, vExt2IR, vExt2MT, .true.)
      call TestPlotVext2( atoms, lathar, input, stars, cell, dimens, sym, harSw, extSw, xcSw, ngdp, paPoX, paPoY, paPoZ, memd_atom,&
        & gdp, mlh_atom, nmem_atom, clnu_atom, vExt2IR, vExt2MT/2., logUnit )
      write(*, *) 'bug: factor of 2 above?'
    end if

    if (.false.) then
      write(*, *) 'scan eii2q'
     call ScanEii2q( atoms, cell, stars, dimens, input, ngdp, gdp )
    end if

    if ( testEii2LatPeriodQSw .or. testEii2PsDensSw ) then
      write(*, *)
      write(*, '(a)') 'Initiating dynamical matrix Hellmann--Feynman contribution test(s)...'
      write(*, '(a)') '---------------------------------------------------------------------'
      write(logUnit, *)
      write(logUnit, '(a)') 'Dynamical matrix Hellmann--Feynman contribution test(s)'
      write(logUnit, '(a)') '*******************************************************'
      write(logUnit, *)
    else
      write(*, '(a)') '----------------------------------------------------------------'
      write(*, '(a)') 'DISABLED dynamical matrix Hellmann-Feynman contribution test(s)!'
      write(*, '(a)') '----------------------------------------------------------------'
      write(logUnit, '(a)')
      write(logUnit, '(a)') 'DISABLED dynamical matrix Hellmann-Feynman contribution test(s)!'
      write(logUnit, '(a)') '****************************************************************'
      return
    end if

    if (testEii2LatPeriodQSw) then
      write(*, '(2x,a)') 'Performing TestEii2LatPeriodQ...'
      call testEii2LatPeriodQ( atoms, stars, input, cell, dimens, qpts, gdp, ngdp, logUnit )
    else
      write(*, '(2x,a)') 'DISABLED TestEii2LatPeriodQ!'
      write (logUnit, '(a)')'Check Eii2 lattice periodicity with resepct to q!'
      write (logUnit, '(a)')'-------------------------------------------------'
      write (logUnit, '(a)')'                                                |__ DISABLED!'
    end if

    if (testEii2PsDensSw) then
      write(*, '(2x,a)') 'Performing TestEii2PsDens...'
      call testEii2PsDens(atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, logUnit, memd_atom, rho0IR, rho0MT, gdp, mlh_atom, nmem_atom, clnu_atom)
    else
      write(*, '(2x,a)') 'DISABLED TestEii2PsDens!'
      write (logUnit, '(a)')'Check Eii2 lattice periodicity with resepct to q!'
      write (logUnit, '(a)')'-------------------------------------------------'
      write (logUnit, '(a)')'                                                |__ DISABLED!'
    end if

  end subroutine testDynMatHF


  subroutine testEii2LatPeriodQ( atoms, stars, input, cell, dimens, qpts, gdp, ngdp, logUnit)

    use m_jpPotDensHelper, only : genPertPotDensGvecs
    use m_jp2ndOrdQuant, only : CalcIIEnerg2MatElem, GenPsDens2ndOrd, CalcIIEnerg2
    use m_jPConstants, only : fpi
    use m_juDFT_NOstopNO, only : juDFT_warn
    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in) :: atoms
    type(t_stars),                 intent(in) :: stars
    type(t_input),                 intent(in) :: input
    type(t_cell),                  intent(in) :: cell
    type(t_dimension),             intent(in) :: dimens
    type(t_kpts),                  intent(in) :: qpts

    ! Scalar parameters
    integer,                       intent(in) :: ngdp
    integer,                       intent(in) :: logUnit

    ! Array parameters
    integer,                       intent(in) :: gdp(:, :)

    ! Scalar variables
    integer                                   :: ngpqdp
    integer                                   :: ngpqdp2km
    integer                                   :: iAtype
    integer                                   :: iBtype
    integer                                   :: iAatom
    integer                                   :: iBatom
    integer                                   :: iAdir
    integer                                   :: iBdir
    integer                                   :: iG
    integer                                   :: iGvar
    integer                                   :: G0index
    logical                                   :: testMode
    logical                                   :: passed = .true.
    real                                      :: minQcomp
    integer                                   :: iqpt
    integer                                   :: idir
    integer                                   :: lcm
    integer                                   :: subIqpt

    ! Array variables
    real                                      :: qpoint(3)
    real                                      :: Gpqext(3)
    integer                                   :: gpqdp2iLim(2, 3)
    integer                                   :: gShift(3)
    integer,           allocatable            :: qIndex(:, :, :)
    integer,           allocatable            :: gpqdp(:, :)
    integer,           allocatable            :: gpqdp2(:, :)
    integer,           allocatable            :: gpqdp2Ind(:, :, :)
    complex,           allocatable            :: E2ndOrdIIat(:, :)
    complex,           allocatable            :: E2ndOrdIIat2(:, :)
    complex,           allocatable            :: psDens2ndOrd(:, :, :, :)
    integer,           allocatable            :: gdpIndex(:, :, :)

    write (logUnit, '(a)')'Check Eii2 lattice periodicity with resepct to q!'
    write (logUnit, '(a)')'-------------------------------------------------'
    allocate( E2ndOrdIIat(3, 3) )
    allocate( E2ndOrdIIat2(3, 3) )
    E2ndOrdIIat = cmplx(0., 0.)
    E2ndOrdIIat2 = cmplx(0., 0.)
    iAtype = 1
    iBtype = 1
    iAatom = 1
    iBatom = 1
    testMode = .false.

    minQcomp = 2.
    do iqpt = 1, qpts%nkpt
      do idir = 1, 3
        if ( (qpts%bk(idir, iqpt) > 1e-12) .and. (qpts%bk(idir, iqpt) < minQcomp)) minQcomp = qpts%bk(idir, iqpt)
      end do ! idir
    end do ! iqpt
    lcm = int(1. / minQcomp)
    allocate(qIndex(0:lcm, 0:lcm, 0:lcm))
    qIndex = 0
    do iqpt = 1, qpts%nkpt
      qIndex(nint(qpts%bk(1, iqpt) * lcm), nint(qpts%bk(2, iqpt) * lcm), nint(qpts%bk(3, iqpt) * lcm)) = iqpt
    end do ! iqpt

    ! Test q = (0 0 0.25)
!    qpoint = [0., 0., 0.25]
!    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
!    call GenPsDens2ndOrd(atoms, cell, dimens, ngpqdp, G0index, gpqdp, qpoint, psDens2ndOrd, testMode)
!!    write(*, *) 'ng', ngpqdp
!!    write(*, '(a)') 'Eii2 for 0.25'
!!    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(1, :)
!!    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(2, :)
!!    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(3, :)
!
!!    do iBdir = 1, 3
!!      do iAdir = 1, 3
!!        do iG = 1, ngpqdp
!!          write(2213, '(6i8,2f15.8)') iBdir, iAdir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), psDens2ndOrd(iG, iAdir, 1, iBdir)
!!        end do ! iG
!!      end do ! iAdir
!!    end do ! iBdir
!
!    Gpqext(:) = 0.
!    E2ndOrdIIat(:, :) = cmplx(0., 0.)
!    do iBdir = 1, 3
!      do iAdir = 1, 3
!        do iG = 1, ngpqdp
!          Gpqext(1:3) = matmul(cell%bmat, gpqdp(1:3, iG) + qpoint(1:3))
!          !todo actually we also need an exp(i (G + q) tau)
!          E2ndOrdIIat(iAdir, iBdir) = E2ndOrdIIat(iAdir, iBdir) + atoms%zatom(1) * fpi / norm2(Gpqext(1:3))**2 &
!                                                                            & * psDens2ndOrd(iG, iAdir, 1, iBdir)
!        end do ! iG
!      end do ! iAdir
!    end do ! iBdir
!!    write(*, '(a)') 'Eii2 for 0.25'
!!    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(1, :)
!!    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(2, :)
!!    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(3, :)
!
!    ! For 0.25
!    !call CalcIIEnerg2MatElem( atoms, cell, dimens, qpoint, ngpqdp, gpqdp, iAtype, iBtype, iAatom, iBatom, E2ndOrdIIat )

    ! todo this does not work for polyatomic systems
    write(*, *) 'todo this does not work for polyatomic systems'
    do iqpt = 1, qpts%nkpt
      do idir = 1, 3
        if (qpts%bk(idir, iqpt) > 1e-12) then
          gShift(idir) = 1
        else
          gshift(idir) = 0
        end if
      end do ! idir

        subIqpt = iqpt
        qpoint(1:3) = qpts%bk(1:3, subIqpt)
      E2ndOrdIIat(:, :) = 0.
!      call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
      call CalcIIEnerg2(atoms, cell, dimens, qpts, stars, input, subIqpt, ngdp, gdp, E2ndOrdIIat)
!      deallocate(gpqdp2Ind, gpqdp)
      if (.false.) then
        write(*, '(a)') 'Eii2ordIIat2'
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(1, :)
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(2, :)
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(3, :)
      end if


     ! deallocate(gpqdp2Ind, psDens2ndOrd)
  !    ngpqdp = 0
  !
  !    ! Test q = (0 0 -0.25)
  !    qpoint = [0., 0., -0.25]
  !    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp2, gpqdp2Ind, gpqdp2iLim )
  !    allocate(gdpIndex(minval(gpqdp2(1, :)) : maxval(gpqdp2(1, :)), minval(gpqdp2(2, :)) : maxval(gpqdp2(2, :)), minval(gpqdp2(3, :))-1 : maxval(gpqdp2(3, :))))
  !    gdpIndex(:, :, :) = 0
  !    do iG = 1, ngpqdp
  !      gdpIndex(gpqdp2(1, iG), gpqdp2(2, iG), gpqdp2(3, iG)) = iG
  !    end do ! iG
  !    call GenPsDens2ndOrd(atoms, cell, dimens, ngpqdp, G0index, gpqdp2, qpoint, psDens2ndOrd, testMode)
  !
  !    do iBdir = 1, 3
  !      do iAdir = 1, 3
  !        do iG = 1, ngpqdp
  !          !iGvar = gdpIndex(-gpqdp(1, iG), -gpqdp(2, iG), -gpqdp(3, iG))
  !          ! todo segfault because the pointer array is not quite correct, but this was only the way to the final test, I guess, so this is not so severe at the moment because calcIIEnerg2 does work
  !!          write(2214, '(6i8,2f15.8)') iBdir, iAdir, iG, -gpqdp2(1, iGvar), -gpqdp2(2, iGvar), -gpqdp2(3, iGvar), psDens2ndOrd(iGvar, iAdir, 1, iBdir)
  !        end do ! iG
  !      end do ! iAdir
  !    end do ! iBdir
  !    call CalcIIEnerg2MatElem( atoms, cell, dimens, qpoint, ngpqdp, gpqdp, iAtype, iBtype, iAatom, iBatom, E2ndOrdIIat )
  !    write(*, *) 'ng', ngpqdp
  !    write(*, '(a)') 'Eii2 for 0.25'
  !    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(1, :)
  !    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(2, :)
  !    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(3, :)

  !    deallocate(gpqdp2Ind, gpqdp2, gdpIndex)

!      ! Test q = (0 0 -0.25)
!      qpoint = [0., 0., 0.75]
!      call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp2, gpqdp2Ind, gpqdp2iLim )
!      allocate(gdpIndex(minval(gpqdp2(1, :)) : maxval(gpqdp2(1, :)), minval(gpqdp2(2, :)) : maxval(gpqdp2(2, :)), minval(gpqdp2(3, :))-1 : maxval(gpqdp2(3, :))))
!      gdpIndex(:, :, :) = 0
!      do iG = 1, ngpqdp
!        gdpIndex(gpqdp2(1, iG), gpqdp2(2, iG), gpqdp2(3, iG)) = iG
!      end do ! iG
!      call GenPsDens2ndOrd(atoms, cell, dimens, ngpqdp, G0index, gpqdp2, qpoint, psDens2ndOrd, testMode)
!
!      do iBdir = 1, 3
!        do iAdir = 1, 3
!          do iG = 1, ngpqdp
!            ! todo segfault because the pointer array is not quite correct, but this was only the way to the final test, I guess, so this is not so severe at the moment because calcIIEnerg2 does work
!  !          iGvar = gdpIndex(-gpqdp(1, iG), -gpqdp(2, iG), -gpqdp(3, iG)-1)
!  !          write(2216, '(6i8,2f15.8)') iBdir, iAdir, iG, -gpqdp2(1, iGvar), -gpqdp2(2, iGvar), -gpqdp2(3, iGvar) -1, psDens2ndOrd(iGvar, iAdir, 1, iBdir)
!          end do ! iG
!        end do ! iAdir
!      end do ! iBdir
!
!      Gpqext(:) = 0.
!      E2ndOrdIIat2(:, :) = cmplx(0., 0.)
!      do iBdir = 1, 3
!        do iAdir = 1, 3
!          do iG = 1, ngpqdp
!            Gpqext(1:3) = matmul(cell%bmat, gpqdp2(1:3, iG) + qpoint(1:3))
!            !todo actually we also need an exp(i (G + q) tau)
!            E2ndOrdIIat2(iAdir, iBdir) = E2ndOrdIIat2(iAdir, iBdir) + atoms%zatom(1) * fpi  / norm2(Gpqext(1:3))**2 &
!                                                                              & * psDens2ndOrd(iG, iAdir, 1, iBdir)
!          end do ! iG
!        end do ! iAdir
!      end do ! iBdir
!  !    write(*, '(a)') 'Eii2 for 0.75'
!  !    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(1, :)
!  !    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(2, :)
!  !    write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat(3, :)
!
!      ! For q = 0.75
!      call CalcIIEnerg2MatElem( atoms, cell, dimens, qpoint, ngpqdp, gpqdp2, iAtype, iBtype, iAatom, iBatom, E2ndOrdIIat2 )

      E2ndOrdIIat2(:, :) = 0.
      ngpqdp = 0
      qPoint(1:3) = nint((-qpts%bk(1:3, iqpt) + gShift(1:3)) * lcm)
      subIqpt = qIndex(qPoint(1), qPoint(2), qPoint(3))
 !     call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
      call CalcIIEnerg2(atoms, cell, dimens, qpts, stars, input, subIqpt, ngdp, gdp, E2ndOrdIIat2)
 !     deallocate(gpqdp2Ind, gpqdp)

      if (.false.) then
        write(*, '(a)') 'Eii2ordIIat'
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat2(1, :)
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat2(2, :)
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat2(3, :)

        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat2(1, :) - E2ndOrdIIat(1, :)
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat2(2, :) - E2ndOrdIIat(2, :)
        write(*, '(3(2(es16.8,1x),3x))') E2ndOrdIIat2(3, :) - E2ndOrdIIat(3, :)
      end if

      if (any ( abs( E2ndOrdIIat2(:, :) - E2ndOrdIIat(:, :) ) >1e-9) ) then
        passed = .false.
        !exit
      end if
    end do ! iqpt

    if (passed) then
      write (logUnit, '(a)')'                                                |__ passed!'
      write(logUnit, *)
    else
      write (logUnit, '(a)')'                                                |__ failed!'
      write(logUnit, *)
      call juDFT_warn('Check Eii2 lattice periodicity with resepct to q!', &
        & calledby='testEii2LatPeriodQ', hint='Debug.')
    end if

  end subroutine testEii2LatPeriodQ

 subroutine testEii2PsDens(atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, logUnit, memd_atom, rho0IR, rho0MT, gdp, mlh_atom, nmem_atom, clnu_atom)

   use m_jpVeff1, only : GenVeff1
   use m_jp2ndOrdQuant, only : GenPsDens2ndOrd
   use m_jpConstants, only : iu, fpi
   use m_jpPotDensHelper, only : genPertPotDensGvecs
   use m_juDFT_NOstopNO, only : juDFT_warn

   implicit none

   ! Type parameters
   type(t_atoms),                  intent(in)  :: atoms
   type(t_stars),                  intent(in)  :: stars
   type(t_cell),                   intent(in)  :: cell
   type(t_sphhar),                 intent(in)  :: lathar
   type(t_dimension),              intent(in)  :: dimens
   type(t_sym),                    intent(in)  :: sym
   type(t_input),                  intent(in)  :: input
   type(t_kpts),                   intent(in)  :: qpts

   ! Scalar parameter
   integer,                        intent(in)  :: ngdp
   integer,                        intent(in)  :: memd_atom
   integer,                        intent(in)  :: logUnit

   ! Array parameter
   complex,                        intent(in)  :: rho0IR(:, :)
   real,                           intent(in)  :: rho0MT(:, 0:, :, :)
   integer,                        intent(in)  :: gdp(:, :)
   integer,                        intent(in)  :: mlh_atom(:, 0:, :)
   integer,                        intent(in)  :: nmem_atom(0:, :)
   complex,                        intent(in)  :: clnu_atom(:, 0:, :)

   ! Scalar variables
   logical                                     :: harSw
   logical                                     :: extSw
   logical                                     :: xcSw
   logical                                     :: vExtFull
   integer                                     :: iDatom
   integer                                     :: iDtype
   integer                                     :: iqpt
   integer                                     :: G0index
   logical                                     :: testMode
   integer                                     :: iG
   integer                                     :: idirC
   integer                                     :: idirR
   real                                        :: noGqext
   integer                                     :: ngpqdp2km
   integer                                     :: ngpqdp
   logical                                     :: passed
   logical                                     :: vHarNum

   ! Array variable
   real                                        :: qpoint(3)
   real                                        :: Gext(3)
   real                                        :: Gpqext(3)
   complex,           allocatable              :: rho1PWDummy(:, :)
   complex,           allocatable              :: rho1MTDummy(:, :, :, :)
   complex,           allocatable              :: grRho0MTDummy(:, :, :, :)
   real,              allocatable              :: dKernMTGPtsDummy(:, :, :)
   complex,           allocatable              :: vxc1IRKernDummy(:)
   complex,           allocatable              :: ylmDummy(:, :)
   real,              allocatable              :: gWghtsDummy(:) ! gaussian weights belonging to gausPts
   complex,           allocatable              :: vEff1IR(:, :)
   complex,           allocatable              :: vEff1MT(:, :, :, :)
   complex,           allocatable              :: psDens2ndOrd(:, :, :, :)
   complex,           allocatable              :: gradVext1(:, :, :)
   integer,           allocatable              :: gpqdp(:, :)
   integer,           allocatable              :: gpqdp2Ind(:, :, :)
   complex,           allocatable              :: rho0IRDummy(:, :)
   complex,           allocatable              :: rho0MTDummy(:, :, :, :)
   integer                                     :: gpqdp2iLim(2, 3)

   write (logUnit, '(a)')'Compare Eii2 pseudodensity with the derivative of the pseudodensity used for Vext1!'
   write (logUnit, '(a)')'-----------------------------------------------------------------------------------'

   allocate( gWghtsDummy(dimens%nspd), &
           & ylmDummy(dimens%nspd, ( atoms%lmaxd + 1 )**2), &
           & dKernMTGPtsDummy(dimens%nspd, atoms%jmtd, atoms%nat), &
           & vxc1IRKernDummy(ngdp), &
           & vEff1IR(ngdp, 3), &
           & vEff1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
           & grRho0MTDummy(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3), &
           & rho1PWDummy(ngdp, 3), &
           & rho1MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat) &
           & )
   allocate(  rho0IRDummy(ngdp, 1) )
   allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )

   passed = .true.

   harSw = .false.
   extSw = .true.
   xcSw = .false.
   vExtFull = .false.
   vHarNum = .false.

   do iqpt = 1, qpts%nkpt

     qpoint = qpts%bk(1:3, iqpt)

     iDtype = 1
     iDatom = 1

     call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpts%bk(1:3, iqpt), gpqdp, gpqdp2Ind, gpqdp2iLim )

     allocate(gradVext1(ngpqdp, 3, 3))

     gWghtsDummy(:) = 0.
     ylmDummy(:, :) = cmplx(0., 0.)
     dKernMTGPtsDummy(:, :, :) = 0.
     vxc1IRKernDummy(:) = cmplx(0., 0.)
     vEff1IR(:, :) = cmplx(0., 0.)
     vEff1MT(:, :, :, :) = cmplx(0., 0.)
     grRho0MTDummy(:, :, :, :) = cmplx(0., 0.)
     rho1MTDummy(:, :, :, :) = cmplx(0., 0.)
     rho1PWDummy(:, :) = cmplx(0., 0.)
     gradVext1(:, :, :) = cmplx(0., 0.)
     rho0IRDummy(:, :) = cmplx(0., 0.)
     rho0MTDummy(:, :, :, :) = cmplx(0., 0.)

     call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRDummy, rho0MTDummy, rho1PWDummy, &
       & rho1MTDummy(:, :, :, :), grRho0MTDummy, gdp, vEff1IR, vEff1MT, vxc1IRKernDummy, ylmDummy, dKernMTGPtsDummy, gWghtsDummy, &
       & iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

     do iG = 1, ngpqdp
       Gpqext(1:3) = matmul(cell%bmat(1:3, 1:3), gpqdp(1:3, iG) + qpoint(1:3))
       do idirC = 1, 3
         do idirR = 1, 3
           gradVext1(iG, idirR, idirC) = iu * Gpqext(idirR) * vEff1IR(iG, idirC)
         end do ! idirR
       end do ! idirC
     end do ! iG

     testMode = .true.
     call GenPsDens2ndOrd( atoms, cell, dimens, ngpqdp, G0index, gpqdp, qpoint, psDens2ndOrd, testMode )

     ! Probably we do not need a shift of the Gset into -q
     do idirC = 1, 3
       do idirR = 1, 3
         do iG = 1, ngpqdp
           Gpqext(1:3) = matmul(cell%bmat(1:3, 1:3), gpqdp(1:3, iG) + qpoint(1:3))
           noGqext = norm2(Gpqext)
           if (abs(-gradVext1(iG, idirR, idirC) / fpi * noGqext**2 - psDens2ndOrd(iG, idirR, 1, idirC)) > 1e-8) passed = .false.
           if (.false.) then
             write(4000, '(3i8,2f15.8)') idirR, idirC, iG, -gradVext1(iG, idirR, idirC) / fpi * noGqext**2
             write(4001, '(3i8,2f15.8)') idirR, idirC, iG, psDens2ndOrd(iG, idirR, 1, idirC)
           end if
         end do ! iG
       end do ! idirR
     end do ! idirC

     deallocate(gradVext1)
   end do ! iqpt

   if (passed) then
     write (logUnit, '(a)')'                                                                                  |__ passed!'
     write(logUnit, *)
   else
     write (logUnit, '(a)')'                                                                                  |__ failed!'
     write(logUnit, *)
     call juDFT_warn('Compare Eii2 pseudodensity with the derivative of the pseudodensity used for Vext1!', &
       & calledby='testEii2PsDens', hint='Debug.')
   end if

 end subroutine testEii2PsDens

end module m_jpTestDynMatHF
