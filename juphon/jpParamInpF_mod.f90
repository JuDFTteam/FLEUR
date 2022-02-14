!----------------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!----------------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Contains IO routines for inputfile
!
!> @author
!> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
!>
!> @brief
!> Contains IO routines for inputfile juPhon.inp
!>
!> @note
!> Will be substituted by the xml input developed for FLEUR in the future. Therefore only minimal effort achieving modern code,
!> readability and documentation will be put into here.
!----------------------------------------------------------------------------------------------------------------------------------------
module m_jpParamInpF

  use mod_juPhonUtils, only : fopen, fclose

  implicit none

  character(8), save :: write_section = ' '

  interface getkey
    module procedure  getkeyScalInt, getkeyArrInt, getkeyScalReal, getkeyArrReal, getkeyLogical, getkeySingChar, getkeyArrChar
  end interface

  contains

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Read parameter from input file Juphon.inp.
  !>
  !> @details
  !> This concept of parsing strings from input file was adapted from ChriNOstopNOh Friedrich. A lot of code from him was recycled.
  !>
  !> @param[out]     addQs: additional phonon wavevectors q, stored sequentially in one-dimensional array
  !> @param[out] kSetShift: shift of k-point set
  !> @param[out]   kSetDim: dimension of k-point set
  !> @param[out]   qSetDim: dimension of q-point set
  !> @param[out]   kModeSw: k-point generation mode switch
  !> @param[out]  noPtsCon: number of points to be checked in continuity tests of densities and potentials
  !> @param[out]     paPoX: x-direction of path through unit cell for path potential-gradient test
  !> @param[out]     paPoY: y-direction of path through unit cell for path potential-gradient test
  !> @param[out]     paPoZ: z-direction of path through unit cell for path potential-gradient test
  !> @param[out]     harSw: Exchange-correlation potential switch for path potential-gradient test
  !> @param[out]     extSw: Exchange-correlation potential switch for path potential-gradient test
  !> @param[out]      xcSw: Exchange-correlation potential switch for path potential-gradient test
  !>
  !> @todo rename variables for inputFile
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine ReadInpFile( addQs, kSetDim, qSetDim, kModeSw, kSetShift, iqpt, imixJP, mixAlpha, noPtsCon, paPoX, paPoY, paPoZ, harSw, &
      & extSw, xcSw, writeKpqArraySw, calcEigenVec, oneSternhCycle, recipLengthUnit, onlyTests, testCompareGrVeff0FleurSw, testVeff1Sw, testUnfoldStarsSw, testRadDerivativeSw,&
      & testGauntCoeffSw, testGradLhExpandFuncSw, testContGrVeff0Sw, testWarpingSw, testSternhHSMEtestSw,                          &
      & testSternhSchroedPertTheoSw, testz1Phi0ContSw, testRho1BasCorrSw, testPlotRho03Dsw, testRadSolSw, testKptsWeightSw,        &
      & testCountValElecSw, testVeff0ContSw, testrho0ContSw, testBackRotMTCoordSysSw, testPsi0ContSw, testOverlapSw, testGradRho0PathSw,   &
      & testEii2LatPeriodQSw, testVarphiHepsVarphiSw, testRho1IRsw, testRho1MTsw, testsActivated,               &
      & test1st2ndPulDynMatEps1, test1st2ndPulDynMatCancel, test3rdPulDynMatCancel, testIntVeff1Rho1Val,     &
      & testGrPsiPsiMatElem, testCompareSurfInt, testSplitMTSurfIntSterh,                    &
      & testVeff1IRMESternh, testEps1q0, testVeff1IRMatqBackFold, testVeff1IRqLatPeriod, testGrMatElemPsiHepsPsiGaussTheo, testPsiHepsTildePsi, testGoldsteinRemaining, testR2orNotWfMtGradNgrNTensGrOvls, testComp3ArgSFIntsSw, testComp2ArgSFIntsSw, testGoldsteinSurfSw, testComp2ArgGrSFIntsSw, testIRIntegralSw, testIR3rdMatElemSw, testActionHgrPhiSw, testXCintegrals, testEii2PsDens )

    implicit none

    ! Array parameters
    real,         allocatable, intent(out) :: addQs(:)
    real,         allocatable, intent(out) :: kSetShift(:)
    integer,      allocatable, intent(out) :: kSetDim(:)
    integer,      allocatable, intent(out) :: qSetDim(:)

    ! Scalar parameters
    logical,                   intent(out) :: kModeSw
    logical,                   intent(out) :: writeKpqArraySw
    logical,                   intent(out) :: calcEigenVec
    logical,                   intent(out) :: oneSternhCycle
    logical,                   intent(out) :: recipLengthUnit
    logical,                   intent(out) :: onlyTests
    logical,                   intent(out) :: testsActivated
    integer,                   intent(out) :: noPtsCon
    integer,                   intent(out) :: iqpt
    integer,                   intent(out) :: imixJP
    real,                      intent(out) :: mixAlpha
    real,                      intent(out) :: paPoX
    real,                      intent(out) :: paPoY
    real,                      intent(out) :: paPoZ
    logical,                   intent(out) :: harSw
    logical,                   intent(out) :: extSw
    logical,                   intent(out) :: xcSw
    logical,                   intent(out) :: testCompareGrVeff0FleurSw
    logical,                   intent(out) :: testVeff1Sw
    logical,                   intent(out) :: testUnfoldStarsSw
    logical,                   intent(out) :: testRadDerivativeSw
    logical,                   intent(out) :: testGauntCoeffSw
    logical,                   intent(out) :: testGradLhExpandFuncSw
    logical,                   intent(out) :: testContGrVeff0Sw
    logical,                   intent(out) :: testWarpingSw
    logical,                   intent(out) :: testSternhHSMEtestSw
    logical,                   intent(out) :: testSternhSchroedPertTheoSw
    logical,                   intent(out) :: testz1Phi0ContSw
    logical,                   intent(out) :: testRho1BasCorrSw
    logical,                   intent(out) :: testPlotRho03Dsw
    logical,                   intent(out) :: testRadSolSw
    logical,                   intent(out) :: testKptsWeightSw
    logical,                   intent(out) :: testCountValElecSw
    logical,                   intent(out) :: testVeff0ContSw
    logical,                   intent(out) :: testrho0ContSw
    logical,                   intent(out) :: testBackRotMTCoordSysSw
    logical,                   intent(out) :: testPsi0ContSw
    logical,                   intent(out) :: testOverlapSw
    logical,                   intent(out) :: testGradRho0PathSw
    logical,                   intent(out) :: testEii2LatPeriodQSw
    logical,                   intent(out) :: testVarphiHepsVarphiSw
    logical,                   intent(out) :: testRho1IRsw
    logical,                   intent(out) :: testRho1MTsw
    logical,                   intent(out) :: test1st2ndPulDynMatEps1
    logical,                   intent(out) :: test1st2ndPulDynMatCancel
    logical,                   intent(out) :: test3rdPulDynMatCancel
    logical,                   intent(out) :: testIntVeff1Rho1Val
    logical,                   intent(out) :: testGrPsiPsiMatElem
    logical,                   intent(out) :: testCompareSurfInt
    logical,                   intent(out) :: testSplitMTSurfIntSterh
    logical,                   intent(out) :: testVeff1IRMESternh
    logical,                   intent(out) :: testEps1q0
    logical,                   intent(out) :: testVeff1IRMatqBackFold
    logical,                   intent(out) :: testVeff1IRqLatPeriod
    logical,                   intent(out) :: testGrMatElemPsiHepsPsiGaussTheo
    logical,                   intent(out) :: testPsiHepsTildePsi
    logical,                   intent(out) :: testGoldsteinRemaining
    logical,                   intent(out) :: testR2orNotWfMtGradNgrNTensGrOvls
    logical,                   intent(out) :: testComp3ArgSFIntsSw
    logical,                   intent(out) :: testGoldsteinSurfSw
    logical,                   intent(out) :: testComp2ArgSFIntsSw
    logical,                   intent(out) :: testComp2ArgGrSFIntsSw
    logical,                   intent(out) :: testIRIntegralSw
    logical,                   intent(out) :: testIR3rdMatElemSw
    logical,                   intent(out) :: testActionHgrPhiSw
    logical,                   intent(out) :: testXCintegrals
    logical,                   intent(out) :: testEii2PsDens

    ! Scalar variable
    integer                                :: inpUnit=101 ! unit for JuPhon.inp


    ! character(1)                           :: dummyChar
    ! character(1), allocatable              :: dummyCharArr(:)

    call fopen ( inpUnit, 'JuPhon.inp', status='old', action='read', form='formatted' )

    ! todo reorder the sections
    call getkeyArrReal( inpUnit,'additionalQs', addQs, writeout=.false., section='PHONON', default=[0., 0., 0.] )
    call getkeyArrInt( inpUnit, 'dimKptSet', kSetDim, writeout=.false., section='K-POINTS', default=[4, 4, 4] )
    call getkeyArrInt( inpUnit, 'dimQptSet', qSetDim, writeout=.false., section='K-POINTS', default=[2, 2, 2] )
    call getkeyLogical( inpUnit,'createKpoints', kModeSw, writeout=.false., section='K-POINTS', default=.false. )
    call getkeyLogical( inpUnit,'writeKpqArray', writeKpqArraySw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'Hartree contribution', harSw, writeout=.false., section='POTENTIAL', default=.false. )
    call getkeyLogical( inpUnit,'external contribution', extSw, writeout=.false., section='POTENTIAL', default=.false. )
    call getkeyLogical( inpUnit,'xc contribution', xcSw, writeout=.false., section='POTENTIAL', default=.false. )
    call getkeyLogical( inpUnit,'calcEigenVec', calcEigenVec, writeout=.false., section='DYNAMICAL MATRIX', default=.false. )
    call getkeyLogical( inpUnit,'recipLengthUnit', recipLengthUnit, writeout=.false., section='DYNAMICAL MATRIX', default=.false. )
    call getkeyLogical( inpUnit,'oneSternhCycle', oneSternhCycle, writeout=.false., section='STERNHEIMER', default=.false. )
    call getkeyLogical( inpUnit,'onlyTests', onlyTests, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testsActivated', testsActivated, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testCompareGrVeff0FleurSw', testCompareGrVeff0FleurSw , writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testVeff1Sw', testVeff1Sw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testUnfoldStarsSw', testUnfoldStarsSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testRadDerivativeSw', testRadDerivativeSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testGauntCoeffSw', testGauntCoeffSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testGradLhExpandFuncSw', testGradLhExpandFuncSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testContGrVeff0Sw', testContGrVeff0Sw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testWarpingSw', testWarpingSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testSternhHSMEtestSw', testSternhHSMEtestSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testSternhSchroedPertTheoSw', testSternhSchroedPertTheoSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testz1Phi0ContSw', testz1Phi0ContSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testRho1BasCorrSw', testRho1BasCorrSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testPlotRho03Dsw ', testPlotRho03Dsw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testRadSolSw', testRadSolSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testKptsWeightSw', testKptsWeightSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testCountValElecSw', testCountValElecSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testVeff0ContSw', testVeff0ContSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testrho0ContSw', testrho0ContSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testBackRotMTCoordSysSw', testBackRotMTCoordSysSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testPsi0ContSw', testPsi0ContSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testOverlapSw', testOverlapSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testGradRho0PathSw', testGradRho0PathSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testEii2LatPeriodQSw', testEii2LatPeriodQSw, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testVarphiHepsVarphiSw', testVarphiHepsVarphiSw , writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testRho1IRsw',  testRho1IRsw , writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testRho1MTsw',  testRho1MTsw , writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'test1st2ndPulDynMatEps1', test1st2ndPulDynMatEps1, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'test1st2ndPulDynMatCancel', test1st2ndPulDynMatCancel, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'test3rdPulDynMatCancel',  test3rdPulDynMatCancel, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testIntVeff1Rho1Val',  testIntVeff1Rho1Val, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testGrPsiPsiMatElem',  testGrPsiPsiMatElem, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testCompareSurfInt',  testCompareSurfInt, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testSplitMTSurfIntSterh', testSplitMTSurfIntSterh, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testVeff1IRMESternh',  testVeff1IRMESternh, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testEps1q0',  testEps1q0, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testVeff1IRMatqBackFold',  testVeff1IRMatqBackFold, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testVeff1IRqLatPeriod',  testVeff1IRqLatPeriod, writeout=.false., section='DEBUG', default=.false. )
    call getkeyLogical( inpUnit,'testGrMatElemPsiHepsPsiGaussTheo', testGrMatElemPsiHepsPsiGaussTheo, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testPsiHepsTildePsi', testPsiHepsTildePsi, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testGoldsteinRemaining', testGoldsteinRemaining, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testR2orNotWfMtGradNgrNTensGrOvls', testR2orNotWfMtGradNgrNTensGrOvls, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testComp3ArgSFIntsSw', testComp3ArgSFIntsSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testComp2ArgSFIntsSw', testComp2ArgSFIntsSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testGoldsteinSurfSw', testGoldsteinSurfSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testComp2ArgGrSFIntsSw', testComp2ArgGrSFIntsSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testIRIntegralSw', testIRIntegralSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testIR3rdMatElemSw', testIR3rdMatElemSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testActionHgrPhiSw', testActionHgrPhiSw, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testXCintegrals', testXCintegrals, writeout=.false., section='DEBUG', default=.false.)
    call getkeyLogical( inpUnit,'testEii2PsDens', testEii2PsDens, writeout=.false., section='DEBUG', default=.false.)
    ! todo kSetShift is not implemented yet
    call getkeyArrReal( inpUnit, 'kPSetShift', kSetShift, writeout=.false., section='K-POINTS', default=[0., 0., 0.] )
    call getkeyScalInt( inpUnit, 'Pts2ChkCont', noPtsCon, writeout=.false., section='DEBUG',  default=10 )
    call getkeyScalReal( inpUnit, 'paPoXCo', paPoX, writeout=.false., section='PLOT POTENTIAL', default=1.0)
    call getkeyScalReal( inpUnit, 'paPoYCo', paPoY, writeout=.false., section='PLOT POTENTIAL', default=0.0)
    call getkeyScalReal( inpUnit, 'paPoZCo', paPoZ, writeout=.false., section='PLOT POTENTIAL', default=0.0)
    call getkeyScalInt( inpUnit, 'method', imixJP, writeout=.false., section='MIXING',  default=1 )
    call getkeyScalReal( inpUnit, 'alpha', mixAlpha, writeout=.false., section='MIXING', default=0.05)
    call getkeyScalInt( inpUnit, 'qIndex', iqpt, writeout=.false., section='PHONON',  default=1 )

    ! todo put new variables also to logfile!
    ! call getkeyArrInt(  inpUnit,'duArrInt', dummyintarr,  writeout=.true.,  section='PHON')
    ! call getkeyScalReal(inpUnit,'dScalRea', dummyReal,    writeout=.true.,  section='PHON')
    ! call getkeyScalInt( inpUnit,'duScalIn', dummyInt,     writeout=.true.,  section='PHON', status=inputStatus, default=100)
    ! call getkeyLogical( inpUnit,'dummyLog', dummyLog,     writeout=.true.,  section='PHON')
    ! call getkeySingChar(inpUnit,'dummyCha', dummyChar,    writeout=.true.,  section='PHON')
    ! call getkeyArrChar( inpUnit,'dummyStr', dummyCharArr, writeout=.true.,  section='PHON')

    call fclose( inpUnit )

  end subroutine ReadInpFile

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses scalar integer from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeyScalInt(iunit, key, arg, section, default, status, writeout, mini, maxi, mine, maxe)

    implicit none

    integer,      intent(in)                :: iunit
    character(*), intent(in)                :: key
    character(*), intent(in),   optional    :: section
    logical,      intent(in),   optional    :: writeout
    integer,      intent(out),  optional    :: status

    integer,      intent(out)               :: arg
    integer,      intent(in),   optional    :: default
    integer,      intent(in),   optional    :: mini, mine
    integer,      intent(in),   optional    :: maxi, maxe

    integer                                 :: ios, ind, ind1, narg, i
    character(100)                          :: line1
    character(65536)                        :: line
    logical                                 :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /= ' ') then
        write(6, '(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    do
      line = ' '
 1    read(iunit, '(A)', iostat=ios) line1
      if(ios /= 0 .or. line1 .eq. 'EXIT') exit
      ind = index(line1, '#')         ! Remove
      if(ind /= 0) line1(ind:) = ' ' ! comments
      if(line1(80:) /= ' ') then
        write (6, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1, len_trim(line1))
      if(len_trim(line) + ind > len(line)) then
        write (6, *) 'getkey: Increase length of line.'
        NOstopNO
      endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line) + 1:) = line1(:ind - 1)    !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line) + 1:) = line1              !
      endif                                            !
      line = adjustl(line) ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8, len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /= section) then
            if(write_section /= ' ') write(6,'(A)') 'END'
            write(6,'(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line == 'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:),'#') + ind - 1 ! truncate comments
        if(ind1 == ind - 1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) == ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6,'(A,$)') '  '
            write(6,'(A)') key
          endif
          goto 2
        endif

        if(present(section)) then
          call getvalScalInt(line(ind:ind1), arg, narg, key, section)
        else
          call getvalScalInt(line(ind:ind1), arg, narg, key)
        endif

        if(present(mini)) call checkScalInt(  arg, mini, -1, key)
        if(present(maxi)) call checkScalInt(  arg, maxi,  1, key)
        if(present(mine)) call check_eScalInt(arg, mine, -1, key)
        if(present(maxe)) call check_eScalInt(arg, maxe,  1, key)

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)

      write(0,'(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6,'(A,$)') ' (section ' // trim(section) // ')'
      write(0,'(A)') ' not found.'
      NOstopNO
    endif

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(section)) write(6, '(A,$)') '  '
      write(6, '(A,$)') key // '       '(:7-len_trim(key))

      write(6, '(''  '',$)')
      call write_int(arg)
      write(6,*)
    endif
  end subroutine getkeyScalInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses integer array from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeyArrInt(iunit, key, arg, section, default, status, writeout, mini, maxi, mine, maxe)

    implicit none

    integer,      intent(in)                  :: iunit
    character(*), intent(in)                  :: key
    character(*), intent(in),   optional      :: section
    logical,      intent(in),   optional      :: writeout
    integer,      intent(out),  optional      :: status

    integer,                    allocatable   :: arg(:)
    integer,      intent(in),   optional      :: default(:)
    integer,      intent(in),   optional      :: mini, mine
    integer,      intent(in),   optional      :: maxi, maxe

    integer                                   :: ios, ind, ind1, narg, i
    character(100)                            :: line1
    character(65536)                          :: line
    logical                                   :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /= ' ') then
        write(6,'(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    narg = size(arg) ; if(narg == 0 .and. allocated(arg)) deallocate (arg)
    if(present(default) .and. allocated(arg)) then
      if(size(default) /= narg) then
        write (6, *) 'getkey: Number of arguments and default values do not agree.'
        NOstopNO
      endif
    endif

    do
      line = ' '
 1    read(iunit, '(A)', iostat=ios) line1
      if(ios /= 0 .or. line1 == 'EXIT') exit
      ind = index(line1, '#')         ! Remove
      if(ind /= 0) line1(ind:) = ' ' ! comments
      if(line1(80:) /= ' ') then
        write (6, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1, len_trim(line1))                     !
      if(len_trim(line) + ind > len(line)) then             !
       write (6, *) 'getkey: Increase length of line.'
       NOstopNO!
      endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line) + 1:) = line1(:ind - 1)        !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line) + 1:) = line1                !
      endif                                            !
      line = adjustl(line) ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8,len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /= section) then
            if(write_section /= ' ') write(6,'(A)') 'END'
            write(6,'(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line == 'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:), '#') + ind - 1 ! truncate comments
        if(ind1 == ind - 1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) == ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6,'(A,$)') '  '
            write(6,'(A)') key
          endif
          goto 2
        endif


        if(present(section)) then
          call getvalArrInt(line(ind:ind1), arg, narg, key, section)
        else
          call getvalArrInt(line(ind:ind1), arg, narg, key)
        endif

        if(present(mini)) call checkArrInt(   arg, mini, -1, key)
        if(present(maxi)) call checkArrInt(   arg, maxi,  1, key)
        if(present(mine)) call check_eArrInt( arg, mine, -1, key)
        if(present(maxe)) call check_eArrInt( arg, maxe,  1, key)

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)

      write(0, '(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6,'(A,$)') ' (section ' // trim(section) // ')'
      write(0, '(A)') ' not found.'
      NOstopNO
    endif

    if(.not.allocated(arg)) allocate(arg(size(default)))

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(section)) write(6,'(A,$)') '  '
      write(6,'(A,$)') key //'       ' (:7-len_trim(key))

      do i = 1,narg
        write(6, '(''  '',$)') ; call write_int(arg(i))
      enddo
      write(6, *)
    endif
  end subroutine getkeyArrInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses scalar real from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeyScalReal(iunit, key, arg, section, default, status, writeout, mini, maxi, mine, maxe)

    implicit none

    integer,      intent(in)                :: iunit
    character(*), intent(in)                :: key
    character(*), intent(in),   optional    :: section
    logical,      intent(in),   optional    :: writeout
    integer,      intent(out),  optional    :: status

    real(8),      intent(out)               :: arg
    real(8),      intent(in),   optional    :: default
    real(8),      intent(in),   optional    :: mini, mine
    real(8),      intent(in),   optional    :: maxi, maxe

    integer                                 :: ios, ind, ind1, narg, i
    character(100)                          :: line1
    character(65536)                        :: line
    logical                                 :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /= ' ') then
        write(6,'(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    do
      line = ' '
 1    read(iunit, '(A)', iostat=ios) line1
      if(ios /= 0 .or. line1 == 'EXIT') exit
      ind = index(line1, '#')                          ! Remove
      if(ind /= 0) line1(ind:) = ' '                   ! comments
      if(line1(80:) /= ' ') then
        write (6, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1, len_trim(line1))
      if(len_trim(line) + ind > len(line)) then
       write (6, *) 'getkey: Increase length of line.'
       NOstopNO!
      endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line) + 1:) = line1(:ind - 1)    !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line) + 1:) = line1              !
      endif                                            !
      line = adjustl(line)                             ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8,len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /=  section) then
            if(write_section /= ' ') write(6,'(A)') 'END'
            write(6,'(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line ==  'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:),'#') + ind - 1 ! truncate comments
        if(ind1 == ind - 1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) .eq. ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6, '(A,$)') '  '
            write(6, '(A)') key
          endif
          goto 2
        endif

        if(present(section)) then
          call getvalScalReal (line(ind:ind1), arg, narg, key, section)
        else
          call getvalScalReal (line(ind:ind1), arg, narg, key)
        endif

        if(present(mini)) call checkScalReal   (arg, mini, -1, key)
        if(present(maxi)) call checkScalReal   (arg, maxi,  1, key)
        if(present(mine)) call check_eScalReal (arg, mine, -1, key)
        if(present(maxe)) call check_eScalReal (arg, maxe,  1, key)

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)

      write(0,'(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6, '(A,$)') ' (section ' // trim(section) // ')'
      write(0,'(A)') ' not found.'
      NOstopNO
    endif

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(section)) write(6, '(A,$)') '  '
      write(6, '(A,$)') key // '       '(:7-len_trim(key))

      write(6, '(''  '',$)') ; call write_real(arg)

      write(6, *)
    endif
  end subroutine getkeyScalReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses real array from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeyArrReal(iunit, key, arg, section, default, status, writeout, mini, maxi, mine, maxe)

    implicit none

    integer,      intent(in)                :: iunit
    character(*), intent(in)                :: key
    character(*), intent(in),   optional    :: section
    logical,      intent(in),   optional    :: writeout
    integer,      intent(out),  optional    :: status

    real(8),                    allocatable :: arg(:)
    real(8),      intent(in),   optional    :: default(:)
    real(8),      intent(in),   optional    :: mini, mine
    real(8),      intent(in),   optional    :: maxi, maxe

    integer                                 :: ios, ind, ind1, narg, i
    character(100)                          :: line1
    character(65536)                        :: line
    logical                                 :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /= ' ') then
        write(6,'(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    narg = size(arg) ; if(narg == 0 .and. allocated(arg)) deallocate (arg)
    if(present(default) .and. allocated(arg)) then
      if(size(default) /= narg) then
        write (*, *) 'getkey: Number of arguments and default values do not agree.'
        NOstopNO
      endif
    endif

    do
      line = ' '
 1    read(iunit,'(A)', iostat=ios) line1
      if(ios /= 0 .or. line1 ==  'EXIT') exit
      ind = index(line1, '#')         ! Remove
      if(ind /= 0) line1(ind:) = ' ' ! comments
      if(line1(80:) /= ' ') then
        write (*, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1, len_trim(line1))                     !
      if(len_trim(line) + ind > len(line)) then         !
       write (6, *) 'getkey: Increase length of line.'
       NOstopNO
      endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line) + 1:) = line1(:ind - 1)        !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line) + 1:) = line1                !
      endif                                            !
      line = adjustl(line) ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8,len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /= section) then
            if(write_section /=  ' ') write(6,'(A)') 'END'
            write(6,'(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line == 'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:),'#') + ind - 1 ! truncate comments
        if(ind1 == ind - 1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) == ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6,'(A,$)') '  '
            write(6,'(A)') key
          endif
          goto 2
        endif

        if(present(section)) then
          call getvalArrReal(line(ind:ind1), arg, narg, key, section)
        else
          call getvalArrReal(line(ind:ind1), arg, narg, key)
        endif

        if(present(mini)) call checkArrReal(  arg, mini, -1, key)
        if(present(maxi)) call checkArrReal(  arg, maxi,  1, key)
        if(present(mine)) call check_eArrReal(arg, mine, -1, key)
        if(present(maxe)) call check_eArrReal(arg, maxe,  1, key)

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)

      write(0, '(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6,'(A,$)') ' (section ' // trim(section) // ')'
      write(0, '(A)') ' not found.'
      NOstopNO
    endif

    if(.not.allocated(arg)) allocate(arg(size(default)))

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(section)) write(6, '(A,$)') '  '
      write(6, '(A,$)') key // '       '(:7-len_trim(key))

      do i = 1,narg
        write(6,'(''  '',$)') ; call write_real(arg(i))
      enddo
      write(6,*)
    endif
  end subroutine getkeyArrReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses logical from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeyLogical(iunit, key, arg, section, default, status, writeout)

    implicit none

    integer,      intent(in)                :: iunit
    character(*), intent(in)                :: key
    character(*), intent(in),   optional    :: section
    logical,      intent(in),   optional    :: writeout
    integer,      intent(out),  optional    :: status

    logical,      intent(out)               :: arg
    logical,      intent(in),   optional    :: default

    integer                                 :: ios, ind, ind1, narg, i
    character(100)                          :: line1
    character(65536)                        :: line
    logical                                 :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /= ' ') then
        write(6,'(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    if(present(default)) then
      arg = default
    else
      arg = .false.
    endif


    do
      line = ' '
 1    read(iunit, '(A)', iostat=ios) line1
      if(ios /= 0 .or. line1 == 'EXIT') exit
      ind = index(line1,'#')         ! Remove
      if(ind  /=  0) line1(ind:) = ' ' ! comments
      if(line1(80:) /= ' ') then
        write (6, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1, len_trim(line1))
      if(len_trim(line) + ind > len(line)) then
       write (6, *) 'getkey: Increase length of line.'
       NOstopNO
      endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line) + 1:) = line1(:ind - 1)        !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line) + 1:) = line1                !
      endif                                            !
      line = adjustl(line) ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8, len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /= section) then
            if(write_section /= ' ') write(6,'(A)') 'END'
            write(6,'(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line == 'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:),'#') + ind - 1 ! truncate comments
        if(ind1 == ind - 1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) == ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6,'(A,$)') '  '
            write(6,'(A)') key
          endif
          goto 2
        endif

        arg = .not.arg

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)
      write(6, * )'getkey: No default given for logical key ' // trim(key) // '. (bug?)'
      write(0,'(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6,'(A,$)') ' (section ' // trim(section) // ')'
      write(0,'(A)') ' not found.'
      NOstopNO
    endif

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(default)) then
        if(arg.neqv.default) then
          if(present(section)) write(6,'(A,$)') '  '
          write(6, '(A)') key
        endif
      else
        if(arg) then
          if(present(section)) write(6,'(A,$)') '  '
          write(6, '(A)') key
        endif
      endif
    endif
  end subroutine getkeyLogical

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses scalar character from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeySingChar(iunit, key, arg, section, default, status, writeout)

    implicit none

    integer,      intent(in)                :: iunit
    character(*), intent(in)                :: key
    character(*), intent(in),   optional    :: section
    logical,      intent(in),   optional    :: writeout
    integer,      intent(out),  optional    :: status

    character(*), intent(out)               :: arg
    character(*), intent(in),   optional    :: default

    integer                                 :: ios, ind, ind1, narg, i
    character(100)                          :: line1
    character(65536)                        :: line
    logical                                 :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /=  ' ') then
        write(6,'(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    do
      line = ' '
 1    read(iunit, '(A)', iostat=ios) line1
      if(ios /=  0 .or. line1 == 'EXIT') exit
      ind = index(line1, '#')         ! Remove
      if(ind /= 0) line1(ind:) = ' ' ! comments
      if(line1(80:) /= ' ') then
        write (6, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1,len_trim(line1))                     !
      if(len_trim(line) + ind > len(line)) then             !
       write (6, *) 'getkey: Increase length of line.'
       NOstopNO!
     endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line)+1:) = line1(:ind-1)        !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line)+1:) = line1                !
      endif                                            !
      line = adjustl(line) ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8,len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /= section) then
            if(write_section /= ' ') write(6,'(A)') 'END'
            write(6, '(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line == 'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:), '#') + ind - 1 ! truncate comments
        if(ind1 == ind-1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) ==  ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6, '(A,$)') '  '
            write(6, '(A)') key
          endif
          goto 2
        endif

        if(present(section)) then
          call getvalSingChar(line(ind:ind1), arg, narg, key, section)
        else
          call getvalSingChar(line(ind:ind1), arg, narg, key)
        endif

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)

      write(0,'(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6,'(A,$)') ' (section ' // trim(section) // ')'
      write(0,'(A)') ' not found.'
      NOstopNO
    endif

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(section)) write(6,'(A,$)') '  '
      write(6,'(A,$)') key // '       '(:7-len_trim(key))
      write(6,'(3X,A,$)') trim(arg)
      write(6,*)
    endif
  end subroutine getkeySingChar

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Parses character array from file with unit iunit.
  !>
  !> @details
  !> The usage is
  !> call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>],[writeout=<writeout>]
  !>             [,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
  !>
  !> Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
  !> <dest> can be a single argument or an array of integer values, real(8) values or character strings.
  !> NOTE: The array must be an allocatable array!
  !>
  !> If <keyword> is not found, the program NOstopNOs with an error message
  !> unless one of <default> and <status> is given.
  !> If <default> is given, then <dest>:=<default>.
  !> If <status>  is given, then <status> takes one of the values
  !> 0 : <keyword> not found,
  !> 1 : <keyword> is found but no arguments are provided,
  !> 2 : <keyword> and arguments are found.
  !> In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
  !>
  !> <dest> can also be logical. Then, the occurrence of <keyword> switches its value
  !> to .true., or if <default> is given to .not.<default>.
  !>
  !> If <section> is given, <keyword> is searched for only between "SECTION <section>"
  !> and "END".
  !>
  !> With mini, maxi, mine and maxe a range for values can be specified:
  !> <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.c
  !>
  !> Everything in the line after "#" is treated as a comment.
  !> Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
  !> keywords another routine should look through the input file and check the keywords.
  !>
  !> @note Input will be done by xml in the long run so no optimal effort is put into this routine
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getkeyArrChar(iunit, key, arg, section, default, status, writeout)

    implicit none

    integer,      intent(in)                :: iunit
    character(*), intent(in)                :: key
    character(*), intent(in),   optional    :: section
    logical,      intent(in),   optional    :: writeout
    integer,      intent(out),  optional    :: status

    character(*),               allocatable :: arg(:)
    character(*), intent(in),   optional    :: default(:)

    integer                                 :: ios, ind, ind1, narg, i
    character(100)                          :: line1
    character(65536)                        :: line
    logical                                 :: insection, searchkey, writeout1

    writeout1 = .false.
    if(present(writeout)) writeout1 = writeout

    if(.not.present(section)) then
      if(write_section /= ' ') then
        write(6,'(A)') 'END'
      endif
    endif

    insection = .false.
    if(present(section)) then
      searchkey = .false.
    else
      searchkey = .true.
    endif

    rewind(iunit)

    narg = size(arg) ; if(narg == 0 .and. allocated(arg)) deallocate (arg)
    if(present(default) .and. allocated(arg)) then
      if(size(default) /= narg) then
        write (6, *) 'getkey: Number of arguments and default values do not agree.'
        NOstopNO
      endif
    endif

    do
      line = ' '
 1    read(iunit, '(A)', iostat=ios) line1
      if(ios /= 0 .or. line1 == 'EXIT') exit
      ind = index(line1, '#')         ! Remove
      if(ind /= 0) line1(ind:) = ' ' ! comments
      if(line1(80:) /= ' ') then
        write (6, *) 'getkey: Line exceeds 80 columns.'
        NOstopNO
      endif
      ind = max(1, len_trim(line1))                     !
      if(len_trim(line) + ind > len(line)) then             !
     write (6, *) 'getkey: Increase length of line.'
     NOstopNO!
   endif
      if(line1(ind:ind) == '\') then                   !
        line(len_trim(line) + 1:) = line1(:ind - 1)        !
        goto 1                                         ! allow for multiple lines
      else                                             ! (continued with '\')
        line(len_trim(line) + 1:) = line1                !
      endif                                            !
      line = adjustl(line) ! Remove leading spaces
      if(line(:8) == 'SECTION ' .and. key(:min(8,len(key))) /= 'SECTION ') then
        insection = .true.
        if(searchkey) searchkey = .false.
        if(present(section)) then
          if(adjustl(line(8:)) == section) searchkey = .true.
          if(writeout1 .and. write_section /= section) then
            if(write_section /= ' ') write(6,'(A)') 'END'
            write(6,'(A)') 'SECTION ' // trim(section)
            write_section = trim(section) // '      '
          endif
        endif
        cycle
      endif
      if(line == 'END') then
        insection = .false.
        if(searchkey) searchkey = .false.
        if(.not.present(section)) searchkey = .true.
        cycle
      endif
      if(.not.searchkey) cycle
      ind = len_trim(key) + 1
      if(line(:ind) == key) then
        ind1 = index(line(ind:), '#') + ind - 1 ! truncate comments
        if(ind1 == ind-1) ind1 = len(line)     !
        if(present(status) .and. line(ind:ind1) == ' ') then
          status = 1
          if(writeout1) then
            if(present(section)) write(6,'(A,$)') '  '
            write(6,'(A)') key
          endif
          goto 2
        endif

        if(present(section)) then
          call getvalArrChar(line(ind:ind1), arg, narg, key, section)
        else
          call getvalArrChar(line(ind:ind1), arg, narg, key)
        endif

        goto 3 ! write out and leave
      endif
    enddo

    if(present(status)) status = 0

 2  continue

    if(.not.present(default)) then
      if(present(status)) return ! arg remains undefined on return (status=0)


      write(0,'(2A,$)') 'getkey: Keyword ', trim(key)
      if(present(section)) write(6, '(A,$)') ' (section ' // trim(section) // ')'
      write(0,'(A)') ' not found.'
      NOstopNO
    endif

    if(.not.allocated(arg)) allocate(arg(size(default)))

    arg = default
    return

 3  if(present(status)) status = 2
    if(writeout1) then

      if(present(section)) write(6, '(A,$)') '  '
      write(6, '(A,$)') key // '       ' (:7-len_trim(key))

      write(6, '('' '',$)')
      do i = 1,narg
        write(6, '(2X,A,$)') trim(arg(i))
      enddo

      write(6, *)
    endif
  end subroutine getkeyArrChar

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyScalInt to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getvalScalInt(line, arg, narg, key, section)

    use mod_juPhonUtils

    implicit none

    character(*), intent(in)               :: line, key
    character(*), intent(in), optional     :: section
    integer                                :: narg
    integer                                :: iword, i, j, ios
    logical                                :: isinteger
    character(80)                          :: word(2000), type

    integer,      intent(out)              :: arg
    type = ' integer'

    word  = ' '
    iword = 1
    do i = 1, len(line)
      if(line(i:i) /= ' ') then
        j                = len_trim(word(iword)) + 1
        if(j > len(word(1))) then
          write (6, *) 'getval: Increase length of word(:).'
          NOstopNO
        endif
        word(iword)(j:j) = line(i:i)
      else if(word(iword) /= ' ') then
        iword = iword + 1
        if(iword > size(word)) then
          write (6, *) 'getval: Increase dimension of word(:).'
          NOstopNO
        endif
      endif
    enddo
    if(word(iword) == ' ') iword = iword - 1
    narg = iword

    if(iword > 1) then
      write(0,'(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
      if(present(section)) write(0,'(A,$)') ' in section ' // trim(section)
      write(0,'(/A)')       '        but needed only one.'
      NOstopNO
    endif
    if(iword < 1) then
      write(0,'(A,$)') 'getval: Could not read' // trim(type) // ' argument after keyword ' // trim(key)
      if(present(section)) then
        write(0,'(A)') ' in section ' // trim(section)
      else
        write(0,*)
      endif
      NOstopNO
    endif

    if(.not.checkisinteger(word(1))) then
      write(0, '(A)') 'getval: argument ' // trim(word(1)) // ' is not an integer.'
      NOstopNO
    endif

    read(word(1),*,iostat=ios) arg
    if(ios /= 0) then
      write(0, '(A)') 'getval: Could not read' // trim(type) // ' argument from ' // trim(word(i))
      write (6, *) 'getval: Read error'
      NOstopNO
    endif
  end subroutine getvalScalInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrInt to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getvalArrInt(line, arg, narg, key, section)

    use mod_juPhonUtils

    implicit none

    character(*), intent(in)               :: line,key
    character(*), intent(in),  optional    :: section
    integer                                :: narg
    integer                                :: iword, i, j, ios
    logical                                :: isinteger
    character(80)                          :: word(2000), type

    integer,                   allocatable :: arg(:)
    type = ' integer'

    word  = ' '
    iword = 1
    do i = 1, len(line)
      if(line(i:i) /= ' ') then
        j                = len_trim(word(iword)) + 1
        if(j > len(word(1))) then
          write (6, *) 'getval: Increase length of word(:).'
          NOstopNO
        endif
        word(iword)(j:j) = line(i:i)
      else if(word(iword) /= ' ') then
        iword = iword + 1
        if(iword >= size(word)) then
          write (6, *) 'getval: Increase dimension of word(:).'
          NOstopNO
        endif
      endif
    enddo
    if(word(iword) == ' ') iword = iword - 1
    narg = iword

    if(allocated(arg)) then
      if(iword > size(arg,1)) then
        write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
        if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
        write(0, '(/A,I3)')    '        but needed only', size(arg,1)
        write (6, *) 'getval: Read more values than needed.'; NOstopNO
      else if(iword < size(arg,1)) then
        write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
        if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
        write(0, '(/A,I2)')    '        but needed', size(arg, 1)
        NOstopNO
      endif
    else
      allocate(arg(iword))
    endif
    do i = 1, size(arg, 1)

      if(.not.checkisinteger(word(i))) then
        write(0,'(A)') 'getval: argument ' // trim(word(i)) // ' is not an integer.'
        NOstopNO
      endif
      read(word(i), *, iostat=ios) arg(i)
      if(ios /= 0) then
        write(0, '(A)') 'getval: Could not read' // trim(type) // ' argument from ' // trim(word(i))
        NOstopNO
      endif
    enddo
  end subroutine getvalArrInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyScalReal to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getvalScalReal(line, arg, narg, key, section)

    implicit none

    character(*), intent(in)            :: line, key
    character(*), intent(in), optional  :: section
    integer                             :: narg
    integer                             :: iword, i, j, ios
    logical                             :: isinteger
    character(80)                       :: word(2000), type

    real(8),      intent(out)           :: arg
    type = ' real'

    word  = ' '
    iword = 1
    do i = 1, len(line)
      if(line(i:i) /= ' ') then
        j                = len_trim(word(iword)) + 1
        if(j > len(word(1))) then
          write (6, *) 'getval: Increase length of word(:).'
          NOstopNO
        endif
        word(iword)(j:j) = line(i:i)
      else if(word(iword) /= ' ') then
        iword = iword + 1
        if(iword > size(word)) then
          write (6, *) 'getval: Increase dimension of word(:).'
          NOstopNO
        endif
      endif
    enddo
    if(word(iword) == ' ') iword = iword - 1
    narg = iword

    if(iword > 1) then
      write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
      if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
      write(0, '(/A)')       '        but needed only one.'
      NOstopNO
    endif
    if(iword < 1) then
      write(0, '(A,$)') 'getval: Could not read' // trim(type) // ' argument after keyword ' // trim(key)
      if(present(section)) then
        write(0, '(A)') ' in section ' // trim(section)
      else
        write(0, *)
      endif
      NOstopNO
    endif
    read(word(1), *, iostat=ios) arg
    if(ios /=  0) then
      write(0, '(A)') 'getval: Could not read' // trim(type) // ' argument from ' // trim(word(i))
      write(6, *) 'getval: Read error'
    endif
  end subroutine getvalScalReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrReal to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getvalArrReal(line, arg, narg, key, section)

    implicit none

    character(*), intent(in)               :: line, key
    character(*), intent(in), optional     :: section
    integer                                :: narg
    integer                                :: iword, i, j, ios
    logical                                :: isinteger
    character(80)                          :: word(2000), type

    real(8),                   allocatable :: arg(:)
    type = ' real'

    word  = ' '
    iword = 1
    do i = 1, len(line)
      if(line(i:i) /= ' ') then
        j                = len_trim(word(iword)) + 1
        if(j > len(word(1))) then
        write (6, *) 'getval: Increase length of word(:).'
        NOstopNO
      endif
        word(iword)(j:j) = line(i:i)
      else if(word(iword) /= ' ') then
        iword = iword + 1
        if(iword > size(word)) then
          write (6, *) 'getval: Increase dimension of word(:).'
          NOstopNO
        endif
      endif
    enddo
    if(word(iword) == ' ') iword = iword - 1
    narg = iword

    if(allocated(arg)) then
      if(iword > size(arg,1)) then
        write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
        if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
        write(0, '(/A,I3)')    '        but needed only', size(arg, 1)
        write (6, *) 'getval: Read more values than needed.'; NOstopNO
      else if(iword < size(arg,1)) then
        write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
        if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
        write(0, '(/A,I2)')    '        but needed', size(arg, 1)
        NOstopNO
      endif
    else
      allocate(arg(iword))
    endif
    do i = 1, size(arg, 1)

      read(word(i), *, iostat=ios) arg(i)
      if(ios /= 0) then
        write(0,'(A)') 'getval: Could not read' // trim(type) // ' argument from ' // trim(word(i))
        NOstopNO
      endif
    enddo
  end subroutine getvalArrReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeySingChar to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getvalSingChar(line, arg, narg, key, section)

    implicit none

    character(*), intent(in)               :: line, key
    character(*), intent(in), optional     :: section
    integer                                :: narg
    integer                                :: iword, i, j, ios
    logical                                :: isinteger
    character(80)                          :: word(2000), type

    character(*), intent(out)              :: arg
    type = ' character'

    word  = ' '
    iword = 1
    do i = 1, len(line)
      if(line(i:i) /= ' ') then
        j                = len_trim(word(iword)) + 1
        if(j > len(word(1))) then
          write (6, *) 'getval: Increase length of word(:).'
          NOstopNO
        endif
        word(iword)(j:j) = line(i:i)
      else if(word(iword) /= ' ') then
        iword = iword + 1
        if(iword > size(word)) then
          write (6, *) 'getval: Increase dimension of word(:).'
          NOstopNO
        endif
      endif
    enddo
    if(word(iword) == ' ') iword = iword - 1
    narg = iword

    if(iword > 1) then
      write(0, '(A,I3,A,$)') 'getval: Read' ,iword, ' arguments after keyword ' // trim(key)
      if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
      write(0, '(/A)')       '        but needed only one.'
      NOstopNO
    endif
    if(iword < 1) then
      write(0, '(A,$)') 'getval: Could not read' // trim(type) // ' argument after keyword ' // trim(key)
      if(present(section)) then
        write(0, '(A)') ' in section ' // trim(section)
      else
        write(0, *)
      endif
      NOstopNO
    endif

    arg = word(1)
  end subroutine getvalSingChar

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrChar to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine getvalArrChar(line, arg, narg, key, section)

    implicit none

    character(*), intent(in)               :: line, key
    character(*), intent(in), optional     :: section
    integer                                :: narg
    integer                                :: iword, i, j, ios
    logical                                :: isinteger
    character(80)                          :: word(2000), type

    character(*), allocatable              :: arg(:)
    type = ' character'

    word  = ' '
    iword = 1
    do i = 1, len(line)
      if(line(i:i) /= ' ') then
        j                = len_trim(word(iword)) + 1
        if(j > len(word(1))) then
          write (6, *) 'getval: Increase length of word(:).'
          NOstopNO
        endif
        word(iword)(j:j) = line(i:i)
      else if(word(iword) /= ' ') then
        iword = iword + 1
        if(iword > size(word)) then
          write (6, *) 'getval: Increase dimension of word(:).'
          NOstopNO
        endif
      endif
    enddo
    if(word(iword) == ' ') iword = iword - 1
    narg = iword

    if(allocated(arg)) then
      if(iword > size(arg,1)) then
        write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
        if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
        write(0, '(/A,I3)')    '        but needed only', size(arg, 1)
        write (6, *) 'getval: Read more values than needed.'; NOstopNO
      else if(iword < size(arg,1)) then
        write(0, '(A,I3,A,$)') 'getval: Read', iword, ' arguments after keyword ' // trim(key)
        if(present(section)) write(0, '(A,$)') ' in section ' // trim(section)
        write(0, '(/A,I2)')    '        but needed', size(arg, 1)
        NOstopNO
      endif
    else
      allocate(arg(iword))
    endif
    do i = 1, size(arg,1)
      arg(i) = word(i)
    enddo
  end subroutine getvalArrChar

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyScalInt to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkScalInt(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    integer                  :: arg, minmax, arg1

    ldum = arg * mode > minmax * mode
    arg1 = arg

    i = (mode + 1) / 2 + 1
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_int(arg1) ; write(0, '(1X,A,$)') '<>'(i:i) // ' ' ; call write_int(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine checkScalInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrInt to parse value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkArrInt(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    integer                  :: arg(:), minmax, arg1

    ldum = any(arg * mode > minmax * mode) ; arg1 = mode * maxval(mode * arg)

    i = (mode + 1) / 2 + 1
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_int(arg1) ; write(0, '(1X,A,$)') '<>'(i:i) // ' ' ; call write_int(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine checkArrInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyScalReal to check value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkScalReal(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    real(8)                  :: arg, minmax, arg1
    ldum = arg * mode > minmax * mode      ; arg1 = arg

    i = (mode + 1) / 2 + 1
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_real(arg1) ; write(0, '(1X,A,$)') '<>'(i:i) // ' ' ; call write_real(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine checkScalReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrReal to check value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkArrReal(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    real(8)                  :: arg(:), minmax, arg1
    ldum = any(arg * mode > minmax * mode) ; arg1 = mode * maxval(mode * arg)

    i = (mode + 1) / 2 + 1
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_real(arg1) ; write(0, '(1X,A,$)') '<>'(i:i) // ' ' ; call write_real(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine checkArrReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyScalInt to check value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_eScalInt(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    integer                  :: arg, minmax, arg1
    ldum = arg * mode >=  minmax * mode      ; arg1 = arg

    i = mode + 2
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_int(arg1) ; write(0, '(1X,A,$)') '<=>='(i:i+1) // ' ' ; call write_int(minmax)
      write(0,'(A)')   ' .'
      NOstopNO
    endif
  end subroutine check_eScalInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrInt to check value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_eArrInt(arg, minmax, mode, key)

    implicit none

    integer, intent(in) :: mode
    logical             :: ldum
    character(*)        :: key
    integer             :: i

    integer             :: arg(:), minmax, arg1
    ldum = any(arg * mode >= minmax * mode) ; arg1 = mode * maxval(mode * arg)

    i = mode + 2
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_int(arg1) ; write(0, '(1X,A,$)') '<=>='(i:i+1) // ' ' ; call write_int(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine check_eArrInt

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyScalReal to check value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_eScalReal(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    real(8)                  :: arg,minmax,arg1
    ldum = arg * mode >= minmax * mode      ; arg1 = arg

    i = mode + 2
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_real(arg1) ; write(0, '(1X,A,$)') '<=>='(i:i+1) // ' ' ; call write_real(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine check_eScalReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of getkeyArrReal to check value for respective key.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine check_eArrReal(arg, minmax, mode, key)

    implicit none

    integer,      intent(in) :: mode
    logical                  :: ldum
    character(*)             :: key
    integer                  :: i

    real(8)                  :: arg(:), minmax, arg1
    ldum = any(arg * mode >= minmax * mode) ; arg1 = mode * maxval(mode * arg)

    i = mode + 2
    if(ldum) then
      write(0, '(A,$)') 'getkey: Value out of range after keyword ' // trim(key)
      write(0, '(A,$)') ': ' ; call write_real(arg1) ; write(0, '(1X,A,$)') '<=>='(i:i+1) // ' ' ; call write_real(minmax)
      write(0, '(A)')   ' .'
      NOstopNO
    endif
  end subroutine check_eArrReal

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of real value parser routines
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_int(arg)

    implicit none

    character(20)              :: form, line
    integer,        intent(in) :: arg
    integer                    :: minusspace

    write(line, *) arg
    if (arg >= 0) then
      minusspace = 1
    else
      minusspace = 0
    end if
    write(form, '(A,I3,A)') '(I', len_trim(adjustl(line)) + minusspace, ',$)'
    write(6, form) arg
  end subroutine write_int

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Auxiliary routine of real value parser routines
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine write_real(arg)

    implicit none

    character(20)             :: form
    character(80)             :: line
    real(8),       intent(in) :: arg
    real(8)                   :: lg
!    write(*,'(A,F10.5,$)') '!!!',arg ; return

    if (nint(arg) == arg) then ! is already an integer
      call write_int(nint(arg))
      return
    endif

    if (arg == 0) then
      lg = 0
    else
      lg = log(abs(arg)) / log(10d0)
    endif

    if(lg >= 0) then
      write(form, '(A,I3.3,A,I3.3,A)') '(F', int(lg) + 3 + 15, '.', 15, ',$)'
    else if(lg >= -3) then
      write(form, '(A,I3.3,A,I3.3,A)') '(F', int(-lg) + 3 + 15, '.', int(-lg) + 15,',$)'
    else
      form = '(ES9.2,$)'
    endif

    line = ' '
    write(line, form) arg
    if(lg >= -3) then
      do while (line(len_trim(line):len_trim(line)) == '0') ! cut trailing zeros after three 0 after comma
        line = line(:len_trim(line)-1)
      enddo
    endif
    write(6,'(A,$)') line(:len_trim(line))
  end subroutine write_real

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
  !>
  !> @brief
  !> Checks the syntax of the input file.
  !>
  !> @details
  !> Does not check whether enough arguments are provided. (This is done in getkey.)
  !> The line length is restricted to 256 characters.
  !> Sections within sections are not possible.
  !> A leading '#' denotes a comment.
  !> Multiple lines can be separated with \
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine checkkeys(iunit)

  implicit none

  integer,           intent(in) :: iunit
  integer,           parameter  :: nkey = 9!put in here number of keys
  !dummy string with length 20: '                    '
  character(len=20),      parameter  :: key(nkey) = (/'additionalQs        ', 'PHONON              ', 'dimKptSet           ', &
                                                     &'dimQptSet           ', 'createKpoints       ', 'kPSetShift          ', &
                                                     &'K-POINTS            ', 'Pts2ChkCont         ', 'DEBUG               '/)
                                           !    &  'PHONON         ', 'K-POINTS       ', 'createKpoints  ', &
                                           !    & 'dimDenseKSet   ', 'dimDenseQSet   ',&
                                           !    &'duArrInt','dScalRea','duScalIn','dummyLog', &
                                           !    &  'dummyCha', 'dummyStr',&
                                           !    & 'kPSetShift     ','SECTION        ','PHON           ' /) !,'SELECT  ','PRECONST',
                                             ! 'ORBITALS','MAXIMIZE','FROZEN  ','DISENTGL' /)
  character(2048)               :: line
  character(256)                :: line1, keyword, section
  integer                       :: ios, ind, i, keyfirst
  logical                       :: insection, found, exist(nkey)

  exist     = .false.
  insection = .false.
  keyfirst  = 1
  rewind(iunit)

  do
    line1 = ' '
 1  read(iunit, '(A)', iostat=ios) line
    if(ios /= 0 .or. line == 'EXIT') exit
    ind = index(line, '#')         ! Remove
    if(ind /= 0) line(ind:) = ' ' ! comments

    ind = max(1, len_trim(line))                   !
    if(line(ind:ind) == '\') then                  !
      line1(len_trim(line1) + 1:) = line(:ind - 1) !
      goto 1                                       ! allow for multiple lines
    else                                           ! (continued with '\')
      line1(len_trim(line1) + 1:) = line           !
      line                        = line1          !
    endif                                          !

    if(line == ' ') cycle

    line    = adjustl(line)   ! Remove leading spaces
    ind     = index(line, ' ') ! Get keyword
    keyword = line(:ind - 1)    !
    line    = adjustl(line(ind + 1:))

    if(keyword == 'END') then
      if(.not.insection) then
        write (6, *) 'checkkeys: Missing SECTION statement before END statement.'
        NOstopNO
      endif
      insection = .false.
      keyfirst  = 1
      write(6,'(A)') 'END'
      cycle
    endif
    if(keyword == 'SECTION' .and. insection) then
      write (6, *) 'checkkeys: Missing END statement after SECTION statement.'
      NOstopNO
    endif

    ! Look for keyword in key-list (between two "SECTION" keys if insection=.true.)
    found = .false.
    i     = keyfirst
    do while(i <= nkey)
      if(trim(key(i)) == trim(keyword)) then
        if(keyword == 'SECTION') then
          i         = i + 1
          section   = line
          if(section == ' ') then
            write (6, *)  'checkkeys: Missing section keyword after SECTION.'
            NOstopNO
          endif
          if(trim(section) == trim(key(i))) then
            insection = .true.
            keyfirst  = i + 1
            found     = .true.
            write(6,'(A)') 'SECTION ' // trim(key(i))
            exit
          endif
        else
          if(exist(i)) then
            write (6, *) 'checkkeys: Keyword ' // trim(keyword) // ' given more than once.'
            NOstopNO
          endif
          found    = .true.
          exist(i) = .true.
          if(insection) write(6,'(''  '',$)')
          write(6,'(A)') keyword(:9) // trim(line)
          exit
        endif
      endif
      if(trim(key(i)) == 'SECTION') exit
      i = i + 1
    enddo

    if(.not.found) then
      if(keyword == 'SECTION') then
        write (6, *) 'checkkeys: Section keyword ' // trim(section) // ' unknown.'; NOstopNO
      else
        write(0,'(3A,$)') 'checkkeys: Keyword ', trim(keyword), ' unknown'
        if(insection) then
          write (6, *) ' (in section ' // trim(section) // ').'; NOstopNO
        else
          write (6, *) ' (outside sections).'; NOstopNO
        endif
      endif
    endif

    enddo

    if(insection) then
      write (6, *) 'checkkeys: Missing END statement after SECTION statement.'
      NOstopNO
    endif
    write(6,*)
  end subroutine checkkeys

end module m_jpParamInpF
