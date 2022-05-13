!----------------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!----------------------------------------------------------------------------------------------------------------------------------------
!
! module m_jpTestInput: Sophisticated tests for m_jpInit module
!
!> @author
!> Christian-Roman Gerhorst and Markus Betzinger
!>
!> @brief
!> Tests input data to be consistent.
!>
!> @details
!> This module contains routines performing more sophisticated tests than just writing data to logFile. Input data is used to calculate
!> known quantities or calculate quantities in different way to check data consistency.
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpTestInput.pdf'>document</a>.
!----------------------------------------------------------------------------------------------------------------------------------------
module m_jpTestInput

  use mod_juPhonUtils

  implicit none

  contains

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Main input-data testing routine calling other routines to perform the tests.
  !>
  !> @details
  !> This routine calls every routine performing a test related to data which is read in or parsed in m_jpInit module. Instead of calling
  !> every test routine from main program this routine summarizes the call of the routines for sake of readability in the main program.
  !>
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] stars      : Stars type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[out] dimens     : Dimension type, see types.f90.
  !> @param[out] lathar     : Lattice harmonics type, see types.f90.
  !> @param[out] sym        : Symmetries type, see types.f90.
  !> @param[out] input      : Input type, see types.f90.
  !> @param[out] kpts       : K-points type, see types.f90.
  !> @param[out] usdus      : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] Veff0      : Type containing unperturbed output potentials of Fleur, see types.f90.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !> @param[out] noPotsCon  : Number of points for which continuity test should be performed.
  !> @param[out] rbas1      : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2      : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] iloTable   : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] nRadFun    : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] gridf      : Initialized grid quantity to use intgrf routine.
  !> @param[out] kveclo     : Basis G-vectors of local orbitals.
  !> @param[out] nv         : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] GbasVec    : G-basis vectors
  !> @param[out] ilst       : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                          pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] ne         : Number of eigenvalues per k-point.
  !> @param[out] z          : Kohn-Sham eigenvectors.
  !> @param[out] rho0IR     : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur.
  !> @param[out] rho0MT     : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
  !> @param[out] mlh_atom   : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom  : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom  : Phase mediating between stars and plane waves.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine DetailedInitializationTests( atoms, stars, cell, dimens, lathar, sym, input, kpts, usdus, Veff0, testRadSolSw,        &
      & testKptsWeightSw, testCountValElecSw, testVeff0ContSw, testrho0ContSw, testBackRotMTCoordSysSw, testPsi0ContSw,            &
      & testOverlapSw, testXCintegrals, logUnit, noPtsCon, rbas1, rbas2, iloTable, nRadFun, gridf, kveclo, nv, GbasVec, ilst, ne, z, rho0IR, rho0MT,&
      & mlh_atom, nmem_atom, clnu_atom, vXC0IR, eXCIR, vXC0MT, eXCMT, ngdp, gdp )

    use m_types
    use m_juDFT_time, only : TimeStart, TimeNOstopNO

#include "cppmacro.h"

    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in) :: atoms
    type(t_stars),                 intent(in) :: stars
    type(t_cell),                  intent(in) :: cell
    type(t_dimension),             intent(in) :: dimens
    type(t_sphhar),                intent(in) :: lathar
    type(t_sym),                   intent(in) :: sym
    type(t_input),                 intent(in) :: input
    type(t_kpts),                  intent(in) :: kpts
    type(t_usdus),                 intent(in) :: usdus
    type(t_potential),             intent(in) :: Veff0

    ! Scalar parameters
    logical,                       intent(in) :: testRadSolSw
    logical,                       intent(in) :: testKptsWeightSw
    logical,                       intent(in) :: testCountValElecSw
    logical,                       intent(in) :: testVeff0ContSw
    logical,                       intent(in) :: testrho0ContSw
    logical,                       intent(in) :: testBackRotMTCoordSysSw
    logical,                       intent(in) :: testPsi0ContSw
    logical,                       intent(in) :: testOverlapSw
    logical,                       intent(in) :: testXCintegrals
    integer,                       intent(in) :: logUnit
    integer,                       intent(in) :: noPtsCon
    integer,                       intent(in) :: ngdp

    ! Array parameters
    real,                          intent(in) :: rbas1(:, :, 0:, :, :)
    real,                          intent(in) :: rbas2(:, :, 0:, :, :)
    integer,                       intent(in) :: iloTable(:, :, :)
    integer,                       intent(in) :: nRadFun(0:, :)
    real,                          intent(in) :: gridf(:, :)
    integer,                       intent(in) :: kveclo(:, :)
    integer,                       intent(in) :: nv(:, :)
    integer,                       intent(in) :: GbasVec(:, :)
    integer,                       intent(in) :: ilst(:, :, :)
    integer,                       intent(in) :: ne(:)
    MCOMPLEX,                      intent(in) :: z(:,:,:,:)
    complex,                       intent(in) :: rho0IR(:,:)
    real,                          intent(in) :: rho0MT(:,:,:,:)
    integer,                       intent(in) :: mlh_atom(:,0:,:)
    integer,                       intent(in) :: nmem_atom(0:, :)
    complex,                       intent(in) :: clnu_atom(:,0:,:)
    complex,                       intent(in) :: vXC0IR(:, :)
    complex,                       intent(in) :: eXCIR(:)
    real,                          intent(in) :: vXC0MT(:, :, :, :)
    real,                          intent(in) :: eXCMT(:, :, :)
    integer,                       intent(in) :: gdp(:, :)

    ! Local scalar variables
    integer                                   :: ispin
    logical                                   :: density_check
    integer                                   :: ikpt

    ! Local array parameters
    integer,           allocatable            :: nobd(:)


   !todo review this test, it should be activated seperately
    if ( testBackRotMTCoordSysSw ) then
      !Calculate radial component of the density for every l and m to compare it to a system where symmetry was artificially broken.
      !For instance, StTiO3, where each O has its own type
      call TimeStart( 'Test: Init: Rotated CS' )
      call CheckBackRotLHLCS ( rho0MT, atoms, input, lathar, clnu_atom, nmem_atom, mlh_atom )
      call TimeNOstopNO( 'Test: Init: Rotated CS' )
    end if

    if ( testRadSolSw .or. testKptsWeightSw .or. testCountValElecSw .or. testVeff0ContSw .or. testrho0ContSw .or. testPsi0ContSw &
      & .or. testOverlapSw .or. testXCintegrals ) then
      write(*, *)
      write(*, '(a)') 'Initiating input data test(s)...'
      write(*, '(a)') '--------------------------------'
      write(logUnit, *)
    write ( logUnit, * )
      write(logUnit, '(a)') 'Input data test(s)'
      write(logUnit, '(a)') '******************'
      write(logUnit, *)
    else
      write ( logUnit, * )
      write ( logUnit, * )
      write(*, *)
      write(*, '(a)') '----------------------------'
      write(*, '(a)') 'DISABLED input data test(s)!'
      write(*, '(a)') '----------------------------'
      write(logUnit, '(a)') 'DISABLED input data tests!'
      write(logUnit, '(a)') '**************************'
      return
    end if

    if ( testRadSolSw ) then
     ! Check rbas1 and rbas2
     write(*, '(2x,a)') 'Performing CheckOrthNormRadSol...'
     call TimeStart( 'Test: Init: radSolut' )
     call CheckOrthNormRadSol( atoms, sym, input, kpts, usdus, logUnit, rbas1, rbas2, nRadFun, gridf, kveclo, nv )
     call TimeNOstopNO( 'Test: Init: radSolut' )
     write ( logUnit, * )
    else
      write(*, '(2x,a)') 'DISABLED CheckOrthNormRadSol!'
      write ( logUnit, '(a)' ) "Orthogonality of <u|u'> and normalization of <u|u>..."
      write ( logUnit, '(a)' ) '-----------------------------------------------------'
      write (logUnit, '(a)')   '                                                    |__ DISABLED!'
    end if

    if ( testKptsWeightSw ) then
      ! Check weights in kpts file
      write(*, '(2x,a)') 'Performing CheckKptsWeights...'
      call TimeStart( 'Test: Init: kpts' )
      call CheckKptsWeights( kpts, logUnit )
      write ( logUnit, * )
      call TimeNOstopNO( 'Test: Init: kpts' )
    else
      write(*, '(2x,a)') 'DISABLED CheckKptsWeights!'
      write(logUnit,*)
      write ( logUnit, '(a)' ) 'Weights in kpts file'
      write ( logUnit, '(a)' ) '--------------------'
      write ( logUnit, '(a)' ) '                   |__DISABLED!'
    end if

    if ( testCountValElecSw ) then
      ! Check number of valence electrons
      write(*, '(2x,a)') 'Performing CheckValenceElectrons...'
      call TimeStart( 'Test: Init: #elec' )
      call CheckValenceElectrons( atoms, cell, kpts, dimens, input, logUnit )
      call TimeNOstopNO( 'Test: Init: #elec' )
      write ( logUnit, * )
    else
      write(*, '(2x,a)') 'DISABLED CheckValenceElectrons!'
      write(logUnit,*)
    write ( logUnit, '(a)' ) 'Number of valence electrons'
    write ( logUnit, '(a)' ) '---------------------------'
    write ( logUnit, '(a)' ) '                          |__DISABLED!'
    end if

   if ( testVeff0ContSw ) then
     write(*, '(2x,a)') 'Performing CheckDensnPot(potential)...'
     call TimeStart( 'Test: Init: potential' )
     density_check = .false.
     do ispin = 1, input%jspins
       !vpw can also be density expansion coefficients or wavefunctions expansion coefficients
       call CheckDensnPot( atoms, sym, stars, cell, lathar, noPtsCon, logUnit, density_check, ispin, Veff0%vpw_uw, Veff0%vr )
     end do
     call TimeNOstopNO( 'Test: Init: potential' )
    else
      write(*, '(2x,a)') 'DISABLED CheckDensnPot(potential)!'
      write(logUnit,*)
      write ( logUnit, '(a)' ) 'Check potential for continuity at MT boundary:'
      write(logUnit,'(a)')     '----------------------------------------------'
      write(logUnit,'(a)') '                                             |_DISABLED'
   end if

   if ( testrho0ContSw ) then
     write(*, '(2x,a)') 'Performing CheckDensnPot(density)...'
     call TimeStart( 'Test: Init: density' )
     density_check = .true. !density is checked
     do ispin = 1, input%jspins
       call CheckDensnPot( atoms, sym, stars, cell, lathar, noPtsCon, logUnit, density_check, ispin, rho0IR, rho0MT )
     end do
     call TimeNOstopNO( 'Test: Init: density' )
    else
      write(*, '(2x,a)') 'DISABLED CheckDensnPot(density)!'
      write(logUnit, *)
      write ( logUnit, '(a)' ) 'Check density for continuity at MT boundary:'
      write(logUnit,'(a)')     '--------------------------------------------'
      write(logUnit,'(a)')     '                                           |_DISABLED'
   end if


   ! Check the continuity of the wavefunctions at the MT boundary at a given number of points.
   ! The wavefunction test costs time therefore it can be switched off
   if ( testPsi0ContSw ) then
     write(*, '(2x,a)') 'Performing CheckWaveFunctions...'
     call TimeStart( 'Test: Init: Wavefunctions' )
     allocate( nobd(kpts%nkpt) )
     nobd = ne

     call CheckWaveFunctions( atoms, cell, kpts, sym, dimens, usdus, noPtsCon, GbasVec, ilst, nv, ne, z, nRadFun, kveclo, rbas1, rbas2, &
       & iloTable, logUnit )
     call TimeNOstopNO( 'Test: Init: Wavefunctions' )
   else
     write(*, '(2x,a)') 'DISABLED CheckWaveFunctions!'
     write(logUnit, *)
     write(logUnit,'(a)') 'Check wave functions for continuity at MT boundary:'
     write(logUnit,'(a)') '---------------------------------------------------'
     write(logUnit,'(a)') '                                                  |__DISABLED!'
   end if
   write(logUnit, *)

   ! Test where the basis functions and the expansion coefficients for the wavefunctions are tested calculating the overlap matrix
   if ( testOverlapSw ) then
     write(*, '(2x,a)') 'Performing OverlapMatrix...'
     write (logUnit, '(a)')   'Orthogonality and normalization test of Kohn--Sham wavefunctions!'
     write ( logUnit, '(a)' ) '-----------------------------------------------------------------'
     do ispin = 1, input%jspins
       do ikpt = 1, kpts%nkpt
         call OverlapMatrix( atoms, cell, input, sym, dimens, kpts, usdus, GbasVec, ilst, rbas1, rbas2, ikpt, nv, kveclo, ne,&
           &z(:, :, ikpt, 1), ispin, gridf, iloTable, nRadFun, logUnit )
       end do
     end do
     write (logUnit, '(a)')   '                                                                |__ passed!'
    else
      write(*, '(2x,a)') 'DISABLED OverlapMatrix!'
     write (logUnit, '(a)')   'Orthogonality and normalization test of Kohn--Sham wavefunctions!'
     write ( logUnit, '(a)' ) '-----------------------------------------------------------------'
     write (logUnit, '(a)')   '                                                                |__ DISABLED!'
   end if

   if ( testXCintegrals ) then
     write(*, '(2x,a)') 'Performing CalcXCintegrals...'
     call CalcXCintegrals( atoms, cell, stars, lathar, ngdp, rho0IR, vXC0IR, vXC0MT, eXCIR, eXCMT, rho0MT, gdp, mlh_atom, nmem_atom, clnu_atom, logUnit )
   else
     write(*, '(2x,a)') 'DISABLED CalcXCintegrals!'
     write(logUnit,'(a)')
     write(logUnit,'(a)') 'Calculate IR and MT integrals with integrand rho V_xc (rhoVxc) and integrand rho E_xc (rhoExc):'
     write(logUnit,'(a)') '-----------------------------------------------------------------------------------------------'
     write(logUnit,'(a)') '                                                                                              |__DISABLED'
   end if

  end subroutine DetailedInitializationTests

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Check the orthonormality of the radial solutions.
  !>
  !> @details
  !> Calculates overlaps of <u|u> = 1, <u|udot> = 0, <u_LO|u_LO> with same ( = 1) and different l ( = 0 ) at same atom
  !>
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] sym        : Symmetries type, see types.f90.
  !> @param[out] input      : Input type, see types.f90.
  !> @param[out] kpts       : K-points type, see types.f90.
  !> @param[out] usdus      : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !> @param[out] rbas1      : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2      : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] nRadFun    : Number of radial functions per orbital quantum number l and atom type.
  !> @param[out] gridf      : Initialized grid quantity to use intgrf routine.
  !> @param[out] kveclo     : Basis G-vectors of local orbitals.
  !> @param[out] nv         : Number of LAPW G-basis vectors for given k-point.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CheckOrthNormRadSol( atoms, sym, input, kpts, usdus, logUnit, rbas1, rbas2, nRadFun, gridf, kveclo, nv )

     
    use m_types, only : t_atoms, t_input, t_kpts, t_usdus, t_sym
    use m_setabc1locdn1
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameters
    type(t_atoms),             intent(in) :: atoms
    type(t_sym),               intent(in) :: sym
    type(t_input),             intent(in) :: input
    type(t_kpts),              intent(in) :: kpts
    type(t_usdus),             intent(in) :: usdus

    ! Scalar parameter
    integer,                   intent(in) :: logUnit

    ! Array parameters
    real,                      intent(in) :: rbas1(:, :, 0:, :, :)
    real,                      intent(in) :: rbas2(:, :, 0:, :, :)
    integer,                   intent(in) :: nRadFun(0:, :)
    real,                      intent(in) :: gridf(:, :)
    integer,                   intent(in) :: kveclo(:, :)
    integer,                   intent(in) :: nv(:, :)

    ! Scalar variables
    integer                               :: itype
    real                                  :: overlap
    real                                  :: expec_overlap
    integer                               :: ispin
    integer                               :: oqn_l
    integer                               :: p
    integer                               :: ilo
    integer                               :: mqn_m
    integer                               :: nvmax
    integer                               :: nmat
    integer                               :: ikpt
    integer                               :: imesh

    ! Array variables
    real,          allocatable            :: intgrnd(:)
!    complex, allocatable                  :: matchDum(:, :, :)
!    complex, allocatable                  :: bascof_lo(:, :, :, :, :, :, :)
    logical,       allocatable            :: enough(:)
    integer,       allocatable            :: nkvec(:, :)
    integer,       allocatable            :: kvec(:, :, :)
    integer,       allocatable            :: nbasf0(:, :)
    real,          allocatable            :: alo1(:, :)
    real,          allocatable            :: blo1(:, :)
    real,          allocatable            :: clo1(:, :)


    write ( logUnit, '(a)' ) "Orthogonality of <u|u'> and normalization of <u|u>..."
    write ( logUnit, '(a)' ) '-----------------------------------------------------'

    allocate( intgrnd(atoms%jmtd) )
    do ispin = 1, input%jspins
      do itype = 1, atoms%ntype
        do oqn_l = 0, atoms%lmax( itype )
          ! Test of normalization of <u|u> = 1
          intgrnd = 0.0
          do imesh = 1, atoms%jri(itype)
            intgrnd(imesh) = rbas1(imesh, 1, oqn_l, itype, ispin) * rbas1(imesh, 1, oqn_l, itype, ispin) + rbas2(imesh, 1, oqn_l, itype, ispin) &
            &* rbas2(imesh, 1, oqn_l, itype, ispin)
          end do
          overlap = Intgrf( intgrnd, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf )
          if  ( abs( overlap - 1 ) > 1e-8 ) then
            write ( logUnit, '(a,i5)' ) 'Atom type, orbital quantum number: ', itype, oqn_l
            write ( logUnit, '(a,f15.8)' ) 'Error for <u|u>:  expected = 1, calculated = ', overlap
            call JuDFT_warn( 'Wrong norm of <u|u> detected', calledby='TestOrthNormCondRadSol', &
              &hint='Check logfile for more information!', file='jpTestInput_mod.F90' )
          end if

          ! Test of orthogonality of <u|u'> = 0
          intgrnd = 0.0
          do imesh = 1, atoms%jri(itype)
            intgrnd(imesh) = rbas1(imesh, 1, oqn_l, itype, ispin) * rbas1(imesh, 2, oqn_l, itype, ispin) + rbas2(imesh, 1, oqn_l, itype, ispin) &
            &* rbas2(imesh, 2, oqn_l, itype, ispin) !cmplx conj! are rbas1 and rbas2 complex when no inversion symmetry
          end do
          overlap = Intgrf( intgrnd, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf )
          if  ( abs( overlap - 0 ) > 1e-8 ) then
            write ( logUnit, '(a,f15.8)' ) 'Atom type, orbital quantum number: ', itype, oqn_l
            write ( logUnit, * ) "Error for <u|u'>. expected = 0, calculated = ", overlap
            call JuDFT_warn('Unorthogonal radial solution detected', calledby='TestOrthNormCondRadSol', &
              &hint='Check logfile for more information!', file='jpTestInput_mod.F90')
          end if
        end do ! oqn_l
      end do ! itype
    end do ! ispin

    ! Check normalization of <u_LO|u_LO>
    if ( any( nRadFun > 2 ) ) then
      allocate( enough(atoms%nat), alo1(atoms%nlod, atoms%ntype), blo1(atoms%nlod, atoms%ntype), clo1(atoms%nlod, atoms%ntype), &
        & nkvec(atoms%nlod, atoms%nat), kvec(2 * (2* atoms%llod + 1), atoms%nlod, atoms%nat), nbasf0(atoms%nlod, atoms%nat) )
      do ikpt =  1, kpts%nkpt
        do ispin = 1, input%jspins
          nvmax = nv(ispin, ikpt)
          nmat = nvmax + atoms%nlotot
          call Setabc1locdn1( atoms%ntype, atoms%nat, atoms%nlod, atoms%llod, atoms%lmaxd, nvmax, nmat, atoms%ntype, atoms%neq, &
            & atoms%nlo, atoms%llo, atoms%l_dulo, atoms%invsat, sym%invsatnr, usdus%ddn, usdus%us, usdus%dus, usdus%uds, usdus%duds, &
            &usdus%dulos, usdus%ulos, usdus%dulon, usdus%uulon, atoms%nlotot, kveclo, enough, nkvec, kvec, nbasf0, alo1, blo1, clo1 )
          do itype = 1, atoms%ntype
            do oqn_l = 0, atoms%llod
              do p = 3, nRadFun(oqn_l, itype)
                intgrnd = rbas1(:, p, oqn_l, itype, ispin) * rbas1(:, p, oqn_l, itype, ispin) &
                &+ rbas2(:, p, oqn_l, itype, ispin) * rbas2(:, p, oqn_l, itype, ispin)
                overlap = Intgrf( intgrnd, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf )
                if ( abs( overlap - 1 ) > 1e-8 ) then
                  write ( logUnit, '(a,4(i3))' ) 'Atom type, LOindex, orbital quantum number, k-point: ', itype, ilo, oqn_l, ikpt
                  write ( logUnit, '(a,f15.8)' ) "Error for <u|u'>. expected = 0, calculated = ", overlap
                  call JuDFT_warn('Wrong norm of <u_LO|u_LO>', calledby='TestOrthNormCondRadSol', &
                    &hint='Check logfile for more information!', file='jpTestInput_mod.F90')
                end if
              end do ! p
            end do ! oqn_l
          end do ! itype
        end do ! ispin
      end do ! ikpt
    endif ! if there are LOs
    write (logUnit, '(a)') '                                                    |__ passed!'

  end subroutine CheckOrthNormRadSol

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Check consistency of kpts-weights.
  !>
  !> @details
  !> All kpts-weights summed up should be 1. This routines checks it.
  !>
  !> @param[out] kpts       : K-points type, see types.f90.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CheckKptsWeights( kpts, logUnit )

    use m_types
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type variable
    type(t_kpts), intent(in) :: kpts

    ! Scalar variable
    integer,      intent(in) :: logUnit

    ! Local scalar
    real                     :: sumOfKpts

    sumOfKpts = sum( kpts%wtkpt )

    write ( logUnit, '(a)' ) 'Weights in kpts file'
    write ( logUnit, '(a)' ) '--------------------'
    if ( sumOfKpts - 1.0 > 9e-12 ) then
      write ( logUnit, '(a)' ) 'Inconsistency detected:'
      write ( logUnit, '(a)' ) '  expected sum = ', 1, '; calculated sum = ', sumOfKpts
      call juDFT_warn( 'Weights in kpts file incorrect!', calledby='checkKptsWeights', hint='Check logfile for more information!', &
        &file='jpTestInput_mod.F90' )
    end if
    write ( logUnit, '(a)' ) '                    |__passed!'

  end subroutine CheckKptsWeights

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Checks if stored number of electrons is consistent.
  !>
  !> @details
  !> Stored number of electrons in zelec is compared to the number of electrons the routine Fermie calculates from the eig file.
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[out] kpts       : K-points type, see types.f90.
  !> @param[out] dimens     : Dimension type, see types.f90.
  !> @param[out] input      : Input type, see types.f90.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CheckValenceElectrons( atoms, cell, kpts, dimens, input, logUnit )

    use m_fermie
    use m_types
    use m_JPConstants, only: tpi

    implicit none

    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: kpts
    type(t_dimension),              intent(in) :: dimens
    type(t_input),                  intent(in) :: input
    integer,                        intent(in) :: logUnit

!   .. Array Arguments ..

    real                                       :: ef
    real                                       :: seigscv
    real                                       :: ts
    real                                       :: zelec_calc

    real,              allocatable             :: w (:, :, :)
    real                                       :: zelec_dummy(1)
    integer                                    :: irecl
    integer                                    :: iunit =  203
    integer                                    :: idum
    integer                                    :: int_len
    integer                                    :: nkpt_dummy(1)

    type(t_Jij)                                :: jij_type
    type(t_noco)                               :: noco_type

    write ( logUnit, '(a)' ) 'Number of valence electrons'
    write ( logUnit, '(a)' ) '---------------------------'
    allocate( w(dimens%neigd, kpts%nkptf, dimens%jspd) )
    allocate( jij_type%eig_l(dimens%neigd +  6, jij_type%nkpt_l) )

    zelec_dummy(1) = input%zelec
    nkpt_dummy(1) = kpts%nkpt
    ! calculate record length of the eig-file
    irecl = atoms%ntype * ( atoms%lmaxd + 1 + atoms%nlod ) + 2
#ifdef CPP_INVERSION
    irecl = 8 * ( irecl + 4 + dimens%neigd * ( dimens%nbasfcn + 1 ) )
#else
    irecl = 8 * ( irecl + 4 + dimens%neigd * ( 2 * dimens%nbasfcn + 1 ) )
#endif
    irecl = irecl + 4 * atoms%nat * atoms%nlod * ( 2 * atoms%llod + 1 )

    idum = size(transfer(tpi,(/1/)))
    if ( idum == 1 ) then
      int_len = 8                 ! idum = 1 if real*8
    else if ( idum == 2 ) then
      int_len = 4                 ! idum = 2 if real*4
    else
      NOstopNO'read_eig: strange length of integers'
    endif
    irecl = irecl + int_len * ( 3 + 3 * dimens%nvd )

    call Fopen( iunit, name='eig', status='old', action='read', form='unformatted', access='direct', recl=irecl )
    !TODO changes in fermie with external subroutine!!!!!!
    call Fermie( dimens%neigd, kpts%nkptf, dimens%jspd, 1, atoms%ntype, atoms%lmaxd, atoms%nlod, input%jspins, 1, atoms%ntype, &
      & nkpt_dummy, iunit, irecl, .false., .false., input%gauss, input%delgau, zelec_dummy, input%tkb, input%tria, noco_type%l_noco, &
      & noco_type%l_ss, ef, seigscv, ts, w, noco_type%qss, jij_type%l_J, jij_type%l_disp, cell%bmat, jij_type%nkpt_l, jij_type%eig_l, &
      & 0, 1)
    call Fclose(iunit)

    zelec_calc =  sum( w ) * 2 ! spin factor
    if ( abs(input%zelec - zelec_calc) > 1e-7) then
      write ( logUnit, '(a)' ) 'Inconsistency in number of valence electrons!'
      write ( logUnit, '(a,f15.8,a,f15.8)' ) 'Expected electrons = ', input%zelec, 'Calculated electrons = ', zelec_calc
    end if
    write ( logUnit, '(a)' ) '                          |__passed!'
  end subroutine CheckValenceElectrons

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Check continuity of density or potential at MT boundary.
  !>
  !> @details
  !> The density (if density=true) or potential(if density = false) at the MT boundary for noPtsCon random points is checked.
  !>
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] sym        : Symmetries type, see types.f90.
  !> @param[out] stars      : Stars type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[out] lathar     : Lattice harmonics type, see types.f90.
  !> @param[out] noPotsCon  : Number of points for which continuity test should be performed.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CheckDensnPot( atoms, sym, stars, cell, lathar, noPtsCon, logUnit, density, ispin, fpw, fr )

    use m_cotra
    use m_starf
    use m_ylm_old
    use mod_juPhonUtils
    use m_types

    implicit none

    ! Type parameters
    type(t_atoms),        intent(in)  :: atoms
    type(t_sym),          intent(in)  :: sym
    type(t_stars),        intent(in)  :: stars
    type(t_cell),         intent(in)  :: cell
    type(t_sphhar),       intent(in)  :: lathar

    ! Scalar parameters
    integer,              intent(in)  :: noPtsCon !number of points to check continuity
    integer                           :: logUnit
    logical,              intent(in)  :: density ! is function density?
    integer,              intent(in)  :: ispin

    ! Array parameters
    complex,              intent(in)  :: fpw(:, :)
    real,                 intent(in)  :: fr(:, 0:, :, :)


    ! Local scalars
    integer                           :: irandPt,ieq !loop variable
    real                              :: euclidNorm !variable to store eculidian L2 norm of vector
    integer                           :: itype !certain atom type
    integer                           :: icoord !loop variable
    integer                           :: imatmulc, imatmulr !loop variable
    integer                           :: istar !loop variable
    integer                           :: iatom
    integer                           :: symAt_temp
    integer                           :: symOpr_temp
    integer                           :: ilath
    real                              :: linCombBas
    integer                           :: l_temp
    integer                           :: imem
    integer                           :: lm_temp
    real                              :: av, dms, rms

    ! Local arrays
    real,   allocatable               :: randPtsCart(:, :, :) !cartesian dimension, number of Pts, itype, random points
    real,   allocatable               :: randPtsGrid(:, :, :) !random points on grid, coordinates, number of points
    real,   allocatable               :: randPtsMTLoc(:, :, :) !random points in local atomic coordinate (coordinates, number of points, ntype)
    complex, allocatable              :: starAtRanPt(:) !value of star at certain point ! ng3 = nq3
    real,   allocatable               :: sumfuncvalI(:,:) !sum of function values coming from interstitial
    real                              :: rvec(3), rvec_int(3), rotatedVector(3) !rotated vector after symmetry operation
    complex, allocatable              :: ylm(:) !TODO eigentlich lmax(itype)
    real,    allocatable              :: sumOfMTFuncVal(:,:)

    allocate( randPtsCart(3, noPtsCon, atoms%nat), randPtsGrid(3, noPtsCon, atoms%ntype), randPtsMTLoc(3, noPtsCon, atoms%ntype) )
    allocate( starAtRanPt(stars%nq3), ylm( ( atoms%lmaxd + 1 )**2 ) )

    allocate( sumfuncvalI(noPtsCon,atoms%nat), sumOfMTFuncVal(noPtsCon, atoms%nat) )
    !Generate noPtsCon random points on the MT sphere
    randPtsCart(:, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieq = 1, atoms%neq(itype)
        iatom = iatom + 1
        call sphpts(randPtsCart(:, :, iatom), noPtsCon, atoms%rmt(itype), atoms%pos(:, iatom))
      end do
    end do

    !Evaluate IR density or potential for every star for each random point at each MT sphere boundary in the unit cell
    sumfuncValI = 0.0
    iatom       = 0
    do itype = 1, atoms%ntype
      do ieq = 1, atoms%neq(itype)
        iatom = iatom + 1
        do irandPt = 1, noPtsCon
          rvec = randPtsCart(:, irandPt, iatom)
          ! transform rvec in internal, real-space coordinates --> rvec_int
          call cotra1(rvec, rvec_int, cell%bmat)
          ! Evaluate stars at random point
          call Starf3(sym%nop, stars%nq3, sym%symor, stars%kv3, sym%mrot, sym%tau, rvec_int, sym%invtab, starAtRanPt)
          sumfuncValI(irandPt, iatom) = 0
          do istar = 1, stars%nq3 !nq3 is number of stars
            sumfuncValI(irandPt, iatom) = sumfuncvalI(irandPt, iatom) + real(fpw(istar, ispin) * starAtRanPt(istar)) * stars%nstr(istar)
            !jsp noch einfügen, warum Multiplikation with G-vectors
          end do  !istar
        end do  !irandPt
      end do  !ieq
    end do  !itype


    iatom = 0; sumOfMTFuncVal = 0
    do itype = 1, atoms%ntype
      !Evaluate MT density or potential for each lattice harmonic at each random point
      randPtsGrid = 0.0
      do ieq = 1, atoms%neq(itype)
        iatom       = iatom + 1
        symAt_temp  = atoms%ntypsy(iatom) ! symmetry of atom iatom
        symOpr_temp = atoms%ngopr (iatom) ! symmetry operation mapping local to global coordinate system
        do irandPt = 1, noPtsCon
          rvec = randPtsCart(:, irandPt, iatom) - atoms%pos(:, iatom)
          if ( symOpr_temp /= 1 ) then
            !Transform random point into internal real space coordinates
            call cotra1( rvec, rvec_int, cell%bmat ) !allocate and deallocate randPtsGrid
            do imatmulr = 1, 3
              rotatedVector(imatmulr) = 0.0
              do imatmulc = 1, 3
                rotatedVector(imatmulr) = rotatedVector(imatmulr) + sym%mrot(imatmulr, imatmulc, symOpr_temp) *rvec_int(imatmulc)
              end do
            end do
            call cotra0( rotatedVector, rvec, cell%amat )
          end if
          call Ylm4( atoms%lmax(itype), rvec, ylm )

          ! Evaluate for this random point for every lattice harmonic
          do ilath = 0, lathar%nlh(symAt_temp)
            linCombBas = 0.0
            l_temp = lathar%llh(ilath, symAt_temp) * (lathar%llh(ilath, symAt_temp) + 1) + 1 !understand this
            do imem = 1, lathar%nmem(ilath, symAt_temp)
              lm_temp = l_temp +  lathar%mlh(imem, ilath, symAt_temp)
              linCombBas = linCombBas +  real(lathar%clnu(imem, ilath, symAt_temp) * ylm(lm_temp))
            end do
              sumOfMTFuncVal(irandPt, iatom) = sumOfMTFuncVal(irandPt, iatom) + fr(atoms%jri(itype), ilath, itype, ispin) * linCombBas
          end do
          sumOfMTFuncVal(irandPt, iatom) = sumOfMTFuncVal(irandPt, iatom)
        end do
      end do
    end do

    ! Logfile output
    if( density) then
      write(logUnit,*)
      write ( logUnit, * ) 'Check density for continuity at MT boundary:'
      write(logUnit,*)     '--------------------------------------------'
    else
      write(logUnit,*)
      write ( logUnit, * ) 'Check potential for continuity at MT boundary:'
      write(logUnit,*)     '----------------------------------------------'
    end if
    iatom = 0
    do itype = 1, atoms%ntype
      do ieq = 1, atoms%neq(itype)
        iatom = iatom + 1
        write ( logUnit, * ) '  Atom:', iatom
        write ( logUnit, * ) '  coordinate(cartesian)          IR          MT'
        do irandPt = 1,noPtsCon
          write( logUnit, '(3f8.4,3x,2f12.8)' ) randPtsCart(:, irandPt, itype), sumfuncValI(irandPt, iatom),&
            &sumOfMTFuncVal(irandPt, iatom)
        end do
        call Fitchk(sumfuncvalI(:, iatom), sumOfMTFuncVal(:,iatom), noPtsCon, av, rms, dms)
        write(logUnit,'(/3x,A,f10.5)'   ) 'average value interstitial:', av
        write(logUnit,'(3x,A,2f10.5,A/)') 'rms, dms:          ', rms, dms, ' in %'
      end do
    end do
  end subroutine CheckDensnPot

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
  subroutine checkWaveFunctions( atoms, cell, kpts, sym, dimens, usdus, noPtsCon, GbasVec, ilst, nv, ne, z, nRadFun, kveclo, rbas1, &
    & rbas2, iloTable, logUnit )

#include "cppmacro.h"

    use m_JPConstants
    use m_types, only : t_noco, t_atoms, t_cell, t_kpts, t_sym, t_dimension, t_usdus
    use m_abcof
    use m_cotra
    use m_ylm_old
     
    use m_juDFT_time, only : TimeStart, TimeNOstopNO

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in) :: atoms
    type(t_cell),              intent(in) :: cell
    type(t_kpts),              intent(in) :: kpts
    type(t_sym),               intent(in) :: sym
    type(t_dimension),         intent(in) :: dimens
    type(t_usdus),             intent(in) :: usdus

    ! Scalar parameter
    integer,                   intent(in) :: noPtsCon !number of points to check continuity
    integer,                   intent(in) :: logUnit

    ! Array parameter
    integer,                   intent(in) :: GbasVec(:, :)
    integer,                   intent(in) :: ilst(:, :, :)
    integer,                   intent(in) :: nv(:, :)
    integer,                   intent(in) :: ne(:)
    MCOMPLEX,                  intent(in) :: z(:,:,:,:)
    integer,                   intent(in) :: nRadFun(0:, :)
    integer,                   intent(in) :: kveclo(:, :)
    real,                      intent(in) :: rbas1(:,:,0:,:,:)
    real,                      intent(in) :: rbas2(:,:,0:,:,:)
    integer,                   intent(in) :: iloTable(:, 0:,: )

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
    integer                               :: ikpt
    integer                               :: nobd
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
    real                                  :: bkpt(3)
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
    real,        allocatable              :: mismatchAvg(:)
    real                                  :: currPnt(3)
    integer, allocatable :: ngoprI(:)


    allocate( noco_type%alph(atoms%ntype), noco_type%beta(atoms%ntype) ) !Up to now those variables are only of dummy character
    allocate( randPtsCart(3, noPtsCon, atoms%ntype) )
    allocate( mismatch(atoms%nat), mismatchAvg(atoms%nat) )
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

    write(logUnit,'(a)') 'Check wave functions for continuity at MT boundary:'
    write(logUnit,'(a)') '---------------------------------------------------'
    write(logUnit,*)
    write(logUnit,'(a)') '   jspin    ikpt    band-averaged abs. discontinuity per atom (index represented by column)'
    ngoprI(:) = 1
    do jspin =  1, dimens%jspd
      mismatchAvg(:) = 0
      do ikpt = 1, kpts%nkpt

        mismatch = 0

        bkpt = kpts%bk(:,ikpt)
        nobd = ne(ikpt)
        nmat = nv(jspin, ikpt) + atoms%nlotot !dimension of Hamiltonian

        allocate(acof(nobd, 0:dimens%lmd, atoms%nat), bcof(nobd, 0:dimens%lmd, atoms%nat), &
          &ccof(-atoms%llod:atoms%llod, nobd, atoms%nlod, atoms%nat))

        call Abcof( atoms%lmaxd, atoms%ntype, dimens%neigd, nobd, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, dimens%lmd, &
          &dimens%nbasfcn, atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, &
          &atoms%neq, atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, bkpt, GbasVec(1, ilst(:nv(jspin, ikpt), ikpt, jspin)), &
          &GbasVec(2, ilst(:nv(jspin, ikpt), ikpt, jspin)), GbasVec(3, ilst(: nv(jspin, ikpt), ikpt, jspin)), nv(:, ikpt), nmat, &
          &ne(ikpt), z(:, :, ikpt, jspin), usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, &
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
                call cotra3(kpts%bk(:,ikpt)+GbasVec(:, ilst(iG, ikpt, jspin)),Gext_temp,cell%bmat)
                expFunc = exp( iu *  dot_product( Gext_temp, currPnt ) )

                do iband = 1,ne(ikpt)
                  ! this is no linear run in storage but instead of ne * nv matrix vector operations and dot products we only do nv
                  ! sum up all contributions of the wavefunctions looping over the IR basis functions
                  sumInterst = sumInterst + z(iG, iband, ikpt, jspin)  * expFunc
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
                    do iband = 1, ne(ikpt)
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
            mismatch(iatom) = mismatch(iatom) + abs( sumOfMTFuncVal - sumInterst ) / sum( ne(:kpts%nkpt) ) / noPtsCon
            deallocate(coeffMT)
          end do ! ieqat
        deallocate(ylm)
      end do ! itype
        deallocate(acof, bcof, ccof)

        !todo output should be reviewed
        write(logUnit,'(2(4x,i4),4x,40es15.8)') jspin,ikpt,mismatch
        mismatchAvg(:) = mismatchAvg(:) + mismatch(:)
      end do  !ikpt
      mismatchAvg(:) = mismatchAvg(:) / kpts%nkpt
      write(logUnit,*)
      write(logUnit,'(a,2x,40es15.8)') 'average mismatch: ', mismatch
      write(logUnit,*)
    end do  !jspin

  end subroutine checkWaveFunctions

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Orthonormality of the wavefunctions is tested.
  !>
  !> @details
  !> The analytical overlap is calculated and multiplied with the zs for a given k-point. The resulting overlapping matrix of the wave
  !> functions should be equals to the unit matrix. This is checked within this routine. See also a picture in the smartphone.
  !>
  !> @note this has not been tested for systems with inversion symmetry.
  !>
  !> @todo insert picture from smartphone
  !>
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[out] input      : Input type, see types.f90.
  !> @param[out] sym        : Symmetries type, see types.f90.
  !> @param[out] dimens     : Dimension type, see types.f90.
  !> @param[out] kpts       : K-points type, see types.f90.
  !> @param[out] usdus      : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[out] GbasVec    : G-basis vectors
  !> @param[out] ilst       : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon. This
  !>                          pointer array contains the right index for GbasVec array to "unfold" G-basis vectors again.
  !> @param[out] rbas1      : Large components of radial solution, its energy derivative and u_LO
  !> @param[out] rbas2      : Small components of radial solution, its energy derivative and u_LO
  !> @param[out] nv         : Number of LAPW G-basis vectors for given k-point.
  !> @param[out] kveclo     : Basis G-vectors of local orbitals.
  !> @param[out] ne         : Number of eigenvalues per k-point.
  !> @param[out] z          : Kohn-Sham eigenvectors.
  !> @param[out] gridf      : Initialized grid quantity to use intgrf routine.
  !> @param[out] iloTable   : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] nRadFun    : Number of radial functions per orbital quantum number l and atom type.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine overlapMatrix( atoms, cell, input, sym, dimens, kpts, usdus,  GbasVec, ilst, rbas1, rbas2, ikpt, nv, kveclo, ne, z, ispin, &
    & gridf, iloTable, nRadFun, logUnit )

    use m_JPConstants, only: fpi, iu, tpi
    use m_intgr
     
    use m_abcof3
    use m_cotra
    use m_types, only : t_atoms, t_cell, t_input, t_sym, t_dimension, t_usdus, t_kpts
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type variables
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_sym),                    intent(in) :: sym
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_usdus),                  intent(in) :: usdus

    ! Scalar variables
    integer,                        intent(in) :: ikpt
    integer,                        intent(in) :: ispin
    integer,                        intent(in) :: logUnit

    ! Array variables
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: kveclo(:, :)
    integer,                        intent(in) :: ne(:)
    MCOMPLEX,                       intent(in) :: z(:,:)
    real,                           intent(in) :: gridf(:, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)

    ! Local type variables
    type (od_inp)                              :: odi_type !Both types, this and the one beneath have to be taken from od_types module
    type (od_sym)                              :: ods_type

    ! Local scalar variables
    integer                                    :: iatom,iatom0
    integer                                    :: itype,iband1,iband2
    integer                                    :: ieqat
    integer                                    :: iG,ng
    integer                                    :: iGprime
    complex                                    :: contrI
    real                                       :: Gdiff_norm, rdum
    complex                                    :: contrMT
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: p1, p2
    real                                       :: argument_temp
    integer                                    :: overlapIndex
    complex                                    :: cdum
    integer                                    :: lm
    integer                                    :: nmat
    integer                                    :: logFileUnit
    integer                                    :: ioqn_l
    integer                                    :: p
    integer                                    :: Glocmax
    integer                                    :: jspin
    integer                                    :: ilo
    integer                                    :: ieq
    logical                                    :: passed
    integer                                    :: l
    integer                                    :: nGLOs
    integer                                    :: iGloc
    integer                                    :: jatom
    integer                                    :: iGlocprime
    integer                                    :: maxlmp
    integer                                    :: locGindex
    integer                                    :: idump
    integer                                    :: jdump
    integer                                    :: pprime
    integer                                    :: loGshift
    integer                                    :: loGprimeShift
    integer                                    :: jlo
    integer                                    :: l_jlo
    integer                                    :: nGLO_jlo
    integer                                    :: mqn_m_jlo
    integer                                    :: idum
    integer                                    :: jdum
    complex                                    :: overlap

    ! Local array variables
    integer                                    :: Gdiff_int(3)
    real                                       :: Gdiff_ext(3)
    complex,           allocatable             :: overlapMat(:)
    complex,           allocatable             :: smat(:,:)
    real,              allocatable             :: f(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,           allocatable             :: matchingCoeff(:, :, :, :)
    complex,           allocatable             :: carr(:)
    real,              allocatable             :: overlapRadFun(:, :, :, :, :)
    complex,           allocatable             :: overlapGGLO(:, :)
    complex,           allocatable             :: overlapGLOGLO(:)
    integer,           allocatable             :: iPtable(:, :)
    integer :: loShift
    integer :: atomShift
    integer :: atomShiftprev
    integer, allocatable :: ngoprI(:)

    allocate(ngoprI(atoms%nat))
    allocate( overlapMat(nv(ispin, ikpt) * ( nv(ispin, ikpt) + 1 ) /  2) )
    ! lmp is index encoding oqn_l, mqn_m and p, this calculates the max value
!       write(*, *) 'nbasfcn', dimens%nbasfcn, nv(ispin, ikpt), atoms%nlotot, dimens%nvd
    maxlmp = maxval( (/ (sum( (/ (( 2 * l + 1 ) * nRadFun(l, itype), l = 0, atoms%lmax(itype)) /) ), itype=1, atoms%ntype) /) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), bascof_lo(3, -atoms%llod:atoms%llod, &
      &4 * atoms%llod + 2, atoms%nlod, atoms%nat))
    allocate( matchingCoeff(dimens%nvd,2, 0:dimens%lmd, atoms%nat) )
    allocate( f(atoms%jmtd) )
    allocate( iPtable(atoms%nlod, atoms%ntype) )

    nmat = nv(ispin, ikpt) + atoms%nlotot !dimension of Hamiltonian
    allocate( smat(nmat, nmat) )
 !   if (ikpt ==64 ) then
 !     write (962, *) nv(ispin, 28)
 !     write (962, *) atoms%nlotot
 !   end if

    atomShiftprev = 0
    !we do not want rotated local systems
    ngoprI(:) = 1

      ! Determine the matching coefficients for the basis functions, i.e. these are not the matching coefficients acof, bcof, ccof!
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, ispin, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, ilst(:nv(ispin, ikpt), ikpt, ispin)), &
        & GbasVec(2, ilst(:nv(ispin, ikpt), ikpt, ispin)), GbasVec(3, ilst(:nv(ispin, ikpt), ikpt, ispin)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi_type, ods_type, a, b, bascof_lo )
      matchingCoeff(:,1, :, :) = a
      matchingCoeff(:,2, :, :) = b
      iatom = 0
      ! Sort the matching coefficients, todo bascof_lo should also come into these matching coefficients. But they are used later so LOs
      ! are built in
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + mqn_m
              matchingCoeff(:nv(ispin, ikpt), 1, lm, iatom) = a(:nv(ispin, ikpt), lm, iatom)
              matchingCoeff(:nv(ispin, ikpt), 2, lm, iatom) = b(:nv(ispin, ikpt), lm, iatom)
            end do
          end do
        end do
      end do

    ! overlapMat muss reel sein fuer inversions-symmetrische Systeme

    ! IR part
    !
    overlapIndex = 0; overlapMat = 0
    do iG = 1, nv(ispin, ikpt)
      do iGprime = 1, iG
        overlapIndex = overlapIndex + 1
        if ( iG /= iGprime ) then
          Gdiff_int  = GbasVec(:, ilst(iG, ikpt, ispin)) - GbasVec(:, ilst(iGprime, ikpt, ispin))
          call cotra3( real(Gdiff_int), Gdiff_ext,cell%bmat )
          Gdiff_norm = norm2(Gdiff_ext)
          contrI = 0d0; iatom = 0
          do itype = 1, atoms%ntype
            argument_temp = Gdiff_norm * atoms%rmt(itype)
            rdum          = (sin(argument_temp) - argument_temp * cos(argument_temp)) / Gdiff_norm**3
            rdum          = fpi * rdum / cell%vol
            do ieqat = 1, atoms%neq(itype)
              iatom  = iatom + 1
              contrI = contrI - rdum *  exp( +iu * tpi * dot_product( atoms%taual(:,iatom), real(Gdiff_int) ) )
            end do
          end do
        else
          contrI = cell%volint / cell%vol
        endif
        overlapMat(overlapIndex) = contrI
      end do
    end do


    !MT part
    iatom0 = 0
    do itype = 1, atoms%ntype
      lm = -1
      ng = atoms%jri(itype)
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = - oqn_l, oqn_l
          lm = lm + 1

          do p1 = 1,2
            do p2 = 1,2
              ! dies ist unabhaengig von m --> daher waere es sogar sinnvoll den m-Loop nach den p1,p2 Loop durchzufuehren
              ! man koennte sogar die Normiertheit von (1/1) und die Orthogonalitaet von (1/2) ausnutzen
              f(:ng) = rbas1(:ng, p1, oqn_l, itype, ispin) * rbas1(:ng, p2, oqn_l, itype, ispin) + rbas2(:ng, p1, oqn_l, itype, ispin) &
                & * rbas2(:ng, p2, oqn_l, itype, ispin)
              rdum = Intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
             ! if( mqn_m .eq. oqn_l ) write(*,*) oqn_l,p1,p2,rdum !todo this can be uncommented
!               if( p1.eq.p2 ) then
!                 if( p1.eq. 1 ) then
!                   rdum = 1d0
!                 else
!                   f(:ng) = rbas1(:ng, p1, oqn_l, itype, ispin) * rbas1(:ng, p2, oqn_l, itype, ispin) + &
!                     & rbas2(:ng, p1, oqn_l, itype, ispin) * rbas2(:ng, p2, oqn_l, itype, ispin)
!                   rdum   = Intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)
!                 endif
!               else
!                 rdum = 0.
!               endif

              iatom = iatom0
              do ieqat = 1,atoms%neq(itype)
                iatom = iatom + 1
                ! speicher wird nicht linear durchlaufen
                overlapIndex = 0
                do iG = 1,nv(ispin, ikpt)
                  cdum = rdum * matchingCoeff(iG, p2, lm, iatom)
                  do iGprime = 1, iG
                    overlapIndex = overlapIndex + 1
                    overlapMat(overlapIndex) = overlapMat(overlapIndex) + conjg(matchingCoeff(iGprime, p1, lm, iatom)) * cdum
                  end do  !iGprime
                end do  !iG

              end do  !ieqat

            end do  !p2
          end do  !p1

         end do  !mqn_m
       end do  !oqn_l

       iatom0 = iatom0 + atoms%neq(itype)
     end do  !itype

    ! calculate overlaps of radial basis functions
    !
    ! den radialen Überlapp wird du vermutlich an vielen Stellen benötigen, daher würde ich
    ! ihn einmal in init oder so berechnen und im Speicher halten

    f = 0
    !radbasFun are real, so overlap is invariant with respect to complex conjugation
    allocate( overlapRadFun( 2 + maxval(nRadFun), 2 + maxval(nRadFun), 0: atoms%lmaxd, atoms%nat, input%jspins) )
    overlapRadFun = 0
    do jspin = 1, input%jspins
      do itype = 1, atoms%ntype
        do ioqn_l = 0, atoms%lmax(atoms%ntype)
          overlapRadFun(1, 1, ioqn_l, itype, jspin) = 1

          overlapRadFun(1, 2, ioqn_l, itype, jspin) = 0
          overlapRadFun(2, 1, ioqn_l, itype, jspin) = 0

          f(:) = rbas1(:, 2, ioqn_l, itype, jspin) * rbas1(:, 2, ioqn_l, itype, jspin) + rbas2(:, 2, ioqn_l, itype, jspin) &
            & * rbas2(:, 2, ioqn_l, itype, jspin) !cmplx conj! are rbas1 and rbas2 complex when no inversion symmetry
          overlapRadFun(2, 2, ioqn_l, itype, jspin) = intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

          do p = 3, nRadFun(ioqn_l, itype)
            iPtable(iloTable(p, ioqn_l, itype), itype) = p
            overlapRadFun(p, p, ioqn_l, iatom, jspin) = 1

            f(:) = rbas1(:, 1, ioqn_l, itype, jspin) * rbas1(:, p, ioqn_l, itype, jspin) + rbas2(:, 1, ioqn_l, itype, jspin) &
              &* rbas2(:, p, ioqn_l, itype, jspin)
            overlapRadFun(1, p, ioqn_l, itype, jspin) = intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

            f(:) = rbas1(:, 2, ioqn_l, itype, jspin) * rbas1(:, p, ioqn_l, itype, jspin) + rbas2(:, 2, ioqn_l, itype, jspin) * &
              & rbas2(:, p, ioqn_l, itype, jspin)
            overlapRadFun(2, p, ioqn_l, itype, jspin) = intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, gridf)

            do pPrime = 3, nRadFun(ioqn_l, itype)
              f(:) = rbas1(:, p, ioqn_l, itype, jspin) * rbas1(:, pPrime, ioqn_l, itype, jspin) + rbas2(:, p, ioqn_l, itype, jspin) * &
                & rbas2(:, pPrime, ioqn_l, itype, jspin) !cmplx conj! are rbas1 and rbas2 complex when no inversion symmetry
              overlapRadFun(pPrime, p, ioqn_l, itype, jspin) = intgrf(f, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, &
                & itype, gridf)
            end do
          end do
        end do
      end do
    end do

    iatom = 0
    allocate( overlapGGLo(atoms%nlotot, nv(ispin, ikpt) ) )
    allocate( overlapGLOGLO(atoms%nlotot *  ( atoms%nlotot + 1 ) / 2) )
    overlapGGLo = 0
    overlapGLOGLO = 0
    overlapIndex = 0 ! maybe this has to be initialized for every atom
    loGshift = 0
    jdum = 0
    idum = 0
    loShift = 0
    atomShift = 0

    ! calculate the analytical LAPW-LO overlap matrix in basis of LAPW basis functions.
    do itype = 1, atoms%ntype
!      idum = 0
      do ieq = 1, atoms%neq(itype)
        iatom = iatom + 1
        atomShift = 0
        !nkvecat = 0
        !if ( atoms%invsat(iatom) == 2 ) then ! if atom is reached by inversion symmetry then skip because it is already considered !todo si with inversion symmetry works not but gaas does
        !  cycle
        !endif
        do ilo = 1, atoms%nlo(itype) ! run over all local orbitals of this atom type
          l = atoms%llo(ilo, itype)
          nGLOs = 2 * l + 1 ! we have 2 l + 1 functions per l
          !if ( atoms%invsat(iatom) == 1 ) then ! here we consider atoms which can reached by inversion symmetry and were skipped before, we have one counteratom which can be reached by inversion symmetry by another parent atom
          !  nGLOs = 2 * nGLOs
          !endif
          ! submatrix GGLO
          do iGloc = 1, nGLOs
            do iG = 1, nv(ispin, ikpt) ! es würde sich anbieten die Loop Reihenfolge von iG und mqn_m zu tauschen;
                                       ! dann könnte man   bascof_lo(1, mqn_m, iGloc, ilo, iatom) &
                                       !                 + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * overlapRadFun(1, 2 + iloTable(l, itype), l, itype, ispin))
                                       ! etc. bereits ausführen
                                       ! loop over G vectors for a certain kpoint and spin
              lm = l * l - 1 ! this formula reaches the adequate lm for a given l
              do mqn_m = -l, l
                lm = lm + 1
                overlapGGLO(iGloc + loShift, iG) = overlapGGLO(iGloc + loShift, iG) & ! contribution from parent atom
                  &+ conjg(matchingCoeff(iG, 1, lm, iatom)) * (bascof_lo(1, mqn_m, iGloc, ilo, iatom) &
                  &+ bascof_lo(3, mqn_m, iGloc, ilo, iatom) * overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) &
                  &+ conjg(matchingCoeff(iG, 2, lm, iatom)) * (bascof_lo(2, mqn_m, iGloc, ilo, iatom) * &
                  &overlapRadFun(2, 2, l, itype, ispin) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
                  &overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin))
                if (atoms%invsat(iatom) == 1) then ! contribution from child atom mapped to by inversion symmetry
                  jatom = sym%invsatnr(iatom)
                  overlapGGLO(iGloc + loShift, iG) = overlapGGLO(iGloc + loShift, iG) &
                    &+ conjg(matchingCoeff(iG, 1, lm, jatom)) * (bascof_lo(1, mqn_m, iGloc, ilo, jatom) &
                    &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) * overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) &
                    &+ conjg(matchingCoeff(iG, 2, lm, jatom)) * (bascof_lo(2, mqn_m, iGloc, ilo, jatom) &
                    &* overlapRadFun(2, 2, l, itype, ispin) + bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
                    &* overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin))
                endif
              end do ! mqn_m
            end do ! iG
          end do ! iGloc
          loShift = loShift + nGLOs
        end do ! ilo

          ! This is the LO LO submatrix in the basis of LAPW basis functions.
        do ilo = 1, atoms%nlo(itype)
          l = atoms%llo(ilo, itype)
          nGLOs = 2 * l + 1 ! we have 2 l + 1 functions per l
          do iGloc = 1, nGLOs
           !go to the next atom in this linear index, LOs only have an overlap within one MT!
            overlapIndex = overlapIndex + atomShiftprev
            !loGprimeShift = loGShift ! understand this, loGprimeshift is reset in the first run
            !overlapIndex = overlapIndex + ( itype - 1 ) * loGshift ! is zero in the first run
            do jlo = 1, (ilo - 1)!this is for submatrix LO LO only triangular matrix is considered ! this order of loops only works for one LO we have to fill up G and Gprime for 1 LO first then for the next one
              l_jlo = atoms%llo(jlo, itype)
              nGLO_jlo = 2 * l_jlo + 1
              !if ( atoms%invsat(iatom) == 1 ) then
              !  nGLO_jlo = 2 * nGLO_jlo
              !endif

              if ( l_jlo /= l ) then !overlap is only for same l
                overlapIndex = overlapIndex + nGLO_jlo
                cycle
              endif

              do iGlocprime = 1, nGLO_jlo
              !  if( iGlocprime + loGprimeShift >  iGloc + idum ) then ! understand this only fill up lower triangular matrix? this induces a lot of unneeded loops
              !    cycle
              !  endif

                overlapIndex = overlapIndex +  1

                do mqn_m = - l, l
                  overlapGLOGLO(overlapIndex) = overlapGLOGLO(overlapIndex) + conjg(bascof_lo(1, mqn_m, iGlocprime, jlo, iatom)) * &
                    &(bascof_lo(1, mqn_m, iGloc, ilo, iatom) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
                    & overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) + conjg(bascof_lo(2, mqn_m, iGlocprime, jlo, iatom)) * &
                    & (bascof_lo(2, mqn_m, iGloc, ilo, iatom) * overlapRadFun(2, 2, l, itype, ispin) + &
                    & bascof_lo(3, mqn_m, iGloc, ilo, iatom) * overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin)) &
                    & + conjg( bascof_lo(3, mqn_m, iGlocprime, jlo, iatom) ) * ( bascof_lo(1, mqn_m, iGloc, ilo, iatom) * &
                    & overlapRadFun(1, iPtable(jlo, itype), l, itype, ispin) + bascof_lo(2, mqn_m, iGloc, ilo, iatom) * &
                    & overlapRadFun(2, iPtable(jlo, itype), l, itype, ispin) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
                    & overlapRadFun(iPtable(ilo, itype), iPtable(jlo, itype),  l, itype, ispin) )
                  if (atoms%invsat(iatom)==1) then
                    jatom = sym%invsatnr(iatom)
                    overlapGLOGLO(overlapIndex) = overlapGLOGLO(overlapIndex) + conjg(bascof_lo(1, mqn_m, iGlocprime, jlo, jatom)) * &
                      &(bascof_lo(1, mqn_m, iGloc, ilo, jatom) + bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
                      &* overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) + conjg(bascof_lo(2, mqn_m, iGlocprime, jlo, jatom)) &
                      &* (bascof_lo(2, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, 2, l, itype, ispin) &
                      &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin)) &
                      &+ conjg(bascof_lo(3, mqn_m, iGlocprime, jlo, jatom)) * (bascof_lo(1, mqn_m, iGloc, ilo, jatom) &
                      &* overlapRadFun(1, iPtable(jlo, itype), l, itype, ispin) &
                      &+ bascof_lo(2, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, iPtable(jlo, itype), l, itype, ispin) &
                      &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
                      &* overlapRadFun(iPtable(ilo, itype), iPtable(jlo, itype), l, itype, ispin))
                  endif
                enddo !mqn
              enddo ! iGlocprime
            enddo !jlo
            !LO with itself
            do iGlocprime = 1, iGloc
            !  if( iGlocprime + loGprimeShift >  iGloc + idum ) then ! understand this only fill up lower triangular matrix? this induces a lot of unneeded loops
            !    cycle
            !  endif

              overlapIndex = overlapIndex +  1

              do mqn_m = - l, l
                overlapGLOGLO(overlapIndex) = overlapGLOGLO(overlapIndex) + conjg(bascof_lo(1, mqn_m, iGlocprime, jlo, iatom)) * &
                  &(bascof_lo(1, mqn_m, iGloc, ilo, iatom) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
                  & overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) + conjg(bascof_lo(2, mqn_m, iGlocprime, jlo, iatom)) * &
                  & (bascof_lo(2, mqn_m, iGloc, ilo, iatom) * overlapRadFun(2, 2, l, itype, ispin) + &
                  & bascof_lo(3, mqn_m, iGloc, ilo, iatom) * overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin)) &
                  & + conjg( bascof_lo(3, mqn_m, iGlocprime, jlo, iatom) ) * ( bascof_lo(1, mqn_m, iGloc, ilo, iatom) * &
                  & overlapRadFun(1, iPtable(jlo, itype), l, itype, ispin) + bascof_lo(2, mqn_m, iGloc, ilo, iatom) * &
                  & overlapRadFun(2, iPtable(jlo, itype), l, itype, ispin) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
                  & overlapRadFun(iPtable(ilo, itype), iPtable(jlo, itype),  l, itype, ispin) )
                if (atoms%invsat(iatom)==1) then !todo for systems with inversion symmetry this has to be tested
                  jatom = sym%invsatnr(iatom)
                  overlapGLOGLO(overlapIndex) = overlapGLOGLO(overlapIndex) + conjg(bascof_lo(1, mqn_m, iGlocprime, jlo, jatom)) * &
                    &(bascof_lo(1, mqn_m, iGloc, ilo, jatom) + bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
                    &* overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) + conjg(bascof_lo(2, mqn_m, iGlocprime, jlo, jatom)) &
                    &* (bascof_lo(2, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, 2, l, itype, ispin) &
                    &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin)) &
                    &+ conjg(bascof_lo(3, mqn_m, iGlocprime, jlo, jatom)) * (bascof_lo(1, mqn_m, iGloc, ilo, jatom) &
                    &* overlapRadFun(1, iPtable(jlo, itype), l, itype, ispin) &
                    &+ bascof_lo(2, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, iPtable(jlo, itype), l, itype, ispin) &
                    &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
                    &* overlapRadFun(iPtable(ilo, itype), iPtable(jlo, itype), l, itype, ispin))
                endif
              enddo !mqn
            enddo ! iGlocprime
          end do !iGloc
          ! Determine how much steps we have to jump for the current atom to get to the next atom in the linear index
          atomShift = atomShift + nGLOs
        end do !ilo
        ! Determine how much steps we have to jump to get to the current atom skipping all the atoms before in the linear index
        atomShiftprev = atomShiftprev + atomShift


            !do iGlocprime = 1, nGLO_jlo
            !!  if( iGlocprime + loGprimeShift >  iGloc + idum ) then ! understand this only fill up lower triangular matrix? this induces a lot of unneeded loops
            !!    cycle
            !!  endif

            !  overlapIndex = overlapIndex +  1
            !  do mqn_m = - l, l
            !    overlapGLOGLO(overlapIndex) = overlapGLOGLO(overlapIndex) + conjg(bascof_lo(1, mqn_m, iGlocprime, jlo, iatom)) * &
            !      &(bascof_lo(1, mqn_m, iGloc, ilo, iatom) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
            !      & overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) + conjg(bascof_lo(2, mqn_m, iGlocprime, jlo, iatom)) * &
            !      & (bascof_lo(2, mqn_m, iGloc, ilo, iatom) * overlapRadFun(2, 2, l, itype, ispin) + &
            !      & bascof_lo(3, mqn_m, iGloc, ilo, iatom) * overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin)) &
            !      & + conjg( bascof_lo(3, mqn_m, iGlocprime, jlo, iatom) ) * ( bascof_lo(1, mqn_m, iGloc, ilo, iatom) * &
            !      & overlapRadFun(1, iPtable(jlo, itype), l, itype, ispin) + bascof_lo(2, mqn_m, iGloc, ilo, iatom) * &
            !      & overlapRadFun(2, iPtable(jlo, itype), l, itype, ispin) + bascof_lo(3, mqn_m, iGloc, ilo, iatom) * &
            !      & overlapRadFun(iPtable(ilo, itype), iPtable(jlo, itype),  l, itype, ispin) )
            !    if (atoms%invsat(iatom)==1) then
            !      jatom = sym%invsatnr(iatom)
            !      overlapGLOGLO(overlapIndex) = overlapGLOGLO(overlapIndex) + conjg(bascof_lo(1, mqn_m, iGlocprime, jlo, jatom)) * &
            !        &(bascof_lo(1, mqn_m, iGloc, ilo, jatom) + bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
            !        &* overlapRadFun(1, iPtable(ilo, itype), l, itype, ispin)) + conjg(bascof_lo(2, mqn_m, iGlocprime, jlo, jatom)) &
            !        &* (bascof_lo(2, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, 2, l, itype, ispin) &
            !        &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, iPtable(ilo, itype), l, itype, ispin)) &
            !        &+ conjg(bascof_lo(3, mqn_m, iGlocprime, jlo, jatom)) * (bascof_lo(1, mqn_m, iGloc, ilo, jatom) &
            !        &* overlapRadFun(1, iPtable(jlo, itype), l, itype, ispin) &
            !        &+ bascof_lo(2, mqn_m, iGloc, ilo, jatom) * overlapRadFun(2, iPtable(jlo, itype), l, itype, ispin) &
            !        &+ bascof_lo(3, mqn_m, iGloc, ilo, jatom) &
            !        &* overlapRadFun(iPtable(ilo, itype), iPtable(jlo, itype), l, itype, ispin))
            !    endif
         !     enddo !mqn
         !   enddo ! iGlocprime
       !   enddo !iGloc
          !nkvecat = nkvecat + nGLOs

         ! idum = idum + nGLOs !this variable is probably to less used
         ! jdum = jdum + nGLOs
        !nkvecprevat = nkvecprevat + nkvecat
      end do !ieq
      !loGshift = loGshift + idum ! this should be within ieq not within itype or not?
      !idum = loGshift
    end do !itype
    !
    ! test: wave function orthogonality
    !

    ! unpack overlapMat
    overlapIndex = 0
    smat = 0
    do iG = 1, nv(ispin, ikpt)
      locGindex = 0
      do iGprime = 1, iG
          overlapIndex = overlapIndex + 1
          smat(iGprime, iG) = overlapMat(overlapIndex)
          smat(iG, iGprime) = conjg( overlapMat(overlapIndex) )
      end do
    end do

    overlapIndex = 0
    !nGLOs = 2 * (2 * atoms%llod + 1)
      do iGprime = nv(ispin, ikpt) + 1, nv(ispin, ikpt) + atoms%nlotot
        do iG = nv(ispin, ikpt) + 1, iGprime
        overlapIndex =  overlapIndex + 1
        smat(iGprime, iG) = conjg(overlapGLOGLO(overlapIndex))
        smat(iG, iGprime) =  overlapGLOGLO(overlapIndex)
      end do
    end do


    nGLOs = 2 * (2 * atoms%llod + 1)
    do iG = 1, nv(ispin, ikpt)
      do iGprime = nv(ispin, ikpt) + 1, nv(ispin, ikpt) + atoms%nlotot
         smat(iGprime, iG) = conjg(overlapGGLO(iGprime - nv(ispin, ikpt), iG))
         smat(iG, iGprime) =  overlapGGLO(iGprime - nv(ispin, ikpt), iG)
      end do
    end do


    !debugging code
!    if (ikpt == 28) then
!      do idump = 1, nv(ispin, ikpt) + atoms%nlotot
!        do jdump = 1, idump
!          write(12346,"((2(i8, 2x), 2(f18.10, 2x)))") idump, jdump, smat(idump, jdump)
!        end do
!      end do
!    endif
    !smat = real(smat) !todo ??????


    !Multiply the zs and test if the result is the unit matrix
    passed = .true.
    allocate( carr(nv(ispin, ikpt) + atoms%nlotot) )
!    write(*, *) 'nv + glos', nv(ispin, ikpt) + atoms%nlotot
    !do itype = 1, atoms%ntype
      do iband1 = 1,ne(ikpt) !to 4 when comparing with Sternheimer matrix element
        carr = matmul( smat, z(:nv(ispin, ikpt) + atoms%nlotot, iband1) )
        do iband2 = 1,ne(ikpt)
         overlap = dot_product( z(:nv(ispin, ikpt) + atoms%nlotot,iband2), carr )
          !todo this has to be conjugated right?
!         write(*,*) iband2, iband1, overlap
          if ( iband1 == iband2 ) then
            if ( ( real(overlap) - 1 > 1e-9 ) .or. ( aimag(overlap) > 1e-9) ) then
              write (*, '(a,x,i2,x,a,x,i3,x,a,x,i3,x,a,x,i1)') 'not orthogonal at k-point =', ikpt, 'with bands', iband1, 'and', &
                & iband2, 'with spin', ispin
              passed = .false.
            end if
          else
            if ( (abs(overlap) > 1e-9) ) then
              write (*, '(a,x,i2,x,a,x,i3,x,a,x,i3,x,a,x,i1)') 'not orthogonal at k-point =', ikpt, 'with bands', iband1, 'and', &
                & iband2, 'with spin', ispin
              passed = .false.
            end if
          end if
          !write (446, '(i4,i4,2(f24.12))') iband2, iband1, overlap
        end do
      end do
    !end do
    if ( .not. passed ) then
      write (logUnit, '(a,i3,a)') 'k-point', ikpt, 'failed'
      call JuDFT_warn( 'OverlapMatrix Input Test failed', calledby='OverlapMatrix', &
        &hint='Debug! Check input data!', file='jpTestInput_mod.F90' )
    end if

  end subroutine OverlapMatrix

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Test of m_jppotdens::genrotllh and m_jppotdens::rotate_clnu.
  !>
  !> @details
  !> Test of m_jppotdens::genrotllh and m_jppotdens::rotate_clnu by printing out the radial spherical harmonic expansion coefficients of
  !> the muffin-tin electronic density. Although convertLH2SphHarm is never called within the code, it serves as a test of the altered
  !> lattice harmonic quantities.
  !>
  !> @note This test requires additional output of FLEUR namely the same system but with an artificial reduction of symmetry. There
  !> shoud be only one atom per atom type.
  !>
  !> @param[out] rho0MT     : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] input      : Input type, see types.f90.
  !> @param[out] lathar     : Lattice harmonics type, see types.f90.
  !> @param[out] mlh_atom   : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom  : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom  : Phase mediating between stars and plane waves.
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CheckBackRotLHLCS ( rho0MT, atoms, input, lathar, clnu_atom, nmem_atom, mlh_atom )

    use m_types
    use m_jpPotDensHelper, only : ConvertLH2SphHarm

    implicit none

    ! Scalar Arguments
    type(t_atoms),                 intent(in) :: atoms
    type(t_input),                 intent(in) :: input
    type(t_sphhar),                intent(in) :: lathar

    ! Array Arguments
    complex,                       intent(in) :: clnu_atom(:, 0:, :)
    integer,                       intent(in) :: nmem_atom(0:, :)
    integer,                       intent(in) :: mlh_atom(:, 0:, :)
    real,                          intent(in) :: rho0MT(:,0:,:,:)

    ! Local Scalar Variables !todo documentation
    integer                                   :: iatom
    integer                                   :: jspin
    integer                                   :: itype
    integer                                   :: ieq
    integer                                   :: oqn_l
    integer                                   :: mqn_m
    integer                                   :: lm
    integer                                   :: igrid

    ! Local Array Variables
    character(len=16)                         :: filenameR
    character(len=16)                         :: filenameI
    complex,          allocatable             :: rho0MT_atom(:, :, :, :)

    allocate( rho0MT_atom( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
    ! Call of convertLH2SphHarm as a test of the changed lattice harmonic quantities and output of the real and imaginary part.
    do jspin = 1, input%jspins
      iatom = 0
      do itype = 1, atoms%ntype
        do ieq = 1, atoms%neq(itype)
          iatom =  iatom + 1
          !todo changed in a recent commit
          call ConvertLH2SphHarm( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, rho0MT_atom )
          write (filenameR, "(A10,I1,A5)") 'radDenslmA', iatom, 'R.ssv'
          write (filenameI, "(A10,I1,A5)") 'radDenslmA', iatom, 'I.ssv'
          call Fopen(1000, name=filenameR, status='replace', action='write', form='formatted')
          call Fopen(1001, name=filenameI, status='replace', action='write', form='formatted')
          write(1000, '(A56,I1,A25,I1,A11,I3)') '# Radial density coefficients for all lms for atom type ', itype, &
            &                                                ' and its equivalent atom ', ieq, '. Columns: ', (atoms%lmax(itype) + 1)**2
          write(1001, '(A56,I1,A25,I1,A11,I3)') '# Radial density coefficients for all lms for atom type ', itype, &
            &                                                ' and its equivalent atom ', ieq, '. Columns: ', (atoms%lmax(itype) + 1)**2
          do igrid = 1, atoms%jri(itype)
            write(1000,'(2x,*(es15.7))') atoms%rmsh(igrid, itype), &
              & ( real( rho0MT_atom(lm, igrid, iatom, 1) ), lm = 1, ( atoms%lmax(itype) + 1 )**2 )
            write(1001,'(2x,*(es15.7))') atoms%rmsh(igrid, itype), &
              & ( aimag( rho0MT_atom(lm, igrid, iatom, 1) ), lm = 1, ( atoms%lmax(itype) + 1 )**2 )
          end do ! igrid
          call Fclose(1000)
          call Fclose(1001)
        end do ! ieq
      end do ! itype
    end do ! jspin
  end subroutine CheckBackRotLHLCS

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Calcualtes integrals from xc-quantities to compare with out file from FLEUR.
  !>
  !> @details
  !> Calculates integrals of density and Exc or density and Vxc, respectively. The output is also printed to the out file of a FLEUR
  !> calculation. Tests also the routines to calculate 2-argument IR and MT integrals.
  !>
  !> @param[out] atoms      : Atoms type, see types.f90.
  !> @param[out] stars      : Stars type, see types.f90.
  !> @param[out] cell       : Unit cell type, see types.f90.
  !> @param[out] lathar     : Lattice harmonics type, see types.f90.
  !> @param[in]  logUnit    : Unit number for juPhon.log.
  !> @param[out] rho0IR     : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur.
  !> @param[out] rho0MT     : Radial coefficients of the unperturbed and converged muffin-tin densities parsed from Fleur.
  !> @param[out] mlh_atom   : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[out] nmem_atom  : Number of lattice harmonic members for every atom.
  !> @param[out] clnu_atom  : Phase mediating between stars and plane waves.
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcXCintegrals( atoms, cell, stars, lathar, ngdp, rho0IRst, vXC0IRst, vXC0MTlh, eXCIRst, eXCMTlh, rho0MTlh, gdp, mlh_atom, nmem_atom, clnu_atom, logUnit )

    use m_types
    use m_jpPotDensHelper, only : WarpIRPot, ConvertStar2G
    use mod_juPhonUtils, only : Calc2ArgIntIR, Calc2ArgCmplxIntMT

    implicit none

    ! Type parameter
    type(t_atoms),               intent(in) :: atoms
    type(t_cell),                intent(in) :: cell
    type(t_stars),               intent(in) :: stars
    type(t_sphhar),              intent(in) :: lathar

    ! Scalar parameter
    integer,                     intent(in) :: logUnit
    integer,                     intent(in) :: ngdp

    ! Array parameter

    complex,                     intent(in) :: vXC0IRst(:, :)
    complex,                     intent(in) :: eXCIRst(:)
    complex,                     intent(in) :: rho0IRst(:, :)
    real,                        intent(in) :: vXC0MTlh(:, 0:, :, :)
    real,                        intent(in) :: eXCMTlh(:, 0:, :)
    real,                        intent(in) :: rho0MTlh(:, 0:, :, :)
    integer,                    intent(in)  :: mlh_atom(:, 0:, :)
    integer,                    intent(in)  :: nmem_atom(0:, :)
    complex,                    intent(in)  :: clnu_atom(:, 0:, :)
    integer,                     intent(in) :: gdp(:, :)

    ! Scalar variables
    integer                                 :: idirR
    complex                                 :: integralIRV
    complex                                 :: integralIRE
    integer                                 :: ptsym
    integer                                 :: ilh
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: lm_pre
    integer                                 :: imesh
    complex                                 :: integralMTv
    complex                                 :: integralMTe
    complex                                 :: integralMTvSum
    complex                                 :: integralMTeSum
    integer                                 :: iatom
    integer                                 :: itype
    integer                                 :: ieqat
    integer                                 :: imem
    integer                                 :: lm

    integer                                 :: iG

    ! Array variables
    complex,        allocatable             :: vXC0IRpwContainer(:, :)
    complex,        allocatable             :: eXCIRpwContainer(:, :)
    complex,        allocatable             :: w_vXC0IRpwContainer(:)
    complex,        allocatable             :: w_eXCIRpwContainer(:)
    complex,        allocatable             :: vXC0MTSpH(:, :, :)
    complex,        allocatable             :: eXCMTSpH(:, :, :)
    complex,        allocatable             :: rho0MTSpH(:, :, :)
    complex,        allocatable             :: rho0IRpw(:)

    !todo this routine can be further optimized
    write(logUnit,'(a)')
    write(logUnit,'(a)') 'Calculate IR and MT integrals with integrand rho V_xc (rhoVxc) and integrand rho E_xc (rhoExc):'
    write(logUnit,'(a)') '-----------------------------------------------------------------------------------------------'

    idirR = 1
    allocate( vXC0IRpwContainer(ngdp, 3), eXCIRpwContainer(ngdp, 3), w_vXC0IRpwContainer(ngdp), w_eXCIRpwContainer(ngdp) )
    allocate( vXC0MTSpH( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat ), eXCMTSpH( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat ), &
            & rho0MTSpH( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat ) )
    allocate( rho0IRpw(ngdp) )

    vXC0IRpwContainer(:, :) = cmplx(0., 0.)
    eXCIRpwContainer(:, :) = cmplx(0., 0.)
    w_vXC0IRpwContainer(:) = cmplx(0., 0.)
    w_eXCIRpwContainer(:) = cmplx(0., 0.)
    ! Convert star expansion coefficients to plane-wave expansion coefficients
    call convertStar2G(vXC0IRst(:, 1), vXC0IRpwContainer(:, idirR), stars, ngdp, gdp)
    call convertStar2G(rho0IRst(:, 1), rho0IRpw, stars, ngdp, gdp)
    call convertStar2G(eXCIRst(:), eXCIRpwContainer(:, idirR), stars, ngdp, gdp)

    ! Warp xc potential and energy density
    call warpIRPot( stars, ngdp, idirR, gdp, vXC0IRpwContainer, w_vXC0IRpwContainer(:) )
    call warpIRPot( stars, ngdp, idirR, gdp, eXCIRpwContainer, w_eXCIRpwContainer(:) )

    ! Calculate integrals
    call Calc2ArgIntIR( cell, ngdp, rho0IRpw, w_vXC0IRpwContainer, integralIRV)
    call Calc2ArgIntIR( cell, ngdp, rho0IRpw, w_eXCIRpwContainer, integralIRE)

    if (.false.) then
      write(logUnit,'(a,1x,2f15.8)') 'rhoVxcIR:', integralIRV
      write(logUnit,'(a,1x,2f15.8)') 'rhoExcIR:', integralIRE
      write(logUnit,'(a)')
    end if

    ! Convert lattice harmonic expansion coefficients to spherical harmonic expansion coefficients
    vXC0MtSpH(:, :, :) = cmplx(0.0, 0.0)
    eXCMTSpH(:, :, :) = cmplx(0.0, 0.0)
    rho0MTSpH(:, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m + 1
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(itype)
              vXC0MtSpH(imesh, lm, iatom) = vXC0MtSpH(imesh, lm, iatom) + vXC0MTlh(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              rho0MTSpH(imesh, lm, iatom) = rho0MTSpH(imesh, lm, iatom) + rho0MTlh(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              eXCMtSpH(imesh, lm, iatom) = eXCMtSpH(imesh, lm, iatom) + eXCMTlh(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

    ! Calculate integrals
    integralMTvSum = cmplx(0., 0.)
    integralMTeSum = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            call Calc2ArgCmplxIntMT( atoms, itype, rho0MTSpH(:, lm, iatom), vXC0MTSpH(:, lm, iatom), integralMTv )
            call Calc2ArgCmplxIntMT( atoms, itype, rho0MTSpH(:, lm, iatom), eXCMTSpH(:, lm, iatom), integralMTe )
            integralMTvSum = integralMTvSum + integralMTv
            integralMTeSum = integralMTeSum + integralMTe
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

    if (.false.) then
      write(logUnit,'(a,1x,2f15.8)') 'rho0VxcMT', integralMTvSum
      write(logUnit,'(a,1x,2f15.8)') 'rhoExcMT', integralMTeSum
    end if

    write(logUnit,'(a,1x,2f15.8)') 'rho0Vxc', integralIRV + integralMTvSum
    write(logUnit,'(a,1x,2f15.8)') 'rhoExc', integralIRE + integralMTeSum
    write(logUnit,*)
    write(logUnit,'(a)') 'Please compare to the output of FLEUR, where the same integrals are calculated!'
    write(logUnit,*)

  end subroutine CalcXCintegrals

end module m_jpTestInput
