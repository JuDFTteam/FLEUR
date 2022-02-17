!----------------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!----------------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Contains IO routines for inputfile and logfile.
!
!> @author
!> Christian-Roman Gerhorst and ChriNOstopNOh Friedrich
!
!> @brief
!> Contains IO routines for inputfile and logfile.
!>
!> @note
!> Will be substituted by the xml input developed for FLEUR in the future. Therefore only minimal effort achieving modern code,
!> readability and documentation will be put into here.
!----------------------------------------------------------------------------------------------------------------------------------------
module m_jpLog

  use mod_juPhonUtils, only : fopen, fclose

  implicit none

  contains

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Generates initial output of juPhon.log.
  !>
  !> @details
  !> This routine writes information about starting-date, starting-time and architecture on which juPhon is run to the beginning of the
  !> log-file juphon.log.
  !>
  !> @param[in]  logFUnit  : Global unit number of juPhon.log
  !> @param[out] startTime : Starting time of juPhon
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine InitLogFile( logFUnit, startTime )

    implicit none

    ! Scalar parameter
    integer, intent (in)  :: logFUnit

    ! Array parameter
    integer, intent (in)  :: startTime(8)

    ! Local array
    character(len=30)     :: hostName     ! architecture on which juPhon is running


    ! Acquire starting-date, starting-time and hostname
    call hostnm( hostName )

    ! Init logfile
    write( logFUnit, '(a)' ) '********************'
    write( logFUnit, '(a)' ) '*  juPhon logfile  *'
    write( logFUnit, '(a)' ) '********************'
    write( logFUnit, * )
    write( logFUnit, '(a)' ) 'juPhon Calculation started!'
    write( logFUnit, '(a)' ) '***************************'
    write( logFUnit, * )
    write( logFUnit, '(a16,a2,a1,a2,a1,i4,a4,a2,a1,a2,a1,a2,a1,a3)' ) 'Date and Time:  ', zero2l(startTime(3)), '.', &
      & zero2l(startTime(2)), '.', startTime(1), ' at ', zero2l(startTime(5)), ':', zero2l(startTime(6)), ':', zero2l(startTime(7)), &
      & ':', zero3l(startTime(8))
    write ( logFUnit, '(a16,a)' ) 'Hostname:       ', trim(hostName)
    write ( logFUnit, * )

#ifdef DEBUG_MODE
    write (logFUnit, '(a)') 'Debug mode activated!'
#endif DEBUG_MODE
#ifdef CPP_INVERSION
    write (logFUnit, '(a)') 'juPhon compiled for systems with inversion symmetry!'
#else
    write (logFUnit, '(a)') 'juPhon compiled without simplifications for inversion symmetry!'
#endif

  end subroutine InitLogFile

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Adds leading zero to get a 2-digit integer
  !>
  !> @param[in] intNr: Global unit number of juPhon.log
  !> @return  zero2l: Formatted string
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  function Zero2l(intNr)

    implicit none

    ! Scalar parameter
    integer, intent(in) :: intNr

    ! Scalar variable
    character(len=2)    :: zero2l

    if ( intNr < 10 ) then
      write(zero2l, '(a1,i1)') '0', intNr

    else
      write(zero2l, '(i2)') intNr

    end if

  end function Zero2l

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Adds leading zero to get a 2-digit integer
  !>
  !> @param[in] intNr: Global unit number of juPhon.log
  !> @returns  zero3l: Formatted string
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
  function Zero3l(intNr)

    implicit none

    ! Scalar parameter
    integer, intent(in) :: intNr

    ! Scalar variable
    character(len=3)    :: zero3l

    if ( intNr < 10 ) then
      write(zero3l, '(a2,i1)') '00', intNr

    else if (intNr < 100) then
      write(zero3l, '(a1,i2)') '0', intNr

    else
      write(zero3l, '(i3)') intNr

    end if

  end function Zero3l

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Writes out read-in data from juPhon.inp to log file juPhon.log
  !>
  !> @param[in] kModeSw:   k-point generation mode switch
  !> @param[in] kSetDim:   dimension of k-point set
  !> @param[in] qSetDim:   dimension of q-point set
  !> @param[in] kSetShift: shift of k-point set
  !> @param[in] addQs:     additional phonon wavevectors q, stored sequentially in one-dimensional array
  !> @param[in] noPtsCon:  number of points to be checked in continuity tests of densities and potentials
  !> @param[in] pathDir:   (non-normed) direction vector of path through unit cell for path potential-gradient test
  !> @param[in] harSw:     Exchange-correlation potential switch for path potential-gradient test
  !> @param[in] extSw:     Exchange-correlation potential switch for path potential-gradient test
  !> @param[in] xcSw:      Exchange-correlation potential switch for path potential-gradient test
  !> @param[in] writeKpqArraySw: Switch if kpq array is written
  !> @param[in] logFUnit:  normed unit number for juPhon.log logfile.
  !>
  !>-------------------------------------------------------------------------------------------------------------------------------------
    subroutine InpParam2Log(kModeSw, kSetDim, qSetDim, kSetShift, addQs, noPtsCon, pathDir, harSw, extSw, xcSw, writeKpqArraySw, &
        & logFUnit)

      implicit none

      !Scalar parameter
      logical, intent(in) :: kModeSw
      integer, intent(in) :: noPtsCon
      logical, intent(in) :: harSw
      logical, intent(in) :: extSw
      logical, intent(in) :: xcSw
      logical, intent(in) :: writeKpqArraySw
      integer, intent(in) :: logFUnit

      !Type parameter
      real,    intent(in) :: pathDir(3)
      integer, intent(in) :: kSetDim(:)
      integer, intent(in) :: qSetDim(:)
      real,    intent(in) :: kSetShift(:)
      real,    intent(in) :: addQs(:)

      ! Scalar variables
      integer             :: nrAddQs ! number of additional q-wave vectors
      integer             :: ii ! loop index

      write(logFUnit, * )
      write(logFUnit, '(a)') 'Parameter from juPhon.inp'
      write(logFUnit, '(a)') '========================='
      write(logFUnit, * )


      write(logFUnit, '(a)')   'k- and q-point set'
      write(logFUnit, '(a)')   '------------------'
      write(logFUnit, '(a,3(i2,1x,a1))') 'Dimension of k-point set: ', kSetDim(1), 'x', kSetDim(2), 'x',kSetDim(3)
      write(logFUnit, '(a,3(i2,1x,a1))') 'Dimension of q-point set: ', qSetDim(1), 'x', qSetDim(2), 'x',qSetDim(3)
      if ( any( addQs > 0 ) ) then
        nrAddQs = size( addQs ) / 3
        do ii = 1, nrAddQs
          write(logFUnit, '(a,3(f5.3,a1,1x),a1)') 'Additional q-point: (', addQs(ii * 1), ',', addQs(ii * 2), ',', addQs(ii * 3), ')'
        end do
      end if
      write(logFUnit, '(a,a1,3(f5.3,a1,1x),a)') 'Shift     of k-point set: ', '(', kSetShift(1), ',', kSetShift(2), ',',kSetShift(3), ')'
      write (logFUnit, '(a)') 'SHIFT NOT IMPLEMENTED YET!'

      if ( .not.kModeSw ) then

#ifdef DEBUG_MODE
        write(logFUnit, *)
        write(logFUnit, '(a)')   'Debug parameter'
        write(logFUnit, '(a)')   '---------------'
        write(logFUnit, '(a)') 'Potential and density continuity check:'
        write(logFUnit, '(2x,a,i3)') 'Number of test points: ', noPtsCon
        write(logFUnit, *)
        write(logFUnit, '(a)') 'Test of unperturbed potential gradient:'
        write(logFUnit, '(2x,a,a1,3(f5.2,a1,1x),a1)') 'Direction of path: ', '(', pathDir(1), ',', pathDir(2), ',', pathDir(3), ')'
        if (harSw) write(logFUnit, '(2x,a)') 'Hartree potential activated'
        if (extSw) write(logFUnit, '(2x,a)') 'External potential activated'
        if (xcSw)  write(logFUnit, '(2x,a)') 'xα exchange-correlation potential activated'
        write(logFUnit, *)
        if ( writeKpqArraySw ) then
          write(logFUnit, '(a)') 'Mapping array k + q -> k is written into file kpqMapArray'
        else
          write(logFUnit, '(a)') 'Mapping array k + q -> k is not written out to file'
        end if
#endif
      end if

    end subroutine InpParam2Log


  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Writes information read in and generated by fleur init.
  !>
  !> @details
  !> 
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine logFleurInit(atoms, input, dimens, cell, lathar, stars, sym, xcpot, logUnit)

    use m_types
    use m_JPConstants, only : namat_const

    implicit none

    ! Type parameters
    type(t_atoms),      intent(in)  :: atoms
    type(t_input),      intent(in)  :: input
    type(t_dimension),  intent(in)  :: dimens
    type(t_cell),       intent(in)  :: cell
    type(t_sphhar),     intent(in)  :: lathar
    type(t_stars),      intent(in)  :: stars
    type(t_sym),        intent(in)  :: sym
    type(t_xcpot),      intent(in)  :: xcpot

    ! Scalar parameters
    integer,            intent(in)  :: logUnit

    ! Scalar variables
    integer                         :: ii
    integer                         :: jj
    integer                         :: m
    integer                         :: iatom
    integer                         :: itype
    integer                         :: ieqat
    integer                         :: isym
    integer                         :: ilh

    write ( logUnit, * )
    write ( logUnit, '(a)') 'Relevant parameters after calling Fleur initialization routine'
    write ( logUnit, '(a)') '=============================================================='
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Atom-related information'
    write ( logUnit, '(a)' ) '------------------------'
    write ( logUnit, '(a,1x,i4)' ) 'Total number of atoms      :', atoms%nat
    write ( logUnit, '(a,1x,i4,/)' ) 'Total number of atom types :', atoms%ntype
    write ( logUnit, '(2x,"index",3x,"type",3x,"name",3x,"atomic number",10x,"cartesian coordinates",15x,"internal coordinates",/)' )
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        write ( logUnit, '(3x,i4,4x,i3,5x,a,8x,f8.4,4x,2("(",3(f9.5,1x),")",4x))' ) iatom, itype, &
          & namat_const(nint( atoms%zatom(itype ))), atoms%zatom(itype ), (atoms%pos(ii,iatom),ii=1,3), (atoms%taual(ii,iatom),ii=1,3)
      end do
    end do
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a,1x,f8.4)' ) 'Number of valence electrons          :', input%zelec
    write ( logUnit, '(a,1x,i4)' )   'Maximal number of spin values        :', input%jspins
    write ( logUnit, '(a,1x,i4)' )   'Maximal orbital quantum number       :', atoms%lmaxd
    write ( logUnit, '(a,1x,i4,/)' ) 'Maximal number of radial grid points :', atoms%jmtd
    write ( logUnit, '(2x,"type",3x,"name",3x,"MT radius",3x,"MT volume",3x,"MT grid points",5x,"MT dx",3x,"lmax",3x,"non-spheric. &
      &lmax",3x,"core states",3x,"number of LOs",/)')
    do itype = 1, atoms%ntype
      write ( logUnit, '(3x,i3,5x,a,4x,f8.5,4x,f8.5,8x,i4,8x,f8.5,3x,i2,11x,i2,15x,i2,12x,i3)' ) itype, &
        & namat_const(nint( atoms%zatom(itype ))), atoms%rmt(itype), atoms%volmts(itype), atoms%jri(itype), atoms%dx(itype),&
        & atoms%lmax(itype), atoms%lnonsph(itype), atoms%ncst(itype), atoms%nlo(itype)
    end do

    write (logUnit, *)
    if ( any( atoms%nlo /= 0) ) then
      write ( logUnit, * )
      write ( logUnit, '(a,1x,i8)' )   'Maximal LO orbital quantum number l  :', atoms%llod
      write ( logUnit, '(a,1x,i8,/)' ) 'Total number of Gs for local orbitals:', atoms%nlotot
      write ( logUnit, '(2x,"type",3x,"name",3x,"number of local orbital(s),",3x, "their orbital quantum number l"/)' )
      do itype = 1, atoms%ntype
        write ( logUnit, '(3x,i3,5x,a,14x,i3,7x)',advance='no' ) itype, namat_const(nint( atoms%zatom(itype ))), atoms%nlo(itype)
        write ( logUnit, *) atoms%llo(:, itype)
      end do
    end if
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a)' )   'Unit cell information'
    write ( logUnit, '(a)' )   '---------------------'
    write ( logUnit, '(a,a)' ) 'Latice name: ', cell%latnam
    write ( logUnit, '(a)' )   'Bravais matrix of cartesian lattice vectors:'
    do ii = 1, 3
      write ( logUnit, '(5x,3(f10.6))' ) ( cell%amat(ii,jj), jj=1, 3 )
    enddo
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Transposed Bravais matrix of reciprocal lattice vectors:'
    do ii = 1, 3
      write ( logUnit, '(5x,3(f10.6))' ) ( cell%bmat(ii,jj), jj=1, 3 )
    enddo
    write ( logUnit, * )
    write ( logUnit, '(a,1x,f15.8)' ) 'Volume of direct unit cell     =', cell%vol
    write ( logUnit, '(a,1x,f15.8)' ) 'Volume of reciprocal unit cell =', cell%omtil
    write ( logUnit, '(a,1x,f15.8)' ) 'Volume of interstitial region  =', cell%volint
    write ( logUnit, '(a,1x,f15.8)' ) 'Volume of all muffin-tins      =', cell%vol - cell%volint
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Symmetry, stars and lattice harmonic information'
    write ( logUnit, '(a)' ) '------------------------------------------------'
    write ( logUnit, '(a,1x,a)' )  'Name of space group                 :', sym%namgrp
    write ( logUnit, '(a,1x,i2)' ) 'Number of found symmetry operations =', sym%nop
    if ( sym%symor ) then
      write ( logUnit, '(2x,a,/)' ) 'All symmetry operations are symmorphic!'
    else
      write ( logUnit, '(2x,a,/)' ) 'There are non-symmorphic operations!'
    end if
    if ( sym%invs ) then
      write ( logUnit, '(a)' ) 'The system shows inversion symmetry!'
      do iatom = 1, atoms%nat
        if ( atoms%invsat(iatom) == 0 ) then
          write ( logUnit, '(2x,a,1x,i3,a,i4)' ) 'atom ', iatom, ': Not mappable to other atom by inversion symmetry! &
          & Atom index (starting from 1) in atom mapping array: ', sym%invsatnr(iatom)
        else if ( atoms%invsat(iatom) == 1 ) then
          write ( logUnit, '(2x,a,1x,i3,a,i4)' ) 'atom ', iatom, ': Parent atom (can be mapped FROM by inversion symmetry)! &
          & Atom index (starting from 1) in atom mapping array: ', sym%invsatnr(iatom)
        else if ( atoms%invsat(iatom) == 2 ) then
          write ( logUnit, '(2x,a,1x,i3,a,i4)' ) 'atom ', iatom, ': Child atom (can be mapped TO by inversion symmetry)! &
          & Atom index (starting from 1) in atom mapping array: ', sym%invsatnr(iatom)
        end if
      end do
    else
      write ( logUnit, '(a)' ) 'There is no inversion symmetry within the system!'
    end if
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a,1x,f12.5)' ) 'Reciprocal cutoff value for interstitial plane waves (kmax)    =', input%rkmax
    write ( logUnit, '(a,1x,i6)')    'Number of stars                                                =', stars%ng3
    write ( logUnit, '(a,1x,f12.5)' ) 'Maximal length of stars (gmax)                                 =', stars%gmax
    write ( logUnit, '(a,1x,f12.5)' ) 'Maximal length of stars for xc-potential (gmaxxc)              =', xcpot%gmaxxc
    write ( logUnit, '(a,1x,i3,1x,"x",1x,i3,1x,"x",i3)' ) 'Dimensions used for Fast-Fourier-Transformation (k1d,k2d,k3d)  =',stars%k1d, stars%k2d, stars%k3d
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a,1x,i2)' ) 'Maximal index of used site symmetries             =', lathar%ntypsd
    write ( logUnit, '(a,1x,i3)' ) 'Maximal number of lattice harmonics for all atoms =', lathar%nlhd + 1
    write ( logUnit, '(a,1x,i3)' ) 'Maximal number of members of sphhar               =', lathar%memd
    write ( logUnit, '(2x,"index",3x,"type",3x,"name",3x,"site symmetry index",3x,"number of lattice harmonics",/)' )
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        write ( logUnit, '(3x,i4,4x,i3,5x,a,10x,i4,22x,i4)' ) iatom, itype, namat_const(nint( atoms%zatom(itype ))), &
          & atoms%ntypsy(iatom), lathar%nlh(atoms%ntypsy(iatom)) + 1
      end do
    end do

#ifdef DEBUG_MODE
    write ( logUnit, * )
    do isym = 1, sym%nsymt
      write ( logUnit, '(/,"Local symmetry",i3,":",i4," lattice harmonics ")' ) isym, lathar%nlh(isym) + 1
      do ilh = 0, lathar%nlh(isym)
        if ( lathar%nmem( ilh, isym) > 1 ) then
          write ( logUnit,'(/,2x,"lattice harmonic",i4,":  l=",i2,",",i3," members:")' ) ilh+1, lathar%llh(ilh,isym), &
            & lathar%nmem(ilh,isym)
        else
          write ( logUnit,'(/,2x,"lattice harmonic",i4,":  l=",i2,",",i3," member:")' ) ilh+1, lathar%llh(ilh,isym), &
            & lathar%nmem(ilh,isym)
        end if
        if ( mod( lathar%nmem(ilh,isym), 2 ) == 1 ) then
          write ( logUnit,'(2x,i5,2f14.8,5x,i5,2f14.8)' ) lathar%mlh(1, ilh, isym), lathar%clnu( 1, ilh, isym )
          if ( lathar%nmem(ilh, isym) > 1 ) then
            write ( logUnit,'(2x,i5,2f14.8,5x,i5,2f14.8)' ) &
              & (lathar%mlh(m, ilh, isym), lathar%clnu(m, ilh,isym), m=2, lathar%nmem(ilh, isym))
          endif
        else
          write ( logUnit, '(2x,i5,2f14.8,5x,i5,2f14.8)' ) &
            & (lathar%mlh(m,ilh,isym), lathar%clnu(m, ilh, isym), m=1, lathar%nmem(ilh, isym))
        end if
      end do
    end do
#endif
    write ( logUnit, * )
    write ( logUnit, * )
    write ( logUnit, '(a)' )       'Information on internal dimensions'
    write ( logUnit, '(a)' )       '----------------------------------'
    write ( logUnit, '(a,1x,i6)' ) 'Number of points for Lebedev quadrature (nspd)                     :', dimens%nspd
    write ( logUnit, '(a,1x,i6)' ) 'Number of basis functions (nvd)                                    :', dimens%nvd
    write ( logUnit, '(a,1x,i6)' ) 'Maximal number of eigenvalues for all k-points (neigd)             :', dimens%neigd
    write ( logUnit, '(a,1x,i6)' ) 'Maximal Weinert parameter (ncvd)                                   :', dimens%ncvd
    write ( logUnit, '(a,1x,i6)' ) 'Maximal value of lm combined index (lmd)                           :', dimens%lmd
    write ( logUnit, '(a,1x,i6)' ) 'Maximal value of combined index lm and lm (lmplm)                  :', dimens%lmplmd
    write ( logUnit, '(a,1x,i6)' ) 'Total number of basis functions including local orbitals (nbasfcn) :', dimens%nbasfcn
    write ( logUnit, * )
    write ( logUnit, * )

  end subroutine logFleurInit

  !>--------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Writes readable information from eig file to log file for testing reason.
  !>
  !> @details
  !> 
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine LogReadFromEigFile(atoms, kpts, dimens, El, nv, ne, logUnit)

    use m_types
    use m_JPConstants

    implicit none

    type(t_atoms),     intent(in) :: atoms
    type(t_kpts),      intent(in) :: kpts
    type(t_dimension), intent(in) :: dimens
    real,              intent(in) :: El(:, 0:, :, :)
    integer,           intent(in) :: nv(:, :)
    integer,           intent(in) :: ne(:)
    integer,           intent(in) :: logUnit

    integer                       :: itype
    integer                       :: oqn_l
    integer                       :: ilo
    integer                       :: ikpt

    write ( logUnit, '(a)' ) 'Read in Fleur eig file'
    write ( logUnit, '(a)' ) '======================'
    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Energy parameters of converged Fleur calculation'
    write ( logUnit, '(a,/)' ) '------------------------------------------------'
    write ( logUnit, '(2x,"type",3x,"name",3x,"LO",3x,"quantum number l",9x,"E_l",/)' )
    do itype = 1, atoms%ntype
      do oqn_l = 0, atoms%lmax(itype)
        write ( logUnit, '(3x,i3,5x,a,3x,a,9x,i2,8x,f15.8)' ) itype, namat_const(nint( atoms%zatom(itype ))), 'no', oqn_l, El(1, oqn_l, itype, 1)
      end do
      do ilo = 1, atoms%nlo(itype)
        write ( logUnit, '(3x,i3,5x,a,3x,a,9x,i2,8x,f15.8)' ) itype, namat_const(nint( atoms%zatom(itype ))), 'yes', atoms%llo(ilo, itype), El(ilo, atoms%llo(ilo, itype), itype, 1)
      end do
    end do

  end subroutine LogReadFromEigFile

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Logs overview of number of basis function, number of all bands and occupied bands per k-point
  !>
  !> @details
  !> Partially this information is already known after reading the eigfile so it should be given out after reading the eig file. But for
  !> sake of readability, it was more logical to accumulate coherent data
  !> 
  !>--------------------------------------------------------------------------------------------------------------------------------------
  !todo Efermie ausgeben und nlotot ausgeben und summe von nvd und nlotot ausgeben
  subroutine LogBandsBFcn( kpts, logUnit, nv, ne, nobd )

    use m_types

    implicit none

    ! Type parameter
    type(t_kpts), intent(in) :: kpts

    ! Scalar parameter
    integer,      intent(in) :: logUnit

    ! Array parameters
    integer,      intent(in) :: nv(:, :)
    integer,      intent(in) :: ne(:)
    integer,      intent(in) :: nobd(:, :)

    ! Scalar variable
    integer                  :: ikpt

    write ( logUnit, * )
    write ( logUnit, '(a)' ) 'Number of regular basis vectors (nv) and eigenvalues (ne) per k-point'
    write ( logUnit, '(a,/)' ) '---------------------------------------------------------------------'
    write ( logUnit, '("k-point index",17x,"k-point",23x,"nv",9x, "ne",5x,"nobd",/)' )
    do ikpt = 1, kpts%nkpt
      write( logUnit, '(4x,i3,8x,"(",1x,3(f10.8,1x),")",5x,i6,5x,i6,5x,i4)' ) ikpt, kpts%bk(:, ikpt), nv(1, ikpt), ne(ikpt), &
        &nobd(ikpt, 1)
    end do

  end subroutine LogBandsBFcn

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Writes log after reading in cnd1 file from Fleur
  !>
  !> @details
  !> nradFun is tested by printing out all orbitals
  !> 
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine LogReadcdn1( atoms, qtotIR, qMTs, qtot, logUnit )

    use m_types
    use m_JPConstants, only : namat_const

    implicit none

    ! Type parameter
    type(t_atoms), intent(in) :: atoms

    ! Scalar parameter
    real,          intent(in) :: qtotIR
    real,          intent(in) :: qtot
    integer,       intent(in) :: logUnit

    ! Array parameter
    real,          intent(in) :: qMTs(:)

    ! Scalar variable
    integer                   :: itype

    write ( logUnit, '(a)' ) 'Read in converged and unperturbed density for intersitital and muffin-tin regions from file cdn1'
    write ( logUnit, '(a)' ) '================================================================================================'
    write ( logUnit, * )
    write ( logUnit, '(a,1x,f12.6)' )   'Total charge in interstitial region :', qtotIR
    write ( logUnit, '(a,1x,f12.6)' )   'Total charge in muffin-tin regions  :', sum(qMts)
    write ( logUnit, '(a,1x,f12.6,/)' ) 'Total charge in unit cell           :', qtot
    write ( logUnit, '(2x,"type",3x,"name",5x,"MT charge",3x,"equivalent atoms",2x,"MT charge per atom",/)' )
    do itype = 1, atoms%ntype
      write ( logUnit, '(3x,i3,5x,a,2x,f12.6,16x,i3,8x,f12.6)' ) itype, namat_const(nint( atoms%zatom(itype ))), qMts(itype), &
        & atoms%neq(itype), qMts(itype) / atoms%neq(itype)
    end do

  end subroutine LogReadcdn1

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Writes log after calculating the radial functions
  !>
  !> @details
  !> nradFun is tested by printing out all orbitals
  !> 
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine LogRadSolutions( atoms, usdus, nRadFun, iloTable, uuilon, duilon, ulouilopn, logUnit )

    use m_types
    use m_JPConstants, only : namat_const

    implicit none

    ! Type parameters
    type(t_atoms), intent(in) :: atoms
    type(t_usdus), intent(in) :: usdus

    ! Scalar parameters
    integer,       intent(in) :: logUnit

    ! Array parameters
    integer,       intent(in) :: nRadFun(0:, :) ! number of radial Functions, i.e. u, udot, uLO for all LOs per l and itype
    integer,       intent(in) :: iloTable(:, 0:, :)
    real,          intent(in) :: uuilon(:, :)
    real,          intent(in) :: duilon(:, :)
    real,          intent(in) :: ulouilopn(:, :, :)

    ! Scalar variables
    integer                   :: itype
    integer                   :: ilo, jlo
    integer                   :: p
    integer                   :: oqn_l


    write( logUnit, * )
    write( logUnit, * )
    write( logUnit, '(a)' ) 'Generate solutions of radial MT problem: u, udot, and uLOs'
    write( logUnit, '(a,/)' ) '=========================================================='

    write ( logUnit, '(2x,"type",3x,"name",3x, "l",10x,"u(R)",11x, "du(R)/dr", 10x,"du/dE",11x,"d2u/drdE",6x, "<du/dE|du/dE>",/)' )
    do itype = 1, atoms%ntype
      do oqn_l = 0, atoms%lmax(itype)
        write ( logUnit, '(3x,i3,5x,a,2x,i2,2x,5(f15.8,2x))' ) itype, namat_const(nint( atoms%zatom(itype ))), oqn_l, &
          & usdus%us(oqn_l, itype, 1), usdus%dus(oqn_l, itype, 1), usdus%uds(oqn_l, itype, 1), usdus%duds(oqn_l, itype, 1), &
          & usdus%ddn(oqn_l, itype, 1)
      end do
    end do

    if ( any( nRadFun > 2 ) ) then
      write ( logUnit, '(2x,"type",3x,"name",3x, "l", "number of LO","u_LO(R)", "du_LO(R)/dr","<u_LO|u>","<u_LO|du/dE>",3x,"uuilon",3x,"duilon",/)' )
      do itype = 1, atoms%ntype
        do oqn_l = 0, atoms%lmax(itype)
          do p = 3, nRadFun(oqn_l, itype)
            write( logUnit, '(3x,i3,5x,a,2x,i2,1x,i2,2x,6(f15.8,2x))' ) itype, namat_const(nint( atoms%zatom(itype ))), oqn_l, &
              & iloTable(p, oqn_l, itype ), usdus%ulos(iloTable(p, oqn_l, itype), itype, 1), &
              & usdus%dulos(iloTable(p, oqn_l, itype), itype, 1), usdus%uulon(iloTable(p, oqn_l, itype), itype, 1), &
              & usdus%dulon(iloTable(p, oqn_l, itype), itype, 1), uuilon(iloTable(p, oqn_l, itype), itype), duilon(iloTable(p, oqn_l, itype), itype)
          end do
        end do
      end do
    end if

    if ( any( nRadFun > 3 ) ) then
      write ( logUnit, '(2x,"type",3x,"name",3x, "l", "number of LO", "number of LOprime",3x, "<u_LO|u_LOprime>",3x,"ulouilopn",/)' )
      do itype = 1, atoms%ntype
        do ilo = 1, atoms%nlo(itype)
          do jlo = 1, atoms%nlo(itype)
             write ( logUnit, '(3x,i3,5x,a,i3,2x,i3,2x,i3,2x,i3,2x,f15.8,2x,f15.8)' ) itype, namat_const(nint( atoms%zatom(itype ))), ilo, &
              &atoms%llo(ilo, itype), jlo, atoms%llo(jlo, itype), usdus%uloulopn(ilo, jlo, itype, 1), ulouilopn(ilo, jlo, itype)
          end do
        end do
      end do
    end if

  end subroutine logRadSolutions

  !?--------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Writes log after back-rotating local coordinate systems and calculating coefficients for lattice harmonics
  !>
  !> @details
  !> nradFun is tested by printing out all orbitals
  !> 
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine LogRotLh( atoms, lathar, memd_atom, logUnit, clnu_atom, nmem_atom, mlh_atom )

    use m_types
    use m_jPConstants, only : namat_const

    implicit none

    ! Type parameters
    type(t_atoms),  intent(in) :: atoms
    type(t_sphhar), intent(in) :: lathar

    ! Scalar parameters
    integer,        intent(in) :: memd_atom
    integer,        intent(in) :: logUnit

    ! Array parameters
    integer,        intent(in) :: mlh_atom(:,0:,:)
    integer,        intent(in) :: nmem_atom(0:, :)
    complex,        intent(in) :: clnu_atom(:,0:,:)

    ! Scalar variables
    integer                    :: iatom
    integer                    :: itype
    integer                    :: ieqat
    integer                    :: itypsym
    integer                    :: ilh
    integer                    :: imem


    write (logUnit, * )
    write (logUnit, * )
    write (logUnit, '(a)') 'Created information to use lattice-harmonic expanded quantities for non-symmetric juPhon Calculation'
    write (logUnit, '(a)') '===================================================================================================='
    write (logUnit, * )


    write ( logUnit, '(2x,"index",3x,"type",3x,"name",3x,"site symmetry index",3x,"lattice harmonic",3x, "magnetic quantum number", &
      &13x,"phase clnu",/)' )
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        itypsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(itypsym)
          do imem = 1, nmem_atom(ilh, itypsym)
            write( logUnit, '(3x,i4,4x,i3,5x,a,10x,i4,17x,i2,18x,i3,12x,2(f15.8,1x))' ) iatom, itype, &
              & namat_const(nint( atoms%zatom(itype ))), itypsym, ilh, mlh_atom(imem, ilh, iatom), clnu_atom(imem, ilh, iatom)
          end do
        end do
      end do
    end do

  end subroutine LogRotLh

  !>--------------------------------------------------------------------------------------------------------------------------------------
  !> @author 
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Finishes logfile
  !>
  !> @details
  !> 
  !>--------------------------------------------------------------------------------------------------------------------------------------
  subroutine finishLogFile( startTime, logFUnit )
    implicit none

    ! Scalar parameter
    integer, intent (in)  :: logFUnit

    ! Array parameter
    integer               :: startTime(8)

    ! Local variables
    integer               :: dateTime(8)
    character(len=30)     :: hostName

    ! Find out starting date and time
    call date_and_time (VALUES=dateTime)

    ! Init logfile
    write (logFUnit, * )
    write (logFUnit, * )
    write (logFUnit, '(a)')   'juPhon calculation successfully completed!'
    write (logFUnit, '(a,/)') '******************************************'
    write (logFUnit, '(a16,a2,a1,a2,a1,i4,a4,a2,a1,a2,a1,a2,a1,a3)') 'Date and Time : ', zero2l(dateTime(3)), '.', &
    & zero2l(dateTime(2)), '.', dateTime(1), ' at ', zero2l(dateTime(5)), ':', zero2l(dateTime(6)), ':', zero2l(dateTime(7)), &
    & ':', zero3l(dateTime(8))
    write (logFUnit, '(a,i2,1x,a,2x,i2,1x,a,2x,i2,1x,a,2x,i3,1x,a,2x,i4,1x,a)')   'Duration      : ', abs(dateTime(3) - startTime(3)), &
      &'d', abs(dateTime(5) -startTime(5)), 'h', abs(dateTime(6) - startTime(6)), 'min', abs(dateTime(7) - startTime(7)), 's', &
      & abs(dateTime(8) - startTime(8)), 'ms'

  end subroutine finishLogFile

end module m_jpLog
