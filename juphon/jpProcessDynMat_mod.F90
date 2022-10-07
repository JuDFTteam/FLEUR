module m_jpProcessDynMat

  use m_types

  implicit none

  contains

  subroutine DiagonalizeDynMat(atoms, qpts, calcEv, dynMat, w, a, iqpt)

    USE m_juDFT_stop
    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in)  :: atoms
    type(t_kpts),                  intent(in)  :: qpts

    ! Scalar parameters
    logical,                       intent(in)  :: calcEv

    ! Array parameters
    complex,                       intent(in)  :: dynMat(:, :)
    real,             allocatable, intent(out) :: w(:)
    complex, allocatable, intent(out) :: a(:, :)

    ! Array parameters

    ! Scalar variables
    ! jobVl : (LAPACK) switch for calculation of left eigenvector
    ! jobVr : (LAPACK) switch for calculation of right eigenvector
    ! n     : (LAPACK) order of the matrix a to diagonalize
    ! ldA   : (LAPACK) leading dimension of the array a
    ! ldVl  : (LAPACK) leading dimension of the array vL
    ! ldVr  : (LAPACK) leading dimension of the array vR
    ! lWork : (LAPACK) dimension of the array work or if equals -1 switch for LAPACK to determine optimal size of work array
    ! info  : (LAPACK) status of calculation
    character                                  :: jobz
    character                                  :: uplo
    integer                                    :: n
    integer                                    :: ldA
    integer                                    :: lWork
    integer                                    :: info
    integer                                    :: ii
    integer                                    :: jj
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: iatom
    integer, intent(in)                       :: iqpt

    ! Array variables
    ! a       : (LAPACK) matrix to diagonalize
    ! w       : (LAPACK) computed eigenvalues
    ! work    : (LAPACK) swap array for LAPACK diagonalization routine
    ! rwork   : (LAPACK) swap array for LAPACK diagonalization routine
    real,             allocatable              :: rwork(:)
    complex, allocatable              :: work(:)
    character(len=:), allocatable             :: filename
    character(len=11)                         :: filenameTemp

    if (iqpt < 10) then
      write(filenameTemp, '("dynMatq=00",i1)') iqpt
    else if (iqpt > 9 .and. iqpt < 100) then
      write(filenameTemp, '("dynMatq=0",i2)') iqpt
    else if (iqpt > 99 .and. iqpt < 1000) then
      write(filenameTemp, '("dynMatq=",i3)') iqpt
    end if
!    filename = trim(filenameTemp)
!    write(*, *) filename
    open( 109, file=filenameTemp, status='replace', action='write', form='formatted')

    ! Set parameter for LAPACK diagonalization routine
    if (calcEv ) then
      jobz = 'V'
    else
      jobz = 'N'
    end if

    ! Evaluate the lower triangle of the Dynamical matrix
    uplo = 'U'

    ! Dimensions of the Dynamical matrix
    n = 3 * atoms%nat
    lda = 3 * atoms%nat

    ! Copy Dynamical matrix into work array as it is destroyed if not calculating eigenvectors
    allocate( a( lda, n ) )
    do jj = 1, n
      do ii = 1, lda
      !todo is this hermitization correct, discuss with Markus
        !!!anfix
        a(ii, jj) = (dynMat(ii, jj) + conjg(dynMat(jj, ii)))/2.0
      end do
    end do

    !todo only works for one atom
    write(*, '(a,3f9.3)') 'q =', qpts%bk(1:3, iqpt)
    write(*, '(a)')       '==================================='
    write(*, '(a)')
    write(*, '(a)') 'Original Dynamical Matrix'
    DO ii = 1, lda
      write(*, '(3(2(es16.8,1x),3x))') a(ii, :)
    END DO

!    write(*, '(a)') 'Deviation from Hermiticity'
!    write(*, '(3(2(es16.8,1x),3x))') a(1, :) - dynMat(1, :)
!    write(*, '(3(2(es16.8,1x),3x))') a(2, :) - dynMat(2, :)
!    write(*, '(3(2(es16.8,1x),3x))') a(3, :) - dynMat(3, :)

    write(109, '(a,3f9.3)') 'q =', qpts%bk(1:3, iqpt)
    write(109, '(a)')       '==================================='
    write(109, '(a)')
    write(109, '(a)') 'Original Dynamical Matrix'
    DO ii = 1, lda
      write(*, '(3(2(es16.8,1x),3x))') a(ii, :)
    END DO

!    write(1000, '(a)') 'Deviation from Hermiticity'
!    write(1000, '(3(2(es16.8,1x),3x))') a(1, :) - dynMat(1, :)
!    write(1000, '(3(2(es16.8,1x),3x))') a(2, :) - dynMat(2, :)
!    write(1000, '(3(2(es16.8,1x),3x))') a(3, :) - dynMat(3, :)

    allocate( w(n))
    w = 0.

    allocate( rwork(3 * n))
    rwork = 0.

    lwork = 2 * n
    allocate(work(lwork))
    work = 0.

    ! Diagonalize dynamical matrix
    call zheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info )

    if ( info < 0 ) then
      write(*, *) 'The ', info, '(st/nd/)th argument had an illegal value'
      CALL juDFT_error("Illegal argument value.",calledby="DiagonalizeDynMat")
    else if ( info > 0 ) then
      write(*, *) 'The diagonalization algorithm failed to converge; ', info, ' off-diagonal elements of an intermediate tridiaonal form&
                                                                                                             & did not converge to zero.'
      CALL juDFT_error("Diagonalization failed.",calledby="DiagonalizeDynMat")
    end if

    write(*, '(a)') 'The eigenvalues of the Dynamical matrix are:'
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        write(*, "(a,i2,a,1x,3(es16.8,1x),',',5x)") 'Atom', iatom, ':', w((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        write(109, "(a,i2,a,1x,3(es16.8,1x),',',5x)") 'Atom', iatom, ':', w((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype

    if ( calcEv ) then
      !todo make the output nicer
      DO ii = 1, lda
        write(*, '(3(2(es16.8,1x),3x))') a(ii, :)
        write(*, *)
        write(109, '(3(2(es16.8,1x),3x))') a(ii, :)
        write(109, *)
      END DO
    else
      a(:, :) = cmplx(0., 0.)
    end if
    close( 109 )

  end subroutine diagonalizeDynMat

  subroutine CalculateFrequencies( atoms, iqpt, eigenVals, eigenFreqs )

    implicit none

    ! Type parameter
    type(t_atoms),              intent(in)  :: atoms

    ! Scalar parameter
    integer,                    intent(in)  :: iqpt

    ! Array parameter
    real,                       intent(in)  :: eigenVals(:)
    complex,          allocatable, intent(out) :: eigenFreqs(:)

    ! Scalar variables
    integer                                 :: itype
    integer                                 :: ieqat
    integer                                 :: iatom
    integer                                 :: idir
    real                                    :: massInElectronMasses
    real                                    :: convFact

    ! Array variables
    character(len=11)                         :: filenameTemp
    REAL                                      :: atomic_mass_array(118)

    if (iqpt < 10) then
      write(filenameTemp, '("dynMatq=00",i1)') iqpt
    else if (iqpt > 9 .and. iqpt < 100) then
      write(filenameTemp, '("dynMatq=0",i2)') iqpt
    else if (iqpt > 99 .and. iqpt < 1000) then
      write(filenameTemp, '("dynMatq=",i3)') iqpt
    end if
    open( 109, file=filenameTemp, status='old', action='write', form='formatted', position='append')

    allocate(eigenFreqs(3*atoms%nat))
    eigenFreqs = 0.
    itype = 1

    atomic_mass_array = [1.01, 4.00, 6.94, 9.01, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, &      ! up to neon
                      & 22.99, 24.31, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95, &                 ! up to argon
                      & 39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.85, 58.93, &          ! up to cobalt
                      & 58.69, 63.55, 65.38, 69.72, 72.63, 74.92, 78.97, 79.90, 83.80, &          ! up to krypton
                      & 85.47, 87.62, 88.91, 91.22, 92.91, 95.95, 97.40, 101.07, 102.91, &        ! up to ruthenium
                      & 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29, & ! up to xenon
                      & 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 146.00, 150.36, 151.96, & ! up to europium
                      & 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, & ! up to hafnium
                      & 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, & ! up to thallium
                      & 207.20, 208.98, 209.98, 210.00, 222.00, 223.00, 226.00, 227.00, 232.04, & ! up to thorium
                      & 231.04, 238.03, 237.00, 244.00, 243.00, 247.00, 247.00, 251.00, 252.00, & ! up to einsteinium
                      & 257.00, 258.00, 259.00, 262.00, 267.00, 269.00, 270.00, 272.00, 273.00, & ! up to hassium
                      & 277.00, 281.00, 281.00, 285.00, 286.00, 289.00, 288.00, 293.00, 294.00, 294.00]

    !if (atoms%nz(itype) == 10) then ! Neon
      !write(*, *) 'Mass for Neon'
      !massInElectronMasses = 20.18 * 1836.15 ! For Neon: 20.18 * 1836.15 m_e
    !else if (atoms%nz(itype) == 13) then ! Aluminium
      !write(*, *) 'Mass for Aluminium'
      !massInElectronMasses = 26.982 * 1836.15 ! For Aluminium : 26.982 * 1836.15 m_e
    !else if (atoms%nz(itype) == 18) then ! Argon
      !write(*, *) 'Mass for Argon'
      !massInElectronMasses = 39.948 * 1836.15 ! For Argon : 39.948 * 1836.15 m_e
    !else if (atoms%nz(itype) == 27) then ! Cobalt
      !write(*, *) 'Mass for Cobalt'
      !massInElectronMasses = 58.933 * 1836.15 ! For Cobalt : 58.933 * 1836.15 m_e
    !else if (atoms%nz(itype) == 29) then ! Copper
      !write(*, *) 'Mass for Copper'
      !massInElectronMasses = 63.546 * 1836.15 ! For Copper : 63.546 * 1836.15 m_e
    !else if (atoms%nz(itype) == 42) then ! Molybdenum
      !write(*, *) 'Mass for Molybdenum'
      !massInElectronMasses = 95.951 * 1836.15 ! For Molybdenum : 95.951 * 1836.15 m_e
    !else if (atoms%nz(itype) == 79) then ! Gold
      !write(*, *) 'Mass for Gold'
      !massInElectronMasses = 196.967 * 1836.15 ! For Gold : 196.967 * 1836.15 m_e
    !end if

    ! TODO: This will all no longer work for systems with different types of atoms;
    !       the mass prefactors need to be accounted for in the dynamical matrix
    !       for their respective displacements!
    massInElectronMasses = atomic_mass_array(atoms%nz(itype)) * 1836.15

    convFact = 4.35974472220e-18 / (5.2918e-11)**2 / 9.1093837015e-31 / massInElectronMasses

    write(*, *)
    write(*, '(a)') 'Eigenfrequencies in THz'
    write(109, *)
    write(109, '(a)') 'Eigenfrequencies in THz'
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          if (abs(eigenVals((iatom - 1) * 3 + idir)) < 1e-5) then
            eigenFreqs((iatom - 1) * 3 + idir) = cmplx(0., 0.)
          else if (eigenVals((iatom - 1) * 3 + idir) < 0 ) then
            eigenFreqs((iatom - 1) * 3 + idir) = cmplx(0., -sqrt(abs(eigenVals((iatom - 1) * 3 + idir)) * convFact) / tpi_const * 1e-12)
          else
            eigenFreqs((iatom - 1) * 3 + idir) = sqrt(eigenVals((iatom - 1) * 3 + idir) * convFact) / tpi_const * 1e-12
          end if
        end do ! idir
        write(*, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        write(109, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype


    write(*, *)
    write(*, '(a)') 'Eigenfrequencies in 1/cm'
    write(109, *)
    write(109, '(a)') 'Eigenfrequencies in 1/cm'
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          eigenFreqs((iatom - 1) * 3 + idir) = eigenFreqs((iatom - 1) * 3 + idir) * 33
        end do ! idir
        write(*, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        write(109, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype
    close( 109 )

  end subroutine CalculateFrequencies

end module m_jpProcessDynMat
