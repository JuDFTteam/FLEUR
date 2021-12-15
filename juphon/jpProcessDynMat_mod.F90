module m_jpProcessDynMat

  use m_types

  implicit none

  contains

  subroutine DiagonalizeDynMat(atoms, qpts, calcEv, dynMat, w, a, iqpt)

    use m_JPConstants, only : iu
    use mod_juPhonUtils, only : fopen, fclose

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
    call fopen( 1000, name=filenameTemp, status='replace', action='write', form='formatted')

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
    write(*, '(3(2(es16.8,1x),3x))') a(1, :)
    write(*, '(3(2(es16.8,1x),3x))') a(2, :)
    write(*, '(3(2(es16.8,1x),3x))') a(3, :)

!    write(*, '(a)') 'Deviation from Hermiticity'
!    write(*, '(3(2(es16.8,1x),3x))') a(1, :) - dynMat(1, :)
!    write(*, '(3(2(es16.8,1x),3x))') a(2, :) - dynMat(2, :)
!    write(*, '(3(2(es16.8,1x),3x))') a(3, :) - dynMat(3, :)

    write(1000, '(a,3f9.3)') 'q =', qpts%bk(1:3, iqpt)
    write(1000, '(a)')       '==================================='
    write(1000, '(a)')
    write(1000, '(a)') 'Original Dynamical Matrix'
    write(1000, '(3(2(es16.8,1x),3x))') a(1, :)
    write(1000, '(3(2(es16.8,1x),3x))') a(2, :)
    write(1000, '(3(2(es16.8,1x),3x))') a(3, :)

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
      NOstopNO
    else if ( info > 0 ) then
      write(*, *) 'The diagonalization algorithm failed to converge; ', info, ' off-diagonal elements of an intermediate tridiaonal form&
                                                                                                             & did not converge to zero.'
      NOstopNO
    end if

    write(*, '(a)') 'The eigenvalues of the Dynamical matrix are:'
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        write(*, "(a,i2,a,1x,3(es16.8,1x),',',5x)") 'Atom', iatom, ':', w((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        write(1000, "(a,i2,a,1x,3(es16.8,1x),',',5x)") 'Atom', iatom, ':', w((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype

    if ( calcEv ) then
      !todo make the output nicer
      write(*, '(3(2(es16.8,1x),3x))') a(1, :)
      write(*, *)
      write(*, '(3(2(es16.8,1x),3x))') a(2, :)
      write(*, *)
      write(*, '(3(2(es16.8,1x),3x))') a(3, :)
      write(*, *)
      write(1000, '(3(2(es16.8,1x),3x))') a(1, :)
      write(1000, *)
      write(1000, '(3(2(es16.8,1x),3x))') a(2, :)
      write(1000, *)
      write(1000, '(3(2(es16.8,1x),3x))') a(3, :)
      write(1000, *)
    else
      a(:, :) = cmplx(0., 0.)
    end if
    call fclose( 1000)

  end subroutine diagonalizeDynMat

  subroutine CalculateFrequencies( atoms, iqpt, eigenVals, eigenFreqs )

    use mod_juPhonUtils, only : fopen, fclose
    use m_JPConstants, only : tpi

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

    if (iqpt < 10) then
      write(filenameTemp, '("dynMatq=00",i1)') iqpt
    else if (iqpt > 9 .and. iqpt < 100) then
      write(filenameTemp, '("dynMatq=0",i2)') iqpt
    else if (iqpt > 99 .and. iqpt < 1000) then
      write(filenameTemp, '("dynMatq=",i3)') iqpt
    end if
    call fopen( 1000, name=filenameTemp, status='old', action='write', form='formatted', position='append')

    allocate(eigenFreqs(3*atoms%nat))
    eigenFreqs = 0.
    itype = 1

    if (atoms%nz(itype) == 10) then ! Neon
      write(*, *) 'Mass for Neon'
      massInElectronMasses = 20.18 * 1836.15 ! For Neon: 20.18 * 1836.15 m_e
    else if (atoms%nz(itype) == 13) then ! Aluminium
      write(*, *) 'Mass for Aluminium'
      massInElectronMasses = 26.982 * 1836.15 ! For Aluminium : 26.982 * 1836.15 m_e
    else if (atoms%nz(itype) == 18) then ! Argon
      write(*, *) 'Mass for Argon'
      massInElectronMasses = 39.948 * 1836.15 ! For Argon : 39.948 * 1836.15 m_e
    else if (atoms%nz(itype) == 27) then ! Cobalt
      write(*, *) 'Mass for Cobalt'
      massInElectronMasses = 58.933 * 1836.15 ! For Cobalt : 58.933 * 1836.15 m_e
    else if (atoms%nz(itype) == 29) then ! Copper
      write(*, *) 'Mass for Copper'
      massInElectronMasses = 63.546 * 1836.15 ! For Copper : 63.546 * 1836.15 m_e
    else if (atoms%nz(itype) == 42) then ! Molybdenum
      write(*, *) 'Mass for Molybdenum'
      massInElectronMasses = 95.951 * 1836.15 ! For Molybdenum : 95.951 * 1836.15 m_e
    else if (atoms%nz(itype) == 79) then ! Gold
      write(*, *) 'Mass for Gold'
      massInElectronMasses = 196.967 * 1836.15 ! For Gold : 196.967 * 1836.15 m_e
    end if

    convFact = 4.35974472220e-18 / (5.2918e-11)**2 / 9.1093837015e-31 / massInElectronMasses

    write(*, *)
    write(*, '(a)') 'Eigenfrequencies in THz'
    write(1000, *)
    write(1000, '(a)') 'Eigenfrequencies in THz'
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          if (abs(eigenVals((iatom - 1) * 3 + idir)) < 1e-5) then
            eigenFreqs((iatom - 1) * 3 + idir) = cmplx(0., 0.)
          else if (eigenVals((iatom - 1) * 3 + idir) < 0 ) then
            eigenFreqs((iatom - 1) * 3 + idir) = cmplx(0., -sqrt(abs(eigenVals((iatom - 1) * 3 + idir)) * convFact) / tpi * 1e-12)
          else
            eigenFreqs((iatom - 1) * 3 + idir) = sqrt(eigenVals((iatom - 1) * 3 + idir) * convFact) / tpi * 1e-12
          end if
        end do ! idir
        write(*, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        write(1000, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype


    write(*, *)
    write(*, '(a)') 'Eigenfrequencies in 1/cm'
    write(1000, *)
    write(1000, '(a)') 'Eigenfrequencies in 1/cm'
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          eigenFreqs((iatom - 1) * 3 + idir) = eigenFreqs((iatom - 1) * 3 + idir) * 33
        end do ! idir
        write(*, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        write(1000, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype
    call fclose( 1000)

  end subroutine CalculateFrequencies

end module m_jpProcessDynMat
!    write(formatString, '(i8)') 2 * 3 * atoms%nat
!    inquire(file='phonDispersion', exist=dispFile)
!    if (dispFile) then
!      call fopen(1000, name='phonDispersion', status='old', action='write', position='append', form='formatted')
!    else
!      call fopen(1000, name='phonDispersion', status='new', action='write', form='formatted')
!    end if
!    write(1000, '(i5,1x,3(es20.10,1x),1x,'// trim(formatString) // '(es20.10,1x))') iqpt, qpts%bk(:, iqpt), cmplx(wR, wI)
!    call fclose(1000)

!    allocate( eigenFreqs(3, atoms%nat) )
!    allocate( rEigenVecs(3, atoms%nat, 3, atoms%nat) )

!    ! Does such a structure make sense if eigenvalues are reordered?
!    iDatom = 0
!    do iDtype = 1, atoms%ntype
!      do iDeqat = 1, atoms%neq(iDtype)
!        iDatom = iatom + 1
!        do iDdir = 1, 3
!          iAlpha = iDdir + 3 * (iDatom - 1)
!          iatom = 0
!          do itype = 1, atoms%ntype
!            do ieqat = 1, atoms%neq(itype)
!              iatom = iatom + 1
!              ! Reformat the eigenvalues
!              eigenFreqs(iDdir, iDatom) = cmplx(wR(iAlpha), wI(iAlpha))
!
!              do idir = 1, 3
!                iBeta = idir + 3 * (iatom - 1)
!                ! Reformat the eigenvectors
!                if ( abs(wi(iAlpha)) > 1e-12 ) then
!                  ! We have a conjugate pair and are at the eigenvalue with positive imaginary part
!                  if ( (iAlpha < 3 * atoms%nat) .and. (wi(iAlpha + 1) > 1e-12) .and. (wi(iAlpha) > 0) .and. (wi(iAlpha + 1) < 0) ) then
!                    rEigenVecs(idir, iatom, iDdir, iDatom) = vR(iBeta, iAlpha) + iu * vR(iBeta, iAlpha)
!                  ! We have a conjugate pair and are at the eigenvalue with negative imaginary part
!                  else if( ( jj > 1 ) .and. ( wi(jj - 1) > 1e-12 ) .and. (wi(jj) < 0) .and. (wi(jj - 1) > 0) ) then
!                    rEigenVecs(idir, iatom, iDdir, iDatom) = vR(iBeta, iAlpha) - iu * vR(iBeta, iAlpha)
!                  end if
!                else
!                  rEigenVecs(idir, iatom, iDdir, iDatom) = vR(iBeta, iAlpha)
!                end if
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
