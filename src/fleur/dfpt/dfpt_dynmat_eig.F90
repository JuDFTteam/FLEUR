!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_dfpt_dynmat_eig
  use m_juDFT_stop
  use m_constants
  use m_xmlOutput
  use m_types_atoms
  implicit none
  private
  public :: diagonalizedynmat, calculatefrequencies, hermitize_scale_dynmat, apply_sumrule

  contains

  subroutine hermitize_scale_dynmat(atoms, dynMat)
    !! Hermitizes the dynamical matrix in place and divides it by the square root of the atomic masses.

    type(t_atoms), intent(in)    :: atoms
    complex,       intent(inout) :: dynMat(:, :)

    integer :: ii, jj
    real    :: atomic_mass_array(118)

    atomic_mass_array = atomicMasses_const * massInElectronMasses
    dynMat = (dynMat + conjg(transpose(dynMat)))/2.0
    do jj = 1, size(dynMat, 2)
      do ii = 1, size(dynMat, 1)
        dynMat(ii, jj) = dynMat(ii, jj)/sqrt(atomic_mass_array(atoms%nz(ceiling(jj/3.0)))*atomic_mass_array(atoms%nz(ceiling(ii/3.0))))
      end do
    end do

  end subroutine hermitize_scale_dynmat

  subroutine apply_sumrule(dynMat, dynMatGamma)
    !! Imposes the acoustic sum rule by projecting the three acoustic modes of the Gamma matrix out of dynMat.

    complex, intent(inout) :: dynMat(:, :)
    complex, intent(in)    :: dynMatGamma(:, :)

    integer :: n, lwork, info, ii, jj, iDir

    real,    allocatable :: w(:), rwork(:)
    complex, allocatable :: work(:), dynMat0(:, :)

    n = size(dynMat, 1)
    lwork = 2 * n
    allocate(w(n), rwork(3 * n), work(lwork))
    w = 0.
    rwork = 0.
    work = 0.
    dynMat0 = dynMatGamma

    call zheev('V', 'U', n, dynMat0, n, w, work, lwork, rwork, info)
    if (info /= 0) call juDFT_error("Diagonalization of the Gamma dynamical matrix failed.", calledby="apply_sumrule")

    do jj = 1, n
      do ii = 1, n
        do iDir = 1, 3
          dynMat(jj, ii) = dynMat(jj, ii) - w(iDir)*dynMat0(jj, iDir)*conjg(dynMat0(ii, iDir))
        end do
      end do
    end do

  end subroutine apply_sumrule

  subroutine DiagonalizeDynMat(atoms, qvec, dynMat, w, a, l_writeOutput)
    !! Diagonalizes the dynamical matrix; w are the eigenvalues, a the eigenvectors.

    type(t_atoms),                  intent(in)  :: atoms
    real,                           intent(in)  :: qvec(3)
    complex,                        intent(in)  :: dynMat(:, :)
    real,              allocatable, intent(out) :: w(:)
    complex,           allocatable, intent(out) :: a(:, :)
    logical,                        intent(in)  :: l_writeOutput

    integer :: n, lwork, info
    integer :: ii, itype, ieqat, iatom

    real,    allocatable :: rwork(:)
    complex, allocatable :: work(:)

    n = size(dynMat, 1)
    a = dynMat

    if (l_writeOutput) then
      write(*, '(a,3f9.3)') 'q =', qvec
      write(*, '(a)')       '==================================='
      write(*, '(a)')
      write(*, '(a)') 'Original Dynamical Matrix [mass corrected]'
      do ii = 1, n
        write(*, '(3(2(es16.8,1x),3x))') a(ii, :)
      end do
    end if

    lwork = 2 * n
    allocate(w(n), rwork(3 * n), work(lwork))
    w = 0.
    rwork = 0.
    work = 0.
    call zheev('V', 'U', n, a, n, w, work, lwork, rwork, info)

    if (info < 0) then
      write(*, *) 'The ', info, '(st/nd/)th argument had an illegal value'
      call juDFT_error("Illegal argument value.", calledby="DiagonalizeDynMat")
    else if (info > 0) then
      write(*, *) 'The diagonalization algorithm failed to converge; ', info, ' off-diagonal elements of an intermediate tridiagonal form did not converge to zero.'
      call juDFT_error("Diagonalization failed.", calledby="DiagonalizeDynMat")
    end if

    if (l_writeOutput) then
      write(*, '(a)') 'The eigenvalues of the Dynamical matrix are:'
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          write(*, "(a,i2,a,1x,3(es16.8,1x),',',5x)") 'Atom', iatom, ':', w((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
        end do
      end do
      do ii = 1, n
        write(*, '(3(2(es16.8,1x),3x))') a(ii, :)
        write(*, *)
      end do
    end if

  end subroutine diagonalizeDynMat

  subroutine CalculateFrequencies( atoms, iqpt, eigenVals, eigenFreqs, qvec ,l_writeOutput)

    implicit none

    ! Type parameter
    type(t_atoms),              intent(in)  :: atoms

    ! Scalar parameter
    integer,                    intent(in)  :: iqpt

    ! Array parameter
    real,                       intent(in)  :: eigenVals(:)
    complex,          allocatable, intent(out) :: eigenFreqs(:)
    real,             intent(in) :: qvec(3)
    logical,          intent(in) :: l_writeOutput

    ! Scalar variables
    integer                                 :: itype
    integer                                 :: ieqat
    integer                                 :: iatom
    integer                                 :: idir
    real                                    :: convFact
    real                                    :: THz_to_cm

    ! Array variables
    CHARACTER(LEN=20)                         :: attributes(5)

    allocate(eigenFreqs(3*atoms%nat))
    eigenFreqs = 0.
    itype = 1


    convFact = 4.35974472220e-18 / (5.2918e-11)**2 / 9.1093837015e-31 
    THz_to_cm = 33.3564

    if (l_writeOutput) then
      write(*, *)
      write(*, '(a)') 'Eigenfrequencies in THz'
    end if
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          !if (abs(eigenVals((iatom - 1) * 3 + idir)) < 1e-5) then
          !  eigenFreqs((iatom - 1) * 3 + idir) = cmplx(0., 0.)
          if (eigenVals((iatom - 1) * 3 + idir) < 0 ) then
            eigenFreqs((iatom - 1) * 3 + idir) = cmplx(0., -sqrt(abs(eigenVals((iatom - 1) * 3 + idir)) * convFact) / tpi_const * 1e-12)
          else
            eigenFreqs((iatom - 1) * 3 + idir) = sqrt(eigenVals((iatom - 1) * 3 + idir) * convFact) / tpi_const * 1e-12
          end if
        end do ! idir
        if (l_writeOutput) write(*, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype


    if (l_writeOutput) then
      write(*, *)
      write(*, '(a)') 'Eigenfrequencies in 1/cm'
    end if
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do idir = 1, 3
          eigenFreqs((iatom - 1) * 3 + idir) = eigenFreqs((iatom - 1) * 3 + idir) * THz_to_cm
        end do ! idir
        if (l_writeOutput) write(*, "(a,i2,a,1x,3(2es16.8,1x),',',5x)") 'Atom', iatom, ':', eigenFreqs((iatom - 1)* 3 + 1:(iatom - 1)* 3 + 3)
      end do ! ieqat
    end do ! itype

    ! write into out.xml
    CALL openXMLElementNoAttributes('phonon')
    CALL openXMLElementNoAttributes('dynMat')
    attributes = ''
    WRITE(attributes(1),'(i0)') iqpt
    WRITE(attributes(2),'(f0.8)') qvec(1)
    WRITE(attributes(3),'(f0.8)') qvec(2)
    WRITE(attributes(4),'(f0.8)') qvec(3)
    attributes(5) = '1/cm'
    CALL writeXMLElementPoly('frequencies',(/'iqpt', 'q_x ' , "q_y " , "q_z ",'unit'/), attributes, REAL(eigenFreqs)+AIMAG(eigenFreqs))
    CALL closeXMLElement('dynMat')
    CALL closeXMLElement('phonon')
  end subroutine CalculateFrequencies

end module m_dfpt_dynmat_eig
