!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  F(R) = <0n|r_a r_b|Rm>, the second moment of the position operator, obtained from the
!>  quantum geometric tensor F_ab(k) = <d_a u_n|d_b u_m> the same way A(R) is obtained from
!>  the Berry connection.
!>
!>  What arrives here is already contracted over the pairs of neighbours and already in the
!>  Wannier gauge, because that gauge has to be applied at two different k-points at once and
!>  a slice that knows one k cannot do it. Only the Fourier sum to R is left.
MODULE m_melem_coeff_f
  USE m_juDFT
  USE m_constants, ONLY : oUnit, tpi_const
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_melem_io, ONLY : melem_write_realspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_fmn
CONTAINS

  ! F(R)_{ab,nm} = (1/N^3) sum_k e^{-ik.R} F_ab,nm(k), the Fourier partner of the geometric
  ! tensor, in the structure of Eq. (84) of Lopez, Vanderbilt, Thonhauser and Souza,
  ! PRB 85, 014435 (2012) with the Hamiltonian left out.
  !
  ! The nine components are stored with b running fastest, c = 3*(a-1) + b, which is the
  ! row-major order the cart2 format writes: xx xy xz yx yy yz zx zy zz.
  !
  ! No conversion of length: the shell weights carry Angstrom squared and the b vectors its
  ! inverse, so what arrives is already Angstrom squared.
  !
  ! Collective: every rank owns a k-slice and the sum over R is reduced across all of them.
  !
  ! NOT VALIDATED against an external reference: no other code writes this file. What it is
  ! checked against are two internal identities, the pair overlap at b1 = b2 and the
  ! hermiticity F_ab = F_ba^dagger.
  SUBROUTINE melem_write_fmn(this, kpts, f0_loc, gk_loc, irvec, nrpts, mpicm, irank, wfpref)
#ifdef CPP_MPI
    use mpi
#endif
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: f0_loc(:, :, :, :, :)   ! (nw,nw,3,3,nk_loc) gauged, per rank
    INTEGER, INTENT(IN) :: gk_loc(:)               ! (nk_loc) global k of each slice entry
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts, mpicm, irank
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref

    INTEGER :: nw, nk, nkl, kl, k, a, b, c, irpt, ierr
    REAL    :: rdotk0
    COMPLEX :: fac
    COMPLEX, ALLOCATABLE :: fr(:, :, :, :)
    CHARACTER(LEN=64) :: fn

    nw = this%num_wann; nk = kpts%nkptf; nkl = SIZE(gk_loc)

    ALLOCATE(fr(nw, nw, nrpts, 9), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO kl = 1, nkl
        k = gk_loc(kl)
        rdotk0 = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
        fac = EXP(CMPLX(0.0, -rdotk0)) / REAL(nk)
        DO a = 1, 3
          DO b = 1, 3
            c = 3*(a - 1) + b
            fr(:, :, irpt, c) = fr(:, :, irpt, c) + fac * f0_loc(:, :, a, b, kl)
          END DO
        END DO
      END DO
    END DO
#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, fr, SIZE(fr), MPI_DOUBLE_COMPLEX, MPI_SUM, mpicm, ierr)
#endif

    IF (irank == 0) THEN
      fn = 'WF1_fmn'
      IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_fmn'
      CALL melem_write_realspace(fr, irvec, [(0, irpt = 1, nrpts)], nrpts, nw, 9, 'cart2', &
                                 TRIM(fn)//'.dat', 0)
      WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (F(R)=<0n|r_a r_b|Rm>, Ang^2)'
    END IF
    DEALLOCATE(fr)
  END SUBROUTINE melem_write_fmn

END MODULE m_melem_coeff_f
