!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  C(R) = <0n|(r_a - R_a) H (r_b - R_b)|Rm>, obtained from C_ab(k) = <d_a u_n|H|d_b u_m>
!>  the same way F(R) is obtained from the geometric tensor. Eq. (84) of Lopez, Vanderbilt,
!>  Thonhauser and Souza, PRB 85, 014435 (2012), with the Hamiltonian left in.
!>
!>  What arrives here is already contracted over the pairs of neighbours and already in the
!>  Wannier gauge, for the same reason F is: the gauge has to be applied at two different
!>  k-points at once and a slice that knows one k cannot do it. Only the Fourier sum to R
!>  is left.
MODULE m_melem_coeff_c
  USE m_juDFT
  USE m_constants, ONLY : oUnit, tpi_const
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_melem_io, ONLY : melem_write_realspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_cmn
CONTAINS

  ! C(R)_{ab,nm} = (1/N^3) sum_k e^{-ik.R} C_ab,nm(k).
  !
  ! Nine components with b running fastest, c = 3*(a-1) + b: the row-major order the cart2
  ! family writes, xx xy xz yx yy yz zx zy zz. Same as F, so the two files are read the
  ! same way.
  !
  ! Units, and this is where C parts company with F. The shell weights carry Angstrom
  ! squared, as they do there, but the Hamiltonian in the middle carries energy, so what
  ! arrives is Hartree times Angstrom squared. It goes out through 'cart2e', which converts
  ! the energy, so that C sits in eV*Ang^2 next to an H(R) and a B(R) that are in eV.
  !
  ! Collective: every rank owns a k-slice and the sum over R is reduced across all of them.
  !
  ! NOT VALIDATED against an external reference: no other code writes this file. What
  ! backs the C(k) reaching here is the identity it is assembled from; the layout itself,
  ! and the convention its consumer reads it with, are unverified exactly as F's are.
  SUBROUTINE melem_write_cmn(this, kpts, c0_loc, gk_loc, irvec, nrpts, mpicm, irank, wfpref)
#ifdef CPP_MPI
    use mpi
#endif
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: c0_loc(:, :, :, :, :)   ! (nw,nw,3,3,nk_loc) gauged, per rank
    INTEGER, INTENT(IN) :: gk_loc(:)               ! (nk_loc) global k of each slice entry
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts, mpicm, irank
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref

    INTEGER :: nw, nk, nkl, kl, k, a, b, c, irpt, ierr
    REAL    :: rdotk0
    COMPLEX :: fac
    COMPLEX, ALLOCATABLE :: cr(:, :, :, :)
    CHARACTER(LEN=64) :: fn

    nw = this%num_wann; nk = kpts%nkptf; nkl = SIZE(gk_loc)

    ALLOCATE(cr(nw, nw, nrpts, 9), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO kl = 1, nkl
        k = gk_loc(kl)
        rdotk0 = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
        fac = EXP(CMPLX(0.0, -rdotk0)) / REAL(nk)
        DO a = 1, 3
          DO b = 1, 3
            c = 3*(a - 1) + b
            cr(:, :, irpt, c) = cr(:, :, irpt, c) + fac * c0_loc(:, :, a, b, kl)
          END DO
        END DO
      END DO
    END DO
#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, cr, SIZE(cr), MPI_DOUBLE_COMPLEX, MPI_SUM, mpicm, ierr)
#endif

    IF (irank == 0) THEN
      fn = 'WF1_cmn'
      IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_cmn'
      CALL melem_write_realspace(cr, irvec, [(0, irpt = 1, nrpts)], nrpts, nw, 9, 'cart2e', &
                                 TRIM(fn)//'.dat', 0)
      WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (C(R)=<0n|r_a H r_b|Rm>, eV*Ang^2)'
    END IF
    DEALLOCATE(cr)
  END SUBROUTINE melem_write_cmn

END MODULE m_melem_coeff_c
