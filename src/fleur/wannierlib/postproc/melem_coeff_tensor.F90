!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The two real-space tensors with two Cartesian indices, in the structure of Eq. (84) of
!>  Lopez, Vanderbilt, Thonhauser and Souza, PRB 85, 014435 (2012):
!>
!>     F(R) = <0n|r_a r_b|Rm>        the second moment of the position operator,
!>                                   the Fourier partner of the quantum geometric
!>                                   tensor F_ab(k) = <d_a u_n|d_b u_m>
!>     C(R) = <0n|r_a H r_b|Rm>      the same one with the Hamiltonian in the middle
!>
!>  They are two different physical objects and keep two entry points, but the Fourier sum
!>  is the same and used to be written twice. What separates them is one line: the length
!>  units are already Angstrom squared in both, and the Hamiltonian in C adds an energy that
!>  the 'cart2e' layout converts to eV, so that C sits next to an H(R) and a B(R) that are
!>  in eV. F needs no conversion and goes out through 'cart2'.
!>
!>  What arrives here is already contracted over the pairs of neighbours and already in the
!>  Wannier gauge, because that gauge has to be applied at two different k-points at once and
!>  a slice that knows one k cannot do it. Only the Fourier sum to R is left.
!>
!>  NEITHER IS VALIDATED against an external reference: no other code writes these files.
!>  What backs F are two internal identities, the pair overlap at b1 = b2 and the hermiticity
!>  F_ab = F_ba^dagger; what backs the C(k) reaching here is the identity it is assembled
!>  from. The layout itself, and the convention its consumer reads it with, are unverified.
MODULE m_melem_coeff_tensor
  USE m_juDFT
  USE m_constants, ONLY : oUnit, tpi_const
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_melem_io, ONLY : melem_write_realspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_fmn, melem_write_cmn
CONTAINS

  !> F(R), in Angstrom squared. No conversion: the shell weights carry Angstrom squared and
  !> the b vectors its inverse, so what arrives is already in those units.
  SUBROUTINE melem_write_fmn(this, kpts, f0_loc, gk_loc, irvec, nrpts, mpicm, irank, wfpref)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: f0_loc(:, :, :, :, :)   ! (nw,nw,3,3,nk_loc) gauged, per rank
    INTEGER, INTENT(IN) :: gk_loc(:)               ! (nk_loc) global k of each slice entry
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts, mpicm, irank
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref

    CALL write_tensor_r(this, kpts, f0_loc, gk_loc, irvec, nrpts, mpicm, irank, &
                        'cart2', '_fmn', 'F(R)=<0n|r_a r_b|Rm>, Ang^2', wfpref)
  END SUBROUTINE melem_write_fmn

  !> C(R), in eV*Angstrom squared. This is where C parts company with F: the Hamiltonian in
  !> the middle carries energy, so what arrives is Hartree times Angstrom squared and the
  !> 'cart2e' layout converts it.
  SUBROUTINE melem_write_cmn(this, kpts, c0_loc, gk_loc, irvec, nrpts, mpicm, irank, wfpref)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: c0_loc(:, :, :, :, :)   ! (nw,nw,3,3,nk_loc) gauged, per rank
    INTEGER, INTENT(IN) :: gk_loc(:)               ! (nk_loc) global k of each slice entry
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts, mpicm, irank
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref

    CALL write_tensor_r(this, kpts, c0_loc, gk_loc, irvec, nrpts, mpicm, irank, &
                        'cart2e', '_cmn', 'C(R)=<0n|r_a H r_b|Rm>, eV*Ang^2', wfpref)
  END SUBROUTINE melem_write_cmn

  !> O(R)_{ab,nm} = (1/N^3) sum_k e^{-ik.R} O_ab,nm(k), for either of the two.
  !>
  !> Nine components with b running fastest, c = 3*(a-1) + b: the row-major order the cart2
  !> family writes, xx xy xz yx yy yz zx zy zz. Both use it, so the two files are read the
  !> same way.
  !>
  !> The pairing of layout and suffix stays here rather than at the call site: the caller
  !> asks for C or for F, not for a format. Getting 'cart2e' onto a file named _fmn would
  !> convert an energy that is not there.
  !>
  !> Collective: every rank owns a k-slice and the sum over R is reduced across all of them.
  SUBROUTINE write_tensor_r(this, kpts, o0_loc, gk_loc, irvec, nrpts, mpicm, irank, &
                            fmt, suffix, what, wfpref)
#ifdef CPP_MPI
    use mpi
#endif
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: o0_loc(:, :, :, :, :)
    INTEGER, INTENT(IN) :: gk_loc(:)
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts, mpicm, irank
    CHARACTER(LEN=*), INTENT(IN) :: fmt      !> 'cart2' or 'cart2e'
    CHARACTER(LEN=*), INTENT(IN) :: suffix   !> '_fmn' or '_cmn'
    CHARACTER(LEN=*), INTENT(IN) :: what     !> for the line written to the output file
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref

    !> The cart2 family writes no degeneracy block, so melem_write_realspace ignores this
    !> argument for these two layouts. Named rather than spelled out at the call site, so
    !> that a reader does not have to check the writer to see that the zeros mean nothing.
    INTEGER, ALLOCATABLE :: ndegen_unused(:)

    INTEGER :: nw, nk, nkl, kl, k, a, b, c, irpt, ierr
    REAL    :: rdotk0
    COMPLEX :: fac
    COMPLEX, ALLOCATABLE :: o_r(:, :, :, :)
    CHARACTER(LEN=64) :: fn

    nw = this%num_wann; nk = kpts%nkptf; nkl = SIZE(gk_loc)

    ALLOCATE(o_r(nw, nw, nrpts, 9), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO kl = 1, nkl
        k = gk_loc(kl)
        rdotk0 = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
        fac = EXP(CMPLX(0.0, -rdotk0)) / REAL(nk)
        DO a = 1, 3
          DO b = 1, 3
            c = 3*(a - 1) + b
            o_r(:, :, irpt, c) = o_r(:, :, irpt, c) + fac * o0_loc(:, :, a, b, kl)
          END DO
        END DO
      END DO
    END DO
#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, o_r, SIZE(o_r), MPI_DOUBLE_COMPLEX, MPI_SUM, mpicm, ierr)
#endif

    IF (irank == 0) THEN
      ALLOCATE(ndegen_unused(nrpts), source=0)
      fn = 'WF1'//suffix
      IF (PRESENT(wfpref)) fn = TRIM(wfpref)//suffix
      CALL melem_write_realspace(o_r, irvec, ndegen_unused, nrpts, nw, 9, fmt, &
                                 TRIM(fn)//'.dat', 0)
      WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat ('//what//')'
      DEALLOCATE(ndegen_unused)
    END IF
    DEALLOCATE(o_r)
  END SUBROUTINE write_tensor_r

END MODULE m_melem_coeff_tensor
