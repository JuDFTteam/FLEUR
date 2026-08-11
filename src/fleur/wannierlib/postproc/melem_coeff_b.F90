!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grunberg Institut, Forschungszentrum Julich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  B(R) = <0n|H (r-R)|Rm>, Eq. (83) of Lopez, Vanderbilt, Thonhauser and Souza,
!>  PRB 85, 014435 (2012).
!>
!>  The Hamiltonian enters only through the ab-initio eigenvalues, so B comes out of the
!>  same neighbour overlaps the wannierisation was given and needs no uHu.
MODULE m_melem_coeff_b
  USE m_juDFT
  USE m_constants, ONLY : oUnit, tpi_const
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_types_melem_bmesh
  USE m_melem_io, ONLY : melem_write_realspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_bmn
CONTAINS

  ! ---------------------------------------------------------------------------
  ! B(R) = <0n| H r_alpha |Rm>, Eq. (83) of Lopez, Vanderbilt, Thonhauser and Souza,
  ! PRB 85, 014435 (2012), in the translationally invariant form (the r_avg centring).
  !
  ! It needs NO uHu. The Hamiltonian enters only through the ab-initio eigenvalues:
  !
  !     M^H(k,b) = V(k)^dagger diag(eig(k)) M^(0,b)(k) V(k+b)
  !
  ! which is the same neighbour overlap the wannierization was given, weighted by the
  ! eigenvalues before the gauge is applied on the left.
  !
  ! The (r - R) of the definition is already carried by the k derivative of the periodic
  ! part, so the sum takes no centring term and needs no Wannier centres.
  !
  ! NOT VALIDATED NUMERICALLY against a reference: the test suite has none for B.
  SUBROUTINE melem_write_bmn(this, kpts, eig, u_matrix, u_opt, mmn_loc, gk_loc, &
                             bmesh, irvec, nrpts, mpicm, irank, wfpref)
#ifdef CPP_MPI
    use mpi
#endif
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)                 ! (nb,nk) ab-initio, Hartree
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)         ! (nw,nw,nk)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)            ! (nb,nw,nk)
    COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)       ! (nb,nb,nntot,nk_loc) this rank's slice
    INTEGER, INTENT(IN) :: gk_loc(:)                 ! (nk_loc) global k of each slice entry
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts, mpicm, irank
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref

    INTEGER :: nb, nw, nk, nnt, nkl, kl, k, kb, nn, a, i, j, m, n, irpt, ierr
    REAL    :: rdotk0
    COMPLEX :: fac
    COMPLEX, ALLOCATABLE :: Vk(:, :), Vkb(:, :), tmp(:, :), mh(:, :, :, :), br(:, :, :, :)
    CHARACTER(LEN=64) :: fn

    nb = this%num_bands; nw = this%num_wann; nk = kpts%nkptf
    nnt = bmesh%nntot; nkl = SIZE(gk_loc)
    IF (nnt < 1 .OR. .NOT. ALLOCATED(bmesh%wb)) THEN
      IF (irank == 0) WRITE(oUnit,'(a)') &
        'wannierlib operators_r: B(R) needs the b-shell weights -> skipped'
      RETURN
    END IF

    ! ---- M^H(k,b) = V(k)^dag diag(eig) M V(k+b) on this rank's k-slice ----
    ALLOCATE(mh(nw, nw, nnt, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(Vk(nb, nw), Vkb(nb, nw), tmp(nb, nw))
    DO kl = 1, nkl
      k = gk_loc(kl)
      Vk = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
      DO i = 1, nb                      ! the eigenvalue weighting acts on the k side
        Vk(i, :) = eig(i, k) * Vk(i, :)
      END DO
      DO nn = 1, nnt
        kb = bmesh%nnlist(k, nn)
        Vkb = MATMUL(u_opt(:, :, kb), u_matrix(:, :, kb))
        tmp = MATMUL(mmn_loc(:, :, nn, kl), Vkb)
        mh(:, :, nn, kl) = MATMUL(CONJG(TRANSPOSE(Vk)), tmp)
      END DO
    END DO
    DEALLOCATE(Vk, Vkb, tmp)

    ! ---- the R sum ----
    ALLOCATE(br(3, nw, nw, nrpts), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO m = 1, nw
        DO n = 1, nw
          DO kl = 1, nkl
            k = gk_loc(kl)
            rdotk0 = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
            fac = EXP(CMPLX(0.0, -rdotk0)) / REAL(nk)
            DO nn = 1, nnt
              DO a = 1, 3
                br(a, n, m, irpt) = br(a, n, m, irpt) + &
                  CMPLX(0.0, bmesh%wb(nn) * bmesh%bk(a, nn, k)) * mh(n, m, nn, kl) * fac
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(mh)
#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, br, SIZE(br), MPI_DOUBLE_COMPLEX, MPI_SUM, mpicm, ierr)
#endif

    IF (irank == 0) THEN
      fn = 'WF1_bmn'
      IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_bmn'
      BLOCK
        COMPLEX, ALLOCATABLE :: b4(:, :, :, :)
        ALLOCATE(b4(nw, nw, nrpts, 3))
        DO a = 1, 3
          DO irpt = 1, nrpts
            b4(:, :, irpt, a) = br(a, :, :, irpt)
          END DO
        END DO
        CALL melem_write_realspace(b4, irvec, [(0, i = 1, nrpts)], nrpts, nw, 3, 'bmn', &
                                   TRIM(fn)//'.dat', 0)
        DEALLOCATE(b4)
      END BLOCK
      WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (B(R)=<0n|H (r-R)|Rm>, eV*Ang) '// &
                         '-- NOT numerically validated'
    END IF
    DEALLOCATE(br)
  END SUBROUTINE melem_write_bmn

END MODULE m_melem_coeff_b
