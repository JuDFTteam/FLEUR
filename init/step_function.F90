!-------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!-------------------------------------------------------------------------------
MODULE m_step_function
   !> Contains revised subroutines to construct the step function \f$\Theta(r)\f$
   !! both on a fine real space grid as well as in \f$G<G_{max}\f$ space.
   USE m_juDFT
   USE m_types

   IMPLICIT NONE

CONTAINS
   SUBROUTINE stepf_analytical(sym, stars, atoms, oneD, input, cell, fmpi, fftgrid, qvec, iDtype, iDir)
      !> Construct the analytical representation of the step function on a big
      !! reciprocal grid.
#include"cpp_double.h"

      USE m_constants
      USE m_od_cylbes

      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_oneD),  INTENT(IN) :: oneD
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_mpi),   INTENT(IN) :: fmpi

      TYPE(t_fftgrid), INTENT(OUT) :: fftgrid

      REAL,    OPTIONAL, INTENT(IN) :: qvec(3)
      INTEGER, OPTIONAL, INTENT(IN) :: iDtype, iDir

      INTEGER :: x, y, z, z2, z_min, z_max
      INTEGER :: layerDim, ifftd, gInd, n, na, nn
      REAL    :: th, inv_omtil, fJ
      REAL    :: g2, absg, absgr, help, fp_omtil, radg, gExtx, gExty
      REAL    :: r_c, r_phs
      COMPLEX :: c_c, c_phs

      REAL    :: gInt(3), gExt(3)

      CALL timestart("New stepf")
      ifftd = 27*stars%mx1*stars%mx2*stars%mx3
      layerDim = 9*stars%mx1*stars%mx2

      ! Initialize fftgrid
      CALL fftgrid%init((/3*stars%mx1,3*stars%mx2,3*stars%mx3/))

      inv_omtil = 1.0/cell%omtil
      fp_omtil = -fpi_const*inv_omtil

      fftgrid%grid = CMPLX(0.0,0.0)

      IF (.NOT.PRESENT(qvec)) THEN
         DO n = 1, atoms%ntype
            fftgrid%grid(0) = fftgrid%grid(0) + atoms%neq(n)*atoms%volmts(n)
         END DO
         fftgrid%grid(0) = 1.0 - fftgrid%grid(0)*inv_omtil
      END IF

      z_min = 0
      z_max = 3*stars%mx3 - 1

      DO z = z_min, z_max
         gInt(3) = REAL(z)
         IF (2*z > 3*stars%mx3) gInt(3) = gInt(3) - 3.0*stars%mx3
         DO y = 0, 3*stars%mx2 - 1
            gInt(2) = REAL(y)
            IF (2*y > 3*stars%mx2) gInt(2) = gInt(2) - 3.0*stars%mx2
            x_loop: DO x = 0, 3*stars%mx1 - 1
               gInt(1) = REAL(x)
               IF (2*x > 3*stars%mx1) gInt(1) = gInt(1) - 3.0*stars%mx1
               gInd = x + 3*stars%mx1*y + layerDim*z
               IF (ALL([x,y,z]==0)) CYCLE x_loop

               gExt = MATMUL(TRANSPOSE(cell%bmat), gInt)
               IF (PRESENT(qvec)) gExt = gExt + MATMUL(TRANSPOSE(cell%bmat), qvec)
               g2 = DOT_PRODUCT(gExt, gExt)
               absg = SQRT(g2)
               help = fp_omtil/g2

               c_c = CMPLX(0.0, 0.0)
               DO n = 1, MERGE(atoms%ntype,iDtype,.NOT.PRESENT(qvec))
                  c_phs = CMPLX(0.0, 0.0)
                  na = MERGE(SUM(atoms%neq(:n - 1)),iDtype,.NOT.PRESENT(qvec))
                  DO nn = 1, atoms%neq(n)
                     th = -tpi_const*DOT_PRODUCT(gInt, atoms%taual(:, na + nn))
                     IF (.NOT.PRESENT(qvec)) THEN
                        c_phs = c_phs + EXP(CMPLX(0, th))
                     ELSE
                        c_phs = c_phs + (-ImagUnit*gExt(iDir))*EXP(CMPLX(0, th))
                     END IF
                  END DO
                  absgr = absg*atoms%rmt(n)
                  c_c = c_c + atoms%rmt(n)*(SIN(absgr)/absgr - COS(absgr))*c_phs
               END DO
               fftgrid%grid(gInd) = help * c_c

               !IF (((z  .EQ. 3*stars%mx3/2) .OR. (y  .EQ. 3*stars%mx2/2)) .OR. (x  .EQ. 3*stars%mx1/2)) THEN
               IF (((2*z  .EQ. 3*stars%mx3) .OR. (2*y  .EQ. 3*stars%mx2)) .OR. (2*x  .EQ. 3*stars%mx1)) THEN
                  fftgrid%grid(gInd) = CMPLX(0.0,0.0)
               END IF

               IF (oneD%odd%d1) THEN ! fmpi version untested, but obsolete anyway
                  IF (gInd<layerDim.AND.gInd.NE.0) THEN
                     gExtx = (cell%bmat(1, 1)*x + cell%bmat(2, 1)*y)
                     gExty = (cell%bmat(1, 2)*x + cell%bmat(2, 2)*y)
                     radg  = SQRT(gExtx**2 + gExty**2)
                     CALL od_cylbes(1, radg*cell%z1, fJ)
                     fftgrid%grid(gInd) = fftgrid%grid(gInd) + 2*cell%vol*fJ/(radg*cell%z1*cell%omtil)
                  END IF
               END IF
            END DO x_loop ! 0, 3*stars%mx1 - 1
         END DO ! 0, 3*stars%mx2 - 1
      END DO ! 0, 3*stars%mx3 - 1

      IF (fmpi%irank == 0) THEN
         IF (input%film) THEN
            fftgrid%grid(0) = fftgrid%grid(0) + cell%vol*inv_omtil - 1.0
            IF (.NOT.oneD%odd%d1) THEN
               DO z2 = 1, 3*stars%mx3 - 1
                  gInt(3) = REAL(z2)
                  IF (gInt(3) > 1.5*stars%mx3) gInt(3) = gInt(3) - 3.0 * stars%mx3
                  th = cell%bmat(3, 3)*gInt(3)*cell%z1
                  fftgrid%grid(z2*layerDim) = fftgrid%grid(z2*layerDim) + cell%vol*inv_omtil*SIN(th)/th
               END DO
            END IF
         END IF
      END IF

      CALL timestop("New stepf")
   END SUBROUTINE

   SUBROUTINE stepf_stars(stars,fftgrid)

      TYPE(t_stars),   INTENT(INOUT) :: stars
      TYPE(t_fftgrid), INTENT(INOUT) :: fftgrid

!#ifdef CPP_MPI
!      INTEGER :: ierr
!      CALL MPI_BCAST(stars%ng3, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
!      CALL MPI_BCAST(stars%ig, stars%ng3, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
!#endif

      ! Save G coefficients to ustep
      CALL fftgrid%takeFieldFromGrid(stars, stars%ustep)
      stars%ustep = stars%ustep * 3 * stars%mx1 * 3 * stars%mx2 * 3 * stars%mx3
      CALL fftgrid%perform_fft(forward=.FALSE.)

      stars%ufft=fftgrid%grid

   END SUBROUTINE
END MODULE m_step_function
