!-------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!-------------------------------------------------------------------------------
MODULE m_step_function
   !! Contains revised subroutines to construct the step function \f$\Theta(r)\f$
   !! both on a fine real space grid as well as in \f$G<G_{max}\f$ space.
   USE m_juDFT
   USE m_types

   IMPLICIT NONE

CONTAINS
   SUBROUTINE stepf_analytical(sym, stars, atoms, input, cell, fmpi, fftgrid, qvec, iDtype, iDir, iOrd, stepf_array)
      !! Construct the analytical representation of the step function on a big
      !! reciprocal grid.

      USE m_constants

      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_cell),  INTENT(IN) :: cell
      TYPE(t_mpi),   INTENT(IN) :: fmpi

      TYPE(t_fftgrid), INTENT(OUT) :: fftgrid

      REAL,    OPTIONAL, INTENT(IN) :: qvec(3)
      INTEGER, OPTIONAL, INTENT(IN) :: iDtype, iDir, iOrd

      COMPLEX, OPTIONAL, INTENT(INOUT) :: stepf_array(0:,:,:)

      INTEGER :: x, y, z, z2, z_min, z_max
      INTEGER :: layerDim, ifftd, gInd, n, na, nn, iDir_loc
      INTEGER :: n_lims(2), dir_lims(2)
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

      n_lims = [1, atoms%ntype]
      dir_lims = [1, 1]

      IF (PRESENT(qvec)) THEN
         n_lims(1) = MERGE(1,          iDtype,PRESENT(stepf_array).AND.iOrd==1)
         n_lims(2) = MERGE(atoms%ntype,iDtype,PRESENT(stepf_array).AND.iOrd==1)

         dir_lims(1) = MERGE(1,iDir,PRESENT(stepf_array))
         dir_lims(2) = MERGE(3,iDir,PRESENT(stepf_array))
      END IF

      DO z = z_min, z_max
         gInt(3) = REAL(z)
         IF (2*z > 3*stars%mx3) gInt(3) = gInt(3) - 3.0*stars%mx3
         DO y = 0, 3*stars%mx2 - 1
            gInt(2) = REAL(y)
            IF (2*y > 3*stars%mx2) gInt(2) = gInt(2) - 3.0*stars%mx2
            x_loop: DO x = 0, 3*stars%mx1 - 1
               ! TODO: Maybe add the possibility to skip calculations for elements
               !       that are already available by inversion; problematic when
               !       parallelized. Also: parallelize!
               gInt(1) = REAL(x)
               IF (2*x > 3*stars%mx1) gInt(1) = gInt(1) - 3.0*stars%mx1
               gInd = x + 3*stars%mx1*y + layerDim*z

               IF (PRESENT(qvec)) THEN
                  IF (ALL([x,y,z]==0).AND.norm2(qvec)<1e-9) CYCLE x_loop
               ELSE
                  IF (ALL([x,y,z]==0)) CYCLE x_loop
               END IF

               gExt = MATMUL(TRANSPOSE(cell%bmat), gInt)
               IF (PRESENT(qvec)) gExt = gExt + MATMUL(TRANSPOSE(cell%bmat), qvec)
               g2 = DOT_PRODUCT(gExt, gExt)
               absg = SQRT(g2)
               help = fp_omtil/g2

               c_c = CMPLX(0.0, 0.0)

               DO n = n_lims(1), n_lims(2)
                  DO iDir_loc = dir_lims(1), dir_lims(2)
                     c_phs = CMPLX(0.0, 0.0)
                     !na = MERGE(SUM(atoms%neq(:n - 1)),iDtype-1,.NOT.PRESENT(qvec))
                     na = SUM(atoms%neq(:n - 1))
                     DO nn = 1, atoms%neq(n)
                        IF (.NOT.PRESENT(qvec)) THEN
                           th = -tpi_const*DOT_PRODUCT(gInt, atoms%taual(:, na + nn))
                           c_phs = c_phs + EXP(CMPLX(0, th))
                        ELSE
                           th = -tpi_const*DOT_PRODUCT(gInt+qvec, atoms%taual(:, na + nn))
                           IF (iOrd==1) THEN
                              c_phs = c_phs + (-ImagUnit*gExt(iDir_loc))*EXP(CMPLX(0, th))
                           ELSE
                              c_phs = c_phs + (-gExt(iDir)*gExt(iDir_loc))*EXP(CMPLX(0, th))
                           END IF
                        END IF
                     END DO
                     absgr = absg*atoms%rmt(n)
                     c_c = c_c + atoms%rmt(n)*(SIN(absgr)/absgr - COS(absgr))*c_phs

                     IF (PRESENT(stepf_array)) THEN
                        stepf_array(gInd, n, iDir_loc) = help * c_c
                        c_c = CMPLX(0.0, 0.0)
                     END IF
                  END DO
               END DO

               IF (.NOT.PRESENT(stepf_array)) THEN
                  fftgrid%grid(gInd) = help * c_c

                  IF (((2*z  .EQ. 3*stars%mx3) .OR. (2*y  .EQ. 3*stars%mx2)) .OR. (2*x  .EQ. 3*stars%mx1)) THEN
                     fftgrid%grid(gInd) = CMPLX(0.0,0.0)
                  END IF
               ELSE
                  IF (((2*z  .EQ. 3*stars%mx3) .OR. (2*y  .EQ. 3*stars%mx2)) .OR. (2*x  .EQ. 3*stars%mx1)) THEN
                     stepf_array(gInd, :, :) = CMPLX(0.0,0.0)
                  END IF
               END IF

            END DO x_loop ! 0, 3*stars%mx1 - 1
         END DO ! 0, 3*stars%mx2 - 1
      END DO ! 0, 3*stars%mx3 - 1

      IF (fmpi%irank == 0) THEN
         IF (input%film) THEN
            fftgrid%grid(0) = fftgrid%grid(0) + cell%vol*inv_omtil - 1.0
            DO z2 = 1, 3*stars%mx3 - 1
               gInt(3) = REAL(z2)
               IF (gInt(3) > 1.5*stars%mx3) gInt(3) = gInt(3) - 3.0 * stars%mx3
               th = cell%bmat(3, 3)*gInt(3)*cell%z1
               fftgrid%grid(z2*layerDim) = fftgrid%grid(z2*layerDim) + cell%vol*inv_omtil*SIN(th)/th
            END DO
         END IF
      END IF

      CALL timestop("New stepf")
   END SUBROUTINE

   SUBROUTINE stepf_stars(stars,fftgrid,qvec)

      TYPE(t_stars),   INTENT(INOUT) :: stars
      TYPE(t_fftgrid), INTENT(INOUT) :: fftgrid
      REAL,    OPTIONAL, INTENT(IN) :: qvec(3)

!#ifdef CPP_MPI
!      INTEGER :: ierr
!      CALL MPI_BCAST(stars%ng3, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
!      CALL MPI_BCAST(stars%ig, stars%ng3, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
!#endif

      ! Save G coefficients to ustep
      CALL fftgrid%takeFieldFromGrid(stars, stars%ustep)
      stars%ustep = stars%ustep * 3 * stars%mx1 * 3 * stars%mx2 * 3 * stars%mx3
      CALL fftgrid%perform_fft(forward=.FALSE.)

      IF (.NOT.PRESENT(qvec)) THEN
         stars%ufft=fftgrid%grid
      ELSE
         stars%ufft1=fftgrid%grid
      END IF

   END SUBROUTINE
END MODULE m_step_function
