MODULE m_stepf
   USE m_juDFT
   USE m_cdn_io
#ifdef CPP_MPI
   use mpi 
#endif
CONTAINS
   SUBROUTINE stepf(sym, stars, atoms, oneD, input, cell, vacuum, fmpi)
      !
      !*********************************************************************
      !     calculates the fourier components of the interstitial step
      !     function for the reciprocal vectors of the star list.
      !           m. weinert  1986
      !*********************************************************************
      !
      !     also set up FFT of U(G) on a (-2G:+2G) grid for convolutions
      !
      !*********************************************************************
#include"cpp_double.h"
      USE m_cfft
      USE m_constants
      USE m_od_cylbes
      USE m_types
      IMPLICIT NONE
      !     ..
      TYPE(t_sym), INTENT(IN)        :: sym
      TYPE(t_stars), INTENT(INOUT)   :: stars
      TYPE(t_atoms), INTENT(IN)      :: atoms
      TYPE(t_oneD), INTENT(IN)       :: oneD
      TYPE(t_input), INTENT(IN)      :: input
      TYPE(t_cell), INTENT(IN)       :: cell
      TYPE(t_vacuum), INTENT(IN)     :: vacuum
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      !     ..
      !     .. Local Scalars ..
      COMPLEX c_c, c_phs
      REAL c, dd, gs, th, inv_omtil, r_phs
      REAL g_rmt, g_sqr, help, g_abs, fp_omtil, r_c, gr, gx, gy
      INTEGER i, k, n, na, nn, i1, i2, i3, ic, ifft2d, ifftd, kk
      INTEGER ic1, ic2, ic3, icc, im1, im2, im3, loopstart
      LOGICAL l_error
      !     ..
      !     .. Local Arrays ..
      COMPLEX, ALLOCATABLE:: sf(:)
      REAL g(3), gm(3), fJ
      REAL, ALLOCATABLE :: bfft(:)
      INTEGER, ALLOCATABLE :: icm(:, :, :)
      INTEGER :: i3_start, i3_end, chunk_size, leftover_size
#ifdef CPP_MPI
      INTEGER ierr
      INTEGER, ALLOCATABLE :: icm_local(:, :, :)
      REAL, ALLOCATABLE :: ufft_local(:), bfft_local(:)

      CALL MPI_BCAST(stars%mx1, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
      CALL MPI_BCAST(stars%mx2, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
      CALL MPI_BCAST(stars%mx3, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
#endif

      ifftd = 27*stars%mx1*stars%mx2*stars%mx3
      !     ..
      !     ..
      !--->    if step function stored on disc, then just read it in
      !
      l_error = .FALSE.
      IF (fmpi%irank == 0) CALL readStepfunction(stars, atoms, cell, vacuum, l_error)
#ifdef CPP_MPI
      CALL MPI_BCAST(l_error, 1, MPI_LOGICAL, 0, fmpi%mpi_comm, ierr)
#endif
      IF (.NOT. l_error) THEN
         RETURN
      END IF

      IF (fmpi%irank == 0) THEN
         ALLOCATE (sf(stars%ng3))
         IF (input%film) THEN
            dd = vacuum%dvac*cell%area/cell%omtil
            IF (oneD%odd%d1) dd = cell%vol/cell%omtil
         ELSE
            dd = 1.0
         END IF
         !--->    G=0 star
         c = 0.0
         DO n = 1, atoms%ntype
            c = c + atoms%neq(n)*atoms%volmts(n)/cell%omtil
         ENDDO
         stars%ustep(1) = CMPLX(dd - c, 0.0)
         !--->    G(parallel)=0  (for film)
         IF (input%film .AND. .NOT. oneD%odd%d1) THEN
            DO k = 2, stars%ng3
               IF (stars%ig2(k) .EQ. 1) THEN
                  th = cell%bmat(3, 3)*stars%kv3(3, k)*cell%z1
                  stars%ustep(k) = CMPLX(cell%vol*SIN(th)/th/cell%omtil, 0.0)
               ELSE
                  stars%ustep(k) = CMPLX(0.0, 0.0)
               END IF
            ENDDO
            !-odim
         ELSEIF (oneD%odd%d1) THEN
            DO k = 2, stars%ng3
               gr = 0.0
               IF (stars%kv3(3, k) .EQ. 0) THEN
                  kk = stars%ig2(k)
                  gr = stars%sk2(kk)
                  CALL od_cylbes(1, gr*cell%z1, fJ)
                  stars%ustep(k) = CMPLX(2.*dd*fJ/(gr*cell%z1), 0.)
               ELSE
                  stars%ustep(k) = CMPLX(0., 0.)
               END IF

            ENDDO
            !+odim
         ELSE
            DO k = 2, stars%ng3
               stars%ustep(k) = CMPLX(0.0, 0.0)
            END DO
         END IF
         !--->    sphere contributions
         na = 0
         DO n = 1, atoms%ntype
            c = 3.*atoms%volmts(n)/cell%omtil
            !-->     structure factors: loop over equivalent atoms
            na = na + 1
            DO k = 2, stars%ng3
               th = -tpi_const*DOT_PRODUCT(stars%kv3(:, k), atoms%taual(:, na))
               sf(k) = CMPLX(COS(th), SIN(th))
            END DO
            DO nn = 2, atoms%neq(n)
               na = na + 1
               DO k = 2, stars%ng3
                  th = -tpi_const*DOT_PRODUCT(stars%kv3(:, k), atoms%taual(:, na))
                  sf(k) = sf(k) + CMPLX(COS(th), SIN(th))
               END DO
            END DO
            !--->    update step function
            DO k = 2, stars%ng3
               gs = stars%sk3(k)*atoms%rmt(n)
               stars%ustep(k) = stars%ustep(k) - (c*(SIN(gs)/gs - COS(gs))/(gs*gs))*sf(k)
            ENDDO
         ENDDO
      ENDIF ! (fmpi%irank == 0)

      !
      ! --> set up stepfunction on fft-grid:
      !
#ifdef CPP_MPI
      ALLOCATE (bfft_local(0:ifftd - 1), ufft_local(0:ifftd - 1))
      bfft_local = 0.0
      ufft_local = 0.0
#endif

      ALLOCATE (bfft(0:ifftd - 1))
      im1 = CEILING(1.5*stars%mx1); im2 = CEILING(1.5*stars%mx2); im3 = CEILING(1.5*stars%mx3)
      ALLOCATE (icm(-im1:im1, -im2:im2, -im3:im3))
      icm = 0
#ifdef CPP_MPI
      ALLOCATE (icm_local(-im1:im1, -im2:im2, -im3:im3))
      icm_local = 0
#endif

      inv_omtil = 1.0/cell%omtil
      fp_omtil = -fpi_const*inv_omtil
      !DO first vector before loop
      stars%ufft = 0.0
      bfft = 0.0

#ifdef CPP_MPI
      IF (fmpi%irank == 0) THEN
         DO n = 1, atoms%ntype
            ufft_local(0) = ufft_local(0) + atoms%neq(n)*atoms%volmts(n)
         ENDDO
         ufft_local(0) = 1.0 - ufft_local(0)*inv_omtil
      ENDIF
#else
      DO n = 1, atoms%ntype
         stars%ufft(0) = stars%ufft(0) + atoms%neq(n)*atoms%volmts(n)
      ENDDO
      stars%ufft(0) = 1.0 - stars%ufft(0)*inv_omtil
#endif

      i3_start = 0
      i3_end = im3
#ifdef CPP_MPI
      chunk_size = (im3 + 1)/fmpi%isize
      leftover_size = modulo(im3 + 1, fmpi%isize)
      IF (fmpi%irank < leftover_size) THEN
         i3_start = fmpi%irank*(chunk_size + 1)
         i3_end = (fmpi%irank + 1)*(chunk_size + 1) - 1
      ELSE
         i3_start = leftover_size*(chunk_size + 1) + (fmpi%irank - leftover_size)*chunk_size
         i3_end = (i3_start + chunk_size) - 1
      ENDIF
#endif
      DO i3 = i3_start, i3_end
         gm(3) = REAL(i3)
         DO i2 = 0, 3*stars%mx2 - 1
            gm(2) = REAL(i2)
            IF (2*i2 > 3*stars%mx2) gm(2) = gm(2) - 3.0*stars%mx2
            IF (i3 == 0 .AND. i2 == 0) THEN
               loopstart = 1
            ELSE
               loopstart = 0
            ENDIF
            DO i1 = loopstart, 3*stars%mx1 - 1
               ic = i1 + 3*stars%mx1*i2 + 9*stars%mx1*stars%mx2*i3
               gm(1) = REAL(i1)
               IF (2*i1 > 3*stars%mx1) gm(1) = gm(1) - 3.0*stars%mx1
               !
               !-> use inversion <-> c.c.
               !
               ic1 = NINT(gm(1)); ic2 = NINT(gm(2)); ic3 = NINT(gm(3))
#ifdef CPP_MPI
               icm_local(ic1, ic2, ic3) = ic
               IF (ic1 == im1) icm_local(-ic1, ic2, ic3) = ic
               IF (ic2 == im2) icm_local(ic1, -ic2, ic3) = ic
               IF ((ic1 == im1) .AND. (ic2 == im2)) icm_local(-ic1, -ic2, ic3) = ic
#else
               icm(ic1, ic2, ic3) = ic
               IF (ic1 == im1) icm(-ic1, ic2, ic3) = ic
               IF (ic2 == im2) icm(ic1, -ic2, ic3) = ic
               IF ((ic1 == im1) .AND. (ic2 == im2)) icm(-ic1, -ic2, ic3) = ic
#endif
               g = MATMUL(TRANSPOSE(cell%bmat), gm)
               g_sqr = DOT_PRODUCT(g, g)
               g_abs = SQRT(g_sqr)
               help = fp_omtil/g_sqr
               IF (sym%invs) THEN
                  r_c = 0.0
                  DO n = 1, atoms%ntype
                     r_phs = 0.0
                     na = SUM(atoms%neq(:n - 1))
                     DO nn = 1, atoms%neq(n)
                        th = -tpi_const*DOT_PRODUCT(gm, atoms%taual(:, na + nn))
                        r_phs = r_phs + COS(th)
                     ENDDO
                     g_rmt = g_abs*atoms%rmt(n)
                     r_c = r_c + atoms%rmt(n)*(SIN(g_rmt)/g_rmt - COS(g_rmt))*r_phs
                  ENDDO
#ifdef CPP_MPI
                  ufft_local(ic) = help*r_c
#else
                  stars%ufft(ic) = help*r_c
#endif
               ELSE
                  c_c = CMPLX(0.0, 0.0)
                  DO n = 1, atoms%ntype
                     c_phs = CMPLX(0.0, 0.0)
                     na = SUM(atoms%neq(:n - 1))
                     DO nn = 1, atoms%neq(n)
                        th = -tpi_const*DOT_PRODUCT(gm, atoms%taual(:, na + nn))
                        c_phs = c_phs + EXP(CMPLX(0, th))
                     ENDDO
                     g_rmt = g_abs*atoms%rmt(n)
                     c_c = c_c + atoms%rmt(n)*(SIN(g_rmt)/g_rmt - COS(g_rmt))*c_phs
                  ENDDO
#ifdef CPP_MPI
                  ufft_local(ic) = help*REAL(c_c)
                  bfft_local(ic) = help*AIMAG(c_c)
#else
                  stars%ufft(ic) = help*REAL(c_c)
                  bfft(ic) = help*AIMAG(c_c)
#endif
               ENDIF

               IF (((i3 .EQ. 3*stars%mx3/2) .OR. (i2 .EQ. 3*stars%mx2/2)) .OR. (i1 .EQ. 3*stars%mx1/2)) THEN
#ifdef CPP_MPI
                  ufft_local(ic) = 0.0
                  bfft_local(ic) = 0.0
#else
                  stars%ufft(ic) = 0.0
                  bfft(ic) = 0.0
#endif
               ENDIF
               !-odim
               IF (oneD%odd%d1) THEN  !!!! fmpi version is not tested yet !!!!
                  IF (ic .LT. 9*stars%mx1*stars%mx2 .AND. ic .NE. 0) THEN
                     gx = (cell%bmat(1, 1)*gm(1) + cell%bmat(2, 1)*gm(2))
                     gy = (cell%bmat(1, 2)*gm(1) + cell%bmat(2, 2)*gm(2))
                     gr = SQRT(gx**2 + gy**2)
                     CALL od_cylbes(1, gr*cell%z1, fJ)
#ifdef CPP_MPI
                     ufft_local(ic) = ufft_local(ic) + 2*cell%vol*fJ/(gr*cell%z1*cell%omtil)
#else
                     stars%ufft(ic) = stars%ufft(ic) + 2*cell%vol*fJ/(gr*cell%z1*cell%omtil)
#endif
                  END IF
               END IF
               !+odim
            ENDDO
         ENDDO
      ENDDO

#ifdef CPP_MPI
      CALL MPI_REDUCE(ufft_local, stars%ufft, ifftd, CPP_MPI_REAL, MPI_SUM, 0, fmpi%mpi_comm, ierr)
      CALL MPI_REDUCE(bfft_local, bfft, ifftd, CPP_MPI_REAL, MPI_SUM, 0, fmpi%mpi_comm, ierr)
      CALL MPI_REDUCE(icm_local, icm, size(icm), MPI_INTEGER, MPI_SUM, 0, fmpi%mpi_comm, ierr)
#endif

      IF (fmpi%irank == 0) THEN
         ic = 9*stars%mx1*stars%mx2*(im3 + 1)
         DO i3 = im3 + 1, 3*stars%mx3 - 1
            gm(3) = REAL(i3)
            gm(3) = gm(3) - 3.0*stars%mx3
            DO i2 = 0, 3*stars%mx2 - 1
               gm(2) = REAL(i2)
               IF (2*i2 > 3*stars%mx2) gm(2) = gm(2) - 3.0*stars%mx2
               DO i1 = 0, 3*stars%mx1 - 1
                  gm(1) = REAL(i1)
                  IF (2*i1 > 3*stars%mx1) gm(1) = gm(1) - 3.0*stars%mx1
                  !
                  !-> use inversion <-> c.c.
                  !
                  ic1 = NINT(gm(1)); ic2 = NINT(gm(2)); ic3 = NINT(gm(3))
                  icc = icm(-ic1, -ic2, -ic3)
                  stars%ufft(ic) = stars%ufft(icc)
                  IF (.NOT. sym%invs) bfft(ic) = -bfft(icc)
                  ic = ic + 1
               ENDDO
            ENDDO
         ENDDO
         !
         ! --> add film-contributions
         !
         IF (input%film .AND. .NOT. oneD%odd%d1) THEN

            ifft2d = 9*stars%mx1*stars%mx2
            stars%ufft(0) = stars%ufft(0) + cell%vol*inv_omtil - 1.0

            DO i3 = 1, 3*stars%mx3 - 1
               gm(3) = REAL(i3)
               IF (gm(3) > 1.5*stars%mx3) gm(3) = gm(3) - 3.0*stars%mx3
               th = cell%bmat(3, 3)*gm(3)*cell%z1
               stars%ufft(i3*ifft2d) = stars%ufft(i3*ifft2d) + cell%vol*inv_omtil*SIN(th)/th
            ENDDO

         ELSEIF (oneD%odd%d1) THEN
            !-odim
            stars%ufft(0) = stars%ufft(0) + cell%vol*inv_omtil - 1.0
            !+odim

         ENDIF
         !
         ! --> make fft
         !
         IF (sym%invs) bfft = 0.0
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx1, 3*stars%mx1, +1)
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx2, 9*stars%mx1*stars%mx2, +1)
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx3, ifftd, +1)

         DEALLOCATE (bfft, icm)
#ifdef CPP_MPI
         DEALLOCATE (bfft_local, ufft_local, icm_local)
#endif

         CALL writeStepfunction(stars)
      ENDIF ! (fmpi%irank == 0)

   END SUBROUTINE stepf
END MODULE m_stepf
