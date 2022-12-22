MODULE m_stepf
   USE m_juDFT
   USE m_cdn_io
#ifdef CPP_MPI
   use mpi
#endif
CONTAINS
   SUBROUTINE stepf(sym, stars, atoms, input, cell, vacuum, fmpi)
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

      USE m_cfft
      USE m_constants
       
      USE m_types
      IMPLICIT NONE
      !     ..
      TYPE(t_sym), INTENT(IN)        :: sym
      TYPE(t_stars), INTENT(INOUT)   :: stars
      TYPE(t_atoms), INTENT(IN)      :: atoms
       
      TYPE(t_input), INTENT(IN)      :: input
      TYPE(t_cell), INTENT(IN)       :: cell
      TYPE(t_vacuum), INTENT(IN)     :: vacuum
      TYPE(t_mpi), INTENT(IN)        :: fmpi

      !     ..
      !     .. Local Scalars ..
      COMPLEX c_c, c_phs, sf
      REAL factorA, dd, gs, th, inv_omtil, r_phs
      REAL g_rmt, g_sqr, help, g_abs, fp_omtil, r_c, gr, gx, gy
      INTEGER i, iType, na, naInit, nn, i1, i2, i3, ic, ifft2d, ifftd, kk
      INTEGER iStar, iStarStart, iStarEnd
      INTEGER ic1, ic2, ic3, icc, im1, im2, im3, loopstart
      LOGICAL l_error
      !     ..
      !     .. Local Arrays ..
      REAL g(3), gm(3), fJ
      REAL, ALLOCATABLE :: bfft(:)
      COMPLEX, ALLOCATABLE :: ustepLocal(:)
      INTEGER, ALLOCATABLE :: icm(:, :, :)
      INTEGER :: i3_start, i3_end, chunk_size, leftover_size
      INTEGER ierr
      INTEGER, ALLOCATABLE :: icm_local(:, :, :)
      REAL, ALLOCATABLE :: ufft_local(:), bfft_local(:)

      ifftd = 27*stars%mx1*stars%mx2*stars%mx3
      !     ..
      !     ..
      !--->    if step function stored on disc, then just read it in
      !
      l_error = .FALSE.
      IF (fmpi%irank == 0) THEN
         CALL readStepfunction(stars, atoms, cell, vacuum, l_error)
      END IF
#ifdef CPP_MPI
      CALL MPI_BCAST(l_error, 1, MPI_LOGICAL, 0, fmpi%mpi_comm, ierr)
#endif
      IF (.NOT. l_error) THEN
         RETURN
      END IF

      ALLOCATE (ustepLocal(stars%ng3))
      stars%ustep(:) = CMPLX(0.0,0.0)
      ustepLocal = CMPLX(0.0,0.0)

      IF (fmpi%irank == 0) THEN
         CALL timestart("ustep")
         IF (input%film) THEN
            dd = vacuum%dvac*cell%area/cell%omtil
         ELSE
            dd = 1.0
         END IF
         !--->    G=0 star
         factorA = 0.0
         DO iType = 1, atoms%ntype
            factorA = factorA + atoms%neq(iType)*atoms%volmts(iType)/cell%omtil
         ENDDO

         ! Treat 1st star separately:

         ustepLocal(1) = CMPLX(dd - factorA, 0.0)

         ! Film contributions:
         !--->    G(parallel)=0  (for film)
         IF (input%film ) THEN
            DO iStar = 2, stars%ng3
               IF (stars%ig2(iStar) .EQ. 1) THEN
                  th = cell%bmat(3, 3)*stars%kv3(3, iStar)*cell%z1
                  ustepLocal(iStar) = CMPLX(cell%vol*SIN(th)/th/cell%omtil, 0.0)
               END IF
            ENDDO
         END IF
      END IF

      ! Now the other stars:

      CALL calcIndexBounds(fmpi,2, stars%ng3, iStarStart, iStarEnd)

#ifdef CPP_MPI
      IF(.NOT.ALLOCATED(stars%sk3)) ALLOCATE(stars%sk3(stars%ng3))
      CALL MPI_BCAST(stars%sk3, stars%ng3, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
      IF(.NOT.ALLOCATED(stars%kv3)) ALLOCATE(stars%kv3(3,stars%ng3))
      CALL MPI_BCAST(stars%kv3, 3*stars%ng3, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
#endif

!$OMP PARALLEL DO DEFAULT(none) & 
!$OMP SHARED(atoms,stars,cell,ustepLocal,iStarStart,iStarEnd) &
!$OMP PRIVATE(iStar,iType,naInit,factorA,th,sf,nn,na,gs)
      DO iStar = iStarStart, iStarEnd
         DO iType = 1, atoms%ntype
            !-->     structure factors: loop over equivalent atoms
            naInit = SUM(atoms%neq(:iType - 1))
            factorA = 3.0 * atoms%volmts(iType) / cell%omtil
            sf = CMPLX(0.0,0.0)

            DO nn = 1, atoms%neq(iType)
               na = naInit + nn
               th = -tpi_const * DOT_PRODUCT(stars%kv3(:,iStar), atoms%taual(:, na))
               sf = sf + CMPLX(COS(th), SIN(th))
            END DO

            gs = stars%sk3(iStar)*atoms%rmt(iType)
            ustepLocal(iStar) = ustepLocal(iStar) - (factorA*(SIN(gs)/gs - COS(gs))/(gs*gs))*sf
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#ifdef CPP_MPI
   ! sum reduce stars%ustep over all MPI ranks
   CALL MPI_REDUCE(ustepLocal, stars%ustep, stars%ng3, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr)
#else
   stars%ustep(:) = ustepLocal(:)
#endif
      DEALLOCATE(ustepLocal)

      IF (fmpi%irank == 0) THEN
         CALL timestop("ustep")
         CALL timestart("Oldstepf")
      END IF ! (fmpi%irank == 0)





      !
      ! --> set up stepfunction on fft-grid:
      !
      ALLOCATE (bfft_local(0:ifftd - 1), ufft_local(0:ifftd - 1))
      bfft_local = 0.0
      ufft_local = 0.0

      ALLOCATE (bfft(0:ifftd - 1))
      im1 = CEILING(1.5*stars%mx1)
      im2 = CEILING(1.5*stars%mx2)
      im3 = CEILING(1.5*stars%mx3)
      ALLOCATE (icm(-im1:im1, -im2:im2, -im3:im3))
      icm = 0
      ALLOCATE (icm_local(-im1:im1, -im2:im2, -im3:im3))
      icm_local = 0

      inv_omtil = 1.0/cell%omtil
      fp_omtil = -fpi_const*inv_omtil
      !DO first vector before loop
      stars%ufft = 0.0
      bfft = 0.0

      IF (fmpi%irank == 0) THEN
         DO iType = 1, atoms%ntype
            ufft_local(0) = ufft_local(0) + atoms%neq(iType)*atoms%volmts(iType)
         ENDDO
         ufft_local(0) = 1.0 - ufft_local(0)*inv_omtil
      ENDIF

      CALL calcIndexBounds(fmpi,0, im3, i3_start, i3_end)

!$OMP PARALLEL DO DEFAULT(none) &
!$OMP SHARED(atoms,stars,cell,sym,bfft_local,ufft_local,icm_local,i3_start,i3_end,im1,im2,im3,fp_omtil,inv_omtil) &
!$OMP PRIVATE(i3,gm,i2,loopstart,i1,ic,na,gs,ic1,ic2,ic3,g,g_sqr,g_abs,help,r_c,iType,r_phs,nn,th,g_rmt,c_c,c_phs)
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

               ic1 = NINT(gm(1))
               ic2 = NINT(gm(2))
               ic3 = NINT(gm(3))
               icm_local(ic1, ic2, ic3) = ic
               IF (ic1 == im1) icm_local(-ic1, ic2, ic3) = ic
               IF (ic2 == im2) icm_local(ic1, -ic2, ic3) = ic
               IF ((ic1 == im1) .AND. (ic2 == im2)) icm_local(-ic1, -ic2, ic3) = ic
               g = MATMUL(TRANSPOSE(cell%bmat), gm)
               g_sqr = DOT_PRODUCT(g, g)
               g_abs = SQRT(g_sqr)
               help = fp_omtil/g_sqr
               IF (sym%invs) THEN
                  r_c = 0.0
                  DO iType = 1, atoms%ntype
                     r_phs = 0.0
                     na = SUM(atoms%neq(:iType - 1))
                     DO nn = 1, atoms%neq(iType)
                        th = -tpi_const*DOT_PRODUCT(gm, atoms%taual(:, na + nn))
                        r_phs = r_phs + COS(th)
                     ENDDO
                     g_rmt = g_abs*atoms%rmt(iType)
                     r_c = r_c + atoms%rmt(iType)*(SIN(g_rmt)/g_rmt - COS(g_rmt))*r_phs
                  ENDDO
                  ufft_local(ic) = help*r_c
               ELSE
                  c_c = CMPLX(0.0, 0.0)
                  DO iType = 1, atoms%ntype
                     c_phs = CMPLX(0.0, 0.0)
                     na = SUM(atoms%neq(:iType - 1))
                     DO nn = 1, atoms%neq(iType)
                        th = -tpi_const*DOT_PRODUCT(gm, atoms%taual(:, na + nn))
                        c_phs = c_phs + EXP(CMPLX(0, th))
                     ENDDO
                     g_rmt = g_abs*atoms%rmt(iType)
                     c_c = c_c + atoms%rmt(iType)*(SIN(g_rmt)/g_rmt - COS(g_rmt))*c_phs
                  ENDDO
                  ufft_local(ic) = help*REAL(c_c)
                  bfft_local(ic) = help*AIMAG(c_c)
               ENDIF

               IF (((2*i3 .EQ. 3*stars%mx3) .OR. (2*i2 .EQ. 3*stars%mx2)) .OR. (2*i1 .EQ. 3*stars%mx1)) THEN
                  ufft_local(ic) = 0.0
                  bfft_local(ic) = 0.0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#ifdef CPP_MPI
      CALL MPI_REDUCE(ufft_local, stars%ufft, ifftd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%mpi_comm, ierr)
      CALL MPI_REDUCE(bfft_local, bfft, ifftd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%mpi_comm, ierr)
      CALL MPI_REDUCE(icm_local, icm, size(icm), MPI_INTEGER, MPI_SUM, 0, fmpi%mpi_comm, ierr)
#else
      stars%ufft(:) = ufft_local(:)
      bfft(:) = bfft_local(:)
      icm(:,:,:) = icm_local(:,:,:)
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
         IF (input%film ) THEN

            ifft2d = 9*stars%mx1*stars%mx2
            stars%ufft(0) = stars%ufft(0) + cell%vol*inv_omtil - 1.0

            DO i3 = 1, 3*stars%mx3 - 1
               gm(3) = REAL(i3)
               IF (gm(3) > 1.5*stars%mx3) gm(3) = gm(3) - 3.0*stars%mx3
               th = cell%bmat(3, 3)*gm(3)*cell%z1
               stars%ufft(i3*ifft2d) = stars%ufft(i3*ifft2d) + cell%vol*inv_omtil*SIN(th)/th
            ENDDO

         ENDIF

         CALL timestop("Oldstepf")
         CALL timestart("FFT - step function")

         !
         ! --> make fft
         !
         IF (sym%invs) bfft = 0.0
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx1, 3*stars%mx1, +1)
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx2, 9*stars%mx1*stars%mx2, +1)
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx3, ifftd, +1)

         DEALLOCATE (bfft, icm)
         DEALLOCATE (bfft_local, ufft_local, icm_local)

         CALL timestop("FFT - step function")

         CALL writeStepfunction(stars)
      ENDIF ! (fmpi%irank == 0)

   END SUBROUTINE stepf
END MODULE m_stepf
