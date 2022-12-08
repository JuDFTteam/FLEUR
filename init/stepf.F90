MODULE m_stepf
   USE m_juDFT
   USE m_cdn_io
#ifdef CPP_MPI
   use mpi
#endif
CONTAINS
   SUBROUTINE stepf(sym, stars, atoms, input, cell, vacuum, fmpi, qvec, iDtype, iDir)
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

      REAL,    OPTIONAL, INTENT(IN) :: qvec(3)
      INTEGER, OPTIONAL, INTENT(IN) :: iDtype, iDir

      !     ..
      !     .. Local Scalars ..
      COMPLEX c_c, c_phs, sf
      REAL factorA, dd, gs, th, inv_omtil, r_phs
      REAL g_rmt, g_sqr, help, g_abs, fp_omtil, r_c, gr, gx, gy
      INTEGER i, iType, na, naInit, nn, i1, i2, i3, ic, ifft2d, ifftd, kk
      INTEGER iStar, iStarStart, iStarEnd
      INTEGER ic1, ic2, ic3, icc, im1, im2, im3, loopstart
      LOGICAL l_error, presentQvec
      !     ..
      !     .. Local Arrays ..
      REAL g(3), gm(3), fJ
      REAL, ALLOCATABLE :: bfft(:)
      INTEGER, ALLOCATABLE :: icm(:, :, :)
      INTEGER :: i3_start, i3_end, chunk_size, leftover_size
#ifdef CPP_MPI
      INTEGER ierr
      INTEGER, ALLOCATABLE :: icm_local(:, :, :)
      REAL, ALLOCATABLE :: ufft_local(:), bfft_local(:)
#endif

      ifftd = 27*stars%mx1*stars%mx2*stars%mx3
      !     ..
      !     ..
      !--->    if step function stored on disc, then just read it in
      !
      l_error = .FALSE.
      IF (fmpi%irank == 0 .AND. .NOT.PRESENT(qvec)) THEN
         CALL readStepfunction(stars, atoms, cell, vacuum, l_error)
      ELSE IF (PRESENT(qvec)) THEN
         l_error = .TRUE.
      END IF
#ifdef CPP_MPI
      CALL MPI_BCAST(l_error, 1, MPI_LOGICAL, 0, fmpi%mpi_comm, ierr)
#endif
      IF (.NOT. l_error) THEN
         RETURN
      END IF

      stars%ustep(:) = CMPLX(0.0,0.0)
      presentQvec = PRESENT(qvec)

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

         IF (.NOT.presentQvec) stars%ustep(1) = CMPLX(dd - factorA, 0.0)
         !--->    G(parallel)=0  (for film)
         IF (input%film ) THEN
            DO iStar = 2, stars%ng3
               IF (stars%ig2(iStar) .EQ. 1) THEN
                  th = cell%bmat(3, 3)*stars%kv3(3, iStar)*cell%z1
                  stars%ustep(iStar) = CMPLX(cell%vol*SIN(th)/th/cell%omtil, 0.0)
               END IF
            ENDDO
         END IF

         !--->    sphere contributions

         ! Treat 1st star separately:

         iStar = 1
         DO iType = MERGE(1,iDtype,.NOT.presentQvec), MERGE(atoms%ntype,iDtype,.NOT.presentQvec)
            !-->     structure factors: loop over equivalent atoms
            naInit = SUM(atoms%neq(:iType - 1)) + 1
            factorA = 3.*atoms%volmts(iType)/cell%omtil

            IF(norm2(stars%center).GT.1.0e-8) THEN
               th = -tpi_const*DOT_PRODUCT(stars%kv3(:,iStar)+stars%center, atoms%taual(:, naInit))
               sf = MERGE(-ImagUnit*stars%gq(iDir,iStar)*CMPLX(COS(th), SIN(th)), CMPLX(COS(th), SIN(th)), presentQvec)
               gs = stars%sk3(iStar)*atoms%rmt(iType)
               stars%ustep(iStar) = stars%ustep(iStar) - (factorA*(SIN(gs)/gs - COS(gs))/(gs*gs))*sf
            END IF
         END DO
      END IF

      iStarStart = 2
      iStarEnd = stars%ng3
!#ifdef CPP_MPI
!      !calculate loop boundaries for the MPI parallelization
!      chunk_size = (stars%ng3-1) / fmpi%isize
!#endif
      IF (fmpi%irank == 0) THEN
         ! Now the other stars:
!$OMP PARALLEL DO DEFAULT(none) & 
!$OMP SHARED(atoms,stars,cell,iStarStart,iStarEnd,iDType,iDir,presentQvec) &
!$OMP PRIVATE(iStar,iType,naInit,factorA,th,sf,nn,na,gs)
         DO iStar = iStarStart, iStarEnd
            DO iType = MERGE(1,iDtype,.NOT.presentQvec), MERGE(atoms%ntype,iDtype,.NOT.presentQvec)
               !-->     structure factors: loop over equivalent atoms
               naInit = SUM(atoms%neq(:iType - 1)) + 1

               factorA = 3.0 * atoms%volmts(iType) / cell%omtil

               th = -tpi_const*DOT_PRODUCT(stars%kv3(:,iStar)+stars%center, atoms%taual(:, naInit))
               sf = MERGE(-ImagUnit*stars%gq(iDir,iStar)*CMPLX(COS(th), SIN(th)), CMPLX(COS(th), SIN(th)), presentQvec)

               DO nn = 2, atoms%neq(iType) !Should automatically be nn = 2, 1 for DFPT
                  na = naInit + nn - 1
                  th = -tpi_const * DOT_PRODUCT(stars%kv3(:,iStar), atoms%taual(:, na))
                  sf = sf + CMPLX(COS(th), SIN(th))
               END DO

               gs = stars%sk3(iStar)*atoms%rmt(iType)
               stars%ustep(iStar) = stars%ustep(iStar) - (factorA*(SIN(gs)/gs - COS(gs))/(gs*gs))*sf
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
!#ifdef CPP_MPI
!   ! sum reduce stars%ustep over all MPI ranks
!#endif
      END IF
      IF (fmpi%irank == 0) THEN
         CALL timestop("ustep")
         CALL timestart("Oldstepf")
      END IF ! (fmpi%irank == 0)



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
         DO iType = 1, atoms%ntype
            ufft_local(0) = ufft_local(0) + atoms%neq(iType)*atoms%volmts(iType)
         ENDDO
         ufft_local(0) = 1.0 - ufft_local(0)*inv_omtil
      ENDIF
#else
      DO n = iType, atoms%ntype
         stars%ufft(0) = stars%ufft(0) + atoms%neq(iType)*atoms%volmts(iType)
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
               ! TODO: (How) do these conditions apply for DFPT?
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
               IF (PRESENT(qvec)) g = g + MATMUL(TRANSPOSE(cell%bmat), qvec)
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
#ifdef CPP_MPI
                  ufft_local(ic) = help*r_c
#else
                  stars%ufft(ic) = help*r_c
#endif
               ELSE
                  c_c = CMPLX(0.0, 0.0)
                  DO iType = 1, MERGE(atoms%ntype,iDtype,.NOT.PRESENT(qvec))
                     c_phs = CMPLX(0.0, 0.0)
                     na = MERGE(SUM(atoms%neq(:iType - 1)),iDtype,.NOT.PRESENT(qvec))
                     DO nn = 1, atoms%neq(iType)
                        th = -tpi_const*DOT_PRODUCT(gm, atoms%taual(:, na + nn))
                        IF (.NOT.PRESENT(qvec)) THEN
                           c_phs = c_phs + EXP(CMPLX(0, th))
                        ELSE
                           c_phs = c_phs + (-ImagUnit*g(iDir))*EXP(CMPLX(0, th))
                        END IF
                     ENDDO
                     g_rmt = g_abs*atoms%rmt(iType)
                     c_c = c_c + atoms%rmt(iType)*(SIN(g_rmt)/g_rmt - COS(g_rmt))*c_phs
                  ENDDO
#ifdef CPP_MPI
                  ufft_local(ic) = help*REAL(c_c)
                  bfft_local(ic) = help*AIMAG(c_c)
#else
                  stars%ufft(ic) = help*REAL(c_c)
                  bfft(ic) = help*AIMAG(c_c)
#endif
               ENDIF

               ! TODO: (How)do these conditions apply for DFPT?
               IF (((2*i3 .EQ. 3*stars%mx3) .OR. (2*i2 .EQ. 3*stars%mx2)) .OR. (2*i1 .EQ. 3*stars%mx1)) THEN
               !IF (.FALSE.) THEN !NEWER; BREAKS A TEST (FePt LO)
#ifdef CPP_MPI
                  ufft_local(ic) = 0.0
                  bfft_local(ic) = 0.0
#else
                  stars%ufft(ic) = 0.0
                  bfft(ic) = 0.0
#endif
               ENDIF
\
            ENDDO
         ENDDO
      ENDDO

#ifdef CPP_MPI
      CALL MPI_REDUCE(ufft_local, stars%ufft, ifftd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%mpi_comm, ierr)
      CALL MPI_REDUCE(bfft_local, bfft, ifftd, MPI_DOUBLE_PRECISION, MPI_SUM, 0, fmpi%mpi_comm, ierr)
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

         !
         ! --> make fft
         !
         IF (sym%invs) bfft = 0.0
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx1, 3*stars%mx1, +1)
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx2, 9*stars%mx1*stars%mx2, +1)
         CALL cfft(stars%ufft, bfft, ifftd, 3*stars%mx3, ifftd, +1)

         IF (PRESENT(qvec)) stars%ufft1 = CMPLX(stars%ufft,bfft)

         DEALLOCATE (bfft, icm)
#ifdef CPP_MPI
         DEALLOCATE (bfft_local, ufft_local, icm_local)
#endif

         IF (.NOT.PRESENT(qvec)) CALL writeStepfunction(stars)
      ENDIF ! (fmpi%irank == 0)

   END SUBROUTINE stepf
END MODULE m_stepf
