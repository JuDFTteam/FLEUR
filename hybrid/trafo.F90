!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_trafo
   use m_judft
   use m_glob_tofrom_loc
   use m_types
   use m_constants
CONTAINS

   SUBROUTINE waveftrafo_symm(cmt_out, z_out, cmt, l_real, z_r, z_c, bandi, ndb, &
                              nk, iop, atoms, mpdata, hybinp, hybdat, kpts, &
                              sym, jsp, lapw)

      USE m_constants
      USE m_wrapper
      USE m_types
      USE m_juDFT
      IMPLICIT NONE

      TYPE(t_mpdata), INTENT(IN)      :: mpdata
      TYPE(t_hybinp), INTENT(IN)      :: hybinp
      TYPE(t_hybdat), INTENT(IN)      :: hybdat
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_lapw), INTENT(IN)        :: lapw

!     - scalars -
      INTEGER, INTENT(IN)      :: nk, jsp, ndb
      INTEGER, INTENT(IN)      ::  bandi, iop

!     - arrays -
      COMPLEX, INTENT(IN)      ::  cmt(:, :, :)
      LOGICAL, INTENT(IN)      ::  l_real
      REAL, INTENT(IN)         ::  z_r(:, :)
      COMPLEX, INTENT(IN)      ::  z_c(:, :)
      COMPLEX, INTENT(INOUT)   ::  cmt_out(hybdat%maxlmindx, atoms%nat, ndb)
      COMPLEX, INTENT(INOUT)   ::  z_out(lapw%nv(jsp), ndb)

!     - local -

!     - scalars -
      INTEGER                 ::  iatom, iatom1, iiatom, itype, igpt, igpt1, ieq, iiop
      INTEGER                 ::  i, l, n, nn, lm0, lm1, lm2
      COMPLEX                 ::  cdum

!     - arrays -
      REAL                    ::  rrot(3, 3), invrrot(3, 3)
      INTEGER                 ::  g(3), g1(3)
      REAL                    ::  tau1(3), rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  cmthlp(2*atoms%lmaxd + 1)
      LOGICAL                 ::  trs

      if (l_real) THEN
         rrot = transpose(1.0*sym%mrot(:, :, sym%invtab(iop)))
         invrrot = transpose(1.0*sym%mrot(:, :, iop))
         trans = sym%tau(:, iop)
      else
         IF (iop <= sym%nop) THEN
            trs = .false.
            rrot = transpose(1.0*sym%mrot(:, :, sym%invtab(iop)))
            invrrot = transpose(1.0*sym%mrot(:, :, iop))
            trans = sym%tau(:, iop)
         ELSE
            trs = .true.
            iiop = iop - sym%nop
            rrot = -transpose(1.0*sym%mrot(:, :, sym%invtab(iiop)))
            invrrot = -transpose(1.0*sym%mrot(:, :, iiop))
            trans = sym%tau(:, iiop)
         END IF
      end if

      rkpt = matmul(rrot, kpts%bkf(:, nk))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g1 = nint(rkpt - rkpthlp)

! MT coefficients
      cmt_out = 0
      iatom = 0
      iiatom = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            iatom1 = hybinp%map(iatom, iop)
            tau1 = hybinp%tvec(:, iatom, iop)

            cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt, tau1))

            lm0 = 0
            DO l = 0, atoms%lmax(itype)
               nn = mpdata%num_radfun_per_l(l, itype)
               DO n = 1, nn
                  lm1 = lm0 + n
                  lm2 = lm0 + n + 2*l*nn
                  DO i = 1, ndb
                     if (l_real) THEN
                        cmt_out(lm1:lm2:nn, iatom1, i) = cdum* &
                                                         matmul(cmt(bandi + i - 1, lm1:lm2:nn, iatom), &
                                                                sym%d_wgn(-l:l, -l:l, l, iop))
                     else
                        IF (trs) THEN
                           cmthlp(:2*l + 1) = CONJG(cmt(bandi + i - 1, lm1:lm2:nn, iatom))
                        ELSE
                           cmthlp(:2*l + 1) = cmt(bandi + i - 1, lm1:lm2:nn, iatom)
                        END IF
                        cmt_out(lm1:lm2:nn, iatom1, i) = cdum*matmul(cmthlp(:2*l + 1), sym%d_wgn(-l:l, -l:l, l, iop))
                     end if
                  END DO
               END DO
               lm0 = lm2
            END DO
         END DO
         iiatom = iiatom + atoms%neq(itype)
      END DO

! PW coefficients
      z_out = 0

      DO igpt = 1, lapw%nv(jsp)
         g = INT(matmul(invrrot, lapw%gvec(:, igpt, jsp) + g1))
!determine number of g
         igpt1 = 0
         DO i = 1, lapw%nv(jsp)
            IF (maxval(abs(g - lapw%gvec(:, i, jsp))) <= 1E-06) THEN
               igpt1 = i
               EXIT
            END IF
         END DO
         IF (igpt1 == 0) THEN
            call judft_error('wavetrafo_symm: rotated G vector not found')
         END IF
         cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt + lapw%gvec(:, igpt, jsp), trans(:)))
         if (l_real) THEN
            z_out(igpt, 1:ndb) = cdum*z_r(igpt1, bandi:bandi + ndb - 1)
         else
            IF (trs) THEN
               z_out(igpt, 1:ndb) = cdum*CONJG(z_c(igpt1, bandi:bandi + ndb - 1))
            ELSE
               z_out(igpt, 1:ndb) = cdum*z_c(igpt1, bandi:bandi + ndb - 1)
            END IF
         end if
      END DO

   END SUBROUTINE waveftrafo_symm

   SUBROUTINE waveftrafo_gen_cmt(cmt, c_phase, l_real, nk, iop, atoms, &
                                 mpdata, hybinp, kpts, sym, nbands, cmt_out)

      use m_juDFT
      USE m_constants
      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      TYPE(t_mpdata), INTENT(IN) :: mpdata
      TYPE(t_hybinp), INTENT(IN) :: hybinp
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)  :: atoms
!     - scalars -
      INTEGER, INTENT(IN)      :: nk, nbands
      INTEGER, INTENT(IN)      ::  iop
      LOGICAL, INTENT(in)      :: l_real
!     - arrays -
      COMPLEX, INTENT(IN)      ::  cmt(:, :, :), c_phase(nbands)

      COMPLEX, INTENT(INOUT)  ::  cmt_out(:, :, :)
!        - local -

!     - scalars -
      INTEGER                 ::  itype, iatom, iatom1, iiatom, ieq, iiop
      INTEGER                 ::  i, l, n, nn, lm0, lm1, lm2
      COMPLEX                 ::  cdum
      LOGICAL                 ::  trs

!     - arrays -
      INTEGER                 ::  rrot(3, 3), invrrot(3, 3)
      INTEGER                 ::  g1(3)
      REAL                    ::  tau1(3), rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  cmthlp(2*atoms%lmaxd + 1)

      call timestart("gen_cmt")
      if (l_real) THEN
         rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
         invrrot = transpose(sym%mrot(:, :, iop))
         trans = sym%tau(:, iop)
      else
         IF (iop <= sym%nop) THEN
            trs = .false.
            rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
            invrrot = transpose(sym%mrot(:, :, iop))
            trans = sym%tau(:, iop)
         ELSE
! in the case of SOC (l_soc=.true.)
! time reversal symmetry is not valid anymore;
! nsym should thus equal nop
            trs = .true.
            iiop = iop - sym%nop
            rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
            invrrot = -transpose(sym%mrot(:, :, iiop))
            trans = sym%tau(:, iiop)
         END IF
      end if

      rkpt = matmul(rrot, kpts%bkf(:, nk))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g1 = nint(rkpt - rkpthlp)

      ! MT coefficients
      cmt_out = 0
      iatom = 0
      iiatom = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            iatom1 = hybinp%map(iatom, iop)
            tau1 = hybinp%tvec(:, iatom, iop)

            cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt, tau1))

            lm0 = 0
            DO l = 0, atoms%lmax(itype)
               nn = mpdata%num_radfun_per_l(l, itype)
               DO n = 1, nn
                  lm1 = lm0 + n
                  lm2 = lm0 + n + 2*l*nn

                  DO i = 1, nbands
                     if (l_real) THEN
                        cmt_out(i, lm1:lm2:nn, iatom1) = cdum*matmul(cmt(i, lm1:lm2:nn, iatom), &
                                                                     hybinp%d_wgn2(-l:l, -l:l, l, iop))
                     else
                        IF (trs) THEN
                           cmthlp(:2*l + 1) = conjg(cmt(i, lm1:lm2:nn, iatom))
                        ELSE
                           cmthlp(:2*l + 1) = cmt(i, lm1:lm2:nn, iatom)
                        END IF
                        cmt_out(i, lm1:lm2:nn, iatom1) = cdum*matmul(cmthlp(:2*l + 1), hybinp%d_wgn2(-l:l, -l:l, l, iop))
                     end if

                  END DO
               END DO
               lm0 = lm2
            END DO
         END DO
         iiatom = iiatom + atoms%neq(itype)
      END DO

      ! If phase and inversion-sym. is true,
      ! define the phase such that z_out is real.
      if (l_real) then
         DO i = 1, nbands
            cmt_out(i, :, :) = cmt_out(i, :, :)/c_phase(i)
         END DO
      end if
      call timestop("gen_cmt")
   END SUBROUTINE waveftrafo_gen_cmt

   SUBROUTINE waveftrafo_genwavf( &
      cmt, z_in, nk, iop, atoms, &
      mpdata, hybinp, kpts, sym, jsp, input, nbands, &
      lapw_nk, lapw_rkpt, cmt_out, z_out)

      use m_juDFT
      USE m_constants
      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      type(t_mat), intent(in)     :: z_in
      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_mpdata), INTENT(IN)    :: mpdata
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_lapw), INTENT(IN)    :: lapw_nk, lapw_rkpt
      type(t_mat), intent(inout)  :: z_out
!     - scalars -
      INTEGER, INTENT(IN)      :: nk, jsp, nbands
      INTEGER, INTENT(IN)      ::  iop
!     - arrays -
      COMPLEX, INTENT(IN)      ::  cmt(:, :, :)

      COMPLEX, INTENT(INOUT)  ::  cmt_out(:, :, :)
!        - local -

!     - scalars -
      INTEGER                 ::  itype, iatom, iatom1, iiatom, igpt, igpt1, ieq, iiop
      INTEGER                 ::  i, l, n, nn, lm0, lm1, lm2
      COMPLEX                 ::  cdum
      LOGICAL                 ::  trs

!     - arrays -
      INTEGER                 ::  rrot(3, 3), invrrot(3, 3)
      INTEGER                 ::  g(3), g1(3)
      REAL                    ::  tau1(3), rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  zhlp(z_in%matsize1, input%neig)
      COMPLEX                 ::  cmthlp(2*atoms%lmaxd + 1)

      call timestart("genwavf")
      if (z_in%l_real) THEN
         rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
         invrrot = transpose(sym%mrot(:, :, iop))
         trans = sym%tau(:, iop)
      else
         IF (iop <= sym%nop) THEN
            trs = .false.
            rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
            invrrot = transpose(sym%mrot(:, :, iop))
            trans = sym%tau(:, iop)
         ELSE
! in the case of SOC (l_soc=.true.)
! time reversal symmetry is not valid anymore;
! nsym should thus equal nop
            trs = .true.
            iiop = iop - sym%nop
            rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
            invrrot = -transpose(sym%mrot(:, :, iiop))
            trans = sym%tau(:, iiop)
         END IF
      end if

      rkpt = matmul(rrot, kpts%bkf(:, nk))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g1 = nint(rkpt - rkpthlp)

      ! MT coefficients
      cmt_out = 0
      iatom = 0
      iiatom = 0

      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            iatom1 = hybinp%map(iatom, iop)
            tau1 = hybinp%tvec(:, iatom, iop)

            cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt, tau1))

            lm0 = 0
            DO l = 0, atoms%lmax(itype)
               nn = mpdata%num_radfun_per_l(l, itype)
               DO n = 1, nn
                  lm1 = lm0 + n
                  lm2 = lm0 + n + 2*l*nn

                  DO i = 1, nbands
                     if (z_in%l_real) THEN
                        cmt_out(i, lm1:lm2:nn, iatom1) = cdum*matmul(cmt(i, lm1:lm2:nn, iatom), &
                                                                     hybinp%d_wgn2(-l:l, -l:l, l, iop))
                     else
                        IF (trs) THEN
                           cmthlp(:2*l + 1) = conjg(cmt(i, lm1:lm2:nn, iatom))
                        ELSE
                           cmthlp(:2*l + 1) = cmt(i, lm1:lm2:nn, iatom)
                        END IF
                        cmt_out(i, lm1:lm2:nn, iatom1) = cdum*matmul(cmthlp(:2*l + 1), hybinp%d_wgn2(-l:l, -l:l, l, iop))
                     end if

                  END DO
               END DO
               lm0 = lm2
            END DO
         END DO
         iiatom = iiatom + atoms%neq(itype)
      END DO

      ! PW coefficients

      zhlp = 0
      DO igpt = 1, lapw_rkpt%nv(jsp)
         g = matmul(invrrot, lapw_rkpt%gvec(:, igpt, jsp) + g1)
         !determine number of g
         igpt1 = 0
         DO i = 1, lapw_nk%nv(jsp)
            IF (all(abs(g - lapw_nk%gvec(:, i, jsp)) <= 1E-06)) THEN
               igpt1 = i
               EXIT
            END IF
         END DO
         IF (igpt1 == 0) CYCLE
         cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt + lapw_rkpt%gvec(:, igpt, jsp), trans))
         if (z_in%l_real) THEN
            zhlp(igpt, :nbands) = cdum*z_in%data_r(igpt1, :nbands)
         else
            IF (trs) THEN
               zhlp(igpt, :nbands) = cdum*conjg(z_in%data_c(igpt1, :nbands))
            ELSE
               zhlp(igpt, :nbands) = cdum*z_in%data_c(igpt1, :nbands)
            END IF
         end if
      END DO

      ! If phase and inversion-sym. is true,
      ! define the phase such that z_out is real.

      DO i = 1, nbands
         if (z_in%l_real) THEN
            cdum = commonphase(zhlp(:, i), z_in%matsize1)

            IF (any(abs(aimag(zhlp(:, i)/cdum)) > 1e-8)) THEN
               WRITE (*, *) maxval(abs(aimag(zhlp(:, i)/cdum)))
               WRITE (*, *) zhlp
               call judft_error('waveftrafo1: Residual imaginary part.')
            END IF
            z_out%data_r(:, i) = real(zhlp(:, i)/cdum)
            cmt_out(i, :, :) = cmt_out(i, :, :)/cdum
         else
            z_out%data_c(:, i) = zhlp(:, i)
         end if
      END DO
      call timestop("genwavf")
   END SUBROUTINE waveftrafo_genwavf

   SUBROUTINE waveftrafo_gen_zmat(z_in, nk, iop, &
                                  kpts, sym, jsp, nbands, &
                                  lapw_nk, lapw_rkpt, z_out, c_phase)

      use m_juDFT
      USE m_constants
      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      type(t_mat), intent(in)     :: z_in
      TYPE(t_sym), INTENT(IN)     :: sym
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_lapw), INTENT(IN)    :: lapw_nk, lapw_rkpt
      type(t_mat), intent(inout)  :: z_out
      complex, intent(inout), optional :: c_phase(:)
!     - scalars -
      INTEGER, INTENT(IN)      :: nk, jsp, nbands
      INTEGER, INTENT(IN)      :: iop

!     - scalars -
      INTEGER                 ::  igpt, igpt1, iiop, i
      COMPLEX                 ::  cdum
      LOGICAL                 ::  trs

!     - arrays -
      INTEGER                 ::  rrot(3, 3), invrrot(3, 3)
      INTEGER                 ::  g(3), g1(3)
      REAL                    ::  rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  zhlp(z_in%matsize1, nbands)

      call timestart("gen_zmat")
      if (present(c_phase)) c_phase = 0

      if (z_in%l_real) THEN
         rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
         invrrot = transpose(sym%mrot(:, :, iop))
         trans = sym%tau(:, iop)
      else
         IF (iop <= sym%nop) THEN
            trs = .false.
            rrot = transpose(sym%mrot(:, :, sym%invtab(iop)))
            invrrot = transpose(sym%mrot(:, :, iop))
            trans = sym%tau(:, iop)
         ELSE
! in the case of SOC (l_soc=.true.)
! time reversal symmetry is not valid anymore;
! nsym should thus equal nop
            trs = .true.
            iiop = iop - sym%nop
            rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
            invrrot = -transpose(sym%mrot(:, :, iiop))
            trans = sym%tau(:, iiop)
         END IF
      end if

      rkpt = matmul(rrot, kpts%bkf(:, nk))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g1 = nint(rkpt - rkpthlp)

      ! PW coefficients

      zhlp = 0
      DO igpt = 1, lapw_rkpt%nv(jsp)
         g = matmul(invrrot, lapw_rkpt%gvec(:, igpt, jsp) + g1)
         !determine number of g
         igpt1 = 0
         DO i = 1, lapw_nk%nv(jsp)
            IF (all(abs(g - lapw_nk%gvec(:, i, jsp)) <= 1E-06)) THEN
               igpt1 = i
               EXIT
            END IF
         END DO
         IF (igpt1 == 0) CYCLE
         cdum = exp(-ImagUnit*tpi_const*dot_product(rkpt + lapw_rkpt%gvec(:, igpt, jsp), trans))
         if (z_in%l_real) THEN
            zhlp(igpt, :nbands) = cdum*z_in%data_r(igpt1, :nbands)
         else
            IF (trs) THEN
               zhlp(igpt, :nbands) = cdum*conjg(z_in%data_c(igpt1, :nbands))
            ELSE
               zhlp(igpt, :nbands) = cdum*z_in%data_c(igpt1, :nbands)
            END IF
         end if
      END DO

      ! If phase and inversion-sym. is true,
      ! define the phase such that z_out is real.

      DO i = 1, nbands
         if (z_in%l_real) THEN
            cdum = commonphase(zhlp(:, i), z_in%matsize1)
            if (present(c_phase)) c_phase(i) = cdum
            if (abs(cdum) < 1e-30) THEN
               call juDFT_error("commonphase can't be 0.")
            end if
            IF (any(abs(aimag(zhlp(:, i)/cdum)) > 1e-8)) THEN
               WRITE (*, *) maxval(abs(aimag(zhlp(:, i)/cdum)))
               call judft_error('waveftrafo1: Residual imaginary part.')
            END IF
            z_out%data_r(:, i) = real(zhlp(:, i)/cdum)
         else
            z_out%data_c(:, i) = zhlp(:, i)
         end if
      END DO
      call timestop("gen_zmat")
   END SUBROUTINE waveftrafo_gen_zmat

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Symmetrizes MT part of input matrix according to inversion symmetry.
   ! This is achieved by a transformation to
   !        1/sqrt(2) * ( exp(ikR) Y_lm(r-R) + (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) )
   ! and                                                                                 if R /=0 or m<0
   !        i/sqrt(2) * ( exp(ikR) Y_lm(r-R) - (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) ) .
   !
   !  or
   !        i*Y_l,0(r)                                                                   if R=0,m=0 and l odd
   ! These functions have the property f(-r)=f(r)* which makes the output matrix real symmetric.
   ! (Array mat is overwritten! )

   SUBROUTINE symmetrize_mpimat(fi, fmpi, mpimat, start_dim, end_dim, imode, lreal, nindxm)
      USE m_types
      use m_constants
      IMPLICIT NONE
      type(t_fleurinput), intent(in)  :: fi
      type(t_mpi), intent(in)         :: fmpi

!     - scalars -
      INTEGER, INTENT(IN)    :: imode, start_dim(2), end_dim(2)
      LOGICAL, INTENT(IN)    ::  lreal

!     - arrays -
      INTEGER, INTENT(IN)    ::  nindxm(0:maxval(fi%hybinp%lcutm1), fi%atoms%ntype)
      COMPLEX, INTENT(INOUT) ::  mpimat(:, :)

!     -local scalars -
      INTEGER               :: i, j, itype, ieq, ic, ic1, l, m, n, nn, ifac, ishift, start_loc, end_loc
      integer               :: i_loc, j_loc, i_pe, j_pe, ierr, len_dim(2)
      REAL                  :: rfac

!     - local arrays -
      COMPLEX               ::  mpicarr(maxval(end_dim)), cfac

      call timestart("symmetrize_mpimat")
      len_dim = end_dim - start_dim + 1
      rfac = sqrt(0.5)
      cfac = sqrt(0.5)*ImagUnit
      ic = 0
      i = 0

      DO itype = 1, fi%atoms%ntype
         nn = sum([((2*l + 1)*nindxm(l, itype), l=0, fi%hybinp%lcutm1(itype))])
         DO ieq = 1, fi%atoms%neq(itype)
            ic = ic + 1
            IF (fi%sym%invsat(ic) == 0) THEN
! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
! invsat(atom) = 0, invsatnr(atom) = 0
! but we need invsatnr(atom) = natom
               ic1 = ic
            ELSE
               ic1 = fi%sym%invsatnr(ic)
            END IF
!ic1 = invsatnr(ic)
            IF (ic1 < ic) THEN
               i = i + nn
               CYCLE
            END IF
!     IF( ic1 .lt. ic ) cycle
            DO l = 0, fi%hybinp%lcutm1(itype)
               ifac = -1
               DO m = -l, l
                  ifac = -ifac
                  ishift = (ic1 - ic)*nn - 2*m*nindxm(l, itype)
                  DO n = 1, nindxm(l, itype)
                     i = i + 1
                     j = i + ishift
                     call glob_to_loc(fmpi, i, i_pe, i_loc)
                     call glob_to_loc(fmpi, j, j_pe, j_loc)
                     IF (ic1 /= ic .or. m < 0) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           call range_from_glob_to_loc(fmpi, start_dim(2), start_loc)
                           call range_to_glob_to_loc(fmpi, end_dim(2), end_loc)
                           mpicarr(start_loc:end_loc)  = mpimat(i,start_loc:end_loc)
                           mpimat(i,start_loc:end_loc) = (mpicarr(start_loc:end_loc) + ifac*mpimat(j,start_loc:end_loc))*rfac
                           mpimat(j,start_loc:end_loc) = (mpicarr(start_loc:end_loc) - ifac*mpimat(j,start_loc:end_loc))*(-cfac)
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           if(i_pe == j_pe .and. fmpi%n_rank == i_pe) then 
                              mpicarr(start_dim(1):end_dim(1)) = mpimat(start_dim(1):end_dim(1), i_loc)
                              mpimat(start_dim(1):end_dim(1),i_loc) &
                                 = (mpimat(start_dim(1):end_dim(1), i_loc) + ifac*mpimat(start_dim(1):end_dim(1), j_loc))*rfac
                              mpimat(start_dim(1):end_dim(1),j_loc) &
                                 = (mpicarr(start_dim(1):end_dim(1)) - ifac*mpimat(start_dim(1):end_dim(1), j_loc))*cfac
#ifdef CPP_MPI
                           else
                              if(fmpi%n_rank == i_pe) then 
                                 call MPI_Send(mpimat(start_dim(1),i_loc), len_dim(1), MPI_DOUBLE_COMPLEX, j_pe, i, fmpi%sub_comm, ierr)
                                 call MPI_Recv(mpicarr(start_dim(1)), len_dim(1), MPI_DOUBLE_COMPLEX, j_pe, j, fmpi%sub_comm, MPI_STATUS_IGNORE, ierr)
                                 mpimat(start_dim(1):end_dim(1),i_loc) = (mpimat(start_dim(1):end_dim(1), i_loc) + ifac*mpicarr(start_dim(1):end_dim(1)))*rfac
                              elseif(fmpi%n_rank == j_pe) then 
                                 call MPI_Recv(mpicarr(start_dim(1)), len_dim(1), MPI_DOUBLE_COMPLEX, i_pe, i, fmpi%sub_comm, MPI_STATUS_IGNORE, ierr)
                                 call MPI_Send(mpimat(start_dim(1),j_loc), len_dim(1), MPI_DOUBLE_COMPLEX, i_pe, j, fmpi%sub_comm, ierr)
                                 mpimat(start_dim(1):end_dim(1),j_loc) = (mpicarr(start_dim(1):end_dim(1)) - ifac*mpimat(start_dim(1):end_dim(1), j_loc))*cfac
                              endif
#endif
                           endif
                        END IF
                     ELSE IF (m == 0 .and. ifac == -1) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           call range_from_glob_to_loc(fmpi, start_dim(2), start_loc)
                           call range_to_glob_to_loc(fmpi, end_dim(2), end_loc)
                           mpimat(i,start_loc:end_loc) = -ImagUnit*mpimat(i,start_loc:end_loc)
                        END IF
                        IF (iand(imode, 2) /= 0 .and. fmpi%n_rank == i_pe) THEN
                           mpimat(start_dim(1):end_dim(1), i_loc) = ImagUnit*mpimat(start_dim(1):end_dim(1), i_loc)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      IF (lreal) THEN
         call juDFT_error("this isn't impemented for mpimat")
! Determine common phase factor and divide by it to make the output matrix real.
         ! cfac = commonphase_mtx(mat, dims(1), dims(2))
         ! do i = 1, dims(1)
         ! do j = 1, dims(2)
         !    mat(i, j) = mat(i, j)/cfac
         !    if (abs(aimag(mat(i, j))) > 1e-8) then
         !       call judft_error('symmetrize: Residual imaginary part. Symmetrization failed.')
         !    end if
         ! end do
         ! end do
      END IF

      call timestop("symmetrize_mpimat")
   END SUBROUTINE symmetrize_mpimat


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ! Symmetrizes MT part of input matrix according to inversion symmetry.
   ! This is achieved by a transformation to
   !        1/sqrt(2) * ( exp(ikR) Y_lm(r-R) + (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) )
   ! and                                                                                 if R /=0 or m<0
   !        i/sqrt(2) * ( exp(ikR) Y_lm(r-R) - (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) ) .
   !
   !  or
   !        i*Y_l,0(r)                                                                   if R=0,m=0 and l odd
   ! These functions have the property f(-r)=f(r)* which makes the output matrix real symmetric.
   ! (Array mat is overwritten! )

   SUBROUTINE symmetrize(mat, dim1, dim2, imode,&
                         atoms, lcutm, maxlcutm, nindxm, sym)
      USE m_types
      use m_constants
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_sym), INTENT(IN)     :: sym

!     - scalars -
      INTEGER, INTENT(IN)    ::  imode, dim1, dim2
      INTEGER, INTENT(IN)    :: maxlcutm

!     - arrays -
      INTEGER, INTENT(IN)    :: lcutm(:)
      INTEGER, INTENT(IN)    ::  nindxm(0:maxlcutm, atoms%ntype)
      COMPLEX, INTENT(INOUT) ::  mat(:, :)

!     -local scalars -
      INTEGER               ::  i, j, itype, ieq, ic, ic1, l, m, n, nn, ifac, ishift
      REAL, parameter       ::  rfac = sqrt(0.5)

!     - local arrays -
      COMPLEX               ::  carr(max(dim1, dim2)), cfac = sqrt(0.5)*ImagUnit

      call timestart("symmetrize")
      ic = 0
      i = 0

      DO itype = 1, atoms%ntype
         nn = sum([((2*l + 1)*nindxm(l, itype), l=0, lcutm(itype))])
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            IF (sym%invsat(ic) == 0) THEN
! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
! invsat(atom) = 0, invsatnr(atom) = 0
! but we need invsatnr(atom) = natom
               ic1 = ic
            ELSE
               ic1 = sym%invsatnr(ic)
            END IF
!ic1 = invsatnr(ic)
            IF (ic1 < ic) THEN
               i = i + nn
               CYCLE
            END IF
!     IF( ic1 .lt. ic ) cycle
            DO l = 0, lcutm(itype)
               ifac = -1
               DO m = -l, l
                  ifac = -ifac
                  ishift = (ic1 - ic)*nn - 2*m*nindxm(l, itype)
                  DO n = 1, nindxm(l, itype)
                     i = i + 1
                     j = i + ishift
                     IF (ic1 /= ic .or. m < 0) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           carr(:dim2) = mat(i, :dim2)
                           mat(i, :dim2) = (carr(:dim2) + ifac*mat(j, :dim2))*rfac
                           mat(j, :dim2) = (carr(:dim2) - ifac*mat(j, :dim2))*(-cfac)
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           carr(:dim1) = mat(:dim1, i)
                           mat(:dim1, i) = (carr(:dim1) + ifac*mat(:dim1, j))*rfac
                           mat(:dim1, j) = (carr(:dim1) - ifac*mat(:dim1, j))*cfac
                        END IF
                     ELSE IF (m == 0 .and. ifac == -1) THEN
                        IF (iand(imode, 1) /= 0) THEN
                           mat(i, :dim2) = -ImagUnit*mat(i, :dim2)
                        END IF
                        IF (iand(imode, 2) /= 0) THEN
                           mat(:dim1, i) = ImagUnit*mat(:dim1, i)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      call timestop("symmetrize")
   END SUBROUTINE symmetrize

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Undoes symmetrization with routine symmetrize.
   SUBROUTINE desymmetrize(mat, dim1, dim2, &
                           atoms, lcutm, maxlcutm, nindxm, sym)

      USE m_types
      IMPLICIT NONE
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars -
      INTEGER, INTENT(IN)      :: dim1, dim2
      INTEGER, INTENT(IN)      :: maxlcutm

!     - arrays -
      INTEGER, INTENT(IN)      :: lcutm(:)
      INTEGER, INTENT(IN)      ::  nindxm(0:maxlcutm, atoms%ntype)
      COMPLEX, INTENT(INOUT)   ::  mat(dim1, dim2)

!     - local scalars -
      INTEGER                 ::  ifac, i, istart, j, itype, ieq, ic, ic1, l, m, n, nn, ishift
      REAL, parameter         ::  rfac1 = sqrt(0.5)
      real                    ::  rfac2
!     - local arrays -
      COMPLEX                 ::  carr(max(dim1, dim2))

      call timestart("desymmetrize")
      ic = 0
      istart = 0
      DO itype = 1, atoms%ntype
         nn = sum([((2*l + 1)*nindxm(l, itype), l=0, lcutm(itype))])
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            ! if the structure is inversion-symmetric, but the equivalent atom belongs to a different unit cell
            ! invsat(atom) = 0, invsatnr(atom) =0
            ! but we need invsatnr(atom) = natom
            ic1 = merge(ic, sym%invsatnr(ic), sym%invsat(ic) == 0)
            !ic1 = invsatnr(ic)
            !IF( ic1 .lt. ic ) cycle
            IF (ic1 < ic) THEN
               istart = istart + nn
            else
               DO l = 0, lcutm(itype)
                  ifac = -1
                  DO m = -l, l
                     ifac = -ifac
                     rfac2 = rfac1*ifac
                     ishift = (ic1 - ic)*nn - 2*m*nindxm(l, itype)
                     IF (ic1 /= ic .or. m < 0) THEN
                        if (ishift <= nindxm(l, itype)) call juDFT_error("if ishift is zero the parallelization is wrong")
                        DO n = 1, nindxm(l, itype)
                           i = istart + n
                           j = i + ishift
                           carr(:dim2) = mat(i, :)
                           mat(i, :) = (carr(:dim2) + ImagUnit*mat(j, :))*rfac1
                           mat(j, :) = (carr(:dim2) - ImagUnit*mat(j, :))*rfac2
                        enddo
                     ELSE IF (m == 0 .and. ifac == -1) THEN
                        DO n = 1, nindxm(l, itype)
                           mat(istart + n, :) = ImagUnit*mat(istart + n, :)
                        enddo
                     endif
                     istart = istart +  nindxm(l, itype)
                  END DO
               END DO
            endif
         END DO
      END DO
      call timestop("desymmetrize")
   END SUBROUTINE desymmetrize

   ! bra_trafo1 rotates cprod at kpts%bkp(ikpt)(<=> not irreducible k-point) to cprod at ikpt (bkp(kpts%bkp(ikpt))), which is the
   ! symmetrie equivalent one
   ! isym maps kpts%bkp(ikpt) on ikpt

   subroutine bra_trafo(fi, mpdata, hybdat, nbands, ikpt, psize, phase, vecin, vecout)
      use m_types
      use m_constants
      use m_judft
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      INTEGER, INTENT(IN)               :: ikpt, nbands, psize
      type(t_mat), INTENT(IN)           :: vecin
      type(t_mat), INTENT(INOUT)        :: vecout
      COMPLEX, INTENT(INOUT)            :: phase(:, :)

      if (vecin%l_real) then
         call bra_trafo_real(fi, mpdata, hybdat, nbands, ikpt, psize, phase, vecin%data_r, vecout%data_r)
      else
         phase = cmplx_1
         call bra_trafo_cmplx(fi, mpdata, hybdat, nbands, ikpt, psize, vecin%data_c, vecout%data_c)
      end if

   end subroutine bra_trafo

   subroutine bra_trafo_real(fi, mpdata, hybdat, nbands, ikpt, psize, phase, matin_r, matout_r)
      use m_types
      use m_constants
      use m_judft
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      INTEGER, INTENT(IN)               :: ikpt, nbands, psize
      REAL, INTENT(IN)                  ::  matin_r(:, :)
      REAL, INTENT(INOUT)               ::  matout_r(:, :)
      COMPLEX, INTENT(INOUT)            ::  phase(:, :)

      COMPLEX, ALLOCATABLE    ::  vecin(:, :), vecout(:, :)
      integer :: ok, i, j, cnt
      integer :: igptm2_list(mpdata%n_g(ikpt))

      phase = cmplx_0
      call timestart("bra trafo real")

      IF (maxval(fi%hybinp%lcutm1) > fi%atoms%lmaxd) call judft_error('bra_trafo: maxlcutm > atoms%lmaxd')   ! very improbable case
      call find_corresponding_g(fi%sym, fi%kpts, mpdata, ikpt, igptm2_list)

!     transform back to unsymmetrized product basis in case of inversion symmetry
      !$OMP parallel default(none) private(i,j, cnt, vecin, vecout, ok) &
      !$OMP shared(nbands, psize, fi, hybdat, mpdata, phase, matin_r, matout_r, ikpt, igptm2_list) 
      allocate (vecin(size(matin_r, dim=1), 1), vecout(size(matin_r, dim=1), 1),  stat=ok, source=cmplx_0)
      IF (ok /= 0) call judft_error('bra_trafo: error allocating vecin or vecout')

      !$OMP do collapse(2)
      DO i = 1, nbands
         DO j = 1, psize
            cnt = (i-1) * psize + j
            vecin(:,1) = matin_r(:,cnt)
            CALL desymmetrize(vecin(:hybdat%nbasp, 1), hybdat%nbasp, 1, &
                              fi%atoms, fi%hybinp%lcutm1, maxval(fi%hybinp%lcutm1), mpdata%num_radbasfn, fi%sym)

            call bra_trafo_core(1, ikpt, 1, fi%sym, mpdata, &
                              fi%hybinp, hybdat, fi%kpts, fi%atoms, igptm2_list, vecin(:,1:1), vecout(:,1:1))

            CALL symmetrize(vecout(:, 1:1), hybdat%nbasm(ikpt), 1, 1, &
                            fi%atoms, fi%hybinp%lcutm1, maxval(fi%hybinp%lcutm1), mpdata%num_radbasfn, fi%sym)

            phase(j, i) = commonphase(vecout(:, 1), hybdat%nbasm(ikpt))
            matout_r(:, cnt) = real(vecout(:, 1)/phase(j, i))
            IF (any(abs(aimag(vecout(:, 1)/phase(j, i))) > 1e-8)) THEN
               WRITE (*, *) vecout(:, 1)/phase(j, i)
               call judft_error('bra_trafo: Residual imaginary part.')
            END IF

         END DO
      END DO
      !$OMP end do

      deallocate (vecout, vecin)
      !$OMP end parallel

      call timestop("bra trafo real")
   end subroutine bra_trafo_real

   subroutine bra_trafo_cmplx(fi, mpdata, hybdat, nbands, ikpt, psize, vecin_c, vecout_c)
      use m_types
      use m_constants
      use m_judft
      implicit none
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(in)        :: mpdata
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      INTEGER, INTENT(IN)               :: ikpt, nbands, psize
      COMPLEX, INTENT(IN)               ::  vecin_c(:, :)
      COMPLEX, INTENT(INOUT)            ::  vecout_c(:, :)

      integer :: igptm2_list(mpdata%n_g(ikpt))

      call timestart("bra trafo cmplx")

      IF (maxval(fi%hybinp%lcutm1) > fi%atoms%lmaxd) call judft_error('bra_trafo: maxlcutm > fi%atoms%lmaxd')   ! very improbable case
      call find_corresponding_g(fi%sym, fi%kpts, mpdata, ikpt, igptm2_list)

      call bra_trafo_core(nbands, ikpt, psize, fi%sym, mpdata, fi%hybinp, hybdat, fi%kpts, fi%atoms, igptm2_list, vecin_c, vecout_c)

      call timestop("bra trafo cmplx")
   end subroutine bra_trafo_cmplx

   subroutine bra_trafo_core(nbands, ikpt, psize, sym, &
                             mpdata, hybinp, hybdat, kpts, atoms, igptm2_list, vecin1, vecout1)
      use m_constants
      implicit none
      type(t_mpdata), intent(in)  :: mpdata
      TYPE(t_hybinp), INTENT(IN)  :: hybinp
      TYPE(t_hybdat), INTENT(IN)  :: hybdat
      TYPE(t_sym), INTENT(IN)     :: sym
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      integer, intent(in)         :: igptm2_list(:)

      INTEGER, INTENT(IN)      ::  ikpt, nbands, psize

      COMPLEX, intent(in)     :: vecin1(:, :)
      complex, intent(inout)  :: vecout1(:, :)

      INTEGER                 :: nrkpt, itype, ic, l, n, i, nn, i1, i2, j1, j2
      INTEGER                 :: igptm, igptm2, igptp, iiop, inviop
      COMPLEX                 :: cexp, cdum

      INTEGER                 :: rrot(3, 3), invrot(3, 3)
      INTEGER                 :: pnt(maxval(mpdata%num_radbasfn), 0:maxval(hybinp%lcutm1), atoms%nat)
      INTEGER                 :: g(3), g1(3)
      REAL                    :: rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 :: dwgn(-maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), &
                                      -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), 0:maxval(hybinp%lcutm1))

      call timestart("bra_trafo_core")
      call timestart("setup")
      IF (kpts%bksym(ikpt) <= sym%nop) THEN
         inviop = sym%invtab(kpts%bksym(ikpt))
         rrot = transpose(sym%mrot(:, :, sym%invtab(kpts%bksym(ikpt))))
         invrot = sym%mrot(:, :, sym%invtab(kpts%bksym(ikpt)))
         trans = sym%tau(:, kpts%bksym(ikpt))

         dwgn(-maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), 0:maxval(hybinp%lcutm1)) &
            = hybinp%d_wgn2(-maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), 0:maxval(hybinp%lcutm1), inviop)

      ELSE
         iiop = kpts%bksym(ikpt) - sym%nop
         inviop = sym%invtab(iiop) + sym%nop
         rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
         invrot = sym%mrot(:, :, sym%invtab(iiop))
         trans = sym%tau(:, iiop)

         dwgn(-maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), 0:maxval(hybinp%lcutm1)) &
            = conjg(hybinp%d_wgn2(-maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1), 0:maxval(hybinp%lcutm1), inviop))
      END IF

      rkpt = matmul(rrot, kpts%bkf(:, kpts%bkp(ikpt)))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g = nint(rkpthlp - rkpt)
      call timestop("setup")

      !test
      call timestart("test")
      nrkpt = 0
      DO i = 1, kpts%nkptf
         IF (maxval(abs(rkpt - kpts%bkf(:, i))) <= 1E-06) THEN
            nrkpt = i
            EXIT
         END IF
      END DO
      IF (nrkpt /= ikpt) THEN
         PRINT *, kpts%bkp(ikpt), ikpt
         PRINT *, kpts%bkf(:, ikpt)
         PRINT *, kpts%bkf(:, kpts%bkp(ikpt))
         PRINT *, rkpt

         call judft_error('bra_trafo: rotation failed')
      END IF
      call timestop("test")

!     Define pointer to first mixed-basis functions (with m = -l)
      call timestart("def pointer to first mpb")
      i = 0
      do ic = 1, atoms%nat
         itype = atoms%itype(ic)
         DO l = 0, hybinp%lcutm1(itype)
            DO n = 1, mpdata%num_radbasfn(l, itype)
               i = i + 1
               pnt(n, l, ic) = i
            END DO
            i = i + mpdata%num_radbasfn(l, itype)*2*l
         END DO
      END DO
      call timestop("def pointer to first mpb")

!     Multiplication
      ! MT
      call timestart("MT part")
      cexp = exp(ImagUnit*tpi_const*dot_product(kpts%bkf(:, ikpt) + g, trans(:)))
      !$OMP parallel do default(none) private(ic, itype, cdum, l, nn, n, i1, i2, j1, j2, i)&
      !$OMP shared(atoms, cexp, hybinp, kpts, mpdata, pnt, dwgn, vecin1, vecout1, ikpt, g, nbands, psize)
      do ic = 1, atoms%nat
         itype = atoms%itype(ic)

         cdum = cexp*exp(-ImagUnit*tpi_const*dot_product(g, atoms%taual(:, hybinp%map(ic, kpts%bksym(ikpt)))))

         DO l = 0, hybinp%lcutm1(itype)
            nn = mpdata%num_radbasfn(l, itype)
            DO n = 1, nn

               i1 = pnt(n, l, ic)
               i2 = i1 + nn*2*l
               j1 = pnt(n, l, hybinp%map(ic, kpts%bksym(ikpt)))
               j2 = j1 + nn*2*l

               DO i = 1, nbands*psize
                  vecout1(i1:i2:nn, i) = cdum*matmul(vecin1(j1:j2:nn, i), dwgn(-l:l, -l:l, l))
               END DO
            END DO
         END DO
      END DO
      !$OMP end parallel do
      call timestop("MT part")

      ! PW
      call timestart("PW part")
      !$OMP parallel do default(none) private(igptm, igptp, g1, igptm2, i, cdum) &
      !$OMP shared(vecout1, vecin1, mpdata, ikpt, igptm2_list, kpts, rrot, g, hybdat, trans, nbands, psize)
      DO igptm = 1, mpdata%n_g(kpts%bkp(ikpt))
         igptp = mpdata%gptm_ptr(igptm, kpts%bkp(ikpt))
         g1 = matmul(rrot, mpdata%g(:, igptp)) + g
         igptm2 = igptm2_list(igptm)                 

         cdum = exp(ImagUnit*tpi_const*dot_product(kpts%bkf(:, ikpt) + g1, trans(:)))
         vecout1(hybdat%nbasp + igptm, :) = cdum*vecin1(hybdat%nbasp + igptm2, :)
      END DO
      !$OMP end parallel do
      call timestop("PW part")
      call timestop("bra_trafo_core")
   end subroutine bra_trafo_core

   subroutine find_corresponding_g(sym, kpts, mpdata, ikpt, igptm2_list)
      implicit none
      type(t_sym), intent(in)    :: sym
      type(t_kpts), intent(in)   :: kpts
      type(t_mpdata), intent(in) :: mpdata
      integer, intent(in)        :: ikpt
      integer, intent(inout)     :: igptm2_list(:)

      integer :: igptm, igptp, g1(3), igptm2, i, iiop
      integer :: g(3), rrot(3, 3)
      REAL    :: rkpt(3), rkpthlp(3)

      call timestart("find correpsonding g")
      call timestart("setup")
      IF (kpts%bksym(ikpt) <= sym%nop) THEN
         rrot = transpose(sym%mrot(:, :, sym%invtab(kpts%bksym(ikpt))))
      ELSE
         iiop = kpts%bksym(ikpt) - sym%nop
         rrot = -transpose(sym%mrot(:, :, sym%invtab(iiop)))
     END IF

      rkpt = matmul(rrot, kpts%bkf(:, kpts%bkp(ikpt)))
      rkpthlp = rkpt
      rkpt = kpts%to_first_bz(rkpt)
      g = nint(rkpthlp - rkpt)
      call timestop("setup")

      !$OMP parallel do default(none) schedule(dynamic, 10) private(igptm, igptp, g1, igptm2) &
      !$OMP shared(kpts, mpdata, ikpt, rrot, g, igptm2_list)
      do igptm = 1, mpdata%n_g(kpts%bkp(ikpt))
         igptp = mpdata%gptm_ptr(igptm, kpts%bkp(ikpt))
         g1 = matmul(rrot, mpdata%g(:, igptp)) + g

         igptm2 = 0
         DO i = 1, mpdata%n_g(ikpt)
            IF (maxval(abs(g1 - mpdata%g(:, mpdata%gptm_ptr(i, ikpt)))) <= 1E-06) THEN
               igptm2 = i
               EXIT
            END IF
         END DO
         IF (igptm2 == 0) THEN
            WRITE (*, *) kpts%bkp(ikpt), ikpt, g1
            WRITE (*, *) mpdata%n_g(kpts%bkp(ikpt)), mpdata%n_g(ikpt)
            WRITE (*, *)
            WRITE (*, *) igptp, mpdata%g(:, igptp)
            WRITE (*, *) g
            WRITE (*, *) rrot
            WRITE (*, *) "Failed tests:", g1
            DO i = 1, mpdata%n_g(ikpt)
               WRITE (*, *) mpdata%g(:, mpdata%gptm_ptr(i, ikpt))
            END DO
            call judft_error('bra_trafo: G-point not found in G-point set.')
         END IF

         igptm2_list(igptm) = igptm2
      enddo
      !$OMP end parallel do
      call timestop("find correpsonding g")
   end subroutine find_corresponding_g

   ! Determines common phase factor (with unit norm)
   function commonphase(carr, n) result(cfac)
      USE m_juDFT
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: n
      COMPLEX, INTENT(IN)      :: carr(n)
      COMPLEX                  :: cfac
      REAL                     :: rdum, rmax
      INTEGER                  :: i

      cfac = 0
      rmax = 0
      DO i = 1, n
         rdum = abs(carr(i))
         IF (rdum > 1e-6) THEN
            cfac = carr(i)/rdum
            EXIT
         ELSE IF (rdum > rmax) THEN
            cfac = carr(i)/rdum
            rmax = rdum
         END IF
      END DO
      IF (abs(cfac) < 1e-10 .and. all(abs(carr) > 1e-10)) THEN
         WRITE (999, *) carr
         call judft_error('commonphase: Could not determine common phase factor. (Wrote carr to fort.999)')
      END IF
   END function commonphase

   function commonphase_mtx(mtx, dim1, dim2) result(cfac)
      implicit none

      COMPLEX, INTENT(IN)      :: mtx(:, :)
      integer, intent(in)      :: dim1, dim2
      COMPLEX                  :: cfac
      REAL                     :: rdum, rmax
      INTEGER                  :: i, j

      do j = 1, dim2
         do i = 1, dim1
            rdum = abs(mtx(i, j))
            IF (rdum > 1e-6) THEN
               cfac = mtx(i, j)/rdum
               EXIT
            ELSE IF (rdum > rmax) THEN
               cfac = mtx(i, j)/rdum
               rmax = rdum
            END IF
         end do
      end do

      IF (abs(cfac) < 1e-10 .and. all(abs(mtx(:dim1, :dim2)) > 1e-10)) THEN
         WRITE (999, *) mtx(:dim1, :dim2)
         call judft_error('commonphase: Could not determine common phase factor. (Wrote carr to fort.999)')
      END IF
   END function commonphase_mtx

   SUBROUTINE bramat_trafo(vecin, igptm_in, ikpt0, iop, writevec, pointer, sym, &
                           rrot, invrrot, mpdata, hybinp, kpts, maxlcutm, atoms, lcutm, nindxm, maxindxm, &
                           dwgn, nbasp, nbasm, vecout, igptm_out)

      USE m_constants
      USE m_util
      USE m_types
      IMPLICIT NONE
      type(t_mpdata), intent(in) :: mpdata
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms

!     - scalars
      INTEGER, INTENT(IN)      ::  ikpt0, igptm_in, iop, maxindxm
      INTEGER, INTENT(IN)      ::  maxlcutm
      INTEGER, INTENT(IN)      ::  nbasp
      LOGICAL, INTENT(IN)      ::  writevec
      INTEGER, INTENT(INOUT)   ::  igptm_out
!     - arrays -
      INTEGER, INTENT(IN)      ::  rrot(:, :), invrrot(:, :)
      INTEGER, INTENT(IN)      :: lcutm(atoms%ntype), &
                                  nindxm(0:maxlcutm, atoms%ntype)
      INTEGER, INTENT(IN)      :: nbasm(:)
      INTEGER, INTENT(IN)      ::  pointer( &
                                  minval(mpdata%g(1, :)) - 1:maxval(mpdata%g(1, :)) + 1, &
                                  minval(mpdata%g(2, :)) - 1:maxval(mpdata%g(2, :)) + 1, &
                                  minval(mpdata%g(3, :)) - 1:maxval(mpdata%g(3, :)) + 1)

      COMPLEX, INTENT(IN)      ::  vecin(:)
      COMPLEX, INTENT(IN)      ::  dwgn(-maxlcutm:maxlcutm, &
                                        -maxlcutm:maxlcutm, &
                                        0:maxlcutm)
      COMPLEX, INTENT(INOUT)     ::  vecout(maxval(nbasm), 1)

!     - private scalars -
      INTEGER                 ::  itype, ieq, ic, l, n, i, nn, i1, i2, j1, j2
      INTEGER                 ::  igptm, igptm2, igptp, isym
      INTEGER                 ::  ikpt1
      LOGICAL                 ::  trs, touch
      COMPLEX                 ::  cexp, cdum
!     - private arrays -
      INTEGER                 ::  pnt(maxindxm, 0:maxlcutm, atoms%nat), g(3), &
                                 g1(3), iarr(maxval(mpdata%n_g))
      REAL                    ::  rkpt(3), rkpthlp(3), trans(3)
      COMPLEX                 ::  vecin1(nbasm(ikpt0))
      COMPLEX                 ::  carr(maxval(mpdata%n_g))

      call timestart("bramat_trafo")
      igptm_out = -1; vecout = CMPLX_NOT_INITALIZED; touch = .false.

      IF (iop <= sym%nop) THEN
         isym = iop
         trs = .false.
         trans = sym%tau(:, isym)
      ELSE
         isym = iop - sym%nop
         trs = .true.
         trans = sym%tau(:, isym)
      END IF

      rkpthlp = matmul(rrot, kpts%bkf(:, ikpt0))
      rkpt = kpts%to_first_bz(rkpthlp)
      g = nint(rkpthlp - rkpt)
      !
      ! determine number of rotated k-point bk(:,ikpt) -> ikpt1
      !
      call timestart("det. kpoint")
      DO i = 1, kpts%nkpt
         IF (maxval(abs(rkpt - kpts%bkf(:, i))) <= 1E-06) THEN
            ikpt1 = i
            EXIT
         END IF
      END DO
      call timestop("det. kpoint")

      call timestart("calc igptm_out")
      DO igptm = 1, mpdata%n_g(ikpt1)
         igptp = mpdata%gptm_ptr(igptm, ikpt1)
         g1 = matmul(invrrot, mpdata%g(:, igptp) - g)
         igptm2 = pointer(g1(1), g1(2), g1(3))
         IF (igptm2 == igptm_in) THEN
            igptm_out = igptm
            touch = .true.
            IF (writevec) THEN
               cdum = exp(ImagUnit*tpi_const*dot_product(kpts%bkf(:, ikpt1) + mpdata%g(:, igptp), trans))
               EXIT
            ELSE
               call timestop("calc igptm_out")
               call timestop("bramat_trafo")
               RETURN
            END IF
         END IF
      END DO
      call timestop("calc igptm_out")

      if (.not. touch) call judft_error("g-point could not be found.")
!     Transform back to unsymmetrized product basis in case of inversion symmetry.

      call timestart("desymm")
      vecout(:nbasm(ikpt0), 1) = vecin(:nbasm(ikpt0))
      if (sym%invs) CALL desymmetrize(vecout, nbasp, 1, &
                                      atoms, lcutm, maxlcutm, nindxm, sym)
      call timestop("desymm")

!     Right-multiplication
      ! PW
      call timestart("right-multiply")
      IF (trs) THEN; vecin1(:nbasm(ikpt0)) = cdum*conjg(vecout(:nbasm(ikpt0), 1))
      ELSE; vecin1(:nbasm(ikpt0)) = cdum*vecout(:nbasm(ikpt0), 1)
      END IF
      call timestop("right-multiply")

!     Define pointer to first mixed-basis functions (with m = -l)
      call timestart("def pointer")
      i = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            DO l = 0, lcutm(itype)
               DO n = 1, nindxm(l, itype)
                  i = i + 1
                  pnt(n, l, ic) = i
               END DO
               i = i + nindxm(l, itype)*2*l
            END DO
         END DO
      END DO
      call timestop("def pointer")

!     Left-multiplication
      ! MT
      call timestart("left multi MT")
      cexp = exp(-ImagUnit*tpi_const*dot_product(kpts%bkf(:, ikpt1) + g, trans))
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            cdum = cexp*exp(ImagUnit*tpi_const*dot_product(g, atoms%taual(:, ic)))
            cdum = conjg(cdum)
            DO l = 0, lcutm(itype)
               nn = nindxm(l, itype)
               DO n = 1, nn

                  i1 = pnt(n, l, ic)
                  i2 = i1 + nn*2*l
                  j1 = pnt(n, l, hybinp%map(ic, sym%invtab(isym)))
                  j2 = j1 + nn*2*l

                  vecout(i1:i2:nn, 1) = cdum*matmul(dwgn(-l:l, -l:l, l), vecin1(j1:j2:nn))

               END DO
            END DO
         END DO
      END DO
      call timestop("left multi MT")

      ! PW
      call timestart("left multi pw")
      DO igptm = 1, mpdata%n_g(ikpt1)
         igptp = mpdata%gptm_ptr(igptm, ikpt1)
         g1 = matmul(invrrot, mpdata%g(:, igptp) - g)
         iarr(igptm) = pointer(g1(1), g1(2), g1(3))
         carr(igptm) = exp(-ImagUnit*tpi_const*dot_product(kpts%bkf(:, ikpt1) + mpdata%g(:, igptp), trans))
      END DO
      DO i1 = 1, mpdata%n_g(ikpt1)
         vecout(nbasp + i1, 1) = carr(i1)*vecin1(nbasp + iarr(i1))
      END DO
      call timestop("left multi pw")

      ! If inversion symmetry is applicable, symmetrize to make the values real.
      call timestart("symmetrize")
      if (sym%invs) CALL symmetrize(vecout, nbasp, 1, 1, &
                                    atoms, lcutm, maxlcutm, nindxm, sym)
      call timestop("symmetrize")
      call timestop("bramat_trafo")
   END SUBROUTINE bramat_trafo

END MODULE m_trafo
