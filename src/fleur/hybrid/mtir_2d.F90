!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Film correction of MT-IR block (2a).  The potential of a slab-restricted plane wave is the
!>bulk 4 pi/Q^2 term plus two null plane waves exp(i Q||.r +- k z), whose projection onto an MT
!>function reduces to moment(n,l) times a solid harmonic, so the correction is additive.
MODULE m_mtir_2d
   USE m_juDFT
   USE m_constants
   USE m_solid_harmonics, ONLY: solidharm_init, solidharm_eval_conj
   USE m_glob_tofrom_loc
   USE m_calc_l_m_from_lm
   USE m_types_fleurinput
   USE m_types_hybdat
   USE m_types_mat
   USE m_types_mpdata
   USE m_types_mpi
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mtir_film_2a_correction

CONTAINS

   PURE REAL FUNCTION dblfac(l)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l
      INTEGER :: j
      dblfac = 1.0
      DO j = 1, 2*l + 1, 2
         dblfac = dblfac*j
      END DO
   END FUNCTION dblfac

   !>Coefficient of moment(n,l) for in-plane |Q||| > 0 (phase and 1/sqrt(omtil) excluded).
   COMPLEX FUNCTION coeff_2a_gpar(l, m, qx, qy, kpar, beta, z1, za)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, m
      REAL, INTENT(IN)    :: qx, qy, kpar, beta, z1, za

      COMPLEX :: kp(3), km(3), cap, cam

      kp = [CMPLX(qx, 0.0), CMPLX(qy, 0.0), CMPLX(0.0, kpar)]
      km = [CMPLX(qx, 0.0), CMPLX(qy, 0.0), CMPLX(0.0, -kpar)]

      cap = EXP(-kpar*(z1 + za))*EXP(-ImagUnit*beta*z1)/CMPLX(kpar, beta)
      cam = EXP(-kpar*(z1 - za))*EXP(ImagUnit*beta*z1)/CMPLX(kpar, -beta)

      coeff_2a_gpar = -(tpi_const/kpar)*(fpi_const/dblfac(l))*ImagUnit**l &
                      *(cap*solidharm_eval_conj(l, m, kp) + cam*solidharm_eval_conj(l, m, km))
   END FUNCTION coeff_2a_gpar

   !>Q|| = 0, Q_z /= 0: only l = 0 and (1,0) change; the divergent 2 pi/g head is omitted.
   COMPLEX FUNCTION coeff_2a_g0(l, m, beta, z1, za)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, m
      REAL, INTENT(IN)    :: beta, z1, za

      IF (m /= 0) THEN
         coeff_2a_g0 = cmplx_0
      ELSE IF (l == 0) THEN
         coeff_2a_g0 = -tpi_const*sfp_const &
                       *(2.0*z1*SIN(beta*z1)/beta + 2.0*COS(beta*z1)/beta**2 &
                         + 2.0*ImagUnit*COS(beta*z1)*za/beta)
      ELSE IF (l == 1) THEN
         coeff_2a_g0 = -tpi_const*SQRT(fpi_const/3.0)*(2.0*ImagUnit*COS(beta*z1)/beta)
      ELSE
         coeff_2a_g0 = cmplx_0
      END IF
   END FUNCTION coeff_2a_g0

   !>G = 0 at Gamma: potential -2 pi (z^2 + z1^2); the bulk moment2 term stays as it is.
   COMPLEX FUNCTION coeff_2a_gamma(l, m, z1, za)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, m
      REAL, INTENT(IN)    :: z1, za

      IF (m /= 0) THEN
         coeff_2a_gamma = cmplx_0
      ELSE IF (l == 0) THEN
         coeff_2a_gamma = -tpi_const*sfp_const*(za**2 + z1**2)
      ELSE IF (l == 1) THEN
         coeff_2a_gamma = -tpi_const*SQRT(fpi_const/3.0)*2.0*za
      ELSE IF (l == 2) THEN
         coeff_2a_gamma = -tpi_const*SQRT(fpi_const/5.0)*(2.0/3.0)
      ELSE
         coeff_2a_gamma = cmplx_0
      END IF
   END FUNCTION coeff_2a_gamma

   SUBROUTINE mtir_film_2a_correction(fi, hybdat, mpdata, fmpi, moment, pgptm1, ngptm1, ikpt, coul)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      REAL, INTENT(IN)               :: moment(:, 0:, :)
      INTEGER, INTENT(IN)            :: pgptm1(:), ngptm1, ikpt
      CLASS(t_mat), INTENT(INOUT)    :: coul

      REAL, PARAMETER :: QMIN = 1e-10

      INTEGER :: igpt0, igpt, igptp, ix, ix_loc, pe_ix, ic, itype, lm, l, m, n, iy, lmx
      REAL    :: q(3), kc(3), ra(3), kpar, beta, z1, svol
      COMPLEX :: cphase, cval
      LOGICAL :: l_gamma

      CALL timestart("MT-IR 2a film")

      z1 = fi%cell%z1
      svol = SQRT(fi%cell%omtil)
      lmx = MAXVAL(fi%hybinp%lcutm1)
      CALL solidharm_init(lmx)

      kc = MATMUL(fi%kpts%bk(:, ikpt), fi%cell%bmat)
      IF (ABS(kc(3)) > 1e-10) CALL juDFT_error( &
         "film hybrid: k-point with a finite k_z", calledby="mtir_film_2a_correction", &
         hint="the MT-IR film branch assumes a 2D k-mesh")

      DO igpt0 = 1, ngptm1
         igpt = pgptm1(igpt0)
         igptp = mpdata%gptm_ptr(igpt, ikpt)
         ix = hybdat%n_mt + igpt
         CALL glob_to_loc(fmpi, ix, pe_ix, ix_loc)
         IF (pe_ix /= fmpi%n_rank) CYCLE

         q = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)
         kpar = NORM2(q(1:2))
         beta = q(3)

         l_gamma = kpar <= QMIN .AND. ABS(beta) <= QMIN

         iy = 0
         DO ic = 1, fi%atoms%nat
            itype = fi%atoms%itype(ic)
            ra = fi%atoms%pos(:, ic)
            cphase = EXP(ImagUnit*((q(1) - kc(1))*ra(1) + (q(2) - kc(2))*ra(2)))

            DO lm = 1, (fi%hybinp%lcutm1(itype) + 1)**2
               CALL calc_l_m_from_lm(lm, l, m)

               IF (l_gamma) THEN
                  cval = coeff_2a_gamma(l, m, z1, ra(3))
               ELSE IF (kpar > QMIN) THEN
                  cval = coeff_2a_gpar(l, m, q(1), q(2), kpar, beta, z1, ra(3))
               ELSE
                  cval = coeff_2a_g0(l, m, beta, z1, ra(3))
               END IF
               cval = cval*cphase/svol

               DO n = 1, mpdata%num_radbasfn(l, itype)
                  iy = iy + 1
                  IF (cval /= cmplx_0) THEN
                     coul%data_c(iy, ix_loc) = coul%data_c(iy, ix_loc) + cval*moment(n, l, itype)
                  END IF
               END DO
            END DO
         END DO
      END DO

      CALL timestop("MT-IR 2a film")
   END SUBROUTINE mtir_film_2a_correction

END MODULE m_mtir_2d
