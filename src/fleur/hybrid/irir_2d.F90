!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Block 3a of the IR-IR Coulomb matrix for films: slab-restricted plane waves with the
!>2D kernel.  Blocks 3b/3c are unchanged; the 2D structure constants enter through 3b.
!>The slab-slab term is diagonal in G|| only; overflow-prone exponentials are grouped.
MODULE m_irir_2d
   USE m_juDFT
   USE m_constants
   USE m_irvac_2d, ONLY: irvac_sphere, irvac_sphere_dgamma, irvac_sphere_d2gamma
   USE m_glob_tofrom_loc, ONLY: glob_to_loc
   USE m_types_fleurinput
   USE m_types_hybdat
   USE m_types_mat
   USE m_types_mpdata
   USE m_types_mpi

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: irir_slab_zz, slab_mt_coeff, irir_slab_mt
   PUBLIC :: assemble_irir_3a_film, irir_slab_zz_g0, irir_slab_mt_g0

CONTAINS

   !>int int exp(-iaz) exp(-q|z-z'|) exp(ibz') over [-z1,z1]^2, in closed form.
   PURE COMPLEX FUNCTION irir_slab_zz(a, b, q, z1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a, b, q, z1

      REAL    :: d, zolap, e2
      COMPLEX :: qa, qb

      d = b - a
      IF (ABS(d)*z1 < 1e-10) THEN
         zolap = 2.0*z1
      ELSE
         zolap = 2.0*SIN(d*z1)/d
      END IF
      qa = CMPLX(q, a)
      qb = CMPLX(q, b)
      e2 = EXP(-2.0*q*z1)

      irir_slab_zz = (2.0*q/(q*q + b*b))*zolap &
                     - EXP(-ImagUnit*b*z1)*(EXP(ImagUnit*a*z1) - e2*EXP(-ImagUnit*a*z1))/(qb*qa) &
                     - EXP(ImagUnit*b*z1)*(EXP(-ImagUnit*a*z1) - e2*EXP(ImagUnit*a*z1))/(CONJG(qb)*CONJG(qa))
   END FUNCTION irir_slab_zz

   !>g = 0: int int exp(-iaz) |z-z'| exp(ibz'); inner integral exact, outer by Simpson.
   COMPLEX FUNCTION irir_slab_zz_g0(a, b, z1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a, b, z1

      INTEGER, PARAMETER :: NQ = 2001
      INTEGER :: i
      REAL    :: z, h, w
      COMPLEX :: gz, mz, gend, mend, zolap, acc

      h = 2.0*z1/(NQ - 1)
      gend = cum_g(z1, b, z1)
      mend = cum_m(z1, b, z1)
      acc = cmplx_0
      DO i = 1, NQ
         z = -z1 + (i - 1)*h
         gz = cum_g(z, b, z1)
         mz = cum_m(z, b, z1)
         zolap = 2.0*z*gz - 2.0*mz + mend - z*gend
         w = MERGE(1.0, MERGE(4.0, 2.0, MOD(i, 2) == 0), i == 1 .OR. i == NQ)
         acc = acc + w*EXP(-ImagUnit*a*z)*zolap
      END DO
      irir_slab_zz_g0 = acc*h/3.0
   END FUNCTION irir_slab_zz_g0

   PURE COMPLEX FUNCTION cum_g(z, b, z1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: z, b, z1

      IF (ABS(b)*z1 < 1e-10) THEN
         cum_g = CMPLX(z + z1, 0.0)
      ELSE
         cum_g = (EXP(ImagUnit*b*z) - EXP(-ImagUnit*b*z1))/(ImagUnit*b)
      END IF
   END FUNCTION cum_g

   PURE COMPLEX FUNCTION cum_m(z, b, z1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: z, b, z1

      COMPLEX :: ib

      IF (ABS(b)*z1 < 1e-10) THEN
         cum_m = CMPLX(0.5*(z*z - z1*z1), 0.0)
      ELSE
         ib = ImagUnit*b
         cum_m = EXP(ib*z)*(z/ib - 1.0/(ib*ib)) - EXP(-ib*z1)*(-z1/ib - 1.0/(ib*ib))
      END IF
   END FUNCTION cum_m

   !>The three exponentials of the slab-PW potential at z' = z_a + u_z.
   PURE SUBROUTINE slab_mt_coeff(g1z, g2z, q, z1, za, gam, cf)
      IMPLICIT NONE
      REAL, INTENT(IN)     :: g1z, g2z, q, z1, za
      COMPLEX, INTENT(OUT) :: gam(3), cf(3)

      COMPLEX :: qa

      qa = CMPLX(q, g1z)

      gam(1) = ImagUnit*(g2z - g1z)
      cf(1) = (2.0*q/(q*q + g1z*g1z))*EXP(gam(1)*za)

      gam(2) = CMPLX(-q, g2z)
      cf(2) = -EXP(-q*(z1 + za))*EXP(ImagUnit*(g1z*z1 + g2z*za))/CONJG(qa)

      gam(3) = CMPLX(q, g2z)
      cf(3) = -EXP(-q*(z1 - za))*EXP(-ImagUnit*(g1z*z1 - g2z*za))/qa
   END SUBROUTINE slab_mt_coeff

   !>-(1/omtil) <slab PW|v|MT part of PW>, summed over atoms.
   COMPLEX FUNCTION irir_slab_mt(g1z, g2z, q, z1, omtil, nat, rmt, za, dl, dph)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: g1z, g2z, q, z1, omtil
      INTEGER, INTENT(IN) :: nat
      REAL, INTENT(IN)    :: rmt(:), za(:), dl
      COMPLEX, INTENT(IN) :: dph(:)

      INTEGER :: ia, iterm
      COMPLEX :: gam(3), cf(3), s

      s = cmplx_0
      DO ia = 1, nat
         CALL slab_mt_coeff(g1z, g2z, q, z1, za(ia), gam, cf)
         DO iterm = 1, 3
            s = s + dph(ia)*cf(iterm)*irvac_sphere(dl, gam(iterm), rmt(ia))
         END DO
      END DO
      irir_slab_mt = -(tpi_const/q)*s/omtil
   END FUNCTION irir_slab_mt

   SUBROUTINE assemble_irir_3a_film(fi, mpdata, hybdat, fmpi, ikpt, ngptm1, pgptm1, coul)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      INTEGER, INTENT(IN)            :: ikpt, ngptm1, pgptm1(:)
      CLASS(t_mat), INTENT(INOUT)    :: coul

      REAL, PARAMETER :: QMIN = 1e-10

      INTEGER :: igpt0, igpt1, igpt2, igptp1, igptp2, ix, iy, ix_loc, pe_ix, ia, nat
      REAL    :: q1(3), q2(3), qn1, qn2, dv(2), dl, z1, omtil, dtil
      COMPLEX :: c
      REAL, ALLOCATABLE    :: rmt(:), za(:)
      COMPLEX, ALLOCATABLE :: dph(:)

      CALL timestart("coulomb matrix 3a film")
      nat = fi%atoms%nat
      z1 = fi%cell%z1
      omtil = fi%cell%omtil
      dtil = fi%cell%amat(3, 3)
      ALLOCATE (rmt(nat), za(nat), dph(nat))
      DO ia = 1, nat
         rmt(ia) = fi%atoms%rmt(fi%atoms%itype(ia))
         za(ia) = fi%atoms%pos(3, ia)
      END DO

      DO igpt0 = 1, ngptm1
         igpt2 = pgptm1(igpt0)
         igptp2 = mpdata%gptm_ptr(igpt2, ikpt)
         ix = hybdat%n_mt + igpt2
         CALL glob_to_loc(fmpi, ix, pe_ix, ix_loc)
         IF (fmpi%n_rank /= pe_ix) CYCLE
         q2 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp2), fi%cell%bmat)
         qn2 = NORM2(q2(1:2))

         DO igpt1 = 1, igpt2
            igptp1 = mpdata%gptm_ptr(igpt1, ikpt)
            iy = hybdat%n_mt + igpt1
            q1 = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp1), fi%cell%bmat)
            qn1 = NORM2(q1(1:2))

            dv = q2(1:2) - q1(1:2)
            dl = NORM2(dv)
            DO ia = 1, nat
               dph(ia) = EXP(ImagUnit*(dv(1)*fi%atoms%pos(1, ia) + dv(2)*fi%atoms%pos(2, ia)))
            END DO

            IF (qn1 > QMIN) THEN
               c = irir_slab_mt(q1(3), q2(3), qn1, z1, omtil, nat, rmt, za, dl, dph)
            ELSE
               c = irir_slab_mt_g0(q1(3), q2(3), z1, omtil, nat, rmt, za, dl, dph)
            END IF
            IF (qn2 > QMIN) THEN
               c = c + CONJG(irir_slab_mt(q2(3), q1(3), qn2, z1, omtil, nat, rmt, za, dl, CONJG(dph)))
            ELSE
               c = c + CONJG(irir_slab_mt_g0(q2(3), q1(3), z1, omtil, nat, rmt, za, dl, CONJG(dph)))
            END IF

            ! slab-slab; at q|| = 0 the divergent head is left to m_gamma_2d
            IF (ALL(mpdata%g(1:2, igptp1) == mpdata%g(1:2, igptp2))) THEN
               IF (qn1 > QMIN) THEN
                  c = c + (tpi_const/(qn1*dtil))*irir_slab_zz(q1(3), q2(3), qn1, z1)
               ELSE
                  c = c - (tpi_const/dtil)*irir_slab_zz_g0(q1(3), q2(3), z1)
               END IF
            END IF
            coul%data_c(iy, ix_loc) = c
         END DO
      END DO

      CALL timestop("coulomb matrix 3a film")
   END SUBROUTINE assemble_irir_3a_film

   !>g = 0 counterpart of irir_slab_mt; the potential is linear (G1z = 0: quadratic) in z'.
   COMPLEX FUNCTION irir_slab_mt_g0(g1z, g2z, z1, omtil, nat, rmt, za, dl, dph)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: g1z, g2z, z1, omtil
      INTEGER, INTENT(IN) :: nat
      REAL, INTENT(IN)    :: rmt(:), za(:), dl
      COMPLEX, INTENT(IN) :: dph(:)

      INTEGER :: ia
      REAL    :: a, r
      COMPLEX :: ge, g0, ke, k1, k0, sph, s

      a = g1z
      g0 = ImagUnit*g2z
      s = cmplx_0
      DO ia = 1, nat
         r = rmt(ia)
         IF (ABS(a) < 1e-10) THEN
            sph = EXP(g0*za(ia))*((za(ia)**2 + z1*z1)*irvac_sphere(dl, g0, r) &
                                  + 2.0*za(ia)*irvac_sphere_dgamma(dl, g0, r) &
                                  + irvac_sphere_d2gamma(dl, g0, r))
         ELSE
            ge = ImagUnit*(g2z - a)
            ke = -2.0/(a*a)
            k1 = -2.0*ImagUnit*COS(a*z1)/a
            k0 = 2.0*z1*SIN(a*z1)/a + 2.0*COS(a*z1)/(a*a)
            sph = ke*EXP(ge*za(ia))*irvac_sphere(dl, ge, r) &
                  + EXP(g0*za(ia))*(k1*(za(ia)*irvac_sphere(dl, g0, r) &
                                        + irvac_sphere_dgamma(dl, g0, r)) &
                                    + k0*irvac_sphere(dl, g0, r))
         END IF
         s = s + dph(ia)*sph
      END DO
      irir_slab_mt_g0 = (tpi_const/omtil)*s
   END FUNCTION irir_slab_mt_g0

END MODULE m_irir_2d
