!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>IR-VAC block of the film Coulomb matrix.  Inside the slab the potential of a vacuum
!>function is a null plane wave, (2 pi/k) mu exp(i k.r||) exp(-k(z1 - sigma z)), so slab and
!>sphere integrals are closed form.  Only the moment carriers couple (at g = 0 two of them).
MODULE m_irvac_2d
   USE m_vac_rows, ONLY: NVAC_MPB
   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_mtir_size
   USE m_coulomb_vac, ONLY: vac_exp_mom, vac_mom_g0
   USE m_vac_rows, ONLY: row_offset, basfn_offset
   USE m_types_mpimat
   USE m_glob_tofrom_loc, ONLY: glob_to_loc
#ifdef CPP_MPI
   USE mpi
#endif

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: irvac_slab_z, irvac_sphere, irvac_sphere_dgamma, irvac_sphere_d2gamma
   PUBLIC :: irvac_g0_slab, irvac_g0_mt
   PUBLIC :: assemble_irvac_dense, copy_irvac_to_sparse

CONTAINS

   !>One element of the dense matrix at global (irow, icol); columns are distributed.
   SUBROUTINE put_dense(coulomb, fmpi, irow, icol, x)
      IMPLICIT NONE
      CLASS(t_mat), INTENT(INOUT) :: coulomb
      TYPE(t_mpi), INTENT(IN)     :: fmpi
      INTEGER, INTENT(IN)         :: irow, icol
      COMPLEX, INTENT(IN)         :: x

      INTEGER :: pe, jloc

      SELECT TYPE (coulomb)
      CLASS IS (t_mpimat)
         CALL glob_to_loc(fmpi, icol, pe, jloc)
         IF (fmpi%n_rank /= pe) RETURN
         coulomb%data_c(irow, jloc) = x
      CLASS IS (t_mat)
         IF (coulomb%l_real) THEN
            coulomb%data_r(irow, icol) = REAL(x)
         ELSE
            coulomb%data_c(irow, icol) = x
         END IF
      END SELECT
   END SUBROUTINE put_dense

   !>g(x) = (sin x - x cos x)/x^3, series for small |x|.
   PURE COMPLEX FUNCTION sph_g(x)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: x

      COMPLEX :: x2, term
      INTEGER :: n

      IF (ABS(x) < 0.2) THEN
         x2 = x*x
         sph_g = cmplx_0
         term = CMPLX(1.0, 0.0)
         DO n = 1, 8
            sph_g = sph_g + ((-1)**(n + 1))*(2*n)*term/factorial(2*n + 1)
            term = term*x2
         END DO
      ELSE
         sph_g = (SIN(x) - x*COS(x))/x**3
      END IF
   END FUNCTION sph_g

   PURE REAL FUNCTION factorial(n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER :: i
      factorial = 1.0
      DO i = 2, n
         factorial = factorial*i
      END DO
   END FUNCTION factorial

   PURE COMPLEX FUNCTION irvac_slab_z(alpha, z1)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: alpha
      REAL, INTENT(IN)    :: z1

      IF (ABS(alpha) < 1e-10) THEN
         irvac_slab_z = CMPLX(2.0*z1, 0.0)
      ELSE
         irvac_slab_z = 2.0*SINH(alpha*z1)/alpha
      END IF
   END FUNCTION irvac_slab_z

   !>int_sphere exp(i Delta.u|| + alpha u_z) d^3u.
   PURE COMPLEX FUNCTION irvac_sphere(dl, alpha, r)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: dl
      COMPLEX, INTENT(IN) :: alpha
      REAL, INTENT(IN)    :: r

      COMPLEX :: kap

      kap = SQRT(CMPLX(dl*dl, 0.0) - alpha*alpha)
      irvac_sphere = fpi_const*r**3*sph_g(kap*r)
   END FUNCTION irvac_sphere

   PURE COMPLEX FUNCTION sph_dg_ox(x)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: x

      COMPLEX :: x2

      IF (ABS(x) < 0.2) THEN
         x2 = x*x
         sph_dg_ox = -1.0/15.0 + x2/210.0 - x2*x2/7560.0 + x2*x2*x2/498960.0
      ELSE
         sph_dg_ox = SIN(x)/x**3 - 3.0*(SIN(x) - x*COS(x))/x**5
      END IF
   END FUNCTION sph_dg_ox

   PURE COMPLEX FUNCTION irvac_sphere_dgamma(dl, gam, r)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: dl
      COMPLEX, INTENT(IN) :: gam
      REAL, INTENT(IN)    :: r

      COMPLEX :: kap

      kap = SQRT(CMPLX(dl*dl, 0.0) - gam*gam)
      irvac_sphere_dgamma = -fpi_const*r**5*gam*sph_dg_ox(kap*r)
   END FUNCTION irvac_sphere_dgamma

   PURE COMPLEX FUNCTION sph_d2g_ox(x)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: x

      COMPLEX :: x2

      IF (ABS(x) < 0.2) THEN
         x2 = x*x
         sph_d2g_ox = 1.0/105.0 - x2/1890.0 + x2*x2/83160.0 - x2*x2*x2/6486480.0
      ELSE
         sph_d2g_ox = COS(x)/x**4 - 6.0*SIN(x)/x**5 + 15.0*(SIN(x) - x*COS(x))/x**7
      END IF
   END FUNCTION sph_d2g_ox

   PURE COMPLEX FUNCTION irvac_sphere_d2gamma(dl, gam, r)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: dl
      COMPLEX, INTENT(IN) :: gam
      REAL, INTENT(IN)    :: r

      COMPLEX :: kap

      kap = SQRT(CMPLX(dl*dl, 0.0) - gam*gam)
      irvac_sphere_d2gamma = -fpi_const*r**5*(sph_dg_ox(kap*r) - gam*gam*r*r*sph_d2g_ox(kap*r))
   END FUNCTION irvac_sphere_d2gamma

   !>IR-VAC into the dense matrix, before apply_inverse_olaps: the IR index needs O^-1.
   SUBROUTINE assemble_irvac_dense(fi, mpdata, hybdat, fmpi, coulomb, ikpt)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      CLASS(t_mat), INTENT(INOUT)    :: coulomb
      INTEGER, INTENT(IN)            :: ikpt

      INTEGER :: ivac, igm, ilen, nn, igpt, igptp, ia, irow, icol, nz
      REAL    :: kvac(3), kir(3), dl2, kk, sigma, z1, dtil, mu, delz, dv(2), gzc, area
      REAL    :: q0v, q1v
      INTEGER :: islot, nslot, nsl
      LOGICAL :: l_g0
      COMPLEX :: alpha, cc, x, ssum

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("IR-VAC coulomb blocks")
      z1 = fi%cell%z1
      dtil = fi%cell%amat(3, 3)
      area = fi%cell%omtil/dtil
      nz = mpdata%nmz_vac
      delz = fi%vacuum%delz

      DO ivac = 1, NVAC_MPB
         sigma = MERGE(1.0, -1.0, ivac == 1)
         DO igm = 1, mpdata%n_g_vac(ikpt)
            ilen = mpdata%glen_ptr_vac(igm, ikpt)
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            kk = mpdata%glen_vac(ilen)
            l_g0 = (kk < 1e-12)
            kvac = MATMUL(fi%kpts%bk(:, ikpt) &
                          + [REAL(mpdata%g_vac(1, mpdata%gptm_ptr_vac(igm, ikpt))), &
                             REAL(mpdata%g_vac(2, mpdata%gptm_ptr_vac(igm, ikpt))), 0.0], &
                          fi%cell%bmat)
            nslot = MERGE(2, 1, l_g0 .AND. nn >= 2)
            DO islot = 1, nslot
               nsl = MERGE(nn, nn - 1, islot == 1)
               IF (l_g0) THEN
                  CALL vac_mom_g0(mpdata%zbasfn_vac(:, nsl, ilen, ivac), nz, delz, q0v, q1v)
                  cc = -tpi_const/SQRT(dtil)
               ELSE
                  mu = vac_exp_mom(mpdata%zbasfn_vac(:, nsl, ilen, ivac), nz, delz, kk)
                  cc = (tpi_const/kk)*mu*EXP(-kk*z1)/SQRT(dtil)
               END IF
               cc = cc*SQRT(area)
               icol = row_offset(mpdata, hybdat, ikpt, ivac) &
                      + basfn_offset(mpdata, ikpt, ivac, igm) + nsl

               DO igpt = 1, mpdata%n_g(ikpt)
                  igptp = mpdata%gptm_ptr(igpt, ikpt)
                  kir = MATMUL(fi%kpts%bk(:, ikpt) + mpdata%g(:, igptp), fi%cell%bmat)
                  dv = kvac(1:2) - kir(1:2)
                  dl2 = SQRT(dv(1)**2 + dv(2)**2)
                  gzc = kir(3)
                  alpha = CMPLX(sigma*kk, -gzc)

                  ssum = cmplx_0
                  IF (l_g0) THEN
                     IF (dl2 < 1e-10) ssum = irvac_g0_slab(gzc, z1, sigma, q0v, q1v)
                  ELSE IF (dl2 < 1e-10) THEN
                     ssum = irvac_slab_z(alpha, z1)
                  END IF

                  DO ia = 1, fi%atoms%nat
                     IF (l_g0) THEN
                        ssum = ssum - EXP(ImagUnit*(dv(1)*fi%atoms%pos(1, ia) &
                                                    + dv(2)*fi%atoms%pos(2, ia))) &
                               *EXP(CMPLX(0.0, -gzc)*fi%atoms%pos(3, ia)) &
                               *irvac_g0_mt(gzc, z1, sigma, fi%atoms%pos(3, ia), q0v, q1v, &
                                            dl2, fi%atoms%rmt(fi%atoms%itype(ia)))/area
                     ELSE
                        ssum = ssum - EXP(ImagUnit*(dv(1)*fi%atoms%pos(1, ia) &
                                                    + dv(2)*fi%atoms%pos(2, ia))) &
                               *EXP(alpha*fi%atoms%pos(3, ia)) &
                               *irvac_sphere(dl2, alpha, fi%atoms%rmt(fi%atoms%itype(ia)))/area
                     END IF
                  END DO

                  x = cc*ssum
                  irow = hybdat%n_mt + igpt
                  CALL put_dense(coulomb, fmpi, irow, icol, x)
                  CALL put_dense(coulomb, fmpi, icol, irow, CONJG(x))
               END DO
            END DO
         END DO
      END DO
      CALL timestop("IR-VAC coulomb blocks")
   END SUBROUTINE assemble_irvac_dense

   SUBROUTINE copy_irvac_to_sparse(fi, mpdata, hybdat, fmpi, coulomb, ikpt)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)  :: hybdat
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      CLASS(t_mat), INTENT(IN)       :: coulomb
      INTEGER, INTENT(IN)            :: ikpt

      INTEGER :: ic, n_noVac, n_g, n_vacf, nb

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      ic = mtir_size(fi, mpdata%n_g, ikpt) - mpdata%n_g(ikpt)
      n_noVac = mtir_size(fi, mpdata%n_g, ikpt)
      n_g = mpdata%n_g(ikpt)
      nb = hybdat%nbasm(ikpt)
      n_vacf = nb - hybdat%n_mt - n_g
      IF (n_vacf <= 0 .OR. n_g <= 0) RETURN

      SELECT TYPE (coulomb)
      CLASS IS (t_mpimat)
#ifdef CPP_SCALAPACK
         CALL unpack_irvac_dist(fi, mpdata, hybdat, fmpi, coulomb, ikpt, ic, n_noVac, n_g, n_vacf)
#else
         CALL juDFT_error("distributed IR-VAC unpack needs ScaLAPACK", calledby="copy_irvac_to_sparse")
#endif
      CLASS IS (t_mat)
         IF (hybdat%coul(ikpt)%l_participate) &
            CALL unpack_irvac_local(fi, mpdata, hybdat, ikpt, ic, n_noVac, coulomb%data_c)
      END SELECT
   END SUBROUTINE copy_irvac_to_sparse

   !>Copy both IR x VAC rectangles into mtir; invs: vacuum index to (M_1 +- M_2)/sqrt(2).
   SUBROUTINE unpack_irvac_rect(fi, mpdata, hybdat, ikpt, ic, n_noVac, n_off, a_ir_vac, b_vac_ir)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)  :: hybdat
      INTEGER, INTENT(IN)            :: ikpt, ic, n_noVac, n_off
      COMPLEX, INTENT(IN)            :: a_ir_vac(:, :), b_vac_ir(:, :)

      INTEGER :: ivac, igm, ilen, nn, igpt, icol_d, irow_s, icol_s
      INTEGER :: islot, nslot, nsl, icol_d2, icol_s2, ngv
      COMPLEX :: a, b, a2, b2, va1, va2, vb1, vb2
      LOGICAL :: l_rot

      ngv = mpdata%n_g_vac(ikpt)
      l_rot = hybdat%coul(ikpt)%mtir%l_real .AND. NVAC_MPB == 2

      DO ivac = 1, NVAC_MPB
         IF (l_rot .AND. ivac == 2) CYCLE
         DO igm = 1, ngv
            ilen = mpdata%glen_ptr_vac(igm, ikpt)
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            nslot = MERGE(2, 1, mpdata%glen_vac(ilen) < 1e-12 .AND. nn >= 2)
            DO islot = 1, nslot
               nsl = MERGE(nn, nn - 1, islot == 1)
               icol_d = row_offset(mpdata, hybdat, ikpt, ivac) &
                        + basfn_offset(mpdata, ikpt, ivac, igm) + nsl - n_off
               IF (islot == 1) THEN
                  icol_s = n_noVac + (ivac - 1)*ngv + igm
               ELSE
                  icol_s = n_noVac + NVAC_MPB*ngv + ivac
               END IF
               IF (l_rot) THEN
                  icol_d2 = row_offset(mpdata, hybdat, ikpt, 2) &
                            + basfn_offset(mpdata, ikpt, 2, igm) + nsl - n_off
                  IF (islot == 1) THEN
                     icol_s2 = n_noVac + ngv + igm
                  ELSE
                     icol_s2 = n_noVac + NVAC_MPB*ngv + 2
                  END IF
               END IF
               DO igpt = 1, mpdata%n_g(ikpt)
                  irow_s = ic + igpt
                  a = a_ir_vac(igpt, icol_d)
                  b = b_vac_ir(icol_d, igpt)
                  IF (l_rot) THEN
                     a2 = a_ir_vac(igpt, icol_d2)
                     b2 = b_vac_ir(icol_d2, igpt)
                     va1 = (a + a2)/sqrt_2
                     va2 = ImagUnit*(a - a2)/sqrt_2
                     vb1 = (b + b2)/sqrt_2
                     vb2 = -ImagUnit*(b - b2)/sqrt_2
                     hybdat%coul(ikpt)%mtir%data_r(irow_s, icol_s) = REAL(va1)
                     hybdat%coul(ikpt)%mtir%data_r(icol_s, irow_s) = REAL(vb1)
                     hybdat%coul(ikpt)%mtir%data_r(irow_s, icol_s2) = REAL(va2)
                     hybdat%coul(ikpt)%mtir%data_r(icol_s2, irow_s) = REAL(vb2)
                  ELSE IF (hybdat%coul(ikpt)%mtir%l_real) THEN
                     hybdat%coul(ikpt)%mtir%data_r(irow_s, icol_s) = REAL(a)
                     hybdat%coul(ikpt)%mtir%data_r(icol_s, irow_s) = REAL(b)
                  ELSE
                     hybdat%coul(ikpt)%mtir%data_c(irow_s, icol_s) = a
                     hybdat%coul(ikpt)%mtir%data_c(icol_s, irow_s) = b
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE unpack_irvac_rect

   SUBROUTINE unpack_irvac_local(fi, mpdata, hybdat, ikpt, ic, n_noVac, src)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)  :: hybdat
      INTEGER, INTENT(IN)            :: ikpt, ic, n_noVac
      COMPLEX, INTENT(IN)            :: src(:, :)

      INTEGER :: n_mt, n_g, n_vacf, r0, c0

      n_mt = hybdat%n_mt
      n_g = mpdata%n_g(ikpt)
      n_vacf = hybdat%nbasm(ikpt) - n_mt - n_g
      r0 = n_mt
      c0 = n_mt + n_g

      CALL unpack_irvac_rect(fi, mpdata, hybdat, ikpt, ic, n_noVac, c0, &
                             src(r0 + 1:r0 + n_g, c0 + 1:c0 + n_vacf), &
                             src(c0 + 1:c0 + n_vacf, r0 + 1:r0 + n_g))
   END SUBROUTINE unpack_irvac_local

#ifdef CPP_SCALAPACK
   !>Distributed case: gather both rectangles on rank 0 with one BLACS grid.
   SUBROUTINE unpack_irvac_dist(fi, mpdata, hybdat, fmpi, coulomb, ikpt, ic, n_noVac, n_g, n_vacf)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(INOUT)  :: hybdat
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      TYPE(t_mpimat), INTENT(IN)     :: coulomb
      INTEGER, INTENT(IN)            :: ikpt, ic, n_noVac, n_g, n_vacf

      TYPE(t_mat) :: locA, locB
      INTEGER     :: ctxt, descA(9), descB(9), umap(1, 1), n_mt

      n_mt = hybdat%n_mt
      CALL locA%alloc(.FALSE., n_g, n_vacf)
      CALL locB%alloc(.FALSE., n_vacf, n_g)

      umap(1, 1) = 0
      CALL BLACS_GET(coulomb%blacsdata%blacs_desc(2), 10, ctxt)
      CALL BLACS_GRIDMAP(ctxt, umap, 1, 1, 1)
      descA = [1, ctxt, n_g, n_vacf, n_g, n_vacf, 0, 0, n_g]
      descB = [1, ctxt, n_vacf, n_g, n_vacf, n_g, 0, 0, n_vacf]

      CALL pzgemr2d(n_g, n_vacf, coulomb%data_c, n_mt + 1, n_mt + n_g + 1, &
                    coulomb%blacsdata%blacs_desc, locA%data_c, 1, 1, descA, &
                    coulomb%blacsdata%blacs_desc(2))
      CALL pzgemr2d(n_vacf, n_g, coulomb%data_c, n_mt + n_g + 1, n_mt + 1, &
                    coulomb%blacsdata%blacs_desc, locB%data_c, 1, 1, descB, &
                    coulomb%blacsdata%blacs_desc(2))

      IF (fmpi%n_rank == 0 .AND. hybdat%coul(ikpt)%l_participate) &
         CALL unpack_irvac_rect(fi, mpdata, hybdat, ikpt, ic, n_noVac, n_mt + n_g, &
                                locA%data_c, locB%data_c)

      CALL locA%free()
      CALL locB%free()
   END SUBROUTINE unpack_irvac_dist
#endif

   !>g = 0: slab integral of the linear potential (z1 - sigma z) Q + M.
   PURE COMPLEX FUNCTION irvac_g0_slab(a, z1, sigma, q0, q1)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a, z1, sigma, q0, q1

      COMPLEX :: j0c, j1c

      IF (ABS(a)*z1 < 1e-10) THEN
         j0c = CMPLX(2.0*z1, 0.0)
         j1c = cmplx_0
      ELSE
         j0c = CMPLX(2.0*SIN(a*z1)/a, 0.0)
         j1c = -2.0*ImagUnit*(SIN(a*z1)/(a*a) - z1*COS(a*z1)/a)
      END IF
      irvac_g0_slab = (z1*q0 + q1)*j0c - sigma*q0*j1c
   END FUNCTION irvac_g0_slab

   COMPLEX FUNCTION irvac_g0_mt(a, z1, sigma, za, q0, q1, dl, r)
      IMPLICIT NONE
      REAL, INTENT(IN) :: a, z1, sigma, za, q0, q1, dl, r

      COMPLEX :: gam
      REAL    :: c0, c1

      gam = CMPLX(0.0, -a)
      c0 = (z1 - sigma*za)*q0 + q1
      c1 = -sigma*q0
      irvac_g0_mt = c0*irvac_sphere(dl, gam, r) + c1*irvac_sphere_dgamma(dl, gam, r)
   END FUNCTION irvac_g0_mt

END MODULE m_irvac_2d
