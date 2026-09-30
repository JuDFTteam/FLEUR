!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>VAC-VAC and MT-VAC blocks of the film Coulomb matrix (IR-VAC: m_irvac_2d).
!>Per (length group, vacuum) only the last z-function carries the exponential moment; it
!>alone couples outside its own vacuum and is stored in mtir, the others in vac1/vac2.
!>The vacuum rows of cprod are <M_I|f>/sqrt(A), hence the factor A (sqrt(A) for cross blocks).
MODULE m_coulomb_vac_blocks
   USE m_vac_rows, ONLY: NVAC_MPB
   USE m_vac_rows, ONLY: vac_mtir_idx, vac_mtir_idx2
   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_types_coul, ONLY: t_coul
   USE m_mtir_size
   USE m_coulomb_vac, ONLY: vac_int_samevac, vac_int_samevac_g0, vac_exp_mom, vac_pref, vac_mom_g0, zquad
   USE m_coulomb_vac, ONLY: vac_elem_gram_cached

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: assemble_vac_blocks
   PUBLIC :: assemble_mtvac_blocks, vac_onsite

CONTAINS

   !>Same-vacuum element, symmetrised.
   REAL FUNCTION vac_onsite(f, h, nz, delz, g, area)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:), h(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz, g, area

      IF (g > 1e-12) THEN
         vac_onsite = 0.5*area*vac_pref(g)*(vac_int_samevac(f, h, nz, delz, g) &
                                            + vac_int_samevac(h, f, nz, delz, g))
      ELSE
         vac_onsite = 0.5*area*vac_pref(g)*(vac_int_samevac_g0(f, h, nz, delz) &
                                            + vac_int_samevac_g0(h, f, nz, delz))
      END IF
   END FUNCTION vac_onsite

   SUBROUTINE assemble_vac_blocks(fi, mpdata, coul, ikpt)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_coul), INTENT(INOUT)    :: coul
      INTEGER, INTENT(IN)            :: ikpt

      INTEGER :: ivac, jvac, ilen, jlen, i, j, nn, mm, nz, ig, jg, irow, icol
      INTEGER :: n_noVac
      REAL    :: g, gj, delz, val, mu_i, mu_j, z1, area, pref
      INTEGER :: ldz
      REAL, ALLOCATABLE :: agram(:, :), az(:, :), vblk(:, :)

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("VAC coulomb blocks")

      nz = mpdata%nmz_vac
      delz = fi%vacuum%delz
      area = fi%cell%omtil/fi%cell%amat(3, 3)
      z1 = fi%cell%z1
      n_noVac = mtir_size(fi, mpdata%n_g, ikpt)

      ! vac1/vac2: moment-free functions of one group, V = Z^T A Z
      ALLOCATE (agram(nz, nz))
      ldz = SIZE(mpdata%zbasfn_vac, 1)
      DO ivac = 1, NVAC_MPB
         DO ilen = 1, SIZE(mpdata%glen_vac)
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            g = mpdata%glen_vac(ilen)

            CALL vac_elem_gram_cached(nz, delz, MERGE(g, 0.0, g > 1e-12), agram)
            ALLOCATE (az(nz, nn), vblk(nn, nn))
            CALL dgemm('N', 'N', nz, nn, nz, 1.0, agram, nz, &
                       mpdata%zbasfn_vac(:, :, ilen, ivac), ldz, 0.0, az, nz)
            CALL dgemm('T', 'N', nn, nn, nz, 1.0, &
                       mpdata%zbasfn_vac(:, :, ilen, ivac), ldz, az, nz, 0.0, vblk, nn)
            pref = 0.5*area*vac_pref(g)

            DO i = 1, nn - 1
               DO j = 1, nn - 1
                  val = pref*(vblk(i, j) + vblk(j, i))
                  IF (fi%sym%invs) THEN
                     coul%vac1_r(i, j, ilen, ivac) = val
                  ELSE
                     coul%vac1_c(i, j, ilen, ivac) = CMPLX(val, 0.0)
                  END IF
               END DO

               val = pref*(vblk(i, nn) + vblk(nn, i))
               IF (fi%sym%invs) THEN
                  coul%vac2_r(i, ilen, ivac) = val
               ELSE
                  coul%vac2_c(i, ilen, ivac) = CMPLX(val, 0.0)
               END IF
            END DO
            DEALLOCATE (az, vblk)
         END DO
      END DO
      DEALLOCATE (agram)

      DO ivac = 1, NVAC_MPB
         DO ig = 1, mpdata%n_g_vac(ikpt)
            ilen = mpdata%glen_ptr_vac(ig, ikpt)
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            g = mpdata%glen_vac(ilen)
            irow = vac_mtir_idx(n_noVac, mpdata%n_g_vac(ikpt), ivac, ig)

            DO jvac = 1, NVAC_MPB
               DO jg = 1, mpdata%n_g_vac(ikpt)
                  IF (jg /= ig) CYCLE
                  jlen = mpdata%glen_ptr_vac(jg, ikpt)
                  mm = mpdata%num_zbasfn_vac(jlen, jvac)
                  IF (mm == 0) CYCLE
                  gj = mpdata%glen_vac(jlen)
                  icol = vac_mtir_idx(n_noVac, mpdata%n_g_vac(ikpt), jvac, jg)

                  ! carriers: same vacuum, or rank one across the slab (g = 0: applied in m_spmm_vac)
                  IF (jvac == ivac) THEN
                     val = vac_onsite(mpdata%zbasfn_vac(:, nn, ilen, ivac), &
                                      mpdata%zbasfn_vac(:, mm, jlen, jvac), nz, delz, g, area)
                  ELSE
                     IF (g > 1e-12) THEN
                        mu_i = vac_exp_mom(mpdata%zbasfn_vac(:, nn, ilen, ivac), nz, delz, g)
                        mu_j = vac_exp_mom(mpdata%zbasfn_vac(:, mm, jlen, jvac), nz, delz, gj)
                        val = area*vac_pref(g)*EXP(-2.0*g*z1)*mu_i*mu_j
                     ELSE
                        val = 0.0
                     END IF
                  END IF

                  IF (coul%mtir%l_real) THEN
                     coul%mtir%data_r(irow, icol) = val
                  ELSE
                     coul%mtir%data_c(irow, icol) = CMPLX(val, 0.0)
                  END IF
               END DO
            END DO
         END DO
      END DO

      ! invs: carriers in the real basis (M_1 +- M_2)/sqrt(2)
      IF (coul%mtir%l_real .AND. NVAC_MPB == 2) THEN
         DO ig = 1, mpdata%n_g_vac(ikpt)
            ilen = mpdata%glen_ptr_vac(ig, ikpt)
            IF (mpdata%num_zbasfn_vac(ilen, 1) == 0) CYCLE
            irow = vac_mtir_idx(n_noVac, mpdata%n_g_vac(ikpt), 1, ig)
            icol = vac_mtir_idx(n_noVac, mpdata%n_g_vac(ikpt), 2, ig)
            val = coul%mtir%data_r(irow, irow)
            g = coul%mtir%data_r(irow, icol)
            coul%mtir%data_r(irow, irow) = val + g
            coul%mtir%data_r(icol, icol) = val - g
            coul%mtir%data_r(irow, icol) = 0.0
            coul%mtir%data_r(icol, irow) = 0.0
         END DO
      END IF

      CALL timestop("VAC coulomb blocks")
   END SUBROUTINE assemble_vac_blocks

   !>MT-VAC: only MT and vacuum moment carriers couple; at g = 0 a second vacuum carrier holds
   !>the first z-moment.  The atom phase uses G|| only, since V(q) is reused for all q in the star.
   SUBROUTINE assemble_mtvac_blocks(fi, mpdata, coul, moment, ikpt)
      USE m_mtvac_2d, ONLY: mtvac_shape
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_coul), INTENT(INOUT)    :: coul
      REAL, INTENT(IN)               :: moment(:, 0:, :)
      INTEGER, INTENT(IN)            :: ikpt

      INTEGER :: icol, icolp, ivac, igm, ilen, nn, nz, i
      LOGICAL :: l_rot
      INTEGER :: n_noVac, ngv, islot, nslot, n_mt_red, ncolv, j
      COMPLEX, ALLOCATABLE :: mtv(:, :)
      REAL    :: kv3(3), kv(2), kc3(3), gv(2), kk, sigma, z1, area, cross(3), mu, delz
      REAL    :: q_i, zm_i
      LOGICAL :: l_g0

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("MT-VAC coulomb blocks")
      cross(1) = fi%cell%amat(2,1)*fi%cell%amat(3,2) - fi%cell%amat(3,1)*fi%cell%amat(2,2)
      cross(2) = fi%cell%amat(3,1)*fi%cell%amat(1,2) - fi%cell%amat(1,1)*fi%cell%amat(3,2)
      cross(3) = fi%cell%amat(1,1)*fi%cell%amat(2,2) - fi%cell%amat(2,1)*fi%cell%amat(1,2)
      area = NORM2(cross)
      z1 = fi%cell%z1
      nz = mpdata%nmz_vac
      delz = fi%vacuum%delz
      n_noVac = mtir_size(fi, mpdata%n_g, ikpt)
      ngv = mpdata%n_g_vac(ikpt)

      l_rot = coul%mtir%l_real .AND. NVAC_MPB == 2
      n_mt_red = n_noVac - mpdata%n_g(ikpt)
      kc3 = MATMUL(fi%kpts%bk(:, ikpt), fi%cell%bmat)
      ncolv = NVAC_MPB*ngv + NVAC_MPB
      IF (l_rot) ALLOCATE (mtv(n_mt_red, ncolv), source=cmplx_0)
      DO ivac = 1, NVAC_MPB
         sigma = MERGE(1.0, -1.0, ivac == 1)
         DO igm = 1, ngv
            ilen = mpdata%glen_ptr_vac(igm, ikpt)
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            kk = mpdata%glen_vac(ilen)
            l_g0 = (kk < 1e-12)
            kv3 = MATMUL(fi%kpts%bk(:, ikpt) &
                         + [REAL(mpdata%g_vac(1, mpdata%gptm_ptr_vac(igm, ikpt))), &
                            REAL(mpdata%g_vac(2, mpdata%gptm_ptr_vac(igm, ikpt))), 0.0], &
                         fi%cell%bmat)
            kv = kv3(1:2)
            gv = kv - kc3(1:2)
            nslot = MERGE(2, 1, l_g0 .AND. nn >= 2)
            DO islot = 1, nslot
               IF (islot == 1) THEN
                  icol = vac_mtir_idx(n_noVac, ngv, ivac, igm)
                  icolp = vac_mtir_idx(n_noVac, ngv, 2, igm)
                  IF (l_g0) THEN
                     CALL vac_mom_g0(mpdata%zbasfn_vac(:, nn, ilen, ivac), nz, delz, q_i, zm_i)
                     mu = 0.0
                  ELSE
                     mu = vac_exp_mom(mpdata%zbasfn_vac(:, nn, ilen, ivac), nz, delz, kk)
                  END IF
               ELSE
                  icol = vac_mtir_idx2(n_noVac, ngv, NVAC_MPB, ivac)
                  icolp = vac_mtir_idx2(n_noVac, ngv, NVAC_MPB, 2)
                  CALL vac_mom_g0(mpdata%zbasfn_vac(:, nn - 1, ilen, ivac), nz, delz, q_i, zm_i)
                  mu = 0.0
               END IF
               IF (l_rot) THEN
                  CALL mtvac_write_col(fi, mpdata, coul, moment, kv, gv, sigma, z1, area, kk, &
                                       l_g0, q_i, zm_i, mu, icol, buf=mtv, jcol=icol - n_noVac)
               ELSE
                  CALL mtvac_write_col(fi, mpdata, coul, moment, kv, gv, sigma, z1, area, kk, &
                                       l_g0, q_i, zm_i, mu, icol)
               END IF
            END DO
         END DO
      END DO
      IF (l_rot) THEN
         block
            USE m_trafo, ONLY: symmetrize
            INTEGER :: ones(0:MAXVAL(fi%hybinp%lcutm1), fi%atoms%ntype)
            INTEGER :: jc1, jc2
            COMPLEX :: w1, w2
            ones = 1
            ! invs: MT index to real harmonics, vacuum index to (M_1 +- M_2)/sqrt(2)
            CALL symmetrize(mtv, n_mt_red, ncolv, 1, fi%atoms, fi%hybinp%lcutm1, &
                            MAXVAL(fi%hybinp%lcutm1), ones, fi%sym)
            DO igm = 1, ngv
               jc1 = igm; jc2 = ngv + igm
               DO i = 1, n_mt_red
                  w1 = (mtv(i, jc1) + mtv(i, jc2))/sqrt_2
                  w2 = ImagUnit*(mtv(i, jc1) - mtv(i, jc2))/sqrt_2
                  mtv(i, jc1) = w1; mtv(i, jc2) = w2
               END DO
            END DO
            jc1 = NVAC_MPB*ngv + 1; jc2 = NVAC_MPB*ngv + 2
            DO i = 1, n_mt_red
               w1 = (mtv(i, jc1) + mtv(i, jc2))/sqrt_2
               w2 = ImagUnit*(mtv(i, jc1) - mtv(i, jc2))/sqrt_2
               mtv(i, jc1) = w1; mtv(i, jc2) = w2
            END DO
            DO j = 1, ncolv
               DO i = 1, n_mt_red
                  coul%mtir%data_r(n_noVac + j, i) = REAL(mtv(i, j))
                  coul%mtir%data_r(i, n_noVac + j) = REAL(mtv(i, j))
               END DO
            END DO
         end block
         DEALLOCATE (mtv)
      END IF

      CALL timestop("MT-VAC coulomb blocks")
   END SUBROUTINE assemble_mtvac_blocks

   SUBROUTINE mtvac_write_col(fi, mpdata, coul, moment, kv, gv, sigma, z1, area, kk, &
                              l_g0, q_i, zm_i, mu, icol, buf, jcol)
      USE m_mtvac_2d, ONLY: mtvac_shape
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_coul), INTENT(INOUT)    :: coul
      REAL, INTENT(IN)               :: moment(:, 0:, :)
      REAL, INTENT(IN)               :: kv(2), gv(2), sigma, z1, area, kk, q_i, zm_i, mu
      LOGICAL, INTENT(IN)            :: l_g0
      INTEGER, INTENT(IN)            :: icol
      COMPLEX, INTENT(INOUT), OPTIONAL :: buf(:, :)
      INTEGER, INTENT(IN), OPTIONAL    :: jcol

      INTEGER :: ia, itype, l, m, lm, mtoff, irow, n
      REAL    :: za, mom, sg
      COMPLEX :: x
      LOGICAL :: l_buf

      l_buf = PRESENT(buf)
      sg = sigma

      mtoff = 0
      DO ia = 1, fi%atoms%nat
         itype = fi%atoms%itype(ia)
         za = fi%atoms%pos(3, ia)
         lm = 0
         DO l = 0, fi%hybinp%lcutm1(itype)
            n = mpdata%num_radbasfn(l, itype)
            mom = moment(n, l, itype)
            DO m = -l, l
               lm = lm + 1
               irow = mtoff + lm
               IF (l_g0) THEN
                  IF (l == 0 .AND. m == 0) THEN
                     ! g = 0: the vacuum potential is linear in z, only (0,0) and (1,0) survive
                     x = -(tpi_const/SQRT(area))*((z1 - sg*za)*q_i + zm_i) &
                         *sfp_const*mom
                  ELSE IF (l == 1 .AND. m == 0) THEN
                     x = (tpi_const/SQRT(area))*sg*q_i*SQRT(fpi_const/3.0)*mom
                  ELSE
                     x = cmplx_0
                  END IF
               ELSE
                  x = mtvac_shape(l, m, kv, sg, area)*mom &
                      *EXP(-kk*(z1 - sg*za))*mu &
                      *EXP(-ImagUnit*(gv(1)*fi%atoms%pos(1, ia) + gv(2)*fi%atoms%pos(2, ia)))
               END IF
               x = x*SQRT(area)

               IF (l_buf) THEN
                  ! <MT|v|vac>, as expected by the rotations
                  buf(irow, jcol) = CONJG(x)
               ELSE IF (coul%mtir%l_real) THEN
                  coul%mtir%data_r(icol, irow) = REAL(x)
                  coul%mtir%data_r(irow, icol) = REAL(x)
               ELSE
                  coul%mtir%data_c(icol, irow) = x
                  coul%mtir%data_c(irow, icol) = CONJG(x)
               END IF
            END DO
         END DO
         mtoff = mtoff + (fi%hybinp%lcutm1(itype) + 1)**2
      END DO
   END SUBROUTINE mtvac_write_col

END MODULE m_coulomb_vac_blocks
