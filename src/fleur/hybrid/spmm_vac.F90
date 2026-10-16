!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Application of the vacuum part of the Coulomb matrix in spmm_inv/spmm_noinv.  The callers'
!>GEMM covers the MT+IR corner of mtir only; the VAC carriers are interleaved with the
!>moment-free functions in the vector layout, so they are gathered and applied here.
MODULE m_spmm_vac
   USE m_vac_rows, ONLY: NVAC_MPB
   USE m_juDFT
   USE m_constants
   USE m_mtir_size
   USE m_coulomb_vac, ONLY: vac_mom_g0
   USE m_types_coul, ONLY: t_coul
   USE m_vac_rows, ONLY: row_offset, basfn_offset, vac_g0_moments
   USE m_vac_rows, ONLY: vac_mtir_idx, vac_mtir_idx2
   USE m_types_fleurinput
   USE m_types_hybdat
   USE m_types_mpdata

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: spmm_vac_r, spmm_vac_c, vac_spmm_layout
   PUBLIC :: apply_mtir_vac_r, apply_mtir_vac_c

CONTAINS

   SUBROUTINE vac_spmm_layout(fi, mpdata, hybdat, ikpt, rbase, nfun, ilenv, nvac, ngv)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      INTEGER, INTENT(IN)            :: ikpt
      INTEGER, ALLOCATABLE, INTENT(OUT) :: rbase(:, :), nfun(:, :), ilenv(:)
      INTEGER, INTENT(OUT)           :: nvac, ngv

      INTEGER :: ivac, igm, ilen

      nvac = NVAC_MPB
      ngv = mpdata%n_g_vac(ikpt)
      ALLOCATE (rbase(nvac, ngv), nfun(nvac, ngv), ilenv(ngv))

      DO igm = 1, ngv
         ilen = mpdata%glen_ptr_vac(igm, ikpt)
         ilenv(igm) = ilen
         DO ivac = 1, nvac
            nfun(ivac, igm) = mpdata%num_zbasfn_vac(ilen, ivac)
            rbase(ivac, igm) = row_offset(mpdata, hybdat, ikpt, ivac) &
                               + basfn_offset(mpdata, ikpt, ivac, igm)
         END DO
      END DO
   END SUBROUTINE vac_spmm_layout

   !>g = 0 coupling across the slab, -A 2 pi (2 z1 Q_v Q_w + M_v Q_w + Q_v M_w): rank two, applied
   !>from the charges Q and first moments M of the input vector.
   SUBROUTINE apply_vac_g0_cross(fi, mpdata, hybdat, ikpt, mat_in, mat_out)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      INTEGER, INTENT(IN)            :: ikpt
      COMPLEX, INTENT(IN)            :: mat_in(:, :)
      COMPLEX, INTENT(INOUT)         :: mat_out(:, :)

      INTEGER :: iv, jv, i, nn, ilen, bi, igm0, n_vec, ivec, nz
      REAL    :: area, z1, delz, q_i, zm_i
      COMPLEX :: qv(NVAC_MPB), mv(NVAC_MPB), sq, sm
      LOGICAL :: l_rot

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN
      area = fi%cell%omtil/fi%cell%amat(3, 3); z1 = fi%cell%z1
      nz = mpdata%nmz_vac; delz = fi%vacuum%delz
      n_vec = SIZE(mat_in, 2)
      l_rot = fi%sym%invs .AND. NVAC_MPB == 2

      DO ivec = 1, n_vec
         CALL vac_g0_moments(mpdata, hybdat, ikpt, NVAC_MPB, nz, delz, &
                             mat_in(:, ivec), qv, mv, igm0)
         IF (igm0 == 0) RETURN
         ilen = mpdata%glen_ptr_vac(igm0, ikpt)
         DO iv = 1, NVAC_MPB
            nn = mpdata%num_zbasfn_vac(ilen, iv)
            IF (nn == 0) CYCLE
            IF (l_rot) THEN
               sq = MERGE(1.0, -1.0, iv == 1)*qv(iv)
               sm = MERGE(1.0, -1.0, iv == 1)*mv(iv)
            ELSE
               sq = cmplx_0; sm = cmplx_0
               DO jv = 1, NVAC_MPB
                  IF (jv == iv) CYCLE
                  sq = sq + qv(jv); sm = sm + mv(jv)
               END DO
            END IF
            bi = row_offset(mpdata, hybdat, ikpt, iv) + basfn_offset(mpdata, ikpt, iv, igm0)
            DO i = 1, nn
               CALL vac_mom_g0(mpdata%zbasfn_vac(:, i, ilen, iv), nz, delz, q_i, zm_i)
               mat_out(bi + i, ivec) = mat_out(bi + i, ivec) &
                  - area*tpi_const*((2.0*z1*q_i + zm_i)*sq + q_i*sm)
            END DO
         END DO
      END DO
   END SUBROUTINE apply_vac_g0_cross

   !>(mtir column, vector row) of every VAC carrier; the g = 0 second carrier is function nn-1.
   SUBROUTINE vac_mtir_pairs(fi, mpdata, hybdat, ikpt, vcol, vrow, npair)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN)    :: fi
      TYPE(t_mpdata), INTENT(IN)        :: mpdata
      TYPE(t_hybdat), INTENT(IN)        :: hybdat
      INTEGER, INTENT(IN)               :: ikpt
      INTEGER, ALLOCATABLE, INTENT(OUT) :: vcol(:), vrow(:)
      INTEGER, INTENT(OUT)              :: npair

      INTEGER :: nvac, ngv, igm, ivac, nn, islot, nslot, n_noVac
      INTEGER, ALLOCATABLE :: rbase(:, :), nfun(:, :), ilenv(:)

      npair = 0
      CALL vac_spmm_layout(fi, mpdata, hybdat, ikpt, rbase, nfun, ilenv, nvac, ngv)
      n_noVac = mtir_size(fi, mpdata%n_g, ikpt)

      ALLOCATE (vcol(2*nvac*ngv), vrow(2*nvac*ngv))

      DO igm = 1, ngv
         DO ivac = 1, nvac
            nn = nfun(ivac, igm)
            IF (nn == 0) CYCLE
            nslot = MERGE(2, 1, mpdata%glen_vac(ilenv(igm)) < 1e-12 &
                          .AND. nn >= 2)
            DO islot = 1, nslot
               npair = npair + 1
               IF (islot == 1) THEN
                  vcol(npair) = vac_mtir_idx(n_noVac, ngv, ivac, igm)
                  vrow(npair) = rbase(ivac, igm) + nn
               ELSE
                  vcol(npair) = vac_mtir_idx2(n_noVac, ngv, nvac, ivac)
                  vrow(npair) = rbase(ivac, igm) + nn - 1
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE vac_mtir_pairs

   !>MT-VAC and IR-VAC couplings; conjugated together with the rest of mtir when required.
   SUBROUTINE apply_mtir_vac_c(fi, mpdata, hybdat, coul, ikpt, ibasm, indx1, conjg_mtir, mat_in, mat_out)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_coul), INTENT(IN)       :: coul
      INTEGER, INTENT(IN)            :: ikpt, ibasm, indx1
      LOGICAL, INTENT(IN)            :: conjg_mtir
      COMPLEX, INTENT(IN)            :: mat_in(:, :)
      COMPLEX, INTENT(INOUT)         :: mat_out(:, :)

      INTEGER :: j, iv, n_vec, ip, npair, vcol, vrow
      INTEGER, ALLOCATABLE :: vcols(:), vrows(:)

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("spmm MT/IR-VAC")
      CALL vac_mtir_pairs(fi, mpdata, hybdat, ikpt, vcols, vrows, npair)
      n_vec = SIZE(mat_in, 2)

      DO iv = 1, n_vec
         DO ip = 1, npair
            vcol = vcols(ip); vrow = vrows(ip)
            IF (conjg_mtir) THEN
               DO j = 1, indx1
                  mat_out(ibasm + j, iv) = mat_out(ibasm + j, iv) &
                                           + CONJG(coul%mtir%data_c(j, vcol))*mat_in(vrow, iv)
                  mat_out(vrow, iv) = mat_out(vrow, iv) &
                                      + CONJG(coul%mtir%data_c(vcol, j))*mat_in(ibasm + j, iv)
               END DO
            ELSE
               DO j = 1, indx1
                  mat_out(ibasm + j, iv) = mat_out(ibasm + j, iv) &
                                           + coul%mtir%data_c(j, vcol)*mat_in(vrow, iv)
                  mat_out(vrow, iv) = mat_out(vrow, iv) &
                                      + coul%mtir%data_c(vcol, j)*mat_in(ibasm + j, iv)
               END DO
            END IF
         END DO
      END DO
      CALL timestop("spmm MT/IR-VAC")
   END SUBROUTINE apply_mtir_vac_c

   SUBROUTINE apply_mtir_vac_r(fi, mpdata, hybdat, coul, ikpt, ibasm, indx1, mat_in, mat_out)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_coul), INTENT(IN)       :: coul
      INTEGER, INTENT(IN)            :: ikpt, ibasm, indx1
      REAL, INTENT(IN)               :: mat_in(:, :)
      REAL, INTENT(INOUT)            :: mat_out(:, :)

      INTEGER :: j, iv, n_vec, ip, npair, vcol, vrow
      INTEGER, ALLOCATABLE :: vcols(:), vrows(:)

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("spmm MT/IR-VAC")
      CALL vac_mtir_pairs(fi, mpdata, hybdat, ikpt, vcols, vrows, npair)
      n_vec = SIZE(mat_in, 2)

      DO iv = 1, n_vec
         DO ip = 1, npair
            vcol = vcols(ip); vrow = vrows(ip)
            DO j = 1, indx1
               mat_out(ibasm + j, iv) = mat_out(ibasm + j, iv) &
                                        + coul%mtir%data_r(j, vcol)*mat_in(vrow, iv)
               mat_out(vrow, iv) = mat_out(vrow, iv) &
                                   + coul%mtir%data_r(vcol, j)*mat_in(ibasm + j, iv)
            END DO
         END DO
      END DO
      CALL timestop("spmm MT/IR-VAC")
   END SUBROUTINE apply_mtir_vac_r

   !>VAC-VAC part: vac1/vac2 within each group, carriers via mtir.  Adds to mat_out.
   SUBROUTINE spmm_vac_r(fi, mpdata, hybdat, coul, ikpt, mat_in, mat_out)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_coul), INTENT(IN)       :: coul
      INTEGER, INTENT(IN)            :: ikpt
      REAL, INTENT(IN)               :: mat_in(:, :)
      REAL, INTENT(INOUT)            :: mat_out(:, :)

      INTEGER :: nvac, ngv, igm, ivac, jvac, ilen, nn, mm, i, j, iv, n_vec, b, n_noVac
      INTEGER, ALLOCATABLE :: rbase(:, :), nfun(:, :), ilenv(:)
      REAL, ALLOCATABLE    :: v1(:, :, :, :), v2(:, :, :), vmc(:, :, :)
      REAL :: s

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("spmm_vac_r")
      CALL vac_spmm_layout(fi, mpdata, hybdat, ikpt, rbase, nfun, ilenv, nvac, ngv)
      n_vec = SIZE(mat_in, 2)
      n_noVac = mtir_size(fi, mpdata%n_g, ikpt)

      ALLOCATE (vmc(nvac, nvac, ngv), source=0.0)
      DO igm = 1, ngv
         DO ivac = 1, nvac
            DO jvac = 1, nvac
               vmc(ivac, jvac, igm) = coul%mtir%data_r( &
                                      vac_mtir_idx(n_noVac, ngv, ivac, igm), &
                                      vac_mtir_idx(n_noVac, ngv, jvac, igm))
            END DO
         END DO
      END DO
      v1 = coul%vac1_r
      v2 = coul%vac2_r

      !$acc data copyin(v1, v2, vmc, rbase, nfun, ilenv) copyin(mat_in) copy(mat_out)
      !$acc kernels
      !$acc loop independent private(igm, ivac, jvac, ilen, nn, mm, i, j, b, s)
      DO iv = 1, n_vec
         !$acc loop seq
         DO igm = 1, ngv
            !$acc loop seq
            DO ivac = 1, nvac
               nn = nfun(ivac, igm)
               IF (nn == 0) CYCLE
               ilen = ilenv(igm)
               b = rbase(ivac, igm)

               DO i = 1, nn - 1
                  s = 0.0
                  DO j = 1, nn - 1
                     s = s + v1(i, j, ilen, ivac)*mat_in(b + j, iv)
                  END DO
                  s = s + v2(i, ilen, ivac)*mat_in(b + nn, iv)
                  mat_out(b + i, iv) = mat_out(b + i, iv) + s
               END DO

               s = 0.0
               DO i = 1, nn - 1
                  s = s + v2(i, ilen, ivac)*mat_in(b + i, iv)
               END DO
               DO jvac = 1, nvac
                  mm = nfun(jvac, igm)
                  IF (mm == 0) CYCLE
                  s = s + vmc(ivac, jvac, igm)*mat_in(rbase(jvac, igm) + mm, iv)
               END DO
               mat_out(b + nn, iv) = mat_out(b + nn, iv) + s
            END DO
         END DO
      END DO
      !$acc end kernels
      !$acc end data

      CALL timestop("spmm_vac_r")
      block
         COMPLEX, ALLOCATABLE :: cin(:, :), cout(:, :)
         ALLOCATE (cin(SIZE(mat_in, 1), SIZE(mat_in, 2)), source=CMPLX(mat_in, 0.0))
         ALLOCATE (cout(SIZE(mat_out, 1), SIZE(mat_out, 2)), source=cmplx_0)
         CALL apply_vac_g0_cross(fi, mpdata, hybdat, ikpt, cin, cout)
         mat_out = mat_out + REAL(cout)
      end block
   END SUBROUTINE spmm_vac_r

   SUBROUTINE spmm_vac_c(fi, mpdata, hybdat, coul, ikpt, conjg_mtir, mat_in, mat_out)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_coul), INTENT(IN)       :: coul
      INTEGER, INTENT(IN)            :: ikpt
      LOGICAL, INTENT(IN)            :: conjg_mtir
      COMPLEX, INTENT(IN)            :: mat_in(:, :)
      COMPLEX, INTENT(INOUT)         :: mat_out(:, :)

      INTEGER :: nvac, ngv, igm, ivac, jvac, ilen, nn, mm, i, j, iv, n_vec, b, n_noVac
      INTEGER, ALLOCATABLE :: rbase(:, :), nfun(:, :), ilenv(:)
      COMPLEX, ALLOCATABLE :: v1(:, :, :, :), v2(:, :, :), vmc(:, :, :)
      COMPLEX :: s

      IF (.NOT. fi%input%film) RETURN
      IF (.NOT. ALLOCATED(mpdata%num_zbasfn_vac)) RETURN

      CALL timestart("spmm_vac_c")
      CALL vac_spmm_layout(fi, mpdata, hybdat, ikpt, rbase, nfun, ilenv, nvac, ngv)
      n_vec = SIZE(mat_in, 2)
      n_noVac = mtir_size(fi, mpdata%n_g, ikpt)

      ALLOCATE (vmc(nvac, nvac, ngv), source=cmplx_0)
      DO igm = 1, ngv
         DO ivac = 1, nvac
            DO jvac = 1, nvac
               vmc(ivac, jvac, igm) = coul%mtir%data_c( &
                                      vac_mtir_idx(n_noVac, ngv, ivac, igm), &
                                      vac_mtir_idx(n_noVac, ngv, jvac, igm))
            END DO
         END DO
      END DO
      v1 = coul%vac1_c
      v2 = coul%vac2_c

      IF (conjg_mtir) THEN
         v1 = CONJG(v1); v2 = CONJG(v2); vmc = CONJG(vmc)
      END IF

      !$acc data copyin(v1, v2, vmc, rbase, nfun, ilenv) copyin(mat_in) copy(mat_out)
      !$acc kernels
      !$acc loop independent private(igm, ivac, jvac, ilen, nn, mm, i, j, b, s)
      DO iv = 1, n_vec
         !$acc loop seq
         DO igm = 1, ngv
            !$acc loop seq
            DO ivac = 1, nvac
               nn = nfun(ivac, igm)
               IF (nn == 0) CYCLE
               ilen = ilenv(igm)
               b = rbase(ivac, igm)

               DO i = 1, nn - 1
                  s = cmplx_0
                  DO j = 1, nn - 1
                     s = s + v1(i, j, ilen, ivac)*mat_in(b + j, iv)
                  END DO
                  s = s + v2(i, ilen, ivac)*mat_in(b + nn, iv)
                  mat_out(b + i, iv) = mat_out(b + i, iv) + s
               END DO

               s = cmplx_0
               DO i = 1, nn - 1
                  s = s + v2(i, ilen, ivac)*mat_in(b + i, iv)
               END DO
               DO jvac = 1, nvac
                  mm = nfun(jvac, igm)
                  IF (mm == 0) CYCLE
                  s = s + vmc(ivac, jvac, igm)*mat_in(rbase(jvac, igm) + mm, iv)
               END DO
               mat_out(b + nn, iv) = mat_out(b + nn, iv) + s
            END DO
         END DO
      END DO
      !$acc end kernels
      !$acc end data

      CALL timestop("spmm_vac_c")
      CALL apply_vac_g0_cross(fi, mpdata, hybdat, ikpt, mat_in, mat_out)
   END SUBROUTINE spmm_vac_c

END MODULE m_spmm_vac
