!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Vacuum rows of cprod: sum over in-plane pairs G||_j' = G||_j + Gm of
!>conjg(C_k) int U_i f_j f_j' dz C_kq, with f = u, udot.  The rows are <M_I|f>/sqrt(A).
MODULE m_wavefproducts_vac
   USE m_vac_rows, ONLY: NVAC_MPB
   USE m_juDFT
   USE m_constants
   USE m_vac_rows, ONLY: row_offset, basfn_offset
   USE m_vac_abcof
   USE m_coulomb_vac, ONLY: vac_zint
   USE m_types_fleurinput
   USE m_types_hybdat
   USE m_types_lapw
   USE m_types_mat
   USE m_types_mpdata
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wavefproducts_vac

CONTAINS

   SUBROUTINE wavefproducts_vac(fi, ik, iq, ikqpt, g_t, jsp, bandoi, bandof, &
                                mpdata, hybdat, lapw, lapw_kq, z_k, z_kq, cprod)
      IMPLICIT NONE
      TYPE(t_fleurinput), INTENT(IN) :: fi
      INTEGER, INTENT(IN)            :: ik, iq, ikqpt, g_t(3), jsp, bandoi, bandof
      TYPE(t_mpdata), INTENT(IN)     :: mpdata
      TYPE(t_hybdat), INTENT(IN)     :: hybdat
      TYPE(t_lapw), INTENT(IN)       :: lapw, lapw_kq
      TYPE(t_mat), INTENT(IN)        :: z_k, z_kq
      COMPLEX, INTENT(INOUT)         :: cprod(:, :)

      INTEGER :: nv2d, nz, ivac, nbnd_k, psize, row0, row
      INTEGER :: nv2_k, nv2_kq, igm, gm(2), ilen, ni, i, j, jp, c, cp, n, np
      INTEGER :: ok

      INTEGER, ALLOCATABLE :: kvac1_k(:), kvac2_k(:), map2_k(:)
      INTEGER, ALLOCATABLE :: kvac1_q(:), kvac2_q(:), map2_q(:)
      INTEGER, ALLOCATABLE :: jpart(:)
      REAL, ALLOCATABLE    :: u_k(:, :), ue_k(:, :), u_q(:, :), ue_q(:, :)
      REAL, ALLOCATABLE    :: t_k(:), dt_k(:), te_k(:), dte_k(:), tei_k(:)
      REAL, ALLOCATABLE    :: t_q(:), dt_q(:), te_q(:), dte_q(:), tei_q(:)
      COMPLEX, ALLOCATABLE :: ac_k(:, :), bc_k(:, :), ac_q(:, :), bc_q(:, :)
      COMPLEX, ALLOCATABLE :: cf_k(:, :)
      COMPLEX, ALLOCATABLE :: cf_q(:, :, :)
      REAL, ALLOCATABLE    :: zint(:, :, :)
      COMPLEX, ALLOCATABLE :: tmat(:, :)
      COMPLEX, ALLOCATABLE :: blk(:, :)

      IF (.NOT. fi%input%film) RETURN

      CALL timestart("wavefproducts_vac")

      nv2d = lapw%dim_nv2d()
      nz = mpdata%nmz_vac
      nbnd_k = hybdat%nbands(ik, jsp)
      psize = bandof - bandoi + 1

      ALLOCATE (kvac1_k(nv2d), kvac2_k(nv2d), map2_k(lapw%dim_nvd()), source=0)
      ALLOCATE (kvac1_q(nv2d), kvac2_q(nv2d), map2_q(lapw%dim_nvd()), source=0)
      CALL vac_map2(lapw, jsp, nv2d, kvac1_k, kvac2_k, map2_k, nv2_k)
      CALL vac_map2(lapw_kq, jsp, nv2d, kvac1_q, kvac2_q, map2_q, nv2_kq)

      ALLOCATE (u_k(fi%vacuum%nmz, nv2d), ue_k(fi%vacuum%nmz, nv2d), source=0.0)
      ALLOCATE (u_q(fi%vacuum%nmz, nv2d), ue_q(fi%vacuum%nmz, nv2d), source=0.0)
      ALLOCATE (t_k(nv2d), dt_k(nv2d), te_k(nv2d), dte_k(nv2d), tei_k(nv2d), source=0.0)
      ALLOCATE (t_q(nv2d), dt_q(nv2d), te_q(nv2d), dte_q(nv2d), tei_q(nv2d), source=0.0)
      ALLOCATE (ac_k(nv2d, nbnd_k), bc_k(nv2d, nbnd_k), source=cmplx_0)
      ALLOCATE (ac_q(nv2d, psize), bc_q(nv2d, psize), source=cmplx_0)
      ALLOCATE (cf_k(2*nv2d, nbnd_k), source=cmplx_0)
      ALLOCATE (cf_q(nv2d, psize, 2), source=cmplx_0)
      ALLOCATE (jpart(nv2d), source=0)
      ALLOCATE (zint(nv2d, 2, 2), source=0.0)
      ALLOCATE (tmat(2*nv2d, psize), blk(nbnd_k, psize), stat=ok)
      IF (ok /= 0) CALL juDFT_error("wavefproducts_vac: allocation failed")

      DO ivac = 1, NVAC_MPB
         CALL vac_uz(fi%vacuum, fi%cell, hybdat%evac_vac(ivac, jsp), &
                     hybdat%vz_vac(:, ivac, jsp), lapw%bkpt(1:2), &
                     kvac1_k, kvac2_k, nv2_k, u_k, ue_k, t_k, dt_k, te_k, dte_k, tei_k)
         CALL vac_uz(fi%vacuum, fi%cell, hybdat%evac_vac(ivac, jsp), &
                     hybdat%vz_vac(:, ivac, jsp), lapw_kq%bkpt(1:2), &
                     kvac1_q, kvac2_q, nv2_kq, u_q, ue_q, t_q, dt_q, te_q, dte_q, tei_q)

         ac_k = cmplx_0; bc_k = cmplx_0
         ac_q = cmplx_0; bc_q = cmplx_0
         CALL vac_abcof(fi%cell, lapw, jsp, ivac, nv2d, nbnd_k, 0, 1.0, map2_k, &
                        t_k, dt_k, te_k, dte_k, z_k, ac_k, bc_k)
         ! z_kq holds only the bands bandoi..bandof of this package, as columns 1..psize
         CALL vac_abcof(fi%cell, lapw_kq, jsp, ivac, nv2d, psize, 0, 1.0, map2_q, &
                        t_q, dt_q, te_q, dte_q, z_kq, ac_q, bc_q)

         cf_k = cmplx_0
         cf_k(1:nv2_k, :) = ac_k(1:nv2_k, :)
         cf_k(nv2_k + 1:2*nv2_k, :) = bc_k(1:nv2_k, :)
         cf_q(:, :, 1) = ac_q; cf_q(:, :, 2) = bc_q

         row0 = row_offset(mpdata, hybdat, iq, ivac)
         DO igm = 1, mpdata%n_g_vac(iq)
            gm = mpdata%g_vac(:, mpdata%gptm_ptr_vac(igm, iq))
            ilen = mpdata%glen_ptr_vac(igm, iq)
            ni = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (ni == 0) CYCLE

            CALL partner_map(kvac1_k, kvac2_k, nv2_k, kvac1_q, kvac2_q, nv2_kq, &
                             gm - g_t(1:2), jpart)

            DO i = 1, ni
               zint = 0.0
               DO j = 1, nv2_k
                  jp = jpart(j)
                  IF (jp == 0) CYCLE
                  zint(j, 1, 1) = z_dot3(mpdata%zbasfn_vac(:, i, ilen, ivac), &
                                        u_k(:, j), u_q(:, jp), nz, fi%vacuum%delz)
                  zint(j, 1, 2) = z_dot3(mpdata%zbasfn_vac(:, i, ilen, ivac), &
                                        u_k(:, j), ue_q(:, jp), nz, fi%vacuum%delz)
                  zint(j, 2, 1) = z_dot3(mpdata%zbasfn_vac(:, i, ilen, ivac), &
                                        ue_k(:, j), u_q(:, jp), nz, fi%vacuum%delz)
                  zint(j, 2, 2) = z_dot3(mpdata%zbasfn_vac(:, i, ilen, ivac), &
                                        ue_k(:, j), ue_q(:, jp), nz, fi%vacuum%delz)
               END DO

               tmat = cmplx_0
               DO c = 1, 2
                  DO j = 1, nv2_k
                     jp = jpart(j)
                     IF (jp == 0) CYCLE
                     DO cp = 1, 2
                        IF (ABS(zint(j, c, cp)) < 1e-30) CYCLE
                        tmat(j + (c - 1)*nv2_k, 1:psize) = tmat(j + (c - 1)*nv2_k, 1:psize) &
                           + zint(j, c, cp)*cf_q(jp, 1:psize, cp)
                     END DO
                  END DO
               END DO

               blk = cmplx_0
               CALL zgemm('C', 'N', nbnd_k, psize, 2*nv2_k, cmplx_1, &
                          cf_k, 2*nv2d, tmat, 2*nv2d, cmplx_0, blk, nbnd_k)

               row = row0 + basfn_offset(mpdata, iq, ivac, igm) + i
               DO n = 1, nbnd_k
                  DO np = 1, psize
                     cprod(row, np + (n - 1)*psize) = cprod(row, np + (n - 1)*psize) + blk(n, np)
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL timestop("wavefproducts_vac")
   END SUBROUTINE wavefproducts_vac

   !>Index in the k+q set of G||_j + gm, 0 if absent.
   SUBROUTINE partner_map(k1, k2, n1, q1, q2, n2, gm, jpart)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: k1(:), k2(:), n1, q1(:), q2(:), n2, gm(2)
      INTEGER, INTENT(OUT) :: jpart(:)

      INTEGER :: j, jp, t1, t2

      jpart = 0
      DO j = 1, n1
         t1 = k1(j) + gm(1)
         t2 = k2(j) + gm(2)
         DO jp = 1, n2
            IF (q1(jp) == t1 .AND. q2(jp) == t2) THEN
               jpart(j) = jp
               EXIT
            END IF
         END DO
      END DO
   END SUBROUTINE partner_map

   REAL FUNCTION z_dot3(f, g, h, nz, delz)
      IMPLICIT NONE
      REAL, INTENT(IN)    :: f(:), g(:), h(:)
      INTEGER, INTENT(IN) :: nz
      REAL, INTENT(IN)    :: delz

      z_dot3 = vac_zint(f(:nz)*g(:nz)*h(:nz), nz, delz)
   END FUNCTION z_dot3

END MODULE m_wavefproducts_vac
