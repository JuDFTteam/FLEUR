!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Index bookkeeping of the VAC block (leaf module: must not use m_types).
!>Mixed-basis rows: MT, IR, then VAC ordered (vacuum, G||, z-function); in mtir the moment
!>carriers follow MT+IR ordered (vacuum, G||), then the two g = 0 second carriers.
MODULE m_vac_rows
   USE m_constants, ONLY: cmplx_0
   USE m_vac_const, ONLY: NVAC_MPB
   USE m_types_mpdata
   USE m_types_hybdat
   USE m_coulomb_vac, ONLY: vac_mom_g0
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: row_offset, basfn_offset, vac_g0_moments
   PUBLIC :: NVAC_MPB, vac_src, vac_mtir_idx, vac_mtir_idx2

CONTAINS

   !>Stored vacuum that represents vacuum ivac.
   PURE INTEGER FUNCTION vac_src(nvac, ivac)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nvac, ivac

      vac_src = MIN(ivac, nvac)
   END FUNCTION vac_src

   PURE INTEGER FUNCTION vac_mtir_idx(n_mtir_noVac, n_g_vac_k, ivac, ig)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_mtir_noVac, n_g_vac_k, ivac, ig

      vac_mtir_idx = n_mtir_noVac + (ivac - 1)*n_g_vac_k + ig
   END FUNCTION vac_mtir_idx

   PURE INTEGER FUNCTION vac_mtir_idx2(n_mtir_noVac, n_g_vac_k, nvac, ivac)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_mtir_noVac, n_g_vac_k, nvac, ivac

      vac_mtir_idx2 = n_mtir_noVac + nvac*n_g_vac_k + ivac
   END FUNCTION vac_mtir_idx2

   PURE INTEGER FUNCTION row_offset(mpdata, hybdat, iq, ivac)
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(IN) :: mpdata
      TYPE(t_hybdat), INTENT(IN) :: hybdat
      INTEGER, INTENT(IN)        :: iq, ivac

      INTEGER :: i, n

      n = hybdat%n_mt + mpdata%n_g(iq)
      IF (ivac == 2) THEN
         DO i = 1, mpdata%n_g_vac(iq)
            n = n + mpdata%num_zbasfn_vac(mpdata%glen_ptr_vac(i, iq), 1)
         END DO
      END IF
      row_offset = n
   END FUNCTION row_offset

   PURE INTEGER FUNCTION basfn_offset(mpdata, iq, ivac, igm)
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(IN) :: mpdata
      INTEGER, INTENT(IN)        :: iq, ivac, igm

      INTEGER :: i, n

      n = 0
      DO i = 1, igm - 1
         n = n + mpdata%num_zbasfn_vac(mpdata%glen_ptr_vac(i, iq), ivac)
      END DO
      basfn_offset = n
   END FUNCTION basfn_offset

   !>Charge and first moment of the g = 0 part of a coefficient vector, per vacuum.
   SUBROUTINE vac_g0_moments(mpdata, hybdat, ikpt, nvac, nz, delz, c, qv, mv, igm0)
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(IN) :: mpdata
      TYPE(t_hybdat), INTENT(IN) :: hybdat
      INTEGER, INTENT(IN)        :: ikpt, nvac, nz
      REAL, INTENT(IN)           :: delz
      COMPLEX, INTENT(IN)        :: c(:)
      COMPLEX, INTENT(OUT)       :: qv(:), mv(:)
      INTEGER, INTENT(OUT)       :: igm0

      INTEGER :: ivac, igm, ilen, nn, i, bi
      REAL    :: q_i, zm_i

      qv = cmplx_0; mv = cmplx_0; igm0 = 0
      DO igm = 1, mpdata%n_g_vac(ikpt)
         ilen = mpdata%glen_ptr_vac(igm, ikpt)
         IF (mpdata%glen_vac(ilen) > 1e-12) CYCLE
         igm0 = igm
         DO ivac = 1, nvac
            nn = mpdata%num_zbasfn_vac(ilen, ivac)
            IF (nn == 0) CYCLE
            bi = row_offset(mpdata, hybdat, ikpt, ivac) + basfn_offset(mpdata, ikpt, ivac, igm)
            DO i = 1, nn
               CALL vac_mom_g0(mpdata%zbasfn_vac(:, i, ilen, ivac), nz, delz, q_i, zm_i)
               qv(ivac) = qv(ivac) + c(bi + i)*q_i
               mv(ivac) = mv(ivac) + c(bi + i)*zm_i
            END DO
         END DO
         EXIT
      END DO
   END SUBROUTINE vac_g0_moments

END MODULE m_vac_rows
