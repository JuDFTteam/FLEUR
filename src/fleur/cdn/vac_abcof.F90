!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>Vacuum representation of LAPW eigenvectors, sum_j [ac u_j(z) + bc udot_j(z)] exp(i(k+G_j).r||):
!>in-plane vector map, the 1D solutions and the matching coefficients.  Shared by vacden and the
!>hybrid film code.
MODULE m_vac_abcof
   USE m_juDFT
   USE m_vacuz
   USE m_vacudz
   USE m_types_cell
   USE m_types_lapw
   USE m_types_mat
   USE m_types_vacuum
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: vac_map2, vac_uz, vac_abcof

CONTAINS

   SUBROUTINE vac_map2(lapw, jspin, nv2d, kvac1, kvac2, map2, nv2)
      IMPLICIT NONE
      TYPE(t_lapw), INTENT(IN) :: lapw
      INTEGER, INTENT(IN)      :: jspin, nv2d
      INTEGER, INTENT(INOUT)   :: kvac1(:), kvac2(:), map2(:)
      INTEGER, INTENT(OUT)     :: nv2

      INTEGER :: k, j, n2

      n2 = 0
      k_loop: DO k = 1, lapw%nv(jspin)
         DO j = 1, n2
            IF (lapw%gvec(1, k, jspin) == kvac1(j) .AND. &
                lapw%gvec(2, k, jspin) == kvac2(j)) THEN
               map2(k) = j
               CYCLE k_loop
            END IF
         END DO
         n2 = n2 + 1
         IF (n2 > nv2d) CALL juDFT_error("vac_map2: too many in-plane vectors", &
                                         calledby="vac_map2")
         kvac1(n2) = lapw%gvec(1, k, jspin)
         kvac2(n2) = lapw%gvec(2, k, jspin)
         map2(k) = n2
      END DO k_loop
      nv2 = n2
   END SUBROUTINE vac_map2

   SUBROUTINE vac_uz(vacuum, cell, evacp, vz, bk2, kvac1, kvac2, nv2, &
                     u, ue, t, dt, te, dte, tei)
      IMPLICIT NONE
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_cell), INTENT(IN)   :: cell
      REAL, INTENT(IN)           :: evacp, vz(:), bk2(2)
      INTEGER, INTENT(IN)        :: kvac1(:), kvac2(:), nv2
      REAL, INTENT(INOUT)        :: u(:, :), ue(:, :)
      REAL, INTENT(INOUT)        :: t(:), dt(:), te(:), dte(:), tei(:)

      INTEGER :: ik, j
      REAL    :: v(3), ev, scale, wronk
      REAL    :: ubuf(vacuum%nmz), uebuf(vacuum%nmz)

      wronk = 2.0
      DO ik = 1, nv2
         v(1) = bk2(1) + kvac1(ik)
         v(2) = bk2(2) + kvac2(ik)
         v(3) = 0.

         ev = evacp - 0.5*DOT_PRODUCT(v, MATMUL(v, cell%bbmat))

         CALL vacuz(ev, vz, vz(vacuum%nmz), vacuum%nmz, vacuum%delz, &
                    t(ik), dt(ik), ubuf)
         CALL vacudz(ev, vz, vz(vacuum%nmz), vacuum%nmz, vacuum%delz, te(ik), &
                     dte(ik), tei(ik), uebuf, dt(ik), ubuf)

         scale = wronk/(te(ik)*dt(ik) - dte(ik)*t(ik))
         te(ik) = scale*te(ik)
         dte(ik) = scale*dte(ik)
         tei(ik) = scale*tei(ik)
         DO j = 1, vacuum%nmz
            uebuf(j) = scale*uebuf(j)
         END DO
         u(:vacuum%nmz, ik) = ubuf
         ue(:vacuum%nmz, ik) = uebuf
      END DO
   END SUBROUTINE vac_uz

   SUBROUTINE vac_abcof(cell, lapw, jspin, ivac, nv2d, ne, koff, fac, map2, &
                        t, dt, te, dte, zMat, ac, bc)
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_lapw), INTENT(IN) :: lapw
      INTEGER, INTENT(IN)      :: jspin, ivac, nv2d, ne, koff
      REAL, INTENT(IN)         :: fac
      INTEGER, INTENT(IN)      :: map2(:)
      REAL, INTENT(IN)         :: t(:), dt(:), te(:), dte(:)
      TYPE(t_mat), INTENT(IN)  :: zMat
      COMPLEX, INTENT(INOUT)   :: ac(:, :), bc(:, :)

      INTEGER :: k, ikG, kspin
      REAL    :: sign, zks, arg, const
      COMPLEX :: c_1, av, bv

      const = 1.0/(SQRT(cell%omtil)*2.0)
      sign = 3.-2.*ivac

      DO k = 1, lapw%nv(jspin)
         kspin = koff + k
         ikG = map2(k)
         zks = lapw%k3(k, jspin)*cell%bmat(3, 3)*sign
         arg = zks*cell%z1
         c_1 = CMPLX(COS(arg), SIN(arg))*const
         av = -c_1*CMPLX(dte(ikG), zks*te(ikG))
         bv = c_1*CMPLX(dt(ikG), zks*t(ikG))
         IF (zMat%l_real) THEN
            ac(ikG, :ne) = ac(ikG, :ne) + fac*zMat%data_r(kspin, :ne)*av
            bc(ikG, :ne) = bc(ikG, :ne) + fac*zMat%data_r(kspin, :ne)*bv
         ELSE
            ac(ikG, :ne) = ac(ikG, :ne) + fac*zMat%data_c(kspin, :ne)*av
            bc(ikG, :ne) = bc(ikG, :ne) + fac*zMat%data_c(kspin, :ne)*bv
         END IF
      END DO
   END SUBROUTINE vac_abcof

END MODULE m_vac_abcof
