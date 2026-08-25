!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The vacuum half of <u_a| e^{-i b.r} |u_b>, the third and last region of the pair
!>  overlap. The muffin-tin and interstitial halves are m_melem_mmkb_sph and
!>  m_melem_mmkb_int; this one accumulates into the same mmnk, so a caller that gains a
!>  film gains one more call and nothing else.
!>
!>  Only parallel G vectors that differ by the parallel part of b contribute: outside the
!>  film the two-dimensional plane waves are still orthogonal, so the whole overlap is a
!>  sum of z integrals over matching pairs. Each pair needs four of them, one per product
!>  of the two z solutions, and the phase e^{-i b_z z} rides along inside the integrand
!>  rather than outside it.
!>
!>  The two sides are independent -- their own t_melem_vacabc, their own spin -- and gb is
!>  the reciprocal lattice vector that brings the ket's k-point back into the zone. That
!>  the caller happens to use them with k and k+b is its business, not this routine's.
!>
!>  BOTH SIDES OF THE FILM, ALWAYS. The loop is over the two slots of t_melem_vacabc, not
!>  over vacuum%nvac: when mirror symmetry leaves only one side stored, the second slot is
!>  its reflection, and this routine does not need to know which case it is in.
!>
!>  The area is the in-plane measure of the integral, and it is explicit here because
!>  t_melem_vacabc leaves it out of ac/bc -- see the note on normalisation in that module
!>  before comparing anything with eigen/hsvac.F90.
MODULE m_melem_mmkb_vac
   USE m_juDFT
   USE m_types
   USE m_intgr, ONLY: intgz0
   USE m_types_melem_vacabc, ONLY: t_melem_vacabc
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_mmkb_vac

CONTAINS

   SUBROUTINE melem_mmkb_vac(vac, vac_b, gb, mmnk, vacchi)
      TYPE(t_melem_vacabc), INTENT(IN) :: vac      !< the bra side, already expanded
      TYPE(t_melem_vacabc), INTENT(IN) :: vac_b    !< the ket side
      INTEGER, INTENT(IN) :: gb(3)
      COMPLEX, INTENT(INOUT) :: mmnk(:, :)
      !> Spin rotation of the vacuum in a non-collinear run; 1 when there is none.
      COMPLEX, INTENT(IN), OPTIONAL :: vacchi

      INTEGER :: nmz, nb, nb_b, islot, iv, iv_b, l, lp, i, j
      REAL :: sgn, zks, delz, z1
      COMPLEX :: chi, uu, ud, du, dd
      REAL, ALLOCATABLE :: re(:), im(:)

      nmz = vac%nmz
      delz = vac%delz
      z1 = vac%z1
      nb = SIZE(mmnk, 1)
      nb_b = SIZE(mmnk, 2)
      chi = CMPLX(1.0, 0.0)
      IF (PRESENT(vacchi)) chi = vacchi

      IF (vac_b%nmz /= nmz .OR. vac_b%delz /= delz .OR. vac_b%z1 /= z1) CALL juDFT_error( &
         'melem_mmkb_vac: the two expansions were built on different z meshes', &
         calledby='melem_mmkb_vac')
      IF (SIZE(vac%ac, 2) < nb .OR. SIZE(vac_b%ac, 2) < nb_b) CALL juDFT_error( &
         'melem_mmkb_vac: fewer expanded states than the overlap asks for', &
         calledby='melem_mmkb_vac')

      CALL timestart('melem_mmkb_vac')
      ALLOCATE (re(nmz), im(nmz))

      DO islot = 1, 2
         sgn = vac%slot_sign(islot)
         iv = vac%slot_vac(islot)
         iv_b = vac_b%slot_vac(islot)
         zks = gb(3)*vac%bmat33*sgn

         DO l = 1, vac%nv2
            DO lp = 1, vac_b%nv2
               !> the ket's parallel G, brought back into the zone, has to match the bra's
               IF (vac%kvac(1, l) /= vac_b%kvac(1, lp) - gb(1)) CYCLE
               IF (vac%kvac(2, l) /= vac_b%kvac(2, lp) - gb(2)) CYCLE

               CALL zint(vac%u(:, l, iv), vac_b%u(:, lp, iv_b), uu)
               CALL zint(vac%u(:, l, iv), vac_b%ue(:, lp, iv_b), ud)
               CALL zint(vac%ue(:, l, iv), vac_b%u(:, lp, iv_b), du)
               CALL zint(vac%ue(:, l, iv), vac_b%ue(:, lp, iv_b), dd)

               DO j = 1, nb_b
                  DO i = 1, nb
                     mmnk(i, j) = mmnk(i, j) + vac%area*chi*( &
                                  vac%ac(l, i, islot)*CONJG(vac_b%ac(lp, j, islot))*uu + &
                                  vac%ac(l, i, islot)*CONJG(vac_b%bc(lp, j, islot))*ud + &
                                  vac%bc(l, i, islot)*CONJG(vac_b%ac(lp, j, islot))*du + &
                                  vac%bc(l, i, islot)*CONJG(vac_b%bc(lp, j, islot))*dd)
                  END DO
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (re, im)
      CALL timestop('melem_mmkb_vac')

   CONTAINS

      !>  int_{z1}^{inf} f(z) g(z) e^{-i b_z z} dz, on the linear mesh of the vacuum.
      !>
      !>  The mesh is reversed before integrating because intgz0 works inward from the
      !>  outermost point and adds the exponential tail beyond it -- that tail is the part
      !>  of the vacuum the mesh does not reach, and dropping it would quietly truncate
      !>  every overlap.
      SUBROUTINE zint(f, g, res)
         REAL, INTENT(IN) :: f(:), g(:)
         COMPLEX, INTENT(OUT) :: res
         INTEGER :: iz
         REAL :: zr, zi, arg

         DO iz = 1, nmz
            arg = (z1 + (iz - 1)*delz)*zks
            re(nmz + 1 - iz) = f(iz)*g(iz)*COS(arg)
            im(nmz + 1 - iz) = f(iz)*g(iz)*SIN(arg)
         END DO
         CALL intgz0(re, delz, nmz, zr, .TRUE.)
         CALL intgz0(im, delz, nmz, zi, .TRUE.)
         res = CMPLX(zr, zi)
      END SUBROUTINE zint

   END SUBROUTINE melem_mmkb_vac

END MODULE m_melem_mmkb_vac
