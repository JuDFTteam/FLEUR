!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The vacuum half of <u_a| e^{-i b.r} p |u_b>, the third region of the momentum. The
!>  muffin-tin and interstitial halves are m_melem_nabla_sph and m_melem_nabla_int; this one
!>  accumulates into the same pnk, so a caller that gains a film gains one more call.
!>
!>  It is deliberately the same routine as m_melem_mmkb_vac with two changes, and a reader
!>  should be able to diff them:
!>
!>    - in plane, the momentum of a two-dimensional plane wave is its own wavevector, so the
!>      x and y components are the overlap integrand times the real constant q_b, which
!>      belongs to THE KET ALONE. Same convention melem_nabla_int states for the
!>      interstitial: the wavevector that multiplies each plane wave is the ket's k + G.
!>
!>    - along z, the momentum differentiates the ket's z function. The vacuum solutions are
!>      real and carry no e^{i k_z z} of their own, so p_z acts on u and ue and on nothing
!>      else. The folding phase e^{-i b_z z} is NOT differentiated: it plays the part the
!>      step function plays in the interstitial, bookkeeping for the two sides being folded
!>      differently, not a factor of the ket.
!>
!>  Returns p and not nabla, like its two siblings, so the caller adds no factor.
!>
!>  THE CONJUGATION. m_melem_mmkb_vac leaves the bra unconjugated and conjugates the ket, so
!>  what it hands over is the complex conjugate of the overlap, and cf_one_k undoes that for
!>  the overlap and the momentum together. This routine is in the same convention, which is
!>  where the +i comes from: p_z = -i d/dz, and conjugating it with real z functions turns
!>  the -i into +i.
!>
!>  SLOT SIGN. The mesh runs away from the film, so on the mirrored slot d/dz is minus
!>  d/d(mesh). slot_sign carries it, the same way it carries the reflection of the phase.
MODULE m_melem_nabla_vac
   USE m_juDFT
   USE m_intgr, ONLY: intgz0
   USE m_types_melem_vacabc, ONLY: t_melem_vacabc
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_nabla_vac

CONTAINS

   SUBROUTINE melem_nabla_vac(vac, vac_b, gb, pnk, vacchi)
      TYPE(t_melem_vacabc), INTENT(IN) :: vac      !< the bra side, already expanded
      TYPE(t_melem_vacabc), INTENT(IN) :: vac_b    !< the ket side
      INTEGER, INTENT(IN) :: gb(3)
      COMPLEX, INTENT(INOUT) :: pnk(:, :, :)       !< (nb, nb_b, 3), cartesian
      !> Spin rotation of the vacuum in a non-collinear run; 1 when there is none.
      COMPLEX, INTENT(IN), OPTIONAL :: vacchi

      INTEGER :: nmz, nb, nb_b, islot, iv, iv_b, l, lp, i, j
      REAL :: sgn, zks, delz, z1, qb(3)
      COMPLEX :: chi, s_uu, s_ud, s_du, s_dd, z_uu, z_ud, z_du, z_dd, cpar, cz
      REAL, ALLOCATABLE :: re(:), im(:)

      nmz = vac%nmz
      delz = vac%delz
      z1 = vac%z1
      nb = SIZE(pnk, 1)
      nb_b = SIZE(pnk, 2)
      chi = CMPLX(1.0, 0.0)
      IF (PRESENT(vacchi)) chi = vacchi

      IF (SIZE(pnk, 3) /= 3) CALL juDFT_error( &
         'melem_nabla_vac: pnk has to carry the three cartesian components', &
         calledby='melem_nabla_vac')
      IF (vac_b%nmz /= nmz .OR. vac_b%delz /= delz .OR. vac_b%z1 /= z1) CALL juDFT_error( &
         'melem_nabla_vac: the two expansions were built on different z meshes', &
         calledby='melem_nabla_vac')
      IF (SIZE(vac%ac, 2) < nb .OR. SIZE(vac_b%ac, 2) < nb_b) CALL juDFT_error( &
         'melem_nabla_vac: fewer expanded states than the momentum asks for', &
         calledby='melem_nabla_vac')
      IF (.NOT. ALLOCATED(vac_b%du) .OR. .NOT. ALLOCATED(vac_b%qpar)) CALL juDFT_error( &
         'melem_nabla_vac: the ket expansion carries no z derivative', &
         hint='t_melem_vacabc builds du/due/qpar; an instance from an older call does not', &
         calledby='melem_nabla_vac')

      CALL timestart('melem_nabla_vac')
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

               !> the overlap integrals, for the in-plane components
               CALL zint(vac%u(:, l, iv), vac_b%u(:, lp, iv_b), s_uu)
               CALL zint(vac%u(:, l, iv), vac_b%ue(:, lp, iv_b), s_ud)
               CALL zint(vac%ue(:, l, iv), vac_b%u(:, lp, iv_b), s_du)
               CALL zint(vac%ue(:, l, iv), vac_b%ue(:, lp, iv_b), s_dd)

               !> the same four with the KET differentiated, for the z component
               CALL zint(vac%u(:, l, iv), vac_b%du(:, lp, iv_b), z_uu)
               CALL zint(vac%u(:, l, iv), vac_b%due(:, lp, iv_b), z_ud)
               CALL zint(vac%ue(:, l, iv), vac_b%du(:, lp, iv_b), z_du)
               CALL zint(vac%ue(:, l, iv), vac_b%due(:, lp, iv_b), z_dd)

               qb = vac_b%qpar(:, lp)

               DO j = 1, nb_b
                  DO i = 1, nb
                     cpar = vac%ac(l, i, islot)*CONJG(vac_b%ac(lp, j, islot))*s_uu + &
                            vac%ac(l, i, islot)*CONJG(vac_b%bc(lp, j, islot))*s_ud + &
                            vac%bc(l, i, islot)*CONJG(vac_b%ac(lp, j, islot))*s_du + &
                            vac%bc(l, i, islot)*CONJG(vac_b%bc(lp, j, islot))*s_dd
                     cz = vac%ac(l, i, islot)*CONJG(vac_b%ac(lp, j, islot))*z_uu + &
                          vac%ac(l, i, islot)*CONJG(vac_b%bc(lp, j, islot))*z_ud + &
                          vac%bc(l, i, islot)*CONJG(vac_b%ac(lp, j, islot))*z_du + &
                          vac%bc(l, i, islot)*CONJG(vac_b%bc(lp, j, islot))*z_dd
                     pnk(i, j, 1) = pnk(i, j, 1) + vac%area*chi*qb(1)*cpar
                     pnk(i, j, 2) = pnk(i, j, 2) + vac%area*chi*qb(2)*cpar
                     pnk(i, j, 3) = pnk(i, j, 3) + vac%area*chi*CMPLX(0.0, sgn)*cz
                  END DO
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (re, im)
      CALL timestop('melem_nabla_vac')

   CONTAINS

      !>  int_{z1}^{inf} f(z) g(z) e^{-i b_z z} dz, the same integral m_melem_mmkb_vac does
      !>  and for the same reason reversed: intgz0 works inward from the outermost point and
      !>  adds the exponential tail beyond it.
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

   END SUBROUTINE melem_nabla_vac

END MODULE m_melem_nabla_vac
