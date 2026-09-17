!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mgga_alpha
   !! MetaGGA diagnostic: extrema of the iso-orbital indicator
   !!
   !!   alpha_sigma = (tau_sigma - tau_W,sigma) / tau_unif,sigma
   !!
   !! with the von Weizsaecker and uniform-gas references
   !!
   !!   tau_W,sigma    = |grad rho_sigma|^2 / (8 rho_sigma)
   !!   tau_unif,sigma = (3/10) (6 pi^2)^(2/3) rho_sigma^(5/3)
   !!
   !! alpha is the quantity SCAN and r2SCAN interpolate on, so it is the sharpest
   !! available probe of the tau construction:
   !!
   !!   alpha = 0  wherever a single orbital determines the density (tau = tau_W),
   !!              which for a one- or two-electron system must hold *everywhere*;
   !!   alpha = 1  for the uniform electron gas.
   !!
   !! A factor-of-two error, a sign error or a missing spin scaling in tau all show
   !! up here immediately, whereas they are nearly invisible in a total energy.
   !!
   !! Note this is evaluated from rho, |grad rho| and tau directly, so it is
   !! independent of the inner-10% GGA fallback that get_exc applies in the muffin
   !! tins - that makes it usable as an exact check.

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: mgga_alpha_extrema

CONTAINS

   SUBROUTINE mgga_alpha_extrema(jspins, rho, sigma, ked, alphaMin, alphaMax)
      USE m_constants
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: jspins
      REAL,    INTENT(IN)    :: rho(:,:)    !! (point, spin) density
      REAL,    INTENT(IN)    :: sigma(:,:)  !! (nsigma, point) contracted gradients
      REAL,    INTENT(IN)    :: ked(:,:)    !! (point, spin) kinetic energy density
      REAL,    INTENT(INOUT) :: alphaMin, alphaMax !! running extrema, updated in place

      INTEGER :: i, js, isig
      REAL    :: rs, taus, grad2, tau_w, tau_unif, alpha
      ! In an exponential density tail tau_unif ~ rho^(5/3) vanishes faster than
      ! (tau - tau_W), so alpha diverges there and its maximum would be dominated by
      ! roundoff in a difference of two nearly equal numbers. Restrict the extrema to
      ! points carrying meaningful density; 1e-6 matches xcpot%apply_cutoffs.
      REAL, PARAMETER :: rho_cut = 1.0e-6
      REAL, PARAMETER :: cunif = 0.3 * (6.0*pi_const**2)**(2.0/3.0)

      DO js = 1, jspins
         ! sigma is stored as (|grad rho_up|^2, grad up . grad dn, |grad rho_dn|^2)
         isig = MERGE(1, 3, js == 1)
         DO i = 1, SIZE(rho, 1)
            IF (jspins == 1) THEN
               ! Spin-scaling: an unpolarized density is two identical half-densities,
               ! so alpha of the total equals alpha of either channel.
               rs    = 0.5 * rho(i, 1)
               taus  = 0.5 * ked(i, 1)
               grad2 = 0.25 * sigma(1, i)
            ELSE
               rs    = rho(i, js)
               taus  = ked(i, js)
               grad2 = sigma(isig, i)
            END IF

            IF (rs <= rho_cut) CYCLE

            tau_w    = grad2 / (8.0 * rs)
            tau_unif = cunif * rs**(5.0/3.0)
            IF (tau_unif <= 0.0) CYCLE

            alpha = (taus - tau_w) / tau_unif

            alphaMin = MIN(alphaMin, alpha)
            alphaMax = MAX(alphaMax, alpha)
         END DO
      END DO

   END SUBROUTINE mgga_alpha_extrema

END MODULE m_mgga_alpha
