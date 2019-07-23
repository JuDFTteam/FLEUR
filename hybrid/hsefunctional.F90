!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

! Provides all procedures needed to implement the HSE-exchange functional
! within a FLAPW framework
! Author: M. Schlipf 2009
MODULE m_hsefunctional
   USE m_judft
   IMPLICIT NONE

#ifdef __PGI
   REAL, EXTERNAL ::erfc
#endif
   ! Constant omega of the HSE exchange functional
   REAL, PARAMETER :: omega_HSE = 0.11

   ! Constant for the maximum number of G points
   INTEGER, PARAMETER :: maxNoGPts = 50

   ! these arrays are calculated once and then reused
   LOGICAL, ALLOCATABLE   :: already_known(:)
   REAL, ALLOCATABLE      :: known_potential(:, :)
#ifdef CPP_INVERSION
   REAL, ALLOCATABLE      :: known_fourier_trafo(:, :, :)
#else
   COMPLEX, ALLOCATABLE   :: known_fourier_trafo(:, :, :)
#endif

   ! variables needed for fourier transformation
   PRIVATE already_known, known_potential, known_fourier_trafo

   ! functions needed internally
   PRIVATE calcYlm, calcSphBes, my_dot_product3x1, my_dot_product3x3, my_dot_product4x1, &
      my_dot_product4x4, my_sum3d, my_sum4d, gPtsSummation, approximateIntegral, generateIntegrals, &
      integrateY1ExpErfc, integrateY3ExpErfc, integrateY5ExpErfc, integrateY7ExpErfc, &
      integrateY9ExpErfc, calculateF, calculateG, calculateH

   INTERFACE my_sum
      MODULE PROCEDURE my_sum3d, my_sum4d
   END INTERFACE

   INTERFACE my_dot_product
      MODULE PROCEDURE my_dot_product4x1, my_dot_product4x4, my_dot_product3x1, my_dot_product3x3
   END INTERFACE

   INTERFACE gPtsSummation
      MODULE PROCEDURE crc_gPtsSummation, rrr_gPtsSummation
   END INTERFACE

CONTAINS

   ! Calculate the enhancement factor of the HSE03 functional
   ! References:
   ! [1] Heyd, Scuseria, Ernzerhof: Hybrid functionals based on a screened Coulomb potential,
   !     J. Chem. Phys. 118 (2003) 8207-8215
   ! [2] Heyd, Scuseria: Assessment and validation of a screened Coulomb hybrid density
   !     functional, J. Chem. Phys. 120 (2004) 7274-7280
   ! [3] Ernzerhof, Perdew: Generalized gradient approximation to the angle- and system-
   !     averaged exchange hole, J. Chem. Phys. 109 (1998) 3313-3319
   ! Input:  kF       - kF = (3 pi^2 rho) ** 1/3
   !         s_inp    - reduced gradient = |grad rho| / 2 kF rho
   ! Output: F_x      - HSE03 enhancement factor
   !         dFx_ds   - derivative of this factor with respect to s
   !         d2Fx_ds2 - second derivative with respect to s
   SUBROUTINE calculateEnhancementFactor(kF, s_inp, F_x, dFx_Ds, d2Fx_Ds2, dFx_dkF, d2Fx_dsdkF)

      IMPLICIT NONE

      REAL, INTENT(IN)  :: kF, s_inp
      REAL, INTENT(OUT) :: F_x, dFx_ds, d2Fx_ds2
      REAL, INTENT(OUT) :: dFx_dkF, d2Fx_dsdkF

      ! Helper variables
      REAL :: r1_kF, &                   ! 1 / kF
              omega_kF, omega_kF_Sqr, &      ! omega_HSE / kF, (omega/kF)^2
              s, s2, &                      ! reduced gradient s and s^2
              s_si2, r2s_si3, r6s_si4, &    ! quotient s_chg / s_inp^n (n = 2,3,4)
              F, G, H, &                    ! results of the functions F, G, and H defined in [3]
              Hs2, dHs2_ds, d2Hs2_ds2, &    ! H * s^2 and derivatives with respect to s
              Fs2, dFs2_ds, d2Fs2_ds2, &    ! F * s^2 and derivatives with respect to s
              dGs2_ds, d2Gs2_ds2, &         ! derivatives with respect to s of G s^2
              C_term, dCt_ds, d2Ct_ds2, &   ! C * (1 + F(s) s^2) and derivatives with respect to s
              E_term, dEt_ds, d2Et_ds2, &   ! E * (1 + G(s) s^2) and derivatives with respect to s
              D_term                        ! D + H s^2

      ! integrals of y^n Exp(-ay^2) Erfc(w/k y) for n = 1, 3, 5, 7, 9
      REAL :: intY1ExpErfc, intY3ExpErfc, intY5ExpErfc, intY7ExpErfc, intY9ExpErfc
      ! integrals of 2/sqrt(pi) y^n Exp(-arg*y^2) for n = 2, 4, 6, 8
      REAL :: r1_arg, intY2Gauss, intY4Gauss, intY6Gauss, intY8Gauss

      ! approximation of the not analytical part of the integral and derivatives
      REAL :: appInt, dAppInt_ds, d2AppInt_ds2
      REAL :: dAppInt_dkF, d2AppInt_dsdkF

      REAL, PARAMETER :: r8_9 = 8.0/9.0  ! 8/9
      REAL, PARAMETER :: &               ! Parameters of the exchange hole as defined in [3]
         B = -0.37170836, C = -0.077215461, &
         D = 0.57786348, E = -0.051955731

      REAL, PARAMETER :: &               ! Correction of the reduced gradient to ensure Lieb-Oxford bound [2]
         s_thresh = 8.3, s_max = 8.572844, s_chg = 18.796223
      LOGICAL :: correction              ! Has the value of s been corrected

      ! If a large value of s would violate the Lieb-Oxford bound, the value of s is reduced,
      ! so that this condition is fullfilled
      correction = s_inp > s_thresh
      IF (correction) THEN
         s_si2 = s_chg/(s_inp*s_inp)
         s = s_max - s_si2
      ELSE
         s = s_inp
      END IF

      ! Calculate different helper variables
      r1_kF = 1.0/kF
      omega_kF = omega_hse*r1_kF !was omega_VHSE()
      omega_kF_Sqr = omega_kF*omega_kF

      ! calculate the functions H and F in [3] and its derivatives
      CALL calculateH(s, H, dHs2_ds, d2Hs2_ds2)
      CALL calculateF(s, H, dHs2_ds, d2Hs2_ds2, F, dFs2_ds, d2Fs2_ds2)
      s2 = s*s
      Hs2 = H*s2
      Fs2 = F*s2
      CALL calculateG(s2, Fs2, dFs2_ds, d2Fs2_ds2, Hs2, dHs2_ds, d2Hs2_ds2, & !Input
                      G, dGs2_ds, d2Gs2_ds2)                                 !Output

      C_term = C*(1 + s2*F)
      dCt_ds = C*dFs2_ds
      d2Ct_ds2 = C*d2Fs2_ds2
      E_term = E*(1 + s2*G)
      dEt_ds = E*dGs2_ds
      d2Et_ds2 = E*d2Gs2_ds2
      D_term = D + Hs2

      ! approximate the integral using an expansion of the error function (cf. [2])
      CALL approximateIntegral(omega_kF, Hs2, D_term, dHs2_ds, d2Hs2_ds2, & !Input
                               appInt, dAppInt_ds, d2AppInt_ds2, dAppInt_dkF, d2AppInt_dsdkF)  !Output

      ! Calculate the integrals
      !
      ! inf
      !   /                2   2      /       \
      !  |     n  -(D+H(s)s ) y      | omega   |
      !  | dy y  e               Erfc| ----- y |
      !  |                           |   k     |
      ! /                             \   F   /
      !  0
      ! for n = 1, 3, 5, 7, and 9
      !
      intY1ExpErfc = integrateY1ExpErfc(omega_kf, omega_kf_Sqr, D_term)
      intY3ExpErfc = integrateY3ExpErfc(omega_kf, omega_kf_Sqr, D_term)
      intY5ExpErfc = integrateY5ExpErfc(omega_kf, omega_kf_Sqr, D_term)
      intY7ExpErfc = integrateY7ExpErfc(omega_kf, omega_kf_Sqr, D_term)
      intY9ExpErfc = integrateY9ExpErfc(omega_kf, omega_kf_Sqr, D_term)

      ! Calculate the integrals
      !
      !       inf
      !         /          /                      2    \
      !   2    |     n    |          2   2   omega   2  |
      ! -----  | dy y  Exp| -(D+H(s)s ) y  - ------ y   |
      !  ____  |          |                   k ^2      |
      ! V Pi  /            \                   F       /
      !        0
      ! for n = 2, 4, 6, 8
      !
      r1_arg = 1.0/(D_term + omega_kF_Sqr)
      intY2Gauss = 0.5*SQRT(r1_arg)*r1_arg
      intY4Gauss = 1.5*intY2Gauss*r1_arg
      intY6Gauss = 2.5*intY4Gauss*r1_arg
      intY8Gauss = 3.5*intY6Gauss*r1_arg

      ! Calculate the integral
      !  inf
      !    /                 /       \
      !   |                 | omega   |
      !   | dy y J(s,y) Erfc| ----- y |
      !   |                 |   k     |
      !  /                   \   F   /
      !   0
      ! where J(s, y) is the exchange hole defined in [3]
      ! the exchange factor is proportional to this integral
      F_x = -r8_9*(appInt + B*intY1ExpErfc + C_term*intY3ExpErfc + E_term*intY5ExpErfc)

      ! Calculate the derivatives with respect to s using that the derivatative of the integral
      ! yields higher orders of the same kind of integral intY1 -> -intY3 -> intY5 ... times
      ! the derivative of the exponent
      dFx_ds = -r8_9*(dAppInt_ds - (B*intY3ExpErfc + C_term*intY5ExpErfc + E_term*intY7ExpErfc)*dHs2_ds &
                      + dCt_ds*intY3ExpErfc + dEt_ds*intY5ExpErfc)
      d2Fx_ds2 = -r8_9*(d2AppInt_ds2 + (B*intY5ExpErfc + C_term*intY7ExpErfc + E_term*intY9ExpErfc)*dHs2_ds**2 &
                        - (B*intY3ExpErfc + C_term*intY5ExpErfc + E_term*intY7ExpErfc)*d2Hs2_ds2 &
                        - 2.0*(dCt_ds*intY5ExpErfc + dEt_ds*intY7ExpErfc)*dHs2_ds &
                        + d2Ct_ds2*intY3ExpErfc + d2Et_ds2*intY5ExpErfc)

      dFx_dkF = -r8_9*r1_kF*(omega_kF*(B*intY2Gauss + C_term*intY4Gauss + E_term*intY6Gauss) + dAppInt_dkF)
      d2Fx_dsdkF = -r8_9*r1_kF*(d2AppInt_dsdkF + omega_kF*(dCt_ds*intY4Gauss + dEt_ds*intY6Gauss &
                                                           - (B*intY4Gauss + C_term*intY6Gauss + E_term*intY8Gauss)*dHs2_ds))

      ! Correction to the derivatives, if s_inp > s_thresh
      IF (correction) THEN
         r2s_si3 = 2.0*s_si2/s_inp
         r6s_si4 = 3.0*r2s_si3/s_inp
         d2Fx_ds2 = d2Fx_ds2*r2s_si3**2 - dFx_ds*r6s_si4
         d2Fx_dsdkF = d2Fx_dsdkF*r2s_si3
         dFx_ds = dFx_ds*r2s_si3
      END IF

   END SUBROUTINE calculateEnhancementFactor

   ! Approximation for the first part of the exchange hole integral
   ! Calculated using the algorithm described in [2].
   ! Additionally the first and second derivative of this approximation with
   ! respect to s are calculated
   ! Input:  omega_kF       - omega / kF
   !         Hs2            - H s^2
   !         D_Hs2          - D + H s^2
   !         dHs2_ds        - d (H s^2) / ds
   !         d2Hs2_ds2      - d^2 (H s^2) / ds^2
   ! Output: appInt         - approximation for one part of the exchange hole integral
   !         dAppInt_ds     - first derivative with respect to s
   !         d2AppInt_ds2   - second derivative with respect to s
   !         dAppInt_dkF    - first derivative with respect to kF
   !         d2AppInt_dsdkF - mixed derivative with respect to s and kF
   SUBROUTINE approximateIntegral(omega_kF, Hs2, D_Hs2, dHs2_ds, d2Hs2_ds2, &
                                  appInt, dAppInt_ds, d2AppInt_ds2, dAppInt_dkF, d2AppInt_dsdkF)

      USE m_exponential_integral, ONLY: calculateExponentialIntegral, gauss_laguerre

      IMPLICIT NONE

      REAL, INTENT(IN)  :: omega_kF, Hs2, D_Hs2, dHs2_ds, d2Hs2_ds2
      REAL, INTENT(OUT) :: appInt, dAppInt_ds, d2AppInt_ds2
      REAL, INTENT(OUT) :: dAppInt_dkF, d2AppInt_dsdkF

      REAL    :: w2, bw, r2bw, bw_Hs2, bw_D_Hs2
      ! variables for temporary storage of the integrals, the prefactors and the dot_product
      REAL    :: integral(0:12), a_omegaI(0:8), aI_omegaI(0:8), dotpr, dotpr2
      INTEGER :: i

      REAL    :: r2w2, r4w2, r2w2_Hs2, r2w2_D_Hs2, arg, exp_e1, dAppInt_dh, d2AppInt_dh2

      ! parameters of the erfc fit as given in [2]
      REAL, PARAMETER :: A_2 = 0.5080572, scale = 1.125/A_2, &
                         a(1:8) = (/-1.128223946706117, 1.452736265762971, &
                                    -1.243162299390327, 0.971824836115601, &
                                    -0.568861079687373, 0.246880514820192, &
                                    -0.065032363850763, 0.008401793031216/), &
                         b = 1.455915450052607, cutoff = 14.0

      ! Calculate helper variables
      w2 = omega_kF**2
      bw = b*w2
      r2bw = 2.0*bw
      bw_Hs2 = bw + Hs2
      bw_D_Hs2 = bw + D_Hs2

      IF (bw_Hs2 < cutoff) THEN
         ! integrate the erfc approximation times the exchange hole
         CALL generateIntegrals(bw_Hs2, bw_D_Hs2, integral)

         ! combine the solutions of the integrals with the appropriate prefactors
         a_omegaI(0) = 1.0
         a_omegaI(1:8) = a*(/(omega_kF**i, i=1, 8)/)
         aI_omegaI = (/(i*a_omegaI(i), i=0, 8)/)
         appInt = DOT_PRODUCT(a_omegaI, integral(0:8))
         dotpr = DOT_PRODUCT(a_omegaI, integral(2:10))
         dAppInt_ds = -dotpr*dHs2_ds
         dAppInt_dkF = r2bw*dotpr - DOT_PRODUCT(aI_omegaI, integral(0:8))
         dotpr2 = DOT_PRODUCT(a_omegaI, integral(4:12))
         d2AppInt_ds2 = dotpr2*dHs2_ds**2 - dotpr*d2Hs2_ds2
         d2AppInt_dsdkF = -(r2bw*dotpr2 - DOT_PRODUCT(aI_omegaI, integral(2:10)))*dHs2_ds
      ELSE
         r2w2 = 2.0*w2
         r4w2 = 4.0*w2
         r2w2_Hs2 = r2w2 + Hs2
         r2w2_D_Hs2 = r2w2 + D_Hs2
         arg = scale*r2w2_Hs2
         exp_e1 = gauss_laguerre(arg)
         appInt = A_2*(LOG(r2w2_Hs2/r2w2_D_Hs2) + exp_e1)
         dAppInt_dh = -A_2/r2w2_D_Hs2 + 1.125*exp_e1
         d2AppInt_dh2 = 1.125/r2w2_Hs2 - A_2/r2w2_D_Hs2**2 - 1.265625/A_2*exp_e1
         dAppInt_ds = dAppInt_dh*dHs2_ds
         dAppInt_dkF = -dAppInt_dh*r4w2
         d2AppInt_ds2 = d2AppInt_dh2*dHs2_ds**2 + dAppInt_dh*d2Hs2_ds2
         d2AppInt_dsdkF = d2AppInt_dh2*dHs2_ds*r4w2
      END IF

   END SUBROUTINE approximateIntegral

   ! Generate the integrals for n = 0, 1, ..., 12
   !
   ! inf
   !   /        /        2                  \      /                        \
   !  |     n  |  A   -Dy           A        |    |    /   omega      2\   2 |
   !  | dy y   | --- e     - --------------  | Exp| - (  b ----- + H s  ) y  |
   !  |        |  y            /    4   2\   |    |    \     k         /     |
   ! /          \            y( 1 + - Ay  ) /      \          F             /
   !  0                        \    9    /
   !
   ! Input:  bw_Hs2   - b (omega/kF)^2 + H s^2
   !         bw_D_Hs2 - b (omega/kF)^2 + D + H s^2
   ! Output: integral - array with the calculated integrals
   ! To simplify the calculation use integral(n+2) = - d(integral(n))/d(b omega/kF)
   SUBROUTINE generateIntegrals(bw_Hs2, bw_D_Hs2, integral)
      USE m_exponential_integral, ONLY: calculateExponentialIntegral, gauss_laguerre, series_laguerre
      USE m_constants
      IMPLICIT NONE

      REAL, INTENT(IN)  :: bw_Hs2, bw_D_Hs2
      REAL, INTENT(OUT) :: integral(0:12)

      ! Helper variables
      REAL :: bw_Hs2_Sqr, bw_Hs2_Cub, sqrt_bw_Hs2, &
              bw_D_Hs2_Sqr, bw_D_Hs2_Cub, bw_D_Hs2_Tet, sqrt_bw_D_Hs2, &
              arg, sqrt_arg, r1_arg, e1_arg, exp_arg, exp_erfc, &
              term1, factor2, term2, arg_n, sum_term, add_term, half_i2_1
      INTEGER :: i

      ! A is an exchange hole parameter [3] and b the fit parameter for erfc [2]
      REAL, PARAMETER :: &
         A = 1.0161144, A_2 = A/2.0, scale = 2.25/A, sqrtA = 1.008025, & !sqrt(A)
         b = 1.455915450052607

      ! Calculate many helper variables
      bw_Hs2_Sqr = bw_Hs2*bw_Hs2
      bw_Hs2_Cub = bw_Hs2*bw_Hs2_Sqr
      sqrt_bw_Hs2 = SQRT(bw_Hs2)
      bw_D_Hs2_Sqr = bw_D_Hs2*bw_D_Hs2
      bw_D_Hs2_Cub = bw_D_Hs2*bw_D_Hs2_Sqr
      bw_D_Hs2_Tet = bw_D_Hs2_Sqr*bw_D_Hs2_Sqr
      sqrt_bw_D_Hs2 = SQRT(bw_D_Hs2)

      ! arg = 9/4 * (b omega/kF + H s^2) / A
      arg = scale*bw_Hs2
      sqrt_arg = SQRT(arg)
      r1_arg = 1.0/arg

      ! calculate e^(arg), E1(arg), and erfc(sqrt(arg))
      exp_arg = EXP(arg)
      exp_erfc = exp_arg*erfc(sqrt_arg)

      IF (arg > series_laguerre) THEN
         term2 = gauss_laguerre(arg)
      ELSE
         CALL calculateExponentialIntegral(arg, e1_arg)
         term2 = exp_arg*e1_arg
      END IF

      ! The n = 0 integral is
      ! A/2 ( ln((b omega/kF + H s^2) / (b omega/kF + D + H s^2))
      !     + e^(arg) E1(arg) )
      integral(0) = A_2*(LOG(bw_Hs2/bw_D_Hs2) + term2)

      ! Calculate now all even n's by successive derivation
      ! The log(...) term gives term proportional to 1/(b omega/kF + D + H s^2)^i
      ! The e^(arg) E1(arg) term reproduces itself with a prefactor and
      ! generates an additional 1/arg term which produces higher 1/arg^i terms
      ! when further deriviated
      term1 = A_2/bw_D_Hs2
      factor2 = -1.125
      arg_n = -1.0/arg
      integral(2) = term1 + factor2*term2

      DO i = 1, 5
         term1 = term1/bw_D_Hs2*REAL(i)
         factor2 = -factor2*scale
         term2 = term2 + arg_n

         integral(2*(i + 1)) = term1 + factor2*term2

         arg_n = -arg_n*REAL(i)/arg
      END DO

      ! The n = 1 integral is
      ! A/2 ( sqrt(pi) / sqrt( b omega/kF + D + H s2 )
      ! - 3/4 sqrt(A) pi e^(arg) erfc(sqrt(arg))
      term1 = A_2*SQRT(PI_const)/sqrt_bw_D_Hs2
      term2 = PI_const*exp_erfc
      factor2 = -0.75*sqrtA

      integral(1) = term1 + factor2*term2

      ! Calculate now all uneven n's by successive derivation
      ! The 1 / sqrt(...) term gives higher orders of 1 / (...)^((2i+1)/2)
      ! The e^(arg) erfc(sqrt(arg)) term reproduces itself with a prefactor
      ! and generates an additional 1/sqrt(arg) term which produces higher
      ! 1/(arg)^((2i+1)/2) terms when further deriviated
      sum_term = -1.125*SQRT(pi_const)/sqrt_bw_Hs2
      add_term = sum_term
      half_i2_1 = -0.5

      DO i = 3, 11, 2
         factor2 = -factor2*scale
         term1 = -term1*half_i2_1/bw_D_Hs2
         integral(i) = term1 + term2*factor2 + sum_term

         add_term = -add_term*half_i2_1/bw_Hs2
         sum_term = -sum_term*scale + add_term
         half_i2_1 = half_i2_1 - 1.0
      ENDDO

   END SUBROUTINE generateIntegrals

   ! Calculate the integral
   !
   ! inf
   !   /                  /       \         /                 ________________ \
   !  |       -a y^2     | omega   |    1  |     omega   /   /     / omega \2   |
   !  | dy y e       Erfc| ----- y | = --- | 1 - -----  /   / a + (  -----  )   |
   !  |                  |   k     |   2 a |       k   /   V       \   k   /    |
   ! /                    \   F   /         \       F                   F      /
   !  0
   !
   ! Input:  omega_kF     - omega / kF
   !         omega_kF_Sqr - (omega / kF)^2
   !         a            - factor of the gaussian function
   ! Return: result of the integration
   REAL FUNCTION integrateY1ExpErfc(omega_kF, omega_kF_Sqr, a)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: omega_kF, omega_kF_Sqr, a

      ! helper variable sqrt(a + (omega/kF)^2)
      REAL :: sqrt_term

      ! calculate helper variable
      sqrt_term = SQRT(a + omega_kF_Sqr)

      ! calculate the integral
      integrateY1ExpErfc = 0.5*(1 - omega_kF/sqrt_term)/a

   END FUNCTION integrateY1ExpErfc

   ! Calculate the integral
   !
   ! inf
   !   /                   /       \           /                   ________________                   ________________3 \
   !  |     3  -a y^2     | omega   |     1   |       omega   /   /     / omega \2       omega   /   /     / omega \2    |
   !  | dy y  e       Erfc| ----- y | = ----- | 2 - 2 -----  /   / a + (  -----  )  - a  -----  /   / a + (  -----  )    |
   !  |                   |   k     |   4 a^2 |         k   /   V       \   k   /          k   /   V       \   k   /     |
   ! /                     \   F   /           \         F                   F              F                   F       /
   !  0
   !
   ! Input:  omega_kF     - omega / kF
   !         omega_kF_Sqr - (omega / kF)^2
   !         a            - factor of the gaussian function
   ! Return: result of the integration
   REAL FUNCTION integrateY3ExpErfc(omega_kF, omega_kF_Sqr, a)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: omega_kF, omega_kF_Sqr, a

      ! helper variables
      REAL :: term, sqrt_term, sqrt_term_3

      ! calculate helper variables
      term = a + omega_kF_Sqr
      sqrt_term = SQRT(term)
      sqrt_term_3 = term*sqrt_term

      ! calculate the integral
      integrateY3ExpErfc = 0.25*(2.0 - a*omega_kF/sqrt_term_3 - 2.0*omega_kF/sqrt_term)/(a*a)

   END FUNCTION integrateY3ExpErfc

   ! Calculate the integral
   !
   ! inf                                     /                                   \
   !   /                   /       \        |                  /              2\  |
   !  |     5  -a y^2     | omega   |    1  |       omega     |      a     3 a  | |
   !  | dy y  e       Erfc| ----- y | = --- | 1 - ----------  | 1 + --- + ----- | |
   !  |                   |   k     |   a^3 |     k  sqrt(x)  |     2 x   8 x^2 | |
   ! /                     \   F   /        |      F           \               /  |
   !  0                                      \                                   /
   ! where
   !          /     \ 2
   !         | omega |
   ! x = a + | ----- |
   !         |   k   |
   !          \   F /
   !
   ! Input:  omega_kF     - omega / kF
   !         omega_kF_Sqr - (omega / kF)^2
   !         a            - factor of the gaussian function
   ! Return: result of the integration
   REAL FUNCTION integrateY5ExpErfc(omega_kF, omega_kF_Sqr, a)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: omega_kF, omega_kF_Sqr, a

      ! helper variables
      REAL :: term, sqrt_term, a_Sqr

      ! calculate helper variables
      term = a + omega_kF_Sqr
      sqrt_term = SQRT(term)
      a_Sqr = a*a

      ! calculate the integral
      integrateY5ExpErfc = (1.0 - (0.375*a_Sqr/(term*term) + a/(2.0*term) + 1.0)*omega_kF/sqrt_term)/(a_Sqr*a)

   END FUNCTION integrateY5ExpErfc

   ! Calculate the integral
   !
   ! inf                                     /                                            \
   !   /                   /       \        |                  /              2       3 \  |
   !  |     7  -a y^2     | omega   |    3  |       omega     |      a     3 a     5 a   | |
   !  | dy y  e       Erfc| ----- y | = --- | 1 - ----------  | 1 + --- + ----- + ------ | |
   !  |                   |   k     |   a^4 |     k  sqrt(x)  |     2 x   8 x^2   16 x^3 | |
   ! /                     \   F   /        |      F           \                        /  |
   !  0                                      \                                            /
   ! where
   !          /     \ 2
   !         | omega |
   ! x = a + | ----- |
   !         |   k   |
   !          \   F /
   !
   ! Input:  omega_kF     - omega / kF
   !         omega_kF_Sqr - (omega / kF)^2
   !         a            - factor of the gaussian function
   ! Return: result of the integration
   REAL FUNCTION integrateY7ExpErfc(omega_kF, omega_kF_Sqr, a)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: omega_kF, omega_kF_Sqr, a

      ! helper variables
      REAL :: term2, term, sqrt_term, a_Sqr

      ! calculate helper variables
      term = a + omega_kF_Sqr
      term2 = term*term
      sqrt_term = SQRT(term)
      a_Sqr = a*a

      ! calculate the integral
      integrateY7ExpErfc = 3.0/(a_Sqr*a_Sqr)*(1.0 - omega_kF/sqrt_term* &
                                              (0.3125*(a_Sqr*a)/(term2*term) + 0.375*a_Sqr/term2 + a/(2.0*term) + 1.0))

   END FUNCTION integrateY7ExpErfc

   ! Calculate the integral
   !
   ! inf                                     /                                                      \
   !   /                   /       \        |                  /              2       3           \  |
   !  |     9  -a y^2     | omega   |    12 |       omega     |      a     3 a     5 a      35 a^4 | |
   !  | dy y  e       Erfc| ----- y | = --- | 1 - ----------  | 1 + --- + ----- + ------ + ------- | |
   !  |                   |   k     |   a^5 |     k  sqrt(x)  |     2 x   8 x^2   16 x^3   128 x^4 | |
   ! /                     \   F   /        |      F           \                                  /  |
   !  0                                      \                                                      /
   ! where
   !          /     \ 2
   !         | omega |
   ! x = a + | ----- |
   !         |   k   |
   !          \   F /
   !
   ! Input:  omega_kF     - omega / kF
   !         omega_kF_Sqr - (omega / kF)^2
   !         a            - factor of the gaussian function
   ! Return: result of the integration
   REAL FUNCTION integrateY9ExpErfc(omega_kF, omega_kF_Sqr, a)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: omega_kF, omega_kF_Sqr, a

      ! helper variables
      REAL :: term2, term, sqrt_term, a_Sqr, a_Tet

      ! calculate helper variables
      term = a + omega_kF_Sqr
      term2 = term*term
      sqrt_term = SQRT(term)
      a_Sqr = a*a
      a_Tet = a_Sqr*a_Sqr

      ! calculate the integral
      integrateY9ExpErfc = 12.0*(1.0 - (0.2734375*a_Tet/(term2*term2) + 0.3125*(a_Sqr*a)/(term2*term) &
                                        + 0.375*a_Sqr/term2 + a/(2.0*term) + 1.0)*omega_kF/sqrt_term)/(a_Tet*a)

   END FUNCTION integrateY9ExpErfc

   ! Calculate the function F(s) given in [3]
   ! Input:  s         - reduced gradient
   !         H         - value of the function H(s)
   !         dHs2_ds   - first derivative of H(s)s^2
   !         dsHs2_ds2 - second derivative of H(s)s^2
   ! Output: F         - value of the function F
   !         dFs2_ds   - first derivative of F(s)s^2
   !         d2Fs2_ds2 - second derivative of F(s)s^2
   SUBROUTINE calculateF(s, H, dHs2_ds, d2Hs2_ds2, F, dFs2_ds, d2Fs2_ds2)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: s, H, dHs2_ds, d2Hs2_ds2
      REAL, INTENT(OUT) :: F, dFs2_ds, d2Fs2_ds2

      REAL, PARAMETER   :: slope = 6.4753871
      REAL, PARAMETER   :: shift = 0.4796583

      ! calculate the function F(s), d(s^2F(s))/ds, and d^2(s^2F(s))/ds^2
      F = slope*H + shift
      dFs2_ds = slope*dHs2_ds + 2.0*s*shift
      d2Fs2_ds2 = slope*d2Hs2_ds2 + 2.0*shift
   END SUBROUTINE calculateF

   ! Calculate the function G(s) given in [3]
   ! Input:  s2        - s^2 where s is the reduced gradient
   !         Fs2       - F(s) s^2
   !         dFs2_ds   - first derivative of F(s)s^2 with respect to s
   !         d2Fs2_ds2 - second derivative of F(s)s^2
   !         Hs2       - H(s) s^2
   !         dHs2_ds   - first derivative of H(s)s^2 with respect to s
   !         d2Hs2_ds2 - second derivative of H(s)s^2
   ! Output: G         - value of the function G
   !         dGs2_ds   - first derivative of G(s)s^2 with respect to s
   !         d2Gs2_ds2 - second derivative of G(s)s^2
   SUBROUTINE calculateG(s2, Fs2, dFs2_ds, d2Fs2_ds2, Hs2, dHs2_ds, d2Hs2_ds2, G, dGs2_ds, d2Gs2_ds2)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: s2, Fs2, dFs2_ds, d2Fs2_ds2, Hs2, dHs2_ds, d2Hs2_ds2
      REAL, INTENT(OUT) :: G, dGs2_ds, d2Gs2_ds2

      ! helper variables
      REAL :: AHs2_1_2, AHs2_3_2, r1_Fs2, D_Hs2, D_Hs2Sqr, D_Hs2Cub, &
              D_Hs2_5_2, D_Hs2_7_2, D_Hs2_9_2, D_Hs2_11_2
      REAL :: part1, dpart1_dh, d2part1_dh2
      REAL :: arg1, arg2, exp_erfc, &
              part2, dpart2_dh, d2part2_dh2
      REAL :: alpha, Ebeta, r3Pi_4_alpha, beta_s2, Ebeta_s2, &
              dalpha_dh, d2alpha_dh2, dalpha_df, d2alpha_dfdh, dbeta_dh, d2beta_dh2

      ! parameters of the exchange hole given in [3]
      REAL, PARAMETER :: PI_75 = 2.356194490192344929, SQRT_PI = 1.77245385090551602729816748334, &
                         A = 1.0161144, sqrtA = 1.008025, r9_4A = 2.25/A, & !sqrt(A), 9/4A
                         B = -0.37170836, &
                         C = -0.077215461, &
                         D = 0.57786348, &
                         E = -0.051955731

      ! calculate the helper variables
      AHs2_1_2 = sqrtA*SQRT(Hs2)
      AHs2_3_2 = AHs2_1_2*A*Hs2
      r1_Fs2 = 1.0 + Fs2
      D_Hs2 = D + Hs2
      D_Hs2Sqr = D_Hs2*D_Hs2
      D_Hs2Cub = D_Hs2*D_Hs2Sqr
      D_Hs2_5_2 = D_Hs2Sqr*SQRT(D_Hs2)
      D_Hs2_7_2 = D_Hs2_5_2*D_Hs2
      D_Hs2_9_2 = D_Hs2_5_2*D_Hs2Sqr
      D_Hs2_11_2 = D_Hs2_5_2*D_Hs2Cub

      ! calculate first part of the term called 'a' in [3] eq. (A2) and its derivatives with respect to H(s)s^2 and F(s)s^2
      !
      !                          2      2           2 2          2 3
      !          __ 15E + 6C(1+Fs )(D+Hs ) + 4B(D+Hs )  + 8A(D+Hs )
      ! part1 = VPi ------------------------------------------------
      !                                    2  7/2
      !                        16 ( D + H s  )
      !
      part1 = SQRT_PI*(15.0*E + 6.0*C*r1_Fs2*D_Hs2 + 4.0*B*D_Hs2Sqr + 8.0*A*D_Hs2Cub)/(16.0*D_Hs2_7_2)
      !
      !                               2      2            2 2          2 3
      ! d part1     __ 105E + 30C(1+Fs )(D+Hs ) + 12B(D+Hs )  + 8A(D+Hs )
      ! ------- = -VPi ---------------------------------------------------
      ! d(Hs^2)                                  2  9/2
      !                              32 ( D + H s  )
      !
      dpart1_dh = -SQRT_PI*(105.0*E + 30.0*C*r1_Fs2*D_Hs2 + 12.0*B*D_Hs2Sqr + 8.0*A*D_Hs2Cub)/(32.0*D_Hs2_9_2)
      !
      !  2                             2      2            2 2           2 3
      ! d  part1    __ 945E + 210C(1+Fs )(D+Hs ) + 60B(D+Hs )  + 24A(D+Hs )
      ! -------- = VPi -----------------------------------------------------
      !        2                                   2  11/2
      ! d(Hs^2)                        64 ( D + H s  )
      !
      d2part1_dh2 = SQRT_PI*(945.0*E + 210.0*C*r1_Fs2*D_Hs2 + 60.0*B*D_Hs2Sqr + 24.0*A*D_Hs2Cub)/(64.0*D_Hs2_11_2)
      !
      ! d part1    __       3 C
      ! ------- = VPi -----------------
      ! d(Fs^2)                  2  5/2
      !               8 ( D + H s  )
      !
      dalpha_df = SQRT_PI*0.375*C/D_Hs2_5_2
      !
      !      2
      !     d part1         __       15 C
      ! --------------- = -VPi ------------------
      ! d(Fs^2) d(Hs^2)                    2  7/2
      !                        16 ( D + H s  )
      !
      d2alpha_dfdh = -2.5*dalpha_df/D_Hs2

      ! calculate second part of the term called 'a' in [3] eq. (A2) and its derivatives
      !                                            ________
      !                       /     2 \       /   /     2   \
      !         3 Pi  ___    | 9 H s   |     |   / 9 H s     |
      ! part2 = ---- V A  Exp| ------- | Erfc|  /  -------   |
      !          4           |   4 A   |     | V     4 A     |
      !                       \       /       \             /
      !
      arg1 = r9_4A*Hs2
      arg2 = SQRT(arg1)
      exp_erfc = EXP(arg1)*erfc(arg2)
      part2 = PI_75*sqrtA*exp_erfc
      !
      !                      /                                \
      ! d part2   3 Pi  ___ |  9                    3          |
      ! ------- = ---- V A  | --- exp_erfc - ----------------  |
      ! d(Hs^2)    4        | 4 A                          2   |
      !                      \               2 sqrt(Pi A Hs ) /
      !
      dpart2_dh = PI_75*sqrtA*(r9_4A*exp_erfc - 1.5/(SQRT_PI*AHs2_1_2))
      !
      !  2                    /                                 2    \
      ! d  part2   3 Pi  ___ |   81                  6 A - 27 Hs      |
      ! -------- = ---- V A  | ------ exp_erfc + -------------------- |
      !        2    4        | 16 A^2                           2 3/2 |
      ! d(Hs^2)               \                  8 sqrt(Pi)(A Hs )   /
      !
      d2part2_dh2 = PI_75*sqrtA*(r9_4A**2*exp_erfc + (6.0*A - 27.0*Hs2)/(8.0*SQRT_PI*AHs2_3_2))

      ! calculate the 'a' and 'b' terms ([3] eq. (A2) and (A3) )
      !
      ! a = part1 - part2
      !
      !             __  2
      !         15 VPi s
      ! b = -----------------
      !                2  7/2
      !     16 (D + H s  )
      !
      alpha = part1 - part2
      beta_s2 = 0.9375*SQRT_PI/D_Hs2_7_2
      Ebeta = E*beta_s2*s2

      ! the derivatives of a with respect to Fs^2 and Hs^2
      !
      !    da      d part1   d part2
      ! -------- = ------- - -------
      ! d (Hs^2)   d(Hs^2)   d(Hs^2)
      !
      dalpha_dh = dpart1_dh - dpart2_dh
      !
      !    2        2          2
      !   d a      d part1    d part2
      ! -------- = -------- - --------
      !        2          2          2
      ! d(Hs^2)    d(Hs^2)    d(Hs^2)
      !
      d2alpha_dh2 = d2part1_dh2 - d2part2_dh2
      !
      !    da      d part1
      ! -------- = -------
      ! d (Fs^2)   d(Fs^2)
      !
      !        2               2
      !       d a             d part1
      ! --------------- = ---------------
      ! d(Fs^2) d(Hs^2)   d(Fs^2) d(Hs^2)
      !
      ! (this expressions are already defined earlier)
      !
      ! and the derivatives of b/s^2 with respect to Hs^2
      !                        __
      ! d(b/s^2)          105 VPi          7  b / s^2
      ! -------- = - ----------------- = - -  --------
      ! d (Hs^2)                2  9/2     2         2
      !              32 (D + H s  )           D + H s
      !
      dbeta_dh = -3.5*beta_s2/D_Hs2
      !
      !  2                     __
      ! d (b/s^2)         945 VPi          9  d(b/s^2)
      ! --------- = ------------------ = - -  --------
      !         2              2  11/2     2  d (Hs^2)
      ! d (Hs^2)    64 (D + H s  )
      !
      d2beta_dh2 = -4.5*dbeta_dh/D_Hs2

      ! calculate the function G(s) and its derivatives
      !
      !       3Pi/4 + a
      ! G = - ---------
      !          E b
      !
      r3Pi_4_alpha = PI_75 + alpha
      Ebeta_s2 = E*beta_s2
      G = -r3Pi_4_alpha/Ebeta
      !
      !             /  /3 Pi    \  d(b/s^2)   /           d a    \  d(Hs^2)       d a   d(Fs^2)
      !            <  ( ---- + a ) -------   /  b/s^2  - -------  > -------  -  ------- -------
      ! d (Gs^2)    \  \ 4      /  d(Hs^2)  /            d(Hs^2) /    d s       d(Fs^2)   ds
      ! -------- = ----------------------------------------------------------------------------
      !    ds                                     E b/s^2
      !
      dGs2_ds = ((r3Pi_4_alpha*dbeta_dh/beta_s2 - dalpha_dh)*dHs2_ds - dalpha_df*dFs2_ds)/Ebeta_s2
      !
      !                                                     /3Pi    \  /             2\           /           \
      !                                                    ( --- + a )( b  h" + b   h' ) + 2b  h'( f'a  + h'a  ) + ...
      !                                            2        \ 4     /  \ h       hh   /      h    \   f      h/
      !  2         -f" a  - h" a  - 2 f' h' a   - h' a   + -----------------------------------------------------------
      ! d (Gs^2)        f       h            fh       hh                                b/s^2
      ! -------- = ---------------------------------------------------------------------------------------------------
      !   ds^2                                               E b/s^2
      !
      !            /3Pi    \   2  2
      !           ( --- + a ) b  h'
      !            \ 4     /   h
      ! ... = - 2 -----------------
      !                  b/s^2
      !
      ! where the primes indicate derivative with respect to s and the subscripts derivatives with respect to the subscript
      ! and f and h are abbrevations for f = Fs^2 and h = Hs^2
      !
      d2Gs2_ds2 = (-d2Fs2_ds2*dalpha_df - d2Hs2_ds2*dalpha_dh &
                   - 2.0*dFs2_ds*dHs2_ds*d2alpha_dfdh - dHs2_ds**2*d2alpha_dh2 &
                   + (r3Pi_4_alpha*(dbeta_dh*d2Hs2_ds2 + d2beta_dh2*dHs2_ds**2) &
                      + 2.0*dbeta_dh*dHs2_ds*(dFs2_ds*dalpha_df + dHs2_ds*dalpha_dh) &
                      - 2.0*r3Pi_4_alpha*(dbeta_dh*dHs2_ds)**2/beta_s2)/beta_s2)/Ebeta_s2
   END SUBROUTINE calculateG

   ! Calculate the function H(s) given in [3]
   ! Input:  s         - reduced gradient
   ! Output: H         - value of the function H
   !         dHs2_ds   - first derivative d(s^2*H(s))/ds
   !         d2Hs2_ds2 - second derivative d^2(s^2H(s))/ds^2
   SUBROUTINE calculateH(s, H, dHs2_ds, d2Hs2_ds2)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: s
      REAL, INTENT(OUT) :: H, dHs2_ds, d2Hs2_ds2

      ! helper variables
      REAL :: s2, s3, s4, s5, s6
      REAL :: numer, denom
      REAL :: dnum_ds, dden_ds
      REAL :: d2num_ds2, d2den_ds2

      ! parameters given in [3]
      REAL, PARAMETER :: &
         a1 = 0.00979681, &
         a2 = 0.0410834, &
         a3 = 0.187440, &
         a4 = 0.00120824, &
         a5 = 0.0347188

      ! calculate helper variables
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      s5 = s3*s2
      s6 = s4*s2

      ! calculate function H(s) with [3] eq. (A5)
      !
      !            2         4
      !        a  s   +  a  s
      !         1         2
      ! H = -------------------------
      !             4       5       6
      !     1 + a  s  + a  s  + a  s
      !          3       4       5
      !
      numer = a1*s2 + a2*s4
      denom = 1.0 + a3*s4 + a4*s5 + a5*s6
      H = numer/denom

      ! calculate the first derivative of s^2 H(s)
      !
      !     /    \
      !  d | f(x) |   f'(x) - f(x)g'(x) / g(x)
      ! -- | ---- | = ------------------------
      ! dx | g(x) |            g(x)
      !     \    /
      !
      numer = numer*s2
      dnum_ds = 4.0*s3*a1 + 6.0*s5*a2
      dden_ds = 4.0*s3*a3 + 5.0*s4*a4 + 6.0*s5*a5
      dHs2_ds = (dnum_ds - (dden_ds*numer)/denom)/denom

      ! calculate the second derivative of s^2 H(s)
      !
      !                                                               2
      !                         2 f'(x)g'(x) + f(x)g"(x) - 2 f(x)g'(x) /g(x)
      !    2  /    \    f"(x) - --------------------------------------------
      !   d  | f(x) |                              g(x)
      ! ---- | ---- | = ----------------------------------------------------
      ! dx^2 | g(x) |                          g(x)
      !       \    /
      !
      d2num_ds2 = 12.0*s2*a1 + 30.0*s4*a2
      d2den_ds2 = 12.0*s2*a3 + 20.0*s3*a4 + 30.0*s4*a5
      d2Hs2_ds2 = (d2num_ds2 - (2.0*dnum_ds*dden_ds + numer*d2den_ds2 - 2.0*numer*dden_ds**2/denom)/denom)/denom
   END SUBROUTINE calculateH

   ! Calculate several Fourier transformations
   ! a) the Fourier transformed potential
   !                                                 /          \
   !          /      | erf |       \     4 Pi       |   |k+G|^2  |
   ! V (k) = ( k + G | --- | k + G  ) = -------  Exp| - -------  |                             [1]
   !  G       \      |  r  |       /    |k+G|^2     |    4 w^2   |
   !                                                 \          /
   !
   ! b) muffin-tin basis function
   !                                                            R
   !                                     -L                      /
   !           /      |      \     4 Pi i         ^     -i G R  |     2
   ! MT (k) = ( k + G | M     ) = ---------- Y ( k+G ) e        | dr r  Phi(r) j ( |k+G| r )  [2]
   !   G,I     \      |  k,I /     Sqrt(Om)   LM                |               L
   !                                                           /
   !                                                            0
   !
   ! c) interstitial basis function
   !                                                      -----
   !          /      |      \                       4 Pi   \     -i(G - G_I)R_a  Sin(|G_I - G| R_a) - |G_I - G| R_a Cos(|G_I - G| R_a)
   ! IN    = ( k + G | M     ) = kronecker(G,G ) - ------   )   e               ------------------------------------------------------- [3]
   !   G,I    \      |  k,I /                 I      Om    /                                          |G_I - G|^3
   !                                                      -----
   !                                                        a
   !                                                                            \_________________________ ___________________________/
   !                                                                                                      V
   !                                                                                                     I_a
   !
   ! In the code:
   ! V_G(k):  potential
   ! MT_G(k): muffintin
   ! IN_G:    interstitial
   ! I_a:     inter_atom
   !
   ! Input:
   ! rmsh        - array of radial mesh points
   ! rmt         - radius of the muffin tin
   ! dx          - logarithmic increment of radial mesh
   ! jri         - number of radial mesh points
   ! jmtd        - dimension of radial mesh points
   ! bk          - one(!) k-vector for which FT is calculated
   ! bmat        - reciprocal lattice vector for all directions
   ! vol         - volume of the unit cell
   ! ntype       - number of atom types
   ! neq         - number of symmetry-equivalent atoms of atom type i
   ! natd        - dimesion of atoms in unit cell
   ! taual       - vector of atom positions (internal coordinates)
   ! lcutm       - l cutoff for mixed basis
   ! maxlcutm    - maximum of all these l cutoffs
   ! nindxm      - number of radial functions of mixed basis
   ! maxindxm    - maximum of these numbers
   ! gptm        - reciprocal lattice vectors of the mixed basis (internal coord.)
   ! ngptm       - number of vectors (for treated k-vector)
   ! pgptm       - pointer to the appropriate g-vector (for treated k-vector)
   ! gptmd       - dimension of gptm
   ! basm        - mixed basis functions (mt + inter) for treated k-vector
   ! lexp        - cutoff of spherical harmonics expansion of plane wave
   ! noGPts      - no g-vectors used for Fourier trafo
   ! Output:
   ! potential   - Fourier transformation of the potential
   ! muffintin   - Fourier transformation of all MT functions
   ! interstital - Fourier transformation of all IR functions
   SUBROUTINE calculate_fourier_transform( &
      ! Input
      rmsh, rmt, dx, jri, jmtd, bk, &
      bmat, vol, ntype, neq, natd, taual, lcutm, maxlcutm, &
      nindxm, maxindxm, gptm, ngptm, pgptm, gptmd, &
      basm, noGPts, irank, &
      ! Output
      potential, muffintin, interstitial)

      USE m_constants
      USE m_types_hybrid, ONLY: gptnorm
      USE m_util, ONLY: sphbessel, pure_intgrf, intgrf_init, intgrf_out, NEGATIVE_EXPONENT_WARNING, NEGATIVE_EXPONENT_ERROR

      IMPLICIT NONE

      ! scalar input
      INTEGER, INTENT(IN)    :: natd, ntype, maxlcutm
      INTEGER, INTENT(IN)    :: jmtd, irank
      INTEGER, INTENT(IN)    :: maxindxm
      INTEGER, INTENT(IN)    :: gptmd, noGPts
      REAL, INTENT(IN)       :: vol

      ! array input
      INTEGER, INTENT(IN)    :: lcutm(ntype)
      INTEGER, INTENT(IN)    :: nindxm(0:maxlcutm, ntype), neq(ntype)
      INTEGER, INTENT(IN)    :: jri(ntype)
      INTEGER, INTENT(IN)    :: gptm(3, gptmd)
      INTEGER, INTENT(IN)    :: ngptm
      INTEGER, INTENT(IN)    :: pgptm(ngptm)

      REAL, INTENT(IN)       :: bk(3)
      REAL, INTENT(IN)       :: rmsh(jmtd, ntype), rmt(ntype), dx(ntype)
      REAL, INTENT(IN)       :: basm(jmtd, maxindxm, 0:maxlcutm, ntype)
      REAL, INTENT(IN)       :: bmat(3, 3)!,amat(3,3)
      REAL, INTENT(IN)       :: taual(3, natd)

      ! array output
      REAL, INTENT(OUT)   :: potential(noGPts)                           ! Fourier transformed potential
      COMPLEX, INTENT(OUT)   :: muffintin(noGPts, maxindxm, &                 ! muffin-tin overlap integral
                                          (maxlcutm + 1)**2, ntype, MAXVAL(neq))
      COMPLEX, INTENT(OUT)   :: interstitial(noGPts, gptmd)                  ! interstistial overlap intergral

      ! private scalars
      INTEGER                :: cg, cg2, ci, cl, cn, cr                          ! counter variables
      REAL                   ::  r2Pi, r4Pi, pi_omega2                   !  2*Pi, 4*Pi, Pi/omega^2
      REAL                   :: sVol, r4Pi_sVol, r4Pi_Vol                   ! sqrt(vol), 4*Pi/sqrt(Vol), 4*Pi/Vol
      REAL                   :: omega, r1_omega2, r1_4omega2
!     REAL,    PARAMETER     :: omega = omega_VHSE()                        ! omega of the HSE functional
!     REAL,    PARAMETER     :: r1_omega2  = 1.0 / omega**2                 ! 1/omega^2
!     REAL,    PARAMETER     :: r1_4omega2 = 0.25 * r1_omega2               ! 1/(4 omega^2)
      COMPLEX, PARAMETER     :: img = (0.0, 1.0)                             ! i

      ! private arrays
      INTEGER                :: gPts(3, noGPts)                              ! g vectors (internal units)
      INTEGER                :: gPts_gptm(3, noGpts, gptmd)                   ! gPts - gptm
      INTEGER                :: natdPtr(ntype + 1)                            ! pointer to all atoms of one type
      REAL, ALLOCATABLE   :: gridf(:, :)                                  ! grid for radial integration
      REAL                   :: k_G(3, noGPts)                               ! k + G
      REAL                   :: AbsK_G(noGPts), AbsK_G2(noGPts)              ! abs(k+G), abs(k+G)^2
      REAL                   :: arg(noGPts)                                 ! abs(k+G)^2 / (4*omega^2)
      REAL                   :: sphbesK_Gr(noGPts, jmtd, 0:maxlcutm, ntype)    ! spherical bessel function of abs(k+G)r
      TYPE(intgrf_out)       :: intgrMT(noGPts, maxindxm, 0:maxlcutm, ntype)   ! integration in muffin-tin
      REAL                   :: abs_dg(noGpts, gptmd)                        ! abs(gPts - gptm)
      COMPLEX                :: imgl(0:maxlcutm)                            ! i^l
      COMPLEX                :: Ylm(noGPts, (maxlcutm + 1)**2)                 ! spherical harmonics for k+G and all lm
      COMPLEX                :: expIGR(noGPts, ntype, MAXVAL(neq))            ! exp(-iGR) for all atom types
      COMPLEX                :: sumInter_atom(noGpts, gptmd)                 ! sum over inter-atom factors

      ! Calculate helper variables
      r2Pi = 2.0*pi_const
      r4Pi = 4.0*pi_const
      sVol = SQRT(vol)
      r4Pi_sVol = r4Pi/sVol
      r4Pi_Vol = r4Pi/Vol
      omega = omega_hse!omega_VHSE()
      r1_omega2 = 1.0/omega**2
      r1_4omega2 = 0.25*r1_omega2
      pi_omega2 = pi_const*r1_omega2

      ! calculate pointers for all atom-types to all atoms
      natdPtr(ntype + 1) = natd
      DO cn = ntype, 1, -1
         natdPtr(cn) = natdPtr(cn + 1) - neq(cn)
      END DO

      ! Define imgl(l) = img**l
      imgl(0) = 1.0
      DO ci = 1, maxlcutm
         imgl(ci) = imgl(ci - 1)*img
      END DO

      ! generate grid for fast radial integration
      CALL intgrf_init(ntype, jmtd, jri, dx, rmsh, gridf)

      ! Set all arrays to 0
      gPts = 0
      k_G = 0
      AbsK_G = 0
      arg = 0
      sphbesK_Gr = 0
      intgrMT%value = 0
      intgrMT%ierror = 0
      Ylm = 0
      expIGR = 0
      gPts_gptm = 0
      abs_dg = 0
      sumInter_atom = 0
      potential = 0
      muffintin = 0
      interstitial = 0

! Calculate the muffin-tin basis function overlap integral
!                                                            R
!                                     -L                      /
!           /      |      \     4 Pi i         ^     -i G R  |     2
! MT (k) = ( k + G | M     ) = ---------- Y ( k+G ) e        | dr r  Phi(r) j ( |k+G| r )  [2]
!   G,I     \      |  k,I /     Sqrt(Om)   LM                |               L
!                                                           /
!                                                            0
      IF (ngptm < noGpts) STOP 'hsefunctional: error calculating Fourier coefficients, noGpts too large'

      gPts(:, :) = gptm(:, pgptm(1:noGPts))
#ifndef __PGI

      gpoints:FORALL (cg=1:noGPts)
      ntypesA:FORALL (cn=1:ntype)
      ! Calculate the phase factor exp(-iGR) for all atoms
      FORALL (ci=1:neq(cn))
      expIGR(cg, cn, ci) = EXP(-r2Pi*img*DOT_PRODUCT(taual(:, natdPtr(cn) + ci), gPts(:, cg)))
      muffintin(cg, :, :, cn, ci) = r4Pi_sVol*expIGR(cg, cn, ci)
      END FORALL

      ! Multiplication with i^-L
      FORALL (cl=0:lcutm(cn))
      muffintin(cg, :, cl*cl + 1:(cl + 1)*(cl + 1), cn, :) = muffintin(cg, :, cl*cl + 1:(cl + 1)*(cl + 1), cn, :)/imgl(cl)
      END FORALL

      END FORALL ntypesA

      ! Calculate the k+G, abs(k+G), and abs(k+G)^2
      k_G(:, cg) = MATMUL(gPts(:, cg) + bk(:), bmat(:, :))
      AbsK_G2(cg) = SUM(k_G(:, cg)**2)
      AbsK_G(cg) = SQRT(AbsK_G2(cg))

      ! Spherical harmonics are calculated for all lm's and k+G's
      Ylm(cg, :) = calcYlm(k_G(:, cg), maxlcutm)

      ! Perform the integration in eq.[2] for all muffin-tins
      ntypesB:FORALL (cn=1:ntype)

      ! Multiplication with the spherical harmonic
      FORALL (cl=1:(lcutm(cn) + 1)**2)
      muffintin(cg, :, cl, cn, :) = muffintin(cg, :, cl, cn, :)*Ylm(cg, cl)
      END FORALL

      ! Calculate the spherical bessel function
      FORALL (cr=1:jri(cn))
      sphbesK_Gr(cg, cr, :, cn) = calcSphBes(AbsK_G(cg)*rmsh(cr, cn), maxlcutm)
      END FORALL
      ! integrate the function and multiplicate to the result
      FORALL (cl=0:lcutm(cn))
      FORALL (ci=1:nindxm(cl, cn))
      intgrMT(cg, ci, cl, cn) = pure_intgrf(rmsh(:, cn)*basm(:, ci, cl, cn)* &
                                            sphbesK_Gr(cg, :, cl, cn), jri, jmtd, rmsh, dx, ntype, cn, gridf)
      muffintin(cg, ci, cl*cl + 1:(cl + 1)*(cl + 1), cn, :) = &
         muffintin(cg, ci, cl*cl + 1:(cl + 1)*(cl + 1), cn, :)*intgrMT(cg, ci, cl, cn)%value
      END FORALL
      END FORALL

      END FORALL ntypesB

! calculate the overlap with the interstitial basis function
!                                                      -----
!          /      |      \                       4 Pi   \     -i(G - G_I)R_a  Sin(|G_I - G| R_a) - |G_I - G| R_a Cos(|G_I - G| R_a)
! IN    = ( k + G | M     ) = kronecker(G,G ) - ------   )   e               ------------------------------------------------------- [3]
!   G,I    \      |  k,I /                 I      Om    /                                          |G_I - G|^3
!                                                      -----
!                                                        a
!                                                                            \_________________________ ___________________________/
!                                                                                                      V
!                                                                                                     I_a
      ! Calculate the difference of the G vectors and its absolute value
      FORALL (cg2=1:gptmd)

      gPts_gptm(:, cg, cg2) = gPts(:, cg) - gptm(:, cg2)
      abs_dg(cg, cg2) = gptnorm(gPts_gptm(:, cg, cg2), bmat)

      END FORALL

      END FORALL gpoints

      ! Check if any of the integrations failed and abort if one did
      IF (ANY(intgrMT%ierror == NEGATIVE_EXPONENT_ERROR)) THEN
         IF (irank == 0) WRITE (6, *) 'intgrf: Warning! Negative exponent x in extrapolation a+c*r**x'
      ELSEIF (ANY(intgrMT%ierror == NEGATIVE_EXPONENT_WARNING)) THEN
         IF (irank == 0) WRITE (6, *) 'intgrf: Negative exponent x in extrapolation a+c*r**x'
         STOP 'intgrf: Negative exponent x in extrapolation a+c*r**x'
      END IF

      ! Calculate the interstitial value with eq.[3] using the limit
      !           3
      !        R_a
      ! I_a = ------
      !         3
      ! if there's no difference of the two G vectors
      WHERE (abs_dg == 0)
      sumInter_atom = DOT_PRODUCT(rmt**3, neq)/3.0
      interstitial = 1.0 - r4Pi_Vol*sumInter_atom
      ELSEWHERE
      sumInter_atom = calculateSummation(abs_dg, gPts_gptm)
      interstitial = -r4Pi_Vol*sumInter_atom
      END WHERE

! Calculate the Fourier transformed potential
!                                                 /          \
!          /      | erf |       \      4 Pi      |   |k+G|^2  |
! V (k) = ( k + G | --- | k + G' ) = -------  Exp| - -------  | kronecker(G,G')             [1]
!  G       \      |  r  |       /    |k+G|^2     |    4 w^2   |
!                                                 \          /
      WHERE (absK_G /= 0)
      ! Calculate (k+G)^2 / 4*w^2
      arg = -r1_4omega2*AbsK_G2

      ! The potential is calculated using eq.[1]
      potential = r4Pi/AbsK_G2*EXP(arg)

! Calculate the Fourier transformed potential for the full(!) potential
!
!          /          | erfc |          \      Pi
! V (k) = ( k + G = 0 | ---- | k + G = 0 ) = -----
!  G       \          |   r  |          /     w^2
!
      ELSEWHERE
      ! the negative sign is added to compensate the sign when the
      ! contribution is added to the Coulomb matrix
      potential = -pi_omega2
      endwhere

      DEALLOCATE (gridf)
#else
      call judft_error("hsefunctional not implemented for PGI")
#endif
   CONTAINS

      ! Calculates the inter-atom parts and the summation over all atoms:
      ! -----
      !  \     -i(G - G_I)R_a  Sin(|G_I - G| R_a) - |G_I - G| R_a Cos(|G_I - G| R_a)
      !   )   e               -------------------------------------------------------
      !  /                                          |G_I - G|^3
      ! -----
      !   a
      ! Input:  abs_dg    - absolute value |G_I - G|
      !         gPts_gptm - vector G - G_I
      PURE FUNCTION calculateSummation(abs_dg, gPts_gptm)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: gPts_gptm(3, noGPts, gptmd)
         REAL, INTENT(IN)    :: abs_dg(noGPts, gptmd)
         COMPLEX             :: calculateSummation(noGPts, gptmd)
         INTEGER             :: cn, ci                                     ! counter variables
         REAL                :: abs_dgR(noGPts, gptmd, ntype)               ! abs(gPts - gptm)*R (R: radius MT)
         REAL                :: inter_atom(noGPts, gptmd, ntype)            ! inter-atom interaction for interstitial
         COMPLEX             :: expIdGR(noGPts, gptmd, ntype, MAXVAL(neq))   ! exp(-i(gPts-gptm)R)
         COMPLEX             :: sumExpIdGR(noGPts, gptmd, ntype)            ! sum over atom of same type

         atoms:FORALL (cn=1:ntype)
         !                                                         -i(G-G_I)R_a
         ! Calculate for all similar atoms (same type) the factor e
         ! use that the G vectors and atomic positions are given in internal units
         FORALL (ci=1:neq(cn))
         expIdGR(:, :, cn, ci) = EXP(-r2pi*img*my_dot_product(gPts_gptm(:, :, :), taual(:, natdPtr(cn) + ci)))
         END FORALL
         ! Add all factors for atoms of the same type
         sumExpIdGR(:, :, cn) = my_sum(expIdGR(:, :, cn, :))

         ! Calculate the inter-atom factor which is the same for all atoms of the same type
         abs_dgR(:, :, cn) = abs_dg*rmt(cn)
         inter_atom(:, :, cn) = (SIN(abs_dgR(:, :, cn)) - abs_dgR(:, :, cn)*COS(abs_dgR(:, :, cn)))/abs_dg**3
         END FORALL atoms

         ! Add the factors of all atoms together
         calculateSummation = my_dot_product(sumExpIdGR, inter_atom)
      END FUNCTION calculateSummation

   END SUBROUTINE calculate_fourier_transform

   ! Calculate several Fourier transformations
   ! a) the Fourier transformed potential
   !                                                 /          \
   !          /      | erf |       \     4 Pi       |   |k+G|^2  |
   ! V (k) = ( k + G | --- | k + G  ) = -------  Exp| - -------  |                             [1]
   !  G       \      |  r  |       /    |k+G|^2     |    4 w^2   |
   !                                                 \          /
   !
   ! b) muffin-tin basis function
   !                                                            R
   !                                     -L                      /
   !           /      |      \     4 Pi i         ^     -i G R  |     2
   ! MT (k) = ( k + G | M     ) = ---------- Y ( k+G ) e        | dr r  Phi(r) j ( |k+G| r )  [2]
   !   G,I     \      |  k,I /     Sqrt(Om)   LM                |               L
   !                                                           /
   !                                                            0
   !
   ! In the code:
   ! V_G(k):  potential
   ! MT_G(k): muffintin
   !
   ! Input:
   ! rmsh        - array of radial mesh points
   ! rmt         - radius of the muffin tin
   ! dx          - logarithmic increment of radial mesh
   ! jri         - number of radial mesh points
   ! jmtd        - dimension of radial mesh points
   ! bk          - one(!) k-vector for which FT is calculated
   ! bmat        - reciprocal lattice vector for all directions
   ! vol         - volume of the unit cell
   ! ntype       - number of atom types
   ! neq         - number of symmetry-equivalent atoms of atom type i
   ! natd        - dimesion of atoms in unit cell
   ! taual       - vector of atom positions (internal coordinates)
   ! lcutm       - l cutoff for mixed basis
   ! maxlcutm    - maximum of all these l cutoffs
   ! nindxm      - number of radial functions of mixed basis
   ! maxindxm    - maximum of these numbers
   ! gptm        - reciprocal lattice vectors of the mixed basis (internal coord.)
   ! ngptm       - number of vectors (for treated k-vector)
   ! pgptm       - pointer to the appropriate g-vector (for treated k-vector)
   ! gptmd       - dimension of gptm
   ! basm        - mixed basis functions (mt + inter) for treated k-vector
   ! lexp        - cutoff of spherical harmonics expansion of plane wave
   ! noGPts      - no g-vectors used for Fourier trafo
   ! Output:
   ! potential   - Fourier transformation of the potential
   ! muffintin   - Fourier transformation of all MT functions
   SUBROUTINE calculate_fourier_transform_once( &
      ! Input
      rmsh, rmt, dx, jri, jmtd, bk, ikpt, nkptf, &
      bmat, vol, ntype, neq, natd, taual, lcutm, maxlcutm, &
      nindxm, maxindxm, gptm, ngptm, pgptm, gptmd, &
      nbasp, basm, noGPts, invsat, invsatnr, irank, &
      ! Output
      potential, fourier_trafo)

      USE m_constants
      USE m_util, ONLY: sphbessel, pure_intgrf, intgrf_init, intgrf_out, NEGATIVE_EXPONENT_WARNING, NEGATIVE_EXPONENT_ERROR
      USE m_trafo, ONLY: symmetrize

      IMPLICIT NONE

      ! scalar input
      INTEGER, INTENT(IN)    :: natd, ntype, maxlcutm
      INTEGER, INTENT(IN)    :: ikpt, nkptf, jmtd
      INTEGER, INTENT(IN)    :: maxindxm, irank
      INTEGER, INTENT(IN)    :: gptmd, noGPts
      REAL, INTENT(IN)       :: vol

      ! array input
      INTEGER, INTENT(IN)    :: lcutm(ntype)
      INTEGER, INTENT(IN)    :: nindxm(0:maxlcutm, ntype), neq(ntype)
      INTEGER, INTENT(IN)    :: jri(ntype)
      INTEGER, INTENT(IN)    :: gptm(3, gptmd)
      INTEGER, INTENT(IN)    :: ngptm, nbasp
      INTEGER, INTENT(IN)    :: pgptm(ngptm)
      INTEGER, INTENT(IN)    :: invsat(natd), invsatnr(natd)

      REAL, INTENT(IN)       :: bk(3)
      REAL, INTENT(IN)       :: rmsh(jmtd, ntype), rmt(ntype), dx(ntype)
      REAL, INTENT(IN)       :: basm(jmtd, maxindxm, 0:maxlcutm, ntype)
      REAL, INTENT(IN)       :: bmat(3, 3)
      REAL, INTENT(IN)       :: taual(3, natd)

      ! array output
      REAL, INTENT(OUT)   :: potential(noGPts)                           ! Fourier transformed potential
#ifdef CPP_INVERSION
      REAL, INTENT(OUT)   :: fourier_trafo(nbasp, noGPts) !muffintin_out(nbasp,noGPts)
#else
      COMPLEX, INTENT(OUT)   :: fourier_trafo(nbasp, noGPts) !muffintin_out(nbasp,noGPts)
#endif

      ! private scalars
      INTEGER                :: cg, cg2, ci, cl, cm, cn, cr                       ! counter variables
      REAL                   :: r2Pi, r4Pi, pi_omega2                   ! Pi, 2*Pi, 4*Pi, Pi/omega^2
      REAL                   :: sVol, r4Pi_sVol, r4Pi_Vol                   ! sqrt(vol), 4*Pi/sqrt(Vol), 4*Pi/Vol
      REAL                   :: omega, r1_omega2, r1_4omega2
!     REAL,    PARAMETER     :: omega = omega_VHSE()                        ! omega of the HSE functional
!     REAL,    PARAMETER     :: r1_omega2  = 1.0 / omega**2                 ! 1/omega^2
!     REAL,    PARAMETER     :: r1_4omega2 = 0.25 * r1_omega2               ! 1/(4 omega^2)
      COMPLEX, PARAMETER     :: img = (0.0, 1.0)                             ! i

      ! private arrays
      INTEGER                :: gPts(3, noGPts)                              ! g vectors (internal units)
      INTEGER                :: natdPtr(ntype + 1)                            ! pointer to all atoms of one type
      INTEGER                :: ptrType(nbasp), ptrEq(nbasp), &               ! pointer from ibasp to corresponding
                                ptrLM(nbasp), ptrN(nbasp)                    ! type, atom, l and m, and n
      REAL, ALLOCATABLE   :: gridf(:, :)                                  ! grid for radial integration
      REAL                   :: k_G(3, noGPts)                               ! k + G
      REAL                   :: AbsK_G(noGPts), AbsK_G2(noGPts)              ! abs(k+G), abs(k+G)^2
      REAL                   :: arg(noGPts)                                 ! abs(k+G)^2 / (4*omega^2)
      REAL                   :: sphbesK_Gr(noGPts, jmtd, 0:maxlcutm, ntype)    ! spherical bessel function of abs(k+G)r
      TYPE(intgrf_out)       :: intgrMT(noGPts, maxindxm, 0:maxlcutm, ntype)   ! integration in muffin-tin
      COMPLEX                :: imgl(0:maxlcutm)                            ! i^l
      COMPLEX                :: Ylm(noGPts, (maxlcutm + 1)**2)                 ! spherical harmonics for k+G and all lm
      COMPLEX                :: expIGR(noGPts, ntype, MAXVAL(neq))            ! exp(-iGR) for all atom types
      COMPLEX                :: muffintin(noGPts, maxindxm, &                 ! muffin-tin overlap integral
                                          (maxlcutm + 1)**2, ntype, MAXVAL(neq))
#ifdef CPP_INVERSION
      COMPLEX                :: sym_muffintin(nbasp, noGPts)                 ! symmetrized muffin tin
#endif

      LOGICAL, SAVE          :: first_entry = .TRUE.                        ! allocate arrays in first entry

      ! allocate arrays in first entry reuse later
      IF (first_entry) THEN
         ALLOCATE (already_known(nkptf), &            ! stores which elements are known
                   known_potential(maxNoGPts, nkptf), &            ! stores the potential for all k-points
                   known_fourier_trafo(nbasp, maxNoGPts, nkptf))            ! stores the fourier transform of the mixed basis
         ! initialization
         already_known = .FALSE.
         known_potential = 0.0
         known_fourier_trafo = 0.0
         ! unset flag as arrays are allocated
         first_entry = .FALSE.
      ELSE
         ! check if size of arrays has changed and stop with error if they did
         IF (SIZE(already_known) /= nkptf) STOP 'hsefunctional: Array size changed!'
         IF (SIZE(known_fourier_trafo, 1) /= nbasp) STOP 'hsefunctional: Array size changed!'
      END IF

      ! if the current k-point was not calculated yet
      IF (.NOT. (already_known(ikpt))) THEN

         ! Calculate helper variables
         r2Pi = 2.0*pi_const
         r4Pi = 4.0*pi_const
         sVol = SQRT(vol)
         r4Pi_sVol = r4Pi/sVol
         r4Pi_Vol = r4Pi/Vol
         omega = omega_hse!omega_VHSE()
         r1_omega2 = 1.0/omega**2
         r1_4omega2 = 0.25*r1_omega2
         pi_omega2 = pi_const*r1_omega2

         ! calculate pointers for all atom-types to all atoms
         natdPtr(ntype + 1) = natd
         DO cn = ntype, 1, -1
            natdPtr(cn) = natdPtr(cn + 1) - neq(cn)
         END DO

         ! Define imgl(l) = img**l
         imgl(0) = 1.0
         DO ci = 1, maxlcutm
            imgl(ci) = imgl(ci - 1)*img
         END DO

         ! generate grid for fast radial integration
         CALL intgrf_init(ntype, jmtd, jri, dx, rmsh, gridf)

         ! Set all arrays to 0
         gPts = 0
         k_G = 0
         AbsK_G = 0
         arg = 0
         sphbesK_Gr = 0
         intgrMT%value = 0
         intgrMT%ierror = 0
         Ylm = 0
         expIGR = 0
         potential = 0
         muffintin = 0

! Calculate the muffin-tin basis function overlap integral
!                                                            R
!                                     -L                      /
!           /      |      \     4 Pi i         ^     -i G R  |     2
! MT (k) = ( k + G | M     ) = ---------- Y ( k+G ) e        | dr r  Phi(r) j ( |k+G| r )  [2]
!   G,I     \      |  k,I /     Sqrt(Om)   LM                |               L
!                                                           /
!                                                            0
         IF (ngptm < noGpts) STOP 'hsefunctional: error calculating Fourier coefficients, noGpts too large'

         gPts(:, :) = gptm(:, pgptm(1:noGPts))

         gpoints:FORALL (cg=1:noGPts)
         ntypesA:FORALL (cn=1:ntype)

         ! Calculate the phase factor exp(-iGR) for all atoms
         FORALL (ci=1:neq(cn))
         expIGR(cg, cn, ci) = EXP(-r2Pi*img*DOT_PRODUCT(taual(:, natdPtr(cn) + ci), gPts(:, cg)))
         muffintin(cg, :, :, cn, ci) = r4Pi_sVol*expIGR(cg, cn, ci)
         END FORALL

         ! Multiplication with i^-L
         FORALL (cl=0:lcutm(cn))
         muffintin(cg, :, cl*cl + 1:(cl + 1)*(cl + 1), cn, :) = muffintin(cg, :, cl*cl + 1:(cl + 1)*(cl + 1), cn, :)/imgl(cl)
         END FORALL

         END FORALL ntypesA

         ! Calculate the k+G, abs(k+G), and abs(k+G)^2
         k_G(:, cg) = MATMUL(gPts(:, cg) + bk(:), bmat(:, :))
         AbsK_G2(cg) = SUM(k_G(:, cg)**2)
         AbsK_G(cg) = SQRT(AbsK_G2(cg))

         ! Spherical harmonics are calculated for all lm's and k+G's
         Ylm(cg, :) = calcYlm(k_G(:, cg), maxlcutm)

         ! Perform the integration in eq.[2] for all muffin-tins
         ntypesB:FORALL (cn=1:ntype)

         ! Multiplication with the spherical harmonic
         FORALL (cl=1:(lcutm(cn) + 1)**2)
         muffintin(cg, :, cl, cn, :) = muffintin(cg, :, cl, cn, :)*Ylm(cg, cl)
         END FORALL

         ! Calculate the spherical bessel function
         FORALL (cr=1:jri(cn))
         sphbesK_Gr(cg, cr, :, cn) = calcSphBes(AbsK_G(cg)*rmsh(cr, cn), maxlcutm)
         END FORALL
         ! integrate the function and multiplicate to the result
         FORALL (cl=0:lcutm(cn))
         FORALL (ci=1:nindxm(cl, cn))
         intgrMT(cg, ci, cl, cn) = pure_intgrf(rmsh(:, cn)*basm(:, ci, cl, cn)* &
                                               sphbesK_Gr(cg, :, cl, cn), jri, jmtd, rmsh, dx, ntype, cn, gridf)
         muffintin(cg, ci, cl*cl + 1:(cl + 1)*(cl + 1), cn, :) = &
            muffintin(cg, ci, cl*cl + 1:(cl + 1)*(cl + 1), cn, :)*intgrMT(cg, ci, cl, cn)%value
         END FORALL
         END FORALL

         END FORALL ntypesB

         END FORALL gpoints

         ! Check if any of the integrations failed and abort if one did
         IF (ANY(intgrMT%ierror == NEGATIVE_EXPONENT_ERROR)) THEN
            IF (irank == 0) WRITE (6, *) 'intgrf: Warning! Negative exponent x in extrapolation a+c*r**x'
         ELSEIF (ANY(intgrMT%ierror == NEGATIVE_EXPONENT_WARNING)) THEN
            IF (irank == 0) WRITE (6, *) 'intgrf: Negative exponent x in extrapolation a+c*r**x'
            STOP 'intgrf: Negative exponent x in extrapolation a+c*r**x'
         END IF

! Calculate the Fourier transformed potential
!                                                 /          \
!          /      | erf |       \      4 Pi      |   |k+G|^2  |
! V (k) = ( k + G | --- | k + G' ) = -------  Exp| - -------  | kronecker(G,G')             [1]
!  G       \      |  r  |       /    |k+G|^2     |    4 w^2   |
!                                                 \          /
         WHERE (absK_G /= 0)
         ! Calculate (k+G)^2 / 4*w^2
         arg = -r1_4omega2*AbsK_G2

         ! The potential is calculated using eq.[1]
         potential = r4Pi/AbsK_G2*EXP(arg)

! Calculate the Fourier transformed potential for the full(!) potential
!
!          /          | erfc |          \      Pi
! V (k) = ( k + G = 0 | ---- | k + G = 0 ) = -----
!  G       \          |   r  |          /     w^2
!
         ELSEWHERE
         ! the negative sign is added to compensate the sign when the
         ! contribution is added to the Coulomb matrix
         potential = -pi_omega2
         endwhere

         DEALLOCATE (gridf)

         !
         ! Create pointer which correlate the position in the array with the
         ! appropriate indices of the MT mixed basis function
         !
         cg = 0
         DO cn = 1, ntype
            DO ci = 1, neq(cn)
               DO cl = 0, lcutm(cn)
                  DO cm = -cl, cl
                     DO cr = 1, nindxm(cl, cn)
                        cg = cg + 1
                        ptrType(cg) = cn
                        ptrEq(cg) = ci
                        ptrLM(cg) = (cl + 1)*cl + cm + 1
                        ptrN(cg) = cr
                     END DO
                  END DO
               END DO
            END DO
         END DO
         IF (nbasp /= cg) STOP 'hsefunctional: wrong array size: nbasp'

#ifdef CPP_INVERSION
         ! Symmetrize muffin tin fourier transform
         DO ci = 1, nbasp
            sym_muffintin(ci, :noGPts) = muffintin(:, ptrN(ci), ptrLM(ci), ptrType(ci), ptrEq(ci))
         END DO
         DO cg = 1, noGPts
            CALL symmetrize(sym_muffintin(:, cg), 1, nbasp, 2, .FALSE., &
                            ntype, ntype, neq, lcutm, maxlcutm, &
                            nindxm, natd, invsat, invsatnr)
         END DO
         ! store the fourier transform of the muffin tin basis
         known_fourier_trafo(:, :, ikpt) = REAL(sym_muffintin)
#else
         ! store the fourier transform of the muffin tin basis
         DO ci = 1, nbasp
            known_fourier_trafo(ci, :noGPts, ikpt) = CONJG(muffintin(:, ptrN(ci), ptrLM(ci), ptrType(ci), ptrEq(ci)))
         END DO
#endif
         ! store the fourier transform of the potential
         known_potential(:noGPts, ikpt) = potential
         ! set the flag so that the fourier transform is not calculated again
         already_known(ikpt) = .TRUE.
         IF (MINVAL(ABS(potential)) > 1e-12) THEN
            WRITE (*, *) 'hsefunctional: Warning! Smallest potential bigger than numerical precision!', MINVAL(ABS(potential))
            WRITE (*, *) 'Perhaps you should increase the number of g-Points used and recompile!'
         END IF

      END IF

      ! return the fourier transform
      potential = known_potential(:noGPts, ikpt)
      fourier_trafo = known_fourier_trafo(:, :noGPts, ikpt)

   END SUBROUTINE calculate_fourier_transform_once

   ! Correct the pure Coulomb Matrix with by subtracting the long-range component
   !
   !  /     |      |      \     /     |       |      \     /     |     |      \
   ! ( M    | V    | M     ) = ( M    | V     | M     ) - ( M    | V   | M     )
   !  \ k,I |  HSE |  k,J /     \ k,I |  coul |  k,J /     \ k,I |  LR |  k,J /
   !
   ! The long-range component is given py the potential
   !
   !         erf( w r )
   ! V (r) = ----------
   !  LR         r
   !
   ! Input:
   ! rmsh     - array of radial mesh points
   ! rmt      - radius of the muffin tin
   ! dx       - logarithmic increment of radial mesh
   ! jri      - number of radial mesh points
   ! jmtd     - dimension of radial mesh points
   ! nkptf    - number of total k-points
   ! nkptd    - dimension of k-points
   ! nkpti    - number of irreducible k-points in window (in KPTS)
   ! bk       - k-vector for all irreduble k-points
   ! bmat     - reciprocal lattice vector for all directions
   ! vol      - volume of unit cell
   ! ntype    - number of atom types
   ! neq      - number of symmetry-equivalent atoms of atom type i
   ! natd     - dimension of atoms in unit cell
   ! taual    - vector of atom positions (internal coordinates)
   ! lcutm    - l cutoff for mixed basis
   ! maxlcutm - maximum of all these l cutoffs
   ! nindxm   - number of radial functions of mixed basis
   ! maxindxm - maximum of these numbers
   ! gptm     - reciprocal lattice vectors of the mixed basis (internal coord.)
   ! ngptm    - number of vectors
   ! pgptm    - pointer to the appropriate g-vector
   ! gptmd    - dimension of gptm
   ! basm     - radial mixed basis functions (mt + inter)
   ! lexp     - cutoff of spherical harmonics expansion of plane wave
   ! maxbasm  - maximum number of mixed basis functions
   ! nbasm    - number of mixed basis function
   ! invsat   - number of inversion-symmetric atom
   ! invsatnr - number of inversion-symmetric atom
   ! Inout:
   ! coulomb  - Coulomb matrix which is changed
   SUBROUTINE change_coulombmatrix( &
      ! Input
      rmsh, rmt, dx, jri, jmtd, nkptf, nkptd, nkpti, bk, &
      bmat, vol, ntype, neq, natd, taual, lcutm, maxlcutm, &
      nindxm, maxindxm, gptm, ngptm, pgptm, gptmd, &
      basm, lexp, maxbasm, nbasm, invsat, invsatnr, irank, &
      ! Input & output
      coulomb)

      USE m_trafo, ONLY: symmetrize
      USE m_wrapper, ONLY: packmat, unpackmat, diagonalize, inverse
      USE m_olap, ONLY: olap_pw

      IMPLICIT NONE

      ! scalar input
      INTEGER, INTENT(IN)    :: natd, ntype, maxlcutm, lexp
      INTEGER, INTENT(IN)    :: jmtd, nkptf, nkptd, nkpti
      INTEGER, INTENT(IN)    :: maxindxm
      INTEGER, INTENT(IN)    :: gptmd, irank
      INTEGER, INTENT(IN)    :: maxbasm
      REAL, INTENT(IN)       :: vol

      ! array input
      INTEGER, INTENT(IN)    :: lcutm(ntype)
      INTEGER, INTENT(IN)    :: nindxm(0:maxlcutm, ntype), neq(ntype)
      INTEGER, INTENT(IN)    :: jri(ntype)
      INTEGER, INTENT(IN)    :: gptm(3, gptmd)
      INTEGER, INTENT(IN)    :: ngptm(nkptf)
      INTEGER, INTENT(IN)    :: pgptm(MAXVAL(ngptm), nkptf)
      INTEGER, INTENT(IN)    :: nbasm(nkptf)
      INTEGER, INTENT(IN)    :: invsat(natd), invsatnr(natd)

      REAL, INTENT(IN)       :: bk(3, nkptd)
      REAL, INTENT(IN)       :: rmsh(jmtd, ntype), rmt(ntype), dx(ntype)
      REAL, INTENT(IN)       :: basm(jmtd, maxindxm, 0:maxlcutm, ntype)
      REAL, INTENT(IN)       :: bmat(3, 3)
      REAL, INTENT(IN)       :: taual(3, natd)

      ! array inout
#ifdef CPP_INVERSION
      REAL, INTENT(INOUT)    :: coulomb(maxbasm*(maxbasm + 1)/2, nkpti)
#else
      COMPLEX, INTENT(INOUT) :: coulomb(maxbasm*(maxbasm + 1)/2, nkpti)
#endif

      ! private scalars
      INTEGER                :: ikpt, itype, ieq, il, im, iindxm, idum, n1, n2, ok    ! counters and other helper variables
      INTEGER                :: noGPts, nbasp                                 ! no used g-vectors, no MT functions

      ! private arrays
      INTEGER, ALLOCATABLE   :: ptrType(:), ptrEq(:), ptrL(:), ptrM(:), ptrN(:)  ! Pointer
      REAL                   :: potential(maxNoGPts)                         ! Fourier transformed potential
      COMPLEX                :: muffintin(maxNoGPts, maxindxm, &               ! muffin-tin overlap integral
                                          (maxlcutm + 1)**2, ntype, MAXVAL(neq))
      COMPLEX                :: interstitial(maxNoGPts, gptmd)                ! interstistial overlap intergral
      COMPLEX, ALLOCATABLE   :: coulmat(:, :)                                 ! helper array to symmetrize coulomb

      ! Check size of arrays
      IF (nkpti > nkptd) STOP 'hsefunctional: missmatch in dimension of arrays'
      nbasp = maxbasm - MAXVAL(ngptm)
      IF (ANY(nbasm - ngptm /= nbasp)) STOP 'hsefunctional: wrong assignment of nbasp'

      !
      ! Create pointer which correlate the position in the array with the
      ! appropriate indices of the MT mixed basis function
      !
      ALLOCATE (ptrType(nbasp), ptrEq(nbasp), ptrL(nbasp), ptrM(nbasp), ptrN(nbasp))
      nbasp = 0
      DO itype = 1, ntype
         DO ieq = 1, neq(itype)
            DO il = 0, lcutm(itype)
               DO im = -il, il
                  DO iindxm = 1, nindxm(il, itype)
                     nbasp = nbasp + 1
                     ptrType(nbasp) = itype
                     ptrEq(nbasp) = ieq
                     ptrL(nbasp) = il
                     ptrM(nbasp) = im
                     ptrN(nbasp) = iindxm
                  END DO
               END DO
            END DO
         END DO
      END DO

      !
      ! Change the Coulomb matrix for all k-points
      !
      DO ikpt = 1, nkpti
         ! use the same g-vectors as in the mixed basis
         ! adjust the limit of the array if necessary
         IF (ngptm(ikpt) < maxNoGPts) THEN
            noGPts = ngptm(ikpt)
         ELSE
            noGPts = maxNoGPts
         END IF

         !
         ! Calculate the Fourier transform of the mixed basis and the potential
         !
         CALL calculate_fourier_transform( &
            ! Input
            rmsh, rmt, dx, jri, jmtd, bk(:, ikpt), &
            bmat, vol, ntype, neq, natd, taual, lcutm, maxlcutm, &
            nindxm, maxindxm, gptm, ngptm(ikpt), pgptm(:, ikpt), gptmd, &
            basm, noGPts, irank, &
            ! Output
            potential, muffintin, interstitial)
         interstitial = CONJG(interstitial)

         ! Helper matrix for temporary storage of the attenuated Coulomb matrix
         ALLOCATE (coulmat(nbasm(ikpt), nbasm(ikpt)), stat=ok)
         IF (ok /= 0) STOP 'hsefunctional: failure at matrix allocation'
         coulmat = 0
         !
         ! Calculate the difference of the Coulomb matrix by the attenuation
         !
         DO n1 = 1, nbasp
            DO n2 = 1, n1
               ! muffin tin - muffin tin contribution
               coulmat(n2, n1) = -gPtsSummation(noGpts, &
                                                muffintin(:, ptrN(n2), (ptrL(n2) + 1)*ptrL(n2) + ptrM(n2) + 1, ptrType(n2), ptrEq(n2)), potential, &
                                                CONJG(muffintin(:, ptrN(n1), (ptrL(n1) + 1)*ptrL(n1) + ptrM(n1) + 1, ptrType(n1), ptrEq(n1))))
               coulmat(n1, n2) = CONJG(coulmat(n2, n1))
            END DO
         END DO
         DO n1 = nbasp + 1, nbasm(ikpt)
            DO n2 = 1, n1
               IF (n2 <= nbasp) THEN
                  ! muffin tin - interstitial contribution
                  coulmat(n2, n1) = -gPtsSummation(noGPts, &
                                                   muffintin(:, ptrN(n2), (ptrL(n2) + 1)*ptrL(n2) + ptrM(n2) + 1, ptrType(n2), ptrEq(n2)), &
                                                   potential, CONJG(interstitial(:, pgptm(n1 - nbasp, ikpt))))
                  coulmat(n1, n2) = CONJG(coulmat(n2, n1))
               ELSE
                  ! interstitial - interstitial contribution
                  coulmat(n2, n1) = -gPtsSummation(noGPts, &
                                                   interstitial(:, pgptm(n2 - nbasp, ikpt)), potential, &
                                                   CONJG(interstitial(:, pgptm(n1 - nbasp, ikpt))))
                  coulmat(n1, n2) = CONJG(coulmat(n2, n1))
               END IF
            END DO
         END DO

#ifdef CPP_INVERSION
         ! symmetrize matrix if system has inversion symmetry
         CALL symmetrize(coulmat, nbasm(ikpt), nbasm(ikpt), 3, .FALSE., &
                         ntype, ntype, neq, lcutm, maxlcutm, &
                         nindxm, natd, invsat, invsatnr)
#endif
         ! add the changes to the Coulomb matrix
         coulomb(:nbasm(ikpt)*(nbasm(ikpt) + 1)/2, ikpt) = packmat(coulmat) + coulomb(:nbasm(ikpt)*(nbasm(ikpt) + 1)/2, ikpt)
         DEALLOCATE (coulmat)

      END DO

      DEALLOCATE (ptrType, ptrEq, ptrL, ptrM, ptrN)

   END SUBROUTINE change_coulombmatrix

   ! Correct the pure Coulomb Matrix with by subtracting the long-range component
   ! during execution time
   !
   !  /     |      |      \     /     |       |      \     /     |     |      \
   ! ( M    | V    | M     ) = ( M    | V     | M     ) - ( M    | V   | M     )
   !  \ k,I |  HSE |  k,J /     \ k,I |  coul |  k,J /     \ k,I |  LR |  k,J /
   !
   ! The long-range component is given py the potential
   !
   !         erf( w r )
   ! V (r) = ----------
   !  LR         r
   !
   ! Input:
   ! rmsh     - array of radial mesh points
   ! rmt      - radius of the muffin tin
   ! dx       - logarithmic increment of radial mesh
   ! jri      - number of radial mesh points
   ! jmtd     - dimension of radial mesh points
   ! bk       - k-vector for this k-point
   ! ikpt     - number of this k-point
   ! nkptf    - number of total k-points
   ! bmat     - reciprocal lattice vector for all directions
   ! vol      - volume of unit cell
   ! ntype    - number of atom types
   ! neq      - number of symmetry-equivalent atoms of atom type i
   ! natd     - dimension of atoms in unit cell
   ! taual    - vector of atom positions (internal coordinates)
   ! lcutm    - l cutoff for mixed basis
   ! maxlcutm - maximum of all these l cutoffs
   ! nindxm   - number of radial functions of mixed basis
   ! maxindxm - maximum of these numbers
   ! gptm     - reciprocal lattice vectors of the mixed basis (internal coord.)
   ! ngptm    - number of vectors
   ! pgptm    - pointer to the appropriate g-vector
   ! gptmd    - dimension of gptm
   ! basm     - radial mixed basis functions (mt + inter)
   ! nbasm    - number of mixed basis function
   ! nobd     - dimension of occupied bands
   ! nbands   - number of bands
   ! nsst     - size of indx
   ! indx     - pointer to bands
   ! invsat   - number of inversion-symmetric atom
   ! invsatnr - number of inversion-symmetric atom
   ! cprod    - scalar product of mixed basis and wavefunction product basis
   ! wl_iks   -
   ! n_q      -
   ! Return:
   ! Change of the Coulomb matrix
   FUNCTION dynamic_hse_adjustment( &
      rmsh, rmt, dx, jri, jmtd, bk, ikpt, nkptf, bmat, vol, &
      ntype, neq, natd, taual, lcutm, maxlcutm, nindxm, maxindxm, &
      gptm, ngptm, pgptm, gptmd, basm, nbasm, &
      nobd, nbands, nsst, ibando, psize, indx, invsat, invsatnr, irank, &
      cprod_r, cprod_c, l_real, wl_iks, n_q)

      USE m_trafo, ONLY: symmetrize
      USE m_olap, ONLY: olap_pw, olap_pwp
      USE m_wrapper, ONLY: diagonalize, dotprod, matvec, packmat, inverse

      IMPLICIT NONE

      ! scalar input
      INTEGER, INTENT(IN)  :: natd, ntype, maxlcutm
      INTEGER, INTENT(IN)  :: jmtd, ikpt, nkptf
      INTEGER, INTENT(IN)  :: maxindxm
      INTEGER, INTENT(IN)  :: gptmd, irank
      INTEGER, INTENT(IN)  :: nbasm, nobd, nbands, ibando, psize
      INTEGER, INTENT(IN)  :: n_q
      REAL, INTENT(IN)     :: vol

      ! array input
      INTEGER, INTENT(IN)  :: lcutm(ntype)
      INTEGER, INTENT(IN)  :: nindxm(0:maxlcutm, ntype), neq(ntype)
      INTEGER, INTENT(IN)  :: jri(ntype)
      INTEGER, INTENT(IN)  :: gptm(3, gptmd)
      INTEGER, INTENT(IN)  :: ngptm
      INTEGER, INTENT(IN)  :: pgptm(ngptm)
      INTEGER, INTENT(IN)  :: nsst(nbands), indx(nbands, nbands)
      INTEGER, INTENT(IN)  :: invsat(natd), invsatnr(natd)
      REAL, INTENT(IN)     :: bk(3)
      REAL, INTENT(IN)     :: rmsh(jmtd, ntype), rmt(ntype), dx(ntype)
      REAL, INTENT(IN)     :: basm(jmtd, maxindxm, 0:maxlcutm, ntype)
      REAL, INTENT(IN)     :: bmat(3, 3)
      REAL, INTENT(IN)     :: taual(3, natd)
      REAL, INTENT(IN)     :: wl_iks(nobd)
      REAL, INTENT(IN)     :: cprod_r(nbasm, psize, nbands)
      COMPLEX, INTENT(IN)  :: cprod_c(nbasm, psize, nbands)
      LOGICAL, INTENT(IN)   :: l_real

      ! return type definition
      COMPLEX              :: dynamic_hse_adjustment(nbands, nbands)

      ! private scalars
      INTEGER              :: noGpts, nbasp
      INTEGER              :: itype, ieq, il, im, iindxm, ibasp
      INTEGER              :: igpt, iobd, iobd0, isst, iband1, iband2
      COMPLEX              :: cdum, gPtsSum

      ! private arrays
      INTEGER, ALLOCATABLE :: ptrType(:), ptrEq(:), ptrL(:), ptrM(:), ptrN(:) ! pointer to reference muffintin-array
      REAL                 :: potential(maxNoGPts)                        ! Fourier transformed potential
      COMPLEX              :: muffintin(maxNoGPts, maxindxm, &              ! muffin-tin overlap integral
                                        (maxlcutm + 1)**2, ntype, MAXVAL(neq))
      COMPLEX              :: interstitial(maxNoGPts, gptmd)               ! interstistial overlap intergral
#ifdef CPP_INVERSION
      REAL                 :: fourier_trafo(nbasm - ngptm, maxNoGPts)        ! Fourier trafo of all mixed basis functions
#else
      COMPLEX              :: fourier_trafo(nbasm - ngptm, maxNoGPts)        ! Fourier trafo of all mixed basis functions
#endif

      REAL                 :: cprod_fourier_trafo_r(maxNoGpts, psize, nbands)  ! Product of cprod and Fourier trafo
      COMPLEX              :: cprod_fourier_trafo_c(maxNoGpts, psize, nbands)  ! Product of cprod and Fourier trafo

      ! Initialisation
      dynamic_hse_adjustment = 0.0
      potential = 0.0
      fourier_trafo = 0.0
      cprod_fourier_trafo_r = 0.0
      cprod_fourier_trafo_c = 0.0
      noGPts = MIN(ngptm, maxNoGPts)
      nbasp = nbasm - ngptm

      !
      ! Calculate the fourier transform of the mixed basis once
      ! If it was already calculated load the old results
      !
      CALL calculate_fourier_transform_once( &
         rmsh, rmt, dx, jri, jmtd, bk, ikpt, nkptf, &
         bmat, vol, ntype, neq, natd, taual, lcutm, maxlcutm, &
         nindxm, maxindxm, gptm, ngptm, pgptm, gptmd, &
         nbasp, basm, noGPts, invsat, invsatnr, irank, &
         potential, fourier_trafo)

      ! Calculate the Fourier transform of the 'normal' basis
      ! by summing over the mixed basis
      DO igpt = 1, noGPts
         DO iobd0 = 1, psize!ibando,min(ibando+psize,nobd)
            iobd = iobd0 + ibando - 1
            IF (iobd > nobd) CYCLE
            DO iband1 = 1, nbands
               if (l_real) THEN
                  cprod_fourier_trafo_r(igpt, iobd0, iband1) = &
                     ! muffin tin contribution
                     dotprod(fourier_trafo(:nbasp, igpt), cprod_r(:nbasp, iobd0, iband1)) &
                     ! interstitial contribution (interstitial is kronecker_G,G')
                     + cprod_r(nbasp + igpt, iobd0, iband1)
               else
                  cprod_fourier_trafo_c(igpt, iobd0, iband1) = &
                     ! muffin tin contribution
                     dotprod(cprod_c(:nbasp, iobd0, iband1), fourier_trafo(:nbasp, igpt)) &
                     ! interstitial contribution (interstitial is kronecker_G,G')
                     + CONJG(cprod_c(nbasp + igpt, iobd0, iband1))
               endif
            END DO
         END DO
      END DO

      ! Summation over all G-vectors
      DO iband1 = 1, nbands
         DO iobd0 = 1, psize!ibando,min(ibando+psize,nobd)
            iobd = iobd0 + ibando - 1
            IF (iobd > nobd) CYCLE
            ! prefactor for all elements
            cdum = wl_iks(iobd)/n_q
            DO isst = 1, nsst(iband1)
               iband2 = indx(isst, iband1)
               ! do summation over G-vectors
               if (l_real) THEN
                  gPtsSum = cdum*gPtsSummation(noGPts, &
                                               cprod_fourier_trafo_r(:, iobd0, iband1), &
                                               potential, cprod_fourier_trafo_r(:, iobd0, iband2))
               ELSE
                  gPtsSum = cdum*gPtsSummation(noGPts, &
                                               conjg(cprod_fourier_trafo_c(:, iobd0, iband1)), &
                                               potential, cprod_fourier_trafo_c(:, iobd0, iband2))
               endif
               dynamic_hse_adjustment(iband2, iband1) = &
                  dynamic_hse_adjustment(iband2, iband1) - gPtsSum
            END DO
         END DO
      END DO

   END FUNCTION dynamic_hse_adjustment

   ! Helper function needed to use forall statements instead of do loop's
   ! calls the harmonicsr subroutine from 'util.F'
   PURE FUNCTION calcYlm(rvec, ll)
      USE m_util, ONLY: harmonicsr
      IMPLICIT NONE
      REAL, INTENT(IN)    :: rvec(3)
      INTEGER, INTENT(IN) :: ll
      COMPLEX             :: calcYlm((ll + 1)**2)
      CALL harmonicsr(calcYlm, rvec, ll)
   END FUNCTION calcYlm

   ! Helper function needed to use forall statements instead of do loop's
   ! calls the sphbessel subroutine from 'util.F'
   PURE FUNCTION calcSphBes(x, l)
      USE m_util, ONLY: sphbessel
      IMPLICIT NONE
      REAL, INTENT(IN)    :: x
      INTEGER, INTENT(IN) :: l
      REAL                :: calcSphBes(0:l)
      CALL sphbessel(calcSphBes, x, l)
   END FUNCTION calcSphBes

   ! Dot_product definition for two arrays (3d x 1d) with common size of first dimension
   PURE FUNCTION my_dot_product3x1(first, second)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: first(:, :, :)
      REAL, INTENT(IN)    :: second(SIZE(first, 1))
      REAL                :: my_dot_product3x1(SIZE(first, 2), SIZE(first, 3))
      INTEGER             :: ci, cj

      FORALL (ci=1:SIZE(first, 2), cj=1:SIZE(first, 3))
      my_dot_product3x1(ci, cj) = DOT_PRODUCT(first(:, ci, cj), second)
      END FORALL
   END FUNCTION my_dot_product3x1

   ! Dot_product definition for two arrays (4d x 1d) with common size of first dimension
   PURE FUNCTION my_dot_product4x1(first, second)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: first(:, :, :, :)
      REAL, INTENT(IN)    :: second(SIZE(first, 1))
      REAL                :: my_dot_product4x1(SIZE(first, 2), SIZE(first, 3), SIZE(first, 4))
      INTEGER             :: ci, cj, ck

      FORALL (ci=1:SIZE(first, 2), cj=1:SIZE(first, 3), ck=1:SIZE(first, 4))
      my_dot_product4x1(ci, cj, ck) = DOT_PRODUCT(first(:, ci, cj, ck), second)
      END FORALL
   END FUNCTION my_dot_product4x1

   ! Dot_product definition for two arrays (3d x 3d) with common size of all dimensions,
   ! The scalar product is evaluated over the last index
   PURE FUNCTION my_dot_product3x3(first, second)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: first(:, :, :)
      REAL, INTENT(IN) :: second(SIZE(first, 1), SIZE(first, 2), SIZE(first, 3))
      COMPLEX             :: my_dot_product3x3(SIZE(first, 1), SIZE(first, 2))
      INTEGER             :: ci, cj

      FORALL (ci=1:SIZE(first, 1), cj=1:SIZE(first, 2))
      my_dot_product3x3(ci, cj) = DOT_PRODUCT(first(ci, cj, :), second(ci, cj, :))
      END FORALL
   END FUNCTION my_dot_product3x3

   ! Dot_product definition for two arrays (4d x 4d) with common size of all dimensions,
   ! The scalar product is evaluated over the last index
   PURE FUNCTION my_dot_product4x4(first, second)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: first(:, :, :, :)
      REAL, INTENT(IN) :: second(SIZE(first, 1), SIZE(first, 2), SIZE(first, 3), SIZE(first, 4))
      COMPLEX             :: my_dot_product4x4(SIZE(first, 1), SIZE(first, 2), SIZE(first, 3))
      INTEGER             :: ci, cj, ck

      FORALL (ci=1:SIZE(first, 1), cj=1:SIZE(first, 2), ck=1:SIZE(first, 3))
      my_dot_product4x4(ci, cj, ck) = DOT_PRODUCT(first(ci, cj, ck, :), second(ci, cj, ck, :))
      END FORALL
   END FUNCTION my_dot_product4x4

   ! Sum over last index of 4d array
   PURE FUNCTION my_sum3d(array)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: array(:, :, :)
      COMPLEX             :: my_sum3d(SIZE(array, 1), SIZE(array, 2))
      INTEGER             :: ci, cj

      FORALL (ci=1:SIZE(array, 1), cj=1:SIZE(array, 2))
      my_sum3d(ci, cj) = SUM(array(ci, cj, :))
      END FORALL
   END FUNCTION my_sum3d

   ! Sum over last index of 4d array
   PURE FUNCTION my_sum4d(array)
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: array(:, :, :, :)
      COMPLEX             :: my_sum4d(SIZE(array, 1), SIZE(array, 2), SIZE(array, 3))
      INTEGER             :: ci, cj, ck

      FORALL (ci=1:SIZE(array, 1), cj=1:SIZE(array, 2), ck=1:SIZE(array, 3))
      my_sum4d(ci, cj, ck) = SUM(array(ci, cj, ck, :))
      END FORALL
   END FUNCTION my_sum4d

   ! Sum product of three vectors
   COMPLEX PURE FUNCTION crc_gPtsSummation(noGpts, vec1, vec2, vec3)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: noGPts
      COMPLEX, INTENT(IN) :: vec1(noGPts)
      REAL, INTENT(IN)    :: vec2(noGPts)
      COMPLEX, INTENT(IN) :: vec3(noGPts)
      COMPLEX             :: temp(noGPts)
      temp = vec1*vec3
      crc_gPtsSummation = DOT_PRODUCT(temp, vec2)
   END FUNCTION crc_gPtsSummation
   REAL PURE FUNCTION rrr_gPtsSummation(noGpts, vec1, vec2, vec3)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: noGPts
      REAL, INTENT(IN)    :: vec1(noGPts)
      REAL, INTENT(IN)    :: vec2(noGPts)
      REAL, INTENT(IN)    :: vec3(noGPts)
      REAL                :: temp(noGPts)
      temp = vec1*vec3
      rrr_gPtsSummation = DOT_PRODUCT(temp, vec2)
   END FUNCTION rrr_gPtsSummation

   ! Include the core valence exchange for the HSE functional
   ! The HSE potential is the Coulomb one attenuated with the complementary
   ! error function which can be expanded in Legendre polynomials (r' > r)
   !
   !                      -----                l+2n
   ! erfc( w |r - r'| )    \                  r           /         \
   ! ------------------ =   )    d  (w r')  ---------  P ( cos(r,r') )
   !    | r - r' |         /      ln          l+2n+1    l \         /
   !                      -----              r'
   !                       l,n
   !
   ! Then an analytical integration of the angular part of the exchange
   ! potential is possible and only the radial part remains
   !
   !          R                        r                              R
   !           /   -----  /  d (w r)    /                              /           d  (w r')  \
   !  4 Pi    |     \    |    ln       |            l+2n+2      l+2n  |             ln         |
   ! -------  | dr   )   |  ---------  | dr' f(r') r'       +  r      | dr' f(r') -----------  |
   ! 2 l + 1  |     /    |    l+2n+1   |                              |              l+2n-1    |
   !         /     -----  \  r        /                              /              r'        /
   !          0      n                 0                              r
   !
   ! The function f is a product of core and valence function
   ! (note: in the code f is defined with an additional 1/r factor)
   !
   SUBROUTINE exchange_vccvHSE(nk, bkpt, nkptd, nkpt, nkpti, ntype, neq, natd, &
                               lmax, lmaxd, nindx, maxindx, lmaxc, nindxc, &
                               maxindxc, core1, core2, lcutm, maxlmindx, &
                               bas1, bas2, jmtd, rmsh, dx, jri, jspd, jsp, &
                               maxfac, fac, sfac, nv, neigd, nbasfcn, nbands, &
                               gridf, nsymop, nsest, indx_sest, irank, &
                               a_ex, nobd, w_iks, &
                               mat_ex, te_hfex_core)

      USE m_constants
      USE m_util
      USE m_wrapper

      IMPLICIT NONE

!   -scalars -
      INTEGER, INTENT(IN)      ::  ntype, jmtd, lmaxd, jspd, jsp, neigd
      INTEGER, INTENT(IN)      ::  maxfac, nbands, maxlmindx, natd, &
                                  maxindxc, nk, maxindx, nbasfcn
      INTEGER, INTENT(IN)      ::  nkptd, nkpt, nkpti, irank
      INTEGER, INTENT(IN)      ::  nsymop
      REAL, INTENT(IN)      ::  a_ex
      REAL, INTENT(INOUT)   ::  te_hfex_core

!   - arrays -
      INTEGER, INTENT(IN)      ::  neq(ntype), lcutm(ntype), lmax(ntype), &
                                  lmaxc(ntype), jri(ntype), nv(jspd), &
                                  nindxc(0:MAXVAL(lmaxc), ntype), &
                                  nindx(0:lmaxd, ntype)
      INTEGER, INTENT(IN)      ::  nsest(nbands), indx_sest(nbands, nbands)
      INTEGER, INTENT(IN)      ::  nobd(nkpt)
      REAL, INTENT(IN)         ::  rmsh(jmtd, ntype), dx(ntype)
      REAL, INTENT(IN)         ::  bas1(jmtd, maxindx, 0:lmaxd, ntype), &
                                  bas2(jmtd, maxindx, 0:lmaxd, ntype)
      REAL, INTENT(IN)         ::  core1(jmtd, maxindxc, 0:MAXVAL(lmaxc), ntype), &
                                  core2(jmtd, maxindxc, 0:MAXVAL(lmaxc), ntype)
      REAL, INTENT(IN)         ::  fac(0:maxfac), sfac(0:maxfac)
      REAL, INTENT(IN)         ::  bkpt(3)
      REAL, INTENT(IN)         ::  gridf(jmtd, ntype)
      REAL, INTENT(IN)         ::  w_iks(neigd, nkptd, jspd)

#ifdef CPP_INVERSION
      REAL, INTENT(INOUT)  ::  mat_ex(nbasfcn*(nbasfcn + 1)/2)
#else
      COMPLEX, INTENT(INOUT)  ::  mat_ex(nbasfcn*(nbasfcn + 1)/2)
#endif
      INTEGER, PARAMETER      ::  ncut = 5                  ! cut-off value of n-summation
      INTEGER                 ::  cn                        ! counter for n-summation
      REAL                    ::  d_ln(jmtd, 0:lmaxd, 0:ncut) ! expansion coefficients of erfc(wr)/r
      ! in Legendre polynomials

!   - local scalars -
      INTEGER                 ::  iatom, ieq, itype, ic, l, l1, l2, &
                                 ll, lm, m, m1, m2, p1, p2, n, n1, n2, nn2, i, j
      INTEGER                 ::  iband1, iband2, ndb1, ndb2, ic1, ic2
      INTEGER                 ::  irecl_cmt

      REAL                    ::  rdum
      REAL                    ::  sum_offdia

      COMPLEX                 ::  cdum
!   - local arrays -
      INTEGER, ALLOCATABLE     ::  larr(:), larr2(:)
      INTEGER, ALLOCATABLE     ::  parr(:), parr2(:)

      REAL                    ::  sum_primf(jmtd), integrand(jmtd)
      REAL                    ::  primf1(jmtd), primf2(jmtd)
      REAL, ALLOCATABLE        ::  fprod(:, :), fprod2(:, :)
      REAL, ALLOCATABLE        ::  integral(:, :)

      COMPLEX                 ::  cmt(neigd, maxlmindx, natd)
      COMPLEX                 ::  exchange(nbands, nbands)
      COMPLEX, ALLOCATABLE     ::  carr(:, :), carr2(:, :), carr3(:, :)

      LOGICAL                 ::  ldum(nbands, nbands)

      ! check if a_ex is consistent
!     IF ( a_ex /= aMix_HSE ) STOP 'hsefunctional: inconsistent mixing!'

      ! read in mt wavefunction coefficients from file cmt
      irecl_cmt = neigd*maxlmindx*natd*16
      OPEN (unit=777, file='cmt', form='unformatted', access='direct', recl=irecl_cmt)
      READ (777, rec=nk) cmt(:, :, :)
      CLOSE (777)

      ALLOCATE (fprod(jmtd, 5), larr(5), parr(5))

      exchange = 0
      iatom = 0
      rdum = 0
      DO itype = 1, ntype
         ! Calculate the expansion coefficients of the potential in Legendre polynomials
         d_ln(:, :lcutm(itype), :) = calculate_coefficients(rmsh(:, itype), lcutm(itype), ncut, fac)
         DO ieq = 1, neq(itype)
            iatom = iatom + 1
            DO l1 = 0, lmaxc(itype)
               DO p1 = 1, nindxc(l1, itype)

                  DO l = 0, lcutm(itype)

                     ! Define core-valence product functions

                     n = 0
                     DO l2 = 0, lmax(itype)
                        IF (l < ABS(l1 - l2) .OR. l > l1 + l2) CYCLE

                        DO p2 = 1, nindx(l2, itype)
                           n = n + 1
                           m = SIZE(fprod, 2)
                           IF (n > m) THEN
                              ALLOCATE (fprod2(jmtd, m), larr2(m), parr2(m))
                              fprod2 = fprod; larr2 = larr; parr2 = parr
                              DEALLOCATE (fprod, larr, parr)
                              ALLOCATE (fprod(jmtd, m + 5), larr(m + 5), parr(m + 5))
                              fprod(:, :m) = fprod2
                              larr(:m) = larr2
                              parr(:m) = parr2
                              DEALLOCATE (fprod2, larr2, parr2)
                           END IF
                           fprod(:, n) = (core1(:, p1, l1, itype) &
                                          *bas1(:, p2, l2, itype) &
                                          + core2(:, p1, l1, itype) &
                                          *bas2(:, p2, l2, itype))/rmsh(:, itype)
                           larr(n) = l2
                           parr(n) = p2
                        END DO
                     END DO

                     ! Evaluate radial integrals (special part of Coulomb matrix : contribution from single MT)

                     ALLOCATE (integral(n, n), carr(n, nbands), &
                               carr2(n, nv(jsp)), carr3(n, nv(jsp)))

                     DO i = 1, n
                        ! Initialization of n Summation
                        sum_primf = 0.0
                        ! Integration over r'
                        DO cn = 0, ncut ! Summation over Legendre polynomials
                           primf1 = 0.0
                           primf2 = 0.0
                           ! Calculate integral for 0 < r' < r
                           CALL primitivef(primf1, fprod(:, i)*rmsh(:, itype)**(l + 2*cn + 1), &
                                           rmsh, dx, jri, jmtd, itype, ntype)
                           ! Calculate integral for r < r' < R
                           CALL primitivef(primf2, d_ln(:, l, cn)*fprod(:, i)/rmsh(:, itype)**(l + 2*cn), &
                                           rmsh, dx, jri, jmtd, -itype, ntype)  ! -itype is to enforce inward integration
                           ! Multiplication with appropriate prefactors
                           primf1 = primf1/rmsh(:, itype)**(l + 2*cn)*d_ln(:, l, cn)
                           primf2 = primf2*rmsh(:, itype)**(l + 2*cn + 1)
                           ! Summation over n
                           sum_primf = sum_primf + primf1 + primf2
                        END DO
                        DO j = 1, n
                           ! Integration over r
                           integrand = fprod(:, j)*sum_primf
                           integral(i, j) = fpi_const/(2*l + 1)* &
                                            intgrf(integrand, jri, jmtd, rmsh, &
                                                   dx, ntype, itype, gridf)
                        END DO

                     END DO

                     ! Add everything up

                     DO m1 = -l1, l1
                        DO m = -l, l
                           m2 = m1 + m

                           carr = 0
                           ic = 0
                           DO n1 = 1, nbands

                              DO i = 1, n
                                 ll = larr(i)
                                 IF (ABS(m2) > ll) CYCLE

                                 lm = SUM((/((2*l2 + 1)*nindx(l2, itype), l2=0, ll - 1)/)) &
                                      + (m2 + ll)*nindx(ll, itype) + parr(i)

                                 carr(i, n1) = cmt(n1, lm, iatom) &
                                               *gaunt(l1, ll, l, m1, m2, m, maxfac, fac, sfac)

                              END DO
                              DO n2 = 1, nsest(n1)!n1
                                 nn2 = indx_sest(n2, n1)
                                 exchange(nn2, n1) = exchange(nn2, n1) &
                                                     + DOT_PRODUCT(carr(:, n1), MATMUL(integral, carr(:, nn2)))
                              END DO
                           END DO
                        END DO
                     END DO

                     DEALLOCATE (integral, carr, carr2, carr3)

                  END DO
               END DO
            END DO
         END DO
      END DO

#ifdef CPP_INVERSION
      IF (ANY(ABS(aimag(exchange)) > 10.0**-10)) THEN
         IF (irank == 0) WRITE (6, '(A)') 'exchangeCore: Warning! Unusually large imaginary component.'
         WRITE (*, *) MAXVAL(ABS(aimag(exchange)))
         STOP 'exchangeCore: Unusually large imaginary component.'
      END IF
#endif

      DO n1 = 1, nobd(nk)
         te_hfex_core = te_hfex_core - a_ex*w_iks(n1, nk, jsp)*exchange(n1, n1)
      END DO

      ! add the core-valence contribution to the exchange matrix mat_ex

      ic = 0
      sum_offdia = 0
      DO n1 = 1, nbands
         DO n2 = 1, n1
            ic = ic + 1
            mat_ex(ic) = mat_ex(ic) + CONJG(exchange(n2, n1))/nsymop
         END DO
      END DO

   END SUBROUTINE exchange_vccvHSE

   ! Include the core core exchange for the HSE functional
   ! The HSE potential is the Coulomb one attenuated with the complementary
   ! error function which can be expanded in Legendre polynomials (r' > r)
   !
   !                      -----                l+2n
   ! erfc( w |r - r'| )    \                  r           /         \
   ! ------------------ =   )    d  (w r')  ---------  P ( cos(r,r') )
   !    | r - r' |         /      ln          l+2n+1    l \         /
   !                      -----              r'
   !                       l,n
   !
   ! Then an analytical integration of the angular part of the exchange
   ! potential is possible and only the radial part remains
   !
   !          R                        r                              R
   !           /   -----  /  d (w r)    /                              /           d  (w r')  \
   !  4 Pi    |     \    |    ln       |            l+2n+2      l+2n  |             ln         |
   ! -------  | dr   )   |  ---------  | dr' f(r') r'       +  r      | dr' f(r') -----------  |
   ! 2 l + 1  |     /    |    l+2n+1   |                              |              l+2n-1    |
   !         /     -----  \  r        /                              /              r'        /
   !          0      n                 0                              r
   !
   ! The function f is a product of core and valence function
   ! (note: in the code f is defined with an additional 1/r factor)
   !
   SUBROUTINE exchange_ccccHSE( &
      ! Input
      nk, nkpti, nw, nwd, ntype, neq, natd, lmaxcd, lmaxc, &
      nindxc, maxindxc, ncst, ncstd, jmtd, jri, &
      rmsh, dx, lmaxd, core1, core2, bkpt, gridf, &
      invsat, invsatnr, wtkpt, maxfac, fac, a_ex, irank, &
      ! Output
      te_hfex)

      USE m_constants
      USE m_util
      USE m_wrapper
      USE m_gaunt
      USE m_trafo

      IMPLICIT NONE

      ! - scalars -
      INTEGER, INTENT(IN)    ::  nk, nkpti, nw, nwd, ntype, natd, ncstd
      INTEGER, INTENT(IN)    ::  lmaxcd, lmaxd, irank
      INTEGER, INTENT(IN)    ::  jmtd, maxindxc, maxfac

      REAL, INTENT(IN)    ::  a_ex
      REAL, INTENT(INOUT) ::  te_hfex

      ! - arays -
      INTEGER, INTENT(IN)    ::  neq(ntype), ncst(ntype), lmaxc(ntype)
      INTEGER, INTENT(IN)    ::  nindxc(0:lmaxcd, ntype)
      INTEGER, INTENT(IN)    ::  jri(ntype)
      INTEGER, INTENT(IN)    ::  invsat(natd), invsatnr(natd)

      REAL, INTENT(IN)    ::  rmsh(jmtd, ntype), dx(ntype)
      REAL, INTENT(IN)    ::  core1(jmtd, maxindxc, 0:lmaxcd, ntype), &
                             core2(jmtd, maxindxc, 0:lmaxcd, ntype)
      REAL, INTENT(IN)    ::  bkpt(3)
      REAL, INTENT(IN)    ::  gridf(jmtd, ntype)
      REAL, INTENT(IN)    ::  wtkpt(nkpti, nwd)
      REAL, INTENT(IN)    ::  fac(0:maxfac)

      ! - local scalars -
      INTEGER               ::  itype, ieq, icst, icst1, icst2, iatom, iatom0
      INTEGER               ::  l1, l2, l, ll, llmax
      INTEGER               ::  m1, m2, m, mm
      INTEGER               ::  n1, n2, n

      REAL                  ::  rdum, rdum1
      ! - local arrays -
      INTEGER               ::  point(maxindxc, -lmaxcd:lmaxcd, 0:lmaxcd, natd)
      REAL                  ::  rprod(jmtd), primf1(jmtd), primf2(jmtd), sum_primf(jmtd), integrand(jmtd)
      COMPLEX               ::  exch(ncstd, ncstd)

      INTEGER, PARAMETER    ::  ncut = 5                     ! cut-off value of n-summation
      INTEGER               ::  cn                           ! counter for n-summation
      REAL                  ::  d_ln(jmtd, 0:2*lmaxcd, 0:ncut) ! expansion coefficients of erfc(wr)/r
      ! in Legendre polynomials
      CHARACTER*100         :: outtext

!   IF ( a_ex /= aMix_HSE ) STOP 'hsefunctional: mixing parameter inconsistent'

!   IF ( irank == 0 ) THEN
!     WRITE(outtext,'(A)') new_line('n') // new_line('n') // '### core-core-core-core exchange ###'
!     CALL writeout(outtext,irank)
!     WRITE(outtext,'(A)') new_line('n') // '        k-point       band    exchange'
!     CALL writeout(outtext,irank)
!   END IF

      ! set up point
      icst = 0
      iatom = 0
      DO itype = 1, ntype
         DO ieq = 1, neq(itype)
            iatom = iatom + 1
            DO l = 0, lmaxc(itype)
               DO m = -l, l
                  DO n = 1, nindxc(l, itype)
                     icst = icst + 1
                     point(n, m, l, iatom) = icst
                  END DO
               END DO
            END DO
         END DO
      END DO

      llmax = 2*lmaxcd
      exch = 0
      iatom0 = 0
      DO itype = 1, ntype
         ! Calculate the expansion coefficients of the potential in Legendre polynomials
         d_ln(:, 0:2*lmaxc(itype), :) = calculate_coefficients(rmsh(:, itype), 2*lmaxc(itype), ncut, fac)

         DO l1 = 0, lmaxc(itype)  ! left core state

            DO l2 = 0, lmaxc(itype)  ! right core state
               DO l = 0, lmaxc(itype)   ! occupied core state

                  DO ll = ABS(l1 - l), l1 + l
                     IF (ll < ABS(l - l2) .OR. ll > l + l2) CYCLE
                     IF (MOD(l + l1 + ll, 2) /= 0) CYCLE
                     IF (MOD(l + l2 + ll, 2) /= 0) CYCLE

                     DO m1 = -l1, l1
                        m2 = m1
                        IF (ABS(m2) > l2) CYCLE
                        DO m = -l, l
                           mm = m - m1
                           IF (ABS(mm) > ll) CYCLE
                           rdum = fpi_const/(2*ll + 1)*gaunt1(l, ll, l1, m, mm, m1, llmax)*gaunt1(l, ll, l2, m, mm, m2, llmax)

                           DO n = 1, nindxc(l, itype)
                              DO n2 = 1, nindxc(l2, itype)
                                 rprod(:) = (core1(:, n, l, itype)*core1(:, n2, l2, itype) &
                                             + core2(:, n, l, itype)*core2(:, n2, l2, itype))/rmsh(:, itype)

                                 ! Initialization of n Summation
                                 sum_primf = 0.0
                                 ! Integration over r'
                                 DO cn = 0, ncut ! Summation over Legendre polynomials
                                    primf1 = 0.0
                                    primf2 = 0.0
                                    ! Calculate integral for 0 < r' < r
                                    CALL primitivef(primf1, rprod(:)*rmsh(:, itype)**(ll + 2*cn + 1), rmsh, dx, jri, jmtd, itype, ntype)
                                    ! Calculate integral for r < r' < R
                                    !-itype is to enforce inward integration
                                    CALL primitivef(primf2, d_ln(:, ll, cn)*rprod(:)/rmsh(:, itype)**(ll + 2*cn), rmsh, dx, jri, jmtd, -itype, ntype)
                                    ! Multiplication with appropriate prefactors
                                    primf1 = primf1/rmsh(:, itype)**(ll + 2*cn)*d_ln(:, ll, cn)
                                    primf2 = primf2*rmsh(:, itype)**(ll + 2*cn + 1)
                                    ! Summation over n
                                    sum_primf = sum_primf + primf1 + primf2
                                 END DO

                                 DO n1 = 1, nindxc(l1, itype)

                                    rprod(:) = (core1(:, n, l, itype)*core1(:, n1, l1, itype) &
                                                + core2(:, n, l, itype)*core2(:, n1, l1, itype))/rmsh(:, itype)

                                    integrand = rprod*sum_primf

                                    rdum1 = rdum*intgrf(integrand, jri, jmtd, rmsh, dx, ntype, itype, gridf)

                                    iatom = iatom0
                                    DO ieq = 1, neq(itype)
                                       iatom = iatom + 1
                                       icst1 = point(n1, m1, l1, iatom)
                                       icst2 = point(n2, m2, l2, iatom)
                                       exch(icst1, icst2) = exch(icst1, icst2) + rdum1
                                    END DO
                                 END DO  !n1

                              END DO  !n2
                           END DO  !n

                        END DO  !m
                     END DO  !m1

                  END DO  !ll

               END DO  !l
            END DO  !l2
         END DO  !l1
         iatom0 = iatom0 + neq(itype)
      END DO  !itype

# ifdef CPP_INVERSION
      CALL symmetrize(exch, ncstd, ncstd, 3, .FALSE., &
                      ntype, ntype, neq, lmaxc, lmaxcd, &
                      nindxc, natd, invsat, invsatnr)
      IF (ANY(ABS(aimag(exch)) > 1E-6)) STOP 'exchange_cccc: exch possesses significant imaginary part'
# endif
!   DO icst = 1,ncstd
!     IF ( irank == 0 ) &
!       WRITE(6,'(    ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,F12.5)')bkpt,icst,REAL(exch(icst,icst))*(-27.211608)
!   END DO

      ! add core exchange contributions to the te_hfex

      DO icst1 = 1, ncstd
         te_hfex = te_hfex - a_ex*wtkpt(nk, nw)*exch(icst1, icst1)
      END DO

   END SUBROUTINE exchange_ccccHSE

   ! The expansion coefficients of the attenuated Coulomb interaction in Legendre polynomials are calculated.
   ! The attenuation is expressed in terms of a complementary error function
   !
   !                    -----              l+2n
   !  erfc(w|r-r'|)      \                r           /           \
   ! --------------- =    )    d  (wr') ---------  P (  cos(r,r')  )
   !     |r-r'|          /      ln        l+2n+1    l \           /
   !                    -----            r'
   !                     l,n
   ! where w is omega of the HSE functional, P_l are the Legendre polynomials and r' > r
   !
   ! The coefficients are given as
   !                                    l-1
   !                                2  -----   i+1  2i+1
   !                       1      -x    \     2    x
   !  d (x) = erfc(x) + -------  e       )   ------------
   !   l,0                ____          /     (2i + 1)!!
   !                     V Pi          -----
   !                                    i=0
   !
   !                                       n-1
   !               n     2                -----      i+1  l+i+1  2l+2n+2i+1
   !           (-1)    -x        2l+1      \     (-1)    2      x
   ! d (x)  = ------- e     -------------   )   ----------------------------
   !  l,n>0     ____         n (2l+2n+1)   /      i! (n-i-1)! (2l+2i+1)!!
   !           V Pi                       -----
   !                                       i=0
   !
   ! Input:  rmsh - grid of radial position
   !         lmax - maximum l value
   !         ncut - maximum value for which coefficients are calculated
   !         fac  - (optional) known factorials
   ! Return: d_ln's for all r in rmsh, all 0 <= l <= lmax and all n <= ncut
   FUNCTION calculate_coefficients(rmsh, lmax, ncut, fac) RESULT(d_ln)

      USE m_constants

      IMPLICIT NONE

      REAL, INTENT(IN)           :: rmsh(:)
      INTEGER, INTENT(IN)        :: lmax, ncut
      REAL, INTENT(IN), OPTIONAL :: fac(0:)

      INTEGER :: ci, cn, cl, fac_size
      REAL    :: r1_spi, r2l_1, ri, r2l_2i_1, rn, r2l_2n_1
      REAL    :: fn(0:MAX(0, ncut - 1))
      REAL    :: rx(SIZE(rmsh)), x2(SIZE(rmsh)), x2n(SIZE(rmsh))
      REAL    :: init(SIZE(rmsh)), addend(SIZE(rmsh))
      REAL    :: erfc_x(SIZE(rmsh)), exp_x2(SIZE(rmsh))
      REAL    :: d_ln(SIZE(rmsh), 0:lmax, 0:ncut)

      ! Initialize factorials
      IF (PRESENT(fac)) THEN
         fac_size = MIN(ncut, SIZE(fac))
         fn(0:fac_size - 1) = fac(0:fac_size - 1)
      ELSE
         fac_size = 1
         fn(0) = 1.0
      END IF
      DO cn = fac_size, ncut - 1
         fn(cn) = fn(cn - 1)*cn
      END DO

      ! Initialize 1 / sqrt(pi)
      r1_spi = 1.0/SQRT(pi_const)

      ! Calculate x and x^2 and the functions erfc(x) and exp(x^2)
      rx = omega_HSe*rmsh
      x2 = rx**2
      erfc_x = erfc(rx)
      exp_x2 = EXP(-x2)

      ! Calculate the first coefficient
      !                                    l-1
      !                                2  -----   i+1  2i+1
      !                       1      -x    \     2    x
      !  d (x) = erfc(x) + -------  e       )   ------------
      !   l,0                ____          /     (2i + 1)!!
      !                     V Pi          -----
      !                                    i=0
      ! Initialize the first coefficient
      r2l_1 = 1.0
      init = 2.0*rx*r1_spi*exp_x2
      addend = init
      IF (lmax >= 0) d_ln(:, 0, 0) = erfc_x
      IF (lmax >= 1) d_ln(:, 1, 0) = addend
      ! Iterative calculation of the other coefficients
      DO cl = 2, lmax
         r2l_1 = r2l_1 + 2.0
         addend = 2.0*addend*x2/r2l_1
         d_ln(:, cl, 0) = d_ln(:, cl - 1, 0) + addend
      END DO
      ! Adding the erfc(x) term which is present in all coefficients
      FORALL (cl=1:lmax)
      d_ln(:, cl, 0) = erfc_x + d_ln(:, cl, 0)
      END FORALL

      ! Calculate all other coefficients
      !                                       n-1
      !               n     2                -----      i+1  l+i+1  2l+2n+2i+1
      !           (-1)    -x        2l+1      \     (-1)    2      x
      ! d (x)  = ------- e     -------------   )   ----------------------------
      !  l,n>0     ____         n (2l+2n+1)   /      i! (n-i-1)! (2l+2i+1)!!
      !           V Pi                       -----
      !                                       i=0
      ! Initialize helper variables
      r2l_1 = -1.0
      DO cl = 0, lmax
         ! Calculation of the l-dependent part
         r2l_1 = r2l_1 + 2.0
         addend = init*x2               ! addend ~ 2^(l+1) x^(2l+1) / (2l-1)!!
         init = 2.0*addend/r2l_1
         ! Summation over all i
         ! i = 0 term
         DO cn = 1, ncut
            d_ln(:, cl, cn) = addend/fn(cn - 1)
         END DO
         ! higher values of i
         ri = 0.0
         r2l_2i_1 = r2l_1
         DO ci = 1, ncut - 1
            ri = ri + 1.0
            r2l_2i_1 = r2l_2i_1 + 2.0
            addend = -2.0*addend*x2/(ri*r2l_2i_1)       ! addend ~ (-1)^i+1 2^(l+i+1) x^(2l+2i+1) (2l+1) / i!(2l+2i+1)!!
            DO cn = ci + 1, ncut
               d_ln(:, cl, cn) = d_ln(:, cl, cn) + addend/fn(cn - ci - 1)
            END DO
         END DO
         r2l_2n_1 = r2l_1
         ! Divide by n and l dependent part
         DO cn = 1, ncut
            r2l_2n_1 = r2l_2n_1 + 2.0
            d_ln(:, cl, cn) = d_ln(:, cl, cn)/r2l_2n_1
         END DO
      END DO
      ! Resolve only-n dependent part
      rn = 1.0
      x2n = x2
      IF (ncut >= 1) THEN
         FORALL (cl=0:lmax)
         d_ln(:, cl, 1) = d_ln(:, cl, 1)*x2
         END FORALL
      END IF
      DO cn = 2, ncut
         rn = rn + 1.0
         x2n = -x2n*x2
         FORALL (cl=0:lmax)
         d_ln(:, cl, cn) = d_ln(:, cl, cn)*x2n/rn
         END FORALL
      END DO

   END FUNCTION calculate_coefficients

END MODULE m_hsefunctional
