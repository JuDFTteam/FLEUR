!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_pwden_kinEnergyDen
   !! Module to compute the interstitial (plane-wave) kinetic energy density.
   !!
   !! For each k-point, the routine computes:
   !!   τ_IS(r) = (1/2) Σ_ν f_ν |∇ψ_ν(r)|²
   !!
   !! using FFT-based gradient evaluation:
   !!   ∂ψ/∂x_j(r) = FFT^{-1}[ (k+G)_j × c_G^ν ]
   !!
   !! The pattern is identical to the force-related KED in pwden.F90 (see the
   !! `input%l_f` block) but without the eigenvalue-weighted density subtraction.
   !!
   !! After accumulating over all bands, the real-space KED is FFT'd back to
   !! reciprocal space and stored as star coefficients in kinEnergyDen%pw.

CONTAINS

   SUBROUTINE pwden_kinEnergyDen(stars, kpts, input, cell, atoms, sym, &
                                  ikpt, jspin, lapw, ne, we, kinEnergyDen, zMat)
      !! Compute the interstitial kinetic energy density for one k-point.
      !!
      !! This is called once per k-point from cdnval_kinEnergyDen.
      !! The result is added to kinEnergyDen%pw(:, jspin).

      USE m_types
      USE m_types_fftGrid
      USE m_fft_interface
      USE m_juDFT

      IMPLICIT NONE

      TYPE(t_lapw), INTENT(IN)        :: lapw
      TYPE(t_input), INTENT(IN)       :: input
      TYPE(t_sym), INTENT(IN)         :: sym
      TYPE(t_stars), INTENT(IN)       :: stars
      TYPE(t_cell), INTENT(IN)        :: cell
      TYPE(t_kpts), INTENT(IN)        :: kpts
      TYPE(t_atoms), INTENT(IN)       :: atoms
      TYPE(t_mat), INTENT(IN)         :: zMat
      TYPE(t_potden), INTENT(INOUT)   :: kinEnergyDen

      REAL, INTENT(IN)       :: we(:)      !! occupation weights (ne)
      INTEGER, INTENT(IN)    :: ne         !! number of occupied states
      INTEGER, INTENT(IN)    :: ikpt       !! k-point index
      INTEGER, INTENT(IN)    :: jspin      !! spin index

      ! Local variables
      TYPE(t_fftGrid) :: stateDeriv, ekinGrid
      INTEGER :: nu, iv, ir, istr, i, j
      REAL    :: s, xk(3), stateRadius, stateFFTRadius
      REAL    :: wtf(ne)
      COMPLEX :: tempState(lapw%nv(jspin))
      COMPLEX, ALLOCATABLE :: cwk(:)
      LOGICAL :: forw

      INTEGER, ALLOCATABLE :: stateIndices(:), fieldSphereIndices(:)

      CALL timestart("pwden_kinEnergyDen")

      ! ---------------------------------------------------------------
      ! Determine FFT grid radius (same logic as pwden)
      ! ---------------------------------------------------------------
      stateRadius = MAXVAL(ABS(lapw%rk(:lapw%nv(jspin), jspin)))
      stateRadius = stateRadius + SQRT(DOT_PRODUCT(kpts%bk(:,ikpt), kpts%bk(:,ikpt)))
      stateFFTRadius = 2.0 * stateRadius

      ! ---------------------------------------------------------------
      ! Initialize FFT grids
      ! ---------------------------------------------------------------
      CALL stateDeriv%init(cell, sym, stateFFTRadius + 0.001)
      CALL ekinGrid%init(cell, sym, stateFFTRadius + 0.001)
      ekinGrid%grid(:) = CMPLX(0.0, 0.0)

      ALLOCATE(stateIndices(lapw%nv(jspin)))
      CALL stateDeriv%fillStateIndexArray(lapw, jspin, stateIndices)
      CALL ekinGrid%fillFieldSphereIndexArray(stars, stateFFTRadius + 0.0008, fieldSphereIndices)

      ! Weights: w_k * f_ν / Ω
      wtf(:ne) = we(:ne) / cell%omtil

      forw = .FALSE.  ! inverse FFT direction (reciprocal → real)

      ! ---------------------------------------------------------------
      ! Loop over occupied states
      ! ---------------------------------------------------------------
      DO nu = 1, ne
         ! Loop over 3 Cartesian directions for the gradient
         DO j = 1, 3
            ! Compute (k+G)_j × z(G, ν) for each PW component
            DO iv = 1, lapw%nv(jspin)
               xk = lapw%gvec(:, iv, jspin) + lapw%bkpt
               ! Convert to Cartesian: (k+G)_j = Σ_i xk(i) * bmat(i,j)
               s = 0.0
               DO i = 1, 3
                  s = s + xk(i) * cell%bmat(i, j)
               END DO
               IF (zMat%l_real) THEN
                  tempState(iv) = s * zMat%data_r(iv, nu)
               ELSE
                  tempState(iv) = s * zMat%data_c(iv, nu)
               END IF
            END DO

            ! Put gradient component on FFT grid and transform to real space
            CALL stateDeriv%putComplexStateOnGrid(lapw, jspin, tempState)
            CALL fft_interface(3, stateDeriv%dimensions(:), stateDeriv%grid, forw, stateIndices)

            ! Accumulate KED: τ += (1/2) w_ν/Ω |∂ψ/∂x_j(r)|²
            DO ir = 0, ekinGrid%gridLength - 1
               ekinGrid%grid(ir) = ekinGrid%grid(ir) &
                                  + wtf(nu) * 0.5 * ABS(stateDeriv%grid(ir))**2
            END DO
         END DO ! j (Cartesian direction)
      END DO ! nu (states)

      ! ---------------------------------------------------------------
      ! FFT back to reciprocal space and collect star coefficients
      ! ---------------------------------------------------------------
      forw = .TRUE.
      CALL fft_interface(3, ekinGrid%dimensions(:), ekinGrid%grid, forw, fieldSphereIndices)

      ALLOCATE(cwk(stars%ng3))
      cwk = CMPLX(0.0, 0.0)
      CALL ekinGrid%takeFieldFromGrid(stars, cwk, stateFFTRadius + 0.0005)

      ! Add to interstitial KED star coefficients
      DO istr = 1, stars%ng3_fft
         kinEnergyDen%pw(istr, jspin) = kinEnergyDen%pw(istr, jspin) + cwk(istr)
      END DO

      DEALLOCATE(cwk)

      CALL timestop("pwden_kinEnergyDen")

   END SUBROUTINE pwden_kinEnergyDen

END MODULE m_pwden_kinEnergyDen
