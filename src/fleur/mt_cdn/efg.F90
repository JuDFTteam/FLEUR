!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_efg

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: calc_efg

   CONTAINS

   SUBROUTINE calc_efg(atoms, sym, sphhar, fmpi, vCoul, efgTensor, l_efgValid)

      !-------------------------------------------------------------------------
      ! Computes the electric field gradient (EFG) tensor at each atomic
      ! nucleus from the l=2 channel of the muffin-tin Coulomb potential.
      !
      ! Near the nucleus the regular l=2 solution of the potential behaves as
      ! V_2m(r) -> v_2m * r^2; v_2m maps directly onto the Cartesian traceless
      ! symmetric EFG tensor V_ij.
      !-------------------------------------------------------------------------

      USE m_types
      USE m_constants
      USE m_lattHarmsSphHarmsConv

      TYPE(t_atoms),  INTENT(IN)  :: atoms
      TYPE(t_sym),    INTENT(IN)  :: sym
      TYPE(t_sphhar), INTENT(IN)  :: sphhar
      TYPE(t_mpi),    INTENT(IN)  :: fmpi
      TYPE(t_potden), INTENT(IN)  :: vCoul
      REAL,           INTENT(OUT) :: efgTensor(3,3,atoms%ntype)
      LOGICAL,        INTENT(OUT) :: l_efgValid(atoms%ntype)

      INTEGER :: iType, nd, m, lm

      COMPLEX, ALLOCATABLE :: vlm(:,:)
      COMPLEX :: v2m(-2:2)

      REAL    :: v20, p, q, u, v
      REAL    :: s5pi, s152pi

      efgTensor = 0.0
      l_efgValid = .FALSE.

      IF (fmpi%irank /= 0) RETURN

      s5pi   = SQRT(5.0/pi_const)
      s152pi = SQRT(15.0/(2.0*pi_const))

      ALLOCATE(vlm(atoms%jmtd,(atoms%lmaxd+1)**2))

      DO iType = 1, atoms%ntype

         IF (atoms%lmax(iType) < 2) CYCLE

         nd = sym%ntypsy(atoms%firstAtom(iType))

         CALL lattHarmsRepToSphHarms(sym, atoms, sphhar, iType, vCoul%mt(:,0:,iType,1), vlm)

         ! Extract the r->0 coefficient v_2m of each l=2 channel via a linear
         ! least-squares fit of vlm(r)/r^2 against r^2, using the first 4
         ! (innermost) points of the logarithmic radial mesh.
         DO m = -2, 2
            lm = 2*3 + m + 1
            v2m(m) = extrapolate_r0(vlm(:,lm), atoms%rmsh(:,iType))
         END DO

         ! vCoul is FLEUR's internal Coulomb potential, in the electron
         ! potential-energy convention (added directly, with a plain "+"
         ! sign, into the KS Hamiltonian; attractive/negative near a
         ! nucleus -- the same convention the core-state solvers in
         ! cored.F90/differ.f rely on to find bound, negative-eigenvalue
         ! states). The classic EFG formula below is written in terms of
         ! the standard electrostatic potential (positive-test-charge
         ! convention, i.e. the negative of vCoul), so the l=2 expansion
         ! coefficients are negated here before use.
         v20 = -REAL(v2m(0))
         p   = -REAL(v2m(1))
         q   = -AIMAG(v2m(1))
         u   = -REAL(v2m(2))
         v   = -AIMAG(v2m(2))

         efgTensor(3,3,iType) = s5pi*v20
         efgTensor(1,1,iType) = -0.5*s5pi*v20 + s152pi*u
         efgTensor(2,2,iType) = -0.5*s5pi*v20 - s152pi*u
         efgTensor(1,2,iType) = -s152pi*v
         efgTensor(2,1,iType) = efgTensor(1,2,iType)
         efgTensor(1,3,iType) = -s152pi*p
         efgTensor(3,1,iType) = efgTensor(1,3,iType)
         efgTensor(2,3,iType) =  s152pi*q
         efgTensor(3,2,iType) = efgTensor(2,3,iType)

         l_efgValid(iType) = .TRUE.

      END DO

      DEALLOCATE(vlm)

   END SUBROUTINE calc_efg

   COMPLEX FUNCTION extrapolate_r0(vals, rmsh)
      ! Linear least-squares fit of vals(i)/rmsh(i)**2 against rmsh(i)**2,
      ! using the first 4 (innermost) radial mesh points, extrapolated to r=0.
      COMPLEX, INTENT(IN) :: vals(:)
      REAL,    INTENT(IN) :: rmsh(:)

      INTEGER, PARAMETER :: nFit = 4
      REAL :: x(nFit), yr(nFit), yi(nFit)
      REAL :: sx, sxx, sx2
      REAL :: ar, br, ai, bi
      INTEGER :: i

      DO i = 1, nFit
         x(i)  = rmsh(i)**2
         yr(i) = REAL(vals(i)) / rmsh(i)**2
         yi(i) = AIMAG(vals(i)) / rmsh(i)**2
      END DO

      sx  = SUM(x)
      sxx = SUM(x*x)
      sx2 = REAL(nFit)*sxx - sx*sx

      br = (REAL(nFit)*SUM(x*yr) - sx*SUM(yr)) / sx2
      ar = (SUM(yr) - br*sx) / REAL(nFit)

      bi = (REAL(nFit)*SUM(x*yi) - sx*SUM(yi)) / sx2
      ai = (SUM(yi) - bi*sx) / REAL(nFit)

      extrapolate_r0 = CMPLX(ar,ai)

   END FUNCTION extrapolate_r0

END MODULE m_efg
