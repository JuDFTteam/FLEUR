!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_tlmplm_vtau
   USE m_judft
   IMPLICIT NONE

CONTAINS
   SUBROUTINE tlmplm_vtau(n, sphhar, atoms, sym, enpara, nococonv, &
       ilSpinPr, ilSpin, iSpinV, fmpi, v, vTau, input, hub1data, td, ud)
      !! Compute V_tau contribution to the local Hamiltonian.
      !!
      !! The MetaGGA potential derivative V_tau = dE_xc/d(tau) acts on the
      !! Hamiltonian through the operator  -(1/2) nabla . (V_tau nabla).
      !! After integration by parts, the matrix elements become:
      !!
      !!   H_tau_{L'L} = (1/2) sum_lh int dV V_tau^{lh}(r) Y_{lh}
      !!                   * nabla[u_{l',order'} Y_{l'}^{m'*}]
      !!                   . nabla[u_{l,order} Y_l^m]
      !!
      !! The gradient separates into radial and angular parts:
      !!   nabla(u_l Y_{lm}) = (du_l/dr) Y_{lm} r_hat + (u_l/r) nabla_Omega Y_{lm}
      !!
      !! The dot product gives:
      !!   radial-radial: (du_{l'}/dr)(du_l/dr) Y_{l'm'}* Y_{lm}
      !!   angular-angular: (u_{l'}/r)(u_l/r) (nabla_Omega Y_{l'm'})* . (nabla_Omega Y_{lm})
      !!
      !! The angular gradient dot product satisfies:
      !!   (nabla_Omega Y_{l'm'})* . (nabla_Omega Y_{lm})
      !!     = (1/2)[l(l+1)+l'(l'+1)-lambda(lambda+1)] * G(l',lambda,l;m',mu,m) Y_{lambda,mu}
      !!
      !! Using R_l(r) = r*u_l(r) and D_l(r) = dR_l/dr - R_l/r = r*(du_l/dr):
      !!   integral = int [D_{l'}*D_l + angfac/r^2 * R_{l'}*R_l] * V_tau^{lh}(r) dr
      !!
      !! This has the same selection rules as the regular tlmplm integrals
      !! (Gaunt coefficient constraints), and the result is ADDED to td%h_loc.

      USE m_constants
      USE m_intgr, ONLY : intgr3
      USE m_genMTBasis
      USE m_gaunt, ONLY: gaunt1
      USE m_gradYlm, ONLY: Derivative
      USE m_types

      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_sphhar),   INTENT(IN)    :: sphhar
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_potden),   INTENT(IN)    :: v     ! Total potential for genMTBasis
      TYPE(t_potden),   INTENT(IN)    :: vTau  ! V_tau in lattice harmonic representation
      TYPE(t_hub1data), INTENT(IN)    :: hub1data
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      TYPE(t_usdus),    INTENT(INOUT) :: ud

      ! Indices
      INTEGER, INTENT(IN) :: n, ilSpinPr, ilSpin, iSpinV

      ! Local arrays for the 4 integral types (with gradient weights)
      REAL, ALLOCATABLE :: uvu(:,:), uvd(:,:), dvu(:,:), dvd(:,:)
      REAL, ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:), flo(:,:,:,:)
      REAL, ALLOCATABLE :: vtau_lh(:,:), x(:)

      ! Derivative arrays: D_f = df/dr - f/r, D_g = dg/dr - g/r
      REAL, ALLOCATABLE :: Df(:,:,:,:), Dg(:,:,:,:)  ! (jri, 2-components, 0:lmaxd, 2-spins)
      REAL, ALLOCATABLE :: dR_dr(:)  ! temporary for derivative computation

      COMPLEX :: cil
      REAL    :: temp, angfac, r_inv
      INTEGER :: i, l, l2, lamda, lh, lm, lmin, lmin0, lmp, lmx, lp, lh0
      INTEGER :: lp1, lpl, mem, mems, mp, mu, nh, na, m, nsym, s, lplmax, jri, comp

      jri = atoms%jri(n)
      lplmax = atoms%lmaxd*(atoms%lmaxd+3)/2

      ALLOCATE(uvu(0:lplmax,0:sphhar%nlhd)); uvu = 0.0
      ALLOCATE(uvd(0:lplmax,0:sphhar%nlhd)); uvd = 0.0
      ALLOCATE(dvu(0:lplmax,0:sphhar%nlhd)); dvu = 0.0
      ALLOCATE(dvd(0:lplmax,0:sphhar%nlhd)); dvd = 0.0

      ALLOCATE(f(atoms%jmtd,2,0:atoms%lmaxd,2), g(atoms%jmtd,2,0:atoms%lmaxd,2), x(atoms%jmtd))
      ALLOCATE(flo(atoms%jmtd,2,atoms%nlod,2))
      ALLOCATE(vtau_lh(SIZE(vTau%mt,1), 0:SIZE(vTau%mt,2)-1))

      ! Derivative arrays
      ALLOCATE(Df(atoms%jmtd, 2, 0:atoms%lmaxd, 2)); Df = 0.0
      ALLOCATE(Dg(atoms%jmtd, 2, 0:atoms%lmaxd, 2)); Dg = 0.0
      ALLOCATE(dR_dr(atoms%jmtd))

      ! Load V_tau lattice harmonic coefficients
      vtau_lh(:jri, 0:) = vTau%mt(:jri, 0:, n, iSpinV)

      ! No lh=0 offset for V_tau (it doesn't get the enpara subtraction)
      lh0 = 0

      ! Generate radial basis functions using total potential (same as in tlmplm)
      DO i = MIN(ilSpinPr,ilSpin), MAX(ilSpinPr,ilSpin)
         CALL genMTBasis(atoms, enpara, v, fmpi, n, i, ud, &
              f(:,:,:,i), g(:,:,:,i), flo(:,:,:,i), hub1data=hub1data)
      END DO

      ! Compute radial derivatives: D(r) = dR/dr - R/r
      ! for f (u functions) and g (udot functions), both components
      DO i = MIN(ilSpinPr,ilSpin), MAX(ilSpinPr,ilSpin)
         DO l = 0, atoms%lmax(n)
            DO comp = 1, 2  ! large and small components
               ! Derivative of f (u)
               CALL Derivative(f(1:jri, comp, l, i), n, atoms, dR_dr(1:jri))
               DO s = 1, jri
                  Df(s, comp, l, i) = dR_dr(s) - f(s, comp, l, i) / atoms%rmsh(s, n)
               END DO

               ! Derivative of g (udot)
               CALL Derivative(g(1:jri, comp, l, i), n, atoms, dR_dr(1:jri))
               DO s = 1, jri
                  Dg(s, comp, l, i) = dR_dr(s) - g(s, comp, l, i) / atoms%rmsh(s, n)
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE(dR_dr)

      na = atoms%firstAtom(n)
      nsym = sym%ntypsy(na)
      nh = sphhar%nlh(nsym)

      ! Generate the gradient-weighted integrals:
      ! <nabla(u_{l'}) | V_tau^{lh} | nabla(u_l)> etc.
      ! for l <= l' [lower triangle], satisfying Gaunt selection rules.

      DO lp = 0, atoms%lmax(n)
         lp1 = (lp*(lp+1))/2
         DO l = 0, lp
            lpl = lp1 + l

            DO lh = lh0, nh
               lamda = sphhar%llh(lh, nsym)
               lmin = lp - l
               lmx = lp + l
               IF ((MOD(lamda+lmx,2).EQ.1) .OR. (lamda.LT.lmin) .OR. (lamda.GT.lmx)) THEN
                  uvu(lpl,lh) = 0.0
                  uvd(lpl,lh) = 0.0
                  dvu(lpl,lh) = 0.0
                  dvd(lpl,lh) = 0.0
               ELSE
                  ! Angular gradient factor: (1/2)[l(l+1) + l'(l'+1) - lambda(lambda+1)]
                  angfac = 0.5 * REAL(l*(l+1) + lp*(lp+1) - lamda*(lamda+1))

                  ! uvu: <nabla(u_{l'}) | V_tau | nabla(u_l)>
                  DO i = 1, jri
                     r_inv = 1.0 / atoms%rmsh(i, n)
                     x(i) = ( (Df(i,1,lp,ilSpinPr)*Df(i,1,l,ilSpin) + Df(i,2,lp,ilSpinPr)*Df(i,2,l,ilSpin)) &
                            + angfac * r_inv * r_inv * &
                              (f(i,1,lp,ilSpinPr)*f(i,1,l,ilSpin) + f(i,2,lp,ilSpinPr)*f(i,2,l,ilSpin)) &
                            ) * vtau_lh(i, lh)
                  END DO
                  CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                  uvu(lpl,lh) = 0.5 * temp  ! Factor 1/2 from operator -(1/2)nabla.Vtau.nabla

                  ! uvd: <nabla(u_{l'}) | V_tau | nabla(udot_l)>
                  DO i = 1, jri
                     r_inv = 1.0 / atoms%rmsh(i, n)
                     x(i) = ( (Df(i,1,lp,ilSpinPr)*Dg(i,1,l,ilSpin) + Df(i,2,lp,ilSpinPr)*Dg(i,2,l,ilSpin)) &
                            + angfac * r_inv * r_inv * &
                              (f(i,1,lp,ilSpinPr)*g(i,1,l,ilSpin) + f(i,2,lp,ilSpinPr)*g(i,2,l,ilSpin)) &
                            ) * vtau_lh(i, lh)
                  END DO
                  CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                  uvd(lpl,lh) = 0.5 * temp

                  ! dvu: <nabla(udot_{l'}) | V_tau | nabla(u_l)>
                  DO i = 1, jri
                     r_inv = 1.0 / atoms%rmsh(i, n)
                     x(i) = ( (Dg(i,1,lp,ilSpinPr)*Df(i,1,l,ilSpin) + Dg(i,2,lp,ilSpinPr)*Df(i,2,l,ilSpin)) &
                            + angfac * r_inv * r_inv * &
                              (g(i,1,lp,ilSpinPr)*f(i,1,l,ilSpin) + g(i,2,lp,ilSpinPr)*f(i,2,l,ilSpin)) &
                            ) * vtau_lh(i, lh)
                  END DO
                  CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                  dvu(lpl,lh) = 0.5 * temp

                  ! dvd: <nabla(udot_{l'}) | V_tau | nabla(udot_l)>
                  DO i = 1, jri
                     r_inv = 1.0 / atoms%rmsh(i, n)
                     x(i) = ( (Dg(i,1,lp,ilSpinPr)*Dg(i,1,l,ilSpin) + Dg(i,2,lp,ilSpinPr)*Dg(i,2,l,ilSpin)) &
                            + angfac * r_inv * r_inv * &
                              (g(i,1,lp,ilSpinPr)*g(i,1,l,ilSpin) + g(i,2,lp,ilSpinPr)*g(i,2,l,ilSpin)) &
                            ) * vtau_lh(i, lh)
                  END DO
                  CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                  dvd(lpl,lh) = 0.5 * temp
               END IF
            END DO
         END DO
      END DO

      ! Assemble the V_tau contribution into td%h_loc using the same
      ! angular coupling (Gaunt coefficients) as in tlmplm.
      s = td%h_loc2(n) ! Offset for udot elements
      ! (lm)' loop:
      DO lp = 0, atoms%lmax(n)
         lp1 = (lp*(lp+1))/2
         DO mp = -lp, lp
            lmp = lp*(lp+1) + mp
            ! lh loop:
            DO lh = lh0, nh
               lamda = sphhar%llh(lh, nsym)
               lmin0 = ABS(lp - lamda)
               IF (lmin0.GT.lp) CYCLE
               ! Ensure l+l'+lamda even
               lmx = lp - MOD(lamda, 2)
               mems = sphhar%nmem(lh, nsym)
               DO mem = 1, mems
                  mu = sphhar%mlh(mem, lh, nsym)
                  m = mp - mu
                  lmin = MAX(lmin0, ABS(m))
                  l2 = ABS(lmx - lmin)
                  lmin = lmin + MOD(l2, 2)
                  DO l = lmin, lmx, 2
                     lm = l*(l+1) + m
                     IF (lm.GT.lmp) CYCLE
                     lpl = lp1 + l
                     cil = ImagUnit**(l-lp) * sphhar%clnu(mem, lh, nsym) &
                         * gaunt1(lp, lamda, l, mp, mu, m, atoms%lmaxd)

                     ! Add V_tau contribution to the local Hamiltonian
                     td%h_loc(lmp,lm,n,ilSpinPr,ilSpin)     = td%h_loc(lmp,lm,n,ilSpinPr,ilSpin)     + cil*uvu(lpl,lh)
                     td%h_loc(lmp,lm+s,n,ilSpinPr,ilSpin)   = td%h_loc(lmp,lm+s,n,ilSpinPr,ilSpin)   + cil*uvd(lpl,lh)
                     td%h_loc(lmp+s,lm,n,ilSpinPr,ilSpin)   = td%h_loc(lmp+s,lm,n,ilSpinPr,ilSpin)   + cil*dvu(lpl,lh)
                     td%h_loc(lmp+s,lm+s,n,ilSpinPr,ilSpin) = td%h_loc(lmp+s,lm+s,n,ilSpinPr,ilSpin) + cil*dvd(lpl,lh)
                     ! Hermitian conjugate for the upper triangle
                     IF (lm.NE.lmp) THEN
                        td%h_loc(lm,lmp,n,ilSpinPr,ilSpin)     = td%h_loc(lm,lmp,n,ilSpinPr,ilSpin)     + CONJG(cil*uvu(lpl,lh))
                        td%h_loc(lm,lmp+s,n,ilSpinPr,ilSpin)   = td%h_loc(lm,lmp+s,n,ilSpinPr,ilSpin)   + CONJG(cil*dvu(lpl,lh))
                        td%h_loc(lm+s,lmp,n,ilSpinPr,ilSpin)   = td%h_loc(lm+s,lmp,n,ilSpinPr,ilSpin)   + CONJG(cil*uvd(lpl,lh))
                        td%h_loc(lm+s,lmp+s,n,ilSpinPr,ilSpin) = td%h_loc(lm+s,lmp+s,n,ilSpinPr,ilSpin) + CONJG(cil*dvd(lpl,lh))
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE tlmplm_vtau
END MODULE m_tlmplm_vtau
