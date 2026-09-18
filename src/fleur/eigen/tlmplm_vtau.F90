!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
      !!
      !! Local orbitals are covered in two parts:
      !!  - The a/b parts of the LO expansion come for free: local_hamiltonian derives
      !!    td%h_loc_LO from td%h_loc *after* this routine runs.
      !!  - The parts multiplying the LO radial function itself are built here, in the
      !!    second half of the routine, using tlo.f90 as the template: the same loops and
      !!    Gaunt assembly, with the integrand replaced by
      !!      (D_a*D_b + angfac/r^2 * R_a*R_b) * vtau(r)
      !!    and the result added to td%h_LO, td%h_LO2 (APW-LO) and td%tuloulo_newer (LO-LO).
      !!    tlo's spherical-Hamiltonian block (its `l_V1` branch) has no V_tau analogue and
      !!    is deliberately absent: V_tau is never absorbed into the radial functions, so
      !!    its lh=0 component is picked up by the lattice-harmonic loop like any other.

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

      ! Local-orbital counterparts of uvu/uvd (APW-LO) and of the LO-LO block
      REAL, ALLOCATABLE :: uvulo(:,:,:), dvulo(:,:,:), ulovulo(:,:)

      ! Derivative arrays: D_f = df/dr - f/r, D_g = dg/dr - g/r
      REAL, ALLOCATABLE :: Df(:,:,:,:), Dg(:,:,:,:)  ! (jri, 2-components, 0:lmaxd, 2-spins)
      REAL, ALLOCATABLE :: Dflo(:,:,:,:)             ! (jri, 2-components, nlod, 2-spins)
      REAL, ALLOCATABLE :: dR_dr(:)  ! temporary for derivative computation

      COMPLEX :: cil
      REAL    :: temp, angfac, r_inv
      INTEGER :: i, l, l2, lamda, lh, lm, lmin, lmin0, lmp, lmx, lp, lh0
      INTEGER :: lp1, lpl, mem, mems, mp, mu, nh, na, m, nsym, s, lplmax, jri, comp
      INTEGER :: lo, lop, loplo, mlo, s_lo, lpmin, lpmax, lpmin0, lpmax0

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
      ALLOCATE(Dflo(atoms%jmtd, 2, MAX(atoms%nlod,1), 2)); Dflo = 0.0
      ALLOCATE(dR_dr(atoms%jmtd))

      ! Local-orbital integrals (same layout as in tlo, but with lh0 = 0)
      ALLOCATE(uvulo(MAX(atoms%nlod,1), 0:atoms%lmaxd, 0:sphhar%nlhd)); uvulo = 0.0
      ALLOCATE(dvulo(MAX(atoms%nlod,1), 0:atoms%lmaxd, 0:sphhar%nlhd)); dvulo = 0.0
      ALLOCATE(ulovulo(MAX(atoms%nlod*(atoms%nlod+1)/2,1), 0:sphhar%nlhd)); ulovulo = 0.0

      CALL timestart("tlmplm_vtau")

      ! Load V_tau lattice harmonic coefficients
      vtau_lh(:jri, 0:) = vTau%mt(:jri, 0:, n, iSpinV)

      ! No lh=0 offset for V_tau (it doesn't get the enpara subtraction)
      lh0 = 0

      CALL timestart("tlmplm_vtau: basis+derivatives")
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

         ! Same for the local-orbital radial functions
         DO lo = 1, atoms%nlo(n)
            DO comp = 1, 2
               CALL Derivative(flo(1:jri, comp, lo, i), n, atoms, dR_dr(1:jri))
               DO s = 1, jri
                  Dflo(s, comp, lo, i) = dR_dr(s) - flo(s, comp, lo, i) / atoms%rmsh(s, n)
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE(dR_dr)
      CALL timestop("tlmplm_vtau: basis+derivatives")

      na = atoms%firstAtom(n)
      nsym = sym%ntypsy(na)
      nh = sphhar%nlh(nsym)

      CALL timestart("tlmplm_vtau: radial integrals")
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
               IF ((MOD(lamda+lmx,2)==1) .OR. (lamda<lmin) .OR. (lamda>lmx)) THEN
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

      CALL timestop("tlmplm_vtau: radial integrals")

      ! Assemble the V_tau contribution into td%h_loc using the same
      ! angular coupling (Gaunt coefficients) as in tlmplm.
      CALL timestart("tlmplm_vtau: assemble h_loc")
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
               IF (lmin0>lp) CYCLE
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
                     IF (lm>lmp) CYCLE
                     lpl = lp1 + l
                     cil = ImagUnit**(l-lp) * sphhar%clnu(mem, lh, nsym) &
                         * gaunt1(lp, lamda, l, mp, mu, m, atoms%lmaxd)

                     ! Add V_tau contribution to the local Hamiltonian
                     td%h_loc(lmp,lm,n,ilSpinPr,ilSpin)     = td%h_loc(lmp,lm,n,ilSpinPr,ilSpin)     + cil*uvu(lpl,lh)
                     td%h_loc(lmp,lm+s,n,ilSpinPr,ilSpin)   = td%h_loc(lmp,lm+s,n,ilSpinPr,ilSpin)   + cil*uvd(lpl,lh)
                     td%h_loc(lmp+s,lm,n,ilSpinPr,ilSpin)   = td%h_loc(lmp+s,lm,n,ilSpinPr,ilSpin)   + cil*dvu(lpl,lh)
                     td%h_loc(lmp+s,lm+s,n,ilSpinPr,ilSpin) = td%h_loc(lmp+s,lm+s,n,ilSpinPr,ilSpin) + cil*dvd(lpl,lh)
                     ! Hermitian conjugate for the upper triangle
                     IF (lm/=lmp) THEN
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

      CALL timestop("tlmplm_vtau: assemble h_loc")

      ! ================================================================
      ! Local orbitals. Structure follows tlo.f90 exactly; only the
      ! integrand differs (gradient-weighted, against V_tau instead of V).
      ! ================================================================
      IF (atoms%nlo(n) > 0) THEN
         CALL timestart("tlmplm_vtau: LO radial integrals")

         ! <nabla(u_{l'}) | V_tau | nabla(u_lo)>  and the same with udot
         DO lo = 1, atoms%nlo(n)
            l = atoms%llo(lo, n)
            DO lp = 0, atoms%lmax(n)
               lmin = ABS(lp - l)
               lmx = lp + l
               DO lh = lh0, nh
                  lamda = sphhar%llh(lh, nsym)
                  IF ((MOD(l+lp+lamda,2)==1) .OR. (lamda<lmin) .OR. (lamda>lmx)) THEN
                     uvulo(lo,lp,lh) = 0.0
                     dvulo(lo,lp,lh) = 0.0
                  ELSE
                     angfac = 0.5 * REAL(l*(l+1) + lp*(lp+1) - lamda*(lamda+1))

                     DO i = 1, jri
                        r_inv = 1.0 / atoms%rmsh(i, n)
                        x(i) = ( (Df(i,1,lp,ilSpinPr)*Dflo(i,1,lo,ilSpin) + Df(i,2,lp,ilSpinPr)*Dflo(i,2,lo,ilSpin)) &
                               + angfac * r_inv * r_inv * &
                                 (f(i,1,lp,ilSpinPr)*flo(i,1,lo,ilSpin) + f(i,2,lp,ilSpinPr)*flo(i,2,lo,ilSpin)) &
                               ) * vtau_lh(i, lh)
                     END DO
                     CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                     uvulo(lo,lp,lh) = 0.5 * temp

                     DO i = 1, jri
                        r_inv = 1.0 / atoms%rmsh(i, n)
                        x(i) = ( (Dg(i,1,lp,ilSpinPr)*Dflo(i,1,lo,ilSpin) + Dg(i,2,lp,ilSpinPr)*Dflo(i,2,lo,ilSpin)) &
                               + angfac * r_inv * r_inv * &
                                 (g(i,1,lp,ilSpinPr)*flo(i,1,lo,ilSpin) + g(i,2,lp,ilSpinPr)*flo(i,2,lo,ilSpin)) &
                               ) * vtau_lh(i, lh)
                     END DO
                     CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                     dvulo(lo,lp,lh) = 0.5 * temp
                  END IF
               END DO
            END DO
         END DO

         ! <nabla(u_lo') | V_tau | nabla(u_lo)> for lo <= lo'
         loplo = 0
         DO lop = 1, atoms%nlo(n)
            lp = atoms%llo(lop, n)
            DO lo = 1, lop
               l = atoms%llo(lo, n)
               loplo = loplo + 1
               IF (loplo > SIZE(ulovulo,1)) CALL juDFT_error("loplo too large!!!", calledby="tlmplm_vtau")
               DO lh = lh0, nh
                  lamda = sphhar%llh(lh, nsym)
                  lmin = ABS(lp - l)
                  lmx = lp + l
                  IF ((MOD(l+lp+lamda,2)==1) .OR. (lamda<lmin) .OR. (lamda>lmx)) THEN
                     ulovulo(loplo,lh) = 0.0
                  ELSE
                     angfac = 0.5 * REAL(l*(l+1) + lp*(lp+1) - lamda*(lamda+1))
                     DO i = 1, jri
                        r_inv = 1.0 / atoms%rmsh(i, n)
                        x(i) = ( (Dflo(i,1,lop,ilSpinPr)*Dflo(i,1,lo,ilSpin) &
                                + Dflo(i,2,lop,ilSpinPr)*Dflo(i,2,lo,ilSpin)) &
                               + angfac * r_inv * r_inv * &
                                 (flo(i,1,lop,ilSpinPr)*flo(i,1,lo,ilSpin) &
                                + flo(i,2,lop,ilSpinPr)*flo(i,2,lo,ilSpin)) &
                               ) * vtau_lh(i, lh)
                     END DO
                     CALL intgr3(x, atoms%rmsh(1,n), atoms%dx(n), jri, temp)
                     ulovulo(loplo,lh) = 0.5 * temp
                  END IF
               END DO
            END DO
         END DO
         CALL timestop("tlmplm_vtau: LO radial integrals")

         CALL timestart("tlmplm_vtau: assemble LO")
         mlo  = SUM(atoms%nlo(:n-1))
         s_lo = td%h_loc2_nonsph(n) ! Offset for udot elements in the LO matrices

         ! APW-LO blocks -> h_LO / h_LO2
         DO lo = 1, atoms%nlo(n)
            l = atoms%llo(lo, n)
            DO m = -l, l
               DO lh = lh0, nh
                  lamda = sphhar%llh(lh, nsym)
                  lpmin0 = ABS(l - lamda)
                  lpmax0 = l + lamda
                  ! Cap l' at the non-spherical expansion of this atom
                  lpmax = MIN(lpmax0, atoms%lnonsph(n))
                  ! Ensure l + lamda + l' even
                  lpmax = lpmax - MOD(l+lamda+lpmax, 2)
                  DO mem = 1, sphhar%nmem(lh, nsym)
                     mu = sphhar%mlh(mem, lh, nsym)
                     mp = m + mu
                     lpmin = MAX(lpmin0, ABS(mp))
                     lpmin = lpmin + MOD(ABS(lpmax-lpmin), 2)
                     DO lp = lpmin, lpmax, 2
                        lmp = lp*(lp+1) + mp
                        cil = ImagUnit**(l-lp) * sphhar%clnu(mem, lh, nsym) &
                            * gaunt1(lp, lamda, l, mp, mu, m, atoms%lmaxd)

                        td%h_LO(lmp,m,lo+mlo,ilSpinPr,ilSpin) = &
                           td%h_LO(lmp,m,lo+mlo,ilSpinPr,ilSpin) + cil*uvulo(lo,lp,lh)
                        td%h_LO(lmp+s_lo,m,lo+mlo,ilSpinPr,ilSpin) = &
                           td%h_LO(lmp+s_lo,m,lo+mlo,ilSpinPr,ilSpin) + cil*dvulo(lo,lp,lh)
                        td%h_LO2(lmp,m,lo+mlo,ilSpinPr,ilSpin) = &
                           td%h_LO2(lmp,m,lo+mlo,ilSpinPr,ilSpin) + CONJG(cil*uvulo(lo,lp,lh))
                        td%h_LO2(lmp+s_lo,m,lo+mlo,ilSpinPr,ilSpin) = &
                           td%h_LO2(lmp+s_lo,m,lo+mlo,ilSpinPr,ilSpin) + CONJG(cil*dvulo(lo,lp,lh))
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! LO-LO block -> tuloulo_newer
         DO lop = 1, atoms%nlo(n)
            lp = atoms%llo(lop, n)
            DO mp = -lp, lp
               DO lh = lh0, nh
                  lamda = sphhar%llh(lh, nsym)
                  DO mem = 1, sphhar%nmem(lh, nsym)
                     mu = sphhar%mlh(mem, lh, nsym)
                     m = mp - mu
                     DO lo = 1, lop
                        l = atoms%llo(lo, n)
                        loplo = ((lop-1)*lop)/2 + lo
                        IF ((ABS(l-lamda)<=lp) .AND. (lp<=(l+lamda)) .AND. &
                            (MOD(l+lp+lamda,2)==0) .AND. (ABS(m)<=l)) THEN
                           cil = ImagUnit**(l-lp) * sphhar%clnu(mem, lh, nsym) &
                               * gaunt1(lp, lamda, l, mp, mu, m, atoms%lmaxd)
                           td%tuloulo_newer(mp,m,lop,lo,n,ilSpinPr,ilSpin) = &
                              td%tuloulo_newer(mp,m,lop,lo,n,ilSpinPr,ilSpin) + cil*ulovulo(loplo,lh)
                           IF (lop/=lo) THEN
                              td%tuloulo_newer(m,mp,lo,lop,n,ilSpinPr,ilSpin) = &
                                 td%tuloulo_newer(m,mp,lo,lop,n,ilSpinPr,ilSpin) + CONJG(cil*ulovulo(loplo,lh))
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CALL timestop("tlmplm_vtau: assemble LO")
      END IF

      CALL timestop("tlmplm_vtau")

   END SUBROUTINE tlmplm_vtau
END MODULE m_tlmplm_vtau
