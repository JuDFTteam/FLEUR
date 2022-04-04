MODULE m_tlmplm
   USE m_judft

   IMPLICIT NONE

CONTAINS
   SUBROUTINE tlmplm(n,sphhar,atoms,sym,enpara,nococonv,&
       ilSpinPr,ilSpin,iSpinV,fmpi,v,vx,input,hub1inp,hub1data,td,ud,alpha_hybrid,lh0,one,v1)
      ! Contruct the local potential matrices
      ! t_{L'L}^{\mu} = \sum_{lh} \int dV u_{l',order'}^{\mu}(r)Y_{l'}^{m'*}(\Omega)
      !                           * V_{lh}(r)Y_{lh}(\Omega)
      !                           * u_{l,order}^{\mu}(r)Y_{l}^{m}(\Omega)
      !                           * i^{l-l'}
      ! of a real valued potential V(\bm{r}). The superindex L is defined as
      ! L := (l,m,order)
      ! with order = 0 refering to radial functions u and order = 1 denoting
      ! their energy derivatives (udot). This construction is not k-dependent
      ! and therefore executed only once each scf iteration.

      USE m_constants
      USE m_intgr, ONLY : intgr3
      USE m_genMTBasis
      USE m_tlo
      USE m_gaunt, ONLY: gaunt1
      USE m_types

      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_sphhar),   INTENT(IN)    :: sphhar
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_potden),   INTENT(IN)    :: v, vx
      TYPE(t_hub1inp),  INTENT(IN)    :: hub1inp
      TYPE(t_hub1data), INTENT(IN)    :: hub1data
      TYPE(t_tlmplm),   INTENT(INOUT) :: td
      TYPE(t_usdus),    INTENT(INOUT) :: ud

      ! Indices of atom type, local spins, accessed spin of V and where to start summing lattice harmonics.
      INTEGER, INTENT(IN) :: n, ilSpinPr, ilSpin, iSpinV, lh0
      REAL,    INTENT(IN) :: alpha_hybrid

      COMPLEX, INTENT(IN) :: one ! 1 for real part of a component, i if there is an imaginary one.

      REAL, OPTIONAL, INTENT(IN) :: v1(:,:) ! DFPT related. Get t of V1, not V.

      REAL, ALLOCATABLE :: uvu(:,:), uvd(:,:), dvu(:,:), dvd(:,:)
      REAL, ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:), flo(:,:,:,:)
      REAL, ALLOCATABLE :: vr0(:,:), x(:)

      COMPLEX :: cil
      REAL    :: temp
      INTEGER :: i,l,l2,lamda,lh,lm,lmin,lmin0,lmp,lmx,lp,info,in
      INTEGER :: lp1,lpl ,mem,mems,mp,mu,nh,na,m,nsym,s,i_u,lplmax
      LOGICAL :: l_remove

      lplmax = atoms%lmaxd*(atoms%lmaxd+3)/2

      ALLOCATE(uvu(0:lplmax,0:sphhar%nlhd)); uvu = 0.0
      ALLOCATE(uvd(0:lplmax,0:sphhar%nlhd)); uvd = 0.0
      ALLOCATE(dvu(0:lplmax,0:sphhar%nlhd)); dvu = 0.0
      ALLOCATE(dvd(0:lplmax,0:sphhar%nlhd)); dvd = 0.0

      ALLOCATE(f(atoms%jmtd,2,0:atoms%lmaxd,2),g(atoms%jmtd,2,0:atoms%lmaxd,2),x(atoms%jmtd))
      ALLOCATE(flo(atoms%jmtd,2,atoms%nlod,2))
      ALLOCATE( vr0(SIZE(v%mt,1),0:SIZE(v%mt,2)-1))

      IF (.NOT.PRESENT(v1)) THEN
         vr0 = v%mt(:,:,n,iSpinV)
         IF (iSpinV<3) THEN
            vr0(:,0)=0.0
            IF (alpha_hybrid.NE.0) vr0=vr0-alpha_hybrid*vx%mt(:,:,n,iSpinV)
         ELSE
            vr0(:,0)=vr0(:,0)-0.5*nococonv%b_con(iSpinV-2,n) !Add constraining field
         END IF
      ELSE
         vr0=v1
      END IF

      DO i = MIN(ilSpinPr,ilSpin),MAX(ilSpinPr,ilSpin)
         CALL genMTBasis(atoms,enpara,v,fmpi,n,i,ud,f(:,:,:,i),g(:,:,:,i),flo(:,:,:,i),hub1data=hub1data)
      END DO

      na = SUM(atoms%neq(:n-1)) + 1
      nsym = sym%ntypsy(na)
      nh = sphhar%nlh(nsym)

      ! Generate the irreducible integrals
      ! <u_{l',order'}^{\mu}|V_{lh}|u_{l,order}^{\mu}>
      ! for l <= l' [lower triangle], but only those that will contribute!

      DO lp = 0,atoms%lmax(n)
         lp1 = (lp*(lp+1))/2
         DO l = 0, lp
            lpl = lp1 + l
            ! Remove non-spherical components for the orbitals treated with DFT+Hubbard-1
            l_remove=.FALSE.
            IF(l.EQ.lp.AND.hub1inp%l_nonsphDC) THEN
               DO i = atoms%n_u+1, atoms%n_u+atoms%n_hia
                  IF (atoms%lda_u(i)%atomType.EQ.n.AND.atoms%lda_u(i)%l.EQ.l) l_remove=.TRUE.
               END DO
            END IF

            ! Loop over the required components of the potential: must satisfy
            ! the triangular conditions and that l'+l+l_{lh} even (Gaunt
            ! coefficient selection rules)
            DO lh = lh0, nh
               lamda = sphhar%llh(lh,nsym)
               lmin = lp - l
               lmx = lp + l
               IF ((MOD(lamda+lmx,2).EQ.1) .OR. (lamda.LT.lmin) .OR. (lamda.GT.lmx) .OR. l_remove) THEN
                  uvu(lpl,lh) = 0.0
                  uvd(lpl,lh) = 0.0
                  dvu(lpl,lh) = 0.0
                  dvd(lpl,lh) = 0.0
               ELSE
                  DO i = 1,atoms%jri(n)
                     x(i) = (f(i,1,lp,ilSpinPr)*f(i,1,l,ilSpin)+f(i,2,lp,ilSpinPr)*f(i,2,l,ilSpin))* vr0(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                  uvu(lpl,lh) = temp

                  DO i = 1,atoms%jri(n)
                     x(i) = (f(i,1,lp,ilSpinPr)*g(i,1,l,ilSpin)+f(i,2,lp,ilSpinPr)*g(i,2,l,ilSpin))* vr0(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                  uvd(lpl,lh) = temp

                  DO i = 1,atoms%jri(n)
                     x(i) = (g(i,1,lp,ilSpinPr)*f(i,1,l,ilSpin)+g(i,2,lp,ilSpinPr)*f(i,2,l,ilSpin))* vr0(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                  dvu(lpl,lh) = temp

                  DO i = 1,atoms%jri(n)
                     x(i) = (g(i,1,lp,ilSpinPr)*g(i,1,l,ilSpin)+g(i,2,lp,ilSpinPr)*g(i,2,l,ilSpin))* vr0(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),temp)
                  dvd(lpl,lh) = temp
               END IF
            END DO
         END DO
      END DO

      ! Generate the various t matrices for lm <= (lm)'
      s = td%h_loc2(n) ! Offset for writing udot elements behind u
      ! (lm)' loop:
      DO lp = 0,atoms%lmax(n)
         lp1 = (lp*(lp+1))/2
         DO mp = -lp,lp
            lmp = lp* (lp+1) + mp
            ! lh loop:
            DO lh = lh0, nh
               lamda = sphhar%llh(lh,nsym)
               lmin0 = ABS(lp-lamda)
               IF (lmin0.GT.lp) CYCLE
               ! Ensure l+l'+lamda even
               lmx = lp - MOD(lamda,2)
               mems = sphhar%nmem(lh,nsym)
               DO mem = 1,mems
                  mu = sphhar%mlh(mem,lh,nsym)
                  m = mp - mu
                  lmin = MAX(lmin0,ABS(m))
                  l2 = ABS(lmx-lmin)
                  lmin = lmin + MOD(l2,2)
                  DO l = lmin,lmx,2
                     lm = l* (l+1) + m
                     IF (lm.GT.lmp) CYCLE
                     lpl = lp1 + l
                     cil = ImagUnit**(l-lp) * sphhar%clnu(mem,lh,nsym) &
                         * gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)

                     td%h_loc(lmp,lm,n,ilSpinPr,ilSpin)     =  td%h_loc(lmp,lm,n,ilSpinPr,ilSpin)     + one*cil*uvu(lpl,lh)
                     td%h_loc(lmp,lm+s,n,ilSpinPr,ilSpin)   =  td%h_loc(lmp,lm+s,n,ilSpinPr,ilSpin)   + one*cil*uvd(lpl,lh)
                     td%h_loc(lmp+s,lm,n,ilSpinPr,ilSpin)   =  td%h_loc(lmp+s,lm,n,ilSpinPr,ilSpin)   + one*cil*dvu(lpl,lh)
                     td%h_loc(lmp+s,lm+s,n,ilSpinPr,ilSpin) =  td%h_loc(lmp+s,lm+s,n,ilSpinPr,ilSpin) + one*cil*dvd(lpl,lh)
                     ! Use the fact that the t matrices are Hermitian by definition to contruct the upper triangle as well.
                     IF (lm.NE.lmp) THEN
                        td%h_loc(lm,lmp,n,ilSpinPr,ilSpin)     =  td%h_loc(lm,lmp,n,ilSpinPr,ilSpin)     + one*CONJG(cil*uvu(lpl,lh))
                        td%h_loc(lm,lmp+s,n,ilSpinPr,ilSpin)   =  td%h_loc(lm,lmp+s,n,ilSpinPr,ilSpin)   + one*CONJG(cil*dvu(lpl,lh))
                        td%h_loc(lm+s,lmp,n,ilSpinPr,ilSpin)   =  td%h_loc(lm+s,lmp,n,ilSpinPr,ilSpin)   + one*CONJG(cil*uvd(lpl,lh))
                        td%h_loc(lm+s,lmp+s,n,ilSpinPr,ilSpin) =  td%h_loc(lm+s,lmp+s,n,ilSpinPr,ilSpin) + one*CONJG(cil*dvd(lpl,lh))
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! If necessary, set up the t-matrices for the local orbitals as well.
      IF (atoms%nlo(n).GE.1) THEN
         CALL tlo(atoms,sym,sphhar,ilSpinPr,ilSpin,iSpinV,n,enpara,lh0,input,vr0,&
            na,flo,f,g,ud, td, one)
      END IF
   END SUBROUTINE tlmplm
END MODULE m_tlmplm
