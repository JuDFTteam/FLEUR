MODULE m_force_a4
CONTAINS
   SUBROUTINE force_a4(atoms,sym,sphhar,input,vr,force)
      !--------------------------------------------------------------------------
      ! Core density force contribution à la Rici et al.
      ! 
      ! Equation A4, Phys. Rev. B 43, 6411
      !--------------------------------------------------------------------------
      USE m_types
      USE m_constants
      USE m_intgr, ONLY : intgr0,intgr3
      USE m_differentiate, ONLY: difcub
      USE m_cdn_io

      IMPLICIT NONE

      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_sym),    INTENT(IN) :: sym
      TYPE(t_sphhar), INTENT(IN) :: sphhar
      TYPE(t_input),  INTENT(IN) :: input

      REAL, INTENT (IN)    :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      REAL, INTENT (INOUT) :: force(3,atoms%ntype,input%jspins)

      ! Local scalars
      COMPLEX, PARAMETER :: czero=(0.0,0.0)
      REAL a4_1, a4_2, qcore, s13, s23, w, xi
      INTEGER i, ir, jsp, lh, lindex, mem, mindex, n, nd, na

      ! Local arrays
      COMPLEX forc_a4(3),gv(3),ycomp1(3,-1:1)
      REAL    rhoaux(atoms%jmtd),rhoc(atoms%jmtd,atoms%ntype,input%jspins)
      REAL    tec(atoms%ntype,input%jspins),qintc(atoms%ntype,input%jspins)

      ! Set Ylm related components.
      s13 = SQRT(1.0/3.0)
      s23 = SQRT(2.0/3.0)
      ycomp1(1,0) = czero
      ycomp1(2,0) = czero
      ycomp1(3,0) = CMPLX(2.0*s13,0.0)
      ycomp1(1,-1) = CMPLX(s23,0.0)
      ycomp1(2,-1) = CMPLX(0.0,-s23)
      ycomp1(3,-1) = czero
      ycomp1(1,1) = CMPLX(-s23,0.0)
      ycomp1(2,1) = CMPLX(0.0,-s23)
      ycomp1(3,1) = czero

      CALL timestart("force_a4")

      ! Read in core density.
      CALL readCoreDensity(input,atoms,rhoc,tec,qintc)

      DO jsp = 1, input%jspins
         DO n = 1, atoms%ntype
            na = atoms%firstAtom(n)
            IF (atoms%l_geo(n)) THEN
               nd = sym%ntypsy(na)

               DO i = 1, 3
                  forc_a4(i) = czero
               END DO

               ! TODO: There is no output for this. Do we want some?

               CALL intgr0(rhoc(:,n,jsp),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),qcore)

8000           FORMAT (' FORCE_A4: core charge=',1p,e16.8)

               DO lh = 1, sphhar%nlh(nd)
                  lindex = sphhar%llh(lh,nd)

                  IF (lindex.GT.1) EXIT

                  DO i = 1, 3
                     gv(i) = czero
                  END DO

                  ! Sum over all m for a particular lattice harmonic.
                  DO mem = 1, sphhar%nmem(lh,nd)
                     mindex = sphhar%mlh(mem,lh,nd)
                     DO i = 1, 3
                        gv(i) = gv(i) + sphhar%clnu(mem,lh,nd)*ycomp1(i,mindex)
                     END DO
                  END DO

                  ! Construct the integrand rhocore*d/dr(v)*r**2.
                  ! Note: rhoc is already multiplied by r**2 and sqrt(4*pi).
                  ! difcub performs the analytic derivative of Lagrangian of 3rd order

                  xi = atoms%rmsh(1,n)
                  rhoaux(1) = difcub(atoms%rmsh(1,n),vr(1,lh,n,jsp),xi)*rhoc(1,n,jsp)

                  DO ir = 2, atoms%jri(n) - 2
                     xi = atoms%rmsh(ir,n)
                     rhoaux(ir) = difcub(atoms%rmsh(ir-1,n),vr(ir-1,lh,n,jsp),xi) * rhoc(ir,n,jsp)
                  END DO

                  ir = atoms%jri(n) - 1
                  xi = atoms%rmsh(ir,n)
                  rhoaux(ir) = difcub(atoms%rmsh(atoms%jri(n)-3,n),vr(atoms%jri(n)-3,lh,n,jsp),xi)*rhoc(ir,n,jsp)

                  ir = atoms%jri(n)
                  xi = atoms%rmsh(ir,n)
                  rhoaux(ir) = difcub(atoms%rmsh(atoms%jri(n)-3,n),vr(atoms%jri(n)-3,lh,n,jsp),xi)*rhoc(ir,n,jsp)

                  CALL intgr3(rhoaux,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),w)

                  a4_1 = 0.5*w/sfp_const

                  ! Construct the integrand rhocore*v*r.
                  ! Note: rhocore is already multiplied by r**2 and sqrt(4*pi).
                  DO ir = 1, atoms%jri(n)
                     rhoaux(ir) = rhoc(ir,n,jsp)/atoms%rmsh(ir,n)*vr(ir,lh,n,jsp)
                  END DO

                  CALL intgr3(rhoaux,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),w)

                  a4_2 = w/sfp_const

                  ! Surface contribution from non-confined core-states
                  ! Klueppelberg Sep'12 (force level 1)
                  IF (input%ctail.AND.(input%f_level.GE.1)) THEN
                     w = rhoc(atoms%jri(n),n,jsp)*vr(atoms%jri(n),lh,n,jsp)
                     w = 0.5*w/sfp_const
                  ELSE
                     w = 0
                  END IF

                  DO i = 1, 3
                     forc_a4(i) = forc_a4(i) - (a4_1+a4_2-w)*gv(i)
                  END DO
               END DO ! lh (0:sphhar%nlh(nd))

               ! Add onto existing forces.
               DO i = 1, 3
                  force(i,n,jsp) = force(i,n,jsp) + REAL(forc_a4(i))
               END DO

               ! Write out result.
               WRITE (oUnit,FMT=8010) n
               WRITE (oUnit,FMT=8020) (forc_a4(i),i=1,3)
8010           FORMAT (' FORCES: EQUATION A4 FOR ATOM TYPE',i4)
8020           FORMAT (' FX_A4=',2f10.6,' FY_A4=',2f10.6,' FZ_A4=',2f10.6)

            END IF ! atoms%l_geo(n)
         END DO ! n (1:atoms%ntype)
      END DO ! jsp (1:input%jpins)

      CALL timestop("force_a4")

   END SUBROUTINE force_a4
END MODULE m_force_a4
