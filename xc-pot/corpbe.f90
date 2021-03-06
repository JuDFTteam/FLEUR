MODULE m_corpbe

!----------------------------------------------------------------------
!  official pbe correlation code. k. burke, may 14, 1996.
!----------------------------------------------------------------------
! references:
! [a] j.p.~perdew, k.~burke, and m.~ernzerhof,
!     {\sl generalized gradient approximation made simple}, sub.
!     to phys. rev.lett. may 1996.
! [b] j. p. perdew, k. burke, and y. wang, {\sl real-space cutoff
!     construction of a generalized gradient approximation:  the pw91
!     density functional}, submitted to phys. rev. b, feb. 1996.
! [c] j. p. perdew and y. wang, phys. rev. b {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE corpbe( &
      l_pbes,rs,zet,t,uu,vv,ww,lgga,lpot, &
      ec,vcup,vcdn,h,dvcup,dvcdn)
!----------------------------------------------------------------------
!  input: rs=seitz radius=(3/4pi rho)^(1/3)
!       : zet=relative spin polarization = (rhoup-rhodn)/rho
!       : t=abs(grad rho)/(rho*2.*ks*g)  -- only needed for pbe
!       : uu=(grad rho)*grad(abs(grad rho))/(rho**2 * (2*ks*g)**3)
!       : vv=(laplacian rho)/(rho * (2*ks*g)**2)
!       : ww=(grad rho)*(grad zet)/(rho * (2*ks*g)**2
!       :  uu,vv,ww, only needed for pbe potential
!       : lgga=flag to do gga (0=>lsd only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=nonlocal part of correlation energy per electron
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------

      USE m_pbecor2
      IMPLICIT NONE
      LOGICAL,INTENT(IN)    :: l_pbes
      INTEGER, INTENT (IN)  :: lgga,lpot
      REAL,    INTENT (IN)  :: rs,zet,t,uu,vv,ww
      REAL,    INTENT (OUT) :: dvcdn,dvcup,ec,h,vcdn,vcup

! thrd*=various multiples of 1/3
! numbers for use in lsd energy spin-interpolation formula, [c](9).
!      gam= 2^(4/3)-2
!      fzz=f''(0)= 8/(9*gam)
! numbers for construction of pbe
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.

      REAL, PARAMETER :: thrd=1.e0/3.e0
      REAL, PARAMETER :: thrdm=-thrd
      REAL, PARAMETER :: thrd2=2.e0*thrd
      REAL, PARAMETER :: sixthm=thrdm/2.e0
      REAL, PARAMETER :: thrd4=4.e0*thrd
      REAL, PARAMETER :: gam=0.5198420997897463295344212145565e0
      REAL, PARAMETER :: fzz=8.e0/ (9.e0*gam)
      REAL, PARAMETER :: gamma=0.03109069086965489503494086371273e0
      REAL, PARAMETER :: eta=1.e-12
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find lsd energy contributions, using [c](10) and table i[c].
! eu=unpolarized lsd correlation energy
! eurs=deu/drs
! ep=fully polarized lsd correlation energy
! eprs=dep/drs
! alfm=-spin stiffness, [c](3).
! alfrsm=-dalpha/drs
! f=spin-scaling factor from [c](9).
! construct ec, using [c](8)

      REAL :: alfm,alfrsm,b,b2,bec,bg,comm,ecrs,eczet,ep,eprs,eu,eurs, &
              f,fac,fact0,fact1,fact2,fact3,fact5,fz,g,g3,g4,gz,hb,hbt,hrs, &
              hrst,ht,htt,hz,hzt,pon,pref,q4,q5,q8,q9,rsthrd,rtrs, &
              t2,t4,t6,z4,delt,bet
!     ..

      IF (l_pbes) THEN ! PBE_sol
         bet=0.046e0
      ELSE
         bet=0.06672455060314922e0
      ENDIF
      delt=bet/gamma

      rtrs = sqrt(rs)
      CALL pbecor2(0.0310907,0.21370,7.5957,3.5876,1.6382, &
      &            0.49294,rtrs,eu,eurs)
      CALL pbecor2(0.01554535,0.20548,14.1189,6.1977,3.3662, &
      &            0.62517,rtrs,ep,eprs)
      CALL pbecor2(0.0168869,0.11125,10.357,3.6231,0.88026, &
      &            0.49671,rtrs,alfm,alfrsm)
      z4 = zet**4
      f = ((1.e0+zet)**thrd4+ (1.e0-zet)**thrd4-2.e0)/gam
      ec = eu* (1.e0-f*z4) + ep*f*z4 - alfm*f* (1.e0-z4)/fzz
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! lsd potential from [c](a1)
! ecrs = dec/drs [c](a2)
! eczet=dec/dzeta [c](a3)
! fz = df/dzeta [c](a4)
      ecrs = eurs* (1.e0-f*z4) + eprs*f*z4 - alfrsm*f* (1.e0-z4)/fzz
      fz = thrd4* ((1.e0+zet)**thrd- (1.e0-zet)**thrd)/gam
      eczet = 4.e0* (zet**3)*f* (ep-eu+alfm/fzz) + &
              fz* (z4*ep-z4*eu- (1.e0-z4)*alfm/fzz)
      comm = ec - rs*ecrs/3.e0 - zet*eczet
      vcup = comm + eczet
      vcdn = comm - eczet
      IF (lgga == 0) RETURN
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! pbe correlation energy
! g=phi(zeta), given after [a](3)
! delt=bet/gamma
! b=a of [a](8)
      g = ((1.e0+zet)**thrd2+ (1.e0-zet)**thrd2)/2.e0
      g3 = g**3
      pon = -ec/ (g3*gamma)
      b = delt/ (exp(pon)-1.e0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      q4 = 1.e0 + b*t2
      q5 = 1.e0 + b*t2 + b2*t4
      h = g3* (bet/delt)*log(1.e0+delt*q4*t2/q5)
      IF (lpot == 0) RETURN
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! energy done. now the potential, using appendix e of [b].
      g4 = g3*g
      t6 = t4*t2
      rsthrd = rs/3.e0
      gz = (((1.e0+zet)**2+eta)**sixthm- ((1.e0-zet)**2+eta)**sixthm)/ &
      &      3.e0
      fac = delt/b + 1.e0
      bg = -3.e0*b2*ec*fac/ (bet*g4)
      bec = b2*fac/ (bet*g3)
      q8 = q5*q5 + delt*q4*q5*t2
      q9 = 1.e0 + 2.e0*b*t2
      hb = -bet*g3*b*t6* (2.e0+b*t2)/q8
      hrs = -rsthrd*hb*bec*ecrs
      fact0 = 2.e0*delt - 6.e0*b
      fact1 = q5*q9 + q4*q9*q9
      hbt = 2.e0*bet*g3*t4* ((q4*q5*fact0-delt*fact1)/q8)/q8
      hrst = rsthrd*t2*hbt*bec*ecrs
      hz = 3.e0*gz*h/g + hb* (bg*gz+bec*eczet)
      ht = 2.e0*bet*g3*q9/q8
      hzt = 3.e0*gz*ht/g + hbt* (bg*gz+bec*eczet)
      fact2 = q4*q5 + b*t2* (q4*q9+q5)
      fact3 = 2.e0*b*q5*q9 + delt*fact2
      htt = 4.e0*bet*g3*t* (2.e0*b/q8- (q9*fact3/q8)/q8)
      comm = h + hrs + hrst + t2*ht/6.e0 + 7.e0*t2*t*htt/6.e0
      pref = hz - gz*t2*ht/g
      fact5 = gz* (2.e0*ht+t*htt)/g
      comm = comm - pref*zet - uu*htt - vv*ht - ww* (hzt-fact5)
      dvcup = comm + pref
      dvcdn = comm - pref

   END SUBROUTINE corpbe
END MODULE m_corpbe
