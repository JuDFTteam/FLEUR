MODULE m_exchpbe
!----------------------------------------------------------------------
!     pbe exchange for a spin-unpolarized electronic system
!     k burke's modification of pw91 codes, may 14, 1996
!     modified again by k. burke, june 29, 1996, with simpler fx(s)
!     inclusion of HSE function by M. Schlipf 2009
!----------------------------------------------------------------------
!     references:
!     [a]j.p.~perdew, k.~burke, and m.~ernzerhof, submiited to prl, may96
!     [b]j.p. perdew and y. wang, phys. rev.  b {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (e).
!     [c] B.~Hammer, L.~B.~Hansen and J.~K.~Norskov PRB 59 7413 (1999)
!     [d] J.~Heyd, G.~E.~Scuseria, M.~Ernzerhof, J. Chem. Phys. {\bf 118},
!     8207 (2003)
!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE exchpbe(xcpot,rho,s,u,v,lgga,lpot, &
                      ex,vx,vx_sr)
      USE m_hsefunctional, ONLY: calculateEnhancementFactor
      USE m_constants,     ONLY: pi_const
      USE m_types_xcpot_data
      USE m_judft
      IMPLICIT NONE

!     .. Arguments
      TYPE(t_xcpot_data), INTENT (IN)  :: xcpot
      INTEGER, INTENT (IN)  :: lgga ! =0=>don't put in gradient corrections, just lda
      INTEGER, INTENT (IN)  :: lpot ! =0=>don't get potential and don't need u and v
      REAL,    INTENT (IN)  :: rho ! density
      REAL,    INTENT (IN)  :: s ! abs(grad rho)/(2*kf*rho), where kf=(3 pi^2 rho)^(1/3)
      REAL,    INTENT (IN)  :: u ! (grad rho)*grad(abs(grad rho))/(rho**2 * (2*kf)**3)
      REAL,    INTENT (IN)  :: v ! (laplacian rho)/(rho*(2*kf)**2) [pw86(24)]
      REAL,    INTENT (OUT) :: ex,vx ! exchange energy per electron (ex) and potential (vx)

!     new variables for the HSE functional
      REAL,    INTENT (OUT) :: vx_sr ! short ranged part of potential
      REAL :: kF,factor,fxhse
      REAL :: dFx_ds,d2Fx_ds2,dFx_dkF,d2Fx_dsdkF

!     .. local variables ..
      REAL :: ul,exunif,fs,fss,fxpbe,p0,s2
      REAL :: xwu,css,dxwu,ddx  ! for wu-cohen
      REAL, PARAMETER :: teo = 10.e0/81.e0 ! for wu-cohen
      REAL, PARAMETER :: cwu = 0.0079325 ! for wu-cohen
      REAL, PARAMETER :: thrd=1.e0/3.e0
      REAL, PARAMETER :: thrd4=4.e0/3.e0
      REAL, PARAMETER :: ax=-0.738558766382022405884230032680836e0 ! -0.75*(3/pi)^(1/3)

!----------------------------------------------------------------------
!     uk, ul defined after [a](13)  (for icorr==7)
!----------------------------------------------------------------------
!      IF (xcpot%is_name("pbe").OR.xcpot%is_name("wc").OR.
!     +     xcpot%is_name("PBEs")) THEN
!         uk=0.8040
!      ELSEIF (xcpot%is_name("rpbe")) THEN
!         uk=1.2450
!      ELSEIF (xcpot%is_name("Rpbe")) THEN  ! changed to [c]
!         uk=0.8040
!      ELSEIF (xcpot%is_name("pbe0").OR.xcpot%is_name("hse")
!     +        .OR. xcpot%is_name("lhse").OR.xcpot%is_name("vhse")) THEN
!         uk=0.8040
!      ELSE
!         CALL judft_error("BUG:Wrong xcpot",calledby="exchpbe.F")
!      ENDIF
!      IF (xcpot%is_name("PBEs")) THEN     ! pbe_sol
!         um=0.123456790123456d0
!      ELSE
!         um=0.2195149727645171e0
!      ENDIF

      ul=xcpot%um/xcpot%uk
!----------------------------------------------------------------------
!     construct lda exchange energy density:
!     e_x[unif] = -0.75 * (3/pi)^(1/3) * rho^(4/3)
!----------------------------------------------------------------------
      exunif = ax*rho**thrd
      IF (lgga == 0) THEN
         ex = exunif
         vx = ex*thrd4
         RETURN
      ENDIF
!----------------------------------------------------------------------
!     construct pbe enhancement factor
!     e_x[pbe]=e_x[unif]*fxpbe(s)
!     fxpbe(s)=1+xcpot%uk-xcpot%uk/(1+ul*s*s)                 [a](13)
!----------------------------------------------------------------------
      s2 = s*s

!     calculate fxpbe
      p0 = 1.e0 + ul*s2
      fxpbe = 1e0 + xcpot%uk - xcpot%uk/p0
      IF (xcpot%is_Rpbe) THEN
         p0 = exp( - ul*s2 )
         fxpbe = 1e0 + xcpot%uk * ( 1e0 - p0 )
      ELSEIF (xcpot%is_wc) THEN
         css = 1+cwu*s2*s2
         IF(s2.GT.100) THEN ! This is introduced because the PGI compiler has problems with calculating exp(-somethingLarge).
            xwu = teo*s2 + log(css)
         ELSE
            xwu = teo*s2 + (xcpot%um-teo)*s2*exp(-s2) + log(css)
         END IF
         p0 = 1.e0 + xwu/xcpot%uk
         fxpbe = 1e0 + xcpot%uk - xcpot%uk/p0
      ENDIF
!     -gu
!     Mixing of short and long range components
      IF (xcpot%is_hse) THEN
         !     ex = (1-a)ex,SR + ex,LR
         !     = (1-a)ex,SR + (ex,PBE - ex,SR)
         !     = (fxpbe - a*fxhse)*exunif
         !     Calculate the enhancement factor fxhse and its derivatives
         !     as integral over the exchange hole (cf. [d])
         kF = (3.0 * pi_const**2 * rho)**thrd
         !! CALL judft_error("HSE not implemented",calledby="exchpbe")
         ! his creates a depency loop
         CALL calculateEnhancementFactor(kF, s, fxhse, dFx_ds, d2Fx_ds2, &
                                         dFx_dkF, d2Fx_dsdkF)
         ex = exunif * (fxpbe - xcpot%exchange_weight * fxhse )
      ELSE
         ex = exunif*fxpbe
      END IF
      IF (lpot == 0) RETURN
!----------------------------------------------------------------------
!     energy done. now the potential:
!     find first and second derivatives of fx w.r.t s.
!     fs=(1/s)*d fxpbe/ ds
!     fss=d fs/ds
!----------------------------------------------------------------------
!     derivatives for the pbe part
      fs = 2.e0*xcpot%uk*ul/ (p0*p0)
      fss = -4.e0*ul*s*fs/p0
      IF (xcpot%is_Rpbe) THEN
         fs = 2.e0*ul*p0
         fss = -2.e0*ul*s*fs
      ELSEIF (xcpot%is_wc) THEN
         IF (s2.GT.100) THEN ! This is introduced because the PGI compiler has problems with calculating exp(-somethingLarge).
            dxwu = 2*teo + 4*cwu*s2/css
            fs = dxwu / (p0*p0)
            ddx = 4*s*(2*cwu*(1-cwu*s2*s2)/css**2)
         ELSE
            dxwu = 2*teo + 2*(xcpot%um-teo)*exp(-s2)*(1-s2) + 4*cwu*s2/css
            fs = dxwu / (p0*p0)
            ddx = 4*s*((xcpot%um-teo)*exp(-s2)*(s2-2)+2*cwu* &
                       (1-cwu*s2*s2)/css**2)
         END IF
         fss = ( ddx - 2*s*dxwu*dxwu/(p0*xcpot%uk) ) / (p0*p0)
      ENDIF
!     -gu
!----------------------------------------------------------------------
!     calculate potential from [b](24)
!----------------------------------------------------------------------
      vx = exunif* (thrd4*fxpbe- (u-thrd4*s2*s)*fss-v*fs)

      IF ( .NOT. (xcpot%is_hse)) RETURN

!----------------------------------------------------------------------
!     short ranged potential (HSE functional)
!----------------------------------------------------------------------
!     calculate fs and fss for the HSE functional
!     where the 1st and 2nd derivative of Fx are known
      fs    = dFx_ds / s
      fss   = (d2Fx_ds2 - fs) / s
      vx_sr = exunif * ( thrd4*fxhse - (u-thrd4*s2*s)*fss - v*fs &
                        + thrd*kF * ( dFx_dkF - d2Fx_dsdkF*s ) )

   END SUBROUTINE exchpbe
END MODULE m_exchpbe
