!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xcvwn
   use m_juDFT
!-----------------------------------------------------------------------
!     Called in case of icorr=4 : spin-polarized exchange-correlation
!       potential from Ceperley-Alder monte carlo results with parametrization
!       of Vosko,Wilk,Nusair, Can. J. Phys. 58, 1200 (1980)
!         ( In case of non-spinpolarization the fit of the correlation
!           energy and potential is essentially equivallent to the
!           parametrization of Perdew,Zunger,PRB 23, 5048 (1981), but
!           in case of spin polarization the magnetization dependent
!           interpolation formula for correlation is better.
!           Part of the subroutine was supplied by M. Manninen. )

!     krla=1: Relativistic correction of exchange energy and potential
!             related to Dirac kinetic energy, according to:
!             A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)

!     be careful: calculation in rydberg!

!     vxc   calculates the XC-potential and
!     exc   calculates the XC-energy

!     based on a subroutine by S. Bluegel;   r.pentcheva 22.01.96
!-----------------------------------------------------------------------

   USE m_constants, ONLY : pi_const
   USE m_relcor
   IMPLICIT NONE

   REAL, PARAMETER, PRIVATE :: cex = 0.91633058742  ! 3/2 * ( 3/(2*pi) )^(2/3)
   REAL, PARAMETER, PRIVATE :: d_15 = 1.e-15
   REAL, PARAMETER, PRIVATE :: one = 1.0 , three = 3.0 , four = 4.0
   REAL, PARAMETER, PRIVATE :: thrd = one/three , two = 2.0
   REAL, PARAMETER, PRIVATE :: fothrd = four * thrd
   REAL, PARAMETER, PRIVATE :: ap = 0.0621814 , xp0 = -0.10498
   REAL, PARAMETER, PRIVATE :: af = 0.0310907 , xf0 = -0.32500
   REAL, PARAMETER, PRIVATE :: al =-0.03377373, xl0 = -0.0047584
   REAL, PARAMETER, PRIVATE :: bp = 3.72744   , cp  = 12.9352
   REAL, PARAMETER, PRIVATE :: bf = 7.06042   , cf  = 18.0578
   REAL, PARAMETER, PRIVATE :: bl = 1.13107   , cl  = 13.0045
   REAL, PARAMETER, PRIVATE :: qp = 6.1519908 , fdd0 = 1.70992093
   REAL, PARAMETER, PRIVATE :: qf = 4.7309269 , ql = 7.123109

CONTAINS
!************************************************************************
   SUBROUTINE vxcvwn( &
      iofile,krla,jspins, &
      mgrid,ngrid,rh, &
      vx,vxc)
!************************************************************************

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points

!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: vxc(mgrid,jspins)     ! xc potential
      REAL, INTENT (OUT) :: vx (mgrid,jspins)     ! x  potential

!     .. Local Scalars ..
      REAL :: s3, s4, c_2, cbrt1, cbrt2, dfds, decdrp, decdrf
      REAL :: dacdr, dbdr, dec1, dec2, cvx, mucp
      REAL :: rho, rh1, rh2 ! total, spin up & spin down  charge density
      REAL :: x, y1, y2, s, thfpi, c_1, rs, beta, bs41, alc
      REAL :: xpx, xfx, xlx, xpx0, xfx0, xlx0, ecp, ecf, fs, ec
      INTEGER :: ir

!     .. Local Arrays ..
      REAL, ALLOCATABLE :: psi(:)       ! relativistic exchange potential corr.

      thfpi  = three / ( four * pi_const )

!-----> evaluate relativistic corrections for exchange

      ALLOCATE ( psi(ngrid) )
      CALL relcor( &
         mgrid,ngrid,jspins,krla, .TRUE. ,rh, &
         psi)

      cvx = fothrd * cex
      IF ( jspins == 2 ) THEN         ! spinpolarized calculation

         c_1 = one / ( two**fothrd - two )
         DO ir = 1,ngrid
            rh1 = max(d_15,rh(ir,1))
            rh2 = max(d_15,rh(ir,jspins))
            rho = rh1 + rh2
            y1 = rh1/rho ; y2 = rh2/rho ; s = y2 - y1  ! s = (rh2 - rh1) / rho
            s3 = s**3 ; s4 = s*s3                      ! spin polarisation
            cbrt1 = (one-s) ** thrd
            cbrt2 = (one+s) ** thrd
            fs = c_1 * ( (one+s)**fothrd + (one-s)**fothrd - two )
            dfds = c_1 * fothrd * ( cbrt1-cbrt2 )
            rs = ( thfpi/rho )**thrd
            x = sqrt(rs)
            xpx = x*x + bp*x + cp
            xfx = x*x + bf*x + cf
            xlx = x*x + bl*x + cl
            xpx0 = xp0*xp0 + bp*xp0 + cp
            xfx0 = xf0*xf0 + bf*xf0 + cf
            xlx0 = xl0*xl0 + bl*xl0 + cl
            ecp  = fec(ap,x,xpx,xp0,xpx0,bp,qp)   ! paramagnetic correlation energy
            ecf  = fec(af,x,xfx,xf0,xfx0,bf,qf)   ! ferromagnetic correlation ener.
            alc  = fec(al,x,xlx,xl0,xlx0,bl,ql)   ! alpha_c
            beta = fdd0* (ecf-ecp)/alc - one      ! beta = (f"(0)*(ecf-ecp))/alc -1
            bs41 = one + beta * s**4

            ec = ecp + alc * fs / fdd0 * bs41        ! total correlation energy
            decdrp = fdedr(rho,ap,x,xpx,xp0,bp,cp)   ! d(ecp)/d(rho)
            decdrf = fdedr(rho,af,x,xfx,xf0,bf,cf)   ! d(ecf)/d(rho)
            dacdr  = fdedr(rho,al,x,xlx,xl0,bl,cl)   ! d(alc)/d(rho)

            dbdr =fdd0* ((decdrf-decdrp)*alc- (ecf-ecp)*dacdr)/alc**2     ! d(beta)/d(rho)
            dec1 =ec + rho* (decdrp+ (dacdr*fs*bs41+alc*fs*dbdr*s4)/fdd0) ! mucp= d(rho*ec)/d(rho)
            dec2 =two*alc/ (fdd0*rho)* (dfds*bs41+four*fs*beta*s3)        ! 2/rho*d(ec)/ds

            c_2 = cvx / rs * psi(ir)                   ! exchange potential muxp=-cvx/rs= 4/3*ex
            vxc(ir,1)     =dec1+dec2*rh2 - c_2*cbrt1   !                        muxp*(2x)**(1/3)
            vxc(ir,jspins)=dec1-dec2*rh1 - c_2*cbrt2   ! calculate exchange correlation potential

            vx (ir,1)     = - c_2*cbrt1
            vx (ir,jspins)= - c_2*cbrt2
         ENDDO

      ELSEIF ( jspins == 1 ) THEN     ! non-spinpolarized calculation

         DO ir = 1,ngrid
            rho = max(d_15,rh(ir,1))
            rs = ( thfpi/rho )**thrd
            x = sqrt(rs)
            xpx = x*x + bp*x + cp
            xpx0 = xp0*xp0 + bp*xp0 + cp

            ecp = fec(ap,x,xpx,xp0,xpx0,bp,qp)     ! correlation energy paramagn.
            decdrp = fdedr(rho,ap,x,xpx,xp0,bp,cp) ! d(ecp)/d(rho)
            mucp = ecp + rho*decdrp                ! d(rho*ec)/d(rho)

            vxc(ir,1)     = mucp - cvx/rs*psi(ir)  ! -1.221774/rs :exchange potential

            vx (ir,1)     = - cvx/rs*psi(ir)
         ENDDO

      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
         CALL juDFT_error("vxcvwn",calledby="xcvwn")
      ENDIF

      DEALLOCATE (psi)
      RETURN

   END SUBROUTINE vxcvwn
!***********************************************************************
   SUBROUTINE excvwn( &
      iofile,krla,jspins, &
      mgrid,ngrid,rh, &
      exc)
!***********************************************************************

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points

!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy

!     .. Local Scalars ..
      REAL ::  c_1,ex
      REAL :: rho, rh1, rh2     ! total, spin up & spin down  charge density
      REAL :: x, y1, y2, s, thfpi, rs, beta, bs41, alc
      REAL :: xpx, xfx, xlx, xpx0, xfx0, xlx0, ecp, ecf, fs, ec
      INTEGER :: ir

!     .. Local Arrays ..
      REAL, ALLOCATABLE :: phi(:)       ! relativistic exchange energy correct.

      thfpi  = three / ( four * pi_const )

      ALLOCATE ( phi(ngrid) )
      CALL relcor( &
         mgrid,ngrid,jspins,krla, .FALSE. ,rh, &
         phi)

      IF ( jspins == 2 ) THEN         ! spinpolarized calculation

         c_1 = one / ( two**fothrd - two )
         DO ir = 1,ngrid
            rh1 = max(d_15,rh(ir,1))
            rh2 = max(d_15,rh(ir,jspins))
            rho = rh1 + rh2
            y1 = rh1/rho ; y2 = rh2/rho ; s = y2 - y1
            fs = c_1 * ( (one+s)**fothrd + (one-s)**fothrd - two )
            rs = ( thfpi/rho )**thrd
            x = sqrt(rs)
            xpx = x*x + bp*x + cp
            xfx = x*x + bf*x + cf
            xlx = x*x + bl*x + cl
            xpx0 = xp0*xp0 + bp*xp0 + cp
            xfx0 = xf0*xf0 + bf*xf0 + cf
            xlx0 = xl0*xl0 + bl*xl0 + cl
            ecp  = fec(ap,x,xpx,xp0,xpx0,bp,qp)   ! paramagnetic correlation energy
            ecf  = fec(af,x,xfx,xf0,xfx0,bf,qf)   ! ferromagnetic correlation ener.
            alc  = fec(al,x,xlx,xl0,xlx0,bl,ql)   ! alpha_c
            beta = fdd0* (ecf-ecp)/alc - one      ! beta = (f"(0)*(ecf-ecp))/alc -1
            bs41 = one + beta * s**4

            ec = ecp + alc * fs / fdd0 * bs41        ! total correlation energy
            ex = -cex/rs* (one + 0.2599210482 * fs)  ! ex = exp + (exf-exp)*f(s)
            ! exf - exp = (2**(1/3)-1) * exp
            !             = 0.25992105 * exp
            exc(ir) = ec + ex*phi(ir)
         ENDDO

      ELSEIF ( jspins == 1 ) THEN     ! non-spinpolarized calculation

         DO ir = 1,ngrid
            rho = max(d_15,rh(ir,1))
            rs = ( thfpi/rho )**thrd
            x = sqrt(rs)
            xpx = x*x + bp*x + cp
            xpx0 = xp0*xp0 + bp*xp0 + cp
            ecp = fec(ap,x,xpx,xp0,xpx0,bp,qp)  ! correlation energy paramagn.
            ex  = -cex / rs                     ! exp: paramag. exchange energy
            ! (like in wigner formula)
            exc(ir) = ecp + ex*phi(ir)
         ENDDO

      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
         CALL juDFT_error("excvwn",calledby="xcvwn")
      ENDIF

      DEALLOCATE (phi)
      RETURN

   END SUBROUTINE excvwn

!--------------------------------------------------------------------
   REAL FUNCTION fec(ai,z,ziz,zi0,ziz0,bi,qi)
      REAL :: ai,z,ziz,zi0,ziz0,bi,qi
      fec = ai* (alog(z*z/ziz)+ two*bi/qi * atan(qi/(two*z+bi)) - &
                 bi*zi0/ziz0* ( alog((z-zi0)**2/ziz) + &
                               two*(bi+ two*zi0)/qi*atan(qi/ (two*z+bi))) )
   END  FUNCTION fec
   REAL FUNCTION fdedr(rho,ai,z,ziz,zi0,bi,ci)
      REAL :: rho,ai,z,ziz,zi0,bi,ci
      fdedr = -ai*z/(three*rho*ziz)*(ci/z-bi*zi0/(z-zi0))
   END FUNCTION fdedr
!--------------------------------------------------------------------

END MODULE m_xcvwn
