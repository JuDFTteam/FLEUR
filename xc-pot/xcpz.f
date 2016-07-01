      MODULE m_xcpz
      use m_juDFT
!-----------------------------------------------------------------------
!     Called in case of icorr=5 : spin-polarized exchange-correlation 
!        from Ceperley-Alder monte carlo results with parametrization
!        of Perdew,Zunger,PRB 23, 5048 (1981)
!
!     krla=1: Relativistic correction of exchange energy and potential 
!             related to Dirac kinetic energy, according to:
!             A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)
!
!     be careful: calculation in rydberg!
!
!     vxcpz   calculates the XC-potential and
!     excpz   calculates the XC-energy
!
!     based on a subroutine by S. Bluegel;   r.pentcheva 22.01.96
!-----------------------------------------------------------------------

      USE m_constants, ONLY : pi_const
      USE m_relcor
      IMPLICIT NONE

      REAL, PARAMETER, PRIVATE :: cvx = 1.221774115422  ! 2 * ( 3/(2*pi) )^(2/3)
      REAL, PARAMETER, PRIVATE :: d_15 = 1.e-15 , c76 = 7.0 / 6.0
      REAL, PARAMETER, PRIVATE :: one = 1.0 , three = 3.0 , four = 4.0
      REAL, PARAMETER, PRIVATE :: half = 0.5 , thrd = one/three
      REAL, PARAMETER, PRIVATE :: thrhalf = three * half , two = 2.0
      REAL, PARAMETER, PRIVATE :: c23 = two * thrd , c43 = four * thrd
      REAL, PARAMETER, PRIVATE :: ap  =  0.0622 , af  =  0.0311
      REAL, PARAMETER, PRIVATE :: bp  = -0.0960 , bf  = -0.0538
      REAL, PARAMETER, PRIVATE :: cp  =  0.0040 , cf  =  0.0014
      REAL, PARAMETER, PRIVATE :: dp  = -0.0232 , df  = -0.0096
      REAL, PARAMETER, PRIVATE :: gp  = -0.2846 , gf  = -0.1686
      REAL, PARAMETER, PRIVATE :: b1p =  1.0529 , b1f =  1.3981
      REAL, PARAMETER, PRIVATE :: b2p =  0.3334 , b2f =  0.2611


      CONTAINS
!************************************************************************
      SUBROUTINE vxcpz(
     >                 iofile,krla,jspins,
     >                 mgrid,ngrid,rh,
     <                 vx,vxc)
!************************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)                   ! charge density
      REAL, INTENT (OUT) :: vx(mgrid,jspins),vxc(mgrid,jspins) ! x/xc potential
!
!     .. Local Scalars ..
      REAL c_2, cbrt1, cbrt2, dfds, decds, dec1, dec2, vcf, vcp
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: psi(:)       ! relativistic exchange potential corr.
!
!-----s Intrinsic Functions
      INTRINSIC max
      REAL :: rho, rh1, rh2 ! total, spin up & spin down  charge density
      REAL :: fothrd, thfpi, c_1, y1, y2, s, fs, rs
      REAL :: ecp, ecf
      INTEGER :: ir

      fothrd = c43
      thfpi  = three / ( four * pi_const )
!
!-----> evaluate relativistic corrections for exchange
!
      ALLOCATE ( psi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.true.,rh,
     <            psi)

      IF ( jspins.EQ.2 ) THEN         ! spinpolarized calculation

        c_1 = one / ( two**fothrd - two )
        DO ir = 1,ngrid
           rh1 = max(d_15,rh(ir,1))
           rh2 = max(d_15,rh(ir,jspins))
           rho = rh1 + rh2
           y1 = rh1/rho ; y2 = rh2/rho ; s = y2 - y1  ! s = (rh2 - rh1) / rho
           cbrt1 = (one-s) ** thrd
           cbrt2 = (one+s) ** thrd
           fs = c_1 * ( (one+s)**fothrd + (one-s)**fothrd - two )
           dfds = c_1 *  fothrd * (cbrt1 - cbrt2) 
           rs = ( thfpi/rho )**thrd
         
           IF (rs.GE.one) THEN
              ecp = fecl(rs,gp,b1p,b2p)  ! correlation energy paramagnetic
              ecf = fecl(rs,gf,b1f,b2f)  ! correlation energy ferromagnetic
              vcp = fvcl(ecp,rs,b1p,b2p) ! d(rho*ecp)/d(rho) = ecp - rs/3*d(ecp)/d(rs)
              vcf = fvcl(ecf,rs,b1f,b2f) ! d(rho*ecf)/d(rho)
           ELSE
              ecp = fecs(rs,ap,bp,cp,dp)
              ecf = fecs(rs,af,bf,cf,df)
              vcp = fvcs(rs,ap,bp,cp,dp)
              vcf = fvcs(rs,af,bf,cf,df)
           ENDIF
           decds = (ecf-ecp)*dfds     ! = d(ec)/d(rho)
           dec1  = vcp + (vcf-vcp)*fs ! = d(rho*ec)/d(rho)
           dec2  = two/rho*decds      ! = 2/rho*d(ec)/ds

           c_2 = cvx / rs * psi(ir)                   ! exchange potential muxp=-cvx/rs= 4/3*ex
           vxc(ir,1)     =dec1+dec2*rh2 - c_2*cbrt1   !                        muxp*(2x)**(1/3)
           vxc(ir,jspins)=dec1-dec2*rh1 - c_2*cbrt2   ! calculate exchange correlation potential
                                         ! vc = ec +vcp + (vcf-vcp)f(s) - (ecf-ecp)df/ds*(s+/-1)
        
           vx (ir,1)     = - c_2*cbrt1 
           vx (ir,jspins)= - c_2*cbrt2 
        ENDDO

      ELSEIF ( jspins.EQ.1 ) THEN     ! non-spinpolarized calculation

        DO ir = 1,ngrid
           rho = max(d_15,rh(ir,1))
           rs = ( thfpi/rho )**thrd
           IF (rs.GE.one) THEN
              ecp = fecl(rs,gp,b1p,b2p)
              vcp = fvcl(ecp,rs,b1p,b2p)
           ELSE
              ecp = fecs(rs,ap,bp,cp,dp)
              vcp = fvcs(rs,ap,bp,cp,dp)
           ENDIF
           vxc(ir,1) = vcp - cvx / rs * psi(ir)
           
           vx (ir,1) = cvx / rs * psi(ir)
        ENDDO

      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
          CALL juDFT_error("vxcpz",calledby="xcpz")
      ENDIF

      DEALLOCATE (psi)
      RETURN

      END SUBROUTINE vxcpz
C***********************************************************************
      SUBROUTINE excpz(
     >                 iofile,krla,jspins,
     >                 mgrid,ngrid,rh,
     <                 exc)
C***********************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy
!
!     .. Local Scalars ..
      REAL ec, ex, cex
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: phi(:)       ! relativistic exchange energy correct.
!
!-----> Intrinsic Functions
      INTRINSIC max
!
      REAL :: rho, rh1, rh2 ! total, spin up & spin down  charge density
      REAL :: fothrd, thfpi, c_1, y1, y2, s, fs, rs
      REAL :: ecp, ecf
      INTEGER :: ir

      fothrd = c43
      thfpi  = three / ( four * pi_const )
      cex = cvx / c43

      ALLOCATE ( phi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.false.,rh,
     <            phi)

      IF ( jspins.EQ.2 ) THEN         ! spinpolarized calculation

        c_1 = one / ( two**fothrd - two )
        DO ir = 1,ngrid
           rh1 = max(d_15,rh(ir,1))
           rh2 = max(d_15,rh(ir,jspins))
           rho = rh1 + rh2
           rs = ( thfpi/rho )**thrd
           y1 = rh1/rho ; y2 = rh2/rho ; s = y2 - y1  ! s = (rh2 - rh1) / rho
           fs = c_1 * ( (one+s)**fothrd + (one-s)**fothrd - two )
           IF (rs.GE.one) THEN
              ecp = fecl(rs,gp,b1p,b2p)  ! correlation energy paramagnetic
              ecf = fecl(rs,gf,b1f,b2f)  ! correlation energy ferromagnetic
           ELSE
              ecp = fecs(rs,ap,bp,cp,dp)
              ecf = fecs(rs,af,bf,cf,df)
           ENDIF

           ec = ecp + (ecf-ecp)*fs                ! total correlation energy
           ex = -cex/rs* (one + 0.2599210482*fs)  ! ex = exp + (exf-exp)*f(s)
                                                  ! exf-exp = (2**(1/3)-1) * exp 
           exc(ir) = ec + ex*phi(ir)        
        ENDDO

      ELSEIF ( jspins.EQ.1 ) THEN     ! non-spinpolarized calculation

        DO ir = 1,ngrid
           rho = max(d_15,rh(ir,1))
           rs = ( thfpi/rho )**thrd
           IF (rs.GE.one) THEN
              ecp = fecl(rs,gp,b1p,b2p)  
           ELSE
              ecp = fecs(rs,ap,bp,cp,dp)
           ENDIF
           ex = -cex/rs
           exc(ir) = ecp + ex*phi(ir)
        ENDDO

      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
          CALL juDFT_error("vxcpz",calledby="xcpz")
      ENDIF

      DEALLOCATE (phi)
      RETURN

      END SUBROUTINE excpz

!--------------------------------------------------------------------
      REAL FUNCTION fecl(r,g,b1,b2)
        REAL r,g,b1,b2
        fecl = g / ( one + b1*sqrt(r) + b2*r )
      END  FUNCTION fecl
      REAL FUNCTION fvcl(ce,r,b1,b2)
        REAL ce,r,b1,b2
        fvcl = ce* (one+c76*b1*sqrt(r)+c43*b2*r)/(one+b1*sqrt(r)+b2*r)
      END FUNCTION fvcl
      REAL FUNCTION fecs(r,a,b,c,d)
        REAL r,a,b,c,d
        INTRINSIC alog
        fecs = a*alog(r) + b + c*r*alog(r) + d*r
      END FUNCTION fecs
      REAL FUNCTION fvcs(r,a,b,c,d)
        REAL r,a,b,c,d
        INTRINSIC alog
        fvcs = a*alog(r) + (b-a/three) + c23*c*r*alog(r) +
     +                    (two*d-c)*r/three
      END FUNCTION fvcs
!--------------------------------------------------------------------

      END MODULE m_xcpz
