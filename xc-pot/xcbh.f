!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_xcbh 
      use m_juDFT
!-----------------------------------------------------------------------
!     Called in case of icorr=2,3 : spin-polarized exchange-correlation 
!       potential of U. von Barth and L. Hedin, J.Phys.C5,1629 (1972)
!       icorr = 2: parametrization of Moruzzi,Janak,Williams
!       icorr = 3: parametrization of von Barth and Hedin
!
!     krla=1: Relativistic correction of exchange energy and potential 
!             related to Dirac kinetic energy, according to:
!             A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)
!
!     be careful: calculation in rydberg!
!
!     vxcbh calculates the XC-potential and
!     excbh calculates the XC-energy
!
!     based on a subroutine by S. Bluegel;   r.pentcheva 22.01.96
!-----------------------------------------------------------------------

      USE m_constants, ONLY : pi_const
      USE m_relcor
      IMPLICIT NONE

      REAL, PARAMETER, PRIVATE :: ff  = 3.847322101863  ! 1 / ( 2^(1/3) - 1 )
      REAL, PARAMETER, PRIVATE :: cvx = 1.221774115422  ! 2 * ( 3/(2*pi) )^(2/3)
      REAL, PARAMETER, PRIVATE :: cpmjw = 0.045  , cfmjw = 0.0225
      REAL, PARAMETER, PRIVATE :: cpvbh = 0.0504 , cfvbh = 0.0254
      REAL, PARAMETER, PRIVATE :: rpmjw = 21.0 , rfmjw = 52.916684096
      REAL, PARAMETER, PRIVATE :: rpvbh = 30.0 , rfvbh = 75.0
      REAL, PARAMETER, PRIVATE :: d_15 = 1.e-15
      REAL, PARAMETER, PRIVATE :: one = 1.0 , three = 3.0 , four = 4.0
      REAL, PARAMETER, PRIVATE :: half = 0.5 , thrd = one/three
      REAL, PARAMETER, PRIVATE :: hfthrd = 0.79370052705 ! 2^(-1/3)
      REAL, PARAMETER, PRIVATE :: thrhalf = three * half
      REAL, PARAMETER, PRIVATE :: fothrd = four * thrd , two = 2.0


      CONTAINS
!************************************************************************
      SUBROUTINE vxcbh
     >                (iofile,xcpot,jspins,
     >                 mgrid,ngrid,rh,
     <                 vx,vxc)
!************************************************************************
      USE m_types
!     
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      TYPE(t_xcpot), INTENT (IN) :: xcpot  !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: vxc(mgrid,jspins)     ! xc potential
      REAL, INTENT (OUT) :: vx (mgrid,jspins)     ! x  potential
!
!     .. Local Scalars ..
      REAL txthrd,tythrd,muxp,mucp,mucf,ecfmp,tauc,mucnm
      REAL :: rho, rh1, rh2 ! total, spin up & spin down charge density
      REAL :: x, y, cp, cf, rp, rf, rs, ecprs, ecfrs
      INTEGER :: ir
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: psi(:)       ! relativistic exchange potential corr.
!
!---- Intrinsic Functions
      INTRINSIC alog,max

!
!-----> evaluate relativistic corrections for exchange
! 
      ALLOCATE ( psi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,xcpot%krla,.true.,rh,
     <            psi)
!
!-----> select exchange correlation potential
!
      IF (xcpot%is_name("mjw")) THEN
         cp = cpmjw ; cf = cfmjw
         rp = rpmjw ; rf = rfmjw
      ELSEIF (xcpot%is_name("bh")) THEN
         cp = cpvbh ; cf = cfvbh
         rp = rpvbh ; rf = rfvbh
      ELSE
         WRITE (iofile,FMT=2000)
          CALL juDFT_error("BUG:vxcbh",calledby="xcbh")
      END IF
 2000 FORMAT (13x,'set key for exchange-correlation potential')
!
!-----> calculate exchange correlation potential
!
      IF ( jspins .EQ. 2) THEN               ! spinpolarized calculation

        DO ir = 1,ngrid                        ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1))
          rh2 = max(d_15,rh(ir,jspins))
          rho = rh1 + rh2
          x = rh1/rho
          y = rh2/rho
          txthrd = (2*x)**thrd
          tythrd = (2*y)**thrd
          rs= (four*pi_const*rho/three)**thrd
          rs = 1/rs

          ecprs = -cp*fc(rs/rp)            ! calculate correlation energy
          ecfrs = -cf*fc(rs/rf)            ! p : paramagnetic, f : ferromagnetic
                                           ! x : exchange,     c : correlation
          muxp = -psi(ir)* (cvx/rs)        ! paramagnetic exchange potential 
                                           !       (psi contains rel. corr.)
          mucp = -cp*alog(one+rp/rs)       ! calculate correlation potential
          mucf = -cf*alog(one+rf/rs)
          ecfmp = fothrd * (ecfrs-ecprs)
          tauc = mucf - mucp - ecfmp
          mucnm = mucp + tauc*fex(x) - ff*ecfmp

          vxc(ir,1)      = mucnm + (muxp+ff*ecfmp)*txthrd  ! collect correlation 
          vxc(ir,jspins) = mucnm + (muxp+ff*ecfmp)*tythrd  ! and exchange parts
        
          vx (ir,1)      = muxp*txthrd
          vx (ir,jspins) = muxp*tythrd
        ENDDO

      ELSEIF ( jspins .EQ. 1 ) THEN        ! non - spinpolarized calculation

        DO ir = 1, ngrid                   ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1))
          rs = (four*pi_const*rh1/three)**thrd
          rs = 1/rs 
          muxp = -psi(ir) * (cvx/rs)       ! paramagnetic exchange potential 
                                           !       (psi contains rel. corr.)
          mucp = -cp* alog(one+rp/rs)      ! calculate correlation potential
          vxc(ir,1)     = mucp + muxp      ! collect correlation & exchange part
        
          vx (ir,1)     = muxp
        ENDDO

      ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
          CALL juDFT_error("vxcbh",calledby="xcbh")
      ENDIF

      DEALLOCATE (psi)
      RETURN   

      END SUBROUTINE vxcbh
C***********************************************************************
      SUBROUTINE excbh
     >                (iofile,xcpot,jspins,
     >                 mgrid,ngrid,rh,
     <                 exc)
C***********************************************************************
      USE m_types
!     
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      TYPE(t_xcpot), INTENT (IN) :: xcpot !  run mode parameters
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy
!
!     .. Local Scalars ..
      REAL thfpi,thrquart,exprs,exfrs,excprs,excfrs
      REAL :: rho, rh1, rh2 ! total, spin up & spin down charge density
      REAL :: x, y, cp, cf, rp, rf, rs, ecprs, ecfrs
      INTEGER :: ir
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: phi(:)       ! relativistic exchange energy correct.
!
!-----> Intrinsic Functions
      INTRINSIC alog,max
!
      thrquart = 0.75
      thfpi = thrquart/pi_const

      ALLOCATE ( phi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,xcpot%krla,.false.,rh,
     <            phi)
!
!-----> select exchange correlation potential
!
      IF (xcpot%is_name("mjw")) THEN
         cp = cpmjw ; cf = cfmjw
         rp = rpmjw ; rf = rfmjw
      ELSEIF (xcpot%is_name("bh")) THEN
         cp = cpvbh ; cf = cfvbh
         rp = rpvbh ; rf = rfvbh
      ELSE
         WRITE (iofile,FMT=2000)
          CALL juDFT_error("excbh",calledby="xcbh")
      END IF
 2000 FORMAT (13x,'set key for exchange-correlation potential')

      IF ( jspins .EQ. 2) THEN       ! spinpolarized calculation

        DO  ir = 1,ngrid                  ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1))
          rh2 = max(d_15,rh(ir,jspins))
          rho = rh1 + rh2
          x = rh1/rho
          rs= (thfpi/rho)**thrd

          exprs = -phi(ir)*thrquart*cvx/rs    ! first exchange energy 
          exfrs = (2.0**thrd)*exprs           ! phi contains rel. corr.

          ecprs = -cp*fc(rs/rp)               ! calculate correlation energy
          ecfrs = -cf*fc(rs/rf)               ! p: paramagnetic, f: ferromagn.

          excprs = exprs + ecprs              ! now add correlation energy
          excfrs = exfrs + ecfrs

          exc(ir) = excprs + (excfrs-excprs)*fex(x) ! collect all terms
        ENDDO 

      ELSEIF ( jspins .EQ. 1 ) THEN  ! non - spinpolarized calculation

        DO ir = 1,ngrid                    ! loop over realspace gridpoints
          rh1 = max(d_15,rh(ir,1)) 
          rs = (thfpi/rh1)**thrd
          exprs = -phi(ir)*thrquart*cvx/rs ! exchange energy ; phi contains 
                                           ! relativistic correctionS
          ecprs = -cp*fc(rs/rp)            ! calculate correlation energy
          exc(ir) = exprs + ecprs          ! add correlation energy
        ENDDO

       ELSE
         WRITE (iofile,'('' error in jspins, jspins ='',i2)') jspins
          CALL juDFT_error("excbh",calledby="xcbh")
       ENDIF

      DEALLOCATE (phi)
      RETURN

      END SUBROUTINE excbh
!--------------------------------------------------------------------
      REAL FUNCTION fc(x)
        REAL x
        fc = (one+(x)*(x)*(x))*alog(one+one/(x))
     +        + half*(x) - (x)*(x) - thrd
      END  FUNCTION fc
      REAL FUNCTION fex(x)
        REAL x
        fex = ff/hfthrd*((x)**fothrd +(1-(x))**fothrd - hfthrd)
      END FUNCTION fex
!--------------------------------------------------------------------

      END MODULE m_xcbh
