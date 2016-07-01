!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_xcxal
      use m_juDFT
!-----------------------------------------------------------------------
!     Called in case of icorr=0 : nonspin-polarized exchange-correlation 
!       potential with the X-alpha method
!
!     We caclulate  V_x = alpha V_{x,Slater}  where alpha = 2/3 
!     i.e. the output corresponds to the Gaspar result for the
!     local approximation to the Hartree-Fock method.
!
!     krla=1: Relativistic correction of exchange energy and potential 
!             related to Dirac kinetic energy, according to:
!             A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)
!
!     be careful: calculation in rydberg!
!
!     vxcxal   calculates the XC-potential and
!     excxal   calculates the XC-energy
!
!     based on a subroutine by S. Bluegel;   r.pentcheva 22.01.96
!-----------------------------------------------------------------------

      USE m_constants, ONLY : pi_const
      USE m_relcor
      IMPLICIT NONE

      REAL, PARAMETER, PRIVATE :: one = 1.0 , three = 3.0 , four = 4.0
      REAL, PARAMETER, PRIVATE :: thrd = one/three , d_15 = 1.e-15 
      REAL, PARAMETER, PRIVATE :: cvx = 1.221774115422  ! 2 * ( 3/(2*pi) )^(2/3)

      REAL, PRIVATE ::  rs, rho, thfpi
      INTEGER, PRIVATE :: ir

      CONTAINS
!************************************************************************
      SUBROUTINE vxcxal(
     >                  krla,jspins,
     >                  mgrid,ngrid,rh,
     <                  vx,vxc)
!************************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspins
      INTEGER, INTENT (IN) :: krla        !  switch for rel. corr.
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: vx (mgrid,jspins)
      REAL, INTENT (OUT) :: vxc(mgrid,jspins)     ! xc potential
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: psi(:)       ! relativistic exchange potential corr.
!
!-----s Intrinsic Functions
      INTRINSIC alog,max

      thfpi  = three / ( four * pi_const )
!
!-----> evaluate relativistic corrections for exchange
!
      ALLOCATE ( psi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.true.,rh,
     <            psi)

      DO ir = 1,ngrid
         IF (jspins.EQ.1) THEN
           rho = max(d_15,rh(ir,1))
         ELSE
           rho = max(d_15,rh(ir,1)) + max(d_15,rh(ir,jspins))
         ENDIF
         rs = (thfpi/rho)**thrd
                                       ! exc:  exchange energy = -0.9163306/rs
         vxc(ir,1) = - psi(ir)*cvx/rs  ! paramagnetic exchange potential = 4/3 exc
         
         vx (ir,1) = - psi(ir)*cvx/rs
      
      ENDDO
      IF ( jspins .EQ. 2 ) THEN        ! spinpolarized calculation
         DO ir = 1,ngrid
            vxc(ir,jspins) = vxc(ir,1)
            
            vx(ir,jspins)  = vx(ir,1)
         ENDDO
      ENDIF

      DEALLOCATE (psi)
      RETURN

      END SUBROUTINE vxcxal
C***********************************************************************
      SUBROUTINE excxal(
     >                  iofile,krla,jspins,
     >                  mgrid,ngrid,rh,
     <                  exc)
C***********************************************************************
!
!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: krla, jspins
      INTEGER, INTENT (IN) :: iofile      !  file number for read and write
      INTEGER, INTENT (IN) :: mgrid,ngrid !  mesh points
!
!     .. Array Arguments ..
      REAL, INTENT (IN)  :: rh(mgrid,jspins)      ! charge density
      REAL, INTENT (OUT) :: exc(mgrid)            ! xc energy
!
!     .. Local Scalars ..
      REAL cex
!
!     .. Local Arrays ..
      REAL, ALLOCATABLE :: phi(:)       ! relativistic exchange energy correct.
!
!-----> Intrinsic Functions
      INTRINSIC alog,max
!
      cex = three * cvx / four
      thfpi  = three / ( four * pi_const )

      ALLOCATE ( phi(ngrid) )
      CALL relcor(
     >            mgrid,ngrid,jspins,krla,.false.,rh,
     <            phi)

      DO ir = 1,ngrid
        IF (jspins.EQ.1) THEN
          rho = max(d_15,rh(ir,1))
        ELSE
          rho = max(d_15,rh(ir,1)) + max(d_15,rh(ir,jspins))
        ENDIF
        rs = (thfpi/rho)**thrd
                                       ! exc:  exchange energy = -0.9163306/rs
        exc(ir)     = - phi(ir)*cex / rs
      ENDDO
      IF ( jspins.EQ.2 ) WRITE (iofile,*)
     +                'WARNING: X-alpha XC & spin-polarized calculation'
      IF ( jspins.GT.2 ) THEN
          WRITE (iofile,'('' or error in jspins, jspins ='',i2)') jspins
           CALL juDFT_error("excxal",calledby="xcxal")
      END IF

      DEALLOCATE (phi)
      RETURN

      END SUBROUTINE excxal
      END MODULE m_xcxal
