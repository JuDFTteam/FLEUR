!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_dwigner
      use m_juDFT

! Calculate the Wigner rotation matrices for complex spherical
! harmonics for all space-group rotations and l=1,2,3. Needed 
! for the calculation of the density matrix in nmat.
!                 gb 2002

c******************************************************
c     interface for real and integer rotation matrices
c       FF, Oct 2006
c*****************************************************
      PRIVATE
      INTERFACE d_wigner
      MODULE PROCEDURE real_wigner, integer_wigner, integer_wigner_1op
      END INTERFACE
      LOGICAL :: written=.false.
      PUBLIC :: d_wigner

                
      CONTAINS
c***************************************************
c        private routine for integer rotation where only
c        one rotation is passed 
c***************************************************
      SUBROUTINE integer_wigner_1op(
     >                          mrot,bmat,lmax,
     <                          d_wgn)

      INTEGER, INTENT(IN)  :: lmax
      INTEGER, INTENT(IN)  :: mrot(3,3)
      REAL,    INTENT(IN)  :: bmat(3,3)
      COMPLEX, INTENT(OUT) :: d_wgn(-lmax:lmax,-lmax:lmax,lmax)
      REAL                    realmrot(3,3,1)
      
      realmrot(:,:,1) = mrot(:,:)
      CALL real_wigner(
     >                 1,realmrot,bmat,lmax,
     <                 d_wgn)

      END SUBROUTINE

c***************************************************
c        private routine for integer rotation
c***************************************************
      SUBROUTINE integer_wigner(
     >                          nop,mrot,bmat,lmax,
     <                          d_wgn)

      INTEGER, INTENT(IN)  :: nop,lmax
      INTEGER, INTENT(IN)  :: mrot(3,3,nop)
      REAL,    INTENT(IN)  :: bmat(3,3)
      COMPLEX, INTENT(OUT) :: d_wgn(-lmax:lmax,-lmax:lmax,lmax,nop)
      REAL                    realmrot(3,3,nop)
      
      realmrot(:,:,:) = mrot(:,:,:)
      CALL real_wigner(
     >                 nop,realmrot,bmat,lmax,
     <                 d_wgn)

      END SUBROUTINE

c**************************************************
c         private routine for real rotation
c**************************************************
      SUBROUTINE real_wigner(
     >                    nop,mrot,bmat,lmax,
     <                    d_wgn)

      USE m_constants, ONLY : pimach
      USE m_inv3

      IMPLICIT NONE

! .. arguments:
      INTEGER, INTENT(IN)  :: nop,lmax
      REAL,    INTENT(IN)  :: mrot(3,3,nop)
      REAL,    INTENT(IN)  :: bmat(3,3)
      COMPLEX, INTENT(OUT) :: d_wgn(-lmax:lmax,-lmax:lmax,lmax,nop)

! .. local variables:
      INTEGER              :: ns,signum
      INTEGER              :: i,j,k,l,m,mp,x_lo,x_up,x,e_c,e_s
      REAL                 :: fac_l_m,fac_l_mp,fac_lmpx,fac_lmx,
     +                        fac_x,fac_xmpm
      REAL                 :: pi,co_bh,si_bh,zaehler,nenner,cp,sp
      REAL                 :: sina,sinb,sinc,cosa,cosb,cosc,determ,dt
      COMPLEX              :: ci,phase_g,phase_a,bas,
     +                        d(-lmax:lmax,-lmax:lmax)

      REAL                 :: alpha(nop),beta(nop),gamma(nop)
      REAL                 :: dmat(3,3),dmati(3,3),det(nop),bmati(3,3)

      INTRINSIC sqrt,max,min


      ci = cmplx(0.0,1.0)
      pi = pimach()
c
c determine the eulerian angles of all the rotations
c
      CALL inv3(bmat,bmati,dt)
      DO ns = 1, nop
c
c first determine the determinant of the rotation
c +1 for a proper rotation
c -1 for a proper rotation times inversion
c
        determ = 0.00
        DO i = 1,3
          DO j = 1,3
            IF (i.NE.j) THEN
              k = 6 - i - j
              signum = 1
              IF ( (i.EQ.(j+1)).OR.(j.EQ.(k+1)) ) signum=-signum
              determ = determ + signum*
     +                 mrot(i,1,ns)*mrot(j,2,ns)*mrot(k,3,ns)
            ENDIF
          ENDDO
        ENDDO
        IF (abs(1.0-abs(determ)).GT.1.0e-5)
     +       CALL juDFT_error("d_wigner: determ/==+-1",calledby
     +       ="d_wigner")
        det(ns) = determ
c
c store the proper part of the rotation in dmati and convert to
c cartesian coordinates
c
        dmati = determ*matmul(bmati,matmul(mrot(:,:,ns),bmat))
c
c the eulerian angles are derived from the inverse of
c dmati, because we use the convention that we rotate functions
c
        CALL inv3(dmati,dmat,dt)
c
c beta follows directly from d33
c
        cosb = dmat(3,3)
        sinb = 1.00 - cosb*cosb
        sinb = max(sinb,0.00)
        sinb = sqrt(sinb)
c
c if beta = 0 or pi , only alpha+gamma or -gamma have a meaning:
c
        IF ( abs(sinb).LT.1.0e-5 ) THEN 

          beta(ns) = 0.0
          IF ( cosb.lt.0.0 ) beta(ns) = pi
          gamma(ns) = 0.0
          cosa = dmat(1,1)/cosb
          sina = dmat(1,2)/cosb
          IF ( abs(sina).LT.1.0e-5 ) THEN
            alpha(ns)=0.0
            IF ( cosa.LT.0.0 ) alpha(ns)=alpha(ns)+pi
          ELSE
            alpha(ns) = 0.5*pi - atan(cosa/sina)
            IF ( sina.LT.0.0 ) alpha(ns)=alpha(ns)+pi
          ENDIF

        ELSE

          beta(ns) = 0.5*pi - atan(cosb/sinb)
c
c determine alpha and gamma from d13 d23 d32 d31
c
          cosa = dmat(3,1)/sinb
          sina = dmat(3,2)/sinb
          cosc =-dmat(1,3)/sinb
          sinc = dmat(2,3)/sinb
          IF ( abs(sina).lt.1.0e-5 ) THEN
            alpha(ns)=0.0
            IF ( cosa.LT.0.0 ) alpha(ns)=alpha(ns)+pi
          ELSE
            alpha(ns) = 0.5*pi - atan(cosa/sina)
            IF ( sina.LT.0.0 ) alpha(ns)=alpha(ns)+pi
          ENDIF
          IF ( abs(sinc).lt.1.0e-5 ) THEN
            gamma(ns) = 0.0
            IF ( cosc.LT.0.0 ) gamma(ns)=gamma(ns)+pi
          ELSE
            gamma(ns) = 0.5*pi - atan(cosc/sinc)
            IF ( sinc.LT.0.0 ) gamma(ns)=gamma(ns)+pi
          ENDIF

        ENDIF

      ENDDO ! loop over nop

#ifndef CPP_MPI
      IF(.NOT.written) THEN
        WRITE (6,8000)
        DO ns = 1, nop
          WRITE (6,8010) ns,alpha(ns),beta(ns),gamma(ns),det(ns)
        ENDDO
        written=.true.
      ENDIF
 8000 FORMAT(//,'   eulerian angles for the rotations ',
     $  //,'   ns   alpha     beta      gamma    determ ')
 8010 FORMAT(i5,4f10.5)
#endif
      DO ns = 1, nop

        co_bh = cos(beta(ns)*0.5)
        si_bh = sin(beta(ns)*0.5)

        DO l = 1, lmax

          DO m = -l,l
            fac_l_m = fac(l+m) * fac(l-m)
            phase_g = exp( - ci * gamma(ns) * m ) 

            DO mp = -l,l
              fac_l_mp = fac(l+mp) * fac(l-mp)

              zaehler = sqrt( real(fac_l_m * fac_l_mp) )
              phase_a = exp( - ci * alpha(ns) * mp ) 
              x_lo = max(0, m-mp)
              x_up = min(l-mp, l+m)

              bas = zaehler * phase_a * phase_g 
              d(m,mp) = cmplx(0.0,0.0)
              DO x = x_lo,x_up
                fac_lmpx = fac(l-mp-x)
                fac_lmx  = fac(l+m-x)
                fac_x    = fac(x)
                fac_xmpm = fac(x+mp-m)
                nenner = fac_lmpx * fac_lmx * fac_x * fac_xmpm
                e_c = 2*l + m - mp - 2*x 
                e_s = 2*x + mp - m
                IF (e_c.EQ.0) THEN
                  cp = 1.0
                ELSE
                  cp = co_bh ** e_c
                ENDIF
                IF (e_s.EQ.0) THEN
                  sp = 1.0
                ELSE
                  sp = si_bh ** e_s
                ENDIF
                d(m,mp) = d(m,mp) + bas * (-1)**x * cp * sp / nenner
              ENDDO

            ENDDO ! loop over mp
          ENDDO   ! loop over m

          DO m = -l,l
            DO mp = -l,l
             d( m,mp ) = d( m,mp ) * (-1)**(m-mp)
            ENDDO
          ENDDO
          DO m = -l,l
            DO mp = -l,l
              IF(abs(det(ns)+1).lt.1e-5) THEN
                d_wgn(m,mp,l,ns) = d( m,mp) * (-1)**l ! adds inversion
              ELSE
                d_wgn(m,mp,l,ns) = d( m,mp)
              ENDIF
            ENDDO
          ENDDO

        ENDDO
      ENDDO

      END SUBROUTINE real_wigner

      ELEMENTAL REAL FUNCTION  fac(n)

      INTEGER, INTENT (IN) :: n
      INTEGER :: i
 
      fac = 0
      IF (n.LT.0) RETURN
      fac = 1
      IF (n.EQ.0) RETURN
      DO i = 2,n
        fac = fac * i
      ENDDO

      END FUNCTION  fac
      
      END MODULE m_dwigner
