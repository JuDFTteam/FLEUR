!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_nmat_rot

! Calculate the Wigner rotation matrices for complex spherical
! harmonics for all space-group rotations and l=1,2,3. Needed 
! for the calculation of the density matrix in nmat.
!
! also allows to use rotated "n_mmp_mat" file by specifying a 
! "n_mmp_rot" file (see its use in u_mix or u_setup) 
!                                                       gb10
      CONTAINS
      SUBROUTINE nmat_rot(
     >                    alpha,beta,gamma,l_in,n_u,jspins,lty,
     X                    n_mmp)

      use m_constants
      USE m_inv3

      IMPLICIT NONE

! .. arguments:
      INTEGER, INTENT(IN)  :: l_in,n_u,jspins,lty(n_u)
      REAL,    INTENT(IN)  :: alpha(n_u),beta(n_u),gamma(n_u)
      COMPLEX, INTENT(INOUT) :: n_mmp(-3:3,-3:3,n_u,jspins)

! .. local variables:
      INTEGER ns,signum,ispin,n
      INTEGER i,j,k,l,m,mp,x_lo,x_up,x,e_c,e_s
      REAL fac_l_m,fac_l_mp,fac_lmpx,fac_lmx,fac_x,fac_xmpm
      REAL co_bh,si_bh,zaehler,nenner,cp,sp
      REAL sina,sinb,sinc,cosa,cosb,cosc,determ,dt
      COMPLEX phase_g,phase_a,bas,d(-l_in:l_in,-l_in:l_in)
      COMPLEX d_wig(-l_in:l_in,-l_in:l_in,l_in,n_u)
      COMPLEX n_tmp(-l_in:l_in,-l_in:l_in)
      COMPLEX nr_tmp(-l_in:l_in,-l_in:l_in)
      LOGICAL, SAVE :: written = .false.

      REAL dmat(3,3),dmati(3,3)

      

      
      DO n = 1, n_u

      co_bh = cos(beta(n)*0.5)
      si_bh = sin(beta(n)*0.5)

      DO l = 1, lty(n)
        d = (0.0,0.0)

        DO m = -l,l
          fac_l_m = fac(l+m) * fac(l-m)
          phase_g = exp( - ImagUnit * gamma(n) * m )

          DO mp = -l,l
            fac_l_mp = fac(l+mp) * fac(l-mp)

            zaehler = sqrt( real(fac_l_m * fac_l_mp) )
            phase_a = exp( - ImagUnit * alpha(n) * mp ) 
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
        d_wig(:,:,l,n) = d(:,:)

      ENDDO ! l
      ENDDO ! n 

      DO ispin = 1, jspins
        DO n = 1, n_u
           n_tmp(:,:) = n_mmp(-l_in:l_in,-l_in:l_in,n,ispin)
           d(:,:) = d_wig(:,:,lty(n),n)

           nr_tmp = matmul( transpose( conjg(d) ) , n_tmp)
           n_tmp =  matmul( nr_tmp, d )

           n_mmp(-l_in:l_in,-l_in:l_in,n,ispin) = n_tmp(:,:)
         ENDDO
      ENDDO
      !write(*,'(14f8.4)') n_mmp

      END SUBROUTINE nmat_rot

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
      
      END MODULE m_nmat_rot
