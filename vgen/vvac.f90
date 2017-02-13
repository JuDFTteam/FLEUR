      MODULE m_vvac
!     ****************************************************************
!     calculates the g(2-dim)=0 part of the vacuum coulomb potential *
!     for general symmetry.          c.l.fu, r.podloucky             *
!     ****************************************************************
      CONTAINS
      SUBROUTINE vvac(&
     &                vacuum,stars,&
     &                cell,sym,input,&
     &                psq,rht,&
     &                vz,rhobar,sig1dh,vz1dh)

      USE m_constants
      USE m_qsf
      USE m_types
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_vacuum),INTENT(IN)  :: vacuum
      TYPE(t_stars),INTENT(IN)   :: stars
      TYPE(t_cell),INTENT(IN)    :: cell
      TYPE(t_sym),INTENT(IN)     :: sym

      COMPLEX, INTENT (OUT):: rhobar
      REAL,    INTENT (OUT):: sig1dh,vz1dh
      TYPE(t_input), INTENT(INOUT) :: input !efield is modified here
!     ..
!     .. Array Arguments ..
      REAL,    INTENT (IN) :: rht(vacuum%nmzd,2,input%jspins) 
      COMPLEX, INTENT (IN) :: psq(stars%ng3) ! pseudo charge density (see psqpw.F)
      REAL,    INTENT (OUT):: vz(vacuum%nmzd,2,input%jspins)
!     ..
!     .. Local Scalar Parameters ..
      COMPLEX, PARAMETER :: ci = (0.0,1.0)
!     ..
!     .. Local Scalars ..
      COMPLEX sumq,vcons
      REAL bj0,bj1,dh,qzh,sigmaa(2)
      INTEGER ig3,imz,ivac,ncsh
!     ..
!     .. Local Arrays ..
      REAL f(vacuum%nmzd),sig(vacuum%nmzd),vtemp(vacuum%nmzd)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC cos,sin

     
      vz(:,1:vacuum%nvac,1) = 0.0  ! initialize potential

      ! obtain mesh point (ncsh) of charge sheet for external electric field
      !                                                 (X.Nie, IFF, 10/97)
      ncsh = input%efield%zsigma/vacuum%delz + 1.01
      sigmaa(1) = (input%efield%sigma+input%efield%sig_b(1))/cell%area
      sigmaa(2) = (input%efield%sigma+input%efield%sig_b(2))/cell%area

      ! g=0 vacuum potential due to neutral charge density
      ! inside slab and zero charge density outside

      vcons = -2.*fpi_const*ci
      dh = cell%z1
      rhobar = -psq(1)
      sumq = cmplx(0.0,0.0)

      DO ig3 = 2,stars%ng3
        IF (stars%ig2(ig3) == 1) THEN           ! select G_|| = 0
          qzh = stars%kv3(3,ig3)*cell%bmat(3,3)*dh
          bj0 = sin(qzh)/qzh
          rhobar = rhobar - psq(ig3)*bj0*stars%nstr(ig3)
          IF (.NOT.(sym%zrfs .OR. sym%invs)) THEN
            bj1 = (sin(qzh) - qzh*cos(qzh)) / (qzh*qzh)
            sumq = sumq + bj1*psq(ig3)*dh*dh
          ENDIF
        ENDIF
      ENDDO

      ivac = 2                        ! lower (ivac=2) vacuum
      IF (vacuum%nvac.EQ.2) THEN
        vz(1:vacuum%nmz,ivac,1) = vcons*sumq
      ENDIF

      ! g=0 vacuum potential due to
      ! negative of rhobar + vacuum (g=0) charge ----> v2(z)

      ivac = 1     ! upper vacuum

      ! =================== DIRICHLET ==============================
      IF (input%efield%dirichlet) THEN
        vz(ncsh+1:vacuum%nmz,ivac,1) = input%efield%sig_b(1)

        CALL qsf(vacuum%delz,rht(1,ivac,1),sig,ncsh,1)
        sig(1:ncsh) = sig(ncsh)- sig(1:ncsh)
        CALL qsf(vacuum%delz,sig,vtemp,ncsh,1)

        DO imz = 1,ncsh
          vz(imz,ivac,1) = -fpi_const* (vtemp(ncsh)-vtemp(imz))&
     &                     + input%efield%sig_b(1)
        ENDDO

        sig1dh = sig(1)
        vz1dh = vz(1,ivac,1)   ! potential on vacuum boundary

        IF (vacuum%nvac == 1) RETURN


        ivac = 2     ! lower vacuum

        CALL qsf(vacuum%delz,rht(1,ivac,1),sig,ncsh,1)

        f(1:ncsh) = sig(1:ncsh) - rhobar*vacuum%dvac + sig1dh

        CALL qsf(vacuum%delz,f,vtemp,ncsh,1)

        DO imz = 1,ncsh
          vz(imz,ivac,1) = -fpi_const* (vtemp(imz)+sig1dh*vacuum%dvac&
     &                     -rhobar*vacuum%dvac*vacuum%dvac/2.) + vz(imz,ivac,1)&
     &                     +vz1dh
        END DO

        ! Force matching on the other side
        input%efield%vslope = (input%efield%sig_b(2)-vz(ncsh,1,1)) / (2*vacuum%delz*(ncsh+1)+vacuum%dvac)
        ivac = 1
        DO imz = 1,ncsh
          vz(imz,ivac,1) = vz(imz,ivac,1) + input%efield%vslope*vacuum%delz*(ncsh-imz+1)
        END DO
        ivac = 2
        DO imz = 1,ncsh
          vz(imz,ivac,1) = vz(imz,ivac,1)&
     &                    + input%efield%vslope*(vacuum%dvac+vacuum%delz*imz+vacuum%delz*ncsh)
        END DO

        vz(ncsh+1:vacuum%nmz,ivac,1) = input%efield%sig_b(2)

      ! =================== NEUMANN ==============================
      ELSE ! Neumann

      CALL qsf(vacuum%delz,rht(1,ivac,1),sig,vacuum%nmz,1)

      sig1dh = sig(vacuum%nmz) - sigmaa(1)  ! need to include contribution from
                                     ! electric field
      sig(1:vacuum%nmz) = sig(vacuum%nmz) - sig(1:vacuum%nmz)

      CALL qsf(vacuum%delz,sig,vtemp,vacuum%nmz,1)

      ! external electric field contribution (X.Nie, IFF, 10/97)
      !                                       corrected 10/99 mw
      DO imz = 1,ncsh
         vz(imz,ivac,1) = -fpi_const* (vtemp(vacuum%nmz)-vtemp(imz)) + vz(imz,ivac,1)&
     &                    -fpi_const*(imz-ncsh)*vacuum%delz*sigmaa(1)
      ENDDO
      DO imz =ncsh+1,vacuum%nmz
         vz(imz,ivac,1) = -fpi_const* (vtemp(vacuum%nmz)-vtemp(imz)) + vz(imz,ivac,1)
      ENDDO

      vz1dh = vz(1,ivac,1)   ! potential on vacuum boundary

      IF (vacuum%nvac.EQ.1) RETURN

      ivac = 2     ! lower vacuum

      CALL qsf(vacuum%delz,rht(1,ivac,1),sig,vacuum%nmz,1)

      f(1:vacuum%nmz) = sig(1:vacuum%nmz) - rhobar*vacuum%dvac + sig1dh

      CALL qsf(vacuum%delz,f,vtemp,vacuum%nmz,1)

      !   external electric field contribution
      DO imz = 1,ncsh
         vz(imz,ivac,1) = -fpi_const* (vtemp(imz)+sig1dh*vacuum%dvac-&
     &                    rhobar*vacuum%dvac*vacuum%dvac/2.) + vz1dh + vz(imz,ivac,1)
      ENDDO
      DO imz =ncsh+1,vacuum%nmz
         vz(imz,ivac,1) = -fpi_const* (vtemp(imz)+sig1dh*vacuum%dvac-&
     &                    rhobar*vacuum%dvac*vacuum%dvac/2.) + vz1dh + vz(imz,ivac,1)&
     &                    +fpi_const*(imz-ncsh)*vacuum%delz*sigmaa(2)
      ENDDO

      END IF ! Dirichlet (vs. Neumann)

      END SUBROUTINE vvac
      END MODULE m_vvac
