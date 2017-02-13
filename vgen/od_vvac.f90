!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_od_vvac
CONTAINS
  SUBROUTINE od_vvac(&
       &     stars,vacuum,cell,&
       &     psq,rht,&
       &     vz)

    !     subroutine which calculates the non warped part of the
    !     vacuum potential (m=0,gz=0)
    !                               Y. Mokrousov
    !     the potential in this subroutine can be defined in two
    !     equivalent ways, which nevertheless give a bit defferent 
    !     results, 2nd one seems to be more precise 

    USE m_qsf
    USE m_od_cylbes
    USE m_types
    USE m_constants
    IMPLICIT NONE
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_cell),INTENT(IN)     :: cell


    COMPLEX, INTENT (IN) :: psq(stars%ng3)
    REAL,    INTENT (IN) :: rht(:,:,:) !(vacuum%nmzd,2,dimension%jspd)
    REAL,    INTENT (OUT) :: vz(:,:,:) !(vacuum%nmzd,2,dimension%jspd)

    COMPLEX  rhobar
    INTEGER  k1,k2,irec3,irec2,i,j,ivac,imz,imz1
    REAL     g2 ,a(vacuum%nmzd)
    REAL     fJ,z,zp,phi
    REAL     rht1(vacuum%nmzd)
    REAL     f2(vacuum%nmzd),f22(vacuum%nmzd)

    INTRINSIC cmplx


    DO i = 1,vacuum%nmz
       f2(i) = 0.
       f22(i) = 0.
       DO ivac = 1,vacuum%nvac
          vz(i,ivac,1) = 0.
       END DO
    END DO


    rhobar = -psq(1)

    DO  k1 = -stars%mx1,stars%mx1
       DO  k2 = -stars%mx2,stars%mx2
          irec3 = stars%ig(k1,k2,0)
          IF (irec3.NE.0) THEN
             irec2 = stars%ig2(irec3)
             IF (irec2.NE.1) THEN
                g2 = stars%sk2(irec2)
                phi = stars%phi2(irec2)
                CALL od_cylbes(1,cell%z1*g2,fJ)
                rhobar = rhobar - 2.*psq(irec3)*CMPLX(fJ/(g2*cell%z1),0.0)

             END IF
          END IF
       ENDDO
    ENDDO
    !----> 1st equivalent way      

    DO  i=1,vacuum%nmz
       rht1(i) = fpi_const*(cell%z1+(i-1)*vacuum%delz)*rht(i,1,1)
    ENDDO
    CALL qsf(vacuum%delz,rht1(1),f2(1),vacuum%nmz,1)

    DO  i = 1,vacuum%nmz
       f2(i) = tpi_const*cell%z1*cell%z1*rhobar-f2(i)
    ENDDO

    DO  i = 1,vacuum%nmz
       DO  j = 1,vacuum%nmz
          IF (j.LT.i) THEN
             f22(j) = 0.0
          ELSE
             f22(j) = f2(j)/(cell%z1+vacuum%delz*(j-1))
          END IF
       ENDDO
       CALL qsf(vacuum%delz,f22(1),a,vacuum%nmz,0)
       DO  ivac =1,vacuum%nvac
          vz(i,ivac,1) = -a(1)
       ENDDO
    ENDDO
    !----> 2nd equivalent way (via the Green function)

    DO imz = 1,vacuum%nmz
       z = cell%z1 + (imz-1)*vacuum%delz
       DO imz1 = 1,vacuum%nmz
          zp = cell%z1 +  (imz1-1)*vacuum%delz
          IF (imz1.LE.imz) THEN
             rht1(imz1) = fpi_const*LOG(z)*zp*rht(imz1,1,1)
          ELSE
             rht1(imz1) = fpi_const*LOG(zp)*zp*rht(imz1,1,1)
          END IF

       END DO
       CALL qsf(vacuum%delz,rht1,a,vacuum%nmz,0)
       vz(imz,1,1) = tpi_const*LOG(z)*(cell%z1*cell%z1)*rhobar - a(1)
    END DO

    RETURN
  END SUBROUTINE od_vvac
END MODULE m_od_vvac


      
   
   






      



      

