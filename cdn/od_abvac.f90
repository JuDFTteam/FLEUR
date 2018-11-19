!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_od_abvac
CONTAINS
  SUBROUTINE od_abvac(&
       &     cell,vacuum,DIMENSION,stars,&
       &     oneD,qssbti,&
       &     n2d_1,&
       &     wronk,evac,bkpt,MM,vM,&
       &     vz,kvac3,nv2,&
       &     uz,duz,u,udz,dudz,ddnv,ud)
    !**************************************************************
    !      determines the nesessary values and derivatives on the 
    !      vacuum cylindrical boundary for finding a and b coefficients
    !      for the construcing vacuum charge density in vacden.F
    !                          Y.Mokrousov, 7th of october 2002
    !*************************************************************** 
    USE m_vacuz
    USE m_vacudz
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell

    !     .. scalar Arguments..
    REAL wronk
    INTEGER, INTENT (in) :: MM ,vM
    INTEGER, INTENT (in) :: n2d_1,nv2
    REAL,    INTENT (in) :: evac
    !     ..array arguments..

    INTEGER, INTENT (in) :: kvac3(DIMENSION%nv2d)
    REAL,    INTENT (in) :: bkpt(3),qssbti 
    REAL,    INTENT (in) :: vz(vacuum%nmz) 
    REAL,    INTENT (out):: udz(DIMENSION%nv2d,-vM:vM)
    REAL,    INTENT (out):: uz(DIMENSION%nv2d,-vM:vM)
    REAL,    INTENT (out):: dudz(DIMENSION%nv2d,-vM:vM)
    REAL,    INTENT (out):: duz(DIMENSION%nv2d,-vM:vM)
    REAL,    INTENT (out):: u(vacuum%nmz,DIMENSION%nv2d,-vM:vM)
    REAL,    INTENT (out):: ud(vacuum%nmz,DIMENSION%nv2d,-vM:vM)
    REAL,    INTENT (out):: ddnv(DIMENSION%nv2d,-vM:vM)
    !     ..local scalars..
    REAL ev,scale,xv,yv,vzero,v1
    INTEGER i,ik,jk,jspin,jsp1,jsp2 ,l,m
    INTEGER i1,i2,i3,ind1,ind3
    !     .. local arrays..
    REAL wdz(DIMENSION%nv2d,-vM:vM),wz(DIMENSION%nv2d,-vM:vM)
    REAL dwdz(DIMENSION%nv2d,-vM:vM),dwz(DIMENSION%nv2d,-vM:vM)
    REAL v(3),x(vacuum%nmz)
    REAL  vr0(vacuum%nmz)
    REAL w(vacuum%nmz,DIMENSION%nv2d,-vM:vM),wd(vacuum%nmz,DIMENSION%nv2d,-vM:vM)

    !     wronksian for the schrodinger equation given by an identity

    wronk = 2.0

    DO  ik = 1,nv2
       DO  m = 0,vM
          v(1) = 0.0
          v(2) = 0.0
          v(3) = bkpt(3) + kvac3(ik) + qssbti
          ev = evac - 0.5*DOT_PRODUCT(v,MATMUL(v,cell%bbmat))

          !     constructing of the 'pseudopotential'

          DO  i=1,vacuum%nmz
             v1 = 1./(8.*((cell%z1+(i-1)*vacuum%delz)**2))&
                  &              -(m*m)/(2.*((cell%z1+(i-1)*vacuum%delz)**2))
             vr0(i) = vz(i)-v1
          ENDDO
          vzero = vr0(vacuum%nmz)

          !     obtaining solutions with the 'pseudopotential'

          CALL vacuz(ev,vr0(1),vzero,vacuum%nmz,vacuum%delz,&
                          wz(ik,m),dwz(ik,m),w(1,ik,m))
          CALL vacudz(ev,vr0(1),vzero,vacuum%nmz,vacuum%delz,&
               wdz(ik,m),dwdz(ik,m),ddnv(ik,m), wd(1,ik,m),dwz(ik,m),w(1,ik,m))

          scale = wronk/(wdz(ik,m)*dwz(ik,m)-&
               &           dwdz(ik,m)*wz(ik,m))
          wdz(ik,m) = scale*wdz(ik,m)
          dwdz(ik,m) = scale*dwdz(ik,m)
          ddnv(ik,m) = scale*ddnv(ik,m)
          IF (m.GT.0) THEN
             wdz(ik,-m) = wdz(ik,m)
             dwdz(ik,-m) = dwdz(ik,m)
             ddnv(ik,-m) = ddnv(ik,m)
          END IF
          DO  i = 1,vacuum%nmz
             wd(i,ik,m) = scale*wd(i,ik,m)
             w(i,ik,m) = scale*w(i,ik,m)
             IF (m.GT.0) THEN
                wd(i,ik,-m) = wd(i,ik,m)
                w(i,ik,-m) = w(i,ik,m)
             END IF
          ENDDO
          !     constructing 'real' solutions

          DO  i=1,vacuum%nmz
             u(i,ik,m)=w(i,ik,m)/SQRT(cell%z1+(i-1)*vacuum%delz)
             ud(i,ik,m)=wd(i,ik,m)/SQRT(cell%z1+(i-1)*vacuum%delz)
             IF (m.GT.0) THEN
                u(i,ik,-m) = u(i,ik,m)
                ud(i,ik,-m) = ud(i,ik,m)
             END IF
          ENDDO
          duz(ik,m)=(-dwz(ik,m))/SQRT(cell%z1)-&
               &           wz(ik,m)/(2.0*((cell%z1)**(1.5)))
          uz(ik,m)=wz(ik,m)/SQRT(cell%z1)
          dudz(ik,m)=(-dwdz(ik,m))/SQRT(cell%z1)-&
               &           wdz(ik,m)/(2.0*((cell%z1)**(1.5)))
          udz(ik,m)=wdz(ik,m)/SQRT(cell%z1)
          IF (m.GT.0) THEN
             duz(ik,-m) = duz(ik,m)
             uz(ik,-m) = uz(ik,m)
             dudz(ik,-m) = dudz(ik,m)
             udz(ik,-m) = udz(ik,m)
          END IF

       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE od_abvac
END MODULE m_od_abvac
