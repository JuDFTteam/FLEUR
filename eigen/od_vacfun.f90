!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_od_vacfun
CONTAINS
  SUBROUTINE od_vacfun(&
       m_cyl,cell,vacuum,DIMENSION,stars,&
       jsp,input,noco,ipot,oneD, n2d_1, ivac,evac,bkpt,MM,vM,&
       vxy,vz,kvac3,nv2, tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv)
    !*********************************************************************
    !     determines the necessary values and derivatives on the cylindrical
    !     vacuum boundary for energy parameter evac. also sets up the
    !     hamiltonian matrices necessary to update the full hamiltonian
    !     matrix.
    !     Y. Mokrousov, June 2002
    !*********************************************************************
    USE m_vacuz
    USE m_vacudz
    USE m_intgr, ONLY : intgz0
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN):: DIMENSION
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_cell),INTENT(IN)     :: cell
    !     ..
    !     .. Scalar Arguments ..

    INTEGER, INTENT (IN) :: jsp   ,ipot,MM ,vM
    INTEGER, INTENT (IN) :: ivac,n2d_1,m_cyl
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nv2(input%jspins)
    INTEGER, INTENT (IN) :: kvac3(DIMENSION%nv2d,input%jspins)
    COMPLEX, INTENT (IN) :: vxy(vacuum%nmzxyd,n2d_1-1)
    COMPLEX, INTENT (OUT):: tddv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX, INTENT (OUT):: tduv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX, INTENT (OUT):: tudv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX, INTENT (OUT):: tuuv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d)
    REAL,    INTENT (IN) :: vz(vacuum%nmzd,2,4) ,evac(2,input%jspins)
    REAL,    INTENT (IN) :: bkpt(3) 
    REAL,    INTENT (OUT):: udz(-vM:vM,DIMENSION%nv2d,input%jspins),uz(-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL,    INTENT (OUT):: dudz(-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL,    INTENT (OUT):: duz(-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL,    INTENT (OUT):: ddnv(-vM:vM,DIMENSION%nv2d,input%jspins)
    !     ..
    !     .. Local Scalars ..
    REAL ev,scale,xv,yv,vzero,v1,wronk
    INTEGER i,ik,jk,np1,jspin,jsp1,jsp2 ,l,m
    INTEGER i1,i2,i3,ind1,ind3
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    REAL wdz(-vM:vM,DIMENSION%nv2d,input%jspins),wz(-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL dwdz(-vM:vM,DIMENSION%nv2d,input%jspins),dwz(-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL u(vacuum%nmzd,-vM:vM,DIMENSION%nv2d,input%jspins),ud(vacuum%nmzd,-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL v(3),x(vacuum%nmzd)
    REAL vr0(vacuum%nmzd,2,4)
    REAL w(vacuum%nmzd,-vM:vM,DIMENSION%nv2d,input%jspins),wd(vacuum%nmzd,-vM:vM,DIMENSION%nv2d,input%jspins)
    REAL qssbti(2)
    !     ..


    tail = .TRUE.
    np1 = vacuum%nmzxy + 1

    !     wronksian for the schrodinger equation given by an identity

    wronk = 2.0

    qssbti(1) = - noco%qss(3)/2.
    qssbti(2) = + noco%qss(3)/2.

    tuuv(:,:,:,:) = CMPLX(0.,0.)
    tudv(:,:,:,:) = CMPLX(0.,0.)
    tduv(:,:,:,:) = CMPLX(0.,0.)
    tddv(:,:,:,:) = CMPLX(0.,0.)

    !     generate basis functions for each 1-d k_z+g_z and m

    DO jspin = 1,input%jspins
       DO  ik = 1,nv2(jspin)
          DO  m = 0,vM
             v(1) = 0.0
             v(2) = 0.0
             v(3) = bkpt(3) + kvac3(ik,jspin) + qssbti(jspin)
             ev = evac(ivac,jspin) - 0.5*DOT_PRODUCT(v,MATMUL(v,cell%bbmat))
             !     constructing of the 'pseudopotential'
             DO  i=1,vacuum%nmzd
                v1 = 1./(8.*((cell%z1+(i-1)*vacuum%delz)**2))&
                     -(m*m)/(2.*((cell%z1+(i-1)*vacuum%delz)**2))
                vr0(i,ivac,jspin) = vz(i,ivac,jspin)-v1
             ENDDO
             vzero = vr0(vacuum%nmzd,ivac,jspin)
             !     obtaining solutions with the 'pseudopotential'

             CALL vacuz(ev,vr0(:,ivac,jspin),vzero,vacuum%nmz,vacuum%delz,&
                  wz(m,ik,jspin),dwz(m,ik,jspin),w(1,m,ik,jspin))
             CALL vacudz(ev,vr0(1,ivac,jspin),vzero,vacuum%nmz,vacuum%delz,&
                  wdz(m,ik,jspin),dwdz(m,ik,jspin),ddnv(m,ik,jspin),&
                  wd(1,m,ik,jspin),dwz(m,ik,jspin),w(1,m,ik,jspin))
             !     make sure the solutions satisfy the wronksian
             scale = wronk/ (wdz(m,ik,jspin)*dwz(m,ik,jspin)-&
                  dwdz(m,ik,jspin)*wz(m,ik,jspin))
             wdz(m,ik,jspin) = scale*wdz(m,ik,jspin)
             dwdz(m,ik,jspin) = scale*dwdz(m,ik,jspin)
             ddnv(m,ik,jspin) = scale*ddnv(m,ik,jspin)
             IF (m.GT.0) THEN
                wdz(-m,ik,jspin) = wdz(m,ik,jspin)
                dwdz(-m,ik,jspin) = dwdz(m,ik,jspin)
                ddnv(-m,ik,jspin) = ddnv(m,ik,jspin)
             END IF
             DO  i = 1,vacuum%nmz
                wd(i,m,ik,jspin) = scale*wd(i,m,ik,jspin)
                w(i,m,ik,jspin) = scale*w(i,m,ik,jspin)
                IF (m.GT.0) THEN
                   wd(i,-m,ik,jspin) = wd(i,m,ik,jspin)
                   w(i,-m,ik,jspin) = w(i,m,ik,jspin)
                END IF
             ENDDO
             !     constructing 'real' solutions
             DO  i=1,vacuum%nmz
                u(i,m,ik,jspin)=w(i,m,ik,jspin)/SQRT(cell%z1+(i-1)*vacuum%delz)
                ud(i,m,ik,jspin)=wd(i,m,ik,jspin)/SQRT(cell%z1+(i-1)*vacuum%delz)
                IF (m.GT.0) THEN
                   u(i,-m,ik,jspin) = u(i,m,ik,jspin)
                   ud(i,-m,ik,jspin) = ud(i,m,ik,jspin)
                END IF
             ENDDO
             duz(m,ik,jspin)=(-dwz(m,ik,jspin))/SQRT(cell%z1)-&
                  wz(m,ik,jspin)/(2.0*((cell%z1)**(1.5)))
             uz(m,ik,jspin)=wz(m,ik,jspin)/SQRT(cell%z1)
             dudz(m,ik,jspin)=(-dwdz(m,ik,jspin))/SQRT(cell%z1)-&
                  wdz(m,ik,jspin)/(2.0*((cell%z1)**(1.5))) 
             udz(m,ik,jspin)=wdz(m,ik,jspin)/SQRT(cell%z1)
             IF (m.GT.0) THEN
                duz(-m,ik,jspin) = duz(m,ik,jspin)
                uz(-m,ik,jspin) = uz(m,ik,jspin)
                dudz(-m,ik,jspin) = dudz(m,ik,jspin)
                udz(-m,ik,jspin) = udz(m,ik,jspin)
             END IF
          ENDDO
       ENDDO
    ENDDO

    !     set up the tuuv, etc. matrices

    IF (noco%l_noco) THEN
       IF (ipot.EQ.1) THEN
          jsp1 = 1
          jsp2 = 1
       ELSEIF (ipot.EQ.2) THEN
          jsp1 = 2
          jsp2 = 2
       ELSEIF (ipot.EQ.3) THEN
          jsp1 = 2
          jsp2 = 1
       ENDIF
    ELSE
       jsp1 = jsp
       jsp2 = jsp
    ENDIF

    DO  ik = 1,nv2(jsp1)
       DO  jk = 1,nv2(jsp2)
          i1 = 0
          i2 = 0
          i3 = kvac3(ik,jsp1) - kvac3(jk,jsp2)
          ind3 = stars%ig(i1,i2,i3)
          IF (ind3.EQ.0) CYCLE
          DO  m = -vM,vM
             DO  l = -vM,vM
                IF (l.EQ.m .OR. (iabs(m).LE.m_cyl&
                     &                             .AND. iabs(l).LE.m_cyl)) THEN
                   !     determine the warping component of the potential
                   ind1 = oneD%ig1(i3,m-l)
                   IF (ind1.NE.0) THEN
                      IF(ind1.NE.1) THEN
                         ind1 = ind1 - 1
                         !     only the warping part
                         !--->             tuuv
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2) *REAL(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2) *AIMAG(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
                         tuuv(m,l,ik,jk) = CMPLX(xv,yv)
                         !--->             tddv
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *REAL(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *AIMAG(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
                         tddv(m,l,ik,jk) = CMPLX(xv,yv)
                         !--->             tudv
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *REAL(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *AIMAG(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
                         tudv(m,l,ik,jk) = CMPLX(xv,yv)
                         !--->             tduv
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2) *REAL(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
                         DO i = 1,vacuum%nmzxy
                            x(np1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2) *AIMAG(vxy(i,ind1))
                         ENDDO
                         CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
                         tduv(m,l,ik,jk) = CMPLX(xv,yv)
                      ELSE
                         !--->          diagonal terms
                         IF ((ipot.EQ.1).OR.(ipot.EQ.2)) THEN
                            tuuv(m,m,ik,ik) = CMPLX(evac(ivac,jsp1),0.0)
                            tddv(m,m,ik,ik) = CMPLX(evac(ivac,jsp1)*&
                                 &                       ddnv(m,ik,jsp1),0.0)
                            tudv(m,m,ik,ik) = CMPLX(0.5,0.0)
                            tduv(m,m,ik,ik) = CMPLX(0.5,0.0)
                         ELSE
                            !--->             tuuv
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2) *vz(i,ivac,3)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = w(i,m,ik,jsp1)*w(i,l,jk,jsp2) *vz(i,ivac,4)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                            tuuv(m,l,ik,jk) = CMPLX(xv,yv)
                            !--->             tddv
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *vz(i,ivac,3)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = wd(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *vz(i,ivac,4)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                            tddv(m,l,ik,jk) = CMPLX(xv,yv)
                            !--->             tudv
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *vz(i,ivac,3)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = w(i,m,ik,jsp1)*wd(i,l,jk,jsp2) *vz(i,ivac,4)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                            tudv(m,l,ik,jk) = CMPLX(xv,yv)
                            !--->             tduv
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2) *vz(i,ivac,3)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                            DO i = 1,vacuum%nmz
                               x(vacuum%nmz+1-i) = wd(i,m,ik,jsp1)*w(i,l,jk,jsp2) *vz(i,ivac,4)
                            ENDDO
                            CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                            tduv(m,l,ik,jk) = CMPLX(xv,yv)

                         ENDIF !ipot
                      END IF !ind1 ne 1
                   ENDIF    ! ind1 ne 0
                END IF
             ENDDO
          ENDDO
       ENDDO
    ENDDO


    RETURN
  END SUBROUTINE od_vacfun
END MODULE m_od_vacfun




      
