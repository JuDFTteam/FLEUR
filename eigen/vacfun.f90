MODULE m_vacfun
  use m_juDFT
CONTAINS
  SUBROUTINE vacfun(&
       vacuum,stars,input,noco,jspin1,jspin2,&
       sym, cell,ivac,evac,bkpt, vxy,vz,kvac,nv2,&
       tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk)
    !*********************************************************************
    !     determines the necessary values and derivatives on the vacuum
    !     boundary (ivac=1 upper vacuum; ivac=2, lower) for energy
    !     parameter evac.  also sets up the 2d hamiltonian matrices
    !     necessary to update the full hamiltonian matrix.
    !               m. weinert
    !*********************************************************************

    USE m_intgr, ONLY : intgz0
    USE m_vacuz
    USE m_vacudz
    USE m_types
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ivac,jspin1,jspin2
    REAL,    INTENT (OUT) :: wronk
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nv2(:)!(input%jspins)
    INTEGER, INTENT (IN) :: kvac(:,:,:)!(2,dimension%nv2d,input%jspins)
    COMPLEX, INTENT (IN) :: vxy(:,:,:,:) !(vacuum%nmzxyd,stars%ng2-1,nvac,:)
    COMPLEX, INTENT (OUT):: tddv(:,:),tduv(:,:)!(dimension%nv2d,dimension%nv2d)
    COMPLEX, INTENT (OUT):: tudv(:,:),tuuv(:,:)!(dimension%nv2d,dimension%nv2d)
    REAL,ALLOCATABLE,INTENT (IN) :: vz(:,:,:) !(vacuum%nmzd,2,4) ,
    REAL,    INTENT (IN) :: evac(:,:)!(2,input%jspins)
    REAL,    INTENT (IN) :: bkpt(3) 
    REAL,    INTENT (OUT):: udz(:,:),uz(:,:)!(dimension%nv2d,input%jspins)
    REAL,    INTENT (OUT):: dudz(:,:),duz(:,:)!(dimension%nv2d,input%jspins)
    REAL,    INTENT (OUT):: ddnv(:,:)!(dimension%nv2d,input%jspins)
    !     ..
    !     .. Local Scalars ..
    REAL ev,scale,xv,yv,vzero,fac
    COMPLEX phase
    INTEGER i,i1,i2,i3,ik,ind2,ind3,jk,np1,jspin,ipot
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    REAL u(vacuum%nmzd,size(duz,1),input%jspins),ud(vacuum%nmzd,size(duz,1),input%jspins)
    REAL v(3),x(vacuum%nmzd), qssbti(2,2)
    !     ..
    fac=MERGE(1.0,-1.0,jspin1>=jspin2)
    ipot=MERGE(jspin1,3,jspin1==jspin2)

    tuuv=0.0;tudv=0.0;tddv=0.0;tduv=0.0
    udz=0.0;duz=0.0;ddnv=0.0;udz=0.;uz=0.
    tail = .true.
    np1 = vacuum%nmzxy + 1
    !--->    wronksian for the schrodinger equation given by an identity
    wronk = 2.0
    !---> setup the spin-spiral q-vector
    qssbti(1:2,1) = - noco%qss(1:2)/2
    qssbti(1:2,2) = + noco%qss(1:2)/2
    !--->    generate basis functions for each 2-d k+g
    DO jspin = MIN(jspin1,jspin2),MAX(jspin1,jspin2)
       DO  ik = 1,nv2(jspin)
          v(1:2) = bkpt(1:2) + kvac(:,ik,jspin) + qssbti(1:2,jspin)
          v(3) = 0.0
          ev = evac(ivac,jspin) - 0.5*dot_product(v,matmul(v,cell%bbmat))
          vzero = vz(vacuum%nmzd,ivac,jspin)
          CALL vacuz(ev,vz(1,ivac,jspin),vzero,vacuum%nmz,vacuum%delz,&
               uz(ik,jspin),duz(ik,jspin),u(1,ik,jspin))
          CALL vacudz(ev,vz(1,ivac,jspin),vzero,vacuum%nmz,vacuum%delz,&
               udz(ik,jspin),dudz(ik,jspin),ddnv(ik,jspin),&
               ud(1,ik,jspin),duz(ik,jspin),u(1,ik,jspin))
          !--->       make sure the solutions satisfy the wronksian
          scale = wronk/ (udz(ik,jspin)*duz(ik,jspin)-&
               &                         dudz(ik,jspin)*uz(ik,jspin))
          udz(ik,jspin)  = scale*udz(ik,jspin)
          dudz(ik,jspin) = scale*dudz(ik,jspin)
          ddnv(ik,jspin) = scale*ddnv(ik,jspin)
          ud(:,ik,jspin) = scale*ud(:,ik,jspin)
       enddo
    ENDDO
    !--->    set up the tuuv, etc. matrices
    DO  ik = 1,nv2(jspin1)
       DO  jk = 1,nv2(jspin2)

          !--->     determine the warping component of the potential
          i1 = kvac(1,ik,jspin1) - kvac(1,jk,jspin2)
          i2 = kvac(2,ik,jspin1) - kvac(2,jk,jspin2)
          i3 = 0
          ind3 = stars%ig(i1,i2,i3)
          IF (ind3.EQ.0) CYCLE
          phase = stars%rgphs(i1,i2,i3)
          ind2 = stars%ig2(ind3)
          IF (ind2.EQ.0) THEN
             WRITE (6,FMT=8000) ik,jk
8000         FORMAT (' **** error in map2 for 2-d stars',2i5)
             CALL juDFT_error("error in map2 for 2-d stars",calledby ="vacfun")
          END IF
          !--->     get the proper warping index (vxy starts with the 2nd star)
          ind2 = ind2 - 1
          IF (ind2.NE.0) THEN
             !--->       only the warping part, 1st star (G=0) is done later

             !--->       obtain the warping matrix elements
             !--->       note that the tail correction (tail=.true.) is included for
             !--->       the integrals, i.e. the integrand is from infinity inward

             !--->       tuuv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*REAL(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tuuv(ik,jk) = phase*cmplx(xv,yv)

             !--->       tddv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*REAL(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) =ud(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tddv(ik,jk) = phase*cmplx(xv,yv)

             !--->       tudv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*real(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tudv(ik,jk) = phase*cmplx(xv,yv)

             !--->       tduv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*real(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tduv(ik,jk) = phase*cmplx(xv,yv)

          ELSE
             !--->       diagonal (film muffin-tin) terms
             IF (jspin1==jspin2) THEN
                tuuv(ik,ik) = cmplx(evac(ivac,jspin1),0.0)
                tddv(ik,ik) = cmplx(evac(ivac,jspin1)*ddnv(ik,jspin1),0.0)
                tudv(ik,ik) = cmplx(0.5,0.0)
                tduv(ik,ik) = cmplx(0.5,0.0)
             ELSE
                !--->          tuuv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*fac*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tuuv(ik,jk) = cmplx(xv,yv)

                !--->          tddv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*fac*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tddv(ik,jk) = cmplx(xv,yv)

                !--->          tudv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*fac*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tudv(ik,jk) = cmplx(xv,yv)

                !--->          tduv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*fac*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tduv(ik,jk) = cmplx(xv,yv)
             ENDIF

          ENDIF
       enddo
    enddo

  END SUBROUTINE vacfun
END MODULE m_vacfun

