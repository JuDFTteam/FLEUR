MODULE m_vacfun
  use m_juDFT
CONTAINS
  SUBROUTINE vacfun(&
       vacuum,DIMENSION,stars, jsp,input,noco,jsp1,jsp2,&
       sym, cell,ivac,evac,bkpt, vxy,vz,kvac1,kvac2,nv2,&
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

    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp ,ivac,jsp1,jsp2
    REAL,    INTENT (OUT) :: wronk
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nv2(input%jspins)
    INTEGER, INTENT (IN) :: kvac1(dimension%nv2d,input%jspins),kvac2(dimension%nv2d,input%jspins)
    COMPLEX, INTENT (IN) :: vxy(vacuum%nmzxyd,stars%ng2-1)
    COMPLEX, INTENT (OUT):: tddv(dimension%nv2d,dimension%nv2d),tduv(dimension%nv2d,dimension%nv2d)
    COMPLEX, INTENT (OUT):: tudv(dimension%nv2d,dimension%nv2d),tuuv(dimension%nv2d,dimension%nv2d)
    REAL,    INTENT (IN) :: vz(vacuum%nmzd,2,4) ,evac(2,input%jspins)
    REAL,    INTENT (IN) :: bkpt(3) 
    REAL,    INTENT (OUT):: udz(dimension%nv2d,input%jspins),uz(dimension%nv2d,input%jspins)
    REAL,    INTENT (OUT):: dudz(dimension%nv2d,input%jspins),duz(dimension%nv2d,input%jspins)
    REAL,    INTENT (OUT):: ddnv(dimension%nv2d,input%jspins)
    !     ..
    !     .. Local Scalars ..
    REAL ev,scale,xv,yv,vzero
    COMPLEX phase
    INTEGER i,i1,i2,i3,ik,ind2,ind3,jk,np1,jspin
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    REAL u(vacuum%nmzd,dimension%nv2d,input%jspins),ud(vacuum%nmzd,dimension%nv2d,input%jspins)
    REAL v(3),x(vacuum%nmzd), qssbti(2,2)
    !     ..
    tuuv=0.0;tudv=0.0;tddv=0.0;tduv=0.0
    udz=0.0;duz=0.0;ddnv=0.0;udz=0.;uz=0.
    tail = .true.
    np1 = vacuum%nmzxy + 1
    !--->    wronksian for the schrodinger equation given by an identity
    wronk = 2.0
    !---> setup the spin-spiral q-vector
    qssbti(1,1) = - noco%qss(1)/2
    qssbti(2,1) = - noco%qss(2)/2
    qssbti(1,2) = + noco%qss(1)/2
    qssbti(2,2) = + noco%qss(2)/2
    !--->    generate basis functions for each 2-d k+g
    DO jspin = 1,input%jspins
       DO  ik = 1,nv2(jspin)
          v(1) = bkpt(1) + kvac1(ik,jspin) + qssbti(1,jspin)
          v(2) = bkpt(2) + kvac2(ik,jspin) + qssbti(2,jspin)
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
          udz(ik,jspin) = scale*udz(ik,jspin)
          dudz(ik,jspin) = scale*dudz(ik,jspin)
          ddnv(ik,jspin) = scale*ddnv(ik,jspin)
          ud(:,ik,jspin) = scale*ud(:,ik,jspin)
       enddo
    ENDDO
    !--->    set up the tuuv, etc. matrices
    DO  ik = 1,nv2(jsp1)
       DO  jk = 1,nv2(jsp2)

          !--->     determine the warping component of the potential
          i1 = kvac1(ik,jsp1) - kvac1(jk,jsp2)
          i2 = kvac2(ik,jsp1) - kvac2(jk,jsp2)
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
                x(np1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*real(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*aimag(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tuuv(ik,jk) = phase*cmplx(xv,yv)

             !--->       tddv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = ud(i,ik,jsp1)*ud(i,jk,jsp2)*real(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) =ud(i,ik,jsp1)*ud(i,jk,jsp2)*aimag(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tddv(ik,jk) = phase*cmplx(xv,yv)

             !--->       tudv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*real(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*aimag(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tudv(ik,jk) = phase*cmplx(xv,yv)

             !--->       tduv
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*real(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
             DO  i = 1,vacuum%nmzxy
                x(np1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*aimag(vxy(i,ind2))
             enddo
             CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
             tduv(ik,jk) = phase*cmplx(xv,yv)

          ELSE

             !--->       diagonal (film muffin-tin) terms
             IF (jsp1==jsp2) THEN
                tuuv(ik,ik) = cmplx(evac(ivac,jsp1),0.0)
                tddv(ik,ik) = cmplx(evac(ivac,jsp1)*ddnv(ik,jsp1),0.0)
                tudv(ik,ik) = cmplx(0.5,0.0)
                tduv(ik,ik) = cmplx(0.5,0.0)
             ELSE

                !--->          tuuv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tuuv(ik,jk) = cmplx(xv,yv)

                !--->          tddv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tddv(ik,jk) = cmplx(xv,yv)

                !--->          tudv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = u(i,ik,jsp1)*ud(i,jk,jsp2)*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tudv(ik,jk) = cmplx(xv,yv)

                !--->          tduv
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,3)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                DO i = 1,vacuum%nmz
                   x(vacuum%nmz+1-i) = ud(i,ik,jsp1)*u(i,jk,jsp2)*vz(i,ivac,4)
                ENDDO
                CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                tduv(ik,jk) = cmplx(xv,yv)
             ENDIF

          ENDIF
       enddo
    enddo

  END SUBROUTINE vacfun
END MODULE m_vacfun

