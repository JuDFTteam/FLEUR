MODULE m_vacfun
  use m_juDFT
#ifdef CPP_MPI
  use mpi
#endif
CONTAINS
  SUBROUTINE vacfun(&
       fmpi,vacuum,stars,input,nococonv,jspin1,jspin2,&
       cell,ivac,evac,bkpt, vxy,vz,kvac,nv2,&
       tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk,&
       bkptq,v1xy,v1z,kvacq,nv2q,uzq,duzq,udzq,dudzq,wronkq)
    !*********************************************************************
    !     determines the necessary values and derivatives on the vacuum
    !     boundary (ivac=1 upper vacuum; ivac=2, lower) for energy
    !     parameter evac.  also sets up the 2d hamiltonian matrices
    !     necessary to update the full hamiltonian matrix.
    !               m. weinert
    !*********************************************************************

    USE m_constants
    USE m_intgr, ONLY : intgz0
    USE m_vacuz
    USE m_vacudz
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_stars),INTENT(IN)       :: stars
    TYPE(t_cell),INTENT(IN)        :: cell
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ivac,jspin1,jspin2
    REAL,    INTENT (OUT) :: wronk
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nv2(:)!(input%jspins)
    INTEGER, INTENT (IN) :: kvac(:,:,:)!(2,lapw%dim_nv2d(),input%jspins)
    COMPLEX, INTENT (IN) :: vxy(:,:,:,:) !(vacuum%nmzxyd,stars%ng2-1,nvac,:)
    COMPLEX, INTENT (OUT):: tddv(:,:),tduv(:,:)!(lapw%dim_nv2d(),lapw%dim_nv2d())
    COMPLEX, INTENT (OUT):: tudv(:,:),tuuv(:,:)!(lapw%dim_nv2d(),lapw%dim_nv2d())
    COMPLEX, INTENT (IN) :: vz(:,:,:) !(vacuum%nmzd,2,4) ,
    REAL,    INTENT (IN) :: evac(:,:)!(2,input%jspins)
    REAL,    INTENT (IN) :: bkpt(3)
    REAL,    INTENT (OUT):: udz(:,:),uz(:,:)!(lapw%dim_nv2d(),input%jspins)
    REAL,    INTENT (OUT):: dudz(:,:),duz(:,:)!(lapw%dim_nv2d(),input%jspins)
    REAL,    INTENT (OUT):: ddnv(:,:)!(lapw%dim_nv2d(),input%jspins)
   ! Optional DFPT stuff
    REAL,    OPTIONAL, INTENT (IN) :: bkptq(3)
    REAL,    OPTIONAL, INTENT (OUT):: udzq(:,:),uzq(:,:)!(lapw%dim_nv2d(),input%jspins)
    REAL,    OPTIONAL, INTENT (OUT):: dudzq(:,:),duzq(:,:)!(lapw%dim_nv2d(),input%jspins)
    REAL,    OPTIONAL, INTENT (OUT) :: wronkq
    INTEGER, OPTIONAL, INTENT (IN) :: kvacq(:,:,:)!(2,lapw%dim_nv2d(),input%jspins)
    INTEGER, OPTIONAL, INTENT (IN) :: nv2q(:)!(input%jspins)
    COMPLEX, OPTIONAL, INTENT (IN) :: v1xy(:,:,:,:) !(vacuum%nmzxyd,stars%ng2-1,nvac,:)
    COMPLEX, OPTIONAL, INTENT (IN) :: v1z(:,:,:) !(vacuum%nmzd,2,4) ,
    !     ..
    !     .. Local Scalars ..
    REAL ev,scale,xv,yv,vzero,fac
    COMPLEX phase
    INTEGER i,i1,i2,i3,ik,ind2,ind3,jk,np1,jspin,ipot,nbuf,ierr,loclen
    INTEGER mat_start,mat_end
    LOGICAL tail, l_dfpt
    !     ..
    !     .. Local Arrays ..
    REAL u(vacuum%nmzd,size(duz,1),input%jspins),ud(vacuum%nmzd,size(duz,1),input%jspins)
    REAL v(3),x(vacuum%nmzd), qssbti(2,2)
    COMPLEX, ALLOCATABLE :: tddv_loc(:,:), tduv_loc(:,:), tudv_loc(:,:), tuuv_loc(:,:)
    COMPLEX, ALLOCATABLE :: tv_gather_buf(:)
    REAL :: ddnvq(SIZE(ddnv,1),input%jspins)
    ! DFPT
    REAL uq(vacuum%nmzd,size(duz,1),input%jspins),udq(vacuum%nmzd,size(duz,1),input%jspins)
    !     ..
    l_dfpt = .FALSE.
    IF (PRESENT(bkptq)) l_dfpt = .TRUE.
    fac=MERGE(1.0,-1.0,jspin1>=jspin2)
    ipot=MERGE(jspin1,3,jspin1==jspin2)

    tuuv=0.0;tudv=0.0;tddv=0.0;tduv=0.0
    udz=0.0;duz=0.0;ddnv=0.0;udz=0.;uz=0.
    tail = .true.
    IF (l_dfpt) THEN
      udzq=0.0;duzq=0.0;ddnvq=0.0;udzq=0.;uzq=0.
    END IF
    np1 = vacuum%nmzxy + 1
    !--->    wronksian for the schrodinger equation given by an identity
    wronk = 2.0
    !---> setup the spin-spiral q-vector
    qssbti(1:2,1) = - nococonv%qss(1:2)/2
    qssbti(1:2,2) = + nococonv%qss(1:2)/2
    !--->    generate basis functions for each 2-d k+g
    DO jspin = MIN(jspin1,jspin2),MAX(jspin1,jspin2)
       DO  ik = 1,nv2(jspin)
          v(1:2) = bkpt(1:2) + kvac(:,ik,jspin) + qssbti(1:2,jspin)
          v(3) = 0.0
          ev = evac(ivac,jspin) - 0.5*dot_product(v,matmul(v,cell%bbmat))
          vzero = vz(vacuum%nmzd,ivac,jspin)
          CALL vacuz(ev,REAL(vz(1:,ivac,jspin)),vzero,vacuum%nmz,vacuum%delz,&
               uz(ik,jspin),duz(ik,jspin),u(1,ik,jspin))
          CALL vacudz(ev,REAL(vz(1:,ivac,jspin)),vzero,vacuum%nmz,vacuum%delz,&
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
       if (l_dfpt) then
         DO  ik = 1,nv2q(jspin)
            v(1:2) = bkptq(1:2) + kvacq(:,ik,jspin) + qssbti(1:2,jspin)
            v(3) = 0.0
            ev = evac(ivac,jspin) - 0.5*dot_product(v,matmul(v,cell%bbmat))
            vzero = vz(vacuum%nmzd,ivac,jspin)
            CALL vacuz(ev,REAL(vz(1:,ivac,jspin)),vzero,vacuum%nmz,vacuum%delz,&
                  uzq(ik,jspin),duzq(ik,jspin),uq(1,ik,jspin))
            CALL vacudz(ev,REAL(vz(1:,ivac,jspin)),vzero,vacuum%nmz,vacuum%delz,&
                  udzq(ik,jspin),dudzq(ik,jspin),ddnvq(ik,jspin),&
                  udq(1,ik,jspin),duzq(ik,jspin),uq(1,ik,jspin))
            !--->       make sure the solutions satisfy the wronksian
            scale = wronk/ (udzq(ik,jspin)*duzq(ik,jspin)-&
                  &                         dudzq(ik,jspin)*uzq(ik,jspin)) !! TODO: Output wronsk always = 2.0?
            udzq(ik,jspin)  = scale*udzq(ik,jspin)
            dudzq(ik,jspin) = scale*dudzq(ik,jspin)
            ddnvq(ik,jspin) = scale*ddnvq(ik,jspin)
            udq(:,ik,jspin) = scale*udq(:,ik,jspin)
         enddo
       end if
    ENDDO
    loclen = size(tddv,2)/fmpi%n_size + 1
    mat_start = fmpi%n_rank*loclen + 1
    mat_end = (fmpi%n_rank+1)*loclen
    ALLOCATE(tddv_loc(size(tddv,1),mat_start:mat_end))
    ALLOCATE(tduv_loc(size(tddv,1),mat_start:mat_end))
    ALLOCATE(tudv_loc(size(tddv,1),mat_start:mat_end))
    ALLOCATE(tuuv_loc(size(tddv,1),mat_start:mat_end))
    tuuv_loc=0.0;tudv_loc=0.0;tddv_loc=0.0;tduv_loc=0.0

    !--->    set up the tuuv, etc. matrices
    DO  jk = mat_start,mat_end
       IF (jk>nv2(jspin2)) EXIT
       IF (.NOT.l_dfpt) THEN
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP& SHARED(tuuv_loc,tddv_loc,tudv_loc,tduv_loc,ddnv,vz,jk) &
       !$OMP& SHARED(stars,jspin1,jspin2,evac,nv2,kvac,vacuum,u,vxy,tail,fac,np1,ivac,ipot,ud) &
       !$OMP& PRIVATE(i1,i2,i3,ind3,phase,ind2,x,xv,yv)
         DO  ik = 1,nv2(jspin1)

            !--->     determine the warping component of the potential
            i1 = fac*(kvac(1,ik,jspin1) - kvac(1,jk,jspin2))
            i2 = fac*(kvac(2,ik,jspin1) - kvac(2,jk,jspin2))
            i3 = 0
            ind3 = stars%ig(i1,i2,i3)
            IF (ind3.EQ.0) CYCLE
            phase = stars%rgphs(i1,i2,i3)
            ind2 = stars%ig2(ind3)
            IF (ind2.EQ.0) THEN
               WRITE (oUnit,FMT=8000) ik,jk
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
               tuuv_loc(ik,jk) = phase*cmplx(xv,yv)

               !--->       tddv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*REAL(vxy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) =ud(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tddv_loc(ik,jk) = phase*cmplx(xv,yv)

               !--->       tudv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*real(vxy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tudv_loc(ik,jk) = phase*cmplx(xv,yv)

               !--->       tduv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*real(vxy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(vxy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tduv_loc(ik,jk) = phase*cmplx(xv,yv)

            ELSE
               !--->       diagonal (film muffin-tin) terms
               IF (jspin1==jspin2) THEN
                  tuuv_loc(ik,ik) = cmplx(evac(ivac,jspin1),0.0)
                  tddv_loc(ik,ik) = cmplx(evac(ivac,jspin1)*ddnv(ik,jspin1),0.0)
                  tudv_loc(ik,ik) = cmplx(0.5,0.0)
                  tduv_loc(ik,ik) = cmplx(0.5,0.0)
               ELSE

                  !--->          tuuv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*vz(i,ivac,3)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(vz(i,ivac,3))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tuuv_loc(ik,jk) = cmplx(xv,yv)

                  !--->          tddv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*vz(i,ivac,3)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(vz(i,ivac,3))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tddv_loc(ik,jk) = cmplx(xv,yv)

                  !--->          tudv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*vz(i,ivac,3)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(vz(i,ivac,3))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tudv_loc(ik,jk) = cmplx(xv,yv)

                  !--->          tduv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*vz(i,ivac,3)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(vz(i,ivac,3))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tduv_loc(ik,jk) = cmplx(xv,yv)
               ENDIF

            ENDIF
         ENDDO
         !$OMP END PARALLEL DO
       ELSE
         !$OMP PARALLEL DO DEFAULT(none) &
         !$OMP& SHARED(tuuv_loc,tddv_loc,tudv_loc,tduv_loc,ddnv,vz,v1z,jk,bkpt,bkptq) &
         !$OMP& SHARED(stars,jspin1,jspin2,evac,nv2,kvac,kvacq,vacuum,u,uq,v1xy,tail,fac,np1,ivac,ipot,ud,udq,l_dfpt) &
         !$OMP& PRIVATE(i1,i2,i3,ind3,phase,ind2,x,xv,yv)
          DO  ik = 1,nv2q(jspin1)
            !--->     determine the warping component of the potential
            i1 = fac*(kvac(1,ik,jspin1) - kvacq(1,jk,jspin2))
            i2 = fac*(kvac(2,ik,jspin1) - kvacq(2,jk,jspin2))
            i3 = 0
            ind3 = stars%ig(i1,i2,i3)
            IF (ind3.EQ.0) CYCLE
            phase = stars%rgphs(i1,i2,i3)
            ind2 = stars%ig2(ind3)
            IF (ind2.EQ.0) THEN
               WRITE (oUnit,FMT=8001) ik,jk
   8001         FORMAT (' **** error in map2 for 2-d stars',2i5)
               CALL juDFT_error("error in map2 for 2-d stars",calledby ="vacfun")
            END IF
            !--->     get the proper warping index (v1xy starts with the 2nd star)
            ind2 = ind2 - 1
            IF ((ind2.NE.0).OR.(norm2(bkptq-bkpt)>1e-8)) THEN
               !--->       tuuv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = uq(i,ik,jspin1)*u(i,jk,jspin2)*REAL(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = uq(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tuuv_loc(ik,jk) = phase*cmplx(xv,yv)

               !--->       tddv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = udq(i,ik,jspin1)*ud(i,jk,jspin2)*REAL(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) =udq(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tddv_loc(ik,jk) = phase*cmplx(xv,yv)

               !--->       tudv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = uq(i,ik,jspin1)*ud(i,jk,jspin2)*real(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = uq(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tudv_loc(ik,jk) = phase*cmplx(xv,yv)

               !--->       tduv
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = udq(i,ik,jspin1)*u(i,jk,jspin2)*real(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,xv,tail)
               DO  i = 1,vacuum%nmzxy
                  x(np1-i) = udq(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(v1xy(i,ind2,ivac,ipot))
               enddo
               CALL intgz0(x,vacuum%delz,vacuum%nmzxy,yv,tail)
               tduv_loc(ik,jk) = phase*cmplx(xv,yv)

            ELSE
               !--->       diagonal (film muffin-tin) terms
                  !--->          tuuv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*v1z(i,ivac,ipot)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(v1z(i,ivac,ipot))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tuuv_loc(ik,jk) = cmplx(xv,yv)

                  !--->          tddv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*v1z(i,ivac,ipot)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(v1z(i,ivac,ipot))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tddv_loc(ik,jk) = cmplx(xv,yv)

                  !--->          tudv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*v1z(i,ivac,ipot)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = u(i,ik,jspin1)*ud(i,jk,jspin2)*fac*AIMAG(v1z(i,ivac,ipot))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tudv_loc(ik,jk) = cmplx(xv,yv)

                  !--->          tduv
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*v1z(i,ivac,ipot)
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,xv,tail)
                  DO i = 1,vacuum%nmz
                     x(vacuum%nmz+1-i) = ud(i,ik,jspin1)*u(i,jk,jspin2)*fac*AIMAG(v1z(i,ivac,ipot))
                  ENDDO
                  CALL intgz0(x,vacuum%delz,vacuum%nmz,yv,tail)
                  tduv_loc(ik,jk) = cmplx(xv,yv)
            ENDIF
         ENDDO
         !$OMP END PARALLEL DO
       END IF
    ENDDO

    IF ( fmpi%n_size == 1 ) THEN
       tuuv = tuuv_loc(:,:size(tuuv,2))
       tduv = tduv_loc(:,:size(tduv,2))
       tudv = tudv_loc(:,:size(tudv,2))
       tddv = tddv_loc(:,:size(tddv,2))
    ELSE
       nbuf = size(tuuv_loc)
#ifdef CPP_MPI
       ALLOCATE(tv_gather_buf(nbuf*fmpi%n_size))
       CALL MPI_ALLGATHER(tuuv_loc,nbuf,MPI_DOUBLE_COMPLEX,tv_gather_buf,nbuf,MPI_DOUBLE_COMPLEX,fmpi%sub_comm,ierr)
       CALL zcopy(size(tuuv),tv_gather_buf,1,tuuv,1)
       CALL MPI_ALLGATHER(tduv_loc,nbuf,MPI_DOUBLE_COMPLEX,tv_gather_buf,nbuf,MPI_DOUBLE_COMPLEX,fmpi%sub_comm,ierr)
       CALL zcopy(size(tduv),tv_gather_buf,1,tduv,1)
       CALL MPI_ALLGATHER(tudv_loc,nbuf,MPI_DOUBLE_COMPLEX,tv_gather_buf,nbuf,MPI_DOUBLE_COMPLEX,fmpi%sub_comm,ierr)
       CALL zcopy(size(tudv),tv_gather_buf,1,tudv,1)
       CALL MPI_ALLGATHER(tddv_loc,nbuf,MPI_DOUBLE_COMPLEX,tv_gather_buf,nbuf,MPI_DOUBLE_COMPLEX,fmpi%sub_comm,ierr)
       CALL zcopy(size(tddv),tv_gather_buf,1,tddv,1)
       DEALLOCATE (tv_gather_buf)
#endif
    ENDIF

    DEALLOCATE(tddv_loc, tduv_loc, tudv_loc, tuuv_loc)

  END SUBROUTINE vacfun
END MODULE m_vacfun
