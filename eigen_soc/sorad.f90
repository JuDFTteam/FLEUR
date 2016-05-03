MODULE m_sorad
  USE m_juDFT
  !*********************************************************************
  !     1. generates radial spin-orbit matrix elements
  !     based on m.weinert's radsra and radsrd subroutines
  !*********************************************************************
CONTAINS
  SUBROUTINE sorad(atoms,input,ntyp,vr,enpara,spav,&
       rsopp,rsopdpd,rsoppd,rsopdp,usdus,&
       rsoplop,rsoplopd,rsopdplo,rsopplo,rsoploplop)

    USE m_constants, ONLY : c_light
    USE m_intgr,     ONLY : intgr0
    USE m_sointg
    USE m_radsra
    USE m_radsrd
    USE m_radsrdn
    USE m_types
    IMPLICIT NONE
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_usdus),INTENT(INOUT)   :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ntyp
    LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,:)!(atoms%jmtd,dimension%jspd),
    REAL,    INTENT (OUT) :: rsopp  (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsoppd (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsopdp (atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsopdpd(atoms%ntypd,atoms%lmaxd,2,2)
    REAL,    INTENT (OUT) :: rsoplop (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsoplopd(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsopdplo(atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsopplo (atoms%ntypd,atoms%nlod,2,2)
    REAL,    INTENT (OUT) :: rsoploplop(atoms%ntypd,atoms%nlod,atoms%nlod,2,2)
    !     ..
    !     .. Local Scalars ..
    REAL ddn1,e ,ulops,dulops,duds1
    INTEGER i,j,ir,jspin,l,noded,nodeu,ilo,ilop
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: p(:,:),pd(:,:),q(:,:),qd(:,:),plo(:,:)
    REAL, ALLOCATABLE :: plop(:,:),glo(:,:),fint(:),pqlo(:,:)
    REAL, ALLOCATABLE :: filo(:,:)
    REAL, ALLOCATABLE :: v0(:),vso(:,:),qlo(:,:)
    !     ..

    IF (atoms%jri(ntyp)>atoms%jmtd)  CALL juDFT_error("atoms%jri(ntyp).GT.atoms%jmtd",calledby ="sorad")
    ALLOCATE ( p(atoms%jmtd,2),pd(atoms%jmtd,2),q(atoms%jmtd,2),plo(atoms%jmtd,2),fint(atoms%jmtd),&
         &   qlo(atoms%jmtd,2),plop(atoms%jmtd,2),qd(atoms%jmtd,2),v0(atoms%jmtd),vso(atoms%jmtd,2) )
    !
    DO l = 0,atoms%lmax(ntyp) 

       DO jspin = 1,input%jspins
          !
          !--->    calculate normalized function at e: p and q 
          !
          e = enpara%el0(l,ntyp,jspin)
          CALL radsra(&
               e,l,vr(:,jspin),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
               usdus%us(l,ntyp,jspin),usdus%dus(l,ntyp,jspin),&
               nodeu,p(:,jspin),q(:,jspin))
          !                     
          !--->    calculate orthogonal energy derivative at e : pd and qd
          !
          CALL radsrd(&
               e,l,vr(:,jspin),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
               usdus%uds(l,ntyp,jspin),usdus%duds(l,ntyp,jspin),&
               usdus%ddn(l,ntyp,jspin),noded,pd(:,jspin),qd(:,jspin),&
               p(:,jspin),q(:,jspin),usdus%dus(l,ntyp,jspin))

       END DO     ! end of spin loop
       !
       !---> in case of jspins=1
       !

       IF (input%jspins.EQ.1) THEN
          DO i = 1,atoms%jri(ntyp)
             p(i,2) =  p(i,1)
             pd(i,2) = pd(i,1)
          ENDDO
       ENDIF
       !
       !---> common spin-orbit integrant V   (average spin directions)
       !                                  SO
       v0(:) = 0.0
       IF (input%jspins.EQ.1) THEN
          v0(1:atoms%jri(ntyp)) = vr(1:atoms%jri(ntyp),1)
          e = enpara%el0(l,ntyp,1)
       ELSE
          DO i = 1,atoms%jri(ntyp)
             v0(i) = (vr(i,1)+vr(i,input%jspins))/2.
          END DO
          e = (enpara%el0(l,ntyp,1)+enpara%el0(l,ntyp,input%jspins))/2.
       END IF

       CALL sointg(ntyp,e,vr,v0,atoms,input,vso)
       IF (spav) THEN
          DO i= 1,atoms%jmtd
             vso(i,1)= (vso(i,1)+vso(i,2))/2.
             vso(i,2)= vso(i,1)
          ENDDO
       ENDIF

       !                        s       s'            .s       s'
       !-->  radial integrals <u  |V  |u  > = rsopp, <u  |V  |u  > = rsopdp etc.
       !                            SO                     SO

       IF (l.GT.0) THEN ! there is no spin-orbit for s-states
          DO i = 1, 2
             DO j = 1, 2
                rsopp(ntyp,l,i,j) = radso( p(:atoms%jri(ntyp),i), p(:atoms%jri(ntyp),j),vso(:atoms%jri(ntyp),i),atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                rsopdp(ntyp,l,i,j) = radso(pd(:atoms%jri(ntyp),i), p(:atoms%jri(ntyp),j),vso(:atoms%jri(ntyp),i),atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                rsoppd(ntyp,l,i,j) = radso( p(:atoms%jri(ntyp),i),pd(:atoms%jri(ntyp),j),vso(:atoms%jri(ntyp),i),atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                rsopdpd(ntyp,l,i,j) = radso(pd(:atoms%jri(ntyp),i),pd(:atoms%jri(ntyp),j),vso(:atoms%jri(ntyp),i),atoms%dx(ntyp),atoms%rmsh(1,ntyp))
             ENDDO
          ENDDO
       ENDIF ! l>0
       !
       !--->  Check for local orbitals with same l
       !
       DO ilo = 1, atoms%nlo(ntyp)
          IF (atoms%llo(ilo,ntyp).EQ.l) THEN

             DO jspin = 1,input%jspins
                e = enpara%ello0(ilo,ntyp,jspin)
                CALL radsra(&
                     e,l,vr(:,jspin),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                     usdus%ulos(ilo,ntyp,jspin),usdus%dulos(ilo,ntyp,jspin),&
                     nodeu,plo(:,jspin),qlo(:,jspin))

                !+apw+lo
                IF (atoms%l_dulo(ilo,ntyp).OR.atoms%ulo_der(ilo,ntyp).GE.1) THEN !  calculate energy derivative (of order atoms%ulo_der) at e
                   ALLOCATE (glo(atoms%jmtd,2),pqlo(atoms%jmtd,2),filo(atoms%jmtd,2))
                   pqlo(1:atoms%jri(ntyp),1)=plo(1:atoms%jri(ntyp),jspin)
                   pqlo(1:atoms%jri(ntyp),2)=qlo(1:atoms%jri(ntyp),jspin)
                   i = atoms%ulo_der(ilo,ntyp)
                   IF(atoms%l_dulo(ilo,ntyp)) i=1
                   CALL radsrdn(&
                        e,l,vr(:,jspin),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                        usdus%ulos(ilo,ntyp,jspin),duds1,ddn1,noded,glo,filo,&!filo is a dummy array&
                        pqlo,usdus%dulos(ilo,ntyp,jspin),i)
                   ddn1 = SQRT(ddn1)
                   IF(atoms%l_dulo(ilo,ntyp)) ddn1=1.0
                   plo(1:atoms%jri(ntyp),jspin) = glo(1:atoms%jri(ntyp),1)/ddn1
                   qlo(1:atoms%jri(ntyp),jspin) = glo(1:atoms%jri(ntyp),2)/ddn1
                   usdus%dulos(ilo,ntyp,jspin) = duds1/ddn1
                   DEALLOCATE (glo,pqlo,filo)
                ENDIF
                !-apw+lo
             ENDDO

             IF (input%jspins.EQ.1) THEN
                plo(1:atoms%jri(ntyp),2) = plo(1:atoms%jri(ntyp),1)
                e = (enpara%ello0(ilo,ntyp,1) + enpara%el0(l,ntyp,1) )/2
             ELSE
                e = (enpara%ello0(ilo,ntyp,1) +  enpara%ello0(ilo,ntyp,input%jspins) +&
                     enpara%el0(l,ntyp,1) + enpara%el0(l,ntyp,input%jspins) )/4
             END IF
             CALL sointg(ntyp,e,vr,v0,atoms,input, vso)
             IF (spav) THEN
                DO i= 1,atoms%jmtd
                   vso(i,1)= (vso(i,1)+vso(i,2))/2.
                   vso(i,2)= vso(i,1)
                ENDDO
             ENDIF

             DO i = 1, 2
                DO j = 1, 2
                   rsoplop (ntyp,ilo,i,j) = radso(plo(:,i),p (:,j),vso(:,i), atoms%rmsh(1,ntyp),atoms%dx(ntyp))
                   rsoplopd(ntyp,ilo,i,j) = radso(plo(:,i),pd(:,j),vso(:,i), atoms%rmsh(1,ntyp),atoms%dx(ntyp))
                   rsopplo (ntyp,ilo,i,j) = radso(p (:,i),plo(:,j),vso(:,i), atoms%rmsh(1,ntyp),atoms%dx(ntyp))
                   rsopdplo(ntyp,ilo,i,j) = radso(pd(:,i),plo(:,j),vso(:,i), atoms%rmsh(1,ntyp),atoms%dx(ntyp))
                ENDDO
             ENDDO

             DO i = 1,input%jspins
                fint(:) = plo(:,i) *  p(:,i) + qlo(:,i) *  q(:,i)
                CALL intgr0(fint,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%uulon(ilo,ntyp,i))
                fint(:) = plo(:,i) * pd(:,i) + qlo(:,i) * qd(:,i)
                CALL intgr0(fint,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%dulon(ilo,ntyp,i))
             ENDDO

             DO ilop = 1, atoms%nlo(ntyp)
                IF (atoms%llo(ilop,ntyp).EQ.l) THEN

                   DO jspin = 1,input%jspins
                      e = enpara%ello0(ilop,ntyp,jspin)
                      CALL radsra(&
                           e,l,vr(:,jspin),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                           ulops,dulops,nodeu,plop(:,jspin),q(:,1))
                      !+apw+lo
                      IF (atoms%l_dulo(ilo,ntyp).OR.atoms%ulo_der(ilo,ntyp).GE.1) THEN ! calculate orthogonal energy derivative at e
                         ALLOCATE (glo(atoms%jmtd,2),pqlo(atoms%jmtd,2),filo(atoms%jmtd,2))
                         pqlo(1:atoms%jri(ntyp),1)=plop(1:atoms%jri(ntyp),jspin)
                         pqlo(1:atoms%jri(ntyp),2)=q(1:atoms%jri(ntyp),1)
                         i = atoms%ulo_der(ilo,ntyp)
                         IF(atoms%l_dulo(ilo,ntyp)) i=1
                         CALL radsrdn(&
                              e,l,vr(:,jspin),atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                              ulops,duds1,ddn1,noded,glo,filo,&!filo is a dummy array&
                              pqlo,dulops,i)
                         plop(1:atoms%jri(ntyp),jspin) = glo(1:atoms%jri(ntyp),1)
                         DEALLOCATE (glo,pqlo,filo)
                      ENDIF
                      !-apw+lo
                   ENDDO

                   IF (input%jspins.EQ.1) THEN
                      plop(1:atoms%jri(ntyp),2) = plop(1:atoms%jri(ntyp),1)
                      e = (enpara%ello0(ilo,ntyp,1) + enpara%ello0(ilop,ntyp,1) )/2
                   ELSE
                      e = (enpara%ello0(ilo,ntyp,1) +  enpara%ello0(ilo,ntyp,input%jspins) +  &
                           enpara%ello0(ilop,ntyp,1) + enpara%ello0(ilop,ntyp,input%jspins) )/4
                   END IF
                   CALL sointg(ntyp,e,vr,v0,atoms,input, vso)
                   IF (spav) THEN
                      DO i= 1,atoms%jmtd
                         vso(i,1)= (vso(i,1)+vso(i,2))/2.
                         vso(i,2)= vso(i,1)
                      ENDDO
                   ENDIF

                   DO i = 1, 2
                      DO j = 1, 2
                         rsoploplop(ntyp,ilo,ilop,i,j) =&
                              radso(plo(:,i),plop(:,j),vso(:,i),atoms%rmsh(1,ntyp),atoms%dx(ntyp))
                      ENDDO
                   ENDDO

                ENDIF
             ENDDO

          ENDIF
       ENDDO ! end of lo-loop

    ENDDO ! end of l-loop

    DEALLOCATE ( p,pd,q,qd,plo,plop,qlo,fint,v0,vso )
    !      rsoplop (:,:,:,:) = 0.0
    !      rsoplopd(:,:,:,:) = 0.0
    !      rsopplo (:,:,:,:) = 0.0
    !      rsopdplo(:,:,:,:) = 0.0
    !      rsoploplop(:,:,:,:,:) = 0.0

  END SUBROUTINE sorad
  !--------------------------------------------------------------------
  REAL FUNCTION radso(a,b,vso,dx,r0)
    !
    !     compute radial spin-orbit integrals
    !
    USE m_intgr, ONLY : intgr0
    IMPLICIT NONE
    !
    !     .. Scalar Arguments ..
    REAL,    INTENT (IN) :: r0,dx
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: a(:),b(:),vso(:)
    !     ..
    !     .. Local Arrays ..
    REAL q(size(a))
    !     ..
    q = a*b*vso
    CALL intgr0(q,r0,dx,size(a),radso)

    RETURN
  END FUNCTION radso
  !--------------------------------------------------------------------
END MODULE m_sorad
