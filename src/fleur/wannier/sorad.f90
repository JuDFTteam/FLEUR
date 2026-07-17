!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_sorad
  USE m_juDFT
  !*********************************************************************
  !     generates radial spin-orbit matrix elements
  !*********************************************************************
   implicit none
CONTAINS
  SUBROUTINE sorad(atoms,input,ntyp,vr,enpara,spav,rsoc,usdus,hub1data)

    USE m_constants, ONLY : c_light
    USE m_intgr,     ONLY : intgr0
    USE m_sointg
    USE m_radsra
    USE m_radsrd
    USE m_radsrdn
    USE m_types_enpara
    USE m_types_input
    USE m_types_atoms
    USE m_types_rsoc
    USE m_types_usdus
    USE m_types_hub1data
    IMPLICIT NONE
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_usdus),INTENT(INOUT) :: usdus
    TYPE(t_rsoc),INTENT(INOUT)  :: rsoc
    TYPE(t_hub1data),OPTIONAL,INTENT(IN)  :: hub1data
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ntyp
    LOGICAL, INTENT (IN) :: spav ! if T, spin-averaged pot is used
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(:,:)!(atoms%jmtd,input%jspins),
    !     ..
    !     .. Local Scalars ..
    REAL ddn1,e ,ulops,dulops,duds1
    INTEGER ir,l,noded,nodeu,ilo,ilop, iSpin, jSpin, iRad, iHia, iDeriv
    LOGICAL l_hia
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: p(:,:),pd(:,:),q(:,:),qd(:,:),plo(:,:)
    REAL, ALLOCATABLE :: plop(:,:),glo(:,:),fint(:),pqlo(:,:)
    REAL, ALLOCATABLE :: filo(:,:)
    REAL, ALLOCATABLE :: v0(:),vso(:,:),qlo(:,:),vrTmp(:)
    !     ..
    
    IF (atoms%jri(ntyp)>atoms%jmtd)  CALL juDFT_error("atoms%jri(ntyp).GT.atoms%jmtd",calledby ="sorad")

    ALLOCATE (p(atoms%jmtd,2),pd(atoms%jmtd,2),q(atoms%jmtd,2),plo(atoms%jmtd,2),fint(atoms%jmtd),&
              qlo(atoms%jmtd,2),plop(atoms%jmtd,2),qd(atoms%jmtd,2),v0(atoms%jmtd),vso(atoms%jmtd,2),vrTmp(atoms%jmtd) )
    p = 0.0
    pd = 0.0
    q = 0.0
    plo = 0.0
    fint = 0.0
    qlo = 0.0
    plop = 0.0
    qd = 0.0
    v0 = 0.0
    vso = 0.0
    vrTmp = 0.0
    

    DO l = 0,atoms%lmax(ntyp) 

       l_hia=.FALSE.
       DO iHia = atoms%n_u+1, atoms%n_u+atoms%n_hia
          IF(atoms%lda_u(iHia)%atomType.EQ.ntyp.AND.atoms%lda_u(iHia)%l.EQ.l) THEN
             l_hia=.TRUE.
          ENDIF
       ENDDO

       DO jspin = 1,input%jspins
          !TODO: here genMTBasis should be used
          vrTmp = enpara%vr(:,ntyp,jspin)

          IF(l_hia.AND.input%jspins.EQ.2) THEN
             IF(PRESENT(hub1data)) THEN
                IF(hub1data%l_performSpinavg) vrTmp = (vr(:,1)+vr(:,2))/2.0
             ENDIF
          ENDIF

          !--->    calculate normalized function at e: p and q 
          e = enpara%el0(l,ntyp,jspin)
          CALL radsra(e,l,vrTmp,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                      usdus%us(l,ntyp,jspin),usdus%dus(l,ntyp,jspin),&
                      nodeu,p(:,jspin),q(:,jspin))

          !--->    calculate orthogonal energy derivative at e : pd and qd
          CALL radsrd(e,l,vrTmp,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                      usdus%uds(l,ntyp,jspin),usdus%duds(l,ntyp,jspin),&
                      usdus%ddn(l,ntyp,jspin),noded,pd(:,jspin),qd(:,jspin),&
                      p(:,jspin),q(:,jspin),usdus%dus(l,ntyp,jspin))

       END DO     ! end of spin loop

       !---> in case of jspins=1
       IF (input%jspins.EQ.1) THEN
          DO iRad = 1,atoms%jri(ntyp)
             p(iRad,2) =  p(iRad,1)
             pd(iRad,2) = pd(iRad,1)
          ENDDO
       ENDIF

       !---> common spin-orbit integrant V   (average spin directions)
       !                                  SO
       v0(:) = 0.0
       IF (input%jspins.EQ.1) THEN
          v0(1:atoms%jri(ntyp)) = vr(1:atoms%jri(ntyp),1)
          e = enpara%el0(l,ntyp,1)
       ELSE
          DO iRad = 1,atoms%jri(ntyp)
             v0(iRad) = (vr(iRad,1) + vr(iRad,input%jspins)) / 2.0
          END DO
          e = (enpara%el0(l,ntyp,1) + enpara%el0(l,ntyp,input%jspins)) / 2.0
       END IF

       CALL sointg(ntyp,e,vr,v0,atoms,input,vso)
       IF (spav) THEN
          DO iRad = 1,atoms%jmtd
             vso(iRad,1)= (vso(iRad,1)+vso(iRad,2))/2.
             vso(iRad,2)= vso(iRad,1)
          ENDDO
       ENDIF

       !                        s       s'            .s       s'
       !-->  radial integrals <u  |V  |u  > = rsopp, <u  |V  |u  > = rsopdp etc.
       !                            SO                     SO

       IF (l.GT.0) THEN ! there is no spin-orbit for s-states
          DO iSpin = 1, 2
             DO jSpin = 1, 2
                ! Legacy rso* arrays removed from t_rsoc - commented out (see rsoc%rso / rad_matrix).
                !rsoc%rsopp(ntyp,l,iSpin,jSpin) = radso( p(:atoms%jri(ntyp),iSpin), p(:atoms%jri(ntyp),jSpin),(vso(:atoms%jri(ntyp),iSpin)+vso(:atoms%jri(ntyp),jSpin))*0.5,atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                !rsoc%rsopdp(ntyp,l,iSpin,jSpin) = radso(pd(:atoms%jri(ntyp),iSpin), p(:atoms%jri(ntyp),jSpin),(vso(:atoms%jri(ntyp),iSpin)+vso(:atoms%jri(ntyp),jSpin))*0.5,atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                !rsoc%rsoppd(ntyp,l,iSpin,jSpin) = radso( p(:atoms%jri(ntyp),iSpin),pd(:atoms%jri(ntyp),jSpin),(vso(:atoms%jri(ntyp),iSpin)+vso(:atoms%jri(ntyp),jSpin))*0.5,atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                !rsoc%rsopdpd(ntyp,l,iSpin,jSpin) = radso(pd(:atoms%jri(ntyp),iSpin),pd(:atoms%jri(ntyp),jSpin),(vso(:atoms%jri(ntyp),iSpin)+vso(:atoms%jri(ntyp),jSpin))*0.5,atoms%dx(ntyp),atoms%rmsh(1,ntyp))
             ENDDO
          ENDDO
       ENDIF ! l>0
       !
       !--->  Check for local orbitals with same l
       !
       DO ilo = 1, atoms%nlo(ntyp)
          IF (atoms%llo(ilo,ntyp).EQ.l) THEN

             DO jspin = 1,input%jspins
                vrTmp = vr(:,jspin)
                if (atoms%l_nonpolbas(ntyp)) vrTmp = (vr(:,1)+vr(:,2))/2.0
                e = enpara%ello0(ilo,ntyp,jspin)
                CALL radsra(e,l,vrtmp,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                            usdus%ulos(ilo,ntyp,jspin),usdus%dulos(ilo,ntyp,jspin),&
                            nodeu,plo(:,jspin),qlo(:,jspin))

                !+apw+lo
                IF (atoms%l_dulo(ilo,ntyp).OR.atoms%ulo_der(ilo,ntyp).GE.1) THEN !  calculate energy derivative (of order atoms%ulo_der) at e
                   ALLOCATE (glo(atoms%jmtd,2),pqlo(atoms%jmtd,2),filo(atoms%jmtd,2))
                   glo = 0.0
                   pqlo = 0.0
                   filo = 0.0
                   pqlo(1:atoms%jri(ntyp),1) = plo(1:atoms%jri(ntyp),jspin)
                   pqlo(1:atoms%jri(ntyp),2) = qlo(1:atoms%jri(ntyp),jspin)
                   iDeriv = atoms%ulo_der(ilo,ntyp)
                   IF(atoms%l_dulo(ilo,ntyp)) iDeriv = 1
                   CALL radsrdn(e,l,vrtmp,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                                usdus%ulos(ilo,ntyp,jspin),duds1,ddn1,noded,glo,filo,&!filo is a dummy array&
                                pqlo,usdus%dulos(ilo,ntyp,jspin),iDeriv)
                   ddn1 = SQRT(ddn1)
                   IF(atoms%l_dulo(ilo,ntyp)) ddn1 = 1.0
                   plo(1:atoms%jri(ntyp),jspin) = glo(1:atoms%jri(ntyp),1) / ddn1
                   qlo(1:atoms%jri(ntyp),jspin) = glo(1:atoms%jri(ntyp),2) / ddn1
                   usdus%dulos(ilo,ntyp,jspin) = duds1 / ddn1
                   usdus%ulos(ilo,ntyp,jspin) = usdus%ulos(ilo,ntyp,jspin) / ddn1
                   DEALLOCATE (glo,pqlo,filo)
                ENDIF
                !-apw+lo
             ENDDO

             IF (input%jspins.EQ.1) THEN
                plo(1:atoms%jri(ntyp),2) = plo(1:atoms%jri(ntyp),1)
                e = (enpara%ello0(ilo,ntyp,1) + enpara%el0(l,ntyp,1)) / 2.0
             ELSE
                e = (enpara%ello0(ilo,ntyp,1) +  enpara%ello0(ilo,ntyp,input%jspins) + &
                     enpara%el0(l,ntyp,1) + enpara%el0(l,ntyp,input%jspins)) / 4.0
             END IF
             CALL sointg(ntyp,e,vr,v0,atoms,input, vso)
             IF (spav) THEN
                DO iRad = 1, atoms%jmtd
                   vso(iRad,1) = (vso(iRad,1) + vso(iRad,2)) / 2.0
                   vso(iRad,2) = vso(iRad,1)
                ENDDO
             ENDIF

             DO iSpin = 1, 2
                DO jSpin = 1, 2
                   ! Legacy rso* arrays removed from t_rsoc - commented out (see rsoc%rso / rad_matrix).
                   !rsoc%rsoplop (ntyp,ilo,iSpin,jSpin) = radso(plo(:atoms%jri(ntyp),iSpin),p (:atoms%jri(ntyp),jSpin),vso(:atoms%jri(ntyp),iSpin),atoms%dx(ntyp),atoms%rmsh(1,ntyp))
                   !rsoc%rsoplopd(ntyp,ilo,iSpin,jSpin) = radso(plo(:atoms%jri(ntyp),iSpin),pd(:atoms%jri(ntyp),jSpin),vso(:atoms%jri(ntyp),iSpin),atoms%dx(ntyp), atoms%rmsh(1,ntyp))
                   !rsoc%rsopplo (ntyp,ilo,iSpin,jSpin) = radso(p (:atoms%jri(ntyp),iSpin),plo(:atoms%jri(ntyp),jSpin),vso(:atoms%jri(ntyp),iSpin),atoms%dx(ntyp), atoms%rmsh(1,ntyp))
                   !rsoc%rsopdplo(ntyp,ilo,iSpin,jSpin) = radso(pd(:atoms%jri(ntyp),iSpin),plo(:atoms%jri(ntyp),jSpin),vso(:atoms%jri(ntyp),iSpin),atoms%dx(ntyp), atoms%rmsh(1,ntyp))
                ENDDO
             ENDDO

             DO iSpin = 1,input%jspins
                fint(:) = plo(:,iSpin) *  p(:,iSpin) + qlo(:,iSpin) *  q(:,iSpin)
                CALL intgr0(fint,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%uulon(ilo,ntyp,iSpin))
                fint(:) = plo(:,iSpin) * pd(:,iSpin) + qlo(:,iSpin) * qd(:,iSpin)
                CALL intgr0(fint,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),usdus%dulon(ilo,ntyp,iSpin))
             ENDDO

             DO ilop = 1, atoms%nlo(ntyp)
                IF (atoms%llo(ilop,ntyp).EQ.l) THEN

                   DO jspin = 1,input%jspins
                      vrTmp = vr(:,jspin)
                      IF (atoms%l_nonpolbas(ntyp)) vrTmp = (vr(:,1) + vr(:,2)) / 2.0
                      e = enpara%ello0(ilop,ntyp,jspin)
                      CALL radsra(e,l,vrtmp,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                                  ulops,dulops,nodeu,plop(:,jspin),q(:,1))
                      !+apw+lo
                      IF (atoms%l_dulo(ilop,ntyp).OR.atoms%ulo_der(ilop,ntyp).GE.1) THEN ! calculate orthogonal energy derivative at e
                         ALLOCATE (glo(atoms%jmtd,2),pqlo(atoms%jmtd,2),filo(atoms%jmtd,2))
                         glo = 0.0
                         pqlo = 0.0
                         filo = 0.0
                         pqlo(1:atoms%jri(ntyp),1) = plop(1:atoms%jri(ntyp),jspin)
                         pqlo(1:atoms%jri(ntyp),2) = q(1:atoms%jri(ntyp),1)
                         iDeriv = atoms%ulo_der(ilop,ntyp)
                         IF(atoms%l_dulo(ilop,ntyp)) iDeriv = 1
                         CALL radsrdn(e,l,vrtmp,atoms%rmsh(1,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),c_light(1.0),&
                                      ulops,duds1,ddn1,noded,glo,filo,&!filo is a dummy array&
                                      pqlo,dulops,iDeriv)
                         ddn1 = SQRT(ddn1)
                         IF(atoms%l_dulo(ilop,ntyp)) ddn1 = 1.0
                         plop(1:atoms%jri(ntyp),jspin) = glo(1:atoms%jri(ntyp),1) / ddn1
                         DEALLOCATE (glo,pqlo,filo)
                      ENDIF
                      !-apw+lo
                   ENDDO

                   IF (input%jspins.EQ.1) THEN
                      plop(1:atoms%jri(ntyp),2) = plop(1:atoms%jri(ntyp),1)
                      e = (enpara%ello0(ilo,ntyp,1) + enpara%ello0(ilop,ntyp,1) ) / 2.0
                   ELSE
                      e = (enpara%ello0(ilo,ntyp,1) +  enpara%ello0(ilo,ntyp,input%jspins) +  &
                           enpara%ello0(ilop,ntyp,1) + enpara%ello0(ilop,ntyp,input%jspins)) / 4.0
                   END IF
                   CALL sointg(ntyp,e,vr,v0,atoms,input, vso)
                   IF (spav) THEN
                      DO iRad = 1,atoms%jmtd
                         vso(iRad,1)= (vso(iRad,1) + vso(iRad,2)) / 2.0
                         vso(iRad,2)= vso(iRad,1)
                      ENDDO
                   ENDIF

                   DO iSpin = 1, 2
                      DO jSpin = 1, 2
                         ! Legacy rso* arrays removed from t_rsoc - commented out (see rsoc%rso / rad_matrix).
                         !rsoc%rsoploplop(ntyp,ilo,ilop,iSpin,jSpin) = &
                         !     radso(plo(:atoms%jri(ntyp),iSpin),plop(:atoms%jri(ntyp),jSpin),vso(:atoms%jri(ntyp),iSpin),atoms%dx(ntyp),atoms%rmsh(1,ntyp))
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
