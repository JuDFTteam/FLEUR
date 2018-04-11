! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnmt
  !***********************************************************************
  !     This subroutine calculates the spherical and non-spherical charge-
  !     density and the orbital moment inside the muffin-tin spheres.
  !     Philipp Kurz 2000-02-03
  !***********************************************************************
CONTAINS
  SUBROUTINE cdnmt(jspd,atoms,sphhar,llpd, noco,l_fmpl,jsp_start,jsp_end, epar,&
       ello,vr,denCoeffs,usdus,&
       orb,denCoeffsOffdiag,&
       chmom,clmom, qa21,rho)
    use m_constants,only: sfp_const
    USE m_rhosphnlo
    USE m_radfun
    USE m_orbmom2
    USE m_types
    USE m_xmlOutput
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(INOUT):: usdus !in fact only the lo part is intent(in)
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: llpd 
    INTEGER, INTENT (IN) :: jsp_start,jsp_end,jspd
    LOGICAL, INTENT (IN) :: l_fmpl
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT    (IN) :: epar(0:atoms%lmaxd,atoms%ntype,jspd)
    REAL, INTENT    (IN) :: vr(atoms%jmtd,atoms%ntype,jspd)
    REAL, INTENT    (IN) :: ello(atoms%nlod,atoms%ntype,jspd)
    REAL, INTENT   (OUT) :: chmom(atoms%ntype,jspd),clmom(3,atoms%ntype,jspd)
    REAL, INTENT (INOUT) :: rho(:,0:,:,:)!(toms%jmtd,0:sphhar%nlhd,atoms%ntype,jspd)
    COMPLEX, INTENT(INOUT) :: qa21(atoms%ntype)
    TYPE (t_orb),              INTENT(IN) :: orb
    TYPE (t_denCoeffs),        INTENT(IN) :: denCoeffs
    TYPE (t_denCoeffsOffdiag), INTENT(IN) :: denCoeffsOffdiag
    !     ..
    !     .. Local Scalars ..
    INTEGER itype,na,nd,l,lp,llp ,lh,j,ispin,noded,nodeu
    INTEGER ilo,ilop,i
    REAL s,wronk,sumlm,qmtt
    COMPLEX cs
    !     ..
    !     .. Local Arrays ..
    REAL qmtl(0:atoms%lmaxd,jspd,atoms%ntype),qmtllo(0:atoms%lmaxd)
    CHARACTER(LEN=20) :: attributes(6)

    !     ..
    !     .. Allocatable Arrays ..
    REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
    COMPLEX, ALLOCATABLE :: rho21(:,:,:)
    !

    CALL timestart("cdnmt")

    IF (noco%l_mperp) THEN
       IF (l_fmpl) THEN
          ALLOCATE ( rho21(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) )
          rho21(:,:,:) = cmplx(0.0,0.0)
       ENDIF
    ENDIF

    !$OMP PARALLEL DEFAULT(none) &

    !$OMP SHARED(usdus,rho,chmom,clmom,qa21,rho21,qmtl) &
    !$OMP SHARED(atoms,jsp_start,jsp_end,epar,vr,denCoeffs,sphhar,ello)&
    !$OMP SHARED(orb,noco,l_fmpl,denCoeffsOffdiag,jspd)&
    !$OMP PRIVATE(itype,na,ispin,l,f,g,nodeu,noded,wronk,i,j,s,qmtllo,qmtt,nd,lh,lp,llp,cs)
    IF (noco%l_mperp) THEN
       ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jspd),g(atoms%jmtd,2,0:atoms%lmaxd,jspd) )
    ELSE
       ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
       ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
    ENDIF

    qmtl = 0
    
    !$OMP DO
    DO itype = 1,atoms%ntype
       na = 1
       DO i = 1, itype - 1
          na = na + atoms%neq(i)
       ENDDO
       !--->    spherical component
       DO ispin = jsp_start,jsp_end
          DO l = 0,atoms%lmax(itype)
             CALL radfun(l,itype,ispin,epar(l,itype,ispin),vr(1,itype,ispin),atoms,&
                  f(1,1,l,ispin),g(1,1,l,ispin),usdus, nodeu,noded,wronk)
             DO j = 1,atoms%jri(itype)
                s = denCoeffs%uu(l,itype,ispin)*( f(j,1,l,ispin)*f(j,1,l,ispin)+f(j,2,l,ispin)*f(j,2,l,ispin) )&
                     +   denCoeffs%dd(l,itype,ispin)*( g(j,1,l,ispin)*g(j,1,l,ispin)+g(j,2,l,ispin)*g(j,2,l,ispin) )&
                     + 2*denCoeffs%du(l,itype,ispin)*( f(j,1,l,ispin)*g(j,1,l,ispin)+f(j,2,l,ispin)*g(j,2,l,ispin) )
                rho(j,0,itype,ispin) = rho(j,0,itype,ispin)+ s/(atoms%neq(itype)*sfp_const)
             ENDDO
          ENDDO

          !--->       add the contribution of the local orbitals and flapw - lo
          !--->       cross-terms to rho, qmtl. the latter are stored in
          !--->       qmtllo. initialize qmtllo
          DO l = 0,atoms%lmaxd
             qmtllo(l) = 0.0
          END DO

          CALL rhosphnlo(itype,atoms,sphhar,&
               usdus%uloulopn(1,1,itype,ispin),usdus%dulon(1,itype,ispin),&
               usdus%uulon(1,itype,ispin),ello(1,itype,ispin),&
               vr(1,itype,ispin),denCoeffs%aclo(1,itype,ispin),denCoeffs%bclo(1,itype,ispin),&
               denCoeffs%cclo(1,1,itype,ispin),denCoeffs%acnmt(0,1,1,itype,ispin),&
               denCoeffs%bcnmt(0,1,1,itype,ispin),denCoeffs%ccnmt(1,1,1,itype,ispin),&
               f(1,1,0,ispin),g(1,1,0,ispin),&
               rho(:,0:,itype,ispin),qmtllo)


          !--->       l-decomposed density for each atom type
          qmtt = 0.0
          DO l = 0,atoms%lmax(itype)
             qmtl(l,ispin,itype) = ( denCoeffs%uu(l,itype,ispin)+denCoeffs%dd(l,itype,ispin)&
                  &              *usdus%ddn(l,itype,ispin) )/atoms%neq(itype) + qmtllo(l)
             qmtt = qmtt + qmtl(l,ispin,itype)
          END DO
          chmom(itype,ispin) = qmtt

          !+soc
          !--->       spherical angular component
          IF (noco%l_soc) THEN
             CALL orbmom2(atoms,itype,ispin,usdus%ddn(0,itype,ispin),&
                          orb,usdus%uulon(1,itype,ispin),usdus%dulon(1,itype,ispin),&
                          usdus%uloulopn(1,1,itype,ispin),clmom(1,itype,ispin))!keep
          ENDIF
          !-soc
          !--->       non-spherical components
          nd = atoms%ntypsy(na)
          DO lh = 1,sphhar%nlh(nd)
             DO l = 0,atoms%lmax(itype)
                DO lp = 0,l
                   llp = (l* (l+1))/2 + lp
                   DO j = 1,atoms%jri(itype)
                      s = denCoeffs%uunmt(llp,lh,itype,ispin)*( &
                           f(j,1,l,ispin)*f(j,1,lp,ispin)+ f(j,2,l,ispin)*f(j,2,lp,ispin) )&
                           + denCoeffs%ddnmt(llp,lh,itype,ispin)*(g(j,1,l,ispin)*g(j,1,lp,ispin)&
                           + g(j,2,l,ispin)*g(j,2,lp,ispin) )&
                           + denCoeffs%udnmt(llp,lh,itype,ispin)*(f(j,1,l,ispin)*g(j,1,lp,ispin)&
                           + f(j,2,l,ispin)*g(j,2,lp,ispin) )&
                           + denCoeffs%dunmt(llp,lh,itype,ispin)*(g(j,1,l,ispin)*f(j,1,lp,ispin)&
                           + g(j,2,l,ispin)*f(j,2,lp,ispin) )
                      rho(j,lh,itype,ispin) = rho(j,lh,itype,ispin)+ s/atoms%neq(itype)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO ! end of spin loop (ispin = jsp_start,jsp_end)

       IF (noco%l_mperp) THEN

          !--->      calculate off-diagonal integrated density
          DO l = 0,atoms%lmax(itype)
             qa21(itype) = qa21(itype) + conjg(&
                  denCoeffsOffdiag%uu21(l,itype) * denCoeffsOffdiag%uu21n(l,itype) +&
                  denCoeffsOffdiag%ud21(l,itype) * denCoeffsOffdiag%ud21n(l,itype) +&
                  denCoeffsOffdiag%du21(l,itype) * denCoeffsOffdiag%du21n(l,itype) +&
                  denCoeffsOffdiag%dd21(l,itype) * denCoeffsOffdiag%dd21n(l,itype) )/atoms%neq(itype)
          ENDDO
          DO ilo = 1, atoms%nlo(itype)
             qa21(itype) = qa21(itype) + conjg(&
                  denCoeffsOffdiag%ulou21(ilo,itype) * denCoeffsOffdiag%ulou21n(ilo,itype) +&
                  denCoeffsOffdiag%ulod21(ilo,itype) * denCoeffsOffdiag%ulod21n(ilo,itype) +&
                  denCoeffsOffdiag%uulo21(ilo,itype) * denCoeffsOffdiag%uulo21n(ilo,itype) +&
                  denCoeffsOffdiag%dulo21(ilo,itype) * denCoeffsOffdiag%dulo21n(ilo,itype) )/&
                  atoms%neq(itype)
             DO ilop = 1, atoms%nlo(itype)
                qa21(itype) = qa21(itype) + conjg(&
                     denCoeffsOffdiag%uloulop21(ilo,ilop,itype) *&
                     denCoeffsOffdiag%uloulop21n(ilo,ilop,itype) )/atoms%neq(itype)
             ENDDO
          ENDDO

          IF (l_fmpl) THEN
             !--->        the following part can be used to calculate the full magnet.
             !--->        density without the atomic sphere approximation for the
             !--->        magnet. density, e.g. for plotting.
             !--->        calculate off-diagonal part of the density matrix
             !--->        spherical component
             DO l = 0,atoms%lmax(itype)
                DO j = 1,atoms%jri(itype)
                   cs = denCoeffsOffdiag%uu21(l,itype)*( f(j,1,l,2)*f(j,1,l,1) +f(j,2,l,2)*f(j,2,l,1) )&
                        + denCoeffsOffdiag%ud21(l,itype)*( f(j,1,l,2)*g(j,1,l,1) +f(j,2,l,2)*g(j,2,l,1) )&
                        + denCoeffsOffdiag%du21(l,itype)*( g(j,1,l,2)*f(j,1,l,1) +g(j,2,l,2)*f(j,2,l,1) )&
                        + denCoeffsOffdiag%dd21(l,itype)*( g(j,1,l,2)*g(j,1,l,1) +g(j,2,l,2)*g(j,2,l,1) )
                   rho21(j,0,itype) = rho21(j,0,itype)+ conjg(cs)/(atoms%neq(itype)*sfp_const)
                ENDDO
             ENDDO

             !--->        non-spherical components
             nd = atoms%ntypsy(na)
             DO lh = 1,sphhar%nlh(nd)
                DO l = 0,atoms%lmax(itype)
                   DO lp = 0,atoms%lmax(itype)
                      llp = lp*(atoms%lmax(itype)+1)+l+1
                      DO j = 1,atoms%jri(itype)
                         cs = denCoeffsOffdiag%uunmt21(llp,lh,itype)*(f(j,1,lp,2)*f(j,1,l,1)&
                              + f(j,2,lp,2)*f(j,2,l,1) )+ denCoeffsOffdiag%udnmt21(llp,lh,itype)*(f(j,1,lp,2)*g(j,1,l,1)&
                              + f(j,2,lp,2)*g(j,2,l,1) )+ denCoeffsOffdiag%dunmt21(llp,lh,itype)*(g(j,1,lp,2)*f(j,1,l,1)&
                              + g(j,2,lp,2)*f(j,2,l,1) )+ denCoeffsOffdiag%ddnmt21(llp,lh,itype)*(g(j,1,lp,2)*g(j,1,l,1)&
                              + g(j,2,lp,2)*g(j,2,l,1) )
                         rho21(j,lh,itype)= rho21(j,lh,itype)+ conjg(cs)/atoms%neq(itype)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

          ENDIF ! l_fmpl
       ENDIF ! noco%l_mperp

    ENDDO ! end of loop over atom types
    !$OMP END DO
    DEALLOCATE ( f,g)
    !$OMP END PARALLEL

    WRITE (6,FMT=8000)
    WRITE (16,FMT=8000)
8000 FORMAT (/,5x,'l-like charge',/,t6,'atom',t15,'s',t24,'p',&
         &     t33,'d',t42,'f',t51,'total')

    DO itype = 1,atoms%ntype
       DO ispin = jsp_start,jsp_end
          WRITE ( 6,FMT=8100) itype, (qmtl(l,ispin,itype),l=0,3),chmom(itype,ispin)
          WRITE (16,FMT=8100) itype, (qmtl(l,ispin,itype),l=0,3),chmom(itype,ispin)
8100      FORMAT (' -->',i3,2x,4f9.5,2x,f9.5)
          attributes = ''
          WRITE(attributes(1),'(i0)') itype
          WRITE(attributes(2),'(f12.7)') chmom(itype,ispin)
          WRITE(attributes(3),'(f12.7)') qmtl(0,ispin,itype)
          WRITE(attributes(4),'(f12.7)') qmtl(1,ispin,itype)
          WRITE(attributes(5),'(f12.7)') qmtl(2,ispin,itype)
          WRITE(attributes(6),'(f12.7)') qmtl(3,ispin,itype)
          CALL writeXMLElementForm('mtCharge',(/'atomType','total   ','s       ','p       ','d       ','f       '/),attributes,&
                                   reshape((/8,5,1,1,1,1,6,12,12,12,12,12/),(/6,2/)))
       ENDDO
    ENDDO

    CALL timestop("cdnmt")

    !---> for testing: to plot the offdiag. part of the density matrix it
    !---> is written to the file rhomt21. This file can read in pldngen.
    IF (l_fmpl) THEN
       OPEN (26,file='rhomt21',form='unformatted',status='unknown')
       WRITE (26) rho21
       CLOSE (26)
       DEALLOCATE ( rho21 )
    ENDIF
    !---> end of test output

  END SUBROUTINE cdnmt
END MODULE m_cdnmt
