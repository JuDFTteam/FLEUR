MODULE m_cdnmt
  !***********************************************************************
  !     This subroutine calculates the spherical and non-spherical charge-
  !     density and the orbital moment inside the muffin-tin spheres.
  !     Philipp Kurz 2000-02-03
  !***********************************************************************
CONTAINS
  SUBROUTINE cdnmt(jspd,atoms,sphhar,llpd, noco,l_fmpl,jsp_start,jsp_end, epar,&
       ello,vr,uu,du,dd,uunmt,udnmt,dunmt,ddnmt, usdus,uloulopn,aclo,bclo,cclo,&
       acnmt,bcnmt,ccnmt, orb,orbl,orblo,mt21,lo21,uloulopn21,uloulop21, uunmt21,&
       ddnmt21,udnmt21,dunmt21, chmom,clmom, qa21,rho)
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
    REAL, INTENT    (IN) :: epar(0:atoms%lmaxd,atoms%ntypd,jspd)
    REAL, INTENT    (IN) :: vr(atoms%jmtd,atoms%ntypd,jspd)
    REAL, INTENT    (IN) :: ello(atoms%nlod,atoms%ntypd,jspd)
    REAL, INTENT    (IN) :: uloulopn(atoms%nlod,atoms%nlod,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: aclo(atoms%nlod,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: bclo(atoms%nlod,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: cclo(atoms%nlod,atoms%nlod,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: acnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: bcnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: ccnmt(atoms%nlod,atoms%nlod,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: uu(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: du(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: dd(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: uunmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: udnmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: dunmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: ddnmt(0:llpd,sphhar%nlhd,atoms%ntypd,jsp_start:jsp_end)
    REAL, INTENT    (IN) :: uloulopn21(atoms%nlod,atoms%nlod,atoms%ntypd)
    COMPLEX, INTENT (IN) :: uloulop21(atoms%nlod,atoms%nlod,atoms%ntypd)
    COMPLEX, INTENT (IN) :: ddnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd)  
    COMPLEX, INTENT (IN) :: dunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd)  
    COMPLEX, INTENT (IN) :: udnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd)  
    COMPLEX, INTENT (IN) :: uunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntypd)  
    REAL, INTENT   (OUT) :: chmom(atoms%ntypd,jspd),clmom(3,atoms%ntypd,jspd)
    REAL, INTENT (INOUT) :: rho(:,0:,:,:)!(toms%jmtd,0:sphhar%nlhd,atoms%ntypd,jspd)
    COMPLEX, INTENT(INOUT) :: qa21(atoms%ntypd)
    TYPE (t_orb),  INTENT (IN) :: orb(0:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end)
    TYPE (t_orbl), INTENT (IN) :: orbl(atoms%nlod,-atoms%llod:atoms%llod,atoms%ntypd,jsp_start:jsp_end)
    TYPE (t_orblo),INTENT (IN) :: orblo(atoms%nlod,atoms%nlod,-atoms%llod:atoms%llod,atoms%ntypd,jsp_start:jsp_end)
    TYPE (t_mt21), INTENT (IN) :: mt21(0:atoms%lmaxd,atoms%ntypd)
    TYPE (t_lo21), INTENT (IN) :: lo21(atoms%nlod,atoms%ntypd)
    !     ..
    !     .. Local Scalars ..
    INTEGER itype,na,nd,l,lp,llp ,lh,j,ispin,noded,nodeu
    INTEGER ilo,ilop
    REAL s,wronk,sumlm,qmtt
    COMPLEX cs
    !     ..
    !     .. Local Arrays ..
    REAL qmtl(0:atoms%lmaxd),qmtllo(0:atoms%lmaxd)
    CHARACTER(LEN=20) :: attributes(6)

    !     ..
    !     .. Allocatable Arrays ..
    REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
    COMPLEX, ALLOCATABLE :: rho21(:,:,:)
    !
    IF (noco%l_mperp) THEN
       ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jspd),g(atoms%jmtd,2,0:atoms%lmaxd,jspd) )
       ALLOCATE ( usdus%us(0:atoms%lmaxd,atoms%ntypd,jspd),usdus%uds(0:atoms%lmaxd,atoms%ntypd,jspd) )
       ALLOCATE ( usdus%dus(0:atoms%lmaxd,atoms%ntypd,jspd),usdus%duds(0:atoms%lmaxd,atoms%ntypd,jspd) )
       ALLOCATE ( usdus%ddn(0:atoms%lmaxd,atoms%ntypd,jspd) )
       IF (l_fmpl) THEN
          ALLOCATE ( rho21(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd) )
          rho21(:,:,:) = cmplx(0.0,0.0)
       ENDIF
    ELSE
       ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
       ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jsp_start:jsp_end) )
       ALLOCATE (   usdus%us(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
       ALLOCATE (  usdus%uds(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
       ALLOCATE (  usdus%dus(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
       ALLOCATE ( usdus%duds(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
       ALLOCATE (  usdus%ddn(0:atoms%lmaxd,atoms%ntypd,jsp_start:jsp_end) )
    ENDIF
    WRITE (6,FMT=8000)
    WRITE (16,FMT=8000)
8000 FORMAT (/,5x,'l-like charge',/,t6,'atom',t15,'s',t24,'p',&
         &     t33,'d',t42,'f',t51,'total')
    CALL timestart("cdnmt")

    CALL openXMLElementNoAttributes('mtCharges')
    na = 1
    DO itype = 1,atoms%ntype
       !--->    spherical component
       DO ispin = jsp_start,jsp_end
          DO l = 0,atoms%lmax(itype)
             CALL radfun(l,itype,ispin,epar(l,itype,ispin),vr(1,itype,ispin),atoms,&
                  f(1,1,l,ispin),g(1,1,l,ispin),usdus, nodeu,noded,wronk)
             DO j = 1,atoms%jri(itype)
                s = uu(l,itype,ispin)*( f(j,1,l,ispin)*f(j,1,l,ispin)+f(j,2,l,ispin)*f(j,2,l,ispin) )&
                     +   dd(l,itype,ispin)*( g(j,1,l,ispin)*g(j,1,l,ispin)+g(j,2,l,ispin)*g(j,2,l,ispin) )&
                     + 2*du(l,itype,ispin)*( f(j,1,l,ispin)*g(j,1,l,ispin)+f(j,2,l,ispin)*g(j,2,l,ispin) )
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
               uloulopn(1,1,itype,ispin),usdus%dulon(1,itype,ispin),&
               usdus%uulon(1,itype,ispin),ello(1,itype,ispin),&
               vr(1,itype,ispin),aclo(1,itype,ispin),bclo(1,itype,ispin),&
               cclo(1,1,itype,ispin),acnmt(0,1,1,itype,ispin),&
               bcnmt(0,1,1,itype,ispin),ccnmt(1,1,1,itype,ispin),&
               f(1,1,0,ispin),g(1,1,0,ispin),&
               rho(:,0:,itype,ispin),qmtllo)


          !--->       l-decomposed density for each atom type
          qmtt = 0.
          DO l = 0,atoms%lmax(itype)
             qmtl(l) = ( uu(l,itype,ispin)+dd(l,itype,ispin)&
                  &              *usdus%ddn(l,itype,ispin) )/atoms%neq(itype) + qmtllo(l)
             qmtt = qmtt + qmtl(l)
          END DO
          chmom(itype,ispin) = qmtt
          WRITE (6,FMT=8100) itype, (qmtl(l),l=0,3),qmtt
          WRITE (16,FMT=8100) itype, (qmtl(l),l=0,3),qmtt
8100      FORMAT (' -->',i2,2x,4f9.5,2x,f9.5)

          attributes = ''
          WRITE(attributes(1),'(i0)') itype
          WRITE(attributes(2),'(f15.10)') qmtt
          WRITE(attributes(3),'(f15.10)') qmtl(0)
          WRITE(attributes(4),'(f15.10)') qmtl(1)
          WRITE(attributes(5),'(f15.10)') qmtl(2)
          WRITE(attributes(6),'(f15.10)') qmtl(3)
          CALL writeXMLElementForm('mtCharge',(/'atomType','total   ','s       ','p       ','d       ','f       '/),attributes,&
                                   reshape((/8,5,1,1,1,1,6,15,15,15,15,15/),(/6,2/)))

          !+soc
          !--->       spherical angular component
          IF (noco%l_soc) THEN
             CALL orbmom2(&
                  atoms,itype,&
                  usdus%ddn(0,itype,ispin),&
                  orb(0,-atoms%lmaxd,itype,ispin),usdus%uulon(1,itype,ispin),&
                  usdus%dulon(1,itype,ispin),uloulopn(1,1,itype,ispin),&
                  orbl(1,-atoms%llod,itype,ispin),orblo(1,1,-atoms%llod,itype,ispin),&
                  clmom(1,itype,ispin))!keep
          ENDIF
          !-soc
          !--->       non-spherical components
          nd = atoms%ntypsy(na)
          DO lh = 1,sphhar%nlh(nd)
             DO l = 0,atoms%lmax(itype)
                DO lp = 0,l
                   llp = (l* (l+1))/2 + lp
                   DO j = 1,atoms%jri(itype)
                      s = uunmt(llp,lh,itype,ispin)*( &
                           f(j,1,l,ispin)*f(j,1,lp,ispin)+ f(j,2,l,ispin)*f(j,2,lp,ispin) )&
                           + ddnmt(llp,lh,itype,ispin)*(g(j,1,l,ispin)*g(j,1,lp,ispin)&
                           + g(j,2,l,ispin)*g(j,2,lp,ispin) )&
                           + udnmt(llp,lh,itype,ispin)*(f(j,1,l,ispin)*g(j,1,lp,ispin)&
                           + f(j,2,l,ispin)*g(j,2,lp,ispin) )&
                           + dunmt(llp,lh,itype,ispin)*(g(j,1,l,ispin)*f(j,1,lp,ispin)&
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
                  mt21(l,itype)%uu * mt21(l,itype)%uun +&
                  mt21(l,itype)%ud * mt21(l,itype)%udn +&
                  mt21(l,itype)%du * mt21(l,itype)%dun +&
                  mt21(l,itype)%dd * mt21(l,itype)%ddn )/atoms%neq(itype)
          ENDDO
          DO ilo = 1, atoms%nlo(itype)
             qa21(itype) = qa21(itype) + conjg(&
                  lo21(ilo,itype)%ulou * lo21(ilo,itype)%uloun +&
                  lo21(ilo,itype)%ulod * lo21(ilo,itype)%ulodn +&
                  lo21(ilo,itype)%uulo * lo21(ilo,itype)%uulon +&
                  lo21(ilo,itype)%dulo * lo21(ilo,itype)%dulon )/&
                  atoms%neq(itype)
             DO ilop = 1, atoms%nlo(itype)
                qa21(itype) = qa21(itype) + conjg(&
                     uloulop21(ilo,ilop,itype) *&
                     uloulopn21(ilo,ilop,itype) )/atoms%neq(itype)
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
                   cs = mt21(l,itype)%uu*( f(j,1,l,2)*f(j,1,l,1) +f(j,2,l,2)*f(j,2,l,1) )&
                        + mt21(l,itype)%ud*( f(j,1,l,2)*g(j,1,l,1) +f(j,2,l,2)*g(j,2,l,1) )&
                        + mt21(l,itype)%du*( g(j,1,l,2)*f(j,1,l,1) +g(j,2,l,2)*f(j,2,l,1) )&
                        + mt21(l,itype)%dd*( g(j,1,l,2)*g(j,1,l,1) +g(j,2,l,2)*g(j,2,l,1) )
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
                         cs = uunmt21(llp,lh,itype)*(f(j,1,lp,2)*f(j,1,l,1)&
                              + f(j,2,lp,2)*f(j,2,l,1) )+ udnmt21(llp,lh,itype)*(f(j,1,lp,2)*g(j,1,l,1)&
                              + f(j,2,lp,2)*g(j,2,l,1) )+ dunmt21(llp,lh,itype)*(g(j,1,lp,2)*f(j,1,l,1)&
                              + g(j,2,lp,2)*f(j,2,l,1) )+ ddnmt21(llp,lh,itype)*(g(j,1,lp,2)*g(j,1,l,1)&
                              + g(j,2,lp,2)*g(j,2,l,1) )
                         rho21(j,lh,itype)= rho21(j,lh,itype)+ conjg(cs)/atoms%neq(itype)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO

          ENDIF ! l_fmpl
       ENDIF ! noco%l_mperp

       na = na + atoms%neq(itype)
    ENDDO ! end of loop over atom types
    CALL closeXMLElement('mtCharges')
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

    DEALLOCATE ( f,g)

  END SUBROUTINE cdnmt
END MODULE m_cdnmt
