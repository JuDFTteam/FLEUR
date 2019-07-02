MODULE m_potdis
CONTAINS
  SUBROUTINE potdis(stars,vacuum,atoms,sphhar, input,cell,sym)
    !
    !     *****************************************************
    !     calculates the root-mean-square distance between input
    !     and output potentials (averaged over unit cell volume)
    !                                 based on code by   c.l.fu
    !     *****************************************************
    !
    USE m_intgr, ONLY : intgr3, intgz0
    USE m_constants, ONLY : fpi_const
    USE m_loddop
    USE m_cfft
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     .. Local Scalars ..
    REAL fact,facv,rhs,sumis,sumz
    COMPLEX phase
    INTEGER i,i1,i2,i3,id2,id3,io,ip,iter,ivac,j,k1,k2,k3,lh,n,&
                 nk12,npz,nt,num,na
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    COMPLEX rhpw(stars%ng3,2,2),rhv1(vacuum%nmzxyd,stars%ng2-1,2,2,2)
    REAL rhsp(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,2,2),rhv0(vacuum%nmzd,2,2,2)
    REAL dis(2),disz(vacuum%nmzd),rh(atoms%jmtd)
    REAL af3((2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1),input%jspins),&
         bf3((2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1),input%jspins), &
         bf2((2*stars%mx1+1)*(2*stars%mx2+1),input%jspins),&
         af2((2*stars%mx1+1)*(2*stars%mx2+1),input%jspins)
    REAL pdis(0:4,0:atoms%ntype,2)
    !     ..
    !

    dis = 0.
    pdis=0.0

    tail = .TRUE.
    npz = vacuum%nmz + 1
    nt = (2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)
    nk12 = (2*stars%mx1+1)*(2*stars%mx2+1)
    fact = cell%omtil/REAL(nt)
    facv = 1.0
    !     ---> reload potentials
    OPEN (9,file='nrp',form='unformatted',status='unknown')
    REWIND 9
    DO  io = 1,2
       CALL loddop(stars,vacuum,atoms,sphhar, input,sym,&
            9, iter,rhsp(:,0:,:,:,io),rhpw(:,:,io), rhv0(:,:,:,io),rhv1(:,:,:,:,io))
    ENDDO
    CLOSE (9)
    IF (input%jspins.EQ.1) THEN
       !       ---> total potential difference
       CALL priv_cdndif(1,1,1,1,1,2,1.0,-1.0,input%film, rhsp,rhpw,rhv0,rhv1)
    ELSE
       !       ---> total input potential
       CALL priv_cdndif(1,1,2,1,1,1,1.0,1.0,input%film, rhsp,rhpw,rhv0,rhv1)
       !       ---> total output potential
       CALL priv_cdndif(1,1,2,2,2,2,1.0,1.0,input%film, rhsp,rhpw,rhv0,rhv1)
       !       ---> input 'spin' potential
       CALL priv_cdndif(2,1,2,1,1,1,1.0,-2.0,input%film, rhsp,rhpw,rhv0,rhv1)
       !       ---> output 'spin' potential
       CALL priv_cdndif(2,1,2,2,2,2,1.0,-2.0,input%film, rhsp,rhpw,rhv0,rhv1)
       !       ---> potential difference
       CALL priv_cdndif(1,1,1,1,1,2,1.0,-1.0,input%film, rhsp,rhpw,rhv0,rhv1)
       !       ---> spin potential difference
       CALL priv_cdndif(2,2,2,1,1,2,1.0,-1.0,input%film, rhsp,rhpw,rhv0,rhv1)
    END IF
    DO  num = 1,input%jspins
       !     ----> m.t. part
       na = 1
       DO n = 1,atoms%ntype
          rh(:atoms%jri(n))=rhsp(:atoms%jri(n),0,n,num,1)*rhsp(:atoms%jri(n),0,n,num,1) *fpi_const
          CALL intgr3(rh,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),pdis(0,n,num))
          pdis(4,n,num)=pdis(0,n,num)
          DO lh = 1,sphhar%nlh(sym%ntypsy(na))
             rh(:atoms%jri(n))=rh(:atoms%jri(n))+rhsp(:atoms%jri(n),lh,n,num,1)*&
                  rhsp(:atoms%jri(n),lh,n,num,1)* atoms%rmsh(:atoms%jri(n),n)*atoms%rmsh(:atoms%jri(n),n)
             CALL intgr3(rh,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),rhs)
             IF (lh<4) pdis(lh,n,num)=rhs
             pdis(4,n,num)=pdis(4,n,num)+rhs
          ENDDO
          pdis(:,n,num) = pdis(:,n,num)*atoms%neq(n)
          dis(num)=dis(num)+SUM(pdis(:,n,num))
          pdis(:,n,num) = SQRT(pdis(:,n,num)/atoms%volmts(n))*1000.
          na = na + atoms%neq(n)
       ENDDO
       !     ----> interstitial part
       !         ---> create density in the real space
       i = 0
       DO  i3 = 0,2*stars%mx3
          k3 = i3
          IF (k3.GT.stars%mx3) k3 = k3 - 2*stars%mx3-1
          DO  i2 = 0,2*stars%mx2
             k2 = i2
             IF (k2.GT.stars%mx2) k2 = k2 - 2*stars%mx2-1
             DO  i1 = 0,2*stars%mx1
                k1 = i1
                IF (k1.GT.stars%mx1) k1 = k1 - 2*stars%mx3-1
                i = i + 1
                id3 = stars%ig(k1,k2,k3)
                phase = stars%rgphs(k1,k2,k3)
                IF (id3.EQ.0) THEN
                   af3(i,1) = 0.
                   bf3(i,1) = 0.
                ELSE
                   af3(i,1) =  REAL(rhpw(id3,num,1))*REAL(phase) &
                           - AIMAG(rhpw(id3,num,1))*AIMAG(phase)
                   bf3(i,1) = AIMAG(rhpw(id3,num,1))*REAL(phase) &
                           +  REAL(rhpw(id3,num,1))*AIMAG(phase)
                END IF
             ENDDO
          ENDDO
       ENDDO

       CALL cfft(af3(1,1),bf3(1,1),nt,stars%mx1*2+1,stars%mx1*2+1,1)
       CALL cfft(af3(1,1),bf3(1,1),nt,stars%mx2*2+1,nk12,1)
       CALL cfft(af3(1,1),bf3(1,1),nt,stars%mx3*2+1,nt,1)
       !         ---> form the dot product
       i = 0
       DO  i3 = 0,stars%mx3*2
          DO  i2 = 0,stars%mx2*2
             DO  i1 = 0,stars%mx3*2
                i = i + 1
                af3(i,1) = af3(i,1)*af3(i,1)
                bf3(i,1) = 0.
             ENDDO
          ENDDO
       ENDDO
       !         ---> back fft
       CALL cfft(af3(1,1),bf3(1,1),nt,stars%mx1*2+1,stars%mx1*2+1,-1)
       CALL cfft(af3(1,1),bf3(1,1),nt,stars%mx2*2+1,nk12,-1)
       CALL cfft(af3(1,1),bf3(1,1),nt,stars%mx3*2+1,nt,-1)
       sumis = 0.
       i = 0
       DO  i3 = 0,stars%mx3*2
          k3 = i3
          IF (k3.GT.stars%mx3) k3 = k3 - stars%mx3*2-1
          DO  i2 = 0,stars%mx2*2
             k2 = i2
             IF (k2.GT.stars%mx2) k2 = k2 - stars%mx2*2-1
             DO  i1 = 0,stars%mx1*2
                k1 = i1
                IF (k1.GT.stars%mx1) k1 = k1 - stars%mx1*2-1
                i = i + 1
                phase = stars%rgphs(k1,k2,k3)
                id3 = stars%ig(k1,k2,k3)
                IF (id3.NE.0) THEN
                   sumis = sumis + REAL(CMPLX(af3(i,1),bf3(i,1))* CONJG(stars%ustep(id3)))*phase*fact
                END IF
             ENDDO
          ENDDO
       ENDDO
       pdis(0,0,num)=SQRT(sumis/cell%volint)*1000.
       dis(num) = dis(num) + sumis
       IF (input%film) THEN
          !     ----> vacuum part
          DO  ivac = 1,vacuum%nvac
             DO  ip = 1,vacuum%nmzxy
                !         ---> create density in the real space
                i = 0
                DO  i2 = 0,stars%mx2*2
                   k2 = i2
                   IF (k2.GT.stars%mx2) k2 = k2 - stars%mx2*2-1
                   DO  i1 = 0,stars%mx1*2
                      k1 = i1
                      IF (k1.GT.stars%mx1) k1 = k1 - stars%mx1*2-1
                      i = i + 1
                      id3 = stars%ig(k1,k2,0)
                      IF (id3.EQ.0) THEN
                         af2(i,1) = 0.
                         bf2(i,1) = 0.
                      ELSE
                         id2 = stars%ig2(id3)
                         phase = stars%rgphs(k1,k2,0)
                         IF (id2.EQ.1) THEN
                            af2(i,1) = rhv0(ip,ivac,num,1)
                            bf2(i,1) = 0.
                         ELSE
                            af2(i,1) = REAL(rhv1(ip,id2-1,ivac,num,1))
                            bf2(i,1) =AIMAG(rhv1(ip,id2-1,ivac,num,1))
                         END IF
                      END IF
                   ENDDO
                ENDDO
                CALL cfft(af2(1,1),bf2(1,1),nk12,stars%mx1*2+1,stars%mx1*2+1,1)
                CALL cfft(af2(1,1),bf2(1,1),nk12,stars%mx2*2+1,nk12,1)
                !         ---> form dot product
                i = 0
                DO  i2 = 0,stars%mx2*2
                   DO  i1 = 0,stars%mx1*2
                      i = i + 1
                      af2(i,1) = af2(i,1)*af2(i,1)
                      bf2(i,1) = 0.
                   ENDDO
                ENDDO
                !         ---> back fft
                CALL cfft(af2(1,1),bf2(1,1),nk12,stars%mx1*2+1,stars%mx1*2+1,-1)
                CALL cfft(af2(1,1),bf2(1,1),nk12,stars%mx2*2+1,nk12,-1)
                disz(npz-ip) = cell%area*af2(1,1)/REAL(nk12)
             ENDDO
             !         ---> beyond warping region
             DO  ip = vacuum%nmzxy + 1,vacuum%nmz
                disz(npz-ip) = cell%area*rhv0(ip,ivac,num,1)
             ENDDO
             CALL intgz0(disz,vacuum%delz,vacuum%nmz,sumz,tail)
             IF (sym%zrfs .OR. sym%invs) facv = 2.0
             dis(num) = dis(num) + facv*sumz
          ENDDO
       END IF
    ENDDO

    dis(:) = SQRT(dis(:)/cell%vol)*1000.

    WRITE (6,FMT=8000) iter,dis(1)
8000 FORMAT (/,'----> distance of  the potential for it=',i3,':',f11.6, ' mhtr/bohr**3')
    WRITE(6,*) "Details of potential differences for each atom type"
    WRITE(6,*) "Atom: total difference: difference of first sphhar"
    DO n=1,atoms%ntype
       WRITE(6,"(i5,' : ',f10.6,' : ',4f10.6)") n,pdis(4,n,1),pdis(0:3,n,1)
    ENDDO
    WRITE(6,*) "Difference of interstitial:",pdis(0,0,1)
    IF (input%jspins.EQ.2) THEN
       WRITE (6,FMT=8010) iter,dis(2)
8010   FORMAT (/,'----> distance of spin potential for it=',i3,':', f11.6,' mhtr/bohr**3')
       WRITE(6,*) "Details of potential differences for each atom type"
       WRITE(6,*) "Atom: total difference: difference of first sphhar"
       DO n=1,atoms%ntype
          WRITE(6,"(i5,' : ',f10.6,' : ',4f10.6)") n,pdis(4,n,2),pdis(0:3,n,2)
       ENDDO
       WRITE(6,*) "Difference of interstitial:",pdis(0,0,2)
    END IF
  CONTAINS
    SUBROUTINE priv_cdndif(m1,m2,m3,n1,n2,n3,s1,s2,film,rhsp,rhpw,rhv0,rhv1)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: m1,m2,m3,n1,n2,n3
      REAL,    INTENT (IN) :: s1,s2
      COMPLEX, INTENT (INOUT) :: rhpw(:,:,:)
      COMPLEX, INTENT (INOUT) :: rhv1(:,:,:,:,:)
      REAL,    INTENT (INOUT) :: rhsp(:,:,:,:,:)
      REAL,    INTENT (INOUT) :: rhv0(:,:,:,:)
      LOGICAL,INTENT (IN)     :: film

      rhsp(:,:,:,m1,n1)=s1*rhsp(:,:,:,m2,n2) +  s2*rhsp(:,:,:,m3,n3)
      rhpw(:,m1,n1)    =s1*rhpw(:,m2,n2)     +  s2*rhpw(:,m3,n3)
      IF (film) THEN
         rhv0(:,:,m1,n1)  =s1*rhv0(:,:,m2,n2)   + s2*rhv0(:,:,m3,n3)
         rhv1(:,:,:,m1,n1)=s1*rhv1(:,:,:,m2,n2) + s2*rhv1(:,:,:,m3,n3)
      ENDIF
    END SUBROUTINE priv_cdndif

  END SUBROUTINE potdis
END MODULE m_potdis
