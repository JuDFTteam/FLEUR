      MODULE m_brzone
      use m_juDFT
!
! This subroutine finds the corner-points, the edges, and the
! faces of the irreducible wedge of the brillouin zone (IBZ).
!
      CONTAINS
      SUBROUTINE brzone(
     >                   rcmt,nsym,idrot,mface,nbsz,nv48,
     =                   cpoint,
     <                   xvec,ncorn,nedge,nface,fnorm,fdist)

      USE m_constants, ONLY : pimach
      IMPLICIT NONE

      INTEGER, PARAMETER :: ibfile = 42

      INTEGER, INTENT (IN) :: mface,nbsz,nv48
      INTEGER, INTENT (IN) :: nsym               ! number of symmetry elements
      REAL,    INTENT (IN) :: rcmt(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
      REAL,    INTENT (IN) :: idrot(3,3,48)      ! rotation matrices in cartesian repr.

      INTEGER, INTENT (OUT) :: ncorn,nedge,nface ! number of corners, faces and edges of the IBZ
      REAL,    INTENT (OUT) :: fnorm(3,mface)    ! normal vector of the planes bordering the IBZ
      REAL,    INTENT (OUT) :: fdist(mface)      ! distance vector of the planes bordering the IBZ
      REAL,    INTENT (OUT) :: cpoint(3,mface)   ! cartesian coordinates of corner points of IBZ
      REAL,    INTENT (OUT) :: xvec(3)           ! arbitrary vector lying in the IBZ
C
C   LOCAL variables
C
      REAL pi
      REAL scale,sum,amin,alpha
      REAL bmin,beta,cmin,cmax,gamma
      REAL sx,xmin
      INTEGER ntl,krecip(3),ntot,ip,i,j,l,n,m,nfp
      INTEGER nmin,mmin,lmin,lmax,nf,ncf
      INTEGER n1,n2,n3,nn,k,ii
C
C     Local working arrays and pointers
C
      REAL epoint(3,mface),fpoint(3,mface),cstart(3,2,mface)
      REAL fvec(3),evec(3),dir(3),c0(3),c1(3),c2(3),csum(3)
      REAL sk(3),yvec(3),ddist(nv48),dvec(3,nv48)
      INTEGER nplane(mface)
C
C----->  Intrinsic Functions
C
      INTRINSIC min,sqrt
C
      OPEN (ibfile,form='formatted',status='scratch')
c
      WRITE (ibfile,'('' brzone '')')
      ntot = (2*nbsz + 1)**3
      WRITE (ibfile,'('' ntot = '',i4,'' nsym = '',i4)') ntot,nsym
      ntl = ntot + nsym - 2
      WRITE (ibfile,'('' ntl = '',i4)') ntl
      WRITE (ibfile,'('' rcmt '',/)')
c     WRITE (ibfile,*) rcmt
      WRITE (ibfile,101) ((rcmt(i,j),j=1,3),i=1,3)
 101  FORMAT(/5x,3(f10.6,3x),2(/5x,3(f10.6,3x)))
C
C construct all boundary-planes
C first the planes that determine the first brillouin zone
C that is, the planes bisecting the line connecting the
C origin with a reciprocal lattice vector ( <> 0 )
C
      pi = pimach()
      DO i = 1,3
       sk(i) = 0.0
        DO j = 1,3
         sk(i)=sk(i)+rcmt(j,i)*rcmt(j,i)
        ENDDO
      ENDDO
      WRITE (ibfile,'('' sk(1) = j=1,3 of rcmt(j,1)**2 '')')
      WRITE (ibfile,97) (sk (ii),ii=1,3)
 97   FORMAT (/5x,'  sk(i) ',3(f13.6,2x))

      scale = sqrt(min(sk(1),sk(2),sk(3)))*0.1
      xvec(1) = scale
      xvec(2) = scale/sqrt(pi)
      xvec(3) = scale/pi
      WRITE (ibfile,98) (xvec(ii),ii=1,3)
 98   FORMAT (/5x,' xvec(i) ',3(f13.6,2x))

      n = 0
      DO n1 = -nbsz,nbsz
        krecip(1) = n1
        DO n2 = -nbsz,nbsz
        krecip(2) = n2
          DO n3 = -nbsz,nbsz
            IF ( .NOT.(n1.EQ.0.AND.n2.EQ.0.AND.n3.EQ.0) ) THEN
            krecip(3) = n3

            n = n + 1
            DO i = 1,3
              dvec(i,n) = 0.0
              DO j = 1,3
                dvec(i,n) = dvec(i,n) + rcmt(i,j)*krecip(j)
              ENDDO
            ENDDO
            WRITE (ibfile,99) n,(dvec(k,n),k=1,3)
 99         FORMAT(/5x,'  dvec(k,',i4,') ',3(f13.6,2x))

            sum = 0.0
            DO i = 1,3
              sum = sum + dvec(i,n)**2
              WRITE (ibfile,'('' sum = dvec**2 = '',f13.6)') sum
            ENDDO
            sum = sqrt(sum)
            ddist(n) = 0.5*sum
            WRITE (ibfile,'(/'' ddist('',i3,'')=(.5*sum**.5) '',f13.6)')
     >                                                        n,ddist(n)
            sum = 1.0/sum
            WRITE (ibfile,'(/'' sum = ( 1/(.5*sum**.5) )'',f13.6)') sum
            DO i = 1,3
              dvec(i,n) = dvec(i,n)*sum
            ENDDO
            WRITE (ibfile,'('' dvec(i,n) * latest sum '')')
            WRITE (ibfile,99) n,(dvec(k,n),k=1,3)

            ENDIF
          ENDDO
        ENDDO
      ENDDO

C
C construct the planes that determine the irreducible wedge
C that is, the planes bisecting the line connecting xvec
C with an element of the star of xvec ( <> xvec )
C
      WRITE (ibfile,'('' working on star of xvec '')')
      WRITE (ibfile,'('' ntot = '',i4,'' ntl = '',i4,/)') ntot,ntl
      DO n = ntot,ntl
        ddist(n) = 0.0
        WRITE (ibfile,'(/)')

        DO i = 1,3
          dvec(i,n)=-xvec(i)
          WRITE (ibfile,'('' dvec('',i3,i4,'')=(here-xvec('',i2,'') '',
     +                                         f13.6)') i,n,i,dvec(i,n)
          DO j=1,3
            dvec(i,n) = dvec(i,n) + idrot(i,j,n+2-ntot)*xvec(j)
            WRITE (ibfile,'('' idrot('',i3,i3,i4,'') = '',f10.6)') 
     +                            i,j,n+2-ntot,idrot(i,j,n+2-ntot)
            WRITE (ibfile,'('' xvec('',i3,'') = '',f6.4)') j,xvec(j)
            WRITE (ibfile,'('' dvec('',i3,i4,'') = '',f13.6,/)') 
     +                                             i,n,dvec(i,n)
          ENDDO
        ENDDO

        sum = 0.0
        DO i = 1,3
          sum =sum + dvec(i,n)**2
          WRITE (ibfile,'('' sum = dvec**2 = '',f13.6)') sum
        ENDDO
        sum = 1.0/sqrt(sum)
        WRITE (ibfile,'(/'' sum = ( 1/(sum**.5) )'',f13.6)') sum
        DO i = 1,3
            dvec(i,n)=dvec(i,n)*sum
        ENDDO
        WRITE (ibfile,'('' dvec(i,n) * latest sum '')')
        WRITE (ibfile,99) n,(dvec(k,n),k=1,3)
      ENDDO
      nn = ntl - ntot + 1
C
C find the point on the line determined by the origin and xvec
C which is on the nearest boundary plane
C
      WRITE (ibfile,'(/,'' find points on nearest boundary plane '')')
      amin = scale*99999.9
      nmin = 0
      DO n = 1,ntl
        sum = 0.0
        DO i = 1,3
          sum = sum + xvec(i)*dvec(i,n)
        ENDDO
        WRITE (ibfile,'('' sum('',i4,'') = '',f13.6)') n,sum
        IF ( abs(sum).GT.1.0e-10 ) THEN
          alpha=ddist(n)/sum
          WRITE (ibfile,'('' alpha('',i4,'') = '',f13.6)') n,alpha
          IF ( .NOT.((alpha.LE.0.0).OR.(alpha.GT.amin)) ) THEN 
            amin = alpha
            nmin = n
            WRITE (ibfile,'('' nmin = '',i4)') n
          ENDIF
        ENDIF
      ENDDO
      IF ( nmin==0 )  CALL juDFT_error("bzone1",calledby ="brzone")
      WRITE (ibfile,'('' amin = '',f13.6)') amin
      DO i = 1,3
        fvec(i) = amin*xvec(i)
      ENDDO
      WRITE (ibfile,'('' fvec('',i3,'') = '',f13.6)') (i,fvec(i),i=1,3)
      nplane(1) = nmin
C
C find the nearest edge in this plane, along the line connecting
C fvec and the center of the plane, given by dvec*ddist
C
      WRITE (ibfile,'(/,'' find nearest edge in this plane '')')
      bmin = scale*99999.9
      mmin = 0
      DO m = 1,ntl
        IF ( m.NE.nmin ) THEN
          sum=0.0
          DO i = 1,3
            sum = sum + dvec(i,m)*(fvec(i)-dvec(i,nmin)*ddist(nmin))
          ENDDO
          WRITE (ibfile,'('' sum('',i4,'') = '',f13.6)') m,sum
          IF ( abs(sum).GT.1.0e-10 ) THEN
            beta = ddist(m)
            WRITE (ibfile,'('' beta('',i4,'') = '',f13.6)') m,ddist(m)
            DO i = 1,3
              beta=beta-fvec(i)*dvec(i,m)
            ENDDO
            WRITE (ibfile,'('' beta-fvec(i)*dvec(i,m) = '',f13.6)') beta
            beta = beta/sum
            IF ( .NOT.((beta.LT.0.0).OR.(beta.GT.bmin)) ) THEN
              bmin=beta
              mmin=m
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ( mmin==0 )  CALL juDFT_error("bzone2",calledby ="brzone")
      DO i = 1,3
        evec(i) = fvec(i) + bmin*(fvec(i)-dvec(i,nmin)*ddist(nmin))
      ENDDO
      WRITE (ibfile,'(/,'' evec('',i3,'') = '',f13.6)') 
     +                                (i,evec(i),i=1,3)
C
C find innermost boundary plane for this edge
C
      WRITE (ibfile,'(/,'' find innermost boundary plane '')')
      xmin = scale*99999.9
      mmin = 0
      DO  m = 1,ntl
        IF ( m.NE.nmin ) THEN
          sum = ddist(m)
          sx  = 0.0
          DO i = 1,3
            sum = sum - dvec(i,m)*evec(i)
            sx  = sx  + dvec(i,m)*dvec(i,nmin)
          ENDDO
          IF ( .NOT.((abs(sum).GT.1.e-10).OR.(sx.GT.xmin)) ) THEN
            xmin = sx
            mmin = m
          ENDIF
        ENDIF
      ENDDO
      IF ( mmin.EQ.0 )  CALL juDFT_error("bzone25",calledby="brzone")
      nplane(2) = mmin
C
C find direction of the edge
C
      dir(1) = dvec(2,nmin)*dvec(3,mmin) - dvec(3,nmin)*dvec(2,mmin)
      dir(2) = dvec(3,nmin)*dvec(1,mmin) - dvec(1,nmin)*dvec(3,mmin)
      dir(3) = dvec(1,nmin)*dvec(2,mmin) - dvec(2,nmin)*dvec(1,mmin)
      WRITE (ibfile,'('' dir('',i3,'') = '',f13.6)') (i,dir(i),i=1,3)
C
C find the corner points on this edge
C
      WRITE (ibfile,'('' find corner points on this edge '')')
      cmin = scale*99999.9
      cmax = -cmin
      lmin = 0
      lmax = 0
      DO l=1,ntl
        IF ( (l.EQ.nmin).OR.(l.EQ.mmin) ) GOTO 2700
        sum = 0.0
        DO i=1,3
          sum = sum + dir(i)*dvec(i,l)
        ENDDO
        IF ( abs(sum).LT.1.0e-10 ) GOTO 2700
        gamma=ddist(l)
        DO i = 1,3
          gamma = gamma - evec(i)*dvec(i,l)
        ENDDO
        gamma = gamma/sum
        IF ( gamma.GE.0.0 ) THEN
          IF (gamma.gt.cmin) GOTO 2700
          cmin=gamma
          lmin=l
          GOTO 2700
        ENDIF
        IF ( gamma.GE.cmax ) THEN
          cmax=gamma
          lmax=l
        ENDIF
 2700   CONTINUE
      ENDDO
      IF ( lmax*lmin.EQ.0 ) CALL juDFT_error("bzone3",calledby="brzone")
      WRITE (ibfile,'('' cmax = '',f13.6)') cmax
      WRITE (ibfile,'('' cmin = '',f13.6)') cmin
      DO i=1,3
         c0(i) = evec(i)+cmax*dir(i)
         WRITE (ibfile,'('' dir('',i3,'') = '',f13.6)') i,dir(i)
         WRITE (ibfile,'('' evec('',i3,'') = '',f13.6)') i,evec(i)
         WRITE (ibfile,'('' c0('',i3,'') = '',f13.6)') i,c0(i)
         c1(i)=evec(i)+cmin*dir(i)
         WRITE (ibfile,'('' c1('',i3,'') = '',f13.6)') i,c1(i)
      ENDDO
C
C prepare the list of corner points, etc, for the
C general scheme of finding the boundaries of the
C irreducible wedge of the first brillouin zone
C
      WRITE (ibfile,'(/,'' prepare list of corner points '')')
      DO i = 1,3
        cstart(i,1,1) = c0(i)
        cstart(i,2,1) = c1(i)
        cstart(i,1,2) = c1(i)
        cstart(i,2,2) = c0(i)
        cpoint(i,1)   = c0(i)
        cpoint(i,2)   = c1(i)
        epoint(i,1)   = 0.5*(c0(i)+c1(i))
        WRITE (ibfile,'('' cstart('',i2,'',1,1) = '',f13.6)') i,c0(i)
        WRITE (ibfile,'('' cstart('',i2,'',2,1) = '',f13.6)') i,c1(i)
        WRITE (ibfile,'('' cstart('',i2,'',1,2) = '',f13.6)') i,c1(i)
        WRITE (ibfile,'('' cstart('',i2,'',2,2) = '',f13.6)') i,c0(i)
        WRITE (ibfile,'('' cpoint('',i2,'',1) = '',f13.6)') i,c0(i)
        WRITE (ibfile,'('' cpoint('',i2,'',2) = '',f13.6)') i,c1(i)
        WRITE (ibfile,'('' epoint('',i2,'',1) = '',f13.6)')i,epoint(i,1)
      ENDDO
      ncorn = 2
      nedge = 1
      nface = 2
      nf = 0
C
C enter general loop which determines all corners and all edges
C of all faces , new faces are added to the list nplane
C
 4000 CONTINUE
      nf  = nf + 1
      nfp = nplane(nf)
C
C we consider face number nf
C start with the corner points of cstart , notice that the order
C of the corner points is important and is determined by the
C order in the outer product of the vectors dvec
C
      DO i=1,3
        c0(i)   = cstart(i,1,nf)
        c1(i)   = cstart(i,1,nf)
        c2(i)   = cstart(i,2,nf)
        csum(i) = cstart(i,1,nf) + cstart(i,2,nf)
      ENDDO
      ncf = 2
 4200 CONTINUE
C
C determine the point fvec
C
      fvec(1) = dvec(2,nfp)*(c2(3)-c1(3))-dvec(3,nfp)*(c2(2)-c1(2))
      fvec(2) = dvec(3,nfp)*(c2(1)-c1(1))-dvec(1,nfp)*(c2(3)-c1(3))
      fvec(3) = dvec(1,nfp)*(c2(2)-c1(2))-dvec(2,nfp)*(c2(1)-c1(1))
c     WRITE (ibfile,'('' pt fvec('',i3,'') = '',f13.6)') (i,fvec(i),i=1,3)
      DO i = 1,3
        fvec(i) = 0.5*(c2(i)+c1(i)) + 0.001*fvec(i)
      ENDDO
C
C determine the edge connected to c2 by moving outwards on c2-c1
C and finding the nearest intersection with a boundary plane
C on the line connecting this point and fvec , which is
C on the correct side of the line c2-c1 , by construction ,
C because of the way we order the corner points
C
      DO i = 1,3
        yvec(i) = c2(i) + 1.0e-5*(c2(i)-c1(i))
      ENDDO
C
C find nearest boundary plane
C
      bmin = scale*99999.9
      mmin = 0
      DO m = 1,ntl
       IF ( m.NE.nplane(nf) ) THEN
         sum = 0.0
         DO i = 1,3
           sum = sum + dvec(i,m)*(yvec(i)-fvec(i))
         ENDDO
         IF ( abs(sum).GE.1.0e-10 ) THEN
           beta=ddist(m)
           DO i = 1,3
             beta = beta - fvec(i)*dvec(i,m)
           ENDDO
           beta = beta/sum
           IF ( .NOT.((beta.LT.0.0).OR.(beta.GT.bmin)) ) THEN
             bmin = beta
             mmin = m
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ( mmin.EQ.0 )  CALL juDFT_error("bzone4",calledby="brzone")
C
C construct direction of this edge
C
      dir(1) = dvec(2,nfp)*dvec(3,mmin) - dvec(3,nfp)*dvec(2,mmin)
      dir(2) = dvec(3,nfp)*dvec(1,mmin) - dvec(1,nfp)*dvec(3,mmin)
      dir(3) = dvec(1,nfp)*dvec(2,mmin) - dvec(2,nfp)*dvec(1,mmin)
      WRITE (ibfile,'(''2 dir('',i3,'') = '',f13.6)') (i,dir(i),i=1,3)
C
C find other corner point on this edge
C
      cmin = scale*99999.9
      lmin = 0
      DO l = 1,ntl
        IF ( .NOT.((l.EQ.nplane(nf)).OR.(l.EQ.mmin)) ) THEN 
          sum = 0.0
          DO i = 1,3
            sum = sum + dir(i)*dvec(i,l)
          ENDDO
          IF ( abs(sum).GE.1.0e-10 )  THEN
            gamma=ddist(l)
            DO i = 1,3
              gamma = gamma - dvec(i,l)*c2(i)
            ENDDO
            gamma = gamma/sum
            IF ( .NOT.((gamma.LT.1.0e-9).OR.(gamma.GT.cmin)) ) THEN
              cmin = gamma
              lmin = l
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF ( lmin.EQ.0 )  CALL juDFT_error("bzone5",calledby="brzone")
C
C move c2 and c1
C
      DO i = 1,3
        c1(i)   = c2(i)
        c2(i)   = c1(i) + cmin*dir(i)
        evec(i) = 0.5*( c1(i)+c2(i) )
      ENDDO
      WRITE (ibfile,'(''corner c1('',i3,'')='',f13.6)') (i,c1(i),i=1,3)
      WRITE (ibfile,'(''corner c2('',i3,'')='',f13.6)') (i,c2(i),i=1,3)
      WRITE (ibfile,'(''evec('',i3,'') = '',f13.6)') (i,evec(i),i=1,3)
C
C find innermost boundary plane for this edge
C
      xmin = scale*99999.9
      mmin = 0
      WRITE (ibfile,'(/,''bzone55 loop ntl='',i4,'' nfp='',i4)') ntl,nfp
      DO m = 1,ntl
        IF ( m.NE.nfp ) THEN
          sum = ddist(m)
          sx  = 0.0
          DO i=1,3
           sum = sum - dvec(i,m)*evec(i)
           sx  = sx  + dvec(i,m)*dvec(i,nfp)
          ENDDO
          IF ( .NOT.((abs(sum).GT.1.0e-6).OR.(sx.GT.xmin)) ) THEN
            xmin = sx
            mmin = m
            WRITE (ibfile,'('' m = '',i4,'' xmin = '',f16.12,'' nfp = ''
     +                                                ,i4)')  m,xmin,nfp
          ENDIF
        ENDIF
      ENDDO
      WRITE (ibfile,'('' m = '',i4,'' xmin = '',f16.12,'' nfp = '',i4)')
     +                                                        m,xmin,nfp
      IF ( mmin.EQ.0 )  CALL juDFT_error("bzone55",calledby="brzone")
C
C check if we have found a new face or not
C
      DO ip = 1,nface
         IF (nplane(ip).EQ.mmin) GOTO 5400
      ENDDO
      nface = nface + 1
      WRITE (ibfile,'('' nface = '',i4)') nface
      nplane(nface) = mmin
      DO i = 1,3
        cstart(i,1,nface) = c2(i)
        cstart(i,2,nface) = c1(i)
        WRITE (ibfile,'('' cstart('',i3,'', 1,'',i3,'') = '',f13.6)')
     +   i,nface,c2(i)
        WRITE (ibfile,'('' cstart('',i3,'', 2,'',i3,'') = '',f13.6)')
     +   i,nface,c1(i)
      ENDDO 
 5400 CONTINUE
C
C check if the new corner and edge points are contained
C in the list of existing points
C
      DO ip = 1,ncorn
        sum = 0.00
        DO i = 1,3
          sum = sum + (c2(i) - cpoint(i,ip))**2
        ENDDO
        IF ( abs(sum).LT.1.0e-10 ) GOTO 6300
      ENDDO 
      ncorn = ncorn + 1
      WRITE (ibfile,'('' ncorn = '',i5)') ncorn
      DO i = 1,3
        cpoint(i,ncorn) = c2(i)
      ENDDO
 6300 CONTINUE
c
      DO ip = 1,nedge
        sum = 0.0
        DO i = 1,3
          sum = sum + (evec(i) - epoint(i,ip))**2
        ENDDO 
      IF ( abs(sum).LT.1.0e-10 ) GOTO 6700
      ENDDO
      nedge = nedge + 1
      WRITE (ibfile,'('' nedge = '',i5)') nedge
      DO i = 1,3
       epoint(i,nedge) = evec(i)
      ENDDO
 6700 CONTINUE
C
C check if we have all points on this face
C
      sum = 0.0
      DO i = 1,3
       sum = sum + ( c2(i) - c0(i) )**2
      ENDDO
      IF ( abs(sum).GT.1.0e-10 ) THEN
        ncf = ncf + 1
        WRITE (ibfile,'('' nface = '',i4)') nface
        GOTO 4200
      ENDIF
C
C we have found all corner points on this face
C determine the center of gravity of this face
C
      DO i = 1,3
        fpoint(i,nf) = csum(i)/ncf
      ENDDO
      IF ( nf.LT.nface ) GOTO 4000
c
      DO ip  =1,nface
        nf = nplane(ip)
        fdist(ip) = ddist(nf)
        DO i=1,3
           fnorm(i,ip) = dvec(i,nf)
        ENDDO
      ENDDO
c

!      WRITE(*,*) 'ncorn', ncorn
!      WRITE(*,*) 'nedge', nedge
!      WRITE(*,*) 'nface', nface
!      WRITE(*,*) 'faces:'
!      DO ip  =1,nface
!         WRITE(*,'(4f20.13)') fnorm(:,ip), fdist(ip)
!      END DO
!      WRITE(*,*) 'coners:'
!      DO ip = 1,ncorn
!         WRITE(*,'(3f20.13)') cpoint(:,ip)
!      END DO

      WRITE (6,7100) ncorn,nedge,nface
      WRITE (ibfile,7100) ncorn,nedge,nface
 7100 FORMAT (///,'  the irreducible wedge of the first brillouin'
     $,' zone has :  ',/,
     $     i10,'     corner points   ',/,
     $     i10,'     edges           ',/,
     $     i10,'     faces           ')
      IF ( (ncorn + nface - nedge)/=2 )  CALL juDFT_error("bzone6"
     +     ,calledby ="brzone")
      WRITE (6,7200) ((cpoint(i,ip),i=1,3),ip=1,ncorn)
      WRITE (ibfile,7200) ((cpoint(i,ip),i=1,3),ip=1,ncorn)
 7200 FORMAT(//,'    corner points in carthesian units ',
     $     99(/,3f10.5))

      CLOSE (ibfile)
      RETURN
      END SUBROUTINE brzone
      END MODULE m_brzone
