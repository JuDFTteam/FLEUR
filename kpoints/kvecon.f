      MODULE m_kvecon
      use m_juDFT
!
! This subroutine determines the k-points with which we
! will calculate the band-structure. The first ncorn
! points are the corners of the irreducible wedge, as
! determined in subroutine bzone. We then generate
! nkpt-ncorn points inside the wedge, distributed in
! such a way that the minimal distance between the points
! is maximal. This subroutine does not find the optimal
! distribution of k-points, it only finds a reasonable one.
!
      CONTAINS
      SUBROUTINE kvecon(
     >                  iofile,ibfile,mkpt,mface,
     >                  nkpt,ncorn,nsym,nface,rcmt,fdist,fnorm,cpoint,
     <                  kvc )

      IMPLICIT NONE

! ... Arguments ...
      INTEGER, INTENT (IN) :: mkpt,mface
      INTEGER, INTENT (IN) :: iofile,ibfile
      INTEGER, INTENT (IN) :: nkpt,ncorn,nsym,nface
      REAL,    INTENT (IN) :: fdist(mface),fnorm(3,mface)
      REAL,    INTENT (IN) :: rcmt(3,3),cpoint(3,mface)
      REAL,    INTENT (OUT) :: kvc(3,mkpt)

! ... Locals ...
      INTEGER nk,i,lmin,l,j,n1,n2,n3,i1,i2,i3,ncd
      REAL d,dm,dist,dmax,xmin,cn,alpha,thrd
      REAL knew(3),kc(3),cnorm(3),xnorm(3),cand(3,48*mkpt)
!
! ... Intrinsic Functions ...
      INTRINSIC abs,sqrt

      IF ( nkpt.LT.ncorn ) THEN
        WRITE (iofile,'(1x,''nkpt='',i4,'' ,ncorn='',i4)') nkpt,ncorn
         CALL juDFT_error("nkpt<ncorn ",calledby="kvecon")
      ENDIF
      IF ( nkpt.GT.mkpt ) THEN
        WRITE (iofile,'(1x,''nkpt='',i4,'' , mkpt='',i4)') nkpt, mkpt
         CALL juDFT_error("nkpt>mkpt ",calledby="kvecon")
      ENDIF
      thrd=1.00/3.00
!
! Construct list of candidate points
!
      d = 1.00
      DO  i = 1, 3
        dist = 0.00
        DO  j=1,3
          dist = dist + rcmt(j,i)**2
        ENDDO
        dist = sqrt(dist)
        xnorm(i) = dist
        d = d * dist
      ENDDO

      d = thrd*(d/(nkpt*nsym))**thrd
      n1 = int(xnorm(1)/d)+1
      n2 = int(xnorm(2)/d)+1
      n3 = int(xnorm(3)/d)+1
      WRITE( iofile,'('' n1 = '',i5,'' n2 = '',i5,'' n3 = '',i5)') 
     +                                                    n1,n2,n3
      ncd=0

      WRITE (iofile,'('' $$$$$ check $$$$$ '')')
      WRITE (iofile,'(''   ncorn = '',i4)') ncorn
      WRITE (ibfile,'(''    cpoints '')')
      DO i = 1, ncorn
        WRITE (ibfile,8901) cpoint(1,i),cpoint(2,i),cpoint(3,i)
      ENDDO
 8901 FORMAT ('    ( ',2(f13.6,','),f13.6,' )',/)
      WRITE (ibfile,'(''   nface = '',i4)') nface
      DO i = 1, nface
        WRITE (ibfile,8902) fdist(i)
        WRITE (ibfile,8901) fnorm(1,i),fnorm(2,i),fnorm(3,i)
      ENDDO
 8902 FORMAT ('  fdist = ',f13.6)
C
      DO i1=-n1,n1
        DO i2=-n2,n2
          DO i3=-n3,n3
            DO i = 1, 3
              knew(i) = i1*rcmt(i,1)/n1 +
     +                  i2*rcmt(i,2)/n2 + 
     +                  i3*rcmt(i,3)/n3
            ENDDO
            DO l = 1, nface
              alpha=0.0e0
              DO i = 1, 3
                alpha = alpha + knew(i)*fnorm(i,l)
              ENDDO
              IF ( alpha.GT.( fdist(l)+1.0e-6*d ) ) GOTO 200
            ENDDO
            ncd = ncd + 1
            IF (ncd>48*mkpt)  CALL juDFT_error("ncd>ncmax",calledby
     +           ="kvecon")
            DO i = 1, 3
              cand(i,ncd)=knew(i)
            ENDDO 
  200       CONTINUE
          ENDDO
        ENDDO
      ENDDO
      IF (ncd<nkpt)  CALL juDFT_error("ncd<nkpt",calledby ="kvecon")

! Initialize the kpoints
      DO nk=1,ncorn
        DO i=1,3
          kvc(i,nk) = cpoint(i,nk)
        ENDDO
      ENDDO
      DO nk=ncorn+1,nkpt
        DO i=1,3
          kvc(i,nk) = cpoint(i,1)
        ENDDO
      ENDDO
      IF ( nkpt.eq.ncorn ) GOTO 2
!
! Enter dynamical cycle. We find a k-point with the shortest
! distance to the other k-points, take it out, and put it in
! a position where it has a maximal minimal distance to all
! the other kpoints
!
    1 continue
C
C FIND MINIMAL DISTANCE K-POINT , EXCLUDE CORNER POINTS
C
      lmin=0
      d=1.0e6
      do 2000 l=ncorn+1,nkpt
      dm=1.0e6
      do 2100 j=1,nkpt
      if ( j.eq.l ) goto 2100
      dist=0.0e0
      do 2200 i=1,3
      dist=dist+(kvc(i,j)-kvc(i,l))**2
 2200 continue
      if ( dist.lt.(dm*(1.0-1e-6))  ) dm=dist
 2100 continue
      if ( d.lt.(dm*(1.0+1e-6))  ) goto 2000
      d=dm
      lmin=l
 2000 continue
      if ( lmin.eq.0 )  CALL juDFT_error(" lmin=0 ",calledby="kvecon")
C
C Find the point which is furthest from all k-points. We
C restrict our points to be on a finite grid determined by
C the corner points.
C
      dmax=0.0e0
      write(iofile,'('' kpoint ('',i3,'')  ncd = '',i8,
     +               '' max. allowed is '',i8)') lmin,ncd,48*mkpt
      do 3000 n1=1,ncd
      do 3100 i=1,3
      knew(i)=cand(i,n1)
 3100 continue
      dm=1.0e6
      do 3200 j=1,nkpt
      if ( j.eq.lmin ) goto 3200
      dist=0.0e0
      do 3300 i=1,3
      dist=dist+(kvc(i,j)-knew(i))**2
 3300 continue
      if ( dist.lt.(dm*(1.0-1e-6))  ) dm=dist
 3200 continue
      if ( dm.lt.(dmax*(1.0+1e-6)) ) goto 3000
      dmax=dm
      do 3400 i=1,3
      kc(i)=knew(i)
 3400 continue
 3000 continue
C
C WE HAVE FOUND A NEW POINT KC . IF THE MINIMAL DISTANCE OF
C THIS POINT IS LESS THAN THE DISTANCE D WE ALREADY HAVE ,
C WE STOP THE DYNAMIC CYCLE .
C
      if ( dmax.lt.1.0001*d ) goto 4000
      do 3500 i=1,3
      kvc(i,lmin)=kc(i)
 3500 continue
      goto 1
 4000 continue
C
C WE HAVE FOUND A REASONABLE DISTRIBUTION . TO AVOID TETRAHEDRA
C WITH A FUNNY SHAPE , WE PUT ALL THE POINTS WHICH ARE NEAR A
C SIDE ON THIS SIDE .
C
      d=0.49*sqrt(d)
      do 4100 l=ncorn+1,nkpt
      xmin=1.0e6
      do n1=1,ncorn-2
       do n2=n1+1,ncorn-1
        do 4200 n3=n2+1,ncorn
        cnorm(1)=(cpoint(2,n1)-cpoint(2,n2))*(cpoint(3,n2)-cpoint(3,n3))
     &          -(cpoint(3,n1)-cpoint(3,n2))*(cpoint(2,n2)-cpoint(2,n3))
        cnorm(2)=(cpoint(3,n1)-cpoint(3,n2))*(cpoint(1,n2)-cpoint(1,n3))
     &          -(cpoint(1,n1)-cpoint(1,n2))*(cpoint(3,n2)-cpoint(3,n3))
        cnorm(3)=(cpoint(1,n1)-cpoint(1,n2))*(cpoint(2,n2)-cpoint(2,n3))
     &          -(cpoint(2,n1)-cpoint(2,n2))*(cpoint(1,n2)-cpoint(1,n3))
        cn=0.0e0
        do i=1,3
          cn=cn+cnorm(i)**2
        end do
        cn=1.0e0/sqrt(cn)
        do i=1,3
          cnorm(i)=cnorm(i)*cn
        end do
        alpha=0.0e0
        do i=1,3
          alpha=alpha+cnorm(i)*(kvc(i,l)-cpoint(i,n1))
        end do
        if ( abs(alpha).gt.abs(xmin) ) goto 4200
        xmin=alpha
        do i=1,3
          xnorm(i)=cnorm(i)
        end do
 4200   continue
       end do
      end do
      if ( abs(xmin).gt.d ) goto 4100
      do 4700 i=1,3
      kvc(i,l)=kvc(i,l)-xmin*xnorm(i)
 4700 continue
 4100 continue
C
C IF A POINT IS TOO CLOSE TO AN EDGE , WE PUT IT ON THAT EDGE .
C
      do 5000 l=ncorn+1,nkpt
      xmin=1.0e6
      do n1=1,ncorn-1
        do 5100 n2=n1+1,ncorn
          cnorm(1)=(cpoint(2,n1)-cpoint(2,n2))*(cpoint(3,n2)-kvc(3,l))
     &            -(cpoint(3,n1)-cpoint(3,n2))*(cpoint(2,n2)-kvc(2,l))
          cnorm(2)=(cpoint(3,n1)-cpoint(3,n2))*(cpoint(1,n2)-kvc(1,l))
     &            -(cpoint(1,n1)-cpoint(1,n2))*(cpoint(3,n2)-kvc(3,l))
          cnorm(3)=(cpoint(1,n1)-cpoint(1,n2))*(cpoint(2,n2)-kvc(2,l))
     &            -(cpoint(2,n1)-cpoint(2,n2))*(cpoint(1,n2)-kvc(1,l))
          cn=0.0e0
          do i=1,3
            cn=cn+cnorm(i)**2
          end do
          dist=0.0e0
          do i=1,3
            dist=dist+(cpoint(i,n1)-cpoint(i,n2))**2
          end do
          cn=sqrt(cn/dist)
          if ( cn.gt.xmin ) goto 5100
          xmin=cn
          cn=0.0e0
          do i=1,3
            cn=cn+(cpoint(i,n2)-cpoint(i,n1))*(kvc(i,l)-cpoint(i,n1))
          end do
          cn=cn/dist
          do i=1,3
            knew(i)=cpoint(i,n1)+cn*(cpoint(i,n2)-cpoint(i,n1))
          end do
 5100   continue
      end do
      if ( xmin.gt.d ) goto 5000
      do 5600 i=1,3
      kvc(i,l)=knew(i)
 5600 continue
 5000 continue
    2 CONTINUE
C
C OUTPUT OF KVECTORS
C
      write(iofile,6000) 2.0e0*d,nkpt
 6000 format(//,'   minimal distance  ',f10.5,
     & //,'    coordinates of',i5,'  kpoints in cart. units ',//)
      write(iofile,6100) ((kvc(i,l),i=1,3),l=1,nkpt)
 6100 format(2(3f13.6,3x),3f13.6)
      RETURN
      END SUBROUTINE kvecon
      END MODULE m_kvecon
