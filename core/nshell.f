      MODULE m_nshell
      use m_juDFT
c-------------------------------------------------------------------
c     Constructs the neighbouring shells up to the given number nsh
c                                                 M. Lezaic '04
c-------------------------------------------------------------------
      CONTAINS
      SUBROUTINE nshell(
     >                  amat,t,nsh,dims,nmax,shmax,film,zcoord,
     <                  nat,R,lenR,nop,mrot,deltaz)
      IMPLICIT NONE
c     ..
c     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: nsh,dims,nmax,shmax,nop
      LOGICAL, INTENT (IN)  :: film
      REAL,    INTENT (IN)   :: zcoord,deltaz
c     ..
c     .. Array Arguments ..
      REAL,    INTENT (IN)   :: amat(3,3),t(3)
      INTEGER, INTENT (IN)   :: mrot(3,3,nop)
      INTEGER, INTENT (OUT)  :: nat(dims)
      REAL,    INTENT (OUT)  :: R(3,shmax,dims),lenR(dims)

c     
c     .. Local Scalars ..
      INTEGER n,i,c,n1,n2,n3,fill,fill1,xmax,ymax,zmax
      REAL lentn,t3
      REAL, PARAMETER:: tol=0.0000001
c
c     .. Local Arrays
      REAL tnC(3),tn(3),Raux(3,nop),Raux1(3,shmax),Rrot(3)
c    
c     .. Intrinsic Functions ..
      INTRINSIC SQRT,ABS,REAL,MIN
c------------------------------------------------------------------
      c=0
      nat(:)=0
      R(:,:,:)=0.
      lenR(:)=0.

      IF (film) THEN  !Added for Film-Jij calculations 10/23/06 B. Hardrat
        zmax = 0
        t3 = deltaz    !Added for films with more than one monolayer 07/10  S.Schroeder
      ELSE             !format of coordinates in file shells in case of film:
       zmax=nmax       !x (relative coord.),y (relative coord.), z (a.u.)
       t3 = t(3)
      END IF

      xmax=nmax
      ymax=nmax

      fst: DO n3=-zmax,zmax
        snd: DO n2=-ymax,ymax
          trd: DO n1=-xmax,xmax
            tn(1)=t(1)-REAL(n1)
            tn(2)=t(2)-REAL(n2)
            tn(3)=t3  -REAL(n3)
            IF ( (ABS(tn(1)).LT.tol).AND.
     &           (ABS(tn(2)).LT.tol).AND.
     &           (ABS(tn(3)).LT.tol) ) CYCLE trd

            tnC(:)=tn(1)*amat(:,1)+tn(2)*amat(:,2)+tn(3)*amat(:,3)
            lentn = SQRT( tnC(1)**2+tnC(2)**2+tnC(3)**2 )
            DO i = 1, c
              IF (ABS(lentn-lenR(i)).LT.tol) THEN
                DO n = 1, nat(i)
                  IF( (ABS(tn(1)-R(1,n,i)).LT.tol).AND.
     &                (ABS(tn(2)-R(2,n,i)).LT.tol).AND.
     &                (ABS(tn(3)-R(3,n,i)).LT.tol) ) CYCLE trd
                ENDDO
                nat(i) = nat(i) + 1
                R(:,nat(i),i) = tn(:)
                  IF (film) THEN
                   R(3,nat(i),i) = zcoord
                  ENDIF
                CYCLE trd
              ENDIF
            ENDDO
            c=c+1
            DO i = 1, c-1
              IF (ABS(min(lentn,lenR(i))-lentn).LT.tol) THEN
                DO n=1,c-i
                  nat(c-n+1) = nat(c-n)
                  lenR(c-n+1) = lenR(c-n)
                  R(:,:,c-n+1) =R (:,:,c-n)
                ENDDO
                nat(i) = 1
                lenR(i) = lentn
                R(:,1,i) = tn(:)
                  IF (film) THEN
                     R(3,1,i) = zcoord
                  ENDIF

                CYCLE trd
              ENDIF
            ENDDO
            nat(c) = 1
            lenR(c) = lentn
            R(:,1,c) = tn(:)
              IF (film) THEN
                R(3,1,c) = zcoord
              ENDIF

          ENDDO trd
        ENDDO snd
      ENDDO fst


c-----------------------------------------------------
c ..
c .. Checking for inequivalent shells with the same lenR(i)
c ..
      dimsl: DO i = 1, dims
        IF (i.GT.c) EXIT dimsl
        Raux(:,1) = R(:,1,i)
        fill = 1

        nopl: DO n = 1, nop 
          Rrot(1)=Raux(1,1)*mrot(1,1,n) + Raux(2,1)*mrot(1,2,n) +
     +            Raux(3,1)*mrot(1,3,n)
          Rrot(2)=Raux(1,1)*mrot(2,1,n) + Raux(2,1)*mrot(2,2,n) +
     +            Raux(3,1)*mrot(2,3,n)
          Rrot(3)=Raux(1,1)*mrot(3,1,n) + Raux(2,1)*mrot(3,2,n) +
     +            Raux(3,1)*mrot(3,3,n)
         
          DO n1 = 1, fill
            IF((ABS(Rrot(1)-Raux(1,n1)).LT.tol).AND.
     &         (ABS(Rrot(2)-Raux(2,n1)).LT.tol).AND.
     &         (ABS(Rrot(3)-Raux(3,n1)).LT.tol)) CYCLE nopl
          ENDDO

          fill=fill+1              
          Raux(:,fill)=Rrot(:)
        
        ENDDO nopl

        IF (fill.LT.nat(i)) THEN
          fill1=0
          eqat: DO n = 1, nat(i)
            DO n1 = 1, fill
              IF(((ABS(R(1,n,i)-Raux(1,n1)).LT.tol).AND.
     &            (ABS(R(2,n,i)-Raux(2,n1)).LT.tol).AND.
     &            (ABS(R(3,n,i)-Raux(3,n1)).LT.tol)).OR.
     &           ((ABS(R(1,n,i)+Raux(1,n1)).LT.tol).AND.
     &            (ABS(R(2,n,i)+Raux(2,n1)).LT.tol).AND.
     &            (ABS(R(3,n,i)+Raux(3,n1)).LT.tol))) CYCLE eqat
            ENDDO

            fill1=fill1+1
            Raux1(:,fill1)=R(:,n,i)
          
          ENDDO eqat

          IF(fill1.GT.0) THEN
            c=c+1
            IF (c.GT.dims)  CALL juDFT_error("nshell:1")
            DO n = 1, c-i-1
              nat(c-n+1) = nat(c-n)
              lenR(c-n+1) = lenR(c-n)
              R(:,:,c-n+1) = R(:,:,c-n)
            ENDDO

            DO n = 1, fill1
              R(:,n,i+1) = Raux1(:,n)
            ENDDO
            lenR(i+1) = lenR(i)
            nat(i+1) = fill1

            DO n = 1, fill
              R(:,n,i) = Raux(:,n)
            ENDDO
            nat(i)=fill
          ENDIF
        ENDIF    ! fill < nat(i)
      ENDDO dimsl
c --------------------------------------------------
      DO i=1,c
        WRITE(117,5005) i,lenR(i),nat(i)
        DO n=1,nat(i)
          WRITE(117,5519) R(1,n,i),R(2,n,i),R(3,n,i)
        ENDDO
        WRITE(117,*) 
      ENDDO 
 5005 FORMAT(i4,1x,f14.10,1x,i4)
 5519 FORMAT(3(1x,f14.10))

      IF (nsh>c) THEN
         CALL juDFT_error("nsh greater than dimensioned",calledby
     +        ="nshell",hint ='increase nmax in jcoff2')
      ENDIF

      END SUBROUTINE nshell
      END MODULE m_nshell
