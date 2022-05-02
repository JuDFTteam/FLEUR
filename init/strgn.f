
c
c     *********************************************************
c     generate two- and three-dimensional stars
c     for slab geometry
c     e. wimmer   nov.1984    c.l.fu  1987
c     implementation of new box-dimension: to treat nonorthogonal
c     lattice systems
c     S. Bl"ugel, IFF, 17.Nov.97
c     *********************************************************
      CONTAINS
      SUBROUTINE strgn(
     >                 k1d,k2d,k3d,n3d,n2d,nop,nn2d,nn3d,natd,
     >                 nmzd,jspd,ntypsd,nmzxyd,ntypd,nlhd,jmtd,
     >                 nmz,jspins,ntypsy,nmzxy,ntype,nlh,jri,film,zrfs,
     >                 invs2,nvac,delz,rmt,dx,zatom,z1,invtab,
     >                 igrd,invs,nop2,bmat,gmax,symor,mrot,tau,
     <                 mx1,mx2,mx3,nq2,nq3,ngz,nk1,nk2,nk3,neq,
     <                 kv2,kv3,sk2,sk3,nstr,nstr2,star2d,ig,
     <                 rgphs,izmin,izmax,phi2,
     <                 kimax,igfft,pgfft,kimax2,igfft2,pgfft2,
     <                 pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy)
c
      USE m_constants, ONLY : pimach
      USE m_set, ONLY : dset
      USE m_spgrot
      USE m_angle
      USE m_types, ONLY : t_2dstar
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: k1d,k2d,k3d,n3d,n2d,nop,nn2d,nn3d,natd
      INTEGER, INTENT (IN) :: nmzd,jspd,ntypsd,nmzxyd,ntypd,nlhd,jmtd
      INTEGER, INTENT (IN) :: igrd,nop2,nmz,jspins,nmzxy
      INTEGER, INTENT (IN) :: nvac,ntype
      LOGICAL, INTENT (IN) :: zrfs,invs,film,invs2,symor
      REAL,    INTENT (IN) :: delz,z1
      REAL,    INTENT (INOUT) :: gmax
      INTEGER, INTENT (OUT):: kimax,kimax2,mx1,mx2,mx3,nq2,nq3
      INTEGER, INTENT (OUT):: ngz,nk1,nk2,nk3,izmin,izmax
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN)  :: jri(ntypd),nlh(ntypsd),invtab(nop)
      INTEGER, INTENT (IN)  :: ntypsy(natd),neq(ntypd),mrot(3,3,nop)
      REAL,    INTENT (IN)  :: rmt(ntypd),dx(ntypd),zatom(ntypd)
      REAL,    INTENT (IN)  :: bmat(3,3),tau(3,nop)
      INTEGER, INTENT (OUT) :: kv2(2,n2d),kv3(3,n3d)
      INTEGER, INTENT (OUT) :: nstr(n3d),nstr2(n2d)
      INTEGER, INTENT (OUT) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      INTEGER, INTENT (OUT) :: igfft(0:nn3d-1,2)
      INTEGER, INTENT (OUT) :: igfft2(0:nn2d-1,2)
      REAL,    INTENT (OUT) :: sk2(n2d),phi2(n2d),sk3(n3d)
      COMPLEX, INTENT (OUT) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      COMPLEX, INTENT (OUT) :: pgfft(0:nn3d-1)
      COMPLEX, INTENT (OUT) :: pgfft2(0:nn2d-1)
      COMPLEX, INTENT (OUT) :: pgft2x(0:) ! 0:nn2d-1 for igrd=1
      COMPLEX, INTENT (OUT) :: pgft2y(0:)
      COMPLEX, INTENT (OUT) :: pgft2xx(0:)
      COMPLEX, INTENT (OUT) :: pgft2yy(0:)
      COMPLEX, INTENT (OUT) :: pgft2xy(0:)
      TYPE (t_2dstar), INTENT (INOUT) :: star2d
!      contains  INTENT (OUT) ::
!      INTEGER ig2(n3d)
!      INTEGER igvac(n2d,-k3d:k3d)
!      INTEGER i2g(-k1d:k1d,-k2d:k2d)
!      REAL r2gphs(-k1d:k1d,-k2d:k2d)
C     ..
C     .. Local Scalars ..
      REAL arltv1,arltv2,arltv3,s,tpi
      REAL gmi,gla,eps,phi,gmax2
      REAL gfx,gfy,pon,pon2
      INTEGER j,k,k1,k2,k3,m0,mxx1,mxx2,n
      INTEGER ned1,nint,kdone,i
      LOGICAL to_ordered,ordered,to_unordered
      LOGICAL new,l_cdn1
      INTEGER kfx,kfy,kfz,kidx,nfftx,nffty,nfftz,kfft
      INTEGER nfftxy,norm,n1,kidx2,k2i
C     ..
C     .. Local Arrays ..
      REAL g(3),gsk3(n3d),phi3(n2d)
      COMPLEX phas(nop)
      INTEGER index(n3d),kr(3,nop),kv(3),kv3rev(n3d,3)
      INTEGER index2(n2d),index3(n3d),igz(n3d)
      CHARACTER*8 space(10),name(10)
c     ..
C     .. External Subroutines ..
      EXTERNAL boxdim,sort
c     ..
C     .. Intrinsic Functions ..
      INTRINSIC iabs,int,min0,sqrt
C     ..
      DATA space/10*'        '/
C
      tpi = 2 * pimach()
c
      nfftx = 3*k1d
      nffty = 3*k2d
      nfftz = 3*k3d
      nfftxy= 9*k1d*k2d
c
c..roa...look if stars are in ordered form in cdn1 file
c     ..
      DO i = 1,10
         name(i) = space(i)
      ENDDO
      ordered = .false.
      INQUIRE (file='cdn1',exist=l_cdn1)
      IF (l_cdn1) THEN
        OPEN (71,file='cdn1',form='unformatted',status='unknown')
        READ (71,end=99,err=99) name
        IF (name(10).EQ.'ordered*') ordered = .true.
  99    CONTINUE
        CLOSE (71)
      ELSE
        ordered = .true.
      ENDIF

      WRITE (6,'(a16,l1)') ' ordered stars: ',ordered
      INQUIRE (file='to_ordered',exist=to_ordered)
      INQUIRE (file='to_unordered',exist=to_unordered)
      IF (to_ordered.AND.ordered)
     >    STOP ' ..strgn1: to_ordered.and.ordered conflicts! '
      IF (to_unordered.AND.(.not.ordered))
     >    STOP ' ..strgn1: to_unordered.and..not.ordered  !!!'
c
c--->    read in information if exists
c
      OPEN (51,file='stars',form='unformatted',status='unknown')
      REWIND 51
c
c...roa..security..
      IF (to_unordered.OR.to_ordered)  GOTO 10
c...roa..security..
c
      READ (51,end=10,err=10) gmax,nq3,nq2,ngz,izmin,izmax,mx1,
     +  mx2,mx3,nk1,nk2,nk3
      IF (igrd.EQ.0) THEN
        READ (51,end=10,err=10) nstr,nstr2,rgphs,sk3,sk2,phi2,kv3,kv2,
     +            ig,star2d%igvac,star2d%ig2,star2d%i2g,star2d%r2gphs,
     +                         kimax,igfft,pgfft,kimax2,igfft2,pgfft2
      ELSE
        READ (51,end=10,err=10) nstr,nstr2,rgphs,sk3,sk2,phi2,kv3,kv2,
     +            ig,star2d%igvac,star2d%ig2,star2d%i2g,star2d%r2gphs,
     +                         kimax,igfft,pgfft,kimax2,igfft2,pgfft2,
     +                           pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy
      ENDIF

      igz = 0
      DO k = 1,nq3
        DO k3 = -k3d, k3d ! ngz
          DO k2 = 1,nq2
            IF (star2d%igvac(k2,k3) == k) THEN
              igz(k) = k3 + k3d
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      GOTO 270

   10 CONTINUE
C
C---> Determine Gmax box of size mx1, mx2, mx3,
c     for which |G(mx1,mx2,mx3)| < Gmax
c     arltv(i) length of reciprical lattice vector along direction (i)
C
c
c-gu   WRITE (16,*) gmax,bmat
!      CALL boxdim(
!     >            bmat,
!     <            arltv1,arltv2,arltv3)

      mx1 = k1d !int( gmax/arltv1 ) + 1
      mx2 = k2d !int( gmax/arltv2 ) + 1
      mx3 = k3d !int( gmax/arltv3 ) + 1
c
      WRITE (16,*) mx1,mx2,mx3
c
      mxx1 = 0
      mxx2 = 0
      nk1 = 2*mx1 + 1
      nk2 = 2*mx2 + 1
      nq2 = 0
      kv(3) = 0

      IF (film) THEN 

      k1_loop: DO k1 = mx1,-mx1,-1
        kv(1) = k1
        k2_loop: DO k2 = mx2,-mx2,-1
          kv(2) = k2
            DO j = 1,2
               g(j) = kv(1)*bmat(1,j) + kv(2)*bmat(2,j)
            ENDDO
            s = sqrt(g(1)**2+g(2)**2)
c--->   determine the angle of the G_{||} needed in odim calculations
c+odim YM cute little 'angle' is written by me to determine the phi2
            phi = angle(g(1),g(2))
c-odim
c
c--->   check if generated vector belongs to a star already
c--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel
c
            IF (s.LT.gmax) THEN      
               CALL spgrot(
     >                     nop,symor,tpi,mrot,tau,invtab,
     >                     kv,
     <                     kr,phas)
               DO n = 1,nop2
                  IF (mxx1.lt.kr(1,n)) mxx1 = kr(1,n)
                  IF (mxx2.lt.kr(2,n)) mxx2 = kr(2,n)
               ENDDO
               DO k = 1,nq2 
                  DO n = 1,nop2
                     IF (kr(1,n).EQ.kv2(1,k) .AND.
     +                   kr(2,n).EQ.kv2(2,k)) CYCLE k2_loop
                  ENDDO
               ENDDO
               nq2 = nq2 + 1            ! new representative found
               IF (nq2.GT.n2d) THEN
                 WRITE (6,8070) nq2,n2d
                 STOP 'nq2.GT.n2d'
               ENDIF
               DO j = 1,2
                  kv2(j,nq2) = kv(j)
               ENDDO
               sk2(nq2) = s
               phi2(nq2) = phi
            ENDIF
        ENDDO k2_loop
      ENDDO k1_loop
 8070 FORMAT ('nq2 = ',i5,' > n2d =',i5)

      IF( mxx1.GT.k1d) THEN
        WRITE (6,'(/'' mxx1.gt.k1d. mxx1,k1d='',2i4)')  mxx1,k1d
        STOP 'mxx1.gt.k1d'
      ENDIF

      IF (mxx2.gt.k2d) THEN
        WRITE (6,'(/'' mxx2.gt.k2d. mxx2,k2d='',2i4)')  mxx2,k2d
        STOP 'mxx2.gt.k2d'
      ENDIF

c--->    sort for increasing length sk2

      CALL sort(nq2,sk2,index)
      DO k = 1,nq2
         kv3rev(k,1) = kv2(1,index(k))
         kv3rev(k,2) = kv2(2,index(k))
         gsk3(k) = sk2(index(k))
         phi3(k) = phi2(index(k))
      ENDDO
      DO k = 1,nq2
         kv2(1,k) = kv3rev(k,1)
         kv2(2,k) = kv3rev(k,2)
         sk2(k) = gsk3(k)
         phi2(k) = phi3(k)
      ENDDO

c--roa   ....
c......sort stars of equal length 2D .....
      IF (ordered.OR.to_ordered) THEN
         i=1
         gla=0.
         gsk3(1)=0.0
         eps=1.e-10
         DO  k = 2,nq2
           IF (sk2(k)-gla.GE.eps) i=i+1
           gla=sk2(k)
           gmi = (mx1+kv2(1,k)) + (mx2+kv2(2,k))*(2*mx1+1)
           gsk3(k) = i * ((2*mx1+1)*(2*mx2+1)+9)+gmi
         ENDDO
         CALL sort(nq2,gsk3,index2)
         DO  k = 1,nq2
            kv3rev(k,1) = kv2(1,index2(k))
            kv3rev(k,2) = kv2(2,index2(k))
            gsk3(k) = sk2(index2(k))
            phi3(k) = phi2(index2(k))
         ENDDO
         DO  k = 1,nq2
            kv2(1,k) = kv3rev(k,1)
            kv2(2,k) = kv3rev(k,2)
            sk2(k) = gsk3(k)
            phi2(k) = phi3(k)
c         if (index2(k).ne.k) write(*,*) ' ic2: ',k,index2(k)
         ENDDO
      ENDIF
c--roa   ....

      WRITE (6,'(/'' nq2='',i4/'' k,kv2(1,2), sk2, phi2''
     &     /(3i4,f10.5,f10.5))')
     &  nq2,(k,kv2(1,k),kv2(2,k),sk2(k),phi2(k),k=1,nq2)

      ELSE
       nq2 = 2
       kv2 = 0 ; sk2 = 0 ; phi2 = 0
      ENDIF ! film
c
c     three dimensional stars
c
      WRITE (6,'(/'' bmat(3,3),mx3='',f10.5,i5)') bmat(3,3),mx3

      IF (mx3.GT.k3d) THEN
         WRITE ( 6,FMT=8000) mx3,k3d
         WRITE (16,FMT=8000) mx3,k3d
         STOP 'mx3.gt.k3d'
      ENDIF
 8000 FORMAT('   mx3.gt.k3d:',2i6)

      nk3 = 2*mx3 + 1
      nq3 = 0
      ig = 0
      rgphs = 0.0
      igfft = 0
      pgfft = 0.0
      gmax2 = gmax * gmax
      izmin =  mx3
      izmax = -mx3

      x_dim: DO k1 = mx1,-mx1,-1
        kv(1) = k1
        y_dim: DO k2 = mx2,-mx2,-1
          kv(2) = k2
          z_dim: DO k3 = mx3,-mx3,-1
            IF ( ig(k1,k2,k3) .NE. 0 ) CYCLE z_dim  ! belongs to another star
            kv(3) = k3

            DO j = 1,3
               g(j) = kv(1)*bmat(1,j) + kv(2)*bmat(2,j) +
     +                kv(3)*bmat(3,j)
            ENDDO
            s = g(1)**2 + g(2)**2 + g(3)**2
            IF (s.LT.gmax2) THEN               ! within gmax-sphere
               nq3 = nq3 + 1
               IF (nq3.GT.n3d) THEN
                 WRITE (6,'("nq3 = ",i5," > n3d =",i5)') nq3,n3d
                 STOP 'strgn: nq3.GT.n3d'
               ENDIF
               DO j = 1,3
                  kv3(j,nq3) = kv(j)
               ENDDO
               sk3(nq3) = sqrt(s)
               IF (k3.LT.izmin) izmin = k3
               IF (k3.GT.izmax) izmax = k3
               CALL spgrot(
     >                     nop,symor,tpi,mrot,tau,invtab,
     >                     kv,
     <                     kr,phas)
               DO n = 1,nop
                  IF (mxx1.lt.kr(1,n)) mxx1 = kr(1,n)
                  IF (mxx2.lt.kr(2,n)) mxx2 = kr(2,n)
                  ig(kr(1,n),kr(2,n),kr(3,n)) = nq3
               ENDDO
            ENDIF

          ENDDO z_dim
        ENDDO y_dim
      ENDDO x_dim
c
c--->    sort for increasing length sk3
c
      CALL sort(nq3,sk3,index)
      DO k = 1,nq3
         kv3rev(k,1) = kv3(1,index(k))
         kv3rev(k,2) = kv3(2,index(k))
         kv3rev(k,3) = kv3(3,index(k))
         gsk3(k) = sk3(index(k))
      ENDDO
      DO k = 1,nq3
         kv3(1,k) = kv3rev(k,1)
         kv3(2,k) = kv3rev(k,2)
         kv3(3,k) = kv3rev(k,3)
         sk3(k) = gsk3(k)
      ENDDO
c
c--roa   ....
c......sort stars of equal length 3D .....
      IF (ordered.OR.to_ordered) THEN
         i=1
         gla=0.
         gsk3(1)=0.
         eps=1.e-10
         DO  k = 2,nq3
           IF (sk3(k)-gla.GE.eps) i=i+1
           gla = sk3(k)
           gmi = (mx1+kv3(1,k)) +
     +           (mx2+kv3(2,k))*(2*mx1+1) +
     +           (mx3+kv3(3,k))*(2*mx1+1)*(2*mx2+1)
           gsk3(k) = i * (9.+(2*mx1+1)*(2*mx2+1)*(2*mx3+1)) + gmi
         ENDDO
         CALL sort(nq3,gsk3,index3)
         DO  k = 1,nq3
            kv3rev(k,1) = kv3(1,index3(k))
            kv3rev(k,2) = kv3(2,index3(k))
            kv3rev(k,3) = kv3(3,index3(k))
            gsk3(k) = sk3(index3(k))
         ENDDO
         DO  k = 1,nq3
            kv3(1,k) = kv3rev(k,1)
            kv3(2,k) = kv3rev(k,2)
            kv3(3,k) = kv3rev(k,3)
            sk3(k) = gsk3(k)
c           if (index3(k).ne.k) write(*,*) ' ic: ',k,index3(k)
         ENDDO
      ENDIF
c--roa   ....
c
c--->  determine true gmax and change old gmax to new gmax
c
      WRITE (6,8060) gmax, sk3(nq3)
      gmax = sk3(nq3)
 8060 FORMAT (/,1x,'old gmax    =',f10.5, '(a.u.)**(-1) ==>  new gmax  '
     >       ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new E_cut   =',
     >            f10.5,' Ry')
!
c--->    generate all star members
c+gu
      kidx=0
c-gu
      rgphs(:,:,:) = 0.0
!
c     set up the pointers and phases for 3d stars
!
      DO  k = 1,nq3
        CALL spgrot(
     >               nop,symor,tpi,mrot,tau,invtab,
     >               kv3(1,k),
     <               kr,phas)
        DO n = 1,nop
c+gu
c -->       set up the igfft(*,3) array as (1d) fft-pointer:
c
c           star ------------> g-vector ------------> fft-grid & phase
c                igfft(*,1)             igfft(*,2)           igfft(*,3)
c
c           size of fft-grid is chosen to be ( 3*k1d x 3*k2d x 3*k3d )
c
c+gu
           new=.true.
           DO n1 = 1,n-1
             norm=(kr(1,n)-kr(1,n1))**2 +
     +            (kr(2,n)-kr(2,n1))**2 +
     +            (kr(3,n)-kr(3,n1))**2
             IF (norm.EQ.0) new = .false.
           ENDDO

           IF (new) THEN

             kfx = kr(1,n)
             kfy = kr(2,n)
             kfz = kr(3,n)
             IF (kfx.LT.0) kfx = kfx+nfftx
             IF (kfy.LT.0) kfy = kfy+nffty
             IF (kfz.LT.0) kfz = kfz+nfftz
!
! -->            store the number of the star, its position 
!                  on fft-grid and phase

             kfft=kfx + kfy*nfftx + kfz*nfftxy
             igfft(kidx,1)=k
             igfft(kidx,2)=kfft
             pgfft(kidx)=phas(n)
             kidx=kidx+1

           ENDIF
c-gu
           ig(kr(1,n),kr(2,n),kr(3,n)) = k
           rgphs(kr(1,n),kr(2,n),kr(3,n)) = 
     +     rgphs(kr(1,n),kr(2,n),kr(3,n)) + phas(n)

        ENDDO
      ENDDO
c
c -->            now for 2d - stars
c
      star2d%i2g = 0
      igfft2 = 0
      pgfft2 = 0.0
      kidx2 = 0

      IF (film) THEN

      DO k = 1,nq2

        kv(1) = kv2(1,k)
        kv(2) = kv2(2,k)
        kv(3) = 0
        CALL spgrot(
     >               nop2,symor,tpi,mrot,tau,invtab,
     >               kv,
     <               kr,phas)

        DO n = 1,nop2
          new=.true.

          DO n1=1,n-1
            norm=(kr(1,n)-kr(1,n1))**2 +
     +           (kr(2,n)-kr(2,n1))**2 +
     +           (kr(3,n)-kr(3,n1))**2
            IF (norm.eq.0) new=.false.
          ENDDO

          IF (new) THEN

            kfx = kr(1,n)
            kfy = kr(2,n)
            gfx = bmat(1,1)*kfx+bmat(2,1)*kfy
            gfy = bmat(1,2)*kfx+bmat(2,2)*kfy
            IF (kfx.LT.0) kfx = kfx+nfftx
            IF (kfy.LT.0) kfy = kfy+nffty
            kfft=kfx + kfy*nfftx

            igfft2(kidx2,1) = k
            igfft2(kidx2,2) = kfft
            pgfft2(kidx2)   = phas(n)
c+guta
            IF (igrd.NE.0) THEN
c                   pgft2x: exp(i*(gfx,gfy,gfz)*tau)*gfx.
c                        y                             y.
c                   pgft2xx: exp(i*(gfx,gfy,gfz)*tau)*gfx*gfx.
c                        yy                             y   y
c                        xy                             x   y

               pgft2x(kidx2)  = phas(n)*gfx
               pgft2y(kidx2)  = phas(n)*gfy
               pgft2xx(kidx2) = phas(n)*gfx*gfx
               pgft2yy(kidx2) = phas(n)*gfy*gfy
               pgft2xy(kidx2) = phas(n)*gfx*gfy
            ENDIF
c-guta
            kidx2 = kidx2 + 1
          ENDIF

          star2d%i2g(kr(1,n),kr(2,n)) = k
          star2d%r2gphs(kr(1,n),kr(2,n)) = 
     +                 star2d%r2gphs(kr(1,n),kr(2,n)) + phas(n)

        ENDDO     ! n -> nop2
      ENDDO       ! k -> nq2 
      
c    
! Array igvac determines for each 2D star and k_z value the 3D star
!
      star2d%igvac = 0
      DO k = 1,nq2
        k1 = kv2(1,k) ; k2 = kv2(2,k)
        DO k3 = -mx3, mx3 
          star2d%igvac(k,k3) = ig(k1,k2,k3) 
        ENDDO
      ENDDO
c
c--->    store number of star with respect to z-index in igz
c
!
! rebuild igz and ig2
!
      star2d%ig2 = 0
      igz = 0
      DO k = 1,nq3
        DO k2 = 1,nq2
          DO k3 = -k3d, k3d ! ngz
            IF (star2d%igvac(k2,k3) == k) THEN
              star2d%ig2(k) = k2
              igz(k) = k3 + k3d
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ELSE
         star2d%ig2 = 0 ; igz = 0 ; star2d%igvac = 0 ; igz = 0
         kimax2= 0 ; igfft2 = 0; pgfft2 = 0.0  ; nstr2 = 0 
         IF (igrd.NE.0) THEN
           pgft2x = 0.0 ; pgft2y = 0.0
           pgft2xx = 0.0 ; pgft2yy = 0.0 ; pgft2xy = 0.0
         ENDIF
      ENDIF ! film

      ngz = izmax - izmin + 1
      kimax=kidx ! -1
      kimax2=kidx2 ! -1
c
c     count number of members for each star
c     nstr2 ... members of 2-dim stars
c
      nstr = 0
      DO k3 = -mx3,mx3
         DO k2 = -mxx2,mxx2
            DO k1 = -mxx1,mxx1
               k = ig(k1,k2,k3)
               IF ( k .NE. 0 ) THEN
                  nstr(k) = nstr(k) + 1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      nstr2 = 0
    
      IF (film) THEN

      DO k2 = -mxx2,mxx2
         DO k1 = -mxx1,mxx1
            k = star2d%i2g(k1,k2)
            IF ( k .NE. 0 ) THEN
               nstr2(k) = nstr2(k) + 1
            ENDIF
         ENDDO
      ENDDO

      ENDIF
!
! normalize phases:
!
      IF (symor) THEN
        rgphs(:,:,:) = 1.0
        star2d%r2gphs(:,:) = 1.0
      ELSE
        pon = 1.0 / nop
        DO k3 = -mx3,mx3
          DO k2 = -mxx2,mxx2
            DO k1 = -mxx1,mxx1

              kfx = k1 ; kfy = k2; kfz = k3
              k = ig(k1,k2,k3)

              IF (k.GT.0) THEN
                IF (kfx.LT.0) kfx = kfx+nfftx
                IF (kfy.LT.0) kfy = kfy+nffty
                IF (kfz.LT.0) kfz = kfz+nfftz
                kfft=kfx + kfy*nfftx + kfz*nfftxy
                rgphs(k1,k2,k3) = rgphs(k1,k2,k3) * nstr(k)*pon
                kidx = -90
                DO i = 1, kimax-1
                  IF ( igfft(i,2) == kfft ) kidx = i
                ENDDO
                IF ( kidx > 0 ) pgfft(kidx)=rgphs(k1,k2,k3)
              ENDIF 

            ENDDO
          ENDDO
        ENDDO
!
! 2D phases
!
        IF (film) THEN
        pon2 = 1.0 / nop2
        DO k2 = -mxx2,mxx2
          DO k1 = -mxx1,mxx1

            kfx = k1 ; kfy = k2
            k = star2d%i2g(k1,k2)

            IF ( k.GT.0 ) THEN
              gfx = bmat(1,1)*kfx+bmat(2,1)*kfy
              gfy = bmat(1,2)*kfx+bmat(2,2)*kfy
              IF (kfx.LT.0) kfx = kfx+nfftx
              IF (kfy.LT.0) kfy = kfy+nffty
              kfft=kfx + kfy*nfftx 
              star2d%r2gphs(k1,k2)=star2d%r2gphs(k1,k2)*nstr2(k)*pon2
              kidx2 = -90
              DO i = 1, kimax2-1
                IF ( igfft2(i,2) == kfft ) kidx2 = i
              ENDDO
              IF ( kidx2 > 0 ) THEN
                pgfft2(kidx2)  = star2d%r2gphs(k1,k2)
                IF (igrd.NE.0) THEN
                  pgft2x(kidx2)  = star2d%r2gphs(k1,k2)*gfx
                  pgft2y(kidx2)  = star2d%r2gphs(k1,k2)*gfy
                  pgft2xx(kidx2) = star2d%r2gphs(k1,k2)*gfx*gfx
                  pgft2yy(kidx2) = star2d%r2gphs(k1,k2)*gfy*gfy
                  pgft2xy(kidx2) = star2d%r2gphs(k1,k2)*gfx*gfy
                ENDIF
              ENDIF
            ENDIF

          ENDDO
        ENDDO
        ENDIF

      ENDIF

      mx1 = mxx1
      mx2 = mxx2
      nk1 = 2*mx1 + 1
      nk2 = 2*mx2 + 1
c
c--->    write /str0/ and /str1/ to unit 51
c
      REWIND 51
      WRITE (51) gmax,nq3,nq2,ngz,izmin,izmax,mx1,mx2,mx3,nk1,nk2,nk3
      IF (igrd.EQ.0) THEN
         WRITE (51) nstr,nstr2,rgphs,sk3,sk2,phi2,kv3,kv2,
     +            ig,star2d%igvac,star2d%ig2,star2d%i2g,star2d%r2gphs,
     +                         kimax,igfft,pgfft,kimax2,igfft2,pgfft2
      ELSE
         WRITE (51) nstr,nstr2,rgphs,sk3,sk2,phi2,kv3,kv2,
     +            ig,star2d%igvac,star2d%ig2,star2d%i2g,star2d%r2gphs,
     +                         kimax,igfft,pgfft,kimax2,igfft2,pgfft2,
     +                           pgft2x,pgft2y,pgft2xx,pgft2yy,pgft2xy
      ENDIF

  270 CONTINUE
      CLOSE (51)
c
c-->  listing
c
      WRITE (16,FMT=8010) gmax,nq3,nq2,ngz,izmin,izmax,nk1,nk2,nk3
 8010 FORMAT (' gmax=',f10.6,/,' nq3=  ',i5,/,' nq2=  ',i5,/,' ngz=  ',
     +       i5,/,' izmin=',i5,/,' izmax=',i5,/,' nk1=  ',i5,/,
     +       ' nk2=  ',i5,/,' nk3=  ',i5,/)
      WRITE (16,FMT=8020) mx1,mx2
 8020 FORMAT (' mx1= ',i5,/,' mx2= ',i5,/)
      WRITE (6,FMT=8030)
 8030 FORMAT (/,/,/,'   s t a r   l i s t',/)

c     k1d,k2d,k3d should be half of nk1,2,3.

      WRITE (6,FMT=8010) gmax,nq3,nq2,ngz,izmin,izmax,nk1,nk2,nk3
      WRITE (6,'('' mx1,mx2,mx3='',3i3)') mx1,mx2,mx3
      WRITE (6,'('' kimax2,kimax='',2i7,'', (start from 0)'')') kimax2,
     &  kimax

      IF (film) THEN
        WRITE (6,FMT=8040)
 8040   FORMAT(/4x,'no.',5x,'kv3',9x,'sk3',9x,'sk2',5x,
     &    'ig2',1x,'igz',1x,'nstr',2x,'nstr2'/)
      ELSE
        WRITE (6,FMT=8041)
 8041   FORMAT(/4x,'no.',5x,'kv3',9x,'sk3',7x,'nstr'/)
      ENDIF

      ned1=9
      nint=30
      DO k = 1,ned1
        IF (film) THEN
        WRITE (6,FMT=8050) k,(kv3(j,k),j=1,3),sk3(k),
     +                     sk2(star2d%ig2(k)),star2d%ig2(k),
     +                     igz(k),nstr(k),nstr2(star2d%ig2(k))
        ELSE
        WRITE (6,FMT=8051) k,(kv3(j,k),j=1,3),sk3(k),nstr(k)
        ENDIF
      ENDDO
 8050 FORMAT (1x,i5,3i4,2f12.6,i4,i3,2i6)
 8051 FORMAT (1x,i5,3i4,f12.6,i6)

      DO k = ned1+1,nq3,nint
        IF (film) THEN
        WRITE (6,FMT=8050) k,(kv3(j,k),j=1,3),sk3(k),
     +                     sk2(star2d%ig2(k)),star2d%ig2(k),
     +                     igz(k),nstr(k),nstr2(star2d%ig2(k))
        ELSE
        WRITE (6,FMT=8051) k,(kv3(j,k),j=1,3),sk3(k),nstr(k)
        ENDIF
        kdone = k
      ENDDO

      IF (kdone.LT.nq3) THEN
        IF (film) THEN
        WRITE (6,FMT=8050) nq3,(kv3(j,nq3),j=1,3),
     +                     sk3(nq3),sk2(star2d%ig2(nq3)),
     +                     star2d%ig2(nq3),igz(nq3),
     +                     nstr(nq3),nstr2(star2d%ig2(nq3))
        ELSE
        WRITE (6,FMT=8051) nq3,(kv3(j,nq3),j=1,3),
     +                     sk3(nq3),nstr(nq3)
        ENDIF
      ENDIF

c.....roa....do sorting for cdn...
      IF (to_ordered.AND..NOT.ordered) THEN
        CALL unor2or(
     >               index2,index3,'ordered*',
     >               n3d,n2d,jspd,jmtd,nlhd,ntypd,nmzd,nmzxyd,natd,
     >               nq3,nq2,jspins,jri,nlh,ntype,nmz,nmzxy,neq,
     >               ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >               nvac,zatom,dx,rmt)
        STOP ' starts.. '
      ENDIf
c
c.....roa....do unsorting for cdn...
      IF (to_unordered.AND.ordered) THEN
c
c..invert..index2......
         DO k = 1,nq2
           gsk3(k)=real(index2(k))
         ENDDO
         CALL sort(nq2,gsk3,index)
         DO k = 1,nq2
           index2(k)=index(k)
         ENDDO
c..invert..index3......
         DO k = 1,nq3
           gsk3(k)=real(index3(k))
         ENDDO
         CALL sort(nq3,gsk3,index)
         DO k = 1,nq3
           index3(k)=index(k)
         ENDDO
         CALL unor2or(
     >                index2,index3,'unorder*',
     >                n3d,n2d,jspd,jmtd,nlhd,ntypd,nmzd,nmzxyd,natd,
     >                nq3,nq2,jspins,jri,nlh,ntype,nmz,nmzxy,neq,
     >                ntypsd,ntypsy,z1,delz,invs,invs2,film,
     >                nvac,zatom,dx,rmt)
         STOP ' starts..unordered '
c.....roa....do sorting for cdn...
      ENDIF

      END SUBROUTINE strgn
!----------------------------------------------------------------
      END MODULE m_strgn
