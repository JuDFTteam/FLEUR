MODULE m_strgn
  USE m_juDFT
  !
  !     *********************************************************
  !     generate two- and three-dimensional stars
  !     for slab geometry
  !     e. wimmer   nov.1984    c.l.fu  1987
  !     implementation of new box-dimension: to treat nonorthogonal
  !     lattice systems
  !     S. Bl"ugel, IFF, 17.Nov.97
  !     *********************************************************
CONTAINS
  SUBROUTINE strgn1(&
       &                  stars,sym,atoms,&
       &                  vacuum,sphhar,input,cell,xcpot)
    !
    USE m_spgrot
    USE m_angle
    USE m_types
    USE m_boxdim
    USE m_sort
    IMPLICIT NONE
    TYPE(t_stars),INTENT(INOUT)  :: stars
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_xcpot),INTENT(IN)     :: xcpot

    !     ..
    !     ..
    !     .. Local Scalars ..
    REAL arltv1,arltv2,arltv3,s
    REAL gmi,gla,eps 
    REAL gfx,gfy,pon,pon2
    INTEGER j,k,k1,k2,k3,m0,mxx1,mxx2,n
    INTEGER ned1,nint,kdone,i,ngz,izmin,izmax
    LOGICAL NEW,l_cdn1
    INTEGER kfx,kfy,kfz,kidx,nfftx,nffty,nfftz,kfft
    INTEGER nfftxy,norm,n1,kidx2,k2i
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: gsk3(:)
    INTEGER, ALLOCATABLE :: ig2p(:),INDEX(:),index3(:),kv3rev(:,:)
    REAL g(3),phi3(stars%ng2),phi
    COMPLEX phas(sym%nop)
    INTEGER kr(3,sym%nop),kv(3)
    INTEGER index2(stars%ng2)
    INTEGER :: igz(stars%ng3)
    !     ..
    !     ..
    !
    !
    nfftx = 3*stars%mx1
    nffty = 3*stars%mx2
    nfftz = 3*stars%mx3
    nfftxy= 9*stars%mx1*stars%mx2

    ALLOCATE (gsk3(stars%ng3),INDEX(stars%ng3),index3(stars%ng3),kv3rev(stars%ng3,3))

    !
    WRITE (*,*) ' stars are always ordered '
    !
    !--->    read in information if exists
    !
    OPEN (51,file='stars',form='unformatted',status='unknown')
    REWIND 51
    !
    READ (51,END=10,err=10) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax,stars%mx1,&
         &  stars%mx2,stars%mx3

    IF (xcpot%igrd.EQ.0) THEN
       READ (51,END=10,err=10) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,&
            &               stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2
    ELSE
       READ (51,END=10,err=10) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,&
            &               stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2,&
            &                         stars%ft2_gfx,stars%ft2_gfy
    ENDIF

    GOTO 270

10  CONTINUE
    !
 
    mxx1 = 0
    mxx2 = 0
    stars%ng2 = 0
    kv(3) = 0

    DO  k1 = stars%mx1,-stars%mx1,-1
       kv(1) = k1
       k2_loop:DO  k2 = stars%mx2,-stars%mx2,-1
          kv(2) = k2
          DO j = 1,2
             g(j) = kv(1)*cell%bmat(1,j) + kv(2)*cell%bmat(2,j)
          ENDDO
          s = SQRT(g(1)**2+g(2)**2)
          !--->   determine the angle of the G_{||} needed in odim calculations
          !+odim YM cute little 'angle' is written by me to determine the phi2
          phi = angle(g(1),g(2))
          !-odim
          !
          !--->   check if generated vector belongs to a star already
          !--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel
          !
          IF (s.LT.stars%gmax) THEN      
             CALL spgrot(&
                  &                     sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
                  &                     kv,&
                  &                     kr)
             DO n = 1,sym%nop2
                IF (mxx1.LT.kr(1,n)) mxx1 = kr(1,n)
                IF (mxx2.LT.kr(2,n)) mxx2 = kr(2,n)
             ENDDO
             DO k = 1,stars%ng2 
                DO n = 1,sym%nop2
                   IF (kr(1,n).EQ.stars%kv2(1,k) .AND.&
                        &                   kr(2,n).EQ.stars%kv2(2,k)) CYCLE k2_loop
                ENDDO
             ENDDO
             !--->    new representative found
             stars%ng2 = stars%ng2 + 1
             IF (stars%ng2.GT.stars%ng2) THEN
                WRITE (6,8070) stars%ng2,stars%ng2
                CALL juDFT_error("nq2.GT.n2d",calledby="strgn")
             ENDIF
             DO j = 1,2
                stars%kv2(j,stars%ng2) = kv(j)
             ENDDO
             stars%sk2(stars%ng2) = s
             stars%phi2(stars%ng2) = phi
          ENDIF
       ENDDO k2_loop
    ENDDO
8070 FORMAT ('nq2 = ',i5,' > n2d =',i5)

    IF( mxx1.GT.stars%mx1) THEN
       WRITE (6,'(/'' mxx1.gt.k1d. mxx1,k1d='',2i4)')  mxx1,stars%mx1
       CALL juDFT_error("mxx1.gt.k1d",calledby="strgn")
    ENDIF

    IF (mxx2.GT.stars%mx2) THEN
       WRITE (6,'(/'' mxx2.gt.k2d. mxx2,k2d='',2i4)')  mxx2,stars%mx2
       CALL juDFT_error("mxx2.gt.k2d",calledby="strgn")
    ENDIF

    !--->    sort for increasing length sk2

    CALL sort(stars%ng2,stars%sk2,index)
    DO k = 1,stars%ng2
       kv3rev(k,1) = stars%kv2(1,INDEX(k))
       kv3rev(k,2) = stars%kv2(2,INDEX(k))
       gsk3(k) = stars%sk2(INDEX(k))
       phi3(k) = stars%phi2(INDEX(k))
    ENDDO
    DO k = 1,stars%ng2
       stars%kv2(1,k) = kv3rev(k,1)
       stars%kv2(2,k) = kv3rev(k,2)
       stars%sk2(k) = gsk3(k)
       stars%phi2(k) = phi3(k)
    ENDDO

    !--roa   ....
    !......sort stars of equal length 2D .....
    i=1
    gla=0.
    gsk3(1)=0.0
    eps=1.e-10
    DO  k = 2,stars%ng2
       IF (stars%sk2(k)-gla.GE.eps) i=i+1
       gla=stars%sk2(k)
       gmi = (stars%mx1+stars%kv2(1,k)) + (stars%mx2+stars%kv2(2,k))*(2*stars%mx1+1)
       gsk3(k) = i * ((2*stars%mx1+1)*(2*stars%mx2+1)+9)+gmi
    ENDDO
    CALL sort(stars%ng2,gsk3,index2)
    DO  k = 1,stars%ng2
       kv3rev(k,1) = stars%kv2(1,index2(k))
       kv3rev(k,2) = stars%kv2(2,index2(k))
       gsk3(k) = stars%sk2(index2(k))
       phi3(k) = stars%phi2(index2(k))
    ENDDO
    DO  k = 1,stars%ng2
       stars%kv2(1,k) = kv3rev(k,1)
       stars%kv2(2,k) = kv3rev(k,2)
       stars%sk2(k) = gsk3(k)
       stars%phi2(k) = phi3(k)
       !         if (index2(k).ne.k) write(*,*) ' ic2: ',k,index2(k)
    ENDDO
    !--roa   ....

    WRITE (6,'(/'' nq2='',i4/'' k,kv2(1,2), sk2, phi2''&
         &     /(3i4,f10.5,f10.5))')&
         &  stars%ng2,(k,stars%kv2(1,k),stars%kv2(2,k),stars%sk2(k),stars%phi2(k),k=1,stars%ng2)
    !
    !     three dimensional stars
    !
    stars%ng3 = 0
    DO k3 = -stars%mx3,stars%mx3
       DO k2 = -stars%mx2,stars%mx2
          DO k1 = -stars%mx1,stars%mx1
             stars%ig(k1,k2,k3) = 0
             stars%rgphs(k1,k2,k3) = cmplx(0.0,0.0)
          ENDDO
       ENDDO
    ENDDO
    !+gu
    stars%igfft2(:,:)=0
    stars%igfft(:,:)=0
    stars%pgfft=0.0
    stars%pgfft2=0.0
    !-gu
    WRITE (6,'(/'' bmat(3,3),mx3='',f10.5,i5)') cell%bmat(3,3),stars%mx3

    IF (stars%mx3.GT.stars%mx3) THEN
       WRITE ( 6,FMT=8000) stars%mx3,stars%mx3
       WRITE (16,FMT=8000) stars%mx3,stars%mx3
       CALL juDFT_error("mx3.gt.k3d",calledby="strgn")
    ENDIF
8000 FORMAT('   mx3.gt.k3d:',2i6)

    m0 = -stars%mx3
    !     zrfs,invs: z-reflection, inversion.
    IF (sym%zrfs .OR. sym%invs) m0 = 0
    izmin = stars%mx3
    izmax = -stars%mx3

    DO  k2 = 1,stars%ng2
       DO  k3 = m0,stars%mx3
          s = SQRT(stars%sk2(k2)**2+ (k3*cell%bmat(3,3))**2)
          !
          !--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel
          IF (s.LT.stars%gmax) THEN      
             !
             stars%ng3 = stars%ng3 + 1
             IF (stars%ng3.GT.stars%ng3) THEN
                WRITE (6,*) stars%ng3,stars%ng3
                CALL juDFT_error("nq3.GT.n3d",calledby="strgn")
             ENDIF
             IF (k3.LT.izmin) izmin = k3
             IF (k3.GT.izmax) izmax = k3
             DO j = 1,2
                stars%kv3(j,stars%ng3) = stars%kv2(j,k2)
             ENDDO
             stars%kv3(3,stars%ng3) = k3
             stars%ig2(stars%ng3) = k2
             stars%sk3(stars%ng3) = s
          ENDIF
       ENDDO
    ENDDO

    !--->    sort for increasing length sk3

    CALL sort(stars%ng3,stars%sk3,index)

    ALLOCATE (ig2p(stars%ng3))

    DO k = 1,stars%ng3
       kv3rev(k,1) = stars%kv3(1,INDEX(k))
       kv3rev(k,2) = stars%kv3(2,INDEX(k))
       kv3rev(k,3) = stars%kv3(3,INDEX(k))
       gsk3(k) = stars%sk3(INDEX(k))
       ig2p(k) = stars%ig2(INDEX(k))
    ENDDO
    DO k = 1,stars%ng3
       stars%kv3(1,k) = kv3rev(k,1)
       stars%kv3(2,k) = kv3rev(k,2)
       stars%kv3(3,k) = kv3rev(k,3)
       stars%sk3(k) = gsk3(k)
       stars%ig2(k) = ig2p(k)
    ENDDO
    !
    !--roa   ....
    !......sort stars of equal length 3D .....
    i=1
    gla=0.
    gsk3(1)=0.
    eps=1.e-10
    DO  k = 2,stars%ng3
       IF (stars%sk3(k)-gla.GE.eps) i=i+1
       gla = stars%sk3(k)
       gmi = (stars%mx1+stars%kv3(1,k)) +&
            &           (stars%mx2+stars%kv3(2,k))*(2*stars%mx1+1) +&
            &           (stars%mx3+stars%kv3(3,k))*(2*stars%mx1+1)*(2*stars%mx2+1)
       gsk3(k) = i * (9.+(2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)) + gmi
    ENDDO
    CALL sort(stars%ng3,gsk3,index3)
    DO  k = 1,stars%ng3
       kv3rev(k,1) = stars%kv3(1,index3(k))
       kv3rev(k,2) = stars%kv3(2,index3(k))
       kv3rev(k,3) = stars%kv3(3,index3(k))
       gsk3(k) = stars%sk3(index3(k))
       ig2p(k) = stars%ig2(index3(k))
    ENDDO
    DO  k = 1,stars%ng3
       stars%kv3(1,k) = kv3rev(k,1)
       stars%kv3(2,k) = kv3rev(k,2)
       stars%kv3(3,k) = kv3rev(k,3)
       stars%sk3(k) = gsk3(k)
       stars%ig2(k) = ig2p(k)
       !           if (index3(k).ne.k) write(*,*) ' ic: ',k,index3(k)
    ENDDO

    DEALLOCATE (ig2p)

    !--roa   ....
    !
    !--->  determine true gmax and change old gmax to new gmax
    !
    WRITE (6,8060) stars%gmax, stars%sk3(stars%ng3)
    stars%gmax = stars%sk3(stars%ng3)
8060 FORMAT (/,1x,'old stars%gmax    =',f10.5, '(a.u.)**(-1) ==>  new stars%gmax  '&
         &       ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new E_cut   =',&
         &            f10.5,' Ry')
    !
    !--->    store number of star with respect to z-index in igz
    !
    !ngz = izmax - izmin + 1
    !DO k = 1,stars%ng3
    !   igz(k) = stars%kv3(3,k) - izmin + 1
    !ENDDO
    igz=0;ngz=0;izmin=0;izmax=0

    !--->    generate all star members
    !+gu
    kidx=0
    kidx2=0
    !-gu
    stars%rgphs(:,:,:) = cmplx(0.0,0.0)
    DO  k = 1,stars%ng3

       CALL spgrot(&
            &               sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
            &               stars%kv3(:,k),&
            &               kr,phas)

       IF (stars%kv3(3,k).EQ.0) THEN

          DO  n = 1,sym%nop2
             !+gu
             ! -->       set up the igfft(*,3) array as (1d) fft-pointer:
             !
             !           star ------------> g-vector ------------> fft-grid & phase
             !                igfft(*,1)             igfft(*,2)           igfft(*,3)
             !
             !           size of fft-grid is chosen to be ( 3*k1d x 3*k2d x 3*k3d )
             !
             NEW=.TRUE.

             DO n1=1,n-1
                norm=(kr(1,n)-kr(1,n1))**2 +&
                     &                (kr(2,n)-kr(2,n1))**2 +&
                     &                (kr(3,n)-kr(3,n1))**2
                IF (norm.EQ.0) NEW=.FALSE.
             ENDDO

             IF (NEW) THEN

                kfx = kr(1,n)
                kfy = kr(2,n)
                kfz = kr(3,n)
                !+guta
                gfx = cell%bmat(1,1)*kfx+cell%bmat(2,1)*kfy+cell%bmat(3,1)*kfz 
                gfy = cell%bmat(1,2)*kfx+cell%bmat(2,2)*kfy+cell%bmat(3,2)*kfz
                !-guta
                IF (kfx.LT.0) kfx = kfx+nfftx
                IF (kfy.LT.0) kfy = kfy+nffty
                IF (kfz.LT.0) kfz = kfz+nfftz
                kfft = kfx + kfy*nfftx + kfz*nfftxy
                !
                ! -->            store the number of the star, its position 
                !c                 on fft-grid and phase
                !
                stars%igfft(kidx,1) = k
                stars%igfft(kidx,2) = kfft
                stars%pgfft(kidx)   = phas(n)
                kidx          = kidx+1
                !
                ! -->            now for 2d - stars
                !
                kfft=kfx + kfy*nfftx
                DO k2 = 1,stars%ng2
                   IF ((stars%kv3(1,k).EQ.stars%kv2(1,k2)).AND.&
                        &                 (stars%kv3(2,k).EQ.stars%kv2(2,k2))) k2i = k2
                ENDDO
                stars%igfft2(kidx2,1) = k2i
                stars%igfft2(kidx2,2) = kfft
                stars%pgfft2(kidx2)   = phas(n)
                !+guta
                IF (xcpot%igrd.NE.0) THEN
                   !!                   pgft2x: exp(i*(gfx,gfy,gfz)*tau)*gfx.
                   !!                        y                             y.
                   !!                   pgft2xx: exp(i*(gfx,gfy,gfz)*tau)*gfx*gfx.
                   !!                        yy                             y   y
                   !!                        xy                             x   y

                   stars%ft2_gfx(kidx2)  = gfx
                   stars%ft2_gfy(kidx2)  = gfy
                   !                    stars%pgft2xx(kidx2) = phas(n)*gfx*gfx
                   !                    stars%pgft2yy(kidx2) = phas(n)*gfy*gfy
                   !                    stars%pgft2xy(kidx2) = phas(n)*gfx*gfy
                ENDIF
                !-guta
                kidx2=kidx2+1

             ENDIF
             !-gu
             stars%ig(kr(1,n),kr(2,n),kr(3,n)) = k
             stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) = &
                  &         stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) + phas(n)

          ENDDO

       ELSE
          !        here: kv3(3,k) =/= 0 

          DO  n = 1,sym%nop
             !+gu
             NEW=.TRUE.
             DO n1 = 1,n-1
                norm=(kr(1,n)-kr(1,n1))**2 +&
                     &                (kr(2,n)-kr(2,n1))**2 +&
                     &                (kr(3,n)-kr(3,n1))**2
                IF (norm.EQ.0) NEW = .FALSE.
             ENDDO

             IF (NEW) THEN

                kfx = kr(1,n)
                kfy = kr(2,n)
                kfz = kr(3,n)
                IF (kfx.LT.0) kfx = kfx+nfftx
                IF (kfy.LT.0) kfy = kfy+nffty
                IF (kfz.LT.0) kfz = kfz+nfftz

                kfft=kfx + kfy*nfftx + kfz*nfftxy
                stars%igfft(kidx,1)=k
                stars%igfft(kidx,2)=kfft
                stars%pgfft(kidx)=phas(n)
                kidx=kidx+1

             ENDIF
             !-gu
             stars%ig(kr(1,n),kr(2,n),kr(3,n)) = k
             stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) = &
                  &         stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) + phas(n)

          ENDDO

       ENDIF

    ENDDO
    !    
    stars%kimax=kidx-1
    stars%kimax2=kidx2-1
    !
    !     count number of members for each star
    !     nstr2 ... members of 2-dim stars
    !
    DO k = 1,stars%ng3
       stars%nstr2(stars%ig2(k)) = 0
       stars%nstr(k) = 0
    ENDDO

    DO k3 = -stars%mx3,stars%mx3
       DO k2 = -mxx2,mxx2
          DO k1 = -mxx1,mxx1
             k = stars%ig(k1,k2,k3)
             IF ( k .NE. 0 ) THEN
                stars%nstr(k) = stars%nstr(k) + 1
                stars%nstr2(stars%ig2(k)) = stars%nstr2(stars%ig2(k)) +(1-min0(iabs(k3),1))
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
    ! normalize phases:
    !
    IF (sym%symor) THEN
       stars%rgphs(:,:,:) = cmplx(1.0,0.0)
    ELSE
       pon = 1.0 / sym%nop
       pon2 = 1.0 / sym%nop2
       DO k3 = -stars%mx3,stars%mx3
          DO k2 = -mxx2,mxx2
             DO k1 = -mxx1,mxx1

                kfx = k1 ; kfy = k2; kfz = k3
                k = stars%ig(k1,k2,k3)

                IF (kfx.LT.0) kfx = kfx+nfftx
                IF (kfy.LT.0) kfy = kfy+nffty
                IF (kfz.LT.0) kfz = kfz+nfftz
                kfft=kfx + kfy*nfftx + kfz*nfftxy
                IF (k.GT.0) THEN
                   IF (stars%kv3(3,k).EQ.0) THEN
                      stars%rgphs(k1,k2,k3) = stars%rgphs(k1,k2,k3) * stars%nstr(k)*pon2
                   ELSE
                      stars%rgphs(k1,k2,k3) = stars%rgphs(k1,k2,k3) * stars%nstr(k)*pon
                   ENDIF
                   kidx = -90
                   DO i = 1, stars%kimax
                      IF ( stars%igfft(i,2) == kfft ) kidx = i
                   ENDDO
                   IF ( kidx > 0 ) stars%pgfft(kidx)=stars%rgphs(k1,k2,k3)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF (stars%mx1 .ne. mxx1 .or.  stars%mx2 .ne. mxx2) CALL judft_error("BUG in strgn")
    !
    !--->    write /str0/ and /str1/ to unit 51
    !
    REWIND 51
    WRITE (51) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax,stars%mx1,stars%mx2,stars%mx3
    IF (xcpot%igrd.EQ.0) THEN
       WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,stars%ig,stars%ig2,igz,&
            &                   stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2
    ELSE
       WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,stars%ig,stars%ig2,igz,&
            &                   stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2,&
            &                    stars%ft2_gfx,stars%ft2_gfy
    ENDIF

270 CONTINUE
    CLOSE (51)
    !
    !-->  listing
    !
    WRITE (16,FMT=8010) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax
8010 FORMAT (' gmax=',f10.6,/,' nq3=  ',i5,/,' nq2=  ',i5,/,' ngz=  ',&
         &       i5,/,' izmin=',i5,/,' izmax=',i5/)
    WRITE (16,FMT=8020) stars%mx1,stars%mx2
8020 FORMAT (' mx1= ',i5,/,' mx2= ',i5,/)
    WRITE (6,FMT=8030)
8030 FORMAT (/,/,/,'   s t a r   l i s t',/)

    WRITE (6,FMT=8010) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax
    WRITE (6,'('' mx1,mx2,mx3='',3i3)') stars%mx1,stars%mx2,stars%mx3
    WRITE (6,'('' kimax2,kimax='',2i7,'', (start from 0)'')') stars%kimax2,&
         &  stars%kimax

    WRITE (6,FMT=8040)
8040 FORMAT(/4x,'no.',5x,'kv3',9x,'sk3',9x,'sk2',5x,&
         &  'ig2',1x,'igz',1x,'nstr',2x,'nstr2'/)

    ned1=9
    nint=30
    DO k = 1,ned1
       WRITE (6,FMT=8050) k,(stars%kv3(j,k),j=1,3),stars%sk3(k),&
            &                     stars%sk2(stars%ig2(k)),&
            &                     stars%ig2(k),igz(k),stars%nstr(k),stars%nstr2(stars%ig2(k))
    ENDDO
8050 FORMAT (1x,i5,3i4,2f12.6,i4,i3,2i6)

    DO k = ned1+1,stars%ng3,nint
       WRITE (6,FMT=8050) k,(stars%kv3(j,k),j=1,3),stars%sk3(k),&
            &                     stars%sk2(stars%ig2(k)),&
            &                     stars%ig2(k),igz(k),stars%nstr(k),stars%nstr2(stars%ig2(k))
       kdone = k
    ENDDO

    IF (kdone.LT.stars%ng3) THEN
       WRITE (6,FMT=8050) stars%ng3,(stars%kv3(j,stars%ng3),j=1,3),&
            &                     stars%sk3(stars%ng3),stars%sk2(stars%ig2(stars%ng3)),&
            &                     stars%ig2(stars%ng3),igz(stars%ng3),stars%nstr(stars%ng3),stars%nstr2(stars%ig2(stars%ng3))
    ENDIF

    DEALLOCATE (gsk3,index,index3,kv3rev)

  END SUBROUTINE strgn1
  !----------------------------------------------------------------
  SUBROUTINE strgn2(&
       &                  stars,sym,atoms,&
       &                  vacuum,sphhar,input,cell,xcpot)
    USE m_boxdim
    USE m_sort
    USE m_spgrot
    USE m_constants
    USE m_types
    IMPLICIT NONE
    TYPE(t_stars),INTENT(INOUT)  :: stars
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
    !     ..
    !     .. Local Scalars ..
    REAL arltv1,arltv2,arltv3,s
    REAL gmi,gla,eps,gmax2,pon
    REAL gfx,gfy
    INTEGER j,k,k1,k2,k3,m0,mxx1,mxx2,mxx3,n
    INTEGER ned1,nint,kdone,i
    LOGICAL NEW,l_cdn1
    INTEGER kfx,kfy,kfz,kidx,nfftx,nffty,nfftz,kfft
    INTEGER nfftxy,norm,n1,kidx2,k2i,izmin,izmax,ngz
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: gsk3(:)
    INTEGER, ALLOCATABLE :: INDEX(:),index3(:),kv3rev(:,:)
    REAL g(3)
    COMPLEX phas(sym%nop)
    INTEGER kr(3,sym%nop),kv(3)
    INTEGER index2(stars%ng2)
    INTEGER :: igz(stars%ng3)

    !     ..
    !
    nfftx = 3*stars%mx1
    nffty = 3*stars%mx2
    nfftz = 3*stars%mx3
    nfftxy= 9*stars%mx1*stars%mx2

    ALLOCATE (gsk3(stars%ng3),INDEX(stars%ng3),index3(stars%ng3),kv3rev(stars%ng3,3))

    !
    WRITE (*,*) ' stars are always ordered '

    !--->    read in information if exists
    !
    OPEN (51,file='stars',form='unformatted',status='unknown')
    REWIND 51
    !
    !
    READ (51,END=10,err=10) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax,stars%mx1,&
         &  stars%mx2,stars%mx3

    IF (xcpot%igrd.EQ.0) THEN
       READ (51,END=10,err=10) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,stars%ig,&
            &               stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2
    ELSE
       READ (51,END=10,err=10) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,stars%ig,&
            &               stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2,&
            &                         stars%ft2_gfx,stars%ft2_gfy
    ENDIF

    GOTO 270

10  CONTINUE
    !

    !
    !
    mxx1 = 0
    mxx2 = 0
    mxx3 = 0
    stars%ng3 = 0
    stars%ig = 0
    gmax2 = stars%gmax * stars%gmax

    x_dim: DO k1 = stars%mx1,-stars%mx1,-1
       kv(1) = k1
       y_dim: DO k2 = stars%mx2,-stars%mx2,-1
          kv(2) = k2
          z_dim: DO k3 = stars%mx3,-stars%mx3,-1
             IF ( stars%ig(k1,k2,k3) .NE. 0 ) CYCLE z_dim
             kv(3) = k3

             DO j = 1,3
                g(j) = kv(1)*cell%bmat(1,j) + kv(2)*cell%bmat(2,j) +&
                     &                kv(3)*cell%bmat(3,j)
             ENDDO
             s = g(1)**2 + g(2)**2 + g(3)**2
             !
             !--->   check if generated vector belongs to a star already
             IF (s.LT.gmax2) THEN      
                !
                !--->    new representative found
                CALL spgrot(&
                     &                     sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
                     &                     kv,&
                     &                     kr)
                stars%ng3 = stars%ng3 + 1
                IF (stars%ng3.GT.stars%ng3) THEN
                   WRITE (6,'("nq3 = ",i5," > n3d =",i5)') stars%ng3,stars%ng3
                   CALL juDFT_error("nq3>n3d",calledby ="strgn")
                ENDIF
                DO j = 1,3
                   stars%kv3(j,stars%ng3) = kv(j)
                ENDDO
                stars%sk3(stars%ng3) = SQRT(s)
                DO n = 1,sym%nop
                   IF (mxx1.LT.kr(1,n)) mxx1 = kr(1,n)
                   IF (mxx2.LT.kr(2,n)) mxx2 = kr(2,n)
                   IF (mxx3.LT.kr(3,n)) mxx3 = kr(3,n)
                   stars%ig(kr(1,n),kr(2,n),kr(3,n)) = stars%ng3
                ENDDO
             ENDIF
          ENDDO z_dim
       ENDDO y_dim
    ENDDO x_dim

    IF( mxx1.GT.stars%mx1) THEN
       WRITE (6,'(/'' mxx1.gt.k1d. mxx1,k1d='',2i4)')  mxx1,stars%mx1
       CALL juDFT_error("mxx1>k1d",calledby ="strgn")
    ENDIF

    IF (mxx2.GT.stars%mx2) THEN
       WRITE (6,'(/'' mxx2.gt.k2d. mxx2,k2d='',2i4)')  mxx2,stars%mx2
       CALL juDFT_error("mxx2>k2d",calledby ="strgn")
    ENDIF

    IF (mxx3.GT.stars%mx3) THEN
       WRITE (6,'(/'' mxx3.gt.k3d. mxx3,k3d='',2i4)')  mxx3,stars%mx3
       CALL juDFT_error("mxx3>k3d",calledby ="strgn")
    ENDIF

    !--->    sort for increasing length sk3
    !
    CALL sort(stars%ng3,stars%sk3,index)
    DO k = 1,stars%ng3
       kv3rev(k,1) = stars%kv3(1,INDEX(k))
       kv3rev(k,2) = stars%kv3(2,INDEX(k))
       kv3rev(k,3) = stars%kv3(3,INDEX(k))
       gsk3(k) = stars%sk3(INDEX(k))
    ENDDO
    DO k = 1,stars%ng3
       stars%kv3(1,k) = kv3rev(k,1)
       stars%kv3(2,k) = kv3rev(k,2)
       stars%kv3(3,k) = kv3rev(k,3)
       stars%sk3(k) = gsk3(k)
    ENDDO
    !
    !--roa   ....
    !......sort stars of equal length 3D .....
    i=1
    gla=0.
    gsk3(1)=0.
    eps=1.e-10
    DO  k = 2,stars%ng3
       IF (stars%sk3(k)-gla.GE.eps) i=i+1
       gla = stars%sk3(k)
       gmi = (stars%mx1+stars%kv3(1,k)) +&
            &           (stars%mx2+stars%kv3(2,k))*(2*stars%mx1+1) +&
            &           (stars%mx3+stars%kv3(3,k))*(2*stars%mx1+1)*(2*stars%mx2+1)
       gsk3(k) = i * (9.+(2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)) + gmi
    ENDDO
    CALL sort(stars%ng3,gsk3,index3)
    DO  k = 1,stars%ng3
       kv3rev(k,1) = stars%kv3(1,index3(k))
       kv3rev(k,2) = stars%kv3(2,index3(k))
       kv3rev(k,3) = stars%kv3(3,index3(k))
       gsk3(k) = stars%sk3(index3(k))
    ENDDO
    DO  k = 1,stars%ng3
       stars%kv3(1,k) = kv3rev(k,1)
       stars%kv3(2,k) = kv3rev(k,2)
       stars%kv3(3,k) = kv3rev(k,3)
       stars%sk3(k) = gsk3(k)
       !           if (index3(k).ne.k) write(*,*) ' ic: ',k,index3(k)
    ENDDO
    !--roa   ....
    !
    !--->  determine true gmax and change old gmax to new gmax
    !
    WRITE (6,8060) stars%gmax, stars%sk3(stars%ng3)
    stars%gmax = stars%sk3(stars%ng3)
8060 FORMAT (/,1x,'gmax    =',f10.5, '(a.u.)**(-1) ==>  new gmax  '&
         &       ,'  =',f10.5,'(a.u.)**(-1) ',/,t38,'==>  new E_cut   =',&
         &            f10.5,' Ry')
    !
    !--->    store number of star with respect to z-index in igz
    !
    !ngz = izmax - izmin + 1
    !DO k = 1,stars%ng3
    !   igz(k) = stars%kv3(3,k) - izmin + 1
    !ENDDO
    igz=0;ngz=0;izmax=0;izmin=0
    !--->    generate all star members
    !+gu
    kidx=0
    kidx2=0
    DO k3 = -mxx3,mxx3
       DO k2 = -mxx2,mxx2
          DO k1 = -mxx1,mxx1
             stars%ig(k1,k2,k3) = 0
          ENDDO
       ENDDO
    ENDDO

    !-gu
    !
    ! sum over phases
    !
    stars%rgphs(:,:,:) = cmplx(0.0,0.0)
    starloop: DO k = 1,stars%ng3

       CALL spgrot(&
            &               sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,&
            &               stars%kv3(:,k),&
            &               kr,phas)
       !
       ! -->    set up the igfft(*,3) array as (1d) fft-pointer:
       !
       !        star ------------> g-vector ------------> fft-grid & phase
       !             igfft(*,1)             igfft(*,2)           igfft(*,3)
       !
       !        size of fft-grid is chosen to be ( 3*k1d x 3*k2d x 3*k3d )
       !
       ops: DO n = 1,sym%nop

          NEW=.TRUE.
          DO n1 = 1,n-1
             norm=(kr(1,n)-kr(1,n1))**2 +&
                  &             (kr(2,n)-kr(2,n1))**2 +&
                  &             (kr(3,n)-kr(3,n1))**2
             IF (norm.EQ.0) NEW = .FALSE.
          ENDDO

          IF (NEW) THEN

             kfx = kr(1,n)
             kfy = kr(2,n)
             kfz = kr(3,n)
             IF (kfx.LT.0) kfx = kfx+nfftx
             IF (kfy.LT.0) kfy = kfy+nffty
             IF (kfz.LT.0) kfz = kfz+nfftz

             kfft=kfx + kfy*nfftx + kfz*nfftxy
             stars%igfft(kidx,1)=k
             stars%igfft(kidx,2)=kfft
             stars%pgfft(kidx)=phas(n)
             kidx=kidx+1

          ENDIF
          stars%ig(kr(1,n),kr(2,n),kr(3,n)) = k
          stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) = &
               &      stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) + phas(n)

       ENDDO ops
    ENDDO starloop
    !    
    stars%kimax=kidx-1
    !
    !     count number of members for each star
    !
    DO k = 1,stars%ng3
       stars%nstr(k) = 0
    ENDDO

    DO k3 = -mxx3,mxx3
       DO k2 = -mxx2,mxx2
          DO k1 = -mxx1,mxx1
             k = stars%ig(k1,k2,k3)
             IF ( k .NE. 0 ) THEN
                stars%nstr(k) = stars%nstr(k) + 1
             ENDIF
             !               DO k = 1,nq3
             !                  IF (ig(k1,k2,k3).eq.k) THEN
             !                     nstr(k) = nstr(k) + 1
             !                  ENDIF
             !               ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    ! normalize phases:
    !
    IF (sym%symor) THEN
       stars%rgphs(:,:,:) = cmplx(1.0,0.0)
    ELSE
       pon = 1.0 / sym%nop
       DO k3 = -mxx3,mxx3
          DO k2 = -mxx2,mxx2
             DO k1 = -mxx1,mxx1
                k = stars%ig(k1,k2,k3)
                IF (k.GT.0) THEN
                   stars%rgphs(k1,k2,k3) = stars%rgphs(k1,k2,k3) * stars%nstr(k)*pon
                   kfx = k1 ; kfy = k2; kfz = k3
                   IF (kfx.LT.0) kfx = kfx+nfftx
                   IF (kfy.LT.0) kfy = kfy+nffty
                   IF (kfz.LT.0) kfz = kfz+nfftz
                   kfft=kfx + kfy*nfftx + kfz*nfftxy
                   kidx = -90
                   DO i = 1, stars%kimax
                      IF ( stars%igfft(i,2) == kfft ) kidx = i
                   ENDDO
                   IF (kidx > 0 ) stars%pgfft(kidx)=stars%rgphs(k1,k2,k3)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    if ( stars%mx1 .ne. mxx1 .or. stars%mx2 .ne. mxx2 .or. stars%mx3 .ne. mxx3 ) call &
         judft_error("BUG 1 in strgen") 
    stars%ng2 = 2 ; stars%kv2 = 0 ; stars%ig2 = 0 ; igz = 0 ; stars%kimax2= 0 ; stars%igfft2 = 0
    stars%sk2 = 0.0 ; stars%pgfft2 = 0.0  ; stars%nstr2 = 0
    IF (xcpot%igrd.NE.0) THEN
       stars%ft2_gfx = 0.0 ; stars%ft2_gfy = 0.0 
    ENDIF
    !
    !--->    write /str0/ and /str1/ to unit 51
    !
    REWIND 51
    WRITE (51) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax,stars%mx1,stars%mx2,stars%mx3
    IF (xcpot%igrd.EQ.0) THEN
       WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,stars%ig,stars%ig2,igz,&
            &                   stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2
    ELSE
       WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,stars%ig,stars%ig2,igz,&
            &                   stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,stars%igfft2,stars%pgfft2,&
            &                    stars%ft2_gfx,stars%ft2_gfy
    ENDIF

270 CONTINUE
    CLOSE (51)
    !
    !-->  listing
    WRITE (16,FMT=8010) stars%gmax,stars%ng3,izmin,izmax
8010 FORMAT (' gmax=',f10.6,/,' nq3=  ',i5,/,&
         &       ' izmin=',i5,/,' izmax=',i5,/)
    WRITE (16,FMT=8020) stars%mx1,stars%mx2,stars%mx3
8020 FORMAT (' mx1= ',i5,/,' mx2= ',i5,' mx3= ',i5,/)
    WRITE (6,FMT=8030)
8030 FORMAT (/,/,/,'   s t a r   l i s t',/)


    WRITE (6,FMT=8010) stars%gmax,stars%ng3,izmin,izmax
    WRITE (6,'('' mx1,mx2,mx3='',3i3)') stars%mx1,stars%mx2,stars%mx3
    WRITE (6,'('' kimax2,kimax='',2i7,'', (start from 0)'')') stars%kimax2,&
         &  stars%kimax

    WRITE (6,FMT=8040)
8040 FORMAT(/4x,'no.',5x,'kv3',9x,'sk3',7x,'nstr'/)

    ned1=9
    nint=30
    DO k = 1,ned1
       WRITE (6,FMT=8050) k,(stars%kv3(j,k),j=1,3),stars%sk3(k),stars%nstr(k)
    ENDDO
8050 FORMAT (1x,i5,3i4,f12.6,i6)

    DO k = ned1+1,stars%ng3,nint
       WRITE (6,FMT=8050) k,(stars%kv3(j,k),j=1,3),stars%sk3(k),stars%nstr(k)
       kdone = k
    ENDDO

    IF (kdone.LT.stars%ng3) THEN
       WRITE (6,FMT=8050) stars%ng3,(stars%kv3(j,stars%ng3),j=1,3),&
            &                     stars%sk3(stars%ng3),stars%nstr(stars%ng3)
    ENDIF


    DEALLOCATE (gsk3,index,index3,kv3rev)

  END SUBROUTINE strgn2
END MODULE m_strgn
