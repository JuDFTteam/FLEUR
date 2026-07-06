!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_stars 
      use m_juDFT
      IMPLICIT NONE
      PRIVATE 
      PUBLIC gf_stargen, gf_write_stars,gf_read_stars,deallocate_stars 
      CONTAINS 
      !<-- S: gf_stargen                                                
                                                                        
      SUBROUTINE gf_stargen(                                            &
     &     sym,bmat,gmax,                                               &
     &     stars,boxed)                                                 
!                                                                       
!     *********************************************************         
!     generate two- and three-dimensional stars                         
!     intrinic call to dimensioning routine                             
!     See strgn1 from FLEUR for more details                            
!     Daniel Wortmann, 2002                                             
!     *********************************************************         
!     for slab geometry                                                 
!     e. wimmer   nov.1984    c.l.fu  1987                              
!     implementation of new box-dimension: to treat nonorthogonal       
!     lattice systems                                                   
!     S. Bl"ugel, IFF, 17.Nov.97                                        
!     *********************************************************c        
!                                                                       
!     Added option boxed to generate stars filling the whole box and not
!     only the usual sphere                                             
!     (last modified: 04-11-04)        D.Wortmann                       
!     *********************************************************c        
      use m_juDFT 
      USE m_gf_math,ONLY:sort 
      USE m_fleur_interface
      USE m_gf_types 
      IMPLICIT NONE 
!     ..                                                                
!     .. Scalar Arguments ..                                            
      REAL,    INTENT (IN) :: gmax 
      TYPE(t_sym),INTENT(IN)  :: sym 
      TYPE(t_stars),INTENT(INOUT):: stars 
      LOGICAL,OPTIONAL,INTENT(IN)::boxed 
!     ..                                                                
!     .. Array Arguments ..                                             
      REAL,    INTENT (IN)  :: bmat(3,3) 
!     ..                                                                
!     .. Local Scalars ..                                               
      REAL s 
      REAL gmi,gla,eps 
      REAL gfx,gfy 
      INTEGER j,k,k1,k2,k3,m0,n 
      INTEGER ned1,n_int,kdone,i 
      LOGICAL new,pe0 
      INTEGER kfx,kfy,kfz,kidx,nfftx,nffty,nfftz,kfft 
      INTEGER nfftxy,norm,n1,kidx2,k2i 
!     ..                                                                
!     .. Local Arrays ..                                                
      REAL g(3) 
      INTEGER kv(3) 
      REAL,ALLOCATABLE    :: gsk3(:),phas(:) 
      INTEGER,ALLOCATABLE :: ig2p(:),index1(:),kr(:,:),kv3rev(:,:) 
      INTEGER,ALLOCATABLE :: index2(:),index3(:) 
                                                                        
                                                                        
      LOGICAL ::l_boxed 
                                                                        
#ifdef CPP_MPI                                                          
      INTEGER ::irank,ierr(3) 
      include 'mpif.h' 
#endif                                                                  
!     ..                                                                
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC iabs,int,min0,sqrt 
!     ..                                                                
!                                                                       
#ifdef CPP_MPI                                                          
      CALL  MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr) 
      pe0=(irank==0) 
#else                                                                   
      pe0=.TRUE. 
#endif                                                                  
      IF (PRESENT(boxed)) THEN 
         IF (boxed) THEN 
            l_boxed=.TRUE. 
         ELSE 
            l_boxed=.FALSE. 
         ENDIF 
      ELSE 
         l_boxed=.FALSE. 
      ENDIF 
      CALL priv_star_dim(gmax,bmat,sym,stars,l_boxed) 
      ALLOCATE(ig2p(stars%nq3),index1(stars%nq3),kr(3,size(sym%tau,2))  &
     &     ,kv3rev(stars%nq3,3),index2(stars%nq2),index3(stars%nq3)     &
     &     ,gsk3(stars%nq3),phas(size(sym%tau,2)))                      
!                                                                       
      nfftx = 3*stars%mx1 
      nffty = 3*stars%mx2 
      nfftz = 3*stars%mx3 
      nfftxy= nfftx*nffty 
                                                                        
      stars%nq2 = 0 
      kv(3) = 0 
                                                                        
      DO  k1 = stars%mx1,-stars%mx1,-1 
         kv(1) = k1 
         k2loop:DO k2 = stars%mx2,-stars%mx2,-1 
             kv(2) = k2 
             g(1:2)=matmul(kv(1:2),bmat(1:2,1:2)) 
             s = sqrt(g(1)**2+g(2)**2) 
!                                                                       
!--->   check if generated vector belongs to a star already             
!--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel   
!                                                                       
             IF (s<=gmax) THEN 
                CALL fleur_spgrot(                                      &
     &               sym,kv,                                            &
     &               kr,phas)                                           
                DO k = 1,stars%nq2 
                   DO n = 1,sym%nop2 
                      IF (all(kr(1:2,n)==stars%kv2(:,k))) CYCLE k2loop 
                   ENDDO 
                ENDDO 
!--->    new representative found                                       
            stars%nq2 = stars%nq2 + 1 
            IF (stars%nq2>size(stars%kv2,2)) THEN 
               WRITE(*,*) k1,k2,stars%mx1,stars%mx2 
               WRITE(*,*) s,gmax 
               WRITE(*,8070) stars%nq2,size(stars%kv2,2) 
                CALL juDFT_error("2D-Dim error in gf_genstars",calledby="gf_stars.F90")
            ENDIF 
            stars%kv2(1:2,stars%nq2) = kv(1:2) 
            stars%sk2(stars%nq2) = s 
         ENDIF 
      ENDDO k2loop 
      ENDDO 
 8070 FORMAT ('stars%nq2 = ',i5,' > n2d =',i5) 
      IF (stars%nq2/=SIZE(stars%kv2,2)) THEN 
         WRITE(*,*) "Some problem with creating 2D-stars" 
         WRITE(*,*) "Dimension different from actual number of stars" 
         WRITE(*,*) stars%nq2,"/=",SIZE(stars%kv2,2) 
         CALL juDFT_error("m_gf_stars:gf_stargen nq2") 
      ENDIF 
                                                                        
!--->    sort for increasing length sk2                                 
                                                                        
      CALL sort(stars%sk2,index1(:stars%nq2)) 
      DO k = 1,stars%nq2 
         kv3rev(k,1:2) = stars%kv2(1:2,index1(k)) 
         gsk3(k) = stars%sk2(index1(k)) 
      ENDDO 
      DO k = 1,stars%nq2 
         stars%kv2(1:2,k) = kv3rev(k,1:2) 
         stars%sk2(k) = gsk3(k) 
      ENDDO 
!--roa   ....                                                           
!......sort stars of equal length 2D .....                              
      i=1 
      gla=0. 
      gsk3(1)=0.0 
      eps=1.E-10 
      DO  k = 2,stars%nq2 
         IF (stars%sk2(k)-gla>=eps) i=i+1 
         gla=stars%sk2(k) 
         gmi = (stars%mx1+stars%kv2(1,k)) + (stars%mx2+stars%kv2(2,k))  &
     &        *(2*stars%mx1+1)                                          
         gsk3(k) = i * ((2*stars%mx1+1)*(2*stars%mx2+1)+9)+gmi 
      ENDDO 
      CALL sort(gsk3(:stars%nq2),index2(:stars%nq2)) 
      DO  k = 1,stars%nq2 
         kv3rev(k,1) = stars%kv2(1,index2(k)) 
         kv3rev(k,2) = stars%kv2(2,index2(k)) 
         gsk3(k) = stars%sk2(index2(k)) 
      ENDDO 
      DO  k = 1,stars%nq2 
         stars%kv2(1,k) = kv3rev(k,1) 
         stars%kv2(2,k) = kv3rev(k,2) 
         stars%sk2(k) = gsk3(k) 
      ENDDO 
                                                                        
      IF (pe0) WRITE (6                                                 &
     &     ,'(/'' stars%nq2='',i4/'' k,kv2(1,2), sk2''/(3i4,f10.5))'    &
     &     )stars%nq2,(k,stars%kv2(1,k),stars%kv2(2,k),stars%sk2(k),k   &
     &     =1,stars%nq2)                                                
!                                                                       
!     three dimensional stars                                           
!                                                                       
      stars%nq3 = 0 
                                                                        
      stars%ig=0 
      stars%rgphs=0.0 
      stars%igfft2=0 
      stars%igfft2=0 
      stars%igfft=0 
      stars%igfft=0 
      stars%pgfft=0.0 
      stars%pgfft2=0.0 
!-gu                                                                    
      IF (pe0) WRITE (6,'(/'' bmat(3,3),stars%mx3='',f10.5,i5)') bmat(3 &
     &     ,3),stars%mx3                                                
                                                                        
                                                                        
      m0 = -stars%mx3 
!     sym%zrfs,sym%invs: z-reflection, inversion.                       
      IF (sym%zrfs .OR. sym%invs) m0 = 0 
      stars%izmin = stars%mx3 
      stars%izmax = -stars%mx3 
                                                                        
      DO 150 k2 = 1,stars%nq2 
         DO 140 k3 = m0,stars%mx3 
            s = sqrt(stars%sk2(k2)**2+ (k3*bmat(3,3))**2) 
!                                                                       
!--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel   
            IF (s<=gmax.OR.l_boxed) THEN 
               stars%nq3 = stars%nq3 + 1 
               IF (stars%nq3>size(stars%kv3,2)) THEN 
                  WRITE(*,*) "stars%nq3=",stars%nq3,                    &
     &                 ">",size(stars%kv3,2)                            
                   CALL juDFT_error("3D-Dim error in gf_genstars",calledby="gf_stars.F90")
               ENDIF 
               stars%kv3(1:2,stars%nq3) = stars%kv2(1:2,k2) 
               stars%kv3(3,stars%nq3) = k3 
               stars%ig2(stars%nq3) = k2 
               stars%sk3(stars%nq3) = s 
            ENDIF 
  140    ENDDO 
  150 ENDDO 
      IF (stars%nq3 /= SIZE(stars%kv3,2)) THEN 
         WRITE(*,*) "Some problem with creating 3D-stars" 
         WRITE(*,*) "Dimension different from actual number of stars" 
         CALL juDFT_error("m_gf_stars:gf_stargen nq3") 
      ENDIF 
                                                                        
      stars%izmin=minval(stars%kv3(3,:)) 
      stars%izmax=maxval(stars%kv3(3,:)) 
!                                                                       
!--->    sort for increasing length sk3                                 
!                                                                       
      CALL sort(stars%sk3,index1) 
      DO k = 1,stars%nq3 
         kv3rev(k,1:3) = stars%kv3(1:3,index1(k)) 
         gsk3(k) = stars%sk3(index1(k)) 
         ig2p(k) = stars%ig2(index1(k)) 
      ENDDO 
      DO k = 1,stars%nq3 
         stars%kv3(1:3,k) = kv3rev(k,1:3) 
         stars%sk3(k) = gsk3(k) 
         stars%ig2(k) = ig2p(k) 
      ENDDO 
                                                                        
!......sort stars of equal length 3D .....                              
                                                                        
      i=1 
      gla=0. 
      gsk3(1)=0. 
      eps=1.E-10 
      DO  k = 2,stars%nq3 
         IF (stars%sk3(k)-gla>=eps) i=i+1 
         gla = stars%sk3(k) 
         gmi = (stars%mx1+stars%kv3(1,k)) +                             &
     &        (stars%mx2+stars%kv3(2,k))*(2*stars%mx1+1) +              &
     &        (stars%mx3+stars%kv3(3,k))*(2*stars%mx1+1)*(2*stars%mx2+1)
         gsk3(k) = i * (9.+(2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3 &
     &        +1)) + gmi                                                
      ENDDO 
      CALL sort(gsk3,index3) 
      DO  k = 1,stars%nq3 
         kv3rev(k,1:3) = stars%kv3(1:3,index3(k)) 
         gsk3(k) = stars%sk3(index3(k)) 
         ig2p(k) = stars%ig2(index3(k)) 
      ENDDO 
      DO  k = 1,stars%nq3 
         stars%kv3(1:3,k) = kv3rev(k,1:3) 
         stars%sk3(k) = gsk3(k) 
         stars%ig2(k) = ig2p(k) 
      ENDDO 
                                                                        
!                                                                       
!--->    store number of star with respect to z-index in igz            
!                                                                       
      stars%ngz = stars%izmax - stars%izmin + 1 
      stars%igz(:) = stars%kv3(3,:) - stars%izmin + 1 
                                                                        
!--->    generate all star members                                      
!+gu                                                                    
      kidx=0 
      kidx2=0 
!-gu                                                                    
      DO 210 k = 1,stars%nq3 
                                                                        
         CALL fleur_spgrot(sym,stars%kv3(:,k),                          &
     &        kr,phas)                                                  
                                                                        
         IF (stars%kv3(3,k)==0) THEN 
                                                                        
            DO 190 n = 1,sym%nop2 
!+gu                                                                    
! -->       set up the stars%igfft(*,3) array as (1d) fft-pointer:      
!                                                                       
!           star ------------> g-vector ------------> fft-grid & phase  
!                stars%igfft(*,1)             stars%igfft(*,2)          
! stars%igfft(*,3)                                                      
!                                                                       
!           size of fft-grid is chosen to be ( 3*k1d x 3*k2d x 3*k3d )  
!                                                                       
               new=.TRUE. 
                                                                        
               DO n1=1,n-1 
                  norm=(kr(1,n)-kr(1,n1))**2 +                          &
     &                 (kr(2,n)-kr(2,n1))**2 +                          &
     &                 (kr(3,n)-kr(3,n1))**2                            
                  IF (norm==0) new=.FALSE. 
               ENDDO 
                                                                        
               IF (new) THEN 
                                                                        
                  kfx = kr(1,n) 
                  kfy = kr(2,n) 
                  kfz = kr(3,n) 
!+guta                                                                  
                  gfx = bmat(1,1)*kfx+bmat(2,1)*kfy+bmat(3,1)*kfz 
                  gfy = bmat(1,2)*kfx+bmat(2,2)*kfy+bmat(3,2)*kfz 
                                                                        
!-guta                                                                  
                  IF (kfx<0) kfx = kfx+nfftx 
                  IF (kfy<0) kfy = kfy+nffty 
                  IF (kfz<0) kfz = kfz+nfftz 
                  kfft = kfx + kfy*nfftx + kfz*nfftxy 
!                                                                       
! -->            store the number of the star, its position             
!c                 on fft-grid and phase                                
!                                                                       
                  stars%igfft(kidx,1) = k 
                  stars%igfft(kidx,2) = kfft 
                  stars%pgfft(kidx)   = phas(n) 
                  kidx          = kidx+1 
                  IF (kidx>size(stars%pgfft))  CALL juDFT_error("DIMERR",calledby="gf_stars.F90")
!                                                                       
! -->            now for 2d - stars                                     
!                                                                       
                  kfft=kfx + kfy*nfftx 
                  DO k2 = 1,stars%nq2 
                     IF ((stars%kv3(1,k)==stars%kv2(1,k2)).AND.       &
     &                    (stars%kv3(2,k)==stars%kv2(2,k2))) k2i = k2 
                  ENDDO 
                  stars%igfft2(kidx2,1) = k2i 
                  stars%igfft2(kidx2,2) = kfft 
                  stars%pgfft2(kidx2)   = phas(n) 
                  !initialize gga-arrays anyway...                      
                  stars%pgft2x(kidx2)  = phas(n)*gfx 
                  stars%pgft2y(kidx2)  = phas(n)*gfy 
                  stars%pgft2xx(kidx2) = phas(n)*gfx*gfx 
                  stars%pgft2yy(kidx2) = phas(n)*gfy*gfy 
                  stars%pgft2xy(kidx2) = phas(n)*gfx*gfy 
                                                                        
                  kidx2=kidx2+1 
                  IF (kidx2>size(stars%pgft2x))  CALL juDFT_error("DIMERR",calledby="gf_stars.F90")
               ENDIF 
!-gu                                                                    
               stars%ig(kr(1,n),kr(2,n),kr(3,n)) = k 
               stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) = phas(n) 
                                                                        
  190       ENDDO 
                                                                        
         ELSE 
!        here: stars%kv3(3,k) =/= 0                                     
                                                                        
            DO 200 n = 1,sym%nop 
!+gu                                                                    
               new=.TRUE. 
               DO n1 = 1,n-1 
                  norm=(kr(1,n)-kr(1,n1))**2 +                          &
     &                 (kr(2,n)-kr(2,n1))**2 +                          &
     &                 (kr(3,n)-kr(3,n1))**2                            
                  IF (norm==0) new = .FALSE. 
               ENDDO 
                                                                        
               IF (new) THEN 
                                                                        
                  kfx = kr(1,n) 
                  kfy = kr(2,n) 
                  kfz = kr(3,n) 
                  IF (kfx<0) kfx = kfx+nfftx 
                  IF (kfy<0) kfy = kfy+nffty 
                  IF (kfz<0) kfz = kfz+nfftz 
                                                                        
                  kfft=kfx + kfy*nfftx + kfz*nfftxy 
                  stars%igfft(kidx,1)=k 
                  stars%igfft(kidx,2)=kfft 
                  stars%pgfft(kidx)=phas(n) 
                  kidx=kidx+1 
                                                                        
               ENDIF 
!-gu                                                                    
                                                                        
               stars%ig(kr(1,n),kr(2,n),kr(3,n)) = k 
               stars%rgphs(kr(1,n),kr(2,n),kr(3,n)) = phas(n) 
                                                                        
  200       ENDDO 
                                                                        
         ENDIF 
                                                                        
  210 ENDDO 
!                                                                       
      stars%kimax=kidx-1 
      stars%kimax2=kidx2-1 
!                                                                       
!     count number of members for each star                             
!     stars%nstr2 ... members of 2-dim stars                            
!                                                                       
      stars%nstr2 = 0 
      stars%nstr = 0 
      DO k3 = -stars%mx3,stars%mx3 
         DO k2 = -stars%mx2,stars%mx2 
            DO k1 = -stars%mx1,stars%mx1 
               k=stars%ig(k1,k2,k3) 
               IF (k<1) CYCLE 
               stars%nstr(k) = stars%nstr(k) + 1 
               stars%nstr2(stars%ig2(k)) = stars%nstr2(stars%ig2(k      &
     &              )) +(1-min0(iabs(k3),1))                            
                                                                        
            ENDDO 
         ENDDO 
      ENDDO 
!                                                                       
!                                                                       
!-->  listing                                                           
!                                                                       
      IF (pe0)WRITE (16,FMT=8010) gmax,stars%nq3,stars%nq2,stars%ngz    &
     &     ,stars%izmin,stars%izmax,2*stars%mx1+1,2*stars%mx2+1,2       &
     &     *stars%mx3+1                                                 
 8010 FORMAT (' gmax=',f10.6,/,' stars%nq3=  ',i5,/,' stars%nq2=  ',i5,/&
     &     , ' ngz=  ',i5,/,' izmin=',i5,/,' izmax=',i5,/,' nk1=  ',i5,/&
     &     ,'nk2=  ',i5,/,' nk3=  ',i5,/)                               
      IF (pe0)WRITE (16,FMT=8020) stars%mx1,stars%mx2 
 8020 FORMAT (' stars%mx1= ',i5,/,' stars%mx2= ',i5,/) 
      IF (pe0)WRITE (6,FMT=8030) 
 8030 FORMAT (/,/,/,'   s t a r   l i s t',/) 
                                                                        
!     k1d,k2d,k3d should be half of nk1,2,3.                            
                                                                        
      IF (pe0)WRITE (6,FMT=8010) gmax,stars%nq3,stars%nq2,stars%ngz     &
     &     ,stars%izmin,stars%izmax,2*stars%mx1+1,2*stars%mx2+1,2       &
     &     *stars%mx3+1                                                 
      IF (pe0)WRITE (6,'('' stars%mx1,stars%mx2,stars%mx3='',3i3)')     &
     &     stars%mx1,stars%mx2,stars%mx3                                
      IF (pe0)WRITE (6                                                  &
     &     ,'('' stars%kimax2,stars%kimax='',2i7,'', (start from 0)'')')&
     &     stars%kimax2,stars%kimax                                     
                                                                        
      IF (pe0)WRITE (6,FMT=8040) 
 8040 FORMAT(/4x,'no.',5x,'kv3',9x,'sk3',9x,'sk2',5x,                   &
     &     'stars%ig2',1x,'igz',1x,'stars%nstr',2x,'stars%nstr2'/)      
                                                                        
      ned1=9 
      n_int=30 
      DO k = 1,ned1 
         IF (pe0)WRITE (6,FMT=8050) k,(stars%kv3(j,k),j=1,3),stars%sk3(k&
     &        ),stars%sk2(stars%ig2(k)),stars%ig2(k),stars%igz(k)       &
     &        ,stars%nstr(k),stars%nstr2(stars%ig2(k))                  
      ENDDO 
 8050 FORMAT (1x,i5,3i4,2f12.6,i4,i3,2i6) 
                                                                        
      DO k = ned1+1,stars%nq3,n_int 
         IF (pe0)WRITE (6,FMT=8050) k,(stars%kv3(j,k),j=1,3),stars%sk3(k&
     &        ),stars%sk2(stars%ig2(k)),stars%ig2(k),stars%igz(k)       &
     &        ,stars%nstr(k),stars%nstr2(stars%ig2(k))                  
         kdone = k 
      ENDDO 
                                                                        
      IF (kdone<stars%nq3) THEN 
         IF (pe0)WRITE (6,FMT=8050) stars%nq3,(stars%kv3(j,stars%nq3),j &
     &        =1,3),stars%sk3(stars%nq3),stars%sk2(stars%ig2(stars%nq3))&
     &        ,stars%ig2(stars%nq3),stars%igz(stars%nq3)                &
     &        ,stars%nstr(stars%nq3),stars%nstr2(stars%ig2(stars%nq3))  
      ENDIF 
      DEALLOCATE(ig2p,index1,kr,kv3rev,index2,index3,gsk3,phas) 
                                                                        
      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:priv_star_dim                                              
                                                                        
      SUBROUTINE priv_star_dim(                                         &
     &     gmax,bmat,sym,                                               &
     &     stars,l_boxed)                                               
!                                                                       
!     *********************************************************         
!     Determines dimensions of stars and allocates pointers             
!                                                                       
!     Algorithm from FLEUR strgn1_dim                                   
!     *********************************************************         
      USE m_fleur_INTERFACE
      USE m_gf_types 
      USE m_constants 
      IMPLICIT NONE 
                                                                        
      REAL, INTENT(IN)          :: gmax,bmat(3,3) 
      TYPE(t_Sym),INTENT(IN)    :: sym 
      TYPE(T_stars),INTENT(INOUT) :: stars 
      LOGICAL,INTENT(IN)        :: l_boxed 
!                                                                       
      REAL arltv1,arltv2,arltv3,s 
      INTEGER kv(3),mxx1,mxx2 
      INTEGER n,k,k1,k2,k3,m0 
      INTEGER nn2d,nn3d,iofile,ksfft 
      INTEGER,ALLOCATABLE:: kv2(:,:) 
      INTEGER kr(3,sym%nop) 
      REAL    phas(sym%nop),g(3) 
                                                                        
                                                                        
      INTRINSIC int 
                                                                        
!---> Determine Gmax box of size stars%mx1, stars%mx2, stars%mx3,       
!     for which |G(stars%mx1,stars%mx2,stars%mx3)| < Gmax               
!     arltv(i) length of reciprical lattice vector along direction (i)  
!                                                                       
      CALL fleur_boxdim(bmat,arltv1,arltv2,arltv3) 
                                                                        
      stars%mx1 = int(gmax/arltv1) + 1 
      stars%mx2 = int(gmax/arltv2) + 1 
      IF (.NOT.l_boxed) THEN 
         stars%mx3 = int(gmax/arltv3) + 1 
      ELSE 
         stars%mx3 = int(gmax/arltv3) + 1 
      ENDIF 
      !<--make box-dimensions compatible with fft                       
      WRITE(6,*) "Adjusting FFT-Box" 
      WRITE(6,*) "Old Box:", stars%mx1, stars%mx2, stars%mx3 
      iofile = 6 
      ksfft  = 1 
      !>                                                                
      ALLOCATE(kv2(2,(2*stars%mx1+1)*(2*stars%mx2+1))) 
!                                                                       
!     two-dimensional stars                                             
!                                                                       
      mxx1=0 
      mxx2=0 
      stars%nq2 = 0 
      kv(3) = 0 
      DO  k1 = stars%mx1,-stars%mx1,-1 
         kv(1) = k1 
         k2loop:DO  k2 = stars%mx2,-stars%mx2,-1 
            kv(2) = k2 
            g(1:2)=matmul(kv(1:2),bmat(1:2,1:2)) 
            s = sqrt(g(1)**2+g(2)**2) 
            IF (s<=gmax) THEN 
               CALL fleur_spgrot(sym,kv,                                &
     &                     kr,phas)
               DO  n = 1,sym%nop2
                   mxx1=max(mxx1,kr(1,n))
                   mxx2=max(mxx2,kr(2,n))
               enddo
               DO  k = 1,stars%nq2 
                  DO  n = 1,sym%nop2 
                     IF (all(kr(1:2,n)==kv2(1:2,k)))CYCLE k2loop 
                ENDDO 
             ENDDO 
!--->    new representative found                                       
               stars%nq2 = stars%nq2 + 1 
               kv2(1:2,stars%nq2) = kv(1:2) 
            END IF 
         ENDDO k2loop 
      ENDDO 
      stars%mx1=mxx1
      stars%mx2=mxx2
      !stars%mx1 = fleur_ifft235(iofile,ksfft,stars%mx1,2.0)
      !stars%mx2 = fleur_ifft235(iofile,ksfft,stars%mx2,2.0)
      !stars%mx3 = fleur_ifft235(iofile,ksfft,stars%mx3,2.0)
      WRITE(6,*) "New Box:", stars%mx1, stars%mx2, stars%mx3
!                                                                       
!                                                                       
!     three dimensional stars                                           
      stars%nq3 = 0 
      m0 = -stars%mx3 
      IF (sym%zrfs .OR. sym%invs) m0 = 0 
      DO k2 = 1,stars%nq2 
         DO k3 = m0,stars%mx3 
            g=matmul((/ kv2(1,k2),kv2(2,k2),k3 /),bmat) 
            s = sqrt(g(1)**2+g(2)**2+g(3)**2) 
            IF (s<=gmax.OR.l_boxed) THEN 
               stars%nq3 = stars%nq3 + 1 
            END IF 
         END DO 
      END DO 
      DEALLOCATE(kv2) 
!                                                                       
!    Now allocate the storage for the arrays                            
!                                                                       
      nn2d= (2*stars%mx1+1)* (2*stars%mx2+1) 
      nn3d= (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1) 
                                                                        
      ALLOCATE(stars%kv3(3,stars%nq3)) 
      ALLOCATE(stars%sk3(stars%nq3)) 
      ALLOCATE(stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,      &
     &     -stars%mx3:stars%mx3))                                       
      ALLOCATE(stars%nstr(stars%nq3)) 
      ALLOCATE(stars%kv2(2,stars%nq2)) 
      ALLOCATE(stars%sk2(stars%nq2)) 
      ALLOCATE(stars%nstr2(stars%nq2)) 
      ALLOCATE(stars%ig2(stars%nq3)) 
      ALLOCATE(stars%igz(stars%nq3)) 
      ALLOCATE(stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,   &
     &     -stars%mx3:stars%mx3))                                       
      ALLOCATE(stars%igfft(0:nn3d-1,2)) 
      ALLOCATE(stars%igfft2(0:nn2d-1,2)) 
      ALLOCATE(stars%pgfft(0:nn3d-1)) 
      ALLOCATE(stars%pgfft2(0:nn2d-1)) 
      ALLOCATE(stars%pgft2x(0:nn2d-1)) 
      ALLOCATE(stars%pgft2xx(0:nn2d-1)) 
      ALLOCATE(stars%pgft2xy(0:nn2d-1)) 
      ALLOCATE(stars%pgft2yy(0:nn2d-1)) 
      ALLOCATE(stars%pgft2y(0:nn2d-1)) 
                                                                        
      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: gf_write_stars                                            
      SUBROUTINE gf_WRITE_stars(fid,stars) 
!-----------------------------------------------                        
!     write the stars to the hdf-file                                   
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      use m_juDFT 
      USE m_gf_types 
      USE m_hdf_tools 
      USE hdf5 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN)  :: fid 
      TYPE(t_stars),INTENT(IN)   :: stars 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T)      :: gid 
      INTEGER             :: hdferr 
      !>                                                                
                                                                        
      IF (io_groupexists(fid,"stars")) RETURN 
                                                                        
      CALL io_gcreate(fid,"stars",gid,hdferr) 
      !write attributes                                                 
      CALL io_WRITE_att(gid, "gmax",stars%gmax) 
      CALL io_WRITE_att(gid, "gmax_inp",stars%gmax_inp) 
      CALL io_WRITE_att(gid, "nq3",stars%nq3) 
      CALL io_WRITE_att(gid, "nq2",stars%nq2) 
      CALL io_WRITE_att(gid, "kimax",stars%kimax) 
      CALL io_WRITE_att(gid, "kimax2",stars%kimax2) 
      CALL io_WRITE_att(gid, "mx1",stars%mx1) 
      CALL io_WRITE_att(gid, "mx2",stars%mx2) 
      CALL io_WRITE_att(gid, "mx3",stars%mx3) 
      CALL io_WRITE_att(gid, "ngz",stars%ngz) 
      CALL io_WRITE_att(gid, "izmax",stars%izmax) 
      CALL io_WRITE_att(gid, "izmin",stars%izmin) 
                                                                        
      !write large variables                                            
      CALL io_WRITE_var(gid,"kv3",stars%kv3) 
      CALL io_WRITE_var(gid,"sk3",stars%sk3) 
      CALL io_WRITE_var(gid,"ig",stars%ig) 
      CALL io_WRITE_var(gid,"nstr",stars%nstr) 
      CALL io_WRITE_var(gid,"kv2",stars%kv2) 
      CALL io_WRITE_var(gid,"sk2",stars%sk2) 
      CALL io_WRITE_var(gid,"nstr2",stars%nstr2) 
      CALL io_WRITE_var(gid,"ig2",stars%ig2) 
      CALL io_WRITE_var(gid,"igz",stars%igz) 
      CALL io_WRITE_var(gid,"rgphs",stars%rgphs) 
      CALL io_WRITE_var(gid,"igfft",stars%igfft) 
      CALL io_WRITE_var(gid,"igfft2",stars%igfft2) 
      CALL io_WRITE_var(gid,"pgfft",stars%pgfft) 
      CALL io_WRITE_var(gid,"pgfft2",stars%pgfft2) 
      CALL io_WRITE_var(gid,"pgft2xy",stars%pgft2xy) 
      CALL io_WRITE_var(gid,"pgft2x",stars%pgft2x) 
      CALL io_WRITE_var(gid,"pgft2y",stars%pgft2y) 
      CALL io_WRITE_var(gid,"pgft2xx",stars%pgft2xx) 
      CALL io_WRITE_var(gid,"pgft2yy",stars%pgft2yy) 
                                                                        
      CALL io_gclose(gid,hdferr) 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: gf_read_stars                                             
                                                                        
      SUBROUTINE gf_read_stars(fid,stars) 
!-----------------------------------------------                        
!     write the stars to the hdf-file                                   
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      use m_juDFT 
      USE m_gf_types 
      USE m_hdf_tools 
      USE hdf5 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER(HID_T),INTENT(IN)   :: fid 
      TYPE(t_stars),INTENT(OUT)   :: stars 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T)      :: gid 
      INTEGER             :: hdferr 
      !>                                                                
                                                                        
      IF (.NOT.io_groupexists(fid,"stars")) THEN 
         CALL juDFT_error("Could not read stars") 
      ENDIF 
                                                                        
      CALL io_gopen(fid,"stars",gid,hdferr) 
      !read attributes                                                  
      CALL io_READ_att(gid, "gmax",stars%gmax) 
      CALL io_READ_att(gid, "gmax_inp",stars%gmax_inp) 
      CALL io_READ_att(gid, "nq3",stars%nq3) 
      CALL io_READ_att(gid, "nq2",stars%nq2) 
      CALL io_READ_att(gid, "kimax",stars%kimax) 
      CALL io_READ_att(gid, "kimax2",stars%kimax2) 
      CALL io_READ_att(gid, "mx1",stars%mx1) 
      CALL io_READ_att(gid, "mx2",stars%mx2) 
      CALL io_READ_att(gid, "mx3",stars%mx3) 
      CALL io_READ_att(gid, "ngz",stars%ngz) 
      CALL io_READ_att(gid, "izmax",stars%izmax) 
      CALL io_READ_att(gid, "izmin",stars%izmin) 
                                                                        
      CALL ALLOCATE_stars(stars) 
                                                                        
      !read large variables                                             
      CALL io_READ_var(gid,"kv3",stars%kv3) 
      CALL io_READ_var(gid,"sk3",stars%sk3) 
      CALL io_READ_var(gid,"ig",stars%ig) 
      CALL io_READ_var(gid,"nstr",stars%nstr) 
      CALL io_READ_var(gid,"kv2",stars%kv2) 
      CALL io_READ_var(gid,"sk2",stars%sk2) 
      CALL io_READ_var(gid,"nstr2",stars%nstr2) 
      CALL io_READ_var(gid,"ig2",stars%ig2) 
      CALL io_READ_var(gid,"igz",stars%igz) 
      CALL io_READ_var(gid,"rgphs",stars%rgphs) 
      CALL io_READ_var(gid,"igfft",stars%igfft) 
      CALL io_READ_var(gid,"igfft2",stars%igfft2) 
      CALL io_READ_var(gid,"pgfft",stars%pgfft) 
      CALL io_READ_var(gid,"pgfft2",stars%pgfft2) 
      CALL io_READ_var(gid,"pgft2xy",stars%pgft2xy) 
      CALL io_READ_var(gid,"pgft2x",stars%pgft2x) 
      CALL io_READ_var(gid,"pgft2y",stars%pgft2y) 
      CALL io_READ_var(gid,"pgft2xx",stars%pgft2xx) 
      CALL io_READ_var(gid,"pgft2yy",stars%pgft2yy) 
                                                                        
      CALL io_gclose(gid,hdferr) 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: allocate_stars(stars)                                     
      SUBROUTINE allocate_stars(stars) 
!-----------------------------------------------                        
! allocates the storage for stars                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(INOUT) :: stars 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: nn2d,nn3d 
      !>                                                                
      nn2d = (2*stars%mx1+1)* (2*stars%mx2+1) 
      nn3d = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1) 
                                                                        
      ALLOCATE(stars%kv3(3,stars%nq3)) 
      ALLOCATE(stars%sk3(stars%nq3)) 
      ALLOCATE(stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,      &
     &     -stars%mx3:stars%mx3))                                       
      ALLOCATE(stars%nstr(stars%nq3)) 
      ALLOCATE(stars%kv2(2,stars%nq2)) 
      ALLOCATE(stars%sk2(stars%nq2)) 
      ALLOCATE(stars%nstr2(stars%nq2)) 
      ALLOCATE(stars%ig2(stars%nq3)) 
      ALLOCATE(stars%igz(stars%nq3)) 
      ALLOCATE(stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,   &
     &     -stars%mx3:stars%mx3))                                       
      ALLOCATE(stars%igfft(0:nn3d-1,2)) 
      ALLOCATE(stars%igfft2(0:nn2d-1,2)) 
      ALLOCATE(stars%pgfft(0:nn3d-1)) 
      ALLOCATE(stars%pgfft2(0:nn2d-1)) 
      ALLOCATE(stars%pgft2x(0:nn2d-1)) 
      ALLOCATE(stars%pgft2xx(0:nn2d-1)) 
      ALLOCATE(stars%pgft2xy(0:nn2d-1)) 
      ALLOCATE(stars%pgft2yy(0:nn2d-1)) 
      ALLOCATE(stars%pgft2y(0:nn2d-1)) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: deallocate_stars(stars)                                   
      SUBROUTINE deallocate_stars(stars) 
!-----------------------------------------------                        
!  deallocate storage from stars-type                                   
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(INOUT) :: stars 
      !>                                                                
      INTEGER      :: err 
      DEALLOCATE(stars%kv3,STAT = err) 
      DEALLOCATE(stars%sk3,STAT=err) 
      DEALLOCATE(stars%ig,STAT=err) 
      DEALLOCATE(stars%nstr,STAT=err) 
      DEALLOCATE(stars%kv2,STAT=err) 
      DEALLOCATE(stars%sk2,STAT=err) 
      DEALLOCATE(stars%nstr2,STAT=err) 
      DEALLOCATE(stars%ig2,STAT=err) 
      DEALLOCATE(stars%igz,STAT=err) 
      DEALLOCATE(stars%rgphs,STAT=err) 
      DEALLOCATE(stars%igfft,STAT=err) 
      DEALLOCATE(stars%igfft2,STAT=err) 
      DEALLOCATE(stars%pgfft,STAT=err) 
      DEALLOCATE(stars%pgfft2,STAT=err) 
      DEALLOCATE(stars%pgft2x,STAT=err) 
      DEALLOCATE(stars%pgft2xx,STAT=err) 
      DEALLOCATE(stars%pgft2xy,STAT=err) 
      DEALLOCATE(stars%pgft2yy,STAT=err) 
      DEALLOCATE(stars%pgft2y,STAT=err) 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
