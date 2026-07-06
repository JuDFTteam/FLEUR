!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_apws 
      use m_juDFT
      use m_juDFT 
      IMPLICIT NONE
                                                                        
      PRIVATE 
      PUBLIC gf_apws, gf_apws_DIM 
      CONTAINS 
      SUBROUTINE gf_apws(                                               &
     &     jspins,jspin,bk,lpe0,noco,gfinp,cell,lapw,layer)             
                                                                        
!*********************************************************************  
!     determines the lapw list such that |k+G|<rkmax.                   
!     bk(i) is the nk k-point given in internal (i.e. b1,b2,b3) units.  
!     the lapw list is stored on unit29 in direct access mode.          
!     (it is the iteration number )                                     
!        m. weinert  1986                                               
!*********************************************************************  
!     modified for explicit USE of z-reflection symmetry in seclr4.f    
!        g. bihlmayer '96                                               
!     subroutine boxdim added to treat non-orthogonal lattice vectors   
!        s.bluegel, IFF, 18.Nov.97                                      
!*********************************************************************  
!     changed subroutine to generate g-vectors out of a zylinder not    
!     the usual sphere.                                                 
!                                                                       
!     added additional index kp of parallel g-vectors to handle         
!     matrixes expanded in terms of the 2D Basis as used in Embedding   
!     New variables:                                                    
!     kp :index                                                         
!     nv2:no of 2D Basis fcts                                           
!            Daniel Wortmann, Tokyo, 2001                               
!*******************************************************************    
      USE m_gf_types 
      USE m_fleur_interface,ONLY:fleur_boxdim 
      IMPLICIT NONE 
!     ..                                                                
!     .. Scalar Arguments ..                                            
      INTEGER, INTENT  (IN)      :: jspins,jspin 
      LOGICAL , INTENT  (IN)     :: lpe0 
      TYPE(t_lapw),INTENT(INOUT) :: lapw 
      TYPE(t_gfinp),INTENT(IN)   :: gfinp 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_noco),INTENT(IN)    :: noco 
      REAL   ,INTENT(IN)         :: bk(3) 
      INTEGER,INTENT(IN)         :: layer 
                                                                        
!     ..                                                                
!     .. Local Scalars ..                                               
      REAL arltv1,arltv2,arltv3,s3,s2,rk2,rkm,t 
      INTEGER i,itt,j,j1,j2,j3,ki,l,m,mk1,mk2,mk3,n,                    &
     &        ispin,jsp_start,jsp_end,n2                                
      LOGICAL l_sorted 
!     ..                                                                
!     .. Local Arrays ..                                                
      REAL s(3) 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC int,sqrt 
!     ..                                                                
      !initialize
      lapw%nv=0
      lapw%nv2=0
                                                                        
                                                                        
      !<-- Calculate the Box                                            
                                                                        
!---> Determine rkmax box of size mk1, mk2, mk3,                        
!     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax                   
!     arltv(i) length of reciprical lattice vector along direction (i)  
!                                                                       
      !cylinder setup is used if the projected Green fct's are          
      !calculated                                                       
                                                                        
                                                                        
      CALL fleur_boxdim(                                                &
     &            cell%bmat,                                            &
     &            arltv1,arltv2,arltv3)                                 
                                                                        
      mk1 = int( lapw%rkmax/arltv1 ) + 4 
      mk2 = int( lapw%rkmax/arltv2 ) + 4 
                                                                        
      IF (lapw%l_cylinder) THEN 
         mk3 = gfinp%napw(layer) 
      ELSE 
         mk3=int( lapw%rkmax/arltv3 ) + 4 
      ENDIF 
                                                                        
      !>                                                                
      ! No noco-setup in this subroutine !!                             
      jsp_start = jspin 
      jsp_end   = jspin 
      rkm = lapw%rkmax 
                                                                        
      DO ispin = jsp_start,jsp_end 
         rk2 = rkm*rkm 
!                                                                       
!--->    2D obtain vectors                                              
!     First calculate only 2D Vectors, sort them and then add z-komponen
!                                                                       
         n2= 0 
         s(3)=0 
         DO  j1 = -mk1,mk1 
            s(1) = bk(1) + j1 + (2*ispin - 3)/2.0*noco%qss(1) 
            DO  j2 = -mk2,mk2 
               s(2) = bk(2) + j2 + (2*ispin - 3)/2.0*noco%qss(2) 
               s2=dot_product(s,matmul(s,cell%bbmat)) 
!               s2 =dotirp(s,s,bbmat)                                   
               IF (s2<=rk2) THEN 
                  n2=n2+1 
                  lapw%kp%rkp(n2,ispin)=sqrt(s2) 
                  lapw%kp%k1p(n2,ispin)=j1 
                  lapw%kp%k2p(n2,ispin)=j2 
               ENDIF 
            ENDDO 
         ENDDO 
                                                                        
         lapw%nv2(ispin)=n2 
         IF (n2>lapw%nv2d)CALL juDFT_error('gf_apws: wrong 2D-dimensions') 
         IF (lapw%l_cylinder)  lapw%nv(ispin) = n2*(2*mk3+1) 
!                                                                       
! Sort the 2D-PW                                                        
!                                                                       
         l_sorted=.FALSE. 
         DO WHILE(.NOT.l_sorted) 
            l_sorted=.TRUE. 
            DO i=1,n2-1 
               !reorderer such that |G(i)|<|G(i+1)| and                 
               !if |G(i)|=|G(i+1)| and G_x(i)<G_x(i+1) and              
               !if |G(i)|=|G(i+1)| and G_x(i)=G_x(i+1) G_y(i)<G_y(i+1)  
               IF (                                                     &
     &              (lapw%kp%rkp(i,ispin)-lapw%kp%rkp(i+1,ispin)        &
     &              )>EPSILON(1.0).OR.((ABS(lapw%kp%rkp(i,ispin)        &
     &              -lapw%kp%rkp(i+1,ispin))<EPSILON(1.0)).AND.         &
     &               ((lapw%kp%k1p(i,ispin)>lapw%kp%k1p(i+1,ispin)).OR. &
     &               (lapw%kp%k1p(i,ispin) == lapw%kp%k1p(i+1,ispin).AND.&
     &               (lapw%kp%k2p(i,ispin)>lapw%kp%k2p(i+1,ispin))))))  &
     &              THEN                                                
                  t = lapw%kp%rkp(i,ispin) 
                  lapw%kp%rkp(i,ispin) = lapw%kp%rkp(i+1,ispin) 
                  lapw%kp%rkp(i+1,ispin) = t 
                  itt = lapw%kp%k1p(i,ispin) 
                  lapw%kp%k1p(i,ispin)=lapw%kp%k1p(i+1,ispin) 
                  lapw%kp%k1p(i+1,ispin)=itt 
                  itt =lapw%kp%k2p(i,ispin) 
                  lapw%kp%k2p(i,ispin) = lapw%kp%k2p(i+1,ispin) 
                  lapw%kp%k2p(i+1,ispin) = itt 
                  l_sorted=.FALSE. 
               ENDIF 
            ENDDO 
         ENDDO 

         ! now generate a map-to the global g-vectors
         lapw%g2map(:)=0
         do n2=1,lapw%nv2(ispin)
            lapw%g2map(n2)=lapw%global2Dmap(lapw%kp%k1p(n2,ispin),lapw%kp%k2p(n2,ispin))
         enddo
         IF (any(lapw%g2map(:lapw%nv2(ispin))==0)) THEN
            DO n2=1,lapw%nv2(ispin)
                 write(*,*) n2,lapw%g2map(n2),lapw%kp%k1p(n2,ispin),lapw%kp%k2p(n2,ispin)
            ENDDO
            CALL juDFT_error("Mapping failed in gf_apws")
         ENDIF

!                                                                       
!     now add third component of g-vector                               
!                                                                       
         n=0 
         DO n2=1, lapw%nv2(ispin) 
            DO j3 = -mk3,mk3 
               s(1) = bk(1)+lapw%kp%k1p(n2,ispin)+(2*ispin-3)/2.0       &
     &              *noco%qss(1)                                        
               s(2) = bk(2)+lapw%kp%k2p(n2,ispin)+(2*ispin-3)/2.0       &
     &              *noco%qss(2)                                        
               s(3) = bk(3) + j3 + (2*ispin - 3)/2.0*noco%qss(3) 
               s3=dot_product(s,matmul(s,cell%bbmat)) 
               !s3 = dotirp(s,s,bbmat)                                  
               IF ((.NOT.lapw%l_cylinder).AND.s3>rk2) CYCLE 
               n = n + 1 
!     new stop added ! no reduction of gmax anymore!  
               IF (n>lapw%nvd)                                       &
     &             CALL juDFT_error('Too many basis functions !',calledby="gf_apws")
               lapw%k%k1(n,ispin) = lapw%kp%k1p(n2,ispin) 
               lapw%k%k2(n,ispin) = lapw%kp%k2p(n2,ispin) 
               lapw%k%k3(n,ispin) = j3 
               lapw%k%kp(n,ispin) = n2 
               lapw%k%rk(n,ispin) = SQRT(s3) 
            ENDDO 
         ENDDO 
         IF (.NOT.lapw%l_cylinder)  lapw%nv(ispin)=n 
!--->    sort by shell-metzner                                          
         m = n 
   81    m = m/2 
         IF (m<=0) GO TO 131 
         ki = n - m 
         j = 1 
   91    i = j 
  101    l = i + m 
         IF (lapw%k%rk(i,ispin)>lapw%k%rk(l,ispin)) GO TO 121 
  111    j = j + 1 
         IF (j>ki) GO TO 81 
         GO TO 91 
  121    t = lapw%k%rk(i,ispin) 
         lapw%k%rk(i,ispin) = lapw%k%rk(l,ispin) 
         lapw%k%rk(l,ispin) = t 
         itt = lapw%k%kp(i,ispin) 
         lapw%k%kp(i,ispin)=lapw%k%kp(l,ispin) 
         lapw%k%kp(l,ispin)=itt 
         itt = lapw%k%k1(i,ispin) 
         lapw%k%k1(i,ispin) = lapw%k%k1(l,ispin) 
         lapw%k%k1(l,ispin) = itt 
         itt = lapw%k%k2(i,ispin) 
         lapw%k%k2(i,ispin) = lapw%k%k2(l,ispin) 
         lapw%k%k2(l,ispin) = itt 
         itt = lapw%k%k3(i,ispin) 
         lapw%k%k3(i,ispin) = lapw%k%k3(l,ispin) 
         lapw%k%k3(l,ispin) = itt 
         i = i - m 
         IF (i<1) GO TO 111 
         GO TO 101 
  131    CONTINUE 

         !Now check for  largest k+g in sphere
         IF (lapw%l_cylinder) THEN
            DO n=lapw%nv(ispin),1,-1
                IF (lapw%k%rk(n,jspin)<lapw%rkmax) exit
            ENDDO
            lapw%nv_sphere(ispin)=n
            write(6,*) "Basis setup:"
            write(6,*) "Cylinder:,",lapw%nv(ispin)," sphere:",lapw%nv_sphere(ispin)
         ELSE
            lapw%nv_sphere(ispin)=lapw%nv(ispin)
         ENDIF
                                                                        
         IF ((.NOT. noco%l_ss) .AND. (jspins==2) ) THEN 
             lapw%nv(jspins-(jspin-1)) =  lapw%nv(jspin) 
            DO i = 1, lapw%nv(jspin) 
               lapw%k%rk(i,jspins-(jspin-1)) = lapw%k%rk(i,jspin) 
               lapw%k%k1(i,jspins-(jspin-1)) = lapw%k%k1(i,jspin) 
               lapw%k%k2(i,jspins-(jspin-1)) = lapw%k%k2(i,jspin) 
               lapw%k%k3(i,jspins-(jspin-1)) = lapw%k%k3(i,jspin) 
               lapw%k%kp(i,jspins-(jspin-1)) = lapw%k%kp(i,jspin) 
            ENDDO 
             lapw%nv2(jspins-(jspin-1)) =  lapw%nv2(jspin)
             lapw%nv_sphere(jspins-(jspin-1)) =  lapw%nv_sphere(jspin)
             DO i = 1, lapw%nv2(jspin) 
                lapw%kp%rkp(i,jspins-(jspin-1))=lapw%kp%rkp(i,jspin) 
                lapw%kp%k1p(i,jspins-(jspin-1))=lapw%kp%k1p(i,jspin) 
                lapw%kp%k2p(i,jspins-(jspin-1))=lapw%kp%k2p(i,jspin) 
             ENDDO 
         ENDIF 
                                                                        
                                                                        
                                !spin loop                              
      ENDDO 
                                                                        
                                                                        
                                                                        
      !<--Write list of 2D-gvectors to out                              
                                 !Do not write this anymore             
      IF (lpe0.AND..FALSE.) THEN 
         WRITE (6,*) 'G_p vectors:' 
         DO i = 1, lapw%nv2(jspin) 
            WRITE (6,*) i, lapw%kp%k1p(i,jspin),lapw%kp%k2p(i,jspin) 
         ENDDO 
      ENDIF 
      !>                                                                
                                                                        
      !<--Set dimensions correctly for noco                             
      lapw%nmat = lapw%nv(jspin)
      lapw%nmat_sphere=  lapw%nv_sphere(jspin)
      IF (noco%L_noco) lapw%nmat = lapw%nv(1)+lapw%nv(2)
      IF (noco%L_noco) lapw%nmat_sphere = lapw%nv_sphere(1)+lapw%nv_sphere(2)
                                                                        
      IF (noco%l_noco) THEN 
         lapw%nv_tot = lapw%nv(1)+lapw%nv(2) 
         lapw%nv_tot_sphere = lapw%nv_sphere(1)+lapw%nv_sphere(2)
         lapw%nv2_tot = lapw%nv2(1)+lapw%nv2(2) 
      ELSE 
         lapw%nv_tot = lapw%nv(jspin) 
         lapw%nv_tot_sphere = lapw%nv_sphere(jspin)
         lapw%nv2_tot = lapw%nv2(jspin) 
      ENDIF 
      !>                                                                
                                                                        
                                                                        
      RETURN 
      END SUBROUTINE gf_apws 
      !<-- S:gf_apws_dim(qss,kpts,gfinp,cell,lapw)                      
                                                                        
      SUBROUTINE gf_apws_DIM(jspins,                                    &
     &     qss,kpts,gfinp,cell,rmtlmax,                                 &
     &     lapw,layer)                                                  
!******************************************                             
!     Calculates the max-no of lapw's                                   
!     D. Wortmann                                                       
!******************************************                             
      USE m_gf_types 
      USE m_fleur_interface,ONLY:fleur_boxdim 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(IN)            :: jspins 
      REAL   ,INTENT(IN)            :: qss(3),rmtlmax(2) 
      TYPE(t_kpts),INTENT(IN)       :: kpts 
      TYPE(t_cell),INTENT(IN)       :: cell 
      TYPE(t_gfinp),INTENT(INOUT)   :: gfinp 
      TYPE(t_lapw),INTENT(INOUT)    :: lapw 
      INTEGER,INTENT(IN)            :: layer 
                                                                        
      !>                                                                
      !<--Locals                                                        
                                                                        
      INTEGER :: mk1,mk2,mk3 
      REAL    :: s(3),rk2,s2,s3 
      INTEGER :: ispin,nk,j1,j2,j3,n,n2,nn 
      INTEGER,ALLOCATABLE :: k1p(:),k2p(:) 
      REAL                :: gkmax,gkmax2 
      !>                                                                
                                                                        
                                                                        
      gfinp%napw(layer)=ABS(gfinp%napw(layer)) 
                                                                        
      !<--Determine maximal boxsize!                                    
                                                                        
      IF (gfinp%napw(layer)/=0) THEN 
         lapw%l_cylinder=.TRUE. 
      ELSE 
         lapw%l_cylinder=.FALSE. 
      ENDIF 
      CALL fleur_boxdim(                                                &
     &     cell%bmat,                                                   &
     &     s(1),s(2),s(3))                                              
      mk1 = int( lapw%rkmax/s(1) ) + 4 
      mk2 = int( lapw%rkmax/s(2) ) + 4 
      IF (lapw%l_cylinder) THEN 
         mk3 = gfinp%napw(layer) 
      ELSE 
         mk3=INT( lapw%rkmax/s(3) ) + 4 
      ENDIF 
      ALLOCATE(k1p((2*mk1+1)*(2*mk2+1)),k2p((2*mk1+1)*(2*mk2+1))) 
                                                                        
      !>                                                                
                                                                        
      !<--loop over all kpts and spins                                  
                                                                        
      lapw%nvd=0 
      lapw%nv2d=0 
      gkmax2 = 0.0 
      gkmax = 0.0 
                   !Only needed for spin-spiral                         
      DO ispin=1,2 
         DO nk=1,kpts%nkpts 
            !<--Calculate 2D-Basis functions!                           
                                                                        
            rk2 = lapw%rkmax*lapw%rkmax 
            n2= 0 
            s(3)=0 
            DO  j1 = -mk1,mk1 
               s(1) = kpts%bk(1,nk) + j1 + (2*ispin - 3)/2.0*qss(1) 
               DO  j2 = -mk2,mk2 
                  s(2) = kpts%bk(2,nk) + j2 + (2*ispin - 3)/2.0*qss(2) 
                  s2=dot_product(s,matmul(s,cell%bbmat)) 
!                  s2 =dotirp(s,s,cell%bbmat)                           
                  IF (s2<=rk2) THEN 
                     gkmax2=max(gkmax2,s2) 
                     n2=n2+1 
                     k1p(n2)=j1 
                     k2p(n2)=j2 
                  ENDIF 
               ENDDO 
            ENDDO 
            lapw%nv2d=max(n2,lapw%nv2d) 
                                                                        
            !>                                                          
            IF (lapw%l_cylinder) THEN 
               lapw%nvd=max(lapw%nvd,lapw%nv2d*(2*gfinp%napw(layer)+1)) 
            ELSE 
               !<--Calculate the 3d-Basis                               
                                                                        
               n=0 
               DO nn=1, n2 
                  DO j3 = -mk3,mk3 
                     s(1) = kpts%bk(1,nk)+k1p(nn)+(2*ispin-3)/2.*qss(1) 
                     s(2) = kpts%bk(2,nk)+k2p(nn)+(2*ispin-3)/2.*qss(2) 
                     s(3) = kpts%bk(3,nk)+j3 +(2*ispin-3)/2. *qss(3) 
                     s3 = dot_PRODUCT(s,MATMUL(s,cell%bbmat)) 
!                     s3 = dotirp(s,s,cell%bbmat)                       
                     IF (s3>rk2) CYCLE 
                     gkmax=max(s3,gkmax) 
                                                                 !Set na
                     gfinp%napw(layer)=max(j3,gfinp%napw(layer)) 
                                                 !maximum if in sc-mode 
                     n = n + 1 
                  ENDDO 
               ENDDO 
               lapw%nvd=max(n,lapw%nvd) 
                                                                        
               !>                                                       
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
                                                                        
                                                                        
      !>                                                                
      !>                                                                
      DEALLOCATE(k1p,k2p) 
      WRITE(6,*) "Basis setup:" 
      WRITE(6,*) "Maximal dimensions:" 
      WRITE(6,*) "nvd:",lapw%nvd 
      WRITE(6,*) "nv2d:",lapw%nv2d 
      WRITE(6,*) "Matching at MT:" 
      IF (gkmax == 0.0)                                                 &
     &     gkmax = sqrt(gkmax2**2+(cell%bmat(3,3)*gfinp%napw(layer))**2)
      WRITE(6,*) "Maximal |k+g|=",gkmax 
      WRITE(6,*) "Maximal r*l  =",rmtlmax(2) 
      WRITE(6,*) "Minimal r*l  =",rmtlmax(1) 
                                                                        
      !<-- Allocate storage for lapw's                                  
      !Basis                                                            
      ALLOCATE(lapw%k%rk(lapw%nvd,jspins));lapw%k%rk = 0. 
      ALLOCATE(lapw%k%k1(lapw%nvd,jspins));lapw%k%k1= 0 
      ALLOCATE(lapw%k%k2(lapw%nvd,jspins));lapw%k%k2= 0 
      ALLOCATE(lapw%k%k3(lapw%nvd,jspins));lapw%k%k3= 0 
      ALLOCATE(lapw%k%kp(lapw%nvd,jspins));lapw%k%kp= 0 
      ALLOCATE(lapw%kp%rkp(lapw%nv2d,jspins));lapw%kp%rkp= 0. 
      ALLOCATE(lapw%kp%k1p(lapw%nv2d,jspins));lapw%kp%k1p= 0 
      ALLOCATE(lapw%kp%k2p(lapw%nv2d,jspins));lapw%kp%k2p=0 
      !>                                                                
                                                                        
      CALL gf_apws_globalmap(cell,lapw) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_apws_globalmap(cell,lapw)                               
                                                                        
      SUBROUTINE gf_apws_globalmap(cell,                                &
     &     lapw)                                                        
!******************************************                             
!   Creates the mapping list of the g-vectors needed to store the       
!   embedding potentials                                                
!******************************************                             
      USE m_gf_types 
      USE m_fleur_interface,ONLY:fleur_boxdim 
      USE m_gf_math,ONLY:sort 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_cell),INTENT(IN)       :: cell 
      TYPE(t_lapw),INTENT(INOUT)    :: lapw 
      !>                                                                
      !<--Locals                                                        
                                                                        
      INTEGER :: mk1,mk2 
      REAL    :: s(3),rk2,s_min 
      INTEGER :: j1,j2,n,n2 
      INTEGER,ALLOCATABLE :: k1p(:),k2p(:),sindex(:) 
      REAL   ,ALLOCATABLE :: kl(:) 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- Determine maximal boxsize!                                   
                                                                        
      CALL fleur_boxdim(                                                &
     &     cell%bmat,                                                   &
     &     s(1),s(2),s(3))                                              
      mk1 = int( lapw%rkmax/s(1) ) + 4 
      mk2 = int( lapw%rkmax/s(2) ) + 4 
                                                                        
      ALLOCATE(sindex((2*mk1+1)*(2*mk2+1)),kl((2*mk1+1)*(2*mk2+1))      &
     &     ,k1p((2*mk1+1)*(2*mk2+1)),k2p((2*mk1+1)*(2*mk2+1)))          
                                                                        
      !>                                                                
                                                                        
      !<-- Calculate 2D-g vectors!                                      
      rk2 = lapw%rkmax*lapw%rkmax 
      n2 = 0 
      s(3)= 0 
      DO  j1 = -mk1,mk1 
         DO  j2 = -mk2,mk2 
            !Try all combinations with maximal k-vectors                
            s(1:2) = (/j1,j2/) 
            s_min  = dot_PRODUCT(s,MATMUL(s,cell%bbmat)) 
            s(1:2) = (/j1+0.5,j2+0.0/) 
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat))) 
            s(1:2) = (/j1+0.5,j2+0.5/) 
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat))) 
            s(1:2) = (/j1+0.0,j2+0.5/) 
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat))) 
            s(1:2) = (/j1-0.5,j2+0.0/) 
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat))) 
            s(1:2) = (/j1-0.5,j2-0.5/) 
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat))) 
            s(1:2) = (/j1+0.0,j2-0.5/) 
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat))) 
            s(1:2) = (/j1+0.5,j2-0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1-0.5,j2+0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            IF (s_min <= rk2) THEN 
               !We found a g-vector                                     
               n2 = n2+1 
               k1p(n2) = j1 
               k2p(n2) = j2 
               kl(n2) = sqrt(1.0*(j1**2+j2**2)) 
            ENDIF 
         ENDDO 
      ENDDO 
      !>                                                                
      !<-- Sort and create map (use sort from gf_tools)                 
      CALL sort(kl(:n2),sindex(:n2)) 
                                                                        
      ALLOCATE(lapw%global2Dlist(2,n2)) 
      ALLOCATE(lapw%g2map(n2))

      lapw%global2Dlist(1,:) = k1p(sindex(:n2)) 
      lapw%global2Dlist(2,:) = k2p(sindex(:n2)) 
      ALLOCATE(lapw%global2Dmap(-mk1:mk1,-mk2:mk2)) 
      lapw%global2Dmap = 0 
      DO n = 1,n2 
         lapw%global2Dmap(lapw%global2Dlist(1,n),lapw%global2Dlist(2,n))&
     &        = n                                                       
      ENDDO 
      !>                                                                
                                                                        
      DO N = 1,N2 
         WRITE(6,*) N,lapw%global2Dlist(:,n) 
         IF (lapw%global2Dmap(lapw%global2Dlist(1,n),lapw%global2Dlist(2&
     &        ,n)) /= n) CALL juDFT_error                                 &
     &        ("Internal error in gf_apws_globalmap")                   
      ENDDO 
      lapw%g_max = (/MAXVAL(lapw%global2Dlist(1,:))                     &
     &     ,MAXVAL(lapw%global2Dlist(2,:))/)                            
      DEALLOCATE(k1p,k2p,kl,sindex) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
