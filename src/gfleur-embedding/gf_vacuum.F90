!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_vacuum 
      use m_juDFT
      IMPLICIT NONE
      PRIVATE 
#ifdef oldstuff
                                                !matrix for coordinate  
      REAL,SAVE               :: coordconv(2,2) 
                                                !transformation, i.e. bm
                                                !kp-vectors             
      INTEGER,POINTER,SAVE    :: k1p(:),k2p(:) 
                                                !k-vectors              
      REAL,Allocatable,SAVE       :: bk(:,:) 
                                                !size of aux-region     
      REAL,SAVE               :: d_aux 
                                                !d of cell              
      REAL,SAVE               :: d 
                                                !vacuum level           
      REAL,SAVE               :: vac_pot 
                                                !position of image plane
      REAL,SAVE               :: imageplane 
                                                !position of image bound
      REAL,SAVE               :: imageboundary 
                                                !EField applied to surfa
      REAL,SAVE               :: EField 
                                                     !do we have a surfa
      LOGICAL,SAVE            :: surface = .FALSE.,l_noco
                                                                        
                                                 !dt of cell            
      REAL,SAVE               :: dt 
                                                 !FFT box in  z directio
      INTEGER,SAVE            :: mx3 
                                                                        
                                                     !do we have an imag
      LOGICAL,SAVE            :: imagepoti = .FALSE. 
                                                                        
                                                  ! #gridpoints of mesh 
      INTEGER,PARAMETER       :: gridpts=1000000 
                                             !stores V(z)               
      REAL, SAVE              :: V(gridpts) 
                                            ! stores x pos. on mesh     
      REAL, SAVE              :: x(gridpts) 
                                                    !mesh spacing       
      REAL, SAVE              :: xdistance,xdist_sq 
                                        !If V=const set                 
      REAL, SAVE              :: V0=0 
                                                                        
                                                                        
      PUBLIC                  :: gf_vacuum_embpot,gf_vacuum_init 
      PUBLIC                  :: gf_vacuum_Poti_init 
!     PUBLIC                  :: gf_vacuum,gf_vacuum_image_embpot       
      PUBLIC                  :: gf_vacuum 
      PUBLIC                  :: gridpts 
      public                  :: gf_vacuum_check_vacpot
#endif
      public gf_simple_vac
                                                                        
      CONTAINS 

      subroutine gf_simple_vac(en,nk,lapw,kpts,cell,vacuum_energy,sigma,efield)
      USE m_gf_types
      USE m_gf_energies
      integer,intent(in)       :: en,nk
      type(t_lapw),intent(in)  :: lapw
      type(t_kpts),intent(in)  :: kpts
      type(t_cell),intent(in)  :: cell
      real,intent(in)          :: vacuum_energy
      complex,intent(out)      :: sigma(:,:)
      real,intent(in),optional :: efield

      integer :: nn,n
      real    :: gk(2),s


      nn=size(sigma,1)
      if (lapw%nv2(1)<nn) nn=nn/2 !noco calculation
      if (present(efield)) THEN
         if (abs(efield)>epsilon(0.0)) CALL vacuum_EField(en,nk,lapw,kpts,cell,vacuum_energy,sigma,efield)
      ELSE
        sigma=cmplx(0.0,0.0)
        DO n=1,nn
          gk=MATMUL(cell%bmat(1:2,1:2),kpts%bk(1:2,nk)+(/lapw%kp%k1p(n,1),lapw%kp%k2p(n,1)/))
          s=dot_PRODUCT(gk,gk)
          sigma(n,n)=cmplx(0,0.5)* SQRT(2*(gf_Z(en,0)-vacuum_energy)-s)
        enddo
      ENDIF
      if (nn.ne.size(sigma,1)) then !noco
        sigma(nn+1:,nn+1:)=sigma(:nn,:nn)
      endif

      end subroutine


       SUBROUTINE vacuum_EField(en,nk,lapw,kpts,cell,vacuum_energy,sigma,efield)
!******************************************
! calculates numericaly the embedding potential for
! linear E-field E*z+V_vac applyed to the system
!       01/2007            A. Hanuschkin
!       modified by DW 2012
!******************************************
      USE m_gf_types
      USE m_gf_energies
      IMPLICIT NONE
      integer,intent(in)       :: en,nk
      type(t_lapw),intent(in)  :: lapw
      type(t_kpts),intent(in)  :: kpts
      type(t_cell),intent(in)  :: cell
      real,intent(in)          :: vacuum_energy
      complex,intent(out)      :: sigma(:,:)
      real,intent(in),optional :: efield
      !-Locals
      REAL    :: gk(2),s
      COMPLEX :: dphi
      INTEGER :: i,n,nn
      INTEGER,PARAMETER       :: gridpts=100000
      REAL                    :: V(gridpts),x(gridpts)
      COMPLEX                 :: phi(gridpts),ksq(gridpts)
      REAL, SAVE              :: xdistance,xdist_sq !mesh spacing


      !Generate potential in vacuum
      xdistance = REAL(100)/REAL(gridpts)
      xdist_sq = xdistance*xdistance

      DO i=1,gridpts
         x(i)=REAL(i-1)*xdistance;
         V(i) = vacuum_energy+(EField*x(i))
      END DO

      !init Embedding potential
      sigma=0.0

      !CALCULATE PHI and DPHI
      nn=size(sigma,1)
      if (lapw%nv2(1)<nn) nn=nn/2 !noco calculation

      DO n=1,nn
           !initial wavefunction
           !inf. in vacuum region V=vac_pot
           gk=MATMUL(cell%bmat(:2,:2),kpts%bk(1:2,nk)+(/lapw%kp%k1p(n,1),lapw%kp%k2p(n,1)/))
           s=dot_PRODUCT(gk,gk)
           !numerical integration ksq=k**2xh**2
                          !set up k^2*h^2 array
           DO i=1,gridpts
               ksq(i) = ( 2.0*(gf_Z(en,0)-V(i))-s)*xdist_sq
           ENDDO
           phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts) )-s)*x(gridpts))
           phi(gridpts-1)= exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts-1) )-s)*x(gridpts-1))
           DO i=gridpts-2,1,-1
             ! normalize!
              IF (ABS(phi(i+1))>1000.) THEN
                 phi(i+1)=phi(i+1)/ABS(phi(i+1))
                 phi(i+2)=phi(i+2)/ABS(phi(i+1))
              ENDIF
              !calculate phi
               phi(i)=((2.-(10./12.)*ksq(i+1))*phi(i+1)-(1.+(1./12.)*ksq(i+2))*phi(i+2))/(1.+(1./12.)*ksq(i))
           ENDDO
           dphi = -0.5 * ((phi(1)-phi(1+1))/(xdistance)+((phi(1+1)-phi(1+2)))/(xdistance)  )
           sigma(n,n)= 0.5*dphi/phi(1)
      ENDDO


      END SUBROUTINE

      !>
#ifdef old_stuff
      !<-- S: gf_vacuum_embpot(en,nk,g2)                                

      SUBROUTINE gf_vacuum_embpot(en,nk,embpot) 
!******************************************                             
!     Construct the embedding potential for a surface calculation.      
!     The vacuum potential must be set before by calling gf_vacuum_init 
!     Please note that the potential in the aux-volume might not be zero
!     See Eq. 5.50 from my Thesis!                                      
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_energies 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en,nk 
      COMPLEX,INTENT(INOUT) ::embpot(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n,nn
      REAL    :: gk(2),s 
      COMPLEX :: kl,kr,I 
      I=cmplx(0.0,1.0) 
      !>                                                                
                                                                        
                               !no surface calculation                  
      IF (.NOT.surface) RETURN 
                          !cal with imagepoti                           
      IF (imagepoti) THEN 
          CALL gf_vacuum_image_embpot(en,nk,embpot) 
          RETURN 
      ENDIF 
                              !cal with E-field                         
      IF (EField/=0.0) THEN 
          CALL gf_vacuum_EField(en,nk,embpot) 
          RETURN 
      ENDIF 
                                                                        
                                                                        
                                                                        
                               !init Embedding potential                
      embpot=0.0 

      !calculate embedding potential
      ! EmbPot vacuum V=0
      nn=size(embpot,1)
	  if (l_noco) nn=nn/2
      DO n=1,nn
          gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/))
          s=dot_PRODUCT(gk,gk)
          embpot(n,n)=I/2* SQRT(2*(gf_Z(en,0)-vac_pot)-s)
      enddo
      if (l_noco) then
      	embpot(nn+1:,nn+1:)=embpot(:nn,:nn)
      endif
!      DO n=1,SIZE(embpot,1)
!         gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/))
!         s=dot_PRODUCT(gk,gk)
!         kl=SQRT(2*gf_Z(en,0)-s)
!         kr=SQRT(2*(gf_Z(en,0)-vac_pot)-s)
!         embpot(n,n)=I/2*kl*(kr*COS(kl*d_aux)+I*kl*SIN(kl*d_aux))/(kl   &
!     &        *COS(kl*d_aux)+I*kr*SIN(kl*d_aux))
!      WRITE (0,*) en, REAL(gf_Z(en)),nk,"S:",embpot(n,n)               
!     &    ,gk,s,REAL(gf_Z(en)),"k:",REAL(kl)                           
!     & ,REAL(kr)                                                       
                                                                        
!       WRITE (0,*) en, REAL(gf_Z(en)),embpot(n,n)                      
!      ENDDO
! "(5i,f8.4,5i,2f12.5,2f8.4,2f8.4)")                                    
                                                                        
      END SUBROUTINE 

      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum_image_embpot(en,nk,g2)                          

      SUBROUTINE gf_vacuum_image_embpot(en,nk,embpot) 
!******************************************                             
!       12/2006            A. Hanuschkin                                
!******************************************                             
      USE m_gf_energies 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en,nk 
      COMPLEX,INTENT(INOUT) ::embpot(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n 
      REAL    :: gk(2),s 
      COMPLEX :: k_r,I 
                                              ! k**2                    
      COMPLEX, ALLOCATABLE :: ksq(:) 
                                              ! derivation of wf        
      COMPLEX :: dphi 
                                              ! wavefunction on mesh    
      COMPLEX , ALLOCATABLE :: phi(:) 
                                              ! save comp. time         
      REAL    :: one12 = (1/12.0) 
      INTEGER :: ii,index_c,zz 
      CHARACTER (LEN=20) :: filename 
                                                                        
      I=cmplx(0.0,1.0) 
      !>                                                                
                                                                        
!     if (.not.surface) return !no surface calculation                  
      !CHECK THAT IMAGEBOUNDARY >= D/2                                  
      !CHECK THAT IMAGEPLANE  < IMAGEBOUNDARY                           
      !CHECK imageplane < embPlane @ reading gf_inp                     
                                                                        
                                                                        
      !allocate and initialize                                          
      ALLOCATE(ksq(gridpts)) 
      ALLOCATE(phi(gridpts)); 
                               !init Embedding potential                
      embpot=0.0 
                                                                        
#ifdef CPP_DEBUG   
                   !not nessecary because it'll be calculated in necessa
                   !but nice if plotted phi(x<c/2)                      
      phi = 0 
                                                                        
      WRITE (*,*) "xdistance,imageboundary,imageplane"                  &
     &            ,xdistance,imageboundary,imageplane                   
#endif                                                                  
                                                                        
                                                                        
       !phi and dphi at z= c/2 =d-d_aux needed                          
       !up to which point integration is needed?                        
                                                        !check imageplan
       index_c = Int((d-d_aux-imageplane)/(xdistance)) 
!       WRITE (*,*) index_c," index_c ",d-d_aux-imageplane              
!     &         ," d-d_aux-imageplane ",d-d_aux," c"                    
                                                                        
       IF ((index_c>gridpts-2).OR.(index_c<=0)) THEN 
         WRITE (*,*) index_c 
         WRITE (*,*) "Imageboundary closer than Embedding plane?" 
          CALL juDFT_error("PROBLEMS WITH SPACING: INDEX_C",calledby="gf_vacuum.F90")
       ENDIF 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      !CALCULATE PHI and DPHI at c/2                                    
      DO n=1,SIZE(embpot,1) 
      !initial wavefunction                                             
      !inf. in vacuum region V=vac_pot                                  
                                                                        
       gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/)) 
       s=dot_PRODUCT(gk,gk) 
                                                                        
                                                                        
#ifdef CPP_DEBUG                                                        
       WRITE (6,*) "s ",s," vac_pot ",vac_pot 
#endif                                                                  
                                                                        
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-vac_pot)-s)            &
     &                    *x(gridpts))                                  
       phi(gridpts-1)= exp(I*sqrt(2.0*(gf_Z(en,0)-vac_pot)-s)           &
     &                    *x(gridpts-1))                                
       ksq(gridpts) = (2.0*(gf_Z(en,0)-vac_pot)-s)*xdist_sq 
       ksq(gridpts-1)= ksq(gridpts) 
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
! WITHOUT POTENTIAL                                                     
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)                 &
     &                    *x(gridpts))                                  
       phi(gridpts-1) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)               &
     &                    *x(gridpts-1))                                
      ksq(gridpts) = (2.0*(gf_Z(en,0))-s)*xdist_sq 
      ksq(gridpts-1)= ksq(gridpts) 
#endif                                                                  
                                                                        
                                                                        
       !numerical integration ksq=k**2xh**2                             
                          !set up k^2*h^2 array                         
       DO ii=1,gridpts-2 
        ksq(ii) = ( 2.0*(gf_Z(en,0)-V(ii)) -S)*xdist_sq 
       ENDDO 
                                                                        
                                   !calculate phi                       
       DO ii=gridpts-2,index_c,-1 
         phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)                 &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
!         IF (ABS(REAL(phi(ii))).gt.1000) THEN  ! normalize!            
!           phi(ii+1)=phi(ii+1)/ABS(Real(phi(ii)))                      
!           phi(ii+2)=phi(ii+2)/ABS(Real(phi(ii)))  !index_x+2 ii+3?    
                                                                        
                                         ! normalize!                   
         IF (ABS(phi(ii))>1000) THEN 
           phi(ii+1)=phi(ii+1)/ABS(phi(ii)) 
                                             !index_x+2 ii+3?           
           phi(ii+2)=phi(ii+2)/ABS(phi(ii)) 
                                                                        
                                                                        
#ifdef CPP_DEBUG_IP_NA                                                  
                                             ! normalize all!           
            DO zz=gridpts,ii+3,-1 
!           phi(zz)=phi(zz)/ABS(Real(phi(ii)))                          
             phi(zz)=phi(zz)/ABS(phi(ii)) 
            END DO 
#endif                                                                  
            phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)              &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
         ENDIF 
       ENDDO 
                                                                        
       !get phi and dphi at z= c/2 =d-d_aux                             
       !calculate simple derivation of dphi with finite differences     
                                                                        
       !phi(index_c) and                                                
       dphi = -0.5 * ((phi(index_c)-phi(index_c+1))/(xdistance)         &
     &               +((phi(index_c+1)-phi(index_c+2)))/(xdistance)  )  
                                                                        
                                                                        
                                                                        
#ifdef CPP_DEBUG_AP    
!write EmbPot to file                            
!       dphi = 0.5 * ( (phi(1)-phi(2))/(xdistance)                      
!     &               +(phi(2)-phi(3))/(xdistance)  )                   
                                                                        
       OPEN(109,FILE="EmbPot",FORM='formatted',                         &
     &         STATUS='OLD',ACCESS='append')                            
!      WRITE (109,'(A)') "# data Re(E),Im(E):Phi:dPhi:Re(Emb):Im(Emb)"  
       WRITE (109,'(8F10.5)')  gf_Z(en,0),phi(index_c),dphi             &
     &                         ,(dphi/phi(index_c))                     
!       WRITE (109,'(7F11.5)')  gf_Z(en,0),REAL(phi(1)),AIMAG(phi(1))   
!     &                        ,REAL(dphi),AIMAG(dphi)                  
!     &                        ,REAL(dphi/phi(1)),AIMAG(dphi/phi(1))    
       CLOSE(109) 
#endif                                                                  
                                                                        
#ifdef CPP_DEBUG_IP    
!write wavefunction to file                      
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,".dat"                                            
        IF (n/=1) CLOSE (66) 
        OPEN(66,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
                                                                        
       WRITE (66,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(phi), Im(phi), V, S)"     
                        !output phi                                     
       DO ii=1,gridpts 
         WRITE (66,'(5F11.5,F7.4)') (imageplane+x(ii)),                 &
     &             (imageplane+x(ii))*(mx3-1)/dt,                       &
     &             phi(ii),V(ii),S                                      
       ENDDO 
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,"-dphi.dat"                                       
        IF (n/=1) CLOSE (77) 
        OPEN(77,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
       WRITE (77,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(dphi), Im(dphi), V, S)"   
         WRITE (77,'(5F11.5,F7.4)') (imageplane+x(ii)),                 &
     &             (imageplane+x(ii))*(mx3-1)/dt,                       &
     &             dphi,V(ii),S                                         
!       WRITE (*,*) "d_aux ",d_aux                                      
#endif                                                                  
                                                                        
                                                                        
                            !saves in case d_aux==0 some cpu time       
       IF (d_aux > 0 ) THEN 
        !move trough delta region  & store EmbPoti                      
                ! phi(index_c) , dphi at left side                      
        k_r = sqrt(2.0*gf_Z(en,0)-s) 
        embpot(n,n)= I/2*k_r*                                           &
     &                (                                                 &
     &                         (dphi/(I*k_r)* COS(k_r*d_aux)            &
     &                         +I*phi(index_c)*SIN(k_r*d_aux))          &
     &                  /                                               &
     &                         ( phi(index_c)*COS(k_r*d_aux)            &
     &                        +I*(dphi/(I*k_r))*SIN(k_r*d_aux) )        &
     &                 )                                                
      ELSE 
        embpot(n,n)= 0.5*dphi/phi(index_c) 
      ENDIF 
                                                                        
!      WRITE (0,*) en, REAL(gf_Z(en,0)),embpot(n,n)                     
             ! n                                                        
      ENDDO 
                                                                        
#ifdef CPP_DEBUG_IP    
!close wf. file                                  
      CLOSE (66) 
      CLOSE (77) 
#endif                                                                  
      RETURN 
                                                                        
      END SUBROUTINE 

      !>                                                                
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum_EField(en,nk,g2)                                

      SUBROUTINE gf_vacuum_EField(en,nk,embpot) 
!******************************************                             
! calculates numericaly the embedding potential for                     
! linear E-field E*z+V_vac applyed to the system                        
!       01/2007            A. Hanuschkin                                
!******************************************                             
      USE m_gf_energies 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en,nk 
      COMPLEX,INTENT(INOUT) ::embpot(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER :: n 
      REAL    :: gk(2),s 
      COMPLEX :: k_r,I 
                                              ! k**2                    
      COMPLEX, ALLOCATABLE :: ksq(:) 
                                              ! derivation of wf        
      COMPLEX :: dphi 
                                              ! wavefunction on mesh    
      COMPLEX , ALLOCATABLE :: phi(:) 
                                              ! save comp. time         
      REAL    :: one12 = (1/12.0) 
      INTEGER :: ii,zz 
      CHARACTER (LEN=20) :: filename 
                                                                        
      I=cmplx(0.0,1.0) 
      !>                                                                
                                                                        
!      if (.not.surface) return !no surface calculation                 
      !CHECK THAT IMAGEBOUNDARY >= D/2                                  
      !CHECK THAT IMAGEPLANE  < IMAGEBOUNDARY                           
      !CHECK imageplane < embPlane @ reading gf_inp                     
                                                                        
                                                                        
      !allocate and initialize                                          
      ALLOCATE(ksq(gridpts)) 
      ALLOCATE(phi(gridpts)); 
                               !init Embedding potential                
      embpot=0.0 
                                                                        
#ifdef CPP_DEBUG   
!not nessecary because it'll be calculated in necessa
                   !but nice if plotted phi(x<c/2)                      
      phi = 0 
      WRITE (*,*) "xdistance,imageboundary,imageplane"                  &
     &            ,xdistance,imageboundary,imageplane                   
#endif                                                                  
                                                                        
       !phi and dphi at z= c/2 =d-d_aux needed                          
       !up to which point integration is needed?                        
!       index_c = Int((d-d_aux)/(xdistance))  !check imageplane < embPla
!       WRITE (*,*) index_c," index_c ",d-d_aux-imageplane              
!     &         ," d-d_aux-imageplane ",d-d_aux," c"                    
                                                                        
!       IF ((index_c.gt.gridpts-2).or.(index_c.le.0)) THEN              
!         WRITE (*,*) index_c                                           
!        WRITE (*,*) "Imageboundary closer than Embedding plane?"       
!          CALL juDFT_error("PROBLEMS WITH SPACING: INDEX_C",calledby="gf_vacuum.F90")
!       ENDIF                                                           
                                                                        
                                                                        
      !CALCULATE PHI and DPHI at c/2                                    
      DO n=1,SIZE(embpot,1) 
      !initial wavefunction                                             
      !inf. in vacuum region V=vac_pot                                  
                                                                        
       gk=MATMUL(coordconv,bk(1:2,nk)+(/k1p(n),k2p(n)/)) 
       s=dot_PRODUCT(gk,gk) 
                                                                        
#ifdef CPP_DEBUG                                                        
       WRITE (6,*) "s ",s," vac_pot ",vac_pot 
#endif                                                                  
                                                                        
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts) )-s)        &
     &                    *x(gridpts))                                  
       phi(gridpts-1)= exp(I*sqrt(2.0*(gf_Z(en,0)-V(gridpts-1) )-s)     &
     &                    *x(gridpts-1))                                
       ksq(gridpts) = (2.0*(gf_Z(en,0)-V(gridpts))-s)*xdist_sq 
       ksq(gridpts-1)= (2.0*(gf_Z(en,0)-V(gridpts-1))-s)*xdist_sq 
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
! WITHOUT POTENTIAL                                                     
       phi(gridpts) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)                 &
     &                    *x(gridpts))                                  
       phi(gridpts-1) = exp(I*sqrt(2.0*(gf_Z(en,0)-V0)-s)               &
     &                    *x(gridpts-1))                                
      ksq(gridpts) = (2.0*(gf_Z(en,0))-s)*xdist_sq 
      ksq(gridpts-1)= ksq(gridpts) 
#endif                                                                  
                                                                        
                                                                        
       !numerical integration ksq=k**2xh**2                             
                          !set up k^2*h^2 array                         
       DO ii=1,gridpts-2 
        ksq(ii) = ( 2.0*(gf_Z(en,0)-V(ii)) -S)*xdist_sq 
       ENDDO 
                                                                        
                             !calculate phi                             
       DO ii=gridpts-2,1,-1 
         phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)                 &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
!         IF (ABS(REAL(phi(ii))).gt.1000) THEN  ! normalize!            
!           phi(ii+1)=phi(ii+1)/ABS(Real(phi(ii)))                      
!           phi(ii+2)=phi(ii+2)/ABS(Real(phi(ii)))  !index_x+2 ii+3?    
                                                                        
                                         ! normalize!                   
         IF (ABS(phi(ii))>1000) THEN 
           phi(ii+1)=phi(ii+1)/ABS(phi(ii)) 
                                             !index_x+2 ii+3?           
           phi(ii+2)=phi(ii+2)/ABS(phi(ii)) 
                                                                        
                                                                        
#ifdef CPP_DEBUG_IP_NA                                                  
                                             ! normalize all!           
            DO zz=gridpts,ii+3,-1 
!           phi(zz)=phi(zz)/ABS(Real(phi(ii)))                          
             phi(zz)=phi(zz)/ABS(phi(ii)) 
            END DO 
#endif                                                                  
            phi(ii)=( (2-(10.0*one12)*ksq(ii+1))*phi(ii+1)              &
     &             -(1+(1*one12)*ksq(ii+2))*phi(ii+2)  )                &
     &               /(1+(1*one12)*ksq(ii))                             
         ENDIF 
       ENDDO 
                                                                        
       !get phi and dphi at z= c/2 =d-d_aux                             
       !calculate simple derivation of dphi with finite differences     
                                                                        
       !phi(index_c) and                                                
       dphi = -0.5 * ((phi(1)-phi(1+1))/(xdistance)                     &
     &               +((phi(1+1)-phi(1+2)))/(xdistance)  )              
                                                                        
                                                                        
                                                                        
#ifdef CPP_DEBUG_AP
    !write EmbPot to file                            
!       dphi = 0.5 * ( (phi(1)-phi(2))/(xdistance)                      
!     &               +(phi(2)-phi(3))/(xdistance)  )                   
                                                                        
       OPEN(109,FILE="EmbPot",FORM='formatted',                         &
     &         STATUS='OLD',ACCESS='append')                            
!      WRITE (109,'(A)') "# data Re(E),Im(E):Phi:dPhi:Re(Emb):Im(Emb)"  
       WRITE (109,'(8F10.5)')  gf_Z(en,0),phi(1),dphi                   &
     &                         ,(dphi/phi(1))                           
!       WRITE (109,'(7F11.5)')  gf_Z(en,0),REAL(phi(1)),AIMAG(phi(1))   
!     &                        ,REAL(dphi),AIMAG(dphi)                  
!     &                        ,REAL(dphi/phi(1)),AIMAG(dphi/phi(1))    
       CLOSE(109) 
#endif                                                                  
                                                                        
#ifdef CPP_DEBUG_IP   
 !write wavefunction to file                      
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,".dat"                                            
        IF (n/=1) CLOSE (66) 
        OPEN(66,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
                                                                        
       WRITE (66,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(phi), Im(phi), V, S)"     
                        !output phi                                     
       DO ii=1,gridpts 
         WRITE (66,'(5F11.5,F7.4)')(x(ii)),                            &
     &             (x(ii))*(mx3-1)/dt,                                  &
     &             phi(ii),V(ii),S                                      
       ENDDO 
      WRITE (filename,'(A,F0.3,A,I0,A)') "wf-",REAL(gf_Z(en,0)),        &
     &          "-",n,"-dphi.dat"                                       
        IF (n/=1) CLOSE (77) 
        OPEN(77,FILE=filename,                                          &
     &              FORM='formatted',STATUS='REPLACE')                  
       WRITE (77,'(A,I5,A)') "# phi for n=",n,                          &
     &                   "  (z, fft point, Re(dphi), Im(dphi), V, S)"   
         WRITE (77,'(5F11.5,F7.4)') (x(ii)),                            &
     &             (x(ii))*(mx3-1)/dt,                                  &
     &             dphi,V(ii),S                                         
!       WRITE (*,*) "d_aux ",d_aux                                      
#endif                                                                  
                                                                        
                                                                        
                            !saves in case d_aux==0 some cpu time       
       IF (d_aux > 0 ) THEN 
        !move trough delta region  & store EmbPoti                      
                ! phi(index_c) , dphi at left side                      
        k_r = sqrt(2.0*gf_Z(en,0)-s) 
        embpot(n,n)= I/2*k_r*                                           &
     &                (                                                 &
     &                         (dphi/(I*k_r)* COS(k_r*d_aux)            &
     &                         +I*phi(1)*SIN(k_r*d_aux))                &
     &                  /                                               &
     &                         ( phi(1)*COS(k_r*d_aux)                  &
     &                        +I*(dphi/(I*k_r))*SIN(k_r*d_aux) )        &
     &                 )                                                
      ELSE 
        embpot(n,n)= 0.5*dphi/phi(1) 
      ENDIF 
                                                                        
!      WRITE (0,*) en, REAL(gf_Z(en,0)),embpot(n,n)                     
             ! n                                                        
      ENDDO 
                                                                        
#ifdef CPP_DEBUG_IP    
!close wf. file                                  
      CLOSE (66) 
      CLOSE (77) 
#endif                                                                  
      RETURN 
                                                                        
      END SUBROUTINE 

      !>

                                                                        
      !<-- S: gf_vacuum_init(jspin,cell,kpts,kp,gfinp,vac_pot)          
      SUBROUTINE gf_vacuum_init(jspin,cell,kpts,lapw,gfinp,vac_pot_in,  &
     &           stars,l_noco_in)
!******************************************                             
!     Initializes the private variables of the module                   
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)             :: jspin 
      TYPE(t_cell),INTENT(IN)        :: cell 
      TYPE(t_kpts),INTENT(IN)        :: kpts 
      TYPE(t_lapw),INTENT(IN)        :: lapw 
      TYPE(t_gfinp),INTENT(IN)       :: gfinp 
      TYPE(t_stars),INTENT(IN)       :: stars 
      REAL   ,INTENT(IN)             :: vac_pot_in
      logical,intent(in)             :: l_noco_in
      !>                                                                
      coordconv = cell%bmat(1:2,1:2) 
      allocate(bk(size(kpts%bk,1),size(kpts%bk,2)))
      bk        = kpts%bk 
      k1p       => lapw%kp%k1p(:,jspin) 
      k2p       => lapw%kp%k2p(:,jspin) 
                                 ! d/2                                  
      d         = cell%z1 
                                     ! dt/2                             
      dt        = cell%amat(3,3)/2.0 
      d_aux     = dt-d 
                                 !fft box                               
      mx3       = stars%mx3 
      surface   = gfinp%l_surface 
      imagepoti = gfinp%l_imagepoti 
      imageplane= gfinp%imageplane 
      imageboundary= gfinp%imageboundary 
      Efield    = gfinp%Efield 
      vac_pot   = vac_pot_in 
      l_noco=l_noco_in
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
      !<-- S: gf_vacuum_Poti_init(                                      
      SUBROUTINE gf_vacuum_Poti_init()
      use m_gf_types
      IMPLICIT NONE
!******************************************                             
!     Setup of potential if EmbPot is calculated numericaly             
!         01/24/2007       A. Hanuschkin                                
!******************************************                             

      IF (imagepoti) THEN 
       CALL gf_vacuum_image_init() 
       RETURN 
      ENDIF 
      IF (EField/=0.0) THEN 
       CALL gf_vacuum_EField_init() 
       RETURN 
      ENDIF 
      END SUBROUTINE 
      !>                                                                
                                                                        

      subroutine gf_vacuum_check_vacpot(vpw,stars,c,vac_pot)
      use m_gf_types
      use m_constants,only:pimach
      use m_gf_fft_singleton
      implicit none
      type(t_stars),intent(in)::stars
      complex,intent(in)::vpw(:,:)
      real,intent(in)::c
      real,intent(inout)::vac_pot

      integer::k,n,in
      complex::v(2*stars%mx3+1,size(vpw,2))
      v=0.0
      do k=-stars%mx3,stars%mx3
      	in=stars%ig(0,0,k)
      	if (in==0) cycle
      	if (k<0) then
      	    n=2*stars%mx3+k+2
      	else
      	    n=k+1
      	endif
      	v(n,:)=vpw(in,:)
      enddo
      v(:,1)=fft(v(:,1),inv=.true.)
      if (size(v,2)>1) v(:,2)=fft(v(:,2),inv=.true.)






      write(6,*) "Vacuum potential calculation"

	  v(1,1)=maxval(real(v(:,1)))
	  if (size(v,2)>1) v(1,2)=maxval(real(v(:,2)))

      write(6,*) "Potential:",real(v(1,:))
      k=min(size(v,2),2)
      vac_pot=sum(real(v(1,:k)))/k
      write(6,*) "Vacpot:",vac_pot

      end subroutine
                                                                        
      !<-- S: gf_vacuum_EField_image_init(                              
      SUBROUTINE gf_vacuum_EField_init() 
      IMPLICIT NONE
!******************************************                             
!     Initializes the E*z+V_vac potential                               
!         01/24/2007       A. Hanuschkin                                
!******************************************                             
      INTEGER :: cutoff,ii 
                                                                        
                                                                        
      !calculate gridspacing                                            
      xdistance = REAL(100)/REAL(gridpts) 
      xdist_sq = xdistance*xdistance 
                                                                        
      !set up potential mesh                                            
!      WRITE (6,*) ''                                                   
!      WRITE (6,*) '#EFie potential meshpoint, V, z (a.u.),z (fftbox u)'
      DO ii=1,gridpts 
         x(ii)=REAL(ii)*xdistance; 
         V(ii) = vac_pot+(EField*x(ii)) 
!        WRITE (6,'(I6,3F11.5)') ii,V(ii),(x(ii))                       
!     &         ,(x(ii))*(mx3-1)/dt  !to out file                       
      END DO 
!      WRITE (6,*) '#EField potential END'                              
!      WRITE (6,*) ''                                                   
                                                                        
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
!    SET POTENTIAL CONSTANT FOR TESTING                                 
      V0=vac_pot 
      V =vac_pot 
      WRITE (6,*) "SET POTENTIAL CONSTANT=V0=",V0," FOR TESTING" 
#endif                                                                  
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum_image_init(                                     
      SUBROUTINE gf_vacuum_image_init() 
      IMPLICIT NONE
!******************************************                             
!     Initializes the 1/z potential                                     
!                          A. Hanuschkin                                
!******************************************                             
      INTEGER :: cutoff,ii 
      REAL :: embplane 
                                                                        
                                                                        
      !calculate gridspacing                                            
      xdistance = REAL(imageboundary-imageplane)/REAL(gridpts) 
      xdist_sq = xdistance*xdistance 
                                                                        
      !set up potential mesh                                            
      IF (EField==0.0) THEN 
       DO ii=1,gridpts 
          x(ii)=REAL(ii)*xdistance; 
          V(ii) = vac_pot-(1/(4.0*x(ii))); 
!         V(1)=vac_pot                                                  
!         V(2)=vac_pot                                                  
       END DO 
      ELSE 
        embplane = d-d_aux 
        DO ii=1,gridpts 
          x(ii)=REAL(ii)*xdistance; 
          V(ii) = vac_pot+EField*(x(ii)+imageplane-embplane)            &
     &            -(1/(4.0*x(ii)))                                      
       END DO 
      ENDIF 
                                                                        
                                                                        
      !set up potential mesh                                            
!      WRITE (6,*) ''                                                   
!      WRITE (6,*) '#IMAGE potential meshpoint, V, z (a.u.),z (fftbox u)
!       DO ii=1,gridpts                                                 
!         WRITE (6,'(I6,3F11.5)') ii,V(ii),(imageplane+x(ii))           
!     &          ,(imageplane+x(ii))*(mx3-1)/dt  !to out file           
!      END DO                                                           
!      WRITE (6,*) '#IMAGE potential END'                               
!      WRITE (6,*) ''                                                   
                                                                        
                                                                        
                                                                        
      !set up mixed potential mesh                                      
!     WRITE (6,*) ''                                                    
!     WRITE (6,*) '#MIXEDIM. potential meshpoint, V, z, fftbox'         
      cutoff = INT(0.16/xdistance) 
      IF (cutoff>gridpts-2) THEN 

          CALL juDFT_error("CUTOFF CHOOSEN TO LARGE !",calledby="gf_vacuum.F90")
      ENDIF 
      DO ii=1,cutoff 
!         V(ii) = (1-IBETA(z))*V_LDA(ii)+IBETA(z)*V(ii);                
          V(ii) = V(cutoff+1) 
!         WRITE (6,'(I6,3F11.5)') ii,V(ii),(imageplane+x(ii))           
!    &         ,(imageplane+x(ii))*(mx3-1)/dt  !to out file             
      END DO 
      DO ii=cutoff,gridpts 
!!        V(ii) = (1-IBETA(z))*V_LDA(ii)+IBETA(z)*V(ii);                
      END DO 
                                                                        
#ifdef CPP_DEBUG_IMAGEPOTI_V0                                           
!    SET POTENTIAL CONSTANT FOR TESTING                                 
      V0=vac_pot 
      V =vac_pot 
      WRITE (6,*) "SET IMAGEPOTENTIAL CONSTANT=V0=",V0," FOR TESTING" 
#endif                                                                  
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S: gf_vacuum(gfinp,lapw,kp(:,jsp),en,nk,jsp)                 
      SUBROUTINE gf_vacuum(gfinp,lapw,en,nk,jsp) 
!******************************************                             
!     Calculate Embedding potential with a potential in Vacuum which    
!     might be different from zero. See Eq. 5.50 from my Thesis!        
!                                                                       
!     (Very ugly hack!!)                                                
!     Two additional parameters are needed:                             
!     The thickness of the auxillary volume ->stored in gfinp%dp1       
!     The value of the vacuum potential ->stored in gfinp%dp2           
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      USE m_gf_energies 
      USE m_gf_io2dmat 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_gfinp),INTENT(IN) :: gfinp 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      INTEGER,INTENT(IN)       :: en,nk,jsp 
      !>                                                                
      !<--Locals                                                        
      COMPLEX               :: embpot(lapw%nv2d,lapw%nv2d) 
      INTEGER               :: n 
      COMPLEX               :: kl,kr 
      COMPLEX               :: I 
      REAL                  :: d 
      I=CMPLX(0.0,1.0) 
      !>                                                                
      d=gfinp%dp1 
      embpot=0.0 
      DO n=1,lapw%nv2(jsp) 
         kl=SQRT(2*gf_Z(en,0)-lapw%kp%rkp(n,jsp)**2) 
         kr=SQRT(2*(gf_Z(en,0)-gfinp%dp2)-lapw%kp%rkp(n,jsp)**2) 
         embpot(n,n)=I/2*kl*(kr*COS(kl*d)+I*kl*SIN(kl*d))/(kl*COS(kl*d) &
     &        +I*kr*SIN(kl*d))                                          
      ENDDO 
      CALL gf_write2dmat(IO2D_EMB,2,2,en,nk,jsp,lapw,embpot) 
                                                                        
      END SUBROUTINE 
      !>                                                                
#endif
      END                                           
