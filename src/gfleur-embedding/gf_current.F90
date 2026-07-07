!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_current 
      USE m_constants, ONLY: oUnit
      use m_juDFT 
      USE m_gf_math 
      IMPLICIT NONE
      !Module contains the various subroutines for the calculation of   
      !the current. Only gf_current is public!                          
      !At present the following modes are supported:                    
      !gfinp%curr = 1:Landauer Equation with two Planes                 
      !gfinp%curr = 2:Landauer Equation on a single Plane               
      !gfinp%curr = 4:Bardeen's Equation                                
      !gfinp%curr=8:new equation on single Plane                        
      !Old modes that are no longer implemented include the T&S-Matrix a
      !Any sum of these values is also possible to evaluate more than   
      !one formula.                                                     
      PRIVATE 
      PUBLIC gf_current,gf_landauer2plane,gf_channels 
      CONTAINS 
      !<-- S:gf_current(Nv2,nk,jspin,sym,cell,gfinp,bkpts)              
                                                                        
      SUBROUTINE gf_current(Nv2,nk,jspin,sym,cell,gfinp,mpi,bkpts,lapw,lapw_gf)
!******************************************                             
!     New subroutine to calculate the current using different possible  
!     formulas                                                          
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_io2dmat 
      USE m_gf_embedding 
      USE m_gf_types 
      USE m_gf_energies,ONLY:gf_noen 
      USE m_gf_phasematrix 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(IN)        :: nv2,nk,jspin 
      TYPE(t_sym),INTENT(IN)    :: sym 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      TYPE(t_gfmpi),INTENT(IN)    :: mpi 
      TYPE(t_lapw),INTENT(IN)   :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      REAL,INTENT(IN)           :: bkpts(:,:) 
                                                                        
      !>                                                                
      !<--Locals                                                        
                                                                        
      INTEGER               :: en,region 
                                                                        
      COMPLEX,ALLOCATABLE   :: embpot(:,:,:,:),gmat(:,:,:) 
                                                                        
      REAL                  :: curr(8+nv2) 
                                                                        
      LOGICAL               :: l_shift 
                                                                        
      !>                                                                
                                                                        
      ALLOCATE(embpot(nv2,nv2,2,2),gmat(nv2,nv2,2)) 
      !<--Output information                                            
                                                                        
      IF (jspin==1.AND.nk==1.AND.mpi%pe0) THEN 
         WRITE(oUnit,*) 'Calculating the current' 
         WRITE(oUnit,*) 'Formulas implemented      used' 
         WRITE(oUnit,*) 'Landauer on two planes  ',btest(gfinp%curr,0) 
         WRITE(oUnit,*) 'Landauer on one plane   ',btest(gfinp%curr,1) 
         WRITE(oUnit,*) 'Bardeen (Transfer H.)   ',btest(gfinp%curr,2) 
         WRITE(oUnit,*) 'Ishida                  ',btest(gfinp%curr,3) 
         WRITE(oUnit,*) 'Folding ImG1ImG2        ',btest(gfinp%curr,4) 
         WRITE(oUnit,*) 'Channel decomposed      ',btest(gfinp%curr,5) 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !<--Loop over all energies                                        
                                                                        
      DO en=1,gf_noen() 
                                                                        
                                                                        
         !<--Read Embedding potentials                                  
                                                                        
         DO region=1,2 
            CALL gf_GETEMB(embpot(:,:,1,region),embpot(:,:,2,region)    &
     &           ,region,en,nk,jspin,lapw,lapw_gf)                              
!            DO side=1,2                                                
!               embpot_ok(side,region)=gf_read2dmat(IO2D_EMB,region,side
!     +              ,nk,jspin,embpot(:,:,side,region))                 
!            ENDDO                                                      
         ENDDO 
                                                                        
         !>                                                             
         curr=0.0 
         !<--Landauer on two planes                                     
                                                                        
         IF(btest(gfinp%curr,0)) THEN 
            IF(.NOT.gf_read2dmat(IO2D_G12,1,1,en,nk                     &
     &           ,jspin,lapw_gf,gmat(:,:,1))) CALL                         &
     &           juDFT_error("gf_current:Data missing for Landauer2")     
            CALL gf_landauer2Plane(gmat(:,:,1),embpot(:,:               &
     &           ,1,1),embpot(:,:,2,1),curr(1),curr(2))                 
         ENDIF 
                                                                        
         !>                                                             
         !<--Landauer on one plane                                      
                                                                        
         IF(btest(gfinp%curr,1)) THEN 
            CALL gf_Landauer1Plane(embpot(:,:,1,1),                     &
     &           embpot(:,:,2,1),curr(3))                               
         ENDIF 
                                                                        
         !>                                                             
         !<--Bardeen formula                                            
                                                                        
         IF(btest(gfinp%curr,2)) THEN 
            INQUIRE(FILE='bard_shift',EXIST=l_shift) 
            IF(.NOT.gf_read2dmat(IO2D_GMAT,1,1,en,nk                    &
     &           ,jspin,lapw_gf,gmat(:,:,1))) CALL                         &
     &           juDFT_error("gf_current:Data missing for Bardeen")       
            IF (.NOT.gf_read2dmat(IO2D_GMAT,2,1,en,nk                   &
     &           ,jspin,lapw_gf,gmat(:,:,2))) CALL                         &
     &           juDFT_error("gf_current:Data missing for Bardeen")       
            IF (l_shift) THEN 
               gmat(:,:,2) = MATMUL(MATMUL((getPhaseMatrix()),          &
     &              gmat(:,:,2)),mat_inverse(getPhaseMatrix()))         
!                embpot(:,:,1,2) = MATMUL(MATMUL((getPhaseMatrix()),    
!     +        embpot(:,:,1,2)),mat_inverse(getPhaseMatrix()))          
            ENDIF 
            CALL gf_Bardeen(embpot(:,:,1,2),                            &
     &           embpot(:,:,1,2),gmat(:,:,1),gmat(:,:,2)                &
     &           ,curr(4))                                              
         ENDIF 
                                                                        
         !>                                                             
         !<--Ishida's new formula                                       
                                                                        
         IF(btest(gfinp%curr,3)) THEN 
            IF(.NOT.gf_read2dmat(IO2D_GMAT,1,1,en,nk                    &
     &           ,jspin,lapw_gf,gmat(:,:,1))) CALL                         &
     &           juDFT_error("gf_current:Data missing for Surface")       
            CALL gf_CurrIshida(embpot(:,:,1,1),gmat(:,:,1)              &
     &           ,curr(5))                                              
         ENDIF 
                                                                        
         !>                                                             
         !<--Simple Bloch-State Matching                                
                                                                        
         IF(btest(gfinp%curr,4)) THEN 
            CALL gf_foldImemb(embpot(:,:,1,1),embpot(:,:,2,1)           &
     &           ,curr(6),curr(7))                                      
         ENDIF 
                                                                        
         !>                                                             
         !<--Conductance Channel decomposition                          
                                                                        
         IF(btest(gfinp%curr,5)) THEN 
            IF(.NOT.gf_read2dmat(IO2D_G12,1,1,en,nk                     &
     &           ,jspin,lapw_gf,gmat(:,:,1))) CALL                         &
     &           juDFT_error("gf_current:Data missing for Channel")       
            CALL gf_channels(embpot(:,:,1,1),embpot(:,:,2,1),gmat(:,:,1)&
     &           ,curr(9:),curr(8),lapw)                                
         ENDIF 
                                                                        
         !>                                                             
         IF(BTEST(gfinp%curr,5)) THEN 
            CALL writetrans(en,nk,jspin,bkpts,sym,cell,8,curr,mpi)
         ELSE 
            CALL writetrans(en,nk,jspin,bkpts,sym,cell,8,curr(1:8),mpi)
         ENDIF 
      ENDDO 
                                                                        
      !>                                                                
      DEALLOCATE(gmat,embpot) 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_landauer2plane(g12,g1,g2,j1,j2)                         
                                                                        
      SUBROUTINE gf_landauer2Plane(g12,g1,g2,j1,j2) 
!*********************************************************************  
!     subroutine to calculate the current using the landauer equation   
!     on two planes                                                     
!      Daniel Wortmann                                                  
!*********************************************************************  
      IMPLICIT NONE 
      COMPLEX, INTENT(IN) ::g12(:,:),g1(:,:),g2(:,:) 
      REAL,INTENT(OUT)    ::j1,j2 
                                                                        
      COMPLEX :: A(SIZE(g1,1),SIZE(g1,1)),B(SIZE(g1,1),SIZE(g1,1)) 
                                                                        
                                                                        
      ! short version                                                   
      ! Tr(Im(G1) G12 Im(G(2)) G12^*)                                   
      ! note:                                                           
      !  G21^*(r,r') = G12^*(r',r)                                      
      !              => G21^* is given by G12^*(g',g)                   
                                                                        
      j1 = 4.*REAL(trace(MATMUL(MATMUL(imag2d(G1),g12),MATMUL(imag2d(G2)&
     &     ,TRANSPOSE(CONJG(g12))))))                                   
                                                                        
      ! longer version                                                  
      ! Re(G1 G12 G2 G21^* - G1 G12 G2^* G21^*)                         
      A=matmul(matmul(g1,g12),g2) 
      A=matmul(A,transpose(conjg(g12))) 
      B=matmul(matmul(g1,g12),transpose(conjg(g2))) 
      B=matmul(B,transpose(conjg(g12))) 
                                                                        
      j2=-2.*REAL(trace(A-B)) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_Landauer1Plane(g1,g2,j1)                                
                                                                        
      SUBROUTINE gf_Landauer1Plane(g1,g2,j1) 
!c********************************************************************* 
!c     subroutine to calculate the current from two embedding           
!c     potentials on the same plane                                     
!c                                                                      
!c                                      Daniel Wortmann                 
!c********************************************************************* 
      IMPLICIT NONE 
      COMPLEX,INTENT(IN)  ::g1(:,:),g2(:,:) 
      REAL, INTENT(OUT)   ::j1 
                                                                        
                                                                        
      COMPLEX ::G(size(g1,1),size(g1,1)) 
                                                                        
      G=mat_inverse(G1+G2) 
      j1 = 2.0*real(trace(matmul(matmul(imag2d(G1),G),matmul(imag2d(g2) &
     &     ,TRANSPOSE(CONJG(G))))))                                     
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_Bardeen(G1,G2,GR,GL,j1)                                 
                                                                        
      SUBROUTINE gf_Bardeen(GeL,GeR,GR,GL,j1) 
!******************************************                             
!     Evaluate Landauer Formula                                         
!     GeL/R Embedding potential of left/right system                    
!     GL/R Green fct of left/rigth system                               
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(INOUT)  :: geL(:,:),geR(:,:) 
      COMPLEX,INTENT(INOUT)  :: gL(:,:),gR(:,:) 
      REAL   ,INTENT(OUT)    :: j1 
      !>                                                                
      !<--Locals                                                        
      COMPLEX :: A(size(geL,1),size(geL,1)) 
      !>                                                                
      WHERE (ABS(GL)>1E20)   GL=0.0 
      WHERE (ABS(GR)>1E20)   GR=0.0 
      WHERE (ABS(geL)>1E20)  GeL=0.0 
      WHERE (ABS(geR)>1E20)  Ger=0.0 
      !Gamma = Tr(GeL+GeR)*Im(GR)+(GeL+GeR)*Im(GL)                      
      !calculate the matrix products                                    
      A=GeL+GeR 
      j1 = 1/4.*trace(MATMUL(MATMUL(A,imag2d(GR)),MATMUL(A,imag2d(GL)))) 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_CurrIshida(G1,G,j1)                                     
                                                                        
      SUBROUTINE gf_CurrIshida(G1,G,j1) 
!******************************************                             
!     Use New Formula by Ishida                                         
!     Gamma=Tr(G*ImG1)                                                  
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)  :: g(:,:),g1(:,:) 
      REAL   ,INTENT(OUT) :: j1 
      !>                                                                
                                                                        
      j1=2.*aimag(trace(matmul(g,imag2d(G1)))) 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: gf_foldImEmb(G1,G,j1,j2)                                   
                                                                        
      SUBROUTINE gf_foldImEmb(G1,G2,j1,j2) 
!******************************************                             
!     Set G = 1 in Landauer Formula                                     
!     i.e. calculate j1 = Tr(ImG1*ImG2)                                 
!     Additionally calculate j2 = max(no_of_channels)                   
!     D. Wortmann                                                       
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)  :: g1(:,:),g2(:,:) 
      REAL   ,INTENT(OUT) :: j1,j2 
      !>                                                                
      !<--Locals                                                        
      COMPLEX :: A(SIZE(g1,1),SIZE(g1,1)) 
      INTEGER :: n1,n2 
      COMPLEX :: ev(SIZE(g1,1),SIZE(g1,1)),ew(SIZE(g1,1)) 
      !>                                                                
      j1=4.*real(trace(MATMUL(IMAG2d(G1),IMAG2d(G2)))) 
                                                                        
      !<--count no of bloch states                                      
      A = IMAG2d(g1) 
      CALL eigenvalues(A,ew,ev) 
      n1 = COUNT(ABS(ew)>1E-3) 
      A = IMAG2d(g2) 
      CALL eigenvalues(A,ew,ev) 
      n2 = COUNT(ABS(ew)>1E-3) 
      j2 = MAX(n1,n2) 
      !>                                                                
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: gf_channels(g1,g2,g12,j,jtotal)                           
      SUBROUTINE gf_channels(g1,g2,g12,j,jtotal,lapw) 
!-----------------------------------------------                        
!  calculate the conductance per channel                                
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)     :: g1(:,:),g2(:,:),g12(:,:) 
      REAL   ,INTENT(OUT)    :: j(:),jtotal 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX             :: sqrtL(SIZE(g1,2),SIZE(g1,2)) 
      COMPLEX             :: sqrtR(SIZE(g1,2),SIZE(g1,2)) 
      COMPLEX             :: T(SIZE(g1,2),SIZE(g1,2)) 
      COMPLEX       :: ew(SIZE(g1,2)),ev(SIZE(g1,1),SIZE(g1,1)) 
                                                                        
      !>                                                                
                                                                        
!      Gamma = MATMUL(MATMUL(imag2d(G1),g12),MATMUL(imag2d(G2)          
!     +     ,TRANSPOSE(CONJG(g12))))                                    
!      CALL mat_eigenvalues(Gamma,ew,ev)                                
!      J = 4*ew                                                         
!      jtotal = SUM(J)                                                  
       !make channel conductance relative                               
!      j = j/jtotal                                                     
                                                                        
      !the transmission matrix is give by                               
      ! T = sqrt(Im(Sigma_L))*G_LR*sqrt(Im(Sigma_R))                    
                                                                        
!      sqrtL = mat_sqrt(imag2d(G1))                                     
!      sqrtR = TRANSPOSE(CONJG(mat_sqrt(imag2d(G2))))                   
                                                                        
!      T = matmul(sqrtL,matmul(g12,sqrtR))                              
                                                                        
      T=mat_sqrt(MATMUL(MATMUL(imag2d(G1),g12),MATMUL(imag2d(G2)        &
     &     ,TRANSPOSE(CONJG(g12)))))                                    
                                                                        
      CALL eigenvalues(T,ew,ev) 
                                                                        
      j = 4*ew*CONJG(ew) 
      jtotal = SUM(j) 
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- S: gf_old_modes  Not in use anymore                          
                                                                        
#ifdef CPP_NEVER                                                        
      !<-- S: gf_barstocurr(nv2,en,nk,jspin,bk,sym,cell,gij,l_noco)     
      SUBROUTINE gf_barstocurr(nv2,en,nk,jspin,bk,sym,cell,gij,l_noco) 
!*********************************************************************  
!     subroutine to calculate the current by the formula of             
!     Baranger+Stones                                                   
!                                                                       
!     Daniel Wortmann                                                   
!     (Jussi's first NOCO version)                                      
!*********************************************************************  
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_types
      USE m_gf_embedding,ONLY:gf_getemb 
      IMPLICIT NONE 
      INTEGER ,INTENT(IN) :: nv2,en,nk ,jspin 
      REAL, INTENT(IN)   :: bk(:,:) 
      TYPE(t_sym),INTENT(IN) :: sym 
      TYPE(t_cell),INTENT(IN) :: cell 
      COMPLEX,DIMENSION(:,:),INTENT(IN) :: gij 
      LOGICAL, INTENT(IN) :: l_noco 
                                                                        
      INTEGER :: i 
      COMPLEX :: c1(1),c2(1) 
      COMPLEX,ALLOCATABLE ::G1(:,:),G2(:,:) 
      COMPLEX,ALLOCATABLE::A(:,:),B(:,:) 
      ALLOCATE(G1(nv2,nv2),g2(nv2,nv2)) 
                                                                        
                                                                        
      IF (l_noco.AND..TRUE.) THEN 
        ALLOCATE( A(nv2/2,nv2/2),B(nv2/2,nv2/2)) 
      ELSE 
         ALLOCATE( A(nv2,nv2),B(nv2,nv2)) 
      ENDIF 
                                                                        
!                                                                       
!   Load embedding potentials                                           
!                                                                       
                                                                        
                                                                        
      CALL gf_getemb(g1,g2,1,en,nk,jspin,l_noco) 
                                                                        
!                                                                       
!    Test different formulas                                            
!                                                                       
      ! short version                                                   
      ! Tr(Im(G1) G12 Im(G(2)) G12^*)                                   
      IF (l_noco.AND..TRUE.) THEN 
! spin up up                                                            
      A = cmplx(0.0,-0.5)*(G1(:nv2/2,:nv2/2)                            &
     &  -transpose(conjg(G1(:nv2/2,:nv2/2))))                           
                                               !Im part of G1           
      A = matmul(A,gij(1:nv2/2,nv2+1:3*nv2/2)) 
                                                                        
      B = cmplx(0.0,-0.5)*(G2(:nv2/2,:nv2/2)                            &
     &   -transpose(conjg(G2(:nv2/2,:nv2/2))))                          
                                                ! Im of G2              
      B = matmul(B,conjg(gij(nv2+1:3*nv2/2,1:nv2/2))) 
                                                                        
      A=matmul(A,B) 
      c1=0.0 
      DO i=1,nv2/2 
         c1=c1+A(i,i) 
      ENDDO 
      CALL writetrans(en,nk,1,bk,sym,cell,(/(4.)*REAL(c1),-2.*REAL(c2)/)&
     &     )                                                            
! spin dn dn                                                            
      A=cmplx(0.0,-0.5)*(G1(1+nv2/2:nv2,1+nv2/2:nv2)                    &
     &   -transpose(conjg(G1(1+nv2/2:nv2,1+nv2/2:nv2))))                
                                                          !Im part of G1
      A=matmul(A,gij(1+nv2/2:nv2,3*nv2/2+1:2*nv2)) 
                                                                        
      B=cmplx(0.0,-0.5)*(G2(1+nv2/2:nv2,1+nv2/2:nv2)                    &
     &  -transpose(conjg(G2(1+nv2/2:nv2,1+nv2/2:nv2))))                 
                                                         ! Im of G2     
      B=matmul(B,conjg(gij(3*nv2/2+1:2*nv2,nv2/2+1:nv2))) 
                                                                        
      A=matmul(A,B) 
      c1=0.0 
      DO i=1,nv2/2 
         c1=c1+A(i,i) 
      ENDDO 
      CALL writetrans(en,nk,2,bk,sym,cell,(/(4.)*REAL(c1),-2.*REAL(c2)/)&
     &     )                                                            
! spin flip (up dn)                                                     
      A=cmplx(0.0,-0.5)*(G1(:nv2/2,1+nv2/2:nv2)                         &
     &   -transpose(conjg(G1(:nv2/2,1+nv2/2:nv2))))                     
                                                     !Im part of G1     
      A=matmul(A,gij(1:nv2/2,3*nv2/2+1:2*nv2)) 
                                                                        
      B=cmplx(0.0,-0.5)*(G2(:nv2/2,1+nv2/2:nv2)                         &
     &   -transpose(conjg(G2(:nv2/2,1+nv2/2:nv2))))                     
                                                     ! Im of G2         
      B=matmul(B,conjg(gij(nv2+1:3*nv2/2,nv2/2+1:nv2))) 
                                                                        
      A=matmul(A,B) 
      c1=0.0 
      DO i=1,nv2/2 
         c1=c1+A(i,i) 
      ENDDO 
      CALL writetrans(en,nk,3,bk,sym,cell,(/(4.)*REAL(c1),-2.*REAL(c2)/)&
     &     )                                                            
      ELSE 
                                                   !Im part of G1       
      A=cmplx(0.0,-0.5)*(G1-transpose(conjg(G1))) 
      A=matmul(A,gij(1:nv2,nv2+1:2*nv2)) 
                                                                        
                                                   ! Im of G2           
      B=cmplx(0.0,-0.5)*(G2-transpose(conjg(G2))) 
      B=matmul(B,transpose(conjg(gij(nv2+1:2*nv2,1:nv2)))) 
                                                                        
      A=matmul(A,B) 
      c1=0.0 
      DO i=1,nv2 
         c1=c1+A(i,i) 
      ENDDO 
                                                                        
      ! longer version                                                  
      ! Re(G1 G12 G2 G21^* - G1 G12 G2^* G21^*)                         
                                                                        
      A=matmul(matmul(g1,gij(1:nv2,nv2+1:2*nv2)),g2) 
      A=matmul(A,transpose(conjg(gij(nv2+1:2*nv2,1:nv2)))) 
      B=matmul(matmul(transpose(g1),gij(1:nv2,nv2+1:2*nv2)),            &
     &   transpose(conjg(g2)))                                          
      B=matmul(B,transpose(conjg(gij(nv2+1:2 *nv2,1:nv2)))) 
                                                                        
      c2=0.0 
      DO i=1,nv2 
         c2=c2+A(i,i)-B(i,i) 
      ENDDO 
                                                                        
                                                                        
      WRITE(*,'(a,i3,1x,i3,1x,e16.10,1x,e16.10)') 'BarSto:',en,nk,      &
     &     (4.)*real(c1),-2*real(c2)                                    
                                                                        
!     WRITE output                                                      
      CALL writetrans(en,nk,jspin,bk,sym,cell,                          &
     &    (/(4.)*REAL(c1),-2.*REAL(c2)/))                               
      ENDIF 
      DEALLOCATE(g1,g2,a,b) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
#endif                                                                  
                                                                        
                                                                        
      !>                                                                
      END                                           
