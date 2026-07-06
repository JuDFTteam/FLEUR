!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_embedding 
      use m_juDFT 
      IMPLICIT NONE
!-----------------------------------------------                        
! This module contains many subroutines needed for                      
! calculations with an embedding potential!!                            
!                                                                       
! setemb: saves the embedding pot to file                               
! getemb: gets embedding potenial from file                             
! addemb: adds embedding potential to Hamiltonian                       
!                                                                       
! checkemb: checks if provided embedding potential fits setup           
! generateEmbPot: generate Embpot from T-Matrix&eigenvectors            
!                                                                       
! added subroutines to deal with embedding-plane descriptors            
!                       (last modified: 04-06-15)                       
!-----------------------------------------------                        
      CONTAINS 
                                                                        
      !<--S: gf_setemb(region,en,nk,jspin,lapw,gij)                     

      SUBROUTINE gf_setemb(region,en,nk,jspin,lapw,lapw_gf,gij)
!********************************************************************** 
!     * This SUBROUTINE saves the embedding potential for the           
!     * region, energy and spin                                         
!     *                                                                 
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001          
!********************************************************************** 
      USE m_gf_math,ONLY:mat_inverse 
      USE m_gf_io2Dmat 
      USE m_gf_types
      IMPLICIT NONE 
      INTEGER,INTENT(IN)       :: region,en,jspin,nk 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      COMPLEX,INTENT(IN)       :: gij(:,:) 
                                                                        
      INTEGER :: nv2,side 
      COMPLEX,ALLOCATABLE::A(:,:) 
                                                                        
      nv2 = size(gij,2)/2 
      IF (2*nv2/=SIZE(gij,1)) CALL juDFT_error                            &
     &     ('Wrong dimensions in gf_setemb')                            
      ALLOCATE(A(nv2,nv2)) 
                                                                        
                                                                        
      DO side = 1,2 
         A = mat_inverse(                                               &
     &        gij(nv2*(side-1)+1:nv2*(side-1)+nv2,                      &
     &        nv2*(side-1)+1:nv2*(side-1)+nv2))                         
         !write(*,*) region,side,a(1,1)                                 
         IF (side==1) A=-1.0*A 
         CALL gf_write2Dmat(IO2D_EMB,region,side,en,nk,jspin,lapw_gf,A) 
                                                                        
                                                                        
      ENDDO 
      DEALLOCATE(A) 
      END SUBROUTINE 

      !>                                                                
      !<--S: gf_getemb(G1,G2,region,en,nk,jspin,altpot)                 
                                                                        
      SUBROUTINE gf_GETEMB(G1,G2,region,en,nk,jspin,lapw,lapw_gf,altpot)
!********************************************************************** 
!     * This SUBROUTINE returns the correct embedding potentials for the
!     * region, energy, kpoint and spin                                 
!     *                                                                 
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001          
!********************************************************************** 
      USE m_gf_types 
      USE m_gf_io2Dmat 
      USE m_gf_vacuum 
      IMPLICIT NONE 
      TYPE(t_lapw),INTENT(IN)       :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      INTEGER,INTENT(IN)            :: en,nk,jspin,region 
      COMPLEX,INTENT(OUT)           :: G1(:,:),G2(:,:) 
      LOGICAL,OPTIONAL,INTENT(IN)   :: altpot 
                                                                        
      INTEGER :: key 
      LOGICAL :: notused 
      key = IO2D_EMB 
      IF(PRESENT(altpot)) THEN 
         IF(altpot) key = IO2D_EMBADV 
      ENDIF 
                                                                        
                                                                        
      notused=gf_read2dmat(key,region,1,en,nk,jspin,lapw_gf,G1) 
      IF(.NOT.gf_read2dmat(key,region,2,en,nk,jspin,lapw_gf,G2)) THEN 
         !on second side we might have a vacuum                         
         !IF (region == 1) CALL gf_vacuum_embpot(en,nk,g2)
         CALL juDFT_warn("UUPS, no vacuum potential in gf_embedding")
          !CALL gf_vacuum_embpot(en,nk,g2)
      ENDIF 
                                                                        
      END SUBROUTINE gf_GETEMB 
      !>                                                                
      !<--S: gf_getemb2(G1,side,region,en,nk,jspin,altpot)              

      SUBROUTINE gf_GETEMB2(G1,side,region,en,nk,jspin,lapw,lapw_gf,altpot)
!********************************************************************** 
!     * This SUBROUTINE returns the correct embedding potentials for the
!     * region, energy, kpoint and spin                                 
!     *                                                                 
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001          
!********************************************************************** 
      USE m_gf_types 
      USE m_gf_io2Dmat 
      USE m_gf_vacuum_hs
      IMPLICIT NONE 
      TYPE(t_lapw),INTENT(IN)       :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      INTEGER,INTENT(IN)            :: en,nk,jspin,region,side 
      COMPLEX,INTENT(OUT)           :: G1(:,:) 
      LOGICAL,OPTIONAL,INTENT(IN)   :: altpot 
                                                                        
      INTEGER :: key 
      LOGICAL :: notused 
      key = IO2D_EMB
      IF(PRESENT(altpot)) THEN 
         IF(altpot) key = IO2D_EMBADV 
      ENDIF 
      CALL timestart("IO: reading embpot") 
                                                                        
      notused=gf_read2dmat(key,region,side,en,nk,jspin,lapw_gf,G1)
      CALL timestop("IO: reading embpot")
      IF (side == 2.AND..NOT.notused) THEN 
         CALL gf_generate_embpot(en,jspin,g1)
      ENDIF 

!      if(.not.notused) CALL juDFT_error("not usable",calledby="gf_embedding.F90")
!      IF(.NOT.gf_read2dmat(key,region,2,en,nk,jspin,lapw_gf,G2)) THEN     
!         !on second side we might have a vacuum                        
!         IF (region == 1) CALL gf_vacuum_embpot(en,nk,g2)              
!      ENDIF                                                            

      END SUBROUTINE gf_GETEMB2 

      !>                                                                
      !<--S: gf_checkemb(region,amat)                                   

      SUBROUTINE gf_checkemb(layers,amat) 
!********************************************************************** 
!     * This subroutine checks the distances to the embedding           
!     * planes are compatible in the slab and in the embedding          
!     potentials                                                        
!     *                                                                 
!     *     Jussi Enkovaara, Juelich 2003                               
!     Added code to check for embedding plane descriptors               
!     Removed code from Jussi after change in io2dmat                   
!        (last modified: 05-11-17) D.Wortmann                           
!*********************************************************************  
      USE m_gf_io2Dmat 
      USE m_gf_types 
      USE m_gf_embdesc 
      IMPLICIT NONE 
      TYPE(t_layers),INTENT(IN) :: layers 
      REAL   ,INTENT(IN)        :: amat(:,:) 
                                                                        
      !<-- Locals                                                       
      TYPE(T_embdesc) :: ds_l,ds_r,de_l,de_r,ds_tmp 
      !>                                                                
                                                                        
                                                                        
      !<-- First read setup-descriptors                                 
                                                                        
      CALL gf_readdesc_setup(1,ds_tmp,ds_l) 
      CALL gf_deallocEmbDesc(ds_tmp) 
      CALL gf_readdesc_setup(layers%num_layers,ds_r,ds_tmp) 
      CALL gf_deallocEmbDesc(ds_tmp) 
      CALL gf_readdescriptor(1,1,amat,de_l) 
      CALL gf_readdescriptor(layers%num_layers,2,amat,de_r) 
                                                                        
      !>                                                                
      IF (de_l%valid) THEN 
         IF (.NOT.gf_checkdesc(de_l,ds_l,amat)) THEN 
            WRITE(*,*) 'Setup of left embedding plane inconsistent' 
            CALL juDFT_error('Embedding-plane descriptors differ') 
         ENDIF 
      ENDIF 
      IF (de_r%valid) THEN 
         IF (.NOT.gf_checkdesc(de_r,ds_r,amat)) THEN 
            WRITE(*,*) 'Setup of right embedding plane inconsistent' 
            CALL juDFT_error('Embedding-plane descriptors differ') 
         ENDIF 
      ENDIF 
      CALL gf_deallocEmbDesc(ds_r) 
      CALL gf_deallocEmbDesc(ds_l) 
      CALL gf_deallocEmbDesc(de_r) 
      CALL gf_deallocEmbDesc(de_l) 
      END SUBROUTINE 

      !>                                                                
                                                                        
      !<-- S: gf_addemb_g1g2(jspin,lapw,bz,z1,z2,H,G1,G2,l_noco,altpot,f
                                                                        
      SUBROUTINE gf_ADDEMB_g1g2(jspin,lapw,lapw_gf,bz,z1,z2,H,G1,G2,    &
     &     l_noco,factor)                                               
                                                                        
!********************************************************************** 
!* This subroutine adds the embedding potential to the Hamiltonian      
!*                                                                      
!*                                                                      
!*                           Daniel Wortmann, Tokyo, 2001             * 
!********************************************************************** 
!   Comments on signs:                                                  
!        -we add Sigma here to e*S-H, so effectively we subtract Sigma  
!           from H                                                      
!        -we add both left and right Sigma, so we use the normal        
!           derivative convention for Sigma                             
!        -we must therefore have the same Sigma on both sides for a     
!            system with z-reflection symmetry                          
!                                                                       
      USE m_constants, ONLY: pi_const, oUnit 
      USE m_gf_types 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)          :: jspin 
      TYPE(t_lapw),INTENT(IN)     :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      REAL,INTENT(IN)             :: bz,z1,z2 
      COMPLEX,INTENT(INOUT)       :: H(:,:) 
      COMPLEX,INTENT(IN)          :: G1(:,:),G2(:,:) 
      LOGICAL, INTENT(IN)         :: l_noco 
      REAL   ,OPTIONAL,INTENT(IN) :: factor 
                                                                        
      INTEGER :: n1,n2 
      COMPLEX ::exp1,exp2,e1,e2 
      REAL ::norm 
                                                                        
!     Now loop over Matrix elements                                     
      norm = 2.0*PIMACH()/bz 
      IF(PRESENT(factor)) norm = norm*factor 
      exp1 = exp(cmplx(0.0,bz*z1)) 
      exp2 = exp(cmplx(0.0,bz*z2)) 
      IF (l_noco) THEN 
         DO n1=1,lapw_gf%nv_tot_sphere/2
            DO n2=1,lapw_gf%nv_tot_sphere/2
               e1=exp1**(lapw%k3(n2,jspin)-lapw%k3(n1,jspin)) 
               e2=exp2**(lapw%k3(n2,jspin)-lapw%k3(n1,jspin)) 
!               up up                                                   
                           !second durface !first surface               
               H(n1,n2) =                                               &
     &              (e1*G1(lapw%kp(n1,jspin),lapw%kp(n2,jspin))     &
     &              +e2*(G2(lapw%kp(n1,jspin),lapw%kp(n2,jspin))))  &
     &              /norm+H(n1,n2)                                      
!               up dn                                                   
                                        !second durface !first surface  
               H(n1+lapw_gf%nv_tot_sphere/2,n2) =                                 &
     &              (e1*G1(lapw%kp(n1,jspin)+lapw_gf%nv2_tot/2           &
     &              ,lapw%kp(n2,jspin))+e2*G2(lapw%kp(n1,jspin)     &
     &              +lapw_gf%nv2_tot/2,lapw%kp(n2,jspin)))/norm+H(n1     &
     &              +lapw_gf%nv_tot_sphere/2,n2)
!               dn up                                                   
               H(n1,n2+lapw_gf%nv_tot_sphere/2)=(e1*                              &
     &              (G1(lapw%kp(n1,jspin),lapw%kp(n2,jspin)         &
     &              +lapw_gf%nv2_tot/2))+e2*(G2(lapw%kp(n1,jspin)        &
     &              ,lapw%kp(n2,jspin)+lapw_gf%nv2_tot/2)))/norm+H(n1,n2 &
     &              +lapw_gf%nv_tot_sphere/2)
!               dn dn                                                   
               H(n1+lapw_gf%nv_tot_sphere/2,n2+lapw%nv_tot/2) = (e1*              &
     &              G1(lapw%kp(n1,jspin)+lapw_gf%nv2_tot/2,lapw%kp(n2  &
     &              ,jspin)+lapw_gf%nv2_tot/2)+e2                          &
     &              *G2(lapw%kp(n1,jspin)+lapw_gf%nv2_tot/2,lapw%kp(n2 &
     &              ,jspin)+lapw_gf%nv2_tot/2))/norm+H(n1                  &
     &              +lapw_gf%nv_tot_sphere/2,n2+lapw_gf%nv_tot_sphere/2)
            ENDDO 
         ENDDO 
      ELSE 
         DO n1=1,lapw_gf%nv_tot_sphere
            DO n2=1,lapw_gf%nv_tot_sphere
               H(n1,n2) = (exp1**(lapw%k3(n2,jspin)-lapw%k3(n1,jspin&
     &              ))*G1(lapw%kp(n1,jspin),lapw%kp(n2,jspin))      &
     &              +exp2**(lapw%k3(n2,jspin)-lapw%k3(n1,jspin      &
     &              ))*G2(lapw%kp(n1,jspin),lapw%kp(n2,jspin)))     &
     &              /norm+H(n1,n2)                                      
                                                                   !seco
            ENDDO 
         ENDDO 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<--S: gf_addemb(lapw,en,nk,jspin,region,bz,z1,z2,H,altpot,factor)
                                                                        
      SUBROUTINE gf_ADDEMB(lapw,lapw_gf,en,nk,jspin,region,bz,z1,z2     &
     &     ,H,l_noco,altpot,factor)                                     
                                                                        
!********************************************************************** 
!* This subroutine adds the embedding potential to the Hamiltonian      
!*                                                                      
!*                                                                      
!*                           Daniel Wortmann, Tokyo, 2001             * 
!********************************************************************** 
!   Comments on signs:                                                  
!        -we add Sigma here to e*S-H, so effectively we subtract Sigma  
!           from H                                                      
!        -we add both left and right Sigma, so we use the normal        
!           derivative convention for Sigma                             
!        -we must therefore have the same Sigma on both sides for a     
!            system with z-reflection symmetry                          
!                                                                       
      USE m_constants, ONLY: pi_const, oUnit 
      USE m_gf_types 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)      :: en,nk,jspin,region 
      REAL,INTENT(IN)         :: bz,z1,z2 
      COMPLEX,INTENT(INOUT)   :: H(:,:) 
      LOGICAL, INTENT(IN)     :: l_noco 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      LOGICAL,OPTIONAL,INTENT(IN) :: altpot 
                                                !Needed in gf_getlessers
      REAL   ,OPTIONAL,INTENT(IN) :: factor 
                                                                        
      COMPLEX,ALLOCATABLE::G1(:,:),G2(:,:) 
                                                                        
                                                                        
      ALLOCATE(G1(lapw_gf%nv2_tot,lapw_gf%nv2_tot),G2(lapw_gf%nv2_tot            &
     &     ,lapw_gf%nv2_tot))                                              
                                                                        
      !<-- First get the embedding potential                            
      IF (PRESENT(altpot)) THEN 
         CALL gf_GETEMB(G1,G2,region,en,nk,jspin,lapw,lapw_gf,altpot) 
      ELSE 
         CALL gf_GETEMB(G1,G2,region,en,nk,jspin,lapw,lapw_gf)
      ENDIF 
      !>                                                                
      !<-- call gf_addemb_g1g2 to actual add the emb-potential          
      IF (PRESENT(altpot)) THEN 
         IF (PRESENT(factor)) THEN 
            CALL gf_addemb_g1g2(jspin,lapw,lapw_gf,bz,z1,z2,H,G1,G2,l_noco      &
     &           ,factor)                                               
         ELSE 
            CALL gf_addemb_g1g2(jspin,lapw,lapw_gf,bz,z1,z2,H,G1,G2,l_noco) 
         ENDIF 
      ELSE 
         CALL gf_addemb_g1g2(jspin,lapw,lapw_gf,bz,z1,z2,H,G1,G2,l_noco) 
      ENDIF 
      !>                                                                
      DEALLOCATE(g1,g2) 
      END SUBROUTINE gf_addemb 
                                                                        
      !>                                                                
      !<--S: gf_generateEmbPot(en,nk,jspin,nv2,l_pe0,ev,ew,T1,T2,curr,la
                                                                        
      SUBROUTINE gf_generateEmbPot(en,nk,jspin,l_pe0,ev,ew,T1,T2        &
     &     ,curr,lapw,lapw_gf,gfinp,layer)                                      
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
!      USE m_gf_CBS_embtest                                             
      USE m_gf_math 
      USE m_gf_phasematrix 
      USE m_gf_types 
      USE m_gf_embdesc 
      USE m_gf_io2dmat 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      COMPLEX,INTENT(INOUT)    :: ev(:,:) 
      COMPLEX,INTENT(IN)       :: ew(:) 
      COMPLEX,INTENT(IN)       :: T1(:,:),T2(:,:) 
      REAL   ,INTENT(IN)       :: curr(:) 
      INTEGER,INTENT(IN)       :: en,nk,jspin 
      LOGICAL,INTENT(IN)       :: l_pe0 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      INTEGER,INTENT(IN)        :: layer 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      LOGICAL,SAVE          :: firstcall = .TRUE. 
      TYPE(t_embdesc),SAVE  :: desc_right,desc_left,desc_tmp 
      REAL                  :: miss 
      COMPLEX               :: sigma1(lapw_gf%nv2_tot,lapw_gf%nv2_tot) 
      COMPLEX               :: sigma2(lapw_gf%nv2_tot,lapw_gf%nv2_tot) 
      !>                                                                
                                                                        
      !<-- Write the descriptors to the files                           
                                                                        
      IF (firstcall.AND.layer == 0) THEN 
         firstcall = .FALSE. 
         CALL gf_readdesc_setup(SIZE(gfinp%napw),desc_right,desc_tmp) 
         CALL gf_deallocEmbDesc(desc_tmp) 
         CALL gf_readdesc_setup(1,desc_tmp,desc_left) 
         CALL gf_deallocEmbDesc(desc_tmp) 
         desc_left%dvec(3) =-1.*desc_left%dvec(3) 
         !Left side: translation vector should be inversed              
         !           z-component stays positive because of z-reflection 
         !           of definition of direction                         
         desc_left%dvec(1:2) = -1.*MATMUL(desc_left%ucell,(/gfinp%dp1   &
     &        ,gfinp%dp2/))                                             
         desc_right%dvec(1:2) = MATMUL(desc_left%ucell,(/gfinp%dp1      &
     &        ,gfinp%dp2/))                                             
         CALL gf_writedescriptor(1,1,desc_left) 
         CALL gf_writedescriptor(1,2,desc_right) 
         CALL gf_deallocEmbDesc(desc_left) 
         CALL gf_deallocEmbDesc(desc_right) 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !<-- First for region II, side 2                                  
      IF (layer == 0) THEN 
         CALL priv_make_emb(curr,ev,ew,2,2,en,nk,jspin,lapw,lapw_gf,gfinp) 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !<-- Transfer T-mat Eigenstates to II/1 = I/1                     
                                                                        
      !Comments on signs:                                               
      !     -the T-matrix is defined in terms of z-derivatives          
      !     -if we go from right to left we must choose the states with 
      !     -the sigma should be defined in terms of normal derivatives,
      !     hence we must calculate Sigma =-0.5*(\partial_z*psi)(psi)^{-
                                                                        
      ev = MATMUL(t2,ev) 
      IF (layer == 0) THEN 
         CALL priv_make_emb(curr,ev,ew,1,1,en,nk,jspin,lapw,lapw_gf,gfinp) 
      ELSE 
         CALL priv_make_emb(curr,ev,ew,layer,1,en,nk,jspin,lapw,lapw_gf,gfinp) 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      !<-- Transfer T-mat Eigenstates to I/2                            
                                                                        
      ev = MATMUL(t1,ev) 
      ! For region I, side 2                                            
      IF (layer == 0) THEN 
         CALL priv_make_emb(curr,ev,ew,1,2,en,nk,jspin,lapw,lapw_gf,gfinp) 
      ELSE 
         CALL priv_make_emb(curr,ev,ew,1,2,en,nk,jspin,lapw,lapw_gf,gfinp) 
      ENDIF 
      !>                                                                
      IF (layer/=0) RETURN 
      !<-- Test transformation of embedding-potential                   
                                                                        
      !Since (\tilde{S}_2) and S_2 are connect by lattice vector the    
      !embedding potential should be equal except for a basis transforma
      !due to the PhaseMatrix  
      IF(gf_read2dmat(IO2D_EMB,1,1,en,nk,jspin,lapw_gf,sigma1).AND.        &
     &     gf_read2dmat(IO2D_EMB,1,1,en,nk,jspin,lapw_gf,sigma2)) THEN     
         sigma2 =MATMUL(getPhaseMatrix(),sigma2)
         sigma2 = MATMUL(sigma2,mat_inverse(getPhaseMatrix()))                    
         miss = MAXVAL(ABS(sigma1-sigma2)) 
                                                                        
         IF (l_pe0) THEN 
            WRITE(oUnit,*) 'Missmatch of EmbPot:',miss 
            IF (miss>0.1) THEN 
               WRITE(*,*)                                               &
     &              "Something wrong with generated embedding potential"
               WRITE(*,*) "Check gf_embedding for details" 
               WRITE(*,*) "Large missmatch in embedding potentials" 
            ENDIF 
         ENDIF 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: priv_make_emb(curr,ev,ew,region,side,en,nk,jspin,lapw,gfin
      SUBROUTINE priv_make_emb(curr,ev,ew,region,side,en,nk,jspin,lapw,lapw_gf  &
     &     ,gfinp)                                                      
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_io2dmat 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)       :: curr(:) 
      COMPLEX,INTENT(IN)       :: ev(:,:),ew(:) 
      INTEGER,INTENT(IN)       :: region,side,en,nk,jspin 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_embinp),INTENT(IN) :: gfinp 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: indx(SIZE(curr)/2) 
      INTEGER             :: indx_full(SIZE(curr)) 
      INTEGER             :: nv2 
      COMPLEX             :: sigma(SIZE(ev,1)/2,SIZE(ev,2)/2) 
      REAL                :: sign 
      !>                                                                
      nv2 = SIZE(curr)/2 
      IF (side==2) THEN 
         sign=1.0 
      ELSE 
         sign =-1.0 
      ENDIF 
                                                                        
      indx     = priv_selectIndex(curr,sign) 
      sigma(:,:) = priv_calculate_sigma(ev(1:nv2,indx),ev(nv2           &
     &     +1:,indx),ew(indx),gfinp,lapw,lapw_gf,sign)                          
                                                                        
      !<-- Improve embpots                                              
                                                                        
         IF (.TRUE.) THEN 
            !Test if embedding potential gives resonable currents and   
            !correct it if necessary !Untested!!!                       
!            CALL gf_CBS_embtest(gpot(:,:,1,s),ev_sorted,ew_sorted      
!     +           ,curr_sorted)                                         
            CALL priv_cutembpot(lapw,lapw_gf,en,jspin,sigma(:,:)) 
         ENDIF 
                                                                        
         !>                                                             
      !<-- Write Embpots to file                                        
                                                                        
      CALL gf_write2dmat(IO2D_EMB,region,side,en,nk,jspin,lapw_gf,sigma) 
                                                                        
                                                                        
      indx_full(:nv2) = priv_selectIndex(curr,1.0) 
      indx_full(nv2+1:) = priv_selectIndex(curr,-1.0) 
                                                                        
      CALL gf_write2dmat(IO2D_CBS,region,side,en,nk,jspin,lapw_gf,ev(:     &
     &     ,indx_full))                                                 
                                                                        
      !>                                                                
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- F: priv_selectIndex(curr,dir)                                

      FUNCTION priv_selectIndex(curr,dir) RESULT(indx) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      REAL   ,INTENT(IN)     :: curr(:) 
      REAL   ,INTENT(IN)     :: dir

      INTEGER                :: indx(SIZE(curr)/2) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER             :: n,in,nv2 
                                                                        
      !>                                                                
      nv2=SIZE(curr)/2 
      in = 0 
      DO n = 1, 2*Nv2 
         IF ( curr(n)*dir>0.0 ) THEN 
            in = in + 1 
            IF (in>nv2) THEN 
                WRITE(*,*) 'CBS: Too many states for EmbPot' 
                RETURN 
            ENDIF
            indx(in) = n 
         ENDIF 
      ENDDO 
                                                                        
      IF (in/=nv2) CALL juDFT_error('CBS: Too few states for EmbPot')
                                                                        
      END FUNCTION 

      !>                                                                
      !<-- S: priv_cutembpot(lapw,en,jspin,sigma)                       
                                                                        
      SUBROUTINE  priv_cutembpot(lapw,lapw_gf,en,jspin,sigma)
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_energies 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_lapw),INTENT(IN) :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      INTEGER,INTENT(IN)      :: en,jspin 
      COMPLEX,INTENT(INOUT)   :: sigma(:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n1,n2 
      REAL                :: s,g_max 
      REAL,PARAMETER      :: fraction = 1.3 
      !>                                                                
                                                                        
                                                                        
      RETURN 
      g_max = MAXVAL(lapw_gf%rkp)*fraction 
      DO n1 = 1,lapw_gf%nv2_tot 
         DO n2 = 1,lapw_gf%nv2_tot 
            s=sqrt(lapw_gf%rkp(n1,jspin)**2+lapw_gf%rkp(n2,jspin)**2) 
            IF (s>g_max) THEN 
               IF (n1==n2) THEN 
                  !On diagonal set analytic value                       
                  sigma(n1,n2) = cmplx(0.,.5)*sqrt(gf_z(en,0)           &
     &                 -lapw_gf%rkp(n1,jspin)**2)                       
               ELSE 
                  sigma(n1,n2)=0.0 
               ENDIF 
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- F: priv_calculate_sigma(psi,delpsi,ew)                       
                                                                        
      FUNCTION priv_calculate_sigma(psi,delpsi,ew,gfinp,lapw,lapw_gf,sign) &
     &     RESULT(sigma)                                                
!-----------------------------------------------                        
!                                                                       
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      USE m_gf_math 
      USE m_gf_types 
!      USE m_gf_cbs_embtest                                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)       :: psi(:,:),delpsi(:,:) 
      COMPLEX,INTENT(IN)       :: ew(:) 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_embinp),INTENT(IN) :: gfinp 
      REAL   ,INTENT(IN)       :: sign 
      COMPLEX                  :: sigma(SIZE(psi,1),SIZE(psi,2)) 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX             :: corr_delpsi(SIZE(psi,1),SIZE(psi,2)) 
      INTEGER             :: n,nn,i 
      !>                                                                
      !<-- if kappa_max is set use different routine                    
      IF (gfinp%kappa_max<1E19) THEN 
!         CALL gf_CBS_restricted(sigma,psi,delpsi,ew,gfinp%kappa_max)   
!         sigma=sign*0.5*sigma                                          
!         RETURN                                                        
           CALL juDFT_error("Kappa_max not implemented",calledby="gf_embedding.F90")
      ENDIF 
      i=0 
      DO n = 1,SIZE(delpsi,1) 
         IF (ABS(AIMAG(ew(n)))>gfinp%kappa_max) THEN 
!            E = -0.5*AIMAG(ew(n))**2-e_shift !Pseudo energy            
            corr_delpsi(:,n)=0.0 
!            DO nn = 1,lapw_gf%nv2_tot                                     
!               corr_delpsi(nn,n) = AIMAG(SQRT(E-lapw_gf%rkp(nn,1)**2/2.
!     $              ))*psi(nn,n)                                       
!            ENDDO                                                      
         ELSE 
            i=i+1 
            corr_delpsi(:,n)=delpsi(:,n) 
         ENDIF 
      ENDDO 
!      WRITE(*,*) "States: ",lapw_gf%nv2_tot," Replaced: ",lapw_gf%nv2_tot-i  
                                                                        
      sigma = mat_inverse(psi) 
      sigma = sign*0.5*MATMUL(delpsi,sigma) 
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      END                                           
