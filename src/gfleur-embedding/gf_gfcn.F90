!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_gfcn 
      use m_juDFT
      IMPLICIT NONE
      CONTAINS 
      !<-- S: gf_GFCN(Nv2,en,nk,jspin,Nv, K,gfinp,l_noco, g,gij,dgij,alt
                                                                        
      SUBROUTINE gf_GFCN(l_layerprep,layer,en,nk,jspin,cell,lapw,lapw_gf,       &
     &     gfinp,l_noco,l_invs,l_sph,                                         &
     &     g,gij,dgij,altpot,real_energy)
!********************************************************************** 
!     This SUBROUTINE takes the Hamiltonian and inverts it to obtain    
!     the Green's funktion.                                             
!     The Surface-Projections of G are calculated                       
!                                                                       
!     Daniel Wortmann, Tokyo, 2001                                      
!           cleaned up, Juelich 2005                                    
!                                                                       
!***********************************************************************
! Addition of several opportunities to circumvent the time-consuming    
! matrix-inversion in special situations:                               
! No embedding potential added => opportunity to use the spectral       
! representation, which is very much faster (Switch: l_spectral).       
! In the spectral-mode an embedding potential may be taken into         
! account later on, by calling the subroutine projembed, which converts 
! the Green's function without embedding                                
! potential into the Green's function with                              
! embedding potential.                                                  
! Only the surface projection of the Green's function is needed => do   
! not obtain the Green's function by matrix inversion but solve a system
! linear equations instead, which is much faster (Routine: surfprojo)   
!                                                                       
!     Frank Freimuth, February 2007                                     
!                                                                       
!********************************************************************** 
!     Switch l_layerprep:                                               
!     .true. => calculate the T-matrix for a subsystem ("layer")        
!     .false.=> consider a complete system                              
!                                                                       
!     Frank Freimuth, November 2007                                     
!                                                                       
!********************************************************************** 
      USE m_gf_energies,ONLY: gf_z 
      USE m_gf_types 
      USE m_gf_proj ,ONLY: gf_gproj, gf_dgproj, gf_gprojnohelpregion 
      USE m_gf_embedding,ONLY: gf_getemb,gf_addemb,gf_setemb 
      USE m_gf_invertmatrix 
      USE m_gf_io2dmat 
      USE m_gf_hsdata 
      USE m_gf_projembed 
      USE m_gf_fnaeproj 
      USE m_gf_surfprojo 
      USE m_gf_addembnhr 
      USE m_gf_curvy2dprojector 
#ifdef CPP_WANNIER                                                      
      USE m_gf_selfenergy 
      USE m_gf_projcorrelate 
#endif                                                                  
      IMPLICIT NONE 
      !<-- Arguments                                                    
                                                                        
      LOGICAL, INTENT(IN)          :: l_layerprep 
      INTEGER, INTENT(IN)          :: layer 
      INTEGER,      INTENT(IN)     :: en,nk,jspin 
      TYPE(t_lapw),    INTENT(IN)  :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_cell),INTENT(IN)      :: cell 
      TYPE(t_embinp),INTENT(IN)     :: gfinp 
      LOGICAL, INTENT(IN)          :: l_noco,l_invs,l_sph
      COMPLEX,ALLOCATABLE, INTENT(INOUT) :: g(:,:) 
      COMPLEX,OPTIONAL,INTENT(OUT) :: gij(:,:,:) 
      COMPLEX,OPTIONAL,INTENT(OUT) :: dgij(:,:,:) 
      LOGICAL,OPTIONAL,INTENT(IN)  :: altpot,real_energy
                                                                        
      !>                                                                
      !<-- locals                                                       
                                                                        
      INTEGER                      :: r1,r2 
      LOGICAL,PARAMETER            :: l_numdelG = .FALSE. 
      ! if l_numdelG = .true. calculate the derivate nummerically. If no
      COMPLEX, ALLOCATABLE         :: g1(:,:),g2(:,:) 
                                                  !this is always region
      INTEGER,PARAMETER            :: region = 1 
      LOGICAL                      :: l_fullgreen, l_addemb, l_nogno    &
     &     ,l_spectral                                                  
                                                                        
      !>                                                                
      l_addemb=gfinp%l_addemb 
      l_fullgreen = gfinp%l_fullgreen 
      l_nogno = gfinp%l_nogno 
      l_spectral = gfinp%l_spectral 
      IF(l_layerprep)THEN 
                          !calculate T-matrix in this mode              
         l_addemb=.FALSE. 
                             !not needed for T-matrix                   
         l_fullgreen=.FALSE. 
                         !not needed for T-matrix                       
         l_nogno=.FALSE. 
      ELSE 
         IF (gfinp%l_dos.AND.l_nogno) THEN 
!            l_nogno = .FALSE.
!            l_fullgreen = .TRUE.
         ENDIF 
      ENDIF 
      ! calculate Green function from spectral representation           
      CALL gf_curvy2dealloc() 
      IF(l_spectral)THEN 
         CALL timestart("gf_fnaeproj") 
         CALL gf_fnaeproj(layer,                                        &
     &        jspin,gfinp,cell,lapw,lapw_gf,en,l_sph,                                 &
     &        l_invs,l_fullgreen,l_nogno,                               &
     &        gfinp%l_nohelpregion,                                     &
     &        gij(:,:,1),                                               &
     &        g,real_energy)
         CALL timestop("gf_fnaeproj") 
                                                                        
#ifdef CPP_WANNIER                                                      
         IF(gfinp%l_addselfen) THEN 
            CALL gf_projcorrelate(                                      &
     &           lapw,en,nk,jspin,cell,l_fullgreen,                     &
     &           l_nogno,                                               &
     &           gij(1:2*lapw_gf%nv2_tot,1:2*lapw_gf%nv2_tot,1),              &
     &           g)                                                     
         ENDIF 
#endif
         IF(l_addemb) THEN
            CALL timestart("gf_projembed") 
                                                                        
            CALL gf_curvy2dprojector(layer,                             &
     &           cell,lapw,lapw_gf,.TRUE.)

            CALL gf_projembed(                                          &
     &           lapw,lapw_gf,en,nk,jspin,cell,l_fullgreen,                     &
     &           l_nogno,gfinp%l_nohelpregion,layer,                    &
     &           gij(1:2*lapw_gf%nv2_tot,1:2*lapw_gf%nv2_tot,1),              &
     &           g)

            CALL timestop("gf_projembed") 
         ENDIF 
               CALL gf_write2dmat(IO2D_GMAT,region,1,en,nk,jspin,lapw_gf            &
     &     ,gij(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region))
      CALL gf_write2dmat(IO2D_GMAT,region,2,en,nk,jspin,lapw_gf            &
     &     ,gij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot                           &
     &     ,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region))
      CALL gf_write2dmat(IO2D_G12,region,1,en,nk,jspin,lapw_gf             &
     &     ,gij(1:lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region))
      CALL gf_write2dmat(IO2D_G12,region,2,en,nk,jspin,lapw_gf             &
     &     ,gij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region))
           !straightforwardly invert (H-z*S) or solve system of lin. eqs
      ELSE 
         CALL timestart("Inversion of Hamiltonian") 
                                                                        
         !<-- read h,S, make H-z*S                                      
                                                                        
         IF (.NOT.ALLOCATED(g)) THEN 
            ALLOCATE(G(lapw_gf%nmat_sphere,lapw_gf%nmat_sphere))
         ELSEIF (SIZE(G,1) /= lapw_gf%nmat_sphere) THEN
            DEALLOCATE(g) 
            ALLOCATE(G(lapw_gf%nmat_sphere,lapw_gf%nmat_sphere))
         ENDIF 
                                                                        
         CALL gf_getLargeH_eS(layer,jspin,nk,en,G) 
                                                                        
         !>                                                             
                                                                        
         !<-- add embedding potential term!                             
                                                                        
         IF (l_addemb) THEN 
          IF(.NOT.gfinp%l_nohelpregion)THEN 
            IF (PRESENT(altpot)) THEN 
               CALL gf_addemb(lapw,lapw_gf,en,nk,jspin,region                   &
     &              ,cell%bmat(3,3),-1.0*cell%z1,cell%z1,G,l_noco       &
     &              ,altpot)                                            
            ELSE 
               CALL gf_addemb(lapw,lapw_gf,en,nk,jspin,region                   &
     &              ,cell%bmat(3,3),-1.0*cell%z1,cell%z1,G,l_noco)      
            ENDIF 
               !l_nohelpregion                                          
          ELSE 
             CALL gf_addembnhr(cell,en,nk,jspin,layer,lapw,lapw_gf,g) 
          ENDIF 
         ENDIF 
                                                                        
         !>                                                             
!calculate only the surface projection gij                              
         IF(.NOT.l_fullgreen)THEN 
            CALL gf_surfprojo(layer,g,gfinp%l_nohelpregion,l_invs,      &
     &           gfinp%l_addemb,                                        &
     &           AIMAG(gf_z(en,layer)) == 0.0,jspin,lapw,lapw_gf,cell,l_noco, &
     &           gij(1:2*lapw_gf%nv2_tot,1:2*lapw_gf%nv2_tot,1))              
!calculate the full green's function by matrix inversion                
         ELSE 
            !<-- Now invert the Matrix!                                 
                                                                        
            CALL gf_invertmatrix(G,l_invs,gfinp%l_addemb,               &
     &           AIMAG(gf_z(en,layer)) == 0.0)                          
                                                                        
            !>                                                          
               !only surface projection needed?                         
         ENDIF 
         CALL timestop("Inversion of Hamiltonian") 
                                                                        
               !how is green calculated?                                
      ENDIF 
                                                                        
                                                                        
                                                                        
      ! Return if no projections are needed                             
      IF ((gfinp%curr==0).AND.(.NOT.gfinp%l_tmat).AND.                  &
     &     (.NOT.gfinp%l_gproj)) RETURN                                 
                                                                        
      !<-- Calculate G11 to G22                                         

      CALL timestart("gf_gproj") 
      IF (.NOT.PRESENT(gij)) CALL                                       &
     &     juDFT_error("Called gf_gfcn without gij")                      
                                                                        
      IF(.NOT.l_spectral.AND.l_fullgreen)THEN 
       IF(.NOT.gfinp%l_nohelpregion)THEN 
         DO r1 = 0,1 
            DO r2 = 0,1 
               CALL gf_gproj(r1,r2,jspin,lapw,lapw_gf,                  &
     &              cell,l_noco,l_sph,G,                                      &
     &              gij(1+lapw_gf%nv2_tot*r1:lapw_gf%nv2_tot+lapw_gf%nv2_tot*r1,1&
     &              +lapw_gf%nv2_tot*r2:lapw_gf%nv2_tot+lapw_gf%nv2_tot*r2,region&
     &              ))                                                  
            ENDDO 
         ENDDO 
       ELSE 
         CALL gf_gprojnohelpregion(layer,cell,lapw,lapw_gf,l_noco,g,gij(:,:,1))
       ENDIF 
      ENDIF 
                                                                        
      CALL timestop("gf_gproj") 

      !>                                                                
             !rest is no longer of use in layer code                    
      RETURN 
      !<-- write out the embedding potentials and the projections       
      IF (gfinp%l_addemb) CALL gf_setemb(region,en,nk,jspin,lapw,       &
     &     lapw_gf,gij(:,:,region))                                             
                                                                        
      CALL gf_write2dmat(IO2D_GMAT,region,1,en,nk,jspin,lapw_gf            &
     &     ,gij(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region))                  
      CALL gf_write2dmat(IO2D_GMAT,region,2,en,nk,jspin,lapw_gf            &
     &     ,gij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot                           &
     &     ,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region))                      
      CALL gf_write2dmat(IO2D_G12,region,1,en,nk,jspin,lapw_gf             &
     &     ,gij(1:lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region))   
      CALL gf_write2dmat(IO2D_G12,region,2,en,nk,jspin,lapw_gf             &
     &     ,gij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region))   
      !>                                                                
                                                                        
      !<-- Calculate derivative of projection                           
                                                                        
      IF (.NOT.PRESENT(dgij)) RETURN 
                                                                        
      CALL timestart("gf_gproj") 
      IF (.NOT.gfinp%l_addemb) THEN 
                    !zero surface derivative!                           
         dgij = 0.0 
      ELSE 
         IF (l_numdelG) THEN 
            !<--old code                                                
            CALL gf_dgproj(g,dgij(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region),&
     &           -1.*CELL%Z1,-1.*CELL%Z1,SIZE(G,1),lapw_gf%nv2_tot         &
     &           ,Lapw%nv_Tot,Cell%bmat(3,3),cell%amat(3,3),-2.0,1.0    &
     &           ,lapw%kp(:,jspin),lapw%k3(:,jspin))                
            CALL gf_dgproj(g,dgij(1:lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2       &
     &           *lapw_gf%nv2_tot,region),-1.*CELL%Z1,CELL%Z1,SIZE(G,1)    &
     &           ,lapw_gf%nv2_tot,Lapw%nv_Tot,Cell%bmat(3,3),cell%amat(3,3)&
     &           ,0.0,1.0,lapw%kp(:,jspin),lapw%k3(:,jspin))        
            CALL gf_dgproj(g,dgij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot         &
     &           ,1:lapw_gf%nv2_tot,region),cell%z1,-1.*CELL%Z1,SIZE(G,1)  &
     &           ,lapw_gf%nv2_tot,Lapw%nv_Tot,Cell%bmat(3,3),cell%amat(3,3)&
     &           ,0.0,1.0,lapw%kp(:,jspin),lapw%k3(:,jspin))        
            CALL gf_dgproj(g,dgij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot         &
     &           ,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region),CELL%Z1,CELL%Z1 &
     &           ,SIZE(G,1),lapw_gf%nv2_tot,Lapw%nv_Tot,Cell%bmat(3,3)     &
     &           ,cell%amat(3,3),2.0,1.0,lapw%kp(:,jspin),lapw%k3(: &
     &           ,jspin))                                               
            !>                                                          
         ELSE 
            !<--get embedding potential and calcalate for all sides     
                                                                        
            ALLOCATE(g1(lapw_gf%nv2_tot,lapw_gf%nv2_tot),g2(lapw_gf%nv2_tot      &
     &           ,lapw_gf%nv2_tot))                                        
            CALL gf_getemb(g1,g2,region,en,nk,jspin,lapw,lapw_gf)
            dgij(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region) = MATMUL(g1      &
     &           ,gij(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region))            
            dgij(1:lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region) = &
     &           MATMUL(g2,gij(1:lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2          &
     &           *lapw_gf%nv2_tot,region))                                 
            dgij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region) = &
     &           MATMUL(g1,gij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot            &
     &           ,1:lapw_gf%nv2_tot,region))                               
            dgij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2         &
     &           *lapw_gf%nv2_tot,region) = MATMUL(g2,gij(lapw_gf%nv2_tot+1:2 &
     &           *lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,region))   
            DEALLOCATE(g1,g2) 
                                                                        
            !>                                                          
         ENDIF 
      ENDIF 
                                                                        
      CALL timestop("gf_gproj") 
                                                                        
      !>                                                                
      END SUBROUTINE gf_gfcn 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_gfcn_analy(lapw,en,nk,jspin,region,                     
                                                                        
      SUBROUTINE gf_gfcn_analy(lapw,lapw_gf,en,nk,jspin,region,                 &
     &    gfinp,pot_aux,gij,dgij,l_noco)                                
!**********************************************************             
!     This subroutine calculates the Green function analytically for    
!     a constant potential in region II                                 
!     Daniel Wortmann, Tokyo 2002                                       
!**********************************************************             
      USE m_gf_energies  ,ONLY:gf_z 
      USE m_gf_types
      USE m_gf_embedding ,ONLY:gf_getemb,gf_setemb 
      USE m_gf_io2dmat 
      USE m_gf_math,ONLY:mat_inverse 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,         INTENT(IN)::en,nk,jspin,region 
      TYPE(t_lapw),INTENT(IN)    :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      COMPLEX,        INTENT(IN) :: pot_aux 
      TYPE(t_embinp), INTENT(IN)  :: gfinp 
      COMPLEX,       INTENT(INOUT) :: gij(:,:,:) 
      COMPLEX,       INTENT(INOUT) :: dgij(:,:,:) 
      LOGICAL,       INTENT(IN)  :: l_noco 
!     locals                                                            
      INTEGER :: n 
      COMPLEX :: k 
      REAL    :: l 
      COMPLEX,ALLOCATABLE::A(:,:),B(:,:) 
                                                                        
!      l=gfinp%d                                                        
      gij(:,:,region)=CMPLX(0.0,0.0) 
                                                                        
      IF (l_noco) THEN 
         DO n = 1,lapw_gf%nv2_tot/2 
            k    = SQRT(2*(gf_z(en,0)-pot_aux)-lapw_gf%rkp(n,jspin)**2) 
            IF (ABS(k)<3*EPSILON(0.0)) THEN 
               k = CMPLX(0.0,3*EPSILON(0.0)) 
            ENDIF 
            gij(n,n,region) = 2.0/k*COS(k*l)/SIN(k*l) 
            gij(n+lapw_gf%nv2_tot,n+lapw_gf%nv2_tot,region) = gij(n,n,region) 
            gij(n+lapw_gf%nv2_tot,n,region) = 2.0/k/SIN(k*l) 
            gij(n,lapw_gf%nv2_tot+n,region) = gij(n+lapw_gf%nv2_tot,n,region) 
                                                                        
            gij(n+lapw_gf%nv2_tot/2,n+lapw_gf%nv2_tot/2,region) = gij(n,n     &
     &           ,region)                                               
            gij(n+3*lapw_gf%nv2_tot/2,n+3*lapw_gf%nv2_tot/2,region) = gij(n,n &
     &           ,region)                                               
            gij(n+3*lapw_gf%nv2_tot/2,n+lapw_gf%nv2_tot/2,region) = gij(n     &
     &           +lapw_gf%nv2_tot,n,region)                                
            gij(n+lapw_gf%nv2_tot/2,3*lapw_gf%nv2_tot/2+n,region) = gij(n     &
     &           +lapw_gf%nv2_tot,n,region)                                
         ENDDO 
      ELSE 
         DO n = 1,lapw_gf%nv2_tot 
            k = SQRT(2*(gf_z(en,0)-pot_aux)-lapw_gf%rkp(n,jspin)**2) 
            IF (ABS(k)<3*EPSILON(0.0)) THEN 
               k = CMPLX(0.0,3*EPSILON(0.0)) 
            ENDIF 
            gij(n,n,region) = 2.0/k*COS(k*l)/SIN(k*l) 
            gij(n+lapw_gf%nv2_tot,n+lapw_gf%nv2_tot,region) = gij(n,n,region) 
            gij(n+lapw_gf%nv2_tot,n,region) = 2.0/k/SIN(k*l) 
            gij(n,lapw_gf%nv2_tot+n,region) = gij(n+lapw_gf%nv2_tot,n,region) 
         ENDDO 
      ENDIF 
      ! The derivative is always assumed to be zero (FIX this later!)   
      dgij(:,:,region) = CMPLX(0.0,0.0) 
                                                                        
                               ! no embedding potential                 
      IF (gfinp%l_addemb) THEN 
                                                                        
!                                                                       
!     Use the Dyson-equation to calculate the Green-function with the   
!     Embedding potential                                               
!                                                                       
         ALLOCATE(A(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot),B(2*lapw_gf%nv2_tot,2   &
     &        *lapw_gf%nv2_tot))                                           
         B = CMPLX(0.0) 
         CALL gf_getemb(B(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot),B(lapw_gf%nv2_tot &
     &        +1:2*lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot),region,en&
     &        ,nk,jspin,lapw,lapw_gf)                                           
         A = MATMUL(gij(:,:,region),B) 
         DO n = 1,2*lapw_gf%nv2_tot 
            A(n,n) = A(n,n)+CMPLX(1.0,0.0) 
         ENDDO 
         A = mat_inverse(A) 
         gij(:,:,region) = MATMUL(A,gij(:,:,region)) 
         DEALLOCATE(A,B) 
                                                                        
         CALL gf_setemb(region,en,nk,jspin,lapw,lapw_gf,gij(:,:,region))
      ENDIF 
      CALL gf_write2dmat(IO2D_GMAT,region,1,en,nk,jspin,lapw_gf            &
     &     ,gij(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot,region))                  
      CALL gf_write2dmat(IO2D_GMAT,region,2,en,nk,jspin,lapw_gf            &
     &     ,gij(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,lapw_gf%nv2_tot+1:2          &
     &     *lapw_gf%nv2_tot,region))                                       
                                                                        
      END SUBROUTINE gf_gfcn_analy 
                                                                        
      !>                                                                
      END                                           
