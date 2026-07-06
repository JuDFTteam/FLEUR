!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_projembed 
      use m_juDFT
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_projembed(                                          &
     &           lapw,lapw_gf,en,nk,jspin,cell,l_fullgreen,             &
     &           l_nogno,l_nohelpregion,layer,                          &
     &           gij,                                                   &
     &           g)                                                     
!**************************************************************         
!     Obtain the surface projection of the Green's function             
!     with embedding potential from the surface projection of           
!     the Green's function without embedding potential.                 
!     Optionally: do the same with the 3D-Green function (controlled    
!     by logical l_fullgreen).                                          
!     Frank Freimuth, February 2007                                     
!**************************************************************         
      USE m_gf_types 
      USE m_gf_embedding 
      USE m_gf_energies,ONLY:direction 
      USE m_gf_spectrum 
      USE m_gf_curvy2dprojector 
                                                                        
                                                                        
      IMPLICIT NONE 
      TYPE(t_lapw),INTENT(IN)::lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      INTEGER,INTENT(IN)::en,nk,jspin 
      TYPE(t_cell),INTENT(IN)::cell 
      LOGICAL,INTENT(IN) :: l_fullgreen,l_nogno,l_nohelpregion 
      INTEGER,INTENT(IN) :: layer 
      COMPLEX,INTENT(INOUT)::gij(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot) 
      COMPLEX,INTENT(INOUT),OPTIONAL::g(:,:) 
                                                                        
      COMPLEX,ALLOCATABLE::semb(:,:) 
      COMPLEX,ALLOCATABLE::gnoughtp(:,:) 
      COMPLEX,ALLOCATABLE::pdaggnought(:,:) 
      COMPLEX,ALLOCATABLE::projtwodright(:,:) 
      COMPLEX,ALLOCATABLE::projtwodleft(:,:) 
      COMPLEX,ALLOCATABLE::gnpg(:,:) 
      COMPLEX,ALLOCATABLE::omgemb(:,:) 
      COMPLEX,ALLOCATABLE::g1(:,:),g2(:,:) 
      COMPLEX,ALLOCATABLE::sigggno(:,:) 
      COMPLEX,ALLOCATABLE::rightsig(:,:) 
      COMPLEX,ALLOCATABLE::sigggnors(:,:) 
      COMPLEX,ALLOCATABLE::rssigggno(:,:) 
      INTEGER,ALLOCATABLE::ipiv(:) 
      COMPLEX alpha,beta 
      INTEGER n1,info,matsize,n2 
      COMPLEX norm,exp1,f1 
                                                                        
                                                                        
      ALLOCATE(omgemb(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot)) 
      ALLOCATE(g1(lapw_gf%nv2_tot,lapw_gf%nv2_tot)) 
      ALLOCATE(g2(lapw_gf%nv2_tot,lapw_gf%nv2_tot)) 
                                                                        
                                                  !left emb-pot         
      CALL gf_getemb2(g1,1,layer,en,nk,jspin,lapw,lapw_gf) 
                                                  !right emb-pot        
      CALL gf_getemb2(g2,2,layer,en,nk,jspin,lapw,lapw_gf) 
                                                                        
!      WRITE(*,*) "G1:",layer,en,nk,SUM(g1)                             
!      WRITE(*,*) "G2:",layer,en,nk,SUM(g2)                             
                                                                        
                                                                        
      alpha=cmplx(1.0,0.0) 
      beta=cmplx(0.0,0.0) 
#ifdef CPP_CUOVLP                                                       
      IF(l_nohelpregion)THEN 
         g1=matmul(basisoverlaps(:,:,1),g1) 
         g2=matmul(basisoverlaps(:,:,2),g2) 
      ENDIF 
#endif                                                                  
      CALL zgemm('N','N',lapw_gf%nv2_tot,lapw_gf%nv2_tot,            &
     &          lapw_gf%nv2_tot,alpha,gij(1,1),2*lapw_gf%nv2_tot,             &
     &          g1,lapw_gf%nv2_tot,beta,omgemb(1,1),2*lapw_gf%nv2_tot)        
      CALL zgemm('N','N',lapw_gf%nv2_tot,lapw_gf%nv2_tot,            &
     &          lapw_gf%nv2_tot,alpha,gij(1,lapw_gf%nv2_tot+1),               &
     &          2*lapw_gf%nv2_tot,g2,lapw_gf%nv2_tot,beta,                    &
     &          omgemb(1,lapw_gf%nv2_tot+1),2*lapw_gf%nv2_tot)                
      CALL zgemm('N','N',lapw_gf%nv2_tot,lapw_gf%nv2_tot,            &
     &          lapw_gf%nv2_tot,alpha,gij(1+lapw_gf%nv2_tot,1),               &
     &          2*lapw_gf%nv2_tot,g1,lapw_gf%nv2_tot,beta,                    &
     &          omgemb(1+lapw_gf%nv2_tot,1),2*lapw_gf%nv2_tot)                
      CALL zgemm('N','N',lapw_gf%nv2_tot,lapw_gf%nv2_tot,            &
     &          lapw_gf%nv2_tot,alpha,                                     &
     &          gij(1+lapw_gf%nv2_tot,1+lapw_gf%nv2_tot),2*lapw_gf%nv2_tot,      &
     &          g2,lapw_gf%nv2_tot,beta,                                   &
     &          omgemb(1+lapw_gf%nv2_tot,1+lapw_gf%nv2_tot),2*lapw_gf%nv2_tot)   
                                                                        
                                                                        
                                                                        
      DO n1=1,2*lapw_gf%nv2_tot 
         omgemb(n1,n1)=omgemb(n1,n1)+cmplx(1.0,0.0) 
      ENDDO 
                                                                        
      ALLOCATE(ipiv(2*lapw_gf%nv2_tot)) 
      CALL zgetrf(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,omgemb,      &
     &          2*lapw_gf%nv2_tot,ipiv,info)                               
      IF(info/=0) CALL juDFT_error("projembed: cgetrf",calledby="gf_projembed.F90")
      IF(.NOT.l_nogno)THEN 
      CALL zgetrs('N',2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,         &
     &        omgemb,2*lapw_gf%nv2_tot,                                    &
     &        ipiv,gij,2*lapw_gf%nv2_tot,info)                             
      IF(info/=0) CALL juDFT_error("projembed: cgetrs",calledby="gf_projembed.F90")
      ENDIF 
                                                                        
      IF(l_fullgreen .OR. l_nogno)THEN 
!*************************************************************          
!     l_fullgreen => calculate the 3D embedded Green function           
!     from the 2D embedded Green function and the 3D Green              
!     function without embedding                                        
!*************************************************************          
                                                                        
!project right hand side of 3D-Green function to 2D-basis               
         matsize=lapw_gf%nmat_sphere
         IF(size(g,1)/=matsize) CALL juDFT_error("size(g,1)",calledby="gf_projembed.F90")
         IF(size(g,2)/=matsize) CALL juDFT_error("size(g,2)",calledby="gf_projembed.F90")
         ALLOCATE(gnoughtp(matsize,2*lapw_gf%nv2_tot)) 
         IF(l_fullgreen)THEN 
          IF(.NOT.l_nohelpregion)THEN 
            gnoughtp(:,:)=cmplx(0.0,0.0) 
            norm=cmplx(1.0/sqrt(cell%amat(3,3)),0.0) 
            exp1=exp(cmplx(0.0,cell%bmat(3,3)*cell%z1)) 
            DO n1=1,lapw_gf%nv_tot_sphere
              f1=norm*exp1**lapw%k3(n1,jspin) 
              DO n2=1,matsize 
               gnoughtp(n2,lapw%kp(n1,jspin))=                        &
     &            gnoughtp(n2,lapw%kp(n1,jspin))+                     &
     &              g(n2,n1)*f1                                         
                    !n2                                                 
              ENDDO 
                  !n1                                                   
            ENDDO 
            exp1=exp(cmplx(0.0,-cell%bmat(3,3)*cell%z1)) 
            DO n1=1,lapw_gf%nv_tot_sphere
              f1=norm*exp1**lapw%k3(n1,jspin) 
              DO n2=1,matsize 
               gnoughtp(n2,lapw_gf%nv2_tot+lapw%kp(n1,jspin))=           &
     &            gnoughtp(n2,lapw_gf%nv2_tot+lapw%kp(n1,jspin))+        &
     &              g(n2,n1)*f1                                         
                    !n2                                                 
              ENDDO 
                  !n1                                                   
            ENDDO 
               !l_nohelpregion                                          
          ELSE 
            ALLOCATE( projtwodright(matsize,2*lapw_gf%nv2_tot) ) 
            CALL gf_get_curvy2dproject(lapw,lapw_gf,projtwodright,.TRUE.) 
!            projtwodright(1:matsize,1:lapw_gf%nv2_tot)=                   
!     =        transpose(conjg(curvyproject(1:lapw_gf%nv2_tot,1:matsize,1))
!            projtwodright(1:matsize,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot)=    
!     =        transpose(conjg(curvyproject(1:lapw_gf%nv2_tot,1:matsize,2))
            CALL zgemm('N','N',matsize,2*lapw_gf%nv2_tot,matsize, &
     &                           cmplx(1.0,0.0),                        &
     &                           g,matsize,projtwodright,matsize,       &
     &                           cmplx(0.0,0.0),gnoughtp,matsize)       
            DEALLOCATE(projtwodright) 
          ENDIF 
         ELSE 
            gnoughtp(:,:)=gright(:,:) 
         ENDIF 
                                                                        
!prepare the calculation of the correction [G0-G]                       
         g1=transpose(g1) 
         g2=transpose(g2) 
         ALLOCATE(semb(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot)) 
         semb(:,:)=cmplx(0.0,0.0) 
         semb(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot)=                           &
     &              g1(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot)                   
         semb(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,                            &
     &                              lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot)=     &
     &              g2(1:lapw_gf%nv2_tot,1:lapw_gf%nv2_tot)                   
         CALL zgetrs('T',2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,      &
     &        omgemb,2*lapw_gf%nv2_tot,                                    &
     &        ipiv,semb,2*lapw_gf%nv2_tot,info)                            
         IF(info/=0) CALL juDFT_error("projembed: cgetrs-full",calledby="gf_projembed.F90")
                                                                        
!multiply with the 2D embedded Green function                           
         ALLOCATE(gnpg(matsize,2*lapw_gf%nv2_tot)) 
         CALL zgemm('N','T',matsize,2*lapw_gf%nv2_tot,            &
     &            2*lapw_gf%nv2_tot,cmplx(1.0,0.0),gnoughtp,matsize,       &
     &            semb,2*lapw_gf%nv2_tot,cmplx(0.0,0.0),gnpg,matsize)      
         DEALLOCATE(gnoughtp) 
         DEALLOCATE(semb) 
                                                                        
!project left hand side of 3D-Green function to 2D-basis                
         ALLOCATE(pdaggnought(2*lapw_gf%nv2_tot,matsize)) 
         IF(l_fullgreen)THEN 
          IF(.NOT.l_nohelpregion)THEN 
            pdaggnought(:,:)=cmplx(0.0,0.0) 
            exp1=exp(cmplx(0.0,-cell%bmat(3,3)*cell%z1)) 
            DO n2=1,matsize 
             DO n1=1,lapw_gf%nv_tot_sphere
              f1=norm*exp1**lapw%k3(n1,jspin) 
               pdaggnought(lapw%kp(n1,jspin),n2)=                     &
     &            pdaggnought(lapw%kp(n1,jspin),n2)+                  &
     &              g(n1,n2)*f1                                         
                   !n1                                                  
             ENDDO 
                  !n2                                                   
            ENDDO 
            exp1=exp(cmplx(0.0,cell%bmat(3,3)*cell%z1)) 
            DO n2=1,matsize 
             DO n1=1,lapw_gf%nv_tot_sphere
              f1=norm*exp1**lapw%k3(n1,jspin) 
               pdaggnought(lapw%kp(n1,jspin)+lapw_gf%nv2_tot,n2)=        &
     &            pdaggnought(lapw%kp(n1,jspin)+lapw_gf%nv2_tot,n2)+     &
     &              g(n1,n2)*f1                                         
                   !n1                                                  
             ENDDO 
                  !n2                                                   
            ENDDO 
               !l_nohelpregion                                          
          ELSE 
            ALLOCATE( projtwodleft(2*lapw_gf%nv2_tot,matsize) ) 
            CALL gf_get_curvy2dproject(lapw,lapw_gf,projtwodleft) 
!            projtwodleft(1:lapw_gf%nv2_tot,1:matsize)=                    
!     =        curvyproject(1:lapw_gf%nv2_tot,1:matsize,1)                 
!            projtwodleft(lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot,1:matsize)=     
!     =        curvyproject(1:lapw_gf%nv2_tot,1:matsize,2)                 
            CALL zgemm('N','N',2*lapw_gf%nv2_tot,matsize,matsize, &
     &                       cmplx(1.0,0.0),projtwodleft,               &
     &                       2*lapw_gf%nv2_tot,g,matsize,                  &
     &                       cmplx(0.0,0.0),pdaggnought,2*lapw_gf%nv2_tot) 
            DEALLOCATE(projtwodleft) 
          ENDIF 
         ELSE 
            pdaggnought(:,:)=gleft(:,:) 
         ENDIF 
                                                                        
!calculate embedded 3D-Green function: G=G0-[G0-G]                      
         CALL zgemm('N','N',matsize,matsize,2*lapw_gf%nv2_tot,    &
     &        cmplx(-1.0,0.0),gnpg,matsize,pdaggnought,2*lapw_gf%nv2_tot,  &
     &        cmplx(1.0,0.0),g,matsize)                                 
                                                                        
         DEALLOCATE(pdaggnought) 
         DEALLOCATE(gnpg) 
            !l_fullgreen.or.l_nogno                                     
      ENDIF 
                                                                        
      DEALLOCATE(ipiv) 
      DEALLOCATE(omgemb) 
      DEALLOCATE(g1,g2) 
      IF(allocated(gright))DEALLOCATE(gright) 
      IF(allocated(gleft))DEALLOCATE(gleft) 
      END SUBROUTINE 
      END                                           
