!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_rt_coefs 
      IMPLICIT NONE
!                                                                       
!     This module collects all different methods to calculate the       
!     Transmission and Reflection coefficients                          
!                                                                       
                                                                        
      PRIVATE 
      PUBLIC gf_wave_match 
      CONTAINS 
                                                                        
      !<-- S: gf_wave_match                                             
                                                                        
      SUBROUTINE gf_wave_match(Nv2,en,nk,jspin,sym,cell,lapw,lapw_gf,matrix     &
     &     ,gfinp,bk,mode,fmpi)
!********************************************************************** 
!     * This SUBROUTINE calculates the current through a INTERFACE      
!     * The CBS of the two bulks has to be given                        
!     * The INTERFACE is characterized by the matrix which might be     
!     * the S-matrix,the T-matrix, or the Projected Greensfunction      
!     *                                                                 
!     *                           Daniel Wortmann, Tokyo, 2001,2002     
!********************************************************************** 
      USE m_gf_writetrans,ONLY:writetrans 
      USE m_gf_types 
      USE m_gf_embedding 
      USE m_gf_io2dmat 
      USE m_gf_energies 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,      INTENT(IN) :: Nv2, jspin 
                                       !for output                      
      INTEGER,      INTENT(IN):: en,nk 
      TYPE(t_gfmpi),intent(in)::fmpi

      TYPE(t_sym),INTENT(IN)  :: sym 
      TYPE(t_cell),INTENT(IN) :: cell 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
                                         !for output                    
      REAL,         INTENT(IN):: bk(:,:) 
      COMPLEX,      INTENT(IN):: matrix(:,:) 
      TYPE(t_embinp),INTENT(IN):: gfinp 
                                      !switches mode                    
      INTEGER,      INTENT(IN):: mode 
                                                                        
!     locals                                                            
      INTEGER  :: n, n_ins,i,n_trans,n_ref,ii 
      COMPLEX  :: cbs1(2*nv2,2*nv2), cbs2(2*nv2,2*nv2) 
      COMPLEX  :: r(nv2),t(nv2) 
      REAL     :: j1(2*nv2), j2(2*nv2),j_trans,j_refl 
                                                                        
                                                                        
!                                                                       
!     READ CBS for both Systems                                         
!                                                                       
      IF (.NOT.gf_read2dmat(IO2D_CBS,1,1,en,nk,jspin,lapw_gf,cbs1)) THEN 
         CALL gf_io2dstatus(IO2D_READ,IO2D_CBS,1,1) 
         CALL juDFT_error('gf_wave_match:Failed to read first CBS') 
      ENDIF 
      IF (.NOT.gf_read2dmat(IO2D_CBS,1,2,en,nk,jspin,lapw_gf,cbs2)) THEN 
         CALL gf_io2dstatus(IO2D_READ,IO2D_CBS,1,2) 
         CALL juDFT_error('gf_wave_match:Failed to read second CBS') 
      ENDIF 
      !<-- Calculate the currents                                       
      DO n = 1, 2*Nv2 
         j1(n) = AIMAG(DOT_PRODUCT(CBS1(1:Nv2,n),CBS1(Nv2+1:2*Nv2,n))) 
         j2(n) = AIMAG(DOT_PRODUCT(CBS2(1:Nv2,n),CBS2(Nv2+1:2*Nv2,n))) 
      ENDDO 
                                                                        
      !>                                                                
      !<-- Normalize the Bloch states to unit current                   
      DO n = 1, 2*Nv2 
         CBS1(nv2+1:,n) =-1*CBS1(nv2+1:,n) 
         IF (ABS(j1(n))>gfinp%eps_current) THEN 
            CBS1(:,n) = CBS1(:,n)/SQRT(ABS(j1(n))) 
         ENDIF 
         IF (ABS(j2(n))>gfinp%eps_current) THEN 
            CBS2(:,n) = CBS2(:,n)/SQRT(ABS(j2(n))) 
         ENDIF 
      ENDDO 
      !>                                                                
      !<-- loop over all incomming states                               
      n_ins = 0 
      j_trans = 0.0 
      j_refl = 0.0 
      DO n   = nv2+1,2*nv2 
         IF (ABS(j1(n))>gfinp%eps_current ) THEN 
            !found an incomming state                                   
            n_ins = n_ins + 1 
            !<-- determine the tramsission&reflection coefficients      
            CALL gf_tmat_rt(nv2,matrix,CBS1(:,:nv2),CBS2(:,nv2+1:)      &
     &           ,cbs1(:,n),r(:),t(:))                                  
            !>                                                          
                                                                        
            !<-- Calculate the currents                                 
            n_trans=0;n_ref=0 
            DO i = 1, nv2 
               IF (ABS(j2(nv2+i))>gfinp%eps_current) THEN 
                  j_trans = j_trans + t(i)*CONJG(t(i)) 
                  n_trans = n_trans+1 
#ifndef CPP_MPI                                                         
                  IF (ABS(t(i))>gfinp%eps_current) THEN 
                     WRITE(81,"(5i4,7(1x,f0.8))") jspin,nk,en,n_ins     &
     &                    ,n_trans,gf_z(en,0),t(i),ABS(t(i))            &
     &                    ,ATAN2(REAL(t(i)),AIMAG(t(i)))                
                  ENDIF 
#endif                                                                  
               ENDIF 
            ENDDO 
            DO i = 1, nv2 
               IF (ABS(j1(i))>gfinp%eps_current) THEN 
                  j_refl = j_refl + r(i)*CONJG(r(i)) 
                  n_ref  = n_ref+1 
#ifndef CPP_MPI                                                         
                  IF (ABS(r(i))>gfinp%eps_current) THEN 
                     WRITE(82,"(5i4,7(1x,f0.8))") jspin,nk,en,n_ins     &
     &                    ,n_ref,gf_z(en,0),r(i),ABS(r(i))              &
     &                    ,ATAN2(REAL(r(i)),AIMAG(r(i)))                
                  ENDIF 
#endif                                                                  
               ENDIF 
            ENDDO 
!            WRITE (*,'(a,i4,a,i4,a,i3,a,f10.8,a,f10.8,a)') 'K:',nk     
!     $           ,' En:',en,'State:',n_ins," J_trans:",j_trans         
!     $           ,' J_refl:',j_refl                                    
                                                                        
         ENDIF 
         !>                                                             
      ENDDO 
      !>                                                                
      j_refl = n_ins-j_refl 
      IF (n_ins>0)    CALL writetrans(en,nk,jspin,bk,sym,cell,2,(/j_trans,j_refl,1.*n_ins,1.*n_trans,1.*n_ref/),fmpi)
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:gf_tmat_rt(Nv2,T_mat,pl,pr,dpl,dpr,n_ins,pin,dpin,r,t)     
      SUBROUTINE gf_tmat_rt(Nv2,T_mat,pl,pr,pin,                        &
     &     r,t)                                                         
                                                                        
!********************************************************************** 
!     subroutine to calculate the reflection and transmission coefs     
!     from the t-matrix                                                 
!     input:                                                            
!     t_mat    : t-matrix                                               
!     pl,pr    : psi_left/right                                         
!     n_ins    : no of incomming waves                                  
!     pin      : boundary values of incomming waves                     
!     output:                                                           
!     r,t      : transmission+reflection coefs                          
!                                  D. Wortmann, Tokyo 2002              
!     ***************************************************************** 
      USE m_gf_math 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,INTENT(IN) :: nv2 
      COMPLEX,INTENT(IN),DIMENSION(:,:)::pr,pl 
      COMPLEX,INTENT(IN)::pin(:) 
      COMPLEX,INTENT(IN)::t_mat(2*nv2,2*nv2) 
      COMPLEX,INTENT(OUT)::r(nv2),t(nv2) 
                                                                        
                           !temp array for r and t                      
      COMPLEX :: rt(2*nv2) 
      COMPLEX :: mat(2*nv2,2*nv2) 
                                                                        
                                                                        
                                                                        
      mat(:,1:nv2) = pr 
                                                                        
      mat(:,nv2+1:2*nv2) = MATMUL(-t_mat,pl) 
                                                                        
      rt = MATMUL(mat_inverse(T_mat),pin) 
                                                                        
      rt = lin_equation(mat,rt) 
                                                                        
      t = rt(1:nv2) 
      r = rt(nv2+1:2*nv2) 
                                                                        
                                                                        
      END SUBROUTINE gf_tmat_rt 
      !>                                                                
                                                                        
      END                                           
