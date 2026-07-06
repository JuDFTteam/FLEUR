!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_embdos 
          IMPLICIT NONE
                                                                        
      CONTAINS 
      !<-- S:gf_embdos(jspins,lapw,kpts,cell,noco,gfinp)                
      SUBROUTINE gf_embdos(jspins,lapw,kpts,cell,noco,gfinp) 
!-----------------------------------------------                        
! DESC:Calculate eigenvalues of imag(Sigma) and interpret this as a kind
!                 Daniel Wortmann, (08-04-17)                           
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_energies 
      USE m_gf_math 
      USE m_gf_embedding 
      USE m_gf_apws 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)         :: jspins 
      TYPE(t_lapw),INTENT(INOUT) :: lapw 
      TYPE(t_kpts),INTENT(IN)    :: kpts 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_noco),INTENT(IN)    :: noco 
      TYPE(t_gfinp),INTENT(IN)   :: gfinp 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: jspin,nk,en,n,in1,in2 
      COMPLEX,ALLOCATABLE :: sigma(:,:),ew(:),ev(:,:) 
      REAL,PARAMETER      :: smallest = 1E-9 
      !>                                                                
      OPEN(99,FILE ="gf_embdos") 
      OPEN(98,FILE ="gf_embdos_S") 
      OPEN(97,FILE ="gf_embdos_AS") 
      DO jspin = 1,jspins 
         DO nk = 1,kpts%nkpts 
            CALL gf_apws(                                               &
     &           jspins,jspin,kpts%bk(:,nk),.FALSE.,noco,gfinp,cell,lapw&
     &           ,1)                                                    
            DO n = 1,lapw_gf%nv2(1) 
               IF (lapw_gf%k1p(n,1) == 1.AND.lapw_gf%k2p(n,1) == 2) in1 &
     &              = n                                                 
               IF (lapw_gf%k1p(n,1) == 2.AND.lapw_gf%k2p(n,1) == 1) in2 &
     &              = n                                                 
            ENDDO 
            ALLOCATE(sigma(lapw_gf%nv2_tot,lapw_gf%nv2_tot)) 
            ALLOCATE(ev(lapw_gf%nv2_tot,lapw_gf%nv2_tot)) 
            ALLOCATE(ew(lapw_gf%nv2_tot)) 
            DO en = 1,gf_Noen() 
               sigma = 0.0 
               CALL gf_getemb2(sigma,1,1,en,nk,jspin,lapw) 
               sigma = imag2d(sigma) 
               CALL eigenvalues(sigma,ew,ev) 
               WRITE(98,"(4i5,999(1x,f0.6))") jspin,nk,en,COUNT((ABS(ew &
     &              )>smallest).AND.ABS(ev(in1,:) - ev(in2,:))<smallest)
               WRITE(97,"(4i5,999(1x,f0.6))") jspin,nk,en,COUNT((ABS(ew &
     &              )>smallest).AND.ABS(ev(in1,:) + ev(in2,:))<smallest)
               WRITE(99,"(4i5,999(1x,f0.6))") jspin,nk,en,COUNT(ABS(ew  &
     &              )>smallest),SUM(PACK(ABS(ew),ABS(ew)>smallest))     &
     &              ,PACK(ABS(ew),ABS(ew)>smallest)                     
            ENDDO 
            DEALLOCATE(sigma,ew,ev) 
         ENDDO 
      ENDDO 
      CLOSE(99) 
      CLOSE(98) 
      CLOSE(97) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
