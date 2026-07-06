!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_makepseudocharge 
      IMPLICIT NONE
                                                                        
!      PRIVATE                                                          
      REAL,PARAMETER :: gausscale = 3.0 
!      PUBLIC gf_makepseudocharge                                       
      CONTAINS 
      !<-- S: gf_makepseudocharge(layer,jspins,atoms,stars,sphhar,cell,s
      SUBROUTINE  gf_makepseudocharge(layer,jspins,atoms,stars,sphhar   &
     &     ,cell,sym,layers,noco,mpi)
!-----------------------------------------------                        
! DESC:Construct the pseudo charge                                      
!           (last modified: 07-11-20) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_iodop 
      USE m_fleur_psqpw 
      USE m_gf_stepsanaly 
      USE m_gf_plot 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER                    :: layer,jspins 
      TYPE(t_atoms),INTENT(IN)   :: atoms 
      TYPE(t_stars),INTENT(IN)   :: stars 
      TYPE(t_sphhar),INTENT(IN)  :: sphhar 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_sym),INTENT(IN)     :: sym 
      TYPE(t_layers),INTENT(IN)  :: layers 
      type(t_noco),INTENT(IN)    :: noco
      TYPE(t_mpi),INTENT(IN)     :: mpi 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX :: vr(1:2*stars%mx3+1), vrr(1:2*stars%mx3+1) 
                                                                        
      COMPLEX,ALLOCATABLE :: qpw(:,:),psq(:),ppsq(:,:) 
      REAL   ,ALLOCATABLE :: rho(:,:,:,:) 
      REAL                :: dt,z1,sigma,psq0,z2 
      INTEGER             :: gz,index,l 
      LOGICAL             :: l_exists 
      INTEGER             :: n,n2,n3
      !>                                                                
                                                                        
      ALLOCATE(qpw(stars%nq3,jspins),psq(stars%nq3)) 
                                                                        
      ALLOCATE(rho(maxval(atoms%jri),0:maxval(sphhar%nlh),atoms%ntype   &
     &     ,jspins))                                                    
                                                                        
      !<-- load the charge                                              
      CALL gf_loddop(GF_CDNFILE,layer,jspins,atoms,stars,sphhar,rho,qpw &
     &     ,noco,.FALSE.)
                                                                        
      !>                                                                
                                                                        
      IF (jspins>1) THEN 
         qpw(:,1)=qpw(:,1)+qpw(:,2) 
         rho(:,:,:,1)=rho(:,:,:,1)+rho(:,:,:,2) 
      ENDIF 
                                                                        
      CALL fleur_psqpw(mpi%self_subcom,atoms,stars,sphhar,cell,sym      &
     &     ,qpw,rho,psq)                                                
!      CALL gf_plot(layer,stars,cell,atoms,sym,1,qpw,GF_PLOT_CHARGE)!res
!     $     ,1/)),GF_PLOT_CHARGE)                                       
      !<-- convolute with step-function                                 
      ALLOCATE(ppsq(2*stars%mx3+1,stars%nq2))
      ppsq=0.0
      DO n=1,size(psq)
         n2=stars%ig2(n)
         n3=stars%kv3(3,n)
         if (n3<0) n3=n3+size(ppsq,1)
         ppsq(n3+1,n2)=psq(n)
      ENDDO


      !CALL  gf_initstepsanaly(stars,1)
      !CALL gf_steps_pot_convolve(layer,stars,psq,ppsq)
      !>                                                                
      !<--If layerfix file exists perhaps we need to fix the charge     
      INQUIRE(FILE ="layerfix",EXIST = l_exists) 
      IF (l_exists) THEN 
         l_exists = .FALSE. 
         OPEN(99,FILE ="layerfix") 
         loop:DO 
            READ(99,*,END = 100) l 
            IF (l == layer.OR.l == 0) THEN 
               l_exists = .TRUE. 
               EXIT loop 
            ENDIF 
         ENDDO loop 
  100    IF (l_exists)  ppsq(1,1) = 0.0 
         CLOSE(99) 
      ENDIF 
      !>                                                                
      ! write pseudo-charge (and charge as dummy potential)             
      WRITE(*,*) "Pseudo:",layer 
      CALL priv_plotPlanars(stars,qpw(:,1),psq,ppsq) 
      CALL gf_iodop_writepseudo(layer,ppsq,mpi%iodop_subcom) 
      DEALLOCATE(qpw,psq,rho,ppsq) 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- F: priv_plotPlanar(stars,pw)results(planar)                  
                                                                        
      SUBROUTINE priv_plotPlanars(stars,qpw,psq,ppsq) 
!-----------------------------------------------                        
!     calculate the planar average of the charge/potential in pw        
!             (last modified: 04-08-13) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_fft_singleton 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN) :: stars 
      COMPLEX,INTENT(IN)       :: qpw(stars%nq3), psq(stars%nq3),ppsq(: &
     &     ,:)                                                          
                                                                        
      !>                                                                
      !<-- Locals                                                       
      COMPLEX :: vz(0:2*stars%mx3,3)
      INTEGER :: n 
      !>                                                                
      vz = 0.0 
      DO n =-stars%mx3,stars%mx3
         IF (stars%ig(0,0,n) == 0) CYCLE 
         IF (n<0) THEN 
            vz(n+2*stars%mx3+1,1) = qpw(stars%ig(0,0,n))
            vz(n+2*stars%mx3+1,2) = psq(stars%ig(0,0,n))
         ELSE 
            vz(n,1)             = qpw(stars%ig(0,0,n)) 
            vz(n,2)             = psq(stars%ig(0,0,n)) 
         ENDIF 
      ENDDO 
      vz(:,1) = fft(vz(:,1),inv = .TRUE.) 
      vz(:,2) = fft(vz(:,2),inv = .TRUE.) 
      vz(:,3) = fft(ppsq(:,1),inv = .TRUE.) 
      WRITE(6,*) "Constructed pseudo-charge" 
      WRITE(6,*) "z      qpw          psq       ppsq" 
      DO n = -stars%mx3,stars%mx3
         IF (n<0) THEN 
            WRITE(6,"(i5,1x,3(f0.7,1x))") n,REAL(vz(n+2*stars%mx3+1,1))   &
     &           ,REAL(vz(n+2*stars%mx3+1,2)),REAL(vz(n+2*stars%mx3+1,3))
         ELSE 
            WRITE(6,"(i5,1x,3(f0.7,1x))") n,REAL(vz(n,1)),REAL(vz(n,2)) &
     &           ,REAL(vz(n,3))                                         
         ENDIF 
      ENDDO 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
