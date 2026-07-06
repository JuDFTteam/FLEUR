!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!     This procedure is much faster than the matrix inversion.          
!                                                                       
!     Frank Freimuth, February 2007                                     
!*************************************************************          
      MODULE m_gf_surfprojo 
      use m_juDFT
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_surfprojo(layer,                                    &
     &        g,                                                        &
     &        l_nohelpregion,l_inversion,l_addemb,                      &
     &        l_realenergy,jspin,lapw,lapw_gf,cell,l_noco,                      &
     &        gij)                                                      
      USE m_gf_types 
      use m_juDFT 
      USE m_gf_curvy2dprojector 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)     :: layer 
      COMPLEX,INTENT(INOUT)::g(:,:) 
      LOGICAL, INTENT(IN):: l_nohelpregion 
      LOGICAL,INTENT(IN)::l_inversion,l_addemb,l_realenergy 
      INTEGER,INTENT(IN)::jspin 
      TYPE(t_lapw),INTENT(IN)     :: lapw 
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_cell),INTENT(IN)     :: cell 
      LOGICAL,INTENT(IN)::l_noco 
      COMPLEX,INTENT(OUT)::gij(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot) 
                                                                        
      INTEGER n1 
      INTEGER matsize 
      INTEGER info,lwork 
      COMPLEX exp1,norm 
      COMPLEX, ALLOCATABLE :: surfgreen(:,:) 
      REAL,    ALLOCATABLE :: b(:,:) 
      COMPLEX, ALLOCATABLE :: cwork(:) 
      REAL,    ALLOCATABLE :: work(:) 
      COMPLEX, ALLOCATABLE :: projtwodright(:,:) 
      COMPLEX, ALLOCATABLE :: projtwodleft(:,:) 
      INTEGER, ALLOCATABLE :: ipiv(:) 
      LOGICAL isreal,issym,blasproj 
      REAL time1,time2 
                                                                        
                                                                        
                                                                        
      !<-- tests                                                        
      isreal  = .FALSE.;issym = .FALSE. 
      IF (l_realenergy) THEN 
         !if no real energy we have a general complex matrix!           
         IF (l_addemb) THEN 
            !if an embedding potential is added, the Green function will
            !be a complex matrix in real space. If we have inversion    
            !symmetry in the system this matrix will be real, hence also
            !hermitian, otherwise the matrix is complex symmetric in    
            !real space and general complex matrix in lapw              
            IF (l_inversion) THEN 
               issym = .TRUE. 
            ENDIF 
         ELSE 
            !if no embedding potential is added and the energy is real  
            !the matrix will be hermitian in any case                   
            issym = .TRUE. 
            !if we additionaly have inversion symmetry, the matrix will 
            !be real                                                    
            IF (l_inversion) THEN 
               isreal = .TRUE. 
            ENDIF 
         ENDIF 
      ENDIF 
      !>                                                                
                                                                        
      gij(:,:)=cmplx(0.0,0.0) 
      matsize=size(g,1) 
      ALLOCATE(projtwodright(matsize,2*lapw_gf%nv2_tot)) 
      projtwodright(:,:)=cmplx(0.0,0.0) 
      exp1=exp(cmplx(0.0,cell%bmat(3,3)*cell%z1)) 
      norm=cmplx(1.0/sqrt(cell%amat(3,3)),0.0) 
!********************************************************               
!set up the projector-matrix for projecting from 3D to 2D               
!********************************************************               
      IF(.NOT.l_nohelpregion)THEN 
       DO n1=1,matsize
         projtwodright(n1,lapw%kp(n1,jspin))=                         &
     &          exp1**lapw%k3(n1,jspin)*norm                          
       ENDDO 
       DO n1=1,matsize
         projtwodright(n1,lapw_gf%nv2_tot+lapw%kp(n1,jspin))=            &
     &          exp1**(-1*lapw%k3(n1,jspin))*norm                     
       ENDDO 
      ELSE 
         CALL gf_curvy2dprojector(layer,cell,lapw,lapw_gf,.TRUE.) 
         CALL gf_get_curvy2dproject(lapw,lapw_gf,projtwodright,.TRUE.) 
!          projtwodright(1:matsize,1:lapw_gf%nv2_tot)=                     
!     =        transpose(conjg(curvyproject(1:lapw_gf%nv2_tot,1:matsize,1))
!          projtwodright(1:matsize,lapw_gf%nv2_tot+1:2*lapw_gf%nv2_tot)=      
!     =        transpose(conjg(curvyproject(1:lapw_gf%nv2_tot,1:matsize,2))
      ENDIF 
       ALLOCATE(ipiv(matsize)) 
       ALLOCATE(projtwodleft(matsize,2*lapw_gf%nv2_tot)) 
       CALL cpu_time(time1) 
       projtwodleft(:,:)=projtwodright(:,:) 
       CALL cpu_time(time2) 
       time2=time2-time1 
                                                                        
!*************************************************                      
!        real matrix                                                    
!*************************************************                      
      IF (isreal) THEN 
         !<-- we have a real-matrix!                                    

         ALLOCATE(B(matsize,matsize)) 
         B = REAL(g) 
         CALL dgetrf(matsize,matsize,b,matsize,ipiv,info) 
         IF(info/=0)  CALL juDFT_error("surfprojo: sgetrf",calledby="gf_surfprojo.F90")
         g=cmplx(b,0.0) 
         DEALLOCATE(B) 
         WRITE(*,*) "surfprojo: real matrix" 
         CALL zgetrs('N',matsize,2*lapw_gf%nv2_tot,g,matsize,   &
     &                    ipiv,projtwodright,matsize,info)              
         IF(info/=0) CALL juDFT_error("surfprojo: cgetrs-real",calledby="gf_surfprojo.F90")

         !>                                                             
         !<-- complex-matrix                                            
!************************************************************           
!        complex matrix                                                 
!************************************************************           
      ELSE 
           CALL zgetrf(matsize,matsize,g,matsize,ipiv,info) 
           IF(info/=0) CALL juDFT_error("surfprojo: cgetrf",calledby="gf_surfprojo.F90")
           WRITE(*,*) "surfprojo: complex matrix" 
           CALL zgetrs('N',matsize,2*lapw_gf%nv2_tot,g,matsize, &
     &                    ipiv,projtwodright,matsize,info)              
           IF(info/=0) CALL juDFT_error("surfprojo: cgetrs",calledby="gf_surfprojo.F90")
         !>                                                             
            !isreal                                                     
      ENDIF 
                                                                        
      DEALLOCATE(ipiv) 
                                                                        
                                                                        
!*****************************************************************      
!     project the left-hand side of matrix from 3D to 2D                
!*****************************************************************      
      CALL cpu_time(time1) 
                                                                        
      INQUIRE(FILE='blasproj',EXIST=blasproj) 
      blasproj=.TRUE. 
      IF(.NOT.blasproj)THEN 
                                                                        
      DO n1=1,matsize
         gij(lapw%kp(n1,jspin),:)=gij(lapw%kp(n1,jspin),:)+         &
     &       exp1**(-1*lapw%k3(n1,jspin))*norm*projtwodright(n1,:)    
      ENDDO 
      DO n1=1,matsize
         gij(lapw%kp(n1,jspin)+lapw_gf%nv2_tot,:)=                       &
     &       gij(lapw%kp(n1,jspin)+lapw_gf%nv2_tot,:)+                   &
     &       exp1**lapw%k3(n1,jspin)*norm*projtwodright(n1,:)         
      ENDDO 
                                                                        
      ELSE 
         CALL zgemm('C','N',2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,     &
     &                   matsize,cmplx(1.0,0.0),projtwodleft,           &
     &                   matsize,projtwodright,matsize,cmplx(0.0,0.0),  &
     &                   gij,2*lapw_gf%nv2_tot)                            
      ENDIF 
      CALL cpu_time(time2) 
      time2=time2-time1 
                                                                        
      DEALLOCATE(projtwodright) 
                                                                        
      END SUBROUTINE gf_surfprojo 
      END                                           
