!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_potproject 
          IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Projectors used to calculate Potential on curvy planes           
!                 Daniel Wortmann, (08-01-28)                           
!-----------------------------------------------                        
      CONTAINS 
      !<-- S: gf_makePotProjector(layer,cell,stars,derivative,projector)
                                                                        
      SUBROUTINE gf_makePotProjector(layer,layers,cell,stars,derivative &
     &     ,curvy2dproj)                                                
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_circlestep 
      USE m_gf_embdesc 
      USE m_gf_types 
      USE m_gf_curvystep 
      use m_juDFT 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_stars),INTENT(IN)  :: stars 
      INTEGER,INTENT(IN)        :: layer 
      LOGICAL,INTENT(IN)        :: derivative 
                                                                        
      COMPLEX,     INTENT(OUT)  :: curvy2dproj(:,:,:) 
                                                                        
      !>                                                                
                                                                        
      !<-- Locals                                                       
                                                                        
      COMPLEX             :: exp1 
      COMPLEX             :: norm(-stars%mx3:stars%mx3) 
      COMPLEX             :: normcircle(-stars%mx3:stars%mx3) 
      INTEGER             :: side,n1,n2 
      COMPLEX,ALLOCATABLE :: circlestep(:,:) 
      COMPLEX,ALLOCATABLE :: curvystep (:,:,:) 
      INTEGER             :: k1max,k2max,k3max,k1,k2,n,index 
                                                                        
      TYPE(t_embdesc),TARGET  :: embdesc_left,embdesc_right 
      TYPE(t_embdesc),POINTER :: embdesc 
                                                                        
      !>                                                                
                                                                        
      CALL gf_readdesc_setup(layer,embdesc_right,embdesc_left) 
      IF (.NOT.(embdesc_left%valid.AND.embdesc_right%valid)) CALL       &
     &     juDFT_error("No valid embedding plane descriptors found")      
                                                                        
                                                                        
      norm = cmplx(1.0,0.0) 
                                                                        
      IF(derivative) THEN 
         !if the dreicative is needed, put additional factor into norm  
         !see code for normcicle below                                  
         DO n =-stars%mx3,stars%mx3 
            norm(n) =norm(n)*cmplx(0.0,n*cell%bmat(3,3)) 
         ENDDO 
      ENDIF 
                                                                        
      k1max = MAXVAL(stars%kv2(1,:size(curvy2dproj,1))) 
      k2max = MAXVAL(stars%kv2(2,:size(curvy2dproj,1))) 
      k3max = stars%mx3 
      ALLOCATE( circlestep(-2*k1max:2*k1max,-2*k2max:2*k2max) ) 
      ALLOCATE(curvystep(-2*k1max:2*k1max,-2*k2max:2*k2max,-k3max:k3max)&
     &     )                                                            
      curvy2dproj(:,:,:)=cmplx(0.0,0.0) 
      DO side = 1,2 
         IF (side == 1) THEN 
            exp1  = EXP(CMPLX(0.0,cell%bmat(3,3)*layers%c(layer)/(2.))) 
            embdesc => embdesc_left 
         ELSE 
            exp1  = EXP(CMPLX(0.0,cell%bmat(3,3)*layers%c(layer)/(-2.))) 
            embdesc => embdesc_right 
         ENDIF 
         !generate the 2d-kronecker delta                               
         DO n2 = 1,size(curvy2dproj,1) 
            DO n =-stars%mx3,stars%mx3 
               index = stars%ig(stars%kv2(1,n2),stars%kv2(2,n2),n) 
               IF (index == 0) CYCLE 
               curvy2dproj(n2,index,side) =curvy2dproj(n2,index,side)+  &
     &              norm(n)/exp1**n                                     
            ENDDO 
         ENDDO 
         IF (embdesc%cut_atoms_in == 0.AND.embdesc%cut_atoms_out == 0)  &
     &        CYCLE                                                     
         !subtract the circle&add the integrals over the sphere surfaces
         normcircle = CMPLX(SQRT(cell%amat(3,3))/cell%omtil,0.0)*norm   &
     &        *SQRT(cell%amat(3,3))                                     
                                                                        
         circlestep = CMPLX(0.0,0.0) 
         CALL gf_circlestep_embdesc(                                    &
     &        cell%bmat,                                                &
     &        embdesc,                                                  &
     &        2*k1max,2*k2max,                                          &
     &        circlestep)                                               
                                                                        
         curvystep = CMPLX(0.0,0.0) 
         CALL gf_curvystep_embdesc(                                     &
     &        cell%bmat,                                                &
     &        embdesc,side==2,                                          &
     &        2*k1max,2*k2max,k3max,                                    &
     &        curvystep)                                                
         DO n1 = 1,size(curvy2dproj,1) 
            DO n2 = 1,size(curvy2dproj,1) 
               k1 = stars%kv2(1,n2)-stars%kv2(1,n1) 
               k2 = stars%kv2(2,n2)-stars%kv2(2,n1) 
               DO n =-stars%mx3,stars%mx3 
                  index = stars%ig(stars%kv2(1,n2),stars%kv2(2,n2),n) 
                  IF (index == 0) CYCLE 
                  curvy2dproj(n1,index,side) =                          &
     &                 curvy2dproj(n1,index,side)-                      &
     &                 normcircle(n)/exp1**n*                           &
     &                 (circlestep(k1,k2)-curvystep(k1,k2,-n))          
                     !n                                                 
               ENDDO 
                     !n1                                                
            ENDDO 
                     !n2                                                
         ENDDO 
            !side                                                       
      ENDDO 
      DEALLOCATE( circlestep,curvystep) 
      CALL gf_deallocEmbDesc(embdesc_left) 
      CALL gf_deallocEmbDesc(embdesc_right) 
                                                                        
                                                                        
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: gf_make_totalPotProj(cell,layers,stars,d_total,d_aux,curvy
      SUBROUTINE  gf_make_totalPotProj(surface,cell,layers,stars,d_total&
     &     ,d_aux,delta_point,curvy2dproj)                              
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_circlestep 
      USE m_gf_embdesc 
      USE m_gf_types 
      USE m_gf_curvystep 
      use m_juDFT 
      USE m_constants 
      IMPLICIT NONE 
      !<--Arguments                                                     
      LOGICAL,INTENT(IN)        :: surface 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_stars),INTENT(IN)  :: stars 
      REAL   ,INTENT(IN)        :: d_total,d_aux,delta_point 
                                                                        
      COMPLEX,     INTENT(OUT)  :: curvy2dproj(:,:,:,:) 
      !>                                                                
                                                                        
      !<-- Locals                                                       
                                                                        
      COMPLEX             :: exp1,exp2 
      COMPLEX             :: normcircle 
      INTEGER             :: side,n1,n2 
      COMPLEX,ALLOCATABLE :: circlestep(:,:) 
      COMPLEX,ALLOCATABLE :: curvystep (:,:,:) 
      INTEGER             :: k1max,k2max,k3max,k1,k2,n,index,mx3 
      INTEGER,ALLOCATABLE :: gz_map(:) 
      REAL                :: bmat(3,3) 
                                                                        
      TYPE(t_embdesc),TARGET  :: embdesc_left,embdesc_right,embdesc_tmp 
      TYPE(t_embdesc),POINTER :: embdesc 
      LOGICAL                 :: ok 
                                                                        
      !>                                                                
                                                                        
      IF(priv_readtotalprojector(curvy2dproj)) RETURN 
                                                                        
      CALL gf_readdesc_setup(1,embdesc_tmp,embdesc_left) 
      CALL gf_deallocEmbDesc(embdesc_tmp) 
      CALL gf_readdesc_setup(layers%num_layers,embdesc_right,embdesc_tmp&
     &     )                                                            
      CALL gf_deallocEmbDesc(embdesc_tmp) 
                                                                        
      IF (.NOT.(embdesc_left%valid.AND.embdesc_right%valid)) CALL       &
     &     juDFT_error("No valid embedding plane descriptors found")      
                                                                        
                                                                        
                                                                        
      bmat = 0.0 
      bmat(1:2,1:2) = cell%bmat(1:2,1:2) 
      bmat(3,3) = 2*pi_const/d_total 
                                                                        
                                                                        
      k1max = MAXVAL(stars%kv2(1,:size(curvy2dproj,1))) 
      k2max = MAXVAL(stars%kv2(2,:size(curvy2dproj,1))) 
      mx3   = SIZE(curvy2dproj,3) 
      k3max = mx3/2+1 
      ALLOCATE(gz_map(mx3)) 
      DO n = 1,mx3 
         IF(n-1<=mx3/2) THEN 
            gz_map(n) = n-1 
         ELSE 
            gz_map(n) = n-mx3-1 
         ENDIF 
      ENDDO 
      ALLOCATE( circlestep(-2*k1max:2*k1max,-2*k2max:2*k2max) ) 
      ALLOCATE(curvystep(-2*k1max:2*k1max,-2*k2max:2*k2max,-k3max:k3max)&
     &     )                                                            
      curvy2dproj(:,:,:,:)=cmplx(0.0,0.0) 
      DO side = 1,2 
         IF (side == 2) THEN 
            exp1  = EXP(CMPLX(0.0,2*pi_const/d_total*(d_total-d_aux) )) 
            embdesc => embdesc_right 
         ELSE 
            exp1  = EXP(CMPLX(0.0,2*pi_const/d_total*d_aux)) 
            embdesc => embdesc_left 
         ENDIF 
         !generate the 2d-kronecker delta                               
         DO n2 = 1,size(curvy2dproj,1) 
            DO n = 1,mx3 
               curvy2dproj(n2,n2,n,side) = 1.0*exp1**gz_map(n) 
            ENDDO 
         ENDDO 
         !subtract the circle&add the integrals over the sphere surfaces
         normcircle = CMPLX(SQRT(cell%amat(3,3))/cell%omtil,0.0) 
                                                                        
                                                                        
         circlestep = CMPLX(0.0,0.0) 
         CALL gf_circlestep_embdesc(                                    &
     &        bmat,                                                     &
     &        embdesc,                                                  &
     &        2*k1max,2*k2max,                                          &
     &        circlestep)                                               
                                                                        
         curvystep = CMPLX(0.0,0.0) 
         CALL gf_curvystep_embdesc(                                     &
     &        bmat,                                                     &
     &        embdesc,side==2,                                          &
     &        2*k1max,2*k2max,k3max,                                    &
     &        curvystep)                                                
         DO n1 = 1,size(curvy2dproj,1) 
            DO n2 = 1,size(curvy2dproj,1) 
               k1 = stars%kv2(1,n2)-stars%kv2(1,n1) 
               k2 = stars%kv2(2,n2)-stars%kv2(2,n1) 
               DO n = 1,mx3 
                  curvy2dproj(n1,n2,n,side) =                           &
     &                 curvy2dproj(n1,n2,n,side)-                       &
     &                 normcircle*exp1**(gz_map(n))*                    &
     &                 (circlestep(k1,k2)-curvystep(k1,k2,-gz_map(n)))  
                     !n                                                 
               ENDDO 
                     !n1                                                
            ENDDO 
                     !n2                                                
         ENDDO 
            !side                                                       
      ENDDO 
                                                                        
      IF (surface) THEN 
         exp1 = EXP(CMPLX(0.0,2.*pi_const/d_total*(d_total-d_aux-7      &
     &        *delta_point)))                                           
         exp2 = EXP(CMPLX(0.0,2.*pi_const/d_total*(d_total-d_aux-8      &
     &        *delta_point)))                                           
         DO n  = 1,mx3 
            curvy2dproj(1,1,n,2) = (exp1**gz_map(n)-exp2**gz_map(n))    &
     &           /delta_point                                           
         ENDDO 
       ENDIF 
                                                                        
      DEALLOCATE( circlestep,curvystep,gz_map) 
      CALL gf_deallocEmbDesc(embdesc_left) 
      CALL gf_deallocEmbDesc(embdesc_right) 
      CALL priv_writetotalprojector(curvy2dproj) 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- F:priv_readtotalprojector(proj)                              
      FUNCTION priv_readtotalprojector(proj) RESULT(ok) 
!-----------------------------------------------                        
!  try to read the totalprojector from gf_setup/layers                  
!           (last modified:08-06-10) D. Wortmann                        
!-----------------------------------------------                        
      USE m_hdf_tools 
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(OUT)    :: proj(:,:,:,:) 
      LOGICAL                :: ok 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: dim(5),hdferr 
      INTEGER(HID_T)     :: ffid,fid,did 
                                                                        
      !>                                                                
      ok = .FALSE. 
      CALL io_hdfopen('gf_setup.hdf',H5F_ACC_RDONLY_F,ffid,hdferr) 
      CALL io_gopen(ffid,"layers",fid,hdferr) 
      IF (.NOT.io_dataexists(fid,"totalprojector")) THEN 
         !Step function already calculated                              
         CALL io_gclose (fid, hdferr) 
         CALL io_hdfclose(ffid,hdferr) 
         RETURN 
      ENDIF 
                                                                        
      CALL io_dopen(fid,"totalprojector",did,hdferr) 
      CALL io_datadim(did,dim) 
                                                                        
      IF (PRODUCT(dim) == 2*SIZE(proj)) THEN 
                                                                        
         CALL io_read(did,(/-1,1,1,1,1/),(/1,size(proj,1),size(proj,2)  &
     &        ,SIZE(proj,3),SIZE(proj,4)/),proj)                        
         ok = .TRUE. 
      ELSE 
                                                                        
         WRITE(*,*) PRODUCT(dim),SIZE(proj) 
      ENDIF 
      CALL io_dclose (did, hdferr) 
      CALL io_gclose (fid, hdferr) 
      CALL io_hdfclose(ffid,hdferr) 
                                                                        
      END FUNCTION 
      !>                                                                
      !<-- S:priv_writetotalprojector(proj,ok)                          
                                                                        
      SUBROUTINE priv_writetotalprojector(proj) 
!-----------------------------------------------                        
!  try to write the totalprojector from gf_setup/layers                 
!           (last modified:08-06-10) D. Wortmann                        
!-----------------------------------------------                        
      USE m_hdf_tools 
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)    :: proj(:,:,:,:) 
                                                                        
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: hdferr 
      INTEGER(HID_T)      :: ffid,fid,did 
                                                                        
      !>                                                                
                                                                        
      CALL io_hdfopen('gf_setup.hdf',H5F_ACC_RDWR_F,ffid,hdferr) 
      CALL io_gopen(ffid,"layers",fid,hdferr) 
      IF (io_dataexists(fid,"totalprojector")) THEN 
         !Projector already present                                     
         CALL io_gclose (fid, hdferr) 
         CALL io_hdfclose(ffid,hdferr) 
         RETURN 
      ENDIF 
                                                                        
      CALL io_createvar(fid,"totalprojector",H5T_NATIVE_DOUBLE,(/2      &
     &     ,SIZE(proj,1),SIZE(proj,2),SIZE(proj,3),SIZE(proj,4)/),did)  
      CALL io_WRITE(did,(/-1,1,1,1,1/),(/1,SIZE(proj,1),SIZE(proj,2)    &
     &     ,SIZE(proj,3),SIZE(proj,4)/),proj)                           
                                                                        
      CALL io_dclose (did, hdferr) 
      CALL io_gclose (fid, hdferr) 
      CALL io_hdfclose(ffid,hdferr) 
                                                                        
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      END                                           
