!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_ioCoulBoundary 
          IMPLICIT NONE
      CONTAINS 
                                                                        
      !<-- S: gf_readCoulBoundary(layer,stars,side,v_bound)             
      SUBROUTINE gf_readCoulBoundary(layer,stars,side,v_bound) 
!-----------------------------------------------                        
!  read coulomb potential boundary values from embpot file              
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      USE m_gf_io2dmat 
      use m_juDFT 
      USE m_constants 
      USE m_gf_types 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)       :: layer,side 
      TYPE(t_stars),INTENT(IN) :: stars 
      COMPLEX,INTENT(OUT)      :: v_bound(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(hid_t)         :: fid,gid,did 
      INTEGER                :: hdferr,n 
      REAL,ALLOCATABLE       :: v_b(:,:) 
      REAL                   :: k(2),shift(2) 
      !>                                                                
                                                                        
                                                                        
      fid = gf_io2dmatFID(IO2D_EMB,layer,side,io2d_READ) 
                                                                        
      IF (fid < 0) CALL juDFT_error(                                      &
     &     "No embedding potential file for boundary values")           
      CALL h5gopen_f(fid,"CoulBoundary",gid,hdferr) 
      IF (hdferr /= 0)  CALL                                            &
     &     juDFT_error("No boundary values found")                        
      CALL io_read_att(gid,"n",n) 
      IF (n/=size(v_bound)) THEN 
         WRITE(6,*) "Different no of boundary values" 
         WRITE(6,*) "In file:",n 
         WRITE(6,*) "In setup:",size(v_bound) 
         n=min(n,size(v_bound)) 
      ENDIF 
      ALLOCATE(v_b(n,2)) 
      CALL h5dopen_f(gid,"v_coul_real",did,hdferr) 
      CALL io_read(did,(/1/),(/n/),v_b(:,1)) 
      CALL h5dclose_f(did,hdferr) 
      CALL h5dopen_f(gid,"v_coul_aimag",did,hdferr) 
      CALL io_read(did,(/1/),(/n/),v_b(:,2)) 
      CALL h5dclose_f(did,hdferr) 
      v_bound = 0.0 
      v_bound(:n) = CMPLX(v_b(:,1),v_b(:,2)) 
      CALL h5gclose_f(gid,hdferr) 
      DEALLOCATE(v_b) 
                                                                        
      shift = gf_io2dmatshift(IO2D_EMB,layer,side,io2d_READ) 
      IF (ANY(shift /= 0.0)) THEN 
         WRITE(*,*) "Shifted Boundary Potential" 
         DO n = 1,stars%nq2 
            k = 1.0*stars%kv2(:,n) 
            v_bound(n) = v_bound(n)*EXP(CMPLX(0.0,-2.*pimach()          &
     &           *dot_PRODUCT(k,shift)))                                
         ENDDO 
      ENDIF 
      END SUBROUTINE 
      !>                                                                
      !<-- S: gf_saveCoulBoundary(layer,side,v_bound)                   
      SUBROUTINE gf_saveCoulBoundary(layer,side,v_bound) 
!-----------------------------------------------                        
!  read coulomb potential boundary values from embpot file              
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      USE m_gf_io2dmat 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)     :: layer,side 
      COMPLEX,INTENT(IN)     :: v_bound(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(hid_t)         :: fid,gid 
      INTEGER                :: hdferr 
      !>                                                                
                                                                        
                                                                        
      fid = gf_io2dmatFID(IO2D_EMB,layer,side,io2d_write) 
                          !nothing to save                              
      IF (fid < 0) RETURN 
      CALL h5gcreate_f(fid,"CoulBoundary",gid,hdferr) 
      IF (hdferr /= 0)  CALL                                            &
     &     juDFT_error("Could not create group for boundary values")      
      CALL io_write_att(gid,"n",size(v_bound)) 
      CALL io_write_var(gid,"v_coul_real",real(v_bound)) 
      CALL io_write_var(gid,"v_coul_aimag",aimag(v_bound)) 
      CALL h5gclose_f(gid,hdferr) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
