!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_embdesc 
      use m_juDFT
      use m_juDFT 
      IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Module that deals with embedding-plane descriptors               
!  These descibe the embedding planes by giving the atoms that cut the  
!  curved surface, the distance to the embedding planes etc.            
!                                                                       
! gf_checkdesc: compares two descriptors                                
! gf_writedescriptor: write descriptor to embedding-pot file            
! gf_writedesc_setup: write two descripter to gf_setup file             
! gf_readdescriptor: read descriptor from embedding-pot file            
! gf_readdesc_setup: read two descripter from gf_setup file             
! gf_deallocEmbDesc: deallocates descriptor                             
!                 Daniel Wortmann, (last modified: 05-01-06)            
!-----------------------------------------------                        
      PRIVATE 
      ! array index: 1=in,2=out,3=all                                   
      !<-- Descriptor of the embedding plane                            
                                                                        
      TYPE t_embdesc 
         LOGICAL         :: valid 
                                                                 !Atoms 
         INTEGER         :: all_atoms,cut_atoms_in,cut_atoms_out 
                                                       !the plane from  
                                                       !both sides      
                                                                        
                                                      !Atomic positions,
         REAL,POINTER    :: atoms_pos(:,:,:) 
                                                       !                
         REAL,POINTER    :: atoms_rmt(:,:) 
         INTEGER,POINTER :: atoms_z(:,:) 
                                                                        
                                                       !2D-unit cell    
         REAL            :: ucell(2,2) 
                                                       !Translation vect
         REAL            :: dvec(3) 
         REAL            :: dist_aux 
      END TYPE 
                                                                        
      !>                                                                
                                                                        
      PUBLIC :: t_embdesc,gf_checkdesc,gf_writedescriptor               &
     &     ,gf_readdescriptor,gf_readdesc_setup,gf_writedesc_setup      &
     &     ,gf_deallocEmbDesc                                           
      CONTAINS 
      !<-- F:gf_checkdesc(d1,d2)                                        
                                                                        
      FUNCTION gf_checkdesc(d1,d2,amat,shift,invert) 
!-----------------------------------------------                        
!    Compares two embedding-plane descriptions                          
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_embdesc),INTENT(IN)  :: d1,d2 
      LOGICAL                     :: gf_checkdesc 
      REAL,INTENT(IN)             :: amat(:,:) 
      REAL   ,INTENT(IN),OPTIONAL :: shift(:) 
      LOGICAL,INTENT(IN),OPTIONAL :: invert 
      !>                                                                
      !<-- Locals                                                       
      LOGICAL             :: same 
      INTEGER             :: n,nn,nnn 
      REAL                :: s(3,25) 
      LOGICAL             :: flip 
      !>                                                                
      flip = .FALSE. 
      IF (PRESENT(invert)) flip = invert 
      IF (d1%cut_atoms_in<0.OR.d2%cut_atoms_in<0) THEN 
         WRITE(*,*)                                                     &
     &        'WARNING: embedding-plane setup could not be checked'     
         gf_checkdesc = .TRUE. 
         RETURN 
      ENDIF 
      gf_checkdesc = .TRUE. 
      !<--Check no of atoms                                             
      IF (flip) THEN 
         IF ((d1%cut_atoms_in /= d2%cut_atoms_out).OR.(d1%cut_atoms_out &
     &        /= d2%cut_atoms_in)) THEN                                 
            gf_checkdesc         =  .FALSE. 
            WRITE(6,*) "No of atoms cutting the plane:" 
            WRITE(6,*) "         D1:", d1%cut_atoms_in,                 &
     &           d1%cut_atoms_out                                       
            WRITE(6,*) "         D2:", d2%cut_atoms_in,                 &
     &           d2%cut_atoms_out                                       
         ENDIF 
      ELSE 
         IF (d1%cut_atoms_in /= d2%cut_atoms_in) THEN 
            gf_checkdesc         =  .FALSE. 
            WRITE(6,*) "No of atoms cutting the plane:" 
            WRITE(6,*) "         from within:", d1%cut_atoms_in,        &
     &           d2%cut_atoms_in                                        
         ENDIF 
         IF (d1%cut_atoms_out /= d2%cut_atoms_out) THEN 
            gf_checkdesc          =  .FALSE. 
            WRITE(6,*) "No of atoms cutting the plane:" 
            WRITE(6,*) "         from outside:", d1%cut_atoms_out,      &
     &           d2%cut_atoms_out                                       
         ENDIF 
      ENDIF 
                                                                        
      !>                                                                
      !<--Check setup                                                   
                                                                        
      IF (ABS(d1%dist_aux - d2%dist_aux) > 2.E-3 )THEN 
         WRITE(6,*) "Distance dist_aux:", d1%dist_aux,d2%dist_aux
         WRITE(*,*) "WARNING: dist_aux differs, might indicate wrong setup"
         WRITE(*,*) "Distance dist_aux:", d1%dist_aux,d2%dist_aux
      ENDIF 
      IF (ANY(ABS(d1%ucell - d2%ucell) > 1E-4).OR.                      &
     &    ANY(ABS(amat(1:2,1:2) - d2%ucell) > 1E-4) ) THEN              
         gf_checkdesc = .FALSE. 
         WRITE(6,*) "2DUnit cell:" 
         WRITE(6,"(2f15.6)") d1%ucell
         WRITE(6,"(2f15.6)") d2%ucell
         write(6,"(2f15.6)") amat(1:2,1:2)
      ENDIF 
                                                                        
      !>                                                                
      !<--MT-radii                                                      
      IF (flip) THEN 
         IF ((ANY(ABS(d1%atoms_rmt(:d1%cut_atoms_in,1) -                &
     &        d2%atoms_rmt(:d2%cut_atoms_out,2)) > 1E-4)).OR.            &
     &        (ANY(ABS(d1%atoms_rmt(:d1%cut_atoms_out,2)               &
     &        -d2%atoms_rmt(:d2%cut_atoms_in,1)) > 1E-4))) THEN         
            gf_checkdesc = .FALSE. 
            WRITE(6,*) "MT-radius of atoms differ" 
         ENDIF 
      ELSE 
         IF (ANY(ABS(d1%atoms_rmt(:d1%cut_atoms_in,1) -                 &
     &        d2%atoms_rmt(:d2%cut_atoms_in,1)) > 1E-4)) THEN           
            gf_checkdesc = .FALSE. 
            WRITE(6,*) "MT-radius of 'in'-atoms differ" 
            WRITE(6,*) d1%atoms_rmt(:d1%cut_atoms_in,1)                 &
     &           ,d2%atoms_rmt(:d2%cut_atoms_in,1)                      
         ENDIF 
         IF (ANY(ABS(d1%atoms_rmt(:d1%cut_atoms_out,2) -                &
     &        d2%atoms_rmt(:d2%cut_atoms_out,2)) > 1E-4)) THEN          
            gf_checkdesc = .FALSE. 
            WRITE(6,*) "MT-radius of 'out'-atoms differ" 
            WRITE(6,*) d1%atoms_rmt(:d1%cut_atoms_out,2)                &
     &           ,d2%atoms_rmt(:d2%cut_atoms_out,2)                     
         ENDIF 
      ENDIF 
      !>                                                                
      !<-- Z                                                            
                                                                        
      IF (flip) THEN 
         IF ((ANY(ABS(d1%atoms_z(:d1%cut_atoms_in,1) -                  &
     &        d2%atoms_z(:d2%cut_atoms_out,2)) > 1E-4)).OR.              &
     &        (ANY(ABS(d1%atoms_z(:d1%cut_atoms_out,2)                 &
     &        -d2%atoms_z(:d2%cut_atoms_in,1)) > 1E-4)))THEN            
            gf_checkdesc = .FALSE. 
            WRITE(6,*) "Z-Values of atoms differ" 
         ENDIF 
      ELSE 
         IF (ANY(ABS(d1%atoms_z(:d1%cut_atoms_in,1) -                   &
     &        d2%atoms_z(:d2%cut_atoms_in,1)) > 1E-4)) THEN             
            gf_checkdesc = .FALSE. 
            WRITE(6,*) "Z-Values of 'in'-atoms differ" 
            WRITE(6,*) d1%atoms_z(:d1%cut_atoms_in,1)                   &
     &           ,d2%atoms_z(:d2%cut_atoms_in,1)                        
         ENDIF 
         IF (ANY(ABS(d1%atoms_z(:d1%cut_atoms_out,2) -                  &
     &        d2%atoms_z(:d2%cut_atoms_out,2)) > 1E-4)) THEN            
            gf_checkdesc = .FALSE. 
            WRITE(6,*) "Z-Values of 'out'-atoms differ" 
            WRITE(6,*) d1%atoms_z(:d1%cut_atoms_out,2)                  &
     &           ,d2%atoms_z(:d2%cut_atoms_out,2)                       
         ENDIF 
      ENDIF 
                                                                        
      !>                                                                
      !<-- Positions                                                    
      nnn = 0 
      DO n =-2,2 
         DO nn =-2,2 
            nnn  = nnn+1 
            IF (PRESENT(shift)) THEN 
               s(:2,nnn) = MATMUL(amat(:2,:2),(/1.0*n,1.0*nn/)+shift) 
            ELSE 
               s(:2,nnn) = MATMUL(amat(:2,:2),(/1.0*n,1.0*nn/)) 
            ENDIF 
            s(3,nnn) = 0.0 
         ENDDO 
      ENDDO 
      IF (flip) THEN 
         DO n = 1,d1%cut_atoms_in 
            same = .FALSE. 
            DO nn = 1,d2%cut_atoms_out 
               DO nnn  = 1,SIZE(s,2) 
                  IF (ALL(ABS(d1%atoms_pos(:,n,1)-d2%atoms_pos(:,nn,2)  &
     &                 -s(:,nnn))<0.03)) same = .TRUE.                  
               ENDDO 
            ENDDO 
            IF (.NOT.same) THEN 
               WRITE(6,*) "Positions of atoms differ" 
               WRITE(6,*) d1%atoms_pos(:,n,1)                           &
     &              ,d2%atoms_pos(:,:d2%cut_atoms_out,2)                
               gf_checkdesc = .FALSE. 
            ENDIF 
         ENDDO 
         DO n = 1,d1%cut_atoms_out 
            same = .FALSE. 
            DO nn = 1,d2%cut_atoms_in 
               DO nnn  = 1,SIZE(s,2) 
                  IF (ALL(ABS(d1%atoms_pos(:,n,2)-d2%atoms_pos(:,nn,1)  &
     &                 -s(:,nnn)) <0.03)) same = .TRUE.                 
               ENDDO 
            ENDDO 
            IF (.NOT.same) THEN 
               WRITE(6,*) "Positions of atoms differ" 
               WRITE(6,*) d1%atoms_pos(:,n,2)                           &
     &              ,d2%atoms_pos(:,:d2%cut_atoms_in,1)                 
               gf_checkdesc = .FALSE. 
            ENDIF 
         ENDDO 
      ELSE 
         DO n = 1,d1%cut_atoms_in 
            same = .FALSE. 
            DO nn = 1,d2%cut_atoms_in 
               DO nnn  = 1,SIZE(s,2) 
                  IF (ALL(ABS(d1%atoms_pos(:,n,1)-d2%atoms_pos(:,nn,1)  &
     &                 -s(:,nnn))<0.03)) same = .TRUE.                  
               ENDDO 
            ENDDO 
            IF (.NOT.same) THEN 
               WRITE(6,*) "Positions of 'in'-atoms differ" 
               WRITE(6,*) d1%atoms_pos(:,n,1)                           &
     &              ,d2%atoms_pos(:,:d2%cut_atoms_in,1)                 
               gf_checkdesc = .FALSE. 
            ENDIF 
         ENDDO 
         DO n = 1,d1%cut_atoms_out 
            same = .FALSE. 
            DO nn = 1,d2%cut_atoms_out 
               DO nnn  = 1,SIZE(s,2) 
                  IF (ALL(ABS(d1%atoms_pos(:,n,2)-d2%atoms_pos(:,nn,2)  &
     &                 -s(:,nnn)) <0.03)) same = .TRUE.                 
               ENDDO 
            ENDDO 
            IF (.NOT.same) THEN 
               WRITE(6,*) "Positions of 'out'-atoms differ" 
               WRITE(6,*) d1%atoms_pos(:,n,2)                           &
     &              ,d2%atoms_pos(:,:d2%cut_atoms_out,2)                
               gf_checkdesc = .FALSE. 
            ENDIF 
         ENDDO 
      ENDIF 
      !>                                                                
                                                                        
      !<-- if not same, print list                                      
      IF (.NOT.gf_checkdesc) THEN 
         WRITE(6,"(a30,2x,a30)") "First embedding plane"                &
     &        ,"Second embedding plane"                                 
         WRITE(6,"(a20,f10.6,2x,f10.6)") "Thickness of Delta volume:"   &
     &        ,d1%dist_aux,d2%dist_aux                                  
         WRITE(6,"(a20,i10,2x,i10)") "Cut atoms(in)",d1%cut_atoms_in    &
     &        ,d2%cut_atoms_in                                          
         DO n = 1,MAX(d1%cut_atoms_in,d2%cut_atoms_in) 
            IF (n <= d1%cut_atoms_in) THEN 
               IF (n <= d2%cut_atoms_in) THEN 
                  WRITE(6,"(2(4(f5.3,1x),i6,2x))") d1%atoms_pos(:,n,1)  &
     &                 ,d1%atoms_rmt(n,1),d1%atoms_z(n,1),d2%atoms_pos(:&
     &                 ,n,1),d2%atoms_rmt(n,1),d2%atoms_z(n,1)          
               ELSE 
                  WRITE(6,("(2(4(f5.3,1x),i6,2x))")) d1%atoms_pos(:,n,1)&
     &                 ,d1%atoms_rmt(n,1),d1%atoms_z(n,1)               
               ENDIF 
            ELSE 
               IF (n <= d2%cut_atoms_in) THEN 
                  WRITE(6,"(32x,4(f5.3,1x),i6,2x)") d2%atoms_pos(:,n,1) &
     &                 ,d2%atoms_rmt(n,1),d2%atoms_z(n,1)               
               ENDIF 
            ENDIF 
         ENDDO 
                                                                        
         WRITE(6,"(a20,i10,2x,i10)") "Cut atoms(out)"                   &
     &        ,d1%cut_atoms_out,d2%cut_atoms_out                        
         DO n = 1,MAX(d1%cut_atoms_out,d2%cut_atoms_out) 
            IF (n <= d1%cut_atoms_out) THEN 
               IF (n <= d2%cut_atoms_out) THEN 
                  WRITE(6,"(2(4(f5.3,1x),i6,2x))") d1%atoms_pos(:,n,2)  &
     &                 ,d1%atoms_rmt(n,2),d1%atoms_z(n,2),d2%atoms_pos(:&
     &                 ,n,2),d2%atoms_rmt(n,2),d2%atoms_z(n,2)          
               ELSE 
                  WRITE(6,"(2(4(f5.3,1x),i6,2x))") d1%atoms_pos(:,n,2)  &
     &                 ,d1%atoms_rmt(n,2),d1%atoms_z(n,2)               
               ENDIF 
            ELSE 
               IF (n <= d2%cut_atoms_out) THEN 
                  WRITE(6,"(32x,4(f5.3,1x),i6,2x)") d2%atoms_pos(:,n,2) &
     &                 ,d2%atoms_rmt(n,2),d2%atoms_z(n,2)               
               ENDIF 
            ENDIF 
         ENDDO 
      ENDIF 
      !>                                                                
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
      !<-- S:gf_writedescriptor(region,side,desc)                       
      SUBROUTINE gf_writedescriptor(region,side,desc) 
!******************************************                             
!     Writes the descriptor of the embedding potential to the           
!     embedding-potential-file                                          
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_gf_io2dmat 
      USE m_gf_types 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)            :: region,side 
      TYPE(t_embdesc),INTENT(INOUT) :: desc 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER        :: hdferr 
      INTEGER(HID_T) :: fid,gid 
                                                                        
      !>                                                                
                                                                        
                                   !invalid descriptor->no write        
      IF (.NOT.desc%valid) RETURN 
      fid = gf_io2dmatfid(IO2D_EMB,region,side,IO2D_WRITE) 
                        !file not opened -> no write                    
      IF (fid<0) RETURN 
      !Check if the group embdesc exists                                
      IF (io_groupexists(fid,"embdesc")) THEN 
         CALL io_gdelete(fid,"embdesc",hdferr) 
      ENDIF 
      CALL io_gcreate(fid,"embdesc",gid,hdferr) 
      !CALL the priv_savedesc to do the IO                              
      CALL priv_savedesc(gid,desc) 
      CALL io_gclose(gid,hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<-- S:gf_readdescriptor(region,side,desc)                        
      SUBROUTINE gf_readdescriptor(region,side,amat,desc) 
!******************************************                             
!     reads the descriptor of the embedding potential                   
!     if an old embedding potential file with no descriptor is found    
!     set desc%cut_atoms_in =-1                                         
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_gf_types 
      USE m_hdf_tools 
      USE m_gf_io2dmat 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)             :: region,side 
      REAL                           :: amat(:,:) 
      TYPE(t_embdesc),INTENT(OUT)    :: desc 
      !>                                                                
      !<--Locals                                                        
      INTEGER        :: hdferr 
      INTEGER(HID_T) :: fid,gid 
      REAL           :: shift(2) 
      !>                                                                
      fid = gf_io2dmatfid(IO2D_EMB,region,side,IO2D_READ) 
      shift = gf_io2dmatshift(IO2D_EMB,region,side,IO2D_READ) 
      desc%valid = .FALSE. 
      IF (fid<0) THEN 
         RETURN 
      ENDIF 
      !Check if the group embdesc exists                                
      IF (.NOT.io_groupexists(fid,"embdesc")) THEN 
         RETURN 
      ENDIF 
      desc%valid = .TRUE. 
      CALL io_gopen(fid,"embdesc",gid,hdferr) 
      !rest is in priv_readdesc()                                       
      CALL priv_readdesc(gid,desc,shift,amat) 
      CALL io_gclose(gid,hdferr) 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_readdesc_setup(desc_right,desc_left)                    
      SUBROUTINE gf_readdesc_setup(layer,desc_right,desc_left) 
!******************************************                             
!     reads the descriptors for the two embedding planes                
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: layer 
      TYPE(t_embdesc),INTENT(OUT) :: desc_right,desc_left 
      !>                                                                
      !<--Locals                                                        
      INTEGER        :: hdferr 
      INTEGER(HID_T) :: fid,gid,gid_r,gid_l,ggid 
      !>                                                                
      CALL io_hdfopen("gf_setup.hdf",  H5F_ACC_RDONLY_F, fid, hdferr) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      !<-- if no descriptors are present return invalid descriptors     
                                                                        
      IF (.NOT.io_groupexists(ggid,"embdesc")) THEN 
         desc_right%valid = .FALSE. 
         desc_left%valid = .FALSE. 
         CALL io_gclose(ggid,hdferr) 
         CALL io_hdfclose(fid,hdferr) 
         RETURN 
      ENDIF 
      desc_right%valid = .TRUE. 
      desc_left%valid  = .TRUE. 
                                                                        
      !>                                                                
      CALL io_gopen(ggid,"embdesc",gid,hdferr) 
      CALL io_gopen(gid,"left",gid_l,hdferr) 
      CALL io_gopen(gid,"right",gid_r,hdferr) 
      CALL priv_readdesc(gid_l,desc_left) 
      CALL priv_readdesc(gid_r,desc_right) 
      CALL io_gclose(gid_l,hdferr) 
      CALL io_gclose(gid_r,hdferr) 
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_hdfclose(fid,hdferr) 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- S:gf_writedesc_setup(desc_right,desc_left)                   
      SUBROUTINE gf_writedesc_setup(layer,desc_right,desc_left) 
!******************************************                             
!     saves the descriptors for the two embedding planes                
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      USE m_hdf_tools 
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: layer 
      TYPE(t_embdesc),INTENT(INOUT) :: desc_right,desc_left 
      !>                                                                
      !<-- Locals                                                       
      INTEGER        :: hdferr 
      INTEGER(HID_T) :: fid,gid,gid_r,gid_l,ggid 
      !>                                                                
      CALL io_hdfopen("gf_setup.hdf",H5F_ACC_RDWR_F, fid, hdferr) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      !<-- if group emb_plane exists, stop!                             
                                                                        
      IF (io_groupexists(ggid,"embdesc")) THEN 
         CALL juDFT_error("Group embdesc does exist already in gf_setup") 
      ENDIF 
                                                                        
      !>                                                                
      CALL io_gcreate(ggid,"embdesc",gid,hdferr) 
      CALL io_gcreate(gid,"left",gid_l,hdferr) 
      CALL io_gcreate(gid,"right",gid_r,hdferr) 
      CALL priv_savedesc(gid_l,desc_left) 
      CALL priv_savedesc(gid_r,desc_right) 
      CALL io_gclose(gid_l,hdferr) 
      CALL io_gclose(gid_r,hdferr) 
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_hdfclose(fid,hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<-- S:gf_deallocEmbDesc(d)                                       
                                                                        
      SUBROUTINE gf_deallocEmbDesc(d) 
!-----------------------------------------------                        
!   deallocates all data from the descriptor                            
!           (last modified: 04-06-15) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_embdesc),INTENT(OUT) :: d 
      !>                                                                
                               !empty description                       
      IF (.NOT.d%valid) RETURN 
      !NAG does not like this....                                       
!      IF (d%cut_atoms_in+d%cut_atoms_out>0)                            
!     +     DEALLOCATE(d%atoms_pos,d%atoms_rmt,d%atoms_z)               
      d%valid = .FALSE. 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:priv_savedesc(gid,desc)                                    
                                                                        
      SUBROUTINE priv_savedesc(gid,desc) 
!******************************************                             
!     Writes the descriptor of the embedding potential to the           
!     hdf-location                                                      
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_gf_types 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)  :: gid 
      TYPE(t_embdesc),INTENT(INOUT) :: desc 
      !>                                                                
      CALL priv_desc_normalize(desc) 
      CALL io_write_att(gid,"atoms",(/desc%cut_atoms_in                 &
     &     ,desc%cut_atoms_out,desc%all_atoms/))                        
      CALL io_WRITE_att(gid,"ucell",desc%ucell) 
      CALL io_WRITE_att(gid,"dist_aux",desc%dist_aux) 
      CALL io_write_att(gid,"dvec",desc%dvec) 
      IF (SIZE(desc%atoms_pos)>0) THEN 
         CALL io_WRITE_att(gid,"atoms_pos",desc%atoms_pos) 
         CALL io_WRITE_att(gid,"atoms_rmt",desc%atoms_rmt) 
         CALL io_WRITE_att(gid,"atoms_z",desc%atoms_z) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:priv_readdesc(gid,desc)                                    
                                                                        
      SUBROUTINE priv_readdesc(gid,desc,shift,amat) 
!******************************************                             
!     reads the descriptor of the embedding plane from hdf-location     
!     if an old embedding potential file with no descriptor is found    
!     set desc%cut_atoms_in =-1                                         
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_gf_types 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN)    :: gid 
      TYPE(t_embdesc),INTENT(OUT)  :: desc 
      REAL,OPTIONAL   ,INTENT(IN)  :: shift(2),amat(:,:) 
      !>                                                                
      !<--Locals                                                        
      INTEGER        :: tmpint(3) 
      REAL           :: s(2) 
      !>                                                                
      CALL io_READ_att(gid,"atoms",tmpint) 
      desc%cut_atoms_in= tmpint(1) 
      desc%cut_atoms_out= tmpint(2) 
      desc%all_atoms = tmpint(3) 
      CALL io_READ_att(gid,"ucell",desc%ucell) 
      CALL io_READ_att(gid,"dist_aux",desc%dist_aux) 
      CALL io_READ_att(gid,"dvec",desc%dvec) 
      !Allocate Storage                                                 
      IF (MAXVAL(tmpint) == 0) THEN 
         ALLOCATE(desc%atoms_pos(1,1,1)) 
      ELSE 
         ALLOCATE(desc%atoms_pos(3,MAXVAL(tmpint),3)) 
         ALLOCATE(desc%atoms_rmt(MAXVAL(tmpint),3)) 
         ALLOCATE(desc%atoms_z(MAXVAL(tmpint),3)) 
         CALL io_READ_att(gid,"atoms_pos",desc%atoms_pos) 
         IF (PRESENT(shift)) THEN 
            s=matmul(amat(:2,:2),shift) 
            desc%atoms_pos(1,:,:) = desc%atoms_pos(1,:,:)+s(1) 
            desc%atoms_pos(2,:,:) = desc%atoms_pos(2,:,:)+s(2) 
         ENDIF 
         CALL io_READ_att(gid,"atoms_rmt",desc%atoms_rmt) 
         CALL io_READ_att(gid,"atoms_z",desc%atoms_z) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:priv_desc_normalize(d)                                     
                                                                        
      SUBROUTINE priv_desc_normalize(d) 
!-----------------------------------------------                        
!     Normalize a descriptor,i.e. change atomic positions               
!     to have coordinates within [0,1] in internal units                
!           (last modified: 08-06-24) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_math 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_embdesc),INTENT(INOUT) :: d 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER             :: n,idx(SIZE(d%atoms_z,1)) 
      REAL                :: a_inv(2,2),x(2) 
                                                                        
      !>                                                                
                                                                        
      a_inv = mat_inverse(d%ucell) 
      !<-- Shift "all" atoms                                            
      DO n  = 1,d%all_atoms 
         x  = MATMUL(a_inv,d%atoms_pos(1:2,n,3)) 
         x = FLOOR(x) 
         d%atoms_pos(1:2,n,3) = d%atoms_pos(1:2,n,3) -MATMUL(d%ucell,x) 
      ENDDO 
      !sort atoms                                                       
      IF (d%all_atoms>0) THEN 
         CALL sort2(d%atoms_pos(1,:d%all_atoms,3),d%atoms_pos(2         &
     &        ,:d%all_atoms,3),idx(:d%all_atoms))                       
         d%atoms_pos(:,:d%all_atoms,3) = d%atoms_pos(:,idx(:d%all_atoms)&
     &        ,3)                                                       
         d%atoms_rmt(:d%all_atoms,3) = d%atoms_rmt(idx(:d%all_atoms),3) 
         d%atoms_z(:d%all_atoms,3)      = d%atoms_z(idx(:d%all_atoms),3) 
      ENDIF 
      !>                                                                
      !<-- Shift "in" atoms                                             
      DO n  = 1,d%cut_atoms_in 
         x  = MATMUL(a_inv,d%atoms_pos(1:2,n,1)) 
         x = FLOOR(x) 
         d%atoms_pos(1:2,n,1) = d%atoms_pos(1:2,n,1) -MATMUL(d%ucell,x) 
      ENDDO 
      IF (d%cut_atoms_in>0) THEN 
         CALL sort2(d%atoms_pos(1,:d%cut_atoms_in,1),d%atoms_pos(2      &
     &        ,:d%cut_atoms_in,1),idx(:d%cut_atoms_in))                 
         d%atoms_pos(:,:d%cut_atoms_in,1) = d%atoms_pos(:               &
     &        ,idx(:d%cut_atoms_in),1)                                  
         d%atoms_rmt(:d%cut_atoms_in,1) =                               &
     &        d%atoms_rmt(idx(:d%cut_atoms_in),1)                       
         d%atoms_z(:d%cut_atoms_in,1) = d%atoms_z(idx(:d%cut_atoms_in),1&
     &        )                                                         
      ENDIF 
      !>                                                                
      !<-- Shift "out" atoms                                            
      DO n = 1,d%cut_atoms_out 
         x  = MATMUL(a_inv,d%atoms_pos(1:2,n,2)) 
         x  = FLOOR(x) 
         d%atoms_pos(1:2,n,2) = d%atoms_pos(1:2,n,2) -MATMUL(d%ucell,x) 
      ENDDO 
      IF (d%cut_atoms_out>0) THEN 
         CALL sort2(d%atoms_pos(1,:d%cut_atoms_out,2),d%atoms_pos(2     &
     &        ,:d%cut_atoms_out,2),idx(:d%cut_atoms_out))               
         d%atoms_pos(:,:d%cut_atoms_out,2) = d%atoms_pos(:              &
     &        ,idx(:d%cut_atoms_out),2)                                 
         d%atoms_rmt(:d%cut_atoms_out,2) =                              &
     &        d%atoms_rmt(idx(:d%cut_atoms_out),2)                      
         d%atoms_z(:d%cut_atoms_out,2) = d%atoms_z(idx(:d%cut_atoms_out)&
     &        ,2)                                                       
      ENDIF 
      !>                                                                
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
