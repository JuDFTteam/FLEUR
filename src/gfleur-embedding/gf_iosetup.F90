!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_iosetup 
      IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Subroutines for the IO of setup-data from and to gf_setup.hdf    
!                 Daniel Wortmann, (08-02-06)                           
!-----------------------------------------------                        
      PRIVATE 
      PUBLIC gf_READ_layerinfo,gf_readatt,gf_writeatt,gf_link_layer 
      CONTAINS 
                                                                        
      !<-- S: gf_read_layerinfo(layers)                                 
                                                                        
      SUBROUTINE gf_READ_layerinfo(                                     &
     &     layers)                                                      
!***********************************************************************
!     Read the group 'layers' in the gf_pot.hdf file.                   
!     This group contains the information of the segmentation of the    
!     entire system into smaller subsystems.                            
!     Frank Freimuth, November 2007                                     
!***********************************************************************
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      TYPE(t_layers),INTENT(OUT) :: layers 
                                                                        
!     To make the code shorter                                          
      INTEGER::n 
!     the hdf-IDs                                                       
      INTEGER(HID_T)::fid,gid 
      INTEGER::hdferr 
      CALL io_hdfopen ('gf_setup.hdf', H5F_ACC_RDONLY_F, fid, hdferr) 
      CALL io_check('gf_iodop%gf_read_layerinfo:',hdferr) 
                                                                        
      CALL io_gopen(fid,'layers',gid,hdferr) 
      CALL io_check('gf_iodop%gf_read_layerinfo:open-group-layers',     &
     &                                                     hdferr)      
      CALL io_read_att(gid,"num_layers",layers%num_layers) 
      ALLOCATE(layers%c(layers%num_layers)) 
      ALLOCATE(layers%d(layers%num_layers)) 
      ALLOCATE(layers%dt(layers%num_layers)) 
      CALL io_read_att(gid,"d",layers%d) 
      CALL io_read_att(gid,"dt",layers%dt) 
      CALL io_read_att(gid,"c",layers%c) 
      CALL io_gclose(gid,hdferr) 
      CALL io_check('gf_iodop%gf_read_layerinfo: reading layers',hdferr) 
                                                                        
      !close file                                                       
      CALL io_hdfclose(fid,hdferr) 
      CALL io_check('gf_iodop$gf_read_layerinfo: close',hdferr) 
                                                                        
                                                                        
      END SUBROUTINE gf_read_layerinfo 
                                                                        
      !>                                                                
      !<-- S:gf_readatt(jspins,c1,c2,d,dt,atoms,cell,sym)               
                                                                        
      SUBROUTINE gf_readatt(layer,jspins,atoms,                         &
     &                      cell,sym)                                   
                                                                        
!***********************************************************************
!     This subroutine reads in the attributes from the potential file.  
!                                                                       
!                                                                       
!                                  Daniel Wortmann, Juelich, 2003       
!***********************************************************************
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      IMPLICIT NONE 
      INTEGER,INTENT(IN) :: layer 
                                                                        
      INTEGER,       INTENT(OUT) ::jspins 
      TYPE(t_atoms), INTENT(OUT) ::atoms 
      TYPE(t_cell), INTENT(OUT) ::cell 
      TYPE(t_sym), INTENT(OUT)  ::sym 
                                                                        
      REAL                      ::c1,c2,d,dt 
!     To make the code shorter                                          
      INTEGER::n 
!     the hdf-IDs                                                       
      INTEGER(HID_T)::fid,gid,did,ffid 
      INTEGER::hdferr 
                                                                        
      CALL io_hdfopen ('gf_setup.hdf', H5F_ACC_RDONLY_F, ffid, hdferr) 
      IF(hdferr /= 0) THEN 
         WRITE(*,*) 'Failed to open gf_setup.hdf' 
         CALL juDFT_error("No setup found") 
      ENDIF 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
                                                                        
!     READ d,dt                                                         
                                                                        
      CALL io_gopen (fid,'misc', gid, hdferr) 
      CALL io_check('gf_iodop%gf_readatt:open-group-misc',hdferr) 
      !Get No of spins                                                  
      CALL io_read_att(gid, "jspins" ,jspins) 
      CALL io_gclose (gid, hdferr) 
      CALL io_check('gf_iodop%gf_readatt:reading misc',hdferr) 
!                                                                       
!-->  read atoms from file                                              
!                                                                       
      CALL io_gopen (fid,'atoms', gid, hdferr) 
      CALL io_check('gf_iodop%gf_readatt:open-group-atoms',hdferr) 
      CALL io_read_att(gid, "ntype" ,atoms%ntype) 
      CALL io_read_att(gid, "nat" ,atoms%nat) 
      !Allocate Memory!                                                 
      n = atoms%ntype 
      ALLOCATE(atoms%nz(n),atoms%neq(n),atoms%jri(n),atoms%ncst(n)      &
     &     ,atoms%rmt(n),atoms%dx(n),atoms%volmts(n),atoms%zatom(n))    
      n=atoms%nat 
      ALLOCATE(atoms%ngopr(n),atoms%ntypsy(n),atoms%nlhtyp(n)           &
     &     ,atoms%invsat(n),atoms%pos(3,n),atoms%taual(3,n))            
      !read int_arrays                                                  
      CALL io_read_att(gid,"nz",atoms%nz) 
      CALL io_read_att(gid,"neq",atoms%neq) 
      CALL io_read_att(gid,"jri",atoms%jri) 
      CALL io_read_att(gid,"ncst",atoms%ncst) 
      CALL io_read_att(gid,"ngopr",atoms%ngopr) 
      CALL io_read_att(gid,"ntypsy",atoms%ntypsy) 
      CALL io_read_att(gid,"invsat",atoms%invsat) 
      !read real_arrays                                                 
      CALL io_READ_att(gid,"rmt",atoms%rmt) 
      CALL io_read_att(gid,"dx",atoms%dx) 
      CALL io_read_att(gid,"volmts",atoms%volmts) 
      CALL io_read_att(gid,"zatom",atoms%zatom) 
      CALL io_read_att(gid,"pos",atoms%pos) 
      CALL io_read_att(gid,"taual",atoms%taual) 
      ALLOCATE(atoms%rmsh(maxval(atoms%jri),atoms%ntype)) 
      !in newer files rmsh is no longer an attribute but a dataset      
      IF (IO_attexists(gid,"rmsh")) THEN 
         CALL io_READ_att(gid,"rmsh",atoms%rmsh) 
      ELSE 
         CALL io_dopen(gid,"rmsh",did,hdferr) 
         CALL io_READ(did,(/1,1/),(/SIZE(atoms%rmsh,1),atoms%ntype/)    &
     &        ,atoms%rmsh)                                              
         CALL io_dclose(did,hdferr) 
      ENDIF 
      CALL io_gclose (gid, hdferr) 
      CALL io_check('gf_iodop$gf_readatt:read atoms',hdferr) 
!                                                                       
!-->  read cell from file                                               
!                                                                       
      CALL io_gopen (fid,'cell', gid, hdferr) 
      CALL io_check('gf_iodop$gf_readatt:open-group cell',hdferr) 
      !read scalars                                                     
      CALL io_read_att(gid,"omtil",cell%omtil) 
      CALL io_read_att(gid,"area",cell%area) 
      CALL io_read_att(gid,"z1",cell%z1) 
      CALL io_read_att(gid,"vol",cell%vol) 
      CALL io_read_att(gid,"volint",cell%volint) 
      !read matrices                                                    
      CALL io_read_att(gid,"amat",cell%amat) 
      CALL io_read_att(gid,"bmat",cell%bmat) 
      CALL io_read_att(gid,"bbmat",cell%bbmat) 
      !read character                                                   
      CALL io_read_att(gid,"latnam",cell%latnam) 
      CALL io_gclose (gid, hdferr) 
      CALL io_check('gf_iodop$gf_readatt:read cell',hdferr) 
!                                                                       
!-->  read sym from file                                                
!                                                                       
      CALL io_gopen (fid,'sym', gid, hdferr) 
      CALL io_check('gf_iodop$gf_readatt:open-group sym',hdferr) 
      !read scalars                                                     
      CALL io_read_att(gid,"symor",sym%symor) 
      CALL io_read_att(gid,"invs2",sym%invs2) 
      CALL io_read_att(gid,"invs",sym%invs) 
      CALL io_read_att(gid,"zrfs",sym%zrfs) 
      CALL io_read_att(gid,"nop",sym%nop) 
      CALL io_read_att(gid,"nop2",sym%nop2) 
      !read matrices                                                    
      ALLOCATE(sym%mrot(3,3,sym%nop),sym%invtab(sym%nop),sym%tau(3      &
     &     ,sym%nop))                                                   
      CALL io_read_att(gid,"mrot",sym%mrot) 
      CALL io_read_att(gid,"invtab",sym%invtab) 
      CALL io_read_att(gid,"tau",sym%tau) 
      !read character                                                   
      CALL io_read_att(gid,"latnam",sym%latnam) 
      CALL io_read_att(gid,"namgrp",sym%namgrp) 
      CALL io_check('gf_iodop$gf_readatt:read sym',hdferr) 
      CALL io_gclose (gid, hdferr) 
      !close file                                                       
      CALL io_gclose(fid,hdferr) 
      CALL io_check('gf_iodop$gf_readatt:close layer',hdferr) 
      CALL io_hdfclose(ffid,hdferr) 
      CALL io_check('gf_iodop$gf_readatt:close file',hdferr) 
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_writeatt(layer,jspins,d,dt,atoms,cell,sym)              
                                                                        
      SUBROUTINE gf_writeatt(layer,jspins,atoms,cell,sym) 
!***********************************************************************
!     This subroutine writes the attributes to the setup file.          
!                                                                       
!                                                                       
!                                  Daniel Wortmann, Juelich, 2003,2008  
!***********************************************************************
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      INTEGER,       INTENT(IN)::jspins,layer 
      TYPE(t_atoms), INTENT(IN)::atoms 
      TYPE(t_cell), INTENT(IN) ::cell 
      TYPE(t_sym), INTENT(IN)  ::sym 
                                                                        
!     the hdf-IDs                                                       
      INTEGER(HID_T) ::fid,gid,did,ffid 
      INTEGER::hdferr 
                                                                        
      CALL io_hdfopen ('gf_setup.hdf', H5F_ACC_RDWR_F, ffid, hdferr) 
      IF(hdferr/=0) THEN 
         WRITE(*,*) 'opening file failed:gf_setup.hdf' 
         WRITE(*,*) 'Will create new file' 
         CALL h5fcreate_f('setup.hdf',H5F_ACC_TRUNC_F,ffid,hdferr) 
      ENDIF 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      CALL io_check('gf_iodop$gf_write:open ',hdferr) 
      !write misc-data                                                  
      CALL io_gcreate (fid,'misc', gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:open-group misc',hdferr) 
      CALL io_write_att(gid, "jspins" ,jspins) 
      CALL io_gclose (gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:write-group misc',hdferr) 
                                                                        
                                                                        
!                                                                       
!-->  put atoms from file                                               
!                                                                       
      CALL io_gcreate (fid,'atoms', gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:open-group atoms',hdferr) 
      !Integers for atoms                                               
      CALL io_write_att(gid, "ntype" ,atoms%ntype) 
      CALL io_write_att(gid, "nat" ,atoms%nat) 
      !One dimensional int-arrays                                       
      CALL io_write_att(gid,"nz",atoms%nz) 
      CALL io_write_att(gid,"neq",atoms%neq) 
      CALL io_write_att(gid,"jri",atoms%jri) 
      CALL io_write_att(gid,"ncst",atoms%ncst) 
      CALL io_write_att(gid,"ngopr",atoms%ngopr) 
      CALL io_write_att(gid,"ntypsy",atoms%ntypsy) 
      CALL io_write_att(gid,"invsat",atoms%invsat) 
      !real_arrays                                                      
      CALL io_write_att(gid,"rmt",atoms%rmt) 
      CALL io_write_att(gid,"dx",atoms%dx) 
      CALL io_write_att(gid,"volmts",atoms%volmts) 
      CALL io_write_att(gid,"zatom",atoms%zatom) 
      CALL io_write_att(gid,"pos",atoms%pos) 
      CALL io_write_att(gid,"taual",atoms%taual) 
      !CALL io_write_att(gid,"rmsh",atoms%rmsh)                         
      !rmsh might be to large for attribute                             
      CALL IO_createvar(gid,"rmsh",H5T_NATIVE_DOUBLE,(/SIZE(atoms%rmsh,1&
     &     ),SIZE(atoms%rmsh,2)/),did)                                  
      CALL IO_write(did,(/1,1/),(/size(atoms%rmsh,1)                    &
     &     ,SIZE(atoms%rmsh,2)/),atoms%rmsh)                            
      CALL io_dclose(did,hdferr) 
      CALL io_gclose (gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:write-group atoms',hdferr) 
!                                                                       
!-->  put cell to file                                                  
!                                                                       
      CALL io_gcreate (fid,'cell', gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:open-group cell',hdferr) 
      ! scalars                                                         
      CALL io_write_att(gid,"omtil",(/cell%omtil/)) 
      CALL io_write_att(gid,"area",cell%area) 
      CALL io_write_att(gid,"z1",cell%z1) 
      CALL io_write_att(gid,"vol",cell%vol) 
      CALL io_write_att(gid,"volint",cell%volint) 
      ! matrices                                                        
      CALL io_write_att(gid,"amat",cell%amat) 
      CALL io_write_att(gid,"bmat",cell%bmat) 
      CALL io_write_att(gid,"bbmat",cell%bbmat) 
      ! character                                                       
      CALL io_write_att(gid,"latnam",cell%latnam) 
      CALL io_gclose (gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:write-group cell',hdferr) 
!                                                                       
!-->  put sym to file                                                   
!                                                                       
      CALL io_gcreate (fid,'sym', gid, hdferr) 
      CALL io_check('gf_iodop$gf_write:open-group sym',hdferr) 
      !scalars                                                          
      CALL io_write_att(gid,"symor",sym%symor) 
      CALL io_write_att(gid,"invs2",sym%invs2) 
      CALL io_write_att(gid,"invs",sym%invs) 
      CALL io_write_att(gid,"zrfs",sym%zrfs) 
      CALL io_write_att(gid,"nop",sym%nop) 
      CALL io_write_att(gid,"nop2",sym%nop2) 
      !matrices                                                         
      CALL io_write_att(gid,"mrot",sym%mrot) 
      CALL io_write_att(gid,"invtab",sym%invtab) 
      CALL io_write_att(gid,"tau",sym%tau) 
      !character                                                        
      CALL io_write_att(gid,"latnam",sym%latnam) 
      CALL io_write_att(gid,"namgrp",sym%namgrp) 
      CALL io_gclose (gid, hdferr) 
      CALL io_gclose (fid, hdferr) 
      CALL io_check('gf_iodop$gf_write:write-group sym',hdferr) 
      !close file                                                       
      CALL io_hdfclose(ffid,hdferr) 
      CALL io_check('gf_iodop$gf_write:close',hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_link_layer(layer,orig_layer)                            
      SUBROUTINE gf_link_layer(layer,orig_layer) 
!-----------------------------------------------                        
!  creates a link to an existing layer in gf_setup file                 
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)     :: layer,orig_layer 
      !>                                                                
      !<-- Locals                                                       
      INTEGER (HID_T)       :: ffid 
      INTEGER               :: hdferr 
                                                                        
      !>                                                                
      CALL io_hdfopen ('gf_setup.hdf', H5F_ACC_RDWR_F, ffid, hdferr) 
      IF (hdferr /= 0) CALL juDFT_error                                   &
     &     ("gf_setup.hdf not found in gf_link_layer")                  
                                                                        
      IF (.NOT.io_groupexists(ffid,io_layername(orig_layer))) CALL      &
     &     juDFT_error(" Original layer not found in gf_link_layer")      
                                                                        
      CALL h5glink_f(ffid,H5G_LINK_HARD_F,io_layername(orig_layer)      &
     &     ,io_layername(layer),hdferr)                                 
                                                                        
      IF (hdferr /= 0) CALL juDFT_error                                   &
     &     ("linking failed in gf_link_layer")                          
                                                                        
      CALL io_hdfclose(ffid,hdferr) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
                                                                        
      END                                           
