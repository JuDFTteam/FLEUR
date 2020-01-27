!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_hdf_tools4 
      USE hdf5
!-----------------------------------------------                        
!     major rewrite of hdf_tools                                        
!     this module contains various subroutines                          
!                                                                       
!-----------------------------------------------                        
      !<--Parameters...                                                 
      INTEGER,PARAMETER::rkind=SELECTED_REAL_KIND(12) 
      INTEGER,PARAMETER::ckind=SELECTED_REAL_KIND(12) 
      !>
      PRIVATE:: io_datadim_did,io_datadim_name


      INTERFACE io_datadim
      MODULE PROCEDURE io_datadim_did,io_datadim_name
      END INTERFACE

      CONTAINS

      SUBROUTINE  io_hdfopen(filename,access_mode,fid,hdferr,access_prp)
      USE hdf5
      IMPLICIT NONE
      character(len=*),intent(in)    :: filename
      INTEGER       ,INTENT(in)      :: access_mode
      INTEGER(HID_T),INTENT(out)     :: fid
      INTEGER,INTENT(OUT),optional   :: hdferr
      INTEGER(HID_T),INTENT(in),optional ::access_prp
      INTEGER:: err

#ifdef CPP_DEBUG
#ifdef CPP_HDFMPI
      include "mpif.h"
      integer:: irank
      call MPI_COMM_RANK (MPI_COMM_WORLD,irank,err)
      write(*,"('PE:',i3,' opened:',a20,' rw:',l1)") irank,filename,access_mode==H5F_ACC_RDWR_F
#else
      write(*,"('Opened:',a20,' rw:',1l)") filename,access_mode==H5F_ACC_RDWR_F
#endif
#endif


      CALL h5fopen_f (filename,access_Mode,fid,err,access_prp)

      IF (present(hdferr)) hdferr=err

      end subroutine

      subroutine io_hdfclose(fid,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(in)    :: fid
      INTEGER,INTENT(OUT),optional :: hdferr

      INTEGER::err
#ifdef CPP_DEBUG
#ifdef CPP_HDFMPI
      include "mpif.h"
      integer:: irank
      character(len=20)::filename
      integer(size_t)::flength
      call h5fget_name_f(fid,filename,flength,err)
      call MPI_COMM_RANK (MPI_COMM_WORLD,irank,err)
      write(*,"('PE:',i3,' closed:',a20)") irank,filename
#else
      character(len=20)::filename
      integer(size_t)::flength
      call H5Fget_name_f(fid,filename,flength,err)
      write(*,"('Closed:',a20)") filename
#endif
#endif
      call h5fclose_f(fid,err)
      IF (present(hdferr)) hdferr=err
      end subroutine

      SUBROUTINE IO_gopen(fid,name,gid,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(in)    :: fid
      character(len=*),intent(in)  :: name
      INTEGER(HID_T),INTENT(out)   :: gid
      INTEGER,INTENT(OUT),optional :: hdferr
      INTEGER::err

      call h5gopen_f(fid,name,gid,err)
      IF (present(hdferr)) hdferr=err
      end subroutine

      SUBROUTINE IO_gcreate(fid,name,gid,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(in)    :: fid
      character(len=*),intent(in)  :: name
      INTEGER(HID_T),INTENT(out)   :: gid
      INTEGER,INTENT(OUT),optional :: hdferr
      INTEGER::err

      call h5gcreate_f(fid,name,gid,err)
      IF (present(hdferr)) hdferr=err
      end subroutine

      SUBROUTINE IO_gdelete(fid,name,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(in)    :: fid
      character(len=*),intent(in)  :: name
      INTEGER,INTENT(OUT),optional :: hdferr
      INTEGER::err

      call h5gunlink_f(fid,name,err)
      IF (present(hdferr)) hdferr=err
      end subroutine


      SUBROUTINE IO_gclose(gid,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)   :: gid
      INTEGER,INTENT(OUT),optional :: hdferr
      INTEGER::err

      call h5gclose_f(gid,err)
      IF (present(hdferr)) hdferr=err
      end subroutine

      SUBROUTINE IO_dopen(fid,name,did,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(in)    :: fid
      character(len=*),intent(in)  :: name
      INTEGER(HID_T),INTENT(out)   :: did
      INTEGER,INTENT(OUT),optional :: hdferr
      INTEGER::err

      call h5dopen_f(fid,name,did,err)
      IF (present(hdferr)) hdferr=err
      end subroutine

      SUBROUTINE IO_dclose(did,hdferr)
      USE hdf5
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)   :: did
      INTEGER,INTENT(OUT),optional :: hdferr
      INTEGER::err

      call h5dclose_f(did,err)
      IF (present(hdferr)) hdferr=err
      end subroutine

      subroutine io_datadim_name(gid,name,dim)
      IMPLICIT NONE
      INTEGER(HID_T),INTENT(IN)::gid
      CHARACTER*(*),INTENT(IN)::name
      INTEGER,INTENT(OUT)::dim(:)

      INTEGER(HID_T)::did
      CALL io_dopen(gid,name,did)
      CALL io_datadim_did(did,dim)
      CALL io_dclose(did)
      END SUBROUTINE

      !<-- S: io_datadim(did,dim)                                       
      SUBROUTINE io_datadim_did(did,dim)
!-----------------------------------------------                        
!  determine the dimesions of a hdf-dataspace                           
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE hdf5 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(in) :: did 
      INTEGER,INTENT(OUT)       :: DIM(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T)      :: sid 
      INTEGER(hsize_T)    :: md(SIZE(dim)),d(SIZE(dim)) 
      INTEGER             :: n,hdferr 
                                                                        
      !>                                                                
      CALL h5dget_space_f(did, sid, hdferr) 
      CALL io_check("Invalid dataset-id in data_dimension",hdferr) 
      CALL h5sget_simple_extent_ndims_f(sid,n,hdferr) 
      CALL io_check("Invalid dataspace in data_dimension",hdferr) 
                                                                        
      IF (n>SIZE(dim)) CALL                                             &
     &     hdf_err("data_dimension called with too small array")        
      d = 0 
      call h5sget_simple_extent_dims_f(sid, d, md, hdferr) 
      CALL io_check("Invalid dataspace in data_dimension",hdferr) 
      dim = d 
      CALL h5sclose_f(sid, hdferr) 
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- F: io_layername(layer)                                       
      FUNCTION io_layername(layer) 
!-----------------------------------------------                        
!  return string for layername                                          
!             (last modified: 07-11-08) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: layer 
      CHARACTER(len = 10)       ::io_layername 
      !>                                                                

                                                                        
      WRITE(io_layername,"(a6,i0)") "layer-",layer
                                                                        
      END FUNCTION 
      !>                                                                
      !<-- init,close library                                           
      SUBROUTINE hdf_init() 
!*****************************************************************      
! DESC:Opens library    No longer needed?!                              
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
      INTEGER :: hdferr 
      CALL h5open_f(hdferr) 
      !Turn automatic error checking on!                                
      CALL h5eset_auto_f(1,hdferr) 
      !CALL h5init_types(hdferr)                                        
#ifndef CPP_AIX                                                         
      CALL checklib() 
#endif                                                                  
      END SUBROUTINE 
      SUBROUTINE hdf_close() 
!*****************************************************************      
! DESC:Closes library   No longer needed?!                              
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
      INTEGER::hdferr 
      !CALL h5close_types(hdferr)                                       
      CALL h5close_f(hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<-- create a variable                                            
!*****************************************************************      
!                                                                       
!     The following subroutines create a var                            
!                                                                       
!                                                                       
!************************************************************           
      SUBROUTINE io_createvar(did,name,TYPE,dims,vid,chunk,fill)
!*****************************************************************      
      USE hdf5 
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(IN)  :: did 
      CHARACTER,INTENT(IN)       :: name*(*) 
      INTEGER(HID_T),INTENT(IN)  ::TYPE 
      INTEGER,INTENT(IN)         ::dims(:) 
      INTEGER(HID_T),INTENT(OUT) :: vid 
      INTEGER,INTENT(IN),OPTIONAL :: chunk 
      logical,intent(in),optional :: fill
                                                                        
      !locals                                                           
      INTEGER(HSIZE_T) ::DIM(SIZE(dims)) 
      INTEGER(HSIZE_T) :: chunk_DIM(SIZE(dims)) 
      INTEGER(HID_t) ::spaceid,prp_id 
      INTEGER       ::hdferr 
      dim = dims 
      WHERE (dim == 0) 
         dim = 1 
      endwhere 

      CALL h5pcreate_f(H5P_DATASET_CREATE_F, prp_id, hdferr)
      IF (PRESENT(chunk)) THEN 
         chunk_dim = dim 
         IF (chunk<SIZE(chunk_dim)) chunk_DIM(chunk+1:) = 1 
         CALL h5pset_chunk_f(prp_id, SIZE(dims),chunk_dim, hdferr) 
      ENDIF
      if (present(fill)) then
        if (.not.fill) call h5pset_fill_time_f(prp_id, H5D_FILL_TIME_NEVER_F, hdferr)
      endif

      CALL h5screate_simple_f(SIZE(dim),dim ,spaceid,hdferr) 
      CALL h5dcreate_f(did,name,TYPE,spaceid,                           &
     &                 vid, hdferr,prp_id)                              
      call h5pclose_f(prp_id,hdferr)
      CALL h5sclose_f(spaceid,hdferr) 
      CALL io_check('io_createvar:'//name,hdferr) 
      END SUBROUTINE 
      !>                                                                
      !<-- F:gettransprop()RESULT(trans)                                
      FUNCTION gettransprop()RESULT(trans) 
!********************************************************************** 
!  local FUNCTION to get default transfer-property                      
!********************************************************************** 
      USE hdf5 
      IMPLICIT NONE 
      INTEGER(HID_T)::trans 
#ifdef CPP_HDFMPI
      INCLUDE 'mpif.h' 
      INTEGER::hdferr 
      LOGICAL::l_mpi
      CALL MPI_INITIALIZED(l_mpi,hdferr)
      IF (l_mpi) THEN
         CALL h5pcreate_f(H5P_DATASET_XFER_F, trans, hdferr) 
         CALL h5pset_dxpl_mpio_f(trans,H5FD_MPIO_INDEPENDENT_F,hdferr)
      ELSE
         trans=H5P_DEFAULT_f 
      ENDIF
#else                                                                   
      trans=H5P_DEFAULT_f 
#endif                                                                  
      END FUNCTION 
!---------------------------------------------------------------------- 
      !>                                                                
      !<-- S:cleartransprop(trans)                                      
      SUBROUTINE cleartransprop(trans) 
!********************************************************************** 
!  local FUNCTION to get default transfer-property                      
!********************************************************************** 
      USE hdf5 
      IMPLICIT NONE 
      INTEGER(HID_T),INTENT(INOUT)::trans 
#ifdef CPP_HDFMPI
      INTEGER::hdferr 
      INCLUDE 'mpif.h' 
      IF (trans==H5P_DEFAULT_f) RETURN 
      CALL h5pclose_f(trans,hdferr) 
#else                                                                   
                             !just to use trans                         
      IF (trans == 1) RETURN 
      RETURN 
#endif                                                                  
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:io_check(text,err)                                         
                                                                        
      SUBROUTINE io_check(text,err,oid) 
!***********************************************************************
!      SUBROUTINE to check IO for error                                 
!                                                                       
!                                  Daniel Wortmann, Juelich, 2002       
!***********************************************************************
      USE hdf5 
      IMPLICIT NONE 
      INTEGER,INTENT(INOUT)              :: err 
      CHARACTER*(*),OPTIONAL             :: text 
      INTEGER(hid_t),INTENT(IN),OPTIONAL :: oid 
                                                                        
      CHARACTER(len = 5)  :: pe ="    :" 
      CHARACTER(len = 500):: object_name 
      INTEGER (hid_t)     :: itype 
      INTEGER             :: hdferr 
      INTEGER (size_t)    :: n,nn 
#ifdef CPP_HDFMPI                                                          
      include 'mpif.h' 
      INTEGER             :: irank,nerr 
      LOGICAL             :: l_mpi
      CALL MPI_INITIALIZED(l_mpi,nerr)
      IF (l_mpi) THEN
         CALL MPI_COMM_rank(MPI_COMM_WORLD,irank,nerr) 
         WRITE(pe,"(i4,a)") irank,":" 
      ENDIF
#endif                                                                  
      n = 500 
      IF (err>=0) RETURN 
                                                                        
      CALL h5eprint_f(err) 
      IF (PRESENT(text)) THEN 
         WRITE(*,*) pe,'IO-Error detected in: ',text 
      ENDIF 
                                                                        
      !Try to get a name for the oid                                    
      IF (PRESENT(oid)) THEN 
         write(*,*) "Offending Object-ID:",oid 
!         CALl h5iget_name_f(oid, object_name,n,nn,Hdferr) 
         hdferr=1
         If (hdferr /= 0) THEN 
            write(*,*) "Name of OID could not be determined" 
         ELSE 
            WRITE(*,*) "Name of OID:",object_name 
            !call h5iget_type_f(oid, itype, hdferr) 
            IF (hdferr /= 0) THEN 
               write(*,*) "Type of OID could not be determined" 
            ELSE 
               IF (iTYPE == H5I_FILE_F)                                 &
     &              WRITE(*,*) "OID is a file"                          
               IF (iTYPE == H5I_group_F)                                &
     &              WRITE(*,*) "OID is a group"                         
               IF (iTYPE == H5I_datatype_F)                             &
     &              WRITE(*,*) "OID is a datatype"                      
               IF (iTYPE == H5I_dataset_F)                              &
     &              WRITE(*,*) "OID is a dataset"                       
               IF (iTYPE == H5I_attr_F)                                 &
     &              WRITE(*,*) "OID is a attribute"                     
            ENDIF 
         ENDIF 
         !CALL h5fget_name_f(oid, object_name, n, hdferr) 
         !IF (hdferr /= 0 ) THEN 
         !   WRITE(*,*) "No File found" 
         !ELSE 
         !   WRITE(*,*) "Filename:", object_name 
         !ENDIF 
      ENDIF 
      !try to generate a execption to get a traceback :-)               
                              !log of a negative real                   
      write(*,*) log(1.0*err) 
                                                                        
      CALL hdf_err('IO-Error in hdf_tools') 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:checklib()                                                 
      SUBROUTINE checklib() 
!***********************************************************************
!      SUBROUTINE to check the basic library functions                  
!                                                                       
!                                  Daniel Wortmann, Juelich, 2002       
!***********************************************************************
      USE hdf5 
      IMPLICIT NONE 
      INTEGER :: hdferr 
      INTEGER(HID_T)::fid,gid 
                                                                        
      !In most cases this check is not needed!                          
             !comment this line if you want to perform the check!       
      return 
                                                                        
      CALL h5fcreate_f("hdftest_tmp.hdf",H5F_ACC_TRUNC_F,fid,hdferr) 
      CALL h5gcreate_f(fid,"testgroup",gid,hdferr) 
      CALL h5gclose_f(gid,hdferr) 
      CALL h5fclose_f(fid,hdferr) 
                                                                        
      CALL h5fopen_f("hdftest_tmp.hdf",H5F_ACC_RDONLY_F,fid,hdferr) 
      CALL h5gopen_f(fid,"testgroup",gid,hdferr) 
      CALL h5gclose_f(gid,hdferr) 
      CALL h5fclose_f(fid,hdferr) 
                                                                        
      OPEN(99,file ="hdftest_tmp.hdf") 
      CLOSE(99) 
                                                                        
      WRITE(*,*) 'HDF library was initialized' 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<--S:hdf_err                                                     
      SUBROUTINE hdf_err(message) 
!-----------------------------------------------                        
!     Version for LINUX compiled with IFC                               
!             (last modified: 05-02-25) D. Wortmann                     
!-----------------------------------------------                        
      use m_juDFT_stop
      IMPLICIT NONE 
      !<-- Arguments                                                    
      CHARACTER*(*)        ::message 
                                                                        
      !>                                                                
      WRITE(*,*) "Error in HDF-io" 
      WRITE(*,*) message 
      call judft_error(message)
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
