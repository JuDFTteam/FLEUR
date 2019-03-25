MODULE m_iomatrix_hdf
  USE m_judft
  USE hdf5
  USE m_hdf_tools
  USE m_types_mat
  USE m_types_mpimat
  IMPLICIT NONE
  PRIVATE
  PUBLIC iomatrix_hdf_close,iomatrix_hdf_open,iomatrix_hdf_write,iomatrix_hdf_read
  
CONTAINS
  SUBROUTINE  iomatrix_hdf_read(mat,nrec,did)
    CLASS(t_Mat),INTENT(INOUT)  :: mat
    INTEGER,INTENT(IN)          :: nrec
    INTEGER(HID_t),INTENT(in)   :: did

    INTEGER::mpi_comm,dim(4)

   
    INTEGER(HID_t) :: memspace,fspace,trans
    INTEGER(hsize_t):: dims(4)
    INTEGER        :: err
    REAL,ALLOCATABLE :: dat(:,:,:,:)
    SELECT TYPE(mat)
    TYPE is (t_mpimat)
       mpi_comm=mat%blacsdata%mpi_com !Only information used from mat intent(in)
    CLASS default
       mpi_comm=0
    END SELECT
    CALL io_datadim(did,dim)
    CALL mat%init(DIM(1)==1,DIM(2),DIM(3),mpi_comm,.TRUE.)
    SELECT TYPE(mat)
    TYPE is (t_mat)
       IF (mat%l_real) THEN
          CALL io_read(did,(/1,1,1,nrec/),(/1,mat%matsize1,mat%matsize2,1/),mat%data_r)
       ELSE
          CALL io_read(did,(/-1,1,1,nrec/),(/1,mat%matsize1,mat%matsize2,1/),mat%data_c)
       END IF
    TYPE is (t_mpimat)
       ALLOCATE(dat(MERGE(1,2,mat%l_real),mat%matsize1,mat%matsize2,1))
       
       CALL h5dget_space_f(did,fspace,err)
       CALL priv_create_hyperslab_from_blacsdesc(mat%l_real,nrec,fspace,mat%blacsdata%blacs_desc)
       dims=SHAPE(dat)
       CALL h5screate_simple_f(4,dims,memspace,err)
       trans=gettransprop()
       CALL h5dread_f(did,H5T_NATIVE_DOUBLE,dat,dims,err,memspace,fspace,trans)
       CALL h5sclose_f(memspace,err)                                  
       CALL h5sclose_f(fspace,err)
       CALL cleartransprop(trans) 
       IF (mat%l_real) THEN 
          mat%data_r=dat(1,:,:,1)
       ELSE
         mat%data_c=CMPLX(dat(1,:,:,1),dat(2,:,:,1))
      ENDIF
    END SELECT
  END SUBROUTINE iomatrix_hdf_read

  SUBROUTINE  iomatrix_hdf_write(mat,rec,did)
    CLASS(t_Mat),INTENT(IN)  :: mat
    INTEGER,INTENT(IN)       :: rec
    INTEGER(HID_t),INTENT(in)::did
    
    INTEGER(HID_t) :: memspace,fspace,trans
    INTEGER(HSIZE_t):: dims(4)
    INTEGER :: err
    REAL,ALLOCATABLE :: dat(:,:,:,:)

    
    SELECT TYPE(mat)
    TYPE is (t_mat)
       IF (mat%l_real) THEN
          CALL io_write(did,(/1,1,1,rec/),(/1,mat%matsize1,mat%matsize2,1/),mat%data_r)
       ELSE
          CALL io_write(did,(/1,1,1,rec/),(/1,mat%matsize1,mat%matsize2,1/),REAL(mat%data_c))
          CALL io_write(did,(/2,1,1,rec/),(/1,mat%matsize1,mat%matsize2,1/),AIMAG(mat%data_c))
       END IF
    TYPE is (t_mpimat)
       ALLOCATE(dat(MERGE(1,2,mat%l_real),mat%matsize1,mat%matsize2,1))
       IF (mat%l_real) THEN 
          dat(1,:,:,1)=mat%data_r
       ELSE
         dat(1,:,:,1)=REAL(mat%data_c)
         dat(2,:,:,1)=REAL(mat%data_c)
      ENDIF
      CALL h5dget_space_f(did,fspace,err)
      CALL priv_create_hyperslab_from_blacsdesc(mat%l_real,rec,fspace,mat%blacsdata%blacs_desc)
      dims=SHAPE(dat)
      CALL h5screate_simple_f(4,dims,memspace,err)
      trans=gettransprop()
      CALL h5dwrite_f(did,H5T_NATIVE_DOUBLE,dat,dims,err,memspace,fspace,trans)
      CALL h5sclose_f(memspace,err)                                  
      CALL h5sclose_f(fspace,err)
      CALL cleartransprop(trans) 
   END SELECT

  END SUBROUTINE iomatrix_hdf_write

  SUBROUTINE iomatrix_hdf_close(fid,did)
    INTEGER(hid_t),INTENT(inout):: fid,did
    INTEGER:: err
    CALL h5dclose_f(did,err)
    CALL h5fclose_f(fid,err)
  END SUBROUTINE iomatrix_hdf_close

  SUBROUTINE iomatrix_hdf_open(l_real,matsize,no_rec,filename,fid,did)
    LOGICAL,INTENT(IN)          :: l_real
    INTEGER,INTENT(in)          :: matsize,no_rec
    CHARACTER(len=*),INTENT(in) :: filename
    INTEGER(hid_t),INTENT(out)  :: fid,did
    
    INTEGER :: dims(4),err
    LOGICAL :: l_exist
    INTEGER(HID_T)  :: access_prp
#if defined(CPP_HDFMPI) && defined(CPP_MPI)
    include 'mpif.h'
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, err)
    CALL h5pset_fapl_mpio_f(access_prp, MPI_COMM_WORLD, MPI_INFO_NULL,err)
#else
    access_prp=H5P_DEFAULT_f
#endif


    INQUIRE(file=filename//'.hdf',exist=l_exist)
    IF (l_exist) THEN
       CALL h5fopen_f(filename//'.hdf',H5F_ACC_RDWR_F,fid,err,access_prp)
    ELSE
       CALL h5fcreate_f(filename//'.hdf',H5F_ACC_TRUNC_F,fid,err,H5P_DEFAULT_f,access_prp)
    ENDIF
    IF (io_dataexists(fid,'Matrix')) THEN
       CALL h5dopen_f(fid,"Matrix", did, err)
    ELSE
       !Create data-space
       dims(1)  = MERGE(1,2,l_real)
       dims(2:3)= matsize
       dims(4)  = no_rec
       call io_createvar(fid,"Matrix",H5T_NATIVE_DOUBLE,dims,did)
    END IF
  END SUBROUTINE iomatrix_hdf_open


  SUBROUTINE priv_create_hyperslab_from_blacsdesc(l_real,nrec,sid,blacsdesc)
    LOGICAL,INTENT(IN) :: l_real
    INTEGER,INTENT(in) :: nrec,blacsdesc(9)
    INTEGER(hid_t),INTENT(in):: sid
    
    INTEGER(hsize_t):: start(4),COUNT(4),stride(4),bloc(4)
    INTEGER         :: nprow,npcol,myrow,mycol,block_row,block_col,matsize,blacs_ctxt,err
    LOGICAL         :: ok
    !For readability get data from blacsdesc
    blacs_ctxt=blacsdesc(2)
    block_row=blacsdesc(5)
    block_col=blacsdesc(6)
    matsize=blacsdesc(4)
    CALL blacs_gridinfo(blacs_ctxt,nprow,npcol,myrow,mycol)

    CALL h5Sselect_none_f(sid,err) !unselect all elements in dataspace
    !Select blocks of blacs-grid
    start=(/0,myrow*block_row,mycol*block_col,nrec-1/)
    count=(/1,FLOOR(1.*matsize/block_row)+1,FLOOR(1.*matsize/block_col)+1,1/)
    stride=(/1,nprow*block_row,npcol*block_col,1/)
    bloc=(/MERGE(1,2,l_real),block_row,block_col,1/)
    CALL h5sselect_hyperslab_f(sid,H5S_SELECT_OR_F,start,count,err,stride,bloc)
    CALL h5sselect_valid_f(sid,ok,err)
    IF (.NOT.ok) THEN
       CALL h5sget_simple_extent_dims_f(sid,start,stride,err)
       !Cut to actual sizes
       start=(/0,0,0,0/)
       count=(/MERGE(1,2,l_real),matsize,matsize,int(stride(4))/)
       CALL h5sselect_hyperslab_f(sid,H5S_SELECT_AND_F,start,count,err) 
       CALL h5sselect_valid_f(sid,ok,err)
       IF (.NOT.ok) CALL judft_error("Writing of matrix failed, BUG in parallel HDF-IO")
    ENDIF
  END SUBROUTINE priv_create_hyperslab_from_blacsdesc
END MODULE m_iomatrix_hdf
