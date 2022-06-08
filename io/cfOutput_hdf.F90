MODULE m_cfOutput_hdf
#ifdef CPP_HDF

   USE hdf5
   USE m_hdf_tools

   IMPLICIT NONE

   PUBLIC opencfFile, closecfFile, writeCFpot, writeCFcdn

   CONTAINS

   SUBROUTINE opencfFile(fileID, atoms, cell, inFilename, l_create)

      USE m_types_atoms
      USE m_types_cell
      USE m_juDFT

      TYPE(t_atoms),                INTENT(IN)  :: atoms
      TYPE(t_cell),                 INTENT(IN)  :: cell
      INTEGER(HID_T),               INTENT(OUT) :: fileID
      CHARACTER(len=:), OPTIONAL, ALLOCATABLE,   INTENT(IN)  :: inFilename
      LOGICAL, OPTIONAL,            INTENT(IN)  :: l_create

      INTEGER          :: version,numCDN, numPOT
      CHARACTER(len=:),ALLOCATABLE :: filename
      LOGICAL          :: l_exist
      LOGICAL          :: l_error,l_createIn
      INTEGER          :: hdfError
      INTEGER(HSIZE_T) :: dims(2)
      INTEGER          :: dimsInt(2)

      INTEGER(HID_T)   :: metaGroupID
      INTEGER(HID_T)   :: generalGroupID
      INTEGER(HID_T)   :: bravaisMatrixSpaceID,bravaisMatrixSetID

      l_createIn = .TRUE.
      IF(PRESENT(l_create)) l_createIn = l_create

      version = 1
      IF(PRESENT(inFilename)) THEN
         filename = inFilename
      ELSE
         filename = "CFdata.hdf"
      ENDIF

      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)

      IF(l_createIn) THEN
         IF(l_exist) THEN
            CALL system('rm '//TRIM(ADJUSTL(filename)))
         ENDIF

         CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

         CALL h5gcreate_f(fileID, '/meta', metaGroupID, hdfError)
         CALL io_write_attint0(metaGroupID,'version',version)

         CALL h5gclose_f(metaGroupID, hdfError)

         !How many potentials and charge densities are written out
         CALL h5gcreate_f(fileID, '/general', generalGroupID, hdfError)
         CALL io_write_attint0(generalGroupID,'numPOT',COUNT(atoms%l_outputCFpot(:)))
         CALL io_write_attint0(generalGroupID,'numCDN',COUNT(atoms%l_outputCFcdn(:)))

         !Write out the Bravais Matrix (important to keep track of phase differences for coefficients with m != 0)
         dims(:2)=(/3,3/)
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),bravaisMatrixSpaceID,hdfError)
         CALL h5dcreate_f(generalGroupID, "bravaisMatrix", H5T_NATIVE_DOUBLE, bravaisMatrixSpaceID, bravaisMatrixSetID, hdfError)
         CALL h5sclose_f(bravaisMatrixSpaceID,hdfError)
         CALL io_write_real2(bravaisMatrixSetID,(/1,1/),dimsInt(:2),"amat",cell%amat)
         CALL h5dclose_f(bravaisMatrixSetID, hdfError)

         CALL h5gclose_f(generalGroupID, hdfError)
      ELSE IF(l_exist) THEN
         !Only open file
         CALL h5fopen_f(filename, H5F_ACC_RDWR_F, fileID, hdfError, H5P_DEFAULT_F)
      ELSE
         CALL juDFT_error("File not found", calledby="opencfFile")
      ENDIF

   END SUBROUTINE opencfFile

   SUBROUTINE closecfFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE closecfFile

   SUBROUTINE writeCFpot(fileID, atoms,input,iType,vlm)

      USE m_types_atoms
      USE m_types_input
      USE m_juDFT

      INTEGER(HID_T),   INTENT(IN)  :: fileID
      TYPE(t_atoms),    INTENT(IN)  :: atoms
      TYPE(t_input),    INTENT(IN)  :: input
      INTEGER,          INTENT(IN)  :: iType
      COMPLEX,          INTENT(IN)  :: vlm(:,:,:)

      INTEGER(HID_T) :: potGroupID, vlmGroupID
      INTEGER(HID_T) :: rmeshDataSpaceID,rmeshDataSetID
      INTEGER(HID_T) :: vlmDataSpaceID,vlmDataSetID

      INTEGER(HSIZE_T)  :: dims(7)
      INTEGER           :: dimsInt(7)
      INTEGER           :: hdfError
      INTEGER           :: l,m,lm
      LOGICAL           :: l_exist
      CHARACTER(len=:), ALLOCATABLE  :: groupName

      groupName = '/pot-'//int2str(iType)

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF(l_exist) THEN
         CALL juDFT_error('Group already exists: '//groupName, calledby="writeCFpot")
      ENDIF

      CALL h5gcreate_f(fileID, groupName, potGroupID, hdfError)

      !Radial Mesh
      CALL io_write_attint0(potGroupID,'atomType',iType)
      CALL io_write_attreal0(potGroupID,'RMT',atoms%rmt(iType))
      dims(:1)=[atoms%jri(iType)]
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),rmeshDataSpaceID,hdfError)
      CALL h5dcreate_f(potGroupID, "rmesh", H5T_NATIVE_DOUBLE, rmeshDataSpaceID, rmeshDataSetID, hdfError)
      CALL h5sclose_f(rmeshDataSpaceID,hdfError)
      CALL io_write_real1(rmeshDataSetID,[1],dimsInt(:1),"rmsh",atoms%rmsh(:atoms%jri(iType),iType))
      CALL h5dclose_f(rmeshDataSetID, hdfError)

      DO l = 2, 6, 2
         DO m = -l,l
            lm = l * (l+1) + m + 1
            CALL h5gcreate_f(potGroupID, 'VKS.'//int2str(l)//'.'//int2str(m), vlmGroupID, hdfError)
            CALL io_write_attint0(vlmGroupID,'l',l)
            CALL io_write_attint0(vlmGroupID,'m',m)

            dims(:3)=[2,atoms%jri(iType),input%jspins]
            dimsInt=dims
            CALL h5screate_simple_f(3,dims(:3),vlmDataSpaceID,hdfError)
            CALL h5dcreate_f(vlmGroupID, "vlm", H5T_NATIVE_DOUBLE, vlmDataSpaceID, vlmDataSetID, hdfError)
            CALL h5sclose_f(vlmDataSpaceID,hdfError)
            CALL io_write_complex2(vlmDataSetID,[-1,1,1],dimsInt(:3),"vlm",vlm(:atoms%jri(iType),lm,:))
            CALL h5dclose_f(vlmDataSetID, hdfError)

            CALL h5gclose_f(vlmGroupID, hdfError)
         ENDDO
      ENDDO
      CALL h5gclose_f(potGroupID, hdfError)

   END SUBROUTINE writeCFpot

   SUBROUTINE writeCFcdn(fileID, atoms,iType, n4f)

      USE m_types_atoms
      USE m_types_input
      USE m_juDFT

      INTEGER(HID_T),   INTENT(IN)  :: fileID
      TYPE(t_atoms),    INTENT(IN)  :: atoms
      INTEGER,          INTENT(IN)  :: iType
      REAL,             INTENT(IN)  :: n4f(:)

      INTEGER(HID_T) :: cdnGroupID
      INTEGER(HID_T) :: rmeshDataSpaceID,rmeshDataSetID
      INTEGER(HID_T) :: cdnDataSpaceID,cdnDataSetID

      INTEGER(HSIZE_T)  :: dims(7)
      INTEGER           :: dimsInt(7)
      INTEGER           :: hdfError
      LOGICAL           :: l_exist
      CHARACTER(len=:),ALLOCATABLE  :: groupName

      groupName = '/cdn-'//int2str(iType)

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF(l_exist) THEN
         CALL juDFT_error('Group already exists: '//groupName, calledby="writeCFcdn")
      ENDIF

      CALL h5gcreate_f(fileID, groupName, cdnGroupID, hdfError)

      !Radial Mesh
      CALL io_write_attint0(cdnGroupID,'atomType',iType)
      CALL io_write_attreal0(cdnGroupID,'RMT',atoms%rmt(iType))
      dims(:1)=[atoms%jri(iType)]
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),rmeshDataSpaceID,hdfError)
      CALL h5dcreate_f(cdnGroupID, "rmesh", H5T_NATIVE_DOUBLE, rmeshDataSpaceID, rmeshDataSetID, hdfError)
      CALL h5sclose_f(rmeshDataSpaceID,hdfError)
      CALL io_write_real1(rmeshDataSetID,[1],dimsInt(:1),"rmsh",atoms%rmsh(:atoms%jri(iType),iType))
      CALL h5dclose_f(rmeshDataSetID, hdfError)

      dims(:1)=[atoms%jri(iType)]
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),cdnDataSpaceID,hdfError)
      CALL h5dcreate_f(cdnGroupID, "cdn", H5T_NATIVE_DOUBLE, cdnDataSpaceID, cdnDataSetID, hdfError)
      CALL h5sclose_f(cdnDataSpaceID,hdfError)
      CALL io_write_real1(cdnDataSetID,[1],dimsInt(:1),"n4f",n4f(:atoms%jri(iType)))
      CALL h5dclose_f(cdnDataSetID, hdfError)

      CALL h5gclose_f(cdnGroupID, hdfError)

   END SUBROUTINE writeCFcdn

#endif
END MODULE m_cfOutput_hdf
