MODULE m_greensf_io

#ifdef CPP_HDF

   USE hdf5
   USE m_hdf_tools

   IMPLICIT NONE

   PUBLIC openGreensFFile, closeGreensFFile, writeGreensFData

   CONTAINS

   SUBROUTINE openGreensFFile(fileID, input, gfinp, atoms, greensf, inFilename)

      USE m_types
      USE m_cdn_io

      TYPE(t_input),                INTENT(IN)  :: input
      TYPE(t_gfinp),                INTENT(IN)  :: gfinp
      TYPE(t_atoms),                INTENT(IN)  :: atoms
      TYPE(t_greensf),              INTENT(IN)  :: greensf
      CHARACTER(len=*), OPTIONAL,   INTENT(IN)  :: inFilename
      INTEGER(HID_T),               INTENT(OUT) :: fileID

      LOGICAL           :: l_exist
      CHARACTER(LEN=30) :: filename
      INTEGER(HID_T)    :: metaGroupID
      INTEGER(HID_T)    :: generalGroupID
      INTEGER(HID_T)    :: energyContourGroupID
      INTEGER(HID_T)    :: energyPointsSpaceID, energyPointsSetID
      INTEGER(HID_T)    :: energyWeightsSpaceID, energyWeightsSetID

      LOGICAL           :: l_error
      INTEGER           :: hdfError
      INTEGER           :: version
      INTEGER           :: dimsInt(7)
      REAL              :: eFermiPrev
      INTEGER(HID_T)    :: dims(7)

      version = 1
      IF(PRESENT(filename_in)) THEN
         filename = TRIM(ADJUSTL(inFilename))
      ELSE
         filename = "greensf.hdf"
      ENDIF

      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
      IF(l_exist) THEN
         CALL system('rm '//TRIM(ADJUSTL(filename)))
      ENDIF

      CALL h5fcreate_f(TRIM(ADJUSTL(filename)), H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

      CALL h5gcreate_f(fileID, '/meta', metaGroupID, hdfError)
      CALL io_write_attint0(metaGroupID,'version',version)
      CALL h5gclose_f(metaGroupID, hdfError)

      CALL readPrevEFermi(eFermiPrev,l_error)
      IF(l_error) THEN
         ! No previous eFermi available
         eFermiPrev = 0.0
      END IF

      CALL h5gcreate_f(fileID, '/general', generalGroupID, hdfError)
      CALL io_write_attint0(generalGroupID,'spins',input%jspins)
      CALL io_write_attreal0(generalGroupID,'FermiEnergy',eFermiPrev)
      CALL io_write_attlog0(generalGroupID,'sphavg',gfinp%l_sphavg)
      CALL io_write_attlog0(generalGroupID,'mperp',gfinp%l_mperp)
      CALL h5gclose_f(generalGroupID, hdfError)

      !Write out the energy contour and integration weights
      CALL h5gcreate_f(fileID, '/energyContour', energyContourGroupID, hdfError)
      CALL io_write_attint0(energyContourGroupID,'nz',greensf%nz)
      CALL io_write_attint0(energyContourGroupID,'shape',gfinp%mode) !Replace with string description

      dims(:2)=[2,greensf%nz]
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),energyPointsSpaceID,hdfError)
      CALL h5dcreate_f(energyContourGroupID, "Points", H5T_NATIVE_DOUBLE, energyPointsSpaceID, energyPointsSetID, hdfError)
      CALL h5sclose_f(energyPointsSpaceID,hdfError)
      CALL io_write_complex1(energyPointsSetID,[-1,1],dimsInt(:2),greensf%e)
      CALL h5dclose_f(energyPointsSetID, hdfError)

      dims(:2)=[2,greensf%nz]
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),energyWeightsSpaceID,hdfError)
      CALL h5dcreate_f(energyContourGroupID, "Weights", H5T_NATIVE_DOUBLE, energyWeightsSpaceID, energyWeightsSetID, hdfError)
      CALL h5sclose_f(energyWeightsSpaceID,hdfError)
      CALL io_write_complex1(energyWeightsSetID,[-1,1],dimsInt(:2),greensf%de)
      CALL h5dclose_f(energyWeightsSetID, hdfError)

      CALL h5gclose_f(energyContourGroupID, hdfError)

   END SUBROUTINE openGreensFFile

   SUBROUTINE closeGreensFFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE closeGreensFFile

   SUBROUTINE writeGreensFData(fileID, input, gfinp, greensf, mmpmat,)

      USE m_types
      USE m_constants

      TYPE(t_input),       INTENT(IN)  :: input
      TYPE(t_gfinp),       INTENT(IN)  :: gfinp
      TYPE(t_greensf),     INTENT(IN)  :: greensf
      COMPLEX,             INTENT(IN)  :: mmpmat(-lmaxU_Const:,-lmaxU_Const:,:,:)

      INTEGER(HID_T),      INTENT(IN)  :: fileID

      INTEGER(HID_T)       :: elementsGroupID
      INTEGER(HID_T)       :: currentelementGroupID
      INTEGER(HID_T)       :: mmpmatSpaceID, mmpmatSetID
      INTEGER(HID_T)       :: sphavgDataSpaceID, sphavgDataSetID
      INTEGER(HID_T)       :: uuDataSpaceID, uuDataSetID
      INTEGER(HID_T)       :: udDataSpaceID, udDataSetID
      INTEGER(HID_T)       :: duDataSpaceID, duDataSetID
      INTEGER(HID_T)       :: ddDataSpaceID, ddDataSetID

      CHARACTER(len=30)    :: elementName
      INTEGER              :: hdfError
      INTEGER              :: dimsInt(7)
      INTEGER              :: i_gf,ispin,m
      INTEGER(HID_T)       :: dims(7)
      REAL                 :: trc(input%jspins)

      CALL h5gcreate_f(fileID, '/elements', elementsGroupID, hdfError)
      CALL io_write_attint0(elementsGroupID,'NumElements',gfinp%n)
      CALL io_write_attint0(elementsGroupID,'maxl',lmaxU_Const)

      DO i_gf = 1, gfinp%n
         WRITE(elementName,200) i_gf
200      FORMAT('element-',i0)

         CALL h5gcreate_f(elementsGroupID, elementName, currentelementGroupID, hdfError)
         CALL io_write_attint0(currentelementGroupID,"l",gfinp%elem(i_gf)%l)
         CALL io_write_attint0(currentelementGroupID,"lp",gfinp%elem(i_gf)%lp)
         CALL io_write_attint0(currentelementGroupID,"atomType",gfinp%elem(i_gf)%atomType)
         CALL io_write_attint0(currentelementGroupID,"atomTypep",gfinp%elem(i_gf)%atomTypep)

         !Trace of occupation matrix
         trc=0.0
         DO ispin = 1, input%jspins
            DO m = -gfinp%elem(i_gf)%l, gfinp%elem(i_gf)%l
               trc(ispin) = trc(ispin) + REAL(mmpmat(m,m,i_gf,ispin))
            ENDDO
         ENDDO
         CALL io_write_attreal0(currentelementGroupID,"SpinUpTrace",trc(1))
         IF(input%jspins.EQ.2) THEN
            CALL io_write_attreal0(currentelementGroupID,"SpinDownTrace",trc(2))
         ENDIF

         !Occupation matrix
         IF(gfinp%l_mperp) THEN
            dims(:4)=[2,2*lmaxU_Const+1,2*lmaxU_Const+1,3]
            dimsInt=dims
            CALL h5screate_simple_f(4,dims(:4),mmpmatSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "mmpmat", H5T_NATIVE_DOUBLE, mmpmatSpaceID, mmpmatSetID, hdfError)
            CALL h5sclose_f(mmpmatSpaceID,hdfError)
            CALL io_write_complex3(mmpmatSetID,[-1,1,1,1],dimsInt(:4),mmpmat(:,:,i_gf,:))
            CALL h5dclose_f(mmpmatSetID, hdfError)
         ELSE
            dims(:4)=[2,2*lmaxU_Const+1,2*lmaxU_Const+1,input%jspins]
            dimsInt=dims
            CALL h5screate_simple_f(4,dims(:4),mmpmatSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "mmpmat", H5T_NATIVE_DOUBLE, mmpmatSpaceID, mmpmatSetID, hdfError)
            CALL h5sclose_f(mmpmatSpaceID,hdfError)
            CALL io_write_complex3(mmpmatSetID,[-1,1,1,1],dimsInt(:4),mmpmat(:,:,i_gf,:input%jspins))
            CALL h5dclose_f(mmpmatSetID, hdfError)
         ENDIF

         !Spherically averaged greensfData
         dims(:6)=[2,greensf%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
         dimsInt=dims
         CALL h5screate_simple_f(6,dims(:6),sphavgDataSpaceID,hdfError)
         CALL h5dcreate_f(currentelementGroupID, "SphAvg", H5T_NATIVE_DOUBLE, sphavgDataSpaceID, sphavgDataSetID, hdfError)
         CALL h5sclose_f(sphavgDataSpaceID,hdfError)
         CALL io_write_complex5(sphavgDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf%gmmpmat(:,:,:,:,:,i_gf))
         CALL h5dclose_f(sphavgDataSetID, hdfError)

         IF(.NOT.gfinp%l_sphavg) THEN

            !uu
            dims(:6)=[2,greensf%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),uuDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "UU", H5T_NATIVE_DOUBLE, uuDataSpaceID, uuDataSetID, hdfError)
            CALL h5sclose_f(uuDataSpaceID,hdfError)
            CALL io_write_complex5(uuDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf%uu(:,:,:,:,:,i_gf))
            CALL h5dclose_f(uuDataSetID, hdfError)

            !ud
            dims(:6)=[2,greensf%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),udDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "UD", H5T_NATIVE_DOUBLE, udDataSpaceID, udDataSetID, hdfError)
            CALL h5sclose_f(udDataSpaceID,hdfError)
            CALL io_write_complex5(udDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf%ud(:,:,:,:,:,i_gf))
            CALL h5dclose_f(udDataSetID, hdfError)

            !du
            dims(:6)=[2,greensf%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),duDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "DU", H5T_NATIVE_DOUBLE, duDataSpaceID, duDataSetID, hdfError)
            CALL h5sclose_f(duDataSpaceID,hdfError)
            CALL io_write_complex5(duDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf%du(:,:,:,:,:,i_gf))
            CALL h5dclose_f(duDataSetID, hdfError)

            !dd
            dims(:6)=[2,greensf%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),ddDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "DD", H5T_NATIVE_DOUBLE, ddDataSpaceID, ddDataSetID, hdfError)
            CALL h5sclose_f(ddDataSpaceID,hdfError)
            CALL io_write_complex5(ddDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf%dd(:,:,:,:,:,i_gf))
            CALL h5dclose_f(ddDataSetID, hdfError)

            !TODO write radial functions
         ENDIF

         CALL h5gclose_f(currentelementGroupID, hdfError)
      ENDDO

      CALL h5gclose_f(elementsGroupID, hdfError)

   END SUBROUTINE writeGreensFData

#endif

END MODULE m_greensf_io