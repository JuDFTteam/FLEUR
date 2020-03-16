MODULE m_greensf_io

#ifdef CPP_HDF

   USE hdf5
   USE m_hdf_tools

   IMPLICIT NONE

   PUBLIC openGreensFFile, closeGreensFFile, writeGreensFData

   PUBLIC GREENSF_GENERAL_CONST, GREENSF_HUBBARD_CONST

   !---------------
   ! Storage Types
   !-------------------------------
   ! GREENS_GENERAL_CONST => All Green's function elements are saved
   ! GREENS_HUBBARD_CONST => Only Hubbard elements are saved with selfenergy (if available)
   INTEGER,    PARAMETER :: GREENSF_GENERAL_CONST = 1
   INTEGER,    PARAMETER :: GREENSF_HUBBARD_CONST = 2

   CONTAINS

   SUBROUTINE openGreensFFile(fileID, input, gfinp, atoms, inFilename)

      USE m_types
      USE m_cdn_io

      TYPE(t_input),                INTENT(IN)  :: input
      TYPE(t_gfinp),                INTENT(IN)  :: gfinp
      TYPE(t_atoms),                INTENT(IN)  :: atoms
      CHARACTER(len=*), OPTIONAL,   INTENT(IN)  :: inFilename
      INTEGER(HID_T),               INTENT(OUT) :: fileID

      LOGICAL           :: l_exist
      CHARACTER(LEN=30) :: filename
      INTEGER(HID_T)    :: metaGroupID
      INTEGER(HID_T)    :: generalGroupID

      LOGICAL           :: l_error
      INTEGER           :: hdfError
      INTEGER           :: version
      REAL              :: eFermiPrev

      version = 1
      IF(PRESENT(inFilename)) THEN
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

   END SUBROUTINE openGreensFFile

   SUBROUTINE closeGreensFFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE closeGreensFFile

   SUBROUTINE writeGreensFData(fileID, input, gfinp, atoms, archiveType, greensf, mmpmat, selfen)

      USE m_types
      USE m_constants
      USE m_juDFT

      INTEGER(HID_T),      INTENT(IN)  :: fileID
      TYPE(t_input),       INTENT(IN)  :: input
      TYPE(t_gfinp),       INTENT(IN)  :: gfinp
      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_greensf),     INTENT(IN)  :: greensf(:)
      INTEGER,             INTENT(IN)  :: archiveType
      COMPLEX,             INTENT(IN)  :: mmpmat(-lmaxU_Const:,-lmaxU_Const:,:,:)
      COMPLEX, OPTIONAL,   INTENT(IN)  :: selfen(:,:,:,:,:) !Only in IO mode for Hubbard 1

      INTEGER(HID_T)    :: elementsGroupID
      INTEGER(HID_T)    :: currentelementGroupID
      INTEGER(HID_T)    :: mmpmatSpaceID, mmpmatSetID
      INTEGER(HID_T)    :: sphavgDataSpaceID, sphavgDataSetID
      INTEGER(HID_T)    :: uuDataSpaceID, uuDataSetID
      INTEGER(HID_T)    :: udDataSpaceID, udDataSetID
      INTEGER(HID_T)    :: duDataSpaceID, duDataSetID
      INTEGER(HID_T)    :: ddDataSpaceID, ddDataSetID
      INTEGER(HID_T)    :: selfenDataSpaceID, selfenDataSetID
      INTEGER(HID_T)    :: energyPointsSpaceID, energyPointsSetID
      INTEGER(HID_T)    :: energyWeightsSpaceID, energyWeightsSetID


      CHARACTER(len=30) :: elementName, groupName
      INTEGER           :: hdfError
      INTEGER           :: dimsInt(7)
      INTEGER           :: i_gf,ispin,m
      INTEGER           :: i_elem,n_elem
      INTEGER(HSIZE_T)  :: dims(7)
      REAL              :: trc(input%jspins)


      SELECT CASE(archiveType)

      CASE(GREENSF_GENERAL_CONST)
         groupName = '/GreensFunctionElements'
         n_elem = gfinp%n
      CASE(GREENSF_HUBBARD_CONST)
         groupName = 'Hubbard1Elements'
         n_elem = atoms%n_hia
      CASE DEFAULT
         CALL juDFT_error("Unknown GF archiveType", calledby="writeGreensFData")
      END SELECT

      CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), elementsGroupID, hdfError)
      CALL io_write_attint0(elementsGroupID,'NumElements',n_elem)
      CALL io_write_attint0(elementsGroupID,'maxl',lmaxU_Const)

      DO i_elem = 1, n_elem

         WRITE(elementName,200) i_elem
200      FORMAT('element-',i0)

         SELECT CASE(archiveType)

         CASE(GREENSF_GENERAL_CONST)
            i_gf = i_elem
         CASE(GREENSF_HUBBARD_CONST)
            i_gf = gfinp%hiaElem(i_elem)
         CASE DEFAULT
            CALL juDFT_error("Unknown GF archiveType", calledby="writeGreensFData")
         END SELECT

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

         CALL io_write_attint0(currentelementGroupID,'nz',greensf(i_gf)%contour%nz)
         !CALL io_write_attint0(currentelementGroupID,'shape',gfinp%mode) !Replace with string description

         dims(:2)=[2,greensf(i_gf)%contour%nz]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),energyPointsSpaceID,hdfError)
         CALL h5dcreate_f(currentelementGroupID, "ContourPoints", H5T_NATIVE_DOUBLE, energyPointsSpaceID, energyPointsSetID, hdfError)
         CALL h5sclose_f(energyPointsSpaceID,hdfError)
         CALL io_write_complex1(energyPointsSetID,[-1,1],dimsInt(:2),greensf(i_gf)%contour%e)
         CALL h5dclose_f(energyPointsSetID, hdfError)
         dims(:2)=[2,greensf(i_gf)%contour%nz]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),energyWeightsSpaceID,hdfError)
         CALL h5dcreate_f(currentelementGroupID, "IntegrationWeights", H5T_NATIVE_DOUBLE, energyWeightsSpaceID, energyWeightsSetID, hdfError)
         CALL h5sclose_f(energyWeightsSpaceID,hdfError)
         CALL io_write_complex1(energyWeightsSetID,[-1,1],dimsInt(:2),greensf(i_gf)%contour%de)
         CALL h5dclose_f(energyWeightsSetID, hdfError)


         !Spherically averaged greensfData
         IF(gfinp%l_sphavg) THEN

            dims(:6)=[2,greensf(i_gf)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),sphavgDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "SphAvg", H5T_NATIVE_DOUBLE, sphavgDataSpaceID, sphavgDataSetID, hdfError)
            CALL h5sclose_f(sphavgDataSpaceID,hdfError)
            CALL io_write_complex5(sphavgDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf(i_gf)%gmmpmat(:,:,:,:,:))
            CALL h5dclose_f(sphavgDataSetID, hdfError)

         ELSE IF(.NOT.gfinp%l_sphavg.AND.archiveType.NE.GREENSF_HUBBARD_CONST) THEN

            !uu
            dims(:6)=[2,greensf(i_gf)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),uuDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "UU", H5T_NATIVE_DOUBLE, uuDataSpaceID, uuDataSetID, hdfError)
            CALL h5sclose_f(uuDataSpaceID,hdfError)
            CALL io_write_complex5(uuDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf(i_gf)%uu(:,:,:,:,:))
            CALL h5dclose_f(uuDataSetID, hdfError)

            !ud
            dims(:6)=[2,greensf(i_gf)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),udDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "UD", H5T_NATIVE_DOUBLE, udDataSpaceID, udDataSetID, hdfError)
            CALL h5sclose_f(udDataSpaceID,hdfError)
            CALL io_write_complex5(udDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf(i_gf)%ud(:,:,:,:,:))
            CALL h5dclose_f(udDataSetID, hdfError)

            !du
            dims(:6)=[2,greensf(i_gf)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),duDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "DU", H5T_NATIVE_DOUBLE, duDataSpaceID, duDataSetID, hdfError)
            CALL h5sclose_f(duDataSpaceID,hdfError)
            CALL io_write_complex5(duDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf(i_gf)%du(:,:,:,:,:))
            CALL h5dclose_f(duDataSetID, hdfError)

            !dd
            dims(:6)=[2,greensf(i_gf)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,MERGE(3,input%jspins,gfinp%l_mperp),2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),ddDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "DD", H5T_NATIVE_DOUBLE, ddDataSpaceID, ddDataSetID, hdfError)
            CALL h5sclose_f(ddDataSpaceID,hdfError)
            CALL io_write_complex5(ddDataSetID,[-1,1,1,1,1,1],dimsInt(:6),greensf(i_gf)%dd(:,:,:,:,:))
            CALL h5dclose_f(ddDataSetID, hdfError)

            !TODO write radial functions
         ENDIF

         IF(archiveType.EQ.GREENSF_HUBBARD_CONST.AND.PRESENT(selfen)) THEN
            dims(:5)=[2,2*(2*lmaxU_Const+1),2*(2*lmaxU_Const+1),greensf(i_gf)%contour%nz,2]
            dimsInt=dims
            CALL h5screate_simple_f(5,dims(:5),selfenDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "selfen", H5T_NATIVE_DOUBLE, selfenDataSpaceID, selfenDataSetID, hdfError)
            CALL h5sclose_f(selfenDataSpaceID,hdfError)
            CALL io_write_complex4(selfenDataSetID,[-1,1,1,1,1],dimsInt(:5),selfen(:,:,:greensf(i_gf)%contour%nz,:,i_elem))
            CALL h5dclose_f(selfenDataSetID, hdfError)
         ENDIF

         CALL h5gclose_f(currentelementGroupID, hdfError)
      ENDDO

      CALL h5gclose_f(elementsGroupID, hdfError)

   END SUBROUTINE writeGreensFData

#endif

END MODULE m_greensf_io

