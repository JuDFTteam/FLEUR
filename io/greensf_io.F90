MODULE m_greensf_io

#ifdef CPP_HDF

   USE hdf5
   USE m_hdf_tools
   USE m_types
   USE m_types_selfen
   USE m_constants
   USE m_juDFT

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

   SUBROUTINE openGreensFFile(fileID, input, gfinp, atoms, kpts, inFilename)

      USE m_types
      USE m_cdn_io

      TYPE(t_input),                INTENT(IN)  :: input
      TYPE(t_gfinp),                INTENT(IN)  :: gfinp
      TYPE(t_atoms),                INTENT(IN)  :: atoms
      TYPE(t_kpts),                 INTENT(IN)  :: kpts
      CHARACTER(len=*), OPTIONAL,   INTENT(IN)  :: inFilename
      INTEGER(HID_T),               INTENT(OUT) :: fileID

      LOGICAL           :: l_exist
      CHARACTER(LEN=30) :: filename
      INTEGER(HID_T)    :: metaGroupID
      INTEGER(HID_T)    :: generalGroupID, kptsGroupID
      INTEGER(HID_T)    :: kptCoordSpaceID, kptCoordSetID
      INTEGER(HID_T)    :: kptWeightSpaceID, kptWeightSetID
      INTEGER(HID_T)    :: kptsSPLabelsSpaceID, kptsSPLabelsSetID
      INTEGER(HID_T)    :: kptsSPIndicesSpaceID, kptsSPIndicesSetID

      INTEGER(HID_T)    :: stringTypeID
      INTEGER(SIZE_T)   :: stringLength

      LOGICAL           :: l_error
      INTEGER           :: hdfError
      INTEGER           :: version
      REAL              :: eFermiPrev
      INTEGER           :: dimsInt(7)
      INTEGER(HSIZE_T)  :: dims(7)

      version = 5
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
      CALL io_write_attlog0(generalGroupID,'mperp',gfinp%l_mperp)

      CALL h5gcreate_f(generalGroupID, 'kpts', kptsGroupID, hdfError)

      CALL io_write_attint0(kptsGroupID,'nkpt',kpts%nkpt)
      CALL io_write_attchar0(kptsGroupID,'kind',TRIM(ADJUSTL(kptsKindString_consts(kpts%kptsKind))))
      CALL io_write_attint0(kptsGroupID,'nSpecialPoints',kpts%numSpecialPoints)

      dims(:2)=(/3,kpts%nkpt/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),kptCoordSpaceID,hdfError)
      CALL h5dcreate_f(kptsGroupID, "coordinates", H5T_NATIVE_DOUBLE, kptCoordSpaceID, kptCoordSetID, hdfError)
      CALL h5sclose_f(kptCoordSpaceID,hdfError)
      CALL io_write_real2(kptCoordSetID,(/1,1/),dimsInt(:2),kpts%bk)
      CALL h5dclose_f(kptCoordSetID, hdfError)

      dims(:1)=(/kpts%nkpt/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),kptWeightSpaceID,hdfError)
      CALL h5dcreate_f(kptsGroupID, "weights", H5T_NATIVE_DOUBLE, kptWeightSpaceID, kptWeightSetID, hdfError)
      CALL h5sclose_f(kptWeightSpaceID,hdfError)
      CALL io_write_real1(kptWeightSetID,(/1/),dimsInt(:1),kpts%wtkpt)
      CALL h5dclose_f(kptWeightSetID, hdfError)

      IF (ALLOCATED(kpts%specialPointIndices)) THEN
         stringLength = LEN(kpts%specialPointNames(:))
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeID, hdfError)
         CALL h5tset_size_f(stringTypeID, stringLength, hdfError)
         CALL h5tset_strpad_f(stringTypeID, H5T_STR_SPACEPAD_F, hdfError)
         CALL h5tset_cset_f(stringTypeID, H5T_CSET_ASCII_F, hdfError)
         dims(:1)=(/kpts%numSpecialPoints/)
         dimsInt=dims
         CALL h5screate_simple_f(1,dims(:1),kptsSPLabelsSpaceID,hdfError)
         CALL h5dcreate_f(kptsGroupID, "specialPointLabels", stringTypeID, kptsSPLabelsSpaceID, kptSPLabelsSetID, hdfError)
         CALL h5tclose_f(stringTypeID,hdfError)
         CALL h5sclose_f(kptsSPLabelsSpaceID,hdfError)
         CALL io_write_string1(kptsSPLabelsSetID,dimsInt(:1),LEN(kpts%specialPointNames(:)),kpts%specialPointNames)
         CALL h5dclose_f(kptsSPLabelsSetID, hdfError)

         dims(:1)=(/kpts%numSpecialPoints/)
         dimsInt=dims
         CALL h5screate_simple_f(1,dims(:1),kptsSPIndicesSpaceID,hdfError)
         CALL h5dcreate_f(kptsGroupID, "specialPointIndices", H5T_NATIVE_INTEGER, kptsSPIndicesSpaceID, kptsSPIndicesSetID, hdfError)
         CALL h5sclose_f(kptsSPIndicesSpaceID,hdfError)
         CALL io_write_integer1(kptsSPIndicesSetID,(/1/),dimsInt(:1),kpts%specialPointIndices)
         CALL h5dclose_f(kptsSPIndicesSetID, hdfError)
      END IF

      CALL h5gclose_f(kptsGroupID, hdfError)
      CALL h5gclose_f(generalGroupID, hdfError)

   END SUBROUTINE openGreensFFile

   SUBROUTINE closeGreensFFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE closeGreensFFile

   SUBROUTINE writeGreensFData(fileID, input, gfinp, atoms, cell, archiveType, greensf, mmpmat, selfen,&
                               u, udot, ulo)

      INTEGER(HID_T),   INTENT(IN) :: fileID
      TYPE(t_input),    INTENT(IN) :: input
      TYPE(t_gfinp),    INTENT(IN) :: gfinp
      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_cell),     INTENT(IN) :: cell
      TYPE(t_greensf),  INTENT(IN) :: greensf(:)
      INTEGER,          INTENT(IN) :: archiveType
      COMPLEX,          INTENT(IN) :: mmpmat(-lmaxU_Const:,-lmaxU_Const:,:,:)

      TYPE(t_selfen),   OPTIONAL, INTENT(IN) :: selfen(:) !Only in IO mode for Hubbard 1
      REAL,             OPTIONAL, INTENT(IN) :: u(:,:,0:,:,:)      !Radial Functions for IO
      REAL,             OPTIONAL, INTENT(IN) :: udot(:,:,0:,:,:)
      REAL,             OPTIONAL, INTENT(IN) :: ulo(:,:,:,:,:)

      INTEGER(HID_T)    :: elementsGroupID,contoursGroupID,radialGroupID
      INTEGER(HID_T)    :: currentelementGroupID,currentcontourGroupID
      INTEGER(HID_T)    :: mmpmatSpaceID, mmpmatSetID
      INTEGER(HID_T)    :: selfenDataSpaceID, selfenDataSetID
      INTEGER(HID_T)    :: energyPointsSpaceID, energyPointsSetID
      INTEGER(HID_T)    :: energyWeightsSpaceID, energyWeightsSetID
      INTEGER(HID_T)    :: uDataSpaceID,uDataSetID
      INTEGER(HID_T)    :: udotDataSpaceID,udotDataSetID
      INTEGER(HID_T)    :: uloDataSpaceID,uloDataSetID
      INTEGER(HID_T)    :: nLODataSpaceID,nLODataSetID
      INTEGER(HID_T)    :: lLODataSpaceID,lLODataSetID
      INTEGER(HID_T)    :: DataSpaceID, DataSetID

      CHARACTER(len=30) :: elementName, groupName, shapeStr
      INTEGER           :: hdfError
      INTEGER           :: dimsInt(7)
      INTEGER           :: ispin,m,jspinsOut,iContour
      INTEGER           :: i_elem,i,iContourOut,nLO
      INTEGER           :: contour_mapping(gfinp%numberContours)
      INTEGER(HSIZE_T)  :: dims(7)
      COMPLEX           :: trc(3)
      LOGICAL           :: l_anyradial
      TYPE(t_greensf)   :: gfOut


      contour_mapping = -1
      jspinsOut = MERGE(3,input%jspins,gfinp%l_mperp)

      !Check dimensions of mmpmat and selfen
      IF(SIZE(mmpmat,3) /= SIZE(greensf)) CALL juDFT_error("Mismatch in sizes: mmpmat", calledby="writeGreensFData")
      IF(PRESENT(selfen)) THEN
         IF(SIZE(selfen) /= SIZE(greensf)) CALL juDFT_error("Mismatch in sizes: selfen", calledby="writeGreensFData")
         IF(archiveType /= GREENSF_HUBBARD_CONST) CALL juDFT_error("Wrong archiveType for selfen", calledby="writeGreensFData")
      ENDIF

      IF(PRESENT(u)) THEN
         IF(archiveType /= GREENSF_GENERAL_CONST) CALL juDFT_error("Wrong archiveType for u", calledby="writeGreensFData")
         IF(.NOT.PRESENT(udot)) CALL juDFT_error("udot not provided (u is present)", calledby="writeGreensFData")
      ENDIF

      !--> Start: Energy Contour Output
      CALL h5gcreate_f(fileID, '/EnergyContours', contoursGroupID, hdfError)

      iContourOut = 0
      DO iContour = 1, gfinp%numberContours
         !Find a greens function element which has this contour (if not skip)
         i_elem = -1
         DO i = 1, SIZE(greensf)
            IF(iContour == greensf(i)%elem%iContour) THEN
               i_elem = i
               EXIT
            ENDIF
         ENDDO

         IF(i_elem==-1) CYCLE

         iContourOut = iContourOut + 1
         WRITE(elementName,100) iContourOut
100      FORMAT('contour-',i0)
         contour_mapping(iContour) = iContourOut

         CALL h5gcreate_f(contoursGroupID, elementName, currentcontourGroupID, hdfError)

         CALL io_write_attint0(currentcontourGroupID,'nz',greensf(i_elem)%contour%nz)

         SELECT CASE (gfinp%contour(iContour)%shape)

         CASE(CONTOUR_RECTANGLE_CONST)
            shapeStr = 'Rectangle'
         CASE(CONTOUR_SEMICIRCLE_CONST)
            shapeStr = 'Semicircle'
         CASE(CONTOUR_DOS_CONST)
            shapeStr = 'DOS'
         CASE DEFAULT
         END SELECT

         CALL io_write_attchar0(currentcontourGroupID,'contourShape',TRIM(ADJUSTL(shapeStr)))
         CALL io_write_attchar0(currentcontourGroupID,'contourLabel',TRIM(ADJUSTL(gfinp%contour(iContour)%label)))

         dims(:2)=[2,greensf(i_elem)%contour%nz]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),energyPointsSpaceID,hdfError)
         CALL h5dcreate_f(currentcontourGroupID, "ContourPoints", H5T_NATIVE_DOUBLE, energyPointsSpaceID, energyPointsSetID, hdfError)
         CALL h5sclose_f(energyPointsSpaceID,hdfError)
         CALL io_write_complex1(energyPointsSetID,[-1,1],dimsInt(:2),greensf(i_elem)%contour%e)
         CALL h5dclose_f(energyPointsSetID, hdfError)
         dims(:2)=[2,greensf(i_elem)%contour%nz]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),energyWeightsSpaceID,hdfError)
         CALL h5dcreate_f(currentcontourGroupID, "IntegrationWeights", H5T_NATIVE_DOUBLE, energyWeightsSpaceID, energyWeightsSetID, hdfError)
         CALL h5sclose_f(energyWeightsSpaceID,hdfError)
         CALL io_write_complex1(energyWeightsSetID,[-1,1],dimsInt(:2),greensf(i_elem)%contour%de)
         CALL h5dclose_f(energyWeightsSetID, hdfError)

         CALL h5gclose_f(currentcontourGroupID, hdfError)
      ENDDO
      CALL io_write_attint0(contoursGroupID,'NumContours',iContourOut)
      CALL h5gclose_f(contoursGroupID, hdfError)
      !--> End: Energy Contour Output


      !--> Start: GF data output
      SELECT CASE(archiveType)

      CASE(GREENSF_GENERAL_CONST)
         groupName = '/GreensFunctionElements'
      CASE(GREENSF_HUBBARD_CONST)
         groupName = '/Hubbard1Elements'
      CASE DEFAULT
         CALL juDFT_error("Unknown GF archiveType", calledby="writeGreensFData")
      END SELECT

      CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), elementsGroupID, hdfError)
      CALL io_write_attint0(elementsGroupID,'NumElements',SIZE(greensf))
      CALL io_write_attint0(elementsGroupID,'maxl',lmaxU_Const)


      l_anyradial = .FALSE.
      DO i_elem = 1, SIZE(greensf)

         WRITE(elementName,200) i_elem
200      FORMAT('element-',i0)

         IF(.NOT.greensf(i_elem)%l_sphavg.AND.gfinp%l_outputSphavg) THEN
            gfOut = greensf(i_elem)%integrateoverMT(atoms,input,gfinp,u,udot,ulo,l_fullRadial=gfinp%l_intFullRadial)
         ELSE
            gfOut = greensf(i_elem)
         ENDIF

         CALL h5gcreate_f(elementsGroupID, elementName, currentelementGroupID, hdfError)
         nLO = greensf(i_elem)%elem%countLOs(atoms)
         IF(nLO>0 .AND..NOT.gfOut%l_sphavg.AND.PRESENT(u).AND..NOT.PRESENT(ulo)) THEN
            CALL juDFT_error("LO Radial Functions needed, but not present", calledby="writeGreensFData")
         ENDIF
         IF(.NOT.gfOut%l_sphavg) l_anyradial = .TRUE.

         !Trace of occupation matrix
         trc=0.0
         DO ispin = 1, jspinsOut
            DO m = -gfOut%elem%l, gfOut%elem%l
               trc(ispin) = trc(ispin) + mmpmat(m,m,i_elem,ispin)
            ENDDO
         ENDDO
         CALL io_write_attreal0(currentelementGroupID,"SpinUpTrace",REAL(trc(1)))
         IF(input%jspins.EQ.2) THEN
            CALL io_write_attreal0(currentelementGroupID,"SpinDownTrace",REAL(trc(2)))
         ENDIF
         IF(gfinp%l_mperp) THEN
            CALL io_write_attreal0(currentelementGroupID,"OffDTrace-x",REAL(trc(3)))
            CALL io_write_attreal0(currentelementGroupID,"OffDTrace-y",AIMAG(trc(3)))
         ENDIF

         CALL writeGreensFElement(currentelementGroupID, gfOut, atoms, cell, jspinsOut, contour_mapping)

         !Occupation matrix
         dims(:4)=[2,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut]
         dimsInt=dims
         CALL h5screate_simple_f(4,dims(:4),mmpmatSpaceID,hdfError)
         CALL h5dcreate_f(currentelementGroupID, "mmpmat", H5T_NATIVE_DOUBLE, mmpmatSpaceID, mmpmatSetID, hdfError)
         CALL h5sclose_f(mmpmatSpaceID,hdfError)
         CALL io_write_complex3(mmpmatSetID,[-1,1,1,1],dimsInt(:4),mmpmat(:,:,i_elem,:jspinsOut))
         CALL h5dclose_f(mmpmatSetID, hdfError)

         IF(archiveType.EQ.GREENSF_HUBBARD_CONST.AND.PRESENT(selfen)) THEN
            dims(:5)=[2,2*(2*lmaxU_Const+1),2*(2*selfen(i_elem)%l+1),greensf(i_elem)%contour%nz,2]
            dimsInt=dims
            CALL h5screate_simple_f(5,dims(:5),selfenDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "selfen", H5T_NATIVE_DOUBLE, selfenDataSpaceID, selfenDataSetID, hdfError)
            CALL h5sclose_f(selfenDataSpaceID,hdfError)
            CALL io_write_complex4(selfenDataSetID,[-1,1,1,1,1],dimsInt(:5),selfen(i_elem)%data)
            CALL h5dclose_f(selfenDataSetID, hdfError)
         ENDIF

         CALL h5gclose_f(currentelementGroupID, hdfError)
      ENDDO

      CALL h5gclose_f(elementsGroupID, hdfError)
      !--> End: GF data output

      !--> Start: Radial Function output
      IF(PRESENT(u).AND.l_anyradial) THEN
         CALL h5gcreate_f(fileID, 'RadialFunctions', radialGroupID, hdfError)

         dims(:2)=[atoms%jmtd,atoms%ntype]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),DataSpaceID,hdfError)
         CALL h5dcreate_f(radialGroupID, "rmsh", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real2(DataSetID,[1,1],dimsInt(:2),atoms%rmsh)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:1)=[atoms%ntype]
         dimsInt=dims
         CALL h5screate_simple_f(1,dims(:1),DataSpaceID,hdfError)
         CALL h5dcreate_f(radialGroupID, "jri", H5T_NATIVE_INTEGER, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_integer1(DataSetID,[1],dimsInt(:1),atoms%jri)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:5)=[atoms%jmtd,2,lmaxU_Const+1,input%jspins,atoms%ntype]
         dimsInt=dims
         CALL h5screate_simple_f(5,dims(:5),uDataSpaceID,hdfError)
         CALL h5dcreate_f(radialGroupID, "u", H5T_NATIVE_DOUBLE, uDataSpaceID, uDataSetID, hdfError)
         CALL h5sclose_f(uDataSpaceID,hdfError)
         CALL io_write_real5(uDataSetID,[1,1,1,1,1],dimsInt(:5),u(:,:,0:lmaxU_Const,:,:))
         CALL h5dclose_f(uDataSetID, hdfError)

         dims(:5)=[atoms%jmtd,2,lmaxU_Const+1,input%jspins,atoms%ntype]
         dimsInt=dims
         CALL h5screate_simple_f(5,dims(:5),udotDataSpaceID,hdfError)
         CALL h5dcreate_f(radialGroupID, "udot", H5T_NATIVE_DOUBLE, udotDataSpaceID, udotDataSetID, hdfError)
         CALL h5sclose_f(udotDataSpaceID,hdfError)
         CALL io_write_real5(udotDataSetID,[1,1,1,1,1],dimsInt(:5),udot(:,:,0:lmaxU_Const,:,:))
         CALL h5dclose_f(udotDataSetID, hdfError)

         IF(PRESENT(ulo)) THEN
            !Mapping array
            dims(:1) = [atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(1,dims(:1),nLODataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "nlo", H5T_NATIVE_INTEGER, nLODataSpaceID, nLODataSetID, hdfError)
            CALL h5sclose_f(nLODataSpaceID,hdfError)
            CALL io_write_integer1(nLODataSetID,[1],dimsInt(:1),atoms%nlo)
            CALL h5dclose_f(nLODataSetID, hdfError)

            dims(:2) = [atoms%nlod,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),lLODataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "llo", H5T_NATIVE_INTEGER, lLODataSpaceID, lLODataSetID, hdfError)
            CALL h5sclose_f(lLODataSpaceID,hdfError)
            CALL io_write_integer2(lLODataSetID,[1,1],dimsInt(:2),atoms%llo)
            CALL h5dclose_f(lLODataSetID, hdfError)

            dims(:5)=[atoms%jmtd,2,atoms%nlod,input%jspins,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(5,dims(:5),uloDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "ulo", H5T_NATIVE_DOUBLE, uloDataSpaceID, uloDataSetID, hdfError)
            CALL h5sclose_f(uloDataSpaceID,hdfError)
            CALL io_write_real5(uloDataSetID,[1,1,1,1,1],dimsInt(:5),ulo)
            CALL h5dclose_f(uloDataSetID, hdfError)
         ENDIF

         CALL h5gclose_f(radialGroupID, hdfError)
      ENDIF
      !--> End: Radial Function output


   END SUBROUTINE writeGreensFData

   SUBROUTINE writeGreensFElement(groupID, g, atoms, cell, jspins, contour_mapping)

      INTEGER(HID_T),   INTENT(IN)  :: groupID
      TYPE(t_greensf),  INTENT(IN)  :: g
      TYPE(t_atoms),    INTENT(IN)  :: atoms
      TYPE(t_cell),     INTENT(IN)  :: cell
      INTEGER,          INTENT(IN)  :: jspins
      INTEGER,          INTENT(IN)  :: contour_mapping(:)


      CHARACTER(len=30) :: groupName, datasetName
      INTEGER :: dimsInt(7)
      INTEGER(HSIZE_T)  :: dims(7)
      INTEGER(HID_T) :: DataSpaceID, DataSetID
      INTEGER(HID_T)    :: loGroupID,currentloGroupID, scalarGroupID
      INTEGER :: hdfError,ikpt
      INTEGER :: nLO, iLO, iLOp

      CALL io_write_attint0(groupID,"l",g%elem%l)
      CALL io_write_attint0(groupID,"lp",g%elem%lp)
      CALL io_write_attint0(groupID,"atomType",g%elem%atomType)
      CALL io_write_attint0(groupID,"atomTypep",g%elem%atomTypep)
      IF(g%elem%atom/=0) THEN
         CALL io_write_attchar0(groupID,"atom",TRIM(ADJUSTL(atoms%label(g%elem%atom))))
         CALL io_write_attchar0(groupID,"atomp",TRIM(ADJUSTL(atoms%label(g%elem%atomp))))
      ELSE
         CALL io_write_attchar0(groupID,"atom",'0')
         CALL io_write_attchar0(groupID,"atomp",'0')
      ENDIF
      CALL io_write_attint0(groupID,'iContour',contour_mapping(g%elem%iContour))
      CALL io_write_attlog0(groupID,'l_onsite',.NOT.g%elem%isOffDiag())
      CALL io_write_attlog0(groupID,'l_sphavg',g%l_sphavg)
      CALL io_write_attlog0(groupID,'l_kresolved',g%elem%l_kresolved)
      CALL io_write_attreal1(groupID,'atomDiff',matmul(cell%amat,g%elem%atomDiff))
      nLO = g%elem%countLOs(atoms)
      CALL io_write_attint0(groupID,'numLOs',nLO)

      IF(g%l_sphavg) THEN

         dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
         dimsInt=dims
         CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "sphavg", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%gmmpmat)
         CALL h5dclose_f(DataSetID, hdfError)

      ELSE IF(g%l_kresolved) THEN
         dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
         dimsInt=dims
         DO ikpt = 1, SIZE(g%gmmpmat_k)
            WRITE(datasetName,201) ikpt
201         FORMAT('kresolved-',i0)
            CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
            CALL h5dcreate_f(groupID, TRIM(ADJUSTL(datasetName)), H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
            CALL h5sclose_f(DataSpaceID,hdfError)
            CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%gmmpmat_k(:,:,:,:,:,ikpt))
            CALL h5dclose_f(DataSetID, hdfError)
         ENDDO
      ELSE

         !--> Start: Radial Coefficients
         dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
         dimsInt=dims
         CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "uu", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%uu)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
         dimsInt=dims
         CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "ud", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%ud)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
         dimsInt=dims
         CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "du", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%du)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
         dimsInt=dims
         CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "dd", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%dd)
         CALL h5dclose_f(DataSetID, hdfError)
         !--> End: Radial Coefficients

         !--> Start: LO Coefficients
         IF(nLO>0) THEN

            CALL h5gcreate_f(groupID, 'LOcontribution', loGroupID, hdfError)
            DO iLO = 1, nLO
               WRITE(groupName,300) iLO
300            FORMAT('lo-',i0)
               CALL h5gcreate_f(loGroupID, groupName, currentloGroupID, hdfError)

               dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
               dimsInt=dims
               CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
               CALL h5dcreate_f(currentloGroupID, "uulo", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
               CALL h5sclose_f(DataSpaceID,hdfError)
               CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%uulo(:,:,:,iLO,:,:))
               CALL h5dclose_f(DataSetID, hdfError)

               dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
               dimsInt=dims
               CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
               CALL h5dcreate_f(currentloGroupID, "ulou", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
               CALL h5sclose_f(DataSpaceID,hdfError)
               CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%ulou(:,:,:,iLO,:,:))
               CALL h5dclose_f(DataSetID, hdfError)

               dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
               dimsInt=dims
               CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
               CALL h5dcreate_f(currentloGroupID, "dulo", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
               CALL h5sclose_f(DataSpaceID,hdfError)
               CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%dulo(:,:,:,iLO,:,:))
               CALL h5dclose_f(DataSetID, hdfError)

               dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
               dimsInt=dims
               CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
               CALL h5dcreate_f(currentloGroupID, "ulod", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
               CALL h5sclose_f(DataSpaceID,hdfError)
               CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%ulod(:,:,:,iLO,:,:))
               CALL h5dclose_f(DataSetID, hdfError)


               DO iLop = 1, nLO
                  WRITE(datasetName,400) iLop
400               FORMAT('uloulop-',i0)

                  dims(:6)=[2,g%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspins,2]
                  dimsInt=dims
                  CALL h5screate_simple_f(6,dims(:6),DataSpaceID,hdfError)
                  CALL h5dcreate_f(currentloGroupID, TRIM(ADJUSTL(datasetName)), H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
                  CALL h5sclose_f(DataSpaceID,hdfError)
                  CALL io_write_complex5(DataSetID,[-1,1,1,1,1,1],dimsInt(:6),g%uloulop(:,:,:,iLO,iLOp,:,:))
                  CALL h5dclose_f(DataSetID, hdfError)
               ENDDO
               CALL h5gclose_f(currentloGroupID, hdfError)

            ENDDO
            CALL h5gclose_f(loGroupID, hdfError)
            !--> End: LO Coefficients
         ENDIF

         !--> Start: Scalar Products
         CALL h5gcreate_f(groupID, 'scalarProducts', scalarGroupID, hdfError)

         dims(:2)=[2,2]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "uun", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real2(DataSetID,[1,1],dimsInt(:2),g%scalarProducts%uun)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:2)=[2,2]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "dun", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real2(DataSetID,[1,1],dimsInt(:2),g%scalarProducts%dun)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:2)=[2,2]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "udn", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real2(DataSetID,[1,1],dimsInt(:2),g%scalarProducts%udn)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:2)=[2,2]
         dimsInt=dims
         CALL h5screate_simple_f(2,dims(:2),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "ddn", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real2(DataSetID,[1,1],dimsInt(:2),g%scalarProducts%ddn)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:3)=[atoms%nlod,2,2]
         dimsInt=dims
         CALL h5screate_simple_f(3,dims(:3),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "uulon", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real3(DataSetID,[1,1,1],dimsInt(:3),g%scalarProducts%uulon)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:3)=[atoms%nlod,2,2]
         dimsInt=dims
         CALL h5screate_simple_f(3,dims(:3),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "uloun", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real3(DataSetID,[1,1,1],dimsInt(:3),g%scalarProducts%uloun)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:3)=[atoms%nlod,2,2]
         dimsInt=dims
         CALL h5screate_simple_f(3,dims(:3),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "dulon", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real3(DataSetID,[1,1,1],dimsInt(:3),g%scalarProducts%dulon)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:3)=[atoms%nlod,2,2]
         dimsInt=dims
         CALL h5screate_simple_f(3,dims(:3),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "ulodn", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real3(DataSetID,[1,1,1],dimsInt(:3),g%scalarProducts%ulodn)
         CALL h5dclose_f(DataSetID, hdfError)

         dims(:4)=[atoms%nlod,atoms%nlod,2,2]
         dimsInt=dims
         CALL h5screate_simple_f(4,dims(:4),DataSpaceID,hdfError)
         CALL h5dcreate_f(scalarGroupID, "uloulopn", H5T_NATIVE_DOUBLE, DataSpaceID, DataSetID, hdfError)
         CALL h5sclose_f(DataSpaceID,hdfError)
         CALL io_write_real4(DataSetID,[1,1,1,1],dimsInt(:4),g%scalarProducts%uloulopn)
         CALL h5dclose_f(DataSetID, hdfError)

         CALL h5gclose_f(scalarGroupID, hdfError)
         !--> End: Scalar Products

      ENDIF

   END SUBROUTINE writeGreensFElement

#endif

END MODULE m_greensf_io

