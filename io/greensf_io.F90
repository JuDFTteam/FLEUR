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
      CALL io_write_attlog0(generalGroupID,'mperp',gfinp%l_mperp)
      CALL h5gclose_f(generalGroupID, hdfError)

   END SUBROUTINE openGreensFFile

   SUBROUTINE closeGreensFFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE closeGreensFFile

   SUBROUTINE writeGreensFData(fileID, input, gfinp, atoms, cell, archiveType, greensf, mmpmat, selfen,&
                               u, udot, ulo, usdus, denCoeffsOffDiag, scalarGF)

      USE m_types
      USE m_types_selfen
      USE m_types_scalarGF
      USE m_constants
      USE m_juDFT

      INTEGER(HID_T),   INTENT(IN) :: fileID
      TYPE(t_input),    INTENT(IN) :: input
      TYPE(t_gfinp),    INTENT(IN) :: gfinp
      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_cell),     INTENT(IN) :: cell
      TYPE(t_greensf),  INTENT(IN) :: greensf(:)
      INTEGER,          INTENT(IN) :: archiveType
      COMPLEX,          INTENT(IN) :: mmpmat(-lmaxU_Const:,-lmaxU_Const:,:,:)

      TYPE(t_selfen),            OPTIONAL, INTENT(IN) :: selfen(:) !Only in IO mode for Hubbard 1
      REAL,                      OPTIONAL, INTENT(IN) :: u(:,:,0:,:,:)      !Radial Functions for IO
      REAL,                      OPTIONAL, INTENT(IN) :: udot(:,:,0:,:,:)
      REAL,                      OPTIONAL, INTENT(IN) :: ulo(:,:,:,:,:)
      TYPE(t_usdus),             OPTIONAL, INTENT(IN) :: usdus
      TYPE(t_denCoeffsOffDiag),  OPTIONAL, INTENT(IN) :: denCoeffsOffDiag
      TYPE(t_scalarGF),          OPTIONAL, INTENT(IN) :: scalarGF(:)

      INTEGER(HID_T)    :: elementsGroupID,contoursGroupID,radialGroupID
      INTEGER(HID_T)    :: currentelementGroupID,currentcontourGroupID
      INTEGER(HID_T)    :: mmpmatSpaceID, mmpmatSetID
      INTEGER(HID_T)    :: sphavgDataSpaceID, sphavgDataSetID
      INTEGER(HID_T)    :: uuDataSpaceID, uuDataSetID
      INTEGER(HID_T)    :: udDataSpaceID, udDataSetID
      INTEGER(HID_T)    :: duDataSpaceID, duDataSetID
      INTEGER(HID_T)    :: ddDataSpaceID, ddDataSetID
      INTEGER(HID_T)    :: selfenDataSpaceID, selfenDataSetID
      INTEGER(HID_T)    :: energyPointsSpaceID, energyPointsSetID
      INTEGER(HID_T)    :: energyWeightsSpaceID, energyWeightsSetID
      INTEGER(HID_T)    :: uDataSpaceID,uDataSetID
      INTEGER(HID_T)    :: udotDataSpaceID,udotDataSetID
      INTEGER(HID_T)    :: uloDataSpaceID,uloDataSetID
      INTEGER(HID_T)    :: loGroupID,currentloGroupID
      INTEGER(HID_T)    :: uuloDataSpaceID,uuloDataSetID
      INTEGER(HID_T)    :: ulouDataSpaceID,ulouDataSetID
      INTEGER(HID_T)    :: duloDataSpaceID,duloDataSetID
      INTEGER(HID_T)    :: ulodDataSpaceID,ulodDataSetID
      INTEGER(HID_T)    :: uloulopDataSpaceID,uloulopDataSetID
      INTEGER(HID_T)    :: nLODataSpaceID,nLODataSetID
      INTEGER(HID_T)    :: lLODataSpaceID,lLODataSetID
      INTEGER(HID_T)    :: ddnDataSpaceID,ddnDataSetID
      INTEGER(HID_T)    :: uulonDataSpaceID,uulonDataSetID
      INTEGER(HID_T)    :: dulonDataSpaceID,dulonDataSetID
      INTEGER(HID_T)    :: uloulopnDataSpaceID,uloulopnDataSetID
      INTEGER(HID_T)    :: uu21nDataSpaceID,uu21nDataSetID
      INTEGER(HID_T)    :: du21nDataSpaceID,du21nDataSetID
      INTEGER(HID_T)    :: ud21nDataSpaceID,ud21nDataSetID
      INTEGER(HID_T)    :: dd21nDataSpaceID,dd21nDataSetID
      INTEGER(HID_T)    :: uulo21nDataSpaceID,uulo21nDataSetID
      INTEGER(HID_T)    :: ulou21nDataSpaceID,ulou21nDataSetID
      INTEGER(HID_T)    :: dulo21nDataSpaceID,dulo21nDataSetID
      INTEGER(HID_T)    :: ulod21nDataSpaceID,ulod21nDataSetID
      INTEGER(HID_T)    :: uloulop21nDataSpaceID,uloulop21nDataSetID


      CHARACTER(len=30) :: elementName, groupName, shapeStr
      INTEGER           :: hdfError
      INTEGER           :: dimsInt(7)
      INTEGER           :: ispin,m,l,lp,atomType,atomTypep,jspinsOut,iContour
      INTEGER           :: i_elem,i,iContourOut,nLO,iLo,iLOp
      INTEGER(HSIZE_T)  :: dims(7)
      REAL              :: trc(MERGE(3,input%jspins,gfinp%l_mperp)),atomDiff(3)
      LOGICAL           :: l_sphavg,l_onsite,l_anyradial
      TYPE(t_greensf)   :: gfOut


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

         CALL h5gcreate_f(contoursGroupID, elementName, currentcontourGroupID, hdfError)

         CALL io_write_attint0(currentcontourGroupID,'nz',greensf(i_elem)%contour%nz)
         CALL io_write_attint0(currentcontourGroupID,'iContour',iContour)

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

         !Get information about the element
         l  = greensf(i_elem)%elem%l
         lp = greensf(i_elem)%elem%lp
         atomType  = greensf(i_elem)%elem%atomType
         atomTypep = greensf(i_elem)%elem%atomTypep
         l_sphavg  = greensf(i_elem)%elem%l_sphavg
         atomDiff  = matmul(cell%amat,greensf(i_elem)%elem%atomDiff)
         iContour  = greensf(i_elem)%elem%iContour

         IF(.NOT.l_sphavg.AND.gfinp%l_outputSphavg) THEN
            gfOut = greensf(i_elem)%integrateoverMT(atoms,input,gfinp,u,udot,ulo,usdus=usdus,denCoeffsOffdiag=denCoeffsOffDiag,&
                                                    scalarGF=scalarGF(i_elem),l_fullRadial=gfinp%l_intFullRadial)
            l_sphavg = .TRUE.
         ELSE
            gfOut = greensf(i_elem)
         ENDIF

         l_onsite = l.EQ.lp.AND.atomType.EQ.atomTypep.AND.ALL(ABS(atomDiff).LT.1e-12)
         CALL h5gcreate_f(elementsGroupID, elementName, currentelementGroupID, hdfError)
         CALL io_write_attint0(currentelementGroupID,"l",l)
         CALL io_write_attint0(currentelementGroupID,"lp",lp)
         CALL io_write_attint0(currentelementGroupID,"atomType",atomType)
         CALL io_write_attint0(currentelementGroupID,"atomTypep",atomTypep)
         CALL io_write_attint0(currentelementGroupID,'iContour',iContour)
         CALL io_write_attlog0(currentelementGroupID,'l_onsite',l_onsite)
         CALL io_write_attlog0(currentelementGroupID,'l_sphavg',l_sphavg)
         CALL io_write_attreal1(currentelementGroupID,'atomDiff',atomDiff)
         nLO = greensf(i_elem)%elem%countLOs(atoms)
         IF(nLO>0 .AND..NOT.l_sphavg.AND.PRESENT(u).AND..NOT.PRESENT(ulo)) THEN
            CALL juDFT_error("LO Radial Functions needed, but not present", calledby="writeGreensFData")
         ENDIF
         CALL io_write_attint0(currentelementGroupID,'numLOs',nLO)

         IF(.NOT.l_sphavg) l_anyradial = .TRUE.

         !Trace of occupation matrix
         trc=0.0
         DO ispin = 1, jspinsOut
            DO m = -l, l
               trc(ispin) = trc(ispin) + REAL(mmpmat(m,m,i_elem,ispin))
            ENDDO
         ENDDO
         CALL io_write_attreal0(currentelementGroupID,"SpinUpTrace",trc(1))
         IF(input%jspins.EQ.2) THEN
            CALL io_write_attreal0(currentelementGroupID,"SpinDownTrace",trc(2))
         ENDIF
         IF(gfinp%l_mperp) THEN
            CALL io_write_attreal0(currentelementGroupID,"OffDTrace",trc(3))
         ENDIF

         !Occupation matrix
         dims(:4)=[2,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut]
         dimsInt=dims
         CALL h5screate_simple_f(4,dims(:4),mmpmatSpaceID,hdfError)
         CALL h5dcreate_f(currentelementGroupID, "mmpmat", H5T_NATIVE_DOUBLE, mmpmatSpaceID, mmpmatSetID, hdfError)
         CALL h5sclose_f(mmpmatSpaceID,hdfError)
         CALL io_write_complex3(mmpmatSetID,[-1,1,1,1],dimsInt(:4),mmpmat(:,:,i_elem,:jspinsOut))
         CALL h5dclose_f(mmpmatSetID, hdfError)

         !Spherically averaged greensfData
         IF(l_sphavg) THEN

            dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),sphavgDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "sphavg", H5T_NATIVE_DOUBLE, sphavgDataSpaceID, sphavgDataSetID, hdfError)
            CALL h5sclose_f(sphavgDataSpaceID,hdfError)
            CALL io_write_complex5(sphavgDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%gmmpmat)
            CALL h5dclose_f(sphavgDataSetID, hdfError)

         ELSE IF(.NOT.l_sphavg.AND.archiveType.NE.GREENSF_HUBBARD_CONST) THEN

            !uu
            dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),uuDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "uu", H5T_NATIVE_DOUBLE, uuDataSpaceID, uuDataSetID, hdfError)
            CALL h5sclose_f(uuDataSpaceID,hdfError)
            CALL io_write_complex5(uuDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%uu)
            CALL h5dclose_f(uuDataSetID, hdfError)

            !ud
            dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),udDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "ud", H5T_NATIVE_DOUBLE, udDataSpaceID, udDataSetID, hdfError)
            CALL h5sclose_f(udDataSpaceID,hdfError)
            CALL io_write_complex5(udDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%ud)
            CALL h5dclose_f(udDataSetID, hdfError)

            !du
            dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),duDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "du", H5T_NATIVE_DOUBLE, duDataSpaceID, duDataSetID, hdfError)
            CALL h5sclose_f(duDataSpaceID,hdfError)
            CALL io_write_complex5(duDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%du)
            CALL h5dclose_f(duDataSetID, hdfError)

            !dd
            dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
            dimsInt=dims
            CALL h5screate_simple_f(6,dims(:6),ddDataSpaceID,hdfError)
            CALL h5dcreate_f(currentelementGroupID, "dd", H5T_NATIVE_DOUBLE, ddDataSpaceID, ddDataSetID, hdfError)
            CALL h5sclose_f(ddDataSpaceID,hdfError)
            CALL io_write_complex5(ddDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%dd)
            CALL h5dclose_f(ddDataSetID, hdfError)

            !--> Start: LO contributions
            IF(nLO>0) THEN
               CALL h5gcreate_f(currentelementGroupID, 'LOcontribution', loGroupID, hdfError)
               DO iLO = 1, nLO
                  WRITE(elementName,300) iLO
300               FORMAT('lo-',i0)
                  CALL h5gcreate_f(loGroupID, elementName, currentloGroupID, hdfError)

                  !uulo
                  dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
                  dimsInt=dims
                  CALL h5screate_simple_f(6,dims(:6),uuloDataSpaceID,hdfError)
                  CALL h5dcreate_f(currentloGroupID, "uulo", H5T_NATIVE_DOUBLE, uuloDataSpaceID, uuloDataSetID, hdfError)
                  CALL h5sclose_f(uuloDataSpaceID,hdfError)
                  CALL io_write_complex5(uuloDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%uulo(:,:,:,iLO,:,:))
                  CALL h5dclose_f(uuloDataSetID, hdfError)

                  !ulou
                  dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
                  dimsInt=dims
                  CALL h5screate_simple_f(6,dims(:6),ulouDataSpaceID,hdfError)
                  CALL h5dcreate_f(currentloGroupID, "ulou", H5T_NATIVE_DOUBLE, ulouDataSpaceID, ulouDataSetID, hdfError)
                  CALL h5sclose_f(ulouDataSpaceID,hdfError)
                  CALL io_write_complex5(ulouDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%ulou(:,:,:,iLO,:,:))
                  CALL h5dclose_f(ulouDataSetID, hdfError)

                  !uulo
                  dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
                  dimsInt=dims
                  CALL h5screate_simple_f(6,dims(:6),duloDataSpaceID,hdfError)
                  CALL h5dcreate_f(currentloGroupID, "dulo", H5T_NATIVE_DOUBLE, duloDataSpaceID, duloDataSetID, hdfError)
                  CALL h5sclose_f(duloDataSpaceID,hdfError)
                  CALL io_write_complex5(duloDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%dulo(:,:,:,iLO,:,:))
                  CALL h5dclose_f(duloDataSetID, hdfError)

                  !ulou
                  dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
                  dimsInt=dims
                  CALL h5screate_simple_f(6,dims(:6),ulodDataSpaceID,hdfError)
                  CALL h5dcreate_f(currentloGroupID, "ulod", H5T_NATIVE_DOUBLE, ulodDataSpaceID, ulodDataSetID, hdfError)
                  CALL h5sclose_f(ulodDataSpaceID,hdfError)
                  CALL io_write_complex5(ulodDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%ulod(:,:,:,iLO,:,:))
                  CALL h5dclose_f(ulodDataSetID, hdfError)

                  DO iLop = 1, nLO
                     WRITE(elementName,400) iLop
400                  FORMAT('uloulop-',i0)
                     !uloulop
                     dims(:6)=[2,greensf(i_elem)%contour%nz,2*lmaxU_Const+1,2*lmaxU_Const+1,jspinsOut,2]
                     dimsInt=dims
                     CALL h5screate_simple_f(6,dims(:6),uloulopDataSpaceID,hdfError)
                     CALL h5dcreate_f(currentloGroupID, elementName, H5T_NATIVE_DOUBLE, uloulopDataSpaceID, uloulopDataSetID, hdfError)
                     CALL h5sclose_f(uloulopDataSpaceID,hdfError)
                     CALL io_write_complex5(uloulopDataSetID,[-1,1,1,1,1,1],dimsInt(:6),gfOut%uloulop(:,:,:,iLO,iLOp,:,:))
                     CALL h5dclose_f(uloulopDataSetID, hdfError)
                  ENDDO
                  CALL h5gclose_f(currentloGroupID, hdfError)
               ENDDO
               CALL h5gclose_f(loGroupID, hdfError)
            ENDIF
            !--> End: LO contributions

         ENDIF

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

         !--> Start: Scalar products
         IF(PRESENT(usdus)) THEN
            dims(:3) = [atoms%lmaxd+1,atoms%ntype,input%jspins]
            dimsInt=dims
            CALL h5screate_simple_f(3,dims(:3),ddnDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "ddn", H5T_NATIVE_DOUBLE, ddnDataSpaceID, ddnDataSetID, hdfError)
            CALL h5sclose_f(ddnDataSpaceID,hdfError)
            CALL io_write_real3(ddnDataSetID,[1,1,1],dimsInt(:3),usdus%ddn)
            CALL h5dclose_f(ddnDataSetID, hdfError)

            dims(:3) = [atoms%nlod,atoms%ntype,input%jspins]
            dimsInt=dims
            CALL h5screate_simple_f(3,dims(:3),uulonDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "uulon", H5T_NATIVE_DOUBLE, uulonDataSpaceID, uulonDataSetID, hdfError)
            CALL h5sclose_f(uulonDataSpaceID,hdfError)
            CALL io_write_real3(uulonDataSetID,[1,1,1],dimsInt(:3),usdus%uulon)
            CALL h5dclose_f(uulonDataSetID, hdfError)

            dims(:3) = [atoms%nlod,atoms%ntype,input%jspins]
            dimsInt=dims
            CALL h5screate_simple_f(3,dims(:3),dulonDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "dulon", H5T_NATIVE_DOUBLE, dulonDataSpaceID, dulonDataSetID, hdfError)
            CALL h5sclose_f(dulonDataSpaceID,hdfError)
            CALL io_write_real3(dulonDataSetID,[1,1,1],dimsInt(:3),usdus%dulon)
            CALL h5dclose_f(dulonDataSetID, hdfError)

            dims(:4) = [atoms%nlod,atoms%nlod,atoms%ntype,input%jspins]
            dimsInt=dims
            CALL h5screate_simple_f(4,dims(:4),uloulopnDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "uloulopn", H5T_NATIVE_DOUBLE, uloulopnDataSpaceID, uloulopnDataSetID, hdfError)
            CALL h5sclose_f(uloulopnDataSpaceID,hdfError)
            CALL io_write_real4(uloulopnDataSetID,[1,1,1,1],dimsInt(:4),usdus%uloulopn)
            CALL h5dclose_f(uloulopnDataSetID, hdfError)
         ENDIF

         IF(PRESENT(denCoeffsOffDiag).AND.gfinp%l_mperp) THEN
            dims(:2) = [atoms%lmaxd+1,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),uu21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "uu21n", H5T_NATIVE_DOUBLE, uu21nDataSpaceID, uu21nDataSetID, hdfError)
            CALL h5sclose_f(uu21nDataSpaceID,hdfError)
            CALL io_write_real2(uu21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%uu21n)
            CALL h5dclose_f(uu21nDataSetID, hdfError)

            dims(:2) = [atoms%lmaxd+1,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),ud21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "ud21n", H5T_NATIVE_DOUBLE, ud21nDataSpaceID, ud21nDataSetID, hdfError)
            CALL h5sclose_f(ud21nDataSpaceID,hdfError)
            CALL io_write_real2(ud21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%ud21n)
            CALL h5dclose_f(ud21nDataSetID, hdfError)

            dims(:2) = [atoms%lmaxd+1,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),du21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "du21n", H5T_NATIVE_DOUBLE, du21nDataSpaceID, du21nDataSetID, hdfError)
            CALL h5sclose_f(du21nDataSpaceID,hdfError)
            CALL io_write_real2(du21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%du21n)
            CALL h5dclose_f(du21nDataSetID, hdfError)

            dims(:2) = [atoms%lmaxd+1,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),dd21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "dd21n", H5T_NATIVE_DOUBLE, dd21nDataSpaceID, dd21nDataSetID, hdfError)
            CALL h5sclose_f(dd21nDataSpaceID,hdfError)
            CALL io_write_real2(dd21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%dd21n)
            CALL h5dclose_f(dd21nDataSetID, hdfError)

            dims(:2) = [atoms%nlod,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),uulo21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "uulo21n", H5T_NATIVE_DOUBLE, uulo21nDataSpaceID, uulo21nDataSetID, hdfError)
            CALL h5sclose_f(uulo21nDataSpaceID,hdfError)
            CALL io_write_real2(uulo21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%uulo21n)
            CALL h5dclose_f(uulo21nDataSetID, hdfError)

            dims(:2) = [atoms%nlod,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),ulou21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "ulou21n", H5T_NATIVE_DOUBLE, ulou21nDataSpaceID, ulou21nDataSetID, hdfError)
            CALL h5sclose_f(ulou21nDataSpaceID,hdfError)
            CALL io_write_real2(ulou21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%ulou21n)
            CALL h5dclose_f(ulou21nDataSetID, hdfError)

            dims(:2) = [atoms%nlod,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),dulo21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "dulo21n", H5T_NATIVE_DOUBLE, dulo21nDataSpaceID, dulo21nDataSetID, hdfError)
            CALL h5sclose_f(dulo21nDataSpaceID,hdfError)
            CALL io_write_real2(dulo21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%dulo21n)
            CALL h5dclose_f(dulo21nDataSetID, hdfError)

            dims(:2) = [atoms%nlod,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(2,dims(:2),ulod21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "ulod21n", H5T_NATIVE_DOUBLE, ulod21nDataSpaceID, ulod21nDataSetID, hdfError)
            CALL h5sclose_f(ulod21nDataSpaceID,hdfError)
            CALL io_write_real2(ulod21nDataSetID,[1,1],dimsInt(:2),denCoeffsOffDiag%ulod21n)
            CALL h5dclose_f(ulod21nDataSetID, hdfError)

            dims(:3) = [atoms%nlod,atoms%nlod,atoms%ntype]
            dimsInt=dims
            CALL h5screate_simple_f(3,dims(:3),uloulop21nDataSpaceID,hdfError)
            CALL h5dcreate_f(radialGroupID, "uloulop21n", H5T_NATIVE_DOUBLE, uloulop21nDataSpaceID, uloulop21nDataSetID, hdfError)
            CALL h5sclose_f(uloulop21nDataSpaceID,hdfError)
            CALL io_write_real3(uloulop21nDataSetID,[1,1,1],dimsInt(:3),denCoeffsOffDiag%uloulop21n)
            CALL h5dclose_f(uloulop21nDataSetID, hdfError)
         ENDIF
         !--> End: Scalar products

         CALL h5gclose_f(radialGroupID, hdfError)
      ENDIF
      !--> End: Radial Function output


   END SUBROUTINE writeGreensFData

#endif

END MODULE m_greensf_io

