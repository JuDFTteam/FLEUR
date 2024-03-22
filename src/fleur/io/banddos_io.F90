!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module is used to write out data that can be used to generate augmented
! DOS or bandstructure plots. For the augmentation additional data, e.g., weights
! for the orbital character of states, is stored.
!
!                                          GM' 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_banddos_io

#ifdef CPP_HDF

   USE hdf5
   USE m_hdf_tools

   IMPLICIT NONE

   PUBLIC openBandDOSFile, closeBandDOSFile, writeBandData, writedosData

   CONTAINS

   SUBROUTINE openBandDOSFile(fileID, input, atoms, cell, kpts, sym, banddos, eFermiPrev)

      USE m_types_input
      USE m_types_atoms
      USE m_types_cell
      USE m_types_kpts
      USE m_types_sym
      USE m_types_banddos

      USE hdf5
      !USE m_cdn_io

      TYPE(t_input),   INTENT(IN)  :: input
      TYPE(t_atoms),   INTENT(IN)  :: atoms
      TYPE(t_cell),    INTENT(IN)  :: cell
      TYPE(t_kpts),    INTENT(IN)  :: kpts
      TYPE(t_sym),     INTENT(IN)  :: sym
      TYPE(t_banddos), INTENT(IN)  :: banddos

      INTEGER(HID_T),       INTENT(OUT) :: fileID
      REAL,INTENT(IN)           :: eFermiPrev


      LOGICAL           :: l_exist
      CHARACTER(LEN=30) :: filename
      INTEGER(HID_T)    :: metaGroupID
      INTEGER(HID_T)    :: generalGroupID
      INTEGER(HID_T)    :: cellGroupID
      INTEGER(HID_T)    :: atomsGroupID
      INTEGER(HID_T)    :: kptsGroupID

      INTEGER(HID_T)    :: stringTypeID
      INTEGER(SIZE_T)   :: stringLength

      INTEGER(HID_T)    :: bravaisMatrixSpaceID, bravaisMatrixSetID
      INTEGER(HID_T)    :: reciprocalCellSpaceID, reciprocalCellSetID
      INTEGER(HID_T)    :: rotMatricesSpaceID, rotMatricesSetID
      INTEGER(HID_T)    :: transVecsSpaceID, transVecsSetID

      INTEGER(HID_T)    :: atomPosSpaceID, atomPosSetID
      INTEGER(HID_T)    :: atomicNumbersSpaceID, atomicNumbersSetID
      INTEGER(HID_T)    :: equivAtomsClassSpaceID, equivAtomsClassSetID

      INTEGER(HID_T)    :: kptCoordSpaceID, kptCoordSetID
      INTEGER(HID_T)    :: kptWeightSpaceID, kptWeightSetID
      INTEGER(HID_T)    :: kptSPLabelsSpaceID, kptSPLabelsSetID
      INTEGER(HID_T)    :: kptsSPIndicesSpaceID, kptsSPIndicesSetID

      INTEGER           :: iType, j, iAtom

      INTEGER           :: hdfError, dimsInt(7)
      INTEGER           :: version

      LOGICAL           :: l_error

      INTEGER           :: atomicNumbers(atoms%nat)
      INTEGER           :: equivAtomsGroup(atoms%nat)

      INTEGER(HSIZE_T)  :: dims(7)

      version = 2
      filename = 'banddos.hdf'

      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
      IF(l_exist) THEN
         CALL system('rm '//TRIM(ADJUSTL(filename)))
      END IF

      CALL h5fcreate_f(TRIM(ADJUSTL(filename)), H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

      CALL h5gcreate_f(fileID, '/meta', metaGroupID, hdfError)
      CALL io_write_attint0(metaGroupID,'version',version)
      CALL h5gclose_f(metaGroupID, hdfError)



      CALL h5gcreate_f(fileID, '/general', generalGroupID, hdfError)
      CALL io_write_attint0(generalGroupID,'spins',input%jspins)
      CALL io_write_attreal0(generalGroupID,'lastFermiEnergy',eFermiPrev)
      CALL io_write_attreal0(generalGroupID,'valenceCharge',input%zelec)
      CALL io_write_attlog0(generalGroupID,'bandUnfolding',banddos%unfoldband)
      CALL h5gclose_f(generalGroupID, hdfError)

      CALL h5gcreate_f(fileID, '/cell', cellGroupID, hdfError)

      CALL io_write_attint0(cellGroupID,'nop',sym%nop)
      CALL io_write_attint0(cellGroupID,'nop2',sym%nop2)
      CALL io_write_attlog0(cellGroupID,'film',input%film)

      dims(:2)=(/3,3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),bravaisMatrixSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "bravaisMatrix", H5T_NATIVE_DOUBLE, bravaisMatrixSpaceID, bravaisMatrixSetID, hdfError)
      CALL h5sclose_f(bravaisMatrixSpaceID,hdfError)
      CALL io_write_real2(bravaisMatrixSetID,(/1,1/),dimsInt(:2),"amat",cell%amat)
      CALL h5dclose_f(bravaisMatrixSetID, hdfError)

      dims(:2)=(/3,3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),reciprocalCellSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "reciprocalCell", H5T_NATIVE_DOUBLE, reciprocalCellSpaceID, reciprocalCellSetID, hdfError)
      CALL h5sclose_f(reciprocalCellSpaceID,hdfError)
      CALL io_write_real2(reciprocalCellSetID,(/1,1/),dimsInt(:2),"bmar",cell%bmat)
      CALL h5dclose_f(reciprocalCellSetID, hdfError)

      dims(:3)=(/3,3,sym%nop/)
      dimsInt=dims
      CALL h5screate_simple_f(3,dims(:3),rotMatricesSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "symRotMatrices", H5T_NATIVE_INTEGER, rotMatricesSpaceID, rotMatricesSetID, hdfError)
      CALL h5sclose_f(rotMatricesSpaceID,hdfError)
      CALL io_write_integer3(rotMatricesSetID,(/1,1,1/),dimsInt(:3),"mrot",sym%mrot(:,:,:sym%nop))
      CALL h5dclose_f(rotMatricesSetID, hdfError)

      dims(:2)=(/3,sym%nop/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),transVecsSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "symTransVecs", H5T_NATIVE_DOUBLE, transVecsSpaceID, transVecsSetID, hdfError)
      CALL h5sclose_f(transVecsSpaceID,hdfError)
      CALL io_write_real2(transVecsSetID,(/1,1/),dimsInt(:2),"tau",sym%tau(:,:sym%nop))
      CALL h5dclose_f(transVecsSetID, hdfError)

      CALL h5gclose_f(cellGroupID, hdfError)

      iAtom = 0
      DO iType = 1, atoms%ntype
         DO j = 1, atoms%neq(iType)
            iAtom = iAtom + 1
            atomicNumbers(iAtom) = atoms%nz(iType)
            equivAtomsGroup(iAtom) = iType
         END DO
      END DO

      CALL h5gcreate_f(fileID, '/atoms', atomsGroupID, hdfError)
      CALL io_write_attint0(atomsGroupID,'nAtoms',atoms%nat)
      CALL io_write_attint0(atomsGroupID,'nTypes',atoms%ntype)

      dims(:2)=(/3,atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),atomPosSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "positions", H5T_NATIVE_DOUBLE, atomPosSpaceID, atomPosSetID, hdfError)
      CALL h5sclose_f(atomPosSpaceID,hdfError)
      CALL io_write_real2(atomPosSetID,(/1,1/),dimsInt(:2),"taual",atoms%taual)
      CALL h5dclose_f(atomPosSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "atomicNumbers", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),"atomicNumbers",atomicNumbers)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),equivAtomsClassSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "equivAtomsGroup", H5T_NATIVE_INTEGER, equivAtomsClassSpaceID, equivAtomsClassSetID, hdfError)
      CALL h5sclose_f(equivAtomsClassSpaceID,hdfError)
      CALL io_write_integer1(equivAtomsClassSetID,(/1/),dimsInt(:1),"equivAtomsGroup",equivAtomsGroup)
      CALL h5dclose_f(equivAtomsClassSetID, hdfError)

      CALL h5gclose_f(atomsGroupID, hdfError)

      CALL h5gcreate_f(fileID, '/kpts', kptsGroupID, hdfError)

      CALL io_write_attint0(kptsGroupID,'nkpt',kpts%nkpt)
      CALL io_write_attint0(kptsGroupID,'nSpecialPoints',kpts%numSpecialPoints)

      dims(:2)=(/3,kpts%nkpt/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),kptCoordSpaceID,hdfError)
      CALL h5dcreate_f(kptsGroupID, "coordinates", H5T_NATIVE_DOUBLE, kptCoordSpaceID, kptCoordSetID, hdfError)
      CALL h5sclose_f(kptCoordSpaceID,hdfError)
      CALL io_write_real2(kptCoordSetID,(/1,1/),dimsInt(:2),"bk",kpts%bk)
      CALL h5dclose_f(kptCoordSetID, hdfError)

      dims(:1)=(/kpts%nkpt/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),kptWeightSpaceID,hdfError)
      CALL h5dcreate_f(kptsGroupID, "weights", H5T_NATIVE_DOUBLE, kptWeightSpaceID, kptWeightSetID, hdfError)
      CALL h5sclose_f(kptWeightSpaceID,hdfError)
      CALL io_write_real1(kptWeightSetID,(/1/),dimsInt(:1),"wtkpt",kpts%wtkpt)
      CALL h5dclose_f(kptWeightSetID, hdfError)

      IF (ALLOCATED(kpts%specialPointIndices)) THEN
         stringLength = LEN(kpts%specialPointNames(:))
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeID, hdfError)
         CALL h5tset_size_f(stringTypeID, stringLength, hdfError)
         CALL h5tset_strpad_f(stringTypeID, H5T_STR_SPACEPAD_F, hdfError)
         CALL h5tset_cset_f(stringTypeID, H5T_CSET_ASCII_F, hdfError)
         dims(:1)=(/kpts%numSpecialPoints/)
         dimsInt=dims
         CALL h5screate_simple_f(1,dims(:1),kptSPLabelsSpaceID,hdfError)
         CALL h5dcreate_f(kptsGroupID, "specialPointLabels", stringTypeID, kptSPLabelsSpaceID, kptSPLabelsSetID, hdfError)
         CALL h5tclose_f(stringTypeID,hdfError)
         CALL h5sclose_f(kptSPLabelsSpaceID,hdfError)
         CALL io_write_string1(kptSPLabelsSetID,dimsInt(:1),LEN(kpts%specialPointNames(:)),kpts%specialPointNames)
         CALL h5dclose_f(kptSPLabelsSetID, hdfError)

         dims(:1)=(/kpts%numSpecialPoints/)
         dimsInt=dims
         CALL h5screate_simple_f(1,dims(:1),kptsSPIndicesSpaceID,hdfError)
         CALL h5dcreate_f(kptsGroupID, "specialPointIndices", H5T_NATIVE_INTEGER, kptsSPIndicesSpaceID, kptsSPIndicesSetID, hdfError)
         CALL h5sclose_f(kptsSPIndicesSpaceID,hdfError)
         CALL io_write_integer1(kptsSPIndicesSetID,(/1/),dimsInt(:1),"indices",kpts%specialPointIndices)
         CALL h5dclose_f(kptsSPIndicesSetID, hdfError)
      END IF

      CALL h5gclose_f(kptsGroupID, hdfError)

   END SUBROUTINE

   SUBROUTINE closeBandDOSFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE

   SUBROUTINE writeBandData(fileID,kpts,name_of_dos,weight_name,weight_eig,eig)
      USE m_types_kpts
      character(len=*),intent(in) :: name_of_dos
      character(len=*),intent(in) :: weight_name
      real,intent(in)             :: weight_eig(:,:,:)
      real,optional,intent(in)    :: eig(:,:,:)
      type(t_kpts),intent(in)     :: kpts


      INTEGER(HID_T),  INTENT(IN) :: fileID
      INTEGER                     :: n

      INTEGER(HID_T)    :: BSGroupID,groupID
      INTEGER           :: hdfError

      if (io_groupexists(fileID,name_of_dos)) THEN
        call io_gopen(fileID,name_of_dos,groupID)
      ELSE
        call h5gcreate_f(fileID, name_of_dos, GroupID, hdfError)
      endif
      if (io_groupexists(GroupID,"BS")) THEN
        call io_gopen(GroupID,"BS",BSgroupID)
      ELSE
        CALL h5gcreate_f(GroupID, "BS", BSGroupID, hdfError)
      endif
      IF (PRESENT(eig)) THEN
         if (.not.io_dataexists(BSGroupID,"eigenvalues")) call io_write_var(BSGroupID,"eigenvalues",eig)
      END IF
      call io_write_var(BSGroupID,weight_name,weight_eig(:,:,:))
      CALL h5gclose_f(BSGroupID, hdfError)
      CALL h5gclose_f(GroupID, hdfError)
   END SUBROUTINE

   SUBROUTINE writedosData(fileID,name_of_dos,e_grid,weight_name,dos)
     character(len=*),intent(in) :: name_of_dos
     character(len=*),intent(in) :: weight_name
     real,intent(in):: e_grid(:),dos(:,:)
      INTEGER(HID_T),  INTENT(IN) :: fileID

      INTEGER                     :: n

      INTEGER(HID_T)    :: DOSGroupID,groupID
      INTEGER           :: hdfError

      if (io_groupexists(fileID,name_of_dos)) THEN
        call io_gopen(fileID,name_of_dos,groupID)
      ELSE
        call h5gcreate_f(fileID, name_of_dos, GroupID, hdfError)
      endif
      if (io_groupexists(groupID,"DOS")) then
        call io_gopen(groupID,"DOS",DOSGroupID)
      ELSE
        CALL h5gcreate_f(GroupID, "DOS", DOSGroupID, hdfError)
      endif
      if (.not.io_dataexists(DOSGroupID,"energyGrid")) call io_write_var(DOSGroupID,"energyGrid",e_grid)
      print *,name_of_dos,weight_name
      call io_write_var(DOSGroupID,weight_name,dos(:,:))
      CALL h5gclose_f(DOSGroupID, hdfError)
      CALL h5gclose_f(GroupID, hdfError)

   END SUBROUTINE

   SUBROUTINE writeEVData(fileID,name_of_dos,weight_name,eig,weight_eig)
      character(len=*),intent(in) :: name_of_dos
      character(len=*),intent(in) :: weight_name
      real,intent(in)             :: eig(:,:,:),weight_eig(:,:,:)

      INTEGER(HID_T),  INTENT(IN) :: fileID
      INTEGER                     :: n

      INTEGER(HID_T)    :: EVGroupID,groupID
      INTEGER           :: hdfError

      if (io_groupexists(fileID,name_of_dos)) THEN
        call io_gopen(fileID,name_of_dos,groupID)
      ELSE
        call h5gcreate_f(fileID, name_of_dos, GroupID, hdfError)
      endif
      if (io_groupexists(GroupID,"EV")) THEN
        call io_gopen(GroupID,"EV",EVgroupID)
      ELSE
        CALL h5gcreate_f(GroupID, "EV", EVGroupID, hdfError)
      endif
      if (.not.io_dataexists(EVGroupID,"eigenvalues")) call io_write_var(EVGroupID,"eigenvalues",eig)
      call io_write_var(EVGroupID,weight_name,weight_eig(:,:,:))
      CALL h5gclose_f(EVGroupID, hdfError)
      CALL h5gclose_f(GroupID, hdfError)
   END SUBROUTINE

   SUBROUTINE io_write_string1(datasetID,dims,stringLength,dataArray)

      USE hdf5
      USE m_hdf_tools4

      IMPLICIT NONE

      INTEGER(HID_T),              INTENT(IN) :: datasetID
      INTEGER,                     INTENT(IN) :: dims(1)
      INTEGER,                     INTENT(IN) :: stringLength
      CHARACTER(LEN=stringLength), INTENT(IN) :: dataArray(:)

      INTEGER          :: hdfError
      INTEGER(HID_T)   :: dataspaceID, memSpaceID
      INTEGER(HID_T)   :: stringTypeID
      INTEGER(HID_t)   :: trans
      INTEGER(HSIZE_t) :: memOffset(1), fncount(1)
      INTEGER(HSIZE_t) :: dimsHDF(1)
      INTEGER(SIZE_T)  :: stringLengthHDF

      stringLengthHDF = stringLength
      dimsHDF(:) = dims(:)
      memOffset(:) = 0
      fnCount(:) = dims(:)

      trans = gettransprop()

      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeID, hdfError)
      CALL h5tset_size_f(stringTypeID, stringLengthHDF, hdfError)
      CALL h5tset_strpad_f(stringTypeID, H5T_STR_SPACEPAD_F, hdfError)
      CALL h5tset_cset_f(stringTypeID, H5T_CSET_ASCII_F, hdfError)

      CALL h5dget_space_f(datasetID,dataspaceID,hdfError)
      CALL h5sselect_hyperslab_f(dataspaceID,H5S_SELECT_SET_F,memOffset,fncount,hdfError)
      CALL h5screate_simple_f(1,dimsHDF,memSpaceID,hdfError)
      CALL h5dwrite_f(datasetID,stringTypeID,dataArray,dimsHDF,hdfError,memSpaceID,dataspaceID,trans)
      CALL h5sclose_f(memSpaceID,hdfError)
      CALL h5sclose_f(dataspaceID,hdfError)
      CALL cleartransprop(trans)

      CALL h5tclose_f(stringTypeID,hdfError)

      CALL io_check("io_write_string1 !",hdfError)

   END SUBROUTINE io_write_string1

#endif

END MODULE m_banddos_io
