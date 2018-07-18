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

   PUBLIC openBandDOSFile, closeBandDOSFile, writeBandDOSData

   CONTAINS

   SUBROUTINE openBandDOSFile(fileID, input, atoms, cell, kpts)

      USE m_types
      USE m_cdn_io

      TYPE(t_input), INTENT(IN)  :: input
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell),  INTENT(IN)  :: cell
      TYPE(t_kpts),  INTENT(IN)  :: kpts

      INTEGER(HID_T),       INTENT(OUT) :: fileID

      LOGICAL           :: l_exist
      CHARACTER(LEN=30) :: filename
      INTEGER(HID_T)    :: metaGroupID
      INTEGER(HID_T)    :: generalGroupID
      INTEGER(HID_T)    :: cellGroupID
      INTEGER(HID_T)    :: atomsGroupID
      INTEGER(HID_T)    :: kptsGroupID

      INTEGER(HID_T)    :: bravaisMatrixSpaceID, bravaisMatrixSetID
      INTEGER(HID_T)    :: reciprocalCellSpaceID, reciprocalCellSetID

      INTEGER(HID_T)    :: atomPosSpaceID, atomPosSetID
      INTEGER(HID_T)    :: atomicNumbersSpaceID, atomicNumbersSetID

      INTEGER(HID_T)    :: kptCoordSpaceID, kptCoordSetID
      INTEGER(HID_T)    :: kptWeightSpaceID, kptWeightSetID

      INTEGER           :: iType, j, iAtom

      INTEGER           :: hdfError, dimsInt(7)
      INTEGER           :: version
      REAL              :: eFermiPrev
      LOGICAL           :: l_error

      INTEGER           :: atomicNumbers(atoms%nat)

      INTEGER(HSIZE_T)  :: dims(7)

      version = 1
      filename = 'banddos.hdf'

      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
      IF(l_exist) THEN
         CALL system('rm '//TRIM(ADJUSTL(filename)))       
      END IF

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
      CALL io_write_attreal0(generalGroupID,'lastFermiEnergy',eFermiPrev)
      CALL h5gclose_f(generalGroupID, hdfError)

      CALL h5gcreate_f(fileID, '/cell', cellGroupID, hdfError)

      dims(:2)=(/3,3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),bravaisMatrixSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "bravaisMatrix", H5T_NATIVE_DOUBLE, bravaisMatrixSpaceID, bravaisMatrixSetID, hdfError)
      CALL h5sclose_f(bravaisMatrixSpaceID,hdfError)
      CALL io_write_real2(bravaisMatrixSetID,(/1,1/),dimsInt(:2),cell%amat)
      CALL h5dclose_f(bravaisMatrixSetID, hdfError)

      dims(:2)=(/3,3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),reciprocalCellSpaceID,hdfError)
      CALL h5dcreate_f(cellGroupID, "reciprocalCell", H5T_NATIVE_DOUBLE, reciprocalCellSpaceID, reciprocalCellSetID, hdfError)
      CALL h5sclose_f(reciprocalCellSpaceID,hdfError)
      CALL io_write_real2(reciprocalCellSetID,(/1,1/),dimsInt(:2),cell%bmat)
      CALL h5dclose_f(reciprocalCellSetID, hdfError)

      CALL h5gclose_f(cellGroupID, hdfError)

      iAtom = 0
      DO iType = 1, atoms%ntype
         DO j = 1, atoms%neq(iType)
            iAtom = iAtom + 1
            atomicNumbers(iAtom) = atoms%nz(iType)
         END DO
      END DO

      CALL h5gcreate_f(fileID, '/atoms', atomsGroupID, hdfError)
      CALL io_write_attint0(atomsGroupID,'nAtoms',atoms%nat)

      dims(:2)=(/3,atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),atomPosSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "positions", H5T_NATIVE_DOUBLE, atomPosSpaceID, atomPosSetID, hdfError)
      CALL h5sclose_f(atomPosSpaceID,hdfError)
      CALL io_write_real2(atomPosSetID,(/1,1/),dimsInt(:2),atoms%taual)
      CALL h5dclose_f(atomPosSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),atomicNumbersSpaceID,hdfError)
      CALL h5dcreate_f(atomsGroupID, "atomicNumbers", H5T_NATIVE_INTEGER, atomicNumbersSpaceID, atomicNumbersSetID, hdfError)
      CALL h5sclose_f(atomicNumbersSpaceID,hdfError)
      CALL io_write_integer1(atomicNumbersSetID,(/1/),dimsInt(:1),atomicNumbers)
      CALL h5dclose_f(atomicNumbersSetID, hdfError)

      CALL h5gclose_f(atomsGroupID, hdfError)

      CALL h5gcreate_f(fileID, '/kpts', kptsGroupID, hdfError)

      CALL io_write_attint0(kptsGroupID,'nkpt',kpts%nkpt)

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

      CALL h5gclose_f(kptsGroupID, hdfError)

   END SUBROUTINE

   SUBROUTINE closeBandDOSFile(fileID)

      INTEGER(HID_T), INTENT(IN)  :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE

   SUBROUTINE writeBandDOSData(fileID,input,atoms,cell,kpts,results,banddos,dos,vacuum)

      USE m_types

      TYPE(t_input),   INTENT(IN) :: input
      TYPE(t_atoms),   INTENT(IN) :: atoms
      TYPE(t_cell),    INTENT(IN) :: cell
      TYPE(t_kpts),    INTENT(IN) :: kpts
      TYPE(t_results), INTENT(IN) :: results
      TYPE(t_banddos), INTENT(IN) :: banddos
      TYPE(t_dos),     INTENT(IN) :: dos
      TYPE(t_vacuum),  INTENT(IN) :: vacuum

      INTEGER                     :: neigd

      INTEGER(HID_T),  INTENT(IN) :: fileID

      INTEGER(HID_T)    :: eigenvaluesGroupID

      INTEGER(HID_T)    :: eigenvaluesSpaceID, eigenvaluesSetID
      INTEGER(HID_T)    :: numFoundEigsSpaceID, numFoundEigsSetID
      INTEGER(HID_T)    :: lLikeChargeSpaceID, lLikeChargeSetID

      INTEGER           :: hdfError, dimsInt(7)

      INTEGER(HSIZE_T)  :: dims(7)

      neigd = MAXVAL(results%neig(:,:))

      CALL h5gcreate_f(fileID, '/eigenvalues', eigenvaluesGroupID, hdfError)

      CALL io_write_attint0(eigenvaluesGroupID,'neigd',neigd)
      CALL io_write_attint0(eigenvaluesGroupID,'maxL',3)

      dims(:2)=(/kpts%nkpt,input%jspins/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),numFoundEigsSpaceID,hdfError)
      CALL h5dcreate_f(eigenvaluesGroupID, "numFoundEigenvals", H5T_NATIVE_INTEGER, numFoundEigsSpaceID, numFoundEigsSetID, hdfError)
      CALL h5sclose_f(numFoundEigsSpaceID,hdfError)
      CALL io_write_integer2(numFoundEigsSetID,(/1,1/),dimsInt(:2),results%neig)
      CALL h5dclose_f(numFoundEigsSetID, hdfError)

      dims(:3)=(/neigd,kpts%nkpt,input%jspins/)
      dimsInt = dims
      CALL h5screate_simple_f(3,dims(:3),eigenvaluesSpaceID,hdfError)
      CALL h5dcreate_f(eigenvaluesGroupID, "eigenvalues", H5T_NATIVE_DOUBLE, eigenvaluesSpaceID, eigenvaluesSetID, hdfError)
      CALL h5sclose_f(eigenvaluesSpaceID,hdfError)
      CALL io_write_real3(eigenvaluesSetID,(/1,1,1/),dimsInt(:3),results%eig(:neigd,:,:))
      CALL h5dclose_f(eigenvaluesSetID, hdfError)

      dims(:5)=(/4,atoms%ntype,neigd,kpts%nkpt,input%jspins/)
      dimsInt = dims
      CALL h5screate_simple_f(5,dims(:5),lLikeChargeSpaceID,hdfError)
      CALL h5dcreate_f(eigenvaluesGroupID, "lLikeCharge", H5T_NATIVE_DOUBLE, lLikeChargeSpaceID, lLikeChargeSetID, hdfError)
      CALL h5sclose_f(lLikeChargeSpaceID,hdfError)
      CALL io_write_real5(lLikeChargeSetID,(/1,1,1,1,1/),dimsInt(:5),dos%qal(0:3,:,:neigd,:,:))
      CALL h5dclose_f(lLikeChargeSetID, hdfError)

      CALL h5gclose_f(eigenvaluesGroupID, hdfError)

   END SUBROUTINE

#endif

END MODULE m_banddos_io
