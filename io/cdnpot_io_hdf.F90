!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnpot_io_hdf

   USE m_types
   USE m_juDFT
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools
#endif
   IMPLICIT NONE

   PRIVATE
#ifdef CPP_HDF
   PUBLIC openCDN_HDF, openPOT_HDF, closeCDNPOT_HDF
   PUBLIC writeStarsHDF, readStarsHDF
   PUBLIC writeLatharmsHDF, readLatharmsHDF
   PUBLIC writeStructureHDF, readStructureHDF
   PUBLIC writeStepfunctionHDF, readStepfunctionHDF, peekStepfunctionHDF
   PUBLIC writeDensityHDF, readDensityHDF
   PUBLIC writePotentialHDF, readPotentialHDF
   PUBLIC writeCoreDensityHDF, readCoreDensityHDF
   PUBLIC writeCDNHeaderData, writePOTHeaderData
   PUBLIC isCoreDensityPresentHDF, deleteDensityEntryHDF
   PUBLIC isDensityEntryPresentHDF, isPotentialEntryPresentHDF
   PUBLIC peekDensityEntryHDF
#endif

   PUBLIC DENSITY_TYPE_UNDEFINED_const
   PUBLIC DENSITY_TYPE_IN_const, DENSITY_TYPE_OUT_const
   PUBLIC DENSITY_TYPE_NOCO_IN_const, DENSITY_TYPE_NOCO_OUT_const
   PUBLIC DENSITY_TYPE_PRECOND_const
   PUBLIC POTENTIAL_TYPE_IN_const, POTENTIAL_TYPE_OUT_const

   INTEGER,          PARAMETER :: DENSITY_TYPE_UNDEFINED_const = 0
   INTEGER,          PARAMETER :: DENSITY_TYPE_IN_const = 1
   INTEGER,          PARAMETER :: DENSITY_TYPE_OUT_const = 2
   INTEGER,          PARAMETER :: DENSITY_TYPE_NOCO_IN_const = 3
   INTEGER,          PARAMETER :: DENSITY_TYPE_NOCO_OUT_const = 4
   INTEGER,          PARAMETER :: DENSITY_TYPE_PRECOND_const = 5

   INTEGER,          PARAMETER :: POTENTIAL_TYPE_IN_const = 1
   INTEGER,          PARAMETER :: POTENTIAL_TYPE_OUT_const = 2

   INTEGER,          PARAMETER :: FILE_FORMAT_VERSION_const = 28

   CONTAINS

#ifdef CPP_HDF

   SUBROUTINE openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

      INTEGER(HID_T), INTENT(OUT) :: fileID
      INTEGER, INTENT(OUT) :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER, INTENT(OUT) :: currentStepfunctionIndex, readDensityIndex,lastDensityIndex

      INTEGER(HID_T) :: generalGroupID
      INTEGER        :: hdfError, fileFormatVersion
      LOGICAL        :: l_exist

      currentStarsIndex = 0
      currentLatharmsIndex = 0
      currentStructureIndex = 0
      currentStepfunctionIndex = 0
      readDensityIndex = 0
      lastDensityIndex = 0
      fileFormatVersion = 0

      INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
      IF(l_exist) THEN ! only open file
         CALL h5fopen_f('cdn.hdf', H5F_ACC_RDWR_F, fileID, hdfError, H5P_DEFAULT_F)

         CALL h5gopen_f(fileID, '/general', generalGroupID, hdfError)
         ! read in primary attributes from the header '/general'
         CALL io_read_attint0(generalGroupID,'currentStarsIndex',currentStarsIndex)
         CALL io_read_attint0(generalGroupID,'currentLatharmsIndex',currentLatharmsIndex)
         CALL io_read_attint0(generalGroupID,'currentStructureIndex',currentStructureIndex)
         CALL io_read_attint0(generalGroupID,'currentStepfunctionIndex',currentStepfunctionIndex)
         CALL io_read_attint0(generalGroupID,'readDensityIndex',readDensityIndex)
         CALL io_read_attint0(generalGroupID,'lastDensityIndex',lastDensityIndex)
         CALL io_read_attint0(generalGroupID,'fileFormatVersion',fileFormatVersion)

         CALL h5gclose_f(generalGroupID, hdfError)
         IF(fileFormatVersion.NE.FILE_FORMAT_VERSION_const) THEN
            WRITE(*,'(a,i4)') 'cdn.hdf has file format version ', fileFormatVersion
            CALL juDFT_error('cdn.hdf file format not readable.' ,calledby ="openCDN_HDF")
         END IF
      ELSE ! create file
         CALL h5fcreate_f('cdn.hdf', H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

         CALL h5gcreate_f(fileID, '/general', generalGroupID, hdfError)
         ! write initial values to primary attributes in the header '/general'
         CALL io_write_attint0(generalGroupID,'currentStarsIndex',currentStarsIndex)
         CALL io_write_attint0(generalGroupID,'currentLatharmsIndex',currentLatharmsIndex)
         CALL io_write_attint0(generalGroupID,'currentStructureIndex',currentStructureIndex)
         CALL io_write_attint0(generalGroupID,'currentStepfunctionIndex',currentStepfunctionIndex)
         CALL io_write_attint0(generalGroupID,'readDensityIndex',readDensityIndex)
         CALL io_write_attint0(generalGroupID,'lastDensityIndex',lastDensityIndex)
         CALL io_write_attint0(generalGroupID,'fileFormatVersion',FILE_FORMAT_VERSION_const)

         CALL h5gclose_f(generalGroupID, hdfError)
      END IF

   END SUBROUTINE openCDN_HDF


   SUBROUTINE openPOT_HDF(fileID,currentStarsIndex,currentLatharmsIndex,&
                          currentStructureIndex,currentStepfunctionIndex)

      INTEGER(HID_T), INTENT(OUT) :: fileID
      INTEGER, INTENT(OUT) :: currentStarsIndex,currentLatharmsIndex
      INTEGER, INTENT(OUT) :: currentStructureIndex, currentStepfunctionIndex

      INTEGER(HID_T) :: generalGroupID
      INTEGER        :: hdfError, fileFormatVersion
      LOGICAL        :: l_exist

      currentStarsIndex = 0
      currentLatharmsIndex = 0
      currentStructureIndex = 0
      currentStepfunctionIndex = 0
      fileFormatVersion = 0

      INQUIRE(FILE='pot.hdf',EXIST=l_exist)
      IF(l_exist) THEN ! only open file
         CALL h5fopen_f('pot.hdf', H5F_ACC_RDWR_F, fileID, hdfError, H5P_DEFAULT_F)

         CALL h5gopen_f(fileID, '/general', generalGroupID, hdfError)
         ! read in primary attributes from the header '/general'
         CALL io_read_attint0(generalGroupID,'currentStarsIndex',currentStarsIndex)
         CALL io_read_attint0(generalGroupID,'currentLatharmsIndex',currentLatharmsIndex)
         CALL io_read_attint0(generalGroupID,'currentStructureIndex',currentStructureIndex)
         CALL io_read_attint0(generalGroupID,'currentStepfunctionIndex',currentStepfunctionIndex)
         CALL io_read_attint0(generalGroupID,'fileFormatVersion',fileFormatVersion)

         CALL h5gclose_f(generalGroupID, hdfError)
         IF(fileFormatVersion.NE.FILE_FORMAT_VERSION_const) THEN
            WRITE(*,'(a,i4)') 'pot.hdf has file format version ', fileFormatVersion
            CALL juDFT_error('pot.hdf file format not readable.' ,calledby ="openPOT_HDF")
         END IF
      ELSE ! create file
         CALL h5fcreate_f('pot.hdf', H5F_ACC_TRUNC_F, fileID, hdfError, H5P_DEFAULT_F, H5P_DEFAULT_F)

         CALL h5gcreate_f(fileID, '/general', generalGroupID, hdfError)
         ! write initial values to primary attributes in the header '/general'
         CALL io_write_attint0(generalGroupID,'currentStarsIndex',currentStarsIndex)
         CALL io_write_attint0(generalGroupID,'currentLatharmsIndex',currentLatharmsIndex)
         CALL io_write_attint0(generalGroupID,'currentStructureIndex',currentStructureIndex)
         CALL io_write_attint0(generalGroupID,'currentStepfunctionIndex',currentStepfunctionIndex)
         CALL io_write_attint0(generalGroupID,'fileFormatVersion',FILE_FORMAT_VERSION_const)

         CALL h5gclose_f(generalGroupID, hdfError)
      END IF

   END SUBROUTINE openPOT_HDF

   SUBROUTINE closeCDNPOT_HDF(fileID)

      INTEGER(HID_T), INTENT(IN) :: fileID

      INTEGER hdfError

      CALL h5fclose_f(fileID, hdfError)

   END SUBROUTINE closeCDNPOT_HDF

   SUBROUTINE writeCDNHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                 currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

      INTEGER(HID_T), INTENT(IN) :: fileID
      INTEGER, INTENT(IN)        :: currentStarsIndex
      INTEGER, INTENT(IN)        :: currentLatharmsIndex
      INTEGER, INTENT(IN)        :: currentStructureIndex
      INTEGER, INTENT(IN)        :: currentStepfunctionIndex
      INTEGER, INTENT(IN)        :: readDensityIndex
      INTEGER, INTENT(IN)        :: lastDensityIndex

      INTEGER(HID_T) :: generalGroupID
      INTEGER        :: hdfError

      CALL h5gopen_f(fileID, '/general', generalGroupID, hdfError)

      CALL io_write_attint0(generalGroupID,'currentStarsIndex',currentStarsIndex)
      CALL io_write_attint0(generalGroupID,'currentLatharmsIndex',currentLatharmsIndex)
      CALL io_write_attint0(generalGroupID,'currentStructureIndex',currentStructureIndex)
      CALL io_write_attint0(generalGroupID,'currentStepfunctionIndex',currentStepfunctionIndex)
      CALL io_write_attint0(generalGroupID,'readDensityIndex',readDensityIndex)
      CALL io_write_attint0(generalGroupID,'lastDensityIndex',lastDensityIndex)
      CALL io_write_attint0(generalGroupID,'fileFormatVersion',FILE_FORMAT_VERSION_const)

      CALL h5gclose_f(generalGroupID, hdfError)

   END SUBROUTINE writeCDNHeaderData

   SUBROUTINE writePOTHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,&
                                 currentStructureIndex,currentStepfunctionIndex)

      INTEGER(HID_T), INTENT(IN) :: fileID
      INTEGER, INTENT(IN)        :: currentStarsIndex
      INTEGER, INTENT(IN)        :: currentLatharmsIndex
      INTEGER, INTENT(IN)        :: currentStructureIndex
      INTEGER, INTENT(IN)        :: currentStepfunctionIndex

      INTEGER(HID_T) :: generalGroupID
      INTEGER        :: hdfError

      CALL h5gopen_f(fileID, '/general', generalGroupID, hdfError)

      CALL io_write_attint0(generalGroupID,'currentStarsIndex',currentStarsIndex)
      CALL io_write_attint0(generalGroupID,'currentLatharmsIndex',currentLatharmsIndex)
      CALL io_write_attint0(generalGroupID,'currentStructureIndex',currentStructureIndex)
      CALL io_write_attint0(generalGroupID,'currentStepfunctionIndex',currentStepfunctionIndex)
      CALL io_write_attint0(generalGroupID,'fileFormatVersion',FILE_FORMAT_VERSION_const)

      CALL h5gclose_f(generalGroupID, hdfError)

   END SUBROUTINE writePOTHeaderData

   SUBROUTINE writeStarsHDF(fileID, starsIndex, stars)

      INTEGER(HID_T), INTENT(IN) :: fileID
      INTEGER,        INTENT(IN) :: starsIndex
      TYPE(t_stars),  INTENT(IN) :: stars

      INTEGER(HID_T)            :: groupID
      INTEGER                   :: hdfError, ft2_gf_dim, dimsInt(7)
      CHARACTER(LEN=30)         :: groupName
      INTEGER(HSIZE_T)          :: dims(7)
      LOGICAL        :: l_exist

      INTEGER(HID_T)                   :: kv3SpaceID, kv3SetID
      INTEGER(HID_T)                   :: kv2SpaceID, kv2SetID
      INTEGER(HID_T)                   :: sk3SpaceID, sk3SetID
      INTEGER(HID_T)                   :: sk2SpaceID, sk2SetID
      INTEGER(HID_T)                   :: igSpaceID, igSetID
      INTEGER(HID_T)                   :: ig2SpaceID, ig2SetID
      INTEGER(HID_T)                   :: nstrSpaceID, nstrSetID
      INTEGER(HID_T)                   :: nstr2SpaceID, nstr2SetID
      INTEGER(HID_T)                   :: phi2SpaceID, phi2SetID
      INTEGER(HID_T)                   :: rgphsSpaceID, rgphsSetID
      INTEGER(HID_T)                   :: igfftSpaceID, igfftSetID
      INTEGER(HID_T)                   :: igfft2SpaceID, igfft2SetID
      INTEGER(HID_T)                   :: pgfftSpaceID, pgfftSetID
      INTEGER(HID_T)                   :: pgfft2SpaceID, pgfft2SetID
      INTEGER(HID_T)                   :: ft2_gfxSpaceID, ft2_gfxSetID
      INTEGER(HID_T)                   :: ft2_gfySpaceID, ft2_gfySetID

      WRITE(groupname,'(a,i0)') '/stars-', starsIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (l_exist) THEN
         CALL juDFT_error('stars entry '//TRIM(ADJUSTL(groupName))//' already exists.' ,calledby ="writeStarsHDF")
      END IF

      CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      ft2_gf_dim = SIZE(stars%ft2_gfx,1)

      CALL io_write_attreal0(groupID,'gmax',stars%gmax)
      CALL io_write_attreal0(groupID,'gmaxInit',stars%gmaxInit)
      CALL io_write_attint0(groupID,'ng3',stars%ng3)
      CALL io_write_attint0(groupID,'ng2',stars%ng2)
      CALL io_write_attint0(groupID,'mx1',stars%mx1)
      CALL io_write_attint0(groupID,'mx2',stars%mx2)
      CALL io_write_attint0(groupID,'mx3',stars%mx3)
      CALL io_write_attint0(groupID,'kimax',stars%kimax)
      CALL io_write_attint0(groupID,'kimax2',stars%kimax2)
      CALL io_write_attint0(groupID,'kq1_fft',stars%kq1_fft)
      CALL io_write_attint0(groupID,'kq2_fft',stars%kq2_fft)
      CALL io_write_attint0(groupID,'kq3_fft',stars%kq3_fft)
      CALL io_write_attint0(groupID,'kmxq_fft',stars%kmxq_fft)
      CALL io_write_attint0(groupID,'kxc1_fft',stars%kxc1_fft)
      CALL io_write_attint0(groupID,'kxc2_fft',stars%kxc2_fft)
      CALL io_write_attint0(groupID,'kxc3_fft',stars%kxc3_fft)
      CALL io_write_attint0(groupID,'ng3_fft',stars%ng3_fft)
      CALL io_write_attint0(groupID,'kmxxc_fft',stars%kmxxc_fft)
      CALL io_write_attint0(groupID,'nxc3_fft',stars%nxc3_fft)
      CALL io_write_attint0(groupID,'ft2_gf_dim',ft2_gf_dim)

      dims(:2)=(/3,stars%ng3/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),kv3SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "kv3", H5T_NATIVE_INTEGER, kv3SpaceID, kv3SetID, hdfError)
      CALL h5sclose_f(kv3SpaceID,hdfError)
      CALL io_write_integer2(kv3SetID,(/1,1/),dimsInt(:2),stars%kv3)
      CALL h5dclose_f(kv3SetID, hdfError)

      dims(:2)=(/2,stars%ng2/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),kv2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "kv2", H5T_NATIVE_INTEGER, kv2SpaceID, kv2SetID, hdfError)
      CALL h5sclose_f(kv2SpaceID,hdfError)
      CALL io_write_integer2(kv2SetID,(/1,1/),dimsInt(:2),stars%kv2)
      CALL h5dclose_f(kv2SetID, hdfError)

      dims(:1)=(/stars%ng3/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),sk3SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "sk3", H5T_NATIVE_DOUBLE, sk3SpaceID, sk3SetID, hdfError)
      CALL h5sclose_f(sk3SpaceID,hdfError)
      CALL io_write_real1(sk3SetID,(/1/),dimsInt(:1),stars%sk3)
      CALL h5dclose_f(sk3SetID, hdfError)

      dims(:1)=(/stars%ng2/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),sk2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "sk2", H5T_NATIVE_DOUBLE, sk2SpaceID, sk2SetID, hdfError)
      CALL h5sclose_f(sk2SpaceID,hdfError)
      CALL io_write_real1(sk2SetID,(/1/),dimsInt(:1),stars%sk2)
      CALL h5dclose_f(sk2SetID, hdfError)

      dims(:3)=(/2*stars%mx1+1,2*stars%mx2+1,2*stars%mx3+1/)
      dimsInt=dims
      CALL h5screate_simple_f(3,dims(:3),igSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ig", H5T_NATIVE_INTEGER, igSpaceID, igSetID, hdfError)
      CALL h5sclose_f(igSpaceID,hdfError)
      CALL io_write_integer3(igSetID,(/1,1,1/),dimsInt(:3),stars%ig)
      CALL h5dclose_f(igSetID, hdfError)

      dims(:1)=(/stars%ng3/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),ig2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ig2", H5T_NATIVE_INTEGER, ig2SpaceID, ig2SetID, hdfError)
      CALL h5sclose_f(ig2SpaceID,hdfError)
      CALL io_write_integer1(ig2SetID,(/1/),dimsInt(:1),stars%ig2)
      CALL h5dclose_f(ig2SetID, hdfError)

      dims(:1)=(/stars%ng3/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),nstrSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "nstr", H5T_NATIVE_INTEGER, nstrSpaceID, nstrSetID, hdfError)
      CALL h5sclose_f(nstrSpaceID,hdfError)
      CALL io_write_integer1(nstrSetID,(/1/),dimsInt(:1),stars%nstr)
      CALL h5dclose_f(nstrSetID, hdfError)

      dims(:1)=(/stars%ng2/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),nstr2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "nstr2", H5T_NATIVE_INTEGER, nstr2SpaceID, nstr2SetID, hdfError)
      CALL h5sclose_f(nstr2SpaceID,hdfError)
      CALL io_write_integer1(nstr2SetID,(/1/),dimsInt(:1),stars%nstr2)
      CALL h5dclose_f(nstr2SetID, hdfError)

      dims(:1)=(/stars%ng2/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),phi2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "phi2", H5T_NATIVE_DOUBLE, phi2SpaceID, phi2SetID, hdfError)
      CALL h5sclose_f(phi2SpaceID,hdfError)
      CALL io_write_real1(phi2SetID,(/1/),dimsInt(:1),stars%phi2)
      CALL h5dclose_f(phi2SetID, hdfError)

      dims(:4)=(/2,2*stars%mx1+1,2*stars%mx2+1,2*stars%mx3+1/)
      dimsInt=dims
      CALL h5screate_simple_f(4,dims(:4),rgphsSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "rgphs", H5T_NATIVE_DOUBLE, rgphsSpaceID, rgphsSetID, hdfError)
      CALL h5sclose_f(rgphsSpaceID,hdfError)
      CALL io_write_complex3(rgphsSetID,(/-1,1,1,1/),dimsInt(:4),stars%rgphs)
      CALL h5dclose_f(rgphsSetID, hdfError)

      dims(:2)=(/stars%kimax+1,2/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),igfftSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "igfft", H5T_NATIVE_INTEGER, igfftSpaceID, igfftSetID, hdfError)
      CALL h5sclose_f(igfftSpaceID,hdfError)
      CALL io_write_integer2(igfftSetID,(/1,1/),dimsInt(:2),stars%igfft(0:stars%kimax,:))
      CALL h5dclose_f(igfftSetID, hdfError)

      dims(:2)=(/stars%kimax2+1,2/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),igfft2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "igfft2", H5T_NATIVE_INTEGER, igfft2SpaceID, igfft2SetID, hdfError)
      CALL h5sclose_f(igfft2SpaceID,hdfError)
      CALL io_write_integer2(igfft2SetID,(/1,1/),dimsInt(:2),stars%igfft2(0:stars%kimax2,:))
      CALL h5dclose_f(igfft2SetID, hdfError)

      dims(:2)=(/2,stars%kimax+1/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),pgfftSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "pgfft", H5T_NATIVE_DOUBLE, pgfftSpaceID, pgfftSetID, hdfError)
      CALL h5sclose_f(pgfftSpaceID,hdfError)
      CALL io_write_complex1(pgfftSetID,(/-1,1/),dimsInt(:2),stars%pgfft(0:stars%kimax))
      CALL h5dclose_f(pgfftSetID, hdfError)

      dims(:2)=(/2,stars%kimax2+1/)
      dimsInt=dims
      CALL h5screate_simple_f(2,dims(:2),pgfft2SpaceID,hdfError)
      CALL h5dcreate_f(groupID, "pgfft2", H5T_NATIVE_DOUBLE, pgfft2SpaceID, pgfft2SetID, hdfError)
      CALL h5sclose_f(pgfft2SpaceID,hdfError)
      CALL io_write_complex1(pgfft2SetID,(/-1,1/),dimsInt(:2),stars%pgfft2(0:stars%kimax2))
      CALL h5dclose_f(pgfft2SetID, hdfError)

      dims(:1)=(/ft2_gf_dim/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),ft2_gfxSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ft2_gfx", H5T_NATIVE_DOUBLE, ft2_gfxSpaceID, ft2_gfxSetID, hdfError)
      CALL h5sclose_f(ft2_gfxSpaceID,hdfError)
      CALL io_write_real1(ft2_gfxSetID,(/1/),dimsInt(:1),stars%ft2_gfx)
      CALL h5dclose_f(ft2_gfxSetID, hdfError)

      dims(:1)=(/ft2_gf_dim/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),ft2_gfySpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ft2_gfy", H5T_NATIVE_DOUBLE, ft2_gfySpaceID, ft2_gfySetID, hdfError)
      CALL h5sclose_f(ft2_gfySpaceID,hdfError)
      CALL io_write_real1(ft2_gfySetID,(/1/),dimsInt(:1),stars%ft2_gfy)
      CALL h5dclose_f(ft2_gfySetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE writeStarsHDF

   SUBROUTINE readStarsHDF(fileID, starsIndex, stars)

      INTEGER(HID_T), INTENT(IN)    :: fileID
      INTEGER,        INTENT(IN)    :: starsIndex
      TYPE(t_stars),  INTENT(INOUT) :: stars

      INTEGER(HID_T)            :: groupID
      INTEGER                   :: hdfError, ft2_gf_dim
      INTEGER                   :: dimsInt(7)
      CHARACTER(LEN=30)         :: groupName
      LOGICAL                   :: l_exist

      INTEGER(HID_T)                   :: kv3SetID
      INTEGER(HID_T)                   :: kv2SetID
      INTEGER(HID_T)                   :: sk3SetID
      INTEGER(HID_T)                   :: sk2SetID
      INTEGER(HID_T)                   :: igSetID
      INTEGER(HID_T)                   :: ig2SetID
      INTEGER(HID_T)                   :: nstrSetID
      INTEGER(HID_T)                   :: nstr2SetID
      INTEGER(HID_T)                   :: phi2SetID
      INTEGER(HID_T)                   :: rgphsSetID
      INTEGER(HID_T)                   :: igfftSetID
      INTEGER(HID_T)                   :: igfft2SetID
      INTEGER(HID_T)                   :: pgfftSetID
      INTEGER(HID_T)                   :: pgfft2SetID
      INTEGER(HID_T)                   :: ft2_gfxSetID
      INTEGER(HID_T)                   :: ft2_gfySetID

      WRITE(groupname,'(a,i0)') '/stars-', starsIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (.NOT.l_exist) THEN
         CALL juDFT_error('stars entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readStarsHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_read_attreal0(groupID,'gmax',stars%gmax)
      CALL io_read_attreal0(groupID,'gmaxInit',stars%gmaxInit)
      CALL io_read_attint0(groupID,'ng3',stars%ng3)
      CALL io_read_attint0(groupID,'ng2',stars%ng2)
      CALL io_read_attint0(groupID,'mx1',stars%mx1)
      CALL io_read_attint0(groupID,'mx2',stars%mx2)
      CALL io_read_attint0(groupID,'mx3',stars%mx3)
      CALL io_read_attint0(groupID,'kimax',stars%kimax)
      CALL io_read_attint0(groupID,'kimax2',stars%kimax2)
      CALL io_read_attint0(groupID,'kq1_fft',stars%kq1_fft)
      CALL io_read_attint0(groupID,'kq2_fft',stars%kq2_fft)
      CALL io_read_attint0(groupID,'kq3_fft',stars%kq3_fft)
      CALL io_read_attint0(groupID,'kmxq_fft',stars%kmxq_fft)
      CALL io_read_attint0(groupID,'kxc1_fft',stars%kxc1_fft)
      CALL io_read_attint0(groupID,'kxc2_fft',stars%kxc2_fft)
      CALL io_read_attint0(groupID,'kxc3_fft',stars%kxc3_fft)
      CALL io_read_attint0(groupID,'ng3_fft',stars%ng3_fft)
      CALL io_read_attint0(groupID,'kmxxc_fft',stars%kmxxc_fft)
      CALL io_read_attint0(groupID,'nxc3_fft',stars%nxc3_fft)
      CALL io_read_attint0(groupID,'ft2_gf_dim',ft2_gf_dim)

      IF(ALLOCATED(stars%kv3)) DEALLOCATE(stars%kv3)
      IF(ALLOCATED(stars%kv2)) DEALLOCATE(stars%kv2)
      IF(ALLOCATED(stars%sk3)) DEALLOCATE(stars%sk3)
      IF(ALLOCATED(stars%sk2)) DEALLOCATE(stars%sk2)
      IF(ALLOCATED(stars%ig)) DEALLOCATE(stars%ig)
      IF(ALLOCATED(stars%ig2)) DEALLOCATE(stars%ig2)
      IF(ALLOCATED(stars%nstr)) DEALLOCATE(stars%nstr)
      IF(ALLOCATED(stars%nstr2)) DEALLOCATE(stars%nstr2)
      IF(ALLOCATED(stars%phi2)) DEALLOCATE(stars%phi2)
      IF(ALLOCATED(stars%rgphs)) DEALLOCATE(stars%rgphs)
      IF(ALLOCATED(stars%igfft)) DEALLOCATE(stars%igfft)
      IF(ALLOCATED(stars%igfft2)) DEALLOCATE(stars%igfft2)
      IF(ALLOCATED(stars%pgfft)) DEALLOCATE(stars%pgfft)
      IF(ALLOCATED(stars%pgfft2)) DEALLOCATE(stars%pgfft2)
      IF(ALLOCATED(stars%ft2_gfx)) DEALLOCATE(stars%ft2_gfx)
      IF(ALLOCATED(stars%ft2_gfy)) DEALLOCATE(stars%ft2_gfy)

      ALLOCATE(stars%kv3(3,stars%ng3))
      ALLOCATE(stars%kv2(2,stars%ng2))
      ALLOCATE(stars%sk3(stars%ng3))
      ALLOCATE(stars%sk2(stars%ng2))
      ALLOCATE(stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
      ALLOCATE(stars%ig2(stars%ng3))
      ALLOCATE(stars%nstr(stars%ng3))
      ALLOCATE(stars%nstr2(stars%ng2))
      ALLOCATE(stars%phi2(stars%ng2))
      ALLOCATE(stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
      ALLOCATE(stars%igfft(0:stars%kimax,2))
      ALLOCATE(stars%igfft2(0:stars%kimax2,2))
      ALLOCATE(stars%pgfft(0:stars%kimax))
      ALLOCATE(stars%pgfft2(0:stars%kimax2))
      ALLOCATE(stars%ft2_gfx(0:ft2_gf_dim-1))
      ALLOCATE(stars%ft2_gfy(0:ft2_gf_dim-1))

      dimsInt(:2)=(/3,stars%ng3/)
      CALL h5dopen_f(groupID, 'kv3', kv3SetID, hdfError)
      CALL io_read_integer2(kv3SetID,(/1,1/),dimsInt(:2),stars%kv3)
      CALL h5dclose_f(kv3SetID, hdfError)

      dimsInt(:2)=(/2,stars%ng2/)
      CALL h5dopen_f(groupID, 'kv2', kv2SetID, hdfError)
      CALL io_read_integer2(kv2SetID,(/1,1/),dimsInt(:2),stars%kv2)
      CALL h5dclose_f(kv2SetID, hdfError)

      dimsInt(:1)=(/stars%ng3/)
      CALL h5dopen_f(groupID, 'sk3', sk3SetID, hdfError)
      CALL io_read_real1(sk3SetID,(/1/),dimsInt(:1),stars%sk3)
      CALL h5dclose_f(sk3SetID, hdfError)

      dimsInt(:1)=(/stars%ng2/)
      CALL h5dopen_f(groupID, 'sk2', sk2SetID, hdfError)
      CALL io_read_real1(sk2SetID,(/1/),dimsInt(:1),stars%sk2)
      CALL h5dclose_f(sk2SetID, hdfError)

      dimsInt(:3)=(/2*stars%mx1+1,2*stars%mx2+1,2*stars%mx3+1/)
      CALL h5dopen_f(groupID, 'ig', igSetID, hdfError)
      CALL io_read_integer3(igSetID,(/1,1,1/),dimsInt(:3),stars%ig)
      CALL h5dclose_f(igSetID, hdfError)

      dimsInt(:1)=(/stars%ng3/)
      CALL h5dopen_f(groupID, 'ig2', ig2SetID, hdfError)
      CALL io_read_integer1(ig2SetID,(/1/),dimsInt(:1),stars%ig2)
      CALL h5dclose_f(ig2SetID, hdfError)

      dimsInt(:1)=(/stars%ng3/)
      CALL h5dopen_f(groupID, 'nstr', nstrSetID, hdfError)
      CALL io_read_integer1(nstrSetID,(/1/),dimsInt(:1),stars%nstr)
      CALL h5dclose_f(nstrSetID, hdfError)

      dimsInt(:1)=(/stars%ng2/)
      CALL h5dopen_f(groupID, 'nstr2', nstr2SetID, hdfError)
      CALL io_read_integer1(nstr2SetID,(/1/),dimsInt(:1),stars%nstr2)
      CALL h5dclose_f(nstr2SetID, hdfError)

      dimsInt(:1)=(/stars%ng2/)
      CALL h5dopen_f(groupID, 'phi2', phi2SetID, hdfError)
      CALL io_read_real1(phi2SetID,(/1/),dimsInt(:1),stars%phi2)
      CALL h5dclose_f(phi2SetID, hdfError)

      dimsInt(:4)=(/2,2*stars%mx1+1,2*stars%mx2+1,2*stars%mx3+1/)
      CALL h5dopen_f(groupID, 'rgphs', rgphsSetID, hdfError)
      CALL io_read_complex3(rgphsSetID,(/-1,1,1,1/),dimsInt(:4),stars%rgphs)
      CALL h5dclose_f(rgphsSetID, hdfError)

      dimsInt(:2)=(/stars%kimax+1,2/)
      CALL h5dopen_f(groupID, 'igfft', igfftSetID, hdfError)
      CALL io_read_integer2(igfftSetID,(/1,1/),dimsInt(:2),stars%igfft(0:stars%kimax,:))
      CALL h5dclose_f(igfftSetID, hdfError)

      dimsInt(:2)=(/stars%kimax2+1,2/)
      CALL h5dopen_f(groupID, 'igfft2', igfft2SetID, hdfError)
      CALL io_read_integer2(igfft2SetID,(/1,1/),dimsInt(:2),stars%igfft2(0:stars%kimax2,:))
      CALL h5dclose_f(igfft2SetID, hdfError)

      dimsInt(:2)=(/2,stars%kimax+1/)
      CALL h5dopen_f(groupID, 'pgfft', pgfftSetID, hdfError)
      CALL io_read_complex1(pgfftSetID,(/-1,1/),dimsInt(:2),stars%pgfft(0:stars%kimax))
      CALL h5dclose_f(pgfftSetID, hdfError)

      dimsInt(:2)=(/2,stars%kimax2+1/)
      CALL h5dopen_f(groupID, 'pgfft2', pgfft2SetID, hdfError)
      CALL io_read_complex1(pgfft2SetID,(/-1,1/),dimsInt(:2),stars%pgfft2(0:stars%kimax2))
      CALL h5dclose_f(pgfft2SetID, hdfError)

      dimsInt(:1)=(/ft2_gf_dim/)
      CALL h5dopen_f(groupID, 'ft2_gfx', ft2_gfxSetID, hdfError)
      CALL io_read_real1(ft2_gfxSetID,(/1/),dimsInt(:1),stars%ft2_gfx)
      CALL h5dclose_f(ft2_gfxSetID, hdfError)

      dimsInt(:1)=(/ft2_gf_dim/)
      CALL h5dopen_f(groupID, 'ft2_gfy', ft2_gfySetID, hdfError)
      CALL io_read_real1(ft2_gfySetID,(/1/),dimsInt(:1),stars%ft2_gfy)
      CALL h5dclose_f(ft2_gfySetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE readStarsHDF

   SUBROUTINE writeStepfunctionHDF(fileID, stepfunctionIndex, starsIndex, stars)

      INTEGER(HID_T), INTENT(IN)    :: fileID
      INTEGER,        INTENT(IN)    :: stepfunctionIndex, starsIndex
      TYPE(t_stars),  INTENT(IN)    :: stars

      INTEGER                   :: ifftd

      INTEGER                   :: hdfError
      INTEGER(HID_T)            :: groupID
      INTEGER(HID_T)            :: ustepSpaceID, ustepSetID
      INTEGER(HID_T)            :: ufftSpaceID, ufftSetID
      CHARACTER(LEN=30)         :: groupName
      INTEGER(HSIZE_T)          :: dims(7)
      INTEGER                   :: dimsInt(7)
      LOGICAL                   :: l_exist

      WRITE(groupname,'(a,i0)') '/stepfunction-', stepfunctionIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (l_exist) THEN
         CALL juDFT_error('stepfunction entry '//TRIM(ADJUSTL(groupName))//' already exists.' ,calledby ="writeStepfunctionHDF")
      END IF

      ifftd = size(stars%ufft)

      CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_write_attint0(groupID,'starsIndex',starsIndex)
      CALL io_write_attint0(groupID,'ifftd',ifftd)
      CALL io_write_attint0(groupID,'ng3',stars%ng3)

      dims(:2)=(/2,stars%ng3/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),ustepSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ustep", H5T_NATIVE_DOUBLE, ustepSpaceID, ustepSetID, hdfError)
      CALL h5sclose_f(ustepSpaceID,hdfError)
      CALL io_write_complex1(ustepSetID,(/-1,1/),dimsInt(:2),stars%ustep)
      CALL h5dclose_f(ustepSetID, hdfError)

      dims(:1)=(/ifftd/)
      dimsInt=dims
      CALL h5screate_simple_f(1,dims(:1),ufftSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ufft", H5T_NATIVE_DOUBLE, ufftSpaceID, ufftSetID, hdfError)
      CALL h5sclose_f(ufftSpaceID,hdfError)
      CALL io_write_real1(ufftSetID,(/1/),dimsInt(:1),stars%ufft)
      CALL h5dclose_f(ufftSetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE writeStepfunctionHDF

   SUBROUTINE readStepfunctionHDF(fileID, stepfunctionIndex, stars)

      INTEGER(HID_T), INTENT(IN)    :: fileID
      INTEGER,        INTENT(IN)    :: stepfunctionIndex
      TYPE(t_stars),  INTENT(INOUT) :: stars

      INTEGER                   :: starsIndex, ng3Temp, ifftd, ifftdStars

      INTEGER(HID_T)            :: groupID
      INTEGER                   :: hdfError
      INTEGER                   :: dimsInt(7)
      CHARACTER(LEN=30)         :: groupName
      LOGICAL                   :: l_exist

      INTEGER(HID_T)            :: ustepSetID
      INTEGER(HID_T)            :: ufftSetID

      WRITE(groupname,'(a,i0)') '/stepfunction-', stepfunctionIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (.NOT.l_exist) THEN
         CALL juDFT_error('stepfunction entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readStepfunctionHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      ifftdStars = 27*stars%mx1*stars%mx2*stars%mx3

      CALL io_read_attint0(groupID,'starsIndex',starsIndex)
      CALL io_read_attint0(groupID,'ng3',ng3Temp)
      CALL io_read_attint0(groupID,'ifftd',ifftd)

      IF((ng3Temp.NE.stars%ng3).OR.(ifftd.NE.ifftdStars)) THEN
         WRITE(*,'(a,i7,a,i7)') 'ng3   (stepfunction): ', ng3Temp, ' ng3   (stars) :', stars%ng3
         WRITE(*,'(a,i7,a,i7)') 'ifftd (stepfunction): ', ifftd,   ' ifftd (stars) :', ifftdStars
         CALL juDFT_error('stepfunction entry '//TRIM(ADJUSTL(groupName))//' does not fit to stars.' ,calledby ="readStepfunctionHDF")
      END IF

      IF(ALLOCATED(stars%ustep)) DEALLOCATE(stars%ustep)
      IF(ALLOCATED(stars%ufft)) DEALLOCATE(stars%ufft)

      ALLOCATE(stars%ustep(stars%ng3))
      ALLOCATE(stars%ufft(0:ifftd-1))

      dimsInt(:2)=(/2,ng3Temp/)
      CALL h5dopen_f(groupID, 'ustep', ustepSetID, hdfError)
      CALL io_read_complex1(ustepSetID,(/-1,1/),dimsInt(:2),stars%ustep)
      CALL h5dclose_f(ustepSetID, hdfError)

      dimsInt(:1)=(/ifftd/)
      CALL h5dopen_f(groupID, 'ufft', ufftSetID, hdfError)
      CALL io_read_real1(ufftSetID,(/1/),dimsInt(:1),stars%ufft(0:ifftd-1))
      CALL h5dclose_f(ufftSetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE readStepfunctionHDF

   SUBROUTINE peekStepfunctionHDF(fileID, stepfunctionIndex, starsIndex)

      INTEGER(HID_T), INTENT(IN)    :: fileID
      INTEGER,        INTENT(IN)    :: stepfunctionIndex
      INTEGER,        INTENT(OUT)   :: starsIndex

      INTEGER(HID_T)            :: groupID
      INTEGER                   :: hdfError
      CHARACTER(LEN=30)         :: groupName
      LOGICAL                   :: l_exist

      WRITE(groupname,'(a,i0)') '/stepfunction-', stepfunctionIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (.NOT.l_exist) THEN
         CALL juDFT_error('stepfunction entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readStepfunctionHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attint0(groupID,'starsIndex',starsIndex)
      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE peekStepfunctionHDF

   SUBROUTINE writeLatharmsHDF(fileID, latharmsIndex, latharms)

      INTEGER(HID_T), INTENT(IN)  :: fileID
      INTEGER,        INTENT(IN)  :: latharmsIndex
      TYPE(t_sphhar), INTENT(IN)  :: latharms

      INTEGER                   :: hdfError
      INTEGER(HID_T)            :: groupID
      INTEGER(HID_T)            :: nlhSpaceID, nlhSetID
      INTEGER(HID_T)            :: llhSpaceID, llhSetID
      INTEGER(HID_T)            :: nmemSpaceID, nmemSetID
      INTEGER(HID_T)            :: mlhSpaceID, mlhSetID
      INTEGER(HID_T)            :: clnuSpaceID, clnuSetID
      CHARACTER(LEN=30)         :: groupName
      INTEGER(HSIZE_T)          :: dims(7)
      INTEGER                   :: dimsInt(7)
      LOGICAL                   :: l_exist

      WRITE(groupname,'(a,i0)') '/latharms-', latharmsIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (l_exist) THEN
         CALL juDFT_error('latharms entry '//TRIM(ADJUSTL(groupName))//' already exists.' ,calledby ="writeLatharmsHDF")
      END IF

      CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_write_attint0(groupID,'ntypsd',latharms%ntypsd)
      CALL io_write_attint0(groupID,'memd',latharms%memd)
      CALL io_write_attint0(groupID,'nlhd',latharms%nlhd)

      dims(:1)=(/latharms%ntypsd/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),nlhSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "nlh", H5T_NATIVE_INTEGER, nlhSpaceID, nlhSetID, hdfError)
      CALL h5sclose_f(nlhSpaceID,hdfError)
      CALL io_write_integer1(nlhSetID,(/1/),dimsInt(:1),latharms%nlh)
      CALL h5dclose_f(nlhSetID, hdfError)

      dims(:2)=(/latharms%nlhd+1,latharms%ntypsd/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),llhSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "llh", H5T_NATIVE_INTEGER, llhSpaceID, llhSetID, hdfError)
      CALL h5sclose_f(llhSpaceID,hdfError)
      CALL io_write_integer2(llhSetID,(/1,1/),dimsInt(:2),latharms%llh)
      CALL h5dclose_f(llhSetID, hdfError)

      dims(:2)=(/latharms%nlhd+1,latharms%ntypsd/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),nmemSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "nmem", H5T_NATIVE_INTEGER, nmemSpaceID, nmemSetID, hdfError)
      CALL h5sclose_f(nmemSpaceID,hdfError)
      CALL io_write_integer2(nmemSetID,(/1,1/),dimsInt(:2),latharms%nmem)
      CALL h5dclose_f(nmemSetID, hdfError)

      dims(:3)=(/latharms%memd,latharms%nlhd+1,latharms%ntypsd/)
      dimsInt = dims
      CALL h5screate_simple_f(3,dims(:3),mlhSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "mlh", H5T_NATIVE_INTEGER, mlhSpaceID, mlhSetID, hdfError)
      CALL h5sclose_f(mlhSpaceID,hdfError)
      CALL io_write_integer3(mlhSetID,(/1,1,1/),dimsInt(:3),latharms%mlh)
      CALL h5dclose_f(mlhSetID, hdfError)

      dims(:4)=(/2,latharms%memd,latharms%nlhd+1,latharms%ntypsd/)
      dimsInt = dims
      CALL h5screate_simple_f(4,dims(:4),clnuSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "clnu", H5T_NATIVE_DOUBLE, clnuSpaceID, clnuSetID, hdfError)
      CALL h5sclose_f(clnuSpaceID,hdfError)
      CALL io_write_complex3(clnuSetID,(/-1,1,1,1/),dimsInt(:4),latharms%clnu)
      CALL h5dclose_f(clnuSetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE writeLatharmsHDF

   SUBROUTINE readLatharmsHDF(fileID, latharmsIndex, latharms)

      INTEGER(HID_T), INTENT(IN)  :: fileID
      INTEGER,        INTENT(IN)  :: latharmsIndex
      TYPE(t_sphhar), INTENT(INOUT) :: latharms

      INTEGER(HID_T)            :: nlhSetID, llhSetID, nmemSetID
      INTEGER(HID_T)            :: mlhSetID, clnuSetID, groupID
      CHARACTER(LEN=30)         :: groupName
      INTEGER                   :: hdfError
      INTEGER                   :: dimsInt(7)
      LOGICAL                   :: l_exist

      WRITE(groupname,'(a,i0)') '/latharms-', latharmsIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (.NOT.l_exist) THEN
         CALL juDFT_error('latharms entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readLatharmsHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_read_attint0(groupID,'ntypsd',latharms%ntypsd)
      CALL io_read_attint0(groupID,'memd',latharms%memd)
      CALL io_read_attint0(groupID,'nlhd',latharms%nlhd)

      IF(ALLOCATED(latharms%nlh)) DEALLOCATE(latharms%nlh)
      IF(ALLOCATED(latharms%llh)) DEALLOCATE(latharms%llh)
      IF(ALLOCATED(latharms%nmem)) DEALLOCATE(latharms%nmem)
      IF(ALLOCATED(latharms%mlh)) DEALLOCATE(latharms%mlh)
      IF(ALLOCATED(latharms%clnu)) DEALLOCATE(latharms%clnu)

      ALLOCATE(latharms%clnu(latharms%memd,0:latharms%nlhd,latharms%ntypsd))
      ALLOCATE(latharms%llh(0:latharms%nlhd,latharms%ntypsd))
      ALLOCATE(latharms%mlh(latharms%memd,0:latharms%nlhd,latharms%ntypsd))
      ALLOCATE(latharms%nlh(latharms%ntypsd),latharms%nmem(0:latharms%nlhd,latharms%ntypsd))

      dimsInt(:1)=(/latharms%ntypsd/)
      CALL h5dopen_f(groupID, 'nlh', nlhSetID, hdfError)
      CALL io_read_integer1(nlhSetID,(/1/),dimsInt(:1),latharms%nlh)
      CALL h5dclose_f(nlhSetID, hdfError)

      dimsInt(:2)=(/latharms%nlhd+1,latharms%ntypsd/)
      CALL h5dopen_f(groupID, 'llh', llhSetID, hdfError)
      CALL io_read_integer2(llhSetID,(/1,1/),dimsInt(:2),latharms%llh)
      CALL h5dclose_f(llhSetID, hdfError)

      dimsInt(:2)=(/latharms%nlhd+1,latharms%ntypsd/)
      CALL h5dopen_f(groupID, 'nmem', nmemSetID, hdfError)
      CALL io_read_integer2(nmemSetID,(/1,1/),dimsInt(:2),latharms%nmem)
      CALL h5dclose_f(nmemSetID, hdfError)

      dimsInt(:3)=(/latharms%memd,latharms%nlhd+1,latharms%ntypsd/)
      CALL h5dopen_f(groupID, 'mlh', mlhSetID, hdfError)
      CALL io_read_integer3(mlhSetID,(/1,1,1/),dimsInt(:3),latharms%mlh)
      CALL h5dclose_f(mlhSetID, hdfError)

      dimsInt(:4)=(/2,latharms%memd,latharms%nlhd+1,latharms%ntypsd/)
      CALL h5dopen_f(groupID, 'clnu', clnuSetID, hdfError)
      CALL io_read_complex3(clnuSetID,(/-1,1,1,1/),dimsInt(:4),latharms%clnu)
      CALL h5dclose_f(clnuSetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE readLatharmsHDF

   SUBROUTINE writeStructureHDF(fileID, input, atoms, cell, vacuum, oneD, structureIndex)

      INTEGER(HID_T), INTENT(IN) :: fileID
      INTEGER, INTENT(IN)        :: structureIndex
      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_oneD),INTENT(IN)    :: oneD

      INTEGER(HID_T)            :: groupID
      INTEGER                   :: hdfError
      CHARACTER(LEN=30)         :: groupName
      INTEGER(HSIZE_T)          :: dims(7)
      INTEGER                   :: dimsInt(7)
      LOGICAL                   :: l_exist

      INTEGER(HID_T)                   :: amatSpaceID, amatSetID
      INTEGER(HID_T)                   :: nzSpaceID, nzSetID
      INTEGER(HID_T)                   :: neqSpaceID, neqSetID
      INTEGER(HID_T)                   :: jriSpaceID, jriSetID
      INTEGER(HID_T)                   :: lmaxSpaceID, lmaxSetID
      INTEGER(HID_T)                   :: ngoprSpaceID, ngoprSetID
      INTEGER(HID_T)                   :: ntypsySpaceID, ntypsySetID
      INTEGER(HID_T)                   :: nlhtypSpaceID, nlhtypSetID
      INTEGER(HID_T)                   :: invsatSpaceID, invsatSetID
      INTEGER(HID_T)                   :: rmtSpaceID, rmtSetID
      INTEGER(HID_T)                   :: dxSpaceID, dxSetID
      INTEGER(HID_T)                   :: volmtsSpaceID, volmtsSetID
      INTEGER(HID_T)                   :: rmshSpaceID, rmshSetID
      INTEGER(HID_T)                   :: zatomSpaceID, zatomSetID
      INTEGER(HID_T)                   :: posSpaceID, posSetID
      INTEGER(HID_T)                   :: taualSpaceID, taualSetID

      WRITE(groupname,'(a,i0)') '/structure-', structureIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (l_exist) THEN
         CALL juDFT_error('structure entry '//TRIM(ADJUSTL(groupName))//' already exists.' ,calledby ="writeStructureHDF")
      END IF

      CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_write_attlog0(groupID,'l_film',input%film)

      CALL io_write_attreal0(groupID,'omtil',cell%omtil)
      CALL io_write_attreal0(groupID,'area',cell%area)
      CALL io_write_attreal0(groupID,'z1',cell%z1)
      CALL io_write_attreal0(groupID,'vol',cell%vol)
      CALL io_write_attreal0(groupID,'volint',cell%volint)

      CALL io_write_attint0(groupID,'ntype',atoms%ntype)
      CALL io_write_attint0(groupID,'nat',atoms%nat)
      CALL io_write_attint0(groupID,'lmaxd',atoms%lmaxd)
      CALL io_write_attint0(groupID,'jmtd',atoms%jmtd)

      CALL io_write_attint0(groupID,'nmz',vacuum%nmz)
      CALL io_write_attint0(groupID,'nmzd',vacuum%nmzd)
      CALL io_write_attint0(groupID,'nmzxy',vacuum%nmzxy)
      CALL io_write_attint0(groupID,'nmzxyd',vacuum%nmzxyd)
      CALL io_write_attint0(groupID,'layerd',vacuum%layerd)
      CALL io_write_attint0(groupID,'layers',vacuum%layers)
      CALL io_write_attint0(groupID,'nvac',vacuum%nvac)
      CALL io_write_attint0(groupID,'nvacd',vacuum%nvacd)
      CALL io_write_attint0(groupID,'nstars',vacuum%nstars)
      CALL io_write_attint0(groupID,'nstm',vacuum%nstm)
      CALL io_write_attreal0(groupID,'delz',vacuum%delz)
      CALL io_write_attreal0(groupID,'dvac',vacuum%dvac)

      CALL io_write_attint0(groupID,'od_nq2',oneD%odi%nq2)

      dims(:2)=(/3,3/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),amatSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "amat", H5T_NATIVE_DOUBLE, amatSpaceID, amatSetID, hdfError)
      CALL h5sclose_f(amatSpaceID,hdfError)
      CALL io_write_real2(amatSetID,(/1,1/),dimsInt(:2),cell%amat)
      CALL h5dclose_f(amatSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),nzSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "nz", H5T_NATIVE_INTEGER, nzSpaceID, nzSetID, hdfError)
      CALL h5sclose_f(nzSpaceID,hdfError)
      CALL io_write_integer1(nzSetID,(/1/),dimsInt(:1),atoms%nz)
      CALL h5dclose_f(nzSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),neqSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "neq", H5T_NATIVE_INTEGER, neqSpaceID, neqSetID, hdfError)
      CALL h5sclose_f(neqSpaceID,hdfError)
      CALL io_write_integer1(neqSetID,(/1/),dimsInt(:1),atoms%neq)
      CALL h5dclose_f(neqSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),jriSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "jri", H5T_NATIVE_INTEGER, jriSpaceID, jriSetID, hdfError)
      CALL h5sclose_f(jriSpaceID,hdfError)
      CALL io_write_integer1(jriSetID,(/1/),dimsInt(:1),atoms%jri)
      CALL h5dclose_f(jriSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),lmaxSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "lmax", H5T_NATIVE_INTEGER, lmaxSpaceID, lmaxSetID, hdfError)
      CALL h5sclose_f(lmaxSpaceID,hdfError)
      CALL io_write_integer1(lmaxSetID,(/1/),dimsInt(:1),atoms%lmax)
      CALL h5dclose_f(lmaxSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),ngoprSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ngopr", H5T_NATIVE_INTEGER, ngoprSpaceID, ngoprSetID, hdfError)
      CALL h5sclose_f(ngoprSpaceID,hdfError)
      CALL io_write_integer1(ngoprSetID,(/1/),dimsInt(:1),atoms%ngopr)
      CALL h5dclose_f(ngoprSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),ntypsySpaceID,hdfError)
      CALL h5dcreate_f(groupID, "ntypsy", H5T_NATIVE_INTEGER, ntypsySpaceID, ntypsySetID, hdfError)
      CALL h5sclose_f(ntypsySpaceID,hdfError)
      CALL io_write_integer1(ntypsySetID,(/1/),dimsInt(:1),atoms%ntypsy)
      CALL h5dclose_f(ntypsySetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),nlhtypSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "nlhtyp", H5T_NATIVE_INTEGER, nlhtypSpaceID, nlhtypSetID, hdfError)
      CALL h5sclose_f(nlhtypSpaceID,hdfError)
      CALL io_write_integer1(nlhtypSetID,(/1/),dimsInt(:1),atoms%nlhtyp)
      CALL h5dclose_f(nlhtypSetID, hdfError)

      dims(:1)=(/atoms%nat/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),invsatSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "invsat", H5T_NATIVE_INTEGER, invsatSpaceID, invsatSetID, hdfError)
      CALL h5sclose_f(invsatSpaceID,hdfError)
      CALL io_write_integer1(invsatSetID,(/1/),dimsInt(:1),atoms%invsat)
      CALL h5dclose_f(invsatSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),rmtSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "rmt", H5T_NATIVE_DOUBLE, rmtSpaceID, rmtSetID, hdfError)
      CALL h5sclose_f(rmtSpaceID,hdfError)
      CALL io_write_real1(rmtSetID,(/1/),dimsInt(:1),atoms%rmt)
      CALL h5dclose_f(rmtSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),dxSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "dx", H5T_NATIVE_DOUBLE, dxSpaceID, dxSetID, hdfError)
      CALL h5sclose_f(dxSpaceID,hdfError)
      CALL io_write_real1(dxSetID,(/1/),dimsInt(:1),atoms%dx)
      CALL h5dclose_f(dxSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),volmtsSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "volmts", H5T_NATIVE_DOUBLE, volmtsSpaceID, volmtsSetID, hdfError)
      CALL h5sclose_f(volmtsSpaceID,hdfError)
      CALL io_write_real1(volmtsSetID,(/1/),dimsInt(:1),atoms%volmts)
      CALL h5dclose_f(volmtsSetID, hdfError)

      dims(:2)=(/atoms%jmtd,atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),rmshSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "rmsh", H5T_NATIVE_DOUBLE, rmshSpaceID, rmshSetID, hdfError)
      CALL h5sclose_f(rmshSpaceID,hdfError)
      CALL io_write_real2(rmshSetID,(/1,1/),dimsInt(:2),atoms%rmsh)
      CALL h5dclose_f(rmshSetID, hdfError)

      dims(:1)=(/atoms%ntype/)
      dimsInt = dims
      CALL h5screate_simple_f(1,dims(:1),zatomSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "zatom", H5T_NATIVE_DOUBLE, zatomSpaceID, zatomSetID, hdfError)
      CALL h5sclose_f(zatomSpaceID,hdfError)
      CALL io_write_real1(zatomSetID,(/1/),dimsInt(:1),atoms%zatom)
      CALL h5dclose_f(zatomSetID, hdfError)

      dims(:2)=(/3,atoms%nat/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),posSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "pos", H5T_NATIVE_DOUBLE, posSpaceID, posSetID, hdfError)
      CALL h5sclose_f(posSpaceID,hdfError)
      CALL io_write_real2(posSetID,(/1,1/),dimsInt(:2),atoms%pos)
      CALL h5dclose_f(posSetID, hdfError)

      dims(:2)=(/3,atoms%nat/)
      dimsInt = dims
      CALL h5screate_simple_f(2,dims(:2),taualSpaceID,hdfError)
      CALL h5dcreate_f(groupID, "taual", H5T_NATIVE_DOUBLE, taualSpaceID, taualSetID, hdfError)
      CALL h5sclose_f(taualSpaceID,hdfError)
      CALL io_write_real2(taualSetID,(/1,1/),dimsInt(:2),atoms%taual)
      CALL h5dclose_f(taualSetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE writeStructureHDF

   SUBROUTINE readStructureHDF(fileID, input, atoms, cell, vacuum, oneD, structureIndex)

      INTEGER(HID_T), INTENT(IN)    :: fileID
      INTEGER, INTENT(IN)           :: structureIndex
      TYPE(t_input),INTENT(INOUT)   :: input
      TYPE(t_atoms), INTENT(INOUT)  :: atoms
      TYPE(t_cell), INTENT(INOUT)   :: cell
      TYPE(t_vacuum), INTENT(INOUT) :: vacuum
      TYPE(t_oneD),INTENT(INOUT)    :: oneD

      INTEGER(HID_T)            :: groupID
      INTEGER                   :: hdfError
      CHARACTER(LEN=30)         :: groupName
      INTEGER                   :: dimsInt(7)
      LOGICAL                   :: l_exist

      INTEGER(HID_T)                   :: amatSetID
      INTEGER(HID_T)                   :: nzSetID
      INTEGER(HID_T)                   :: neqSetID
      INTEGER(HID_T)                   :: jriSetID
      INTEGER(HID_T)                   :: lmaxSetID
      INTEGER(HID_T)                   :: ngoprSetID
      INTEGER(HID_T)                   :: ntypsySetID
      INTEGER(HID_T)                   :: nlhtypSetID
      INTEGER(HID_T)                   :: invsatSetID
      INTEGER(HID_T)                   :: rmtSetID
      INTEGER(HID_T)                   :: dxSetID
      INTEGER(HID_T)                   :: volmtsSetID
      INTEGER(HID_T)                   :: rmshSetID
      INTEGER(HID_T)                   :: zatomSetID
      INTEGER(HID_T)                   :: posSetID
      INTEGER(HID_T)                   :: taualSetID

      WRITE(groupname,'(a,i0)') '/structure-', structureIndex

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      IF (.NOT.l_exist) THEN
         CALL juDFT_error('structure entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readStructureHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_read_attlog0(groupID,'l_film',input%film)

      CALL io_read_attreal0(groupID,'omtil',cell%omtil)
      CALL io_read_attreal0(groupID,'area',cell%area)
      CALL io_read_attreal0(groupID,'z1',cell%z1)
      CALL io_read_attreal0(groupID,'vol',cell%vol)
      CALL io_read_attreal0(groupID,'volint',cell%volint)

      CALL io_read_attint0(groupID,'ntype',atoms%ntype)
      CALL io_read_attint0(groupID,'nat',atoms%nat)
      CALL io_read_attint0(groupID,'lmaxd',atoms%lmaxd)
      CALL io_read_attint0(groupID,'jmtd',atoms%jmtd)

      CALL io_read_attint0(groupID,'nmz',vacuum%nmz)
      CALL io_read_attint0(groupID,'nmzd',vacuum%nmzd)
      CALL io_read_attint0(groupID,'nmzxy',vacuum%nmzxy)
      CALL io_read_attint0(groupID,'nmzxyd',vacuum%nmzxyd)
      CALL io_read_attint0(groupID,'layerd',vacuum%layerd)
      CALL io_read_attint0(groupID,'layers',vacuum%layers)
      CALL io_read_attint0(groupID,'nvac',vacuum%nvac)
      CALL io_read_attint0(groupID,'nvacd',vacuum%nvacd)
      CALL io_read_attint0(groupID,'nstars',vacuum%nstars)
      CALL io_read_attint0(groupID,'nstm',vacuum%nstm)
      CALL io_read_attreal0(groupID,'delz',vacuum%delz)
      CALL io_read_attreal0(groupID,'dvac',vacuum%dvac)

      CALL io_read_attint0(groupID,'od_nq2',oneD%odi%nq2)

      IF(ALLOCATED(atoms%nz)) DEALLOCATE(atoms%nz)
      IF(ALLOCATED(atoms%neq)) DEALLOCATE(atoms%neq)
      IF(ALLOCATED(atoms%jri)) DEALLOCATE(atoms%jri)
      IF(ALLOCATED(atoms%lmax)) DEALLOCATE(atoms%lmax)
      IF(ALLOCATED(atoms%ngopr)) DEALLOCATE(atoms%ngopr)
      IF(ALLOCATED(atoms%ntypsy)) DEALLOCATE(atoms%ntypsy)
      IF(ALLOCATED(atoms%nlhtyp)) DEALLOCATE(atoms%nlhtyp)
      IF(ALLOCATED(atoms%invsat)) DEALLOCATE(atoms%invsat)
      IF(ALLOCATED(atoms%rmt)) DEALLOCATE(atoms%rmt)
      IF(ALLOCATED(atoms%dx)) DEALLOCATE(atoms%dx)
      IF(ALLOCATED(atoms%volmts)) DEALLOCATE(atoms%volmts)
      IF(ALLOCATED(atoms%rmsh)) DEALLOCATE(atoms%rmsh)
      IF(ALLOCATED(atoms%zatom)) DEALLOCATE(atoms%zatom)
      IF(ALLOCATED(atoms%pos)) DEALLOCATE(atoms%pos)
      IF(ALLOCATED(atoms%taual)) DEALLOCATE(atoms%taual)

      ALLOCATE(atoms%nz(atoms%ntype))
      ALLOCATE(atoms%neq(atoms%ntype))
      ALLOCATE(atoms%jri(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype))
      ALLOCATE(atoms%ngopr(atoms%nat))
      ALLOCATE(atoms%ntypsy(atoms%nat))
      ALLOCATE(atoms%nlhtyp(atoms%ntype))
      ALLOCATE(atoms%invsat(atoms%nat))
      ALLOCATE(atoms%rmt(atoms%ntype))
      ALLOCATE(atoms%dx(atoms%ntype))
      ALLOCATE(atoms%volmts(atoms%ntype))
      ALLOCATE(atoms%rmsh(atoms%jmtd,atoms%ntype))
      ALLOCATE(atoms%zatom(atoms%ntype))
      ALLOCATE(atoms%pos(3,atoms%nat))
      ALLOCATE(atoms%taual(3,atoms%nat))

      dimsInt(:2)=(/3,3/)
      CALL h5dopen_f(groupID, 'amat', amatSetID, hdfError)
      CALL io_read_real2(amatSetID,(/1,1/),dimsInt(:2),cell%amat)
      CALL h5dclose_f(amatSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'nz', nzSetID, hdfError)
      CALL io_read_integer1(nzSetID,(/1/),dimsInt(:1),atoms%nz)
      CALL h5dclose_f(nzSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'neq', neqSetID, hdfError)
      CALL io_read_integer1(neqSetID,(/1/),dimsInt(:1),atoms%neq)
      CALL h5dclose_f(neqSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'jri', nzSetID, hdfError)
      CALL io_read_integer1(jriSetID,(/1/),dimsInt(:1),atoms%jri)
      CALL h5dclose_f(jriSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'lmax', lmaxSetID, hdfError)
      CALL io_read_integer1(lmaxSetID,(/1/),dimsInt(:1),atoms%lmax)
      CALL h5dclose_f(lmaxSetID, hdfError)

      dimsInt(:1)=(/atoms%nat/)
      CALL h5dopen_f(groupID, 'ngopr', ngoprSetID, hdfError)
      CALL io_read_integer1(ngoprSetID,(/1/),dimsInt(:1),atoms%ngopr)
      CALL h5dclose_f(ngoprSetID, hdfError)

      dimsInt(:1)=(/atoms%nat/)
      CALL h5dopen_f(groupID, 'ntypsy', ntypsySetID, hdfError)
      CALL io_read_integer1(ntypsySetID,(/1/),dimsInt(:1),atoms%ntypsy)
      CALL h5dclose_f(ntypsySetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'nlhtyp', nlhtypSetID, hdfError)
      CALL io_read_integer1(nlhtypSetID,(/1/),dimsInt(:1),atoms%nlhtyp)
      CALL h5dclose_f(nlhtypSetID, hdfError)

      dimsInt(:1)=(/atoms%nat/)
      CALL h5dopen_f(groupID, 'invsat', invsatSetID, hdfError)
      CALL io_read_integer1(invsatSetID,(/1/),dimsInt(:1),atoms%invsat)
      CALL h5dclose_f(invsatSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'rmt', rmtSetID, hdfError)
      CALL io_read_real1(rmtSetID,(/1/),dimsInt(:1),atoms%rmt)
      CALL h5dclose_f(rmtSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'dx', dxSetID, hdfError)
      CALL io_read_real1(dxSetID,(/1/),dimsInt(:1),atoms%dx)
      CALL h5dclose_f(dxSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'volmts', volmtsSetID, hdfError)
      CALL io_read_real1(volmtsSetID,(/1/),dimsInt(:1),atoms%volmts)
      CALL h5dclose_f(volmtsSetID, hdfError)

      dimsInt(:2)=(/atoms%jmtd,atoms%ntype/)
      CALL h5dopen_f(groupID, 'rmsh', rmshSetID, hdfError)
      CALL io_read_real2(rmshSetID,(/1,1/),dimsInt(:2),atoms%rmsh)
      CALL h5dclose_f(rmshSetID, hdfError)

      dimsInt(:1)=(/atoms%ntype/)
      CALL h5dopen_f(groupID, 'zatom', zatomSetID, hdfError)
      CALL io_read_real1(zatomSetID,(/1/),dimsInt(:1),atoms%zatom)
      CALL h5dclose_f(zatomSetID, hdfError)

      dimsInt(:2)=(/3,atoms%nat/)
      CALL h5dopen_f(groupID, 'pos', posSetID, hdfError)
      CALL io_read_real2(posSetID,(/1,1/),dimsInt(:2),atoms%pos)
      CALL h5dclose_f(posSetID, hdfError)

      dimsInt(:2)=(/3,atoms%nat/)
      CALL h5dopen_f(groupID, 'taual', taualSetID, hdfError)
      CALL io_read_real2(taualSetID,(/1,1/),dimsInt(:2),atoms%taual)
      CALL h5dclose_f(taualSetID, hdfError)

      CALL h5gclose_f(groupID, hdfError)

   END SUBROUTINE readStructureHDF

   SUBROUTINE writeDensityHDF(input, fileID, archiveName, densityType, previousDensityIndex,&
                              starsIndex, latharmsIndex, structureIndex, stepfunctionIndex,&
                              date,time,distance,fermiEnergy,l_qfix,iter,&
                              fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

      TYPE(t_input),    INTENT(IN) :: input
      INTEGER(HID_T),   INTENT(IN) :: fileID
      INTEGER,          INTENT(IN) :: densityType, previousDensityIndex
      INTEGER,          INTENT(IN) :: starsIndex, latharmsIndex, structureIndex
      INTEGER,          INTENT(IN) :: stepfunctionIndex
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER, INTENT (IN)         :: date, time, iter
      REAL,    INTENT (IN)         :: fermiEnergy, distance
      LOGICAL, INTENT (IN)         :: l_qfix

      REAL,    INTENT (IN)         :: fr(:,:,:,:)
      REAL,    INTENT (IN)         :: fz(:,:,:)
      COMPLEX, INTENT (IN)         :: fpw(:,:)
      COMPLEX, INTENT (IN)         :: fzxy(:,:,:,:)
      COMPLEX, INTENT (IN)         :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)

      INTEGER                      :: ntype,jmtd,nmzd,nmzxyd,nlhd,ng3,ng2
      INTEGER                      :: nmz, nvac, od_nq2, nmzxy
      INTEGER                      :: hdfError
      LOGICAL                      :: l_film, l_exist
      INTEGER(HID_T)               :: archiveID, groupID
      CHARACTER(LEN=30)            :: groupName, densityTypeName
      INTEGER(HSIZE_T)             :: dims(7)
      INTEGER                      :: dimsInt(7)

      INTEGER(HID_T)               :: frSpaceID, frSetID
      INTEGER(HID_T)               :: fpwSpaceID, fpwSetID
      INTEGER(HID_T)               :: fzSpaceID, fzSetID
      INTEGER(HID_T)               :: fzxySpaceID, fzxySetID
      INTEGER(HID_T)               :: cdomSpaceID, cdomSetID
      INTEGER(HID_T)               :: cdomvzSpaceID, cdomvzSetID
      INTEGER(HID_T)               :: cdomvxySpaceID, cdomvxySetID

      WRITE(groupname,'(a,i0)') '/structure-', structureIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Structure entry "//TRIM(ADJUSTL(groupName))//" does not exist",calledby ="writeDensityHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attlog0(groupID,'l_film',l_film)
      CALL io_read_attint0(groupID,'ntype',ntype)
      CALL io_read_attint0(groupID,'jmtd',jmtd)
      CALL io_read_attint0(groupID,'nmzd',nmzd)
      CALL io_read_attint0(groupID,'nmzxyd',nmzxyd)
      CALL io_read_attint0(groupID,'nmzxy',nmzxy)
      CALL io_read_attint0(groupID,'nmz',nmz)
      CALL io_read_attint0(groupID,'nvac',nvac)
      CALL io_read_attint0(groupID,'od_nq2',od_nq2)
      CALL h5gclose_f(groupID, hdfError)

      WRITE(groupname,'(a,i0)') '/latharms-', latharmsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Latharms entry "//TRIM(ADJUSTL(groupName))//" does not exist",calledby ="writeDensityHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attint0(groupID,'nlhd',nlhd)
      CALL h5gclose_f(groupID, hdfError)

      WRITE(groupname,'(a,i0)') '/stars-', starsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Stars entry "//TRIM(ADJUSTL(groupName))//" does not exist",calledby ="writeDensityHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attint0(groupID,'ng3',ng3)
      CALL io_read_attint0(groupID,'ng2',ng2)
      CALL h5gclose_f(groupID, hdfError)

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))

      SELECT CASE (densityType)
         CASE(DENSITY_TYPE_IN_const)
            densityTypeName = '/in'
         CASE(DENSITY_TYPE_OUT_const)
            densityTypeName = '/out'
         CASE(DENSITY_TYPE_NOCO_IN_const)
            densityTypeName = '/noco_in'
         CASE(DENSITY_TYPE_NOCO_OUT_const)
            densityTypeName = '/noco_out'
         CASE(DENSITY_TYPE_PRECOND_const)
            densityTypeName = '/precond'
         CASE DEFAULT
            CALL juDFT_error("Unknown density type selected",calledby ="writeDensityHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))

      IF(l_exist) THEN
         CALL h5gopen_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)

         CALL io_write_attint0(archiveID,'previousDensityIndex',previousDensityIndex)
         CALL io_write_attint0(archiveID,'starsIndex',starsIndex)
         CALL io_write_attint0(archiveID,'latharmsIndex',latharmsIndex)
         CALL io_write_attint0(archiveID,'structureIndex',structureIndex)
         CALL io_write_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
         CALL io_write_attint0(archiveID,'spins',input%jspins)
         CALL io_write_attint0(archiveID,'iter',iter)
         IF (distance.GE.-1e-10) THEN
            CALL io_write_attreal0(archiveID,'distance',distance)
         END IF

         l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

         IF(l_exist) THEN
            CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

            CALL io_write_attreal0(groupID,'fermiEnergy',fermiEnergy)
            CALL io_write_attlog0(groupID,'l_qfix',l_qfix)

            dimsInt(:4)=(/jmtd,nlhd+1,ntype,input%jspins/)
            CALL h5dopen_f(groupID, 'fr', frSetID, hdfError)
            CALL io_write_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
            CALL h5dclose_f(frSetID, hdfError)

            dimsInt(:3)=(/2,ng3,input%jspins/)
            CALL h5dopen_f(groupID, 'fpw', fpwSetID, hdfError)
            CALL io_write_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
            CALL h5dclose_f(fpwSetID, hdfError)

            IF (l_film) THEN
               dimsInt(:3)=(/nmzd,2,input%jspins/)
               CALL h5dopen_f(groupID, 'fz', fzSetID, hdfError)
               CALL io_write_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
               CALL h5dclose_f(fzSetID, hdfError)

               dimsInt(:5)=(/2,nmzxyd,ng2-1,2,input%jspins/)
               CALL h5dopen_f(groupID, 'fzxy', fzxySetID, hdfError)
               CALL io_write_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
               CALL h5dclose_f(fzxySetID, hdfError)
            END IF

            IF((densityType.EQ.DENSITY_TYPE_NOCO_IN_const).OR.&
               (densityType.EQ.DENSITY_TYPE_NOCO_OUT_const)) THEN

               dimsInt(:2)=(/2,ng3/)
               CALL h5dopen_f(groupID, 'cdom', cdomSetID, hdfError)
               CALL io_write_complex1(cdomSetID,(/-1,1/),dimsInt(:2),cdom)
               CALL h5dclose_f(cdomSetID, hdfError)

               IF (l_film) THEN
                  dimsInt(:3)=(/2,nmz,nvac/)
                  CALL h5dopen_f(groupID, 'cdomvz', cdomvzSetID, hdfError)
                  CALL io_write_complex2(cdomvzSetID,(/-1,1,1/),dimsInt(:3),cdomvz)
                  CALL h5dclose_f(cdomvzSetID, hdfError)

                  dimsInt(:4)=(/2,nmzxy,od_nq2-1,nvac/)
                  CALL h5dopen_f(groupID, 'cdomvxy', cdomvxySetID, hdfError)
                  CALL io_write_complex3(cdomvxySetID,(/-1,1,1,1/),dimsInt(:4),cdomvxy)
                  CALL h5dclose_f(cdomvxySetID, hdfError)
               END IF
            END IF

            CALL h5gclose_f(groupID, hdfError)
         ELSE
            CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

            CALL io_write_attreal0(groupID,'fermiEnergy',fermiEnergy)
            CALL io_write_attlog0(groupID,'l_qfix',l_qfix)

            dims(:4)=(/jmtd,nlhd+1,ntype,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(4,dims(:4),frSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fr", H5T_NATIVE_DOUBLE, frSpaceID, frSetID, hdfError)
            CALL h5sclose_f(frSpaceID,hdfError)
            CALL io_write_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
            CALL h5dclose_f(frSetID, hdfError)

            dims(:3)=(/2,ng3,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(3,dims(:3),fpwSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fpw", H5T_NATIVE_DOUBLE, fpwSpaceID, fpwSetID, hdfError)
            CALL h5sclose_f(fpwSpaceID,hdfError)
            CALL io_write_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
            CALL h5dclose_f(fpwSetID, hdfError)

            IF (l_film) THEN
               dims(:3)=(/nmzd,2,input%jspins/)
               dimsInt = dims
               CALL h5screate_simple_f(3,dims(:3),fzSpaceID,hdfError)
               CALL h5dcreate_f(groupID, "fz", H5T_NATIVE_DOUBLE, fzSpaceID, fzSetID, hdfError)
               CALL h5sclose_f(fzSpaceID,hdfError)
               CALL io_write_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
               CALL h5dclose_f(fzSetID, hdfError)

               dims(:5)=(/2,nmzxyd,ng2-1,2,input%jspins/)
               dimsInt = dims
               CALL h5screate_simple_f(5,dims(:5),fzxySpaceID,hdfError)
               CALL h5dcreate_f(groupID, "fzxy", H5T_NATIVE_DOUBLE, fzxySpaceID, fzxySetID, hdfError)
               CALL h5sclose_f(fzxySpaceID,hdfError)
               CALL io_write_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
               CALL h5dclose_f(fzxySetID, hdfError)
            END IF

            IF((densityType.EQ.DENSITY_TYPE_NOCO_IN_const).OR.&
               (densityType.EQ.DENSITY_TYPE_NOCO_OUT_const)) THEN

               dims(:2)=(/2,ng3/)
               dimsInt = dims
               CALL h5screate_simple_f(2,dims(:2),cdomSpaceID,hdfError)
               CALL h5dcreate_f(groupID, "cdom", H5T_NATIVE_DOUBLE, cdomSpaceID, cdomSetID, hdfError)
               CALL h5sclose_f(cdomSpaceID,hdfError)
               CALL io_write_complex1(cdomSetID,(/-1,1/),dimsInt(:2),cdom)
               CALL h5dclose_f(cdomSetID, hdfError)

               IF (l_film) THEN
                  dims(:3)=(/2,nmz,nvac/)
                  dimsInt = dims
                  CALL h5screate_simple_f(3,dims(:3),cdomvzSpaceID,hdfError)
                  CALL h5dcreate_f(groupID, "cdomvz", H5T_NATIVE_DOUBLE, cdomvzSpaceID, cdomvzSetID, hdfError)
                  CALL h5sclose_f(cdomvzSpaceID,hdfError)
                  CALL io_write_complex2(cdomvzSetID,(/-1,1,1/),dimsInt(:3),cdomvz)
                  CALL h5dclose_f(cdomvzSetID, hdfError)

                  dims(:4)=(/2,nmzxy,od_nq2-1,nvac/)
                  dimsInt = dims
                  CALL h5screate_simple_f(4,dims(:4),cdomvxySpaceID,hdfError)
                  CALL h5dcreate_f(groupID, "cdomvxy", H5T_NATIVE_DOUBLE, cdomvxySpaceID, cdomvxySetID, hdfError)
                  CALL h5sclose_f(cdomvxySpaceID,hdfError)
                  CALL io_write_complex3(cdomvxySetID,(/-1,1,1,1/),dimsInt(:4),cdomvxy)
                  CALL h5dclose_f(cdomvxySetID, hdfError)
               END IF
            END IF

            CALL h5gclose_f(groupID, hdfError)
         END IF

         CALL h5gclose_f(archiveID, hdfError)
      ELSE
         CALL h5gcreate_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)

         CALL io_write_attint0(archiveID,'previousDensityIndex',previousDensityIndex)
         CALL io_write_attint0(archiveID,'starsIndex',starsIndex)
         CALL io_write_attint0(archiveID,'latharmsIndex',latharmsIndex)
         CALL io_write_attint0(archiveID,'structureIndex',structureIndex)
         CALL io_write_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
         CALL io_write_attint0(archiveID,'spins',input%jspins)
         CALL io_write_attint0(archiveID,'iter',iter)
         CALL io_write_attint0(archiveID,'date',date)
         CALL io_write_attint0(archiveID,'time',time)
         CALL io_write_attreal0(archiveID,'distance',distance)

         CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

         CALL io_write_attreal0(groupID,'fermiEnergy',fermiEnergy)
         CALL io_write_attlog0(groupID,'l_qfix',l_qfix)

         dims(:4)=(/jmtd,nlhd+1,ntype,input%jspins/)
         dimsInt = dims
         CALL h5screate_simple_f(4,dims(:4),frSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "fr", H5T_NATIVE_DOUBLE, frSpaceID, frSetID, hdfError)
         CALL h5sclose_f(frSpaceID,hdfError)
         CALL io_write_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
         CALL h5dclose_f(frSetID, hdfError)

         dims(:3)=(/2,ng3,input%jspins/)
         dimsInt = dims
         CALL h5screate_simple_f(3,dims(:3),fpwSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "fpw", H5T_NATIVE_DOUBLE, fpwSpaceID, fpwSetID, hdfError)
         CALL h5sclose_f(fpwSpaceID,hdfError)
         CALL io_write_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
         CALL h5dclose_f(fpwSetID, hdfError)

         IF (l_film) THEN
            dims(:3)=(/nmzd,2,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(3,dims(:3),fzSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fz", H5T_NATIVE_DOUBLE, fzSpaceID, fzSetID, hdfError)
            CALL h5sclose_f(fzSpaceID,hdfError)
            CALL io_write_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
            CALL h5dclose_f(fzSetID, hdfError)

            dims(:5)=(/2,nmzxyd,ng2-1,2,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(5,dims(:5),fzxySpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fzxy", H5T_NATIVE_DOUBLE, fzxySpaceID, fzxySetID, hdfError)
            CALL h5sclose_f(fzxySpaceID,hdfError)
            CALL io_write_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
            CALL h5dclose_f(fzxySetID, hdfError)
         END IF

         IF((densityType.EQ.DENSITY_TYPE_NOCO_IN_const).OR.&
            (densityType.EQ.DENSITY_TYPE_NOCO_OUT_const)) THEN

            dims(:2)=(/2,ng3/)
            dimsInt = dims
            CALL h5screate_simple_f(2,dims(:2),cdomSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "cdom", H5T_NATIVE_DOUBLE, cdomSpaceID, cdomSetID, hdfError)
            CALL h5sclose_f(cdomSpaceID,hdfError)
            CALL io_write_complex1(cdomSetID,(/-1,1/),dimsInt(:2),cdom)
            CALL h5dclose_f(cdomSetID, hdfError)

            IF (l_film) THEN
               dims(:3)=(/2,nmz,nvac/)
               dimsInt = dims
               CALL h5screate_simple_f(3,dims(:3),cdomvzSpaceID,hdfError)
               CALL h5dcreate_f(groupID, "cdomvz", H5T_NATIVE_DOUBLE, cdomvzSpaceID, cdomvzSetID, hdfError)
               CALL h5sclose_f(cdomvzSpaceID,hdfError)
               CALL io_write_complex2(cdomvzSetID,(/-1,1,1/),dimsInt(:3),cdomvz)
               CALL h5dclose_f(cdomvzSetID, hdfError)

               dims(:4)=(/2,nmzxy,od_nq2-1,nvac/)
               dimsInt = dims
               CALL h5screate_simple_f(4,dims(:4),cdomvxySpaceID,hdfError)
               CALL h5dcreate_f(groupID, "cdomvxy", H5T_NATIVE_DOUBLE, cdomvxySpaceID, cdomvxySetID, hdfError)
               CALL h5sclose_f(cdomvxySpaceID,hdfError)
               CALL io_write_complex3(cdomvxySetID,(/-1,1,1,1/),dimsInt(:4),cdomvxy)
               CALL h5dclose_f(cdomvxySetID, hdfError)
            END IF
         END IF

         CALL h5gclose_f(groupID, hdfError)

         CALL h5gclose_f(archiveID, hdfError)
      END IF

   END SUBROUTINE writeDensityHDF

   SUBROUTINE writePotentialHDF(input, fileID, archiveName, potentialType,&
                                starsIndex, latharmsIndex, structureIndex,stepfunctionIndex,&
                                iter,fr,fpw,fz,fzxy)

      TYPE(t_input),    INTENT(IN) :: input
      INTEGER(HID_T),   INTENT(IN) :: fileID
      INTEGER,          INTENT(IN) :: potentialType
      INTEGER,          INTENT(IN) :: starsIndex, latharmsIndex, structureIndex
      INTEGER,          INTENT(IN) :: stepfunctionIndex
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER, INTENT (IN)         :: iter

      REAL,    INTENT (IN)         :: fr(:,:,:,:)
      REAL,    INTENT (IN)         :: fz(:,:,:)
      COMPLEX, INTENT (IN)         :: fpw(:,:)
      COMPLEX, INTENT (IN)         :: fzxy(:,:,:,:)

      INTEGER                      :: ntype,jmtd,nmzd,nmzxyd,nlhd,ng3,ng2
      INTEGER                      :: nmz, nvac, od_nq2, nmzxy
      INTEGER                      :: hdfError
      LOGICAL                      :: l_film, l_exist
      INTEGER(HID_T)               :: archiveID, groupID
      CHARACTER(LEN=30)            :: groupName, potentialTypeName
      INTEGER(HSIZE_T)             :: dims(7)
      INTEGER                      :: dimsInt(7)

      INTEGER(HID_T)               :: frSpaceID, frSetID
      INTEGER(HID_T)               :: fpwSpaceID, fpwSetID
      INTEGER(HID_T)               :: fzSpaceID, fzSetID
      INTEGER(HID_T)               :: fzxySpaceID, fzxySetID

      WRITE(groupname,'(a,i0)') '/structure-', structureIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Structure entry "//TRIM(ADJUSTL(groupName))//" does not exist",calledby ="writePotentialHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attlog0(groupID,'l_film',l_film)
      CALL io_read_attint0(groupID,'ntype',ntype)
      CALL io_read_attint0(groupID,'jmtd',jmtd)
      CALL io_read_attint0(groupID,'nmzd',nmzd)
      CALL io_read_attint0(groupID,'nmzxyd',nmzxyd)
      CALL io_read_attint0(groupID,'nmzxy',nmzxy)
      CALL io_read_attint0(groupID,'nmz',nmz)
      CALL io_read_attint0(groupID,'nvac',nvac)
      CALL io_read_attint0(groupID,'od_nq2',od_nq2)
      CALL h5gclose_f(groupID, hdfError)

      WRITE(groupname,'(a,i0)') '/latharms-', latharmsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Latharms entry "//TRIM(ADJUSTL(groupName))//" does not exist",calledby ="writePotentialHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attint0(groupID,'nlhd',nlhd)
      CALL h5gclose_f(groupID, hdfError)

      WRITE(groupname,'(a,i0)') '/stars-', starsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Stars entry "//TRIM(ADJUSTL(groupName))//" does not exist",calledby ="writePotentialHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)
      CALL io_read_attint0(groupID,'ng3',ng3)
      CALL io_read_attint0(groupID,'ng2',ng2)
      CALL h5gclose_f(groupID, hdfError)

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))

      SELECT CASE (potentialType)
         CASE(POTENTIAL_TYPE_IN_const)
            potentialTypeName = '/in'
         CASE(POTENTIAL_TYPE_OUT_const)
            potentialTypeName = '/out'
         CASE DEFAULT
            CALL juDFT_error("Unknown potential type selected",calledby ="writePotentialHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(potentialTypeName))

      IF(l_exist) THEN
         CALL h5gopen_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)

         CALL io_write_attint0(archiveID,'starsIndex',starsIndex)
         CALL io_write_attint0(archiveID,'latharmsIndex',latharmsIndex)
         CALL io_write_attint0(archiveID,'structureIndex',structureIndex)
         CALL io_write_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
         CALL io_write_attint0(archiveID,'spins',input%jspins)
         CALL io_write_attint0(archiveID,'iter',iter)

         l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

         IF(l_exist) THEN
            CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

            dimsInt(:4)=(/jmtd,nlhd+1,ntype,input%jspins/)
            CALL h5dopen_f(groupID, 'fr', frSetID, hdfError)
            CALL io_write_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
            CALL h5dclose_f(frSetID, hdfError)

            dimsInt(:3)=(/2,ng3,input%jspins/)
            CALL h5dopen_f(groupID, 'fpw', fpwSetID, hdfError)
            CALL io_write_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
            CALL h5dclose_f(fpwSetID, hdfError)

            IF (l_film) THEN
               dimsInt(:3)=(/nmzd,2,input%jspins/)
               CALL h5dopen_f(groupID, 'fz', fzSetID, hdfError)
               CALL io_write_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
               CALL h5dclose_f(fzSetID, hdfError)

               dimsInt(:5)=(/2,nmzxyd,ng2-1,2,input%jspins/)
               CALL h5dopen_f(groupID, 'fzxy', fzxySetID, hdfError)
               CALL io_write_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
               CALL h5dclose_f(fzxySetID, hdfError)
            END IF

            CALL h5gclose_f(groupID, hdfError)
         ELSE
            CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

            dims(:4)=(/jmtd,nlhd+1,ntype,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(4,dims(:4),frSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fr", H5T_NATIVE_DOUBLE, frSpaceID, frSetID, hdfError)
            CALL h5sclose_f(frSpaceID,hdfError)
            CALL io_write_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
            CALL h5dclose_f(frSetID, hdfError)

            dims(:3)=(/2,ng3,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(3,dims(:3),fpwSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fpw", H5T_NATIVE_DOUBLE, fpwSpaceID, fpwSetID, hdfError)
            CALL h5sclose_f(fpwSpaceID,hdfError)
            CALL io_write_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
            CALL h5dclose_f(fpwSetID, hdfError)

            IF (l_film) THEN
               dims(:3)=(/nmzd,2,input%jspins/)
               dimsInt = dims
               CALL h5screate_simple_f(3,dims(:3),fzSpaceID,hdfError)
               CALL h5dcreate_f(groupID, "fz", H5T_NATIVE_DOUBLE, fzSpaceID, fzSetID, hdfError)
               CALL h5sclose_f(fzSpaceID,hdfError)
               CALL io_write_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
               CALL h5dclose_f(fzSetID, hdfError)

               dims(:5)=(/2,nmzxyd,ng2-1,2,input%jspins/)
               dimsInt = dims
               CALL h5screate_simple_f(5,dims(:5),fzxySpaceID,hdfError)
               CALL h5dcreate_f(groupID, "fzxy", H5T_NATIVE_DOUBLE, fzxySpaceID, fzxySetID, hdfError)
               CALL h5sclose_f(fzxySpaceID,hdfError)
               CALL io_write_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
               CALL h5dclose_f(fzxySetID, hdfError)
            END IF

            CALL h5gclose_f(groupID, hdfError)
         END IF

         CALL h5gclose_f(archiveID, hdfError)
      ELSE
         CALL h5gcreate_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)

         CALL io_write_attint0(archiveID,'starsIndex',starsIndex)
         CALL io_write_attint0(archiveID,'latharmsIndex',latharmsIndex)
         CALL io_write_attint0(archiveID,'structureIndex',structureIndex)
         CALL io_write_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
         CALL io_write_attint0(archiveID,'spins',input%jspins)
         CALL io_write_attint0(archiveID,'iter',iter)

         CALL h5gcreate_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

         dims(:4)=(/jmtd,nlhd+1,ntype,input%jspins/)
         dimsInt = dims
         CALL h5screate_simple_f(4,dims(:4),frSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "fr", H5T_NATIVE_DOUBLE, frSpaceID, frSetID, hdfError)
         CALL h5sclose_f(frSpaceID,hdfError)
         CALL io_write_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
         CALL h5dclose_f(frSetID, hdfError)

         dims(:3)=(/2,ng3,input%jspins/)
         dimsInt = dims
         CALL h5screate_simple_f(3,dims(:3),fpwSpaceID,hdfError)
         CALL h5dcreate_f(groupID, "fpw", H5T_NATIVE_DOUBLE, fpwSpaceID, fpwSetID, hdfError)
         CALL h5sclose_f(fpwSpaceID,hdfError)
         CALL io_write_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
         CALL h5dclose_f(fpwSetID, hdfError)

         IF (l_film) THEN
            dims(:3)=(/nmzd,2,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(3,dims(:3),fzSpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fz", H5T_NATIVE_DOUBLE, fzSpaceID, fzSetID, hdfError)
            CALL h5sclose_f(fzSpaceID,hdfError)
            CALL io_write_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
            CALL h5dclose_f(fzSetID, hdfError)

            dims(:5)=(/2,nmzxyd,ng2-1,2,input%jspins/)
            dimsInt = dims
            CALL h5screate_simple_f(5,dims(:5),fzxySpaceID,hdfError)
            CALL h5dcreate_f(groupID, "fzxy", H5T_NATIVE_DOUBLE, fzxySpaceID, fzxySetID, hdfError)
            CALL h5sclose_f(fzxySpaceID,hdfError)
            CALL io_write_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
            CALL h5dclose_f(fzxySetID, hdfError)
         END IF

         CALL h5gclose_f(groupID, hdfError)

         CALL h5gclose_f(archiveID, hdfError)
      END IF

   END SUBROUTINE writePotentialHDF

   SUBROUTINE readDensityHDF(fileID, archiveName, densityType,&
                             fermiEnergy,l_qfix,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

      INTEGER(HID_T), INTENT(IN)   :: fileID
      INTEGER, INTENT(IN)          :: densityType
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER, INTENT (OUT)        :: iter
      REAL,    INTENT (OUT)        :: fermiEnergy
      LOGICAL, INTENT (OUT)        :: l_qfix

      REAL,    INTENT (OUT)        :: fr(:,:,:,:)
      REAL,    INTENT (OUT)        :: fz(:,:,:)
      COMPLEX, INTENT (OUT)        :: fpw(:,:)
      COMPLEX, INTENT (OUT)        :: fzxy(:,:,:,:)
      COMPLEX, INTENT (OUT)        :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)

      INTEGER              :: starsIndex, latharmsIndex, structureIndex, stepfunctionIndex
      INTEGER              :: previousDensityIndex, jspins
      INTEGER              :: ntype,jmtd,nmzd,nmzxyd,nlhd,ng3,ng2
      INTEGER              :: nmz, nvac, od_nq2, nmzxy
      INTEGER              :: localDensityType
      LOGICAL              :: l_film, l_exist
      INTEGER(HID_T)       :: archiveID, groupID, groupBID
      INTEGER              :: hdfError
      CHARACTER(LEN=30)    :: groupName, groupBName, densityTypeName
      INTEGER              :: dimsInt(7)

      INTEGER(HID_T)              :: frSetID
      INTEGER(HID_T)              :: fpwSetID
      INTEGER(HID_T)              :: fzSetID
      INTEGER(HID_T)              :: fzxySetID
      INTEGER(HID_T)              :: cdomSetID
      INTEGER(HID_T)              :: cdomvzSetID
      INTEGER(HID_T)              :: cdomvxySetID

      cdom = CMPLX(0.0,0.0)
      cdomvz = CMPLX(0.0,0.0)
      cdomvxy = CMPLX(0.0,0.0)

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error('density archive '//TRIM(ADJUSTL(archiveName))//' does not exist.' ,calledby ="readDensityHDF")
      END IF

      localDensityType = densityType
      SELECT CASE (densityType)
         CASE(DENSITY_TYPE_IN_const)
            densityTypeName = '/in'
         CASE(DENSITY_TYPE_OUT_const)
            densityTypeName = '/out'
         CASE(DENSITY_TYPE_NOCO_IN_const)
            densityTypeName = '/noco_in'
            groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
            l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
            IF(.NOT.l_exist) THEN
               localDensityType = DENSITY_TYPE_IN_const
               densityTypeName = '/in'
            END IF
         CASE(DENSITY_TYPE_NOCO_OUT_const)
            densityTypeName = '/noco_out'
            groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
            l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
            IF(.NOT.l_exist) THEN
               localDensityType = DENSITY_TYPE_OUT_const
               densityTypeName = '/out'
            END IF
         CASE(DENSITY_TYPE_PRECOND_const)
            densityTypeName = '/precond'
         CASE DEFAULT
            CALL juDFT_error("Unknown density type selected",calledby ="readDensityHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error('Density entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readDensityHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_read_attint0(archiveID,'previousDensityIndex',previousDensityIndex)
      CALL io_read_attint0(archiveID,'starsIndex',starsIndex)
      CALL io_read_attint0(archiveID,'latharmsIndex',latharmsIndex)
      CALL io_read_attint0(archiveID,'structureIndex',structureIndex)
      CALL io_read_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
      CALL io_read_attint0(archiveID,'spins',jspins)
      CALL io_read_attint0(archiveID,'iter',iter)

      WRITE(groupBName,'(a,i0)') '/structure-', structureIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupBName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Structure entry "//TRIM(ADJUSTL(groupBName))//" does not exist",calledby ="readDensityHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupBName)), groupBID, hdfError)
      CALL io_read_attlog0(groupBID,'l_film',l_film)
      CALL io_read_attint0(groupBID,'ntype',ntype)
      CALL io_read_attint0(groupBID,'jmtd',jmtd)
      CALL io_read_attint0(groupBID,'nmzd',nmzd)
      CALL io_read_attint0(groupBID,'nmzxyd',nmzxyd)
      CALL io_read_attint0(groupBID,'nmzxy',nmzxy)
      CALL io_read_attint0(groupBID,'nmz',nmz)
      CALL io_read_attint0(groupBID,'nvac',nvac)
      CALL io_read_attint0(groupBID,'od_nq2',od_nq2)
      CALL h5gclose_f(groupBID, hdfError)

      WRITE(groupBName,'(a,i0)') '/latharms-', latharmsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupBName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Latharms entry "//TRIM(ADJUSTL(groupBName))//" does not exist",calledby ="readDensityHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupBName)), groupBID, hdfError)
      CALL io_read_attint0(groupBID,'nlhd',nlhd)
      CALL h5gclose_f(groupBID, hdfError)

      WRITE(groupBName,'(a,i0)') '/stars-', starsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupBName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Stars entry "//TRIM(ADJUSTL(groupBName))//" does not exist",calledby ="readDensityHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupBName)), groupBID, hdfError)
      CALL io_read_attint0(groupBID,'ng3',ng3)
      CALL io_read_attint0(groupBID,'ng2',ng2)
      CALL h5gclose_f(groupBID, hdfError)

      CALL io_read_attreal0(groupID,'fermiEnergy',fermiEnergy)
      CALL io_read_attlog0(groupID,'l_qfix',l_qfix)

      dimsInt(:4)=(/jmtd,nlhd+1,ntype,jspins/)
      CALL h5dopen_f(groupID, 'fr', frSetID, hdfError)
      CALL io_read_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
      CALL h5dclose_f(frSetID, hdfError)

      dimsInt(:3)=(/2,ng3,jspins/)
      CALL h5dopen_f(groupID, 'fpw', fpwSetID, hdfError)
      CALL io_read_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
      CALL h5dclose_f(fpwSetID, hdfError)

      IF (l_film) THEN
         dimsInt(:3)=(/nmzd,2,jspins/)
         CALL h5dopen_f(groupID, 'fz', fzSetID, hdfError)
         CALL io_read_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
         CALL h5dclose_f(fzSetID, hdfError)

         dimsInt(:5)=(/2,nmzxyd,ng2-1,2,jspins/)
         CALL h5dopen_f(groupID, 'fzxy', fzxySetID, hdfError)
         CALL io_read_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
         CALL h5dclose_f(fzxySetID, hdfError)
      END IF

      IF((localDensityType.EQ.DENSITY_TYPE_NOCO_IN_const).OR.&
         (localDensityType.EQ.DENSITY_TYPE_NOCO_OUT_const)) THEN

         dimsInt(:2)=(/2,ng3/)
         CALL h5dopen_f(groupID, 'cdom', cdomSetID, hdfError)
         CALL io_read_complex1(cdomSetID,(/-1,1/),dimsInt(:2),cdom)
         CALL h5dclose_f(cdomSetID, hdfError)

         IF (l_film) THEN
            dimsInt(:3)=(/2,nmz,nvac/)
            CALL h5dopen_f(groupID, 'cdomvz', cdomvzSetID, hdfError)
            CALL io_read_complex2(cdomvzSetID,(/-1,1,1/),dimsInt(:3),cdomvz)
            CALL h5dclose_f(cdomvzSetID, hdfError)

            dimsInt(:4)=(/2,nmzxy,od_nq2-1,nvac/)
            CALL h5dopen_f(groupID, 'cdomvxy', cdomvxySetID, hdfError)
            CALL io_read_complex3(cdomvxySetID,(/-1,1,1,1/),dimsInt(:4),cdomvxy)
            CALL h5dclose_f(cdomvxySetID, hdfError)
         END IF
      END IF

      CALL h5gclose_f(groupID, hdfError)
      CALL h5gclose_f(archiveID, hdfError)

   END SUBROUTINE readDensityHDF

   SUBROUTINE readPotentialHDF(fileID, archiveName, potentialType,&
                               iter,fr,fpw,fz,fzxy)

      INTEGER(HID_T), INTENT(IN)   :: fileID
      INTEGER, INTENT(IN)          :: potentialType
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER, INTENT (OUT)        :: iter

      REAL,    INTENT (OUT)        :: fr(:,:,:,:)
      REAL,    INTENT (OUT)        :: fz(:,:,:)
      COMPLEX, INTENT (OUT)        :: fpw(:,:)
      COMPLEX, INTENT (OUT)        :: fzxy(:,:,:,:)

      INTEGER              :: starsIndex, latharmsIndex, structureIndex, stepfunctionIndex
      INTEGER              :: jspins
      INTEGER              :: ntype,jmtd,nmzd,nmzxyd,nlhd,ng3,ng2
      INTEGER              :: nmz, nvac, od_nq2, nmzxy
      LOGICAL              :: l_film, l_exist
      INTEGER(HID_T)       :: archiveID, groupID, groupBID
      INTEGER              :: hdfError
      CHARACTER(LEN=30)    :: groupName, groupBName, potentialTypeName
      INTEGER              :: dimsInt(7)

      INTEGER(HID_T)              :: frSetID
      INTEGER(HID_T)              :: fpwSetID
      INTEGER(HID_T)              :: fzSetID
      INTEGER(HID_T)              :: fzxySetID


      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error('density archive '//TRIM(ADJUSTL(archiveName))//' does not exist.' ,calledby ="readPotentialHDF")
      END IF

      SELECT CASE (potentialType)
         CASE(POTENTIAL_TYPE_IN_const)
            potentialTypeName = '/in'
         CASE(POTENTIAL_TYPE_OUT_const)
            potentialTypeName = '/out'
         CASE DEFAULT
            CALL juDFT_error("Unknown potential type selected",calledby ="readPotentialHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(potentialTypeName))
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error('Potential entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="readPotentialHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_read_attint0(archiveID,'starsIndex',starsIndex)
      CALL io_read_attint0(archiveID,'latharmsIndex',latharmsIndex)
      CALL io_read_attint0(archiveID,'structureIndex',structureIndex)
      CALL io_read_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
      CALL io_read_attint0(archiveID,'spins',jspins)
      CALL io_read_attint0(archiveID,'iter',iter)

      WRITE(groupBName,'(a,i0)') '/structure-', structureIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupBName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Structure entry "//TRIM(ADJUSTL(groupBName))//" does not exist",calledby ="readPotentialHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupBName)), groupBID, hdfError)
      CALL io_read_attlog0(groupBID,'l_film',l_film)
      CALL io_read_attint0(groupBID,'ntype',ntype)
      CALL io_read_attint0(groupBID,'jmtd',jmtd)
      CALL io_read_attint0(groupBID,'nmzd',nmzd)
      CALL io_read_attint0(groupBID,'nmzxyd',nmzxyd)
      CALL io_read_attint0(groupBID,'nmzxy',nmzxy)
      CALL io_read_attint0(groupBID,'nmz',nmz)
      CALL io_read_attint0(groupBID,'nvac',nvac)
      CALL io_read_attint0(groupBID,'od_nq2',od_nq2)
      CALL h5gclose_f(groupBID, hdfError)

      WRITE(groupBName,'(a,i0)') '/latharms-', latharmsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupBName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Latharms entry "//TRIM(ADJUSTL(groupBName))//" does not exist",calledby ="readPotentialHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupBName)), groupBID, hdfError)
      CALL io_read_attint0(groupBID,'nlhd',nlhd)
      CALL h5gclose_f(groupBID, hdfError)

      WRITE(groupBName,'(a,i0)') '/stars-', starsIndex
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupBName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("Stars entry "//TRIM(ADJUSTL(groupBName))//" does not exist",calledby ="readPotentialHDF")
      END IF
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupBName)), groupBID, hdfError)
      CALL io_read_attint0(groupBID,'ng3',ng3)
      CALL io_read_attint0(groupBID,'ng2',ng2)
      CALL h5gclose_f(groupBID, hdfError)

      dimsInt(:4)=(/jmtd,nlhd+1,ntype,jspins/)
      CALL h5dopen_f(groupID, 'fr', frSetID, hdfError)
      CALL io_read_real4(frSetID,(/1,1,1,1/),dimsInt(:4),fr)
      CALL h5dclose_f(frSetID, hdfError)

      dimsInt(:3)=(/2,ng3,jspins/)
      CALL h5dopen_f(groupID, 'fpw', fpwSetID, hdfError)
      CALL io_read_complex2(fpwSetID,(/-1,1,1/),dimsInt(:3),fpw)
      CALL h5dclose_f(fpwSetID, hdfError)

      IF (l_film) THEN
         dimsInt(:3)=(/nmzd,2,jspins/)
         CALL h5dopen_f(groupID, 'fz', fzSetID, hdfError)
         CALL io_read_real3(fzSetID,(/1,1,1/),dimsInt(:3),fz)
         CALL h5dclose_f(fzSetID, hdfError)

         dimsInt(:5)=(/2,nmzxyd,ng2-1,2,jspins/)
         CALL h5dopen_f(groupID, 'fzxy', fzxySetID, hdfError)
         CALL io_read_complex4(fzxySetID,(/-1,1,1,1,1/),dimsInt(:5),fzxy)
         CALL h5dclose_f(fzxySetID, hdfError)
      END IF

      CALL h5gclose_f(groupID, hdfError)
      CALL h5gclose_f(archiveID, hdfError)

   END SUBROUTINE readPotentialHDF


   SUBROUTINE peekDensityEntryHDF(fileID, archiveName, densityType,&
                                  iter, starsIndex, latharmsIndex, structureIndex,&
                                  stepfunctionIndex, previousDensityIndex, jspins,&
                                  date, time, distance, fermiEnergy, l_qfix)

      INTEGER(HID_T), INTENT(IN)   :: fileID
      INTEGER, INTENT(IN)          :: densityType
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER, INTENT(OUT)          :: date, time, iter
      INTEGER, INTENT(OUT)          :: starsIndex, latharmsIndex, structureIndex, stepfunctionIndex
      INTEGER, INTENT(OUT)          :: previousDensityIndex, jspins
      REAL,    INTENT(OUT)          :: fermiEnergy, distance
      LOGICAL, INTENT(OUT)          :: l_qfix

      INTEGER              :: localDensityType
      LOGICAL              :: l_exist
      INTEGER(HID_T)       :: archiveID, groupID
      INTEGER              :: hdfError
      CHARACTER(LEN=30)    :: groupName, densityTypeName

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error('density archive '//TRIM(ADJUSTL(archiveName))//' does not exist.' ,calledby ="peekDensityEntryHDF")
      END IF

      localDensityType = densityType
      SELECT CASE (densityType)
         CASE(DENSITY_TYPE_UNDEFINED_const)
            densityTypeName = ''
         CASE(DENSITY_TYPE_IN_const)
            densityTypeName = '/in'
         CASE(DENSITY_TYPE_OUT_const)
            densityTypeName = '/out'
         CASE(DENSITY_TYPE_NOCO_IN_const)
            densityTypeName = '/noco_in'
            groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
            l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
            IF(.NOT.l_exist) THEN
               localDensityType = DENSITY_TYPE_IN_const
               densityTypeName = '/in'
            END IF
         CASE(DENSITY_TYPE_NOCO_OUT_const)
            densityTypeName = '/noco_out'
            groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
            l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
            IF(.NOT.l_exist) THEN
               localDensityType = DENSITY_TYPE_OUT_const
               densityTypeName = '/out'
            END IF
         CASE(DENSITY_TYPE_PRECOND_const)
            densityTypeName = '/precond'
         CASE DEFAULT
            CALL juDFT_error("Unknown density type selected",calledby ="peekDensityEntryHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
      IF(.NOT.l_exist) THEN
         CALL juDFT_error('density entry '//TRIM(ADJUSTL(groupName))//' does not exist.' ,calledby ="peekDensityEntryHDF")
      END IF

      CALL h5gopen_f(fileID, TRIM(ADJUSTL(archiveName)), archiveID, hdfError)
      CALL h5gopen_f(fileID, TRIM(ADJUSTL(groupName)), groupID, hdfError)

      CALL io_read_attint0(archiveID,'previousDensityIndex',previousDensityIndex)
      CALL io_read_attint0(archiveID,'starsIndex',starsIndex)
      CALL io_read_attint0(archiveID,'latharmsIndex',latharmsIndex)
      CALL io_read_attint0(archiveID,'structureIndex',structureIndex)
      CALL io_read_attint0(archiveID,'stepfunctionIndex',stepfunctionIndex)
      CALL io_read_attint0(archiveID,'spins',jspins)
      CALL io_read_attint0(archiveID,'iter',iter)
      CALL io_read_attint0(archiveID,'date',date)
      CALL io_read_attint0(archiveID,'time',time)
      CALL io_read_attreal0(archiveID,'distance',distance)

      IF (densityType.NE.DENSITY_TYPE_UNDEFINED_const) THEN
         CALL io_read_attreal0(groupID,'fermiEnergy',fermiEnergy)
         CALL io_read_attlog0(groupID,'l_qfix',l_qfix)
      END IF

      CALL h5gclose_f(groupID, hdfError)
      CALL h5gclose_f(archiveID, hdfError)

   END SUBROUTINE peekDensityEntryHDF

   SUBROUTINE writeCoreDensityHDF(fileID,input,atoms,dimension,rhcs,tecs,qints)

      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_input),    INTENT(IN) :: input
      TYPE(t_dimension),INTENT(IN) :: DIMENSION

      INTEGER(HID_T), INTENT(IN) :: fileID
      REAL,           INTENT(IN) :: rhcs(atoms%jmtd,atoms%ntype,DIMENSION%jspd)
      REAL,           INTENT(IN) :: tecs(atoms%ntype,DIMENSION%jspd)
      REAL,           INTENT(IN) :: qints(atoms%ntype,DIMENSION%jspd)

      INTEGER hdfError
      INTEGER(HID_T) cdncGroupID, rhcsSpaceID, rhcsSetID
      INTEGER(HID_T) tecsSpaceID, tecsSetID, qintsSpaceID, qintsSetID
      LOGICAL l_exist
      INTEGER(HSIZE_T) :: dims(7)
      INTEGER          :: dimsInt(7)

      l_exist = io_groupexists(fileID,'/cdnc')

      IF(l_exist) THEN ! replace current core density
         CALL h5gopen_f(fileID, '/cdnc', cdncGroupID, hdfError)

         CALL io_write_attint0(cdncGroupID,'jspd',DIMENSION%jspd)
         CALL io_write_attint0(cdncGroupID,'ntype',atoms%ntype)
         CALL io_write_attint0(cdncGroupID,'jmtd',atoms%jmtd)

         dimsInt(:3)=(/atoms%jmtd,atoms%ntype,DIMENSION%jspd/)
         CALL h5dopen_f(cdncGroupID, 'rhcs', rhcsSetID, hdfError)
         CALL io_write_real3(rhcsSetID,(/1,1,1/),dimsInt(:3),rhcs)
         CALL h5dclose_f(rhcsSetID, hdfError)

         dimsInt(:2)=(/atoms%ntype,DIMENSION%jspd/)
         CALL h5dopen_f(cdncGroupID, 'tecs', tecsSetID, hdfError)
         CALL io_write_real2(tecsSetID,(/1,1/),dimsInt(:2),tecs)
         CALL h5dclose_f(tecsSetID, hdfError)

         dimsInt(:2)=(/atoms%ntype,DIMENSION%jspd/)
         CALL h5dopen_f(cdncGroupID, 'qints', qintsSetID, hdfError)
         CALL io_write_real2(qintsSetID,(/1,1/),dimsInt(:2),qints)
         CALL h5dclose_f(qintsSetID, hdfError)

         CALL h5gclose_f(cdncGroupID, hdfError)
      ELSE ! write new core density
         CALL h5gcreate_f(fileID, '/cdnc', cdncGroupID, hdfError)

         CALL io_write_attint0(cdncGroupID,'jspd',DIMENSION%jspd)
         CALL io_write_attint0(cdncGroupID,'ntype',atoms%ntype)
         CALL io_write_attint0(cdncGroupID,'jmtd',atoms%jmtd)

         dims(:3)=(/atoms%jmtd,atoms%ntype,DIMENSION%jspd/)
         dimsInt = dims
         CALL h5screate_simple_f(3,dims(:3),rhcsSpaceID,hdfError)
         CALL h5dcreate_f(cdncGroupID, "rhcs", H5T_NATIVE_DOUBLE, rhcsSpaceID, rhcsSetID, hdfError)
         CALL h5sclose_f(rhcsSpaceID,hdfError)
         CALL io_write_real3(rhcsSetID,(/1,1,1/),dimsInt(:3),rhcs)
         CALL h5dclose_f(rhcsSetID, hdfError)

         dims(:2)=(/atoms%ntype,DIMENSION%jspd/)
         dimsInt = dims
         CALL h5screate_simple_f(2,dims(:2),tecsSpaceID,hdfError)
         CALL h5dcreate_f(cdncGroupID, "tecs", H5T_NATIVE_DOUBLE, tecsSpaceID, tecsSetID, hdfError)
         CALL h5sclose_f(tecsSpaceID,hdfError)
         CALL io_write_real2(tecsSetID,(/1,1/),dimsInt(:2),tecs)
         CALL h5dclose_f(tecsSetID, hdfError)

         dims(:2)=(/atoms%ntype,DIMENSION%jspd/)
         dimsInt = dims
         CALL h5screate_simple_f(2,dims(:2),qintsSpaceID,hdfError)
         CALL h5dcreate_f(cdncGroupID, "qints", H5T_NATIVE_DOUBLE, qintsSpaceID, qintsSetID, hdfError)
         CALL h5sclose_f(qintsSpaceID,hdfError)
         CALL io_write_real2(qintsSetID,(/1,1/),dimsInt(:2),qints)
         CALL h5dclose_f(qintsSetID, hdfError)

         CALL h5gclose_f(cdncGroupID, hdfError)
      END IF

   END SUBROUTINE writeCoreDensityHDF

   SUBROUTINE readCoreDensityHDF(fileID,input,atoms,dimension,rhcs,tecs,qints)

      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_input),    INTENT(IN) :: input
      TYPE(t_dimension),INTENT(IN) :: DIMENSION

      INTEGER(HID_T), INTENT(IN) :: fileID
      REAL,    INTENT(OUT)       :: rhcs(atoms%jmtd,atoms%ntype,DIMENSION%jspd)
      REAL,    INTENT(OUT)       :: tecs(atoms%ntype,DIMENSION%jspd)
      REAL,    INTENT(OUT)       :: qints(atoms%ntype,DIMENSION%jspd)

      INTEGER hdfError
      INTEGER(HID_T) cdncGroupID, rhcsSetID, tecsSetID, qintsSetID
      INTEGER jspdTemp, ntypeTemp, jmtdTemp
      LOGICAL l_exist
      INTEGER :: dimsInt(7)

      l_exist = io_groupexists(fileID,'/cdnc')
      IF(.NOT.l_exist) THEN
         CALL juDFT_error("no core density found",calledby ="readCoreDensityHDF")
      END IF

      CALL h5gopen_f(fileID, '/cdnc', cdncGroupID, hdfError)

      CALL io_read_attint0(cdncGroupID,'jspd',jspdTemp)
      CALL io_read_attint0(cdncGroupID,'ntype',ntypeTemp)
      CALL io_read_attint0(cdncGroupID,'jmtd',jmtdTemp)

      IF(jspdTemp.NE.DIMENSION%jspd) CALL juDFT_error("jspd is inconsistent",calledby ="readCoreDensityHDF")
      IF(ntypeTemp.NE.atoms%ntype) CALL juDFT_error("ntype is inconsistent",calledby ="readCoreDensityHDF")
      IF(jmtdTemp.NE.atoms%jmtd) CALL juDFT_error("jmtd is inconsistent",calledby ="readCoreDensityHDF")

      dimsInt(:3)=(/atoms%jmtd,atoms%ntype,DIMENSION%jspd/)
      CALL h5dopen_f(cdncGroupID, 'rhcs', rhcsSetID, hdfError)
      CALL io_read_real3(rhcsSetID,(/1,1,1/),dimsInt(:3),rhcs)
      CALL h5dclose_f(rhcsSetID, hdfError)

      dimsInt(:2)=(/atoms%ntype,DIMENSION%jspd/)
      CALL h5dopen_f(cdncGroupID, 'tecs', tecsSetID, hdfError)
      CALL io_read_real2(tecsSetID,(/1,1/),dimsInt(:2),tecs)
      CALL h5dclose_f(tecsSetID, hdfError)

      dimsInt(:2)=(/atoms%ntype,DIMENSION%jspd/)
      CALL h5dopen_f(cdncGroupID, 'qints', qintsSetID, hdfError)
      CALL io_read_real2(qintsSetID,(/1,1/),dimsInt(:2),qints)
      CALL h5dclose_f(qintsSetID, hdfError)

      CALL h5gclose_f(cdncGroupID, hdfError)

   END SUBROUTINE readCoreDensityHDF

   LOGICAL FUNCTION deleteDensityEntryHDF(fileID,archiveName)

      INTEGER(HID_T), INTENT(IN)   :: fileID
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER                      :: hdfError
      LOGICAL                      :: l_exist

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))
      IF(.NOT.l_exist) THEN
         deleteDensityEntryHDF = .FALSE.
         RETURN
      END IF

      CALL h5ldelete_f(fileID, archiveName, hdfError)

      deleteDensityEntryHDF = .TRUE.

   END FUNCTION deleteDensityEntryHDF

   LOGICAL FUNCTION isCoreDensityPresentHDF()

      INTEGER(HID_T) :: fileID
      INTEGER        :: currentStarsIndex,currentLatharmsIndex
      INTEGER        :: currentStructureIndex, currentStepfunctionIndex
      INTEGER        :: readDensityIndex, lastDensityIndex
      LOGICAL        :: l_exist

      INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
      IF(.NOT.l_exist) THEN
         isCoreDensityPresentHDF = .FALSE.
         RETURN
      END IF
      CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                       currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
      isCoreDensityPresentHDF = io_groupexists(fileID,'/cdnc')
      CALL closeCDNPOT_HDF(fileID)
   END FUNCTION

   LOGICAL FUNCTION isDensityEntryPresentHDF(fileID,archiveName,densityType)

      INTEGER(HID_T), INTENT(IN)   :: fileID
      INTEGER , INTENT(IN)         :: densityType
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      INTEGER                      :: localDensityType
      CHARACTER(LEN=30)            :: groupName, densityTypeName
      LOGICAL                      :: l_exist

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))
      IF(.NOT.l_exist) THEN
         isDensityEntryPresentHDF = .FALSE.
         RETURN
      END IF

      localDensityType = densityType
      SELECT CASE (densityType)
         CASE(DENSITY_TYPE_UNDEFINED_const)
            densityTypeName = ''
         CASE(DENSITY_TYPE_IN_const)
            densityTypeName = '/in'
         CASE(DENSITY_TYPE_OUT_const)
            densityTypeName = '/out'
         CASE(DENSITY_TYPE_NOCO_IN_const)
            densityTypeName = '/noco_in'
            groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
            l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
            IF(.NOT.l_exist) THEN
               localDensityType = DENSITY_TYPE_IN_const
               densityTypeName = '/in'
            END IF
         CASE(DENSITY_TYPE_NOCO_OUT_const)
            densityTypeName = '/noco_out'
            groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
            l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))
            IF(.NOT.l_exist) THEN
               localDensityType = DENSITY_TYPE_OUT_const
               densityTypeName = '/out'
            END IF
         CASE(DENSITY_TYPE_PRECOND_const)
            densityTypeName = '/precond'
         CASE DEFAULT
            CALL juDFT_error("Unknown density type selected",calledby ="isDensityEntryPresentHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(densityTypeName))
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      isDensityEntryPresentHDF = l_exist
   END FUNCTION

   LOGICAL FUNCTION isPotentialEntryPresentHDF(fileID,archiveName,potentialType)

      INTEGER(HID_T), INTENT(IN)   :: fileID
      INTEGER , INTENT(IN)         :: potentialType
      CHARACTER(LEN=*), INTENT(IN) :: archiveName

      CHARACTER(LEN=30)            :: groupName, potentialTypeName
      LOGICAL                      :: l_exist

      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(archiveName)))
      IF(.NOT.l_exist) THEN
         isPotentialEntryPresentHDF = .FALSE.
         RETURN
      END IF

      SELECT CASE (potentialType)
         CASE(POTENTIAL_TYPE_IN_const)
            potentialTypeName = '/in'
         CASE(POTENTIAL_TYPE_OUT_const)
            potentialTypeName = '/out'
         CASE DEFAULT
            CALL juDFT_error("Unknown potential type selected",calledby ="isPotentialEntryPresentHDF")
      END SELECT

      groupName = TRIM(ADJUSTL(archiveName))//TRIM(ADJUSTL(potentialTypeName))
      l_exist = io_groupexists(fileID,TRIM(ADJUSTL(groupName)))

      isPotentialEntryPresentHDF = l_exist
   END FUNCTION

#endif

END MODULE m_cdnpot_io_hdf
