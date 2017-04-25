!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This module is a wrapper for the charge density I/O
!!!
!!!                             GM'17
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_cdn_io

   USE m_types
   USE m_juDFT
   USE m_loddop
   USE m_wrtdop
   USE m_cdnpot_io_hdf
   USE m_cdnpot_io_common
#ifdef CPP_HDF
   USE hdf5
#endif
   IMPLICIT NONE

   PRIVATE
   PUBLIC printDensityFileInfo
   PUBLIC readDensity, writeDensity
   PUBLIC isDensityFilePresent, isCoreDensityPresent
   PUBLIC readCoreDensity, writeCoreDensity
   PUBLIC readStars, writeStars
   PUBLIC readStepfunction, writeStepfunction
   PUBLIC setStartingDensity, readPrevEFermi, deleteDensities
   PUBLIC storeStructureIfNew
   PUBLIC CDN_INPUT_DEN_const, CDN_OUTPUT_DEN_const
   PUBLIC CDN_ARCHIVE_TYPE_CDN1_const, CDN_ARCHIVE_TYPE_NOCO_const
   PUBLIC CDN_ARCHIVE_TYPE_CDN_const

   INTEGER,          PARAMETER :: CDN_INPUT_DEN_const = 1
   INTEGER,          PARAMETER :: CDN_OUTPUT_DEN_const = 2

   INTEGER,          PARAMETER :: CDN_ARCHIVE_TYPE_CDN1_const = 1
   INTEGER,          PARAMETER :: CDN_ARCHIVE_TYPE_NOCO_const = 2
   INTEGER,          PARAMETER :: CDN_ARCHIVE_TYPE_CDN_const  = 3

   INTEGER,          PARAMETER :: CDN_DIRECT_MODE = 1
   INTEGER,          PARAMETER :: CDN_STREAM_MODE = 2
   INTEGER,          PARAMETER :: CDN_HDF5_MODE   = 3

   CONTAINS

   SUBROUTINE printDensityFileInfo()

      INTEGER            :: mode, i
      LOGICAL            :: l_exist

#ifdef CPP_HDF
      INTEGER(HID_T)    :: fileID
#endif
      INTEGER           :: currentStarsIndex,currentLatharmsIndex
      INTEGER           :: currentStructureIndex,currentStepfunctionIndex
      INTEGER           :: readDensityIndex, lastDensityIndex
      CHARACTER(LEN=30) :: archiveName

      INTEGER           :: dateTemp, timeTemp
      INTEGER           :: iterTemp, starsIndexTemp, latharmsIndexTemp 
      INTEGER           :: structureIndexTemp,stepfunctionIndexTemp
      INTEGER           :: previousDensityIndex, jspinsTemp
      REAL              :: fermiEnergyTemp, distanceTemp
      LOGICAL           :: l_qfixTemp
      CHARACTER(LEN=10) :: dateString
      CHARACTER(LEN=10) :: timeString
      CHARACTER(LEN=19) :: timeStampString
      CHARACTER(LEN=15) :: distanceString

      CALL getMode(mode)

      WRITE(*,*) 'Available densities info:'
      WRITE(*,*)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
            WRITE(*,*) 'densityIndex   iteration   prevDensity   prevDistance        timeStamp'
            DO i = 1, lastDensityIndex
               archiveName = ''
               WRITE(archiveName,'(a,i0)') '/cdn-', i

               l_exist = isDensityEntryPresentHDF(fileID,archiveName,DENSITY_TYPE_UNDEFINED_const)
               IF(.NOT.l_exist) THEN
                  CYCLE
               END IF

               CALL peekDensityEntryHDF(fileID, archiveName, DENSITY_TYPE_UNDEFINED_const,&
                                        iterTemp, starsIndexTemp, latharmsIndexTemp, structureIndexTemp,&
                                        stepfunctionIndexTemp,previousDensityIndex, jspinsTemp,&
                                        dateTemp, timeTemp, distanceTemp, fermiEnergyTemp, l_qfixTemp)

               WRITE(dateString,'(i8)') dateTemp
               WRITE(timeString,'(i6)') timeTemp

               distanceString = ''
               IF (distanceTemp.GE.-1e-10) THEN
                  WRITE(distanceString,'(f15.8)') distanceTemp
               END IF

               WRITE(timeStampString,'(a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2)') &
                  dateString(1:4),'/',dateString(5:6),'/',dateString(7:8),&
                  timeString(1:2),':',timeString(3:4),':',timeString(5:6)

               WRITE(*,'(1x,i7,6x,i7,7x,i7,4x,a15,3x,a)') i, iterTemp, previousDensityIndex, distanceString,&
                                                            TRIM(ADJUSTL(timeStampString))
            END DO
            CALL closeCDNPOT_HDF(fileID)
         ELSE
            WRITE(*,'(a)') "No cdn.hdf file found. Density file info is not available."
         END IF
#else
         WRITE(*,'(a)') "Fleur is not compiled with HDF5 support. Density file info is not available."
#endif
      ELSE
         WRITE(*,'(a)') "Density file info is only available if '-hdf_cdn' switch is used."
      END IF

   END SUBROUTINE printDensityFileInfo


   SUBROUTINE readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                          relCdnIndex,fermiEnergy,l_qfix,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_oneD),INTENT(IN)   :: oneD

      INTEGER, INTENT (IN)      :: inOrOutCDN
      INTEGER, INTENT (IN)      :: relCdnIndex
      INTEGER, INTENT (OUT)     :: iter
      INTEGER, INTENT (IN)      :: archiveType
      REAL,    INTENT (OUT)     :: fermiEnergy
      LOGICAL, INTENT (OUT)     :: l_qfix

      !     ..
      !     .. Array Arguments ..
      COMPLEX, INTENT (OUT) :: fpw(stars%ng3,input%jspins), fzxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
      COMPLEX, INTENT (OUT) :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)
      REAL,    INTENT (OUT) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins), fz(vacuum%nmzd,2,input%jspins)

      ! local variables
      INTEGER            :: mode, datend, k, i, iVac, j, iUnit, l
      LOGICAL            :: l_exist, l_rhomatFile, l_DimChange
      CHARACTER(LEN=30)  :: filename

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER           :: currentStarsIndex,currentLatharmsIndex
      INTEGER           :: currentStructureIndex,currentStepfunctionIndex
      INTEGER           :: readDensityIndex, lastDensityIndex
      INTEGER           :: previousDensityIndex, densityType
      CHARACTER(LEN=30) :: archiveName
      TYPE(t_cell)      :: cellTemp

      fermiEnergy = 0.0
      l_qfix = .FALSE.

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF

         densityType = 0
         archiveName = ''

         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

            IF (archiveType.EQ.CDN_ARCHIVE_TYPE_CDN_const) THEN
               archiveName = 'cdn'
            ELSE
               WRITE(archiveName,'(a,i0)') '/cdn-', readDensityIndex
            END IF

            SELECT CASE (inOrOutCDN)
               CASE (CDN_INPUT_DEN_const)
                  IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
                     densityType = DENSITY_TYPE_NOCO_IN_const
                  ELSE
                     densityType = DENSITY_TYPE_IN_const
                  END IF
               CASE (CDN_OUTPUT_DEN_const)
                  IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
                     densityType = DENSITY_TYPE_NOCO_OUT_const
                  ELSE
                     densityType = DENSITY_TYPE_OUT_const
                  END IF
               CASE DEFAULT
                  WRITE(*,*) 'inOrOutCDN = ', inOrOutCDN
                  CALL juDFT_error("Invalid inOrOutCDN selected.",calledby ="readDensity")
            END SELECT
            l_exist = isDensityEntryPresentHDF(fileID,archiveName,densityType)
            CALL closeCDNPOT_HDF(fileID)
         END IF

         IF (l_exist) THEN
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

            CALL readDensityHDF(fileID, input, stars, sphhar, atoms, vacuum, oneD, archiveName, densityType,&
                                fermiEnergy,l_qfix,l_DimChange,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

            CALL closeCDNPOT_HDF(fileID)

            IF(l_DimChange) THEN
               CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                           1,-1.0,fermiEnergy,l_qfix,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)
            END IF
         ELSE
            WRITE(*,*) 'cdn.hdf file or relevant density entry not found.'
            WRITE(*,*) 'Falling back to stream access file cdn.str.'
            mode = CDN_STREAM_MODE
         END IF
#endif
      END IF

      IF(mode.EQ.CDN_STREAM_MODE) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from cdn.str and exit subroutine

         ELSE
            WRITE(*,*) 'cdn.str file not found.'
            WRITE(*,*) 'Falling back to direct access file cdn1.'
            mode = CDN_DIRECT_MODE
         END IF
      END IF

      IF (mode.EQ.CDN_DIRECT_MODE) THEN
         filename = 'cdn1'
         l_rhomatFile = .FALSE.
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
            INQUIRE(file="rhomat_inp",EXIST=l_exist)
            IF (l_exist) filename = 'rhomat_inp'
            IF (inOrOutCDN.EQ.CDN_OUTPUT_DEN_const) THEN
               INQUIRE(file="rhomat_out",EXIST=l_exist)
               IF (l_exist) filename = 'rhomat_out'
            END IF
            IF(l_exist) l_rhomatFile = .TRUE.
         END IF
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_CDN_const) THEN
            filename = 'cdn'
         END IF

         INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=l_exist)
         IF(.NOT.l_exist) THEN
            CALL juDFT_error("charge density file "//TRIM(ADJUSTL(filename))//" missing",calledby ="readDensity")
         END IF

         iUnit = 93
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),FORM='unformatted',STATUS='old')

         IF ((inOrOutCDN.EQ.CDN_OUTPUT_DEN_const).AND.(archiveType.NE.CDN_ARCHIVE_TYPE_NOCO_const)) THEN
            ! call loddop to move the file position to the output density
            CALL loddop(stars,vacuum,atoms,sphhar,input,sym,&
                        iUnit,iter,fr,fpw,fz,fzxy)
         END IF

         ! read in the density
         CALL loddop(stars,vacuum,atoms,sphhar,input,sym,&
                     iUnit,iter,fr,fpw,fz,fzxy)

         ! read in additional data if l_noco and data is present
         IF ((archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const).AND.l_rhomatFile) THEN
            READ (iUnit,iostat=datend) (cdom(k),k=1,stars%ng3)
            IF (datend == 0) THEN
               IF (input%film) THEN
                  READ (iUnit) ((cdomvz(i,iVac),i=1,vacuum%nmz),iVac=1,vacuum%nvac)
                  READ (iUnit) (((cdomvxy(i,j-1,iVac),i=1,vacuum%nmzxy),j=2,oneD%odi%nq2), iVac=1,vacuum%nvac)
               END IF
            ELSE
               ! (datend < 0)  =>  no off-diagonal magnetisation stored
               !                   in "rhomat_inp"
               IF (datend > 0) THEN
                  WRITE(*,*) 'datend = ', datend
                  CALL juDFT_error("density file has illegal state.",calledby ="readDensity")
               END IF
               cdom = CMPLX(0.0,0.0)
               IF (input%film) THEN
                  cdomvz = CMPLX(0.0,0.0)
                  cdomvxy = CMPLX(0.0,0.0)
               END IF
            END IF
         ELSE IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
            cdom = CMPLX(0.0,0.0)
            IF (input%film) THEN
               cdomvz = CMPLX(0.0,0.0)
               cdomvxy = CMPLX(0.0,0.0)
            END IF
         END IF
         CLOSE(iUnit)
      END IF

   END SUBROUTINE readDensity

   SUBROUTINE writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                           relCdnIndex,distance,fermiEnergy,l_qfix,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_oneD),INTENT(IN)   :: oneD

      INTEGER, INTENT (IN)      :: inOrOutCDN
      INTEGER, INTENT (IN)      :: relCdnIndex, iter
      INTEGER, INTENT (IN)      :: archiveType
      REAL,    INTENT (IN)      :: fermiEnergy, distance
      LOGICAL, INTENT (IN)      :: l_qfix
      !     ..
      !     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: fpw(stars%ng3,input%jspins), fzxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
      COMPLEX, INTENT (IN) :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)
      REAL,    INTENT (IN) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins), fz(vacuum%nmzd,2,input%jspins)

      TYPE(t_stars)        :: starsTemp
      TYPE(t_vacuum)       :: vacuumTemp
      TYPE(t_atoms)        :: atomsTemp
      TYPE(t_sphhar)       :: sphharTemp
      TYPE(t_input)        :: inputTemp
      TYPE(t_sym)          :: symTemp
      TYPE(t_cell)         :: cellTemp
      TYPE(t_oneD)         :: oneDTemp

      COMPLEX, ALLOCATABLE :: fpwTemp(:,:), fzxyTemp(:,:,:,:)
      REAL, ALLOCATABLE    :: frTemp(:,:,:,:), fzTemp(:,:,:)

      INTEGER           :: mode, iterTemp, k, i, iVac, j, iUnit
      INTEGER           :: d1, d10, asciioffset, iUnitTemp
      LOGICAL           :: l_exist, l_storeIndices, l_writeNew, l_same
      LOGICAL           :: l_writeAll
      CHARACTER(len=30) :: filename
      CHARACTER(len=5)  :: cdnfile

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER           :: currentStarsIndex,currentLatharmsIndex
      INTEGER           :: currentStructureIndex,currentStepfunctionIndex
      INTEGER           :: readDensityIndex, writeDensityIndex, lastDensityIndex
      INTEGER           :: previousDensityIndex, densityType
      INTEGER           :: starsIndexTemp, latharmsIndexTemp, structureIndexTemp
      INTEGER           :: stepfunctionIndexTemp
      INTEGER           :: jspinsTemp
      INTEGER           :: date, time, dateTemp, timeTemp
      REAL              :: fermiEnergyTemp, distanceTemp
      LOGICAL           :: l_qfixTemp
      CHARACTER(LEN=30) :: archiveName
      CHARACTER(LEN=8)  :: dateString
      CHARACTER(LEN=10) :: timeString
      CHARACTER(LEN=10) :: zone

      CALL getMode(mode)
      CALL DATE_AND_TIME(dateString,timeString,zone)
      READ(dateString,'(i8)') date
      READ(timeString,'(i6)') time

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         CALL checkAndWriteMetadataHDF(fileID, input, atoms, cell, vacuum, oneD, stars, sphhar, sym,&
                                       currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                       currentStepfunctionIndex,l_storeIndices)

         previousDensityIndex = readDensityIndex
         writeDensityIndex = readDensityIndex
         IF(relCdnIndex.NE.0) THEN
            writeDensityIndex = lastDensityIndex+relCdnIndex
            lastDensityIndex = writeDensityIndex
            readDensityIndex = writeDensityIndex
            l_storeIndices = .TRUE.
         END IF

         archiveName = ''
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_CDN_const) THEN
            archiveName = 'cdn'
         ELSE
            WRITE(archiveName,'(a,i0)') '/cdn-', writeDensityIndex
         END IF

         densityType = 0
         SELECT CASE (inOrOutCDN)
            CASE (CDN_INPUT_DEN_const)
               IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
                  densityType = DENSITY_TYPE_NOCO_IN_const
               ELSE
                  densityType = DENSITY_TYPE_IN_const
               END IF
            CASE (CDN_OUTPUT_DEN_const)
               IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
                  densityType = DENSITY_TYPE_NOCO_OUT_const
               ELSE
                  densityType = DENSITY_TYPE_OUT_const
               END IF
            CASE DEFAULT
               WRITE(*,*) 'inOrOutCDN = ', inOrOutCDN
               CALL juDFT_error("Invalid inOrOutCDN selected.",calledby ="writeDensity")
         END SELECT

         IF(relCdnIndex.EQ.0) THEN
            l_exist = isDensityEntryPresentHDF(fileID,archiveName,DENSITY_TYPE_UNDEFINED_const)
            IF(l_exist) THEN
               CALL peekDensityEntryHDF(fileID, archiveName, DENSITY_TYPE_UNDEFINED_const,&
                                        iterTemp, starsIndexTemp, latharmsIndexTemp, structureIndexTemp,&
                                        stepfunctionIndexTemp,previousDensityIndex, jspinsTemp,&
                                        dateTemp, timeTemp, distanceTemp, fermiEnergyTemp, l_qfixTemp)
            END IF
         END IF

         ALLOCATE (fzxyTemp(vacuum%nmzxyd,stars%ng2-1,2,input%jspins))
         ALLOCATE (fzTemp(vacuum%nmzd,2,input%jspins))
         fzTemp(:,:,:) = fz(:,:,:)
         fzxyTemp(:,:,:,:) = fzxy(:,:,:,:)
         IF(vacuum%nvac.EQ.1) THEN
            fzTemp(:,2,:)=fzTemp(:,1,:)
            IF (sym%invs) THEN
               fzxyTemp(:,:,2,:) = CONJG(fzxyTemp(:,:,1,:))
            ELSE
               fzxyTemp(:,:,2,:) = fzxyTemp(:,:,1,:)
            END IF
         END IF

         CALL writeDensityHDF(input, fileID, archiveName, densityType, previousDensityIndex,&
                              currentStarsIndex, currentLatharmsIndex, currentStructureIndex,&
                              currentStepfunctionIndex,date,time,distance,fermiEnergy,l_qfix,iter+relCdnIndex,&
                              fr,fpw,fzTemp,fzxyTemp,cdom,cdomvz,cdomvxy)

         DEALLOCATE(fzTemp,fzxyTemp)

         IF(l_storeIndices) THEN
            CALL writeCDNHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                    currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
         END IF

         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         ! Write density to cdn.str file
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE
         filename = 'cdn1'
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
            filename = 'rhomat_inp'
            IF(inOrOutCDN.EQ.CDN_OUTPUT_DEN_const) filename = 'rhomat_out'
         END IF
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_CDN_const) THEN
            filename = 'cdn'
         END IF

         IF ((relCdnIndex.EQ.1).AND.(archiveType.EQ.CDN_ARCHIVE_TYPE_CDN1_const).AND.(iter.EQ.0)) THEN
            INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=l_exist)
            IF(l_exist) THEN
               CALL juDFT_error("Trying to generate starting density while a density exists.",calledby ="writeDensity")
            END IF
         END IF

         iUnit = 93
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),FORM='unformatted',STATUS='unknown')

         IF ((relCdnIndex.EQ.1).AND.(archiveType.EQ.CDN_ARCHIVE_TYPE_CDN1_const).AND.(iter.GE.1)) THEN
            inputTemp%jspins = input%jspins
            vacuumTemp%nmzxyd = vacuum%nmzxyd
            atomsTemp%jmtd = atoms%jmtd
            sphharTemp%nlhd = sphhar%nlhd
            vacuumTemp%nmzd = vacuum%nmzd
            atomsTemp%ntype = atoms%ntype
            ALLOCATE (sphharTemp%nlh(SIZE(sphhar%nlh)))
            sphharTemp%nlh(:) = sphhar%nlh(:)
            ALLOCATE (atomsTemp%ntypsy(SIZE(atoms%ntypsy)))
            atomsTemp%ntypsy(:) = atoms%ntypsy(:)
            ALLOCATE (atomsTemp%jri(SIZE(atoms%jri)))
            atomsTemp%jri(:) = atoms%jri(:)
            ALLOCATE (atomsTemp%neq(SIZE(atoms%neq)))
            atomsTemp%neq(:) = atoms%neq(:)
            ALLOCATE (atomsTemp%zatom(SIZE(atoms%zatom)))
            atomsTemp%zatom(:) = atoms%zatom(:)
            ALLOCATE (atomsTemp%rmt(SIZE(atoms%rmt)))
            atomsTemp%rmt(:) = atoms%rmt(:)
            ALLOCATE (atomsTemp%dx(SIZE(atoms%dx)))
            atomsTemp%dx(:) = atoms%dx(:)
            starsTemp%ng3 = stars%ng3
            symTemp%invs = sym%invs
            inputTemp%film = input%film
            vacuumTemp%nvac = vacuum%nvac
            starsTemp%ng2 = stars%ng2
            symTemp%invs2 = sym%invs2
            ALLOCATE (fpwTemp(stars%ng3,input%jspins))
            ALLOCATE (fzxyTemp(vacuum%nmzxyd,stars%ng2-1,2,input%jspins))
            ALLOCATE (frTemp(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins))
            ALLOCATE (fzTemp(vacuum%nmzd,2,input%jspins))

            !--->    generate name of file to hold the results of this iteration
            d1 = MOD(iter,10)
            d10 = MOD(INT((iter+0.5)/10),10)
            asciioffset = IACHAR('1')-1
            IF ( d10.GE.10 ) asciioffset = IACHAR('7')
            cdnfile = 'cdn'//ACHAR(d10+asciioffset)//ACHAR(d1+IACHAR('1')-1)

            iUnitTemp = 72
            OPEN (iUnitTemp,file=cdnfile,form='unformatted',status='unknown')
            REWIND iUnitTemp


            CALL loddop(starsTemp,vacuumTemp,atomsTemp,sphharTemp,inputTemp,symTemp,&
                        iUnit,iterTemp,frTemp,fpwTemp,fzTemp,fzxyTemp)
            CALL wrtdop(starsTemp,vacuumTemp,atomsTemp,sphharTemp,inputTemp,symTemp,&
                        iUnitTemp,iterTemp,frTemp,fpwTemp,fzTemp,fzxyTemp)
            CALL loddop(starsTemp,vacuumTemp,atomsTemp,sphharTemp,inputTemp,symTemp,&
                        iUnit,iterTemp,frTemp,fpwTemp,fzTemp,fzxyTemp)
            CALL wrtdop(starsTemp,vacuumTemp,atomsTemp,sphharTemp,inputTemp,symTemp,&
                        iUnitTemp,iterTemp,frTemp,fpwTemp,fzTemp,fzxyTemp)

            CLOSE(iUnitTemp)
            REWIND iUnit

            DEALLOCATE (fzTemp, frTemp, fzxyTemp, fpwTemp)
            DEALLOCATE (atomsTemp%neq, atomsTemp%jri, atomsTemp%zatom, atomsTemp%ntypsy, sphharTemp%nlh)
            DEALLOCATE (atomsTemp%rmt, atomsTemp%dx)
         END IF

         IF ((inOrOutCDN.EQ.CDN_OUTPUT_DEN_const).AND.(archiveType.NE.CDN_ARCHIVE_TYPE_NOCO_const)) THEN

            ! Generate data in temp arrays and variables to be able to perform loddop call.
            ! loddop is called to move the file position to the output density position.
            inputTemp%jspins = input%jspins
            vacuumTemp%nmzxyd = vacuum%nmzxyd
            atomsTemp%jmtd = atoms%jmtd
            sphharTemp%nlhd = sphhar%nlhd
            vacuumTemp%nmzd = vacuum%nmzd
            atomsTemp%ntype = atoms%ntype
            ALLOCATE (sphharTemp%nlh(SIZE(sphhar%nlh)))
            sphharTemp%nlh(:) = sphhar%nlh(:)
            ALLOCATE (atomsTemp%ntypsy(SIZE(atoms%ntypsy)))
            atomsTemp%ntypsy(:) = atoms%ntypsy(:)
            ALLOCATE (atomsTemp%jri(SIZE(atoms%jri)))
            atomsTemp%jri(:) = atoms%jri(:)
            ALLOCATE (atomsTemp%neq(SIZE(atoms%neq)))
            atomsTemp%neq(:) = atoms%neq(:)
            starsTemp%ng3 = stars%ng3
            symTemp%invs = sym%invs
            inputTemp%film = input%film
            vacuumTemp%nvac = vacuum%nvac
            starsTemp%ng2 = stars%ng2
            symTemp%invs2 = sym%invs2
            ALLOCATE (fpwTemp(stars%ng3,input%jspins))
            ALLOCATE (fzxyTemp(vacuum%nmzxyd,stars%ng2-1,2,input%jspins))
            ALLOCATE (frTemp(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins))
            ALLOCATE (fzTemp(vacuum%nmzd,2,input%jspins))

            CALL loddop(starsTemp,vacuumTemp,atomsTemp,sphharTemp,inputTemp,symTemp,&
                        iUnit,iterTemp,frTemp,fpwTemp,fzTemp,fzxyTemp)

            DEALLOCATE (fzTemp, frTemp, fzxyTemp, fpwTemp)
            DEALLOCATE (atomsTemp%neq, atomsTemp%jri, atomsTemp%ntypsy, sphharTemp%nlh)
         END IF

         ! Write the density
         CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym,&
                     iUnit,iter+relCdnIndex,fr,fpw,fz,fzxy)

         ! Write additional data if l_noco
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
            WRITE (iUnit) (cdom(k),k=1,stars%ng3)
            IF (input%film) THEN
               WRITE (iUnit) ((cdomvz(i,iVac),i=1,vacuum%nmz),iVac=1,vacuum%nvac)
               WRITE (iUnit) (((cdomvxy(i,j-1,iVac),i=1,vacuum%nmzxy),j=2,oneD%odi%nq2), iVac=1,vacuum%nvac)
            END IF
         END IF

         CLOSE(iUnit)
      END IF

   END SUBROUTINE writeDensity

   SUBROUTINE readPrevEFermi(eFermiPrev,l_error)

      REAL,    INTENT(OUT) :: eFermiPrev
      LOGICAL, INTENT(OUT) :: l_error

      INTEGER        :: mode
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER        :: currentStarsIndex,currentLatharmsIndex
      INTEGER        :: currentStructureIndex,currentStepfunctionIndex
      INTEGER        :: readDensityIndex, lastDensityIndex

      INTEGER           :: starsIndex, latharmsIndex, structureIndex
      INTEGER           :: stepfunctionIndex
      INTEGER           :: date, time, iter, jspins, previousDensityIndex
      REAL              :: fermiEnergy, distance
      LOGICAL           :: l_qfix, l_exist
      CHARACTER(LEN=30) :: archiveName

      CALL getMode(mode)

      eFermiPrev = 0.0
      l_error = .FALSE.
      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
         WRITE(archiveName,'(a,i0)') '/cdn-', readDensityIndex
         CALL peekDensityEntryHDF(fileID, archiveName, DENSITY_TYPE_UNDEFINED_const,&
                                  iter, starsIndex, latharmsIndex, structureIndex, stepfunctionIndex,&
                                  previousDensityIndex, jspins, date, time, distance, fermiEnergy, l_qfix)
         archiveName = ''
         WRITE(archiveName,'(a,i0)') '/cdn-', previousDensityIndex
         l_exist = isDensityEntryPresentHDF(fileID,archiveName,DENSITY_TYPE_NOCO_OUT_const)
         IF(l_exist) THEN
            CALL peekDensityEntryHDF(fileID, archiveName, DENSITY_TYPE_NOCO_OUT_const,&
                                     iter, starsIndex, latharmsIndex, structureIndex, stepfunctionIndex,&
                                     previousDensityIndex, jspins, date, time, distance, fermiEnergy, l_qfix)
            eFermiPrev = fermiEnergy
         ELSE
            l_error = .TRUE.
         END IF
         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         STOP 'cdn.str not yet implemented!'
      ELSE
         l_error = .TRUE.
      END IF

   END SUBROUTINE

   SUBROUTINE readCoreDensity(input,atoms,dimension,rhcs,tecs,qints)

      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_dimension),INTENT(IN) :: DIMENSION

      REAL, INTENT(OUT) :: rhcs(atoms%jmtd,atoms%ntype,DIMENSION%jspd)
      REAL, INTENT(OUT) :: tecs(atoms%ntype,DIMENSION%jspd)
      REAL, INTENT(OUT) :: qints(atoms%ntype,DIMENSION%jspd)

      INTEGER            :: mode, iUnit, iSpin, iAtom, i
      LOGICAL            :: l_exist
      CHARACTER(LEN=30)  :: filename

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER        :: currentStarsIndex,currentLatharmsIndex
      INTEGER        :: currentStructureIndex,currentStepfunctionIndex
      INTEGER        :: readDensityIndex, lastDensityIndex

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         l_exist = isCoreDensityPresentHDF()
         IF (l_exist) THEN
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
            CALL readCoreDensityHDF(fileID,input,atoms,dimension,rhcs,tecs,qints)
            CALL closeCDNPOT_HDF(fileID)
            RETURN
         ELSE
            WRITE(*,*) 'No core density is available in HDF5 format.'
            WRITE(*,*) 'Falling back to stream access file cdn.str.'
            mode = CDN_STREAM_MODE
         END IF
#endif
      END IF

      IF(mode.EQ.CDN_STREAM_MODE) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from cdn.str and exit subroutine

            RETURN
         ELSE
            WRITE(*,*) 'cdn.str file not found.'
            WRITE(*,*) 'Falling back to direct access file cdnc.'
            mode = CDN_DIRECT_MODE
         END IF
      END IF

      IF (mode.EQ.CDN_DIRECT_MODE) THEN
         filename = 'cdnc'
         INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=l_exist)
         IF(.NOT.l_exist) THEN
            CALL juDFT_error("core charge density file "//TRIM(ADJUSTL(filename))//" missing",calledby ="readCoreDensity")
         END IF
         iUnit = 17
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),form='unformatted',status='unknown')
         DO iSpin = 1,input%jspins
            DO iAtom = 1,atoms%ntype
               READ (iUnit) (rhcs(i,iAtom,iSpin),i=1,atoms%jri(iAtom))
               READ (iUnit) tecs(iAtom,iSpin)
            END DO
            READ (iUnit) (qints(iAtom,iSpin),iAtom=1,atoms%ntype)
         END DO
         CLOSE (iUnit)
      END IF

   END SUBROUTINE readCoreDensity

   SUBROUTINE writeCoreDensity(input,atoms,dimension,rhcs,tecs,qints)

      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_dimension),INTENT(IN) :: DIMENSION

      REAL, INTENT(IN) :: rhcs(atoms%jmtd,atoms%ntype,DIMENSION%jspd)
      REAL, INTENT(IN) :: tecs(atoms%ntype,DIMENSION%jspd)
      REAL, INTENT(IN) :: qints(atoms%ntype,DIMENSION%jspd)

      INTEGER :: mode, iUnit, iSpin, iAtom, i

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER        :: currentStarsIndex,currentLatharmsIndex
      INTEGER        :: currentStructureIndex,currentStepfunctionIndex
      INTEGER        :: readDensityIndex, lastDensityIndex

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
         CALL writeCoreDensityHDF(fileID,input,atoms,dimension,rhcs,tecs,qints)
         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         ! Write core density to cdn.str file
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE
         iUnit = 17
         OPEN (iUnit,file='cdnc',form='unformatted',status='unknown')
         DO iSpin = 1,input%jspins
            DO iAtom = 1,atoms%ntype
               WRITE (iUnit) (rhcs(i,iAtom,iSpin),i=1,atoms%jri(iAtom))
               WRITE (iUnit) tecs(iAtom,iSpin)
            END DO
            WRITE (iUnit) (qints(iAtom,iSpin),iAtom=1,atoms%ntype)
         END DO
         CLOSE (iUnit)
      END IF

   END SUBROUTINE writeCoreDensity

   SUBROUTINE storeStructureIfNew(input, atoms, cell, vacuum, oneD, sym)

      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_oneD),INTENT(IN)    :: oneD
      TYPE(t_sym),INTENT(IN)     :: sym

      TYPE(t_input)              :: inputTemp
      TYPE(t_atoms)              :: atomsTemp
      TYPE(t_cell)               :: cellTemp
      TYPE(t_vacuum)             :: vacuumTemp
      TYPE(t_oneD)               :: oneDTemp
      TYPE(t_sym)                :: symTemp

      INTEGER        :: mode
      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
      LOGICAL        :: l_writeStructure, l_same
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         l_writeStructure = .FALSE.
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         IF (currentStructureIndex.EQ.0) THEN
            currentStructureIndex = currentStructureIndex + 1
            l_writeStructure = .TRUE.
         ELSE
            CALL readStructureHDF(fileID, inputTemp, atomsTemp, cellTemp, vacuumTemp, oneDTemp, symTemp, currentStructureIndex)
            CALL compareStructure(atoms, vacuum, cell, sym, atomsTemp, vacuumTemp, cellTemp, symTemp, l_same)
            IF(.NOT.l_same) THEN
               currentStructureIndex = currentStructureIndex + 1
               l_writeStructure = .TRUE.
            END IF
         END IF

         IF (l_writeStructure) THEN
            CALL writeStructureHDF(fileID, input, atoms, cell, vacuum, oneD, sym, currentStructureIndex)
            CALL writeCDNHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                    currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
         END IF

         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         ! Write stars to stars file
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE
         ! In direct access mode no structure information is written to any file.
      END IF

   END SUBROUTINE storeStructureIfNew

   SUBROUTINE writeStars(stars,l_xcExtended,l_ExtData)

      TYPE(t_stars),INTENT(IN)   :: stars
      LOGICAL, INTENT(IN)        :: l_xcExtended, l_ExtData

      INTEGER        :: mode, ngz, izmin, izmax

      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      INTEGER :: igz(stars%ng3)

      ngz = 0
      izmin = 0
      izmax = 0
      igz = 0

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         currentStarsIndex = currentStarsIndex + 1
         CALL writeStarsHDF(fileID, currentStarsIndex, currentStructureIndex, stars)

         CALL writeCDNHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                 currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         ! Write stars to stars file
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE
         OPEN (51,file='stars',form='unformatted',status='unknown')
         WRITE (51) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax,stars%mx1,stars%mx2,stars%mx3
         IF(l_ExtData) THEN
            IF (.NOT.l_xcExtended) THEN
               WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,&
                          stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                          stars%igfft2,stars%pgfft2
            ELSE
               WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,&
                          stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                          stars%igfft2,stars%pgfft2,stars%ft2_gfx,stars%ft2_gfy
            END IF
         ELSE
            IF (.NOT.l_xcExtended) THEN
               WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,&
                          stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                          stars%igfft2,stars%pgfft2
            ELSE
               WRITE (51) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,&
                          stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                          stars%igfft2,stars%pgfft2,stars%ft2_gfx,stars%ft2_gfy
            END IF
         END IF
         CLOSE (51)
      END IF
   END SUBROUTINE writeStars

   SUBROUTINE readStars(stars,l_xcExtended,l_ExtData,l_error)

      TYPE(t_stars),INTENT(INOUT) :: stars
      LOGICAL, INTENT(IN)         :: l_xcExtended,l_ExtData
      LOGICAL, INTENT(OUT)        :: l_error


      TYPE(t_stars)               :: starsTemp
      INTEGER                     :: mode, ioStatus, ngz,izmin,izmax
      LOGICAL                     :: l_exist, l_same

      INTEGER        :: structureIndexTemp
      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      INTEGER :: igz(stars%ng3)

      l_error = .FALSE.
      ngz = 0
      izmin = 0
      izmax = 0
      igz = 0

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
#ifdef CPP_HDF
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
            IF (currentStarsIndex.LT.1) THEN
               mode = CDN_DIRECT_MODE ! (no stars entry found in cdn.hdf file)
            ELSE
               CALL peekStarsHDF(fileID, currentStarsIndex, structureIndexTemp)
               l_same = structureIndexTemp.EQ.currentStructureIndex
               IF(l_same) THEN
                  CALL readStarsHDF(fileID, currentStarsIndex, starsTemp)
                  CALL compareStars(stars, starsTemp, l_same)
               END IF
               IF(l_same) THEN
                  CALL readStarsHDF(fileID, currentStarsIndex, stars)
               ELSE
                 mode = CDN_DIRECT_MODE ! (no adequate stars entry found in cdn.hdf file)
               END IF
            END IF
            CALL closeCDNPOT_HDF(fileID)
#endif
         END IF
         IF(.NOT.l_exist) THEN
            mode = CDN_STREAM_MODE
         END IF
      END IF

      IF(mode.EQ.CDN_STREAM_MODE) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF (l_exist) THEN
            STOP 'cdn.str code path not yet implemented!'
         END IF
         IF (.NOT.l_exist) THEN
            mode = CDN_DIRECT_MODE
         END IF
      END IF

      IF (mode.EQ.CDN_DIRECT_MODE) THEN
         INQUIRE(FILE='stars',EXIST=l_exist)
         IF(.NOT.l_exist) THEN
            l_error = .TRUE.
            RETURN
         END IF
         OPEN (51,file='stars',form='unformatted',status='unknown')

         READ (51,IOSTAT=ioStatus) stars%gmax,stars%ng3,stars%ng2,ngz,izmin,izmax,stars%mx1,stars%mx2,stars%mx3
         IF (ioStatus.NE.0) THEN
            l_error = .TRUE.
            RETURN
         END IF

         IF (l_ExtData) THEN
            IF (.NOT.l_xcExtended) THEN
               READ (51,IOSTAT=ioStatus) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,&
                                         stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                                         stars%igfft2,stars%pgfft2
               stars%ft2_gfx = 0.0
               stars%ft2_gfy = 0.0
            ELSE
               READ (51,IOSTAT=ioStatus) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%phi2,stars%kv3,stars%kv2,&
                                         stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                                         stars%igfft2,stars%pgfft2,stars%ft2_gfx,stars%ft2_gfy
            END IF
         ELSE
            IF (.NOT.l_xcExtended) THEN
               READ (51,IOSTAT=ioStatus) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,&
                                         stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                                         stars%igfft2,stars%pgfft2
               stars%ft2_gfx = 0.0
               stars%ft2_gfy = 0.0
            ELSE
               READ (51,IOSTAT=ioStatus) stars%nstr,stars%nstr2,stars%rgphs,stars%sk3,stars%sk2,stars%kv3,stars%kv2,&
                                         stars%ig,stars%ig2,igz,stars%kimax,stars%igfft,stars%pgfft,stars%kimax2,&
                                         stars%igfft2,stars%pgfft2,stars%ft2_gfx,stars%ft2_gfy
            END IF
         END IF

         IF (ioStatus.NE.0) THEN
            l_error = .TRUE.
            RETURN
         END IF

         CLOSE (51)
      END IF

   END SUBROUTINE readStars

   SUBROUTINE writeStepfunction(stars)

      TYPE(t_stars),INTENT(IN) :: stars

      INTEGER                  :: mode, ifftd, i

      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      ifftd=size(stars%ufft)

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         currentStepfunctionIndex = currentStepfunctionIndex + 1
         CALL writeStepfunctionHDF(fileID, currentStepfunctionIndex, currentStarsIndex, currentStructureIndex, stars)
         CALL writeCDNHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                 currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         ! Write stars to stars file
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE
         OPEN (14,file='wkf2',form='unformatted',status='unknown')

         WRITE (14) stars%ng3,ifftd
         WRITE (14) (stars%ustep(i),i=1,stars%ng3)
         WRITE (14) (stars%ufft(i),i=0,ifftd-1)

         CLOSE (14)
      END IF

   END SUBROUTINE writeStepfunction

   SUBROUTINE readStepfunction(stars, atoms, cell, vacuum, l_error)

      TYPE(t_stars),INTENT(INOUT)   :: stars
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_cell), INTENT(IN)      :: cell
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      LOGICAL, INTENT(OUT)          :: l_error

      TYPE(t_stars)                 :: starsTemp
      TYPE(t_input)                 :: inputTemp
      TYPE(t_atoms)                 :: atomsTemp
      TYPE(t_cell)                  :: cellTemp
      TYPE(t_vacuum)                :: vacuumTemp
      TYPE(t_oneD)                  :: oneDTemp

      INTEGER        :: mode
      INTEGER        :: ifftd, ng3Temp, ifftdTemp, ioStatus, i, starsIndexTemp
      INTEGER        :: structureIndexTemp
      LOGICAL        :: l_exist, l_same, l_sameTemp

      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      l_error = .FALSE.
      ioStatus = 0
      ifftd = 27*stars%mx1*stars%mx2*stars%mx3

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
#ifdef CPP_HDF
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
            IF(currentStepfunctionIndex.GT.0) THEN
               CALL peekStepfunctionHDF(fileID, currentStepfunctionIndex, starsIndexTemp, structureIndexTemp)
               IF((starsIndexTemp.EQ.currentStarsIndex).AND.(structureIndexTemp.EQ.currentStructureIndex)) THEN
                  CALL readStepfunctionHDF(fileID, currentStepfunctionIndex, stars)
               ELSE
                  mode = CDN_STREAM_MODE ! No adequate stepfunction entry found. Fall back to other IO modes.
               END IF
            ELSE
               mode = CDN_STREAM_MODE ! No adequate stepfunction entry found. Fall back to other IO modes.
            END IF
            CALL closeCDNPOT_HDF(fileID)
#endif
         END IF
         IF(.NOT.l_exist) THEN
            mode = CDN_STREAM_MODE
         END IF
      END IF

      IF(mode.EQ.CDN_STREAM_MODE) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF (l_exist) THEN
            STOP 'cdn.str code path not yet implemented!'
         END IF
         IF (.NOT.l_exist) THEN
            mode = CDN_DIRECT_MODE
         END IF
      END IF

      IF (mode.EQ.CDN_DIRECT_MODE) THEN
         INQUIRE(FILE='wkf2',EXIST=l_exist)
         IF(.NOT.l_exist) THEN
            l_error = .TRUE.
            RETURN
         END IF
         OPEN (14,file='wkf2',form='unformatted',status='unknown')
         ng3temp=0;ifftdTemp=0
         READ (14,IOSTAT=ioStatus) ng3Temp, ifftdTemp
         IF (ng3Temp.NE.stars%ng3) ioStatus = 1
         IF (ifftdTemp.NE.ifftd) ioStatus = 1
         IF (ioStatus.NE.0) THEN
            l_error = .TRUE.
            CLOSE (14)
            RETURN
         END IF
         READ (14) (stars%ustep(i),i=1,stars%ng3)
         READ (14) (stars%ufft(i),i=0,ifftd-1)

         CLOSE (14)
      END IF

   END SUBROUTINE readStepfunction

   SUBROUTINE setStartingDensity(l_noco)

      LOGICAL,INTENT(IN) :: l_noco

#ifdef CPP_HDF
      INTEGER(HID_T)    :: fileID
#endif
      INTEGER           :: currentStarsIndex,currentLatharmsIndex
      INTEGER           :: currentStructureIndex,currentStepfunctionIndex
      INTEGER           :: readDensityIndex, lastDensityIndex
      INTEGER           :: sdIndex, ioStatus, mode
      INTEGER           :: densityType
      CHARACTER(LEN=20) :: numberString
      CHARACTER(LEN=30) :: archiveName
      LOGICAL           :: l_exist

      IF (.NOT.juDFT_was_argument("-sd")) THEN
         RETURN
      END IF
      numberString = juDFT_string_for_argument("-sd")
      IF (TRIM(ADJUSTL(numberString)).EQ.'') THEN
         CALL juDFT_error("No number for starting density set in command line.",calledby ="setStartingDensity")
      END IF
      ioStatus = 0
      READ(numberString,'(i8)',iostat=ioStatus) sdIndex
      IF(ioStatus.NE.0) THEN
         CALL juDFT_error("Could not convert starting density index string to number.",calledby ="setStartingDensity")
      END IF

      WRITE(*,'(a,i0,a)') 'Using density ', sdIndex, ' as starting density.'

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
         densityType = DENSITY_TYPE_IN_const
         IF(l_noco) THEN
            densityType = DENSITY_TYPE_NOCO_IN_const
         END IF
         archiveName = ''
         WRITE(archiveName,'(a,i0)') '/cdn-', sdIndex
         l_exist = isDensityEntryPresentHDF(fileID,archiveName,densityType)
         IF(.NOT.l_exist) THEN
            WRITE(*,*) 'archiveName: ', TRIM(ADJUSTL(archiveName))
            CALL juDFT_error("For selected starting density index no in-density is present.",calledby ="setStartingDensity")
         END IF
         CALL writeCDNHeaderData(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                 currentStepfunctionIndex,sdIndex,lastDensityIndex)
         CALL closeCDNPOT_HDF(fileID)
#endif
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE      
         WRITE(*,*) 'Explicit setting of starting density in direct access mode'
         WRITE(*,*) 'not implemented.'
         WRITE(*,*) ''
         WRITE(*,*) 'Ignoring -sd command line argument.'
      END IF

   END SUBROUTINE setStartingDensity

   SUBROUTINE deleteDensities()

#ifdef CPP_HDF
      INTEGER(HID_T)    :: fileID
#endif
      INTEGER           :: currentStarsIndex,currentLatharmsIndex
      INTEGER           :: currentStructureIndex,currentStepfunctionIndex
      INTEGER           :: readDensityIndex, lastDensityIndex
      INTEGER           :: ioStatus, mode, i
      INTEGER           :: startNumber, endNumber, separatorIndex
      CHARACTER(LEN=20) :: ddString
      CHARACTER(LEN=30) :: archiveName
      LOGICAL           :: l_exist, l_deleted

      IF (.NOT.juDFT_was_argument("-delden")) THEN
         RETURN
      END IF
      ddString = juDFT_string_for_argument("-delden")
      IF (TRIM(ADJUSTL(ddString)).EQ.'') THEN
         CALL juDFT_error("Densities to be deleted not specified.",calledby ="deleteDensities")
      END IF

      separatorIndex = -1
      startNumber = -1
      endNumber = -1
      DO i = 1, LEN(TRIM(ADJUSTL(ddString)))
         IF(VERIFY(ddString(i:i),'1234567890').NE.0) THEN
            IF ((ddString(i:i).EQ.'-').AND.(separatorIndex.EQ.-1)) THEN
               separatorIndex = i
            ELSE
               CALL juDFT_error("density deletion string format error",calledby ="deleteDensities")
            END IF
         END IF
      END DO

      IF(separatorIndex.NE.-1) THEN
         READ(ddString(1:separatorIndex-1),'(i0)') startNumber
         READ(ddString(separatorIndex+1:LEN(TRIM(ADJUSTL(ddString)))),'(i0)') endNumber
      ELSE
         READ(ddString(1:LEN(TRIM(ADJUSTL(ddString)))),'(i0)') startNumber
         READ(ddString(1:LEN(TRIM(ADJUSTL(ddString)))),'(i0)') endNumber
      END IF

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
#ifdef CPP_HDF
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

            DO i = startNumber, endNumber
               archiveName = ''
               WRITE(archiveName,'(a,i0)') '/cdn-', i

               l_exist = isDensityEntryPresentHDF(fileID,archiveName,DENSITY_TYPE_UNDEFINED_const)
               IF(.NOT.l_exist) THEN
                  CYCLE
               END IF
               
               l_deleted = deleteDensityEntryHDF(fileID,archiveName)
               IF (l_deleted) THEN
                  WRITE(*,*) 'deleted density entry ', TRIM(ADJUSTL(archiveName))
               END IF
            END DO

            CALL closeCDNPOT_HDF(fileID)
#endif
            WRITE(*,*) 'Please note:'
            WRITE(*,*) 'The deletion of the densities does not free the associated disk space.'
            WRITE(*,*) 'To do this you have to repack the cdn.hdf file.'
            WRITE(*,*) 'It can be done by using the tool h5repack, e.g., by invoking'
            WRITE(*,*) 'h5repack -i cdn.hdf -o cdn-packed.hdf'
            WRITE(*,*) 'mv cdn-packed.hdf cdn.hdf'
         ELSE
            WRITE(*,*) "No cdn.hdf file found. No density entry deleted."
         END IF
      ELSE IF(mode.EQ.CDN_STREAM_MODE) THEN
         STOP 'CDN_STREAM_MODE not yet implemented!'
      ELSE      
         WRITE(*,*) 'Explicit deletion of densities in direct access mode'
         WRITE(*,*) 'not implemented.'
         WRITE(*,*) ''
         WRITE(*,*) 'Ignoring -dd command line argument.'
      END IF

      CALL juDFT_error("Densities deleted.")
      
   END SUBROUTINE deleteDensities

   SUBROUTINE getMode(mode)
      INTEGER, INTENT(OUT) :: mode

      mode = CDN_DIRECT_MODE
      IF (juDFT_was_argument("-stream_cdn")) THEN
         mode=CDN_STREAM_MODE
      END IF
      IF (juDFT_was_argument("-hdf_cdn")) THEN
#ifdef CPP_HDF
         mode=CDN_HDF5_MODE
#else
         WRITE(*,*) 'Code not compiled with HDF5 support.'
         WRITE(*,*) 'Falling back to direct access.'
#endif
      END IF
   END SUBROUTINE getMode

   LOGICAL FUNCTION isDensityFilePresent(archiveType)

      INTEGER, INTENT(IN) :: archiveType

      LOGICAL             :: l_exist
      INTEGER             :: mode

      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      CALL getMode(mode)

      IF (mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF(l_exist) THEN
#ifdef CPP_HDF
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
            CALL closeCDNPOT_HDF(fileID)
            IF(readDensityIndex.GT.0) THEN
               isDensityFilePresent = .TRUE.
               RETURN
            END IF
#endif
         END IF
      END IF

      IF ((mode.EQ.CDN_STREAM_MODE).OR.(mode.EQ.CDN_HDF5_MODE)) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF(l_exist) THEN
            isDensityFilePresent = l_exist
            RETURN
         END IF
      END IF

      !cdn1 or rhomat_inp should be enough for any mode...
      INQUIRE(FILE='cdn1',EXIST=l_exist)
      IF (archiveType.EQ.CDN_ARCHIVE_TYPE_CDN1_const) THEN
         isDensityFilePresent = l_exist
         RETURN
      END IF
      IF (archiveType.NE.CDN_ARCHIVE_TYPE_NOCO_const) THEN
         CALL juDFT_error("Illegal archive type selected.",calledby ="isDensityFilePresent")
      END IF
      IF (l_exist) THEN
         isDensityFilePresent = l_exist
         RETURN
      END IF
      INQUIRE(FILE='rhomat_inp',EXIST=l_exist)
      isDensityFilePresent = l_exist
   END FUNCTION isDensityFilePresent

   LOGICAL FUNCTION isCoreDensityPresent()

      LOGICAL             :: l_exist
      INTEGER             :: mode

      CALL getMode(mode)

      IF (mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         l_exist = isCoreDensityPresentHDF()
         IF(l_exist) THEN
            isCoreDensityPresent = l_exist
            RETURN
         END IF
#endif
      END IF

      IF ((mode.EQ.CDN_STREAM_MODE).OR.(mode.EQ.CDN_HDF5_MODE)) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF(l_exist) THEN
            STOP 'Not yet implemented!'
            RETURN
         END IF
      END IF

      !cdnc should be enough for any mode...
      INQUIRE(FILE='cdnc',EXIST=l_exist)
      IF (l_exist) THEN
         isCoreDensityPresent = l_exist
         RETURN
      END IF
      isCoreDensityPresent = .FALSE.
   END FUNCTION isCoreDensityPresent


END MODULE m_cdn_io
