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
   USE m_constants
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
   PUBLIC storeStructureIfNew,transform_by_moving_atoms
   PUBLIC getIOMode
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

      CALL getIOMode(mode)

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
                          relCdnIndex,fermiEnergy,l_qfix,den,inFilename)

      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_cell), INTENT(IN)     :: cell
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_potden),INTENT(INOUT) :: den

      INTEGER, INTENT (IN)      :: inOrOutCDN
      INTEGER, INTENT (IN)      :: relCdnIndex
      INTEGER, INTENT (IN)      :: archiveType
      REAL,    INTENT (OUT)     :: fermiEnergy
      LOGICAL, INTENT (OUT)     :: l_qfix

      CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: inFilename

      ! local variables
      INTEGER            :: mode, datend, k, i, iVac, j, iUnit, l, numLines, ioStatus, iofl
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

      COMPLEX, ALLOCATABLE :: cdomvz(:,:)

      fermiEnergy = 0.0
      l_qfix = .FALSE.

      CALL getIOMode(mode)

#ifndef CPP_HDF
      filename = 'cdn.hdf'
      IF(PRESENT(inFilename)) filename = TRIM(ADJUSTL(inFilename))//'.hdf'
      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
      IF (l_exist) THEN
         CALL juDFT_warn('Fleur not compiled for HDF5, but '//TRIM(ADJUSTL(filename))//' present',calledby='readDensity')
      END IF
#endif

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF

         densityType = 0
         archiveName = ''

         filename = 'cdn.hdf'
         IF(PRESENT(inFilename)) filename = TRIM(ADJUSTL(inFilename))//'.hdf'

         INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
         IF (l_exist) THEN
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex,inFilename)

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
                             currentStepfunctionIndex,readDensityIndex,lastDensityIndex,inFilename)

            CALL readDensityHDF(fileID, input, stars, sphhar, atoms, vacuum, oneD, archiveName, densityType,&
                                fermiEnergy,l_qfix,l_DimChange,den)

            CALL closeCDNPOT_HDF(fileID)

            IF(l_DimChange) THEN
               CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                           1,-1.0,fermiEnergy,l_qfix,den)
            END IF
         ELSE
            INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=l_exist)
            IF(l_exist) THEN
               WRITE(*,*) 'densityType is ', densityType
               WRITE(*,*) 'Relevant density entry is '//TRIM(ADJUSTL(archiveName))//'.'
               WRITE(*,*) 'Entry not found in '//TRIM(ADJUSTL(filename))//'.'
            ELSE
               WRITE(*,*) TRIM(ADJUSTL(filename))//' file not found.'
            END IF
            WRITE(*,*) 'Falling back to stream access.'
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

         IF(PRESENT(inFilename)) filename = inFilename

         INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=l_exist)
         IF(.NOT.l_exist) THEN
            CALL juDFT_error("charge density file "//TRIM(ADJUSTL(filename))//" missing",calledby ="readDensity")
         END IF

         iUnit = 93
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),FORM='unformatted',STATUS='old')

         IF ((inOrOutCDN.EQ.CDN_OUTPUT_DEN_const).AND.(archiveType.NE.CDN_ARCHIVE_TYPE_NOCO_const)) THEN
            ! call loddop to move the file position to the output density
            CALL loddop(stars,vacuum,atoms,sphhar,input,sym,&
                        iUnit,den%iter,den%mt,den%pw,den%vacz,den%vacxy)
         END IF

         ! read in the density
         CALL loddop(stars,vacuum,atoms,sphhar,input,sym,&
                     iUnit,den%iter,den%mt,den%pw,den%vacz,den%vacxy)

         ! read in additional data if l_noco and data is present
         IF ((archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const).AND.l_rhomatFile) THEN
            READ (iUnit,iostat=datend) (den%pw(k,3),k=1,stars%ng3)
            IF (datend == 0) THEN
               IF (input%film) THEN
                  ALLOCATE(cdomvz(vacuum%nmz,vacuum%nvac))
                  READ (iUnit) ((cdomvz(i,iVac),i=1,vacuum%nmz),iVac=1,vacuum%nvac)
                  DO iVac = 1, vacuum%nvac
                     DO i = 1, vacuum%nmz
                        den%vacz(i,iVac,3) = REAL(cdomvz(i,iVac))
                        den%vacz(i,iVac,4) = AIMAG(cdomvz(i,iVac))
                     END DO
                  END DO
                  DEALLOCATE(cdomvz)
                  READ (iUnit) (((den%vacxy(i,j-1,iVac,3),i=1,vacuum%nmzxy),j=2,stars%ng2), iVac=1,vacuum%nvac)
               END IF
            ELSE
               ! (datend < 0)  =>  no off-diagonal magnetisation stored
               !                   in "rhomat_inp"
               IF (datend > 0) THEN
                  WRITE(*,*) 'datend = ', datend
                  CALL juDFT_error("density file has illegal state.",calledby ="readDensity")
               END IF
               den%pw(:,3) = CMPLX(0.0,0.0)
               IF (input%film) THEN
                  den%vacz(:,:,3:4) = 0.0
                  den%vacxy(:,:,:,3) = CMPLX(0.0,0.0)
               END IF
            END IF
         ELSE IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
            den%pw(:,3) = CMPLX(0.0,0.0)
            IF (input%film) THEN
               den%vacz(:,:,3:4) = 0.0
               den%vacxy(:,:,:,3) = CMPLX(0.0,0.0)
            END IF
         END IF
         CLOSE(iUnit)
      END IF

      INQUIRE(FILE='n_mmp_mat',EXIST=l_exist)
      IF(l_exist.AND.atoms%n_u.GT.0) THEN
         OPEN (69,file='n_mmp_mat',status='unknown',form='formatted')
         READ (69,'(7f20.13)',IOSTAT=ioStatus) den%mmpMat
         REWIND(69)
         numLines = 0
         DO
            READ (69,*,iostat=iofl)
            IF (iofl < 0) EXIT
            numLines = numLines + 1
         END DO
         IF (MOD(numLines,14*input%jspins).NE.0) THEN
            WRITE(*,*) 'The n_mmp_mat file could not be read.'
            WRITE(*,*) 'Was it an old style file with linear mixing parameter specification'
            WRITE(*,*) 'in the last line? Linear mixing for the density matrix can now be'
            WRITE(*,*) 'activated and specified in the inp.xml file.'
            CALL juDFT_error("strange n_mmp_mat-file...", calledby = "readDensity")
         END IF
         CLOSE(69)
      END IF

   END SUBROUTINE readDensity

   SUBROUTINE writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                           relCdnIndex,distance,fermiEnergy,l_qfix,den,inFilename)

      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_cell), INTENT(IN)     :: cell
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_potden),INTENT(INOUT) :: den

      INTEGER, INTENT (IN)      :: inOrOutCDN
      INTEGER, INTENT (IN)      :: relCdnIndex
      INTEGER, INTENT (IN)      :: archiveType
      REAL,    INTENT (IN)      :: fermiEnergy, distance
      LOGICAL, INTENT (IN)      :: l_qfix

      CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: inFilename

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
      LOGICAL           :: l_qfixTemp, l_CheckBroyd
      CHARACTER(LEN=30) :: archiveName
      CHARACTER(LEN=8)  :: dateString
      CHARACTER(LEN=10) :: timeString
      CHARACTER(LEN=10) :: zone

      COMPLEX, ALLOCATABLE :: cdomvz(:,:)

      CALL getIOMode(mode)
      CALL DATE_AND_TIME(dateString,timeString,zone)
      READ(dateString,'(i8)') date
      READ(timeString,'(i6)') time

      l_CheckBroyd = .TRUE.
      IF(PRESENT(inFilename)) l_CheckBroyd = .FALSE.

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex,inFilename)

         CALL checkAndWriteMetadataHDF(fileID, input, atoms, cell, vacuum, oneD, stars, sphhar, sym,&
                                       currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                       currentStepfunctionIndex,l_storeIndices,l_CheckBroyd)

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

         IF(vacuum%nvac.EQ.1) THEN
            den%vacz(:,2,:)=den%vacz(:,1,:)
            IF (sym%invs) THEN
               den%vacxy(:,:,2,:) = CONJG(den%vacxy(:,:,1,:))
            ELSE
               den%vacxy(:,:,2,:) = den%vacxy(:,:,1,:)
            END IF
         END IF

         CALL writeDensityHDF(input, fileID, archiveName, densityType, previousDensityIndex,&
                              currentStarsIndex, currentLatharmsIndex, currentStructureIndex,&
                              currentStepfunctionIndex,date,time,distance,fermiEnergy,l_qfix,&
                              den%iter+relCdnIndex,den)

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

         IF(PRESENT(inFilename)) filename = inFilename

         IF ((relCdnIndex.EQ.1).AND.(archiveType.EQ.CDN_ARCHIVE_TYPE_CDN1_const).AND.(den%iter.EQ.0)) THEN
            INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=l_exist)
            IF(l_exist) THEN
               CALL juDFT_error("Trying to generate starting density while a density exists.",calledby ="writeDensity")
            END IF
         END IF

         iUnit = 93
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),FORM='unformatted',STATUS='unknown')

         IF ((relCdnIndex.EQ.1).AND.(archiveType.EQ.CDN_ARCHIVE_TYPE_CDN1_const).AND.(den%iter.GE.1)) THEN
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
            vacuumTemp%nmzxy = vacuum%nmzxy
            vacuumTemp%nmz = vacuum%nmz
            vacuumTemp%dvac = vacuum%dvac
            vacuumTemp%delz = vacuum%delz
            starsTemp%ng2 = stars%ng2
            symTemp%invs2 = sym%invs2
            ALLOCATE (fpwTemp(stars%ng3,input%jspins))
            ALLOCATE (fzxyTemp(vacuum%nmzxyd,stars%ng2-1,2,input%jspins))
            ALLOCATE (frTemp(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins))
            ALLOCATE (fzTemp(vacuum%nmzd,2,input%jspins))

            !--->    generate name of file to hold the results of this iteration
            d1 = MOD(den%iter,10)
            d10 = MOD(INT((den%iter+0.5)/10),10)
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
            vacuumTemp%nmzxy = vacuum%nmzxy
            vacuumTemp%nmz = vacuum%nmz
            vacuumTemp%dvac = vacuum%dvac
            vacuumTemp%delz = vacuum%delz
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
                     iUnit,den%iter+relCdnIndex,den%mt,den%pw,den%vacz,den%vacxy)

         ! Write additional data if l_noco
         IF (archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const) THEN
            WRITE (iUnit) (den%pw(k,3),k=1,stars%ng3)
            IF (input%film) THEN
               ALLOCATE(cdomvz(vacuum%nmz,vacuum%nvac))
               DO iVac = 1, vacuum%nvac
                  DO i = 1, vacuum%nmz
                     cdomvz(i,iVac) = CMPLX(den%vacz(i,iVac,3),den%vacz(i,iVac,4))
                  END DO
               END DO
               WRITE (iUnit) ((cdomvz(i,iVac),i=1,vacuum%nmz),iVac=1,vacuum%nvac)
               WRITE (iUnit) (((den%vacxy(i,j-1,iVac,3),i=1,vacuum%nmzxy),j=2,stars%ng2), iVac=1,vacuum%nvac)
               DEALLOCATE(cdomvz)
            END IF
         END IF

         CLOSE(iUnit)
      END IF

      !write density matrix to n_mmp_mat_out file
      IF((inOrOutCDN.EQ.CDN_INPUT_DEN_const).AND.(relCdnIndex.EQ.1).AND.&
         ((archiveType.EQ.CDN_ARCHIVE_TYPE_CDN1_const).OR.(archiveType.EQ.CDN_ARCHIVE_TYPE_NOCO_const))) THEN
         IF(atoms%n_u.GT.0) THEN
            filename = 'n_mmp_mat'
            IF (mode.EQ.CDN_HDF5_MODE) THEN
               filename = 'n_mmp_mat_out'
            END IF
            IF(ANY(den%mmpMat(:,:,:,:).NE.0.0)) THEN
               IF ((mode.EQ.CDN_HDF5_MODE).AND..NOT.(input%ldauLinMix.AND.(input%ldauMixParam.EQ.0.0))) THEN
                  INQUIRE(file='n_mmp_mat',exist=l_exist)
                  IF(l_exist) THEN
                     CALL system('mv n_mmp_mat n_mmp_mat_old')
                     PRINT *,"n_mmp_mat moved to n_mmp_mat_old"
                  END IF
               END IF
               OPEN (69,file=TRIM(ADJUSTL(filename)),status='replace',form='formatted')
               WRITE (69,'(7f20.13)') den%mmpMat(:,:,:,:)
               CLOSE (69)
            END IF
         END IF
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

      CALL getIOMode(mode)

      eFermiPrev = 0.0
      l_error = .FALSE.
      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
         WRITE(archiveName,'(a,i0)') '/cdn-', readDensityIndex
         CALL peekDensityEntryHDF(fileID, archiveName, DENSITY_TYPE_NOCO_IN_const,&
                                  iter, starsIndex, latharmsIndex, structureIndex, stepfunctionIndex,&
                                  previousDensityIndex, jspins, date, time, distance, fermiEnergy, l_qfix)
         eFermiPrev = fermiEnergy
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

      REAL, INTENT(OUT) :: rhcs(:,:,:)!(atoms%jmtd,atoms%ntype,input%jspins)
      REAL, INTENT(OUT) :: tecs(:,:)!(atoms%ntype,input%jspins)
      REAL, INTENT(OUT) :: qints(:,:)!(atoms%ntype,input%jspins)

      INTEGER            :: mode, iUnit, iSpin, iAtom, i
      LOGICAL            :: l_exist
      CHARACTER(LEN=30)  :: filename

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER        :: currentStarsIndex,currentLatharmsIndex
      INTEGER        :: currentStructureIndex,currentStepfunctionIndex
      INTEGER        :: readDensityIndex, lastDensityIndex

      CALL getIOMode(mode)

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

      REAL, INTENT(IN) :: rhcs(:,:,:)!(atoms%jmtd,atoms%ntype,input%jspins)
      REAL, INTENT(IN) :: tecs(:,:)!(atoms%ntype,input%jspins)
      REAL, INTENT(IN) :: qints(:,:)!(atoms%ntype,input%jspins)

      INTEGER :: mode, iUnit, iSpin, iAtom, i

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif
      INTEGER        :: currentStarsIndex,currentLatharmsIndex
      INTEGER        :: currentStructureIndex,currentStepfunctionIndex
      INTEGER        :: readDensityIndex, lastDensityIndex

      CALL getIOMode(mode)

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

   SUBROUTINE storeStructureIfNew(input,stars, atoms, cell, vacuum, oneD, sym,mpi,sphhar,noco)

      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_oneD),INTENT(IN)    :: oneD
      TYPE(t_sym),INTENT(IN)     :: sym
      TYPE(t_mpi),INTENT(IN)      :: mpi
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_noco),INTENT(IN)     :: noco
      TYPE(t_stars),INTENT(IN)    :: stars

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

      CALL getIOMode(mode)

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
            CALL compareStructure(input, atoms, vacuum, cell, sym, inputTemp, atomsTemp, vacuumTemp, cellTemp, symTemp, l_same)
            IF(.NOT.l_same) THEN
               currentStructureIndex = currentStructureIndex + 1
               l_writeStructure = .TRUE.
            END IF
         END IF

 
         IF (l_writeStructure) THEN
            CALL writeStructureHDF(fileID, input, atoms, cell, vacuum, oneD, sym, currentStructureIndex,.TRUE.)
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

   SUBROUTINE transform_by_moving_atoms(mpi,stars,atoms,vacuum, cell, sym, sphhar,input,oned,noco)
     USE m_types
     USE m_qfix
     use m_fix_by_gaussian
     IMPLICIT NONE
     TYPE(t_mpi),INTENT(IN)      :: mpi
     TYPE(t_atoms),INTENT(IN)    :: atoms
     TYPE(t_sym),INTENT(IN)      :: sym
     TYPE(t_vacuum),INTENT(IN)   :: vacuum
     TYPE(t_sphhar),INTENT(IN)   :: sphhar
     TYPE(t_input),INTENT(IN)    :: input
     TYPE(t_oneD),INTENT(IN)     :: oneD
     TYPE(t_cell),INTENT(IN)     :: cell
     TYPE(t_noco),INTENT(IN)     :: noco
     TYPE(t_stars),INTENT(IN)    :: stars
     !Locals
     INTEGER :: archiveType
     LOGICAL :: l_qfix
     REAL    :: fermiEnergy,fix
     REAL    :: shifts(3,atoms%nat)

     TYPE(t_potden):: den

     TYPE(t_input)              :: inputTemp
     TYPE(t_atoms)              :: atomsTemp
     TYPE(t_cell)               :: cellTemp
     TYPE(t_vacuum)             :: vacuumTemp
     TYPE(t_oneD)               :: oneDTemp
     TYPE(t_sym)                :: symTemp
     

     INTEGER :: mode,currentStarsIndex,currentLatharmsIndex,currentStructureIndex
     INTEGER :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex,structureindex
     LOGICAL :: l_same,l_structure_by_shift

#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
      character(len=50) :: archivename
#endif
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
      integer :: ierr
#endif
      l_same=.TRUE.;l_structure_by_shift=.TRUE.
          
      CALL getIOMode(mode)
      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         IF (mpi%irank==0) THEN
            CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                 currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
            IF (currentStructureIndex>0.AND.lastdensityindex>0) THEN
               !Determine structure of last density
               WRITE(archivename,"(a,i0)") "cdn-",lastdensityindex
               CALL peekDensityEntryHDF(fileID, archivename, DENSITY_TYPE_IN_const, structureIndex=structureIndex)
               !Read that structure
               CALL readStructureHDF(fileID, inputTemp, atomsTemp, cellTemp, vacuumTemp, oneDTemp, symTemp, StructureIndex)
               CALL compareStructure(input, atoms, vacuum, cell, sym, inputTemp, atomsTemp, vacuumTemp, cellTemp, symTemp, l_same,l_structure_by_shift)
            ENDIF
            CALL closeCDNPOT_HDF(fileID)
         ENDIF
#ifdef CPP_MPI
         CALL mpi_bcast(l_same,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
         CALL mpi_bcast(l_structure_by_shift,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif
         IF (l_same.OR..NOT.l_structure_by_shift) RETURN ! nothing to do
         
         IF (mpi%irank==0) THEN
            WRITE(6,*) "Atomic movement detected, trying to adjust charge density" 
            
            !Calculate shifts
            shifts=atomsTemp%taual-atoms%taual
            
            !Determine type of charge
            archiveType = MERGE(CDN_ARCHIVE_TYPE_NOCO_const,CDN_ARCHIVE_TYPE_CDN1_const,noco%l_noco)
            !read the current density
            CALL den%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)
            CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
                 0,fermiEnergy,l_qfix,den)
         ENDIF
         !Now fix the density
         SELECT CASE(input%qfix)
         CASE (0,1) !just qfix the density
            IF (mpi%irank==0) WRITE(6,*) "Using qfix to adjust density"
            if (mpi%irank==0) CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,&
                 den,noco%l_noco,mpi%isize==1,force_fix=.TRUE.,fix=fix)
         CASE(2,3)
            if (mpi%irank==0) CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,&
                 den,noco%l_noco,mpi%isize==1,force_fix=.TRUE.,fix=fix,fix_pw_only=.true.)
         CASE(4,5)
            if (mpi%irank==0) CALL fix_by_gaussian(shifts,atoms,stars,mpi,sym,vacuum,sphhar,input,oned,cell,noco,den)
         CASE default
            CALL judft_error("Wrong choice of qfix in input")
         END SELECT
         !Now write the density to file
         IF (mpi%irank==0) CALL writedensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,CDN_INPUT_DEN_const,&
              0,-1.0,fermiEnergy,l_qfix,den)
         
#endif
      END IF
    END SUBROUTINE transform_by_moving_atoms

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

      CALL getIOMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         currentStarsIndex = currentStarsIndex + 1
         CALL writeStarsHDF(fileID, currentStarsIndex, currentStructureIndex, stars,.TRUE.)

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

      CALL getIOMode(mode)

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

      CALL getIOMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
#ifdef CPP_HDF
         CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                          currentStepfunctionIndex,readDensityIndex,lastDensityIndex)

         currentStepfunctionIndex = currentStepfunctionIndex + 1
         CALL writeStepfunctionHDF(fileID, currentStepfunctionIndex, currentStarsIndex, currentStructureIndex, stars,.TRUE.)
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

      CALL getIOMode(mode)

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
      CHARACTER(LEN=1000) :: numberString
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

      CALL getIOMode(mode)

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
      CHARACTER(LEN=1000) :: ddString
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
      IF (TRIM(ADJUSTL(ddString)).EQ.'allbutlast') THEN
         CALL getIOMode(mode)

         IF(mode.EQ.CDN_HDF5_MODE) THEN
            INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
            IF (l_exist) THEN
#ifdef CPP_HDF
               CALL openCDN_HDF(fileID,currentStarsIndex,currentLatharmsIndex,currentStructureIndex,&
                                currentStepfunctionIndex,readDensityIndex,lastDensityIndex)
               CALL closeCDNPOT_HDF(fileID)
               startNumber = 1
               endNumber = lastDensityIndex - 1
#endif
            END IF
         END IF
      ELSE
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
            READ(ddString(1:separatorIndex-1),'(i7)') startNumber
            READ(ddString(separatorIndex+1:LEN(TRIM(ADJUSTL(ddString)))),'(i7)') endNumber
         ELSE
            READ(ddString(1:LEN(TRIM(ADJUSTL(ddString)))),'(i7)') startNumber
            READ(ddString(1:LEN(TRIM(ADJUSTL(ddString)))),'(i7)') endNumber
         END IF
      END IF

      CALL getIOMode(mode)

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
         WRITE(*,*) 'Ignoring -delden command line argument.'
      END IF

      CALL juDFT_end("Selected densities deleted.")
      
   END SUBROUTINE deleteDensities

   SUBROUTINE getIOMode(mode)
      INTEGER, INTENT(OUT) :: mode

      mode = CDN_DIRECT_MODE
      !IF (juDFT_was_argument("-stream_cdn")) THEN
      !   mode=CDN_STREAM_MODE
      !END IF
      IF (.NOT.juDFT_was_argument("-no_cdn_hdf")) THEN !juDFT_was_argument("-hdf_cdn")) THEN
#ifdef CPP_HDF
         mode=CDN_HDF5_MODE
#else
!         WRITE(*,*) 'Code not compiled with HDF5 support.'
!         WRITE(*,*) 'Falling back to direct access.'
#endif
      END IF
   END SUBROUTINE getIOMode

   LOGICAL FUNCTION isDensityFilePresent(archiveType)

      INTEGER, INTENT(IN) :: archiveType

      LOGICAL             :: l_exist
      INTEGER             :: mode

      INTEGER        :: currentStarsIndex,currentLatharmsIndex,currentStructureIndex
      INTEGER        :: currentStepfunctionIndex,readDensityIndex,lastDensityIndex
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
#endif

      CALL getIOMode(mode)

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

      CALL getIOMode(mode)

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
