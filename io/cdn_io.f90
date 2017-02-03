!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
   IMPLICIT NONE

   PRIVATE
   PUBLIC readDensity, writeDensity
   PUBLIC isDensityFilePresent, isCoreDensityPresent
   PUBLIC readCoreDensity, writeCoreDensity
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

   SUBROUTINE readDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                          relCdnIndex,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_oneD),INTENT(IN)   :: oneD

      INTEGER, INTENT (IN)      :: inOrOutCDN
      INTEGER, INTENT (IN)      :: relCdnIndex
      INTEGER, INTENT (OUT)     :: iter
      INTEGER, INTENT (IN)      :: archiveType

      !     ..
      !     .. Array Arguments ..
      COMPLEX, INTENT (OUT) :: fpw(stars%n3d,input%jspins), fzxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)
      COMPLEX, INTENT (OUT) :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)
      REAL,    INTENT (OUT) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins), fz(vacuum%nmzd,2,input%jspins)

      ! local variables
      INTEGER            :: mode, datend, k, i, iVac, j, iUnit
      LOGICAL            :: l_exist, l_rhomatFile
      CHARACTER(LEN=30)  :: filename

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from cdn.hdf and exit subroutine

            RETURN
         ELSE
            WRITE(*,*) 'cdn.hdf file not found.'
            WRITE(*,*) 'Falling back to stream access file cdn.str.'
            mode = CDN_STREAM_MODE
         END IF
      END IF

      IF(mode.EQ.CDN_STREAM_MODE) THEN
         INQUIRE(FILE='cdn.str',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from cdn.str and exit subroutine

            RETURN
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
      END IF

      CLOSE(iUnit)

   END SUBROUTINE readDensity

   SUBROUTINE writeDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,archiveType,inOrOutCDN,&
                           relCdnIndex,iter,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_oneD),INTENT(IN)   :: oneD

      INTEGER, INTENT (IN)      :: inOrOutCDN
      INTEGER, INTENT (IN)      :: relCdnIndex, iter
      INTEGER, INTENT (IN)      :: archiveType
      !     ..
      !     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: fpw(stars%n3d,input%jspins), fzxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)
      COMPLEX, INTENT (IN) :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)
      REAL,    INTENT (IN) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins), fz(vacuum%nmzd,2,input%jspins)

      TYPE(t_stars)        :: starsTemp
      TYPE(t_vacuum)       :: vacuumTemp
      TYPE(t_atoms)        :: atomsTemp
      TYPE(t_sphhar)       :: sphharTemp
      TYPE(t_input)        :: inputTemp
      TYPE(t_sym)          :: symTemp

      COMPLEX, ALLOCATABLE :: fpwTemp(:,:), fzxyTemp(:,:,:,:)
      REAL, ALLOCATABLE    :: frTemp(:,:,:,:), fzTemp(:,:,:)

      INTEGER           :: mode, iterTemp, k, i, iVac, j, iUnit
      INTEGER           :: d1, d10, asciioffset, iUnitTemp
      LOGICAL           :: l_exist
      CHARACTER(len=30) :: filename
      CHARACTER(len=5)  :: cdnfile

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         ! Write density to cdn.hdf file
         STOP 'CDN_HDF5_MODE not yet implemented!'
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
            starsTemp%n3d = stars%n3d
            inputTemp%jspins = input%jspins
            vacuumTemp%nmzxyd = vacuum%nmzxyd
            starsTemp%n2d = stars%n2d
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
            ALLOCATE (fpwTemp(stars%n3d,input%jspins))
            ALLOCATE (fzxyTemp(vacuum%nmzxyd,stars%n2d-1,2,input%jspins))
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
            starsTemp%n3d = stars%n3d
            inputTemp%jspins = input%jspins
            vacuumTemp%nmzxyd = vacuum%nmzxyd
            starsTemp%n2d = stars%n2d
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
            ALLOCATE (fpwTemp(stars%n3d,input%jspins))
            ALLOCATE (fzxyTemp(vacuum%nmzxyd,stars%n2d-1,2,input%jspins))
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

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from cdn.hdf and exit subroutine

            RETURN
         ELSE
            WRITE(*,*) 'cdn.hdf file not found.'
            WRITE(*,*) 'Falling back to stream access file cdn.str.'
            mode = CDN_STREAM_MODE
         END IF
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

      CALL getMode(mode)

      IF(mode.EQ.CDN_HDF5_MODE) THEN
         ! Write core density to cdn.hdf file
         STOP 'CDN_HDF5_MODE not yet implemented!'
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


   SUBROUTINE getMode(mode)
      INTEGER, INTENT(OUT) :: mode

      mode = CDN_DIRECT_MODE
      IF (juDFT_was_argument("-stream_cdn")) mode=CDN_STREAM_MODE
      IF (juDFT_was_argument("-hdf_cdn")) mode=CDN_HDF5_MODE
   END SUBROUTINE getMode

   LOGICAL FUNCTION isDensityFilePresent(archiveType)

      INTEGER, INTENT(IN) :: archiveType

      LOGICAL             :: l_exist
      INTEGER             :: mode

      CALL getMode(mode)

      IF (mode.EQ.CDN_HDF5_MODE) THEN
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF(l_exist) THEN
            isDensityFilePresent = l_exist
            RETURN
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
         INQUIRE(FILE='cdn.hdf',EXIST=l_exist)
         IF(l_exist) THEN
            STOP 'Not yet implemented!'
            RETURN
         END IF
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
