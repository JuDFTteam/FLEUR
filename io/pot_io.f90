!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This module is a wrapper for the potential I/O
!!!
!!!                             GM'17
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_pot_io

   USE m_types
   USE m_juDFT
   USE m_loddop
   USE m_wrtdop
   IMPLICIT NONE

   PRIVATE
   PUBLIC readPotential, writePotential
   PUBLIC POT_ARCHIVE_TYPE_TOT_const, POT_ARCHIVE_TYPE_COUL_const
   PUBLIC POT_ARCHIVE_TYPE_X_const

   INTEGER,          PARAMETER :: POT_ARCHIVE_TYPE_TOT_const = 1
   INTEGER,          PARAMETER :: POT_ARCHIVE_TYPE_COUL_const = 2
   INTEGER,          PARAMETER :: POT_ARCHIVE_TYPE_X_const = 3

   INTEGER,          PARAMETER :: POT_DIRECT_MODE = 1
   INTEGER,          PARAMETER :: POT_STREAM_MODE = 2
   INTEGER,          PARAMETER :: POT_HDF5_MODE   = 3

   CONTAINS

   SUBROUTINE readPotential(stars,vacuum,atoms,sphhar,input,sym,archiveType,&
                            iter,fr,fpw,fz,fzxy)

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym

      INTEGER, INTENT (OUT)     :: iter
      INTEGER, INTENT (IN)      :: archiveType

      !     ..
      !     .. Array Arguments ..
      COMPLEX, INTENT (OUT) :: fpw(stars%n3d,input%jspins), fzxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)
      REAL,    INTENT (OUT) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins), fz(vacuum%nmzd,2,input%jspins)

      ! local variables
      INTEGER           :: mode, iUnit
      LOGICAL           :: l_exist
      CHARACTER(len=30) :: filename

      CALL getMode(mode)

      IF(mode.EQ.POT_HDF5_MODE) THEN
         INQUIRE(FILE='pot.hdf',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from pot.hdf and exit subroutine

            RETURN
         ELSE
            WRITE(*,*) 'pot.hdf file not found.'
            WRITE(*,*) 'Falling back to stream access file pot.str.'
            mode = POT_STREAM_MODE
         END IF
      END IF

      IF(mode.EQ.POT_STREAM_MODE) THEN
         INQUIRE(FILE='pot.str',EXIST=l_exist)
         IF (l_exist) THEN
            !load density from cdn.str and exit subroutine

            RETURN
         ELSE
            WRITE(*,*) 'pot.str file not found.'
            WRITE(*,*) 'Falling back to direct access file.'
            mode = POT_DIRECT_MODE
         END IF
      END IF

      IF (mode.EQ.POT_DIRECT_MODE) THEN
         filename = 'illegalPotentialArchive'
         IF (archiveType.EQ.POT_ARCHIVE_TYPE_TOT_const) THEN
            filename = 'pottot'
         END IF
         IF (archiveType.EQ.POT_ARCHIVE_TYPE_COUL_const) THEN
            filename = 'potcoul'
         END IF
         IF (archiveType.EQ.POT_ARCHIVE_TYPE_X_const) THEN
            filename = 'potx'
         END IF

         INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=l_exist)
         IF(.NOT.l_exist) THEN
            CALL juDFT_error("potential file "//TRIM(ADJUSTL(filename))//" missing",calledby ="readPotential")
         END IF

         iUnit = 11
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),form='unformatted',status='unknown')

         CALL loddop(stars,vacuum,atoms,sphhar,input,sym,&
                     iUnit,iter,fr,fpw,fz,fzxy)
         CLOSE(iUnit)

      END IF

   END SUBROUTINE readPotential

   SUBROUTINE writePotential(stars,vacuum,atoms,sphhar,input,sym,archiveType,&
                             iter,fr,fpw,fz,fzxy)

      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sym),INTENT(IN)    :: sym

      INTEGER, INTENT (IN)      :: iter
      INTEGER, INTENT (IN)      :: archiveType
      !     ..
      !     .. Array Arguments ..
      COMPLEX, INTENT (IN) :: fpw(stars%n3d,input%jspins), fzxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)
      REAL,    INTENT (IN) :: fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins), fz(vacuum%nmzd,2,input%jspins)

      ! local variables
      INTEGER           :: mode, iUnit
      LOGICAL           :: l_exist
      CHARACTER(len=30) :: filename

      CALL getMode(mode)

      IF(mode.EQ.POT_HDF5_MODE) THEN
         ! Write potential to pot.hdf file
         STOP 'POT_HDF5_MODE not yet implemented!'
      ELSE IF(mode.EQ.POT_STREAM_MODE) THEN
         ! Write potential to pot.str file
         STOP 'POT_STREAM_MODE not yet implemented!'
      ELSE
         ! Direct mode
         filename = 'illegalPotentialArchive'
         IF (archiveType.EQ.POT_ARCHIVE_TYPE_TOT_const) THEN
            filename = 'pottot'
         END IF
         IF (archiveType.EQ.POT_ARCHIVE_TYPE_COUL_const) THEN
            filename = 'potcoul'
         END IF
         IF (archiveType.EQ.POT_ARCHIVE_TYPE_X_const) THEN
            filename = 'potx'
         END IF

         iUnit = 11
         OPEN (iUnit,file=TRIM(ADJUSTL(filename)),form='unformatted',status='unknown')
         CALL wrtdop(stars,vacuum,atoms,sphhar,input,sym,&
                     iUnit,iter,fr,fpw,fz,fzxy)
         CLOSE(iUnit)
      END IF

   END SUBROUTINE writePotential

   SUBROUTINE getMode(mode)
      INTEGER, INTENT(OUT) :: mode

      mode = POT_DIRECT_MODE
      IF (juDFT_was_argument("-stream_pot")) mode=POT_STREAM_MODE
      IF (juDFT_was_argument("-hdf_pot")) mode=POT_HDF5_MODE
   END SUBROUTINE getMode

END MODULE m_pot_io
