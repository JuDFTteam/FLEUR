!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_dfpt_cdn_io
   use m_juDFT
   use m_constants, only: POTDEN_TYPE_DEN
   use m_cdn_io, only: readDensity, writeDensity, isDensityArchivePresent, CDN_INPUT_DEN_const, &
                       CDN_ARCHIVE_TYPE_CDN1_const, CDN_ARCHIVE_TYPE_NOCO_const, CDN_ARCHIVE_TYPE_FFN_const
   use m_types_atoms
   use m_types_cell
   use m_types_input
   use m_types_noco
   use m_types_potden
   use m_types_sphhar
   use m_types_stars
   use m_types_sym
   use m_types_vacuum

   implicit none
   private
   public :: dfptDensityArchiveName, dfptArchiveType, writeDfptDensity, readDfptDensity, isDfptDensityPresent
   public :: storeDfptDensity, clearDfptWorkFile

   character(len=*), parameter :: DFPT_CDN_FILE = 'dfpt_cdn'
   character(len=*), parameter :: DFPT_WORK_FILE = 'dfpt_cdn_current'

contains

   function dfptDensityArchiveName(kind, iQ, iDtype, iDir) result(archiveName)
      !! Archive name /dfpt/<kind>/q_NNNNNN_b_NNN_j_N of one density perturbation.

      character(len=*), intent(in) :: kind
      integer,          intent(in) :: iQ, iDtype, iDir
      character(len=64)            :: archiveName

      if (iQ < 0 .or. iQ > 999999 .or. iDtype < 0 .or. iDtype > 999 .or. iDir < 0 .or. iDir > 9) then
         call judft_error("DFPT perturbation index out of range for the archive name.", calledby="dfptDensityArchiveName")
      end if
      write(archiveName, '(a,i6.6,a,i3.3,a,i1)') '/dfpt/'//trim(kind)//'/q_', iQ, '_b_', iDtype, '_j_', iDir
   end function dfptDensityArchiveName

   integer function dfptArchiveType(noco)
      !! Archive type of the DFPT density perturbations for the given magnetism setup.

      type(t_noco), intent(in) :: noco

      dfptArchiveType = CDN_ARCHIVE_TYPE_CDN1_const
      if (any(noco%l_unrestrictMT)) then
         dfptArchiveType = CDN_ARCHIVE_TYPE_FFN_const
      else if (noco%l_noco) then
         dfptArchiveType = CDN_ARCHIVE_TYPE_NOCO_const
      end if
   end function dfptArchiveType

   subroutine writeDfptDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,archiveType,archiveName,qvec,distance,fermiEnergy,den)
      !! Writes the density perturbation of the running perturbation into the working file dfpt_cdn_current.hdf.

      type(t_stars),    intent(in)    :: stars
      type(t_noco),     intent(in)    :: noco
      type(t_vacuum),   intent(in)    :: vacuum
      type(t_atoms),    intent(in)    :: atoms
      type(t_cell),     intent(in)    :: cell
      type(t_sphhar),   intent(in)    :: sphhar
      type(t_input),    intent(in)    :: input
      type(t_sym),      intent(in)    :: sym
      integer,          intent(in)    :: archiveType
      character(len=*), intent(in)    :: archiveName
      real,             intent(in)    :: qvec(3)
      real,             intent(in)    :: distance, fermiEnergy
      type(t_potden),   intent(inout) :: den

      logical :: l_backup

      call backupFile(DFPT_WORK_FILE, l_backup)
      call writeDensity(stars, noco, vacuum, atoms, cell, sphhar, input, sym, archiveType, CDN_INPUT_DEN_const, 1, distance, fermiEnergy, &
                        -1.0, -1.0, .true., den, inFilename=DFPT_WORK_FILE, inArchive=archiveName, qvec=qvec)
      if (l_backup) call system('rm '//DFPT_WORK_FILE//'_old.hdf')
   end subroutine writeDfptDensity

   subroutine readDfptDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,archiveType,archiveName,qvec,fermiEnergy,lastDistance,den)
      !! Reads the density perturbation of q-vector qvec, from the working file if it holds it, otherwise from dfpt_cdn.hdf.

      type(t_stars),    intent(in)    :: stars
      type(t_noco),     intent(in)    :: noco
      type(t_vacuum),   intent(in)    :: vacuum
      type(t_atoms),    intent(in)    :: atoms
      type(t_cell),     intent(in)    :: cell
      type(t_sphhar),   intent(in)    :: sphhar
      type(t_input),    intent(in)    :: input
      type(t_sym),      intent(in)    :: sym
      integer,          intent(in)    :: archiveType
      character(len=*), intent(in)    :: archiveName
      real,             intent(in)    :: qvec(3)
      real,             intent(out)   :: fermiEnergy, lastDistance
      type(t_potden),   intent(inout) :: den

      character(len=len(DFPT_WORK_FILE)) :: filename
      logical                            :: l_qfix

      filename = DFPT_CDN_FILE
      if (isDensityArchivePresent(DFPT_WORK_FILE, archiveType, archiveName, qvec)) then
         filename = DFPT_WORK_FILE
      else if (isDensityArchivePresent(DFPT_CDN_FILE, archiveType, archiveName) .and. &
               .not. isDensityArchivePresent(DFPT_CDN_FILE, archiveType, archiveName, qvec)) then
         call judft_error("Density perturbation "//trim(archiveName)//" in "//DFPT_CDN_FILE//".hdf belongs to a different q-vector.", &
                          calledby="readDfptDensity")
      end if
      call readDensity(stars, noco, vacuum, atoms, cell, sphhar, input, sym, archiveType, CDN_INPUT_DEN_const, 0, fermiEnergy, &
                       lastDistance, l_qfix, den, inFilename=trim(filename), inArchive=archiveName, qvec=qvec)
   end subroutine readDfptDensity

   logical function isDfptDensityPresent(archiveType, archiveName, qvec)
      !! Checks whether the working file or dfpt_cdn.hdf holds the density perturbation archiveName of q-vector qvec.

      integer,          intent(in) :: archiveType
      character(len=*), intent(in) :: archiveName
      real,             intent(in) :: qvec(3)

      isDfptDensityPresent = isDensityArchivePresent(DFPT_WORK_FILE, archiveType, archiveName, qvec)
      if (.not. isDfptDensityPresent) isDfptDensityPresent = isDensityArchivePresent(DFPT_CDN_FILE, archiveType, archiveName, qvec)
   end function isDfptDensityPresent

   subroutine storeDfptDensity(stars,vacuum,atoms,cell,sphhar,input,sym,noco,archiveName,qvec)
      !! Copies a converged density perturbation from the working file into dfpt_cdn.hdf.

      type(t_stars),    intent(in) :: stars
      type(t_vacuum),   intent(in) :: vacuum
      type(t_atoms),    intent(in) :: atoms
      type(t_cell),     intent(in) :: cell
      type(t_sphhar),   intent(in) :: sphhar
      type(t_input),    intent(in) :: input
      type(t_sym),      intent(in) :: sym
      type(t_noco),     intent(in) :: noco
      character(len=*), intent(in) :: archiveName
      real,             intent(in) :: qvec(3)

      type(t_potden) :: den
      integer        :: archiveType
      real           :: fermiEnergy, lastDistance
      logical        :: l_qfix, l_backup

      archiveType = dfptArchiveType(noco)
      if (.not. isDensityArchivePresent(DFPT_WORK_FILE, archiveType, archiveName, qvec)) return

      call den%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN, l_dfpt=.true.)
      call readDensity(stars, noco, vacuum, atoms, cell, sphhar, input, sym, archiveType, CDN_INPUT_DEN_const, 0, fermiEnergy, &
                       lastDistance, l_qfix, den, inFilename=DFPT_WORK_FILE, inArchive=archiveName, qvec=qvec)

      call backupFile(DFPT_CDN_FILE, l_backup)
      call writeDensity(stars, noco, vacuum, atoms, cell, sphhar, input, sym, archiveType, CDN_INPUT_DEN_const, 0, lastDistance, fermiEnergy, &
                        -1.0, -1.0, l_qfix, den, inFilename=DFPT_CDN_FILE, inArchive=archiveName, qvec=qvec)
      if (l_backup) call system('rm '//DFPT_CDN_FILE//'_old.hdf')
   end subroutine storeDfptDensity

   subroutine clearDfptWorkFile()
      !! Removes the working file of the finished perturbation.

      logical :: l_exist

      inquire(file=DFPT_WORK_FILE//'.hdf', exist=l_exist)
      if (l_exist) call system('rm '//DFPT_WORK_FILE//'.hdf')
      inquire(file=DFPT_WORK_FILE//'_old.hdf', exist=l_exist)
      if (l_exist) call system('rm '//DFPT_WORK_FILE//'_old.hdf')
   end subroutine clearDfptWorkFile

   subroutine backupFile(filename, l_backup)
      !! Moves filename.hdf to filename_old.hdf and copies it back, so that a complete file survives a crash during the next write.

      character(len=*), intent(in)  :: filename
      logical,          intent(out) :: l_backup

      logical :: l_exist

      inquire(file=filename//'_old.hdf', exist=l_exist)
      if (l_exist) call system('rm '//filename//'_old.hdf')
      inquire(file=filename//'.hdf', exist=l_backup)
      ! the copy gives the file to be written a new inode, so stale HDF5 file locks on the old one do not block it
      if (l_backup) call system('mv '//filename//'.hdf '//filename//'_old.hdf && cp '//filename//'_old.hdf '//filename//'.hdf')
   end subroutine backupFile

end module m_dfpt_cdn_io
